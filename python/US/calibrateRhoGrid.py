r""" Calibrate the US model across a grid of CRRA parameter values.

Run:  .venv\Scripts\python.exe python\US\calibrateRhoGrid.py [options]
      ... --lo 0.5 --hi 2.0 --step 0.1        the default sweep
      ... --commonX                           the common-X variant (its own csv and pickle directory)
      ... --force                             re-solve points already in the csv
      ... --rhos 0.5 1 2                      the grid as a list (instead of --lo/--hi/--step)
      ... --theta VARIANT [--root DIR]        theta imposed as db['θ0'] = the U.S. row of data/oecd26_theta.csv
                                              under VARIANT (ModelUS.getθ, oecd26.thetaFromFile); csv and pickles
                                              under DIR/calibration/, DIR the variant's directory (default
                                              results/oecd26_theta_<variant>/, oecd26.thetaDir),
                                              never results/calibration/; its rows carry thetaSource and θ0, and
                                              a csv of another θ0 is refused unless --force.
      ... --theta VARIANT --main              the pipeline's main variant (python/paper/config.py US['theta'],
                                              oecd26.thetaRoot): the same rows into results/calibration/.
      ... --theta index                       θ from the workbook's RR0 (no thetaSource column); without --theta that
                                              run is refused into results/ (oecd26.indexTheta): use --out and --pkldir.
A relative --out, --pkldir or --root is the caller's working directory's (the script itself runs in python/US).

The march is anchored at rho=1, where the LOG solver applies and no warm start is needed, and walks
outward in both directions (gridsearch.continuation.marchGrid). Each CRRA point is a full backward
recursion inside an outer 2-D root, so the sweep is minutes rather than the hours the Argentina one
costs -- but it is built to survive being interrupted all the same:

  * results/calibration/US_rhoGrid.csv is the record, rewritten after *every* point. A run started over
    reads it back, and any rho already present is returned from it without re-solving -- which also means
    those points still seed the extrapolation, so a resumed march is no worse warm-started than an
    uninterrupted one. A row is resumed only at the run's settings (staleSettings, finding #13): one that
    records another n, ns (CRRA rows), smoothKnots or interpKind is re-solved.
  * results/calibration/instancesUS/*.pkl is a cache, not the record: one pickled model per point, holding
    the converged db. Reconstructing the same thing from the csv costs one _calSetPars and one solve, so
    nothing is lost if these are deleted. Full solutions (policy functions) are deliberately NOT stored:
    CRRA's are not pickleable (gridsearch.interp returns closures) and they regenerate in seconds.

    The pickle directory is per-variant and separate from the Argentina sweep's `instances/`, whose
    filenames are the rho alone: pointing two sweeps at one directory silently overwrites wherever the
    grids share a value. Same reason --out and --pkldir must be changed together.

--verify re-evaluates the outer residual at the converged parameters on a finer inner grid and stores it
beside the point. It is what separates a converged calibration from a resolved one: a point whose
verifyResidual is not small is converged but not resolved, and should not be read as a result.

Two things to know before reading that column. The rho=1 row is NOT informative -- the LOG solver has no
inner state grid to refine, and solveVectorized does not touch the tau grid either, so its verifyResidual
comes back equal to its residual by construction (measured: 1.24e-13 for both). And on the CRRA rows it
OVERSTATES the uncertainty in the answer: at rho=0.5 it sits at ~7e-4 across every ns from 100 to 300,
while beta itself is settled to 1e-4 relative over ns>=150. It measures how much the residual moves under
a 1.5x grid change at fixed parameters, which does not vanish just because the parameters have converged.
Read it as a level across points -- a point an order of magnitude above its neighbours is the one to look
at -- not as an error bar on beta.

Note what is NOT swept here. beta and omega are the only calibrated parameters; theta is closed-form from
the replacement-rate ratio (or imposed, --theta) and does not depend on rho, and under --commonX the hours-unit parameter X is
recovered after the root rather than inside it (docs eq:calibration:Xsolve). So the columns that move
across rho are beta, omega and the untargeted savings rate -- R and tau are targets and must not move.
"""
import os, sys, argparse, pickle, subprocess, datetime
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                 # where a relative --out, --pkldir or --root points: the chdir below moves to HERE

# Progress has to appear as it happens rather than when the process exits: whenever stdout is a pipe or a
# log file it is block-buffered by default, which hides the whole run. Encoding too, not only buffering --
# every progress line carries Greek parameter names, and a redirected stdout defaults to the ANSI codepage
# on Windows, which raises UnicodeEncodeError on the first point and takes the march down with it.
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)                                        # test.py reads data/ relative to its own location

import test as testmod
from model import ModelUS
from gridsearch import continuation
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'calibration')

# Ordered so the csv reads as: what was solved, what it cost, whether to believe it, then the answer.
COLUMNS = ['ρ', 'preferences', 'requested', 'residual', 'verifyResidual',
           'β', 'ω', 'X',
           'R', 'τ', 'sr', 'h', 'hbar', 'nfev', 'time', 'n', 'ns', 'smoothKnots', 'interpKind',
           'nEqMax', 'nCandMax', 'nFallback',
           'x0', 'x1', 'commit', 'timestamp']
THETACOLUMNS = ['thetaSource', 'θ0']      # appended under --theta only (oecd26.usBeta reads them)
# R and τ are targets, so they are constants down the csv by construction -- kept anyway, because a column
# that is supposed to be constant is the cheapest possible check that a point converged to the right thing
# rather than merely converged. sr/h/hbar are the untargeted quantities and are what actually vary.


def gitCommit():
    """ Short hash of the code that produced a row -- a sweep outlives several edits of the solver. """
    try:
        return subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd = REPO, capture_output = True,
                              text = True, timeout = 10).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def toRow(rec, requested, commit):
    """ calibratePoint's record -> one flat csv row, with the tax counts of the PEE solve at the calibrated
    point (the record's 'multiplicity' as nEqMax, nCandMax, nFallback; policy.multiplicityColumns). """
    g = rec['gridSettings']
    row = {k: rec.get(k) for k in ('ρ', 'preferences', 'residual', 'verifyResidual', 'β', 'ω', 'X',
                                   'R', 'τ', 'sr', 'h', 'hbar', 'nfev', 'time')}
    row |= {'requested': requested, 'n': g.get('n'), 'ns': g.get('ns'),
            'smoothKnots': g.get('smoothKnots'), 'interpKind': g.get('interpKind'), 'commit': commit,
            'timestamp': datetime.datetime.now().isoformat(timespec = 'seconds')}
    row |= {f'x{i}': v for i, v in enumerate(rec['x'])}
    return row | multiplicityColumns(rec.get('multiplicity'))


def readDone(path):
    """ {rho: x} for the points already solved, as the resume cache. Keyed on the rounded value so a grid
    rebuilt from the same --lo/--step matches it exactly rather than by floating-point luck. """
    if not os.path.exists(path):
        return {}, pd.DataFrame(columns = COLUMNS)
    df = pd.read_csv(path)
    done = {round(float(r['ρ']), 6): np.array([r['x0'], r['x1']], dtype = float)
            for _, r in df.iterrows() if np.isfinite(r.get('residual', np.nan))}
    return done, df


def rowSettings(a):
    """ The settings a row records and a resumed row must match (toRow's n, ns, smoothKnots, interpKind), from the
    parsed arguments: smoothKnots 0 is None (the adaptive smoother). """
    return {'n': a.n, 'ns': a.ns, 'smoothKnots': a.smoothKnots or None, 'interpKind': a.interpKind}


def _sameSetting(v, w):
    """ A recorded value v (from the csv) against the run's w: None is NaN, numbers to 1e-12 relative. """
    if w is None:
        return bool(pd.isna(v))
    if pd.isna(v):
        return False
    if isinstance(w, str):
        return str(v) == w
    return abs(float(v) - float(w)) <= 1e-12*max(1., abs(float(w)))


def staleSettings(df, want, grid, legacy = ()):
    """ The resume check of a sweep csv (finding #13) over the solved rows (finite residual) at a rho of `grid`.
    Returns ({rho: [the columns that differ]}, [the rhos of rows that predate a recorded setting]).

    `want` maps a column to the run's value, or to a callable of rho. A row records the grid where its 'n' is set:
    then n, smoothKnots, interpKind and, on a CRRA row, ns (the LOG solver has no state grid; its rows record the
    policy's default) must equal want's. A `legacy` column, added to the rows after the csvs were first written, is
    checked where the row has a value. A row without the grid or without a legacy column's value is accepted. """
    if df is None or not len(df):
        return {}, []
    on = set(np.round(np.asarray(grid, dtype = float), 6))
    stale, old = {}, []
    for r in df.to_dict('records'):
        ρ = round(float(r['ρ']), 6)
        if ρ not in on or not np.isfinite(pd.to_numeric(r.get('residual', np.nan), errors = 'coerce')):
            continue
        cols = [c for c in want if c in legacy or not pd.isna(r.get('n', np.nan))]
        cols = [c for c in cols if not (c == 'ns' and r.get('preferences') != 'CRRA')]
        unrecorded = pd.isna(r.get('n', np.nan)) or any(pd.isna(r.get(c, np.nan)) for c in legacy)
        diff = [c for c in cols if not (c in legacy and pd.isna(r.get(c, np.nan)))
                and not _sameSetting(r.get(c, np.nan), want[c](ρ) if callable(want[c]) else want[c])]
        if diff:
            stale[ρ] = diff
        elif unrecorded:
            old.append(ρ)
    return stale, old


def resumeAt(done, df, want, grid, legacy = ()):
    """ Drop from `done` the rows staleSettings finds at other settings (they are re-solved) and print what it found:
    the stale rows, and one notice for the rows that predate a recorded setting (resumed unchecked). """
    stale, old = staleSettings(df, want, grid, legacy)
    for ρ, cols in stale.items():
        done.pop(ρ, None)
        print('  rho={:<6} stale: recorded at another {} -- re-solved'.format(ρ, ', '.join(cols)))
    if old:
        print('notice: {} row(s) (rho {}) predate the recorded settings {}: resumed without the check'.format(
            len(old), ', '.join('{:g}'.format(v) for v in old), ', '.join(legacy) if legacy else 'n, ns, smoothKnots, interpKind'))
    return stale, old


def thetaRun(variant, commonX, out = None, pkldir = None, root = None, main = False):
    """ --theta VARIANT: (θ0, csv path, pickle directory), the last two by default under root/calibration/, root
    the variant's directory (oecd26.thetaDir(variant) by default; a relative path is the repo's), results/ under
    `main` (oecd26.thetaRoot: VARIANT must be the pipeline's main variant). SystemExit with --main and --root, where
    the variant has no U.S. θ, or, without --main, for a path under results/calibration/. """
    import oecd26
    if main and root is not None:
        raise SystemExit('--main writes results/calibration/: not with --root')
    try:
        θ0 = oecd26.thetaFromFile('USA', variant)
        top = oecd26.thetaRoot(variant, main) if root is None else os.path.join(REPO, root)
    except (ValueError, FileNotFoundError) as e:          # NoTheta is a ValueError
        raise SystemExit('--theta {}: {}'.format(variant, e))
    root = os.path.join(top, 'calibration')
    out = out or os.path.join(root, 'US' + ('_rhoGridCommonX' if commonX else '_rhoGrid') + '.csv')
    pkldir = pkldir or os.path.join(root, 'instancesUS' + ('CommonX' if commonX else ''))
    for q in (out, pkldir):
        if not main and os.path.normcase(os.path.abspath(q)).startswith(os.path.normcase(OUTDIR + os.sep)):
            raise SystemExit('--theta {}: {} is under results/calibration/, which holds the main θ variant '
                             '(--main) only'.format(variant, q))
    return θ0, out, pkldir


def rhoGrid(a):
    """ The sweep's grid: --rhos when given, else --lo to --hi by --step; rounded to 6 dp (the csv's key). """
    if a.rhos:
        return np.round(np.array(sorted(set(a.rhos)), dtype = float), 6)
    return np.round(np.arange(a.lo, a.hi + 0.5*a.step, a.step), 6)


def main():
    p = argparse.ArgumentParser(description = 'Calibrate the US model across a grid of rho.')
    p.add_argument('--lo', type = float, default = 0.5)
    p.add_argument('--hi', type = float, default = 2.0)
    p.add_argument('--step', type = float, default = 0.1)
    p.add_argument('--anchor', type = float, default = 1.0)
    p.add_argument('--commonX', action = 'store_true',
                   help = 'the common-X calibration variant (targets average hours as well)')
    p.add_argument('--n', type = int, default = 101, help = 'nodes on the tau grid')
    # 150, not policy.py's solve default of 50: the outer root is far more demanding of the inner grid
    # than a single PEE solve. At rho=0.5 beta settles to 1e-4 relative from ns=150 on and is 2.4% off at
    # ns=50 -- see policy.py's _gridSettings for the measured sequence.
    p.add_argument('--ns', type = int, default = 150, help = 'nodes on the CRRA state grid (LOG ignores it)')
    p.add_argument('--verify', type = int, default = 225,
                   help = 'refined ns for the resolution check; 0 to skip')
    p.add_argument('--verifyN', type = int, default = 151, help = 'refined tau nodes for the same check')
    # Both are well-posedness settings rather than tuning knobs -- see policy.py's _gridSettings for the
    # measurement behind the defaults. Exposed so a sweep can reproduce the pre-fix behaviour on demand
    # (--smoothKnots 0 selects the adaptive smoother), and recorded in the csv either way.
    p.add_argument('--smoothKnots', type = int, default = 4,
                   help = 'fixed-knot policy smoother: interior knot every m-th node. 0 selects the '
                          'adaptive smoother, whose knot COUNT is chosen from the data and makes the '
                          'outer residual jitter rather than converge under refinement.')
    p.add_argument('--interpKind', default = 'linear', choices = ('linear', 'cubic', 'pchip'),
                   help = "continuation interpolant. 'cubic' is NOT recommended here: it does not change "
                          'the refinement behaviour and fails to converge at some (rho, ns).')
    p.add_argument('--degree', type = int, default = 1, help = 'extrapolation degree for the warm start')
    p.add_argument('--maxHalvings', type = int, default = 2)
    p.add_argument('--force', action = 'store_true', help = 're-solve points already present in the csv')
    p.add_argument('--theta', default = None,
                   help = "a variant of data/oecd26_theta.csv: θ imposed at the U.S. row's value (db['θ0']); 'index' the "
                          "workbook's RR0 (required for that into results/)")
    p.add_argument('--root', default = None,
                   help = "the variant's directory (default results/oecd26_theta_<variant>); writes calibration/")
    p.add_argument('--main', action = 'store_true',
                   help = "with --theta config.US['theta']: results/calibration/ (the paper's), oecd26.thetaRoot")
    p.add_argument('--rhos', type = float, nargs = '*', default = None,
                   help = 'the grid as a list of rho (overrides --lo/--hi/--step); must contain --anchor')
    p.add_argument('--out', default = None)
    p.add_argument('--pkldir', default = None, help = 'where the per-point pickled instances go')
    a = p.parse_args()
    a.out, a.pkldir, a.root = [None if x is None else os.path.abspath(os.path.join(CWD, x))
                               for x in (a.out, a.pkldir, a.root)]

    import oecd26
    if a.theta in (None, oecd26.INDEX):
        if a.root is not None or a.main:
            raise SystemExit('--root and --main need --theta VARIANT')
        out = a.out or os.path.join(OUTDIR, 'US' + ('_rhoGridCommonX' if a.commonX else '_rhoGrid') + '.csv')
        pkldir = a.pkldir or os.path.join(OUTDIR, 'instancesUS' + ('CommonX' if a.commonX else ''))
        a.theta = oecd26.indexTheta(a.theta, [out, pkldir], 'calibrateRhoGrid.py')
        pars, columns, extra = dict(testmod.pars), COLUMNS, {}
        os.makedirs(OUTDIR, exist_ok = True)
    else:
        θ0, out, pkldir = thetaRun(a.theta, a.commonX, a.out, a.pkldir, a.root, a.main)
        pars, columns, extra = dict(testmod.pars) | {'θ0': θ0}, COLUMNS + THETACOLUMNS, {'thetaSource': a.theta, 'θ0': θ0}
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok = True)
    os.makedirs(pkldir, exist_ok = True)

    grid = rhoGrid(a)
    if round(a.anchor, 6) not in set(grid):
        raise SystemExit('the anchor rho={} is not on the grid {}'.format(a.anchor, list(grid)))
    commit = gitCommit()
    # Keyed by solver: 'ns' sizes a state grid the LOG solver does not have, and passing it there would
    # be silently ignored rather than rejected. 'n' is shared -- both solvers grid tau. smoothKnots and
    # interpKind go to BOTH: they are well-posedness settings, and leaving the LOG anchor on different
    # ones would solve rho=1 under a different residual from every other point -- and since the anchor's
    # x is the warm start the whole march is seeded from, that difference would propagate.
    wp = {'smoothKnots': a.smoothKnots if a.smoothKnots else None, 'interpKind': a.interpKind}
    gridSettings = {'CRRA': {'n': a.n, 'ns': a.ns} | wp, 'LOG': {'n': a.n} | wp}
    verify = ({'CRRA': {'n': a.verifyN, 'ns': a.verify} | wp, 'LOG': {'n': a.verifyN} | wp}
              if a.verify else None)

    done, df = ({}, pd.DataFrame(columns = COLUMNS)) if a.force else readDone(out)
    if extra and len(df):                 # a resumed variant csv must be at this θ0 (#13)
        same = (df['thetaSource'] == a.theta) & (pd.to_numeric(df['θ0']) - θ0).abs().le(1e-12) \
            if set(THETACOLUMNS) <= set(df.columns) else pd.Series(False, index = df.index)
        if not same.all():
            raise SystemExit('{}: rows at another thetaSource or θ0 than {} {!r}: rerun with --force'.format(
                out, a.theta, θ0))
    if not extra and len(df) and 'thetaSource' in df.columns and df['thetaSource'].notna().any():
        raise SystemExit('{}: rows at thetaSource {}, not the index θ: rerun with --force or another --out'.format(
            out, sorted(set(df['thetaSource'].dropna().astype(str)))))
    if not a.force:
        resumeAt(done, df, rowSettings(a), grid)
    # Keyed by rho, not a list: marchGrid reports one record per *attempt*, so a point that fails on the
    # extrapolated start and succeeds on the carried one arrives twice, and a resumed run revisits values
    # already in the csv. Keying makes every write idempotent and lets a success overwrite a failure.
    rows = {} if a.force else {round(float(r['ρ']), 6): r for r in df.to_dict('records')}

    print('grid: {} points, {}'.format(len(grid), ', '.join('{:g}'.format(v) for v in grid)))
    print('variant: {},  anchor rho={},  tau nodes {},  state nodes {},  verify at {}/{}'.format(
        'common X' if a.commonX else 'vector X_i', a.anchor, a.n, a.ns, a.verifyN, a.verify or 'off'))
    print('targets: tau0={:.6f}  RR0={:.6f}  Ushare0={:.6f}'.format(pars['τ0'], pars['RR0'], pars['Ushare0']))
    if extra:
        print('theta imposed: θ0={!r} (variant {!r} of {}; RR0 not read)'.format(θ0, a.theta, oecd26.FILES['theta']))
    print('out: {}'.format(os.path.relpath(out, REPO)))
    if done:
        print('resuming: {} of {} points already present'.format(
            sum(round(float(v), 6) in done for v in grid), len(grid)))

    m = ModelUS(pars = pars | {'ρ': float(a.anchor)}, commonX = a.commonX, **testmod.kwargs)

    def write():
        (pd.DataFrame(list(rows.values())).reindex(columns = columns).sort_values('ρ')
           .to_csv(out, index = False))

    def solve(ρ, x0):
        """ marchGrid's callback. A point already in the csv is returned from it -- not re-solved, but
        still installed into db, so the next point's *inner* warm start is not left at a stale parameter
        point while its outer one is up to date.

        Under --commonX the reinstall covers (β,ω) only, not X: X is not part of the search vector, and it
        enters no aggregate (docs eq:calibration:yxCommonX), so the next point's warm start cannot see it.
        A cached point keeps whatever X its own row already records. """
        key = round(float(ρ), 6)
        if key in done:
            m.db.update(m.adjPar('ρ', float(ρ)))
            m._calSetPars(m._calFromX(done[key]))
            print('  rho={:<6} cached'.format(key))
            return {'x': done[key], 'cached': True}
        rec = m.calibratePoint(ρ, x0 = x0, gridSettings = gridSettings, verify = verify)
        print('  rho={:<6} {:<4} max|res|={:.2e}  verify={:.2e}  nfev={:<3} {:5.1f}s  '
              'β={:.5f} ω={:.5f}{}  sr={:.4f}'.format(
                  key, rec['preferences'], rec['residual'], rec.get('verifyResidual', np.nan),
                  rec['nfev'], rec['time'], rec['β'], rec['ω'],
                  ' X={:.4f}'.format(rec['X']) if 'X' in rec else '', rec['sr']))
        return rec

    def onPoint(r):
        """ Persist after every successful point. Failures are only *reported* here -- a failed attempt is
        not a failed point, since the march still has the carried start and step-halving to try. """
        key = round(float(r['value']), 6)
        if not r['ok']:
            print('  rho={:<6} attempt failed ({}): {}'.format(key, r['x0Source'], r['error'][:100]))
            return
        if r['result'].get('cached'):
            return
        rows[key] = toRow(r['result'], r['requested'], commit) | extra
        write()
        with open(os.path.join(pkldir, 'rho_{:.4f}.pkl'.format(r['value'])), 'wb') as f:
            pickle.dump(m, f)

    # marchGrid directly rather than m.calibrateGrid: the resume-from-csv and progress-printing live in
    # `solve` above, and calibrateGrid would substitute its own callback and silently discard both.
    res = continuation.marchGrid(grid, solve, anchor = a.anchor, degree = a.degree,
                                 maxHalvings = a.maxHalvings, onPoint = onPoint)

    # A failed *attempt* is not a failed point -- the retry ladder and step-halving may have rescued it.
    # Only values with no successful record anywhere are written as failures.
    solvedValues = {round(float(r['value']), 6) for r in res['records'] if r['ok']}
    hardFails = [r for r in res['failures'] if round(float(r['value']), 6) not in solvedValues]
    for r in hardFails:
        rows[round(float(r['value']), 6)] = {
            'ρ': round(float(r['value']), 6), 'requested': r['requested'], 'commit': commit,
            'residual': np.nan, 'preferences': 'FAILED: ' + r['error'][:120],
            'timestamp': datetime.datetime.now().isoformat(timespec = 'seconds')} | extra
    write()

    ok = len(grid)-len(hardFails)
    print('\n{} of {} points solved -> {}'.format(ok, len(grid), os.path.relpath(out, REPO)))
    if hardFails:
        print('failed: ' + ', '.join('{:g}'.format(r['value']) for r in hardFails))
    return 1 if hardFails else 0


if __name__ == '__main__':
    sys.exit(main())

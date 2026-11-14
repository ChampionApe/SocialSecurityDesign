r""" Calibrate the France / UK models (modelFR.ModelFR) across a grid of CRRA parameter values.

Run:  .venv\Scripts\python.exe python\US\calibrateRhoGridEU.py --country FR
      ... --country UK                        the UK on its own income groups
      ... --commonX                           the common-X variant (its own csv and pickle directory)
      ... --force                             re-solve points already in the csv
      ... --rhos 0.5 1 2                      the grid as a list (calibrateRhoGrid.rhoGrid)
      ... --theta VARIANT [--root DIR | --main]   theta imposed as db['θ0'] = the country's row of
                                              data/oecd26_theta.csv under VARIANT (oecd26.countryθ0: FRA, GBR),
                                              against the U.S. sweep of the same variant (calibrateRhoGrid.py
                                              --theta, at the variant's U.S. θ0 at every rho, oecd26.usBeta);
                                              csv and pickles under DIR/calibration/ (DIR the variant's directory,
                                              oecd26.thetaDir), results/calibration/ under --main (the pipeline's
                                              main variant, oecd26.thetaRoot); rows carry thetaSource and θ0.
      ... --theta index                       the workbook's θ; without --theta that run is refused into results/
                                              (oecd26.indexTheta): use --out and --pkldir.

Same march as calibrateRhoGrid.py (read that first -- the anchor, the resume-from-csv contract, the
pickle-is-a-cache rule and the note on how to read verifyResidual all carry over unchanged). Three
differences, all of them consequences of the protocol in modelFR.py:

  * The search is ONE-dimensional. beta is imposed from the US calibration at the same rho, so only omega
    is solved for, against tau_{t0}. The csv therefore has a single x0 column, not x0/x1, and beta is a
    RECORDED column rather than a solved one -- worth reading down the column all the same, since it must
    reproduce the US sweep's beta exactly at every point.

  * Every visited rho must exist in the US sweep csv. usRefFromCsv matches rows exactly and refuses to
    interpolate: an interpolated beta is not the US calibration's beta at that rho, and quietly importing
    one would make the whole point of "impose the US discount factor" false. The grid is checked against
    the reference csv BEFORE the march starts, so a missing value fails in a second rather than halfway
    through a sweep.

    This is also why --maxHalvings defaults to 0 here. marchGrid's step-halving inserts intermediate
    values (0.55, 0.65, ...) and KEEPS them in the record; with no US row at those values they cannot be
    solved. If a point resists, the fix is to run the US sweep on the finer grid first, not to soften the
    reference.

  * A resumed row must record this run's settings (calibrateRhoGrid.resumeAt, finding #13): the grid, and in the
    columns `retAge` (France's or the UK's own ν correction, testEU.retAgeHost) and `βUS` (the U.S. reference β at
    its rho) the two inputs from outside the workbook; rows written before those two columns are resumed unchecked.

  * Two extra columns. `lambda` is the proportional X_i rescaling that puts average hours on the
    US-referenced target -- under vector X it is the second calibrated quantity, recovered in closed form
    after the root (docs eq:us:model:scaleInvariance). It is absent under --commonX, where X plays that
    role instead. `hoursDrift` is how far tau/R/hbar moved when the rescaled model was re-solved; see
    ModelFR.hoursDriftTol for what level to expect and why it is recorded rather than asserted away.

A relative --out, --pkldir or --root is the caller's working directory's (the script itself runs in python/US).

What must be constant down the csv: tau (the one target) and beta (imposed, and equal to the US sweep's).
What moves: omega, lambda, the untargeted savings rate, and R -- which is a prediction here, not a target,
and is the most interesting column in the file for exactly that reason.
"""
import os, sys, argparse, pickle, subprocess, datetime
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                 # where a relative --out, --pkldir or --root points: the chdir below moves to HERE

# See calibrateRhoGrid.py: line buffering so progress appears as it happens, UTF-8 so the Greek parameter
# names in every progress line do not take the march down on a redirected Windows stdout.
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)                                        # testEU.py reads data/ relative to its own location

import testEU
from gridsearch import continuation
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'calibration')
THETACOLUMNS = ['thetaSource', 'θ0']      # appended under --theta only

# Ordered as the US csv is: what was solved, what it cost, whether to believe it, then the answer.
COLUMNS = ['ρ', 'preferences', 'requested', 'residual', 'verifyResidual', 'hoursDrift',
           'β', 'ω', 'λ', 'X',
           'R', 'τ', 'sr', 'h', 'hbar', 'nfev', 'time', 'n', 'ns', 'smoothKnots', 'interpKind', 'retAge', 'βUS',
           'nEqMax', 'nCandMax', 'nFallback',
           'x0', 'commit', 'timestamp']
INPUTCOLUMNS = ('retAge', 'βUS')          # recorded since 2026-11-04; a row without them is resumed unchecked


def gitCommit():
    try:
        return subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd = REPO, capture_output = True,
                              text = True, timeout = 10).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def toRow(rec, requested, commit):
    """ calibratePoint's record -> one flat csv row, with the tax counts of the PEE solve at the calibrated
    point (the record's 'multiplicity' as nEqMax, nCandMax, nFallback; policy.multiplicityColumns). """
    g = rec['gridSettings']
    row = {k: rec.get(k) for k in ('ρ', 'preferences', 'residual', 'verifyResidual', 'hoursDrift',
                                   'β', 'ω', 'λ', 'X', 'R', 'τ', 'sr', 'h', 'hbar', 'nfev', 'time')}
    row |= {'requested': requested, 'n': g.get('n'), 'ns': g.get('ns'),
            'smoothKnots': g.get('smoothKnots'), 'interpKind': g.get('interpKind'), 'commit': commit,
            'timestamp': datetime.datetime.now().isoformat(timespec = 'seconds')}
    row |= {f'x{i}': v for i, v in enumerate(rec['x'])}
    return row | multiplicityColumns(rec.get('multiplicity'))


def readDone(path):
    """ {rho: x} for the points already solved, as the resume cache. One-element x -- see the docstring. """
    if not os.path.exists(path):
        return {}, pd.DataFrame(columns = COLUMNS)
    df = pd.read_csv(path)
    done = {round(float(r['ρ']), 6): np.array([r['x0']], dtype = float)
            for _, r in df.iterrows() if np.isfinite(r.get('residual', np.nan))}
    return done, df


def thetaRun(variant, country, commonX, out = None, pkldir = None, root = None, main = False, rhos = ()):
    """ --theta VARIANT: (θ0 of `country`, csv path, pickle directory, the U.S. sweep it is referenced to), all
    under root/calibration/ (oecd26.thetaRoot(variant, main) unless `root`, a relative path the repo's). The U.S.
    sweep's rows at every rho of `rhos` must be at the variant's U.S. θ0 (oecd26.usBeta). SystemExit for --main with
    --root, a variant without the country's or the U.S.'s θ, a U.S. sweep missing or at another θ0, or, without
    --main, a path under results/calibration/. """
    import oecd26
    if main and root is not None:
        raise SystemExit('--main writes results/calibration/: not with --root')
    try:
        θ0 = oecd26.countryθ0(country, variant)
        top = oecd26.thetaRoot(variant, main) if root is None else os.path.join(REPO, root)
        cal = os.path.join(top, 'calibration')
        usCsv = os.path.join(cal, 'US_rhoGrid' + ('CommonX' if commonX else '') + '.csv')
        for ρ in rhos:
            oecd26.usBeta(usCsv, float(ρ), thetaSource = variant)
    except (ValueError, FileNotFoundError, KeyError) as e:          # NoTheta is a ValueError
        raise SystemExit('--theta {}: {}'.format(variant, e))
    out = out or os.path.join(cal, country + '_rhoGrid' + ('CommonX' if commonX else '') + '.csv')
    pkldir = pkldir or os.path.join(cal, 'instances' + country + ('CommonX' if commonX else ''))
    for q in (out, pkldir):
        if not main and os.path.normcase(os.path.abspath(q)).startswith(os.path.normcase(OUTDIR + os.sep)):
            raise SystemExit('--theta {}: {} is under results/calibration/, which holds the main θ variant '
                             '(--main) only'.format(variant, q))
    return θ0, out, pkldir, usCsv


def main():
    p = argparse.ArgumentParser(description = 'Calibrate the France/UK model across a grid of rho.')
    p.add_argument('--country', default = 'FR', choices = ('FR', 'UK'))
    p.add_argument('--lo', type = float, default = 0.5)
    p.add_argument('--hi', type = float, default = 2.0)
    p.add_argument('--step', type = float, default = 0.1)
    p.add_argument('--anchor', type = float, default = 1.0)
    p.add_argument('--commonX', action = 'store_true')
    p.add_argument('--n', type = int, default = 101, help = 'nodes on the tau grid')
    p.add_argument('--ns', type = int, default = 150, help = 'nodes on the CRRA state grid (LOG ignores it)')
    p.add_argument('--verify', type = int, default = 225,
                   help = 'refined ns for the resolution check; 0 to skip')
    p.add_argument('--verifyN', type = int, default = 151)
    p.add_argument('--smoothKnots', type = int, default = 4)
    p.add_argument('--interpKind', default = 'linear', choices = ('linear', 'cubic', 'pchip'))
    p.add_argument('--degree', type = int, default = 1)
    # 0, not the US sweep's 2: an inserted intermediate rho has no row in the US reference csv. See the
    # module docstring.
    p.add_argument('--maxHalvings', type = int, default = 0)
    p.add_argument('--force', action = 'store_true')
    p.add_argument('--out', default = None)
    p.add_argument('--pkldir', default = None)
    p.add_argument('--rhos', type = float, nargs = '*', default = None,
                   help = 'the grid as a list of rho (overrides --lo/--hi/--step); must contain --anchor')
    p.add_argument('--theta', default = None,
                   help = "a variant of data/oecd26_theta.csv: θ imposed at the country's row (db['θ0']); 'index' the "
                          "workbook's θ (required for that into results/)")
    p.add_argument('--root', default = None, help = "the variant's directory (default results/oecd26_theta_<variant>)")
    p.add_argument('--main', action = 'store_true',
                   help = "with --theta config.US['theta']: results/calibration/ (the paper's)")
    a = p.parse_args()
    a.out, a.pkldir, a.root = [None if x is None else os.path.abspath(os.path.join(CWD, x))
                               for x in (a.out, a.pkldir, a.root)]

    from calibrateRhoGrid import rhoGrid, rowSettings, resumeAt
    import oecd26
    grid = rhoGrid(a)
    if round(a.anchor, 6) not in set(grid):
        raise SystemExit('the anchor rho={} is not on the grid {}'.format(a.anchor, list(grid)))
    if a.theta in (None, oecd26.INDEX):
        if a.root is not None or a.main:
            raise SystemExit('--root and --main need --theta VARIANT')
        tag = a.country + '_rhoGrid' + ('CommonX' if a.commonX else '')
        out = a.out or os.path.join(OUTDIR, tag + '.csv')
        pkldir = a.pkldir or os.path.join(OUTDIR, 'instances' + a.country + ('CommonX' if a.commonX else ''))
        a.theta = oecd26.indexTheta(a.theta, [out, pkldir], 'calibrateRhoGridEU.py')
        θ0, usCsv, columns, extra = None, None, COLUMNS, {}
        os.makedirs(OUTDIR, exist_ok = True)
    else:
        θ0, out, pkldir, usCsv = thetaRun(a.theta, a.country, a.commonX, a.out, a.pkldir, a.root, a.main, grid)
        columns, extra = COLUMNS + THETACOLUMNS, {'thetaSource': a.theta, 'θ0': θ0}
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok = True)
    os.makedirs(pkldir, exist_ok = True)

    commit = gitCommit()
    wp = {'smoothKnots': a.smoothKnots if a.smoothKnots else None, 'interpKind': a.interpKind}
    gridSettings = {'CRRA': {'n': a.n, 'ns': a.ns} | wp, 'LOG': {'n': a.n} | wp}
    verify = ({'CRRA': {'n': a.verifyN, 'ns': a.verify} | wp, 'LOG': {'n': a.verifyN} | wp}
              if a.verify else None)

    # Fail before the march, not during it: every rho on the grid needs a row in the US sweep csv.
    ref = testEU.usReference(commonX = a.commonX, csv = usCsv)
    missing = []
    for ρ in grid:
        try:
            ref(float(ρ))
        except KeyError:
            missing.append(float(ρ))
    if missing:
        print('The US reference sweep has no row at: ' + ', '.join('{:g}'.format(v) for v in missing))
        print('Run python/US/calibrateRhoGrid.py{} over those values first -- the reference is matched '
              'exactly and is deliberately not interpolated.'.format(' --commonX' if a.commonX else ''))
        return 2

    done, df = ({}, pd.DataFrame(columns = COLUMNS)) if a.force else readDone(out)
    if extra and len(df):                 # a resumed variant csv must be at this θ0 (#13)
        same = ((df['thetaSource'] == a.theta) & (pd.to_numeric(df['θ0']) - θ0).abs().le(1e-12)
                if set(THETACOLUMNS) <= set(df.columns) else pd.Series(False, index = df.index))
        if not same.all():
            raise SystemExit('{}: rows at another thetaSource or θ0 than {} {!r}: rerun with --force'.format(
                out, a.theta, θ0))
    if not extra and len(df) and 'thetaSource' in df.columns and df['thetaSource'].notna().any():
        raise SystemExit('{}: rows at thetaSource {}, not the index θ: rerun with --force or another --out'.format(
            out, sorted(set(df['thetaSource'].dropna().astype(str)))))
    retAge = testEU.retAgeHost(a.country)
    βUS = lambda ρ: float(ref(float(ρ))['β'])
    if not a.force:
        resumeAt(done, df, rowSettings(a) | {'retAge': retAge, 'βUS': βUS}, grid, INPUTCOLUMNS)
    rows = {} if a.force else {round(float(r['ρ']), 6): r for r in df.to_dict('records')}

    print('country: {},  variant: {}'.format(a.country, 'common X' if a.commonX else 'vector X_i'))
    print('grid: {} points, {}'.format(len(grid), ', '.join('{:g}'.format(v) for v in grid)))
    if extra:
        print('theta imposed: θ0={!r} (variant {!r}; RR0 not read), U.S. reference {}'.format(
            θ0, a.theta, os.path.relpath(usCsv, REPO)))
    print('anchor rho={},  tau nodes {},  state nodes {},  verify at {}/{}'.format(
        a.anchor, a.n, a.ns, a.verifyN, a.verify or 'off'))
    print('out: {}'.format(os.path.relpath(out, REPO)))
    if done:
        print('resuming: {} of {} points already present'.format(
            sum(round(float(v), 6) in done for v in grid), len(grid)))

    m = testEU.model(a.country, commonX = a.commonX, ρ = float(a.anchor), usCsv = usCsv, θ0 = θ0)

    def write():
        (pd.DataFrame(list(rows.values())).reindex(columns = columns).sort_values('ρ')
           .to_csv(out, index = False))

    def solve(ρ, x0):
        """ marchGrid's callback. As the US sweep's: a cached point is reinstalled into db rather than
        re-solved, so the next point's inner warm start moves with its outer one. setUSRef is called on
        the cached path too -- beta is a function of rho here, and leaving the previous point's beta in
        db would warm-start the next solve from the wrong model. """
        key = round(float(ρ), 6)
        if key in done:
            m.db.update(m.adjPar('ρ', float(ρ)))
            m.setUSRef()
            m._calSetPars(m._calFromX(done[key]))
            print('  rho={:<6} cached'.format(key))
            return {'x': done[key], 'cached': True}
        rec = m.calibratePoint(ρ, x0 = x0, gridSettings = gridSettings, verify = verify)
        print('  rho={:<6} {:<4} max|res|={:.2e}  verify={:.2e}  drift={:.1e}  nfev={:<3} {:5.1f}s  '
              'β={:.5f} ω={:.5f} {}  R={:.4f} sr={:.4f}'.format(
                  key, rec['preferences'], rec['residual'], rec.get('verifyResidual', np.nan),
                  rec['hoursDrift'], rec['nfev'], rec['time'], rec['β'], rec['ω'],
                  'X={:.4f}'.format(rec['X']) if 'X' in rec else 'λ={:.5f}'.format(rec['λ']),
                  rec['R'], rec['sr']))
        return rec

    def onPoint(r):
        key = round(float(r['value']), 6)
        if not r['ok']:
            print('  rho={:<6} attempt failed ({}): {}'.format(key, r['x0Source'], r['error'][:100]))
            return
        if r['result'].get('cached'):
            return
        rows[key] = toRow(r['result'], r['requested'], commit) | {'retAge': retAge, 'βUS': βUS(key)} | extra
        write()
        with open(os.path.join(pkldir, 'rho_{:.4f}.pkl'.format(r['value'])), 'wb') as f:
            pickle.dump(m, f)

    res = continuation.marchGrid(grid, solve, anchor = a.anchor, degree = a.degree,
                                 maxHalvings = a.maxHalvings, onPoint = onPoint)

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

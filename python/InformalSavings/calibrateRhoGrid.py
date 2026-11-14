r""" Calibrate the informalSavings model across a grid of CRRA parameter values.

Run:  .venv\Scripts\python.exe python\InformalSavings\calibrateRhoGrid.py [options]
      ... --lo 0.5 --hi 2.0 --step 0.1        the default sweep
      ... --lo 0.5 --hi 2.0 --step 0.25       a coarse pass first
      ... --force                             re-solve points already in the csv

The march is anchored at rho=1, where the LOG solver applies and no warm start is needed, and walks
outward in both directions (gridsearch.continuation.marchGrid). Each calibration costs ~26 PEE solves at
~21 s each, so the sweep is measured in hours and is built to survive being interrupted:

  * results/calibration/informalSavings_rhoGrid.csv is the record, rewritten after *every* point. A run
    started over reads it back, and any rho already present is returned from it without re-solving --
    which also means those points still seed the extrapolation, so a resumed march is no worse warm-started
    than an uninterrupted one.
  * results/calibration/instances/*.pkl is a cache, not the record: one pickled model per point (~62 kB),
    holding the converged db. Reconstructing the same thing from the csv costs one _calSetPars and one
    solve, so nothing is lost if these are deleted -- they exist to make a solved point cheap to open
    interactively. Full solutions (policy functions) are deliberately NOT stored: CRRA's are not
    pickleable at all (gridsearch.interp's griddedInterp2D returns a closure) and they regenerate in ~21 s.

The inner grid defaults to nι=ns=45, which is NOT the PEE solve's default and must not be reduced -- the
CRRA calibration needs a finer inner grid than the CRRA solve (finding #3). Whether 45 is still enough
far from rho=1 has never been tested, which is what --verify exists to answer: it re-evaluates the outer
residual at the converged parameters on a finer grid and stores it beside the point. A point whose
verifyResidual is not small is converged but not resolved, and should not be read as a result.
"""
import os, sys, argparse, pickle, subprocess, datetime
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))

# A point takes minutes and the sweep takes hours, so progress has to appear as it happens rather than
# when the process exits. Whenever stdout is a pipe or a log file it is block-buffered by default, which
# hides the entire run; set here rather than relying on the caller remembering `python -u`.
# encoding too, not only buffering: every progress line carries Greek parameter names, and when
# stdout is a pipe or a log file it defaults to the ANSI codepage on Windows -- which raises
# UnicodeEncodeError on the first point and takes the whole march down with it.
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)                                        # relative paths resolve here; test.py finds data/ itself

import test as testmod
from model import ModelInformalSavings
from gridsearch import continuation

OUTDIR = os.path.join(REPO, 'results', 'calibration')
CSV    = os.path.join(OUTDIR, 'informalSavings_rhoGrid.csv')
PKLDIR = os.path.join(OUTDIR, 'instances')

# Ordered so the csv reads as: what was solved, what it cost, whether to believe it, then the answer.
COLUMNS = ['ρ', 'preferences', 'requested', 'residual', 'verifyResidual', 'occupancyι', 'occupancys',
           'nEqMax', 'nCandMax', 'nFallback',
           'β', 'ω', 'η0', 'X0', 'X', 'hbar', 'commonX', 'epsPre0',
           'KY', 'sr', 'τ', 'ι', 'yInf', 'spendGDP', 'srGDP',
           'nRoots', 'nfev', 'time', 'nι', 'ns', 'nτ', 'x0', 'x1', 'x2', 'x3',
           'commit', 'timestamp']
# yInf/spendGDP/srGDP: informal output over GDP, pension spending over GDP (the datum the tax target is
# derived from, model.taxTarget) and savings over GDP at t0, GDP = formal plus informal output.
# occupancy* sit beside verifyResidual because they answer the same kind of question and neither is
# asserted: verifyResidual catches a point that is converged but not resolved, occupancy catches a state
# grid whose nodes are mostly somewhere the dynamics never go. A grid that is too NARROW announces itself
# (states outside it are reported infeasible, never clipped); one that is too WIDE is silent, and this is
# the column that makes it visible. toRow reindexes to this list, so a record key absent here is dropped.
# nEqMax/nCandMax/nFallback are the tax candidates' counts of the point's final solve (num_robustroot.tex):
# nEqMax > 1 means the equilibrium selection bound somewhere, nFallback > 0 that some state had no
# equilibrium on the grid. Non-key, like every column but ρ (finding #13); -1 on a failed point.


def gitCommit():
    """ Short hash of the code that produced a row -- the sweep outlives several edits of the solver. """
    try:
        return subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd = REPO, capture_output = True,
                              text = True, timeout = 10).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def toRow(rec, requested, commit):
    """ calibratePoint's record -> one flat csv row. """
    g = rec['gridSettings']
    row = {k: rec.get(k) for k in ('ρ', 'preferences', 'residual', 'verifyResidual', 'occupancyι',
                                   'occupancys', 'β', 'ω', 'η0', 'X0', 'X', 'hbar', 'commonX', 'epsPre0',
                                   'KY', 'sr', 'τ', 'ι', 'yInf', 'spendGDP', 'srGDP', 'nRoots', 'nfev',
                                   'time')}
    row |= {k: rec.get(k, -1) for k in ('nEqMax', 'nCandMax', 'nFallback')}
    row |= {'requested': requested, 'nι': g.get('nι'), 'ns': g.get('ns'), 'nτ': g.get('n'),
            'commit': commit, 'timestamp': datetime.datetime.now().isoformat(timespec = 'seconds')}
    row |= {f'x{i}': v for i, v in enumerate(rec['x'])}
    return row


def readDone(path):
    """ {rho: x} for the points already solved, as the resume cache. Keyed on the rounded value so a
    grid rebuilt from the same --lo/--step matches it exactly rather than by floating-point luck. """
    if not os.path.exists(path):
        return {}, pd.DataFrame(columns = COLUMNS)
    df = pd.read_csv(path)
    done = {round(float(r['ρ']), 6): np.array([r[f'x{i}'] for i in range(4)], dtype = float)
            for _, r in df.iterrows() if np.isfinite(r.get('residual', np.nan))}
    return done, df


def main():
    p = argparse.ArgumentParser(description = 'Calibrate informalSavings across a grid of rho.')
    p.add_argument('--lo', type = float, default = 0.5)
    p.add_argument('--hi', type = float, default = 2.0)
    p.add_argument('--step', type = float, default = 0.1)
    p.add_argument('--anchor', type = float, default = 1.0)
    p.add_argument('--nι', '--niota', dest = 'nι', type = int, default = 45)
    p.add_argument('--ns', type = int, default = 45)
    p.add_argument('--verify', type = int, default = 60, help = 'refined nι=ns for the resolution check; 0 to skip')
    p.add_argument('--verifyLOG', type = int, default = 75,
                   help = "refined nι for the LOG anchor's resolution check (its working nι is 50)")
    # Default 4, policy.py's own. It must NOT be None here: initGS merges this dict over the class
    # defaults, so a None would select the adaptive smoother rather than request the default. 0 selects
    # the adaptive smoother.
    p.add_argument('--smoothKnots', type = int, default = 4,
                   help = 'fixed-knot policy smoother: interior knot every m-th node. 0 selects the '
                          'adaptive smoother, which re-chooses its knot COUNT from the data and makes the '
                          'outer residual discontinuous -- only for reproducing pre-2026-08-19 results.')
    p.add_argument('--interpKind', default = 'cubic', choices = ('linear', 'cubic', 'pchip'),
                   help = "continuation interpolant, BOTH solvers; 'linear' does not converge away "
                          "from rho=1 under CRRA, and does not converge in nι at all under LOG")
    p.add_argument('--degree', type = int, default = 1, help = 'extrapolation degree for the warm start')
    p.add_argument('--maxHalvings', type = int, default = 2)
    # interpKind is shared, so crossing rho=1 changes the solver, the state dimension and nι (50->45) but
    # not the interpolant. --common additionally matches the grid SIZES, to isolate the two recursions
    # from every setting at once (a diagnostic: the LOG answer is converged in nι at 'cubic').
    p.add_argument('--common', action = 'store_true',
                   help = 'also give the LOG anchor the CRRA grid SIZES, so nothing but the recursion '
                          'changes at rho=1 (diagnostic; interpKind is already shared)')
    p.add_argument('--force', action = 'store_true', help = 're-solve points already present in the csv')
    # The common-X calibration variant (model.py __init__): its own csv and instance directory, since
    # the two variants answer the same rho with different (eta_i, X_i, eta_0, X_0) (finding #13).
    p.add_argument('--commonX', action = 'store_true',
                   help = 'the common-X variant: informalSavings_rhoGridCommonX.csv, instancesCommonX/')
    # The pre-reform eps = 0 variant (model.py __init__), likewise its own csv and instance directory; it
    # composes with --commonX (the tags concatenate: ...CommonXEpsPre0).
    p.add_argument('--epsPre0', action = 'store_true',
                   help = 'pre-reform eps = 0: informalSavings_rhoGrid[CommonX]EpsPre0.csv, '
                          'instances[CommonX]EpsPre0/')
    p.add_argument('--out', default = None, help = "default: the variant's own csv")
    # A sweep written to a non-default --out must also get its own --pkldir, or it silently overwrites the
    # canonical sweep's instances wherever the two grids share a value (the filename is the rho alone).
    # That matters even when the settings agree: the code is bitwise reproducible within a process but not
    # across them, so the overwrite is not a no-op and shockUniversal.py reads these by name.
    p.add_argument('--pkldir', default = None, help = "where the per-point pickled instances go (default: the variant's own)")
    # The anchor's starting parameters. Every other point is seeded from the march history, but the anchor
    # has only test.py's defaults (beta=0.6, omega=2), and the iota state grid is built from the steady
    # state at whatever parameters the residual is FIRST evaluated at -- at a capital share of 0.35 those
    # defaults put the informal steady state outside the net-saver region and the grid is degenerate
    # before the root takes a single step. Give the anchor a guess near its solution instead.
    p.add_argument('--x0', type = float, nargs = 4, metavar = ('BETA', 'OMEGA', 'ETA0', 'X0'),
                   default = None, help = "the anchor's starting (beta, omega, eta0, X0)")
    args = p.parse_args()
    tag = ('CommonX' if args.commonX else '') + ('EpsPre0' if args.epsPre0 else '')
    if args.out is None:
        args.out = os.path.join(OUTDIR, f'informalSavings_rhoGrid{tag}.csv')
    if args.pkldir is None:
        args.pkldir = os.path.join(OUTDIR, f'instances{tag}')

    pkldir = args.pkldir
    os.makedirs(OUTDIR, exist_ok = True)
    os.makedirs(pkldir, exist_ok = True)
    grid = np.round(np.arange(args.lo, args.hi + 0.5*args.step, args.step), 6)
    commit = gitCommit()
    # Keyed by solver. The grid SIZES are the CRRA calibration's alone (the LOG anchor keeps its own nι=50):
    # they are resolution choices and legitimately differ. interpKind and smoothKnots go to BOTH solvers:
    # they are well-posedness choices (findings #4, #5, #7), and the anchor's x seeds the whole march, so
    # solving rho=1 under a different interpolant or smoother would propagate to every point.
    # interpKind='cubic' is required rather than preferred: at 'linear' the CRRA calibration away from rho=1
    # does not converge and the LOG solve does not converge in nι (test_calibrationGrid.py's docstring;
    # archive/notes/informalSavings_resolvedIssues.md).
    knots = args.smoothKnots if args.smoothKnots else None      # 0 -> the adaptive branch
    crraGS = {'nι': args.nι, 'ns': args.ns, 'interpKind': args.interpKind, 'smoothKnots': knots}
    logGS = {'smoothKnots': knots, 'interpKind': args.interpKind}
    gridSettings = {'CRRA': crraGS, 'LOG': dict(crraGS) if args.common else logGS}
    # The resolution check runs at both solvers; 75 is LOG's refinement rung (1.5x its nι=50), matching
    # test_calibrationGrid.py's VERIFY.
    verify = ({'CRRA': {'nι': args.verify, 'ns': args.verify, 'interpKind': args.interpKind},
               'LOG': {'nι': args.verifyLOG}}
              if args.verify else None)

    done, df = ({}, pd.DataFrame(columns = COLUMNS)) if args.force else readDone(args.out)
    # Keyed by rho, not a list: marchGrid reports one record per *attempt*, so a point that fails on the
    # extrapolated start and succeeds on the carried one arrives twice, and a resumed run revisits values
    # already in the csv. Keying makes every write idempotent and lets a success overwrite an earlier
    # failure at the same rho.
    rows = {} if args.force else {round(float(r['ρ']), 6): r for r in df.to_dict('records')}
    print('grid: {} points, {} to {} step {}'.format(len(grid), grid[0], grid[-1], args.step))
    print('anchor rho={}, inner grid {}x{}, verify at {}/{}, smoothKnots={}, interpKind={} (both '
          'solvers), LOG grid {}'.format(
              args.anchor, args.nι, args.ns, args.verify or 'off', args.verifyLOG,
              knots or 'adaptive', args.interpKind, 'CRRA sizes' if args.common else 'own nι=50'))
    if done:
        print('resuming: {} of {} points already in {}'.format(
            sum(round(float(v), 6) in done for v in grid), len(grid), os.path.relpath(args.out, REPO)))

    m = ModelInformalSavings(pars = testmod.pars | {'ρ': float(args.anchor)}, commonX = args.commonX,
                             epsPre0 = args.epsPre0, **testmod.kwargs)
    print('calibration variant: {}, pre-reform eps {}'.format(
        'common X (formal workweek targeted)' if args.commonX else 'vector X_i (relative hours as data)',
        '= 0' if args.epsPre0 else 'the basic pension (getEps)'))

    def write():
        (pd.DataFrame(list(rows.values())).reindex(columns = COLUMNS).sort_values('ρ')
           .to_csv(args.out, index = False))

    def solve(ρ, x0):
        """ marchGrid's callback. A point already in the csv is returned from it -- not re-solved, but
        still installed into db, so the next point's *inner* warm start is not left at a stale parameter
        point while its outer one is up to date. """
        key = round(float(ρ), 6)
        if key in done:
            x = done[key]
            m.db.update(m.adjPar('ρ', float(ρ)))
            m._calSetPars(m._calFromX(x))
            print('  rho={:<6} cached'.format(key))
            return {'x': x, 'cached': True}
        rec = m.calibratePoint(ρ, x0 = x0, gridSettings = gridSettings, verify = verify)
        print('  rho={:<6} {:<4} max|res|={:.2e}  verify={:.2e}  nfev={:<3} {:.0f}s  '
              'β={:.5f} ω={:.5f} η0={:.5f} X0={:.5f}{}'.format(
                  key, rec['preferences'], rec['residual'], rec.get('verifyResidual', np.nan),
                  rec['nfev'], rec['time'], rec['β'], rec['ω'], rec['η0'], rec['X0'],
                  '  X={:.4f} hbar={:.4f}'.format(rec['X'], rec['hbar']) if 'X' in rec else ''))
        return rec

    def onPoint(r):
        """ Persist after every successful point. Failures are only *reported* here -- a failed attempt
        does not mean a failed point, since the march still has the carried start and step-halving to try,
        so which values genuinely failed is not known until the march returns. """
        key = round(float(r['value']), 6)
        if not r['ok']:
            print('  rho={:<6} attempt failed ({}): {}'.format(key, r['x0Source'], r['error'][:100]))
            return
        if r['result'].get('cached'):
            return
        rows[key] = toRow(r['result'], r['requested'], commit)
        write()
        with open(os.path.join(pkldir, 'rho_{:.4f}.pkl'.format(r['value'])), 'wb') as f:
            pickle.dump(m, f)

    x0 = None
    if args.x0 is not None:
        x0Pars = dict(zip(('β', 'ω', 'η0', 'X0'), args.x0))
        m._calSetPars(x0Pars)
        x0 = m._calToX(x0Pars)
        print('anchor starting point: ' + ', '.join('{}={:.4f}'.format(k, v) for k, v in x0Pars.items()))
    out = continuation.marchGrid(grid, solve, x0 = x0, anchor = args.anchor, degree = args.degree,
                                 maxHalvings = args.maxHalvings, onPoint = onPoint)

    # A failed *attempt* is not a failed point -- the retry ladder and step-halving may have rescued it.
    # Only values with no successful record anywhere are written as failures.
    solvedValues = {round(float(r['value']), 6) for r in out['records'] if r['ok']}
    hardFails = [r for r in out['failures'] if round(float(r['value']), 6) not in solvedValues]
    for r in hardFails:
        rows[round(float(r['value']), 6)] = {
            'ρ': round(float(r['value']), 6), 'requested': r['requested'], 'commit': commit,
            'residual': np.nan, 'preferences': 'FAILED: ' + r['error'][:120],
            'nEqMax': -1, 'nCandMax': -1, 'nFallback': -1,
            'timestamp': datetime.datetime.now().isoformat(timespec = 'seconds')}
    write()

    solved = [r for r in out['records'] if r['ok']]
    print('\n{} solved ({} requested, {} inserted by step-halving), {} failed'.format(
        len(solved), sum(r['requested'] for r in solved), sum(not r['requested'] for r in solved),
        len(hardFails)))
    print('csv: ' + os.path.relpath(args.out, REPO))
    if hardFails:
        print('failed at rho: ' + ', '.join('{:.4f}'.format(r['value']) for r in hardFails))
    final = pd.DataFrame(list(rows.values())).reindex(columns = COLUMNS).sort_values('ρ')
    bad = final[final['verifyResidual'] > 1e-3]
    if len(bad):
        print('\nWARNING: {} point(s) converged but did not survive grid refinement -- do not read these '
              'as results, re-run them with a larger --nι/--ns:'.format(len(bad)))
        print(bad[['ρ', 'residual', 'verifyResidual']].to_string(index = False))
    return 0 if not hardFails else 1


if __name__ == '__main__':
    sys.exit(main())

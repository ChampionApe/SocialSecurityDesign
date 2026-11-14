r""" Comparative statics of the politico-economic equilibrium on a full CARTESIAN GRID of (eps, theta) at
the calibration year t0. Not a test -- an experiment script. Produces the data behind the paper figure
that plots each variable against eps with one line per theta (shaded between adjacent theta lines, theta
on a colorbar), which needs a genuine product grid rather than a cross through the calibrated point.

Run:  .venv\Scripts\python.exe python\InformalSavings\sweepEpsThetaGrid.py
      ... --nEps 41 --nTheta 21             finer grids
      ... --eps 0.1 0.3 0.5 --theta 0.7 0.9 explicit grids
      ... --force                           re-solve rows already in the csv

NOT A RECALIBRATION. The calibrated rho=1 instance is loaded and (beta, omega, eta0, X0) stay pinned at
their calibrated values at every grid point; eps and theta move and the full politico-economic equilibrium
is re-solved over the whole horizon at each. Both calibrated values are forced into their grids, so the
grid carries exactly one row with BOTH at their calibrated values (statusQuo) which must reproduce the
baseline solve -- that row is the check that the installation routines left nothing stale.

EPS AND THETA ARE TWO INDEPENDENT AXES HERE, BY CONSTRUCTION. model.getEps makes eps a decreasing function
of theta, so a recalibration would move eps whenever theta moves; on a product grid that would collapse
the grid onto a curve and there would be no theta-family of eps-curves to plot. theta is therefore
installed on its own and eps is installed independently afterwards. The other reading -- re-applying the
calibration formula for eps at each new theta -- is what a recalibration would do, and has no meaning on a
product grid; it is not offered.

WHAT INSTALLING EACH PARAMETER REQUIRES

  eps: shockUniversal.installEps, never a bare `eps=` argument to solvePEE_*. kappa_t(eps_{t+1}) is
    consumed everywhere through a CACHED db['kappa'], so passing a new eps while leaving db alone solves
    the households' problem against one government budget and the government's against another. Nothing
    raises -- a mutually inconsistent (eps, kappa) violates no equilibrium condition, since kappa enters
    everywhere as a given (README, "Open items").

  theta: installTheta, which writes db['theta'] plus its lead/lag through adjPar. theta has no cached
    descendant EXCEPT through eps, and that chain is deliberately cut here.

ORDER: theta first, then eps. installTheta touches only db['theta']; installEps rewrites db['eps'] and
db['kappa'] (aux_kappa reads db['eps[t+1]'] and nothing else), so neither can clobber the other. solvePoint asserts the (eps, theta) actually sitting in db are the intended pair, and both are
recorded on every row.

REPORTED AT ONE YEAR. Every column is read at t0 = db['t'][db['t0']] (model index 3), the calibration
year, the figure's "year 2010". The pickled instances carry neither db['dates'] nor db['workweek'], so h
is the model's own aggregate hours, not a workweek in hours.

RESUMABLE, AND THE ROW DOES NOT RECORD WHAT PRODUCED IT. The csv is rewritten after every point and a
re-run returns any (eps, theta) already present in it without re-solving. The key is the parameter point
alone, so anything else that moved is invisible to the resume and needs a new --out or --force:

  * a RECALIBRATION -- the loaded instance carries (beta, omega, eta0, X0) and no column names them. This
    is the one that has actually bitten: the 2026-08-24 K/Y retarget moved the calibrated eps, so a
    resumed run ADDED the new calibrated column and kept 378 rows solved at the old beta, leaving a
    two-calibration csv that is a complete rectangle and passes every shape check but one.
  * changed settings -- a different --interpKind/--smoothKnots/--nι/--ns.

paper/datasets.epsThetaGrid now rejects a csv with two statusQuo rows, or one whose statusQuo does not
match the current calibration, which is what makes the first case loud rather than a kink in a figure.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd
from gridsearch.testing import utf8Stdout

utf8Stdout()                                    # progress lines and --help carry Greek
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
os.chdir(HERE)                                  # relative paths resolve here; test.py finds data/ itself

from shockUniversal import loadCalibrated, installEps, solvePEE, universalEps, PKLDIR
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'sweeps')
# Ordered so the csv reads as: the (eps, theta) the row is keyed on, the system it defines and which two
# distinguished points it is, then the equilibrium it produced, then what it cost. nEqMax/nCandMax/
# nFallback are the tax candidates' counts of the row's solve (num_robustroot.tex): non-key (#13), and
# NaN on a row resumed from a csv written before they existed.
COLUMNS = ['eps', 'theta', 'ρ', 'statusQuo', 'universalEps',
           'τ', 'sr', 'srGDP', 'h', 'ι', 's', 's_', 'Y', 'Y0', 'c10', 'c20', 'nRoots', 'nEqMax', 'nCandMax',
           'nFallback', 'time']
# sr is savings over formal output (Base.savingsRate); srGDP is savings over GDP, s/(Y+Y0), with Y and Y0
# formal and informal output (EE_report), the concept of the tax target (model.taxTarget).


def installTheta(m, θ):
    """ Write a theta path into db, and nothing behind it. theta is a 1D parameter, so adjPar gives
    db['theta'] its lead and lag; nothing else caches theta directly.

    The one derived chain, theta -> eps (model.getEps) -> db['kappa'], is deliberately NOT followed --
    see the module docstring. eps is installed separately by the caller, and installEps refreshes
    db['kappa'] itself, so the cached kappa is consistent with the eps that ends up in db. """
    m.db.update(m.adjPar('θ', np.full(m.T, float(θ))))


def atT0(m, sol, t0):
    """ The reported equilibrium at the single year t0, as the csv's columns. s_/h/c10/c20 report on
    db['t'] and iota on db['txE'] (README, "Conventions and traps"), so both are looked up by label at t0
    rather than positionally. The tax candidates' counts cover the whole recursion, not t0 alone. """
    r = sol['report']
    s, s_, h, Y, Y0 = (float(r[k].xs(t0)) for k in ('s', 's_', 'h', 'Y', 'Y0'))
    return {'τ': float(sol['τ'].xs(t0)), 'sr': float(m.B.savingsRate(s, s_, h, t0)), 'srGDP': s/(Y+Y0),
            'h': h, 'ι': float(r['ι'].xs(t0)), 's': s, 's_': s_, 'Y': Y, 'Y0': Y0,
            'c10': float(r['c10'].xs(t0)), 'c20': float(r['c20'].xs(t0)),
            'nRoots': (sol['init'] or {}).get('nRoots')} | multiplicityColumns(sol.get('multiplicity'))


def buildGrid(explicit, lo, hi, n, pinned):
    """ A grid that contains `pinned` exactly. The pinned values (the calibrated parameter, and for eps
    the universal target) are inserted rather than approximated: statusQuo must land on a row, and the
    eps=eps^U row is what the universalisation shock is compared against. """
    g = np.asarray(explicit, dtype = float) if explicit else np.linspace(lo, hi, n)
    return np.unique(np.concatenate([g, np.asarray(pinned, dtype = float)]))


def solvePoint(ε, θ, ρ, settings, pkldir, cal, εU):
    """ One grid point: a freshly loaded calibrated instance, both characteristics installed, the whole
    PEE re-solved. Reloading per point rather than mutating one instance keeps every point independent of
    the order they were solved in -- solvePEE_* caches warm starts on the model (m.x0), so a shared
    instance would make a row depend on its predecessor. """
    t = time.time()
    m, preferences = loadCalibrated(ρ, settings, pkldir = pkldir)
    installTheta(m, θ)                          # must not track eps: the two axes are independent here
    installEps(m, np.full(m.T, float(ε)))
    εIn, θIn = float(m.db['eps'].iloc[0]), float(m.db['θ'].iloc[0])
    assert np.isclose(εIn, ε, rtol = 0, atol = 1e-12) and np.isclose(θIn, θ, rtol = 0, atol = 1e-12), \
        'installed (ε,θ)=({}, {}) but asked for ({}, {})'.format(εIn, θIn, ε, θ)
    t0 = m.db['t'][m.db['t0']]
    row = {'eps': εIn, 'theta': θIn, 'ρ': float(ρ),
           'statusQuo': bool(np.isclose(εIn, cal['eps']) and np.isclose(θIn, cal['theta'])),
           'universalEps': bool(np.isclose(εIn, εU))}
    row |= atT0(m, solvePEE(m, preferences), t0)
    row['time'] = time.time() - t
    return row


def readDone(path, force):
    """ {(rounded eps, rounded theta): row} for the points already in the csv. Rounded so a grid rebuilt
    from the same --lo/--hi/--n matches it exactly rather than by floating-point luck. """
    if force or not os.path.exists(path):
        return {}
    df = pd.read_csv(path)
    return {(round(float(r['eps']), 10), round(float(r['theta']), 10)): r
            for r in df.to_dict('records')}


def main():
    p = argparse.ArgumentParser(description = 'Politico-economic equilibrium on a cartesian (eps, theta) '
                                              'grid at the calibrated rho=1 instance.')
    p.add_argument('--rho', type = float, default = 1.0)
    p.add_argument('--eps', type = float, nargs = '+', default = None,
                   help = 'explicit eps grid; overrides --epsLo/--epsHi/--nEps')
    p.add_argument('--epsLo', type = float, default = 0.02)
    p.add_argument('--epsHi', type = float, default = 0.65)
    p.add_argument('--nEps', type = int, default = 25,
                   help = 'before the calibrated eps and the universal eps^U are inserted')
    p.add_argument('--theta', type = float, nargs = '+', default = None,
                   help = 'explicit theta grid; overrides --thetaLo/--thetaHi/--nTheta')
    p.add_argument('--thetaLo', type = float, default = 0.5)
    p.add_argument('--thetaHi', type = float, default = 1.0)
    p.add_argument('--nTheta', type = int, default = 13,
                   help = 'before the calibrated theta is inserted')
    p.add_argument('--nι', '--niota', dest = 'nι', type = int, default = 45)
    p.add_argument('--ns', type = int, default = 45)
    p.add_argument('--interpKind', default = 'cubic', choices = ('linear', 'cubic', 'pchip'))
    p.add_argument('--smoothKnots', type = int, default = 4)
    p.add_argument('--pkldir', default = PKLDIR, help = 'directory of pickled calibrated instances')
    p.add_argument('--force', action = 'store_true', help = 're-solve rows already present in the csv')
    p.add_argument('--out', default = 'epsThetaGrid_rho{ρ:.4f}.csv',
                   help = 'relative to results/sweeps/, or an absolute path')
    a = p.parse_args()
    settings = {'nι': a.nι, 'ns': a.ns, 'interpKind': a.interpKind, 'smoothKnots': a.smoothKnots}

    os.makedirs(OUTDIR, exist_ok = True)
    out = a.out.format(ρ = a.rho)
    path = out if os.path.isabs(out) else os.path.join(OUTDIR, out)

    # The calibrated pair and eps^U are read off the instance, not hard-coded: they are what the grids are
    # pinned to and what statusQuo/universalEps are decided against, and all three move whenever the
    # calibration is re-run.
    m0, preferences = loadCalibrated(a.rho, settings, pkldir = a.pkldir)
    cal = {'eps': float(m0.db['eps'].iloc[0]), 'theta': float(m0.db['θ'].iloc[0])}
    εU = float(universalEps(m0)[0])
    εGrid = buildGrid(a.eps, a.epsLo, a.epsHi, a.nEps, [cal['eps'], εU])
    θGrid = buildGrid(a.theta, a.thetaLo, a.thetaHi, a.nTheta, [cal['theta']])

    print('ρ={}  {} solver, T={}, t0 = db index {}   grid {}'.format(
        a.rho, preferences, m0.T, m0.db['t'][m0.db['t0']], settings))
    print('calibrated: ε={:.6f}  θ={:.6f}   universal ε^U={:.6f}'.format(cal['eps'], cal['theta'], εU))
    print('cartesian grid: {} ε in [{:.4f}, {:.4f}] x {} θ in [{:.4f}, {:.4f}] = {} points'.format(
        len(εGrid), εGrid[0], εGrid[-1], len(θGrid), θGrid[0], θGrid[-1], len(εGrid)*len(θGrid)))

    rows = readDone(path, a.force)
    todo = [(ε, θ) for θ in θGrid for ε in εGrid]
    if rows:
        print('resuming: {} of {} points already in {}'.format(
            sum((round(float(ε), 10), round(float(θ), 10)) in rows for ε, θ in todo), len(todo),
            os.path.relpath(path, REPO)))

    def write():
        (pd.DataFrame(list(rows.values())).reindex(columns = COLUMNS)
           .sort_values(['theta', 'eps']).to_csv(path, index = False))

    for n, (ε, θ) in enumerate(todo, start = 1):
        key = (round(float(ε), 10), round(float(θ), 10))
        if key in rows:
            continue
        row = solvePoint(ε, θ, a.rho, settings, a.pkldir, cal, εU)
        rows[key] = row
        write()
        print('  [{:>4}/{}] ε={:.6f} θ={:.6f}{}  τ={:.6f}  sr={:.6f}  h={:.6f}  ι={:.6f}  '
              'nRoots={}  {:.1f}s'.format(n, len(todo), row['eps'], row['theta'],
                                          ' *SQ*' if row['statusQuo'] else
                                          (' *εU*' if row['universalEps'] else '     '),
                                          row['τ'], row['sr'], row['h'], row['ι'],
                                          row['nRoots'], row['time']))
    write()

    df = pd.DataFrame(list(rows.values())).reindex(columns = COLUMNS).sort_values(['theta', 'eps'])
    print('\n{} rows written: {}'.format(len(df), os.path.relpath(path, REPO)))
    bad = df[df['nRoots'] != 1]
    if len(bad):
        print('WARNING: {} row(s) with nRoots != 1:'.format(len(bad)))
        print(bad[['eps', 'theta', 'τ', 'nRoots']].to_string(index = False))
    bound = df[(df['nEqMax'] > 1) | (df['nFallback'] > 0)]
    if len(bound):
        print('WARNING: {} row(s) where the equilibrium selection bound (nEqMax > 1) or fell back:'.format(len(bound)))
        print(bound[['eps', 'theta', 'τ', 'nEqMax', 'nCandMax', 'nFallback']].to_string(index = False))
    sq = df[df['statusQuo'] == True]
    print('status quo row (must reproduce the baseline solve):')
    print(sq[['eps', 'theta', 'τ', 'sr', 'h', 'ι']].to_string(index = False))
    return 0


if __name__ == '__main__':
    sys.exit(main())

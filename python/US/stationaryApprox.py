r""" How far is a stationary policy function from the exact date-specific one along the US demographic path?

Run:  .venv\Scripts\python.exe python\US\stationaryApprox.py --theta C --main   # rho = 0.5, 2 (config.US['ρStationary'])
      ... --rho 2 --Tstat 15 --ns 150 --n 101 [--vectorX]

The question (sec:numerical's literature paragraph). The numerical politico-economic literature computes a
TIME-INVARIANT policy function tau(s) as a fixed point, the right object in a stationary environment. Along a
demographic transition the Markov policy function is date-specific: the continuation at t is the policy of
the economy at t+1, with its own nu_{t+1}. Evaluating a stationary function along the projected nu_t path
holds the continuation at its stationary form, i.e. a steady-state approximation. This script measures what
that approximation would cost in OUR model, under CRRA (under LOG the periods decouple and the two coincide).

Two objects per rho, at the calibrated (beta, omega) of the headline sweep:
  exact   the paper's solution: CRRA.solveBackward over the 1960-2200 horizon, walked forward from the
          steady-state initial state (solvePEE_CRRA).
  stat    for each distinct nu_t on the path, the stationary policy function: the backward recursion over a
          long horizon (Tstat periods) with nu frozen at nu_t, read at its FIRST period (its convergence is
          reported as the max gap between the first two periods' policies on the state grid). The
          approximate path walks tau_t = tau^stat_{nu_t}(s_{t-1}), s_t = s^stat_{nu_t}(s_{t-1}) from the
          SAME initial state as the exact path.
Also reported: the policy-function gap at the exact path's own states, tau^stat_{nu_t}(s^exact_{t-1}) -
tau^exact_t(s^exact_{t-1}), which isolates the error in the function from its propagation through the state.

The exact solution's last periods (constant nu from 2080 on, three appended steady-state periods) carry the
terminal condition rather than demographic change: read the 1960-2080 rows for the non-stationarity cost and
the tail for the finite-horizon effect. Both solves share one state grid (the baseline's default sGrid) and
grid settings, so the gaps are not interpolation artefacts of different grids.

Writes results/numerical/US_stationaryApprox.csv, one row per (rho, t), each with the exact solve's tax counts
(nEqMax, nCandMax, nFallback; policy.multiplicityColumns), merged into the csv on (rho, t) and written atomically
(writeRows): a run at some rho keeps the other rho's rows at its θ and drops the rows at another θ (one θ per csv).
A relative --out is the caller's working directory's. --theta C --main is the pipeline's; --theta index (or an --out
outside results/) the workbook's θ, which without --theta is refused into results/ (oecd26.indexTheta).
"""
import os, sys, argparse, time
import numpy as np, pandas as pd
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                 # where a relative --out points: the chdir below moves to HERE
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import test as testmod
from model import ModelUS
from policy import multiplicityColumns

OUT = os.path.join(REPO, 'results', 'numerical', 'US_stationaryApprox.csv')
# --theta VARIANT --main (the pipeline's main θ variant, oecd26.thetaRoot): θ0 = the variant's U.S. row of
# data/oecd26_theta.csv imposed on every model here (db['θ0'], ModelUS.getθ), the sweep row read at that θ0
# (oecd26.usBeta), and thetaSource, θ0 on every output row. Set by main.
THETA = {'variant': None, 'θ0': None}


def usRow(ρ, commonX):
    csv = os.path.join(REPO, 'results', 'calibration', 'US_rhoGridCommonX.csv' if commonX else 'US_rhoGrid.csv')
    df = pd.read_csv(csv)
    hit = df.loc[(df['ρ'] - ρ).abs() < 1e-9]
    if hit.empty:
        raise KeyError(f'No US calibration at ρ={ρ} in {csv}.')
    if THETA['variant'] is not None:
        import oecd26
        oecd26.usBeta(csv, ρ, thetaSource = THETA['variant'])          # ValueError at another θ0
    return hit.iloc[-1]


def build(row, ρ, commonX, gs, ν = None, T = None):
    """ The calibrated US model at ρ; with ν a scalar and T an int, the same economy with ν frozen at that
    value over a T-period horizon (the stationary problem). """
    pars = (dict(testmod.pars) | {'ρ': float(ρ), 'β': float(row['β']), 'ω': float(row['ω'])}
            | ({} if THETA['θ0'] is None else {'θ0': THETA['θ0']}))
    kwargs = dict(testmod.kwargs)
    if ν is not None:
        pars['ν'] = np.full(T, float(ν))
        kwargs['T'] = T
    m = ModelUS(pars = pars, commonX = commonX, **kwargs)
    if ν is None:
        m.db['dates'] = testmod.dates
    m.CRRA.initGS(gs)
    if commonX:
        m.initProductivity_commonX(X = float(row['X']))
        m.updateAuxPars()
    return m


def stationaryPolicy(row, ρ, commonX, gs, ν, T, sGrid, sRange):
    """ The stationary policy function at ν: first period of a T-period backward recursion with ν frozen.
    Returns (report dict of that period, convergence gap between periods 1 and 2 over the nodes of sGrid
    inside sRange, τ*). The gap is confined to sRange (the states the exact path visits, padded) because
    the recursion does not settle at the very bottom of the grid (s ~ 1e-4, an order of magnitude below any
    visited state), where the gap stays ~1e-4 however long the horizon; in the visited range it is ~1e-6. τ* is
    the steady-state politico-economic tax at ν: the fixed point at which the stationary policy function,
    evaluated at the steady-state savings under a constant τ*, returns τ* (model.steadyStatePEE_CRRA's
    construction). τ* is what a steady-state comparison across demographic scenarios reports. """
    ms = build(row, ρ, commonX, gs, ν = ν, T = T)
    θ, ε = ms.db['θ'].values, ms.db['eps'].values
    sols = ms.CRRA.solveBackward(θ, ε, sGrid = sGrid)
    t0, t1 = ms.db['t'][0], ms.db['t'][1]
    vis = (sGrid >= sRange[0]) & (sGrid <= sRange[1])
    conv = float(np.nanmax(np.abs(sols[t0]['τ'].values - sols[t1]['τ'].values)[vis]))
    l, u = ms.CRRA.GS['PEE']['gridSettings']['l'], ms.CRRA.GS['PEE']['gridSettings']['u']
    residual = lambda τ: τ - np.clip(sols[t0]['τPolicy'](ms.steadyState_CRRA_solve(τ, θ[0], t = t0)['s']), l, u)
    τStar = float(optimize.brentq(residual, l, u))
    return sols[t0], conv, τStar


def runRho(ρ, commonX, gs, Tstat):
    row = usRow(ρ, commonX)
    m = build(row, ρ, commonX, gs)
    tIdx = m.db['t']
    θ, ε = m.db['θ'].values, m.db['eps'].values
    ν = m.db['ν'].values.astype(float)
    l, u = m.CRRA.GS['PEE']['gridSettings']['l'], m.CRRA.GS['PEE']['gridSettings']['u']
    sGrid = m.CRRA.defaultSGrid(θ[-1], tIdx[-1], n = gs['ns'])

    tic = time.time()
    exact = m.solvePEE_CRRA(backwardKwargs = {'sGrid': sGrid})
    τE, rep = exact['τ'].values, exact['report']
    sE_ = rep['s_'].values                                  # s_{t-1} entering each period, length T
    s0 = float(sE_[0])
    print(f'rho={ρ}: exact solve {time.time()-tic:.0f}s, tau_2020={τE[m.db["t0"]]:.5f}')

    stat, conv, τStar = {}, {}, {}
    for νt in sorted(set(np.round(ν, 12))):
        tic = time.time()
        stat[νt], conv[νt], τStar[νt] = stationaryPolicy(row, ρ, commonX, gs, νt, Tstat, sGrid,
                                                         (0.8 * sE_.min(), 1.2 * sE_.max()))
        print(f'   stationary at nu={νt:.4f}: {time.time()-tic:.0f}s, first-vs-second-period gap {conv[νt]:.2e}, '
              f'steady-state PEE tau*={τStar[νt]:.5f}')
    νLR = np.round(ν[-1], 12)                               # the long-run (terminal) demographics

    # walk the stationary functions forward from the exact initial state
    τA, sA_ = np.empty(m.T), np.empty(m.T)
    s_ = s0
    for pos in range(m.T):
        νt = np.round(ν[pos], 12)
        sA_[pos] = s_
        τA[pos] = float(np.clip(stat[νt]['τPolicy'](s_), l, u))
        if pos < m.T - 1:
            s_ = float(stat[νt]['sPolicy'](s_))
    # the stationary function at the exact path's own states
    τAatE = np.array([float(np.clip(stat[np.round(ν[pos], 12)]['τPolicy'](sE_[pos]), l, u)) for pos in range(m.T)])
    # one stationary function, computed at the long-run demographics, applied at every date (at the exact states)
    τLR = np.array([float(np.clip(stat[νLR]['τPolicy'](sE_[pos]), l, u)) for pos in range(m.T)])
    # the steady-state comparison: the stationary equilibrium tax at each date's own nu
    τSS = np.array([τStar[np.round(ν[pos], 12)] for pos in range(m.T)])

    # exact economic equilibrium at the approximate tax path, for allocations comparable to the paper's
    solA = m.EE_CRRA_solve(τA, θ, ε, s0)
    repA = m.EE_report(solA, τA, θ, ε, s0)
    def srOverY(r):
        """ s/Y at every t, the paper's savings-rate unit (shocks.readout's 'srOverY'). """
        return np.array([float(m.B.savingsRate(r['s'].xs(t), r['s_'].xs(t), r['h'].xs(t), t)) for t in tIdx])

    rows = []
    for pos, t in enumerate(tIdx):
        rows.append({'ρ': ρ, 't': int(t), 'date': int(testmod.dates[pos]), 'ν': ν[pos],
                     'τ_exact': τE[pos], 'τ_stat': τA[pos], 'τ_statAtExactState': τAatE[pos],
                     'τ_ssPEE': τSS[pos], 'τ_statLongRun': τLR[pos],
                     'gap_pp': 100 * (τA[pos] - τE[pos]), 'gapAtExactState_pp': 100 * (τAatE[pos] - τE[pos]),
                     'gapSS_pp': 100 * (τSS[pos] - τE[pos]), 'gapLongRun_pp': 100 * (τLR[pos] - τE[pos]),
                     's_exact': sE_[pos], 's_stat': sA_[pos],
                     'srY_exact': srOverY(rep)[pos], 'srY_stat': srOverY(repA)[pos],
                     'R_exact': rep['R'].values[pos], 'R_stat': repA['R'].values[pos],
                     'statConv': conv[np.round(ν[pos], 12)],
                     'inGrid_stat': bool(sGrid[0] <= sA_[pos] <= sGrid[-1])}
                    | multiplicityColumns(exact.get('multiplicity')))
    return pd.DataFrame(rows)


KEY = ['ρ', 't']                  # the csv's row key (python/paper/datasets.stationaryApprox)


def writeRows(df, path):
    """ Merge this run's rows `df` into csv `path` and write it through runESC.writeCsv (atomic): rows on disk at a
    key (ρ, t) of `df` are replaced, the others kept if at the run's θ (oecd26.rowsAtTheta at THETA) and dropped
    otherwise. Carried rows are read round trip, so their digits survive. Returns the written frame, sorted on KEY. """
    import oecd26, runESC
    os.makedirs(os.path.dirname(path), exist_ok = True)
    out = df
    if os.path.exists(path):
        old = pd.read_csv(path, float_precision = 'round_trip')
        keep = oecd26.rowsAtTheta(old, THETA['variant'], THETA['θ0'])
        if len(keep) < len(old):
            print(f'  {os.path.relpath(path, REPO)}: {len(old) - len(keep)} rows at another θ dropped')
        k = lambda d: [(round(float(ρ), 9), int(t)) for ρ, t in zip(d['ρ'], d['t'])]
        mine = set(k(df))
        out = pd.concat([keep[[x not in mine for x in k(keep)]], df], ignore_index = True)
    out = out.sort_values(KEY, kind = 'stable').reset_index(drop = True)
    runESC.writeCsv(out, path)
    return out


def main():
    p = argparse.ArgumentParser(description = 'Stationary-policy approximation vs exact date-specific policies (US, CRRA).')
    p.add_argument('--rho', type = float, nargs = '*', default = [0.5, 2.0])      # config.US['ρStationary']
    p.add_argument('--Tstat', type = int, default = 40, help = 'horizon of the frozen-nu recursion')
    p.add_argument('--n', type = int, default = 101)
    p.add_argument('--ns', type = int, default = 150)
    p.add_argument('--smoothKnots', type = int, default = 4)
    p.add_argument('--vectorX', action = 'store_true', help = 'the vector-X calibration variant (headline is common X)')
    p.add_argument('--out', default = OUT)
    p.add_argument('--theta', default = None,
                   help = "the pipeline's main θ variant (with --main); 'index' the workbook's θ")
    p.add_argument('--main', action = 'store_true', help = "with --theta config.US['theta']")
    a = p.parse_args()
    a.out = os.path.abspath(os.path.join(CWD, a.out))
    import oecd26
    if a.theta in (None, oecd26.INDEX) and not a.main:
        a.theta = oecd26.indexTheta(a.theta, [a.out], 'stationaryApprox.py')
    if a.theta is not None or a.main:
        if not (a.theta is not None and a.main):
            raise SystemExit('--theta needs --main here (and --main --theta)')
        try:
            oecd26.thetaRoot(a.theta, main = True)
            THETA.update({'variant': a.theta, 'θ0': oecd26.thetaFromFile('USA', a.theta)})
        except (ValueError, FileNotFoundError) as e:
            raise SystemExit('--theta {}: {}'.format(a.theta, e))
    gs = {'n': a.n, 'ns': a.ns, 'smoothKnots': a.smoothKnots or None, 'interpKind': 'linear'}
    dfs = [runRho(ρ, not a.vectorX, gs, a.Tstat) for ρ in a.rho]
    df = pd.concat(dfs, ignore_index = True)
    if THETA['variant'] is not None:
        df = df.assign(thetaSource = THETA['variant'], θ0 = THETA['θ0'])
    writeRows(df, a.out)
    pd.set_option('display.width', 200)
    print('\n' + df[['ρ', 'date', 'ν', 'τ_exact', 'τ_stat', 'τ_ssPEE', 'τ_statLongRun', 'gap_pp', 'gapAtExactState_pp',
                     'gapSS_pp', 'gapLongRun_pp', 'srY_exact', 'srY_stat', 'statConv', 'inGrid_stat']].round(5).to_string(index = False))
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()

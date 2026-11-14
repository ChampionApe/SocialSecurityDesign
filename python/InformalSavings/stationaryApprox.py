r""" How far is a stationary policy function from the exact date-specific one along Argentina's demographic path?

Run:  .venv\Scripts\python.exe python\InformalSavings\stationaryApprox.py                # rho = 0.5, 1, 2
      ... --rho 1 --pkldir results\calibration\instancesCommonX --out ARG_stationaryApprox_commonX.csv
      grid flags as shockUniversal.py (--nι --ns --interpKind --smoothKnots); the instance is re-solved on the
      settings it was calibrated at (loadCalibrated's split between the two solvers).

The Argentina counterpart of python/US/stationaryApprox.py (see its docstring for the question). Here even
LOG carries a state, the informal savings ratio iota_{t-1}, so the stationary approximation is not exact at
rho = 1, and the demographic change is large: nu falls from 1.65 to 0.96 across the horizon.

Objects per rho, on the pickled calibrated instance:
  exact  the paper's solution: solvePEE_LOG/CRRA (backward recursion over the 1920-2190 horizon, initial
         state from the steady-state fixed point, exact economic equilibrium at the walked tax path).
  stat   for each distinct nu_t, the stationary policy function: the same recursion on a deep copy with nu
         frozen at nu_t over the whole horizon, read at its FIRST period, on the BASELINE's state grids
         (so the two solutions are interpolated on the same nodes). Its convergence is reported as the gap
         between the first two periods' policies over the states the exact path visits. The approximate
         path walks tau_t = tau^stat_{nu_t}(state_{t-1}) from the exact initial state, with the state
         transition re-solved on the frozen model (as approximatePEE(exact=True) does).
  ssPEE  the steady-state politico-economic tax at each nu_t (initialStatePEE on the frozen model): what a
         comparison of steady states across demographic scenarios reports.
  statLongRun  one stationary function, at the terminal nu, evaluated at the exact path's own states.

The exact solution's last periods carry the terminal condition rather than demographic change (nu is
constant from 2070 on and the last three periods are appended steady-state periods); read the dated rows
for the non-stationarity cost and the tail for the finite-horizon effect.

Writes results/numerical/ARG_stationaryApprox{,_commonX}.csv, one row per (rho, t), with the tax
candidates' counts nEqMax/nCandMax/nFallback over every recursion behind a rho (num_robustroot.tex).
"""
import os, sys, argparse, time
from copy import deepcopy
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import test as testmod
from shockUniversal import loadCalibrated, solvePEE, PKLDIR
from policy import multiplicitySummary, multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'numerical')


def frozen(m, ν):
    """ A deep copy of m with nu frozen at `nu` over the whole horizon and the warm-start caches cleared. """
    ms = deepcopy(m)
    ms.x0 = {}
    ms.db.update(ms.adjPar('ν', np.full(ms.T, float(ν))))
    ms.updateAuxPars()
    return ms


def gridsOf(sol, prefs):
    """ The state and candidate grids a solved period used, as solveBackward keyword arguments. """
    if prefs == 'LOG':
        return {'ιGrid': sol['ι_'], 'ιCandGrid': sol['ιCand']}
    return {'sGrid': sol['s_'], 'ιGrid': sol['ι_'], 'sCandGrid': sol['sCand'], 'ιCandGrid': sol['ιCand']}


def convergence(sol0, sol1, prefs, sRange, ιRange):
    """ max |tau_0 - tau_1| over the grid nodes inside the visited ranges. """
    if prefs == 'LOG':
        ι = sol0['τ'].index.values
        vis = (ι >= ιRange[0]) & (ι <= ιRange[1])
        return float(np.nanmax(np.abs(sol0['τ'].values - sol1['τ'].values)[vis]))
    s, ι = sol0['τ'].index.values, sol0['τ'].columns.values
    vis = ((s >= sRange[0]) & (s <= sRange[1]))[:, None] & ((ι >= ιRange[0]) & (ι <= ιRange[1]))[None, :]
    return float(np.nanmax(np.abs(sol0['τ'].values - sol1['τ'].values)[vis]))


def stationary(m, ν, prefs, grids, sRange, ιRange):
    """ The frozen-nu model, its policy functions read at the first period, the convergence gap, the
    steady-state politico-economic tax tau* at nu, and the tax candidates' counts of its recursion. """
    ms = frozen(m, ν)
    θ, ε = ms.db['θ'].values, ms.db['eps'].values
    sols = getattr(ms, prefs).solveBackward(θ, ε, **grids)
    t0, t1 = ms.db['t'][0], ms.db['t'][1]
    conv = convergence(sols[t0], sols[t1], prefs, sRange, ιRange)
    try:
        τStar = float(ms.initialStatePEE(sols, θ, ε, prefs)['τ'])
    except RuntimeError as e:
        print('   steady-state PEE at nu={:.4f} failed: {}'.format(ν, str(e).split('\n')[0][:120]))
        τStar = np.nan
    return {'m': ms, 'sols': sols, 'conv': conv, 'τStar': τStar, 'mult': multiplicitySummary(sols)}


def runRho(ρ, settings, pkldir):
    m, prefs = loadCalibrated(ρ, settings, pkldir = pkldir)
    tIdx, T = m.db['t'], m.T
    ν = m.db['ν'].values.astype(float)
    pol = getattr(m, prefs)
    l, u = pol.GS['PEE']['gridSettings']['l'], pol.GS['PEE']['gridSettings']['u']
    isLOG = prefs == 'LOG'

    tic = time.time()
    exact = solvePEE(m, prefs)
    τE, rep, init = exact['τ'].values, exact['report'], exact['init']
    sE_ = rep['s_'].values                                            # s_{t-1} entering each period
    ιE_ = np.concatenate([[init['ι']], rep['ι'].values])   # iota entering each period: init, then ι_t (txE, T-1)
    print(f'rho={ρ} ({prefs}): exact solve {time.time()-tic:.0f}s, tau_t0={τE[m.db["t0"]]:.5f}')
    grids = gridsOf(exact['sols'][tIdx[0]], prefs)
    sRange, ιRange = (0.8*sE_.min(), 1.2*sE_.max()), (0.8*ιE_.min(), 1.2*ιE_.max())

    stat = {}
    for νt in sorted(set(np.round(ν, 12))):
        tic = time.time()
        stat[νt] = stationary(m, νt, prefs, grids, sRange, ιRange)
        print(f'   stationary at nu={νt:.4f}: {time.time()-tic:.0f}s, first-vs-second-period gap '
              f'{stat[νt]["conv"]:.2e}, steady-state PEE tau*={stat[νt]["τStar"]:.5f}')
    νLR = np.round(ν[-1], 12)

    def τAt(st, s_, ι_):
        pol0 = st['sols'][st['m'].db['t'][0]]['τPolicy']
        return float(np.clip(pol0(ι_) if isLOG else pol0(s_, ι_), l, u))

    # walk the stationary functions forward from the exact initial state, re-solving the transition on
    # the frozen model of that date
    τA, sA_, ιA_ = np.empty(T), np.empty(T), np.empty(T)
    s_, ι_ = float(sE_[0]), float(ιE_[0])
    for pos in range(T):
        st = stat[np.round(ν[pos], 12)]
        ms, sols = st['m'], st['sols']
        t0, t1 = ms.db['t'][0], ms.db['t'][1]
        θs, εs = ms.db['θ'].values, ms.db['eps'].values
        sA_[pos], ιA_[pos] = s_, ι_
        τA[pos] = τAt(st, s_, ι_)
        if pos < T - 1:
            if isLOG:
                r, _ = ms.LOG.solveStateApprox_t(np.array([τA[pos]]), sols[t0]['ιCand'], t0, θs[1], εs[1], sols[t1])
                ι_ = float(r[0])
            else:
                sSol, ιSol, _, _ = ms.CRRA.solveStateApprox_t(np.array([τA[pos]]), np.array([s_]), sols[t0]['sCand'],
                                                              sols[t0]['ιCand'], t0, θs[1], εs[1], sols[t1])
                s_, ι_ = float(sSol[0, 0]), float(ιSol[0, 0])
    τAatE = np.array([τAt(stat[np.round(ν[pos], 12)], sE_[pos], ιE_[pos]) for pos in range(T)])
    τLR = np.array([τAt(stat[νLR], sE_[pos], ιE_[pos]) for pos in range(T)])
    τSS = np.array([stat[np.round(ν[pos], 12)]['τStar'] for pos in range(T)])

    # exact economic equilibrium at the approximate tax path, for allocations in the paper's units
    θ, ε, s0 = m.db['θ'].values, m.db['eps'].values, float(sE_[0])
    if isLOG:
        repA = m.EE_report(m.EE_LOG_solve(τA, θ, ε, s0), τA, θ, ε, s0)
    else:
        pathE = exact['path']                    # the exact solve's own warm start layout (solvePEE_CRRA)
        x0 = np.concatenate([pathE['Γs'], pathE['h'], pathE['s']])
        repA = m.EE_report(m.EE_CRRA_solve(τA, θ, ε, s0, x0 = x0), τA, θ, ε, s0)

    def srOverY(r):
        return np.array([float(m.B.savingsRate(r['s'].xs(t), r['s_'].xs(t), r['h'].xs(t), t)) for t in tIdx])
    srE, srA = srOverY(rep), srOverY(repA)
    dates = list(testmod.dates)
    # the tax candidates' counts over every political solve this rho's rows draw on: the exact recursion
    # and the stationary one at every nu (the walk passes through all of them); non-key columns
    counts = multiplicityColumns(exact.get('multiplicity'), *(st['mult'] for st in stat.values()))
    print('   tax candidates (exact and every stationary recursion): ' + ', '.join(
        '{}={}'.format(k, v) for k, v in counts.items()))
    rows = []
    for pos, t in enumerate(tIdx):
        νt = np.round(ν[pos], 12)
        rows.append({'ρ': ρ, 'preferences': prefs, 't': int(t), 'date': int(dates[pos]), 'ν': ν[pos],
                     'τ_exact': τE[pos], 'τ_stat': τA[pos], 'τ_statAtExactState': τAatE[pos],
                     'τ_ssPEE': τSS[pos], 'τ_statLongRun': τLR[pos],
                     'gap_pp': 100*(τA[pos] - τE[pos]), 'gapAtExactState_pp': 100*(τAatE[pos] - τE[pos]),
                     'gapSS_pp': 100*(τSS[pos] - τE[pos]), 'gapLongRun_pp': 100*(τLR[pos] - τE[pos]),
                     's_exact': sE_[pos], 's_stat': sA_[pos], 'ι_exact': ιE_[pos], 'ι_stat': ιA_[pos],
                     'srY_exact': srE[pos], 'srY_stat': srA[pos],
                     'statConv': stat[νt]['conv']} | counts)
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description = 'Stationary-policy approximation vs exact date-specific policies (Argentina).')
    p.add_argument('--rho', type = float, nargs = '*', default = [0.5, 1.0, 2.0])
    p.add_argument('--nι', '--niota', dest = 'nι', type = int, default = 45)
    p.add_argument('--ns', type = int, default = 45)
    p.add_argument('--interpKind', default = 'cubic', choices = ('linear', 'cubic', 'pchip'))
    p.add_argument('--smoothKnots', type = int, default = 4)
    p.add_argument('--pkldir', default = PKLDIR, help = 'directory of pickled calibrated instances')
    p.add_argument('--out', default = 'ARG_stationaryApprox.csv', help = 'file name under results/numerical')
    a = p.parse_args()
    settings = {'nι': a.nι, 'ns': a.ns, 'interpKind': a.interpKind, 'smoothKnots': a.smoothKnots}
    pkldir = a.pkldir if os.path.isabs(a.pkldir) else os.path.join(REPO, a.pkldir)   # cwd is HERE by now
    dfs = [runRho(ρ, settings, pkldir) for ρ in a.rho]
    df = pd.concat(dfs, ignore_index = True)
    os.makedirs(OUTDIR, exist_ok = True)
    out = os.path.join(OUTDIR, a.out)
    df.to_csv(out, index = False)
    pd.set_option('display.width', 220)
    print('\n' + df[['ρ', 'date', 'ν', 'τ_exact', 'τ_stat', 'τ_ssPEE', 'gap_pp', 'gapAtExactState_pp', 'gapSS_pp',
                     'gapLongRun_pp', 'ι_exact', 'ι_stat', 'srY_exact', 'srY_stat', 'statConv']].round(5).to_string(index = False))
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()

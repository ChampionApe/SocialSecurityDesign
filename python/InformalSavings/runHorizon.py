r""" The horizon, the terminal condition and the accuracy of the tax policy functions, Argentina (appendix D,
app:T).

Run:  .venv\Scripts\python.exe python\InformalSavings\runHorizon.py --stage horizon --rho 0.5 1 2
      ... --stage residuals --rho 0.5 1 2 [--refine 4]

HORIZON. The pickled calibrated instances of the headline variant (results/calibration/instancesCommonX/,
common X, pre-reform eps of the workbook), re-solved on the settings they were calibrated at
(shockUniversal.loadCalibrated) with the backward induction of app:T started at T = 2130, 2190 (the paper's:
test.py's t_ss = 3), 2280, and under the steady-state terminal condition: nu constant from the last
projection date (2100) on and the policy functions stationary from T on, computed as the finite-horizon
solution with --kSS periods appended (statGap: max |tau_2190(s, iota) - tau_2220(s, iota)| over the grid nodes
inside the range the path visits, padded 20%: the stationarity of the policy at 2190). The horizon is changed
on the instance itself (resize): every time-indexed parameter is extended by, or truncated to, its constant
tail, which the function asserts; each run first checks that the resize equals a model built afresh at that
horizon from the workbook and the calibrated (beta, omega, eta_0, X_0, X). Read at every dated period
(1920-2100); the paper reads 2010 and 2040:

  baseline  tau, savings over GDP s/(Y + Y0) (GDP formal plus informal output, model.taxTarget's), the
            informal savings ratio iota.
  reform    the universal pension of section 5 (shockUniversal.py's 'match' against type 1: the copy from
            2010, eps raised to eps^U, seeded with the baseline's state and the retirees' inherited shares),
            reported as the change against the baseline of the same horizon.

The 2190 rows must reproduce universal_match_rho*_commonX.csv's tau_2010 baseline and reform (1e-8).

RESIDUALS. As python/US/runHorizon.py's: at every period of the 2190 horizon the period's own selection is
re-run, without smoothing, at the nodes of the state grid (iota_{t-1} under LOG; (s_{t-1}, iota_{t-1}) under
CRRA), at the midpoints of its cells (geometric along iota, whose grid is geometric) and at the state the path
enters t with, against the solved continuation (LOG/CRRA.solveBackward_t, solveTerminal, with the candidate
grids the recursion used). The residual is tau_FOC(x) - tau_hat_t(x), the root of the first order condition
(eq:zdecomposition, eq:PEEcrraTerms) at x against the policy function the path uses, in units of the tax.
'' is the condition as the solver evaluates it (state interpolation and policy smoothing); 'Ref' repeats it on
a tax grid refined --refine-fold, the residual against the condition itself: the tau-derivatives
(policy._gradProfile) are already an interpolating spline in ln(1 - tau), so only the grid is refined. Under
CRRA a point is 'node' (both coordinates nodes), 'cell' (both midpoints) or 'edge' (one of each). Summary
statistics of |residual| per (rho, period, set, region), region 'all' or 'visited'.

Writes results/numerical/ARG_horizon_commonX.csv (key rho, terminal, tss, date) and
ARG_residuals_commonX.csv. Every horizon row carries the tax counts of the baseline and reform solves
(policy.multiplicityColumns).
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
from model import ModelInformalSavings
from shockUniversal import loadCalibrated, solvePEE, universalEps, installEps
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'numerical')
PKLDIR = os.path.join(REPO, 'results', 'calibration', 'instancesCommonX')
SHOCKCSV = os.path.join(REPO, 'results', 'shocks', 'universal_match_rho{:.4f}_commonX.csv')
SETTINGS = {'interpKind': 'cubic', 'smoothKnots': 4, 'nι': 45, 'ns': 45}   # config.ARG['gridSettings']
BASE_TSS = 3                    # test.py's t_ss: the paper's horizon, T = 2190
HORIZONS = (1, 3, 6)            # T = 2130, 2190, 2280
KEY = ['ρ', 'terminal', 'tss', 'date']


# ------------------------------------------------------------------------------------------------- horizon
def horizonList(kSS):
    """ The paper's horizon first (its self-check), the other two, the steady-state terminal condition. """
    return [('finite', BASE_TSS)] + [('finite', k) for k in HORIZONS if k != BASE_TSS] + [('steadyState', kSS)]


def horizonDates(tss):
    return testmod.datesLog.union(pd.Index([testmod.datesLog[-1] + 30*i for i in range(1, tss + 1)]))


def resize(m, T):
    """ A deep copy of m on a T-period horizon: every time-indexed db entry extended by its last value or
    truncated, the indices rebuilt, the warm starts cleared. Exact only where every time-indexed parameter is
    constant over its last three periods (so that x, x[t-1] and x[t+1] all extend or truncate to their own
    values); raises otherwise. """
    mT = deepcopy(m)
    db, T0 = mT.db, m.T
    tNew = pd.Index(range(T), name = 't')
    for k, v in list(db.items()):
        if isinstance(v, (pd.Series, pd.DataFrame)) and isinstance(v.index, pd.MultiIndex):
            raise NotImplementedError(f'resize: {k} is indexed by a MultiIndex.')
        if isinstance(v, (pd.Series, pd.DataFrame)) and v.index.name == 't':
            a = np.asarray(v, dtype = float)
        elif isinstance(v, np.ndarray) and v.shape[:1] == (T0,):
            a = v.astype(float)
        else:
            continue
        if not (np.array_equal(a[-1], a[-2], equal_nan = True) and np.array_equal(a[-2], a[-3], equal_nan = True)):
            raise ValueError(f'resize: {k} is not constant over its last three periods.')
        if isinstance(v, np.ndarray):
            db[k] = np.concatenate([a[:T], np.repeat(a[-1:], max(T - T0, 0), axis = 0)])
        else:
            db[k] = v.reindex(tNew).ffill()
    db['t'], db['txE'] = tNew, pd.Index(range(T - 1), name = 't')
    db['tj'] = pd.MultiIndex.from_product([db['t'], db['j']])
    db['ti'] = pd.MultiIndex.from_product([db['t'], db['i']])
    mT.T = mT.LOG.T = mT.CRRA.T = T
    mT.x0 = {}
    mT.addNamespaces()
    return mT


def fresh(m, T):
    """ The model built afresh at horizon T from the workbook (test.py) and m's calibrated (beta, omega,
    eta_0, X_0), plus X under common X, as calibrate leaves them (_calPostRoot). The check on resize. """
    tss = T - testmod.TLog
    ν = np.hstack([testmod.νLog, np.full(tss, testmod.νLog[-1])])
    t0 = m.db['t'][m.db['t0']]
    pars = dict(testmod.pars) | {'ρ': float(m.db['ρ'].xs(t0)), 'ν': ν}
    mf = ModelInformalSavings(pars = pars, T = T, nj = testmod.kwargs['nj'], commonX = m.commonX, epsPre0 = m.epsPre0)
    if m.commonX:
        mf.initProductivity_commonX(X = float(m.db['Xj'].values[0, 1]), keepInformal = True)
        mf.updateAuxPars()
    mf._calSetPars(m.calibrationPars)
    return mf


def dbGap(a, b):
    """ max relative difference over the time-indexed db entries two models share. """
    gap = 0.0
    for k, v in a.db.items():
        w = b.db.get(k)
        if isinstance(v, (pd.Series, pd.DataFrame, np.ndarray)) and w is not None and not isinstance(v, pd.Index):
            x, y = np.asarray(v, dtype = float), np.asarray(w, dtype = float)
            if x.shape != y.shape:
                return np.inf
            with np.errstate(divide = 'ignore', invalid = 'ignore'):
                d = np.abs(x - y)/np.maximum(np.abs(x), 1e-300)
            gap = max(gap, float(np.nanmax(np.where(x == y, 0, d))) if d.size else 0.0)
    return gap


def gdpSavings(m, report):
    """ s_t/(Y_t + Y0_t) per period (NaN where s_t is not positive: the terminal period), as
    python/paper/datasets.savingsRatePath. """
    s = report['s'].reindex(m.db['t']).values.astype(float)
    Y, Y0 = report['Y'].values.astype(float), report['Y0'].values.astype(float)
    with np.errstate(divide = 'ignore', invalid = 'ignore'):
        sr = s/(Y + Y0)
    return np.where(s > 0, sr, np.nan)


def statGap(sols, tIdx, pos, path, prefs):
    """ max |tau_t - tau_{t+1}| at the period in position pos over the grid nodes inside the visited range
    (padded 20%). """
    if pos + 1 >= len(tIdx):
        return np.nan
    a, b = sols[tIdx[pos]]['τ'], sols[tIdx[pos + 1]]['τ']
    ιv = np.asarray(path['ι_'], dtype = float)
    if prefs == 'LOG':
        ι = a.index.values
        vis = (ι >= 0.8*ιv.min()) & (ι <= 1.2*ιv.max())
    else:
        sv = np.asarray(path['s_'], dtype = float)
        s, ι = a.index.values, a.columns.values
        vis = (((s >= 0.8*sv.min()) & (s <= 1.2*sv.max()))[:, None]
               & ((ι >= 0.8*ιv.min()) & (ι <= 1.2*ιv.max()))[None, :])
    return float(np.nanmax(np.abs(a.values - b.values)[vis]))


def mergeWrite(path, rows, key):
    new = pd.DataFrame(rows)
    if os.path.exists(path):
        old = pd.read_csv(path)
        k = lambda df: df[key].astype(str).agg('|'.join, axis = 1)
        old = old[~k(old).isin(set(k(new)))]
        new = pd.concat([old, new], ignore_index = True)
    os.makedirs(os.path.dirname(path), exist_ok = True)
    new.to_csv(path, index = False)


def stageHorizon(ρs, kSS, out):
    for ρ in ρs:
        m0, prefs = loadCalibrated(ρ, SETTINGS, pkldir = PKLDIR, epsPre0 = False)
        posT0 = int(m0.db['t0'])
        ref = pd.read_csv(SHOCKCSV.format(ρ)).set_index('t')
        rows = []
        for terminal, tss in horizonList(kSS):
            tic = time.time()
            T = testmod.TLog + tss
            m = resize(m0, T)
            gap = dbGap(m, fresh(m0, T))
            if gap > 1e-12:
                raise RuntimeError(f'ρ={ρ}, T={T}: the resized instance differs from a fresh build by {gap:.1e}.')
            tIdx, t0 = m.db['t'], m.db['t'][posT0]
            base = solvePEE(m, prefs)
            τB = base['τ'].values
            seed = m.stateAtT0(base['report'], t0, init = base['init'])
            mt0 = m.createCopyFromt0(t0)
            installEps(mt0, universalEps(mt0, rule = 'match', refType = 1))
            shock = solvePEE(mt0, prefs, **seed)
            τR = np.concatenate([np.full(posT0, np.nan), shock['τ'].values])
            srB, srR = gdpSavings(m, base['report']), np.concatenate([np.full(posT0, np.nan),
                                                                     gdpSavings(mt0, shock['report'])])
            ιB = base['report']['ι'].reindex(tIdx).values.astype(float)
            ιR = np.concatenate([np.full(posT0, np.nan), shock['report']['ι'].reindex(mt0.db['t']).values.astype(float)])
            if (terminal, tss) == ('finite', BASE_TSS):
                d = max(abs(τB[posT0] - float(ref.loc[t0, 'τ_base'])), abs(τR[posT0] - float(ref.loc[t0, 'τ_reform'])))
                print(f'ρ={ρ} ({prefs}) self-check against {os.path.basename(SHOCKCSV.format(ρ))}: τ_2010 base '
                      f'{τB[posT0]:.8f} vs {float(ref.loc[t0, "τ_base"]):.8f}, reform {τR[posT0]:.8f} vs '
                      f'{float(ref.loc[t0, "τ_reform"]):.8f}: max diff {d:.1e}')
                if d > 1e-8:
                    raise RuntimeError(f'ρ={ρ}: the 2190 horizon does not reproduce the published reform ({d:.1e}).')
            sg = statGap(base['sols'], tIdx, testmod.TLog + BASE_TSS - 1, base['path'], prefs) if terminal == 'steadyState' else np.nan
            dates = horizonDates(tss)
            counts = multiplicityColumns(base.get('multiplicity'), shock.get('multiplicity'))
            for pos in range(testmod.TLog):
                rows.append({'ρ': ρ, 'preferences': prefs, 'terminal': terminal, 'tss': tss, 'T': int(dates[-1]),
                             'pos': pos, 'date': int(dates[pos]), 'ν': float(m.db['ν'].iloc[pos]),
                             'τ': τB[pos], 'srGDP': srB[pos], 'ι': ιB[pos],
                             'τ_reform': τR[pos], 'srGDP_reform': srR[pos], 'ι_reform': ιR[pos],
                             'dτ': τR[pos] - τB[pos], 'dsrGDP': srR[pos] - srB[pos],
                             'statGap': sg, 'seconds': time.time() - tic} | counts)
            print(f'  ρ={ρ} {terminal:<11} T={int(dates[-1])}: τ_2010={τB[posT0]:.6f} τ_2040={τB[posT0+1]:.6f} '
                  f'Δτ_2010={100*(τR[posT0]-τB[posT0]):+.4f} p.p. Δτ_2040={100*(τR[posT0+1]-τB[posT0+1]):+.4f} p.p. '
                  f'statGap={sg:.1e} resize-vs-fresh {gap:.0e}  ({time.time()-tic:.0f}s)')
            mergeWrite(out, rows, KEY)
    return out


# ------------------------------------------------------------------------------------------------- residuals
def refined(τGrid, k):
    a, b = τGrid[:-1], τGrid[1:]
    return np.sort(np.concatenate([τGrid] + [a + (b - a)*j/k for j in range(1, k)]))


def stats(r):
    a = np.abs(np.asarray(r, dtype = float).ravel())
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {'n': 0, 'max': np.nan, 'p99': np.nan, 'median': np.nan}
    return {'n': int(a.size), 'max': float(a.max()), 'p99': float(np.percentile(a, 99)), 'median': float(np.median(a))}


def stageResiduals(ρs, refine, out):
    rows = []
    for ρ in ρs:
        tic = time.time()
        m, prefs = loadCalibrated(ρ, SETTINGS, pkldir = PKLDIR, epsPre0 = False)
        pol = getattr(m, prefs)
        base = solvePEE(m, prefs)
        sols, tIdx, path = base['sols'], m.db['t'], base['path']
        θ, ε = m.db['θ'].values, m.db['eps'].values
        l, u = pol.GS['PEE']['gridSettings']['l'], pol.GS['PEE']['gridSettings']['u']
        τGrid = pol.GS['PEE']['solGrids']['τ'].copy()
        ιGrid = np.asarray(sols[tIdx[0]]['ι_'], dtype = float)
        ιCells = np.sqrt(ιGrid[:-1]*ιGrid[1:])
        ιPath = np.asarray(path['ι_'], dtype = float)
        ιVis = (0.8*ιPath.min(), 1.2*ιPath.max())
        if prefs == 'CRRA':
            sGrid = np.asarray(sols[tIdx[0]]['s_'], dtype = float)
            sCells = 0.5*(sGrid[:-1] + sGrid[1:])
            sPath = np.asarray(path['s_'], dtype = float)
            sVis = (0.8*sPath.min(), 1.2*sPath.max())
        pool = {}
        for pos, t in enumerate(tIdx):
            last = pos == len(tIdx) - 1
            tLag = tIdx[pos - 1] if pos > 0 else m.B.tFirst
            xι = np.unique(np.concatenate([ιGrid, ιCells, [ιPath[pos]]]))
            kι = np.where(np.isin(xι, ιGrid), 'node', np.where(np.isin(xι, ιCells), 'cell', 'path'))
            if prefs == 'LOG':
                τHat = np.clip(sols[t]['τPolicy'](xι), l, u)
                kind, vis = kι, (xι >= ιVis[0]) & (xι <= ιVis[1])
            else:
                xs = np.unique(np.concatenate([sGrid, sCells, [sPath[pos]]]))
                ks = np.where(np.isin(xs, sGrid), 'node', np.where(np.isin(xs, sCells), 'cell', 'path'))
                S, I = np.meshgrid(xs, xι, indexing = 'ij')
                τHat = np.clip(sols[t]['τPolicy'](S, I), l, u)
                KS, KI = np.meshgrid(ks, kι, indexing = 'ij')
                kind = np.where((KS == 'path') & (KI == 'path'), 'path',
                                np.where((KS == 'path') | (KI == 'path'), 'other',
                                         np.where(KS == KI, KS, 'edge')))
                vis = (((xs >= sVis[0]) & (xs <= sVis[1]))[:, None] & ((xι >= ιVis[0]) & (xι <= ιVis[1]))[None, :])
            res = {}
            for tag, grid in (('', τGrid), ('Ref', refined(τGrid, refine))):
                pol.GS['PEE']['solGrids']['τ'] = grid
                try:
                    if prefs == 'LOG':
                        rep = (pol.solveTerminal(θ[pos], ε[pos], t = t, ιGrid = xι) if last else
                               pol.solveBackward_t(sols[tIdx[pos + 1]], t, tLag, θ[pos], θ[pos + 1], ε[pos], ε[pos + 1],
                                                   xι, sols[t]['ιCand'], smooth = 0))
                    else:
                        rep = (pol.solveTerminal(θ[pos], ε[pos], t = t, ιGrid = xι, sGrid = xs) if last else
                               pol.solveBackward_t(sols[tIdx[pos + 1]], t, tLag, θ[pos], θ[pos + 1], ε[pos], ε[pos + 1],
                                                   xs, xι, sols[t]['sCand'], sols[t]['ιCand'], smooth = 0))
                finally:
                    pol.GS['PEE']['solGrids']['τ'] = τGrid
                res[tag] = np.asarray(rep['τ'], dtype = float) - τHat
            date = int(testmod.dates[pos])
            for setName in ('node', 'cell', 'edge', 'path'):
                for region, mask in (('all', np.ones_like(vis)), ('visited', vis)):
                    sel = (kind == setName) & mask
                    if not sel.any() or (setName == 'path' and region == 'all'):
                        continue
                    for tag in ('', 'Ref'):
                        pool.setdefault((setName, region, tag), []).append(res[tag][sel])
                    rows.append({'model': 'ARG', 'preferences': prefs, 'ρ': ρ, 'nTau': len(τGrid), 'refine': refine,
                                 'date': date, 'set': setName, 'region': region}
                                | stats(res[''][sel]) | {f'{k}Ref': v for k, v in stats(res['Ref'][sel]).items()})
            print(f'  ρ={ρ} {date}: max |r| {100*np.nanmax(np.abs(res[""])):.2e} p.p., reference {100*np.nanmax(np.abs(res["Ref"])):.2e}'
                  f'  ({time.time()-tic:.0f}s)')
        for (setName, region, tag), v in pool.items():
            if tag:
                continue
            st, stF = stats(np.concatenate(v)), stats(np.concatenate(pool[(setName, region, 'Ref')]))
            rows.append({'model': 'ARG', 'preferences': prefs, 'ρ': ρ, 'nTau': len(τGrid), 'refine': refine,
                         'date': 'all', 'set': setName, 'region': region} | st | {f'{k}Ref': w for k, w in stF.items()})
            print(f'  ρ={ρ} {setName:<4} {region:<7}: max {100*st["max"]:.2e} p.p., p99 {100*st["p99"]:.2e}; '
                  f'reference FOC max {100*stF["max"]:.2e}, p99 {100*stF["p99"]:.2e}')
        mergeWrite(out, rows, ['model', 'ρ', 'date', 'set', 'region'])
    return out


def main():
    p = argparse.ArgumentParser(description = 'Horizon, terminal condition and FOC residuals (Argentina).')
    p.add_argument('--stage', nargs = '+', default = ['horizon', 'residuals'], choices = ('horizon', 'residuals'))
    p.add_argument('--rho', type = float, nargs = '+', default = [1.0, 0.5, 2.0])
    p.add_argument('--kSS', type = int, default = 40, help = 'periods appended for the steady-state terminal condition')
    p.add_argument('--refine', type = int, default = 4, help = 'residuals: tax-grid refinement factor')
    a = p.parse_args()
    for st in a.stage:
        tic = time.time()
        print(f'\n=== stage {st} ===')
        if st == 'horizon':
            stageHorizon(a.rho, a.kSS, os.path.join(OUTDIR, 'ARG_horizon_commonX.csv'))
        else:
            stageResiduals(a.rho, a.refine, os.path.join(OUTDIR, 'ARG_residuals_commonX.csv'))
        print(f'=== stage {st} done ({time.time()-tic:.0f}s) ===')
    print('DONE')
    return 0


if __name__ == '__main__':
    sys.exit(main())

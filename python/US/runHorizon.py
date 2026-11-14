r""" The horizon, the terminal condition and the accuracy of the tax policy functions (appendix D, app:T).

Run:  .venv\Scripts\python.exe python\US\runHorizon.py --stage horizon --rho 0.5 2     # CRRA, theta given
      ... --stage design --rho 0.5 2          # CRRA, theta chosen (LeadedCRRA2D, alg esc:crra2D), lambda calibrated
      ... --stage log                         # rho = 1: LOG, theta given and chosen
      ... --stage residuals --rho 0.5 2 --ns 150 50
      ... --stage jumps --rho 2 --ns 50 150   # tau_2020 along a fine scan of omega (finding #5's tell)
      ... --stage scheme --rho 0.5 2          # whole solves under converged tau-derivatives
      ... --theta C --main                    # every stage at theta variant C (the paper's), as below

THE THETA VARIANT. --theta VARIANT --main (VARIANT = python/paper/config.py's US['theta'], oecd26.thetaRoot): every
model at theta0 = the variant's U.S. row of data/oecd26_theta.csv (db['θ0'], ModelUS.getθ), every calibration row
read (the sweep, the cost, the design path of the self-checks) kept at that thetaSource and theta0 (_at), every row
written carrying both. A csv holds one theta: a merge drops the rows on disk at another (mergeWrite). Without --theta
the run is refused (every output is under results/numerical/, oecd26.indexTheta); --theta index runs the workbook's
theta there.

HORIZON. The calibrated economy (results/calibration/US_rhoGridCommonX.csv; the cost lambda from
results/esc/escCalibration{,CRRA}.csv) with every parameter held at its calibrated value, re-solved with the
backward induction of app:T started at T = 2140, 2200 (the paper's: test.py's t_ss = 3), 2290, and under the
steady-state terminal condition: nu constant from the last projection date (2110) on and the policy functions
stationary from T on, the infinite-horizon Markov equilibrium of that economy. The last is computed as the
finite-horizon solution with --kSS periods appended (2110 + 30 kSS), its policy at 2200 checked for
stationarity (statGap: max |tau_2200(s) - tau_2230(s)| over the states the path visits, padded 20%). Each
horizon is a fresh model at that T (ModelUS / ModelESC with nu extended by its last value), so nothing but
T differs. Rows are the dated periods 1960-2110; the paper reads 2020 and 2050.

  horizon  theta given (the paper's CRRA rows, rho = 0.5, 2): tau, s/Y, R, h.
  design   theta chosen one period ahead under the paper's cost 'size' (eq esc:AB), the exact recursion of
           runESCcrra.py --exact: (beta, omega, X) recalibrated with the cost on at the 2200 horizon
           (ModelESC.calibrate, as runESCcrra.py does) and installed at every horizon; free path
           (pinAtT0 = False), so theta_2020 is the design in force, chosen in 1990.
  log      rho = 1. The LOG tax condition decouples across dates (eq PEELOG:decoupling), so the tax at every date
           does not depend on the horizon, nor does the design chosen at every date but T-1 (eq
           esc:separability): the rows must agree to solver precision at every date up to 2110, both readings.

RESIDUALS. The first order condition z_t of the CRRA terms (app:T, eq:terminalPEECRRA), at states off the
grid the policy functions were solved on, on the 2200 horizon. At every period t the period's own selection
(the consistent condition's crossings tested at frozen shares, num_robustroot.tex) is re-run, without
smoothing, at the nodes of the state grid S, the midpoints of its cells and the state the equilibrium path
enters t with, against the solved continuation tau_{t+1}(.) (CRRA.solveBackward_t / solveTerminal). The
residual is tau_FOC(s) - tau_hat_t(s): the root of the first order condition at s against the policy
function the recursion and the path use (smoothed, interpolated, clipped to [l, u]), in units of the tax.
Two references: '' the condition as the solver evaluates it (its tax grid of 101 nodes, its tau-derivatives:
CRRA._gradProfile's least-squares cubic in ln(1-tau), knots gradKnotSpacing = 0.04 apart in tau), which isolates
the state interpolation and the policy smoothing; 'Ref' the same derivatives on a tax grid refined --refine-fold
(801 nodes at the default 8, the knots still 0.04 apart: CRRA._gradKnots), where they have converged (to 1.8e-5 in
the 2020 tax of a whole solve from 101 to 801 nodes; an interpolating spline does not converge here, python/US
RESEARCH_LOG.md, 2026-10-24), so 'Ref' is the residual against the first order condition itself. Written as summary
statistics (max, 99th percentile, median of |residual|) per (rho, ns, period, set, region): set 'node'
(smoothing alone), 'cell' (interpolation and smoothing), 'path'; region 'all' or 'visited' (within the
path's range of s, padded 20%).

SCHEME. The whole solve at the calibrated parameters under the reference of RESIDUALS (n = 801, the knots
held) against the published one (n = 101): tau and s/Y at 2020 and 2050, and at theta = 0 and theta = 1 (shocks.py),
whose gap in 2020 is section 6's design effect.

JUMPS. tau_2020 over a fine scan of omega around its calibrated value (beta, X fixed): isolated spikes in the
successive differences are discrete selection changes (finding #5); their size against the residuals says
whether the residuals see them.

Writes results/numerical/US_horizon.csv (stages horizon, log; key block, preferences, rho, terminal, tss,
date), US_horizonDesign.csv (stage design), US_residuals.csv, US_jumps.csv, US_derivativeScheme.csv. Every
horizon row carries the counts of its solve (policy.multiplicityColumns; the design rows also the design
layer's).
"""
import os, sys, argparse, time
from contextlib import contextmanager
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import test as testmod
import shocks as sh
import oecd26
from model import ModelUS
from modelESC import ModelESC
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'numerical')
CALDIR = os.path.join(REPO, 'results', 'calibration')
ESCDIR = os.path.join(REPO, 'results', 'esc')
GS = {'n': 101, 'ns': 150, 'smoothKnots': 4, 'interpKind': 'linear'}   # runShocksUS.py's / runESCcrra.GSC
BASE_TSS = 3                    # test.py's t_ss: the paper's horizon, T = 2200
HORIZONS = (1, 3, 6)            # T = 2140, 2200, 2290
KEY = ['block', 'preferences', 'ρ', 'terminal', 'tss', 'date']
WEDGE = {'spec': 'size', 'phi': 0.5}           # the paper's cost (config.US['esc'])
N_CAND_2D = 41                                 # config.US['esc']['nCand2D']
# The run's theta variant (main's --theta V --main): theta0 imposed on every model, rows read and written at it.
THETA = {'variant': None, 'θ0': None}


def _at(df):
    """ The rows of df at the run's theta (oecd26.rowsAtTheta at THETA's variant and θ0). """
    return oecd26.rowsAtTheta(df, THETA['variant'], THETA['θ0'])


def _tag():
    """ The columns every row written carries under a variant: thetaSource, θ0 ({} without one). """
    return {} if THETA['variant'] is None else {'thetaSource': THETA['variant'], 'θ0': THETA['θ0']}


def setTheta(theta, main):
    """ --theta V --main: V must be config.US['theta'] (results/numerical/ holds the paper's variant only). SystemExit
    for one flag without the other or another variant; the sweep rows are checked where they are read (sweepRow). """
    if theta is None and not main:
        return
    if theta is None or not main:
        raise SystemExit('--theta needs --main here (and --main --theta)')
    try:
        oecd26.thetaRoot(theta, main = True)
        THETA.update({'variant': theta, 'θ0': oecd26.thetaFromFile('USA', theta)})
    except (ValueError, FileNotFoundError) as e:
        raise SystemExit('--theta {}: {}'.format(theta, e))


# ------------------------------------------------------------------------------------------------- set-up
def horizon(tss):
    """ (T, nu (T,), dates (T,)) of the horizon with tss periods appended after the last projection date. """
    T = testmod.TLog + tss
    ν = np.hstack([testmod.νLog, np.full(T - testmod.TLog, testmod.νLog[-1])])
    dates = testmod.datesLog.union(pd.Index([testmod.datesLog[-1] + 30*i for i in range(1, tss + 1)]))
    return T, ν, dates


def sweepRow(ρ):
    """ The common-X U.S. sweep's last row at rho, at the run's theta (_at). """
    df = _at(pd.read_csv(os.path.join(CALDIR, 'US_rhoGridCommonX.csv')))
    hit = df.loc[(df['ρ'] - ρ).abs() < 1e-9]
    if hit.empty:
        raise KeyError(f'No common-X US calibration at ρ={ρ} at θ {THETA}.')
    return hit.iloc[-1]


def build(ρ, β, ω, X, tss, gs = GS, cls = ModelUS, **kw):
    """ The common-X US model at (rho, beta, omega, X) on the horizon tss, every other input the workbook's. """
    T, ν, dates = horizon(tss)
    pars = (dict(testmod.pars) | {'ρ': float(ρ), 'β': float(β), 'ω': float(ω), 'ν': ν}
            | ({} if THETA['θ0'] is None else {'θ0': THETA['θ0']}))
    m = cls(pars = pars, commonX = True, T = T, nj = testmod.kwargs['nj'], **kw)
    m.db['dates'], m.db['workweek'] = dates, testmod.workweek
    m.CRRA.initGS(gs)
    m.LOG.initGS({k: v for k, v in gs.items() if k != 'ns'})
    m.initProductivity_commonX(X = float(X))
    m.updateAuxPars()
    return m


def commonX(m):
    """ The calibrated scalar X of a common-X model (Xj's formal columns). """
    return float(m.db['Xj'].values[0, 1])


def readRows(m, τ, report, θ, ids, extra = None):
    """ One row per dated period (1960 to the last projection date): tau, theta, s/Y (the paper's unit), s/(wh),
    R, h. """
    rows = []
    for pos in range(testmod.TLog):
        t = m.db['t'][pos]
        s, s_, h = report['s'].xs(t), report['s_'].xs(t), report['h'].xs(t)
        sY = float(m.B.savingsRate(s, s_, h, t))
        α = float(m.db['α'].xs(t))
        rows.append(ids | {'pos': pos, 'date': int(m.db['dates'][pos]), 'ν': float(m.db['ν'].xs(t)),
                           'τ': float(np.asarray(τ, dtype = float)[pos]), 'θ': float(np.asarray(θ, dtype = float)[pos]),
                           'srOverY': sY, 'sr': sY/(1-α), 'R': float(report['R'].xs(t)), 'h': float(h)}
                    | (extra or {}) | _tag())
    return rows


def statGap(m, sols, report, pos):
    """ max |tau_t(s) - tau_{t+1}(s)| at the period in position pos, over the grid nodes inside the range of s
    the path visits (padded 20%): the stationarity of the policy there. NaN if pos+1 is not on the horizon. """
    tIdx = m.db['t']
    if pos + 1 >= len(tIdx):
        return np.nan
    s_ = report['s_'].values.astype(float)
    a, b = sols[tIdx[pos]]['τ'], sols[tIdx[pos + 1]]['τ']
    x = a.index.values
    vis = (x >= 0.8*s_.min()) & (x <= 1.2*s_.max())
    return float(np.nanmax(np.abs(a.values - b.values)[vis]))


@contextmanager
def _locked(path, wait = 600.):
    """ An exclusive lock file beside `path` (path + '.lock'), so that processes merging into one csv (the design
    stage per rho) never interleave a read and a write. RuntimeError after `wait` seconds. """
    lock, tic = path + '.lock', time.time()
    while True:
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            break
        except FileExistsError:
            if time.time() - tic > wait:
                raise RuntimeError(f'{lock} held for {wait:.0f}s: remove it if no run is writing {path}')
            time.sleep(0.5)
    try:
        yield
    finally:
        os.close(fd)
        os.remove(lock)


def mergeWrite(path, rows, key):
    """ Replace the rows of `path` whose key is among `rows`' and write; rows on disk at another theta than the
    run's (_at) are dropped, so a csv holds one theta. """
    new = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(path), exist_ok = True)
    with _locked(path):
        if os.path.exists(path):
            old = pd.read_csv(path)
            keep = _at(old)
            if len(keep) < len(old):
                print(f'  {os.path.relpath(path, REPO)}: {len(old) - len(keep)} rows at another θ dropped')
            k = lambda df: df[key].astype(str).agg('|'.join, axis = 1)
            new = pd.concat([keep[~k(keep).isin(set(k(new)))], new], ignore_index = True)
        new.to_csv(path, index = False)


def horizonList(kSS):
    """ [(terminal, tss)]: the paper's horizon first (its self-check stops a run that is not the published
    solve before anything else is spent), the other two finite horizons, the steady-state terminal condition. """
    return [('finite', BASE_TSS)] + [('finite', k) for k in HORIZONS if k != BASE_TSS] + [('steadyState', kSS)]


# ------------------------------------------------------------------------------------------------- stages
def stageHorizon(ρs, kSS, out, gs = GS):
    """ theta given, CRRA. The 2200 row must reproduce the sweep's tau_2020 (1e-6, runShocksUS.py's check). """
    for ρ in ρs:
        row = sweepRow(ρ)
        rows = []
        for terminal, tss in horizonList(kSS):
            tic = time.time()
            m = build(ρ, row['β'], row['ω'], row['X'], tss, gs)
            outp = m.solvePEE_CRRA()
            τ = outp['τ'].values
            if terminal == 'finite' and tss == BASE_TSS and abs(τ[m.db['t0']] - float(row['τ'])) > 1e-6:
                raise RuntimeError(f'ρ={ρ}: the 2200 horizon gives τ_2020={τ[m.db["t0"]]:.8f}, the sweep '
                                   f'recorded {float(row["τ"]):.8f}: not the published solve.')
            gap = statGap(m, outp['sols'], outp['report'], testmod.TLog + BASE_TSS - 1) if terminal == 'steadyState' else np.nan
            ids = {'block': 'tax', 'preferences': 'CRRA', 'ρ': ρ, 'terminal': terminal, 'tss': tss,
                   'T': int(m.db['dates'][-1]), 'ns': gs['ns']}
            rows += readRows(m, τ, outp['report'], m.db['θ'].values, ids,
                             {'statGap': gap, 'seconds': time.time() - tic} | multiplicityColumns(outp['multiplicity']))
            print(f'  ρ={ρ} {terminal:<11} T={ids["T"]}: τ_2020={τ[2]:.6f} τ_2050={τ[3]:.6f}  '
                  f'statGap={gap:.1e}  ({time.time()-tic:.0f}s)')
            mergeWrite(out, rows, KEY)
    return out


def escLambda(ρ, method):
    """ The converged calibrated cost lambda of the paper's spec, common X: LOG from escCalibration.csv, the exact
    CRRA recursion's from escCalibrationCRRA.csv (method = 'exact'). """
    if method == 'LOG':
        df = _at(pd.read_csv(os.path.join(ESCDIR, 'escCalibration.csv')))
        hit = df[(df['spec'] == WEDGE['spec']) & (df['phi'] == WEDGE['phi']) & df['commonX'].astype(bool)
                 & df['converged'].astype(bool)]
    else:
        df = _at(pd.read_csv(os.path.join(ESCDIR, 'escCalibrationCRRA.csv')))
        hit = df[(df['spec'] == WEDGE['spec']) & (df['phi'] == WEDGE['phi']) & df['commonX'].astype(bool)
                 & df['converged'].astype(bool) & ((df['ρ'] - ρ).abs() < 1e-9) & (df['method'] == method)]
    if hit.empty:
        raise KeyError(f'No {method} calibration of {WEDGE} at ρ={ρ} at θ {THETA}.')
    return hit.iloc[-1]


def escPathRow(ρ, method):
    """ The published design path's 2020 row (escPath.csv / escPathCRRA.csv), for the self-check. """
    f = 'escPath.csv' if method == 'LOG' else 'escPathCRRA.csv'
    df = _at(pd.read_csv(os.path.join(ESCDIR, f)))
    hit = df[(df['spec'] == WEDGE['spec']) & (df['phi'] == WEDGE['phi']) & df['commonX'].astype(bool) & (df['pos'] == 2)]
    if method != 'LOG':
        hit = hit[((hit['ρ'] - ρ).abs() < 1e-9) & (hit['method'] == method)]
    return hit.iloc[-1]


def stageDesign(ρs, kSS, out, gs = GS, only = None):
    """ theta chosen, CRRA, the exact 2-D recursion. Calibrated at 2200 with the cost on, then every horizon at
    those parameters. The 2200 row must reproduce escPathCRRA.csv's exact theta_2020 and tau_2020. """
    for ρ in ρs:
        cal = escLambda(ρ, 'exact')
        λ = float(cal['p'])
        kw = dict(wedge = WEDGE | {'p': λ}, nθCandCRRA = 13, nθCand2D = N_CAND_2D, designRule = 'root', Ma = 5)
        tic = time.time()
        mB = build(ρ, cal['β'], cal['ω'], sweepRow(ρ)['X'], BASE_TSS, gs, cls = ModelESC, **kw)
        mB.calibrate()
        β, ω, X = mB.simpleβinv(), float(mB.db['ω'].iloc[0]), commonX(mB)
        print(f'ρ={ρ}: λ={λ:.6f}, recalibrated with the cost on: β={β:.6f} ω={ω:.6f} X={X:.6f} '
              f'(csv β={float(cal["β"]):.6f} ω={float(cal["ω"]):.6f})  ({time.time()-tic:.0f}s)')
        ref = escPathRow(ρ, 'exact')
        rows = []
        for terminal, tss in horizonList(kSS):
            if only and terminal not in only and (terminal, tss) != ('finite', BASE_TSS):
                continue                                   # the 2200 path always runs: it is the self-check
            tic = time.time()
            m = mB if (terminal, tss) == ('finite', BASE_TSS) else build(ρ, β, ω, X, tss, gs, cls = ModelESC, **kw)
            led = m.solveLeaded2D(pinAtT0 = False)
            τ, θ = led['τ'].values, led['θ'].values
            if (terminal, tss) == ('finite', BASE_TSS):
                d = max(abs(θ[2] - float(ref['θ'])), abs(τ[2] - float(ref['τ'])))
                print(f'  self-check against escPathCRRA.csv (exact): θ_2020 {θ[2]:.6f} vs {float(ref["θ"]):.6f}, '
                      f'τ_2020 {τ[2]:.6f} vs {float(ref["τ"]):.6f}: max diff {d:.1e}')
                if d > 1e-4:
                    raise RuntimeError(f'ρ={ρ}: the 2200 design path does not reproduce the published one ({d:.1e}).')
            ids = {'block': 'design', 'preferences': 'CRRA', 'ρ': ρ, 'terminal': terminal, 'tss': tss,
                   'T': int(m.db['dates'][-1]), 'ns': gs['ns'], 'spec': WEDGE['spec'], 'p': λ}
            rows += readRows(m, τ, led['report'], θ, ids,
                             {'seconds': time.time() - tic} | multiplicityColumns(led['multiplicity'], design = True))
            print(f'  ρ={ρ} {terminal:<11} T={ids["T"]}: θ_2020={θ[2]:.6f} θ_2050={θ[3]:.6f} τ_2020={τ[2]:.6f} '
                  f'τ_2050={τ[3]:.6f}  ({time.time()-tic:.0f}s)')
            mergeWrite(out, rows, KEY)
    return out


def stageLog(kSS, out, gs = GS):
    """ rho = 1. theta given (solvePEE_LOG) and chosen (LeadedLOG, the paper's cost at its LOG lambda). """
    row = sweepRow(1.0)
    cal = escLambda(1.0, 'LOG')
    λ = float(cal['p'])
    mB = build(1.0, cal['β'], cal['ω'], row['X'], BASE_TSS, gs, cls = ModelESC, wedge = WEDGE | {'p': λ})
    mB.calibrate()
    β, ω, X = mB.simpleβinv(), float(mB.db['ω'].iloc[0]), commonX(mB)
    ref = escPathRow(1.0, 'LOG')
    rows = []
    for terminal, tss in horizonList(kSS):
        tic = time.time()
        m = build(1.0, row['β'], row['ω'], row['X'], tss, gs)
        outp = m.solvePEE_LOG()
        τ = outp['τ'].values
        if (terminal, tss) == ('finite', BASE_TSS) and abs(τ[2] - float(row['τ'])) > 1e-6:
            raise RuntimeError(f'ρ=1: τ_2020={τ[2]:.8f}, the sweep recorded {float(row["τ"]):.8f}.')
        ids = {'block': 'tax', 'preferences': 'LOG', 'ρ': 1.0, 'terminal': terminal, 'tss': tss,
               'T': int(m.db['dates'][-1])}
        rows += readRows(m, τ, outp['report'], m.db['θ'].values, ids,
                         {'seconds': time.time() - tic} | multiplicityColumns(outp['multiplicity']))
        tic = time.time()
        mE = mB if (terminal, tss) == ('finite', BASE_TSS) else build(1.0, β, ω, X, tss, gs, cls = ModelESC,
                                                                      wedge = WEDGE | {'p': λ})
        led = mE.solveLeaded(pinAtT0 = False)
        τE, θE = led['τ'].values, led['θ'].values
        if (terminal, tss) == ('finite', BASE_TSS):
            d = max(abs(θE[2] - float(ref['θ'])), abs(τE[2] - float(ref['τ'])))
            print(f'  self-check against escPath.csv: θ_2020 {θE[2]:.6f} vs {float(ref["θ"]):.6f}, τ_2020 '
                  f'{τE[2]:.6f} vs {float(ref["τ"]):.6f}: max diff {d:.1e}')
            if d > 1e-6:
                raise RuntimeError(f'ρ=1: the 2200 LOG design path does not reproduce escPath.csv ({d:.1e}).')
        ids = ids | {'block': 'design', 'spec': WEDGE['spec'], 'p': λ}
        rows += readRows(mE, τE, led['report'], θE, ids,
                         {'seconds': time.time() - tic} | multiplicityColumns(led['multiplicity']))
        print(f'  ρ=1 {terminal:<11} T={ids["T"]}: τ_2020={τ[2]:.8f} τ_2050={τ[3]:.8f} | chosen θ_2020={θE[2]:.6f} '
              f'θ_2050={θE[3]:.6f} θ_2110={θE[5]:.6f}')
        mergeWrite(out, rows, KEY)
    return out


def refined(τGrid, k):
    """ The tax grid with k-1 points inserted in every cell (its nodes kept). """
    a, b = τGrid[:-1], τGrid[1:]
    inner = [a + (b - a)*j/k for j in range(1, k)]
    return np.sort(np.concatenate([τGrid] + inner))


@contextmanager
def referenceDerivatives(m, k, spacing = None):
    """ The CRRA tax condition with converged tau-derivatives: m's tax grid refined k-fold with the derivatives'
    knots held `spacing` apart in tau (gridSettings' gradKnotSpacing, the production 0.04 by default), so that
    CRRA._gradKnots puts k times as many nodes under each knot; both restored on exit. """
    settings = m.CRRA.GS['PEE']['gridSettings']
    grid, spacing0 = m.CRRA.GS['PEE']['solGrids']['τ'].copy(), settings.get('gradKnotSpacing')
    m.CRRA.GS['PEE']['solGrids']['τ'] = refined(grid, k)
    if spacing is not None:
        settings['gradKnotSpacing'] = float(spacing)
    try:
        yield
    finally:
        m.CRRA.GS['PEE']['solGrids']['τ'] = grid
        if spacing0 is None:
            settings.pop('gradKnotSpacing', None)
        else:
            settings['gradKnotSpacing'] = spacing0


def stats(r):
    a = np.abs(np.asarray(r, dtype = float))
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {'n': 0, 'max': np.nan, 'p99': np.nan, 'median': np.nan}
    return {'n': int(a.size), 'max': float(a.max()), 'p99': float(np.percentile(a, 99)), 'median': float(np.median(a))}


def stageResiduals(ρs, nss, refine, out):
    """ The tax policy functions' first-order-condition residuals off the grid, 2200 horizon. """
    rows = []
    for ρ in ρs:
        row = sweepRow(ρ)
        for ns in nss:
            tic = time.time()
            gs = GS | {'ns': ns}
            m = build(ρ, row['β'], row['ω'], row['X'], BASE_TSS, gs)
            outp = m.solvePEE_CRRA()
            sols, tIdx = outp['sols'], m.db['t']
            θ, ε = m.db['θ'].values, m.db['eps'].values
            l, u = m.CRRA.GS['PEE']['gridSettings']['l'], m.CRRA.GS['PEE']['gridSettings']['u']
            sGrid = sols[tIdx[0]]['τ'].index.values.astype(float)
            sPath = outp['report']['s_'].values.astype(float)
            lo, hi = 0.8*sPath.min(), 1.2*sPath.max()
            τGrid = m.CRRA.GS['PEE']['solGrids']['τ'].copy()
            pool = {}
            for pos, t in enumerate(tIdx):
                cells = 0.5*(sGrid[:-1] + sGrid[1:])
                x = np.unique(np.concatenate([sGrid, cells, [sPath[pos]]]))
                kind = np.where(np.isin(x, sGrid), 'node', np.where(np.isin(x, cells), 'cell', 'path'))
                τHat = np.clip(sols[t]['τPolicy'](x), l, u)

                def rootAt():
                    if pos == len(tIdx) - 1:
                        return m.CRRA.solveTerminal(θ[pos], ε[pos], t = t, sGrid = x)['τ'].values.astype(float)
                    tLag = tIdx[pos - 1] if pos > 0 else m.B.tFirst
                    return m.CRRA.solveBackward_t(sols[tIdx[pos + 1]], t, tLag, tIdx[pos + 1], θ[pos], θ[pos + 1],
                                                  ε[pos], ε[pos + 1], x, sGrid, smooth = 0)['τ'].values.astype(float)
                res = {'': rootAt() - τHat}
                with referenceDerivatives(m, refine):
                    res['Ref'] = rootAt() - τHat
                vis = (x >= lo) & (x <= hi)
                date = int(m.db['dates'][pos])
                for setName in ('node', 'cell', 'path'):
                    for region, mask in (('all', np.ones_like(vis)), ('visited', vis)):
                        sel = (kind == setName) & mask
                        if setName == 'path' and region == 'all':
                            continue
                        for tag in ('', 'Ref'):
                            pool.setdefault((setName, region, tag), []).append(res[tag][sel])
                        rows.append({'model': 'US', 'preferences': 'CRRA', 'ρ': ρ, 'ns': ns, 'nTau': len(τGrid),
                                     'refine': refine, 'date': date, 'set': setName, 'region': region}
                                    | {f'{k}': v for k, v in stats(res[''][sel]).items()}
                                    | {f'{k}Ref': v for k, v in stats(res['Ref'][sel]).items()} | _tag())
            for (setName, region, tag), v in pool.items():
                if tag:
                    continue
                st, stF = stats(np.concatenate(v)), stats(np.concatenate(pool[(setName, region, 'Ref')]))
                rows.append({'model': 'US', 'preferences': 'CRRA', 'ρ': ρ, 'ns': ns, 'nTau': len(τGrid),
                             'refine': refine, 'date': 'all', 'set': setName, 'region': region}
                            | st | {f'{k}Ref': w for k, w in stF.items()} | _tag())
                print(f'  ρ={ρ} ns={ns} {setName:<4} {region:<7}: max {100*st["max"]:.2e} p.p., p99 '
                      f'{100*st["p99"]:.2e}; reference FOC max {100*stF["max"]:.2e}, p99 {100*stF["p99"]:.2e}')
            print(f'  ρ={ρ} ns={ns}: {time.time()-tic:.0f}s')
            mergeWrite(out, rows, ['model', 'ρ', 'ns', 'date', 'set', 'region'])
    return out


def stageJumps(ρs, nss, nScan, width, out):
    """ tau_2020 along omega in [omega*(1-width), omega*(1+width)], nScan points, beta and X fixed. """
    rows = []
    for ρ in ρs:
        row = sweepRow(ρ)
        for ns in nss:
            tic = time.time()
            for k, ω in enumerate(float(row['ω'])*(1 + np.linspace(-width, width, nScan))):
                m = build(ρ, row['β'], ω, row['X'], BASE_TSS, GS | {'ns': ns})
                outp = m.solvePEE_CRRA()
                rows.append({'ρ': ρ, 'ns': ns, 'k': k, 'ω': ω, 'ωCal': float(row['ω']), 'τ2020': float(outp['τ'].values[2]),
                             'τ2050': float(outp['τ'].values[3])} | multiplicityColumns(outp['multiplicity']) | _tag())
            d = pd.DataFrame([r for r in rows if r['ρ'] == ρ and r['ns'] == ns])
            d2 = np.diff(d['τ2020'].values, 2)
            print(f'  ρ={ρ} ns={ns}: max |second difference of τ_2020| {np.abs(d2).max():.2e}, median '
                  f'{np.median(np.abs(d2)):.2e}  ({time.time()-tic:.0f}s)')
            mergeWrite(out, rows, ['ρ', 'ns', 'k'])            # an integer key: a float read back need not print alike
    return out


def stageScheme(ρs, refine, out):
    """ The whole CRRA solve at the calibrated parameters, published tax grid against the reference
    (referenceDerivatives on every period): the baseline and the theta = 0 / theta = 1 paths read at 2020/2050.
    Parameters are held at the calibration, so the difference is the tax grid's alone. """
    rows = []
    for ρ in ρs:
        row = sweepRow(ρ)
        for scheme in ('published', 'reference'):
            tic = time.time()
            for scen in ('baseline', 'theta0', 'theta1'):
                m = build(ρ, row['β'], row['ω'], row['X'], BASE_TSS, GS)
                if scen != 'baseline':
                    m, _ = sh.shockedCopy(m, scen, {})
                if scheme == 'published':
                    outp = m.solvePEE_CRRA()
                else:
                    with referenceDerivatives(m, refine):
                        outp = m.solvePEE_CRRA()
                ids = {'ρ': ρ, 'scheme': scheme, 'nTau': 101 if scheme == 'published' else 100*refine + 1,
                       'scenario': scen}
                for r in readRows(m, outp['τ'].values, outp['report'], m.db['θ'].values, ids,
                                  multiplicityColumns(outp['multiplicity'])):
                    if r['date'] in (2020, 2050):
                        rows.append(r)
                print(f'  ρ={ρ} {scheme:<9} {scen:<8}: τ_2020={outp["τ"].values[2]:.6f} τ_2050={outp["τ"].values[3]:.6f}  '
                      f'({time.time()-tic:.0f}s)')
            mergeWrite(out, rows, ['ρ', 'scheme', 'scenario', 'date'])
    return out


def main():
    p = argparse.ArgumentParser(description = 'Horizon, terminal condition and FOC residuals (US).')
    p.add_argument('--stage', nargs = '+', default = ['log', 'horizon', 'residuals'],
                   choices = ('log', 'horizon', 'design', 'residuals', 'jumps', 'scheme'))
    p.add_argument('--rho', type = float, nargs = '+', default = [0.5, 2.0])
    p.add_argument('--kSS', type = int, default = 40, help = 'periods appended for the steady-state terminal condition')
    p.add_argument('--terminal', nargs = '*', default = None, choices = ('finite', 'steadyState'),
                   help = 'design stage: run only these terminal conditions')
    p.add_argument('--ns', type = int, nargs = '+', default = [150, 50], help = 'residuals/jumps: state grids')
    p.add_argument('--refine', type = int, default = 8, help = 'residuals, scheme: the reference tax-grid refinement')
    p.add_argument('--nScan', type = int, default = 41)
    p.add_argument('--width', type = float, default = 1e-3)
    p.add_argument('--theta', default = None,
                   help = "the pipeline's main θ variant, config.US['theta'] (with --main); 'index' the workbook's θ")
    p.add_argument('--main', action = 'store_true', help = "with --theta config.US['theta']")
    a = p.parse_args()
    if a.theta in (None, oecd26.INDEX) and not a.main:
        a.theta = oecd26.indexTheta(a.theta, [OUTDIR], 'runHorizon.py')
    setTheta(a.theta, a.main)
    print(f"θ variant {THETA['variant']!r}: θ0 = {THETA['θ0']!r}")
    for st in a.stage:
        tic = time.time()
        print(f'\n=== stage {st} ===')
        if st == 'horizon':
            stageHorizon(a.rho, a.kSS, os.path.join(OUTDIR, 'US_horizon.csv'))
        elif st == 'log':
            stageLog(a.kSS, os.path.join(OUTDIR, 'US_horizon.csv'))
        elif st == 'design':
            stageDesign(a.rho, a.kSS, os.path.join(OUTDIR, 'US_horizonDesign.csv'), only = a.terminal)
        elif st == 'residuals':
            stageResiduals(a.rho, a.ns, a.refine, os.path.join(OUTDIR, 'US_residuals.csv'))
        elif st == 'jumps':
            stageJumps(a.rho, a.ns, a.nScan, a.width, os.path.join(OUTDIR, 'US_jumps.csv'))
        elif st == 'scheme':
            stageScheme(a.rho, a.refine, os.path.join(OUTDIR, 'US_derivativeScheme.csv'))
        print(f'=== stage {st} done ({time.time()-tic:.0f}s) ===')
    print('DONE')
    return 0


if __name__ == '__main__':
    sys.exit(main())

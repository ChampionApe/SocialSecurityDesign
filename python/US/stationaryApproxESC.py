r""" Stationary vs date-specific policy functions for the ENDOGENOUS design (leaded choice under the wedge).

Run:  .venv\Scripts\python.exe python\US\stationaryApproxESC.py --commonX            # rho = 1 (LOG), 0.5, 2 (CRRA, exact 2-D)
      ... --rho 1 --spec scale --phi 0.5 --ns 150 --nCand2D 41

The endogenous-theta counterpart of stationaryApprox.py. The object is now a pair of policy functions per
date, the tax tau_t(.) and the design chosen for t+1, Theta_t(.), over the inherited design (LOG) or over
(s_{t-1}, theta_t) (CRRA, LeadedCRRA2D -- the published method). The design choice is forward-looking by
construction: the electorate at t weighs the continuation at t+1 and t+2, which carry nu_{t+1}, nu_{t+2}. A
stationary policy function computed at nu_t holds those at nu_t. Under LOG the tax function itself is
static in (tau_t, theta_t), so the whole approximation error sits in the design; under CRRA both move.

exact  the paper's free path (pinAtT0 = False, the choice binding from the first period), on the model
       calibrated under the wedge at the calibrated cost p read from results/esc/escCalibration{,CRRA}.csv.
stat   for each distinct nu_t, the same recursion on a deep copy with nu frozen at nu_t, read at its first
       period (on the baseline's savings grid under CRRA). Convergence: the gap between the first two
       periods' design and tax policies over the visited states. The approximate path walks the
       stationary functions from the exact initial state (theta*, and s_0 under CRRA).
Also: the stationary function evaluated at the exact path's own states (the error in the function alone).

Writes results/numerical/US_ESC_stationaryApprox.csv, one row per (rho, t), each with the exact solve's counts
(policy.multiplicityColumns): the tax rule's (nEqMax, nCandMax, nFallback) and the root design layer's of the
2-D recursion (nEqθMax, nBrθMax, nFallbackθ; -1 under LOG, whose design is a grid argmax).

Under CRRA the cost is the converged exact calibration of escCalibrationCRRA.csv at the run's θ; with only a
path-iteration row at that (rho, spec, phi) the run stops, unless --allowPathP lets the path p stand in (a smoke run).
A relative --out is the caller's working directory's. --theta C --main is the pipeline's; --theta index (or an --out
outside results/) the workbook's θ, which without --theta is refused into results/ (oecd26.indexTheta).
"""
import os, sys, argparse, time
from copy import deepcopy
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                 # where a relative --out points: the chdir below moves to HERE
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import runESC, runESCcrra
from runESC import pickCalib
from runESCcrra import readCalibratedP, GSC
from policy import multiplicityColumns

ESCDIR = os.path.join(REPO, 'results', 'esc')
OUT = os.path.join(REPO, 'results', 'numerical', 'US_ESC_stationaryApprox.csv')


def frozen(m, ν):
    """ A deep copy of m with nu frozen at `nu` over the whole horizon and the warm-start caches cleared. """
    ms = deepcopy(m)
    ms.x0, ms.LOG.x0, ms.CRRA.x0 = {}, {}, {}
    ms.db.update(ms.adjPar('ν', np.full(ms.T, float(ν))))
    ms.updateAuxPars()
    return ms


def rows(ρ, method, m, ν, θE, θA, θAatE, τE, τA, τAatE, sE_, sA_, conv, mult = None):
    dates = m.db['dates']
    out = []
    for pos, t in enumerate(m.db['t']):
        νt = np.round(ν[pos], 12)
        out.append({'ρ': ρ, 'method': method, 't': int(t), 'date': int(dates[pos]) if pos < len(dates) else np.nan,
                    'ν': ν[pos], 'θ_exact': θE[pos], 'θ_stat': θA[pos], 'θ_statAtExactState': θAatE[pos],
                    'gapθ': θA[pos] - θE[pos], 'gapθAtExactState': θAatE[pos] - θE[pos],
                    'τ_exact': τE[pos], 'τ_stat': τA[pos], 'τ_statAtExactState': τAatE[pos],
                    'gapτ_pp': 100*(τA[pos] - τE[pos]), 'gapτAtExactState_pp': 100*(τAatE[pos] - τE[pos]),
                    's_exact': sE_[pos] if sE_ is not None else np.nan,
                    's_stat': sA_[pos] if sA_ is not None else np.nan,
                    'statConvθ': conv[νt][0], 'statConvτ': conv[νt][1]} | multiplicityColumns(mult, design = True))
    return out


def runLOG(spec, phi, commonX):
    path = os.path.join(ESCDIR, 'escCalibration.csv')
    calib = runESC.checkCalib(pd.read_csv(path), path, [spec], [phi], commonX)      # under a θ variant, its rows
    hit = pickCalib(calib, spec, phi, commonX)
    if hit.empty:
        raise KeyError(f'no converged LOG wedge calibration for ({spec}, {phi}, commonX={commonX})')
    p = float(hit.iloc[0]['p'])
    m = runESC.buildUS({'spec': spec, 'phi': phi, 'p': p}, ρ = 1.0, commonX = commonX)
    m.calibrate()
    tIdx, T = m.db['t'], m.T
    ν = m.db['ν'].values.astype(float)
    θStar = float(m.db['θ'].xs(m.t0Year))
    tic = time.time()
    led = m.solveLeaded(pinAtT0 = False)
    θE, τE = led['θ'].values, led['τ'].values
    print(f'rho=1 (LOG, p={p:.4f}): exact solve {time.time()-tic:.0f}s, theta path '
          + '  '.join(f'{x:.4f}' for x in θE[:6]))

    stat, conv = {}, {}
    for νt in sorted(set(np.round(ν, 12))):
        tic = time.time()
        ms = frozen(m, νt)
        sols = ms.ESC.solveBackward()
        t0, t1 = tIdx[0], tIdx[1]
        conv[νt] = (float(np.nanmax(np.abs(sols[t0]['θNext'] - sols[t1]['θNext']))),
                    float(np.nanmax(np.abs(sols[t0]['τ'] - sols[t1]['τ']))))
        stat[νt] = {'m': ms, 'sols': sols}
        print(f'   stationary at nu={νt:.4f}: {time.time()-tic:.0f}s, first-vs-second-period gap '
              f'theta {conv[νt][0]:.2e}, tau {conv[νt][1]:.2e}; design chosen {float(np.mean(sols[t0]["θNext"])):.4f}')

    def choice(νt, θt):
        s0 = stat[νt]['sols'][tIdx[0]]
        return float(np.interp(θt, s0['θGrid'], s0['θNext']))

    θA, τA = np.empty(T), np.empty(T)
    θA[0] = θStar
    θAatE = np.full(T, np.nan); θAatE[0] = θStar
    τAatE = np.empty(T)
    for pos in range(T):
        νt = np.round(ν[pos], 12)
        τA[pos] = stat[νt]['m'].ESC.τAt(tIdx[0], θA[pos])
        τAatE[pos] = stat[νt]['m'].ESC.τAt(tIdx[0], θE[pos])
        if pos < T - 1:
            θA[pos+1] = choice(νt, θA[pos])
            θAatE[pos+1] = choice(νt, θE[pos])
    return rows(1.0, 'LOG', m, ν, θE, θA, θAatE, τE, τA, τAatE, None, None, conv, led['multiplicity'])


def runCRRA(ρ, spec, phi, commonX, ns, nCand2D, allowPathP = False):
    fCal = os.path.join(ESCDIR, 'escCalibrationCRRA.csv')
    p, from_ = readCalibratedP(fCal, ρ, spec, phi, commonX, 'exact', fallback = allowPathP)
    if not np.isfinite(p):
        raise KeyError(f'no exact CRRA wedge calibration for (rho={ρ}, {spec}, {phi}, commonX={commonX}); '
                       '--allowPathP uses the path-iteration p')
    if from_ != 'exact':
        print(f'  WARNING: no exact calibration on file at rho={ρ}; using the {from_} p = {p:.6f} (--allowPathP)')
    gs = GSC | {'ns': ns}
    m = runESCcrra.buildUS(ρ, {'spec': spec, 'phi': phi, 'p': p}, commonX = commonX, gs = gs, nθCand2D = nCand2D)
    m.calibrate()
    tIdx, T = m.db['t'], m.T
    ν = m.db['ν'].values.astype(float)
    θStar = float(m.db['θ'].xs(m.t0Year))
    tic = time.time()
    led = m.solveLeaded2D(pinAtT0 = False, verbose = True)
    θE, τE, sE_ = led['θ'].values, led['τ'].values, led['report']['s_'].values
    sGrid = led['sols'][tIdx[0]]['sGrid']
    print(f'rho={ρ} (CRRA 2-D, p={p:.4f}): exact solve {time.time()-tic:.0f}s, theta path '
          + '  '.join(f'{x:.4f}' for x in θE[:6]))
    vis = (sGrid >= 0.8*sE_.min()) & (sGrid <= 1.2*sE_.max())

    stat, conv = {}, {}
    for νt in sorted(set(np.round(ν, 12))):
        tic = time.time()
        ms = frozen(m, νt)
        sols = ms.ESCC2.solvePolicies(θStar, pinPos = None, sGrid = sGrid, verbose = False)
        t0, t1 = tIdx[0], tIdx[1]
        conv[νt] = (float(np.nanmax(np.abs(sols[t0]['θNext'] - sols[t1]['θNext'])[vis])),
                    float(np.nanmax(np.abs(sols[t0]['τ'] - sols[t1]['τ'])[vis])))
        stat[νt] = {'m': ms, 'sols': sols}
        print(f'   stationary at nu={νt:.4f}: {time.time()-tic:.0f}s, first-vs-second-period gap '
              f'theta {conv[νt][0]:.2e}, tau {conv[νt][1]:.2e}')

    θA, τA, sA_ = np.empty(T), np.empty(T), np.empty(T)
    θAatE, τAatE = np.full(T, np.nan), np.empty(T)
    θA[0] = θAatE[0] = θStar
    s_ = float(led['s0'])
    for pos in range(T):
        sol0 = stat[np.round(ν[pos], 12)]['sols'][tIdx[0]]
        sA_[pos] = s_
        τA[pos] = float(sol0['τPolicy'](s_, θA[pos]))
        τAatE[pos] = float(sol0['τPolicy'](sE_[pos], θE[pos]))
        if pos < T - 1:
            θA[pos+1] = float(sol0['θPolicy'](s_, θA[pos]))
            θAatE[pos+1] = float(sol0['θPolicy'](sE_[pos], θE[pos]))
            s_ = float(sol0['sPolicy'](s_, θA[pos]))
    return rows(ρ, 'CRRA2D', m, ν, θE, θA, θAatE, τE, τA, τAatE, sE_, sA_, conv, led['multiplicity'])


def main():
    p = argparse.ArgumentParser(description = 'Stationary vs date-specific policies for the endogenous design.')
    p.add_argument('--rho', type = float, nargs = '*', default = [1.0, 0.5, 2.0])
    p.add_argument('--spec', default = 'scale')
    p.add_argument('--phi', type = float, default = 0.5)
    p.add_argument('--commonX', action = 'store_true')
    p.add_argument('--ns', type = int, default = GSC['ns'])
    p.add_argument('--nCand2D', type = int, default = 41)
    p.add_argument('--out', default = OUT)
    p.add_argument('--allowPathP', action = 'store_true',
                   help = 'under CRRA, use the path-iteration p where no exact calibration is on file (a smoke run)')
    p.add_argument('--theta', default = None,
                   help = "the pipeline's main θ variant (with --main): runESC.setTheta's θ0, rows carry it")
    p.add_argument('--main', action = 'store_true', help = "with --theta config.US['theta']")
    a = p.parse_args()
    a.out = os.path.abspath(os.path.join(CWD, a.out))
    import oecd26
    if a.theta in (None, oecd26.INDEX) and not a.main:
        a.theta = oecd26.indexTheta(a.theta, [a.out], 'stationaryApproxESC.py')
    if a.theta is not None or a.main:
        if not (a.theta is not None and a.main):
            raise SystemExit('--theta needs --main here (and --main --theta)')
        runESC.setTheta(a.theta, None, ('calib',), False, a.commonX, a.rho[0], main = True)
        try:
            for ρ in a.rho:
                oecd26.usBeta(runESC.usSweep(a.commonX), ρ, thetaSource = a.theta)
        except (FileNotFoundError, ValueError, KeyError) as e:
            raise SystemExit('--theta {}: {}'.format(a.theta, e))
    allRows = []
    for ρ in a.rho:
        allRows += [r | runESC.thetaColumns() for r in (runLOG(a.spec, a.phi, a.commonX) if np.isclose(ρ, 1)
                    else runCRRA(ρ, a.spec, a.phi, a.commonX, a.ns, a.nCand2D, a.allowPathP))]
        df = pd.DataFrame(allRows)
        os.makedirs(os.path.dirname(a.out), exist_ok = True)
        df.to_csv(a.out, index = False)                     # after every rho: a long run is readable as it goes
    pd.set_option('display.width', 220)
    print('\n' + df[['ρ', 'method', 'date', 'ν', 'θ_exact', 'θ_stat', 'θ_statAtExactState', 'gapθ', 'gapθAtExactState',
                     'τ_exact', 'τ_stat', 'gapτ_pp', 'gapτAtExactState_pp', 'statConvθ', 'statConvτ']].round(5).to_string(index = False))
    print(f'\nwrote {a.out}')


if __name__ == '__main__':
    main()

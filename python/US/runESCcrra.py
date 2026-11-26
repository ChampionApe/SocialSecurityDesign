r""" Endogenous system characteristics under CRRA: calibrate the wedge at rho != 1, solve the design path,
and run the ageing counterfactual.

Run:  .venv\Scripts\python.exe python\US\runESCcrra.py                       # rho = 2, spec 'scale', phi = .5
      ... --rho 2.0 0.5  --spec scale  --phi 0.5  --commonX
      ... --stage calib path shocks
      ... --exact --stage calib path shocks --ns 150                        # the published method
      ... --exact --rho 2.0 --spec size --phi 0.5 --commonX --bracket L U     # the paper's wedge, lambda in p

TWO METHODS, one `method` column. Every row written here carries method = 'exact' or 'path', and the
column is part of every merge key, so the two vintages coexist in one csv and the paper pipeline selects
by it (python/paper/datasets.py, config.US['esc']['exact']).

  --exact   LeadedCRRA2D through ModelESC.solveLeaded2D: backward-solved policy functions over the
            two-dimensional state (s_{t-1}, theta_t), no held-fixed future. The PUBLISHED method for the
            CRRA ESC leg. One full recursion per calibration trial and per chosen-design reading
            (~6 min at ns = 150, ~95 s at ns = 50; the two grids agree to 3e-5 in the design, the
            theta-state grid is what the accuracy hangs on -- test_escCRRA.py). The wedge calibration
            therefore SCANS on the coarse grid (--nsScan) and refines the bracketed root on --ns, so it
            needs no warm start from the path iteration; a path p on file only narrows the scanned
            bracket to half..double of it.
            --designRule picks its design layer: 'root' (the default, alg esc:crra2D of num_esc.tex, --Ma
            tabulation nodes) or 'legacy' (the earlier layer, for comparisons).
  default   LeadedCRRA's path iteration, the development stand-in certified against the exact recursion
            to ~0.01 in the design (test_escCRRA.py). No paper output reads it.

Every row also carries designRule, Ma and the counts of the solve behind it (solverColumns), as non-key
columns: nEqMax, nCandMax, nFallback (the tax rule) and nEqθMax, nBrθMax, nFallbackθ (the root design
layer), -1 where that solve counted nothing.

Why this is a separate driver from runESC.py: under CRRA nothing about the leaded choice is cheap. The tau
FOC depends on s_{t-1}, so there is no static tauPolicy(theta); W_t is not additively separable, so the
choice depends on both states; and every candidate design costs a full solvePEE_CRRA (~3 s) instead of a
closed form. The path iteration budgets:

    calibrate p     ~70 s per trial value  (a 35 s recalibration of (beta, omega) + 13 candidate solves)
    the design path ~2 min per iteration
    each shock      one path solve per reading

STAGES
  calib   p such that the equilibrium design IN FORCE at t0 reproduces theta*, per (rho, spec, phi) --
          ModelESC.leadedDesignAtT0_CRRA (path) or leadedDesignAtT0_2D (exact), the design chosen one
          period BEFORE 2020, which is what the tables read. The bracket is scanned an order of magnitude
          BELOW the LOG one: at rho = 2 the LOG-calibrated p = 0.41 already puts the choice at the
          theta = 1 corner (a higher EIS needs far less wedge to reach an interior choice).
  path    the design path at the calibrated p, plus targetDrift.
  shocks  counterfactuals, theta pinned vs chosen, each a NEW EQUILIBRIUM PATH (shocked parameters over
          the whole horizon, own steady state, the choice binding from the first period) and REPORTED AT
          t0. Default scenario set: acute ageing plus the French characteristics -- income distribution,
          leisure, voting, their income+voting combination and all three at once -- the same set
          runESC.py's LOG stage runs, selectable via --scenarios, plus France's own calibrated path as
          the endpoint row. The French data are frenchData's (runShocksUS.py), calibrated at the SAME rho
          under CRRA, computed once per rho and shared across specs.
  sequential  the COSTLESS sequential FOC of eq:esc:seqFOC (ModelESC.sequentialFOC) on the solved
              no-wedge baseline at every dated period, escSequentialCRRA.csv: its sign over theta in
              [0,1] is the "three timings, one corner" claim under CRRA (the LOG case is test_esc.py's).

--theta VARIANT [--root DIR | --main] (stages calib, path, shocks; runESC.setTheta): the U.S. at θ0 = the
variant's U.S. row of data/oecd26_theta.csv, France and the UK at theirs (runESC.countryθ0), seeded off the variant's
sweeps in DIR/calibration/ (every rho asked for at the variant: oecd26.usBeta; the UK's for --host UK), the csvs in
DIR/esc/ (results/ under --main, the pipeline's main variant), every row carrying thetaSource, θ0 and grid
(solverColumns); a path or shocks stage reads only a converged calibration row at the same three (readCalibratedP).
--theta index: the workbook's θ into results/esc/; without --theta the run is refused (oecd26.indexTheta).

--spendingTarget gross|net (runESC.SPENDING): the calibrations' spending target under a cost, a key column of the three
csvs (blank = 'gross'); readCalibratedP reads only the run's. targetDrift's 'τ' (τDrift) is the drift of the targeted
quantity, tau f under 'net'; the calib row's τ0/fStar are at the calibrated tax (ModelESC.calibratedTax0); path rows
carry τNet, shocks rows τNet_t0, τNet_t1 (benefits paid over labour income, ModelESC.netTax).

--outdir DIR: the csvs' directory instead of the one above (runESC.py's --outdir); every csv this driver reads is one
it writes there (escCalibrationCRRA{host}{tag}.csv: a path stage's cost, the exact scan's narrowed bracket, the untagged
fallback), the sweeps stay --theta's. The cost of results/esc/ is runESC.PAPERCOST: a run under another SSD_TAUPOWER or
SSD_CARRIER is refused there (runESC.refuseCost), and its rows carry tauPower and carrier (solverColumns). --bracket
is in the run's λ units, which scale by about τ0^(1-k) (BRACKET_CRRA's default carries TAU_SCALE).

A stage that reads its cost from file (no calib stage in the run) under --exact needs an exact row: with only a
path-iteration row at that (rho, spec, phi) the combination is a recorded failure, unless --allowPathP lets the path p
stand in (a smoke run). Failures a run survives (a calibration, a path or a scenario that raised) are printed, the run
goes on, and main exits 1 with the list (runESC.noteFailure). A relative --root is the caller's working directory's.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd
from copy import deepcopy

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                 # where a relative --root points: the chdir below moves to HERE
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import test as testmod
import shocks as sh
import runESC
from modelESC import ModelESC
from runESC import SHOCKS_ESC, mergeWrite, buildEU, readout as escReadout, wedgeReadout, piColumn
from runShocksUS import frenchData, hostRow
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'esc')
GSC = {'n': 101, 'ns': 150, 'smoothKnots': 4, 'interpKind': 'linear'}

# The default scan bracket per spec when --bracket is not given (see the argument's comment): p for the
# theta^p forms, lambda for 'size' (about 8.6 at rho = 1 under LOG and expected to fall in rho, as p does).
from modelESC import TAU_SCALE      # the lambda brackets scale with tau0^(1-k)
BRACKET_CRRA = {'scale': (0.01, 3.0), 'flat': (0.01, 3.0), 'size': (0.25*TAU_SCALE, 200.*TAU_SCALE),
                'sizePiR': (0.5*TAU_SCALE, 400.*TAU_SCALE)}

# Row-identity keys, rho-indexed counterparts of runESC's. commonX is in all three -- see mergeWrite.
# method too: an exact row and a path row answer the same question by different solvers and must not
# overwrite each other (the pipeline selects by method).
# spendingTarget too, blank at 'gross' (runESC.spendingColumns, runESC.mergeWrite).
KEYCAL = ['ρ', 'spec', 'phi', 'commonX', 'method', 'spendingTarget']
KEYPATH = ['ρ', 'spec', 'phi', 'commonX', 'method', 'spendingTarget', 'pos']
KEYSHK = ['ρ', 'spec', 'phi', 'commonX', 'method', 'spendingTarget', 'scenario', 'θpinned']


def buildUS(ρ, wedge = None, nθCandCRRA = 13, commonX = False, gs = None, nθCand2D = 21, designRule = 'root',
            Ma = 5):
    """ runESC.buildUS's CRRA counterpart: (beta, omega) seeded from the matching variant's sweep, the
    calibration itself redone by the caller. gs: the CRRA grid settings (GSC by default); LeadedCRRA2D
    borrows them at solve time, so --ns reaches the exact solver through here. nθCand2D: the exact
    recursion's candidate grid for θ_{t+1} -- the objective is flat near its maximum (~1e-5 in W over
    ±0.01 in θ at the frVoting choice), so the grid sets the resolution of the reported design.
    designRule, Ma: LeadedCRRA2D's design layer and the root layer's tabulation nodes. The sweep and θ0 are the
    run's θ variant's (runESC.usSweep, runESC.THETA; results/calibration/ and the workbook's θ without one). """
    gs = GSC if gs is None else gs
    row = pd.read_csv(runESC.usSweep(commonX))
    row = row.loc[(row['ρ'] - ρ).abs() < 1e-9].iloc[-1]
    θ0 = runESC.THETA['θ0']
    m = ModelESC(pars = testmod.pars | {'ρ': float(ρ), 'β': float(row['β']), 'ω': float(row['ω'])}
                 | ({} if θ0 is None else {'θ0': float(θ0)}),
                 wedge = wedge, nθCandCRRA = nθCandCRRA, nθCand2D = nθCand2D, commonX = commonX,
                 designRule = designRule, Ma = Ma, spendingTarget = runESC.SPENDING['target'], **testmod.kwargs)
    m.db['dates'], m.db['workweek'] = testmod.dates, testmod.workweek
    m.CRRA.initGS(gs)
    m.LOG.initGS({k: v for k, v in gs.items() if k != 'ns'})
    return m


def buildHost(host, ρ, wedge = None, nθCandCRRA = 13, commonX = False, gs = None, nθCand2D = 21,
              designRule = 'root', Ma = 5):
    """ buildUS for host 'US'; for 'UK' the UK's ModelESCFR (beta imposed from the US, omega calibrated),
    omega seeded from the UK's no-wedge sweep at this rho, with the same grid and design-layer settings. """
    if host == 'US':
        return buildUS(ρ, wedge, nθCandCRRA = nθCandCRRA, commonX = commonX, gs = gs, nθCand2D = nθCand2D,
                       designRule = designRule, Ma = Ma)
    gs = GSC if gs is None else gs
    m = buildEU(host, wedge, ρ = ρ, commonX = commonX, ω = float(hostRow(ρ, commonX, host, runESC.THETA['calibDir'])['ω']),
                nθCandCRRA = nθCandCRRA, nθCand2D = nθCand2D, designRule = designRule, Ma = Ma)
    m.CRRA.initGS(gs)
    m.LOG.initGS({k: v for k, v in gs.items() if k != 'ns'})
    return m


def franceRowCRRA(m, ρ, spec, phi, p, hbarRef, commonX = False, gs = None):
    """ runESC.franceRow at rho != 1: France's own calibrated equilibrium at 2020 under the same wedge,
    exogenous theta, in escShocksCRRA's row schema. See runESC.franceRow for what the row means and why
    its workweek is a target rather than a prediction. Exogenous theta, so the row is method-free; the
    caller stamps the method it is filed under. """
    mFR = buildEU('FR', {'spec': spec, 'phi': phi, 'p': p}, ρ = ρ, commonX = commonX)
    mFR.CRRA.initGS(GSC if gs is None else gs)
    mFR.calibrate(preferences = 'CRRA')
    pos = mFR.db['t0']
    out = mFR.solvePEE_CRRA()
    ww = float(m.db['workweek'])
    r0 = escReadout(mFR, out['τ'], out['report'], hbarRef, pos, workweekData = ww)
    r1 = escReadout(mFR, out['τ'], out['report'], hbarRef, pos+1, workweekData = ww)
    θ = mFR.db['θ'].values.astype(float)
    return {'ρ': ρ, 'spec': spec, 'phi': phi, 'commonX': commonX, 'p': p, 'scenario': 'France',
            'θpinned': True,
            'θ_tm1': float(θ[pos-1]), 'θ_t0': float(θ[pos]), 'θ_t1': float(θ[pos+1]),
            'τ_t0': r0['τ'], 'sr_t0': r0['sr'], 'ww_t0': r0['workweek'],
            'τ_t1': r1['τ'], 'sr_t1': r1['sr'], 'ww_t1': r1['workweek']} | netColumns(mFR, out['τ'], θ, pos) | {
            'multiplicity': out.get('multiplicity')}


def netColumns(m, τ, θ, pos):
    """ {'τNet_t0', 'τNet_t1'}: benefits paid over labour income at positions pos and pos+1 of a solved path
    (ModelESC.netTax at db's pi, the solve's own; tau itself without a cost). τ, θ: the path's (T,). """
    τ, θ = np.asarray(τ, dtype = float), np.asarray(θ, dtype = float)
    return {f'τNet_t{k}': float(m.netTax(float(τ[pos+k]), float(θ[pos+k]), m.db['t'][pos+k])) for k in (0, 1)}


def stageSequential(ρs, out, commonX = False, gs = None, nθ = 201):
    """ The costless sequential FOC on the solved no-wedge baseline at every non-terminal dated period
    (position 0 has no predetermined ratio). One row per (rho, pos): the largest and smallest value over
    the theta grid and whether it is negative throughout -- the corner claim is `negative` at every row --
    and the baseline solve's tax counts (policy.multiplicityColumns). """
    rows, θg = [], np.linspace(0., 1., nθ)
    for ρ in ρs:
        m = buildUS(ρ, None, commonX = commonX, gs = gs)
        m.calibrate()
        sol = m.solvePEE_LOG() if ρ == 1.0 else m.solvePEE_CRRA()
        dates = m.db['dates']
        for pos in range(1, len(m.db['t']) - 1):
            foc = m.sequentialFOC(θg, sol, pos = pos)
            k = int(np.argmax(foc))
            rows.append({'ρ': ρ, 'commonX': commonX, 'pos': pos,
                         'date': dates[pos] if pos < len(dates) else np.nan,
                         'τ': float(np.asarray(sol['τ'])[pos]),
                         'focMax': float(foc[k]), 'θAtMax': float(θg[k]), 'focMin': float(np.min(foc)),
                         'negative': bool(np.all(foc < 0.))} | multiplicityColumns(sol.get('multiplicity')))
            print('  rho={:<4} pos={} ({}): max FOC {:+.4e} at θ={:.3f}, min {:+.4e}  -> {}'.format(
                ρ, pos, rows[-1]['date'], foc[k], θg[k], np.min(foc),
                'negative throughout' if rows[-1]['negative'] else 'NOT NEGATIVE'))
        mergeWrite(out, rows, ['ρ', 'commonX', 'pos'])
    return pd.DataFrame(rows)


def readCalibratedP(fCal, ρ, spec, phi, commonX, method, fallback = False):
    """ (p, method of the row) for this (rho, spec, phi, variant) from the csv: the last converged row at `method`.
    With `fallback` and no row at 'exact', the last converged row of another method (the path iteration's p, which
    agrees with the exact one to ~0.01 in the design: a stand-in for a smoke run, not for the paper). (NaN, None)
    when nothing matches. """
    if not os.path.exists(fCal):
        return np.nan, None
    hit = runESC.atSpending(pd.read_csv(fCal))          # the run's spending target only
    hit = hit[(hit['ρ'] == ρ) & (hit['spec'] == spec) & (hit['phi'] == phi)]
    if 'converged' in hit.columns:
        hit = hit[hit['converged'].fillna(False).astype(bool)]
    if 'commonX' in hit.columns:
        hit = hit[hit['commonX'].astype(bool) == bool(commonX)]
    want = runESC.thetaColumns()          # under a θ variant only a row at its thetaSource, θ0 and grid
    if want:
        hit = (hit[(hit['thetaSource'] == want['thetaSource']) & (hit['grid'] == want['grid'])
                   & (pd.to_numeric(hit['θ0']) - want['θ0']).abs().le(1e-12)]
               if all(c in hit.columns for c in want) else hit.iloc[:0])
    if 'method' not in hit.columns:
        hit = hit.assign(method = 'path')
    hit['method'] = hit['method'].fillna('path')
    own = hit[hit['method'] == method]
    if not own.empty:
        return float(own.iloc[-1]['p']), method
    if fallback and method == 'exact' and not hit.empty:
        return float(hit.iloc[-1]['p']), str(hit.iloc[-1]['method'])
    return np.nan, None


def solverColumns(designRule, Ma, mult):
    """ The non-key columns every row of this driver carries: the run's design-layer settings (designRule,
    Ma; LeadedCRRA2D under --exact) and the counts of the solve behind the row (policy.multiplicityColumns:
    nEqMax, nCandMax, nFallback for the tax rule, nEqθMax, nBrθMax, nFallbackθ for the root design layer,
    -1 where that solve counted nothing; nRootsSMax the state fixed point's root count), and runESC.tagColumns (the
    θ variant's, the spending target and the cost where not the defaults). Never part of a merge key (finding #13)
    but for spendingTarget (KEYCAL). """
    return {'designRule': designRule, 'Ma': int(Ma)} | multiplicityColumns(mult, design = True) | runESC.tagColumns()


def multiplicityOf(led):
    """ The multiplicity summary of a solved design path: solveLeaded2D's own, or the path iteration's
    last equilibrium solve (tax counts only); None if neither carries one. """
    if 'multiplicity' in led:
        return led['multiplicity']
    return led['out'].get('multiplicity') if 'out' in led else None


def pathRowsFrom(m, led, base, ρ, spec, phi, pCal, commonX, method, designRule = 'root', Ma = 5):
    """ escPathCRRA rows from a solved design path (solveLeadedCRRA or solveLeaded2D return), each with
    solverColumns. """
    solver = solverColumns(designRule, Ma, multiplicityOf(led))
    t0 = m.t0Year
    hbarRef = float(m.B.avgHours(base['report']['h'].xs(t0), t0))
    dates = m.db['dates']
    rows = []
    for pos, t in enumerate(m.db['t'][:-1]):
        r = sh.readout(m, led['out']['τ'] if 'out' in led else led['τ'],
                       led['out']['report'] if 'out' in led else led['report'],
                       float(m.db['workweek']), hbarRef, pos = pos)
        rb = sh.readout(m, base['τ'], base['report'], float(m.db['workweek']), hbarRef, pos = pos)
        τLed = led['out']['τ'] if 'out' in led else led['τ']
        rows.append({'ρ': ρ, 'spec': spec, 'phi': phi, 'commonX': commonX, 'method': method,
                     'p': pCal, 'pos': pos,
                     'date': dates[pos] if pos < len(dates) else np.nan,
                     'ν': float(m.db['ν'].xs(t)), 'θ': float(led['θ'].xs(t)),
                     'τ': r['τ'], 'sr': r['sr'], 'workweek': r['workweek'],
                     'τNet': float(m.netTax(float(τLed.xs(t)), float(led['θ'].xs(t)), t)),
                     'τExo': rb['τ'], 'srExo': rb['sr'],
                     'converged': led.get('converged', True), 'step': led.get('step', np.nan),
                     'τDrift': led['targetDrift']['τ'], 'RDrift': led['targetDrift']['R']}
                    | piColumn(m, t) | solver)
    return rows


def main():
    p = argparse.ArgumentParser(description = 'Endogenous theta under CRRA.')
    p.add_argument('--rho', type = float, nargs = '*', default = [2.0])
    p.add_argument('--spec', nargs = '*', default = ['scale'], choices = ('scale', 'flat', 'size', 'sizePiR'),
                   help = "'size' is the paper's wedge (f(theta, tau), lambda in the p slot); 'sizePiR' the "
                          "paper's with lambda*Pi_t/R_t (ModelESC; --exact only)")
    p.add_argument('--commonX', action = 'store_true',
                   help = 'the common-X calibration variant (the paper leads with it). Written as a '
                          'column and part of the merge key -- see runESC.mergeWrite.')
    p.add_argument('--phi', type = float, nargs = '*', default = [0.5])
    p.add_argument('--stage', nargs = '*', default = ['calib', 'path', 'shocks'],
                   choices = ('calib', 'path', 'shocks', 'sequential'),
                   help = "add 'sequential' for the costless sequential FOC")
    p.add_argument('--scenarios', nargs = '*',
                   default = ['baseline', 'acute', 'frIncome', 'frLeisure', 'frVoting', 'frBoth', 'frAll'],
                   help = 'shock stage scenarios; any key of runESC.SHOCKS_ESC plus "baseline"')
    p.add_argument('--exact', action = 'store_true',
                   help = 'solve the leaded choice with the exact 2-D recursion (LeadedCRRA2D) instead of '
                          'the path iteration; rows are written with method = "exact"')
    p.add_argument('--allowPathP', action = 'store_true',
                   help = 'under --exact, a stage without calib that finds no exact calibration at a (rho, spec, phi) '
                          'uses the path-iteration p (a smoke run); otherwise that combination fails')
    p.add_argument('--ns', type = int, default = GSC['ns'],
                   help = 'the CRRA savings-state grid (default {}); the 2-D recursion borrows it, and '
                          '50 is ~4x faster at no measurable cost in the choice'.format(GSC['ns']))
    p.add_argument('--nsScan', type = int, default = 50,
                   help = 'under --exact, the grid the wedge calibration SCANS on before refining at '
                          '--ns (the warm start; 0 = scan at --ns too)')
    p.add_argument('--nCand2D', type = int, default = 21,
                   help = 'under --exact, the candidate grid for θ_{t+1} (see buildUS)')
    p.add_argument('--designRule', default = 'root', choices = ('root', 'legacy'),
                   help = "under --exact, LeadedCRRA2D's design layer: 'root' (alg esc:crra2D) or 'legacy' "
                          '(the earlier layer, for comparisons)')
    p.add_argument('--Ma', type = int, default = 5,
                   help = "under --exact and designRule 'root', the nodes on which the residual in a is "
                          'tabulated before its sign changes are closed')
    # The required cost falls steeply in the intertemporal elasticity (p is about 0.95, 0.41, 0.09 at rho = 0.5,
    # 1, 2), so the default bracket spans all three and nScan keeps the resolution per decade. Under --exact the
    # default is half..double the path iteration's p when one is on file (readCalibratedP), each node being a full
    # recursion. The default is per spec (BRACKET_CRRA): under 'size' the parameter is lambda, an order of
    # magnitude larger than p; one explicit --bracket applies to every spec in the run.
    p.add_argument('--bracket', type = float, nargs = 2, default = None)
    p.add_argument('--nScan', type = int, default = None, help = 'default 14 (path), 4 (exact)')
    p.add_argument('--xtol', type = float, default = None,
                   help = 'brentq tolerance on p; default 1e-6 (path), 2e-5 (exact)')
    p.add_argument('--nCand', type = int, default = 13)
    p.add_argument('--maxIter', type = int, default = 4)
    p.add_argument('--tag', default = '')
    p.add_argument('--host', default = 'US', choices = ('US', 'UK'),
                   help = "'UK': calibrate the UK's OWN cost (the deep test), the calib stage only; the csv "
                          'name carries the host (escCalibrationCRRAUK.csv)')
    p.add_argument('--theta', default = None,
                   help = "a variant of data/oecd26_theta.csv: the U.S., France and the UK at its θ (db['θ0'])")
    p.add_argument('--root', default = None, help = "the variant's directory (default results/oecd26_theta_<variant>)")
    p.add_argument('--main', action = 'store_true',
                   help = "with --theta config.US['theta']: read results/calibration/, write results/esc/")
    p.add_argument('--spendingTarget', default = 'gross', choices = runESC.SPENDING_TARGETS,
                   help = "runESC.py's: 'net' calibrates every cost to the benefits paid, tau_t0 f = tau0 "
                          "(targetDrift then the drift of tau f); a key column, blank at 'gross'")
    p.add_argument('--outdir', default = None,
                   help = "the csvs' directory (read and written) instead of --theta's (results/esc/ under --main); the "
                          "sweeps are still --theta's. A relative path is the working directory's.")
    a = p.parse_args()
    a.root = None if a.root is None else os.path.abspath(os.path.join(CWD, a.root))
    runESC.setSpendingTarget(a.spendingTarget)
    if a.host != 'US' and set(a.stage) != {'calib'}:
        p.error('--host {} runs the calib stage only (--stage calib)'.format(a.host))
    OUTDIR = setThetaCRRA(a)
    OUTDIR = OUTDIR if a.outdir is None else os.path.abspath(os.path.join(CWD, a.outdir))
    runESC.refuseCost(OUTDIR)
    os.makedirs(OUTDIR, exist_ok = True)
    hs = '' if a.host == 'US' else a.host
    fCal = os.path.join(OUTDIR, f'escCalibrationCRRA{hs}{a.tag}.csv')
    fPath = os.path.join(OUTDIR, f'escPathCRRA{hs}{a.tag}.csv')
    fShk = os.path.join(OUTDIR, f'escShocksCRRA{hs}{a.tag}.csv')
    method = 'exact' if a.exact else 'path'
    gs = GSC | {'ns': a.ns}
    nScan = a.nScan if a.nScan is not None else (6 if a.exact else 14)
    xtol = a.xtol if a.xtol is not None else (2e-5 if a.exact else 1e-6)
    gsScan = gs | {'ns': a.nsScan} if (a.exact and a.nsScan > 0) else gs
    print('method = {}, ns = {}{}'.format(method, a.ns, f', scan at ns = {a.nsScan}, designRule = '
                                          f'{a.designRule}, Ma = {a.Ma}' if a.exact else ''))
    layer = {'designRule': a.designRule, 'Ma': a.Ma}

    if 'sequential' in a.stage:
        print('=== the costless sequential FOC on the no-wedge baseline ===')
        stageSequential(a.rho, os.path.join(OUTDIR, f'escSequentialCRRA{a.tag}.csv'),
                        commonX = a.commonX, gs = gs)
    if set(a.stage) <= {'sequential'}:
        print()
        print('-> {}'.format(os.path.relpath(OUTDIR, REPO)))
        return 0

    calRows, pathRows, shkRows = [], [], []
    frDataCache = {}
    for ρ in a.rho:
        for spec in a.spec:
            for phi in a.phi:
                tag = f'rho={ρ}, {spec}, phi={phi}'
                # ---------------------------------------------------- calibrate p
                pCal = np.nan
                if 'calib' in a.stage:
                    tic = time.time()
                    print(f'\n=== [{tag}] calibrating p under CRRA ({method}) ===')
                    bracket, nScanHere = a.bracket, nScan
                    if bracket is None:
                        bracket = BRACKET_CRRA[spec]
                        if a.exact:
                            # a path-iteration p on file narrows the scan; without one the wide bracket
                            # is scanned on the coarse grid, which is what --nsScan is for
                            pRef, from_ = readCalibratedP(fCal, ρ, spec, phi, a.commonX, 'path')
                            if np.isfinite(pRef):
                                bracket, nScanHere = (0.5*pRef, 2.0*pRef), min(nScan, 4)
                                print('  bracket {:.4f}..{:.4f} around the {} p = {:.4f}'.format(
                                    *bracket, from_, pRef))
                            else:
                                nScanHere = max(nScan, 14)
                                print('  no path-iteration p on file: full bracket, nScan = {}'
                                      .format(nScanHere))
                    m = buildHost(a.host, ρ, {'spec': spec, 'phi': phi, 'p': 0.2}, nθCandCRRA = a.nCand,
                                commonX = a.commonX, gs = gs, nθCand2D = a.nCand2D, **layer)
                    try:
                        rec = m.calibrateWedge(spec = spec, phi = phi,
                                               preferences = 'CRRA2D' if a.exact else 'CRRA',
                                               bracket = tuple(bracket), nScan = nScanHere, xtol = xtol,
                                               beforeScan = lambda: m.CRRA.initGS(gsScan),
                                               beforeRefine = lambda: m.CRRA.initGS(gs))
                        pCal = rec['p']
                        # the wedge read at the calibrated point's tax: tau0, the calibration target, which
                        # the path stage's targetDrift measures the solved path against (under 'net' the tax
                        # whose benefits paid are tau0, ModelESC.calibratedTax0)
                        wr = (wedgeReadout(m, m.calibratedTax0()) if rec['converged']
                              else {'Vtilde': float(m.B.Vtilde(m.t0Year))})
                        calRows.append({'ρ': ρ, 'spec': spec, 'phi': phi, 'commonX': a.commonX, 'host': a.host,
                                        'method': method, 'p': rec['p'],
                                        'converged': rec['converged'], 'θStar': rec['θ'],
                                        'residual': rec['residual'], 'message': rec['message'],
                                        'β': m.simpleβinv(), 'ω': float(m.db['ω'].xs(m.t0Year)),
                                        'ns': a.ns, 'nsScan': a.nsScan if a.exact else np.nan,
                                        'nCand2D': a.nCand2D if a.exact else np.nan,
                                        'nScan': len(rec['scan']),
                                        'seconds': time.time()-tic} | wr
                                       # calibrateWedge's last recursion is the one at the calibrated p
                                       | solverColumns(a.designRule, a.Ma,
                                                       getattr(m.ESCC2, 'lastMultiplicity', None)
                                                       if a.exact and rec['converged'] else None))
                        print('  -> p={}  ({})  [{:.0f}s]'.format(rec['p'], rec['message'], time.time()-tic))
                    except Exception as e:
                        runESC.noteFailure(f'{tag} calibration', e)
                        calRows.append({'ρ': ρ, 'spec': spec, 'phi': phi, 'commonX': a.commonX, 'host': a.host,
                                        'method': method, 'p': np.nan,
                                        'converged': False, 'message': f'{type(e).__name__}: {e}'}
                                       | solverColumns(a.designRule, a.Ma, None))
                    mergeWrite(fCal, calRows, KEYCAL)
                else:
                    # a --tag run without its own calibration reads the untagged one (smoke runs)
                    for fc in (fCal, os.path.join(OUTDIR, f'escCalibrationCRRA{hs}.csv')):
                        pCal, from_ = readCalibratedP(fc, ρ, spec, phi, a.commonX, method, fallback = a.allowPathP)
                        if np.isfinite(pCal):
                            break
                    if np.isfinite(pCal) and from_ != method:
                        print(f'  [{tag}] WARNING: no {method} calibration on file; using the {from_} '
                              f'p = {pCal:.6f} as a stand-in (--allowPathP).')
                    if not np.isfinite(pCal):
                        runESC.noteFailure(tag, LookupError('no converged {} calibration on file (--stage calib, or '
                                                            '--allowPathP under --exact)'.format(method)))
                if not np.isfinite(pCal):
                    print(f'  [{tag}] no calibrated p -- skipping the remaining stages.')
                    continue

                # ---------------------------------------------------- the design path
                m = buildHost(a.host, ρ, {'spec': spec, 'phi': phi, 'p': pCal}, nθCandCRRA = a.nCand,
                            commonX = a.commonX, gs = gs, nθCand2D = a.nCand2D, **layer)
                m.calibrate()
                t0 = m.t0Year
                θStar = float(m.db['θ'].xs(t0))
                base = m.solvePEE_CRRA()
                hbarRef = float(m.B.avgHours(base['report']['h'].xs(t0), t0))
                ledBase = None          # the free path on the calibrated model: shared by path and shocks
                if 'path' in a.stage:
                    tic = time.time()
                    print(f'\n=== [{tag}] design path (p={pCal:.4f}, θ*={θStar:.4f}, {method}) ===')
                    try:
                        # pinAtT0 explicitly, not by default: the reported path is the free one, so that
                        # theta_{t0} is the equilibrium design the shocks stage and the tables read.
                        ledBase = (m.solveLeaded2D(pinAtT0 = False, verbose = True) if a.exact
                                   else m.solveLeadedCRRA(maxIter = a.maxIter, pinAtT0 = False))
                        pathRows += pathRowsFrom(m, ledBase, base, ρ, spec, phi, pCal, a.commonX, method,
                                                 **layer)
                        print('  θ path: {}   (converged={}, {:.0f}s)'.format(
                            '  '.join('{:.4f}'.format(x) for x in ledBase['θ'].values[:7]),
                            ledBase.get('converged', True), time.time()-tic))
                        mergeWrite(fPath, pathRows, KEYPATH)
                    except Exception as e:
                        runESC.noteFailure(f'{tag} design path', e)
                        ledBase = None

                # ---------------------------------------------------- counterfactuals
                if 'shocks' in a.stage:
                    pos0 = m.db['t0']
                    frData = None
                    if any(n.startswith('fr') for n in a.scenarios):
                        if ρ not in frDataCache:
                            print('  calibrating France at rho={} for the French scenarios ...'.format(ρ))
                            frDataCache[ρ] = frenchData(
                                m, ρ, 'CRRA', gs = gs, commonX = a.commonX, θ0FR = runESC.countryθ0('FR'),
                                usCsv = None if runESC.THETA['calibDir'] is None else runESC.usSweep(a.commonX))
                        frData = frDataCache[ρ]
                    for name in a.scenarios:
                        for pin in (True, False):
                            tic = time.time()
                            try:
                                # New equilibrium path, as runESC.leadedNewPath: the shock holds over the
                                # whole horizon, the economy starts at its own steady state (no s0 seed),
                                # and the design is chosen from the first period when pin is False.
                                mt, _ = sh.shockedCopy(m, name, frData, SHOCKS_ESC) \
                                        if name != 'baseline' else (deepcopy(m), None)
                                if name == 'baseline':
                                    mt.x0, mt.LOG.x0, mt.CRRA.x0 = {}, {}, {}
                                if pin:
                                    θp = mt.db['θ'].values.astype(float)
                                    o = mt.solvePEE_CRRA(θ = θp, ε = mt.db['eps'].values.astype(float))
                                    θPath, out = pd.Series(θp, index = mt.db['t']), o
                                    mult = o.get('multiplicity')
                                elif name == 'baseline' and ledBase is not None:
                                    θPath = ledBase['θ']          # the path stage already solved it
                                    out = ledBase['out'] if 'out' in ledBase else ledBase
                                    mult = multiplicityOf(ledBase)
                                elif a.exact:
                                    if mt._hasΠR():
                                        # pi of the shocked economy's pinned path: a closer start for the
                                        # 2-D recursion's pi iteration than the host's free path
                                        mt.setΠR(mt.solvePEE_CRRA()['ΠRpath'])
                                        mt._ΠRstart = {k: v for k, v in mt._ΠRstart.items()
                                                       if not k.startswith('solveLeaded2D')}
                                    rec = mt.solveLeaded2D(pinAtT0 = False)
                                    θPath, out = rec['θ'], rec
                                    mult = multiplicityOf(rec)
                                else:
                                    rec = mt.solveLeadedCRRA(maxIter = a.maxIter, verbose = False,
                                                             pinAtT0 = False)
                                    θPath, out = rec['θ'], rec['out']
                                    mult = multiplicityOf(rec)
                                r0 = sh.readout(mt, out['τ'], out['report'], float(mt.db['workweek']),
                                                hbarRef, pos = pos0)
                                r1 = sh.readout(mt, out['τ'], out['report'], float(mt.db['workweek']),
                                                hbarRef, pos = pos0+1)
                                shkRows.append({'ρ': ρ, 'spec': spec, 'phi': phi,
                                                'commonX': a.commonX, 'host': a.host, 'method': method,
                                                'p': pCal,
                                                'scenario': name, 'θpinned': pin,
                                                'θ_tm1': float(θPath.iloc[pos0-1]),
                                                'θ_t0': float(θPath.iloc[pos0]),
                                                'θ_t1': float(θPath.iloc[pos0+1]),
                                                'τ_t0': r0['τ'], 'sr_t0': r0['sr'], 'ww_t0': r0['workweek'],
                                                'τ_t1': r1['τ'], 'sr_t1': r1['sr'], 'ww_t1': r1['workweek']}
                                               | netColumns(mt, out['τ'], θPath, pos0)
                                               | piColumn(mt, mt.t0Year) | solverColumns(a.designRule, a.Ma, mult))
                                print('  {:<9} pin={:<5} θ_t0={:.4f}  τ_t0={:.4f} sr_t0={:.4f} '
                                      'ww_t0={:.2f}  (θ_t1={:.4f}, {:.0f}s)'.format(
                                          name, str(pin), float(θPath.iloc[pos0]), r0['τ'], r0['sr'],
                                          r0['workweek'], float(θPath.iloc[pos0+1]), time.time()-tic))
                                mergeWrite(fShk, shkRows, KEYSHK)
                            except Exception as e:
                                runESC.noteFailure(f'{tag} {name} pin={pin}', e)

                    # France's own calibrated path at this rho -- the endpoint the French-characteristics
                    # rows are read against (runESC.franceRow's CRRA counterpart, exogenous theta). Run with
                    # the French scenarios only, so a run of other scenarios leaves the France rows on disk.
                    if not any(n.startswith('fr') for n in a.scenarios):
                        continue
                    try:
                        tic = time.time()
                        f = franceRowCRRA(m, ρ, spec, phi, pCal, hbarRef, commonX = a.commonX, gs = gs)
                        f['method'], f['host'] = method, a.host
                        f |= solverColumns(a.designRule, a.Ma, f.pop('multiplicity'))
                        shkRows.append(f)
                        print('  {:<9} pin={:<5} θ_t0={:.4f}  τ_t0={:.4f} sr_t0={:.4f} ww_t0={:.2f}'
                              '  [{:.0f}s]'.format('France', 'True', f['θ_t0'], f['τ_t0'], f['sr_t0'],
                                                   f['ww_t0'], time.time()-tic))
                        mergeWrite(fShk, shkRows, KEYSHK)
                    except Exception as e:
                        runESC.noteFailure(f'{tag} France', e)
    print('\n-> {}'.format(os.path.relpath(OUTDIR, REPO)))
    if runESC.FAILURES:
        print('{} FAILED:'.format(len(runESC.FAILURES)) + ''.join('\n  ' + x for x in runESC.FAILURES))
    return 1 if runESC.FAILURES else 0


def setThetaCRRA(a):
    """ runESC.setTheta for this driver's --theta/--root/--main, and the output directory (OUTDIR without a
    variant). Under a variant: stages calib, path, shocks only, the U.S. sweep at the variant at every rho of
    --rho (oecd26.usBeta), and for --host UK the UK's sweep rows too. """
    import oecd26
    if a.theta in (None, oecd26.INDEX):
        runESC.setTheta(None, a.root, main = a.main)        # SystemExit for --root or --main alone
        a.theta = oecd26.indexTheta(a.theta, [OUTDIR], 'runESCcrra.py')
        return OUTDIR
    if set(a.stage) - {'calib', 'path', 'shocks'}:
        raise SystemExit('--theta: stages calib, path, shocks only')
    out = runESC.setTheta(a.theta, a.root, ('calib',), False, a.commonX, a.rho[0], a.main)
    try:
        for ρ in a.rho:
            oecd26.usBeta(runESC.usSweep(a.commonX), ρ, thetaSource = a.theta)
            if a.host != 'US':
                r = hostRow(ρ, a.commonX, a.host, runESC.THETA['calibDir'])
                if r.get('thetaSource') != a.theta or not abs(float(r.get('θ0', np.nan))
                                                            - runESC.countryθ0(a.host)) <= 1e-12:
                    raise ValueError('the {} sweep at rho = {} is not at variant {!r}'.format(a.host, ρ, a.theta))
    except (FileNotFoundError, ValueError, KeyError) as e:
        raise SystemExit('--theta {}: {}'.format(a.theta, e))
    return out


if __name__ == '__main__':
    sys.exit(main())

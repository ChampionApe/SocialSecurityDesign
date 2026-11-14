r""" The endogenous design under the 'size' cost at other Frisch elasticities: log preferences (rho = 1),
the model recalibrated at each xi.

Run:  .venv\Scripts\python.exe python\US\runESCxi.py --theta C --main      # xi = 0.3 and the self-check, then 0.2, 0.4
      ... --theta C --main --xi 0.2 0.4 --stage calib path acute      # no self-check
      ... --theta C --main --xi 0.3 --stage check                     # the self-check on the rows on disk

xi is the workbook's 'Labor supply elasticity' (data/USMain_test.xlsx, 0.3), replaced in the parameters the
model is built from (buildUSxi); the workbook is never written. The cost is Eq (esc:AB) of
writing/US/model_esc.tex, f = exp(-lambda tau Vtilde (1-theta)^2 / 2) with Vtilde of Eq (esc:Vtilde), lambda in
the 'p' slot and column, phi the dummy 0.5; the calibration variant is config.US['commonX'] (common X by
default, --no-commonX for vector X). Stages, per xi:

  calib   lambda such that the design in force in 2020 on the freely simulated leaded path is theta*
          (Eq esc:calibration of num_esc.tex; ModelESC.calibrateWedge -> leadedDesignAtT0: a log-grid scan of
          ModelESC.WEDGE_BRACKET['size'], then a bracketed root), with (beta, omega) recalibrated to
          (R_t0, tau_t0) and X to hbar_t0 at every trial lambda (Eqs calibration, calibration:Xsolve).
          theta* is read off the replacement-rate ratio and relative incomes (Eqs esc:RR, calibration:theta)
          and is asserted equal to the workbook-xi model's, before and after calibration. Row: lambda,
          theta*, beta, omega, X, Vtilde, f(theta*, tau_t0), f(0, tau_t0), tau_t0, the (R, tau) target drift
          of the endogenous path, and the scan as 'lambda:residual' pairs.
  path    the leaded equilibrium (Alg esc:leaded, Eq esc:leadedObjective) at the calibrated lambda, the choice
          binding from the first period, as runESC.stagePath: theta_t, tau_t, savings rate, workweek and R per
          date, beside the same at the exogenous design theta* (columns *Exo).
  acute   nu_t = 1 over the whole horizon (shocks.shockAgeing), a new equilibrium path read at 2020
          (runESC.leadedNewPath, as runESC.stageShocks): theta pinned at theta*, and chosen.
  check   the rows at the workbook xi against the published ones (escCalibration.csv, escPath.csv,
          escExperiments.csv; read only), CHECK_TOL relative for lambda, absolute otherwise. With the workbook
          xi in --xi it runs first and a failure stops the run before any other xi: rows from a code path that
          does not reproduce the headline are not comparable with it.

--theta VARIANT --main are required (runESC.setTheta; VARIANT = python/paper/config.py's US['theta']): every U.S.
model imposes theta = the variant's U.S. row of data/oecd26_theta.csv as db['θ0'], seeded off results/calibration/'s
sweep at that θ0, so theta* is the same at every xi; the published rows the self-check reads (escCalibration.csv,
escPath.csv, escExperiments.csv) and the rows of this file are kept at that thetaSource and θ0 (oecd26.rowsAtTheta), and
every row written carries thetaSource, θ0 and grid (runESC.thetaColumns).

Writes results/esc/escXiRobustness.csv, long format, merged through runESC.mergeWrite on KEYXI after dropping the rows
on file at another thetaSource or θ0 (a csv holds one θ): per xi one
'calibration' row and one 'acute pinned' and one 'acute chosen' row at 2020, and one 'path' row per date.
The path and acute stages read lambda from the calibration row on disk when calib is not run with them;
every row carries its lambda in 'p', and the summary flags a row whose 'p' differs from its xi's calibration
row (a stale row after a recalibration, crossCuttingFindings #13). A xi whose scan finds no sign change
writes its calibration row (converged False, the scan) and skips path and acute; a design at 0 or 1 is a
result and is written as one.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                     # --out resolves against the caller's directory, not HERE
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import test as testmod
import oecd26
import runESC
from policy import multiplicityColumns

ESCDIR = os.path.join(REPO, 'results', 'esc')
OUT = os.path.join(ESCDIR, 'escXiRobustness.csv')
PUBLISHED = {'calib': os.path.join(ESCDIR, 'escCalibration.csv'),
             'path': os.path.join(ESCDIR, 'escPath.csv'),
             'experiments': os.path.join(ESCDIR, 'escExperiments.csv')}
XI0 = float(testmod.pars['ξ'])        # the workbook's value
SPEC, PHI = 'size', 0.5               # the paper's cost; phi is a dummy key under 'size'
KEYXI = ['xi', 'spec', 'phi', 'commonX', 'kind', 'pos']
THETA_TOL = 1e-10                     # theta* at xi against the workbook-xi model
CHECK_TOL = 1e-6


def buildUSxi(wedge, ξ, commonX = True, ρ = 1.0):
    """ runESC.buildUS with the labour supply elasticity set to ξ in the workbook parameters it builds from.

    buildUS reads testmod.pars at call time, so ξ is swapped there for the duration of the call and restored.
    It has to reach the constructor: η_i and X_i are derived from db['ξ'] there (initProductivity; Eq
    calibration:etaCommonX under common X, where _calPostRoot re-derives them after every calibration), so
    writing db['ξ'] on a built model would leave η_i at the workbook ξ. """
    saved = testmod.pars['ξ']
    testmod.pars['ξ'] = float(ξ)
    try:
        m = runESC.buildUS(wedge, ρ = ρ, commonX = commonX)
    finally:
        testmod.pars['ξ'] = saved
    if not np.allclose(m.db['ξ'].values.astype(float), ξ, rtol = 0, atol = 0):
        raise RuntimeError(f'buildUSxi: db[ξ] is not {ξ} at every date.')
    return m


def θStarAt(ξ, commonX):
    """ θ* (Eq calibration:theta) on the uncalibrated model at ξ. """
    m = buildUSxi(None, ξ, commonX)
    return float(m.db['θ'].xs(m.t0Year))


def _assertθ(θ, θ0, ξ, where):
    if not abs(θ - θ0) <= THETA_TOL:
        raise AssertionError(f'θ* at ξ = {ξ} ({where}) is {θ!r}, the workbook-ξ model has {θ0!r}: '
                             f'gap {θ-θ0:+.3e} > {THETA_TOL}.')


def _ids(ξ, commonX, kind, pos, date, λ):
    return {'xi': ξ, 'spec': SPEC, 'phi': PHI, 'commonX': commonX, 'kind': kind, 'pos': int(pos),
            'date': date, 'p': λ} | runESC.thetaColumns()


def atθ(df):
    """ The rows of csv frame `df` at the run's θ (oecd26.rowsAtTheta at runESC.THETA's variant and θ0). """
    return oecd26.rowsAtTheta(df, runESC.THETA['variant'], runESC.THETA['θ0'])


def write(out, rows):
    """ runESC.mergeWrite on KEYXI, after the rows of `out` at another thetaSource or θ0 than the run's are dropped. """
    if os.path.exists(out):
        old = pd.read_csv(out, float_precision = 'round_trip')
        keep = atθ(old)
        if len(keep) < len(old):
            print(f'  {os.path.relpath(out, REPO)}: {len(old) - len(keep)} rows at another θ dropped')
            runESC.writeCsv(keep, out)
    runESC.mergeWrite(out, rows, KEYXI)


# ---------------------------------------------------------------- stages

def stageCalib(ξ, out, commonX, θ0, bracket = None):
    """ The wedge calibration at ξ; returns its row. θ0: θ* of the workbook-ξ model. """
    tic = time.time()
    m = buildUSxi({'spec': SPEC, 'phi': PHI, 'p': 1.0}, ξ, commonX)
    t0, pos = m.t0Year, m.db['t0']
    θξ = float(m.db['θ'].xs(t0))
    _assertθ(θξ, θ0, ξ, 'built')
    print(f'\n[ξ={ξ}] θ* = {θξ:.12f} (workbook ξ: {θ0:.12f}, gap {θξ-θ0:+.1e}); calibrating λ ...')
    try:
        rec = m.calibrateWedge(spec = SPEC, phi = PHI, bracket = bracket)
    except Exception as e:
        rec = {'p': np.nan, 'residual': np.nan, 'scan': [], 'converged': False,
               'message': f'{type(e).__name__}: {e}'}
    v = np.array([s['residual'] for s in rec['scan']], dtype = float)
    ok = np.isfinite(v)
    nSign = int(sum(1 for k in range(len(v)-1) if ok[k] and ok[k+1] and v[k]*v[k+1] < 0))
    r = _ids(ξ, commonX, 'calibration', pos, m.db['dates'][pos], rec['p']) | {
        'converged': bool(rec['converged']), 'θStar': θξ, 'residual': rec['residual'],
        'nSignChanges': nSign, 'Vtilde': float(m.B.Vtilde(t0))}
    if rec['converged']:
        _assertθ(float(m.db['θ'].xs(t0)), θ0, ξ, 'calibrated')
        led = m.solveLeaded(pinAtT0 = False)
        r |= {'β': m.simpleβinv(), 'ω': float(m.db['ω'].xs(t0)),
              'X': float(m.db['Xi'].xs(t0).iloc[0]) if commonX else np.nan,
              'τDrift': led['targetDrift']['τ'], 'RDrift': led['targetDrift']['R'],
              'choice': float(led['θ'].iloc[pos]), 'choiceAtT0': float(led['θ'].iloc[pos+1])}
        r |= runESC.wedgeReadout(m, float(led['τ'].xs(t0)))
        r |= multiplicityColumns(led['multiplicity'])
        print('  -> λ={p:.6f}  θ*={θStar:.4f}  β={β:.4f} ω={ω:.4f} X={X:.4f}  Ṽ={Vtilde:.4f} f(θ*)={fStar:.4f} '
              'f(0)={f0:.4f}  sign changes {nSignChanges}  ({s:.0f}s)'.format(s = time.time()-tic, **r))
    r |= {'message': rec['message'],
          'scan': ';'.join('{:.4f}:{:+.5f}'.format(s['p'], s['residual']) for s in rec['scan'])}
    if not rec['converged']:
        r |= multiplicityColumns(None)
    if not rec['converged']:
        print(f"  -> NOT CONVERGED ({r['message']}); scan: {r['scan']}")
    write(out, [r])
    return r


def stagePath(ξ, λ, out, commonX):
    """ runESC.stagePath at ξ: the leaded path and the exogenous-θ path, every dated position. """
    tic = time.time()
    m = buildUSxi({'spec': SPEC, 'phi': PHI, 'p': λ}, ξ, commonX)
    m.calibrate()
    base, hbarRef = runESC.baselineRefs(m)
    led = m.solveLeaded(pinAtT0 = False)
    dates = m.db['dates']
    rows = []
    for pos, t in enumerate(m.db['t']):
        if pos + 1 >= len(m.db['t']):
            break
        r = runESC.readout(m, led['τ'], led['report'], hbarRef, pos)
        rb = runESC.readout(m, base['τ'], base['report'], hbarRef, pos)
        rows.append(_ids(ξ, commonX, 'path', pos, dates[pos] if pos < len(dates) else np.nan, λ) | {
            'ν': float(m.db['ν'].xs(t)), 'θ': float(led['θ'].xs(t)), 'τ': r['τ'], 'sr': r['sr'],
            'workweek': r['workweek'], 'R': r['R'],
            'τExo': rb['τ'], 'srExo': rb['sr'], 'workweekExo': rb['workweek']} | multiplicityColumns(led['multiplicity']))
    write(out, rows)
    print('[ξ={}] θ path: {}   τ path: {}  ({:.0f}s)'.format(
        ξ, '  '.join('{:.4f}'.format(x) for x in led['θ'].values[:8]),
        '  '.join('{:.4f}'.format(x) for x in led['τ'].values[:8]), time.time()-tic))


def stageAcute(ξ, λ, out, commonX):
    """ runESC.stageShocks' acute-ageing scenario at ξ, pinned and chosen, read at t0 (and t0 ± 1). """
    tic = time.time()
    m = buildUSxi({'spec': SPEC, 'phi': PHI, 'p': λ}, ξ, commonX)
    m.calibrate()
    _, hbarRef = runESC.baselineRefs(m)
    apply = runESC.SHOCKS_ESC['acute'][1]
    pos = m.db['t0']
    rows = []
    for pin in (True, False):
        r = runESC.leadedNewPath(m, hbarRef, apply, None, pin = pin)
        mt = r['m']
        rows.append(_ids(ξ, commonX, 'acute ' + ('pinned' if pin else 'chosen'), pos, m.db['dates'][pos], λ) | {
            'ν': float(mt.db['ν'].xs(mt.db['t'][pos])), 'θ': r['θ0'], 'τ': r['t0']['τ'],
            'sr': r['t0']['sr'], 'workweek': r['t0']['workweek'], 'R': r['t0']['R'],
            'θ_tm1': r['θ_'], 'θ_t1': r['θ1'], 'τ_t1': r['t1']['τ'], 'sr_t1': r['t1']['sr'],
            'ww_t1': r['t1']['workweek']} | multiplicityColumns(r['multiplicity']))
        print('[ξ={}] acute {:<6}: θ_t0={:.4f}  τ_t0={:.4f} sr_t0={:.4f} ww_t0={:.2f}  (θ_t1={:.4f})'.format(
            ξ, 'pinned' if pin else 'chosen', r['θ0'], r['t0']['τ'], r['t0']['sr'], r['t0']['workweek'], r['θ1']))
    write(out, rows)
    print(f'  ({time.time()-tic:.0f}s)')


# ---------------------------------------------------------------- the self-check and the summary

def _rows(out, ξ, commonX):
    """ This file's rows at ξ under (SPEC, PHI, commonX) and the run's θ; an empty frame if none. """
    if not os.path.exists(out):
        return pd.DataFrame()
    d = atθ(pd.read_csv(out))
    return d[np.isclose(d['xi'], ξ) & (d['spec'] == SPEC) & np.isclose(d['phi'], PHI)
             & (d['commonX'].astype(bool) == bool(commonX))]


def calibrationRow(out, ξ, commonX):
    d = _rows(out, ξ, commonX)
    d = d[d['kind'] == 'calibration'] if len(d) else d
    return None if d.empty else d.iloc[0].to_dict()


def selfCheck(out, commonX, tol = CHECK_TOL):
    """ The workbook-ξ rows of `out` against the published ones, both at the run's θ (atθ). True iff at least one
    pair was compared and every gap is within tol (relative for λ, absolute otherwise). Pairs whose rows are
    missing on either side are listed, not failed. """
    d = _rows(out, XI0, commonX)
    cal = runESC.pickCalib(atθ(pd.read_csv(PUBLISHED['calib'])), SPEC, PHI, commonX)
    path = atθ(pd.read_csv(PUBLISHED['path']))
    path = path[(path['spec'] == SPEC) & np.isclose(path['phi'], PHI) & (path['commonX'].astype(bool) == bool(commonX))]
    exp = atθ(pd.read_csv(PUBLISHED['experiments']))
    exp = exp[np.isclose(exp['ρ'], 1.0) & (exp['preferences'] == 'LOG') & (exp['spec'] == SPEC)
              & np.isclose(exp['phi'], PHI) & (exp['commonX'].astype(bool) == bool(commonX))
              & (exp['scenario'] == 'acute')]
    pairs, missing = [], []

    def pair(label, mine, pub, rel = False):
        if len(mine) == 0 or len(pub) == 0:
            missing.append(label)
            return
        a, b = float(mine.iloc[0]), float(pub.iloc[0])
        gap = (a - b)/abs(b) if rel else a - b
        pairs.append((label, a, b, gap, abs(gap) <= tol))

    c = d[d['kind'] == 'calibration'] if len(d) else d
    for k in ('p', 'θStar', 'β', 'ω', 'X', 'Vtilde', 'fStar', 'f0', 'choice'):
        pair(f'calibration {k}', c[k] if (len(c) and k in c) else [], cal[k] if k in cal else [], rel = k == 'p')
    for pos in (2, 3, 4, 5):
        mine = d[(d['kind'] == 'path') & (d['pos'] == pos)] if len(d) else d
        pub = path[path['pos'] == pos]
        for k in ('θ', 'τ'):
            pair(f'path {k} {int(pub["date"].iloc[0]) if len(pub) else pos}', mine[k] if len(mine) else [],
                 pub[k] if len(pub) else [])
    for pin in (True, False):
        kind = 'acute ' + ('pinned' if pin else 'chosen')
        mine = d[d['kind'] == kind] if len(d) else d
        pub = exp[exp['θpinned'].astype(bool) == pin]
        for k, kp in (('θ', 'θ_t0'), ('τ', 'τ_t0')):
            pair(f'{kind} {k} 2020', mine[k] if len(mine) else [], pub[kp] if len(pub) else [])

    print(f'\nSELF-CHECK at the workbook ξ = {XI0} against results/esc/escCalibration, escPath, escExperiments '
          f'(tol {tol:g}):')
    for label, a, b, gap, good in pairs:
        print('  {:<24} {:>18.12f} {:>18.12f}  gap {:+.2e}  {}'.format(label, a, b, gap, 'ok' if good else 'FAIL'))
    if missing:
        print('  not compared (row missing): ' + ', '.join(missing))
    passed = bool(pairs) and all(p[-1] for p in pairs)
    print('SELF-CHECK ' + ('PASSED' if passed else 'FAILED'))
    return passed


def summary(out, xis, commonX):
    """ One block per ξ from the csv; flags rows whose λ differs from their calibration row's. """
    print('\nSUMMARY ({})'.format(os.path.relpath(out, REPO)))
    for ξ in xis:
        d = _rows(out, ξ, commonX)
        c = d[d['kind'] == 'calibration'] if len(d) else d
        if c.empty:
            print(f'  ξ={ξ}: no calibration row')
            continue
        c = c.iloc[0]
        if not bool(c['converged']):
            print(f"  ξ={ξ}: λ NOT calibrated ({c['message']}); scan {c['scan']}")
            continue
        stale = d[~np.isclose(d['p'], c['p'], rtol = 0, atol = 0)]
        print('  ξ={}: λ={:.4f} θ*={:.10f} Ṽ={:.4f} f(θ*,τ2020)={:.4f} f(0,τ2020)={:.4f} τ2020={:.4f} '
              'X={:.4f} β={:.4f} ω={:.4f} sign changes {:d} RDrift={:.1e} τDrift={:.1e}'.format(
                  ξ, c['p'], c['θStar'], c['Vtilde'], c['fStar'], c['f0'], c['τ0'], c['X'], c['β'], c['ω'],
                  int(c['nSignChanges']), c['RDrift'], c['τDrift']))
        pth = d[d['kind'] == 'path'].sort_values('pos')
        if len(pth):
            sel = pth[pth['date'].isin([1990, 2020, 2050, 2080, 2110])]
            print('        path  ' + '  '.join('{:.0f}: θ={:.4f} τ={:.4f} (τExo={:.4f})'.format(
                r['date'], r['θ'], r['τ'], r['τExo']) for _, r in sel.iterrows()))
        for kind in ('acute pinned', 'acute chosen'):
            a = d[d['kind'] == kind]
            if len(a):
                a = a.iloc[0]
                print('        {:<12}  2020: θ={:.4f} τ={:.4f} sr={:.4f} ww={:.2f}   (θ 1990={:.4f}, 2050={:.4f})'.format(
                    kind, a['θ'], a['τ'], a['sr'], a['workweek'], a['θ_tm1'], a['θ_t1']))
        if len(stale):
            print('        STALE: {} row(s) carry a λ other than the calibration row\'s ({})'.format(
                len(stale), ', '.join(sorted(set(stale['kind'])))))


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description = 'The endogenous design under the size cost at other Frisch '
                                               'elasticities (rho = 1).')
    ap.add_argument('--xi', type = float, nargs = '+', default = [XI0, 0.2, 0.4],
                    help = 'the Frisch elasticities; the workbook value runs first, then the check')
    ap.add_argument('--stage', nargs = '+', default = ['calib', 'path', 'acute', 'check'],
                    choices = ('calib', 'path', 'acute', 'check'))
    ap.add_argument('--commonX', action = argparse.BooleanOptionalAction, default = True,
                    help = "the calibration variant, config.US['commonX'] (common X, the headline); part of the "
                           "merge key")
    ap.add_argument('--bracket', type = float, nargs = 2, default = None,
                    help = "the λ scan bracket; default ModelESC.WEDGE_BRACKET['size']")
    ap.add_argument('--out', default = OUT)
    ap.add_argument('--theta', default = None, help = "the pipeline's main θ variant, config.US['theta'] (with --main)")
    ap.add_argument('--main', action = 'store_true', help = "with --theta: read and write results/esc/")
    a = ap.parse_args()
    a.out = os.path.abspath(os.path.join(CWD, a.out))
    if a.theta is None or not a.main:
        raise SystemExit('--theta VARIANT --main are required: the rows are at config.US[\'theta\'], the design the '
                         'paper reads')
    runESC.setTheta(a.theta, None, ('calib',), False, a.commonX, 1.0, True)
    xis = [round(float(x), 10) for x in a.xi]
    xis = [x for x in xis if np.isclose(x, XI0)][:1] + [x for x in dict.fromkeys(xis) if not np.isclose(x, XI0)]
    run = [s for s in ('calib', 'path', 'acute') if s in a.stage]
    os.makedirs(os.path.dirname(a.out), exist_ok = True)
    print('runESCxi: ξ = {}, stages {}, commonX = {}, -> {}'.format(xis, a.stage, a.commonX,
                                                                  os.path.relpath(a.out, REPO)))
    θ0 = θStarAt(XI0, a.commonX)
    for ξ in xis:
        cal = None
        if 'calib' in run:
            cal = stageCalib(ξ, a.out, a.commonX, θ0, bracket = None if a.bracket is None else tuple(a.bracket))
        elif run:
            cal = calibrationRow(a.out, ξ, a.commonX)
        if run and (cal is None or not bool(cal['converged'])):
            print(f'[ξ={ξ}] no converged λ -- path and acute skipped')
        elif cal is not None:
            if 'path' in run:
                stagePath(ξ, float(cal['p']), a.out, a.commonX)
            if 'acute' in run:
                stageAcute(ξ, float(cal['p']), a.out, a.commonX)
        if 'check' in a.stage and np.isclose(ξ, XI0):
            if not selfCheck(a.out, a.commonX):
                summary(a.out, xis, a.commonX)
                print('DONE (self-check failed)')
                raise SystemExit(f'self-check FAILED at ξ = {XI0}: stopped before {xis[1:]}.')
    if 'check' in a.stage and not any(np.isclose(x, XI0) for x in xis):
        print(f'check skipped: the workbook ξ = {XI0} is not in --xi')
    summary(a.out, xis, a.commonX)
    print('DONE')
    return 0


if __name__ == '__main__':
    sys.exit(main())

r""" The U.S. tax of 2020 against a design pinned on the whole path, with and without the cost of redistribution.

Run:  .venv\Scripts\python.exe python\US\runTauTheta.py --theta C --main --commonX
      ... --theta index                  the workbook's θ into results/esc/ (without --theta: refused, oecd26.indexTheta)
      ... --grid 0 0.1 ... 1            the pinned designs (default 0, 0.1, ..., 1; θ0 is always added)
      ... --refine 0.01                 also the designs within one grid cell of the cost model's grid maximum, at this
                                        step (0 for none; default 0.01)
      ... --force                       recompute rows already on disk

Two U.S. economies at ρ = 1 (log preferences), each calibrated at its observed design θ0 (the U.S. row of
data/oecd26_theta.csv under --theta, runESC.setTheta):
  'none'  without the cost: section 6's calibration (β, ω, X recalibrated as in runESC's calib stage 'none' row);
  'size'  under the paper's cost eq (esc:AB), f(θ, τ) = exp(-½ λ τ Ṽ (1-θ)²), at the λ of escCalibration.csv
          (spec 'size', phi 0.5, the run's commonX and θ variant), (β, ω, X) recalibrated at that λ.
In each, θ_t = θ for every date (shocks.shockTheta) and the LOG politico-economic equilibrium with θ given is solved
from the economy's own steady state (runESC.leadedNewPath, pin = True: the tax first order condition of
eq (PEELOG:decoupling), with the cost's factor eq (esc:dlnc2i) under 'size'); β, ω, X and λ stay at their calibrated
values. Each row reports at 2020: τ, spending 100 τ (1 - α) in % of GDP, the savings rate, the workweek against the
economy's own baseline, R, and the tax rule's counts (policy.multiplicityColumns), with the economy's calibrated β and
ω, on which a row's reuse is keyed.

Writes results/esc/escTauTheta.csv (under --main; a variant's root/esc otherwise), merged on (spec, phi, commonX,
θpin) with thetaColumns; a row on disk at the same key, θ variant, λ, β and ω is reused unless --force.
Self-checks printed at the end: the 'none' rows at θ = 0 and 1 against results/shocks/US_shocksCommonX.csv's θ = 0 and
θ = 1 rows at ρ = 1, and the 'size' row at θ0 against escShocks.csv's pinned baseline; a gap above SELFCHECKTOL in τ,
or no pair to compare, exits with an error after the rows are written. A relative --root is the caller's working
directory's.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()                     # where a relative --root points: runESC's import changes the directory to HERE
sys.path.insert(0, HERE)
import runESC as R                    # changes the working directory to HERE
import shocks as sh
from policy import multiplicityColumns

KEY = ['spec', 'phi', 'commonX', 'θpin']
SELFCHECKTOL = 1e-6                   # the largest gap in τ against the published rows
PHI = 0.5
GRID = [round(0.1*k, 10) for k in range(11)]


def _onDisk(path, spec, θ, commonX, p, cal):
    """ The row of `path` at (spec, θ, commonX), the run's θ variant, cost parameter `p` and calibration `cal`
    ({'β', 'ω'} of the economy just calibrated: a row of an earlier calibration is recomputed, #13), or None. """
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path, float_precision = 'round_trip')
    want = R.thetaColumns()
    hit = d[(d['spec'] == spec) & np.isclose(d['θpin'], θ, rtol = 0, atol = 1e-12)
            & (d['commonX'].astype(bool) == bool(commonX))]
    for k, v in want.items():
        hit = hit[(np.isclose(pd.to_numeric(hit[k]), v, rtol = 0, atol = 1e-12)) if k == 'θ0' else (hit[k] == v)]
    if spec != 'none':
        hit = hit[np.isclose(hit['p'], p, rtol = 1e-12, atol = 0)]
    for k, v in cal.items():
        if k not in hit.columns:
            return None
        hit = hit[np.isclose(pd.to_numeric(hit[k]), v, rtol = 1e-10, atol = 0)]
    return None if hit.empty else hit.iloc[-1].to_dict()


def economy(spec, calib, commonX):
    """ (model, λ, hbarRef, {'β', 'ω'}): the U.S. calibrated at θ0 without the cost ('none') or at the calibrated λ
    ('size'). """
    if spec == 'none':
        m, p = R.buildUS(None, commonX = commonX), np.nan
    else:
        hit = R.pickCalib(calib, spec, PHI, commonX)
        if hit.empty:
            raise SystemExit('no converged {} calibration in escCalibration.csv at the run\'s θ variant'.format(spec))
        p = float(hit.iloc[0]['p'])
        m = R.buildUS({'spec': spec, 'phi': PHI, 'p': p}, commonX = commonX)
    m.calibrate()
    _, hbarRef = R.baselineRefs(m)
    return m, p, hbarRef, {'β': float(m.simpleβinv()), 'ω': float(m.db['ω'].xs(m.t0Year))}


def pinned(m, θ, hbarRef):
    """ The readouts at 2020 of the path with θ_t = θ at every date. """
    r = R.leadedNewPath(m, hbarRef, apply = lambda mt, d: sh.shockTheta(mt, θ), pin = True)
    t0 = r['m'].t0Year
    α = float(r['m'].db['α'].xs(t0))
    rt = r['t0']
    return {'τ': rt['τ'], 'spending': 100*rt['τ']*(1 - α), 'sr': rt['sr'], 'workweek': rt['workweek'], 'R': rt['R'],
            'θ_t0': r['θ0']} | multiplicityColumns(r['multiplicity'])


def selfCheck(df, outdir, commonX, tol = SELFCHECKTOL):
    """ The rows that reproduce published numbers: 'none' at θ = 0, 1 against the shock csv, 'size' at θ0 against
    escShocks.csv's pinned baseline. Printed; returns the largest absolute gap in τ, and raises RuntimeError when
    it exceeds `tol` or when no pair could be compared (NaN). """
    sys.path.insert(0, os.path.join(R.REPO, 'python', 'paper'))
    import config as PC
    gaps = []
    shk = PC.atTheta(pd.read_csv(os.path.join(R.REPO, 'results', 'shocks',
                                              'US_shocksCommonX.csv' if commonX else 'US_shocks.csv')))
    shk = shk[np.isclose(shk['ρ'], 1.0) & (shk['effect'] == 'full')]
    for θ, scen in ((0., r'$\theta = 0$'), (1., r'$\theta = 1$')):
        mine = df[(df['spec'] == 'none') & np.isclose(df['θpin'], θ)]
        ref = shk[shk['scenario'] == scen]
        if len(mine) and len(ref):
            g = float(mine['τ'].iloc[0]) - float(ref['τ'].iloc[0])
            gaps.append(abs(g))
            print('self-check none θ = {:g}: τ {:.10f} against US_shocks {:.10f} (gap {:.2e})'.format(
                θ, float(mine['τ'].iloc[0]), float(ref['τ'].iloc[0]), g))
    esc = os.path.join(outdir, 'escShocks.csv')
    if os.path.exists(esc):
        e = pd.read_csv(esc)
        if R.THETA['variant'] is not None:
            e = e[(e.get('thetaSource') == R.THETA['variant'])]
        e = e[(e['spec'] == 'size') & (e['scenario'] == 'baseline') & (e['θpinned'].astype(str) == 'True')
              & (e['commonX'].astype(bool) == bool(commonX))]
        mine = df[(df['spec'] == 'size') & np.isclose(df['θpin'], float(R.THETA['θ0'] or np.nan))]
        if len(e) and len(mine):
            g = float(mine['τ'].iloc[0]) - float(e['τ_t0'].iloc[0])
            gaps.append(abs(g))
            print('self-check size θ0: τ {:.10f} against escShocks {:.10f} (gap {:.2e})'.format(
                float(mine['τ'].iloc[0]), float(e['τ_t0'].iloc[0]), g))
    none = df[df['spec'] == 'none'].set_index('θpin')['τ']
    if 0.0 in none.index and 1.0 in none.index:
        print('without the cost, τ(1) - τ(0) = {:+.4f} p.p.'.format(100*(none[1.0] - none[0.0])))
    gap = max(gaps) if gaps else np.nan
    if not gap <= tol:
        raise RuntimeError('self-check: largest gap in τ against the published rows is {} (tolerance {:g}; NaN: no '
                           'row to compare)'.format(gap, tol))
    return gap


def main():
    p = argparse.ArgumentParser(description = __doc__.split('\n')[0])
    p.add_argument('--theta', default = None)
    p.add_argument('--main', action = 'store_true')
    p.add_argument('--root', default = None)
    p.add_argument('--commonX', action = 'store_true')
    p.add_argument('--grid', type = float, nargs = '+', default = GRID)
    p.add_argument('--refine', type = float, default = 0.01)
    p.add_argument('--specs', nargs = '+', default = ['none', 'size'], choices = ('none', 'size'))
    p.add_argument('--force', action = 'store_true')
    a = p.parse_args()
    a.root = None if a.root is None else os.path.abspath(os.path.join(CWD, a.root))
    import oecd26
    theta = a.theta if a.root is not None else oecd26.indexTheta(a.theta, [R.OUTDIR], 'runTauTheta.py')
    outdir = R.setTheta(theta, a.root, ('calib',), False, a.commonX, 1.0, a.main)
    out = os.path.join(outdir, 'escTauTheta.csv')
    calibPath = os.path.join(outdir, 'escCalibration.csv')
    calib = R.checkCalib(pd.read_csv(calibPath), calibPath, ['size'], [PHI], a.commonX)
    θ0 = R.THETA['θ0']
    rows = []
    for spec in a.specs:
        tic = time.time()
        m, lam, hbarRef, cal = economy(spec, calib, a.commonX)
        print('[{}] calibrated at θ0 = {}, λ = {}, β = {:.6f}, ω = {:.6f} ({:.0f}s)'.format(
            spec, θ0, lam, cal['β'], cal['ω'], time.time() - tic))
        grid = sorted(set(a.grid) | ({float(θ0)} if θ0 is not None else set()))

        def solve(θ, kind):
            old = None if a.force else _onDisk(out, spec, θ, a.commonX, lam, cal)
            if old is not None:
                rows.append(old)
                return old
            t = time.time()
            r = {'spec': spec, 'phi': PHI if spec != 'none' else np.nan, 'commonX': a.commonX, 'θpin': float(θ),
                 'p': lam, 'kind': kind} | cal | pinned(m, θ, hbarRef) | R.thetaColumns()
            rows.append(r)
            R.mergeWrite(out, [x for x in rows], KEY)
            print('  [{}] θ = {:.4f} ({}): τ = {:.6f}, spending {:.4f}% of GDP, nEqMax {} ({:.1f}s)'.format(
                spec, θ, kind, r['τ'], r['spending'], r['nEqMax'], time.time() - t))
            return r

        res = {θ: solve(θ, 'grid') for θ in grid}
        if spec != 'none' and a.refine > 0:
            on = sorted(x for x in res if x in a.grid)
            τs = [res[x]['τ'] for x in on]
            k = int(np.argmax(τs))
            if 0 < k < len(on) - 1:
                lo, hi = on[k - 1], on[k + 1]
                n = int(round((hi - lo)/a.refine))
                for θ in np.round(lo + a.refine*np.arange(1, n), 10):
                    if not any(abs(θ - x) < 1e-12 for x in res):
                        res[float(θ)] = solve(float(θ), 'refine')
    df = pd.read_csv(out, float_precision = 'round_trip')
    if R.THETA['variant'] is not None:
        df = df[df['thetaSource'] == R.THETA['variant']]
    df = df[df['commonX'].astype(bool) == bool(a.commonX)]
    gap = selfCheck(df, outdir, a.commonX)
    for spec in a.specs:
        s = df[df['spec'] == spec].sort_values('θpin')
        k = s['τ'].idxmax()
        print('[{}] max τ {:.6f} at θ = {:.4f}; τ(0) {:.6f}, τ(1) {:.6f}'.format(
            spec, s.at[k, 'τ'], s.at[k, 'θpin'], s['τ'].iloc[0], s['τ'].iloc[-1]))
    print('\n-> {}  (largest self-check gap {:.2e})'.format(os.path.relpath(out, R.REPO), gap))
    return 0


if __name__ == '__main__':
    sys.exit(main())

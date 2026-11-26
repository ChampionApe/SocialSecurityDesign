r""" Run the US counterfactuals at one or more rho and write the three paper tables' worth of numbers.

Run:  .venv\Scripts\python.exe python\US\runShocksUS.py                      # rho = 0.5, 1, 2
      ... --rho 1                                the LOG case alone
      ... --family theta                         one family (theta | ageing | french | all)
      ... --commonX                              the common-X calibration variant
      ... --theta VARIANT [--root DIR] [--coarse]  the US at theta imposed (thetaRun)
      ... --theta VARIANT --main                 the pipeline's main variant into results/shocks/

Three families. Each scenario is a NEW EQUILIBRIUM PATH -- the shocked parameters hold over the whole
1960-2200 horizon and the economy starts at its own steady state, not the US's -- read at 2020, and
reported as a full effect (tau re-optimised) and an economic-equilibrium-only effect (tau held at the
baseline path):

    theta    theta = 0 and theta = 1
    ageing   mild (nu -> (1+nu)/2) and acute (nu -> 1), throughout
    french   France's income distribution, leisure preferences, voting, and all three at once

Read by the paper's figure 4 (python/paper/figuresUS.usOverview) and the supplemental appendix's counterfactuals
(python/paper/siteData.oecdCounterfactuals), which also carries the economic-equilibrium effect.

The income-distribution scenario holds theta at the US design, so that row is about inequality alone and the
theta family carries pension design.

Plus a FRANCE reference row (--noFrance to skip): France's own calibrated path, read at its own 2020 with
the US workweek reference, so the French-characteristics rows can be read against the country they are
borrowed from. That row is the point of the new-path convention -- an unanticipated 2020 reform is not
commensurable with a country's own equilibrium path, a counterfactual country that has always had these
characteristics is. Note the France row's WORKWEEK is a calibration target (ModelFR targets hours
relative to the US), not a prediction; only its tau and savings rate are results.

See shocks.py for what each shock does and for the reporting conventions (savings rate is s/(w*h), the
workweek is expressed relative to the US baseline, and everything is read at db['t0']). The French
characteristics are read from data/FRMain.xlsx via testEU.py -- France's income groups are cut at US
percentiles precisely so that these swaps are well defined.

The calibration at each rho comes from the US sweep csv (results/calibration/US_rhoGrid{,CommonX}.csv),
so this script never re-runs the outer root: it installs (beta, omega) and solves. That also means every
rho asked for must already be in the sweep.

THE θ VARIANT. --theta VARIANT imposes theta on the US as db['θ0'] = the US row of data/oecd26_theta.csv under
VARIANT (ModelUS.getθ, oecd26.thetaFromFile): the baseline and every US-based scenario carry it, the income row
holds theta at it, and the theta family still sets 0 and 1. The calibration is the variant's own sweep,
DIR/calibration/US_rhoGrid{,CommonX}.csv (calibrateRhoGrid.py --theta, at that θ0: oecd26.usBeta), which is also
France's US reference (beta imposed, hours); France takes its θ0 from the same file (the FRA row,
oecd26.countryθ0) in the reference row and in the calibration its characteristics are read from (frenchData), with
omega recalibrated at it. DIR is the variant's directory (oecd26.thetaDir, results/oecd26_theta_<variant>/;
results/ under --main, the pipeline's main variant, oecd26.thetaRoot), the csv DIR/shocks/US_shocks{,CommonX}.csv
with thetaSource, θ0 and grid columns; --coarse (variant only) puts the tax grid at 51 nodes. --theta index is the
workbook's θ (no θ columns); without --theta a run into results/ is refused (oecd26.indexTheta).

Writes results/shocks/US_shocks{,CommonX}.csv -- one row per (rho, family, scenario, effect), each with the
tax counts of its solve (nEqMax, nCandMax, nFallback; policy.multiplicityColumns): the baseline's, France's and
each scenario's full effect (shocks.runOne's 'multiplicity'); -1 on the 'ee' rows, which select no tax. A run
replaces the rows of its own (rho, family, scenario, effect) in the csv and keeps the others (`writeRows`), so a
subset (--rho 1, --family theta) leaves the rest; rows at another θ than the run's (thetaSource, θ0) are dropped, a
csv holds one θ. The file is written through a temporary file and renamed. France's reference row is the one failure
a run survives: it is printed, the run goes on, and main exits 1 (runESC.noteFailure). A relative --out or --root is
the caller's working directory's.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CWD = os.getcwd()                 # where a relative --out or --root points: the chdir below moves to HERE
sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
sys.path.insert(0, HERE)
os.chdir(HERE)

import test as testmod
import testEU
import shocks as sh
import oecd26
import runESC
from policy import multiplicityColumns
from model import ModelUS

OUTDIR = os.path.join(REPO, 'results', 'shocks')

FAMILIES = {'theta': ('theta0', 'theta1'), 'ageing': ('mild', 'acute'),
            'french': ('frIncome', 'frLeisure', 'frVoting', 'frAll')}


def hostRow(ρ, commonX, host = 'US', calibDir = None):
    """ A sweep csv's row at this rho: `host` is the sweep's stem ('US', or 'UK' for runESCcrra.py's UK cost
    calibration); calibDir the sweep's directory, results/calibration by default. """
    csv = os.path.join(calibDir or os.path.join(REPO, 'results', 'calibration'),
                       host + '_rhoGrid' + ('CommonX' if commonX else '') + '.csv')
    df = pd.read_csv(csv)
    hit = df.loc[(df['ρ'] - ρ).abs() < 1e-9]
    if hit.empty:
        raise KeyError(f'No {host} calibration at ρ={ρ} in {csv}. Run the sweep there first.')
    return hit.iloc[-1]


def hostModel(ρ, commonX, row, gs, θ0 = None):
    """ The US's calibrated model at rho, rebuilt from its sweep row rather than unpickled: the csv is the
    record of the calibration. ModelUS with (beta, omega) and, under common X, the calibrated X. The baseline
    solve that follows is the check: tau must come back at the workbook target. θ0: theta imposed (db['θ0'],
    ModelUS.getθ), None for the workbook's. """
    preferences = 'LOG' if ρ == 1 else 'CRRA'
    m = ModelUS(pars = testmod.pars | {'ρ': float(ρ), 'β': float(row['β']), 'ω': float(row['ω'])}
                | ({} if θ0 is None else {'θ0': float(θ0)}), commonX = commonX, **testmod.kwargs)
    m.db['dates'] = testmod.dates
    getattr(m, preferences).initGS(gs)
    if commonX:
        m.initProductivity_commonX(X = float(row['X']))
        m.updateAuxPars()
    return m


def frenchData(m, ρ, preferences, gs = None, commonX = False, usCsv = None, θ0FR = None):
    """ France's characteristics, in the form the shocks want, plus the host quantities they are relative to.
    usCsv: the US sweep France is referenced to instead of results/calibration's (a θ variant's); θ0FR France's
    design imposed (db['θ0'], a θ variant's FRA row), None for its RR0's.

    zetaFR and muFR are data and do not depend on rho. xbarRatio does: France's X_i carry the calibrated
    hours-unit level lambda, which moves across rho (0.8835 -> 0.9020 over the sweep), so France is
    calibrated at the SAME rho -- and the SAME calibration variant -- as the US model being shocked,
    rather than once at rho = 1 under vector X.

    zetaFR is France's relative income on ITS OWN groups -- which are the US percentiles, so gamma matches
    the US vector type by type and the swap is a like-for-like replacement. xbarRatio is the ratio of
    population-weighted mean X_i, France's calibrated value over the US's: the row the paper's calibration
    table reports as `X`, and the only thing a pure scale can be matched on.
    """
    parsFR, kwargsFR, _, _ = testEU.load('FR')
    # commonX must match the host model being shocked. Under commonX, France's X_i is a single scalar and
    # its eta_i inverts z^eta directly, so both characteristics differ from the vector-X calibration --
    # taking them from the wrong variant would impose a France that does not exist in this run's units.
    mFR = testEU.model('FR', ρ = float(ρ), commonX = commonX, usCsv = usCsv, θ0 = θ0FR)
    if gs:
        getattr(mFR, preferences).initGS(gs)
    mFR.calibrate(preferences = preferences)
    t0FR = mFR.db['t'][mFR.db['t0']]
    γFR = mFR.db['γi'].xs(t0FR).values.astype(float)
    XbarFR = float((γFR * mFR.db['Xi'].xs(t0FR).values.astype(float)).sum())

    t0 = m.db['t'][m.db['t0']]
    γUS = m.db['γi'].xs(t0).values.astype(float)
    XbarUS = float((γUS * m.db['Xi'].xs(t0).values.astype(float)).sum())
    # gamma is not bit-identical across the two workbooks -- the same nominal percentile cuts land on
    # slightly different shares in CPS and LIS -- so this reports the gap rather than asserting equality.
    ηFR = mFR.db['ηi'].xs(t0FR).values.astype(float)
    # ηScale: the level the income row puts France's eta at (Gamma_h = 1 at the US X, shocks.ηLevel).
    # The leisure row scales X_i by ηScale*xbarRatio so the two rows compose to France's own (eta, X).
    return {'ηFR': ηFR, 'μFR': parsFR['μj'], 'ηScale': sh.ηLevel(m, ηFR),
            'xbarRatio': XbarFR/XbarUS, 'XbarFR': XbarFR, 'XbarUS': XbarUS,
            'γgap': float(np.max(np.abs(γFR - γUS))),
            'θUS': float(m.db['θ'].xs(t0))}


def effectRows(r, ids):
    """ The two csv rows of one shocks.runOne result: `ids` (the row's identifiers) plus the readout of
    each effect, with the full effect's tax counts (r['multiplicity'], policy.multiplicityColumns) and -1
    on the 'ee' row, whose taxes are held at the baseline path. """
    return [{**ids, 'effect': eff, **r[eff]} | multiplicityColumns(r['multiplicity'] if eff == 'full' else None)
            for eff in ('full', 'ee')]


def franceReference(ρ, preferences, gs = None, commonX = False, workweekData = None, hbarRef = None,
                    usCsv = None, θ0FR = None):
    """ France's own calibrated equilibrium at 2020, in the same three columns the host rows carry.

    The endpoint of the "mostly US, partly France" scale: the supplemental appendix reads the French-characteristics
    rows against this (scenario 'france' of siteData.oecdCounterfactuals). It is NOT a shock on the US model -- France has its own eta, X, mu, nu AND its own
    calibrated omega (ModelFR imposes the US beta and searches omega alone), so the gap between the
    all-characteristics row and this one is exactly what the observable characteristics fail to explain.

    hbarRef is the US baseline's hbar and workweekData the US workbook's own workweek, so the reported
    number is France's hours on the US scale. That is the comparable object under vector X, where the
    LEVEL of hbar is not identified -- and it is also ModelFR's own hours target, so this cell reproduces
    France's observed workweek by construction rather than predicting it. Say so wherever it is printed.
    usCsv, θ0FR: the US sweep France is referenced to and France's imposed design, as frenchData.
    """
    mFR = testEU.model('FR', ρ = float(ρ), commonX = commonX, usCsv = usCsv, θ0 = θ0FR)
    if gs:
        getattr(mFR, preferences).initGS(gs)
    mFR.calibrate(preferences = preferences)
    out = getattr(mFR, f'solvePEE_{preferences}')()
    r = sh.readout(mFR, out['τ'], out['report'], workweekData, hbarRef, pos = mFR.db['t0'])
    r['θ'] = float(mFR.db['θ'].xs(mFR.db['t'][mFR.db['t0']]))
    r['ω'] = float(mFR.db['ω'].xs(mFR.db['t'][mFR.db['t0']]))
    return r | multiplicityColumns(out.get('multiplicity'))


KEYROW = ['ρ', 'family', 'scenario', 'effect']      # the identity of a row of the csv


def writeRows(path, rows, tag):
    """ Csv `path` after this run's `rows`: the rows on disk at the run's θ (oecd26.rowsAtTheta at tag's thetaSource
    and θ0; the workbook's without a variant) whose KEYROW the run did not produce, then `rows`; written atomically
    (runESC.writeCsv). Called after every scenario with the rows so far. """
    new = pd.DataFrame(rows)
    if os.path.exists(path):
        old = pd.read_csv(path, float_precision = 'round_trip')
        keep = oecd26.rowsAtTheta(old, tag.get('thetaSource'), tag.get('θ0'))
        keep = keep[~runESC._keyOf(keep, KEYROW).isin(runESC._keyOf(new, KEYROW))]
        if len(keep):
            new = pd.concat([keep, new], ignore_index = True)
    runESC.writeCsv(new, path)


def thetaRun(theta, root = None, commonX = True, rhos = (1.0,), main = False):
    """ --theta VARIANT: (θ0, the variant's calibration directory, the default csv, France's θ0). root: the
    variant's directory (oecd26.thetaDir by default; a relative path is the repo's; results/ under `main`,
    oecd26.thetaRoot). The sweep at every rho must be at θ0 (oecd26.usBeta). SystemExit for --main with --root or
    another variant than config.US['theta'], for a variant without a US or French θ, a missing or other-θ0 sweep, or
    results/shocks/ as the output without --main. """
    if main and root is not None:
        raise SystemExit('--main reads results/calibration/ and writes results/shocks/: not with --root')
    try:
        θ0, θ0FR = oecd26.thetaFromFile('USA', theta), oecd26.countryθ0('FR', theta)
        root = oecd26.thetaRoot(theta, main) if root is None else os.path.join(REPO, root)
    except (ValueError, FileNotFoundError) as e:          # NoTheta is a ValueError
        raise SystemExit('--theta {}: {}'.format(theta, e))
    calibDir = os.path.join(root, 'calibration')
    out = os.path.join(root, 'shocks', 'US_shocks' + ('CommonX' if commonX else '') + '.csv')
    if not main and os.path.normcase(os.path.dirname(os.path.abspath(out))) == os.path.normcase(OUTDIR):
        raise SystemExit('--theta {}: results/shocks/ holds the main variant\'s runs (--main) only; use another '
                         '--root'.format(theta))
    for ρ in rhos:
        try:
            oecd26.usBeta(os.path.join(calibDir, 'US_rhoGrid' + ('CommonX' if commonX else '') + '.csv'), ρ,
                          thetaSource = theta)
        except (FileNotFoundError, ValueError, KeyError) as e:
            raise SystemExit('--theta {}: {}'.format(theta, e))
    return θ0, calibDir, out, θ0FR


def main():
    p = argparse.ArgumentParser(description = 'US counterfactuals: theta, ageing, French characteristics.')
    p.add_argument('--rho', type = float, nargs = '*', default = [0.5, 1.0, 2.0])
    p.add_argument('--family', default = 'all', choices = ('all', 'theta', 'ageing', 'french'))
    p.add_argument('--commonX', action = 'store_true')
    p.add_argument('--n', type = int, default = 101)
    p.add_argument('--ns', type = int, default = 150)
    p.add_argument('--smoothKnots', type = int, default = 4)
    p.add_argument('--interpKind', default = 'linear')
    p.add_argument('--noFrance', action = 'store_true',
                   help = "skip France's own calibrated path, the reference row the French-characteristics "
                          'counterfactuals are read against. See franceReference.')
    p.add_argument('--theta', default = None,
                   help = "a variant of data/oecd26_theta.csv: the US at that θ (db['θ0']); see THE θ VARIANT")
    p.add_argument('--root', default = None, help = "the variant's directory (default results/oecd26_theta_<variant>)")
    p.add_argument('--coarse', action = 'store_true', help = 'the tax grid at 51 nodes, under --theta only')
    p.add_argument('--main', action = 'store_true',
                   help = "with --theta config.US['theta']: read results/calibration/, write results/shocks/")
    p.add_argument('--out', default = None)
    a = p.parse_args()
    a.out, a.root = [None if x is None else os.path.abspath(os.path.join(CWD, x)) for x in (a.out, a.root)]

    names = sum((list(FAMILIES[f]) for f in (FAMILIES if a.family == 'all' else [a.family])), [])
    famOf = {n: f for f, ns in FAMILIES.items() for n in ns}
    if a.theta in (None, oecd26.INDEX):
        if a.root is not None or a.coarse or a.main:
            raise SystemExit('--root, --coarse and --main need --theta VARIANT')
        θ0, θ0FR, calibDir, usCsv, tagθ = None, None, None, None, {}
        out = a.out or os.path.join(OUTDIR, 'US_shocks' + ('CommonX' if a.commonX else '') + '.csv')
        a.theta = oecd26.indexTheta(a.theta, [out], 'runShocksUS.py')
        os.makedirs(OUTDIR, exist_ok = True)
    else:
        θ0, calibDir, outθ, θ0FR = thetaRun(a.theta, a.root, a.commonX, a.rho, a.main)
        usCsv = os.path.join(calibDir, 'US_rhoGrid' + ('CommonX' if a.commonX else '') + '.csv')
        tagθ = {'thetaSource': a.theta, 'θ0': θ0, 'grid': 'coarse' if a.coarse else 'fine'}
        out = a.out or outθ
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok = True)
        if a.coarse:
            a.n = 51
        print('θ variant {!r}: θ0 = {!r} (France {!r}), US sweep {}'.format(a.theta, θ0, θ0FR, usCsv))

    workweek = float(testmod.workweek)
    rows = []
    print('host: US')
    print('scenarios: {}'.format(', '.join(names)))
    print('workweek (data, {}): {:.5f}'.format(testmod.dfc['Calibration year'], workweek))

    for ρ in a.rho:
        row = hostRow(ρ, a.commonX, 'US', calibDir)
        preferences = 'LOG' if ρ == 1 else 'CRRA'
        wp = {'smoothKnots': a.smoothKnots or None, 'interpKind': a.interpKind}
        gs = ({'n': a.n, 'ns': a.ns} if preferences == 'CRRA' else {'n': a.n}) | wp
        m = hostModel(ρ, a.commonX, row, gs, θ0 = θ0)

        t0 = m.db['t'][m.db['t0']]
        tic = time.time()
        base = m.solvePEE_LOG() if preferences == 'LOG' else m.solvePEE_CRRA()
        hbarRef = float(m.B.avgHours(base['report']['h'].xs(t0), t0))
        b = sh.readout(m, base['τ'], base['report'], workweek, hbarRef, pos = m.db['t0'])
        print('\nrho={:<4} {:<4} baseline: tau={:.4f}  sr={:.4f}  workweek={:.2f}  theta={:.4f}'
              '   ({:.0f}s to solve)'.format(ρ, preferences, b['τ'], b['sr'], b['workweek'],
                                             float(m.db['θ'].xs(t0)), time.time()-tic))
        # The calibration is a fixed point of these numbers, so they are also a check on the reinstall:
        # tau must equal the target and the workweek must be the data value by construction of hbarRef.
        # Against the sweep's OWN recorded tau, not the workbook target: on a vector-X ModelFR under CRRA the
        # recorded tau carries the hours-rescaling drift (modelFR.hoursDriftTol, <= 1e-3 from the target).
        if abs(b['τ'] - float(row['τ'])) > 1e-6:
            raise RuntimeError('US baseline at ρ={} gives τ={:.6f} but the sweep recorded {:.6f}: the reinstalled '
                               'calibration does not reproduce the sweep'.format(ρ, b['τ'], float(row['τ'])))
        rows.append({'ρ': ρ, 'preferences': preferences, 'family': 'baseline', 'scenario': 'Baseline',
                     'effect': 'baseline', **b} | multiplicityColumns(base.get('multiplicity')) | tagθ)

        gsFR = gs
        data = (frenchData(m, ρ, preferences, gsFR, a.commonX, usCsv = usCsv, θ0FR = θ0FR)
                if any(famOf[n] == 'french' for n in names) else {})
        if data:
            print('  French: Xbar_FR/Xbar_host={:.4f} ({:.2f} vs {:.2f}), eta level c={:.4f}, '
                  'max|gamma_FR-gamma_host|={:.1e}'.format(data['xbarRatio'], data['XbarFR'],
                                                        data['XbarUS'], data['ηScale'], data['γgap']))

        for name in names:
            tic = time.time()
            r = sh.runOne(m, base, name, data, preferences, workweek, hbarRef)
            rows += [e | tagθ for e in effectRows(r, {'ρ': ρ, 'preferences': preferences, 'family': famOf[name],
                                                      'scenario': r['label']})]
            print('  {:<22} full: tau={:.4f} sr={:.4f} ww={:.2f}   |   EE-only: tau={:.4f} sr={:.4f} '
                  'ww={:.2f}   {:.0f}s'.format(
                      name, r['full']['τ'], r['full']['sr'], r['full']['workweek'],
                      r['ee']['τ'], r['ee']['sr'], r['ee']['workweek'], time.time()-tic))
            writeRows(out, rows, tagθ)

        if not a.noFrance:
            tic = time.time()
            try:
                f = franceReference(ρ, preferences, gsFR, a.commonX, workweek, hbarRef, usCsv = usCsv,
                                    θ0FR = θ0FR)
                rows.append({'ρ': ρ, 'preferences': preferences, 'family': 'french',
                             'scenario': 'France (own calibration)', 'effect': 'full', **f} | tagθ)
                print('  {:<22} full: tau={:.4f} sr={:.4f} ww={:.2f}   |   theta={:.4f} omega={:.4f}'
                      '   {:.0f}s'.format('France', f['τ'], f['sr'], f['workweek'], f['θ'], f['ω'],
                                          time.time()-tic))
            except Exception as e:
                runESC.noteFailure(f'ρ={ρ} France reference', e)
            writeRows(out, rows, tagθ)

    writeRows(out, rows, tagθ)
    print('\n-> {}'.format(os.path.relpath(out, REPO)))
    if runESC.FAILURES:
        print('{} FAILED:'.format(len(runESC.FAILURES)) + ''.join('\n  ' + x for x in runESC.FAILURES))
    return 1 if runESC.FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())

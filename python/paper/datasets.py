r""" Loaders for everything stage (iii) reads. The only module that knows the results/ file layout and
column names, so a schema change is a one-file edit rather than a hunt through the builders.

Nothing here imports the model. A loader whose input is missing raises MissingInput, which build.py
turns into a skipped output and a line in the report -- an unrun experiment must never look like a
built table.
"""
import os, json, functools
import numpy as np, pandas as pd

import config as C


class MissingInput(Exception):
    """ A required results/ file does not exist. Carries the paths so the report can name them. """
    def __init__(self, paths):
        self.paths = [paths] if isinstance(paths, str) else list(paths)
        super().__init__('missing input(s): '
                         + ', '.join(os.path.relpath(p, C.REPO) for p in self.paths))


class StaleInput(ValueError):
    """ A results/ file exists but its rows were computed at other inputs than the current ones (an anchor of an
    earlier U.S. calibration or grid, a layer row of another anchor; crossCuttingFindings #13). A ValueError, so a
    caller that treats it as fatal need not know it. """


def _need(path):
    if not os.path.exists(path):
        raise MissingInput(path)
    return path


def _argCX(commonX):
    return C.ARG['commonX'] if commonX is None else bool(commonX)


def calibrationSummary(commonX = None):
    """ The one record per variant runCalibration.py writes, at the requested variant (the headline,
    config.ARG['commonX'], by default). Vector entries are JSON in the csv and come back as lists of
    float. A summary written before the variant column existed is a vector-X summary. """
    cx = _argCX(commonX)
    df = pd.read_csv(_need(os.path.join(C.PAPERDIR, 'calibrationSummary.csv')))
    if 'commonX' not in df.columns:
        df['commonX'] = False
    df = df[df['commonX'].fillna(False).astype(bool) == cx]
    if df.empty:
        raise MissingInput(os.path.join(C.PAPERDIR, 'calibrationSummary.csv')
                           + ' (no {} row)'.format('common-X' if cx else 'vector-X'))
    rec = df.iloc[-1].to_dict()
    for k in ('ν', 'ηi', 'Xi', 'γi', 'zxi', 'zxPredicted'):
        if isinstance(rec.get(k), str):
            rec[k] = json.loads(rec[k])
    # The calibrated tax hits the workbook's derived target (config.argTaxTarget, model.taxTarget) to the
    # calibration's tolerance; a summary at another target is a calibration of another economy.
    if not abs(float(rec['τ']) - C.calendar()['τ0']) < 1e-6:
        raise MissingInput(os.path.join(C.PAPERDIR, 'calibrationSummary.csv')
                           + ' (tau {:.7f} is not the workbook target {:.7f}: a stale calibration; rerun '
                             'stage (i))'.format(float(rec['τ']), C.calendar()['τ0']))
    return rec


def rhoGrid(commonX = None):
    """ The calibration sweep of one variant, one row per rho (a repeated rho is an error, #8). 'sr' is savings
    over GDP (the sweep's srGDP, GDP = formal plus informal output, model.taxTarget's concept); the sweep's own
    'sr', savings over formal output, is returned as 'srFormal'. A sweep without srGDP predates that concept and is
    MissingInput. """
    path = _need(C.argSweepCsv(_argCX(commonX)))
    df = pd.read_csv(path)
    if 'srGDP' not in df.columns:
        raise MissingInput(path + ' (no srGDP column: the sweep predates the GDP concept; rerun stage (i))')
    df = df.rename(columns = {'sr': 'srFormal'}).assign(sr = lambda d: d['srGDP'])
    return _noDuplicates(df, ['ρ'], path).sort_values('ρ').reset_index(drop = True)


def _shockCsv(ρ, scenario, rule, commonX):
    return os.path.join(C.SHOCKDIR, C.argShockTemplate(scenario, rule, _argCX(commonX)).format(ρ = ρ))


def shockPath(ρ, scenario = 'reform', rule = None, commonX = None):
    """ One rho's full response path, indexed on the model's own t with t0 first. """
    return pd.read_csv(_need(_shockCsv(ρ, scenario, rule, commonX)), index_col = 0)


def shockAtPeriod(period = 0, scenario = 'reform', rule = None, ρGrid = None, commonX = None):
    """ One row per rho, at `period` periods after t0 (0 = the impact period, the paper's year 2010).

    The row is taken POSITIONALLY (iloc), not by label: each csv is already indexed on the model's own
    t starting at t0, so period is an offset into that. Every rho in the grid must be present -- a
    partial sweep is raised on rather than silently plotted as a shorter curve. """
    ρGrid = C.ARG['ρGrid'] if ρGrid is None else ρGrid
    paths = {ρ: _shockCsv(ρ, scenario, rule, commonX) for ρ in ρGrid}
    lack = [p for p in paths.values() if not os.path.exists(p)]
    if lack:
        raise MissingInput(lack)
    rows = []
    for ρ, p in paths.items():
        df = pd.read_csv(p, index_col = 0)
        if period >= len(df):
            raise MissingInput(p + ' (has {} periods, need offset {})'.format(len(df), period))
        r = df.iloc[period].to_dict()
        r['ρ'] = ρ
        rows.append(r)
    return pd.DataFrame(rows).sort_values('ρ').reset_index(drop = True)


def epsThetaGrid(ρ = None, commonX = None):
    """ The Cartesian (eps, theta) grid at the calibration year, one row per pair.

    Two shape checks, because the producer is resumable on its own csv and so has two distinct ways of
    handing back a file that is the wrong ANSWER rather than a missing one.

    COMPLETE RECTANGLE. figures.argLogFourInOne pivots this into an eps x theta matrix and fills
    between adjacent theta columns; a missing pair becomes a NaN in the matrix and fill_between drops
    that span silently, so a partially-resumed sweep would render as a figure with holes.

    ONE statusQuo ROW, AT THE CURRENT CALIBRATION. sweepEpsThetaGrid.py resumes on the (eps, theta)
    key and pins the calibrated eps into its grid, so after a recalibration a resumed run ADDS the new
    calibrated column and keeps every row solved at the old parameters. The rectangle check cannot see
    that -- staleness adds a column rather than leaving a hole -- but the csv then carries two rows
    flagged statusQuo, and a wholly stale csv carries one that does not match calibrationSummary. Both
    are caught here.

    'sr' is savings over GDP (the csv's srGDP), the csv's own 'sr' (over formal output) 'srFormal', as
    rhoGrid. """
    ρ = C.ARG['ρBaseline'] if ρ is None else ρ
    path = _need(C.argEpsThetaCsv(ρ, _argCX(commonX)))
    df = pd.read_csv(path)
    if 'srGDP' not in df.columns:
        raise MissingInput(path + ' (no srGDP column: the sweep predates the GDP concept; rerun it with --force)')
    df = df.rename(columns = {'sr': 'srFormal'}).assign(sr = lambda d: d['srGDP'])
    nε, nθ = df['eps'].nunique(), df['theta'].nunique()
    if len(df) != nε * nθ:
        raise MissingInput('{} (incomplete grid: {} rows for {} eps x {} theta)'.format(
            path, len(df), nε, nθ))
    sq = df[df['statusQuo'] == True]
    if len(sq) != 1:
        raise MissingInput('{} ({} statusQuo rows, expected 1 -- the sweep resumed across a '
                           'recalibration and the csv mixes parameter sets; re-run it with --force)'
                           .format(path, len(sq)))
    cal = calibrationSummary(commonX)
    if np.isclose(cal['ρ'], ρ) and not (np.isclose(float(sq['eps'].iloc[0]), cal['ε'])
                                        and np.isclose(float(sq['theta'].iloc[0]), cal['θ'])):
        raise MissingInput('{} (status quo at (eps, theta) = ({:.6f}, {:.6f}), but the calibration is '
                           'at ({:.6f}, {:.6f}) -- the whole sweep predates it; re-run it with --force)'
                           .format(path, float(sq['eps'].iloc[0]), float(sq['theta'].iloc[0]),
                                   cal['ε'], cal['θ']))
    return df


# ---------------------------------------------------------------------------------------------------
# Derived quantities the csvs do not carry
# ---------------------------------------------------------------------------------------------------
def seedSavings(ρ, rule = None, tol = 1e-6, commonX = None):
    """ s_{t0-1}, the savings level entering the reform year. Common to every scenario, because the
    reform is unanticipated: the seed state is the baseline's in all of them.

    No shock csv carries it -- it is a state, not a row -- so it is recovered two ways and they check
    each other:

      * from the eeOnly csv's `s__base`, which reports the seed directly, when that file exists;
      * by inverting eq (calibration) at t0, where the baseline savings rate over formal output is
        reported by the calibration: sr = s/((s_/nu)^alpha h^(1-alpha)) solved for s_. The achieved
        rate is read off the sweep csv (rhoGrid's srFormal), so this stays exact at the solver's own
        precision.

    The inversion is what keeps every figure buildable before the eeOnly experiment has been run; the
    agreement check is what keeps it honest once it has. """
    ref  = shockPath(ρ, 'reform', rule, commonX)
    grid = rhoGrid(commonX)
    row  = grid.loc[np.isclose(grid['ρ'], ρ)]
    if row.empty:
        raise MissingInput('rho={} absent from {}'.format(ρ, os.path.basename(C.argSweepCsv(_argCX(commonX)))))
    cal  = calibrationSummary(commonX)
    α, ν = cal['α'], np.asarray(cal['ν'], dtype = float)[C.calendar()['t0']]
    sr, s, h = float(row['srFormal'].iloc[-1]), float(ref['s_base'].iloc[0]), float(ref['h_base'].iloc[0])
    inverted = ν * (s / (sr * h**(1-α)))**(1/α)

    path = _shockCsv(ρ, 'ee', rule, commonX)
    if os.path.exists(path):
        reported = float(pd.read_csv(path, index_col = 0)['s__base'].iloc[0])
        if abs(reported/inverted - 1) > tol:
            raise ValueError('seed savings disagree at rho={}: eeOnly csv {!r} vs eq(calibration) '
                             'inversion {!r} ({:.2e} relative). One of the two runs is stale.'
                             .format(ρ, reported, inverted, abs(reported/inverted - 1)))
        return reported
    return inverted


def baselineHours(ρ, rule = None, commonX = None):
    """ The calibrated baseline's aggregate hours at t0, for this rho: the reference point that
    config.workweekHours normalises against. Per rho, because each rho is separately calibrated and h is
    not one of the four targets -- so its level differs across the grid even though tau and the savings
    rate do not. Normalising each rho against its own baseline is what makes the pre-reform workweek
    42.54 everywhere and the reform rows comparable. """
    return float(shockPath(ρ, 'reform', rule, commonX)['h_base'].iloc[0])


def _gdpPath(ρ, which, rule, commonX, tol = 1e-8):
    """ (path frame, Y, Y0) of one scenario: formal and informal output per formal young household, the
    csv's own Y_<which>/Y0_<which> (EE_report), so GDP = Y + Y0 (model.taxTarget's concept). The formal
    output is checked against the technology (s_{t-1}/nu_t)^alpha h_t^(1-alpha) at the path's own lag, the
    seed state at t0 (seedSavings): the check ties the csv's columns to the state the reform starts from.
    A csv without the columns predates the GDP concept and is MissingInput. """
    path = _need(_shockCsv(ρ, 'ee' if which == 'ee' else 'reform', rule, commonX))
    df   = pd.read_csv(path, index_col = 0)
    if 'Y_' + which not in df.columns or 'Y0_' + which not in df.columns:
        raise MissingInput(path + ' (no Y/Y0 columns: the run predates the GDP concept; rerun stage (ii))')
    cal  = calibrationSummary(commonX)
    α, ν = cal['α'], np.asarray(cal['ν'], dtype = float)
    t0i  = C.calendar()['t0']
    s, h = df['s_' + which].values.astype(float), df['h_' + which].values.astype(float)
    Y, Y0 = df['Y_' + which].values.astype(float), df['Y0_' + which].values.astype(float)
    s_ = np.concatenate([[seedSavings(ρ, rule, commonX = commonX)], s[:-1]])
    dev = np.max(np.abs((s_/ν[t0i:t0i+len(s)])**α * h**(1-α)/Y - 1))
    if not dev <= tol:
        raise ValueError('{}: formal output Y_{} disagrees with the technology at the seed state by {:.2e} '
                         'relative; the csv and the calibration are from different runs.'.format(path, which, dev))
    return df, Y, Y0


def savingsRatePath(ρ, which = 'reform', rule = None, commonX = None):
    """ Savings over GDP, s_t/(Y_t + Y0_t), over the whole reported path (_gdpPath). `which` is a column
    suffix: 'base', 'reform' (shockUniversal's csv) or 'ee' (shockEEOnly's). The terminal period is NaN by
    construction (s_T = 0), and is left NaN rather than zero-filled: a zero there would plot as a real
    collapse. """
    df, Y, Y0 = _gdpPath(ρ, which, rule, commonX)
    s = df['s_' + which].values.astype(float)
    with np.errstate(divide = 'ignore', invalid = 'ignore'):
        sr = s/(Y + Y0)
    sr[s <= 0] = np.nan                       # terminal period: s_T = 0 is degenerate, not a datum
    return pd.Series(sr, index = df.index)


def reformSavingsRate(ρ, period = 0, rule = None, commonX = None):
    """ Savings over GDP on the full-effect path at `period` periods after t0. """
    return float(savingsRatePath(ρ, 'reform', rule, commonX).iloc[period])


def spendingPath(ρ, which = 'reform', rule = None, commonX = None):
    """ Pension spending over GDP, tau_t (1-alpha) Y_t/(Y_t + Y0_t), over the whole reported path: the tax
    is levied on formal labour income (1-alpha) Y_t, GDP includes informal output (model.taxTarget). No builder
    calls it: it is the entry point REPLICATION.md names for the prose numbers of section 5. """
    df, Y, Y0 = _gdpPath(ρ, which, rule, commonX)
    α = calibrationSummary(commonX)['α']
    return pd.Series(df['τ_' + which].values.astype(float)*(1-α)*Y/(Y + Y0), index = df.index)


def formalBenefitChange(ρ, which = 'reform', period = 0, rule = None, commonX = None, tol = 1e-10):
    """ The change in formal benefits against the baseline, b_t^i(which)/b_t^i(base) - 1, at `period`
    periods after t0. The same for every formal type (eq:governmentBudget: b_t^i is a type-specific
    parameter bracket times h_{t-1} bbar_t), so it is also the change in their average biAvg: the
    shockUniversal csv's d_biAvg, checked against its reference type's d_bi; shockEEOnly's csv carries d_bi
    only (its EE_report has no lagged state at t0, which b^i does not depend on). """
    path = _need(_shockCsv(ρ, 'ee' if which == 'ee' else 'reform', rule, commonX))
    r = pd.read_csv(path, index_col = 0).iloc[period]
    if 'd_biAvg' in r.index:
        if abs(r['d_biAvg'] - r['d_bi']) > tol:
            raise ValueError('{}: d_biAvg {} and d_bi {} differ at period {}; formal benefits did not move '
                             'uniformly'.format(path, r['d_biAvg'], r['d_bi'], period))
        return float(r['d_biAvg'])
    if which == 'reform':
        raise MissingInput(path + ' (no d_biAvg column: the run predates the formal-benefit readout)')
    return float(r['d_bi'])


# ---------------------------------------------------------------------------------------------------
# US / France / UK
# ---------------------------------------------------------------------------------------------------

def _readUS(path, **kw):
    """ A U.S.-arm csv (results/{calibration,esc,shocks,oecd26}, the U.S. stationary checks, the summary) at the θ
    variant config.US['theta'] only (config.atTheta; a csv without thetaSource reads as 'index') and, for a csv of a
    driver with a cost (config.spendingKeyed), at the spending target config.US['esc']['spendingTarget'] only
    (config.atSpending; a csv without spendingTarget reads as 'gross'): MissingInput when the file is absent or holds no
    row at them. Every reader below goes through this, so a row of another θ or target is never read
    (crossCuttingFindings #13). """
    df = pd.read_csv(_need(path), **kw)
    out = C.atTheta(df)
    if out.empty:
        raise MissingInput('{} (no row at θ variant {})'.format(path, C.US['theta']))
    if C.spendingKeyed(path):
        out = C.atSpending(out)
        if out.empty:
            raise MissingInput('{} (no row at spending target {!r})'.format(path, C.US['esc']['spendingTarget']))
    return out.copy()


def usCalibrationSummary(commonX = None):
    """ runCalibrationUS.py's one-row-per-country record, keyed by country. Vector entries are JSON in
    the csv and come back as lists of float.

    commonX = None means the variant the paper leads with, config.US['commonX']. Every US loader here
    takes it that way, so a builder that says nothing gets the headline calibration and a builder that
    wants the robustness twin has to ask for it by name. """
    commonX = C.US['commonX'] if commonX is None else commonX
    df = _readUS(os.path.join(C.PAPERDIR, 'usCalibrationSummary.csv'))
    df = df[df['commonX'].astype(bool) == bool(commonX)]
    if df.empty:
        raise MissingInput(os.path.join(C.PAPERDIR, 'usCalibrationSummary.csv'))
    out = {}
    for rec in df.to_dict('records'):
        for k in ('ηi', 'Xi', 'γi', 'μi', 'ν'):
            if isinstance(rec.get(k), str):
                rec[k] = json.loads(rec[k])
        out[rec['country']] = rec
    return out


@functools.lru_cache(maxsize = None)
def usAlpha(commonX = None):
    """ The capital share of the US/FR/UK model, read from usCalibrationSummary.csv (the US row of the requested
    variant) and checked against the workbook's `Capital income share` cell, which is where the summary got it.
    MissingInput without the summary or its U.S. row, ValueError if the two disagree: one number, from data, never
    typed here and never taken from the workbook alone. Cached (build.py --map clears it, via clearCaches). """
    fromWb = C.usCalendar('US')['α']
    path = os.path.join(C.PAPERDIR, 'usCalibrationSummary.csv')
    rec = usCalibrationSummary(commonX).get('US')
    if rec is None or 'α' not in rec or rec['α'] != rec['α']:
        raise MissingInput(path + ' (no α for the US)')
    fromCsv = float(rec['α'])
    if not np.isclose(fromCsv, fromWb):
        raise ValueError('capital share disagrees: usCalibrationSummary.csv {!r} vs workbook {!r}'
                         .format(fromCsv, fromWb))
    return fromCsv


def escSavingsOverY(sr, commonX = None):
    """ ESC csv savings rate s/(w h) -> the paper's s/Y.

    The ESC csvs (escExperiments/escShocks*/escCalibration*/escCountry) carry `sr_*` as savings over
    gross labour income, s/(w h); the US shock csv carries both that (`sr`) and s/Y (`srOverY`), and
    the paper reports s/Y everywhere. With Cobb--Douglas production and a competitive labour market
    w h = (1-alpha) Y identically, so s/Y = (1-alpha) * s/(w h) is EXACT in every period and scenario,
    not an approximation -- which is why converting here rather than re-running the ESC leg loses
    nothing. `sr` may be a scalar or an array. """
    return np.asarray(sr, dtype = float) * (1 - usAlpha(commonX))


US_SR = 'srOverY'    # the US shock csv column the tables read: s/Y. `sr` is s/(w h) and is NOT read.


def usShocks(commonX = None, host = 'US'):
    """ python/US/runShocksUS.py's long csv at the θ variant: one row per (rho, family, scenario, effect). `host` is
    the economy the characteristics were imposed on, the U.S. (config.usShockCsv).

    `effect` is 'baseline' | 'full' | 'ee'. tau, sr and srOverY are fractions -- `sr` is s/(w h),
    `srOverY` is s/Y, the paper's savings rate (US_SR); `workweek` is already in hours,
    normalised against that rho's own baseline inside the experiment script -- stage (iii) must NOT
    re-derive it from hbar, which has no scale under vector X (see python/US/shocks.py).

    commonX = None is the headline variant (config.US['commonX']). """
    commonX = C.US['commonX'] if commonX is None else commonX
    path = C.usShockCsv(host, commonX)
    df = _readUS(path)
    df.attrs['source'] = 'results/shocks/' + os.path.basename(path)
    return df


def usShockRow(df, ρ, scenario, effect):
    """ One (rho, scenario, effect) row, raising rather than returning an empty frame -- a silently
    missing scenario would render as a blank table cell instead of a skipped output. """
    hit = df[(np.isclose(df['ρ'], ρ)) & (df['scenario'] == scenario) & (df['effect'] == effect)]
    if hit.empty:
        raise MissingInput('{} (ρ={}, {}) in {}'
                           .format(scenario, ρ, effect, df.attrs.get('source', 'the US shock csv')))
    if len(hit) > 1:
        raise ValueError('{}: {} rows for ρ={}, {}, {}'.format(df.attrs.get('source', 'the US shock csv'), len(hit),
                                                              ρ, scenario, effect))
    return hit.iloc[0]


def usBaseline(df, ρ):
    return usShockRow(df, ρ, 'Baseline', 'baseline')


def usShocksXi():
    """ python/US/runShocksXi.py's csv: usShocks' columns at rho = 1, common X, plus the Frisch elasticity
    'ξ', the US and France recalibrated at each. """
    path = os.path.join(C.SHOCKDIR, 'US_shocksXiCommonX.csv')
    df = _readUS(path)
    df.attrs['source'] = 'results/shocks/US_shocksXiCommonX.csv'
    return _noDuplicates(df, ['ξ', 'ρ', 'scenario', 'effect'], path)


def usSweep(country, commonX = None):
    """ One country's calibration sweep at the θ variant, one row per rho (a repeated rho is an error, #8). """
    commonX = C.US['commonX'] if commonX is None else commonX
    path = C.usSweepCsv(country, commonX)
    return _noDuplicates(_readUS(path), ['ρ'], path).sort_values('ρ').reset_index(drop = True)


# ---------------------------------------------------------------------------------------------------
# Endogenous system characteristics (results/esc/)
# ---------------------------------------------------------------------------------------------------
def escCornerSlope(commonX = None, spec = None):
    """ python/US/runESC.stageCorner's csv: per economy (FR, US, UK), lambda (0, the U.S. value, 100) and
    electorate (`vote`, the db['t'] position; `date`), the slope of the leaded design objective at
    theta_{t+1} = 1 and the anticipated tax there (`τ1Corner`); at the requested variant, spec and phi.
    Each row is recalibrated at its lambda (omega moves wherever the observed design is below the corner). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    path = os.path.join(C.ESCDIR, 'escCornerSlope.csv')
    df = _readUS(path)
    keep = [i for i, rec in enumerate(df.to_dict('records'))
            if _escVariant(rec, commonX) and rec['spec'] == spec
            and np.isclose(float(rec['phi']), C.US['esc']['phi'])]
    return _noDuplicates(df.iloc[keep].copy(), ['country', 'λ', 'vote'], path)


def escCountry(commonX = None, spec = None, k = 2):
    """ python/US/runESC.stageCountry's csv: France and the UK under the ESC cost, at the requested variant,
    spec (config.US['esc']['spec'] by default) and phi; rows of the UK regrouped at US income percentiles
    ('UKUS'), on disk from earlier runs, are returned with the rest and no builder prints them. Two readings
    per economy, `wedgeFrom` 'US' (the US-calibrated parameter imposed; `choice` is what that electorate then
    picks, `τ` its observed 2020 tax rate) and 'own' (the parameter calibrated to re-elect the economy's own
    design; only converged rows are kept, so an economy whose observed design no finite value reaches has no
    'own' row). `Vtilde` is present on rows written since the
    'size' spec. The stage runs under LOG only, and the csv carries no rho column, so rho = 1 is stamped
    on here rather than read -- a CRRA country stage would have to write the column itself. The spec of
    config.US['esc']['costVariants'] ('sizePiR') is read from its tagged csv (escCountryCosts.csv), where
    runESC.py --tag writes it; every other spec from escCountry.csv. Each economy at its θ0 under the θ variant; one
    row per (country, wedgeFrom), a repeated pair being an error. `k`: the cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    cv = C.US['esc']['costVariants']
    name = 'escCountry' + (cv['tag'] if spec == cv['spec'] else '') + '.csv'
    df, _ = _readESC(name, k)
    path = os.path.join(escDir(k), name)
    keep = [i for i, rec in enumerate(df.to_dict('records'))
            if _escVariant(rec, commonX) and rec['spec'] == spec
            and np.isclose(float(rec['phi']), C.US['esc']['phi'])
            and (rec['wedgeFrom'] != 'own' or (rec['converged'] == rec['converged']
                                               and bool(rec['converged'])))]
    df = _noDuplicates(df.iloc[keep].copy(), ['country', 'wedgeFrom'], path)
    if 'ρ' not in df.columns:
        df['ρ'] = 1.0
    return df


def _escVariant(rec, commonX):
    """ Is this ESC row from the calibration variant we are reading?

    The csvs carry a `commonX` column and the drivers key their merge on it, so both variants can sit in
    one file. A row from before that column existed is a vector-X row -- that is the variant everything
    ran under then -- which is what the missing-column fallback means, not "applies to either".

    A MISSING value reads as False for the same reason, and it has to be handled explicitly: pandas
    returns NaN for a blank cell, and `bool(nan)` is True, so the obvious `bool(rec.get(...))` would file
    every untagged row under common X -- the one variant it certainly is not. """
    v = rec.get('commonX', False)
    v = False if v is None or v != v else bool(v)      # v != v catches NaN
    return v == bool(commonX)


def escMethod():
    """ The CRRA ESC solver the paper prints: 'exact' (LeadedCRRA2D) or 'path' (the path iteration),
    config.US['esc']['exact']. Every CRRA row in results/esc carries a `method` column; a row written
    before the column existed is a path-iteration row, which is what a missing value means. """
    return 'exact' if C.US['esc'].get('exact', True) else 'path'


def _escMethodOK(rec):
    v = rec.get('method', None)
    v = 'path' if v is None or v != v else str(v)          # v != v catches NaN
    return v == escMethod()


def escDir(k = 2):
    """ The directory of the U.S. design leg under the cost with tau^k: results/esc/ at k = 2, tau1Dir('esc') at
    config.US['tau1']'s power. """
    if k == 2:
        return C.ESCDIR
    if k != C.US['tau1']['tauPower']:
        raise ValueError("escDir: k = {} (2 or config.US['tau1']['tauPower'])".format(k))
    return tau1Dir('esc')


def _readESC(name, k = 2):
    """ _readUS of escDir(k)/<name>; at k other than 2 every row must carry config.US['tau1']'s power and carrier
    (_tau1Rows). Returns (frame, its path relative to the repository, forward slashes). """
    path = os.path.join(escDir(k), name)
    df = _readUS(path)
    return (df if k == 2 else _tau1Rows(df, path)), os.path.relpath(path, C.REPO).replace(os.sep, '/')


def escCalibration(commonX = None, k = 2):
    """ The calibrated deadweight cost, {(rho, spec): record} with p (the cost parameter under either
    spec: lambda under 'size'), thetaStar, beta, omega, and -- on rows written since the 'size' spec --
    fStar, f0, Vtilde and τ0 (escWedge reads them).

    Two files because two solvers produced them: escCalibration.csv is the LOG case (rho = 1; also
    carries the no-wedge row under spec 'none', keyed here as (1.0, 'none')), escCalibrationCRRA.csv the
    CRRA rows with their own rho column, at the published method (escMethod) only. Only converged rows at
    config.US['esc']['phi'] AND at the requested calibration variant are returned -- an unconverged,
    wrong-variant or wrong-method calibration must surface as a missing key, not as a cell of a number
    that was never calibrated for this economy by the published solver. `k`: the cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    phi = C.US['esc']['phi']
    out = {}
    log, _ = _readESC('escCalibration.csv', k)
    for rec in log.to_dict('records'):
        if not _escVariant(rec, commonX):
            continue
        if rec['spec'] == 'none':
            out[(1.0, 'none')] = rec
        elif bool(rec['converged']) and np.isclose(float(rec['phi']), phi):
            out[(1.0, rec['spec'])] = rec
    crra, _ = _readESC('escCalibrationCRRA.csv', k)
    for rec in crra.to_dict('records'):
        if (_escVariant(rec, commonX) and _escMethodOK(rec) and bool(rec['converged'])
                and np.isclose(float(rec['phi']), phi)):
            out[(float(rec['ρ']), rec['spec'])] = rec
    return out


def escCalibrationLOG(spec, commonX = None, k = 2):
    """ The converged LOG (rho = 1) cost calibration of the U.S. under `spec`, one record (p, θStar, β, ω, Vtilde,
    ...; 'ΠR' under 'sizePiR'): from escCalibration.csv, or for config.US['esc']['costVariants']'s spec from its
    tagged csv (escCalibrationCosts.csv), where runESC.py --tag writes it. 'none' is the row without the cost. `k`: the
    cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    cv = C.US['esc']['costVariants']
    name = 'escCalibration' + (cv['tag'] if spec == cv['spec'] else '') + '.csv'
    df, _ = _readESC(name, k)
    path = os.path.join(escDir(k), name)
    keep = [rec for rec in df.to_dict('records')
            if _escVariant(rec, commonX) and rec['spec'] == spec
            and (spec == 'none' or (bool(rec['converged']) and np.isclose(float(rec['phi']), C.US['esc']['phi'])))]
    if len(keep) != 1:
        raise MissingInput('{} ({} converged {} rows at commonX={})'.format(path, len(keep), spec, commonX))
    return keep[0]


def _col(rec, k):
    """ A float column that may be absent, or blank (NaN) on a row written before it existed: None then. """
    v = rec.get(k, None)
    return None if v is None or v != v else float(v)


def escWedge(rec, spec):
    """ f(θ*), f(0) and Ṽ of one cost-calibration row (escCalibration{,CRRA}.csv), {'fStar', 'f0',
    'Vtilde'} as floats, None where the row does not carry the quantity.

    Rows written since the 'size' spec carry `fStar` = f(θ*, τ_t0), `f0` = f(0, τ_t0) and `Vtilde` under
    every spec and are read as they are. An older row has none of them; under 'scale' the two f values
    are then f(θ) = φ + (1-φ) θ^p from the row's own phi and p (so f(0) = φ) and Ṽ is None. Under any
    other spec the form needs Ṽ and τ, which such a row does not carry, so the row is MissingInput rather
    than a number computed from the wrong form. """
    out = {k: _col(rec, k) for k in ('fStar', 'f0', 'Vtilde')}
    if out['fStar'] is None or out['f0'] is None:
        if spec != 'scale':
            raise MissingInput('columns fStar, f0, Vtilde of the {!r} rows in results/esc/escCalibration{{,CRRA}}.csv'
                               .format(spec))
        φ, p, θ = float(rec['phi']), float(rec['p']), float(rec['θStar'])
        out['fStar'], out['f0'] = φ + (1-φ)*θ**p, φ
    return out


def escCalibrationHost(host = 'UK', commonX = None, spec = None, k = 2):
    """ A non-US host's OWN calibrated cost per rho, {rho: record} with p (lambda under 'size'): rho = 1
    from escCountry.csv (wedgeFrom 'own'), the other rho from escCalibrationCRRA<host>.csv at the
    published method. Converged rows at the variant, spec and phi only, as escCalibration. `k`: the cost's power of
    the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    phi = C.US['esc']['phi']
    out = {}
    cty = escCountry(commonX = commonX, spec = spec, k = k)
    own = cty[(cty['country'] == host) & (cty['wedgeFrom'] == 'own')]
    if not own.empty:
        out[1.0] = own.iloc[0].to_dict()
    path = os.path.join(escDir(k), 'escCalibrationCRRA' + host + '.csv')
    if os.path.exists(path):
        crra = C.atSpending(C.atTheta(pd.read_csv(path)))
        for rec in (crra if k == 2 else _tau1Rows(crra, path)).to_dict('records'):
            if (_escVariant(rec, commonX) and _escMethodOK(rec) and bool(rec['converged'])
                    and rec['spec'] == spec and np.isclose(float(rec['phi']), phi)):
                out[float(rec['ρ'])] = rec
    return out


def escExperiments(k = 2):
    """ collectESCexperiments.py's merged long csv: one row per (rho, spec, scenario, reading), with the
    design, tax, savings rate and workweek at t0 AND t0+1. `θpinned` True is the exogenous-theta reading
    (design held at the shocked model's own exogenous value), False the endogenous one.

    The t0 columns are what the tables read. Every scenario is a NEW EQUILIBRIUM PATH -- the shocked
    parameters hold over the whole horizon and the political choice binds from the first period -- so
    theta_{t0} is an equilibrium outcome and 2020 already carries the design response. The scenario 'France' is
    France's own calibrated path, not a shock on the US model, and carries only the pinned reading. `k`: the cost's
    power of the tax (escDir). """
    df, src = _readESC('escExperiments.csv', k)
    df.attrs['source'] = src
    return df


def escRow(df, ρ, spec, scenario, pinned, commonX = None):
    """ One (rho, spec, scenario, reading) row, raising rather than returning an empty frame -- same
    contract as usShockRow. The calibration variant is part of the key for the reason _escVariant gives.
    """
    commonX = C.US['commonX'] if commonX is None else commonX
    cx = (df['commonX'].fillna(False).astype(bool) if 'commonX' in df.columns
          else pd.Series(False, index = df.index))
    # The LOG rows (rho = 1) are exact by construction and pass under either setting; the CRRA rows must
    # carry the published method (a merged csv from before the column is all path-iteration rows).
    meth = (df['method'].fillna('path').astype(str) if 'method' in df.columns
            else pd.Series('path', index = df.index))
    isLOG = (df['preferences'] == 'LOG') if 'preferences' in df.columns else np.isclose(df['ρ'], 1.0)
    hit = df[np.isclose(df['ρ'], ρ) & (df['spec'] == spec) & (df['scenario'] == scenario)
             & (df['θpinned'].astype(bool) == bool(pinned)) & (cx == bool(commonX))
             & np.isclose(df['phi'], C.US['esc']['phi']) & (isLOG | (meth == escMethod()))]
    if hit.empty:
        raise MissingInput('{} (ρ={}, {}, pinned={}, commonX={}, method={}) in results/esc/escExperiments{{,UK}}.csv'
                           .format(scenario, ρ, spec, pinned, commonX, escMethod()))
    if len(hit) > 1:
        raise ValueError('{}: {} rows for ρ={}, {}, pinned={}, commonX={}, method={}'.format(
            df.attrs.get('source', 'the ESC csv'), len(hit), ρ, scenario, pinned, commonX, escMethod()))
    return hit.iloc[0]


# ---------------------------------------------------------------------------------------------------
# The supplemental appendix's exhibits (tablesOA.py, figuresOA.py): figure 1's data, the calibrations across
# rho, the endogenous design's baseline path, xi and timing checks, and the numerical checks.
# ---------------------------------------------------------------------------------------------------
OECDFIG1        = os.path.join(C.DATA, 'oecdFigure1.csv')           # stage (0): oecdFigure1.main
OECDFIG1SOURCES = os.path.join(C.DATA, 'oecdFigure1_sources.csv')
OECDCORR        = os.path.join(C.PAPERDIR, 'oecdCorrelations.csv')  # written by the OECDdata build


def oecdFigure1Data():
    """ Figure 1's cross-section, one row per country indexed on ISO3, in the footnote's order: every
    concept oecdFigure1 fetches, each with `<column>_year` where the year it is read at varies by
    country. """
    return pd.read_csv(_need(OECDFIG1), index_col = 'iso3')


def oecdFigure1Sources():
    """ One row per data column (or group of columns sharing a source) of oecdFigure1.csv. Every column is
    text, and a blank cell comes back as '' rather than NaN. """
    return pd.read_csv(_need(OECDFIG1SOURCES), dtype = str, keep_default_na = False)


def oecdCorrelations():
    """ Pairwise Pearson and Spearman correlations of figure 1's four variables per concept set (`variant`
    'named' | 'asEPS' | 'epsReadings' | 'figure'), pairwise complete, as the OECDdata build writes them; the role
    'inequality' is the Gini in the first three, config.INEQUALITY's measure in 'figure'. """
    return pd.read_csv(_need(OECDCORR))


# Page D1's section on inequality: the csvs the OECD_InequalityMeasures and OECD26_Layer1Inequality builds write into
# config.PAPERDIR (read at call time, so that a test's PAPERDIR holds them).
OECDINEQCORR    = 'oecdInequalityCorrelations.csv'
OECD26INEQSLOPE = 'oecd26Layer1InequalitySlopes.csv'


def paperCsv(name):
    """ results/paper/<name> under the current config.PAPERDIR. """
    return os.path.join(C.PAPERDIR, name)


def oecdInequalityCorrelations():
    """ The Gini against the model's measures of income inequality on figure 1's countries (tablesOA.
    inequalityCorrelations), as the OECD_InequalityMeasures build writes it. """
    return pd.read_csv(_need(paperCsv(OECDINEQCORR)))


def oecd26InequalitySlopes():
    """ Figure 7's spending and index on each measure of income inequality, data and model at the paper's variant
    (tablesOECD26.inequalitySlopes), as the OECD26_Layer1Inequality build writes it. """
    return pd.read_csv(_need(paperCsv(OECD26INEQSLOPE)))


# The cross-country test: python/US/runOECD26.py writes results/oecd26/, stage (0)'s oecd26Sample.py the sample.
OECD26DIR    = os.path.join(C.RESULTS, 'oecd26')
OECD26SAMPLE = os.path.join(C.DATA, 'oecd26_sample.csv')


def oecd26Sample():
    """ One row per country of figure 1, indexed on ISO3: `sample26` (the sample rule, applied in stage (0) and
    nowhere else) and the informality datum behind it. """
    return pd.read_csv(_need(OECD26SAMPLE), index_col = 'iso3')


def thetaVariant(variant = None):
    """ data/oecd26_theta.csv under `variant` (config.US['theta'] by default), indexed on ISO3: `theta`, `inSample`
    (bool), `spending`, and `index` = 1/(2 - θ), the design on figure 1.1's scale (NaN where θ is). `sample` is the
    variant's sample: inSample and a finite θ. """
    variant = C.US['theta'] if variant is None else variant
    d = pd.read_csv(_need(C.THETAFILE), dtype = {'iso3': str, 'variant': str}, float_precision = 'round_trip')
    d = d[d['variant'] == variant].set_index('iso3')
    if d.empty:
        raise MissingInput('{} (no variant {!r})'.format(C.THETAFILE, variant))
    d['theta'] = pd.to_numeric(d['theta'], errors = 'coerce')
    d['inSample'] = d['inSample'].astype(str).str.strip().str.lower().isin(('true', '1', '1.0'))
    d['index'] = 1/(2 - d['theta'])
    d['sample'] = d['inSample'] & d['theta'].notna()
    return d


def oecd26SampleOf(variant = None):
    """ The cross-section's countries under the θ variant, in figure 1.1's order: thetaVariant's `sample`. """
    t = thetaVariant(variant)
    return [c for c in oecdFigure1Data().index if c in t.index and bool(t.at[c, 'sample'])]


# The income groups the model is fed (python/US/oecd26.readData), the source of config.INEQUALITY's model measures.
OECD29DIST    = os.path.join(C.DATA, 'oecd29_distribution.csv')
OECD29SOURCES = os.path.join(C.DATA, 'oecd29_sources.csv')
GROUPTOL      = 1e-6                                   # the csv's eight digits: Σγy = 1 and Ṽ hold to ~2e-8


def incomeGroups(isos = None):
    r""" config.INCOMEGROUPS' rows of data/oecd29_distribution.csv, one row per country of `isos` (the rows' countries
    when None), indexed on ISO3 in that order: `y1`-`y3`, `gamma1`-`gamma3`, `Vtilde` = Σ_i γ_i (y_i - 1)²/y_i (Eq
    esc:budget's \tilde V), `yHighOverLow` = y_3/y_1 and `year`. Raises unless every country of `isos` has exactly the
    three groups, with finite y and γ, Σγ = Σγy = 1, one year, and the file's Vtilde = Σγ(y - 1)²/y. """
    d = pd.read_csv(_need(OECD29DIST))
    g = C.INCOMEGROUPS
    d = d[(d['role'] == g['role']) & (d['source'] == g['source']) & (d['convention'] == g['convention'])]
    isos = list(dict.fromkeys(d['iso3'])) if isos is None else list(isos)
    d = d[d['iso3'].isin(isos)]
    n = d.groupby('iso3')['group'].apply(lambda s: sorted(s) == [1, 2, 3]).reindex(isos, fill_value = False)
    if not n.all():
        raise ValueError('{}: not the three income groups of {} for {}'.format(
            os.path.relpath(OECD29DIST, C.REPO), g, ', '.join(n.index[~n])))
    y = d.pivot(index = 'iso3', columns = 'group', values = 'y').loc[isos].astype(float)
    γ = d.pivot(index = 'iso3', columns = 'group', values = 'gamma').loc[isos].astype(float)
    first = d.groupby('iso3').first().loc[isos]
    V = (γ*(y - 1)**2/y).sum(axis = 1)
    bad = ~(np.isfinite(y).all(axis = 1) & np.isfinite(γ).all(axis = 1)
            & ((γ.sum(axis = 1) - 1).abs() <= GROUPTOL) & (((γ*y).sum(axis = 1) - 1).abs() <= GROUPTOL)
            & ((first['Vtilde'].astype(float) - V).abs() <= GROUPTOL)
            & (d.groupby('iso3')['year'].nunique().reindex(isos) == 1))
    if bad.any():
        raise ValueError('{}: the groups of {} are not y, γ with Σγ = Σγy = 1 and Vtilde = Σγ(y - 1)²/y for {}'
                         .format(os.path.relpath(OECD29DIST, C.REPO), g, ', '.join(bad.index[bad])))
    out = pd.DataFrame({'y{}'.format(k): y[k] for k in (1, 2, 3)} | {'gamma{}'.format(k): γ[k] for k in (1, 2, 3)},
                       index = pd.Index(isos, name = 'iso3'))
    out['Vtilde'] = first['Vtilde'].astype(float)
    out['yHighOverLow'] = y[3]/y[1]
    out['year'] = first['year'].astype(int)
    return out


def oecd29Sources():
    """ data/oecd29_sources.csv (stage (0)'s oecd29Data.py): one row per data column or group of columns, as text. """
    return pd.read_csv(_need(OECD29SOURCES), dtype = str, keep_default_na = False)


def oecd26Dir(theta = None, root = None):
    """ Where the cross-section's rows of thetaSource `theta` (config.US['theta'] by default) are: results/oecd26/ for
    the main variant and for 'index' (the rows of figure 1.1's index, the paper's measure before the variant),
    results/oecd26_theta_<theta>/ for a data variant (python/US/runOECD26.py --theta without --main); `root`, a
    directory, overrides both. """
    if root is not None:
        return root
    theta = C.US['theta'] if theta is None else theta
    return OECD26DIR if theta in (C.US['theta'], 'index') else os.path.join(C.RESULTS, 'oecd26_theta_' + theta)


# A turnout whose U.S. inputs are another's stands on that one's anchor (python/US/runOECD26.ANCHORTURNOUT, read here
# without importing model code): Luxembourg's flat profile leaves the U.S. as it is.
ANCHORTURNOUT = {'csesFlat': 'cses'}


def oecd26OmegaRule(layer, omegaRule = None):
    """ Layer `layer`'s ω rule: `omegaRule`, else config.US['oecd26']['omegaRule'] for layer 1; layer 2 calibrates each
    country's ω and is always 'anchor' (a layer-1 rule asked of it is an error). """
    if layer == 2:
        if omegaRule not in (None, 'anchor'):
            raise ValueError("layer 2 calibrates each country's ω: omegaRule 'anchor' only, got {!r}".format(omegaRule))
        return 'anchor'
    return C.US['oecd26']['omegaRule'] if omegaRule is None else omegaRule


def oecd26DirRel(theta = None):
    """ oecd26Dir(theta) as the banners of the generated files name it: 'results/oecd26/' or
    'results/oecd26_theta_<V>/'. """
    return 'results/' + os.path.basename(oecd26Dir(theta)) + '/'


def oecd26Path(layer, theta = None, root = None):
    """ anchor.csv (layer 0) or layer<1|2>.csv of oecd26Dir(theta, root). """
    return os.path.join(oecd26Dir(theta, root), 'anchor.csv' if layer == 0 else 'layer{}.csv'.format(layer))


def _oecd26Variant(turnout = None, spec = None, betaRule = None):
    """ (turnout, spec, betaRule), each config.US['oecd26']'s where None. """
    v = C.US['oecd26']
    return (v['turnout'] if turnout is None else turnout, v['spec'] if spec is None else spec,
            v['betaRule'] if betaRule is None else betaRule)


def _oecd26Rows(layer, turnout = None, spec = None, betaRule = None, theta = None, root = None, omegaRule = None,
                retAge = None, usInput = None):
    """ The rows of anchor.csv (layer 0) or layer<1|2>.csv at one turnout, spec and β rule (config.US['oecd26']'s
    where None), at config phi, calibration variant and spending target (config.US['esc']['spendingTarget']), at
    thetaSource `theta` (config.US['theta'] by default; oecd26Dir, or `root`) and, for a layer, at `omegaRule`
    (oecd26OmegaRule: layer 1's ω, 'anchor' or 'mean', config's where None; layer 2 is always 'anchor'), at ν's
    retirement-age correction `retAge` (config.US['oecd26']['retAge'] where None; the anchor has none, the U.S. gap
    being 0) and at `usInput` (layer 1's input at its U.S. value, python/US/runOECD26.py --usInput; 'none', every
    country's own inputs, where None). The anchor of a turnout in ANCHORTURNOUT is that turnout's. A csv without the
    betaRule column, or a blank cell in it, is a 'common' row: the one rule before the column existed
    (runOECD26.LEGACY); one without thetaSource an 'index' row, one without omegaRule an 'anchor' row, one without
    retAge a 'none' row, one without spendingTarget a 'gross' row, one without usInput a 'none' row. """
    turnout, spec, betaRule = _oecd26Variant(turnout, spec, betaRule)
    if layer == 0:
        turnout = ANCHORTURNOUT.get(turnout, turnout)
    else:
        omegaRule = oecd26OmegaRule(layer, omegaRule)
    retAge = C.US['oecd26']['retAge'] if retAge is None else retAge
    theta = C.US['theta'] if theta is None else theta
    path = oecd26Path(layer, theta, root)
    df = C.atTheta(pd.read_csv(_need(path), dtype = {'inputsKey': str, 'anchorKey': str}), theta)
    if df.empty:
        raise MissingInput('{} (no row at θ variant {})'.format(path, theta))
    rule = df['betaRule'].fillna('common') if 'betaRule' in df.columns else pd.Series('common', index = df.index)
    sel = ((df['turnout'] == turnout) & (df['spec'] == spec) & (rule == betaRule)
           & np.isclose(df['phi'].astype(float), C.US['esc']['phi'])
           & (df['commonX'].astype(bool) == bool(C.US['commonX'])))
    sel &= df.index.isin(C.atSpending(df).index)
    if layer != 0:
        col = lambda k, v: df[k].fillna(v) if k in df.columns else pd.Series(v, index = df.index)
        sel &= ((col('omegaRule', 'anchor') == omegaRule) & (col('retAge', 'none') == retAge)
                & (col('usInput', 'none') == ('none' if usInput is None else usInput)))
    return _noDuplicates(df[sel], ['turnout'] + ([] if layer == 0 else ['iso3']), path), path


def oecd26Anchor(turnout = None, spec = None, theta = None, root = None, carrier = None):
    """ The calibrated anchor (the U.S. at the OECD inputs) of one turnout reading and cost spec (config.US['oecd26']'s
    where None; the anchor's β rule is always 'common') at thetaSource `theta` (config.US['theta'] by default), as a
    dict: β, ω, `p` (λ under 'size', the λ of λ·Π/R under 'sizePiR' (κ of 2κπ(1-π) under config.CARRIER 'convex',
    whose anchor must carry that `carrier`, else StaleInput; the layers follow through anchorKey), NaN under 'none',
    which has none), its readouts,
    `inputsKey` and `grid`. Refused unless it was solved on the fine grid and at the β the U.S. calibration without
    the cost has now at ρ = 1 at that θ (usSweep for the main variant, the variant's own sweep in oecd26Dir(theta)/
    calibration/ for a robustness one): an anchor computed before a recalibration of the U.S. is stale, and so is
    everything layered on it (#13): StaleInput. The 'index' rows (figure 1.1's index) have no live sweep since the
    main variant replaced it in results/calibration/; their β is not checked. `root`: oecd26Dir's (a directory run on
    the main variant's sweep: the stylized readings). `carrier`: the carrier the anchor must carry, config.CARRIER where
    None (another for a directory of another cost, config.US['tau1']). """
    turnout, spec, _ = _oecd26Variant(turnout, spec)
    want = C.CARRIER if carrier is None else carrier
    theta = C.US['theta'] if theta is None else theta
    df, path = _oecd26Rows(0, turnout, spec, 'common', theta, root)
    ok = df[df['status'] == 'ok']
    if ok.empty:
        raise MissingInput(path + ' (no calibrated {} {} anchor)'.format(spec, turnout))
    a = ok.iloc[0].to_dict()
    carrier = a.get('carrier', 'pi')
    carrier = 'pi' if carrier is None or carrier != carrier else carrier          # a blank cell is a 'pi' row
    if spec == 'sizePiR' and carrier != want:
        raise StaleInput('{}: the {} anchor was solved under the carrier {!r}, the reader wants {!r} (config.CARRIER, '
                         'SSD_CARRIER): rerun python/US/runOECD26.py under it'.format(os.path.relpath(path, C.REPO),
                                                                                    turnout, carrier, want))
    if theta == 'index':
        if a['grid'] != 'fine':
            raise StaleInput('{}: the {} anchor was solved on the {} grid'.format(path, turnout, a['grid']))
        return a
    if theta == C.US['theta']:
        sweep = usSweep('US')
    else:
        sp = os.path.join(oecd26Dir(theta, root), 'calibration', 'US_rhoGridCommonX.csv')
        sweep = C.atTheta(pd.read_csv(_need(sp)), theta)
        if sweep.empty:
            raise MissingInput('{} (no row at θ variant {})'.format(sp, theta))
    β = float(sweep.loc[np.isclose(sweep['ρ'], 1.0), 'β'].iloc[0])
    if a['grid'] != 'fine' or not np.isclose(float(a['β']), β, rtol = 1e-12, atol = 0):
        raise StaleInput('{}: the {} anchor was solved at β = {!r} on the {} grid; the U.S. calibration has β = {!r}: '
                         'rerun python/US/runOECD26.py'.format(os.path.relpath(path, C.REPO), turnout, a['β'],
                                                               a['grid'], β))
    return a


def oecd26Layer(layer, turnout = None, spec = None, betaRule = None, theta = None, root = None, omegaRule = None,
                retAge = None, usInput = None, carrier = None):
    """ Layer 1 or 2 of one turnout reading, cost spec and β rule (config.US['oecd26']'s where None) at thetaSource
    `theta` (config.US['theta'] by default), ν's retirement-age correction `retAge` (config's where None) and, layer 1,
    at `omegaRule` (config's where None; oecd26OmegaRule) and `usInput` ('none' where None; _oecd26Rows), from
    oecd26Dir(theta, root), indexed on ISO3. Every row must have been computed against the anchor of its spec on disk
    (its `anchorKey` is the anchor's `inputsKey`) on the anchor's grid: a row of an earlier anchor is StaleInput, never
    mixed in (crossCuttingFindings #13). `carrier`: oecd26Anchor's. """
    oecd26OmegaRule(layer, omegaRule)
    if layer != 1 and usInput not in (None, 'none'):
        raise ValueError('usInput {!r}: layer 1 only'.format(usInput))
    a = oecd26Anchor(turnout, spec, theta, root, carrier)
    df, path = _oecd26Rows(layer, turnout, spec, betaRule, theta, root, omegaRule, retAge, usInput)
    stale = (df['anchorKey'] != a['inputsKey']) | (df['grid'] != a['grid'])
    if stale.any():
        raise StaleInput('{}: rows computed against another anchor or grid: {}'.format(
            os.path.relpath(path, C.REPO), ', '.join(df.loc[stale, 'iso3'])))
    return df.set_index('iso3')


def tau1Dir(kind):
    """ results/<config.US['tau1'][kind]>/: the U.S. ('esc') or the cross-section ('oecd26') under the cost linear in the
    tax. """
    return os.path.join(C.RESULTS, C.US['tau1'][kind])


def _tau1Rows(df, path):
    """ The rows of a csv of tau1Dir at config.US['tau1']'s tauPower and carrier; StaleInput where one is another. """
    t = C.US['tau1']
    k = pd.to_numeric(df['tauPower'], errors = 'coerce') if 'tauPower' in df.columns else pd.Series(2., index = df.index)
    car = df['carrier'].fillna('pi') if 'carrier' in df.columns else pd.Series('convex', index = df.index)
    if not (np.isclose(k, t['tauPower']).all() and (car == t['carrier']).all()):
        raise StaleInput('{}: rows at another tauPower or carrier than config.US[\'tau1\'] ({}, {!r})'.format(
            os.path.relpath(path, C.REPO), t['tauPower'], t['carrier']))
    return df


def escTauPower(commonX = None, spec = None):
    """ The U.S. under log preferences (rho = 1) at the paper's cost and at the cost linear in the tax (config.US['tau1']),
    each calibrated to the same targets: {k: {'cal': the converged calibration record (escCalibration.csv), 'shocks':
    escShocks.csv's rows at the spec, indexed on (scenario, θpinned)}} for k = 2 (results/esc/) and k = 1 (tau1Dir('esc')),
    at the θ variant and spending target (_readUS), commonX and phi of config. """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    out = {}
    for k, d in ((2, C.ESCDIR), (int(C.US['tau1']['tauPower']), tau1Dir('esc'))):
        pc, ps = os.path.join(d, 'escCalibration.csv'), os.path.join(d, 'escShocks.csv')
        cal, shk = _escSelect(_readUS(pc), spec, commonX), _escSelect(_readUS(ps), spec, commonX)
        if k != 2:
            cal, shk = _tau1Rows(cal, pc), _tau1Rows(shk, ps)
        cal = cal[cal['converged'].astype(bool)]
        if len(cal) != 1:
            raise MissingInput('{} ({} converged {} rows)'.format(pc, len(cal), spec))
        shk = _noDuplicates(shk, ['scenario', 'θpinned'], ps)
        shk = shk.assign(θpinned = shk['θpinned'].astype(str) == 'True').set_index(['scenario', 'θpinned'])
        out[k] = {'cal': cal.iloc[0].to_dict(), 'shocks': shk, 'paths': (pc, ps)}
    return out


def argReading(which):
    """ Argentina's surprise reform in 2010 (the first row, t0, of shockUniversal.py's csv at rho = 1, common X) under a
    reading of the pre-reform economy: 'paper' (the basic pension, spending of the workbook), 'epsPre0' (no informal
    retiree drew a pension, REPLICATION.md section 3) or 'spending' (the basic pension, calibrated to
    config.ARG['spendingBaseline']['spending'] of GDP). Pension spending over GDP S = tau (1-alpha) Y/(Y + Y0) and the
    informal retirees' part of it S gamma_0 b0/(gamma_0 b0 + biAvg) (eq:governmentBudget; REPLICATION.md section 5), before
    (`_base`) and after (`_reform`) the reform, the tax `τ_*`, the change in formal benefits `dBenefits`, and the
    calibration's β, ω and spending datum `spendGDP`; `paths` the two csvs. """
    ρ = C.ARG['ρBaseline']
    tmpl = C.argShockTemplate('reform', commonX = True).format(ρ = ρ)
    if which == 'paper':
        shock, cal = os.path.join(C.SHOCKDIR, tmpl), C.argSweepCsv(True)
    elif which == 'epsPre0':
        shock = os.path.join(C.SHOCKDIR, tmpl.replace('.csv', '_epsPre0.csv'))
        cal = os.path.join(C.CALIBDIR, 'informalSavings_rhoGridCommonXEpsPre0.csv')
    elif which == 'spending':
        root = os.path.join(C.RESULTS, C.ARG['spendingBaseline']['dir'])
        shock = os.path.join(root, 'shocks', tmpl)
        cal = os.path.join(root, 'calibration', os.path.basename(C.argSweepCsv(True)))
    else:
        raise ValueError("which: 'paper', 'epsPre0' or 'spending', got {!r}".format(which))
    r = pd.read_csv(_need(shock), index_col = 0).iloc[0]
    c = pd.read_csv(_need(cal))
    c = c[np.isclose(c['ρ'].astype(float), ρ) & np.isfinite(pd.to_numeric(c['residual'], errors = 'coerce'))]
    if len(c) != 1:
        raise MissingInput('{} ({} solved rows at rho = {})'.format(cal, len(c), ρ))
    c = c.iloc[0]
    sm = calibrationSummary(True)
    α, γ0 = float(sm['α']), float(sm['γ0'])
    out = {'paths': (shock, cal), 'β': float(c['β']), 'ω': float(c['ω']), 'spendGDP': float(c['spendGDP']),
           'dBenefits': float(r['d_biAvg'])}
    for w in ('base', 'reform'):
        S = float(r['τ_' + w])*(1 - α)*float(r['Y_' + w])/(float(r['Y_' + w]) + float(r['Y0_' + w]))
        out['τ_' + w], out['S_' + w] = float(r['τ_' + w]), S
        out['Sinf_' + w] = S*γ0*float(r['b0_' + w])/(γ0*float(r['b0_' + w]) + float(r['biAvg_' + w]))
    if which != 'epsPre0' and not abs(out['S_base'] - out['spendGDP']) <= 1e-6:
        raise StaleInput('{}: pre-reform spending {:.6f} against the calibration\'s datum {:.6f}: the reform csv is not '
                         'of this calibration'.format(os.path.relpath(shock, C.REPO), out['S_base'], out['spendGDP']))
    return out


def _escSelect(df, spec, commonX, method = False):
    """ The rows of an ESC csv at the calibration variant (_escVariant), `spec` and config phi, and at the
    published CRRA method when `method`. A column the csv does not carry does not filter. """
    keep = []
    for i, rec in enumerate(df.to_dict('records')):
        if not _escVariant(rec, commonX):
            continue
        if spec is not None and 'spec' in rec and rec['spec'] != spec:
            continue
        if 'phi' in rec and not (rec['phi'] == rec['phi'] and np.isclose(float(rec['phi']), C.US['esc']['phi'])):
            continue
        if method and not _escMethodOK(rec):
            continue
        keep.append(i)
    return df.iloc[keep].copy()


def _noDuplicates(df, key, path):
    """ A duplicate key is an error, never resolved by first/last-wins (crossCuttingFindings #8). """
    if df.duplicated(key).any():
        raise ValueError('{}: duplicate rows on {}'.format(os.path.relpath(path, C.REPO), key))
    return df


def escPath(ρ, commonX = None, spec = None, k = 2):
    """ The baseline path of the leaded design choice at `rho`, one row per date from the first period to
    the last before the terminal one, sorted on `pos`: the design in force `θ`, the tax `τ`, `sr` (s/(w h);
    escSavingsOverY converts) and `workweek` on the path where the choice binds from the first period, and
    `τExo`, `srExo` on the same model with the design held at θ* throughout. rho = 1 (log preferences)
    from escPath.csv, any other rho from escPathCRRA.csv at the published method (escMethod). Headline
    variant and the paper's spec by default; `k` the cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    if np.isclose(ρ, 1.0):
        df, _ = _readESC('escPath.csv', k)
        path = os.path.join(escDir(k), 'escPath.csv')
        df = _escSelect(df, spec, commonX)
    else:
        df, _ = _readESC('escPathCRRA.csv', k)
        path = os.path.join(escDir(k), 'escPathCRRA.csv')
        df = _escSelect(df[np.isclose(df['ρ'], ρ)], spec, commonX, method = True)
    if df.empty:
        raise MissingInput('{} (ρ={}, {}, commonX={}, method={})'.format(path, ρ, spec, commonX, escMethod()))
    return _noDuplicates(df, ['pos'], path).sort_values('pos').reset_index(drop = True)


def escPathHost(commonX = None, spec = None, k = 2):
    """ python/US/runESC.stagePathHost's csv: escPath's columns for France and the UK (log preferences), per
    `country` and `wedgeFrom` ('US': the US-calibrated cost parameter; 'own': the country's own, which France
    does not have), plus the country's θStar and ω. Headline variant and the paper's spec by default; `k` the
    cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    df, _ = _readESC('escPathHost.csv', k)
    path = os.path.join(escDir(k), 'escPathHost.csv')
    df = _escSelect(df, spec, commonX)
    if df.empty:
        raise MissingInput('{} ({}, commonX={})'.format(path, spec, commonX))
    return _noDuplicates(df, ['country', 'wedgeFrom', 'pos'], path).sort_values(['country', 'wedgeFrom', 'pos'])


def escXi(commonX = None, spec = None, k = 2):
    """ python/US/runESCxi.py's csv at the variant and spec, log preferences: per Frisch elasticity `xi`
    one 'calibration' row (p = lambda, θStar, β, ω, X, converged), one 'path' row per date (the columns of
    escPath) and the rows 'acute pinned' and 'acute chosen' (ν_t = 1 throughout, a new equilibrium path
    read at 2020: θ, τ, sr = s/(w h), workweek). `k`: the cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    df, _ = _readESC('escXiRobustness.csv', k)
    path = os.path.join(escDir(k), 'escXiRobustness.csv')
    df = _escSelect(df, spec, commonX)
    if df.empty:
        raise MissingInput('{} ({}, commonX={})'.format(path, spec, commonX))
    return _noDuplicates(df, ['xi', 'kind', 'pos'], path)


def escTauTheta(commonX = None, spec = None, k = 2):
    """ python/US/runTauTheta.py's csv (results/esc/escTauTheta.csv) at the θ variant and calibration variant: the U.S.
    at ρ = 1 with θ pinned on the whole path at `θpin`, without the cost (spec 'none') and under `spec`
    (config.US['esc']['spec'] by default) at its calibrated λ, as {'none': frame, spec: frame}, each sorted on θpin
    (τ, spending in % of GDP, `kind` 'grid' | 'refine'). Refused (StaleInput) where the costly rows' λ is not
    escCalibration's at ρ = 1, or the rows without the cost were calibrated at another β than the U.S. sweep's (#13).
    `k`: the cost's power of the tax (escDir). """
    commonX = C.US['commonX'] if commonX is None else commonX
    spec = C.US['esc']['spec'] if spec is None else spec
    df, _ = _readESC('escTauTheta.csv', k)
    path = os.path.join(escDir(k), 'escTauTheta.csv')
    df = df[df['commonX'].astype(bool) == bool(commonX)]
    out = {}
    for s in ('none', spec):
        d = df[df['spec'] == s]
        if s != 'none':
            d = d[np.isclose(d['phi'].astype(float), C.US['esc']['phi'])]
        if d.empty:
            raise MissingInput('{} (no {} rows, commonX={})'.format(path, s, commonX))
        out[s] = _noDuplicates(d, ['θpin'], path).sort_values('θpin').reset_index(drop = True)
    lam = float(escCalibrationLOG(spec, commonX, k)['p'])
    if not np.allclose(out[spec]['p'].astype(float), lam, rtol = 1e-12, atol = 0):
        raise StaleInput('{}: the {} rows were run at λ = {}, escCalibration.csv has {!r}: rerun python/US/runTauTheta.py'
                         .format(os.path.relpath(path, C.REPO), spec, sorted(set(out[spec]['p'])), lam))
    sw = usSweep('US', commonX)
    β = float(sw.loc[np.isclose(sw['ρ'], 1.0), 'β'].iloc[0])
    if not np.allclose(out['none']['β'].astype(float), β, rtol = 1e-10, atol = 0):
        raise StaleInput('{}: the rows without the cost were calibrated at β = {}, the U.S. sweep has {!r}: rerun '
                         'python/US/runTauTheta.py'.format(os.path.relpath(path, C.REPO), sorted(set(out['none']['β'])), β))
    return out


def _variantRows(path, commonX, key):
    commonX = C.US['commonX'] if commonX is None else commonX
    df = _readUS(path)
    df = df.iloc[[i for i, rec in enumerate(df.to_dict('records')) if _escVariant(rec, commonX)]]
    if df.empty:
        raise MissingInput('{} (commonX={})'.format(path, commonX))
    return _noDuplicates(df, key, path).reset_index(drop = True)


def stationaryApprox(which):
    """ The stationary-against-date-specific check of `which`: 'US' (the tax with the design given, CRRA),
    'US_ESC' (the design and the tax with the design chosen; `method` LOG or CRRA2D) or 'ARG' (Argentina's
    tax). One row per (rho, t) with its `date`; `gap_pp` / `gapτ_pp` is the tax on the path walked by the
    stationary policy functions minus the date-specific one, in p.p., `gapθ` the same for the design. Run
    on the headline variant only (python/paper/runShocksUS.py and runShocks.py, entry 'stationary'). """
    name = {'US': 'US_stationaryApprox.csv', 'US_ESC': 'US_ESC_stationaryApprox.csv',
            'ARG': 'ARG_stationaryApprox' + C.argVariantTag(C.ARG['commonX']) + '.csv'}[which]
    path = os.path.join(C.NUMDIR, name)
    df = pd.read_csv(_need(path)) if which == 'ARG' else _readUS(path)
    return _noDuplicates(df, ['ρ', 't'], path).sort_values(['ρ', 't']).reset_index(drop = True)


HORIZONCSVS = {'US': ('US_horizon.csv', 'US_horizonDesign.csv'),
               'ARG': ('ARG_horizon' + C.argVariantTag(C.ARG['commonX']) + '.csv',)}
RESIDUALCSVS = {'US': 'US_residuals.csv', 'ARG': 'ARG_residuals' + C.argVariantTag(C.ARG['commonX']) + '.csv'}


def horizon(arm):
    """ The horizon and terminal-condition runs of appendix D (python/US/runHorizon.py stages log, horizon,
    design; python/InformalSavings/runHorizon.py), at the calibrated parameters of the headline variant. One
    row per (block, preferences, rho, terminal, tss, date), dates 1960-2110 (U.S.) or 1920-2100 (Argentina):
    `terminal` 'finite' with `tss` periods after the last projection date (`T` the last year; tss 3 is the
    paper's horizon) or 'steadyState'; `block` 'tax' (the design given) or 'design' (chosen; U.S. only);
    Argentina's rows are block 'tax' and carry the universal reform (`τ_reform`, `dτ`, `srGDP_reform`, `dsrGDP`).
    Savings: `srOverY` s/Y (U.S.), `srGDP` s/(Y + Y0) (Argentina). """
    frames = []
    for f in HORIZONCSVS[arm]:
        p = os.path.join(C.NUMDIR, f)
        d = _readUS(p) if arm == 'US' else pd.read_csv(_need(p))     # Argentina's csvs carry no thetaSource
        if arm == 'ARG':
            d = d.assign(block = 'tax')
        frames.append(d)
    df = pd.concat(frames, ignore_index = True)
    path = os.path.join(C.NUMDIR, HORIZONCSVS[arm][0])
    return _noDuplicates(df, ['block', 'preferences', 'ρ', 'terminal', 'tss', 'date'], path)


def residuals(arm):
    """ The first-order-condition residuals of appendix D at states off the policy functions' grids (the
    runHorizon.py stage residuals), summary statistics of |tau_FOC - tau_hat| in units of the tax: per (rho,
    ns [U.S.], date, set, region), `date` 'all' pooling the periods; `max`, `p99`, `median`, `n` against the
    condition as the solver evaluates it, the same with suffix `Ref` against the condition with converged tax
    derivatives. `set`: 'node', 'cell', 'edge' (Argentina's 2-D grid), 'path'; `region`: 'all', 'visited'. """
    path = os.path.join(C.NUMDIR, RESIDUALCSVS[arm])
    df = _readUS(path) if arm == 'US' else pd.read_csv(_need(path))
    if 'ns' not in df.columns:
        df['ns'] = np.nan
    df['date'] = df['date'].astype(str)
    return _noDuplicates(df, ['ρ', 'ns', 'date', 'set', 'region'], path)


def taxPolicyState():
    """ The U.S. 2020 tax policy function against inherited savings (python/US/c5Numbers.py a13): per (rho,
    `s_2019/s*`), `τ2020` the policy at the savings the 2020 economy inherits, as a multiple of their value on the
    equilibrium path, and `τ2020Equilibrium` the path's tax; design given, the calibration of each rho, NaN off
    the policy's state grid. """
    path = os.path.join(C.NUMDIR, 'C5_taxPolicyState.csv')
    return _noDuplicates(_readUS(path), ['ρ', 's_2019/s*'], path).sort_values(['ρ', 's_2019/s*'])


def strategicTerm():
    """ The strategic term of the U.S. 2020 choice (python/US/c5Numbers.py strategic): per rho, `dsdτ` the
    economic equilibrium's ds_2020/dtau_2020 at the inherited savings of the path, `dτ1ds` the slope of the
    2050 policy function there, `strategic` their product, and the three as elasticities; under log only the
    zeros of prop:LOG:PEE (the rest NaN). """
    path = os.path.join(C.NUMDIR, 'C5_strategicTerm.csv')
    return _noDuplicates(_readUS(path), ['ρ'], path).sort_values('ρ')


def criticalTurnout():
    """ The sign of eq:esc:seqFOC without the cost, the design chosen with the tax (python/US/c5Numbers.py a1;
    prop:esc:corner, app:esc:corner): per (rho, date, tax), on the calibrated equilibrium of each dated period,
    `negativeEverywhere` over theta_t in [0, 1], `turnoutRatio` the observed mu_H/mu_L and the critical ratio at
    which the derivative first turns non-negative somewhere, weights mu_i = (mu_i^CPS)^m (`ratioCrit_μm`) or
    y_i^k (`ratioCrit_yk`). `tax` 'calibrated' evaluates at the path's tax of that date (the paper's numbers);
    'theta0', at t0 only, at the tax of the theta = 0 counterfactual (US_shocksCommonX.csv). """
    path = os.path.join(C.NUMDIR, 'C5_criticalTurnout.csv')
    return _noDuplicates(_readUS(path), ['ρ', 'date', 'tax'], path).sort_values(['ρ', 'date', 'tax'])


def numDir(k = 2):
    """ The directory of C5's numbers (python/US/c5Numbers.py) under the cost with tau^k: results/numerical/ at k = 2,
    tau1Dir('numerical') at config.US['tau1']'s power. """
    if k == 2:
        return C.NUMDIR
    if k != C.US['tau1']['tauPower']:
        raise ValueError("numDir: k = {} (2 or config.US['tau1']['tauPower'])".format(k))
    return tau1Dir('numerical')


def _readNum(name, k = 2):
    """ _readUS of numDir(k)/<name>, at k other than 2 at config.US['tau1']'s power and carrier only (_tau1Rows).
    Returns (frame, its path). """
    path = os.path.join(numDir(k), name)
    df = _readUS(path)
    return (df if k == 2 else _tau1Rows(df, path)), path


def designDecomp(k = 2):
    """ The design choice of 2020 under log preferences (python/US/c5Numbers.py design; eq:esc:leadFOC at the 1990
    electorate). `kind` 'foc': per (`case` US/FR, `route`, `piece`), `value` dW_t/dtheta_{t+1} at the U.S. choice
    `θ1`, route 'total', 'design' (the benefit formula's theta, the share lost and the tax held), 'cost' (the share
    lost's theta alone, eq:esc:ftheta) or 'tax' (the anticipated tax T_{t+1}(theta)); piece 'W' the sum of the four
    blocs. `kind` 'choice': per `experiment`, `θ2020` the design in force in 2020 on the free path: 'US', 'FR' (France's
    income distribution), 'X1' (FR with the cost at the U.S. size), 'X2' (at the U.S. size and dispersion), 'X3' (the
    U.S. with France's dispersion in the cost), each with the `gap` of the cost to its target. `k`: the cost's power of
    the tax (numDir). """
    df, path = _readNum('C5_designDecomp.csv', k)
    return _noDuplicates(df, ['kind', 'case', 'experiment', 'route', 'piece'], path)


def designDecompCRRA(ρ, k = 2):
    """ c5Numbers.py designCRRA at `rho` (the exact 2-D recursion): designDecomp's 'choice' rows X1 and X2. `k`: the
    cost's power of the tax (numDir). """
    df, path = _readNum('C5_designDecompCRRA_rho{:g}.csv'.format(ρ), k)
    return _noDuplicates(df, ['experiment'], path)


def lambdaElasticity(k = 2):
    """ The calibrated lambda against participation and taxable-income elasticities (c5Numbers.py lambda; appendix
    app:esc:cost), per rho: `πBaseline2020` = Pi/R in 2020 without the cost, `λ`. At k = 2 (eps(T) = kappa T) `τ2020`,
    `tbar2020` the wedge in place eq:esc:wedge, the suffix _A under the carrier lambda = 2 kappa pi (1 - pi), _C under
    the fiscal externality's 3 kappa pi (1 - pi): `εLocal_*` the local elasticity at the wedge in place, `baseOutside_*`
    the share of earnings outside the base, `λ_ε<e>_A` the lambda an elasticity e implies, `multiple_ε<e>_*` lambda
    over it. At k = 1 (a constant elasticity, lambda = eps pi) the suffix _ΠR: `εImplied_ΠR` = lambda/pi, `λ_ε<e>_ΠR` =
    e pi, `multiple_ε<e>_ΠR` lambda over it. `k`: the cost's power of the tax (numDir). """
    df, path = _readNum('C5_lambdaElasticity.csv', k)
    return _noDuplicates(df, ['ρ'], path).sort_values('ρ')


def selectionCsvs():
    """ The results csvs the registered outputs read that a politico-economic solve wrote, for the
    equilibrium counts (NUM_Selection): [(arm, label, [paths])], one entry per file or per rho-grid family.
    The endogenous-design experiments appear through the csvs escExperiments.csv merges, since the merge
    drops the counts (python/US/collectESCexperiments.py). The cross-country test's csvs (results/oecd26, labelled
    'oecd26/<file>') carry -1 on the rows without a solve (status 'noTurnout': Luxembourg at own turnout), which
    countSummary leaves uncounted, as it does a failed cost calibration's. """
    arg, us = [], []
    for cx in (C.ARG['commonX'], not C.ARG['commonX']):
        arg.append(('ARG', os.path.basename(C.argSweepCsv(cx)), [C.argSweepCsv(cx)]))
        for scen in ('reform', 'ee'):
            tmpl = C.argShockTemplate(scen, commonX = cx)
            arg.append(('ARG', tmpl.replace('{ρ:.4f}', '*'),
                        [os.path.join(C.SHOCKDIR, tmpl.format(ρ = ρ)) for ρ in C.ARG['ρGrid']]))
        arg.append(('ARG', os.path.basename(C.argEpsThetaCsv(C.ARG['ρBaseline'], cx)),
                    [C.argEpsThetaCsv(C.ARG['ρBaseline'], cx)]))
    arg.append(('ARG', 'ARG_stationaryApprox' + C.argVariantTag(C.ARG['commonX']) + '.csv',
                [os.path.join(C.NUMDIR, 'ARG_stationaryApprox' + C.argVariantTag(C.ARG['commonX']) + '.csv')]))
    for f in HORIZONCSVS['ARG']:
        arg.append(('ARG', f, [os.path.join(C.NUMDIR, f)]))
    sb = C.ARG['spendingBaseline']['dir']                  # page A1's reform from another pre-reform spending
    tmpl = C.argShockTemplate('reform', commonX = True).format(ρ = C.ARG['ρBaseline'])
    for sub, f in (('calibration', os.path.basename(C.argSweepCsv(True))), ('shocks', tmpl)):
        arg.append(('ARG', sb + '/' + f, [os.path.join(C.RESULTS, sb, sub, f)]))
    for cx in (C.US['commonX'], not C.US['commonX']):
        for country in ('US', 'UK', 'FR'):
            us.append(('US', os.path.basename(C.usSweepCsv(country, cx)), [C.usSweepCsv(country, cx)]))
        us.append(('US', os.path.basename(C.usShockCsv('US', cx)), [C.usShockCsv('US', cx)]))
    us.append(('US', 'US_shocksXiCommonX.csv', [os.path.join(C.SHOCKDIR, 'US_shocksXiCommonX.csv')]))
    escFiles = ('escCalibration.csv', 'escCalibrationCosts.csv', 'escCalibrationCRRA.csv', 'escCalibrationCRRAUK.csv',
                'escCountry.csv', 'escShocks.csv', 'escShocksCosts.csv', 'escShocksCRRA.csv', 'escPath.csv',
                'escPathCRRA.csv', 'escPathHost.csv', 'escTauTheta.csv', 'escXiRobustness.csv')
    for f in escFiles:
        us.append(('US', f, [os.path.join(C.ESCDIR, f)]))
    for layer in (0, 1, 2):
        p = oecd26Path(layer)
        us.append(('US', 'oecd26/' + os.path.basename(p), [p]))
    for V in C.US['oecd26']['thetaVariants']:
        for layer in (0, 1, 2):
            p = oecd26Path(layer, V)
            us.append(('US', 'oecd26_theta_{}/{}'.format(V, os.path.basename(p)), [p]))
    for f in ('US_stationaryApprox.csv', 'US_ESC_stationaryApprox.csv') + HORIZONCSVS['US'] + ('US_derivativeScheme.csv',):
        us.append(('US', f, [os.path.join(C.NUMDIR, f)]))
    for f in escFiles:                                          # the cost linear in the tax (config.US['tau1'])
        us.append(('US', C.US['tau1']['esc'] + '/' + f, [os.path.join(tau1Dir('esc'), f)]))
    for f in ('anchor.csv', 'layer1.csv', 'layer2.csv'):
        us.append(('US', C.US['tau1']['oecd26'] + '/' + f, [os.path.join(tau1Dir('oecd26'), f)]))
    return arg + us


COUNTCOLUMNS = {'tax':    ('nEqMax', 'nCandMax', 'nFallback'),
                'design': ('nEqθMax', 'nBrθMax', 'nFallbackθ')}
# The root counts of the state's own fixed point (eq:stateResidual; python/US/policy.rootSummary,
# python/InformalSavings/policy.ROOTS): the savings level s_t under CRRA, the informal share ι_t in Argentina; -1 where
# the solve has no such fixed point (the U.S. arm under log), blank on a row written before the columns.
ROOTCOLUMNS = ('nRootsSMax', 'nRootsιMax')


def countSummary(paths, layer = 'tax'):
    """ The equilibrium counts over the csvs `paths`: {'files', 'rows', 'counted', 'uncounted', 'unrecorded',
    'eqMax', 'candMax', 'fallbacks'} for the tax rule (layer 'tax': nEqMax, nCandMax, nFallback,
    writing/US/num_robustroot.tex) or the design layer of the exact CRRA recursion (layer 'design': nEqθMax,
    nBrθMax, nFallbackθ, writing/US/num_esc.tex). A row is counted where its three columns are present and
    non-negative; 'uncounted' rows carry -1, what a solve that selected nothing at that layer writes
    (python/US/policy.multiplicityColumns); 'unrecorded' rows are blank, written before the columns existed.
    With nothing counted the three statistics are None, not zero. 'rootsMax' is the largest root count of the
    state's fixed point (ROOTCOLUMNS) over the rows that record a non-negative one, whatever the layer, 'rootsRows' the
    number of those rows; None and 0 where no row does. """
    lack = [p for p in paths if not os.path.exists(p)]
    if lack:
        raise MissingInput(lack)
    cols, rows, parts, neg, blank = COUNTCOLUMNS[layer], 0, [], 0, 0
    roots = []
    for p in paths:
        d = _countFrame(p)
        rows += len(d)
        rc = [c for c in ROOTCOLUMNS if c in d.columns]
        if rc:
            r = d[rc].apply(pd.to_numeric, errors = 'coerce').max(axis = 1)
            roots.append(r[r >= 0])
        if all(c in d.columns for c in cols):
            sub = d[list(cols)].apply(pd.to_numeric, errors = 'coerce')
            isBlank = sub.isna().any(axis = 1)
            isNeg = ~isBlank & (sub < 0).any(axis = 1)
            neg, blank = neg + int(isNeg.sum()), blank + int(isBlank.sum())
            parts.append(sub[~isBlank & ~isNeg])
        else:
            blank += len(d)
    sub = pd.concat(parts) if parts else pd.DataFrame(columns = list(cols))
    none = sub.empty
    roots = pd.concat(roots) if roots else pd.Series(dtype = float)
    return {'files': len(paths), 'rows': rows, 'counted': len(sub), 'uncounted': neg, 'unrecorded': blank,
            'eqMax': None if none else int(sub[cols[0]].max()),
            'candMax': None if none else int(sub[cols[1]].max()),
            'fallbacks': None if none else int(sub[cols[2]].sum()),
            'rootsMax': None if roots.empty else int(roots.max()), 'rootsRows': int(len(roots))}


COUNTIDS = ('ρ', 'spec', 'scenario', 'θpinned', 'method', 'family', 'effect', 'terminal', 'tss', 'date')


def thetaKeyed(path):
    """ Is `path` a U.S.-arm csv of a run at the θ variant (the U.S., France and UK sweeps and shocks, results/esc,
    results/oecd26, the U.S. stationary checks), whose rows are read at config.US['theta'] only? """
    p, name = os.path.abspath(path), os.path.basename(path)
    under = lambda d: os.path.normcase(p).startswith(os.path.normcase(os.path.abspath(d)) + os.sep)
    return ((under(C.CALIBDIR) and name.startswith(('US_', 'FR_', 'UK_'))) or under(C.ESCDIR) or under(OECD26DIR)
            or (under(C.SHOCKDIR) and name.startswith('US_'))
            or name in ('US_stationaryApprox.csv', 'US_ESC_stationaryApprox.csv'))


def _countFrame(p):
    """ The rows of csv `p` the counts read: at the θ variant where the file is thetaKeyed, at its own variant in a
    robustness variant's directory (results/oecd26_theta_<V>/), all of them otherwise; and at the spending target where
    the file is config.spendingKeyed. """
    d = pd.read_csv(_need(p))
    top = os.path.basename(os.path.dirname(os.path.abspath(p)))
    if top.startswith('oecd26_theta_'):
        d = C.atTheta(d, top[len('oecd26_theta_'):])
    elif thetaKeyed(p):
        d = C.atTheta(d)
    return C.atSpending(d) if C.spendingKeyed(p) else d


def countExceptions(paths, layer = 'tax'):
    """ The counted rows of `paths` with more than one equilibrium or a fallback at `layer`: [(file, {its
    identifying columns among COUNTIDS}, nEq max, fallbacks)]. """
    cols, out = COUNTCOLUMNS[layer], []
    for p in paths:
        d = _countFrame(p)
        if not all(c in d.columns for c in cols):
            continue
        sub = d[list(cols)].apply(pd.to_numeric, errors = 'coerce')
        hit = sub.notna().all(axis = 1) & (sub >= 0).all(axis = 1) & ((sub[cols[0]] > 1) | (sub[cols[2]] > 0))
        for i in d.index[hit]:
            ids = {k: d.at[i, k] for k in COUNTIDS if k in d.columns and pd.notna(d.at[i, k])}
            out.append((os.path.basename(p), ids, int(sub.at[i, cols[0]]), int(sub.at[i, cols[2]])))
    return out

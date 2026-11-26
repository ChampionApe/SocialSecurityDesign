r""" The supplemental appendix's datasets: one tidy csv per interactive figure of the site
(writing/OnlineAppendix/assets/oa.js), also offered for download beside it.

Each builder reads results/ through datasets.py only and returns a DataFrame; build.py writes it to
results/paper/Data/<name>.csv (kind 'data', ONLINE_ONLY) and `build.py --site` copies it into the site.

Units are the paper's: tax and savings rates in percent and their changes in percentage points, savings over
GDP, the workweek in hours, the design theta and workers per retiree as levels. A change is against the
baseline of the same economy, rho and xi; with the design chosen, against the baseline of the same reading,
pinned against pinned and chosen against chosen (figure 6). Time paths run from the reform year (Argentina) or
from one period before the calibration year (the U.S., the UK, France) to the last dated period of the
calendar: the first period, the steady-state tail and the terminal period are left out.
"""
import os
import numpy as np, pandas as pd

import config as C
import datasets as D
from tablesUS import ROWKEY

# The ESC csvs' scenario names -> the site's keys (ROWKEY's, for the rows the two exercises share).
ESCKEY = {'acute': 'acute', 'mild': 'mild', 'frIncome': 'income', 'frVoting': 'voting', 'frLeisure': 'leisure',
          'frBoth': 'incomeVoting', 'frAll': 'frall', 'France': 'france'}
# A consistency check between two files that hold the same reading (levels; theta, tau and s/Y as fractions).
CHECK_TOL = 1e-6


def _note(text):
    """ A deliberate omission, printed with the build's other lines. """
    print('  siteData: ' + text)


def _workbookXi():
    """ The Frisch elasticity of every run without an explicit xi: the workbook's, read off the calibration
    summary rather than typed. """
    xs = {float(r['ξ']) for r in D.usCalibrationSummary().values()}
    if len(xs) != 1:
        raise ValueError('usCalibrationSummary.csv: more than one ξ across countries: {}'.format(sorted(xs)))
    return xs.pop()


# ---------------------------------------------------------------------------------------------------
# Argentina
# ---------------------------------------------------------------------------------------------------
def argReformPaths():
    """ ARG_ReformPaths: the reform of table 2 at every rho of the grid. Per period, the pre-reform and the
    reform path of the tax rate, savings over GDP, the average workweek (normalised as table 2, so the
    pre-reform path reads the observed workweek in the reform year at every rho) and the informal savings
    ratio. Headline calibration variant. """
    cal = C.calendar()
    years = [y for k, y in sorted(cal['dates'].items()) if k >= cal['t0']]
    rows = []
    for ρ in C.ARG['ρGrid']:
        p = D.shockPath(ρ, 'reform')
        if len(p) < len(years):
            raise D.MissingInput('{} periods in the reform path at ρ={}, {} dated'.format(len(p), ρ, len(years)))
        hRef = D.baselineHours(ρ)
        sr = {w: D.savingsRatePath(ρ, w).values for w in ('base', 'reform')}
        for k, year in enumerate(years):
            r = p.iloc[k]
            vals = {'tau': (100*r['τ_base'], 100*r['τ_reform']),
                    'sr': (100*sr['base'][k], 100*sr['reform'][k]),
                    'workweek': tuple(float(C.workweekHours(r['h_' + w], hRef)) for w in ('base', 'reform')),
                    'iota': (r['ι_base'], r['ι_reform'])}
            for v, (pre, post) in vals.items():
                rows.append({'rho': round(float(ρ), 1), 'year': int(year), 'variable': v,
                             'pre': float(pre), 'reform': float(post)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------
# The rich economies with the design given
# ---------------------------------------------------------------------------------------------------
# (site key, shock csv column, scale to the site's unit)
OECDOUTCOMES = (('tau', 'τ', 100.), ('sr', D.US_SR, 100.), ('workweek', 'workweek', 1.))


def _oecdBlock(df, source, host, calibration, ξ):
    """ The long rows of one shock csv block (one rho, or one xi at rho = 1): every scenario and effect
    against the block's baseline. """
    rows = []
    for ρ in sorted(df['ρ'].unique()):
        sub = df[np.isclose(df['ρ'], ρ)]
        base = D.usBaseline(sub, ρ)
        for _, r in sub[sub['effect'] != 'baseline'].iterrows():
            for key, col, scale in OECDOUTCOMES:
                lev, b = scale*float(r[col]), scale*float(base[col])
                rows.append({'source': source, 'host': host, 'calibration': calibration, 'xi': ξ,
                             'rho': round(float(ρ), 1), 'scenario': ROWKEY[r['scenario']], 'effect': r['effect'],
                             'outcome': key, 'level': lev, 'baseline': b, 'change': lev - b})
    return rows


OECDHOSTS = ('US',)


def oecdCounterfactuals():
    """ OECD_Counterfactuals: every counterfactual of section 6 on the U.S., under both calibrations of the taste for
    leisure and at the three rho of the paper ('main'); and on the U.S. at three Frisch elasticities, rho = 1, common
    X ('xi'; python/US/runShocksXi.py). Both effects: 'full' and 'ee' (the tax held at its baseline path). """
    ξ0 = _workbookXi()
    rows = []
    for host in OECDHOSTS:
        for cx in (True, False):
            rows += _oecdBlock(D.usShocks(commonX = cx, host = host), 'main', host,
                               'commonX' if cx else 'vectorX', ξ0)
    xi = D.usShocksXi()
    for ξ in sorted(xi['ξ'].unique()):
        rows += _oecdBlock(xi[np.isclose(xi['ξ'], ξ)], 'xi', 'US', 'commonX', round(float(ξ), 2))
    out = pd.DataFrame(rows)
    # the xi file's workbook-xi block is the main file's rho = 1 block (runShocksXi.check); keep it honest
    a = out[(out['source'] == 'xi') & np.isclose(out['xi'], ξ0)].set_index(['scenario', 'effect', 'outcome'])['level']
    b = out[(out['source'] == 'main') & (out['host'] == 'US') & (out['calibration'] == 'commonX')
            & np.isclose(out['rho'], 1.0)].set_index(['scenario', 'effect', 'outcome'])['level'].reindex(a.index)
    if not (a.notna().all() and b.notna().all()):
        raise D.StaleInput('US_shocksXiCommonX.csv at ξ = {} and US_shocksCommonX.csv at ρ = 1 do not hold the same '
                           'rows ({} of {} compared)'.format(ξ0, int((a.notna() & b.notna()).sum()), len(a)))
    gap = (a - b).abs().max()
    if not gap <= 100*CHECK_TOL:
        raise ValueError('US_shocksXiCommonX.csv at ξ = {} differs from US_shocksCommonX.csv at ρ = 1 by {}'.format(ξ0, gap))
    return out


# ---------------------------------------------------------------------------------------------------
# The endogenous design
# ---------------------------------------------------------------------------------------------------
def _escLevels(r, θ = 'θ_t0', τ = 'τ_t0', sr = 'sr_t0', ww = 'ww_t0'):
    """ The four site outcomes of one ESC reading as levels in the site's units; savings over GDP. """
    return {'theta': float(r[θ]), 'tau': 100*float(r[τ]), 'sr': 100*float(D.escSavingsOverY(float(r[sr]))),
            'workweek': float(r[ww])}


def _escRows(source, cost, host, spec, ξ, ρ, scenario, pinned, lev, base):
    return [{'source': source, 'cost': cost, 'host': host, 'spec': spec, 'xi': ξ, 'rho': round(float(ρ), 1),
             'scenario': scenario, 'reading': 'pinned' if pinned else 'chosen', 'outcome': o,
             'level': lev[o], 'baseline': base[o], 'change': lev[o] - base[o]} for o in lev]


# Section 7's two costs (config.COSTKEYS, the datasets' `cost` column), by the power k of the tax rate: every
# endogenous-design dataset carries the rows of both, each read from its own directory (datasets.escDir).
ESCCOSTS = (1, 2)


def _escCounterfactualsAt(k, ξ0):
    """ ESC_Counterfactuals' rows under the cost with tau^k (its `cost` config.COSTKEYS[k]). """
    cost = C.COSTKEYS[k]
    rows, left = [], []
    for host in ('US',):
        df = D.escExperiments(k)
        for spec in sorted(df['spec'].unique()):
            for ρ in C.US['esc']['ρTable']:
                for pinned in (True, False):
                    try:
                        base = _escLevels(D.escRow(df, ρ, spec, 'baseline', pinned))
                    except D.MissingInput:
                        left.append((spec, ρ, 'baseline', pinned))
                        continue
                    for scen, key in ESCKEY.items():
                        try:
                            r = D.escRow(df, ρ, spec, scen, pinned)
                        except D.MissingInput:          # France has no chosen reading; mild ageing is rho = 1 only
                            left.append((spec, ρ, scen, pinned))
                            continue
                        rows += _escRows('main', cost, host, spec, ξ0, ρ, key, pinned, _escLevels(r), base)
    if left:
        _note('ESC_Counterfactuals ({}): no escExperiments.csv row for {} of (spec, rho, scenario, reading); left out'
              .format(cost, len(left)))
    xi = D.escXi(k = k)                     # required: page C3's prose states the xi results under both costs
    year0 = C.usCalendar()['year0']
    for ξ in sorted(xi['xi'].unique()):
        sub = xi[np.isclose(xi['xi'], ξ)]
        cal = sub[sub['kind'] == 'calibration'].iloc[0]
        p20 = sub[(sub['kind'] == 'path') & (sub['date'] == year0)].iloc[0]
        base = {False: _escLevels(p20, θ = 'θ', τ = 'τ', sr = 'sr', ww = 'workweek'),
                True: _escLevels(pd.Series({'θ': cal['θStar'], 'τ': p20['τExo'], 'sr': p20['srExo'],
                                            'ww': p20['workweekExo']}), θ = 'θ', τ = 'τ', sr = 'sr', ww = 'ww')}
        for pinned in (True, False):
            r = sub[sub['kind'] == 'acute ' + ('pinned' if pinned else 'chosen')].iloc[0]
            lev = _escLevels(r, θ = 'θ', τ = 'τ', sr = 'sr', ww = 'workweek')
            rows += _escRows('xi', cost, 'US', C.US['esc']['spec'], round(float(ξ), 2), 1.0, 'acute', pinned, lev,
                             base[pinned])
    return rows


def escCounterfactuals():
    """ ESC_Counterfactuals: section 7's counterfactuals on the U.S. under both cost specifications, each read
    with the design pinned and chosen, at the three rho of the paper ('main'); and acute ageing on the U.S. at
    three Frisch elasticities, rho = 1 ('xi'; escXiRobustness.csv, whose baselines are its path rows in 2020;
    python/US/runESCxi.py, stage (ii)'s entry escXi, by hand under the linear cost). Each under both of section 7's
    costs (`cost`, ESCCOSTS); the cost on the design does not depend on k and is in both. The xi block is required
    under both (MissingInput otherwise): page C3 and the figure's note state its results. The design is a level; the
    tax rate, savings over GDP and the workweek are also changes against the baseline of the same reading. """
    ξ0 = _workbookXi()
    out = pd.DataFrame(sum((_escCounterfactualsAt(k, ξ0) for k in ESCCOSTS), []))
    # the xi file at the workbook xi is the main file at rho = 1 (runESCxi.py's check); keep the two honest
    key = ['scenario', 'reading', 'outcome']
    for k in ESCCOSTS:
        at = out[out['cost'] == C.COSTKEYS[k]]
        a = at[(at['source'] == 'xi') & np.isclose(at['xi'], ξ0)].set_index(key)[['level', 'baseline']]
        if a.empty:
            raise D.MissingInput('{}/escXiRobustness.csv (no rows at the workbook ξ = {})'
                                 .format(os.path.basename(D.escDir(k)), ξ0))
        b = at[(at['source'] == 'main') & (at['host'] == 'US') & (at['spec'] == C.US['esc']['spec'])
               & np.isclose(at['rho'], 1.0)].set_index(key)[['level', 'baseline']].reindex(a.index)
        if not (a.notna().all().all() and b.notna().all().all()):
            raise D.StaleInput('escXiRobustness.csv at ξ = {} and escExperiments.csv at ρ = 1 do not hold the same '
                               'rows ({})'.format(ξ0, C.COSTKEYS[k]))
        gap = float((a - b).abs().max().max())
        if not gap <= 100*CHECK_TOL:
            raise ValueError('escXiRobustness.csv at ξ = {} differs from escExperiments.csv at ρ = 1 by {} ({})'
                             .format(ξ0, gap, C.COSTKEYS[k]))
    return out


def _pathYears():
    """ The dated years from one period before the calibration year to the last dated one. """
    cal = C.usCalendar()
    return [y for k, y in sorted(cal['dates'].items()) if k >= cal['t0'] - 1]


def _pathRows(df, source, cost, country, wedges, spec, ρ, ξ, θStar, years):
    """ The long rows of one baseline path (escPath's columns): per year and variable, the reading with the
    design chosen and with it pinned at theta*. `wedges`: the wedgeFrom labels the rows go under; `cost` the
    config.COSTKEYS key of the cost the path is under. """
    rows = []
    for _, r in df[df['date'].isin(years)].iterrows():
        sr = lambda c: 100*float(D.escSavingsOverY(float(r[c])))
        vals = {'nu': (float(r['ν']), float(r['ν'])), 'theta': (float(r['θ']), float(θStar)),
                'tau': (100*float(r['τ']), 100*float(r['τExo'])), 'sr': (sr('sr'), sr('srExo')),
                'workweek': (float(r['workweek']),
                             float(r['workweekExo']) if 'workweekExo' in r and r['workweekExo'] == r['workweekExo']
                             else np.nan)}
        for w in wedges:
            for v, (ch, pin) in vals.items():
                rows.append({'source': source, 'cost': cost, 'country': country, 'wedgeFrom': w, 'spec': spec,
                             'rho': round(float(ρ), 1), 'xi': ξ, 'year': int(r['date']), 'variable': v,
                             'chosen': ch, 'pinned': pin})
    return rows


def _escPathsAt(k, years, ξ0):
    """ ESC_Paths' rows under the cost with tau^k (its `cost` config.COSTKEYS[k], every csv of datasets.escDir(k)). """
    cost = C.COSTKEYS[k]
    calib = D.escCalibration(k = k)
    rows = []
    for spec in (C.US['esc']['spec'], C.US['esc']['comparisonSpec']):
        for ρ in C.US['esc']['ρTable']:
            if spec != C.US['esc']['spec'] and (float(ρ), spec) not in calib:
                continue                                # the comparison cost is calibrated under log only
            θStar = float(calib[(float(ρ), spec)]['θStar'])
            rows += _pathRows(D.escPath(ρ, spec = spec, k = k), 'main', cost, 'US', ('US', 'own'), spec, ρ, ξ0, θStar,
                              years)
    xi = D.escXi(k = k)                     # required: page C3's prose states the xi paths under both costs
    # A stale row (crossCuttingFindings #13): the workbook-xi path must be escPath.csv's rho = 1 path.
    own = xi[np.isclose(xi['xi'], ξ0) & (xi['kind'] == 'path')].sort_values('pos')
    if own.empty:
        raise D.MissingInput('{}/escXiRobustness.csv (no path at the workbook ξ = {})'
                             .format(os.path.basename(D.escDir(k)), ξ0))
    base = D.escPath(1.0, k = k)
    for c in ('θ', 'τ', 'τExo'):
        if len(own) != len(base) or not np.all(np.abs(own[c].values - base[c].values) <= 1e-8):
            raise ValueError('{}/escXiRobustness.csv at the workbook ξ = {} does not reproduce escPath.csv ({}); '
                             're-run python/US/runESCxi.py'.format(os.path.basename(D.escDir(k)), ξ0, c))
    for ξ in sorted(xi['xi'].unique()):
        sub = xi[np.isclose(xi['xi'], ξ)]
        θStar = float(sub[sub['kind'] == 'calibration'].iloc[0]['θStar'])
        rows += _pathRows(sub[sub['kind'] == 'path'], 'xi', cost, 'US', ('US', 'own'), C.US['esc']['spec'], 1.0,
                          round(float(ξ), 2), θStar, years)
    host = D.escPathHost(k = k)
    for (country, w), sub in host.groupby(['country', 'wedgeFrom']):
        rows += _pathRows(sub, 'host', cost, country, (w,), C.US['esc']['spec'], 1.0, ξ0, float(sub['θStar'].iloc[0]),
                          years)
    return rows


def escPaths():
    """ ESC_Paths: the baseline path with the design chosen and with it pinned at the observed design. The
    U.S. at the three rho of the paper under the cost and under the comparison cost where it is calibrated
    (rho = 1; 'main') and at three Frisch elasticities, rho = 1 ('xi', escXiRobustness.csv, required as in
    escCounterfactuals); the UK and France at rho = 1 under the U.S. cost parameter and under their own where they have one
    ('host'; escPathHost.csv). Each under both of section 7's costs (`cost`, ESCCOSTS). The U.S. rows go under both
    wedgeFrom labels, its own parameter being the U.S. one. Under CRRA the pinned workweek is not stored
    (escPathCRRA.csv): NaN. """
    years, ξ0 = _pathYears(), _workbookXi()
    return pd.DataFrame(sum((_escPathsAt(k, years, ξ0) for k in ESCCOSTS), []))

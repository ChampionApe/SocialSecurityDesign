r""" Workbook loaders for the France and UK variants -- the ModelFR counterpart of test.py.

    import testEU
    m = testEU.model('FR')                 # France
    m = testEU.model('UK')                 # UK, its own income groups

Same sheet layout as data/USMain_test.xlsx with two differences, both of which are the calibration
protocol showing up in the data (see modelFR.py):

  * There is NO '30y interest' column. R is not a target for these countries -- beta is imposed from the
    US -- so the number was never collected. db['R0'] is set to NaN rather than left at ModelUS's default
    2.443: nothing in ModelFR reads it, and a stray ModelUS.calibrate on one of these workbooks must fail
    loudly rather than quietly target the US interest rate.
  * The population sheet's 'Worker-to-retiree' is the column the model uses. France's is adjusted for its
    lower retirement age (1.5 years below the US in 1990 and 2020, 1 from 2050; PaG 2021, appendix F) and so differs from
    the raw census ratio in the neighbouring column; the UK's does not need the adjustment. Under config.US['retAgeFR']
    'exit' France's is instead its WPP ratio corrected for the gap in effective exit ages, and under
    config.US['retAgeUK'] 'exit' the UK's (hostNu).

France's 'Replacement rate' is 1, and getTheta (getθ) returns exactly 1 at RR0 = 1 for ANY income
grouping -- which is why France's groups can be cut at US percentiles. Nothing imposes theta = 1 in code. Under a θ
variant `model(..., θ0 = ...)` imposes θ as db['θ0'] instead (ModelUS.getθ; oecd26.countryθ0), which does not
depend on the grouping either, so France stays at the US cuts with its ω recalibrated there.

Only the 'population', 'heterogeneity' and 'calibration' sheets are read. The regrouped pairs in the
workbooks (the UK at US percentiles, 'heterogeneityUS'/'calibrationUS'; France at the UK's cuts,
'heterogeneityUK'/'calibrationUK') are not (data/README.md).
"""
import os
import numpy as np, pandas as pd

DIR = os.path.dirname(os.path.abspath(__file__))
PATHS = {'FR': os.path.join(DIR, '..', '..', 'data', 'FRMain.xlsx'),
         'UK': os.path.join(DIR, '..', '..', 'data', 'UKMain.xlsx')}
PAPERCONFIG = os.path.join(DIR, '..', 'paper', 'config.py')
t_ss = 3   # number of periods in steady state, as in test.py
RETAGESFR = ('nra', 'exit')
HOSTKEY = {'FR': 'retAgeFR', 'UK': 'retAgeUK'}       # config.US's key of each country's own ν correction
HOSTISO = {'FR': 'FRA', 'UK': 'GBR'}
HOSTWPP = {'FR': 'Relative cohort sizes', 'UK': 'Worker-to-retiree'}   # the population sheet's WPP ratio


def retAgeHost(country):
    """ python/paper/config.py's US['retAgeFR'] or US['retAgeUK'] ('nra' where absent): which correction `country`'s
    own ν carries. """
    import importlib.util
    sp = importlib.util.spec_from_file_location('_paperConfigHost', PAPERCONFIG)
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    r = mod.US.get(HOSTKEY[country], 'nra')
    if r not in RETAGESFR:
        raise ValueError("config.US[{!r}] must be one of {}, got {!r}".format(HOSTKEY[country], RETAGESFR, r))
    return r


def hostNu(country, dft, retAge):
    """ `country`'s ν at the population sheet's dates (dft, indexed on t) under `retAge`: 'nra' the sheet's
    'Worker-to-retiree', bitwise (France's appendix-F column; the UK's WPP ratio, its normal retirement ages being the
    U.S. ones); 'exit' the sheet's WPP ratio (HOSTWPP) corrected by oecd26.retAgeNu at the country's gap to the U.S.
    men's effective exit age (one gap at every date, oecd26.retAgeGaps), the cross-section's row under retAge 'exit'. """
    if retAge == 'nra':
        return dft['Worker-to-retiree'].values.astype(float)
    import oecd26                                   # lazily: oecd26 imports the U.S. workbook and the model
    g = float(oecd26.retAgeGaps(HOSTISO[country], 'exit')[0])
    return oecd26.retAgeNu(dft[HOSTWPP[country]].values.astype(float), np.full(len(dft), g))


def load(country, t_ss = t_ss):
    """ Read a workbook into (pars, kwargs, dates, workweek), in the form ModelFR's __init__ wants. """
    wb = pd.read_excel(PATHS[country], sheet_name = None, header = None)

    def sheet(name):
        d = pd.DataFrame(wb[name].values)
        return pd.DataFrame(d.iloc[1:, ].values, columns = d.iloc[0, :])

    dft = sheet('population').set_index('t')
    datesLog = dft.index
    νLog = hostNu(country, dft, retAgeHost(country))
    TLog = len(datesLog)
    dates = datesLog.union(pd.Index([datesLog[-1]+30*i for i in range(1, t_ss+1)]))

    dfj = sheet('heterogeneity')
    dfc = sheet('calibration')
    dfc = pd.Series(dfc.iloc[0, :].values, index = dfc.columns)
    workweek = dfc['Average workweek']

    kwargs = {'T': TLog+t_ss, 'nj': dfj.shape[0]}
    pars = {'α': dfc['Capital income share'], 'ξ': dfc['Labor supply elasticity'],
            'ν': np.hstack([νLog, np.full(kwargs['T']-TLog, νLog[-1])]),
            'τ0': dfc['Pension tax'], 'RR0': dfc['Replacement rate'],
            'RRGroups': tuple(eval(str(dfc['Replacement rate groups']))),
            't0': dates.get_loc(dfc['Calibration year']),
            'γj': dfj['Population shares'].values.astype(float),
            'μj': dfj['Voting shares'].values.astype(float),
            'Xj': 1,
            'zxj': (dfj['Hours'].values/(dfj['Hours'].values.mean())).astype(float),
            'zηj': (dfj['Income'].values/(dfj['Income'].values.mean())).astype(float)}
    pars['h0'] = workweek / (7*12)      # share of time spent working -- the SAME unit as the US's h0,
                                        # which is what makes the h̄ ratio in ModelFR.hbarTarget meaningful
    pars['Ushare0'] = dfc['Universal share']
    pars['R0'] = np.nan                 # not a target here, and not in the workbook -- see the docstring
    return pars, kwargs, dates, workweek


def usReference(commonX = False, h0US = None, csv = None):
    """ ModelFR's US reference as a callable of rho, read off the US rho sweep. h0US defaults to the US
    workbook's own average-hours input, loaded from test.py. csv: the sweep's path instead of
    results/calibration/US_rhoGrid[CommonX].csv (a θ variant's own, calibrateRhoGrid.py --theta). """
    from modelFR import usRefFromCsv
    if h0US is None:
        import test as testmod
        h0US = testmod.pars['h0']
    if csv is None:
        csv = os.path.join(DIR, '..', '..', 'results', 'calibration',
                           'US_rhoGrid' + ('CommonX' if commonX else '') + '.csv')
    return usRefFromCsv(csv, h0US = h0US)


def model(country, commonX = False, ρ = 1., ω = 2., usCsv = None, θ0 = None, **over):
    """ A ready ModelFR on `country`'s workbook, with the US reference wired to the sweep of the matching
    variant (or the sweep at usCsv; see usReference). beta is NOT set here: setUSRef
    installs it from the reference at db's rho when calibrate runs. θ0: the design imposed as db['θ0']
    (ModelUS.getθ); None leaves it to the workbook's RR0.
    """
    from modelFR import ModelFR
    pars, kwargs, dates, workweek = load(country)
    pars.update({'ρ': ρ, 'ω': ω} | ({} if θ0 is None else {'θ0': float(θ0)}))
    m = ModelFR(pars = pars | over, commonX = commonX,
                usRef = usReference(commonX = commonX, csv = usCsv), **kwargs)
    m.db['dates'], m.db['workweek'], m.db['country'] = dates, workweek, country
    return m

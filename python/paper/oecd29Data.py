r""" Stage (0) inputs for calibrating the log model to the 29 OECD members before 2000 (figure 1's sample):
workers per retiree nu_t by model period, the paper's three income groups, and secondary inputs. Not wired
into build.py, REPLICATION.md or the tests.

Run:  .venv\Scripts\python.exe -I python\paper\oecd29Data.py [--force] [--refetch]
      -> data/oecd29.csv, data/oecd29_population.csv, data/oecd29_distribution.csv, data/oecd29_sources.csv

`main` is the only code here that touches the network: `fetch` downloads every raw file missing from
logs/oecd29_raw/<source>/ (gitignored), then `collect` reads them. main skips when the four outputs exist
unless --force; raw files are reused unless --refetch. Raw files are untrusted data: never executed, zip
members read in memory, nothing written into their directories. A raw file's mtime is its retrieval date.

nu_t = P(t)/P(t-30): UN WPP 2024 total population on 1 July, medium variant. It reproduces the workbooks'
'Worker-to-retiree' (US, UK) and France's 'Relative cohort sizes' to 5e-7 at 1990-2080. 1960 repeats 1990
and 2110 repeats 2080, as in the workbooks.

Income groups follow appendix F: 'own' cuts the bottom group where its mean is half the population mean and
the next where its mean equals it; 'usPct' cuts at the US cumulative shares 0.582266/0.757993. y_i is the
group mean over the population mean; Vtilde = sum_i gamma_i (y_i - 1)^2 / y_i. A distribution is a
piecewise-linear Lorenz curve (income constant within a reported bin) or a two-piece lognormal fitted to
quantile ratios; 'own' is NaN when no such group exists at the data's resolution (never filled).

Secondary inputs: CSES turnout by household-income quintile, with mu_i its average over each group's
percentile range, flat within quintiles (FRMain.xlsx's and UKMain.xlsx's voting shares are this rule on the
file's own elections, France at the US cuts and the UK at its workbook's cuts; data/README.md); PaG 2021 normal retirement ages, also behind the nu variant
wpp2024_nraAdj; PWT 11.0; OECD hours worked; ILO SDG 8.3.1 (informality) and 1.3.1 (pension coverage).
"""
import os, io, re, sys, zipfile, argparse, datetime
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import oecdFigure1 as F1                     # COUNTRIES, the User-Agent, _get, _RemoteFile, the WID url

ISO3, NAME = F1.ISO3, {c[1]: c[0] for c in F1.COUNTRIES}
C29 = '+'.join(ISO3)
RAWDIR = os.path.join(C.REPO, 'logs', 'oecd29_raw')
OUT = {k: os.path.join(C.DATA, 'oecd29{}.csv'.format(s)) for k, s in
       (('main', ''), ('pop', '_population'), ('dist', '_distribution'), ('src', '_sources'))}
FIG1 = os.path.join(C.DATA, 'oecdFigure1.csv')
WORKBOOKS = {'USA': ('USMain_test.xlsx', 'heterogeneity'), 'GBR': ('UKMain.xlsx', 'heterogeneity'),
             'FRA': ('FRMain.xlsx', 'heterogeneity'), 'FRA@UK': ('FRMain.xlsx', 'heterogeneityUK')}

PERIODS = [1990, 2020, 2050, 2080]           # nu_t = P(t)/P(t-30); the workbooks copy 1990 to 1960, 2080 to 2110
FR_GAPS = {1990: 2.5, 2020: 2.5, 2050: 2.0, 2080: 2.0}          # appendix F until 2026-10-26: France 2.5 years before
                                                                # the US/UK, 2 from 2050 (now the nraAdj gaps)
FR_WORKBOOK = {1990: 1.026299, 2020: 0.966884, 2050: 0.903744, 2080: 0.870591}   # FRMain.xlsx population!D until then
DIST_YEAR, DIST_WINDOW = 2019, 2             # last pre-pandemic year; ESP, IRL, ITA earnings ratios jump in 2020-21
SURVEY = {'USA': 2018, 'GBR': 2006, 'FRA': 2018}   # CPS March 2019, BHPS wave 17 (2006-07), LIS 2018
US_CUTS = (0.582266, 0.757993)               # data/USMain_test.xlsx heterogeneity, cumulative shares
UK_CUTS = (0.291632, 0.872847)               # data/UKMain.xlsx heterogeneity, cumulative shares
Z9 = norm.ppf(0.9)

WPP = 'https://population.un.org/wpp/assets/Excel%20Files/'
ILO = 'https://rplumber.ilo.org/'
SDMX = F1.SDMX
# key: (raw subdirectory, file name, url, zip member or None). The archive zips of WPP 2019/2022 (0.3/3.5 GB)
# are read member by member through HTTP range requests, as oecdFigure1 does for WID.
RAW = {
    'wpp2024':   ('wpp2024', 'WPP2024_TotalPopulationBySex.csv.gz',
                  WPP + '1_Indicator%20(Standard)/CSV_FILES/WPP2024_TotalPopulationBySex.csv.gz', None),
    'wpp2024up': ('wpp2024update', 'WPP2024_CSV_files_update.zip',
                  WPP + '1_Indicator%20(Standard)/WPP2024_CSV_files_update.zip', None),
    'wpp2022':   ('wpp2022', 'WPP2022_TotalPopulationBySex.zip', WPP + '5_Archive/WPP2022-CSV-data.zip',
                  'CSV_FILES/WPP2022_TotalPopulationBySex.zip'),
    'wpp2019':   ('wpp2019', 'WPP2019_TotalPopulationBySex.csv', WPP + '5_Archive/WPP2019-CSV-data.zip',
                  'WPP2019_TotalPopulationBySex.csv'),
    'decI':      ('oecd_earn', 'dec_i_29.csv', SDMX + 'OECD.ELS.SAE,DSD_EARNINGS@DEC_I,1.0/' + C29
                  + '.......?startPeriod=2004&endPeriod=2024&format=csvfilewithlabels', None),
    'iloLid':    ('ilo_lid', 'lid_29.csv', ILO + 'data/indicator/?id=LAP_2LID_QTL_RT_A&ref_area=' + C29
                  + '&timefrom=2004&timeto=2025&format=.csv', None),
    'iloToc':    ('ilo_catalog', 'toc_indicator.csv', ILO + 'metadata/toc/indicator/?lang=en&format=.csv', None),
    'eag':       ('oecd_eag', 'eag_distr_29.csv', SDMX + 'OECD.EDU.IMEP,DSD_EAG_LSO_EA@DF_LSO_EARN_DISTR_MEDIAN,1.0/'
                  + C29 + '.M+_T.Y25T64._T.............?startPeriod=2005&endPeriod=2024&format=csvfilewithlabels', None),
    'idd':       ('oecd_idd2', 'idd_USA_FRA_GBR.csv', SDMX + 'OECD.WISE.INE,DSD_WISE_IDD@DF_IDD,1.0/USA+FRA+GBR.A.'
                  'D5_1_INC_DISP+D9_5_INC_DISP......?startPeriod=2004&endPeriod=2021&format=csvfilewithlabels', None),
    'silc':      ('eurostat_silc2', 'ilc_di01_FR_UK.csv', 'https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/'
                  'data/ilc_di01/A..SHARE.EUR.FR+UK?format=SDMX-CSV&startPeriod=2005&endPeriod=2020', None),
    **{'wid' + c2: ('wid', 'WID_data_{}.csv'.format(c2), F1.URL['wid'], 'WID_data_{}.csv'.format(c2))
       for c2 in ('US', 'FR', 'GB')},
    # secondary inputs. CSES: plain HTTP, no login; module 5 (2023-07-25), module 6 second advance release
    # (2025-12-16), the Integrated Module Dataset of modules 1-4 (2024-02-27; 655 MB csv in the zip).
    'csesM5':    ('cses_m5', 'cses5_csv.zip', 'https://cses.org/wp-content/uploads/2023/07/cses5_csv.zip', None),
    'csesM6':    ('cses_m6', 'cses6_csv.zip', 'https://cses.org/wp-content/uploads/2025/12/cses6_csv.zip', None),
    'csesImd':   ('cses_imd', 'cses_imd_csv.zip', 'https://cses.org/wp-content/uploads/2024/02/cses_imd_csv.zip', None),
    'pagNra':    ('pag2021_ages', 'statlinks/2cvxya.xlsx', 'https://stat.link/2cvxya', None),   # PaG 2021 figure 3.7
    'pwt':       ('pwt110', 'pwt110.xlsx', 'https://dataverse.nl/api/access/datafile/554105', None),
    'hours':     ('oecd_hours', 'AVG_ANN_HRS_WKD_oecd29_total.csv', SDMX + 'OECD.ELS.SAE,DSD_HW@DF_AVG_ANN_HRS_WKD,1.0/'
                  + C29 + '.HW.H_Y_PS._Z._Z.EMP.A.ACTUAL._T._Z.MEAN._Z._T?startPeriod=2015&endPeriod=2022'
                  '&format=csvfilewithlabels', None),
    'ilo831':    ('ilostat', 'SDG_0831_SEX_ECO_RT_A.csv', ILO + 'data/indicator/?id=SDG_0831_SEX_ECO_RT_A&ref_area='
                  + C29 + '&timefrom=2010&timeto=2025&format=.csv', None),
    'ilo131':    ('ilostat', 'SDG_0131_SEX_SOC_RT_A.csv', ILO + 'data/indicator/?id=SDG_0131_SEX_SOC_RT_A&ref_area='
                  + C29 + '&timefrom=2010&timeto=2025&format=.csv', None)}


def _path(key):
    d, f = RAW[key][:2]
    return os.path.join(RAWDIR, d, *f.split('/'))


def _retrieved(key):
    p = _path(key)
    return datetime.date.fromtimestamp(os.path.getmtime(p)).isoformat() if os.path.exists(p) else ''


# --- fetching (network; called from main only) -----------------------------------------------------------
def fetch(refetch = False):
    """ Download every raw file that is missing (all of them with refetch). """
    for key, (d, f, url, member) in RAW.items():
        p = _path(key)
        if os.path.exists(p) and not refetch:
            continue
        os.makedirs(os.path.dirname(p), exist_ok = True)
        print('  fetching', os.path.relpath(p, C.REPO), flush = True)
        if member is None:
            data = F1._get(url)
            with open(p, 'wb') as fh:
                fh.write(data)
        else:
            z = zipfile.ZipFile(io.BufferedReader(F1._RemoteFile(url), buffer_size = 1 << 20))
            with z.open(member) as src, open(p, 'wb') as dst:
                while (chunk := src.read(1 << 22)):
                    dst.write(chunk)


# --- population: nu_t ---------------------------------------------------------------------------------------
def wpp(key, locid = None):
    """ Total population (thousands, 1 July, medium variant), index year, columns the 29 ISO3. WPP 2019 has no
    ISO3 column: its LocID (ISO 3166 numeric) maps through `locid`, read off WPP 2024. """
    p = _path(key)
    if key == 'wpp2022':
        z = zipfile.ZipFile(p)
        d = pd.read_csv(z.open([m for m in z.namelist() if m.endswith('.csv')][0]), low_memory = False,
                        usecols = ['ISO3_code', 'Variant', 'Time', 'PopTotal'])
    elif key == 'wpp2019':
        d = pd.read_csv(p, usecols = ['LocID', 'Variant', 'Time', 'PopTotal'])
        d['ISO3_code'] = d['LocID'].map(locid)
    else:
        d = pd.read_csv(p, low_memory = False, usecols = ['ISO3_code', 'LocID', 'Variant', 'Time', 'PopTotal'])
    d = d[(d['Variant'] == 'Medium') & d['ISO3_code'].isin(ISO3)]
    pop = d.pivot(index = 'Time', columns = 'ISO3_code', values = 'PopTotal')[ISO3]
    if key == 'wpp2024':
        ids = d.drop_duplicates('ISO3_code').set_index('LocID')['ISO3_code']
        return pop, ids.to_dict()
    return pop


def wppUpdateCountries():
    """ Countries revised by the WPP 2024 CSV update (January 2026); none of the 29 may be among them. """
    z = zipfile.ZipFile(_path('wpp2024up'))
    d = pd.read_csv(z.open([m for m in z.namelist() if 'TotalPopulation' in m][0]), low_memory = False)
    return sorted(set(d['ISO3_code'].dropna()))


def nu(pop):
    """ nu_t = P(t)/P(t-30) at PERIODS, index period, columns ISO3. """
    return pd.DataFrame({t: pop.loc[t]/pop.loc[t - 30] for t in PERIODS}).T


def population(nra):
    """ Long table country x period x definition, plus the WPP 2024 nu table and the update check. `nra`:
    pagNra(), for the retirement-age adjustment of every country relative to the US. """
    p24, locid = wpp('wpp2024')
    pops = {'wpp2024': p24, 'wpp2022': wpp('wpp2022'), 'wpp2019': wpp('wpp2019', locid)}
    rows = []
    for dfn, pop in pops.items():
        n = nu(pop)
        for c in ISO3:
            for t in PERIODS:
                rows.append({'iso3': c, 'country': NAME[c], 'definition': dfn, 'period': t, 'nu': n.at[t, c],
                             'pop': pop.at[t, c], 'pop_lag30': pop.at[t - 30, c], 'note': ''})
            if dfn == 'wpp2024':
                for t, src in ((1960, 1990), (2110, 2080)):
                    rows.append({'iso3': c, 'country': NAME[c], 'definition': dfn, 'period': t,
                                 'nu': n.at[src, c], 'pop': np.nan, 'pop_lag30': np.nan,
                                 'note': 'copy of {} (workbook convention)'.format(src)})
    n24 = nu(p24)
    for t, v in nraRule(n24['FRA'], FR_GAPS).items():
        rows.append({'iso3': 'FRA', 'country': NAME['FRA'], 'definition': 'wpp2024_frRule', 'period': t, 'nu': v,
                     'pop': p24.at[t, 'FRA'], 'pop_lag30': p24.at[t - 30, 'FRA'],
                     'note': 'appendix F rule until 2026-10-26, gap {} years; FRMain.xlsx population!D had {:.6f} '
                             '(wpp2024_nraAdj since)'.format(
                         FR_GAPS[t], FR_WORKBOOK[t])})
    us = nra.loc['USA']
    for c in ISO3:
        gaps = {t: us[col] - nra.at[c, col] for t, col in
                ((t, 'nra_current_men' if t <= 2020 else 'nra_future_men') for t in PERIODS)}
        for t, v in nraRule(n24[c], gaps).items():
            rows.append({'iso3': c, 'country': NAME[c], 'definition': 'wpp2024_nraAdj', 'period': t, 'nu': v,
                         'pop': p24.at[t, c], 'pop_lag30': p24.at[t - 30, c],
                         'note': "men's normal retirement age {} years below the US's (PaG 2021, {})".format(
                             round(gaps[t], 2), 'current' if t <= 2020 else 'future')})
    fig1 = pd.read_csv(FIG1, index_col = 'iso3')
    for c in ISO3:
        rows.append({'iso3': c, 'country': NAME[c], 'definition': 'oecdHist_fig1', 'period': 2020,
                     'nu': fig1.at[c, 'popGrowth'], 'pop': fig1.at[c, 'pop2020_oecd']/1e3,
                     'pop_lag30': fig1.at[c, 'pop1990_oecd']/1e3, 'note': 'figure 1 popGrowth (data/oecdFigure1.csv)'})
    up = wppUpdateCountries()
    if set(up) & set(ISO3):
        raise ValueError('the WPP 2024 update revises {}: read it over the main file'.format(set(up) & set(ISO3)))
    return pd.DataFrame(rows), nu(p24), up


# --- income groups --------------------------------------------------------------------------------------------
def lorenz(edges, means):
    """ (P, L) at the bin edges, L(1) = 1, from population edges (K+1,) and bin means (K,). """
    edges, means = np.asarray(edges, float), np.asarray(means, float)
    cum = np.concatenate([[0.], np.cumsum(np.diff(edges)*means)])
    return edges, cum/cum[-1]


def _bisect(f, a, b, tol = 1e-12):
    fa = f(a)
    while b - a > tol:
        m = 0.5*(a + b)
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5*(a + b)


def cutsOwn(P, L, n = 200001):
    """ The paper's cuts. Grid search for the bracket, then bisection: L(p)/p - 1/2 rises in p, and
    L(p) - L(p1) - (p - p1) falls then rises to p1 - L(p1) > 0 at p = 1. (nan, nan) or (p1, nan) when the
    group does not exist, e.g. the bin just above p1 already averages more than the mean. """
    Lf = lambda p: np.interp(p, P, L)
    grid = np.linspace(0, 1, n)[1:]
    r = Lf(grid)/grid - 0.5
    k = int(np.argmax(r >= 0))
    if k == 0 or r[k] < 0:
        return np.nan, np.nan
    p1 = _bisect(lambda p: Lf(p)/p - 0.5, grid[k - 1], grid[k])
    g2 = grid[grid > p1]
    h = Lf(g2) - Lf(p1) - (g2 - p1)
    j0 = int(np.argmin(h))
    if h[j0] >= 0:
        return p1, np.nan
    j = j0 + int(np.argmax(h[j0:] >= 0))
    p2 = _bisect(lambda p: Lf(p) - Lf(p1) - (p - p1), g2[j - 1], g2[j])
    return p1, p2


def groups(P, L, cuts):
    """ (gamma (3,), y (3,), Vtilde) at the cumulative shares `cuts`; NaN throughout if a cut is NaN. """
    if np.isnan(cuts).any():
        return np.full(3, np.nan), np.full(3, np.nan), np.nan
    c = np.array([0., cuts[0], cuts[1], 1.])
    g = np.diff(c)
    y = np.diff(np.interp(c, P, L))/g
    return g, y, float(np.sum(g*(y - 1)**2/y))


def lognormal2(sL, sH, n = 200000):
    """ Two-piece lognormal, ln y = s z with s = sL below the median and sH above, on n equal bins. """
    z = norm.ppf((np.arange(n) + 0.5)/n)
    return lorenz(np.linspace(0, 1, n + 1), np.exp(np.where(z < 0, sL, sH)*z))


def decileRatioFit(r51, r95):
    """ sL, sH of the two-piece lognormal through P50/P10 and P90/P50. """
    return np.log(r51)/Z9, np.log(r95)/Z9


def eagFit(F):
    """ Two-piece lognormal (mu, sL, sH) through the cumulative shares F at 0.5, 1, 1.5, 2 times the median
    of all earners, least squares in normal scores; returns (sL, sH, max abs residual). """
    x, zF = np.log([0.5, 1, 1.5, 2]), norm.ppf(np.clip(F, 1e-6, 1 - 1e-6))
    res = least_squares(lambda t: (x - t[0])/np.where(x < t[0], t[1], t[2]) - zF, [0., 0.6, 0.6],
                        bounds = ([-2, 0.05, 0.05], [2, 3, 3]))
    return res.x[1], res.x[2], float(np.abs(res.fun).max())


def quintilesToGroups(q, cuts):
    """ Group averages of a step function of the percentile rank, flat within each quintile (q: 5 values in
    rank order), over the groups [0, c1), [c1, c2), [c2, 1]: FRMain.xlsx heterogeneityUK's rule. Assumes a
    voter's household-income rank is the rank of the group's earnings. """
    q = np.asarray(q, float)
    if np.isnan(cuts).any() or np.isnan(q).any():
        return np.full(3, np.nan)
    e = np.linspace(0, 1, 6)
    F = lambda p: np.interp(p, e, np.concatenate([[0.], np.cumsum(q/5)]))
    c = np.array([0., cuts[0], cuts[1], 1.])
    return np.diff(F(c))/np.diff(c)


def _year(years, target, window = None):
    """ The year nearest target (within window if given), ties to the earlier year; None if there is none. """
    ys = [int(y) for y in years if window is None or abs(int(y) - target) <= window]
    return min(ys, key = lambda t: (abs(t - target), t)) if ys else None


# The candidate sources. Each loader returns {iso3: {year: (P, L, fit-string)}} for the countries it has.
def decI(countries = ISO3):
    """ OECD earnings database, decile ratios of gross earnings, full-time employees, men. """
    d = pd.read_csv(_path('decI'))
    d = d[d['SEX'] == 'M'].pivot_table(index = ['REF_AREA', 'TIME_PERIOD'], columns = 'AGGREGATION_OPERATION',
                                       values = 'OBS_VALUE').dropna(subset = ['D5_1', 'D9_5'])
    out = {}
    for (c, y), r in d.iterrows():
        if c in countries:
            sL, sH = decileRatioFit(r['D5_1'], r['D9_5'])
            out.setdefault(c, {})[int(y)] = ('sigmaL={:.4f};sigmaH={:.4f};P50/P10={:.4f};P90/P50={:.4f}'.format(
                sL, sH, r['D5_1'], r['D9_5']), lambda sL = sL, sH = sH: lognormal2(sL, sH))
    return out


def _ilo():
    d = pd.read_csv(_path('iloLid'))
    d.columns = [x.strip('﻿').strip('"') for x in d.columns]
    return d


def iloLid(countries = ISO3):
    """ ILO modelled labour income distribution by decile (share of labour income, % of total). """
    d, out = _ilo(), {}
    for (c, y), g in d.groupby(['ref_area', 'time']):
        s = g.sort_values('classif1')['obs_value'].to_numpy(float)
        if c in countries and len(s) == 10:
            flag = ''.join(sorted(set(g['obs_status'].dropna())))
            out.setdefault(c, {})[int(y)] = ('deciles;obs_status={}'.format(flag or 'none'),
                                             lambda s = s: lorenz(np.linspace(0, 1, 11), s*10))
    return out


def wid(var):
    """ WID g-percentile averages (127 bins) of `var` for the US, France and the UK. """
    out = {}
    for c, c2 in (('USA', 'US'), ('FRA', 'FR'), ('GBR', 'GB')):
        d = pd.read_csv(_path('wid' + c2), sep = ';', usecols = ['variable', 'percentile', 'year', 'value'])
        d = d[d['variable'] == var]
        for y, g in d.groupby('year'):
            bins = []
            for p, v in zip(g['percentile'], g['value']):
                m = re.fullmatch(r'p([\d.]+)p([\d.]+)', p)
                lo, w = float(m[1]), float(m[2]) - float(m[1])
                step = 1 if lo < 99 else 0.1 if lo < 99.9 else 0.01 if lo < 99.99 else 0.001
                if abs(w - step) < 1e-9:
                    bins.append((lo, v))
            if len(bins) == 127:
                bins.sort()
                e = np.array([b[0] for b in bins] + [100.])/100
                out.setdefault(c, {})[int(y)] = ('g-percentiles(127)', lambda e = e, v = np.array([b[1] for b in bins]):
                                                 lorenz(e, v))
    return out


def eag(countries = ('USA', 'FRA', 'GBR')):
    """ OECD Education at a Glance, men 25-64 with earnings (full and part time), shares at <= 0.5, 1, 1.5,
    2 times the median of all earners; two-piece lognormal fit. """
    d = pd.read_csv(_path('eag'))
    d = d[(d['ATTAINMENT_LEV'] == '_T') & (d['SEX'] == 'M') & (d['AGE'] == 'Y25T64')]
    p = d.pivot_table(index = ['REF_AREA', 'TIME_PERIOD'], columns = 'INCOME', values = 'OBS_VALUE')/100
    out = {}
    for (c, y), r in p.iterrows():
        if c not in countries:
            continue
        F = np.cumsum([r['MILE_50'], r['MIGT_50_LE100'], r['MIGT_100_LE150'], r['MIGT_150_LE200']])
        if np.isnan(F).any():
            continue
        sL, sH, err = eagFit(F)
        out.setdefault(c, {})[int(y)] = ('sigmaL={:.4f};sigmaH={:.4f};maxResid={:.3f}'.format(sL, sH, err),
                                         lambda sL = sL, sH = sH: lognormal2(sL, sH))
    return out


def idd():
    """ OECD IDD, P50/P10 and P90/P50 of equivalised disposable income (total population), two-piece lognormal. """
    d = pd.read_csv(_path('idd'))
    d = d[(d['AGE'] == '_T') & (d['DEFINITION'] == 'D_CUR')]
    p = d.pivot_table(index = ['REF_AREA', 'TIME_PERIOD'], columns = 'MEASURE', values = 'OBS_VALUE')
    out = {}
    for (c, y), r in p.dropna().iterrows():
        sL, sH = decileRatioFit(r['D5_1_INC_DISP'], r['D9_5_INC_DISP'])
        out.setdefault(c, {})[int(y)] = ('sigmaL={:.4f};sigmaH={:.4f}'.format(sL, sH),
                                         lambda sL = sL, sH = sH: lognormal2(sL, sH))
    return out


def silc():
    """ Eurostat EU-SILC ilc_di01, decile shares of equivalised disposable income, France and the UK. """
    d = pd.read_csv(_path('silc'))
    d = d[(d['statinfo'] == 'SHARE') & (d['unit'] == 'EUR') & d['quant_inc'].str.fullmatch(r'D\d+')]
    out = {}
    for (geo, y), g in d.groupby(['geo', 'TIME_PERIOD']):
        s = g.set_index('quant_inc')['OBS_VALUE'].reindex(['D{}'.format(k) for k in range(1, 11)]).to_numpy(float)
        if not np.isnan(s).any():
            out.setdefault({'FR': 'FRA', 'UK': 'GBR'}[geo], {})[int(y)] = (
                'deciles', lambda s = s: lorenz(np.linspace(0, 1, 11), s*10))
    return out


SOURCES_DIST = {   # name: (concept, loader, role) -- 'data' sources are built for all 29 countries
    'oecd_earnings_men_ft': ('gross earnings of full-time employees, men (OECD earnings database, P10/P50/P90), '
                             'two-piece lognormal', decI, 'data'),
    'ilo_labour_income':    ('labour income of all employed persons, both sexes (ILO modelled estimates), deciles',
                             iloLid, 'data'),
    'oecd_eag_men':         ('earnings of employed men 25-64, full and part time (OECD Education at a Glance), '
                             'five bins, two-piece lognormal', eag, 'validation'),
    'wid_ptinc':            ('pre-tax national income, equal-split adults 20+ (WID ptinc j992)',
                             lambda: wid('aptincj992'), 'validation'),
    'wid_diinc':            ('post-tax disposable income, equal-split adults 20+ (WID diinc j992)',
                             lambda: wid('adiincj992'), 'validation'),
    'oecd_idd_disposable':  ('equivalised household disposable income (OECD IDD, P50/P10 and P90/P50), '
                             'two-piece lognormal', idd, 'validation'),
    'eusilc_disposable':    ('equivalised household disposable income (EU-SILC ilc_di01), deciles', silc, 'validation')}
RECOMMENDED = 'oecd_earnings_men_ft'


def _rows(c, role, source, concept, year, PL, fit, conventions, note = ''):
    out = []
    for conv, cuts in conventions:
        if cuts is None:
            cuts = cutsOwn(*PL)
        g, y, V = groups(*PL, cuts)
        edges = [0., cuts[0], cuts[1], 1.]
        for i in range(3):
            out.append({'iso3': c.split('@')[0], 'country': NAME[c.split('@')[0]], 'role': role, 'source': source,
                        'concept': concept, 'year': year, 'convention': conv, 'group': i + 1,
                        'cut_lo': edges[i], 'cut_hi': edges[i + 1], 'gamma': g[i], 'y': y[i], 'Vtilde': V,
                        'fit': fit, 'note': note})
    return out


def workbookTargets():
    """ The three calibrated groupings (and France at the UK's cuts), read from the workbooks. """
    rows = []
    label = {'USA': ('own', 'CPS March Supplement 2019'), 'GBR': ('own', 'BHPS wave 17 (2006-07)'),
             'FRA': ('usPct', 'LIS 2018'), 'FRA@UK': ('ukCuts', 'LIS 2018, at the UK cuts')}
    for c, (f, sheet) in WORKBOOKS.items():
        d = pd.read_excel(os.path.join(C.DATA, f), sheet_name = sheet, header = None)
        d = pd.DataFrame(d.iloc[1:].to_numpy(), columns = d.iloc[0].to_numpy())
        d = d[d['Group'] != 0]
        s, inc = d['Population shares'].to_numpy(float), d['Income'].to_numpy(float)
        P, L = lorenz(np.concatenate([[0.], np.cumsum(s/s.sum())]), inc)
        conv, survey = label[c]
        cuts = tuple(P[1:3])
        rows += _rows(c, 'target', 'workbook', survey + '; labour income, men 25-60', SURVEY[c.split('@')[0]],
                      (P, L), '{}!{}'.format(f, sheet), [(conv, cuts)])
    return rows


def distribution():
    """ The long distribution table: workbook targets, every candidate at the survey years of the three
    calibrated countries, and the 'data' sources for all 29 at DIST_YEAR (+-DIST_WINDOW). """
    rows = workbookTargets()
    for name, (concept, loader, role) in SOURCES_DIST.items():
        data = loader()
        for c in ('USA', 'GBR', 'FRA'):
            if c not in data:
                continue
            y = _year(data[c], SURVEY[c])
            fit, PL = data[c][y]
            conv = [('own', None), ('usPct', US_CUTS)] + ([('ukCuts', UK_CUTS)] if c == 'FRA' else [])
            rows += _rows(c, 'validation', name, concept, y, PL(), fit, conv,
                          '' if y == SURVEY[c] else 'nearest year to the survey ({})'.format(SURVEY[c]))
        if role == 'data':
            for c in ISO3:
                y = _year(data.get(c, {}), DIST_YEAR, DIST_WINDOW)
                if y is None:
                    continue
                fit, PL = data[c][y]
                rows += _rows(c, 'data', name, concept, y, PL(), fit, [('own', None), ('usPct', US_CUTS)])
    return pd.DataFrame(rows)


def validationTable(dist):
    """ Vtilde and group shares of every candidate under the convention of each workbook target. """
    t = dist[dist['role'] == 'target']
    v = dist[dist['role'].isin(['validation', 'target'])]
    out = []
    for (c, conv), g in t.groupby(['iso3', 'convention']):
        for src, h in v[(v['iso3'] == c) & (v['convention'] == conv)].groupby('source', sort = False):
            h = h.sort_values('group')
            out.append({'iso3': c, 'convention': conv, 'source': src, 'year': int(h['year'].iloc[0]),
                        'gamma': '/'.join('{:.3f}'.format(x) for x in h['gamma']),
                        'y': '/'.join('{:.3f}'.format(x) for x in h['y']), 'Vtilde': h['Vtilde'].iloc[0],
                        'dV': h['Vtilde'].iloc[0] - g['Vtilde'].iloc[0]})
    return pd.DataFrame(out)


# --- secondary inputs -----------------------------------------------------------------------------------------
# CSES harmonised names per file: study, ISO3, election year, sample and demographic weights, household-income
# quintile (1-5; 6-9 study-specific or missing), turnout in the main election (0/1; 'not registered' coded
# 93 in modules 5/6 and 9999993 in the IMD).
CSES = {'csesM5':  dict(member = 'cses5.csv', study = 'E1004', iso = 'E1006_UNALPHA3', year = 'E1008',
                        w_sample = 'E1010_1', w_demo = 'E1010_2', income = 'E2010', vote = 'E3012', notreg = 93),
        'csesM6':  dict(member = 'cses6.csv', study = 'F1004', iso = 'F1006_UNALPHA3', year = 'F1009',
                        w_sample = 'F1101_1', w_demo = 'F1101_2', income = 'F2010_1', vote = 'F3010', notreg = 93),
        'csesImd': dict(member = 'cses_imd.csv', study = 'IMD1004', iso = 'IMD1006_UNALPHA3', year = 'IMD1008_YEAR',
                        w_sample = 'IMD1010_1', w_demo = 'IMD1010_2', income = 'IMD2006', vote = 'IMD3001',
                        notreg = 9999993)}
CSES_WINDOW, CSES_MINQ = (2015, 2022), 30
IMD_MODS = ['IMD1008_MOD_{}'.format(i) for i in range(1, 6)]
FIELDS = ('study', 'iso', 'year', 'w_sample', 'w_demo', 'income', 'vote')


def _csesRead(key):
    """ One CSES file as harmonised respondent rows for the 29. The IMD's module-5 rows carry no income
    quintile (all coded 9), so they are dropped in favour of the module-5 file. """
    s = CSES[key]
    with zipfile.ZipFile(_path(key)) as z, z.open(s['member']) as f:
        d = pd.read_csv(f, usecols = [s[c] for c in FIELDS] + (IMD_MODS if key == 'csesImd' else []),
                        encoding = 'latin1', low_memory = False)
    d = d[d[s['iso']].isin(ISO3)].copy()
    if key == 'csesImd':
        d = d[d['IMD1008_MOD_5'] != 1]
        d['module'] = d[IMD_MODS].to_numpy().argmax(axis = 1) + 1
    else:
        d['module'] = {'csesM5': 5, 'csesM6': 6}[key]
    d = d.rename(columns = {s[c]: c for c in FIELDS})
    d['notreg'] = d['vote'] == s['notreg']
    return d[['module'] + list(FIELDS) + ['notreg']]


def _csesWeight(g):
    """ Per study: the demographic weight if it varies, else the sample weight if it varies, else none;
    normalised to mean 1 within the study. """
    for col, lab in (('w_demo', 'demographic'), ('w_sample', 'sample')):
        if g[col].nunique() > 1:
            w = g[col].astype(float)
            return w/w[np.isfinite(w) & (w > 0)].mean(), lab
    return pd.Series(1.0, index = g.index), 'none'


def cses():
    """ Turnout by household-income quintile, one election per country: among elections with >= CSES_MINQ
    valid respondents in every quintile, the one nearest 2020 inside CSES_WINDOW (ties to the earlier year,
    then the later module), else the nearest outside it (cses_inWindow False). A volunteered 'not registered'
    counts as not voting: turnout among adults, not registered voters. Studies of one election in one module
    are pooled (Belgium's Flanders and Wallonia). """
    d = pd.concat([_csesRead(k) for k in ('csesM5', 'csesM6', 'csesImd')], ignore_index = True)
    d.loc[d['notreg'], 'vote'] = 0
    rows = []
    for (c, y, mod), e in d.groupby(['iso', 'year', 'module']):
        w = [_csesWeight(g) for _, g in e.groupby('study')]
        e = e.assign(w = pd.concat([x[0] for x in w]))
        v = e[e['income'].between(1, 5) & e['vote'].isin([0, 1])]
        vw = v[np.isfinite(v['w']) & (v['w'] > 0)]
        if v.empty:
            continue
        r = {'iso3': c, 'cses_year': int(y), 'cses_module': int(mod), 'cses_study': '+'.join(sorted(e['study'].unique())),
             'cses_weight': '+'.join(sorted({x[1] for x in w})), 'cses_n': len(v), 'cses_nNotRegistered': int(e['notreg'].sum())}
        for q in range(1, 6):
            x, xw = v[v['income'] == q], vw[vw['income'] == q]
            r['nq{}'.format(q)] = len(x)
            r['turnout_q{}_unw'.format(q)] = x['vote'].mean() if len(x) else np.nan
            r['turnout_q{}'.format(q)] = np.average(xw['vote'], weights = xw['w']) if len(xw) else np.nan
        rows.append(r)
    el = pd.DataFrame(rows)
    el = el[el[['nq{}'.format(q) for q in range(1, 6)]].min(axis = 1) >= CSES_MINQ]
    keep = []
    for c, g in el.groupby('iso3'):
        inside = g[g['cses_year'].between(*CSES_WINDOW)]
        y = _year((inside if len(inside) else g)['cses_year'], 2020)
        r = g[g['cses_year'] == y].sort_values('cses_module').iloc[-1].copy()
        r['cses_inWindow'] = len(inside) > 0
        keep.append(r)
    out = pd.DataFrame(keep).set_index('iso3').drop(columns = ['nq{}'.format(q) for q in range(1, 6)])
    return out.reindex(ISO3)


PAG_NAMES = {**{n: c for n, c, _ in F1.COUNTRIES}, 'Korea': 'KOR'}


def pagNra():
    """ PaG 2021 figure 3.7 (sheet Data): normal retirement age after a full career from age 22, men and women;
    current = retiring in 2020, future = entering the labour market in 2020 (legislated, including links to
    life expectancy). Fractional years. """
    d = pd.read_excel(_path('pagNra'), sheet_name = 'Data', header = None)
    head = ' '.join(str(v) for v in d.iloc[3].tolist() + d.iloc[4].tolist())
    if 'Current (retiring 2020)' not in head or list(d.iloc[4, 1:7]) != ['Men', 'Women', 'Gap', 'Men', 'Women', 'Gap']:
        raise ValueError('PaG 2021 figure 3.7 layout changed: ' + head)
    d = d.iloc[5:, :6]
    d.columns = ['name', 'nra_current_men', 'nra_current_women', 'gap', 'nra_future_men', 'nra_future_women']
    d = d.assign(iso3 = d['name'].astype(str).str.strip().map(PAG_NAMES)).dropna(subset = ['iso3']).set_index('iso3')
    return d[['nra_current_men', 'nra_current_women', 'nra_future_men', 'nra_future_women']].astype(float).reindex(ISO3)


def pwt(years = (2019, 2020)):
    """ PWT 11.0: labsh (labour compensation over GDP), K/Y = rnna/rgdpna (constant national prices) and
    cn/cgdpo (current PPPs), avh (annual hours per person engaged). """
    d = pd.read_excel(_path('pwt'), sheet_name = 'Data',
                      usecols = ['countrycode', 'year', 'labsh', 'rnna', 'rgdpna', 'cn', 'cgdpo', 'avh'])
    d = d[d['countrycode'].isin(ISO3) & d['year'].isin(years)].set_index(['countrycode', 'year'])
    v = pd.DataFrame({'labsh': d['labsh'], 'ky_na': d['rnna']/d['rgdpna'], 'ky_pp': d['cn']/d['cgdpo'], 'avh': d['avh']})
    return pd.DataFrame({'pwt_{}_{}'.format(k, y): v[k].xs(y, level = 'year') for k in v for y in years}).reindex(ISO3)


def oecdHours(years = (2019, 2020)):
    """ OECD average annual hours actually worked per person in employment (all workers). """
    d = pd.read_csv(_path('hours'))
    d = d[(d['MEASURE'] == 'HW') & (d['WORKER_STATUS'] == '_T') & (d['HOURS_TYPE'] == 'ACTUAL')]
    p = d.pivot(index = 'REF_AREA', columns = 'TIME_PERIOD', values = 'OBS_VALUE').reindex(ISO3)
    return pd.DataFrame({'hours_{}'.format(y): p.get(y) for y in years}, index = ISO3).astype(float)


def iloSdg(key, classif, prefix, target = 2020, window = 2):
    """ ILOSTAT SDG indicator in percent: total, men and women at the year nearest target within +-window
    (chosen on the total), with that year, its obs_status (B break, U unreliable) and source id. """
    d = pd.read_csv(_path(key), encoding = 'utf-8-sig')
    d = d[d['classif1'] == classif]
    out = pd.DataFrame(index = ISO3, columns = [prefix + s for s in ('', '_men', '_women', '_year', '_status')],
                       dtype = object)
    for c in ISO3:
        t = d[(d['ref_area'] == c) & (d['sex'] == 'SEX_T')].dropna(subset = ['obs_value'])
        y = _year(t['time'], target, window)
        if y is None:
            continue
        for sex, s in (('SEX_T', ''), ('SEX_M', '_men'), ('SEX_F', '_women')):
            x = d[(d['ref_area'] == c) & (d['sex'] == sex) & (d['time'] == y)]
            out.at[c, prefix + s] = float(x['obs_value'].iloc[0]) if len(x) else np.nan
        st = t.loc[t['time'] == y, 'obs_status'].iloc[0]
        out.at[c, prefix + '_year'], out.at[c, prefix + '_status'] = y, st if isinstance(st, str) else ''
    for s in ('', '_men', '_women', '_year'):
        out[prefix + s] = pd.to_numeric(out[prefix + s])
    return out


def secondary():
    """ The secondary inputs, one row per country. """
    return pd.concat([cses(), pagNra(), pwt(), oecdHours(),
                      iloSdg('ilo831', 'ECO_SECTOR_TOTAL', 'informal'),
                      iloSdg('ilo131', 'SOC_CONTIG_PENSION', 'pensionCoverage')], axis = 1).reindex(ISO3)


def nraRule(nuSeries, gaps):
    """ Appendix F's rule for a country retiring g years before the reference, nu (30 - g)/(30 + nu g); for
    one retiring |g| years later, the mirror image: (30 nu + |g|)/(30 - |g|). """
    return {t: (v*(30 - g)/(30 + v*g) if g >= 0 else (30*v - g)/(30 + g)) for (t, v), g in
            ((tv, gaps[tv[0]]) for tv in nuSeries.items())}


# --- the outputs ----------------------------------------------------------------------------------------------
def wide(popLong, nu24, dist, extra = None):
    """ One row per country: nu (WPP 2024), the recommended grouping, figure 1's spending and replacement
    rates (copied, not refetched), and the secondary inputs in `extra`. """
    df = pd.DataFrame(index = pd.Index(ISO3, name = 'iso3'))
    df['country'] = [NAME[c] for c in ISO3]
    for t in PERIODS:
        df['nu_{}'.format(t)] = nu24.loc[t]
    fig1 = pd.read_csv(FIG1, index_col = 'iso3')
    df['popGrowth_fig1'] = fig1['popGrowth']
    d = dist[(dist['role'] == 'data') & (dist['source'] == RECOMMENDED)]
    for c in ISO3:
        for conv, sfx in (('own', ''), ('usPct', '_usPct')):
            h = d[(d['iso3'] == c) & (d['convention'] == conv)].sort_values('group')
            if h.empty:
                continue
            df.at[c, 'dist_year'] = h['year'].iloc[0]
            if conv == 'own':
                df.at[c, 'dist_fit'] = h['fit'].iloc[0]
                df.at[c, 'cut1'], df.at[c, 'cut2'] = h['cut_hi'].iloc[0], h['cut_hi'].iloc[1]
                for i in range(3):
                    df.at[c, 'gamma{}'.format(i + 1)] = h['gamma'].iloc[i]
            for i in range(3):
                df.at[c, 'y{}{}'.format(i + 1, sfx)] = h['y'].iloc[i]
            df.at[c, 'Vtilde' + sfx] = h['Vtilde'].iloc[0]
    for col in ('pensionSpending', 'pensionSpending_cashInKind', 'grr050_men', 'grr100_men', 'bbIndex'):
        df[col] = fig1[col]
    if extra is not None:
        df = df.join(extra)
        q = extra[['turnout_q{}'.format(k) for k in range(1, 6)]]
        for c in ISO3:
            for cuts, sfx in (((df.at[c, 'cut1'], df.at[c, 'cut2']), ''), (US_CUTS, '_usPct')):
                for i, m in enumerate(quintilesToGroups(q.loc[c], np.array(cuts, float))):
                    df.at[c, 'mu{}{}'.format(i + 1, sfx)] = m
    for col in [c for c in df.columns if c.endswith('_year')] + ['cses_module', 'cses_n', 'cses_nNotRegistered']:
        if col in df:
            df[col] = pd.to_numeric(df[col]).round().astype('Int64')
    return df


def _byYear(s):
    """ '2019: AUT, BEL; 2020: ...; none: LUX' for a series of years indexed by ISO3. """
    s = pd.to_numeric(s.reindex(ISO3))
    parts = ['{}: {}'.format(int(y), ', '.join(s.index[s == y])) for y in sorted(set(s.dropna()))]
    return '; '.join(parts + (['none: ' + ', '.join(s.index[s.isna()])] if s.isna().any() else []))


def _secondarySources(sec):
    toc = pd.read_csv(_path('iloToc')).set_index('id')['last.update']
    n = lambda col: '{} of 29'.format(int(sec[col].notna().sum()))
    cses = 'https://cses.org/wp-content/uploads/'
    return [
        ('cses_year, cses_module, cses_study, cses_weight, cses_n, cses_nNotRegistered, cses_inWindow, turnout_q1-5, '
         'turnout_q1-5_unw; mu1-3, mu1-3_usPct', 'CSES', 'Comparative Study of Electoral Systems: module 5 (full '
         'release 2023-07-25), module 6 (second advance release 2025-12-16), Integrated Module Dataset of modules 1-4 '
         '(2024-02-27)', 'turnout in the main election (E3012, F3010, IMD3001; a volunteered "not registered" counts '
         'as not voting) by household-income quintile (E2010, F2010_1, IMD2006); turnout_q* weighted (demographic '
         'weight, else sample weight, as cses_weight), *_unw unweighted; mu1-3 average the quintiles over each income '
         "group's percentile range at cut1/cut2, mu1-3_usPct at the US cuts",
         'see dataset', n('turnout_q1') + '; election year ' + _byYear(sec['cses_year']), _retrieved('csesM5'),
         cses + '2023/07/cses5_csv.zip; ' + cses + '2025/12/cses6_csv.zip; ' + cses + '2024/02/cses_imd_csv.zip',
         'One election per country, nearest 2020 within {}-{} with >= {} valid respondents per quintile, else the '
         'nearest outside (cses_inWindow False). mu takes a voter\'s household-income rank for the rank of the '
         "group's earnings. The U.S. workbook's voting shares are the CPS's; France's and the UK's are these "
         'elections (the UK at its own cuts); see data/README.md.'.format(*CSES_WINDOW, CSES_MINQ)),
        ('nra_current_men, nra_current_women, nra_future_men, nra_future_women; oecd29_population.csv definition '
         'wpp2024_nraAdj', 'OECD', 'Pensions at a Glance 2021, figure 3.7 (StatLink xlsx, sheet Data)',
         'normal retirement age after a full career from age 22: current = retiring in 2020, future = entering the '
         'labour market in 2020 (legislated, including links to life expectancy)', 'PaG 2021',
         n('nra_current_men') + '; rules of 2020', _retrieved('pagNra'), RAW['pagNra'][2],
         'wpp2024_nraAdj applies appendix F\'s rule to every country with the gap to the US men\'s age (current for '
         '1990 and 2020, future for 2050 and 2080), mirrored for later retirement. France here: 64.5 and 66 against '
         'the US 66 and 67 (gaps 1.5 and 1), as in appendix F and FRMain.xlsx population!D since 2026-10-26; until '
         'then a French age of 63.5 with no source (gaps 2.5 and 2, definition wpp2024_frRule).'),
        ('pwt_labsh_2019, pwt_labsh_2020, pwt_ky_na_*, pwt_ky_pp_*, pwt_avh_*', 'Groningen Growth and Development '
         'Centre', 'Penn World Table 11.0 (DataverseNL, doi:10.34894/FABVLR), sheet Data',
         'labsh; ky_na = rnna/rgdpna (constant national prices); ky_pp = cn/cgdpo (current PPPs); avh', 'PWT 11.0',
         '2019, 2020 (' + n('pwt_labsh_2020') + ')', _retrieved('pwt'), RAW['pwt'][2],
         'K/Y in annual output units. The paper sets alpha = 0.30 for the US, France and the UK.'),
        ('hours_2019, hours_2020', 'OECD', 'Average annual hours actually worked, '
         'OECD.ELS.SAE:DSD_HW@DF_AVG_ANN_HRS_WKD(1.0)', 'hours actually worked per year per person in employment, '
         'all workers (WORKER_STATUS _T)', 'Data Explorer', '2019, 2020 (' + n('hours_2020') + ')',
         _retrieved('hours'), RAW['hours'][2], 'Hours levels move no rate under log preferences.'),
        ('informal, informal_men, informal_women, informal_year, informal_status', 'ILO',
         'ILOSTAT, SDG 8.3.1 (SDG_0831_SEX_ECO_RT_A)', 'informal employment, % of total employment, all economic '
         'activities (ECO_SECTOR_TOTAL)', 'last.update ' + str(toc.get('SDG_0831_SEX_ECO_RT_A', '')),
         n('informal') + '; ' + _byYear(sec['informal_year']), _retrieved('ilo831'), RAW['ilo831'][2],
         'Nearest 2020 within 2018-2022. informal_status: B break in series, U unreliable.'),
        ('pensionCoverage, pensionCoverage_men, pensionCoverage_women, pensionCoverage_year, pensionCoverage_status',
         'ILO', 'ILOSTAT, SDG 1.3.1 (SDG_0131_SEX_SOC_RT_A)', 'old age (SOC_CONTIG_PENSION): % of the population above '
         'statutory pensionable age receiving an old-age pension; capped at 100 by the ILO',
         'last.update ' + str(toc.get('SDG_0131_SEX_SOC_RT_A', '')), n('pensionCoverage') + '; '
         + _byYear(sec['pensionCoverage_year']), _retrieved('ilo131'), RAW['ilo131'][2], 'Nearest 2020 within 2018-2022.')]


def sourcesTable(nu24, dist, upd, sec):
    """ One row per column (or group of columns sharing a source), as data/oecdFigure1_sources.csv. """
    d = dist[(dist['role'] == 'data') & (dist['group'] == 1) & (dist['convention'] == 'own')]
    yrs = lambda s: _yearsText(d[d['source'] == s].set_index('iso3')['year'])
    toc = pd.read_csv(_path('iloToc'))
    iloUpd = toc.loc[toc['id'] == 'LAP_2LID_QTL_RT_A', 'last.update'].iloc[0]
    rows = [
        ('nu_1990, nu_2020, nu_2050, nu_2080; oecd29_population.csv definition wpp2024', 'UN DESA Population Division',
         'World Population Prospects 2024, TotalPopulationBySex', 'PopTotal (thousands, 1 July), Variant Medium; '
         'nu_t = P(t)/P(t-30)', 'WPP 2024 (July 2024; file of 13 Dec 2024)', '1960-2080 (29 of 29)',
         _retrieved('wpp2024'), RAW['wpp2024'][2],
         'Reproduces data/USMain_test.xlsx and UKMain.xlsx population!C (Worker-to-retiree) and FRMain.xlsx '
         'population!C (Relative cohort sizes) to 5e-7. The January 2026 CSV update revises {} only.'.format(
             ', '.join(upd) or 'no country')),
        ('oecd29_population.csv definition wpp2022', 'UN DESA Population Division',
         'World Population Prospects 2022, TotalPopulationBySex (archive zip member)', 'as wpp2024', 'WPP 2022',
         '1960-2080 (29 of 29)', _retrieved('wpp2022'), RAW['wpp2022'][2] + ' :: ' + RAW['wpp2022'][3],
         'Vintage sensitivity; does not reproduce the workbooks.'),
        ('oecd29_population.csv definition wpp2019', 'UN DESA Population Division',
         'World Population Prospects 2019, TotalPopulationBySex (archive zip member)', 'as wpp2024; countries by LocID',
         'WPP 2019', '1960-2080 (29 of 29)', _retrieved('wpp2019'), RAW['wpp2019'][2] + ' :: ' + RAW['wpp2019'][3],
         'Vintage sensitivity; does not reproduce the workbooks.'),
        ('popGrowth_fig1; oecd29_population.csv definition oecdHist_fig1', 'OECD', 'see data/oecdFigure1_sources.csv',
         'popGrowth: OECD historical population 2020/1990', '', '2020 (29 of 29)', '', 'data/oecdFigure1.csv',
         'Copied from figure 1, not refetched. France differs most from nu_2020 (overseas departments added from 1991).'),
        ('dist_year, dist_fit, cut1, cut2, gamma1-3, y1-3, Vtilde, y1-3_usPct, Vtilde_usPct; '
         'oecd29_distribution.csv source oecd_earnings_men_ft', 'OECD',
         'Earnings database, Decile ratios of gross earnings, OECD.ELS.SAE:DSD_EARNINGS@DEC_I(1.0)',
         'men (SEX=M), full-time employees, D5_1 = P50/P10 and D9_5 = P90/P50; two-piece lognormal through them',
         'Data Explorer', yrs('oecd_earnings_men_ft'), _retrieved('decI'), RAW['decI'][2],
         'Recommended source (data/README.md). Validation: see the validation rows of oecd29_distribution.csv.'),
        ('oecd29_distribution.csv source ilo_labour_income', 'ILO', 'ILOSTAT, ILO modelled estimates (ILOEST)',
         'LAP_2LID_QTL_RT_A: labour income distribution by decile (%), all employed persons',
         'last.update ' + str(iloUpd), yrs('ilo_labour_income'), _retrieved('iloLid'), RAW['iloLid'][2],
         'Alternative. obs_status M = model-based extrapolation, I = imputation (KOR, NZL all years), given in fit.'),
        ('oecd29_distribution.csv source oecd_eag_men (validation only)', 'OECD',
         'Education at a Glance, OECD.EDU.IMEP:DSD_EAG_LSO_EA@DF_LSO_EARN_DISTR_MEDIAN(1.0)',
         'men 25-64 with earnings, all work-time arrangements; shares at <=0.5, 0.5-1, 1-1.5, 1.5-2, >2 x median',
         'Data Explorer', 'USA, GBR, FRA nearest the survey year', _retrieved('eag'), RAW['eag'][2],
         'Not used for the 29: no ISL, JPN only 2012, GRC only 2017-18, and level breaks (e.g. GBR 2019, FRA 2020, '
         'BEL 2021, HUN 2019, NOR 2016).'),
        ('oecd29_distribution.csv sources wid_ptinc, wid_diinc (validation only)', 'World Inequality Lab',
         'World Inequality Database, bulk download (members WID_data_US/FR/GB.csv)',
         'aptincj992 / adiincj992: g-percentile averages, equal-split adults 20+', 'bulk zip', 'survey years',
         _retrieved('widUS'), F1.URL['wid'], 'Household-income concepts; WID has labour income (pllin) for France only.'),
        ('oecd29_distribution.csv source oecd_idd_disposable (validation only)', 'OECD',
         'Income Distribution Database, OECD.WISE.INE:DSD_WISE_IDD@DF_IDD(1.0)',
         'D5_1_INC_DISP, D9_5_INC_DISP (total population, current definition); two-piece lognormal', 'Data Explorer',
         'survey years', _retrieved('idd'), RAW['idd'][2], 'The flow has no decile shares.'),
        ('oecd29_distribution.csv source eusilc_disposable (validation only)', 'Eurostat', 'EU-SILC, ilc_di01',
         'decile shares (SHARE, D1-D10) of equivalised disposable income', 'API', 'FRA 2018, GBR 2006',
         _retrieved('silc'), RAW['silc'][2], 'No US.'),
        ('oecd29_distribution.csv role target', 'this repository', 'data/USMain_test.xlsx, UKMain.xlsx, FRMain.xlsx',
         'heterogeneity sheets (FRMain heterogeneityUK for ukCuts): Population shares, Income', '', '', '',
         'data/', 'Read with pandas, never written.'),
        ('pensionSpending, pensionSpending_cashInKind, grr050_men, grr100_men, bbIndex', 'OECD',
         'see data/oecdFigure1_sources.csv', 'copied from data/oecdFigure1.csv', '', '2020 (29 of 29)', '',
         'data/oecdFigure1.csv', 'Not refetched. Complete for the 29.')]
    rows += _secondarySources(sec)
    return pd.DataFrame(rows, columns = ['column', 'provider', 'dataset', 'series', 'vintage', 'years', 'retrieved',
                                         'url', 'note'])


def _yearsText(s, target = DIST_YEAR):
    s = s.reindex(ISO3)
    parts = ['{} ({} of {})'.format(target, int((s == target).sum()), len(s))]
    for y in sorted(set(s.dropna()) - {target}):
        parts.append('{}: {}'.format(int(y), ', '.join(s.index[s == y])))
    if s.isna().any():
        parts.append('none in {}-{}: {}'.format(target - DIST_WINDOW, target + DIST_WINDOW, ', '.join(s.index[s.isna()])))
    return '; '.join(parts)


def collect():
    sec = secondary()
    popLong, nu24, upd = population(sec)
    dist = distribution()
    return popLong, nu24, upd, dist, sec


def main():
    p = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--force', action = 'store_true', help = 'rebuild the outputs even if they exist')
    p.add_argument('--refetch', action = 'store_true', help = 'download the raw files again')
    a = p.parse_args()
    if all(os.path.exists(f) for f in OUT.values()) and not a.force:
        print('up to date: data/oecd29*.csv (--force to rebuild)')
        return
    fetch(a.refetch)
    popLong, nu24, upd, dist, sec = collect()
    main_ = wide(popLong, nu24, dist, sec)
    src = sourcesTable(nu24, dist, upd, sec)
    popLong.to_csv(OUT['pop'], index = False, float_format = '%.8g')
    dist.to_csv(OUT['dist'], index = False, float_format = '%.8g')
    main_.to_csv(OUT['main'], float_format = '%.8g')
    src.to_csv(OUT['src'], index = False)
    pd.set_option('display.width', 200)
    print(validationTable(dist).to_string(index = False, float_format = '{:.3f}'.format))
    for f in OUT.values():
        print('wrote', os.path.relpath(f, C.REPO))


if __name__ == '__main__':
    main()

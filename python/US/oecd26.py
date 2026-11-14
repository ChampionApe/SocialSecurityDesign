r""" The model of section 7 (ModelESC, a proportional cost of Eq esc:AB or none, log preferences) built from one
OECD country's data: the loader of the cross-country test (python/US/runOECD26.py -> results/oecd26/).

    import oecd26
    inp = oecd26.countryInputs('GBR', 'csesFlat')     # the turnout is always named: 'csesFlat' (the paper's), 'cses', ...
    inp = oecd26.countryInputs('GBR', 'csesFlat', thetaSource = 'C')   # θ from data/oecd26_theta.csv's variant 'C'
    β, ωSeed, src = oecd26.usBeta()
    m = oecd26.build(inp, β, ω = 1.5, λ = 10.)        # a ModelOECD, LOG grid installed

Inputs per country (`countryInputs`):

  ν_t     data/oecd29.csv nu_1990..nu_2080 (UN WPP 2024, P(t)/P(t-30)), on the workbooks' calendar
          (`DATES`): 1960 repeats 1990, 2110 repeats 2080, the steady-state periods repeat the last (`nuPath`).
          retAge 'exit' | 'nra' (not 'none', the default): appendix F's ν(30 - g)/(30 + ν g) at each WPP date,
          its mirror image for g < 0 (`retAgeNu`), g = A_US - A_c (`retAgeGaps`): A the men's average effective age of labour-market exit
          (data/oecd29_exitAge.csv, PaG 2021 figure 6.14, 2020) at every date, or nra_current_men at 1990, 2020
          and nra_future_men at 2050, 2080 (= data/oecd29_population.csv's wpp2024_nraAdj). The U.S. gap is 0
          and leaves its ν bitwise, so the anchor does not depend on retAge.
  groups  three income groups at the U.S. percentiles (cumulative shares 0.582266/0.757993): γ_i and
          relative income y_i (Σ_i γ_i y_i = 1) of data/oecd29_distribution.csv, role 'data', source
          oecd_earnings_men_ft, convention usPct. y_i enters as z^η_i of the common-X calibration
          (model.initProductivity_commonX), so Base.hηRatio = y_i and Base.Vtilde (Eq esc:Vtilde) is the
          file's Vtilde_usPct.
  θ       thetaSource 'index' (the default): eq calibration:thetaSimple of writing/US, θ = 2 - 1/I, with I =
          grr100_men/grr050_men (figure 1.1's index, PaG 2021), clipped to [0, 1] (`thetaFromIndex`; its inverse
          I = 1/(2 - θ) is exact under the benefit formula and any proportional cost, archive/notes/reviewResponse_2026-10/theory_closedForm.md section 2). Not getθ on the
          group means (eq calibration:theta): ModelOECD.getθ re-derives the same value at every refresh, so no
          updateAuxPars can undo it (crossCuttingFindings #9).
          thetaSource a variant name: θ is the variant's row of data/oecd26_theta.csv (columns iso3, variant,
          theta, inSample, spending; `thetaFromFile`), the datum I = 1/(2 - θ), so that getθ returns θ to rounding,
          and θRaw = θ. A missing row, or θ NaN or outside [0, 1], raises NoTheta.
  τ0      spending over the labour share, spending/100/(1 - α): under 'index' figure 1.1's cash pension spending
          (pensionSpending); under a variant the file's `spending` column (`spendingFromFile`), which for the
          variants of SPENDINGASFIGURE must equal pensionSpending to 1e-9 (ValueError otherwise); B's is public plus
          mandatory and quasi-mandatory private spending.
  α, ξ    the U.S. workbook's (data/USMain_test.xlsx through test.py), and its workweek for h0 (hours levels
          move no rate under log preferences).
  μ_i     by `turnout`, always named by the caller: 'cses': the country's own CSES turnout averaged over the same
          groups (mu1-3_usPct; none for Luxembourg, which raises NoTurnout); 'cps' (not the paper's; runOECD26 refuses
          it into results/oecd26/): the U.S. workbook's CPS propensities to vote, in every country (the U.S. percentiles
          are the CPS groups' cuts). 'pooled' (appendix H's stylized reading, not in
          TURNOUT): the U.S. at its own CSES turnout; the countries of COMPULSORY flat, μ = (1, 1, 1); every other
          country at `pooledProfile`, the mean of μ_c/μ_c,1 over the voluntary non-U.S. countries of the θ source's
          sample with CSES turnout (`turnoutSample`), so Luxembourg enters. 'csesFlat' (not in TURNOUT): 'cses',
          but a COMPULSORY country without CSES (Luxembourg) flat; every other country, the U.S. included, has its
          'cses' inputs but for the label. Only the ratios of μ across groups
          matter: μ_j multiplies retirees and workers alike in the political objective W_t of def:PEE
          (writing/Paper/Sections/Model.tex), and the rows at μ and 0.9μ agree to the solver's tolerance.
  β       betaRule 'common': imposed, the U.S. calibration without the cost at ρ = 1 in the common-X variant
          (`usBeta`, results/calibration/US_rhoGridCommonX.csv), as appendix G.1 imposes it on the UK and France.
          Under a θ variant the sweep is the variant's own (`sweepPath`, calibrateRhoGrid.py --theta; under --main,
          `thetaRoot`, results/calibration/'s), the U.S. at θ0 = the variant's U.S. θ, the same number as the
          anchor's. France and the UK take their θ0 from the same file (`countryθ0`).
          betaRule 'commonR': calibrated so that R_{t0} is db['R0'], the anchor's 2020 return (build's R0).

ModelOECD calibrates, with θ exogenous at the datum and the cost installed, under betaRule 'common' ω alone
against τ_{t0} (model.py eq:calibration with R dropped, as modelFR.ModelFR); under 'commonR' (β, ω) against
(R_{t0}, τ_{t0}) (model.py eq:calibration), or β alone against R_{t0} when ω is held (calibrateω False). The cost
parameter (the 'p' slot: λ under 'size', the λ of λ·Π_t/R_t under 'sizePiR') is calibrated by
ModelESC.calibrateWedge to eq esc:calibration (the design in force in 2020 on the freely simulated path is the
datum).
`closedFormTax` is the log tax at θ ≡ 1, archive/notes/reviewResponse_2026-10/theory_closedForm.md eq (1.4): at θ = 1 the cost vanishes with
its τ-derivative, so it holds under any λ and any μ, β -> 0 in the terminal period.

Only the common-X calibration is defined: relative hours are not among the inputs.
"""
import os, json, hashlib
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(REPO, 'data')
FILES = {'main': os.path.join(DATA, 'oecd29.csv'), 'pop': os.path.join(DATA, 'oecd29_population.csv'),
         'dist': os.path.join(DATA, 'oecd29_distribution.csv'), 'sample': os.path.join(DATA, 'oecd26_sample.csv'),
         'theta': os.path.join(DATA, 'oecd26_theta.csv'), 'exit': os.path.join(DATA, 'oecd29_exitAge.csv')}


def shown(path):
    """ A file as the messages name it: relative to the repository, '/'-separated. The messages of NoTheta and
    NoTurnout are written into the results csvs, which must not carry a local path. """
    try:
        rel = os.path.relpath(path, REPO)
    except ValueError:                     # another drive
        return os.path.basename(path)
    return os.path.basename(path) if rel.startswith('..') else rel.replace(os.sep, '/')


USSWEEP = os.path.join(REPO, 'results', 'calibration', 'US_rhoGridCommonX.csv')
RESULTS = os.path.join(REPO, 'results')
PAPERCONFIG = os.path.join(REPO, 'python', 'paper', 'config.py')
INDEX = 'index'                           # thetaSource of figure 1.1's index; reserved, no variant of the file
SPENDINGASFIGURE = ('paper', 'A', 'A2pt', 'C')     # variants whose `spending` is figure 1.1's pensionSpending
SPENDINGTOL = 1e-9
ISO3 = {'US': 'USA', 'FR': 'FRA', 'UK': 'GBR'}     # the workbooks' countries in data/oecd26_theta.csv

import test as testmod                    # the U.S. workbook: α, ξ, the CPS turnout, the workweek
from modelESC import ModelESC

DISTSOURCE, DISTCONV = 'oecd_earnings_men_ft', 'usPct'
NUYEARS = (1990, 2020, 2050, 2080)        # the four WPP periods of data/oecd29.csv
T_SS = 3                                  # steady-state periods after 2110, as test.py
DATES = pd.Index([1960, 1990, 2020, 2050, 2080, 2110] + [2110 + 30*k for k in range(1, T_SS + 1)])
T0YEAR = 2020
GS = {'n': 101, 'smoothKnots': 4, 'interpKind': 'linear'}     # runESC.GS
TURNOUT = ('cps', 'cses')                 # the readings the loader's tests loop over; the paper's is 'csesFlat' (config.US['oecd26'])
TURNOUTS = TURNOUT + ('pooled', 'csesFlat')   # every reading countryInputs accepts (module docstring, μ_i)
RETAGES = ('none', 'exit', 'nra')         # the retirement-age correction of ν (module docstring, ν_t)
PERIOD = 30.                              # years per model period, the 30 of retAgeNu
MU = ['mu1_usPct', 'mu2_usPct', 'mu3_usPct']
# Enforced compulsory voting around 2020 (International IDEA, Compulsory Voting database): flat turnout under
# 'pooled'. Greece's compulsory voting is not enforced (voluntary here); Turkey lies outside every sample.
COMPULSORY = ('AUS', 'BEL', 'LUX')
BETARULES = ('common', 'commonR')


class NoTurnout(ValueError):
    """ The country has no CSES election in data/oecd29.csv (Luxembourg). """


class NoTheta(ValueError):
    """ The country has no θ in [0, 1] under the variant of data/oecd26_theta.csv (no row, or NaN or out of range). """


def readTheta(path = None):
    """ data/oecd26_theta.csv (FILES['theta'] unless `path`): iso3, variant, theta (float), inSample (bool). A
    duplicate (iso3, variant) or a variant named 'index' raises (#8). """
    path = FILES['theta'] if path is None else path
    if not os.path.exists(path):
        raise FileNotFoundError('{}: no θ variants (the stage (0) file is written separately)'.format(path))
    d = pd.read_csv(path, dtype = {'iso3': str, 'variant': str}, float_precision = 'round_trip')   # θ as written
    lack = {'iso3', 'variant', 'theta', 'inSample'} - set(d.columns)
    if lack:
        raise ValueError('{}: no column {}'.format(path, sorted(lack)))
    d['theta'] = pd.to_numeric(d['theta'], errors = 'coerce')
    d['inSample'] = d['inSample'].astype(str).str.strip().str.lower().isin(('true', '1', '1.0'))
    if (d['variant'] == INDEX).any():
        raise ValueError("{}: the variant name {!r} is reserved for figure 1.1's index".format(path, INDEX))
    dup = d[d.duplicated(['iso3', 'variant'], keep = False)]
    if len(dup):
        raise ValueError('{}: duplicate rows {}'.format(path, sorted(set(zip(dup['iso3'], dup['variant'])))))
    return d


def thetaVariants(path = None):
    """ The variant names of data/oecd26_theta.csv, in the file's order. """
    return list(dict.fromkeys(readTheta(path)['variant']))


def thetaFromFile(iso3, variant, path = None):
    """ θ of `iso3` under `variant` of data/oecd26_theta.csv. An unknown variant raises ValueError; a missing row,
    or θ NaN or outside [0, 1], NoTheta. """
    d = readTheta(path)
    if variant not in set(d['variant']):
        raise ValueError('{}: no variant {!r} (variants: {})'.format(shown(path or FILES['theta']), variant,
                                                                    list(dict.fromkeys(d['variant']))))
    r = d[(d['iso3'] == iso3) & (d['variant'] == variant)]
    if r.empty:
        raise NoTheta('{}: no row for {} under variant {!r}'.format(shown(path or FILES['theta']), iso3, variant))
    θ = float(r['theta'].iloc[0])
    if not 0. <= θ <= 1.:                 # NaN fails too
        raise NoTheta('{}: θ = {!r} for {} under variant {!r}, not in [0, 1]'.format(
            shown(path or FILES['theta']), θ, iso3, variant))
    return θ


def spendingFromFile(iso3, variant, path = None):
    """ Pension spending, % of GDP, of `iso3` under `variant` of data/oecd26_theta.csv (its `spending` column). A file
    without the column raises ValueError; a missing row or a NaN, NoTheta. """
    d = readTheta(path)
    if 'spending' not in d.columns:
        raise ValueError('{}: no column spending'.format(shown(path or FILES['theta'])))
    r = d[(d['iso3'] == iso3) & (d['variant'] == variant)]
    if r.empty or not np.isfinite(pd.to_numeric(r['spending'], errors = 'coerce').iloc[0]):
        raise NoTheta('{}: no spending for {} under variant {!r}'.format(shown(path or FILES['theta']), iso3, variant))
    return float(r['spending'].iloc[0])


def rowsAtTheta(df, variant, θ0):
    """ The rows of csv frame `df` at a run's θ: thetaSource `variant` and θ0 within 1e-12; without a variant
    (None) the rows with no thetaSource, or INDEX, the workbook's θ. A missing column reads as NaN. """
    ts = df['thetaSource'] if 'thetaSource' in df.columns else pd.Series(np.nan, index = df.index)
    if variant is None:
        return df[ts.isna() | (ts == INDEX)]
    θ = pd.to_numeric(df['θ0'], errors = 'coerce') if 'θ0' in df.columns else pd.Series(np.nan, index = df.index)
    return df[(ts == variant) & (θ - θ0).abs().le(1e-12)]


def countryθ0(country, variant, path = None):
    """ θ0 of a workbook country ('US', 'FR', 'UK'; ISO3) under `variant`: thetaFromFile of its ISO3 row. """
    return thetaFromFile(ISO3[country], variant, path)


def paperConfig():
    """ python/paper/config.py as a module, loaded from the file (python/US imports nothing from python/paper). """
    import importlib.util
    sp = importlib.util.spec_from_file_location('_paperConfigθ', PAPERCONFIG)
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


def paperTheta():
    """ python/paper/config.py's US['theta'], the pipeline's main θ variant, read from the file. """
    return paperConfig().US.get('theta')


def thetaRoot(variant, main = False):
    """ The directory a θ variant's run reads and writes (calibration/, esc/, shocks/ below it): results/ under
    `main`, which must be the pipeline's main variant (paperTheta; ValueError otherwise), thetaDir(variant)
    without it. """
    if main:
        if variant is None or variant == INDEX or variant != paperTheta():
            raise ValueError("--main writes results/ for the pipeline's main θ variant, config.US['theta'] = {!r}; "
                             'got --theta {!r}'.format(paperTheta(), variant))
        return RESULTS
    return thetaDir(variant)


def underResults(path):
    """ Is `path` results/ or below it? """
    q, r = os.path.normcase(os.path.abspath(path)), os.path.normcase(RESULTS)
    return q == r or q.startswith(r + os.sep)


def indexTheta(theta, paths, driver):
    """ A U.S. driver's --theta where it may mean the index θ (the workbooks' RR0, figure 1.1's index): returns None
    for None and INDEX (the driver's no-variant path), any other value unchanged. Without --theta (None) a run that
    would write any of `paths` under results/ exits: the standard paths hold the pipeline's main variant
    (--theta config.US['theta'] --main), and --theta index writes the index θ there on purpose. """
    if theta == INDEX:
        return None
    if theta is None:
        bad = [p for p in paths if p is not None and underResults(p)]
        if bad:
            raise SystemExit("{}: without --theta the rows are at the index θ, and {} is under results/, which holds "
                             "the main variant's (--theta config.US['theta'] --main): pass --theta {} to write the "
                             'index θ there, or an output outside results/'.format(
                                 driver, os.path.relpath(bad[0], REPO), INDEX))
    return theta


def thetaSample(variant, path = None):
    """ The countries with inSample True under `variant`, in data/oecd29.csv's order, the U.S. always included. """
    d = readTheta(path)
    d = d[d['variant'] == variant]
    if d.empty:
        raise ValueError('{}: no variant {!r}'.format(shown(path or FILES['theta']), variant))
    keep = set(d.loc[d['inSample'], 'iso3']) | {'USA'}
    return [c for c in countries() if c in keep]


def thetaDir(variant):
    """ results/oecd26_theta_<variant>/, a θ variant's own directory: the layers' csvs, calibration/ (its U.S.
    sweep) and esc/ (its cost calibration). """
    return os.path.join(REPO, 'results', 'oecd26_theta_' + variant)


def sweepPath(root):
    """ The U.S. sweep of a variant's directory `root`: calibration/US_rhoGridCommonX.csv. """
    return os.path.join(root, 'calibration', 'US_rhoGridCommonX.csv')


def thetaFromIndex(I):
    """ θ = 2 - 1/I, clipped to [0, 1] (I = 1/(2 - θ) inverted): Australia (I < 1/2) and Ireland sit at 0, the
    countries with I = 1 at 1. """
    return float(np.clip(2. - 1./float(I), 0., 1.))


def indexFromTheta(θ):
    """ I = 1/(2 - θ): a design on the scale of figure 1.1's index. """
    return 1./(2. - np.asarray(θ, dtype = float))


def nuPath(nu4):
    """ ν over DATES from the four WPP values (1990, 2020, 2050, 2080). """
    n = [float(nu4[k]) for k in range(4)]
    return np.array([n[0], n[0], n[1], n[2], n[3], n[3]] + [n[3]]*T_SS)


def readExitAge(path = None):
    """ data/oecd29_exitAge.csv's exitAge_men (men's average effective age of labour-market exit, 2020), on ISO3. """
    d = pd.read_csv(FILES['exit'] if path is None else path, index_col = 'iso3', float_precision = 'round_trip')
    return d['exitAge_men'].astype(float)


def retAgeGaps(iso3, retAge, main = None):
    """ g = A_US - A_c in years at NUYEARS, (4,): 'exit' readExitAge's at every date; 'nra' nra_current_men at 1990
    and 2020, nra_future_men at 2050 and 2080 (`main` readData()[0]); 'none' zeros. A missing age raises. """
    if retAge not in RETAGES:
        raise ValueError('retAge must be one of {}, got {!r}'.format(RETAGES, retAge))
    if retAge == 'none':
        return np.zeros(4)
    if retAge == 'exit':
        A = readExitAge()
        g = np.full(4, float(A['USA']) - float(A[iso3]))
    else:
        main = readData()[0] if main is None else main
        cols = ['nra_current_men']*2 + ['nra_future_men']*2
        g = main.loc['USA', cols].to_numpy(float) - main.loc[iso3, cols].to_numpy(float)
    if not np.isfinite(g).all():
        raise ValueError('{}: no retirement age under retAge {!r}'.format(iso3, retAge))
    return g


def retAgeNu(nu4, g):
    """ Appendix F's ν(30 - g)/(30 + ν g) at each date of nu4 (PERIOD = 30) for a country retiring g > 0 years before
    the U.S.; for one retiring |g| years later its mirror image (30ν + |g|)/(30 - |g|), the extra years worked by the
    old cohort (python/paper/oecd29Data.nraRule, so 'nra' is data/oecd29_population.csv's wpp2024_nraAdj). A zero
    gap leaves ν bitwise. """
    ν, g = np.asarray(nu4, dtype = float), np.asarray(g, dtype = float)
    return np.where(g == 0., ν, np.where(g > 0., ν*(PERIOD - g)/(PERIOD + ν*g), (PERIOD*ν - g)/(PERIOD + g)))


def readData():
    """ The three data files the loader reads, and the sample of the data step (python/paper/oecd26Sample.py). """
    main = pd.read_csv(FILES['main'], index_col = 'iso3')
    dist = pd.read_csv(FILES['dist'])
    dist = dist[(dist['role'] == 'data') & (dist['source'] == DISTSOURCE) & (dist['convention'] == DISTCONV)]
    sample = pd.read_csv(FILES['sample'], index_col = 'iso3')
    return main, dist, sample


def countries():
    """ The 29 ISO3 codes in data/oecd29.csv's order. """
    return list(pd.read_csv(FILES['main'], usecols = ['iso3'])['iso3'])


def turnoutSample(thetaSource = INDEX, data = None):
    """ The countries `pooledProfile` averages over: the θ source's sample, thetaSample under a variant and
    data/oecd26_sample.csv's sample26 under 'index', in data/oecd29.csv's order. """
    if thetaSource != INDEX:
        return thetaSample(thetaSource)
    sample = (readData() if data is None else data)[2]
    return [c for c in countries() if bool(sample.at[c, 'sample26'])]


def pooledProfile(thetaSource = INDEX, data = None):
    """ (μ (3,), countries): the mean of μ_c/μ_c,1 (mu1-3_usPct) over the countries of turnoutSample other than the
    U.S. and COMPULSORY with CSES turnout; μ_1 = 1. """
    main = (readData() if data is None else data)[0]
    keep = [c for c in turnoutSample(thetaSource, data)
            if c != 'USA' and c not in COMPULSORY and main.loc[c, MU].notna().all()]
    μ = main.loc[keep, MU].to_numpy(float)
    return (μ/μ[:, :1]).mean(axis = 0), keep


def countryInputs(iso3, turnout, data = None, thetaSource = INDEX, retAge = 'none'):
    """ One country's inputs (module docstring) as a dict; `data` a readData() triple to reuse. thetaSource
    'index' (figure 1.1's index) or a variant of data/oecd26_theta.csv (NoTheta where it has no θ). turnout one of
    TURNOUTS; under 'pooled' μ depends on thetaSource's sample (pooledProfile). retAge one of RETAGES (ν's
    correction, retAgeNu). """
    if turnout not in TURNOUTS:
        raise ValueError('turnout must be one of {}, got {!r}'.format(TURNOUTS, turnout))
    main, dist, sample = readData() if data is None else data
    d = main.loc[iso3]
    g = dist[dist['iso3'] == iso3].sort_values('group')
    if len(g) != 3:
        raise ValueError('{}: {} rows of {}/{} in {}'.format(iso3, len(g), DISTSOURCE, DISTCONV, FILES['dist']))
    γ, y = g['gamma'].to_numpy(float), g['y'].to_numpy(float)
    if turnout == 'cps':
        μ = np.asarray(testmod.pars['μj'], dtype = float)[1:]
    elif turnout == 'pooled' and iso3 in COMPULSORY:
        μ = np.ones(3)
    elif turnout == 'pooled' and iso3 != 'USA':
        μ = pooledProfile(thetaSource, (main, dist, sample))[0]
    elif turnout == 'csesFlat' and iso3 in COMPULSORY and d[MU].isna().any():
        μ = np.ones(3)
    else:                                 # 'cses', 'csesFlat', and the U.S. under 'pooled'
        μ = d[MU].to_numpy(float)
        if np.isnan(μ).any():
            raise NoTurnout('{}: no CSES turnout in {}'.format(iso3, shown(FILES['main'])))
    nu4 = retAgeNu(d[['nu_{}'.format(t) for t in NUYEARS]].to_numpy(float), retAgeGaps(iso3, retAge, main))
    spending = float(d['pensionSpending'])
    if thetaSource == INDEX:
        index, θRaw = float(d['grr100_men'])/float(d['grr050_men']), None
    else:
        θRaw = thetaFromFile(iso3, thetaSource)
        index = float(indexFromTheta(θRaw))
        sFile = spendingFromFile(iso3, thetaSource)
        if thetaSource in SPENDINGASFIGURE and not abs(sFile - spending) <= SPENDINGTOL:
            raise ValueError('{}: spending {!r} under variant {!r} is not figure 1.1\'s pensionSpending {!r}'.format(
                iso3, sFile, thetaSource, spending))
        spending = sFile
    inp = inputsFrom(iso3, d['country'], nuPath(nu4), γ, y, μ, spending, index,
                     turnout = turnout, sample26 = bool(sample.at[iso3, 'sample26']), thetaSource = thetaSource)
    if θRaw is not None:
        inp['θRaw'] = θRaw
    return inp


def inputsFrom(iso3, country, ν, γ, y, μ, spending, index, turnout, sample26 = None, thetaSource = INDEX):
    """ The inputs dict of countryInputs from its parts: ν over DATES (T,), γ_i and y_i (3,), μ_i (3,), spending
    in % of GDP and the replacement-rate ratio I; α, ξ and the workweek from the U.S. workbook. thetaSource is
    recorded only (countryInputs reads θ). """
    α, ξ = float(testmod.pars['α']), float(testmod.pars['ξ'])
    γ, y, ν = np.array(γ, dtype = float), np.array(y, dtype = float), np.array(ν, dtype = float)   # copies
    if ν.shape != (len(DATES),):
        raise ValueError('{}: ν must have one value per date of DATES ({}), got {}'.format(iso3, len(DATES), ν.shape))
    return {'iso3': iso3, 'country': country, 'turnout': turnout, 'thetaSource': thetaSource, 'sample26': sample26,
            'ν': ν, 'γ': γ, 'y': y,
            'μ': np.array(μ, dtype = float), 'Vtilde': float(np.sum(γ*(y - 1)**2/y)), 'spending': float(spending),
            'τ0': float(spending)/100/(1 - α), 'index': float(index), 'θ': thetaFromIndex(index),
            'θRaw': 2. - 1./float(index), 'α': α, 'ξ': ξ, 'workweek': float(testmod.workweek)}


def usBeta(path = USSWEEP, ρ = 1.0, thetaSource = INDEX):
    """ (β, ω, source) of the U.S. calibration without the cost at `ρ` (the sweep's row, finite residual): β is
    imposed on every economy here, ω only seeds a root. `source` names the file and the row's commit and
    timestamp, which the drivers record so a row computed before a recalibration is recognisable (#13).
    Under a θ variant `path` is the variant's own sweep (sweepPath; FileNotFoundError while it is absent), whose
    row must carry that thetaSource and θ0 = the variant's U.S. θ (ValueError otherwise); `source` adds both. """
    if thetaSource != INDEX and not os.path.exists(path):
        raise FileNotFoundError('{}: no U.S. sweep for θ variant {!r}: run calibrateRhoGrid.py --commonX --theta {} '
                                '--lo 1 --hi 1 first'.format(path, thetaSource, thetaSource))
    df = pd.read_csv(path)
    df = df[np.isclose(df['ρ'].astype(float), ρ) & np.isfinite(pd.to_numeric(df['residual'], errors = 'coerce'))]
    if df.empty:
        raise KeyError('{}: no converged row at rho = {}'.format(path, ρ))
    r = df.iloc[-1]
    src = '{} rho={:g} commit={} {}'.format(os.path.relpath(path, REPO).replace(os.sep, '/'), ρ,
                                            r.get('commit', ''), r.get('timestamp', ''))
    if thetaSource != INDEX:
        θ = thetaFromFile('USA', thetaSource)
        if r.get('thetaSource') != thetaSource or not abs(float(r.get('θ0', np.nan)) - θ) <= 1e-12:
            raise ValueError('{}: the row at rho = {} was calibrated at thetaSource {!r}, θ0 {!r}; variant {!r} has '
                             'U.S. θ {!r}: rerun calibrateRhoGrid.py --theta {} --force'.format(
                                 path, ρ, r.get('thetaSource'), r.get('θ0'), thetaSource, θ, thetaSource))
        src += ' theta={} θ0={!r}'.format(thetaSource, θ)
    return float(r['β']), float(r['ω']), src


class ModelOECD(ModelESC):
    """ ModelESC on one country's OECD inputs: θ from the datum RR0 = I (figure 1.1's index, or 1/(2 - θ) of a θ
    variant), never from db['θ0'] (ModelUS.θOverride), which getθ refuses. The cost must be proportional ('size',
    'sizePiR' or 'scale', or none), under which it cancels from the replacement-rate ratio and I = 1/(2 - θ)
    holds; under 'flat' it does not.

    betaRule 'common': β imposed by `build`, _calPars ('ω',) against τ_{t0} (R is a prediction). 'commonR':
    _calPars ('β', 'ω') against (R_{t0}/R0 - 1, τ_{t0} - τ0), ModelUS's eq:calibration, or ('β',) against
    R_{t0}/R0 - 1 alone when calibrateω is False (ω held at the value `build` installed). Under 'sizePiR'
    ModelESC.calibrate wraps either root in the π iteration. """

    _calBounds = {k: (0., np.inf) for k in ('β', 'ω')}

    def __init__(self, *args, betaRule = 'common', calibrateω = True, **kwargs):
        if betaRule not in BETARULES:
            raise ValueError('ModelOECD: betaRule must be one of {}, got {!r}'.format(BETARULES, betaRule))
        self.betaRule = betaRule
        self._calPars = (('β',) if betaRule == 'commonR' else ()) + (('ω',) if calibrateω else ())
        super().__init__(*args, **kwargs)

    def getθ(self):
        if self.db.get('wedgeSpec') == 'flat':
            raise NotImplementedError("ModelOECD: θ = 2 - 1/I needs a proportional cost; 'flat' is not one.")
        if self.θOverride() is not None:
            raise ValueError("ModelOECD: θ enters through RR0 = 1/(2 - θ) (countryInputs' thetaSource); "
                             "db['θ0'] is not read here.")
        return thetaFromIndex(self.db['RR0'])

    def _calSetPars(self, pars):
        """ The parameters of _calPars, each through adjPar: under 'common' β is installed once by `build` and
        never written from `pars`. """
        for k in self._calPars:
            self.db.update(self.adjPar(k, pars[k]))
        self.updateAuxPars()

    def _calResidual(self, report):
        """ One entry per element of _calPars, same order: R_{t0}/R0 - 1 for β, τ_{t0} - τ0 for ω. """
        return np.array(([report['R']/self.db['R0'] - 1] if 'β' in self._calPars else [])
                        + ([report['τ'] - self.db['τ0']] if 'ω' in self._calPars else []))


def build(inp, β, ω, λ = None, spec = 'size', phi = 0.5, gs = None, nθ = 41, nθCand = 121, betaRule = 'common',
          R0 = None, calibrateω = True):
    """ A ModelOECD on `inp` (countryInputs) at β (imposed under betaRule 'common', the root's start under
    'commonR'), with ω and the cost parameter λ (the 'p' slot; None = no cost), the LOG grid `gs` (GS by default)
    installed. R0: the 2020 return 'commonR' calibrates β to, required there; NaN otherwise, so that a stray R
    target fails. calibrateω: ModelOECD's. Calendar DATES, t0 = 2020. """
    if betaRule == 'commonR' and R0 is None:
        raise ValueError("build: betaRule 'commonR' needs the target R0")
    T = len(DATES)
    pars = {'α': inp['α'], 'ξ': inp['ξ'], 'ν': np.array(inp['ν'], dtype = float), 'τ0': inp['τ0'],
            'RR0': inp['index'], 'R0': np.nan if R0 is None else float(R0), 'Ushare0': 0.,
            't0': DATES.get_loc(T0YEAR),
            'γj': np.r_[0., inp['γ']], 'μj': np.r_[1., inp['μ']], 'Xj': 1, 'zxj': np.ones(4),
            'zηj': np.r_[inp['y'][0], inp['y']],    # the zero-mass slot cloned off type 1: finite, unused
            'h0': inp['workweek']/(7*12), 'ρ': 1., 'β': float(β), 'ω': float(ω)}
    wedge = None if λ is None else {'spec': spec, 'phi': phi, 'p': float(λ)}
    m = ModelOECD(pars = pars, wedge = wedge, nθ = nθ, nθCand = nθCand, commonX = True, nj = 4, T = T,
                  betaRule = betaRule, calibrateω = calibrateω)
    m.db['dates'], m.db['workweek'], m.db['country'] = DATES, inp['workweek'], inp['iso3']
    m.LOG.initGS(GS if gs is None else gs)
    return m


def closedFormTax(α, β, ξ, ω, ν, terminal = False):
    """ The log tax at θ_t = 1, archive/notes/reviewResponse_2026-10/theory_closedForm.md eq (1.4); β -> 0 in the terminal period. """
    Ψ = (1 + ξ)*(1 + (0. if terminal else α*β*(1 + ξ)/(1 + α*ξ)))
    return max(0., (ω*(1 - α) - α*ν*Ψ)/(ω*(1 - α)*(1 + ξ) + (1 - α)*ν*Ψ))


def closedFormCheck(m):
    """ The model's tax path at θ ≡ 1 against closedFormTax at each date's ν: (max over the dates of |gap|, the
    number of dates at the corner τ = 0, which the solver returns exactly). """
    out = m.solvePEE_LOG(θ = np.ones(m.T))
    τ = out['τ'].values.astype(float)
    t0 = m.t0Year
    α, ξ, ω = float(m.db['α'].xs(t0)), float(m.db['ξ'].xs(t0)), float(m.db['ω'].xs(t0))
    cf = np.array([closedFormTax(α, m.simpleβinv(), ξ, ω, ν, terminal = (pos == m.T - 1))
                   for pos, ν in enumerate(m.db['ν'].values.astype(float))])
    return float(np.max(np.abs(τ - cf))), int(np.sum(cf == 0.))


def inputsKey(*parts):
    """ 'h' and 12 hex digits hashing everything a row depends on (inputs, imposed and anchor parameters, grid
    settings), numbers at 12 significant digits; the letter keeps a csv reader from taking it for a number. A
    resumable row is reused only under the same key (#13). """
    def canon(v):
        if isinstance(v, dict):
            return {str(k): canon(v[k]) for k in sorted(v, key = str)}
        if isinstance(v, (list, tuple, np.ndarray, pd.Series)):
            return [canon(x) for x in list(np.asarray(v, dtype = object).ravel())]
        if isinstance(v, (bool, np.bool_)) or v is None:
            return v if v is None else bool(v)
        if isinstance(v, (int, float, np.integer, np.floating)):
            return '{:.12g}'.format(float(v))
        return str(v)
    return 'h' + hashlib.sha1(json.dumps(canon(list(parts)), sort_keys = True).encode()).hexdigest()[:12]


def modelInputs(inp):
    """ The part of `inp` a solve depends on, for inputsKey. """
    return {k: inp[k] for k in ('iso3', 'turnout', 'ν', 'γ', 'y', 'μ', 'τ0', 'index', 'α', 'ξ')}

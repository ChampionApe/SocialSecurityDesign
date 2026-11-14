r""" Paths, the Argentina specification, and the unit conversions shared by all three pipeline stages.

This module imports nothing from the model. Stage (iii) must stay able to rebuild every table and figure
from `results/` alone, in seconds and without unpickling a model instance -- see README.md.
"""
import os, sys, functools
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))

DATA      = os.path.join(REPO, 'data')
RESULTS   = os.path.join(REPO, 'results')
CALIBDIR  = os.path.join(RESULTS, 'calibration')
SHOCKDIR  = os.path.join(RESULTS, 'shocks')
SWEEPDIR  = os.path.join(RESULTS, 'sweeps')
ESCDIR    = os.path.join(RESULTS, 'esc')
NUMDIR    = os.path.join(RESULTS, 'numerical')      # pre-publication numerical checks (stationary vs
                                                    # date-specific policy functions), no paper output
PAPERDIR  = os.path.join(RESULTS, 'paper')          # stage (iii) writes here first
PAPERTEX  = os.path.join(REPO, 'writing', 'Paper')  # ... then copies here unless --no-copy

MODELDIR  = os.path.join(REPO, 'python', 'InformalSavings')
USDIR     = os.path.join(REPO, 'python', 'US')
PYTHON    = sys.executable                           # the stage's own interpreter runs its children, on any OS
# The power k on the tax in the cost f (python/US/modelESC.TAU_POWER reads the same variable; the
# branch's default is 2). TAUS_*: tau0^(1-k), the factor the calibrated lambda moves by, for the scan brackets.
TAU_POWER = float(os.environ.get('SSD_TAUPOWER', 2.))
TAUS_US, TAUS_UK = 0.105229**(1. - TAU_POWER), 0.111943**(1. - TAU_POWER)
# The carrier of the Pi-inclusive cost 'sizePiR' (python/US/modelESC.WEDGE_CARRIER reads the same
# variable; the branch's default is 'convex'): 'pi' lambda_{n,t} = epsilon pi_{n,t}, 'convex' lambda_{n,t} =
# 2 kappa pi_{n,t}(1 - pi_{n,t}), pi = Pi/R. PIRSYMBOL: the anchor's `p` as the builders print it (math mode, no
# dollars); PIRCOST: the cost a country faces as the notes print it, dollars included, subscript n,t (PIRCOSTN: n
# alone, for the captions). Under 'pi' each is the text the builders printed before the carrier existed.
CARRIER = os.environ.get('SSD_CARRIER', 'convex')
CARRIERS = ('pi', 'convex')
if CARRIER not in CARRIERS:
    raise ValueError('SSD_CARRIER={!r}: one of {}'.format(CARRIER, CARRIERS))
PIRSYMBOL = {'pi': r'\varepsilon', 'convex': r'\kappa'}[CARRIER]
PIRCOST = {'pi': r'$\lambda_{n,t} = \varepsilon\,\Pi_{n,t}/R_{n,t}$',
           'convex': r'$\lambda_{n,t} = 2\kappa\,\pi_{n,t}(1-\pi_{n,t})$, $\pi_{n,t} = \Pi_{n,t}/R_{n,t}$'}[CARRIER]
PIRCOSTN = {'pi': r'$\lambda_n = \varepsilon\,\Pi_n/R_n$',
            'convex': r'$\lambda_n = 2\kappa\,\pi_n(1-\pi_n)$, $\pi_n = \Pi_n/R_n$'}[CARRIER]


def carrierFactor(π):
    """ g(π) of CARRIER, as python/US/base.carrierFactor (stage (iii) imports no model code): π under 'pi',
    2π(1-π) under 'convex'. """
    return π if CARRIER == 'pi' else 2*π*(1 - π)

# ---------------------------------------------------------------------------------------------------
# The Argentina specification: what the paper's numbers are, as a declaration rather than as CLI history.
# Stages (i) and (ii) pass these to the experiment scripts; changing a paper number starts here.
# ---------------------------------------------------------------------------------------------------
ARG = {
    'workbook':   'ArgentinaTest.xlsx',
    'ρGrid':      [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0],
    'ρAnchor':    1.0,      # the sole LOG point; calibrateGrid marches outward from it
    'ρBaseline':  1.0,      # the paper's headline specification
    'ρTable':     [0.5, 1.0, 2.0],   # rows the funcOfRho table prints uncommented; the rest are
                                     # emitted commented. Matches US['ρTable'], so both arms of the
                                     # paper show the same three points of their rho grids
    'rule':       'match',  # b^0 = b^refType
    'refType':    1,
    # The anchor's starting (beta, omega, eta0, X0). The march seeds every other point from its history,
    # but the anchor starts from test.py's defaults (beta = 0.6, omega = 2), and at alpha = 0.35 those put
    # the informal steady state outside the net-saver region, so the iota state grid is degenerate before
    # the root takes a step. This is the vector-X solution at rho = 1 for alpha = 0.35, K/Y = 3.23 and the
    # GDP-based tax target (model.taxTarget, tau0 = 0.1181), z_0 relative to the formal average
    # (test_calibrationGrid.py, 2026-10-23); retune it if the capital share, the spending share, the K/Y
    # target or the household-survey targets move.
    'anchorGuess': {'β': 0.668641, 'ω': 1.496254, 'η0': 0.285587, 'X0': 0.373999},
    # Grid settings. calibrateRhoGrid.py gives BOTH solvers interpKind/smoothKnots and only the grid
    # SIZES to CRRA; LOG keeps its own documented nι=50. Anything re-solving a calibrated instance must
    # mirror that split or it solves under a different interpolant than it was fitted under
    # (archive/notes/informalSavings_resolvedIssues.md). loadCalibrated() enforces it; do not bypass.
    'gridSettings': {'interpKind': 'cubic', 'smoothKnots': 4, 'nι': 45, 'ns': 45},
    # WHICH CALIBRATION VARIANT THE ARGENTINA OUTPUTS LEAD WITH (as US['commonX']). False: vector X_i
    # from relative income and relative hours; True: one scalar X across the formal types pinned by the
    # formal workweek, relative formal hours a prediction. X enters no aggregate, so tau, K/Y, the
    # savings rate and every counterfactual coincide across the two to solver precision (measured at
    # <= 4e-12 over the whole rho grid, archive/notes/argentina_commonX_vs_vectorX.md); what differs is the
    # calibration table (eta_i, X_i, eta_0, X_0) and the relative-hours diagnostic. The headline outputs
    # read this; their twins pass commonX = not this.
    # Common X, matching the OECD arm: it identifies one leisure parameter
    # rather than four, puts the workweek level on a target, and turns relative formal hours into a
    # prediction to be checked against the survey.
    'commonX': True,
}


def argVariantTag(commonX = False):
    """ The file-name suffix of the Argentina variant's results: '' (vector X) or '_commonX'. """
    return '_commonX' if commonX else ''


def argSweepCsv(commonX = False):
    return os.path.join(CALIBDIR, 'informalSavings_rhoGrid' + ('CommonX' if commonX else '') + '.csv')


def argInstanceDir(commonX = False):
    return os.path.join(CALIBDIR, 'instances' + ('CommonX' if commonX else ''))


def argShockTemplate(scenario, rule = None, commonX = False):
    """ The shock csv name template of one scenario ('reform' | 'ee' | 'flat'), with {ρ} left open. """
    rule = rule or ARG['rule']
    stem = {'reform': 'universal_' + rule, 'ee': 'eeOnly_' + rule, 'flat': 'universal_flat'}[scenario]
    return stem + '_rho{ρ:.4f}' + argVariantTag(commonX) + '.csv'


def argEpsThetaCsv(ρ, commonX = False):
    return os.path.join(SWEEPDIR, 'epsThetaGrid_rho{:.4f}{}.csv'.format(ρ, argVariantTag(commonX)))

# ---------------------------------------------------------------------------------------------------
# Calendar and units
# ---------------------------------------------------------------------------------------------------
# The pickled instances do NOT carry db['dates'] or db['workweek'] -- test.py sets them on its own
# mLOG, and calibratePoint pickles an instance that never saw them. Both are re-read from the workbook.


@functools.lru_cache(maxsize = None)
def calendar():
    """ {model t index: calendar year} and the calibration year's index, from the workbook.
    Mirrors test.py's construction: the population sheet's dates, extended by t_ss=3 further 30-year
    steps. Only the dated part is returned -- the steady-state tail has no calendar meaning. """
    wb = pd.read_excel(os.path.join(DATA, ARG['workbook']), sheet_name = None, header = None)
    dft = pd.DataFrame(wb['population'].values)
    dates = pd.DataFrame(dft.iloc[1:, ].values, columns = dft.iloc[0, :]).set_index('t').index
    dfc = pd.DataFrame(wb['calibration'].values)
    dfc = pd.Series(dfc.iloc[1, :].values, index = dfc.iloc[0, :].values)
    τ0, y = argTaxTarget(wb)
    return {'dates': {i: int(d) for i, d in enumerate(dates)},
            't0': int(list(dates).index(dfc['Calibration year'])),
            'year0': int(dfc['Calibration year']),
            'workweek': float(dfc['Average workweek']),
            # Derived as in test.py (model.taxTarget), never a literal tax rate.
            'τ0': τ0, 'yInformal': y,
            'α': float(dfc['Capital income share']),
            'spending': float(dfc['Pension spending']),
            's0': float(dfc['Savings rate'])}


def argTaxTarget(wb):
    """ (tau0, y) from the workbook `wb` (pd.read_excel's dict of sheets): eq (calibration:tau)'s target, pension
    spending over GDP with GDP formal plus informal output and y the informal share of GDP. The closed form of
    python/InformalSavings/model.py's taxTarget, restated because this module imports nothing from the model,
    on test.py's inputs: a = gamma_0 z_0^eta (1-alpha), z_0^eta the informal income over the population-weighted
    formal mean, tau0 = spending (1+a)/((1-alpha) + spending a), y = a(1-tau0)/(1+a(1-tau0)).
    datasets.calibrationSummary checks tau0 against the calibrated tax. """
    dfc = pd.DataFrame(wb['calibration'].values)
    dfc = pd.Series(dfc.iloc[1, :].values, index = dfc.iloc[0, :].values)
    dfj = pd.DataFrame(wb['heterogeneity'].values)
    dfj = pd.DataFrame(dfj.iloc[1:, ].values, columns = dfj.iloc[0, :])
    γ, inc = dfj['Population shares'].values.astype(float), dfj['Income'].values.astype(float)
    z0η = inc[0]/np.average(inc[1:], weights = γ[1:])
    S, α = float(dfc['Pension spending']), float(dfc['Capital income share'])
    a = γ[0]*z0η*(1-α)
    τ0 = S*(1+a)/((1-α) + S*a)
    return τ0, a*(1-τ0)/(1 + a*(1-τ0))


# ---------------------------------------------------------------------------------------------------
# The US/France/UK specification. Same role as ARG above: a declaration of what the paper's rich-OECD
# numbers are, not a second implementation. See python/US/README.md for the model and the protocol.
# ---------------------------------------------------------------------------------------------------
US = {
    'countries':  ('US', 'FR', 'UK'),        # the columns of OECD_Calibration (the supplemental appendix's)
    'workbooks':  {'US': 'USMain_test.xlsx', 'FR': 'FRMain.xlsx', 'UK': 'UKMain.xlsx'},
    # THE PAPER'S MEASURE OF THE DESIGN: a variant of data/oecd26_theta.csv (notes/thetaTau_alternatives.md),
    # imposed as theta0 on the U.S., France and the UK (their USA, FRA, GBR rows) and taken by every OECD country.
    # Every U.S./France/UK/ESC/OECD command of stages (i) and (ii) carries --theta <this> --main, which writes the
    # standard results/ paths (python/US/oecd26.thetaRoot); every reader of those paths keeps the rows at this
    # thetaSource only (datasets.atTheta), a csv written without the column reading as 'index'.
    'theta':      'C',
    # FRANCE'S nu IN ITS OWN CALIBRATION (python/US/testEU.load; every France run of stages (i) and (ii)):
    # 'nra' the population sheet's 'Worker-to-retiree' (appendix F: the WPP ratio corrected for the gap to the
    # U.S. men's normal retirement age) | 'exit' the WPP ratio ('Relative cohort sizes') corrected for the gap to the
    # U.S. men's effective exit age (oecd26.retAgeNu), as the cross-section's FRA row under US['oecd26']['retAge'].
    # 'retAgeUK' the same for the UK: 'nra' its sheet's column, the WPP ratio (its normal retirement ages are the
    # U.S. ones, a gap of 0) | 'exit' that ratio at its gap in effective exit ages.
    'retAgeFR':   'nra',
    'retAgeUK':   'nra',
    'ρGrid':      [0.5, 1.0, 2.0],            # the sweeps' points (calibrateRhoGrid*.py --rhos)
    'ρAnchor':    1.0,                        # the sole LOG point, and the sweeps' march anchor
    'ρBaseline':  1.0,                        # the headline specification: the LOG tables
    'ρTable':     [0.5, 1.0, 2.0],            # the rows the CRRA tables print
    # WHICH CALIBRATION VARIANT THE PAPER LEADS WITH. True = the common scalar X of the docs' common-X calibration,
    # where the hours unit is a calibration target and relative hours become a prediction; False = the
    # vector X_i calibration, where relative hours are data and the LEVEL of hbar is not identified.
    # beta, omega, tau, R, the savings rate and aggregate h are the SAME under both (block recursivity,
    # measured to ~1e-13 down the sweeps) -- what the choice changes is eta, X, and therefore every
    # counterfactual defined THROUGH eta or X: the French income-distribution and leisure rows.
    # The headline tables and figures read this; their vector-X twins pass commonX = not this.
    'commonX': True,
    # gradKnotSpacing is recorded, not forwarded (python/US/policy.CRRA's default): the distance in tau between
    # the knots of the least-squares cubic in ln(1 - tau) the CRRA tau-derivatives are taken from. Appendix
    # D's residual table (tablesHorizon.numResiduals) prints it.
    'gridSettings': {'interpKind': 'linear', 'smoothKnots': 4, 'n': 101, 'ns': 150,
                     'verify': 225, 'verifyN': 151, 'gradKnotSpacing': 0.04},
    # --- Endogenous system characteristics (sec:esc): the leaded choice of theta under a deadweight
    # cost on the pension system (python/US/base.py fWedge, writing/US/model_esc.tex). Every spec is
    # PROPORTIONAL -- the cost scales the whole benefit and cancels from the replacement-rate ratio, so
    # theta* stays the data's own -- and its one parameter is calibrated per rho so the design IN FORCE in
    # 2020 on a freely simulated path is the observed one (ModelESC.leadedDesignAtT0,
    # results/esc/escCalibration{,CRRA}.csv). The counterfactual tables report at t0: every scenario is a
    # new equilibrium path whose political choice binds from the first period, so 2020's design is an
    # outcome and already carries the response. The ESC leg runs under the headline variant, US['commonX'].
    #   'size'  (the paper's): f(theta, tau) = exp(-1/2 lambda tau^2 Vtilde (1-theta)^2), a deadweight cost
    #           quadratic in the implicit tax the flat component levies on each type and in the size of the
    #           system, lambda = 2 kappa pi (1 - pi) held constant (modelESC.TAU_POWER = 2; tau^1 is paper/tau1); Vtilde = sum_i gamma_i (y_i-1)^2/y_i. Nothing is lost when nothing is
    #           redistributed. lambda lives in the csvs' `p` column (the paper prints it as lambda).
    #   'scale' (the previous wedge, the supplemental appendix's second cost specification): f(theta) =
    #           phi + (1-phi) theta^p, the cost attached to the design label rather than to the transfer.
    #   'flat'  (only the flat component carries the cost): no run since the permanent stage was retired.
    # The rho points of the pre-publication stationary-vs-date-specific check (runShocksUS.py 'stationary';
    # rho = 1 is exact by the LOG decoupling and is not run): the non-log points of ρGrid. sec:numerical's
    # footnote quotes the maximal gap over these points.
    'ρStationary': [0.5, 2.0],
    'esc': {
        'spec':           'size',      # the paper's cost specification
        'comparisonSpec': 'scale',     # the previous wedge, the supplemental appendix's second cost (page C2)
        # phi is a DUMMY KEY under 'size': f does not use it, but every ESC csv is merge-keyed on
        # (spec, phi, commonX) and the readers (runESC.pickCalib, datasets.escRow/escCalibration) filter on
        # it, so the 'size' rows are written and read at 0.5. Under 'scale' it is the imposed floor f(0).
        'phi':            0.5,
        'ρTable':         [0.5, 1.0, 2.0],
        # The scan bracket of the EXACT CRRA calibration (runESCcrra.py --bracket), per spec and rho: the
        # required parameter falls steeply in rho, so one bracket cannot serve both. A rho absent from a
        # spec's dict leaves the bracket to runESCcrra.py's own default for that spec; a rho present with
        # None is a PLACEHOLDER and stage (i) refuses to run it (config.escBracket). 'size' (the paper's cost, the
        # one exact calibration run): [lambda/20, 3 lambda] at rho = 2 and [lambda/3, 5 lambda] at rho = 0.5 of
        # lambda = 11.563954 (the LOG lambda at the pension-only workbook). At US['theta'] the LOG lambda is 6.247
        # (results/esc/escCalibration.csv), and the exact roots expected there, ~13.8 and ~1.46, lie inside both.
        # Lambda scales by about tau0^(1-k) (U.S. tau0 0.105229, the UK's 0.111943), so do the brackets.
        'bracket':        {'size':    {0.5: (3.855*TAUS_US, 57.82*TAUS_US), 2.0: (0.5782*TAUS_US, 34.69*TAUS_US)}},
        # brentq's xtol on the cost parameter where runESCcrra.py's default (2e-5) is too coarse, per spec and rho.
        'xtol':           {},
        # THE PUBLISHED CRRA METHOD. True: every CRRA ESC output is built from rows with method = 'exact'
        # (LeadedCRRA2D, python/US/runESCcrra.py --exact, the pre-publication part of stages (i)/(ii))
        # and a missing exact row is MissingInput -- never a fallback to the path iteration's rows, which
        # stay on disk under method = 'path' as the cross-check they are. False builds from the path rows.
        'exact':     True,
        'ns2D':      150,      # savings-state grid of the exact recursion (its choice is insensitive to it)
        'nsScan':    50,       # the coarse grid the wedge calibration SCANS on before refining at ns2D
        # The candidate grid for θ_{t+1}. The objective is flat near its maximum (1e-5 in W over ±0.01 in
        # θ at the frVoting choice, rho = 2), so this grid sets the resolution of the printed design: 13
        # nodes gave 0.285, 21 gave 0.273.
        'nCand2D':   41,
        # The design layer of the exact recursion (runESCcrra.py --designRule/--Ma): algorithm esc:crra2D of
        # writing/US/num_esc.tex, the root in a, adopted on the pilot (python/US/RESEARCH_LOG.md, 2026-10-02).
        'designRule': 'root',
        'Ma':         5,
        # The UK's OWN cost parameter per rho (the deep test, UK_ESC_Calibration), calibrated as the US one is:
        # LOG from escCountry.csv ('own'), CRRA by the exact recursion into escCalibrationCRRAUK.csv (stage (i)
        # --prepub, runESCcrra.py --host UK --stage calib). The brackets are scaled to the UK's LOG lambda at its
        # θ0 of US['theta'] (0.181; 3.628 on the coarse grid), the exact lambda being 2.08 and 0.38 times the LOG
        # one at rho = 0.5 and 2: ~7.8 and ~1.45 expected, inside the brackets, scanned at nScan nodes. nsScan 0 at
        # rho = 0.5 scans at ns2D: the UK's (beta-imposed) calibration does not converge on the ns = 50 grid there.
        'uk': {'bracket':   {0.5: (3.5*TAUS_UK, 16.0*TAUS_UK), 2.0: (0.55*TAUS_UK, 3.2*TAUS_UK)},
               'nsScan':    {0.5: 0, 2.0: 50},
               'nScan':     5},
        # THE COST VARIANT OF SECTION 7 (A8 of archive/notes/reviewResponse_2026-10/paper_reviewResponse.md): the
        # Pi-inclusive cost 'sizePiR' (python/US/modelESC.py: lambda*Pi_t/R_t, the factor the paper's cost absorbs
        # into lambda, solved to the path's own Pi_t/R_t), calibrated like the paper's cost under LOG into the
        # `tag`ged csv escCalibrationCosts.csv (its epsilon is the cross-section's), and acute ageing against its
        # baseline under it at rho = 1 (escShocksCosts.csv); collectESCexperiments.py --costs merges it with the
        # paper's csvs into escAgeingCosts.csv.
        'costVariants': {'spec': 'sizePiR', 'tag': 'Costs', 'scenarios': ['baseline', 'acute']},
        # escExperiments.csv scenario keys -> the labels the appendix tables print.
        # 'frLeisure' is run and merged but no longer printed (dropped from the paper 2026-09-11).
        'scenarios': {'acute': 'Acute ageing', 'frIncome': 'Income distribution',
                      'frVoting': 'Voting', 'frAll': 'All French characteristics'},
    },
    # THE CROSS-COUNTRY TEST (python/US/runOECD26.py, results/oecd26/), the one place its paper variant lives:
    # the builders (tablesOECD26.PAPER, figuresOECD26) and the readers (datasets.oecd26*) default to it, and
    # runOECD26.py's --spec/--betaRule/--retAge defaults are read from here.
    #   spec      the cost transported from the U.S. anchor: 'sizePiR' (lambda_n = epsilon Pi_n/R_n, or under the
    #             carrier 'convex' 2 kappa pi_n (1 - pi_n) with pi_n = Pi_n/R_n (CARRIER above); epsilon or kappa
    #             the U.S. one), 'size' (lambda itself), 'none' (no cost)
    #   betaRule  'common' (the U.S. beta in every country) | 'commonR' (each country's beta set so that its
    #             2020 return is the U.S. anchor's)
    #   turnout   'csesFlat' (each country's own CSES propensities to vote; Luxembourg, which has no CSES survey and
    #             enforces compulsory voting, flat: 21 countries under C; it stands on the 'cses' anchor) | 'cses'
    #             (Luxembourg left out, 20 countries)
    #   omegaRule layer 1's omega: 'mean' (the mean of layer 2's; beta and the cost parameter the anchor's) |
    #             'anchor' (the U.S. anchor's)
    #   reading   the model's design: 'inForce' (in force in 2020, chosen in 1990) | 'chosen' (chosen in 2020)
    #   spendingReading  the model's spending in section 7's test (OECD26_Layer1, OECD26_Slopes, OECD26_Variants):
    #             'chosen' (100 tauInForce (1 - alpha), the tax on the path whose design in force is the chosen one,
    #             in both layers) | 'observed' (100 tau2020 (1 - alpha), the tax at the observed design, which layer
    #             2's omega matches). Section 6's spending test (OECD26_Spending) is always 'observed'.
    #   spendingSpec  the cost of section 6's spending test at the observed designs
    #   retAge    nu's retirement-age correction in every country but the U.S. (python/US/oecd26.retAgeNu, appendix F's
    #             rule): 'nra' (the gap to the U.S. men's normal retirement age, current at 1990 and 2020, future at
    #             2050 and 2080; France's nu is then appendix F's) | 'none' (the WPP ratio) | 'exit' (the effective exit
    #             age). A key column of the layers' rows (runOECD26.LEGACY 'none'); the anchor does not depend on it
    # `variants`: the runs stage (ii)'s entry 'oecd26' makes into results/oecd26/, one runOECD26.py command per retAge
    # of each (--theta US['theta'] --main --inSample), in this order (an anchor and layer 2 before the mean omega that
    # reads them); rows already current are reused, so a rerun costs seconds. They are the paper's specification with
    # and without the cost, appendix G's population variants (`populationVariants`) and the supplemental appendix's
    # readings (the last two entries, each the parts that differ from the paper's: Luxembourg left out, layer 1 at
    # the U.S. omega). `thetaVariants`: appendix G's data variants, measures of data/oecd26_theta.csv (THETANAMES),
    # each run by stage (ii)'s entry 'oecd26Variants' at the paper's specification (a U.S. sweep at rho = 1, then the
    # cross-section) into its own results/oecd26_theta_<V>/.
    'oecd26': {
        'spec': 'sizePiR', 'betaRule': 'common', 'turnout': 'csesFlat', 'omegaRule': 'mean', 'reading': 'inForce',
        'spendingSpec': 'none', 'spendingReading': 'chosen', 'retAge': 'nra',
        'variants': [
            {'spec': ['sizePiR', 'none'], 'betaRule': ['common'], 'turnout': ['csesFlat'],
             'stage': ['anchor', 'layer2', 'layer1'], 'omegaRule': 'mean', 'retAge': ['nra']},
            {'spec': ['sizePiR'], 'betaRule': ['common'], 'turnout': ['csesFlat'],
             'stage': ['layer2', 'layer1'], 'omegaRule': 'mean', 'retAge': ['none', 'exit']},
            {'spec': ['sizePiR'], 'betaRule': ['common'], 'turnout': ['cses'],
             'stage': ['anchor', 'layer2', 'layer1'], 'omegaRule': 'mean', 'retAge': ['nra']},
            {'spec': ['sizePiR'], 'betaRule': ['common'], 'turnout': ['cses'],
             'stage': ['layer2', 'layer1'], 'omegaRule': 'mean', 'retAge': ['none', 'exit']},
            {'spec': ['sizePiR'], 'betaRule': ['common'], 'turnout': ['csesFlat'],
             'stage': ['layer1'], 'omegaRule': 'anchor', 'retAge': ['nra']},
        ],
        'populationVariants': ['none', 'exit'],
        'thetaVariants': ['A', 'B'],
    },
    # THE TAX AGAINST A PINNED DESIGN (python/US/runTauTheta.py, results/esc/escTauTheta.csv; US_ESC_TauTheta): the
    # U.S. at rho = 1 without the cost and under US['esc']['spec'] at its calibrated lambda, theta pinned on the whole
    # path at each point of `grid` and at theta0, and at `refine` steps within one grid cell of the costly economy's
    # grid maximum.
    'tauTheta': {'grid': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 'refine': 0.01},
}

THETAFILE = os.path.join(DATA, 'oecd26_theta.csv')
ISO3 = {'US': 'USA', 'FR': 'FRA', 'UK': 'GBR'}
# US['theta']'s measure as every generated note states it, after "$\theta$": one wording for all of them
# (python/paper/oecd26Theta.py, notes/thetaTau_alternatives.md, variant C).
THETADEF = (r'the earnings-linked share of the least-squares line through the public pension at 0.5, 1 and 2 times '
            r'average earnings ($c/(a+c)$ of $p(y) = a + c\,y$), times one minus the share of public spending on '
            r'minimum top-ups and tested benefits')
# The names the paper gives the measures of the design (US['theta'], the paper's, and appendix G's data variants
# US['oecd26']['thetaVariants']) and the readings of nu (US['oecd26']['retAge'] and appendix G's population variants).
THETANAMES = {'C': 'public tier', 'A': 'public tier, 23 countries', 'B': 'all mandatory tiers'}
RETAGENAMES = {'nra': 'statutory age', 'none': 'uncorrected', 'exit': 'exit age'}


def thetaFlags():
    """ The flags every U.S./France/UK/ESC/OECD command of stages (i) and (ii) carries: --theta US['theta'] --main. """
    return ['--theta', US['theta'], '--main']


@functools.lru_cache(maxsize = None)
def _thetaFile():
    return pd.read_csv(THETAFILE, dtype = {'iso3': str, 'variant': str}, float_precision = 'round_trip')


def theta0(country = 'US', variant = None):
    """ θ0 of a workbook country ('US', 'FR', 'UK') or an ISO3 code under `variant` (US['theta'] by default): its row of
    data/oecd26_theta.csv, as python/US/oecd26.thetaFromFile reads it. KeyError where there is none. """
    variant = US['theta'] if variant is None else variant
    d = _thetaFile()
    r = d[(d['iso3'] == ISO3.get(country, country)) & (d['variant'] == variant)]
    if r.empty or not np.isfinite(float(r['theta'].iloc[0])):
        raise KeyError('{}: no θ for {} under variant {!r}'.format(THETAFILE, country, variant))
    return float(r['theta'].iloc[0])


def atTheta(df, variant = None, country = None, tol = 1e-12):
    """ The rows of `df` at thetaSource `variant` (US['theta'] by default; a csv or cell without the column reads as
    'index'), and with `country` also at its θ0 (the θ0 column, theta0). Every reader of a U.S.-arm csv under
    results/{calibration,esc,shocks,oecd26} goes through this, so a row of another θ is never read (#13). """
    variant = US['theta'] if variant is None else variant
    ts = df['thetaSource'].fillna('index') if 'thetaSource' in df.columns else pd.Series('index', index = df.index)
    keep = ts == variant
    if country is not None and variant != 'index':
        θ = pd.to_numeric(df['θ0'], errors = 'coerce') if 'θ0' in df.columns else pd.Series(np.nan, index = df.index)
        keep &= (θ - theta0(country, variant)).abs().le(tol)
    return df[keep]


def escBracket(ρ, spec = None):
    """ The exact CRRA calibration's scan bracket at `rho` under `spec` (US['esc']['spec'] by default):
    (lo, hi) to pass as runESCcrra.py --bracket, or None to leave the child its own default. Raises on a
    placeholder (rho listed with None) so an unset bracket stops stage (i) before a 45-minute scan. """
    spec = US['esc']['spec'] if spec is None else spec
    table = US['esc']['bracket'].get(spec, {})
    hit = [table[k] for k in table if np.isclose(float(k), ρ)]
    if not hit:
        return None
    if hit[0] is None:
        raise ValueError("config.US['esc']['bracket'][{!r}][{}] is a placeholder: set it from the LOG "
                         "calibration before running the exact CRRA calibration".format(spec, ρ))
    lo, hi = hit[0]
    return float(lo), float(hi)


# Which sweep csv belongs to which country and calibration variant. The US is calibrated on its own
# (calibrateRhoGrid.py); FR/UK impose the US beta and are swept by calibrateRhoGridEU.py.
def usSweepCsv(country, commonX = False):
    tag = ('US_rhoGrid' if country == 'US' else country + '_rhoGrid') + ('CommonX' if commonX else '')
    return os.path.join(CALIBDIR, tag + '.csv')


def usInstanceDir(country, commonX = False):
    return os.path.join(CALIBDIR, 'instances' + ('US' if country == 'US' else country)
                        + ('CommonX' if commonX else ''))


def usShockCsv(host = 'US', commonX = False):
    """ python/US/runShocksUS.py's long csv for `host` ('US') and variant. """
    return os.path.join(SHOCKDIR, host + '_shocks' + ('CommonX' if commonX else '') + '.csv')


@functools.lru_cache(maxsize = None)
def usCalendar(country = 'US'):
    """ {model t index: calendar year}, the calibration year, and the observed workweek, from `country`'s
    workbook. Mirrors python/US/test.py and testEU.py: the population sheet's dates extended by t_ss = 3
    further 30-year steps, of which only the dated part is returned.

    Read from the workbook rather than from a pickled instance for the same reason ARG's calendar is:
    calibratePoint pickles an instance that never saw db['dates']. On a shock COPY db['dates'] is worse
    than absent -- it is present and stale (python/US/test_createCopyFromt0.py), which is the other
    reason nothing in this pipeline reads it.

    """
    wb = pd.read_excel(os.path.join(DATA, US['workbooks'][country]), sheet_name = None, header = None)
    dft = pd.DataFrame(wb['population'].values)
    dates = pd.DataFrame(dft.iloc[1:, ].values, columns = dft.iloc[0, :]).set_index('t').index
    dfc = pd.DataFrame(wb['calibration'].values)
    dfc = pd.Series(dfc.iloc[1, :].values, index = dfc.iloc[0, :].values)
    return {'dates': {i: int(d) for i, d in enumerate(dates)},
            't0': int(list(dates).index(dfc['Calibration year'])),
            'year0': int(dfc['Calibration year']),
            'workweek': float(dfc['Average workweek']),
            'τ0': float(dfc['Pension tax']),
            'RR0': float(dfc['Replacement rate']),
            'α': float(dfc['Capital income share'])}


def clearCaches():
    """ Forget the cached workbook and θ-file reads (calendar, usCalendar, _thetaFile). build.py --map traces each
    output's file reads, which a read cached by an earlier output would hide. """
    for f in (calendar, usCalendar, _thetaFile):
        f.cache_clear()


# ---------------------------------------------------------------------------------------------------
# The two US calibration variants, as they appear in an output's identity and in its table note.
#
# The variant US['commonX'] names is the HEADLINE: its outputs keep the plain names and tex labels the
# paper already \ref{}s, so switching which variant leads never renames the outputs the draft cites --
# it changes what they contain. The other variant is the robustness twin, and its name and label carry
# ITS OWN variant rather than the word "alternative", so a file on disk says what is in it.
#
# The printed note is keyed on commonX, not on which variant leads: the paper's tables and figures are
# the common-X calibration and say nothing about it, and only a vector-X table names itself.
# ---------------------------------------------------------------------------------------------------
def isLead(commonX, arm = 'US'):
    """ Is `commonX` the headline variant of `arm` ('US' | 'ARG')? """
    return bool(commonX) == bool((US if arm == 'US' else ARG)['commonX'])


def variantSuffix(commonX, arm = 'US'):
    r""" '' for the headline variant, '_commonX'/'_vectorX' for the twin. Appended to an output's name,
    its tex filename and its \label. """
    return '' if isLead(commonX, arm) else ('_commonX' if commonX else '_vectorX')


def variantCaption(commonX, arm = 'US'):
    """ The caption tail that marks a twin. Empty for the headline -- its caption is the paper's own. """
    return '' if isLead(commonX, arm) else (r' (common $X$ calibration)' if commonX
                                            else r' (vector $X_i$ calibration)')


def variantNote(commonX, full = False, arm = 'US'):
    r""" One sentence naming the calibration variant, appended to a VECTOR-X table's note.

    Empty under common X, which is what the paper is: every table and figure in the main text is the
    common-X calibration, so none of them says so, and only the vector-X twins in the appendices carry a
    label. The vector-X calibration table (`full = True`) spells the variant out; every other vector-X
    table names it and points at that table's note. """
    if commonX:
        return ''
    calib = r'table:US:Calib' if arm == 'US' else r'table:Arg:Calib'
    if not full:
        return (r' Vector-$X_i$ calibration: see the note to '
                r'table~\ref{' + calib + variantSuffix(commonX, arm) + '}.')
    return (r' Vector-$X_i$ calibration: $X_i$ is identified from relative hours, which are data here, '
            r'and the level of $\bar h$ is then not identified --- only its ratio to the baseline is. '
            r'$\beta$, $\omega$, the tax rate and the savings rate are common to the two variants; what '
            r'differs is anything defined through $\eta$ or $X$.')


def workweekHours(h, hRef):
    """ Aggregate hours h_t -> an average workweek in hours, by NORMALISATION against `hRef`.

    `h` has no well-defined scale in the model, so no expression converts it to hours on its own. The
    observed average workweek is a REFERENCE POINT instead: the calibrated model's hours at the
    calibration year are defined to be `workweek` hours, and any other h is reported as
    `workweek * h/hRef`. Pass the calibrated baseline's own h at t0 as `hRef` -- per rho, since each rho
    is separately calibrated. The baseline then reads 42.54 at every rho by construction, and a shocked
    h is read as the change in hours it implies.

    NOT `h * 7 * 12`. That inverts test.py's `h0 = workweek/(7*12)`, which is how the PRE-DETERMINED
    period's hours are fed in as a model input; it is not a scale the solved h_t inherits, and using it
    to report makes the baseline miss 42.54 and manufactures a spread across rho out of a free
    normalisation. """
    return calendar()['workweek'] * np.asarray(h, dtype = float) / hRef


def pct(x, digits = 2):
    r""" 0.1256 -> '12.56\%'. The escaped percent is what goes into a tex cell. """
    return r'{:.{d}f}\%'.format(100*x, d = digits)


def pp(x, digits = 2):
    r""" A CHANGE in a rate, 0.0014 -> '$+0.14$ p.p.', -0.0158 -> '$-1.58$ p.p.'. Math mode so the
    minus is a minus rather than a hyphen. A change that rounds to zero prints unsigned ('0.00 p.p.')
    rather than as '-0.00': the sign of a solver residual is not a result. """
    v = round(100*x, digits)
    if v == 0:
        return '{:.{d}f} p.p.'.format(0., d = digits)
    return '${:+.{d}f}$ p.p.'.format(v, d = digits)


def num(x, digits = 2):
    return '{:.{d}f}'.format(x, d = digits)


def vec(x, digits = 1):
    """ [1.03, 1.45] -> '[1.0, 1.5]', the form the calibration table already uses. """
    return '[' + ', '.join('{:.{d}f}'.format(v, d = digits) for v in x) + ']'

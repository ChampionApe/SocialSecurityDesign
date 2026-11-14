r""" Endogenous system characteristics, A+B: run the leaded choice of theta under the deadweight wedge and
write everything the appendix needs.

Run:  .venv\Scripts\python.exe python\US\runESC.py                     # every stage
      ... --stage calib                        the wedge calibration table alone
      ... --stage path                         the endogenous design path
      ... --stage shocks                       the counterfactuals
      ... --stage country                      France and the UK
      ... --spec scale flat  --phi 0.25 0.5 0.75
      ... --spec size --phi 0.5 --commonX --stage calib path shocks country     # the paper's wedge

Four stages, each writing its own csv under results/esc/:

  calib    for every (spec, phi): p such that the equilibrium design IN FORCE at 2020 is the observed one
           (ModelESC.leadedDesignAtT0), with (beta, omega) recalibrated inside. Plus the no-wedge row,
           which is the corner the whole exercise is trying to escape. Every row also carries Vtilde,
           fStar = f(theta*, tau0), f0 = f(0, tau0) and tau0 (wedgeReadout), which the paper's calibration
           table prints; under 'size' the 'p' column IS lambda and 'phi' is the dummy 0.5.
  path     the equilibrium path of theta_t under the calibrated wedge, with the choice binding from the
           first period -- the ageing prediction, and the closest thing the model has to figure 1.1's
           cross-section.
  shocks   the counterfactuals of shocks.py (ageing, French income distribution / leisure / voting, and
           all three at once), plus a France reference row, each run twice: theta PINNED at the shocked
           model's exogenous design, and theta chosen. Every one is a NEW EQUILIBRIUM PATH -- the shocked
           parameters hold over the whole horizon, the economy starts at its own steady state, the
           political choice binds from the first period -- and is REPORTED AT t0 (2020). t0+1 is written
           too, but the tables read t0: because theta_{t0} is now an equilibrium outcome rather than a
           pinned datum, 2020 already carries the design response.
  country  France and the UK under (i) the US-calibrated wedge and (ii) their own. (i) asks whether one
           common technology of redistribution puts the three countries in the observed order; (ii) is the
           appendix's own convention (a separately calibrated p per country).
  pathHost the path stage on France and the UK, under the US-calibrated p and under their own
           (escCountry.csv) -> escPathHost.csv.
  corner   'size' only: the slope of the leaded design objective at the Bismarckian corner, theta_{t+1} = 1,
           and the tax there, for France, the U.S. and the UK at lambda = 0, the U.S. lambda and 100, the
           1990 and the 2020 electorate (app:esc:cost: under log both are free of lambda, so a negative slope
           means no lambda re-elects the corner) -> escCornerSlope.csv. Each row is recalibrated at its
           lambda, so the slope is lambda-free only where the observed design is the corner (France at its
           workbook θ = 1): an economy observed below it (the U.S. and the UK; France too under a θ variant, at
           countryθ0) has omega, and with it the slope, move with the cost.

--theta VARIANT [--root DIR | --main] [--coarse] (stages calib, path, pathHost, shocks, country, corner; `setTheta`):
every U.S. model imposes theta as db['θ0'] = the U.S. row of data/oecd26_theta.csv under VARIANT (ModelUS.getθ,
oecd26.thetaFromFile), so the calib stage puts that design in force in 2020, (beta, omega) recalibrated inside as
always; France and the UK impose theirs, the FRA and GBR rows (oecd26.countryθ0, `buildEU`, `frenchData`), with omega
recalibrated at it (France at the U.S. income cuts). The variant's directory DIR (default
results/oecd26_theta_<variant>/, oecd26.thetaDir; results/ under --main, the pipeline's main variant,
oecd26.thetaRoot) holds its U.S. sweep, DIR/calibration/ (calibrateRhoGrid.py --theta: required, at that θ0), which
seeds the U.S. and is France's and the UK's U.S. reference (beta imposed, hours). The csvs go to DIR/esc/ (results/esc/
only under --main), their rows adding thetaSource, θ0 (the U.S.'s) and grid, and a stage that reads
escCalibration.csv or escCountry.csv refuses the rows it reads at another thetaSource, θ0 or grid (checkCalib).
--coarse (variant only): GRIDS['coarse']. --theta index: the workbook's θ into results/esc/; without --theta main
refuses to run (results/esc/ is its only output, oecd26.indexTheta).

Spec 'sizePiR' (ModelESC: the 'size' cost with lambda*Pi_t/R_t, pi solved to the path's own) runs the calib,
path, shocks and country stages; its rows carry 'ΠR', pi at the row's date (t0 for calib, shocks and country).
Run it with a --tag: it is an experiment beside the paper's cost, not a row of the paper's csvs.

Everything is LOG (rho = 1). The CRRA case needs a two-dimensional (s, theta) state and is not this file.

A failure a stage survives (a scenario, France's row, a country's cost) is printed, the run goes on, and main exits 1
with the list at the end (`noteFailure`); a calibration that does not converge is a row with converged False, not a
failure. A relative --root is the caller's working directory's (the script itself runs in python/US).
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
import testEU
import shocks as sh
from modelESC import ModelESC
from modelFR import ModelFR
from policy import multiplicityColumns

OUTDIR = os.path.join(REPO, 'results', 'esc')
GS = {'n': 101, 'smoothKnots': 4, 'interpKind': 'linear'}
# The tax grid's nodes and LeadedLOG's design state grid and candidates of buildUS/buildEU/frenchData.
GRIDS = {'fine': {'n': 101, 'nθ': 41, 'nθCand': 121}, 'coarse': {'n': 51, 'nθ': 21, 'nθCand': 41}}
# The run's θ variant, set by setTheta (main's --theta) and read by every builder here, so no stage can build a
# U.S. at another θ: θ0 imposed on the U.S. (db['θ0']), calibDir the variant's U.S. sweep (buildUS's seed, France's
# and the UK's U.S. reference), grid a GRIDS key. The defaults are the workbook's θ, results/calibration/, 'fine'.
THETA = {'variant': None, 'θ0': None, 'calibDir': None, 'grid': 'fine', 'main': False}


def countryθ0(country):
    """ θ0 of France or the UK under the run's variant (oecd26.countryθ0), None without one. """
    if THETA['variant'] is None:
        return None
    import oecd26
    return oecd26.countryθ0(country, THETA['variant'])


FAILURES = []                      # the failures this run survived (noteFailure); main exits 1 when there are any


def noteFailure(where, e):
    """ Print and record a failure a stage survives: the run goes on, main exits 1 with the list. """
    FAILURES.append('{}: {}: {}'.format(where, type(e).__name__, e))
    print('  FAILED ' + FAILURES[-1])


def writeCsv(df, path):
    """ `df` to csv `path` through a temporary file in its directory and os.replace, so that an interrupted write
    leaves the file as it was. """
    tmp = path + '.tmp'
    try:
        df.to_csv(tmp, index = False)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def mergeWrite(path, rows, key):
    """ Write `rows` to csv `path`, PRESERVING rows on disk whose `key` combination this run did not
    produce. The stage loops write after every point, and the paper pipeline re-runs exactly the missing
    (spec, phi) or (rho, spec, phi) combinations -- a plain to_csv would clobber everything else in the
    file on such a partial run. NaN keys (the no-wedge row's phi) are compared as equal.

    `commonX` IS PART OF EVERY KEY HERE, and it has to be: the two calibration variants answer the same
    (rho, spec, phi) question about two different economies, so without it a common-X run would silently
    overwrite the vector-X rows -- or, worse, leave them in place to be read as if they were current. A
    resumable producer must be keyed on the question AND on what is answering it
    (notes/crossCuttingFindings.md #13).

    Keys are compared as CANONICAL strings (_keyOf): a key column that mixes bools with NaN -- legacy rows
    kept beside new ones -- comes back from the csv as 1.0/NaN, not True/NaN, and a naive str() would then
    never match the run's own rows again and duplicate them on every write. """
    df = pd.DataFrame(rows)
    if os.path.exists(path) and len(df):
        # round_trip: the default parser keeps 17 digits counting the zeros after the decimal point, so a rewrite
        # would truncate the rows it merely carries (a 1e-4 diagnostic lost its last three digits).
        old = pd.read_csv(path, float_precision = 'round_trip')
        # A csv written before one of the key columns existed has no way to match a new row, and must
        # not be dropped on the strength of the columns it does share -- so the missing column is added
        # as NaN, which the '__na__' fill turns into a key no current row can equal. The old rows are
        # then preserved verbatim and stay distinguishable from the new ones.
        for k in key:
            if k not in old.columns:
                old[k] = np.nan
        df = pd.concat([old[~_keyOf(old, key).isin(_keyOf(df, key))], df], ignore_index = True)
    writeCsv(df, path)


def _keyCell(v):
    """ One key value as a string that survives the csv round trip: NaN/NA/'' -> '__na__'; bools, ints
    and floats through repr(float(v)) (True == 1.0 == 1, 3 == 3.0); everything else str(v). """
    if isinstance(v, (bool, np.bool_, int, np.integer, float, np.floating)):
        f = float(v)
        return '__na__' if np.isnan(f) else repr(f)
    if v is None or v is pd.NA or (isinstance(v, str) and v == ''):
        return '__na__'
    return str(v)


def _keyOf(d, key):
    """ The row keys of frame d over the columns `key`, canonicalised cell by cell (_keyCell). """
    return pd.MultiIndex.from_arrays([[_keyCell(v) for v in d[k].tolist()] for k in key])


# Every row also carries the tax rule's counts of the solve behind it, nEqMax, nCandMax, nFallback, as
# non-key columns (policy.multiplicityColumns; the permanent timing's over its candidate designs and the
# equilibrium at the choice, modelESC.solvePermanent; -1 where there is no solve: a calibration that did not
# converge).

# The row-identity keys of the merged csvs, one per stage. spec and commonX are in every one -- see
# mergeWrite -- so the 'size' rows coexist with the 'scale' rows (the comparison arm) in one file and a
# rerun replaces only its own rows. No csv here carries rho: this driver is LOG at rho = 1; runESCcrra's
# keys add it.
KEYCAL = ['spec', 'phi', 'commonX']
KEYSHK = ['spec', 'phi', 'commonX', 'scenario', 'θpinned']
KEYPATH = ['spec', 'phi', 'commonX', 'pos']
KEYPATHHOST = ['spec', 'phi', 'commonX', 'country', 'wedgeFrom', 'pos']
KEYCTY = ['spec', 'phi', 'commonX', 'country', 'wedgeFrom']
KEYCORNER = ['spec', 'phi', 'commonX', 'country', 'λ', 'vote']


def pickCalib(calib, spec, phi, commonX):
    """ The converged wedge calibration for one (spec, phi, variant), or an empty frame.

    Filtering on commonX is not optional: escCalibration.csv carries both variants, and reading a
    vector-X p into a common-X run would report a wedge that was never calibrated for that economy. """
    hit = calib[(calib['spec'] == spec) & (calib['phi'] == phi) & calib['converged']]
    if 'commonX' in calib.columns:
        hit = hit[hit['commonX'].astype(bool) == bool(commonX)]
    return hit


class ModelESCFR(ModelESC, ModelFR):
    """ France/UK under the wedge. ModelESC first in the MRO so its __init__ strips `wedge` before
    ModelFR's strips `usRef`; calibrate() then resolves to ModelFR's (beta imposed from the US, omega the
    only searched parameter) and getθ to ModelESC's wedge-aware one. """


# ---------------------------------------------------------------- builders

def usSweep(commonX = False):
    """ The U.S. sweep of this run: results/calibration/US_rhoGrid{,CommonX}.csv, or THETA['calibDir']'s. """
    return os.path.join(THETA['calibDir'] or os.path.join(REPO, 'results', 'calibration'),
                        'US_rhoGridCommonX.csv' if commonX else 'US_rhoGrid.csv')


def _usRef(commonX):
    """ France's and the UK's U.S. reference (testEU.usReference) off this run's U.S. sweep. """
    return testEU.usReference(commonX = commonX, csv = None if THETA['calibDir'] is None else usSweep(commonX))


def thetaColumns():
    """ The columns a θ variant's rows add: thetaSource, θ0 and grid ({} without a variant). """
    return {} if THETA['variant'] is None else {'thetaSource': THETA['variant'], 'θ0': THETA['θ0'],
                                                 'grid': THETA['grid']}


def buildUS(wedge = None, ρ = 1.0, nθ = None, nθCand = None, commonX = False, θ0 = None):
    """ The US model under the wedge, with (beta, omega) SEEDED from the no-wedge sweep at this rho (usSweep).
    θ0: the design imposed as db['θ0'] (ModelUS.getθ), THETA's by default; None derives it from the workbook's
    RR0. nθ, nθCand and the tax grid: GRIDS[THETA['grid']] unless given.

    Only seeded: every stage that uses this calls calibrate() (or calibrateWedge, which calibrates
    inside), so the numbers on the csv are a starting point, not the calibration. What the variant
    selects is therefore the sweep file to seed from and -- through commonX -- the meaning of eta and X.
    Under commonX the hours unit is a target rather than free, so the calibration adds h0 as a third
    condition and returns a scalar X; nothing here has to install it, calibrate() does.
    """
    g = GRIDS[THETA['grid']]
    θ0 = THETA['θ0'] if θ0 is None else θ0
    row = pd.read_csv(usSweep(commonX))
    row = row.loc[(row['ρ'] - ρ).abs() < 1e-9].iloc[-1]
    m = ModelESC(pars = testmod.pars | {'ρ': float(ρ), 'β': float(row['β']), 'ω': float(row['ω'])}
                 | ({} if θ0 is None else {'θ0': float(θ0)}),
                 wedge = wedge, nθ = g['nθ'] if nθ is None else nθ, nθCand = g['nθCand'] if nθCand is None else nθCand,
                 commonX = commonX, **testmod.kwargs)
    m.db['dates'] = testmod.dates
    m.db['workweek'] = testmod.workweek
    m.LOG.initGS(GS | {'n': g['n']})
    return m


def buildEU(country, wedge = None, nθ = None, nθCand = None, ρ = 1.0, commonX = False, ω = 2., **escKw):
    """ France or the UK under the wedge. commonX must MATCH the US model it is compared against: under
    commonX France's X_i is a single scalar and its eta_i inverts z^eta directly, so both characteristics
    differ from the vector-X calibration and a mismatched pair would put a France that does not exist in
    this run's units next to the US rows (runShocksUS.frenchData says the same about the shock data).
    ω only seeds calibrate(); escKw goes to ModelESC (nθCandCRRA, nθCand2D). θ is the country's own: its RR0, or
    under a θ variant db['θ0'] = countryθ0(country); the U.S. reference is this run's sweep (_usRef), the grids
    GRIDS[THETA['grid']] unless given. """
    g = GRIDS[THETA['grid']]
    pars, kwargs, dates, workweek = testEU.load(country)
    θ0 = countryθ0(country)
    pars.update({'ρ': float(ρ), 'ω': float(ω)} | ({} if θ0 is None else {'θ0': θ0}))
    # usReference must be read off the SAME variant's sweep: it carries hbar_US, and hbar differs
    # between vector X and common X by construction (python/US/README.md).
    m = ModelESCFR(pars = pars, wedge = wedge, usRef = _usRef(commonX), nθ = g['nθ'] if nθ is None else nθ,
                   nθCand = g['nθCand'] if nθCand is None else nθCand, commonX = commonX, **escKw, **kwargs)
    m.db['dates'], m.db['workweek'], m.db['country'] = dates, workweek, country
    m.LOG.initGS(GS | {'n': g['n']})
    return m


# ---------------------------------------------------------------- reporting

def readout(m, τ, report, hbarRef, pos, workweekData = None):
    """ shocks.readout at one position of the path.

    workweekData defaults to the model's own datum, which is right for every US-based counterfactual. It
    must be passed explicitly for a NON-US model (franceRow): the workweek is reported as
    workweekData * hbar/hbarRef with hbarRef the US baseline's, so pairing France's hbar with France's own
    datum would scale the ratio twice and put its row at ~31.9 instead of 35.4. """
    return sh.readout(m, τ, report,
                      float(m.db['workweek'] if workweekData is None else workweekData),
                      hbarRef, pos = pos)


def wedgeReadout(m, τ0):
    """ The wedge read at t0 on a calibrated model, as four csv columns: Vtilde (the young's vintage at t0),
    fStar = f(theta*, tau0), f0 = f(0, tau0) and tau0 -- tau_{t0} on the calibrated path, which the
    calibration's own target db['tau0'] reproduces to tauDrift. Under 'scale'/'flat' f ignores tau (and
    under no wedge f = 1); tau0 is written regardless so the column means one thing in every row. f is
    that of the t0 retirees (lag '[t-1]': their Vtilde, equal to the young's in every calibration here, and
    under 'sizePiR' pi_t0, written as a fifth column 'ΠR', and carrierColumn's). """
    t0 = m.t0Year
    θStar = float(m.db['θ'].xs(t0))
    with m.B.cacheParams():
        out = {'Vtilde': float(m.B.Vtilde(t0)), 'fStar': float(m.B.fWedge(θStar, τ0, t0, '[t-1]')),
               'f0': float(m.B.fWedge(0., τ0, t0, '[t-1]')), 'τ0': float(τ0)}
        if m.db.get('wedgeSpec') == 'sizePiR':
            out['ΠR'] = float(m.B.wedgeΠR(t0, '[t-1]'))
    return out | carrierColumn(m)


def carrierColumn(m):
    """ {'carrier': g's name} under 'sizePiR' at a carrier other than 'pi' (base.carrierFactor), {}
    otherwise, so the rows of the 'pi' carrier keep their columns. """
    c = m.db.get('wedgeCarrier', 'pi')
    return {'carrier': c} if m.db.get('wedgeSpec') == 'sizePiR' and c != 'pi' else {}


def piColumn(m, t):
    """ {'ΠR': pi_t} under 'sizePiR' (and carrierColumn's), {} otherwise: the return-factor ratio of a path row. """
    return ({'ΠR': float(m.db['wedgeΠR'].xs(t))} | carrierColumn(m)) if m.db.get('wedgeSpec') == 'sizePiR' else {}


def baselineRefs(m):
    """ The exogenous-theta baseline on the SAME model (same wedge, same calibration): its hbar at t0 is
    the workweek's reference point, and its (tau, sr, workweek) the row every endogenous number is read
    against. """
    base = m.solvePEE_LOG()
    t0 = m.t0Year
    hbarRef = float(m.B.avgHours(base['report']['h'].xs(t0), t0))
    return base, hbarRef


# ---------------------------------------------------------------- stage: calibration

def stageCalib(specs, phis, out, ρ = 1.0, commonX = False, bracket = None):
    """ bracket: the scan bracket for the wedge parameter, None = ModelESC.WEDGE_BRACKET[spec]. """
    rows = []
    tagθ = thetaColumns()
    m = buildUS(None, ρ = ρ, commonX = commonX)
    m.calibrate()
    t0 = m.t0Year
    sols = m.ESC.solveBackward()
    rows.append({'spec': 'none', 'phi': np.nan, 'commonX': commonX, 'p': np.nan, 'converged': True,
                 'θStar': float(m.db['θ'].xs(t0)), 'choice': m.leadedDesignAtT0(sols),
                 'β': m.simpleβinv(), 'ω': float(m.db['ω'].xs(t0)), 'τDrift': np.nan, 'RDrift': np.nan}
                | wedgeReadout(m, float(m.db['τ0'])) | multiplicityColumns(m.ESC.lastMultiplicity) | tagθ)
    print('no wedge: θ*={θStar:.4f} -> choice {choice:.4f}  (Vtilde={Vtilde:.4f})'.format(**rows[-1]))
    mergeWrite(out, rows, KEYCAL)

    for spec in specs:
        for phi in phis:
            tic = time.time()
            print(f'\n[{spec}, φ={phi}] calibrating p ...')
            m = buildUS({'spec': spec, 'phi': phi, 'p': 1.0}, ρ = ρ, commonX = commonX)
            try:
                rec = m.calibrateWedge(spec = spec, phi = phi, bracket = bracket)
            except Exception as e:
                print(f'  FAILED: {type(e).__name__}: {e}')
                rows.append({'spec': spec, 'phi': phi, 'commonX': commonX, 'p': np.nan,
                             'converged': False} | multiplicityColumns(None) | tagθ)
                mergeWrite(out, rows, KEYCAL); continue
            r = {'spec': spec, 'phi': phi, 'commonX': commonX, 'p': rec['p'],
                 'converged': rec['converged'],
                 'θStar': rec['θ'], 'residual': rec['residual'],
                 'β': m.simpleβinv(), 'ω': float(m.db['ω'].xs(m.t0Year)),
                 'Vtilde': float(m.B.Vtilde(m.t0Year))}
            if rec['converged']:
                led = m.solveLeaded(pinAtT0 = False)
                # 'choice' is the calibration's own target read back off the solved path: the design in
                # FORCE at t0, not the one chosen there. At a converged p it equals θ* to `residual`.
                r |= {'τDrift': led['targetDrift']['τ'], 'RDrift': led['targetDrift']['R'],
                      'choice': float(led['θ'].iloc[m.db['t0']]),
                      'choiceAtT0': float(led['θ'].iloc[m.db['t0']+1])}
                r |= wedgeReadout(m, float(led['τ'].xs(m.t0Year)))
            r |= multiplicityColumns(led['multiplicity'] if rec['converged'] else None) | tagθ
            rows.append(r)
            print('  -> p={p}  θ*={θStar:.4f}  β={β:.4f} ω={ω:.4f}  Vtilde={Vtilde:.4f} f(θ*)={fStar} f(0)={f0}  ({:.0f}s)'
                  .format(time.time()-tic, **({'fStar': np.nan, 'f0': np.nan} | r)))
            mergeWrite(out, rows, KEYCAL)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- stage: the design path

def stagePath(specs, phis, calib, out, ρ = 1.0, commonX = False):
    rows = []
    for spec in specs:
        for phi in phis:
            hit = pickCalib(calib, spec, phi, commonX)
            if hit.empty:
                continue
            p = float(hit.iloc[0]['p'])
            m = buildUS({'spec': spec, 'phi': phi, 'p': p}, ρ = ρ, commonX = commonX)
            m.calibrate()
            base, hbarRef = baselineRefs(m)
            led = m.solveLeaded(pinAtT0 = False)
            dates = m.db['dates']
            for pos, t in enumerate(m.db['t']):
                if pos + 1 >= len(m.db['t']):
                    break
                r = readout(m, led['τ'], led['report'], hbarRef, pos)
                rb = readout(m, base['τ'], base['report'], hbarRef, pos)
                rows.append({'spec': spec, 'phi': phi, 'commonX': commonX, 'p': p, 'pos': pos,
                             'date': dates[pos] if pos < len(dates) else np.nan,
                             'ν': float(m.db['ν'].xs(t)),
                             'θ': float(led['θ'].xs(t)), 'τ': r['τ'], 'sr': r['sr'],
                             'workweek': r['workweek'], 'R': r['R'],
                             'τExo': rb['τ'], 'srExo': rb['sr'], 'workweekExo': rb['workweek']}
                            | piColumn(m, t) | multiplicityColumns(led['multiplicity']) | thetaColumns())
            mergeWrite(out, rows, KEYPATH)
            print('[{}, φ={}] θ path: {}'.format(spec, phi,
                  '  '.join('{:.3f}'.format(x) for x in led['θ'].values[:8])))
    return pd.DataFrame(rows)


def stagePathHost(specs, phis, calib, country, out, ρ = 1.0, commonX = False):
    """ stagePath's design path for France and the UK, each on its own calibration and income groups (buildEU),
    under the US-calibrated cost ('US', escCalibration.csv) and under the country's own ('own',
    escCountry.csv), where one exists: at its workbook design theta = 1 France has none, since no value re-elects
    it (under a θ variant France's design is countryθ0's and may have one). The workweek is
    the country's own (its datum at the exogenous baseline in 2020); the *Exo columns pin the design at the
    country's observed theta*. `calib`, `country`: the two csvs' frames. """
    rows = []
    for spec in specs:
        for phi in phis:
            hitUS = pickCalib(calib, spec, phi, commonX)
            if hitUS.empty:
                continue
            for label in ('UK', 'FR'):
                hitOwn = pickOwn(country, label, spec, phi, commonX)
                wedges = [('US', float(hitUS.iloc[0]['p']))] + ([('own', float(hitOwn.iloc[0]['p']))]
                                                               if not hitOwn.empty else [])
                for wedgeFrom, p in wedges:
                    m = buildEU(label, {'spec': spec, 'phi': phi, 'p': p}, ρ = ρ, commonX = commonX)
                    m.calibrate()
                    base, hbarRef = baselineRefs(m)
                    led = m.solveLeaded(pinAtT0 = False)
                    dates = m.db['dates']
                    for pos, t in enumerate(m.db['t']):
                        if pos + 1 >= len(m.db['t']):
                            break
                        r = readout(m, led['τ'], led['report'], hbarRef, pos)
                        rb = readout(m, base['τ'], base['report'], hbarRef, pos)
                        rows.append({'spec': spec, 'phi': phi, 'commonX': commonX, 'country': label,
                                     'wedgeFrom': wedgeFrom, 'p': p, 'pos': pos,
                                     'date': dates[pos] if pos < len(dates) else np.nan,
                                     'ν': float(m.db['ν'].xs(t)), 'θStar': float(m.db['θ'].xs(m.t0Year)),
                                     'ω': float(m.db['ω'].xs(m.t0Year)),
                                     'θ': float(led['θ'].xs(t)), 'τ': r['τ'], 'sr': r['sr'],
                                     'workweek': r['workweek'], 'R': r['R'],
                                     'τExo': rb['τ'], 'srExo': rb['sr'], 'workweekExo': rb['workweek']}
                                    | multiplicityColumns(led['multiplicity']) | thetaColumns())
                    mergeWrite(out, rows, KEYPATHHOST)
                    print('[{}, φ={}] {} under the {} cost p={:.4f}: θ path {}'.format(
                        spec, phi, label, wedgeFrom, p, '  '.join('{:.3f}'.format(x) for x in led['θ'].values[:8])))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- stage: counterfactuals

def leadedNewPath(m, hbarRef, apply = None, data = None, pin = False):
    """ One counterfactual as a NEW EQUILIBRIUM PATH: the shocked parameters hold over the whole horizon,
    the economy starts at its own steady state, and everything is read at t0 (2020). shocks.shockedCopy
    is the same construction on the exogenous-theta side.

    pin=True holds theta at the shocked model's own exogenous design for the whole horizon -- the
    exogenous-theta reading. pin=False lets the political choice bind from the first period the horizon
    allows, so theta_{t0} is itself an equilibrium outcome rather than an inherited datum. That is what
    makes a table at 2020 informative: under the old copy-from-t0 convention theta was pinned as history
    through 2020 and the two readings agreed there by construction, so the tables had to be read at 2050.

    Returns readouts at t0 and t0+1 (t0 is what the tables print; t0+1 is kept because the design path's
    drift is the ageing prediction) plus the design at t0-1, t0, t0+1 and the solve's multiplicity. """
    pos = m.db['t0']
    mt = deepcopy(m)
    # as shocks.shockedCopy: the warm starts hold the BASELINE's solution, so clearing them is what keeps
    # each experiment independent of which one ran before it.
    mt.x0, mt.LOG.x0, mt.CRRA.x0 = {}, {}, {}
    if apply is not None:
        apply(mt, data)
    ε = mt.db['eps'].values.astype(float)
    if pin:
        θ = mt.db['θ'].values.astype(float)
        pol = mt.solvePEE_LOG(θ = θ, ε = ε)
        τ, report, θPath, mult = pol['τ'], pol['report'], pd.Series(θ, index = mt.db['t']), pol.get('multiplicity')
    else:
        led = mt.solveLeaded(pinAtT0 = False)
        τ, report, θPath, mult = led['τ'], led['report'], led['θ'], led['multiplicity']
    return {'t0': readout(mt, τ, report, hbarRef, pos), 't1': readout(mt, τ, report, hbarRef, pos+1),
            'θ_': float(θPath.iloc[pos-1]), 'θ0': float(θPath.iloc[pos]),
            'θ1': float(θPath.iloc[pos+1]), 'm': mt, 'multiplicity': mult}


def stageShocks(specs, phis, calib, out, ρ = 1.0, commonX = False, scenarios = None):
    """ The counterfactuals on the US at the US-calibrated cost parameter (escCalibration.csv's frame,
    `calib`). `scenarios` runs a subset (rows of other scenarios on disk are kept by mergeWrite); France's
    own row is run with any French scenario. Every row carries host 'US'. """
    rows = []
    host = 'US'
    for spec in specs:
        for phi in phis:
            hit = pickCalib(calib, spec, phi, commonX)
            if hit.empty:
                print(f'  [{spec}, φ={phi}] no calibrated cost -- skipped')
                continue
            p = float(hit.iloc[0]['p'])
            m = buildUS({'spec': spec, 'phi': phi, 'p': p}, ρ = ρ, commonX = commonX)
            m.calibrate()
            _, hbarRef = baselineRefs(m)
            names = scenarios or ESC_SCENARIOS
            data = frenchData(m, commonX = commonX) if any(n.startswith('fr') for n in names) else None
            print(f'\n[{host}: {spec}, φ={phi}, p={p:.4f}] counterfactuals')
            for name in names:
                apply, d = (None, None) if name == 'baseline' else (SHOCKS_ESC[name][1], data)
                for pin in (True, False):
                    tic = time.time()
                    try:
                        r = leadedNewPath(m, hbarRef, apply, d, pin = pin)
                    except Exception as e:
                        noteFailure(f'{spec} φ={phi} {name} pin={pin}', e)
                        continue
                    rows.append({'spec': spec, 'phi': phi, 'commonX': commonX, 'host': host, 'p': p,
                                 'scenario': name,
                                 'θpinned': pin, 'θ_tm1': r['θ_'], 'θ_t0': r['θ0'], 'θ_t1': r['θ1'],
                                 'τ_t0': r['t0']['τ'], 'sr_t0': r['t0']['sr'], 'ww_t0': r['t0']['workweek'],
                                 'τ_t1': r['t1']['τ'], 'sr_t1': r['t1']['sr'], 'ww_t1': r['t1']['workweek']}
                                | piColumn(r['m'], r['m'].t0Year) | multiplicityColumns(r['multiplicity'])
                                | thetaColumns())
                    print('  {:<10} pin={:<5} θ_t0={:.4f}  τ_t0={:.4f} sr_t0={:.4f} ww_t0={:.2f}'
                          '  (θ_t1={:.4f}, {:.0f}s)'
                          .format(name, str(pin), r['θ0'], r['t0']['τ'], r['t0']['sr'],
                                  r['t0']['workweek'], r['θ1'], time.time()-tic))
                    mergeWrite(out, rows, KEYSHK)

            # France's own calibrated path under the same wedge -- the endpoint the French-characteristics
            # rows are read against. Exogenous theta only: France's design is a datum here, and a leaded
            # France would be a different exercise (that is stageCountry's question).
            if not any(n.startswith('fr') for n in names):
                continue
            try:
                f = franceRow(m, spec, phi, p, hbarRef, ρ = ρ, commonX = commonX)
                rows.append(f | {'host': host} | thetaColumns())
                print('  {:<10} pin={:<5} θ_t0={:.4f}  τ_t0={:.4f} sr_t0={:.4f} ww_t0={:.2f}'
                      .format('France', 'True', f['θ_t0'], f['τ_t0'], f['sr_t0'], f['ww_t0']))
                mergeWrite(out, rows, KEYSHK)
            except Exception as e:
                noteFailure(f'{spec} φ={phi} France', e)
    return pd.DataFrame(rows)


def franceRow(m, spec, phi, p, hbarRef, ρ = 1.0, commonX = False):
    """ France's own calibrated equilibrium at 2020 under the same wedge, in escShocks' row schema.

    Not a shock on the host model: France carries its own eta, X, mu, nu AND its own calibrated omega, so
    what separates it from the all-French-characteristics row is precisely what the observable
    characteristics do not explain. The workweek is on the host's scale (hbarRef, and the host's datum
    passed explicitly to readout) and is ModelFR's own hours target, so that cell reproduces the datum
    rather than predicting it. """
    mFR = buildEU('FR', {'spec': spec, 'phi': phi, 'p': p}, ρ = ρ, commonX = commonX)
    mFR.calibrate()
    pos = mFR.db['t0']
    out = mFR.solvePEE_LOG()
    ww = float(m.db['workweek'])
    r0 = readout(mFR, out['τ'], out['report'], hbarRef, pos, workweekData = ww)
    r1 = readout(mFR, out['τ'], out['report'], hbarRef, pos+1, workweekData = ww)
    θ = mFR.db['θ'].values.astype(float)
    return {'spec': spec, 'phi': phi, 'commonX': commonX, 'p': p, 'scenario': 'France',
            'θpinned': True,
            'θ_tm1': float(θ[pos-1]), 'θ_t0': float(θ[pos]), 'θ_t1': float(θ[pos+1]),
            'τ_t0': r0['τ'], 'sr_t0': r0['sr'], 'ww_t0': r0['workweek'],
            'τ_t1': r1['τ'], 'sr_t1': r1['sr'], 'ww_t1': r1['workweek']} | multiplicityColumns(out.get('multiplicity'))


def shockIncomeAndVoting(mt0, d):
    """ France's income distribution AND its voting profile together.

    Neither shocks.py nor the paper runs the pair, and the pair is what the cross-section actually
    confronts: the two move in OPPOSITE directions for the design (more equality pushes theta up, a
    flatter voting profile pushes it down), so whether the model implies a theta-inequality relation at
    all depends on their net effect -- the question figure 1.1 poses.

    theta is re-installed after both, for the reason shocks.shockFrenchAll gives: shockVoting's
    updateAuxPars re-derives it from France's eta and would undo the income shock's pin. """
    θ0 = float(mt0.db['θ'].xs(mt0.db['t'][0])) if d.get('θUS') is None else float(d['θUS'])
    sh.shockIncomeDistribution(mt0, d['ηFR'], θPin = θ0)
    sh.shockVoting(mt0, d['μFR'])
    mt0.db.update(mt0.adjPar('θ', θ0))
    return float(mt0.db['θ'].xs(mt0.db['t'][mt0.db['t0']]))


SHOCKS_ESC = dict(sh.SHOCKS) | {
    'frBoth': ('Income distribution + voting', shockIncomeAndVoting),
}

# The order the shocks stage runs and the tables read. 'frAll' (shocks.shockFrenchAll) is the far end of
# the "mostly US, partly France" scale and franceRow is its endpoint; 'frBoth' stays because it isolates
# the two characteristics that move the DESIGN in opposite directions.
ESC_SCENARIOS = ('baseline', 'mild', 'acute', 'frIncome', 'frLeisure', 'frVoting', 'frBoth', 'frAll')

def pickOwn(country, host, spec, phi, commonX):
    """ `host`'s OWN converged cost calibration from escCountry.csv's frame (`country`), or an empty frame. """
    hit = country[(country['country'] == host) & (country['wedgeFrom'] == 'own') & (country['spec'] == spec)
                  & (country['phi'] == phi) & (country['converged'].astype(str) == 'True')]
    return hit[hit['commonX'].astype(bool) == bool(commonX)]


def frenchData(m, commonX = False):
    """ France's characteristics in the form shocks.py wants (runShocksUS.frenchData, LOG only), France's
    households cut at the US percentiles. 'θUS' is the US design.

    commonX must match the US model being shocked, for the reason buildEU gives. Under a θ variant France is
    calibrated at its own θ0 (countryθ0); 'θUS' is still the U.S. one, which the income row holds. """
    parsFR, kwFR, datesFR, wwFR = testEU.load('FR')
    θ0FR = countryθ0('FR')
    mFR = ModelFR(pars = parsFR | {'ρ': 1., 'ω': 2.} | ({} if θ0FR is None else {'θ0': θ0FR}), usRef = _usRef(commonX),
                  commonX = commonX, **kwFR)
    mFR.db['dates'], mFR.db['workweek'] = datesFR, wwFR
    mFR.LOG.initGS(GS | {'n': GRIDS[THETA['grid']]['n']})
    mFR.calibrate()
    t0FR, t0US = mFR.db['t'][mFR.db['t0']], m.t0Year
    XbarFR = float((mFR.db['γi'].xs(t0FR)*mFR.db['Xi'].xs(t0FR)).sum())
    XbarUS = float((m.db['γi'].xs(t0US)*m.db['Xi'].xs(t0US)).sum())
    # muFR is the WORKBOOK's mu_j -- nj entries, including the zero-mass slot -- because shockVoting
    # installs it through adjPar('μj'). Passing the ni-length mu_i instead raises a shape error, which is
    # how this was caught. eta is the other way round: shockIncomeDistribution takes eta_i.
    ηFR = mFR.db['ηi'].xs(t0FR).values.astype(float)
    return {'ηFR': ηFR, 'μFR': parsFR['μj'],
            'ηScale': sh.ηLevel(m, ηFR),   # as runShocksUS.frenchData: the leisure row needs the income row's scale
            'xbarRatio': XbarFR/XbarUS, 'θUS': float(m.db['θ'].xs(t0US))}


# ---------------------------------------------------------------- stage: France and the UK

def usWedgeRow(m, spec, phi, commonX, label, pUS):
    """ The country stage's row of `m` (France or the UK built at the U.S. parameter pUS, calibrated): the design in
    force in 2020 on the freely simulated path (leadedDesignAtT0 on its own policies: under 'sizePiR' solveLeaded
    at its own fixed point in pi), the counts of that recursion, and the 2020 tax on the theta*-pinned path
    (solvePEE_LOG, at its own pi, which wedgeReadout's 'ΠR' then reads). """
    t0 = m.t0Year
    θStar = float(m.db['θ'].xs(t0))
    ch = m.leadedDesignAtT0()
    mult = m.ESC.lastMultiplicity
    τc = float(m.solvePEE_LOG()['τ'].xs(t0))
    return ({'spec': spec, 'phi': phi, 'commonX': commonX, 'country': label, 'wedgeFrom': 'US',
             'p': pUS, 'θStar': θStar, 'choice': ch, 'ω': float(m.db['ω'].xs(t0)), 'τ': τc}
            | wedgeReadout(m, τc) | multiplicityColumns(mult))


def stageCountry(specs, phis, calib, out, commonX = False):
    rows = []
    for spec in specs:
        for phi in phis:
            hit = pickCalib(calib, spec, phi, commonX)
            if hit.empty:
                continue
            pUS = float(hit.iloc[0]['p'])
            for country in ('FR', 'UK'):
                label = country
                # (i) the US-calibrated wedge, imposed
                try:
                    m = buildEU(country, {'spec': spec, 'phi': phi, 'p': pUS}, commonX = commonX)
                    m.calibrate()
                    rows.append(usWedgeRow(m, spec, phi, commonX, label, pUS) | thetaColumns())
                    θStar, ch = rows[-1]['θStar'], rows[-1]['choice']
                    print('[{}, φ={}] {:<5} US wedge p={:.4f}: θ*={:.4f} -> choice {:.4f}  (Vtilde={:.4f})'
                          .format(spec, phi, label, pUS, θStar, ch, rows[-1]['Vtilde']))
                except Exception as e:
                    noteFailure(f'{spec} φ={phi} {label} US wedge', e)
                # (ii) its own calibrated p
                try:
                    m = buildEU(country, {'spec': spec, 'phi': phi, 'p': pUS}, commonX = commonX)
                    rec = m.calibrateWedge(spec = spec, phi = phi, verbose = False)
                    own = (wedgeReadout(m, m.ESC.τAt(m.t0Year, rec['θ'])) if rec['converged']
                           else {'Vtilde': float(m.B.Vtilde(m.t0Year))})
                    rows.append({'spec': spec, 'phi': phi, 'commonX': commonX, 'country': label,
                                 'wedgeFrom': 'own',
                                 'p': rec['p'], 'θStar': rec['θ'], 'choice': rec['θ'] + rec['residual'],
                                 'converged': rec['converged'],
                                 'ω': float(m.db['ω'].xs(m.t0Year))} | own
                                # calibrateWedge's last recursion is the one at the calibrated p
                                | multiplicityColumns(m.ESC.lastMultiplicity if rec['converged'] else None)
                                | thetaColumns())
                    print('[{}, φ={}] {:<5} own wedge: p={} θ*={:.4f} ({})'
                          .format(spec, phi, label, rec['p'], rec['θ'], rec['message']))
                except Exception as e:
                    noteFailure(f'{spec} φ={phi} {label} own wedge', e)
                mergeWrite(out, rows, KEYCTY)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- stage: the slope at the Bismarckian corner

CORNERH = 1e-4                       # the step of the one-sided difference at theta_{t+1} = 1
CORNERλ = (0.0, None, 100.0)         # None: the U.S. lambda of escCalibration


def cornerSlope(m, vote):
    """ d W_t/d theta_{t+1} at theta_{t+1} = 1 for the electorate at db['t'] position `vote` (one-sided,
    second order, step CORNERH), the design in force and its tax at t given, theta_{t+2} = 1; and the
    anticipated tax there, T_{t+1}(1). Under log the choice is separable (prop:esc:separability), so the
    continuation at theta_{t+2} = 1 does not move the slope. """
    tIdx = m.db['t']
    t, t1, t2 = tIdx[vote], tIdx[vote + 1], tIdx[vote + 2]
    tLag = tIdx[vote - 1] if vote > 0 else m.B.tFirst
    θt = float(m.db['θ'].xs(t))
    θ1 = 1.0 - CORNERH*np.arange(3)
    τt = m.ESC.τOfθ(t, np.array([θt]), tLag, terminal = False)[0][0]
    τ1 = m.ESC.τOfθ(t1, θ1, t, terminal = False)[0]
    τ2 = m.ESC.τOfθ(t2, np.array([1.0]), t1, terminal = (t2 == tIdx[-1]))[0][0]
    cont = {'τ1': τ1, 'θ2': np.ones_like(θ1), 'τ2': np.full_like(θ1, τ2), 'terminal1': False}
    W, _ = m.ESC.objective(t, tLag, t1, np.full_like(θ1, τt), np.full_like(θ1, θt), θ1, cont)
    return {'date': m.db['dates'][vote], 'ω': float(m.db['ω'].xs(t)),
            'τ1Corner': float(τ1[0]), 'slope': float((3*W[0] - 4*W[1] + W[2])/(2*CORNERH))}


def stageCorner(specs, phis, calib, out, commonX = False):
    rows = []
    for spec in [s for s in specs if s == 'size']:
        for phi in phis:
            hit = pickCalib(calib, spec, phi, commonX)
            if hit.empty:
                continue
            for country in ('FR', 'US', 'UK'):
                for λ in CORNERλ:
                    λv = float(hit.iloc[0]['p']) if λ is None else λ
                    wedge = None if λv == 0 else {'spec': spec, 'phi': phi, 'p': λv}
                    m = (buildUS(wedge, commonX = commonX) if country == 'US'
                         else buildEU(country, wedge, commonX = commonX))
                    try:
                        m.calibrate()
                    except RuntimeError as e:                 # a calibration that fails at this λ: no row
                        print('[{}] {} λ={:8.3f}: calibration FAILED, no row ({})'.format(spec, country, λv, str(e)[:80]))
                        continue
                    p0 = int(m.db['t'].get_loc(m.t0Year))
                    for vote in (p0 - 1, p0):
                        r = cornerSlope(m, vote)
                        rows.append({'spec': spec, 'phi': phi, 'commonX': commonX, 'country': country,
                                     'λ': λv, 'λFrom': 'US' if λ is None else 'fixed', 'vote': vote,
                                     'θStar': float(m.db['θ'].xs(m.t0Year))} | r | thetaColumns())
                        print('[{}] {} λ={:8.3f} vote {}: slope at θ=1 {:+.6f}, T(1) {:.6f}'
                              .format(spec, country, λv, r['date'], r['slope'], r['τ1Corner']))
                    mergeWrite(out, rows, KEYCORNER)
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description = 'Endogenous theta: leaded choice under the deadweight wedge.')
    p.add_argument('--stage', nargs = '*', default = ['calib', 'path', 'shocks', 'country'],
                   choices = ('calib', 'path', 'pathHost', 'shocks', 'country', 'corner'),
                   help = "pathHost: the design path of France and the UK under the US and their own cost "
                          "(escPathHost.csv; needs escCalibration.csv and escCountry.csv)")
    p.add_argument('--scenarios', nargs = '*', default = None,
                   help = 'the shocks stage runs these scenarios only (default: all); rows of the others '
                          'on disk are kept')
    p.add_argument('--spec', nargs = '*', default = ['scale', 'flat'],
                   choices = ('scale', 'flat', 'size', 'sizePiR'),
                   help = "'size' is the paper's wedge (f(theta, tau), lambda in the p slot); 'scale' the "
                          "previous one and the comparison arm; 'flat' MGE's variant; 'sizePiR' the paper's "
                          "with lambda*Pi_t/R_t (A8's Pi-inclusive cost)")
    p.add_argument('--phi', type = float, nargs = '*', default = [0.5, 0.25, 0.75])
    p.add_argument('--rho', type = float, default = 1.0)
    p.add_argument('--commonX', action = 'store_true',
                   help = 'the common-X calibration variant (the paper leads with it). Recorded as a '
                          'column in every csv written here and part of the merge key, so the two '
                          'variants coexist in one file instead of overwriting each other.')
    p.add_argument('--bracket', type = float, nargs = 2, default = None,
                   help = 'the scan bracket for the wedge parameter (p, or lambda under size); default per '
                          'spec, ModelESC.WEDGE_BRACKET. Applies to every spec asked for in one run.')
    p.add_argument('--tag', default = '')
    p.add_argument('--theta', default = None,
                   help = "a variant of data/oecd26_theta.csv: every U.S. model at that θ (db['θ0']); stages "
                          + ', '.join(THETASTAGES) + "; 'index' the workbook's θ into results/esc/ (required for it)")
    p.add_argument('--root', default = None,
                   help = "the variant's directory (default results/oecd26_theta_<variant>): reads calibration/, "
                          "writes esc/. A relative path is the working directory's.")
    p.add_argument('--coarse', action = 'store_true', help = "GRIDS['coarse'], under --theta only")
    p.add_argument('--main', action = 'store_true',
                   help = "with --theta config.US['theta']: read results/calibration/, write results/esc/")
    a = p.parse_args()
    a.root = None if a.root is None else os.path.abspath(os.path.join(CWD, a.root))
    import oecd26
    theta = a.theta if a.root is not None else oecd26.indexTheta(a.theta, [OUTDIR], 'runESC.py')
    outdir = setTheta(theta, a.root, a.stage, a.coarse, a.commonX, a.rho, a.main)
    os.makedirs(outdir, exist_ok = True)
    f = lambda n: os.path.join(outdir, f'{n}{a.tag}.csv')

    calib = None
    if 'calib' in a.stage:
        calib = stageCalib(a.spec, a.phi, f('escCalibration'), ρ = a.rho, commonX = a.commonX,
                           bracket = None if a.bracket is None else tuple(a.bracket))
    else:
        calib = checkCalib(pd.read_csv(f('escCalibration')), f('escCalibration'), a.spec, a.phi, a.commonX)
    if 'path' in a.stage:
        stagePath(a.spec, a.phi, calib, f('escPath'), ρ = a.rho, commonX = a.commonX)
    if 'pathHost' in a.stage:
        cty = pd.read_csv(f('escCountry'))
        cty = checkCalib(cty[cty['country'].isin(['UK', 'FR'])], f('escCountry'), a.spec, a.phi, a.commonX)
        stagePathHost(a.spec, a.phi, calib, cty, f('escPathHost'), ρ = a.rho, commonX = a.commonX)
    if 'shocks' in a.stage:
        stageShocks(a.spec, a.phi, calib, f('escShocks'), ρ = a.rho, commonX = a.commonX,
                    scenarios = a.scenarios)
    if 'country' in a.stage:
        stageCountry(a.spec, a.phi, calib, f('escCountry'), commonX = a.commonX)
    if 'corner' in a.stage:
        stageCorner(a.spec, a.phi, calib, f('escCornerSlope'), commonX = a.commonX)
    print('\n-> {}'.format(os.path.relpath(outdir, REPO)))
    if FAILURES:
        print('{} FAILED:'.format(len(FAILURES)) + ''.join('\n  ' + x for x in FAILURES))
    return 1 if FAILURES else 0


THETASTAGES = ('calib', 'path', 'pathHost', 'shocks', 'country', 'corner')     # the stages a θ variant runs


def setTheta(theta = None, root = None, stages = ('calib',), coarse = False, commonX = True, ρ = 1.0, main = False):
    """ Install a run's θ variant in THETA and return its output directory: results/esc without one (THETA left at
    its defaults); under VARIANT, root/esc with root the variant's directory (oecd26.thetaDir by default; a
    relative path is the repo's; results/ under `main`, oecd26.thetaRoot), θ0 its U.S. θ and calibDir
    root/calibration, whose U.S. sweep at ρ must be at that θ0 (oecd26.usBeta). SystemExit for --root, --coarse or
    --main without a variant, --main with --root or with another variant than config.US['theta'], a stage outside
    THETASTAGES, a variant without a U.S. θ, a missing or other-θ0 sweep, or results/esc/ as the output without
    --main. """
    if theta is None:
        if root is not None or coarse or main:
            raise SystemExit('--root, --coarse and --main need --theta VARIANT')
        return OUTDIR
    import oecd26
    if set(stages) - set(THETASTAGES):
        raise SystemExit('--theta {}: stages {} only, not {}'.format(theta, THETASTAGES,
                                                                    sorted(set(stages) - set(THETASTAGES))))
    if main and root is not None:
        raise SystemExit('--main reads results/calibration/ and writes results/esc/: not with --root')
    try:
        θ0 = oecd26.thetaFromFile('USA', theta)
        root = oecd26.thetaRoot(theta, main) if root is None else os.path.join(REPO, root)
    except (ValueError, FileNotFoundError) as e:          # NoTheta is a ValueError
        raise SystemExit('--theta {}: {}'.format(theta, e))
    outdir = os.path.join(root, 'esc')
    if not main and os.path.normcase(os.path.abspath(outdir)) == os.path.normcase(OUTDIR):
        raise SystemExit('--theta {}: results/esc/ holds the main variant\'s runs (--main) only; use another '
                         '--root'.format(theta))
    calibDir = os.path.join(root, 'calibration')
    try:
        oecd26.usBeta(os.path.join(calibDir, 'US_rhoGridCommonX.csv' if commonX else 'US_rhoGrid.csv'), ρ,
                      thetaSource = theta)
    except (FileNotFoundError, ValueError, KeyError) as e:
        raise SystemExit('--theta {}: {}'.format(theta, e))
    THETA.update({'variant': theta, 'θ0': θ0, 'calibDir': calibDir, 'grid': 'coarse' if coarse else 'fine',
                  'main': bool(main)})
    print('θ variant {!r}: θ0 = {!r}, U.S. sweep {}, grid {}'.format(theta, θ0, usSweep(commonX), THETA['grid']))
    return outdir


def checkCalib(calib, path, specs = None, phis = None, commonX = None):
    """ escCalibration.csv's (or escCountry.csv's) frame, refused under a θ variant unless every row a stage reads --
    spec in `specs`, phi in `phis`, the run's commonX (all rows where None) -- is at THETA's thetaSource, θ0 and grid
    (#13: the cost a later stage imposes was calibrated at that θ). Rows of other keys, which mergeWrite keeps beside
    a run's own, are not read and not checked. """
    want = thetaColumns()
    if not want:
        return calib
    read = pd.Series(True, index = calib.index)
    if specs is not None:
        read &= calib['spec'].isin(list(specs))
    if phis is not None:
        read &= np.isclose(pd.to_numeric(calib['phi'], errors = 'coerce').values[:, None],
                           np.asarray(list(phis), dtype = float)[None, :]).any(axis = 1)
    if commonX is not None and 'commonX' in calib.columns:
        read &= calib['commonX'].astype(bool) == bool(commonX)
    rows = calib[read]
    ok = all(c in calib.columns for c in want) and bool(
        ((rows['thetaSource'] == want['thetaSource']) & (rows['grid'] == want['grid'])
         & (pd.to_numeric(rows['θ0']) - want['θ0']).abs().le(1e-12)).all())
    if not ok:
        raise SystemExit('{}: rows at another thetaSource, θ0 or grid than {}: rerun --stage calib (or country)'
                         .format(path, want))
    return calib


if __name__ == '__main__':
    sys.exit(main())

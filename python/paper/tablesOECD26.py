r""" The cross-country test of section 7's model (python/US/runOECD26.py, results/oecd26/): the panel data shared
with figuresOECD26, the slopes of figure 1.1's relations in the data and in the model (all five, OECD26_Slopes; the
two of section 6's spending test at the observed designs without the cost, spendingSlopes), the calibrated political
weights (OECD26_Omega), the five slopes under the supplemental appendix's population and data variants (OECD26_Variants,
page C4), and the predicted against the observed design (thetaFit, figuresOECD26's OECD26_ThetaFit*).

NAMES: layer 1 (layer1.csv, `layer` 1) is printed as "the model" (header `Model`), layer 2 (layer2.csv, each country's
own political weight) as "own ω" (`Own $\omega$`); no output prints "layer". OECD26_Slopes at the paper's variant
shows the model alone (table 5); own ω's slopes there are OECD26_Variants' starred rows.

THE DESIGN'S MEASURE (`theta`, config.US['theta'] by default): the data side's index is I = 1/(2 - θ) of that variant
of data/oecd26_theta.csv and its spending the variant's, on its countries (datasets.oecd26SampleOf: inSample and a
finite θ); the model's rows are the variant's (datasets.oecd26Dir). `theta` 'paper' is figure 1.1's index on
data/oecd26_sample.csv's sample26, with the 'index' rows of results/oecd26/.

Every builder takes
  turnout  'csesFlat' (each economy at its own CSES propensities to vote; Luxembourg, which has none and enforces
           compulsory voting, flat) | 'cses' (Luxembourg left out of data and model alike);
  sample   '26', the sample of the design's measure (above);
  reading  'inForce' (the design in force in 2020 on the freely simulated path, chosen in 1990: the reading of the
           cost parameter's target) | 'chosen' (the design chosen in 2020 at the observed design, in force from 2050);
  spec     the cost transported from the U.S. anchor: 'sizePiR' (λ_n = ε Π_n/R_n for country n, the anchor's ε in
           the csvs' `p`; under config.CARRIER 'convex' λ_n = 2κ π_n(1 - π_n), π_n = Π_n/R_n, κ in `p`), 'size' (the anchor's λ) | 'none' (no cost; the model then reads no design);
  betaRule 'common' (the anchor's β in every country) | 'commonR' (each country's β set so that its 2020 return is
           the anchor's);
  spending the model's spending: 'chosen' (100 τInForce (1 - α), the 2020 tax on the path whose design in force is the
           chosen one, = spending2020 τInForce/τ2020; reading 'inForce' and a cost only) | 'observed' (spending2020, the
           tax at the observed design, which own ω (layer 2) matches by construction); under a cost and the spending
           target 'net' (config.US['esc']['spendingTarget'], netSpending) each is the benefits paid, 100 τ f (1 - α):
           spendingNetInForce and spendingNet2020;
  omegaRule the model's (layer 1's) ω: 'mean' (the mean of the countries' own ω) | 'anchor' (the U.S. anchor's);
  retAge   ν's retirement-age correction: 'nra' (statutory age) | 'none' | 'exit' (effective exit age);
  theta    the design's measure (config.THETANAMES): config.US['theta'] or a data variant (results/oecd26_theta_<V>/);
each None for the paper's variant, PAPER, which is config.US['oecd26']'s (sample '26', spending its
`spendingReading`, theta config.US['theta']; under spec 'none' spending is 'observed'). An output's name and label
carry the parts that differ from the paper's variant (`suffix`: '_cses', '_omegaUS', '_nuNone', '_thetaB', ...). The
model's design enters every regression on the index's scale, I = 1/(2 - θ),
the inverse of the loader's θ = 2 - 1/I; spending is 100 τ (1 - α), % of GDP, the inverse of the calibration's τ0. Both
sides of every slope use figure 1.1's regressors (OECD population growth 2020/1990, the measure of income inequality
config.INEQUALITY, the role 'inequality': y_3/y_1 of the income groups the model is fed, datasets.incomeGroups) and the
same countries. Page D1's inequalitySlopes set figure 7.2's panels on inequality on the Gini and both model measures.
"""
import os, re
import numpy as np, pandas as pd

import config as C
import datasets as D
import oecdFigure1 as F1
from tables import rowKey, keyed
from tablesUS import _xwrap

# Figure 1.1's relations, (y role, x role), in the order the tables print them: the four on the exogenous
# dimensions first, spending on the index last. The x role 'index' is the design on the index scale.
RELATIONS = (('spending', 'growth'), ('spending', 'inequality'), ('index', 'growth'), ('index', 'inequality'),
             ('spending', 'index'))
# Section 6's spending test (spendingSlopes, figuresOECD26.oecd26Spending), a subset of RELATIONS.
SPENDING = (('spending', 'growth'), ('spending', 'inequality'))
READINGS = {'chosen': 'θChosen', 'inForce': 'θInForce'}
SPENDINGS = ('chosen', 'observed')
TURNOUTS = ('csesFlat', 'cses')
OMEGARULES = ('mean', 'anchor')
RETAGES = ('nra', 'none', 'exit')
# The measures of the design: the paper's, then the supplemental appendix's data variants (page C4).
MEASURES = (C.US['theta'],) + tuple(v for v in C.US['oecd26']['thetaVariants'] if v != C.US['theta'])
VARIANT = ('turnout', 'sample', 'reading', 'spec', 'betaRule', 'spending', 'omegaRule', 'retAge', 'theta')  # suffix's order
PAPER = ({k: C.US['oecd26'][k] for k in ('turnout', 'reading', 'spec', 'betaRule', 'omegaRule', 'retAge')}
         | {'sample': '26', 'spending': C.US['oecd26']['spendingReading'], 'theta': C.US['theta']})
# suffix's part for a value that alone would be ambiguous ('none' is a spec and a retAge) or opaque.
SUFFIXPART = {'omegaRule': {'anchor': 'omegaUS', 'mean': 'omegaMean'},
              'retAge': {'none': 'nuNone', 'exit': 'nuExit', 'nra': 'nuNra'}}
YTEX = {'spending': 'Pension spending', 'index': 'Bismarckian index'}
XTEX = {'growth': 'population growth', 'inequality': C.INEQUALITYLABELS[C.INEQUALITY]['text'],
        'index': 'the Bismarckian index'}
# Page D1's measures of income inequality, in the order inequalitySlopes and figuresOECD26.oecd26Layer1Inequality
# print them: the Gini first, the comparison, then the model's measures (panelData's columns of the same names).
INEQUALITYROWS = ('gini', 'yHighOverLow', 'Vtilde')
# The data files every slope reads, for the banners: figure 1.1's, the design's measure and the income groups.
DATASRC = 'data/oecdFigure1.csv, oecd26_theta.csv, oecd29_distribution.csv'
# What a caption adds for a cost, β rule, turnout or spending other than the paper's.
SPECCAPTION = {'size': r'the U.S.\ $\lambda$ in every country',
               'sizePiR': r'the U.S.\ $' + C.PIRSYMBOL + r'$ in every country, ' + C.PIRCOSTN,
               'none': 'without the cost'}
BETACAPTION = {'common': r'the U.S.\ $\beta$ in every country',
               'commonR': r'each country\textquotesingle s $\beta$ set to the U.S.\ 2020 return'}
TURNOUTCAPTION = {'csesFlat': 'each country at its own turnout', 'cses': 'without Luxembourg'}
SPENDINGCAPTION = {'chosen': 'spending at the chosen designs', 'observed': 'spending at the observed designs'}
OMEGACAPTION = {'anchor': r'every country at the U.S.\ $\omega$', 'mean': r'every country at the mean $\omega$'}
RETAGECAPTION = {'none': r'$\nu$ uncorrected', 'exit': r'$\nu$ corrected at the effective exit age',
                 'nra': r'$\nu$ corrected at the normal retirement age'}
THETACAPTION = {m: r'the design measured as ' + C.THETANAMES.get(m, m) for m in MEASURES}
# A layer as a header or row label prints it (NAMES above).
LAYERNAME = {1: 'Model', 2: r'Own $\omega$'}
# Section 7.3's two costs, by the power k of the tax rate in f (Eq esc:budget at k = 1, esc:budgetSquared at k = 2): k = 2 in
# datasets.oecd26Dir(), k = 1 in datasets.tau1Dir('oecd26') under config.US['tau1']'s carrier (costKw). Table 5 and
# figures 7.2 and 7.3 show the model under both at the paper's variant; every other exhibit is k = 2's.
COSTS = (1, 2)
COSTKEY = C.COSTKEYS                                                   # mark series and `% row:` parts
COSTNAME = {1: 'linear cost', 2: 'squared cost'}                         # legends and table headers, after "model, "
COSTTITLE = {1: 'Cost linear in the tax rate', 2: 'Cost in the tax rate squared'}       # panel titles
COSTTEX = {1: 'the cost linear in the tax rate', 2: 'the cost in the tax rate squared'}


def costKw(k):
    """ panelData's `root` and `carrier` of the cost with tau^k: {} at k = 2 (the defaults), config.US['tau1']'s at 1. """
    if k == 2:
        return {}
    t = C.US['tau1']
    if k != t['tauPower']:
        raise ValueError('cost k = {}: one of {}'.format(k, COSTS))
    return {'root': D.tau1Dir('oecd26'), 'carrier': t['carrier']}


def costsSrc(layers = 'layer1.csv, anchor.csv'):
    """ The %% Source: of an exhibit of both costs. """
    return 'results/' + C.US['tau1']['oecd26'] + '/, ' + D.oecd26DirRel() + ' ' + layers + '; ' + DATASRC


# The notes' sentence on the propensities to vote by income group.
TURNOUTNOTE = {'csesFlat': r'Turnout: each country\textquotesingle s own propensities to vote by income group, from the '
                           r'CSES; Luxembourg, which has no CSES survey and enforces compulsory voting, at a flat profile.',
               'cses': r'Turnout: each country\textquotesingle s own propensities to vote by income group, from the CSES.'}


# What the notes of the paper's variant say in words (its cost \varepsilon, own turnout with Luxembourg flat, the design in
# force in 2020 chosen in 1990, spending on the path of that design) rather than from PAPER.
PAPERWORDING = {'spec': 'sizePiR', 'turnout': 'csesFlat', 'reading': 'inForce', 'spending': 'chosen'}


def paperWording():
    """ Raises unless config.US['oecd26'] is the variant PAPERWORDING describes, the one the paper's notes (figures
    OECD26_Layer1 and OECD26_ThetaFit, table OECD26_Omega) state in words. """
    bad = {k: PAPER[k] for k, v in PAPERWORDING.items() if PAPER[k] != v}
    if bad:
        raise ValueError("the notes of the paper's cross-country test are worded for {}; config has {}: reword them"
                         .format(PAPERWORDING, bad))


def turnoutNote(turnout):
    """ The notes' sentence on the propensities to vote by income group under `turnout`. """
    return TURNOUTNOTE[turnout]
COSTSYMBOL = {'size': r'$\lambda$', 'sizePiR': '$' + C.PIRSYMBOL + '$'}         # the anchor's `p` (config.CARRIER)
# The design's measure as the notes state it: config.US['theta']'s (config.THETADEF), and the index on figure 1.1's
# scale.
MEASURENOTE = r'$\theta$, ' + C.THETADEF + r', and $I = 1/(2-\theta)$, the Bismarckian index.'


def measureNote(v):
    r""" The notes' sentence on θ and the index at the variant `v`: config.US['theta']'s measure in full (MEASURENOTE),
    a data variant by its name (config.THETANAMES), which the supplemental appendix defines (page C4's section
    `esc-oecd26-measure`, cited as \oa: the paper has no label for it). """
    if v['theta'] == C.US['theta']:
        return MEASURENOTE
    return (r'$\theta$ measured as ' + C.THETANAMES.get(v['theta'], v['theta']) + r' (\oa{esc-oecd26-measure}), '
            r'and $I = 1/(2-\theta)$.')


def variant(turnout = None, sample = None, reading = None, spec = None, betaRule = None, spending = None,
            omegaRule = None, retAge = None, theta = None):
    """ The variant as a dict over VARIANT, PAPER's value where None; spending 'observed' under spec 'none', which
    chooses no design. """
    v = dict(zip(VARIANT, (turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta)))
    v = {k: PAPER[k] if v[k] is None else v[k] for k in VARIANT}
    if spending is None and v['spec'] == 'none':
        v['spending'] = 'observed'
    return v


def suffix(turnout = None, sample = None, reading = None, spec = None, betaRule = None, spending = None,
           omegaRule = None, retAge = None, theta = None):
    r""" '' for the paper's variant (PAPER); otherwise the parts that differ, e.g. '_cses', '_size_commonR',
    '_nuExit', '_thetaB' (SUFFIXPART, 'theta' + the measure). Appended to an output's name and \label. Under spec
    'none' the spending reading is not a part. """
    v = variant(turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta)
    part = lambda k: 'theta' + v[k] if k == 'theta' else SUFFIXPART.get(k, {}).get(v[k], v[k])
    return ''.join('_' + part(k) for k in VARIANT
                   if v[k] != PAPER[k] and not (k == 'spending' and v['spec'] == 'none'))


def variantCaption(turnout = None, sample = None, reading = None, spec = None, betaRule = None, spending = None,
                   omegaRule = None, retAge = None, theta = None):
    """ What a caption adds for a variant other than the paper's: '' for the paper's own. """
    v = variant(turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta)
    parts = ([r'every country of figure \ref{fig:US:OECDdata}'] if v['sample'] != PAPER['sample'] else []) \
        + ([TURNOUTCAPTION[v['turnout']]] if v['turnout'] != PAPER['turnout'] else []) \
        + (['the design chosen in 2020'] if v['reading'] != PAPER['reading'] else []) \
        + ([SPECCAPTION[v['spec']]] if v['spec'] != PAPER['spec'] else []) \
        + ([BETACAPTION[v['betaRule']]] if v['betaRule'] != PAPER['betaRule'] else []) \
        + ([SPENDINGCAPTION[v['spending']]] if v['spending'] != PAPER['spending'] and v['spec'] != 'none' else []) \
        + ([OMEGACAPTION[v['omegaRule']]] if v['omegaRule'] != PAPER['omegaRule'] else []) \
        + ([RETAGECAPTION[v['retAge']]] if v['retAge'] != PAPER['retAge'] else []) \
        + ([THETACAPTION[v['theta']]] if v['theta'] != PAPER['theta'] else [])
    return ''.join(', ' + p for p in parts)


def _thetaSource(v):
    """ The thetaSource of the variant's model rows: 'index' for figure 1.1's index ('paper'), else the measure. """
    return 'index' if v['theta'] == 'paper' else v['theta']


def _check(v):
    if (v['turnout'] not in TURNOUTS or v['sample'] != '26' or v['reading'] not in READINGS
            or v['spec'] not in SPECCAPTION or v['betaRule'] not in BETACAPTION or v['spending'] not in SPENDINGS
            or v['omegaRule'] not in OMEGARULES or v['retAge'] not in RETAGES or v['theta'] not in MEASURES + ('paper',)):
        raise ValueError('turnout {}, sample 26, reading chosen|inForce, spec size|sizePiR|none, betaRule common|commonR, '
                         'spending chosen|observed, omegaRule {}, retAge {}, theta {}; got {}'.format(
                             '|'.join(TURNOUTS), '|'.join(OMEGARULES), '|'.join(RETAGES), '|'.join(MEASURES), v))
    if v['spending'] == 'chosen' and (v['reading'] != 'inForce' or v['spec'] == 'none'):
        raise ValueError("spending 'chosen' is the tax on the path of the design in force in 2020: reading 'inForce' "
                         "and a cost only; got {}".format(v))


def netSpending(v):
    """ Is the model's spending of variant `v` the benefits paid, 100 τ f (1 - α) (`spendingNet2020`,
    `spendingNetInForce` of the rows), rather than 100 τ (1 - α)? Under a cost and config.US['esc']['spendingTarget']
    'net', where every ω is calibrated to it; without a cost f = 1 and the two coincide. """
    return v['spec'] != 'none' and C.US['esc']['spendingTarget'] == 'net'


def spendingTex(v):
    r""" The model's spending as the notes print it at variant `v`: the benefits paid under netSpending, the tax
    times labour's share otherwise. """
    return r'$\tau f(1-\alpha)$, the benefits paid,' if netSpending(v) else r'$\tau(1-\alpha)$'


def dataSide(theta = None):
    """ (countries, index, spending) of the data side under the design's measure `theta` (module docstring): figure
    1.1's order, I and spending (% of GDP) as Series on ISO3. """
    theta = C.US['theta'] if theta is None else theta
    fig = D.oecdFigure1Data()
    if theta == 'paper':
        smp = D.oecd26Sample()
        return ([c for c in fig.index if bool(smp.at[c, 'sample26'])], fig[F1.NAMED['index']].astype(float),
                fig[F1.NAMED['spending']].astype(float))
    t = D.thetaVariant(theta)
    return D.oecd26SampleOf(theta), t['index'].astype(float), t['spending'].astype(float)


def panelData(layer, turnout = None, sample = None, reading = None, spec = None, betaRule = None, theta = None,
              spending = None, omegaRule = None, root = None, retAge = None, inequality = None, usInput = None,
              carrier = None):
    r""" One row per country of the sample with a model row: figure 1.1's x variables (`growth`, `inequality` = the
    measure `inequality` of INEQUALITYROWS, config.INEQUALITY's where None), every measure of INEQUALITYROWS (`gini`,
    `yHighOverLow`, `Vtilde`), the data's
    `spendingData`, `indexData` and `θData` = 2 - 1/I under the design's measure `theta` (dataSide), and the model's
    `spendingModel` (at the variant's spending reading), `spendingObserved` (spending2020), `θModel`, `indexModel` =
    1/(2 - θ), from layer `layer` (1 | 2) of that measure at the variant's retAge; layer 1 at the variant's omegaRule
    ('mean' | 'anchor'), both from datasets.oecd26Dir(theta, root). Countries whose model row is 'noTurnout' (cses) are
    dropped from both sides and listed in `.attrs['dropped']`; a sample country with no row is MissingInput, one
    without its three income groups (datasets.incomeGroups) or a failed row an error, a model row whose Vtilde is not
    its income groups' StaleInput. `usInput` (layer 1): the rows with that input of every country at its U.S. value
    (python/US/runOECD26.py --usInput; under 'income' every row's Vtilde is the U.S. groups'); `carrier`: the carrier the
    anchor of `root` carries (datasets.oecd26Anchor; config.CARRIER where None). """
    v = variant(turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta)
    _check(v)
    theta = v['theta']
    fig = D.oecdFigure1Data()
    smp = D.oecd26Sample()
    model = D.oecd26Layer(layer, v['turnout'], v['spec'], v['betaRule'], _thetaSource(v), root,
                          v['omegaRule'] if layer == 1 else 'anchor', v['retAge'], usInput if layer == 1 else None,
                          carrier)
    isos, index, spend = dataSide(theta)
    lack = [c for c in isos if c not in model.index]
    if lack:
        raise D.MissingInput('{} (layer {}, {} {} {}: no row for {})'.format(
            D.oecd26Path(layer, None if theta == 'paper' else theta, root), layer, v['turnout'], v['spec'],
            v['betaRule'], ', '.join(lack)))
    m = model.loc[isos]
    bad = m.index[~m['status'].isin(['ok', 'noTurnout'])]
    if len(bad):
        raise ValueError('layer {} ({}): failed rows for {}'.format(layer, v['turnout'], ', '.join(bad)))
    moved = [c for c in isos if bool(m.at[c, 'sample26']) != bool(smp.at[c, 'sample26'])]
    if moved:
        raise ValueError('layer {} ({}): sample26 differs from data/oecd26_sample.csv for {}: rerun the layer'
                         .format(layer, v['turnout'], ', '.join(moved)))
    groups = D.incomeGroups(isos)                                  # every sample country has its three groups
    dropped = list(m.index[m['status'] == 'noTurnout'])
    keep = [c for c in isos if c not in dropped]
    ref = (pd.Series(float(D.incomeGroups(['USA']).loc['USA', 'Vtilde']), index = keep) if usInput == 'income'
           else groups.loc[keep, 'Vtilde'])
    gap = float((m.loc[keep, 'Vtilde'].astype(float) - ref).abs().max())
    if not gap <= D.GROUPTOL:
        raise D.StaleInput('{}: the rows\' Vtilde differs from the income groups of {} by {:.2e}: rows computed at other '
                           'income groups; rerun python/US/runOECD26.py'.format(
                               D.oecd26Path(layer, None if theta == 'paper' else theta, root),
                               os.path.relpath(D.OECD29DIST, C.REPO), gap))
    θ = m.loc[keep, READINGS[v['reading']]].astype(float)
    if netSpending(v):
        lack = [k for k in ('spendingNet2020', 'spendingNetInForce') if k not in m.columns]
        if lack:
            raise D.StaleInput('{}: no {} (rows written before the net spending target): rerun python/US/runOECD26.py'
                               .format(D.oecd26Path(layer, None if theta == 'paper' else theta, root), ', '.join(lack)))
        observed = m.loc[keep, 'spendingNet2020'].astype(float)
        model = observed if v['spending'] == 'observed' else m.loc[keep, 'spendingNetInForce'].astype(float)
    else:
        observed = m.loc[keep, 'spending2020'].astype(float)
        model = (observed if v['spending'] == 'observed' else
                 observed*m.loc[keep, 'τInForce'].astype(float)/m.loc[keep, 'τ2020'].astype(float))
    inequality = C.INEQUALITY if inequality is None else inequality
    if inequality not in INEQUALITYROWS:
        raise ValueError('inequality {!r}: one of {}'.format(inequality, INEQUALITYROWS))
    out = pd.DataFrame({'country': fig.loc[keep, 'country'],
                        'growth': fig.loc[keep, F1.NAMED['growth']].astype(float),
                        'gini': fig.loc[keep, 'gini'].astype(float),
                        'yHighOverLow': groups.loc[keep, 'yHighOverLow'], 'Vtilde': groups.loc[keep, 'Vtilde'],
                        'spendingData': spend.loc[keep], 'indexData': index.loc[keep],
                        'θData': 2 - 1/index.loc[keep],
                        'spendingModel': model, 'spendingObserved': observed,
                        'θModel': θ, 'indexModel': 1/(2 - θ)}, index = pd.Index(keep, name = 'iso3'))
    out.insert(2, 'inequality', out[inequality])
    out.attrs['dropped'] = dropped
    return out


def ols(x, y):
    """ Least squares of y on a constant and x: slope, its classical standard error, intercept, Pearson's r
    and n. """
    x, y = np.asarray(x, dtype = float), np.asarray(y, dtype = float)
    n = len(x)
    xc = x - x.mean()
    b = float(np.sum(xc*(y - y.mean()))/np.sum(xc**2))
    a = float(y.mean() - b*x.mean())
    e = y - a - b*x
    se = float(np.sqrt(np.sum(e**2)/(n - 2)/np.sum(xc**2)))
    return {'slope': b, 'se': se, 'intercept': a, 'r': float(np.corrcoef(x, y)[0, 1]), 'n': n}


def _xy(p, yRole, xRole, side):
    """ (x, y) of one relation on `side` 'Data' | 'Model': the regressors are the data's, except the index,
    which is each side's own. """
    y = p[yRole + side]
    x = p[xRole + side] if xRole == 'index' else p[xRole]
    return x, y


def slopes(turnout = None, sample = None, reading = None, spec = None, betaRule = None, relations = RELATIONS,
           layers = (1, 2), theta = None, spending = None, omegaRule = None, root = None, keep = None, retAge = None,
           inequality = None, usInput = None, carrier = None):
    """ One row per relation of `relations` (RELATIONS by default): OLS of the data and of each layer of `layers`
    (`data`, `layer1`, `layer2`, each an ols dict, None for a layer not asked for), on the countries every layer
    asked for covers (and among `keep`, if given), under the design's measure `theta` and at `retAge` (panelData;
    layer 1 at `omegaRule`, both from `root`; income inequality the measure `inequality`, config.INEQUALITY's where
    None). Layer 2 matches spending at the observed design by construction, so
    under spending 'observed' its spending-on-regressor rows are flagged `targeted2`. Under spec 'none' no row has a
    design: ask for the spending relations only (SPENDING). `usInput`, `carrier`: panelData's. """
    ps = {k: panelData(k, turnout, sample, reading, spec, betaRule, theta, spending, omegaRule, root, retAge,
                       inequality, usInput, carrier) for k in layers}
    common = ps[layers[0]].index
    for k in layers[1:]:
        common = common.intersection(ps[k].index)
    if keep is not None:
        common = common.intersection(pd.Index(keep))
    ps = {k: p.loc[common] for k, p in ps.items()}
    first = ps[layers[0]]
    obs = variant(turnout, sample, reading, spec, betaRule, spending)['spending'] == 'observed'
    rows = []
    for yRole, xRole in relations:
        rows.append({'y': yRole, 'x': xRole, 'data': ols(*_xy(first, yRole, xRole, 'Data'))}
                    | {'layer{}'.format(k): ols(*_xy(ps[k], yRole, xRole, 'Model')) if k in ps else None
                       for k in (1, 2)}
                    | {'targeted2': obs and yRole == 'spending' and xRole != 'index'})
    out = pd.DataFrame(rows)
    out.attrs.update({'n': len(common), 'dropped': first.attrs['dropped'], 'countries': list(common)})
    return out


def spendingSlopes(turnout = None, sample = None, **kw):
    """ Section 6's spending test: the two relations of SPENDING in the data and in layer 1 at config.US['oecd26']'s
    spendingSpec (no cost: no design is read) and the β rule 'common', at the observed designs; `kw` slopes()'s
    omegaRule, retAge, theta, root. """
    return slopes(turnout, sample, None, C.US['oecd26']['spendingSpec'], 'common', relations = SPENDING,
                  layers = (1,), spending = 'observed', **kw)


def layer1Omega(v, root = None, carrier = None):
    """ Layer 1's ω̄ at the variant `v` (a dict of variant()) under omegaRule 'mean': the `ωMean` its rows carry (the
    mean of layer 2's ω, runOECD26.meanOmega); None under 'anchor', where every country takes the anchor's. `root`,
    `carrier`: panelData's (costKw). """
    if v['omegaRule'] != 'mean':
        return None
    rows = D.oecd26Layer(1, v['turnout'], v['spec'], v['betaRule'], _thetaSource(v), root, 'mean', v['retAge'], None,
                         carrier)
    ω = rows.loc[rows['status'] == 'ok', 'ωMean'].astype(float).unique()
    if len(ω) != 1:
        raise ValueError('layer 1 at the mean ω ({}): {} values of ωMean'.format(v, len(ω)))
    return float(ω[0])


def thetaFit(layer, turnout = None, units = 'theta', omegaRule = None, root = None, retAge = None, theta = None,
             carrier = None):
    """ The predicted against the observed design of layer `layer` at the paper's variant but `turnout`, `omegaRule`,
    `retAge` and the measure `theta`, from `root` (oecd26Dir by default; `carrier` its anchor's, costKw): `panel` (panelData) and, over the countries
    other than the U.S. (the anchor), Pearson's `r` of model on data, `rmse` of the model's prediction, `rmseMean` that
    of predicting every country at the sample mean, `sdModel` and `sdData` the standard deviations (population, ddof
    0: sdData equals rmseMean) of predicted and observed, and `fit` (ols of model on data). `units` 'theta' (θ) or
    'index' (I = 1/(2 - θ)). `n` counts the U.S. in, `nStat` out. """
    p = panelData(layer, turnout, omegaRule = omegaRule, root = root, retAge = retAge, theta = theta, carrier = carrier)
    col = {'theta': ('θData', 'θModel'), 'index': ('indexData', 'indexModel')}[units]
    h = p.drop(index = 'USA')
    x, y = h[col[0]].to_numpy(float), h[col[1]].to_numpy(float)
    return {'panel': p, 'x': p[col[0]], 'y': p[col[1]], 'n': len(p), 'nStat': len(h),
            'r': float(np.corrcoef(x, y)[0, 1]), 'rmse': float(np.sqrt(np.mean((y - x)**2))),
            'rmseMean': float(np.sqrt(np.mean((x - x.mean())**2))), 'sdData': float(np.std(x)), 'sdModel': float(np.std(y)),
            'fit': ols(x, y), 'dropped': p.attrs['dropped']}


def spread(**kw):
    """ Why the model's slope of spending on the index is steeper than the data's (section 7), on layer 1's countries
    at the paper's variant (`kw`: panelData's): for `data` and `model`, the standard deviations (population, ddof 0,
    as thetaFit) of the index and of spending, `sdIndex` and `sdSpending`, their correlation `r`, and `r2Growth`, the
    share of the index's variance that population growth accounts for (the R2 of ols of the index on growth). """
    p = panelData(1, **kw)
    out = {}
    for side, tag in (('data', 'Data'), ('model', 'Model')):
        I, S = p['index' + tag].to_numpy(float), p['spending' + tag].to_numpy(float)
        out[side] = {'sdIndex': float(np.std(I)), 'sdSpending': float(np.std(S)),
                     'r': float(np.corrcoef(I, S)[0, 1]), 'r2Growth': ols(p['growth'], I)['r']**2}
    return out


def slopeDigits(y, x):
    """ The digits a slope of relation (y, x) prints with: the index on inequality config.INEQUALITY's `digits`, every
    other 2. """
    return C.INEQUALITYLABELS[C.INEQUALITY]['digits'] if (y, x) == ('index', 'inequality') else 2


def _cell(o, digits = 2, stack = False):
    r""" A slope and its standard error; `stack` puts the error on a line of its own (\newline), so that every cell of
    a narrow column breaks alike. """
    return '$' + C.num(o['slope'], digits) + '$' + (r'\newline' if stack else ' ') + '$(' + C.num(o['se'], digits) + ')$'


def _dashNote(body):
    r""" The note's key to `--`, only when a cell of `body` (rows of a table) prints one. """
    return ' --: not computed.' if any(re.search(r'(&\s*--\s*(&|\\\\))|\{--\}', b) for b in body) else ''


def _names(isos):
    """ ISO3 codes as figure 1.1's country names, joined by commas and 'and'. """
    names = [str(D.oecdFigure1Data().at[c, 'country']) for c in isos]
    return names[0] if len(names) == 1 else ', '.join(names[:-1]) + ' and ' + names[-1]


def sampleText(n, dropped):
    """ '20 countries (without Luxembourg, which has no CSES survey)' or '21 countries'. """
    return '{} countries'.format(n) + (' (without {}, which {} no CSES survey)'.format(
        _names(dropped), 'has' if len(dropped) == 1 else 'have') if dropped else '')


def _shared(a, v):
    """ The anchor's parameters every layer-1 country takes, as tex items: β (common rule), ω, the cost parameter. """
    return (([r'$\beta = ' + C.num(float(a['β']), 2) + '$'] if v['betaRule'] == 'common' else [])
            + [r'$\omega = ' + C.num(float(a['ω']), 2) + '$']
            + ([COSTSYMBOL[v['spec']][:-1] + ' = ' + C.num(float(a['p']), 2) + '$'] if v['spec'] in COSTSYMBOL else []))


def _and(items):
    return items[0] if len(items) == 1 else ', '.join(items[:-1]) + ' and ' + items[-1]


def _meanWords(v):
    """ The model's ω under omegaRule 'mean' in words, for a note that does not define own ω beside it. """
    return (r'the mean of the political weights that match each country\textquotesingle s spending'
            + (r' at its observed design' if v['spending'] == 'chosen' else ''))


def _layersText(a, v, matches = True, layers = (1, 2), ω = None):
    r""" The note's sentence on the model (layer 1: what every country shares with the anchor `a`) and, where `layers`
    holds 2, own ω (layer 2: what each country calibrates), at the variant `v`; `matches` adds that own ω then matches
    spending (under spending 'observed'); `ω` the model's ω under omegaRule 'mean' (layer1Omega(v) where None), the
    mean of the countries' own ω where own ω is defined beside it, else in words (_meanWords). """
    ω = layer1Omega(v) if ω is None else ω
    items = _shared(a, v)
    if ω is not None:
        items = [i for i in items if not i.startswith(r'$\omega')]
    both = 2 in layers
    mean = '' if ω is None else (r', at $\omega = ' + C.num(ω, 2) + r'$, '
                                 + (r'the mean of the countries\textquotesingle{} own $\omega$' if both else _meanWords(v)))
    tail = {'size': '', 'none': ', without the cost',
            'sizePiR': r', each at its cost ' + C.PIRCOST}[v['spec']] \
        + (r', with the $\beta$ that sets its 2020 return to the U.S.\ one' if v['betaRule'] == 'commonR' else '')
    calibrated = r'$\omega$' if v['betaRule'] == 'common' else r'$\beta$ and $\omega$'
    one = r'The model: every country' + (r' at the U.S.\ ' + _and(items) if items else '') + mean + tail
    two = (r'own $\omega$: each country\textquotesingle s ' + calibrated + r' calibrated to its spending'
           + (r' at its observed design' if v['spending'] == 'chosen' else '')
           + (r' and the U.S.\ 2020 return' if v['betaRule'] == 'commonR' else '')
           + (', which it then matches' if matches and v['spending'] == 'observed' else ''))
    if not both:
        return one + '.'
    return one + '; ' + two + '.'


def specText(a, v, matches = True, layers = (1, 2)):
    r""" The notes' account of the model, own ω (where `layers` holds 2) and the turnout. The paper's variant points to
    figure 7.2 (`OECD26_Layer1`), whose note and section 7.3 define the model, and to figure 7.3 (`OECD26_ThetaFit`),
    whose note defines own ω (notes/paper_styleGuide.md, "What a note holds"); a variant states them in full, since
    they differ from figures 7.2's and 7.3's. """
    if suffix(**v) == '':
        return (r'The model as in figure \ref{fig:OECD26:layer1}'
                + (r', own $\omega$ as in figure \ref{fig:OECD26:thetaFit}' if 2 in layers else '') + '.')
    return modelText(v) + ' ' + _layersText(a, v, matches, layers) + ' ' + turnoutNote(v['turnout'])


def modelText(v):
    r""" The note's sentence on what the model's index and spending are at the variant `v`. """
    design = ('chosen in 2020, which binds from 2050' if v['reading'] == 'chosen' else
              'in force in 2020, chosen in 1990')
    spend = (r'in 2020 on the path of the design it chooses' if v['spending'] == 'chosen' else
             r'in 2020 at the observed design')
    return (r'The model\textquotesingle s index is $1/(2-\theta)$ of the design ' + design
            + r', and its spending ' + spendingTex(v) + ' ' + spend + '.')


def _slopesCosts():
    r""" Table 5 (OECD26_Slopes at the paper's variant): the data beside the model under each cost of COSTS, each with
    its own anchor and mean ω (costKw). Rows keyed `<y>_<x>`. """
    ss = {k: slopes(layers = (1,), **costKw(k)) for k in COSTS}
    if len({tuple(s.attrs['countries']) for s in ss.values()}) != 1:
        raise ValueError('OECD26_Slopes: the two costs cover different countries')
    first = ss[COSTS[0]]
    body = []
    for i in range(len(first)):
        r = first.iloc[i]
        d = slopeDigits(r['y'], r['x'])
        body.append(' & '.join([YTEX[r['y']] + ' on ' + XTEX[r['x']], _cell(r['data'], d, stack = True)]
                               + [_cell(ss[k].iloc[i]['layer1'], d, stack = True) for k in COSTS])
                    + r' \\ % row: ' + rowKey(r['y'] + '_' + r['x']))
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, ' + '{} countries. '.format(first.attrs['n'])
            + r'The model under each cost as in figure \ref{fig:OECD26:layer1}.')
    return _xwrap('OECD26_Slopes', costsSrc(),
                  r'Figure \ref{fig:US:OECDdata}\textquotesingle s relations in the data and in the model',
                  'table:OECD26:slopes', 'l' + 'Y'*(1 + len(COSTS)),
                  [r'\textbf{Relation}', r'\textbf{Data}'] + [r'\textbf{Model, ' + COSTNAME[k] + '}' for k in COSTS],
                  '\n'.join(body), note, width = r'\textwidth')


def slopesLayers(v):
    """ The layers OECD26_Slopes shows at a variant `v` other than the paper's (table 5 is _slopesCosts): the model and
    own ω. """
    return (1,) if suffix(**v) == '' else (1, 2)


def oecd26Slopes(turnout = None, sample = None, reading = None, spec = None, betaRule = None, spending = None,
                 omegaRule = None, retAge = None, theta = None):
    r""" The slopes of figure 1.1's relations: the data on the sample and the model (layer 1: the U.S. cost and one
    political weight, the mean of the countries' own or the U.S. one), and at a variant other than the paper's also
    with each country's own political weight (own ω, layer 2; slopesLayers). """
    v = variant(turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta)
    if suffix(**v) == '':
        return _slopesCosts()
    layers = slopesLayers(v)
    s = slopes(**v, layers = layers)
    a = D.oecd26Anchor(v['turnout'], v['spec'], _thetaSource(v))
    body = []
    for _, r in s.iterrows():
        label = YTEX[r['y']] + ' on ' + XTEX[r['x']]
        d = slopeDigits(r['y'], r['x'])
        cells = ([label, _cell(r['data'], d, stack = True), _cell(r['layer1'], d, stack = True)]
                 + ([r'\textit{targeted}' if r['targeted2'] else _cell(r['layer2'], d, stack = True)] if 2 in layers
                    else []))
        body.append(' & '.join(cells) + r' \\ % row: ' + rowKey(r['y'] + '_' + r['x']))
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, ' + '{} countries. '.format(s.attrs['n'])
            + specText(a, v, layers = layers))
    sfx = suffix(**v)
    return _xwrap('OECD26_Slopes' + sfx,
                  D.oecd26DirRel(_thetaSource(v)) + ''.join('layer{}.csv, '.format(k) for k in layers)
                  + 'anchor.csv; ' + DATASRC,
                  r'Figure \ref{fig:US:OECDdata}\textquotesingle s relations in the data and in the model'
                  + variantCaption(**v),
                  'table:OECD26:slopes' + sfx, 'l' + 'Y'*(1 + len(layers)),
                  [r'\textbf{Relation}', r'\textbf{Data}'] + [r'\textbf{' + LAYERNAME[k] + '}' for k in layers],
                  '\n'.join(body), note, width = r'\textwidth')


def oecd26Omega(turnout = None, sample = None, reading = None, spec = None, betaRule = None, retAge = None,
                theta = None):
    r""" Own ω (layer 2) by country: the inputs that differ across countries (workers per retiree in 2020, the
    dispersion of relative incomes, spending, the index and the design it implies), the calibrated political weight and
    the design predicted at the U.S. cost, in the order of figure 1.1's footnote. The model's ω rule is no part of it. """
    v = variant(turnout, sample, reading, spec, betaRule, retAge = retAge, theta = theta)
    _check(v)
    p = panelData(2, **v)
    m = D.oecd26Layer(2, v['turnout'], v['spec'], v['betaRule'], _thetaSource(v), None, 'anchor',
                      v['retAge']).loc[p.index]
    a = D.oecd26Anchor(v['turnout'], v['spec'], _thetaSource(v))
    gap = max(float((m[k].astype(float) - p[k].astype(float)).abs().max()) for k in ('spendingData', 'indexData', 'θData'))
    if not gap <= 1e-12:
        raise D.StaleInput('{}: the data columns of layer 2 differ from the data side (data/oecd26_theta.csv) by {:.2e}: '
                           'rows calibrated to other data; rerun python/US/runOECD26.py'.format(
                               D.oecd26Path(2, _thetaSource(v)), gap))
    smp = D.oecd26Sample()
    outside = [c for c in p.index if not bool(smp.at[c, 'sample26'])]
    body = []
    for c in p.index:
        r = m.loc[c]
        tag = r'$^{\ast}$' if c == 'USA' else r'$^{\dagger}$' if c in outside else ''
        cells = [str(p.at[c, 'country']) + tag, C.num(float(r['ν2020']), 2),
                 C.num(float(r['Vtilde']), 3), C.num(float(p.at[c, 'spendingData']), 1),
                 C.num(float(p.at[c, 'indexData']), 3),
                 C.num(float(p.at[c, 'θData']), 3), C.num(float(r['ω']), 3), C.num(float(p.at[c, 'θModel']), 3),
                 C.num(float(p.at[c, 'indexModel']), 3)]
        body.append(' & '.join(cells) + r' \\ % row: ' + c.lower())
    design = ('chosen in 2020, in force from 2050' if v['reading'] == 'chosen' else
              'in force in 2020 on the freely simulated path, chosen in 1990')
    head, join = ((r'$\omega$ is calibrated to spending at the observed design at $\beta = ' + C.num(float(a['β']), 2)
                   + '$', ' and ') if v['betaRule'] == 'common' else
                  (r'$\beta$ and $\omega$ are calibrated to the U.S.\ 2020 return and to spending at the observed '
                   r'design', ' at '))
    cost = {'size': join + r'the U.S.\ $\lambda = ' + C.num(float(a['p']), 2) + '$',
            'sizePiR': join + r'the U.S.\ $' + C.PIRSYMBOL + ' = ' + C.num(float(a['p']), 2) + '$, ' + C.PIRCOST,
            'none': ', without the cost'}[v['spec']]
    anchor = (r' $^{\ast}$The U.S.\ is the anchor, where $\omega$ and ' + COSTSYMBOL[v['spec']] + ' are calibrated.'
              if v['spec'] in COSTSYMBOL else r' $^{\ast}$The U.S.\ is the anchor.')
    dagger = (r' $^{\dagger}$Outside the paper\textquotesingle s sample: informal employment above '
              + '{:g}'.format(float(smp['informalMax'].iloc[0])) + r'\% of employment.' if outside else '')
    if suffix(**v) == '':
        paperWording()
        note = (r'\textit{Note:} ' + '{} countries'.format(len(p)) + r'. $\nu_{2020}$, workers per retiree; '
                r'$\tilde V$, the dispersion of relative incomes at the U.S.\ income percentiles; spending in \% of '
                r'GDP; $\theta$ and $I = 1/(2-\theta)$ as in appendix \ref{app:oecd26:data}. $\omega$ is calibrated to '
                r'spending at the observed design (own $\omega$, figure \ref{fig:OECD26:thetaFit}); the last two columns are '
                r'the model\textquotesingle s design for 2020, chosen in 1990, and its index. $^{\ast}$The anchor, '
                r'where $\omega$ and $' + C.PIRSYMBOL + r'$ are calibrated.' + dagger)
    else:
        note = (r'\textit{Note:} ' + sampleText(len(p), p.attrs['dropped']) + r'. $\nu_{2020}$, workers per retiree; '
                r'$\tilde V$, the dispersion of relative incomes at the U.S.\ income percentiles; spending in \% of GDP; '
                + measureNote(v) + ' ' + head + cost + r'; the last two columns are the model\textquotesingle s design, '
                + design + r', and its index.' + anchor + dagger + ' ' + turnoutNote(v['turnout']))
    # two header rows, the data's columns and the model's under a spanner each, so that nine columns fit
    header = (r' & \multicolumn{5}{c}{\textbf{Data}} & \multicolumn{3}{c}{\textbf{Model, own $\omega$}} \\' + '\n'
              + r'\cmidrule(lr){2-6}\cmidrule(lr){7-9}' + '\n'
              + r'\textbf{Country} & $\nu_{2020}$ & $\tilde V$ & spending & $I$ & $\theta$ & $\omega$ & $\theta$ & $I$')
    sfx = suffix(**v)
    return _xwrap('OECD26_Omega' + sfx,
                  D.oecd26DirRel(_thetaSource(v)) + 'layer2.csv, anchor.csv; data/oecdFigure1.csv, oecd26_theta.csv',
                  r'The countries of the cross-country test, each calibrated to its pension spending'
                  + variantCaption(**v),
                  'table:OECD26:omega' + sfx, 'l' + 'Y'*8, header,
                  '\n'.join(body), note, width = r'\textwidth', size = r'\footnotesize\setlength{\tabcolsep}{3pt}')



# ---------------------------------------------------------------------------------------------------
# The variants of the supplemental appendix's page C4 (OECD26_Variants, ONLINE_ONLY): the paper's specification beside
# its population variants (ν uncorrected or at the exit age, config.US['oecd26']['populationVariants']) and its data
# variants (the measures of the design of config.US['oecd26']['thetaVariants'], at the paper's ν).
RETAGEORDER = ('none', 'nra', 'exit')
RELTEX = {('spending', 'growth'): 'Spending on population growth',
          ('spending', 'inequality'): 'Spending on ' + XTEX['inequality'],
          ('spending', 'index'): r'Spending on $I$', ('index', 'growth'): r'$I$ on population growth',
          ('index', 'inequality'): r'$I$ on ' + XTEX['inequality']}


def populationRows():
    """ The readings of ν page C4 sets side by side, in RETAGEORDER: the paper's and config's populationVariants. """
    have = {C.US['oecd26']['retAge']} | set(C.US['oecd26']['populationVariants'])
    return [r for r in RETAGEORDER if r in have]


def _variantSlopes(paper = False, **kw):
    """ slopes() of RELATIONS at the paper's specification but `kw`, or None while a variant's rows are missing or stale
    (MissingInput, StaleInput): `--`, never a stale value. The paper's own rows (`paper`) and anything else raise. """
    try:
        return slopes(**kw)
    except (D.MissingInput, D.StaleInput):
        if paper:
            raise
        return None


def oecd26Variants():
    r""" The slopes of figure 1.1's relations (OECD26_Slopes's five) under page C4's variants, one column per
    relation; the supplemental appendix's alone, so its label, table:OECD26:variants, resolves on the site only. Panel A, the population data: the data (the same at every ν) and the model and own ω (layers 1 and 2,
    rows `Model`, `Own $\omega$`) at each reading of ν (populationRows, the paper's starred), then the uncorrected ν
    without Luxembourg (turnout 'cses'), with its own data on its 20 countries. Panel B, the measure of the design: per
    measure (MEASURES, the paper's first) its own data, on its own countries, and the model and own ω, the data variants
    from results/oecd26_theta_<V>/ (stage (ii)'s oecd26Variants), at the paper's ν. A variant whose rows are missing or
    stale prints `--`; the paper's own (starred) raise; the starred own-ω rows are the paper's only own-ω slopes. Rows keyed
    `population_data`, `population_<retAge>_<layer1|layer2>`, `population_none_cses_<data|layer1|layer2>`,
    `measure_<θ>_<data|layer1|layer2>`. """
    pop = {r: _variantSlopes(r == PAPER['retAge'], retAge = r) for r in populationRows()}
    noLux = _variantSlopes(retAge = 'none', turnout = 'cses')         # Luxembourg is what moves the uncorrected ν
    meas = {m: _variantSlopes(m == PAPER['theta'], theta = m) for m in MEASURES}
    ncol = 2 + len(RELATIONS)

    def cells(s, side):
        if s is None:
            return ['--']*len(RELATIONS)
        t = s.set_index(['y', 'x'])
        return [r'\textit{targeted}' if side == 'layer2' and t.loc[k, 'targeted2'] else
                _cell(t.loc[k, side], slopeDigits(*k), stack = True) for k in RELATIONS]

    star = lambda on: r'$^{\ast}$' if on else ''
    cap = lambda s: s[0].upper() + s[1:]
    first = pop[PAPER['retAge']]
    body = [r'\multicolumn{%d}{l}{\textit{A. The population data, $\nu$}} \\' % ncol,
            keyed(' & '.join(['Data', ''] + cells(first, 'data')) + r' \\', 'population_data')]
    for r, s in pop.items():
        for side, lab in (('layer1', LAYERNAME[1]), ('layer2', LAYERNAME[2])):
            name = cap(C.RETAGENAMES[r]) + star(r == PAPER['retAge']) if side == 'layer1' else ''
            body.append(keyed(' & '.join([name, lab] + cells(s, side)) + r' \\', 'population_{}_{}'.format(r, side)))
    for side, lab in (('data', 'Data'), ('layer1', LAYERNAME[1]), ('layer2', LAYERNAME[2])):
        name = cap(C.RETAGENAMES['none']) + r', without Luxembourg' if side == 'data' else ''
        body.append(keyed(' & '.join([name, lab] + cells(noLux, side)) + r' \\', 'population_none_cses_{}'.format(side)))
    body.append(r'\multicolumn{%d}{l}{\textit{B. The measure of the design, $\theta$}} \\' % ncol)
    for m, s in meas.items():
        for side, lab in (('data', 'Data'), ('layer1', LAYERNAME[1]), ('layer2', LAYERNAME[2])):
            name = cap(C.THETANAMES[m]) + star(m == PAPER['theta']) if side == 'data' else ''
            body.append(keyed(' & '.join([name, lab] + cells(s, side)) + r' \\', 'measure_{}_{}'.format(m, side)))
    n = first.attrs['n']
    ns = ', '.join('{} {}'.format(C.THETANAMES[m], s.attrs['n']) for m, s in meas.items() if s is not None)
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, the U.S.\ included. Panel A: $\nu$ '
            r'corrected for each country\textquotesingle s gap to the U.S.\ normal retirement age, uncorrected, or '
            r'corrected for its gap in the effective age of labor-market exit'
            + r', on the {} countries of the paper\textquotesingle s measure'.format(n)
            + r'; the data do not depend on $\nu$. The last three rows leave out Luxembourg, which has no CSES survey'
            + (r', on {} countries'.format(noLux.attrs['n']) if noLux is not None else '')
            + r'. Panel B: each measure of $\theta$ is its own cross-section, the data '
            r'and the model calibrated to it, on its countries' + (' (' + ns + ')' if ns else '')
            + r', at the paper\textquotesingle s $\nu$. ' + specText(None, variant())
            + r' $^{\ast}$The paper\textquotesingle s specification; its data and model rows are those of table '
            r'\ref{table:OECD26:slopes}.' + _dashNote(body))
    header = r'\textbf{Variant} & & ' + ' & '.join(RELTEX[k] for k in RELATIONS)
    return _xwrap('OECD26_Variants',
                  'results/oecd26/, results/oecd26_theta_<V>/ layer1.csv, layer2.csv, anchor.csv; ' + DATASRC,
                  r'Figure \ref{fig:US:OECDdata}\textquotesingle s relations under other population data and other '
                  r'measures of the design',
                  'table:OECD26:variants', r'>{\raggedright\arraybackslash}p{2.5cm}l' + 'Y'*len(RELATIONS), header,
                  '\n'.join(body), note, width = r'\textwidth', size = r'\footnotesize\setlength{\tabcolsep}{3pt}')


# ---------------------------------------------------------------------------------------------------
# Page C4's further readings: section 6's test with partial slopes (OECD26_GivenDesign), the model with one input at its
# U.S. value (OECD26_USInputs), table 5 without the U.S. (OECD26_SlopesNoUS) and under the cost linear in the tax
# (OECD26_SlopesTau1). All ONLINE_ONLY.


def olsControl(x, z, y):
    """ Least squares of y on a constant, x and the control z: the coefficient on x and its classical standard error,
    and n. """
    x, z, y = (np.asarray(a, dtype = float) for a in (x, z, y))
    X = np.column_stack([np.ones_like(x), x, z])
    b, *_ = np.linalg.lstsq(X, y, rcond = None)
    e = y - X @ b
    V = np.sum(e**2)/(len(y) - X.shape[1])*np.linalg.inv(X.T @ X)
    return {'slope': float(b[1]), 'se': float(np.sqrt(V[1, 1])), 'n': len(y)}


GIVEN = (('spending', 'index'), ('spending', 'growth'), ('spending', 'inequality'))


def oecd26GivenDesign():
    r""" Section 6's test with the design given and no cost (config.US['oecd26']['spendingSpec'], the model's spending at
    each country's observed design; layer 1 at the paper's turnout, ω rule and ν): spending on the index I, on population
    growth and on income inequality, in the data and in the model, alone and controlling for population growth (the
    coefficient in the regression on the regressor and growth; none for growth itself). The model's index is the data's,
    the design being given. Rows keyed `<y>_<x>`. """
    v = variant(spec = C.US['oecd26']['spendingSpec'])
    p = panelData(1, spec = v['spec'], spending = 'observed')
    x = lambda role: p['indexData'] if role == 'index' else p[role]
    body = []
    for yR, xR in GIVEN:
        cells = [RELTEX[(yR, xR)] if (yR, xR) in RELTEX else 'Spending on ' + XTEX[xR]]
        for side in ('Data', 'Model'):
            y = p['spending' + side]
            cells.append(_cell(ols(x(xR), y), 2, stack = True))
            cells.append('--' if xR == 'growth' else _cell(olsControl(x(xR), p['growth'], y), 2, stack = True))
        body.append(' & '.join(cells) + r' \\ % row: ' + rowKey(yR + '_' + xR))
    a = D.oecd26Anchor(v['turnout'], v['spec'])
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, ' + sampleText(len(p), p.attrs['dropped'])
            + r'. Alone: the regression on the regressor; with growth: the coefficient on it in the regression on it and '
            r'population growth. The model is that of the test of section \ref{sec:oecd} (figure '
            r'\ref{fig:OECD26:spending}), without the cost: '
            r'every country at its observed design, so that its index is the data\textquotesingle s, and its spending '
            r'$\tau(1-\alpha)$ in 2020 there. ' + _layersText(a, v, layers = (1,)) + ' ' + turnoutNote(v['turnout'])
            + _dashNote(body))
    header = (r' & \multicolumn{2}{c}{\textbf{Data}} & \multicolumn{2}{c}{\textbf{Model}} \\' + '\n'
              + r'\cmidrule(lr){2-3}\cmidrule(lr){4-5}' + '\n'
              + r'\textbf{Relation} & alone & with growth & alone & with growth')
    return _xwrap('OECD26_GivenDesign', D.oecd26DirRel() + 'layer1.csv, anchor.csv; ' + DATASRC,
                  r'Pension spending in the data and in the model with the design given, without the cost',
                  'table:OA:OECD26_GivenDesign', r'>{\raggedright\arraybackslash}p{4.2cm}YYYY', header,
                  '\n'.join(body), note, width = r'\textwidth')


USINPUTNAMES = {'demography': 'demography', 'income': 'income distribution', 'turnout': 'turnout'}
USINPUTTEXT = {'demography': r'workers per retiree, after the correction for the retirement age',
               'income': r'the income groups, and with them $\tilde V$',
               'turnout': r'the propensities to vote by income group'}


def oecd26USInputs(k = 2):
    r""" The model at the paper's specification (layer 1) with one input of every country at its U.S. value
    (config.US['oecd26']['usInputs']; python/US/runOECD26.py --usInput), beside the data and the paper's model, under the
    cost with tau^k (costKw; OECD26_USInputs + config.costSuffix(k), its label too): the five slopes of OECD26_Slopes, and
    spending on the index controlling for population growth. A reading whose rows are missing or stale prints `--`; the
    paper's own raise. Rows keyed `<y>_<x>` and `spending_index_growth`. """
    v = variant()
    kw = costKw(k)
    sfx = C.costSuffix(k)
    symbol = C.CARRIERSYMBOL[kw.get('carrier', C.CARRIER)]
    cols = [('data', None), ('model', None)] + [('model', u) for u in C.US['oecd26']['usInputs']]
    panels = {}
    for _, u in cols:
        if u in panels:
            continue
        try:
            panels[u] = panelData(1, usInput = u, **kw)
        except (D.MissingInput, D.StaleInput):
            if u is None:
                raise
            panels[u] = None
    base = panels[None]
    common = base.index
    for u, pu in panels.items():
        if pu is not None:
            common = common.intersection(pu.index)

    def cell(side, u, yR, xR, control = False):
        pu = panels[u if side == 'model' else None]
        if pu is None:
            return '--'
        pu = pu.loc[common]
        tag = 'Data' if side == 'data' else 'Model'
        x, y = _xy(pu, yR, xR, tag)
        o = olsControl(x, pu['growth'], y) if control else ols(x, y)
        return _cell(o, slopeDigits(yR, xR), stack = True)

    body = []
    for yR, xR in RELATIONS:
        body.append(' & '.join([RELTEX[(yR, xR)]] + [cell(sd, u, yR, xR) for sd, u in cols])
                    + r' \\ % row: ' + rowKey(yR + '_' + xR))
    body.append(' & '.join([r'Spending on $I$, with population growth'] + [cell(sd, u, 'spending', 'index', True)
                                                                          for sd, u in cols])
                + r' \\ % row: spending_index_growth')
    D.oecd26Anchor(v['turnout'], v['spec'], **kw)
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, ' + '{} countries'.format(len(common))
            + r'; the last row is the coefficient on $I$ in the regression on $I$ and population growth. The model as in '
            r'figure \ref{fig:OECD26:layer1}, under ' + COSTTEX[k] + r'; in the last three columns every country takes '
            r'one input at its U.S.\ value, '
            + _and([r'\textit{' + USINPUTNAMES[u] + '}: ' + USINPUTTEXT[u] for u in C.US['oecd26']['usInputs']])
            + r', every other input, its spending datum and observed design among them, its own, and $\beta$, $\omega$ '
            r'and $' + symbol + r'$ those of the model under that cost. The regressors are the '
            r'data\textquotesingle s, except the index, which is each column\textquotesingle s own.' + _dashNote(body))
    header = (r' & & & \multicolumn{3}{c}{\textbf{Model, one input at its U.S.\ value}} \\' + '\n'
              + r'\cmidrule(lr){4-6}' + '\n'
              + r'\textbf{Relation} & \textbf{Data} & \textbf{Model} & '
              + ' & '.join(USINPUTNAMES[u] for u in C.US['oecd26']['usInputs']))
    src = (D.oecd26DirRel() if k == 2 else 'results/' + C.US['tau1']['oecd26'] + '/') + 'layer1.csv, anchor.csv; ' + DATASRC
    return _xwrap('OECD26_USInputs' + sfx, src,
                  r'Figure \ref{fig:US:OECDdata}\textquotesingle s relations in the model with one input at its U.S.\ value, '
                  + COSTTEX[k], 'table:OA:OECD26_USInputs' + sfx, r'>{\raggedright\arraybackslash}p{3.6cm}YYYYY', header,
                  '\n'.join(body), note, width = r'\textwidth', size = r'\footnotesize')


def oecd26SlopesNoUS():
    r""" Table 5 (OECD26_Slopes at the paper's variant: the data and the model) on the sample and without the U.S., the
    anchor where the cost parameter, β and the U.S. ω are calibrated. Rows keyed `<y>_<x>`. """
    v = variant()
    every = slopes(layers = (1,))
    noUS = slopes(layers = (1,), keep = [c for c in every.attrs['countries'] if c != 'USA'])
    body = []
    for (_, r), (_, q) in zip(every.iterrows(), noUS.iterrows()):
        d = slopeDigits(r['y'], r['x'])
        body.append(' & '.join([RELTEX[(r['y'], r['x'])], _cell(r['data'], d, stack = True),
                                _cell(q['data'], d, stack = True), _cell(r['layer1'], d, stack = True),
                                _cell(q['layer1'], d, stack = True)]) + r' \\ % row: ' + rowKey(r['y'] + '_' + r['x']))
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, on the ' + str(every.attrs['n'])
            + r' countries of table \ref{table:OECD26:slopes} and on the ' + str(noUS.attrs['n'])
            + r' without the U.S., the anchor of the model. ' + specText(None, v, layers = (1,)))
    header = (r' & \multicolumn{2}{c}{\textbf{Data}} & \multicolumn{2}{c}{\textbf{Model}} \\' + '\n'
              + r'\cmidrule(lr){2-3}\cmidrule(lr){4-5}' + '\n'
              + r'\textbf{Relation} & all & without the U.S. & all & without the U.S.')
    return _xwrap('OECD26_SlopesNoUS', D.oecd26DirRel() + 'layer1.csv, anchor.csv; ' + DATASRC,
                  r'Figure \ref{fig:US:OECDdata}\textquotesingle s relations without the U.S.',
                  'table:OA:OECD26_SlopesNoUS', r'>{\raggedright\arraybackslash}p{4.2cm}YYYY', header, '\n'.join(body),
                  note, width = r'\textwidth')


def oecd26SlopesTau1():
    r""" Table 5's slopes under the paper's cost (quadratic in the tax, k = 2) and under the cost linear in the tax
    (config.US['tau1'], its own anchor and own ω in datasets.tau1Dir('oecd26')), the model at the paper's specification
    otherwise, beside the data. Rows keyed `<y>_<x>`. """
    t = C.US['tau1']
    k2 = slopes(layers = (1,))
    k1 = slopes(layers = (1,), root = D.tau1Dir('oecd26'), carrier = t['carrier'])
    if list(k1.attrs['countries']) != list(k2.attrs['countries']):
        raise ValueError('OECD26_SlopesTau1: the two costs cover different countries')
    a1 = D.oecd26Anchor(PAPER['turnout'], PAPER['spec'], root = D.tau1Dir('oecd26'), carrier = t['carrier'])
    body = []
    for (_, r), (_, q) in zip(k2.iterrows(), k1.iterrows()):
        d = slopeDigits(r['y'], r['x'])
        body.append(' & '.join([RELTEX[(r['y'], r['x'])], _cell(r['data'], d, stack = True),
                                _cell(r['layer1'], d, stack = True), _cell(q['layer1'], d, stack = True)])
                    + r' \\ % row: ' + rowKey(r['y'] + '_' + r['x']))
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, ' + '{} countries. '.format(k2.attrs['n'])
            + r'$k = 2$: the model as in figure \ref{fig:OECD26:layer1}, its cost the tax rate squared. $k = 1$: the cost '
            r'linear in the tax rate, $f = \exp(-\frac{1}{2}\lambda_{n,t}\tau\tilde V(1-\theta)^2)$ with $\lambda_{n,t} = '
            r'\varepsilon\,\pi_{n,t}$, $\varepsilon = ' + C.num(float(a1['p']), 2) + r'$ calibrated on the U.S.\ anchor and '
            r'every country at the mean of its own $\omega$ under that cost; the rest as at $k = 2$.')
    header = [r'\textbf{Relation}', r'\textbf{Data}', r'\textbf{Model}, $\bm{k = 2}$', r'\textbf{Model}, $\bm{k = 1}$']
    return _xwrap('OECD26_SlopesTau1', D.oecd26DirRel() + ', results/' + t['oecd26'] + '/ layer1.csv, anchor.csv; '
                  + DATASRC, r'Figure \ref{fig:US:OECDdata}\textquotesingle s relations under the cost linear in the tax',
                  'table:OA:OECD26_SlopesTau1', r'>{\raggedright\arraybackslash}p{4.2cm}YYY', header, '\n'.join(body),
                  note, width = r'\textwidth')


# ---------------------------------------------------------------------------------------------------
# Page D1's section on inequality (OECD26_Layer1Inequality, ONLINE_ONLY): figure 7.2's spending and index on each
# measure of INEQUALITYROWS, the data and the model at the paper's variant.


def inequalitySlopes(p = None):
    """ Least squares (ols) of the data's and the model's spending and index on each measure of INEQUALITYROWS, on
    panelData(1)'s countries at the paper's variant (`p`, read when None): one row per (measure, y, side), side 'data'
    | 'model', with slope, se, intercept, r and n. """
    p = panelData(1) if p is None else p
    rows = []
    for m in INEQUALITYROWS:
        for y in ('spending', 'index'):
            for side in ('Data', 'Model'):
                rows.append({'measure': m, 'y': y, 'side': side.lower()} | ols(p[m].astype(float), p[y + side]))
    return pd.DataFrame(rows)


def oecd26Layer1InequalitySlopes():
    r""" OECD26_Layer1InequalitySlopes (page D1): the slopes of figure 7.2's panels on inequality on the Gini and on the
    model's two measures, the data and the model at the paper's variant, with Pearson's r (results/paper/
    oecd26Layer1InequalitySlopes.csv, written by the OECD26_Layer1Inequality build). Rows keyed `<y>_<measure>`. The
    csv is checked against the current results (inequalitySlopes): one of other rows or data is StaleInput. """
    df = D.oecd26InequalitySlopes()
    cur = inequalitySlopes()
    key = ['measure', 'y', 'side']
    have = df.set_index(key).reindex(cur.set_index(key).index)
    num = ['slope', 'se', 'r', 'n']
    if have[num].isna().any().any() or len(df) != len(cur) or not np.allclose(
            cur[num].values.astype(float), have[num].values.astype(float), rtol = 0, atol = 1e-9):
        raise D.StaleInput('results/paper/oecd26Layer1InequalitySlopes.csv does not hold the slopes of the current '
                           'results (results/oecd26/layer1.csv, ' + DATASRC + '): rebuild OECD26_Layer1Inequality')
    t = df.set_index(key)
    n = int(t['n'].iloc[0])
    body = []
    for m in INEQUALITYROWS:
        for y in ('spending', 'index'):
            d = 3 if y == 'index' else 2
            cells = [YTEX[y] + ' on ' + C.INEQUALITYLABELS[m]['text']
                     + (r'$^{\ast}$' if m == C.INEQUALITY else '')]
            for side in ('data', 'model'):
                o = t.loc[(m, y, side)]
                cells += [_cell(o, d), '$' + C.num(float(o['r']), 2) + '$']
            body.append(' & '.join(cells) + r' \\ % row: ' + rowKey(y + '_' + m.lower()))
    a = D.oecd26Anchor()
    note = (r'\textit{Note:} Least-squares slopes, standard errors in parentheses, and Pearson\textquotesingle s $r$, '
            + '{} countries. '.format(n) + specText(a, variant(), layers = (1,))
            + r' The Gini is that of disposable income; ' + C.INEQUALITYLABELS['yHighOverLow']['text']
            + r', $y_3/y_1$, and the dispersion $\tilde V = \sum_i \gamma_i (y_i - 1)^2/y_i$ are of the income groups '
            r'the model is fed. $^{\ast}$The measure of figure \ref{fig:OECD26:layer1} and table '
            r'\ref{table:OECD26:slopes}.')
    header = (r' & \multicolumn{2}{c}{\textbf{Data}} & \multicolumn{2}{c}{\textbf{' + LAYERNAME[1] + r'}} \\' + '\n'
              + r'\cmidrule(lr){2-3}\cmidrule(lr){4-5}' + '\n'
              + r'\textbf{Relation} & slope & $r$ & slope & $r$')
    return _xwrap('OECD26_Layer1InequalitySlopes', 'results/paper/oecd26Layer1InequalitySlopes.csv',
                  r'Figure \ref{fig:OECD26:layer1}\textquotesingle s relations with income inequality under three '
                  r'measures of it', 'table:OA:OECD26_Layer1InequalitySlopes',
                  r'>{\raggedright\arraybackslash}p{5.2cm}YYYY', header, '\n'.join(body), note, width = r'\textwidth')

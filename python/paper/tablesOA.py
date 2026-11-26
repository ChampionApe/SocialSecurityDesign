r""" Table builders for the supplemental appendix (writing/OnlineAppendix; archive/notes/brief_onlineAppendix_2026-10-06.md,
section 2.3). Each returns the complete tex body of one file in results/paper/Tables/, which is not copied
into writing/Paper. Layout as tablesUS: tables.BANNER and tablesUS._xwrap, label `table:OA:<Name>`, a note
in the house form.

Every body row ends with ` % row: <key>`, the last thing on its line; block titles, sub-headers and rules
carry none (brief section 2.2). `_body` refuses a duplicate key.

Units are the paper tables': savings over GDP throughout (the ESC csvs carry s/(w h) and go through
datasets.escSavingsOverY), tax and savings rates in percent to two decimals, their changes in p.p.
(config.pp), the design theta to three decimals. The cell formatters the figures reuse for their marks' text
(ARGRHOCOLS, OECDRHOCOLS, escPathCells) live here.
"""
import os, re
import numpy as np

import config as C
import datasets as D
import oecdFigure1 as F1
from tables import LQ, RQ
from tablesUS import _xwrap, COUNTRYNAME


# ---------------------------------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------------------------------
def _row(cells, key, tail = ''):
    r""" One body row: the cells, the row end and any spacing `tail`, then the key comment last. """
    return ' & '.join(cells) + r' \\' + tail + ' % row: ' + key


def _title(text, ncol):
    r""" A block title spanning the table: no key. """
    return r'\multicolumn{' + str(ncol) + r'}{l}{\textit{' + text + r'}} \\'


def _body(lines):
    """ The joined body, refusing a duplicate row key; a block gap on the last row is dropped. """
    keys = [l.split(' % row: ')[1] for l in lines if ' % row: ' in l]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        raise ValueError('duplicate row keys: {}'.format(dup))
    lines = list(lines)
    lines[-1] = lines[-1].replace(r' \\[.5em] % row: ', r' \\ % row: ')
    return '\n'.join(lines)


def _signed(x, digits = 2):
    """ A change as a signed math-mode number, unsigned when it rounds to zero (as config.pp). """
    v = round(float(x), digits)
    return C.num(0., digits) if v == 0 else '${:+.{d}f}$'.format(v, d = digits)


def _dash(x, fmt):
    """ `fmt(x)`, or '--' for a value the csv does not hold (NaN). """
    return '--' if x is None or x != x else fmt(x)


def _rhoKey(ρ):
    return C.num(float(ρ), 1)


def _tex(s):
    r""" Plain text from a csv into a tex cell: the special characters escaped, straight double quotes as
    tex quotes, apostrophes as \textquotesingle (the house form). """
    s = str(s)
    out = []
    for ch in s:
        out.append({'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_',
                    '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
                    '|': r'$|$', '<': r'$<$', '>': r'$>$', "'": r'\textquotesingle{}'}.get(ch, ch))
    s = ''.join(out)
    return re.sub(r'"([^"]*)"', LQ + r'\1' + RQ, s)


def _keyOf(s):
    """ An ASCII row key from a name: lowercase, runs of anything outside [a-z0-9] as one '_'. """
    return re.sub(r'[^a-z0-9]+', '_', str(s).lower()).strip('_')


def _fileKey(name):
    """ A file name's row key, its extension dropped: 'US_shocks.csv' -> 'us_shocks'. """
    return _keyOf(os.path.splitext(name)[0])


def _label(name):
    return 'table:OA:' + name


# ---------------------------------------------------------------------------------------------------
# Figure 1's data (data.qmd)
# ---------------------------------------------------------------------------------------------------
ROLES = F1.ROLES                                                # oecdFigure1's variable roles, in order
ROLEFMT = {'spending': lambda v: C.num(v, 1), 'growth': lambda v: C.num(v, 3),
           'index': lambda v: C.num(v, 3),
           'inequality': lambda v: C.INEQUALITYLABELS[C.INEQUALITY]['fmt'].format(v)}


def oecdCountries():
    r""" Figure 1's cross-section, one row per country it plots (oecdFigure1.plotted: config.US['theta']'s measure's
    countries): the four variables oecdFigure1.PLOT draws, each followed by the year it is read at where that varies
    by country (data/oecdFigure1.csv; the index from data/oecd26_theta.csv). Rows keyed on lowercase ISO3; the
    economies section 6 calibrates (oecdFigure1.HIGHLIGHT) set in bold. """
    df = F1.plotted(D.oecdFigure1Data())
    cols = []                                         # (header, csv column, formatter)
    for role in ROLES:
        c = F1.PLOT[role]
        cols.append((r'\textbf{' + _texLabel(c) + '}', c, ROLEFMT[role]))
        if c + '_year' in df.columns:
            cols.append((r'\textbf{Year}', c + '_year', lambda v: str(int(round(float(v))))))
    lines = []
    for iso, r in df.iterrows():
        name = _tex(r['country'])
        if iso in F1.HIGHLIGHT:
            name = r'\textbf{' + name + '}'
        lines.append(_row([name] + [_dash(r[c], f) for _, c, f in cols], _keyOf(iso)))
    note = (r'\textit{Note:} The ' + str(len(df)) + r' OECD members plotted in figure '
            r'\ref{fig:US:OECDdata}; the Bismarckian index is $1/(2-\theta)$ of the measure of the design in its note. '
            r'Each variable is read in ' + str(F1.YEAR)
            + r' or the nearest year within ' + str(F1.WINDOW) + r' years of it (the year columns); population '
            r'growth is the population of ' + str(F1.YEAR) + r' over that of ' + str(F1.BASE) + r'. Sources in '
            r'table \ref{' + _label('OECD_Sources') + r'}. In bold the economies the paper calibrates.')
    return _xwrap('OECD_Countries', 'data/oecdFigure1.csv, data/oecd26_theta.csv' + F1.distSource(),
                  r'The cross-section of figure \ref{fig:US:OECDdata}, by country',
                  _label('OECD_Countries'), 'l' + 'Y'*len(cols), ['\\textbf{Country}'] + [h for h, _, _ in cols],
                  _body(lines), note, width = r'\textwidth')


EPS = re.compile(r'\bEPS\b|eps_')


def _texLabel(column):
    """ oecdFigure1.LABEL[column] as a tex cell: the text escaped (_tex), its math ($...$) kept. """
    return ''.join(part if k % 2 else _tex(part) for k, part in enumerate(re.split(r'(\$[^$]*\$)', F1.LABEL[column])))


def _withoutEPS(note):
    """ A sources note without its sentences on the figure's previous version (the EPS and its eps_* readings),
    which data/README.md documents. """
    return ' '.join(s for s in re.split(r'(?<=\.)\s+(?=[A-Z])', str(note)) if not EPS.search(s))


def _plottedColumns():
    """ The csv columns figure 1 draws. """
    return set(F1.PLOT.values())


def _thetaSourceRow():
    r""" OECD_Sources' row of the design θ and its index 1/(2 - θ): config.US['theta']'s measure, variant C of
    data/oecd26_theta.csv (app:oecd26:data), its source, the earnings levels and the countries read off the data
    (oecdFigure1.thetaSource). The two sentences the note types are checked: the sample rule against
    oecd26Theta.PRIVMAX, the one exception against the measure's notes; another measure is a ValueError. """
    from oecd26Theta import PRIVMAX
    if C.US['theta'] != 'C' or PRIVMAX != 0.25:
        raise ValueError('OECD_Sources describes variant C at PRIVMAX 0.25; config.US[\'theta\'] is {!r}, PRIVMAX {}'
                         .format(C.US['theta'], PRIVMAX))
    s = F1.thetaSource()
    if set(s['inForce']) != {'GBR'}:
        raise ValueError('OECD_Sources names the UK as the one θ not read off the benefit formula; the measure\'s notes '
                         'name {}'.format(sorted(s['inForce'])))
    names = D.oecdFigure1Data()['country']
    low = ', '.join('{} at {:g}'.format(names[c], v) for c, v in s['low'].items())
    source = (_tex(s['source']) + '; retrieved ' + _tex(s['retrieved']) + r'. Top-ups and tested benefits: OECD SOCX '
              r'2020 and national sources (data/oecd26\_topups.csv; python/paper/oecd26Theta.py, S\_TOPUP)')
    series = (r'$\theta$ is ' + C.THETADEF + r': the gross public pension of men with a full career under the rules '
              r'of 2020' + (' (lowest earnings: ' + low + ')' if low else '') + ', public spending as pensionSpending. '
              r'The index is $1/(2-\theta)$. Years: 2020 (' + '{n} of {n})'.format(n = s['n']))
    note = (r'The paper\textquotesingle s measure (appendix \ref{app:oecd26:data}), variant C of '
            r'data/oecd26\_theta.csv, on the ' + str(s['n']) + r' countries where funded tiers pay less than a '
            r'quarter of mandatory benefits in 2020 and the public pension is not tested against them. The '
            r'UK\textquotesingle s $\theta$ is not read off the benefit formula: it is the State Second '
            r'Pension\textquotesingle s share of spending on the State Pensions and Pension Credit in 2019--20.')
    return [r'$\theta$, the design, and its index', source, series, note]


def oecdSources():
    r""" The sources of figure 1's data, one row per entry of data/oecdFigure1_sources.csv, one for the design θ and
    its index (_thetaSourceRow) and one per column oecdFigure1.plotted derives (oecdFigure1.derivedSources:
    config.INEQUALITY's measure of the model's income groups): the columns it sources, provider and dataset with
    vintage and retrieval date, the series, the years and the note. The entries whose columns figure 1 draws come
    first. URLs are in the csvs, not printed. Rows keyed on the entry's first column, the design's `theta`. """
    import pandas as pd
    src = D.oecdFigure1Sources()
    src = src[[not all(c.strip().startswith('eps_') for c in s.split(',')) for s in src['column']]]
    derived = F1.derivedSources()
    src = pd.concat([src, derived], ignore_index = True)
    plotted = _plottedColumns()
    first = lambda s: s.split(',')[0].strip()
    drawn = [any(c.strip() in plotted for c in s.split(',')) for s in src['column']]
    blocks = [('Plotted in the figure of the paper', src[drawn]),
              ('Further concepts in data/oecdFigure1.csv', src[[not d for d in drawn]])]
    lines = []
    for title, sub in blocks:
        if sub.empty:
            continue
        lines.append(_title(_tex(title), 4))
        if title.startswith('Plotted'):
            lines.append(_row(_thetaSourceRow(), 'theta'))
        for _, r in sub.iterrows():
            source = '; '.join(_tex(x) for x in (r['provider'] + ': ' + r['dataset'], r['vintage'],
                                                  ('retrieved ' + r['retrieved']) if r['retrieved'] else '') if x)
            series = _tex(r['series']) + ('. Years: ' + _tex(r['years']) if r['years'] else '')
            lines.append(_row([_tex(r['column']), source, series, _tex(_withoutEPS(r['note'])) or '--'],
                              _keyOf(first(r['column']))))
    extra = [first(c) for c in derived['column']]
    note = (r'\textit{Note:} The columns are those of data/oecdFigure1.csv, but for $\theta$ and its index, from '
            r'data/oecd26\_theta.csv (data/README.md)'
            + (', and ' + ', '.join(_tex(c) for c in extra) + r', derived from the income groups of '
               r'data/oecd29\_distribution.csv' if extra else '')
            + r'; data/oecdFigure1\_sources.csv, data/oecd26\_pensionTiers\_harmonised.csv'
            + (r', data/oecd29\_sources.csv' if extra else '') + r' and the files the $\theta$ row names carry the URL '
            r'of every source. bbIndex, the ratio of the replacement rates at average and at half of average '
            r'earnings, is not the paper\textquotesingle s index.')
    rag = r'>{\raggedright\arraybackslash}'
    text = _xwrap('OECD_Sources', 'data/oecdFigure1_sources.csv, data/oecd26_theta.csv, '
                  'data/oecd26_pensionTiers_harmonised.csv' + (', data/oecd29_sources.csv' if extra else ''),
                  r'The sources of figure \ref{fig:US:OECDdata}',
                  _label('OECD_Sources'), rag + 'p{3.8cm}' + (rag + 'X')*3,
                  [r'\textbf{Columns}', r'\textbf{Source}', r'\textbf{Series and years}', r'\textbf{Note}'],
                  _body(lines), note, width = r'\textwidth')
    return text.replace('\\centering\n', '\\centering\n\\footnotesize\n', 1)


def oecdCorrelations():
    r""" The pairwise correlations of figure 1's variables in the concepts and on the countries it plots
    (results/paper/oecdCorrelations.csv, written by the OECDdata build, variant 'figure' = oecdFigure1.PLOT on
    oecdFigure1.plotted): n, Pearson's r and Spearman's rank correlation with their p-values. The csv's other concept
    sets, on the 29 (the named concepts, the figure's previous version), are not printed (data/README.md). Rows
    keyed `named_<x role>_<y role>`. The csv is checked against the figure's current data (n, r and rho per pair on
    oecdFigure1.plotted): a csv of another sample or of an earlier data/oecdFigure1.csv is StaleInput. """
    df = D.oecdCorrelations()
    sub = df[df['variant'] == 'figure']
    if sub.empty:
        raise D.MissingInput('results/paper/oecdCorrelations.csv (variant figure)')
    cur = F1.figureCorrelations(F1.plotted()).set_index(['x', 'y'])
    have = sub.set_index(['x', 'y'])
    if (sorted(cur.index) != sorted(have.index) or (cur['n'] != have['n'].reindex(cur.index)).any()
            or not np.allclose(cur[['pearson', 'spearman']].values,
                               have[['pearson', 'spearman']].reindex(cur.index).values.astype(float), rtol = 0, atol = 1e-9)):
        raise D.StaleInput('results/paper/oecdCorrelations.csv does not hold the correlations of the current figure '
                           '(data/oecdFigure1.csv, data/oecd26_theta.csv): rebuild OECDdata')
    role = {v: k for k, v in F1.PLOT.items()}
    label = _texLabel
    p = lambda v: r'$<0.001$' if float(v) < 0.0005 else C.num(float(v), 3)
    lines = [_row([label(r['x']) + ' and ' + label(r['y']), str(int(r['n'])),
                   '$' + C.num(float(r['pearson']), 2) + '$', p(r['pearson_p']),
                   '$' + C.num(float(r['spearman']), 2) + '$', p(r['spearman_p'])],
                  'named_' + role[r['x']] + '_' + role[r['y']]) for _, r in sub.iterrows()]
    note = (r'\textit{Note:} The concepts plotted in figure \ref{fig:US:OECDdata}. Pairwise complete '
            r'observations; the $p$-values are those of the two-sided tests of zero correlation. Variables as in '
            r'table \ref{' + _label('OECD_Countries') + r'} and their sources in table \ref{'
            + _label('OECD_Sources') + '}.')
    return _xwrap('OECD_Correlations', 'results/paper/oecdCorrelations.csv', r'Correlations in the cross-section of figure \ref{fig:US:OECDdata}',
                  _label('OECD_Correlations'), r'>{\raggedright\arraybackslash}p{6.4cm}YYYYY',
                  [r'\textbf{Variables}', '$n$', r'\textbf{Pearson} $r$', '$p$', r'\textbf{Spearman}', '$p$'],
                  _body(lines), note, width = r'\textwidth')


# ---------------------------------------------------------------------------------------------------
# Page D1's section on inequality: the Gini against the model's two measures (OECD_InequalityMeasures, its correlations
# OECD_InequalityCorrelations; figure 7's panels on the three measures are tablesOECD26.inequalitySlopes').
INEQUALITYMEASURES = ('yHighOverLow', 'Vtilde')          # the model's measures, set against the Gini


def inequalityMeasures():
    """ One row per country figure 1 plots (oecdFigure1.plotted), on ISO3: `country`, the Gini (`gini`, `gini_year`)
    and the model's measures of INEQUALITYMEASURES with the groups' `year` (datasets.incomeGroups). """
    fig = F1.plotted()
    g = D.incomeGroups(fig.index)
    out = fig[['country', 'gini', 'gini_year']].copy()
    for m in INEQUALITYMEASURES:
        out[m] = g[m]
    out['year'] = g['year']
    return out


def inequalityCorrelations(m = None):
    """ The Gini against each measure of INEQUALITYMEASURES on inequalityMeasures' countries `m`, all of them (`sample`
    'all') and without the U.S. ('noUSA'): n, Pearson's r and Spearman's rank correlation with their two-sided p-values,
    and the least-squares line of the measure on the Gini (tablesOECD26.ols: slope, its classical standard error,
    intercept). One row per (measure, sample). """
    from scipy import stats
    from tablesOECD26 import ols
    import pandas as pd
    m = inequalityMeasures() if m is None else m
    rows = []
    for col in INEQUALITYMEASURES:
        for sample, keep in (('all', m.index), ('noUSA', m.index[m.index != 'USA'])):
            x, y = m.loc[keep, 'gini'].astype(float), m.loc[keep, col].astype(float)
            pr, sr, o = stats.pearsonr(x, y), stats.spearmanr(x, y), ols(x, y)
            rows.append({'measure': col, 'sample': sample, 'n': len(keep), 'pearson': float(pr[0]),
                         'pearson_p': float(pr[1]), 'spearman': float(sr[0]), 'spearman_p': float(sr[1]),
                         'slope': o['slope'], 'se': o['se'], 'intercept': o['intercept']})
    return pd.DataFrame(rows)


SAMPLETEX = {'all': 'All countries', 'noUSA': r'Without the U.S.'}


def oecdInequalityCorrelations():
    r""" OECD_InequalityCorrelations (page D1): the Gini against the model's measures of income inequality on figure 1's
    countries, with and without the U.S. (results/paper/oecdInequalityCorrelations.csv, written by the
    OECD_InequalityMeasures build): n, Pearson's r and Spearman's rank correlation with their p-values, the
    least-squares slope of the measure on the Gini. Rows keyed `<measure>_<all|nousa>`. The csv is checked against
    the current data (inequalityCorrelations): one of other data or countries is StaleInput. """
    df = D.oecdInequalityCorrelations()
    cur = inequalityCorrelations()
    key = ['measure', 'sample']
    have = df.set_index(key).reindex(cur.set_index(key).index)
    num = ['n', 'pearson', 'spearman', 'slope', 'se']
    if have[num].isna().any().any() or not np.allclose(cur[num].values.astype(float), have[num].values.astype(float),
                                                         rtol = 0, atol = 1e-9) or len(df) != len(cur):
        raise D.StaleInput('results/paper/oecdInequalityCorrelations.csv does not hold the correlations of the current '
                           'data (data/oecdFigure1.csv, oecd29_distribution.csv, oecd26_theta.csv): rebuild '
                           'OECD_InequalityMeasures')
    p = lambda v: r'$<0.001$' if float(v) < 0.0005 else C.num(float(v), 3)
    lines = []
    for col in INEQUALITYMEASURES:
        d = 3 if col == 'Vtilde' else 2
        lines.append(_title(_texLabel(col), 7))
        for _, r in df[df['measure'] == col].iterrows():
            lines.append(_row([SAMPLETEX[r['sample']], str(int(r['n'])),
                               '$' + C.num(float(r['pearson']), 2) + '$', p(r['pearson_p']),
                               '$' + C.num(float(r['spearman']), 2) + '$', p(r['spearman_p']),
                               '$' + C.num(float(r['slope']), d) + r'$ $(' + C.num(float(r['se']), d) + ')$'],
                              _keyOf(col) + '_' + _keyOf(r['sample'])))
    note = (r'\textit{Note:} The countries of figure \ref{fig:US:OECDdata}. Correlations of the Gini of disposable income '
            r'with each measure of the income groups the model is fed, at the U.S.\ shares of the population, and the '
            r'least-squares slope of the measure on the Gini, its standard error in parentheses; the $p$-values are those '
            r'of the two-sided tests of zero correlation. The measures are those of the figure of the Gini against the model\textquotesingle s measures, above.')
    return _xwrap('OECD_InequalityCorrelations', 'results/paper/oecdInequalityCorrelations.csv',
                  r'The Gini and the model\textquotesingle s measures of income inequality',
                  _label('OECD_InequalityCorrelations'), r'>{\raggedright\arraybackslash}p{3.2cm}YYYYYY',
                  [r'\textbf{Countries}', '$n$', r'\textbf{Pearson} $r$', '$p$', r'\textbf{Spearman}', '$p$',
                   r'\textbf{Slope on the Gini}'],
                  _body(lines), note, width = r'\textwidth')


# ---------------------------------------------------------------------------------------------------
# Argentina (argentina.qmd)
# ---------------------------------------------------------------------------------------------------
# (csv column, header, formatter) of the Argentine calibration across rho; ARG_RhoGrid's marks print with
# the formatters.
ARGRHOCOLS = [('β',  r'$\beta$',    lambda v: C.num(v, 3)),
              ('ω',  r'$\omega$',   lambda v: C.num(v, 3)),
              ('η0', r'$\eta_0$',   lambda v: C.num(v, 3)),
              ('X0', r'$X_0$',      lambda v: C.num(v, 3)),
              ('X',  r'$X$',        lambda v: C.num(v, 2)),
              ('sr', r'\textbf{Savings rate}', lambda v: C.pct(v)),
              ('ι',  r'$\iota$',    lambda v: C.num(v, 3))]


def argRhoGrid(commonX = None):
    """ The Argentine calibration sweep of one variant over config.ARG['ρGrid'], raising if a rho is missing. """
    g = D.rhoGrid(commonX)
    return _onGrid(g, C.ARG['ρGrid'], C.argSweepCsv(C.ARG['commonX'] if commonX is None else commonX))


def _onGrid(g, grid, path):
    """ The rows of sweep `g` at the rhos of `grid`, raising MissingInput naming `path` if one is absent. """
    miss = [ρ for ρ in grid if not np.isclose(g['ρ'], ρ).any()]
    if miss:
        raise D.MissingInput('{} (ρ = {})'.format(path, miss))
    return g[g['ρ'].apply(lambda r: any(np.isclose(r, v) for v in grid))].reset_index(drop = True)


# ---------------------------------------------------------------------------------------------------
# The rich economies with the design given (oecd.qmd)
# ---------------------------------------------------------------------------------------------------
OECDHOSTS = (('US', 'usa'), ('UK', 'gbr'), ('FR', 'fra'))      # (sweep key, mark-id part), the paper's order
OECDRHOCOLS = [('β',  r'$\beta$',  lambda v: C.num(v, 3)),
               ('ω',  r'$\omega$', lambda v: C.num(v, 3)),
               ('X',  r'$X$',      lambda v: C.num(v, 2)),
               ('R',  r'$R$',      lambda v: C.num(v, 3)),
               ('sr', r'\textbf{Savings rate}', lambda v: C.pct(v))]


def oecdSweeps(commonX = None):
    """ {country: its calibration sweep over config.US['ρGrid']} for the U.S., the UK and France, raising if a
    rho is missing from any. """
    cx = C.US['commonX'] if commonX is None else commonX
    return {c: _onGrid(D.usSweep(c, cx), C.US['ρGrid'], C.usSweepCsv(c, cx)) for c, _ in OECDHOSTS}


# ---------------------------------------------------------------------------------------------------
# The endogenous design (esc.qmd)
# ---------------------------------------------------------------------------------------------------
def escPathCells(r):
    """ One baseline-path row as the tables print it: the design chosen, the tax chosen and pinned, their
    difference. US_ESC_Path's marks print the first three. """
    return [C.num(r['θ'], 3), C.pct(r['τ']), C.pct(r['τExo']), C.pp(r['τ'] - r['τExo'])]


COSTNAME = {'scale': r'The cost on the design'}         # config.US['esc']['spec'] is named by its k (_costName)


def _costName(k, capital = True):
    """ 'The linear cost' / 'the squared cost' (config.COSTWORDS) of the cost with tau^k. """
    s = C.COSTWORDS[k]
    return s[0].upper() + s[1:] if capital else s


def _escRel(k, names):
    """ 'results/<escDir(k)>/<names>' for a banner. """
    return 'results/' + os.path.basename(D.escDir(k)) + '/' + names


def escCostCalibration(k = 2):
    r""" C2.2: the costs of redistribution calibrated under log preferences (rho = 1, the only rho at which
    config.US['esc']['comparisonSpec'] is calibrated; escCalibration.csv of datasets.escDir(k) through
    datasets.escCalibrationLOG): the cost f = exp(-lambda tau^k Vtilde (1-theta)^2 / 2) with tau^k (k = 1
    eq:esc:fk at k = 2; config.COSTEQ) and the cost on the design f = phi + (1-phi) theta^p,
    which does not depend on k, each with its parameter and the shares of revenue lost in 2020 at the observed design
    and at theta = 0 (datasets.escWedge), at the 2020 tax. At k = config.US['esc']['tauPower'] the first row prints as
    the rho = 1 row of table:US_ESC:calibration (lambda at two decimals, the shares at one; test_onlineAppendix compares
    them). Built as ESC_CostCalibration + config.costSuffix(k). Rows keyed `spec@<spec>`. """
    ρ = C.US['ρAnchor']
    spec, alt = C.US['esc']['spec'], C.US['esc']['comparisonSpec']
    name = 'ESC_CostCalibration' + C.costSuffix(k)
    lines, θs, phi = [], [], None
    for s, symbol, digits in ((spec, r'\lambda', 2), (alt, 'p', 3)):
        r = D.escCalibrationLOG(s, k = k)
        w = D.escWedge(r, s)
        θs.append(float(r['θStar']))
        if s == alt:
            phi = float(r['phi'])
        lines.append(_row([_costName(k) if s == spec else COSTNAME[s],
                           '$' + symbol + ' = ' + C.num(float(r['p']), digits) + '$',
                           C.pct(1 - w['fStar'], 1), C.pct(1 - w['f0'], 1)], 'spec@' + s))
    if max(θs) - min(θs) > 1e-6:
        raise ValueError('{}: theta* differs between the two costs: {}'.format(name, θs))
    year0 = str(C.usCalendar()['year0'])
    header = [r'\textbf{Cost}', r'\textbf{Parameter}', r'$1 - f(\theta^{\ast})$', r'$1 - f(0)$']
    paper = k == C.US['esc']['tauPower']
    note = (r'\textit{Note:} Log preferences ($\rho = ' + C.num(ρ, 0) + r'$). ' + _costName(k) + r' is that of '
            r'equation \eqref{' + C.COSTEQ[k] + '}' + (r', as in table \ref{table:US_ESC:calibration}' if paper else '')
            + r'; the cost on the design is $f(\theta) = \phi + (1-\phi)\theta^{p}$ with $\phi = ' + C.num(phi, 1)
            + r'$ imposed. Each parameter is calibrated so that the U.S.\ electorate re-elects its observed design, '
            r'$\theta^{\ast} = ' + C.num(θs[0], 3) + r'$, in ' + year0 + r'. $1 - f(\theta^{\ast})$ and $1 - f(0)$ are '
            r'the shares of revenue lost in ' + year0 + r' at the observed design and at a flat benefit, both at the '
            + year0 + r' tax rate.' + C.variantNote(C.US['commonX']))
    return _xwrap(name, _escRel(k, 'escCalibration.csv'),
                  r'The costs of redistribution calibrated for the U.S.: ' + _costName(k, False) + ' and the cost on '
                  r'the design', _label(name), 'lYYY', header, _body(lines), note, width = r'.8\textwidth')


def escTiming():
    r""" C2.1: the design chosen without the cost under the two timings of app:esc:corner.

      * With the tax (prop:esc:corner), per rho of C5_criticalTurnout.csv (datasets.criticalTurnout, its rows at
        tax 'calibrated': the calibrated equilibrium of every dated period, the paper's numbers): theta = 0 where
        eq:esc:seqFOC is negative at every date and every theta_t of the grid on [0, 1], '--' otherwise; the
        observed ratio of the high- to the low-income group's weight mu_H/mu_L; the critical ratio in t0 and its
        lowest over the dates, under weights (mu_i^CPS)^m. The critical ratio under weights y_i^k must agree to the
        printed digits (the high/low pair sets the sign, app:esc:corner), else ValueError.
      * One period ahead (eq:esc:leadFOC), log preferences only: theta_{t0}, the design in force at t0 on the free
        path without the cost, the 'none' row of escCalibration.csv (`choice`, ModelESC.leadedDesignAtT0).

    Rows keyed `withTax@<rho>`, `ahead@<rho>`. """
    ct = D.criticalTurnout()
    ct = ct[ct['tax'] == 'calibrated']
    if ct.empty:
        raise D.MissingInput('results/numerical/C5_criticalTurnout.csv (no row at the calibrated tax)')
    cal = C.usCalendar()
    year0, yearBefore = cal['year0'], cal['dates'][cal['t0'] - 1]
    isTrue = lambda s: s.astype(str).isin(['True', 'true', '1', '1.0'])
    ratio = lambda v: _dash(v, lambda x: C.num(x, 2))
    lines, every, dates = [], True, sorted(int(d) for d in ct['date'].unique())
    for k, (ρ, sub) in enumerate(ct.groupby('ρ')):
        neg = bool(isTrue(sub['negativeEverywhere']).all())
        every &= neg
        at0 = sub[sub['date'] == year0]
        if len(at0) != 1:
            raise D.MissingInput('results/numerical/C5_criticalTurnout.csv (rho {}, {})'.format(ρ, year0))
        at0 = at0.iloc[0]
        if np.ptp(sub['turnoutRatio'].astype(float)) > 1e-12:
            raise ValueError('ESC_Timing: the observed turnout ratio moves across dates at rho {}'.format(ρ))
        crit0, low = float(at0['ratioCrit_μm']), float(sub['ratioCrit_μm'].min())
        for a, b in ((crit0, float(at0['ratioCrit_yk'])), (low, float(sub['ratioCrit_yk'].min()))):
            if ratio(a) != ratio(b):
                raise ValueError('ESC_Timing: the critical ratio depends on the weights\' form at rho {}: {} vs {}'
                                 .format(ρ, a, b))
        lines.append(_row(['With the tax' if k == 0 else '', C.num(ρ, 1), C.num(0., 3) if neg else '--',
                           ratio(float(at0['turnoutRatio'])), ratio(crit0), ratio(low)], 'withTax@' + _rhoKey(ρ)))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')
    ρLOG = C.US['ρAnchor']
    none = D.escCalibrationLOG('none')
    lines.append(_row(['One period ahead', C.num(ρLOG, 1), C.num(float(none['choice']), 3), '--', '--', '--'],
                      'ahead@' + _rhoKey(ρLOG)))
    header = (r' & & & \multicolumn{3}{c}{\textbf{High- to low-income weight}} \\' + '\n'
              r'\cmidrule(lr){4-6}' + '\n'
              r'\textbf{Timing} & $\rho$ & $\theta$ \textbf{chosen} & \textbf{Observed} & \textbf{Critical}, '
              + str(year0) + r' & \textbf{Critical}, lowest')
    note = (r'\textit{Note:} No cost of redistribution; $\rho = ' + C.num(ρLOG, 0) + r'$ is log preferences. With the '
            r'tax: the sign of the derivative \eqref{eq:esc:seqFOC} on a grid of $\theta_t$ on $[0,1]$, at the '
            r'calibrated equilibrium of every date from ' + str(dates[0]) + ' to ' + str(dates[-1])
            + (r'; it is negative at every date and design, so the choice is the corner $\theta_t = 0$. '
               if every else '. ')
            + r'High- to low-income weight: the ratio of the two groups\textquotesingle{} weights in the political '
              r'objective; observed, the ratio of their propensities to vote in the calibration; critical, the ratio '
              r'at which the derivative first turns non-negative at some design as the weights rise more steeply '
              r'with income, in ' + str(year0) + r' and at its lowest over those dates. It is the same, to the digits '
              r'shown, whether the weights are the observed propensities raised to a common power or a power of '
              r'relative income. One period ahead: $\theta_{' + str(year0) + r'}$, the design the electorate of '
            + str(yearBefore) + r' chooses by \eqref{eq:esc:leadFOC}, on the path where the choice binds from the '
              r'first period.')
    return _xwrap('ESC_Timing', 'results/numerical/C5_criticalTurnout.csv, results/esc/escCalibration.csv',
                  'The design chosen without the cost, under two timings', _label('ESC_Timing'), 'lYYYYY',
                  header, _body(lines), note)


# Panel A's routes of eq:esc:leadFOC (datasets.designDecomp), in the table's order; the last two offset the first.
DECOMPROUTES = (('design', r'The pull: the design, the share lost and the tax held'),
                ('cost', r'The saving on the loss'),
                ('tax', r'The response of the anticipated tax'),
                ('total', r'Sum: the design condition'))
# Panel B's steps, (key, label, (to, from)) as differences of the 'choice' experiments, the U.S. ('US') first.
DECOMPSTEPS = (('step1', 'Everywhere but in the cost', ('X2', 'US')),
               ('step2', 'The dispersion in the cost', ('X1', 'X2')),
               ('step3', 'The tax rate in the cost', ('FR', 'X1')),
               ('total', 'Total', ('FR', 'US')))
# The relative gap of the cost to its target below which c5Numbers.py's design (log) and designCRRA stop.
DECOMPGAP = {'LOG': 1e-10, 'CRRA': 2e-3}
LAMBDAEPS = (0.1, 0.4)          # the range of participation and taxable-income elasticities, c5Numbers.EPS's ends


def _numRel(k, names):
    """ 'results/<numDir(k)>/<names>' for a banner or a message. """
    return 'results/' + os.path.basename(D.numDir(k)) + '/' + names


def _decompChoices(ρ, exp, k = 2):
    r""" Panel B's designs at `rho` under the cost with tau^k: {'US', 'FR', 'X1', 'X2'[, 'X3']: theta_2020} and the
    largest gap of the X rows. US, FR: escExperiments.csv's free baseline and frIncome rows (`exp`, datasets.escDir(k));
    the X rows from c5Numbers.py's design (rho = 1, whose US and FR rows must be escExperiments', else StaleInput) or
    designCRRA (datasets.numDir(k)). """
    spec = C.US['esc']['spec']
    θ = {e: float(D.escRow(exp, ρ, spec, sc, False)['θ_t0']) for e, sc in (('US', 'baseline'), ('FR', 'frIncome'))}
    LOG = np.isclose(ρ, C.US['ρAnchor'])
    if LOG:
        dd = D.designDecomp(k)
        ch = dd[(dd['kind'] == 'choice') & dd['experiment'].isin(['US', 'FR', 'X1', 'X2', 'X3'])].set_index('experiment')
        for e in ('US', 'FR'):
            if e not in ch.index or abs(float(ch.loc[e, 'θ2020']) - θ[e]) > 1e-6:
                raise D.StaleInput('{}: theta_2020 of {} is not escExperiments.csv\'s ({})'.format(
                    _numRel(k, 'C5_designDecomp.csv'), e, θ[e]))
        xs = ('X1', 'X2', 'X3')
    else:
        ch = D.designDecompCRRA(ρ, k).set_index('experiment')
        xs = ('X1', 'X2')
    miss = [x for x in xs if x not in ch.index]
    if miss:
        raise D.MissingInput('{} (experiments {})'.format(
            _numRel(k, 'C5_designDecomp{}.csv'.format('' if LOG else 'CRRA_rho{:g}'.format(ρ))), miss))
    gap = max(float(ch.loc[x, 'gap']) for x in xs)
    if not gap < DECOMPGAP['LOG' if LOG else 'CRRA']:
        raise ValueError('ESC_Decompositions: the cost misses its target by {:.2e} at rho {}'.format(gap, ρ))
    θ.update({x: float(ch.loc[x, 'θ2020']) for x in xs})
    return θ, gap


def escDecompositions(k = 2):
    r""" C2.3: the numbers section 7 and appendix C.2 state of the design choice under the cost with tau^k, one column
    per rho of config.US['ρGrid'] (c5Numbers.py's csvs in datasets.numDir(k), the ESC csvs in datasets.escDir(k)).
    Built as ESC_Decompositions + config.costSuffix(k).

      A. Log preferences only: dW_t/dtheta_{t+1} of eq:esc:leadFOC at the U.S. choice by route (datasets.designDecomp,
         case US, piece W; DECOMPROUTES), the routes' sum (the total, zero at the choice) and the shares of the offset,
         the saving on the loss ('cost') and the tax's response ('tax') each over their sum. The routes must sum to
         the total, and the total vanish, to 1e-6 of the pull, else ValueError.
      B. theta_2020 as France's income distribution is imposed on the U.S. in three steps (DECOMPSTEPS: the cost first
         at the U.S. size and dispersion, X2, then at the U.S. size, X1, then France's own) and in total; under log
         also X3 - US, France's dispersion in the cost alone (_decompChoices).
      C. App:esc:cost's magnitude (datasets.lambdaElasticity; its lambda must be escCalibration's, else StaleInput):
         pi_2020 and lambda; at k = 2 (the carrier 'convex', eps(T) = kappa T) the wedge in place eq:esc:wedge, the
         local elasticity under lambda = 2 kappa pi (1 - pi) (columns _A), the lambda that the elasticities LAMBDAEPS
         imply and lambda as a multiple of each, the share of earnings outside the base, and the local elasticity and
         that share under the fiscal externality's 3 kappa pi (1 - pi) (columns _C); at k = 1 (the carrier 'pi', a
         constant elasticity) the elasticity lambda/pi that lambda = eps pi implies, the lambda of each elasticity of
         LAMBDAEPS and the multiples (columns _ΠR).

    Rows keyed `<panel>@<item>`. """
    if not C.US['commonX']:
        raise D.StaleInput('results/numerical/C5_*.csv are at common X (python/US/c5Numbers.py), not the headline')
    name = 'ESC_Decompositions' + C.costSuffix(k)
    carrier = 'convex' if k == 2 else C.US['tau1']['carrier']
    ρs = [float(ρ) for ρ in C.US['ρGrid']]
    ρLOG = C.US['ρAnchor']
    ncol = 1 + len(ρs)
    cal = C.usCalendar()
    year0, yearBefore = cal['year0'], cal['dates'][cal['t0'] - 1]
    lines = []
    atLOG = lambda s: [s if np.isclose(ρ, ρLOG) else '--' for ρ in ρs]

    # A. the design condition at the U.S. choice, log preferences
    dd = D.designDecomp(k)
    foc = dd[(dd['kind'] == 'foc') & (dd['case'] == 'US') & (dd['piece'] == 'W')].set_index('route')
    if set(r for r, _ in DECOMPROUTES) - set(foc.index):
        raise D.MissingInput(_numRel(k, 'C5_designDecomp.csv') + ' (the routes of the U.S. design condition)')
    v = {r: float(foc.loc[r, 'value']) for r, _ in DECOMPROUTES}
    pull = abs(v['design'])
    if abs(v['design'] + v['cost'] + v['tax'] - v['total']) > 1e-6*pull or abs(v['total']) > 1e-6*pull:
        raise ValueError('{}: the routes of the design condition do not sum to zero at the choice: {}'.format(name, v))
    θ1 = float(foc['θ1'].iloc[0])
    lines.append(_title(r'A. The design condition at the U.S.\ choice, log preferences', ncol))
    for r, label in DECOMPROUTES:
        lines.append(_row([label] + atLOG(_signed(v[r], 4)), 'A@' + r))
    offset = v['cost'] + v['tax']
    for r, label in (('cost', 'the saving on the loss'), ('tax', 'the response of the tax')):
        lines.append(_row(['Share of the offset: ' + label] + atLOG(C.pct(v[r]/offset, 0)), 'A@share_' + r))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')

    # B. France's income distribution in three steps
    exp = D.escExperiments(k)
    θ, gaps = {}, {}
    for ρ in ρs:
        θ[ρ], gaps[ρ] = _decompChoices(ρ, exp, k)
    lines.append(_title(r'B. France\textquotesingle s income distribution in three steps, the change in $\theta_{'
                        + str(year0) + '}$', ncol))
    for key, label, (a, b) in DECOMPSTEPS:
        lines.append(_row([label] + [_signed(θ[ρ][a] - θ[ρ][b], 3) for ρ in ρs], 'B@' + key))
    lines.append(_row([r'The dispersion in the cost alone']
                      + [_signed(θ[ρ]['X3'] - θ[ρ]['US'], 3) if 'X3' in θ[ρ] else '--' for ρ in ρs], 'B@alone'))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')

    # C. the size of lambda
    le = D.lambdaElasticity(k).set_index('ρ')
    e0, e1 = LAMBDAEPS
    if k == 2:
        sfx = 'A'
        rowsC = [('pi', r'$\pi = \Pi/R$', 'πBaseline2020', lambda x: C.num(x, 3)),
                 ('tbar', r'The wedge in place, $\bar{t}$', 'tbar2020', lambda x: C.num(x, 3)),
                 ('lambda', r'$\lambda$', 'λ', lambda x: C.num(x, 2)),
                 ('eps', r'The local elasticity, $\varepsilon(\bar{t})$', 'εLocal_A', lambda x: C.num(x, 2))]
    else:
        sfx = 'ΠR'
        rowsC = [('pi', r'$\pi = \Pi/R$', 'πBaseline2020', lambda x: C.num(x, 3)),
                 ('lambda', r'$\lambda$', 'λ', lambda x: C.num(x, 2)),
                 ('eps', r'The elasticity of the base, $\varepsilon = \lambda/\pi$', 'εImplied_ΠR',
                  lambda x: C.num(x, 2))]
    lamDigits = 2 if k == 2 else 3
    rowsC += [('lambda_eps{:g}'.format(e0), r'$\lambda$ at $\varepsilon = ' + C.num(e0, 1) + '$',
               'λ_ε{:g}_{}'.format(e0, sfx), lambda x: C.num(x, lamDigits)),
              ('lambda_eps{:g}'.format(e1), r'$\lambda$ at $\varepsilon = ' + C.num(e1, 1) + '$',
               'λ_ε{:g}_{}'.format(e1, sfx), lambda x: C.num(x, lamDigits)),
              ('multiple_eps{:g}'.format(e1), r'$\lambda$ as a multiple of that at $\varepsilon = ' + C.num(e1, 1) + '$',
               'multiple_ε{:g}_{}'.format(e1, sfx), lambda x: C.num(x, 1)),
              ('multiple_eps{:g}'.format(e0), r'$\lambda$ as a multiple of that at $\varepsilon = ' + C.num(e0, 1) + '$',
               'multiple_ε{:g}_{}'.format(e0, sfx), lambda x: C.num(x, 1))]
    if k == 2:
        rowsC += [('outside', 'Earnings outside the base', 'baseOutside_A', lambda x: C.pct(x, 1)),
                  ('eps_fiscal', r'Fiscal externality: $\varepsilon(\bar{t})$', 'εLocal_C', lambda x: C.num(x, 2)),
                  ('outside_fiscal', r'Fiscal externality: earnings outside the base', 'baseOutside_C',
                   lambda x: C.pct(x, 1))]
    need = [col for _, _, col, _ in rowsC]
    miss = [c for c in need if c not in le.columns]
    if miss or any(not np.isclose(le.index, ρ).any() for ρ in ρs):
        raise D.MissingInput('{} (columns {}, rho {})'.format(_numRel(k, 'C5_lambdaElasticity.csv'), miss, ρs))
    if 'carrier' in le.columns and set(le['carrier'].dropna()) != {carrier}:
        raise D.StaleInput('{}: carrier {}, the mapping of k = {} is {!r}'.format(
            _numRel(k, 'C5_lambdaElasticity.csv'), set(le['carrier']), k, carrier))
    escCal = D.escCalibration(k = k)
    row = lambda ρ: le.loc[le.index[np.isclose(le.index, ρ)][0]]
    for ρ in ρs:
        λcal = float(escCal[(ρ, C.US['esc']['spec'])]['p'])
        if abs(float(row(ρ)['λ']) - λcal) > 1e-6*λcal:
            raise D.StaleInput('{}: lambda {} at rho {}, calibrated {}'.format(
                _numRel(k, 'C5_lambdaElasticity.csv'), float(row(ρ)['λ']), ρ, λcal))
    cells = lambda col, fmt: [fmt(float(row(ρ)[col])) for ρ in ρs]
    lines.append(_title(r'C. The size of $\lambda$, ' + str(year0), ncol))
    for key, label, col, fmt in rowsC:
        lines.append(_row([label] + cells(col, fmt), 'C@' + key))

    gapCRRA = max(g for ρ, g in gaps.items() if not np.isclose(ρ, ρLOG))
    header = [''] + [r'$\rho = {:g}$'.format(ρ) for ρ in ρs]
    T = r'\mathcal{T}_{' + str(year0) + '}'
    paper = k == C.US['esc']['tauPower']
    lamRef = r'$\lambda$ of table \ref{table:US_ESC:calibration}' if paper else r'$\lambda$ the calibrated cost parameter'
    if k == 2:
        panelC = (r'$\pi = \Pi/R$ on the path without the cost, $\bar{t} = \tau(1-\pi)$ the wedge in place of eq.\ '
                  r'\eqref{eq:esc:wedge} at the calibrated tax rate, ' + lamRef + r', and $\varepsilon(\bar{t}) = '
                  r'\kappa\bar{t} = \lambda\tau/(2\pi)$ the local elasticity of the base under $\lambda = 2\kappa\pi(1-\pi)$. '
                  r'$\lambda$ at $\varepsilon$ is the cost parameter that a local elasticity $\varepsilon$ implies, '
                  r'$2\pi\varepsilon/\tau$, at the ends of the range of participation and taxable-income elasticities '
                  r'\parencite{ChettyGMW11a1,SaezSG12a1}, and the multiples are the calibrated $\lambda$ over it. Earnings '
                  r'outside the base: the share $\frac{1}{2}\kappa\bar{t}^{2}$ that the wedge in place keeps outside. '
                  r'Fiscal externality: the budget charged with the participation tax of every unit of base that leaves, '
                  r'$\lambda = 3\kappa\pi(1-\pi)$.')
    else:
        panelC = (r'$\pi = \Pi/R$ on the path without the cost, ' + lamRef + r', and $\varepsilon = \lambda/\pi$ the '
                  r'elasticity of the base that $\lambda = \varepsilon\pi$ implies under a constant elasticity. $\lambda$ '
                  r'at $\varepsilon$ is the cost parameter $\varepsilon\pi$ that an elasticity $\varepsilon$ implies, at '
                  r'the ends of the range of participation and taxable-income elasticities '
                  r'\parencite{ChettyGMW11a1,SaezSG12a1}, and the multiples are the calibrated $\lambda$ over it.')
    total = (r'the response to France\textquotesingle s income distribution in figure \ref{fig:US_ESC:overview}' if paper
             else r'the response to France\textquotesingle s income distribution under this cost '
                  r'(\oa{esc-costs-design})')
    note = (r'\textit{Note:} ' + _costName(k) + r', equation \eqref{' + C.COSTEQ[k] + r'}; $\rho = ' + C.num(ρLOG, 0)
            + r'$ is log preferences. Panel A: the derivative of the '
            r'political objective of the electorate of ' + str(yearBefore) + r' in the design it chooses for '
            + str(year0) + r', eq.\ \eqref{eq:esc:leadFOC}, at the U.S.\ choice $\theta_{' + str(year0) + '} = '
            + C.num(θ1, 3) + r'$, by route: the design with the share of revenue lost and the anticipated tax held, '
            r'which pulls the design towards flat benefits; the dependence of the share lost on the design (eq.\ '
            r'\eqref{eq:esc:ftheta}), the saving on the loss; and the response of the anticipated tax $' + T
            + r'(\theta_{' + str(year0) + r'})$. The routes sum to zero at the choice; each share is one of the last '
            r'two routes over their sum, which offsets the pull. Panel B: $\theta_{' + str(year0) + r'}$ on the path '
            r'where the choice binds from the first period, as France\textquotesingle s income distribution is imposed '
            r'on the U.S., first everywhere but in the cost, whose share of revenue lost at each design is held at its '
            r'U.S.\ value, then in the dispersion of relative incomes $\tilde{V}$ in the cost, and last in the tax rate '
            r'in the cost; the first step is the combined effect of France\textquotesingle s distribution outside the '
            r'cost, and the total is ' + total + r'. In the first two steps $\lambda$ is multiplied at each date by a '
            r'factor, taken as given by voters, that equates the share lost at the tax rate of the path solved with the '
            r'share at the U.S.\ tax rate of that date, so that at a candidate tax off the path the share lost still '
            r'changes ' + ('with the square of the tax rate' if k == 2 else 'in proportion to the tax rate')
            + r'. Under CRRA preferences the share held is met to a relative gap of at most ' + C.pct(gapCRRA, 2) + '. '
            r'The dispersion in the cost alone is the U.S.\ with France\textquotesingle s $\tilde{V}$ in the cost, '
            r'computed under log preferences only. Panel C: the magnitudes of appendix \ref{app:esc:cost}, ' + panelC)
    src = (_numRel(k, 'C5_designDecomp.csv, C5_designDecompCRRA_rho*.csv, C5_lambdaElasticity.csv') + ', '
           + _escRel(k, 'escExperiments.csv, escCalibration{,CRRA}.csv'))
    return _xwrap(name, src, r'The cost and the design: the decompositions of section \ref{sec:esc}, '
                  + _costName(k, False), _label(name), '>{\\raggedright\\arraybackslash}p{6.2cm}' + 'Y'*len(ρs), header,
                  _body(lines), note)


def _stationaryBlocks():
    """ [(arm title, [(row label, key block, frame, gap column, formatter)])] of NUM_Stationary, each frame
    cut to the dated periods before the terminal one. """
    def cut(df):
        last = df.groupby('ρ')['t'].transform('max')
        return df[df['t'] < last]
    us, esc, arg = (cut(D.stationaryApprox(w)) for w in ('US', 'US_ESC', 'ARG'))
    pp = lambda v: _signed(v, 2)
    th = lambda v: _signed(v, 3)
    return [('United States', [('Tax rate, design given (p.p.)', 'us', us, 'gap_pp', pp),
                               ('Tax rate, design chosen (p.p.)', 'esctax', esc, 'gapτ_pp', pp),
                               (r'Design chosen ($\theta$)', 'escdesign', esc, 'gapθ', th)]),
            ('Argentina', [('Tax rate (p.p.)', 'arg', arg, 'gap_pp', pp)])]


def numStationary():
    r""" Stationary against date-specific policy functions (results/numerical/*_stationaryApprox*.csv): the
    path walked by a time-invariant policy function solved at each date's own nu_t, minus the paper's
    date-specific solution, by date and rho. The tax in p.p. (the U.S. with the design given, under CRRA; the
    U.S. with the design chosen, log and CRRA; Argentina, log and CRRA) and the chosen design in units of
    theta. The terminal period is not printed. Rows keyed `<block>@<rho>`. The csvs' other columns (the gap
    at the exact path's states, the steady-state and long-run comparisons) are not printed: the technical
    documentation does not define them. """
    blocks = _stationaryBlocks()
    ncol = 2 + max(len(sorted(df['date'].unique())) for _, rows in blocks for _, _, df, _, _ in rows)
    lines = []
    for arm, rows in blocks:
        dates = sorted(int(d) for d in rows[0][2]['date'].unique())
        for _, _, df, _, _ in rows:
            if sorted(int(d) for d in df['date'].unique()) != dates:
                raise ValueError('NUM_Stationary: the {} blocks have different dates'.format(arm))
        lines.append(r'\multicolumn{2}{l}{\textit{' + arm + r'}} & ' + ' & '.join(str(d) for d in dates)
                     + ' &'*(ncol - 2 - len(dates)) + r' \\')
        lines.append(r'\cmidrule(lr){3-' + str(2 + len(dates)) + '}')
        for label, block, df, col, fmt in rows:
            for k, (ρ, sub) in enumerate(df.groupby('ρ')):
                sub = sub.set_index('date')
                cells = [label if k == 0 else '', C.num(ρ, 1)] + [fmt(sub.loc[d, col]) for d in dates]
                cells += [''] * (ncol - len(cells))
                lines.append(_row(cells, block + '@' + _rhoKey(ρ)))
            lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')
    note = (r'\textit{Note:} Each cell is the policy on the path walked by time-invariant policy functions, '
            r'each solved at its own date' + "'" + r's $\nu_t$ and applied from the same initial state, minus '
            r'the date-specific policy of the paper' + "'" + r's solution: the tax rate in percentage points, '
            r'the chosen design in units of $\theta$. Under log preferences ($\rho = 1$) with the design given '
            r'the two coincide and the U.S.\ block has no row. The design chosen is at the paper' + "'" + r's cost '
            r'specification. The terminal period of the finite horizon is not shown; the last period shown '
            r'carries the date-specific solution\textquotesingle s approach to it (appendix \ref{app:T}) rather '
            r'than the demographic transition.')
    return _xwrap('NUM_Stationary', 'results/numerical/{US,US_ESC}_stationaryApprox.csv, ARG_stationaryApprox'
                  + C.argVariantTag(C.ARG['commonX']) + '.csv',
                  'Stationary against date-specific policy functions', _label('NUM_Stationary'),
                  r'>{\raggedright\arraybackslash}p{2.6cm}' + 'Y'*(ncol - 1),
                  r'\textbf{Policy} & $\rho$ & \multicolumn{' + str(ncol - 2)
                  + r'}{c}{\textbf{Stationary minus date-specific}}', _body(lines), note, width = r'\textwidth')


COUNTWORDS = {'spec': {'size': 'the paper' + "'" + 's cost', 'scale': 'the cost on the design',
                        'sizePiR': r'the $\Pi$-inclusive cost'},
              'method': {'exact': 'exact recursion', 'path': 'path iteration'}}


def _exceptionText(f, ids, nEq, fb, layer):
    """ One counted row with several equilibria or a fallback, in words: its file and identifying columns. """
    parts = []
    for k in ('ρ', 'spec', 'scenario', 'θpinned', 'method', 'effect', 'terminal', 'date'):
        if k not in ids:
            continue
        v = ids[k]
        if k == 'θpinned':
            parts.append('design pinned' if bool(v) else 'design chosen')
        elif k == 'scenario':
            parts.append(_tex(C.US['esc']['scenarios'].get(v, v)))
        elif k == 'ρ':
            parts.append(r'$\rho = ' + C.num(float(v), 1) + '$')
        else:
            parts.append(COUNTWORDS.get(k, {}).get(v, _tex(v)))
    where = 'the tax rule' if layer == 'tax' else 'the design layer'
    what = ('the solver fell back at {} state{} of {}'.format(fb, '' if fb == 1 else 's', where) if fb
            else '{} equilibria passed in {}'.format(nEq, where))
    what = what.replace(' at 1 state ', ' at one state ')
    return r'{} in {} ({})'.format(what, _tex(f), ', '.join(parts))


def _countVerdict(entries, stats):
    """ The closing sentence of both count tables: one equilibrium everywhere, or every counted exception named
    (datasets.countExceptions). """
    exc = [_exceptionText(*e, layer) for _, _, ps in entries for layer in ('tax', 'design')
           for e in D.countExceptions(ps, layer)]
    eqMax = max(s['eqMax'] for s in stats.values() if s['eqMax'] is not None)
    if eqMax == 1:      # the exceptions are then all fallbacks: states where no candidate passed
        return r'Every counted row has one equilibrium at every state and period' + (
            ', except where no candidate passed: ' + '; '.join(exc) + '.' if exc else ', and no fallback.')
    return r'The largest number of equilibria at any state is ' + str(eqMax) + ': ' + '; '.join(exc) + '.'


SOLVEGROUPS = {'calibration': 'Calibrations', 'shocks': 'Experiments', 'sweeps': 'Equilibria over designs',
               'esc': 'Endogenous design', 'oecd26': 'Cross-country test',
               **{'oecd26_theta_' + v: r'Cross-country test, $\theta$ measured as ' + C.THETANAMES.get(v, v)
                  for v in C.US['oecd26']['thetaVariants']},
               C.US['tau1']['esc']: r'Endogenous design, the cost linear in the tax',
               C.US['tau1']['oecd26']: r'Cross-country test, the cost linear in the tax',
               'numerical': 'Numerical checks'}                                  # results/<dir>, in print order
ARMWORDS = (('ARG', 'Argentina'), ('US', 'The U.S., the UK and France'))


# The notes' sentence on the roots column (datasets.ROOTCOLUMNS, countSummary's rootsMax).
ROOTSNOTE = (r'Roots: the largest number of solutions of the fixed point for next period\textquotesingle s state (the '
             r'savings of the formal households under CRRA preferences; in Argentina also the informal households\textquotesingle{} '
             r'share of savings), over the rows that record it, whatever their tax count; under log preferences the '
             r'U.S. arm solves no such fixed point, and rows written before the count existed do not record it. ')


def numSelection():
    r""" The count table by solve (page D2; appendix D states it in words): datasets.selectionCsvs' files grouped
    by arm and by the results directory that holds them (SOLVEGROUPS), the tax rule first and then the design layer
    of the exact CRRA recursion over every file that records it; per group the files, the rows, the rows that
    record counts, and over them the largest number of equilibria and of candidates, the largest number of roots of
    the state's fixed point (rootsMax, over the rows that record it) and the total of fallback states
    (datasets.countSummary over the group's files). The split of the other rows and the file-by-file list are
    numSelectionFiles'. Rows keyed `<layer>_<arm>_<dir>`. """
    entries = D.selectionCsvs()
    stats = {(lab, layer): D.countSummary(ps, layer) for _, lab, ps in entries for layer in ('tax', 'design')}
    num = lambda v: '--' if v is None else str(v)
    cells = lambda s: [str(s['files']), str(s['rows']), str(s['counted']), num(s['eqMax']), num(s['candMax']),
                       num(s['rootsMax']), num(s['fallbacks'])]
    group = lambda ps: os.path.basename(os.path.dirname(ps[0]))
    unknown = sorted({group(ps) for _, _, ps in entries} - set(SOLVEGROUPS))
    if unknown:
        raise ValueError('NUM_Selection: no label for results/{}'.format(', '.join(unknown)))
    lines = []
    for arm, words in ARMWORDS:
        lines.append(_title('Tax rule: ' + words.replace('The U.S.', 'the U.S.'), 8))
        for g, lab in SOLVEGROUPS.items():
            ps = [p for a, _, pp in entries if a == arm and group(pp) == g for p in pp]
            if ps:
                lines.append(_row([lab] + cells(D.countSummary(ps, 'tax')), 'tax_{}_{}'.format(arm.lower(), g)))
        lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')
    lines.append(_title('Design layer of the exact CRRA recursion', 8))
    for g, lab in SOLVEGROUPS.items():
        ps = [p for _, lab2, pp in entries if group(pp) == g and stats[(lab2, 'design')]['counted'] for p in pp]
        if ps:
            lines.append(_row([lab] + cells(D.countSummary(ps, 'design')), 'design_' + g))
    note = (r'\textit{Note:} The files of results that the exhibits read, grouped by the solves that wrote them: '
            r'per group the files, their rows, the rows that record the counts, and over these the largest number '
            r'of equilibria at any state and period, the largest number of candidates, and the number of states at '
            r'which no candidate passed the equilibrium test and the solver fell back. The other rows select '
            r'nothing at that layer (the economic-equilibrium rows, with taxes held at the baseline path; the tax '
            r'of the exact recursion, selected inside its design layer; the '
            r'cross-country test\textquotesingle s rows of Luxembourg at its own turnout, which has no CSES '
            r'election and so no solve; the rows with the design pinned) or were '
            r'written before the counts existed; the supplemental appendix lists every file. The tax rule\textquotesingle '
            r's candidates are the two corners and the crossings of the first order condition; the design '
            r'layer\textquotesingle s are the brackets of its consistency residual. ' + ROOTSNOTE
            + _countVerdict(entries, stats))
    header = (r' & & \multicolumn{2}{c}{\textbf{Rows}} & \multicolumn{3}{c}{\textbf{Largest number of}} & \\'
              + '\n' + r'\cmidrule(lr){3-4}\cmidrule(lr){5-7}' + '\n'
              r'\textbf{Solves} & \textbf{Files} & all & counted & equilibria & candidates & roots & \textbf{Fallbacks}')
    return _xwrap('NUM_Selection', 'results/{calibration,shocks,sweeps,esc,oecd26,numerical}/*.csv',
                  'Equilibrium counts in the results the exhibits read', _label('NUM_Selection'),
                  r'>{\raggedright\arraybackslash}p{3.0cm}YYYYYYY', header, _body(lines), note, width = r'\textwidth',
                  size = r'\footnotesize\setlength{\tabcolsep}{5pt}\renewcommand{\arraystretch}{1.1}')


def numSelectionFiles():
    r""" The equilibrium counts per results csv behind the outputs (datasets.selectionCsvs, the horizon runs of
    appendix D among them): rows; rows that record counts, rows whose solve selects nothing at that layer
    (-1), rows written before the counts existed (blank); over the counted rows the largest number of
    equilibria, the largest number of candidates and the total of fallbacks. First for the tax rule of
    writing/US/num_robustroot.tex (nEqMax, nCandMax, nFallback), then for the design layer of the exact CRRA
    recursion of writing/US/num_esc.tex (nEqθMax, nBrθMax, nFallbackθ), whose candidates are brackets. A csv
    that records no count of a layer prints '--' in its statistics; the design block lists only the csvs that
    record that layer. The note names every counted row with more than one equilibrium or a fallback
    (datasets.countExceptions). Rows keyed `<layer>_<file>`. The supplemental appendix's, beside the grouped
    numSelection. """
    entries = D.selectionCsvs()
    stats = {(lab, layer): D.countSummary(ps, layer) for _, lab, ps in entries for layer in ('tax', 'design')}
    num = lambda v: '--' if v is None else str(v)
    # a file name may break after an underscore and inside a camel-case word (\allowbreak, which the site drops),
    # so that seven count columns fit beside it
    fname = lambda lab, n: (re.sub(r'(?<=[a-z])(?=[A-Z])', r'\\allowbreak{}', _tex(lab).replace(r'\_', r'\_\allowbreak{}'))
                            + ('' if n == 1 else r' (' + str(n) + ' files)'))
    cells = lambda s: [str(s['rows']), str(s['counted']), str(s['uncounted']), str(s['unrecorded']),
                       num(s['eqMax']), num(s['candMax']), num(s['rootsMax']), num(s['fallbacks'])]
    lines = []
    for arm, title in (('ARG', 'Tax rule: Argentina'), ('US', 'Tax rule: the U.S., the UK and France')):
        lines.append(_title(title, 9))
        for a, lab, ps in entries:
            if a == arm:
                lines.append(_row([fname(lab, stats[(lab, 'tax')]['files'])] + cells(stats[(lab, 'tax')]),
                                  'tax_' + _fileKey(lab)))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')
    lines.append(_title('Design layer of the exact CRRA recursion: brackets as candidates', 9))
    for a, lab, ps in entries:
        s = stats[(lab, 'design')]
        if s['counted'] == 0:
            continue
        lines.append(_row([fname(lab, s['files'])] + cells(s), 'design_' + _fileKey(lab)))
    last = _countVerdict(entries, stats)
    note = (r'\textit{Note:} Per csv, the rows; the rows that record the counts; the rows whose solve selects '
            r'nothing at that layer and records $-1$ (no selection); the rows written before the counts were '
            r'recorded (not recorded); and over the counted rows the largest number of equilibria at any state and '
            r'period, the largest number of candidates and the total number of states at which the solver fell '
            r'back. The tax rule is the selection of the technical documentation, which falls back on the tax '
            r'that maximizes its integral criterion where no candidate passes the equilibrium test. No selection '
            r'in the tax rule: the economic-equilibrium rows, with taxes held at the baseline path; the rows of '
            r'the exact recursion under CRRA preferences with the design chosen, whose tax is selected inside the '
            r'design layer and counted there; and the '
            r'cross-country test\textquotesingle s rows of Luxembourg at its own turnout, which has no CSES election '
            r'and so no solve. '
            r'The design layer is that of the exact recursion, which compares the candidate designs with '
            r'the retirees\textquotesingle{} savings shares held fixed and then makes the shares consistent: its '
            r'candidates are the brackets of that consistency residual, and where none closes to an equilibrium it '
            r'falls back on comparing each design at its own consistent shares; no selection there are the rows '
            r'with the design pinned, and the log rows of the stationary check. Not recorded: the rows of the '
            r'earlier cost on the design and of the path iteration under CRRA preferences, the cross-check the '
            r'exhibits do not read, written before the counts existed and not rerun at the current calibration. '
            r'escExperiments.csv, which the endogenous-design exhibits read, merges the escShocks* files, whose '
            r'counts are the rows above. ' + ROOTSNOTE + last)
    return _xwrap('NUM_SelectionFiles', 'results/{calibration,shocks,sweeps,esc,oecd26,numerical}/*.csv',
                  'Equilibrium counts in the results the exhibits read, file by file', _label('NUM_SelectionFiles'),
                  r'>{\raggedright\arraybackslash}p{3.0cm}YYYYYYYY',
                  r'\textbf{File} & \multicolumn{4}{c}{\textbf{Rows}} & \multicolumn{3}{c}{\textbf{Largest number of}} '
                  r'& \textbf{Fallback} \\' + '\n' + r'\cmidrule(lr){2-5}\cmidrule(lr){6-8}' + '\n'
                  r' & all & counted & no selection & not recorded & equilibria & candidates & roots & states',
                  _body(lines), note, width = r'\textwidth', size = r'\footnotesize\setlength{\tabcolsep}{4pt}')


# ---------------------------------------------------------------------------------------------------
# Page A1: the reform of section 5 under three readings of the pre-reform economy (ARG_SpendingBaseline).
# ---------------------------------------------------------------------------------------------------
ARGREADINGS = ('paper', 'epsPre0', 'spending')


def argSpendingBaseline():
    r""" The reform of section 5 in 2010 (datasets.argReading) under three readings of the pre-reform economy, one column
    each: the paper's (the basic pension, pension spending of the workbook), no informal retiree drawing a pension
    ($\epsilon = 0$ before the reform) and the paper's at pension spending of config.ARG['spendingBaseline'] of GDP. Rows:
    the calibrated $\beta$ and $\omega$, the pre-reform tax rate and spending, its part paid to informal retirees, and the
    reform's changes in 2010 in the tax rate (p.p.), spending and spending on informal retirees (p.p. of GDP), and
    formal benefits (%). Rows keyed on the quantity. """
    r = {w: D.argReading(w) for w in ARGREADINGS}
    share = lambda x: C.num(100*x, 2)
    rows = [('beta', r'$\beta$', lambda d: C.num(d['β'], 3)),
            ('omega', r'$\omega$', lambda d: C.num(d['ω'], 3)),
            ('tax', r'Tax rate before the reform (\%)', lambda d: share(d['τ_base'])),
            ('spending', r'Pension spending before the reform (\% of GDP)', lambda d: share(d['S_base'])),
            ('informal', r'\quad paid to informal retirees', lambda d: share(d['Sinf_base'])),
            ('dtax', r'Change in the tax rate (p.p.)', lambda d: _signed(100*(d['τ_reform'] - d['τ_base']))),
            ('dspending', r'Change in pension spending (p.p.\ of GDP)', lambda d: _signed(100*(d['S_reform'] - d['S_base']))),
            ('dinformal', r'\quad paid to informal retirees',
             lambda d: _signed(100*(d['Sinf_reform'] - d['Sinf_base']))),
            ('dbenefits', r'Change in formal benefits (\%)', lambda d: _signed(100*d['dBenefits']))]
    lines = [_title('Calibration', 4)]
    for k, lab, f in rows:
        if k == 'dtax':
            lines[-1] = lines[-1].replace(r' \\ % row: ', r' \\[.5em] % row: ')
            lines.append(_title('The reform in {}'.format(C.calendar()['year0']), 4))
        lines.append(_row([lab] + [f(r[w]) for w in ARGREADINGS], k))
    S = 100*C.ARG['spendingBaseline']['spending']
    note = (r'\textit{Note:} The reform of section \ref{sec:argentina} at $\rho = ' + C.num(C.ARG['ρBaseline'], 0)
            + r'$, each column calibrated anew to the targets of appendix \ref{app:EPH}, pension spending as its '
            r'header gives it. Basic '
            r'pension: informal retirees draw the basic pension before the reform, the paper\textquotesingle s '
            r'reading; $\epsilon = 0$: none draws a pension before the reform. Spending is pension spending relative to '
            r'GDP, formal plus informal output, in ' + str(C.calendar()['year0']) + r'; the part paid to informal '
            r'retirees is their share of the outlay. Changes are reform minus pre-reform path in '
            + str(C.calendar()['year0']) + r', the formal benefits\textquotesingle{} the change of every formal '
            r'retiree\textquotesingle s benefit.')
    header = (r' & \multicolumn{2}{c}{\textbf{Spending ' + C.num(100*r['paper']['spendGDP'], 1) + r'\% of GDP}} & '
              r'\textbf{Spending ' + C.num(S, 1) + r'\% of GDP} \\' + '\n' + r'\cmidrule(lr){2-3}\cmidrule(lr){4-4}' + '\n'
              r' & \textbf{Basic pension} & $\bm{\epsilon = 0}$ & \textbf{Basic pension}')
    src = ', '.join(os.path.relpath(p, C.REPO).replace('\\', '/') for w in ARGREADINGS for p in r[w]['paths'])
    return _xwrap('ARG_SpendingBaseline', src, r'The reform under three readings of the pre-reform economy, Argentina',
                  _label('ARG_SpendingBaseline'), r'>{\raggedright\arraybackslash}p{6.2cm}YYY', header, _body(lines),
                  note, width = r'\textwidth')


# ---------------------------------------------------------------------------------------------------
# Page C2: the U.S. under the linear and the squared cost side by side (ESC_TauPower).
# ---------------------------------------------------------------------------------------------------
def escTauPower():
    r""" The U.S. under the cost linear in the tax ($k = 1$, config.US['tau1'], the cost of sections 7.1-7.2) and under
    the cost with the tax squared ($k = 2$), each calibrated to the same targets: the cost parameter at each rho of
    config.US['esc']['ρTable'] (datasets.escCalibration), and under log preferences (datasets.escTauPower) the design
    in force in 2020 in the baseline and under acute ageing and the 2020 tax rate in each, with the design chosen and
    pinned. Columns k = 1 and k = 2; rows keyed on the quantity (`lambda` the rho = 1 row, `lambda@<rho>` the others). """
    d = D.escTauPower()
    ks = sorted(d)
    cal = {k: D.escCalibration(k = k) for k in ks}
    ρ1, spec = C.US['ρAnchor'], C.US['esc']['spec']
    for k in ks:
        if abs(float(cal[k][(ρ1, spec)]['p']) - float(d[k]['cal']['p'])) > 1e-9*float(d[k]['cal']['p']):
            raise D.StaleInput('ESC_TauPower: the rho = 1 lambda of k = {} differs between its readers'.format(k))
    lam = lambda ρ: (lambda k: C.num(float(cal[k][(float(ρ), spec)]['p']), 2) if (float(ρ), spec) in cal[k] else '--')
    get = lambda k, scen, pin, col: float(d[k]['shocks'].loc[(scen, pin), col])
    rows = [(('lambda' if np.isclose(ρ, ρ1) else 'lambda@' + _rhoKey(ρ)), r'$\lambda$, $\rho = ' + C.num(ρ, 1) + '$',
             lam(ρ)) for ρ in C.US['esc']['ρTable']]
    rows += [
            ('theta_base', r'Design in force in 2020, baseline', lambda k: C.num(get(k, 'baseline', False, 'θ_t0'), 3)),
            ('theta_acute', r'\quad under acute ageing', lambda k: C.num(get(k, 'acute', False, 'θ_t0'), 3)),
            ('dtheta', r'\quad change', lambda k: _signed(get(k, 'acute', False, 'θ_t0') - get(k, 'baseline', False, 'θ_t0'), 3)),
            ('tax_base', r'Tax rate in 2020, baseline (\%)', lambda k: C.num(100*get(k, 'baseline', False, 'τ_t0'), 2)),
            ('dtax_chosen', r'\quad change under acute ageing, design chosen (p.p.)',
             lambda k: _signed(100*(get(k, 'acute', False, 'τ_t0') - get(k, 'baseline', False, 'τ_t0')))),
            ('dtax_pinned', r'\quad change under acute ageing, design pinned (p.p.)',
             lambda k: _signed(100*(get(k, 'acute', True, 'τ_t0') - get(k, 'baseline', True, 'τ_t0'))))]
    lines = [_row([lab] + [f(k) for k in ks], key) for key, lab, f in rows]
    n = len(C.US['esc']['ρTable'])
    lines[n - 1] = lines[n - 1].replace(r' \\ % row:', r' \\[.5em] % row:')
    year0 = str(C.usCalendar()['year0'])
    note = (r'\textit{Note:} The U.S.\ under the cost of section \ref{sec:esc}, $f = \exp(-\frac{1}{2}\lambda\tau^k'
            r'\tilde V(1-\theta)^2)$, linear in the tax rate ($k = 1$, equation \eqref{' + C.COSTEQ[1] + r'}, as in table '
            r'\ref{table:US_ESC:calibration}) and with the tax rate squared ($k = 2$, equation \eqref{' + C.COSTEQ[2]
            + r'}), each $\lambda$ calibrated at each $\rho$ so that the design in force in ' + year0 + r' is the '
            r'observed one, with $\omega$ and $\beta$, so that the benefits paid match spending. The other rows are under '
            r'log preferences ($\rho = 1$). Acute ageing holds workers per retiree at one over the whole horizon; each '
            r'reading is a new equilibrium path read in ' + year0 + '.')
    src = ', '.join(os.path.relpath(p, C.REPO).replace('\\', '/') for k in ks for p in d[k]['paths']
                    + (os.path.join(D.escDir(k), 'escCalibrationCRRA.csv'),))
    return _xwrap('ESC_TauPower', src, r'The U.S.\ under the linear and the squared cost',
                  _label('ESC_TauPower'), r'>{\raggedright\arraybackslash}p{7.4cm}YY',
                  [''] + [r'\textbf{' + _costName(k)[4:].capitalize() + r'}, $\bm{k = ' + str(int(k)) + '}$'
                          for k in ks],
                  _body(lines), note, width = r'\textwidth')


# ---------------------------------------------------------------------------------------------------
# The two calibrations of the taste for leisure side by side (pages A3, B3) and the two cost specifications
# (C2). Only what differs across the columns compared is printed; each builder checks that
# what it leaves out is the same, and the note says so.
# ---------------------------------------------------------------------------------------------------
VARIANTHEAD = (r'\textbf{Common} $\bm{X}$', r'\textbf{Vector} $\bm{X_i}$')


def argLeisureCalibration():
    r""" A3.1: the rows of table \ref{table:Arg:Calib} that differ between the common-X and the vector-X_i
    calibration (rho = 1), with what identifies each. A row expected to differ that does not raises, and so
    does a row left out that differs. """
    c, v = D.calibrationSummary(True), D.calibrationSummary(False)
    ww = C.num(C.calendar()['workweek'], 1)
    rows = [(r'$X$, $X_i$', '$' + C.num(c['Xi'][0], 2) + '$', '$' + C.vec(v['Xi'], 2) + '$',
             'Average formal workweek of ' + ww + r' hours (common); relative hours (vector)', 'X'),
            (r'$\eta_i$', '$' + C.vec(c['ηi'], 2) + '$', '$' + C.vec(v['ηi'], 2) + '$', 'Income distribution', 'eta'),
            (r'$\eta_0$', '$' + C.num(c['η0'], 3) + '$', '$' + C.num(v['η0'], 3) + '$', 'Informal relative income',
             'eta0'),
            (r'$X_0$', '$' + C.num(c['X0'], 3) + '$', '$' + C.num(v['X0'], 3) + '$', 'Informal relative hours', 'X0')]
    for _, a, b, _, k in rows:
        if a == b:
            raise ValueError('ARG_LeisureCalibration: {} is the same under the two calibrations'.format(k))
    differ = [k for k in ('ε', 'θ', 'α', 'ξ', 'β', 'γ0', 'ω')
              if not np.isclose(float(c[k]), float(v[k]), rtol = 1e-9, atol = 0)]
    if differ:
        raise ValueError('ARG_LeisureCalibration: {} differ between the calibrations'.format(differ))
    note = (r'\textit{Note:} $\rho = ' + C.num(c['ρ'], 0) + r'$. The other parameters of table '
            r'\ref{table:Arg:Calib} ($\epsilon$, $\theta$, $\alpha$, $\nu_t$, $\xi$, $\beta$, $\gamma_0$, '
            r'$\omega$) are the same under the two calibrations. Under the common $X$ its level is the unit of '
            r'hours, fixed by the observed formal workweek; under the vector $X_i$ relative hours are data and '
            r'the unit is a normalization.')
    return _xwrap('ARG_LeisureCalibration', 'results/paper/calibrationSummary.csv',
                  'The Argentine calibration under the two treatments of the taste for leisure',
                  _label('ARG_LeisureCalibration'), r'lcc>{\raggedright\arraybackslash}X',
                  [r'\textbf{Parameter}', VARIANTHEAD[0], VARIANTHEAD[1], r'\textbf{Identified by}'],
                  _body([_row([lab, a, b, t], k) for lab, a, b, t, k in rows]), note, width = r'\textwidth')


def argLeisureRho():
    r""" A3.2: the parameters solved with the equilibrium path that differ between the two calibrations, at the
    rho of the paper's figures (config.ARG['ρTable']): X under the common X, eta_0 and X_0 under both. The
    vector-X_i X_i are data, the same at every rho; beta and omega are the same under the two at every rho,
    which this builder checks rather than prints. """
    gc, gv = D.rhoGrid(True), D.rhoGrid(False)
    lines = []
    for ρ in C.ARG['ρTable']:
        a, b = gc[np.isclose(gc['ρ'], ρ)].iloc[-1], gv[np.isclose(gv['ρ'], ρ)].iloc[-1]
        for k in ('β', 'ω'):
            if not np.isclose(float(a[k]), float(b[k]), rtol = 1e-8, atol = 0):
                raise ValueError('ARG_LeisureRho: {} differs between the calibrations at ρ = {}'.format(k, ρ))
        lines.append(_row([C.num(ρ, 1), C.num(float(a['X']), 3), C.num(float(a['η0']), 3), C.num(float(b['η0']), 3),
                           C.num(float(a['X0']), 3), C.num(float(b['X0']), 3)], 'rho@' + _rhoKey(ρ)))
    v = D.calibrationSummary(False)
    header = (r' & \textbf{Common} & \multicolumn{2}{c}{$\bm{\eta_0}$} & \multicolumn{2}{c}{$\bm{X_0}$} \\' + '\n'
              r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}' + '\n'
              r'\textbf{CRRA} ($\rho$) & $\bm{X}$ & ' + ' & '.join([VARIANTHEAD[0], VARIANTHEAD[1]]*2))
    note = (r'\textit{Note:} Each $\rho$ is calibrated to the same targets. Under the vector $X_i$ the formal $X_i$ '
            r'are identified from the survey\textquotesingle s relative hours and are the same at every $\rho$, '
            '$' + C.vec(v['Xi'], 2) + r'$; $\beta$ and $\omega$ are the same under the two calibrations at every '
            r'$\rho$.')
    return _xwrap('ARG_LeisureRho', 'results/calibration/informalSavings_rhoGrid{,CommonX}.csv',
                  r'The leisure and productivity parameters across $\rho$, Argentina', _label('ARG_LeisureRho'),
                  'YYYYYY', header, _body(lines), note, width = r'.85\textwidth')


OECDCOLS = ('US', 'UK', 'FR')       # table 3's countries, in its order


def oecdLeisureCalibration():
    r""" B3.1: the rows of table \ref{table:US:CalibUSUKFR} that differ between the two calibrations: X (the common
    value; under vector X_i the population-weighted mean of X_i) and the productivity spread. The other rows
    are checked to be the same. """
    c, v = D.usCalibrationSummary(True), D.usCalibrationSummary(False)
    for k in OECDCOLS:
        for col in ('θ', 'ω', 'ν2020', 'β', 'ρ', 'ξ', 'α'):
            if not np.isclose(float(c[k][col]), float(v[k][col]), rtol = 1e-8, atol = 0):
                raise ValueError('OECD_LeisureCalibration: {} of {} differs between the calibrations'.format(col, k))
    cells = lambda col: [C.num(float(r[k][col]), 2) for k in OECDCOLS for r in (c, v)]
    lines = [_row([r'$X$'] + cells('Xbar'), 'X'), _row([r'$\eta_H/\eta_L$'] + cells('ηHηL'), 'etaratio')]
    header = (' & ' + ' & '.join(r'\multicolumn{2}{c}{\textbf{' + COUNTRYNAME[k] + '}}' for k in OECDCOLS) + r' \\'
              + '\n' + ''.join(r'\cmidrule(lr){%d-%d}' % (2 + 2*i, 3 + 2*i) for i in range(len(OECDCOLS))) + '\n'
              + r'\textbf{Parameter} & ' + ' & '.join([r'\textbf{Common}', r'\textbf{Vector}']*len(OECDCOLS)))
    note = (r'\textit{Note:} $\rho = 1$. Under the vector $X_i$, $X$ is the population-weighted mean of the $X_i$. '
            r'The other parameters of table \ref{table:US:CalibUSUKFR} are the same under the two calibrations.')
    return _xwrap('OECD_LeisureCalibration', 'results/paper/usCalibrationSummary.csv',
                  'The U.S., the UK and France under the two treatments of the taste for leisure',
                  _label('OECD_LeisureCalibration'), 'l' + 'Y'*(2*len(OECDCOLS)), header, _body(lines), note,
                  width = r'.9\textwidth')


def oecdLeisureGroups():
    r""" B3.2: the taste for leisure X_i and the productivity eta_i by income group in the three countries under
    the two calibrations (rho = 1). The population shares and the propensities to vote are data and the same
    under both, which this builder checks. """
    c, v = D.usCalibrationSummary(True), D.usCalibrationSummary(False)
    lines = []
    for k in OECDCOLS:
        for col in ('γi', 'μi'):
            if not np.allclose(np.asarray(c[k][col], float), np.asarray(v[k][col], float), rtol = 0, atol = 1e-12):
                raise ValueError('OECD_LeisureGroups: {} of {} differs between the calibrations'.format(col, k))
        lines.append(_title(COUNTRYNAME[k], 7))
        for col, lab, key in (('Xi', r'$X_i$', 'X'), ('ηi', r'$\eta_i$', 'eta')):
            lines.append(_row([lab] + [C.num(float(x), 2) for r in (c, v) for x in r[k][col]],
                              '{}@{}'.format(key, k.lower())))
    header = (' & ' + ' & '.join(r'\multicolumn{3}{c}{' + h + '}' for h in VARIANTHEAD) + r' \\' + '\n'
              + r'\cmidrule(lr){2-4}\cmidrule(lr){5-7}' + '\n'
              + r'\textbf{Parameter} & ' + ' & '.join([r'\textbf{Low}', r'\textbf{Medium}', r'\textbf{High}']*2))
    note = (r'\textit{Note:} $\rho = 1$. The income groups are those of tables \ref{table:a_US:CalibUS}, '
            r'\ref{table:a_US:CalibFR} and \ref{table:a_US:CalibUK}; the population shares $\gamma_i$ and the '
            r'propensities to vote $\mu_i$ are data and the same under the two calibrations.')
    return _xwrap('OECD_LeisureGroups', 'results/paper/usCalibrationSummary.csv',
                  'Leisure and productivity by income group under the two calibrations',
                  _label('OECD_LeisureGroups'), 'lYYYYYY', header, _body(lines), note, width = r'.9\textwidth')




def usStrategic():
    r""" B2.1: the strategic term of the U.S. 2020 choice per rho (datasets.strategicTerm): the elasticity of the
    savings the 2020 young carry into 2050 to the 2020 tax, the elasticity of the 2050 tax to those savings (the
    slope of the 2050 policy function) and their product as the change of the 2050 tax per unit change of the
    2020 tax. Design given, each rho at its own calibration, at the equilibrium path. """
    df = D.strategicTerm()
    year0 = C.usCalendar()['year0']
    year1 = year0 + 30
    lines = []
    for _, r in df.iterrows():
        ρ = float(r['ρ'])
        cells = [C.num(ρ, 1) + (' (log)' if np.isclose(ρ, 1.0) else ''),
                 _dash(r['elasticity s,τ'], lambda v: _signed(v, 3)),
                 _signed(r['elasticity τ1,s'], 3), _signed(r['strategic'], 3)]
        lines.append(_row(cells, 'rho@' + _rhoKey(ρ)))
    header = (r'\textbf{CRRA} ($\rho$) & Savings ' + str(year0) + r', $\frac{d\ln s_{' + str(year0) + r'}}{d\ln\tau_{'
              + str(year0) + r'}}$ & Policy ' + str(year1) + r', $\frac{d\ln\tau_{' + str(year1) + r'}}{d\ln s_{'
              + str(year0) + r'}}$ & Strategic term, $\frac{d\tau_{' + str(year1) + r'}}{d\tau_{' + str(year0) + r'}}$')
    note = (r'\textit{Note:} At the equilibrium path, the design given. The savings response is the economic '
            r'equilibrium\textquotesingle s at the savings inherited in ' + str(year0) + r', given the policies of '
            r'the later electorates; the policy slope is that of the ' + str(year1) + r' tax policy at the '
            r'savings carried into ' + str(year1) + r'. The strategic term, their product in levels, is the change '
            r'in the ' + str(year1) + r' tax rate, in percentage points, when the ' + str(year0) + r' electorate '
            r'raises its tax rate by one percentage point. Under log preferences the tax rate does not depend on '
            r'the state and the term is zero.')
    return _xwrap('US_Strategic', 'results/numerical/C5_strategicTerm.csv',
                  r'The strategic term of the U.S.\ tax choice in ' + str(year0), _label('US_Strategic'),
                  'lYYY', header, _body(lines), note, width = r'.9\textwidth')


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


def oecdSources():
    r""" The sources of figure 1's data, one row per entry of data/oecdFigure1_sources.csv and per column
    oecdFigure1.plotted derives (oecdFigure1.derivedSources: config.INEQUALITY's measure of the model's income groups):
    the columns it sources, provider and dataset with vintage and retrieval date, the series, the years and the note.
    The entries whose columns figure 1 draws come first. URLs are in the csvs, not printed. Rows keyed on the entry's
    first column. """
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
        for _, r in sub.iterrows():
            source = '; '.join(_tex(x) for x in (r['provider'] + ': ' + r['dataset'], r['vintage'],
                                                  ('retrieved ' + r['retrieved']) if r['retrieved'] else '') if x)
            series = _tex(r['series']) + ('. Years: ' + _tex(r['years']) if r['years'] else '')
            lines.append(_row([_tex(r['column']), source, series, _tex(_withoutEPS(r['note'])) or '--'],
                              _keyOf(first(r['column']))))
    extra = [first(c) for c in derived['column']]
    note = (r'\textit{Note:} The columns are those of data/oecdFigure1.csv'
            + (', but for ' + ', '.join(_tex(c) for c in extra) + r', derived from the income groups of '
               r'data/oecd29\_distribution.csv' if extra else '')
            + r'; data/oecdFigure1\_sources.csv' + (r' and data/oecd29\_sources.csv carry' if extra else ' carries')
            + r' the URL of every source. The figure\textquotesingle s Bismarckian index, $1/(2-\theta)$, is '
            r'not among them: $\theta$ is from data/oecd26\_theta.csv (data/README.md).')
    rag = r'>{\raggedright\arraybackslash}'
    text = _xwrap('OECD_Sources', 'data/oecdFigure1_sources.csv' + (', oecd29_sources.csv' if extra else ''),
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


COSTNAME = {'size': r'The paper\textquotesingle s cost', 'scale': r'The cost on the design'}


def escCostCalibration():
    r""" C2.2: the two costs of redistribution calibrated under log preferences (rho = 1, the only rho at which
    config.US['esc']['comparisonSpec'] is calibrated; escCalibration.csv through datasets.escCalibrationLOG): the
    paper's cost f = exp(-lambda tau^2 Vtilde (1-theta)^2 / 2) of eq:esc:budget and the cost on the design
    f = phi + (1-phi) theta^p, each with its parameter and the shares of revenue lost in 2020 at the observed
    design and at theta = 0 (datasets.escWedge), at the 2020 tax. The paper's row prints as the rho = 1 row of
    table:US_ESC:calibration (lambda at two decimals, the shares at one; test_onlineAppendix compares them).
    Rows keyed `spec@<spec>`. """
    ρ = C.US['ρAnchor']
    spec, alt = C.US['esc']['spec'], C.US['esc']['comparisonSpec']
    lines, θs, phi = [], [], None
    for s, symbol, digits in ((spec, r'\lambda', 2), (alt, 'p', 3)):
        r = D.escCalibrationLOG(s)
        w = D.escWedge(r, s)
        θs.append(float(r['θStar']))
        if s == alt:
            phi = float(r['phi'])
        lines.append(_row([COSTNAME[s], '$' + symbol + ' = ' + C.num(float(r['p']), digits) + '$',
                           C.pct(1 - w['fStar'], 1), C.pct(1 - w['f0'], 1)], 'spec@' + s))
    if max(θs) - min(θs) > 1e-6:
        raise ValueError('ESC_CostCalibration: theta* differs between the two costs: {}'.format(θs))
    year0 = str(C.usCalendar()['year0'])
    header = [r'\textbf{Cost}', r'\textbf{Parameter}', r'$1 - f(\theta^{\ast})$', r'$1 - f(0)$']
    note = (r'\textit{Note:} Log preferences ($\rho = ' + C.num(ρ, 0) + r'$). The paper\textquotesingle s cost is '
            r'that of equation \eqref{eq:esc:budget}, as in table \ref{table:US_ESC:calibration}; the cost on the '
            r'design is $f(\theta) = \phi + (1-\phi)\theta^{p}$ with $\phi = ' + C.num(phi, 1) + r'$ imposed. Each '
            r'parameter is calibrated so that the U.S.\ electorate re-elects its observed design, $\theta^{\ast} = '
            + C.num(θs[0], 3) + r'$, in ' + year0 + r'. $1 - f(\theta^{\ast})$ and $1 - f(0)$ are the shares of '
            r'revenue lost in ' + year0 + r' at the observed design and at a flat benefit, both at the ' + year0
            + r' tax rate.' + C.variantNote(C.US['commonX']))
    return _xwrap('ESC_CostCalibration', 'results/esc/escCalibration.csv',
                  r'The two costs of redistribution calibrated for the U.S.', _label('ESC_CostCalibration'),
                  'lYYY', header, _body(lines), note, width = r'.8\textwidth')


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


def _decompChoices(ρ, exp):
    r""" Panel B's designs at `rho`: {'US', 'FR', 'X1', 'X2'[, 'X3']: theta_2020} and the largest gap of the X rows.
    US, FR: escExperiments.csv's free baseline and frIncome rows under the paper's cost; the X rows from
    c5Numbers.py's design (rho = 1, whose US and FR rows must be escExperiments', else StaleInput) or designCRRA. """
    spec = C.US['esc']['spec']
    θ = {k: float(D.escRow(exp, ρ, spec, sc, False)['θ_t0']) for k, sc in (('US', 'baseline'), ('FR', 'frIncome'))}
    LOG = np.isclose(ρ, C.US['ρAnchor'])
    if LOG:
        dd = D.designDecomp()
        ch = dd[(dd['kind'] == 'choice') & dd['experiment'].isin(['US', 'FR', 'X1', 'X2', 'X3'])].set_index('experiment')
        for k in ('US', 'FR'):
            if k not in ch.index or abs(float(ch.loc[k, 'θ2020']) - θ[k]) > 1e-6:
                raise D.StaleInput('results/numerical/C5_designDecomp.csv: theta_2020 of {} is not escExperiments.csv\'s '
                                   '({})'.format(k, θ[k]))
        xs = ('X1', 'X2', 'X3')
    else:
        ch = D.designDecompCRRA(ρ).set_index('experiment')
        xs = ('X1', 'X2')
    miss = [x for x in xs if x not in ch.index]
    if miss:
        raise D.MissingInput('results/numerical/C5_designDecomp{}.csv (experiments {})'
                             .format('' if LOG else 'CRRA_rho{:g}'.format(ρ), miss))
    gap = max(float(ch.loc[x, 'gap']) for x in xs)
    if not gap < DECOMPGAP['LOG' if LOG else 'CRRA']:
        raise ValueError('ESC_Decompositions: the cost misses its target by {:.2e} at rho {}'.format(gap, ρ))
    θ.update({x: float(ch.loc[x, 'θ2020']) for x in xs})
    return θ, gap


def escDecompositions():
    r""" C2.3: the numbers section 7 and appendix C.2 state of the design choice, one column per rho of
    config.US['ρGrid'].

      A. Log preferences only: dW_t/dtheta_{t+1} of eq:esc:leadFOC at the U.S. choice by route (datasets.designDecomp,
         case US, piece W; DECOMPROUTES), the routes' sum (the total, zero at the choice) and the shares of the offset,
         the saving on the loss ('cost') and the tax's response ('tax') each over their sum. The routes must sum to
         the total, and the total vanish, to 1e-6 of the pull, else ValueError.
      B. theta_2020 as France's income distribution is imposed on the U.S. in three steps (DECOMPSTEPS: the cost first
         at the U.S. size and dispersion, X2, then at the U.S. size, X1, then France's own) and in total; under log
         also X3 - US, France's dispersion in the cost alone (_decompChoices).
      C. App:esc:cost's magnitude at k = 2 (datasets.lambdaElasticity; its lambda must be escCalibration's, else
         StaleInput): pi_2020, the wedge in place eq:esc:wedge, lambda, the local elasticity under the paper's carrier
         2 kappa pi (1 - pi) (columns _A), the lambda that the elasticities LAMBDAEPS imply and lambda as a multiple of
         each, the share of earnings outside the base; the local elasticity and that share under the fiscal
         externality's 3 kappa pi (1 - pi) (columns _C).

    Rows keyed `<panel>@<item>`. """
    if not C.US['commonX']:
        raise D.StaleInput('results/numerical/C5_*.csv are at common X (python/US/c5Numbers.py), not the headline')
    if C.TAU_POWER != 2 or C.CARRIER != 'convex':
        raise D.StaleInput('ESC_Decompositions: appendix C.2\'s mapping is that of k = 2 and the carrier \'convex\'')
    ρs = [float(ρ) for ρ in C.US['ρGrid']]
    ρLOG = C.US['ρAnchor']
    ncol = 1 + len(ρs)
    cal = C.usCalendar()
    year0, yearBefore = cal['year0'], cal['dates'][cal['t0'] - 1]
    lines = []
    atLOG =lambda s: [s if np.isclose(ρ, ρLOG) else '--' for ρ in ρs]

    # A. the design condition at the U.S. choice, log preferences
    dd = D.designDecomp()
    foc = dd[(dd['kind'] == 'foc') & (dd['case'] == 'US') & (dd['piece'] == 'W')].set_index('route')
    if set(r for r, _ in DECOMPROUTES) - set(foc.index):
        raise D.MissingInput('results/numerical/C5_designDecomp.csv (the routes of the U.S. design condition)')
    v = {r: float(foc.loc[r, 'value']) for r, _ in DECOMPROUTES}
    pull = abs(v['design'])
    if abs(v['design'] + v['cost'] + v['tax'] - v['total']) > 1e-6*pull or abs(v['total']) > 1e-6*pull:
        raise ValueError('ESC_Decompositions: the routes of the design condition do not sum to zero at the choice: '
                         '{}'.format(v))
    θ1 = float(foc['θ1'].iloc[0])
    lines.append(_title(r'A. The design condition at the U.S.\ choice, log preferences', ncol))
    for r, label in DECOMPROUTES:
        lines.append(_row([label] + atLOG(_signed(v[r], 4)), 'A@' + r))
    offset = v['cost'] + v['tax']
    for r, label in (('cost', 'the saving on the loss'), ('tax', 'the response of the tax')):
        lines.append(_row(['Share of the offset: ' + label] + atLOG(C.pct(v[r]/offset, 0)), 'A@share_' + r))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')

    # B. France's income distribution in three steps
    exp = D.escExperiments()
    θ, gaps = {}, {}
    for ρ in ρs:
        θ[ρ], gaps[ρ] = _decompChoices(ρ, exp)
    lines.append(_title(r'B. France\textquotesingle s income distribution in three steps, the change in $\theta_{'
                        + str(year0) + '}$', ncol))
    for key, label, (a, b) in DECOMPSTEPS:
        lines.append(_row([label] + [_signed(θ[ρ][a] - θ[ρ][b], 3) for ρ in ρs], 'B@' + key))
    lines.append(_row([r'The dispersion in the cost alone']
                      + [_signed(θ[ρ]['X3'] - θ[ρ]['US'], 3) if 'X3' in θ[ρ] else '--' for ρ in ρs], 'B@alone'))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')

    # C. the size of lambda
    le = D.lambdaElasticity().set_index('ρ')
    need = (['πBaseline2020', 'tbar2020', 'λ', 'εLocal_A', 'εLocal_C', 'baseOutside_A', 'baseOutside_C']
            + ['{}_ε{:g}_{}'.format(c, e, s) for e in LAMBDAEPS for c, s in (('λ', 'A'), ('multiple', 'A'))])
    miss = [c for c in need if c not in le.columns]
    if miss or any(not np.isclose(le.index, ρ).any() for ρ in ρs):
        raise D.MissingInput('results/numerical/C5_lambdaElasticity.csv (columns {}, rho {})'.format(miss, ρs))
    if 'carrier' in le.columns and set(le['carrier'].dropna()) != {C.CARRIER}:
        raise D.StaleInput('results/numerical/C5_lambdaElasticity.csv: carrier {}'.format(set(le['carrier'])))
    escCal = D.escCalibration()
    row = lambda ρ: le.loc[le.index[np.isclose(le.index, ρ)][0]]
    for ρ in ρs:
        λcal = float(escCal[(ρ, C.US['esc']['spec'])]['p'])
        if abs(float(row(ρ)['λ']) - λcal) > 1e-6*λcal:
            raise D.StaleInput('results/numerical/C5_lambdaElasticity.csv: lambda {} at rho {}, calibrated {}'
                               .format(float(row(ρ)['λ']), ρ, λcal))
    cells = lambda col, fmt: [fmt(float(row(ρ)[col])) for ρ in ρs]
    e0, e1 = LAMBDAEPS
    lines.append(_title(r'C. The size of $\lambda$, ' + str(year0), ncol))
    for key, label, col, fmt in (
            ('pi', r'$\pi = \Pi/R$', 'πBaseline2020', lambda x: C.num(x, 3)),
            ('tbar', r'The wedge in place, $\bar{t}$', 'tbar2020', lambda x: C.num(x, 3)),
            ('lambda', r'$\lambda$', 'λ', lambda x: C.num(x, 2)),
            ('eps', r'The local elasticity, $\varepsilon(\bar{t})$', 'εLocal_A', lambda x: C.num(x, 2)),
            ('lambda_eps{:g}'.format(e0), r'$\lambda$ at $\varepsilon = ' + C.num(e0, 1) + '$',
             'λ_ε{:g}_A'.format(e0), lambda x: C.num(x, 2)),
            ('lambda_eps{:g}'.format(e1), r'$\lambda$ at $\varepsilon = ' + C.num(e1, 1) + '$',
             'λ_ε{:g}_A'.format(e1), lambda x: C.num(x, 2)),
            ('multiple_eps{:g}'.format(e1), r'$\lambda$ as a multiple of that at $\varepsilon = ' + C.num(e1, 1) + '$',
             'multiple_ε{:g}_A'.format(e1), lambda x: C.num(x, 1)),
            ('multiple_eps{:g}'.format(e0), r'$\lambda$ as a multiple of that at $\varepsilon = ' + C.num(e0, 1) + '$',
             'multiple_ε{:g}_A'.format(e0), lambda x: C.num(x, 1)),
            ('outside', 'Earnings outside the base', 'baseOutside_A', lambda x: C.pct(x, 1)),
            ('eps_fiscal', r'Fiscal externality: $\varepsilon(\bar{t})$', 'εLocal_C', lambda x: C.num(x, 2)),
            ('outside_fiscal', r'Fiscal externality: earnings outside the base', 'baseOutside_C',
             lambda x: C.pct(x, 1))):
        lines.append(_row([label] + cells(col, fmt), 'C@' + key))

    gapCRRA = max(g for ρ, g in gaps.items() if not np.isclose(ρ, ρLOG))
    header = [''] + [r'$\rho = {:g}$'.format(ρ) for ρ in ρs]
    T = r'\mathcal{T}_{' + str(year0) + '}'
    note = (r'\textit{Note:} $\rho = ' + C.num(ρLOG, 0) + r'$ is log preferences. Panel A: the derivative of the '
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
            r'in the cost; the total is the response to France\textquotesingle s income distribution in figure '
            r'\ref{fig:US_ESC:overview}, and the split depends on the order of the steps. '
            r'Under CRRA preferences the share held is met to a relative gap of at most ' + C.pct(gapCRRA, 2) + '. '
            r'The dispersion in the cost alone is the U.S.\ with France\textquotesingle s $\tilde{V}$ in the cost, '
            r'computed under log preferences only. Panel C: the magnitudes of appendix \ref{app:esc:cost}, $\pi = \Pi/R$ '
            r'on the path without the cost, $\bar{t} = \tau(1-\pi)$ the wedge in place of eq.\ \eqref{eq:esc:wedge} at '
            r'the calibrated tax rate, $\lambda$ of table \ref{table:US_ESC:calibration}, and $\varepsilon(\bar{t}) = '
            r'\kappa\bar{t} = \lambda\tau/(2\pi)$ the local elasticity of the base under $\lambda = 2\kappa\pi(1-\pi)$. '
            r'$\lambda$ at $\varepsilon$ is the '
            r'cost parameter that a local elasticity $\varepsilon$ implies, $2\pi\varepsilon/\tau$, at the ends of the '
            r'range of participation and taxable-income elasticities \parencite{ChettyGMW11a1,SaezSG12a1}, and the '
            r'multiples are the calibrated $\lambda$ over it. Earnings outside the base: the share '
            r'$\frac{1}{2}\kappa\bar{t}^{2}$ that the wedge in place keeps outside. Fiscal externality: the budget '
            r'charged with the participation tax of every unit of base that leaves, $\lambda = 3\kappa\pi(1-\pi)$.')
    src = ('results/numerical/C5_designDecomp.csv, C5_designDecompCRRA_rho*.csv, C5_lambdaElasticity.csv, '
           'results/esc/escExperiments.csv, escCalibration{,CRRA}.csv')
    return _xwrap('ESC_Decompositions', src, r'The cost and the design: the decompositions of section \ref{sec:esc}',
                  _label('ESC_Decompositions'), '>{\\raggedright\\arraybackslash}p{6.2cm}' + 'Y'*len(ρs), header,
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
               'numerical': 'Numerical checks'}                                  # results/<dir>, in print order
ARMWORDS = (('ARG', 'Argentina'), ('US', 'The U.S., the UK and France'))


def numSelection():
    r""" The count table by solve (page D2; appendix D states it in words): datasets.selectionCsvs' files grouped
    by arm and by the results directory that holds them (SOLVEGROUPS), the tax rule first and then the design layer
    of the exact CRRA recursion over every file that records it; per group the files, the rows, the rows that
    record counts, and over them the largest number of equilibria and of candidates and the total of fallback
    states (datasets.countSummary over the group's files). The split of the other rows and the file-by-file list
    are numSelectionFiles'. Rows keyed `<layer>_<arm>_<dir>`. """
    entries = D.selectionCsvs()
    stats = {(lab, layer): D.countSummary(ps, layer) for _, lab, ps in entries for layer in ('tax', 'design')}
    num = lambda v: '--' if v is None else str(v)
    cells = lambda s: [str(s['files']), str(s['rows']), str(s['counted']), num(s['eqMax']), num(s['candMax']),
                       num(s['fallbacks'])]
    group = lambda ps: os.path.basename(os.path.dirname(ps[0]))
    unknown = sorted({group(ps) for _, _, ps in entries} - set(SOLVEGROUPS))
    if unknown:
        raise ValueError('NUM_Selection: no label for results/{}'.format(', '.join(unknown)))
    lines = []
    for arm, words in ARMWORDS:
        lines.append(_title('Tax rule: ' + words.replace('The U.S.', 'the U.S.'), 7))
        for g, lab in SOLVEGROUPS.items():
            ps = [p for a, _, pp in entries if a == arm and group(pp) == g for p in pp]
            if ps:
                lines.append(_row([lab] + cells(D.countSummary(ps, 'tax')), 'tax_{}_{}'.format(arm.lower(), g)))
        lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')
    lines.append(_title('Design layer of the exact CRRA recursion', 7))
    for g, lab in SOLVEGROUPS.items():
        ps = [p for _, lab2, pp in entries if group(pp) == g and stats[(lab2, 'design')]['counted'] for p in pp]
        if ps:
            lines.append(_row([lab] + cells(D.countSummary(ps, 'design')), 'design_' + g))
    note = (r'\textit{Note:} The files of results that the exhibits read, grouped by the solves that wrote them: '
            r'per group the files, their rows, the rows that record the counts, and over these the largest number '
            r'of equilibria at any state and period, the largest number of candidates, and the number of states at '
            r'which no candidate passed the equilibrium test and the solver fell back. The other rows select '
            r'nothing at that layer (the economic-equilibrium rows, with taxes held at the baseline path; the tax '
            r'of the exact recursion, selected inside its design layer; the cost calibrations that failed; the '
            r'cross-country test\textquotesingle s rows of Luxembourg at its own turnout, which has no CSES '
            r'election and so no solve; the rows with the design pinned) or were '
            r'written before the counts existed; the supplemental appendix lists every file. The tax rule\textquotesingle '
            r's candidates are the two corners and the crossings of the first order condition; the design '
            r'layer\textquotesingle s are the brackets of its consistency residual. '
            + _countVerdict(entries, stats))
    header = (r' & & \multicolumn{2}{c}{\textbf{Rows}} & \multicolumn{2}{c}{\textbf{Largest number of}} & \\'
              + '\n' + r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}' + '\n'
              r'\textbf{Solves} & \textbf{Files} & all & counted & equilibria & candidates & \textbf{Fallbacks}')
    return _xwrap('NUM_Selection', 'results/{calibration,shocks,sweeps,esc,oecd26,numerical}/*.csv',
                  'Equilibrium counts in the results the exhibits read', _label('NUM_Selection'),
                  r'>{\raggedright\arraybackslash}p{3.0cm}YYYYYY', header, _body(lines), note, width = r'\textwidth',
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
                       num(s['eqMax']), num(s['candMax']), num(s['fallbacks'])]
    lines = []
    for arm, title in (('ARG', 'Tax rule: Argentina'), ('US', 'Tax rule: the U.S., the UK and France')):
        lines.append(_title(title, 8))
        for a, lab, ps in entries:
            if a == arm:
                lines.append(_row([fname(lab, stats[(lab, 'tax')]['files'])] + cells(stats[(lab, 'tax')]),
                                  'tax_' + _fileKey(lab)))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.5em] % row:')
    lines.append(_title('Design layer of the exact CRRA recursion: brackets as candidates', 8))
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
            r'design layer and counted there; and cost calibrations that failed, France\textquotesingle{}s '
            r'observed design being out of the cost\textquotesingle{}s reach, and under the cost on the design '
            r'that of the UK at the U.S.\ income groups, a row escCountry.csv keeps and no exhibit prints; and the '
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
            r'counts are the rows above. ' + last)
    return _xwrap('NUM_SelectionFiles', 'results/{calibration,shocks,sweeps,esc,oecd26,numerical}/*.csv',
                  'Equilibrium counts in the results the exhibits read, file by file', _label('NUM_SelectionFiles'),
                  r'>{\raggedright\arraybackslash}p{3.0cm}YYYYYYY',
                  r'\textbf{File} & \multicolumn{4}{c}{\textbf{Rows}} & \multicolumn{2}{c}{\textbf{Largest number of}} '
                  r'& \textbf{Fallback} \\' + '\n' + r'\cmidrule(lr){2-5}\cmidrule(lr){6-7}' + '\n'
                  r' & all & counted & no selection & not recorded & equilibria & candidates & states',
                  _body(lines), note, width = r'\textwidth', size = r'\footnotesize\setlength{\tabcolsep}{4pt}')


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


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
ROLES = ('spending', 'growth', 'index', 'gini')                 # oecdFigure1's variable roles, in order
ROLEFMT = {'spending': lambda v: C.num(v, 1), 'growth': lambda v: C.num(v, 3),
           'index': lambda v: C.num(v, 3), 'gini': lambda v: C.num(v, 3)}


def oecdCountries():
    r""" Figure 1's cross-section, one row per country it plots (oecdFigure1.plotted: config.US['theta']'s measure's
    countries): the four variables oecdFigure1.PLOT draws, each followed by the year it is read at where that varies
    by country (data/oecdFigure1.csv; the index from data/oecd26_theta.csv). Rows keyed on lowercase ISO3; the
    economies section 6 calibrates (oecdFigure1.HIGHLIGHT) set in bold. """
    df = F1.plotted(D.oecdFigure1Data())
    cols = []                                         # (header, csv column, formatter)
    for role in ROLES:
        c = F1.PLOT[role]
        cols.append((r'\textbf{' + _tex(_varLabel(role, c)) + '}', c, ROLEFMT[role]))
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
    return _xwrap('OECD_Countries', 'data/oecdFigure1.csv', r'The cross-section of figure \ref{fig:US:OECDdata}, by country',
                  _label('OECD_Countries'), 'l' + 'Y'*len(cols), ['\\textbf{Country}'] + [h for h, _, _ in cols],
                  _body(lines), note, width = r'\textwidth')


EPS = re.compile(r'\bEPS\b|eps_')


def _varLabel(role, column):
    return F1.LABEL[column]


def _withoutEPS(note):
    """ A sources note without its sentences on the figure's previous version (the EPS and its eps_* readings),
    which data/README.md documents. """
    return ' '.join(s for s in re.split(r'(?<=\.)\s+(?=[A-Z])', str(note)) if not EPS.search(s))


def _plottedColumns():
    """ The csv columns figure 1 draws. """
    return set(F1.PLOT.values())


def oecdSources():
    r""" The sources of figure 1's data, one row per entry of data/oecdFigure1_sources.csv: the columns it
    sources, provider and dataset with vintage and retrieval date, the series, the years and the note. The
    entries whose columns figure 1 draws come first. URLs are in the csv, not printed. Rows keyed on the
    entry's first column. """
    src = D.oecdFigure1Sources()
    src = src[[not all(c.strip().startswith('eps_') for c in s.split(',')) for s in src['column']]]
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
    note = (r'\textit{Note:} The columns are those of data/oecdFigure1.csv; data/oecdFigure1\_sources.csv '
            r'carries the URL of every source. The figure\textquotesingle s Bismarckian index, $1/(2-\theta)$, is '
            r'not among them: $\theta$ is from data/oecd26\_theta.csv (data/README.md).')
    rag = r'>{\raggedright\arraybackslash}'
    text = _xwrap('OECD_Sources', 'data/oecdFigure1_sources.csv', r'The sources of figure \ref{fig:US:OECDdata}',
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
    label = lambda c: _tex(_varLabel(role[c], c))
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


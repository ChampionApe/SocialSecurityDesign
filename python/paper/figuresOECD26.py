r""" Figure 1.1's panels with the model beside the data (python/US/runOECD26.py, results/oecd26/): the data's
countries and fitted line in the first series colour, the model's in the second, one panel per relation of
tablesOECD26.RELATIONS. Layer 1 draws all five panels; layer 2 draws the two panels of the index. The arguments
(turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta) are tablesOECD26's, each None for the
paper's variant. Section 6's spending test (`oecd26Spending`) draws the two spending panels of layer 1 without the cost
(config.US['oecd26']['spendingSpec']) at the observed designs. The predicted against the observed design
(`oecd26ThetaFit`; the supplemental appendix's `oecd26ThetaFitRows`, the population and the data variants of page C4)
draws tablesOECD26.thetaFit's panels; page D1's `oecd26Layer1Inequality` layer 1's panels on inequality under each
measure of tablesOECD26.INEQUALITYROWS. Each figure also writes `<name>_note.tex`, its note, which build.py copies beside
the pdf for a figure the paper inputs and the site shows for one it does not (onlineAppendix.generatedNote).
Layer 1 prints as "the model" (legend `model`, panel title `Model`), layer 2 as "own ω" (`Own $\omega$`), whose first
definition is OECD26_ThetaFit's note (tablesOECD26's NAMES).
"""
import os
import numpy as np
from matplotlib.transforms import Bbox

import config as C
import datasets as D
import oecdFigure1 as F1
from figures import plt, SERIES, INK, _panel, _save, markId, mark
from tables import BANNER
import tablesOECD26 as T
from tablesOECD26 import (RELATIONS, SPENDING, panelData, ols, suffix, variant, spendingSlopes, _xy, _layersText,
                          modelText, sampleText, turnoutNote)

MARKFMT = dict(F1.MARKFMT)
AXIS = {'spending': (0, 20), 'index': (0.44, 1.07)}
EDGE = '#fcfcfb'
LEGEND = {1: 'model', 2: r'model, own $\omega$'}                 # the model's legend entry, by layer


def _src(theta = None, layers = 'layer1.csv, layer2.csv, anchor.csv'):
    """ A note's %% Source: line: the csvs of the measure `theta`'s results directory (datasets.oecd26Dir) and the data. """
    return D.oecd26DirRel(theta) + layers + '; ' + T.DATASRC


def _note(name, src, text):
    """ Writes results/paper/Figs/<name>_note.tex (the banner and `text`) and returns its path. """
    path = os.path.join(C.PAPERDIR, 'Figs', name + '_note.tex')
    with open(path, 'w', encoding = 'utf-8', newline = '\n') as f:
        f.write(BANNER.format(name = name, src = src) + text + '\n')
    return path


def oecd26Figure(layer = 1, turnout = None, sample = None, reading = None, spec = None, betaRule = None,
                 spending = None, omegaRule = None, retAge = None, theta = None, slopesTable = True):
    r""" OECD26_Layer1 (layer 1) or its layer-2 counterpart OECD26_Layer2: 2x3 as figure 1.1, the legend in
    the sixth cell (layer 1), or 1x3 with the index panels (layer 2). Every dot and line is a mark; a dot is
    linked to its country's row of OECD26_Omega (lowercase ISO3), a line carries its slope and, `slopesTable`, the row
    of the variant's OECD26_Slopes (False where that table is not built). Layer 1 also writes its note
    (`OECD26_Layer1<suffix>_note.tex`): the sample, the model's index and spending, its values, the turnout; it names
    no layer and, at the paper's variant, not own ω either (figure 7.3 defines it). """
    v = variant(turnout, sample, reading, spec, betaRule, spending, omegaRule, retAge, theta)
    p = panelData(layer, **v)
    sfx = suffix(**v)
    name = 'OECD26_Layer{}'.format(layer) + sfx
    cells = ([(0, 0, 'spending', 'growth'), (0, 1, 'spending', 'inequality'), (0, 2, 'spending', 'index'),
              (1, 0, 'index', 'growth'), (1, 1, 'index', 'inequality')] if layer == 1 else
             [(0, 0, 'index', 'growth'), (0, 1, 'index', 'inequality')])
    rows = 2 if layer == 1 else 1
    fig, axes = plt.subplots(rows, 3, figsize = (5.91, 3.9 if layer == 1 else 2.15), constrained_layout = True,
                             squeeze = False)
    marks = []
    for r, c, yv, xv in cells:
        if (yv, xv) not in RELATIONS:
            raise ValueError('{}: ({}, {}) is not one of tablesOECD26.RELATIONS'.format(name, yv, xv))
        ax = axes[r, c]
        _panel(ax, '', '', titlesize = 8, labelsize = 7)
        (xn, xf), (yn, yf) = MARKFMT[xv], MARKFMT[yv]
        for side, colour, z in (('Data', SERIES[0], 3), ('Model', SERIES[1], 4)):
            x, y = _xy(p, yv, xv, side)
            for iso in p.index:
                pt, = ax.plot(x[iso], y[iso], linestyle = 'none', marker = 'o' if side == 'Data' else 'D',
                              markersize = 3.6, color = colour, markeredgecolor = EDGE,
                              markeredgewidth = 0.5, zorder = z)
                pt.set_gid(markId(name, yv + '_' + xv, side.lower(), iso.lower()))
                marks.append(mark(pt.get_gid(), p.at[iso, 'country'], y[iso],
                                  '{}: {}; {}: {}'.format(xn, xf.format(x[iso]), yn, yf.format(y[iso])),
                                  series = side.lower(), panel = '{} against {}'.format(yn, xn),
                                  table = 'OECD26_Omega' + sfx if layer == 2 else None,
                                  row = iso.lower() if layer == 2 else None))
            o = ols(x, y)
            xs = np.linspace(float(x.min()), float(x.max()), 2)
            ln, = ax.plot(xs, o['intercept'] + o['slope']*xs, color = colour, linewidth = 1.2, zorder = z - 1)
            ln.set_gid(markId(name, yv + '_' + xv, side.lower(), 'fit'))
            marks.append(mark(ln.get_gid(), 'fitted line', o['slope'],
                              'slope {:.{d}f} (s.e. {:.{d}f})'.format(o['slope'], o['se'], d = T.slopeDigits(yv, xv)),
                              series = side.lower(),
                              panel = '{} against {}'.format(yn, xn),
                              table = 'OECD26_Slopes' + sfx if slopesTable else None,
                              row = yv + '_' + xv if slopesTable else None))
        if c == 0:
            ax.set_ylabel(F1.LABEL[F1.PLOT[yv]], color = INK['secondary'], fontsize = 7.5)
        if r == rows - 1 or (layer == 1 and c == 2):
            ax.set_xlabel(F1.LABEL[F1.PLOT[xv]], color = INK['secondary'], fontsize = 7.5)
        if yv == 'index':
            ax.set_ylim(*AXIS['index'])
        else:
            ax.set_ylim(0, max(AXIS['spending'][1], float(np.nanmax(p[['spendingData', 'spendingModel']].values)) + 1))
        if xv == 'index':
            ax.set_xlim(*AXIS['index'])
    legendAx = axes[rows - 1, 2]
    legendAx.axis('off')
    dot = lambda colour, marker: plt.Line2D([], [], linestyle = '-', linewidth = 1.2, marker = marker,
                                            markersize = 3.6, color = colour, markeredgecolor = EDGE,
                                            markeredgewidth = 0.5)
    n = len(p)
    leg = legendAx.legend([dot(SERIES[0], 'o'), dot(SERIES[1], 'D')],
                          ['data', LEGEND[layer]], loc = 'center', frameon = False, fontsize = 7,
                          labelcolor = INK['secondary'], title = '{} countries, fitted lines'.format(n),
                          title_fontsize = 7)
    leg.get_title().set_color(INK['secondary'])
    out = _save(fig, name, marks)
    if layer != 1:
        return out
    a = D.oecd26Anchor(v['turnout'], v['spec'], T._thetaSource(v))
    if sfx == '':
        # the anchor of the OECD26 notes: the later exhibits point here (tablesOECD26.specText)
        T.paperWording()
        ω = T.layer1Omega(v)
        shared = (r'every country at the U.S.\ $\beta$ and at $\omega = ' + C.num(ω, 2)
                  + '$, ' + T._meanWords(v) + r' (appendix \ref{app:oecd26:layer2})' if ω is not None else
                  r'every country at the U.S.\ $\beta$ and $\omega$')
        note = (r'\textit{Note:} ' + '{} countries'.format(n) + r', as in figure \ref{fig:OECD26:spending}, in the data '
                r'and in the model (section \ref{sec:esc:crossSection}): ' + shared + r', at its cost '
                r'$\lambda_{n,t}$ and at its own CSES propensities to vote, Luxembourg, which has none, at a flat '
                r'profile. The model\textquotesingle s index is $1/(2-\theta)$ of the design chosen in 1990 for '
                r'2020, its spending $\tau(1-\alpha)$ in 2020 on the path of that design. The lines are least-squares '
                r'fits; table \ref{table:OECD26:slopes} reports their slopes.')
    else:
        note = (r'\textit{Note:} ' + sampleText(n, p.attrs['dropped']) + r' in the data and in the model, '
                r'with the least-squares lines whose slopes table~\ref{table:OECD26:slopes' + sfx + '} reports. '
                + modelText(v) + ' ' + _layersText(a, v, layers = (1,)) + ' ' + turnoutNote(v['turnout']))
    return out + [_note(name, _src(T._thetaSource(v)), note)]


def oecd26Spending(turnout = None, sample = None):
    r""" OECD26_Spending (section 6's spending test): pension spending on population growth and on income inequality
    (config.INEQUALITY), the data
    and layer 1 at config.US['oecd26']['spendingSpec'] (no cost) and the β rule 'common', at the observed designs;
    1x3 as layer 2's figure, the legend in the third cell. Also writes `OECD26_Spending_note.tex`, the figure's note
    (the sample, the anchor's β and ω without the cost, the turnout; the slopes are in section 6's text), which build.py copies beside the pdf for the paper's tablenotes to \input.
    Every dot and line is a mark, a line carrying its slope. """
    spec = C.US['oecd26']['spendingSpec']
    p = panelData(1, turnout, sample, None, spec, 'common', spending = 'observed')
    s = spendingSlopes(turnout, sample).set_index(['y', 'x'])
    sfx = suffix(turnout, sample)
    name = 'OECD26_Spending' + sfx
    fig, axes = plt.subplots(1, 3, figsize = (5.91, 2.15), constrained_layout = True, squeeze = False)
    marks = []
    for c, (yv, xv) in enumerate(SPENDING):
        ax = axes[0, c]
        _panel(ax, '', '', titlesize = 8, labelsize = 7)
        (xn, xf), (yn, yf) = MARKFMT[xv], MARKFMT[yv]
        for side, colour, z in (('Data', SERIES[0], 3), ('Model', SERIES[1], 4)):
            x, y = _xy(p, yv, xv, side)
            for iso in p.index:
                pt, = ax.plot(x[iso], y[iso], linestyle = 'none', marker = 'o' if side == 'Data' else 'D',
                              markersize = 3.6, color = colour, markeredgecolor = EDGE,
                              markeredgewidth = 0.5, zorder = z)
                pt.set_gid(markId(name, yv + '_' + xv, side.lower(), iso.lower()))
                marks.append(mark(pt.get_gid(), p.at[iso, 'country'], y[iso],
                                  '{}: {}; {}: {}'.format(xn, xf.format(x[iso]), yn, yf.format(y[iso])),
                                  series = side.lower(), panel = '{} against {}'.format(yn, xn)))
            o = s.loc[(yv, xv), 'data' if side == 'Data' else 'layer1']
            xs = np.linspace(float(x.min()), float(x.max()), 2)
            ln, = ax.plot(xs, o['intercept'] + o['slope']*xs, color = colour, linewidth = 1.2, zorder = z - 1)
            ln.set_gid(markId(name, yv + '_' + xv, side.lower(), 'fit'))
            marks.append(mark(ln.get_gid(), 'fitted line', o['slope'],
                              'slope {:.{d}f} (s.e. {:.{d}f})'.format(o['slope'], o['se'], d = T.slopeDigits(yv, xv)),
                              series = side.lower(),
                              panel = '{} against {}'.format(yn, xn)))
        if c == 0:
            ax.set_ylabel(F1.LABEL[F1.PLOT[yv]], color = INK['secondary'], fontsize = 7.5)
        ax.set_xlabel(F1.LABEL[F1.PLOT[xv]], color = INK['secondary'], fontsize = 7.5)
        ax.set_ylim(0, max(AXIS['spending'][1], float(np.nanmax(p[['spendingData', 'spendingModel']].values)) + 1))
    legendAx = axes[0, 2]
    legendAx.axis('off')
    dot = lambda colour, marker: plt.Line2D([], [], linestyle = '-', linewidth = 1.2, marker = marker,
                                            markersize = 3.6, color = colour, markeredgecolor = EDGE,
                                            markeredgewidth = 0.5)
    leg = legendAx.legend([dot(SERIES[0], 'o'), dot(SERIES[1], 'D')], ['data', 'model'], loc = 'center',
                          frameon = False, fontsize = 7, labelcolor = INK['secondary'],
                          title = '{} countries, fitted lines'.format(len(p)), title_fontsize = 7)
    leg.get_title().set_color(INK['secondary'])
    out = _save(fig, name, marks)
    v = variant(turnout, sample, None, spec, 'common', 'observed')
    a = D.oecd26Anchor(v['turnout'], spec)
    ω = T.layer1Omega(v)
    # the mean ω here is that of the weights calibrated without the cost (section 6), not figure 7.2's
    weight = (r'$ and at $\omega = ' + C.num(ω, 2) + r'$, the mean of the political weights that reproduce each '
              r'country\textquotesingle s spending without the cost' if ω is not None
              else r'$ and $\omega = ' + C.num(float(a['ω']), 2) + '$, without the cost')
    note = (r'\textit{Note:} ' + sampleText(len(p), p.attrs['dropped']) + r'. The model\textquotesingle s spending is '
            r'$\tau(1-\alpha)$ in 2020 at the observed design, every country at the U.S.\ $\beta = '
            + C.num(float(a['β']), 2) + weight + '. ' + turnoutNote(v['turnout'])
            + ' The lines are least-squares fits.')
    return out + [_note(name, _src(layers = 'layer1.csv, anchor.csv'), note)]


# ---------------------------------------------------------------------------------------------------
# Page D1's section on inequality: figure 7.2's panels on inequality under each measure of tablesOECD26.INEQUALITYROWS.
INEQUALITYSOURCE = {'gini': 'WIID', 'yHighOverLow': 'OECD earnings, men full-time',
                    'Vtilde': 'OECD earnings, men full-time'}


def yearsText(s):
    """ '2019' or '2019 (POL 2018)': the modal year of `s` (on ISO3) and the countries read at another. """
    s = s.astype(int)
    mode = int(s.mode().iloc[0])
    other = s[s != mode]
    return str(mode) + (' ({})'.format(', '.join('{} {}'.format(c, v) for c, v in other.items())) if len(other) else '')


def oecd26Layer1Inequality():
    r""" OECD26_Layer1Inequality (page D1, ONLINE_ONLY): figure 7.2's spending and index panels on income inequality,
    one row per measure of tablesOECD26.INEQUALITYROWS (the Gini, y_3/y_1, Ṽ), the data and the model at the paper's
    variant, the least-squares line of each, and in each panel the data's and the model's slope (s.e.) and Pearson r.
    Writes results/paper/oecd26Layer1InequalitySlopes.csv (tablesOECD26.inequalitySlopes, which
    OECD26_Layer1InequalitySlopes prints) and its note. Every dot and line is a mark, a line linked to its row of that
    table. """
    name = 'OECD26_Layer1Inequality'
    p = panelData(1)
    s = T.inequalitySlopes(p)
    os.makedirs(C.PAPERDIR, exist_ok = True)
    csv = D.paperCsv(D.OECD26INEQSLOPE)
    s.to_csv(csv, index = False)
    years = {'gini': D.oecdFigure1Data().loc[p.index, 'gini_year'], 'yHighOverLow': D.incomeGroups(p.index)['year']}
    years['Vtilde'] = years['yHighOverLow']
    rows = T.INEQUALITYROWS
    H = 2.6*len(rows) + 0.9
    fig, axes = plt.subplots(len(rows), 2, figsize = (5.91, H))
    fig.subplots_adjust(left = 0.10, right = 0.985, top = 1 - 0.32/H, bottom = 0.95/H, wspace = 0.22, hspace = 0.48)
    marks = []
    for r, m in enumerate(rows):
        x = p[m].astype(float)
        lab = C.INEQUALITYLABELS[m]
        xn, xf = lab['mark'], lab['fmt']
        for c, yv in enumerate(('spending', 'index')):
            ax = axes[r, c]
            _panel(ax, '', '', labelsize = 6.5)
            ax.tick_params(length = 2.5)
            if c == 0:
                ax.set_title('({}) {}{}\n'.format('abc'[r], lab['axis'], " (the paper's)" if m == C.INEQUALITY else ''),
                             color = INK['primary'], fontsize = 8, loc = 'left', pad = 3)
                ax.text(0, 1.0, INEQUALITYSOURCE[m] + ', ' + yearsText(years[m]), transform = ax.transAxes,
                        ha = 'left', va = 'bottom', fontsize = 6.3, color = INK['secondary'])
            ax.set_ylabel(F1.LABEL[F1.PLOT[yv]], color = INK['secondary'], fontsize = 6.5)
            ax.set_xlabel(lab['axis'], color = INK['secondary'], fontsize = 6.5, labelpad = 1)
            yn, yf = MARKFMT[yv]
            txt = []
            for side, colour, marker, z in (('Data', SERIES[0], 'o', 3), ('Model', SERIES[1], 'D', 4)):
                y = p[yv + side]
                for iso in p.index:
                    pt, = ax.plot(x[iso], y[iso], ls = 'none', marker = marker, ms = 3.2, color = colour, mec = EDGE,
                                  mew = 0.5, zorder = z)
                    pt.set_gid(markId(name, m, yv, side.lower(), iso.lower()))
                    marks.append(mark(pt.get_gid(), p.at[iso, 'country'], y[iso],
                                      '{}: {}; {}: {}'.format(xn, xf.format(x[iso]), yn, yf.format(y[iso])),
                                      series = side.lower(), panel = '{} against {}'.format(yn, xn)))
                o = s[(s['measure'] == m) & (s['y'] == yv) & (s['side'] == side.lower())].iloc[0]
                xs = np.array([float(x.min()), float(x.max())])
                ln, = ax.plot(xs, o['intercept'] + o['slope']*xs, color = colour, lw = 1.1, zorder = z - 1)
                ln.set_gid(markId(name, m, yv, side.lower(), 'fit'))
                d = 3 if yv == 'index' else 2                       # as OECD26_Layer1InequalitySlopes prints them
                marks.append(mark(ln.get_gid(), 'fitted line', o['slope'],
                                  'slope {:.{d}f} (s.e. {:.{d}f}), r {:.2f}'.format(o['slope'], o['se'], o['r'], d = d),
                                  series = side.lower(), panel = '{} against {}'.format(yn, xn),
                                  table = 'OECD26_Layer1InequalitySlopes', row = yv + '_' + m.lower()))
                txt.append((colour, '{} {:.{d}f} ({:.{d}f}), r {:.2f}'.format(LEGEND[1] if side == 'Model' else 'data',
                                                                             o['slope'], o['se'], o['r'], d = d)))
            if yv == 'index':
                ax.set_ylim(AXIS['index'][0], AXIS['index'][1] + 0.2)
            else:
                ax.set_ylim(0, max(AXIS['spending'][1], float(np.nanmax(p[['spendingData', 'spendingModel']].values)) + 1)
                            + 7)
            for k, (colour, t) in enumerate(txt):
                ax.text(0.02, 0.97 - 0.07*k, t, transform = ax.transAxes, ha = 'left', va = 'top', fontsize = 5.8,
                        color = colour, zorder = 7)
    h = [plt.Line2D([], [], color = SERIES[0], lw = 1.1, marker = 'o', ms = 3.2, mec = EDGE, mew = 0.5),
         plt.Line2D([], [], color = SERIES[1], lw = 1.1, marker = 'D', ms = 3.2, mec = EDGE, mew = 0.5)]
    leg = fig.legend(h, ['data', LEGEND[1]], loc = 'lower center', ncol = 2, frameon = False, fontsize = 6.5,
                     labelcolor = INK['secondary'], title = '{} countries, fitted lines; in each panel: slope (s.e.), '
                                                            'Pearson r'.format(len(p)), title_fontsize = 6.5)
    leg.get_title().set_color(INK['secondary'])
    out = _save(fig, name, marks)
    paper = lambda k: r', the measure of figure \ref{fig:OECD26:layer1}' if k == C.INEQUALITY else ''
    note = (r'\textit{Note:} ' + '{} countries'.format(len(p)) + r', the data and the model of figure '
            r'\ref{fig:OECD26:layer1}: spending and the index against three measures of income inequality, the Gini of '
            r'disposable income (WIID)' + paper('gini') + r'; high over low income, $y_3/y_1$' + paper('yHighOverLow')
            + r'; and the dispersion $\tilde V = \sum_i \gamma_i (y_i - 1)^2/y_i$' + paper('Vtilde')
            + r', the last two of the income groups the model is fed. The lines are least-squares fits; table '
            r'\ref{table:OA:OECD26_Layer1InequalitySlopes} reports their slopes.')
    return out + [csv, _note(name, _src(layers = 'layer1.csv, anchor.csv'), note)]


# ---------------------------------------------------------------------------------------------------
# The predicted against the observed design (OECD26_ThetaFit and its index version; page C4's
# OECD26_ThetaFitPopulation and OECD26_ThetaFitMeasure, ONLINE_ONLY).
UNITS = {'theta': (r'$\theta$', (-0.04, 1.04), 'θ'), 'index': (r'index $1/(2-\theta)$', (0.47, 1.02), 'I')}
NUTITLE = {'none': r'$\nu$ uncorrected', 'nra': r'$\nu$ at the normal age', 'exit': r'$\nu$ at the exit age'}


def layerTitle(layer, omegaRule = None):
    """ A panel's title: layer 2 `Own ω`; layer 1 `Model`, with its ω rule only where that is not
    config.US['oecd26']'s (the paper's model; None for it). """
    if layer == 2:
        return T.LAYERNAME[2]
    rule = C.US['oecd26']['omegaRule'] if omegaRule is None else omegaRule
    if rule == C.US['oecd26']['omegaRule']:
        return T.LAYERNAME[1]
    return T.LAYERNAME[1] + (r': mean $\omega$' if rule == 'mean' else r': U.S. $\omega$')


def fitRows(kind):
    """ The rows of page C4's figure `kind` ('population' | 'measure'): (title, thetaFit's keywords), the paper's
    marked: the readings of ν (tablesOECD26.populationRows) or the measures of the design (tablesOECD26.MEASURES). """
    paper = lambda on: " (the paper's)" if on else ''
    if kind == 'population':
        return [(NUTITLE[r] + paper(r == T.PAPER['retAge']), {'retAge': r}) for r in T.populationRows()]
    return [(C.THETANAMES[m][0].upper() + C.THETANAMES[m][1:] + paper(m == T.PAPER['theta']), {'theta': m})
            for m in T.MEASURES]


def _overlap(a, b):
    return max(0., min(a.x1, b.x1) - max(a.x0, b.x0))*max(0., min(a.y1, b.y1) - max(a.y0, b.y0))


def _isoLabels(ax, xs, ys, texts, bold = (), fontsize = 5.6):
    """ Greedy placement of point labels: the first of twelve offsets that hits no placed label, no other marker and
    stays inside the axes, else the least overlapping one; the most crowded points are placed first. """
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.c_[xs, ys])
    rad = 2.6*fig.dpi/72
    dots = [Bbox.from_extents(px - rad, py - rad, px + rad, py + rad) for px, py in pts]
    box = ax.get_window_extent(rend)
    offs = [(3, 2, 'left', 'bottom'), (3, -2, 'left', 'top'), (-3, 2, 'right', 'bottom'), (-3, -2, 'right', 'top'),
            (0, 4, 'center', 'bottom'), (0, -4, 'center', 'top'), (5, 0, 'left', 'center'), (-5, 0, 'right', 'center'),
            (8, 6, 'left', 'bottom'), (-8, -6, 'right', 'top'), (8, -6, 'left', 'top'), (-8, 6, 'right', 'bottom')]
    placed = []
    order = np.argsort([-sum(_overlap(dots[i].expanded(4, 4), m) > 0 for m in dots) for i in range(len(xs))],
                       kind = 'stable')
    for i in order:
        best, score0 = None, None
        for dx, dy, ha, va in offs:
            t = ax.annotate(texts[i], (xs[i], ys[i]), xytext = (dx, dy), textcoords = 'offset points', ha = ha,
                            va = va, fontsize = fontsize, color = INK['primary'] if texts[i] in bold else INK['secondary'],
                            fontweight = 'bold' if texts[i] in bold else 'normal', zorder = 6)
            t.set_in_layout(False)
            bb = t.get_window_extent(rend).expanded(1.04, 1.08)
            score = sum(_overlap(bb, q) for q in placed) + sum(_overlap(bb, m) for j, m in enumerate(dots) if j != i)
            if not (box.x0 <= bb.x0 and bb.x1 <= box.x1 and box.y0 <= bb.y0 and bb.y1 <= box.y1):
                score += 1e5
            if best is None or score < score0:
                if best is not None:
                    best.remove()
                best, score0 = t, score
            else:
                t.remove()
            if score == 0:
                break
        placed.append(best.get_window_extent(rend))


def _statText(f):
    return 'r = {:.2f}\nRMSE {:.2f}  (at the mean: {:.2f})'.format(f['r'], f['rmse'], f['rmseMean'])


def _fitPanel(ax, f, name, part, units, lim, ms = 3.6, msUS = 5.2):
    """ One predicted-against-observed panel of thetaFit's `f`: the 45-degree line, the fitted line over the countries
    other than the U.S., the countries, the U.S. open, the statistics. Returns the marks. """
    sym = UNITS[units][2]
    _panel(ax, '', '', titlesize = 7.5, labelsize = 6.5)
    ax.tick_params(length = 2.5)
    ax.grid(True, color = INK['grid'], linewidth = 0.6)
    ax.plot(lim, lim, color = INK['muted'], lw = 0.8, ls = (0, (3, 2)), zorder = 1)
    x, y, p = f['x'], f['y'], f['panel']
    nu = [c for c in p.index if c != 'USA']
    o = f['fit']
    xs = np.array([float(x[nu].min()), float(x[nu].max())])
    marks = []
    ln, = ax.plot(xs, o['intercept'] + o['slope']*xs, color = SERIES[1], lw = 1.0, zorder = 2)
    ln.set_gid(markId(name, part, 'fit'))
    marks.append(mark(ln.get_gid(), 'fitted line, without the U.S.', o['slope'], 'slope {:.2f}'.format(o['slope']),
                      series = 'model', panel = part))
    for iso in p.index:
        us = iso == 'USA'
        pt, = ax.plot(x[iso], y[iso], 'o' if us else 'D', ms = msUS if us else ms, ls = 'none',
                      zorder = 4 if us else 3,
                      **({'mfc': 'white', 'mec': INK['primary'], 'mew': 0.9} if us else
                         {'color': SERIES[1], 'mec': EDGE, 'mew': 0.5}))
        pt.set_gid(markId(name, part, iso.lower()))
        marks.append(mark(pt.get_gid(), p.at[iso, 'country'] + (' (anchor)' if us else ''), y[iso],
                          'observed {s} {:.3f}; predicted {s} {:.3f}'.format(x[iso], y[iso], s = sym),
                          series = 'model', panel = part))
    tx = ax.text(0.03, 0.97, _statText(f), transform = ax.transAxes, ha = 'left', va = 'top', fontsize = 6.2,
                 color = INK['secondary'], linespacing = 1.35, zorder = 7)
    tx.set_gid(markId(name, part, 'stats'))
    marks.append(mark(tx.get_gid(), 'Pearson r, without the U.S.', f['r'],
                      'r {:.2f}; RMSE {:.2f}, at the mean {:.2f}'.format(f['r'], f['rmse'], f['rmseMean']),
                      series = 'model', panel = part))
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect('equal', adjustable = 'box')
    return marks


def oecd26ThetaFit(units = 'theta'):
    r""" OECD26_ThetaFit (section 7.3; `units` 'index': OECD26_ThetaFitIndex, built, not input): the design each
    country's model predicts, in force in 2020, against the observed one, the model | own ω (layers 1 | 2), at the
    paper's variant; the 45-degree line, the fitted line over the countries other than the U.S., the U.S. open (the
    anchor, on the line by construction), ISO labels, and r and the RMSE with that of predicting the sample mean, the
    U.S. left out (tablesOECD26.thetaFit). Writes its note, which defines own ω; the statistics are on the panels. """
    name = 'OECD26_ThetaFit' + ('' if units == 'theta' else 'Index')
    unit, lim, _ = UNITS[units]
    T.paperWording()
    fs = {k: T.thetaFit(k, units = units) for k in (1, 2)}
    fig, axes = plt.subplots(1, 2, figsize = (5.91, 3.05))
    fig.subplots_adjust(left = 0.085, right = 0.99, bottom = 0.25, top = 0.93, wspace = 0.2)
    marks = []
    for ax, k in zip(axes, (1, 2)):
        f = fs[k]
        part = 'layer{}'.format(k)
        marks += _fitPanel(ax, f, name, part, units, lim)
        ax.set_title(layerTitle(k), fontsize = 7.5, color = INK['primary'], loc = 'left')
        ax.set_xlabel('Observed ' + unit, color = INK['secondary'], fontsize = 7)
        if k == 1:
            ax.set_ylabel('Predicted ' + unit, color = INK['secondary'], fontsize = 7)
        _isoLabels(ax, f['x'].to_numpy(float), f['y'].to_numpy(float), list(f['panel'].index), bold = ('USA',))
    h = [plt.Line2D([], [], ls = 'none', marker = 'D', ms = 3.6, color = SERIES[1], mec = EDGE, mew = 0.5),
         plt.Line2D([], [], ls = 'none', marker = 'o', ms = 5.2, mfc = 'white', mec = INK['primary'], mew = 0.9),
         plt.Line2D([], [], color = SERIES[1], lw = 1.0), plt.Line2D([], [], color = INK['muted'], lw = 0.8,
                                                                    ls = (0, (3, 2)))]
    mean = C.US['oecd26']['omegaRule'] == 'mean'           # the U.S. then leaves the 45-degree line in the model
    fig.legend(h, ['country, model', 'U.S. (anchor{})'.format(r'; on the line with own $\omega$ by construction' if mean
                                                             else ', on the line by construction'),
                   'fitted line', '45° line'],
               loc = 'lower center', ncol = 4, frameon = False, fontsize = 6.3, labelcolor = INK['secondary'],
               handlelength = 1.8, columnspacing = 1.4)
    out = _save(fig, name, marks)
    what = r'$\theta$' if units == 'theta' else r'index $1/(2-\theta)$'
    # own ω is defined here, where the paper first shows it (tablesOECD26's NAMES)
    note = (r'\textit{Note:} The design ' + what + r' of each country\textquotesingle s model for 2020, chosen in 1990, against '
            r'the observed one. Left, the model of figure \ref{fig:OECD26:layer1}; right, the model with each '
            r'country\textquotesingle s own political weight (own $\omega$), calibrated to its spending at its observed '
            r'design (appendix \ref{app:oecd26:layer2}). The U.S.\ is the anchor'
            + (r'; in the right panel it lies on the $45^{\circ}$ line by construction' if mean else
               r' and lies on the $45^{\circ}$ line by construction')
            + r'. The solid line is the least-squares fit over the other ' + str(fs[1]['nStat'])
            + r' countries, and each panel reports, without the U.S., the correlation $r$ and the RMSE of the '
            r'prediction, beside the RMSE of a prediction at the sample mean.')
    return out + [_note(name, _src(), note)]


# The measures are defined on the site (page C4's section esc-oecd26-measure), which the paper has no label for: \oa.
FITNOTE = {
    'population': (r'with $\nu$ corrected for each country\textquotesingle s gap to the U.S.\ normal retirement age '
                   r'(the paper\textquotesingle s), uncorrected, or corrected for its gap in the effective age of '
                   r'labor-market exit'),
    'measure':    (r'under each measure of $\theta$ (\oa{esc-oecd26-measure}), each its own cross-section with '
                   r'the U.S.\ anchor recalibrated under it, at the paper\textquotesingle s $\nu$')}


def oecd26ThetaFitRows(kind, units = 'theta'):
    r""" Page C4's OECD26_ThetaFitPopulation (`kind` 'population') or OECD26_ThetaFitMeasure ('measure'):
    thetaFit's panels without labels on common axes, one row per reading of fitRows(kind), the paper's among them,
    the model | own ω (layers 1 | 2) at the paper's other settings (a data variant from its results/oecd26_theta_<V>/).
    Writes its note; the statistics are on the panels. """
    name = 'OECD26_ThetaFit' + kind.capitalize() + ('' if units == 'theta' else 'Index')
    unit, lim, _ = UNITS[units]
    rows = fitRows(kind)
    H = 0.35 + 2.7*len(rows)
    fig, axes = plt.subplots(len(rows), 2, figsize = (5.91, H), sharex = True, sharey = True, squeeze = False)
    fig.subplots_adjust(left = 0.12, right = 0.99, bottom = 0.45/H, top = 1 - 0.2/H, wspace = 0.1, hspace = 0.16)
    marks = []
    for r, (title, kw) in enumerate(rows):
        fs = {k: T.thetaFit(k, units = units, **kw) for k in (1, 2)}
        part = '_'.join(str(x) for x in kw.values())
        for c, k in enumerate((1, 2)):
            ax, f = axes[r, c], fs[k]
            marks += _fitPanel(ax, f, name, '{}_layer{}'.format(part, k), units, lim, ms = 3.4, msUS = 5)
            ax.set_title(layerTitle(k) + ('   (n = {})'.format(f['n']) if c == 0 else ''), fontsize = 7,
                         color = INK['primary'], loc = 'left')
            if r == len(rows) - 1:
                ax.set_xlabel('Observed ' + unit, color = INK['secondary'], fontsize = 7)
            if c == 0:                              # the row's reading on the left, over the axis label
                ax.set_ylabel(title + '\nPredicted ' + unit, color = INK['secondary'], fontsize = 7)
    out = _save(fig, name, marks)
    note = (r'\textit{Note:} As figure~\ref{fig:OECD26:thetaFit}' + ('' if units == 'theta' else r', in units of the index')
            + ', ' + FITNOTE[kind] + r'. The U.S.\ is open; the solid lines are least-squares fits over the other '
            r'countries.')
    return out + [_note(name, _src() + ('; results/oecd26_theta_<V>/' if kind == 'measure' else ''), note)]

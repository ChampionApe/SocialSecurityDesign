r""" Figure builders for the supplemental appendix: the calibrations across rho, the baseline path of the endogenous design
(escPathPaper, US_ESC_Path), the tax against a pinned design (escTauTheta). House style from figures.py (_panel, SERIES, INK); drawn at the measure (5.91 in).

Each writes results/paper/Figs/<name>.{pdf,png,svg} and <name>.marks.json through figures._save. Every data
point is its own artist, tagged with figures.markId; the text of its mark is the value as the tables print it
(tablesOA's formatters), with tex escapes removed. Lines are drawn untagged beneath the points.
"""
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import config as C
import datasets as D
import tablesOA as T
from tablesUS import COUNTRYNAME
from figures import SERIES, INK, _panel, _save, markId
from tables import BANNER

MARKERS = ('o', 's', '^')
SURFACE = '#fcfcfb'
WIDTH = 5.91                       # \linewidth of the paper, in inches


def _plain(text):
    """ A tex cell as a mark's text: escapes and math shifts removed ('15.06\\%' -> '15.06%'). """
    return text.replace('\\%', '%').replace('$', '')


def _point(ax, x, y, gid, colour, marker = 'o', ms = 3.6, z = 4):
    """ One tagged data mark, a Line2D of its own so that the svg carries its id exactly once. """
    (pt,) = ax.plot([x], [y], linestyle = 'none', marker = marker, markersize = ms, color = colour,
                    markeredgecolor = SURFACE, markeredgewidth = 0.5, zorder = z)
    pt.set_gid(gid)
    return pt


def _mark(gid, label, series, panel, value, text):
    return {'id': gid, 'label': label, 'series': series, 'panel': panel, 'value': float(value),
            'text': text, 'table': None, 'row': None}


def _style(ax, title, ylabel = '', xlabel = None):
    _panel(ax, title, ylabel, titlesize = 8.5, labelsize = 7)
    ax.yaxis.label.set_fontsize(7.5)
    if xlabel:
        ax.set_xlabel(xlabel, color = INK['secondary'], fontsize = 7.5)


def _key(handles, labels, ax, **kw):
    leg = ax.legend(handles, labels, frameon = False, fontsize = 7, labelcolor = INK['secondary'], **kw)
    return leg


# ---------------------------------------------------------------------------------------------------
# The calibrations across rho
# ---------------------------------------------------------------------------------------------------
def _rhoPanel(ax, xs, ys, gids, colour, marker, highlight = None):
    """ A line through (xs, ys), untagged, and one tagged point per value. `highlight`: the x drawn in ink
    (the paper's rho). """
    ax.plot(xs, ys, color = colour, linewidth = 1.2, zorder = 3)
    for x, y, gid in zip(xs, ys, gids):
        ink = highlight is not None and np.isclose(x, highlight)
        _point(ax, x, y, gid, INK['primary'] if ink else colour, marker, ms = 4.2 if ink else 3.4,
               z = 5 if ink else 4)


ARGPANELS = [('β',  r'Discount factor $\beta$',             'Discount factor β',             '',               1.),
             ('ω',  r'Political weight of the old $\omega$', 'Political weight of the old ω', '',               1.),
             ('sr', 'Savings rate, % of GDP',               'Savings rate {year}, % of GDP', '',               100.),
             ('ι',  r'Informal savings ratio $\iota$',       'Informal savings ratio ι',      '',               1.),
             ('η0', r'Productivity $\eta_0$',                'Informal productivity η0',      '',               1.),
             ('X0', r'Taste for leisure $X_0$',              'Informal taste for leisure X0', '',               1.)]
ARGASCII = {'β': 'beta', 'ω': 'omega', 'sr': 'sr', 'ι': 'iota', 'η0': 'eta0', 'X0': 'X0'}    # mark-id parts


def argRhoGrid(commonX = None):
    r""" The Argentine calibration across the rho grid (results/calibration/informalSavings_rhoGrid*): beta,
    omega, the savings rate and the informal savings ratio of the calibration year, eta_0 and X_0, one panel
    each against rho, the paper's rho in ink. The common X is not drawn: the vector-X sweep does not record
    it. """
    commonX = C.ARG['commonX'] if commonX is None else commonX
    sfx = C.variantSuffix(commonX, 'ARG')
    name = 'ARG_RhoGrid' + sfx
    g = T.argRhoGrid(commonX)
    fmt = {c: f for c, _, f in T.ARGRHOCOLS}
    year = C.calendar()['year0']
    ρ1 = C.ARG['ρBaseline']
    fig, axes = plt.subplots(2, 3, figsize = (WIDTH, 3.9), constrained_layout = True, sharex = True)
    marks = []
    for k, (ax, (col, title, plainTitle, ylabel, scale)) in enumerate(zip(axes.ravel(), ARGPANELS)):
        _style(ax, title.format(year = year), ylabel, '$\\rho$' if k >= 3 else None)
        xs, ys = g['ρ'].values.astype(float), scale*g[col].values.astype(float)
        gids = [markId(name, ARGASCII[col], 'rho' + C.num(x, 1)) for x in xs]
        _rhoPanel(ax, xs, ys, gids, SERIES[0], 'o', highlight = ρ1)
        for x, v, gid in zip(xs, g[col].values.astype(float), gids):
            marks.append(_mark(gid, 'ρ = ' + C.num(x, 1), '', plainTitle.format(year = year), scale*v,
                               _plain(fmt[col](v))))
    dot = lambda c, ms: Line2D([], [], linestyle = 'none', marker = 'o', markersize = ms, color = c,
                               markeredgecolor = SURFACE, markeredgewidth = 0.5)
    fig.legend([dot(SERIES[0], 3.4), dot(INK['primary'], 4.2)],
               ['calibration at each $\\rho$', '$\\rho = ' + C.num(ρ1, 0) + '$, the paper\'s specification'],
               loc = 'outside lower center', ncol = 2, frameon = False, fontsize = 7, labelcolor = INK['secondary'])
    return _save(fig, name, marks)


OECDPANELS = [('β',  r'Discount factor $\beta$ (all three)', 'Discount factor β',            '',               1.),
              ('ω',  r'Political weight of the old $\omega$', 'Political weight of the old ω', '',               1.),
              ('X',  r'Taste for leisure $X$',               'Taste for leisure X',           '',               1.),
              ('R',  r'Interest factor $R$, 30 years',       'Interest factor R, 30 years',   '',               1.),
              ('sr', 'Savings rate, % of GDP',               'Savings rate {year}, % of GDP', '',               100.)]
ASCII = {'β': 'beta', 'ω': 'omega', 'X': 'X', 'R': 'R', 'sr': 'sr'}


def oecdRhoGrid(commonX = None):
    r""" The U.S., UK and French calibrations across the rho grid (results/calibration/{US,UK,FR}_rhoGrid*):
    beta (calibrated for the U.S. and imposed on the other two, so drawn once), omega, the common X, the
    30-year interest factor R and the savings rate of the calibration year, one panel each against rho, a
    series per country, the paper's rho in ink. Under vector X the sweeps record no X and its panel is left
    out. """
    commonX = C.US['commonX'] if commonX is None else commonX
    sfx = C.variantSuffix(commonX)
    name = 'OECD_RhoGrid' + sfx
    sw = T.oecdSweeps(commonX)
    fmt = {c: f for c, _, f in T.OECDRHOCOLS}
    year = C.usCalendar()['year0']
    ρ1 = C.US['ρBaseline']
    panels = [p for p in OECDPANELS if commonX or p[0] != 'X']
    fig, axes = plt.subplots(2, 3, figsize = (WIDTH, 4.1), constrained_layout = True, sharex = True)
    slots = axes.ravel()
    marks = []
    colours = dict(zip([c for c, _ in T.OECDHOSTS], (SERIES[0], SERIES[1], INK['secondary'])))
    markers = dict(zip([c for c, _ in T.OECDHOSTS], MARKERS))
    for k, (ax, (col, title, plainTitle, ylabel, scale)) in enumerate(zip(slots, panels)):
        # an axis gets the rho label where no panel sits below it
        below = k >= 3 or k + 3 >= len(panels)
        _style(ax, title.format(year = year), ylabel, '$\\rho$' if below else None)
        ax.tick_params(labelbottom = below)
        hosts = T.OECDHOSTS[:1] if col == 'β' else T.OECDHOSTS
        for c, key in hosts:
            g = sw[c]
            xs, vs = g['ρ'].values.astype(float), g[col].values.astype(float)
            gids = [markId(name, ASCII[col], key, 'rho' + C.num(x, 1)) for x in xs]
            _rhoPanel(ax, xs, scale*vs, gids, colours[c], markers[c], highlight = ρ1)
            for x, v, gid in zip(xs, vs, gids):
                marks.append(_mark(gid, 'ρ = ' + C.num(x, 1), COUNTRYNAME[c], plainTitle.format(year = year),
                                   scale*v, _plain(fmt[col](v))))
    for ax in slots[len(panels):]:
        ax.axis('off')
    handles = [Line2D([], [], color = colours[c], marker = markers[c], markersize = 3.4, linewidth = 1.2,
                      markeredgecolor = SURFACE, markeredgewidth = 0.5) for c, _ in T.OECDHOSTS]
    handles.append(Line2D([], [], linestyle = 'none', marker = 'o', markersize = 4.2, color = INK['primary'],
                          markeredgecolor = SURFACE, markeredgewidth = 0.5))
    _key(handles, [COUNTRYNAME[c] for c, _ in T.OECDHOSTS] + ['$\\rho = ' + C.num(ρ1, 0) + '$, the paper\'s'],
         slots[len(panels)], loc = 'center')
    return _save(fig, name, marks)


# ---------------------------------------------------------------------------------------------------
# The 2020 tax against inherited savings
# ---------------------------------------------------------------------------------------------------
def usTaxState():
    r""" US_TaxState: the U.S. tax policy of 2020 against the savings the 2020 economy inherits, as a multiple of
    their equilibrium value, one line per rho of the paper (datasets.taxPolicyState; design given, each rho at
    its own calibration). The equilibrium point, where every line passes through the calibrated tax, is the
    vertical rule. Colours as the supplemental appendix's rho series. """
    name = 'US_TaxState'
    df = D.taxPolicyState()
    year0 = C.usCalendar()['year0']
    colours = {0.5: SERIES[0], 1.0: SERIES[1], 2.0: INK['secondary']}
    fig, ax = plt.subplots(1, 1, figsize = (WIDTH*0.62, 2.5), constrained_layout = True)
    _style(ax, '', 'Tax rate in {}, percent'.format(year0),
           'Savings inherited in {}, relative to the equilibrium path'.format(year0))
    ax.axvline(1., color = INK['secondary'], linewidth = 0.6, linestyle = (0, (2, 2)), zorder = 1)
    marks, handles, labels = [], [], []
    for k, ρ in enumerate(sorted(df['ρ'].unique())):
        if float(ρ) not in colours:
            raise ValueError('US_TaxState: no colour for ρ = {}'.format(ρ))
        sub = df[np.isclose(df['ρ'], ρ) & np.isfinite(df['τ2020'])]
        xs, ys = sub['s_2019/s*'].values.astype(float), 100*sub['τ2020'].values.astype(float)
        colour = colours[float(ρ)]
        ax.plot(xs, ys, color = colour, linewidth = 1.2, zorder = 3)
        for x, y in zip(xs, ys):
            gid = markId(name, 'rho' + C.num(ρ, 1), 's' + C.num(x, 2))
            _point(ax, x, y, gid, colour, MARKERS[k], ms = 3.2)
            marks.append(_mark(gid, 'inherited savings × ' + C.num(x, 2), 'ρ = ' + C.num(ρ, 1),
                               'Tax rate in {}'.format(year0), y, C.num(y, 2) + '%'))
        handles.append(Line2D([], [], color = colour, marker = MARKERS[k], markersize = 3.2, linewidth = 1.2,
                              markeredgecolor = SURFACE, markeredgewidth = 0.5))
        labels.append(r'$\rho = ' + C.num(ρ, 1) + '$' + (' (log)' if np.isclose(ρ, 1.0) else ''))
    _key(handles, labels, ax, loc = 'upper right')
    return _save(fig, name, marks)


# ---------------------------------------------------------------------------------------------------
# The endogenous design's baseline path
# ---------------------------------------------------------------------------------------------------
def escPathPaper():
    r""" The paper's figure of the baseline path (sec:esc, fig:US_ESC:path), at the baseline rho (log
    preferences) and from the calibration year on, where the calibration pins the design and everything after
    is the model's: workers per retiree, the design chosen against theta* pinned, and the tax rate with the
    design chosen and with it pinned (datasets.escPath). The path from 1990, the other rhos and the Frisch
    elasticities are the supplemental appendix's (siteData.escPaths). """
    name = 'US_ESC_Path'
    ρ1, spec = C.US['ρBaseline'], C.US['esc']['spec']
    year0 = C.usCalendar()['year0']
    df = D.escPath(ρ1)
    df = df[df['date'].astype(int) >= year0].reset_index(drop = True)
    θstar = float(D.escCalibration()[(ρ1, spec)]['θStar'])
    xs = df['date'].values.astype(int)
    chosen, pinned = SERIES[0], SERIES[1]
    dashed = (0, (4, 2))
    fig, axes = plt.subplots(1, 3, figsize = (WIDTH, 2.3), constrained_layout = True)
    panels = [('nu', r'Workers per retiree $\nu_t$', 'Workers per retiree ν_t', ''),
              ('theta', r'Pension design $\theta_t$', 'Pension design θ_t', ''),
              ('tau', r'Tax rate $\tau_t$', 'Tax rate τ_t', 'percent')]
    marks = []
    for ax, (pid, title, plainTitle, ylabel) in zip(axes, panels):
        _style(ax, title, ylabel, 'year')
        # a tick per period, labelled every other one: six four-digit years do not fit a third of the measure
        ax.set_xticks(xs)
        ax.set_xticklabels([str(x) if k % 2 == 0 else '' for k, x in enumerate(xs)])
        if pid == 'nu':
            series = [('ν', INK['secondary'], '', 'o', 'solid', 1., lambda r, v: C.num(v, 2))]
        elif pid == 'theta':
            ax.axhline(θstar, color = pinned, linewidth = 1.2, linestyle = dashed, zorder = 2)
            series = [('θ', chosen, 'chosen', 'o', 'solid', 1., lambda r, v: _plain(T.escPathCells(r)[0]))]
        else:
            series = [('τExo', pinned, 'pinned', 's', dashed, 100., lambda r, v: _plain(T.escPathCells(r)[2])),
                      ('τ', chosen, 'chosen', 'o', 'solid', 100., lambda r, v: _plain(T.escPathCells(r)[1]))]
        for col, colour, label, marker, ls, scale, text in series:
            vs = scale*df[col].values.astype(float)
            ax.plot(xs, vs, color = colour, linewidth = 1.2, linestyle = ls, zorder = 3)
            for (_, r), x, v in zip(df.iterrows(), xs, vs):
                gid = markId(name, pid, label or 'data', x)
                _point(ax, x, v, gid, colour, marker)
                marks.append(_mark(gid, str(x), label, plainTitle, v, text(r, v / scale)))
    # room under the pinned line, and its label above it at the right, where the chosen path has left it
    lo, hi = axes[1].get_ylim()
    axes[1].set_ylim(min(lo, θstar - 0.012), hi)
    axes[1].annotate(r'pinned, $\theta^{\ast} = ' + C.num(θstar, 3) + '$', xy = (1., θstar),
                     xycoords = ('axes fraction', 'data'), xytext = (-2, 2), textcoords = 'offset points',
                     ha = 'right', va = 'bottom', fontsize = 7, color = INK['secondary'])
    handles = [Line2D([], [], color = chosen, marker = 'o', markersize = 3.6, linewidth = 1.2,
                      markeredgecolor = SURFACE, markeredgewidth = 0.5),
               Line2D([], [], color = pinned, marker = 's', markersize = 3.6, linewidth = 1.2, linestyle = dashed,
                      markeredgecolor = SURFACE, markeredgewidth = 0.5)]
    leg = axes[2].legend(handles, [r'$\theta$ chosen', r'$\theta$ pinned at $\theta^{\ast}$'], loc = 'lower right',
                         fontsize = 7, labelcolor = INK['secondary'], frameon = False)
    leg.set_zorder(6)
    return _save(fig, name, marks)


# ---------------------------------------------------------------------------------------------------
# The tax against a pinned design
# ---------------------------------------------------------------------------------------------------
def tauThetaPeaks(d):
    """ For one frame of datasets.escTauTheta: (θ at the largest τ, that τ, τ at θ = 0, τ at θ = 1, whether τ falls
    at every step), τ in percent. """
    τ = 100*d['τ'].values.astype(float)
    k = int(np.argmax(τ))
    return float(d['θpin'].iloc[k]), float(τ[k]), float(τ[0]), float(τ[-1]), bool(np.all(np.diff(τ) < 0))


def escTauTheta():
    r""" US_ESC_TauTheta: the U.S. tax rate of 2020 against a design pinned on the whole path (datasets.escTauTheta,
    python/US/runTauTheta.py), without the cost (section 6's calibration) and under the paper's cost at its calibrated
    λ, each economy calibrated at the observed θ* (vertical rule); every point is a pinned design (the grid, θ*, and
    the steps around the costly economy's peak). Also writes `US_ESC_TauTheta_note.tex`: the peak and the corners of
    both curves. """
    name = 'US_ESC_TauTheta'
    spec = C.US['esc']['spec']
    tt = D.escTauTheta()
    year0 = C.usCalendar()['year0']
    θstar = float(C.theta0())
    fig, ax = plt.subplots(1, 1, figsize = (WIDTH*0.62, 2.6), constrained_layout = True)
    _style(ax, '', 'Tax rate in {}, percent'.format(year0), r'Pension design $\theta$, pinned at every date')
    ax.axvline(θstar, color = INK['secondary'], linewidth = 0.6, linestyle = (0, (2, 2)), zorder = 1)
    marks, handles, labels = [], [], []
    series = (('none', SERIES[0], 'o', 'without the cost'), (spec, SERIES[1], 's', 'with the cost'))
    for which, colour, marker, label in series:
        d = tt[which]
        xs, ys = d['θpin'].values.astype(float), 100*d['τ'].values.astype(float)
        ax.plot(xs, ys, color = colour, linewidth = 1.2, zorder = 3)
        for x, y, kind in zip(xs, ys, d['kind']):
            gid = markId(name, which, 'theta' + '{:.4f}'.format(x))
            _point(ax, x, y, gid, colour, marker, ms = 3.2 if kind == 'grid' else 2.0)
            marks.append(_mark(gid, 'θ = {:.3f}'.format(x), label, 'Tax rate in {}'.format(year0), y,
                               C.num(y, 2) + '%'))
        handles.append(Line2D([], [], color = colour, marker = marker, markersize = 3.2, linewidth = 1.2,
                              markeredgecolor = SURFACE, markeredgewidth = 0.5))
        labels.append(label)
    handles.append(Line2D([], [], color = INK['secondary'], linewidth = 0.6, linestyle = (0, (2, 2))))
    labels.append(r'observed $\theta^{\ast}$')
    _key(handles, labels, ax, loc = 'upper right')
    out = _save(fig, name, marks)
    θ0n, peak0, lo0, hi0, falls0 = tauThetaPeaks(tt['none'])
    θ1, peak1, lo1, hi1, falls1 = tauThetaPeaks(tt[spec])
    lam = float(tt[spec]['p'].iloc[0])
    pc = lambda v: C.num(v, 2) + r'\%'
    first = ((r'Without the cost the tax rate falls with the design, from ' + pc(lo0) + r' at $\theta = 0$ to '
              + pc(hi0) + r' at $\theta = 1$ (' + C.pp((hi0 - lo0)/100) + ').') if falls0 else
             (r'Without the cost the tax rate is ' + pc(lo0) + r' at $\theta = 0$, peaks at ' + pc(peak0)
              + r' at $\theta = ' + C.num(θ0n, 2) + '$ and is ' + pc(hi0) + r' at $\theta = 1$.'))
    second = (r' With the cost it is ' + pc(lo1) + r' at $\theta = 0$, peaks at ' + pc(peak1) + r' at $\theta = '
              + C.num(θ1, 2) + '$, against ' + pc(100*float(tt[spec].loc[np.isclose(tt[spec]['θpin'], θstar), 'τ'].iloc[0]))
              + r' at $\theta^{\ast}$, and falls to ' + pc(hi1) + r' at $\theta = 1$.')
    step = tt[spec].loc[tt[spec]['kind'] == 'refine', 'θpin']
    grid = C.US['tauTheta']['grid']                              # the note names it as the first, second and last point
    if not np.allclose(np.diff(grid), grid[1] - grid[0]):
        raise ValueError("config.US['tauTheta']['grid'] is not equally spaced: {}".format(grid))
    note = (r'\textit{Note:} The U.S.\ in ' + str(year0) + r' under log preferences ($\rho = 1$), with the design '
            r'pinned at the same $\theta$ at every date and the economy started from its own steady state; $\beta$, '
            r'$\omega$ and $X$ stay at their calibration at the observed $\theta^{\ast} = ' + C.num(θstar, 3)
            + r'$ (dashed line), without the cost as in section \ref{sec:oecd} and with the paper\textquotesingle s '
            r'cost at $\lambda = ' + C.num(lam, 2) + r'$ (table \ref{table:US_ESC:calibration}). ' + first + second
            + (r' The designs are $\theta = {:g}, {:g}, \dots, {:g}$, $\theta^{{\ast}}$ and, around the peak, steps of '.format(
                   grid[0], grid[1], grid[-1])
               + C.num(float(np.diff(np.sort(step.values)).min()), 2) + '.' if len(step) > 1 else ''))
    path = os.path.join(C.PAPERDIR, 'Figs', name + '_note.tex')
    with open(path, 'w', encoding = 'utf-8', newline = '\n') as f:
        f.write(BANNER.format(name = name, src = 'results/esc/escTauTheta.csv, escCalibration.csv') + note + '\n')
    return out + [path]

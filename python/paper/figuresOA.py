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


# ---------------------------------------------------------------------------------------------------
# Page D1's section on inequality: the Gini against the model's measures
# ---------------------------------------------------------------------------------------------------
# Label offsets in points, with a cost (pixels^2) that keeps the far ones for when every near one collides.
OFFSETS = [(3, 2, 'left', 'bottom', 0), (3, -2, 'left', 'top', 0), (-3, 2, 'right', 'bottom', 0),
           (-3, -2, 'right', 'top', 0), (0, 4, 'center', 'bottom', 0), (0, -4, 'center', 'top', 0),
           (5, 0, 'left', 'center', 0), (-5, 0, 'right', 'center', 0), (7, 5, 'left', 'bottom', 25),
           (-7, -5, 'right', 'top', 25), (7, -5, 'left', 'top', 25), (-7, 5, 'right', 'bottom', 25),
           (0, 8, 'center', 'bottom', 25), (0, -8, 'center', 'top', 25)]


def _isoLabels(ax, xs, ys, texts, bold = (), avoid = (), fontsize = 5.6, ambiguous = 40., passes = 3):
    """ Point labels at the first of OFFSETS with no overlap, else the least bad: the score is the offset's cost plus
    the overlap (pixels^2) with placed labels (x4), other markers and the boxes `avoid`, plus `ambiguous` (x10 at a far
    offset) where another marker is nearer the label than its own, plus 1e5 outside the axes. Crowded points first,
    then `passes` passes that re-place each label against all the others. Needs the final layout. """
    from matplotlib.transforms import Bbox
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.c_[xs, ys])
    rad = 2.2*fig.dpi/72
    dots = [Bbox.from_extents(px - rad, py - rad, px + rad, py + rad) for px, py in pts]
    box = ax.get_window_extent(rend)
    ov = lambda a, b: max(0., min(a.x1, b.x1) - max(a.x0, b.x0))*max(0., min(a.y1, b.y1) - max(a.y0, b.y0))
    dist = lambda b, p: np.hypot(max(b.x0 - p[0], 0, p[0] - b.x1), max(b.y0 - p[1], 0, p[1] - b.y1))

    def score(i, bb, others, far):
        s = 4*sum(ov(bb, q) for q in others) + sum(ov(bb, d) for j, d in enumerate(dots) if j != i)
        s += sum(ov(bb, a) for a in avoid)
        own = dist(bb, pts[i])
        s += ambiguous*(10 if far else 1)*any(dist(bb, pts[j]) < own for j in range(len(pts)) if j != i)
        if not (box.x0 <= bb.x0 and bb.x1 <= box.x1 and box.y0 <= bb.y0 and bb.y1 <= box.y1):
            s += 1e5
        return s

    def place(i, others):
        best = None
        for dx, dy, ha, va, cost in OFFSETS:
            t = ax.annotate(texts[i], (xs[i], ys[i]), xytext = (dx, dy), textcoords = 'offset points', ha = ha,
                            va = va, fontsize = fontsize, zorder = 6,
                            color = INK['primary'] if texts[i] in bold else INK['secondary'],
                            fontweight = 'bold' if texts[i] in bold else 'normal')
            bb = t.get_window_extent(rend).expanded(1.04, 1.08)
            s = cost + score(i, bb, others, cost > 0)
            if best is None or s < best[0]:
                if best is not None:
                    best[1].remove()
                best = (s, t, bb)
            else:
                t.remove()
            if s == 0:
                break
        return best

    crowd = [-sum(ov(dots[i].expanded(4, 4), d) > 0 for d in dots) for i in range(len(xs))]
    placed = {}
    for i in np.argsort(crowd, kind = 'stable'):
        placed[i] = place(i, [b for _, _, b in placed.values()])
    for _ in range(passes):
        for i in list(placed):
            placed[i][1].remove()
            placed[i] = place(i, [b for j, (_, _, b) in placed.items() if j != i])


def oecdInequalityMeasures():
    r""" OECD_InequalityMeasures (page D1): on the countries of figure 1 (tablesOA.inequalityMeasures), the Gini of
    disposable income against the two measures of the income groups the model is fed (tablesOA.INEQUALITYMEASURES:
    y_3/y_1 and Ṽ), one panel each, with the least-squares line over all of them and over all but the U.S., and Pearson's
    r and Spearman's rank correlation with and without the U.S. Writes results/paper/oecdInequalityCorrelations.csv
    (tablesOA.inequalityCorrelations, which OECD_InequalityCorrelations prints) and its note. Every dot and line is a
    mark, a dot linked to its country's row of OECD_Countries, a line and the statistics to their row of
    OECD_InequalityCorrelations. """
    import oecdFigure1 as F1
    from figuresOECD26 import EDGE, yearsText, INEQUALITYSOURCE
    name = 'OECD_InequalityMeasures'
    m = T.inequalityMeasures()
    corr = T.inequalityCorrelations(m)
    os.makedirs(C.PAPERDIR, exist_ok = True)
    csv = D.paperCsv(D.OECDINEQCORR)
    corr.to_csv(csv, index = False)
    cols = T.INEQUALITYMEASURES
    H = 2.75*len(cols) + 0.7
    fig, axes = plt.subplots(len(cols), 1, figsize = (WIDTH, H), sharex = True, squeeze = False)
    # a fixed layout: constrained_layout counts the source lines and collapses the axes
    fig.subplots_adjust(left = 0.09, right = 0.985, top = 1 - 0.35/H, bottom = 0.75/H, hspace = 0.42)
    x = m['gini'].astype(float)
    lo, hi = float(x.min()), float(x.max())
    pad = 0.06*(hi - lo)
    G = C.INEQUALITYLABELS['gini']
    marks, statBox = [], {}
    for k, (ax, col) in enumerate(zip(axes[:, 0], cols)):
        lab = C.INEQUALITYLABELS[col]
        panel = '{} against {}'.format(lab['plain'], G['plain'])
        c = corr[(corr['measure'] == col) & (corr['sample'] == 'all')].iloc[0]
        cn = corr[(corr['measure'] == col) & (corr['sample'] == 'noUSA')].iloc[0]
        _panel(ax, '', '', labelsize = 6.5)
        ax.set_ylabel(lab['axis'].split(', ')[-1], color = INK['secondary'], fontsize = 7.5)
        # the title's second line is blank: the source line is drawn into it in secondary ink
        ax.set_title('({}) {}\n'.format('abc'[k], lab['axis']), color = INK['primary'], fontsize = 8, loc = 'left',
                     pad = 3)
        ax.text(0, 1.0, INEQUALITYSOURCE[col] + ', ' + yearsText(m['year']), transform = ax.transAxes, ha = 'left',
                va = 'bottom', fontsize = 6.3, color = INK['secondary'])
        ax.tick_params(labelbottom = True, length = 2.5)
        y = m[col].astype(float)
        for iso in m.index:
            pt, = ax.plot(x[iso], y[iso], 'o', ms = 3.8, ls = 'none', zorder = 4 if iso == 'USA' else 3,
                          color = SERIES[0], mec = EDGE, mew = 0.5)
            pt.set_gid(markId(name, col, iso.lower()))
            marks.append({'id': pt.get_gid(), 'label': m.at[iso, 'country'], 'series': 'data', 'panel': panel,
                          'value': float(y[iso]),
                          'text': '{}: {}; {}: {}'.format(G['mark'], G['fmt'].format(x[iso]), lab['mark'],
                                                          lab['fmt'].format(y[iso])),
                          'table': 'OECD_Countries', 'row': iso.lower()})
        nu = np.asarray(m.index != 'USA')
        for o, keep, style, part in ((c, np.ones(len(m), bool), {'color': SERIES[0], 'lw': 1.1}, 'all'),
                                     (cn, nu, {'color': INK['muted'], 'lw': 0.9, 'ls': (0, (3, 2))}, 'nousa')):
            xs = np.array([float(x[keep].min()), float(x[keep].max())])
            ln, = ax.plot(xs, o['intercept'] + o['slope']*xs, zorder = 2, **style)
            ln.set_gid(markId(name, col, 'fit', part))
            marks.append({'id': ln.get_gid(), 'series': 'data', 'panel': panel, 'value': float(o['slope']),
                          'label': 'least squares, ' + ('all' if part == 'all' else 'without the U.S.'),
                          'text': 'slope {} (s.e. {})'.format(lab['fmt'].format(o['slope']), lab['fmt'].format(o['se'])),
                          'table': 'OECD_InequalityCorrelations', 'row': T._keyOf(col) + '_' + part})
        # headroom: the statistics box spans the top left of the axes, so every point under it (left 70% of the
        # axis) stays in the lower 70% of the axis, room for a label above it
        ymin, ymax = float(y.min()), float(y.max())
        bottom = ymin - 0.08*(ymax - ymin)
        yLeft = float(y[x <= lo - pad + 0.70*(hi - lo + 2*pad)].max())
        ax.set_ylim(bottom, max(ymax + 0.10*(ymax - ymin), (yLeft - 0.30*bottom)/0.70))
        ax.set_xlim(lo - pad, hi + pad)
        statBox[ax] = ax.text(0.03, 0.97, 'r = {:.2f}   (without the U.S. {:.2f})\nSpearman {:.2f}   ({:.2f}),  n = {}'.format(
            c['pearson'], cn['pearson'], c['spearman'], cn['spearman'], int(c['n'])), transform = ax.transAxes,
            ha = 'left', va = 'top', fontsize = 6.2, color = INK['secondary'], linespacing = 1.35, zorder = 7,
            bbox = {'facecolor': EDGE, 'edgecolor': 'none', 'alpha': 0.9, 'pad': 1.5})
        statBox[ax].set_gid(markId(name, col, 'stats'))
        marks.append({'id': statBox[ax].get_gid(), 'label': 'Pearson r', 'series': 'data', 'panel': panel,
                      'value': float(c['pearson']),
                      'text': 'r {:.2f}, without the U.S. {:.2f}; Spearman {:.2f}, without the U.S. {:.2f}'.format(
                          c['pearson'], cn['pearson'], c['spearman'], cn['spearman']),
                      'table': 'OECD_InequalityCorrelations', 'row': T._keyOf(col) + '_all'})
    axes[-1, 0].set_xlabel('{} (WIID, {})'.format(F1.LABEL['gini'], yearsText(m['gini_year'])),
                           color = INK['secondary'], fontsize = 7.5)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for ax, col in zip(axes[:, 0], cols):
        _isoLabels(ax, x.to_numpy(float), m[col].to_numpy(float), list(m.index), bold = ('USA',),
                   avoid = [statBox[ax].get_window_extent(rend).expanded(1.06, 1.25)])
    h = [Line2D([], [], color = SERIES[0], lw = 1.1), Line2D([], [], color = INK['muted'], lw = 0.9, ls = (0, (3, 2)))]
    fig.legend(h, ['least squares, all {}'.format(len(m)), 'least squares, without the U.S.'],
               loc = 'lower center', ncol = 2, frameon = False, fontsize = 6.5, labelcolor = INK['secondary'],
               handlelength = 1.8, columnspacing = 1.4)
    out = _save(fig, name, marks)
    g = D.incomeGroups(m.index)
    shares = ', '.join(r'{:.0f}\%'.format(100*float(g['gamma{}'.format(k)].iloc[0])) for k in (1, 2, 3))
    note = (r'\textit{Note:} The ' + str(len(m)) + r' countries of figure \ref{fig:US:OECDdata}. The Gini is that of '
            r'disposable income (WIID); high over low income, $y_3/y_1$, and the dispersion $\tilde V = \sum_i \gamma_i '
            r'(y_i - 1)^2/y_i$ are of the three income groups the model is fed, at the U.S.\ shares of the population ('
            + shares + r', from the bottom), in the gross earnings of men working full time (OECD). The solid line is '
            r'the least-squares fit over all countries, the dashed one over all but the U.S.; each panel gives '
            r'Pearson\textquotesingle s $r$ and Spearman\textquotesingle s rank correlation, with the U.S.\ and, in '
            r'parentheses, without it (table \ref{table:OA:OECD_InequalityCorrelations}).')
    path = os.path.join(C.PAPERDIR, 'Figs', name + '_note.tex')
    with open(path, 'w', encoding = 'utf-8', newline = '\n') as f:
        f.write(BANNER.format(name = name, src = 'data/oecdFigure1.csv, oecd29_distribution.csv, oecd26_theta.csv')
                + note + '\n')
    return out + [csv, path]

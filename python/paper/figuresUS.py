r""" Figure builders for the US/France/UK arm. figures.py is the Argentina arm and owns the house style;
this module imports it rather than restating it, so the two arms stay one visual family.

Style rules that apply here specifically:
  * The counterfactuals are a CATEGORICAL set, not a continuous parameter, so they get the categorical
    pair -- never THETA_RAMP, which is reserved for a parameter that varies continuously.
  * The baseline is a reference line, not a series: it is drawn in ink, not in a series colour, so it
    cannot be mistaken for one of the scenarios being compared against it.
  * Every panel plots the DEVIATION from the baseline, never the level. The levels differ by a few
    points on a 14% base and a 39-hour base, which a level axis compresses into indistinguishable bars;
    the deviation is also the only thing the workweek column means, since under vector X the level of
    hbar is not identified and only its ratio to the baseline is a result. The one exception is the
    pension design in the ESC figure, where the corners 0 and 1 are meaningful and the level is read.
  * Savings are reported over GDP, s/Y, in both figures -- the shock csv's `srOverY`, and in the ESC
    figure (1 - alpha) times the csv's s/(wh), since wh = (1 - alpha) Y. alpha is read from the
    calibration summary, never typed.
  * NOT config.pct anywhere in here. It escapes the percent sign for tex and matplotlib renders the
    backslash literally, so a legend built with it reads "14.4\%". Anything drawn INTO a figure needs
    plain formatting; config.pct is for tex cells only.
"""
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

import config as C
import datasets as D
from figures import SERIES, INK, _panel, _save, markId, mark, signed, plain
from tablesUS import ROWKEY


# The rows of the overview, in the order the paper's discussion moves in: pension design, then ageing (the
# two it identifies as the main determinants), then France's characteristics.
# Reading the finished figure top-to-bottom should reproduce that ranking, so the order is fixed here
# rather than sorted by effect size: a figure whose ordering changes with the data cannot be referred to
# in prose. Each row is (row label, [(csv scenario, mark label), ...], how a second scenario is drawn). A row with two scenarios is one bar per rho (_barPanel): 'tone' draws both from the baseline,
# the FIRST in the series colour and the second in its light tone (theta = 0 and theta = 1, on either side
# of the US design); 'line' draws the first as the bar and marks the second by a line across it (mild
# ageing inside acute ageing). The figure carries no key for either: the paper's note explains them. The
# composite rows (all three French characteristics, France's own path) are not drawn.
OVERVIEW = [('Pension design',       [(r'$\theta = 0$', r'$\theta = 0$'), (r'$\theta = 1$', r'$\theta = 1$')],
             'tone'),
            ('Ageing',               [('Acute ageing', 'Acute ageing'), ('Mild ageing', 'Mild ageing')], 'line'),
            ('French income distr.', [('Income distribution', 'French income distr.')], None),
            ('French voting',        [('Voting', 'French voting')], None),
            ('French leisure',       [('Leisure preferences', 'French leisure')], None)]

# The light tone of a 'tone' bar: the series colour at this share over the surface, which keeps the
# lightest of the three (orange) at about 2:1 against the surface, the floor at which it still reads as a
# mark in print. Opaque rather than an alpha, so the shorter tone drawn over the longer keeps its colour.
TINT, SURFACE = 0.6, '#ffffff'
# The width of a 'line' mark, in points. It is drawn in the surface colour where it falls inside its bar,
# and in the series colour where it does not, so it never vanishes against the surface.
MARKLW = 2.0
# A row every one of whose bars rounds to zero at the paper's precision (one decimal, style guide §3) says
# so in words: an empty bar group reads as a missing run. Only the leisure rows qualify: in the tax and
# savings panels of section 6's figure, in the design and tax panels of section 7's (THETAZERO, half the
# design's last printed digit, for its level).
ZEROROW = 0.05
THETAZERO = 5e-4

# The three reported quantities, as (csv column, panel title, axis label, how to scale a deviation).
# tau and srOverY are fractions on the csv and are read in percentage points; the workweek is already
# in hours (normalised inside the experiment script against that rho's own baseline) and its deviation
# is in hours. Do NOT re-derive the workweek from hbar here -- see tablesUS.py.
PANELS = [('τ',        'Equilibrium tax rate', 'Percentage points',        100.),
          ('srOverY',  'Savings over GDP',     'Percentage points of GDP', 100.),
          ('workweek', 'Average workweek',     'Hours per week',           1.)]

# Colours for the three rho values: the categorical pair plus one ink shade -- three is one more than
# the validated pair carries, so the extra slot is deliberately NEUTRAL rather than a third hue guessed
# by eye. Fixed order, never cycled.
RHOCOLOURS = [SERIES[0], SERIES[1], INK['secondary']]

# A mark's id part and the unit its text carries, per plotted column (figures.signed). None: a level of
# theta, printed to three decimals.
MARKUNIT = {'τ': ('tau', 'p.p.'), 'srOverY': ('sr', 'p.p. of GDP'), 'workweek': ('ww', 'hours'),
            'θ_t0': ('theta', None), 'τ_t0': ('tau', 'p.p.'), 'sr_t0': ('sr', 'p.p. of GDP'),
            'ww_t0': ('ww', 'hours')}


def _rhoPart(ρ):
    return 'rho{:.1f}'.format(ρ)


def _markText(v, unit):
    return '{:.3f}'.format(v) if unit is None else signed(v, unit = unit)


def _tint(colour):
    """ The light tone of `colour`: TINT of it over the surface, as an opaque rgb. """
    return tuple(TINT*np.array(mcolors.to_rgb(colour)) + (1 - TINT)*np.array(mcolors.to_rgb(SURFACE)))


def _barPanel(ax, title, xlabel, labels, series, colours, gids = None, kinds = None):
    r""" One horizontal grouped-bar panel: `series` is a list of (label, values) over `labels`, drawn as
    one bar per series within each scenario group. `gids[k][j]` tags series k's bar for label j.

    `kinds[j]` (OVERVIEW's last field) makes row j's values pairs, with `gids[k][j]` a pair too. 'tone':
    two scenarios in one bar, both measured from the baseline, the first in the series colour and the
    second in its light tone (_tint); the longer is drawn first, so where the two lie on the same side of
    zero the shorter stays visible over it, parted by a surface gap where it ends, and on opposite sides
    the zero line parts them. 'line': the first is the bar and the second a line across it (MARKLW). A row
    whose bars all round to zero (ZEROROW) is labelled 'unchanged'. """
    y = np.arange(len(labels))
    height = 0.8/len(series)
    kinds = kinds or [None]*len(labels)
    # _panel FIRST: it sets tick_params, which would otherwise recolour the scenario labels to the muted
    # ink meant for numeric ticks. These are the figure's row headings and belong in primary ink.
    _panel(ax, title, '', titlesize = 12, labelsize = 10)
    for k, (lab, vals) in enumerate(series):
        yk = y + (k - (len(series)-1)/2)*height
        for j, (v, kind) in enumerate(zip(vals, kinds)):
            g = gids[k][j] if gids else None
            if kind == 'tone':
                tones = list(zip(v, (colours[k], _tint(colours[k])), g or (None, None)))
            elif kind == 'line':
                tones = [(v[0], colours[k], g[0] if g else None)]
            else:
                tones = [(v, colours[k], g)]
            for x, c, gid in sorted(tones, key = lambda t: -abs(t[0])):
                # One legend entry per series, from its first row's full tone.
                key = lab if j == 0 and c is colours[k] else '_nolegend_'
                patch = ax.barh(yk[j], x, height = height, color = c, edgecolor = 'none', zorder = 3,
                                label = key).patches[0]
                if gid:
                    patch.set_gid(gid)
            if kind == 'tone' and v[0]*v[1] > 0:
                inner = np.sign(v[0])*min(abs(v[0]), abs(v[1]))
                ax.vlines(inner, yk[j] - height/2, yk[j] + height/2, colors = SURFACE, linewidth = 1.2,
                          zorder = 3.5)
            elif kind == 'line':
                inside = v[0]*v[1] > 0 and abs(v[1]) < abs(v[0])
                line, = ax.plot([v[1], v[1]], [yk[j] - height/2, yk[j] + height/2], linewidth = MARKLW,
                                color = SURFACE if inside else colours[k], solid_capstyle = 'butt',
                                zorder = 3.6)
                if g:
                    line.set_gid(g[1])
    for j in range(len(labels)):
        row = [x for _, vals in series for x in (vals[j] if isinstance(vals[j], tuple) else (vals[j],))]
        if all(abs(x) < ZEROROW for x in row):
            ax.annotate('unchanged', xy = (0, y[j]), xytext = (5, 0), textcoords = 'offset points',
                        ha = 'left', va = 'center', fontsize = 10, style = 'italic', color = INK['muted'])
    ax.axvline(0, color = INK['primary'], linewidth = 1.0, zorder = 4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize = 11, color = INK['primary'])
    ax.set_xlabel(xlabel, color = INK['secondary'], fontsize = 10)
    ax.grid(axis = 'y', visible = False)


def _topDown(axes):
    """ Put the first scenario at the TOP, once for the whole row of panels.

    NOT inside the per-panel helpers: these axes are created with sharey=True, so they are one y-axis
    wearing several faces and every invert_yaxis() call flips it again. Called per panel it reverses
    the scenario order on any EVEN number of panels and leaves it alone on any odd number. Invert once,
    on one axis, and the count stops mattering. """
    axes[0].invert_yaxis()


def _figLegend(fig, handles, labels, bottom = 0.075, ncol = None):
    r""" One legend for the whole figure, below the panels. Figure-level rather than per-panel: the
    series mean the same thing in every panel, so repeating the key would be redundant.

    The baseline the deviations are measured from is NOT drawn here: it is a \tablenotes under the
    \includegraphics in the paper, which keeps its levels out of a second place they would have to be
    kept in step (they are printed in the tables) and leaves the panels the height. """
    # ncol: a five-entry key does not fit one row at this width and is clipped at both ends rather than
    # shrunk, so the caller that has one wraps it.
    fig.tight_layout(rect = (0, bottom, 1, 1))
    fig.legend(handles, labels, loc = 'lower center', ncol = ncol or len(handles), frameon = False,
               fontsize = 11, bbox_to_anchor = (0.5, 0.01))


def _rowsPresent(df, ρs, rows, effect = 'full', drop = False):
    r""" `rows` that have, at every rho, each of the row's scenarios (entry[1]: one csv scenario, or OVERVIEW's list of
    (scenario, label)). A scenario missing at any rho raises MissingInput, the output being skipped: a bar group with a
    hole in it reads as a zero. `drop`: leave such a row out of the figure instead (a twin calibration's). """
    out = []
    for entry in rows:
        names = [entry[1]] if isinstance(entry[1], str) else [s for s, _ in entry[1]]
        try:
            for ρ in ρs:
                for s in names:
                    D.usShockRow(df, ρ, s, effect)
        except D.MissingInput:
            if not drop:
                raise
            continue
        out.append(entry)
    return out


def usOverview(commonX = None):
    r""" Figure \ref{fig:US:overview}: every single-characteristic counterfactual's effect on the tax
    rate, savings over GDP and the average workweek, in deviations from the baseline, at each rho in
    config.US['ρTable'].

    Read as: which characteristics move each outcome, and does that ranking survive the IES. Ageing and
    pension design should dominate the tax panel and the French characteristics should be visibly minor
    there, while leisure preferences move the workweek alone. Putting the three panels side by side on a
    SHARED scenario axis is what makes the rankings comparable.

    The design and ageing rows are two scenarios each, one bar per rho (OVERVIEW, _barPanel): theta = 0
    and theta = 1 in two tones on either side of the baseline, so the bar's length is the whole effect of
    the design, and acute ageing as the bar with mild ageing a line across it. The legend keys rho only;
    the paper's figure note explains the two rows.
    """
    commonX = C.US['commonX'] if commonX is None else commonX
    sfx = C.variantSuffix(commonX)
    name, marks = 'US_overview' + sfx, []
    df = D.usShocks(commonX = commonX)
    ρs = C.US['ρTable']
    rows = _rowsPresent(df, ρs, OVERVIEW, drop = not C.isLead(commonX))
    labels = [lab for lab, _, _ in rows]
    kinds = [kind for _, _, kind in rows]
    colours = RHOCOLOURS[:len(ρs)]

    fig, axes = plt.subplots(1, len(PANELS), figsize = (8.4, 0.46*len(labels) + 2.2), sharey = True)
    for ax, (col, title, xlabel, scale) in zip(axes, PANELS):
        part, unit = MARKUNIT[col]
        series, gids = [], []
        for ρ in ρs:
            base = D.usBaseline(df, ρ)[col]
            vals, ids = [], []
            for _, scen, _ in rows:
                vs = [scale*(D.usShockRow(df, ρ, s, 'full')[col] - base) for s, _ in scen]
                gs = [markId(name, ROWKEY[s], _rhoPart(ρ), part) for s, _ in scen]
                for (s, lab), v, g in zip(scen, vs, gs):
                    marks.append(mark(g, plain(lab), v, _markText(v, unit),
                                      series = 'ρ = ' + C.num(ρ, 1), panel = title))
                vals.append(tuple(vs) if len(vs) > 1 else vs[0])
                ids.append(tuple(gs) if len(gs) > 1 else gs[0])
            series.append((r'$\rho = ' + C.num(ρ, 1) + '$', vals))
            gids.append(ids)
        _barPanel(ax, title, xlabel, labels, series, colours, gids, kinds)
    _topDown(axes)
    handles, labs = axes[0].get_legend_handles_labels()
    _figLegend(fig, handles, labs)
    return _save(fig, name, marks)


# ---------------------------------------------------------------------------------------------------
# Endogenous system characteristics (app:ESC)
# ---------------------------------------------------------------------------------------------------
# The rows of section 7's figure: the categories of section 6's figure (OVERVIEW) that have a counterpart
# when the design is chosen. No pension-design row, since the design is the outcome; no mild ageing, which
# runESCcrra.py does not run (it exists at rho = 1 only). France's taste for leisure leaves the design, the
# tax and savings where they are and moves hours alone.
ESCSCENARIOS = [('Acute ageing',           'acute'),
                ('French income distr.',   'frIncome'),
                ('French voting',          'frVoting'),
                ('French leisure',         'frLeisure')]


def _escRowsPresent(df, spec, ρs, scenarios):
    """ `scenarios` after checking that each carries BOTH readings at every rho (escRow raises MissingInput, the
    output being skipped): a row of the figure is never dropped. """
    for entry in scenarios:
        for ρ in ρs:
            for pinned in (True, False):
                D.escRow(df, ρ, spec, entry[1], pinned)
    return list(scenarios)


def _escPairs(name, df, spec, ρs, scen, col, title, scale, level, idHost = ()):
    r""" One dumbbell panel's data: `pairs` for _dumbbellPanel (each reading against the baseline of the
    SAME reading at each rho, pinned against pinned and chosen against chosen, or the level for `level`),
    their `gids` and the `marks`. `idHost`: extra id parts naming the host. """
    part, unit = MARKUNIT[col]
    pairs, gids, marks = [], [], []
    for ρ in ρs:
        base = {pinned: 0. if level else float(D.escRow(df, ρ, spec, 'baseline', pinned)[col])
                for pinned in (True, False)}
        vals = [(scale*(float(D.escRow(df, ρ, spec, s, True)[col]) - base[True]),
                 scale*(float(D.escRow(df, ρ, spec, s, False)[col]) - base[False]))
                for _, s in scen]
        pairs.append((r'$\rho = ' + C.num(ρ, 1) + '$', vals))
        gids.append([])
        for (lab, s), pv in zip(scen, vals):
            ids = [markId(name, *idHost, s, _rhoPart(ρ), part, reading) for reading in ('pinned', 'chosen')]
            gids[-1].append(ids)
            marks += [mark(g, plain(lab), v, _markText(v, unit),
                           series = 'ρ = {}, θ {}'.format(C.num(ρ, 1), reading), panel = title)
                      for g, v, reading in zip(ids, pv, ('pinned', 'chosen'))]
    return pairs, gids, marks


def _dumbbellPanel(ax, title, xlabel, labels, pairs, colours, ref = 0., extraRefs = (), gids = None,
                   zero = ZEROROW):
    r""" One horizontal dumbbell panel. `pairs` is a list over series (rho values) of
    (label, [(x_pinned, x_chosen) per scenario]); each pair is drawn as an open marker at the pinned
    reading, a filled marker at the chosen one and a connector in the series colour, so the DISTANCE
    between the two readings -- the design response -- is the mark itself rather than a gap the reader
    has to measure between a bar end and a marker. `ref` is the ink reference line (0 for a deviation
    panel, the US design for the level panel); `extraRefs` are muted lines (the corners of theta).
    `gids[k][j]`: the (pinned, chosen) ids of series k's pair for scenario j. `zero`: the distance from
    `ref` within which a whole row is labelled 'unchanged' (ZEROROW for p.p. and hours, THETAZERO for the
    design's level). """
    y = np.arange(len(labels))
    off = np.linspace(-0.27, 0.27, len(pairs)) if len(pairs) > 1 else [0.]
    _panel(ax, title, '', titlesize = 12, labelsize = 10)
    ax.axvline(ref, color = INK['primary'], linewidth = 1.0, zorder = 2)
    for x in extraRefs:
        ax.axvline(x, color = INK['muted'], linewidth = 0.8, zorder = 2)
    for k, (_, vals) in enumerate(pairs):
        yy = y + off[k]
        x0 = [v[0] for v in vals]
        x1 = [v[1] for v in vals]
        for a, b, yk in zip(x0, x1, yy):
            ax.plot([a, b], [yk, yk], color = colours[k], linewidth = 1.8, zorder = 3,
                    solid_capstyle = 'round')
        # One scatter per point, so each marker carries its own gid: same z-order and draw order as one
        # scatter per reading, so not a pixel moves.
        for j, (a, yk) in enumerate(zip(x0, yy)):
            sc = ax.scatter([a], [yk], s = 34, facecolors = 'white', edgecolors = colours[k],
                            linewidths = 1.4, zorder = 5)
            if gids:
                sc.set_gid(gids[k][j][0])
        for j, (b, yk) in enumerate(zip(x1, yy)):
            sc = ax.scatter([b], [yk], s = 34, facecolors = colours[k], edgecolors = 'white',
                            linewidths = 0.8, zorder = 6)
            if gids:
                sc.set_gid(gids[k][j][1])
    for yk in y[:-1]:
        ax.axhline(yk + 0.5, color = INK['grid'], linewidth = 0.8, zorder = 1)
    # As _barPanel: a row whose every marker rounds to the reference is labelled 'unchanged'.
    for j in range(len(labels)):
        if all(abs(x - ref) < zero for _, vals in pairs for x in vals[j]):
            ax.annotate('unchanged', xy = (ref, y[j]), xytext = (7, 0), textcoords = 'offset points',
                        ha = 'left', va = 'center', fontsize = 10, style = 'italic', color = INK['muted'])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize = 11, color = INK['primary'])
    ax.set_xlabel(xlabel, color = INK['secondary'], fontsize = 10)
    ax.grid(axis = 'y', visible = False)


def escOverview():
    r""" Figure \ref{fig:US_ESC:overview}: what each counterfactual does when the pension design is
    chosen politically, against what it does at a fixed design -- as dumbbells from the pinned reading
    (open marker) to the chosen one (filled marker), one per (scenario, rho).

    Four panels in a 2 x 2 grid. The design itself first, as a LEVEL: every dumbbell there starts at
    the US design, so its length is the design response and the muted lines at 0 and 1 are the corners
    the text refers to. Then the tax rate, savings over GDP and the average workweek as deviations from
    the baseline of the same reading (_escPairs), so each marker is the counterfactual's effect in its own
    economy and the dumbbell is what choosing the design adds to it. Do not difference both readings
    against one baseline: savings and hours in 2020 respond to the design chosen for 2050, which the chosen
    baseline raises above theta*, so the two baselines differ and every pinned marker would carry the gap
    (France's taste for leisure, which leaves every rate unchanged, would read off zero). A row whose
    markers all round to the reference says 'unchanged'.

    2 x 2 rather than one row of four: at \linewidth a row of four narrow panels shrinks the type below
    legibility and squeezes the tax panel's 12-point range into an inch. The grid keeps every panel at
    half the width, and the scenario labels are repeated on the left of each row so a panel can be read
    without tracing back across the page. All four share one scenario axis (sharey), so inverting it
    once (_topDown) puts the first scenario at the top everywhere.

    The workweek is in hours, normalised against that rho's own baseline inside the ESC driver (the
    csv's ww_t0), for the reason tablesUS.py gives: under vector X the level of hbar is not identified,
    so stage (iii) never re-derives it.

    Every row is a new equilibrium path whose political choice binds from the first period, so the
    design in force in 2020 is itself an outcome. All panels are read at 2020.
    """
    df = D.escExperiments()
    spec, ρs = C.US['esc']['spec'], C.US['esc']['ρTable']
    α = D.usAlpha()
    colours = RHOCOLOURS[:len(ρs)]
    scen = _escRowsPresent(df, spec, ρs, ESCSCENARIOS)
    labels = [lab for lab, _ in scen]
    θstar = float(D.escRow(df, C.US['ρBaseline'], spec, 'baseline', False)['θ_t0'])

    # (column, title, axis label, deviation scale, is a level). sr_t0 is s/(wh); (1 - alpha) turns it
    # into s/Y before the percentage-point scaling.
    panels = [('θ_t0',  'Pension design in 2020', '$\\theta_{2020}$', 1., True),
              ('τ_t0',  'Equilibrium tax rate', 'Percentage points', 100., False),
              ('sr_t0', 'Savings over GDP', 'Percentage points of GDP', 100.*(1 - α), False),
              ('ww_t0', 'Average workweek', 'Hours per week', 1., False)]

    fig, grid = plt.subplots(2, 2, figsize = (8.0, 2*(0.5*len(labels) + 1.35) + 0.9), sharey = True)
    axes = grid.ravel()
    marks = []
    for ax, (col, title, xlabel, scale, level) in zip(axes, panels):
        pairs, gids, mk = _escPairs('US_ESC_overview', df, spec, ρs, scen, col, title, scale, level)
        marks += mk
        _dumbbellPanel(ax, title, xlabel, labels, pairs, colours,
                       ref = θstar if level else 0., extraRefs = (0., 1.) if level else (), gids = gids,
                       zero = THETAZERO if level else ZEROROW)
    axes[0].set_xlim(-0.04, 1.06)
    _topDown(axes)
    # sharey hides the inner columns' tick labels; the second ROW's left panel keeps its own, which is
    # what lets the lower panels be read without tracing back up the page.
    for ax in grid[:, 0]:
        ax.tick_params(axis = 'y', labelleft = True)

    handles = [Line2D([], [], color = c, linewidth = 1.8, label = r'$\rho = ' + C.num(ρ, 1) + '$')
               for c, ρ in zip(colours, ρs)]
    handles += [Line2D([], [], marker = 'o', linestyle = 'none', markerfacecolor = 'white',
                       markeredgecolor = INK['primary'], markersize = 6,
                       label = '$\\theta$ pinned at the U.S. design'),
                Line2D([], [], marker = 'o', linestyle = 'none', markerfacecolor = INK['primary'],
                       markeredgecolor = 'white', markersize = 6,
                       label = '$\\theta$ chosen by the electorate')]
    _figLegend(fig, handles, [h.get_label() for h in handles], bottom = 0.12, ncol = 2)
    return _save(fig, 'US_ESC_overview', marks)

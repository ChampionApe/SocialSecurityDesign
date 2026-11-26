r""" Figure builders. Each writes one pdf (and a png alongside for quick viewing) into results/paper/Figs.

House style, applied by `_panel`:
  * two categorical hues in a FIXED order, blue then orange, never cycled. The pair is validated for
    colour-vision deficiency (worst-case adjacent dE 24.7 protan / 32.7 tritan, normal 33.6) -- do not
    substitute by eye.
  * one y-axis per panel, never two. Panels carry different units on purpose; a shared axis across
    unlike measures is what a small-multiple layout exists to avoid.
  * a legend whenever a panel has two series, so identity is never carried by colour alone.
  * recessive grid and axes; text in ink colours, never in a series colour.
"""
import os, re, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

import config as C
import datasets as D

SERIES = ('#2a78d6', '#eb6834')     # categorical slots 1 and 2, in fixed order
SERIES3 = '#00a888'                # slot 3, only where a figure has three series (figure 7); with its own marker
INK    = {'primary': '#0b0b0b', 'secondary': '#52514e', 'muted': '#898781', 'grid': '#e1e0d9'}
LW, MS = 2.0, 5.0

# Sequential ramp for a CONTINUOUS parameter, one hue light->dark -- never a rainbow, and never the
# categorical pair above. Its middle step is SERIES[0], so a shaded figure and a two-series figure read
# as the same family. The lightest step is the ordinal floor against a light surface (2.06:1), so the
# palest curve stays legible in print.
THETA_RAMP = mcolors.LinearSegmentedColormap.from_list(
    'thetaBlue', ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281', '#0d366b'])


def _panel(ax, title, ylabel, titlesize = 10, labelsize = 8):
    """ `titlesize`/`labelsize` default to the Argentina figures' sizes. The US bar and dumbbell figures
    pass larger ones: they are drawn much wider than \\linewidth and lose most of their type to the
    downscaling, so they need to be set larger to arrive at the same size on the page. """
    ax.set_title(title, color = INK['primary'], fontsize = titlesize, loc = 'left', pad = 6)
    ax.set_ylabel(ylabel, color = INK['secondary'], fontsize = 9)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_color(INK['muted'])
    ax.tick_params(colors = INK['muted'], labelsize = labelsize)
    ax.grid(True, color = INK['grid'], linewidth = 0.8)
    ax.set_axisbelow(True)


def markId(*parts):
    """ The SVG id of one tagged mark: the parts joined by '--', every character outside [A-Za-z0-9_.-]
    replaced by '_'. Prefix with the figure's name so ids stay unique on a page that inlines several. """
    return re.sub(r'[^A-Za-z0-9_.-]', '_', '--'.join(str(p) for p in parts))


def mark(gid, label, value, text, series = '', panel = '', table = None, row = None):
    """ One entry of `_save`'s `marks`: the artist's gid, the row the reader sees, the number in the
    panel's own units and as the paper prints it (`signed`), and the table and row key
    (tables.rowKey) that print it, or None for both. Labels are plain unicode (`plain`). """
    return {'id': gid, 'label': label, 'series': series, 'panel': panel, 'value': float(value),
            'text': text, 'table': table, 'row': row}


MINUS = '−'


def signed(v, digits = 1, unit = ''):
    """ A change as the paper prints it (style guide §3): '+6.2 p.p.', '−1.9 hours', '+0.087'. A change
    that rounds to zero at `digits` gets up to two more digits rather than reading '0.0'; one that is
    zero at that precision too prints unsigned. """
    for d in range(digits, digits + 3):
        s = '{:.{d}f}'.format(abs(v), d = d)
        if float(s) != 0:
            break
    sign = '' if float(s) == 0 else ('+' if v > 0 else MINUS)
    return sign + s + (' ' + unit if unit else '')


_PLAIN = [('$', ''), ('\\theta', 'θ'), ('\\rho', 'ρ'), ('\\tau', 'τ'), ('\\epsilon', 'ε'), ('\\iota', 'ι'),
          ('\\xi', 'ξ'), ('\\lambda', 'λ'), ('^U', 'ᵁ'), ('\n', ' ')]


def plain(tex):
    """ A matplotlib label in the TeX subset the figures use, as plain unicode for a mark. """
    for a, b in _PLAIN:
        tex = tex.replace(a, b)
    return tex


def _save(fig, name, marks = None):
    """ Writes `name`.pdf (the paper's), .png and .svg (the supplemental appendix's) into results/paper/Figs.

    `marks`: the supplemental appendix's handle on the figure, a list of dicts, one per mark the reader can
    hover or click, each with at least 'id' (the gid the builder set on the artist, see markId), 'label',
    'value' (float) and 'text' (the value formatted as the paper prints it), and where a table carries the
    number, 'table' (its output name in build.OUTPUTS) and 'row' (the `% row:` key of that table's row).
    Written to `name`.marks.json. Every id must appear exactly once in the svg: a tag lost on the way is
    a mark the site cannot reach, so it raises. """
    os.makedirs(os.path.join(C.PAPERDIR, 'Figs'), exist_ok = True)
    out = []
    for ext in ('pdf', 'png', 'svg'):
        p = os.path.join(C.PAPERDIR, 'Figs', name + '.' + ext)
        if ext == 'svg':
            with matplotlib.rc_context({'svg.hashsalt': name, 'svg.fonttype': 'path'}):
                fig.savefig(p, metadata = {'Date': None})
        else:
            fig.savefig(p, dpi = 300)
        out.append(p)
    plt.close(fig)
    if marks is not None:
        with open(out[-1], encoding = 'utf-8') as f:
            svg = f.read()
        ids = [m['id'] for m in marks]
        if len(set(ids)) != len(ids):
            raise ValueError('{}: duplicate mark ids'.format(name))
        lost = [i for i in ids if svg.count('id="{}"'.format(i)) != 1]
        if lost:
            raise ValueError('{}: {} mark id(s) not found once in the svg, e.g. {}'.format(name, len(lost), lost[:3]))
        p = os.path.join(C.PAPERDIR, 'Figs', name + '.marks.json')
        with open(p, 'w', encoding = 'utf-8') as f:
            json.dump({'figure': name, 'marks': marks}, f, ensure_ascii = False, indent = 1)
        out.append(p)
    return out


# ---------------------------------------------------------------------------------------------------
def argLogFourInOne(commonX = None):
    r""" Figure \ref{fig:Argentina_functionOfParams}: the tax rate, the savings rate, the average
    workweek and the informal savings ratio over the (epsilon, theta) plane, at the calibration year.

    Each panel plots the variable against epsilon with ONE CURVE PER THETA and the span between adjacent
    curves shaded, so the whole surface is visible rather than two slices through it. Reads as: position
    along x is epsilon, position within the band is theta, and the band's WIDTH at a given epsilon is
    how much theta matters there.

    theta is a continuous magnitude, so it gets a ONE-HUE light-to-dark ramp and a colourbar, never a
    categorical or multi-hue scale. THETA_RAMP's middle step is SERIES[0], which keeps this figure in
    the same family as the repo's others.

    The calibrated theta is drawn over the band in INK, not as a further value of the ramp: it is an
    annotation, and giving it a hue would make the reader hunt for it on the colourbar. The reform is
    the horizontal move along that curve from the calibrated epsilon to epsilon^U. """
    commonX = C.ARG['commonX'] if commonX is None else commonX
    grid = D.epsThetaGrid(commonX = commonX)
    hRef = float(grid.loc[grid['statusQuo'], 'h'].iloc[0])
    θs   = np.sort(grid['theta'].unique())
    norm = mcolors.Normalize(float(θs[0]), float(θs[-1]))
    θBase = float(grid.loc[grid['statusQuo'], 'theta'].iloc[0])
    εBase = float(grid.loc[grid['statusQuo'], 'eps'].iloc[0])
    # eps^U is solved at every theta; the marked one is on the calibrated theta's curve.
    εU = float(grid.loc[grid['universalEps'] & np.isclose(grid['theta'], θBase), 'eps'].iloc[0])
    jBase = int(np.argmin(np.abs(θs - θBase)))

    # (title, axis label, value, mark-id part, how a mark prints the level)
    panels = [('Tax rate $\\tau$', 'percent', lambda d: 100*d['τ'], 'tau', '{:.1f}%'),
              ('Savings rate', 'percent of GDP', lambda d: 100*d['sr'], 'sr', '{:.1f}% of GDP'),
              ('Average workweek', 'hours', lambda d: C.workweekHours(d['h'], hRef), 'ww', '{:.1f} hours'),
              ('Informal savings ratio $\\iota = s^0/s$', 'ratio', lambda d: d['ι'], 'iota', '{:.3f}')]
    name = 'ARG_LOG_FourInOne' + C.variantSuffix(commonX, 'ARG')
    marks = []

    fig, axes = plt.subplots(2, 2, figsize = (9.0, 6.6), constrained_layout = True)
    for ax, (title, ylabel, f, part, fmt) in zip(axes.ravel(), panels):
        _panel(ax, title, ylabel)
        piv = grid.assign(y = f(grid)).pivot(index = 'eps', columns = 'theta', values = 'y')
        x = piv.index.values
        # Bands carry the surface; the hairlines only give it grain and mark where the solved thetas
        # actually are, so a reader does not mistake the shading for an interpolated continuum.
        for j in range(len(θs) - 1):
            ax.fill_between(x, piv.iloc[:, j], piv.iloc[:, j+1], linewidth = 0, zorder = 2,
                            color = THETA_RAMP(norm(0.5*(θs[j] + θs[j+1]))))
        for j in range(len(θs)):
            ax.plot(x, piv.iloc[:, j], color = THETA_RAMP(norm(θs[j])), linewidth = 0.4, zorder = 3)

        # A surface-coloured halo, because the calibrated theta sits in the ramp's DARK end on three of
        # the four panels and bare ink on navy is unreadable in print.
        yBase = piv.iloc[:, jBase]
        ax.plot(x, yBase, color = INK['primary'], linewidth = 1.3, zorder = 4,
                path_effects = [pe.Stroke(linewidth = 3.2, foreground = '#fcfcfb'), pe.Normal()])
        for ε, marker, which, lab in ((εBase, 'o', 'pre', 'pre-reform ε = {:.2f}'.format(εBase)),
                                      (εU, 'D', 'post', 'post-reform εᵁ = {:.2f}'.format(εU))):
            v = float(yBase.loc[ε])
            ln, = ax.plot(ε, v, marker = marker, markersize = MS, linestyle = 'none',
                          color = INK['primary'], markeredgecolor = '#fcfcfb', markeredgewidth = 1.5,
                          zorder = 5)
            ln.set_gid(markId(name, part, which))
            # The grid is a csv no table prints (the supplemental appendix offers it for download).
            marks.append(mark(ln.get_gid(), lab, v, fmt.format(v), series = 'θ = {:.2f}'.format(θBase),
                              panel = plain(title)))
        ax.set_xlabel('$\\epsilon$', color = INK['secondary'], fontsize = 9)
        ax.set_xlim(float(x[0]), float(x[-1]))

    cbar = fig.colorbar(plt.cm.ScalarMappable(cmap = THETA_RAMP, norm = norm), ax = axes,
                        location = 'top', shrink = 0.55, aspect = 32, pad = 0.02)
    cbar.set_label('$\\theta$', color = INK['secondary'], fontsize = 9)
    cbar.ax.tick_params(colors = INK['muted'], labelsize = 8)
    cbar.outline.set_visible(False)

    # Distinct marker SHAPES, not two inks: the pre/post pair must stay separable in greyscale.
    ink = lambda **kw: plt.Line2D([], [], color = INK['primary'], **kw)
    fig.legend([ink(linewidth = 1.3),
                ink(marker = 'o', linestyle = 'none', markersize = MS),
                ink(marker = 'D', linestyle = 'none', markersize = MS)],
               ['calibrated $\\theta = {:.2f}$'.format(θBase),
                'pre-reform $\\epsilon = {:.2f}$'.format(εBase),
                'post-reform $\\epsilon^U = {:.2f}$'.format(εU)],
               loc = 'outside lower center', ncol = 3, frameon = False, fontsize = 9,
               labelcolor = INK['secondary'])
    return _save(fig, name, marks)


# ---------------------------------------------------------------------------------------------------
def argCrraLog(longRun = 1, commonX = None):
    r""" Figure \ref{fig:ARG:EffectOfCRRA}: the short- and long-run effect of the reform, as a function
    of rho.

    `longRun` is a positional offset from t0 along the shock path, not a model index; the default is
    t0+1, i.e. year 2040. Any reading landing on or past the TERMINAL period is refused by the guard
    below: there s_T = 0 makes the savings rate and iota degenerate rather than small, and the period
    before it already carries the terminal boundary's influence.

    Each panel is a CHANGE in the natural unit of its variable -- percentage points for the three rates
    and hours for the workweek -- rather than a common percent scale. """
    commonX = C.ARG['commonX'] if commonX is None else commonX
    nPeriods = len(D.shockPath(C.ARG['ρBaseline'], 'reform', commonX = commonX))
    if longRun >= nPeriods - 1:
        raise ValueError('longRun={} lands on or past the terminal period (path has {} periods, and '
                         'the last is degenerate: s_T=0).'.format(longRun, nPeriods))

    short = D.shockAtPeriod(0, 'reform', commonX = commonX)
    long_ = D.shockAtPeriod(longRun, 'reform', commonX = commonX)
    year  = C.calendar()['year0']
    srShort = np.array([100*(D.reformSavingsRate(ρ, 0, commonX = commonX)
                             - D.savingsRatePath(ρ, 'base', commonX = commonX).iloc[0]) for ρ in short['ρ']])
    srLong  = np.array([100*(D.reformSavingsRate(ρ, longRun, commonX = commonX)
                             - D.savingsRatePath(ρ, 'base', commonX = commonX).iloc[longRun]) for ρ in long_['ρ']])

    # Hours are normalised against the CALIBRATED baseline at t0 -- one reference point per rho, fixed
    # across periods. Using each period's own h_base instead would re-anchor the scale every period and
    # report a change against a moving unit.
    hRef = short['h_base'].values.astype(float)
    ww = C.calendar()['workweek']

    # (title, axis label, change, precomputed change, mark-id part, unit of a mark's text)
    panels = [('Tax rate $\\tau$', 'change, percentage points',
               lambda d: 100*(d['τ_reform'] - d['τ_base']), None, 'tau', 'p.p.'),
              ('Savings rate', 'change, percentage points', None, (srShort, srLong), 'sr', 'p.p. of GDP'),
              ('Average workweek', 'change, hours',
               lambda d: ww*(d['h_reform'].values - d['h_base'].values)/hRef, None, 'ww', 'hours'),
              ('Informal savings ratio $\\iota$', 'change, percentage points',
               lambda d: 100*(d['ι_reform'] - d['ι_base']), None, 'iota', 'p.p.')]
    sfx = C.variantSuffix(commonX, 'ARG')
    name, marks = 'ARG_CRRA_LOG' + sfx, []

    fig, axes = plt.subplots(2, 2, figsize = (9.0, 6.4), constrained_layout = True)
    labels = ('short run (year {})'.format(year), 'long run (year {})'.format(year + 30*longRun))
    for ax, (title, ylabel, f, precomputed, part, unit) in zip(axes.ravel(), panels):
        _panel(ax, title, ylabel)
        ys = precomputed if precomputed is not None else (f(short), f(long_))
        for y, colour, label, run in zip(ys, SERIES, labels, ('short', 'long')):
            # The line, then each point as its own artist so it can be tagged: the same draw order as
            # one line with markers, so not a pixel moves.
            ax.plot(short['ρ'], y, color = colour, linewidth = LW, solid_capstyle = 'round', zorder = 3)
            for ρ, v in zip(short['ρ'], np.asarray(y, dtype = float)):
                pt, = ax.plot(ρ, v, color = colour, marker = 'o', markersize = MS, linestyle = 'none',
                              zorder = 3)
                pt.set_gid(markId(name, part, run, 'rho{:.1f}'.format(ρ)))
                marks.append(mark(pt.get_gid(), 'ρ = {:.1f}'.format(ρ), v, signed(v, unit = unit),
                                  series = label, panel = plain(title)))
        ax.axhline(0, color = INK['grid'], linewidth = 1, zorder = 0)
        ax.set_xlabel('$\\rho$', color = INK['secondary'], fontsize = 9)

    handles = [plt.Line2D([], [], color = colour, linewidth = LW, marker = 'o', markersize = MS,
                          solid_capstyle = 'round') for colour in SERIES]
    fig.legend(handles, labels, loc = 'outside lower center', ncol = 2, frameon = False,
               fontsize = 9, labelcolor = INK['secondary'])
    return _save(fig, name, marks)

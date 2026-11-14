r""" Appendix D's tables (app:T): the horizon and the terminal condition (NUM_Horizon), and the first order
condition's residuals off the policy functions' grids (NUM_Residuals); and the supplemental appendix's date-by-date
detail of the first (NUM_HorizonPaths). The candidate counts are tablesOA.numSelection (NUM_Selection), on the
supplemental appendix's page D2; appendix D states them in words.

Inputs: results/numerical/{US_horizon,US_horizonDesign,US_residuals}.csv
(python/US/runHorizon.py) and ARG_horizon_commonX.csv, ARG_residuals_commonX.csv
(python/InformalSavings/runHorizon.py), read through datasets.horizon/residuals (US_derivativeScheme.csv: appendix D's prose). Layout as
tablesOA: tablesUS._xwrap, a note in the house form, every body row ending in ` % row: <key>`. Rates in
percent to two decimals, their changes across horizons in p.p. to three, the design theta to three.
"""
import numpy as np

import config as C
import datasets as D
from tablesUS import _xwrap
from tablesOA import _row, _body, _signed, _rhoKey

BASE_TSS = 3                                     # the paper's horizon: three periods after the last projection
ALTERNATIVES = (('finite', 1), ('finite', 6), ('steadyState', None))
RHOS = (0.5, 1.0, 2.0)
TIGHT = (r'\footnotesize\setlength{\tabcolsep}{4pt}'      # ten columns within the measure, and the rows within a
         r'\renewcommand{\arraystretch}{1.0}')            # page: overrides _xwrap's 1.25


def _label(name):
    return 'table:T:' + name


def _pick(df, block, ρ, terminal, tss):
    """ The rows of one (block, rho, horizon), indexed by date; tss None takes the steady-state rows. """
    sub = df[(df['block'] == block) & np.isclose(df['ρ'], ρ) & (df['terminal'] == terminal)]
    if tss is not None:
        sub = sub[sub['tss'] == tss]
    if sub['tss'].nunique() != 1:
        raise ValueError('NUM_Horizon: {} rows for ({}, rho = {}, {}, {})'.format(
            'no' if sub.empty else 'several horizons of', block, ρ, terminal, tss))
    return sub.set_index('date')


def _horizonPanel(df, title, years, outcomes):
    """ One panel: a title row with the two dates under every horizon, then a row per (outcome, rho): the
    level at the paper's horizon, and each alternative's change against it at the same date. """
    t0, t1 = years
    lines = [r'\multicolumn{2}{l}{\textit{' + title + r'}} & ' + ' & '.join([str(t0), str(t1)]*4) + r' \\',
             r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}']
    for label, block, col, level, change, key in outcomes:
        for k, ρ in enumerate(RHOS):
            base = _pick(df, block, ρ, 'finite', BASE_TSS)
            cells = [label if k == 0 else '', C.num(ρ, 1)] + [level(base.loc[t, col]) for t in years]
            for terminal, tss in ALTERNATIVES:
                alt = _pick(df, block, ρ, terminal, tss)
                cells += [change(alt.loc[t, col] - base.loc[t, col]) for t in years]
            lines.append(_row(cells, key + '@' + _rhoKey(ρ)))
        lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.4em] % row:')
    return lines


def numHorizon():
    r""" Table app:T's horizon check: the tax, the savings rate and (in the U.S.) the chosen design at the
    calibration year and one period later, at the paper's horizon (T = 2200 in the U.S., 2190 in Argentina)
    and the change against it at two periods shorter, three periods longer and under the steady-state terminal
    condition, every parameter at its calibrated value. Rows keyed `<arm>_<outcome>@<rho>`. """
    us, arg = D.horizon('US'), D.horizon('ARG')
    pct = lambda v: C.num(100*v, 2)
    dpp = lambda v: _signed(100*v, 3)
    th, dth = (lambda v: C.num(v, 3)), (lambda v: _signed(v, 3))
    yUS, yARG = C.usCalendar()['year0'], C.calendar()['year0']
    lines = _horizonPanel(us, 'United States', (yUS, yUS + 30), [
        ('Tax rate (\\%), design given', 'tax', 'τ', pct, dpp, 'us_tax'),
        ('Savings rate (\\% of GDP), design given', 'tax', 'srOverY', pct, dpp, 'us_savings'),
        (r'Design $\theta$, chosen', 'design', 'θ', th, dth, 'us_design'),
        ('Tax rate (\\%), design chosen', 'design', 'τ', pct, dpp, 'us_taxchosen'),
        ('Savings rate (\\% of GDP), design chosen', 'design', 'srOverY', pct, dpp, 'us_savingschosen')])
    lines += _horizonPanel(arg, 'Argentina', (yARG, yARG + 30), [
        ('Tax rate (\\%)', 'tax', 'τ', pct, dpp, 'arg_tax'),
        ('Savings rate (\\% of GDP)', 'tax', 'srGDP', pct, dpp, 'arg_savings'),
        ('Universal pension: change in the tax rate (p.p.)', 'tax', 'dτ', lambda v: _signed(100*v, 2), dpp,
         'arg_reform')])
    Tus = int(_pick(us, 'tax', 0.5, 'finite', BASE_TSS)['T'].iloc[0])
    Targ = int(_pick(arg, 'tax', 1.0, 'finite', BASE_TSS)['T'].iloc[0])
    note = (r'\textit{Note:} The first two columns are levels at the paper' + "'" + r's horizon, $T = ' + str(Tus)
            + r'$ in the U.S.\ and $' + str(Targ) + r'$ in Argentina; the others are changes against them at the same '
            r'date, with the backward induction started two periods earlier, three periods later, or from the steady '
            r'state of appendix \ref{app:T}. Tax and savings rates in percent, their changes in percentage points, '
            r'the design in units of $\theta$. The design is chosen one period ahead at the calibrated cost; in '
            r'Argentina the universal pension is that of section \ref{sec:argentina}, as a change against the '
            r'baseline of the same horizon.')
    head = (r'\textbf{Outcome} & $\bm{\rho}$ & \multicolumn{2}{c}{\textbf{Paper' + "'" + r's $\bm{T}$}} & '
            r'\multicolumn{2}{c}{\textbf{Shorter}} & \multicolumn{2}{c}{\textbf{Longer}} & '
            r'\multicolumn{2}{c}{\textbf{Steady state}}')
    return _xwrap('NUM_Horizon', 'results/numerical/{US_horizon,US_horizonDesign,ARG_horizon'
                  + C.argVariantTag(C.ARG['commonX']) + '}.csv',
                  'The horizon and the terminal condition', _label('Horizon'),
                  r'>{\raggedright\arraybackslash}p{3.0cm}c' + 'Y'*8, head, _body(lines), note,
                  width = r'\textwidth', size = TIGHT)


ALTLABEL = {('finite', 1): 'Two periods shorter', ('finite', 6): 'Three periods longer',
            ('steadyState', None): 'Steady state'}
ALTKEY = {('finite', 1): 'shorter', ('finite', 6): 'longer', ('steadyState', None): 'steady'}
LOGTOL = 1e-9          # the U.S. log rows: horizon-free (prop:LOG:PEE) to the solver's tolerance, or refused


def _projectedYears(cal):
    """ The dated years from the calibration year to the last dated one. """
    return [y for k, y in sorted(cal['dates'].items()) if k >= cal['t0']]


def _pathsPanel(df, title, years, outcomes, ρs, arm):
    """ One panel of NUM_HorizonPaths: a title row with the dates, then per (outcome, rho, alternative) the
    change against the paper's horizon at every date. """
    lines = [r'\multicolumn{3}{l}{\textit{' + title + r'}} & ' + ' & '.join(str(y) for y in years) + r' \\']
    for label, block, col, change, key in outcomes:
        for i, ρ in enumerate(ρs):
            base = _pick(df, block, ρ, 'finite', BASE_TSS)
            for j, (terminal, tss) in enumerate(ALTERNATIVES):
                alt = _pick(df, block, ρ, terminal, tss)
                cells = [label if i == j == 0 else '', C.num(ρ, 1) if j == 0 else '', ALTLABEL[(terminal, tss)]]
                cells += [change(alt.loc[t, col] - base.loc[t, col]) if np.isfinite(alt.loc[t, col]) else '--'
                          for t in years]
                lines.append(_row(cells, '{}_{}@{}@{}'.format(arm, key, _rhoKey(ρ), ALTKEY[(terminal, tss)])))
        lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.4em] % row:')
    return lines


def numHorizonPaths():
    r""" The supplemental appendix's detail of NUM_Horizon: the change against the paper's horizon at every projected
    date from the calibration year on, under each alternative horizon. The U.S. under CRRA preferences (the
    tax with the design given; the design and the tax with it chosen), whose log rows do not depend on the
    horizon and are checked to LOGTOL rather than printed; Argentina's tax and the universal pension's change
    in it at every rho. Rows keyed `<arm>_<outcome>@<rho>@<alternative>`. """
    us, arg = D.horizon('US'), D.horizon('ARG')
    dpp = lambda v: _signed(100*v, 3)
    dth = lambda v: _signed(v, 4)
    calUS, calARG = C.usCalendar(), C.calendar()
    yUS, yARG = _projectedYears(calUS), _projectedYears(calARG)
    usOut = [('Tax rate, design given', 'tax', 'τ', dpp, 'tax'),
             (r'Design $\theta$, chosen', 'design', 'θ', dth, 'design'),
             ('Tax rate, design chosen', 'design', 'τ', dpp, 'taxchosen')]
    logGap = 0.
    for _, block, col, _, _ in usOut:
        base = _pick(us, block, 1.0, 'finite', BASE_TSS)
        for terminal, tss in ALTERNATIVES:
            alt = _pick(us, block, 1.0, terminal, tss)
            logGap = max(logGap, float(np.max(np.abs(alt.loc[yUS, col].values - base.loc[yUS, col].values))))
    if not logGap <= LOGTOL:
        raise ValueError('NUM_HorizonPaths: the U.S. log rows move with the horizon by {:.2e}'.format(logGap))
    crra = sorted(float(r) for r in us[us['preferences'] == 'CRRA']['ρ'].unique())
    lines = _pathsPanel(us, 'United States, CRRA preferences', yUS, usOut, crra, 'us')
    lines += _pathsPanel(arg, 'Argentina', yARG, [
        ('Tax rate', 'tax', 'τ', dpp, 'tax'),
        ('Universal pension: change in the tax rate', 'tax', 'dτ', dpp, 'reform')], RHOS, 'arg')
    Tus = {k: int(_pick(us, 'tax', crra[0], *k)['T'].iloc[0]) for k in (('finite', BASE_TSS), ('finite', 1), ('finite', 6))}
    Targ = {k: int(_pick(arg, 'tax', 1.0, *k)['T'].iloc[0]) for k in (('finite', BASE_TSS), ('finite', 1), ('finite', 6))}
    note = (r'\textit{Note:} The change in each outcome against the paper' + "'" + r's horizon ($T = '
            + str(Tus[('finite', BASE_TSS)]) + r'$ in the United States, $' + str(Targ[('finite', BASE_TSS)])
            + r'$ in Argentina) at the same date, every parameter at its calibrated value: with the backward '
            r'induction started two periods earlier ($T = ' + str(Tus[('finite', 1)]) + r'$ and $'
            + str(Targ[('finite', 1)]) + r'$), three periods later ($T = ' + str(Tus[('finite', 6)]) + r'$ and $'
            + str(Targ[('finite', 6)]) + r'$), and under the steady-state terminal condition of appendix '
            r'\ref{app:T}. The dates run from the calibration year to the last dated period. Changes in the tax '
            r'rates in percentage points, in the design in units of $\theta$. '
            r'Under log preferences ($\rho = 1$) the U.S.\ outcomes do not depend on the horizon: their largest '
            r'change at these dates is ' + _sci(logGap) + r', and they are not shown. The levels at the '
            r'calibration year and one period later are in the horizon table of appendix \ref{app:T}.')
    head = (r'\textbf{Outcome} & $\bm{\rho}$ & \textbf{Horizon} & \multicolumn{' + str(len(yUS)) + r'}{c}{'
            r'\textbf{Change against the paper' + "'" + r's} $\bm{T}$}')
    if len(yUS) != len(yARG):
        raise ValueError('NUM_HorizonPaths: {} U.S. and {} Argentine dates'.format(len(yUS), len(yARG)))
    return _xwrap('NUM_HorizonPaths', 'results/numerical/{US_horizon,US_horizonDesign,ARG_horizon'
                  + C.argVariantTag(C.ARG['commonX']) + '}.csv',
                  'The horizon and the terminal condition, date by date', 'table:OA:NUM_HorizonPaths',
                  r'>{\raggedright\arraybackslash}p{3.2cm}c>{\raggedright\arraybackslash}p{2.6cm}' + 'Y'*len(yUS),
                  head, _body(lines), note, width = r'\textwidth', size = r'\footnotesize')


def _sci(v):
    """ A small positive number as $a \times 10^{b}$ (one decimal); zero as $0$. """
    if not np.isfinite(v) or v == 0:
        return '$0$'
    e = int(np.floor(np.log10(abs(v))))
    return r'${:.1f} \times 10^{{{}}}$'.format(v/10**e, e)


def _pp2(v):
    """ |residual| in p.p. of the tax to two significant digits; '--' where there is none. """
    if v is None or not np.isfinite(v):
        return '--'
    x = 100*float(v)
    return '$0$' if x == 0 else '{:#.2g}'.format(x).rstrip('.') if x >= 1e-3 else _sci(x)


RESIDUALCOLS = (('path', 'visited', 'max'), ('cell', 'visited', 'max'), ('cell', 'visited', 'p99'),
                ('cell', 'all', 'max'))


def _residualCells(df, sel):
    """ The four statistics of RESIDUALCOLS for both references (the solver's condition, then the condition
    with converged tax derivatives) from the pooled rows (date 'all') of one (model, rho, ns). """
    pooled = df[sel & (df['date'] == 'all')].set_index(['set', 'region'])
    out = []
    for suffix in ('', 'Ref'):
        for s, r, stat in RESIDUALCOLS:
            out.append(_pp2(pooled.loc[(s, r), stat + suffix]) if (s, r) in pooled.index else '--')
    return out


def numResiduals():
    r""" Table app:T's accuracy check: |tau_FOC - tau_hat| in p.p. of the tax at states off the grids the
    policy functions were solved on, the 2200 (2190) horizon, pooled over periods: at the state the path
    enters each period with, at the centres of the grid's cells within the range the path visits (largest
    and 99th percentile) and anywhere on the grid (largest); against the first order condition as the solver
    evaluates it and against the condition with converged tax derivatives. The U.S. at both state grids the
    paper uses (150 nodes for every solve, 50 for the cost calibration's scan). Rows keyed
    `<arm>@<rho>[@<ns>]`. """
    us, arg = D.residuals('US'), D.residuals('ARG')
    lines = [r'\multicolumn{10}{l}{\textit{United States, CRRA preferences}} \\']
    for ns in sorted(us['ns'].dropna().unique(), reverse = True):
        for ρ in sorted(us['ρ'].unique()):
            sel = np.isclose(us['ρ'], ρ) & (us['ns'] == ns)
            lines.append(_row([C.num(ρ, 1), str(int(ns))] + _residualCells(us, sel),
                              'us@' + _rhoKey(ρ) + '@' + str(int(ns))))
    lines[-1] = lines[-1].replace(r' \\ % row:', r' \\[.4em] % row:')
    lines.append(r'\multicolumn{10}{l}{\textit{Argentina}} \\')
    for ρ in sorted(arg['ρ'].unique()):
        sub = arg[np.isclose(arg['ρ'], ρ)]
        grid = '45' + (r'$\times$45' if (sub['preferences'] == 'CRRA').all() else '')
        lines.append(_row([C.num(ρ, 1), grid] + _residualCells(arg, np.isclose(arg['ρ'], ρ)), 'arg@' + _rhoKey(ρ)))
    note = (r'\textit{Note:} Each cell is the largest (or 99th percentile) absolute difference, in percentage '
            r'points of the tax rate, between the root of the first order condition at a state and the policy '
            r'function there, over the periods of the paper' + "'" + r's horizon: at the state the equilibrium '
            r'path enters each period with, at the centers of the state grid' + "'" + r's cells within the range '
            r'the path visits (padded by a fifth), and at the centers of all its cells. The grid column gives the '
            r'nodes of the state grid. The first four columns evaluate the condition as the solver does, the last '
            r'four with its derivatives in the tax converged (appendix \ref{app:T}). Log preferences in the U.S.\ '
            r'have no state.')
    head = (r'$\bm{\rho}$ & \textbf{Grid} & \multicolumn{4}{c}{\textbf{As solved}} & \multicolumn{4}{c}{\textbf{'
            r'Converged derivatives}} \\' + '\n'
            r'\cmidrule(lr){3-6}\cmidrule(lr){7-10}' + '\n'
            r' & & Path & Visited, max & Visited, 99\% & All, max & Path & Visited, max & Visited, 99\% & '
            r'All, max')
    return _xwrap('NUM_Residuals', 'results/numerical/{US_residuals,ARG_residuals'
                  + C.argVariantTag(C.ARG['commonX']) + '}.csv',
                  'First order condition residuals off the grid', _label('Residuals'),
                  r'cc' + 'Y'*8, head, _body(lines), note, width = r'\textwidth', size = r'\footnotesize')

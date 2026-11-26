r""" Table builders for the US/France/UK arm. tables.py is the Argentina arm; the two are separate
modules for the same reason runCalibration.py and runCalibrationUS.py are.

Each function returns the complete tex body of one file in writing/Paper/Tables/, reproducing the
STRUCTURE of the hand-written table it replaces (tabularx with Y columns, the same caption, label and
rules), so the paper's \ref{}s keep resolving and a diff shows moved numbers rather than a re-layout.

Every US builder takes `commonX`, defaulting to the headline variant (config.US['commonX']), which is the
only one build.py registers; config.variantSuffix names the other variant's file and \label, so the
headline's stay exactly as the draft cites them.

THREE CONVENTIONS of stage (iii), for every builder that reads the experiment csvs (figuresUS, siteData):

  * `Savings rate` is s/Y, savings relative to GDP -- the same quantity as the Argentina tables. The US
    shock csv carries it as `srOverY` (datasets.US_SR); its `sr` column is s/(w h) and is not read. The
    ESC csvs carry only s/(w h) and go through datasets.escSavingsOverY, which is exact (w h = (1-alpha) Y).
  * A counterfactual's savings rate is reported as the CHANGE against that rho's own baseline, in
    percentage points (config.pp); only baseline rows carry a level. The baseline savings rate is a
    prediction that moves with rho (beta is identified by R), so a scenario at rho differenced against
    the rho = 1 baseline is not an effect. (The leisure row, a pure scale that cannot move savings, is
    the check: 0.00 at every rho.) Tax rates and workweeks stay as levels.
  * `Avg. workweek` is normalised against each rho's OWN baseline, inside the experiment script. Under
    vector X the level of hbar is not identified, so there is no expression that converts it to hours;
    the observed workweek is a reference point, not a unit. Stage (iii) therefore reads `workweek`
    straight out of the csv and must never re-derive it from hbar.
"""
import os
import numpy as np

import config as C
import datasets as D
from tables import BANNER, LQ, RQ, notesBlock, rowKey, keyed

# The short key of each shock-csv scenario: figuresUS's mark ids and siteData's scenario column.
ROWKEY = {'Baseline': 'baseline', r'$\theta = 0$': 'theta0', r'$\theta = 1$': 'theta1',
          'Mild ageing': 'mild', 'Acute ageing': 'acute', 'Income distribution': 'income',
          'Leisure preferences': 'leisure', 'Voting': 'voting', 'All French characteristics': 'frall',
          'France (own calibration)': 'france'}


def _xwrap(name, src, caption, label, colspec, header, body, note = None, width = r'.9\textwidth', size = None):
    """ One threeparttable around a tabularx, matching the hand-written US tables' layout. `header` is a
    list of cells for one row, or a pre-formatted string when a table needs more than one header row.

    `width`: the four-column shock tables sit at .9\\textwidth. The six-column ESC tables need the full
    measure -- at .9 each Y column is 2.2cm and the headers break mid-word. `size` (a font-size command,
    e.g. r'\\footnotesize') sets the tabular's type, not the caption's or the notes'; None leaves it at the
    paper's. """
    tn = notesBlock(note)
    head = header if isinstance(header, str) else ' & '.join(header)
    return (BANNER.format(name = name, src = src)
            + '\\begin{table}[!htb]\n\\centering\n\\begin{threeparttable}\n'
            + '\\caption{' + caption + '}\n\\label{' + label + '}\n'
            + '\\renewcommand{\\arraystretch}{1.25}\n'
            + ('' if size is None else size + '\n')
            + '\\begin{tabularx}{' + width + '}{' + colspec + '}\n\\toprule\n'
            + head + ' \\\\\n\\midrule \n'
            + body + '\n\\bottomrule\n\\end{tabularx}\n' + tn
            + '\\end{threeparttable}\n\\end{table}\n')


# ---------------------------------------------------------------------------------------------------
# LaTeX: column headers, captions, block titles. figuresOA draws the US, UK and FR entries as plain text.
COUNTRYNAME = {'US': 'U.S.', 'UK': 'UK', 'FR': 'France'}
CAPTIONNAME = dict(COUNTRYNAME, US = 'the U.S.', UK = 'the UK')    # "Household heterogeneity -- the U.S."


# The calibration table's countries, in its column order, and as its note names the ones other than the U.S.
CALIBCOUNTRIES = ('US', 'UK', 'FR')
NOTENAME = {'UK': 'the UK', 'FR': 'France'}
COUNTWORD = {2: 'two', 3: 'three'}


def _andList(items):
    return items[0] if len(items) == 1 else ', '.join(items[:-1]) + ' and ' + items[-1]


def usukfrCalibration(commonX = None, countries = CALIBCOUNTRIES, name = None, label = None):
    r""" The headline calibration of `countries` (the U.S. first, then some of CALIBCOUNTRIES in its order), built as
    output `name` with \label `label` (both before the variant's suffix; by default named and labelled by the countries,
    'USFRCalibration' and table:US:CalibUSFR). build.py registers two: the paper's USFRCalibration (table:US:CalibUSFR)
    and the supplemental appendix's three countries, OECD_Calibration (table:US:CalibUSUKFR).

    `X` is the POPULATION-WEIGHTED MEAN of X_i -- the only summary of a vector whose level IS the hours
    unit, and the one the leisure counterfactual is matched on. The lower panel holds the parameters the
    countries share, each in one cell across the country columns: beta, calibrated on the US and
    imposed on the others at the same rho, and rho, xi and alpha, imposed. A value that differs
    across the countries raises rather than printing one of them. """
    commonX = C.US['commonX'] if commonX is None else commonX
    if countries[0] != 'US' or any(k not in CALIBCOUNTRIES for k in countries) or len(countries) < 2:
        raise ValueError('usukfrCalibration: the U.S. first, then some of {}; got {}'.format(CALIBCOUNTRIES[1:],
                                                                                         countries))
    c = D.usCalibrationSummary(commonX)
    cols = [k for k in CALIBCOUNTRIES if k in countries]     # the hand-written column order
    lack = [k for k in cols if k not in c]
    if lack:
        raise D.MissingInput(os.path.join(C.PAPERDIR, 'usCalibrationSummary.csv') + ' (no {} row)'.format(lack))
    name = ''.join(cols) + 'Calibration' if name is None else name
    label = 'table:US:Calib' + ''.join(cols) if label is None else label
    others = [NOTENAME[k] for k in ('FR', 'UK') if k in cols]   # the note's order
    year0 = C.usCalendar()['year0']
    # Booktabs horizontals, plus the one deviation: a hairline at 25% black on each side of the value
    # block, which is what separates the numbers from the prose column now that there is no spanner.
    # No \addlinespace with it -- the rules are drawn per row, so a gap between rows breaks them into
    # dashes; \arraystretch carries the air instead. A merged cell carries the right-hand rule in its
    # own spec, since \multicolumn replaces the spec of the columns it spans.
    rule = '!{\\color{black!25}\\vrule width 0.5pt}'

    # Parameter, the three values, then what identifies them. The values are the table's content and sit
    # in the middle under one spanner rule, which is what separates them from the prose column; the
    # identifying phrase is a gloss and closes the row.
    def row(key, label, target, fn):
        return keyed(' & '.join([label] + [fn(c[k]) for k in cols] + [target]) + r' \\', key)

    def shared(key, label, target, col, fmt = C.num):
        v = [float(c[k][col]) for k in cols]
        if not np.allclose(v, v[0], rtol = 0, atol = 1e-12):
            raise ValueError('{} is not common to {}: {}'.format(col, cols, v))
        cell = r'\multicolumn{%d}{c%s}{%s}' % (len(cols), rule, fmt(v[0]))
        return keyed(' & '.join([label, cell, target]) + r' \\', key)

    rows = [
        # θ is the replacement-rate ratio's closed form without a measure of the design (RR0), the measure's own
        # value under config.US['theta'] (imposed as θ0)
        row('theta', r'$\theta$', 'Replacement rate dispersion' if C.US.get('theta') in (None, 'index')
            else 'Earnings-linked share of the public pension', lambda r: C.num(r['θ'], 3)),
        row('omega', r'$\omega$',
            ', '.join(r'$\tau^{' + k + '} = ' + C.pct(c[k]['τ0'], 1) + '$' for k in cols),
            lambda r: C.num(r['ω'])),
        # Two decimals: with the hours unit normalised to μ = 1 (model.addEigenVectors), the vector-X
        # X_i and their mean sit on an O(1) scale where one decimal is two significant figures.
        row('X', '$X$', 'Average workweek', lambda r: C.num(r['Xbar'], 2)),
        row('nu', r'$\nu_{%d}$' % year0, '30-year gross population growth',
            lambda r: C.num(r['ν2020'])),
        row('etaratio', r'$\eta_{H}/\eta_L$', 'Relative productivity, high to low income',
            lambda r: C.num(r['ηHηL'])),
        r'\midrule',
        # Not "imposed on the other two" as well: the note says so, and the cell then fits one line.
        shared('beta', r'$\beta$', '30-year interest rate (U.S.)', 'β'),
        shared('rho', r'$\rho$', 'Log-GHH preferences', 'ρ', fmt = '{:g}'.format),
        shared('xi', r'$\xi$', 'Frisch elasticity of labor supply', 'ξ'),
        shared('alpha', r'$\alpha$', 'Capital income share', 'α'),
    ]
    header = ([r'\textbf{Parameter}'] + [r'\textbf{' + COUNTRYNAME[k] + '}' for k in cols]
              + [r'\textbf{Identified by}'])
    # The one note that spells the variant out; every other US table points here (config.variantNote).
    n = len(cols)
    note = (r'\item \textit{Note:} $X$ is the '
            r'population-weighted mean of $X_i$; its level is the hours unit, pinned for ' + _andList(others)
            + r' by targeting average hours relative to the U.S.\ rather than in levels. The lower panel is '
            r'common to the ' + COUNTWORD[n] + r' countries; $\beta$ is calibrated for the U.S.\ and imposed on '
            + ('the other ' + COUNTWORD[n - 1] if n > 2 else others[0]) + '.'
            + C.variantNote(commonX, full = True))
    return (BANNER.format(name = name + C.variantSuffix(commonX),
                          src = 'results/paper/usCalibrationSummary.csv')
            + '\\begin{table}[!htb]\n\\centering\n\\begin{threeparttable}\n'
            + '\\caption{Calibration, ' + _andList([CAPTIONNAME[k] if k != 'UK' else 'the UK' for k in cols])
            + C.variantCaption(commonX) + '}\n'
            + '\\label{' + label + C.variantSuffix(commonX) + '}\n'
            + '\\renewcommand{\\arraystretch}{1.3}\n'
            + '\\begin{tabularx}{\\textwidth}{l' + rule + 'C{1.4cm}'*len(cols) + rule
            + '>{\\raggedright\\arraybackslash}X}\n\\toprule\n'
            + ' & '.join(header) + ' \\\\ \\midrule\n'
            + '\n'.join(rows) + '\n\\bottomrule\n\\end{tabularx}\n'
            + notesBlock(note)
            + '\\end{threeparttable}\n\\end{table}\n')


def _householdHeterogeneity(country, name, label, commonX = None):
    r""" One country's per-group table: gamma_i, X_i, eta_i, mu_i.

    Under the common-X calibration the X_i row is one number repeated, and the hours row is a prediction
    rather than the target it is under vector X -- so the `Target` column is variant-dependent. """
    commonX = C.US['commonX'] if commonX is None else commonX
    summary = D.usCalibrationSummary(commonX)
    if country not in summary:
        raise D.MissingInput(os.path.join(C.PAPERDIR, 'usCalibrationSummary.csv') + ' (no {} row)'.format(country))
    c = summary[country]
    spec = [(r'$\gamma_i$', 'γi', 2, 'Income percentiles.', 'gamma'),
            ('$X_i$',       'Xi', 2, 'Average hours worked.' if commonX else 'Hours worked.', 'X'),
            (r'$\eta_i$',   'ηi', 2, 'Income distribution.', 'eta'),
            (r'$\mu_i$',    'μi', 2, 'Voting propensity.', 'mu')]
    rows = [keyed(' & '.join([lab] + [C.num(v, d) for v in c[key]] + [target]) + r' \\', rk)
            for lab, key, d, target, rk in spec]
    return (BANNER.format(name = name + C.variantSuffix(commonX),
                          src = 'results/paper/usCalibrationSummary.csv')
            + '\\begin{table}[!htb]\n\\centering\n\\begin{threeparttable}\n'
            + '\\caption{Household heterogeneity -- ' + CAPTIONNAME[country]
            + C.variantCaption(commonX) + '}\n'
            + '\\label{' + label + C.variantSuffix(commonX) + '}\n'
            + '\\renewcommand{\\arraystretch}{1.5}\n'
            + '\\begin{tabularx}{\\textwidth}{Y|YYY|p{5cm}}\n\\hline\n'
            + '& \\multicolumn{3}{c|}{\\textbf{Income group}} & \\\\ \\cline{2-4}\n'
            + ' & '.join([r'\multicolumn{1}{c|}{\textbf{Parameter}}',
                          r'\textbf{Low}', r'\textbf{Medium}', r'\textbf{High}', r'\textbf{Target}'])
            + ' \\\\ \\hline\n' + '\n'.join(rows) + '\n\\hline\n\\end{tabularx}\n'
            + '\\end{threeparttable}\n\\end{table}\n')


def usHouseholdHeterogeneity(commonX = None):
    return _householdHeterogeneity('US', 'US_householdheterogeneity', 'table:a_US:CalibUS', commonX)


def frHouseholdHeterogeneity(commonX = None):
    return _householdHeterogeneity('FR', 'FR_householdheterogeneity', 'table:a_US:CalibFR', commonX)


def ukHouseholdHeterogeneity(commonX = None):
    return _householdHeterogeneity('UK', 'UK_householdheterogeneity', 'table:a_US:CalibUK', commonX)


# ---------------------------------------------------------------------------------------------------
# Endogenous system characteristics (section sec:esc, appendix app:escResults): the cost calibration of
# the U.S. and the UK and the cross-country test.
# ---------------------------------------------------------------------------------------------------
def _wedgeSymbol(spec = None):
    """ The calibrated cost parameter's symbol under `spec` (config.US['esc']['spec'] by default): lambda
    under 'size', p under 'scale'. Both live in the csvs' `p` column. """
    spec = C.US['esc']['spec'] if spec is None else spec
    return r'\lambda' if spec == 'size' else 'p'


def _costRows(recs, spec):
    """ The rows of a cost-calibration table, {rho: record} -> (body, thetaStar, Vtilde): rho, the cost
    parameter and omega at two decimals, and the shares of revenue lost at the observed design and at theta = 0
    (datasets.escWedge), in percent at one decimal. theta* and Vtilde are the same at every rho -- the cost
    cancels from the replacement-rate ratio, and Vtilde is a function of relative incomes alone -- so the
    note prints them once; a row that disagrees is an error in the csvs, not a cell to print. """
    rows, θs, Vs = [], [], []
    for ρ in sorted(recs):
        r = recs[ρ]
        w = D.escWedge(r, spec)
        θs.append(float(r['θStar']))
        Vs.append(w['Vtilde'])
        rows.append(keyed(' & '.join([C.num(ρ, 1), C.num(float(r['p']), 2), C.num(float(r['ω']), 2),
                                      C.pct(1 - w['fStar'], 1), C.pct(1 - w['f0'], 1)]) + r' \\',
                          rowKey('rho', ρ)))
    if max(θs) - min(θs) > 1e-6 or (None not in Vs and max(Vs) - min(Vs) > 1e-6):
        raise ValueError('theta* or Vtilde differs across rho in the cost calibration: {} {}'.format(θs, Vs))
    return '\n'.join(rows), θs[0], (None if None in Vs else Vs[0])


def _costHeader(spec):
    return [r'$\rho$', '$' + _wedgeSymbol(spec) + '$', r'$\omega$', r'$1 - f(\theta^{\ast})$', '$1 - f(0)$']


def escCalibrationTable(spec = None):
    r""" Table \ref{table:US_ESC:calibration} (section sec:esc): the calibrated cost parameter per rho (lambda
    under 'size', p under 'scale'), omega recalibrated with it, and the shares of revenue lost in 2020
    at the observed design and at theta = 0, both at the 2020 tax rate, under the cost of sections 7.1 and 7.2
    (config.US['esc']['tauPower'], datasets.escDir). `spec` defaults to the paper's (config.US['esc']['spec']);
    passing the comparison spec exercises the builder on its rows. """
    spec = C.US['esc']['spec'] if spec is None else spec
    k = C.US['esc']['tauPower']
    cal = D.escCalibration(k = k)
    recs = {}
    for ρ in C.US['esc']['ρTable']:
        if (ρ, spec) not in cal:
            raise D.MissingInput('escCalibration ({}, {})'.format(ρ, spec))
        recs[ρ] = cal[(ρ, spec)]
    body, θ, V = _costRows(recs, spec)
    if spec == 'size':
        cost = (r'the linear cost of equation \eqref{eq:esc:budget}' if k == 1
                else r'the squared cost of section \ref{sec:esc:cost}, $\tau_t^2$ in place of $\tau_t$ in \eqref{eq:esc:budget}')
        form = (r'The cost is ' + cost + r', with $\theta^{\ast} = ' + C.num(θ, 3)
                + r'$ and $\tilde V = ' + C.num(V, 3) + r'$ at every $\rho$')
    else:
        form = (r'The cost is $f(\theta) = \phi + (1-\phi)\theta^{p}$ with $\phi = ' + C.num(C.US['esc']['phi'], 1)
                + r'$ imposed, and $\theta^{\ast} = ' + C.num(θ, 3) + r'$ at every $\rho$')
    note = (r'\textit{Note:} ' + form + r'. $1 - f(\theta^{\ast})$ and $1 - f(0)$ are the shares of revenue lost '
            r'in 2020 at the observed design and at a flat benefit, both at the 2020 tax rate.'
            + C.variantNote(C.US['commonX']))
    return _xwrap('US_ESC_Calibration', 'results/' + os.path.basename(D.escDir(k)) + '/escCalibration{,CRRA}.csv',
                  'The calibrated cost of redistribution for the U.S.',
                  'table:US_ESC:calibration', 'YYYYY', _costHeader(spec), body, note,
                  width = r'.8\textwidth')


# escCountry.csv's economy keys, in the order the table prints them, their labels and row keys. The UK at U.S.
# income groups ('UKUS') stays in the csvs and is not printed.
COUNTRYROWS = (('UK', 'UK', 'gbr'), ('FR', 'France', 'fra'))
# The deep test's transports of the U.S. cost parameter, in the table's order, and the parameter's symbol: the
# paper's 'size' (λ); config.US['esc']['costVariants']'s 'sizePiR' (ε, the cost λ_n = ε Π_n/R_n) when asked for and its
# country stage (escCountryCosts.csv) is on file. ε is appendix C.2's Harberger elasticity and prints as \varepsilon:
# \epsilon is the paper's universal pension. Under the carrier 'convex' (config.CARRIER) 'sizePiR' prints κ.
TRANSPORTS = ('size',)
COSTSYMBOL = {'size': r'\lambda', 'sizePiR': C.PIRSYMBOL, 'scale': 'p'}


def _frenchCornerFalls(tol = 1e-6):
    """ The note's claim, read off escCornerSlope.csv: the 1990 electorate's slope at theta_2020 = 1 is
    negative for France at every lambda on file, and neither it nor the tax there moves with lambda. """
    df = D.escCornerSlope()
    fr = df[df['country'] == 'FR']
    fr = fr[fr['vote'] == fr['vote'].min()]
    if len(fr) < 2 or not np.isclose(fr['λ'], 0).any():
        raise D.MissingInput('escCornerSlope.csv (France at lambda = 0 and at least one other)')
    if not (fr['slope'] < 0).all():
        raise ValueError("US_ESC_Country: France's slope at the corner is not negative: {}".format(fr['slope'].tolist()))
    for col in ('slope', 'τ1Corner'):
        v = fr[col].values.astype(float)
        if np.ptp(v) > tol*max(1.0, np.abs(v).max()):
            raise ValueError("US_ESC_Country: France's {} at the corner moves with lambda: {}".format(col, v.tolist()))


def _costTwin(k):
    """ What a table under the cost with tau^k adds to its name and label (config.costSuffix) and to its caption. """
    return C.costSuffix(k), ', ' + C.COSTWORDS[k]


def escCountryTable(specs = TRANSPORTS, k = 2):
    r""" Table \ref{table:US_ESC:country}: the deep test of the cost, under the cost with tau^k (datasets.escDir;
    US_ESC_Country + config.costSuffix(k), its label too). One row per economy of COUNTRYROWS -- the UK
    and France, each at its observed design (its θ0 under the θ variant) -- with that design and Ṽ, and per transport
    of the U.S. parameter in `specs` (TRANSPORTS: its λ; its ε when asked for) the design its electorate chooses
    under it and the parameter at which it re-elects its own design (`--` where no finite value does). Where France's
    observed design is the corner theta = 1, which no proportional cost reaches, the note says so and checks it
    against escCornerSlope.csv under 'size' (_frenchCornerFalls); an economy observed below the corner may have its own
    parameter, which is then printed. Vtilde is printed when the csvs carry it (rows
    written since the 'size' spec). With ε among `specs`, each economy's Π/R of 2020 on its path under the U.S. ε
    (escCountryCosts.csv's `ΠR`, the 'US' row; the 'size' rows carry none), and the note's U.S. ratio from the ε
    calibration (datasets.escCalibrationLOG), at which the U.S. ε and λ must be one cost (checked, 1e-3 relative).
    rho = 1 (LOG) only: the country stage has no CRRA counterpart. """
    specs = list(specs)
    sfx, capTail = _costTwin(k)
    pir = C.US['esc']['costVariants']['spec']
    dfs = {spec: D.escCountry(spec = spec, k = k) for spec in specs}
    us = {spec: df[df['wedgeFrom'] == 'US'] for spec, df in dfs.items()}
    hasV = all('Vtilde' in df.columns for df in dfs.values()) and all(bool(u['Vtilde'].notna().any())
                                                                       for u in us.values())
    hasR = pir in specs
    if hasR and ('ΠR' not in us[pir].columns or us[pir]['ΠR'].isna().any()):
        raise D.MissingInput('escCountry ({}): no ΠR on a U.S.-parameter row'.format(pir))
    rows, pUS, ownFR, θFR, noOwn = [], {}, False, None, []
    for key, name, rk in COUNTRYROWS:
        cells, θs, Vs = [], [], []
        for spec in specs:
            u = us[spec][us[spec]['country'] == key]
            if u.empty:
                raise D.MissingInput('escCountry ({}, US {}, {}, phi={})'
                                     .format(key, COSTSYMBOL[spec], spec, C.US['esc']['phi']))
            u = u.iloc[0]                     # escCountry: one row per (country, wedgeFrom)
            pUS[spec] = float(u['p'])
            θs.append(float(u['θStar']))
            Vs.append(float(u['Vtilde']) if hasV else np.nan)
            own = dfs[spec][(dfs[spec]['country'] == key) & (dfs[spec]['wedgeFrom'] == 'own')]
            ownFR |= key == 'FR' and not own.empty
            if own.empty and name not in noOwn:
                noOwn.append(name)
            cells += [C.num(float(u['choice']), 3), '--' if own.empty else C.num(float(own.iloc[0]['p']), 2)]
        if np.ptp(θs) > 1e-6 or (hasV and np.ptp(Vs) > 1e-6):
            raise ValueError('US_ESC_Country: {} has another observed design or Ṽ under {}: {} {}'.format(
                key, specs, θs, Vs))
        if key == 'FR':
            θFR = θs[0]
        lead = [name, C.num(θs[0], 3)] + (['--' if Vs[0] != Vs[0] else C.num(Vs[0], 3)] if hasV else [])
        if hasR:
            lead.append(C.num(float(us[pir][us[pir]['country'] == key].iloc[0]['ΠR']), 3))
        rows.append(keyed(' & '.join(lead + cells) + r' \\', rk))
    corner = θFR is not None and abs(θFR - 1.) < 1e-12           # France observed at the Bismarckian corner
    if corner and ownFR:
        raise ValueError("US_ESC_Country: France has an 'own' cost parameter at the corner, which no value gives")
    ratio = ''
    if hasR:
        cal = D.escCalibrationLOG(pir, k = k)
        if abs(float(cal['p']) - pUS[pir]) > 1e-6*pUS[pir]:
            raise ValueError('US_ESC_Country: the U.S. {} {} of escCountry differs from its calibration {}: a stale '
                             'country stage'.format(pir, pUS[pir], cal['p']))
        rUS = float(cal['ΠR'])
        if 'size' in specs and abs(pUS[pir]*C.carrierFactor(rUS) - pUS['size']) > 1e-3*pUS['size']:
            raise ValueError('US_ESC_Country: the U.S. p·g(Π/R) {} is not its λ {}, which the note says'
                             .format(pUS[pir]*C.carrierFactor(rUS), pUS['size']))
        ratio = (r' $\Pi/R$ is the ratio of the pay-as-you-go return factor to the market return in 2020 on the '
                 r'economy\textquotesingle s own path under the U.S.\ $' + COSTSYMBOL[pir] + r'$, at which it faces '
                 + C.PIRCOST + r'; the U.S.\ ratio is $'
                 + C.num(rUS, 3) + '$' + (r', at which the U.S.\ $' + COSTSYMBOL[pir] + r'$ and $\lambda$ are the same '
                                         r'cost' if 'size' in specs else '') + '.')
    ncol = 2 + hasV + hasR
    # \mbox: in the narrow Y column the header otherwise breaks as "ob-served"
    header = (' & '*ncol + ' & '.join(r'\multicolumn{2}{c}{\textbf{U.S.} $' + COSTSYMBOL[spec] + '$}'
                                      for spec in specs) + r' \\' + '\n'
              + ''.join(r'\cmidrule(lr){%d-%d}' % (ncol + 1 + 2*i, ncol + 2 + 2*i) for i in range(len(specs))) + '\n'
              + ' & '.join([r'\textbf{Economy}', r'$\theta^{\ast}$ \mbox{\textbf{observed}}']
                           + ([r'$\tilde V$'] if hasV else []) + ([r'$\Pi/R$'] if hasR else [])
                           + sum(([r'$\theta$ \textbf{chosen}', r'\textbf{Own} $' + COSTSYMBOL[spec] + '$']
                                  for spec in specs), [])))
    pars = [r'U.S.\ $' + COSTSYMBOL[spec] + ' = ' + C.num(pUS[spec], 2) + '$' for spec in specs]
    which = 'either parameter' if len(specs) > 1 else '$' + COSTSYMBOL[specs[0]] + '$'
    if corner:
        gap = (r' No value of ' + which + r' makes France re-elect its observed design, the Bismarckian corner '
               r'$\theta = 1$')
        if 'size' in specs:
            _frenchCornerFalls()
            gap += (r': under log preferences neither the slope of the design objective at that corner nor the tax '
                    r'rate there depends on the cost parameter, and France\textquotesingle s slope there is negative.')
        else:
            gap += '.'
    else:
        gap = (r' --: no value of ' + which + ' makes ' + ' or '.join(noOwn) + r' re-elect its observed design.'
               if noOwn else '')
    note = (r'\textit{Note:} ' + C.COSTWORDS[k][0].upper() + C.COSTWORDS[k][1:] + r', equation \eqref{'
            + C.COSTEQ[k] + r'}. Each economy is its own calibration at $\rho = 1$, with the ' + ' or the '.join(pars)
            + '.' + ratio + ' ' + LQ + 'Chosen' + RQ
            + r' is the design in 2020 on the economy\textquotesingle s own freely simulated path under the '
            r'U.S.\ parameter, to be read against the observed one; ' + LQ + 'own' + RQ + r' is the value at '
            r'which that path re-elects the observed design.' + gap + C.variantNote(C.US['commonX']))
    return _xwrap('US_ESC_Country' + sfx, 'results/' + os.path.basename(D.escDir(k)) + '/escCountry.csv'
                  + (', escCountryCosts.csv' if hasR else '') + (', escCornerSlope.csv' if corner and 'size' in specs else ''),
                  r'The UK and France under the U.S.\ cost of redistribution' + capTail,
                  'table:US_ESC:country' + sfx, 'l' + 'Y'*(ncol - 1 + 2*len(specs)), header, '\n'.join(rows), note,
                  width = r'\textwidth', size = r'\setlength{\tabcolsep}{3pt}')


def ukEscCalibrationTable(spec = None, k = 2):
    r""" Table \ref{table:UK_ESC:calibration} (appendix app:escResults:calibration): the UK's own calibrated
    cost per rho at its own income groups (datasets.escCalibrationHost: rho = 1 from escCountry.csv, the
    other rho from escCalibrationCRRAUK.csv), in the layout of escCalibrationTable, under the cost with tau^k
    (datasets.escDir; UK_ESC_Calibration + config.costSuffix(k), its label too). beta is the US value without the cost
    at each rho, imposed as in appendix F (app:US); omega is recalibrated with the cost parameter. """
    spec = C.US['esc']['spec'] if spec is None else spec
    sfx, capTail = _costTwin(k)
    cal = D.escCalibrationHost('UK', spec = spec, k = k)
    recs = {}
    for ρ in C.US['esc']['ρTable']:
        if ρ not in cal:
            raise D.MissingInput("the UK's own cost calibration at rho = {} ({})".format(ρ, spec))
        recs[ρ] = cal[ρ]
    body, θ, V = _costRows(recs, spec)
    sym = _wedgeSymbol(spec)
    paper = k == C.US['esc']['tauPower']
    note = (r'\textit{Note:} ' + ('As' if paper else 'In the layout of') + r' table \ref{table:US_ESC:calibration}, for '
            r'the UK at its own income groups' + ('' if paper else r' and under ' + C.COSTWORDS[k] + r' of equation '
                                                  r'\eqref{' + C.COSTEQ[k] + '}') + r', where '
            r'$\theta^{\ast} = ' + C.num(θ, 3) + '$' + ('' if V is None else r' and $\tilde V = ' + C.num(V, 3) + '$')
            + r' at every $\rho$; $' + sym + r'$ is the value at which the UK electorate re-elects its observed '
            r'design. $\beta$ is the U.S.\ value of the calibration without the cost, imposed as in appendix '
            r'\ref{app:US}, and $\omega$ is recalibrated with $' + sym + r'$.' + C.variantNote(C.US['commonX']))
    return _xwrap('UK_ESC_Calibration' + sfx, 'results/' + os.path.basename(D.escDir(k))
                  + '/escCountry.csv, escCalibrationCRRAUK.csv', 'The calibrated cost of redistribution in the UK' + capTail,
                  'table:UK_ESC:calibration' + sfx, 'YYYYY', _costHeader(spec), body, note, width = r'.8\textwidth')

r""" The ECONOMIC-EQUILIBRIUM-ONLY reading of the universalisation reform: the new permanent eps with
taxes held at the BASELINE path, i.e. the middle row of the paper's table:Argentina:Universal. Not a test
-- an experiment script, and the companion of shockUniversal.py (which produces the "Full effect" row,
where tau is re-optimised politically).

Run:  .venv\Scripts\python.exe python\InformalSavings\shockEEOnly.py [--rho 1.0] [--control]

    Baseline              the calibrated pre-reform path                     shockUniversal.py's `_base`
    Economic Equilibrium  new permanent eps, tau fixed at the baseline path  THIS script
    Full effect           new permanent eps, tau re-optimised                shockUniversal.py's `_reform`

The decomposition is the point: the paper reports that for labour supply the two effects reinforce, while
for the savings rate the pure economic-equilibrium effect is POSITIVE and the full effect negative,
because the tax rise dominates. Each run prints the three-row table at t0 and writes the middle row.

WHY THIS IS CHEAP. Taxes are exogenous here, so there is no political problem to solve: no backward PEE
recursion, no state grids, no forward walk -- one call to EE_LOG_solve/EE_CRRA_solve (model.py §3) at the
baseline tau with the reformed eps. Seconds against the ~300 s/rho that solvePEE_CRRA costs.

THE SETUP IS shockUniversal.runShock's, ONE STEP SHORTER. The reform is dated at t0 = db['t'][db['t0']]
(the calibration year) and runs on createCopyFromt0(t0) seeded by stateAtT0, so the EE-only path starts
from exactly the seed state the full-effect path starts from and the two rows are comparable. What
changes relative to the full-effect run is only that solvePEE_* is replaced by EE_*_solve at a given tau.
eps is installed with shockUniversal.installEps -- NOT by passing eps= to the solver, which would leave
the cached db['kappa'] (and db['kappa[t-1]'] at the copy's first period) on the status quo and so solve
the household problem against the reform and the government budget against the baseline, silently.

WHERE THE BASELINE COMES FROM. By default the baseline tau path is read off the full-effect run's own CSV
(results/shocks/universal_<rule>_rho<rho>.csv, column tau_base), which is the baseline PEE's tau over the
copy's whole horizon by construction -- shockUniversal.py writes `frame(base['report'], tau=base['tau'])
.loc[t0:]`. Verified numerically at rho=1 against a re-solve (--resolveBaseline): max|dtau| = 0. Reading
it also makes the `_base` columns written here bitwise identical to the full-effect file's, which is what
lets the two rows of the table be differenced. The `_base` columns do not depend on --rule/--refType/
--scale at all, so any universal_* file for this rho supplies them.

The one thing the CSV does not carry is s_{t0-1}, the state EE_*_solve needs (it starts at t0). It is
recovered by inverting eq:factorPrices, s_{t-1} = nu h (R/alpha)^{1/(alpha-1)}, on the baseline's own
R_t0/h_t0 -- exact, not an approximation, and asserted against stateAtT0's seed under --resolveBaseline.

OUTPUT (results/shocks/eeOnly_<rule>_rho<rho>.csv) mirrors the full-effect file's layout with `_reform`
replaced by `_ee`, and adds `s_` (the lagged savings level) and `sr` (the savings rate) to both blocks:
t, rho, <series>_base, <series>_ee, d_<series> = ee/base - 1, and the tax candidates' counts nEqMax/
nCandMax/nFallback (-1 unless --resolveBaseline: no political problem is solved here). tau_ee equals
tau_base exactly by construction and is asserted on. s__base(t0) is the seed state s_{t0-1}, common to all three scenarios,
and is the number the full effect's own savings rate has to be built from.

--control solves the copy at the baseline tau AND the baseline eps before installing the reform. That is
the same EE problem the baseline path already solved over t0.., so every series must come back at the
baseline's own values; it is what catches a wrong seed state, a misaligned tau, or a mis-sliced copy, and
it costs one extra EE solve. It reproduces the baseline to <=3e-9 on every series and every period except
two entries at t0 that no copy can reproduce -- see PROXYROW0, which the full-effect CSV shares.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
os.chdir(HERE)                                  # relative paths resolve here; test.py finds data/ itself

import shockUniversal as su                     # installEps/universalEps/loadCalibrated/frame: one code path
from policy import multiplicityColumns
from gridsearch.testing import utf8Stdout

utf8Stdout()                                    # this module prints Greek; see gridsearch/testing.py

CSV, PKLDIR, OUTDIR = su.CSV, su.PKLDIR, su.OUTDIR
# The full-effect file's series, plus the lagged savings level and the savings rate. s_ is written out
# because savingsRate(s, s_, h, t) needs it and at t0 it is the SEED STATE, not a row of any file -- so a
# reader of these CSVs cannot recover it by shifting a column. It is common to all three scenarios (the
# reform is unanticipated), which is what lets the full effect's savings rate be built from s__base here
# and s_reform/h_reform out of shockUniversal.py's own CSV.
NAMES = su.SERIES + ('s_',)                     # what EE_report is asked for
SERIES = NAMES + ('sr',)                        # the written column blocks, in order
CONTROLTOL = 1e-6                               # max|Δ| the no-shock control is allowed on any series
# EE_report takes the lagged objects at its FIRST period from initialState_solve unless it is handed them
# (`lagged`), and this script does not hand them over -- the baseline csv carries no lagged state -- so on
# the copy those are the model's own initial-state proxies, not the baseline's values entering t0
# (stateAtT0's 'inherited'). Two reported series see them, at that one row only -- measured with the
# no-shock control below:
#   c20   via ι_{t0-1} (ι_/si_s_):  ~5.5% of level, under BOTH solvers.
#   bbar  via h_{t0-1} (hFromS at init['Γs']):  ~5e-4-7e-4 of level under CRRA, EXACT under LOG, where
#         Γ_{s,t-1} is a function of (τ_t,θ_t) alone and the proxy is therefore the true value.
# b0/bi do not: h_{t-1} cancels between them and bbar. Everything else reproduces to ≤3e-9 (CRRA's own
# EE root tolerance) and ≤1e-15 under LOG. The full-effect CSV's reform rows take the baseline's lagged
# state (shockUniversal passes stateAtT0's seed), so its d_c20/d_bbar at t0 are exact and this file's are
# not; no paper output reads either. Both are compared from the second row on, with the t0 gap reported.
PROXYROW0 = ('c20', 'bbar')


def readCsv(path):
    """ A shocks CSV, indexed on t, parsed exactly (see baselineFromCsv). """
    return pd.read_csv(path, index_col = 0, float_precision = 'round_trip')


def baselineFromCsv(path):
    """ The `_base` block of a shockUniversal.py CSV, suffix stripped, indexed on the model's own t (t0
    first). Column set and order match shockUniversal.SERIES.

    float_precision='round_trip' is required, not cosmetic: read_csv's default C parser is accurate only
    to ~1 ulp, so without it the baseline block written back out here would differ from shockUniversal's
    in the last digit and the two files could not be differenced exactly. """
    df = readCsv(path)
    cols = {c[:-len('_base')]: df[c] for c in df.columns if c.endswith('_base')}
    return pd.DataFrame(cols)[[k for k in su.SERIES if k+'_base' in df.columns]]


def s0FromR(m, R, h, t):
    """ s_{t-1} from eq:factorPrices, R_t = alpha (s_{t-1}/(nu_t h_t))^{alpha-1}, inverted. Used to
    recover the seed state from a baseline CSV, which starts at t0 and so has no s_{t0-1} column. """
    α, ν = m.B.get('α', t), m.B.get('ν', t)
    return ν * h * (R/α)**(1/(α-1))


def savingsRatePath(m, s, s_, h, index):
    """ Eq (calibration)'s savings rate at every t of `index`, which is the model's own CALENDAR (t0
    first), so the per-year parameter lookup must be given t explicitly rather than left to default. """
    return pd.Series([m.B.savingsRate(s[i], s_[i], h[i], t) for i, t in enumerate(index)], index = index)


def withLagAndRate(m, df, s0):
    """ Add 's_' (=s_{t-1}; s0 at the first row, since that period's lag is the seed state and not in the
    frame) and 'sr' to a frame of solved series. Where the frame already carries EE_report's own s_ the
    shift identity is CHECKED rather than overwritten -- that column is built independently, over the
    baseline's full horizon, so the two agreeing is what pins s_ at t0 against stateAtT0's seed. """
    s = df['s'].values
    lag = np.append(s0, s[:-1])
    if 's_' in df:
        dev = float(np.max(np.abs(df['s_'].values - lag)))
        if dev != 0:
            raise RuntimeError("s_ disagrees with the shifted s by {:.3e} -- seed state or time "
                               "alignment is wrong".format(dev))
    else:
        df['s_'] = lag
    df['sr'] = savingsRatePath(m, s, df['s_'].values, df['h'].values, df.index)
    return df


def eeSolve(m, preferences, τ, θ, ε, s0, **kwargs):
    """ Economic equilibrium at a GIVEN policy path, expanded into the full object set. No political
    problem is solved -- tau is data here. Returns (sol, report). """
    sol = (m.EE_LOG_solve(τ, θ, ε, s0) if preferences == 'LOG'
           else m.EE_CRRA_solve(τ, θ, ε, s0, **kwargs))
    return sol, m.EE_report(sol, τ, θ, ε, s0)


def decomposition(m, df, univ, s0):
    """ The paper's three-row table at t0 -- tau, savings rate, aggregate h -- for Baseline / Economic
    Equilibrium / Full effect, then the two signs the paper's decomposition claims. The full-effect row is
    read off shockUniversal.py's CSV; its savings rate is rebuilt here from that path's own s/h (the file
    predates the sr columns), seeded by the same s_{t0-1}, since the reform path enters t0 in the same
    state. """
    s, h = univ['s_reform'].values, univ['h_reform'].values
    srFull = savingsRatePath(m, s, np.append(s0, s[:-1]), h, df.index)
    rows = [('Baseline', df['τ_base'], df['sr_base'], df['h_base']),
            ('Economic Equilibrium', df['τ_ee'], df['sr_ee'], df['h_ee']),
            ('Full effect', univ['τ_reform'], srFull, univ['h_reform'])]
    print('\ntable:Argentina:Universal at t0 (t={}):'.format(df.index[0]))
    print('  {:<22}{:>12}{:>16}{:>16}'.format('Scenario', 'Tax rate', 'Savings rate', 'Aggregate h'))
    for name, τ, sr, hh in rows:
        print('  {:<22}{:>11.4f}%{:>15.4f}%{:>16.6f}'.format(
            name, 100*float(np.asarray(τ)[0]), 100*float(np.asarray(sr)[0]), float(np.asarray(hh)[0])))
    print('  impact-period signs:  d_sr  EE-only {:+.4%}  full {:+.4%}     d_h  EE-only {:+.4%}  full '
          '{:+.4%}'.format(df['d_sr'].iloc[0], srFull.iloc[0]/df['sr_base'].iloc[0] - 1,
                           df['d_h'].iloc[0], float(univ['d_h'].iloc[0])))


def runEEOnly(ρ, settings, rule, refType, scale, control, resolveBaseline, out, baseCsv,
              commonSettings = False, pkldir = None):
    print('\n' + '='*94)
    print('rho={}   universal pensions at t0, ECONOMIC EQUILIBRIUM ONLY (tau fixed at the baseline path):'
          '   rule={} ({})   scale={}'.format(
              ρ, rule, 'b^0 = b^{}'.format(refType) if rule == 'match' else 'eps = 1-theta', scale))
    print('='*94)
    m, preferences = su.loadCalibrated(ρ, settings, commonSettings, pkldir)
    t0 = m.db['t'][m.db['t0']]
    univPath = os.path.join(OUTDIR, baseCsv.format(ρ = ρ, rule = rule))
    univ = readCsv(univPath) if os.path.exists(univPath) else None
    print('{} solver, T={}, t0 = db index {} (calibration year)'.format(preferences, m.T, t0))

    if resolveBaseline:
        tStart = time.time()
        base = su.solvePEE(m, preferences)
        seed = m.stateAtT0(base['report'], t0, init = base['init'])
        s0 = float(seed['s0'])
        b = su.frame(base['report'], refType, names = NAMES, τ = base['τ']).loc[t0:]
        print('baseline PEE re-solved ({:.0f}s).  s entering t0 = {:.10f}'.format(time.time()-tStart, s0))
        print('  seed check: |s__base(t0) - stateAtT0 s0| = {:.2e}'.format(abs(float(b['s_'].iloc[0])-s0)))
        if univ is not None:
            print('  vs the full-effect CSV:  max|Δτ_base| = {:.2e}   |Δs_t0| = {:.2e}'.format(
                float(np.max(np.abs(b['τ'].values - univ['τ_base'].values))),
                abs(s0 - s0FromR(m, univ['R_base'].iloc[0], univ['h_base'].iloc[0], t0))))
    else:
        if univ is None:
            raise FileNotFoundError('no baseline CSV at {} -- run shockUniversal.py first, or pass '
                                    '--resolveBaseline'.format(os.path.relpath(univPath, REPO)))
        b = baselineFromCsv(univPath)
        s0 = float(s0FromR(m, b['R'].iloc[0], b['h'].iloc[0], t0))
        print('baseline read from {}.  s entering t0 = {:.10f} (inverted from R_t0)'.format(
            os.path.relpath(univPath, REPO), s0))
    b = withLagAndRate(m, b, s0)
    names = tuple(k for k in NAMES if k in b.columns)   # a baseline csv may predate a series (Y, Y0)
    print('s_ on the baseline: s__base(t0) = {:.10f} (the seed state), and s__base(t) = s_base(t-1) '
          'for t>t0 to {:.1e}'.format(float(b['s_'].iloc[0]),
                                      float(np.max(np.abs(b['s_'].values[1:] - b['s'].values[:-1])))))

    mt0 = m.createCopyFromt0(t0)
    τ, θ = b['τ'].values, mt0.db['θ'].values
    if len(τ) != mt0.T:
        raise ValueError('baseline tau has {} periods, the copy has {}'.format(len(τ), mt0.T))

    x0 = None
    if control:
        _, ctrl = eeSolve(mt0, preferences, τ, θ, mt0.db['eps'].values, s0)
        c = withLagAndRate(m, su.frame(ctrl, refType, names = names, τ = τ).set_index(b.index), s0)
        d = {k: np.abs(c[k].values - b[k].values) for k in SERIES if k in c and k in b}
        dev = {k: float(np.nanmax(v[1:] if k in PROXYROW0 else v)) for k, v in d.items()}
        print('control (copy, eps UNCHANGED, baseline tau) reproduces the baseline:  ' +
              '  '.join('max|Δ{}|={:.1e}'.format(k, v) for k, v in sorted(dev.items(), key = lambda kv: -kv[1])[:4]))
        print('  at t0 only, on initialState_solve\'s proxies (see PROXYROW0):  ' +
              '  '.join('|Δ{}|={:.1e} ({:+.2%})'.format(k, d[k][0], d[k][0]/abs(b[k].iloc[0]))
                        for k in PROXYROW0 if k in d))
        worst, wv = max(dev.items(), key = lambda kv: kv[1])
        if not (wv <= CONTROLTOL):
            raise RuntimeError('no-shock control failed: max|Δ{}| = {:.3e} > {:.1e}'.format(worst, wv, CONTROLTOL))
        x0 = mt0.x0.get('EE_CRRA')              # the no-shock solution warm-starts the reformed one

    εOld = float(mt0.db['eps'].iloc[0])
    εU = su.universalEps(mt0, rule = rule, refType = refType, scale = scale)
    su.installEps(mt0, εU)
    print('eps: {:.6f} -> {:.6f} ({:+.1f}%)   kappa: {:.6f} -> {:.6f}'.format(
        εOld, εU[0], 100*(εU[0]/εOld-1), float(m.db['κ'].iloc[t0]), float(mt0.db['κ'].iloc[0])))

    tStart = time.time()
    _, ee = eeSolve(mt0, preferences, τ, θ, mt0.db['eps'].values, s0,
                    **({} if x0 is None else {'x0': x0}))
    # The reform's defining identity on the SOLVED path (shockUniversal.runShock's check): b^0/b^ref must
    # equal the ratio the installed eps implies -- 1 under 'match'. Fails loudly if installEps left
    # db['kappa'] or db['eps[t+1]'] out of step, or refType was read against the wrong column.
    target = εU/su.relBenefit(mt0, refType)
    ratio = np.asarray(ee['b0'])/np.asarray(ee['bi'].loc[:, refType])
    print('EE solved ({:.2f}s).  b^0/b^{} = {:.6f}, target {:.6f} (max|dev| {:.2e})'.format(
        time.time()-tStart, refType, ratio[0], target[0], np.max(np.abs(ratio-target))))

    r = withLagAndRate(m, su.frame(ee, refType, names = names, τ = τ).set_index(b.index), s0)
    df = b.join(r, lsuffix = '_base', rsuffix = '_ee')
    for k in SERIES:
        if k+'_base' in df:
            df['d_'+k] = df[k+'_ee']/df[k+'_base'] - 1
    df.insert(0, 'ρ', ρ)
    df.index.name = 't'
    # The tax candidates' counts (num_robustroot.tex) of the political solves behind the row: the EE-only
    # solve has none, so only a re-solved baseline (--resolveBaseline) supplies them; -1 otherwise.
    for k, v in multiplicityColumns(*([base.get('multiplicity')] if resolveBaseline else [])).items():
        df[k] = v
    dτ = float(np.max(np.abs(df['τ_ee'].values - df['τ_base'].values)))
    if dτ != 0:                                 # tau is EXOGENOUS here: the two columns are one array
        raise RuntimeError('tau_ee != tau_base (max|Δ| = {:.3e}) -- the tax path was not held fixed'.format(dτ))
    print('tau held fixed: max|τ_ee - τ_base| = 0')

    show = ['τ_base', 'sr_base', 'sr_ee', 'd_sr', 'd_s', 'd_h', 'd_ι', 'd_c10', 'd_c20']
    print('\nrelative change (EE-only/baseline - 1), tau and the savings rate in levels:')
    print(df[[c for c in show if c in df]].to_string(float_format = lambda v: '{:9.5f}'.format(v)))

    if univ is not None:
        decomposition(m, df, univ, s0)

    if out:
        os.makedirs(OUTDIR, exist_ok = True)
        path = os.path.join(OUTDIR, out.format(ρ = ρ, rule = rule))
        df.to_csv(path)
        print('\nwritten: ' + os.path.relpath(path, REPO))
    return df


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--rho', type = float, nargs = '+', default = None,
                   help = 'default: every rho in --csv with both a pickled instance and a full-effect CSV')
    p.add_argument('--csv', default = CSV)
    p.add_argument('--rule', default = 'match', choices = ('match', 'flat'),
                   help = "as shockUniversal.py -- 'match' (default): b^0 = b^refType.  'flat': eps = 1-theta")
    p.add_argument('--refType', type = int, default = 1)
    p.add_argument('--scale', type = float, default = 1.0, help = 'multiplier on eps^U (1 = exact equality)')
    p.add_argument('--control', action = 'store_true',
                   help = 'solve the copy at the baseline eps first and require it to reproduce the '
                          'baseline (also warm-starts the reformed CRRA solve)')
    p.add_argument('--resolveBaseline', action = 'store_true',
                   help = 're-solve the baseline PEE instead of reading tau_base off the full-effect CSV '
                          '(minutes per rho under CRRA, against a file read)')
    p.add_argument('--nι', '--niota', dest = 'nι', type = int, default = 45)
    p.add_argument('--ns', type = int, default = 45)
    p.add_argument('--interpKind', default = 'cubic', choices = ('linear', 'cubic', 'pchip'))
    p.add_argument('--smoothKnots', type = int, default = 4)
    p.add_argument('--commonSettings', action = 'store_true',
                   help = 'also give LOG the grid SIZES (--nι/--ns). Diagnostic only')
    p.add_argument('--pkldir', default = PKLDIR, help = 'directory of pickled calibrated instances')
    p.add_argument('--baseCsv', default = 'universal_{rule}_rho{ρ:.4f}.csv',
                   help = "shockUniversal.py's output for this rho: the baseline block, and the "
                          "full-effect row of the printed decomposition")
    p.add_argument('--out', default = 'eeOnly_{rule}_rho{ρ:.4f}.csv', help = "'' to skip writing")
    a = p.parse_args()
    settings = {'nι': a.nι, 'ns': a.ns, 'interpKind': a.interpKind, 'smoothKnots': a.smoothKnots}
    ρs = a.rho
    if ρs is None:
        grid = pd.read_csv(a.csv)
        ρs = [float(v) for v in grid['ρ']
              if os.path.exists(os.path.join(a.pkldir, 'rho_{:.4f}.pkl'.format(float(v))))
              and os.path.exists(os.path.join(OUTDIR, a.baseCsv.format(ρ = float(v), rule = a.rule)))]
    for ρ in ρs:
        runEEOnly(ρ, settings, a.rule, a.refType, a.scale, a.control, a.resolveBaseline, a.out,
                  a.baseCsv, a.commonSettings, a.pkldir)

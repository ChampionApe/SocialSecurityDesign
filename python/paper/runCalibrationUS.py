r""" Stage (i) of the paper pipeline, US/France/UK arm: produce the calibration the rich-OECD numbers
rest on. The Argentina arm is runCalibration.py; the two are separate entry points because they delegate
to different sweep scripts, and share config.py and results/paper/.

Run:  .venv\Scripts\python.exe python\paper\runCalibrationUS.py       the MAIN part: check what exists, solve what does not
      ... --prepub                                                     the PRE-PUBLICATION part only (the exact CRRA wedge)
      ... --all                                                        both parts
      ... --force                                                      re-solve every point
      ... --summaryOnly                                                skip solving; just rebuild the summary
      ... --commonX                                                    also sweep the common-X variant
      ... --dry                                                        print the commands that would run and exit

TWO PARTS, as runShocksUS.py: the main part is the sweeps, the LOG wedge calibration and the summary; the
pre-publication part is the EXACT CRRA wedge calibration of the paper's cost and the UK's own (runESCcrra.py
--exact, the published method, ~2.5 h per rho from scratch).

THE θ VARIANT. Every command carries config.thetaFlags() (--theta config.US['theta'] --main): the U.S., France and
the UK at the variant's θ0 (data/oecd26_theta.csv), written to the standard results/ paths. A sweep point or a cost
calibration counts as present only at that thetaSource and θ0 (config.atTheta); a sweep csv holding rows of another
θ is re-solved whole (--force to the child, which refuses to resume across a θ0). A sweep point also counts only at
the settings sweepCmd passes, as the row records them (offSettings: n, ns on CRRA rows, smoothKnots, interpKind; for
France and the UK retAge and the U.S. β, βUS), so a change of config.US['gridSettings'] or ['retAgeFR'/'retAgeUK']
runs the child, which re-solves those rows (calibrateRhoGrid.resumeAt); a row without the columns counts as current
(one notice per csv). --dry checks France and the UK against the U.S. sweep on disk.

Two things happen, and only the first is expensive:

  1. The rho sweeps at config.US['ρGrid'] (--rhos). Delegated to python/US/calibrateRhoGrid.py (the US, which
     calibrates beta and omega against R and tau) and python/US/calibrateRhoGridEU.py (France and the UK,
     which IMPOSE the US beta at the same rho and calibrate omega alone). This file is the DECLARATION of
     what the paper calibrates, not a second implementation of it. A few minutes per country per variant;
     all are resumable and return already-solved rho from their own csv, so a re-run with nothing to do costs
     seconds.

     ORDER MATTERS HERE, unlike in the Argentina arm. calibrateRhoGridEU.py reads the US sweep csv for
     beta and hbar at each rho and refuses to interpolate, so the US sweep must be complete over the
     grid before any European sweep starts. That dependency is enforced, not just documented.

  2. The endogenous-theta (sec:esc) cost calibrations: the cost parameter per (rho, spec) such that the
     leaded choice at 2020 reproduces the observed design, at config.US['esc']'s phi. Delegated to
     python/US/runESC.py (LOG, rho = 1 -- the paper's spec AND the comparison spec, also the no-wedge
     corner row, and the cost variant of section 7, config.US['esc']['costVariants'], the Pi-inclusive cost,
     into its --tag csv) and python/US/runESCcrra.py (CRRA, the other rho in the esc table grid, the paper's
     spec and the UK's own, with config.escBracket's scan bracket per rho). EXPENSIVE where missing (~2.5 h
     per CRRA rho: each trial value recalibrates (beta, omega) and runs the exact recursion), which is why
     the check is per-(rho, spec) and the drivers merge into their csvs rather than overwriting them.

  3. The summary, results/paper/usCalibrationSummary.csv -- one row per country. Stage (iii) builds the
     calibration and household-heterogeneity tables from this file rather than from a pickled instance,
     so that stage stays a pure csv -> tex step with no model import (README.md). It is the only place
     in this arm that opens a pickle.

WHY THE SUMMARY IS NOT JUST THE SWEEP CSVS: they carry the calibrated scalars (beta, omega, lambda, X)
but not the per-type vectors eta_i, X_i, gamma_i, mu_i, nor theta -- and those are the whole of the two
household-heterogeneity tables and three rows of the calibration tables (USFRCalibration, OECD_Calibration). They live on the instance.
"""
import os, sys, argparse, subprocess, pickle, datetime, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

# US first -- see the docstring
SWEEPS = ('US',) + C.US['countries'][1:]


def sweepCmd(country, commonX = False, force = False):
    """ The command that produces `country`'s sweep csv, at config.US's grid and settings and its θ variant.
    --force also where the csv holds rows of another θ (stale), which the child would refuse to resume. """
    g, ρ = C.US['gridSettings'], C.US['ρGrid']
    common = (['--rhos'] + [str(r) for r in ρ]
              + ['--anchor', str(C.US['ρAnchor']), '--n', str(g['n']), '--ns', str(g['ns']),
                 '--verify', str(g['verify']), '--verifyN', str(g['verifyN']),
                 '--interpKind', g['interpKind'], '--smoothKnots', str(g['smoothKnots'])] + C.thetaFlags())
    if country == 'US':
        cmd = [C.PYTHON, os.path.join(C.USDIR, 'calibrateRhoGrid.py')] + common
    else:
        cmd = [C.PYTHON, os.path.join(C.USDIR, 'calibrateRhoGridEU.py'), '--country', country] + common
    return cmd + (['--commonX'] if commonX else []) + (['--force'] if force or stale(country, commonX) else [])


def _sweep(country, commonX = False):
    """ `country`'s sweep csv, or None while absent. """
    path = C.usSweepCsv(country, commonX)
    return pd.read_csv(path) if os.path.exists(path) else None


def stale(country, commonX = False):
    """ Does `country`'s sweep csv hold a row of another θ than config.US['theta'] and its θ0? """
    df = _sweep(country, commonX)
    return df is not None and len(C.atTheta(df, country = country)) < len(df)


GRIDCOLUMNS = ('n', 'ns', 'smoothKnots', 'interpKind')     # python/US/calibrateRhoGrid.rowSettings
INPUTCOLUMNS = ('retAge', 'βUS')                          # python/US/calibrateRhoGridEU.INPUTCOLUMNS (France, UK)
RETAGEKEY = {'FR': 'retAgeFR', 'UK': 'retAgeUK'}          # python/US/testEU.retAgeHost's config key ('nra' where absent)
_NOTICED = set()                                          # the csvs whose legacy notice was printed


def _finite(df):
    """ The rows of a sweep frame with a finite residual (the solved ones). """
    return df[np.isfinite(pd.to_numeric(df['residual'], errors = 'coerce'))]


def sweepSettings(country, commonX = False):
    """ The settings sweepCmd passes, as a sweep row records them: n, ns, smoothKnots (0 is None, the adaptive
    smoother), interpKind (calibrateRhoGrid.rowSettings), and for France and the UK retAge (config.US['retAgeFR'],
    ['retAgeUK']) and βUS, a callable of rho: the U.S. sweep's β at that rho, its last solved row at the θ variant
    (modelFR.usRefFromCsv; None where the U.S. sweep has no such row). """
    g = C.US['gridSettings']
    want = {'n': g['n'], 'ns': g['ns'], 'smoothKnots': g['smoothKnots'] or None, 'interpKind': g['interpKind']}
    if country == 'US':
        return want
    us, βUS = _sweep('US', commonX), {}
    if us is not None:
        us = _finite(C.atTheta(us, country = 'US'))
        βUS = {round(float(r), 6): float(b) for r, b in zip(us['ρ'], us['β'])}      # the last row at a rho wins
    return want | {'retAge': C.US.get(RETAGEKEY[country], 'nra'), 'βUS': lambda ρ: βUS.get(round(float(ρ), 6))}


def _same(v, w):
    """ A recorded value v (from the csv) against config's w (calibrateRhoGrid._sameSetting): None is NaN, a string
    by equality, a number to 1e-12 relative. """
    if w is None:
        return bool(pd.isna(v))
    if pd.isna(v):
        return False
    if isinstance(w, str):
        return str(v) == w
    return abs(float(v) - float(w)) <= 1e-12*max(1., abs(float(w)))


def offSettings(country, commonX = False):
    """ The sweep's own resume check (calibrateRhoGrid.staleSettings) read from the csv: over the solved rows at the θ
    variant at a rho of config.US['ρGrid'], ({rho: [the columns that differ from sweepSettings]}, [the rhos of rows
    that predate a recorded setting]). A row with n records the grid: n, smoothKnots, interpKind and, on a CRRA row,
    ns are compared; France's and the UK's retAge and βUS where the row has them (βUS not where the U.S. sweep has
    no row at that rho, which `missing('US')` reports). A row without the grid or without retAge/βUS counts as
    current, with one notice per csv. """
    df = _sweep(country, commonX)
    if df is None:
        return {}, []
    want, legacy = sweepSettings(country, commonX), (() if country == 'US' else INPUTCOLUMNS)
    on = {round(ρ, 6) for ρ in C.US['ρGrid']}
    off, old, unrecorded = {}, [], []
    for r in _finite(C.atTheta(df, country = country)).to_dict('records'):
        ρ = round(float(r['ρ']), 6)
        if ρ not in on:
            continue
        grid = not pd.isna(r.get('n', np.nan))
        cols = [c for c in GRIDCOLUMNS if grid and not (c == 'ns' and r.get('preferences') != 'CRRA')]
        cols += [c for c in legacy if not pd.isna(r.get(c, np.nan))]
        w = {c: want[c](ρ) if callable(want[c]) else want[c] for c in cols}
        diff = [c for c in cols if not (c == 'βUS' and w[c] is None) and not _same(r.get(c, np.nan), w[c])]
        lack = ([] if grid else list(GRIDCOLUMNS)) + [c for c in legacy if pd.isna(r.get(c, np.nan))]
        if diff:
            off[ρ] = diff
        elif lack:
            old.append(ρ)
            unrecorded += [c for c in lack if c not in unrecorded]
    path = C.usSweepCsv(country, commonX)
    if old and path not in _NOTICED:
        _NOTICED.add(path)
        try:
            shown = os.path.relpath(path, C.REPO)
        except ValueError:                    # another drive
            shown = path
        print('notice: {}: {} row(s) (rho {}) predate the recorded settings {}: counted as current'.format(
            shown, len(old), ', '.join('{:g}'.format(v) for v in old), ', '.join(unrecorded)))
    return off, old


def missing(country, commonX = False, settings = True):
    """ rho in the paper's grid that `country`'s sweep csv does not already carry at the θ variant, and (`settings`)
    the rho whose row records other settings than config passes (offSettings), which the child re-solves
    (calibrateRhoGrid.resumeAt). """
    df = _sweep(country, commonX)
    if df is None:
        return list(C.US['ρGrid'])
    have = set(np.round(_finite(C.atTheta(df, country = country))['ρ'], 6))
    off = set(offSettings(country, commonX)[0]) if settings else set()
    return [ρ for ρ in C.US['ρGrid'] if round(ρ, 6) not in have or round(ρ, 6) in off]


def _escLogRows(phis, name = 'escCalibration.csv'):
    """ {(spec, phi)} of converged LOG wedge calibrations in results/esc/`name`, under the headline variant and
    the θ variant (config.atTheta at the U.S. θ0, which every ESC row records). """
    path = os.path.join(C.ESCDIR, name)
    if not os.path.exists(path):
        return set()
    df = C.atTheta(pd.read_csv(path), country = 'US')
    ok = df[df['converged'].astype(bool)]
    if 'commonX' in ok.columns:
        ok = ok[ok['commonX'].astype(bool) == bool(C.US['commonX'])]
    return {(s, round(float(p), 6)) for s, p in zip(ok['spec'], pd.to_numeric(ok['phi'], errors = 'coerce'))
            if p == p}


def _exactRows(path, phi):
    """ {(rho, spec)} of converged exact CRRA cost calibrations in `path` at `phi`, headline variant. """
    if not os.path.exists(path):
        return set()
    df = C.atTheta(pd.read_csv(path), country = 'US')
    if 'method' not in df.columns:
        df['method'] = 'path'
    ok = df[df['converged'].astype(bool) & np.isclose(df['phi'], phi) & (df['method'].fillna('path') == 'exact')]
    if 'commonX' in ok.columns:
        ok = ok[ok['commonX'].astype(bool) == bool(C.US['commonX'])]
    return {(round(float(r), 6), s) for r, s in zip(ok['ρ'], ok['spec'])}


def escMissing(part = 'main', force = False, strict = True):
    """ The ESC cost-calibration commands whose rows results/esc does not already carry (converged, at
    config.US['esc']'s phi, headline variant). One command per missing combination -- the drivers merge
    into their csvs (runESC.mergeWrite), so partial re-runs are safe.

    'main': the LOG calibrations at esc['phi'] under the paper's spec and the comparison spec (about 80 s
    each; also the no-wedge row), and the cost variant's (esc['costVariants'], its own --tag csv).
    'prepub': the EXACT CRRA calibrations of the paper's spec at every rho the ESC tables print (method = 'exact'
    rows -- the path iteration's rows in the same csv do not count), one command per rho with its own scan bracket
    from config.escBracket, and the UK's own (the deep test). `force` lists them all. Every command carries
    config.thetaFlags() and only rows at the θ variant count. A placeholder bracket raises unless `strict` is False
    (the --dry listing), which prints it as `--bracket <unset> <unset>` instead. """
    esc = C.US['esc']
    spec, phi = esc['spec'], esc['phi']
    specsL = [spec] + ([esc['comparisonSpec']] if esc.get('comparisonSpec') else [])
    cv = esc['costVariants']
    # The ESC leg runs under the headline calibration variant only -- see paper/runShocksUS.ESCVARIANT.
    variant = (['--commonX'] if C.US['commonX'] else []) + C.thetaFlags()
    cmds = []
    haveL = set() if force else _escLogRows([phi])
    lackL = [s for s in specsL if (s, round(phi, 6)) not in haveL]
    if part in ('main', 'all') and lackL:
        cmds.append([C.PYTHON, os.path.join(C.USDIR, 'runESC.py'), '--stage', 'calib',
                     '--spec'] + lackL + ['--phi', str(phi)] + variant)
    if part in ('main', 'all') and (force or (cv['spec'], round(phi, 6)) not in
                                    _escLogRows([phi], 'escCalibration' + cv['tag'] + '.csv')):
        cmds.append([C.PYTHON, os.path.join(C.USDIR, 'runESC.py'), '--stage', 'calib', '--spec', cv['spec'],
                     '--phi', str(phi), '--tag', cv['tag']] + variant)
    if part in ('prepub', 'all'):
        exact = [C.PYTHON, os.path.join(C.USDIR, 'runESCcrra.py'), '--exact', '--stage', 'calib']
        # the exact recursion's design layer (config.US['esc']), on every runESCcrra.py command below
        layer = ['--designRule', esc['designRule'], '--Ma', str(esc['Ma'])]
        # the paper's spec into escCalibrationCRRA.csv
        for sp, tag in [(spec, '')]:
            have = set() if force else _exactRows(os.path.join(C.ESCDIR, 'escCalibrationCRRA' + tag + '.csv'), phi)
            lackρ = [ρ for ρ in esc['ρTable'] if ρ != C.US['ρAnchor'] and (round(ρ, 6), sp) not in have]
            scan = ['--nsScan', str(esc['nsScan'])]
            opts = (['--spec', sp, '--phi', str(phi), '--ns', str(esc['ns2D'])] + scan
                    + ['--nCand2D', str(esc['nCand2D'])] + (['--tag', tag] if tag else []) + variant + layer)
            # A rho with its own bracket gets its own command (the bracket is one pair); the rest share one
            # command and runESCcrra.py's own default scan for the spec.
            brackets = {}
            for ρ in lackρ:
                try:
                    brackets[ρ] = C.escBracket(ρ, sp)
                except ValueError:
                    if strict:
                        raise
                    brackets[ρ] = ('<unset>', '<unset>')
            shared = [ρ for ρ in lackρ if brackets[ρ] is None]
            if shared:
                cmds.append(exact + ['--rho'] + [str(r) for r in shared] + opts)
            for ρ in lackρ:
                if brackets[ρ] is not None:
                    xt = [v for k, v in esc.get('xtol', {}).get(sp, {}).items() if np.isclose(float(k), ρ)]
                    cmds.append(exact + ['--rho', str(ρ), '--bracket'] + [str(b) for b in brackets[ρ]] + opts
                                + (['--xtol', str(xt[0])] if xt and not tag else []))
        # The UK's own exact CRRA cost (the deep test, UK_ESC_Calibration), one command per missing rho, its
        # bracket and scan grid from esc['uk'].
        uk = esc['uk']
        pathUK = os.path.join(C.ESCDIR, 'escCalibrationCRRAUK.csv')
        haveUK = set()
        if os.path.exists(pathUK) and not force:
            df = C.atTheta(pd.read_csv(pathUK), country = 'US')
            ok = df[df['converged'].astype(bool) & (df['method'] == 'exact') & (df['spec'] == spec)
                    & (df['commonX'].astype(bool) == bool(C.US['commonX']))]
            haveUK = {round(float(r), 6) for r in ok['ρ']}
        for ρ in esc['ρTable']:
            if ρ == C.US['ρAnchor'] or round(ρ, 6) in haveUK:
                continue
            cmds.append([C.PYTHON, os.path.join(C.USDIR, 'runESCcrra.py'), '--exact', '--host', 'UK',
                         '--stage', 'calib', '--rho', str(ρ), '--bracket']
                        + [str(b) for b in uk['bracket'][ρ]]
                        + ['--nScan', str(uk['nScan']), '--spec', spec, '--phi', str(phi),
                           '--ns', str(esc['ns2D']), '--nsScan', str(uk['nsScan'][ρ]),
                           '--nCand2D', str(esc['nCand2D'])] + variant + layer)
    return cmds


def laneOf(cmd):
    """ The csv a cost-calibration command of escMissing writes (its stem): escCalibration<tag> (runESC.py) or
    escCalibrationCRRA<host><tag> (runESCcrra.py). Commands of one lane must run in turn. """
    opt = lambda k: cmd[cmd.index(k) + 1] if k in cmd else ''
    host = opt('--host')
    crra = os.path.basename(cmd[1]) == 'runESCcrra.py'
    return 'escCalibration' + ('CRRA' + ('' if host in ('', 'US') else host) if crra else '') + opt('--tag')


def summarise(country, ρ = None, commonX = False):
    """ One country's calibration at `rho`, as one flat record.

    Unpickling needs ModelUS/ModelFR importable and the workbook path resolvable, which is why this
    chdir's into python/US exactly as the experiment scripts do. Held to this one function so nothing
    else in the pipeline inherits the requirement.
    """
    ρ = C.US['ρBaseline'] if ρ is None else ρ
    cwd = os.getcwd()
    sys.path.insert(0, C.USDIR)
    os.chdir(C.USDIR)
    try:
        with open(os.path.join(C.usInstanceDir(country, commonX), 'rho_{:.4f}.pkl'.format(ρ)), 'rb') as f:
            m = pickle.load(f)
        t0 = m.db['t'][m.db['t0']]
        atT0 = lambda k: float(m.db[k].xs(t0)) if hasattr(m.db[k], 'xs') else float(np.asarray(m.db[k])[0])
        # Vectors go through the csv as JSON, not as a python repr -- see runCalibration.summarise for
        # the numpy-2 'np.float64(...)' reading that made this a live bug rather than a hypothetical one.
        vec = lambda a: json.dumps([float(v) for v in np.asarray(a, dtype = float).ravel()])
        ηi = m.db['ηi'].xs(t0).values.astype(float)
        Xi = m.db['Xi'].xs(t0).values.astype(float)
        γi = m.db['γi'].xs(t0).values.astype(float)
        μi = m.db['μi'].xs(t0).values.astype(float)
        νt = m.db['ν'].values.astype(float)
        rec = {'country': country, 'ρ': ρ, 'commonX': bool(commonX),
               'preferences': m._calPreferences(), 't0': int(m.db['t0']), 'T': int(m.T),
               'θ': atT0('θ'), 'α': atT0('α'), 'ξ': atT0('ξ'), 'Γh': atT0('Γh'),
               'ν2020': float(νt[m.db['t0']]),
               # The population-weighted mean X_i is the row the calibration tables report as `X`, and the
               # only summary of X_i a pure scale can be matched on (python/US/shocks.shockLeisure).
               'Xbar': float((γi*Xi).sum()),
               'ηHηL': float(ηi[-1]/ηi[0]),
               'ηi': vec(ηi), 'Xi': vec(Xi), 'γi': vec(γi), 'μi': vec(μi), 'ν': vec(νt)}
    finally:
        os.chdir(cwd)
        sys.path.remove(C.USDIR)

    # beta/omega/lambda come from the sweep csv, not the pickle: the csv is the record of the calibration
    # and is what every other paper number is read against, so the table must not be able to disagree
    # with it. theta is on the instance only -- it is derived, not searched.
    df = C.atTheta(pd.read_csv(C.usSweepCsv(country, commonX)), country = country)
    row = df.loc[np.isclose(df['ρ'], ρ)]
    if row.empty:
        raise SystemExit('rho={} at θ variant {} is not in {} -- run the sweep first.'
                         .format(ρ, C.US['theta'], os.path.relpath(C.usSweepCsv(country, commonX), C.REPO)))
    row = row.iloc[-1]
    if not np.isclose(rec['θ'], float(row['θ0']), rtol = 0, atol = 1e-12):
        raise SystemExit('{}: the pickle at rho={} has θ {} but its sweep row θ0 {}: pickle and csv from different '
                         'runs'.format(country, ρ, rec['θ'], row['θ0']))
    for k in ('β', 'ω', 'λ', 'X', 'R', 'τ', 'sr', 'h', 'hbar', 'residual', 'verifyResidual',
              'hoursDrift', 'thetaSource', 'θ0', 'commit', 'timestamp'):
        rec[k] = row[k] if k in row.index else np.nan
    # The workweek and the tax target are inputs, and tau is a target: a mismatch here means the pickle
    # and the csv row came from different runs, or that the workbook moved under the sweep.
    cal = C.usCalendar(country)
    rec['workweek'] = cal['workweek']
    rec['τ0'] = cal['τ0']
    if not np.isclose(float(rec['τ']), cal['τ0'], atol = 1e-3):
        raise SystemExit('{}: solved tau={} but the workbook targets {}. Sweep and workbook disagree.'
                         .format(country, rec['τ'], cal['τ0']))
    rec['builtAt'] = datetime.datetime.now().replace(microsecond = 0).isoformat()
    return rec


def main():
    p = argparse.ArgumentParser(description = __doc__.split('\n')[1])
    p.add_argument('--force', action = 'store_true', help = 're-solve every point, not only the missing')
    p.add_argument('--summaryOnly', action = 'store_true', help = 'rebuild the summary from what exists')
    p.add_argument('--commonX', action = 'store_true',
                   help = 'also SWEEP the common-X variant. The summary always covers both variants -- '
                          'the paper builds a headline table and its twin from them and summarising is '
                          'only an unpickle -- so this flag is about the expensive step alone.')
    p.add_argument('--dry', action = 'store_true', help = 'print the sweep commands and exit')
    p.add_argument('--rho', type = float, default = None, help = 'summarise a rho other than the baseline')
    p.add_argument('--prepub', action = 'store_true',
                   help = 'the pre-publication part only: the exact CRRA wedge, no sweeps')
    p.add_argument('--all', dest = 'all_', action = 'store_true', help = 'both parts')
    p.add_argument('--lanes', default = None, metavar = 'LOGDIR',
                   help = 'run the cost calibrations that write different csvs side by side (one process per csv, '
                          "each csv's commands in turn; laneOf), each lane logging to LOGDIR/escLane_<csv>.log")
    p.add_argument('--keepGoing', action = 'store_true',
                   help = 'after a failed sweep or cost calibration, go on with the others and the summary; exit 1 '
                          'at the end if any failed')
    a = p.parse_args()
    part = 'all' if a.all_ else ('prepub' if a.prepub else 'main')
    failed = []

    def fail(msg):
        if not a.keepGoing:
            raise SystemExit(msg)
        print('FAILED ' + msg)
        failed.append(msg)

    variants = [False] + ([True] if a.commonX else [])
    if a.dry:
        if part != 'prepub':
            for cx in variants:
                for c in SWEEPS:
                    if a.force or missing(c, cx):
                        print(' '.join(sweepCmd(c, cx, a.force)))
                    else:
                        print('# sweep {}{}: all {} rho current, not run'.format(
                            c, ' (common X)' if cx else '', len(C.US['ρGrid'])))
        for cmd in escMissing(part, a.force, strict = False):
            print(' '.join(cmd))
        return

    if not a.summaryOnly and part != 'prepub':
        for cx in variants:
            for c in SWEEPS:
                todo = C.US['ρGrid'] if a.force else missing(c, cx)
                label = c + (' (common X)' if cx else '')
                if not todo:
                    print('{:<14} all {} rho already solved.'.format(label, len(C.US['ρGrid'])))
                    continue
                # The US sweep is the reference every European one reads beta and hbar from, and
                # calibrateRhoGridEU refuses to interpolate it -- so an incomplete US sweep must stop the
                # run here rather than fail per-point halfway through a march.
                if c != 'US' and missing('US', cx):
                    fail('The US sweep is incomplete at {}; France/UK impose its beta and '
                         'cannot be swept first ({}).'.format(missing('US', cx), label))
                    continue
                print('\n{}: calibrating {} point(s)'.format(label, len(todo)))
                cmd = sweepCmd(c, cx, a.force)
                print('  ' + ' '.join(cmd))
                r = subprocess.run(cmd, cwd = C.REPO)
                if r.returncode:
                    fail('sweep for {} exited {}'.format(label, r.returncode))

    if not a.summaryOnly:
        # --- the ESC cost calibrations (sec:esc). Per-(rho, spec, method) check, so a complete
        # results/esc costs nothing here; a missing exact CRRA rho costs ~45 min (pre-publication part).
        # A placeholder bracket (config.escBracket) raises here, before any command runs.
        cmds = escMissing(part, a.force)
        if not cmds:
            print('ESC cost ({}): every combination already calibrated.'.format(part))
        if a.lanes:
            # one lane per csv the commands write (laneOf), each lane's commands in turn, the lanes side by side:
            # no two processes ever merge into one csv (runESC.mergeWrite reads, then writes)
            lanes = {}
            for cmd in cmds:
                lanes.setdefault(laneOf(cmd), []).append(cmd)
            os.makedirs(a.lanes, exist_ok = True)
            script = ('import subprocess, sys, json\ncmds = json.loads(sys.argv[1])\nbad = 0\n'
                      'for c in cmds:\n    print("ESC cost: " + " ".join(c), flush = True)\n'
                      '    bad += subprocess.run(c).returncode != 0\nsys.exit(1 if bad else 0)\n')
            procs = {}
            for lane, cs in lanes.items():
                log = open(os.path.join(a.lanes, 'escLane_{}.log'.format(lane)), 'w', encoding = 'utf-8')
                print('ESC cost lane {} ({} command(s)) -> {}'.format(lane, len(cs), log.name))
                procs[lane] = (subprocess.Popen([C.PYTHON, '-c', script, json.dumps(cs)], cwd = C.REPO,
                                                stdout = log, stderr = subprocess.STDOUT), log)
            for lane, (pr, log) in procs.items():
                if pr.wait():
                    fail('ESC calibration lane {} exited {} ({})'.format(lane, pr.returncode, log.name))
                log.close()
            cmds = []
        for cmd in cmds:
            print('\nESC cost: ' + ' '.join(cmd))
            r = subprocess.run(cmd, cwd = C.REPO)
            if r.returncode:
                fail('ESC calibration exited {}: {}'.format(r.returncode, ' '.join(cmd[1:])))

    # Both variants, always: the paper's headline outputs read one and their robustness twins the
    # other, so a summary carrying only one of them blocks half the build (config.US['commonX']).
    # A sweep that does not exist yet is left out of the summary rather than failing it; its tables
    # report MissingInput in stage (iii). The summary reads what is on file (missing without the settings check).
    recs = [summarise(c, a.rho, cx) for cx in (False, True) for c in SWEEPS
            if os.path.exists(C.usSweepCsv(c, cx))
            and (C.US['ρBaseline'] if a.rho is None else a.rho) not in missing(c, cx, settings = False)]
    os.makedirs(C.PAPERDIR, exist_ok = True)
    out = os.path.join(C.PAPERDIR, 'usCalibrationSummary.csv')
    pd.DataFrame(recs).to_csv(out, index = False)
    print('\nwritten: ' + os.path.relpath(out, C.REPO))
    # Printed for the variant the paper leads with. beta, omega and theta are common to both by block
    # recursivity; Xbar and the eta ratio are not, which is the whole content of the variant.
    print('  the {} calibration (config.US[\'commonX\'] = {})'.format(
        'common-X' if C.US['commonX'] else 'vector-X', C.US['commonX']))
    print('  {:<8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(
        'country', 'θ', 'ω', 'β', 'Xbar', 'ηH/ηL', 'ν2020'))
    for r in recs:
        if bool(r['commonX']) != bool(C.US['commonX']):
            continue
        print('  {:<8} {:8.4f} {:8.4f} {:8.4f} {:8.2f} {:8.3f} {:8.3f}'.format(
            r['country'], r['θ'], r['ω'], r['β'], r['Xbar'], r['ηHηL'], r['ν2020']))
    if failed:
        raise SystemExit('failed: ' + ' | '.join(failed))


if __name__ == '__main__':
    main()

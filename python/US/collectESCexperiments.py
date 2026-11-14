r""" Collect the endogenous-theta counterfactual runs across rho into one table.

Run:  .venv\Scripts\python.exe python\US\collectESCexperiments.py --theta V
      ... --costs        the ageing counterfactuals under the three costs -> escAgeingCosts.csv (see costs())
      ... --theta V      only the rows at thetaSource V (a csv or row without the column reads as 'index'; --theta
                         index keeps those), which the merged csv then carries (atTheta). Without it the run is
                         refused, since results/esc/ holds the main variant's rows (oecd26.underResults); an ESC
                         directory outside results/ merges every row, as written.

Inputs (results/esc/): escShocks.csv (rho = 1, LOG, runESC.py) and escShocksCRRA.csv (rho != 1, CRRA,
runESCcrra.py --stage shocks). Writes results/esc/escExperiments.csv (long form, one row per
(rho, spec, method, scenario, reading)) and prints a per-spec pivot at t0 (2020): the design in force
(theta_t0) and the tax/savings outcomes there, chosen vs pinned. Reported at t0 because every scenario
is a new equilibrium path whose political choice binds from the first period, so 2020's design is an
outcome rather than an inherited datum (runESC.py's shocks-stage docstring). The paper pipeline's
stage (ii) declares this script (python/paper/runShocksUS.py) and stage (iii) reads only its output, so
a change to the merge belongs here, not in the pipeline.

`method` says which solver produced a CRRA row: 'exact' (LeadedCRRA2D, the published method) or 'path'
(the path iteration). LOG rows are stamped 'exact' -- the LOG backward recursion is the exact solution
-- and a CRRA row written before the column existed is a path-iteration row. Both vintages are carried
through; python/paper/datasets.escRow selects by config.US['esc']['exact'].
"""
import os, sys, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
ESC = os.path.join(REPO, 'results', 'esc')

COLUMNS = ('ρ', 'preferences', 'spec', 'phi', 'commonX', 'method', 'p', 'scenario', 'θpinned',
           'θ_tm1', 'θ_t0', 'θ_t1', 'τ_t0', 'sr_t0', 'ww_t0', 'τ_t1', 'sr_t1', 'ww_t1', 'thetaSource', 'θ0')


def atTheta(d, theta):
    """ The rows of frame `d` at thetaSource `theta` (a missing column or cell reads as 'index'); `d` itself when
    theta is None. """
    if theta is None or d is None:
        return d
    ts = d['thetaSource'].fillna('index') if 'thetaSource' in d.columns else pd.Series('index', index = d.index)
    return d[ts == theta]


def merge(log, crra):
    """ The merged long table from the LOG frame (escShocks.csv) and the CRRA frame (escShocksCRRA.csv).
    Pure: no file access, so test_esc.py can drive it on toy frames. """
    log = log.copy()
    crra = crra.copy()
    log['ρ'], log['preferences'] = 1.0, 'LOG'
    log['method'] = 'exact'
    crra['preferences'] = 'CRRA'
    if 'method' not in crra.columns:
        crra['method'] = 'path'
    crra['method'] = crra['method'].fillna('path')
    d = pd.concat([log, crra], ignore_index = True)
    # commonX identifies the calibration variant the row was produced under and is carried through to
    # the merged csv, where datasets.escExperiments filters on it. A row written before the column
    # existed is a vector-X row (the variant everything ran under then), so that is what a missing value
    # means.
    if 'commonX' not in d.columns:
        d['commonX'] = False
    d['commonX'] = d['commonX'].fillna(False).astype(bool)
    d = d[[c for c in COLUMNS if c in d.columns]]
    return d.sort_values(['commonX', 'method', 'spec', 'ρ', 'scenario', 'θpinned']).reset_index(drop = True)


COSTS = ('size', 'scale', 'sizePiR')          # the paper's cost, the cost on the design, the Pi-inclusive one
AGEING = ('baseline', 'mild', 'acute')
OUTCOMES = ('θ_t0', 'τ_t0', 'sr_t0', 'ww_t0', 'θ_t1', 'τ_t1')


def costs(frames):
    """ The ageing counterfactuals (mild: nu -> (1+nu)/2, acute: nu -> 1; shocks.shockAgeing) under the three
    costs of Eq esc:AB at every rho, both readings, as one long table. `frames` maps file stems to frames:
    the paper's csvs ('escShocks', 'escShocksCRRA') and the cost-variant runs ('escShocksCosts',
    'escShocksCRRACosts': 'sizePiR' at every rho; 'escShocksCRRACostsMild': mild ageing under 'size' and
    'scale' at rho != 1). One
    row per (rho, spec, scenario, θpinned) at the published method (LOG rows are exact), common X, phi 0.5,
    with the outcomes at t0 and t0+1, their changes against the baseline of the SAME reading and spec (d*),
    and pi = Pi/R at t0 where the spec carries it (with its carrier where not 'pi'). A key found in
    two files is an error. Pure. """
    parts = []
    for stem, d in frames.items():
        if d is None or d.empty:
            continue
        d = d.copy()
        if 'ρ' not in d.columns:
            d['ρ'] = 1.0
        if 'method' not in d.columns:
            d['method'] = np.nan
        d['method'] = d['method'].fillna('path' if 'CRRA' in stem else 'exact')
        d = d[(d['method'] == 'exact') & d['spec'].isin(COSTS) & d['scenario'].isin(AGEING)
              & d['commonX'].fillna(False).astype(bool) & np.isclose(d['phi'], 0.5)]
        d['source'] = stem
        parts.append(d)
    d = pd.concat(parts, ignore_index = True)
    key = ['ρ', 'spec', 'scenario', 'θpinned']
    dup = d.duplicated(key, keep = False)
    if dup.any():
        raise ValueError('escAgeingCosts: {} rows on {} come from two files: {}'.format(
            int(dup.sum()), key, sorted(set(d.loc[dup, 'source']))))
    base = d[d['scenario'] == 'baseline'].set_index(['ρ', 'spec', 'θpinned'])
    at = pd.MultiIndex.from_frame(d[['ρ', 'spec', 'θpinned']])
    for c in OUTCOMES:
        d['d' + c] = d[c].values - base[c].reindex(at).values
    cols = (['ρ', 'spec', 'p', 'scenario', 'θpinned', 'method', 'source'] + list(OUTCOMES)
            + ['d' + c for c in OUTCOMES] + [c for c in ('ΠR', 'carrier') if c in d.columns]
            + [c for c in ('thetaSource', 'θ0') if c in d.columns])
    return d[cols].sort_values(['spec', 'ρ', 'scenario', 'θpinned']).reset_index(drop = True)


def main():
    sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
    ap = argparse.ArgumentParser()
    ap.add_argument('--costs', action = 'store_true',
                    help = "write escAgeingCosts.csv from the paper's csvs and the --tag Costs runs instead")
    ap.add_argument('--theta', default = None, help = "keep the rows at this thetaSource only (atTheta)")
    args = ap.parse_args()
    target = os.path.join(ESC, 'escAgeingCosts.csv' if args.costs else 'escExperiments.csv')
    if args.theta is None:
        import oecd26
        if oecd26.underResults(target):
            raise SystemExit("collectESCexperiments.py: without --theta every row on file, whatever its θ, would be "
                             "merged into {}, which holds the main variant's (--theta config.US['theta']): pass "
                             "--theta V (--theta index: the rows without thetaSource)".format(os.path.relpath(target, REPO)))
    if args.costs:
        stems = ('escShocks', 'escShocksCRRA', 'escShocksCosts', 'escShocksCRRACosts', 'escShocksCRRACostsMild')
        frames = {st: atTheta(pd.read_csv(os.path.join(ESC, st + '.csv')), args.theta) for st in stems
                  if os.path.exists(os.path.join(ESC, st + '.csv'))}
        d = costs(frames)
        out = os.path.join(ESC, 'escAgeingCosts.csv')
        d.to_csv(out, index = False)
        print(f'-> {os.path.relpath(out, REPO)}  ({len(d)} rows from {", ".join(frames)})\n')
        for pin, lab in ((True, 'theta PINNED'), (False, 'theta CHOSEN')):
            sub = d[(d['θpinned'] == pin) & (d['scenario'] != 'baseline')]
            piv = sub.pivot_table(index = ['scenario', 'spec'], columns = 'ρ', values = ['dθ_t0', 'dτ_t0'],
                                  aggfunc = 'first')
            print('--- {}: changes at 2020 against the baseline of the same reading ---'.format(lab))
            print(piv.round(4).to_string(), '\n')
        return
    d = merge(atTheta(pd.read_csv(os.path.join(ESC, 'escShocks.csv')), args.theta),
              atTheta(pd.read_csv(os.path.join(ESC, 'escShocksCRRA.csv')), args.theta))
    out = os.path.join(ESC, 'escExperiments.csv')
    d.to_csv(out, index = False)
    print(f'-> {os.path.relpath(out, REPO)}  ({len(d)} rows)\n')

    order = ['baseline', 'mild', 'acute', 'frIncome', 'frLeisure', 'frVoting', 'frBoth',
             'frAll', 'France']
    for cx in sorted(d['commonX'].unique()):
        for method in sorted(d['method'].unique()):
            for spec in sorted(d.loc[d['commonX'] == cx, 'spec'].unique()):
                ds = d[(d['spec'] == spec) & (d['commonX'] == cx) & (d['method'] == method)]
                if ds.empty:
                    continue
                print('=' * 100)
                print('spec = {}, phi = 0.5, {}, method = {}   (theta_t0 = design in force at 2020; '
                      'tau/sr at 2020)'.format(spec, 'common X' if cx else 'vector X', method))
                print('=' * 100)
                for pin, lab in ((False, 'theta CHOSEN'), (True, 'theta PINNED at its exogenous value')):
                    sub = ds[ds['θpinned'] == pin]
                    if sub.empty:
                        continue
                    piv = sub.pivot_table(index = 'scenario', columns = 'ρ',
                                          values = ['θ_t0', 'τ_t0', 'sr_t0'], aggfunc = 'first')
                    piv = piv.reindex([s for s in order if s in piv.index])
                    print('\n--- {} ---'.format(lab))
                    print(piv.round(4).to_string())
                print()


if __name__ == '__main__':
    main()

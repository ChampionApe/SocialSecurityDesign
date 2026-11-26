r""" Stage (ii) of the paper pipeline, US/France/UK arm: run the counterfactuals the rich-OECD tables and
figure are read off. Nothing is calibrated here -- every experiment reads (beta, omega) from the sweep
csvs, so stage (i) (runCalibrationUS.py) must have run first.

Run:  .venv\Scripts\python.exe python\paper\runShocksUS.py               the MAIN part: whatever is missing
      ... --prepub                                                       the PRE-PUBLICATION part only
      ... --all                                                          both parts
      ... --only shocks                                                  one experiment
      ... --list                                                         what exists, what is missing
      ... --force                                                        re-run even where output exists
      ... --dry                                                          print the commands and exit

As in the Argentina arm, each entry declares WHAT the paper needs and which script produces it;
python/US/runShocksUS.py is the implementation and keeps its own CLI. This file records the settings the
published numbers were produced at.

THE θ VARIANT. Every command carries config.thetaFlags() (--theta config.US['theta'] --main: the U.S., France and
the UK at the variant's θ0, the standard results/ paths), and an output counts as present only where it holds rows
at that thetaSource (`_atTheta`). The data variants of the supplemental appendix's page C4,
config.US['oecd26']['thetaVariants'], run without --main into their own results/oecd26_theta_<V>/ (`oecd26Variants`).
Every command of a driver with a cost also carries config.spendingFlags() (--spendingTarget
config.US['esc']['spendingTarget']), and its outputs count only where they hold rows at that target (`_atTheta`,
config.spendingKeyed).

TWO PARTS. `part = 'main'` is what a routine rebuild runs: the exogenous-theta shocks in both variants and
the LOG leg of the endogenous-theta appendix, all cheap, and the cross-section. `part = 'prepub'` holds the
heavy runs whose results move only when the model does, run once before submission: the exact CRRA ESC leg
(the published method, config.US['esc']['exact']) and the stationary checks. The merges (`escAgeingCosts`,
`escExperiments`) belong to both. The two parts are selected on the command line, never mixed by default, so
that a development rebuild can never turn into an overnight solve.

ONE SCRIPT PRODUCES ALL THREE TABLES' WORTH OF ROWS. python/US/runShocksUS.py writes one long csv --
one row per (rho, family, scenario, effect) -- covering the theta, ageing and French-characteristic
families at every rho in config.US['ρTable']. That is why the whole set costs one entry here rather than
three: the families share a baseline solve per rho, and splitting them would re-solve it three times.

Cost of the main part: a few minutes for the exogenous shocks and the LOG ESC leg, ~2.5 h for the
cross-section from scratch (its readings in results/oecd26/; seconds when its rows are current), ~40 min per data
variant on its U.S. sweep and anchor.
"""
import os, sys, argparse, subprocess
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

for _s in (sys.stdout, sys.stderr):
    _s.reconfigure(encoding = 'utf-8')

THETA = C.thetaFlags()
# Every command of a driver that calibrates or reads a cost carries the spending target (config.spendingFlags); the
# exogenous-theta shocks have none.
SPENDINGFLAG = C.spendingFlags()
RHOFLAGS = ['--rho'] + [str(r) for r in C.US['ρTable']]
G = C.US['gridSettings']
GRIDFLAGS = ['--n', str(G['n']), '--ns', str(G['ns']),
             '--interpKind', G['interpKind'], '--smoothKnots', str(G['smoothKnots'])]
# The ESC leg runs under the headline calibration variant only. It is the most expensive stage in the
# pipeline, and unlike the exogenous-theta shocks -- where both variants are cheap and both are built --
# a second variant here would double a multi-hour run for an appendix the paper does not read twice.
ESCVARIANT = (['--commonX'] if C.US['commonX'] else []) + THETA + SPENDINGFLAG
ESC = C.US['esc']
# The paper's cost specification, for the CRRA legs and the stationary check (one spec each; the
# stationary check follows the paper's spec so sec:numerical's footnote quotes the paper's own design).
ESCSPEC = ['--spec', ESC['spec'], '--phi', str(ESC['phi'])]
# The LOG leg runs the comparison arm too: it costs seconds, the supplemental appendix's second cost reads
# its rows, and runESC.stagePath and stageCountry OVERWRITE escPath.csv and escCountry.csv with the specs
# of the run rather than merging (stageCalib and stageShocks merge), so a single-spec LOG run would drop
# the other spec's path and country rows.
ESCSPEC_LOG = (['--spec', ESC['spec']] + ([ESC['comparisonSpec']] if ESC.get('comparisonSpec') else [])
               + ['--phi', str(ESC['phi'])])
ESCRHO = [str(r) for r in ESC['ρTable'] if r != C.US['ρAnchor']]
# The exact recursion's design layer, forwarded to every runESCcrra.py command (config.US['esc']).
ESCLAYER = ['--designRule', ESC['designRule'], '--Ma', str(ESC['Ma'])]
# The cost variant of section 7 (config.US['esc']['costVariants']): the Pi-inclusive cost, in its tagged csvs.
CV = ESC['costVariants']
CVSPEC = ['--spec', CV['spec'], '--phi', str(ESC['phi']), '--tag', CV['tag']]
OECD = C.US['oecd26']
TAUTHETA = C.US['tauTheta']


def _read(path):
    return pd.read_csv(path, dtype = {'inputsKey': str, 'anchorKey': str}) if os.path.exists(path) else None


def _atTheta(path, variant = None):
    """ Output check: the csv exists AND holds rows at thetaSource `variant` (config.US['theta'] by default) and, for
    a csv of a driver with a cost (config.spendingKeyed), at config.US['esc']['spendingTarget']. The file's mere
    existence says nothing: a csv written before the variant holds index rows only, one before the target gross
    rows only. """
    df = _read(path)
    if df is None:
        return False
    df = C.atTheta(df, variant)
    return len(C.atSpending(df) if C.spendingKeyed(path) else df) > 0


def _hasExactRows(path, ρs = None, spec = None):
    """ Output check for the exact CRRA leg: the csv exists AND carries method = 'exact' rows at every rho
    the tables print, under the headline variant and the θ variant (and `spec`, if given). The file's mere
    existence says nothing -- the path iteration's rows live in the same csv under method = 'path'. """
    df = _read(path)
    if df is None or 'method' not in df.columns:
        return False
    df = C.atSpending(C.atTheta(df))
    df = df[df['method'].fillna('path') == 'exact']
    if spec is not None:
        df = df[df['spec'] == spec]
    if 'commonX' in df.columns:
        df = df[df['commonX'].fillna(False).astype(bool) == bool(C.US['commonX'])]
    have = set(round(float(r), 6) for r in df['ρ'])
    return all(round(float(r), 6) in have for r in (ρs or [float(r) for r in ESCRHO]))


def _usScript(name):
    return [C.PYTHON, os.path.join(C.USDIR, name)]


def _variantDir(V):
    return os.path.join(C.RESULTS, 'oecd26_theta_' + V)


EXPERIMENTS = {
    # Every family, every rho in the table grid, both readings (full effect and economic-equilibrium-only).
    # The VECTOR-X variant: with config.US['commonX'] = True this is the appendix's robustness twin, and
    # 'shocksCommonX' below carries the headline tables and figure. Both are built either way -- they cost
    # minutes each and the pair is a check on itself (see 'shocksCommonX').
    'shocks': {
        'part':    'main',
        'script':  'runShocksUS.py',
        'args':    RHOFLAGS + GRIDFLAGS + THETA,
        'outputs': lambda: [os.path.join(C.SHOCKDIR, 'US_shocks.csv')],
        'requires': lambda: [C.usSweepCsv('US', False)],
        'note':    'theta / ageing / French, all rho, both effects',
    },
    # The same set under the common-X calibration variant -- the HEADLINE one. Cheap and worth
    # having as a check as much as an output: theta, ageing and voting must come back IDENTICAL to the
    # vector-X run (they touch neither eta nor X, and beta/omega/h agree across variants), while income
    # distribution and leisure must differ -- those two are defined through eta and X, whose meaning is
    # exactly what the variant changes. A commonX run that matched on all seven would mean the variant
    # was not being applied.
    'shocksCommonX': {
        'part':    'main',
        'script':  'runShocksUS.py',
        'args':    RHOFLAGS + GRIDFLAGS + ['--commonX'] + THETA,
        'outputs': lambda: [os.path.join(C.SHOCKDIR, 'US_shocksCommonX.csv')],
        'requires': lambda: [C.usSweepCsv('US', True)],
        'note':    'the same set under common X',
    },
    # The US counterfactuals at other Frisch elasticities, rho = 1, common X: the US and France recalibrated at
    # each xi inside the script (python/US/runShocksXi.py, whose first xi reproduces US_shocksCommonX.csv).
    'shocksXi': {
        'part':    'main',
        'script':  'runShocksXi.py',
        'args':    THETA,
        'outputs': lambda: [os.path.join(C.SHOCKDIR, 'US_shocksXiCommonX.csv')],
        'requires': lambda: [os.path.join(C.SHOCKDIR, 'US_shocksCommonX.csv')],
        'note':    'the US counterfactuals at xi = 0.2, 0.3, 0.4 (LOG, common X)',
    },
    # --- Endogenous system characteristics (sec:esc), the LOG leg (rho = 1), the paper's spec and the
    # comparison spec (ESCSPEC_LOG). Requires the LOG cost calibrations of both from stage (i)
    # (runCalibrationUS.py's escMissing step). All cheap.
    # The equilibrium design path under the calibrated cost: the 2020 -> 2110 drift the text quotes
    # (results/esc/escPath.csv).
    'escPath': {
        'part':    'main',
        'script':  'runESC.py',
        'args':    ['--stage', 'path'] + ESCSPEC_LOG + ESCVARIANT,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escPath.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibration.csv')],
        'note':    'endogenous-theta design path, LOG (rho = 1)',
    },
    # The counterfactuals on the LEADED-choice model at the calibrated cost, each run twice (theta
    # pinned at the calibrated design vs chosen).
    'escShocks': {
        'part':    'main',
        'script':  'runESC.py',
        'args':    ['--stage', 'shocks'] + ESCSPEC_LOG + ESCVARIANT,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escShocks.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibration.csv')],
        'note':    'endogenous-theta counterfactuals, LOG (rho = 1)',
    },
    # The cost variant (section 7, A8): acute ageing and its baseline under the Pi-inclusive cost, LOG, at its own
    # calibrated epsilon (escCalibrationCosts.csv, stage (i)).
    'escCosts': {
        'part':    'main',
        'script':  'runESC.py',
        'args':    ['--stage', 'shocks'] + CVSPEC + ['--scenarios'] + CV['scenarios'] + ESCVARIANT,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escShocksCosts.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibrationCosts.csv')],
        'note':    'the Pi-inclusive cost: acute ageing against its baseline, LOG (rho = 1)',
    },
    # The UK and France under the US parameter and under their own (Tables/US_ESC_Country.tex), each at its θ0. The
    # csv also holds rows of the UK at US income groups ('UKUS'), which the stage no longer writes and nothing prints.
    'escCountry': {
        'part':    'main',
        'script':  'runESC.py',
        'args':    ['--stage', 'country'] + ESCSPEC_LOG + ESCVARIANT,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escCountry.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibration.csv')],
        'note':    'the UK and France under the ESC cost, LOG (rho = 1)',
    },
    # The U.S. tax of 2020 against a design pinned on the whole path, without the cost and at the calibrated lambda,
    # LOG (python/US/runTauTheta.py; config.US['tauTheta']; US_ESC_TauTheta). Its rows are keyed on the run's
    # calibration (beta, omega, lambda), so the child decides what to recompute: always run (seconds).
    'escTauTheta': {
        'part':    'main',
        'always':  True,
        'script':  'runTauTheta.py',
        'args':    (['--grid'] + [repr(float(x)) for x in TAUTHETA['grid']] + ['--refine', repr(TAUTHETA['refine'])]
                    + ESCVARIANT),
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escTauTheta.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibration.csv')],
        'note':    'the U.S. tax against a pinned design, with and without the cost, LOG (rho = 1)',
    },
    # The same path for France and the UK under the US-calibrated cost and their own (escCountry.csv): the
    # supplemental appendix's paths by country (runESC.stagePathHost).
    'escPathHost': {
        'part':    'main',
        'script':  'runESC.py',
        'args':    ['--stage', 'pathHost'] + ESCSPEC + ESCVARIANT,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escPathHost.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibration.csv'), os.path.join(C.ESCDIR, 'escCountry.csv')],
        'note':    'endogenous-theta design path of France and the UK, LOG (rho = 1)',
    },
    # The cross-country test (python/US/runOECD26.py, results/oecd26/): the U.S. at the OECD inputs (anchor), every
    # in-sample country of the θ variant at its omega and cost parameter (layer 1) and with its own omega (layer 2),
    # LOG; one command per variant of config.US['oecd26']['variants'] (costs, turnout, layer 1's omega, stages) and per
    # retAge of it, in its order: the paper's specification, the supplemental appendix's (page C4) population variants
    # and readings. Its rows are keyed on their inputs, so the child decides what to recompute: always run.
    'oecd26': {
        'part':    'main',
        'always':  True,
        'script':  'runOECD26.py',
        'argSets': lambda: [['--spec'] + v['spec'] + ['--betaRule'] + v['betaRule'] + ['--turnout'] + v['turnout']
                            + ['--stage'] + v['stage'] + ['--inSample', '--omegaRule', v['omegaRule'], '--retAge', r]
                            + THETA + SPENDINGFLAG
                            for v in OECD['variants'] for r in v['retAge']],
        'outputs': lambda: [os.path.join(C.RESULTS, 'oecd26', f + '.csv') for f in ('anchor', 'layer1', 'layer2')],
        'requires': lambda: [C.usSweepCsv('US', True)],
        'note':    'the cross-country test, LOG (rho = 1)',
    },
    # Page C4's data variants (config.US['oecd26']['thetaVariants']), each at its own U.S. sweep at rho = 1
    # (calibrateRhoGrid.py --theta V, which resumes at its θ0) and the paper's specification of the cross-section
    # (turnout, layer 1's omega, retAge), into results/oecd26_theta_<V>/ (OECD26_Variants). Always run, as 'oecd26'.
    'oecd26Variants': {
        'part':    'main',
        'always':  True,
        'commands': lambda: sum(([_usScript('calibrateRhoGrid.py') + ['--commonX', '--lo', '1', '--hi', '1',
                                                                      '--theta', V],
                                  _usScript('runOECD26.py') + ['--theta', V, '--inSample', '--spec', OECD['spec'],
                                                               '--betaRule', OECD['betaRule'],
                                                               '--turnout', OECD['turnout'],
                                                               '--stage', 'anchor', 'layer2', 'layer1',
                                                               '--omegaRule', OECD['omegaRule'],
                                                               '--retAge', OECD['retAge']] + SPENDINGFLAG]
                                 for V in OECD['thetaVariants']), []),
        'outputs': lambda: [os.path.join(_variantDir(V), f + '.csv') for V in OECD['thetaVariants']
                            for f in ('anchor', 'layer1', 'layer2')],
        'atTheta': lambda p: os.path.basename(os.path.dirname(p))[len('oecd26_theta_'):],
        'note':    'the cross-country test under the robustness measures of theta',
    },
    # --- PRE-PUBLICATION. The CRRA leg by the EXACT 2-D recursion (LeadedCRRA2D), the published method:
    # the design path and the counterfactuals at the exact-calibrated p from stage (i) --prepub. One full
    # recursion per chosen-design reading, ~6 min each at ns = 150; ~1 h for the two rho. Rows carry
    # method = 'exact' and sit beside the path iteration's rows (method = 'path'), which no paper output
    # reads (config.US['esc']['exact']); the output check therefore looks for exact rows, not the file.
    'escShocksCRRA': {
        'part':    'prepub',
        'script':  'runESCcrra.py',
        'args':    (['--exact', '--stage', 'path', 'shocks', '--rho'] + ESCRHO + ESCSPEC
                    + ['--ns', str(ESC['ns2D']), '--nsScan', str(ESC['nsScan']),
                       '--nCand2D', str(ESC['nCand2D'])] + ESCVARIANT + ESCLAYER),
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escShocksCRRA.csv'),
                            os.path.join(C.ESCDIR, 'escPathCRRA.csv')],
        'complete': lambda: all(_hasExactRows(os.path.join(C.ESCDIR, f), spec = ESC['spec'])
                                for f in ('escShocksCRRA.csv', 'escPathCRRA.csv')),
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibrationCRRA.csv')],
        'note':    'endogenous-theta path and counterfactuals, CRRA, EXACT 2-D recursion',
    },
    # --- PRE-PUBLICATION. Stationary vs date-specific policy functions, the check behind sec:numerical's
    # literature paragraph: the exact date-specific recursion against a stationary policy function computed
    # at each date's own nu and walked along the demographic path (plus the steady-state comparison and one
    # long-run function). No paper table reads these; the footnote quotes the maximal gaps. They move only
    # when the model or the calibration does, hence prepub. The exogenous-theta check under CRRA (LOG is
    # exact by the decoupling) and the endogenous-theta check under LOG and the exact 2-D CRRA recursion
    # (~1 h at ns2D, one recursion per distinct nu), under the paper's cost spec (ESCSPEC). The Argentina
    # counterpart is runShocks.py --prepub.
    'stationary': {
        'part':    'prepub',
        'script':  'stationaryApprox.py',
        'args':    (['--rho'] + [str(r) for r in C.US['ρStationary']]
                    + ['--n', str(G['n']), '--ns', str(G['ns']), '--smoothKnots', str(G['smoothKnots'])]
                    + ([] if C.US['commonX'] else ['--vectorX']) + THETA),
        'outputs': lambda: [os.path.join(C.NUMDIR, 'US_stationaryApprox.csv')],
        'requires': lambda: [C.usSweepCsv('US', C.US['commonX'])],
        'note':    'stationary vs date-specific tax policies, CRRA (~2 min)',
    },
    'stationaryESC': {
        'part':    'prepub',
        'script':  'stationaryApproxESC.py',
        'args':    (['--rho'] + [str(r) for r in ESC['ρTable']] + ESCSPEC
                    + ['--ns', str(ESC['ns2D']), '--nCand2D', str(ESC['nCand2D'])] + ESCVARIANT),
        'outputs': lambda: [os.path.join(C.NUMDIR, 'US_ESC_stationaryApprox.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escCalibrationCRRA.csv')],
        'note':    'stationary vs date-specific design policies, LOG and exact CRRA (~1 h)',
    },
    # The cost variant's ageing table (section 7, A8): the paper's csvs and the Costs-tagged one merged, at the θ
    # variant only (collectESCexperiments.py --theta).
    'escAgeingCosts': {
        'part':    'both',
        'always':  True,
        'script':  'collectESCexperiments.py',
        'args':    ['--costs', '--theta', C.US['theta']] + SPENDINGFLAG,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escAgeingCosts.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, 'escShocksCosts.csv')],
        'note':    'merge the ageing counterfactuals under the three costs',
    },
    # The merge stage (iii) reads. Listed last so a run rebuilds it after the producers, in BOTH parts,
    # and NEVER skipped ('always'): it costs a second, and a merge left over from before a producer ran
    # is exactly the stale-but-present output that blocked the build on 2026-09-11.
    'escExperiments': {
        'part':    'both',
        'always':  True,
        'script':  'collectESCexperiments.py',
        'args':    ['--theta', C.US['theta']] + SPENDINGFLAG,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escExperiments.csv')],
        'note':    'merge the LOG and CRRA legs into one long table',
    },
    # The ESC counterfactuals at other Frisch elasticities (python/US/runESCxi.py; the site's xi blocks of
    # ESC_Counterfactuals and ESC_Paths): the cost recalibrated at each xi, the workbook xi first as a self-check.
    # After escExperiments, whose merge the self-check reads.
    'escXi': {
        'part':    'main',
        'script':  'runESCxi.py',
        'args':    THETA + SPENDINGFLAG,
        'outputs': lambda: [os.path.join(C.ESCDIR, 'escXiRobustness.csv')],
        'requires': lambda: [os.path.join(C.ESCDIR, f) for f in ('escCalibration.csv', 'escPath.csv', 'escExperiments.csv')],
        'note':    'endogenous-theta counterfactuals at xi = 0.2, 0.3, 0.4, LOG (rho = 1)',
    },
}


def selected(prepub = False, all_ = False):
    """ The entry names of the requested part(s), in declaration order. """
    parts = {'main', 'prepub', 'both'} if all_ else ({'prepub', 'both'} if prepub else {'main', 'both'})
    return [n for n, e in EXPERIMENTS.items() if e['part'] in parts]


def status(name):
    """ (script exists, outputs present at the θ variant, outputs lacking). An 'always' entry is re-run whatever is
    on disk (`run`), which says nothing about what is on disk: its present outputs still count as present. """
    e = EXPERIMENTS[name]
    script = os.path.join(C.USDIR, e['script']) if 'script' in e else None
    outs = e['outputs']()
    variantOf = e.get('atTheta', lambda p: None)
    have = [p for p in outs if _atTheta(p, variantOf(p))]
    lack = [p for p in outs if p not in have]
    if not lack and 'complete' in e and not e['complete']():
        lack = [p + ' (no exact rows)' for p in outs]
    return script is None or os.path.exists(script), have, lack


def run(name, force = False):
    """ Does `name` run: anything lacking, --force, or an 'always' entry. """
    _, _, lack = status(name)
    return bool(lack) or force or bool(EXPERIMENTS[name].get('always'))


def commands(name):
    """ The entry's commands: its 'commands' (a callable of full commands), else one per argument list of its
    'argSets' (a callable), else the one of its 'args', each on its 'script'. """
    e = EXPERIMENTS[name]
    if 'commands' in e:
        return e['commands']()
    sets = e['argSets']() if 'argSets' in e else [e['args']]
    return [[C.PYTHON, os.path.join(C.USDIR, e['script'])] + list(a) for a in sets]


def main():
    p = argparse.ArgumentParser(description = __doc__.split('\n')[1])
    p.add_argument('--only', nargs = '+', choices = list(EXPERIMENTS), default = None,
                   help = 'a subset, regardless of part')
    p.add_argument('--prepub', action = 'store_true', help = 'the pre-publication part instead of the main one')
    p.add_argument('--all', dest = 'all_', action = 'store_true', help = 'both parts')
    p.add_argument('--list', action = 'store_true', help = 'report what exists and exit')
    p.add_argument('--force', action = 'store_true', help = 're-run even where output already exists')
    p.add_argument('--dry', action = 'store_true', help = 'print commands and exit')
    p.add_argument('--keepGoing', action = 'store_true',
                   help = 'after a failed entry, go on with the others (each still skips while its inputs lack rows '
                          'at the θ variant); exit 1 at the end if any failed')
    p.add_argument('--exclude', nargs = '+', choices = list(EXPERIMENTS), default = [],
                   help = 'leave these entries out of the selection (e.g. oecd26Variants, run as its own step)')
    a = p.parse_args()
    names = [n for n in (a.only if a.only else selected(a.prepub, a.all_)) if n not in a.exclude]

    if a.list or a.dry:
        print('θ variant {} (config.US[\'theta\']): an output is "had" where it holds rows at it'.format(C.US['theta']))
        for name in names:
            ok, have, lack = status(name)
            n = len(EXPERIMENTS[name]['outputs']())
            print('{:<7} {:<18} {:<66} script:{}  have {}/{}{}{}'.format(
                EXPERIMENTS[name]['part'], name, EXPERIMENTS[name]['note'], 'yes' if ok else 'MISSING',
                len(have), n, '  (no exact rows yet)' if any('no exact' in l for l in lack) else '',
                '  (always re-run)' if EXPERIMENTS[name].get('always') else ''))
            if a.dry and ok:
                for cmd in commands(name):
                    print('    ' + ' '.join(cmd))
        return

    failed = []
    for name in names:
        ok, have, lack = status(name)
        if not ok:
            print('SKIP {}: {} does not exist yet.'.format(name, EXPERIMENTS[name]['script']))
            continue
        if not run(name, a.force):
            print('SKIP {}: all {} output(s) present.'.format(name, len(have)))
            continue
        # an input counts once it holds rows at the θ variant: a stale csv of another θ is not one
        need = [p for p in EXPERIMENTS[name].get('requires', lambda: [])() if not _atTheta(p)]
        if need:
            print('SKIP {}: needs {} at θ variant {} first.'.format(
                name, ', '.join(os.path.relpath(p, C.REPO) for p in need), C.US['theta']))
            if a.keepGoing:
                failed.append(name + ' (inputs)')
            continue
        for cmd in commands(name):
            print('\n' + '='*94 + '\n{} [{}]: {}\n  {}\n'.format(name, EXPERIMENTS[name]['part'],
                                                                 EXPERIMENTS[name]['note'], ' '.join(cmd))
                  + '='*94)
            r = subprocess.run(cmd, cwd = C.REPO)
            if r.returncode:
                if not a.keepGoing:
                    raise SystemExit('{} exited {}'.format(os.path.basename(cmd[1]), r.returncode))
                print('FAILED {}: {} exited {}'.format(name, os.path.basename(cmd[1]), r.returncode))
                failed.append(name)
                break
    if failed:
        raise SystemExit('failed or skipped for lack of inputs: {}'.format(', '.join(failed)))


if __name__ == '__main__':
    main()

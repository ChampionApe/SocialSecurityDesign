r""" Replicate the paper: every step that produces its results, in order of dependence, at the paper's settings.

    python python/paper/replicate.py --dry            # print every command and exit
    python python/paper/replicate.py --list           # what exists, step by step
    python python/paper/replicate.py --quick          # step 'build' only: every exhibit from the committed results
    python python/paper/replicate.py                  # every step below but 'data' (cold: about a day)
    python python/paper/replicate.py --only c5 build  |  --from shocksUS  |  --prepub  |  --data  |  --force

Steps, with the cold run times measured on the authors' machine where they were:
  data        stage (0), only with --data (network): Penn World Table targets, figure 1's cross-section, the 29
              countries' inputs, the sample rule, pension tiers (PowerShell, Windows only), top-ups, the measure θ
  calibARG    stage (i), Argentina: runCalibration.py --commonX                                 ~3 h
  calibUS     stage (i), the U.S., France and the UK: runCalibrationUS.py --commonX             ~40 min
  shocksARG   stage (ii), Argentina: runShocks.py --commonX                                     ~1 h
  shocksUS    stage (ii), the U.S. arm, the endogenous design and the cross-section             ~3 h, + 1.5-6 h for
              (runShocksUS.py)                                                                  page C4's data variants
  argEps0     the ε = 0 Argentine calibration (calibARG's sweep at ε = 0) and its reform
  argDiag     the reform's diagnostics at the headline calibration (reformDiagnostics.py)
  c5          the smaller U.S. numbers the text quotes (python/US/c5Numbers.py)               designCRRA ~1.5 h per ρ
  horizonUS   appendix D's checks on the U.S. (python/US/runHorizon.py)                         ~40 min, + ~4 h design
  horizonARG  appendix D's checks on Argentina (python/InformalSavings/runHorizon.py)           ~45 min
  build       stage (iii): every table and figure into results/paper (build.py --no-copy)       seconds
--prepub adds the pre-publication parts of calibUS, shocksARG and shocksUS (about 9 hours), which no exhibit of
the paper reads. The stage drivers skip work whose output exists and --force makes them redo it; every other step
is skipped when all its outputs exist, unless --force. Commands run from the repository root under this
interpreter, with PYTHONUTF8=1.
"""
import os, sys, glob, time, argparse, datetime, subprocess
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config as C

for _s in (sys.stdout, sys.stderr):              # commands carry '--nι'
    _s.reconfigure(encoding = 'utf-8')

G = C.ARG['gridSettings']
GRIDFLAGS = ['--interpKind', G['interpKind'], '--smoothKnots', str(G['smoothKnots']),
             '--nι', str(G['nι']), '--ns', str(G['ns'])]
THETA = C.thetaFlags()                            # --theta C --main: the paper's measure of the design
SPENDING = C.spendingFlags()                      # --spendingTarget: the drivers that read a cost (c5, horizonUS)
ARGRULE = ['--rule', C.ARG['rule'], '--refType', str(C.ARG['refType'])]


def py(*parts):
    return os.path.join(C.REPO, 'python', *parts)


def shock(name):
    return os.path.join(C.SHOCKDIR, name)


def num(name):
    return os.path.join(C.NUMDIR, name)


class Cmd:
    """ One command. `outputs` (paths or globs) let the step skip it when they all exist and, `atSpending`, each
    holds rows at the paper's spending target (config.atSpending); None: the script decides itself what to redo (a
    stage driver or a resuming sweep), and `force` is what --force appends to it. """
    def __init__(self, args, outputs = None, force = (), atSpending = False):
        self.args, self.outputs, self.force, self.atSpending = list(args), outputs, list(force), atSpending

    def done(self):
        if self.outputs is None or not all(glob.glob(p) for p in self.outputs):
            return False
        return not self.atSpending or all(len(C.atSpending(pd.read_csv(f))) for p in self.outputs for f in glob.glob(p))


def steps(a):
    part = ['--all'] if a.prepub else []
    keep = ['--keepGoing'] if a.keepGoing else []
    eps0Csv = os.path.join(C.CALIBDIR, 'informalSavings_rhoGridCommonXEpsPre0.csv')
    eps0Dir = os.path.join(C.CALIBDIR, 'instancesCommonXEpsPre0')
    reform0 = C.argShockTemplate('reform', commonX = True).replace('.csv', '_epsPre0.csv')
    ee0 = C.argShockTemplate('ee', commonX = True).replace('.csv', '_epsPre0.csv')
    ρs = C.ARG['ρGrid']
    c5out = {'a1': 'C5_criticalTurnout.csv', 'a9': 'C5_ranking.csv', 'loss': 'C5_marginalLoss.csv',
             'lambda': 'C5_lambdaElasticity.csv', 'index': 'C5_indexUnits.csv', 'design': 'C5_designDecomp.csv',
             'a13': 'C5_taxPolicyState.csv', 'strategic': 'C5_strategicTerm.csv'}
    return [
        ('data', 'stage (0): the inputs in data/ (network)', [
            Cmd([py('paper', s)]) for s in ('dataTargets.py', 'oecdFigure1.py', 'oecd29Data.py', 'oecd26Sample.py',
                                            'oecd26PensionTiers.py', 'oecd26Topups.py', 'oecd26Theta.py')]),
        ('calibARG', 'stage (i), Argentina', [
            Cmd([py('paper', 'runCalibration.py'), '--commonX'], force = ['--force'])]),
        ('calibUS', 'stage (i), the U.S., France and the UK', [
            Cmd([py('paper', 'runCalibrationUS.py'), '--commonX'] + part + keep, force = ['--force'])]),
        ('shocksARG', 'stage (ii), Argentina', [
            Cmd([py('paper', 'runShocks.py'), '--commonX'] + part, force = ['--force'])]),
        ('shocksUS', 'stage (ii), the U.S. arm, the endogenous design and the cross-section', [
            Cmd([py('paper', 'runShocksUS.py')] + part + keep, force = ['--force'])]),
        ('argEps0', 'the ε = 0 Argentine calibration (it resumes) and its reform, every ρ', [
            Cmd([py('InformalSavings', 'calibrateRhoGrid.py'), '--lo', '0.5', '--hi', '2.0', '--step', '0.1',
                 '--anchor', '1.0', '--commonX', '--epsPre0', '--x0', '0.693156', '1.663414', '0.292216', '0.411944']
                + GRIDFLAGS, force = ['--force']),
            Cmd([py('InformalSavings', 'shockUniversal.py'), '--epsPre0'] + ARGRULE + GRIDFLAGS
                + ['--pkldir', eps0Dir, '--csv', eps0Csv, '--out', reform0],
                outputs = [shock(reform0.format(ρ = ρ)) for ρ in ρs]),
            Cmd([py('InformalSavings', 'shockEEOnly.py')] + ARGRULE + GRIDFLAGS
                + ['--pkldir', eps0Dir, '--csv', eps0Csv, '--baseCsv', reform0, '--out', ee0],
                outputs = [shock(ee0.format(ρ = ρ)) for ρ in ρs])]),
        ('argDiag', "the reform's diagnostics at the headline calibration", [
            Cmd([py('InformalSavings', 'reformDiagnostics.py'), '--rho', '0.5', '0.6', '1.0', '2.0',
                 '--pkldir', C.argInstanceDir(True), '--refine', '60', '--refineRho', '0.5', '0.6'],
                outputs = [shock('ARG_reform{}_commonX.csv'.format(k)) for k in ('FOC', 'B5', 'HouseholdSavings', 'Refine')])]),
        ('c5', 'the smaller U.S. numbers the text quotes', [
            Cmd([py('US', 'c5Numbers.py'), c] + THETA + SPENDING, outputs = [num(f)], atSpending = True)
            for c, f in c5out.items()] + [
            Cmd([py('US', 'c5Numbers.py'), 'designCRRA', '--rho', r] + THETA + SPENDING,
                outputs = [num('C5_designDecompCRRA_rho{:g}.csv'.format(float(r)))], atSpending = True)
            for r in ('0.5', '2.0')]),
        ('horizonUS', "appendix D's checks, the U.S.", [
            Cmd([py('US', 'runHorizon.py'), '--stage', 'residuals', 'jumps', 'scheme'] + THETA,
                outputs = [num(f) for f in ('US_residuals.csv', 'US_jumps.csv', 'US_derivativeScheme.csv')]),
            Cmd([py('US', 'runHorizon.py'), '--stage', 'log', 'horizon'] + THETA + SPENDING,
                outputs = [num('US_horizon.csv')], atSpending = True),
            Cmd([py('US', 'runHorizon.py'), '--stage', 'design'] + THETA + SPENDING, outputs = [num('US_horizonDesign.csv')],
                atSpending = True)]),
        ('horizonARG', "appendix D's checks, Argentina", [
            Cmd([py('InformalSavings', 'runHorizon.py')],
                outputs = [num('ARG_horizon_commonX.csv'), num('ARG_residuals_commonX.csv')])]),
        ('build', 'stage (iii): every table and figure into results/paper', [
            Cmd([py('paper', 'build.py'), '--no-copy'])]),
    ]


def show(args):
    rel = lambda s: ('.' + s[len(C.REPO):]).replace('\\', '/') if s.startswith(C.REPO) else s
    return ' '.join(rel(x) if ' ' not in x else '"{}"'.format(rel(x)) for x in ['python'] + args)


def stamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def run(args):
    env = dict(os.environ, PYTHONUTF8 = '1')
    return subprocess.run([sys.executable] + args, cwd = C.REPO, env = env).returncode


def main():
    p = argparse.ArgumentParser(description = __doc__.split('\n')[1])
    p.add_argument('--dry', action = 'store_true', help = 'print every command and exit')
    p.add_argument('--list', action = 'store_true', help = 'report what exists, step by step, and exit')
    p.add_argument('--quick', action = 'store_true', help = "only the step 'build'")
    p.add_argument('--only', nargs = '+', default = None, help = 'run only these steps')
    p.add_argument('--from', dest = 'start', default = None, help = 'start at this step')
    p.add_argument('--data', action = 'store_true', help = 'include stage (0), which reaches the network')
    p.add_argument('--prepub', action = 'store_true', help = 'add the pre-publication parts of the stage drivers')
    p.add_argument('--force', action = 'store_true', help = 'redo work whose output exists')
    p.add_argument('--keepGoing', action = 'store_true', help = 'go on after a failed step')
    a = p.parse_args()
    allSteps = steps(a)
    names = [n for n, _, _ in allSteps]
    for n in (a.only or []) + ([a.start] if a.start else []):
        if n not in names:
            p.error('unknown step {!r}; the steps are {}'.format(n, ', '.join(names)))
    if a.quick:
        chosen = ['build']
    elif a.only:
        chosen = a.only
    else:
        chosen = names[names.index(a.start):] if a.start else names
        if not a.data and 'data' not in (a.only or []):
            chosen = [n for n in chosen if n != 'data']
    todo = [s for s in allSteps if s[0] in chosen]

    if a.dry or a.list:
        for name, note, cmds in todo:
            print('[{}] {}'.format(name, note))
            for c in cmds:
                args = c.args + (c.force if a.force else [])
                if a.list:
                    state = ('done' if c.done() else 'missing') if c.outputs is not None else 'the script decides'
                    print('  {:<18} {}'.format(state, show(args)))
                else:
                    print('  ' + show(args))
        if a.list:
            print('\nStage drivers report their own state: runShocks.py/runShocksUS.py --list, runCalibration*.py --dry,'
                  ' build.py --list.')
        return 0

    failed, t00 = [], time.time()
    for name, note, cmds in todo:
        t0 = time.time()
        print('{}  [{}] start: {}'.format(stamp(), name, note), flush = True)
        status = 'OK'
        for c in cmds:
            if c.done() and not a.force:
                print('{}  [{}] skip (outputs exist): {}'.format(stamp(), name, show(c.args)), flush = True)
                continue
            args = c.args + (c.force if a.force else [])
            print('{}  [{}] run: {}'.format(stamp(), name, show(args)), flush = True)
            if run(args) != 0:
                status = 'FAILED'
                break
        print('{}  [{}] {} after {:.1f} min'.format(stamp(), name, status, (time.time() - t0) / 60), flush = True)
        if status != 'OK':
            failed.append(name)
            if not a.keepGoing:
                break
    print('{}  done in {:.1f} h; failed: {}'.format(stamp(), (time.time() - t00) / 3600, ', '.join(failed) or 'none'))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())

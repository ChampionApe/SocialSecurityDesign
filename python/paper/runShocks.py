r""" Stage (ii) of the paper pipeline: run the experiments whose output the paper's tables and figures
are read off. Nothing is calibrated here -- every experiment starts from a pickled calibrated instance,
so stage (i) must have run first.

Run:  .venv\Scripts\python.exe python\paper\runShocks.py                 run whatever is missing
      ... --commonX                                                     also the common-X variant
      ... --only universal                                              one experiment
      ... --list                                                        what exists, what is missing
      ... --force                                                       re-run even where output exists,
                                                                        and re-do work inside a script
                                                                        that resumes from its own csv
      ... --dry                                                         print the commands and exit

Each entry of EXPERIMENTS declares WHAT the paper needs and which script produces it. The scripts are
the implementations and keep their own CLIs; this file is the record of the settings the published
numbers were produced at, so that reproducing them is a run rather than an archaeology.

TWO CALIBRATION VARIANTS (config.ARG['commonX'] names the headline): every entry runs per variant,
reading that variant's instances (config.argInstanceDir) and writing csvs with its suffix
(config.argVariantTag). The scripts take the variant purely through --pkldir/--out/--csv/--baseCsv.

Cost, cold, on the 16-point rho grid: `universal` ~2.5 h (a full backward PEE recursion per rho),
`flat` ~10 min (anchor only), `eeOnly` ~minutes (taxes are exogenous, so it is one EE solve per rho and
no policy recursion), `epsThetaGrid` ~15 min (378 points at ~2.4 s; the older ~5 min was measured under
the pre-2026-08-24 calibration). Every entry is skipped when its output already
exists, so a re-run with nothing to do costs seconds.
"""
import os, sys, argparse, subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

for _s in (sys.stdout, sys.stderr):             # --dry prints GRIDFLAGS, which carry '--nι'
    _s.reconfigure(encoding = 'utf-8')

G = C.ARG['gridSettings']
GRIDFLAGS = ['--interpKind', G['interpKind'], '--smoothKnots', str(G['smoothKnots']),
             '--nι', str(G['nι']), '--ns', str(G['ns'])]


def _rhoFiles(template, ρGrid):
    return {ρ: os.path.join(C.SHOCKDIR, template.format(ρ = ρ)) for ρ in ρGrid}


def _variantFlags(cx):
    """ The instance directory and sweep csv of one variant, as the scripts' own flags. """
    return ['--pkldir', C.argInstanceDir(cx), '--csv', C.argSweepCsv(cx)]


EXPERIMENTS = {
    # The headline reform: b^0 = b^1 at t0, taxes re-optimised politically. Tables funcOfRho and the
    # "Full effect" row of ArgentinaUniversal, and figure ARG_CRRA_LOG.
    'universal': {
        'script':  'shockUniversal.py',
        'args':    lambda cx: (['--rule', C.ARG['rule'], '--refType', str(C.ARG['refType'])] + GRIDFLAGS
                               + _variantFlags(cx) + ['--out', C.argShockTemplate('reform', commonX = cx)]),
        'outputs': lambda cx: _rhoFiles(C.argShockTemplate('reform', commonX = cx), C.ARG['ρGrid']),
        'note':    'full effect, all rho',
    },
    # The other reading of "universal" (eps = 1-theta, the non-contributive component only). It falls on
    # the OTHER side of the status quo, so the two bracket the reform rather than differing in degree.
    # Anchor only -- extending it to the whole grid is an open item (InformalSavings/README.md).
    'flat': {
        'script':  'shockUniversal.py',
        'args':    lambda cx: (['--rule', 'flat'] + GRIDFLAGS + _variantFlags(cx)
                               + ['--out', C.argShockTemplate('flat', commonX = cx)]),
        'rho':     [C.ARG['ρAnchor']],
        'outputs': lambda cx: _rhoFiles(C.argShockTemplate('flat', commonX = cx), [C.ARG['ρAnchor']]),
        'note':    'bracketing reading, anchor only',
    },
    # The same reform with taxes held at the baseline path: the "Economic Equilibrium" row of
    # ArgentinaUniversal, which is what separates the pure equilibrium response from the policy response.
    'eeOnly': {
        'script':  'shockEEOnly.py',
        'args':    lambda cx: (['--rule', C.ARG['rule'], '--refType', str(C.ARG['refType'])] + GRIDFLAGS
                               + _variantFlags(cx)
                               + ['--baseCsv', C.argShockTemplate('reform', commonX = cx),
                                  '--out', C.argShockTemplate('ee', commonX = cx)]),
        'outputs': lambda cx: _rhoFiles(C.argShockTemplate('ee', commonX = cx), C.ARG['ρGrid']),
        'note':    'economic-equilibrium effect only (taxes fixed)',
    },
    # Comparative statics in the two system characteristics at the baseline rho, on a full CARTESIAN
    # (eps, theta) grid: one curve in eps per theta, which is what figure ARG_LOG_FourInOne shades.
    # eps does NOT track theta here -- the two are independent axes by construction (see the script).
    'epsThetaGrid': {
        'script':  'sweepEpsThetaGrid.py',
        'args':    lambda cx: (['--rho', str(C.ARG['ρBaseline'])] + GRIDFLAGS
                               + ['--pkldir', C.argInstanceDir(cx),
                                  '--out', 'epsThetaGrid_rho{ρ:.4f}' + C.argVariantTag(cx) + '.csv']),
        'outputs': lambda cx: {'sweep': C.argEpsThetaCsv(C.ARG['ρBaseline'], cx)},
        # The ONLY entry here whose script resumes from its own csv: it is keyed on (eps, theta) and
        # returns rows solved under any earlier calibration or grid setting without a word. So --force
        # has to REACH it. The two shock scripts overwrite their csvs outright and need nothing.
        'force':   ['--force'],
        'note':    'tau/sr/h/iota over the eps x theta plane',
    },
    # --- PRE-PUBLICATION (`part` = 'prepub'; every entry above is 'main'). Stationary vs date-specific
    # policy functions, the check behind sec:numerical's literature paragraph -- the Argentina counterpart
    # of runShocksUS.py's 'stationary'. Even LOG carries a state here (iota), so rho = 1 is included. Run
    # on the HEADLINE variant only: X enters no aggregate, so the two variants would print the same numbers.
    # ~2 h (a full CRRA recursion per distinct nu at rho = 0.5 and 2). No paper output reads it.
    'stationary': {
        'part':     'prepub',
        'headline': True,
        'script':   'stationaryApprox.py',
        'args':     lambda cx: (['--rho'] + [str(r) for r in C.ARG['ρTable']] + GRIDFLAGS
                                + ['--pkldir', C.argInstanceDir(cx),
                                   '--out', 'ARG_stationaryApprox' + C.argVariantTag(cx) + '.csv']),
        'outputs':  lambda cx: {'csv': os.path.join(C.NUMDIR, 'ARG_stationaryApprox' + C.argVariantTag(cx) + '.csv')},
        'note':     'stationary vs date-specific policies, LOG and CRRA (prepub, headline variant)',
    },
}


def selected(prepub = False, all_ = False):
    """ The entry names of the requested part(s), in declaration order; an entry without `part` is 'main'. """
    parts = {'main', 'prepub'} if all_ else ({'prepub'} if prepub else {'main'})
    return [n for n, e in EXPERIMENTS.items() if e.get('part', 'main') in parts]


def status(name, cx = False):
    """ (script exists, {key: path} present, {key: path} missing) for one experiment and variant. """
    e = EXPERIMENTS[name]
    have, lack = {}, {}
    for k, path in e['outputs'](cx).items():
        (have if os.path.exists(path) else lack)[k] = path
    return os.path.exists(os.path.join(C.MODELDIR, e['script'])), have, lack


def command(name, cx = False, ρ = None, force = False):
    """ `force` appends the entry's own re-do flag, for the scripts that resume from their own output.
    Skipping our own outputs is this file's job; making a resumable script redo work is the script's,
    and only it knows the flag -- so it is declared per entry rather than assumed to be '--force'. """
    e = EXPERIMENTS[name]
    cmd = [C.PYTHON, os.path.join(C.MODELDIR, e['script'])] + list(e['args'](cx))
    ρ = ρ if ρ is not None else e.get('rho')
    if ρ:
        cmd += ['--rho'] + [str(v) for v in ρ] if isinstance(ρ, (list, tuple)) else ['--rho', str(ρ)]
    return cmd + (list(e.get('force', [])) if force else [])


def main():
    p = argparse.ArgumentParser(description = __doc__.split('\n')[1])
    p.add_argument('--only', nargs = '+', choices = list(EXPERIMENTS), default = None,
                   help = 'a subset, regardless of part')
    p.add_argument('--commonX', action = 'store_true', help = 'also run every experiment for the common-X variant')
    p.add_argument('--prepub', action = 'store_true', help = 'the pre-publication part instead of the main one')
    p.add_argument('--all', dest = 'all_', action = 'store_true', help = 'both parts')
    p.add_argument('--list', action = 'store_true', help = 'report what exists and exit')
    p.add_argument('--force', action = 'store_true', help = 're-run even where output already exists')
    p.add_argument('--dry', action = 'store_true', help = 'print commands and exit')
    a = p.parse_args()
    names = a.only if a.only else selected(a.prepub, a.all_)
    variants = [False] + ([True] if a.commonX else [])
    # a 'headline' entry runs once, on the headline variant, whatever --commonX says
    variantsOf = lambda name: [C.ARG['commonX']] if EXPERIMENTS[name].get('headline') else variants

    if a.list or a.dry:
        for name in names:
            for cx in variantsOf(name):
                ok, have, lack = status(name, cx)
                print('{:<9} {:<13} {:<48} script:{}  have {}/{}'.format(
                    'commonX' if cx else 'vectorX', name, EXPERIMENTS[name]['note'],
                    'yes' if ok else 'MISSING', len(have), len(have)+len(lack)))
                if a.dry and ok:
                    print('    ' + ' '.join(command(name, cx, force = a.force)))
                if lack and not a.dry:
                    keys = sorted(lack, key = str)
                    print('    missing: ' + ', '.join(str(k) for k in keys[:8])
                          + (' ...' if len(keys) > 8 else ''))
        return

    for name in names:
        for cx in variantsOf(name):
            ok, have, lack = status(name, cx)
            label = '{} [{}]'.format(name, 'common X' if cx else 'vector X')
            if not ok:
                print('SKIP {}: {} does not exist yet.'.format(label, EXPERIMENTS[name]['script']))
                continue
            if not lack and not a.force:
                print('SKIP {}: all {} output(s) present.'.format(label, len(have)))
                continue
            # The scripts are individually resumable and skip what they already have, so the whole set is
            # handed over rather than only the missing keys -- one process, one warm-started march.
            cmd = command(name, cx, list(lack) if (lack and not a.force
                                                   and all(isinstance(k, float) for k in lack)) else None,
                          force = a.force)
            print('\n' + '='*94 + '\n{}: {}\n  {}\n'.format(label, EXPERIMENTS[name]['note'], ' '.join(cmd))
                  + '='*94)
            r = subprocess.run(cmd, cwd = C.REPO)
            if r.returncode:
                raise SystemExit('{} exited {}'.format(EXPERIMENTS[name]['script'], r.returncode))


if __name__ == '__main__':
    main()

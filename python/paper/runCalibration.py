r""" Stage (i) of the paper pipeline: produce the calibration the paper's numbers rest on.

Run:  .venv\Scripts\python.exe python\paper\runCalibration.py            check what exists, solve what does not
      ... --commonX                                                      also sweep the common-X variant
      ... --force                                                        re-solve every point
      ... --summaryOnly                                                  skip solving; just rebuild the summary
      ... --dry                                                          print the command and exit

Two things happen here, and only the first is expensive:

  1. The rho sweep. Delegated to python/InformalSavings/calibrateRhoGrid.py at config.ARG's grid and
     settings -- this file is the DECLARATION of what the paper calibrates, not a second implementation
     of it. ~75 minutes cold for the 16 points; the sweep is resumable and returns already-solved rho
     from its own csv, so a re-run with nothing to do costs seconds. Two calibration variants (vector
     X_i, the default, and common X with --commonX), each with its own csv and instance directory
     (config.argSweepCsv / argInstanceDir).

  2. The summary, results/paper/calibrationSummary.csv, one row per variant that has an instance.
     Stage (iii) builds the calibration TABLE from this file rather than from a pickled instance, so
     that stage stays a pure csv -> tex step that needs no model import (README.md). It is the one point
     in the pipeline that opens a pickle.

WHY THE SUMMARY IS NOT JUST THE SWEEP CSV: informalSavings_rhoGrid.csv carries the four CALIBRATED
parameters (beta, omega, eta0, X0) but not the ones that are fixed or derived -- eps, theta, gamma0,
alpha, xi, nu, eta_i, X_i -- and those are half of the paper's calibration table. They live on the
instance, which is where this reads them from.
"""
import os, sys, argparse, subprocess, pickle, datetime, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C


def sweepCmd(commonX = False, force = False):
    g = C.ARG['gridSettings']
    ρ = C.ARG['ρGrid']
    step = round(ρ[1] - ρ[0], 10)
    cmd = [C.PYTHON, os.path.join(C.MODELDIR, 'calibrateRhoGrid.py'),
           '--lo', str(min(ρ)), '--hi', str(max(ρ)), '--step', str(step),
           '--anchor', str(C.ARG['ρAnchor']),
           '--nι', str(g['nι']), '--ns', str(g['ns']),
           '--interpKind', g['interpKind'], '--smoothKnots', str(g['smoothKnots'])]
    x0 = C.ARG.get('anchorGuess')
    if x0:
        cmd += ['--x0'] + [str(x0[k]) for k in ('β', 'ω', 'η0', 'X0')]
    return cmd + (['--commonX'] if commonX else []) + (['--force'] if force else [])


def missing(commonX = False):
    """ rho in the paper's grid that the variant's sweep csv does not already carry. """
    path = C.argSweepCsv(commonX)
    if not os.path.exists(path):
        return list(C.ARG['ρGrid'])
    have = set(np.round(pd.read_csv(path)['ρ'].values, 6))
    return [ρ for ρ in C.ARG['ρGrid'] if round(ρ, 6) not in have]


def summarise(ρ = None, commonX = False):
    """ The calibration table's content, at the paper's baseline rho, as one flat record.

    Unpickling needs ModelInformalSavings importable and test.py's data path resolvable, which is why
    this chdir's into the model directory exactly as the experiment scripts do. Held to this one
    function so that nothing else in the pipeline inherits the requirement.
    """
    ρ = C.ARG['ρBaseline'] if ρ is None else ρ
    cwd = os.getcwd()
    sys.path.insert(0, C.MODELDIR)
    os.chdir(C.MODELDIR)
    try:
        with open(os.path.join(C.argInstanceDir(commonX), 'rho_{:.4f}.pkl'.format(ρ)), 'rb') as f:
            m = pickle.load(f)
        t0 = m.db['t0']
        first = lambda k: float(np.asarray(m.db[k])[0]) if np.ndim(m.db[k]) else float(m.db[k])
        atT0  = lambda k: float(m.db[k].xs(m.db['t'][t0])) if hasattr(m.db[k], 'xs') else first(k)
        # Vectors go through the csv as JSON, not as a python repr: under numpy 2 the repr of a list of
        # np.float64 is 'np.float64(1.64...)', whose literal text a number-scraping reader mines a
        # spurious 64.0 out of. json.dumps of plain floats has no such reading.
        vec = lambda a: json.dumps([float(v) for v in np.asarray(a, dtype = float).ravel()])
        rec = {'ρ': ρ, 'commonX': bool(getattr(m, 'commonX', False)),
               'preferences': m._calPreferences(), 't0': t0, 'T': m.T,
               'ε': atT0('eps'), 'θ': atT0('θ'), 'γ0': atT0('γ0'),
               'α': atT0('α'), 'ξ': atT0('ξ'), 'η0': atT0('η0'), 'X0': atT0('X0'),
               'ν': vec(m.db['ν'].values),
               'ηi': vec(m.db['ηi'].xs(m.db['t'][t0])),
               'Xi': vec(m.db['Xi'].xs(m.db['t'][t0])),
               'γi': vec(m.db['γi'].xs(m.db['t'][t0])),
               # relative formal hours: data under vector X, the model's prediction under common X
               'zxi': vec(m.db['zxi'].xs(m.db['t'][t0])),
               'zxPredicted': vec(m.predictedRelativeHours()) if hasattr(m, 'predictedRelativeHours')
                              else json.dumps([])}
        if bool(getattr(m, 'commonX', False)) != bool(commonX):
            raise SystemExit('instance at {} is the {} variant, not the requested one.'.format(
                C.argInstanceDir(commonX), 'common-X' if rec['commonX'] else 'vector-X'))
    finally:
        os.chdir(cwd)
        sys.path.remove(C.MODELDIR)

    # beta/omega come from the sweep csv, not the pickle: the csv is the record of the calibration and
    # is what every other paper number is read against, so the table must not be able to disagree with it.
    df = pd.read_csv(C.argSweepCsv(commonX))
    row = df.loc[np.isclose(df['ρ'], ρ)]
    if row.empty:
        raise SystemExit('rho={} is not in {} -- run the sweep first.'.format(
            ρ, os.path.relpath(C.argSweepCsv(commonX), C.REPO)))
    row = row.iloc[-1]
    for k in ('β', 'ω', 'KY', 'sr', 'τ', 'ι', 'residual', 'verifyResidual', 'commit', 'timestamp'):
        rec[k] = row[k]
    for k in ('X', 'hbar'):
        rec[k] = row[k] if k in row.index else np.nan
    # eta0/X0 are calibrated AND on the instance; they must agree, and a mismatch means the pickle and
    # the csv row came from different runs.
    for k in ('η0', 'X0'):
        if not np.isclose(rec[k], row[k], rtol = 1e-8):
            raise SystemExit('{}: instance {!r} != csv {!r}. Pickle and csv row are from different runs.'
                             .format(k, rec[k], row[k]))
    rec['builtAt'] = datetime.datetime.now().replace(microsecond = 0).isoformat()
    return rec


def main():
    p = argparse.ArgumentParser(description = __doc__.split('\n')[1])
    p.add_argument('--force', action = 'store_true', help = 're-solve every point, not only the missing')
    p.add_argument('--summaryOnly', action = 'store_true', help = 'rebuild the summary from what exists')
    p.add_argument('--commonX', action = 'store_true',
                   help = 'also SWEEP the common-X variant. The summary always covers every variant '
                          'that has an instance -- summarising is only an unpickle.')
    p.add_argument('--dry', action = 'store_true', help = 'print the sweep command and exit')
    p.add_argument('--rho', type = float, default = None, help = 'summarise a rho other than the baseline')
    a = p.parse_args()

    variants = [False] + ([True] if a.commonX else [])
    if a.dry:
        for cx in variants:
            print(' '.join(sweepCmd(cx, a.force)))
        return

    if not a.summaryOnly:
        for cx in variants:
            label = 'common X' if cx else 'vector X'
            todo = C.ARG['ρGrid'] if a.force else missing(cx)
            if todo:
                print('{}: calibrating {} point(s): {}'.format(label, len(todo), todo))
                cmd = sweepCmd(cx, a.force)
                print('  ' + ' '.join(cmd))
                r = subprocess.run(cmd, cwd = C.REPO)
                if r.returncode:
                    raise SystemExit('calibrateRhoGrid.py exited {}'.format(r.returncode))
            else:
                print('{}: all {} rho already solved; nothing to calibrate.'.format(label, len(C.ARG['ρGrid'])))

    # Every variant with an instance at the baseline rho, the headline first.
    ρ = C.ARG['ρBaseline'] if a.rho is None else a.rho
    recs = [summarise(a.rho, cx) for cx in (False, True)
            if os.path.exists(os.path.join(C.argInstanceDir(cx), 'rho_{:.4f}.pkl'.format(ρ)))]
    os.makedirs(C.PAPERDIR, exist_ok = True)
    out = os.path.join(C.PAPERDIR, 'calibrationSummary.csv')
    pd.DataFrame(recs).to_csv(out, index = False)
    print('\nwritten: ' + os.path.relpath(out, C.REPO))
    for rec in recs:
        print('  --- {} ---'.format('common X' if rec['commonX'] else 'vector X'))
        for k in ('ρ', 'β', 'ω', 'η0', 'X0', 'X', 'hbar', 'ε', 'θ', 'KY', 'sr', 'τ', 'residual'):
            print('  {:<6} {}'.format(k, rec[k]))


if __name__ == '__main__':
    main()

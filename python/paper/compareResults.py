r""" Row-by-row comparison of the csvs under results/ against the committed ones (reading a run,
archive/notes/todo_finalRun_2026-10-02.md): which files changed, which columns are new, and per numeric column
the largest absolute difference on rows aligned by the file's key columns.

    .venv\Scripts\python.exe python\paper\compareResults.py                 # every changed csv
    .venv\Scripts\python.exe python\paper\compareResults.py --ref HEAD~3    # against another commit
    .venv\Scripts\python.exe python\paper\compareResults.py --only esc      # path substring filter
    .venv\Scripts\python.exe python\paper\compareResults.py --tol 1e-12     # report columns above tol

Alignment: the columns common to both versions among KEYS, in order; a file with none of them is aligned on
position when the row counts agree and reported as 'positional'. Rows present on one side only are counted. Rows
that repeat a key are never dropped: the n-th repeat of a key in one version is compared with the n-th repeat in the
other, and the number of rows that share a key is reported per version. The committed version is read with
`git show <ref>:<path>`; nothing is written.
"""
import os, io, argparse, subprocess
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
KEYS = ['ρ', 'rho', 'spec', 'phi', 'commonX', 'method', 'host', 'scenario', 'θpinned', 'pos', 'country',
        'variant', 'experiment', 'shock', 'name', 'year', 't', 'preferences', 'timing', 'ξ', 'xi',
        'pinning', 'nu', 'wedge', 'thetaSource', 'turnout', 'omegaRule', 'retAge', 'betaRule', 'iso3', 'effect',
        'family', 'kind', 'wedgeFrom', 'date', 'ns', 'block', 'tss', 'terminal']
# eps and theta are keys only on the (eps, theta) grid sweeps; elsewhere θ is an output (the design path)
GRIDKEYS = ['eps', 'theta', 'θ']


def changed(ref):
    out = subprocess.run(['git', '-C', REPO, 'diff', '--name-only', ref, '--', 'results/'],
                         text = True, capture_output = True, encoding = 'utf-8').stdout.split()
    new = subprocess.run(['git', '-C', REPO, 'ls-files', '--others', '--exclude-standard', 'results/'],
                         text = True, capture_output = True, encoding = 'utf-8').stdout.split()
    return [p for p in out if p.endswith('.csv')], [p for p in new if p.endswith('.csv')]


def committed(ref, path):
    r = subprocess.run(['git', '-C', REPO, 'show', f'{ref}:{path}'], capture_output = True)
    if r.returncode != 0:
        return None
    return pd.read_csv(io.BytesIO(r.stdout))


def compare(old, new, tol, path = ''):
    allKeys = KEYS + (GRIDKEYS if 'epsThetaGrid' in path else [])
    keys = [k for k in allKeys if k in old.columns and k in new.columns]
    added = [c for c in new.columns if c not in old.columns]
    dropped = [c for c in old.columns if c not in new.columns]
    common = [c for c in new.columns if c in old.columns and c not in keys]
    if keys:
        # float keys are rounded so that 0.30000000000000004 and 0.3 align (a grid value written by two runs)
        old, new = old.copy(), new.copy()
        for k in keys:
            if pd.api.types.is_float_dtype(old[k]) and pd.api.types.is_float_dtype(new[k]):
                old[k], new[k] = old[k].round(12), new[k].round(12)
        dupOld, dupNew = int(old.duplicated(keys, keep = False).sum()), int(new.duplicated(keys, keep = False).sum())
        for d in (old, new):                       # the n-th repeat of a key meets the n-th repeat of it
            d['#repeat'] = d.groupby(keys, dropna = False).cumcount()
        o = old.set_index(keys + ['#repeat'], drop = True)
        n = new.set_index(keys + ['#repeat'], drop = True)
        both = o.index.intersection(n.index)
        onlyOld, onlyNew = len(o.index.difference(n.index)), len(n.index.difference(o.index))
        o, n = o.loc[both], n.loc[both]
        how = 'keys ' + ','.join(keys)
    elif len(old) == len(new):
        o, n, onlyOld, onlyNew, how, dupOld, dupNew = old, new, 0, 0, 'positional', 0, 0
    else:
        return {'how': 'unaligned', 'rowsOld': len(old), 'rowsNew': len(new), 'added': added,
                'dropped': dropped, 'diffs': {}}
    diffs = {}
    for c in common:
        a, b = o[c], n[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            d = (a.astype(float) - b.astype(float)).abs()
            m = float(np.nanmax(d.values)) if d.size and np.isfinite(d.values).any() else 0.0
            nanMismatch = int((a.isna() != b.isna()).sum())
            if m > tol or nanMismatch:
                diffs[c] = (m, nanMismatch)
        else:
            # pandas 3's astype(str) keeps a missing value missing, and NaN != NaN: two blanks are not a change
            ne = int(((a.astype(str) != b.astype(str)) & ~(a.isna() & b.isna())).sum())
            if ne:
                diffs[c] = ('text', ne)
    return {'how': how, 'rows': len(o), 'onlyOld': onlyOld, 'onlyNew': onlyNew, 'dupOld': dupOld, 'dupNew': dupNew,
            'added': added, 'dropped': dropped, 'diffs': diffs}


def main():
    p = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ref', default = 'HEAD')
    p.add_argument('--only', default = None, help = 'substring of the path')
    p.add_argument('--tol', type = float, default = 0.0)
    a = p.parse_args()
    mod, new = changed(a.ref)
    if a.only:
        mod = [x for x in mod if a.only in x]
        new = [x for x in new if a.only in x]
    print(f'{len(mod)} changed csv(s), {len(new)} new, against {a.ref}')
    for path in new:
        print(f'\nNEW   {path}')
    for path in mod:
        old = committed(a.ref, path)
        cur = pd.read_csv(os.path.join(REPO, path))
        if old is None:
            print(f'\n{path}: not in {a.ref}')
            continue
        r = compare(old, cur, a.tol, path)
        print(f'\n{path}  [{r["how"]}]')
        if r['how'] == 'unaligned':
            print(f'  rows {r["rowsOld"]} -> {r["rowsNew"]}, no key columns to align on')
        else:
            print(f'  rows aligned {r["rows"]}; only committed {r["onlyOld"]}; only new {r["onlyNew"]}')
            if r['dupOld'] or r['dupNew']:
                print(f'  rows sharing a key with another row: committed {r["dupOld"]}, new {r["dupNew"]}')
        if r['added']:
            print(f'  columns added: {", ".join(r["added"])}')
        if r['dropped']:
            print(f'  columns dropped: {", ".join(r["dropped"])}')
        if not r['diffs']:
            print(f'  every common column identical (tol {a.tol:g})')
        for c, (m, k) in sorted(r['diffs'].items(), key = lambda kv: -float(kv[1][0]) if kv[1][0] != 'text' else 0):
            if m == 'text':
                print(f'  {c:<24} text differs on {k} rows')
            else:
                print(f'  {c:<24} max|Δ| {m:.3e}' + (f', NaN pattern differs on {k} rows' if k else ''))


if __name__ == '__main__':
    main()

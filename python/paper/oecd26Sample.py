r""" Stage (0), no network: the sample of the 26-country test, the 29 OECD members of figure 1 less those whose
ILO informal employment exceeds INFORMAL_MAX percent (Mexico, Turkey and Korea).

Run:  .venv\Scripts\python.exe -I python\paper\oecd26Sample.py [--force]   -> data/oecd26_sample.csv

Reads data/oecd29.csv (oecd29Data.py) and writes one row per country of the 29: the informality datum it is
judged on (ILOSTAT SDG 8.3.1, % of total employment, the year nearest 2020 within 2018-2022, with its status
flag), `sample26` and the reason. A country with no informality series (AUS, CAN, JPN, NZL, USA) is kept: the
rule drops on an observed share above the threshold, never on a missing one. The model drivers
(python/US/runOECD26.py) copy `sample26` onto their rows and the stage (iii) builders filter on this file, so
the rule is applied here and nowhere else.
"""
import os, sys, argparse
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

SRC = os.path.join(C.DATA, 'oecd29.csv')
OUT = os.path.join(C.DATA, 'oecd26_sample.csv')
INFORMAL_MAX = 25.0                      # percent of total employment


def sample(d):
    """ The sample table from data/oecd29.csv's frame `d` (indexed by iso3, in its order). """
    inf = pd.to_numeric(d['informal'], errors = 'coerce')
    keep = ~(inf > INFORMAL_MAX)                       # NaN compares False: no series, kept
    note = []
    for c in d.index:
        if np.isnan(inf[c]):
            note.append('no ILO informality series; kept')
        else:
            note.append('informal employment {:.1f}% ({}{}) {} {:g}%'.format(
                inf[c], int(d.at[c, 'informal_year']),
                ', status ' + d.at[c, 'informal_status'] if isinstance(d.at[c, 'informal_status'], str)
                and d.at[c, 'informal_status'] else '', '>' if not keep[c] else '<=', INFORMAL_MAX))
    out = pd.DataFrame({'country': d['country'], 'informal': inf,
                        'informal_year': pd.to_numeric(d['informal_year']).round().astype('Int64'),
                        'informal_status': d['informal_status'].fillna(''), 'informalMax': INFORMAL_MAX,
                        'sample26': keep.astype(bool), 'note': note}, index = d.index)
    out.index.name = 'iso3'
    return out


def main():
    p = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--force', action = 'store_true', help = 'rebuild the output even if it exists')
    a = p.parse_args()
    if os.path.exists(OUT) and not a.force:
        print('up to date: {} (--force to rebuild)'.format(os.path.relpath(OUT, C.REPO)))
        return
    s = sample(pd.read_csv(SRC, index_col = 'iso3'))
    s.to_csv(OUT, float_format = '%.8g')
    print('wrote {}: {} of {} in the sample; dropped {}'.format(
        os.path.relpath(OUT, C.REPO), int(s['sample26'].sum()), len(s), ', '.join(s.index[~s['sample26']])))


if __name__ == '__main__':
    main()

r""" Stage (0): build the Argentina calibration target that comes from outside the repo.

Run:  .venv\Scripts\python.exe python\paper\dataTargets.py [--force] [--lo 1980] [--hi 2010]

The one target this produces is Argentina's **capital-output ratio** $K/Y$ in annual output units, which
is what identifies $\beta$ (`eq:calibration:KY`). The target is the ratio **at the calibration year**,
which is where the model's other aggregate targets are measured; the mean over the thirty years ending
there -- one model period of data -- is written out beside it as the sensitivity, and `--target window`
selects it instead. See `archive/notes/argentina_calibrationTarget.md` for why the calibration targets this
ratio at all.

Source: Penn World Table 11.0, capital stock and GDP at constant national prices (`rnna`, `rgdpna`), as
mirrored by FRED. Both series are in the same units, so the ratio needs no deflator and no PPP
conversion. The United States is fetched alongside Argentina purely as a reference column -- the US arm
of the model targets the same object through its 30-year interest rate (`R0`), and $K/Y = 30\alpha/R$.

Written outputs (both under `data/`, both committed, so nothing downstream needs the network):
  * `argentina_capitalOutput.csv`      the annual series, i.e. the evidence
  * `argentina_calibrationTargets.csv` the target, its window, source and retrieval date, and
                                       the calibration year on its own as a sensitivity

`python/InformalSavings/test.py` reads the second file. Like the rest of the pipeline this script skips
work whose output already exists; `--force` refetches.
"""
import os, io, argparse, datetime, urllib.request
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(REPO, 'data')

SERIES = os.path.join(DATA, 'argentina_capitalOutput.csv')
TARGETS = os.path.join(DATA, 'argentina_calibrationTargets.csv')

# PWT 11.0 via FRED. K and Y for one country are in the same constant national prices, so K/Y is a pure
# ratio; do NOT mix a K from one of these with a Y from another source.
FRED = {'K_ARG': 'RKNANPARA666NRUG', 'Y_ARG': 'RGDPNAARA666NRUG',
        'K_USA': 'RKNANPUSA666NRUG', 'Y_USA': 'RGDPNAUSA666NRUG'}
SOURCE = ('Penn World Table 11.0 (rnna/rgdpna, constant national prices) via FRED: '
          + ', '.join(sorted(FRED.values())))
# args.hi is the calibration year AND the end of the averaging window: the target is read at the former,
# the sensitivity over the latter. Argentina's K/Y swings between 3.17 and 4.28 inside this window through
# the denominator alone (up through the 1989-90 inflation crisis, up again through the 2001-02 collapse,
# down through the 2003-07 recovery), so the two readings are 13% apart and which to target is a
# judgement, not a measurement: the calibration year wins because every other target in eq:calibration --
# the tax rate, the replacement-rate ratio, the coverage share, the household survey -- is measured at or
# around 2010, and a target averaged over a different span would not describe the same economy.
WINDOW = (1980, 2010)


def fetch():
    """ The four FRED series as one annual frame indexed by year. """
    url = 'https://fred.stlouisfed.org/graph/fredgraph.csv?id=' + ','.join(FRED.values())
    raw = urllib.request.urlopen(url, timeout = 60).read().decode()
    df = pd.read_csv(io.StringIO(raw), parse_dates = ['observation_date'])
    df = df.rename(columns = {v: k for k, v in FRED.items()})
    df['year'] = df['observation_date'].dt.year
    return df.set_index('year')[list(FRED)]


def capitalOutput(df):
    """ K/Y per country, in annual output units -- the object the calibration targets. """
    return pd.DataFrame({'KY_ARG': df['K_ARG']/df['Y_ARG'], 'KY_USA': df['K_USA']/df['Y_USA']})


def main():
    p = argparse.ArgumentParser(description = __doc__,
                               formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lo', type = int, default = WINDOW[0], help = 'first year of the averaging window')
    p.add_argument('--hi', type = int, default = WINDOW[1], help = 'last year of the averaging window')
    p.add_argument('--target', choices = ('year', 'window'), default = 'year',
                   help = "which reading is THE target: the calibration year (--hi) or the window mean. "
                          "The other is written out beside it and read by nothing.")
    p.add_argument('--force', action = 'store_true', help = 'refetch even if the outputs exist')
    args = p.parse_args()

    if os.path.exists(SERIES) and os.path.exists(TARGETS) and not args.force:
        t = pd.read_csv(TARGETS)
        print('up to date: {} (--force to refetch)'.format(os.path.relpath(TARGETS, REPO)))
        print(t.to_string(index = False))
        return

    ky = capitalOutput(fetch()).dropna()
    ky.round(6).to_csv(SERIES)
    window = ky.loc[args.lo:args.hi, 'KY_ARG']
    if len(window) != args.hi - args.lo + 1:
        raise ValueError('window {}-{} has {} observations, not {}'.format(
            args.lo, args.hi, len(window), args.hi - args.lo + 1))

    # Both readings are written every time; --target only decides which one is named capitalOutputRatio,
    # i.e. which one the models read. The other keeps its own name and is read by nothing, so the
    # sensitivity to the choice is on the record rather than in a note.
    reading = {'year':   (float(ky.loc[args.hi, 'KY_ARG']), str(args.hi)),
               'window': (float(window.mean()), '{}-{}'.format(args.lo, args.hi))}
    other = 'window' if args.target == 'year' else 'year'
    common = {'units': 'annual K/Y', 'source': SOURCE,
              'retrieved': datetime.date.today().isoformat()}
    rows = [{'target': 'capitalOutputRatio', 'value': round(reading[args.target][0], 4),
             'window': reading[args.target][1],
             'note': 'THE TARGET ({}): identifies beta via eq:calibration:KY'.format(args.target),
             **common},
            {'target': 'capitalOutputRatio_' + other, 'value': round(reading[other][0], 4),
             'window': reading[other][1],
             'note': 'reference only: the {} reading, read by nothing'.format(other), **common}]
    pd.DataFrame(rows).to_csv(TARGETS, index = False)

    print('{}: {} rows, {}-{}'.format(os.path.relpath(SERIES, REPO), len(ky), ky.index[0], ky.index[-1]))
    print('K/Y Argentina {}-{}: {:.4f}   min {:.3f} ({})  max {:.3f} ({})'.format(
        args.lo, args.hi, window.mean(), window.min(), window.idxmin(), window.max(), window.idxmax()))
    print('K/Y Argentina {} alone:  {:.4f}   ({:+.1f}% against the window mean)'.format(
        args.hi, ky.loc[args.hi, 'KY_ARG'], 100*(ky.loc[args.hi, 'KY_ARG']/window.mean() - 1)))
    print('K/Y United States {}-{}: {:.4f}   (reference: the arm that targets R)'.format(
        args.lo, args.hi, ky.loc[args.lo:args.hi, 'KY_USA'].mean()))
    print('{}: capitalOutputRatio = {}  ({} reading)'.format(
        os.path.relpath(TARGETS, REPO), rows[0]['value'], args.target))


if __name__ == '__main__':
    main()

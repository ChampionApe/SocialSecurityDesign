r""" Stage (0): the pension design θ (and the spending behind τ) of the 26 countries under the variants of
notes/thetaTau_alternatives.md, for python/US/oecd26.py's thetaSource (variant C is the paper's, config.US['theta'])
and the data part of the technical note.

Run:  .venv\Scripts\python.exe -I python\paper\oecd26Theta.py [--force]
      -> data/oecd26_theta.csv, one row per country x variant: iso3, variant, theta, thetaRaw, s,
         spending, privShare, inSample, note, tau

Variants (1 - θ is the share of spending not linked to own earnings; θ the slope of relative benefits on
relative earnings):
  paper   the earlier replacement-rate index, not the paper's measure (figure 1.1 plots 1/(2 - θ_C)):
          θ = 2 - 1/I, I = RR(1 AW)/RR(0.5 AW), mandatory schemes, clipped to [0, 1]
  A       public tier: the least-squares slope c/(a+c) of the pension a + c·y on y at 0.5, 1, 2 AW (PaG 2021
          table 4.2, men; NZL's low earner at 0.6 AW, the table's footnote) times (1 - s); GBR in force
  A2pt    as A on 0.5 and 1 AW only
  B       all mandatory and quasi-mandatory tiers: the slope on table 4.1 times (1 - s·S/S_B); spending S_B
  C       A, with inSample False where funded tiers pay at least a quarter of mandatory benefits in 2020
          (privShare) or A's slope is negative (a public pension tested against the private one)
θ outside [0, 1] is NaN in `theta` (thetaRaw keeps it), except B, clipped to 1 (SWE's occupational tier above
the public ceiling), and paper, clipped to [0, 1]. Spending S is SOCX public old-age and survivors cash 2020 (figure 1.1's), τ = S/(1 - α).

s (S_TOPUP) is 2020 spending on minimum-pension top-ups and tested benefits that PaG's full-career worker at
0.5 AW does not reach, and that SOCX books inside public old-age and survivors cash, over S. Values in % of
GDP from data/oecd26_topups.csv (ITA FRA AUT DEU BEL PRT HUN) and the scheme files
data/oecd26_pensionSchemes_A.csv, _B.csv (ESP FIN USA SWE JPN CHE), all hand-collected drafts, typed in here;
ITA's integrazione and PRT's complementos include their invalidity part (upper bounds), BEL lacks its
minimum-pension supplement (lower bound), LUX is a gap (0); NOR's and POL's minimum is reached at 0.5 AW (0).
"""
import os, sys, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(DATA, 'oecd26_theta.csv')
ALPHA = 0.30
PRIVMAX = 0.25

S_TOPUP = {   # % of GDP, 2020 (SWE 2018); source in the module docstring
    'ITA': 0.4216 + 0.0903 + 0.2911 + 0.0049,   # integrazione, maggiorazioni, assegno sociale, importo aggiuntivo
    'FRA': 0.3679 + 0.1523,                     # minima de pension (all schemes), minimum vieillesse
    'AUT': 0.2265,                              # Ausgleichszulage on old-age and survivors pensions (SOCX lines)
    'DEU': 0.0521,                              # SOCX old age, other cash benefits (social assistance)
    'BEL': 0.1572,                              # GRAPA
    'PRT': 0.9586 + 0.1055 + 0.0587,            # complementos sociais, CSI, pensões sociais
    'HUN': 0.0083 + 0.0125,                     # minimum top-up (upper bound), SOCX means-tested
    'ESP': 0.62, 'FIN': 0.74, 'USA': 0.06, 'SWE': 0.28, 'JPN': 0.04, 'CHE': 0.45,
}
GBR_INFORCE = 0.81/(2.96 + 0.81 + 0.48 + 0.22)  # S2P / (basic SP + S2P + new SP + Pension Credit), 2019-20
QUASI = ('NLD', 'SWE')                          # quasi-mandatory occupational tiers SOCX books as voluntary


def slope(rr, lo):
    """ θ = c/(a + c) of the least-squares line a + c·y through the pensions rr·y at y = (lo, 1, 2). """
    y = np.array([lo, 1., 2.])
    a, c = np.linalg.lstsq(np.c_[np.ones(3), y], np.asarray(rr, float)/100*y, rcond = None)[0]
    return c/(a + c)


def build():
    d = pd.read_csv(os.path.join(DATA, 'oecd26_pensionTiers_harmonised.csv'), comment = '#')
    d['year'] = d['year'].astype(str)
    fig = pd.read_csv(os.path.join(DATA, 'oecdFigure1.csv'), index_col = 'iso3')
    smp = pd.read_csv(os.path.join(DATA, 'oecd26_sample.csv'), index_col = 'iso3')
    isos = list(smp.index[smp['sample26']])

    def rr(item, tier, iso):
        x = d[(d.item == item) & (d.tier == tier) & (d.iso3 == iso) & (d.year == '2020') & (d.sex == 'M')]
        return x.set_index('earnings_multiple')['value'].reindex([0.5, 1., 2.]).to_numpy(float)

    sel = lambda it, tier: d[(d.item == it) & (d.tier == tier) & (d.year == '2020')].set_index('iso3')['value']
    mp, vp = sel('socx_TP01_cash', 'mandatoryPrivate'), sel('socx_TP01_cashInKind', 'voluntaryPrivate')
    rows = []
    for iso in isos:
        lo = 0.6 if iso == 'NZL' else 0.5
        S = float(fig.at[iso, 'pensionSpending'])
        top = S_TOPUP.get(iso, 0.)
        SB = S + np.nan_to_num(mp.get(iso, 0.)) + (np.nan_to_num(vp.get(iso, 0.)) if iso in QUASI else 0.)
        priv = 1 - S/SB
        pub = rr('grr_tier', 'public_mandatory', iso)
        b = pub/100*np.array([lo, 1., 2.])
        two = (b[1] - b[0])/((1 - lo)*b[1]) if b[1] > 0 else -np.inf
        aRaw = GBR_INFORCE if iso == 'GBR' else slope(pub, lo)*(1 - top/S)
        a2Raw = GBR_INFORCE if iso == 'GBR' else two*(1 - top/S)
        bRaw = slope(rr('grr_mandatory', 'mandatory', iso), lo)*(1 - top/SB)
        I = float(fig.at[iso, 'bbIndex'])
        ok = lambda t: min(max(t, 0.), 1.) if -1e-9 <= t <= 1. + 1e-9 else np.nan   # IRL, NLD: 0 to rounding
        inC = bool(priv < PRIVMAX and aRaw >= -1e-9)
        note = 'GBR in force (S2P share)' if iso == 'GBR' else ''
        base = {'iso3': iso, 's': top/S, 'privShare': priv}
        rows += [base | {'variant': 'paper', 'theta': min(max(2 - 1/I, 0.), 1.), 'thetaRaw': 2 - 1/I,
                         'spending': S, 'inSample': True, 'note': 'figure 1.1 index'},
                 base | {'variant': 'A', 'theta': ok(aRaw), 'thetaRaw': aRaw, 'spending': S, 'inSample': True,
                         'note': note},
                 base | {'variant': 'A2pt', 'theta': ok(a2Raw), 'thetaRaw': a2Raw, 'spending': S,
                         'inSample': True, 'note': note},
                 base | {'variant': 'B', 'theta': min(max(bRaw, 0.), 1.) if bRaw >= -1e-9 else np.nan, 'thetaRaw': bRaw,
                         'spending': SB, 'inSample': True, 'note': ''},
                 base | {'variant': 'C', 'theta': ok(aRaw), 'thetaRaw': aRaw, 'spending': S, 'inSample': inC,
                         'note': note}]
    out = pd.DataFrame(rows)[['iso3', 'variant', 'theta', 'thetaRaw', 's', 'spending', 'privShare', 'inSample',
                              'note']]
    out['tau'] = out['spending']/100/(1 - ALPHA)
    return out


def main():
    p = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--force', action = 'store_true')
    a = p.parse_args()
    if os.path.exists(OUT) and not a.force:
        print('{} exists (--force rebuilds)'.format(os.path.relpath(OUT, REPO)))
        return 0
    out = build()
    out.to_csv(OUT, index = False, float_format = '%.6f')
    print(out.pivot(index = 'iso3', columns = 'variant', values = 'theta').round(3).to_string())
    print('C sample:', int(out[(out.variant == 'C') & out.inSample].shape[0]), 'countries')
    return 0


if __name__ == '__main__':
    sys.exit(main())

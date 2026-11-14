r""" Diagnostics of the universalisation reform at the calibration year t0 (2010), per rho: what the paper's
appendix and section 5 quote about the reform beyond its headline numbers. Not a test -- an experiment script
on the calibrated instances, sharing shockUniversal.py's code path (loadCalibrated, stateAtT0, createCopyFromt0,
installEps, universalEps, solvePEE) and shockEEOnly.py's eeSolve.

Run:  .venv\Scripts\python.exe python\InformalSavings\reformDiagnostics.py --rho 0.5 1.0 2.0
          --pkldir <instances> [--refine 60 --refineRho 0.5 0.6] [--tag _commonX]

Per rho it solves the baseline PEE, the no-shock copy and the reform copy from the baseline's state entering t0
(eq:shockState, the retirees' formal shares held), the economic equilibrium at the baseline tax with the
reformed eps, and writes four csvs into results/shocks/ (suffix --tag):

  ARG_reformFOC<tag>.csv   The t0 first order condition z_t (eq:zdecomposition; eq:PEELOG under LOG,
      eq:PEEcrraTerms with eq:dv20 under CRRA) at the seed state (s_{t0-1}, iota_{t0-1}), split by generation and
      type: the old's terms into the pension part of dln c_2/dtau (sigma/tau, where eps enters directly through
      eq:governmentBudget) and the labour part ((1-alpha) dln h_t/dtau_t), each times its level factor c_2^{1-1/rho}
      under CRRA, and the young's numerical terms dv_1. Rows per group: z on the no-shock copy and on the reform
      copy at the two FOC roots tauB and tauR, the reform's push zR(tauB) - zB(tauB), the secant slope
      (zR(tauR) - zR(tauB))/(tauR - tauB), so that tauR - tauB = sum(push)/(-sum(slope)) exactly.
  ARG_reformB5<tag>.csv    Both sides of the paper's condition eq:LOG:epsCondition (appendix B) on the solved path at
      t0, before (the baseline, pre-reform eps) and after (the reform, eps^U): sigma^0 = b^0/c_2^0 and sigma^i =
      b^i/c_2^i (eq:EE:c0, eq:EE:ci), phi_t = 1 + (1-alpha) tau_t dln h_t/dtau_t from the FOC's own derivative,
      lhsB5 = (phi/eps) sigma^0 (1-sigma^0), lhsSavings = sigma^0 (1-sigma^0)/eps (the informal-savings model's own
      coefficient: its informal retirees' savings income moves with h_t as formal retirees' does), rhs =
      sum_i gamma_i sigma^i (1-sigma^i). The condition is the log model's; at rho != 1 the row evaluates it on the
      CRRA calibration and is not the CRRA condition.
  ARG_reformHouseholdSavings<tag>.csv   Household savings rates of the young at t0, savings over disposable
      labour income: informal s_{t,0}/(w^0 eta_0 h_0) with s_{t,0} = iota_t s_t (eq:EE:s0_s), formal s_{t,i}/((1-tau)
      w eta_i h_{t,i}) (eq:EE:si_s) and their gamma-weighted average s_t/((1-tau) w h_t); baseline, economic
      equilibrium (tax at the baseline path) and full effect, and the difference-in-differences informal minus
      formal.
  ARG_reformRefine<tag>.csv  (--refine) The reform's t0 tax response re-solved on a finer state grid (nι = ns =
      --refine; nι only at rho = 1), beside the working grid's.
"""
import os, sys, argparse, time
import numpy as np, pandas as pd
from scipy import interpolate, optimize

sys.stdout.reconfigure(encoding = 'utf-8', line_buffering = True)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
os.chdir(HERE)                                  # relative paths resolve here; test.py finds data/ itself

import shockUniversal as su
import shockEEOnly as see
from gridsearch import CartesianGrid

OUTDIR = os.path.join(REPO, 'results', 'shocks')
GROUPS = ('oldInformalPension', 'oldInformalLabour', 'oldFormalPension', 'oldFormalLabour',
          'youngFormal', 'youngInformal')


def focComponents(mt, preferences, sols, seed):
    """ z_t of the copy's first period on the tax grid at the seed state, by group (GROUPS), plus dln h_t/dtau_t
    and the solver's own z (the check that the groups sum to it). The first period's objects are rebuilt as
    solveBackward_t builds them, at the seed (s0 is a single-point savings grid under CRRA), with the inherited
    formal shares held (eq:shockState). Returns (tauGrid, {group: (M,)}, dlnh (M,), zSolver (M,)). """
    P = getattr(mt, preferences)
    BG = P.BG
    tIdx = mt.db['t']
    t, t1 = tIdx[0], tIdx[1]
    solp = sols[t1]
    θ, ε = mt.db['θ'].values, mt.db['eps'].values
    τG = P.GS['PEE']['solGrids']['τ']
    si = np.asarray(seed['inherited']['si_s'], dtype = float)
    ι0, s0 = float(seed['ι0']), float(seed['s0'])
    α = BG.get('α', t)
    with BG.cacheParams():
        if preferences == 'LOG':
            ι, _ = P.solveStateApprox_t(τG, sols[t]['ιCand'], t, θ[1], ε[1], solp)
            d = P.stateGrid_t(τG, ι, t, t, θ[0], θ[1], ε[1], solp, si)
            parts = P.zbarParts_t(d)
            lev2i, lev20 = 1.0, 1.0
            zSolver = P._zState(P.zbarAtShares(θ[0], t, d['si_s_'], parts), τG, ε[0], np.array([ι0]), t)[:, 0]
        else:
            sG = np.array([s0])
            g = CartesianGrid(τ = τG, s_ = sG)
            sSol, ιSol, _, _ = P.solveStateApprox_t(τG, sG, sols[t]['sCand'], sols[t]['ιCand'], t, θ[1], ε[1], solp)
            d = P.stateGrid_t(g.flat['τ'], sSol.reshape(-1), ιSol.reshape(-1), g.flat['s_'], t, t, θ[0], θ[1],
                              ε[1], solp, si)
            parts = P.zbarParts_t(d, g, t)
            p = parts['p']
            lev2i = BG.c2i(parts['h'], parts['s_'], τG, θ[0], d['si_s_'], t)**p
            lev20 = BG.c20(parts['h'], parts['s_'], ε[0], τG, ι0, t)**p
            zbar = P.zbarAtShares(θ[0], t, d['si_s_'], parts)
            zSolver = P._zStateCRRA(zbar, d, ε[0], np.array([ι0]), t, g).reshape(-1)
        dlnh = np.asarray(parts['dlnh'], dtype = float)
        lab = (1-α)*dlnh
        pen2i = BG.dlnc2i_dτ(dlnh, τG, θ[0], d['si_s_'], t) - lab[:, None]
        A = BG.auxInf1_(ε[0], t)
        pen20 = A/(ι0 + A*τG)
        w2i, w20 = BG.get('γi[t-1]', t)*BG.ω2i(t), BG.get('γ0[t-1]', t)*BG.ω20(t)
        ν = BG.get('ν', t)
        w1i, w10 = ν*BG.get('γi', t)*BG.ω1i(t), ν*BG.get('γ0', t)*BG.ω10(t)
        comp = {'oldInformalPension': w20*lev20*pen20, 'oldInformalLabour': w20*lev20*lab,
                'oldFormalPension': (w2i*lev2i*pen2i).sum(axis = -1),
                'oldFormalLabour': (w2i*lev2i*lab[:, None]).sum(axis = -1),
                'youngFormal': (w1i*np.asarray(parts['dv1i'])).sum(axis = -1),
                'youngInformal': w10*np.asarray(parts['dv10'])}
    total = sum(comp.values())
    ok = np.isfinite(total) & np.isfinite(zSolver)
    if not np.allclose(total[ok], zSolver[ok], rtol = 1e-10, atol = 1e-12):
        raise RuntimeError('the groups do not sum to the solver\'s z_t: max|Δ| {:.2e}'.format(
            np.max(np.abs(total[ok] - zSolver[ok]))))
    return τG, comp, dlnh, zSolver


def _spline(x, y):
    ok = np.isfinite(y)
    return interpolate.CubicSpline(x[ok], y[ok])


def focRoot(τG, z, near):
    """ The crossing of z along the tax grid nearest `near` (the path's own t0 tax), on a cubic spline. """
    f = _spline(τG, z)
    ok = np.isfinite(z)
    x, y = τG[ok], z[ok]
    k = np.flatnonzero(np.sign(y[:-1]) != np.sign(y[1:]))
    if not k.size:
        raise RuntimeError('z_t has no crossing on the tax grid')
    k = k[np.argmin(np.abs(x[k] - near))]
    return float(optimize.brentq(f, x[k], x[k+1], xtol = 1e-14))


def b5Row(m, rep, t, τ, φ, label, ρ):
    """ eq:LOG:epsCondition's two sides on a solved report at t. """
    γ = np.asarray(m.B.get('γi[t-1]', t), dtype = float)
    ε = float(m.B.get('eps', t))
    σ0 = float(rep['b0'].xs(t)/rep['c20'].xs(t))
    σi = np.asarray(rep['bi'].xs(t), dtype = float)/np.asarray(rep['c2i'].xs(t), dtype = float)
    row = {'ρ': ρ, 'point': label, 'ε': ε, 'τ': τ, 'φ': φ, 'σ0': σ0, 'sumγ': float(γ.sum())}
    row |= {'σ{}'.format(j): v for j, v in zip(m.db['i'], σi)}
    row |= {'γ{}'.format(j): v for j, v in zip(m.db['i'], γ)}
    row['lhsB5'] = φ/ε*σ0*(1-σ0)
    row['lhsSavings'] = σ0*(1-σ0)/ε
    row['rhs'] = float((γ*σi*(1-σi)).sum())
    row['marginB5'] = row['lhsB5'] - row['rhs']
    return row


def householdSavings(m, rep, t, τ, scenario, ρ):
    """ Savings over disposable labour income of the young at t, by type (eq:EE:s0_s, eq:EE:si_s). """
    s, w, h = float(rep['s'].xs(t)), float(rep['w'].xs(t)), float(rep['h'].xs(t))
    ι = float(rep['ι'].xs(t))
    η0, ηi = float(m.B.get('η0', t)), np.asarray(m.B.get('ηi', t), dtype = float)
    inf = ι*s/(float(rep['w0'].xs(t))*η0*float(rep['h0'].xs(t)))
    si = np.asarray(rep['si_s'].xs(t), dtype = float)*s
    formal = si/((1-τ)*w*ηi*np.asarray(rep['hi'].xs(t), dtype = float))
    row = {'ρ': ρ, 'scenario': scenario, 'τ': τ, 'informal': inf, 'formalAvg': s/((1-τ)*w*h)}
    row |= {'formal{}'.format(j): v for j, v in zip(m.db['i'], formal)}
    return row


def runRho(ρ, settings, pkldir, epsPre0, refine):
    print('\n' + '='*94 + '\nrho={}'.format(ρ))
    tStart = time.time()
    m, pref = su.loadCalibrated(ρ, settings, pkldir = pkldir, epsPre0 = epsPre0)
    t0 = m.db['t'][m.db['t0']]
    base = su.solvePEE(m, pref)
    seed = m.stateAtT0(base['report'], t0, init = base['init'])
    τBpath = float(base['τ'].xs(t0))
    mB = m.createCopyFromt0(t0)
    outB = su.solvePEE(mB, pref, **seed)
    mR = m.createCopyFromt0(t0)
    εU = su.universalEps(mR)
    su.installEps(mR, εU)
    outR = su.solvePEE(mR, pref, **seed)
    τRpath = float(outR['τ'].iloc[0])
    print('{} solves ({:.0f}s): tau_t0 baseline {:.6f}, no-shock copy {:.6f}, reform {:.6f}, dtau {:+.4f} pp'.format(
        pref, time.time()-tStart, τBpath, float(outB['τ'].iloc[0]), τRpath, 100*(τRpath-τBpath)))

    # W15: the first order condition by group
    τG, cB, dlnhB, zB = focComponents(mB, pref, outB['sols'], seed)
    _, cR, dlnhR, zR = focComponents(mR, pref, outR['sols'], seed)
    τB, τR = focRoot(τG, zB, τBpath), focRoot(τG, zR, τRpath)
    foc = []
    for g in GROUPS + ('total',):
        yB = sum(cB.values()) if g == 'total' else cB[g]
        yR = sum(cR.values()) if g == 'total' else cR[g]
        fB, fR = _spline(τG, yB), _spline(τG, yR)
        foc.append({'ρ': ρ, 'group': g, 'zB_τB': float(fB(τB)), 'zR_τB': float(fR(τB)), 'zR_τR': float(fR(τR)),
                    'push': float(fR(τB) - fB(τB)), 'slopeR': float((fR(τR) - fR(τB))/(τR - τB)),
                    'slopeB_τB': float(fB(τB, 1)), 'τB': τB, 'τR': τR, 'τBpath': τBpath, 'τRpath': τRpath})
    tot = foc[-1]
    print('  FOC roots at the seed state: tauB {:.6f}, tauR {:.6f}; push {:+.4e}, slope {:+.4e} -> dtau {:+.4f} pp'
          .format(τB, τR, tot['push'], tot['slopeR'], -100*tot['push']/tot['slopeR']))
    for r in foc[:-1]:
        print('    {:<20} push {:+.4e}  slopeR {:+.4e}  z at tauB {:+.4e}'.format(r['group'], r['push'], r['slopeR'],
                                                                                 r['zB_τB']))

    # A12: eq:LOG:epsCondition on the solved paths at t0
    α = float(m.B.get('α', t0))
    φB = 1 + (1-α)*τBpath*float(_spline(τG, dlnhB)(τBpath))
    φR = 1 + (1-α)*τRpath*float(_spline(τG, dlnhR)(τRpath))
    b5 = [b5Row(m, base['report'], t0, τBpath, float(φB), 'pre', ρ),
          b5Row(mR, outR['report'], mR.db['t'][0], τRpath, float(φR), 'post', ρ)]
    for r in b5:
        print('  (epsCondition) {:<4} eps {:.4f} phi {:.4f} sigma0 {:.4f} lhsB5 {:.4f} lhsSavings {:.4f} rhs {:.4f}'.format(
            r['point'], r['ε'], r['φ'], r['σ0'], r['lhsB5'], r['lhsSavings'], r['rhs']))

    # #2: household savings rates of the young at t0
    τpath = base['τ'].loc[t0:].values
    mE = m.createCopyFromt0(t0)
    su.installEps(mE, su.universalEps(mE))
    _, repE = see.eeSolve(mE, pref, τpath, mE.db['θ'].values, mE.db['eps'].values, float(seed['s0']))
    hh = [householdSavings(m, base['report'], t0, τBpath, 'base', ρ),
          householdSavings(mE, repE, mE.db['t'][0], τBpath, 'ee', ρ),
          householdSavings(mR, outR['report'], mR.db['t'][0], τRpath, 'reform', ρ)]
    for r in hh:
        print('  savings rate {:<6} informal {:.4f}  formal avg {:.4f}'.format(r['scenario'], r['informal'], r['formalAvg']))

    ref = []
    if refine and any(np.isclose(ρ, v) for v in refine[1]):
        n = refine[0]
        fine = dict(settings) | {'nι': n, 'ns': n}
        gs = lambda mm: getattr(mm, pref).GS['PEE']['gridSettings']
        ref.append({'ρ': ρ, 'grid': 'working', 'nι': gs(m)['nι'], 'ns': gs(m)['ns'] if pref == 'CRRA' else np.nan,
                    'τB': τBpath, 'τR': τRpath, 'dτ': τRpath - τBpath})
        mf, _ = su.loadCalibrated(ρ, fine, pkldir = pkldir, epsPre0 = epsPre0, commonSettings = pref == 'LOG')
        bf = su.solvePEE(mf, pref)
        sf = mf.stateAtT0(bf['report'], t0, init = bf['init'])
        mfR = mf.createCopyFromt0(t0)
        su.installEps(mfR, su.universalEps(mfR))
        rf = su.solvePEE(mfR, pref, **sf)
        τb, τr = float(bf['τ'].xs(t0)), float(rf['τ'].iloc[0])
        ref.append({'ρ': ρ, 'grid': 'fine', 'nι': gs(mf)['nι'], 'ns': gs(mf)['ns'] if pref == 'CRRA' else np.nan,
                    'τB': τb, 'τR': τr, 'dτ': τr - τb})
        print('  refinement: ' + '; '.join('{} dtau {:+.4f} pp'.format(r['grid'], 100*r['dτ']) for r in ref))
    print('  rho={} done ({:.0f}s)'.format(ρ, time.time()-tStart))
    return foc, b5, hh, ref


if __name__ == '__main__':
    p = argparse.ArgumentParser(description = 'Diagnostics of the universalisation reform at t0.')
    p.add_argument('--rho', type = float, nargs = '+', required = True)
    p.add_argument('--pkldir', required = True, help = 'calibrated instances; a relative path is read from the '
                                                       'repository root (the script chdirs to its own directory)')
    p.add_argument('--epsPre0', action = 'store_true', help = 'the instances are the eps = 0 variant')
    p.add_argument('--tag', default = '_commonX', help = 'suffix of the four csv names')
    p.add_argument('--outdir', default = OUTDIR, help = 'where the csvs go (default results/shocks)')
    p.add_argument('--refine', type = int, default = 0, help = 'finer nι = ns for the refinement rows (0: none)')
    p.add_argument('--refineRho', type = float, nargs = '*', default = [])
    p.add_argument('--nι', '--niota', dest = 'nι', type = int, default = 45)
    p.add_argument('--ns', type = int, default = 45)
    p.add_argument('--interpKind', default = 'cubic', choices = ('linear', 'cubic', 'pchip'))
    p.add_argument('--smoothKnots', type = int, default = 4)
    a = p.parse_args()
    a.pkldir = a.pkldir if os.path.isabs(a.pkldir) else os.path.join(REPO, a.pkldir)
    settings = {'nι': a.nι, 'ns': a.ns, 'interpKind': a.interpKind, 'smoothKnots': a.smoothKnots}
    out = {k: [] for k in ('FOC', 'B5', 'HouseholdSavings', 'Refine')}
    for ρ in a.rho:
        foc, b5, hh, ref = runRho(ρ, settings, a.pkldir, a.epsPre0, (a.refine, a.refineRho) if a.refine else None)
        out['FOC'] += foc
        out['B5'] += b5
        out['HouseholdSavings'] += hh
        out['Refine'] += ref
    os.makedirs(a.outdir, exist_ok = True)
    for k, rows in out.items():
        if rows:
            path = os.path.join(a.outdir, 'ARG_reform{}{}.csv'.format(k, a.tag))
            pd.DataFrame(rows).to_csv(path, index = False)
            print('written: ' + os.path.relpath(path, REPO))

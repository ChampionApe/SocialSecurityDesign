r""" ModelESC: the US model with a deadweight wedge on redistributive benefits AND the leaded political
choice of theta (app:ESC). The "A+B" combination.

    m = ModelESC(pars = ..., wedge = {'spec': 'scale', 'phi': 0.5, 'p': 0.4}, **kwargs)
    m.calibrate()                       # (beta, omega) as ModelUS, with the wedge installed
    out = m.solveLeaded()               # the sequence of policy functions + the equilibrium path
    rec = m.calibrateWedge()            # p such that the leaded choice reproduces the observed design

THREE SPECS, and they differ in more than a formula (Base.wedgeA/wedgeB, Eq esc:AB):

    'size'    b^i = f(theta_t, tau_t)[theta*eta_i*h_i + (1-theta)h]*bbar,
              f = exp(-lambda*tau^k*Vtilde*(1-theta)^2/2), k = TAU_POWER     the paper's spec at k = 2: the
              Harberger loss of the flat component's implicit taxes on top of the wedge in place, elasticity
              kappa*T, lambda = 2 kappa pi (1 - pi) held constant; quadratic in the redistribution performed
              (Vtilde = sum_i gamma_i (y_i-1)^2/y_i, Base.Vtilde) and in the size of the system (tau^2).
              k = 1 is the tau-linear cost of the tag paper/tau1.
              lambda lives in the 'p' slot (db['wedgeP']); 'phi' stays in the db as a dummy (0.5) because
              every csv key carries it, and f never reads it.
    'sizePiR' as 'size' with lambda*g(pi_t) in place of lambda, pi_t = Pi_t/R_t (Base.wedgeΠR): keeps the
              factor the U.S. results hold constant; the paper's cross-section, with g = 2 pi (1 - pi). pi_t is the equilibrium Pi_t/R_t of the path being solved, taken as
              given by households and voters: every solve under it iterates pi to a fixed point (_ΠRFixedPoint).
              lambda*g(pi_t) with g the carrier db['wedgeCarrier'] (WEDGE_CARRIER, env SSD_CARRIER):
              'pi' g = pi (p is epsilon), 'convex' g = 2 pi (1 - pi) (p is kappa); pi itself, its dating and its
              fixed point are the same under both
    'scale'   b^i = f(theta)[theta*eta_i*h_i + (1-theta)h]*bbar,  f = phi+(1-phi)theta^p   the previous spec
    'flat'    b^i = [theta*eta_i*h_i + (1-theta)f(theta)h]*bbar                              MGE's variant

Under 'size' and 'scale' f cancels out of the OECD replacement-rate ratio that identifies theta, so theta is
the same number as in the exogenous-theta model. Under 'flat' it does NOT cancel: the same RR0 datum implies a
different theta for every p, so theta and the wedge are jointly identified and getTheta becomes a scalar
root. That is why the calibration below targets the DATUM (RR0, through theta*(p)) rather than a number for
theta -- see notes/crossCuttingFindings.md #9, which is exactly this trap in its earlier incarnation.

THE CALIBRATION, and what is approximate about it. Three targets, three parameters:

    R_{t0}, tau_{t0}     ->  beta, omega     ModelUS.calibrate, wedge installed, theta EXOGENOUS at theta*
    theta_{t0} = theta*  ->  p               the outer root here (lambda under 'size'; phi is imposed)

where theta_{t0} is the design in FORCE at the baseline year on a freely simulated path -- chosen in 1990,
not in 2020 (leadedDesignAtT0). So the observed design is what the model's own political history delivers
by 2020, and the counterfactual tables, which are new paths read at 2020, have a baseline row that sits on
the datum. Everything after 2020 is then a prediction: as nu_t moves the design moves with it, and that
drift is the object of the exercise.

The target is not thetaPolicy_{t0}(theta*) = theta*, the choice MADE at 2020, which is the right object
only when the design is pinned as history through 2020 and the tables are read at 2050; the two differ by
one period of nu_t.

The approximation is in the inner block: (beta, omega) are calibrated on the path with theta held at
theta*, not on the endogenous-theta path. At the calibrated p those two paths agree at t0 by construction
and drift only later, so the error is second order -- but it is an error, not an identity, so solveLeaded
reports `targetDrift`, the distance the endogenous path lands from the (R, tau) targets it was calibrated
to. Read it before trusting a result: it is the diagnostic that says whether the fixed-point approximation
held at this parameter point, and it costs nothing.
"""
import os
from copy import deepcopy
import numpy as np, pandas as pd
from scipy import optimize
from model import ModelUS
from base import CARRIERS
from policy import mergeMultiplicity
from policyESC import LeadedLOG, LeadedCRRA, LeadedCRRA2D, PermanentLOG, PermanentCRRA, SOLVE_FAILURES

WEDGE_SPECS = (None, 'scale', 'flat', 'size', 'sizePiR')
# calibrateWedge's default scan bracket per spec (a log grid): p in [0.05, 3] for the theta^p forms; lambda
# in [0.5, 200] for 'size' at rho = 1 under LOG, where the residual runs from the theta -> 0 corner at small
# lambda to the theta -> 1 corner at large, with one sign change. The CRRA driver has its own defaults
# (runESCcrra.BRACKET_CRRA, --bracket): the required cost falls steeply in rho. The 'sizePiR' lambda is the
# 'size' one over pi_t0 (about 0.55), so its bracket is doubled.
# f carries tau^TAU_POWER (base.fWedge). The calibrated lambda then scales by about tau0^(1-k),
# tau0 = 0.105229 the U.S. 2020 tax every calibration pins, so the 'size'/'sizePiR' brackets carry TAU_SCALE.
TAU_POWER = float(os.environ.get('SSD_TAUPOWER', 2.))
TAU_SCALE = 0.105229**(1. - TAU_POWER)
# The carrier g(pi) that multiplies the 'p' slot in f under 'sizePiR' (base.carrierFactor, db
# 'wedgeCarrier'): 'pi' g = pi (p is epsilon), 'convex' g = 2 pi (1 - pi) (p is kappa, the slope of the elasticity
# in the wedge; lambda_t = 2 kappa pi_t (1 - pi_t)). kappa at the U.S. is about lambda/(2 * 0.247), inside the
# 'sizePiR' brackets below at k = 1 and 2.
WEDGE_CARRIER = os.environ.get('SSD_CARRIER', 'convex')
WEDGE_BRACKET = {'scale': (0.05, 3.0), 'flat': (0.05, 3.0), 'size': (0.5*TAU_SCALE, 200.*TAU_SCALE),
                 'sizePiR': (1.0*TAU_SCALE, 400.*TAU_SCALE)}
# The pi = Pi/R path a 'sizePiR' solve starts from when it has no earlier one (U.S. 2020,
# archive/notes/reviewResponse_2026-10/theory_closedForm.md 3.4), the fixed point's tolerance (max |Δpi| over the path; the map contracts by
# about 0.03 per iteration at the U.S. calibration) and iteration cap. The exact 2-D recursion (~6 min per
# iteration at ns = 150) stops at ΠR_TOL2D: a pi error of 1e-5 moves the design by ~4e-6. A CRRA tax solve
# stops at ΠR_TOLCRRA; a pi error of 1e-8 moves tau by ~1e-11. A fixed point that stalls on a jump of the
# CRRA solve (_ΠRFixedPoint) is accepted below ΠR_STALLTOL, 2e-4 in pi (~1e-4 in the 2020 tax).
ΠR_START, ΠR_TOL, ΠR_TOLCRRA, ΠR_TOL2D, ΠR_STALLTOL, ΠR_MAXITER = 0.55, 1e-10, 1e-8, 1e-5, 2e-4, 40


class ModelESC(ModelUS):
    _defaultWedge = {'spec': None, 'phi': 0.5, 'p': 1.0}

    def __init__(self, *args, wedge = None, nθ = 41, nθCand = 121, nθCandCRRA = 13,
                 nθCandPerm = 21, nθ2D = 13, nθCand2D = 21, designRule = 'root', Ma = 5, **kwargs):
        """ wedge: {'spec': None|'scale'|'flat'|'size'|'sizePiR', 'phi': float, 'p': float} ('p' is lambda
        under 'size', ε under 'sizePiR' at carrier 'pi', κ at 'convex'). spec=None reproduces ModelUS exactly (Base.wedgeA/wedgeB are then the
        identity), which is what test_esc.py pins. designRule, Ma: LeadedCRRA2D's design layer ('root', alg
        esc:crra2D, or 'legacy') and the root layer's tabulation nodes. nθ2D: LeadedCRRA2D's θ-state grid,
        its policies cubic along it (test_gridConvergence.py checks 13 nodes against 97; finding #19). """
        self._wedge0 = self._defaultWedge | (wedge or {})
        if self._wedge0['spec'] not in WEDGE_SPECS:
            raise ValueError(f"ModelESC: unknown wedge spec {self._wedge0['spec']!r}; one of {WEDGE_SPECS}.")
        super().__init__(*args, **kwargs)
        self.ESC = LeadedLOG(self, nθ = nθ, nθCand = nθCand)
        self.ESCC = LeadedCRRA(self, nθCand = nθCandCRRA)
        self.ESCC2 = LeadedCRRA2D(self, nθ = nθ2D, nθCand = nθCand2D, designRule = designRule, Ma = Ma)
        self.ESCP = PermanentLOG(self, nθ = nθ, nθCand = nθCand)
        self.ESCPC = PermanentCRRA(self, nθCand = nθCandPerm)

    def initPars(self, pars = None):
        """ Write the wedge parameters BEFORE the parent fills db: updateAuxPars runs inside the parent's
        __init__ and re-derives theta through getθ, which is wedge-aware under 'flat'. """
        if WEDGE_CARRIER not in CARRIERS:
            raise ValueError(f'SSD_CARRIER={WEDGE_CARRIER!r}: one of {CARRIERS}.')
        self.db.update({'wedgeSpec': self._wedge0['spec'], 'wedgePhi': self._wedge0['phi'],
                        'wedgeP': self._wedge0['p'],
                        'wedgeTauPower': TAU_POWER,      # tau^k in f
                        'wedgeCarrier': WEDGE_CARRIER})  # g(pi) under 'sizePiR'
        super().initPars(pars = pars)

    # ------------------------------------------------------------------ wedge parameters
    @property
    def wedge(self):
        return {'spec': self.db['wedgeSpec'], 'phi': self.db['wedgePhi'], 'p': self.db['wedgeP']}

    # ------------------------------------------------------------------ 'sizePiR': pi_t = Pi_t/R_t
    def _hasΠR(self):
        return self.db.get('wedgeSpec') == 'sizePiR'

    def setΠR(self, π):
        """ Install pi_t (T,) as db['wedgeΠR'] with its lead and lag, which Base.wedgeΠR reads. """
        s = pd.Series(np.broadcast_to(np.asarray(π, dtype = float), (self.T,)).copy(), index = self.db['t'])
        self.db.update(self.addLeadAndLags('wedgeΠR', s))

    def ΠRPath(self, report, τ, θ):
        """ Pi_t/R_t on a solved path, (T,): Pi_t = nu_t w_t h_t/(w_{t-1} h_{t-1}), the pay-as-you-go return
        factor, and R_t = alpha (s_{t-1}/(nu_t h_t))^(alpha-1) (report['R']). The period before the horizon is
        the steady state the path starts from: s_{-2} = s_{-1} = s0 and h_{-1} as EE_report backs it out.
        τ, θ: the path's (T,) policy (its first-period values fix that steady state). """
        tF = self.B.tFirst
        τ, θ = np.asarray(τ, dtype = float), np.asarray(θ, dtype = float)
        s0 = float(report['s_'].values[0])
        init = self.initialState_solve(τ[tF], θ[tF])
        h_ = float(self.B.hFromS(s0, init['Γs'], tF))
        w_ = float(self.B.w(s0, h_, tF))
        wh = report['w'].values.astype(float)*report['h'].values.astype(float)
        Π = self.db['ν'].values.astype(float)*wh/np.append(w_*h_, wh[:-1])
        return Π/report['R'].values.astype(float)

    def _ΠRFixedPoint(self, solve, θOf, what, tol = ΠR_TOL):
        """ `solve()` under 'sizePiR' iterated to pi = ΠRPath(its own path): solve at a starting pi, recompute
        pi, until max|Δpi| < tol. The start is the fixed point the last solve of the same kind (`what`)
        reached on this instance or the one it was copied from, else db's pi, else ΠR_START: only the
        iteration count depends on it. θOf(out) -> the solve's (T,) design path. Returns the last solve,
        whose pi is within tol of its own path's, with 'ΠRiter', 'ΠRpath' (ΠRPath of the returned path, (T,))
        and 'ΠRgap' = max|ΠRpath - pi solved at|; db keeps the pi the returned solve was solved at, the
        warm start (_ΠRstart) the path's. Under every other spec, and inside calibrate's root (_ΠRfrozen),
        solve() once at db's pi.

        A CRRA solve can jump (by ~1e-4 in tau) as pi crosses a point where its selection changes, and the
        map then has no fixed point but a 2-cycle across the jump. When the gap stops halving, the update is
        damped (pi <- (pi + Pi/R)/2), which converges to the jump; if the gap then stalls below ΠR_STALLTOL
        the solve is returned as it stands, its 'ΠRgap' saying how far it is from consistent. """
        if not self._hasΠR() or getattr(self, '_ΠRfrozen', False):
            return solve()
        start = getattr(self, '_ΠRstart', {}).get(what)
        if start is not None and len(start) == self.T:
            self.setΠR(start)
        elif 'wedgeΠR' not in self.db:
            self.setΠR(ΠR_START)
        gaps, damp = [], False
        for k in range(1, ΠR_MAXITER + 1):
            out = solve()
            πIn = self.db['wedgeΠR'].values.astype(float)
            π = self.ΠRPath(out['report'], np.asarray(out['τ'], dtype = float), θOf(out))
            gap = float(np.max(np.abs(π - πIn)))
            gaps.append(gap)
            if gap < tol or (damp and len(gaps) > 6 and gap < ΠR_STALLTOL and gap > 0.9*min(gaps[-6:-1])):
                self._ΠRstart = dict(getattr(self, '_ΠRstart', {})) | {what: π.copy()}
                return out | {'ΠRiter': k, 'ΠRgap': gap, 'ΠRpath': π}
            if not damp and len(gaps) > 2 and gaps[-1] > 0.5*gaps[-2] and gaps[-2] > 0.5*gaps[-3]:
                damp = True
            self.setΠR(0.5*(πIn + π) if damp else π)
        raise RuntimeError(f'{what}: pi = Pi/R did not converge in {ΠR_MAXITER} iterations (last gap {gap:.2e}).')

    def calibrate(self, *args, **kwargs):
        """ ModelUS.calibrate; under 'sizePiR' an outer loop around it: (beta, omega) are calibrated at a
        fixed pi (_ΠRfrozen, so the root sees no feedback), pi is recomputed from the calibrated equilibrium,
        and the two alternate until max|Δpi| < ΠR_TOLCRRA (damped and accepted below ΠR_STALLTOL as in
        _ΠRFixedPoint). Returns the last calibration, with 'ΠRiter', 'ΠRpath' and 'ΠRgap' as _ΠRFixedPoint's;
        db and the warm start keep the pi it was solved at. """
        if not self._hasΠR():
            return super().calibrate(*args, **kwargs)
        start = getattr(self, '_ΠRstart', {}).get('calibrate')
        if start is not None and len(start) == self.T:
            self.setΠR(start)
        elif 'wedgeΠR' not in self.db:
            self.setΠR(ΠR_START)
        gaps, damp = [], False
        for k in range(1, ΠR_MAXITER + 1):
            self._ΠRfrozen = True
            try:
                cal = super().calibrate(*args, **kwargs)
            finally:
                self._ΠRfrozen = False
            pee = cal['report']['PEE']
            πIn = self.db['wedgeΠR'].values.astype(float)
            π = self.ΠRPath(pee['report'], np.asarray(pee['τ'], dtype = float), self.db['θ'].values)
            gap = float(np.max(np.abs(π - πIn)))
            gaps.append(gap)
            if gap < ΠR_TOLCRRA or (damp and len(gaps) > 6 and gap < ΠR_STALLTOL and gap > 0.9*min(gaps[-6:-1])):
                self._ΠRstart = dict(getattr(self, '_ΠRstart', {})) | {'calibrate': πIn.copy()}
                return cal | {'ΠRiter': k, 'ΠRgap': gap, 'ΠRpath': π}
            if not damp and len(gaps) > 2 and gaps[-1] > 0.5*gaps[-2] and gaps[-2] > 0.5*gaps[-3]:
                damp = True
            self.setΠR(0.5*(πIn + π) if damp else π)
        raise RuntimeError(f'calibrate: pi = Pi/R did not converge in {ΠR_MAXITER} iterations (last gap {gap:.2e}).')

    def _refuseΠR(self, what):
        if self._hasΠR():
            raise NotImplementedError(f"{what} does not iterate pi = Pi/R; 'sizePiR' runs through solvePEE_*, "
                                      'solveLeaded, leadedDesignAtT0 and solveLeaded2D only.')

    def solvePEE_LOG(self, θ = None, ε = None, s0 = None, solver = 'Robust', **kwargs):
        """ ModelUS.solvePEE_LOG; under 'sizePiR' at the fixed point in pi (_ΠRFixedPoint). """
        parent, θp = super().solvePEE_LOG, (self.db['θ'].values if θ is None else θ)
        return self._ΠRFixedPoint(lambda: parent(θ, ε, s0, solver, **kwargs), lambda out: θp, 'solvePEE_LOG')

    def solvePEE_CRRA(self, θ = None, ε = None, s0 = None, warmStart = True, backwardKwargs = None,
                      solveKwargs = None):
        """ ModelUS.solvePEE_CRRA; under 'sizePiR' at the fixed point in pi (_ΠRFixedPoint). """
        parent, θp = super().solvePEE_CRRA, (self.db['θ'].values if θ is None else θ)
        return self._ΠRFixedPoint(lambda: parent(θ, ε, s0, warmStart, backwardKwargs, solveKwargs),
                                  lambda out: θp, 'solvePEE_CRRA', tol = ΠR_TOLCRRA)

    def setWedge(self, spec = None, phi = None, p = None, update = True, carrier = None):
        """ Install wedge parameters and refresh the auxiliary ones. `update` re-derives theta (and Γh,
        eps, kappa) -- necessary under 'flat', where theta is a function of p, and harmless under 'scale'.
        Pass update=False only if you are about to overwrite db['θ'] yourself. `carrier`: 'sizePiR''s g(π)
        (CARRIERS; WEDGE_CARRIER where never set). """
        if carrier is not None:
            if carrier not in CARRIERS:
                raise ValueError(f'setWedge: carrier {carrier!r}, one of {CARRIERS}.')
            self.db['wedgeCarrier'] = carrier
        if spec is not None:
            self.db['wedgeSpec'] = None if spec == 'none' else spec
        if phi is not None:
            self.db['wedgePhi'] = float(phi)
        if p is not None:
            self.db['wedgeP'] = float(p)
        if update:
            self.updateAuxPars()
        return self.wedge

    # ------------------------------------------------------------------ theta identification
    def getθ(self):
        """ theta from the OECD replacement-rate ratio (eq:calibration:theta), wedge-aware.

        The datum pins the RATIO of the benefit formula's two coefficients, B/A = (1-RR0)/(RR0*h1-h2),
        and nothing else -- so the parent's closed form already delivers it as (1-theta)/theta. Under
        'scale' A/B carry a common f(theta) which cancels there, and theta is unchanged. Under 'flat'
        B/A = (1-theta)f(theta)/theta and theta is the root of that in (0,1).

        The root is bracketed by a scan rather than handed straight to brentq: (1-x)/x is decreasing and
        f(x) increasing, so the product need not be monotone, and a silently-picked second root would be
        a different calibration wearing the same number.

        db['θ0'] (ModelUS.θOverride) is the design itself and is returned under every spec but 'flat', which
        refuses it: there the datum is B/A, not θ. """
        θ = super().getθ()
        if self.db.get('wedgeSpec') != 'flat':
            return θ
        if self.θOverride() is not None:
            raise NotImplementedError("getθ: db['θ0'] imposes θ, but under 'flat' the datum is B/A = "
                                      "(1-θ)f(θ)/θ; the override is defined for the proportional costs only.")
        if not (0. < θ < 1.):
            raise ValueError(f'getθ: the no-wedge inversion gives θ={θ}, outside (0,1); B/A is not identified.')
        r = (1-θ)/θ
        g = lambda x: (1-x)*self.B.fWedge(x)/x - r
        grid = np.linspace(1e-6, 1-1e-9, 2001)
        v = g(grid)
        sign = np.where(np.diff(np.sign(v)) != 0)[0]
        if len(sign) == 0:
            raise ValueError(f'getθ: no root of (1-θ)f(θ)/θ = {r:.6g} on (0,1) at {self.wedge}.')
        if len(sign) > 1:
            raise ValueError(f'getθ: {len(sign)} roots of (1-θ)f(θ)/θ = {r:.6g} on (0,1) at {self.wedge}; '
                             'the wedge shape makes θ ambiguous at this (φ,p).')
        k = sign[0]
        return float(optimize.brentq(g, grid[k], grid[k+1], xtol = 1e-14))

    # ------------------------------------------------------------------ the leaded solve
    @property
    def t0Year(self):
        """ The db['t'] label of the calibration baseline year (db['t0'] is its POSITION -- see
        ModelUS.default0DParams). """
        return self.db['t'][self.db['t0']]

    def solveLeaded(self, s0 = None, pinAtT0 = True, sols = None, ε = None):
        """ Solve the leaded politico-economic equilibrium end to end.

        pinAtT0: hold theta at its calibrated value up to and including t0, so the political choice binds
        from t0 onward (LeadedLOG.simulate's tPin) -- the timing calibrateWedge targets. False lets the
        choice bind from the first period of the horizon, which makes theta_{t0} a prediction rather than
        the datum, and is only what you want for a diagnostic.

        Returns {'sols', 'θ', 'τ', 'sol', 'report', 'targetDrift'} -- the same keys solvePEE_LOG returns
        plus the design path and the fixed-point diagnostic (see the module docstring). Under 'sizePiR' the
        whole solve iterates pi = Pi/R (_ΠRFixedPoint), so `sols` cannot be passed in. """
        if self._hasΠR():
            if sols is not None:
                raise ValueError("solveLeaded: under 'sizePiR' the policies depend on pi; do not pass sols.")
            return self._ΠRFixedPoint(lambda: self._solveLeaded(s0, pinAtT0, None, ε),
                                      lambda out: out['θ'].values, f'solveLeaded{pinAtT0}')
        return self._solveLeaded(s0, pinAtT0, sols, ε)

    def _solveLeaded(self, s0, pinAtT0, sols, ε):
        sols = self.ESC.solveBackward() if sols is None else sols
        t0 = self.t0Year
        θ0 = float(self.db['θ'].xs(t0))
        θPath, τPath = self.ESC.simulate(sols, θ0, tPin = t0 if pinAtT0 else None)
        ε = self.db['eps'].values if ε is None else ε
        θ, τ = θPath.values, τPath.values
        if s0 is None:
            s0 = self.steadyState_LOG_solve(τ[self.B.tFirst], θ[self.B.tFirst], t = self.B.tFirst)['s']
        sol = self.EE_LOG_solve(τ, θ, ε, s0)
        report = self.EE_report(sol, τ, θ, ε, s0)
        drift = {'R': float(report['R'].xs(t0)) - float(self.db['R0']),
                 'τ': float(τPath.xs(t0)) - float(self.db['τ0'])}
        return {'sols': sols, 'θ': θPath, 'τ': τPath, 'sol': sol, 'report': report,
                'targetDrift': drift, 'multiplicity': LeadedLOG.multiplicity(sols)}

    def leadedChoiceAtT0(self, sols = None):
        """ thetaPolicy_{t0}(theta*) -- the design the t0 electorate picks FOR t0+1, given theta* in force.

        This is not what calibrateWedge targets any more; leadedDesignAtT0 is. Kept because it is the
        object the pinned timing reports and what test_esc.py's separability checks are written against. """
        sols = self.ESC.solveBackward() if sols is None else sols
        return self.ESC.choiceAt(sols, self.t0Year, float(self.db['θ'].xs(self.t0Year)))

    def leadedDesignAtT0(self, sols = None):
        """ theta_{t0} on the FREELY simulated path -- the design in FORCE at the baseline year when the
        political choice binds from the first period of the horizon. What calibrateWedge targets, and
        what every counterfactual table reads.

        The distinction from leadedChoiceAtT0 is one period and it matters. theta_t is a state chosen at
        t-1, so the design in force in 2020 was picked in 1990 under nu_{1990}; the choice made in 2020
        binds in 2050. They coincide only if the policy function is time-invariant, and nu_t is not. At
        p = 0.4022 (scale, phi = 0.5, rho = 1) the two differ by 0.011 -- small, but enough to put the
        baseline row off the observed design in a table read at 2020, which is why the target moved here.

        Under LOG this is exactly thetaPolicy_{t0-1} evaluated anywhere: the policy has no state at all
        (test_esc.py measures it), so the inherited theta_{t0-1} drops out and simulate's forward sweep is
        a formality. It is written through simulate anyway, because that identity is a LOG fact and this
        method must stay right if it stops holding. Under 'sizePiR' the free path's own pi is needed, so it
        is solveLeaded(pinAtT0 = False) at its fixed point. """
        if self._hasΠR():
            if sols is not None:
                raise ValueError("leadedDesignAtT0: under 'sizePiR' the policies depend on pi; do not pass sols.")
            return float(self.solveLeaded(pinAtT0 = False)['θ'].xs(self.t0Year))
        sols = self.ESC.solveBackward() if sols is None else sols
        θPath, _ = self.ESC.simulate(sols, float(self.db['θ'].xs(self.t0Year)), tPin = None)
        return float(θPath.xs(self.t0Year))

    def createCopyFromt0(self, t0):
        """ ModelUS.createCopyFromt0, with the leaded solver's horizon refreshed. self.ESC is deepcopied
        with everything else and its db/B/BG/BT already point at the copy's own objects, but its cached T
        would still describe the original horizon. """
        mt0 = super().createCopyFromt0(t0)
        mt0.ESC.T = mt0.ESCC.T = mt0.ESCC2.T = mt0.ESCP.T = mt0.ESCPC.T = mt0.T
        return mt0

    def _leadedGridChoice_CRRA(self, votePos, s0 = None, solveKwargs = None, name = ''):
        """ The design chosen by the electorate at position `votePos` for votePos+1, under CRRA, with the
        rest of the design path held at theta*.

        Holding theta at theta* away from the one candidate slot is not a shortcut AT THE CALIBRATED
        POINT: the target is precisely that the choice reproduces theta*, so at a converged p the held
        path IS the equilibrium path (up to the drift that time-varying nu_t induces later, which
        solveLeadedCRRA then reports). Away from the calibrated p it is a genuine approximation, and it
        only has to be good enough to locate the root.

        13 candidate solves rather than a full path iteration, which is what makes calibrating p under
        CRRA affordable at all.

        A candidate whose solve fails (policyESC.SOLVE_FAILURES) is not a maximiser; any other exception
        propagates. The last call's failures are in self.lastGridChoiceFailures, (candidate, error) pairs. """
        self._refuseΠR('_leadedGridChoice_CRRA')
        θStar = float(self.db['θ'].xs(self.t0Year))
        Ws = np.empty(self.ESCC.nθCand)
        failures = []
        for k, cand in enumerate(self.ESCC.θCand):
            θTry = np.full(self.T, θStar)
            θTry[votePos+1] = cand
            try:
                Ws[k] = self.ESCC.W(self.ESCC.solveθPath(θTry, s0 = s0, **(solveKwargs or {})), votePos)
            except SOLVE_FAILURES as e:
                Ws[k] = -np.inf
                failures.append((float(cand), f'{type(e).__name__}: {e}'))
        self.lastGridChoiceFailures = failures
        if not np.any(np.isfinite(Ws)):
            raise RuntimeError(f'{name or "_leadedGridChoice_CRRA"}: every candidate failed to solve '
                               f'({len(failures)} solver failures).')
        return self.ESCC._argmax(self.ESCC.θCand, Ws)[0]

    def leadedChoiceAtT0_CRRA(self, s0 = None, solveKwargs = None):
        """ theta_{t0+1} chosen at t0 under CRRA -- the CRRA counterpart of leadedChoiceAtT0. """
        return self._leadedGridChoice_CRRA(self.db['t0'], s0, solveKwargs, 'leadedChoiceAtT0_CRRA')

    def leadedDesignAtT0_CRRA(self, s0 = None, solveKwargs = None):
        """ theta_{t0} under CRRA: the design in FORCE at the baseline year, chosen one period earlier.
        The CRRA counterpart of leadedDesignAtT0, and what calibrateWedge targets at rho != 1.

        Same 13-candidate cost as leadedChoiceAtT0_CRRA -- only the slot moves, from t0+1 to t0, and the
        objective from the t0 electorate's to the t0-1 one's. A full free-path iteration inside the outer
        root on p would be the exact object and is not affordable; this is the same approximation the
        calibration has always made, relocated one period. """
        pos = self.db['t0']
        if pos == 0:
            raise ValueError('leadedDesignAtT0_CRRA: t0 is the first period of the horizon, so the design '
                             'in force there is inherited, not chosen.')
        return self._leadedGridChoice_CRRA(pos-1, s0, solveKwargs, 'leadedDesignAtT0_CRRA')

    def solveLeaded2D(self, s0 = None, pinAtT0 = True, sols = None, verbose = False):
        """ The TRUE leaded CRRA equilibrium (LeadedCRRA2D): backward-solved 2-D policy functions, the
        simulated design/tax path, the exact economic equilibrium at that path, and the same targetDrift
        diagnostic the other solvers report. pinAtT0 as solveLeaded: the design is history up to and
        including t0 -- and here the pinning has to reach the RECURSION, not just the simulation, because
        τ_t responds to θ_{t+1} under CRRA (see LeadedCRRA2D's docstring). 'multiplicity': the tax rule's
        and the design layer's counts over the recursion (LeadedCRRA2D.multiplicity; -1 for a set no period
        counted). Under 'sizePiR' the recursion iterates pi = Pi/R (_ΠRFixedPoint), so `sols` cannot be
        passed in. """
        if self._hasΠR():
            if sols is not None:
                raise ValueError("solveLeaded2D: under 'sizePiR' the policies depend on pi; do not pass sols.")
            return self._ΠRFixedPoint(lambda: self._solveLeaded2D(s0, pinAtT0, None, verbose),
                                      lambda out: out['θ'].values, f'solveLeaded2D{pinAtT0}', tol = ΠR_TOL2D)
        return self._solveLeaded2D(s0, pinAtT0, sols, verbose)

    def _solveLeaded2D(self, s0, pinAtT0, sols, verbose):
        θStar = float(self.db['θ'].xs(self.t0Year))
        pinPos = int(self.db['t0']) if pinAtT0 else None
        if sols is None:
            sols = self.ESCC2.solvePolicies(θStar, pinPos = pinPos, verbose = verbose)
        if s0 is None:
            s0 = self.ESCC2.s0FixedPoint(sols, θStar)
        θPath, τPath, sP, hP, ΓsP = self.ESCC2.simulate(sols, θStar, s0, pinPos = pinPos)
        ε = self.db['eps'].values.astype(float)
        sol = self.EE_CRRA_solve(τPath.values, θPath.values, ε, s0,
                                 x0 = np.concatenate([ΓsP, hP, sP]))
        report = self.EE_report(sol, τPath.values, θPath.values, ε, s0)
        t0 = self.t0Year
        drift = {'R': float(report['R'].xs(t0)) - float(self.db['R0']),
                 'τ': float(τPath.xs(t0)) - float(self.db['τ0'])}
        return {'sols': sols, 'θ': θPath, 'τ': τPath, 'sol': sol, 'report': report,
                'targetDrift': drift, 's0': s0, 'multiplicity': LeadedCRRA2D.multiplicity(sols)}

    def leadedDesignAtT0_2D(self, verbose = False, withMultiplicity = False):
        """ θ_{t0} on the FREELY simulated 2-D path (LeadedCRRA2D, every period chooses): the design in
        force at the baseline year, the exact CRRA counterpart of leadedDesignAtT0 and what calibrateWedge
        targets under preferences = 'CRRA2D'. One full recursion per call. withMultiplicity: return
        (θ_{t0}, solveLeaded2D's combined multiplicity) instead of θ_{t0} alone. """
        out = self.solveLeaded2D(pinAtT0 = False, verbose = verbose)
        θ0 = float(out['θ'].xs(self.t0Year))
        return (θ0, out['multiplicity']) if withMultiplicity else θ0

    def leadedChoiceAtT0_2D(self, out = None, verbose = False):
        """ θPolicy_{t0}(s_{t0-1}, θ*) from the 2-D solver -- the exact counterpart of what
        leadedChoiceAtT0_CRRA approximates with a held-fixed future, evaluated at the pinned path's own
        state entering t0. Pass a solveLeaded2D return to reuse its recursion. """
        if out is None:
            out = self.solveLeaded2D(pinAtT0 = True, verbose = verbose)
        t0 = self.t0Year
        θStar = float(self.db['θ'].xs(t0))
        s_ = float(out['report']['s_'].xs(t0))
        return self.ESCC2.choiceAt(out['sols'], t0, s_, θStar)

    def solveLeadedCRRA(self, s0 = None, maxIter = 6, tol = 1e-3, verbose = True, solveKwargs = None,
                        pinAtT0 = False):
        """ The full CRRA design path (LeadedCRRA.solvePath) plus the same targetDrift diagnostic
        solveLeaded reports.

        pinAtT0 as solveLeaded, and the default is False here because the counterfactuals are new paths:
        the design is chosen from the first period the horizon allows, so theta_{t0} is an equilibrium
        outcome and can be read at 2020. True holds it at theta* through t0 -- the unanticipated-reform
        reading, kept as a diagnostic and as the object leadedChoiceAtT0_CRRA is the one-shot version of.
        """
        self._refuseΠR('solveLeadedCRRA')
        θStar = float(self.db['θ'].xs(self.t0Year))
        rec = self.ESCC.solvePath(θStar, pinPos = self.db['t0'] if pinAtT0 else 0, maxIter = maxIter,
                                  tol = tol, s0 = s0, verbose = verbose, solveKwargs = solveKwargs)
        t0 = self.t0Year
        rec['targetDrift'] = {'R': float(rec['out']['report']['R'].xs(t0)) - float(self.db['R0']),
                              'τ': float(rec['out']['τ'].xs(t0)) - float(self.db['τ0'])}
        return rec

    # ------------------------------------------------------------------ the sequential choice
    def sequentialFOC(self, θGrid, out, pos = None):
        """ ∂W_t/∂θ_t of Eq (esc:seqFOC), COSTLESS (A = θ, B = 1-θ, A' = 1, B' = -1), on a solved path at
        horizon position `pos` (db['t0'] default), for every θ_t on θGrid. Returns (T_grid,).

        out: a solvePEE_LOG/CRRA return; its report supplies s_{t-1,i}/s_{t-1} (the predetermined ratio,
        'si_s' at t-1), h_t and s_{t-1}, and its τ the tax at t. Weights are LeadedBase.weights' old bloc,
        γ_{t-1,i} ω p_{t-1} μ_{t-1,i}. Under CRRA the summand carries (c_{2,t}^i)^{1-1/ρ}, evaluated at the
        candidate θ_t through Base.c2i's own bracket; at ρ = 1 that factor is 1 and the sum is the doc's
        formula. The sign over [0,1] is the object: negative everywhere means the costless sequential
        choice is the Beveridgean corner. pos = 0 has no solved t-1 ratio and raises; so does an installed
        wedge, since c2i would carry f(θ, τ) while the summand stays costless. """
        if self.db.get('wedgeSpec') is not None:
            raise NotImplementedError(f"sequentialFOC is the costless Eq (esc:seqFOC); wedge spec "
                                      f"{self.db.get('wedgeSpec')!r} is installed.")
        pos = int(self.db['t0']) if pos is None else int(pos)
        if pos == 0:
            raise ValueError('sequentialFOC: no predetermined ratio at the first period.')
        tIdx = self.db['t']
        t, tLag = tIdx[pos], tIdx[pos-1]
        B, rep = self.B, out['report']
        with B.cacheParams():
            α, ρ = float(B.get('α', t)), float(B.get('ρ', t))
            coef = (1-α)/α * float(B.get('p[t-1]', t)) * float(np.asarray(out['τ'])[pos]) \
                / float(B.get('κ[t-1]', t))
            hη = np.asarray(B.hηRatio(t, lag = '[t-1]'), dtype = float).ravel()
            si_ = rep['si_s'].xs(tLag).values.astype(float)
            old = np.asarray(self.ESC.weights(t)[0], dtype = float).ravel()
            h, s_ = float(rep['h'].xs(t)), float(rep['s_'].xs(t))
            q = 1 - 1/ρ
            θg = np.asarray(θGrid, dtype = float)
            foc = np.empty(len(θg))
            for k, θ in enumerate(θg):
                inner = si_ + coef*(θ*hη + (1-θ))                  # the costless bracket
                num = coef*(hη - 1)                                # A' y/Γh + B' with A' = 1, B' = -1
                c2 = np.asarray(B.c2i(h, s_, float(np.asarray(out['τ'])[pos]), θ, si_, t),
                                dtype = float).ravel() if q != 0 else np.ones_like(inner)
                foc[k] = float((old * c2**q * num/inner).sum())
        return foc

    # ------------------------------------------------------------------ the permanent choice
    def predeterminedSiRatio(self, base = None, preferences = 'LOG'):
        """ s_{t0-1,i}/s_{t0-1} from the INCUMBENT equilibrium -- the ratio the permanent choice would hold
        fixed if the reform were unanticipated. The default timing is the anticipated one, which pins it at
        the CHOSEN design instead (PermanentLOG.solveFixedPoint); this is the diagnostic alternative, and
        the seed the fixed point starts from.

        Read off the baseline report at the period BEFORE t0: EE_report's 'si_s' is computed at vintage t
        and gives s_{t,i}/s_t, so the row at t0-1 is exactly what c_{2,t0}^i consumes. base: a solved
        solvePEE_* return; solved here if not supplied. """
        t0 = self.t0Year
        pos = self.db['t'].get_loc(t0)
        if pos == 0:
            raise ValueError('predeterminedSiRatio: t0 is the first period of the horizon, so there is no '
                             'predetermined ratio to read -- use initialState_solve instead.')
        if base is None:
            base = getattr(self, f'solvePEE_{preferences}')()
        return base['report']['si_s'].xs(self.db['t'][pos-1]).values.astype(float)

    def θPathPermanent(self, θ):
        """ The design path a permanent reform at t0 actually is: the INCUMBENT design before t0, theta
        from t0 on. Before t0 the design is exogenous -- the vote has not happened -- so applying the
        chosen theta over the whole horizon would report a counterfactual past.

        It does not change the choice: eq (EE:si_s) at vintage t0-1 depends on date-t0 policy alone, so
        the predetermined ratio is the same either way (checked in test_esc.py). It changes the reported
        path before t0, which is a different object. """
        out = self.db['θ'].values.astype(float).copy()
        out[self.db['t0']:] = float(θ)
        return out

    def solvePermanent(self, preferences = 'LOG', base = None, pinning = 'fixedPoint', diagnostics = True,
                       s0 = None, θCand = None, verbose = False):
        """ The permanent design chosen at t0, and the equilibrium it implies.

        pinning selects what s_{t0-1,i}/s_{t0-1} is held at while W_{t0} is maximised. It is a TIMING
        assumption, not a numerical setting -- see PermanentLOG's docstring:

          'fixedPoint' (default)  the vote is ANTICIPATED, so the savings made at t0-1 were made against
                                  the design that wins: theta* = argmax W(theta ; siRatio(theta*)).
          'incumbent'             the ratio from the incumbent equilibrium -- the UNANTICIPATED reform.
                                  Identical to 'fixedPoint' wherever the choice reproduces the incumbent
                                  design, which is what calibrateWedge targets.
          'moving'                recomputed at each candidate. The WRONG object under any timing: savings
                                  are sunk at the vote, and z_t = 0 is built holding the ratio fixed, so
                                  a moving ratio is inconsistent with the tau it is paired with. Exposed
                                  so the size of the error stays on the record. LOG only.

        diagnostics: also report the choice under the other two readings ('θIncumbent', 'θMoving').

        s0 (CRRA): the state before the first period of the horizon (solvePEE_CRRA's s0) for every
        candidate's solve and for the equilibrium at the chosen design; None, the default, is each solve's
        own steady-state seed (steadyStatePEE_CRRA). Unused under LOG.

        Returns {'θ', 'atBound', 'τ', 'sol', 'report', 'W', 'θCand', 'converged', 'multiplicity', ...},
        'multiplicity' the tax candidates' counts over the choice's tax selections (every candidate design)
        and the equilibrium at the chosen design (policy.mergeMultiplicity; num_robustroot.tex). """
        self._refuseΠR('solvePermanent')
        t0 = self.t0Year
        if base is None:
            base = getattr(self, f'solvePEE_{preferences}')()
        θInc = float(self.db['θ'].xs(t0))
        siInc = self.predeterminedSiRatio(base = base)
        if preferences == 'LOG':
            P, kw = self.ESCP, {'θCand': θCand}
        else:
            P, kw = self.ESCPC, {'θCand': θCand, 's0': s0, 'verbose': verbose}
        if pinning == 'fixedPoint':
            rec = P.solveFixedPoint(t0, θInc, **kw)
        elif pinning == 'incumbent':
            rec = P.solve(t0, siInc, **kw)
        elif pinning == 'moving':
            if preferences != 'LOG':
                raise NotImplementedError("solvePermanent: pinning='moving' is LOG only.")
            rec = self.ESCP.solve(t0, None, θCand = θCand)
        else:
            raise ValueError(f"solvePermanent: unknown pinning {pinning!r}.")
        rec.setdefault('converged', True)     # vacuously so for the two single-pass pinnings
        if diagnostics:
            kwInc = (kw | {'cache': rec['cache']}) if 'cache' in rec else kw   # CRRA: reuse the solves
            rec['θIncumbent'] = P.solve(t0, siInc, **kwInc)['θ']
            rec['θMoving'] = (self.ESCP.solve(t0, None, θCand = θCand)['θ']
                              if preferences == 'LOG' else np.nan)

        # the equilibrium at the chosen permanent design: incumbent before t0, chosen from t0 on
        θ = self.θPathPermanent(rec['θ'])
        ε = self.db['eps'].values.astype(float)
        θSave = self.db['θ'].values.astype(float).copy()
        try:
            self.db.update(self.adjPar('θ', θ))
            if preferences == 'LOG':
                τ, multPath = self.ESCP.τPath(θ, withMultiplicity = True)
                s0e = self.steadyState_LOG_solve(τ[self.B.tFirst], θ[self.B.tFirst], t = self.B.tFirst)['s']
                sol = self.EE_LOG_solve(τ, θ, ε, s0e)
                report = self.EE_report(sol, τ, θ, ε, s0e)
            else:
                outC = self.solvePEE_CRRA(θ = θ, ε = ε, s0 = s0)
                τ, multPath = outC['τ'].values.astype(float), outC.get('multiplicity')
                sol, report = outC['sol'], outC['report']
        finally:
            self.db.update(self.adjPar('θ', θSave))
        return rec | {'τ': pd.Series(τ, index = self.db['t']), 'sol': sol, 'report': report, 'θPath': θ,
                      'multiplicity': mergeMultiplicity([rec.get('multiplicity'), multPath])}

    def permanentChoiceAtT0(self, preferences = 'LOG', base = None, verbose = False):
        """ The permanent choice evaluated at the INCUMBENT ratio -- what calibrateWedge's residual is
        built from, and deliberately not solvePermanent's fixed point.

        This is not a second convention. calibrateWedge asks for the p at which the electorate re-elects
        the system it has, i.e. choice(theta*) = theta*, and THAT EQUATION IS the fixed-point condition:
        a p where the incumbent-pinned choice reproduces the incumbent design is a p where the incumbent
        design IS the fixed point. So both readings have the same root, and this one reaches it with one
        pass per trial instead of an inner iteration.

        What the two do NOT share is the residual away from the root, so a scan must not mix them. Under
        CRRA every candidate starts from its own steady-state seed (solvePermanent's s0 = None). """
        t0 = self.t0Year
        if base is None:
            base = getattr(self, f'solvePEE_{preferences}')()
        siRatio_ = self.predeterminedSiRatio(base = base)
        if preferences == 'LOG':
            return self.ESCP.solve(t0, siRatio_)['θ']
        return self.ESCPC.solve(t0, siRatio_, verbose = verbose)['θ']

    # ------------------------------------------------------------------ calibrating the wedge
    def wedgeResidual(self, p, spec, phi, calKwargs = None, preferences = 'LOG'):
        """ (the equilibrium design at t0) - theta* at wedge parameter p, with (beta, omega) recalibrated
        to (R_{t0}, tau_{t0}) at that p. theta* is db['θ'] AFTER setWedge, so it moves with p under
        'flat'. preferences selects which solver evaluates the design.

        For the LEADED timings the target is the design in FORCE at t0 on a freely simulated path
        (leadedDesignAtT0), not the choice made at t0. That is the design the tables report, so it is the
        one the baseline has to reproduce; see leadedDesignAtT0 for the one-period distinction and what
        it costs. For the PERMANENT timings the choice at t0 IS the design at t0 -- the reform is dated
        t0 -- so those two entries are unchanged. """
        self.setWedge(spec = spec, phi = phi, p = float(p))
        self.calibrate(**(calKwargs or {}))
        θStar = float(self.db['θ'].xs(self.t0Year))
        design = {'LOG':      lambda: self.leadedDesignAtT0(),
                  'CRRA':     lambda: self.leadedDesignAtT0_CRRA(),
                  'CRRA2D':   lambda: self.leadedDesignAtT0_2D(),
                  'permLOG':  lambda: self.permanentChoiceAtT0('LOG'),
                  'permCRRA': lambda: self.permanentChoiceAtT0('CRRA')}[preferences]()
        return design - θStar

    def calibrateWedge(self, spec = 'scale', phi = 0.5, bracket = None, nScan = 12,
                       calKwargs = None, xtol = 1e-6, verbose = True, preferences = 'LOG',
                       beforeScan = None, beforeRefine = None):
        """ Calibrate p (lambda under 'size'; phi imposed) so the equilibrium design at t0 is the observed
        one (wedgeResidual). bracket: None picks WEDGE_BRACKET[spec]. preferences: 'LOG', 'CRRA' (path
        iteration), 'CRRA2D' (the exact 2-D recursion; the published CRRA method), 'permLOG', 'permCRRA'.
        beforeScan/beforeRefine: callables run before the scan and before the bracketed root --
        runESCcrra.py uses them to scan on a coarse savings grid and refine on the fine one (the 2-D
        choice is insensitive to that grid, its cost is not).

        Scans `bracket` on a log grid for a sign change before bracketing, because the residual is NOT
        guaranteed monotone in p and, more importantly, is FLAT AT A CORNER: wherever the choice is at
        theta = 0 or 1 the residual does not move with p at all, and a root finder handed such a bracket
        would report a spurious convergence at whichever endpoint it happened to test. The scan reports
        what it saw, so a run that fails to find an interior crossing says so rather than returning a
        number. Returns {'p', 'residual', 'scan', 'converged', 'wedge', 'θ'}. """
        bracket = WEDGE_BRACKET[spec] if bracket is None else bracket
        grid = np.exp(np.linspace(np.log(bracket[0]), np.log(bracket[1]), nScan))
        scan = []
        if beforeScan is not None:
            beforeScan()
        warm = {}                                    # node -> the warm starts it left (calibrate's x0, pi path)
        for p in grid:
            try:
                r = self.wedgeResidual(p, spec, phi, calKwargs, preferences)
                warm[len(scan)] = (deepcopy(self.x0.get('calibration')),
                                   deepcopy(getattr(self, '_ΠRstart', None)))
            except Exception as e:                   # a failed inner calibration is a datum, not a stop
                r = np.nan
                if verbose:
                    print(f'    p={p:.4f}: inner solve failed ({type(e).__name__}: {e})')
            scan.append({'p': float(p), 'residual': float(r),
                         'θ': float(self.db['θ'].xs(self.t0Year))})
            if verbose and np.isfinite(r):
                print(f"    p={p:.4f}: θ*={scan[-1]['θ']:.4f}  choice-θ* = {r:+.5f}")
        v = np.array([s['residual'] for s in scan])
        # A node whose inner calibration failed is skipped, not a wall: the sign change is looked for between
        # finite neighbours (on the scan grid the CRRA tax solve can jump by ~1e-5 and stop the (beta, omega)
        # root at a node the refine grid solves).
        fin = [k for k in range(len(v)) if np.isfinite(v[k])]
        idx = [(fin[j], fin[j+1]) for j in range(len(fin)-1) if v[fin[j]]*v[fin[j+1]] < 0]
        if not idx:
            return {'p': np.nan, 'residual': np.nan, 'scan': scan, 'converged': False,
                    'wedge': self.wedge, 'θ': np.nan,
                    'message': 'no sign change in the scanned bracket -- the choice never crosses θ*'}
        k0, k1 = idx[0]
        # The refine starts from the warm start of the bracket node nearer the root, not of the last node
        # scanned: a wide bracket's top nodes sit at the theta = 1 corner, (beta, omega) far from the root's.
        kw = k0 if abs(v[k0]) <= abs(v[k1]) else k1
        x0w, πw = warm[kw]
        if x0w is not None:
            self.x0['calibration'] = deepcopy(x0w)
        if πw is not None:
            self._ΠRstart = deepcopy(πw)
        if beforeRefine is not None:
            beforeRefine()
        f = lambda p: self.wedgeResidual(p, spec, phi, calKwargs, preferences)
        p = optimize.brentq(f, grid[k0], grid[k1], xtol = xtol)
        res = self.wedgeResidual(p, spec, phi, calKwargs, preferences)  # leaves db at the calibrated point
        return {'p': float(p), 'residual': float(res), 'scan': scan, 'converged': abs(res) < 1e-4,
                'wedge': self.wedge, 'θ': float(self.db['θ'].xs(self.t0Year)),
                'message': 'calibrated'}

r""" Endogenous system characteristics: the LEADED choice of theta (app:ESC, "Leaded choices of theta_t"),
under LOG (LeadedLOG) and CRRA (LeadedCRRA) preferences.

theta_{t+1} is chosen at t by the same probabilistic-voting objective that chooses tau_t, and is a STATE
at t+1. What makes the LOG case cheap is a property of the US model (no informal household), verified in
test_esc.py rather than assumed:

    z_t -- the FOC for tau_t -- depends on (tau_t, theta_t) ALONE.

theta_{t+1} enters z_t only through Theta_{h,t} -- as do tau_{t+1} and, under the 'size' wedge, f(theta_{t+1},
tau_{t+1}) -- which reaches the FOC only through dv20 (the old informal household), and that term carries
weight gamma_0 = 0 here. So tau_t = tauPolicy_t(theta_t) is a STATIC scalar problem, the two choices at t
are separable, and the state is one-dimensional.

Two further properties, both MEASURED (test_esc.py) rather than assumed, and both specific to LOG:

  * the choice does not depend on s_{t-1} -- ln(s_{t-1}) enters every term of W_t additively, which is the
    appendix's own normalisation ("simply assume that s_{t-1} = 1");
  * the choice does not depend on theta_t EITHER. ln(h_t), ln(ctilde_{1,t}^i) and ln(R_{t+1}) are each
    additively separable in tau_t and theta_{t+1}, so W_t = A(tau_t) + B(theta_{t+1}); theta_t reaches W_t
    only through tau_t, so it cannot move the argmax. thetaPolicy_t comes back constant across the state
    grid to machine precision.

Under CRRA both fail: the powers do not separate. LeadedCRRA therefore solves the PATH rather than a
policy function -- see its own docstring for what that costs and how the residual state-dependence is
measured rather than hoped away.

THE OBJECTIVE, in both cases:

    W_t = sum_i gamma_{t-1,i} omega p_{t-1} mu_{t-1,i} v_{2,t}^i + nu_t sum_i gamma_{t,i} mu_{t,i} v_{1,t}^i

with (model_PEE.tex) v_1^i = (1+beta_i)ln(ctilde_1^i) + beta_i ln(R_{t+1}) and v_2^i = ln(c_2^i) under LOG,
and v_1^i = (1+B_{t+1}^i)(ctilde_1^i)^{1-1/rho}/(1-1/rho), v_2^i = (c_2^i)^{1-1/rho}/(1-1/rho) under CRRA.
The weights are Base.FOC's own, so the objective and the tau-FOC cannot be weighted differently.

WHY A GRID AND NOT A FOC. The whole point of the exercise is to find out whether the choice is interior,
so the solver must be able to return a corner as a corner. Both solvers evaluate W_t on a grid of
candidates, take the argmax, and refine it parabolically when interior; a corner is reported as such.
Differentiating W_t would buy speed and cost exactly the property being measured.
"""
import time
import numpy as np, pandas as pd
from scipy import optimize
from gridsearch import roots1d, CartesianGrid, griddedInterp1D, griddedSmooth1D
from policy import CRRA, multiplicitySummary, mergeMultiplicity, DESIGN_KEYS, DESIGN_NAMES, requireZeroMass

# What a candidate's equilibrium solve raises when it fails, and the grid searches over candidate designs
# catch: non-convergence (model._checkConverged, brentq; RuntimeError), a bracket without a sign change or
# degenerate scipy/numpy input, LinAlgError included (ValueError), float overflow or division
# (ArithmeticError). Any other exception is a bug and propagates.
SOLVE_FAILURES = (RuntimeError, ValueError, ArithmeticError)


class LeadedBase:
    """ What the LOG and CRRA leaded solvers share: the political weights, the zero-mass guard, and the
    corner-preserving argmax. """

    def weights(self, t):
        """ W_t's bloc weights: old_i = gamma_{t-1,i}*omega*p_{t-1}*mu_{t-1,i}, young_i = nu_t*gamma_{t,i}*mu_{t,i}.
        Exactly Base.FOC's combination. The informal blocs are absent: gamma_0 = 0 in this model. """
        old = self.BG.get('γi[t-1]', t) * self.BG.ω2i(t)
        young = float(self.BG.get('ν', t)) * self.BG.get('γi', t) * self.BG.ω1i(t)
        return old, young

    def _requireZeroMass(self):
        """ The whole design rests on gamma_0 = 0 (see LeadedLOG.z and the module header). Checked where
        it is relied on rather than asserted in a comment. """
        requireZeroMass(self.db, 'The leaded and permanent design solvers (with mass, z_t depends on '
                                 'theta_{t+1} through dv20 and tau_t, theta_{t+1} stop being separable)')

    @staticmethod
    def _argmax(x, y):
        """ argmax of y over the grid x, refined by the parabola through the three points around it when
        the maximum is interior. Returns (x*, atBound). A corner is returned AS a corner -- whether the
        choice is interior is the object of the exercise, so it must never be interpolated away. """
        k = int(np.argmax(y))
        if k == 0 or k == len(x)-1:
            return float(x[k]), True
        y0, y1, y2 = y[k-1], y[k], y[k+1]
        denom = y0 - 2*y1 + y2
        if denom == 0 or not np.isfinite(denom):
            return float(x[k]), False
        shift = float(np.clip(0.5*(y0 - y2)/denom, -1., 1.))     # in units of the grid step
        return float(x[k] + shift*(x[k+1]-x[k])), False


class LeadedLOG(LeadedBase):
    """ Sequence of policy functions (tauPolicy_t, thetaPolicy_t) over the state theta_t, LOG case.

    The state grid is kept even though thetaPolicy_t is measurably constant along it (module docstring):
    it costs almost nothing, it is what MAKES that property measurable, and it stops being true the moment
    anything breaks the log separability. """

    def __init__(self, m, nθ = 41, nθCand = 121):
        self.m = m
        self.B, self.BG, self.BT = m.B, m.BG, m.BT
        self.db = m.db
        self.ni, self.T = m.ni, m.T
        self.nθ, self.nθCand = nθ, nθCand
        self.θGrid = np.linspace(0., 1., nθ)          # the state grid, and where policies are tabulated
        self.θCand = np.linspace(0., 1., nθCand)      # candidates searched over

    # ------------------------------------------------------------------ tau given theta
    def z(self, t, τ, θ, tLag, terminal):
        """ The LOG political FOC z_t at (tau, theta), both broadcastable to a common shape (M,). Thin
        wrapper over LOG.stateGrid/focGrid -- the SAME code path solvePEE_LOG uses, so the leaded solver
        cannot drift from the exogenous-theta one.

        (tau_{t+1}, theta_{t+1}) are passed as (tau, theta) rather than their true values, and that is
        SOUND rather than approximate: they enter stateGrid only through the LEVEL of Theta_{h,t}, which
        reaches the FOC only through dv20 -- the old informal household -- whose weight is gamma_{t-1,0},
        exactly zero in this model. That covers every wedge spec: under 'size' f(theta_{t+1}, tau_{t+1})
        rides in the same level. The placeholder keeps that term finite (0*NaN would poison the whole
        FOC, README's "zero-mass slot"). test_esc.py drives the placeholder over the whole unit square,
        under 'scale' and under 'size', and asserts z_t does not move; if the informal type is ever given
        mass, that test fails and this shortcut has to go. """
        return self.stateAndZ(t, τ, θ, tLag, terminal)[1]

    def stateAndZ(self, t, τ, θ, tLag, terminal):
        """ z as above, together with the stateGrid dict it was evaluated on (the selection rule
        re-evaluates z_t on that dict with the retirees' shares frozen). """
        τ = np.atleast_1d(np.asarray(τ, dtype = float))
        θ = np.broadcast_to(np.asarray(θ, dtype = float), τ.shape)
        d = self.m.LOG.stateGrid(τ, t, θ, tLag, terminal, τ1 = τ, θ1 = θ)
        return d, self.m.LOG.focGrid(d, t, θ, float(self.db['eps'].xs(t)), terminal)

    def τOfθ(self, t, θ, tLag, terminal, tol = 0.0, polish = True):
        """ tauPolicy_t evaluated on a vector of states theta (shape (K,)): for each, the maximiser of the
        political objective among {l, u} and the downward crossings of z_t (roots1d.selectMax, as
        LOG.solveBackward_t). Vectorised over the (theta, tau) mesh in one z() call -- z_t is elementwise
        in both arguments, so the flattened mesh costs one evaluation, not K.

        polish: brentq the selected interior crossing to solver tolerance inside one grid cell. The grid
        selection alone is only O(spacing^2) accurate (~1e-4 at the default n), and tau_{t0} is a
        CALIBRATION TARGET matched at 1e-8 -- without this the wedge calibration would chase grid noise.
        Corners are never polished: z_t != 0 there is the correct answer, not a residual.

        The candidates are tested and ranked at the savings shares each implies, held fixed along the
        grid (policy.LOG._select / roots1d.selectMaxFrozen, num_robustroot.tex); the polish re-locates
        the selected crossing on z_t itself, whose root is the equilibrium.

        Returns (tau, atBound, nMax, extra), the first three (K,) and extra = {'nCand', 'nEq',
        'fallback', 'W'}, each (K,): the candidates tested, the equilibria among them, whether the
        integral criterion had to stand in, and the selected equilibrium's own objective. """
        θ = np.atleast_1d(np.asarray(θ, dtype = float))
        τGrid = self.m.LOG.GS['PEE']['solGrids']['τ']
        K, M = θ.size, τGrid.size
        τMesh = np.repeat(τGrid, K)                    # (M*K,), τ-first: node m, state k at m*K + k
        θMesh = np.tile(θ, M)                          # (M*K,)
        d, zF = self.stateAndZ(t, τMesh, θMesh, tLag, terminal)
        zz = zF.reshape(M, K)                          # columns = the K states
        εt = float(self.db['eps'].xs(t))
        frozen = lambda cand: self.m.LOG.objectiveFrozen(cand, τGrid, d, θ, t, εt, tLag, terminal)
        sel = self.m.LOG._select(τGrid, zz, frozen, tol = tol)
        τ = np.array(sel['x'], dtype = float).reshape(K)
        atBound = np.asarray(sel['atBound'], dtype = bool).reshape(K)
        nMax = np.asarray(sel['nMax'], dtype = int).reshape(K)
        for k in range(K):
            if polish and not atBound[k]:
                τ[k] = self._polish(t, τ[k], θ[k], tLag, terminal, τGrid)
        extra = {kk: np.asarray(sel[kk]).reshape(K) for kk in ('nCand', 'nEq', 'fallback', 'W')}
        return τ, atBound, nMax, extra

    def _polish(self, t, τ, θt, tLag, terminal, τGrid):
        """ Sharpen one interior crossing to solver tolerance. Brackets by the grid cell containing it and
        widens by whole cells if the selection sat on a node; returns the unpolished value if no sign
        change can be bracketed (which selectMax's own tolerance can allow at a tangential crossing). """
        step = τGrid[1] - τGrid[0]
        f = lambda x: float(self.z(t, x, θt, tLag, terminal)[0])
        for k in (1, 2, 4):
            a, b = max(τ - k*step, τGrid[0]), min(τ + k*step, τGrid[-1])
            fa, fb = f(a), f(b)
            if np.isfinite(fa) and np.isfinite(fb) and fa*fb < 0:
                return float(optimize.brentq(f, a, b, xtol = 1e-12, rtol = 8.9e-16))
        return float(τ)

    def τAt(self, t, θt):
        """ tauPolicy_t(theta_t) at ONE state, solved rather than interpolated. simulate() uses this: the
        policy function is tabulated on self.θGrid and interpolating it back costs O(spacing^2), which at
        the default grid is ~5e-6 in tau -- invisible next to any economic effect, but tau_{t0} is a
        calibration target matched at 1e-8, so the realised path should not carry avoidable grid error. """
        tIdx = self.db['t']
        pos = tIdx.get_loc(t)
        tLag = tIdx[pos-1] if pos > 0 else self.B.tFirst
        τ = self.τOfθ(t, np.atleast_1d(θt), tLag, terminal = (t == tIdx[-1]))[0]
        return float(τ[0])

    # ------------------------------------------------------------------ the objective
    def objective(self, t, tLag, t1, τt, θt, θ1, cont, s_ = 1., siRatio_ = None):
        """ W_t over a mesh of (state theta_t, candidate theta1), all arguments flattened to (M,).

        cont: {'τ1': tau_{t+1}(theta1), 'θ2': theta_{t+2}(theta1), 'τ2': tau_{t+2}(theta2),
               'terminal1': whether t+1 is the terminal period} -- the continuation, already evaluated at
        theta1 by solveBackward. s_ = s_{t-1}; the argmax does not depend on it (module docstring), and it
        defaults to the appendix's own normalisation.

        siRatio_: s_{t-1,i}/s_{t-1}, shape (ni,) or (M,ni). None (the leaded default) recomputes it from
        the candidate (tau_t, theta_t), matching LOG.stateGrid. PermanentLOG passes it explicitly, because
        there theta_t is the object being chosen and letting the predetermined state move with it would
        fold a channel into the choice that the policy maker takes as given (see base.dlnc2i_dτ).

        Returns (W, parts) with W shape (M,). """
        BG = self.BG
        βi, βi1 = BG.get('βi', t), BG.get('βi', t1)
        τ1 = cont['τ1']

        # --- period t equilibrium at the candidate
        Γs = BG.Γs(βi, τ1, θ1, t)
        Θh = BG.Θh(τt, τ1, θ1, Γs, t)
        h = BG.h(Θh, s_, t)
        s = BG.s(BG.Θs(Θh, Γs, t), s_, t)

        # --- period t+1, needed only for R_{t+1}
        if cont['terminal1']:
            Θh1 = BG.ΘhTerminal(τ1, t1)
        else:
            Γs1 = BG.Γs(βi1, cont['τ2'], cont['θ2'], t1)
            Θh1 = BG.Θh(τ1, cont['τ2'], cont['θ2'], Γs1, t1)
        h1 = BG.h(Θh1, s, t1)
        R1 = BG.Rlead(s, h1, t)

        # --- indirect utilities. The old: c_{2,t}^i at the predetermined s_{t-1,i}/s_{t-1}, a function of
        # (tau_t, theta_t) -- NOT of theta1, so it shifts W_t's level without touching the argmax.
        # Computed exactly anyway: it costs nothing and makes W_t the actual objective.
        if siRatio_ is None:
            Γs_ = BG.Γs(BG.get('βi', tLag), τt, θt, tLag)
            siRatio_ = BG.si_s(BG.get('βi', tLag), τt, θt, Γs_, tLag)
        c2 = BG.c2i(h, s_, τt, θt, siRatio_, t)
        tc1 = BG.tildec1i(h, βi, τ1, θ1, Γs, t)

        v2 = np.log(c2)
        v1 = (1+βi)*np.log(tc1) + βi*np.log(R1)[:, None]
        old, young = self.weights(t)
        W = (old*v2).sum(axis = -1) + (young*v1).sum(axis = -1)
        return W, {'h': h, 's': s, 'R1': R1, 'Γs': Γs, 'Θh': Θh}

    # ------------------------------------------------------------------ one period
    def taxes_t(self, t, tLag, terminal, cand = True):
        """ Step 2 of alg esc:leaded at period t: tauPolicy_t solved (τOfθ, polished) on the state grid and,
        with cand, at every candidate design -- the anticipated tax function T_t(θ') that period t-1's
        choice reads at its candidates. Solved there rather than interpolated off the state grid: between
        nodes a piecewise-linear table has the secant's slope, not T_t', and the design objective is flat
        enough at its maximum for that slope to move the argmax (finding #19). Returns the period-dict
        entries 'τ', 'atBoundτ', 'nMaxτ', 'nCandτ', 'nEqτ', 'fallbackτ' (nθ,) on θGrid and 'τCand',
        'nCandτCand', 'nEqτCand', 'fallbackτCand' (nθCand,) on θCand ('τCand' None without cand). """
        K = self.nθ
        θs = np.concatenate([self.θGrid, self.θCand]) if cand else self.θGrid
        τ, atBound, nMax, ex = self.τOfθ(t, θs, tLag, terminal = terminal)
        out = {'τ': τ[:K], 'atBoundτ': atBound[:K], 'nMaxτ': nMax[:K], 'nCandτ': ex['nCand'][:K],
               'nEqτ': ex['nEq'][:K], 'fallbackτ': ex['fallback'][:K], 'τCand': None}
        if cand:
            out |= {'τCand': τ[K:], 'nCandτCand': ex['nCand'][K:], 'nEqτCand': ex['nEq'][K:],
                    'fallbackτCand': ex['fallback'][K:]}
        return out

    def solveBackward_t(self, t, tLag, t1, cont, terminal1, cand = True):
        """ One period of the recursion: tauPolicy_t on the state grid (and at the candidates, taxes_t),
        then thetaPolicy_t by maximising W_t over self.θCand at each state node.

        cont: {'τ1Cand': (nθCand,) T_{t+1} at every candidate design, solved (taxes_t of t+1), 'θPolicy1':
        callable or None (None iff t+1 is terminal), 'τPolicy2': callable theta2 -> tau_{t+2} or None}.
        theta_{t+2} = thetaPolicy_{t+1}(theta') is the same at every candidate (the choice has no state,
        module docstring), so tau_{t+2} read off t+2's table moves W_t's level alone, not its argmax. """
        tax = self.taxes_t(t, tLag, terminal = False, cand = cand)

        K, C = self.nθ, self.nθCand
        θtM = np.repeat(self.θGrid, C)                 # (K*C,)
        τtM = np.repeat(tax['τ'], C)
        θ1M = np.tile(self.θCand, K)

        τ1 = np.tile(np.asarray(cont['τ1Cand'], dtype = float), K)
        if terminal1:
            θ2 = τ2 = None
        else:
            θ2 = cont['θPolicy1'](θ1M)
            τ2 = cont['τPolicy2'](θ2)
        W, parts = self.objective(t, tLag, t1, τtM, θtM, θ1M,
                                  {'τ1': τ1, 'θ2': θ2, 'τ2': τ2, 'terminal1': terminal1})
        W = W.reshape(K, C)

        θNext = np.empty(K); atBoundθ = np.empty(K, dtype = bool)
        for k in range(K):
            θNext[k], atBoundθ[k] = self._argmaxCubic(self.θCand, W[k])
        return tax | {'θGrid': self.θGrid.copy(), 'θNext': θNext, 'atBoundθ': atBoundθ, 'W': W,
                      'stateSpread': float(θNext.max() - θNext.min())}

    @staticmethod
    def _argmaxCubic(x, y):
        """ argmax of y over the uniform grid x, refined by the cubic through four nodes when the maximum is
        interior: the maximising node, its two neighbours, and the next node on the side of the higher
        neighbour. The parabola through three (LeadedBase._argmax) misplaces the vertex by O(h^2 y'''/y''),
        the cubic by O(h^3). A corner is returned as a corner; where the four nodes do not fit on the grid
        or the cubic has no maximum within a cell of the node, the parabola's answer. Returns (x*, atBound). """
        k = int(np.argmax(y))
        n = len(x)
        if k == 0 or k == n-1:
            return float(x[k]), True
        xp, _ = LeadedBase._argmax(x, y)
        j0 = k-2 if y[k-1] >= y[k+1] else k-1
        if j0 < 0 or j0 + 4 > n or not np.all(np.isfinite(y[j0:j0+4])):
            return xp, False
        h = float(x[k+1] - x[k])
        u = (np.asarray(x[j0:j0+4], dtype = float) - x[k])/h          # nodes in units of the step
        a3, a2, a1, _ = np.polyfit(u, np.asarray(y[j0:j0+4], dtype = float) - y[k], 3)
        # the stationary points, 3 a3 v^2 + 2 a2 v + a1 = 0, in the form that stays accurate as a3 -> 0
        A, Bq, Cq = 3*a3, 2*a2, a1
        disc = Bq*Bq - 4*A*Cq
        if disc < 0:
            r = []
        else:
            q = -0.5*(Bq + np.copysign(np.sqrt(disc), Bq))
            r = ([q/A] if A != 0 else []) + ([Cq/q] if q != 0 else [])
        r = [v for v in r if -1. <= v <= 1. and 6*a3*v + 2*a2 < 0]
        return (float(x[k] + r[0]*h), False) if r else (xp, False)

    # ------------------------------------------------------------------ the recursion
    def solveBackward(self, tol = 0.0):
        """ The whole sequence of policy functions, backwards from the terminal period.

        The terminal period T-1 has no theta_T to choose and its tau uses the terminal FOC (beta = 0, as
        LOG.solveBackward_t's `terminal`). Period T-2 is the first with a theta choice, and its
        continuation reaches only into T-1, whose Theta_h is the terminal formula.

        Returns {t: report}. The anticipated tax function T_{t+1}(θ') is solved at every candidate
        (taxes_t); the policies read between state nodes (thetaPolicy_{t+1}, tauPolicy_{t+2}, and the
        simulation's thetaPolicy_t) are piecewise linear (np.interp) -- notes/crossCuttingFindings.md #4/#5:
        this module's lineage has been bitten by adaptive knots and by cubic overshoot in exactly this
        role, and a policy function about to be maximised over must not carry interpolation wiggles. """
        self._requireZeroMass()
        tIdx = self.db['t']
        sols = {}
        with self.BG.cacheParams():
            tT = tIdx[-1]
            posT = tIdx.get_loc(tT)
            tLagT = tIdx[posT-1] if posT > 0 else self.B.tFirst
            sols[tT] = {'θGrid': self.θGrid.copy(), 'θNext': None, 'atBoundθ': None, 'terminal': True} \
                | self.taxes_t(tT, tLagT, terminal = True, cand = posT > 0)

            for pos in range(len(tIdx)-2, -1, -1):
                t, t1 = tIdx[pos], tIdx[pos+1]
                tLag = tIdx[pos-1] if pos > 0 else self.B.tFirst
                terminal1 = (t1 == tIdx[-1])
                s1 = sols[t1]
                cont = {'τ1Cand': s1['τCand'],
                        'θPolicy1': None if terminal1 else self._interp(s1['θGrid'], s1['θNext']),
                        'τPolicy2': None}
                if not terminal1:
                    s2 = sols[tIdx[pos+2]]
                    cont['τPolicy2'] = self._interp(s2['θGrid'], s2['τ'])
                sols[t] = self.solveBackward_t(t, tLag, t1, cont, terminal1, cand = pos > 0) | {'terminal': False}
        self.lastMultiplicity = self.multiplicity(sols)
        return sols

    @staticmethod
    def multiplicity(sols):
        """ The tax rule's counts over a solution, policy.multiplicitySummary's six entries: the selections on
        the state grid (nEqτ, nCandτ, fallbackτ) merged with those at the candidate designs (nEqτCand, ...;
        taxes_t). """
        return mergeMultiplicity([multiplicitySummary(sols, keys = ('nEqτ', 'nCandτ', 'fallbackτ')),
                                  multiplicitySummary(sols, keys = ('nEqτCand', 'nCandτCand', 'fallbackτCand'))])

    @staticmethod
    def _interp(x, y):
        """ Piecewise-linear interpolant, clamped outside [x0, xN] (np.interp's own default). """
        return lambda q: np.interp(q, x, y)

    # ------------------------------------------------------------------ forward simulation
    def simulate(self, sols, θ0, tPin = None):
        """ The equilibrium paths of theta and tau implied by the policy functions, from an inherited
        design theta = theta0.

        theta_t is the STATE at t (chosen at t-1), so theta[pos+1] = thetaPolicy_t(theta[pos]) and
        tau[pos] = tauPolicy_t(theta[pos]).

        tPin: hold theta at theta0 for every period up to AND INCLUDING this one, letting the choice bind
        only from tPin onward -- i.e. the design is history (data) until tPin and a political outcome
        after it. That is the timing the wedge calibration targets (thetaPolicy_{tPin}(theta0) = theta0),
        and it makes the baseline and endogenous paths agree exactly up to tPin, which is what lets the
        two be compared at tPin at all. None = the choice binds from the first period. """
        tIdx = self.db['t']
        θ = np.empty(len(tIdx)); τ = np.empty(len(tIdx))
        θ[0] = θ0
        pin = None if tPin is None else tIdx.get_loc(tPin)
        for pos, t in enumerate(tIdx):
            if pin is not None and pos <= pin:
                θ[pos] = θ0
            s = sols[t]
            τ[pos] = self.τAt(t, θ[pos])
            if pos < len(tIdx)-1:
                θ[pos+1] = np.interp(θ[pos], s['θGrid'], s['θNext'])
        return pd.Series(θ, index = tIdx), pd.Series(τ, index = tIdx)

    def choiceAt(self, sols, t, θt):
        """ thetaPolicy_t(theta_t) at one state -- what the calibration of the wedge targets. """
        s = sols[t]
        return float(np.interp(θt, s['θGrid'], s['θNext']))


class LeadedCRRA(LeadedBase):
    r""" The leaded choice under CRRA, solved as a PATH rather than as a policy function.

    Under CRRA neither LOG simplification survives: v = c^{1-1/rho}/(1-1/rho) does not turn products into
    sums, so W_t is not additively separable in (tau_t, theta_{t+1}) and s_{t-1} does not drop out. The
    honest Markov object is a policy function over the two-dimensional state (s_{t-1}, theta_t), i.e. a
    2-D version of policy.py's CRRA grid. This class does something cheaper and says exactly what it
    assumes:

        Iterate on the equilibrium PATH {theta_t}. At each t, evaluate W_t over a grid of candidate
        theta_{t+1}, re-solving the WHOLE CRRA equilibrium (model.solvePEE_CRRA, unchanged) for each
        candidate while holding theta_{t+2}, ... at the current iterate. Take the argmax, sweep forward,
        repeat until the path stops moving.

    WHAT THAT ASSUMES. Holding theta_{t+2} fixed while theta_{t+1} varies is exactly right if the choice
    at t+1 does not respond to the design it inherits. Under LOG that is not an approximation at all --
    thetaPolicy is provably constant in theta_t (module docstring). Under CRRA it is an approximation, and
    stateSensitivity() measures it directly: it re-runs the t+1 choice at two different theta_{t+1} and
    reports how far the response moves. Report that number alongside any result from this class; if it is
    not small, the 2-D grid is required and this shortcut is not good enough.

    Re-solving for each candidate is also what makes the envelope logic right: tau_t is re-optimised at
    every candidate, so its own response to theta_{t+1} contributes nothing to first order, and W_t is
    evaluated at the equilibrium tau rather than at a stale one. """

    def __init__(self, m, nθCand = 13):
        self.m = m
        self.B, self.BG, self.BT = m.B, m.BG, m.BT
        self.db = m.db
        self.ni, self.T = m.ni, m.T
        self.nθCand = nθCand
        self.θCand = np.linspace(0., 1., nθCand)

    # ------------------------------------------------------------------ one equilibrium
    def solveθPath(self, θ, s0 = None, **kwargs):
        """ The CRRA politico-economic equilibrium at a GIVEN design path (model.solvePEE_CRRA, which
        already accepts a time-varying theta -- CRRA.solveBackward indexes theta[pos] per period). """
        ε = self.db['eps'].values.astype(float)
        return self.m.solvePEE_CRRA(θ = np.asarray(θ, dtype = float), ε = ε, s0 = s0, **kwargs)

    def W(self, out, pos):
        """ W_t at position `pos` of a solved equilibrium (solveθPath's return). CRRA indirect utilities,
        same blocs and weights as the LOG case. """
        t = self.db['t'][pos]
        rep = out['report']
        ρ = float(self.BG.get('ρ', t))
        q = 1 - 1/ρ
        c1 = rep['tildec1i'].xs(t).values.astype(float)
        c2 = rep['c2i'].xs(t).values.astype(float)
        B = rep['B'].xs(t).values.astype(float)
        v1 = (1+B)*c1**q/q
        v2 = c2**q/q
        old, young = self.weights(t)
        return float((old*v2).sum() + (young*v1).sum())

    # ------------------------------------------------------------------ the path iteration
    def solvePath(self, θ0, pinPos = 0, maxIter = 6, tol = 1e-4, s0 = None, verbose = True,
                  solveKwargs = None):
        """ Iterate the design path to a fixed point.

        theta[pos] for pos <= pinPos is held at theta0 (the inherited design -- history, not a choice, as
        LeadedLOG.simulate's tPin). Positions pinPos+1 .. T-1 are chosen, swept forward, Gauss-Seidel
        (each choice sees the updated earlier entries and the previous iterate's later ones). A candidate
        whose solve fails (SOLVE_FAILURES) is not a maximiser; any other exception propagates.

        Returns {'θ', 'out', 'iterations', 'converged', 'atBound', 'history', 'step', 'nFailed',
        'failures'}: the failed candidate solves, as (iteration, pos, candidate, error). """
        self._requireZeroMass()
        kw = dict(solveKwargs or {})
        T = len(self.db['t'])
        θ = np.full(T, float(θ0))
        history = [θ.copy()]
        atBound = np.zeros(T, dtype = bool)
        step = np.inf
        failures = []
        for it in range(maxIter):
            θNew = θ.copy()
            for pos in range(pinPos, T-1):
                Ws = np.empty(self.nθCand)
                for k, cand in enumerate(self.θCand):
                    θTry = θNew.copy()
                    θTry[pos+1] = cand
                    try:
                        Ws[k] = self.W(self.solveθPath(θTry, s0 = s0, **kw), pos)
                    except SOLVE_FAILURES as e:
                        Ws[k] = -np.inf        # an infeasible candidate is not a maximiser
                        failures.append((it, pos, float(cand), f'{type(e).__name__}: {e}'))
                if not np.any(np.isfinite(Ws)):
                    raise RuntimeError(f'LeadedCRRA: every candidate failed to solve at pos={pos}.')
                θNew[pos+1], atBound[pos+1] = self._argmax(self.θCand, Ws)
            step = float(np.max(np.abs(θNew - θ)))
            θ = θNew
            history.append(θ.copy())
            if verbose:
                print('    iter {}: max|dθ|={:.5f}  θ={}'.format(
                    it, step, ' '.join('{:.4f}'.format(x) for x in θ[:6])))
            if step < tol:
                break
        out = self.solveθPath(θ, s0 = s0, **kw)
        return {'θ': pd.Series(θ, index = self.db['t']), 'out': out, 'iterations': it+1,
                'converged': step < tol, 'atBound': atBound, 'history': history, 'step': step,
                'nFailed': len(failures), 'failures': failures}

    # ------------------------------------------------------------------ the assumption, measured
    def stateSensitivity(self, θ, pos, δ = 0.05, s0 = None, solveKwargs = None):
        """ How much the choice at pos+1 responds to the design it inherits -- the quantity solvePath
        assumes away (see the class docstring). Perturbs theta_{pos+1} by +/-delta, re-optimises
        theta_{pos+2} at each, and returns d(theta_{pos+2})/d(theta_{pos+1}).

        Zero under LOG by the separability argument; whatever it is under CRRA is the error term of the
        path iteration, and belongs in the write-up next to the result. A candidate whose solve fails
        (SOLVE_FAILURES) is not a maximiser, counted in 'nFailed'; any other exception propagates. """
        kw = dict(solveKwargs or {})
        θ = np.asarray(θ, dtype = float)
        got = {}
        nFailed = 0
        for s, lab in ((-δ, 'lo'), (+δ, 'hi')):
            θP = θ.copy()
            θP[pos+1] = float(np.clip(θ[pos+1] + s, 0., 1.))
            Ws = np.empty(self.nθCand)
            for k, cand in enumerate(self.θCand):
                θTry = θP.copy()
                θTry[pos+2] = cand
                try:
                    Ws[k] = self.W(self.solveθPath(θTry, s0 = s0, **kw), pos+1)
                except SOLVE_FAILURES:
                    Ws[k] = -np.inf
                    nFailed += 1
            got[lab] = (self._argmax(self.θCand, Ws)[0], θP[pos+1])
        (chLo, θLo), (chHi, θHi) = got['lo'], got['hi']
        slope = (chHi - chLo)/(θHi - θLo) if θHi != θLo else np.nan
        return {'slope': float(slope), 'choiceLo': chLo, 'choiceHi': chHi, 'θLo': θLo, 'θHi': θHi,
                'nFailed': nFailed}


class Interp2D:
    """ A policy table over the 2-D state (s_, θ): piecewise-linear and linearly EXTRAPOLATING along s
    (griddedInterp1D per θ column, NaN nodes dropped per column), and along θ, clamped to the grid (which
    spans the unit interval, so no θ query can leave it), by `kind`:
      'cubic'   the Catmull-Rom cubic through the four nodes around the query (central-difference tangents,
                C1, exact for quadratics; at an end cell the ghost node continues the end secant). The
                design objective reads the continuation at candidate designs between θ nodes and is flat at
                its maximum, so a piecewise-linear θ interpolant's kinks set the argmax (finding #19).
                Catmull-Rom rather than a monotone cubic: PCHIP's tangents vanish at every turning point of
                the nodal data, and the anticipated tax turns near the observed design. Where an outer node
                is NaN the cell falls back to linear, so a query is NaN exactly when the linear one is.
      'linear'  piecewise linear.
    bounds: (lo, hi) the values are clipped to (a cubic can overshoot next to a kink or a flat stretch at a
    bound), None for none. A column with fewer than two finite nodes evaluates to NaN. """

    KINDS = ('cubic', 'linear')

    def __init__(self, sGrid, θGrid, tab, kind = 'cubic', bounds = None):
        """ tab: (ns, nθ) values on (sGrid, θGrid). """
        if kind not in self.KINDS:
            raise ValueError(f'Interp2D: kind {kind!r}, one of {self.KINDS}.')
        self.θGrid = np.asarray(θGrid, dtype = float)
        self.kind, self.bounds = kind, bounds
        sGrid = np.asarray(sGrid, dtype = float)
        self.cols = []
        for k in range(tab.shape[1]):
            col = np.asarray(tab[:, k], dtype = float)
            ok = np.isfinite(col)
            self.cols.append(griddedInterp1D(sGrid[ok], col[ok]) if ok.sum() >= 2 else None)

    def __call__(self, s, θ):
        s, θ = np.broadcast_arrays(np.asarray(s, dtype = float), np.asarray(θ, dtype = float))
        shape = s.shape
        sf, θf = s.reshape(-1), θ.reshape(-1)
        gθ = self.θGrid
        n = len(gθ)
        j = np.clip(np.searchsorted(gθ, θf) - 1, 0, n - 2)
        w = np.clip((θf - gθ[j])/(gθ[j+1] - gθ[j]), 0., 1.)
        V = np.stack([np.full_like(sf, np.nan) if c is None else np.asarray(c(sf), dtype = float)
                      for c in self.cols], axis = 1)
        idx = np.arange(len(sf))
        p1, p2 = V[idx, j], V[idx, j+1]
        out = (1 - w)*p1 + w*p2
        if self.kind == 'cubic':
            p0 = np.where(j >= 1, V[idx, np.maximum(j - 1, 0)], 2*p1 - p2)
            p3 = np.where(j <= n - 3, V[idx, np.minimum(j + 2, n - 1)], 2*p2 - p1)
            cub = p1 + 0.5*w*((p2 - p0) + w*((2*p0 - 5*p1 + 4*p2 - p3) + w*(3*(p1 - p2) + p3 - p0)))
            out = np.where(np.isfinite(cub), cub, out)
        if self.bounds is not None:
            out = np.clip(out, *self.bounds)
        return out.reshape(shape)


# ------------------------------------------------------------------ the frozen state (alg esc:crra2D)
# The inherited distribution is the one-parameter family of eq:esc:aDef: with a discount factor common to
# all types, s_{t-1,i}/s_{t-1} = 1 + a(y_i - 1). B_t at an off-grid tax is formed from h_t interpolated
# along τ (B_t = B(s_{t-1}, h_t)), as CRRA.objectiveFrozen forms the shares of a tax candidate. The step-1
# grid is CartesianGrid(τ, s_, θ1), flat in C-order; per θ_t the root layer works on (ns_, nθ1) arrays
# (state rows, design columns).

def aOf(BG, B, τ, θt, tLag):
    """ Eq esc:aDef: a_t with s_{t-1,i}/s_{t-1} = 1 + a_t(y_i - 1), i.e. 1 minus the type-free term of
    base.si_s at vintage tLag (its second term; the first and third are proportional to y_i once B_t is
    common). B: B_t per type, (..., ni); τ, θt broadcast against B[..., 0]. Returns (...,). Raises unless
    βi is the same for every type at tLag, which is what makes the family one-parameter. """
    βi = np.asarray(BG.get('βi', tLag), dtype = float)
    if βi.size > 1 and np.ptp(βi) != 0:
        raise ValueError(f'aOf: βi differs across types at {tLag} ({βi}); the inherited distribution is '
                         'then not the one-parameter family of eq:esc:aDef.')
    Bc = np.asarray(B, dtype = float)[..., 0]
    τ = np.asarray(τ, dtype = float)
    α, p, κ = BG.get('α', tLag), BG.get('p', tLag), BG.get('κ', tLag)
    return 1 + (1/(1+Bc)) * ((1-α)/α * p*BG.wedgeB(θt, τ, tLag)/κ*τ)


def sharesFrom(BG, a, tLag):
    """ s_{t-1,i}/s_{t-1} = 1 + a(y_i - 1) (eq:esc:aDef), y_i = hηRatio at vintage tLag, the y_i of
    base.si_s. a (...,) -> (..., ni). """
    y = np.asarray(BG.hηRatio(tLag), dtype = float)
    return 1 + np.asarray(a, dtype = float)[..., None]*(y - 1)


def _alongτ(x, Y, c):
    """ Y (M, C, ...) read at c (C,) along its first axis: linear between nodes (roots1d.interpAlong),
    the node value itself where c sits on a node, so that a point on the edge of the feasible sub-grid
    does not pick up its infeasible neighbour. NaN c gives NaN. """
    x = np.asarray(x, dtype = float)
    c = np.asarray(c, dtype = float)
    out = roots1d.interpAlong(x, Y, c[None, :])[0]
    k = np.clip(np.searchsorted(x, np.where(np.isfinite(c), c, x[0])), 0, x.size - 1)
    on = np.isfinite(c) & (x[k] == c)
    if on.any():
        out[on] = np.asarray(Y)[k[on], np.flatnonzero(on)]
    return out


def _interpRows(x, Y, q):
    """ np.interp(q[j], x, Y[j]) per row of Y (n, M); NaN where q is NaN. """
    out = np.full(Y.shape[0], np.nan)
    for j in np.flatnonzero(np.isfinite(q)):
        out[j] = np.interp(q[j], x, Y[j])
    return out


def _argmaxRows(x, Y):
    """ LeadedBase._argmax per row of Y (n, M), non-finite entries excluded, together with the value at
    the refined point (the parabola through the three nodes around the maximum when interior, the node
    value at a corner). Returns (x*, atBound, value), each (n,); NaN rows with no finite entry. """
    n, M = Y.shape
    xs, vs = np.full(n, np.nan), np.full(n, np.nan)
    ab = np.zeros(n, dtype = bool)
    Ym = np.where(np.isfinite(Y), Y, -np.inf)
    for j in range(n):
        if not np.isfinite(Y[j]).any():
            continue
        xs[j], ab[j] = LeadedBase._argmax(x, Ym[j])
        k = int(np.argmax(Ym[j]))
        if ab[j] or not 0 < k < M - 1:
            vs[j] = Ym[j, k]
            continue
        y0, y1, y2 = Ym[j, k-1], Ym[j, k], Ym[j, k+1]
        denom = y0 - 2*y1 + y2
        if denom == 0 or not np.isfinite(denom):
            vs[j] = y1
            continue
        u = (xs[j] - x[k])/(x[k+1] - x[k])
        vs[j] = y1 + 0.5*(y2 - y0)*u + 0.5*denom*u*u
    return xs, ab, vs


def _frozenCache(E, core):
    """ What frozenTaxPass needs that does not move with a or θ_t, once per period on core: the young's
    term of W_t on the step-1 grid, Σ_i young_i hatc1iPow_i/q as (nτ, ns_·nθ1), the retirees' weights,
    q = 1 - 1/ρ and the τ grid. """
    if '_frozen' not in core:
        t, g = core['t'], core['g']
        nτ, ns_, nθ1 = g.shape
        q = float(1 - 1/E.BG.get('ρ', t))
        old, young = E.weights(t)
        core['_frozen'] = {'q': q, 'old': old, 'τGrid': g.values('τ'),
                           'Wy': ((young*core['d']['hatc1iPow']).sum(axis = -1)/q).reshape(nτ, ns_*nθ1)}
    return core['_frozen']


def _cellCrossing(x, z, j, cols):
    """ The crossing of the piecewise-linear interpolant of z (n, C) next to node j, for the columns cols
    (j and cols (K,), 0 < j < n-1): in [x_{j-1}, x_j] or [x_j, x_{j+1}] wherever z changes sign there
    (the one nearer x_j if both, the left one on a tie), at roots1d.allRoots' linear crossing; x_j itself
    where neither cell does (a tangential maximum). At the maximising node of the integrated profile z
    changes sign in one of the two cells unless z_j = 0, and a crossing further out is never nearer x_j,
    so this is the root of the column nearest x_j (test_designChoicePilot.py T8). Returns (K,). """
    x = np.asarray(x, dtype = float)
    zl, z0, zr = z[j - 1, cols], z[j, cols], z[j + 1, cols]
    xl, x0, xr = x[j - 1], x[j], x[j + 1]
    with np.errstate(divide = 'ignore', invalid = 'ignore'):
        rl = xl - zl*(x0 - xl)/(z0 - zl)
        rr = x0 - z0*(xr - x0)/(zr - z0)
    inL, inR = zl*z0 < 0, z0*zr < 0
    left = inL & (~inR | (x0 - rl <= rr - x0))
    return np.where(left, rl, np.where(inR, rr, x0))


def frozenTaxPass(E, core, θt, a, fields = None):
    """ The tax best response at a frozen distribution (eq:esc:stateEq:tax, step 3 of alg esc:crra2D) at
    every (s_, θ1) of the core, for one θ_t. a: (ns_,) or (ns_, nθ1), constant along τ.

    z_t is re-evaluated on the grid with the shares at D(a) (CRRA.zAtShares) and integrated along τ per
    (s_, θ1) column (roots1d.cumtrapzColumns, feasible sub-grid). τ̂ is the corner where the profile's
    maximum is at an end node, else the crossing of the piecewise-linear interpolant of that z_t in a cell
    next to the maximising node (_cellCrossing). No equilibrium test: nothing moves along the grid at a
    frozen a. V = W_t(τ̂, θ1; a) read from the grid by the quadratic through the three nodes nearest τ̂
    (roots1d._quadAt; the node value at a corner), W_t = Σ young·hatc1iPow/q + Σ old·c_2(D(a))^q/q; Vy
    its young's part.

    fields: {name: flat (N,) or (N, ni) grid array} read at τ̂ as well. Returns {'τ', 'atBound', 'V', 'Vy',
    'h', 'B'} (ns_, nθ1), 'B' (ns_, nθ1, ni) formed from h at τ̂, plus the fields; NaN where a column has
    fewer than two feasible nodes. """
    BG, g, d, parts = E.BG, core['g'], core['d'], core['parts']
    t, tLag = core['t'], core['tLag']
    fc = _frozenCache(E, core)
    τGrid, q, old, Wy = fc['τGrid'], fc['q'], fc['old'], fc['Wy']
    nτ, ns_, nθ1 = g.shape
    C = ns_*nθ1
    a = np.asarray(a, dtype = float)
    aSC = np.broadcast_to(a[:, None] if a.ndim == 1 else a, (ns_, nθ1))
    aF = np.broadcast_to(aSC[None], (nτ, ns_, nθ1)).reshape(-1)
    D = sharesFrom(BG, aF, tLag)
    with BG.cacheParams(), np.errstate(divide = 'ignore', invalid = 'ignore'):
        z = np.asarray(E.zAtShares(d, parts, θt, D, t), dtype = float).reshape(nτ, C)
        c2 = BG.c2i(d['h'], d['s_'], d['τ'], θt, D, t)
        Wg = Wy + ((old*c2**q).sum(axis = -1)/q).reshape(nτ, C)
    P = roots1d.cumtrapzColumns(τGrid, z)
    fin = np.isfinite(z)
    usable = fin.sum(axis = 0) >= 2
    first = np.where(usable, fin.argmax(axis = 0), 0)
    last = np.where(usable, nτ - 1 - fin[::-1].argmax(axis = 0), 0)
    jmax = np.argmax(np.where(fin, P, -np.inf), axis = 0)
    corner = (jmax == first) | (jmax == last)
    τh = τGrid[jmax].astype(float)
    inner = np.flatnonzero(usable & ~corner)
    if inner.size:
        τh[inner] = _cellCrossing(τGrid, z, jmax[inner], inner)
    τh[~usable] = np.nan
    cols = np.arange(C)
    V = roots1d._quadAt(τGrid, Wg[None], τh[None])[0]
    Vy = roots1d._quadAt(τGrid, Wy[None], τh[None])[0]
    cc = usable & corner
    V[cc], Vy[cc] = Wg[jmax[cc], cols[cc]], Wy[jmax[cc], cols[cc]]
    h = _alongτ(τGrid, d['h'].reshape(nτ, C), τh)
    s_C = np.repeat(core['sGrid'], nθ1)
    with BG.cacheParams(), np.errstate(divide = 'ignore', invalid = 'ignore'):
        B = BG.B(s_C, h, tLag)
    out = {'τ': τh.reshape(ns_, nθ1), 'atBound': (corner & usable).reshape(ns_, nθ1),
           'V': V.reshape(ns_, nθ1), 'Vy': Vy.reshape(ns_, nθ1), 'h': h.reshape(ns_, nθ1),
           'B': B.reshape(ns_, nθ1, -1)}
    for name, Y in (fields or {}).items():
        Y = np.asarray(Y, dtype = float)
        out[name] = _alongτ(τGrid, Y.reshape((nτ, C) + Y.shape[1:]), τh).reshape((ns_, nθ1) + Y.shape[1:])
    return out


def _legacyKeys(shape):
    """ The tax-rule entries of the period dict where the legacy layer did not run: NaN / -1. """
    return {'τStar': np.full(shape, np.nan), 'atBoundτ': np.zeros(shape, dtype = bool),
            'nEqτ': np.full(shape, -1), 'nCandτ': np.full(shape, -1),
            'fallbackτ': np.zeros(shape, dtype = bool), 'W': np.full(shape, np.nan)}


class LeadedCRRA2D(CRRA, LeadedBase):
    r""" The TRUE leaded choice under CRRA: a sequence of policy functions over the two-dimensional state
    (s_{t-1}, θ_t), identified by backward iteration -- the "honest Markov object" LeadedCRRA's docstring
    names and does not compute. No path iteration and no held-fixed future: the choice of θ_{t+1} at t
    sees continuation policies τ_{t+1}(s_t, θ_{t+1}) and h_{t+1}(s_t, θ_{t+1}) that already embed the
    choice at t+1 responding to the design it inherits.

    STRUCTURE, per period t (mirrors policy.CRRA's alg:CRRA:grid with one extra layer):

      1. For every (τ_t, s_{t-1}, θ_{t+1}-candidate): resolve the equilibrium state s_t (the same fixed
         point as the exogenous solver -- _rootS -- with the continuation read off 2-D interpolants).
      2. z_t on that grid. Everything except the old formal households' dv2i is independent of the
         inherited design θ_t (θ_t enters only through the benefit split of the CURRENT old), so the
         grid carries no θ_t axis: dv2i alone is recomputed per θ_t node, with the numerical
         τ-derivatives (the expensive splines) shared across all of them.
      3-4. The choice (τ_t, θ_{t+1}) at every state (s_{t-1}, θ_t), by designRule:
         'root' (the default; alg esc:crra2D, the equilibrium at a state eq:esc:stateEq). The retirees'
           shares are the one-parameter family of eq:esc:aDef (aOf, sharesFrom). At a frozen a: the tax
           best response on the grid (frozenTaxPass), the design θ̂'(a) by LeadedBase._argmax over θ1 of
           V(θ'; a) = W_t(τ̂(θ'; a), θ'; a), τ̂(a) and h_t by np.interp along θ1, and the residual
           r(a) = a - a_t(τ̂(a), h_t) of eq:esc:aResidual (_bestAt). r is tabulated on Ma nodes spanning
           the consistent a_t of every grid node of the state plus one cell at each end; every sign change
           is closed to aTolBracket (aClose: 'secant', _closeSecant; 'bisection', _closeBisection) and is
           an equilibrium if |r| <= aTolResidual there, a jump otherwise. Among several equilibria the
           highest V(θ̂'; a); none: fallbackθ and the 'legacy' answer at that state.
         'legacy': τ*(s_, θ_t, θ1) by the selection rule along τ (num_robustroot.tex) at every candidate
           design, W_t there (young: Σ young·hatc1iPow/q; old: Σ old·c_2^{1-1/ρ}/q) with the retirees'
           shares consistent with each candidate's own τ* and h_t, argmax over θ1 with parabolic
           refinement, corners preserved, τ at the refined θ1 by linear interpolation. It values a
           candidate design at shares the electorate does not control; kept for comparisons, and runs the
           pinned periods (choose = False) under either rule.
      5. Tabulate τ/θNext/s/h/Γs at the chosen policies, smooth τ along s with PINNED knots
         (crossCuttingFindings #5), rebuild the tables' interpolants for period t-1.

    Steps 1-2 are _periodCore, 3-4 _choose (_chooseRoot or _chooseLegacy), 5 _handBack.

    Period-dict entries of the root layer, (nθ, ns_) unless noted: nBrθ (brackets), nEqθ (closed
    brackets), fallbackθ, nIterθ (closing iterations), aStar, Vstar, rStar, τChosen (τ̂(a*) before the
    smoothing of step 5); Vθ (nθ, ns_, nθ1) = V(θ'; a*); nPass (frozen passes in the period). Its
    tax-rule entries (τStar3, nEqτ, nCandτ, fallbackτ, W) are NaN / -1 unless a fallback ran the legacy
    layer in that period.

    The recursion is DIRECT -- terminal condition plus one backward pass -- so unlike LeadedCRRA's path
    iteration it needs no warm start and no convergence tolerance; the path iteration is kept as the
    cheap cross-check (validated against this class, not the other way around).

    PINNING (the timing the wedge calibration targets) is built into the recursion rather than only the
    simulation: for positions pos < pinPos the choice of θ_{pos+1} is FORCED to θPin (the candidate grid
    collapses to that one point), because under CRRA τ_t genuinely depends on θ_{t+1} -- pinning only at
    simulation time would evaluate τ off the pinned continuation. Under LOG the distinction vanishes,
    which is why LeadedLOG can pin in simulate() alone.

    Grid settings are borrowed from the model's exogenous solver (self.GS = m.CRRA.GS at solve time), so
    a driver that tunes m.CRRA.initGS(...) tunes this class with it. """

    DESIGN_RULES = ('root', 'legacy')
    _BEST = ('θ', 'atBound', 'V', 'τ', 'r', 'Vθ')        # _bestAt's per-state readings
    _MAXCLOSE = 200                                       # closing iterations per bracket, a guard only

    def __init__(self, m, nθ = 13, nθCand = 21, designRule = 'root', Ma = 5, aTolBracket = 1e-9,
                 aTolResidual = 1e-6, aClose = 'secant', θInterp = 'cubic', **kwargs):
        """ nθ, θInterp: the θ-state grid and the policies' interpolant along it (Interp2D's kind). The design
        objective is flat at its maximum, and a piecewise-linear continuation sets its argmax at the kinks
        (finding #19): 'cubic' at 13 nodes, checked against 97 by test_gridConvergence.py. A recursion's
        cost grows about in proportion to nθ; the cubic costs ~1.15 times the linear at a given nθ.
        designRule: 'root' or 'legacy' (class docstring). Ma (>= 4), aTolBracket, aTolResidual, aClose: the
        root layer's tabulation nodes, closing tolerances and closing method. """
        super().__init__(m, **kwargs)
        if designRule not in self.DESIGN_RULES:
            raise ValueError(f'LeadedCRRA2D: designRule {designRule!r}, one of {self.DESIGN_RULES}.')
        if Ma < 4:
            raise ValueError('LeadedCRRA2D: Ma >= 4 (one cell beyond the range at each end).')
        if aClose not in ('secant', 'bisection'):
            raise ValueError(f"LeadedCRRA2D: aClose {aClose!r}, 'secant' or 'bisection'.")
        if θInterp not in Interp2D.KINDS:
            raise ValueError(f'LeadedCRRA2D: θInterp {θInterp!r}, one of {Interp2D.KINDS}.')
        self.θInterp = θInterp
        self.nθ, self.nθCand = nθ, nθCand
        self.θGrid = np.linspace(0., 1., nθ)
        self.θCand = np.linspace(0., 1., nθCand)
        self.designRule, self.Ma, self.aClose = designRule, int(Ma), aClose
        self.aTolBracket, self.aTolResidual = float(aTolBracket), float(aTolResidual)

    # ------------------------------------------------------------------ the one genuine override
    def stateApprox_t(self, τ, s, t, θ1, solp):
        """ CRRA.stateApprox_t with a 2-D continuation: solp's τPolicy/hPolicy take (s_t, θ_{t+1}) and
        θ1 is a flat array over the mesh rather than a scalar. Everything downstream broadcasts. """
        BG = self.BG
        τ1, h1 = solp['τPolicy'](s, θ1), solp['hPolicy'](s, θ1)
        B1 = BG.B(s, h1, t)
        Γs = BG.Γs(B1, τ1, θ1, t)
        Θh = BG.Θh(τ, τ1, θ1, Γs, t)
        Θs = BG.Θs(Θh, Γs, t)
        return {'τ1': τ1, 'h1': h1, 'B1': B1, 'Γs': Γs, 'Θh': Θh, 'Θs': Θs}

    # ------------------------------------------------------------------ state fixed point on the grid
    def _solveStateGrid(self, τGrid, s_Grid, θ1Grid, sCand, t, sol1):
        """ s_t over the (τ, s_, θ1) product, flat in C-order. Θs is evaluated on (τ, sCand, θ1) once and
        broadcast across s_ (only the residual's (s_/ν)^σ factor involves it), exactly as the parent's
        solveStateApprox_t broadcasts across its 1-D state. """
        gA = CartesianGrid(τ = τGrid, s = sCand, θ1 = θ1Grid)
        Θs = gA.reshape(self.stateApprox_t(gA.flat['τ'], gA.flat['s'], t, gA.flat['θ1'], sol1)['Θs'])
        nτ, nsC, nθ1 = Θs.shape
        ns_ = len(s_Grid)
        big = np.ascontiguousarray(np.broadcast_to(Θs.transpose(1, 0, 2)[:, :, None, :],
                                                   (nsC, nτ, ns_, nθ1))).reshape(nsC, -1)
        s_flat = np.ascontiguousarray(np.broadcast_to(s_Grid[None, :, None],
                                                      (nτ, ns_, nθ1))).reshape(-1)
        return self._rootS(big, sCand, s_flat, t)

    # ------------------------------------------------------------------ economics at resolved states
    def _econCore(self, τ, s, s_, θ1, t, tLag, t1, ε, ε1, sol1):
        """ stateGrid_t minus its θ_t block (Γs_/si_s_/c2i live in _zAtθ/_Wgrid instead): everything the
        FOC and the objective need that does NOT depend on the inherited design. Flat arrays throughout. """
        BG = self.BG
        self._requireCRRA(t)
        d = self.stateApprox_t(τ, s, t, θ1, sol1)
        d['τ'], d['s'], d['s_'] = τ, s, s_
        d['h'] = BG.h(d['Θh'], s_, t)
        d['B'] = BG.B(s_, d['h'], tLag)
        d['hatc1iPow'] = BG.hatc1iPow(d['h'], d['B1'], d['τ1'], θ1, d['Γs'], t)
        d['lnhatc1i'] = BG.lnhatc1i(d['h'], d['B1'], d['τ1'], θ1, d['Γs'], t)
        d['tc20'] = BG.tildec20(d['h'], s_, ε, τ, t)
        d['tc20_1'] = BG.tildec20(d['h1'], s, ε1, d['τ1'], t1)
        return d

    def _focParts(self, d, g, t, ε):
        """ The θ_t-independent pieces of z_t, including every numerical τ-derivative -- CRRA.focParts_t,
        computed once and shared across the θ_t loop (the splines dominate the period's cost). """
        return self.focParts_t(d, g, t, ε)

    def _zAtθ(self, d, parts, θt, t, tLag):
        """ z_t at one inherited design θt: only the retirees' term moves (see the class docstring),
        rebuilt at the shares consistent with each node (CRRA.zAtShares). """
        BG, τ = self.BG, d['τ']
        Γs_ = BG.Γs(d['B'], τ, θt, tLag)
        si_s_ = BG.si_s(d['B'], τ, θt, Γs_, tLag)
        return self.zAtShares(d, parts, θt, si_s_, t)

    def _econAt(self, τF, s_F, θtF, θ1F, sCand, t, tLag, t1, ε, ε1, sol1):
        """ Re-solve the state fixed point and evaluate (W, s, h, Γs) at ARBITRARY flat points
        (τ, s_, θt, θ1) -- used at the selected τ* (step 4) and at the final smoothed tables (step 5),
        where τ no longer sits on the grid. """
        BG = self.BG
        N, nsC = τF.size, len(sCand)
        τR = np.ascontiguousarray(np.broadcast_to(τF, (nsC, N))).reshape(-1)
        θ1R = np.ascontiguousarray(np.broadcast_to(θ1F, (nsC, N))).reshape(-1)
        sR = np.ascontiguousarray(np.broadcast_to(sCand[:, None], (nsC, N))).reshape(-1)
        Θs = self.stateApprox_t(τR, sR, t, θ1R, sol1)['Θs'].reshape(nsC, N)
        sPt, nR = self._rootS(Θs, sCand, s_F, t)
        d = self._econCore(τF, sPt, s_F, θ1F, t, tLag, t1, ε, ε1, sol1)
        q = float(1 - 1/self.BG.get('ρ', t))
        old, young = self.weights(t)
        with np.errstate(divide = 'ignore', invalid = 'ignore'):
            Γs_ = BG.Γs(d['B'], τF, θtF, tLag)
            si_s_ = BG.si_s(d['B'], τF, θtF, Γs_, tLag)
            c2 = BG.c2i(d['h'], s_F, τF, θtF, si_s_, t)
            W = (young*d['hatc1iPow']).sum(axis = -1)/q + (old*c2**q).sum(axis = -1)/q
        return {'W': W, 's': sPt, 'h': d['h'], 'Γs': d['Γs'], 'nRoots': nR}

    # ------------------------------------------------------------------ one period
    _CHOICE_KEYS = ('τStar', 'atBoundτ', 'nEqτ', 'nCandτ', 'fallbackτ', 'W', 'θNext', 'atBoundθ', 'τSel')

    def solveBackward_t2D(self, sol1, t, tLag, t1, ε, ε1, sGrid, sCand, θ1Grid, choose):
        """ One period of the recursion (steps 1-5 of the class docstring): _periodCore (step 1-2's
        θ_t-free grid), _choose (steps 3-4) and _handBack (step 5). θ1Grid: the candidate grid for
        θ_{t+1}; a pinned period passes the single forced value and choose = False. Returns the period
        dict with (ns, nθ) tables and their Interp2D interpolants, plus the wall times tCore and tChoose
        of the first two parts. A subclass replaces the choice layer by overriding _choose alone; keys
        its _choose returns beyond _CHOICE_KEYS are carried into the period dict as they are. """
        with self.BG.cacheParams():
            tic = time.perf_counter()
            core = self._periodCore(sol1, t, tLag, t1, ε, ε1, sGrid, sCand, θ1Grid)
            tCore = time.perf_counter() - tic
            tic = time.perf_counter()
            ch = self._choose(core, θ1Grid, choose)
            tChoose = time.perf_counter() - tic
            extra = {'τStar3': ch['τStar'], 'W': ch['W'], 'atBoundτ': ch['atBoundτ'],
                     'atBoundθ': ch['atBoundθ'].T, 'nEqτ': ch['nEqτ'], 'nCandτ': ch['nCandτ'],
                     'fallbackτ': ch['fallbackτ']}
            extra |= {k: v for k, v in ch.items() if k not in self._CHOICE_KEYS}
            extra |= {'tCore': tCore, 'tChoose': tChoose}
            return self._handBack(core, ch['θNext'], ch['τSel'], choose, extra)

    def _periodCore(self, sol1, t, tLag, t1, ε, ε1, sGrid, sCand, θ1Grid):
        """ Step 1 and the θ_t-free part of step 2: s_t on the (τ, s_, θ1) grid (_solveStateGrid), the
        economics there (_econCore) and the share-free pieces of z_t with every numerical τ-derivative
        (_focParts). Returns {'g', 's', 'nRoots', 'd', 'parts'} and the arguments, flat over g in
        C-order (τ, s_, θ1). Neither θ_t nor the retirees' shares enter it, so every choice layer of
        the period runs on one core. """
        τGrid = self.GS['PEE']['solGrids']['τ']
        g = CartesianGrid(τ = τGrid, s_ = sGrid, θ1 = θ1Grid)
        with self.BG.cacheParams():
            s, nRoots = self._solveStateGrid(τGrid, sGrid, θ1Grid, sCand, t, sol1)
            d = self._econCore(g.flat['τ'], s, g.flat['s_'], g.flat['θ1'], t, tLag, t1, ε, ε1, sol1)
            parts = self._focParts(d, g, t, ε)
        return {'g': g, 's': s, 'nRoots': nRoots, 'd': d, 'parts': parts, 'sol1': sol1, 't': t,
                'tLag': tLag, 't1': t1, 'ε': ε, 'ε1': ε1, 'sGrid': sGrid, 'sCand': sCand,
                'θ1Grid': θ1Grid}

    def _choose(self, core, θ1Grid, choose):
        """ Steps 3-4 by designRule: _chooseRoot under 'root' when the period chooses on at least three
        candidate designs, _chooseLegacy otherwise (a pinned period under either rule). """
        if self.designRule == 'root' and choose and len(θ1Grid) >= 3:
            return self._chooseRoot(core, θ1Grid)
        return self._chooseLegacy(core, θ1Grid, choose)

    def _chooseLegacy(self, core, θ1Grid, choose):
        """ Steps 3-4 under designRule 'legacy': per θ_t node the tax at every (s_, θ1) by the selection
        rule (_selectND with the objectiveFrozen callback), W_t at that tax by re-solving the state
        (_econAt, the retirees at the shares consistent with each candidate's own τ* and h_t), and the
        design by _argmax over θ1Grid. Returns {'τStar', 'atBoundτ', 'nEqτ', 'nCandτ', 'fallbackτ', 'W'} of
        shape (nθ, ns_, nθ1) and {'θNext', 'atBoundθ', 'τSel'} of shape (nθ, ns_). A pinned period
        (choose = False) keeps θ1Grid[0] and the tax there. """
        g, d, parts = core['g'], core['d'], core['parts']
        t, tLag, t1, ε, ε1 = core['t'], core['tLag'], core['t1'], core['ε'], core['ε1']
        sGrid, sCand, sol1 = core['sGrid'], core['sCand'], core['sol1']
        nθ, ns_, nθ1 = self.nθ, len(sGrid), len(θ1Grid)
        with self.BG.cacheParams():
            τStar = np.empty((nθ, ns_, nθ1))
            atBoundτ = np.zeros((nθ, ns_, nθ1), dtype = bool)
            nEqτ = np.zeros((nθ, ns_, nθ1), dtype = int)
            nCandτ = np.zeros((nθ, ns_, nθ1), dtype = int)
            fallbackτ = np.zeros((nθ, ns_, nθ1), dtype = bool)
            for it, θt in enumerate(self.θGrid):
                frozen = lambda cand, θt = float(θt): self.objectiveFrozen(cand, g, d, parts, θt, t, tLag)
                sel = self._selectND(g, self._zAtθ(d, parts, θt, t, tLag), 'τ', frozen)
                τStar[it], atBoundτ[it] = sel['x'], sel['atBound']
                nEqτ[it], nCandτ[it], fallbackτ[it] = sel['nEq'], sel['nCand'], sel['fallback']

            # --- W at the selected τ*, per (θt, s_, θ1)
            θtF = np.repeat(self.θGrid, ns_*nθ1)
            s_F = np.tile(np.repeat(sGrid, nθ1), nθ)
            θ1F = np.tile(θ1Grid, nθ*ns_)
            at = self._econAt(τStar.reshape(-1), s_F, θtF, θ1F, sCand, t, tLag, t1, ε, ε1, sol1)
            W = at['W'].reshape(nθ, ns_, nθ1)

            # --- the choice
            θNext = np.full((nθ, ns_), np.nan)
            atBoundθ = np.zeros((nθ, ns_), dtype = bool)
            τSel = np.full((nθ, ns_), np.nan)
            if choose:
                Wm = np.where(np.isfinite(W), W, -np.inf)
                for it in range(nθ):
                    for j in range(ns_):
                        if not np.any(np.isfinite(W[it, j])):
                            continue
                        θNext[it, j], atBoundθ[it, j] = self._argmax(θ1Grid, Wm[it, j])
                        τSel[it, j] = np.interp(θNext[it, j], θ1Grid, τStar[it, j])
            else:
                θNext[:] = θ1Grid[0]
                τSel = τStar[:, :, 0].copy()
        return {'τStar': τStar, 'atBoundτ': atBoundτ, 'nEqτ': nEqτ, 'nCandτ': nCandτ,
                'fallbackτ': fallbackτ, 'W': W, 'θNext': θNext, 'atBoundθ': atBoundθ, 'τSel': τSel}

    def _bestAt(self, core, θt, a):
        """ Steps 3-4 of alg esc:crra2D at a frozen a (ns_,): frozenTaxPass, the design θ̂'(a) by
        LeadedBase._argmax over θ1 of V(θ'; a) with the parabola's value at an interior maximum, τ̂(a) and
        h_t by np.interp along θ1 at θ̂'(a), and the residual r(a) = a - a_t(τ̂(a), h_t) of
        eq:esc:aResidual. Returns {'θ', 'atBound', 'V', 'τ', 'r'} (ns_,) and 'Vθ' (ns_, nθ1). """
        BG, tLag = self.BG, core['tLag']
        x = np.asarray(core['θ1Grid'], dtype = float)
        fp = frozenTaxPass(self, core, θt, a)
        θh, ab, Vh = _argmaxRows(x, fp['V'])
        τh = _interpRows(x, fp['τ'], θh)
        hh = _interpRows(x, fp['h'], θh)
        with BG.cacheParams(), np.errstate(divide = 'ignore', invalid = 'ignore'):
            r = a - aOf(BG, BG.B(core['sGrid'], hh, tLag), τh, θt, tLag)
        return {'θ': θh, 'atBound': ab, 'V': Vh, 'τ': τh, 'r': r, 'Vθ': fp['V']}

    def _chooseRoot(self, core, θ1Grid):
        """ Steps 3-4 under designRule 'root' (alg esc:crra2D), per θ_t node and vectorised over s_: the
        range of a, r tabulated on Ma nodes (_bestAt), the brackets (sign changes between neighbouring
        nodes, and exact zeros at a node), their closing (_closeSecant or _closeBisection), the selection
        among the equilibria. Returns _chooseLegacy's keys -- the tax-rule entries NaN / -1 unless a
        fallback ran _chooseLegacy -- and the per-state report of the class docstring. """
        BG, d, g, tLag = self.BG, core['d'], core['g'], core['tLag']
        nτ, ns_, nθ1 = g.shape
        nθ, Ma = self.nθ, self.Ma
        close = self._closeSecant if self.aClose == 'secant' else self._closeBisection
        θNext, τSel, aStar, Vstar, rStar = (np.full((nθ, ns_), np.nan) for _ in range(5))
        atBoundθ = np.zeros((nθ, ns_), dtype = bool)
        nBr, nEq, nIter = (np.zeros((nθ, ns_), dtype = int) for _ in range(3))
        fallback = np.zeros((nθ, ns_), dtype = bool)
        Vθ = np.full((nθ, ns_, nθ1), np.nan)
        nPass = 0
        for it, θt in enumerate(self.θGrid):
            θt = float(θt)
            with BG.cacheParams(), np.errstate(divide = 'ignore', invalid = 'ignore'):
                aG = aOf(BG, d['B'], d['τ'], θt, tLag).reshape(nτ, ns_, nθ1)
                lo = np.nanmin(np.where(np.isfinite(aG), aG, np.inf), axis = (0, 2))
                hi = np.nanmax(np.where(np.isfinite(aG), aG, -np.inf), axis = (0, 2))
            live = np.isfinite(lo) & np.isfinite(hi)
            cell = np.where(live, np.maximum((hi - lo)/(Ma - 3), 1e-12), 1.)
            A = np.where(live, lo - cell, 1.)[None, :] + cell[None, :]*np.arange(Ma)[:, None]   # (Ma, ns_)
            tab = [self._bestAt(core, θt, A[k]) for k in range(Ma)]
            nPass += Ma
            T = {key: np.stack([np.asarray(b[key]) for b in tab]) for key in self._BEST}       # (Ma, ns_, ...)
            R = T['r']
            br = [[] for _ in range(ns_)]                 # (node lo, node hi) per bracket, per state
            for j in np.flatnonzero(live):
                for k in range(Ma):
                    if R[k, j] == 0:
                        br[j].append((k, k))
                    elif (k < Ma - 1 and np.isfinite(R[k, j]) and np.isfinite(R[k+1, j])
                          and R[k, j]*R[k+1, j] < 0):
                        br[j].append((k, k + 1))
            nBr[it] = [len(b) for b in br]
            found = [[] for _ in range(ns_)]              # (V, a, θ, atBound, τ, r, Vθ) per closed root
            for slot in range(int(nBr[it].max(initial = 0))):
                js = np.flatnonzero(nBr[it] > slot)
                kLo = np.array([br[j][slot][0] for j in js])
                kHi = np.array([br[j][slot][1] for j in js])
                res = close(core, θt, A[Ma//2], js, A[kLo, js], A[kHi, js],
                            {key: T[key][kLo, js] for key in self._BEST},
                            {key: T[key][kHi, js] for key in self._BEST})
                nPass += res['nPass']
                nIter[it, js] += res['nIter']
                for n, j in enumerate(js):
                    r_, V_ = res['r'][n], res['V'][n]
                    if np.isfinite(r_) and abs(r_) <= self.aTolResidual and np.isfinite(V_):
                        found[j].append((V_, res['a'][n], res['θ'][n], res['atBound'][n], res['τ'][n], r_,
                                         res['Vθ'][n]))
            for j in range(ns_):
                nEq[it, j] = len(found[j])
                if not found[j]:
                    fallback[it, j] = bool(live[j])
                    continue
                V_, a_, θ_, ab_, τ_, r_, Vθ_ = max(found[j], key = lambda f: f[0])
                θNext[it, j], atBoundθ[it, j], τSel[it, j] = θ_, ab_, τ_
                aStar[it, j], Vstar[it, j], rStar[it, j], Vθ[it, j] = a_, V_, r_, Vθ_
        legacy = _legacyKeys((nθ, ns_, nθ1))
        if fallback.any():
            legacy = self._chooseLegacy(core, θ1Grid, True)
            θNext[fallback], τSel[fallback] = legacy['θNext'][fallback], legacy['τSel'][fallback]
            atBoundθ[fallback] = legacy['atBoundθ'][fallback]
        return {k: legacy[k] for k in ('τStar', 'atBoundτ', 'nEqτ', 'nCandτ', 'fallbackτ', 'W')} | {
            'θNext': θNext, 'atBoundθ': atBoundθ, 'τSel': τSel, 'τChosen': τSel.copy(),
            'nBrθ': nBr, 'nEqθ': nEq, 'fallbackθ': fallback, 'nIterθ': nIter, 'aStar': aStar,
            'Vstar': Vstar, 'rStar': rStar, 'Vθ': Vθ, 'nPass': nPass}

    def _closeSecant(self, core, θt, fill, js, lo, hi, bLo, bHi):
        """ Close the brackets [lo, hi] (n,) of r (eq:esc:aResidual) at the states js by a bracketed
        secant: false position with the Illinois modification (the value kept at an end retained twice in
        a row is halved); bisection where the secant point is not finite or leaves the bracket, and where
        the bracket has not halved in three iterations; every trial at least aTolBracket/4 inside the
        bracket, so that next to the root a step straddles it. One frozen pass (_bestAt at every state,
        `fill` where no bracket is open) per iteration, until every bracket is narrower than aTolBracket;
        a NaN or exact-zero reading collapses its bracket onto the trial point. bLo, bHi: the _bestAt
        readings at the two ends. Returns the reading at the end with the smaller |r|, {'a'} and _BEST
        (n, ...), with 'nIter' (n,), the iterations each bracket was open, and 'nPass'. """
        tol, δ = self.aTolBracket, 0.25*self.aTolBracket
        lo, hi = np.array(lo, dtype = float), np.array(hi, dtype = float)
        bLo = {k: np.array(v) for k, v in bLo.items()}
        bHi = {k: np.array(v) for k, v in bHi.items()}
        fLo, fHi = bLo['r'].astype(float), bHi['r'].astype(float)
        n = lo.size
        side, stall, nIter = (np.zeros(n, dtype = int) for _ in range(3))
        ref = hi - lo
        nPass = 0
        for _ in range(self._MAXCLOSE):
            live = hi - lo >= tol
            if not live.any():
                break
            with np.errstate(divide = 'ignore', invalid = 'ignore'):
                x = hi - fHi*(hi - lo)/(fHi - fLo)
                bis = ~np.isfinite(x) | (x <= lo) | (x >= hi) | (stall >= 3)
                x = np.clip(np.where(bis, 0.5*(lo + hi), x), lo + δ, hi - δ)
            aFull = np.array(fill, dtype = float)
            aFull[js[live]] = x[live]
            b = self._bestAt(core, θt, aFull)
            nPass += 1
            nIter += live
            new = {k: np.asarray(b[k])[js] for k in self._BEST}
            rx = new['r']
            bad = live & (~np.isfinite(rx) | (rx == 0))
            toLo = live & ~bad & (np.sign(rx) == np.sign(bLo['r']))
            toHi = live & ~bad & ~toLo
            fHi = np.where(toLo & (side == -1), 0.5*fHi, fHi)       # Illinois: hi retained twice in a row
            fLo = np.where(toHi & (side == 1), 0.5*fLo, fLo)        # ... lo retained twice in a row
            fLo, fHi = np.where(toLo, rx, fLo), np.where(toHi, rx, fHi)
            lo, hi = np.where(toLo | bad, x, lo), np.where(toHi | bad, x, hi)
            for k in self._BEST:
                bLo[k][toLo | bad], bHi[k][toHi | bad] = new[k][toLo | bad], new[k][toHi | bad]
            side = np.where(toLo, -1, np.where(toHi, 1, side))
            w = hi - lo
            reset = (w <= 0.5*ref) | bis
            ref, stall = np.where(reset, w, ref), np.where(reset, 0, stall + live)
        pick = np.abs(bHi['r']) < np.abs(bLo['r'])
        out = {k: np.where(pick.reshape((-1,) + (1,)*(bLo[k].ndim - 1)), bHi[k], bLo[k]) for k in self._BEST}
        return out | {'a': np.where(pick, hi, lo), 'nIter': nIter, 'nPass': nPass}

    def _closeBisection(self, core, θt, fill, js, lo, hi, bLo, bHi):
        """ _closeSecant's brackets closed by bisection, all bisected together until every one is narrower
        than aTolBracket, then read at the midpoint (one more pass): the reference _closeSecant is checked
        against (test_designChoicePilot.py T9). Same arguments and return. """
        aLo, aHi = np.array(lo, dtype = float), np.array(hi, dtype = float)
        rLo = np.array(bLo['r'], dtype = float)
        nIter = np.zeros(aLo.size, dtype = int)
        nPass = 0
        while np.any(aHi - aLo >= self.aTolBracket):
            nIter += aHi - aLo >= self.aTolBracket
            am = 0.5*(aLo + aHi)
            aFull = np.array(fill, dtype = float)
            aFull[js] = am
            rm = self._bestAt(core, θt, aFull)['r'][js]
            nPass += 1
            same = np.sign(rm) == np.sign(rLo)
            bad = ~np.isfinite(rm) | (rm == 0)
            aLo = np.where(bad, am, np.where(same, am, aLo))
            aHi = np.where(bad, am, np.where(same, aHi, am))
            rLo = np.where(same, rm, rLo)
        am = 0.5*(aLo + aHi)
        aFull = np.array(fill, dtype = float)
        aFull[js] = am
        b = self._bestAt(core, θt, aFull)
        nPass += 1
        return {k: np.asarray(b[k])[js] for k in self._BEST} | {'a': am, 'nIter': nIter, 'nPass': nPass}

    def _handBack(self, core, θNext, τSel, choose, extra):
        """ Step 5: smooth τSel (nθ, ns_) along s per θt column (pinned knots, #5), then retabulate the
        equilibrium at the final (τ, θNext) so the tables the previous period interpolates are
        self-consistent. Returns the period dict; `extra` (the choice layer's per-state report) is
        merged into it. """
        sGrid, sCand, θ1Grid = core['sGrid'], core['sCand'], core['θ1Grid']
        t, tLag, t1, ε, ε1, sol1 = core['t'], core['tLag'], core['t1'], core['ε'], core['ε1'], core['sol1']
        nθ, ns_ = self.nθ, len(sGrid)
        with self.BG.cacheParams():
            knots = self.GS['PEE']['gridSettings']['smoothKnots']
            τTab = griddedSmooth1D(sGrid, τSel.T, s = 1e-5, knots = knots)     # (ns, nθ)
            θTab = θNext.T                                                     # (ns, nθ)
            fin = self._econAt(τTab.T.reshape(-1), np.tile(sGrid, nθ), np.repeat(self.θGrid, ns_),
                               θTab.T.reshape(-1), sCand, t, tLag, t1, ε, ε1, sol1)
            sTab, hTab = fin['s'].reshape(nθ, ns_).T, fin['h'].reshape(nθ, ns_).T
            ΓsTab = fin['Γs'].reshape(nθ, ns_).T

        return ({'sGrid': sGrid, 'θGrid': self.θGrid.copy(), 'θ1Grid': np.asarray(θ1Grid, dtype = float),
                 'τ': τTab, 'θNext': θTab, 's': sTab, 'h': hTab, 'Γs': ΓsTab}
                | extra
                | {'choose': choose, 'terminal': False,
                   'θSpread_s': np.nan if not choose else float(np.nanmax(θTab, axis = 0).max()
                                                               - np.nanmin(θTab, axis = 0).min())}
                | self._interpolants(sGrid, {'τ': τTab, 'h': hTab, 's': sTab, 'Γs': ΓsTab, 'θ': θTab}))

    def _interpolants(self, sGrid, tabs):
        """ The period's Interp2D policies over (s_, θ) along θ by self.θInterp: 'τPolicy' clipped to the tax
        grid's [l, u] and 'θPolicy' to [0, 1] (a cubic can overshoot a bound), the others unclipped. tabs:
        {'τ', 'h', and optionally 's', 'Γs', 'θ'} (ns, nθ). """
        gs = self.GS['PEE']['gridSettings']
        bounds = {'τ': (gs['l'], gs['u']), 'θ': (0., 1.)}
        return {f'{k}Policy': Interp2D(sGrid, self.θGrid, v, kind = self.θInterp, bounds = bounds.get(k))
                for k, v in tabs.items()}

    def solveTerminal2D(self, ε, sGrid, t, tLag):
        """ The terminal period has no design to choose -- θ_T is its state. One inherited (closed-form)
        CRRA.solveTerminal per θ node, stacked into 2-D tables. """
        ns_ = len(sGrid)
        τTab = np.empty((ns_, self.nθ)); hTab = np.empty((ns_, self.nθ))
        for k, θ in enumerate(self.θGrid):
            rep = self.solveTerminal(float(θ), ε, t = t, sGrid = sGrid)
            τTab[:, k], hTab[:, k] = rep['τ'].values, rep['h'].values
        return {'sGrid': sGrid, 'θGrid': self.θGrid.copy(), 'τ': τTab, 'h': hTab,
                'θNext': None, 'terminal': True, 'choose': False} | self._interpolants(sGrid, {'τ': τTab, 'h': hTab})

    # ------------------------------------------------------------------ the recursion
    def solvePolicies(self, θPin, pinPos = None, sGrid = None, verbose = False):
        """ The whole sequence of 2-D policy functions, backwards from the terminal period.

        θPin: the inherited design -- the value forced on pinned periods AND the reference the default
        state grid is sized at. pinPos: positions pos < pinPos have the choice of θ_{pos+1} forced to
        θPin (the design is history until pinPos, a political outcome from pinPos on -- LeadedCRRA
        .solvePath's own convention); None = every period chooses. Returns {t: period dict}. """
        self._requireZeroMass()
        self.GS = self.m.CRRA.GS                     # borrow the driver-tuned grid settings
        tIdx = self.db['t']
        ε = self.db['eps'].values.astype(float)
        if sGrid is None:
            sGrid = self.defaultSGrid(float(θPin), tIdx[-1], n = self.GS['PEE']['gridSettings']['ns'])
        sGrid = np.asarray(sGrid, dtype = float)

        posT = len(tIdx) - 1
        tT = tIdx[-1]
        tLagT = tIdx[posT-1] if posT > 0 else self.B.tFirst
        sols = {tT: self.solveTerminal2D(ε[posT], sGrid, tT, tLagT)}
        for pos in range(posT-1, -1, -1):
            t, t1 = tIdx[pos], tIdx[pos+1]
            tLag = tIdx[pos-1] if pos > 0 else self.B.tFirst
            choose = (pinPos is None) or (pos >= pinPos)
            θ1Grid = self.θCand if choose else np.array([float(θPin)])
            tic = time.time()
            sols[t] = self.solveBackward_t2D(sols[t1], t, tLag, t1, ε[pos], ε[pos+1],
                                             sGrid, sGrid, θ1Grid, choose)
            if verbose:
                print('    t={}: {}  [{:.1f}s]'.format(
                    t, 'choice' if choose else 'pinned', time.time() - tic))
        self.lastMultiplicity = self.multiplicity(sols)
        return sols

    @staticmethod
    def multiplicity(sols):
        """ The tax rule's counts (nEqτ, nCandτ, fallbackτ, summarised as nEqMax, nCandMax, ...) and the root
        design layer's (policy.DESIGN_KEYS as policy.DESIGN_NAMES: nEqθMax, nBrθMax, nStatesMultipleθ,
        nFallbackθ) over a solution, in one dict. A set no period counted is -1: the design counts under
        'legacy', the tax counts where every choosing period ran the root layer without a fallback and no
        period was pinned. """
        return (multiplicitySummary(sols, keys = ('nEqτ', 'nCandτ', 'fallbackτ'))
                | multiplicitySummary(sols, keys = DESIGN_KEYS, names = DESIGN_NAMES))

    # ------------------------------------------------------------------ forward simulation
    def s0FixedPoint(self, sols, θ0):
        """ The stationary state of the FIRST period's transition at the inherited design -- the CRRA2D
        counterpart of the steady-state seed the other solvers use for s0. brentq on the state grid when
        its ends bracket a sign change; otherwise iteration from the grid's middle, which raises unless it
        converges (|Δs| < 1e-12 within 200 steps) to a state on the grid: off it, sPolicy is an
        extrapolation, not a tabulated equilibrium. """
        tIdx = self.db['t']
        f = lambda s: float(sols[tIdx[0]]['sPolicy'](s, θ0)) - s
        sg = sols[tIdx[0]]['sGrid']
        a, b = float(sg[0]), float(sg[-1])
        fa, fb = f(a), f(b)
        if np.isfinite(fa) and np.isfinite(fb) and fa*fb < 0:
            return float(optimize.brentq(f, a, b, xtol = 1e-12))
        s, gap = 0.5*(a+b), np.inf
        for _ in range(200):
            sNew = float(sols[tIdx[0]]['sPolicy'](s, θ0))
            gap = abs(sNew - s)
            if gap < 1e-12:
                break
            s = sNew
        if not (gap < 1e-12 and a <= s <= b):
            raise RuntimeError(f's0FixedPoint: no sign change of sPolicy(s) - s on the state grid [{a:.6g}, '
                               f'{b:.6g}], and the iteration from its middle ended at s = {s:.6g} with '
                               f'|Δs| = {gap:.2e} (it must converge on the grid).')
        return s

    def simulate(self, sols, θ0, s0, pinPos = None):
        """ The equilibrium paths implied by the policy functions, from inherited design θ0 and initial
        state s0. Pinning here only overrides the REALISED θ; a pinned recursion (solvePolicies' pinPos)
        already tabulated its policies at the forced design, so passing the same pinPos twice is
        consistent, and passing it here alone is the LOG-style approximation. Returns (θ, τ, s, h, Γs)
        -- the last three length T-1, for EE_CRRA_solve's warm start. """
        tIdx = self.db['t']
        T = len(tIdx)
        θ = np.empty(T); τ = np.empty(T)
        sPath = np.empty(T-1); hPath = np.empty(T-1); ΓsPath = np.empty(T-1)
        θ[0] = float(θ0)
        s_ = float(s0)
        for pos, t in enumerate(tIdx):
            sol = sols[t]
            τ[pos] = float(sol['τPolicy'](s_, θ[pos]))
            if pos < T-1:
                nxt = float(sol['θPolicy'](s_, θ[pos]))
                θ[pos+1] = float(θ0) if (pinPos is not None and pos < pinPos) else nxt
                hPath[pos] = float(sol['hPolicy'](s_, θ[pos]))
                ΓsPath[pos] = float(sol['ΓsPolicy'](s_, θ[pos]))
                s_ = float(sol['sPolicy'](s_, θ[pos]))
                sPath[pos] = s_
        return (pd.Series(θ, index = tIdx), pd.Series(τ, index = tIdx), sPath, hPath, ΓsPath)

    def choiceAt(self, sols, t, s_, θt):
        """ θPolicy_t(s_, θ_t) at one state. """
        return float(sols[t]['θPolicy'](s_, θt))


class PermanentLOG(LeadedLOG):
    r""" The PERMANENT choice of theta (app:ESC, "Permanent choice of theta"), LOG case.

    At the reform date t0 the electorate chooses a theta expected to hold forever: theta_t = theta for
    every t >= t0. It therefore enters THREE channels at once, where the leaded choice enters one:

        the sequential channel   theta_{t0} re-splits the CURRENT retirees' benefits (the appendix's
                                 E_{2,t}^{i,theta}), which is what corners the sequential choice at zero;
        the leaded channel       theta_{t0+1} moves h_{t0}, hence every period-t0 quantity;
        the continuation         theta_{t0+k} moves tau_{t0+k} and the whole future path.

    TWO THINGS MAKE THIS CHEAPER THAN THE APPENDIX'S RECIPE, and both are worth stating because they are
    not obvious from the write-up, which proposes a two-dimensional grid over (tau_{t0}, theta).

    1. THE JOINT CHOICE CONCENTRATES. dW/dtau = 0 is the ordinary tau first-order condition evaluated at
       theta_t = theta -- the permanent choice adds nothing to it, since theta is not a function of tau.
       So the optimal tau at any candidate theta is just tauPolicy_{t0}(theta), already available, and what
       is left is a ONE-dimensional maximisation over theta. The 2-D grid is not wrong, only redundant.

    2. tau_t FOR t > t0 IS THE ORDINARY PEE AT CONSTANT theta. Nothing about the continuation is special:
       once theta is fixed forever, later periods face exactly the exogenous-theta problem tauPolicy
       already solves. There is no recursion to run.

    THE ONE THING THAT MUST NOT BE GOT WRONG. s_{t0-1,i}/s_{t0-1} is a PREDETERMINED state, and here theta
    enters it (through theta_{t0}) in a way it never does in the leaded problem. Savings at t0-1 are SUNK
    when the vote is taken at t0: a voter comparing two candidates does not get a different s_{t0-1,i}
    under each. So the ratio is one number, passed in as an argument, and never recomputed per candidate --
    the same error base.dlnc2i_dτ's docstring forbids for tau. Letting it move is worth 0.13 in theta, and
    is also internally inconsistent with the tau it is paired with: si_s depends on tau_t0 as well as
    theta_t0, while tau_t0 comes from z_t = 0, which is BUILT holding the ratio fixed. It would break the
    concentration argument above.

    WHICH number it is pinned at is a separate question, and it is the timing that answers it. The vote is
    ANTICIPATED: households arrive at t0 knowing a design will be chosen, so the savings they made at t0-1
    were made against the design that actually wins. Rational expectations therefore make the choice a
    FIXED POINT -- solveFixedPoint below -- and not the incumbent's ratio, which is the unanticipated
    reading. The two coincide exactly wherever the chosen design equals the incumbent one, which is what
    the wedge calibration targets, so the calibrated p is the same under both; they separate only away
    from it. ModelESC.solvePermanent exposes all three readings and reports the gaps. """

    def τPath(self, θ, tFrom = None, withMultiplicity = False):
        """ tau_t for every t >= tFrom given a design path: the ordinary PEE, period by period, since z_t
        depends on (tau_t, theta_t) alone -- so a design that is one number before t0 and another from t0
        on needs no more work than a constant one. theta: scalar (constant) or a length-T path. Periods
        before tFrom are left NaN -- they are not part of this experiment and must not be read.
        withMultiplicity: return (tau, the tax candidates' counts over the periods solved,
        policy.multiplicitySummary) instead of tau alone. Each period is τAt's τOfθ call. """
        tIdx = self.db['t']
        θPath = np.full(len(tIdx), float(θ)) if np.isscalar(θ) else np.asarray(θ, dtype = float)
        out = np.full(len(tIdx), np.nan)
        start = 0 if tFrom is None else tIdx.get_loc(tFrom)
        counts = {}
        for pos, t in enumerate(tIdx):
            if pos >= start:
                tLag = tIdx[pos-1] if pos > 0 else self.B.tFirst
                τ, _, _, counts[t] = self.τOfθ(t, np.atleast_1d(θPath[pos]), tLag, terminal = (t == tIdx[-1]))
                out[pos] = float(τ[0])
        return (out, multiplicitySummary(counts)) if withMultiplicity else out

    def siRatioAt(self, t0, θ):
        """ s_{t0-1,i}/s_{t0-1} implied by a permanent design theta in force from t0 onward, shape (ni,).

        Eq (EE:si_s) at vintage t0-1 is a function of DATE-t0 policy alone, so the kinked path the reform
        actually is -- the incumbent design before t0, theta from t0 -- gives the same ratio as a constant
        theta path, and the fixed point below stays a scalar problem. Verified against the solved baseline
        report in test_esc.py. """
        tIdx = self.db['t']
        tLag = tIdx[tIdx.get_loc(t0)-1]
        τ0, θv = np.atleast_1d(self.τAt(t0, float(θ))), np.atleast_1d(float(θ))
        with self.BG.cacheParams():
            β_ = self.BG.get('βi', tLag)
            Γs_ = self.BG.Γs(β_, τ0, θv, tLag)
            return np.asarray(self.BG.si_s(β_, τ0, θv, Γs_, tLag), dtype = float).ravel()

    def objectiveOverθ(self, t0, siRatio_, θCand = None, s_ = 1.):
        """ W_{t0} over a grid of candidate permanent designs, fully vectorised.

        siRatio_: the pinned predetermined state, shape (ni,). See the class docstring for why it is an
        argument rather than something recomputed here. 'multiplicity': the tax candidates' counts of the
        selections at t0, t0+1, t0+2 over every candidate design (τOfθ, policy.multiplicitySummary). """
        th = self.θCand if θCand is None else np.asarray(θCand, dtype = float)
        tIdx = self.db['t']
        pos = tIdx.get_loc(t0)
        tLag = tIdx[pos-1] if pos > 0 else self.B.tFirst
        t1 = tIdx[pos+1]
        terminal1 = (t1 == tIdx[-1])

        counts = {}
        τ0, _, _, counts[t0] = self.τOfθ(t0, th, tLag, terminal = False)
        τ1, _, _, counts[t1] = self.τOfθ(t1, th, t0, terminal = terminal1)
        if terminal1:
            τ2 = None
        else:
            t2 = tIdx[pos+2]
            τ2, _, _, counts[t2] = self.τOfθ(t2, th, t1, terminal = (t2 == tIdx[-1]))
        cont = {'τ1': τ1, 'θ2': th, 'τ2': τ2, 'terminal1': terminal1}
        W, parts = self.objective(t0, tLag, t1, τ0, th, th, cont, s_ = s_, siRatio_ = siRatio_)
        return {'θ': th, 'W': W, 'τ0': τ0, 'τ1': τ1, 'parts': parts, 'multiplicity': multiplicitySummary(counts)}

    def solve(self, t0, siRatio_, θCand = None, s_ = 1.):
        """ The permanent design chosen at t0, given the predetermined s_{t0-1,i}/s_{t0-1}.

        nTurning counts sign changes in the gradient of W over the candidate grid, so a multi-peaked
        objective is visible rather than silently resolved by argmax. 'multiplicity': objectiveOverθ's. """
        self._requireZeroMass()
        with self.BG.cacheParams():
            d = self.objectiveOverθ(t0, siRatio_, θCand = θCand, s_ = s_)
        thStar, atBound = self._argmax(d['θ'], d['W'])
        turns = int(np.sum(np.diff(np.sign(np.diff(d['W']))) != 0))
        return {'θ': thStar, 'atBound': atBound, 'W': d['W'], 'θCand': d['θ'],
                'τAtChoice': float(np.interp(thStar, d['θ'], d['τ0'])), 'nTurning': turns,
                'multiplicity': d['multiplicity']}

    def solveFixedPoint(self, t0, θ0, θCand = None, s_ = 1., tol = 1e-9, maxIter = 50):
        """ The permanent design under the ANTICIPATED vote (the class docstring's timing):

            theta* = argmax_theta W_{t0}( theta ; siRatio(theta*) ),

        solved as a best response on the predetermined ratio. theta0 seeds it -- the incumbent design is
        the natural seed, and wherever the wedge is calibrated it is already the fixed point, so the
        iteration costs one extra pass and confirms it.

        `converged` is REPORTED, not asserted: this is a best response and not a contraction, so a cycle
        has to be visible rather than hidden behind a maxIter. 'iterates' carries the whole sequence. """
        θ, hist = float(θ0), [float(θ0)]
        with self.BG.cacheParams():
            for _ in range(maxIter):
                rec = self.solve(t0, self.siRatioAt(t0, θ), θCand = θCand, s_ = s_)
                hist.append(rec['θ'])
                done = abs(rec['θ'] - θ) < tol
                θ = rec['θ']
                if done:
                    break
            return rec | {'iterates': np.array(hist), 'converged': done,
                          'siRatio': self.siRatioAt(t0, θ)}


class PermanentCRRA(LeadedBase):
    r""" The permanent choice under CRRA.

    Simpler than LeadedCRRA, not harder: a permanent theta means a CONSTANT design path, so each candidate
    is one ordinary solvePEE_CRRA call and there is no path to iterate. PermanentLOG's concentration
    argument carries over unchanged -- tau is whatever the CRRA politico-economic equilibrium delivers at
    that design -- so the whole problem is a one-dimensional grid search over theta.

    s_{t0-1,i}/s_{t0-1} is pinned for the same reason as in the LOG case, and here it must come from a
    SOLVED equilibrium: under CRRA it is not a closed form. That is why the candidate solves are cached
    (_grid) rather than folded into the maximisation -- the anticipated-vote fixed point re-weights the
    same equilibria against a different pinned ratio, so it costs iterations of arithmetic rather than
    iterations of solvePEE_CRRA. """

    def __init__(self, m, nθCand = 21):
        self.m = m
        self.B, self.BG, self.BT = m.B, m.BG, m.BT
        self.db = m.db
        self.ni, self.T = m.ni, m.T
        self.nθCand = nθCand
        self.θCand = np.linspace(0., 1., nθCand)

    def _grid(self, t0, ths, s0 = None, solveKwargs = None, verbose = True):
        """ One ordinary CRRA equilibrium per candidate permanent design, reduced to what W needs.

        The candidate is passed to solvePEE_CRRA and also written into db['θ'], so that anything reading
        the design in force from db (solvePEE_*'s θ = None default) sees the candidate; no updateAuxPars
        runs in between, which would re-derive θ from the data (#9). db['θ'] is restored in a finally
        block, so a failed candidate cannot leave db on the wrong design.

        'si' is each candidate's OWN s_{t0-1,i}/s_{t0-1}. W never uses it -- solveFixedPoint does, to
        find the ratio the anticipated vote implies. 'mult' is each candidate's tax counts (solvePEE_CRRA's
        'multiplicity'; None where it failed). A candidate whose solve fails (SOLVE_FAILURES) keeps ok=False
        and NaN rows, its error in 'errors' (None where it solved) and the count in 'nFailed'; any other
        exception propagates. """
        self._requireZeroMass()
        t = self.db['t'][self.db['t'].get_loc(t0)]
        ε = self.db['eps'].values.astype(float)
        n = len(ths)
        c = {'θ': ths, 'ok': np.zeros(n, dtype = bool), 'τ': np.full(n, np.nan),
             'h': np.full(n, np.nan), 's_': np.full(n, np.nan),
             'c1': np.full((n, self.ni), np.nan), 'B': np.full((n, self.ni), np.nan),
             'si': np.full((n, self.ni), np.nan), 'mult': [None]*n, 'errors': [None]*n}
        thSave = self.db['θ'].values.astype(float).copy()
        try:
            for k, cand in enumerate(ths):
                self.db.update(self.m.adjPar('θ', float(cand)))
                try:
                    out = self.m.solvePEE_CRRA(θ = np.full(self.T, float(cand)), ε = ε, s0 = s0,
                                               **(solveKwargs or {}))
                    rep = out['report']
                    tPrev = self.db['t'][self.db['t'].get_loc(t0)-1]
                    c['c1'][k] = rep['tildec1i'].xs(t).values.astype(float)
                    c['B'][k] = rep['B'].xs(t).values.astype(float)
                    c['si'][k] = rep['si_s'].xs(tPrev).values.astype(float)
                    c['h'][k], c['s_'][k] = float(rep['h'].xs(t)), float(rep['s_'].xs(t))
                    c['τ'][k] = float(out['τ'].xs(t0))
                    c['mult'][k] = out.get('multiplicity')
                    c['ok'][k] = True
                except SOLVE_FAILURES as e:
                    c['errors'][k] = f'{type(e).__name__}: {e}'
                    if verbose:
                        print('      theta={:.3f}: failed ({}: {})'.format(cand, type(e).__name__, e))
        finally:
            self.db.update(self.m.adjPar('θ', thSave))
        c['nFailed'] = int(n - c['ok'].sum())
        if not c['ok'].any():
            raise RuntimeError('PermanentCRRA: every candidate failed to solve.')
        return c

    def W(self, cache, siRatio_, t0):
        """ W_{t0} at every cached candidate, with c_{2,t0}^i rebuilt at the PINNED predetermined ratio
        rather than read off each candidate's own report, whose si_s moved with the candidate.

        Evaluated with db back on the incumbent design, which is safe because c2i takes theta and tau
        explicitly and nothing in its chain reads db['θ'] -- wedgeA/wedgeB are functions of their (θ, τ)
        arguments (and, under 'size', of the income distribution through Vtilde, not of the design). """
        t = self.db['t'][self.db['t'].get_loc(t0)]
        q = 1 - 1/float(self.BG.get('ρ', t))
        old, young = self.weights(t)
        out = np.full(len(cache['θ']), -np.inf)
        with self.BG.cacheParams():
            for k in np.flatnonzero(cache['ok']):
                c2 = np.asarray(self.BG.c2i(cache['h'][k], cache['s_'][k], cache['τ'][k],
                                            float(cache['θ'][k]), siRatio_, t), dtype = float)
                v1 = (1+cache['B'][k])*cache['c1'][k]**q/q
                out[k] = float((old*(c2**q/q)).sum() + (young*v1).sum())
        return out

    def solve(self, t0, siRatio_, s0 = None, θCand = None, solveKwargs = None, verbose = True,
              cache = None):
        """ Grid-search the permanent design under CRRA at a given pinned ratio. One full solvePEE_CRRA
        per candidate, unless a cache from a previous _grid call is supplied. 'multiplicity': the tax
        counts over every candidate's solve (policy.mergeMultiplicity of the cache's 'mult'); 'nFailed':
        the candidates whose solve failed (_grid). """
        ths = self.θCand if θCand is None else np.asarray(θCand, dtype = float)
        if cache is None:
            cache = self._grid(t0, ths, s0 = s0, solveKwargs = solveKwargs, verbose = verbose)
        Ws = self.W(cache, siRatio_, t0)
        thStar, atBound = self._argmax(cache['θ'], Ws)
        ok = cache['ok']
        return {'θ': thStar, 'atBound': atBound, 'W': Ws, 'θCand': cache['θ'], 'cache': cache,
                'τAtChoice': float(np.interp(thStar, cache['θ'][ok], cache['τ'][ok])),
                'multiplicity': mergeMultiplicity(cache.get('mult', [])),
                'nFailed': int(cache.get('nFailed', int((~ok).sum())))}

    def solveFixedPoint(self, t0, θ0, s0 = None, θCand = None, solveKwargs = None, verbose = True,
                        tol = 1e-6, maxIter = 50):
        """ PermanentLOG.solveFixedPoint's timing under CRRA: theta* = argmax W(theta ; siRatio(theta*)).

        The candidate equilibria are solved ONCE and the iteration re-weights them, so the fixed point
        costs no extra solvePEE_CRRA calls. siRatio(theta) is interpolated linearly across candidates --
        the same piecewise-linear rule the rest of the module uses, and for the same reason: an object
        feeding an argmax must not carry interpolation artefacts.

        tol is looser than the LOG version's (1e-6 against 1e-9) because each candidate here is a numerical
        equilibrium solve, not a closed form. """
        ths = self.θCand if θCand is None else np.asarray(θCand, dtype = float)
        cache = self._grid(t0, ths, s0 = s0, solveKwargs = solveKwargs, verbose = verbose)
        ok = cache['ok']
        siOf = lambda x: np.array([np.interp(x, cache['θ'][ok], cache['si'][ok, i])
                                   for i in range(self.ni)])
        θ, hist = float(θ0), [float(θ0)]
        for _ in range(maxIter):
            rec = self.solve(t0, siOf(θ), θCand = ths, cache = cache)
            hist.append(rec['θ'])
            done = abs(rec['θ'] - θ) < tol
            θ = rec['θ']
            if done:
                break
        return rec | {'iterates': np.array(hist), 'converged': done, 'siRatio': siOf(θ)}

import numpy as np, pandas as pd
from contextlib import contextmanager
def noneInit(x, fallBackValue):
    return fallBackValue if x is None else x


class Base:
    """ Building blocks

    Convention: parameters (α, ξ, ν, γ, η, X, β, p, κ, Γh, ...) are read from the database via
    self()/self.get(). Objects that may change during a solve or be chosen by policy -- taxes τ,
    pension characteristics θ/ε, aggregate states s, h and their lags/leads, discount factors B --
    are always explicit function arguments, never read from db, so this class stays agnostic to
    how those objects end up being solved for (fixed path, one-off choice, no-commitment policy, ...). """
    def __init__(self, m, t = None):
        self.m = m # associated Model class instance
        self.db = m.db # main database
        self.tFirst = self.db['t'][0] # model's first active/endogenous period (docs' t=1) -- NOT the
        # pre-determined state before it (docs' t=0, code's `s0`) and NOT db['t0'] (calibration baseline
        # year index). See README, "Timing and data conventions".
        self.t = t # current year
        self._cache = None # parameter cache; None = disabled. See cacheParams().

    #######################################################################
    ##########   Parameter caching (cacheParams)                  #########
    #######################################################################
    # Opt-in and block-scoped, not always-on -- a cache surviving a db rewrite (model.py's calibration)
    # would return stale parameters silently.
    @contextmanager
    def cacheParams(self):
        """ Memoise db parameter reads for the duration of the block (see the §header comment above).

        Use around any loop that holds the year fixed while evaluating repeatedly -- a period of a policy
        grid search, say. Nests safely: an inner block reuses the outer block's cache rather than starting
        a second one, and only the outermost exit clears it. """
        if self._cache is not None:
            yield self # already inside a block -- reuse it, and let that block own the teardown
            return
        self._cache = {}
        try:
            yield self
        finally:
            self._cache = None

    def _year(self, t = None):
        """ The year a lookup actually resolves to -- the cache key's time component. Clamped at tFirst,
        matching __call__, so `t=None`, an explicit t, and any t before the horizon all share one entry. """
        return max(noneInit(t, self.t), self.tFirst)

    def _memo(self, key, fn):
        """ Return fn(), memoised under `key` while a cacheParams() block is active (otherwise just
        fn()). Keys must include the resolved year and any `lag` argument -- anything that changes the
        result. """
        c = self._cache
        if c is None:
            return fn()
        if key not in c:
            c[key] = fn()
        return c[key]

    def __call__(self, k, t = None):
        """ Return symbol 'k' from database. If t is not provided, rely on self.t attribute as backup. """
        return self._memo(('call', k, self._year(t)),
                          lambda: self.db[f'{k}'].xs(self._year(t)))

    def get(self, k, t = None):
        """ return numpy version of __call__ method """
        def _v():
            s = self(k, t = t)
            return s.values if isinstance(s, (pd.Series, pd.DataFrame)) else s
        return self._memo(('get', k, self._year(t)), _v)

    #######################################################################
    ##########                    0. Aux methods                ###########
    #######################################################################
    # Political weights:
    def ω2i(self, t = None):
        """ Political weight on types 2i, p_{t-1,i}·ω_t·μ_{t-1,i}: the voting weight is the old generation's
        own, set when it was young (μ attaches to a generation). """
        return (self('pi[t-1]', t).mul(self('ω',t), axis = 0).mul(self('μi[t-1]',t), axis = 0)).values
    def ω20(self, t = None):
        """ Political weight on types 20. Uses μ_{t-1,0} -- see ω2i's docstring. """
        return self.get('p0[t-1]', t) * self.get('ω',t) * self.get('μ0[t-1]',t)
    def ω1i(self, t = None):
        """ Political weight on types 1i """
        return self.get('μi',t)
    def ω10(self, t = None):
        """ Political weight on types 10 """
        return self.get('μ0',t)

    # Power functions:
    def power_s(self, t = None):
        """ α(1+ξ)/(1+αξ) """
        return self.get('α',t)*(1+self.get('ξ',t))/(1+self.get('α',t)*self.get('ξ',t))
    def power_h(self, t = None):
        """ αξ/(1+αξ) """
        return self.get('α',t)*self.get('ξ',t)/(1+self.get('α',t)*self.get('ξ',t))

    # Auxiliary parameters over household types:
    def auxProd(self, t = None, lag = ''):
        """ ηi^(1+ξ)/Xi^ξ. Pass lag = '[t-1]' or '[t+1]' to use a shifted vintage of η,X,ξ
        (needed e.g. in the pension benefit formula, which references last period's productivity distribution).
        Memoised inside cacheParams(). """
        return self._memo(('auxProd', lag, self._year(t)), lambda:
            (self(f'ηi{lag}', t).pow(1+self(f'ξ{lag}',t), axis = 0)/self(f'Xi{lag}',t).pow(self(f'ξ{lag}',t), axis = 0)).values)

    def Γh(self, t = None, lag = ''):
        """ ∑i γi * ηi^(1+ξ)/Xi^ξ. Memoised inside cacheParams(). Deliberately still *computed* rather
        than read from db['Γh'] (holds the same value via model.py's paramsFromFuncs): that's only
        refreshed by updateAuxPars, so it would go stale mid-calibration when ηi/Xi are being solved for. """
        return self._memo(('Γh', lag, self._year(t)), lambda:
            (self(f'γi{lag}', t) * self.auxProd(t, lag = lag)).sum())

    @staticmethod
    def _bcast(x):
        """ Reshape a (t,)-shaped array/scalar to (t,1) so it broadcasts against (t,i)-shaped arrays.
        No-op for true scalars / single-year values (Base, BaseGrid). """
        x = np.asarray(x)
        return x[:, None] if x.ndim else x

    # Two distinct ratios that differ by a factor η_{t,i} -- keep them apart:
    #   hRatio  = h_{t,i}/h_t          = (η_{t,i}/X_{t,i})^ξ / Γ_{h,t}     (eq EE:hi)
    #   hηRatio = h_{t,i}η_{t,i}/h_t   = (η_{t,i}^{1+ξ}/X_{t,i}^ξ) / Γ_{h,t}
    # hηRatio is what the doc's si_s/c2i/dv2i formulas carry (they all write η^{1+ξ}/X^ξ over Γh);
    # hRatio is what individual labour supply and the benefit formula need. Both are parameters only
    # (no solution state), so both memoise inside cacheParams() like their ingredients.
    def hRatio(self, t = None, lag = ''):
        """ Eq (EE:hi): h_{t,i}/h_t = (η_{t,i}/X_{t,i})^ξ/Γ_{h,t}. Satisfies ∑_i γ_{t,i}η_{t,i}·hRatio_i=1. """
        return self._memo(('hRatio', lag, self._year(t)), lambda:
            (self(f'ηi{lag}', t)/self(f'Xi{lag}', t)).pow(self(f'ξ{lag}', t), axis = 0).values
            / self._bcast(self.Γh(t, lag = lag)))

    def hηRatio(self, t = None, lag = ''):
        """ h_{t,i}η_{t,i}/h_t = (η_{t,i}^{1+ξ}/X_{t,i}^ξ)/Γ_{h,t}. Satisfies ∑_i γ_{t,i}·hηRatio_i = 1. """
        return self._memo(('hηRatio', lag, self._year(t)), lambda:
            self.auxProd(t, lag = lag) / self._bcast(self.Γh(t, lag = lag)))

    #######################################################################
    ##########   1. Factor prices (eq:factorPrices, eq:factorPrices0)  #####
    #######################################################################
    def R(self, s_, h, t = None):
        """ Eq (factorPrices): R_t = α (s_{t-1}/(ν_t h_t))^(α-1). s_ = s_{t-1}, h = h_t (explicit -- solution objects). """
        α = self.get('α', t)
        return α * (s_/(self.get('ν', t)*h))**(α-1)

    def w(self, s_, h, t = None):
        """ Eq (factorPrices): w_t = (1-α)(s_{t-1}/(ν_t h_t))^α. TFP A_t normalized to 1. s_ = s_{t-1}, h = h_t. """
        α = self.get('α', t)
        return (1-α) * (s_/(self.get('ν', t)*h))**α

    def w0(self, s_, t = None):
        """ Eq (factorPrices0): informal wage w_t^0 = ((1-α)/Γh^α)^(1/(1+αξ)) (s_{t-1}/ν_t)^(α/(1+αξ)).
        s_ = s_{t-1}. No χ^R scaling here -- that applies to the informal *return* only (see R0). """
        α, ξ = self.get('α', t), self.get('ξ', t)
        return ((1-α)/self.Γh(t)**α)**(1/(1+α*ξ)) * (s_/self.get('ν', t))**(α/(1+α*ξ))

    def R0(self, s_, h, t = None):
        """ Eq (factorPrices0): informal return R_t^0 = R_t·χ_t^R. s_ = s_{t-1}, h = h_t.

        χ^R carries a *period* index, not a generation index: R_t^0 is the return earned between t-1 and t,
        so the discounting of b_{t+1}^0 inside s_{t,0} uses χ^R_{t+1} (see auxInf1) while c_{2,t}^0 uses
        χ^R_t (auxInf1_). """
        return self.R(s_, h, t) * self.get('χR', t)

    #######################################################################
    ##########   2. Individual labour supply (eq:EE:hi)           ##########
    #######################################################################
    def hi(self, h, t = None):
        """ h_{t,i} = h_t * (h_{t,i}/h_t). h = aggregate h_t (explicit -- solution object). """
        return self._bcast(h) * self.hRatio(t)

    #######################################################################
    ##########   3. Pension system / government budget (eq:governmentBudget)  ###
    ##########      τ, θ, ε are always explicit -- to be endogenized later     ###
    #######################################################################
    def κ(self, ε1, t = None):
        """ Eq (governmentBudget:bbar): κ_t(ε_{t+1}) = p_t + ε_{t+1}γ_{t,0}p_{t,0}, so that generation t's
        retirees per formal young cost Σ_j γ_{t,j}p_{t,j}b_{t+1}^j = κ_t h_t bbar_{t+1}. n_t and ν_t count
        FORMAL young (eq:aggregateStates, the factor prices' s_{t-1}/(ν_t h_t), FOC's young weight ν_t), so
        the budget n_t w_t τ_t h_t = n_{t-1} Σ_j γ_{t-1,j}p_{t-1,j}b_t^j carries no (1+γ_0) factor, also
        when γ_0 varies over t. ε1 = ε_{t+1} (explicit -- mirrors the b0/bi convention, since ε may be
        endogenized later; γ0/p/p0 are read from db as primitives, not solution objects). """
        p, p0, γ0 = self.get('p', t), self.get('p0', t), self.get('γ0', t)
        return p + ε1*γ0*p0

    def bbar(self, τ, w, h, h_, t = None):
        """ Eq (governmentBudget:bbar): bbar_t = ν_t w_t h_t τ_t / (h_{t-1} κ_{t-1}).
        τ, w, h = period-t tax rate/wage/aggregate hours; h_ = h_{t-1} (aggregate, lagged). All explicit. """
        return self.get('ν', t)*w*h*τ / (h_*self.get('κ[t-1]', t))

    def bi(self, θ, bbar, h_, t = None):
        """ Eq (governmentBudget): b_t^i = [θ_t h_{t-1,i} η_{t-1,i} + (1-θ_t) h_{t-1}] bbar_t.
        θ, bbar = period-t contributive-incentive parameter / benefit level (explicit -- θ will be endogenized).
        h_ = h_{t-1} (aggregate, lagged). """
        hiη_ = self._bcast(h_) * self.hηRatio(t, lag = '[t-1]')
        bracket = self._bcast(θ)*hiη_ + self._bcast((1-θ)*h_)
        return bracket * self._bcast(bbar)

    def b0(self, ε, bbar, h_):
        """ Eq (governmentBudget): b_t^0 = ε_t h_{t-1} bbar_t.
        ε, bbar = period-t universal-pension parameter / benefit level (explicit -- ε will be endogenized).
        h_ = h_{t-1} (aggregate, lagged). """
        return ε * h_ * bbar

    # auxInf1/auxInf1_: the informal household's discounted pension coefficient, i.e. what multiplies
    # τ_{t+1}·s_t once b_{t+1}^0/(R_{t+1}^0/p_{t,0}) is written out in equilibrium. Used by s0_s/c10/
    # tildec10 (forward, χ^R_{t+1}) and c20 (lagged, χ^R_t) -- see R0's docstring on the χ^R timing.
    def auxInf1(self, ε1, t = None):
        """ (1-α)/(χ^R_{t+1} α) · p_{t,0} ε_{t+1}/κ_t. ε1 = ε_{t+1} (explicit). """
        α = self.get('α', t)
        return (1-α)/(self.get('χR[t+1]', t)*α) * self.get('p0', t)*ε1/self.get('κ', t)

    def auxInf1_(self, ε, t = None):
        """ (1-α)/(χ^R_t α) · p_{t-1,0} ε_t/κ_{t-1} -- auxInf1 shifted back one period. ε = ε_t (explicit). """
        α = self.get('α', t)
        return (1-α)/(self.get('χR', t)*α) * self.get('p0[t-1]', t)*ε/self.get('κ[t-1]', t)

    #######################################################################
    ##########   4. Numerical auxiliary functions (eq:auxiliary)   #########
    ##########      -- the core building blocks used to solve the   #########
    ##########      economic equilibrium given a policy path        #########
    #######################################################################
    def Γs(self, B, τ1, θ1, t = None):
        """ Eq (auxiliary:Gammas): Γ_{s,t}(B_{t+1}, τ_{t+1}, θ_{t+1}).
        B = B_{t+1}^i (per type, explicit); τ1 = τ_{t+1}; θ1 = θ_{t+1} (explicit -- to be endogenized). """
        α, ξ, p, κ = self.get('α', t), self.get('ξ', t), self.get('p', t), self.get('κ', t)
        γi, auxProd = self.get('γi', t), self.auxProd(t)
        Bratio = B/(1+B)
        num = (γi*auxProd*Bratio).sum(axis = -1)/(1+ξ)
        denom = 1 + (1-α)/α * p*τ1/κ * (θ1 + (1-θ1)*(γi/(1+B)).sum(axis = -1))
        return num/denom

    def Θh(self, τ, τ1, θ1, Γs, t = None):
        """ Eq (auxiliary:Thetah): Θ_{h,t}(τ_t, τ_{t+1}, θ_{t+1}, Γ_{s,t}). """
        α, ξ, p, κ, Γh = self.get('α', t), self.get('ξ', t), self.get('p', t), self.get('κ', t), self.Γh(t)
        return Γh**((1+ξ)/(1+α*ξ)) * ((1-α)*(1-τ)/(Γh - (1-α)/α*p*θ1*τ1/κ*Γs))**(ξ/(1+α*ξ))

    def ΘhTerminal(self, τ, t = None):
        """ Eq (auxiliary:ThetahT): terminal-period Θ_{h,T}(τ_T). """
        α, ξ, Γh = self.get('α', t), self.get('ξ', t), self.Γh(t)
        return Γh**(1/(1+α*ξ)) * ((1-α)*(1-τ))**(ξ/(1+α*ξ))

    def Θs(self, Θh, Γs, t = None):
        """ Eq (auxiliary:Thetas): Θ_{s,t}(Θ_{h,t}, Γ_{s,t}). """
        ξ, Γh = self.get('ξ', t), self.Γh(t)
        return (Θh/Γh)**((1+ξ)/ξ) * Γs

    def s(self, Θs, s_, t = None):
        """ Eq (auxiliary:s): s_t(Θ_{s,t}, s_{t-1}). s_ = s_{t-1} (explicit). """
        return Θs * (s_/self.get('ν', t))**self.power_s(t)

    def h(self, Θh, s_, t = None):
        """ Eq (auxiliary:h): h_t(Θ_{h,t}, s_{t-1}). s_ = s_{t-1} (explicit). """
        return Θh * (s_/self.get('ν', t))**self.power_h(t)

    def ΘhFromH(self, h, s_, t = None):
        """ Eq (auxiliary:h), inverted: Θ_{h,t}(h_t, s_{t-1}). Recovers Θh from an already-solved
        equilibrium regardless of which formula (general/terminal) produced h -- used by §10. """
        return h / (s_/self.get('ν', t))**self.power_h(t)

    def sFromH(self, h, Γs, t = None):
        """ Eq (auxiliary:sFromH): alternative s_t(h_t, Γ_{s,t}), given h_t directly instead of s_{t-1}. """
        ξ = self.get('ξ', t)
        return (h/self.Γh(t))**((1+ξ)/ξ) * Γs

    def hFromS(self, s, Γs, t = None):
        """ Eq (auxiliary:sFromH), inverted: h_t(s_t, Γ_{s,t}) -- recovers aggregate hours from a *given*
        savings level and Γ_{s,t}, the reverse direction of sFromH. Used by model.py's initialState_solve
        to back out h_{-1} (docs' t=0's aggregate hours) from the actual given s0 and a Γ_{s,-1} computed
        from db['t'][0]'s own (τ,θ) -- see that method's docstring. Round-trips exactly with sFromH
        (verified numerically). """
        ξ = self.get('ξ', t)
        return self.Γh(t) * (s/Γs)**(ξ/(1+ξ))

    def Rlead(self, s, h1, t = None):
        """ Eq (auxiliary:R): R_{t+1}(s_t, h_{t+1}). s = s_t, h1 = h_{t+1} (both explicit). """
        α = self.get('α', t)
        return α * (s/(self.get('ν[t+1]', t)*h1))**(α-1)

    def B(self, s, h1, t = None):
        """ Eq (auxiliary:B): B_{t+1}^i(s_t, h_{t+1}) = (β_{t,i})^ρ [R_{t+1}(s_t,h_{t+1})/p_t]^(ρ-1). """
        ρc = self._bcast(self.get('ρ', t))
        Rp = self._bcast(self.Rlead(s, h1, t)/self.get('p', t))
        return self.get('βi', t)**ρc * Rp**(ρc-1)

    def Rlead0(self, s, h1, t = None):
        """ Eq (factorPrices0) one period ahead: R_{t+1}^0(s_t, h_{t+1}) = R_{t+1}·χ^R_{t+1}. """
        return self.Rlead(s, h1, t) * self.get('χR[t+1]', t)

    def B0(self, s, h1, t = None):
        """ Eq (auxiliary:B), informal counterpart: B_{t+1}^0 = (β_{t,0})^ρ [R_{t+1}^0/p_{t,0}]^(ρ-1).
        Scalar per period (one informal type), so no _bcast. Collapses to the primitive β_{t,0} at ρ=1. """
        ρ = self.get('ρ', t)
        return self.get('β0', t)**ρ * (self.Rlead0(s, h1, t)/self.get('p0', t))**(ρ-1)

    def lnRleadΘ(self, Θs, Θh1, t = None):
        """ Eq (auxiliary:R) in coefficient form: ln(R_{t+1}) up to an additive constant along τ_t,
            ln(R_{t+1}) = (α-1)[(1-power_h)ln(Θ_{s,t}) - ln(Θ_{h,t+1})] + const.
        Substituting s_t = Θ_{s,t}(s_{t-1}/ν_t)^power_s and h_{t+1} = Θ_{h,t+1}(s_t/ν_{t+1})^power_h leaves
        s_{t-1} only inside the dropped constant, which is all a LOG first order condition needs
        (docs eq:logsep). This is why the political recursion carries Θ_{h,t+1} rather than h_{t+1}: the
        level also depends on s_t and so is not a function of the state (docs §PEELOG).

        α/power_h are read at t, following Rlead's own convention (both are t+1 objects exactly;
        immaterial unless α/ξ vary over t). """
        α = self.get('α', t)
        return (α-1)*((1-self.power_h(t))*np.log(Θs) - np.log(Θh1))

    #######################################################################
    ##########   5. Individual savings ratio (eq:EE:si_s)          #########
    #######################################################################
    def si_s(self, B, τ1, θ1, Γs, t = None):
        """ Eq (EE:si_s): s_{t,i}/s_t(B_{t+1}, τ_{t+1}, θ_{t+1}, Γ_{s,t}).
        B = B_{t+1}^i (per type); τ1 = τ_{t+1}; θ1 = θ_{t+1}; Γs = Γ_{s,t} (all explicit -- solution objects).
        Rewritten in terms of Γ_{s,t}'s own numerator/denominator (Γs = num/denom, see eq:auxiliary:Gammas)
        so we avoid recomputing the aggregate sum ∑_j γ_j(η_j^{1+ξ}/X_j^ξ)B_{t+1}^j/(1+B_{t+1}^j) a second time. """
        α, p, κ = self.get('α', t), self.get('p', t), self.get('κ', t)
        auxProd, Bratio = self.auxProd(t), B/(1+B)
        term1 = Bratio*auxProd / self._bcast((1+self.get('ξ', t))*Γs)
        term2 = -(1/(1+B)) * self._bcast((1-α)/α * p*(1-θ1)/κ*τ1)
        term3 = -self.hηRatio(t) * self._bcast((1-α)/α * p*θ1/κ*τ1)
        return term1 + term2 + term3

    def s0_s(self, B0, Θs, τ1, ε1, t = None):
        """ Eq (EE:s0_s)/(auxiliary:s0_s): ι_t = s_{t,0}/s_t (B_{t+1}^0, Θ_{s,t}, τ_{t+1}, ε_{t+1}).

        Denominator is *aggregate formal* savings s_t -- informal savings stay out of the formal capital
        stock, so this is a ratio across two different pools, not a share. Unlike si_s (a function of
        τ_{t+1} alone), ι_t depends on τ_t too, through Θ_{s,t}: that is what makes it an endogenous state
        of the political problem. """
        α, ξ = self.get('α', t), self.get('ξ', t)
        smooth = B0/(1+B0) * self.auxProd0(t)/(1+ξ) * ((1-α)/self.Γh(t)**α)**((1+ξ)/(1+α*ξ)) / Θs
        pension = -self.auxInf1(ε1, t)*τ1 / (1+B0)
        return smooth + pension

    #######################################################################
    ##########   6. Closed-form consumption -- formal (eq:EE:ci)    ########
    #######################################################################
    def c1i(self, h, s, B, τ1, θ1, t = None):
        """ Eq (EE:ci): c_{1,t}^i(h_t, s_t, B_{t+1}, τ_{t+1}, θ_{t+1}). """
        ξ, α, p, κ = self.get('ξ', t), self.get('α', t), self.get('p', t), self.get('κ', t)
        auxProd, Bratio = self.auxProd(t), B/(1+B)
        smooth = auxProd * self._bcast((h/self.Γh(t))**((1+ξ)/ξ)) * (1 - Bratio/self._bcast(1+ξ))
        pension = self._bcast(s)/(1+B) * self._bcast((1-α)/α*p*τ1*(1-θ1)/κ)
        return smooth + pension

    def tildec1i(self, h, B, τ1, θ1, Γs, t = None):
        """ Eq (EE:ci): tilde-c_{1,t}^i(h_t, B_{t+1}, τ_{t+1}, θ_{t+1}, Γ_{s,t}) -- s_t eliminated using Γ_{s,t}. """
        ξ, α, p, κ = self.get('ξ', t), self.get('α', t), self.get('p', t), self.get('κ', t)
        auxProd = self.auxProd(t)
        hΓh = self._bcast((h/self.Γh(t))**((1+ξ)/ξ))
        bracket = auxProd/self._bcast(1+ξ) + self._bcast(Γs*(1-α)/α*p*τ1*(1-θ1)/κ)
        return hΓh/(1+B) * bracket

    # ĉ_{1,t}^i ≡ (1+B_{t+1}^i)^{1/(1-1/ρ)}·tilde-c_{1,t}^i (docs eq:hatc1) folds B_{t+1}^i's own τ_t
    # -dependence (through s_t) into the consumption level, so υ_{1,t}^i = ĉ^{1-1/ρ}/(1-1/ρ) and one
    # numerical dln(ĉ)/dτ_t picks up both channels.
    #
    # ĉ itself is never computed -- only (ĉ)^{1-1/ρ} and ln(ĉ) are needed, and both avoid the
    # (1/(1-1/ρ))-power that the literal definition requires:
    #       (ĉ)^{1-1/ρ} = (1+B)·(tilde-c)^{1-1/ρ},        ln(ĉ) = ln(1+B)/(1-1/ρ) + ln(tilde-c).
    # DO NOT replace these with the literal form: it overflows float64 as ρ->1 (exponent 1001 at ρ=1.001,
    # overflowing once B>1.03), silently yielding inf. Neither form above overflows for any ρ≠1.
    #
    # ρ=1 itself is genuinely undefined (indirect utility is (1+B)ln(tilde-c), not a power of any level)
    # -- not a gap, since ρ=1 is exactly what policy.py's LOG class solves in closed form.
    def hatc1iPow(self, h, B, τ1, θ1, Γs, t = None):
        """ (ĉ_{1,t}^i)^{1-1/ρ} = (1+B_{t+1}^i)·(tilde-c_{1,t}^i)^{1-1/ρ} -- the level factor of
        dυ_{1,t}^i/dτ_t. See the §header comment above for why this form rather than the literal one. """
        p = 1 - 1/self.get('ρ', t)
        return (1+B) * self.tildec1i(h, B, τ1, θ1, Γs, t)**p

    def lnhatc1i(self, h, B, τ1, θ1, Γs, t = None):
        """ ln(ĉ_{1,t}^i) = ln(1+B_{t+1}^i)/(1-1/ρ) + ln(tilde-c_{1,t}^i) -- differentiate this along τ_t
        to get dln(ĉ_{1,t}^i)/dτ_t. See the §header comment above. """
        p = 1 - 1/self.get('ρ', t)
        return np.log1p(B)/p + np.log(self.tildec1i(h, B, τ1, θ1, Γs, t))

    def c2i(self, h, s_, τ, θ, siRatio_, t = None):
        """ Eq (EE:ci): c_{2,t}^i(h_t, s_{t-1}, τ_t, θ_t, s_{t-1,i}/s_{t-1}).
        siRatio_ = s_{t-1,i}/s_{t-1}, itself a solution object (si_s() evaluated at t-1) -- predetermined at t,
        explicit rather than read from db (mirrors θ/τ convention). """
        α, ν, p_, κ_ = self.get('α', t), self.get('ν', t), self.get('p[t-1]', t), self.get('κ[t-1]', t)
        A = (1-α)/α * p_*τ/κ_
        inner = siRatio_ + self._bcast(A) * (1 + self._bcast(θ)*(self.hηRatio(t, lag = '[t-1]') - 1))
        outer = α * (ν/p_) * h**(1-α) * (s_/ν)**α
        return self._bcast(outer) * inner

    def tildec2i(self, h, s_, τ, θ, siRatio_, t = None):
        """ Eq (EE:ci): tilde-c_{2,t}^i = c_{2,t}^i -- formal retirees don't supply labour when old. """
        return self.c2i(h, s_, τ, θ, siRatio_, t = t)

    #######################################################################
    ##########   7. Closed-form consumption -- informal (eq:EE:c0)  ########
    #######################################################################
    def auxProd0(self, t = None, lag = ''):
        """ η0^(1+ξ)/X0^ξ -- informal analogue of auxProd (scalar type-0 productivity/disutility; ξ is
        always current-period). Memoised inside cacheParams(). """
        return self._memo(('auxProd0', lag, self._year(t)), lambda:
            self.get(f'η0{lag}', t)**(1+self.get('ξ', t)) / self.get(f'X0{lag}', t)**self.get('ξ', t))

    def h0(self, s_, t = None):
        """ Eq (informalOpt): h_{t,0}(s_{t-1}) = (η_{t,0} w_t^0 / X_{t,0})^ξ. """
        return (self.get('η0', t)*self.w0(s_, t)/self.get('X0', t))**self.get('ξ', t)

    def c10(self, s_, s, B0, τ1, ε1, t = None):
        """ Eq (EE:c0): c_{1,t}^0(s_{t-1}, s_t, B_{t+1}^0, τ_{t+1}, ε_{t+1}). Mirrors c1i's argument order;
        s_ enters through the informal wage w_t^0, s_t through the discounted universal pension. """
        ξ = self.get('ξ', t)
        smooth = self.auxProd0(t) * self.w0(s_, t)**(1+ξ) * (1 - (B0/(1+B0))/(1+ξ))
        pension = s/(1+B0) * self.auxInf1(ε1, t)*τ1
        return smooth + pension

    def tildec10(self, s_, s, B0, τ1, ε1, t = None):
        """ Eq (EE:c0): tilde-c_{1,t}^0 -- c10 with the GHH labour-disutility term netted out. """
        ξ = self.get('ξ', t)
        smooth = self.auxProd0(t) * self.w0(s_, t)**(1+ξ) / ((1+ξ)*(1+B0))
        pension = s/(1+B0) * self.auxInf1(ε1, t)*τ1
        return smooth + pension

    def tildec10Θ(self, Θs, B0, τ1, ε1, t = None):
        """ Eq (EE:c0)/(EE:sigma_ci): tilde-c_{1,t}^0's coefficient function, i.e.
        tilde-c_{1,t}^0 = tildec10Θ · (s_{t-1}/ν_t)^power_s. Both of tildec10's terms carry that identical
        power of s_{t-1}: w_t^0(s_{t-1})^{1+ξ} = ((1-α)/Γ_h^α)^{(1+ξ)/(1+αξ)}(s_{t-1}/ν_t)^power_s (the
        same coefficient s0_s carries), and s_t = Θ_{s,t}(s_{t-1}/ν_t)^power_s. Needed because the LOG
        political FOC evaluates ln(tilde-c_{1,t}^0) up to a τ_t-constant (docs eq:logsep, eq:v1LOG); the
        formal twin needs no such method, since tildec1i(Θ_{h,t},...) already *is* its own coefficient
        function (h enters it only through (h/Γ_h)^{(1+ξ)/ξ}). """
        α, ξ = self.get('α', t), self.get('ξ', t)
        smooth = self.auxProd0(t) * ((1-α)/self.Γh(t)**α)**((1+ξ)/(1+α*ξ)) / ((1+ξ)*(1+B0))
        pension = Θs/(1+B0) * self.auxInf1(ε1, t)*τ1
        return smooth + pension

    def c20(self, h, s_, ε, τ, ι_, t = None):
        """ Eq (EE:c0): c_{2,t}^0(h_t, s_{t-1}, ε_t, τ_t, s_{t-1,0}/s_{t-1}). ι_ = ι_{t-1} (explicit --
        predetermined at t, the informal counterpart of c2i's siRatio_, and the political problem's
        endogenous state). """
        α, ν, p0_ = self.get('α', t), self.get('ν', t), self.get('p0[t-1]', t)
        outer = self.get('χR', t) * α * (ν/p0_) * h**(1-α) * (s_/ν)**α
        return outer * (ι_ + self.auxInf1_(ε, t)*τ)

    def tildec20(self, h, s_, ε, τ, ι_, t = None):
        """ Eq (EE:c0): tilde-c_{2,t}^0 = c_{2,t}^0 -- informal retirees supply no labour in this model
        variant (unlike the analytical one, where their old-age endowment was modelled as labour supply). """
        return self.c20(h, s_, ε, τ, ι_, t = t)

    # ĉ_{1,t}^0 ≡ (1+B_{t+1}^0)^{1/(1-1/ρ)}·tilde-c_{1,t}^0 -- informal twin of hatc1iPow/lnhatc1i (see
    # the §6 header comment for why the level itself is never formed: it overflows float64 as ρ->1).
    def hatc10Pow(self, s_, s, B0, τ1, ε1, t = None):
        """ (ĉ_{1,t}^0)^{1-1/ρ} = (1+B_{t+1}^0)·(tilde-c_{1,t}^0)^{1-1/ρ}. """
        p = 1 - 1/self.get('ρ', t)
        return (1+B0) * self.tildec10(s_, s, B0, τ1, ε1, t)**p

    def lnhatc10(self, s_, s, B0, τ1, ε1, t = None):
        """ ln(ĉ_{1,t}^0) = ln(1+B_{t+1}^0)/(1-1/ρ) + ln(tilde-c_{1,t}^0). """
        p = 1 - 1/self.get('ρ', t)
        return np.log1p(B0)/p + np.log(self.tildec10(s_, s, B0, τ1, ε1, t))

    #######################################################################
    ##########   8. Steady state (eq:steadystate_LOG, eq:steadystate_CRRA)  #
    #######################################################################
    # Fixed point s_t=s_{t-1}=s* under constant (τ,θ) at t=self.tFirst's parameters. See README, "Timing
    # and data conventions".
    def sSteadyState(self, Θs, t = None):
        """ Eq (steadystate_LOG:s), generalized: the steady-state savings level s* solving the fixed point
        s* = Θs*(s*/ν_t)^power_s -- i.e. s_t = s_{t-1} = s* under a constant Θs_t = Θs. Algebraically
        identical to eq (steadystate_LOG:s) once Θs is expanded via Θh/Θs (verified by hand), but written
        directly in terms of Θs so the same expression also serves the CRRA steady state (eq
        steadystate_CRRA), which only differs in *how* Γs (and hence Θh/Θs) was obtained. """
        ν, powerS = self.get('ν', t), self.power_s(t)
        return (Θs/ν**powerS)**(1/(1-powerS))

    def RpSteadyState(self, Γs, τ, θ, t = None):
        """ Eq (steadystate_CRRA:Bi)'s bracket: the steady-state gross return per survivor R/p_t implied by
        a candidate Γs under constant policy (τ, θ). """
        α, p, κ, ν, Γh = self.get('α', t), self.get('p', t), self.get('κ', t), self.get('ν', t), self.Γh(t)
        return (α/p) * (Γh - (1-α)/α * p*θ*τ/κ * Γs) / ((1-α)*(1-τ)) * (ν/Γs)

    def BSteadyState(self, Γs, τ, θ, t = None):
        """ Eq (steadystate_CRRA:Bi): steady-state B^i(Γ_s; τ, θ) -- the generalized discount factor
        consistent with a steady state at Γs, given constant policy τ, θ. Γs, τ, θ explicit (Γs is the
        candidate being root-found in steadyState_CRRA_residual; τ, θ are the constant steady-state policy).
        Used to close the Γs fixed point: a candidate Γs is a steady state iff
        self.Γs(self.BSteadyState(Γs, τ, θ), τ, θ) == Γs (see steadyState_CRRA_residual in model.py). At
        ρ=1 this collapses to the primitive βi regardless of Γs -- matching the LOG case's B^i=β_i, so
        BSteadyState is only actually needed for the CRRA (ρ≠1) solve. """
        ρc = self._bcast(self.get('ρ', t))
        return self.get('βi', t)**ρc * self._bcast(self.RpSteadyState(Γs, τ, θ, t))**(ρc-1)

    def ΓsCap(self, τ, θ, t = None):
        """ The Γs at which Θ_{h,t}'s denominator (eq:auxiliary:Thetah) hits zero,

            Γs_cap = Γh·α·κ / ((1-α)·p·θ·τ),      = inf when θτ = 0.

        Above it Θh and RpSteadyState are negative, so a fractional power returns NaN and a bracketing root
        finder dies rather than reporting a bad bracket. A hard feasibility limit: the steady-state root is
        always strictly below it. Scales with α/(1-α), so a constant upper bound that was safe at α = 0.43
        (0.75 against a cap of ≈1.1) is inside the infeasible region at α = 0.35 for τ large enough --
        see model.steadyState_CRRA_bounds and notes/crossCuttingFindings.md #7. """
        denom = (1-self.get('α', t))*self.get('p', t)*θ*τ
        cap = self.Γh(t)*self.get('α', t)*self.get('κ', t)/np.where(denom == 0, np.nan, denom)
        return np.where(np.isnan(cap), np.inf, cap)[()]

    def B0SteadyState(self, Γs, τ, θ, t = None):
        """ Informal counterpart of BSteadyState: B^0 = (β_0)^ρ (R^0/p_{t,0})^(ρ-1), with R^0 = χ^R·R and
        R = p_t·RpSteadyState. Collapses to the primitive β_0 at ρ=1. Not part of the Γs fixed point
        (informal savings never enter Γs) -- only needed to give model.py's initialState_solve an ι_{-1}. """
        ρ = self.get('ρ', t)
        Rp0 = self.get('χR', t) * self.RpSteadyState(Γs, τ, θ, t) * self.get('p', t)/self.get('p0', t)
        return self.get('β0', t)**ρ * Rp0**(ρ-1)

    #######################################################################
    ######   9. Political first-order condition (eq:focbounded, eq:PEELOG)  #
    #######################################################################
    # FOC is the preference-agnostic combiner; `_LOG` methods below are LOG-specific closed forms. CRRA's
    # t<T equivalents need numerical derivatives instead (docs §PEE) -- see dlnc2i_dτ/dlnc20_dτ below.
    def FOC(self, dv1i, dv10, dv2i, dv20, t = None):
        """ Eq (focbounded, the general FOC underlying eq:PEELOG): z_t, the marginal
        political objective w.r.t. τ_t. Combines the young (dv1i,dv10) and old (dv2i,dv20) marginal
        indirect utilities via the political weights (ω1i/ω10/ω2i/ω20, §0) and population shares γ.
        Preference-agnostic: this combination is identical regardless of whether dv1i/dv10/dv2i/dv20 came
        from the LOG closed forms below or a future CRRA version -- only how the dv's themselves are
        computed differs by preference case. The old-generation term uses γ_{t-1,j} (their population share
        when they were young); the young-generation term uses current γ_{t,j} -- matches the doc's z_t
        definition exactly (not a typo). dv1i, dv2i: per-type, shape (...,ni). dv10, dv20: aggregate
        informal, shape (...,). """
        γi, γi_ = self.get('γi', t), self.get('γi[t-1]', t)
        γ0, γ0_ = self.get('γ0', t), self.get('γ0[t-1]', t)
        ω1i, ω10, ω2i, ω20 = self.ω1i(t), self.ω10(t), self.ω2i(t), self.ω20(t)
        ν = self.get('ν', t)
        old = (γi_*ω2i*dv2i).sum(axis = -1) + γ0_*ω20*dv20
        young = ν*((γi*ω1i*dv1i).sum(axis = -1) + γ0*ω10*dv10)
        return old + young

    # --- Old generations' marginal utilities. dln(h_t)/dτ_t is always an explicit argument: it is closed
    # form ONLY in the terminal period (dlnΘhTerminal_dτ). For t<T, Θ_{h,t} depends on τ_{t+1} =
    # τ^{t+1}(ι_t(τ_t)), an interpolant with no closed-form derivative, so policy.py obtains it as a grid
    # derivative of ln(Θ_{h,t}) instead (docs §PEELOG, "Which derivatives are taken on the grid"). Passing
    # it in rather than computing it here is what keeps a terminal-only formula from being callable at
    # t<T. Both are preference-agnostic (no _LOG suffix): CRRA reuses them with a c^{1-1/ρ} weight.
    def dlnc2i_dτ(self, dlnh_dτ, τ, θ, siRatio_, t = None):
        """ Eq (PEE)/(PEELOG): dln(c_{2,t}^i)/dτ_t, given dln(h_t)/dτ_t. siRatio_ = s_{t-1,i}/s_{t-1}
        (explicit -- predetermined at t, but itself a function of τ_t via si_s at the t-1 vintage).

        This closed form is mandatory, not merely cheaper (docs §PEELOG): siRatio_ varies along a grid of
        candidate τ_t, and the policy maker takes it as predetermined, so a grid derivative would fold in
        a channel that does not belong in the first order condition. """
        α, p_, κ_ = self.get('α', t), self.get('p[t-1]', t), self.get('κ[t-1]', t)
        A0 = (1-α)/α * p_/κ_
        bracket = 1 + self._bcast(θ)*(self.hηRatio(t, lag = '[t-1]') - 1)
        num = self._bcast(A0) * bracket
        denom = siRatio_ + self._bcast(A0*τ)*bracket
        return self._bcast((1-α)*dlnh_dτ) + num/denom

    def dlnc20_dτ(self, dlnh_dτ, τ, ε, ι_, t = None):
        """ Eq (dv20): dln(c_{2,t}^0)/dτ_t = (1-α)dln(h_t)/dτ_t + A_t/(ι_{t-1}+A_tτ_t), with
        A_t = auxInf1_(ε_t). ι_ = ι_{t-1} (explicit -- the political problem's endogenous state).

        Contrast dlnc2i_dτ: here the grid derivative would be legitimate, since ι_{t-1} is a grid
        coordinate held fixed as τ_t varies -- promoting the ratio to a state is what removes the hazard.
        The closed form is used because it is exact and cheaper. The pole at ι_{t-1}+A_tτ_t = 0 is kept
        outside the evaluation region by l_ι>0 with τ_t≥0 (docs §PEELOG), so no masking is needed. """
        A = self.auxInf1_(ε, t)
        return (1-self.get('α', t))*dlnh_dτ + A/(ι_ + A*τ)

    # --- Terminal period (eq:terminalPEELOG). Every object entering z_T is closed form in (τ_T, ι_{T-1}):
    # the terminal young neither save nor face a continuation policy, so no numerical differentiation is
    # needed at t=T at all. dυ_{1,T}^0/dτ_T = 0 has no method -- policy.py passes zeros.
    def dlnΘhTerminal_dτ(self, τ, t = None):
        """ Eq (terminalPEELOG): dln(h_T)/dτ_T = dln(Θ_{h,T})/dτ_T = -ξ/(1+αξ)·1/(1-τ_T). Feeds
        dlnc2i_dτ/dlnc20_dτ's first argument in the terminal period. """
        ξ, α = self.get('ξ', t), self.get('α', t)
        return -1/(1-τ) * ξ/(1+α*ξ)

    def dv1iTerminal_dτ_LOG(self, τ, t = None):
        """ Eq (terminalPEELOG): dυ_{1,T}^i/dτ_T = -(1+ξ)/(1+αξ)·1/(1-τ_T) -- identical across types i.
        Returned as (M,1) (not (M,)) so it broadcasts against the (ni,) political weights inside FOC. """
        ξ, α = self.get('ξ', t), self.get('α', t)
        return self._bcast(-1/(1-τ) * (1+ξ)/(1+α*ξ))

    # --- t<T young generations (eq:v1LOG). Both profiles are known only up to an additive constant along
    # τ_t and are built from the Θ's rather than from consumption levels: by eq:EE:sigma_ci every
    # equilibrium consumption function is Θ_x·(s_{t-1}/ν_t)^{σ_x} with a policy-independent exponent, so
    # ln(x_t) - ln(Θ_x) does not vary with τ_t (eq:logsep) and s_{t-1} drops out entirely. Assemble each
    # profile and differentiate it in ONE numerical step -- differentiating ln(c̃_1) and ln(R_{t+1})
    # separately costs a second spline fit for no gain.
    def v1iProfile_LOG(self, Θh, Θs, Θh1, β, τ1, θ1, Γs, t = None):
        """ Eq (v1LOG): υ_{1,t}^i up to an additive τ_t-constant,
            (1+β_{t,i})ln(tilde-c_{1,t}^i) + β_{t,i}ln(R_{t+1}).
        tildec1i(Θ_{h,t},...) is c̃_{1,t}^i's coefficient function -- see tildec10Θ's docstring. Shapes:
        Θh/Θs/Θh1/τ1 are (M,), β is (ni,) -> returns (M,ni). """
        return (1+β)*np.log(self.tildec1i(Θh, β, τ1, θ1, Γs, t)) \
               + β*self._bcast(self.lnRleadΘ(Θs, Θh1, t))

    def v10Profile_LOG(self, Θs, Θh1, β0, τ1, ε1, t = None):
        """ Eq (v1LOG): υ_{1,t}^0 up to an additive τ_t-constant,
            (1+β_{t,0})ln(tilde-c_{1,t}^0) + β_{t,0}ln(R_{t+1}).
        The informal household's own return enters as ln(R_{t+1}^0) = ln(R_{t+1}) + ln(χ^R_{t+1}) with the
        second term exogenous, hence the same lnRleadΘ as the formal line. Returns (M,). """
        return (1+β0)*np.log(self.tildec10Θ(Θs, β0, τ1, ε1, t)) + β0*self.lnRleadΘ(Θs, Θh1, t)

    #######################################################################
    ##########   10. Calibration targets (eq:calibration)          #########
    #######################################################################
    # Evaluated at the single baseline year db['t0'] -- pass that year as `t`, never rely on the default.
    # calibrationη0/X0 read z_0^η/z_0^x (db['zη0']/db['zx0']) but NOT db['η0']/db['X0']: they compute the
    # *implied* η0/X0 that model.py's calibration loop compares against the current ones. z_0^x is hours
    # relative to AVERAGE FORMAL HOURS ∑γ_i h_{t,i} = h_t·hoursUnitRatio, which equals the aggregate h_t
    # only under eq (calibration:yNorm)'s second normalisation (vector X, hoursUnitRatio = 1). Under the
    # common-X variant that normalisation is spent on the workweek target and the ratio carries the
    # difference (docs eq:calibration:eta0/X0, the factor M).
    def hoursUnitRatio(self, t = None):
        """ M_t = ∑_i γ_{t,i}(η_{t,i}/X_{t,i})^ξ / Γ_{h,t}: average formal hours over the aggregate h_t
        (eq EE:hi summed with population weights). Exactly 1 under eq (calibration:yNorm); under common X
        it is X^{-ξ/(1+ξ)}·∑γ_i(z^η_i)^{ξ/(1+ξ)}. """
        return (self.get('γi', t) * self.hRatio(t)).sum(axis = -1)

    def avgHoursFormal(self, h, t = None):
        """ Average formal hours h̄_t = ∑_i γ_{t,i}h_{t,i} = h_t·M_t -- the object the observed workweek is
        comparable to and the common-X calibration pins (NOT the aggregate h_t, which is in efficiency
        units and does not respond to the hours unit). """
        return h * self.hoursUnitRatio(t)
    def savingsRate(self, s, s_, h, t = None):
        """ Eq (calibration:sr): s_t / ((s_{t-1}/ν_t)^α h_t^{1-α}), savings over FORMAL output. s,
        s_=s_{t-1}, h explicit. Reported, not targeted -- capitalOutputRatio identifies β (see it). The
        GDP-based rate is s_t/(outputFormal + outputInformal). """
        α = self.get('α', t)
        return s / ((s_/self.get('ν', t))**α * h**(1-α))

    # GDP = formal plus informal output, both per formal young household (Σ_i γ_{t,i} = 1): the concept of
    # the tax target (model.taxTarget) and of every "% of GDP" readout.
    def outputFormal(self, s_, h, t = None):
        """ Y_t/n_t = (s_{t-1}/ν_t)^α h_t^{1-α}: the Cobb-Douglas technology on eq:aggregateStates' K_t and
        H_t, per formal young household. s_ = s_{t-1}, h = h_t (explicit). Equals w_t h_t/(1-α). """
        α = self.get('α', t)
        return (s_/self.get('ν', t))**α * h**(1-α)

    def outputInformal(self, s_, t = None):
        """ γ_{t,0} w_t^0 η_{t,0} h_{t,0} = γ_{t,0}(η_0^{1+ξ}/X_0^ξ)(w_t^0)^{1+ξ}: the informal type's linear
        technology (eq:factorPrices0) at its labour supply (eq:informalOpt), per formal young household.
        s_ = s_{t-1} (explicit). """
        return self.get('γ0', t) * self.auxProd0(t) * self.w0(s_, t)**(1+self.get('ξ', t))

    def capitalOutputRatio(self, s_, h, t = None):
        """ Eq (calibration:KY): K_t/Y_t in ANNUAL output units, = yearsPerPeriod·(K_t/h_t)^{1-α} with
        K_t = s_{t-1}/ν_t. s_ = s_{t-1}, h = h_t.

        The factor is not cosmetic: Y_t is one 30-year period's output, so K_t/Y_t is a *thirtieth* of the
        capital-output ratio the data report and the two are not comparable without it. Equivalent to the
        US arm's R_{t0} target by K/Y = yearsPerPeriod·α/R (Cobb-Douglas), which is why this target needs
        no s_t: it prices the predetermined stock, not the flow into it. """
        α = self.get('α', t)
        return self.db['yearsPerPeriod'] * (s_/self.get('ν', t)/h)**(1-α)

    def calibrationη0(self, Θh, τ, t = None):
        """ Eq (calibration:eta0): η_0 implied by z_0^η/z_0^x at a given Θ_{h,t}, τ_t, with z_0^x read
        against average formal hours h_t·M_t (hoursUnitRatio; M = 1 under vector X). """
        α, ξ, Γh = self.get('α', t), self.get('ξ', t), self.Γh(t)
        M = self.hoursUnitRatio(t)
        return (self.get('zη0', t)/self.get('zx0', t)) * (1-α)*(1-τ) / (M * Θh**α * ((1-α)/Γh**α)**(1/(1+α*ξ)))

    def calibrationX0(self, η0, Θh, t = None):
        """ Eq (calibration:X0): X_0 implied by a given η_0, Θ_{h,t}. Feed it calibrationη0's output.
        Informal hours are z_0^x·M_t·h_t, so M sits inside the 1/ξ power. """
        α, ξ, Γh = self.get('α', t), self.get('ξ', t), self.Γh(t)
        M = self.hoursUnitRatio(t)
        return η0 * ((1-α)/Γh**α)**(1/(1+α*ξ)) / (Θh*self.get('zx0', t)*M)**(1/ξ)


class BaseGrid(Base):
    """ Assumes parameters from database are selected for a single year,
        but inputs/arguments are grids of the same length. """

class BaseTime(Base):
    """ Assume parameters from database are functions of time and
    so are inputs/arguments. """

    def _year(self, t = None):
        """ No year dimension here -- symbols are returned as whole time paths, so a lookup's result does
        not depend on t at all and every cache key collapses to the symbol (plus lag). """
        return None

    def __call__(self, k, t = None):
        """ Return the full time-indexed symbol 'k' (no year-slicing) -- parameters and inputs are vectors over t here. """
        return self._memo(('call', k, None), lambda: self.db[f'{k}'])

    def Γh(self, t = None, lag = ''):
        """ Matches Base.Γh's return type (plain ndarray, not a pandas Series) -- several formulas
        (Θh, ΘhTerminal, Θs, w0, c10/tildec10, c20/tildec20) use self.Γh(t) directly without going
        through self._bcast, so leaving this as a Series would silently break positional indexing
        (e.g. `result[-1]`) on their output whenever called via BaseTime. Memoised inside cacheParams()
        as in Base -- with its own key, since the value differs (a full path, not one year's scalar). """
        return self._memo(('Γh', lag, None), lambda:
            (self(f'γi{lag}', t) * self.auxProd(t, lag = lag)).sum(axis = 1).values)

    #######################################################################
    ##########   8. Finite-horizon terminal stacking (FH_*)         ########
    #######################################################################
    # Several quantities use a *different* formula at the terminal period (t=T-1) than for t<T-1 (e.g.
    # labour supply: Θh, needing Γs_t/τ_{t+1}/θ_{t+1}, vs. ΘhTerminal). These FH_* methods own that
    # stacking once rather than repeating it at every call site (EE_LOG_solve, EE_CRRA_*, EE_report).
    # Convention: Γs/B are length T-1 (genuinely undefined at T, no period T+1 to look into); h/s/τ/τ1/θ1
    # are length T. Where the terminal formula is a special case of the general one (c1i, tildec1i --
    # feeding B=0/Γs=0 alongside s_T=0 makes the smoothing/pension terms vanish, verified algebraically),
    # padding replaces branching; where it isn't (FH_h), the terminal entry is computed separately.
    def FH_h(self, τ, τ1, θ1, Γs, s_, t = None):
        """ Labour supply h_t, stacking Θh (eq:auxiliary:Thetah, t<T-1) with ΘhTerminal
        (eq:auxiliary:ThetahT, t=T-1). Γs: length T-1. τ, τ1, θ1, s_ (=s_{t-1} for every t): length T.
        Returns length-T array. """
        Γs_pad = np.append(Γs, 0) # dummy -- Θh's terminal-index output is overwritten below, never used
        Θh = self.Θh(τ, τ1, θ1, Γs_pad, t)
        Θh[-1] = self.ΘhTerminal(τ, t)[-1]
        return self.h(Θh, s_, t)

    def FH_c1i(self, h, s, B, τ1, θ1, t = None):
        """ Eq (EE:ci) stacked with the terminal formula (§Finite Horizon Terminal state) at t=T-1:
        c1_T^i = auxProd_i(h_T/Γh_T)^{(1+ξ)/ξ}, with no pension/smoothing term (no period T+1 to smooth
        into). Padding B with 0 and relying on the finite-horizon s[-1]=0 makes c1i's general formula
        collapse to exactly this. h, s, τ1, θ1: length T. B: length T-1. """
        B_pad = np.vstack([B, np.zeros((1, B.shape[1]))])
        return self.c1i(h, s, B_pad, τ1, θ1, t)

    def FH_tildec1i(self, h, B, τ1, θ1, Γs, t = None):
        """ As FH_c1i, for tildec1i (see FH_c1i's docstring for the terminal-collapse argument -- here
        driven by padding both B and Γs with 0). h, τ1, θ1: length T. B, Γs: length T-1. """
        B_pad = np.vstack([B, np.zeros((1, B.shape[1]))])
        Γs_pad = np.append(Γs, 0)
        return self.tildec1i(h, B_pad, τ1, θ1, Γs_pad, t)

    def FH_c10(self, s_, s, B0, τ1, ε1, t = None):
        """ Eq (EE:c0) stacked with the terminal formula (§Finite Horizon Terminal state) at t=T-1:
        c_{1,T}^0 = auxProd0·(w_T^0)^{1+ξ} (informal young no longer save either, so no pension/smoothing
        term). Padding B0 with 0 alongside the finite-horizon s[-1]=0 makes c10's general formula collapse
        to exactly this. s_, s, τ1, ε1: length T. B0: length T-1. """
        return self.c10(s_, s, np.append(B0, 0), τ1, ε1, t)

    def FH_tildec10(self, s_, s, B0, τ1, ε1, t = None):
        """ As FH_c10, for tildec10 (terminal value c_{1,T}^0/(1+ξ)). """
        return self.tildec10(s_, s, np.append(B0, 0), τ1, ε1, t)

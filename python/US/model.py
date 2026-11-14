import numpy as np, pandas as pd, functools, time
from scipy import optimize
from copy import deepcopy
from symMaps import SimpleSys, Lag
from base import Base, BaseGrid, BaseTime
from policy import LOG, CRRA
from gridsearch import continuation


def _shiftT(idx, t0, tName):
    """ Row-wise: subtract t0 from the `tName` axis/level of a pandas Index/MultiIndex, leaving any
    other level's values untouched. """
    if isinstance(idx, pd.MultiIndex):
        arrays = [idx.get_level_values(n) - t0 if n == tName else idx.get_level_values(n) for n in idx.names]
        return pd.MultiIndex.from_arrays(arrays, names = idx.names)
    return idx - t0


def _sliceDb(db, t0, tName = 't'):
    """ Restrict every db entry indexed (wholly or via one level) by `tName` to periods from t0 onward,
    AND renumber that axis back to 0-based (subtract t0), mutating db IN PLACE.

    Renumbering (not just restricting) is required, not cosmetic: model.py's EE_LOG_solve/
    EE_CRRA_solve/EE_report/initialState_solve index the plain ndarrays a caller passes in (τ, θ, ε)
    positionally via self.B.tFirst (e.g. `τ[self.B.tFirst]`), which is only correct when db['t'] is the
    native 0-based range() a fresh instance's initIdxs() builds -- so tFirst must come back out as 0 on
    the copy. Correspondingly, a caller's own τ/θ/ε arrays for the copy must themselves start at
    position 0 = the original t0 (i.e. sliced off the front) -- solvePEE_LOG/CRRA's own θ=None/ε=None
    defaults already get this for free, since db['θ']/db['eps'] are sliced+renumbered here too.

    createCopyFromt0 relies on self.B/self.BG/self.BT/self.LOG/self.CRRA all sharing the same dict
    object post-deepcopy -- rebinding db to a fresh dict here would silently orphan that aliasing, so
    this mutates the existing dict's entries in place rather than returning a new one. Scalars,
    type-indexed (j/i/u) objects, and eigenvector-calibration arrays pass through unchanged. Covers
    db['t']/db['txE'] (plain Index), db['tj']/db['ti'] (MultiIndex), and every 1D/2D parameter plus its
    [t+1]/[t-1] siblings (Series/DataFrame) in one pass. """
    for k, v in db.items():
        if isinstance(v, pd.MultiIndex):
            if tName in v.names:
                db[k] = _shiftT(v[v.get_level_values(tName) >= t0], t0, tName)
        elif isinstance(v, pd.Index):
            if v.name == tName:
                db[k] = v[v >= t0] - t0
        elif isinstance(v, (pd.Series, pd.DataFrame)):
            idx = v.index
            if isinstance(idx, pd.MultiIndex):
                if tName in idx.names:
                    sub = v[idx.get_level_values(tName) >= t0]
                    db[k] = sub.set_axis(_shiftT(sub.index, t0, tName))
            elif idx.name == tName:
                sub = v[idx >= t0]
                db[k] = sub.set_axis(sub.index - t0)


class ModelUS:
    _BaseClass = Base
    _BaseGridClass = BaseGrid
    _BaseTimeClass = BaseTime
    _PolicyLogClass = LOG
    _PolicyCRRAClass = CRRA

    def __init__(self, nj = 4, T = 10, pars = None, commonX = False, **kwargs):
        """ US model (docs writing/US/): formal households only.

        The informal type j=0 has NO mass here (γ_0 = 0) -- the docs' model has no such household at all.
        The slot is nonetheless retained throughout db, base.py and policy.py, because the data are laid
        out that way (data/README.md) and every array is j-indexed. With γ_0 = 0 the
        informal terms are computed and then multiplied by zero: κ_t collapses to p_t (Base.κ), and the
        informal blocs drop out of the political FOC (Base.FOC). Synthetic η_0/X_0 keep those dead terms
        finite -- they must stay finite, since 0·NaN is NaN and would poison the whole FOC.

        Parameters
        ----------
        nj: Integer - number of household types, INCLUDING the zero-mass slot j=0.
        T : Integer - years.
        commonX: bool - which calibration variant of the docs' §calibration to use.
            False (default, "vector X", the paper's robustness check): a full vector X_i, identified from relative-hours data z^x_i
                via the eigenvector system. The hours unit is left arbitrary, so the average workweek
                h̄ has no meaningful level (see calibrate's docstring).
            True ("common X", the paper's, python/paper/config.py US['commonX']): a common scalar X across types, identified from the level of average
                hours db['h0']. Relative hours then become a model *prediction*, not data.

        db : dict
            Keys = Variable names (string).
            Values = Data either stored as scalars, pandas series, or numpy arrays.
        """
        self.nj, self.T = nj, T
        self.commonX = commonX
        self.ni = self.nj-1
        self.db = {} # dictionary database
        self.x0 = {} # cache of last-used/default initial guesses for numerical problems, keyed by problem name
        self.addProperty('paramsFromFuncs', ['Γh','θ','eps','κ']) # list of auxiliary parameters to be recalculated when deep parameters are updated
        self.parTypes = self._parTypes.copy() # Mapping of parameter names to dimensions
        self.initIdxs() # Initialize relevant pandas indices
        self.addNamespaces() # aux classes that help nagivate stacked vectors + lag/lead symbols.
        self.B = self._BaseClass(self) # Class used for scalar computations
        self.BG = self._BaseGridClass(self)  # Class used for vectorized computations across states
        self.BT = self._BaseTimeClass(self)  # Class used for vectorized computations across time
        self.initPars(pars = pars)
        self.initProductivity()
        self.updateAuxPars()
        self.LOG = self._PolicyLogClass(self) 
        self.CRRA = self._PolicyCRRAClass(self)

    #######################################################################
    ##########                0.1. Auxiliary methods            ###########
    #######################################################################

    # Some basic methods for navigating symbols:
    def leadSym(self, symbol, lead = -1, opt = None, ns = 'exo'):
        """ 'lead' uses the same sign convention as Lag.series' 'lag' argument, i.e. result.loc[t] = symbol.loc[t-lead].
        bfill = 'ss' clamps the out-of-domain boundary period to the nearest in-domain value (replaces the old useLoc='nn' behavior). """
        o = {'bfill': 'ss'} if opt is None else opt
        if isinstance(symbol, pd.Series):
            return Lag.series(symbol, lead, **o)
        elif isinstance(symbol, pd.DataFrame):
            return Lag.series(symbol.stack(), lead, **o).unstack()
        elif isinstance(symbol, np.ndarray):
            if symbol.ndim == 1:
                return np.hstack([symbol[1:], symbol[-1]])
            elif symbol.ndim == 2:
                return np.vstack([symbol[1:], symbol[-1]])

    def __call__(self, x, name, ns = 'exo', **kwargs):
        return self.ns[ns](x, name, **kwargs)

    def get(self, x, name, ns = 'exo'):
        return self.ns[ns].get(x, name)

    ### Specify parameter structure - what parameters are 1D, 2D - what auxiliary parameters do we have and how do we recompute them? 
    @property
    def _parTypes(self):
        return {'0D': list(self.default0DParams), '1D': list(self.default1Dparams)+self.aux1DParams+self.paramsFromFuncs, '2D': list(self.default2Dparams)+self.aux2DParams}
    @property
    def default0DParams(self):
        """ scalars. Calibration targets, all matched at db['t0']:
        - 'τ0'      the equivalent social security tax rate implied by observed pension spending.
        - 'R0'      the 30-year gross interest rate (docs eq:calibration:R). This REPLACES the Argentina
                    models' savings-rate target 's0'; the savings rate is still reported, just not targeted.
        - 'h0'      the average workweek as a share of discretionary time (docs eq:calibration:Xsolve).
                    Only a target when commonX=True; ignored otherwise. Set by test.py, not defaulted here.
        - 'RR0'     the OECD ratio of replacement rates at mean vs half-mean income, which identifies θ_P,
                    the design of the pension formula.
        - 'Ushare0' the share of transfers to the old paid as a flat benefit outside the pension formula:
                    θ = (1 - Ushare0)·θ_P (getθ). 0 in every workbook. Also the inert ε (getEps).
        The defaults are USMain_test.xlsx's (pension cash spending 7.366% of GDP over 1-α = 0.7; 39.2/49.6).

        NOTE 't0' here is the *index into db['t'] of the calibration baseline year* (2020 for the US;
        defaults to 2, an arbitrary early index) -- NOT self.B.tFirst (=db['t'][0], the model's first
        active/endogenous period, docs' t=1) and NOT the pre-determined period before the horizon (docs'
        t=0). All three are genuinely different things. """
        return {'τ0': .105229, 'RR0': 0.790323, 'R0': 2.443, 'Ushare0': 0., 'RRGroups': (1,2), 't0': 2}

    @property
    def default1Dparams(self):
        """ 1d vectors - defined over t """
        return {'α': .43, 'ν': 1, 'ξ' : .3, 'ρ': 1, 'ω': 1.25, 'α0': 1, 'χ': 1}
    @property
    def default2Dparams(self):
        """ 2d matrices - defined over (t, j)"""
        return {'γj': np.full(self.nj, 1/self.ni)} | {k: np.full(self.nj, 1) for k in ('pj','μj','Xj','ηj','zxj','zηj','βj')}
    @property
    def aux2DParams(self):
        """ 2d matrices - defined over (t, i)"""
        return [f"{k[:-1]}i" for k in self.default2Dparams]
    @property
    def aux1DParams(self):
        """ 1d vectors - defined over t """
        return [f"{k[:-1]}0" for k in self.default2Dparams]

    ### Methods to compute auxiliary parameters used throughout the code:
    @property
    def aux_θ(self):
        return pd.Series(self.getθ(), index = self.db['t'])
    @property
    def aux_eps(self):
        return pd.Series(self.getEps(), index = self.db['t'])
    @property
    def aux_κ(self):
        return self.BT.κ(self.db['eps[t+1]'])
    @property
    def aux_Γh(self):
        return self.BT.Γh()

    def addProperty(self, key, value):
        """ default dynamic property method """
        setattr(self, f'_{key}', value)
        setattr(type(self), key, property(fget = lambda self: getattr(self, f'_{key}'), fset = lambda self, value: setattr(self, f'_{key}', value)))

    def __setstate__(self, state):
        self.__dict__.update(state)
        value = getattr(self, '_paramsFromFuncs', getattr(self, '_paramFromFuncs', None))
        if value is not None:
            self.addProperty('paramsFromFuncs', value)

    def initIdxs(self):
        self.db['t'] = pd.Index(range(self.T), name = 't')
        self.db['txE'] = pd.Index(range(self.T-1), name = 't') # Time index without terminal period
        self.db['j'] = pd.Index(range(self.nj), name = 'j')
        self.db['i'] = self.db['j'][1:]
        self.db['u'] = self.db['j'][0:1]
        self.db['tj'] = pd.MultiIndex.from_product([self.db['t'], self.db['j']])
        self.db['ti'] = pd.MultiIndex.from_product([self.db['t'], self.db['i']])

    def addNamespaces(self):
        self.ns = {}
        self.ns['exo'] = SimpleSys(v = (dict.fromkeys(self.default1Dparams, self.db['t']) |
                                       dict.fromkeys(self.default2Dparams, self.db['tj'])|
                                       dict.fromkeys(self.aux1DParams, self.db['t']) |
                                       dict.fromkeys(self.aux2DParams, self.db['ti']) |
                                       dict.fromkeys(self.paramsFromFuncs, self.db['t'])))
        # One namespace per numerical problem (§3): here just for stacking/unstacking EE_CRRA's unknown
        # vector x=(Γs,h,s) (all txE-domain) -- NOT used for its x0()/db-array machinery (that's tied to
        # symMaps' own SimpleDB, not our plain-dict self.db), only for __call__/get/unloadSol's name<->
        # position mapping, so this coexists cleanly with keeping solution objects out of self.db.
        self.ns['EE_CRRA'] = SimpleSys(v = dict.fromkeys(('Γs', 'h', 's'), self.db['txE']))
        [ns.compile() for ns in self.ns.values()];

    def initPars(self, pars = None):
        self.db.update(self.defaultParameters) # default parameters and targets
        self.addDefaultHeterogeneity # default heterogeneity
        if pars is not None:
            [self.db.update(self.adjPar(k,v)) for k,v in pars.items()];

    ### Methods used to navigate 0D/1D/2D parameters including lags/leads - and how to update full dimension parameters based inputs with fewer dimensions. 
    @property
    def defaultParameters(self):
        """ Return full dictionary and default values
        and with dimensions broadcasted to appropriate dimensions for 0D/1D """
        return functools.reduce(lambda x,y: x|y, [self.adjPar(k,v) for k,v in self.default0DParams.items()] + [self.adjPar(k,v) for k,v in self.default1Dparams.items()])
    @property
    def addDefaultHeterogeneity(self):
        """ Return dictionary with defautlt values 
        and with dimensions broadcasted to appropriate dimensions for 2D. """
        [self.db.update(self.adjPar(k,v)) for k,v in self.default2Dparams.items()]; 

    def updateAuxPars(self):
        [self.db.update(self.addLeadAndLags(k, getattr(self, f'aux_{k}'))) for k in self.paramsFromFuncs]

    def adjPar(self, k, vals):
        """ The db entries for parameter k at vals, broadcast over the whole horizon (a scalar, a (nj,)/(ni,)
        type vector or a full (T,)/(T, n) array), with their [t+1]/[t-1] siblings. Returns new objects;
        db itself is not touched. """
        if k == 'β':
            return self.adj2Dj_ll('βj', self.simpleβ(vals))
        elif k == 'pj':
            return self.adjpj('pj', vals)
        elif k in self.parTypes['0D']:
            return {k: vals}
        elif k in self.parTypes['1D']:
            return self.adj1D_ll(k, vals)
        elif k in self.aux2DParams:
            return self.adj2Di_ll(k, vals)
        elif k in self.parTypes['2D']:
            return self.adj2Dj_ll(k, vals)
        else:
            return {k:vals}
    def simpleβ(self, β):
        """ Override full βj with scalar """
        return β * self.db['pj']
    def simpleβinv(self):
        """ Get scalar β from full dataframe (assumes constant across t,j) """
        return self.db['βj'].iloc[0,0]/self.db['pj'].iloc[0,0]
    def adjpj(self, k, vals):
        """ Specialized method for p - we have to enforce constant p across pi"""
        d = self.adjDf_tj(k, vals)
        d['p'] = (self.db['γi'] * d['pi']).sum(axis=1)
        return self.addLeadAndLags(d)

    def adj1D_ll(self, k, vals):
        return self.addLeadAndLags(k, self.adjVec_t(vals))
    def adj2Dj_ll(self, k , vals):
        return self.addLeadAndLags(self.adjDf_tj(k, vals))
    def adj2Di_ll(self, k , vals):
        return self.addLeadAndLags(self.adjDf_ti(k, vals))

    def addLeadAndLags(self, k, s = None):
        if isinstance(k, dict):
            return functools.reduce(lambda x,y: x|y, [self.addLeadAndLags(i,v) for i,v in k.items()])
        else:
            return {k: s, f'{k}[t+1]': self.leadSym(s, lead = {'t':-1} if isinstance(s, pd.DataFrame) else -1), f'{k}[t-1]': self.leadSym(s, lead = {'t': 1} if isinstance(s, pd.DataFrame) else  1)}

    def adjVec_t(self, vals):
        return pd.Series(vals, index = self.db['t'])
    def adjDf_ti(self, k, vals):
        if np.isscalar(vals):
            xi = pd.DataFrame(vals, index = self.db['t'], columns = self.db['i'])
        elif vals.ndim == 1:
            xi = pd.DataFrame(np.tile(vals, (self.T,1)), index = self.db['t'], columns = self.db['i'])
        else:
            xi = pd.DataFrame(vals, index = self.db['t'], columns = self.db['i'])
        return {f'{k[:-1]}j': pd.concat([self.db[f'{k[:-1]}0'], xi], axis = 1), k: xi}
    def adjDf_tj(self, k, vals):
        if np.isscalar(vals):
            xj = pd.DataFrame(vals, index = self.db['t'], columns = self.db['j'])
        elif vals.ndim == 1:
            xj =  pd.DataFrame(np.tile(vals, (self.T, 1)), index = self.db['t'], columns = self.db['j'])
        else:
            xj = pd.DataFrame(vals, index = self.db['t'], columns = self.db['j'])
        return {k: xj, f'{k[:-1]}i': xj[self.db['i']], f'{k[:-1]}0': xj[0]}



    #######################################################################
    ##########            2. Simple calibration methods         ###########
    #######################################################################

    def θOverride(self):
        """ db['θ0'], the design imposed as a datum (a variant of data/oecd26_theta.csv, oecd26.thetaFromFile),
        as a float; None where the key is absent or not finite. Outside [0, 1] it raises. """
        v = self.db.get('θ0')
        if v is None or not np.isfinite(float(v)):
            return None
        if not 0. <= float(v) <= 1.:
            raise ValueError(f"db['θ0'] = {v!r} is outside [0, 1]")
        return float(v)

    def getθ(self):
        """ θ from the replacement-rate ratio RR0 of the groups RRGroups (calibration:theta): θ_P, the design
        of the pension formula, times 1 - Ushare0. A flat benefit paid outside the formula, a share m of all
        transfers to the old, makes the transfers' design (1 - m)θ_P exactly under the benefit formula
        (eq:governmentBudget). Ushare0 = 0 returns θ_P bit for bit.
        With db['θ0'] present and finite (θOverride) θ is db['θ0'] itself, which Ushare0 does not scale; RR0 is
        then not read. Every refresh re-installs it (#9). """
        θ0 = self.θOverride()
        if θ0 is not None:
            return θ0
        i,ii = self.db['RRGroups'][0], self.db['RRGroups'][1]
        ξ= self.db['ξ'].xs(self.db['t0'])
        h1,h2 = (self.db['Γh']*self.db['Xi'][i]**ξ/self.db['ηi'][i]**(1+ξ)).xs(self.db['t0']), (self.db['Γh']*self.db['Xi'][ii]**ξ/self.db['ηi'][ii]**(1+ξ)).xs(self.db['t0'])
        θP = (self.db['RR0']*h1-h2)/(1-h2-self.db['RR0']*(1-h1))
        return (1-self.db['Ushare0'])*θP
    def getEps(self):
        """ ε, the universal-pension parameter of the (zero-mass) informal type.

        INERT in this model. ε reaches the equilibrium through exactly two channels -- κ_t (Base.κ, where
        it appears only as ε·γ_0·p_0) and the informal old household's marginal utility (Base.dlnc20_dτ,
        whose contribution to the FOC carries the factor γ_{t-1,0}) -- and γ_0 = 0 kills both. We set it
        to db['Ushare0'] rather than to zero because that is what the name means in this codebase, and a
        zero would be indistinguishable from "not yet decided" if the type were ever given mass.

        Which is the case this guards against: db['Ushare0'] is the universal (non-earnings-related) share
        of US transfers to the old (it scales θ, getθ), and it is NOT the same object as the Argentina models' ε (a coverage-adjusted,
        discounted benefit ratio for a group that actually receives pensions). If γ_0 is ever made
        positive -- the platform-work extension the paper flags -- ε has to be derived properly first,
        so fail loudly here rather than let a plausible-looking number in through the back door. """
        if float(np.max(self.db['γ0'])) != 0.:
            raise NotImplementedError(
                "getEps: γ_0 > 0, but db['Ushare0'] is the universal share of transfers, not a calibrated ε. "
                "Derive ε for the informal type before giving it mass (see writing/US/model_calibration.tex).")
        return self.db['Ushare0']

    # Initialize productivity distribution:
    def initProductivity(self):
        """ Get ηi, Xi from data on relative income (and, under vector X, relative hours).

        Both variants normalize Γ_h = 1 (docs eq:scaleInvariance -- that normalization is free, and is the
        same one the Argentina models use). They differ only in how the hours unit (eq:hoursUnit) is fixed;
        see __init__'s commonX. Everything is keyed off the baseline year db['t0'] and broadcast across t,
        exactly as in the Argentina models. """
        if self.commonX:
            self.initProductivity_commonX()
        else:
            self.initProductivity_vectorX()

    def initProductivity_vectorX(self):
        """ Vector X: a full vector X_i, identified from relative income AND relative hours via the two
        eigenvector systems (docs eq:calibration:etai/Xi). Relative-hours data cannot speak to the hours
        unit ϖ = ∑γ_i y^x_i (eq:hoursUnit), so addEigenVectors sets it to 1 -- h̄_{t0} = h_{t0} in the
        unit the calibration reports, the same convention in every country. """
        self.addEigenVectors()
        ηi = self.getηi()
        η0 = 0.3 * ηi[0] * self.db['zη0'].xs(self.db['t0'])/self.db['zηi'].xs(self.db['t0'])[1] # initial guess for η0 based on the rest of the vector
        self.db.update(self.adjPar('ηj', np.hstack([η0,ηi])))
        xj = self.db['zxj'].xs(self.db['t0'])
        yx = np.hstack([self.db['yx'][0]*xj[0]/xj[1], self.db['yx']]) # inital guess for yx vector.
        self.db.update(self.adjPar('Xj', self.db['ηj'].values/(np.tile(yx, (self.T,1))**(1/self.db['ξ'].values.reshape(self.T,1)))))

    def zηiNormalized(self):
        """ z^η_i rescaled so that ∑_i γ_i z^η_i = 1, i.e. income relative to the POPULATION-WEIGHTED
        average. test.py divides raw income by its unweighted mean (and over j, including the zero-mass
        slot), which is not the same normalization the model's z^η means. Vector X never notices -- the
        eigenvector system is scale-free and getηi divides the scale back out -- but common X inverts
        z^η directly and would silently mis-scale η without this. Weighted at db['t0']. """
        t0 = self.db['t0']
        zηi, γi = self.db['zηi'].xs(t0).values, self.db['γi'].xs(t0).values
        return zηi/(γi*zηi).sum()

    def initProductivity_commonX(self, X = 1.):
        """ Common X: a common scalar X across types. With Γ_h = 1 the income condition inverts in closed
        form (docs eq:calibration:etaCommonX), η_{t,i} = (z^η_i)^{1/(1+ξ_t)}·X^{ξ_t/(1+ξ_t)}, and no
        eigenvector solve is needed -- with a common X the two systems collapse into one.

        Γ_h = 1 then holds for ANY X (because ∑γ_i z^η_i = 1 by construction of zηiNormalized), which is
        exactly why X is free to be calibrated to the average workweek afterwards without disturbing
        anything else. Relative hours are no longer data here: the model predicts h_i/h ∝ (z^η_i)^{ξ/(1+ξ)},
        so db['zxi'] goes unused and is worth checking against as a diagnostic (see calibrate).

        The zero-mass slot j=0 is cloned off type 1 -- any finite positive value does, since γ_0 = 0
        multiplies it away; what matters is only that it is finite (0·NaN would poison the FOC). """
        ξ = self.db['ξ'].values.reshape(self.T, 1)
        ηi = self.zηiNormalized().reshape(1, self.ni)**(1/(1+ξ)) * X**(ξ/(1+ξ))
        self.db.update(self.adjPar('ηj', np.hstack([ηi[:, :1], ηi])))
        self.db.update(self.adjPar('Xj', np.full((self.T, self.nj), X)))

    def hbarTarget(self):
        """ The value h̄_{t0} is calibrated to. db['h0'] as given -- the observed workweek in the model's
        own time unit. Read through this method rather than db['h0'] directly, because a level of h̄ is
        only meaningful under commonX: modelFR.ModelFR overrides it with a target defined RELATIVE to a
        US reference, which is the only form that means anything under vector X (docs eq:us:model:hoursUnit:
        h̄ is reported against a reference, never as a level). """
        return self.db['h0']

    def solveCommonX(self, h):
        """ Docs eq:calibration:Xsolve: the X at which average hours at db['t0'] hit hbarTarget(),
        given the aggregate h_{t0} an already-solved equilibrium produced.

        Closed form, and exact -- not an approximation to be iterated. X enters no aggregate (see
        initProductivity_commonX), so h_{t0} on the right hand side does not move when X is applied. """
        t0 = self.db['t0']
        ξ, γi = self.db['ξ'].xs(t0), self.db['γi'].xs(t0).values
        S = (γi * self.zηiNormalized()**(ξ/(1+ξ))).sum()
        return (h * S/self.hbarTarget())**((1+ξ)/ξ)

    def rescaleX(self, λ):
        """ The scale invariance (docs eq:us:model:scaleInvariance) applied as a calibration step:
        X_{t,i} -> λ^{-1/ξ} X_{t,i} at fixed η, so Γ_h -> λΓ_h and every LEVEL (h, h̄, h_i, s, c, Y)
        scales by λ while R, w, τ, θ, s_i/s and the savings rate do not move at all.

        Note it moves all X_i by the SAME proportion (λ^{-1/ξ}), so it is not the hours-unit rescaling
        of eq:us:model:hoursUnit -- that one holds y^η fixed and leaves aggregate h alone. Here Γ_h ends
        at λ, so the Γ_h = 1 normalisation both initProductivity_* impose is deliberately given up; that
        is the price of putting h̄ on a scale, and it costs nothing because no target reads Γ_h.

        Identical in effect to test_invariance.py's rescale(λ), which is what asserts the invariance. """
        ξ = self.db['ξ'].values.reshape(self.T, 1)
        self.db.update(self.adjPar('Xj', self.db['Xj'].values * λ**(-1/ξ)))
        self.updateAuxPars() # Γh moves to λ; θ is invariant but is recomputed from η/X, so refresh it too

    def addEigenVectors(self):
        """ The eigenvectors y^η, y^x of Eq (calibration:etai/Xi), each scaled to γ·y = 1: Γ_h = 1 AND the
        hours unit ϖ = ∑γ_i y^x_i = 1 (eq:hoursUnit). The matrix z·γᵀ has rank one, so its eigenvector is
        z itself (eigenvalue γ·z) and the scaled eigenvector is y = z/(γ·z), in closed form.

        The y^x scale is not cosmetic: it is the hours unit of eq (hoursUnit), the same in every country
        only under this normalisation. Equilibrium-neutral (h̄ and h_i are the only objects that respond
        to μ), but the reported η_i, X_i and X̄ are not. """
        t0 = self.db['t0']
        γi = self.db['γi'].xs(t0).values
        zx, zη = self.db['zxi'].xs(t0).values, self.db['zηi'].xs(t0).values
        self.db['yx'], self.db['yη'] = zx/(γi*zx).sum(), zη/(γi*zη).sum()
    def getηi(self):
        return self.db['yη']/(self.db['yx']*sum(self.db['γi'].xs(self.db['t0']).values*self.db['yη']))

    #######################################################################
    ##########   3. Economic equilibrium (EE) solve, given policy   ########
    #######################################################################
    # Each problem is a residual, a solve and a report. Γs/B/si_s live on db['txE'] (length T-1, no period
    # T+1 to look into), everything else on db['t'] (length T). EE_LOG_solve/EE_CRRA_solve/EE_report work
    # on plain ndarrays and wrap into pd.Series/DataFrame only in the final return (_wrapVars, keyed off
    # _t2vars/_txE2vars below) -- not inside base.py's own methods (incl. FH_*), since those run inside
    # EE_CRRA_solve's root-finding loop where pandas' per-call overhead would be pure cost.

    # Which pandas index a symbol reports over -- a fixed, once-per-model mapping.
    _t2vars = ('s', 'h', 's_', 'R', 'w', 'w0', 'hi', 'bbar', 'bi', 'b0',
               'c1i', 'tildec1i', 'c2i', 'tildec2i', 'c10', 'tildec10', 'c20', 'tildec20')
    _txE2vars = ('Γs', 'B', 'si_s')

    def _wrapEE(self, x, idx):
        """ Wrap an EE_* result as pd.Series (1D) or pd.DataFrame (2D, columns=db['i']) indexed by idx
        (db['t'] for T-length, db['txE'] for T-1-length quantities). """
        return pd.DataFrame(x, index = idx, columns = self.db['i']) if x.ndim == 2 else pd.Series(x, index = idx)

    def _wrapVars(self, d):
        """ Wrap {name: raw ndarray} into pd.Series/DataFrame via _wrapEE, looking up each name's index
        (db['t'] vs db['txE']) from _t2vars/_txE2vars. Any name in neither is a bug -- fails loudly
        rather than silently mis-indexing a symbol someone forgot to register. """
        t, txE = self.db['t'], self.db['txE']
        def idx(k):
            if k in self._t2vars: return t
            elif k in self._txE2vars: return txE
            else: raise KeyError(f"'{k}' is not registered in _t2vars/_txE2vars -- add it there.")
        return {k: self._wrapEE(v, idx(k)) for k, v in d.items()}

    def _checkConverged(self, residual, tol = 1e-8, name = '', scipyRes = None):
        """ Verify convergence via the actual residual norm rather than trusting a solver's own success
        flag (scipy's res.success behaves inconsistently across methods -- treat as informational only).
        scipyRes: optional OptimizeResult, included in the error message for diagnostic context. Not
        scipy-specific (just needs a residual array), so it generalizes to any numerical problem here. """
        maxResid = np.max(np.abs(residual))
        if not (maxResid <= tol):  # not `maxResid > tol`: NaN comparisons are always False either way, so
                                    # that form would silently let a NaN residual pass as "converged".
            scipyInfo = f" (scipy success={scipyRes.success}, message={scipyRes.message!r})" if scipyRes is not None else ""
            raise RuntimeError(f"{name} did not converge: max|residual|={maxResid:.3e} > tol={tol:.1e}{scipyInfo}")
        return maxResid

    def EE_LOG_solve(self, τ, θ, ε, s0 = None):
        """ Closed-form economic equilibrium given LOG preferences (ρ=1, so B_{t+1}^i=β_{t,i} is a pure
        primitive -- no root-finding). τ, θ, ε: full length-T policy paths. s0: the state before db['t'][0]
        (docs' s_0); defaults to the LOG steady state at db['t'][0]'s (τ,θ) (§4's steadyState_LOG_solve).
        Γs/Θh/Θs are fully exogenous given β; the only genuine recursion is s_t (nonlinear in s_{t-1}), so
        that's the only loop -- h then follows in one vectorized FH_h call. Returns {'s','h','Γs','B'}
        (Γs/B length T-1, h/s length T). No self.x0 caching -- closed-form, nothing to warm-start. """
        if s0 is None:
            s0 = self.steadyState_LOG_solve(τ[self.B.tFirst], θ[self.B.tFirst], t = self.B.tFirst)['s']
        τ1, θ1 = self.leadSym(τ), self.leadSym(θ)
        Γs_full = self.BT.Γs(self.db['βi'].values, τ1, θ1)
        Θh_full = self.BT.Θh(τ, τ1, θ1, Γs_full)
        Θh_full[-1] = self.BT.ΘhTerminal(τ)[-1]
        Θs = self.BT.Θs(Θh_full, Γs_full)
        ν, powerS = self.db['ν'].values, self.BT.power_s()
        s = np.empty(self.T)
        s_ = s0
        for i in range(self.T-1):
            s[i] = Θs[i] * (s_/ν[i])**powerS[i]
            s_ = s[i]
        s[-1] = 0
        s_lag = np.append(s0, s[:-1])
        h = self.BT.FH_h(τ, τ1, θ1, Γs_full[:-1], s_lag)
        return self._wrapVars({'s': s, 'h': h, 'Γs': Γs_full[:-1], 'B': self.db['βi'].values[:-1]})

    def defaultX0_EE_CRRA(self, τ, θ, ε, s0):
        """ Fallback initial guess for EE_CRRA_solve when neither an explicit x0 nor a cached
        self.x0['EE_CRRA'] is available: the LOG closed-form (Γs,h,s) path, txE-sliced and stacked into
        the layout EE_CRRA_residual expects -- cheap, and exact at ρ=1. """
        nx = self.T - 1
        sol = self.EE_LOG_solve(τ, θ, ε, s0)
        return np.concatenate([sol['Γs'].values, sol['h'].values[:nx], sol['s'].values[:nx]])

    def EE_CRRA_residual(self, x, τ, θ, ε, s0):
        """ Residual for the CRRA economic equilibrium. Unknowns are the full (Γs_t, h_t, s_t) triple for
        t=0..T-2 (nx=T-1 each, 3*nx total, matching the docs' square system) -- evaluated fully vectorized
        via self.BT with zero Python loops, since each defining equation (auxiliary:h, auxiliary:sFromH,
        auxiliary:Gammas) only checks *consistency* of the candidate. Extraction uses self.ns['EE_CRRA']
        rather than hardcoded slicing. Uses __call__ (raw ndarray), not .get(): this runs inside the
        root-finder's hot loop. Returns the stacked (Γs_resid, h_resid, s_resid), length 3*nx. """
        nx = self.T - 1
        ns = self.ns['EE_CRRA']
        Γs_u, h_u, s_u = ns(x, 'Γs'), ns(x, 'h'), ns(x, 's')
        τ1, θ1 = self.leadSym(τ), self.leadSym(θ)

        s_lag = np.append(s0, s_u)          # s_{t-1}, t=0..T-1
        Γs_pad = np.append(Γs_u, 0)         # dummy -- sFromH's terminal-index output discarded by [:nx]
        s_full = np.append(s_u, 0) 

        h_full = self.BT.FH_h(τ, τ1, θ1, Γs_u, s_lag)
        h_resid = h_u - h_full[:nx]
        s_resid = s_u - self.BT.sFromH(h_full, Γs_pad)[:nx]

        h_lead = np.append(h_full[1:], h_full[-1])                   # h_{t+1}, only [:nx] used
        with np.errstate(divide = 'ignore', invalid = 'ignore'): # s_full[-1]=0 makes B's terminal entry inf/nan; discarded by [:nx] below
            Γs_new = self.BT.Γs(self.BT.B(s_full, h_lead), τ1, θ1)
        Γs_resid = Γs_u - Γs_new[:nx]

        return np.concatenate([Γs_resid, h_resid, s_resid])

    def EE_CRRA_solve(self, τ, θ, ε, s0 = None, x0 = None, update = True, tol = 1e-8, **kwargs):
        """ Root-find the CRRA economic equilibrium. s0: state before db['t'][0] (docs' s_0); defaults to
        the CRRA steady state at db['t'][0]'s (τ,θ) (§4's steadyState_CRRA_solve). x0 defaults to
        self.x0.get('EE_CRRA'), else defaultX0_EE_CRRA. update: cache the solved x for warm-starting later
        calls. Returns {'s','h','Γs','B'} (Γs/B length T-1, h/s length T -- h/s need the FH_h/terminal-
        append reconstruction below since the solver only searches their txE-domain values). """
        if s0 is None:
            s0 = self.steadyState_CRRA_solve(τ[self.B.tFirst], θ[self.B.tFirst], t = self.B.tFirst)['s']
        nx = self.T - 1
        ns = self.ns['EE_CRRA']
        if x0 is None:
            x0 = self.x0['EE_CRRA'] if 'EE_CRRA' in self.x0 else self.defaultX0_EE_CRRA(τ, θ, ε, s0)
        res = optimize.root(self.EE_CRRA_residual, x0, args = (τ, θ, ε, s0), **kwargs)
        self._checkConverged(res.fun, tol = tol, name = 'EE_CRRA_solve', scipyRes = res)
        if update:
            self.x0['EE_CRRA'] = res.x
        Γs_u, h_u, s_u = ns(res.x, 'Γs'), ns(res.x, 'h'), ns(res.x, 's')
        τ1, θ1 = self.leadSym(τ), self.leadSym(θ)
        s_lag = np.append(s0, s_u)
        s = np.append(s_u, 0)
        h = self.BT.FH_h(τ, τ1, θ1, Γs_u, s_lag)
        h_lead = np.append(h[1:], h[-1])
        with np.errstate(divide = 'ignore', invalid = 'ignore'): # s[-1]=0 makes B's terminal entry inf/nan; discarded below
            B_full = self.BT.B(s, h_lead)
        return self._wrapVars({'s': s, 'h': h, 'Γs': Γs_u, 'B': B_full[:-1]})

    def EE_report(self, sol, τ, θ, ε, s0):
        """ Expand a solved {'s','h','Γs','B'} dict (from either solver) into the full set of equilibrium
        objects via base.py, vectorized over the whole time path (self.BT) -- shared by both since
        everything downstream of (s,h,Γs,B) is identical regardless of origin. Γs/B/si_s come back at
        length T-1; everything else (incl. 's_'=s_{t-1} and 'h_', backed out via §5's
        initialState_solve) is length T. """
        s, h, Γs, B = sol['s'].values, sol['h'].values, sol['Γs'].values, sol['B'].values
        s_ = np.append(s0, s[:-1])
        τ1, θ1 = self.leadSym(τ), self.leadSym(θ)
        init = self.initialState_solve(τ[self.B.tFirst], θ[self.B.tFirst])
        h_ = np.append(self.B.hFromS(s0, init['Γs'], self.B.tFirst), h[:-1])

        R, w, w0 = self.BT.R(s_, h), self.BT.w(s_, h), self.BT.w0(s_)
        hi = self.BT.hi(h)
        bbar = self.BT.bbar(τ, w, h, h_)
        bi = self.BT.bi(θ, bbar, h_, τ = τ)      # τ_t reaches b^i through f(θ_t, τ_t) under the size wedge
        b0 = self.BT.b0(ε, bbar, h_)

        # si_s has no terminal-period counterpart at all (s_{T-1}=0 makes the ratio 0/0), so we pad with a
        # nonzero dummy purely to dodge the resulting divide-by-zero warning, then drop that entry (unlike
        # FH_c1i/FH_tildec1i's 0-padding, this dummy carries no economic meaning -- it's discarded, not used).
        B_dummy, Γs_dummy = np.vstack([B, np.ones((1, self.ni))]), np.append(Γs, 1)
        si_s = self.BT.si_s(B_dummy, τ1, θ1, Γs_dummy)[:-1]
        si_s_ = np.vstack([init['si_s'][None, :], si_s])

        c1i = self.BT.FH_c1i(h, s, B, τ1, θ1)
        tc1i = self.BT.FH_tildec1i(h, B, τ1, θ1, Γs)
        c2i = self.BT.c2i(h, s_, τ, θ, si_s_)
        tc2i = self.BT.tildec2i(h, s_, τ, θ, si_s_)
        c10, tc10 = self.BT.c10(s_), self.BT.tildec10(s_)
        c20, tc20 = self.BT.c20(h, s_, ε, τ), self.BT.tildec20(h, s_, ε, τ)

        return self._wrapVars({'s': s, 'h': h, 's_': s_, 'R': R, 'w': w, 'w0': w0, 'hi': hi, 'bbar': bbar, 'bi': bi, 'b0': b0,
                                'c1i': c1i, 'tildec1i': tc1i, 'c2i': c2i, 'tildec2i': tc2i,
                                'c10': c10, 'tildec10': tc10, 'c20': c20, 'tildec20': tc20,
                                'Γs': Γs, 'B': B, 'si_s': si_s})

    #######################################################################
    ##########   4. Steady state solve (docs §2.1, eq:steadystate_LOG, eq:steadystate_CRRA)  ######
    #######################################################################
    # The default s0 source for EE_LOG_solve/EE_CRRA_solve.
    def steadyState_LOG_solve(self, τ, θ, t = None):
        """ Closed-form steady state given LOG preferences (ρ=1, so B^i=β_i is a pure primitive -- no
        root-finding, matching EE_LOG_solve). τ, θ: constant steady-state policy (scalars). """
        t = self.B.tFirst if t is None else t
        B = self.B.get('βi', t)
        Γs = self.B.Γs(B, τ, θ, t)
        return self.steadyState_report(Γs, B, τ, θ, t)

    def steadyState_CRRA_residual(self, Γs, τ, θ, t = None):
        """ Residual for the CRRA steady state (eq:steadystate_CRRA:Gammas): Γs is a steady state iff it
        reproduces itself through self.B.Γs once B^i is made consistent with it via self.B.BSteadyState.
        Scalar in/out, matching scipy.optimize.brentq's signature. """
        return Γs - self.B.Γs(self.B.BSteadyState(Γs, τ, θ, t), τ, θ, t)

    def steadyState_CRRA_bounds(self, τ, θ, t = None, upper = 0.75, margin = 0.99, lower = 1e-6,
                                maxExpand = 40):
        """ Bracket for steadyState_CRRA_solve's Γs search: (lower, min(upper, margin·Base.ΓsCap)),
        widened upward if that does not actually bracket a root.

        A constant `upper` alone (0.75) is NOT safe at US parameters, and it fails in BOTH directions:

        Too high. The Θh denominator turns negative at Base.ΓsCap, which at this calibration falls below
        0.75 once τ ≳ 0.78, and steadyStatePEE_CRRA searches τ over the whole of [l,u]. brentq then
        evaluates a NaN and raises "function value at x=0.75 is NaN" -- a failure that points at the
        solver rather than at the infeasible bracket that caused it. min(upper, margin·ΓsCap) fixes that.

        Too low. At θ = 0 -- the θ = 0 counterfactual, over the whole horizon -- ΓsCap is infinite
        (eq:auxiliary:Thetah's denominator carries θτ), so NOTHING ties the upper end to the model and
        0.75 is a bare guess. At ρ = 2 and τ near 1 the root sits above it: the residual is negative at
        both ends and brentq raises "f(a) and f(b) must have different signs". Hence the expansion, which
        runs only when the default bracket has failed, so it can never change a call that already worked.
        Expanding rather than retuning the constant is the same choice as above -- a bigger constant would
        only move the (ρ, θ, τ) at which it reappears. See notes/crossCuttingFindings.md #7. """
        t = self.B.tFirst if t is None else t
        cap = self.B.ΓsCap(τ, θ, t)
        hi = min(upper, margin*cap)
        f = lambda x: self.steadyState_CRRA_residual(x, τ, θ, t)
        brackets = lambda a, b: np.isfinite(a) and np.isfinite(b) and a*b < 0
        fLo = f(lower)
        if brackets(fLo, f(hi)):
            return (lower, hi)
        for _ in range(maxExpand):
            hi = min(2*hi, margin*cap)
            if brackets(fLo, f(hi)):
                return (lower, hi)
            if hi >= margin*cap:                 # nothing left to expand into
                break
        return (lower, min(upper, margin*cap))   # give up: let brentq report the original failure

    def steadyState_CRRA_solve(self, τ, θ, t = None, bounds = None, tol = 1e-11, **kwargs):
        """ Root-find the CRRA steady state Γs via brentq (bounded scalar search, per the doc's own
        recommendation). bounds default to steadyState_CRRA_bounds -- see there for why a constant upper
        bound is not enough here. """
        t = self.B.tFirst if t is None else t
        bounds = self.steadyState_CRRA_bounds(τ, θ, t) if bounds is None else bounds
        Γs = optimize.brentq(self.steadyState_CRRA_residual, *bounds, args = (τ, θ, t), **kwargs)
        self._checkConverged(self.steadyState_CRRA_residual(Γs, τ, θ, t), tol = tol, name = 'steadyState_CRRA_solve')
        B = self.B.BSteadyState(Γs, τ, θ, t)
        return self.steadyState_report(Γs, B, τ, θ, t)

    def steadyState_report(self, Γs, B, τ, θ, t = None):
        """ Expand a solved (Γs, B) pair into the full steady-state core: Θh/Θs (at constant
        τ_t=τ_{t+1}=τ, θ_{t+1}=θ), then s*/h* -- closed-form regardless of LOG vs. CRRA, since Γs/B
        already encode whichever preference case produced them. """
        t = self.B.tFirst if t is None else t
        Θh = self.B.Θh(τ, τ, θ, Γs, t)
        Θs = self.B.Θs(Θh, Γs, t)
        s = self.B.sSteadyState(Θs, t)
        h = self.B.h(Θh, s, t)
        return {'Γs': Γs, 'B': B, 's': s, 'h': h, 'Θs': Θs}

    #######################################################################
    ##########   5. Initial (pre-determined) state (docs' t=0)        ######
    #######################################################################
    # Γ_{s,-1}, B_{-1}^i and s_{-1,i}/s_{-1} (the generation already old at db['t'][0]) are the CRRA steady
    # state's at db['t'][0]'s own (τ,θ). initialState_solve does NOT take s0 -- only h_{-1} needs it (via
    # Base.hFromS), computed separately in EE_report.
    def initialState_solve(self, τ, θ, t = None, **kwargs):
        """ Identify Γ_{s,-1}, B_{-1}^i, s_{-1,i}/s_{-1} -- see §5 header. τ, θ: scalars, db['t'][0]'s own
        policy. kwargs passed to steadyState_CRRA_solve. Returns {'Γs','B','si_s'} -- no h_{-1} (needs the
        actual s0; see EE_report's Base.hFromS call). """
        t = self.B.tFirst if t is None else t
        ss = self.steadyState_CRRA_solve(τ, θ, t = t, **kwargs)
        Γs, B = ss['Γs'], ss['B']
        si_s = self.B.si_s(B, τ, θ, Γs, t)
        return {'Γs': Γs, 'B': B, 'si_s': si_s}

    #######################################################################
    ##########   6. Politico-economic equilibrium (PEE) solve, LOG    ######
    #######################################################################
    # End-to-end orchestrator: identify τ via policy.py's LOG class, then solve/report the full economic
    # equilibrium it implies. Deliberately thin -- all the actual work is delegated: self.LOG.solveVectorized
    # (no s0 needed, see policy.py) for τ, then EE_LOG_solve/EE_report (§3) for the equilibrium.
    def solvePEE_LOG(self, θ = None, ε = None, s0 = None, solver = 'Robust', **kwargs):
        """ Solve the LOG politico-economic equilibrium end to end. θ, ε: full length-T paths; default to
        db['θ']/db['eps']. s0: if None, the LOG steady state at db['t'][0], evaluated at the *solved*
        τ[db['t'][0]]. Returns {'policy','sol','report'}.

        solver picks which policy.py method identifies τ:
        - 'Robust'     (default) LOG.solveRobust(check = True) -- gradient solve (alg:LOG:fast), then the
                       full-grid search (alg:LOG:gridsearch) on EVERY solve, not as a fallback: the gradient
                       solution is kept if within one grid cell of the selected equilibrium, else the grid
                       selection is polished by a warm-started gradient solve; the candidate counts are
                       returned under 'multiplicity'.
        - 'Vectorized' LOG.solveVectorized alone. Fails loudly rather than falling back.
        - 'Backward'   LOG.solveBackward alone. Capped at grid resolution -- for diagnosing the FOC. """
        if θ is None:
            θ = self.db['θ'].values
        if ε is None:
            ε = self.db['eps'].values
        if solver not in ('Robust', 'Vectorized', 'Backward'):
            raise ValueError(f"solver must be 'Robust', 'Vectorized' or 'Backward', got {solver!r}.")
        policy = getattr(self.LOG, f'solve{solver}')(θ, ε, **kwargs)
        τ = policy['τ'].values
        if s0 is None:
            s0 = self.steadyState_LOG_solve(τ[self.B.tFirst], θ[self.B.tFirst], t = self.B.tFirst)['s']
        sol = self.EE_LOG_solve(τ, θ, ε, s0)
        report = self.EE_report(sol, τ, θ, ε, s0)
        return {'policy': policy, 'τ': policy['τ'], 'sol': sol, 'report': report,
                'multiplicity': policy.get('multiplicity')}   # the tax candidates' counts (num_robustroot.tex)

    #######################################################################
    ##########   7. Politico-economic equilibrium (PEE) solve, CRRA   ######
    #######################################################################
    # solvePEE_LOG's s0 default doesn't survive to CRRA: it reads τ[db['t'][0]] off the *solved* τ path,
    # which only works because LOG's τ never depends on s. CRRA's τ_{tFirst}=τPolicy_{tFirst}(s0) is
    # itself a function of the state being determined, so "solve τ, then read off τ[tFirst]" is circular.
    # steadyStatePEE_CRRA replaces it.
    def steadyStatePEE_CRRA(self, sols, θ, t = None, bounds = None, tol = 1e-10, **kwargs):
        """ solvePEE_CRRA's default s0: the (τ*, s*) fixed point where τ* is both the constant-policy
        steady-state tax (self.steadyState_CRRA_solve) and what the solved policy function would choose
        at that steady state:
        Returns the full steady-state report ({'Γs','B','s','h','Θs'}), not just 's'. """
        t = self.B.tFirst if t is None else t
        θt = θ[t]                                            # scalar -- θ is the full length-T path here
        l, u = self.CRRA.GS['PEE']['gridSettings']['l'], self.CRRA.GS['PEE']['gridSettings']['u']
        bounds = (l, u) if bounds is None else bounds
        τPolicy = sols[t]['τPolicy']
        def residual(τ):
            s = self.steadyState_CRRA_solve(τ, θt, t = t)['s']
            return τ - np.clip(τPolicy(s), l, u)
        τStar = optimize.brentq(residual, *bounds, **kwargs)
        self._checkConverged(residual(τStar), tol = tol, name = 'steadyStatePEE_CRRA')
        return self.steadyState_CRRA_solve(τStar, θt, t = t)

    def solvePEE_CRRA(self, θ = None, ε = None, s0 = None, warmStart = True,
                      backwardKwargs = None, solveKwargs = None):
        """ Solve the CRRA politico-economic equilibrium end to end (docs §PEE/alg:CRRA:grid): identify
        the backward-solved policy functions (CRRA.solveBackward), forward-simulate the tax path they
        imply from an initial state (CRRA.approximatePEE), then solve the equilibrium given that
        path (EE_CRRA_solve/EE_report).

        θ, ε: full length-T paths, as solvePEE_LOG. s0: if None, steadyStatePEE_CRRA's fixed point.
        warmStart: build EE_CRRA_solve's x0 from the forward-simulated (Γs,h,s) path instead of its own
        default (a LOG closed-form proxy, exact only at ρ=1) -- scoped to this call only; passing
        solveKwargs={'x0': ...} overrides it regardless.

        Returns {'sols','τ','sol','report'}: sols is CRRA.solveBackward's {t: report dict} ; τ is
        the path fed to EE_CRRA_solve; sol/report are its usual outputs. ('τ'/'sol'/'report' are named
        to match solvePEE_LOG's return, so §8's calibration can drive either one interchangeably.) """
        if θ is None:
            θ = self.db['θ'].values
        if ε is None:
            ε = self.db['eps'].values
        sols = self.CRRA.solveBackward(θ, ε, **(backwardKwargs or {}))
        if s0 is None:
            s0 = self.steadyStatePEE_CRRA(sols, θ)['s']
        path = self.CRRA.approximatePEE(sols, s0)
        kwargs = dict(solveKwargs or {})
        if warmStart and 'x0' not in kwargs:
            kwargs['x0'] = np.concatenate([path['Γs'], path['h'], path['s']])
        sol = self.EE_CRRA_solve(path['τ'].values, θ, ε, s0, **kwargs)
        report = self.EE_report(sol, path['τ'].values, θ, ε, s0)
        return {'sols': sols, 'τ': path['τ'], 'sol': sol, 'report': report,
                'multiplicity': getattr(self.CRRA, 'lastMultiplicity', None)}

    #######################################################################
    ##########   8. Calibration (docs §calibration, eq:calibration)   ######
    #######################################################################
    # Nested fixed point: outer root over (β,ω) against eq:calibration's two targets (R_{t0}, τ_{t0}); each
    # residual evaluation refreshes the auxiliary parameters and re-solves the whole PEE (inner loop).
    # ηi/Xi/θ are NOT here -- fixed once in §2, independent of (β,ω).
    #
    # Two targets, R_{t0} and τ_{t0}. The commonX variant adds a third target (average hours) but NOT a
    # third dimension: X enters no aggregate, so it is recovered in closed form after the root, never
    # inside it. See calibrate.

    _calPars = ('β', 'ω')
    # Unbounded reparameterization below: finite upper bound -> logit, infinite -> log. Only positivity is
    # enforced -- β is NOT capped at 1 (simpleβ sets βj=β·p_j, the actual discount factor), so that the
    # search fails rather than stalling on a bound if a target needs β>1.
    _calBounds = {k: (0., np.inf) for k in _calPars}

    @property
    def calibrationPars(self):
        """ (β,ω) as db currently holds them. Inverse of _calSetPars. """
        t0 = self.db['t'][self.db['t0']]
        return {'β': self.simpleβinv(), 'ω': self.db['ω'].xs(t0)}

    def _calToX(self, pars):
        """ Calibration-parameter dict -> unbounded search vector. """
        def f(k):
            l, u = self._calBounds[k]
            v = pars[k]
            return np.log(v-l) if np.isinf(u) else np.log((v-l)/(u-v))
        return np.array([f(k) for k in self._calPars])

    def _calFromX(self, x):
        """ Inverse of _calToX. """
        def f(k, xi):
            l, u = self._calBounds[k]
            return l+np.exp(xi) if np.isinf(u) else l+(u-l)/(1+np.exp(-xi))
        return {k: f(k, xi) for k, xi in zip(self._calPars, x)}

    def _calSetPars(self, pars):
        """ Write (β,ω) into db and refresh dependent auxiliary parameters (κ). """
        self.db.update(self.adjPar('β', pars['β']))
        self.db.update(self.adjPar('ω', pars['ω']))
        self.updateAuxPars()

    def calibration_report(self, pars, preferences, solveKwargs = None):
        """ Install `pars`, solve the PEE, evaluate the target quantities at db['t0'].

        Returns the two targeted quantities (R, τ), plus the savings rate and average hours, which are
        reported but not targeted -- 'hbar' because under commonX it is what solveCommonX is about to be
        pointed at.

        solveKwargs must not pin θ/ε: they need to stay at their None defaults so solvePEE_* re-reads
        db['θ']/db['eps'] after _calSetPars refreshes them. """
        self._calSetPars(pars)
        out = getattr(self, f'solvePEE_{preferences}')(**(solveKwargs or {}))
        t0 = self.db['t'][self.db['t0']]
        rep, τ = out['report'], out['τ'].xs(t0)
        h, s, s_, R = rep['h'].xs(t0), rep['s'].xs(t0), rep['s_'].xs(t0), rep['R'].xs(t0)
        return {'R': R, 'τ': τ, 'h': h, 'hbar': self.B.avgHours(h, t0),
                'sr': self.B.savingsRate(s, s_, h, t0), 'Θh': self.B.ΘhFromH(h, s_, t0), 'PEE': out}

    def _calResidual(self, report):
        """ eq:calibration's targets as a residual on an already-evaluated calibration_report. One entry
        per element of _calPars, same order.

        Split out from calibration_residual because `calibrate` re-forms the residual at the converged
        point, so a subclass that changes the targets (modelFR.ModelFR drops R) must be able to change
        them in exactly one place. R enters relatively (R≈2.4) and τ as a level gap (τ≈0.14), so the two
        carry the same O(0.1) magnitude and the root finder cannot trade accuracy in one for the other. """
        return np.array([report['R']/self.db['R0'] - 1, report['τ'] - self.db['τ0']])

    def calibration_residual(self, x, preferences, solveKwargs = None):
        """ eq:calibration as a residual on the unbounded vector. """
        return self._calResidual(self.calibration_report(self._calFromX(x), preferences, solveKwargs))

    def calibrate(self, preferences = None, x0 = None, tol = 1e-8, update = True,
                  solveKwargs = None, verifyCommonX = True, **kwargs):
        """ Solve eq:calibration for (β,ω): one 2-D root, each evaluation a full PEE solve.
        preferences defaults to 'LOG' iff ρ=1 (CRRA's recursion divides by 1-1/ρ). x0 defaults to
        self.x0['calibration'] if cached, else db's current parameters. kwargs -> scipy.optimize.root.

        Under commonX, X is calibrated to db['h0'] AFTER the root, in closed form (solveCommonX), and the
        equilibrium is then re-solved so that db and the returned report hold the final X. This is
        block-recursive rather than a third search dimension because X enters no aggregate: the re-solve
        reproduces the same (R,τ) it was handed. verifyCommonX checks exactly that, and it is checking an
        identity, not an approximation -- any drift beyond `tol` means an η/X has leaked into an aggregate
        through something other than η^(1+ξ)/X^ξ, so it raises rather than warns.

        db is mutated as the search proceeds; restored to its entry state on any failure (a shallow copy
        suffices -- adjPar/updateAuxPars always rebind db keys rather than mutate in place). On success db
        is left holding the converged parameters and equilibrium. """
        if preferences is None:
            preferences = 'LOG' if self.db['ρ'].xs(self.db['t'][self.db['t0']]) == 1 else 'CRRA'
        if x0 is None:
            x0 = self.x0.get('calibration', self._calToX(self.calibrationPars))
        if self.commonX and 'h0' not in self.db:
            raise KeyError("calibrate: commonX=True needs the average-hours target db['h0'] (see test.py).")
        snapshot = dict(self.db)
        try:
            res = optimize.root(self.calibration_residual, x0, args = (preferences, solveKwargs), **kwargs)
            pars = self._calFromX(res.x)
            report = self.calibration_report(pars, preferences, solveKwargs)
            residual = self._calResidual(report)
            self._checkConverged(residual, tol = tol, name = 'calibrate', scipyRes = res)
            pars, report = self._calPostRoot(pars, report, preferences, solveKwargs, tol, verifyCommonX)
        except Exception:
            self.db.clear() # self.B/BG/BT hold a reference to this dict, so restore in place
            self.db.update(snapshot)
            raise
        if update:
            self.x0['calibration'] = res.x
        return {'pars': pars, 'x': res.x, 'residual': residual, 'report': report, 'scipyRes': res}

    def _calPostRoot(self, pars, report, preferences, solveKwargs, tol, verify):
        """ The block-recursive half of eq:calibration: parameters that hit a target WITHOUT entering the
        root, because they enter no aggregate the residual can see. Returns (pars, report) with the final
        parameters installed in db and the equilibrium re-solved so db and the report agree with them.

        Here that is commonX's X (docs eq:calibration:Xsolve). `verify` re-checks that the re-solve
        reproduced the targeted quantities; it is checking an IDENTITY, not an approximation -- any drift
        beyond tol means an η/X has leaked into an aggregate through something other than η^(1+ξ)/X^ξ --
        so it raises rather than warns. modelFR.ModelFR extends this with the vector-X hours rescaling. """
        if not self.commonX:
            return pars, report
        X = self.solveCommonX(report['h'])
        self.initProductivity_commonX(X = X)
        self.updateAuxPars() # Γh/θ are functions of ηi/Xi -- θ via getθ, so it must be refreshed here
        pars = pars | {'X': X}
        reportX = self.calibration_report({k: pars[k] for k in self._calPars}, preferences, solveKwargs)
        if verify:
            drift = np.array([reportX['R']/report['R'] - 1, reportX['τ'] - report['τ'],
                              reportX['hbar']/self.hbarTarget() - 1])
            self._checkConverged(drift, tol = max(tol, 1e-8), name = 'calibrate (commonX rescaling)')
        return pars, reportX

    #######################################################################
    ##########   8b. Calibration across a parameter grid              ######
    #######################################################################
    # A march over ρ, anchored at ρ=1: the one value where the LOG solver applies, so it is both the
    # cheapest point on the grid and the only one needing no warm start. Every other point is seeded from
    # its neighbours in the UNBOUNDED coordinate, where extrapolation cannot step β or ω across zero.

    def _calPreferences(self):
        """ Which solver eq:calibration's inner loop uses at db's current ρ. Mirrors calibrate's own
        default -- CRRA's recursion divides by 1-1/ρ and its guard refuses ρ=1 -- and is needed separately
        because the *grid settings* have to be installed on that solver before calibrate is called. """
        return 'LOG' if self.db['ρ'].xs(self.db['t'][self.db['t0']]) == 1 else 'CRRA'

    @staticmethod
    def _calGridSettings(gridSettings, preferences):
        """ Resolve calibratePoint's gridSettings for the solver actually in use. A flat dict applies to
        whichever solver runs; a dict keyed by 'LOG'/'CRRA' selects per solver, which is what a march
        across ρ generally wants -- 'ns' sizes a state grid the LOG solver does not have. """
        if gridSettings and set(gridSettings) & {'LOG', 'CRRA'}:
            return gridSettings.get(preferences)
        return gridSettings

    def _calVerify(self, x, preferences, policy, verify):
        """ Re-evaluate the outer residual at the converged x on a finer inner grid (the installed grid
        settings updated by `verify`), and restore the installed settings afterwards. Returns
        max|residual|, or NaN if the refined solve fails.

        This is the one diagnostic that separates a converged calibration from a *resolved* one: the outer
        root can settle on a point that is displaced rather than merely imprecise, and the symptom is that
        refining the inner grid at fixed parameters leaves the residual on a plateau instead of decaying.
        Read the LEVEL across points, not the shape within one -- the rung at the calibration's own grid
        is ~0 by construction. Recorded per point, never asserted: no grid size has standing at every ρ. """
        working = dict(policy.GS['PEE']['gridSettings'])
        try:
            policy.initGS(working | verify)
            return float(np.max(np.abs(self.calibration_residual(x, preferences))))
        except Exception:
            return np.nan
        finally:
            policy.initGS(working)

    def calibratePoint(self, value, x0 = None, par = 'ρ', gridSettings = None, verify = None,
                       calKwargs = None):
        """ Set db[par] = value, then calibrate there. Returns one flat record: the calibrated parameters,
        the unbounded x (which seeds the next point), max|residual|, the target quantities, and the
        diagnostics a sweep has to be readable by afterwards (wall time, scipy's nfev, the inner grid
        actually used, and 'multiplicity': the tax candidates' counts of the PEE solve at the calibrated
        parameters, policy.multiplicitySummary, num_robustroot.tex eq:equilibriumTest).

        x0 is the *unbounded* vector, not a parameter dict. gridSettings (flat, or keyed by 'LOG'/'CRRA',
        _calGridSettings) is installed on the solver in use; None keeps the settings already installed.
        db is left holding the converged point, which is deliberate: it is the next point's inner warm
        start. """
        self.db.update(self.adjPar(par, value))
        preferences = self._calPreferences()
        policy = getattr(self, preferences)
        settings = self._calGridSettings(gridSettings, preferences)
        if settings is not None:
            policy.initGS(settings)
        tStart = time.time()
        cal = self.calibrate(preferences = preferences, x0 = x0, **(calKwargs or {}))
        elapsed = time.time()-tStart
        rep = cal['report']
        rec = {par: float(value), 'preferences': preferences, 'x': cal['x'], 'time': elapsed,
               'residual': float(np.max(np.abs(cal['residual']))),
               'R': float(rep['R']), 'τ': float(rep['τ']), 'sr': float(rep['sr']),
               'h': float(rep['h']), 'hbar': float(rep['hbar']),
               'nfev': int(getattr(cal['scipyRes'], 'nfev', -1)),
               'gridSettings': dict(policy.GS['PEE']['gridSettings']),
               'multiplicity': rep['PEE'].get('multiplicity')}
        rec.update({k: float(v) for k, v in cal['pars'].items()})
        verifySettings = self._calGridSettings(verify, preferences) if verify else None
        if verifySettings:
            rec['verifyResidual'] = self._calVerify(cal['x'], preferences, policy, verifySettings)
        return rec

    def calibrateGrid(self, grid, par = 'ρ', anchor = 1.0, gridSettings = None, verify = None,
                      x0 = None, degree = 1, maxHalvings = 2, onPoint = None, calKwargs = None,
                      **kwargs):
        """ Calibrate at every value of `grid`, marching outward from `anchor` (gridsearch.continuation's
        marchGrid -- see it for the extrapolation and failure-recovery contract; kwargs pass through).

        The march runs on this one instance, so each point also inherits db's converged parameters as its
        inner warm start; that is separate from, and additional to, the extrapolated x0.

        onPoint(record) fires after every point including failures -- pass a writer here, since a sweep is
        one crash away from losing everything solved so far. """
        def solve(value, xStart):
            return self.calibratePoint(value, x0 = xStart, par = par, gridSettings = gridSettings,
                                       verify = verify, calKwargs = calKwargs)
        return continuation.marchGrid(grid, solve, x0 = x0, anchor = anchor, degree = degree,
                                      maxHalvings = maxHalvings, onPoint = onPoint, **kwargs)

    #######################################################################
    ##########   9. Model copies for shock experiments (docs' t0)     ######
    #######################################################################
    # createCopyFromt0 produces a new, independent model instance whose horizon is restricted to
    # db['t'] >= t0 -- for "unexpected shock" experiments: solve the baseline PEE over the full horizon,
    # build the copy, then re-solve PEE on it with s0 (stateAtT0) read off the baseline's own report at
    # t0. The copy's own db['t'] is renumbered to start at 0 again (see _sliceDb's docstring for why --
    # model.py's own EE_LOG_solve/EE_CRRA_solve/EE_report/initialState_solve index a caller's τ/θ/ε
    # ndarrays positionally via self.B.tFirst, which only works when tFirst is 0), so self.B.tFirst on
    # the copy is 0, not t0 -- the original calendar year t0 is not retained anywhere on the copy.
    # State seeding is deliberately NOT done here -- s0 stays an explicit argument to
    # solvePEE_LOG/solvePEE_CRRA on the copy, exactly as for any other instance; see stateAtT0.
    def createCopyFromt0(self, t0):
        """ Return an independent copy of this model with its time horizon restricted to db['t'] >= t0
        and renumbered to start at 0 (db['t'] runs 0..T-t0-1 on the copy, not t0..T-1).
        Warm-start caches (self.x0, self.LOG.x0, self.CRRA.x0) are cleared (stale/wrong length for the
        new horizon); self.LOG.GS/self.CRRA.GS (state-space grids, not time-indexed) are left as-is.
        db['t0'] (the calibration-baseline-year *position* -- unrelated to this t0; see
        default0DParams' docstring) is shifted by -t0 if the calibration year still falls inside the new
        horizon, else set to None: a stale position would otherwise silently resolve to the wrong year
        rather than failing loudly. Recalibrating a copy needs a caller-supplied db['t0']. """
        if t0 not in self.db['t']:
            raise ValueError(f"t0={t0!r} is not in db['t'] (={list(self.db['t'])}).")
        mt0 = deepcopy(self)
        _sliceDb(mt0.db, t0)
        mt0.T = len(mt0.db['t'])
        for baseIns in (mt0.B, mt0.BG, mt0.BT):
            baseIns.tFirst = mt0.db['t'][0]
        mt0.LOG.T = mt0.CRRA.T = mt0.T
        mt0.x0, mt0.LOG.x0, mt0.CRRA.x0 = {}, {}, {}
        mt0.addNamespaces()
        mt0.db['t0'] = (self.db['t0'] - t0) if (self.db['t0'] is not None and self.db['t0'] >= t0) else None
        return mt0

    def stateAtT0(self, report, t0):
        """ The state entering t0 in an already-solved report (solvePEE_LOG/solvePEE_CRRA's 'report'),
        for seeding createCopyFromt0(t0)'s own solvePEE_LOG(s0=...)/solvePEE_CRRA(s0=...). """
        return {'s0': report['s_'].xs(t0)}

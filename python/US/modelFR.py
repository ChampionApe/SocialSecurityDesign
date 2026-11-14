r""" ModelFR -- the US model calibrated against a US reference rather than on its own.

France and the UK are the same model as `ModelUS` (see model.py and writing/US/); only the calibration
protocol differs, in the two ways the paper's calibration of the rich OECD countries states:

  1. beta is IMPOSED at the value the US calibration produced at the same rho, not searched. The interest
     rate stops being a target with it, so eq:calibration collapses from a 2-D root over (beta, omega) to
     a 1-D root over omega against tau_{t0}. R is still reported -- it is now a prediction.

  2. Average hours are targeted in BOTH calibration variants, not only under commonX, and the target is
     defined relative to the US:

         hbar_FR = hbar_US * workweek_FR / workweek_US

     with hbar_US the calibrated US model's average hours at the same rho and the workweeks from data.
     Under vector X this is the only form of hours target that identifies anything: nothing else in that
     calibration pins the level of hbar (docs eq:us:model:hoursUnit; see the README's note that hbar must
     be reported against a reference point). Under commonX it is not a new target at all:
     hbar_US = h0_US exactly there, so it reduces to the observed workweek_FR the US variant already uses.

Under vector X the target is hit by moving every X_i by the same proportion -- the scale invariance
(eq:us:model:scaleInvariance), applied through ModelUS.rescaleX. Consequences worth stating plainly:

  * Gamma_h = 1 NO LONGER HOLDS on a calibrated ModelFR. initProductivity_vectorX still imposes it, but
    the post-root rescaling leaves Gamma_h = lambda. Note what this is and is not: Gamma_h = 1 is a
    NORMALISATION in ModelUS, where nothing pins the level of hours; here the hours target pins it, so
    Gamma_h comes out DATA-DETERMINED rather than free. It still sets the level of h, s, c and Y (it is
    unidentified only jointly with TFP, which base.py normalises to A = 1), so it is not a quantity to
    treat as irrelevant -- only one that no longer equals 1. A check asserting Gamma_h = 1 is wrong here.
  * lambda is block-recursive to the omega root, exactly as commonX's X is: the scale invariance leaves
    tau, R, w, theta, s_i/s and the savings rate pointwise unchanged, so the rescaling cannot disturb the
    target omega was just solved for. That is asserted, not assumed, by the `verify` drift check.

The design needs no code of its own: France's θ0 is the datum of the paper's θ variant ('C' of
data/oecd26_theta.csv, 0.822; the UK's 0.181), imposed as db['θ0'] (oecd26.countryθ0) and returned by
ModelUS.getθ for any income grouping, so calibrate re-derives nothing from RR0 and initProductivity/updateAuxPars
re-install θ0 (#9). So this class serves the UK unchanged; only the workbook differs.

Needs, beyond a ModelUS's inputs: db['h0'] (the observed workweek in the model's time unit, as in
test.py) and a US reference -- see setUSRef and usRefFromCsv.
"""
import numpy as np
from model import ModelUS


class USReference:
    """ Callable rho -> {'β', 'hbar', 'h0'}, backed by a table of solved US points.

    A plain class rather than a closure so that it PICKLES: calibrateRhoGridEU.py stores one ModelFR
    instance per swept point, and the instance holds its reference source.

    Lookup is exact -- keyed on rho rounded to 6 dp, the same convention the sweep drivers use for their
    grid keys -- and deliberately never interpolated. A beta read off a neighbouring rho is not the US
    calibration's beta at this rho, which would make "impose the US discount factor" quietly false. """

    def __init__(self, table, h0):
        self.table = {round(float(k), 6): {'β': float(v['β']), 'hbar': float(v['hbar'])}
                      for k, v in dict(table).items()}
        self.h0 = float(h0)

    def __call__(self, ρ):
        key = round(float(ρ), 6)
        if key not in self.table:
            raise KeyError(f"USReference: no US calibration at ρ={key} (have {sorted(self.table)}). "
                           "Run python/US/calibrateRhoGrid.py over that value first.")
        return self.table[key] | {'h0': self.h0}

    def __repr__(self):
        return f'USReference({len(self.table)} points, h0={self.h0:.6f})'


def usRefFromCsv(path, h0US, rho = None):
    """ A USReference built from a rho-sweep csv written by calibrateRhoGrid.py
    (results/calibration/US_rhoGrid{,CommonX}.csv). Returns the reference, or the resolved dict if `rho`
    is given.

    h0US is the US model's own db['h0'] (workweek_US/(7*12)) -- it is not in the csv, being an input
    rather than a result. Pass the csv of the SAME variant as the model being calibrated: hbar differs
    between vector X and commonX by construction, and it is hbar that carries the reference. Rows whose
    calibration failed (NaN residual) are skipped rather than loaded as if they were answers. """
    import pandas as pd
    df = pd.read_csv(path)
    df = df[np.isfinite(pd.to_numeric(df['residual'], errors = 'coerce'))]
    ref = USReference({r['ρ']: r for _, r in df.iterrows()}, h0US)
    return ref if rho is None else ref(rho)


class ModelFR(ModelUS):
    """ ModelUS with the US-referenced calibration described in this module's docstring. Everything
    outside eq:calibration -- the equilibrium, both PEE solvers, the shock machinery -- is inherited
    unchanged, because none of it differs. """

    # omega alone. beta is installed by setUSRef before the search and never enters it; _calBounds must
    # be redefined with it, since ModelUS builds its own from ModelUS._calPars at class-body time.
    _calPars = ('ω',)
    _calBounds = {k: (0., np.inf) for k in _calPars}

    # Tolerance for the post-rescaling drift check, per solver. The scale invariance is EXACT, so under
    # LOG the drift is machine precision and the check is asserted at 1e-8.
    #
    # Under CRRA policy.py's defaultSGrid floors the state grid at an ABSOLUTE 1e-4 (s_{T-1}=0 makes
    # Rlead/Bi/si_s undefined and every model here can reach exact zero; the floor must NOT be made to
    # scale with the model), while its upper bound 1.25*sMax moves with the model. Rescaling X therefore
    # shifts the floor relative to the grid and τ moves by O(1e-4), a plateau in ns rather than a
    # resolution error, and small only while 1e-4 is a small share of sMax: re-check at a calibration whose
    # sMax is an order of magnitude below the US's. 1e-3 passes that while still catching what the check
    # is for -- an eta/X leaking into an aggregate through something other than y^eta, an O(1) error. The
    # measured drift is recorded on the report as 'hoursDrift'; read it, do not assume it.
    hoursDriftTol = {'LOG': 1e-8, 'CRRA': 1e-3}

    def __init__(self, *args, usRef = None, **kwargs):
        """ usRef: either a dict {'β', 'hbar', 'h0'} (a fixed reference, for a single-rho calibration) or
        a callable rho -> that dict (for a march over rho, where beta and hbar both move with rho -- see
        usRefFromCsv). Installed lazily, at calibrate time, so db['rho'] is whatever the caller has by
        then set it to. """
        super().__init__(*args, **kwargs)
        self.usRefSource = usRef
        self.usRef = None

    ##########  The US reference  ##########

    def setUSRef(self, **override):
        """ Resolve the US reference at db's current rho and install beta from it. Called by calibrate;
        call it directly only to solve (rather than calibrate) at the imposed beta.

        override replaces individual entries of the resolved reference -- the escape hatch for a
        one-off, and what the tests use to drive the workweek ratio. """
        src = self.usRefSource
        if src is None and not override:
            raise ValueError("ModelFR: no US reference. Pass usRef= to __init__ (a dict {'β','hbar','h0'} "
                             "or a callable of ρ -- see usRefFromCsv), or call setUSRef(**ref).")
        ρ = self.db['ρ'].xs(self.db['t'][self.db['t0']])
        # dict(...) on BOTH branches: a callable source is free to hand back a cached/shared dict, and
        # applying `override` in place would then edit the caller's reference table for every later ρ.
        ref = dict(src(ρ) if callable(src) else src) if src is not None else {}
        ref.update(override)
        missing = {'β', 'hbar', 'h0'} - set(ref)
        if missing:
            raise KeyError(f"ModelFR: US reference is missing {sorted(missing)}; got {sorted(ref)}.")
        self.usRef = ref
        self.db.update(self.adjPar('β', ref['β']))
        self.updateAuxPars()
        return ref

    def hbarTarget(self):
        """ hbar_US * workweek_FR/workweek_US, with the workweek ratio read off db['h0']/usRef['h0']
        (both are workweek/(7*12), so the units cancel and only the ratio survives -- which is the whole
        point: under vector X only a RATIO of average hours is meaningful). """
        if self.usRef is None:
            raise ValueError("ModelFR.hbarTarget: call setUSRef() first (calibrate does it for you).")
        return self.usRef['hbar'] * self.db['h0']/self.usRef['h0']

    ##########  eq:calibration  ##########

    def _calSetPars(self, pars):
        """ omega only. beta is not a calibration parameter here and must NOT be written from `pars`:
        it is installed once by setUSRef, and a stray write is exactly how an imposed value silently
        becomes a searched one. """
        self.db.update(self.adjPar('ω', pars['ω']))
        self.updateAuxPars()

    def _calResidual(self, report):
        """ One target: tau_{t0}. R is reported but free (see the module docstring). """
        return np.array([report['τ'] - self.db['τ0']])

    def calibrate(self, *args, **kwargs):
        """ ModelUS.calibrate, with the US reference resolved at db's current rho first (which is what
        installs beta) and the average-hours target required in both variants, not just under commonX.

        Productivity is reset to the Gamma_h = 1 baseline before the root. That is what makes lambda a
        LEVEL rather than an increment: rescaleX multiplies the X_i already in db, so without the reset a
        second calibrate -- or the next point of a march, which reuses one instance -- would start from
        the previous solution's X and report the small step from there instead of the total rescaling.
        The equilibrium is identical either way (the target is absolute, so the rescaling lands on it from
        wherever it starts), so this only ever corrupted the reported lambda; with the reset,
        lambda == Gamma_h exactly, which is the identified quantity.

        initProductivity is deterministic in the workbook data alone -- it reads zx/zeta, not beta, omega
        or rho -- so the reset costs no warm start. Note it happens BEFORE super()'s db snapshot, so a
        failed calibrate restores db to the baseline productivity rather than to whatever X it was handed;
        that state is re-derivable from data, so nothing is lost. """
        self.setUSRef()
        if 'h0' not in self.db:
            raise KeyError("ModelFR.calibrate: needs the average-hours target db['h0'] in BOTH variants "
                           "(vector X hits it by rescaling X_i; see the module docstring).")
        self.initProductivity()
        self.updateAuxPars()   # Γh back to 1, and θ is a function of η/X
        return super().calibrate(*args, **kwargs)

    def _calPostRoot(self, pars, report, preferences, solveKwargs, tol, verify):
        """ ModelUS's commonX step, plus -- under vector X -- the proportional X_i rescaling that puts
        hbar_{t0} on hbarTarget().

        lambda is closed-form and exact, not iterated: the scale invariance scales hbar by exactly
        lambda, so lambda = target/hbar of the solved equilibrium. It is applied AFTER the root for the
        same reason commonX's X is -- it moves no quantity the residual reads. The drift check asserts
        that: tau and R must come back pointwise unchanged and hbar must land exactly on the target.

        beta is folded into `pars` here so that a sweep record carries the value that was imposed --
        it is absent from _calPars and would otherwise never be written down. """
        pars, report = super()._calPostRoot(pars, report, preferences, solveKwargs, tol, verify)
        pars = pars | {'β': self.simpleβinv()}
        if self.commonX:
            self._hoursDrift = np.nan   # X carries the hours unit there; there is no rescaling to drift
            return pars, report
        λ = self.hbarTarget()/report['hbar']
        self.rescaleX(λ)
        pars = pars | {'λ': λ}
        reportλ = self.calibration_report({k: pars[k] for k in self._calPars}, preferences, solveKwargs)
        drift = np.array([reportλ['R']/report['R'] - 1, reportλ['τ'] - report['τ'],
                          reportλ['hbar']/self.hbarTarget() - 1])
        reportλ['hoursDrift'] = self._hoursDrift = float(np.max(np.abs(drift)))
        if verify:
            self._checkConverged(drift, name = 'calibrate (hours-unit rescaling)',
                                 tol = max(tol, self.hoursDriftTol[preferences]))
        return pars, reportλ

    ##########  Sweeps  ##########

    def calibratePoint(self, value, *args, par = 'ρ', **kwargs):
        """ As ModelUS's, except that a march over rho also moves the reference: beta and hbar_US are
        both functions of rho, so usRefSource must be a callable when par == 'ρ'. Resolved here (not
        only inside calibrate) so that setUSRef fires against the NEW rho, which db[par] = value has
        just installed. """
        self.db.update(self.adjPar(par, value))
        if par == 'ρ' and not callable(self.usRefSource):
            raise TypeError("ModelFR.calibratePoint: marching over ρ needs a usRef callable of ρ -- β and "
                            "h̄_US both move with ρ, and a fixed reference would import one ρ's β at every "
                            "point. See usRefFromCsv.")
        self.setUSRef()
        rec = super().calibratePoint(value, *args, par = par, **kwargs)
        # Carried into the sweep csv rather than left in the report: under CRRA it is the measured cost of
        # defaultSGrid's absolute floor (see hoursDriftTol), and a level that moves across a march is the
        # thing worth seeing.
        rec['hoursDrift'] = float(getattr(self, '_hoursDrift', np.nan))
        return rec

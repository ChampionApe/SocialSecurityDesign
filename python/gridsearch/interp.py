""" Interpolation, smoothing and differentiation of gridded data.

Three jobs, all model-agnostic like the rest of this package:

  griddedInterp1D   -- a callable interpolant, e.g. a solved policy tau(s)/h(s) to be evaluated later at
                       an arbitrary candidate state, not only at the nodes the grid search visited.
  griddedInterp2D   -- the same for a policy over two states, evaluated on paired coordinates.
  griddedSmooth1D   -- denoise a solved profile before interpolating it.
  griddedGradient1D -- a derivative along a grid axis, where no closed form exists.

The latter two are built on a *smoothing* spline rather than finite differences, deliberately and
following the precedent this replaces: their inputs are themselves outputs of interpolated policy
functions, so they carry small kinks that plain differencing would amplify into the derivative. The
smoothing factor is what controls that, and it is exposed rather than fixed.
"""
import numpy as np
from scipy import interpolate, ndimage


def _validityMask(values, ndim):
    """ Per-node validity for a values array whose first `ndim` axes are grid axes. A node counts as valid
    only if every trailing entry at it is finite, so a partially-NaN node cannot be half-used. """
    finite = np.isfinite(np.asarray(values, dtype = float))
    return finite if finite.ndim == ndim else finite.reshape(finite.shape[:ndim] + (-1,)).all(axis = -1)


def _fillInvalid(values, valid):
    """ Replace invalid nodes with their nearest valid neighbour, so a spline can be *constructed* over a
    grid carrying infeasible cells. This is scaffolding for the fit only -- the caller must still mask the
    result back to NaN through _maskedCall, or the filled values would be reported as if they were policy.
    The fill still shapes the spline in the valid cells near an invalid node ('cubic' is a global spline,
    'pchip' reads the neighbouring nodes' slopes), so there the result is not a fit over the valid nodes
    alone. An all-invalid grid is returned as zeros; the mask then makes every evaluation NaN anyway. """
    out = np.array(values, dtype = float)
    if valid.all():
        return out
    if not valid.any():
        return np.zeros_like(out)
    idx = ndimage.distance_transform_edt(~valid, return_distances = False, return_indices = True)
    return out[tuple(idx)]


def _maskedCall(f, fNaN):
    """ Wrap an interpolant so that it returns NaN exactly where the 'linear' kind over the unfilled values
    would. fNaN is that linear interpolant of an array that is 0 at valid nodes and NaN at invalid ones, so
    the same interval lookup decides: a point in a cell with an invalid corner, a point extrapolated from
    such an edge cell, and a point exactly on a valid node whose lookup interval has an invalid corner
    (interp1d assigns a node to the interval on its left, RegularGridInterpolator to the one on its right)
    all return NaN. An infeasible region has to surface as a non-finite path rather than as a plausible
    interpolated policy. """
    def call(*args):
        out = np.asarray(f(*args), dtype = float)
        bad = np.isnan(np.asarray(fNaN(*args), dtype = float))
        return np.where(bad.reshape(bad.shape + (1,)*(out.ndim-bad.ndim)), np.nan, out)
    return call


def griddedInterp1D(x, y, kind = 'linear'):
    """ x: (M,) grid nodes -- sorted here defensively, since only the axis a search runs over (roots1d's
    own precondition) is guaranteed sorted, not every state axis a grid happens to carry. y: (M,) or
    (M, ...) values at those nodes; trailing axes (e.g. a per-household-type quantity) pass through
    unchanged.

    Returns a callable f(xNew) -> yNew that extrapolates past [min(x), max(x)] rather than raising or
    clamping to the boundary value: a candidate state produced by an outer root-solve (e.g. an
    endogenous continuation state implied by a candidate policy) need not land inside the grid this
    interpolant was built from, and clamping would misrepresent the policy at exactly the states where the
    boundary is informative.

    kind: 'linear' (piecewise linear, NaN wherever the lookup interval has an invalid node), any other
    interp1d kind ('cubic': a not-a-knot cubic spline), or 'pchip' (monotone C1 cubic, no overshoot where
    a policy is flat at a bound). The non-linear kinds cannot be constructed over NaN, so invalid nodes are
    filled for the fit and masked on evaluation (_maskedCall): they return NaN exactly where 'linear'
    does; elsewhere near an invalid node their values depend on the fill (_fillInvalid). """
    order = np.argsort(x)
    xs, ys = np.asarray(x)[order], np.asarray(y, dtype = float)[order]
    if kind == 'linear':
        return interpolate.interp1d(xs, ys, kind = kind, axis = 0,
                                    bounds_error = False, fill_value = 'extrapolate')
    valid = _validityMask(ys, 1)
    filled = _fillInvalid(ys, valid)
    if kind == 'pchip':
        spline = interpolate.PchipInterpolator(xs, filled, axis = 0, extrapolate = True)
    else:
        spline = interpolate.interp1d(xs, filled, kind = kind, axis = 0,
                                      bounds_error = False, fill_value = 'extrapolate')

    def f(xNew):
        """ As in the 2-D case, a non-finite coordinate must return NaN rather than reach the spline. """
        a = np.asarray(xNew, dtype = float)
        fin = np.isfinite(a)
        if fin.all():
            return spline(a)
        probe = spline(np.where(fin, a, xs[0]))
        return np.where(fin.reshape(fin.shape + (1,)*(probe.ndim-a.ndim)), probe, np.nan)

    if valid.all():
        return f
    fNaN = interpolate.interp1d(xs, np.where(valid, 0., np.nan), kind = 'linear',
                                bounds_error = False, fill_value = 'extrapolate')
    return _maskedCall(f, fNaN)


def griddedInterp2D(x, y, z, kind = 'linear'):
    """ The two-dimensional counterpart of griddedInterp1D, for a policy that is a function of two states.

    x: (M,) nodes of axis 0; y: (N,) nodes of axis 1; z: (M, N) or (M, N, ...) values, with trailing axes
    passed through unchanged. Both axes are sorted here defensively (z with them), since a state axis is
    not guaranteed sorted by any caller.

    Returns a callable f(xNew, yNew) -> values, evaluated *elementwise* on broadcast-compatible arrays --
    not on their outer product -- because the caller's two coordinates are paired (one candidate state per
    point), never a grid to be crossed.

    Extrapolates past the rectangle rather than clamping, for the same reason as griddedInterp1D: a
    candidate state implied by an outer root-solve need not land inside the grid the interpolant was built
    from, and clamping would silently return the boundary policy for an interior question. How far that
    extrapolation may be trusted is a question for the caller's own reachability check, not for this
    function.

    kind: 'linear' or any RegularGridInterpolator method -- 'pchip' for a monotone C1 surface, 'cubic'
    for a smoother but overshooting one (a tensor-product spline, global along each axis). A policy
    surface carries NaN at its infeasible nodes and the spline methods refuse to be constructed over them,
    so invalid nodes are filled for the fit and masked on evaluation (_maskedCall): the result is NaN
    exactly where 'linear' is, and elsewhere near an invalid node it depends on the fill (_fillInvalid).
    The mask is not a nicety: a path through an infeasible cell must go non-finite *without* leaving the
    rectangle, which is what the CRRA path solve's `strict` check reads. """
    xo, yo = np.argsort(x), np.argsort(y)
    grid = (np.asarray(x)[xo], np.asarray(y)[yo])
    values = np.asarray(z, dtype = float)
    values = values[np.ix_(xo, yo)] if values.ndim == 2 else values[xo][:, yo]

    def build(vals, method):
        g = interpolate.RegularGridInterpolator(grid, vals, method = method,
                                                bounds_error = False, fill_value = None)
        def call(xNew, yNew):
            a, b = np.broadcast_arrays(np.asarray(xNew, dtype = float), np.asarray(yNew, dtype = float))
            av, bv = a.ravel(), b.ravel()
            # A non-finite *coordinate* (an infeasible candidate state, which the callers do produce) must
            # come back NaN. 'linear' does that on its own; the spline methods build their splines lazily
            # per axis at call time, so a NaN coordinate makes axis 0 return NaN and axis 1 then raises on
            # constructing a spline over it. Evaluate the finite coordinates only.
            fin = np.isfinite(av) & np.isfinite(bv)
            if fin.all():
                out = g(np.stack([av, bv], axis = -1))
            else:
                probe = g(np.stack([np.where(fin, av, grid[0][0]),
                                    np.where(fin, bv, grid[1][0])], axis = -1))
                out = np.where(fin.reshape(fin.shape + (1,)*(probe.ndim-1)), probe, np.nan)
            return out.reshape(a.shape + out.shape[1:])
        return call

    if kind == 'linear':
        return build(values, kind)
    valid = _validityMask(values, 2)
    f = build(_fillInvalid(values, valid), kind)
    if valid.all():
        return f
    return _maskedCall(f, build(np.where(valid, 0., np.nan), 'linear'))


def _fixedKnots(xv, k, every):
    """ Interior knots at every `every`-th valid node, for the fixed-knot branch of _splineAlongAxis0.

    Depends only on the node positions and the validity mask -- never on the values. That is the whole
    point: with the knots pinned, the least-squares fit is a LINEAR map of y, so at a fixed mask the
    smoothed profile moves continuously with whatever parameter produced y. Placing them on the valid
    nodes also satisfies Schoenberg-Whitney by construction (each knot interval holds `every` data
    points), which fixed positions in grid coordinates would not once a column carries NaNs.

    Returns an empty array when the column is too short to support interior knots; LSQUnivariateSpline
    then fits a single degree-k polynomial, which is still linear in y. """
    idx = np.arange(every, len(xv)-every, every)
    if len(idx) == 0 or len(xv) < len(idx)+k+1:
        return np.empty(0)
    return xv[idx]


def _splineAlongAxis0(x, y, s, k, derivative, knots = None):
    """ Fit a spline in x along axis 0 of y, once per trailing index, and evaluate either it or its
    derivative back at x. Shared by griddedSmooth1D/griddedGradient1D so the two differ only in that
    final choice.

    knots selects HOW the spline is determined, and the two branches differ in a way that matters far
    beyond smoothing quality:
      None      -- UnivariateSpline(s=s): FITPACK chooses the number and placement of interior knots to
                   meet the residual bound s. That count is an integer chosen from the DATA, so it flips
                   discontinuously as y moves, and the fitted profile jumps with it.
      int m     -- LSQUnivariateSpline with interior knots at every m-th valid node (_fixedKnots). The
                   knots no longer depend on y, so the map is linear and the output is continuous.
    Use the fixed branch wherever the result is differentiated with respect to a model parameter, or
    feeds a root problem in one: the adaptive branch's knot flips read as genuine discontinuities there.

    NaNs are dropped per column rather than poisoning the fit. This is not defensive padding: the callers'
    grids carry genuinely infeasible cells (a candidate state for which no economic equilibrium exists),
    and the honest treatment is to fit the surviving points and return NaN where the input was NaN --
    matching the precedent's own dropna handling. A column with too few valid points to support degree k
    is returned as all-NaN rather than silently fitted at a lower degree. """
    x = np.asarray(x, dtype = float)
    y = np.asarray(y, dtype = float)
    if y.shape[:1] != x.shape:
        raise ValueError(f"y's first axis {y.shape[:1]} must match x {x.shape}.")
    order = np.argsort(x)
    xs = x[order]
    flat = y[order].reshape(len(xs), -1)
    out = np.full(flat.shape, np.nan)
    for j in range(flat.shape[1]):
        col = flat[:, j]
        ok = ~np.isnan(col)
        if ok.sum() <= k:                       # both spline classes need > k points
            continue
        if knots is None:
            sp = interpolate.UnivariateSpline(xs[ok], col[ok], s = s, k = k)
        else:
            sp = interpolate.LSQUnivariateSpline(xs[ok], col[ok],
                                                 _fixedKnots(xs[ok], k, knots), k = k)
        out[ok, j] = (sp.derivative() if derivative else sp)(xs[ok])
    inv = np.argsort(order)                     # undo the sort applied above
    return out.reshape(y.shape)[inv]


def griddedSmooth1D(x, y, s = 1e-5, k = 3, knots = None):
    """ Smoothing spline through gridded (x, y), evaluated back at x. y may be (M,) or (M, ...); each
    trailing index is smoothed independently along axis 0. NaNs pass through as NaN (see
    _splineAlongAxis0). s: smoothing factor -- larger smooths harder; 0 would interpolate exactly.

    knots: None keeps the adaptive-knot smoothing spline (s applies); an int m pins the knots at every
    m-th valid node and fits by least squares, and s is then IGNORED: every value of s gives the same fit,
    so a caller's nonzero s only switches the smoothing on (the policy smoothers' `smooth`). See
    _splineAlongAxis0 -- the fixed branch is the one whose output is continuous in the data, and it is also
    ~2.4x faster, since FITPACK no longer iterates over knot counts. """
    return _splineAlongAxis0(x, y, s, k, derivative = False, knots = knots)


def griddedGradient1D(x, y, s = 1e-4, k = 3, knots = None):
    """ dy/dx along a grid axis, as the derivative of a smoothing spline through (x, y). y may be (M,) or
    (M, ...) -- e.g. (M, ni) for a per-household-type quantity -- with each trailing index differentiated
    independently along axis 0. NaNs pass through as NaN (see _splineAlongAxis0).

    Use where no closed form exists. Where one does, prefer it: a grid derivative silently includes every
    channel that varies along the grid, which is not always what a first-order condition wants (see
    base.py's dlnc2i_dτ for a case where that distinction is load-bearing rather than cosmetic).

    knots: as griddedSmooth1D. Callers that pass s=0 (an interpolating spline, whose knots ARE the data
    points and so already do not depend on the values) need nothing here. """
    return _splineAlongAxis0(x, y, s, k, derivative = True, knots = knots)

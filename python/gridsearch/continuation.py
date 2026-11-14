""" Numerical continuation over a grid of parameter values.

Model-agnostic like the rest of this package: nothing here knows what the parameter *is* or what `solve`
does with it. The problem it addresses is a sequence of expensive, ill-conditioned solves indexed by a
parameter, where each solve needs a starting point and the previous solves are the only information about
where that starting point should be.

  marchGrid -- visit a grid of parameter values in continuation order, extrapolating each solve's starting
               point from the ones already completed, and recovering from a failure by stepping to an
               intermediate value first.

Three properties the callers here depend on:

  * **Anchored, bidirectional order.** `anchor` names the value to solve first; the grid is then walked
    *outward from it in both directions*, each direction carrying its own history. This is not a
    convenience -- a caller whose problem is only well-posed at one interior point (e.g. a preference
    parameter at which a different, cheaper solver applies) has nowhere else to start.
  * **Extrapolation in the parameter value, not the grid index.** Step-halving inserts values off the
    requested grid, so index spacing stops being meaningful the first time a point resists.
  * **Failure is recorded, not raised** (unless `stopOnFail`). A sweep of this kind runs for hours; one
    unsolvable value must not discard the points already bought.

The starting point passed to `solve` should be an *unbounded* coordinate wherever the caller has one.
Extrapolating a positivity- or interval-constrained parameter in levels can step across its bound, and
the resulting starting point is not merely poor but invalid; the same extrapolation in the log/logit
coordinate cannot. This function does not enforce that -- it cannot see the constraint -- but every
docstring here assumes it.
"""
import numpy as np


def extrapolateX(history, value, degree = 1):
    """ Starting point for `value` from solves already completed. history: list of (value, x) pairs in
    visit order, most recent last; x arrays of a common shape (k,). Returns (x0, source) with source one
    of 'none' (empty history), 'carry' (the most recent x, used whenever the degree available is 0) or
    'extrap'.

    The polynomial through the last degree+1 points is evaluated at `value` componentwise. Nodes are
    centred on `value` before the Vandermonde solve, so the answer is read off the constant term and the
    system stays well conditioned no matter how far the nodes sit from the origin -- which matters because
    step-halving can put two nodes arbitrarily close together. Nodes closer than 1e-12 fall back to
    'carry' rather than solving a singular system. """
    if not history:
        return None, 'none'
    n = min(int(degree), len(history)-1)
    if n < 1:
        return np.asarray(history[-1][1], dtype = float), 'carry'
    pts = history[-(n+1):]
    u = np.array([p[0] for p in pts], dtype = float) - value
    if np.min(np.abs(u[:, None] - u[None, :]) + np.eye(len(u))) < 1e-12:
        return np.asarray(history[-1][1], dtype = float), 'carry'
    xs = np.array([np.asarray(p[1], dtype = float) for p in pts])   # (n+1, k)
    coef = np.linalg.solve(np.vander(u, n+1), xs)                   # rows: u^n ... u^0
    return coef[-1], 'extrap'                                       # the polynomial evaluated at u=0


def _visitOrder(grid, anchor):
    """ (anchorValue, leftValues, rightValues) for a sorted grid. anchor: a value in the grid (matched to
    the nearest node) or None for the leftmost. left runs *outward*, i.e. descending. """
    g = np.asarray(grid, dtype = float)
    if np.any(np.diff(g) <= 0):
        raise ValueError('grid must be strictly increasing (got {})'.format(g))
    i = 0 if anchor is None else int(np.argmin(np.abs(g-anchor)))
    return g[i], list(g[:i][::-1]), list(g[i+1:])


def marchGrid(grid, solve, x0 = None, anchor = None, degree = 1, maxHalvings = 2, minStep = 1e-9,
              onPoint = None, stopOnFail = False):
    """ Solve at every value of `grid`, in continuation order from `anchor` outward.

    solve(value, x0) -> result: must return a dict carrying key 'x', the coordinate that seeds the next
        solve (see the module docstring on which coordinate that should be). Raising signals failure;
        no exception type is privileged, since an inner solver may report non-convergence any way it likes.
    x0: starting point for the anchor only. Every other point gets one from the history.
    anchor: grid value to solve first (nearest node); None -> the leftmost, i.e. a plain left-to-right pass.
    degree: extrapolation degree. 1 (linear, through the last two solves) is the default because the
        second derivative of the solution path is exactly what a continuation problem has no estimate of
        until it is already in trouble; 0 reuses the previous x unchanged.
    maxHalvings: how many intermediate values may be inserted between the last success and a resisting
        target before it is declared failed. Intermediates are solved and *kept*: they enter the history
        (so the retry extrapolates from closer in) and are recorded with requested=False.
    minStep: an intermediate closer than this to the last success is not attempted -- the point is
        genuinely unsolvable rather than merely far away.
    onPoint(record): called after every point, successful or not, before the march continues. This is the
        persistence hook: a sweep long enough to need this function is long enough to crash, and writing
        from here caps the loss at one point.
    stopOnFail: abandon the current direction on the first failed *requested* value. Off by default -- the
        other direction is unaffected either way, and a failure need not mean every value beyond it fails.

    Returns {'records', 'failures', 'history'}: records one dict per attempted value in visit order
    (keys 'value', 'ok', 'requested', 'x0Source', 'halvings', and on success 'x'/'result', on failure
    'error'); failures the subset with ok=False; history the (value, x) pairs kept, sorted by value. """
    anchorValue, left, right = _visitOrder(grid, anchor)
    records, failures = [], []

    def record(r):
        records.append(r)
        if not r['ok']:
            failures.append(r)
        if onPoint is not None:
            onPoint(r)

    def attempt(value, x0Cand, source, requested, halvings):
        """ One call to solve. Returns the result dict, or None having recorded the failure. """
        try:
            res = solve(value, x0Cand)
        except Exception as e:                      # any inner solver's non-convergence, not just one type
            record({'value': value, 'ok': False, 'requested': requested, 'x0Source': source,
                    'halvings': halvings, 'error': '{}: {}'.format(type(e).__name__, e)})
            return None
        record({'value': value, 'ok': True, 'requested': requested, 'x0Source': source,
                'halvings': halvings, 'x': np.asarray(res['x'], dtype = float), 'result': res})
        return res

    def candidates(history, value):
        """ Starting points to try at `value`, best first: the extrapolation, then the previous x
        unchanged. The second is not redundant -- extrapolating across a kink in the solution path
        overshoots, and the un-extrapolated point is then both closer and cheaper than halving the step.

        With no history there is nothing to extrapolate from, and x0=None is passed through so `solve`
        applies its own default. Returning nothing instead would leave the anchor unattempted. """
        if not history:
            return [(None, 'default')]
        out, seen = [], []
        for x0Cand, source in (extrapolateX(history, value, degree),
                               extrapolateX(history, value, 0)):
            if x0Cand is None or any(np.allclose(x0Cand, s) for s in seen):
                continue
            seen.append(x0Cand)
            out.append((x0Cand, source))
        return out

    def solvePoint(value, history, requested, halvings = 0):
        """ Solve at `value`, halving the step toward it on failure. Appends to history on success. """
        for x0Cand, source in candidates(history, value):
            res = attempt(value, x0Cand, source, requested, halvings)
            if res is not None:
                history.append((value, np.asarray(res['x'], dtype = float)))
                return True
        if halvings >= maxHalvings or not history:
            return False
        mid = 0.5*(history[-1][0] + value)
        if min(abs(mid-history[-1][0]), abs(value-mid)) < minStep:
            return False
        if not solvePoint(mid, history, False, halvings+1):
            return False
        return solvePoint(value, history, requested, halvings+1)

    anchorHistory = []
    if x0 is not None:
        res = attempt(anchorValue, np.asarray(x0, dtype = float), 'user', True, 0)
        if res is not None:
            anchorHistory.append((anchorValue, np.asarray(res['x'], dtype = float)))
    else:
        solvePoint(anchorValue, anchorHistory, True)
    if not anchorHistory:
        raise RuntimeError('the anchor value {} failed to solve; the march has no starting point '
                           '(records: {})'.format(anchorValue, records))

    history = list(anchorHistory)
    for values in (right, left):
        h = list(anchorHistory)                     # each direction extrapolates along its own path only
        for value in values:
            if not solvePoint(value, h, True) and stopOnFail:
                break
        history += [p for p in h if p[0] != anchorValue]
    return {'records': records, 'failures': failures, 'history': sorted(history, key = lambda p: p[0])}

""" Cartesian product of named 1d grids, and the flat <-> ND mapping that goes with it.

Motivation. A grid search over a choice variable at each point of a state grid wants two things at once:
to evaluate the model's equations *once*, on flat vectors covering every (choice, state) combination
(base.py's methods are vectorized over a leading "points" axis, so this is one numpy pass), and then to
look along the choice dimension alone to locate roots/maxima (roots1d, which consumes an (M, N) matrix
of N independent functions sampled on one shared M-node grid).

The inspiration code did this with a pandas MultiIndex, recovering the matrix via `unstack(level=...)`.
That is not needed: when the product is built here, in C-order, the flat layout is fixed by construction,
so `reshape` inverts it *exactly* and cheaply. No sorting, no label lookups, no pandas in the hot path.
Nothing here forbids attaching labels later -- reporting can build an index from `.grids` -- but the
solve path stays plain ndarrays end to end.

Ordering convention. Axes follow the order the grids were given, and `flat` is C-order (last axis varies
fastest), matching `np.reshape`. Nothing depends on the choice variable being first: `asColumns` moves
whichever axis is named to the front. Putting it first is still the natural reading order.

    g = CartesianGrid(tau = tauGrid, s_ = sGrid)
    z = someModelEquation(**g.flat)     # (g.size,)  -- one vectorized pass
    zc = g.asColumns(z, 'tau')          # (M_tau, N_states)  -- exact reshape
    ...                                 # roots1d.selectMax(tauGrid, zc) -> (N_states,)
    sol = sel['x'].reshape(g.stateShape('tau'))
"""
import numpy as np


class CartesianGrid:
    """ The Cartesian product of named 1d grids.

    Construct with a dict or with keyword arguments (`CartesianGrid({'tau': a, 's_': b})` and
    `CartesianGrid(tau = a, s_ = b)` are equivalent -- the dict form exists because grid names are built
    programmatically as often as they are typed, and need not be valid Python identifiers).

    Grids are not required to be sorted here: only the axis actually searched over has to be, and that is
    roots1d's precondition to enforce (it does, loudly), not this class's -- a state axis may legitimately
    be unsorted or even non-numeric. """

    def __init__(self, grids = None, **kwargs):
        g = dict(grids or {}) | kwargs
        if not g:
            raise ValueError("CartesianGrid needs at least one grid.")
        self.grids = {k: np.asarray(v) for k, v in g.items()}
        for k, v in self.grids.items():
            if v.ndim != 1:
                raise ValueError(f"grid {k!r} must be 1-dimensional, got shape {v.shape}.")
            if v.size == 0:
                raise ValueError(f"grid {k!r} is empty.")
        self._flat = None # built on first use; see `flat`

    # ---- structure -------------------------------------------------------------------------------
    @property
    def names(self):
        """ Axis names, in axis order. """
        return tuple(self.grids)

    @property
    def shape(self):
        """ ND shape of the product, in axis order. """
        return tuple(v.size for v in self.grids.values())

    @property
    def size(self):
        """ Total number of points (= len of every vector in `flat`). """
        return int(np.prod(self.shape))

    def __len__(self):
        return self.size

    def __repr__(self):
        dims = ', '.join(f'{k}={v.size}' for k, v in self.grids.items())
        return f'CartesianGrid({dims}; size={self.size})'

    def axis(self, name):
        """ Position of `name` among the axes. """
        try:
            return self.names.index(name)
        except ValueError:
            raise KeyError(f"{name!r} is not a grid here; have {self.names}.") from None

    def values(self, name):
        """ The 1d grid stored under `name` (the nodes, not the expanded flat vector). """
        if name not in self.grids:
            raise KeyError(f"{name!r} is not a grid here; have {self.names}.")
        return self.grids[name]

    def stateShape(self, name):
        """ The ND shape left once `name`'s axis is removed -- i.e. the shape of one solution per state
        combination, which is what a search along `name` returns. Empty tuple when `name` is the only
        axis (the result is then a scalar-shaped array, consistent with numpy's own conventions). """
        a = self.axis(name)
        return self.shape[:a] + self.shape[a+1:]

    # ---- flat <-> ND -----------------------------------------------------------------------------
    @property
    def flat(self):
        """ {name: (size,) ndarray}: every axis broadcast over the full product, C-order.

        Feed straight into any function vectorized over a leading points axis --
        `f(**g.flat)` evaluates the whole product in one pass. Built once and cached, since a grid is
        immutable here and the same flat vectors are reused across every evaluation of a solve. """
        if self._flat is None:
            mesh = np.meshgrid(*self.grids.values(), indexing = 'ij')
            self._flat = {k: m.ravel() for k, m in zip(self.grids, mesh)}
        return self._flat

    def reshape(self, y):
        """ (size, ...) -> shape + (...): undo the flattening, keeping any trailing axes.

        Trailing axes are carried through untouched so per-type quantities work unchanged: base.py returns
        (points, ni) arrays for anything indexed by household type, and those become shape + (ni,). """
        y = np.asarray(y)
        if y.shape[:1] != (self.size,):
            raise ValueError(f"y's first axis ({y.shape[:1]}) must be the grid size ({self.size},).")
        return y.reshape(self.shape + y.shape[1:])

    def asColumns(self, y, name):
        """ (size,) -> (M, N): `name`'s axis moved to the front, every other axis collapsed into columns.

        This is the shape roots1d works in -- M nodes of the searched-over variable, N independent
        functions sampled on them, one per state combination. Column order is C-order over the remaining
        axes *in their original relative order*, which is exactly what `stateShape` describes, so a
        per-column result reshapes back with `result.reshape(g.stateShape(name))`.

        1d input only: a genuine (M, N) matrix is the contract roots1d validates against, and silently
        folding trailing per-type axes into the columns would produce a matrix whose columns no longer
        correspond one-to-one with state combinations. Use `reshape` for those. """
        y = np.asarray(y)
        if y.ndim != 1:
            raise ValueError(f"asColumns expects a flat 1d array, got shape {y.shape}; use reshape() for "
                             "quantities with trailing axes.")
        a = np.moveaxis(self.reshape(y), self.axis(name), 0)
        return a.reshape(a.shape[0], -1)

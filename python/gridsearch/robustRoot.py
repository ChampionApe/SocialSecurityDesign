""" Bounded-root reparameterization (the note's "Robust policy search with bounds", eq:root).

Generic, model-agnostic: wraps any first-order-condition function f (only meaningful/safe to evaluate on
[l,u], e.g. a tax rate confined to [0,1]) into an unconstrained residual h(z,τ) that a gradient-based
root-finder (e.g. scipy.optimize.root/newton) can search over all of ℝ without ever calling f outside its
safe domain. Interior candidates (τ∈[l,u]) evaluate f directly and h=f. Candidates outside [l,u] are
extended with an artificial linear penalty anchored at f(l)/f(u): the sign of the penalty term is set up
so h still points back toward the relevant bound even when f itself has no interior root or isn't
monotonic there (see the doc's two boundary-detection examples). Recover the actual bounded policy from a
solved τ̃ via clip(τ̃, l, u). """
import numpy as np


def clip(τ, l = 0, u = 1):
    """ [τ]_{[l,u]} ≡ min(u, max(l,τ)) -- the bounded interior value τ is actually evaluated at. """
    return np.minimum(u, np.maximum(l, τ))


def penalty(τ, l = 0, u = 1, a0 = 1, a1 = 1):
    """ g(τ) (eq:root): 0 for τ∈[l,u], else a linear penalty-slope driver growing with distance outside
    the bound (a0 below l, a1 above u). a0, a1>0 are technical parameters (steeper = faster to detect a
    boundary solution, but too steep can hurt the root-finder's own step-size heuristics). """
    return a0*np.minimum(τ-l, 0) + a1*np.maximum(τ-u, 0)


def boundedResidual(f, l = 0, u = 1, a0 = 1, a1 = 1):
    """ Wrap f into h(z,τ) (eq:root). f is only ever evaluated at clip(τ,l,u), so it never has to be
    defined/safe outside [l,u]. Returns a callable h(τ) suitable as a root-finder's residual; τ (and f's
    return value) may be scalar or an array -- clip/penalty are elementwise, so this vectorizes for free
    over e.g. a whole time path of τ_t's, matching this codebase's "zero Python loops" style elsewhere. """
    def h(τ):
        z = f(clip(τ, l, u))
        return z - np.abs(z)*penalty(τ, l, u, a0, a1)
    return h

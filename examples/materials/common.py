"""Shared helpers for the materials engineering examples (woven preform, PU foam).

Nothing here changes the solver. These functions turn a PoreWise result into the engineering
quantities the examples report, and check the fields the solver returned:

    sample_permeability      k per sample cross-section (not per enclosing box)
    series_corrected         remove free-fluid buffer layers from an apparent permeability
    plane_fluxes             volume flux through every plane normal to the flow (conservation)
    flow_regime              Mach and pore-scale Reynolds number in lattice units
    trace_streamlines        wall-safe streamlines for figures (vectorised form of the rule in
                             benchmarks/rock_streamlines.py: interpolate only where all eight
                             surrounding cells are pore, check sub-steps, stop otherwise)
    points_in_solid          nearest-voxel test of figure points against the computational grid

Array convention as in PoreWise: (z, y, x); velocity components ux, uy, uz on the same grid.
"""
import itertools
import json
import math

import numpy as np

AXES = {"x": 2, "y": 1, "z": 0}


# ---------------------------------------------------------------------------- permeability
def sample_permeability(k_box, sample_area_fraction):
    """Permeability per sample cross-section.

    PoreWise reports k = nu <u> / F with <u> averaged over every cell of the grid, i.e. the
    superficial velocity over the whole box cross-section. When the sample occupies only part of
    that cross-section (a cylindrical core inside a square grid, the rest solid), the Darcy flux
    of the sample is the same volume flow divided by the sample area, so
        k_sample = k_box / (A_sample / A_box).
    """
    if not 0 < sample_area_fraction <= 1:
        raise ValueError("sample_area_fraction must be in (0, 1]")
    return k_box / sample_area_fraction


def series_corrected(k_apparent, length_total, length_sample, k_buffer=math.inf):
    """Permeability of a sample in series with free-fluid buffer layers (body-force driven).

    With a uniform body force and periodic pressure, the apparent permeability of layers in
    series is the length-weighted harmonic mean: L_tot / k_app = L_s / k_s + L_b / k_b.
    """
    lb = length_total - length_sample
    if lb < 0:
        raise ValueError("sample longer than the domain")
    inv = length_total / k_apparent - (lb / k_buffer if math.isfinite(k_buffer) else 0.0)
    if inv <= 0:
        raise ValueError("buffer resistance exceeds the total resistance")
    return length_sample / inv


# ---------------------------------------------------------------------------- field checks
def plane_fluxes(u, axis):
    """Volume flux through every plane normal to `axis` (sum of the normal velocity)."""
    other = tuple(a for a in range(3) if a != axis)
    q = np.asarray(u, dtype=np.float64).sum(axis=other)
    return dict(mean=float(q.mean()), cv=float(q.std() / abs(q.mean())) if q.mean() else math.inf,
                max_rel_dev=float(np.abs(q / q.mean() - 1).max()) if q.mean() else math.inf)


def flow_regime(fields, blocked, nu, pore_length_lu):
    """Mach number and pore-scale Reynolds number (interstitial mean speed x pore length / nu)."""
    speed = np.sqrt(sum(np.asarray(fields[c], dtype=np.float64) ** 2 for c in ("ux", "uy", "uz")))
    pore = ~np.asarray(blocked, bool)
    mean_pore = float(speed[pore].mean())
    return dict(mach_max=float(speed.max() * math.sqrt(3)), u_mean_pore_lu=mean_pore,
                reynolds_pore=mean_pore * pore_length_lu / nu, pore_length_lu=pore_length_lu)


# ---------------------------------------------------------------------------- streamlines
def _velocity(points, U, blocked):
    """Trilinear velocity at points (N, 3) in index coordinates (z, y, x).

    `ok` is False where any of the eight surrounding cells is solid or outside the grid: there
    the interpolant would mix wall and fluid values, so the point is not traced.
    """
    shape = np.array(blocked.shape)
    base = np.floor(points).astype(np.int64)
    ok = np.all(base >= 0, axis=1) & np.all(base + 1 < shape, axis=1)
    b = np.clip(base, 0, shape - 2); f = points - b
    vel = np.zeros_like(points)
    for corner in itertools.product((0, 1), repeat=3):
        i = b + np.array(corner)
        ok &= ~blocked[i[:, 0], i[:, 1], i[:, 2]]
        w = np.prod(np.where(np.array(corner, bool), f, 1 - f), axis=1)
        vel += w[:, None] * U[:, i[:, 0], i[:, 1], i[:, 2]].T
    return vel, ok


def trace_streamlines(fields, blocked, seeds, *, step=0.2, max_steps=4000, direction=1.0, sub_step=0.05):
    """Integrate streamlines from `seeds` (N, 3), index coordinates (z, y, x), arc-length Euler.

    Returns a list of (n_i, 3) arrays and the stop reason of each line. A line stops where the
    interpolation stencil touches a solid cell or the grid edge, at stagnation, or at max_steps;
    every sub-step of every step (spacing <= `sub_step` voxel) is checked before it is kept, so
    no retained point has a solid cell among its eight neighbours.
    """
    blocked = np.asarray(blocked, bool)
    U = np.stack([np.asarray(fields[c], dtype=np.float64) for c in ("uz", "uy", "ux")])
    p = np.asarray(seeds, dtype=np.float64).copy(); n = len(p)
    active = np.ones(n, bool); reason = np.array(["maximum_steps"] * n, dtype=object)
    path = [p.copy()]; alive = [active.copy()]
    subs = np.linspace(sub_step / step, 1.0, max(1, int(math.ceil(step / sub_step))))
    _, ok0 = _velocity(p, U, blocked); reason[~ok0] = "seed_at_wall"; active &= ok0
    for _ in range(max_steps):
        if not active.any():
            break
        idx = np.flatnonzero(active)
        v, ok = _velocity(p[idx], U, blocked)
        speed = np.linalg.norm(v, axis=1)
        stop = ~ok | (speed < 1e-18)
        delta = np.zeros_like(v); good = ~stop
        delta[good] = direction * step * v[good] / speed[good, None]
        for a in subs:                                            # check every sub-step first
            _, oks = _velocity(p[idx] + a * delta, U, blocked)
            stop |= ~oks
        reason[idx[stop & (speed < 1e-18)]] = "stagnation"
        reason[idx[stop & ~(speed < 1e-18)]] = "wall_or_edge"
        keep = idx[~stop]
        p[keep] += delta[~stop]
        active[idx[stop]] = False
        path.append(p.copy()); alive.append(active.copy())
    P = np.stack(path); A = np.stack(alive)                        # (steps+1, n, 3), (steps+1, n)
    lines = []
    for k in range(n):
        if reason[k] == "seed_at_wall":                            # the seed itself fails the stencil rule: no line
            lines.append(np.empty((0, 3))); continue
        last = int(np.flatnonzero(A[:, k]).max()) if A[:, k].any() else 0
        lines.append(P[: last + 1, k])
    return lines, list(reason)


def trace_both_ways(fields, blocked, seeds, **kw):
    """Streamlines through each seed, upstream and downstream joined into one line."""
    fwd, rf = trace_streamlines(fields, blocked, seeds, direction=1.0, **kw)
    bwd, rb = trace_streamlines(fields, blocked, seeds, direction=-1.0, **kw)
    return [np.concatenate([b[::-1], f[1:]]) if len(f) else f for f, b in zip(fwd, bwd)], list(zip(rb, rf))


def points_in_solid(points, blocked):
    """Fraction-ready boolean: is the nearest cell of each point (z, y, x) solid?"""
    i = np.rint(np.asarray(points)).astype(np.int64)
    ok = np.all((i >= 0) & (i < np.array(blocked.shape)), axis=1)
    out = np.zeros(len(i), bool); out[ok] = np.asarray(blocked)[i[ok, 0], i[ok, 1], i[ok, 2]]
    return out


def stencil_touches_solid(points, blocked):
    """True where any of the eight cells around a point is solid (the tracer's stop rule)."""
    _, ok = _velocity(np.asarray(points, dtype=np.float64), np.zeros((3,) + blocked.shape), np.asarray(blocked, bool))
    return ~ok


# ---------------------------------------------------------------------------- output
def run_summary(run):
    """The scalar provenance of a PoreWise run, JSON-ready (fields removed)."""
    keep = {}
    for k, v in run.items():
        if k in ("ux", "uy", "uz", "rho", "convergence_history"):
            continue
        if isinstance(v, (bool, int, float, str)) or v is None:
            keep[k] = v
        elif isinstance(v, dict):
            keep[k] = json.loads(json.dumps(v, default=float))
        elif isinstance(v, (list, tuple)):
            keep[k] = json.loads(json.dumps(list(v), default=float))
    return keep

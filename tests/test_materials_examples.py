"""Regression tests for the materials engineering examples (examples/materials)."""
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MAT = ROOT / "examples" / "materials"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod; spec.loader.exec_module(mod); return mod


common = _load("materials_common", MAT / "common.py")
geometry = _load("woven_geometry", MAT / "woven_preform" / "geometry.py")
prepare = _load("foam_prepare", MAT / "pu_foam" / "prepare.py")
rock = _load("rock_streamlines_ref", ROOT / "benchmarks" / "rock_streamlines.py")


# ----------------------------------------------------------------------------- woven preform
@pytest.mark.parametrize("state", ["relaxed", "compacted"])
def test_plain_weave_topology_and_mould_contact(state):
    wv = geometry.STATES[state]
    assert geometry.plain_weave_alternates(wv)
    solid, labels, info = geometry.voxelise(wv, 0.04)
    assert info["tow_free_layers"] == 0                       # no slot between preform and mould
    assert abs(info["mould_offset_mm"]) <= 0.02 + 1e-12        # mould within half a voxel of the crests
    assert info["overlap_before_correction_fraction"] < 0.005
    assert set(np.unique(labels)) <= {0, 1, 2, 3, 4, 255}
    assert labels[0].min() == 255 and labels[-1].min() == 255


@pytest.mark.parametrize("state", ["relaxed", "compacted"])
def test_fibre_bookkeeping_is_consistent(state):
    wv = geometry.STATES[state]
    _, _, info = geometry.voxelise(wv, 0.02)
    # fibre volume fraction = intra-tow fraction x tow volume fraction, to within the voxelisation error
    assert info["fibre_volume_fraction"] == pytest.approx(wv.intra_tow_vf * info["tow_volume_fraction"], rel=0.03)
    assert 195 < info["areal_weight_g_m2"] < 215                # a 3K, 5 ends/cm carbon plain weave


def test_weave_is_periodic_in_plane():
    wv = geometry.STATES["relaxed"]; L = wv.cell_mm
    rng = np.random.default_rng(0); x, y, z = rng.uniform(0, L, 500), rng.uniform(0, L, 500), rng.uniform(-0.2, 0.2, 500)
    a = geometry.tow_margins(x, y, z, wv, 0.0, (0.0, 0.0)); b = geometry.tow_margins(x + L, y - L, z, wv, 0.0, (0.0, 0.0))
    for k in ("warp", "weft"):
        np.testing.assert_allclose(a[k][0], b[k][0], atol=1e-12)


def test_interference_correction_removes_overlap():
    wv = geometry.Weave("sinusoid", spacing_mm=2.0, width_mm=1.85, intra_tow_vf=0.72, beta=0.0)   # interferes
    g = np.linspace(0, wv.cell_mm, 120, endpoint=False); z = np.linspace(-0.12, 0.12, 40)
    Z, Y, X = np.meshgrid(z, g, g, indexing="ij")
    mg = geometry.tow_margins(X, Y, Z, wv, 0.0, (0.0, 0.0))
    warp, weft, both = geometry.corrected_inside(mg)
    assert both.any()                                          # the sinusoid does interfere
    assert not (warp & weft).any()                             # after correction no point is in both
    assert ((mg["warp"][0] > 0) | (mg["weft"][0] > 0)).sum() == (warp | weft).sum()   # union unchanged


# ----------------------------------------------------------------------------- normalisation
def test_sample_permeability_and_series_correction():
    assert common.sample_permeability(2.0, 0.5) == pytest.approx(4.0)
    with pytest.raises(ValueError):
        common.sample_permeability(1.0, 0.0)
    # two layers in series: 3 units of k=1 and 1 unit of k=4 -> harmonic mean 4/(3/1+1/4)
    k_app = 4 / (3 / 1 + 1 / 4)
    assert common.series_corrected(k_app, 4, 3, k_buffer=4) == pytest.approx(1.0)
    assert common.series_corrected(k_app * 0 + 4 / 3, 4, 3) == pytest.approx(1.0)


def test_pipe_in_sleeve_matches_poiseuille():
    """A straight circular pipe in a square box: k per pipe area tends to R^2/8."""
    from porewise import lbm_stokes_3d
    n, R = 26, 10.0
    inside = prepare.disc_mask((n, n), ((n - 1) / 2, (n - 1) / 2), R)
    solid = np.broadcast_to(~inside[None], (4, n, n)).copy()
    tau = 1.0; nu = (tau - 0.5) / 3
    run = lbm_stokes_3d(solid, F_x=0.0, F_z=1e-6, backend="numpy", collision="trt", tau=tau, n_steps_max=20000,
                        conv_tol=1e-7, verbose=False)
    k_pipe = common.sample_permeability(run["u_z_mean_total"] * nu / 1e-6, inside.mean())
    r_eff = math.sqrt(inside.sum() / math.pi)                  # staircase disc of the same area
    assert k_pipe == pytest.approx(r_eff ** 2 / 8, rel=0.08)


# ----------------------------------------------------------------------------- streamlines
def _channel_field(nz=6, ny=12, nx=40):
    blocked = np.zeros((nz, ny, nx), bool); blocked[:, :2] = blocked[:, -2:] = True
    y = np.arange(ny) - (ny - 1) / 2; prof = np.clip(1 - (y / 4.5) ** 2, 0, None)
    ux = np.broadcast_to(prof[None, :, None], blocked.shape).copy(); ux[blocked] = 0
    return blocked, dict(ux=ux, uy=np.zeros_like(ux), uz=np.zeros_like(ux))


def test_vectorised_tracer_matches_reference_rule():
    blocked, f = _channel_field()
    seeds = np.array([[2.5, 4.0, 1.0], [2.5, 6.0, 1.0], [2.5, 7.2, 2.0]])
    lines, reasons = common.trace_streamlines(f, blocked, seeds, step=0.2, max_steps=500)
    for s, line in zip(seeds, lines):
        ref, _ = rock.integrate(s, f, blocked, step=0.2, max_steps=500)
        assert len(line) == len(ref)
        np.testing.assert_allclose(line, ref, atol=1e-9)


def test_traced_points_never_touch_solid():
    blocked, f = _channel_field()
    f["uy"] = 0.15 * np.ones_like(f["ux"]); f["uy"][blocked] = 0          # drift towards the wall
    seeds = np.array([[2.5, y, 1.0] for y in np.linspace(0.5, 10.5, 21)])          # some seeds at or in the wall
    lines, reasons = common.trace_streamlines(f, blocked, seeds, step=0.2, max_steps=2000)
    assert any(r == "seed_at_wall" for r in reasons) and all(len(l) == 0 for l, r in zip(lines, reasons) if r == "seed_at_wall")
    pts = np.concatenate([l for l in lines if len(l)])
    assert not common.stencil_touches_solid(pts, blocked).any()
    assert not common.points_in_solid(pts, blocked).any()
    assert all(r in ("wall_or_edge", "maximum_steps", "stagnation", "seed_at_wall") for r in reasons)


def test_plane_fluxes_detects_conservation():
    u = np.ones((5, 4, 6)); assert common.plane_fluxes(u, 2)["cv"] == pytest.approx(0.0)
    u[:, :, 3] *= 1.1; assert common.plane_fluxes(u, 2)["max_rel_dev"] > 0.05


# ----------------------------------------------------------------------------- foam domain
def _synthetic_core(nz=12, n=40, seed=1):
    rng = np.random.default_rng(seed); centre = ((n - 1) / 2, (n - 1) / 2); r = 16.0
    inside = prepare.disc_mask((n, n), centre, r)
    skel = (rng.random((nz, n, n)) < 0.06) & inside[None]
    skel[:, int(centre[0]), int(centre[1] - r + 0.5)] = True; skel[:, int(centre[0]), int(centre[1] + r - 0.5)] = True
    skel[:, int(centre[0] - r + 0.5), int(centre[1])] = True; skel[:, int(centre[0] + r - 0.5), int(centre[1])] = True
    return skel, inside


def test_cylinder_fit_and_domain_ends():
    skel, inside = _synthetic_core()
    c, r = prepare.cylinder(skel)
    assert c == pytest.approx(((40 - 1) / 2, (40 - 1) / 2), abs=0.6) and r == pytest.approx(16.0, abs=0.6)
    s_b, i_b = prepare.build_domain(skel, c, r, ends="buffer", buffer_voxels=3)
    assert s_b.shape[0] == skel.shape[0] + 6 and i_b["sample_slices"] == skel.shape[0]
    assert not s_b[:3][:, prepare.disc_mask(skel.shape[1:], c, r)].any()      # buffers are free fluid in the sleeve
    s_m, i_m = prepare.build_domain(skel, c, r, ends="mirror")
    np.testing.assert_array_equal(s_m[:skel.shape[0]], s_m[skel.shape[0]:][::-1])
    assert (s_m[:, ~prepare.disc_mask(skel.shape[1:], c, r)]).all()             # sleeve is solid


def test_porosity_preserving_coarsening():
    skel, inside = _synthetic_core(nz=16)
    coarse, thr = prepare.coarsen(skel, 2, "porosity", inside)
    ins_c = inside[::2, ::2]
    assert coarse[:, ins_c].mean() == pytest.approx(skel[:, inside].mean(), abs=0.01)
    quarter, _ = prepare.coarsen(skel, 2, "quarter", inside)
    assert quarter[:, ins_c].mean() >= coarse[:, ins_c].mean()

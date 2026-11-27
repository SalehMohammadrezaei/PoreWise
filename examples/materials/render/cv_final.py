"""Blender scene specs for the final portfolio images of the two materials examples.

    viz=python        # an environment with numpy, scipy, pyvista and matplotlib
    $viz examples/materials/render/cv_final.py woven A out_dir [--res 1000 750 --samples 48]
    $viz examples/materials/render/cv_final.py foam  A out_dir
    LD_LIBRARY_PATH=/usr/lib/wsl/lib blender -b --factory-startup -P examples/materials/render/bl_scene.py -- out_dir/<name>.json
    python examples/materials/render/compose.py out_dir/<name>_rgba.png out_dir/<name>.png cv

Woven: two plies of the compacted plain weave over 2 x 2 unit cells (8 mm), true scale. The upper
ply is removed for x > 4 mm in the picture only (the flow was computed with both plies in place).
Streamlines of the 10 um solution, colour |u|/U_D: lower-ply channels over the whole sample (hidden
where the upper ply covers them), upper-ply channels where the upper ply is drawn.

Foam: the simulated section of the micro-CT core (full 52 mm cross-section, the central 15 mm
of its length, scan resolution), struts as segmented; a 90-degree sector is left out of the drawing
and holds streamlines of the computed field, colour |u|/U_D.

All streamline points are checked with the wall-safe rule (no point in or next to a solid).
"""
import argparse, json, math, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent; MAT = HERE.parent; ROOT = MAT.parents[1]
sys.path.insert(0, str(MAT)); sys.path.insert(0, str(MAT / "woven_preform"))
from common import trace_both_ways, points_in_solid, stencil_touches_solid  # noqa: E402

FIELDS = ROOT / "results" / "materials"; CACHE = FIELDS / "render_cache" / "cv_final"
GREY_WARP, GREY_WEFT, GREY_FOAM = [0.53, 0.55, 0.58], [0.66, 0.68, 0.70], [0.66, 0.68, 0.70]


def _cmap(v):
    import matplotlib
    return matplotlib.colormaps["viridis"](np.clip(v, 0, 1))[:, :3]


def _save_lines(name, P, S, vmax):
    CACHE.mkdir(parents=True, exist_ok=True); f = CACHE / f"{name}_lines.npz"
    np.savez(f, pts=np.concatenate(P), sizes=np.array([len(p) for p in P]), col=np.concatenate([_cmap(s / vmax) for s in S]))
    return str(f)


def _segments(lines, keep, min_extent, axis_cols):
    """Split each line where keep(points) is False; keep runs whose extent along axis_cols >= min_extent (index units)."""
    for ln in lines:
        if len(ln) < 2:
            continue
        m = keep(ln); e = np.flatnonzero(np.diff(np.r_[0, m.astype(int), 0]))
        for b0, b1 in zip(e[::2], e[1::2]):
            seg = ln[b0:b1]
            if len(seg) > 2 and np.ptp(seg[:, axis_cols]) >= min_extent:
                yield seg


# ---------------------------------------------------------------- woven
def woven_meshes(state, step=0.005, cells=2, cut_x=4.0):
    import pyvista as pv
    from geometry import STATES, stack_layout, tow_margins, corrected_inside
    wv = STATES[state]; lay = stack_layout(wv); W = cells * wv.cell_mm
    xs = np.arange(-step, W + 1.5 * step, step)
    zs = np.arange(lay["mould_bottom_mm"] - 2 * step, lay["mould_top_mm"] + 2 * step, step)
    Z, Y, X = np.meshgrid(zs, xs, xs, indexing="ij")
    outside = (X < 0) | (X > W) | (Y < 0) | (Y > W)
    out = []; top = len(lay["ply_z0_mm"]) - 1
    for p, (z0, sh) in enumerate(zip(lay["ply_z0_mm"], lay["ply_shift"])):
        mg = tow_margins(X, Y, Z, wv, z0, sh); warp, weft, _ = corrected_inside(mg)
        for kind in ("warp", "weft"):
            m = np.where(np.isfinite(mg[kind][0]), mg[kind][0], -0.05)
            other = warp if kind == "weft" else weft
            f = np.where(other & (m > 0), -1e-3, m).astype(np.float32)
            f[outside] = -1.0                                                 # closed (cut) tow ends at the sample edge
            if p == top and cut_x is not None:
                f[X > cut_x] = -1.0                                           # upper ply left out of the picture here
            g = pv.ImageData(dimensions=(len(xs), len(xs), len(zs)), spacing=(step,) * 3, origin=(xs[0], xs[0], zs[0])); g["f"] = f.ravel()
            mesh = g.contour([0.0], scalars="f").smooth_taubin(n_iter=20, pass_band=0.08).compute_normals(auto_orient_normals=True)
            fn = CACHE / f"woven_{state}_p{p}_{kind}.ply"; CACHE.mkdir(parents=True, exist_ok=True); mesh.save(fn)
            out.append((str(fn), kind, p)); print(fn.name, mesh.n_cells, "triangles")
    return out, wv, lay


def woven_lines(state, wv, lay, cells=2, cut_x=4.0, n_low=40, n_up=26, dx_um=10.0, vmax=30.0):
    d = np.load(FIELDS / "woven" / f"{state}_mould_dx{dx_um:g}um_x.npz"); solid, lab = d["solid"], d["labels"]; dx = dx_um * 1e-3
    nzg, n = solid.shape[0], solid.shape[1]; z0 = lay["mould_bottom_mm"] + (lay["cavity_mm"] - (nzg - 2) * dx) / 2 - dx
    U_D = float(d["ux"][1:-1].astype(np.float64).mean()); T = cells + 2
    S = np.tile(solid, (1, T, T)); F = {c: np.tile(d[c], (1, T, T)) for c in ("ux", "uy", "uz")}
    zc = z0 + (np.arange(nzg) + 0.5) * dx; lo, hi = n, n + cells * n                     # drawn region, tiled indices
    to_mm = lambda q: np.c_[(q[:, 2] - n + 0.5) * dx, (q[:, 1] - n + 0.5) * dx, z0 + (q[:, 0] + 0.5) * dx]
    top_low = lay["ply_z0_mm"][0] + wv.thickness_mm + 0.01
    rng = np.random.default_rng(11); P, Sp, allp = [], [], []
    for ply, code, count in ((0, 1, n_low), (1, 3, n_up)):
        kz = np.flatnonzero((zc > lay["ply_z0_mm"][ply]) & (zc < lay["ply_z0_mm"][ply] + wv.thickness_mm))
        open_y = ~(lab[kz] == code).any(axis=(0, 2))                        # channel columns of this ply
        cand = np.argwhere(~S[kz[0]:kz[-1] + 1, lo:hi, lo:hi]) + [kz[0], lo, lo]
        cand = cand[open_y[cand[:, 1] % n]]
        if ply == 1:
            cand = cand[(cand[:, 2] - n + 0.5) * dx < cut_x - 0.3]
        k = min(10 * count, len(cand)); seeds = cand[rng.choice(len(cand), k, replace=False)] + rng.uniform(0.2, 0.8, (k, 3))
        lines, _ = trace_both_ways(F, S, seeds, step=0.4, max_steps=10000); allp += [l for l in lines if len(l)]
        if ply == 0:   # in the cut-back half, nothing above the lower ply (that space held the upper ply)
            keep = lambda q: (q[:, 1] >= lo) & (q[:, 1] < hi - 1) & (q[:, 2] >= lo) & (q[:, 2] < hi - 1) & \
                ((((q[:, 2] - n + 0.5) * dx) < cut_x) | ((z0 + (q[:, 0] + 0.5) * dx) < top_low))
        else:
            keep = lambda q: (q[:, 1] >= lo) & (q[:, 1] < hi - 1) & (q[:, 2] >= lo) & (((q[:, 2] - n + 0.5) * dx) < cut_x)
        got = 0
        for seg in _segments(lines, keep, 1.5 / dx, [2]):
            i = tuple(np.rint(seg).astype(int).T)
            P.append(to_mm(seg)); Sp.append(np.sqrt(sum(np.asarray(F[c])[i].astype(float) ** 2 for c in F)) / U_D); got += 1
            if got >= count:
                break
    allp = np.concatenate(allp)
    check = dict(points=int(len(allp)), in_solid=int(points_in_solid(allp, S).sum()), stencil_solid=int(stencil_touches_solid(allp, S).sum()), lines=len(P))
    return _save_lines(f"woven_{state}", P, Sp, vmax), check


WOVEN_CAMS = {"A": dict(az=-24, el=17, dist=12.5, lens=45, target=(4.6, 3.6, 0.0)),      # low, looking along the channels (flow towards the viewer's left)
              "B": dict(az=-64, el=27, dist=10.2, lens=45, target=(4.3, 4.0, 0.0)),      # low front view: cut tow sections, step between the plies
              "A0": dict(az=-70, el=28, dist=19.5, lens=55, target=(4.0, 4.0, 0.0)),
              "B0": dict(az=-35, el=31, dist=19.5, lens=55, target=(4.0, 4.0, 0.0))}


# ---------------------------------------------------------------- foam
def foam_scene(tag="t0_full_buffer_crop200_tau1.5", variant="A", n_lines=230, vmax=3.0, az=-60.0):
    import pyvista as pv
    from scipy import ndimage
    d = np.load(FIELDS / "foam" / f"{tag}.npz"); solid = d["solid"]; dx = float(d["dx"]) * 1e3; buf = int(d["buffer"])
    nz = solid.shape[0]; cy, cx = (float(v) for v in d["centre"]); R = float(d["radius"]); zs = slice(buf, nz - buf)
    yy, xx = np.mgrid[:solid.shape[1], :solid.shape[2]]; ry, rx = yy - cy, xx - cx
    ins = np.hypot(ry, rx) <= R
    U_D = float(np.asarray(d["uz"][zs], dtype=np.float64)[:, ins].mean())
    if variant in ("A", "C"):    # the whole simulated section, a 90-degree sector towards the camera left out (C: upper half only)
        rel = (np.degrees(np.arctan2(ry, rx)) - az + 180) % 360 - 180
        open2d = ins & (np.abs(rel) <= 45); keep2d = ins & ~open2d; region2d = ins
    else:                 # an 18 mm block from the centre, the quadrant towards the camera left out
        h = 9.0 / dx; blk = (np.abs(ry) <= h) & (np.abs(rx) <= h)
        ux_, uy_ = math.cos(math.radians(az)), math.sin(math.radians(az))
        open2d = blk & (rx * np.sign(ux_) > 0) & (ry * np.sign(uy_) > 0); keep2d = blk & ~open2d; region2d = blk
    strut = solid[zs] & keep2d[None]
    if variant == "C":   # stepped cutaway: the sector keeps its struts below z_step
        z_step = int(0.45 * strut.shape[0]); strut[:z_step] |= solid[zs][:z_step] & open2d[None]
    f = ndimage.gaussian_filter(np.pad(strut, 1).astype(np.float32), 0.6); del strut
    g = pv.ImageData(dimensions=f.shape[::-1], spacing=(dx,) * 3, origin=(-(cx + 1) * dx, -(cy + 1) * dx, -dx)); g["f"] = f.ravel(); del f
    m = g.contour([0.5], scalars="f"); del g
    m = m.decimate(0.5) if m.n_cells > 8e6 else m
    m = m.smooth_taubin(n_iter=20, pass_band=0.08).compute_normals()
    CACHE.mkdir(parents=True, exist_ok=True); mfile = CACHE / f"foam_{variant}.ply"; m.save(mfile); print("foam mesh", m.n_cells)
    rng = np.random.default_rng(5); cand = np.argwhere(~solid[zs] & open2d[None]); cand[:, 0] += buf
    k = min(10 * n_lines, len(cand)); seeds = cand[rng.choice(len(cand), k, replace=False)] + rng.uniform(0.2, 0.8, (k, 3))
    F = {c: d[c] for c in ("ux", "uy", "uz")}
    lines, _ = trace_both_ways(F, solid, seeds, step=0.5, max_steps=5000)
    allp = np.concatenate([l for l in lines if len(l)])
    check = dict(points=int(len(allp)), in_solid=int(points_in_solid(allp, solid).sum()), stencil_solid=int(stencil_touches_solid(allp, solid).sum()))

    def keep(q):
        iy = np.clip(np.rint(q[:, 1]).astype(int), 0, solid.shape[1] - 1); ix = np.clip(np.rint(q[:, 2]).astype(int), 0, solid.shape[2] - 1)
        return (q[:, 0] >= buf) & (q[:, 0] <= nz - buf - 1) & open2d[iy, ix]
    P, Sp = [], []
    for seg in _segments(lines, keep, 5.0 / dx, [0]):
        i = tuple(np.rint(seg).astype(int).T)
        P.append(np.c_[(seg[:, 2] - cx) * dx, (seg[:, 1] - cy) * dx, (seg[:, 0] - buf) * dx])
        Sp.append(np.sqrt(sum(np.asarray(F[c])[i].astype(float) ** 2 for c in F)) / U_D)
        if len(P) >= n_lines:
            break
    check["lines"] = len(P)
    L = (nz - 2 * buf) * dx; half = R * dx if variant in ("A", "C") else 9.0
    return str(mfile), _save_lines(f"foam_{variant}", P, Sp, vmax), check, L, half


def camera(target, az, el, dist, lens):
    a, e = math.radians(az), math.radians(el)
    loc = [target[0] + dist * math.cos(e) * math.cos(a), target[1] + dist * math.cos(e) * math.sin(a), target[2] + dist * math.sin(e)]
    return {"loc": loc, "target": list(target), "lens": lens}


def look(spec, res, samples):
    spec.update(res=res, samples=samples, world=[1, 1, 1], world_strength=0.5, exposure=-1.3,
                view_transform="Standard", look="None", lights={"key": 0.6, "fill": 0.15, "rim": 0.2})
    return spec


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("variant"); ap.add_argument("outdir")
    ap.add_argument("--res", nargs=2, type=int, default=[1000, 750]); ap.add_argument("--samples", type=int, default=48)
    ap.add_argument("--reuse", action="store_true", help="reuse cached meshes and lines")
    ap.add_argument("--lines", type=int, default=170)
    a = ap.parse_args(); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    if a.what == "woven":
        from geometry import STATES, stack_layout
        wv = STATES["compacted"]; lay = stack_layout(wv)
        if a.reuse and (CACHE / "woven_compacted_lines.npz").exists():
            meshes = [(str(CACHE / f"woven_compacted_p{p}_{k}.ply"), k, p) for p in (0, 1) for k in ("warp", "weft")]
            lines, check = str(CACHE / "woven_compacted_lines.npz"), json.load(open(CACHE / "woven_check.json"))
        else:
            meshes, wv, lay = woven_meshes("compacted"); lines, check = woven_lines("compacted", wv, lay)
            (CACHE / "woven_check.json").write_text(json.dumps(check))
        c = WOVEN_CAMS[a.variant]
        spec = {"out": str(out / f"woven_{a.variant}_rgba.png"),
                "meshes": [{"ply": p, "color": GREY_WARP if k == "warp" else GREY_WEFT, "rough": 0.7, "specular": 0.25} for p, k, _ in meshes],
                "curves": [{"npz": lines, "radius": 0.011, "rough": 0.45, "emit": 0.25, "specular": 0.3}],
                "ground": {"z": lay["mould_bottom_mm"] - 0.003, "size": 300, "x": 4, "y": 4},
                "camera": camera(c["target"], c["az"], c["el"], c["dist"], c["lens"]), "tracer_check": check}
    else:
        if a.reuse and (CACHE / f"foam_{a.variant}_lines.npz").exists():
            mfile, lines = str(CACHE / f"foam_{a.variant}.ply"), str(CACHE / f"foam_{a.variant}_lines.npz")
            check, L, half = json.load(open(CACHE / f"foam_{a.variant}_check.json")).values()
        else:
            mfile, lines, check, L, half = foam_scene(variant=a.variant, n_lines=a.lines)
            (CACHE / f"foam_{a.variant}_check.json").write_text(json.dumps(dict(check=check, L=L, half=half)))
        cam = camera((0, 0, L * 0.30), -60.0, 34 if a.variant != "B" else 28, 2.75 * half + 1.6 * L if a.variant != "B" else 4.2 * half + 1.6 * L, 55)
        spec = {"out": str(out / f"foam_{a.variant}_rgba.png"),
                "meshes": [{"ply": mfile, "color": GREY_FOAM, "rough": 0.7, "specular": 0.25}],
                "curves": [{"npz": lines, "radius": 0.055, "rough": 0.45, "emit": 0.25, "specular": 0.3}],
                "ground": {"z": -0.01, "size": 600}, "camera": cam, "tracer_check": check}
        spec = look(spec, a.res, a.samples); spec.update(world_strength=0.75, exposure=-1.4, lights={"key": 0.42, "fill": 0.15, "rim": 0.2})
    spec = spec if "view_transform" in spec else look(spec, a.res, a.samples)
    f = out / f"{a.what}_{a.variant}.json"; f.write_text(json.dumps(spec, indent=1)); print(f, json.dumps(spec["tracer_check"]))

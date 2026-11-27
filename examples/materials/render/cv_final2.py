"""Blender scene specs for the gallery images of the two materials examples (second design).

    viz=python        # an environment with numpy, scipy, pyvista and matplotlib
    $viz examples/materials/render/cv_final2.py foam_main|foam_detail|woven_meshes|woven_fabric|woven_geometry|woven_exploded OUT [--res W H --samples N]
    LD_LIBRARY_PATH=/usr/lib/wsl/lib blender -b --factory-startup -P examples/materials/render/bl_scene.py -- OUT/<name>.json
    python examples/materials/render/compose.py OUT/<name>_rgba.png OUT/<name>.png cv

Foam: the whole simulated section (52 mm across, 15 mm long, scan resolution), flow upwards. The struts
are drawn opaque up to 70 % of the length, fade out between 70 % and 90 % and are not drawn above
(display only: the flow was computed through the whole section). Streamlines are seeded on an even
grid in the outlet buffer and traced back through the foam; colour |u|/U_D.

Woven (relaxed state, 5 um field): the fabric with a carbon-like material; the tow architecture with
warp and weft coloured; an exploded view with the upper ply lifted and the computed |u|/U_D on the
plane between the plies (z = 0.107 mm) drawn between them (display only: computed with the plies nested).
"""
import argparse, json, math, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent; MAT = HERE.parent; ROOT = MAT.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(MAT)); sys.path.insert(0, str(MAT / "woven_preform"))
from common import trace_both_ways, points_in_solid, stencil_touches_solid  # noqa: E402
import cv_final as c1  # noqa: E402

FIELDS = ROOT / "results" / "materials"; CACHE = FIELDS / "render_cache" / "cv_final2"; CACHE.mkdir(parents=True, exist_ok=True)
FOAM_TAG = "t0_full_buffer_crop200_tau1.5"; FOAM_VMAX = 2.6; FOAM_CMAP = "turbo"; WOVEN_VMAX = 10.0


def cmap(v, name="viridis"):
    import matplotlib
    return matplotlib.colormaps[name](np.clip(v, 0, 1))[:, :3]


def save_lines(name, P, S, vmax, name_cmap="viridis"):
    f = CACHE / f"{name}_lines.npz"
    np.savez(f, pts=np.concatenate(P), sizes=np.array([len(p) for p in P]), col=np.concatenate([cmap(s / vmax, name_cmap) for s in S]))
    return str(f)


# ------------------------------------------------------------------ foam
def foam_load(tag=FOAM_TAG):
    d = np.load(FIELDS / "foam" / f"{tag}.npz"); solid = d["solid"]; dx = float(d["dx"]) * 1e3; buf = int(d["buffer"])
    nz = solid.shape[0]; cy, cx = (float(v) for v in d["centre"]); R = float(d["radius"])
    yy, xx = np.mgrid[:solid.shape[1], :solid.shape[2]]; ins = np.hypot(yy - cy, xx - cx) <= R
    U_D = float(np.asarray(d["uz"][buf:nz - buf], dtype=np.float64)[:, ins].mean())
    return d, solid, dx, buf, nz, cy, cx, R, ins, U_D


def strut_mesh(mask, dx, origin, name, decimate=None):
    import pyvista as pv
    from scipy import ndimage
    f = ndimage.gaussian_filter(np.pad(mask, 1).astype(np.float32), 0.6)
    g = pv.ImageData(dimensions=f.shape[::-1], spacing=(dx,) * 3, origin=tuple(o - dx for o in origin)); g["f"] = f.ravel(); del f
    m = g.contour([0.5], scalars="f"); del g
    if decimate:
        m = m.decimate(decimate)
    m = m.smooth_taubin(n_iter=20, pass_band=0.08).compute_normals()
    fn = CACHE / f"{name}.ply"; m.save(fn); print(name, m.n_cells, "triangles"); return str(fn)


def foam_main(spacing_mm=4.0, tag=FOAM_TAG, name="foam_full", decimate=0.4, min_frac=0.8):
    d, solid, dx, buf, nz, cy, cx, R, ins, U_D = foam_load(tag); L = (nz - 2 * buf) * dx
    mfile = strut_mesh(solid[buf:nz - buf] & ins[None], dx, (-cx * dx, -cy * dx, 0.0), name, decimate=decimate)
    # seeds: jittered hexagonal grid in the outlet buffer (all fluid there), traced back through the foam
    rng = np.random.default_rng(2); s = spacing_mm / dx; pts = []
    for j, y in enumerate(np.arange(-R, R, s * 0.866)):
        for x in np.arange(-R, R, s) + (0.5 * s if j % 2 else 0):
            if math.hypot(x, y) < R - 1.2 / dx:
                pts.append((nz - buf + 0.3 * buf + rng.uniform(0, 0.3 * buf), cy + y + rng.uniform(-0.2, 0.2) * s, cx + x + rng.uniform(-0.2, 0.2) * s))
    seeds = np.array(pts); F = {c: d[c] for c in ("ux", "uy", "uz")}
    lines, _ = trace_both_ways(F, solid, seeds, step=0.5, max_steps=6000)
    allp = np.concatenate([l for l in lines if len(l)])
    check = dict(points=int(len(allp)), in_solid=int(points_in_solid(allp, solid).sum()), stencil_solid=int(stencil_touches_solid(allp, solid).sum()), seeds=len(seeds))
    P, Sp = [], []
    for ln in lines:
        q = ln[(ln[:, 0] >= buf) & (ln[:, 0] <= nz - 2)]                    # from the foam's inlet face to the end of the outlet buffer
        if len(q) < 10 or np.ptp(q[:, 0]) * dx < min_frac * L:
            continue
        i = tuple(np.rint(q).astype(int).T)
        P.append(np.c_[(q[:, 2] - cx) * dx, (q[:, 1] - cy) * dx, (q[:, 0] - buf) * dx])
        Sp.append(np.sqrt(sum(np.asarray(F[c])[i].astype(float) ** 2 for c in F)) / U_D)
    check["lines"] = len(P)
    return mfile, save_lines(name, P, Sp, FOAM_VMAX, FOAM_CMAP), check, L, R * dx


def foam_detail(half_mm=5.0, n_lines=75):
    d, solid, dx, buf, nz, cy, cx, R, ins, U_D = foam_load()
    h = int(round(half_mm / dx)); zc = buf + 100; cyi, cxi = int(round(cy)) + 60, int(round(cx)) - 40
    K0, Y0, X0 = zc - h, cyi - h, cxi - h; n = 2 * h; S = solid[K0:K0 + n, Y0:Y0 + n, X0:X0 + n]
    mfile = strut_mesh(S, dx, (0.5 * dx,) * 3, "foam_detail")
    rng = np.random.default_rng(9); cand = np.argwhere(~S[n // 2 - 2:n // 2 + 2, 8:-8, 8:-8]) + [K0 + n // 2 - 2, Y0 + 8, X0 + 8]
    seeds = cand[rng.choice(len(cand), 8 * n_lines, replace=False)] + rng.uniform(0.2, 0.8, (8 * n_lines, 3))
    F = {c: d[c] for c in ("ux", "uy", "uz")}
    lines, _ = trace_both_ways(F, solid, seeds, step=0.4, max_steps=3000)
    allp = np.concatenate([l for l in lines if len(l)])
    check = dict(points=int(len(allp)), in_solid=int(points_in_solid(allp, solid).sum()), stencil_solid=int(stencil_touches_solid(allp, solid).sum()))
    P, Sp = [], []
    for ln in lines:
        q = ln - [K0, Y0, X0]; q = q[np.all((q >= 0) & (q < n - 1), axis=1)]
        if len(q) < 10 or np.ptp(q[:, 0]) < 0.9 * n:
            continue
        i = tuple(np.rint(q + [K0, Y0, X0]).astype(int).T)
        P.append(np.c_[(q[:, 2] + 0.5) * dx, (q[:, 1] + 0.5) * dx, (q[:, 0] + 0.5) * dx])
        Sp.append(np.sqrt(sum(np.asarray(F[c])[i].astype(float) ** 2 for c in F)) / U_D)
        if len(P) >= n_lines:
            break
    check["lines"] = len(P)
    return mfile, save_lines("foam_detail", P, Sp, FOAM_VMAX, FOAM_CMAP), check, n * dx


# ------------------------------------------------------------------ woven
def woven_setup():
    from geometry import STATES, stack_layout
    wv = STATES["relaxed"]; lay = stack_layout(wv); return wv, lay


def woven_meshes(cells, step, tag):
    """Analytic tow surfaces over cells x cells unit cells, closed at the sample edges: {(ply, kind): path}."""
    import pyvista as pv
    from geometry import tow_margins, corrected_inside
    wv, lay = woven_setup(); W = cells * wv.cell_mm
    xs = np.arange(-step, W + 1.5 * step, step); zs = np.arange(lay["mould_bottom_mm"] - 2 * step, lay["mould_top_mm"] + 2 * step, step)
    Z, Y, X = np.meshgrid(zs, xs, xs, indexing="ij"); outside = (X < 0) | (X > W) | (Y < 0) | (Y > W); out = {}
    for p, (z0, sh) in enumerate(zip(lay["ply_z0_mm"], lay["ply_shift"])):
        mg = tow_margins(X, Y, Z, wv, z0, sh); warp, weft, _ = corrected_inside(mg)
        for kind in ("warp", "weft"):
            m = np.where(np.isfinite(mg[kind][0]), mg[kind][0], -0.05); other = warp if kind == "weft" else weft
            f = np.where(other & (m > 0), -1e-3, m).astype(np.float32); f[outside] = -1.0
            g = pv.ImageData(dimensions=(len(xs), len(xs), len(zs)), spacing=(step,) * 3, origin=(xs[0], xs[0], zs[0])); g["f"] = f.ravel()
            mesh = g.contour([0.0], scalars="f").smooth_taubin(n_iter=20, pass_band=0.08).compute_normals(auto_orient_normals=True)
            fn = CACHE / f"woven_{tag}_p{p}_{kind}.ply"; mesh.save(fn); out[f"{p}_{kind}"] = str(fn); print(fn.name, mesh.n_cells)
    (CACHE / f"woven_{tag}.json").write_text(json.dumps(out)); return out


def woven_slice(cells=2, k=50, sub=2, lift=0.0, shift=(0.0, 0.0)):
    """|u|/U_D on the plane between the plies as a vertex-coloured mesh (resin only), plus direction arrows."""
    import pyvista as pv
    wv, lay = woven_setup(); dx = 0.005
    d = np.load(FIELDS / "woven" / "relaxed_mould_dx5um_x.npz"); nzg = d["solid"].shape[0]
    z0 = lay["mould_bottom_mm"] + (lay["cavity_mm"] - (nzg - 2) * dx) / 2 - dx; zk = z0 + (k + 0.5) * dx
    U_D = float(d["ux"][1:-1].astype(np.float64).mean())
    so = np.tile(d["solid"][k, ::sub, ::sub], (cells, cells))
    u = np.tile(d["ux"][k, ::sub, ::sub].astype(np.float32), (cells, cells)) / U_D; v = np.tile(d["uy"][k, ::sub, ::sub].astype(np.float32), (cells, cells)) / U_D
    w = np.tile(d["uz"][k, ::sub, ::sub].astype(np.float32), (cells, cells)) / U_D; sp = np.sqrt(u * u + v * v + w * w)
    h = dx * sub; ny, nx = so.shape
    g = pv.ImageData(dimensions=(nx, ny, 1), spacing=(h, h, 1), origin=(0.5 * h + shift[0], 0.5 * h + shift[1], zk + lift))
    g["rgb"] = (cmap(sp.ravel() / WOVEN_VMAX, "turbo") * 255).astype(np.uint8); g["s"] = so.ravel().astype(np.float32)
    m = g.threshold(0.5, scalars="s", invert=True, all_scalars=True).extract_surface().triangulate()
    fs = CACHE / "woven_slice.ply"; m.save(fs, texture="rgb")
    st = 26; yy, xx = np.mgrid[st // 2:ny:st, st // 2:nx:st]; ok = ~so[yy, xx] & (sp[yy, xx] > 0.8)
    pts = np.c_[(xx[ok] + 0.5) * h + shift[0], (yy[ok] + 0.5) * h + shift[1], np.full(ok.sum(), zk + lift + 0.014)]
    a = pv.PolyData(pts); vec = np.c_[u[yy, xx][ok], v[yy, xx][ok], np.zeros(ok.sum())]; a["v"] = vec / np.linalg.norm(vec, axis=1)[:, None]
    a["len"] = 0.14 + 0.22 * np.sqrt(np.clip(sp[yy, xx][ok] / WOVEN_VMAX, 0, 1))          # arrow length (mm) grows with speed
    gl = a.glyph(orient="v", scale="len", factor=1.0, geom=pv.Arrow(tip_length=0.45, tip_radius=0.17, shaft_radius=0.06, tip_resolution=14, shaft_resolution=10))
    gl = gl.extract_surface().triangulate(); gl["rgb"] = np.full((gl.n_points, 3), 255, np.uint8)
    fa = CACHE / "woven_arrows.ply"; gl.save(fa, texture="rgb")
    return str(fs), str(fa), zk, int(ok.sum())


def woven_plane_lines(cells=2, k=50, sub=2, dsep_mm=0.24, min_len_mm=0.7):
    """Evenly spaced streamlines of the in-plane velocity (u_x, u_y) on the plane between the plies, as white tubes
    with one arrowhead each. They show direction on that plane; the out-of-plane component is not drawn."""
    import pyvista as pv
    wv, lay = woven_setup(); dx = 0.005
    d = np.load(FIELDS / "woven" / "relaxed_mould_dx5um_x.npz"); nzg = d["solid"].shape[0]
    z0 = lay["mould_bottom_mm"] + (lay["cavity_mm"] - (nzg - 2) * dx) / 2 - dx; zk = z0 + (k + 0.5) * dx
    U_D = float(d["ux"][1:-1].astype(np.float64).mean())
    so = np.tile(d["solid"][k, ::sub, ::sub], (cells, cells))
    u = np.tile(d["ux"][k, ::sub, ::sub].astype(np.float64), (cells, cells)) / U_D; v = np.tile(d["uy"][k, ::sub, ::sub].astype(np.float64), (cells, cells)) / U_D
    h = dx * sub; ny, nx = so.shape; ds = dsep_mm / h; occ = np.zeros((int(ny / ds) + 1, int(nx / ds) + 1), bool)

    def vel(p):
        j, i = int(round(p[0])), int(round(p[1]))
        if j < 1 or i < 1 or j >= ny - 1 or i >= nx - 1 or so[j, i]:
            return None
        w = np.array([v[j, i], u[j, i]]); n = np.hypot(*w)
        return None if n < 0.25 else w / n

    def trace(p0, sign):
        pts = []; p = np.array(p0, float); own = set()
        for _ in range(6000):
            w = vel(p)
            if w is None:
                break
            p = p + sign * 0.5 * w; c = (int(p[0] / ds), int(p[1] / ds))
            if c[0] < 0 or c[1] < 0 or c[0] >= occ.shape[0] or c[1] >= occ.shape[1] or (occ[c] and c not in own):
                break
            own.add(c); pts.append(p.copy())
        return pts, own
    rng = np.random.default_rng(1); sp = np.hypot(u, v); lines = []
    cand = np.argwhere(~so & (sp > 0.6)); cand = cand[rng.permutation(len(cand))[:40000]]
    cand = cand[np.argsort(-sp[cand[:, 0], cand[:, 1]])]                    # fast channels get their lines first
    for c in cand:
        cc = (int(c[0] / ds), int(c[1] / ds))
        if occ[cc]:
            continue
        fw, o1 = trace(c, 1.0); bw, o2 = trace(c, -1.0); ln = np.array(bw[::-1] + [c.astype(float)] + fw)
        if len(ln) * 0.5 * h < min_len_mm:
            continue
        for q in o1 | o2 | {cc}:
            occ[q] = True
        lines.append(ln)
    tubes = pv.PolyData(); cones = pv.PolyData(); zl = zk + 0.014
    for ln in lines:
        P = np.c_[(ln[::4, 1] + 0.5) * h, (ln[::4, 0] + 0.5) * h, np.full(len(ln[::4]), zl)]
        if len(P) < 3:
            continue
        tubes += pv.Spline(P, max(len(P) * 2, 8)).tube(radius=0.011, n_sides=8)
        m = int(len(P) * 0.6); t = P[min(m + 1, len(P) - 1)] - P[m - 1]; t /= np.linalg.norm(t)
        cones += pv.Cone(center=P[m], direction=t, height=0.15, radius=0.05, resolution=14)
    out = []
    for name, mesh in (("woven_plane_lines", tubes), ("woven_plane_heads", cones)):
        mesh = mesh.triangulate(); mesh["rgb"] = np.full((mesh.n_points, 3), 255, np.uint8); fn = CACHE / f"{name}.ply"; mesh.save(fn, texture="rgb"); out.append(str(fn))
    print("in-plane lines", len(lines)); return out, len(lines)


def look(spec, res, samples, exposure=-1.3, world=0.5, lights=(0.6, 0.15, 0.2)):
    spec.update(res=res, samples=samples, world=[1, 1, 1], world_strength=world, exposure=exposure, view_transform="Standard", look="None",
                lights={"key": lights[0], "fill": lights[1], "rim": lights[2]})
    return spec


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("outdir")
    ap.add_argument("--res", nargs=2, type=int, default=[1000, 750]); ap.add_argument("--samples", type=int, default=48)
    ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args(); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True); w = a.what
    GREY = [0.55, 0.52, 0.47]
    if w == "foam_main":
        ck = CACHE / "foam_main.json"
        if a.reuse and ck.exists():
            mfile, lines, check, L, Rm = json.load(open(ck))
        else:
            mfile, lines, check, L, Rm = foam_main(); ck.write_text(json.dumps([mfile, lines, check, L, Rm]))
        spec = {"out": str(out / "foam_main_rgba.png"),
                "meshes": [{"ply": mfile, "color": GREY, "rough": 0.7, "specular": 0.25, "fade": [0.60 * L, 0.86 * L]}],
                "curves": [{"npz": lines, "radius": 0.085, "rough": 0.45, "emit": 0.3, "specular": 0.3}],
                "ground": {"z": -0.01, "size": 600},
                "camera": c1.camera((0, 0, L * 0.52), -62.0, 21.0, 2.75 * Rm + 1.6 * L, 55), "tracer_check": check}
        spec = look(spec, a.res, a.samples, exposure=-1.4, world=0.75, lights=(0.42, 0.15, 0.2))
    elif w == "foam_whole":                                 # the entire scanned sample (44 mm long), 2x coarsened grid as simulated
        ck = CACHE / "foam_whole.json"
        if a.reuse and ck.exists():
            mfile, lines, check, L, Rm = json.load(open(ck))
        else:
            mfile, lines, check, L, Rm = foam_main(spacing_mm=2.1, tag="t0_coarse_porosity_buffer_tau1.5", name="foam_whole", decimate=None, min_frac=0.38); ck.write_text(json.dumps([mfile, lines, check, L, Rm]))
        spec = {"out": str(out / "foam_whole_rgba.png"),
                "meshes": [{"ply": mfile, "color": GREY, "rough": 0.7, "specular": 0.25, "fade": [0.62 * L, 0.86 * L]}],
                "curves": [{"npz": lines, "radius": 0.105, "rough": 0.45, "emit": 0.3, "specular": 0.3}],
                "ground": {"z": -0.01, "size": 900},
                "camera": c1.camera((0, 0, L * 0.50), -62.0, 19.0, 156.0, 55), "tracer_check": check}
        spec = look(spec, a.res, a.samples, exposure=-1.4, world=0.75, lights=(0.42, 0.15, 0.2))
    elif w == "foam_detail":
        ck = CACHE / "foam_detail.json"
        if a.reuse and ck.exists():
            mfile, lines, check, E = json.load(open(ck))
        else:
            mfile, lines, check, E = foam_detail(); ck.write_text(json.dumps([mfile, lines, check, E]))
        spec = {"out": str(out / "foam_detail_rgba.png"),
                "meshes": [{"ply": mfile, "color": GREY, "rough": 0.7, "specular": 0.25}],
                "curves": [{"npz": lines, "radius": 0.05, "rough": 0.45, "emit": 0.3, "specular": 0.3}],
                "camera": c1.camera((E / 2, E / 2, E / 2), -58.0, 17.0, 3.0 * E, 55), "tracer_check": check}
        spec = look(spec, a.res, a.samples, exposure=-1.4, world=0.75, lights=(0.42, 0.15, 0.2))
    elif w == "woven_meshes":
        woven_meshes(3, 0.008, "c3"); woven_meshes(2, 0.008, "c2"); sys.exit()
    else:
        wv, lay = woven_setup(); L = wv.cell_mm
        if w == "woven_fabric":
            M = json.load(open(CACHE / "woven_c3.json")); W = 3 * L
            spec = {"out": str(out / "woven_fabric_rgba.png"),
                    "meshes": [{"ply": f, "color": [0.035, 0.037, 0.042], "rough": 0.38, "specular": 0.6,
                                "fibre": {"across": "Y" if k.endswith("warp") else "X", "scale": 60, "strength": 0.7, "aniso": 0.6, "distortion": 1.5}}
                               for k, f in M.items()],
                    "ground": {"z": lay["mould_bottom_mm"] - 0.004, "size": 400, "x": W / 2, "y": W / 2},
                    "camera": c1.camera((W / 2, W / 2 - 0.3, 0.0), -62.0, 33.0, 33.0, 55)}
            spec = look(spec, a.res, a.samples, exposure=-0.9, world=0.9, lights=(1.3, 0.35, 0.6))
        elif w == "woven_geometry":
            M = json.load(open(CACHE / "woven_c2.json")); W = 2 * L
            spec = {"out": str(out / "woven_geometry_rgba.png"),
                    "meshes": [{"ply": f, "color": [0.22, 0.38, 0.58] if k.endswith("warp") else [0.80, 0.62, 0.34], "rough": 0.7, "specular": 0.25} for k, f in M.items()],
                    "ground": {"z": lay["mould_bottom_mm"] - 0.004, "size": 400, "x": W / 2, "y": W / 2},
                    "camera": c1.camera((W / 2, W / 2 - 0.2, 0.0), -62.0, 26.0, 22.5, 55)}
            spec = look(spec, a.res, a.samples)
        elif w == "woven_exploded":
            # lower ply with the computed flow map on it; the upper ply lifted and set aside (+x) for display
            M = json.load(open(CACHE / "woven_c2.json")); W = 2 * L; lift = 1.6; side = 9.4
            fs, fa, zk, na = woven_slice(2)
            meshes = [{"ply": f, "color": [0.60, 0.62, 0.66] if k.endswith("warp") else [0.80, 0.81, 0.83], "rough": 0.75, "specular": 0.2,
                       "offset": [side, 0, lift] if k.startswith("1") else [0, 0, 0], "no_shadow": k.startswith("1")} for k, f in M.items()]
            meshes.append({"ply": fs, "vertex_colors": True, "emit": 2 ** 1.3, "smooth": False})
            meshes.append({"ply": fa, "color": [0.03, 0.03, 0.04], "rough": 0.6})
            zb = lay["mould_bottom_mm"]
            bars = [{"start": [x, y, 0.12], "end": [x + side, y, lift], "radius": 0.012, "name": f"g{i}"} for i, (x, y) in enumerate(((W, 0), (W, W)))]
            spec = {"out": str(out / "woven_exploded_rgba.png"), "meshes": meshes, "bars": bars,
                    "ground": {"z": zb - 0.004, "size": 400, "x": W / 2, "y": W / 2},
                    "camera": c1.camera(((W + side) / 2 + 0.25, W / 2 - 0.4, 0.3), -82.0, 43.0, 35.5, 55), "n_arrows": na, "slice_z_mm": zk}
            spec = look(spec, a.res, a.samples)
        elif w == "woven_flow":
            # the lower ply with the computed velocity map and in-plane flow lines on the plane between the plies;
            # the upper ply lifted and set back, for display only
            M = json.load(open(CACHE / "woven_c2.json")); W = 2 * L; lift = 2.2; back = 5.4
            fs, fa, zk, nl = woven_slice(2)
            meshes = [{"ply": f, "color": [0.60, 0.62, 0.66] if k.endswith("warp") else [0.80, 0.81, 0.83], "rough": 0.75, "specular": 0.2,
                       "offset": [0, back, lift] if k.startswith("1") else [0, 0, 0], "no_shadow": k.startswith("1")} for k, f in M.items()]
            meshes.append({"ply": fs, "vertex_colors": True, "emit": 2 ** 1.3, "smooth": False})
            meshes.append({"ply": fa, "vertex_colors": True, "emit": 2 ** 1.3})
            zb = lay["mould_bottom_mm"]
            bars = [{"start": [x, y, 0.12], "end": [x, y + back, lift], "radius": 0.012, "name": f"g{i}"} for i, (x, y) in enumerate(((0, 0), (W, 0), (W, W), (0, W)))]
            spec = {"out": str(out / "woven_flow_rgba.png"), "meshes": meshes, "bars": bars,
                    "ground": {"z": zb - 0.004, "size": 400, "x": W / 2, "y": W / 2},
                    "camera": c1.camera((W / 2, W / 2 + back * 0.46, lift * 0.3), -70.0, 47.0, 30.0, 55), "n_arrows": nl, "slice_z_mm": zk}
            spec = look(spec, a.res, a.samples)
    f = out / f"{w}.json"; f.write_text(json.dumps(spec, indent=1)); print(f, json.dumps(spec.get("tracer_check", {})))

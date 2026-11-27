"""Figures for the woven-preform example (docs/figures/materials).

    python examples/materials/woven_preform/figures.py

    woven_resolution.png   k against voxel size for both states and both stackings, with the
                           observed order of convergence and a Richardson estimate
    woven_section.png      u_x / U_D in a y-z section through the 5 um field, tow outlines
"""
import glob, json, math, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
OUT = ROOT / "docs" / "figures" / "materials"; RES = HERE / "results"; FIELDS = ROOT / "results" / "materials" / "woven"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def load_series(state, stack):
    """(voxel um, k_xx, fibre Vf) at tau 1.5 from summary.json (cavity-normalised permeabilities)."""
    S = json.load(open(HERE.parent / "summary.json"))["woven"]; rows = []
    for v in S.values():
        if v["state"] == state and v["stack"] == stack and abs(v["tau"] - 1.5) < 1e-9 and "x" in v["loads"]:
            rows.append((v["dx_um"], v["loads"]["x"]["k_m2"], v["fibre_volume_fraction"]))
    return sorted(rows)


def richardson(rows):
    """Observed order and extrapolated k from three voxel sizes halving each time."""
    if len(rows) < 3:
        return None
    (h1, k1, _), (h2, k2, _), (h3, k3, _) = rows[:3]                        # finest first
    if not (abs(h2 / h1 - 2) < 1e-6 and abs(h3 / h2 - 2) < 1e-6) or (k2 - k1) == 0:
        return None
    r = (k3 - k2) / (k2 - k1)
    if r <= 0:
        return dict(order=None, k_extrapolated=None, ratio=r)
    p = math.log(r) / math.log(2)
    return dict(order=p, k_extrapolated=k1 + (k1 - k2) / (2 ** p - 1), ratio=r)


def resolution_figure():
    fig, ax = plt.subplots(figsize=(6.2, 3.8)); summary = {}
    for state, colour in (("relaxed", "#3b6ea5"), ("compacted", "#c2571a")):
        for stack, ls, mk in (("mould", "-", "o"), ("periodic", "--", "s")):
            rows = load_series(state, stack)
            if not rows:
                continue
            h = [r[0] for r in rows]; k = [r[1] for r in rows]
            vf = rows[0][2]
            ax.plot(h, k, ls, marker=mk, color=colour, label=f"{state}, {'2 plies in mould' if stack == 'mould' else 'infinite laminate'} (Vf {vf:.2f})")
            rc = richardson(rows); summary[f"{state}_{stack}"] = dict(rows=rows, richardson=rc)
            if rc and rc["k_extrapolated"]:
                ax.plot([0], [rc["k_extrapolated"]], marker=mk, mfc="white", color=colour)
    ax.set_xlabel("voxel size (µm)"); ax.set_ylabel("in-plane permeability $k_{xx}$ (m²)"); ax.set_yscale("log"); ax.set_xlim(-1, 22)
    ax.set_title("Resolution study (open symbols at 0: Richardson estimate)", fontsize=9); ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(); OUT.mkdir(parents=True, exist_ok=True); fig.savefig(OUT / "woven_resolution.png"); plt.close(fig)
    (RES / "resolution_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    return summary


def section_figure(state="compacted"):
    f = FIELDS / f"{state}_mould_dx5um_x.npz"
    if not f.exists():
        return
    d = np.load(f); solid, labels, ux = d["solid"], d["labels"], d["ux"].astype(np.float64)
    U_D = ux[1:-1].mean(); dx = 5e-3
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.2), sharex=True)
    for ax, xi, title in ((axs[0], 0, "section at a weft crossing (x = 0)"), (axs[1], solid.shape[2] // 4, "section between crossings (x = s/2)")):
        sl = ux[:, :, xi] / U_D; m = np.ma.masked_where(solid[:, :, xi], sl)
        im = ax.imshow(m, origin="lower", cmap="viridis", vmin=0, vmax=6, extent=[0, solid.shape[1] * dx, 0, solid.shape[0] * dx], aspect="auto", interpolation="nearest")
        lab = labels[:, :, xi].astype(float)
        ax.imshow(np.ma.masked_where(~solid[:, :, xi], lab), origin="lower", cmap="Greys", vmin=-2, vmax=6, alpha=0.55,
                  extent=[0, solid.shape[1] * dx, 0, solid.shape[0] * dx], aspect="auto", interpolation="nearest")
        ax.set_title(title, fontsize=8, loc="left"); ax.set_ylabel("z (mm)")
    axs[1].set_xlabel("y (mm)")
    cb = fig.colorbar(im, ax=axs, shrink=0.85); cb.set_label("$u_x / U_D$  (flow into the page)")
    fig.suptitle(f"Where the resin flows: {state} state, 5 µm field, load along x (warp)", fontsize=9)
    fig.savefig(OUT / f"woven_section_{state}.png", bbox_inches="tight"); plt.close(fig)


def geometry_figure(dx=0.005):
    """Sections through the voxel model of both states: along a warp tow (crimp, nesting) and across the tows."""
    sys.path.insert(0, str(HERE))
    from geometry import STATES, voxelise
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#ffffff", "#5b6f8c", "#9aa9c0", "#8c6a5b", "#c0a99a"] + ["#ffffff"] * 250 + ["#333333"])
    fig, axs = plt.subplots(2, 2, figsize=(9.5, 3.6), sharex="col")
    for row, name in enumerate(("relaxed", "compacted")):
        wv = STATES[name]; solid, lab, info = voxelise(wv, dx)
        nz = lab.shape[0]; ext_z = nz * dx
        sec_xz = lab[:, 0, :]                                         # along x through warp tow 0 (y = 0)
        sec_yz = lab[:, :, lab.shape[2] // 4]                         # across the tows, between crossings
        for col, (sec, xl) in enumerate(((sec_xz, "x (mm), along a warp tow"), (sec_yz, "y (mm), across the warp tows"))):
            ax = axs[row, col]
            ax.imshow(sec, origin="lower", cmap=cmap, vmin=0, vmax=255, interpolation="nearest",
                      extent=[0, sec.shape[1] * dx, 0, ext_z], aspect=2.0)
            ax.set_ylabel("z (mm)" if col == 0 else ""); ax.set_title(f"{name}: fibre Vf {info['fibre_volume_fraction']:.2f}, tow h {wv.thickness_mm*1e3:.0f} µm, gap {wv.spacing_mm*1e3-wv.width_mm*1e3:.0f} µm" if col == 0 else "", fontsize=8, loc="left")
            if row == 1:
                ax.set_xlabel(xl)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.text(0.01, 0.01, "Blue: ply 1 (dark warp, light weft); brown: ply 2; black: mould plates; white: resin flow space. Vertical scale x2.", fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(OUT / "woven_geometry.png"); plt.close(fig)


if __name__ == "__main__":
    print(json.dumps(resolution_figure(), indent=1, default=float)[:2000])
    for s in ("relaxed", "compacted"):
        section_figure(s)
    geometry_figure()

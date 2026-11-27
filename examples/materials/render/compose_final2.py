"""Compose the two materials gallery figures from the Blender renders (cv_final2.py).

    python examples/materials/render/compose_final2.py RENDER_DIR OUT_DIR
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image

R, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
INK, MUTED = "#1f2328", "#59616b"


def img(name):
    return np.asarray(Image.open(R / f"{name}.png").convert("RGB"))


def crop(a, pad=40, thr=247):
    m = (a < thr).any(axis=2); r = np.flatnonzero(m.any(axis=1)); c = np.flatnonzero(m.any(axis=0))
    return a[max(r[0] - pad, 0):r[-1] + pad, max(c[0] - pad, 0):c[-1] + pad]


def cbar(fig, rect, vmax, cmap, label, ticks):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(0, vmax), cmap), cax=cax, orientation="horizontal", ticks=ticks)
    cb.outline.set_linewidth(0.5); cb.ax.tick_params(labelsize=8.5, width=0.5, length=2.5, pad=2, colors=INK)
    cb.ax.set_title(label, fontsize=9.5, pad=4, color=INK)


def foam():
    fig = plt.figure(figsize=(9.333, 7.0), dpi=300); fig.patch.set_facecolor("white")
    m = img("foam_main").astype(np.float32); hh, ww = m.shape[:2]; e = 160                 # feather the render's edges into the white page
    ry = np.clip(np.minimum(np.arange(hh), hh - 1 - np.arange(hh)) / e, 0, 1)[:, None]; rx = np.clip(np.minimum(np.arange(ww), ww - 1 - np.arange(ww)) / e, 0, 1)[None, :]
    wgt = (np.minimum(ry, rx) ** 0.8)[..., None]; m = (255 - (255 - m) * wgt).astype(np.uint8)
    ax = fig.add_axes([0.0, 0.065, 0.86, 0.86]); ax.imshow(m); ax.set_axis_off()
    fig.text(0.035, 0.962, "Open-cell polyurethane foam: pore-scale flow", fontsize=13, color=INK, weight="bold", va="top")
    fig.text(0.035, 0.925, "micro-CT geometry (74.8 µm voxels, 52 mm sample) · lattice-Boltzmann solution, 88.7 M fluid cells",
             fontsize=8.4, color=MUTED, va="top")
    cbar(fig, [0.035, 0.825, 0.20, 0.015], 2.6, "turbo", "local speed / mean (Darcy) speed", [0, 1, 2])
    d = crop(img("foam_detail"), 25); x0, y0, w, h = 0.745, 0.60, 0.235, 0.315          # detail inset, top right
    ia = fig.add_axes([x0, y0, w, h]); ia.imshow(d); ia.set_xticks([]); ia.set_yticks([]); ia.set_anchor("C")
    ia.set_facecolor("white"); ia.set_aspect("equal", adjustable="datalim")
    for s in ia.spines.values():
        s.set_linewidth(0.8); s.set_color("#8a9099")
    fig.text(x0, y0 + h + 0.008, "detail, 10 mm: flow through cell windows", fontsize=7.4, color=INK, va="bottom")
    ax.add_patch(FancyArrowPatch((0.035, 0.30), (0.035, 0.46), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=16, lw=1.6, color=INK))
    fig.text(0.035, 0.475, "flow", fontsize=10, color=INK, va="bottom", ha="center")
    fig.text(0.035, 0.02, "Struts are faded out over the last part of the sample for display only; the flow was computed through the whole sample.\n"
             "Scan: Ruf & Steeb (2024), DaRUS, doi:10.18419/DARUS-3010, CC BY 4.0.", fontsize=7.2, color=MUTED, va="bottom", linespacing=1.5)
    fig.savefig(OUT / "polyurethane_foam_cv.png", dpi=300); plt.close(fig); print("saved foam")


def woven():
    fig = plt.figure(figsize=(9.333, 7.0), dpi=300); fig.patch.set_facecolor("white")
    fig.text(0.03, 0.972, "Woven carbon preform: from fabric to resin-flow simulation", fontsize=13, color=INK, weight="bold", va="top")
    tops = (("woven_fabric", 0.03, "1   Fabric: two plies of 3K plain weave (idealised model)"),
            ("woven_geometry", 0.515, "2   Tow architecture: warp and weft, two nested plies"))
    for name, x0, title in tops:
        a = crop(img(name), 30); ax = fig.add_axes([x0, 0.615, 0.455, 0.27]); ax.imshow(a); ax.set_axis_off(); ax.set_anchor("C")
        fig.text(x0, 0.905, title, fontsize=9.8, color=INK, va="bottom")
    for k, (col, lab) in enumerate((("#39608f", "warp (flow direction)"), ("#c9a05c", "weft"))):
        fig.patches.append(Rectangle((0.84, 0.652 - 0.027 * k), 0.016, 0.014, transform=fig.transFigure, color=col))
        fig.text(0.862, 0.659 - 0.027 * k, lab, fontsize=7.8, color=INK, va="center")
    a = crop(img("woven_exploded"), 30); ax = fig.add_axes([0.03, 0.10, 0.94, 0.42]); ax.imshow(a); ax.set_axis_off(); ax.set_anchor("C")
    fig.text(0.03, 0.555, "3   CFD: computed resin velocity in the channels between the tows (upper ply lifted aside)", fontsize=9.8, color=INK, va="bottom")
    fig.text(0.115, 0.085, "lower ply with the velocity field on the plane between the plies", fontsize=8.2, color=INK)
    fig.text(0.66, 0.085, "upper ply, lifted for display", fontsize=8.2, color=INK)
    cbar(fig, [0.40, 0.505, 0.20, 0.014], 12, "turbo", "local speed / mean (Darcy) speed", [0, 4, 8, 12])
    fig.text(0.03, 0.022, "Lattice-Boltzmann solution on a 5 µm grid (800 × 800 × 105 cells per unit cell), flow along the warp.\n"
             "Tows are modelled as impermeable; the flow was computed with the plies nested as in view 2.", fontsize=7.2, color=MUTED, va="bottom", linespacing=1.5)
    fig.savefig(OUT / "woven_composite_cv.png", dpi=300); plt.close(fig); print("saved woven")


if __name__ == "__main__":
    foam(); woven()

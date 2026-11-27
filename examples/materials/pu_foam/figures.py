"""Figures for the PU-foam example (docs/figures/materials).

    python examples/materials/pu_foam/figures.py

    foam_compression.png    k (per sample area) and porosity against axial strain, all states
    foam_sensitivity.png    t0: resolution, coarsening rule, end treatment and tau
"""
import json, glob, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]; RES = HERE / "results"; OUT = ROOT / "docs" / "figures" / "materials"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def rec(tag):
    f = RES / f"{tag}.json"
    return json.load(open(f)) if f.exists() else None


def compression_figure():
    survey = {r["state"]: r for r in json.load(open(RES / "survey.json"))}
    rows = []
    for st in ("t0", "t1", "t2", "t3", "t4", "t5", "t6"):
        full = rec(f"{st}_full_buffer_crop320_tau1.5") or rec(f"{st}_full_buffer_tau1.5")
        coarse = rec(f"{st}_coarse_porosity_buffer_tau1.5")
        rows.append(dict(state=st, strain=survey[st]["axial_strain_vs_t0"], porosity=survey[st]["porosity"],
                         k_full=full["k_sample_m2"] if full else None, k_coarse=coarse["k_sample_m2"] if coarse else None))
    fig, ax = plt.subplots(figsize=(6.2, 3.8)); ax2 = ax.twinx(); ax2.spines["top"].set_visible(False)
    load = [r for r in rows if r["state"] != "t6"]
    e = [r["strain"] for r in load]
    ax.plot(e, [r["k_full"] for r in load], "o-", color="#3b6ea5", label="k, 74.8 µm (scan resolution)")
    ax.plot(e, [r["k_coarse"] for r in load], "o--", color="#3b6ea5", mfc="white", label="k, 149.6 µm (2x coarsened)")
    t6 = rows[-1]
    if t6["k_full"]:
        ax.plot([t6["strain"]], [t6["k_full"]], "D", color="#2a8a4a", label="t6: unloaded after the test")
    ax2.plot(e, [r["porosity"] for r in load], "s:", color="#888888", label="porosity")
    ax.set_xlabel("axial strain (1 - L / L$_0$)"); ax.set_ylabel("axial permeability (m²)"); ax2.set_ylabel("porosity", color="#666666")
    ax.set_yscale("log"); h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, frameon=False, loc="lower left")
    ax.set_title("PU foam (10 ppi) under uniaxial compression: DaRUS-3010, states t0-t6", fontsize=9)
    fig.tight_layout(); OUT.mkdir(parents=True, exist_ok=True); fig.savefig(OUT / "foam_compression.png"); plt.close(fig)
    (RES / "compression_summary.json").write_text(json.dumps(rows, indent=1)); return rows


def sensitivity_figure():
    cases = [("coarse, quarter rule, periodic, tau 1.0", "t0_coarse_quarter_periodic_tau1"),
             ("coarse, quarter rule, periodic", "t0_coarse_quarter_periodic_tau1.5"),
             ("coarse, porosity rule, periodic", "t0_coarse_porosity_periodic_tau1.5"),
             ("coarse, porosity rule, mirror", "t0_coarse_porosity_mirror_tau1.5"),
             ("coarse, porosity rule, buffer", "t0_coarse_porosity_buffer_tau1.5"),
             ("coarse, porosity rule, buffer, crop 160", "t0_coarse_porosity_buffer_crop160_tau1.5"),
             ("scan resolution, buffer, crop 320", "t0_full_buffer_crop320_tau1.5"),
             ("scan resolution, periodic, crop 320", "t0_full_periodic_crop320_tau1.5")]
    rows = [(lab, rec(t)) for lab, t in cases]; rows = [(l, r) for l, r in rows if r]
    fig, ax = plt.subplots(figsize=(6.6, 0.42 * len(rows) + 1.0))
    y = np.arange(len(rows))[::-1]
    ax.barh(y, [r["k_sample_m2"] * 1e7 for _, r in rows], color=["#3b6ea5" if "scan" in l else "#9db4cf" for l, _ in rows])
    for yi, (l, r) in zip(y, rows):
        ax.text(r["k_sample_m2"] * 1e7, yi, f" {r['k_sample_m2']*1e7:.3f}", va="center", fontsize=7)
    ax.set_yticks(y); ax.set_yticklabels([l for l, _ in rows], fontsize=7); ax.set_xlabel("axial permeability (10$^{-7}$ m²), per sample area")
    ax.set_title("Undeformed foam t0: sensitivity to resolution, coarsening, end treatment, tau", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "foam_sensitivity.png"); plt.close(fig)
    return [(l, r["k_sample_m2"]) for l, r in rows]


if __name__ == "__main__":
    print(json.dumps(compression_figure(), indent=1)); print(json.dumps(sensitivity_figure(), indent=1))

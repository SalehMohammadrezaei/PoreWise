"""PU foam: axial permeability of the scanned core, per compression state.

    python examples/materials/pu_foam/run.py --state t0 --resolution full|coarse [--rule porosity|quarter]
           [--ends buffer|periodic|mirror] [--crop N] [--buffer-mm 1.5] [--tau 1.5] [--save-fields]

Flow along the core axis (z) driven by a body force, TRT collision, sealed sleeve. Reported:
    k_box_m2        PoreWise's k (mean velocity over the whole square grid)
    k_sample_m2     per sample cross-section (k_box / area fraction), buffer layers removed as a
                    series resistance; this is the permeability of the foam core
plus convergence and mass diagnostics, flux conservation through every z plane, flow regime,
and the share of the flow in the outer ring one nominal pore wide (wall effect).
"""
import argparse, json, math, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parents[2]))
from prepare import (STATES, VOXEL_M, DEFAULT_DATA, load, cylinder, disc_mask, coarsen, build_domain)  # noqa: E402
from common import sample_permeability, series_corrected, plane_fluxes, flow_regime, run_summary  # noqa: E402
from porewise import lbm_stokes_3d  # noqa: E402

FIELDS = HERE.parents[2] / "results" / "materials" / "foam"
PORE_M = 25.4e-3 / 10                     # nominal cell size of a 10 ppi foam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True, choices=list(STATES)); ap.add_argument("--resolution", default="full", choices=["full", "coarse"])
    ap.add_argument("--rule", default="porosity", choices=["porosity", "quarter"]); ap.add_argument("--ends", default="buffer", choices=["buffer", "periodic", "mirror"])
    ap.add_argument("--crop", type=int, default=0, help="central slices kept, at the run's resolution (0 = all)")
    ap.add_argument("--buffer-mm", type=float, default=1.5); ap.add_argument("--tau", type=float, default=1.5)
    ap.add_argument("--conv-tol", type=float, default=1e-6); ap.add_argument("--save-fields", action="store_true")
    ap.add_argument("--data", default=str(DEFAULT_DATA)); ap.add_argument("--out", default=str(HERE / "results"))
    a = ap.parse_args()
    t0 = time.time()
    skel = load(a.state, a.data); centre, radius = cylinder(skel)
    inside_full = disc_mask(skel.shape[1:], centre, radius)
    rec = dict(state=a.state, displacement_mm=STATES[a.state], resolution=a.resolution, ends=a.ends, tau=a.tau,
               porosity_full_resolution=float(1 - skel[:, inside_full].mean()), dataset="doi:10.18419/DARUS-3010 (CC BY 4.0)")
    if a.resolution == "coarse":
        skel, thr = coarsen(skel, 2, a.rule, inside_full); dx = 2 * VOXEL_M
        centre = ((centre[0] - 0.5) / 2, (centre[1] - 0.5) / 2); radius = radius / 2
        rec.update(coarsening_rule=a.rule, coarsening_threshold=thr)
    else:
        dx = VOXEL_M
    buf = int(round(a.buffer_mm * 1e-3 / dx)) if a.ends == "buffer" else 0
    solid, info = build_domain(skel, centre, radius, ends=a.ends, buffer_voxels=buf, crop=a.crop or None)
    del skel
    rec.update(geometry=info, dx_m=dx, fluid_nodes=int((~solid).sum()), grid_zyx=list(solid.shape))
    print(f"{a.state} {a.resolution} {a.ends}: grid {solid.shape}, porosity in cylinder {info['porosity_in_cylinder']:.4f}, fluid nodes {rec['fluid_nodes']/1e6:.1f} M", flush=True)
    nu = (a.tau - 0.5) / 3
    run = lbm_stokes_3d(solid, F_x=0.0, F_z=1e-6, backend="cuda-sparse", collision="trt", tau=a.tau, precision="float32",
                        n_steps_max=300000, conv_tol=a.conv_tol, return_fields=True, verbose=False, wall_timeout_s=8 * 3600)
    k_box_lu = run["u_z_mean_total"] * nu / 1e-6
    k_app_lu = sample_permeability(k_box_lu, info["sample_area_fraction"])
    if a.ends == "buffer":
        r = info["radius_voxels"]
        k_lu = series_corrected(k_app_lu, info["domain_slices"], info["sample_slices"], k_buffer=r * r / 8)   # open pipe in the sleeve
    else:
        k_lu = k_app_lu
    uz = run["uz"]
    zs = slice(buf, buf + info["sample_slices"]) if a.ends == "buffer" else slice(None)
    ins = disc_mask(solid.shape[1:], centre, info["radius_voxels"])
    rr = np.hypot(*np.meshgrid(np.arange(solid.shape[1]) - centre[0], np.arange(solid.shape[2]) - centre[1], indexing="ij"))
    ring = ins & (rr > info["radius_voxels"] - PORE_M / dx)
    q = np.asarray(uz[zs], dtype=np.float64).mean(axis=0)
    rec.update(valid=run.get("valid_for_permeability"), termination=run.get("termination_reason"), iterations=run.get("iterations"),
               k_box_m2=k_box_lu * dx * dx, k_sample_apparent_m2=k_app_lu * dx * dx, k_sample_m2=k_lu * dx * dx,
               diagnostics=run_summary(run).get("diagnostics"), plane_flux=plane_fluxes(uz, 0),
               regime=flow_regime(run, solid, nu, pore_length_lu=PORE_M / dx),
               wall_ring=dict(area_share=float(ring.sum() / ins.sum()), flow_share=float(q[ring].sum() / q[ins].sum())),
               wall_s=time.time() - t0)
    print(f"  k_sample {rec['k_sample_m2']:.4e} m2 (box {rec['k_box_m2']:.3e}), {rec['iterations']} it, {rec['termination']}, "
          f"flux CV {rec['plane_flux']['cv']:.1e}, {rec['wall_s']:.0f} s", flush=True)
    tag = f"{a.state}_{a.resolution}" + (f"_{a.rule}" if a.resolution == "coarse" else "") + f"_{a.ends}" + (f"_crop{a.crop}" if a.crop else "") + f"_tau{a.tau:g}"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); (out / f"{tag}.json").write_text(json.dumps(rec, indent=1, default=float))
    if a.save_fields:
        FIELDS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(FIELDS / f"{tag}.npz", solid=solid, ux=run["ux"].astype(np.float32), uy=run["uy"].astype(np.float32),
                            uz=run["uz"].astype(np.float32), centre=np.array(centre), radius=info["radius_voxels"], buffer=buf, dx=dx)


if __name__ == "__main__":
    main()

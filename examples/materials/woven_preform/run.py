"""Woven preform: in-plane permeability, two compaction states, resolution study.

    python examples/materials/woven_preform/run.py [--states relaxed compacted] [--dx 0.02 0.01 0.005]
           [--tau 1.0] [--save-fields 0.005] [--out examples/materials/woven_preform/results]

For each state and voxel size: build the geometry (geometry.py), drive the flow along x (warp)
and along y (weft) with a body force, TRT collision, and record per load: k (m2) and the
in-plane tensor, convergence and mass diagnostics, flux conservation through every plane, the
flow regime, and where the flow goes (next to the mould, in straight line-of-sight channels).
Fields of the runs at --save-fields voxel size go to results/materials/woven (not in git).
"""
import argparse, json, math, os, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parents[2]))
from geometry import STATES, stack_layout, voxelise, voxelise_periodic, periodic_stack_layout, plain_weave_alternates  # noqa: E402
from common import plane_fluxes, flow_regime, run_summary  # noqa: E402
from porewise import lbm_stokes_3d, k_lu_to_m2  # noqa: E402

FIELDS = HERE.parents[2] / "results" / "materials" / "woven"


def where_the_flow_goes(u, solid, axis):
    """Flow shares: first/last cavity layer (against the mould) and straight open channels."""
    q = u.sum()
    cav = slice(1, solid.shape[0] - 1)
    layer = np.asarray(u, dtype=np.float64).sum(axis=(1, 2))
    open_line = ~solid.any(axis=axis)                       # line of sight along the flow axis
    flow_in_open = np.asarray(u, dtype=np.float64).sum(axis=axis)[open_line].sum()
    return dict(mould_layer_bottom_share=float(layer[1] / q), mould_layer_top_share=float(layer[-2] / q),
                mould_layers_volume_share=float(2 / (solid.shape[0] - 2)),
                line_of_sight_area_share=float(open_line[cav].mean()), line_of_sight_flow_share=float(flow_in_open / q))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--states", nargs="+", default=list(STATES)); ap.add_argument("--dx", nargs="+", type=float, default=[0.02, 0.01, 0.005])
    ap.add_argument("--tau", type=float, default=1.5); ap.add_argument("--save-fields", type=float, default=0.005)
    ap.add_argument("--conv-tol", type=float, default=1e-6); ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--directions", nargs="+", default=["x", "y"]); ap.add_argument("--fields-tau", type=float, default=1.5)
    ap.add_argument("--stack", default="mould", choices=["mould", "periodic"], help="two plies in a closed mould, or an infinite laminate")
    a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True); FIELDS.mkdir(parents=True, exist_ok=True)
    for name in a.states:
        wv = STATES[name]; lay = stack_layout(wv) if a.stack == "mould" else periodic_stack_layout(wv)
        for dx in a.dx:
            solid, labels, info = (voxelise if a.stack == "mould" else voxelise_periodic)(wv, dx, lay)
            rec = dict(state=name, stack=a.stack, weave=wv.describe(), geometry=info, plain_weave_alternates=plain_weave_alternates(wv),
                       tau=a.tau, collision="trt", precision="float64", loads={})
            nu = (a.tau - 0.5) / 3
            u_mean = {}
            for d, force in (("x", dict(F_x=1e-6)), ("y", dict(F_x=0.0, F_y=1e-6))):
                if d not in a.directions:
                    continue
                t0 = time.time()
                run = lbm_stokes_3d(solid, backend="cuda-sparse", collision="trt", tau=a.tau, precision="float64",
                                    n_steps_max=400000, conv_tol=a.conv_tol, return_fields=True, verbose=False,
                                    wall_timeout_s=4 * 3600, **force)
                ax = {"x": 2, "y": 1}[d]; u = run["u" + d]
                # PoreWise averages over the whole grid; the preform's Darcy velocity is over the cavity only,
                # so the two mould-plate layers are taken out (no plates in the periodic stack)
                cav_factor = solid.shape[0] / (solid.shape[0] - 2) if a.stack == "mould" else 1.0
                kd = {c: k_lu_to_m2(run[f"u_{c}_mean_total"] * nu / 1e-6, dx * 1e-3) * cav_factor for c in "xy"}   # column of the tensor
                u_mean[d] = kd
                rec["loads"][d] = dict(
                    valid=run.get("valid_for_permeability"), termination=run.get("termination_reason"), iterations=run.get("iterations"),
                    k_m2=kd[d], k_cross_m2=kd["y" if d == "x" else "x"], normalisation="cavity", cavity_factor=cav_factor, wall_s=time.time() - t0,
                    diagnostics=run_summary(run).get("diagnostics"), plane_flux=plane_fluxes(u, ax),
                    regime=flow_regime(run, solid, nu, pore_length_lu=(wv.spacing_mm - wv.width_mm) / dx),
                    flow_paths=where_the_flow_goes(u, solid, ax) if a.stack == "mould" else None)
                print(f"{name} dx {dx*1e3:.1f} um load {d}: k {kd[d]:.4e} m2 (cross {rec['loads'][d]['k_cross_m2']:.2e}), "
                      f"{run.get('iterations')} it, {rec['loads'][d]['termination']}, flux CV {rec['loads'][d]['plane_flux']['cv']:.1e}, "
                      f"{time.time()-t0:.0f} s", flush=True)
                if abs(dx - a.save_fields) < 1e-12 and a.tau == a.fields_tau:
                    np.savez_compressed(FIELDS / f"{name}_{a.stack}_dx{dx*1e3:g}um_{d}.npz", solid=solid, labels=labels,
                                        ux=run["ux"].astype(np.float32), uy=run["uy"].astype(np.float32), uz=run["uz"].astype(np.float32))
            if "x" in u_mean and "y" in u_mean:
                K = np.array([[u_mean["x"]["x"], u_mean["y"]["x"]], [u_mean["x"]["y"], u_mean["y"]["y"]]])
                rec["tensor_inplane_m2"] = K.tolist()
                rec["tensor_asymmetry"] = float(abs(K[0, 1] - K[1, 0]) / np.abs(np.diag(K)).mean())
                ev = np.linalg.eigvalsh((K + K.T) / 2); rec["principal_m2"] = ev.tolist(); rec["anisotropy_ratio"] = float(ev.max() / ev.min())
            fn = out / (f"{name}_{a.stack}_dx{dx*1e3:g}um_tau{a.tau:g}" + ("" if set(a.directions) == {"x", "y"} else "_" + "".join(a.directions)) + ".json")
            fn.write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    main()

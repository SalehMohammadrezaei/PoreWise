"""Collect the materials-example results into one table of reported numbers.

    python examples/materials/summarise.py            # prints and writes examples/materials/summary.json

Every number quoted in docs/examples/*.md and in the README section comes from this file, which
is built only from the run records in examples/materials/*/results/*.json.
"""
import glob, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _j(p):
    return json.load(open(p)) if Path(p).exists() else None


def woven():
    out = {}
    for f in sorted(glob.glob(str(HERE / "woven_preform" / "results" / "*_tau*.json"))):
        d = json.load(open(f)); g = d["geometry"]; key = Path(f).stem
        stack = d.get("stack", "mould")
        # records written before run.py normalised by the cavity: the grid includes the two mould-plate layers
        fac = {k: (1.0 if v.get("normalisation") == "cavity" or stack != "mould" else g["grid_zyx"][0] / (g["grid_zyx"][0] - 2)) for k, v in d["loads"].items()}
        tens = d.get("tensor_inplane_m2")
        if tens and stack == "mould" and not all(v.get("normalisation") == "cavity" for v in d["loads"].values()):
            ff = g["grid_zyx"][0] / (g["grid_zyx"][0] - 2); tens = [[x * ff for x in r] for r in tens]
        loads = {k: dict(k_m2=v["k_m2"] * fac[k], k_cross_m2=v["k_cross_m2"] * fac[k], cavity_factor_applied=fac[k], iterations=v["iterations"], termination=v["termination"],
                         flux_cv=v["plane_flux"]["cv"], mach=v["regime"]["mach_max"], re_pore=v["regime"]["reynolds_pore"],
                         mass_drift=(v.get("diagnostics") or {}).get("mass_drift"), flow_paths=v.get("flow_paths"))
                 for k, v in d["loads"].items()}
        out[key] = dict(state=d["state"], stack=stack, tau=d["tau"], dx_um=g["dx_mm"] * 1e3, grid=g["grid_zyx"],
                        tow_volume_fraction=g["tow_volume_fraction"], fibre_volume_fraction=g["fibre_volume_fraction"],
                        intra_tow_vf=g["intra_tow_vf"], areal_weight_g_m2=g["areal_weight_g_m2"], crimp=g["crimp_factor"],
                        cavity_or_period_um=(g.get("cavity_mm") or g.get("period_mm")) * 1e3,
                        overlap_before_correction=g.get("overlap_before_correction_fraction"), tow_free_layers=g.get("tow_free_layers"),
                        weave=d["weave"], loads=loads, tensor=tens, asym=d.get("tensor_asymmetry"),
                        anisotropy=d.get("anisotropy_ratio"))
    return out


def foam():
    out = {"survey": _j(HERE / "pu_foam" / "results" / "survey.json")}
    for f in sorted(glob.glob(str(HERE / "pu_foam" / "results" / "t*_*.json"))):
        d = json.load(open(f)); key = Path(f).stem
        out[key] = dict(state=d["state"], resolution=d["resolution"], ends=d["ends"], tau=d["tau"], rule=d.get("coarsening_rule"),
                        crop=d["geometry"].get("crop_slices"), sample_slices=d["geometry"]["sample_slices"],
                        porosity_cylinder=d["geometry"]["porosity_in_cylinder"], porosity_full=d["porosity_full_resolution"],
                        k_sample_m2=d["k_sample_m2"], k_box_m2=d["k_box_m2"], iterations=d["iterations"], termination=d["termination"],
                        flux_cv=d["plane_flux"]["cv"], mach=d["regime"]["mach_max"], re_pore=d["regime"]["reynolds_pore"],
                        mass_drift=(d.get("diagnostics") or {}).get("mass_drift"), wall_ring=d["wall_ring"], fluid_nodes=d["fluid_nodes"])
    return out


if __name__ == "__main__":
    s = dict(woven=woven(), foam=foam())
    (HERE / "summary.json").write_text(json.dumps(s, indent=1, default=float))
    for k, v in s["woven"].items():
        print(f"woven {k:40s} Vf {v['fibre_volume_fraction']:.3f} " + " ".join(f"k_{d}={l['k_m2']:.4e}" for d, l in v["loads"].items()))
    for k, v in s["foam"].items():
        if k != "survey":
            print(f"foam  {k:45s} phi {v['porosity_cylinder']:.4f} k {v['k_sample_m2']:.4e} ({v['termination']}, {v['iterations']} it)")

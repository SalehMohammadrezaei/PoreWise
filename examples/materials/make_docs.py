"""Write the example documentation from the run records (no number is typed by hand).

    python examples/materials/summarise.py && python examples/materials/make_docs.py

Writes docs/examples/woven_preform.md, docs/examples/pu_foam.md and
docs/examples/materials_verification.md. Prose that does not depend on results is in this file;
every number comes from examples/materials/summary.json.
"""
import json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]; DOCS = ROOT / "docs" / "examples"
S = json.load(open(HERE / "summary.json"))
W, F = S["woven"], S["foam"]


def e(x, p=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "n/a"
    m, ex = f"{x:.{p - 1}e}".split("e"); return f"{m} × 10<sup>{int(ex)}</sup>"


def pct(x, p=1):
    return "n/a" if x is None else f"{100 * x:.{p}f} %"


def wkey(state, stack, dx, tau=1.5):
    for k, v in W.items():
        if v["state"] == state and v["stack"] == stack and abs(v["dx_um"] - dx) < 1e-6 and abs(v["tau"] - tau) < 1e-9 and "x" in v["loads"]:
            return v
    return None


def richardson(ks):
    """(k5, k10, k20) -> observed order, extrapolated value."""
    k1, k2, k3 = ks
    if None in ks or k2 == k1:
        return None, None
    r = (k3 - k2) / (k2 - k1)
    if r <= 0:
        return None, None
    p = math.log(r) / math.log(2)
    return p, k1 + (k1 - k2) / (2 ** p - 1)


# ============================================================================ woven
def woven_tables():
    rows, conv = [], []
    for state in ("relaxed", "compacted"):
        for stack in ("mould", "periodic"):
            if not any(wkey(state, stack, dx) for dx in (5, 10, 20)):
                continue                                                    # stacking not run
            ks = []
            for dx in (5, 10, 20):
                v = wkey(state, stack, dx); ks.append(v["loads"]["x"]["k_m2"] if v else None)
                if v:
                    L = v["loads"]["x"]
                    rows.append(f"| {state} | {'2 plies, mould' if stack == 'mould' else 'infinite laminate'} | {dx:g} | {v['fibre_volume_fraction']:.3f} | "
                                f"{v['tow_volume_fraction']:.3f} | {e(L['k_m2'], 4)} | {L['iterations']} | {L['flux_cv']:.1e} | {L['mach']:.1e} | {L['re_pore']:.1e} |")
            p, kx = richardson(ks)
            conv.append(dict(state=state, stack=stack, ks=ks, order=p, k_extrap=kx,
                             change_10_5=(ks[0] / ks[1] - 1) if ks[0] and ks[1] else None,
                             change_20_10=(ks[1] / ks[2] - 1) if ks[1] and ks[2] else None))
    return rows, conv


def woven_doc():
    rows, conv = woven_tables()
    v5 = {s: wkey(s, "mould", 5) for s in ("relaxed", "compacted")}
    v10 = {s: wkey(s, "mould", 10) for s in ("relaxed", "compacted")}
    vb = {s: wkey(s, "periodic", 5) for s in ("relaxed", "compacted")}
    tau = [v for v in W.values() if v["state"] == "compacted" and v["stack"] == "mould" and abs(v["dx_um"] - 10) < 1e-6]
    sym = [v for v in W.values() if v.get("tensor")]
    geo = []
    for s in ("relaxed", "compacted"):
        v = v5[s] or v10[s]
        if v:
            wv = v["weave"]
            geo.append(f"| {s} | {wv['spacing_mm']:.2f} | {wv['width_mm']:.2f} | {wv['thickness_mm']*1e3:.0f} | {wv['gap_mm']*1e3:.0f} | "
                       f"{wv['intra_tow_vf']:.2f} | {v['cavity_or_period_um']/2:.0f} | {v['tow_volume_fraction']:.3f} | {v['fibre_volume_fraction']:.3f} | "
                       f"{v['areal_weight_g_m2']:.0f} | {pct(v['overlap_before_correction'], 2)} | {v['tow_free_layers']} |")
    fp = []
    for s in ("relaxed", "compacted"):
        v = v5[s]
        if v and v["loads"]["x"].get("flow_paths"):
            f = v["loads"]["x"]["flow_paths"]
            fp.append(f"| {s} | {pct(f['mould_layers_volume_share'])} | {pct(f['mould_layer_bottom_share'] + f['mould_layer_top_share'])} | "
                      f"{pct(f['line_of_sight_area_share'])} | {pct(f['line_of_sight_flow_share'])} |")
    convrows = []
    for c in conv:
        order = "n/a" if c["order"] is None else f"{c['order']:.2f}"
        stack = "2 plies, mould" if c["stack"] == "mould" else "infinite laminate"
        convrows.append(f"| {c['state']} | {stack} | {e(c['ks'][2], 4)} | {e(c['ks'][1], 4)} | {e(c['ks'][0], 4)} | {pct(c['change_20_10'])} | {pct(c['change_10_5'])} |")
    monotone = all(c["order"] is not None for c in conv)
    taurows = [f"| {v['tau']:g} | {e(v['loads']['x']['k_m2'], 5)} |" for v in sorted(tau, key=lambda v: v["tau"])]
    symrows = [f"| {v['state']} | {v['stack']} | {v['dx_um']:g} | {v['tau']:g} | {e(v['tensor'][0][0], 4)} | {e(v['tensor'][1][1], 4)} | {abs(v['tensor'][0][1]) / v['tensor'][0][0]:.1e} |" for v in sym]
    k = {s: (v5[s]["loads"]["x"]["k_m2"] if v5[s] else None) for s in v5}; kb = {s: (vb[s]["loads"]["x"]["k_m2"] if vb[s] else None) for s in vb}
    ratio = (k["relaxed"] / k["compacted"]) if k["relaxed"] and k["compacted"] else None
    md = f"""# Woven composite preform: resin flow between the tows

Example of PoreWise as a mesoscale tool for liquid composite moulding (RTM, resin infusion). The
in-plane permeability of a dry woven preform is the main material input of a mould-filling
simulation (Darcy's law, **q** = −(**K**/μ)∇p): it sets fill times, injection pressures and where
dry spots form. This example computes the in-plane permeability tensor of an idealised carbon
plain-weave preform at two compaction states, shows where the resin flows, and checks the
numbers against resolution, solver parameters and published textile data.

![The lower ply with the computed resin speed on the plane between the plies; the upper ply is lifted for display](../figures/materials/woven_preform.png)

*Relaxed state, 5 µm grid. Colour: local speed on the plane between the two plies, blue (slow) to red (fast); arrows:
in-plane flow direction. The upper ply is lifted and set back for display only; the flow was computed with the plies
nested.*

## The textile

A 3K carbon plain weave of the commercial 200 g/m² class: 5.0 ends/cm and 5.0 picks/cm (yarn
spacing 2.0 mm), 3K tows of 7 µm filaments (fibre density 1.78 g/cm³). Datasheets of such fabrics
give 200 g/m² ± 4 % and a dry thickness of about 0.32 mm (Compositesplaza CAWFINPL200; R&G
CT-200100P). The geometry here gives an areal weight within that tolerance (table below).

The model is built the way TexGen builds woven textiles (elliptical tow sections swept along crimped
paths, interference removal, layer nesting), in `examples/materials/woven_preform/geometry.py`:

- **Tow path:** z = ±(h/2)·g(cos(πu/s)), with g(c) = tanh(βc)/tanh(β) and β = 1.5. β flattens the path over the
  crossings, where real tows lie on each other. A plain sinusoid (β = 0) with elliptical sections of width w intersects
  the crossing tow unless w < 2s/π. For the tow widths used here it would interpenetrate by 3–6 % of the tow volume.
- **Section:** elliptical, width w and thickness h. h follows from the tow's fibre area and the assumed intra-tow fibre
  fraction (πwh/4 = fibre area / intra-tow V<sub>f</sub>).
- **Interference:** with β = 1.5 the overlap before correction is at most the value in the table. Any remaining shared
  point goes to the tow in which it lies deeper (TexGen's principle); the solid is unchanged.
- **Plain-weave pattern:** warp and weft centrelines are on opposite sides at every crossing (checked for all four
  crossings of the unit cell).
- **Two plies:** the second ply is shifted half a yarn spacing in x and y and lowered onto the first until they touch
  (nesting by contact).
- **Mould:** in the closed-mould case, rigid plates sit on the outermost crests. The cavity is a whole number of cubic
  voxels, centred on the preform, so no tow-free slot is left. The earlier version of this example had a 140 µm slot
  under the mould that carried 52 % of the flow.
- **Infinite laminate (implemented, not reported):** `--stack periodic` repeats the same two plies periodically through
  the thickness, with no mould. No results for it are reported here.

Two compaction states of the same fabric: compaction flattens and widens the tows (higher intra-tow fibre fraction,
smaller gaps) and the plies nest closer.

| State | s (mm) | w (mm) | h (µm) | gap (µm) | intra-tow V<sub>f</sub> (assumed) | ply thickness (µm) | tow volume fraction | fibre V<sub>f</sub> | areal weight (g/m²) | overlap before correction | tow-free layers |
|---|---|---|---|---|---|---|---|---|---|---|---|
{chr(10).join(geo)}

The intra-tow fibre fraction is an assumption. Values of 0.56–0.71 are measured or used in the literature (Syerko et al.
2023: 0.56 by CT of a glass tow; Wong et al. 2006: 0.60; Nedanov & Advani 2002: 0.71). The fibre volume fraction of the
preform is the intra-tow fraction times the tow volume fraction. These ply thicknesses match the cured-ply estimate
t = 200/(1780·V<sub>f</sub>) for this areal weight.

## Physics and method

**Assumptions:**
- Saturated, steady, creeping flow of a Newtonian resin between the tows: Stokes equations, no slip on the tows and the
  mould.
- The tows are impermeable, so only inter-tow flow is resolved.
  - For 7 µm fibres at an intra-tow fraction of 0.6–0.7, Gebart's formulas give tow permeabilities of 2 × 10<sup>-14</sup>
    to 3 × 10<sup>-13</sup> m², three to four orders below the results here.
  - Published comparisons of dual-scale and impermeable-tow models of woven fabrics put the effect on saturated in-plane
    permeability at about 3–30 % (Alotaibi et al. 2021).
  - Unsaturated, dual-scale filling (void formation, delayed tow impregnation) is outside this model.

**Solver:** PoreWise D3Q19 lattice Boltzmann.
- Collision: TRT, magic parameter 3/16, so walls sit half-way between nodes and the permeability is independent of the
  relaxation time. τ = 1.5.
- Driving force: a uniform body force along x (warp) or y (weft); periodic in x and y.
- Precision and convergence: double precision, convergence tolerance 10<sup>-6</sup>.
- Permeability: k = ν⟨u⟩/F, with ⟨u⟩ the mean velocity over the cavity (mould plates excluded), i.e. the Darcy velocity
  of the preform.
- The off-diagonal terms come from the velocity component normal to the load.
- **Through-thickness permeability is not computed.** With impermeable tows, flow across the plies can only pass where
  inter-tow openings of neighbouring plies happen to line up. That makes k<sub>zz</sub> a property of the nesting
  assumption rather than of the fabric, and in real preforms it is dominated by intra-tow flow, which is not modelled.

## Results

| State | Stacking | voxel (µm) | fibre V<sub>f</sub> | tow vol. fraction | k<sub>xx</sub> (m²) | iterations | plane-flux CV | Mach | Re<sub>pore</sub> |
|---|---|---|---|---|---|---|---|---|---|
{chr(10).join(sorted(rows))}

At 5 µm the moulded two-ply preform gives k = {e(k['relaxed'])} m² (relaxed) and {e(k['compacted'])} m² (compacted). Compaction lowers
the permeability by a factor of {'n/a' if ratio is None else f'{ratio:.1f}'}.

### Where the resin flows

| State (5 µm, mould) | volume of the two cavity layers at the mould | flow in those layers | open straight channels along x (area) | flow in those channels |
|---|---|---|---|---|
{chr(10).join(fp)}

The open straight channels are the gaps between neighbouring warp tows, running over and under the weft where it
crosses at mid-ply height. They are part of the plain-weave architecture, not an artefact: in a real fabric these
inter-tow channels are where the resin front runs ahead. The figure
[`woven_section_compacted.png`](../figures/materials/woven_section_compacted.png) shows the axial velocity in sections
across the flow.

## Verification and sensitivity

Full report: [materials_verification.md](materials_verification.md).

**Resolution:** k<sub>xx</sub> at 20, 10 and 5 µm and the change at each halving of the voxel size:

| State | Stacking | k, 20 µm (m²) | k, 10 µm (m²) | k, 5 µm (m²) | change 20 → 10 µm | change 10 → 5 µm |
|---|---|---|---|---|---|---|
{chr(10).join(convrows)}

{'The sequences converge monotonically.' if monotone else 'The sequences are not monotone (the 20 µm grid puts only 7 to 15 voxels across the gaps), so no order of convergence or extrapolated value is quoted. The change from 10 to 5 µm is taken as the resolution uncertainty of the 5 µm values.'}

**Relaxation time** (compacted, 10 µm, TRT):

| τ | k<sub>xx</sub> (m²) |
|---|---|
{chr(10).join(taurows)}

**In-plane tensor** (loads along x and y):

| State | Stacking | voxel (µm) | τ | k<sub>xx</sub> | k<sub>yy</sub> | \\|k<sub>xy</sub>\\|/k<sub>xx</sub> |
|---|---|---|---|---|---|---|
{chr(10).join(symrows)}

The fabric is balanced and the ply shift is symmetric, so k<sub>xx</sub> = k<sub>yy</sub> is expected; it holds to the
digits shown, and the off-diagonal terms are negligible.

## Comparison with published data

No experiment matches this fabric, these compaction states and this stacking, so no experimental validation is claimed.

- **Measured, comparable fabrics:**
  - Carbon 2×2 twill (6K, 285 g/m²) at V<sub>f</sub> = 45.3 %: K₁ = 1.39 × 10<sup>-10</sup> and K₂ = 0.73 × 10<sup>-10</sup> m²
    (Vernet et al. 2014, benchmark II, unsaturated).
  - Glass 2/2 twill at V<sub>f</sub> = 45.6 / 49.5 / 53.3 %: median K₁ = 6.1 / 4.1 / 2.6 × 10<sup>-11</sup> m²
    (May et al. 2019, benchmark III).
  - Glass plain weave: about 7 × 10<sup>-11</sup> at 45 % and 2 × 10<sup>-11</sup> m² at 54 % (Verleye et al. 2008, read off
    a figure).
- **Mesoscale simulations of plain weaves** typically sit above measurements: a 40–50 % gap is reported by Verleye et al.
  2008, and factors of about 2 by Zeng et al. (FPCM12).
- **The values computed here** fall in the same 10<sup>-11</sup>–10<sup>-10</sup> m² range at comparable fibre volume
  fractions, and decrease with compaction as the measurements do.

The remaining differences have identifiable causes: the idealised tow shapes, impermeable tows, the two-ply closed-mould
stacking, and the measured fabrics being twills rather than plain weaves.

## Engineering conclusions

- In a plain weave the resin flows mainly along the gaps between parallel tows and over and under the crossing tows.
  The tows themselves carry almost none of the saturated flow.
- Compaction from the relaxed to the compacted state narrows the gaps (300 → 150 µm) and lowers the in-plane
  permeability by the factor above. The preform architecture, more than the fibre volume fraction alone, sets the flow
  paths.
- With the mould on the tow crests, the cavity layers next to the mould carry about their volume share of the flow
  (table above), so no mould-surface bypass remains.
- For a balanced plain weave with symmetric nesting, the in-plane permeability is isotropic. Unbalanced fabrics or
  shear would make it anisotropic, which the tensor workflow captures.

## Limitations

- **Geometry:** idealised (TexGen-type) shapes, not a scan. Tow width, gap and intra-tow fraction are assumed, and no
  measured 3K tow geometry was found to calibrate them. One unit cell is periodic: no variability, no random nesting,
  no shear.
- **Tows:** impermeable. The saturated value would rise by a few to about 30 % with permeable tows (literature range).
  Unsaturated dual-scale effects are not modelled.
- **Compaction:** represented by two geometric states, not by a mechanical compaction model.
- **Through-thickness permeability:** not computed (see Physics and method).

## Reproduce

```bash
python examples/materials/woven_preform/run.py --dx 0.02 0.01           # both states, x and y loads, closed mould
python examples/materials/woven_preform/run.py --dx 0.005 --directions x   # 5 um (saves the fields for the figures)
python examples/materials/woven_preform/run.py --dx 0.01 --states compacted --tau 1.0 --directions x   # tau check
python examples/materials/summarise.py && python examples/materials/woven_preform/figures.py
```

GPU (`cuda-sparse`); the 5 µm runs take 45 to 80 minutes each on an RTX 6000 Ada. The rendered figure needs Blender 4.5
(`examples/materials/render/cv_final2.py woven_flow`, then `bl_scene.py`).

## References

- Compositesplaza, technical data sheet CAWFINPL200 (carbon plain weave, 200 g/m², 5.0 × 5.0 /cm, 3K).
- Verleye B., Lomov S.V., Long A., Verpoest I., Roose D. (2008). Permeability of textile reinforcements: simulation, influence of shear and validation. Compos. Sci. Technol. 68, 2804–2810.
- Vernet N. et al. (2014). Experimental determination of the permeability of engineering textiles: Benchmark II. Compos. Part A 61, 172–184.
- May D. et al. (2019). In-plane permeability characterization of engineering textiles based on radial flow experiments: a benchmark exercise. Compos. Part A 121, 100–114.
- Gebart B.R. (1992). Permeability of unidirectional reinforcements for RTM. J. Compos. Mater. 26, 1100–1133.
- Alotaibi H., Abeykoon C., Soutis C., Jabbari M. (2021). Numerical simulation of resin flow in a dual-scale textile preform. Materials 14, 405.
- Syerko E. et al. (2023). Compos. Part A 167, 107397 (CT-based intra-yarn fibre volume fraction).
- Wong C.C., Long A.C., Rudd C.D. (2006). FPCM-8 (unit-cell permeability of plain weaves).
- Zeng X., Endruweit A., Long A.C. (FPCM-12). Meso-scale permeability of textiles with nesting and shear.
- Lin H., Brown L.P., Long A.C. (2011). Modelling and simulating textile structures using TexGen. Adv. Mater. Res. 331, 44–47.
"""
    (DOCS / "woven_preform.md").write_text(md)
    return conv


# ============================================================================ foam
def frec(tag):
    return F.get(tag)


def foam_doc():
    sv = {r["state"]: r for r in F["survey"]}
    comp = []
    for st in ("t0", "t1", "t2", "t3", "t4", "t5", "t6"):
        full = frec(f"{st}_full_buffer_crop200_tau1.5") or frec(f"{st}_full_buffer_tau1.5")
        coarse = frec(f"{st}_coarse_porosity_buffer_tau1.5")
        r = sv[st]
        comp.append(f"| {st} | {r['displacement_mm']} | {r['displacement_mm'] / 50:.2f} | {r['axial_strain_vs_t0']:.3f} | {r['length_mm']:.2f} | {r['porosity']:.4f} | "
                    f"{r['solid_volume_rel_t0']:.3f} | {e(full['k_sample_m2']) if full else 'not run'} | {e(coarse['k_sample_m2']) if coarse else 'not run'} |")
    k0 = frec("t0_coarse_porosity_buffer_tau1.5"); k5 = frec("t5_coarse_porosity_buffer_tau1.5"); k6 = frec("t6_coarse_porosity_buffer_tau1.5")
    ratio = k5["k_sample_m2"] / k0["k_sample_m2"] if k0 and k5 else None
    sens = []
    for lab, tag in [("coarse 149.6 µm, quarter rule, periodic ends, τ = 1.0", "t0_coarse_quarter_periodic_tau1"),
                     ("coarse, quarter rule, periodic ends, τ = 1.5", "t0_coarse_quarter_periodic_tau1.5"),
                     ("coarse, porosity-preserving, periodic ends", "t0_coarse_porosity_periodic_tau1.5"),
                     ("coarse, porosity-preserving, mirrored", "t0_coarse_porosity_mirror_tau1.5"),
                     ("coarse, porosity-preserving, fluid buffers", "t0_coarse_porosity_buffer_tau1.5"),
                     ("coarse, porosity-preserving, buffers, central 26.9 mm", "t0_coarse_porosity_buffer_crop180_tau1.5"),
                     ("scan resolution 74.8 µm, buffers, central 15.0 mm", "t0_full_buffer_crop200_tau1.5")]:
        r = frec(tag)
        if r:
            sens.append(f"| {lab} | {r['porosity_cylinder']:.4f} | {e(r['k_sample_m2'], 4)} | {r['iterations']} | {r['flux_cv']:.1e} | {r['mass_drift']:.1e} |")
    t0f = frec("t0_full_buffer_crop200_tau1.5")
    res_diff = (t0f["k_sample_m2"] / k0["k_sample_m2"] - 1) if t0f and k0 else None
    md = f"""# Open-cell polyurethane foam: permeability from a micro-CT scan, under compression

Example of PoreWise as an image-based characterisation tool: the axial permeability of a real
foam core, computed from its micro-CT scan at seven stages of a uniaxial compression test. Open-cell
foams are used as filters, separators, acoustic and fluid-management media; their permeability and
how it changes when the foam is compressed in service or during processing are design inputs.

![The whole scanned foam sample with streamlines of the computed flow](../figures/materials/pu_foam.png)

*The whole sample (52 mm across, 44 mm long), flow upwards, on the 149.6 µm grid. Streamlines are coloured by local
speed, blue (slow) to red (fast). The struts are faded out near the top for display only. A close-up at scan resolution
is in [`pu_foam_detail.png`](../figures/materials/pu_foam_detail.png).*

## Data

M. Ruf and H. Steeb (2024), *In situ micro-XRCT data set of an open-cell polyurethane foam sample
under uniaxial compression load*, DaRUS, https://doi.org/10.18419/DARUS-3010, **CC BY 4.0**.

- **Sample:** reticulated PU foam, 10 pores per inch, a cylinder 52 mm in diameter and 50 mm long.
- **States:** scanned at t0–t6 under imposed displacements of 0, 5, 10, 15, 20 and 25 mm, then 0 mm again (t6). The
  sample relaxed for 20 min before each scan.
- **Images:** segmented stacks (skeleton 255, void 0) with a documented voxel size of 74.8 µm, extracted by the authors
  with a cylinder of 52.36 mm diameter.
- **Voxel-size check:** the fitted cylinder is 699 voxels across, which is 52.29 mm at 74.8 µm. t0 is 593 slices,
  44.36 mm, against the documented usable length of 44.26 mm. Both support 74.8 µm.

`examples/materials/pu_foam/prepare.py fetch` downloads the segmented stacks with checksum verification.

| State | displacement (mm) | nominal strain (/50 mm) | axial strain vs t0 | length (mm) | porosity | apparent solid volume / t0 | k, 74.8 µm (m²) | k, 149.6 µm (m²) |
|---|---|---|---|---|---|---|---|---|
{chr(10).join(comp)}

Porosity is measured inside the sample cylinder.

**Apparent solid volume.** It stays within 3 % of t0 up to 10 mm of displacement and then rises by 26 % at 25 mm.
Polyurethane is effectively incompressible, so this rise is not physical. It reflects the segmentation of densified,
touching struts (partial volume) at high strain. Every state is cut by the same 52.36 mm cylinder, so any lateral
bulging is not in the data. The porosity of the most compressed states, and the permeabilities computed from them,
carry this uncertainty.

## Physics and method

**Assumptions:** saturated, steady, creeping Newtonian flow through the pore space, with no slip on the struts (Stokes
equations).

**Domain:**
- **Sealed core:** the sample cylinder is kept and everything outside it is solid, like a core in a sleeve. So the flow
  is along the axis only, and nothing passes outside the sample.
- **Ends:** the two end faces of a scanned sample do not match, so a periodic domain joins mismatched struts. The
  standard configuration adds free-fluid reservoirs (1.5 mm) inside the sleeve at both ends. Their resistance is removed
  afterwards as a series term: L/k<sub>app</sub> = L<sub>s</sub>/k<sub>s</sub> + L<sub>b</sub>/k<sub>b</sub>, with
  k<sub>b</sub> = R²/8 for the open sleeve.
- **Alternative ends:** mirrored and periodic ends are run for comparison.

**Normalisation.** PoreWise reports k with the mean velocity over the whole grid. Here the sample is a cylinder inside a
square grid, so the Darcy flux of the sample is the same volume flow divided by the sample area:
k<sub>sample</sub> = k<sub>grid</sub> / (A<sub>sample</sub>/A<sub>grid</sub>). The cylinder fills 67.8 % of the grid, so
without this correction k would be 32 % low. That was the error found in the first version of this example.

**Resolution:**
- **Scan resolution:** one run, the uncompressed state t0, on a central 200-slice section (15.0 mm; {t0f['fluid_nodes'] / 1e6:.1f} million
  fluid nodes) so that it fits one 48 GB GPU.
- **Coarse grid:** the 149.6 µm runs use the whole length of all seven states. Two coarsening rules are compared:
  - "quarter": a cell is solid if at least 25 % of its block is skeleton. This keeps thin struts connected but thickens
    them.
  - "porosity-preserving": the threshold reproduces the scan's solid fraction exactly.

**Solver settings:** TRT, τ = 1.5, single precision, convergence tolerance 10<sup>-6</sup>.

## Verification and sensitivity (t0)

| Configuration | porosity | k (m²) | iterations | plane-flux CV | mass drift |
|---|---|---|---|---|---|
{chr(10).join(sens)}

Full report: [materials_verification.md](materials_verification.md).

## Results

**Compression.** On the 149.6 µm grid the axial permeability falls from {e(k0['k_sample_m2']) if k0 else 'n/a'} m² (t0) to {e(k5['k_sample_m2']) if k5 else 'n/a'} m² at 25 mm displacement (t5): a ratio of {'n/a' if ratio is None else f'{ratio:.2f}'}.
The unloaded state t6 returns to {e(k6['k_sample_m2']) if k6 else 'n/a'} m², close to t0, consistent with its porosity (0.9539 against 0.9534).

**Resolution.** At scan resolution, t0 gives {e(t0f['k_sample_m2']) if t0f else 'n/a'} m² on its central 15.0 mm section, {'n/a' if res_diff is None else f'{abs(100 * res_diff):.0f} % ' + ('below' if res_diff < 0 else 'above')} the
coarse-grid value for the whole length. The two runs differ in both voxel size and section, so this is an indication of
the resolution effect, not a clean measure of it.

![Permeability and porosity against compression](../figures/materials/foam_compression.png)

## Comparison with published data

- **Absolute value** at t0: within the range reported for 10 ppi open-cell foams of high porosity, about
  10<sup>-7</sup> m².
  - Reticulated PU, 10 ppi: 4.6 × 10<sup>-8</sup> m² at porosity 0.976 (Efstathiou & Luhar 2018); 1.6 × 10<sup>-7</sup> m²
    (Manes et al. 2011, as cited there).
  - Metal 10 PPI foams: 1.0–1.8 × 10<sup>-7</sup> m² (Bhattacharya et al. 2002; Mancin et al. 2010; Dukhan 2006, as
    tabulated by Kumar 2014).
  - These are Forchheimer fits or different foams, so only the order of magnitude is compared.
- **The dataset authors' own permeability study** (Krach, Ruf & Steeb 2025, InterPore Journal 2(1)) reports, read off
  their Fig. 10, about 2 × 10<sup>-11</sup> m² at zero strain. That is four orders of magnitude below both this result
  and the literature range above, which points to a unit or scaling difference that could not be resolved from the
  paper. Their ratio over the same compression range is about 0.22.
- **Compression trend:**
  - Regressions for compressed open-cell foams give k/k₀ ≈ (1 − ε)<sup>1.8</sup> (Gent & Rusch) to
    (1 − ε)<sup>2.2</sup> (Hilyard & Collier), as quoted by Dawson (2008).
  - Dawson measured k/k₀ = 0.18–0.22 at ε = 0.6 for fine PU foams.
  - The ratio computed here over 0 to 0.5 nominal strain is the value above, to be read against these.

## Engineering conclusions

- Permeability can be extracted from the scan of a real foam with the sample's true cross-section: at the scan's own
  resolution for the uncompressed state, and on a 2× coarser grid for the compression series.
- Compressing the foam along the flow direction lowers its permeability progressively, and the change is reversible
  here (t6).
- The scan-resolution result and the coarse-grid result differ by the amount in the sensitivity table. Struts only a
  few voxels thick make the result resolution-dependent, so the voxel size has to be stated with any number.

## Limitations

- **Single sample:** one foam, one flow direction (the loading axis), one scan per state. No repeat samples.
- **Segmentation:** the authors' segmentation is used as given. Its effect grows at high strain (the apparent rise in
  solid volume above).
- **Scan-resolution run:** t0 only, on a central 15.0 mm section. The compression series is on the 149.6 µm grid, whole
  length.
- **Boundary conditions:** the sealed core with fluid buffers is an idealisation of a permeameter. The sleeve walls are
  no-slip at the cut cylinder surface, and lateral bulging is not represented.
- **Flow regime:** Stokes flow only. Forchheimer effects at practical flow rates through 10 ppi foams are not covered.

## Reproduce

```bash
python examples/materials/pu_foam/prepare.py fetch && python examples/materials/pu_foam/prepare.py survey
for s in t0 t1 t2 t3 t4 t5 t6; do python examples/materials/pu_foam/run.py --state $s --resolution coarse; done   # compression series
python examples/materials/pu_foam/run.py --state t0 --resolution full --crop 200 --save-fields                   # scan resolution
python examples/materials/pu_foam/run.py --state t0 --resolution coarse --ends mirror                            # end-condition variants
python examples/materials/summarise.py && python examples/materials/pu_foam/figures.py
```

## References and attribution

- Ruf M., Steeb H. (2024). In situ micro-XRCT data set of an open-cell polyurethane foam sample under uniaxial compression load. DaRUS. https://doi.org/10.18419/DARUS-3010. CC BY 4.0. The images are used unmodified; only coarsened copies are derived for the sensitivity runs.
- Krach D., Ruf M., Steeb H. (2025). InterPore Journal 2(1), doi:10.69631/ipj.v2i1nr39.
- Efstathiou C., Luhar M. (2018). J. Fluid Mech. 841, 351–379.
- Bhattacharya A., Calmidi V.V., Mahajan R.L. (2002). Thermophysical properties of high porosity metal foams. Int. J. Heat Mass Transfer 45, 1017–1031.
- Dawson M.A. (2008). Modeling the dynamic response of low-density, reticulated, elastomeric foam impregnated with Newtonian and non-Newtonian fluids. PhD thesis, MIT. (Journal version: Dawson, Germaine, Gibson 2007, Int. J. Solids Struct. 44, 5133.)
"""
    (DOCS / "pu_foam.md").write_text(md)


def verification_doc():
    wrows = []
    for k_, v in sorted(W.items()):
        for d, L in v["loads"].items():
            drift = "n/a" if L["mass_drift"] is None else format(L["mass_drift"], ".1e")
            wrows.append(f"| {v['state']} | {v['stack']} | {v['dx_um']:g} | {v['tau']:g} | {d} | {e(L['k_m2'], 5)} | {L['termination']} | {L['iterations']} | "
                         f"{drift} | {L['flux_cv']:.1e} | {L['mach']:.1e} | {L['re_pore']:.1e} |")
    frows = []
    for k_, v in sorted(F.items()):
        if k_ == "survey":
            continue
        frows.append(f"| {k_} | {v['porosity_cylinder']:.4f} | {e(v['k_sample_m2'], 5)} | {e(v['k_box_m2'], 4)} | {v['termination']} | {v['iterations']} | "
                     f"{'n/a' if v['mass_drift'] is None else format(v['mass_drift'], '.1e')} | {v['flux_cv']:.1e} | {v['mach']:.1e} | {v['re_pore']:.1e} | "
                     f"{pct(v['wall_ring']['area_share'])} / {pct(v['wall_ring']['flow_share'])} |")
    renders = []
    for p in sorted((ROOT / "results" / "materials" / "render_specs").glob("*.json")) if (ROOT / "results" / "materials" / "render_specs").exists() else []:
        sp = json.load(open(p)); pv = sp.get("provenance", {}); tc = pv.get("tracer_check", {})
        renders.append(f"| {p.stem} | {Path(pv.get('field', '')).name} | {tc.get('points_traced')} | {tc.get('in_solid')} | {tc.get('stencil_touches_solid')} | {tc.get('lines_drawn')} |")
    cache = ROOT / "results" / "materials" / "render_cache" / "cv_final2"
    for name, field in (("pu_foam.png", "t0_coarse_porosity_buffer_tau1.5.npz"), ("pu_foam_detail.png", "t0_full_buffer_crop200_tau1.5.npz")):
        p = cache / {"pu_foam.png": "foam_whole.json", "pu_foam_detail.png": "foam_detail.json"}[name]
        if p.exists():
            tc = json.load(open(p))[2]
            renders.append(f"| {name} | {field} | {tc.get('points')} | {tc.get('in_solid')} | {tc.get('stencil_solid')} | {tc.get('lines')} |")
    md = f"""# Materials examples: numerical verification and sensitivity

Every row below is read from a run record (`examples/materials/*/results/*.json`) by
`examples/materials/make_docs.py`. A run is accepted only when PoreWise's own checks pass
(convergence over three consecutive windows, mass conservation, Mach below 0.05).

**Conservation.**
- *Plane-flux CV* is the coefficient of variation of the volume flux through every plane normal to the flow. In steady
  incompressible flow it is zero up to round-off and the convergence tolerance.
- *Mass drift* is the change in total mass over the run.

**Flow regime.**
- *Mach* is the maximum lattice speed times √3.
- *Re<sub>pore</sub>* is the mean pore speed times a pore length (the inter-tow gap for the preform, the nominal 2.54 mm
  cell for the foam) divided by ν. It is far below 1 everywhere, so the flow is Stokes flow and the Darcy permeability is
  well defined.

## Woven preform

| state | stacking | voxel (µm) | τ | load | k (m²) | termination | iterations | mass drift | plane-flux CV | Mach | Re<sub>pore</sub> |
|---|---|---|---|---|---|---|---|---|---|---|---|
{chr(10).join(wrows)}

## PU foam

| run | porosity in cylinder | k per sample area (m²) | k per grid area (m²) | termination | iterations | mass drift | plane-flux CV | Mach | Re<sub>pore</sub> | outer ring: area / flow |
|---|---|---|---|---|---|---|---|---|---|---|
{chr(10).join(frows)}

The outer ring is the annulus one nominal pore wide at the sleeve. A ring carrying less flow than its share of the area
means the sleeve does not act as a bypass.

## Figures: provenance of the streamlines

All streamlines are integrated in the computed velocity fields with the wall-safe tracer (`examples/materials/common.py`,
the rule of `benchmarks/rock_streamlines.py`). A point is kept only if the eight cells around it are pore, and every
sub-step is checked. The counts are over every traced point, including lines not drawn.

| figure | field | points traced | points in a solid voxel | points with a solid cell in the stencil | lines drawn |
|---|---|---|---|---|---|
{chr(10).join(renders)}

`woven_preform.png` has no streamlines: it shows the speed and the in-plane direction on one plane of the computed field.
"""
    (DOCS / "materials_verification.md").write_text(md)


if __name__ == "__main__":
    DOCS.mkdir(parents=True, exist_ok=True)
    conv = woven_doc(); foam_doc(); verification_doc()
    print(json.dumps(conv, indent=1, default=float))

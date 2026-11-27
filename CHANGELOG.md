# Changelog

## Unreleased

### Materials engineering examples (`examples/materials/`)
Two worked examples that apply PoreWise beyond rock images. The solver and the public API are
unchanged; everything below is example code, tests and documentation.

- **Woven composite preform** (`woven_preform/`): an idealised 3K carbon plain-weave
  generator.
  - Elliptical tows on crimped paths, interference removal, nesting by contact, mould on the
    tow crests (no tow-free slot), and an infinite-laminate variant.
  - In-plane permeability tensor at two compaction states for the closed-mould stacking, with a
    resolution study (20, 10, 5 um), a tau check and a flow-path analysis.
- **Open-cell PU foam** (`pu_foam/`): the micro-CT scans of Ruf & Steeb (DaRUS-3010,
  CC BY 4.0).
  - Download with checksums, cylinder fitting, a sealed-core domain, porosity-preserving
    coarsening, and fluid-buffer, mirrored or periodic ends.
  - Axial permeability per sample cross-section for seven compression states on a 2x coarsened
    grid, and for the uncompressed state at scan resolution.
- **Shared helpers** (`common.py`):
  - Permeability per sample cross-section (PoreWise's k is per grid cross-section) and a series
    correction for buffer layers.
  - Plane-flux conservation and flow-regime checks.
  - A vectorised form of the wall-safe streamline rule of `benchmarks/rock_streamlines.py`.
- **Figures and docs:** analysis figures, Blender scene scripts for the rendered figures
  (`render/`, optional), documentation generated from the run records (`make_docs.py`,
  `docs/examples/`) and a Materials examples section in the README.
- **Tests:** `tests/test_materials_examples.py`.
  - Weave topology, mould contact, periodicity, interference correction, fibre bookkeeping.
  - Normalisation (including a pipe-in-sleeve Poiseuille check), tracer equivalence with the
    reference rule, and no traced point next to a solid.
  - Foam domain construction and porosity-preserving coarsening.

## 0.1.0

First public release. The project was developed as LBM-Permeability and renamed PoreWise;
`import lbm_permeability` still works and warns.

### Capabilities
- Absolute permeability in one direction, or the full tensor with principal values and
  directions, from 2D (D2Q9) and 3D (D3Q19) binary images, periodic body-force-driven flow
- BGK and two-relaxation-time (magic parameter 3/16) collision
- Backends: `numpy`, `cupy-array`, `cuda` (dense), `cuda-sparse` and `numba-sparse`
  (pore voxels only, half-way bounce-back, in-place streaming, float32 or float64 storage)
- Acceptance checks on convergence, mass, density and Mach number; JSON output with
  provenance; VTI field export; porosity, connectivity and local-thickness tools
- Command line interface `python -m porewise`

### Verification
See `docs/verification.md`: plane channel, pipes, simple cubic sphere arrays (Zick and Homsy
1982), inclined-slit tensor test, the micromodel benchmark of Wagner et al. (2021) and the
1024-cube images of Saxena et al. (2017).

### Limitations
Periodic boundaries only; single-phase creeping flow; binary images with isotropic voxels;
one GPU (no MPI or multi-GPU); pore-only backends are 3D only; the 2D pressure-driven solver
is experimental.

### Results in the accompanying paper
The verification ladder and the rock benchmarks were run between commits `b89d47f` and
`a2efed9` (before and after the rename); the boundary and refinement study at `067e80a` and
later. Each result file records its own `source_commit`. The numerical kernels have not changed
since `8b70ed6`, and a regression at this tag reproduces the headline cases
(`benchmarks/regression.py`).

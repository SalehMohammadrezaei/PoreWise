# Materials examples: numerical verification and sensitivity

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
| compacted | mould | 10 | 1.5 | x | 1.9571 × 10<sup>-11</sup> | converged | 81000 | 2.1e-12 | 8.3e-06 | 3.7e-05 | 1.2e-04 |
| compacted | mould | 10 | 1.5 | y | 1.9571 × 10<sup>-11</sup> | converged | 81000 | 2.1e-12 | 8.3e-06 | 3.7e-05 | 1.2e-04 |
| compacted | mould | 10 | 1 | x | 1.9571 × 10<sup>-11</sup> | converged | 51500 | 2.6e-12 | 7.9e-06 | 7.5e-05 | 4.9e-04 |
| compacted | mould | 10 | 2.5 | x | 1.9571 × 10<sup>-11</sup> | converged | 122500 | 1.9e-12 | 1.0e-05 | 1.9e-05 | 3.1e-05 |
| compacted | mould | 20 | 1 | x | 1.7165 × 10<sup>-11</sup> | converged | 21500 | 1.2e-12 | 4.2e-06 | 1.7e-05 | 6.0e-05 |
| compacted | mould | 20 | 1 | y | 1.7165 × 10<sup>-11</sup> | converged | 21500 | 1.2e-12 | 4.2e-06 | 1.7e-05 | 6.0e-05 |
| compacted | mould | 20 | 1.5 | x | 1.7165 × 10<sup>-11</sup> | converged | 35500 | 1.2e-12 | 4.3e-06 | 8.6e-06 | 1.5e-05 |
| compacted | mould | 20 | 1.5 | y | 1.7165 × 10<sup>-11</sup> | converged | 35500 | 1.2e-12 | 4.3e-06 | 8.6e-06 | 1.5e-05 |
| compacted | mould | 5 | 1.5 | x | 1.9246 × 10<sup>-11</sup> | converged | 144000 | 3.9e-12 | 1.7e-05 | 1.5e-04 | 9.6e-04 |
| relaxed | mould | 10 | 1.5 | x | 4.6821 × 10<sup>-11</sup> | converged | 76000 | 2.8e-12 | 1.0e-05 | 7.2e-05 | 4.8e-04 |
| relaxed | mould | 10 | 1.5 | y | 4.6821 × 10<sup>-11</sup> | converged | 76000 | 2.8e-12 | 1.0e-05 | 7.2e-05 | 4.8e-04 |
| relaxed | mould | 20 | 1.5 | x | 5.7258 × 10<sup>-11</sup> | converged | 50500 | 1.9e-12 | 4.6e-06 | 1.9e-05 | 7.0e-05 |
| relaxed | mould | 20 | 1.5 | y | 5.7258 × 10<sup>-11</sup> | converged | 50500 | 1.9e-12 | 4.6e-06 | 1.9e-05 | 7.0e-05 |
| relaxed | mould | 5 | 1.5 | x | 5.0362 × 10<sup>-11</sup> | converged | 73000 | 2.7e-12 | 1.9e-05 | 3.1e-04 | 4.0e-03 |

## PU foam

| run | porosity in cylinder | k per sample area (m²) | k per grid area (m²) | termination | iterations | mass drift | plane-flux CV | Mach | Re<sub>pore</sub> | outer ring: area / flow |
|---|---|---|---|---|---|---|---|---|---|---|
| t0_coarse_porosity_buffer_crop180_tau1.5 | 0.9545 | 3.2200 × 10<sup>-7</sup> | 2.439 × 10<sup>-7</sup> | converged | 6000 | 1.9e-12 | 1.7e-05 | 3.3e-04 | 2.7e-03 | 18.5 % / 14.1 % |
| t0_coarse_porosity_buffer_tau1.5 | 0.9534 | 3.1313 × 10<sup>-7</sup> | 2.280 × 10<sup>-7</sup> | converged | 6500 | 2.9e-12 | 1.4e-05 | 3.0e-04 | 2.5e-03 | 18.5 % / 14.1 % |
| t0_coarse_porosity_mirror_tau1.5 | 0.9534 | 3.1802 × 10<sup>-7</sup> | 2.169 × 10<sup>-7</sup> | converged | 7500 | 1.1e-12 | 1.2e-05 | 2.9e-04 | 2.4e-03 | 18.5 % / 14.2 % |
| t0_coarse_porosity_periodic_tau1.5 | 0.9534 | 3.1064 × 10<sup>-7</sup> | 2.119 × 10<sup>-7</sup> | converged | 8500 | 4.6e-12 | 7.7e-06 | 2.8e-04 | 2.3e-03 | 18.5 % / 14.2 % |
| t0_coarse_quarter_periodic_tau1 | 0.9335 | 2.3927 × 10<sup>-7</sup> | 1.632 × 10<sup>-7</sup> | converged | 5000 | 4.5e-13 | 8.4e-06 | 5.2e-04 | 7.5e-03 | 18.5 % / 14.7 % |
| t0_coarse_quarter_periodic_tau1.5 | 0.9335 | 2.3927 × 10<sup>-7</sup> | 1.632 × 10<sup>-7</sup> | converged | 11000 | 7.8e-12 | 8.4e-06 | 2.6e-04 | 1.9e-03 | 18.5 % / 14.7 % |
| t0_full_buffer_crop200_tau1.5 | 0.9550 | 2.8926 × 10<sup>-7</sup> | 2.360 × 10<sup>-7</sup> | converged | 14000 | 3.4e-12 | 3.3e-05 | 1.3e-03 | 2.1e-02 | 18.5 % / 14.3 % |
| t1_coarse_porosity_buffer_tau1.5 | 0.9491 | 2.6717 × 10<sup>-7</sup> | 1.960 × 10<sup>-7</sup> | converged | 6500 | 6.0e-12 | 2.7e-05 | 2.7e-04 | 2.2e-03 | 18.5 % / 14.3 % |
| t2_coarse_porosity_buffer_tau1.5 | 0.9389 | 2.1420 × 10<sup>-7</sup> | 1.587 × 10<sup>-7</sup> | converged | 8000 | 2.7e-11 | 2.2e-05 | 2.5e-04 | 1.8e-03 | 18.5 % / 15.5 % |
| t3_coarse_porosity_buffer_tau1.5 | 0.9266 | 1.6891 × 10<sup>-7</sup> | 1.269 × 10<sup>-7</sup> | converged | 8000 | 2.2e-11 | 1.8e-05 | 2.1e-04 | 1.5e-03 | 18.5 % / 16.7 % |
| t4_coarse_porosity_buffer_tau1.5 | 0.9059 | 1.2376 × 10<sup>-7</sup> | 9.463 × 10<sup>-8</sup> | converged | 9500 | 5.9e-11 | 1.7e-05 | 1.7e-04 | 1.1e-03 | 18.5 % / 17.8 % |
| t5_coarse_porosity_buffer_tau1.5 | 0.8677 | 8.0149 × 10<sup>-8</sup> | 6.301 × 10<sup>-8</sup> | converged | 10500 | 3.4e-11 | 1.5e-05 | 1.4e-04 | 7.8e-04 | 18.5 % / 18.7 % |
| t6_coarse_porosity_buffer_tau1.5 | 0.9539 | 3.1500 × 10<sup>-7</sup> | 2.293 × 10<sup>-7</sup> | converged | 6500 | 3.7e-12 | 1.4e-05 | 3.1e-04 | 2.5e-03 | 18.5 % / 14.1 % |

The outer ring is the annulus one nominal pore wide at the sleeve. A ring carrying less flow than its share of the area
means the sleeve does not act as a bypass.

## Figures: provenance of the streamlines

All streamlines are integrated in the computed velocity fields with the wall-safe tracer (`examples/materials/common.py`,
the rule of `benchmarks/rock_streamlines.py`). A point is kept only if the eight cells around it are pore, and every
sub-step is checked. The counts are over every traced point, including lines not drawn.

| figure | field | points traced | points in a solid voxel | points with a solid cell in the stencil | lines drawn |
|---|---|---|---|---|---|
| pu_foam.png | t0_coarse_porosity_buffer_tau1.5.npz | 125253 | 0 | 0 | 202 |
| pu_foam_detail.png | t0_full_buffer_crop200_tau1.5.npz | 316402 | 0 | 0 | 75 |

`woven_preform.png` has no streamlines: it shows the speed and the in-plane direction on one plane of the computed field.

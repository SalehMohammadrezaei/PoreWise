"""Open-cell polyurethane foam (Ruf & Steeb, DaRUS-3010): data download and domain preparation.

Dataset: M. Ruf, H. Steeb (2024). In situ micro-XRCT data set of an open-cell polyurethane foam
sample under uniaxial compression load. DaRUS, https://doi.org/10.18419/DARUS-3010, CC BY 4.0.
Reticulated PU foam, 10 pores per inch, cylinder 52 mm x 50 mm, scanned at seven states t0..t6
under uniaxial displacements {0, 5, 10, 15, 20, 25, 0} mm (t6: unloaded again). The segmented
stacks (foam skeleton 255, void 0) have a documented voxel size of 74.8 um and were extracted by
the authors with a cylinder of 52.36 mm diameter.

    python examples/materials/pu_foam/prepare.py fetch  [--data DIR]      # segmented stacks + UTM data
    python examples/materials/pu_foam/prepare.py survey [--data DIR]      # porosity, length, solid volume

Domain conventions (z = sample axis = flow direction, (z, y, x) arrays):
    sleeve      everything outside the fitted sample cylinder is solid (a sealed core)
    coarsening  2x2x2 blocks; 'quarter': solid if >= 25 % of the block is skeleton (keeps thin
                struts connected, thickens them); 'porosity': the threshold that reproduces the
                full-resolution solid fraction inside the cylinder
    ends        'periodic' (the two end faces meet), 'buffer' (free-fluid reservoirs inside the
                sleeve at both ends, removed afterwards as a series resistance) or 'mirror'
                (sample + its mirror image along z, periodic by construction)
"""
import argparse, glob, hashlib, json, math, os, sys, urllib.request
from pathlib import Path

import numpy as np

VOXEL_M = 74.8e-6
STATES = {"t0": 0, "t1": 5, "t2": 10, "t3": 15, "t4": 20, "t5": 25, "t6": 0}   # imposed displacement, mm
DARUS = "https://darus.uni-stuttgart.de/api"
DOI = "doi:10.18419/DARUS-3010"
DEFAULT_DATA = Path(__file__).resolve().parents[3] / "results" / "materials" / "foam_data"


def fetch(data=DEFAULT_DATA):
    """Download the segmented stacks and the UTM records, checking the published MD5 sums."""
    import tarfile
    data = Path(data); data.mkdir(parents=True, exist_ok=True)
    meta = json.load(urllib.request.urlopen(f"{DARUS}/datasets/:persistentId/?persistentId={DOI}", timeout=60))
    files = meta["data"]["latestVersion"]["files"]
    (data / "dataset.json").write_text(json.dumps(meta, indent=1))
    for f in files:
        df = f["dataFile"]; name = df["filename"]
        if not (name.endswith("_segmented.tar.gz") or name.endswith(".tab") or name.endswith(".csv")):
            continue
        target = data / name
        if not target.exists():
            urllib.request.urlretrieve(f"{DARUS}/access/datafile/{df['id']}" + ("?format=original" if name.endswith(".tab") else ""), target)
        md5 = df.get("md5") or (df.get("checksum") or {}).get("value")
        if md5 and name.endswith(".tar.gz"):
            got = hashlib.md5(target.read_bytes()).hexdigest()
            if got != md5:
                raise RuntimeError(f"checksum mismatch for {name}: {got} != {md5}")
        if name.endswith("_segmented.tar.gz"):
            st = name.split("_")[0]
            if not (data / st).exists():
                with tarfile.open(target) as t:
                    t.extractall(data / st)
    print("data in", data)


def load(state, data=DEFAULT_DATA):
    """Skeleton (True) of one state, (z, y, x), from the segmented TIFF stack."""
    from PIL import Image
    files = sorted(glob.glob(str(Path(data) / state / f"{state}_segmented" / "*.tif")))
    if not files:
        raise FileNotFoundError(f"no segmented stack for {state} under {data}; run 'prepare.py fetch'")
    return np.stack([np.asarray(Image.open(f)) > 0 for f in files])


def cylinder(skel):
    """Sample cylinder fitted to the skeleton's projected extent: centre (y, x), radius (voxels)."""
    proj = skel.any(axis=0); ys, xs = np.nonzero(proj)
    return ((ys.min() + ys.max()) / 2, (xs.min() + xs.max()) / 2), (xs.max() - xs.min() + 1) / 2


def disc_mask(shape_yx, centre, radius):
    yy, xx = np.mgrid[:shape_yx[0], :shape_yx[1]]
    return (yy - centre[0]) ** 2 + (xx - centre[1]) ** 2 <= radius ** 2


def coarsen(skel, factor=2, rule="porosity", inside=None):
    """Block coarsening by `factor`. Returns (coarse skeleton, threshold used)."""
    nz, ny, nx = (n // factor * factor for n in skel.shape)
    if rule == "porosity":
        # block fractions of a binary image take only factor^3 + 1 values, so a threshold cannot hit
        # a target porosity when many blocks tie; a light Gaussian (sigma = factor / 2 voxels) before
        # block averaging makes them continuous (anti-aliasing), the threshold then matches exactly
        from scipy import ndimage
        src = ndimage.gaussian_filter(skel[:nz, :ny, :nx].astype(np.float32), factor / 2)
    else:
        src = skel[:nz, :ny, :nx]
    frac = src.reshape(nz // factor, factor, ny // factor, factor, nx // factor, factor).mean(axis=(1, 3, 5))
    if rule == "quarter":
        thr = 0.25
    elif rule == "porosity":
        target = skel[:, inside].mean() if inside is not None else skel.mean()
        ins_c = inside[:ny:factor, :nx:factor] if inside is not None else np.ones(frac.shape[1:], bool)
        vals = np.sort(frac[:, ins_c].ravel())[::-1]
        k = int(round(target * vals.size)); thr = float(vals[max(k - 1, 0)]) if k else 1.0
    else:
        raise ValueError(rule)
    return frac >= thr, float(thr)


def build_domain(skel, centre, radius, *, ends="buffer", buffer_voxels=0, crop=None):
    """Sealed-core domain: solid sleeve outside the cylinder, chosen end treatment.

    Returns (solid, info) with the sample's share of the cross-section and the lengths needed to
    normalise the permeability.
    """
    if crop:
        z0 = (skel.shape[0] - crop) // 2; skel = skel[z0:z0 + crop]
    inside = disc_mask(skel.shape[1:], centre, radius)
    sample = skel | ~inside[None]
    if ends == "mirror":
        solid = np.concatenate([sample, sample[::-1]], axis=0)
    elif ends == "buffer":
        pad = np.broadcast_to(~inside[None], (buffer_voxels,) + inside.shape)
        solid = np.concatenate([pad, sample, pad], axis=0)
    elif ends == "periodic":
        solid = sample
    else:
        raise ValueError(ends)
    info = dict(ends=ends, crop_slices=crop, sample_slices=int(skel.shape[0]), domain_slices=int(solid.shape[0]),
                buffer_voxels=int(buffer_voxels if ends == "buffer" else 0), sample_area_fraction=float(inside.mean()),
                porosity_in_cylinder=float(1 - skel[:, inside].mean()), radius_voxels=float(radius))
    return np.ascontiguousarray(solid), info


def survey(data=DEFAULT_DATA):
    rows = []
    for st, disp in STATES.items():
        sk = load(st, data); c, r = cylinder(sk); ins = disc_mask(sk.shape[1:], c, r)
        phi_s = float(sk[:, ins].mean()); L = sk.shape[0] * VOXEL_M * 1e3; area = math.pi * (r * VOXEL_M * 1e3) ** 2
        rows.append(dict(state=st, displacement_mm=disp, slices=int(sk.shape[0]), length_mm=L, diameter_mm=2 * r * VOXEL_M * 1e3,
                         porosity=1 - phi_s, solid_volume_mm3=phi_s * area * L))
    v0 = rows[0]["solid_volume_mm3"]; L0 = rows[0]["length_mm"]
    for r in rows:
        r["solid_volume_rel_t0"] = r["solid_volume_mm3"] / v0; r["axial_strain_vs_t0"] = 1 - r["length_mm"] / L0
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("action", choices=["fetch", "survey"]); ap.add_argument("--data", default=str(DEFAULT_DATA))
    a = ap.parse_args()
    if a.action == "fetch":
        fetch(a.data)
    else:
        rows = survey(a.data)
        out = Path(__file__).resolve().parent / "results" / "survey.json"; out.parent.mkdir(exist_ok=True)
        out.write_text(json.dumps(rows, indent=1)); print(json.dumps(rows, indent=1))

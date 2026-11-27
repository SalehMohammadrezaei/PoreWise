"""Idealised plain-weave preform: geometry generation and bookkeeping.

The construction follows the usual mesoscale textile model (as in TexGen): each tow is a
constant cross-section swept along a crimped centreline.

    centreline      z_c(u) = z0 +/- A g(cos(pi u / s)), A = h/2, s = yarn spacing
                    g(c) = tanh(beta c) / tanh(beta)   (beta = 0: plain sinusoid); beta > 0
                    flattens the path over the crossings, where real tows lie on each other
    cross-section   elliptical, width w, thickness h: half-thickness t(d) = h/2 sqrt(1-(2d/w)^2)
    interlacing     warp j (along x, centred at y = j s) and weft i (along y, centred at x = i s)
                    have opposite signs at every crossing (i s, j s): the plain-weave pattern

Interference. An elliptical tow swept along a crimped path intersects the tow it crosses
beside the crossing point unless the path curvature is below the section curvature (for a
sinusoid with A = h/2 that needs w < 2s/pi). As TexGen does with its interference correction,
the overlap is removed: a point claimed by both tows belongs to the tow in which it lies deeper
(larger margin / half-thickness), giving flat contact faces. The solid union is unchanged; the
correction only makes tow volumes and tow surfaces physical. The overlap volume before
correction is reported.

Plies. Ply 2 is shifted in plane by `shift` x s and lowered onto ply 1 until the two just touch
(nesting by contact, column by column). Mould plates are placed on the outermost crests, so no
tow-free slot is left between preform and mould.

Fibre content. 3K carbon tows (3000 filaments, 7 um, 1.78 g/cm3): the fibre area of a tow is
fixed; the elliptical section area pi w h / 4 then sets the intra-tow fibre fraction, and the
cavity thickness sets the fibre volume fraction of the preform. Both are outputs, reported
separately from the tow volume fraction.
"""
from dataclasses import asdict, dataclass
import math

import numpy as np


@dataclass
class Fibre:
    diameter_um: float = 7.0
    filaments: int = 3000
    density_g_cm3: float = 1.78

    @property
    def tow_area_mm2(self):
        return self.filaments * math.pi * (self.diameter_um * 1e-3 / 2) ** 2


@dataclass
class Weave:
    name: str
    spacing_mm: float           # s: yarn spacing (ends and picks per mm = 1/s)
    width_mm: float             # w: tow width
    intra_tow_vf: float         # fibre fraction inside a tow; sets the tow thickness
    beta: float = 1.0           # path flattening over the crossings (0 = sinusoid)
    plies: int = 2
    shift: tuple = (0.5, 0.5)   # in-plane shift of ply 2, in yarn spacings
    fibre: Fibre = None

    def __post_init__(self):
        self.fibre = self.fibre or Fibre()

    @property
    def thickness_mm(self):     # h from the section area pi w h / 4 = fibre area / intra-tow Vf
        return 4 * self.fibre.tow_area_mm2 / (self.intra_tow_vf * math.pi * self.width_mm)

    @property
    def amplitude_mm(self):
        return self.thickness_mm / 2

    @property
    def cell_mm(self):          # plain-weave unit cell: two yarns each way
        return 2 * self.spacing_mm

    def describe(self):
        d = asdict(self); d.update(thickness_mm=self.thickness_mm, amplitude_mm=self.amplitude_mm, cell_mm=self.cell_mm,
                                   tow_fibre_area_mm2=self.fibre.tow_area_mm2, gap_mm=self.spacing_mm - self.width_mm)
        return d


def _path(u, wv):
    c = np.cos(np.pi * u / wv.spacing_mm)
    return c if wv.beta <= 0 else np.tanh(wv.beta * c) / math.tanh(wv.beta)


def _half_thickness(d, wv):
    return wv.thickness_mm / 2 * np.sqrt(np.clip(1 - (2 * d / wv.width_mm) ** 2, 0, None))


def tow_margins(x, y, z, wv, z0, shift):
    """Margin t - |z - z_c| (> 0 inside) and half-thickness t, for the warp and weft of one ply.

    Coordinates in mm, periodic in x and y with period 2 s. Returns dict of (margin, t) arrays.
    """
    s, L = wv.spacing_mm, wv.cell_mm
    xs, ys = (x - shift[0] * s) % L, (y - shift[1] * s) % L
    shape = np.broadcast(x, y, z).shape
    out = {}
    for kind in ("warp", "weft"):
        m = np.full(shape, -np.inf); tt = np.zeros(shape)
        for j in range(2):
            if kind == "warp":   # along x, centred at y = j s, sign (-1)^j at x = 0
                d = (ys - j * s + L / 2) % L - L / 2; zc = z0 + wv.amplitude_mm * _path(xs, wv) * (-1) ** j
            else:                # along y, centred at x = j s, opposite sign at every crossing
                d = (xs - j * s + L / 2) % L - L / 2; zc = z0 - wv.amplitude_mm * _path(ys, wv) * (-1) ** j
            t = _half_thickness(d, wv)
            mj = np.where(t > 0, t - np.abs(z - zc), -np.inf)
            better = mj > m
            m = np.where(better, mj, m); tt = np.where(better, t, tt)
        out[kind] = (m, tt)
    return out


def corrected_inside(margins):
    """Interference-corrected warp / weft membership: shared points go to the deeper tow."""
    (mw, tw), (mf, tf) = margins["warp"], margins["weft"]
    iw, if_ = mw > 0, mf > 0
    both = iw & if_
    with np.errstate(invalid="ignore", divide="ignore"):
        deeper_w = (mw / np.where(tw > 0, tw, 1)) >= (mf / np.where(tf > 0, tf, 1))
    warp = iw & (~both | deeper_w); weft = if_ & (~both | ~deeper_w)
    return warp, weft, both


def _column_extent(wv, z0, shift, n=1200):
    """Top and bottom surface of one ply on an n x n column grid over the unit cell (mm)."""
    g = (np.arange(n) + 0.5) * wv.cell_mm / n
    X, Y = np.meshgrid(g, g, indexing="xy")
    top = np.full(X.shape, -np.inf); bot = np.full(X.shape, np.inf)
    s, L = wv.spacing_mm, wv.cell_mm
    xs, ys = (X - shift[0] * s) % L, (Y - shift[1] * s) % L
    for j in range(2):
        for kind in ("warp", "weft"):
            if kind == "warp":
                d = (ys - j * s + L / 2) % L - L / 2; zc = z0 + wv.amplitude_mm * _path(xs, wv) * (-1) ** j
            else:
                d = (xs - j * s + L / 2) % L - L / 2; zc = z0 - wv.amplitude_mm * _path(ys, wv) * (-1) ** j
            t = _half_thickness(d, wv); has = t > 0
            top = np.where(has, np.maximum(top, zc + t), top); bot = np.where(has, np.minimum(bot, zc - t), bot)
    return top, bot


def stack_layout(wv, n=1200):
    """Ply offsets (nesting by contact) and mould positions, in mm. Ply 1 mid-plane at z = 0."""
    z0 = [0.0]; shifts = [(0.0, 0.0)]
    top_prev, bot_first = _column_extent(wv, 0.0, (0.0, 0.0), n)
    for p in range(1, wv.plies):
        sh = (wv.shift[0] * p % 2, wv.shift[1] * p % 2)
        top_rel, bot_rel = _column_extent(wv, 0.0, sh, n)
        dz = float(np.nanmax(np.where(np.isfinite(top_prev) & np.isfinite(bot_rel), top_prev - bot_rel, -np.inf)))
        z0.append(dz); shifts.append(sh)
        top_prev = np.maximum(top_prev, top_rel + dz)
    bottom = float(np.nanmin(np.where(np.isfinite(bot_first), bot_first, np.inf)))
    top = float(np.nanmax(np.where(np.isfinite(top_prev), top_prev, -np.inf)))
    return dict(ply_z0_mm=z0, ply_shift=shifts, mould_bottom_mm=bottom, mould_top_mm=top, cavity_mm=top - bottom)


def voxelise(wv, dx_mm, layout=None):
    """Voxel image (z, y, x) with one mould-plate layer at the bottom and top.

    Cell centres sample the corrected tow shapes. Returns (solid, labels, info); labels:
    0 pore, 1..2*plies tows (warp of ply p = 2p-1, weft = 2p), 255 mould.
    """
    layout = layout or stack_layout(wv)
    nxy = int(round(wv.cell_mm / dx_mm))
    if abs(nxy * dx_mm - wv.cell_mm) > 1e-9 * wv.cell_mm:
        raise ValueError(f"dx must divide the unit cell ({wv.cell_mm} mm)")
    # cubic voxels: the cavity is the nearest whole number of layers, centred on the preform, so
    # the mould sits within half a voxel of the outermost crests (clearance or contact < dx/2)
    nz = max(1, int(round(layout["cavity_mm"] / dx_mm)))
    z_start = layout["mould_bottom_mm"] + (layout["cavity_mm"] - nz * dx_mm) / 2
    xc = (np.arange(nxy) + 0.5) * dx_mm
    zc = z_start + (np.arange(nz) + 0.5) * dx_mm
    Z, Y, X = np.meshgrid(zc, xc, xc, indexing="ij")
    labels = np.zeros((nz + 2, nxy, nxy), np.uint8); overlap = 0
    for p, (z0, sh) in enumerate(zip(layout["ply_z0_mm"], layout["ply_shift"])):
        mg = tow_margins(X, Y, Z, wv, z0, sh)
        warp, weft, both = corrected_inside(mg)
        overlap += int(both.sum())
        lab = labels[1:-1]
        lab[warp & (lab == 0)] = 2 * p + 1; lab[weft & (lab == 0)] = 2 * p + 2
    labels[0] = labels[-1] = 255
    solid = labels > 0
    cav = solid[1:-1]
    tow_cells = int(cav.sum())
    info = dict(dx_mm=dx_mm, grid_zyx=list(solid.shape), cavity_layers=nz, cavity_mm=layout["cavity_mm"],
                cavity_voxelised_mm=nz * dx_mm, mould_offset_mm=(nz * dx_mm - layout["cavity_mm"]) / 2,
                tow_volume_fraction=tow_cells / cav.size,
                overlap_before_correction_fraction=overlap / max(tow_cells, 1),
                tow_free_layers=int((~cav.reshape(nz, -1).any(axis=1)).sum()), layout=layout)
    info.update(fibre_bookkeeping(wv, layout))
    return solid, labels, info


def periodic_stack_layout(wv, n=1200):
    """Two plies repeating through the thickness (an infinite laminate, no mould).

    Ply 2 rests on ply 1 (offset d12, as in stack_layout) and the next copy of ply 1 rests on
    ply 2 (offset d21, found the same way), so the period is d12 + d21.
    """
    lay = stack_layout(wv, n)
    sh2 = lay["ply_shift"][1]
    top2, _ = _column_extent(wv, 0.0, sh2, n); _, bot1 = _column_extent(wv, 0.0, (0.0, 0.0), n)
    d21 = float(np.nanmax(np.where(np.isfinite(top2) & np.isfinite(bot1), top2 - bot1, -np.inf)))
    return dict(ply_z0_mm=lay["ply_z0_mm"], ply_shift=lay["ply_shift"], period_mm=lay["ply_z0_mm"][1] + d21, d12_mm=lay["ply_z0_mm"][1], d21_mm=d21)


def voxelise_periodic(wv, dx_mm, layout=None):
    """Voxel image of the periodic two-ply stack (periodic in x, y and z; no mould plates)."""
    layout = layout or periodic_stack_layout(wv)
    # cubic voxels: the period is rounded DOWN to whole voxels, which can only press the plies into
    # contact by less than a voxel (never opens a slot between them)
    nxy = int(round(wv.cell_mm / dx_mm)); nz = max(1, int(math.floor(layout["period_mm"] / dx_mm + 1e-9))); P = nz * dx_mm
    xc = (np.arange(nxy) + 0.5) * dx_mm; zc = -wv.thickness_mm + (np.arange(nz) + 0.5) * dx_mm
    Z, Y, X = np.meshgrid(zc, xc, xc, indexing="ij")
    labels = np.zeros((nz, nxy, nxy), np.uint8)
    for p, (z0, sh) in enumerate(zip(layout["ply_z0_mm"], layout["ply_shift"])):
        for img in (-P, 0.0, P):                                           # periodic images through the thickness
            warp, weft, _ = corrected_inside(tow_margins(X, Y, Z, wv, z0 + img, sh))
            labels[warp & (labels == 0)] = 2 * p + 1; labels[weft & (labels == 0)] = 2 * p + 2
    solid = labels > 0
    fb = fibre_bookkeeping(wv, dict(cavity_mm=P))
    info = dict(dx_mm=dx_mm, grid_zyx=list(solid.shape), period_mm=P, period_exact_mm=layout["period_mm"], tow_volume_fraction=float(solid.mean()),
                layout=layout, **fb)
    return solid, labels, info


def fibre_bookkeeping(wv, layout, n=4000):
    """Fibre volume fraction of the cavity, implied areal weight, crimp."""
    u = (np.arange(n) + 0.5) * wv.cell_mm / n
    dz = np.gradient(wv.amplitude_mm * _path(u, wv), u)
    crimp = float(np.mean(np.sqrt(1 + dz ** 2)))                           # path length / projected length
    fibre_thickness = 2 * wv.fibre.tow_area_mm2 / wv.spacing_mm * crimp    # mm of fibre per ply (both directions)
    return dict(crimp_factor=crimp, fibre_thickness_per_ply_mm=fibre_thickness,
                areal_weight_g_m2=fibre_thickness * 1e-3 * wv.fibre.density_g_cm3 * 1e6,
                fibre_volume_fraction=wv.plies * fibre_thickness / layout["cavity_mm"],
                intra_tow_vf=wv.intra_tow_vf)


def plain_weave_alternates(wv):
    """True if warp and weft centrelines are on opposite sides at all four crossings of the cell."""
    s = wv.spacing_mm; ok = True
    for i in range(2):
        for j in range(2):
            zw = wv.amplitude_mm * _path(np.array(i * s), wv) * (-1) ** j
            zf = -wv.amplitude_mm * _path(np.array(j * s), wv) * (-1) ** i
            ok &= bool(np.sign(zw) == -np.sign(zf) and zw != 0)
    return ok


STATES = {
    # Two compaction states of the same fabric (3K carbon, 5 ends/cm): compaction flattens and
    # widens the tows (higher intra-tow fibre fraction, smaller gaps) and the plies nest closer.
    "relaxed": Weave("relaxed", spacing_mm=2.0, width_mm=1.70, intra_tow_vf=0.62, beta=1.5),
    "compacted": Weave("compacted", spacing_mm=2.0, width_mm=1.85, intra_tow_vf=0.72, beta=1.5),
}

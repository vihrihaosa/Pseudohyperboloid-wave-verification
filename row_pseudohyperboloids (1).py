# -*- coding: utf-8 -*-
"""
FINAL SCRIPT.

Construction of open n-th order pseudohyperboloids by the recursive scheme

    B_2 = {d},
    B_{k+1} = {R_k + r, R_k - r : r in B_k},
    N_n = 2^(n-2).

IN THIS VERSION:
1) the script builds not just a set of branches, but the FINAL UNION VOLUME;
2) if individual pseudohyperboloidal toroids intersect or are nested,
   only their common part is removed as duplicated volume;
3) the toroids themselves are NOT removed as generating objects;
4) the final figure is the union of all volumes;
5) if the toroids do not intersect, they are shown separately;
6) the rule applies to all orders and all rows.

Supported:
- vertical type;
- horizontal type;
- all orders from 2 to n;
- 2D meridional sections of the union volume;
- 3D surfaces of the union volume of revolution;
- a stack of identical copies on a common rotation axis.
"""

from __future__ import annotations

import math
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# DISPLAY SETTINGS
# ------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 300,
    "font.size": 11,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "font.family": "DejaVu Sans",
})

SURFACE_ALPHA = 0.24
MERGE_TOL = 1.0e-9
NPHI = 220
# By default, 2D areas are NOT filled so that
# artificial straight closing segments from fill_betweenx do not appear.
FILL_2D_AREAS = False


# ------------------------------------------------------------
# BASE GEOMETRY
# ------------------------------------------------------------
def c_focus(a: float, b: float) -> float:
    return math.sqrt(a * a + b * b)


def L_open(a: float, b: float, R: float) -> float:
    """Longitudinal existence limit of the open vertical-type profile."""
    return a * math.sqrt(1.0 + (R / b) ** 2)


def eta_hyperbola(abs_x: np.ndarray, a: float, b: float) -> np.ndarray:
    abs_x = np.asarray(abs_x, dtype=float)
    val = (abs_x / a) ** 2 - 1.0
    val = np.maximum(val, 0.0)
    return b * np.sqrt(val)


def rho_vertical(abs_s: np.ndarray, a: float, b: float, R: float) -> np.ndarray:
    """
    Open vertical-type profile:
        rho(|s|) = R - b*sqrt((|s|/a)^2 - 1),   a <= |s| <= L.
    """
    out = R - eta_hyperbola(abs_s, a, b)
    return np.maximum(out, 0.0)


def x_horizontal(u: np.ndarray, a: float, b: float, R: float) -> np.ndarray:
    """
    Horizontal type:
        x_h(u) = a * sqrt(1 + ((R - |u|)/b)^2),   |u| <= R.
    """
    u = np.asarray(u, dtype=float)
    return a * np.sqrt(1.0 + ((R - np.abs(u)) / b) ** 2)


# ------------------------------------------------------------
# BASE GRIDS FOR THE TWO TYPES
# ------------------------------------------------------------
def _ensure_zero_in_axis(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    if np.any(np.isclose(arr, 0.0, atol=1e-14)):
        return arr
    arr2 = np.concatenate([arr, np.array([0.0])])
    arr2.sort()
    return arr2


def vertical_base_grid(a: float, b: float, R: float, npts: int = 900) -> Dict[str, np.ndarray]:
    L = L_open(a, b, R)
    s_left = np.linspace(-L, -a, npts)
    s_right = np.linspace(a, L, npts)
    d_left = rho_vertical(np.abs(s_left), a, b, R)
    d_right = rho_vertical(np.abs(s_right), a, b, R)
    return {
        "axis_left": s_left,
        "axis_right": s_right,
        "d_left": d_left,
        "d_right": d_right,
        "L": L,
    }


def horizontal_base_grid(a: float, b: float, R: float, npts: int = 900) -> Dict[str, np.ndarray]:
    u = np.linspace(-R, R, npts)
    u = _ensure_zero_in_axis(u)
    d = x_horizontal(u, a, b, R)
    return {
        "axis": u,
        "d": d,
        "L": R,
    }


# ------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------
def morphology_class(R: float, offsets: Sequence[float]) -> str:
    if not offsets:
        return "base 2nd order only"
    R_star = max(offsets)
    if R > R_star:
        return "with internal intersection"
    if math.isclose(R, R_star, rel_tol=1e-12, abs_tol=1e-12):
        return "symmetric"
    return "non-intersecting annular"


def validate_offsets(base_distance: np.ndarray, offsets: Sequence[float]) -> List[str]:
    """
    Non-blocking check. The construction can still be performed.
    Messages are informational.
    """
    msgs = []
    current_max = float(np.max(base_distance))
    for i, offset in enumerate(offsets, start=1):
        if offset < current_max:
            msgs.append(
                f"At step R{i}: R{i} = {offset:.6f} < {current_max:.6f}. "
                f"This means intersection / nesting of components and requires taking the union of volumes."
            )
        elif math.isclose(offset, current_max, rel_tol=1e-12, abs_tol=1e-12):
            msgs.append(
                f"At step R{i}: R{i} = {offset:.6f} equals {current_max:.6f}. "
                f"This is the limiting case of component contact."
            )
        current_max = offset + current_max
    return msgs


def stack_shifts(axis_half_length: float, h: float, m: int) -> np.ndarray:
    """
    Returns the center shifts of m identical copies along the common rotation axis.

    The full length of one copy along the common axis is 2 * axis_half_length.
    Step between centers of neighboring copies:
        step = 2 * axis_half_length + h
    """
    if int(m) != m or m < 1:
        raise ValueError("Parameter m must be an integer >= 1")
    m = int(m)
    step = 2.0 * float(axis_half_length) + float(h)
    return -np.arange(m, dtype=float) * step


def _ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def _save_or_show(fig, path: Path | None, show: bool):
    if path is not None:
        fig.savefig(path, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


# ------------------------------------------------------------
# RECURSIVE RADIAL INTERVALS
# ------------------------------------------------------------
def build_recursive_interval_arrays(base_distance: np.ndarray,
                                    offsets: Sequence[float]) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Builds radial-volume intervals for one order.

    Order 2:
        [0, d]

    Next step for each interval [lo, hi]:
        [R - hi, R - lo] and [R + lo, R + hi]

    This corresponds to the union of the volumes of the generating toroids,
    without deleting the toroids themselves as construction objects.
    """
    zero = np.zeros_like(base_distance, dtype=float)
    intervals = [(zero.copy(), np.asarray(base_distance, dtype=float).copy())]
    for Rk in offsets:
        nxt: List[Tuple[np.ndarray, np.ndarray]] = []
        for lo, hi in intervals:
            a = np.maximum(Rk - hi, 0.0)
            b = np.maximum(Rk - lo, 0.0)
            c = np.maximum(Rk + lo, 0.0)
            d = np.maximum(Rk + hi, 0.0)
            nxt.append((a, b))
            nxt.append((c, d))
        intervals = nxt
    return intervals


@dataclass
class Piece:
    axis: np.ndarray
    lo: np.ndarray
    hi: np.ndarray


# ------------------------------------------------------------
# BUILDING THE UNION VOLUME FOR A GIVEN ORDER
# ------------------------------------------------------------
def _interp_piece_to_global_axis(piece: Piece, global_axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = piece.axis
    lo = piece.lo
    hi = piece.hi
    lo_i = np.interp(global_axis, x, lo)
    hi_i = np.interp(global_axis, x, hi)
    mask = (global_axis >= x.min() - 1e-12) & (global_axis <= x.max() + 1e-12)
    lo_i[~mask] = np.nan
    hi_i[~mask] = np.nan
    return lo_i, hi_i


def _merge_scalar_intervals(intervals: List[Tuple[float, float]], tol: float = MERGE_TOL) -> List[Tuple[float, float]]:
    if not intervals:
        return []
    data = [(float(min(a, b)), float(max(a, b))) for a, b in intervals if np.isfinite(a) and np.isfinite(b) and b > a + tol]
    if not data:
        return []
    data.sort(key=lambda t: (t[0], t[1]))
    merged = [list(data[0])]
    for a, b in data[1:]:
        if a <= merged[-1][1] + tol:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


def _build_merged_components_from_pieces(pieces: List[Piece], n_global: int | None = None) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    if not pieces:
        raise ValueError("No components available for building the union volume")

    # Global axis: union of all real nodes.
    all_axis = np.concatenate([p.axis for p in pieces])
    global_axis = np.unique(np.round(all_axis, 12))
    global_axis.sort()

    if n_global is not None and n_global > global_axis.size:
        dense = np.linspace(global_axis.min(), global_axis.max(), n_global)
        global_axis = np.unique(np.concatenate([global_axis, dense]))
        global_axis.sort()

    # IMPORTANT: if there is a real empty gap between two rows / components,
    # lines or surfaces must not connect the edges of this gap.
    # Therefore, control points are added at the centers of all empty gaps
    # between component projections onto the common axis. At these points there are no intervals,
    # and the plot breaks by itself without adding anything.
    cover = sorted((float(np.min(p.axis)), float(np.max(p.axis))) for p in pieces)
    merged_cover: List[List[float]] = []
    for a0, b0 in cover:
        if not merged_cover or a0 > merged_cover[-1][1] + MERGE_TOL:
            merged_cover.append([a0, b0])
        else:
            merged_cover[-1][1] = max(merged_cover[-1][1], b0)
    gap_points = []
    for left, right in zip(merged_cover[:-1], merged_cover[1:]):
        if right[0] > left[1] + MERGE_TOL:
            gap_points.append(0.5 * (left[1] + right[0]))
    if gap_points:
        global_axis = np.unique(np.concatenate([global_axis, np.array(gap_points, dtype=float)]))
        global_axis.sort()

    interpolated: List[Tuple[np.ndarray, np.ndarray]] = [
        _interp_piece_to_global_axis(piece, global_axis) for piece in pieces
    ]

    merged_per_sample: List[List[Tuple[float, float]]] = []
    max_count = 0
    for j in range(global_axis.size):
        scalars: List[Tuple[float, float]] = []
        for lo_i, hi_i in interpolated:
            lo = lo_i[j]
            hi = hi_i[j]
            if np.isfinite(lo) and np.isfinite(hi) and hi > lo + MERGE_TOL:
                scalars.append((max(0.0, lo), max(0.0, hi)))
        merged = _merge_scalar_intervals(scalars, tol=MERGE_TOL)
        merged_per_sample.append(merged)
        max_count = max(max_count, len(merged))

    lo_components = [np.full(global_axis.shape, np.nan, dtype=float) for _ in range(max_count)]
    hi_components = [np.full(global_axis.shape, np.nan, dtype=float) for _ in range(max_count)]

    for j, merged in enumerate(merged_per_sample):
        for k, (a, b) in enumerate(merged):
            lo_components[k][j] = a
            hi_components[k][j] = b

    return global_axis, lo_components, hi_components


def _contiguous_segments(mask: np.ndarray, axis: np.ndarray | None = None) -> List[Tuple[int, int]]:
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0:
        return []

    # Additional protection: if two neighboring active points are separated
    # by too large an axial step, this indicates a real empty gap
    # between components. It must not be connected by a straight line.
    if axis is not None and axis.size == mask.size and axis.size >= 3:
        dif = np.diff(axis)
        pos = dif[np.isfinite(dif) & (dif > 0)]
        if pos.size:
            med = float(np.median(pos))
            jump_thr = max(1e-12, 4.0 * med)
            jump_breaks = np.where(dif > jump_thr)[0]
            for j in jump_breaks:
                if 0 <= j < mask.size - 1:
                    mask[j] = False
                    mask[j + 1] = False

    segments: List[Tuple[int, int]] = []
    start = None
    for i, ok in enumerate(mask):
        if ok and start is None:
            start = i
        if (not ok) and start is not None:
            if i - start >= 2:
                segments.append((start, i))
            start = None
    if start is not None and mask.size - start >= 2:
        segments.append((start, mask.size))
    return segments


def _split_on_artificial_jumps(radius: np.ndarray) -> List[Tuple[int, int]]:
    """
    Additionally splits the curve where the point sequence makes
    a non-geometric jump due to a change in the merged-interval index.
    This removes artificial straight connecting segments.
    """
    radius = np.asarray(radius, dtype=float)
    n = radius.size
    if n < 2:
        return [(0, n)] if n > 0 else []
    dr = np.abs(np.diff(radius))
    finite = np.isfinite(dr)
    valid = dr[finite]
    if valid.size == 0:
        return [(0, n)]
    med = float(np.median(valid))
    mx = float(np.max(valid))
    if med <= 1e-12:
        jump_thr = max(1e-8, 0.05 * mx)
    else:
        jump_thr = max(1e-8, 8.0 * med)
    break_ids = np.where(dr > jump_thr)[0] + 1
    if break_ids.size == 0:
        return [(0, n)]
    segments = []
    st = 0
    for br in break_ids:
        if br - st >= 2:
            segments.append((st, int(br)))
        st = int(br)
    if n - st >= 2:
        segments.append((st, n))
    return segments


def _plot_boundary_curve(ax, radius: np.ndarray, axis: np.ndarray, lw: float):
    for st, en in _split_on_artificial_jumps(radius):
        ax.plot(radius[st:en], axis[st:en], lw=lw)




# ------------------------------------------------------------
# REFERENCE LINES, GENERATRICES AND GEOMETRIC FOCI
# ------------------------------------------------------------
def _desc_add_const(C: float, desc):
    if desc[0] == "const":
        return ("const", float(C) + float(desc[1]))
    return ("hyp", float(C) + float(desc[1]), float(desc[2]))


def _desc_sub_const(C: float, desc):
    if desc[0] == "const":
        return ("const", float(C) - float(desc[1]))
    # C - (A + sigma*d) = (C-A) - sigma*d
    return ("hyp", float(C) - float(desc[1]), -float(desc[2]))


def build_recursive_boundary_descriptors(offsets: Sequence[float]):
    """
    Exact descriptors of generating hyperbolic boundaries.

    A boundary of the form ("hyp", C, sigma) means:
        r(x) = C + sigma*d(x),
    where d(x) is the base radial generatrix.

    Merge creates neither new hyperbolas nor new foci.
    """
    intervals = [(('const', 0.0), ('hyp', 0.0, +1.0))]
    for Rk in offsets:
        nxt = []
        for lo, hi in intervals:
            nxt.append((_desc_sub_const(Rk, hi), _desc_sub_const(Rk, lo)))
            nxt.append((_desc_add_const(Rk, lo), _desc_add_const(Rk, hi)))
        intervals = nxt
    return intervals


def _unique_hyperbolic_descriptors(offsets: Sequence[float]):
    intervals = build_recursive_boundary_descriptors(offsets)
    seen = set()
    out = []
    for lo, hi in intervals:
        for desc in (lo, hi):
            if desc[0] != 'hyp':
                continue
            key = (round(float(desc[1]), 10), round(float(desc[2]), 10))
            if key not in seen:
                seen.add(key)
                out.append(desc)
    return out


def _curve_matches_visible_boundary(global_axis: np.ndarray,
                                    curve_axis: np.ndarray,
                                    curve_r: np.ndarray,
                                    lo_components: List[np.ndarray],
                                    hi_components: List[np.ndarray],
                                    tol: float = 1.0e-4) -> bool:
    """Checks whether the exact generating curve is visible on the 2D Merge boundary."""
    ga = np.asarray(global_axis, dtype=float)
    xa = np.asarray(curve_axis, dtype=float)
    rr = np.asarray(curve_r, dtype=float)
    if xa.size < 8:
        return False

    r_glob = np.interp(ga, xa, rr, left=np.nan, right=np.nan)
    in_range = (ga >= xa.min() - 1e-12) & (ga <= xa.max() + 1e-12)
    r_glob[~in_range] = np.nan

    visible_arrays: List[np.ndarray] = []
    for arr in hi_components:
        visible_arrays.append(np.asarray(arr, dtype=float))
    for arr in lo_components:
        if np.any(np.isfinite(arr) & (arr > MERGE_TOL)):
            visible_arrays.append(np.asarray(arr, dtype=float))

    scale = max(1.0, float(np.nanmax(np.abs(rr))) if np.any(np.isfinite(rr)) else 1.0)
    atol = max(tol, 2.0e-3 * scale)

    for vis in visible_arrays:
        mask = np.isfinite(r_glob) & np.isfinite(vis) & (r_glob > MERGE_TOL) & (vis > MERGE_TOL)
        if np.count_nonzero(mask) < 6:
            continue
        close = np.zeros_like(mask, dtype=bool)
        close[mask] = np.abs(r_glob[mask] - vis[mask]) <= atol
        for st, en in _contiguous_segments(close, ga):
            if en - st >= 6:
                span = float(ga[en - 1] - ga[st])
                if span > 0.02 * max(1.0, float(xa.max() - xa.min())):
                    return True
    return False


def _plot_focus_points_2d(ax, points: List[Tuple[float, float]], scale_ref: float):
    if not points:
        return
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    # The white marker fill reduces the false visual impression
    # that the generatrix line passes directly through a focus.
    ax.scatter(xs, ys, s=28, marker='o', facecolors='white', edgecolors='darkgreen',
               linewidths=1.25, zorder=8,
               label='Geometric foci of generating hyperbolas')
    dx = 0.018 * max(scale_ref, 1.0)
    for x, y in points:
        ha = 'left' if x >= 0 else 'right'
        xt = x + dx if x >= 0 else x - dx
        ax.text(xt, y, 'F', fontsize=7, color='darkgreen', va='center', ha=ha)


def _visible_hyperbolic_foci_vertical(a: float,
                                      b: float,
                                      R: float,
                                      offsets: Sequence[float],
                                      global_axis: np.ndarray,
                                      lo_components: List[np.ndarray],
                                      hi_components: List[np.ndarray],
                                      m: int,
                                      h: float,
                                      npts: int = 900) -> List[Tuple[float, float]]:
    """Foci of real visible hyperbolic generatrices of the vertical type."""
    base = vertical_base_grid(a, b, R, npts=npts)
    shifts = stack_shifts(float(base['L']), h, m)
    c = c_focus(a, b)
    descs = _unique_hyperbolic_descriptors(offsets)
    points: List[Tuple[float, float]] = []

    for shift in shifts:
        for axis_base, d_base in ((base['axis_left'], base['d_left']), (base['axis_right'], base['d_right'])):
            axis_abs = axis_base + float(shift)
            for _, C, sigma in descs:
                C = float(C)
                sigma = float(sigma)
                r_curve = C + sigma * d_base
                mask = np.isfinite(r_curve) & (r_curve > MERGE_TOL)
                if np.count_nonzero(mask) < 8:
                    continue
                if not _curve_matches_visible_boundary(global_axis, axis_abs[mask], r_curve[mask], lo_components, hi_components):
                    continue
                rF = C + sigma * R
                for sF in (-c, +c):
                    points.append((rF, float(shift) + sF))
                    if abs(rF) > 1e-10:
                        points.append((-rF, float(shift) + sF))

    unique = sorted({(round(float(x), 8), round(float(y), 8)) for x, y in points}, key=lambda p: (p[1], p[0]))
    return [(float(x), float(y)) for x, y in unique]


def _visible_hyperbolic_foci_horizontal(a: float,
                                        b: float,
                                        R: float,
                                        offsets: Sequence[float],
                                        global_axis: np.ndarray,
                                        lo_components: List[np.ndarray],
                                        hi_components: List[np.ndarray],
                                        shifts: np.ndarray,
                                        npts: int = 900) -> List[Tuple[float, float]]:
    """Foci of real visible hyperbolic generatrices of the horizontal type."""
    base = horizontal_base_grid(a, b, R, npts=npts)
    c = c_focus(a, b)
    descs = _unique_hyperbolic_descriptors(offsets)
    points: List[Tuple[float, float]] = []

    for shift in np.asarray(shifts, dtype=float):
        axis_abs = base['axis'] + float(shift)
        for _, C, sigma in descs:
            C = float(C)
            sigma = float(sigma)
            r_curve = C + sigma * base['d']
            mask = np.isfinite(r_curve) & (r_curve > MERGE_TOL)
            if np.count_nonzero(mask) < 8:
                continue
            if not _curve_matches_visible_boundary(global_axis, axis_abs[mask], r_curve[mask], lo_components, hi_components):
                continue
            # A horizontal hyperbola has foci along the radial coordinate C±c.
            for axisF in (float(shift) - R, float(shift) + R):
                for rF in (C - c, C + c):
                    points.append((rF, axisF))
                # Mirror side of the full 2D section.
                for rF in (-C - c, -C + c):
                    points.append((rF, axisF))

    unique = sorted({(round(float(x), 8), round(float(y), 8)) for x, y in points}, key=lambda p: (p[1], p[0]))
    return [(float(x), float(y)) for x, y in unique]


def _annotate_vertical_order2_2d(ax, R: float, rmax: float, axis_min: float, axis_max: float):
    ax.axvline(0.0, color='black', linestyle='--', linewidth=1.05, zorder=1)
    ax.text(0.03, 0.97, 'Rotation axis', transform=ax.transAxes,
            ha='left', va='top', fontsize=9)
    ax.axvline(R, color='black', linestyle='--', linewidth=1.0, zorder=1)
    ax.text(R + 0.035 * max(rmax, 1.0), axis_min + 0.72 * (axis_max-axis_min),
            'Focal axis\nof the generatrix', ha='left', va='center', fontsize=9)
    yR = axis_min + 0.16 * (axis_max - axis_min)
    ax.annotate('', xy=(R, yR), xytext=(0.0, yR),
                arrowprops=dict(arrowstyle='<->', lw=1.15, color='black'))
    ax.text(0.5 * R, yR + 0.03 * (axis_max-axis_min), 'R',
            ha='center', va='bottom', fontsize=12)


def _annotate_horizontal_order2_2d(ax, R: float, rmax: float, axis_min: float, axis_max: float):
    ax.axvline(0.0, color='black', linestyle='--', linewidth=1.0, zorder=1)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=1.05, zorder=1)
    ax.text(0.98, 0.52, 'Central horizontal\nsymmetry axis y = 0',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
    ax.axhline(R, color='black', linestyle='--', linewidth=1.0, zorder=1)
    ax.text(0.98, 0.88, 'Focal axis\nof the generatrix y = R',
            transform=ax.transAxes, ha='right', va='top', fontsize=9)
    xR = -0.82 * max(rmax, 1.0)
    ax.annotate('', xy=(xR, R), xytext=(xR, 0.0),
                arrowprops=dict(arrowstyle='<->', lw=1.15, color='black'))
    ax.text(xR - 0.03 * max(rmax, 1.0), 0.5 * R, 'R',
            ha='right', va='center', fontsize=12)


def plot_base_generatrices(a: float,
                           b: float,
                           R: float,
                           outpath: Path | None = None,
                           show: bool = True,
                           npts: int = 1200):
    """Base hyperbolic generatrices of the 2nd order with reference lines and R."""
    c = c_focus(a, b)
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Vertical type.
    ax = axes[0]
    base_v = vertical_base_grid(a, b, R, npts=npts)
    ax.plot(base_v['d_left'], base_v['axis_left'], color='navy', lw=2.0, label='Generatrix')
    ax.plot(base_v['d_right'], base_v['axis_right'], color='navy', lw=2.0)
    ax.axvline(0.0, color='black', linestyle='--', lw=1.05)
    ax.axvline(R, color='black', linestyle='--', lw=1.0)
    ax.text(0.03, 0.97, 'Rotation axis', transform=ax.transAxes, ha='left', va='top', fontsize=9)
    ax.text(R + 0.03 * max(R, 1.0), 0.55 * base_v['L'], 'Focal axis\nof the generatrix', ha='left', va='center', fontsize=9)
    yR = -0.88 * base_v['L']
    ax.annotate('', xy=(R, yR), xytext=(0.0, yR), arrowprops=dict(arrowstyle='<->', lw=1.15, color='black'))
    ax.text(0.5 * R, yR + 0.035 * base_v['L'], 'R', ha='center', va='bottom', fontsize=12)
    for sf in (-c, +c):
        ax.plot(R, sf, 'o', color='darkgreen', ms=5)
        ax.text(R + 0.03 * max(R, 1.0), sf, 'F', color='darkgreen', fontsize=8, va='center')
    ax.set_title(f'Base generatrix: vertical type, a={a:g}, b={b:g}, R={R:g}')
    ax.set_xlabel('Radial coordinate')
    ax.set_ylabel('Coordinate along the common axis')
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='best', fontsize=8)

    # Horizontal type.
    ax = axes[1]
    base_h = horizontal_base_grid(a, b, R, npts=npts)
    ax.plot(base_h['d'], base_h['axis'], color='darkred', lw=2.0, label='Generatrix')
    ax.axvline(0.0, color='black', linestyle='--', lw=1.0)
    ax.axhline(0.0, color='black', linestyle='--', lw=1.05)
    ax.axhline(R, color='black', linestyle='--', lw=1.0)
    ax.text(0.98, 0.52, 'Central horizontal\nsymmetry axis y = 0', transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
    ax.text(0.98, 0.88, 'Focal axis\nof the generatrix y = R', transform=ax.transAxes, ha='right', va='top', fontsize=9)
    xR = -0.20 * max(base_h['d'])
    ax.annotate('', xy=(xR, R), xytext=(xR, 0.0), arrowprops=dict(arrowstyle='<->', lw=1.15, color='black'))
    ax.text(xR - 0.035 * max(base_h['d']), 0.5 * R, 'R', ha='right', va='center', fontsize=12)
    for yf in (-R, +R):
        for xf in (-c, +c):
            ax.plot(xf, yf, 'o', color='darkgreen', ms=5)
            ax.text(xf + 0.02 * max(base_h['d']), yf, 'F', color='darkgreen', fontsize=8, va='center')
    ax.set_title(f'Base generatrix: horizontal type, a={a:g}, b={b:g}, R={R:g}')
    ax.set_xlabel('Radial coordinate')
    ax.set_ylabel('Coordinate along the common axis')
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='best', fontsize=8)

    fig.tight_layout()
    _save_or_show(fig, outpath, show)


# ------------------------------------------------------------
# GENERATING COMPONENTS FOR ONE ORDER
# ------------------------------------------------------------
def vertical_order_union_components(a: float,
                                    b: float,
                                    R: float,
                                    offsets: Sequence[float],
                                    npts: int,
                                    m: int,
                                    h: float) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray], float]:
    base = vertical_base_grid(a, b, R, npts=npts)
    shifts = stack_shifts(float(base["L"]), h, m)

    interval_left = build_recursive_interval_arrays(base["d_left"], offsets)
    interval_right = build_recursive_interval_arrays(base["d_right"], offsets)

    pieces: List[Piece] = []
    for shift in shifts:
        for lo, hi in interval_left:
            pieces.append(Piece(base["axis_left"] + shift, lo, hi))
        for lo, hi in interval_right:
            pieces.append(Piece(base["axis_right"] + shift, lo, hi))

    global_axis, lo_components, hi_components = _build_merged_components_from_pieces(pieces, n_global=max(2500, 3 * npts))
    return global_axis, lo_components, hi_components, float(base["L"])


def horizontal_order_union_components(a: float,
                                      b: float,
                                      R: float,
                                      offsets: Sequence[float],
                                      npts: int,
                                      m: int,
                                      h: float) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray], float, np.ndarray]:
    base = horizontal_base_grid(a, b, R, npts=npts)
    shifts = stack_shifts(float(base["L"]), h, m)

    interval_base = build_recursive_interval_arrays(base["d"], offsets)

    pieces: List[Piece] = []
    for shift in shifts:
        for lo, hi in interval_base:
            pieces.append(Piece(base["axis"] + shift, lo, hi))

    global_axis, lo_components, hi_components = _build_merged_components_from_pieces(pieces, n_global=max(2500, 3 * npts))
    return global_axis, lo_components, hi_components, float(base["L"]), shifts


# ------------------------------------------------------------
# SURFACES OF REVOLUTION
# ------------------------------------------------------------
def revolve_about_x(axis_coord: np.ndarray, radius: np.ndarray, nphi: int = NPHI):
    phi = np.linspace(0.0, 2.0 * np.pi, nphi)
    S, Phi = np.meshgrid(axis_coord, phi, indexing="ij")
    RR = np.tile(radius[:, None], (1, phi.size))
    X = S
    Y = RR * np.cos(Phi)
    Z = RR * np.sin(Phi)
    return X, Y, Z


def revolve_about_y(axis_coord: np.ndarray, radius: np.ndarray, nphi: int = NPHI):
    phi = np.linspace(0.0, 2.0 * np.pi, nphi)
    U, Phi = np.meshgrid(axis_coord, phi, indexing="ij")
    RR = np.tile(radius[:, None], (1, phi.size))
    X = RR * np.cos(Phi)
    Y = U
    Z = RR * np.sin(Phi)
    return X, Y, Z


@dataclass
class PlotBounds3D:
    xmin: float = math.inf
    xmax: float = -math.inf
    ymin: float = math.inf
    ymax: float = -math.inf
    zmin: float = math.inf
    zmax: float = -math.inf

    def update(self, X: np.ndarray, Y: np.ndarray, Z: np.ndarray):
        self.xmin = min(self.xmin, float(np.nanmin(X)))
        self.xmax = max(self.xmax, float(np.nanmax(X)))
        self.ymin = min(self.ymin, float(np.nanmin(Y)))
        self.ymax = max(self.ymax, float(np.nanmax(Y)))
        self.zmin = min(self.zmin, float(np.nanmin(Z)))
        self.zmax = max(self.zmax, float(np.nanmax(Z)))


def set_axes_equal_real_3d(ax, bounds: PlotBounds3D):
    xmid = 0.5 * (bounds.xmin + bounds.xmax)
    ymid = 0.5 * (bounds.ymin + bounds.ymax)
    zmid = 0.5 * (bounds.zmin + bounds.zmax)
    half = 0.5 * max(bounds.xmax - bounds.xmin, bounds.ymax - bounds.ymin, bounds.zmax - bounds.zmin)
    if half <= 0:
        half = 1.0
    ax.set_xlim(xmid - half, xmid + half)
    ax.set_ylim(ymid - half, ymid + half)
    ax.set_zlim(zmid - half, zmid + half)
    try:
        ax.set_box_aspect((1.0, 1.0, 1.0))
    except Exception:
        pass
    try:
        ax.set_proj_type("ortho")
    except Exception:
        pass


def _plot_surface_piece_x(ax, axis_seg: np.ndarray, radius_seg: np.ndarray, bounds: PlotBounds3D, stride: int = 3):
    X, Y, Z = revolve_about_x(axis_seg, radius_seg)
    ax.plot_surface(
        X[::stride, ::stride], Y[::stride, ::stride], Z[::stride, ::stride],
        linewidth=0, edgecolor="none", alpha=SURFACE_ALPHA, antialiased=False, shade=True
    )
    bounds.update(X, Y, Z)


def _plot_surface_piece_y(ax, axis_seg: np.ndarray, radius_seg: np.ndarray, bounds: PlotBounds3D,
                          split_positions: Sequence[float] = (), stride: int = 3):
    # For the horizontal type, split at row-shift positions
    # so bends in central zones are not smoothed.
    split_ids = set()
    for pos in split_positions:
        idx = int(np.argmin(np.abs(axis_seg - pos)))
        if 1 <= idx < axis_seg.size - 1 and np.isclose(axis_seg[idx], pos, atol=1e-10):
            split_ids.add(idx)
    starts = [0] + sorted(split_ids)
    ends = sorted(split_ids) + [axis_seg.size - 1]
    for st, en in zip(starts, ends):
        seg_axis = axis_seg[st:en + 1]
        seg_rad = radius_seg[st:en + 1]
        if seg_axis.size < 2:
            continue
        X, Y, Z = revolve_about_y(seg_axis, seg_rad)
        ax.plot_surface(
            X[::stride, ::stride], Y[::stride, ::stride], Z[::stride, ::stride],
            linewidth=0, edgecolor="none", alpha=SURFACE_ALPHA, antialiased=False, shade=True
        )
        bounds.update(X, Y, Z)


# ------------------------------------------------------------
# 2D DRAWINGS OF THE UNION VOLUME
# ------------------------------------------------------------
def _plot_component_2d(ax, axis: np.ndarray, lo: np.ndarray, hi: np.ndarray):
    """
    Plots only real boundaries of the union volume.

    IMPORTANT:
    - the script does not add anything;
    - fill_betweenx / fill_between are not used;
    - if component numbering jumps after interval merging,
      the curve is split and is not connected by an artificial straight segment.
    """
    mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
    for st, en in _contiguous_segments(mask, axis):
        yy = axis[st:en]
        lo_seg = lo[st:en]
        hi_seg = hi[st:en]

        _plot_boundary_curve(ax, hi_seg, yy, lw=1.8)
        _plot_boundary_curve(ax, -hi_seg, yy, lw=1.8)

        if np.nanmax(lo_seg) > MERGE_TOL:
            _plot_boundary_curve(ax, lo_seg, yy, lw=1.5)
            _plot_boundary_curve(ax, -lo_seg, yy, lw=1.5)


def plot_union_section_vertical(order: int,
                                global_axis: np.ndarray,
                                lo_components: List[np.ndarray],
                                hi_components: List[np.ndarray],
                                a: float,
                                b: float,
                                R: float,
                                offsets: Sequence[float],
                                m: int,
                                h: float,
                                outpath: Path | None = None,
                                show: bool = True):
    fig, ax = plt.subplots(figsize=(9, 10))
    for lo, hi in zip(lo_components, hi_components):
        _plot_component_2d(ax, global_axis, lo, hi)

    rmax = 0.0
    for hi in hi_components:
        if np.any(np.isfinite(hi)):
            rmax = max(rmax, float(np.nanmax(hi)))
    if rmax <= 0:
        rmax = 1.0

    axis_min = float(np.nanmin(global_axis))
    axis_max = float(np.nanmax(global_axis))
    # Show the rotation axis as a dashed line on all 2D plots.
    ax.axvline(0.0, color='black', linestyle='--', linewidth=1.0, zorder=1)
    focus_points = _visible_hyperbolic_foci_vertical(
        a=a, b=b, R=R, offsets=offsets, global_axis=global_axis,
        lo_components=lo_components, hi_components=hi_components,
        m=m, h=h, npts=max(900, int(global_axis.size // max(1, m)))
    )
    _plot_focus_points_2d(ax, focus_points, scale_ref=max(rmax, axis_max-axis_min))
    if order == 2:
        _annotate_vertical_order2_2d(ax, R=R, rmax=rmax, axis_min=axis_min, axis_max=axis_max)

    ax.set_xlim(-1.10 * max(rmax, R), 1.18 * max(rmax, R) if order == 2 else 1.10 * max(rmax, R))
    ax.set_title(f"Vertical type, order {order}: union volume, 2D section, m={m}, h={h}")
    ax.set_xlabel("Radial coordinate")
    ax.set_ylabel("Coordinate along the common axis")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


def plot_union_section_horizontal(order: int,
                                  global_axis: np.ndarray,
                                  lo_components: List[np.ndarray],
                                  hi_components: List[np.ndarray],
                                  a: float,
                                  b: float,
                                  R: float,
                                  offsets: Sequence[float],
                                  shifts: np.ndarray,
                                  m: int,
                                  h: float,
                                  outpath: Path | None = None,
                                  show: bool = True):
    fig, ax = plt.subplots(figsize=(10, 8))
    for lo, hi in zip(lo_components, hi_components):
        _plot_component_2d(ax, global_axis, lo, hi)

    rmax = 0.0
    for hi in hi_components:
        if np.any(np.isfinite(hi)):
            rmax = max(rmax, float(np.nanmax(hi)))
    if rmax <= 0:
        rmax = 1.0

    axis_min = float(np.nanmin(global_axis))
    axis_max = float(np.nanmax(global_axis))
    # Show the rotation axis as a dashed line on all 2D plots.
    ax.axhline(0.0, color='black', linestyle='--', linewidth=1.0, zorder=1)
    focus_points = _visible_hyperbolic_foci_horizontal(
        a=a, b=b, R=R, offsets=offsets, global_axis=global_axis,
        lo_components=lo_components, hi_components=hi_components,
        shifts=shifts, npts=max(900, int(global_axis.size // max(1, m)))
    )
    _plot_focus_points_2d(ax, focus_points, scale_ref=max(rmax, axis_max-axis_min))
    if order == 2:
        _annotate_horizontal_order2_2d(ax, R=R, rmax=rmax, axis_min=axis_min, axis_max=axis_max)

    ax.set_xlim(-1.10 * max(rmax, R), 1.10 * max(rmax, R))
    ax.set_title(f"Horizontal type, order {order}: union volume, 2D section, m={m}, h={h}")
    ax.set_xlabel("Radial coordinate")
    ax.set_ylabel("Coordinate along the common axis")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


# ------------------------------------------------------------
# BRIGHT 2D SECTION ON 3D FOR CHECKING THE FULL INTERNAL VOLUME
# ------------------------------------------------------------
def _plot_bright_2d_section_on_3d_vertical(ax, global_axis: np.ndarray, lo_components: List[np.ndarray], hi_components: List[np.ndarray]):
    """
    Thin black meridional 2D section in the Y=0 plane for the vertical type.
    IMPORTANT: the profile is split at jumps so artificial straight connections do not appear.
    """
    y_plane = 0.0
    y_eps = 1.0e-6
    color = 'black'

    def draw_split(axis_arr: np.ndarray, rad_arr: np.ndarray, lw_main: float, lw_aux: float):
        for a0, b0 in _split_on_artificial_jumps(rad_arr):
            xx = axis_arr[a0:b0]
            rr = rad_arr[a0:b0]
            if xx.size < 2:
                continue
            ax.plot(xx, np.full_like(xx, y_plane), rr, color=color, linewidth=lw_main)
            ax.plot(xx, np.full_like(xx, y_eps), rr, color=color, linewidth=lw_aux, alpha=0.9)

    for lo, hi in zip(lo_components, hi_components):
        mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
        for st, en in _contiguous_segments(mask, global_axis):
            axis_seg = global_axis[st:en]
            lo_seg = lo[st:en]
            hi_seg = hi[st:en]
            # outer boundary
            draw_split(axis_seg,  hi_seg, 1.4, 1.0)
            draw_split(axis_seg, -hi_seg, 1.4, 1.0)
            # inner boundary
            if np.nanmax(lo_seg) > MERGE_TOL:
                draw_split(axis_seg,  lo_seg, 1.2, 0.9)
                draw_split(axis_seg, -lo_seg, 1.2, 0.9)


def _plot_bright_2d_section_on_3d_horizontal(ax, global_axis: np.ndarray, lo_components: List[np.ndarray], hi_components: List[np.ndarray]):
    """
    Thin black meridional 2D section in the Z=0 plane for the horizontal type.
    IMPORTANT: the profile is split at jumps so artificial straight connections do not appear.
    """
    z_plane = 0.0
    z_eps = 1.0e-6
    color = 'black'

    def draw_split(rad_arr: np.ndarray, axis_arr: np.ndarray, lw_main: float, lw_aux: float):
        for a0, b0 in _split_on_artificial_jumps(rad_arr):
            rr = rad_arr[a0:b0]
            yy = axis_arr[a0:b0]
            if rr.size < 2:
                continue
            ax.plot(rr, yy, np.full_like(yy, z_plane), color=color, linewidth=lw_main)
            ax.plot(rr, yy, np.full_like(yy, z_eps), color=color, linewidth=lw_aux, alpha=0.9)

    for lo, hi in zip(lo_components, hi_components):
        mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
        for st, en in _contiguous_segments(mask, global_axis):
            axis_seg = global_axis[st:en]
            lo_seg = lo[st:en]
            hi_seg = hi[st:en]
            # outer boundary
            draw_split( hi_seg, axis_seg, 1.4, 1.0)
            draw_split(-hi_seg, axis_seg, 1.4, 1.0)
            # inner boundary
            if np.nanmax(lo_seg) > MERGE_TOL:
                draw_split( lo_seg, axis_seg, 1.2, 0.9)
                draw_split(-lo_seg, axis_seg, 1.2, 0.9)


# ------------------------------------------------------------
# 3D SURFACES OF THE UNION VOLUME
# ------------------------------------------------------------
def plot_union_surface_vertical(order: int,
                                global_axis: np.ndarray,
                                lo_components: List[np.ndarray],
                                hi_components: List[np.ndarray],
                                m: int,
                                h: float,
                                outpath: Path | None = None,
                                show: bool = True,
                                stride: int = 3):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    bounds = PlotBounds3D()

    for lo, hi in zip(lo_components, hi_components):
        mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
        for st, en in _contiguous_segments(mask, global_axis):
            axis_seg = global_axis[st:en]
            lo_seg = lo[st:en]
            hi_seg = hi[st:en]
            _plot_surface_piece_x(ax, axis_seg, hi_seg, bounds, stride=stride)
            if np.nanmax(lo_seg) > MERGE_TOL:
                _plot_surface_piece_x(ax, axis_seg, lo_seg, bounds, stride=stride)

    _plot_bright_2d_section_on_3d_vertical(ax, global_axis, lo_components, hi_components)
    set_axes_equal_real_3d(ax, bounds)
    ax.grid(True)
    ax.set_title(f"Vertical type, order {order}: 3D surface of the union volume + bright 2D section, m={m}, h={h}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


def plot_union_surface_horizontal(order: int,
                                  global_axis: np.ndarray,
                                  lo_components: List[np.ndarray],
                                  hi_components: List[np.ndarray],
                                  shifts: np.ndarray,
                                  m: int,
                                  h: float,
                                  outpath: Path | None = None,
                                  show: bool = True,
                                  stride: int = 3):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    bounds = PlotBounds3D()

    for lo, hi in zip(lo_components, hi_components):
        mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
        for st, en in _contiguous_segments(mask, global_axis):
            axis_seg = global_axis[st:en]
            lo_seg = lo[st:en]
            hi_seg = hi[st:en]
            _plot_surface_piece_y(ax, axis_seg, hi_seg, bounds, split_positions=shifts, stride=stride)
            if np.nanmax(lo_seg) > MERGE_TOL:
                _plot_surface_piece_y(ax, axis_seg, lo_seg, bounds, split_positions=shifts, stride=stride)

    _plot_bright_2d_section_on_3d_horizontal(ax, global_axis, lo_components, hi_components)
    set_axes_equal_real_3d(ax, bounds)
    ax.grid(True)
    ax.set_title(f"Horizontal type, order {order}: 3D surface of the union volume + bright 2D section, m={m}, h={h}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


# ------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------
def run_pseudohyperboloids(a: float = 0.5,
                           b: float = 1.0,
                           R: float = 15.0,
                           offsets: Sequence[float] = (14.0, 30.0),
                           geometry_type: str = "both",   # vertical / horizontal / both
                           mode: str = "all",             # section / surface / all
                           outdir: str = "pseudohyperboloids_notebook_output",
                           npts: int = 900,
                           h: float = -2.0,
                           m: int = 3,
                           show: bool = True,
                           make_base_generatrices: bool = True) -> Dict[str, object]:
    """
    Main function.

    The FINAL UNION VOLUME is built for each order from 2 to n,
    including component intersections and nesting.
    """
    offsets = list(offsets)
    outdir_path = Path(outdir)
    _ensure_dir(outdir_path)

    if geometry_type not in {"vertical", "horizontal", "both"}:
        raise ValueError("geometry_type must be 'vertical', 'horizontal' or 'both'")
    if mode not in {"section", "surface", "all"}:
        raise ValueError("mode must be 'section', 'surface' or 'all'")
    if a <= 0 or b <= 0 or R <= 0:
        raise ValueError("Parameters a, b, R must be positive")
    if int(m) != m or m < 1:
        raise ValueError("Parameter m must be an integer >= 1")

    result = {
        "a": a,
        "b": b,
        "R": R,
        "offsets": offsets,
        "n_order": len(offsets) + 2,
        "focus_c": c_focus(a, b),
        "morphology": morphology_class(R, offsets),
        "m": int(m),
        "h": float(h),
        "volume_mode": "union_of_all_toroidal_components",
        "saved_files": [],
        "warnings": [],
    }

    print(f"a = {a}")
    print(f"b = {b}")
    print(f"R = {R}")
    print(f"offsets = {offsets}")
    print(f"Order n = {len(offsets) + 2}")
    print(f"Focal distance c = {c_focus(a, b):.6f}")
    print(f"Morphological class = {result['morphology']}")
    print(f"m = {m}")
    print(f"h = {h}")
    print("Construction mode = union volume as the union of all components")
    print("Rule: only the common overlapping part of the volume is removed,")
    print("while the pseudohyperboloidal toroids themselves are not deleted as generating objects.")
    print(f"Output folder = {outdir_path.resolve()}")

    if make_base_generatrices:
        path = outdir_path / f"base_generatrices_a{a}_b{b}_R{R}.png"
        plot_base_generatrices(a=a, b=b, R=R, outpath=path, show=show, npts=npts)
        result["saved_files"].append(str(path))

    if geometry_type in {"vertical", "both"}:
        base_v = vertical_base_grid(a, b, R, npts=npts)
        result["warnings"].extend([f"vertical-left: {w}" for w in validate_offsets(base_v["d_left"], offsets)])
        result["warnings"].extend([f"vertical-right: {w}" for w in validate_offsets(base_v["d_right"], offsets)])

        for p in range(2, len(offsets) + 3):
            used_offsets = offsets[:max(0, p - 2)]
            global_axis, lo_components, hi_components, _ = vertical_order_union_components(
                a, b, R, used_offsets, npts=npts, m=m, h=h
            )
            if mode in {"section", "all"}:
                path = outdir_path / f"vertical_order_{p:02d}_section_m{m}_h{h}.png"
                plot_union_section_vertical(
                    order=p,
                    global_axis=global_axis,
                    lo_components=lo_components,
                    hi_components=hi_components,
                    a=a,
                    b=b,
                    R=R,
                    offsets=used_offsets,
                    m=m,
                    h=h,
                    outpath=path,
                    show=show,
                )
                result["saved_files"].append(str(path))
            if mode in {"surface", "all"}:
                path = outdir_path / f"vertical_order_{p:02d}_surface_m{m}_h{h}.png"
                plot_union_surface_vertical(
                    order=p,
                    global_axis=global_axis,
                    lo_components=lo_components,
                    hi_components=hi_components,
                    m=m,
                    h=h,
                    outpath=path,
                    show=show,
                )
                result["saved_files"].append(str(path))

    if geometry_type in {"horizontal", "both"}:
        base_h = horizontal_base_grid(a, b, R, npts=npts)
        result["warnings"].extend([f"horizontal: {w}" for w in validate_offsets(base_h["d"], offsets)])

        for p in range(2, len(offsets) + 3):
            used_offsets = offsets[:max(0, p - 2)]
            global_axis, lo_components, hi_components, _, shifts = horizontal_order_union_components(
                a, b, R, used_offsets, npts=npts, m=m, h=h
            )
            if mode in {"section", "all"}:
                path = outdir_path / f"horizontal_order_{p:02d}_section_m{m}_h{h}.png"
                plot_union_section_horizontal(
                    order=p,
                    global_axis=global_axis,
                    lo_components=lo_components,
                    hi_components=hi_components,
                    a=a,
                    b=b,
                    R=R,
                    offsets=used_offsets,
                    shifts=shifts,
                    m=m,
                    h=h,
                    outpath=path,
                    show=show,
                )
                result["saved_files"].append(str(path))
            if mode in {"surface", "all"}:
                path = outdir_path / f"horizontal_order_{p:02d}_surface_m{m}_h{h}.png"
                plot_union_surface_horizontal(
                    order=p,
                    global_axis=global_axis,
                    lo_components=lo_components,
                    hi_components=hi_components,
                    shifts=shifts,
                    m=m,
                    h=h,
                    outpath=path,
                    show=show,
                )
                result["saved_files"].append(str(path))

    if result["warnings"]:
        print("\nWARNINGS:")
        for w in result["warnings"]:
            print(" -", w)

    print("\nSAVED FILES:")
    for f in result["saved_files"]:
        print(" -", f)

    return result


# ------------------------------------------------------------
# DEFAULT PARAMETERS FOR THIS TOPIC
# ------------------------------------------------------------
if __name__ == "__main__":
    a = 0.5
    b = 1.0
    R = 15.0
    offsets = [14.0, 30.0]    # R1, R2, ..., R_{n-2}
    geometry_type = "both"   # "vertical", "horizontal", "both"
    mode = "all"             # "section", "surface", "all"
    outdir = "pseudohyperboloids_notebook_output"
    npts = 900
    m = 3
    h = -2.0
    show = True

    results = run_pseudohyperboloids(
        a=a,
        b=b,
        R=R,
        offsets=offsets,
        geometry_type=geometry_type,
        mode=mode,
        outdir=outdir,
        npts=npts,
        m=m,
        h=h,
        show=show,
    )

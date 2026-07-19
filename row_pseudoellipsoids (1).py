# -*- coding: utf-8 -*-
"""
Row n-th order pseudoellipsoids.

NEW CONCEPT
---------------
1. Only two types are considered:
   - horizontal: K < 1
   - vertical:   K > 1
   The value K = 1 is not considered.

2. For both types, the base condition is:
   - h1 = 0
   - h  = 0

   That is, the base generatrix is built from two quarter-ellipse segments
   without an equatorial window and without end windows.

3. The number of constructions is aligned with pseudohyperboloids and
   pseudoparaboloids:
   - 2nd order
   - 3rd order
   - 4th order
   - rows (m copies, step is defined by h_row)

4. On the generatrices and on the 2D view for the 2nd order, reference lines
   and parameter R are shown by analogy with pseudoparaboloids.

5. In 3D, a bright meridional 2D section along the rotation axis is shown
   to visually verify that the internal volume has not been lost.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 220,
    "font.size": 10,
    "axes.grid": True,
    "grid.alpha": 0.24,
    "font.family": "DejaVu Sans",
})

SURFACE_ALPHA = 0.28
MERGE_TOL = 1.0e-9
NPHI = 90


# ============================================================
# BASE GEOMETRY OF TWO QUARTER ELLIPSES
# ============================================================

def y_left(x: np.ndarray, a: float, b: float) -> np.ndarray:
    val = 1.0 - ((x + a) ** 2) / (a ** 2)
    return b * np.sqrt(np.clip(val, 0.0, None))


def y_right(x: np.ndarray, a: float, b: float, h1: float) -> np.ndarray:
    val = 1.0 - ((x - a - h1) ** 2) / (a ** 2)
    return b * np.sqrt(np.clip(val, 0.0, None))


def build_signed_profile(a: float, b: float, h1: float, n: int = 1800):
    """
    Base upper generatrix from two quarter ellipses.
    For h1=0 they meet at x=0.
    """
    xmin = min(-a, h1)
    xmax = max(0.0, a + h1)
    x = np.linspace(xmin, xmax, n)

    yL = np.full_like(x, np.nan, dtype=float)
    yR = np.full_like(x, np.nan, dtype=float)

    maskL = (x >= -a) & (x <= 0.0)
    maskR = (x >= h1) & (x <= a + h1)

    yL[maskL] = y_left(x[maskL], a, b)
    yR[maskR] = y_right(x[maskR], a, b, h1)

    y = np.full_like(x, np.nan, dtype=float)

    if h1 >= 0:
        onlyL = maskL & (~maskR)
        onlyR = maskR & (~maskL)
        y[onlyL] = yL[onlyL]
        y[onlyR] = yR[onlyR]

        gap = (x > 0.0) & (x < h1)
        y[gap] = 0.0
        y[np.isclose(x, 0.0)] = 0.0
        if not np.isclose(h1, 0.0):
            y[np.isclose(x, h1)] = 0.0
    else:
        both = maskL & maskR
        onlyL = maskL & (~maskR)
        onlyR = maskR & (~maskL)
        y[onlyL] = yL[onlyL]
        y[onlyR] = yR[onlyR]
        y[both] = np.maximum(yL[both], yR[both])

    return x, y, yL, yR, maskL, maskR


def compute_foci(a: float, b: float, h1: float, R: float, eps: float = 1e-12):
    if b > a + eps:
        c = np.sqrt(b ** 2 - a ** 2)
        foci_1d = [(-a, -c, "F1"), (-a, +c, "F2"),
                   (a + h1, -c, "F3"), (a + h1, +c, "F4")]
        foci_2d = foci_1d + [
            (-a, 2 * R + c, "F5"), (-a, 2 * R - c, "F6"),
            (a + h1, 2 * R + c, "F7"), (a + h1, 2 * R - c, "F8")]
    elif a > b + eps:
        c = np.sqrt(a ** 2 - b ** 2)
        foci_1d = [(-a - c, 0.0, "F1"), (-a + c, 0.0, "F2"),
                   (a + h1 - c, 0.0, "F3"), (a + h1 + c, 0.0, "F4")]
        foci_2d = foci_1d + [
            (-a - c, 2 * R, "F5"), (-a + c, 2 * R, "F6"),
            (a + h1 - c, 2 * R, "F7"), (a + h1 + c, 2 * R, "F8")]
    else:
        c = 0.0
        foci_1d = [(-a, 0.0, "F1=F2"), (a + h1, 0.0, "F3=F4")]
        foci_2d = foci_1d + [(-a, 2 * R, "F5=F6"), (a + h1, 2 * R, "F7=F8")]
    return c, foci_1d, foci_2d


def _base_focus_pairs_labeled(a: float, K: float, h1: float = 0.0) -> List[Tuple[float, float, str]]:
    """
    Real foci of the source ellipses with unified numbering:
    - for the horizontal ellipse: F1 on the left, F2 on the right;
    - for the vertical ellipse: F1 at the bottom, F2 at the top.

    The same numbering direction is used for each ellipse.
    """
    b, h1, _R = params_from_K(a, K, h1=h1, h=0.0)
    eps = 1.0e-12
    if b > a + eps:
        c = float(np.sqrt(b ** 2 - a ** 2))
        return [(-a, -c, 'F1'), (-a, +c, 'F2'),
                (a + h1, -c, 'F1'), (a + h1, +c, 'F2')]
    if a > b + eps:
        c = float(np.sqrt(a ** 2 - b ** 2))
        return [(-a - c, 0.0, 'F1'), (-a + c, 0.0, 'F2'),
                (a + h1 - c, 0.0, 'F1'), (a + h1 + c, 0.0, 'F2')]
    return [(-a, 0.0, 'F1=F2'), (a + h1, 0.0, 'F1=F2')]

def classify_K(K: float, eps: float = 1e-12) -> str:
    if math.isclose(K, 1.0, rel_tol=eps, abs_tol=eps):
        raise ValueError("K = 1 is not considered")
    return "Horizontal" if K < 1.0 else "Vertical"


def params_from_K(a: float, K: float, h1: float = 0.0, h: float = 0.0) -> Tuple[float, float, float]:
    """
    As in the source pseudoellipsoid script:
    - b = a*K (when a=1 this is b=K);
    - h1 enters the geometry of the right ellipse: the right ellipse center is x = a + h1;
    - h enters only the level of the rotation axis: R = b + h
      (when a=1 this is R = K + h).
    """
    typ = classify_K(K)
    b = float(a) * float(K)
    h1 = float(h1)
    h = float(h)
    R = b + h
    return b, h1, R


def ellipsoid_base_grid(a: float, K: float, h1: float = 0.0, h: float = 0.0, npts: int = 900) -> Dict[str, np.ndarray]:
    b, h1, R = params_from_K(a, K, h1=h1, h=h)
    x, y, yL, yR, maskL, maskR = build_signed_profile(a, b, h1, n=npts)
    y_act = np.where(np.isnan(y), np.nan, np.minimum(y, R))
    d = np.where(np.isnan(y_act), np.nan, R - y_act)
    d = np.maximum(d, 0.0)

    finite = np.isfinite(d)
    if not np.any(finite):
        raise ValueError("The base profile is empty")

    return {
        "axis": x,
        "d": d,
        "y_profile": y,
        "y_left": yL,
        "y_right": yR,
        "mask_left": maskL,
        "mask_right": maskR,
        "xmin": float(np.nanmin(x[finite])),
        "xmax": float(np.nanmax(x[finite])),
        "width": float(np.nanmax(x[finite]) - np.nanmin(x[finite])),
        "max_distance": float(np.nanmax(d[finite])),
        "a": float(a),
        "b": float(b),
        "h1": float(h1),
        "R": float(R),
        "type": classify_K(K),
    }


# ============================================================
# MERGE AND ORDER RECURSION
# ============================================================

@dataclass
class Piece:
    axis: np.ndarray
    lo: np.ndarray
    hi: np.ndarray


def build_recursive_interval_arrays(base_distance: np.ndarray,
                                    offsets: Sequence[float]) -> List[Tuple[np.ndarray, np.ndarray]]:
    d0 = np.asarray(base_distance, dtype=float)
    zero = np.zeros_like(d0)
    intervals = [(zero.copy(), d0.copy())]

    for Rk in offsets:
        nxt: List[Tuple[np.ndarray, np.ndarray]] = []
        for lo, hi in intervals:
            a0 = np.maximum(Rk - hi, 0.0)
            b0 = np.maximum(Rk - lo, 0.0)
            c0 = np.maximum(Rk + lo, 0.0)
            d1 = np.maximum(Rk + hi, 0.0)
            nxt.append((a0, b0))
            nxt.append((c0, d1))
        intervals = nxt

    return intervals


def stack_shifts_by_width(axis_width: float, h_row: float, m: int) -> np.ndarray:
    if int(m) != m or m < 1:
        raise ValueError("m must be an integer >= 1")
    m = int(m)
    step = float(axis_width) + float(h_row)
    return -np.arange(m, dtype=float) * step


def _interp_piece_to_global_axis(piece: Piece, global_axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = piece.axis
    lo = piece.lo
    hi = piece.hi
    finite = np.isfinite(lo) & np.isfinite(hi)
    if np.count_nonzero(finite) < 2:
        return np.full_like(global_axis, np.nan), np.full_like(global_axis, np.nan)
    x_f = x[finite]
    lo_f = lo[finite]
    hi_f = hi[finite]
    lo_i = np.interp(global_axis, x_f, lo_f)
    hi_i = np.interp(global_axis, x_f, hi_f)
    mask = (global_axis >= x_f.min() - 1e-12) & (global_axis <= x_f.max() + 1e-12)
    lo_i[~mask] = np.nan
    hi_i[~mask] = np.nan
    return lo_i, hi_i


def _merge_scalar_intervals(intervals: List[Tuple[float, float]], tol: float = MERGE_TOL) -> List[Tuple[float, float]]:
    if not intervals:
        return []
    data = [
        (float(min(a0, b0)), float(max(a0, b0)))
        for a0, b0 in intervals
        if np.isfinite(a0) and np.isfinite(b0) and b0 > a0 + tol
    ]
    if not data:
        return []
    data.sort(key=lambda t: (t[0], t[1]))
    merged = [list(data[0])]
    for a0, b0 in data[1:]:
        if a0 <= merged[-1][1] + tol:
            merged[-1][1] = max(merged[-1][1], b0)
        else:
            merged.append([a0, b0])
    return [(a0, b0) for a0, b0 in merged]


def _build_merged_components_from_pieces(pieces: List[Piece],
                                         n_global: int | None = None) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    all_axis = np.concatenate([p.axis[np.isfinite(p.axis)] for p in pieces])
    global_axis = np.unique(np.round(all_axis, 12))
    global_axis.sort()

    if n_global is not None and n_global > global_axis.size:
        dense = np.linspace(global_axis.min(), global_axis.max(), n_global)
        global_axis = np.unique(np.concatenate([global_axis, dense]))
        global_axis.sort()

    interpolated = [_interp_piece_to_global_axis(piece, global_axis) for piece in pieces]
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
        for k, (a0, b0) in enumerate(merged):
            lo_components[k][j] = a0
            hi_components[k][j] = b0

    return global_axis, lo_components, hi_components


def ellipsoid_order_union_components(a: float,
                                     K: float,
                                     h1: float,
                                     h: float,
                                     offsets: Sequence[float],
                                     npts: int,
                                     m: int,
                                     h_row: float):
    base = ellipsoid_base_grid(a=a, K=K, h1=h1, h=h, npts=npts)
    shifts = stack_shifts_by_width(float(base["width"]), h_row, m)
    interval_base = build_recursive_interval_arrays(base["d"], offsets)

    pieces: List[Piece] = []
    for shift in shifts:
        for lo, hi in interval_base:
            pieces.append(Piece(base["axis"] + shift, lo, hi))

    global_axis, lo_components, hi_components = _build_merged_components_from_pieces(
        pieces,
        n_global=max(4000, 6 * npts, 1200 * (len(offsets) + 1)),
    )
    return global_axis, lo_components, hi_components, shifts, base


# ============================================================
# DRAWING HELPERS
# ============================================================

def _ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def _save_or_show(fig, path: Optional[Path], show: bool):
    """
    Notebook mode:
    - if path is specified, the figure is additionally saved to a file;
    - if show=True, the figure is displayed directly in the notebook via plt.show();
    - if show=False, the figure is closed without interactive display.
    """
    if path is not None:
        fig.savefig(path, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def _axis_limits_from_components(axis: np.ndarray, hi_components: List[np.ndarray]) -> Tuple[float, float, float]:
    rmax = 0.0
    for hi in hi_components:
        if np.any(np.isfinite(hi)):
            rmax = max(rmax, float(np.nanmax(hi)))
    if rmax <= 0:
        rmax = 1.0
    return rmax, float(np.nanmin(axis)), float(np.nanmax(axis))


def _contiguous_segments(mask: np.ndarray, axis: Optional[np.ndarray] = None) -> List[Tuple[int, int]]:
    mask = np.asarray(mask, dtype=bool).copy()
    if mask.size == 0:
        return []
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
    radius = np.asarray(radius, dtype=float)
    n = radius.size
    if n < 2:
        return [(0, n)] if n > 0 else []
    dr = np.abs(np.diff(radius))
    valid = dr[np.isfinite(dr)]
    if valid.size == 0:
        return [(0, n)]
    finite_r = radius[np.isfinite(radius)]
    if finite_r.size == 0:
        return [(0, n)]
    r_range = float(np.nanmax(finite_r) - np.nanmin(finite_r))
    med = float(np.median(valid))
    jump_thr = max(1e-8, 12.0 * med, 0.08 * r_range)
    break_ids = np.where(dr > jump_thr)[0] + 1
    if break_ids.size == 0:
        return [(0, n)]
    segments: List[Tuple[int, int]] = []
    st = 0
    for br in break_ids:
        if br - st >= 2:
            segments.append((st, int(br)))
        st = int(br)
    if n - st >= 2:
        segments.append((st, n))
    return segments


def _plot_boundary_curve(ax, axis: np.ndarray, radius: np.ndarray, lw: float, color=None, alpha=1.0):
    for st, en in _split_on_artificial_jumps(radius):
        if en - st >= 2:
            ax.plot(axis[st:en], radius[st:en], lw=lw, color=color, alpha=alpha)


def _plot_component_2d(ax, axis: np.ndarray, lo: np.ndarray, hi: np.ndarray):
    mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
    for st, en in _contiguous_segments(mask, axis):
        xx = axis[st:en]
        lo_seg = lo[st:en]
        hi_seg = hi[st:en]
        _plot_boundary_curve(ax, xx, hi_seg, lw=1.8)
        _plot_boundary_curve(ax, xx, -hi_seg, lw=1.8)
        if np.nanmax(lo_seg) > MERGE_TOL:
            _plot_boundary_curve(ax, xx, lo_seg, lw=1.4)
            _plot_boundary_curve(ax, xx, -lo_seg, lw=1.4)


# ============================================================
# EXACT FOCI OF GENERATING ELLIPSES
# ============================================================
def _desc_add_const(C: float, desc):
    if desc[0] == "const":
        return ("const", float(C + desc[1]))
    _, c0, sigma = desc
    return ("ellip", float(C + c0), float(sigma))


def _desc_sub_const(C: float, desc):
    if desc[0] == "const":
        return ("const", float(C - desc[1]))
    _, c0, sigma = desc
    return ("ellip", float(C - c0), float(-sigma))


def build_recursive_boundary_descriptors(offsets: Sequence[float]):
    """
    Exact descriptors of recursion boundaries.

    A boundary of the form ("ellip", C, sigma) means:
        r(x) = C + sigma * d(x),
    where d(x) = R - y_ellipse(x).

    IMPORTANT:
    - Merge creates no new ellipses and no new foci;
    - foci are taken only from these generating elliptical boundaries;
    - the visible envelope after Merge is not used to compute foci.
    """
    intervals = [(("const", 0.0), ("ellip", 0.0, +1.0))]
    for Rk in offsets:
        nxt = []
        for lo, hi in intervals:
            nxt.append((_desc_sub_const(Rk, hi), _desc_sub_const(Rk, lo)))
            nxt.append((_desc_add_const(Rk, lo), _desc_add_const(Rk, hi)))
        intervals = nxt
    return intervals


def _unique_elliptic_descriptors(offsets: Sequence[float]):
    intervals = build_recursive_boundary_descriptors(offsets)
    seen = set()
    out = []
    for lo, hi in intervals:
        for desc in (lo, hi):
            if desc[0] != "ellip":
                continue
            key = (round(float(desc[1]), 10), round(float(desc[2]), 10))
            if key not in seen:
                seen.add(key)
                out.append(desc)
    return out


def _descriptor_has_active_elliptic_segment(desc, d_base: np.ndarray) -> bool:
    """
    Discards only completely nonexistent boundaries.
    Boundaries partially clipped by the r=0 axis are kept because their
    visible part still has the focus of the source ellipse.
    """
    _, C, sigma = desc
    rad = float(C) + float(sigma) * np.asarray(d_base, dtype=float)
    return bool(np.any(np.isfinite(rad) & (rad > MERGE_TOL)))


def _exact_focus_points_for_elliptic_order(a: float,
                                           K: float,
                                           offsets: Sequence[float],
                                           shifts: np.ndarray) -> List[Tuple[float, float]]:
    """
    Foci of all real elliptical generating boundaries of this order.

    If r = C + d(x) = C + R - y(x), then the focus is:
        r_F = C + R - y_F.

    If r = C - d(x) = C - R + y(x), then the focus is:
        r_F = C - R + y_F.

    For the full 2D section, the mirror side is also shown:
        r_F -> -r_F.
    """
    b, h1, R = params_from_K(a, K)
    foci_1d = _base_focus_pairs_labeled(a, K)
    base = ellipsoid_base_grid(a=a, K=K, npts=1200)
    d_base = base["d"]

    descs = [
        d for d in _unique_elliptic_descriptors(offsets)
        if _descriptor_has_active_elliptic_segment(d, d_base)
    ]

    points: List[Tuple[float, float, str]] = []
    for x_shift in np.asarray(shifts, dtype=float):
        for _, C, sigma in descs:
            C = float(C)
            sigma = float(sigma)
            for fx, fy, _ in foci_1d:
                fx_abs = float(fx) + float(x_shift)
                if sigma > 0:
                    r_abs = C + R - float(fy)
                else:
                    r_abs = C - R + float(fy)

                points.append((fx_abs, r_abs))
                if abs(r_abs) > 1e-10:
                    points.append((fx_abs, -r_abs))

    unique = sorted(
        {(round(float(x), 8), round(float(y), 8)) for x, y in points},
        key=lambda p: (p[0], p[1])
    )
    return [(float(x), float(y)) for x, y in unique]


def _curve_matches_visible_boundary(global_axis: np.ndarray,
                                  curve_axis: np.ndarray,
                                  curve_r: np.ndarray,
                                  lo_components: List[np.ndarray],
                                  hi_components: List[np.ndarray],
                                  tol: float = 1.0e-4) -> bool:
    """
    Checks whether the exact elliptical boundary coincides with the VISIBLE post-Merge
    2D boundary on at least a nonzero segment. This is the criterion that
    the foci of the corresponding ellipse should actually be shown.
    """
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

    scale = max(1.0, float(np.nanmax(rr)) if np.any(np.isfinite(rr)) else 1.0)
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
                if span > 0.03 * max(1.0, float(xa.max() - xa.min())):
                    return True
    return False


def _visible_exact_focus_points_for_elliptic_order(a: float,
                                                   K: float,
                                                   h1: float,
                                                   h: float,
                                                   offsets: Sequence[float],
                                                   shifts: np.ndarray,
                                                   global_axis: np.ndarray,
                                                   lo_components: List[np.ndarray],
                                                   hi_components: List[np.ndarray]) -> List[Tuple[float, float, str]]:
    """
    Foci only of those REAL generating ellipses that are actually
    visible on the final 2D boundary after Merge.

    This removes false foci that occurred for internal / completely
    hidden ellipses during row overlaps and intersections.
    """
    b, h1, R = params_from_K(a, K, h1=h1, h=h)
    foci_1d = _base_focus_pairs_labeled(a, K, h1=h1)
    base = ellipsoid_base_grid(a=a, K=K, h1=h1, h=h, npts=1600)
    x_base = np.asarray(base['axis'], dtype=float)
    d_base = np.asarray(base['d'], dtype=float)

    descs = [
        d for d in _unique_elliptic_descriptors(offsets)
        if _descriptor_has_active_elliptic_segment(d, d_base)
    ]

    points: List[Tuple[float, float]] = []

    for x_shift in np.asarray(shifts, dtype=float):
        x_abs = x_base + float(x_shift)
        for _, C, sigma in descs:
            C = float(C)
            sigma = float(sigma)
            r_curve = C + sigma * d_base
            mask = np.isfinite(r_curve) & (r_curve > MERGE_TOL)
            if np.count_nonzero(mask) < 8:
                continue

            if not _curve_matches_visible_boundary(
                global_axis=global_axis,
                curve_axis=x_abs[mask],
                curve_r=r_curve[mask],
                lo_components=lo_components,
                hi_components=hi_components,
            ):
                continue

            for fx, fy, flab in foci_1d:
                fx_abs = float(fx) + float(x_shift)
                if sigma > 0:
                    r_abs = C + R - float(fy)
                else:
                    r_abs = C - R + float(fy)
                points.append((fx_abs, r_abs, flab))
                if abs(r_abs) > 1e-10:
                    points.append((fx_abs, -r_abs, flab))

    unique = sorted({(round(float(x), 8), round(float(y), 8), str(lab)) for x, y, lab in points}, key=lambda p: (p[0], p[1], p[2]))
    return [(float(x), float(y), lab) for x, y, lab in unique]


def _plot_focus_points_2d(ax, points: List[Tuple[float, float, str]], scale_ref: float):
    if not points:
        return
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.scatter(xs, ys, s=20, marker='o', color='darkgreen', zorder=7,
               label='Exact foci of visible real ellipses')
    dx = 0.015 * max(scale_ref, 1.0)
    for x, y, lab in points:
        ha = 'left' if x >= 0 else 'right'
        xt = x + dx if x >= 0 else x - dx
        ax.text(xt, y, lab, fontsize=7, color='darkgreen', va='center', ha=ha)



# ============================================================
# ANNOTATIONS FOR GENERATRICES AND 2ND-ORDER 2D DRAWINGS
# ============================================================

def _annotate_generatrix(ax, a: float, b: float, R: float, typ: str, foci_1d):
    ax.axhline(R, color='black', linestyle='--', linewidth=1.05)
    ax.text(0.98, 0.92, 'Rotation axis y = R', transform=ax.transAxes,
            ha='right', va='top', fontsize=9)
    ax.axvline(0.0, color='black', linestyle='--', linewidth=1.0)
    ax.text(0.5, 0.02, 'Symmetry axis x = 0', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=8)

    # Axis / axes of the foci of the generating ellipses.
    if b > a + 1e-12:
        ax.axvline(-a, color='black', linestyle=':', linewidth=1.0)
        ax.axvline(+a, color='black', linestyle=':', linewidth=1.0)
        ax.text(-a, 1.02 * R, 'Axis of foci', ha='center', va='bottom', fontsize=8)
        ax.text(+a, 1.02 * R, 'Axis of foci', ha='center', va='bottom', fontsize=8)
    else:
        ax.axhline(0.0, color='black', linestyle=':', linewidth=1.0)
        ax.text(0.98, 0.06, 'Axis of foci y = 0', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=8)

    xR = -1.10 * a
    ax.annotate('', xy=(xR, R), xytext=(xR, 0.0),
                arrowprops=dict(arrowstyle='<->', lw=1.15, color='black'))
    ax.text(xR - 0.08 * a, 0.5 * R, 'R', ha='right', va='center', fontsize=12)

    for fx, fy, fl in foci_1d:
        ax.plot(fx, fy, 'D', color='red', ms=4, zorder=10)
        ax.annotate(fl, (fx, fy), fontsize=5, color='red',
                    xytext=(3, 3), textcoords='offset points')


def _annotate_order2_2d(ax, a: float, R: float, x_left: float):
    ax.axhline(0.0, color='black', linestyle='--', linewidth=1.05)
    ax.text(0.98, 0.52, 'Rotation axis / symmetry axis r = 0',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
    ax.axvline(0.0, color='black', linestyle='--', linewidth=1.0, alpha=0.8)
    ax.text(0.5, 0.03, 'Reference axis x = 0', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=8)

    ax.axhline(R, color='black', linestyle='--', linewidth=1.0)
    ax.axhline(-R, color='black', linestyle='--', linewidth=1.0)
    ax.text(0.98, 0.88, 'Boundary R', transform=ax.transAxes, ha='right', va='top', fontsize=8)

    xR = x_left
    ax.annotate('', xy=(xR, R), xytext=(xR, 0.0),
                arrowprops=dict(arrowstyle='<->', lw=1.15, color='black'))
    ax.text(xR - 0.06 * max(a, 1.0), 0.5 * R, 'R', ha='right', va='center', fontsize=12)



# ============================================================
# BASE GENERATRICES
# ============================================================

def plot_base_generatrix(a: float,
                         K: float,
                         h1: float = 0.0,
                         h: float = 0.0,
                         outpath: Optional[Path] = None,
                         show: bool = True,
                         npts: int = 1200):
    typ = classify_K(K)
    b, h1, R = params_from_K(a, K, h1=h1, h=h)
    x, y, yL, yR, *_ = build_signed_profile(a, b, h1, n=npts)
    _, foci_1d, _ = compute_foci(a, b, h1, R)

    fig, ax = plt.subplots(figsize=(8.6, 6.2))

    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(-a + a * np.cos(t), b * np.sin(t), '--', color='gray', alpha=0.28, lw=0.8)
    ax.plot((a + h1) + a * np.cos(t), b * np.sin(t), '--', color='gray', alpha=0.28, lw=0.8)

    vmL = ~np.isnan(yL)
    vmR = ~np.isnan(yR)
    if np.any(vmL):
        ax.plot(x[vmL], yL[vmL], '-', color='steelblue', lw=2.0)
    if np.any(vmR):
        ax.plot(x[vmR], yR[vmR], '-', color='darkorange', lw=2.0)
    vm = ~np.isnan(y)
    ax.plot(x[vm], y[vm], '-', color='#2ca02c', lw=2.8)

    _annotate_generatrix(ax, a=a, b=b, R=R, typ=typ, foci_1d=foci_1d)

    ax.set_xlim(-1.25 * a, 1.25 * a)
    ax.set_ylim(min(-0.15 * max(a, b), -0.25), 1.22 * R)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Pseudoellipsoid of the 2nd order: base generatrix; {typ} type; K={K:g}, a={a:g}, b={b:g}, h1={h1:g}, h={h:g}')
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


# ============================================================
# 2D SECTIONS
# ============================================================

def plot_union_section(order: int,
                       a: float,
                       K: float,
                       h1: float,
                       h: float,
                       global_axis: np.ndarray,
                       lo_components: List[np.ndarray],
                       hi_components: List[np.ndarray],
                       shifts: np.ndarray,
                       offsets: Sequence[float],
                       m: int,
                       h_row: float,
                       outpath: Optional[Path] = None,
                       show: bool = True):
    typ = classify_K(K)
    b, h1, R = params_from_K(a, K, h1=h1, h=h)

    fig, ax = plt.subplots(figsize=(9.5, 7.0))
    for lo, hi in zip(lo_components, hi_components):
        _plot_component_2d(ax, global_axis, lo, hi)

    rmax, amin, amax = _axis_limits_from_components(global_axis, hi_components)

    # Foci are shown for all orders, but only as EXACT foci
    # of real generating elliptical boundaries of the recursion. Merge creates
    # no new foci, so only descriptors C ± d(x) are used.
    focus_points = _visible_exact_focus_points_for_elliptic_order(
        a=a,
        K=K,
        h1=h1,
        h=h,
        offsets=offsets,
        shifts=shifts,
        global_axis=global_axis,
        lo_components=lo_components,
        hi_components=hi_components,
    )
    _plot_focus_points_2d(ax, focus_points, scale_ref=max(rmax, amax - amin))

    if order == 2:
        _annotate_order2_2d(ax, a=a, R=R, x_left=amin - 0.10 * (amax - amin + 1e-9))
    else:
        ax.axhline(0.0, color='black', linestyle='--', linewidth=0.9, alpha=0.75)

    ax.set_xlim(amin - 0.10 * (amax - amin + 1e-9), amax + 0.10 * (amax - amin + 1e-9))
    ax.set_ylim(-1.08 * rmax, 1.08 * rmax)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x-axis')
    ax.set_ylabel('Radial coordinate')
    ax.set_title(
        f'Pseudoellipsoid: {typ} type; order {order}; K={K:g}\n'
        f'a={a:g}, b={b:g}, h1={h1:g}, h={h:g}, R={R:g}, offsets={list(offsets)}, m={m}, h_row={h_row:g}'
    )
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


# ============================================================
# 3D SURFACES
# ============================================================

def revolve_about_x(axis_coord: np.ndarray, radius: np.ndarray, nphi: int = NPHI):
    phi = np.linspace(0.0, 2.0 * np.pi, nphi)
    X, Phi = np.meshgrid(axis_coord, phi, indexing='ij')
    RR = np.tile(radius[:, None], (1, phi.size))
    Y = RR * np.cos(Phi)
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
        ax.set_proj_type('ortho')
    except Exception:
        pass


def _surface_subsegments(axis_seg: np.ndarray, radius_seg: np.ndarray) -> List[Tuple[np.ndarray, np.ndarray]]:
    axis_seg = np.asarray(axis_seg, dtype=float)
    radius_seg = np.asarray(radius_seg, dtype=float)
    if axis_seg.size < 2:
        return []
    out = []
    for a0, b0 in _split_on_artificial_jumps(radius_seg):
        if b0 - a0 >= 2:
            out.append((axis_seg[a0:b0], radius_seg[a0:b0]))
    return out


def _plot_surface_piece(ax, axis_seg: np.ndarray, radius_seg: np.ndarray, bounds: PlotBounds3D, stride: int = 3):
    X, Y, Z = revolve_about_x(axis_seg, radius_seg)
    ax.plot_surface(
        X[::stride, ::stride], Y[::stride, ::stride], Z[::stride, ::stride],
        linewidth=0, edgecolor='none', alpha=SURFACE_ALPHA,
        antialiased=False, shade=True,
    )
    bounds.update(X, Y, Z)


def _plot_surface_boundary(ax, axis_seg: np.ndarray, radius_seg: np.ndarray, bounds: PlotBounds3D, stride: int = 3):
    for aa, rr in _surface_subsegments(axis_seg, radius_seg):
        _plot_surface_piece(ax, aa, rr, bounds, stride=stride)


def _plot_bright_2d_section_on_3d(ax, global_axis: np.ndarray, lo_components: List[np.ndarray], hi_components: List[np.ndarray]):
    y_plane = 0.0
    y_eps = 1.0e-6  # micro-shift for stable display of internal lines without noticeable deviation from the rotation axis
    outer_color = 'magenta'
    inner_color = 'lime'
    for lo, hi in zip(lo_components, hi_components):
        mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
        for st, en in _contiguous_segments(mask, global_axis):
            xx = global_axis[st:en]
            lo_seg = lo[st:en]
            hi_seg = hi[st:en]
            for a0, b0 in _split_on_artificial_jumps(hi_seg):
                x = xx[a0:b0]
                z = hi_seg[a0:b0]
                if x.size >= 2:
                    for yp in (y_plane, y_eps):
                        y = np.full_like(x, yp)
                        ax.plot(x, y, z, color=outer_color, linewidth=3.4, alpha=1.0)
                        ax.plot(x, y, -z, color=outer_color, linewidth=3.4, alpha=1.0)
            if np.nanmax(lo_seg) > MERGE_TOL:
                for a0, b0 in _split_on_artificial_jumps(lo_seg):
                    x = xx[a0:b0]
                    z = lo_seg[a0:b0]
                    if x.size >= 2:
                        for yp in (y_plane, y_eps):
                            y = np.full_like(x, yp)
                            ax.plot(x, y, z, color=inner_color, linewidth=3.1, alpha=1.0)
                            ax.plot(x, y, -z, color=inner_color, linewidth=3.1, alpha=1.0)



def _plot_all_generating_elliptic_profiles_on_3d(ax,
                                                 a: float,
                                                 K: float,
                                                 h1: float,
                                                 h: float,
                                                 offsets: Sequence[float],
                                                 m: int,
                                                 h_row: float,
                                                 npts: int = 1400):
    """
    Additional control visualization: shows ALL generating
    elliptical profiles of this order before Merge leaves only
    the outer boundary of the union volume.

    This does not change the geometry of the volume. It is a bright check overlay so that on
    3D 4th-order plots the internal elliptical generatrices remain visible.
    """
    base = ellipsoid_base_grid(a=a, K=K, h1=h1, h=h, npts=npts)
    x_base = np.asarray(base["axis"], dtype=float)
    d_base = np.asarray(base["d"], dtype=float)
    shifts = stack_shifts_by_width(float(base["width"]), h_row, m)
    descs = [
        d for d in _unique_elliptic_descriptors(offsets)
        if _descriptor_has_active_elliptic_segment(d, d_base)
    ]

    y_plane = 0.0
    y_eps = 2.0e-6
    color = "cyan"

    for x_shift in np.asarray(shifts, dtype=float):
        x = x_base + float(x_shift)
        for _, C, sigma in descs:
            r = np.maximum(float(C) + float(sigma) * d_base, 0.0)
            mask = np.isfinite(r) & (r > MERGE_TOL)
            for st, en in _contiguous_segments(mask, x):
                xx = x[st:en]
                rr = r[st:en]
                for a0, b0 in _split_on_artificial_jumps(rr):
                    if b0 - a0 >= 2:
                        xs = xx[a0:b0]
                        zs = rr[a0:b0]
                        for yp in (y_plane, y_eps):
                            yy = np.full_like(xs, yp)
                            ax.plot(xs, yy, zs, color=color, linewidth=1.9, alpha=1.0)
                            ax.plot(xs, yy, -zs, color=color, linewidth=1.9, alpha=1.0)


def plot_union_surface(order: int,
                       a: float,
                       K: float,
                       h1: float,
                       h: float,
                       global_axis: np.ndarray,
                       lo_components: List[np.ndarray],
                       hi_components: List[np.ndarray],
                       offsets: Sequence[float],
                       m: int,
                       h_row: float,
                       outpath: Optional[Path] = None,
                       show: bool = True,
                       stride: int = 3):
    typ = classify_K(K)
    b, h1, R = params_from_K(a, K, h1=h1, h=h)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    bounds = PlotBounds3D()

    for lo, hi in zip(lo_components, hi_components):
        mask = np.isfinite(lo) & np.isfinite(hi) & (hi > lo + MERGE_TOL)
        for st, en in _contiguous_segments(mask, global_axis):
            axis_seg = global_axis[st:en]
            lo_seg = lo[st:en]
            hi_seg = hi[st:en]
            _plot_surface_boundary(ax, axis_seg, hi_seg, bounds, stride=stride)
            if np.nanmax(lo_seg) > MERGE_TOL:
                _plot_surface_boundary(ax, axis_seg, lo_seg, bounds, stride=stride)

    _plot_bright_2d_section_on_3d(ax, global_axis, lo_components, hi_components)
    _plot_all_generating_elliptic_profiles_on_3d(ax, a=a, K=K, h1=h1, h=h, offsets=offsets, m=m, h_row=h_row)
    set_axes_equal_real_3d(ax, bounds)
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(
        f'Pseudoellipsoid: {typ} type; order {order}; K={K:g}\n'
        f'a={a:g}, b={b:g}, h1={h1:g}, h={h:g}, R={R:g}, offsets={list(offsets)}, m={m}, h_row={h_row:g}'
    )
    fig.tight_layout()
    _save_or_show(fig, outpath, show)


# ============================================================
# MAIN RUN
# ============================================================

def run_pseudoellipsoids(a: float = 1.0,
                         K_horizontal: float = 0.5,
                         K_vertical: float = 1.5,
                         h1: float = 0.0,
                         h: float = 0.0,
                         offsets: Sequence[float] = (0.8, 1.6),
                         mode: str = 'all',
                         outdir: str | None = None,
                         npts: int = 700,
                         m: int = 3,
                         h_row: float = -0.25,
                         show: bool = True,
                         make_base_generatrices: bool = True,
                         save_files: bool = False) -> Dict[str, object]:
    """
    Builds the same set of constructions for the two pseudoellipsoid types:
    horizontal (K<1) and vertical (K>1).

    Notebook mode:
    - by default show=True, figures are displayed directly in the cell;
    - by default save_files=False, separate PNG files are not created;
    - to save PNG files, set save_files=True and specify outdir.

    h1 and h work as in the source script: build_signed_profile(a,b,h1), R=b+h.
    offsets of length 2 -> orders 2, 3, 4.
    """
    if mode not in {'section', 'surface', 'all'}:
        raise ValueError("mode must be 'section', 'surface' or 'all'")

    classify_K(K_horizontal)
    classify_K(K_vertical)
    if not (K_horizontal < 1.0 and K_vertical > 1.0):
        raise ValueError('Need K_horizontal < 1 and K_vertical > 1')

    offsets = list(float(v) for v in offsets)
    if len(offsets) != 2:
        raise ValueError('For consistency with the other scripts, two offsets are expected: for the 3rd and 4th orders')

    outdir_path = Path(outdir) if outdir is not None else None
    if save_files:
        if outdir_path is None:
            outdir_path = Path('pseudoellipsoids_output')
        _ensure_dir(outdir_path)

    saved_files: List[str] = []

    cases = [
        ('horizontal', K_horizontal),
        ('vertical', K_vertical),
    ]

    if make_base_generatrices:
        for name, K in cases:
            path = (outdir_path / f'{name}_base_generatrix_K{K:g}.png') if save_files else None
            plot_base_generatrix(a=a, K=K, h1=h1, h=h, outpath=path, show=show, npts=max(1200, npts))
            if path is not None:
                saved_files.append(str(path))

    for name, K in cases:
        for order in (2, 3, 4):
            used_offsets = offsets[:max(0, order - 2)]
            global_axis, lo_components, hi_components, shifts, base = ellipsoid_order_union_components(
                a=a, K=K, h1=h1, h=h, offsets=used_offsets, npts=npts, m=m, h_row=h_row
            )
            prefix = f'{name}_order_{order:02d}_m{m}_h1{h1}_h{h}_hrow{h_row}'

            if mode in {'section', 'all'}:
                path = (outdir_path / f'{prefix}_section.png') if save_files else None
                plot_union_section(
                    order=order, a=a, K=K, h1=h1, h=h,
                    global_axis=global_axis,
                    lo_components=lo_components,
                    hi_components=hi_components,
                    shifts=shifts,
                    offsets=used_offsets,
                    m=m, h_row=h_row,
                    outpath=path, show=show,
                )
                if path is not None:
                    saved_files.append(str(path))

            if mode in {'surface', 'all'}:
                path = (outdir_path / f'{prefix}_surface.png') if save_files else None
                plot_union_surface(
                    order=order, a=a, K=K, h1=h1, h=h,
                    global_axis=global_axis,
                    lo_components=lo_components,
                    hi_components=hi_components,
                    offsets=used_offsets,
                    m=m, h_row=h_row,
                    outpath=path, show=show,
                )
                if path is not None:
                    saved_files.append(str(path))

    return {
        'a': float(a),
        'K_horizontal': float(K_horizontal),
        'K_vertical': float(K_vertical),
        'h1': float(h1),
        'h': float(h),
        'offsets': offsets,
        'm': int(m),
        'h_row': float(h_row),
        'saved_files': saved_files,
    }


if __name__ == '__main__':
    # Default notebook call: figures are displayed directly in the cell,
    # without mandatory saving to separate files.
    run_pseudoellipsoids(
        a=1.0,
        K_horizontal=0.5,
        K_vertical=1.5,
        h1=0.0,
        h=0.0,
        offsets=(0.8, 1.6),
        mode='all',
        outdir=None,
        npts=700,
        m=3,
        h_row=0.0,
        show=True,
        make_base_generatrices=True,
        save_files=False,
    )

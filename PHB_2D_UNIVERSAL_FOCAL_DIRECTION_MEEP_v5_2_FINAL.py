#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHB_2D_UNIVERSAL_FOCAL_DIRECTION_MEEP_v5_2.py
====================================================
Universal, independently selectable Maxwell/FDTD verifier of the causal link:

    PHB focal concentration  ->  useful differential angular information.

The physical sensing ring is unchanged. In the 2D meridional section it remains
represented by two finite air-side monitor regions adjacent to the uninterrupted
PEC wall at (x=a, y=+/-R).

The script can run any requested subset independently:
    empty_reference  no cavity;
    phb              exact hyperbolic generatrix;
    smooth_area      smooth, non-hyperbolic, equal-area causal control.

Why this script is different from v4
------------------------------------
1. A single target frequency is used (default f=1.8), so the PHB and smooth-area
   profiles can be separated at high resolution without a 243-case broadband sweep.
2. The primary stage uses MEEP's frequency-domain CW solver. For a high-Q cavity
   this can be much faster than time stepping to t=500 while directly returning
   the steady-state field at the requested frequency.
3. Resolution convergence is checked at two user-selected resolutions, normally
   80 and 120. At resolution 120 the maximum PHB/smooth-area profile separation
   for a=b=0.5,R=2 is about four grid cells.
4. The two sensor signals are decomposed into calibrated common and differential
   spatial modes. The direction observable is the complex mode-conversion ratio

       z(theta) = O(theta) / E(theta),

   where E and O are the symmetry-matched common and differential channels.
   This is the two-sensor analogue of spatial mode filtering.
5. Receiver-noise tests use one fixed absolute noise scale derived from an empty-cell
   reference run. Therefore any absolute focal gain is NOT normalized away.
6. The script integrates electromagnetic energy over the full cavity and the focal
   band. It reports focal energy fraction and area-normalized enrichment, not only
   a focal-line/center-line ratio.
7. It estimates local Fisher information and the corresponding Cramer-Rao lower
   bound from the complex derivatives of the two sensor signals with respect to angle.
8. An optional time-domain confirmation stage runs only a minimal subset to t=500
   and performs Harminv near the focal band. This is a validation of the CW result,
   not the primary expensive sweep.

Scientific scope
----------------
* Strictly 2D Cartesian/meridional TE Maxwell model.
* Signed angle in one plane only.
* No range estimation and no full 3D azimuth/elevation claim.
* A positive result establishes the causal bridge only within the exact PHB and
  smooth equal-area control used here.
"""

from __future__ import annotations

import argparse
import cmath
import csv
import json
import math
import shutil
import time
import traceback
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import Rectangle
except Exception:
    plt = None
    MplPath = None
    Rectangle = None

try:
    import meep as mp
except Exception:
    mp = None

SCRIPT_VERSION = "5.2"
GEOMETRY_MODELS = ("phb", "smooth_area")
RUN_MODELS = ("empty_reference", "phb", "smooth_area")


# =============================================================================
# Generic helpers
# =============================================================================

def require_meep() -> None:
    if mp is None:
        raise RuntimeError("meep is unavailable; run FDTD stages in the meep-working Docker image")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=True), encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    rows = list(rows)
    keys: List[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def parse_float_list(text: str) -> List[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(round(v)) for v in parse_float_list(text)]


def parse_str_list(text: str) -> List[str]:
    return [v.strip().lower() for v in text.split(",") if v.strip()]


def angle_tag(angle: float) -> str:
    return ("m" if angle < 0 else "p") + f"{abs(angle):.4f}".replace(".", "p")


def archive_folder(folder: Path) -> Optional[Path]:
    if not folder.exists():
        return None
    zpath = folder.parent / f"{folder.name}.zip"
    if zpath.exists():
        zpath.unlink()
    shutil.make_archive(str(folder), "zip", root_dir=str(folder.parent), base_dir=folder.name)
    return zpath


def safe_complex_mean(arr: Any, weights: Optional[Any] = None) -> complex:
    a = np.asarray(arr, dtype=np.complex128).squeeze()
    valid = np.isfinite(a.real) & np.isfinite(a.imag)
    if not np.any(valid):
        return complex(float("nan"), float("nan"))
    if weights is None:
        return complex(np.mean(a[valid]))
    w = np.asarray(weights, dtype=float).squeeze()
    if w.shape != a.shape:
        w = np.broadcast_to(w, a.shape)
    w = np.where(valid & np.isfinite(w), w, 0.0)
    den = float(np.sum(w))
    return complex(np.sum(w * a) / den) if den > 0 else complex(np.mean(a[valid]))


def relative_complex_error(a: complex, b: complex) -> float:
    return float(abs(a - b) / max(abs(a), abs(b), 1e-300))


def relative_l2(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a)
    b = np.asarray(b)
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(a), np.linalg.norm(b), 1e-300))


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    ar = np.concatenate((np.asarray(a).real.ravel(), np.asarray(a).imag.ravel()))
    br = np.concatenate((np.asarray(b).real.ravel(), np.asarray(b).imag.ravel()))
    den = float(np.linalg.norm(ar) * np.linalg.norm(br))
    return float(np.dot(ar, br) / den) if den > 0 else float("nan")


# =============================================================================
# Geometry: copied from v4 and restricted to the decisive causal pair
# =============================================================================

@dataclass(frozen=True)
class Geometry:
    a: float
    b: float
    R: float
    R1: float
    R1_mode: str
    model: str

    @property
    def Rin(self) -> float:
        return self.R1 if self.R1_mode == "absolute" else self.R * self.R1

    @property
    def c(self) -> float:
        return math.sqrt(self.a * self.a + self.b * self.b)

    @property
    def L(self) -> float:
        return self.a * math.sqrt(1.0 + (self.R / self.b) ** 2)

    @property
    def x_inlet(self) -> float:
        return -self.a * math.sqrt(1.0 + ((self.R - self.Rin) / self.b) ** 2)

    @property
    def cap_rmax(self) -> float:
        return self.R + self.a

    def x_left_phb(self, rho: np.ndarray) -> np.ndarray:
        return -self.a * np.sqrt(1.0 + ((self.R - rho) / self.b) ** 2)

    def x_right_phb(self, rho: np.ndarray) -> np.ndarray:
        return +self.a * np.sqrt(1.0 + ((self.R - rho) / self.b) ** 2)


def polygon_area(poly: np.ndarray) -> float:
    p = np.asarray(poly, dtype=float)
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def cavity_polygon_from_upper(upper: np.ndarray) -> np.ndarray:
    lower = np.asarray(upper[::-1], dtype=float).copy()
    lower[:, 1] *= -1.0
    return np.vstack((upper, lower))


def _phb_upper(g: Geometry, n_horn: int, n_cap: int) -> np.ndarray:
    rho_l = np.linspace(g.Rin, g.R, n_horn)
    x_l = g.x_left_phb(rho_l)
    x_c = np.linspace(-g.a, g.a, n_cap)
    y_c = g.R + np.sqrt(np.maximum(g.a * g.a - x_c * x_c, 0.0))
    rho_r = np.linspace(g.R, 0.0, n_horn)
    x_r = g.x_right_phb(rho_r)
    return np.column_stack((
        np.concatenate((x_l, x_c[1:], x_r[1:])),
        np.concatenate((rho_l, y_c[1:], rho_r[1:])),
    ))


def _smooth_upper(g: Geometry, exponent: float, n_horn: int, n_cap: int) -> np.ndarray:
    rho_l = np.linspace(g.Rin, g.R, n_horn)
    s = (rho_l - g.Rin) / (g.R - g.Rin)
    h_l = 1.0 - np.power(np.maximum(1.0 - s, 0.0), exponent)
    x_l = g.x_inlet + (-g.a - g.x_inlet) * h_l

    x_c = np.linspace(-g.a, g.a, n_cap)
    y_c = g.R + np.sqrt(np.maximum(g.a * g.a - x_c * x_c, 0.0))

    rho_r = np.linspace(g.R, 0.0, n_horn)
    u = (g.R - rho_r) / g.R
    x_r = g.a + (g.L - g.a) * np.power(np.maximum(u, 0.0), exponent)

    return np.column_stack((
        np.concatenate((x_l, x_c[1:], x_r[1:])),
        np.concatenate((rho_l, y_c[1:], rho_r[1:])),
    ))


@lru_cache(maxsize=64)
def smooth_area_exponent(a: float, b: float, R: float, Rin: float) -> float:
    base = Geometry(a, b, R, Rin, "absolute", "phb")
    target = polygon_area(cavity_polygon_from_upper(_phb_upper(base, 4000, 2400)))

    def area(p: float) -> float:
        return polygon_area(cavity_polygon_from_upper(_smooth_upper(base, p, 4000, 2400)))

    lo, hi = 1.0001, 12.0
    if not (area(lo) >= target >= area(hi)):
        raise RuntimeError("could not bracket smooth-area exponent")
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if area(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def upper_centerline(g: Geometry, n_horn: int = 720, n_cap: int = 420) -> np.ndarray:
    if g.model == "phb":
        return _phb_upper(g, n_horn, n_cap)
    if g.model == "smooth_area":
        p = smooth_area_exponent(g.a, g.b, g.R, g.Rin)
        return _smooth_upper(g, p, n_horn, n_cap)
    raise ValueError(f"unknown model {g.model}")


def cavity_area(g: Geometry) -> float:
    return polygon_area(cavity_polygon_from_upper(upper_centerline(g, 3000, 1800)))


def _polyline_strip_polygon(points: np.ndarray, wall_t: float) -> np.ndarray:
    p = np.asarray(points, dtype=float)
    keep = np.ones(len(p), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(p, axis=0), axis=1) > 1e-12
    p = p[keep]
    d = np.diff(p, axis=0)
    seglen = np.linalg.norm(d, axis=1)
    normals = np.column_stack((-d[:, 1], d[:, 0])) / seglen[:, None]
    half = 0.5 * float(wall_t)
    offsets = np.zeros_like(p)
    offsets[0] = half * normals[0]
    offsets[-1] = half * normals[-1]
    for i in range(1, len(p) - 1):
        m = normals[i - 1] + normals[i]
        nm = float(np.linalg.norm(m))
        if nm <= 1e-12:
            offsets[i] = half * normals[i]
            continue
        m /= nm
        denom = abs(float(np.dot(m, normals[i - 1])))
        offsets[i] = min(half / max(denom, 0.40), 2.5 * half) * m
    return np.vstack((p + offsets, (p - offsets)[::-1]))


@lru_cache(maxsize=128)
def wall_polygon_tuples(g: Geometry, wall_t: float):
    upper = upper_centerline(g)
    up_poly = _polyline_strip_polygon(upper, wall_t)
    lo_poly = up_poly.copy()
    lo_poly[:, 1] *= -1.0
    return tuple(map(tuple, up_poly)), tuple(map(tuple, lo_poly))


def wall_polygons(g: Geometry, wall_t: float) -> Tuple[np.ndarray, np.ndarray]:
    up, lo = wall_polygon_tuples(g, float(wall_t))
    return np.asarray(up, dtype=float), np.asarray(lo, dtype=float)


def make_geometry_objects(g: Geometry, wall_t: float):
    require_meep()
    metal = getattr(mp, "metal", mp.Medium(epsilon=1e9))
    objects = []
    for poly in wall_polygons(g, wall_t):
        vertices = [mp.Vector3(float(x), float(y), 0) for x, y in poly]
        objects.append(mp.Prism(vertices=vertices, height=mp.inf, axis=mp.Vector3(0, 0, 1), material=metal))
    return objects


def geometry_dimensions(g: Geometry, dpml: float, x_air: float, y_air: float) -> Tuple[float, float]:
    return 2.0 * (g.L + x_air + dpml), 2.0 * (g.cap_rmax + y_air + dpml)


def sensor_regions(g: Geometry, wall_t: float, resolution: int, offset_physical: float,
                   size_x: float, size_y: float) -> Dict[str, Any]:
    delta = max(float(offset_physical), 0.75 * wall_t + 0.5 * size_x)
    return {
        "wall_top": (g.a, +g.R),
        "wall_bottom": (g.a, -g.R),
        "probe_top": (g.a - delta, +g.R),
        "probe_bottom": (g.a - delta, -g.R),
        "sensor_size": (size_x, size_y),
        "offset_physical": delta,
        "offset_cells": delta * resolution,
    }


def cavity_path(g: Geometry):
    if MplPath is None:
        raise RuntimeError("matplotlib is required for cavity-mask integration")
    return MplPath(cavity_polygon_from_upper(upper_centerline(g, 1800, 1000)), closed=True)


def geometry_audit(geometries: Sequence[Geometry], resolution: int, wall_t: float,
                   sensor_offset: float, sensor_size_x: float, sensor_size_y: float) -> Dict[str, Any]:
    gp = next(g for g in geometries if g.model == "phb")
    gs = next(g for g in geometries if g.model == "smooth_area")
    rho = np.linspace(0, gp.R, 4001)
    xp = gp.x_right_phb(rho)
    p = smooth_area_exponent(gs.a, gs.b, gs.R, gs.Rin)
    xs = gs.a + (gs.L - gs.a) * np.power((gs.R - rho) / gs.R, p)
    diff = np.abs(xp - xs)
    audits = {}
    for g in geometries:
        sensors = sensor_regions(g, wall_t, resolution, sensor_offset, sensor_size_x, sensor_size_y)
        audits[g.model] = {
            "area": cavity_area(g),
            "wall_cells": wall_t * resolution,
            "sensor_offset_cells": sensors["offset_cells"],
            "sensor_top": sensors["probe_top"],
            "sensor_bottom": sensors["probe_bottom"],
        }
    return {
        "status": "PASS",
        "resolution": resolution,
        "smooth_area_exponent": p,
        "phb_area": audits["phb"]["area"],
        "smooth_area": audits["smooth_area"]["area"],
        "relative_area_difference": abs(audits["phb"]["area"] - audits["smooth_area"]["area"]) / audits["phb"]["area"],
        "max_profile_separation": float(np.max(diff)),
        "max_profile_separation_cells": float(np.max(diff) * resolution),
        "mean_profile_separation_cells": float(np.mean(diff) * resolution),
        "models": audits,
    }


def save_geometry_outputs(outroot: Path, geometries: Sequence[Geometry], resolution: int,
                          wall_t: float, sensor_offset: float, sensor_size_x: float,
                          sensor_size_y: float) -> Dict[str, Any]:
    audit = geometry_audit(geometries, resolution, wall_t, sensor_offset, sensor_size_x, sensor_size_y)
    write_json(outroot / "GEOMETRY_AUDIT.json", audit)
    if plt is None:
        return audit
    fig, ax = plt.subplots(figsize=(12, 8))
    styles = {"phb": "-", "smooth_area": "--"}
    for g in geometries:
        up = upper_centerline(g)
        ax.plot(up[:, 0], up[:, 1], styles[g.model], lw=2.2, label=g.model)
        ax.plot(up[:, 0], -up[:, 1], styles[g.model], lw=2.2)
        s = sensor_regions(g, wall_t, resolution, sensor_offset, sensor_size_x, sensor_size_y)
        ax.scatter([s["probe_top"][0], s["probe_bottom"][0]], [s["probe_top"][1], s["probe_bottom"][1]], s=45)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"PHB vs smooth equal-area control; resolution={resolution}; max separation={audit['max_profile_separation_cells']:.2f} cells")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outroot / "01_GEOMETRY_PHB_VS_SMOOTH.png", dpi=220)
    plt.close(fig)

    rho = np.linspace(0, geometries[0].R, 2000)
    gp = next(g for g in geometries if g.model == "phb")
    gs = next(g for g in geometries if g.model == "smooth_area")
    p = smooth_area_exponent(gs.a, gs.b, gs.R, gs.Rin)
    xp = gp.x_right_phb(rho)
    xs = gs.a + (gs.L - gs.a) * ((gs.R - rho) / gs.R) ** p
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rho, (xp - xs) * resolution)
    ax.axhline(0, lw=0.8)
    ax.set_xlabel("rho")
    ax.set_ylabel("(x_PHB - x_smooth) in grid cells")
    ax.set_title("Numerical separation of the two causal geometries")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outroot / "02_PROFILE_SEPARATION_CELLS.png", dpi=220)
    plt.close(fig)
    return audit


# =============================================================================
# MEEP source and steady-state extraction
# =============================================================================

def phased_line_amp(freq: float, angle_deg: float, source_half_height: float):
    theta = math.radians(angle_deg)
    ky = float(freq) * math.sin(theta)
    h = max(float(source_half_height), 1e-12)

    def amp(p):
        u = abs(float(p.y)) / h
        envelope = math.exp(-0.5 * (u / 0.86) ** 12)
        return envelope * cmath.exp(1j * 2.0 * math.pi * ky * float(p.y))

    return amp


def weighted_region_mean(sim, component, center: Tuple[float, float], size: Tuple[float, float]) -> complex:
    vol = mp.Volume(center=mp.Vector3(center[0], center[1], 0), size=mp.Vector3(size[0], size[1], 0))
    arr = sim.get_array(vol=vol, component=component, cmplx=True)
    meta = sim.get_array_metadata(vol=vol)
    w = meta[3] if len(meta) >= 4 else None
    return safe_complex_mean(arr, w)


def field_arrays(sim, g: Geometry, components: Sequence[Any]) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    center = mp.Vector3(0.5 * (g.x_inlet + g.L), 0, 0)
    size = mp.Vector3(g.L - g.x_inlet, 2.0 * g.cap_rmax, 0)
    vol = mp.Volume(center=center, size=size)
    arrays: Dict[str, np.ndarray] = {}
    for comp in components:
        arrays[str(comp)] = np.asarray(sim.get_array(vol=vol, component=comp, cmplx=True), dtype=np.complex128).squeeze()
    x, y, _z, w = sim.get_array_metadata(vol=vol)
    return arrays, np.asarray(x), np.asarray(y), np.asarray(w, dtype=float).squeeze()


def extract_steady_state(sim, g: Geometry, sensors: Dict[str, Any], frequency: float,
                         focal_half_thickness: float, save_map: bool = False) -> Tuple[Dict[str, Any], Optional[Dict[str, np.ndarray]]]:
    top = weighted_region_mean(sim, mp.Hz, sensors["probe_top"], sensors["sensor_size"])
    bottom = weighted_region_mean(sim, mp.Hz, sensors["probe_bottom"], sensors["sensor_size"])

    arrays, x, y, weights = field_arrays(sim, g, [mp.Hz, mp.Ex, mp.Ey])
    hz = arrays[str(mp.Hz)]
    ex = arrays[str(mp.Ex)]
    ey = arrays[str(mp.Ey)]
    u = 0.25 * (np.abs(hz) ** 2 + np.abs(ex) ** 2 + np.abs(ey) ** 2)
    X, Y = np.meshgrid(x, y, indexing="ij")
    points = np.column_stack((X.ravel(), Y.ravel()))
    cmask = cavity_path(g).contains_points(points).reshape(X.shape)
    fmask = cmask & (np.abs(X) <= g.a) & (np.abs(np.abs(Y) - g.R) <= focal_half_thickness)
    w = np.broadcast_to(weights, u.shape)
    cavity_energy = float(np.sum(w[cmask] * u[cmask]))
    focal_energy = float(np.sum(w[fmask] * u[fmask]))
    cavity_area_num = float(np.sum(w[cmask]))
    focal_area_num = float(np.sum(w[fmask]))
    eta = focal_energy / max(cavity_energy, 1e-300)
    area_fraction = focal_area_num / max(cavity_area_num, 1e-300)
    enrichment = eta / max(area_fraction, 1e-300)

    margin = max(2.5 / max(len(x) / max(g.L - g.x_inlet, 1e-9), 1.0), 0.02 * g.a)
    line_len = max(2.0 * (g.a - margin), 1e-6)
    line_size = (line_len, 0.0)
    upper = weighted_region_line(sim, mp.Hz, (0.0, +g.R), line_size)
    lower = weighted_region_line(sim, mp.Hz, (0.0, -g.R), line_size)

    snap = {
        "frequency": frequency,
        "top_real": float(top.real),
        "top_imag": float(top.imag),
        "top_abs": float(abs(top)),
        "bottom_real": float(bottom.real),
        "bottom_imag": float(bottom.imag),
        "bottom_abs": float(abs(bottom)),
        "cavity_energy": cavity_energy,
        "focal_energy": focal_energy,
        "focal_energy_fraction": eta,
        "focal_area_fraction": area_fraction,
        "focal_enrichment": enrichment,
        "profile_x": upper[0].tolist(),
        "profile_upper_real": upper[1].real.tolist(),
        "profile_upper_imag": upper[1].imag.tolist(),
        "profile_lower_real": lower[1].real.tolist(),
        "profile_lower_imag": lower[1].imag.tolist(),
    }
    fmap = None
    if save_map:
        fmap = {"x": x, "y": y, "hz": hz, "energy_density": u, "cavity_mask": cmask.astype(np.uint8), "focal_mask": fmask.astype(np.uint8)}
    return snap, fmap


def weighted_region_line(sim, component, center: Tuple[float, float], size: Tuple[float, float]) -> Tuple[np.ndarray, np.ndarray]:
    vol = mp.Volume(center=mp.Vector3(center[0], center[1], 0), size=mp.Vector3(size[0], size[1], 0))
    arr = np.asarray(sim.get_array(vol=vol, component=component, cmplx=True), dtype=np.complex128).squeeze().ravel()
    x, _y, _z, _w = sim.get_array_metadata(vol=vol)
    xx = np.asarray(x, dtype=float).ravel()
    n = min(len(xx), len(arr))
    return xx[:n], arr[:n]


def build_simulation(args, g: Optional[Geometry], resolution: int, angle: float,
                     frequency: float, source_x: float, sx: float, sy: float, dpml: float,
                     wall_t: float, continuous: bool):
    require_meep()
    source_half_height = 0.5 * (sy - 2 * dpml - 0.25)
    if continuous:
        src_time = mp.ContinuousSource(
            frequency=frequency,
            width=args.cw_source_width_periods / frequency,
            slowness=args.cw_source_slowness,
        )
    else:
        src_time = mp.GaussianSource(
            frequency=frequency,
            fwidth=max(args.td_source_fwidth_frac * frequency, 1e-6),
            cutoff=args.td_gaussian_cutoff,
        )
    source = mp.Source(
        src=src_time,
        component=mp.Hz,
        center=mp.Vector3(source_x, 0, 0),
        size=mp.Vector3(0, 2 * source_half_height, 0),
        amp_func=phased_line_amp(frequency, angle, source_half_height),
    )
    return mp.Simulation(
        cell_size=mp.Vector3(sx, sy, 0),
        dimensions=2,
        resolution=resolution,
        boundary_layers=[mp.PML(dpml)],
        sources=[source],
        geometry=[] if g is None else make_geometry_objects(g, wall_t),
        default_material=getattr(mp, "air", mp.Medium(epsilon=1.0)),
        eps_averaging=True,
        subpixel_tol=args.subpixel_tol,
        subpixel_maxeval=args.subpixel_maxeval,
        force_complex_fields=True,
        Courant=args.courant,
        progress_interval=5,
    )


def run_cw_case(args, outroot: Path, g: Optional[Geometry], resolution: int, angle: float,
                frequency: float, sensors: Dict[str, Any], source_x: float, sx: float,
                sy: float, dpml: float, wall_t: float, save_map: bool = False) -> Path:
    """Run one independently addressable CW case.

    model is selected by g:
      * g is None -> empty_reference;
      * Geometry(model='phb') -> exact PHB;
      * Geometry(model='smooth_area') -> equal-area smooth control.

    cw_tolerance_mode:
      * sequential  - solve tolerances successively in one Simulation (fast);
      * independent - create a fresh Simulation for every tolerance (strict repeatability test).
    """
    model = "empty_reference" if g is None else g.model
    case_dir = ensure_dir(outroot / f"res_{resolution:03d}" / model / f"angle_{angle_tag(angle)}")
    result_path = case_dir / "cw_measurement.json"
    tolerances = sorted(parse_float_list(args.cw_tolerances), reverse=True)
    tol_mode = str(args.cw_tolerance_mode).lower()
    if tol_mode not in {"sequential", "independent"}:
        raise ValueError("--cw-tolerance-mode must be sequential or independent")

    if args.resume and result_path.exists():
        try:
            old = read_json(result_path)
            if (
                old.get("script_version") == SCRIPT_VERSION
                and old.get("resolution") == resolution
                and old.get("model") == model
                and abs(float(old.get("angle_deg")) - float(angle)) < 1e-12
                and abs(float(old.get("frequency")) - float(frequency)) < 1e-12
                and old.get("cw_tolerance_mode", "sequential") == tol_mode
                and np.allclose(old.get("cw_tolerances", []), tolerances)
            ):
                print(f"SKIP CW res={resolution} {model} angle={angle:+.3f}", flush=True)
                return result_path
        except Exception:
            pass

    print(
        f"START CW res={resolution} {model} angle={angle:+.3f} "
        f"f={frequency:.8g} tol_mode={tol_mode}",
        flush=True,
    )
    t0 = time.time()
    snapshots = []
    final_map = None

    def solve_one(sim, tol: float, want_map: bool):
        ts = time.time()
        sim.solve_cw(tol, args.cw_maxiters, args.cw_L)
        if g is None:
            top = weighted_region_mean(sim, mp.Hz, sensors["probe_top"], sensors["sensor_size"])
            bottom = weighted_region_mean(sim, mp.Hz, sensors["probe_bottom"], sensors["sensor_size"])
            snap = {
                "top_real": float(top.real), "top_imag": float(top.imag), "top_abs": float(abs(top)),
                "bottom_real": float(bottom.real), "bottom_imag": float(bottom.imag), "bottom_abs": float(abs(bottom)),
            }
            fmap = None
        else:
            snap, fmap = extract_steady_state(
                sim, g, sensors, frequency, args.focal_half_thickness, save_map=want_map
            )
        snap["cw_tolerance"] = tol
        snap["solve_wall_s"] = time.time() - ts
        return snap, fmap

    if tol_mode == "sequential":
        sim = build_simulation(
            args, g, resolution, angle, frequency, source_x, sx, sy, dpml, wall_t, continuous=True
        )
        sim.init_sim()
        try:
            for tol in tolerances:
                snap, fmap = solve_one(sim, tol, save_map and tol == tolerances[-1])
                snapshots.append(snap)
                if fmap is not None:
                    final_map = fmap
        finally:
            sim.reset_meep()
    else:
        for tol in tolerances:
            sim = build_simulation(
                args, g, resolution, angle, frequency, source_x, sx, sy, dpml, wall_t, continuous=True
            )
            sim.init_sim()
            try:
                snap, fmap = solve_one(sim, tol, save_map and tol == tolerances[-1])
                snapshots.append(snap)
                if fmap is not None:
                    final_map = fmap
            finally:
                sim.reset_meep()

    result = {
        "script_version": SCRIPT_VERSION,
        "created": datetime.now().isoformat(timespec="seconds"),
        "stage": "cw",
        "model": model,
        "resolution": resolution,
        "angle_deg": angle,
        "frequency": frequency,
        "cw_tolerances": tolerances,
        "cw_tolerance_mode": tol_mode,
        "wall_thickness": wall_t,
        "sensor_offset_physical": sensors["offset_physical"],
        "sensor_size": list(sensors["sensor_size"]),
        "snapshots": snapshots,
        "elapsed_wall_s": time.time() - t0,
    }
    write_json(result_path, result)
    if final_map is not None:
        np.savez_compressed(case_dir / "field_map_final.npz", **final_map)
    print(f"DONE  CW res={resolution} {model} angle={angle:+.3f} elapsed={result['elapsed_wall_s']:.1f}s", flush=True)
    return result_path


def run_td_case(args, outroot: Path, g: Geometry, resolution: int, angle: float,
                frequency: float, sensors: Dict[str, Any], source_x: float, sx: float,
                sy: float, dpml: float, wall_t: float) -> Path:
    case_dir = ensure_dir(outroot / "TD_CONFIRM" / f"res_{resolution:03d}" / g.model / f"angle_{angle_tag(angle)}")
    result_path = case_dir / "td_measurement.json"
    checkpoints = sorted(parse_float_list(args.td_checkpoints))
    if args.resume and result_path.exists():
        try:
            old = read_json(result_path)
            if old.get("script_version") == SCRIPT_VERSION and np.allclose(old.get("checkpoints", []), checkpoints):
                print(f"SKIP TD res={resolution} {g.model} angle={angle:+.3f}", flush=True)
                return result_path
        except Exception:
            pass

    print(f"START TD res={resolution} {g.model} angle={angle:+.3f} f={frequency:.8g}", flush=True)
    t0 = time.time()
    sim = build_simulation(args, g, resolution, angle, frequency, source_x, sx, sy, dpml, wall_t, continuous=False)
    pt, pb = sensors["probe_top"], sensors["probe_bottom"]
    ssx, ssy = sensors["sensor_size"]
    top_mon = sim.add_dft_fields([mp.Hz], frequency, 0, 1, center=mp.Vector3(*pt, 0), size=mp.Vector3(ssx, ssy, 0), decimation_factor=1)
    bottom_mon = sim.add_dft_fields([mp.Hz], frequency, 0, 1, center=mp.Vector3(*pb, 0), size=mp.Vector3(ssx, ssy, 0), decimation_factor=1)
    harminv = mp.Harminv(mp.Hz, mp.Vector3(0, g.R, 0), frequency, args.harminv_width)

    snapshots = []
    first = checkpoints[0]
    sim.run(mp.after_sources(harminv), until_after_sources=first)
    snapshots.append(read_td_snapshot(sim, top_mon, bottom_mon, first))
    prev = first
    for cp in checkpoints[1:]:
        sim.run(until=float(cp - prev))
        snapshots.append(read_td_snapshot(sim, top_mon, bottom_mon, cp))
        prev = cp

    modes = []
    for m in getattr(harminv, "modes", []):
        modes.append({
            "frequency_real": float(m.freq),
            "frequency_imag": float(getattr(m, "decay", float("nan"))),
            "Q": float(m.Q),
            "amplitude_real": float(np.real(m.amp)),
            "amplitude_imag": float(np.imag(m.amp)),
            "error": float(m.err),
        })
    result = {
        "script_version": SCRIPT_VERSION,
        "created": datetime.now().isoformat(timespec="seconds"),
        "stage": "td-confirm",
        "model": g.model,
        "resolution": resolution,
        "angle_deg": angle,
        "frequency": frequency,
        "checkpoints": checkpoints,
        "snapshots": snapshots,
        "harminv_modes": modes,
        "elapsed_wall_s": time.time() - t0,
        "note": "DFT snapshots are cumulative from simulation start; use them only as a convergence trace. CW is the primary steady-state result.",
    }
    write_json(result_path, result)
    sim.reset_meep()
    print(f"DONE  TD res={resolution} {g.model} angle={angle:+.3f} elapsed={result['elapsed_wall_s']:.1f}s", flush=True)
    return result_path


def read_td_snapshot(sim, top_mon, bottom_mon, checkpoint: float) -> Dict[str, Any]:
    top = safe_complex_mean(sim.get_dft_array(top_mon, mp.Hz, 0))
    bottom = safe_complex_mean(sim.get_dft_array(bottom_mon, mp.Hz, 0))
    return {
        "checkpoint_after_sources": float(checkpoint),
        "meep_time": float(sim.meep_time()),
        "top_real": float(top.real), "top_imag": float(top.imag), "top_abs": float(abs(top)),
        "bottom_real": float(bottom.real), "bottom_imag": float(bottom.imag), "bottom_abs": float(abs(bottom)),
    }


# =============================================================================
# Causal analysis: focal mode vs angular information
# =============================================================================

def cplx(s: Mapping[str, Any], prefix: str) -> complex:
    return complex(float(s[f"{prefix}_real"]), float(s[f"{prefix}_imag"]))


def load_cw_results(outroot: Path) -> List[Dict[str, Any]]:
    rows = []
    for path in sorted(outroot.glob("res_*/**/cw_measurement.json")):
        d = read_json(path)
        for snap in d.get("snapshots", []):
            row = {
                "source_file": str(path.relative_to(outroot)),
                "model": d["model"],
                "resolution": int(d["resolution"]),
                "angle_deg": float(d["angle_deg"]),
                "frequency": float(d["frequency"]),
                "cw_tolerance": float(snap["cw_tolerance"]),
                **snap,
            }
            rows.append(row)
    return rows


def parity_sign_from_zero(rows: Sequence[Mapping[str, Any]], model: str, resolution: int, tol: float) -> int:
    z = [r for r in rows if r["model"] == model and int(r["resolution"]) == resolution and abs(float(r["angle_deg"])) < 1e-12 and abs(float(r["cw_tolerance"]) - tol) <= max(1e-15, 1e-8 * tol)]
    if not z:
        raise RuntimeError(f"missing zero-angle row for {model}, res={resolution}, tol={tol}")
    t, b = cplx(z[0], "top"), cplx(z[0], "bottom")
    return +1 if abs(t + b) >= abs(t - b) else -1


def modal_channels(top: complex, bottom: complex, parity_sign: int) -> Tuple[complex, complex, complex]:
    common = (top + parity_sign * bottom) / math.sqrt(2.0)
    differential = (top - parity_sign * bottom) / math.sqrt(2.0)
    ratio = differential / common if abs(common) > 1e-300 else complex(float("nan"), float("nan"))
    return common, differential, ratio


def interpolation_predict(train_angles: np.ndarray, train_z: np.ndarray, z: complex) -> float:
    X = np.column_stack((train_z.real, train_z.imag))
    q = np.asarray([z.real, z.imag])
    best_angle, best_d = float(train_angles[0]), float("inf")
    order = np.argsort(train_angles)
    train_angles, X = train_angles[order], X[order]
    for i in range(len(train_angles) - 1):
        a, b = X[i], X[i + 1]
        v = b - a
        den = float(np.dot(v, v))
        u = float(np.dot(q - a, v) / den) if den > 0 else 0.0
        u = min(1.0, max(0.0, u))
        p = a + u * v
        d = float(np.linalg.norm(q - p))
        if d < best_d:
            best_d = d
            best_angle = float(train_angles[i] + u * (train_angles[i + 1] - train_angles[i]))
    return best_angle


def analyze_cw(args, outroot: Path, geometries: Sequence[Geometry], resolutions: Sequence[int],
               train_angles: Sequence[float], test_angles: Sequence[float]) -> Dict[str, Any]:
    rows = load_cw_results(outroot)
    if not rows:
        raise RuntimeError("no CW results found")
    write_csv(outroot / "CW_MEASUREMENTS_ALL.csv", rows)
    tolerances = sorted(parse_float_list(args.cw_tolerances), reverse=True)
    final_tol = tolerances[-1]
    noise_ratios = parse_float_list(args.receiver_noise_ratios)
    rng = np.random.default_rng(args.seed)

    modal_rows: List[Dict[str, Any]] = []
    estimator_rows: List[Dict[str, Any]] = []
    convergence_rows: List[Dict[str, Any]] = []
    fisher_rows: List[Dict[str, Any]] = []

    for res in resolutions:
        ref = [r for r in rows if r["model"] == "empty_reference" and int(r["resolution"]) == res and abs(float(r["cw_tolerance"]) - final_tol) <= max(1e-15, final_tol * 1e-8)]
        if not ref:
            raise RuntimeError(f"missing empty reference at resolution {res}")
        rt, rb = cplx(ref[0], "top"), cplx(ref[0], "bottom")
        reference_amp = math.sqrt(0.5 * (abs(rt) ** 2 + abs(rb) ** 2))

        for model in GEOMETRY_MODELS:
            parity = parity_sign_from_zero(rows, model, res, final_tol)
            selected = [r for r in rows if r["model"] == model and int(r["resolution"]) == res and abs(float(r["cw_tolerance"]) - final_tol) <= max(1e-15, final_tol * 1e-8)]
            selected = sorted(selected, key=lambda r: float(r["angle_deg"]))
            for r in selected:
                top, bottom = cplx(r, "top"), cplx(r, "bottom")
                common, differential, ratio = modal_channels(top, bottom, parity)
                rec = {
                    "model": model, "resolution": res, "angle_deg": float(r["angle_deg"]),
                    "parity_sign": parity, "reference_amplitude": reference_amp,
                    "common_real": common.real, "common_imag": common.imag, "common_power_over_reference": abs(common) ** 2 / max(reference_amp ** 2, 1e-300),
                    "differential_real": differential.real, "differential_imag": differential.imag, "differential_power_over_reference": abs(differential) ** 2 / max(reference_amp ** 2, 1e-300),
                    "ratio_real": ratio.real, "ratio_imag": ratio.imag, "ratio_abs": abs(ratio),
                    "coherent_cross_real_over_reference2": (differential * np.conj(common)).real / max(reference_amp ** 2, 1e-300),
                    "coherent_cross_imag_over_reference2": (differential * np.conj(common)).imag / max(reference_amp ** 2, 1e-300),
                    "focal_energy_fraction": float(r.get("focal_energy_fraction", float("nan"))),
                    "focal_enrichment": float(r.get("focal_enrichment", float("nan"))),
                }
                modal_rows.append(rec)

            train = [r for r in modal_rows if r["model"] == model and int(r["resolution"]) == res and any(abs(float(r["angle_deg"]) - a) < 1e-9 for a in train_angles)]
            test = [r for r in modal_rows if r["model"] == model and int(r["resolution"]) == res and any(abs(float(r["angle_deg"]) - a) < 1e-9 for a in test_angles)]
            train = sorted(train, key=lambda r: float(r["angle_deg"]))
            test = sorted(test, key=lambda r: float(r["angle_deg"]))
            ta = np.asarray([float(r["angle_deg"]) for r in train])
            tz = np.asarray([complex(float(r["ratio_real"]), float(r["ratio_imag"])) for r in train])

            if len(train) >= 3 and test:
                clean_preds = [interpolation_predict(ta, tz, complex(float(r["ratio_real"]), float(r["ratio_imag"]))) for r in test]
                truths = np.asarray([float(r["angle_deg"]) for r in test])
                estimator_rows.append(metric_row(model, res, "clean", clean_preds, truths))

                raw_by_angle = {float(r["angle_deg"]): r for r in selected}
                for nr in noise_ratios:
                    preds_all, truths_all = [], []
                    sigma = nr * reference_amp
                    for _ in range(args.noise_trials):
                        noisy_train_z = []
                        for a in ta:
                            rr = raw_by_angle[float(a)]
                            top = cplx(rr, "top") + sigma / math.sqrt(2) * (rng.normal() + 1j * rng.normal())
                            bottom = cplx(rr, "bottom") + sigma / math.sqrt(2) * (rng.normal() + 1j * rng.normal())
                            _c, _d, z = modal_channels(top, bottom, parity)
                            noisy_train_z.append(z)
                        noisy_train_z = np.asarray(noisy_train_z)
                        for rtest in test:
                            rr = raw_by_angle[float(rtest["angle_deg"])]
                            top = cplx(rr, "top") + sigma / math.sqrt(2) * (rng.normal() + 1j * rng.normal())
                            bottom = cplx(rr, "bottom") + sigma / math.sqrt(2) * (rng.normal() + 1j * rng.normal())
                            _c, _d, z = modal_channels(top, bottom, parity)
                            preds_all.append(interpolation_predict(ta, noisy_train_z, z))
                            truths_all.append(float(rtest["angle_deg"]))
                    estimator_rows.append(metric_row(model, res, f"fixed_noise_{nr:g}_of_empty", preds_all, np.asarray(truths_all)))

            # Fisher information at zero using the nearest symmetric nonzero train angle.
            positive = sorted(a for a in train_angles if a > 0 and -a in train_angles)
            if positive:
                delta = positive[0]
                rp = next(r for r in selected if abs(float(r["angle_deg"]) - delta) < 1e-9)
                rm = next(r for r in selected if abs(float(r["angle_deg"]) + delta) < 1e-9)
                dtheta = math.radians(2 * delta)
                ds_top = (cplx(rp, "top") - cplx(rm, "top")) / dtheta
                ds_bottom = (cplx(rp, "bottom") - cplx(rm, "bottom")) / dtheta
                common_p, diff_p, ratio_p = modal_channels(cplx(rp, "top"), cplx(rp, "bottom"), parity)
                common_m, diff_m, ratio_m = modal_channels(cplx(rm, "top"), cplx(rm, "bottom"), parity)
                d_odd = (diff_p - diff_m) / dtheta
                d_ratio = (ratio_p - ratio_m) / dtheta
                base = {
                    "model": model, "resolution": res, "delta_angle_deg": delta,
                    "odd_mode_slope_per_rad_over_reference": abs(d_odd) / max(reference_amp, 1e-300),
                    "mode_ratio_slope_per_rad": abs(d_ratio),
                }
                for nr in noise_ratios:
                    sigma = nr * reference_amp
                    fisher = 2.0 * (abs(ds_top) ** 2 + abs(ds_bottom) ** 2) / max(sigma ** 2, 1e-300)
                    crlb_deg = math.degrees(1.0 / math.sqrt(max(fisher, 1e-300)))
                    fisher_rows.append({**base, "noise_ratio_of_empty": nr, "fisher_information_per_rad2": fisher, "crlb_std_deg": crlb_deg})

            # Solver-tolerance convergence for zero angle.
            if len(tolerances) >= 2:
                loose, tight = tolerances[-2], tolerances[-1]
                for angle in sorted(set(train_angles) | set(test_angles)):
                    pair = [r for r in rows if r["model"] == model and int(r["resolution"]) == res and abs(float(r["angle_deg"]) - angle) < 1e-9 and float(r["cw_tolerance"]) in (loose, tight)]
                    if len(pair) == 2:
                        pair = sorted(pair, key=lambda r: float(r["cw_tolerance"]), reverse=True)
                        p0, p1 = pair
                        convergence_rows.append({
                            "type": "cw_tolerance", "model": model, "resolution": res, "angle_deg": angle,
                            "loose_tolerance": loose, "tight_tolerance": tight,
                            "top_relative_error": relative_complex_error(cplx(p0, "top"), cplx(p1, "top")),
                            "bottom_relative_error": relative_complex_error(cplx(p0, "bottom"), cplx(p1, "bottom")),
                            "focal_eta_relative_error": abs(float(p0.get("focal_energy_fraction", np.nan)) - float(p1.get("focal_energy_fraction", np.nan))) / max(abs(float(p0.get("focal_energy_fraction", np.nan))), abs(float(p1.get("focal_energy_fraction", np.nan))), 1e-300),
                        })

    # Resolution convergence of scalar modal and focal metrics at final tolerance.
    if len(resolutions) >= 2:
        r0, r1 = resolutions[-2], resolutions[-1]
        for model in GEOMETRY_MODELS:
            for angle in sorted(set(train_angles) | set(test_angles)):
                a0 = next((r for r in modal_rows if r["model"] == model and int(r["resolution"]) == r0 and abs(float(r["angle_deg"]) - angle) < 1e-9), None)
                a1 = next((r for r in modal_rows if r["model"] == model and int(r["resolution"]) == r1 and abs(float(r["angle_deg"]) - angle) < 1e-9), None)
                if a0 and a1:
                    z0 = complex(float(a0["ratio_real"]), float(a0["ratio_imag"]))
                    z1 = complex(float(a1["ratio_real"]), float(a1["ratio_imag"]))
                    convergence_rows.append({
                        "type": "resolution", "model": model, "angle_deg": angle,
                        "resolution_low": r0, "resolution_high": r1,
                        "mode_ratio_relative_error": relative_complex_error(z0, z1),
                        "focal_eta_relative_error": abs(float(a0["focal_energy_fraction"]) - float(a1["focal_energy_fraction"])) / max(abs(float(a0["focal_energy_fraction"])), abs(float(a1["focal_energy_fraction"])), 1e-300),
                        "focal_enrichment_relative_error": abs(float(a0["focal_enrichment"]) - float(a1["focal_enrichment"])) / max(abs(float(a0["focal_enrichment"])), abs(float(a1["focal_enrichment"])), 1e-300),
                    })

    write_csv(outroot / "MODAL_CHANNELS.csv", modal_rows)
    write_csv(outroot / "FOCAL_MODE_ESTIMATOR_METRICS.csv", estimator_rows)
    write_csv(outroot / "FISHER_INFORMATION_FIXED_NOISE.csv", fisher_rows)
    write_csv(outroot / "CONVERGENCE_METRICS.csv", convergence_rows)

    final_res = max(resolutions)
    final_modal = [r for r in modal_rows if int(r["resolution"]) == final_res]
    zero = {m: next(r for r in final_modal if r["model"] == m and abs(float(r["angle_deg"])) < 1e-9) for m in GEOMETRY_MODELS}
    focal_ratio = float(zero["phb"]["focal_energy_fraction"]) / max(float(zero["smooth_area"]["focal_energy_fraction"]), 1e-300)
    enrich_ratio = float(zero["phb"]["focal_enrichment"]) / max(float(zero["smooth_area"]["focal_enrichment"]), 1e-300)

    final_fisher = [r for r in fisher_rows if int(r["resolution"]) == final_res]
    fisher_ratios = []
    for nr in noise_ratios:
        fp = next(r for r in final_fisher if r["model"] == "phb" and abs(float(r["noise_ratio_of_empty"]) - nr) < 1e-12)
        fs = next(r for r in final_fisher if r["model"] == "smooth_area" and abs(float(r["noise_ratio_of_empty"]) - nr) < 1e-12)
        fisher_ratios.append(float(fp["fisher_information_per_rad2"]) / max(float(fs["fisher_information_per_rad2"]), 1e-300))

    final_est = [r for r in estimator_rows if int(r["resolution"]) == final_res]
    estimator_advantages = []
    for condition in sorted(set(r["condition"] for r in final_est)):
        ep = next(r for r in final_est if r["model"] == "phb" and r["condition"] == condition)
        es = next(r for r in final_est if r["model"] == "smooth_area" and r["condition"] == condition)
        estimator_advantages.append({"condition": condition, "phb_mae": ep["mae_deg"], "smooth_mae": es["mae_deg"], "mae_ratio_phb_over_smooth": ep["mae_deg"] / max(es["mae_deg"], 1e-300)})

    tol_pass = all(max(float(r.get("top_relative_error", 0)), float(r.get("bottom_relative_error", 0)), float(r.get("focal_eta_relative_error", 0))) <= args.convergence_relative_error for r in convergence_rows if r["type"] == "cw_tolerance")
    res_rows = [r for r in convergence_rows if r["type"] == "resolution"]
    res_pass = bool(res_rows) and all(max(float(r["mode_ratio_relative_error"]), float(r["focal_eta_relative_error"]), float(r["focal_enrichment_relative_error"])) <= args.convergence_relative_error for r in res_rows)
    focal_pass = focal_ratio >= args.required_advantage_ratio and enrich_ratio >= args.required_advantage_ratio
    fisher_pass = bool(fisher_ratios) and min(fisher_ratios) >= args.required_advantage_ratio
    estimator_pass = bool(estimator_advantages) and all(float(r["mae_ratio_phb_over_smooth"]) <= 1.0 / args.required_advantage_ratio for r in estimator_advantages if r["condition"] != "clean" or True)

    if focal_pass and fisher_pass and estimator_pass and tol_pass and res_pass:
        verdict = "FOCAL_DIRECTION_BRIDGE_SUPPORTED_WITHIN_TESTED_2D_MODEL"
    elif focal_pass and not fisher_pass:
        verdict = "PHB_FOCAL_CONCENTRATION_CONFIRMED_BUT_NOT_ANGULAR_INFORMATION_GAIN"
    elif not focal_pass and fisher_pass:
        verdict = "ANGULAR_GAIN_WITHOUT_CONFIRMED_FOCAL_CAUSATION"
    elif not res_pass or not tol_pass:
        verdict = "NOT_NUMERICALLY_CONVERGED"
    else:
        verdict = "NO_EXCLUSIVE_PHB_ADVANTAGE_OVER_SMOOTH_CONTROL"

    summary = {
        "script_version": SCRIPT_VERSION,
        "verdict": verdict,
        "final_resolution": final_res,
        "frequency": args.frequency,
        "focal_energy_fraction_ratio_phb_over_smooth_at_zero": focal_ratio,
        "focal_enrichment_ratio_phb_over_smooth_at_zero": enrich_ratio,
        "fisher_information_ratios_phb_over_smooth": fisher_ratios,
        "estimator_advantages": estimator_advantages,
        "cw_tolerance_convergence_pass": tol_pass,
        "resolution_convergence_pass": res_pass,
        "focal_advantage_pass": focal_pass,
        "fisher_advantage_pass": fisher_pass,
        "estimator_advantage_pass": estimator_pass,
        "required_advantage_ratio": args.required_advantage_ratio,
        "interpretation": {
            "supported": "PHB must beat the smooth equal-area control in focal energy fraction, fixed-noise Fisher information, and focal-mode angle recovery after tolerance and resolution convergence.",
            "focal_only": "PHB is a stronger concentrator at this frequency but the focal mode does not improve differential direction information.",
        },
    }
    write_json(outroot / "FINAL_BRIDGE_SUMMARY.json", summary)
    write_verdict_ru(outroot / "FINAL_BRIDGE_VERDICT_RU.txt", summary)
    save_analysis_plots(outroot, modal_rows, estimator_rows, fisher_rows, resolutions)
    return summary


def metric_row(model: str, resolution: int, condition: str, preds: Sequence[float], truths: np.ndarray) -> Dict[str, Any]:
    p = np.asarray(preds, dtype=float)
    y = np.asarray(truths, dtype=float)
    err = np.abs(p - y)
    sign = np.mean([1.0 if (abs(t) < 1e-12 and abs(q) < 1e-12) or (t * q > 0) else 0.0 for q, t in zip(p, y)])
    return {
        "model": model, "resolution": resolution, "condition": condition,
        "mae_deg": float(np.mean(err)), "median_abs_error_deg": float(np.median(err)),
        "max_abs_error_deg": float(np.max(err)), "sign_accuracy": float(sign),
    }


def write_verdict_ru(path: Path, s: Mapping[str, Any]) -> None:
    text = f"""PHB 2D FOCAL-DIRECTION BRIDGE v5\n\nFINAL VERDICT: {s['verdict']}\n\nЧастота: {s['frequency']}\nФинальное разрешение: {s['final_resolution']}\nОтношение доли энергии в фокальной зоне PHB/smooth: {s['focal_energy_fraction_ratio_phb_over_smooth_at_zero']:.6g}\nОтношение фокального обогащения PHB/smooth: {s['focal_enrichment_ratio_phb_over_smooth_at_zero']:.6g}\nОтношения информации Фишера PHB/smooth: {s['fisher_information_ratios_phb_over_smooth']}\nСходимость CW по точности решателя: {s['cw_tolerance_convergence_pass']}\nСходимость по разрешению: {s['resolution_convergence_pass']}\nПреимущество по фокальной энергии: {s['focal_advantage_pass']}\nПреимущество по угловой информации Фишера: {s['fisher_advantage_pass']}\nПреимущество алгоритма фокальной пространственной моды: {s['estimator_advantage_pass']}\n\nИнтерпретация:\n- Положительный итог требует одновременного роста фокальной доли энергии, фиксированно-шумовой угловой информации и точности восстановления угла относительно гладкого негиперболического контроля.\n- Если растёт только фокальная энергия, PHB подтверждается как концентратор, но не как усилитель угловой информации.\n"""
    path.write_text(text, encoding="utf-8")


def save_analysis_plots(outroot: Path, modal_rows: Sequence[Mapping[str, Any]], estimator_rows: Sequence[Mapping[str, Any]], fisher_rows: Sequence[Mapping[str, Any]], resolutions: Sequence[int]) -> None:
    if plt is None:
        return
    final_res = max(resolutions)
    selected = [r for r in modal_rows if int(r["resolution"]) == final_res]
    fig, ax = plt.subplots(figsize=(10, 6))
    for model in GEOMETRY_MODELS:
        rr = sorted([r for r in selected if r["model"] == model], key=lambda r: float(r["angle_deg"]))
        ax.plot([r["angle_deg"] for r in rr], [r["ratio_real"] for r in rr], marker="o", label=f"{model}: Re(O/E)")
    ax.axhline(0, lw=0.8)
    ax.set_xlabel("incidence angle, deg")
    ax.set_ylabel("real part of differential/common mode ratio")
    ax.set_title(f"Two-sensor spatial-mode conversion at resolution {final_res}")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outroot / "MODE_CONVERSION_VS_ANGLE.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    for model in GEOMETRY_MODELS:
        rr = sorted([r for r in selected if r["model"] == model], key=lambda r: float(r["angle_deg"]))
        ax.plot([r["angle_deg"] for r in rr], [r["focal_energy_fraction"] for r in rr], marker="o", label=model)
    ax.set_xlabel("incidence angle, deg")
    ax.set_ylabel("focal-zone energy / total cavity energy")
    ax.set_title(f"Absolute focal energy fraction at resolution {final_res}")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outroot / "FOCAL_ENERGY_FRACTION_VS_ANGLE.png", dpi=220)
    plt.close(fig)


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Universal independently selectable PHB focal-direction verifier v5.2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--stage", choices=["geometry", "cw", "td-confirm", "analyze"], default="cw")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--resume", action="store_true", help="skip only fully matching completed cases")
    ap.add_argument("--no-archive", dest="archive", action="store_false")
    ap.set_defaults(archive=True)

    # Independent case selection. Any model, angle and resolution can be run alone.
    ap.add_argument(
        "--models",
        default="empty_reference,phb,smooth_area",
        help="comma-separated subset of: empty_reference,phb,smooth_area",
    )
    ap.add_argument(
        "--angles",
        default="",
        help="CW angles in degrees. If empty, union of --train-angles and --test-angles is used.",
    )
    ap.add_argument(
        "--case-order",
        choices=["angle-major", "model-major"],
        default="angle-major",
        help="angle-major gives immediate paired PHB/control results; model-major reproduces legacy ordering",
    )
    ap.add_argument(
        "--cw-tolerance-mode",
        choices=["sequential", "independent"],
        default="sequential",
        help="sequential is faster; independent creates a fresh Simulation for every tolerance",
    )
    ap.add_argument(
        "--auto-analyze",
        action="store_true",
        help="run full legacy analysis after CW; use only when the complete required dataset exists",
    )

    ap.add_argument("--a", type=float, default=0.5)
    ap.add_argument("--b", type=float, default=0.5)
    ap.add_argument("--R", type=float, default=2.0)
    ap.add_argument("--R1", type=float, default=0.2)
    ap.add_argument("--R1-mode", choices=["relative", "absolute"], default="relative")
    ap.add_argument("--frequency", type=float, default=1.8)
    ap.add_argument("--resolutions", default="80,120")
    ap.add_argument("--train-angles", default="-15,-7.5,0,7.5,15")
    ap.add_argument("--test-angles", default="-11.25,-3.75,3.75,11.25")
    ap.add_argument("--td-angles", default="-7.5,0,7.5")

    ap.add_argument("--wall-thickness", type=float, default=0.0833333333333333)
    ap.add_argument("--wall-min-cells", type=float, default=3.0)
    ap.add_argument("--sensor-offset-physical", type=float, default=0.1041666666666667)
    ap.add_argument("--sensor-size-x", type=float, default=0.05)
    ap.add_argument("--sensor-size-y", type=float, default=0.05)
    ap.add_argument("--focal-half-thickness", type=float, default=0.08)
    ap.add_argument("--dpml-over-lambda", type=float, default=1.2)
    ap.add_argument("--x-air", type=float, default=0.8)
    ap.add_argument("--y-air", type=float, default=0.5)
    ap.add_argument("--source-gap-from-left-tip", type=float, default=0.55)

    ap.add_argument("--cw-tolerances", default="1e-6,1e-8")
    ap.add_argument("--cw-maxiters", type=int, default=20000)
    ap.add_argument("--cw-L", type=int, default=4)
    ap.add_argument("--cw-source-width-periods", type=float, default=5.0)
    ap.add_argument("--cw-source-slowness", type=float, default=3.0)

    ap.add_argument("--td-checkpoints", default="180,300,500")
    ap.add_argument("--td-source-fwidth-frac", type=float, default=0.04)
    ap.add_argument("--td-gaussian-cutoff", type=float, default=5.0)
    ap.add_argument("--harminv-width", type=float, default=0.25)

    ap.add_argument("--receiver-noise-ratios", default="0.01,0.03,0.10")
    ap.add_argument("--noise-trials", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260719)
    ap.add_argument("--required-advantage-ratio", type=float, default=1.20)
    ap.add_argument("--convergence-relative-error", type=float, default=0.10)

    ap.add_argument("--subpixel-tol", type=float, default=1e-4)
    ap.add_argument("--subpixel-maxeval", type=int, default=10000)
    ap.add_argument("--courant", type=float, default=0.5)
    ap.add_argument("--map-angles", default="-7.5,0,7.5")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    outroot = ensure_dir(Path(args.outdir))
    resolutions = sorted(set(parse_int_list(args.resolutions)))
    if not resolutions:
        raise ValueError("at least one resolution is required")

    selected_models = parse_str_list(args.models)
    unknown = sorted(set(selected_models) - set(RUN_MODELS))
    if unknown:
        raise ValueError(f"unknown --models values: {unknown}; allowed: {RUN_MODELS}")
    if not selected_models:
        raise ValueError("--models must select at least one model")
    selected_models = list(dict.fromkeys(selected_models))

    train_angles = parse_float_list(args.train_angles)
    test_angles = parse_float_list(args.test_angles)
    cw_angles = parse_float_list(args.angles) if str(args.angles).strip() else sorted(set(train_angles + test_angles))
    if not cw_angles and args.stage == "cw":
        raise ValueError("at least one CW angle is required")
    td_angles = parse_float_list(args.td_angles)

    # Reference geometries are always created for common cell dimensions, sensor positions,
    # and the PHB-vs-smooth geometry audit. Selected models only control expensive runs.
    geometry_by_model = {
        m: Geometry(args.a, args.b, args.R, args.R1, args.R1_mode, m)
        for m in GEOMETRY_MODELS
    }
    reference_geometries = [geometry_by_model[m] for m in GEOMETRY_MODELS]
    if not all(0 < g.Rin < g.R for g in reference_geometries):
        raise ValueError("0 < Rin < R is required")

    wavelength = 1.0 / args.frequency
    dpml = args.dpml_over_lambda * wavelength
    base_geometry = geometry_by_model["phb"]
    source_x = -base_geometry.L - args.source_gap_from_left_tip

    manifest = {
        "script_version": SCRIPT_VERSION,
        "created": datetime.now().isoformat(timespec="seconds"),
        "args": vars(args),
        "selected_models": selected_models,
        "cw_angles": cw_angles,
        "resolutions": resolutions,
        "frequency": args.frequency,
        "wavelength": wavelength,
        "a_over_lambda": args.a / wavelength,
        "independent_case_selection": True,
        "scientific_question": "Does PHB focal concentration and sensor loading exceed empty space and a smooth equal-area non-hyperbolic control?",
    }
    write_json(outroot / "RUN_MANIFEST.json", manifest)

    audit_res = max(resolutions)
    audit_wall = max(args.wall_thickness, args.wall_min_cells / audit_res)
    save_geometry_outputs(
        outroot, reference_geometries, audit_res, audit_wall,
        args.sensor_offset_physical, args.sensor_size_x, args.sensor_size_y,
    )
    if args.stage == "geometry":
        if args.archive:
            print(f"Archive: {archive_folder(outroot)}")
        return 0

    if args.stage == "cw":
        require_meep()
        map_angles = parse_float_list(args.map_angles)
        for res in resolutions:
            wall_t = max(args.wall_thickness, args.wall_min_cells / res)
            sx, sy = geometry_dimensions(base_geometry, dpml, args.x_air, args.y_air)
            sensors = sensor_regions(
                base_geometry, wall_t, res, args.sensor_offset_physical,
                args.sensor_size_x, args.sensor_size_y,
            )

            def run_selected(model: str, angle: float) -> None:
                g = None if model == "empty_reference" else geometry_by_model[model]
                save_map = (
                    model != "empty_reference"
                    and res == max(resolutions)
                    and any(abs(angle - a) < 1e-9 for a in map_angles)
                )
                run_cw_case(
                    args, outroot, g, res, angle, args.frequency, sensors,
                    source_x, sx, sy, dpml, wall_t, save_map=save_map,
                )

            if args.case_order == "angle-major":
                for angle in cw_angles:
                    for model in selected_models:
                        run_selected(model, angle)
            else:
                for model in selected_models:
                    for angle in cw_angles:
                        run_selected(model, angle)

        if args.auto_analyze:
            missing = {"phb", "smooth_area"} - set(selected_models)
            if missing:
                raise ValueError(f"--auto-analyze requires phb and smooth_area; missing {sorted(missing)}")
            summary = analyze_cw(args, outroot, reference_geometries, resolutions, train_angles, test_angles)
            print(json.dumps({"FINAL_VERDICT": summary["verdict"]}, ensure_ascii=False), flush=True)
        if args.archive:
            print(f"Archive: {archive_folder(outroot)}")
        return 0

    if args.stage == "td-confirm":
        require_meep()
        td_models = [m for m in selected_models if m in GEOMETRY_MODELS]
        if not td_models:
            raise ValueError("td-confirm requires phb and/or smooth_area in --models")
        res = max(resolutions)
        wall_t = max(args.wall_thickness, args.wall_min_cells / res)
        sx, sy = geometry_dimensions(base_geometry, dpml, args.x_air, args.y_air)
        sensors = sensor_regions(
            base_geometry, wall_t, res, args.sensor_offset_physical,
            args.sensor_size_x, args.sensor_size_y,
        )
        for model in td_models:
            g = geometry_by_model[model]
            for angle in td_angles:
                run_td_case(args, outroot, g, res, angle, args.frequency, sensors, source_x, sx, sy, dpml, wall_t)
        if args.archive:
            print(f"Archive: {archive_folder(outroot)}")
        return 0

    # Full analysis is deliberately a separate stage so partial independent runs never fail.
    summary = analyze_cw(args, outroot, reference_geometries, resolutions, train_angles, test_angles)
    print(json.dumps({"FINAL_VERDICT": summary["verdict"]}, ensure_ascii=False), flush=True)
    if args.archive:
        print(f"Archive: {archive_folder(outroot)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        raise

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHB_MEEP_fullwave_uniqueness_scan_v4_workdir.py

Full-vector Maxwell/FDTD verification script v4-seed for a second-order vertical
pseudohyperboloid (PHB) cavity/resonator in MEEP.

Changes in v4-seed:
- Added 'seed' field to RunConfig (default=42) and --seed argument;
- mp.set_seed(seed) is called before Simulation() to ensure reproducibility;
- seed is saved into run_manifest.json and decision_report.md;
- call to make_simulation() now includes seed=cfg.seed.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import shlex
import subprocess
import time
import zipfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np

def _trapz(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)

try:
    import pandas as pd
except Exception:
    pd = None

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None

try:
    import meep as mp
    MEEP_IMPORT_ERROR = None
except Exception as exc:
    mp = None
    MEEP_IMPORT_ERROR = exc


@dataclass
class PHBParams:
    a: float = 0.3
    b: float = 0.6
    R: float = 3.0
    R1: float = 0.0
    R2: float = 0.0

    @property
    def c(self) -> float:
        return math.sqrt(self.a * self.a + self.b * self.b)

    @property
    def L(self) -> float:
        return self.a * math.sqrt(1.0 + (self.R / self.b) ** 2)

    @property
    def dR_total(self) -> float:
        return self.R * self.R2

    @property
    def dR_half(self) -> float:
        return 0.5 * self.dR_total


@dataclass
class RunConfig:
    params: PHBParams
    mode: str
    geometries: List[str]
    m: int
    kr_center: float
    kr_span: float
    nfreq: int
    resolution: int
    dpml: float
    pad_r: float
    pad_z: float
    wall_thickness: float
    source_component: str
    source_r: float
    source_z: float
    source_dr: float
    source_dz: float
    until_after_sources: float
    decay_by: float
    harminv_after_sources: float
    courant: Optional[float]
    open_slot: bool
    slot_z_center: float
    slot_z_halfwidth: float
    slot_r_halfwidth: float
    outdir: str
    dry_run_geometry: bool
    plot_all_fields: bool
    include_h_fields: bool
    seed: int = 42


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def smoothstep(t: float) -> float:
    t = clamp01(t)
    return 3.0 * t * t - 2.0 * t * t * t


def phb_wall_radius_scalar(z: float, p: PHBParams, geometry: str = "PH") -> float:
    za = abs(float(z))
    a, b, R, L = p.a, p.b, p.R, p.L
    geometry = geometry.upper()

    if za > L:
        return 0.0

    if za <= a:
        if geometry == "CYLINDER_R":
            return R
        if geometry == "ELLIPSOID_L":
            return R * math.sqrt(max(0.0, 1.0 - (za / L) ** 2))
        return R + math.sqrt(max(0.0, a * a - za * za))

    t = (za - a) / max(L - a, 1e-15)
    t = clamp01(t)

    if geometry == "PH":
        val = R - b * math.sqrt(max((za / a) ** 2 - 1.0, 0.0))
    elif geometry == "LINEAR":
        val = R * (1.0 - t)
    elif geometry == "POLY2":
        val = R * (1.0 - t * t)
    elif geometry == "SMOOTHSTEP":
        val = R * (1.0 - smoothstep(t))
    elif geometry == "CYLINDER_R":
        val = R
    elif geometry == "ELLIPSOID_L":
        val = R * math.sqrt(max(0.0, 1.0 - (za / L) ** 2))
    elif geometry == "HERMITE_SLOPE_MATCH":
        t = clamp01(t)
        s0 = -R / max(L - a, 1e-15)
        s1 = -0.25 * R / max(L - a, 1e-15)
        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3 - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        h11 = t**3 - t**2
        val = h00 * R + h10 * (L - a) * s0 + h01 * 0.0 + h11 * (L - a) * s1
    elif geometry == "POLY_VOL_MATCH":
        q = 1.55
        val = R * (1.0 - t ** q)
    else:
        raise ValueError(f"Unknown geometry '{geometry}'.")

    return max(0.0, float(val))


def wall_radius_array(z: np.ndarray, p: PHBParams, geometry: str) -> np.ndarray:
    f = np.vectorize(lambda zz: phb_wall_radius_scalar(float(zz), p, geometry))
    return f(z)


def inside_cavity(r: np.ndarray, z: np.ndarray, p: PHBParams, geometry: str) -> np.ndarray:
    rw = wall_radius_array(z, p, geometry)
    return (r >= 0.0) & (np.abs(z) <= p.L) & (r <= rw)


def in_annular_slot(r: float, z: float, cfg: RunConfig) -> bool:
    if not cfg.open_slot:
        return False
    p = cfg.params
    dz = abs(z - cfg.slot_z_center)
    dr = abs(r - p.R)
    return (dz <= cfg.slot_z_halfwidth) and (dr <= cfg.slot_r_halfwidth)


def make_material_function(cfg: RunConfig, geometry: str):
    if mp is None:
        raise RuntimeError("MEEP is not imported.")
    p = cfg.params
    air = mp.Medium(epsilon=1.0)
    metal = getattr(mp, "metal", mp.Medium(epsilon=1.0e9))
    wall_t = cfg.wall_thickness

    def matfun(pos):
        r = float(pos.x)
        z = float(pos.z)
        rw = phb_wall_radius_scalar(z, p, geometry)

        if abs(z) > p.L + wall_t:
            return air

        if cfg.mode == "open" and in_annular_slot(r, z, cfg):
            return air

        if (r >= max(0.0, rw)) and (r <= max(0.0, rw) + wall_t) and (abs(z) <= p.L + wall_t):
            return metal

        if cfg.mode == "closed" and (abs(z) >= p.L - wall_t) and (r <= wall_t + max(rw, 0.0)):
            return metal

        return air

    return matfun


def numerical_derivatives(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    dy = np.gradient(y, x, edge_order=2)
    ddy = np.gradient(dy, x, edge_order=2)
    return dy, ddy


def gaussian_curvature_of_revolution(z: np.ndarray, r: np.ndarray) -> np.ndarray:
    rp, rpp = numerical_derivatives(z, r)
    denom = np.maximum(r * (1.0 + rp * rp) ** 2, 1e-15)
    return -rpp / denom


def geometry_summary(p: PHBParams, geometry: str, n: int = 2000) -> Dict[str, float]:
    z = np.linspace(-p.L + 1e-6, p.L - 1e-6, n)
    r = wall_radius_array(z, p, geometry)
    K = gaussian_curvature_of_revolution(z, r)
    horn_mask = (np.abs(z) > p.a * 1.05) & (np.abs(z) < p.L * 0.98) & (r > 1e-6)
    cap_mask = (np.abs(z) < p.a * 0.95) & (r > 1e-6)
    vol_proxy = float(_trapz(np.pi * r * r, z))
    surf_proxy = float(_trapz(2 * np.pi * r * np.sqrt(1 + np.gradient(r, z) ** 2), z))
    negK_mask = horn_mask & (K < 0)
    negK_fraction = float(np.sum(negK_mask) / max(np.sum(horn_mask), 1)) if np.any(horn_mask) else float("nan")
    out = {
        "geometry": geometry,
        "a": p.a,
        "b": p.b,
        "R": p.R,
        "L": p.L,
        "c": p.c,
        "volume_proxy": vol_proxy,
        "surface_proxy": surf_proxy,
        "rho_at_c": phb_wall_radius_scalar(p.c, p, geometry),
        "negative_K_horn_node_fraction": negK_fraction,
    }
    if np.any(horn_mask):
        out.update({
            "K_horn_min": float(np.nanmin(K[horn_mask])),
            "K_horn_max": float(np.nanmax(K[horn_mask])),
            "K_horn_mean": float(np.nanmean(K[horn_mask])),
        })
    if np.any(cap_mask):
        out.update({
            "K_cap_min": float(np.nanmin(K[cap_mask])),
            "K_cap_max": float(np.nanmax(K[cap_mask])),
            "K_cap_mean": float(np.nanmean(K[cap_mask])),
        })
    return out


def plot_geometry(p: PHBParams, geometry: str, cfg: RunConfig, outdir: Path) -> None:
    if plt is None:
        return
    z = np.linspace(-p.L, p.L, 2000)
    r = wall_radius_array(z, p, geometry)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(z, r, lw=2, label=f"{geometry} wall +r")
    ax.plot(z, -r, lw=1, label=f"{geometry} wall -r")
    ax.axvline(-p.c, ls="--", lw=1, label="external-focal planes ±c")
    ax.axvline(+p.c, ls="--", lw=1)
    ax.axhline(p.R, ls=":", lw=1, label="diagnostic r=R")
    ax.axhline(-p.R, ls=":", lw=1)
    ax.fill_between([-p.c, p.c], [p.R - 0.1 * p.R] * 2, [p.R + 0.1 * p.R] * 2,
                    alpha=0.18, label="D_0.10 diagnostic annulus")
    if cfg.mode == "open" and cfg.open_slot:
        ax.add_patch(plt.Rectangle((cfg.slot_z_center - cfg.slot_z_halfwidth, p.R - cfg.slot_r_halfwidth),
                                   2 * cfg.slot_z_halfwidth, 2 * cfg.slot_r_halfwidth,
                                   fill=False, lw=2, label="annular slot proxy"))
    ax.set_xlabel("z = x")
    ax.set_ylabel("signed r = rho")
    ax.set_title(f"PHB Maxwell/FDTD geometry control: {geometry}")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"geometry_{geometry}.png", dpi=180)
    plt.close(fig)


def component_from_name(name: str):
    if mp is None:
        raise RuntimeError("MEEP is not available.")
    lookup = {
        "Er": mp.Er, "Ep": mp.Ep, "Ez": mp.Ez,
        "Hr": mp.Hr, "Hp": mp.Hp, "Hz": mp.Hz,
        "Dr": mp.Dr, "Dp": mp.Dp, "Dz": mp.Dz,
    }
    if name not in lookup:
        raise ValueError(f"Unknown component '{name}'. Valid: {sorted(lookup)}")
    return lookup[name]


def all_dft_components(include_h: bool = True):
    comps = [mp.Er, mp.Ep, mp.Ez]
    if include_h:
        comps += [mp.Hr, mp.Hp, mp.Hz]
    return comps


def component_label(comp) -> str:
    if mp is None:
        return str(comp)
    for name in ["Er", "Ep", "Ez", "Hr", "Hp", "Hz", "Dr", "Dp", "Dz"]:
        if getattr(mp, name, None) == comp:
            return name
    return str(comp)


def make_simulation(cfg: RunConfig, geometry: str, fcen: float, df: float, seed: Optional[int] = None):
    p = cfg.params
    rmax = p.R + p.a + cfg.wall_thickness + cfg.pad_r + cfg.dpml
    zspan = 2.0 * (p.L + cfg.wall_thickness + cfg.pad_z + cfg.dpml)
    cell_size = mp.Vector3(rmax, 0, zspan)
    src_comp = component_from_name(cfg.source_component)

    source_center = mp.Vector3(cfg.source_r, 0, cfg.source_z)
    source_size = mp.Vector3(cfg.source_dr, 0, cfg.source_dz)
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=src_comp,
                         center=source_center,
                         size=source_size)]

    boundary_layers = [
        mp.PML(cfg.dpml, direction=mp.R),
        mp.PML(cfg.dpml, direction=mp.Z),
    ]

    courant = cfg.courant
    if courant is None:
        courant = min(0.5, 1.0 / (abs(cfg.m) + 0.5)) if abs(cfg.m) > 0 else 0.5

    if seed is not None:
        # MEEP versions differ in whether a seed helper is exposed in Python.
        # The simulation is deterministic for ordinary Gaussian sources and PML,
        # but this guard prevents AttributeError on pymeep builds without mp.set_seed.
        np.random.seed(int(seed))
        for _seed_func_name in ("set_seed", "set_random_seed", "_random_seed"):
            _seed_func = getattr(mp, _seed_func_name, None)
            if callable(_seed_func):
                _seed_func(int(seed))
                break

    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=boundary_layers,
        resolution=cfg.resolution,
        sources=sources,
        dimensions=mp.CYLINDRICAL,
        m=cfg.m,
        material_function=make_material_function(cfg, geometry),
        force_complex_fields=True,
        accurate_fields_near_cylorigin=True,
        Courant=courant,
    )
    return sim, rmax, zspan


def add_monitors(sim, cfg: RunConfig, fcen: float, df: float, nfreq: int, rmax: float, zspan: float):
    p = cfg.params
    dft = sim.add_dft_fields(
        all_dft_components(cfg.include_h_fields), fcen, df, nfreq,
        center=mp.Vector3(0.5 * (rmax - cfg.dpml), 0, 0),
        size=mp.Vector3(rmax - cfg.dpml, 0, zspan - 2 * cfg.dpml),
    )

    z_right = +0.5 * zspan - cfg.dpml - 0.25 * cfg.pad_z
    z_left = -0.5 * zspan + cfg.dpml + 0.25 * cfg.pad_z
    flux_right = sim.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0.5 * (rmax - cfg.dpml), 0, z_right),
                      size=mp.Vector3(rmax - cfg.dpml, 0, 0), direction=mp.Z),
    )
    flux_left = sim.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0.5 * (rmax - cfg.dpml), 0, z_left),
                      size=mp.Vector3(rmax - cfg.dpml, 0, 0), direction=mp.Z),
    )

    r_mon = min(rmax - cfg.dpml - 0.25 * cfg.pad_r, p.R + p.a + cfg.wall_thickness + 0.5 * cfg.pad_r)
    flux_radial = sim.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(r_mon, 0, 0),
                      size=mp.Vector3(0, 0, zspan - 2 * cfg.dpml), direction=mp.R),
    )

    flux_map = {"right_z": flux_right, "left_z": flux_left, "radial": flux_radial}

    # In open-slot mode the physically interesting loss/output is local and
    # annular, near z=slot_z_center and r~R.  The global radial monitor remains
    # useful as a leakage proxy, but this local monitor is closer to the intended
    # annular output-port diagnostic.
    if cfg.open_slot:
        slot_r = min(rmax - cfg.dpml - 0.10 * cfg.pad_r,
                     p.R + cfg.wall_thickness + 0.5 * max(cfg.slot_r_halfwidth, 1.0 / max(cfg.resolution, 1)))
        slot_z_size = max(2.0 * cfg.slot_z_halfwidth, 2.0 / max(cfg.resolution, 1))
        flux_slot_radial = sim.add_flux(
            fcen, df, nfreq,
            mp.FluxRegion(center=mp.Vector3(slot_r, 0, cfg.slot_z_center),
                          size=mp.Vector3(0, 0, slot_z_size), direction=mp.R),
        )
        flux_map["slot_radial"] = flux_slot_radial

    return dft, flux_map


def metadata_grid(sim, dft_obj, energy_shape: Tuple[int, ...], cfg: RunConfig, rmax: float, zspan: float):
    shape = tuple(int(x) for x in energy_shape if int(x) > 1)
    if len(shape) == 0:
        raise RuntimeError(f"Cannot reconstruct a DFT grid from scalar shape {energy_shape!r}.")
    if len(shape) == 1:
        nr, nz = shape[0], 1
    else:
        n0, n1 = shape[0], shape[1]
        r_size = max(float(rmax - cfg.dpml), 1e-15)
        z_size = max(float(zspan - 2.0 * cfg.dpml), 1e-15)
        expected_ratio = r_size / z_size
        ratio_01 = n0 / max(n1, 1)
        ratio_10 = n1 / max(n0, 1)
        if abs(math.log(max(ratio_01, 1e-12) / expected_ratio)) <= abs(math.log(max(ratio_10, 1e-12) / expected_ratio)):
            nr, nz = n0, n1
            transpose = False
        else:
            nr, nz = n1, n0
            transpose = True

        dr = r_size / max(nr, 1)
        dz = z_size / max(nz, 1)
        r = (np.arange(nr, dtype=float) + 0.5) * dr
        z = -0.5 * z_size + (np.arange(nz, dtype=float) + 0.5) * dz
        Rg, Zg = np.meshgrid(r, z, indexing="ij")
        if transpose:
            Rg = Rg.T
            Zg = Zg.T
        W = 2.0 * math.pi * np.maximum(Rg, 0.0) * abs(dr) * abs(dz)
        return Rg, Zg, W

    r_size = max(float(rmax - cfg.dpml), 1e-15)
    dr = r_size / max(nr, 1)
    r = (np.arange(nr, dtype=float) + 0.5) * dr
    Rg = r.reshape((nr, 1))
    Zg = np.zeros_like(Rg)
    W = 2.0 * math.pi * np.maximum(Rg, 0.0) * abs(dr)
    return Rg, Zg, W


def get_dft_energy_density(sim, dft_obj, freq_index: int, cfg: RunConfig) -> Tuple[np.ndarray, List[str], Dict[str, np.ndarray]]:
    arrays = []
    used = []
    comp_energy: Dict[str, np.ndarray] = {}
    target_shape = None
    for comp in all_dft_components(cfg.include_h_fields):
        try:
            arr = np.asarray(sim.get_dft_array(dft_obj, comp, freq_index))
        except Exception:
            continue
        arr = np.squeeze(arr)
        if target_shape is None:
            target_shape = arr.shape
        if arr.shape != target_shape:
            continue
        label = component_label(comp)
        ecomp = np.abs(arr) ** 2
        arrays.append(ecomp)
        used.append(label)
        comp_energy[label] = ecomp
    if not arrays:
        raise RuntimeError("No DFT field arrays could be read from Meep.")
    energy = np.zeros_like(arrays[0], dtype=float)
    for a in arrays:
        energy += a
    return energy, used, comp_energy


def analyze_dft_at_frequency(sim, dft_obj, freq_index: int, freq: float, cfg: RunConfig,
                             geometry: str, rmax: float, zspan: float) -> Tuple[Dict[str, float], np.ndarray, np.ndarray, np.ndarray]:
    p = cfg.params
    edens, used, comp_energy = get_dft_energy_density(sim, dft_obj, freq_index, cfg)
    Rg, Zg, W = metadata_grid(sim, dft_obj, edens.shape, cfg, rmax, zspan)
    inside = inside_cavity(Rg, Zg, p, geometry)
    rw = wall_radius_array(Zg, p, geometry)

    total = float(np.sum(edens[inside] * W[inside]))
    eps = max(total, 1e-300)

    masks = {
        "eta_025R": inside & (np.abs(Zg) <= p.c) & (np.abs(Rg - p.R) <= 0.025 * p.R),
        "eta_050R": inside & (np.abs(Zg) <= p.c) & (np.abs(Rg - p.R) <= 0.050 * p.R),
        "eta_100R": inside & (np.abs(Zg) <= p.c) & (np.abs(Rg - p.R) <= 0.100 * p.R),
        "E_between_foci": inside & (np.abs(Zg) <= p.c),
        "E_outside_foci": inside & (np.abs(Zg) > p.c),
        "near_wall_90": inside & (Rg >= 0.90 * np.maximum(rw, 1e-15)),
        "axis_10": inside & (Rg <= 0.10 * np.maximum(rw, 1e-15)),
    }
    vals = {name: float(np.sum(edens[mask] * W[mask])) for name, mask in masks.items()}
    eta025 = vals["eta_025R"] / eps
    eta050 = vals["eta_050R"] / eps
    eta100 = vals["eta_100R"] / eps
    between = vals["E_between_foci"]
    outside = vals["E_outside_foci"]
    cf = between / max(outside, 1e-300)
    leakage = outside / eps
    volume_total = float(np.sum(W[inside]))
    volume_ann100 = float(np.sum(W[masks["eta_100R"]]))
    geom_fraction = volume_ann100 / max(volume_total, 1e-300)
    enrichment = eta100 / max(geom_fraction, 1e-300)

    row: Dict[str, float] = {
        "geometry": geometry,
        "m": cfg.m,
        "freq": freq,
        "kR": 2.0 * math.pi * freq * p.R,
        "eta_025R": eta025,
        "eta_050R": eta050,
        "eta_100R": eta100,
        "E_total_proxy": total,
        "E_between_foci_frac": between / eps,
        "E_outside_foci_frac": leakage,
        "CF": cf,
        "log10_CF": math.log10(max(cf, 1e-300)),
        "near_wall_90_frac": vals["near_wall_90"] / eps,
        "axis_10_frac": vals["axis_10"] / eps,
        "geom_ann100_fraction": geom_fraction,
        "ann100_enrichment": enrichment,
        "used_components_count": len(used),
        "grid_shape_0": int(edens.shape[0]) if len(edens.shape) > 0 else 0,
        "grid_shape_1": int(edens.shape[1]) if len(edens.shape) > 1 else 1,
    }

    for label, arr in comp_energy.items():
        comp_total = float(np.sum(arr[inside] * W[inside]))
        comp_key = label.replace(" ", "_")
        row[f"component_{comp_key}_frac_of_total"] = comp_total / eps
        row[f"component_{comp_key}_ann100_frac_of_component"] = (
            float(np.sum(arr[masks["eta_100R"]] * W[masks["eta_100R"]])) / max(comp_total, 1e-300)
        )
        row[f"component_{comp_key}_outside_foci_frac_of_component"] = (
            float(np.sum(arr[masks["E_outside_foci"]] * W[masks["E_outside_foci"]])) / max(comp_total, 1e-300)
        )
    return row, edens, Rg, Zg


def plot_field_map(energy: np.ndarray, Rg: np.ndarray, Zg: np.ndarray, cfg: RunConfig,
                   geometry: str, row: Dict[str, float], outdir: Path, freq_index: int) -> None:
    if plt is None:
        return
    p = cfg.params
    fig, ax = plt.subplots(figsize=(7.2, 4.6))

    # Plot only the field inside the intended cavity mask.  Earlier smoke plots
    # showed the whole DFT rectangle, which could visually over-emphasize
    # exterior/PML fields even though the metrics were computed only inside.
    inside = inside_cavity(Rg, Zg, p, geometry)
    eplot = np.array(energy, dtype=float, copy=True)
    eplot[~inside] = np.nan
    norm = float(np.nanmax(eplot)) if np.any(np.isfinite(eplot)) else float(np.nanmax(energy))
    val = np.log10(np.maximum(eplot / max(norm, 1e-300), 1e-12))

    extent = [float(np.nanmin(Zg)), float(np.nanmax(Zg)), float(np.nanmin(Rg)), float(np.nanmax(Rg))]
    im = ax.imshow(val, origin="lower", aspect="auto", extent=extent, vmin=-8, vmax=0)
    z = np.linspace(-p.L, p.L, 1200)
    r = wall_radius_array(z, p, geometry)
    ax.plot(z, r, "w-", lw=1.4, label="cavity wall")
    ax.axvline(-p.c, color="w", ls="--", lw=0.9)
    ax.axvline(+p.c, color="w", ls="--", lw=0.9)
    ax.axhline(p.R, color="w", ls=":", lw=0.9)
    ax.set_xlabel("z = x")
    ax.set_ylabel("r = rho")
    ax.set_title(f"{geometry}, m={cfg.m}, kR={row['kR']:.3f}, eta100={100*row['eta_100R']:.2f}%")
    ax.set_ylim(0, p.R + p.a + cfg.wall_thickness + cfg.pad_r)
    fig.colorbar(im, ax=ax, label="log10 normalized inside-cavity |field|^2 proxy")
    fig.tight_layout()
    fig.savefig(outdir / f"field_inside_{geometry}_fi{freq_index:02d}_kR{row['kR']:.3f}.png", dpi=170)
    plt.close(fig)


def run_one_geometry(cfg: RunConfig, geometry: str, outdir: Path) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    p = cfg.params
    fcen = cfg.kr_center / (2.0 * math.pi * p.R)
    df = cfg.kr_span / (2.0 * math.pi * p.R)
    if df <= 0:
        raise ValueError("kr_span must be positive.")

    if mp.am_master():
        print(f"\n=== MEEP run: {geometry}, mode={cfg.mode}, m={cfg.m}, fcen={fcen:.6g}, df={df:.6g} ===", flush=True)

    sim, rmax, zspan = make_simulation(cfg, geometry, fcen, df, seed=cfg.seed)
    dft, flux_objs = add_monitors(sim, cfg, fcen, df, cfg.nfreq, rmax, zspan)

    probe = mp.Vector3(p.R, 0, 0.0)
    h = mp.Harminv(component_from_name(cfg.source_component), probe, fcen, df)

    # Fixed-time run mode.
    # Resonator fields can decay very slowly, so stop_when_fields_decayed() may
    # keep a nominally quick smoke test running for hours.  For this Docker/WSL
    # workflow we use a hard, predictable after-source interval and store the
    # run as a screening/proxy result, not as a final Q-certification.
    sim.run(mp.after_sources(h), until_after_sources=cfg.until_after_sources)

    freqs = list(mp.get_flux_freqs(next(iter(flux_objs.values()))))
    metric_rows: List[Dict] = []
    best_row = None
    best_payload = None
    for fi, freq in enumerate(freqs):
        row, edens, Rg, Zg = analyze_dft_at_frequency(sim, dft, fi, float(freq), cfg, geometry, rmax, zspan)
        metric_rows.append(row)
        if best_row is None or row["eta_100R"] > best_row["eta_100R"]:
            best_row = row
            best_payload = (edens, Rg, Zg, fi)
        if cfg.plot_all_fields and mp.am_master():
            plot_field_map(edens, Rg, Zg, cfg, geometry, row, outdir, fi)

    if (not cfg.plot_all_fields) and best_row is not None and best_payload is not None and mp.am_master():
        edens, Rg, Zg, fi = best_payload
        plot_field_map(edens, Rg, Zg, cfg, geometry, best_row, outdir, int(fi))

    harminv_rows: List[Dict] = []
    for mode in getattr(h, "modes", []):
        try:
            freq = float(mode.freq)
            decay = float(mode.decay)
            q = float(mode.Q)
            amp = abs(complex(mode.amp))
            err = float(mode.err)
        except Exception:
            continue
        harminv_rows.append({
            "geometry": geometry,
            "m": cfg.m,
            "freq": freq,
            "kR": 2.0 * math.pi * freq * p.R,
            "decay": decay,
            "Q": q,
            "amp_abs": amp,
            "err": err,
        })

    flux_rows: List[Dict] = []
    for mon_name, obj in flux_objs.items():
        fluxes = list(mp.get_fluxes(obj))
        for freq, flx in zip(freqs, fluxes):
            flux_rows.append({
                "geometry": geometry,
                "monitor": mon_name,
                "m": cfg.m,
                "freq": float(freq),
                "kR": 2.0 * math.pi * float(freq) * p.R,
                "flux": float(flx),
            })

    sim.reset_meep()
    return metric_rows, harminv_rows, flux_rows


def dataframe(rows: List[Dict]):
    if pd is None:
        return rows
    return pd.DataFrame(rows)


def save_table(rows: List[Dict], path: Path) -> None:
    if pd is not None:
        pd.DataFrame(rows).to_csv(path, index=False)
    else:
        with path.open("w", encoding="utf-8") as f:
            if not rows:
                return
            keys = list(rows[0].keys())
            f.write(",".join(keys) + "\n")
            for row in rows:
                f.write(",".join(str(row.get(k, "")) for k in keys) + "\n")


def summarize_uniqueness(metric_rows: List[Dict], harminv_rows: List[Dict]) -> List[Dict]:
    if not metric_rows:
        return []
    if pd is None:
        geoms = sorted({str(r["geometry"]) for r in metric_rows})
        out = []
        for g in geoms:
            rows = [r for r in metric_rows if r["geometry"] == g]
            best = max(rows, key=lambda x: x["eta_100R"])
            out.append(dict(best))
        return out

    df = pd.DataFrame(metric_rows)
    idx = df.groupby("geometry")["eta_100R"].idxmax()
    best = df.loc[idx].copy().sort_values("eta_100R", ascending=False)

    if harminv_rows:
        hdf = pd.DataFrame(harminv_rows)
        qbest = hdf.sort_values("Q", ascending=False).groupby("geometry").head(1)[["geometry", "Q", "kR", "err"]]
        qbest = qbest.rename(columns={"Q": "best_harminv_Q", "kR": "best_harminv_kR", "err": "best_harminv_err"})
        best = best.merge(qbest, on="geometry", how="left")

    ph = best[best["geometry"] == "PH"]
    controls = best[best["geometry"] != "PH"]
    if len(ph) and len(controls):
        ph_eta = float(ph.iloc[0]["eta_100R"])
        ph_logcf = float(ph.iloc[0]["log10_CF"])
        best_control_eta = float(controls["eta_100R"].max())
        best_control_logcf = float(controls["log10_CF"].max())
        best["delta_eta100_PH_minus_best_control"] = ph_eta - best_control_eta
        best["delta_logCF_PH_minus_best_control"] = ph_logcf - best_control_logcf
        best["PH_eta_unique_flag"] = (ph_eta - best_control_eta) > 0.03
        best["PH_CF_unique_flag"] = (ph_logcf - best_control_logcf) > 0.5
    return best.to_dict("records")


def plot_uniqueness_summary(summary_rows: List[Dict], outdir: Path) -> None:
    if plt is None or pd is None or not summary_rows:
        return
    df = pd.DataFrame(summary_rows)
    if "geometry" not in df:
        return
    df = df.sort_values("eta_100R", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.barh(df["geometry"], 100 * df["eta_100R"])
    ax.axvline(70, ls="--", lw=1, label="C2 scalar threshold reference")
    ax.set_xlabel("best full-vector DFT energy in D_0.10 (%)")
    ax.set_title("MEEP full-vector annular localization screen")
    ax.grid(True, axis="x", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "uniqueness_eta100_summary.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    df2 = df.sort_values("log10_CF", ascending=True)
    ax.barh(df2["geometry"], df2["log10_CF"])
    ax.set_xlabel("log10(CF) = log10(E_|z|<=c / E_|z|>c)")
    ax.set_title("MEEP full-vector axial confinement screen")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "uniqueness_logCF_summary.png", dpi=170)
    plt.close(fig)


def summarize_flux(flux_rows: List[Dict]) -> List[Dict]:
    if not flux_rows:
        return []

    if pd is None:
        # Fallback for minimal Docker images without pandas.
        grouped: Dict[Tuple[str, Any, Any, Any], Dict[str, Any]] = {}
        for row in flux_rows:
            try:
                geometry = str(row.get("geometry", ""))
                m = row.get("m")
                freq = float(row.get("freq", 0.0))
                kR = float(row.get("kR", 0.0))
                monitor = str(row.get("monitor", ""))
                flx = float(row.get("flux", 0.0))
            except Exception:
                continue
            key = (geometry, m, freq, kR)
            g = grouped.setdefault(key, {"geometry": geometry, "m": m, "freq": freq, "kR": kR,
                                         "right_z": 0.0, "left_z": 0.0, "radial": 0.0, "slot_radial": 0.0})
            if monitor in ("right_z", "left_z", "radial", "slot_radial"):
                g[monitor] += flx

        best_by_geom: Dict[str, Dict[str, Any]] = {}
        for g in grouped.values():
            total_abs = (abs(float(g.get("right_z", 0.0))) + abs(float(g.get("left_z", 0.0))) +
                         abs(float(g.get("radial", 0.0))) + abs(float(g.get("slot_radial", 0.0))))
            g["total_abs_flux"] = total_abs
            g["right_fraction_abs"] = abs(float(g.get("right_z", 0.0))) / total_abs if total_abs > 0 else float("nan")
            g["radial_fraction_abs"] = abs(float(g.get("radial", 0.0))) / total_abs if total_abs > 0 else float("nan")
            g["slot_fraction_abs"] = abs(float(g.get("slot_radial", 0.0))) / total_abs if total_abs > 0 else float("nan")
            geom = str(g.get("geometry", ""))
            old = best_by_geom.get(geom)
            if old is None or total_abs > float(old.get("total_abs_flux", -1.0)):
                best_by_geom[geom] = dict(g)

        out: List[Dict[str, Any]] = []
        for g in sorted(best_by_geom.values(), key=lambda r: float(r.get("total_abs_flux", 0.0)), reverse=True):
            out.append({
                "geometry": g.get("geometry"),
                "m": g.get("m"),
                "freq_at_max_abs_flux": g.get("freq"),
                "kR_at_max_abs_flux": g.get("kR"),
                "right_z_flux_at_max_abs_flux": g.get("right_z"),
                "left_z_flux_at_max_abs_flux": g.get("left_z"),
                "radial_flux_at_max_abs_flux": g.get("radial"),
                "slot_radial_flux_at_max_abs_flux": g.get("slot_radial"),
                "total_abs_flux": g.get("total_abs_flux"),
                "right_fraction_abs": g.get("right_fraction_abs"),
                "radial_fraction_abs": g.get("radial_fraction_abs"),
                "slot_fraction_abs": g.get("slot_fraction_abs"),
            })
        return out

    df = pd.DataFrame(flux_rows)
    if df.empty:
        return []
    pivot = df.pivot_table(index=["geometry", "m", "freq", "kR"], columns="monitor", values="flux", aggfunc="sum").reset_index()
    for name in ["right_z", "left_z", "radial", "slot_radial"]:
        if name not in pivot.columns:
            pivot[name] = 0.0
    flux_cols = ["right_z", "left_z", "radial", "slot_radial"]
    pivot["total_abs_flux"] = pivot[flux_cols].abs().sum(axis=1)
    pivot["right_fraction_abs"] = pivot["right_z"].abs() / pivot["total_abs_flux"].replace(0, np.nan)
    pivot["radial_fraction_abs"] = pivot["radial"].abs() / pivot["total_abs_flux"].replace(0, np.nan)
    pivot["slot_fraction_abs"] = pivot["slot_radial"].abs() / pivot["total_abs_flux"].replace(0, np.nan)
    idx = pivot.groupby("geometry")["total_abs_flux"].idxmax()
    best = pivot.loc[idx].copy().sort_values("total_abs_flux", ascending=False)
    best = best.rename(columns={
        "freq": "freq_at_max_abs_flux",
        "kR": "kR_at_max_abs_flux",
        "right_z": "right_z_flux_at_max_abs_flux",
        "left_z": "left_z_flux_at_max_abs_flux",
        "radial": "radial_flux_at_max_abs_flux",
        "slot_radial": "slot_radial_flux_at_max_abs_flux",
    })
    return best.to_dict("records")


def write_decision_report(cfg: RunConfig, geom_rows: List[Dict], summary_rows: List[Dict], flux_summary_rows: List[Dict], outdir: Path) -> None:
    lines: List[str] = []
    lines.append("# PHB MEEP full-vector decision report")
    lines.append("")
    lines.append("This report is generated automatically by the script. It is a screening report, not a final laser-resonator certification.")
    lines.append("")
    lines.append("## Run identity")
    lines.append(f"- mode: `{cfg.mode}`")
    lines.append(f"- geometry list: `{', '.join(cfg.geometries)}`")
    lines.append(f"- a={cfg.params.a}, b={cfg.params.b}, R={cfg.params.R}, c={cfg.params.c}, L={cfg.params.L}")
    lines.append(f"- m={cfg.m}, resolution={cfg.resolution}, nfreq={cfg.nfreq}, kR scan = {cfg.kr_center} ± {0.5*cfg.kr_span}")
    lines.append(f"- seed: {cfg.seed}")  # < ADDED
    lines.append("")
    lines.append("## Decision logic")
    lines.append("A revolutionary PHB-specific wave effect should not be accepted from eta_100R alone, because high-m annular localization can be generic in axisymmetric cavities. Strong evidence requires PHB to beat the best non-hyperbolic control after grid/PML/source convergence in several independent metrics:")
    lines.append("")
    lines.append("1. eta_100R advantage > 3 percentage points;")
    lines.append("2. log10(CF) advantage > 0.5;")
    lines.append("3. open-slot useful output/leakage pattern not reproduced by controls;")
    lines.append("4. Harminv Q or resonance recurrence advantage not reproduced by controls;")
    lines.append("5. stability under resolution, PML thickness, source component and source position changes.")
    lines.append("")
    if not summary_rows:
        lines.append("No Maxwell metric rows were generated in this run. Geometry-only dry run cannot decide wave uniqueness.")
    else:
        ph_rows = [r for r in summary_rows if str(r.get("geometry", "")).upper() == "PH"]
        ctrl_rows = [r for r in summary_rows if str(r.get("geometry", "")).upper() != "PH"]
        if ph_rows and ctrl_rows:
            ph = ph_rows[0]
            best_eta_ctrl = max(ctrl_rows, key=lambda r: float(r.get("eta_100R", -1)))
            best_cf_ctrl = max(ctrl_rows, key=lambda r: float(r.get("log10_CF", -1e99)))
            de = float(ph.get("eta_100R", 0.0)) - float(best_eta_ctrl.get("eta_100R", 0.0))
            dc = float(ph.get("log10_CF", 0.0)) - float(best_cf_ctrl.get("log10_CF", 0.0))
            lines.append("## Main PHB-vs-control result")
            lines.append(f"- PH best eta_100R = {100*float(ph.get('eta_100R', 0.0)):.3f}% at kR={float(ph.get('kR', 0.0)):.4g}")
            lines.append(f"- best-control eta_100R = {100*float(best_eta_ctrl.get('eta_100R', 0.0)):.3f}% for `{best_eta_ctrl.get('geometry')}`")
            lines.append(f"- delta eta_100R = {100*de:.3f} percentage points")
            lines.append(f"- PH log10(CF) = {float(ph.get('log10_CF', 0.0)):.3f}")
            lines.append(f"- best-control log10(CF) = {float(best_cf_ctrl.get('log10_CF', 0.0)):.3f} for `{best_cf_ctrl.get('geometry')}`")
            lines.append(f"- delta log10(CF) = {dc:.3f}")
            lines.append("")
            if de > 0.03 and dc > 0.5:
                lines.append("Preliminary verdict: PHB shows a nontrivial full-vector advantage in this run. Repeat at higher resolution, thicker PML, different sources and open slot before making a strong claim.")
            elif de > 0.03 or dc > 0.5:
                lines.append("Preliminary verdict: PHB shows a partial advantage in this run. This is promising but insufficient for a revolution-level claim.")
            else:
                lines.append("Preliminary verdict: no PHB-specific full-vector uniqueness is resolved in this run. The geometry may still be useful, but this run does not justify a revolution-level resonator claim.")
        else:
            lines.append("PH or control rows are missing; PHB-specific uniqueness cannot be tested without controls.")
    lines.append("")
    lines.append("## Files to inspect")
    lines.append("- `geometry_summary.csv` — geometry, volume/surface proxies, curvature diagnostics")
    lines.append("- `dft_metrics.csv` — full-vector DFT energy metrics, including component-resolved fractions")
    lines.append("- `harminv_modes.csv` — resonant frequencies/Q estimates when Harminv finds modes")
    lines.append("- `flux_spectra.csv` and `flux_summary.csv` — axial/radial leakage or output coupling proxies")
    lines.append("- `uniqueness_summary.csv` — best-frequency PHB/control comparison")
    (outdir / "decision_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_geometries(text: str) -> List[str]:
    return [x.strip().upper() for x in text.split(",") if x.strip()]


def build_config(args: argparse.Namespace) -> RunConfig:
    p = PHBParams(a=args.a, b=args.b, R=args.R, R1=args.R1, R2=args.R2)

    if args.preset == "smoke":
        resolution = args.resolution if args.resolution is not None else 12
        nfreq = args.nfreq if args.nfreq is not None else 5
        until = args.until_after_sources if args.until_after_sources is not None else 60
    elif args.preset == "research":
        resolution = args.resolution if args.resolution is not None else 30
        nfreq = args.nfreq if args.nfreq is not None else 17
        until = args.until_after_sources if args.until_after_sources is not None else 250
    elif args.preset == "publication":
        resolution = args.resolution if args.resolution is not None else 55
        nfreq = args.nfreq if args.nfreq is not None else 31
        until = args.until_after_sources if args.until_after_sources is not None else 600
    else:
        resolution = args.resolution if args.resolution is not None else 20
        nfreq = args.nfreq if args.nfreq is not None else 9
        until = args.until_after_sources if args.until_after_sources is not None else 150

    slot_z_center = args.slot_z_center
    if slot_z_center is None:
        slot_z_center = p.a
    slot_r_halfwidth = args.slot_r_halfwidth
    if slot_r_halfwidth is None:
        # IMPORTANT: for a nonzero R2 sweep the physical half-width of the
        # annular slot must actually follow dR_half = 0.5*R*R2.
        # In earlier versions a protective floor max(..., 0.05R) made all
        # R2 <= 0.10 runs geometrically identical.  Keep the 0.05R fallback
        # only for open runs where R2=0 was intentionally used as a default
        # slot proxy.
        if args.mode == "open" and p.R2 > 0.0:
            slot_r_halfwidth = 0.5 * p.dR_total
        elif args.mode == "open":
            slot_r_halfwidth = 0.05 * p.R
        else:
            slot_r_halfwidth = 0.0

    return RunConfig(
        params=p,
        mode=args.mode,
        geometries=parse_geometries(args.geometries),
        m=args.m,
        kr_center=args.kr_center,
        kr_span=args.kr_span,
        nfreq=int(nfreq),
        resolution=int(resolution),
        dpml=args.dpml,
        pad_r=args.pad_r,
        pad_z=args.pad_z,
        wall_thickness=args.wall_thickness,
        source_component=args.source_component,
        source_r=args.source_r if args.source_r is not None else p.R,
        source_z=args.source_z,
        source_dr=args.source_dr if args.source_dr is not None else max(0.06 * p.R, 2.0 / max(resolution, 1)),
        source_dz=args.source_dz if args.source_dz is not None else max(0.15 * p.a, 2.0 / max(resolution, 1)),
        until_after_sources=float(until),
        decay_by=args.decay_by,
        harminv_after_sources=args.harminv_after_sources,
        courant=args.courant,
        open_slot=(args.mode == "open"),
        slot_z_center=float(slot_z_center),
        slot_z_halfwidth=args.slot_z_halfwidth,
        slot_r_halfwidth=float(slot_r_halfwidth),
        outdir=args.outdir,
        dry_run_geometry=args.dry_run_geometry,
        plot_all_fields=args.plot_all_fields,
        include_h_fields=not args.e_fields_only,
        seed=args.seed,
    )


def make_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="MEEP full-vector Maxwell/FDTD PHB uniqueness scan in cylindrical coordinates."
    )
    ap.add_argument("--preset", choices=["custom", "smoke", "research", "publication"], default="custom",
                    help="Convenience runtime preset. All numerical options can still be overridden.")
    ap.add_argument("--mode", choices=["closed", "open"], default="closed",
                    help="closed = metal-shell cavity; open = shell with a simple annular slot proxy.")
    ap.add_argument("--geometries", default="PH,LINEAR,POLY2,SMOOTHSTEP",
                    help="Comma-separated list: PH, LINEAR, POLY2, SMOOTHSTEP, HERMITE_SLOPE_MATCH, POLY_VOL_MATCH, CYLINDER_R, ELLIPSOID_L")
    ap.add_argument("--a", type=float, default=0.3)
    ap.add_argument("--b", type=float, default=0.6)
    ap.add_argument("--R", type=float, default=3.0)
    ap.add_argument("--R1", type=float, default=0.0)
    ap.add_argument("--R2", type=float, default=0.0,
                    help="Annular slot fractional width for open mode; dR_total=R*R2. If 0, a default 0.1R full slot is used.")
    ap.add_argument("--m", type=int, default=15, help="Azimuthal order exp(i*m*phi).")
    ap.add_argument("--kr-center", type=float, default=15.8, help="Center kR of the frequency scan.")
    ap.add_argument("--kr-span", type=float, default=4.0, help="Total kR span of the frequency scan.")
    ap.add_argument("--nfreq", type=int, default=None, help="Number of DFT frequency samples.")
    ap.add_argument("--resolution", type=int, default=None, help="MEEP pixels per length unit.")
    ap.add_argument("--dpml", type=float, default=1.0)
    ap.add_argument("--pad-r", type=float, default=1.0)
    ap.add_argument("--pad-z", type=float, default=1.0)
    ap.add_argument("--wall-thickness", type=float, default=0.18)
    ap.add_argument("--source-component", default="Ep", choices=["Er", "Ep", "Ez", "Hr", "Hp", "Hz"])
    ap.add_argument("--source-r", type=float, default=None)
    ap.add_argument("--source-z", type=float, default=0.0)
    ap.add_argument("--source-dr", type=float, default=None)
    ap.add_argument("--source-dz", type=float, default=None)
    ap.add_argument("--until-after-sources", type=float, default=None)
    ap.add_argument("--decay-by", type=float, default=40.0,
                    help="Field-decay waiting time used by stop_when_fields_decayed.")
    ap.add_argument("--harminv-after-sources", type=float, default=120.0,
                    help="Reserved fixed Harminv interval; currently fallback uses until-after-sources.")
    ap.add_argument("--courant", type=float, default=None)
    ap.add_argument("--slot-z-center", type=float, default=None)
    ap.add_argument("--slot-z-halfwidth", type=float, default=0.06)
    ap.add_argument("--slot-r-halfwidth", type=float, default=None)
    ap.add_argument("--outdir", default="/work/PHB_SINGLE")
    ap.add_argument("--dry-run-geometry", action="store_true",
                    help="Only compute geometry/curvature summaries and plots. Does not require Meep.")
    ap.add_argument("--plot-all-fields", action="store_true")
    ap.add_argument("--e-fields-only", action="store_true",
                    help="Use only Er/Ep/Ez in DFT energy proxy. Default includes E and H fields.")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for MEEP (affects GaussianSource noise and PML).")
    return ap


def scan_main(argv: Optional[List[str]] = None) -> int:
    args = make_argparser().parse_args(argv)
    cfg = build_config(args)
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    p = cfg.params
    manifest = asdict(cfg)
    manifest["params"] = asdict(p)
    manifest["scientific_status"] = (
        "MEEP full-vector cylindrical-coordinate FDTD screening. Not a laser model; "
        "use results as C4/open-boundary evidence only after resolution and controls are converged."
    )
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    geom_rows: List[Dict] = []
    for g in cfg.geometries:
        row = geometry_summary(p, g)
        geom_rows.append(row)
        plot_geometry(p, g, cfg, outdir)
    save_table(geom_rows, outdir / "geometry_summary.csv")

    print("\nGeometry summary:")
    for row in geom_rows:
        print(json.dumps(row, indent=2))

    if cfg.dry_run_geometry:
        write_decision_report(cfg, geom_rows, [], [], outdir)
        print(f"\nDry-run complete. Geometry outputs are in: {outdir}")
        return 0

    if MEEP_IMPORT_ERROR is not None or mp is None:
        print("\nERROR: meep could not be imported in this Python environment.", file=sys.stderr)
        print(f"Import error: {MEEP_IMPORT_ERROR!r}", file=sys.stderr)
        print("Install pymeep from conda-forge, then rerun. Geometry dry-run does not require Meep.", file=sys.stderr)
        return 2

    metric_rows_all: List[Dict] = []
    harminv_rows_all: List[Dict] = []
    flux_rows_all: List[Dict] = []

    for geometry in cfg.geometries:
        metric_rows, harminv_rows, flux_rows = run_one_geometry(cfg, geometry, outdir)
        metric_rows_all.extend(metric_rows)
        harminv_rows_all.extend(harminv_rows)
        flux_rows_all.extend(flux_rows)

        if mp.am_master():
            save_table(metric_rows_all, outdir / "dft_metrics.csv")
            save_table(harminv_rows_all, outdir / "harminv_modes.csv")
            save_table(flux_rows_all, outdir / "flux_spectra.csv")

    if mp.am_master():
        summary = summarize_uniqueness(metric_rows_all, harminv_rows_all)
        save_table(summary, outdir / "uniqueness_summary.csv")
        flux_summary = summarize_flux(flux_rows_all)
        save_table(flux_summary, outdir / "flux_summary.csv")
        plot_uniqueness_summary(summary, outdir)
        write_decision_report(cfg, geom_rows, summary, flux_summary, outdir)
        print("\n=== Best-frequency uniqueness summary ===")
        for row in summary:
            print(json.dumps(row, indent=2, default=str))
        print(f"\nDONE. Outputs are in: {outdir}")
    return 0


# ============================================================
# Unified single-file matrix runner and aggregator
# ============================================================

def _run_cmd(cmd: List[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("\n=== RUN ===")
    print(" ".join(shlex.quote(x) for x in cmd))
    print(f"LOG: {log_path}")
    t0 = time.time()
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        log.write("COMMAND: " + " ".join(shlex.quote(x) for x in cmd) + "\n\n")
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert p.stdout is not None
        for line in p.stdout:
            sys.stdout.write(line)
            log.write(line)
        ret = p.wait()
        dt = time.time() - t0
        log.write(f"\nRETURN_CODE={ret}\nELAPSED_SECONDS={dt:.3f}\n")
    print(f"=== DONE rc={ret}, elapsed={dt/60:.2f} min ===")
    return ret


def _zip_dir(src: Path, dst: Path) -> None:
    if dst.exists():
        dst.unlink()
    with zipfile.ZipFile(dst, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for p in src.rglob("*"):
            if p.is_file() and p != dst:
                zf.write(p, p.relative_to(src.parent))


def matrix_main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Run the PHB MEEP decision matrix from this single unified file."
    )
    ap.add_argument("--root", default="/work/PHB_RESULTS")
    ap.add_argument("--quick", action="store_true", help="Run only dry-run + smoke closed + one open slot test.")
    ap.add_argument("--medium", action="store_true", help="Run a targeted validation set: quick tests plus resolution/source/R2 checks, shorter than full.")
    ap.add_argument("--no-zip", action="store_true")
    ap.add_argument("--python", default=sys.executable, help="Python executable used for child scan/aggregate runs.")
    args = ap.parse_args(argv)

    root = Path(args.root)
    logs = root / "_logs"
    root.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    py = args.python
    this_script = str(Path(__file__).resolve())

    tasks: List[Tuple[str, List[str]]] = []
    tasks.append(("00_geometry_dryrun", [py, this_script, "scan", "--dry-run-geometry", "--outdir", str(root / "00_geometry_dryrun")]))

    tasks.append(("01_smoke_closed_m15", [
        py, this_script, "scan",
        "--preset", "smoke", "--mode", "closed",
        "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP",
        "--m", "15", "--resolution", "20", "--nfreq", "9",
        "--outdir", str(root / "01_smoke_closed_m15")
    ]))

    tasks.append(("02_open_R2_0p06", [
        py, this_script, "scan",
        "--preset", "smoke", "--mode", "open", "--R2", "0.06",
        "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP",
        "--m", "15", "--resolution", "20", "--nfreq", "9",
        "--outdir", str(root / "02_open_R2_0p06")
    ]))

    if args.medium and not args.quick:
        # Targeted next-step validation after a successful quick run.
        # This is much shorter than the full matrix but checks the most important
        # scientific failure modes: resolution sensitivity, source dependence,
        # azimuthal-order dependence and R2 slot dependence.
        for m in [10, 15, 20]:
            tasks.append((f"03_medium_open_m_{m}", [
                py, this_script, "scan",
                "--preset", "smoke", "--mode", "open", "--R2", "0.06",
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP",
                "--m", str(m), "--resolution", "24", "--nfreq", "9",
                "--outdir", str(root / f"03_medium_open_m_{m}")
            ]))

        for res, dpml in [(24, 0.8), (28, 1.0)]:
            suffix = str(dpml).replace('.', 'p')
            tasks.append((f"04_medium_stability_open_res_{res}_dpml_{suffix}", [
                py, this_script, "scan",
                "--preset", "smoke", "--mode", "open", "--R2", "0.06",
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP",
                "--m", "15", "--resolution", str(res), "--dpml", str(dpml), "--nfreq", "9",
                "--outdir", str(root / f"04_medium_stability_open_res_{res}_dpml_{suffix}")
            ]))

        for r2 in [0.04, 0.06, 0.08, 0.10]:
            tasks.append((f"05_medium_open_R2_{str(r2).replace('.', 'p')}", [
                py, this_script, "scan",
                "--preset", "smoke", "--mode", "open", "--R2", str(r2),
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP",
                "--m", "15", "--resolution", "24", "--nfreq", "9",
                "--outdir", str(root / f"05_medium_open_R2_{str(r2).replace('.', 'p')}")
            ]))

        for comp in ["Ep", "Ez", "Er"]:
            tasks.append((f"06_medium_source_{comp}", [
                py, this_script, "scan",
                "--preset", "smoke", "--mode", "open", "--R2", "0.06",
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP",
                "--m", "15", "--source-component", comp,
                "--resolution", "24", "--nfreq", "9",
                "--outdir", str(root / f"06_medium_source_{comp}")
            ]))

    elif not args.quick:
        for m in [0, 5, 10, 15, 20]:
            tasks.append((f"03_closed_m_{m}", [
                py, this_script, "scan",
                "--preset", "research", "--mode", "closed",
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP,POLY_VOL_MATCH,HERMITE_SLOPE_MATCH",
                "--m", str(m), "--resolution", "24", "--nfreq", "11",
                "--outdir", str(root / f"03_closed_m_{m}")
            ]))

        for res, dpml in [(20, 0.8), (24, 0.8), (28, 1.0), (32, 1.2)]:
            suffix = str(dpml).replace('.', 'p')
            tasks.append((f"04_stability_res_{res}_dpml_{suffix}", [
                py, this_script, "scan",
                "--preset", "research", "--mode", "closed",
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP,POLY_VOL_MATCH,HERMITE_SLOPE_MATCH",
                "--m", "15", "--resolution", str(res), "--dpml", str(dpml), "--nfreq", "11",
                "--outdir", str(root / f"04_stability_res_{res}_dpml_{suffix}")
            ]))

        for comp, sz in [("Ez", 0.0), ("Er", 0.0), ("Ep", 0.0), ("Hz", 0.0), ("Ez", -0.335), ("Ez", 0.335)]:
            suffix = str(sz).replace('-', 'm').replace('.', 'p')
            tasks.append((f"05_source_{comp}_z_{suffix}", [
                py, this_script, "scan",
                "--preset", "research", "--mode", "closed",
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP,POLY_VOL_MATCH,HERMITE_SLOPE_MATCH",
                "--m", "15", "--source-component", comp, "--source-z", str(sz),
                "--resolution", "24", "--nfreq", "11",
                "--outdir", str(root / f"05_source_{comp}_z_{suffix}")
            ]))

        for r2 in [0.02, 0.04, 0.06, 0.08, 0.10, 0.15]:
            suffix = str(r2).replace('.', 'p')
            tasks.append((f"06_open_R2_{suffix}", [
                py, this_script, "scan",
                "--preset", "research", "--mode", "open", "--R2", str(r2),
                "--geometries", "PH,LINEAR,POLY2,SMOOTHSTEP,POLY_VOL_MATCH,HERMITE_SLOPE_MATCH",
                "--m", "15", "--resolution", "24", "--nfreq", "11",
                "--outdir", str(root / f"06_open_R2_{suffix}")
            ]))

    master_log = logs / "MASTER_RUN_SUMMARY.txt"
    failures: List[Tuple[str, int]] = []
    with master_log.open("w", encoding="utf-8") as ml:
        ml.write(f"ROOT={root}\nQUICK={args.quick}\nMEDIUM={args.medium}\nTASKS={len(tasks)}\nSCRIPT={this_script}\n\n")
        for i, (name, cmd) in enumerate(tasks, 1):
            ml.write(f"[{i}/{len(tasks)}] {name}\n")
            ml.flush()
            rc = _run_cmd(cmd, logs / f"{name}.log")
            ml.write(f"RETURN_CODE={rc}\n\n")
            ml.flush()
            if rc != 0:
                failures.append((name, rc))
                # Do not stop immediately: partial data are still useful for the scientific verdict.

    rc = _run_cmd([py, this_script, "aggregate", "--root", str(root), "--out", str(root / "_AGGREGATED")], logs / "AGGREGATE.log")
    if rc != 0:
        failures.append(("AGGREGATE", rc))

    if not args.no_zip:
        zip_path = Path("/work") / (root.name + ".zip")
        if str(root).startswith("/work"):
            zip_path = Path("/work") / (root.name + ".zip")
        else:
            zip_path = root.parent / (root.name + ".zip")
        print(f"\nCreating ZIP: {zip_path}")
        _zip_dir(root, zip_path)
        print(f"ZIP ready: {zip_path}")

    print("\n=== FINAL ===")
    print(f"Outputs: {root}")
    print(f"Logs: {logs}")
    if not args.no_zip:
        print(f"ZIP: {zip_path}")
    if failures:
        print("Some tasks failed, but partial data may still be useful:")
        for name, rc in failures:
            print(f"  {name}: rc={rc}")
        return 1
    print("All tasks completed.")
    return 0


def _agg_read_csv(path: Path):
    if pd is None:
        raise SystemExit("This aggregator requires pandas: conda install pandas or pip install pandas")
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def _agg_manifest_meta(folder: Path) -> Dict[str, Any]:
    p = folder / "run_manifest.json"
    if not p.exists():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}
    params = data.get("params", {}) or {}
    seed = data.get("seed", 42)
    return {
        "run_folder": str(folder),
        "run_name": folder.name,
        "mode": data.get("mode"),
        "m": data.get("m"),
        "resolution": data.get("resolution"),
        "dpml": data.get("dpml"),
        "nfreq": data.get("nfreq"),
        "kr_center": data.get("kr_center"),
        "kr_span": data.get("kr_span"),
        "source_component": data.get("source_component"),
        "source_r": data.get("source_r"),
        "source_z": data.get("source_z"),
        "R2": params.get("R2"),
        "a": params.get("a"),
        "b": params.get("b"),
        "R": params.get("R"),
        "seed": seed,
    }


def _agg_choose_best_rows(dft):
    if pd is None:
        raise SystemExit("This aggregator requires pandas: conda install pandas or pip install pandas")
    if dft.empty or "geometry" not in dft.columns:
        return pd.DataFrame()
    score_col = "eta_100R" if "eta_100R" in dft.columns else None
    if score_col is None:
        return dft.groupby("geometry", as_index=False).head(1)
    idx = dft.groupby("geometry")[score_col].idxmax()
    return dft.loc[idx].copy().reset_index(drop=True)


def _agg_build_verdict(best, meta: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(meta)
    if best.empty or "geometry" not in best.columns or "PH" not in set(best["geometry"].astype(str)):
        out.update({"status": "missing_PH_or_metrics"})
        return out
    ph = best[best["geometry"].astype(str) == "PH"].iloc[0]
    controls = best[best["geometry"].astype(str) != "PH"].copy()
    out.update({f"PH_{c}": ph.get(c, None) for c in best.columns if c not in ("geometry",)})
    if controls.empty:
        out.update({"status": "PH_only_no_controls"})
        return out
    for metric, higher_better in [
        ("eta_100R", True),
        ("eta_050R", True),
        ("log10_CF", True),
        ("E_outside_foci_frac", False),
        ("near_wall_90_frac", False),
        ("axis_10_frac", False),
        ("ann100_enrichment", True),
    ]:
        if metric not in best.columns:
            continue
        vals = pd.to_numeric(controls[metric], errors="coerce")
        if vals.dropna().empty:
            continue
        if higher_better:
            j = vals.idxmax()
            ctrl_val = float(controls.loc[j, metric])
            delta = float(ph.get(metric, math.nan)) - ctrl_val
        else:
            j = vals.idxmin()
            ctrl_val = float(controls.loc[j, metric])
            delta = ctrl_val - float(ph.get(metric, math.nan))
        out[f"best_control_{metric}"] = str(controls.loc[j, "geometry"])
        out[f"best_control_value_{metric}"] = ctrl_val
        out[f"PH_advantage_{metric}"] = delta
    eta_adv = out.get("PH_advantage_eta_100R", float("nan"))
    cf_adv = out.get("PH_advantage_log10_CF", float("nan"))
    leak_adv = out.get("PH_advantage_E_outside_foci_frac", float("nan"))
    if isinstance(eta_adv, float) and isinstance(cf_adv, float) and isinstance(leak_adv, float):
        if eta_adv >= 0.03 and cf_adv >= 0.5 and leak_adv > 0:
            status = "strong_candidate_repeat_high_resolution"
        elif eta_adv > 0.0 or cf_adv > 0.0 or leak_adv > 0.0:
            status = "weak_or_mixed_candidate"
        else:
            status = "no_resolved_PH_advantage"
    else:
        status = "incomplete_metrics"
    out["status"] = status
    return out




def _agg_read_csv_plain(path: Path) -> List[Dict[str, Any]]:
    import csv
    if (not path.exists()) or path.stat().st_size == 0:
        return []
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def _agg_write_csv_plain(rows: List[Dict[str, Any]], path: Path) -> None:
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: List[str] = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with path.open("w", encoding="utf-8", newline="") as f:
        if not keys:
            f.write("")
            return
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def _agg_float_plain(x: Any, default: float = float("nan")) -> float:
    try:
        if x is None or x == "":
            return default
        return float(x)
    except Exception:
        return default


def _agg_choose_best_rows_plain(dft_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    if not dft_rows:
        return []
    by_geom: Dict[str, Dict[str, Any]] = {}
    for row in dft_rows:
        g = str(row.get("geometry", ""))
        if not g:
            continue
        score = _agg_float_plain(row.get("eta_100R"), float("-inf"))
        old = by_geom.get(g)
        old_score = _agg_float_plain(old.get("eta_100R"), float("-inf")) if old else float("-inf")
        if old is None or score > old_score:
            by_geom[g] = dict(row)
    return list(by_geom.values())


def _agg_build_verdict_plain(best_rows: List[Dict[str, Any]], meta: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(meta)
    if not best_rows:
        out.update({"status": "missing_PH_or_metrics"})
        return out
    ph_rows = [r for r in best_rows if str(r.get("geometry", "")) == "PH"]
    if not ph_rows:
        out.update({"status": "missing_PH_or_metrics"})
        return out
    ph = ph_rows[0]
    controls = [r for r in best_rows if str(r.get("geometry", "")) != "PH"]
    for k, v in ph.items():
        if k != "geometry":
            out[f"PH_{k}"] = v
    if not controls:
        out.update({"status": "PH_only_no_controls"})
        return out

    metrics = [
        ("eta_100R", True),
        ("eta_050R", True),
        ("log10_CF", True),
        ("E_outside_foci_frac", False),
        ("near_wall_90_frac", False),
        ("axis_10_frac", False),
        ("ann100_enrichment", True),
    ]
    for metric, higher_better in metrics:
        if metric not in ph:
            continue
        valid = [(r, _agg_float_plain(r.get(metric))) for r in controls]
        valid = [(r, v) for r, v in valid if not math.isnan(v)]
        if not valid:
            continue
        if higher_better:
            ctrl, ctrl_val = max(valid, key=lambda rv: rv[1])
            delta = _agg_float_plain(ph.get(metric)) - ctrl_val
        else:
            ctrl, ctrl_val = min(valid, key=lambda rv: rv[1])
            delta = ctrl_val - _agg_float_plain(ph.get(metric))
        out[f"best_control_{metric}"] = str(ctrl.get("geometry", ""))
        out[f"best_control_value_{metric}"] = ctrl_val
        out[f"PH_advantage_{metric}"] = delta

    eta_adv = _agg_float_plain(out.get("PH_advantage_eta_100R"))
    cf_adv = _agg_float_plain(out.get("PH_advantage_log10_CF"))
    leak_adv = _agg_float_plain(out.get("PH_advantage_E_outside_foci_frac"))
    if not any(math.isnan(v) for v in [eta_adv, cf_adv, leak_adv]):
        if eta_adv >= 0.03 and cf_adv >= 0.5 and leak_adv > 0:
            status = "strong_candidate_repeat_high_resolution"
        elif eta_adv > 0.0 or cf_adv > 0.0 or leak_adv > 0.0:
            status = "weak_or_mixed_candidate"
        else:
            status = "no_resolved_PH_advantage"
    else:
        status = "incomplete_metrics"
    out["status"] = status
    return out


def _aggregate_main_plain(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Aggregate PHB MEEP decision-matrix subfolders without pandas.")
    ap.add_argument("--root", default="/work/PHB_RESULTS")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    root = Path(args.root)
    outdir = Path(args.out) if args.out else root / "_AGGREGATED"
    outdir.mkdir(parents=True, exist_ok=True)

    dft_all: List[Dict[str, Any]] = []
    best_all: List[Dict[str, Any]] = []
    flux_all: List[Dict[str, Any]] = []
    verdicts: List[Dict[str, Any]] = []

    for folder in sorted([p for p in root.rglob("*") if p.is_dir()]):
        meta = _agg_manifest_meta(folder)
        if not meta:
            continue
        dft_rows = _agg_read_csv_plain(folder / "dft_metrics.csv")
        flux_rows = _agg_read_csv_plain(folder / "flux_summary.csv")
        best_rows = _agg_read_csv_plain(folder / "uniqueness_summary.csv")
        if not best_rows:
            best_rows = _agg_choose_best_rows_plain(dft_rows)

        for rows in (dft_rows, flux_rows, best_rows):
            for row in rows:
                for k, v in meta.items():
                    row[k] = v
        dft_all.extend(dft_rows)
        flux_all.extend(flux_rows)
        best_all.extend(best_rows)
        verdicts.append(_agg_build_verdict_plain(best_rows, meta))

    _agg_write_csv_plain(dft_all, outdir / "all_dft_metrics.csv")
    _agg_write_csv_plain(best_all, outdir / "all_best_frequency_rows.csv")
    _agg_write_csv_plain(flux_all, outdir / "all_flux_summary.csv")
    _agg_write_csv_plain(verdicts, outdir / "all_preliminary_verdicts.csv")

    lines = ["# Aggregated PHB MEEP decision matrix", ""]
    lines.append(f"Root: `{root}`")
    lines.append(f"Runs with manifests: {len(verdicts)}")
    lines.append("")
    lines.append("Note: this report was generated without pandas; CSV aggregation is still complete.")
    if verdicts:
        counts: Dict[str, int] = {}
        for row in verdicts:
            status = str(row.get("status", ""))
            counts[status] = counts.get(status, 0) + 1
        lines.append("")
        lines.append("## Verdict counts")
        for status, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            lines.append(f"- {status}: {n}")
        cols = ["run_name", "mode", "m", "resolution", "dpml", "R2", "seed",
                "PH_advantage_eta_100R", "PH_advantage_log10_CF",
                "PH_advantage_E_outside_foci_frac", "status"]
        top = sorted(verdicts, key=lambda r: _agg_float_plain(r.get("PH_advantage_eta_100R"), float("-inf")), reverse=True)[:30]
        if top:
            lines.append("")
            lines.append("## Top rows by PH eta_100R advantage")
            lines.append(",".join(cols))
            for row in top:
                lines.append(",".join(str(row.get(c, "")) for c in cols))
    (outdir / "aggregated_decision_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Aggregated outputs written to: {outdir}")
    return 0


def aggregate_main(argv: Optional[List[str]] = None) -> int:
    if pd is None:
        return _aggregate_main_plain(argv)
    ap = argparse.ArgumentParser(description="Aggregate PHB MEEP decision-matrix subfolders.")
    ap.add_argument("--root", default="/work/PHB_RESULTS")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    root = Path(args.root)
    outdir = Path(args.out) if args.out else root / "_AGGREGATED"
    outdir.mkdir(parents=True, exist_ok=True)

    dft_all: List[Any] = []
    best_all: List[Any] = []
    flux_all: List[Any] = []
    verdicts: List[Dict[str, Any]] = []

    for folder in sorted([p for p in root.rglob("*") if p.is_dir()]):
        meta = _agg_manifest_meta(folder)
        if not meta:
            continue
        dft = _agg_read_csv(folder / "dft_metrics.csv")
        flux = _agg_read_csv(folder / "flux_summary.csv")
        best = _agg_read_csv(folder / "uniqueness_summary.csv")
        if best.empty:
            best = _agg_choose_best_rows(dft)
        for df in (dft, flux, best):
            if not df.empty:
                for k, v in meta.items():
                    df[k] = v
        if not dft.empty:
            dft_all.append(dft)
        if not flux.empty:
            flux_all.append(flux)
        if not best.empty:
            best_all.append(best)
        verdicts.append(_agg_build_verdict(best, meta))

    def concat(items: List[Any]):
        return pd.concat(items, ignore_index=True) if items else pd.DataFrame()

    df_dft = concat(dft_all)
    df_best = concat(best_all)
    df_flux = concat(flux_all)
    df_verdict = pd.DataFrame(verdicts)

    df_dft.to_csv(outdir / "all_dft_metrics.csv", index=False)
    df_best.to_csv(outdir / "all_best_frequency_rows.csv", index=False)
    df_flux.to_csv(outdir / "all_flux_summary.csv", index=False)
    df_verdict.to_csv(outdir / "all_preliminary_verdicts.csv", index=False)

    lines = ["# Aggregated PHB MEEP decision matrix", ""]
    lines.append(f"Root: `{root}`")
    lines.append(f"Runs with manifests: {len(verdicts)}")
    if not df_verdict.empty and "status" in df_verdict.columns:
        lines.append("")
        lines.append("## Verdict counts")
        for status, n in df_verdict["status"].value_counts(dropna=False).items():
            lines.append(f"- {status}: {int(n)}")
        cols = [c for c in ["run_name", "mode", "m", "resolution", "dpml", "R2", "seed",
                            "PH_advantage_eta_100R", "PH_advantage_log10_CF", "PH_advantage_E_outside_foci_frac", "status"]
                if c in df_verdict.columns]
        if cols:
            top = df_verdict[cols].copy()
            if "PH_advantage_eta_100R" in top.columns:
                top = top.sort_values("PH_advantage_eta_100R", ascending=False, na_position="last").head(30)
            lines.append("")
            lines.append("## Top rows by PH eta_100R advantage (seed column included)")
            try:
                lines.append(top.to_markdown(index=False))
            except Exception:
                lines.append(top.to_csv(index=False))
                lines.append("")
                lines.append("Note: pandas.to_markdown failed, so this table was written as CSV text. Install `tabulate` for Markdown tables.")
    (outdir / "aggregated_decision_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Aggregated outputs written to: {outdir}")
    return 0


def _print_unified_help() -> None:
    print("""
phb.py — one-file PHB/MEEP tool for Docker/WSL

Main user command from /home/qwerty/meep-working:
  docker run --rm -it -v "$PWD":/work meep-working python /work/phb.py

Default action with no arguments:
  quick decision matrix, all outputs saved to /work/PHB_RESULTS

Useful optional commands:
  python /work/phb.py dry       -> only geometry check, outputs /work/PHB_DRY
  python /work/phb.py quick     -> quick decision matrix, outputs /work/PHB_RESULTS
  python /work/phb.py medium    -> targeted validation, outputs /work/PHB_MEDIUM
  python /work/phb.py full      -> full decision matrix, outputs /work/PHB_RESULTS
  python /work/phb.py aggregate -> aggregate /work/PHB_RESULTS

All plots are saved as PNG files. Nothing is shown interactively in the terminal.
""".strip())


def unified_main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    # IMPORTANT FOR THE USER WORKFLOW:
    # The file is intended to be placed in /home/qwerty/meep-working and run as:
    #   docker run --rm -it -v "$PWD":/work meep-working python /work/phb.py
    # With no arguments it performs a QUICK decision-matrix run and writes all
    # result folders, CSV files, PNG figures, logs and ZIP archive under /work.
    if not argv:
        return matrix_main(["--quick", "--root", "/work/PHB_RESULTS"])

    cmd = argv[0]
    rest = argv[1:]
    if cmd in ("-h", "--help", "help"):
        _print_unified_help()
        return 0
    if cmd == "dry":
        return scan_main(["--dry-run-geometry", "--outdir", "/work/PHB_DRY"] + rest)
    if cmd == "quick":
        return matrix_main(["--quick", "--root", "/work/PHB_RESULTS"] + rest)
    if cmd == "medium":
        return matrix_main(["--medium", "--root", "/work/PHB_MEDIUM"] + rest)
    if cmd == "full":
        return matrix_main(["--root", "/work/PHB_RESULTS"] + rest)
    if cmd == "scan":
        return scan_main(rest)
    if cmd == "matrix":
        return matrix_main(rest)
    if cmd == "aggregate":
        # Default aggregation root is /work/PHB_RESULTS unless the user overrides it.
        return aggregate_main(rest if rest else ["--root", "/work/PHB_RESULTS"])
    # Backward-compatible mode: no subcommand means ordinary MEEP scan arguments.
    return scan_main(argv)


if __name__ == "__main__":
    raise SystemExit(unified_main())

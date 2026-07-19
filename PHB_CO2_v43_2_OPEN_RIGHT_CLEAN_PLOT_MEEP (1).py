#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHB_CO2_v43_2_OPEN_RIGHT_CLEAN_PLOT_MEEP.py

Full-vector axisymmetric MEEP/FDTD verifier for the OPEN V39 telescope-output
PHB geometry.

Purpose
-------
This script is NOT a Fox-Li scalar resonator model and is NOT a closed halfring
continuation of v47.  It transfers the corrected V39 ray-J geometry into a
2.5D full-vector Maxwell/FDTD problem:

    * two separated open hyperbolic funnels;
    * empty central distance between the mouths z=-a and z=+a;
    * V39 meaning of R2: R2 reduces the RIGHT funnel radius; v41.3 allows R2=0;
    * right-funnel edge radius: R_right = R - R2;
    * fully open right half-space above the right-funnel edge; no outer PEC stop or support;
    * primary diagnostic annular interval at z≈+a: output_r_min=R-R2 if R2>0, otherwise R; output_r_max=R+R3; air continues for r>output_r_max;
    * R3: computed left funnel diameter/continuation scale, R -> R+R3; it also contributes to the computed output annulus outer radius;
    * default historical candidate: a=b=R=1, R2=0.10, outer_radius=1.10;
    * diagnostic ring width is output_r_max - output_r_min, not necessarily R2.
    * IMPORTANT: a is the fixed equatorial axial half-gap / mouth coordinate.
      It is NOT split into a_hyp and gap and is NOT changed when R2 changes.

The script writes per-case folders and a combined *_ALL_RESULTS.zip archive in
the same folder, in the same workflow style as V39.

Critical v40.3 correction
-------------------------
The right hyperbolic funnel is a FULL reflecting funnel, not an empty region and
not a short rim.  For R2>=0, R2 sets the reduced mouth radius of the right funnel:
    R_right = R - R2.
The actual reflecting right PHB surface is
    rho_right(z) = R_right - b*sqrt((z/a)^2 - 1)
which is algebraically identical to the V39 reduced-surface law
    rho_right(z) = rho_base(z) - R2,
but its end is where rho_right becomes zero:
    z_right_end = a*sqrt(1 + ((R-R2)/b)^2).
For a=b=R=1 and R2=0.10 this gives z_right_end≈1.34536, i.e. the right funnel is a real extended mirror.
For R2=0 the right funnel remains full-radius and ends at L.
The axial mouth coordinate remains z=+a.  R2 changes only the right-funnel mouth radius, not the PHB gap.

Coordinate convention
---------------------
MEEP cylindrical 2.5D coordinates are used:
    r = pos.x  (cylindrical radius)
    z = pos.z  (resonator axis; V39 x-axis)

The material is a perfect electric conductor (mp.metal if available).  This is
appropriate as a first CO2 metallic-resonator screening model.  Later versions
can replace PEC by dispersive metal if needed.

Scientific interpretation
-------------------------
A positive result requires more than a bright near-field spot.  The key metrics
are the forward Poynting flux through the V39 output annulus and the Poynting-
derived angular distributions on M1 and optional M2 planes:
    theta50/theta80/theta90/theta95, L5/L10/L15, S_phi fraction, phase RMS.

Author: OpenAI assistant for Vladimir I. Khaustov workflow
Version: v43.1-left-mouth-grid-connected, 2026-07-13
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sys
import time
import traceback
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

try:
    import meep as mp
except Exception:  # plan/geometry stages work without MEEP
    mp = None


# -----------------------------
# Data model
# -----------------------------

@dataclass(frozen=True)
class V39OpenGeometry:
    model: str = "phb"             # phb or linear
    a: float = 1.0                  # FIXED V39 PHB parameter: mouth coordinate / equatorial half-gap; do not split into gap and a_hyp
    b: float = 1.0
    R: float = 1.0
    R1: float = 0.0                # v43.1: absolute radius of the physical input aperture in the LEFT funnel. If R1>0, the left funnel is truncated at rho=R1 and waves are injected through that aperture.
    R2: float = 0.10               # V39 right-funnel reduction; now allowed to be exactly zero
    R3: float = 0.0                # v43.1: left mouth-side radial continuation/window scale; contiguous at z=-a, no detached focal-plane screen
    R4: float = 0.0                # v43.1: reserved; must remain 0. No right PEC stop or support is created.
    outer_radius: Optional[float] = None  # legacy computational-extent hint only; it creates no material
    output_r_min: Optional[float] = None  # optional diagnostic interval inner radius; no outer screen is created
    output_r_max: Optional[float] = None  # optional diagnostic interval outer radius; air continues beyond it

    @property
    def dR(self) -> float:
        # v41.3: R2 is a signed input parameter with physical domain R2>=0.
        # Do not silently abs() negative values, because that hides command-line mistakes.
        return float(self.R2)

    @property
    def R_right(self) -> float:
        return self.R - self.dR

    @property
    def input_window_radius(self) -> float:
        """Effective absolute radius of the physical input aperture in the left funnel.

        R1 is NOT a fraction of R.  If R1>0, the left funnel is physically
        truncated at the axial coordinate where the left-wall radius equals R1.
        Thus the open input cross-section is 0<=r<=R1 at z=z_R1< -a.
        This is different from the v42.2 mistake where R1 was drawn at z=-a.
        """
        if self.R1 <= 0:
            return 0.0
        return max(0.0, min(float(self.R1), float(self.R)))

    @property
    def has_left_input_window(self) -> bool:
        return self.R1 > 0.0

    @property
    def left_mouth_open_radius(self) -> float:
        """Open radial interval for the effective left-side entrance.

        Historical mode: R1=0 keeps the full left horn down to its axis-end.
        R1 mode: the left funnel is truncated at rho=R1 and the cross-section
        0<=r<=R1 is the physical input aperture.
        """
        return self.input_window_radius if self.has_left_input_window else float(self.R)

    def abs_z_for_left_wall_radius(self, radius: float) -> float:
        """Positive |z| coordinate where the left wall has the specified radius."""
        rr = max(0.0, min(float(radius), float(self.R)))
        if self.model == "linear":
            return self.a + ((self.R - rr) / max(self.R, 1e-300)) * (self.L - self.a)
        return self.a * math.sqrt(1.0 + ((self.R - rr) / max(self.b, 1e-300)) ** 2)

    @property
    def left_input_z(self) -> float:
        """Axial coordinate of the R1 physical input aperture plane.

        If R1>0 this is the plane inside the left funnel where rho_wall=R1.
        If R1=0 it coincides with the natural left axial end z=-L.
        """
        if self.has_left_input_window:
            return -self.abs_z_for_left_wall_radius(self.input_window_radius)
        return -self.L

    @property
    def left_z_start(self) -> float:
        """Start of the actually present left reflecting funnel."""
        return self.left_input_z if self.has_left_input_window else -self.L

    @property
    def output_min(self) -> float:
        """Inner radius of the primary right-side diagnostic annulus.

        User-fixed V43.1 rule:
            if R2 > 0: output_r_min = R - R2
            if R2 = 0: output_r_min = R

        Manual --output-r-min is kept only as an explicit control override.
        """
        if self.output_r_min is not None:
            return float(self.output_r_min)
        if self.R2 > 0.0:
            return self.R - self.dR
        return self.R

    @property
    def output_max(self) -> float:
        """Outer radius of the primary right-side diagnostic annulus.

        User-fixed V43.1 rule:
            output_r_max = R + R3

        R4 is not used to compute the output-window outer radius in this rule.
        Manual --output-r-max is kept only as an explicit control override.
        --outer-radius is legacy and does not define material or the diagnostic annulus.
        """
        if self.output_r_max is not None:
            return float(self.output_r_max)
        return self.R + max(0.0, float(self.R3))

    @property
    def output_width(self) -> float:
        return self.output_max - self.output_min

    @property
    def left_extension_outer_radius(self) -> float:
        return self.R + max(0.0, float(self.R3))

    @property
    def right_extension_outer_radius(self) -> float:
        """Reserved compatibility value.

        V43.1 has a completely open right half-space above the right-funnel
        edge.  R4 has no geometric or material effect and must remain zero.
        """
        return float(self.R)

    @property
    def left_mouth_pec_outer_radius(self) -> float:
        """Outer radius of the physical left mouth-side PEC continuation."""
        return self.left_extension_outer_radius

    def left_mouth_bridge_inner_radius(self, z: float) -> float:
        """Inner PEC radius of the grid-robust left-mouth bridge.

        The PHB wall has an infinite meridional slope at z=-a.  A finite-thickness
        shell and a separate radial rim that only starts at r=R can therefore
        touch only on the ideal mathematical line (z=-a,r=R) while remaining
        disconnected on the sampled FDTD grid.  For z on the horn side of the
        mouth slab, the bridge starts at the actual wall radius; on the central-
        gap side it starts at R.  This closes the unintended annular seam without
        reducing the intended open mouth r<R.
        """
        zz = float(z)
        if zz <= -self.a:
            rw = self.wall_rho(zz)
            if math.isfinite(rw):
                return max(0.0, min(float(self.R), float(rw)))
        return float(self.R)

    @property
    def left_extension_z(self) -> float:
        return -self.a

    @property
    def right_extension_z(self) -> float:
        return +self.a

    @property
    def radial_extent_radius(self) -> float:
        vals = [self.R, self.R_right, self.output_max, self.left_extension_outer_radius]
        if self.outer_radius is not None:
            vals.append(float(self.outer_radius))
        return max(vals)

    @property
    def L(self) -> float:
        return self.a * math.sqrt(1.0 + (self.R / self.b) ** 2)

    @property
    def c_focus(self) -> float:
        return math.sqrt(self.a * self.a + self.b * self.b)

    @property
    def right_z_end(self) -> float:
        """End of the full reduced right funnel where rho_right becomes zero.

        For R2>0 the right funnel mouth radius is R_right=R-R2.
        The reflecting right surface is not truncated at rho=R_right on the
        original funnel; it is a complete reduced-radius funnel:
            rho_right(z) = (R-R2) - b*sqrt((z/a)^2 - 1).
        It ends where rho_right=0.
        """
        rr = max(0.0, min(self.R, self.R_right))
        if self.model == "linear":
            return self.a + (rr / max(self.R, 1e-300)) * (self.L - self.a)
        return self.a * math.sqrt(1.0 + (rr / max(self.b, 1e-300)) ** 2)

    def base_rho(self, abs_z: float) -> float:
        if abs_z < self.a or abs_z > self.L:
            return float("nan")
        if self.model == "linear":
            return max(0.0, self.R * (self.L - abs_z) / max(self.L - self.a, 1e-300))
        val = self.R - self.b * math.sqrt(max((abs_z / self.a) ** 2 - 1.0, 0.0))
        return max(0.0, val)

    def wall_rho(self, z: float) -> float:
        """Actual reflecting wall radius at axial coordinate z.

        Left horn: base PHB/linear profile.
        Right horn: full reduced V39 profile rho_right=rho_base-dR = (R-R2)-b*sqrt(...).
        Central gap |z|<a: no wall.
        """
        z = float(z)
        if z <= -self.a and z >= self.left_z_start:
            return self.base_rho(abs(z))
        if z >= self.a and z <= self.right_z_end:
            return max(0.0, self.base_rho(abs(z)) - self.dR)
        return float("nan")


@dataclass
class RunConfig:
    stage: str
    outroot: str
    models: List[str]
    m_list: List[int]
    a_over_lambda: float
    min_a_over_lambda: float
    resolution: int
    dpml_over_lambda: float
    dpml: Optional[float]
    wall_thickness: float
    aperture_stop_thickness: Optional[float]
    source_mode: str
    source_components: str
    nsrc: int
    seed: int
    fwidth_frac: float
    after_sources: float
    aperture_offset_cells: float
    enable_m2: bool
    m2_distance: float
    m2_capture_angle_deg: float
    far_zone_safety: float
    allow_near_m2: bool
    save_profiles: bool
    save_geometry_png: bool
    archive: bool
    skip_existing: bool
    stop_on_error: bool
    narrow_theta95_deg: float
    useful_flux_min: float
    min_output_cells: float
    allow_underresolved: bool
    min_free_gb: float
    # v41 modal / laser-resonator diagnostics
    harminv_component: str
    harminv_points: str
    modal_after_sources: float
    modal_fwidth_frac: float
    verify_frequency: float
    verify_fwidth_frac: float
    m1_subwindows: str
    # v41.5 internal ring-accumulator diagnostics
    enable_internal_diagnostics: bool
    save_internal_map: bool
    internal_r_min: Optional[float]
    internal_r_max: Optional[float]
    internal_z_min: Optional[float]
    internal_z_max: Optional[float]
    eq_r_min: Optional[float]
    eq_r_max: Optional[float]
    eq_z_min: Optional[float]
    eq_z_max: Optional[float]
    axis_r_max: float
    outer_r_min: Optional[float]
    save_density_maps: bool
    density_map_full_section: bool
    match_model_seeds: bool


# -----------------------------
# Utilities
# -----------------------------

def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def parse_int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in str(s).replace(";", ",").split(",") if x.strip()]


def parse_models(s: str) -> List[str]:
    out: List[str] = []
    aliases = {"hyperbolic": "phb", "linear_phb": "linear", "straight": "linear", "cone": "linear"}
    for part in str(s).replace(";", ",").split(","):
        x = part.strip().lower()
        if not x:
            continue
        x = aliases.get(x, x)
        if x not in ("phb", "linear"):
            raise ValueError(f"Unknown model {part!r}; use phb,linear,both")
        if x not in out:
            out.append(x)
    return out or ["phb"]


def write_json(path: Path, obj: object) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: Optional[List[str]] = None) -> None:
    ensure_dir(path.parent)
    if fieldnames is None:
        keys: List[str] = []
        for row in rows:
            for k in row.keys():
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def archive_folder(folder: Path) -> Optional[Path]:
    if not folder.exists():
        return None
    zip_base = folder.with_suffix("")
    zip_path = folder.parent / f"{folder.name}_ALL_RESULTS.zip"
    if zip_path.exists():
        zip_path.unlink()
    shutil.make_archive(str(zip_base) + "_ALL_RESULTS", "zip", root_dir=str(folder.parent), base_dir=folder.name)
    return zip_path


def free_gb(path: Path) -> float:
    usage = shutil.disk_usage(str(path if path.exists() else path.parent))
    return usage.free / (1024 ** 3)


def wavelength_from_a(a: float, a_over_lambda: float) -> float:
    return float(a) / float(a_over_lambda)

def complex_parts(x, default: float = float("nan")) -> Tuple[float, float, float]:
    """Return (real, imag, abs) for Harminv values that may be real or complex.

    MEEP/Harminv can return several attributes as complex Python numbers.
    Direct float(complex) raises TypeError, which caused v41 modal runs to fail
    after the FDTD part had already completed.  This helper keeps the raw
    information without crashing the CSV/JSON post-processing.
    """
    try:
        z = complex(x)
        return float(np.real(z)), float(np.imag(z)), float(abs(z))
    except Exception:
        return float(default), 0.0, float(default)


def real_float(x, default: float = float("nan")) -> float:
    """Safely convert a scalar that may be complex to its real part."""
    try:
        return float(np.real(x))
    except Exception:
        return float(default)


def abs_float(x, default: float = float("nan")) -> float:
    """Safely convert a scalar that may be complex to absolute value."""
    try:
        return float(abs(x))
    except Exception:
        return float(default)


def weighted_percentile(values: np.ndarray, weights: np.ndarray, q01: float) -> float:
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return float("nan")
    v = values[mask]
    w = weights[mask]
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    c = np.cumsum(w)
    target = float(q01) * c[-1]
    idx = int(np.searchsorted(c, target, side="left"))
    return float(v[min(max(idx, 0), len(v) - 1)])


def phase_coherence(phase: np.ndarray, weights: np.ndarray) -> float:
    """Weighted circular coherence |<exp(i phase)>| in [0,1]."""
    phase = np.asarray(phase, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    mask = np.isfinite(phase) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return float("nan")
    p = phase[mask]
    w = weights[mask]
    z = np.sum(w * np.exp(1j * p)) / max(np.sum(w), 1e-300)
    return float(abs(z))


def circular_phase_rms(phase: np.ndarray, weights: np.ndarray) -> float:
    phase = np.asarray(phase, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    mask = np.isfinite(phase) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return float("nan")
    p = phase[mask]
    w = weights[mask]
    z = np.sum(w * np.exp(1j * p)) / max(np.sum(w), 1e-300)
    mean_phase = math.atan2(float(np.imag(z)), float(np.real(z)))
    wrapped = np.angle(np.exp(1j * (p - mean_phase)))
    return float(math.sqrt(np.sum(w * wrapped * wrapped) / max(np.sum(w), 1e-300)))


def phase_polyfit_residual_rms(r: np.ndarray, phase: np.ndarray, weights: np.ndarray, degree: int) -> float:
    """Weighted residual RMS after removing piston/tilt/curvature from radial phase.

    degree=0 is equivalent to removing only the piston phase.
    degree=1 removes radial tilt; degree=2 also removes quadratic curvature.
    The residual is wrapped back to [-pi,pi], so the value is robust to 2*pi
    phase wraps.  This is a screening diagnostic for phase-front quality.
    """
    r = np.asarray(r, dtype=float).ravel()
    phase = np.asarray(phase, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    mask = np.isfinite(r) & np.isfinite(phase) & np.isfinite(weights) & (weights > 0)
    if np.sum(mask) < max(2, degree + 2):
        return float("nan")
    rr = r[mask]
    pp = phase[mask]
    ww = weights[mask]
    order = np.argsort(rr)
    rr, pp, ww = rr[order], pp[order], ww[order]
    pp = np.unwrap(pp)
    x = rr - float(np.average(rr, weights=ww))
    try:
        if degree <= 0:
            fit = np.full_like(pp, float(np.average(pp, weights=ww)))
        else:
            coeff = np.polyfit(x, pp, degree, w=np.sqrt(ww / max(np.max(ww), 1e-300)))
            fit = np.polyval(coeff, x)
        residual = np.angle(np.exp(1j * (pp - fit)))
        return float(math.sqrt(np.sum(ww * residual * residual) / max(np.sum(ww), 1e-300)))
    except Exception:
        return float("nan")


def add_phase_front_metrics(metrics: Dict[str, object], label: str, r_grid: np.ndarray, phase: np.ndarray, weights: np.ndarray, component_label: str) -> None:
    """Append phase-front metrics for a selected complex field component."""
    prefix = f"{label}_phase_{component_label}"
    metrics[f"{prefix}_coherence"] = phase_coherence(phase, weights)
    metrics[f"{prefix}_rms_raw_rad"] = circular_phase_rms(phase, weights)
    metrics[f"{prefix}_rms_after_piston_rad"] = phase_polyfit_residual_rms(r_grid, phase, weights, 0)
    metrics[f"{prefix}_rms_after_linear_fit_rad"] = phase_polyfit_residual_rms(r_grid, phase, weights, 1)
    metrics[f"{prefix}_rms_after_quadratic_fit_rad"] = phase_polyfit_residual_rms(r_grid, phase, weights, 2)


# -----------------------------
# Geometry validation and plots
# -----------------------------

def validate_geometry(g: V39OpenGeometry, cfg: RunConfig) -> Dict[str, object]:
    problems: List[str] = []
    warnings: List[str] = []
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)

    if g.a <= 0 or g.b <= 0 or g.R <= 0:
        problems.append("a,b,R must be positive")
    if g.R2 < 0:
        problems.append("R2 must be >= 0. v41.3 explicitly allows R2=0 but does not allow negative right-funnel reduction.")
    if g.R3 < 0:
        problems.append("R3 must be >= 0. It is a radial left-focal-plane extension length.")
    if abs(g.R4) > 1e-15:
        problems.append("R4 is reserved in V43.1 and must be 0. The right half-space is fully open; no right PEC stop/support is permitted.")
    if g.R_right <= 0:
        problems.append(f"R_right=R-R2={g.R_right:g} <= 0")
    if g.output_min < -1e-12:
        problems.append(f"output_min={g.output_min:g} must be non-negative")
    if g.output_max <= g.output_min:
        problems.append(f"output_max={g.output_max:g} must be > output_min={g.output_min:g}")
    # If the diagnostic M1 interval is decoupled from R2, it should not start
    # inside the right mirror edge unless a future script explicitly supports
    # cutting the right mirror.  The main use case is output_min >= R_right.
    if g.output_min < g.R_right - 1e-12:
        problems.append(
            f"output_min={g.output_min:g} is below right-funnel edge R_right={g.R_right:g}; "
            "this would cut into the reflecting right horn. Use output_min>=R_right."
        )

    if cfg.a_over_lambda < cfg.min_a_over_lambda:
        problems.append(f"a/lambda={cfg.a_over_lambda:g} < required {cfg.min_a_over_lambda:g}")
    output_cells = g.output_width * cfg.resolution
    if output_cells < cfg.min_output_cells and not cfg.allow_underresolved:
        problems.append(f"output annulus width={g.output_width:g} has only {output_cells:.3g} cells < {cfg.min_output_cells:g}; increase resolution")
    if g.R2 == 0:
        warnings.append("R2=0 enabled: right funnel is full-radius; output annulus starts at R unless --output-r-min overrides it.")
    if g.output_r_min is not None or g.output_r_max is not None:
        warnings.append(
            "Independent diagnostic M1 interval is enabled: R2 defines the right-funnel geometry, "
            "while --output-r-min/--output-r-max define the open aperture."
        )
    if g.R2 != 0.10 or abs(g.output_max - 1.10 * g.R) > 1e-12:
        warnings.append("Default historical candidate was R2=0.10R and output_max=1.10R; current values are allowed for optimization sweeps")
    right_horn_cells = (g.right_z_end - g.a) * cfg.resolution
    if right_horn_cells < 4.0 and not cfg.allow_underresolved:
        warnings.append(
            f"Right full/reduced funnel has {right_horn_cells:.3g} axial cells at resolution={cfg.resolution}. "
            "Increase resolution for final publication-grade conclusions."
        )
    if g.R3 > 0:
        r3_cells = g.R3 * cfg.resolution
        if r3_cells < 4.0 and not cfg.allow_underresolved:
            warnings.append(
                f"R3 left focal-plane extension has only {r3_cells:.3g} radial cells; "
                "use higher resolution for final conclusions."
            )
    if g.model == "linear":
        warnings.append("linear model is a negative control: same output-window protocol but straight conical generators")
    if cfg.enable_m2:
        D = 2.0 * g.output_width
        fraunhofer = 2.0 * D * D / max(lam, 1e-300)
        if cfg.m2_distance < cfg.far_zone_safety * fraunhofer and not cfg.allow_near_m2:
            warnings.append(
                f"M2 distance={cfg.m2_distance:g} is below safety*Fraunhofer estimate={cfg.far_zone_safety * fraunhofer:g}; "
                "M2 must be interpreted as downstream/intermediate-field, not strict far-field."
            )
    return {
        "ok": not problems,
        "problems": problems,
        "warnings": warnings,
        "lambda": lam,
        "a_over_lambda": cfg.a_over_lambda,
        "R_over_lambda": g.R / lam,
        "R2_v39_right_reduction": g.R2,
        "R3_left_focal_extension": g.R3,
        "R4_reserved_no_effect": g.R4,
        "R_right": g.R_right,
        "output_r_min": g.output_min,
        "output_r_max": g.output_max,
        "output_width": g.output_width,
        "output_width_cells": output_cells,
        "left_focal_extension_z": -g.c_focus,
        "left_R3_mouth_continuation_r_range": [g.R, g.left_extension_outer_radius] if g.R3 > 0 else [],
        "left_horn_z_range": [-g.L, -g.a],
        "a_is_constant_mouth_coordinate_and_half_gap": True,
        "central_empty_gap_z_range": [-g.a, g.a],
        "right_reduced_horn_z_range": [g.a, g.right_z_end],
        "right_reduced_horn_axial_length": g.right_z_end - g.a,
        "right_reduced_horn_axial_cells": (g.right_z_end - g.a) * cfg.resolution,
        "external_focal_ring_meridional_points": [[-g.c_focus, g.R], [g.c_focus, g.R]],
    }


def plot_geometry(out_png: Path, g: V39OpenGeometry, cfg: RunConfig, rmax: Optional[float] = None, z_m2: Optional[float] = None) -> None:
    """Full V39 geometry plot only.

    v41 deliberately avoids a separate misleading right-horn-only publication
    figure.  Physical reflecting surfaces are solid.  Dashed/dotted lines are
    reference levels, output-window boundaries, or monitors only; they are not
    walls.
    """
    if plt is None or not cfg.save_geometry_png:
        return
    ensure_dir(out_png.parent)
    z = np.linspace(-g.L, max(g.L, g.right_z_end) + 0.10, 1800)
    left_r = np.full_like(z, np.nan, dtype=float)
    right_r = np.full_like(z, np.nan, dtype=float)
    for i, zz in enumerate(z):
        if -g.L <= zz <= -g.a:
            left_r[i] = g.wall_rho(float(zz))
        if g.a <= zz <= g.right_z_end:
            right_r[i] = g.wall_rho(float(zz))

    fig, ax = plt.subplots(figsize=(11.5, 6.6), dpi=160)
    ax.plot(z, left_r, lw=3.0, label="LEFT PHB reflecting funnel: mouth radius R")
    ax.plot(z, -left_r, lw=3.0)
    ax.plot(z, right_r, lw=3.0, label="RIGHT PHB reflecting funnel: mouth radius R-R2")
    ax.plot(z, -right_r, lw=3.0)

    wall_t = effective_wall_t(cfg, wavelength_from_a(g.a, cfg.a_over_lambda))
    if g.left_mouth_pec_outer_radius > g.R:
        ax.plot([-g.a, -g.a], [g.R, g.left_mouth_pec_outer_radius], lw=4.0,
                label="connected left PEC continuation R→R+R3")
        ax.plot([-g.a, -g.a], [-g.left_mouth_pec_outer_radius, -g.R], lw=4.0)
    if g.output_min > g.R_right:
        ax.plot([g.a, g.a], [g.R_right, g.output_min], lw=4.0, label="right inner aperture diaphragm R_right→output_min")
        ax.plot([g.a, g.a], [-g.output_min, -g.R_right], lw=4.0)

    # V43.2 publication-safe convention:
    # Do NOT draw an open-air interval as a thick solid segment, because that
    # visually resembles a PEC screen.  The only right-side geometry is the
    # reflecting funnel itself.  M1 is shown solely as a thin dashed virtual
    # monitor displaced into the air by aperture_offset_cells.
    ax.scatter([g.a, g.a], [g.R_right, -g.R_right], s=80, marker="s", label="right funnel edge R-R2")

    z_m1 = g.a + cfg.aperture_offset_cells / max(cfg.resolution, 1)
    ax.plot([z_m1, z_m1], [g.output_min, g.output_max], lw=1.8, ls="--", color="tab:cyan",
            label="M1 virtual monitor in air (NO material)")
    ax.plot([z_m1, z_m1], [-g.output_max, -g.output_min], lw=1.8, ls="--", color="tab:cyan")
    ax.annotate("OPEN AIR for all r > R_right",
                xy=(g.a + 0.04 * g.R, g.output_max + 0.08 * g.R),
                xytext=(g.a + 0.18 * g.R, g.output_max + 0.28 * g.R),
                arrowprops=dict(arrowstyle="->", lw=1.2), fontsize=8.5)
    if z_m2 is not None:
        ax.axvline(z_m2, ls="--", lw=1.8, label="M2 downstream monitor")

    ax.axhline(g.R, ls=":", lw=1.2, color="0.35", label="rho=R reference / focal-ring radius, not wall")
    ax.axhline(-g.R, ls=":", lw=1.2, color="0.35")
    ax.axhline(g.R_right, ls="--", lw=1.0, color="0.35", label="rho=R-R2 right-mouth radius")
    ax.axhline(-g.R_right, ls="--", lw=1.0, color="0.35")
    ax.axhline(g.output_max, ls="--", lw=1.0, color="0.55", label="rho=output_max boundary")
    ax.axhline(-g.output_max, ls="--", lw=1.0, color="0.55")
    ax.axvspan(-g.a, g.a, alpha=0.06, label="empty central gap; a is constant")
    ax.scatter([-g.c_focus, g.c_focus], [g.R, g.R], s=50, label="external focal points/rings")
    ax.scatter([-g.c_focus, g.c_focus], [-g.R, -g.R], s=50)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x axis / MEEP z")
    ax.set_ylabel("rho / MEEP r")
    ax.grid(True, alpha=0.25)
    ax.set_title(
        f"V41 open V39 geometry: {g.model}, a=b=R={g.R:g}, R2={g.R2:g}; "
        f"right mouth R-R2={g.R_right:g}; output {g.output_min:g}…{g.output_max:g}; R3={g.R3:g}; R4={g.R4:g}\n"
        f"Two physical PHB funnels are solid; dashed/dotted lines are references or monitors, not reflecting surfaces."
    )
    ax.set_xlim(-g.L - 0.10, max(g.L, g.right_z_end) + 0.25)
    rr = rmax if rmax is not None else g.radial_extent_radius + 0.35 * g.R
    ax.set_ylim(-min(rr, g.radial_extent_radius + 0.55 * g.R), min(rr, g.radial_extent_radius + 0.55 * g.R))
    ax.legend(fontsize=7.3, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)

def write_geometry_audit(cdir: Path, g: V39OpenGeometry, cfg: RunConfig) -> None:
    """Save V39 geometry samples, especially the right reduced horn.

    This is a guard against exactly the error where the right funnel is visually
    or numerically mistaken for empty space.
    """
    rows: List[Dict[str, object]] = []
    for side, z0, z1 in (("left_present_after_R1_truncation", g.left_z_start, -g.a), ("right_reduced", g.a, g.right_z_end)):
        n = 41 if side == "right_reduced" else 101
        for z in np.linspace(z0, z1, n):
            rows.append({
                "side": side,
                "x_or_meep_z": float(z),
                "rho_wall_actual": float(g.wall_rho(float(z))),
                "rho_base_reference": float(g.base_rho(abs(float(z)))) if side == "right_reduced" else float(g.wall_rho(float(z))),
                "is_reflecting_wall": True,
                "note": "RIGHT PHB MIRROR, NOT EMPTY" if side == "right_reduced" else "left PHB mirror",
            })
    if g.has_left_input_window:
        for r in np.linspace(0.0, g.input_window_radius, 21):
            rows.append({
                "side": "left_R1_input_window_open",
                "x_or_meep_z": float(g.left_input_z),
                "rho_wall_actual": float(r),
                "rho_base_reference": float("nan"),
                "is_reflecting_wall": False,
                "note": "v43.1 open input aperture: left funnel is truncated where rho=R1; source injected from left; internal sources off",
            })
        for r in np.linspace(g.input_window_radius, max(g.input_window_radius, g.left_mouth_pec_outer_radius), 21):
            rows.append({
                "side": "left_R1_input_diaphragm_closed",
                "x_or_meep_z": float(g.left_input_z),
                "rho_wall_actual": float(r),
                "rho_base_reference": float("nan"),
                "is_reflecting_wall": True,
                "note": "v43.1 PEC rim around the R1 input aperture; not an extra optical screen",
            })

    if g.R3 > 0:
        for r in np.linspace(g.R, g.left_extension_outer_radius, 21):
            rows.append({
                "side": "left_R3_mouth_continuation",
                "x_or_meep_z": float(g.left_extension_z),
                "rho_wall_actual": float(r),
                "rho_base_reference": float(g.R),
                "is_reflecting_wall": True,
                "note": "v43.1 physical left-mouth PEC continuation R..R+R3 at z=-a; no focal-plane screen",
            })
    write_csv(cdir / "geometry_audit_v39_walls.csv", rows)
    write_json(cdir / "geometry_audit_v39_key_numbers.json", {
        "V39_rule": "R2 changes right-funnel radius; right reflecting horn rho_right(x)=rho_base(x)-R2; output is over/outside its edge",
        "a": g.a,
        "b": g.b,
        "R": g.R,
        "R1": g.R1,
        "R1_left_input_window_radius": g.input_window_radius,
        "R1_left_mouth_open_radius": g.left_mouth_open_radius,
        "R1_internal_sources_forced_off": bool(g.has_left_input_window),
        "R2": g.R2,
        "R3": g.R3,
        "R4": g.R4,
        "R_right": g.R_right,
        "outer_radius_legacy": g.outer_radius,
        "output_annulus": [g.output_min, g.output_max],
        "left_R3_mouth_continuation": {"z": g.left_extension_z, "r_min": g.R, "r_max": g.left_extension_outer_radius, "detached_from_wall": False, "grid_bridge_rule": "on horn-side slab planes PEC begins at wall_rho(z); on central-gap side it begins at R"},
        "right_side_above_funnel": {"z": g.right_extension_z, "r_min": g.R_right, "r_max": "computational_boundary", "material": "air", "PEC_stop_present": False},
        "right_reflecting_horn_x_range": [g.a, g.right_z_end],
        "right_reflecting_horn_axial_length": g.right_z_end - g.a,
        "right_reflecting_horn_axial_cells_at_resolution": (g.right_z_end - g.a) * cfg.resolution,
        "warning": "v43.1: a is fixed; R1>0 creates a physical left input aperture; the right funnel is full and reflecting; R2 may be 0; R3 defines the connected left-mouth PEC continuation and the primary annular measurement interval; the right half-space above the funnel is open air; R4 is reserved and must be 0.",
    })


# -----------------------------
# MEEP model
# -----------------------------

def require_meep() -> None:
    if mp is None:
        raise RuntimeError("MEEP is not available. Run fast/confirm/full inside your meep-working Docker image.")


def effective_wall_t(cfg: RunConfig, lam: float) -> float:
    return max(float(cfg.wall_thickness), 2.0 / max(cfg.resolution, 1), 0.035 * lam)


def aperture_stop_t(cfg: RunConfig, lam: float) -> float:
    if cfg.aperture_stop_thickness is not None:
        return max(float(cfg.aperture_stop_thickness), 1.0 / max(cfg.resolution, 1))
    return max(effective_wall_t(cfg, lam), 4.0 / max(cfg.resolution, 1))


def make_material_function(g: V39OpenGeometry, cfg: RunConfig, rmax: float):
    require_meep()
    air = mp.Medium(epsilon=1.0)
    metal = getattr(mp, "metal", mp.Medium(epsilon=1.0e9))
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    wall_t = effective_wall_t(cfg, lam)
    stop_t = aperture_stop_t(cfg, lam)

    def matfun(pos):
        r = float(pos.x)
        z = float(pos.z)
        if r < 0:
            return air

        # Left input aperture / PEC rim.
        # Historical mode (R1=0): the full left funnel is present down to z=-L;
        # only the left mouth at z=-a is open to the central gap.
        # R1 mode (R1>0): the left funnel is truncated at z=g.left_input_z,
        # where rho_wall=R1. The aperture 0<=r<=R1 is open; the annular rim
        # r>=R1 at that same plane is PEC. This is not an additional screen;
        # it is the physical cut edge of the truncated funnel.
        if g.has_left_input_window:
            # The truncated left funnel has a finite annular PEC cut edge.  Its
            # radial reach is controlled only by the left-side scale R3; a
            # right-side R4 value must never silently enlarge the left geometry.
            left_outer = max(g.input_window_radius, g.left_mouth_pec_outer_radius)
            if abs(z - g.left_input_z) <= 0.5 * stop_t and (g.input_window_radius <= r <= left_outer + wall_t):
                return metal
        else:
            # Grid-robust open-left-mouth mode.  The ideal PHB wall and the
            # R3 continuation meet at (z=-a,r=R), but the PHB slope is singular
            # there.  Starting the whole finite-thickness slab only at r=R
            # leaves an unintended annular seam on nearby FDTD planes.  Bridge
            # from the actual horn-wall radius on the horn side, and from R on
            # the central-gap side.  The intended mouth r<R remains open.
            left_outer = g.left_mouth_pec_outer_radius
            if abs(z + g.a) <= 0.5 * stop_t:
                left_inner = g.left_mouth_bridge_inner_radius(z)
                if left_inner <= r <= left_outer + wall_t:
                    return metal

        # V43.2 open-right correction: no PEC screen, stop, support or
        # diaphragm is created outside the right-funnel edge.  The interval
        # output_min..output_max is a diagnostic M1 annulus only; air continues
        # for all r>output_max up to the PML.  An optional inner diaphragm is
        # retained solely for an explicit manual output_min>R_right control.
        if abs(z - g.a) <= 0.5 * stop_t:
            if g.output_min > g.R_right + 0.25 * wall_t and (g.R_right <= r <= g.output_min):
                return metal

        # V39 reflecting horn walls.
        # This includes BOTH mirrors:
        #   left horn: rho = rho_base(|x|) for x in [-L,-a];
        #   right horn: rho = rho_base(x)-R2 for x in [+a,+x_right_end].
        # Therefore the right funnel is NOT empty; waves see a PEC shell at
        # the exact V39 reduced-horn surface.
        rw = g.wall_rho(z)
        if math.isfinite(rw) and rw >= 0.0:
            if max(0.0, rw) <= r <= max(0.0, rw) + wall_t:
                return metal
        return air

    return matfun


def component_from_name(name: str):
    require_meep()
    table = {"Er": mp.Er, "Ep": mp.Ep, "Ez": mp.Ez, "Hr": mp.Hr, "Hp": mp.Hp, "Hz": mp.Hz}
    key = name.strip()
    if key not in table:
        raise ValueError(f"Unknown component {name!r}; use Ez,Er,Ep")
    return table[key]


def sample_source_point(rng: np.random.Generator, g: V39OpenGeometry, cfg: RunConfig) -> Tuple[float, float]:
    # Active medium proxy: mostly inside the empty central gap, with radius below
    # the right funnel edge to avoid direct source-on-aperture artifacts.
    r_min = 0.05 * g.R if abs(cfg.current_m) > 0 else 0.0  # type: ignore[attr-defined]
    r_max = max(r_min + 1e-6, 0.82 * g.R_right)
    z = float(rng.uniform(-0.85 * g.a, 0.85 * g.a))
    r = float(rng.uniform(r_min, r_max))
    return r, z


def make_sources(g: V39OpenGeometry, cfg: RunConfig, m: int):
    require_meep()
    object.__setattr__(cfg, "current_m", m)  # dynamic helper for sampling only
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    fcen = float(cfg.verify_frequency) if (cfg.stage == "mode_verify" and float(cfg.verify_frequency) > 0.0) else 1.0 / lam
    if cfg.stage == "modal":
        fwidth = cfg.modal_fwidth_frac * fcen
    elif cfg.stage == "mode_verify":
        fwidth = cfg.verify_fwidth_frac * fcen
    else:
        fwidth = cfg.fwidth_frac * fcen
    comp_names = [x.strip() for x in cfg.source_components.split(",") if x.strip()]
    comps = [component_from_name(x) for x in (comp_names or ["Ez"])]
    model_offset = 0 if cfg.match_model_seeds else 1000 * (0 if g.model == "phb" else 1)
    rng = np.random.default_rng(int(cfg.seed + 10000 * m + model_offset))
    records: List[Dict[str, object]] = []

    def rec(i: int, kind: str, comp: str, r: float, z: float, amp: complex, phase: float, size_r: float = 0.0, size_z: float = 0.0) -> None:
        records.append({
            "source_index": i, "source_kind": kind, "component": comp,
            "center_r": r, "center_z": z, "size_r": size_r, "size_z": size_z,
            "amplitude_real": float(np.real(amp)), "amplitude_imag": float(np.imag(amp)),
            "phase_rad": phase, "frequency": fcen, "fwidth": fwidth,
            "m": m, "model": g.model,
        })

    if cfg.source_mode == "coherent":
        r0 = 0.08 * g.R if abs(m) > 0 else 0.0
        r1 = 0.82 * g.R_right
        center = mp.Vector3(0.5 * (r0 + r1), 0, 0.0)
        size = mp.Vector3(r1 - r0, 0, 1.65 * g.a)
        rec(0, "coherent_extended_active_gap", comp_names[0], float(center.x), 0.0, 1.0 + 0j, 0.0, float(size.x), float(size.z))
        return [mp.Source(mp.GaussianSource(fcen, fwidth=fwidth), component=comps[0], center=center, size=size)], fcen, fwidth, records

    if cfg.source_mode == "single":
        r, z = sample_source_point(rng, g, cfg)
        rec(0, "single_point_active_gap", comp_names[0], r, z, 1.0 + 0j, 0.0)
        return [mp.Source(mp.GaussianSource(fcen, fwidth=fwidth), component=comps[0], center=mp.Vector3(r, 0, z))], fcen, fwidth, records

    if cfg.source_mode != "random":
        raise ValueError("--source-mode must be random, coherent, or single")

    sources = []
    nsrc = max(1, int(cfg.nsrc))
    for i in range(nsrc):
        r, z = sample_source_point(rng, g, cfg)
        phase = float(rng.uniform(0, 2 * math.pi))
        amp = complex(math.cos(phase), math.sin(phase)) / math.sqrt(nsrc)
        comp_i = comps[i % len(comps)]
        comp_name = comp_names[i % len(comp_names)]
        rec(i, "random_phase_point_active_gap", comp_name, r, z, amp, phase)
        sources.append(mp.Source(mp.GaussianSource(fcen, fwidth=fwidth), component=comp_i, center=mp.Vector3(r, 0, z), amplitude=amp))
    return sources, fcen, fwidth, records


def estimate_cell(g: V39OpenGeometry, cfg: RunConfig) -> Tuple[float, float, float, Optional[float]]:
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    dpml = cfg.dpml if cfg.dpml is not None else cfg.dpml_over_lambda * lam
    z_m1 = g.a + cfg.aperture_offset_cells / max(cfg.resolution, 1)
    z_m2 = z_m1 + cfg.m2_distance if cfg.enable_m2 else None

    r_needed = g.radial_extent_radius + max(0.40 * g.R, 0.6 * lam, 0.25)
    if cfg.enable_m2 and z_m2 is not None:
        r0, r1 = m2_interval(g, cfg, r_needed + dpml + 1.0)
        r_needed = max(r_needed, r1 + 0.20 * g.R)
    rmax = r_needed + dpml + max(0.25 * g.R, 0.6 * lam, 0.25)
    z_pos_needed = max(g.right_z_end + 0.25 * g.R, z_m1 + 0.25 * g.R)
    if z_m2 is not None:
        z_pos_needed = max(z_pos_needed, z_m2 + 0.25 * g.R)
    z_neg_needed = max(g.L, g.c_focus if g.R3 > 0 else g.L) + 0.25 * g.R
    z_half = max(z_pos_needed, z_neg_needed) + dpml + max(0.30 * g.R, 0.6 * lam, 0.25)
    return rmax, 2.0 * z_half, z_m1, z_m2


def m2_interval(g: V39OpenGeometry, cfg: RunConfig, rmax: float) -> Tuple[float, float]:
    center = 0.5 * (g.output_min + g.output_max)
    half_m1 = 0.5 * g.output_width
    half_angle = cfg.m2_distance * math.tan(math.radians(cfg.m2_capture_angle_deg))
    half = max(half_m1, half_angle)
    return max(0.0, center - half), min(rmax, center + half)


def make_simulation(g: V39OpenGeometry, cfg: RunConfig, m: int):
    require_meep()
    rmax, zspan, z_m1, z_m2 = estimate_cell(g, cfg)
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    cfg.dpml = cfg.dpml if cfg.dpml is not None else cfg.dpml_over_lambda * lam
    sources, fcen, fwidth, source_records = make_sources(g, cfg, m)
    matfun = make_material_function(g, cfg, rmax)
    courant = min(0.5, 1.0 / (abs(m) + 0.8)) if abs(m) > 0 else 0.5
    sim = mp.Simulation(
        cell_size=mp.Vector3(rmax, 0, zspan),
        boundary_layers=[mp.PML(cfg.dpml, direction=mp.R), mp.PML(cfg.dpml, direction=mp.Z)],
        resolution=cfg.resolution,
        dimensions=mp.CYLINDRICAL,
        m=int(m),
        sources=sources,
        material_function=matfun,
        force_complex_fields=True,
        accurate_fields_near_cylorigin=True,
        Courant=courant,
    )
    return sim, rmax, zspan, z_m1, z_m2, fcen, fwidth, source_records


# -----------------------------
# Diagnostics
# -----------------------------

def add_diagnostics(sim, g: V39OpenGeometry, cfg: RunConfig, rmax: float, zspan: float, z_m1: float, z_m2: Optional[float], fcen: float):
    require_meep()
    comps = [mp.Er, mp.Ep, mp.Ez, mp.Hr, mp.Hp, mp.Hz]
    r0, r1 = g.output_min, g.output_max
    center = 0.5 * (r0 + r1)
    width = max(r1 - r0, 1.0 / max(cfg.resolution, 1))
    dft_m1 = sim.add_dft_fields(comps, fcen, 0, 1, center=mp.Vector3(center, 0, z_m1), size=mp.Vector3(width, 0, 0))
    flux = {
        "M1_output_annulus_z": sim.add_flux(
            fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(center, 0, z_m1), size=mp.Vector3(width, 0, 0), direction=mp.Z),
        )
    }
    dpml = float(cfg.dpml or 0.0)
    margin = max(0.10 * g.R, 0.25 * wavelength_from_a(g.a, cfg.a_over_lambda))
    z_right = +0.5 * zspan - dpml - margin
    z_left = -0.5 * zspan + dpml + margin
    r_boundary = rmax - dpml - margin
    if z_right > z_m1:
        flux["right_boundary_all_z"] = sim.add_flux(
            fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(0.5 * r_boundary, 0, z_right), size=mp.Vector3(r_boundary, 0, 0), direction=mp.Z),
        )
    if z_left < -g.a:
        flux["left_boundary_all_z"] = sim.add_flux(
            fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(0.5 * r_boundary, 0, z_left), size=mp.Vector3(r_boundary, 0, 0), direction=mp.Z),
        )
    if r_boundary > g.radial_extent_radius:
        flux["outer_boundary_r"] = sim.add_flux(
            fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(r_boundary, 0, 0), size=mp.Vector3(0, 0, zspan - 2.0 * dpml - 2.0 * margin), direction=mp.R),
        )

    dft_internal = None
    if cfg.enable_internal_diagnostics:
        zb = zone_bounds(g, cfg)
        ir0, ir1 = zb["internal_r_min"], zb["internal_r_max"]
        iz0, iz1 = zb["internal_z_min"], zb["internal_z_max"]
        if ir1 > ir0 and iz1 > iz0:
            dft_internal = sim.add_dft_fields(
                comps, fcen, 0, 1,
                center=mp.Vector3(0.5 * (ir0 + ir1), 0, 0.5 * (iz0 + iz1)),
                size=mp.Vector3(ir1 - ir0, 0, iz1 - iz0),
            )

    dft_m2 = None
    if cfg.enable_m2 and z_m2 is not None and z_m2 < z_right:
        mr0, mr1 = m2_interval(g, cfg, rmax)
        if mr1 > mr0:
            dft_m2 = sim.add_dft_fields(comps, fcen, 0, 1, center=mp.Vector3(0.5 * (mr0 + mr1), 0, z_m2), size=mp.Vector3(mr1 - mr0, 0, 0))
            flux["M2_downstream_z"] = sim.add_flux(
                fcen, 0, 1,
                mp.FluxRegion(center=mp.Vector3(0.5 * (mr0 + mr1), 0, z_m2), size=mp.Vector3(mr1 - mr0, 0, 0), direction=mp.Z),
            )
    return dft_m1, dft_m2, dft_internal, flux


def dft_array(sim, dft, comp) -> np.ndarray:
    arr = sim.get_dft_array(dft, comp, 0)
    return np.asarray(arr).squeeze().astype(complex).ravel()


def dft_array_raw(sim, dft, comp) -> np.ndarray:
    """Return raw DFT array without flattening for 2D internal maps."""
    arr = sim.get_dft_array(dft, comp, 0)
    return np.asarray(arr).squeeze().astype(complex)


def zone_bounds(g: V39OpenGeometry, cfg: RunConfig) -> Dict[str, float]:
    """Default diagnostic zones for the V2 ring-accumulator hypothesis.

    The internal diagnostic window is the open central gap by default,
    z in [-a,+a], and r from the axis to the largest relevant aperture radius.
    The default equatorial ring zone is the primary diagnostic radial interval
    projected through the central gap.  All bounds can be overridden from CLI.
    """
    ir0 = 0.0 if cfg.internal_r_min is None else float(cfg.internal_r_min)
    ir1 = max(g.radial_extent_radius, g.output_max, g.R + max(g.R3, 0.0)) if cfg.internal_r_max is None else float(cfg.internal_r_max)
    iz0 = -g.a if cfg.internal_z_min is None else float(cfg.internal_z_min)
    iz1 = +g.a if cfg.internal_z_max is None else float(cfg.internal_z_max)
    er0 = g.output_min if cfg.eq_r_min is None else float(cfg.eq_r_min)
    er1 = g.output_max if cfg.eq_r_max is None else float(cfg.eq_r_max)
    ez0 = iz0 if cfg.eq_z_min is None else float(cfg.eq_z_min)
    ez1 = iz1 if cfg.eq_z_max is None else float(cfg.eq_z_max)
    ar1 = float(cfg.axis_r_max)
    or0 = er1 if cfg.outer_r_min is None else float(cfg.outer_r_min)
    return {
        "internal_r_min": ir0, "internal_r_max": ir1,
        "internal_z_min": iz0, "internal_z_max": iz1,
        "eq_r_min": er0, "eq_r_max": er1,
        "eq_z_min": ez0, "eq_z_max": ez1,
        "axis_r_max": ar1, "outer_r_min": or0,
    }


def component_name_from_mp(comp) -> str:
    require_meep()
    table = {mp.Er: "Er", mp.Ep: "Ep", mp.Ez: "Ez", mp.Hr: "Hr", mp.Hp: "Hp", mp.Hz: "Hz"}
    return table.get(comp, str(comp))


def analyze_internal_ring_accumulator(sim, dft, g: V39OpenGeometry, cfg: RunConfig, label: str = "INTERNAL") -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    """Axis-safe legacy internal-map diagnostic.

    v42.1 critical correction
    -------------------------
    Older V41.5/V41.6 code assumed that a 2D MEEP DFT array was stored as
    [r,z].  In the problematic internal maps MEEP returned the array as [z,r],
    which moved an on-axis right-side maximum into a false left/R3 location.
    This function now orients every raw component to the canonical internal
    convention used by the rest of the code:

        array.shape == (nr, nz),  first index -> r,  second index -> z.

    The orientation decision is based on the requested physical map size and the
    resolution, and is written into the returned metrics.  The V42 full-map
    diagnostic below is still preferred for publication because it additionally
    separates metal/buffer/safe-vacuum masks.
    """
    require_meep()
    if dft is None:
        return {f"{label}_valid": False, f"{label}_reason": "no_internal_dft"}, []

    zb = zone_bounds(g, cfg)
    raw = {
        "Er": dft_array_raw(sim, dft, mp.Er),
        "Ep": dft_array_raw(sim, dft, mp.Ep),
        "Ez": dft_array_raw(sim, dft, mp.Ez),
        "Hr": dft_array_raw(sim, dft, mp.Hr),
        "Hp": dft_array_raw(sim, dft, mp.Hp),
        "Hz": dft_array_raw(sim, dft, mp.Hz),
    }
    if any(np.asarray(v).squeeze().ndim < 2 for v in raw.values()):
        shapes = {k: tuple(np.asarray(v).squeeze().shape) for k, v in raw.items()}
        return {f"{label}_valid": False, f"{label}_reason": f"internal_dft_not_2d_shapes_{shapes}"}, []

    nr_expected = max(2, int(round((zb["internal_r_max"] - zb["internal_r_min"]) * cfg.resolution)) + 1)
    nz_expected = max(2, int(round((zb["internal_z_max"] - zb["internal_z_min"]) * cfg.resolution)) + 1)

    oriented: Dict[str, np.ndarray] = {}
    orient_rows: List[Dict[str, object]] = []
    for name, arr0 in raw.items():
        arr, info = _v42_orient_array_with_info(arr0, nr_expected, nz_expected, name)
        oriented[name] = arr
        orient_rows.append(info)

    nr = min(v.shape[0] for v in oriented.values())
    nz = min(v.shape[1] for v in oriented.values())
    for k in oriented:
        oriented[k] = oriented[k][:nr, :nz]

    r = np.linspace(zb["internal_r_min"], zb["internal_r_max"], nr)
    z = np.linspace(zb["internal_z_min"], zb["internal_z_max"], nz)
    Rg, Zg = np.meshgrid(r, z, indexing="ij")

    Er, Ep, Ez = oriented["Er"], oriented["Ep"], oriented["Ez"]
    Hr, Hp, Hz = oriented["Hr"], oriented["Hp"], oriented["Hz"]
    E2 = np.abs(Er) ** 2 + np.abs(Ep) ** 2 + np.abs(Ez) ** 2
    H2 = np.abs(Hr) ** 2 + np.abs(Hp) ** 2 + np.abs(Hz) ** 2
    U = E2 + H2

    S_r = 0.5 * np.real(Ep * np.conj(Hz) - Ez * np.conj(Hp))
    S_phi = 0.5 * np.real(Ez * np.conj(Hr) - Er * np.conj(Hz))
    S_z = 0.5 * np.real(Er * np.conj(Hp) - Ep * np.conj(Hr))

    # Cylindrical integration weight; constants 2*pi, dr, dz cancel for fractions.
    W = np.maximum(Rg, 0.0)
    total = float(np.sum(U * W))
    if not np.isfinite(total) or total <= 0:
        return {f"{label}_valid": False, f"{label}_reason": "zero_internal_energy"}, []

    axis_mask = Rg <= zb["axis_r_max"]
    eq_mask = (Rg >= zb["eq_r_min"]) & (Rg <= zb["eq_r_max"]) & (Zg >= zb["eq_z_min"]) & (Zg <= zb["eq_z_max"])
    outer_mask = Rg >= zb["outer_r_min"]
    inner_nonaxis_mask = (Rg > zb["axis_r_max"]) & (~eq_mask) & (~outer_mask)

    def efrac(mask: np.ndarray) -> float:
        return float(np.sum(U[mask] * W[mask]) / total) if np.any(mask) else 0.0

    def zone_stats(mask: np.ndarray, name: str) -> Dict[str, object]:
        if not np.any(mask):
            return {f"{label}_{name}_valid": False}
        ww = U[mask] * W[mask]
        denom = float(np.sum(ww))
        if denom <= 0:
            return {f"{label}_{name}_valid": False}
        sr = S_r[mask]; sp = S_phi[mask]; sz = S_z[mask]
        sabs = np.sqrt(sr*sr + sp*sp + sz*sz)
        return {
            f"{label}_{name}_valid": True,
            f"{label}_{name}_energy_fraction": float(np.sum(ww) / total),
            f"{label}_{name}_r_mean": float(np.sum(Rg[mask] * ww) / denom),
            f"{label}_{name}_z_mean": float(np.sum(Zg[mask] * ww) / denom),
            f"{label}_{name}_S_r_mean": float(np.sum(sr * ww) / denom),
            f"{label}_{name}_S_phi_mean": float(np.sum(sp * ww) / denom),
            f"{label}_{name}_S_z_mean": float(np.sum(sz * ww) / denom),
            f"{label}_{name}_S_r_abs_fraction": float(np.sum(np.abs(sr) * W[mask]) / max(np.sum(sabs * W[mask]), 1e-300)),
            f"{label}_{name}_S_phi_abs_fraction": float(np.sum(np.abs(sp) * W[mask]) / max(np.sum(sabs * W[mask]), 1e-300)),
            f"{label}_{name}_S_z_abs_fraction": float(np.sum(np.abs(sz) * W[mask]) / max(np.sum(sabs * W[mask]), 1e-300)),
            f"{label}_{name}_forward_Sz_fraction": float(np.sum(np.maximum(sz,0.0) * W[mask]) / max(np.sum(np.abs(sz) * W[mask]), 1e-300)),
        }

    peak_i, peak_j = np.unravel_index(int(np.nanargmax(U)), U.shape)
    first = orient_rows[0] if orient_rows else {}
    metrics: Dict[str, object] = {
        f"{label}_valid": True,
        f"{label}_WARNING": "legacy internal-map diagnostic corrected for axis orientation in v42.1; V42 FULL_MAP safe-vacuum metrics are preferred for publication",
        f"{label}_axis_convention": "array_oriented_to_(r,z)_after_reading_MEEP_raw_array",
        f"{label}_expected_nr_from_r_range": int(nr_expected),
        f"{label}_expected_nz_from_z_range": int(nz_expected),
        f"{label}_raw_shape_first_component": str(first.get("raw_shape", "")),
        f"{label}_orientation_first_component": str(first.get("orientation", "")),
        f"{label}_transposed_first_component": bool(first.get("transposed", False)),
        f"{label}_orientation_ambiguous_first_component": bool(first.get("ambiguous", False)),
        f"{label}_nr": int(nr),
        f"{label}_nz": int(nz),
        f"{label}_r_min": zb["internal_r_min"],
        f"{label}_r_max": zb["internal_r_max"],
        f"{label}_z_min": zb["internal_z_min"],
        f"{label}_z_max": zb["internal_z_max"],
        f"{label}_eq_r_min": zb["eq_r_min"],
        f"{label}_eq_r_max": zb["eq_r_max"],
        f"{label}_eq_z_min": zb["eq_z_min"],
        f"{label}_eq_z_max": zb["eq_z_max"],
        f"{label}_energy_total_weighted": total,
        f"{label}_eta_axis": efrac(axis_mask),
        f"{label}_eta_eq": efrac(eq_mask),
        f"{label}_eta_outer": efrac(outer_mask),
        f"{label}_eta_inner_nonaxis_other": efrac(inner_nonaxis_mask),
        f"{label}_energy_peak_r": float(r[peak_i]),
        f"{label}_energy_peak_z": float(z[peak_j]),
        f"{label}_energy_peak_U": float(U[peak_i, peak_j]),
        f"{label}_energy_peak_E2": float(E2[peak_i, peak_j]),
        f"{label}_energy_peak_H2": float(H2[peak_i, peak_j]),
    }
    metrics.update(zone_stats(axis_mask, "axis"))
    metrics.update(zone_stats(eq_mask, "eq"))
    metrics.update(zone_stats(outer_mask, "outer"))

    rows: List[Dict[str, object]] = []
    step_r = 1 if cfg.save_internal_map else max(1, nr // 80)
    step_z = 1 if cfg.save_internal_map else max(1, nz // 120)
    for i in range(0, nr, step_r):
        for j in range(0, nz, step_z):
            rows.append({
                "r": float(r[i]), "z": float(z[j]),
                "array_i_r": int(i), "array_j_z": int(j),
                "E2": float(E2[i, j]), "H2": float(H2[i, j]), "U": float(U[i, j]),
                "S_r": float(S_r[i, j]), "S_phi": float(S_phi[i, j]), "S_z": float(S_z[i, j]),
                "zone_axis": bool(axis_mask[i, j]), "zone_eq": bool(eq_mask[i, j]), "zone_outer": bool(outer_mask[i, j]),
            })
    return metrics, rows

def phb_envelope_rho_for_map(g: V39OpenGeometry, z: float) -> float:
    """Return the intended air-cavity envelope radius for map masking.

    This is only for plotting/saving density maps. It never changes the MEEP
    geometry.  It keeps the full open PHB meridional interior:
      * central gap |z|<a: r<=R;
      * left horn: r<=rho_left(z);
      * right horn: r<=rho_right(z)=rho_base(z)-R2.
    """
    z = float(z)
    if abs(z) < g.a:
        return g.R
    rw = g.wall_rho(z)
    if math.isfinite(rw):
        return max(0.0, rw)
    return 0.0


def phb_inside_mask_for_map(g: V39OpenGeometry, Rg: np.ndarray, Zg: np.ndarray, include_extensions: bool = True) -> np.ndarray:
    rho = np.vectorize(lambda zz: phb_envelope_rho_for_map(g, float(zz)))(Zg)
    mask = (Rg <= rho + 1e-12)
    if include_extensions:
        tol = max(2.5 / max(1, int(getattr(g, '_resolution_for_map', 96))), 0.025)
        if g.R3 > 0:
            mask = mask | ((np.abs(Zg + g.a) <= tol) & (Rg >= g.R) & (Rg <= g.left_extension_outer_radius + 1e-12))
    return mask


def save_density_map_products(cdir: Path, g: V39OpenGeometry, cfg: RunConfig, rows: List[Dict[str, object]], label: str = "INTERNAL") -> None:
    """Save publication/debug density-map products directly from the script.

    Outputs:
      density_raw_window_rz.png/csv                -- full DFT diagnostic window, r-z
      density_inside_PHB_rz.png/csv                -- only inside PHB cavity/envelope
      density_full_2D_signed_radius_inside_PHB.png -- full meridional section −r...+r vs z
      density_exit_outside_region_rz.png/csv       -- right/output-side portion of saved map

    White/blank zones in masked PNGs mean either outside PHB or outside saved
    DFT data. This function does not invent data outside the DFT window.
    """
    if plt is None or not rows:
        return
    ensure_dir(cdir)
    df = None
    try:
        import pandas as pd  # type: ignore
        df = pd.DataFrame(rows)
    except Exception:
        return
    if df.empty or not {"r", "z", "U"}.issubset(df.columns):
        return

    r_vals = np.array(sorted(df["r"].unique()), dtype=float)
    z_vals = np.array(sorted(df["z"].unique()), dtype=float)
    nr, nz = len(r_vals), len(z_vals)
    if nr < 2 or nz < 2:
        return
    # Pivot all important scalar fields.
    piv = df.pivot(index="z", columns="r", values="U").sort_index().sort_index(axis=1)
    U = piv.values.astype(float)  # shape nz x nr
    Rg, Zg = np.meshgrid(r_vals, z_vals, indexing="xy")
    object.__setattr__(g, "_resolution_for_map", cfg.resolution)
    inside = phb_inside_mask_for_map(g, Rg, Zg, True)
    exit_region = (Zg >= g.a) & (Rg >= 0.0)
    output_annulus_region = (Zg >= g.a - 0.5 / max(1, cfg.resolution)) & (Rg >= g.output_min) & (Rg <= g.output_max)

    df2 = df.copy()
    # Add masks to CSV rows without assuming row ordering.
    inside_flat = []
    exit_flat = []
    outann_flat = []
    for _, row in df2.iterrows():
        rr = float(row["r"]); zz = float(row["z"])
        ins = phb_inside_mask_for_map(g, np.array([[rr]]), np.array([[zz]]), True)[0, 0]
        inside_flat.append(bool(ins))
        exit_flat.append(bool(zz >= g.a))
        outann_flat.append(bool((zz >= g.a - 0.5 / max(1, cfg.resolution)) and (g.output_min <= rr <= g.output_max)))
    df2["inside_PHB_map_mask"] = inside_flat
    df2["exit_side_z_ge_a"] = exit_flat
    df2["output_annulus_r_window"] = outann_flat
    df2.to_csv(cdir / "density_map_data_with_masks.csv", index=False)

    def draw_map(arr: np.ndarray, out_png: Path, title: str, *, full_signed: bool = False) -> None:
        A = np.asarray(arr, dtype=float)
        if full_signed:
            r_signed = np.concatenate([-r_vals[r_vals > 0][::-1], r_vals])
            A = np.hstack([A[:, r_vals > 0][:, ::-1], A])
            extent = [float(r_signed.min()), float(r_signed.max()), float(z_vals.min()), float(z_vals.max())]
            xlabel = "signed radius ρ  (−r ... +r)"
        else:
            extent = [float(r_vals.min()), float(r_vals.max()), float(z_vals.min()), float(z_vals.max())]
            xlabel = "r"
        Lg = np.log10(A + 1.0)
        finite = np.isfinite(Lg)
        if not finite.any():
            return
        fig, ax = plt.subplots(figsize=(10, 7), dpi=170)
        im = ax.imshow(
            Lg,
            origin="lower",
            extent=extent,
            aspect="auto",
            vmin=float(np.nanpercentile(Lg, 2)),
            vmax=float(np.nanpercentile(Lg, 99.5)),
        )
        # Geometry references only, not fake field.
        if full_signed:
            zz_left = np.linspace(-g.L, -g.a, 500)
            rr_left = np.array([g.wall_rho(float(zz)) for zz in zz_left])
            zz_right = np.linspace(g.a, g.right_z_end, 500)
            rr_right = np.array([g.wall_rho(float(zz)) for zz in zz_right])
            ax.plot(rr_left, zz_left, lw=2.0)
            ax.plot(-rr_left, zz_left, lw=2.0)
            ax.plot(rr_right, zz_right, lw=2.0)
            ax.plot(-rr_right, zz_right, lw=2.0)
            ax.axvline(0.0, lw=0.9)
        else:
            # draw the upper wall curve in r-z coordinates
            zz_left = np.linspace(-g.L, -g.a, 500)
            rr_left = np.array([g.wall_rho(float(zz)) for zz in zz_left])
            zz_right = np.linspace(g.a, g.right_z_end, 500)
            rr_right = np.array([g.wall_rho(float(zz)) for zz in zz_right])
            ax.plot(rr_left, zz_left, lw=2.0)
            ax.plot(rr_right, zz_right, lw=2.0)
        ax.axhline(0.0, lw=0.8)
        ax.axhline(-g.c_focus, lw=0.8)
        ax.axhline(g.c_focus, lw=0.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label="log10(U+1)")
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)

    # Raw full saved window.
    draw_map(U, cdir / "density_raw_window_rz.png", f"{g.model} m-density: raw saved DFT window, r-z")
    # Inside PHB only.
    U_inside = U.copy()
    U_inside[~inside] = np.nan
    draw_map(U_inside, cdir / "density_inside_PHB_rz.png", f"{g.model}: energy density only inside PHB/cavity mask, r-z")
    draw_map(U_inside, cdir / "density_full_2D_signed_radius_inside_PHB.png", f"{g.model}: full 2D section −r…+r, inside PHB only", full_signed=True)
    # Exit/right side and output annulus maps.
    U_exit = U.copy()
    U_exit[~exit_region] = np.nan
    draw_map(U_exit, cdir / "density_exit_side_rz.png", f"{g.model}: exit/right-side saved energy map z≥a, r-z")
    U_outann = U.copy()
    U_outann[~output_annulus_region] = np.nan
    draw_map(U_outann, cdir / "density_output_annulus_region_rz.png", f"{g.model}: output annulus/corridor map, r={g.output_min:g}…{g.output_max:g}")

    # Save compact grid arrays for reproducible plotting.
    np.savez_compressed(
        cdir / "density_map_grids.npz",
        r=r_vals,
        z=z_vals,
        U=U,
        inside_PHB_mask=inside,
        exit_side_mask=exit_region,
        output_annulus_mask=output_annulus_region,
    )


def analyze_z_profile(sim, dft, g: V39OpenGeometry, cfg: RunConfig, r_min: float, r_max: float, z_plane: float, label: str) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    require_meep()
    Er = dft_array(sim, dft, mp.Er)
    Ep = dft_array(sim, dft, mp.Ep)
    Ez = dft_array(sim, dft, mp.Ez)
    Hr = dft_array(sim, dft, mp.Hr)
    Hp = dft_array(sim, dft, mp.Hp)
    Hz = dft_array(sim, dft, mp.Hz)
    n = min(len(Er), len(Ep), len(Ez), len(Hr), len(Hp), len(Hz))
    if n <= 0:
        return {f"{label}_valid": False}, []
    Er, Ep, Ez, Hr, Hp, Hz = Er[:n], Ep[:n], Ez[:n], Hr[:n], Hp[:n], Hz[:n]
    r_grid = np.linspace(r_min, r_max, n)

    # Complex time-average Poynting components in cylindrical coordinates:
    # S = 0.5 Re(E x H*)
    S_r = 0.5 * np.real(Ep * np.conj(Hz) - Ez * np.conj(Hp))
    S_phi = 0.5 * np.real(Ez * np.conj(Hr) - Er * np.conj(Hz))
    S_z = 0.5 * np.real(Er * np.conj(Hp) - Ep * np.conj(Hr))
    E2 = np.abs(Er) ** 2 + np.abs(Ep) ** 2 + np.abs(Ez) ** 2
    H2 = np.abs(Hr) ** 2 + np.abs(Hp) ** 2 + np.abs(Hz) ** 2
    area_w = np.maximum(r_grid, 0.0)  # 2*pi and dr cancel in normalized quantities
    forward_w = np.maximum(S_z, 0.0) * area_w
    abs_w = np.sqrt(S_r * S_r + S_phi * S_phi + S_z * S_z) * area_w

    theta_poloidal = np.degrees(np.arctan2(np.abs(S_r), np.maximum(S_z, 1e-300)))
    theta_3d = np.degrees(np.arctan2(np.sqrt(S_r * S_r + S_phi * S_phi), np.maximum(S_z, 1e-300)))
    theta_poloidal = np.where(S_z > 0, theta_poloidal, 180.0)
    theta_3d = np.where(S_z > 0, theta_3d, 180.0)
    phase_Er = np.angle(Er)
    phase_Ep = np.angle(Ep)
    phase_Ez = np.angle(Ez)

    # Dominant electric component on the output ring.  Ez alone can be a poor
    # phase-front proxy near Ez nodes or for m>0 hybrid fields, so v41.5
    # records both Ez and the dominant E-component phase metrics.
    comp_power = {
        "Er": float(np.sum(np.abs(Er) ** 2 * forward_w)),
        "Ep": float(np.sum(np.abs(Ep) ** 2 * forward_w)),
        "Ez": float(np.sum(np.abs(Ez) ** 2 * forward_w)),
    }
    dominant_comp_name = max(comp_power, key=comp_power.get) if comp_power else "Ez"
    dominant_phase = {"Er": phase_Er, "Ep": phase_Ep, "Ez": phase_Ez}.get(dominant_comp_name, phase_Ez)

    total_forward = float(np.sum(forward_w))
    total_abs = float(np.sum(abs_w))
    # Local axial-tube diagnostics: useful for deciding whether the emitted
    # energy is a tubular beam parallel to the rotation axis or a radial leak.
    fw_mask = np.isfinite(forward_w) & (forward_w > 0)
    if np.any(fw_mask):
        idx_min_pol = int(np.nanargmin(np.where(fw_mask, theta_poloidal, np.nan)))
        idx_min_3d = int(np.nanargmin(np.where(fw_mask, theta_3d, np.nan)))
    else:
        idx_min_pol = idx_min_3d = -1
    sign_Sr = np.sign(S_r)
    sign_changes = int(np.sum((sign_Sr[:-1] * sign_Sr[1:] < 0) & np.isfinite(S_r[:-1]) & np.isfinite(S_r[1:]))) if len(S_r) > 1 else 0

    metrics: Dict[str, object] = {
        f"{label}_valid": True,
        f"{label}_r_min": float(r_min),
        f"{label}_r_max": float(r_max),
        f"{label}_z": float(z_plane),
        f"{label}_samples": int(n),
        f"{label}_forward_weight_sum": total_forward,
        f"{label}_abs_poynting_weight_sum": total_abs,
        f"{label}_forward_fraction_of_abs": float(total_forward / max(total_abs, 1e-300)),
        f"{label}_E2_weighted_mean": float(np.sum(E2 * area_w) / max(np.sum(area_w), 1e-300)),
        f"{label}_E2_peak": float(np.nanmax(E2)) if E2.size else float("nan"),
        f"{label}_E2_peak_r": float(r_grid[int(np.nanargmax(E2))]) if E2.size else float("nan"),
        f"{label}_dominant_E_component": dominant_comp_name,
        f"{label}_dominant_E_power_Er": comp_power.get("Er", float("nan")),
        f"{label}_dominant_E_power_Ep": comp_power.get("Ep", float("nan")),
        f"{label}_dominant_E_power_Ez": comp_power.get("Ez", float("nan")),
        # Backward-compatible raw Ez phase RMS used in v41.3.
        f"{label}_phase_Ez_rms_rad": circular_phase_rms(phase_Ez, forward_w),
        f"{label}_S_phi_abs_fraction": float(np.sum(np.abs(S_phi) * area_w) / max(np.sum((np.abs(S_r) + np.abs(S_phi) + np.abs(S_z)) * area_w), 1e-300)),
        f"{label}_S_r_abs_fraction": float(np.sum(np.abs(S_r) * area_w) / max(np.sum((np.abs(S_r) + np.abs(S_phi) + np.abs(S_z)) * area_w), 1e-300)),
        f"{label}_S_z_abs_fraction": float(np.sum(np.abs(S_z) * area_w) / max(np.sum((np.abs(S_r) + np.abs(S_phi) + np.abs(S_z)) * area_w), 1e-300)),
        f"{label}_S_r_zero_crossings": sign_changes,
        f"{label}_min_theta_poloidal_deg": float(theta_poloidal[idx_min_pol]) if idx_min_pol >= 0 else float("nan"),
        f"{label}_min_theta_poloidal_r": float(r_grid[idx_min_pol]) if idx_min_pol >= 0 else float("nan"),
        f"{label}_min_theta_3d_deg": float(theta_3d[idx_min_3d]) if idx_min_3d >= 0 else float("nan"),
        f"{label}_min_theta_3d_r": float(r_grid[idx_min_3d]) if idx_min_3d >= 0 else float("nan"),
        f"{label}_axial_parallel_fraction_abs": float(np.sum(np.abs(S_z) * area_w) / max(np.sum(np.sqrt(S_r*S_r + S_phi*S_phi + S_z*S_z) * area_w), 1e-300)),
    }
    # v41.5 phase-front diagnostics.  A raw RMS >~1 rad is a strong warning;
    # if the RMS remains >~1 rad after linear/quadratic removal, the output
    # is not a clean coherent tubular beam candidate.
    add_phase_front_metrics(metrics, label, r_grid, phase_Ez, forward_w, "Ez")
    add_phase_front_metrics(metrics, label, r_grid, dominant_phase, forward_w, "dominantE")

    for p in (50, 80, 90, 95):
        metrics[f"{label}_theta{p}_poloidal_deg"] = weighted_percentile(theta_poloidal, forward_w, p / 100.0)
        metrics[f"{label}_theta{p}_3d_deg"] = weighted_percentile(theta_3d, forward_w, p / 100.0)
    for th in (1, 5, 10, 15, 20, 25):
        metrics[f"{label}_L{th}_poloidal_fraction"] = float(np.sum(forward_w[theta_poloidal <= th]) / max(total_forward, 1e-300))
        metrics[f"{label}_L{th}_3d_fraction"] = float(np.sum(forward_w[theta_3d <= th]) / max(total_forward, 1e-300))

    # Dimensionless tube merit in [0,~1]: high only when the beam is forward,
    # 3D-narrow, and weakly azimuthal.  This is a screening metric, not a
    # replacement for theta95/L10/L15/flux tables.
    try:
        metrics[f"{label}_tube_merit_L15_3d_forward_lowSwirl"] = float(
            metrics[f"{label}_L15_3d_fraction"]
            * metrics[f"{label}_forward_fraction_of_abs"]
            * max(0.0, 1.0 - metrics[f"{label}_S_phi_abs_fraction"])
        )
    except Exception:
        metrics[f"{label}_tube_merit_L15_3d_forward_lowSwirl"] = float("nan")

    rows: List[Dict[str, object]] = []
    for i in range(n):
        rows.append({
            "i": i,
            "r": float(r_grid[i]),
            "z": float(z_plane),
            "E2": float(E2[i]),
            "H2": float(H2[i]),
            "S_r": float(S_r[i]),
            "S_phi": float(S_phi[i]),
            "S_z": float(S_z[i]),
            "theta_poloidal_deg": float(theta_poloidal[i]),
            "theta_3d_deg": float(theta_3d[i]),
            "phase_Er_rad": float(phase_Er[i]),
            "phase_Ep_rad": float(phase_Ep[i]),
            "phase_Ez_rad": float(phase_Ez[i]),
            "forward_weight": float(forward_w[i]),
        })
    return metrics, rows


def parse_m1_subwindows(s: str, g: V39OpenGeometry) -> List[Tuple[str, float, float]]:
    """Parse diagnostic subwindows for the M1 output annulus.

    Syntax: "name:r0:r1,name2:r0:r1" or "r0:r1,r0:r1".
    Values are absolute radii in the same dimensionless units as R.
    The full V39 annulus is always available as M1 metrics; these subwindows
    test whether the narrow poloidal core occupies only part of 0.90..1.10.
    """
    text = str(s or "").strip()
    if not text:
        # Geometry-adaptive defaults.  Fixed historical radii such as
        # 0.95..1.05 are invalid for geometries with R != 1.
        r0 = float(g.output_min)
        r1 = float(g.output_max)
        w = r1 - r0
        mid = 0.5 * (r0 + r1)
        text = (
            f"inner_half:{r0}:{mid},outer_half:{mid}:{r1},"
            f"central_50:{r0+0.25*w}:{r1-0.25*w},"
            f"inner_25:{r0}:{r0+0.25*w},outer_25:{r1-0.25*w}:{r1}"
        )
    out: List[Tuple[str, float, float]] = []
    for item in text.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        parts = [x.strip() for x in item.split(":") if x.strip()]
        if len(parts) == 2:
            r0, r1 = float(parts[0]), float(parts[1])
            name = f"sub_{str(r0).replace('.', 'p')}_{str(r1).replace('.', 'p')}"
        elif len(parts) == 3:
            name, r0, r1 = parts[0], float(parts[1]), float(parts[2])
        else:
            raise ValueError(f"Bad --m1-subwindows item {item!r}; use name:r0:r1 or r0:r1")
        r0 = max(float(g.output_min), float(r0))
        r1 = min(float(g.output_max), float(r1))
        if r1 > r0:
            safe = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)
            out.append((safe, r0, r1))
    return out


def analyze_profile_rows_subwindow(rows: List[Dict[str, object]], r0: float, r1: float, label: str) -> Dict[str, object]:
    """Compute angle/flux metrics from an already sampled M1 profile interval."""
    if not rows:
        return {f"{label}_valid": False, f"{label}_reason": "empty_parent_profile"}
    rr = np.array([float(x.get("r", float("nan"))) for x in rows], dtype=float)
    mask = np.isfinite(rr) & (rr >= r0 - 1e-12) & (rr <= r1 + 1e-12)
    if not np.any(mask):
        return {f"{label}_valid": False, f"{label}_r_min": float(r0), f"{label}_r_max": float(r1), f"{label}_samples": 0}
    theta = np.array([float(x.get("theta_poloidal_deg", x.get("theta_poloidal_outward_deg", float("nan")))) for x in rows], dtype=float)[mask]
    theta3 = np.array([float(x.get("theta_3d_deg", x.get("theta_3d_outward_deg", float("nan")))) for x in rows], dtype=float)[mask]
    w = np.array([float(x.get("forward_weight", x.get("outward_weight", 0.0))) for x in rows], dtype=float)[mask]
    E2 = np.array([float(x.get("E2", float("nan"))) for x in rows], dtype=float)[mask]
    Sr = np.array([float(x.get("S_r", 0.0)) for x in rows], dtype=float)[mask]
    Sp = np.array([float(x.get("S_phi", 0.0)) for x in rows], dtype=float)[mask]
    Sz = np.array([float(x.get("S_z", 0.0)) for x in rows], dtype=float)[mask]
    rr_m = rr[mask]
    phase_ez_all = np.array([float(x.get("phase_Ez_rad", float("nan"))) for x in rows], dtype=float)[mask]
    total = float(np.sum(w))
    parent_total = float(sum(float(x.get("forward_weight", x.get("outward_weight", 0.0))) for x in rows))
    denom = np.sum((np.abs(Sr)+np.abs(Sp)+np.abs(Sz)) * np.maximum(rr_m, 0.0))
    metrics: Dict[str, object] = {
        f"{label}_valid": True,
        f"{label}_r_min": float(r0),
        f"{label}_r_max": float(r1),
        f"{label}_samples": int(np.sum(mask)),
        f"{label}_forward_weight_sum": total,
        f"{label}_parent_M1_forward_fraction": float(total / max(parent_total, 1e-300)),
        f"{label}_E2_peak": float(np.nanmax(E2)) if E2.size else float("nan"),
        f"{label}_E2_peak_r": float(rr_m[int(np.nanargmax(E2))]) if E2.size else float("nan"),
        f"{label}_S_phi_abs_fraction": float(np.sum(np.abs(Sp) * np.maximum(rr_m, 0.0)) / max(denom, 1e-300)),
        f"{label}_S_r_abs_fraction": float(np.sum(np.abs(Sr) * np.maximum(rr_m, 0.0)) / max(denom, 1e-300)),
        f"{label}_S_z_abs_fraction": float(np.sum(np.abs(Sz) * np.maximum(rr_m, 0.0)) / max(denom, 1e-300)),
        f"{label}_phase_Ez_rms_rad": circular_phase_rms(phase_ez_all, w),
    }
    add_phase_front_metrics(metrics, label, rr_m, phase_ez_all, w, "Ez")
    for p in (50, 80, 90, 95):
        metrics[f"{label}_theta{p}_poloidal_deg"] = weighted_percentile(theta, w, p / 100.0)
        metrics[f"{label}_theta{p}_3d_deg"] = weighted_percentile(theta3, w, p / 100.0)
    for th in (1, 5, 10, 15, 20, 25):
        metrics[f"{label}_L{th}_poloidal_fraction"] = float(np.sum(w[theta <= th]) / max(total, 1e-300))
        metrics[f"{label}_L{th}_3d_fraction"] = float(np.sum(w[theta3 <= th]) / max(total, 1e-300))
    return metrics


def write_subwindow_summary(path: Path, row_metrics: Dict[str, object], subwindows: List[Tuple[str, float, float]]) -> None:
    rows: List[Dict[str, object]] = []
    for name, r0, r1 in subwindows:
        prefix = f"M1W_{name}"
        rows.append({
            "name": name,
            "r_min": r0,
            "r_max": r1,
            "valid": row_metrics.get(f"{prefix}_valid"),
            "samples": row_metrics.get(f"{prefix}_samples"),
            "parent_M1_forward_fraction": row_metrics.get(f"{prefix}_parent_M1_forward_fraction"),
            "theta50_poloidal_deg": row_metrics.get(f"{prefix}_theta50_poloidal_deg"),
            "theta90_poloidal_deg": row_metrics.get(f"{prefix}_theta90_poloidal_deg"),
            "theta95_poloidal_deg": row_metrics.get(f"{prefix}_theta95_poloidal_deg"),
            "L10_poloidal_fraction": row_metrics.get(f"{prefix}_L10_poloidal_fraction"),
            "L15_poloidal_fraction": row_metrics.get(f"{prefix}_L15_poloidal_fraction"),
            "theta95_3d_deg": row_metrics.get(f"{prefix}_theta95_3d_deg"),
            "S_phi_abs_fraction": row_metrics.get(f"{prefix}_S_phi_abs_fraction"),
            "phase_Ez_rms_raw_rad": row_metrics.get(f"{prefix}_phase_Ez_rms_raw_rad"),
            "phase_Ez_rms_after_linear_fit_rad": row_metrics.get(f"{prefix}_phase_Ez_rms_after_linear_fit_rad"),
            "phase_Ez_rms_after_quadratic_fit_rad": row_metrics.get(f"{prefix}_phase_Ez_rms_after_quadratic_fit_rad"),
            "E2_peak_r": row_metrics.get(f"{prefix}_E2_peak_r"),
        })
    write_csv(path, rows)


def plot_subwindow_bars(csv_path: Path, png_path: Path, title: str) -> None:
    if plt is None or not csv_path.exists():
        return
    rows = []
    with csv_path.open("r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    if not rows:
        return
    labels = [r["name"] for r in rows]
    vals = [float(r.get("theta95_poloidal_deg") or "nan") for r in rows]
    l15 = [float(r.get("L15_poloidal_fraction") or "nan") for r in rows]
    fig, ax1 = plt.subplots(figsize=(max(8, 0.9*len(labels)), 4.8), dpi=150)
    x = np.arange(len(labels))
    ax1.bar(x, vals, alpha=0.7, label="theta95 poloidal")
    ax1.axhline(25, ls=":", lw=1)
    ax1.set_ylabel("theta95, deg")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=35, ha="right")
    ax1.grid(True, axis="y", alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(x, l15, marker="o", label="L15")
    ax2.set_ylabel("L15 fraction")
    ax2.set_ylim(0, 1.05)
    ax1.set_title(title)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def plot_profile(profile_csv: Path, png_path: Path, title: str) -> None:
    if plt is None or not profile_csv.exists():
        return
    rows = []
    with profile_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    if not rows:
        return
    r = np.array([float(x["r"]) for x in rows])
    theta = np.array([float(x["theta_poloidal_deg"]) for x in rows])
    theta3 = np.array([float(x["theta_3d_deg"]) for x in rows])
    Sz = np.array([float(x["S_z"]) for x in rows])
    E2 = np.array([float(x["E2"]) for x in rows])
    fig, ax1 = plt.subplots(figsize=(9, 5), dpi=150)
    ax1.plot(r, np.clip(theta, 0, 80), label="theta poloidal")
    ax1.plot(r, np.clip(theta3, 0, 80), label="theta 3D")
    ax1.set_xlabel("r")
    ax1.set_ylabel("angle, deg")
    ax1.grid(True, alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(r, Sz, alpha=0.45, label="S_z")
    ax2.plot(r, E2 / max(float(np.nanmax(E2)), 1e-300), alpha=0.45, label="E2 norm")
    ax2.set_ylabel("S_z / normalized E2")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, fontsize=8, loc="upper right")
    ax1.set_title(title)
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def classify(metrics: Dict[str, object], cfg: RunConfig) -> str:
    theta95 = float(metrics.get("M1_theta95_poloidal_deg", float("nan")))
    theta95_3d = float(metrics.get("M1_theta95_3d_deg", float("nan")))
    flux = abs(float(metrics.get("flux_M1_output_annulus_z", 0.0)))
    l15 = float(metrics.get("M1_L15_poloidal_fraction", 0.0))
    l15_3d = float(metrics.get("M1_L15_3d_fraction", 0.0))
    sphi = float(metrics.get("M1_S_phi_abs_fraction", 1.0))
    phase_raw = float(metrics.get("M1_phase_Ez_rms_raw_rad", metrics.get("M1_phase_Ez_rms_rad", float("nan"))))
    phase_quad = float(metrics.get("M1_phase_Ez_rms_after_quadratic_fit_rad", phase_raw))
    phase_bad = (math.isfinite(phase_quad) and phase_quad > 1.0) or (math.isfinite(phase_raw) and phase_raw > 1.5)
    if not math.isfinite(theta95):
        return "NO_PROFILE"
    if flux < cfg.useful_flux_min:
        return "WEAK_FLUX"
    if theta95 <= cfg.narrow_theta95_deg and l15 >= 0.85 and math.isfinite(theta95_3d) and theta95_3d <= cfg.narrow_theta95_deg and l15_3d >= 0.70 and sphi <= 0.20:
        return "PHASE_DISORDERED_NARROW_ANGLES" if phase_bad else "PROMISING_TUBULAR_AXIAL_M1_3D"
    if theta95 <= cfg.narrow_theta95_deg and l15 >= 0.85:
        return "PHASE_DISORDERED_POLOIDAL_LAYER" if phase_bad else "PROMISING_POLOIDAL_LAYER_WITH_3D_CHECK_NEEDED"
    if theta95 <= 35.0 and l15 >= 0.60:
        return "PARTIAL_CORE_WITH_TAIL"
    return "BROAD_OR_UNORDERED"


def classify_v42(metrics: Dict[str, object], cfg: UniversalRunConfig) -> str:
    """Classify the dedicated physical M1 output-annulus result.

    V43.1 deliberately uses the signed outward flux through the real annulus;
    a negative/backward flux is never converted to a positive useful output by
    taking an absolute value. Angular metrics are secondary to the existence
    of a resolved forward flux.
    """
    def mf(name: str, default: float = float("nan")) -> float:
        try:
            return float(metrics.get(name, default))
        except Exception:
            return float(default)

    if not bool(metrics.get("M1_valid", metrics.get("M1_OUTPUT_ANNULUS_valid", False))):
        return "NO_M1_OUTPUT_ANNULUS_PROFILE"

    flux = mf("flux_M1_output_annulus_z", mf("flux_M1_OUTPUT_ANNULUS_outward_signed", 0.0))
    positive_proxy = mf("M1_positive_outward_flux_proxy", 0.0)
    theta95 = mf("M1_theta95_poloidal_deg")
    theta95_3d = mf("M1_theta95_3d_deg")
    l15 = mf("M1_L15_poloidal_fraction", 0.0)
    l15_3d = mf("M1_L15_3d_fraction", 0.0)
    sphi = mf("M1_S_phi_abs_fraction", 1.0)
    phase_raw = mf("M1_phase_Ez_rms_raw_rad")
    phase_quad = mf("M1_phase_Ez_rms_after_quadratic_fit_rad", phase_raw)
    output_cells = mf("validation_output_width_cells", float("nan"))

    if flux <= cfg.useful_flux_min or positive_proxy <= cfg.useful_flux_min:
        return "WEAK_OR_REVERSED_M1_FLUX"
    if not math.isfinite(theta95):
        return "M1_FORWARD_FLUX_WITHOUT_VALID_ANGLE_PROFILE"

    phase_bad = ((math.isfinite(phase_quad) and phase_quad > 1.0) or
                 (math.isfinite(phase_raw) and phase_raw > 1.5))
    underresolved = math.isfinite(output_cells) and output_cells < cfg.min_output_cells

    if (theta95 <= cfg.narrow_theta95_deg and l15 >= 0.85 and
            math.isfinite(theta95_3d) and theta95_3d <= cfg.narrow_theta95_deg and
            l15_3d >= 0.70 and sphi <= 0.20):
        if underresolved:
            return "PROMISING_M1_3D_BUT_OUTPUT_UNDERRESOLVED"
        return "PHASE_DISORDERED_NARROW_M1_3D" if phase_bad else "PROMISING_M1_3D_NARROW_LOW_SWIRL"
    if theta95 <= cfg.narrow_theta95_deg and l15 >= 0.85:
        if underresolved:
            return "PROMISING_M1_POLOIDAL_BUT_OUTPUT_UNDERRESOLVED"
        return "PHASE_DISORDERED_M1_POLOIDAL_LAYER" if phase_bad else "PROMISING_M1_POLOIDAL_LAYER_3D_CHECK_NEEDED"
    if theta95 <= 35.0 and l15 >= 0.60:
        return "M1_PARTIAL_NARROW_CORE_WITH_TAIL"
    return "M1_BROAD_OR_UNORDERED"


# -----------------------------
# v41 modal / ring-down diagnostics
# -----------------------------

def parse_harminv_points(s: str, g: V39OpenGeometry) -> List[Tuple[float, float]]:
    """Parse points as 'r,z;r,z;...' in MEEP cylindrical coordinates.

    Default points are inside the active central gap and below the right-mouth
    radius, so they test resonator ring-down rather than the output monitor.
    """
    text = str(s).strip()
    if not text:
        # Multi-plane defaults reduce the risk of missing a mode whose field
        # has a node at z=0.  All probes remain inside the empty central gap.
        radii = (0.20 * g.R_right, 0.50 * g.R_right, 0.80 * g.R_right)
        zplanes = (-0.35 * g.a, 0.0, +0.35 * g.a)
        return [(r, z) for z in zplanes for r in radii]
    pts: List[Tuple[float, float]] = []
    for part in text.replace("|", ";").split(";"):
        part = part.strip()
        if not part:
            continue
        vals = [float(x.strip()) for x in part.replace(":", ",").split(",") if x.strip()]
        if len(vals) != 2:
            raise ValueError(f"Bad --harminv-points item {part!r}; use r,z;r,z")
        pts.append((vals[0], vals[1]))
    return pts


def cluster_modal_rows(rows: List[Dict[str, object]], fcen: float) -> List[Dict[str, object]]:
    """Cluster reliable Harminv detections across points and components.

    Frequency coincidence at several spatial probes and/or field components is
    treated as stronger evidence than a single high-amplitude detection.  The
    ranking penalizes large Harminv errors and broad within-cluster frequency
    spans.  It remains a screening rank, not a substitute for mode_verify.
    """
    finite: List[Dict[str, object]] = []
    for r in rows:
        try:
            f = float(r.get('frequency', float('nan')))
            q = float(r.get('Q', float('nan')))
            a = max(float(r.get('amp_abs', 0.0)), 0.0)
            e = float(r.get('err_abs', float('nan')))
            if math.isfinite(f) and f > 0 and math.isfinite(q) and q > 0:
                finite.append({**r, 'frequency': f, 'Q': q, 'amp_abs': a, 'err_abs': e})
        except Exception:
            pass
    if not finite:
        return []

    finite.sort(key=lambda x: float(x['frequency']))
    tol = max(0.003, 0.0015 * max(float(fcen), 1.0))
    clusters: List[List[Dict[str, object]]] = []
    for row in finite:
        if not clusters:
            clusters.append([row])
            continue
        center = float(np.mean([float(x['frequency']) for x in clusters[-1]]))
        if abs(float(row['frequency']) - center) > tol:
            clusters.append([row])
        else:
            clusters[-1].append(row)

    out: List[Dict[str, object]] = []
    for i, cl in enumerate(clusters):
        freqs = np.array([float(x['frequency']) for x in cl], dtype=float)
        amps = np.array([max(float(x.get('amp_abs', 0.0)), 0.0) for x in cl], dtype=float)
        qs = np.array([max(float(x.get('Q', 0.0)), 0.0) for x in cl], dtype=float)
        errs = np.array([float(x.get('err_abs', float('nan'))) for x in cl], dtype=float)
        spatial_probes = sorted(set(int(x.get('monitor_index', -1)) for x in cl))
        components = sorted(set(str(x.get('component', '')) for x in cl if str(x.get('component', ''))))
        detectors = sorted(set((str(x.get('component', '')), int(x.get('monitor_index', -1))) for x in cl))
        amp_sum = float(np.sum(amps))
        freq_w = float(np.sum(freqs * amps) / max(amp_sum, 1e-300)) if amp_sum > 0 else float(np.mean(freqs))
        q_amp_w = float(np.sum(qs * amps) / max(amp_sum, 1e-300)) if amp_sum > 0 else float(np.mean(qs))
        span = float(np.max(freqs) - np.min(freqs))
        finite_errs = errs[np.isfinite(errs)]
        err_median = float(np.median(finite_errs)) if finite_errs.size else float('nan')
        err_quality = 1.0 / (1.0 + max(err_median, 0.0)) if math.isfinite(err_median) else 0.5
        span_quality = 1.0 / (1.0 + span / max(tol, 1e-300))
        detector_count = len([d for d in detectors if d[1] >= 0])
        score = float(amp_sum * max(q_amp_w, 0.0) * math.sqrt(max(detector_count, 1)) * err_quality * span_quality)
        out.append({
            'cluster_index': i,
            'frequency_mean': float(np.mean(freqs)),
            'frequency_amp_weighted': freq_w,
            'frequency_min': float(np.min(freqs)),
            'frequency_max': float(np.max(freqs)),
            'frequency_span': span,
            'cluster_tolerance': tol,
            'occurrence_count': int(len(cl)),
            'spatial_probe_count': int(len([p for p in spatial_probes if p >= 0])),
            'probe_count': int(len([p for p in spatial_probes if p >= 0])),
            'probe_indices': ';'.join(str(p) for p in spatial_probes),
            'component_count': int(len(components)),
            'components': ';'.join(components),
            'detector_count_component_x_probe': int(detector_count),
            'Q_mean': float(np.mean(qs)),
            'Q_amp_weighted': q_amp_w,
            'Q_max': float(np.max(qs)),
            'amp_sum': amp_sum,
            'amp_max': float(np.max(amps)) if amps.size else float('nan'),
            'err_min': float(np.min(finite_errs)) if finite_errs.size else float('nan'),
            'err_median': err_median,
            'err_max': float(np.max(finite_errs)) if finite_errs.size else float('nan'),
            'error_quality_factor': err_quality,
            'span_quality_factor': span_quality,
            'score_Q_amp_detector_error_span': score,
            'score_Q_amp_probe': score,
        })
    out.sort(key=lambda x: float(x.get('score_Q_amp_detector_error_span', 0.0)), reverse=True)
    for j, row in enumerate(out):
        row['rank_by_score'] = j + 1
    return out

def run_modal_ringdown(sim, g: V39OpenGeometry, cfg: RunConfig, cdir: Path,
                       fcen: float, fwidth: float, m: int) -> Dict[str, object]:
    require_meep()
    comp_names = [x.strip() for x in str(cfg.harminv_component).replace(';', ',').split(',') if x.strip()]
    if not comp_names:
        comp_names = ['Ez']
    components = [(name, component_from_name(name)) for name in comp_names]
    pts = parse_harminv_points(cfg.harminv_points, g)

    detectors: List[Tuple[object, str, int, float, float]] = []
    for cname, comp in components:
        for ip, (r, z) in enumerate(pts):
            h = mp.Harminv(comp, mp.Vector3(r, 0, z), fcen, fwidth)
            detectors.append((h, cname, ip, float(r), float(z)))

    # Broad pulse followed by passive ring-down.  Multiple components and
    # multiple z planes reduce blind spots due to field nodes or polarization.
    sim.run(mp.after_sources(*[x[0] for x in detectors]),
            until_after_sources=float(cfg.modal_after_sources))

    rows: List[Dict[str, object]] = []
    for ih, (h, cname, ip, r, z) in enumerate(detectors):
        for im, mode in enumerate(getattr(h, 'modes', [])):
            freq_raw = getattr(mode, 'freq', float('nan'))
            Q_raw = getattr(mode, 'Q', float('nan'))
            decay_raw = getattr(mode, 'decay', float('nan'))
            amp_raw = getattr(mode, 'amp', 0.0)
            err_raw = getattr(mode, 'err', float('nan'))

            freq = real_float(freq_raw)
            Q = real_float(Q_raw)
            decay = real_float(decay_raw)
            amp_real, amp_imag, amp_abs = complex_parts(amp_raw)
            err_real, err_imag, err_abs = complex_parts(err_raw)
            in_search_band = (fcen - 0.5*fwidth) <= freq <= (fcen + 0.5*fwidth) if math.isfinite(freq) else False
            reliable_basic = bool(math.isfinite(freq) and freq > 0 and math.isfinite(Q) and Q > 0 and in_search_band)

            rows.append({
                'harminv_detector_index': ih,
                'monitor_index': ip,
                'mode_index': im,
                'component': cname,
                'r': r,
                'z': z,
                'frequency': freq,
                'frequency_imag_raw': float(np.imag(complex(freq_raw))) if np.isfinite(real_float(freq_raw, 0.0)) else float('nan'),
                'search_band_min': fcen - 0.5*fwidth,
                'search_band_max': fcen + 0.5*fwidth,
                'inside_search_band': in_search_band,
                'basic_reliable_detection': reliable_basic,
                'wavelength_equiv': float(1.0 / freq) if freq and math.isfinite(freq) and freq > 0 else float('nan'),
                'a_over_lambda_equiv': float(g.a * freq) if math.isfinite(freq) else float('nan'),
                'Q': Q,
                'decay': decay,
                'amp_real': amp_real,
                'amp_imag': amp_imag,
                'amp_abs': amp_abs,
                'err_real': err_real,
                'err_imag': err_imag,
                'err_abs': err_abs,
                'err': err_abs,
                'model_m': m,
            })

    write_csv(cdir / 'modal_modes.csv', rows)
    clusters = cluster_modal_rows(rows, fcen)
    write_csv(cdir / 'modal_cluster_summary.csv', clusters)

    finite = [r for r in rows if bool(r.get('basic_reliable_detection', False))]
    for r in finite:
        err = float(r.get('err_abs', float('nan')))
        err_quality = 1.0 / (1.0 + max(err, 0.0)) if math.isfinite(err) else 0.5
        r['score_Q_amp_error'] = float(max(float(r.get('Q', 0.0)), 0.0) *
                                             max(float(r.get('amp_abs', 0.0)), 0.0) * err_quality)
    best = max(finite, key=lambda x: float(x.get('score_Q_amp_error', 0.0)), default=None)
    best_cluster = clusters[0] if clusters else None
    metrics: Dict[str, object] = {
        'modal_components': ','.join(comp_names),
        'modal_component': ','.join(comp_names),
        'modal_fcen': fcen,
        'modal_fwidth': fwidth,
        'modal_search_band_min': fcen - 0.5*fwidth,
        'modal_search_band_max': fcen + 0.5*fwidth,
        'modal_after_sources': cfg.modal_after_sources,
        'modal_spatial_point_count': len(pts),
        'modal_component_count': len(comp_names),
        'modal_detector_count': len(detectors),
        'modal_point_count': len(pts),
        'modal_mode_count': len(rows),
        'modal_basic_reliable_mode_count': len(finite),
        'modal_cluster_count': len(clusters),
    }
    if best_cluster is not None:
        metrics.update({
            'modal_best_cluster_frequency': best_cluster.get('frequency_amp_weighted'),
            'modal_best_cluster_a_over_lambda': g.a * float(best_cluster.get('frequency_amp_weighted')),
            'modal_best_cluster_Q_amp_weighted': best_cluster.get('Q_amp_weighted'),
            'modal_best_cluster_Q_max': best_cluster.get('Q_max'),
            'modal_best_cluster_amp_sum': best_cluster.get('amp_sum'),
            'modal_best_cluster_probe_count': best_cluster.get('probe_count'),
            'modal_best_cluster_component_count': best_cluster.get('component_count'),
            'modal_best_cluster_detector_count': best_cluster.get('detector_count_component_x_probe'),
            'modal_best_cluster_frequency_span': best_cluster.get('frequency_span'),
            'modal_best_cluster_err_median': best_cluster.get('err_median'),
            'modal_best_cluster_score': best_cluster.get('score_Q_amp_detector_error_span'),
        })
    if best is not None:
        metrics.update({
            'modal_best_frequency': best['frequency'],
            'modal_best_a_over_lambda': g.a * float(best['frequency']),
            'modal_best_wavelength_equiv': best['wavelength_equiv'],
            'modal_best_Q': best['Q'],
            'modal_best_decay': best['decay'],
            'modal_best_amp_abs': best['amp_abs'],
            'modal_best_err_abs': best['err_abs'],
            'modal_best_component': best['component'],
            'modal_best_monitor_index': best['monitor_index'],
            'modal_best_r': best['r'],
            'modal_best_z': best['z'],
        })
    else:
        metrics['modal_warning'] = (
            'Harminv returned no positive-Q finite mode inside the requested search band. '
            'Increase modal_after_sources, adjust fwidth, or inspect modal_modes.csv.'
        )
    write_json(cdir / 'modal_metrics.json', metrics)
    return metrics

def run_modal_v42(sim, g: V39OpenGeometry, cfg: UniversalRunConfig,
                  cdir: Path, fcen: float, fwidth: float, m: int,
                  tag: str = "") -> Dict[str, object]:
    """V42 modal/Harminv entry point.

    V42.7 called run_modal_v42(), but the function itself was missing after
    refactoring from the legacy run_modal_ringdown().  This caused the modal
    run to stop before the actual ringdown started.  The wrapper below keeps
    the V42 call path and delegates to the tested Harminv implementation.
    """
    metrics = run_modal_ringdown(sim, g, cfg, cdir, fcen, fwidth, m)
    metrics["modal_v42_entrypoint"] = "run_modal_v42_wrapper_to_run_modal_ringdown"
    if tag:
        metrics["case_tag"] = tag
    write_json(cdir / "modal_metrics.json", metrics)
    return metrics

# -----------------------------
# Run one case
# -----------------------------

def run_one_case(args, cfg: RunConfig, model: str, m: int) -> Dict[str, object]:
    g = V39OpenGeometry(
        model=model, a=args.a, b=args.b, R=args.R, R1=args.R1, R2=args.R2, R3=args.R3, R4=args.R4,
        outer_radius=args.outer_radius, output_r_min=args.output_r_min, output_r_max=args.output_r_max
    )
    outroot = Path(cfg.outroot)
    tag = (
        f"{model}_m{m}_R2_{g.R2:.4g}_R3_{g.R3:.4g}_R4_{g.R4:.4g}_"
        f"out_{g.output_min:.4g}_{g.output_max:.4g}_aol_{cfg.a_over_lambda:.4g}_res{cfg.resolution}"
    ).replace(".", "p")
    cdir = outroot / tag
    done = cdir / "metrics.json"
    if cfg.skip_existing and done.exists():
        old = json.loads(done.read_text(encoding="utf-8"))
        old["skipped_existing"] = True
        return old
    ensure_dir(cdir)

    started = time.time()
    row: Dict[str, object] = {
        "timestamp": now_iso(),
        "script": Path(__file__).name,
        "version": "v41.5-ring-accumulator-diagnostics",
        "stage": cfg.stage,
        "model": model,
        "m": m,
        "a": g.a,
        "b": g.b,
        "R": g.R,
        "R2_v39_right_reduction": g.R2,
        "R3_left_focal_extension": g.R3,
        "R4_reserved_no_effect": g.R4,
        "R_right": g.R_right,
        "outer_radius_legacy": g.outer_radius,
        "output_r_min": g.output_min,
        "output_r_max": g.output_max,
        "output_width": g.output_width,
        "L": g.L,
        "c_focus": g.c_focus,
        "right_z_end": g.right_z_end,
        "a_over_lambda": cfg.a_over_lambda,
        "resolution": cfg.resolution,
        "source_mode": cfg.source_mode,
        "nsrc": cfg.nsrc,
    }
    validation = validate_geometry(g, cfg)
    row.update({f"validation_{k}": json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in validation.items()})
    write_json(cdir / "config.json", {"geometry": asdict(g), "config": asdict(cfg), "validation": validation})

    # Geometry plot can be written before MEEP simulation.
    try:
        est_rmax, est_zspan, est_m1, est_m2 = estimate_cell(g, cfg)
        plot_geometry(cdir / "geometry_FULL_OPEN_V39_const_a.png", g, cfg, est_rmax, est_m2)
        write_geometry_audit(cdir, g, cfg)
    except Exception as e:
        row["geometry_plot_error"] = repr(e)

    if not validation["ok"]:
        row["candidate_status"] = "INVALID_GEOMETRY"
        row["error"] = "; ".join(validation.get("problems", []))
        write_json(cdir / "metrics.json", row)
        return row

    if cfg.stage in ("plan", "geometry"):
        row["candidate_status"] = "GEOMETRY_ONLY"
        row["elapsed_seconds"] = time.time() - started
        write_json(cdir / "metrics.json", row)
        return row

    require_meep()
    if free_gb(outroot) < cfg.min_free_gb:
        raise RuntimeError(f"Not enough free disk space near {outroot}: {free_gb(outroot):.2f} GB < {cfg.min_free_gb} GB")

    sim = None
    try:
        sim, rmax, zspan, z_m1, z_m2, fcen, fwidth, source_records = make_simulation(g, cfg, m)
        row.update({"rmax": rmax, "zspan": zspan, "z_M1": z_m1, "z_M2": z_m2 if z_m2 is not None else "", "fcen": fcen, "fwidth": fwidth})
        write_csv(cdir / "sources.csv", source_records)
        print(f"RUN {tag}: rmax={rmax:.4g}, zspan={zspan:.4g}, z_M1={z_m1:.4g}, z_M2={z_m2}, fcen={fcen:.4g}, fwidth={fwidth:.4g}", flush=True)

        if cfg.stage == "modal":
            modal_metrics = run_modal_ringdown(sim, g, cfg, cdir, fcen, fwidth, m)
            row.update(modal_metrics)
            row["candidate_status"] = "MODAL_RINGDOWN_DONE" if int(modal_metrics.get("modal_mode_count", 0)) > 0 else "MODAL_NO_MODES_FOUND"
            row["elapsed_seconds"] = time.time() - started
            write_json(cdir / "metrics.json", row)
            return row

        dft_m1, dft_m2, dft_internal, flux_objs = add_diagnostics(sim, g, cfg, rmax, zspan, z_m1, z_m2, fcen)
        sim.run(until_after_sources=float(cfg.after_sources))

        for name, obj in flux_objs.items():
            try:
                row[f"flux_{name}"] = float(mp.get_fluxes(obj)[0])
            except Exception as e:
                row[f"flux_{name}_error"] = repr(e)
        # v41 energy-balance proxy: separate useful output annulus from other monitored losses.
        monitored_abs_loss = 0.0
        for k, v in list(row.items()):
            if str(k).startswith("flux_") and not str(k).endswith("_error"):
                try:
                    monitored_abs_loss += abs(float(v))
                except Exception:
                    pass
        try:
            row["useful_M1_fraction_of_monitored_abs_flux"] = abs(float(row.get("flux_M1_output_annulus_z", 0.0))) / max(monitored_abs_loss, 1e-300)
        except Exception:
            row["useful_M1_fraction_of_monitored_abs_flux"] = float("nan")
        m1_metrics, m1_rows = analyze_z_profile(sim, dft_m1, g, cfg, g.output_min, g.output_max, z_m1, "M1")
        row.update(m1_metrics)
        # v41.2: diagnostic subwindows inside the same physical V39 output annulus.
        # This does NOT change the geometry or monitor; it tests whether the
        # narrow poloidal core occupies only a radial layer inside 0.90..1.10.
        try:
            subwins = parse_m1_subwindows(cfg.m1_subwindows, g)
            for sw_name, sw_r0, sw_r1 in subwins:
                sw_metrics = analyze_profile_rows_subwindow(m1_rows, sw_r0, sw_r1, f"M1W_{sw_name}")
                row.update(sw_metrics)
            if subwins:
                write_subwindow_summary(cdir / "M1_subwindows_summary.csv", row, subwins)
                plot_subwindow_bars(cdir / "M1_subwindows_summary.csv", cdir / "M1_subwindows_summary.png", f"M1 subwindow diagnostics: {tag}")
        except Exception as e:
            row["M1_subwindows_error"] = repr(e)
        if cfg.save_profiles:
            write_csv(cdir / "M1_output_profile.csv", m1_rows)
            plot_profile(cdir / "M1_output_profile.csv", cdir / "M1_output_profile.png", f"M1 output profile: {tag}")
        if dft_m2 is not None and z_m2 is not None:
            mr0, mr1 = m2_interval(g, cfg, rmax)
            m2_metrics, m2_rows = analyze_z_profile(sim, dft_m2, g, cfg, mr0, mr1, z_m2, "M2")
            row.update(m2_metrics)
            if cfg.save_profiles:
                write_csv(cdir / "M2_downstream_profile.csv", m2_rows)
                plot_profile(cdir / "M2_downstream_profile.csv", cdir / "M2_downstream_profile.png", f"M2 downstream profile: {tag}")
        if cfg.enable_internal_diagnostics and dft_internal is not None:
            try:
                int_metrics, int_rows = analyze_internal_ring_accumulator(sim, dft_internal, g, cfg, "INTERNAL")
                row.update(int_metrics)
                write_csv(cdir / "internal_zone_summary.csv", [int_metrics])
                if cfg.save_profiles:
                    write_csv(cdir / "internal_energy_map.csv", int_rows)
            except Exception as e:
                row["INTERNAL_error"] = repr(e)
        row["candidate_status"] = classify(row, cfg)
    except Exception as e:
        row["candidate_status"] = "ERROR"
        row["error"] = repr(e)
        row["traceback"] = traceback.format_exc()
        if cfg.stop_on_error:
            raise
    finally:
        try:
            if sim is not None:
                sim.reset_meep()
        except Exception:
            pass

    row["elapsed_seconds"] = time.time() - started
    write_json(cdir / "metrics.json", row)
    print("SUMMARY", row.get("candidate_status"), "theta95_M1=", row.get("M1_theta95_poloidal_deg"), "flux_M1=", row.get("flux_M1_output_annulus_z"), flush=True)
    return row


def collect_metrics(outroot: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for p in sorted(outroot.glob("*/metrics.json")):
        try:
            rows.append(json.loads(p.read_text(encoding="utf-8")))
        except Exception:
            pass
    return rows


def write_summary(outroot: Path, rows: List[Dict[str, object]]) -> None:
    write_csv(outroot / "summary.csv", rows)
    write_json(outroot / "summary.json", rows)
    if not rows or plt is None:
        return
    # Simple status plot: M1 theta95 vs model/m.
    plot_rows = [r for r in rows if r.get("M1_theta95_poloidal_deg") not in (None, "")]
    if not plot_rows:
        return
    labels = [f"{r.get('model')}-m{r.get('m')}" for r in plot_rows]
    vals = [float(r.get("M1_theta95_poloidal_deg", float("nan"))) for r in plot_rows]
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(labels)), 4.8), dpi=150)
    ax.bar(np.arange(len(vals)), vals)
    ax.axhline(20, ls=":", lw=1)
    ax.set_ylabel("M1 theta95 poloidal, deg")
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_title("V41 open V39 modal verifier: M1 angular tail")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outroot / "summary_theta95_M1.png")
    plt.close(fig)


def write_methodology(outroot: Path, args, cfg: RunConfig) -> None:
    text = f"""# PHB v41.5 open PHB: modal search + phase-front verification

This folder was produced by `PHB_CO2_v41_4_MODAL_PHASE_FRONT_MEEP.py`.

## Why v41 is new

Earlier v47-style diagnostics were mainly for a strict aperture/halfring or
aperture-side near-field question.  This v41 script transfers the corrected V39
ray-J geometry to an open full-vector Maxwell/FDTD setting.

## V39 geometry used here

- Open PHB: no cylinder closure and no halfring closure.
- Left horn: reflective hyperbolic/linear control horn from z=-L to z=-a.
- Empty central distance: z in [-a,+a].
- Right horn: explicitly present reduced telescope-like reflecting horn with edge radius R_right=R-R2.
  It follows the V39 law rho_right(x)=rho_base(x)-R2 from x=+a to x=x_right_end.
- Primary M1 annular interval: by default r in [R_right, R+R3] at z≈+a. It is diagnostic only; the right half-space remains open for larger r.
- R2 is allowed to be exactly zero: the right funnel then remains full-radius.
- R3 is treated as a continuous left mouth-side radial scale starting at the funnel mouth r=R,z=-a; it is not a detached focal-plane PEC screen.
- Historical diagnostic interval example: a=b=R=1, R2=0.10, R_right=0.90, output_max=1.10, interval width=0.20.

## Parameters of this run

```text
a={args.a}
b={args.b}
R={args.R}
R2={args.R2}
R3={args.R3}, R4={args.R4} (v43.1: R3 fixes the connected left wall and primary annular measurement interval; R4 is reserved and must be 0)
outer_radius={args.outer_radius}
output_r_min={args.output_r_min}
output_r_max={args.output_r_max}
a_over_lambda={cfg.a_over_lambda}
resolution={cfg.resolution}
models={cfg.models}
m_list={cfg.m_list}
source_mode={cfg.source_mode}
source_components={cfg.source_components}
nsrc={cfg.nsrc}
after_sources={cfg.after_sources}
M2 enabled={cfg.enable_m2}
```

## Modal and output metrics

- `flux_M1_output_annulus_z`: MEEP flux through the V39 output annulus.
- `M1_theta50/80/90/95_poloidal_deg`: local meridional Poynting angle percentiles.
- `M1_theta50/80/90/95_3d_deg`: full 3D Poynting angle percentiles including S_phi.
- `M1_L5/L10/L15...`: fraction of forward flux inside selected cones.
- `M1_S_phi_abs_fraction`: diagnostic of azimuthal/vortex-like flow.
- `M1_phase_Ez_rms_raw_rad`, `M1_phase_Ez_rms_after_linear_fit_rad`, `M1_phase_Ez_rms_after_quadratic_fit_rad`: v41.5 phase-front quality diagnostics.  A value above about 1 rad is a strong warning for coherent narrow-beam formation.
- `modal_cluster_summary.csv`: v41.5 Harminv frequency clusters across probe points; use this file to choose frequencies for the subsequent `mode_verify` runs.
- `M1W_*`: subwindow metrics inside the same diagnostic M1 annulus; these do not change geometry and are used to locate the radial core.
- `M1_tube_merit_L15_3d_forward_lowSwirl`: screening metric for a forward, 3D-narrow, low-swirl tubular energy layer.
- `M2_*`: downstream propagation check, when enabled.
- `INTERNAL_eta_eq`, `INTERNAL_eta_axis`, `INTERNAL_energy_peak_r/z`: v41.5 internal equatorial ring-accumulator diagnostics, when `--enable-internal-diagnostics` is used.

## Interpretation rule

A bright field map is not enough.  A useful wave-verification candidate must show
both non-negligible output flux and low angular tails at M1, and preferably must
remain ordered at M2.  PHB must also be compared with the `linear` negative
control under the same geometry and source protocol.
"""
    (outroot / "methodology.md").write_text(text, encoding="utf-8")


# -----------------------------
# CLI
# -----------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--stage", choices=["plan", "geometry", "modal", "mode_verify", "fast", "confirm", "full", "summary"], default="plan")
    ap.add_argument("--outroot", default="/work/PHB_V40_OPEN_TELESCOPE_FULLWAVE")
    ap.add_argument("--models", default="phb,linear", help="phb, linear, or phb,linear")
    ap.add_argument("--m-list", default="0,1,2")
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--b", type=float, default=1.0)
    ap.add_argument("--R", type=float, default=1.0)
    ap.add_argument("--R2", type=float, default=0.10, help="V39 right-funnel radius reduction: R_right=R-R2. v41.3 allows R2=0.")
    ap.add_argument("--R3", type=float, default=0.0, help="left-funnel continuation scale; output_r_max is computed as R+R3")
    ap.add_argument("--R4", type=float, default=0.0, help="reserved compatibility parameter; must be 0 in V43.1; creates no geometry")
    ap.add_argument("--outer-radius", type=float, default=None, help="legacy optional computational radial extent hint only; it does NOT define the output window in V43.1")
    ap.add_argument("--output-r-min", type=float, default=None, help="optional physical output aperture inner radius; default is R_right=R-R2")
    ap.add_argument("--output-r-max", type=float, default=None, help="optional physical output aperture outer radius; default is --outer-radius")
    ap.add_argument("--a-over-lambda", type=float, default=3.0, help="dimensionless scale a/lambda")
    ap.add_argument("--min-a-over-lambda", type=float, default=1.0)
    ap.add_argument("--resolution", type=int, default=80, help="grid cells per unit length")
    ap.add_argument("--dpml", type=float, default=None)
    ap.add_argument("--dpml-over-lambda", type=float, default=1.2)
    ap.add_argument("--wall-thickness", type=float, default=0.035)
    ap.add_argument("--aperture-stop-thickness", type=float, default=None)
    ap.add_argument("--source-mode", choices=["random", "coherent", "single"], default="random")
    ap.add_argument("--source-components", default="Ez,Er,Ep")
    ap.add_argument("--nsrc", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--fwidth-frac", type=float, default=0.18)
    ap.add_argument("--after-sources", type=float, default=160.0)
    ap.add_argument("--aperture-offset-cells", type=float, default=3.0, help="M1 monitor offset after z=+a in grid cells")
    ap.add_argument("--disable-m2", action="store_true", help="disable downstream M2 diagnostic monitor")
    ap.add_argument("--m2-distance", type=float, default=8.0)
    ap.add_argument("--m2-capture-angle-deg", type=float, default=25.0)
    ap.add_argument("--far-zone-safety", type=float, default=1.0)
    ap.add_argument("--allow-near-m2", action="store_true")
    ap.add_argument("--no-profiles", action="store_true")
    ap.add_argument("--no-geometry-png", action="store_true")
    ap.add_argument("--no-archive", action="store_true")
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--stop-on-error", action="store_true")
    ap.add_argument("--narrow-theta95-deg", type=float, default=25.0)
    ap.add_argument("--useful-flux-min", type=float, default=1e-12)
    ap.add_argument("--min-output-cells", type=float, default=8.0)
    ap.add_argument("--allow-underresolved", action="store_true")
    ap.add_argument("--min-free-gb", type=float, default=1.0)
    ap.add_argument("--harminv-component", default="Ez", help="field component for ring-down modal extraction")
    ap.add_argument("--harminv-points", default="", help="modal probe points as r,z;r,z;...; default uses three points inside the central gap")
    ap.add_argument("--modal-after-sources", type=float, default=350.0, help="ring-down time after broad pulse for --stage modal")
    ap.add_argument("--modal-fwidth-frac", type=float, default=0.70, help="broad pulse bandwidth fraction for modal search")
    ap.add_argument("--verify-frequency", type=float, default=0.0, help="frequency to verify in --stage mode_verify; if 0, uses 1/lambda from a_over_lambda")
    ap.add_argument("--verify-fwidth-frac", type=float, default=0.025, help="narrow bandwidth fraction for --stage mode_verify")
    ap.add_argument("--m1-subwindows", default="core_095_105:0.95:1.05,core_097_103:0.97:1.03,upper_100_104:1.00:1.04,lower_090_100:0.90:1.00,upper_100_110:1.00:1.10", help="diagnostic radial subwindows inside M1 output annulus; syntax name:r0:r1,name:r0:r1. Does not change geometry.")
    ap.add_argument("--enable-internal-diagnostics", action="store_true", help="v41.5: add a 2D DFT map in the internal PHB region and compute axis/equatorial/outer energy fractions")
    ap.add_argument("--save-internal-map", action="store_true", help="save full internal_energy_map.csv; otherwise a decimated map is saved when profiles are enabled")
    ap.add_argument("--internal-r-min", type=float, default=None)
    ap.add_argument("--internal-r-max", type=float, default=None)
    ap.add_argument("--internal-z-min", type=float, default=None)
    ap.add_argument("--internal-z-max", type=float, default=None)
    ap.add_argument("--eq-r-min", type=float, default=None, help="equatorial ring-zone inner radius for eta_eq; default is output_r_min")
    ap.add_argument("--eq-r-max", type=float, default=None, help="equatorial ring-zone outer radius for eta_eq; default is output_r_max")
    ap.add_argument("--eq-z-min", type=float, default=None, help="equatorial ring-zone lower z; default is internal_z_min")
    ap.add_argument("--eq-z-max", type=float, default=None, help="equatorial ring-zone upper z; default is internal_z_max")
    ap.add_argument("--axis-r-max", type=float, default=0.25, help="axis-zone radius for eta_axis")
    ap.add_argument("--outer-r-min", type=float, default=None, help="outer-zone inner radius for eta_outer; default is eq_r_max")
    ap.add_argument("--save-density-maps", action="store_true", help="save PNG density maps from the full internal DFT map: raw window, inside PHB, full 2D signed section, and exit/outside region")
    ap.add_argument("--density-map-full-section", action="store_true", help="force full signed-radius 2D section density PNGs; implied by --save-density-maps")
    ap.add_argument("--different-model-seeds", action="store_true", help="use different random source seeds for phb and linear; default v41.3 is matched seeds for fair control")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    models = parse_models_v42(args.models)
    ms = parse_int_list(args.m_list)
    cfg = RunConfig(
        stage=args.stage,
        outroot=args.outroot,
        models=models,
        m_list=ms,
        a_over_lambda=args.a_over_lambda,
        min_a_over_lambda=args.min_a_over_lambda,
        resolution=args.resolution,
        dpml_over_lambda=args.dpml_over_lambda,
        dpml=args.dpml,
        wall_thickness=args.wall_thickness,
        aperture_stop_thickness=args.aperture_stop_thickness,
        source_mode=args.source_mode,
        source_components=args.source_components,
        nsrc=args.nsrc,
        seed=args.seed,
        fwidth_frac=args.fwidth_frac,
        after_sources=args.after_sources,
        aperture_offset_cells=args.aperture_offset_cells,
        enable_m2=not args.disable_m2,
        m2_distance=args.m2_distance,
        m2_capture_angle_deg=args.m2_capture_angle_deg,
        far_zone_safety=args.far_zone_safety,
        allow_near_m2=args.allow_near_m2,
        save_profiles=not args.no_profiles,
        save_geometry_png=not args.no_geometry_png,
        archive=not args.no_archive,
        skip_existing=args.skip_existing,
        stop_on_error=args.stop_on_error,
        narrow_theta95_deg=args.narrow_theta95_deg,
        useful_flux_min=args.useful_flux_min,
        min_output_cells=args.min_output_cells,
        allow_underresolved=args.allow_underresolved,
        min_free_gb=args.min_free_gb,
        harminv_component=args.harminv_component,
        harminv_points=args.harminv_points,
        modal_after_sources=args.modal_after_sources,
        modal_fwidth_frac=args.modal_fwidth_frac,
        verify_frequency=args.verify_frequency,
        verify_fwidth_frac=args.verify_fwidth_frac,
        m1_subwindows=args.m1_subwindows,
        enable_internal_diagnostics=args.enable_internal_diagnostics,
        save_internal_map=args.save_internal_map,
        internal_r_min=args.internal_r_min,
        internal_r_max=args.internal_r_max,
        internal_z_min=args.internal_z_min,
        internal_z_max=args.internal_z_max,
        eq_r_min=args.eq_r_min,
        eq_r_max=args.eq_r_max,
        eq_z_min=args.eq_z_min,
        eq_z_max=args.eq_z_max,
        axis_r_max=args.axis_r_max,
        outer_r_min=args.outer_r_min,
        save_density_maps=args.save_density_maps,
        density_map_full_section=args.density_map_full_section or args.save_density_maps,
        match_model_seeds=not args.different_model_seeds,
    )
    outroot = ensure_dir(Path(cfg.outroot))
    write_methodology(outroot, args, cfg)

    print("============================================================", flush=True)
    print("PHB v41.5 open PHB ring-accumulator diagnostics MEEP verifier", flush=True)
    print(f"outroot={outroot}", flush=True)
    print(f"models={models}; m={ms}; stage={cfg.stage}", flush=True)
    gprint = V39OpenGeometry(model="phb", a=args.a, b=args.b, R=args.R, R2=args.R2, R3=args.R3, R4=args.R4, outer_radius=args.outer_radius, output_r_min=args.output_r_min, output_r_max=args.output_r_max)
    print(f"V43.1 geometry: R2={args.R2}; R3={args.R3}, R4={args.R4} (v43.1: R3 fixes the connected left wall and primary annular measurement interval; R4 is reserved and must be 0); R_right=R-R2={gprint.R_right}; output=[{gprint.output_min},{gprint.output_max}], width={gprint.output_width}", flush=True)
    print("NO halfring. NO cylinder. Central gap is open air. RIGHT V39 reduced reflecting horn is included.", flush=True)
    print("============================================================", flush=True)

    if cfg.stage == "summary":
        rows = collect_metrics(outroot)
        write_summary(outroot, rows)
        if cfg.archive:
            z = archive_folder(outroot)
            if z:
                print("Archive written:", z, flush=True)
        return 0

    rows: List[Dict[str, object]] = []
    for model in models:
        for m in ms:
            print(f"\n=== CASE model={model}, m={m}, stage={cfg.stage} ===", flush=True)
            row = run_one_case(args, cfg, model, m)
            rows.append(row)
            # live summary after each case
            all_rows = collect_metrics(outroot)
            write_summary(outroot, all_rows)
            if row.get("candidate_status") == "ERROR" and cfg.stop_on_error:
                break

    all_rows = collect_metrics(outroot)
    write_summary(outroot, all_rows)
    if cfg.archive:
        z = archive_folder(outroot)
        if z:
            print("Archive written:", z, flush=True)
    bad = [r for r in all_rows if str(r.get("candidate_status", "")).startswith("ERROR")]
    return 1 if bad else 0


# Legacy v41 main retained above for reference; V42 main is defined below.


# ============================================================================
# V42.0 UNIVERSAL FULL-WAVE OVERRIDE LAYER
# ============================================================================
# The original v41.6 implementation is intentionally retained above so that
# earlier data products and formulas remain auditable.  All execution from this
# point uses the V42 functions below.

from dataclasses import dataclass as _v42_dataclass
from typing import Any as _V42Any


@_v42_dataclass
class UniversalRunConfig:
    stage: str
    outroot: str
    models: List[str]
    m_list: List[int]
    a_over_lambda: float
    min_a_over_lambda: float
    resolution: int
    dpml_over_lambda: float
    dpml: Optional[float]
    wall_thickness: float
    aperture_stop_thickness: Optional[float]

    # Excitation
    source_mode: str
    source_waveform: str
    source_duration: float
    source_components: str
    nsrc: int
    seed: int
    source_csv: str
    source_r: Optional[float]
    source_z: float
    source_pair_z: float
    source_r_min: Optional[float]
    source_r_max: Optional[float]
    source_z_min: Optional[float]
    source_z_max: Optional[float]
    source_size_r: float
    source_size_z: float
    input_source_offset_cells: float
    fwidth_frac: float
    after_sources: float

    # Monitor geometry.  The old names are kept as attributes because several
    # legacy plotting helpers use them.
    aperture_offset_cells: float
    enable_m2: bool
    m2_distance: float
    m2_capture_angle_deg: float
    far_zone_safety: float
    allow_near_m2: bool
    enable_left_near: bool
    enable_left_far: bool
    enable_right_near: bool
    enable_right_far: bool
    enable_top: bool
    left_near_z: Optional[float]
    left_far_z: Optional[float]
    right_near_z: Optional[float]
    right_far_z: Optional[float]
    far_distance: float
    axial_monitor_r_min: float
    axial_monitor_r_max: Optional[float]
    top_r: Optional[float]
    top_offset: float
    top_z_min: Optional[float]
    top_z_max: Optional[float]
    top_z_margin: float

    # Output and storage
    save_profiles: bool
    save_geometry_png: bool
    archive: bool
    skip_existing: bool
    stop_on_error: bool
    min_free_gb: float
    save_complex_fields: bool
    save_npz: bool

    # Resolution warnings.  None of these values stop a run in V42.
    narrow_theta95_deg: float
    useful_flux_min: float
    min_output_cells: float
    min_wall_cells: float
    min_extension_cells: float
    min_pml_cells: float
    allow_underresolved: bool

    # Modal diagnostics
    harminv_component: str
    harminv_points: str
    modal_after_sources: float
    modal_fwidth_frac: float
    verify_frequency: float
    verify_fwidth_frac: float
    m1_subwindows: str

    # Full map diagnostics
    enable_internal_diagnostics: bool
    save_internal_map: bool
    internal_r_min: Optional[float]
    internal_r_max: Optional[float]
    internal_z_min: Optional[float]
    internal_z_max: Optional[float]
    map_domain: str
    map_r_margin: float
    map_z_margin: float
    map_stride: int
    exclude_metal_buffer_cells: float
    eq_r_min: Optional[float]
    eq_r_max: Optional[float]
    eq_z_min: Optional[float]
    eq_z_max: Optional[float]
    axis_r_max: float
    outer_r_min: Optional[float]
    focal_band_halfwidth: float
    save_density_maps: bool
    density_map_full_section: bool
    match_model_seeds: bool

    # Live terminal output
    progress_interval: float
    meep_verbosity: int


V42_VERSION = 'v43.1-corrected-output-annulus-and-masks'



def parse_models_v42(s: str) -> List[str]:
    text = str(s).strip().lower()
    if text in ('both', 'all', 'phb+linear'):
        return ['phb', 'linear']
    return parse_models(s)

def _v42_safe_float(x: _V42Any, default: float = float('nan')) -> float:
    try:
        return float(x)
    except Exception:
        return float(default)


def _v42_resolution_grade(cells: float) -> str:
    if not math.isfinite(cells):
        return 'unknown'
    if cells < 4:
        return 'CRITICAL_EXPLORATORY_ONLY'
    if cells < 8:
        return 'UNDERRESOLVED_SCREENING'
    if cells < 16:
        return 'ACCEPTABLE_CONFIRMATION_MINIMUM'
    return 'PUBLICATION_GRADE_TARGET'


def compute_monitor_layout_v42(g: V39OpenGeometry, cfg: UniversalRunConfig) -> Dict[str, float]:
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    dpml = float(cfg.dpml) if cfg.dpml is not None else float(cfg.dpml_over_lambda) * lam
    offset = float(cfg.aperture_offset_cells) / max(int(cfg.resolution), 1)

    # Universal output planes are anchored to the external-focal planes by default.
    # In R1 mode, LEFT_MOUTH is the physical R1 input-aperture plane, not z=-a.
    left_mouth = (g.left_input_z - offset) if g.has_left_input_window else (-g.a - offset)
    right_mouth = +g.a + offset
    left_near = float(cfg.left_near_z) if cfg.left_near_z is not None else -g.c_focus - offset
    right_near = float(cfg.right_near_z) if cfg.right_near_z is not None else +g.c_focus + offset
    left_far = float(cfg.left_far_z) if cfg.left_far_z is not None else left_near - float(cfg.far_distance)
    right_far = float(cfg.right_far_z) if cfg.right_far_z is not None else right_near + float(cfg.far_distance)

    if cfg.stage == 'modal':
        # Modal extraction does not need long propagation planes.
        effective_left_far = left_near
        effective_right_far = right_near
    else:
        effective_left_far = left_far if cfg.enable_left_far else left_near
        effective_right_far = right_far if cfg.enable_right_far else right_near

    monitor_r_min = max(0.0, float(cfg.axial_monitor_r_min))
    if cfg.axial_monitor_r_max is not None:
        monitor_r_max = float(cfg.axial_monitor_r_max)
    else:
        expansion = 0.0
        if cfg.stage != 'modal' and (cfg.enable_left_far or cfg.enable_right_far):
            expansion = float(cfg.far_distance) * math.tan(math.radians(float(cfg.m2_capture_angle_deg)))
        monitor_r_max = max(
            g.radial_extent_radius + max(0.30 * g.R, 0.5 * lam),
            g.output_max + expansion,
        )

    top_r = float(cfg.top_r) if cfg.top_r is not None else (
        g.radial_extent_radius + max(float(cfg.top_offset), 0.20 * g.R, 0.5 * lam)
    )
    top_z_min = float(cfg.top_z_min) if cfg.top_z_min is not None else g.left_z_start - float(cfg.top_z_margin)
    top_z_max = float(cfg.top_z_max) if cfg.top_z_max is not None else max(g.L, g.right_z_end) + float(cfg.top_z_margin)

    radial_needed = max(monitor_r_max, top_r, g.radial_extent_radius)
    radial_margin = max(0.25 * g.R, 0.6 * lam, 4.0 / max(cfg.resolution, 1))
    rmax = radial_needed + dpml + radial_margin

    z_needed_min = min(g.left_z_start, effective_left_far, top_z_min, g.left_extension_z if g.R3 > 0 else g.left_z_start)
    z_needed_max = max(g.right_z_end, effective_right_far, top_z_max)
    axial_margin = max(0.30 * g.R, 0.6 * lam, 4.0 / max(cfg.resolution, 1))
    z_half = max(abs(z_needed_min), abs(z_needed_max)) + dpml + axial_margin
    zspan = 2.0 * z_half

    boundary_margin = max(0.10 * g.R, 0.25 * lam, 2.0 / max(cfg.resolution, 1))
    boundary_r = rmax - dpml - boundary_margin
    boundary_left_z = -z_half + dpml + boundary_margin
    boundary_right_z = +z_half - dpml - boundary_margin

    return {
        'lambda': lam,
        'dpml': dpml,
        'rmax': rmax,
        'zspan': zspan,
        'left_mouth_z': left_mouth,
        'right_mouth_z': right_mouth,
        'left_near_z': left_near,
        'left_far_z': left_far,
        'right_near_z': right_near,
        'right_far_z': right_far,
        'monitor_r_min': monitor_r_min,
        'monitor_r_max': min(monitor_r_max, boundary_r),
        'top_r': min(top_r, boundary_r),
        'top_z_min': max(top_z_min, boundary_left_z),
        'top_z_max': min(top_z_max, boundary_right_z),
        'boundary_r': boundary_r,
        'boundary_left_z': boundary_left_z,
        'boundary_right_z': boundary_right_z,
        'estimated_nr_cells': int(math.ceil(rmax * cfg.resolution)),
        'estimated_nz_cells': int(math.ceil(zspan * cfg.resolution)),
        'estimated_total_cells_2d': int(math.ceil(rmax * cfg.resolution) * math.ceil(zspan * cfg.resolution)),
    }


def validate_geometry_v42(g: V39OpenGeometry, cfg: UniversalRunConfig, layout: Dict[str, float]) -> Dict[str, object]:
    problems: List[str] = []
    warnings: List[str] = []
    checks: List[Dict[str, object]] = []
    lam = layout['lambda']

    if g.a <= 0 or g.b <= 0 or g.R <= 0:
        problems.append('a, b and R must be positive.')
    if g.R1 < 0:
        problems.append('R1 must be >= 0. R1 is the absolute input-aperture radius in the left funnel.')
    if g.R1 > g.R:
        warnings.append(f'R1={g.R1:g} is larger than R={g.R:g}; the effective left input window is clipped to R.')
    if g.R2 < 0:
        problems.append('R2 must be >= 0.')
    if g.R3 < 0:
        problems.append('R3 must be >= 0.')
    if abs(g.R4) > 1e-15:
        problems.append('R4 is reserved in V43.1 and must be 0. The right half-space above the right funnel is fully open air.')
    if g.R_right <= 0:
        problems.append(f'R_right=R-R2={g.R_right:g} must be positive.')
    if g.output_min < 0:
        problems.append('output_r_min must be non-negative.')
    if g.output_max <= g.output_min:
        problems.append(
            'The computed output radial interval is empty. In V43.1 the default output is '
            '[R-R2 if R2>0 else R, R+R3]. Increase R2 or R3, or use an explicit output override only for a special control run.'
        )
    if g.output_min < g.R_right - 1e-12:
        problems.append(
            f'output_r_min={g.output_min:g} cuts into the right reflecting horn edge R_right={g.R_right:g}.'
        )
    if layout['monitor_r_max'] <= layout['monitor_r_min']:
        problems.append('The axial monitor radial interval is empty.')
    if layout['top_z_max'] <= layout['top_z_min']:
        problems.append('The top radial monitor z interval is empty.')

    def check(name: str, physical_size: float, recommended: float, note: str) -> None:
        cells = float(physical_size) * float(cfg.resolution)
        grade = _v42_resolution_grade(cells)
        checks.append({
            'feature': name,
            'physical_size': physical_size,
            'cells': cells,
            'recommended_min_cells': recommended,
            'grade': grade,
            'note': note,
        })
        if cells < recommended:
            warnings.append(
                f'UNDERRESOLVED WARNING: {name} has {cells:.3g} cells; recommended minimum is {recommended:g}. '
                f'The calculation WILL CONTINUE, but the result is {grade}.'
            )

    if g.has_left_input_window:
        check('R1_left_input_window_radius', g.input_window_radius, cfg.min_output_cells,
              'Controls representation of the physical input aperture where the left funnel wall radius equals R1.')
    if g.output_width > 0:
        check('output_annulus_width', g.output_width, cfg.min_output_cells,
              'Controls angular and phase-profile reliability on the output annulus. V43.1 computes this as [R-R2 if R2>0 else R, R+R3] unless explicitly overridden.')
    else:
        checks.append({
            'feature': 'output_annulus_width',
            'physical_size': g.output_width,
            'cells': g.output_width * cfg.resolution,
            'recommended_min_cells': cfg.min_output_cells,
            'grade': 'INVALID_ZERO_WIDTH',
            'note': 'V43.1 computed output width is zero: output=[R-R2 if R2>0 else R, R+R3]. Increase R2 or R3.'
        })
    check('right_horn_axial_length', max(g.right_z_end - g.a, 0.0), 4.0,
          'Controls representation of the complete right reflecting funnel.')
    check('PEC_wall_thickness', effective_wall_t(cfg, lam), cfg.min_wall_cells,
          'Controls representation of the conducting wall.')
    check('PML_thickness', layout['dpml'], cfg.min_pml_cells,
          'Controls suppression of reflections from the computational boundary.')
    if g.R3 > 0:
        check('R3_left_extension', g.R3, cfg.min_extension_cells,
              'Controls the left mouth-side continuous radial extension/window scale R->R+R3.')

    if cfg.a_over_lambda < cfg.min_a_over_lambda:
        warnings.append(
            f'SCALE WARNING: a/lambda={cfg.a_over_lambda:g} is below the requested reference '
            f'{cfg.min_a_over_lambda:g}. The run continues as an exploratory calculation.'
        )
    if cfg.allow_underresolved:
        warnings.append(
            '--allow-underresolved is retained only for command compatibility. In V42, low-cell checks never stop a run.'
        )
    if int(layout.get('estimated_total_cells_2d', 0)) > 2_000_000:
        warnings.append(
            f'MEMORY/TIME WARNING: estimated cylindrical grid has about '
            f'{int(layout["estimated_total_cells_2d"]):,} 2D cells before DFT storage. '
            'Far monitors and full complex maps can substantially increase RAM and runtime.'
        )

    # In mode_verify the requested resonant frequency and a/lambda must describe
    # the same dimensionless Maxwell problem: f = (a/lambda)/a.
    if cfg.stage == 'mode_verify':
        if cfg.verify_frequency <= 0:
            problems.append(
                '--stage mode_verify requires a positive --verify-frequency selected from the modal results.'
            )
        else:
            implied_a_over_lambda = g.a * float(cfg.verify_frequency)
            rel = abs(implied_a_over_lambda - cfg.a_over_lambda) / max(abs(cfg.a_over_lambda), 1e-300)
            if rel > 1e-4:
                problems.append(
                    f'INCONSISTENT FREQUENCY SCALE: --verify-frequency={cfg.verify_frequency:g} with a={g.a:g} '
                    f'implies a/lambda={implied_a_over_lambda:g}, but --a-over-lambda={cfg.a_over_lambda:g}. '
                    'Set --a-over-lambda equal to a*verify_frequency.'
                )
            elif rel > 1e-8:
                warnings.append(
                    f'ROUNDING WARNING: verify-frequency implies a/lambda={implied_a_over_lambda:.12g}, '
                    f'while the command uses {cfg.a_over_lambda:.12g}. The relative difference {rel:.3g} is accepted.'
                )

    if g.model == 'linear':
        warnings.append(
            'The linear geometry is a straight-generator control with the same mouth/end coordinates and aperture protocol; '
            'its enclosed volume is not forced to equal the PHB volume.'
        )

    return {
        'ok': not problems,
        'problems': problems,
        'warnings': warnings,
        'resolution_checks': checks,
        'lambda': lam,
        'a_over_lambda': cfg.a_over_lambda,
        'R_over_lambda': g.R / max(lam, 1e-300),
        'R1_input_window_radius': g.input_window_radius,
        'R1_input_z': g.left_input_z,
        'R1_left_z_start': g.left_z_start,
        'R1_input_window_cells': g.input_window_radius * cfg.resolution,
        'R1_internal_sources_forced_off': bool(g.has_left_input_window),
        'output_width_cells': g.output_width * cfg.resolution,
        'output_computation_rule': 'output_min=R-R2 if R2>0 else R; output_max=R+R3 unless manually overridden',
        'computed_output_outer_from_R3': g.output_max,
        'right_open_half_space_no_PEC_stop': True,
        'R3_cells': g.R3 * cfg.resolution,
        'R4_reserved_value': g.R4,
        'monitor_layout': layout,
        'underresolved_policy': 'WARN_AND_CONTINUE',
    }


def _v42_temporal_source(cfg: UniversalRunConfig, fcen: float, fwidth: float):
    require_meep()
    if cfg.source_waveform == 'continuous':
        duration = max(float(cfg.source_duration), 1.0 / max(fcen, 1e-12))
        return mp.ContinuousSource(frequency=fcen, width=max(fwidth, 1e-12), end_time=duration)
    return mp.GaussianSource(fcen, fwidth=max(fwidth, 1e-12))


def _v42_source_bounds(g: V39OpenGeometry, cfg: UniversalRunConfig, m: int) -> Tuple[float, float, float, float]:
    r0 = float(cfg.source_r_min) if cfg.source_r_min is not None else (0.05 * g.R if abs(m) > 0 else 0.0)
    r1 = float(cfg.source_r_max) if cfg.source_r_max is not None else max(r0 + 1e-6, 0.82 * g.R_right)
    z0 = float(cfg.source_z_min) if cfg.source_z_min is not None else -0.85 * g.a
    z1 = float(cfg.source_z_max) if cfg.source_z_max is not None else +0.85 * g.a
    return r0, r1, z0, z1


def make_sources_v42(g: V39OpenGeometry, cfg: UniversalRunConfig, m: int):
    require_meep()
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    fcen = float(cfg.verify_frequency) if (cfg.stage == 'mode_verify' and cfg.verify_frequency > 0) else 1.0 / lam
    if cfg.stage == 'modal':
        fwidth = cfg.modal_fwidth_frac * fcen
    elif cfg.stage == 'mode_verify':
        fwidth = cfg.verify_fwidth_frac * fcen
    else:
        fwidth = cfg.fwidth_frac * fcen

    comp_names = [x.strip() for x in cfg.source_components.split(',') if x.strip()] or ['Ez']
    comps = [component_from_name(x) for x in comp_names]
    model_offset = 0 if cfg.match_model_seeds else (0 if g.model == 'phb' else 1000)
    rng = np.random.default_rng(int(cfg.seed + 10000 * m + model_offset))
    r0, r1, z0, z1 = _v42_source_bounds(g, cfg, m)
    default_r = float(cfg.source_r) if cfg.source_r is not None else 0.5 * (r0 + r1)
    records: List[Dict[str, object]] = []
    sources: List[object] = []

    aliases = {
        'coherent': 'coherent_volume',
        'single': 'single_ring',
        'ring': 'single_ring',
        'symmetric': 'symmetric_even',
        'antisymmetric': 'symmetric_odd',
        'mirrored_random': 'mirrored_random_even',
    }
    mode = aliases.get(cfg.source_mode, cfg.source_mode)

    def add_source(index: int, kind: str, comp_index: int, r: float, z: float,
                   amp: complex = 1.0 + 0j, size_r: float = 0.0, size_z: float = 0.0,
                   pair_id: str = '') -> None:
        comp_index = int(comp_index) % len(comps)
        comp = comps[comp_index]
        comp_name = comp_names[comp_index]
        center = mp.Vector3(float(r), 0, float(z))
        size = mp.Vector3(max(float(size_r), 0.0), 0, max(float(size_z), 0.0))
        sources.append(mp.Source(_v42_temporal_source(cfg, fcen, fwidth), component=comp, center=center, size=size, amplitude=complex(amp)))
        records.append({
            'source_index': index,
            'source_kind': kind,
            'pair_id': pair_id,
            'component': comp_name,
            'center_r': float(r),
            'center_z': float(z),
            'size_r': float(size_r),
            'size_z': float(size_z),
            'amplitude_real': float(np.real(amp)),
            'amplitude_imag': float(np.imag(amp)),
            'amplitude_abs': float(abs(amp)),
            'phase_rad': float(np.angle(amp)),
            'frequency': fcen,
            'fwidth': fwidth,
            'waveform': cfg.source_waveform,
            'm': m,
            'model': g.model,
        })

    internal_modes = {
        'random','coherent_volume','single_ring','symmetric_even','symmetric_odd',
        'mirrored_random_even','mirrored_random_odd','radial_sheet','axial_sheet',
        'left_sheet','right_sheet'
    }
    if g.has_left_input_window and mode in internal_modes:
        # Hard safety rule requested by the user: if R1>0 the old internal
        # active-gap sources are disabled.  A generic/default source command is
        # automatically reinterpreted as injection through the physical R1
        # input window at the left mouth.
        mode = 'left_input_sheet'

    if mode in ('left_input_sheet','left_input_ring','left_input_center','left_input_edge','left_input_random'):
        if not g.has_left_input_window:
            raise ValueError('left-input source modes require --R1 > 0')
        input_r0 = max(0.0, float(cfg.source_r_min) if cfg.source_r_min is not None else 0.0)
        input_r1 = float(cfg.source_r_max) if cfg.source_r_max is not None else g.input_window_radius
        input_r1 = min(max(input_r1, input_r0 + 1e-9), g.input_window_radius)
        z_in = g.left_input_z - float(cfg.input_source_offset_cells) / max(int(cfg.resolution), 1)
        sheet_size_r = cfg.source_size_r if cfg.source_size_r > 0 else (input_r1 - input_r0)
        sheet_size_z = cfg.source_size_z
        if mode == 'left_input_sheet':
            add_source(0, 'R1_left_input_sheet_forced_internal_off', 0,
                       0.5 * (input_r0 + input_r1), z_in, 1.0 + 0j,
                       sheet_size_r, sheet_size_z, 'R1')
        elif mode == 'left_input_center':
            rr = float(cfg.source_r) if cfg.source_r is not None else (0.0 if abs(m) == 0 else max(0.05 * input_r1, 1.5 / max(int(cfg.resolution), 1)))
            rr = min(max(rr, input_r0), input_r1)
            add_source(0, 'R1_left_input_center_forced_internal_off', 0, rr, z_in, 1.0 + 0j,
                       cfg.source_size_r, cfg.source_size_z, 'R1')
        elif mode == 'left_input_ring':
            rr = float(cfg.source_r) if cfg.source_r is not None else 0.5 * (input_r0 + input_r1)
            rr = min(max(rr, input_r0), input_r1)
            add_source(0, 'R1_left_input_ring_forced_internal_off', 0, rr, z_in, 1.0 + 0j,
                       cfg.source_size_r, cfg.source_size_z, 'R1')
        elif mode == 'left_input_edge':
            rr = float(cfg.source_r) if cfg.source_r is not None else (input_r0 + 0.85 * (input_r1 - input_r0))
            rr = min(max(rr, input_r0), input_r1)
            add_source(0, 'R1_left_input_edge_forced_internal_off', 0, rr, z_in, 1.0 + 0j,
                       cfg.source_size_r, cfg.source_size_z, 'R1')
        elif mode == 'left_input_random':
            nsrc = max(1, int(cfg.nsrc))
            norm = 1.0 / math.sqrt(nsrc)
            for i in range(nsrc):
                rr = float(rng.uniform(input_r0, input_r1))
                # keep the source on/just outside the aperture plane; optional
                # source_size_z gives small user-controlled thickness through
                # the input opening.
                zz = z_in
                phase = float(rng.uniform(0, 2 * math.pi))
                amp = norm * complex(math.cos(phase), math.sin(phase))
                add_source(i, 'R1_left_input_random_forced_internal_off', i % len(comps), rr, zz, amp,
                           0.0, cfg.source_size_z, 'R1')

    elif mode == 'custom_csv':
        if not cfg.source_csv:
            raise ValueError('--source-mode custom_csv requires --source-csv PATH')
        with open(cfg.source_csv, newline='', encoding='utf-8-sig') as f:
            rows = list(csv.DictReader(f))
        if not rows:
            raise ValueError('The custom source CSV is empty.')
        for i, row in enumerate(rows):
            cname = str(row.get('component', comp_names[i % len(comp_names)])).strip()
            if cname not in comp_names:
                comp_names.append(cname)
                comps.append(component_from_name(cname))
            ci = comp_names.index(cname)
            rr = float(row.get('r', row.get('center_r', default_r)))
            zz = float(row.get('z', row.get('center_z', cfg.source_z)))
            if row.get('amplitude_real') not in (None, '') or row.get('amplitude_imag') not in (None, ''):
                amp = complex(float(row.get('amplitude_real', 0.0)), float(row.get('amplitude_imag', 0.0)))
            else:
                mag = float(row.get('amplitude', 1.0))
                phase = float(row.get('phase_rad', 0.0))
                amp = mag * complex(math.cos(phase), math.sin(phase))
            add_source(i, 'custom_csv', ci, rr, zz, amp,
                       float(row.get('size_r', 0.0)), float(row.get('size_z', 0.0)),
                       str(row.get('pair_id', '')))

    elif mode == 'coherent_volume':
        add_source(0, mode, 0, 0.5 * (r0 + r1), 0.5 * (z0 + z1), 1.0 + 0j,
                   cfg.source_size_r if cfg.source_size_r > 0 else (r1 - r0),
                   cfg.source_size_z if cfg.source_size_z > 0 else (z1 - z0))

    elif mode == 'single_ring':
        add_source(0, mode, 0, default_r, cfg.source_z, 1.0 + 0j,
                   cfg.source_size_r, cfg.source_size_z)

    elif mode in ('symmetric_even', 'symmetric_odd'):
        pair_z = min(abs(float(cfg.source_pair_z)), 0.95 * g.a)
        sign = 1.0 if mode == 'symmetric_even' else -1.0
        norm = 1.0 / math.sqrt(2.0)
        add_source(0, mode, 0, default_r, -pair_z, norm + 0j, cfg.source_size_r, cfg.source_size_z, 'pair0')
        add_source(1, mode, 0, default_r, +pair_z, sign * norm + 0j, cfg.source_size_r, cfg.source_size_z, 'pair0')

    elif mode in ('mirrored_random_even', 'mirrored_random_odd'):
        pair_count = max(1, int(math.ceil(cfg.nsrc / 2)))
        sign = 1.0 if mode.endswith('even') else -1.0
        norm = 1.0 / math.sqrt(2.0 * pair_count)
        for p in range(pair_count):
            rr = float(rng.uniform(r0, r1))
            zz = float(rng.uniform(max(0.02 * g.a, 0.0), max(abs(z0), abs(z1))))
            phase = float(rng.uniform(0, 2 * math.pi))
            amp = norm * complex(math.cos(phase), math.sin(phase))
            ci = p % len(comps)
            add_source(2*p, mode, ci, rr, -zz, amp, 0.0, 0.0, f'pair{p}')
            add_source(2*p+1, mode, ci, rr, +zz, sign * amp, 0.0, 0.0, f'pair{p}')

    elif mode in ('radial_sheet', 'left_sheet', 'right_sheet'):
        zz = cfg.source_z
        if mode == 'left_sheet':
            zz = -abs(cfg.source_pair_z)
        elif mode == 'right_sheet':
            zz = +abs(cfg.source_pair_z)
        add_source(0, mode, 0, 0.5 * (r0 + r1), zz, 1.0 + 0j,
                   cfg.source_size_r if cfg.source_size_r > 0 else (r1 - r0),
                   cfg.source_size_z)

    elif mode == 'axial_sheet':
        add_source(0, mode, 0, default_r, 0.5 * (z0 + z1), 1.0 + 0j,
                   cfg.source_size_r,
                   cfg.source_size_z if cfg.source_size_z > 0 else (z1 - z0))

    elif mode == 'random':
        nsrc = max(1, int(cfg.nsrc))
        norm = 1.0 / math.sqrt(nsrc)
        for i in range(nsrc):
            rr = float(rng.uniform(r0, r1))
            zz = float(rng.uniform(z0, z1))
            phase = float(rng.uniform(0, 2 * math.pi))
            amp = norm * complex(math.cos(phase), math.sin(phase))
            add_source(i, 'random_phase_ring', i % len(comps), rr, zz, amp)
    else:
        raise ValueError(
            'Unknown --source-mode. Use random, coherent_volume, single_ring, symmetric_even, '
            'symmetric_odd, mirrored_random_even, mirrored_random_odd, radial_sheet, axial_sheet, '
            'left_sheet, right_sheet, left_input_sheet, left_input_ring, left_input_center, left_input_edge, left_input_random, or custom_csv.'
        )

    return sources, fcen, fwidth, records


def make_simulation_v42(g: V39OpenGeometry, cfg: UniversalRunConfig, m: int):
    require_meep()
    layout = compute_monitor_layout_v42(g, cfg)
    cfg.dpml = layout['dpml']
    sources, fcen, fwidth, source_records = make_sources_v42(g, cfg, m)
    matfun = make_material_function(g, cfg, layout['rmax'])
    courant = min(0.5, 1.0 / (abs(m) + 0.8)) if abs(m) > 0 else 0.5
    try:
        mp.verbosity(int(cfg.meep_verbosity))
    except Exception:
        pass
    sim = mp.Simulation(
        cell_size=mp.Vector3(layout['rmax'], 0, layout['zspan']),
        boundary_layers=[mp.PML(layout['dpml'], direction=mp.R), mp.PML(layout['dpml'], direction=mp.Z)],
        resolution=cfg.resolution,
        dimensions=mp.CYLINDRICAL,
        m=int(m),
        sources=sources,
        material_function=matfun,
        force_complex_fields=True,
        accurate_fields_near_cylorigin=True,
        Courant=courant,
    )
    return sim, layout, fcen, fwidth, source_records


def _v42_map_bounds(g: V39OpenGeometry, cfg: UniversalRunConfig, layout: Dict[str, float]) -> Dict[str, float]:
    br = layout['boundary_r']
    bz0 = layout['boundary_left_z']
    bz1 = layout['boundary_right_z']
    if cfg.map_domain == 'full_cell':
        r0, r1, z0, z1 = 0.0, br, bz0, bz1
    elif cfg.map_domain == 'custom':
        if None in (cfg.internal_r_min, cfg.internal_r_max, cfg.internal_z_min, cfg.internal_z_max):
            raise ValueError('--map-domain custom requires all four --internal-*-min/max bounds.')
        r0, r1 = float(cfg.internal_r_min), float(cfg.internal_r_max)
        z0, z1 = float(cfg.internal_z_min), float(cfg.internal_z_max)
    else:  # resonator
        r0 = 0.0
        r1 = min(br, g.radial_extent_radius + max(cfg.map_r_margin, 0.25 * g.R))
        z0 = max(bz0, -g.L - max(cfg.map_z_margin, 0.25 * g.R))
        z1 = min(bz1, max(g.L, g.right_z_end) + max(cfg.map_z_margin, 0.25 * g.R))
        if cfg.internal_r_min is not None:
            r0 = float(cfg.internal_r_min)
        if cfg.internal_r_max is not None:
            r1 = min(br, float(cfg.internal_r_max))
        if cfg.internal_z_min is not None:
            z0 = max(bz0, float(cfg.internal_z_min))
        if cfg.internal_z_max is not None:
            z1 = min(bz1, float(cfg.internal_z_max))
    if r1 <= r0 or z1 <= z0:
        raise ValueError(f'Invalid map bounds r=[{r0},{r1}], z=[{z0},{z1}]')
    return {'r_min': r0, 'r_max': r1, 'z_min': z0, 'z_max': z1}


def add_full_diagnostics_v42(sim, g: V39OpenGeometry, cfg: UniversalRunConfig,
                             layout: Dict[str, float], fcen: float):
    """Create physically distinct output and leakage diagnostics.

    The narrow annular output M1 is never conflated with the full right-mouth
    plane.  This separation is essential for the V43.1 scientific protocol:

      * M1_OUTPUT_ANNULUS: only output_min <= r <= output_max;
      * RIGHT_MOUTH_FULL: the wider right-mouth plane used to quantify leakage;
      * RIGHT_NEAR/RIGHT_FAR: downstream propagation planes;
      * LEFT_* and TOP: other escape channels;
      * BOUNDARY_*: closed flux budget immediately inside PML.
    """
    require_meep()
    comps = [mp.Er, mp.Ep, mp.Ez, mp.Hr, mp.Hp, mp.Hz]
    dfts: Dict[str, object] = {}
    fluxes: Dict[str, object] = {}
    meta: Dict[str, Dict[str, object]] = {}
    full_r0, full_r1 = layout['monitor_r_min'], layout['monitor_r_max']

    def axial_interval(label: str, z: float, r0: float, r1: float,
                       outward_sign: int, enabled: bool, role: str = '') -> None:
        if not enabled:
            return
        r0 = max(0.0, float(r0))
        r1 = float(r1)
        if r1 <= r0:
            raise ValueError(f'Empty axial monitor {label}: r=[{r0},{r1}]')
        rc, rw = 0.5 * (r0 + r1), r1 - r0
        dfts[label] = sim.add_dft_fields(
            comps, fcen, 0, 1,
            center=mp.Vector3(rc, 0, z), size=mp.Vector3(rw, 0, 0)
        )
        fluxes[label] = sim.add_flux(
            fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(rc, 0, z), size=mp.Vector3(rw, 0, 0), direction=mp.Z)
        )
        meta[label] = {
            'kind': 'axial', 'role': role, 'z': float(z),
            'r_min': r0, 'r_max': r1,
            'outward_sign': int(outward_sign),
            'outward_axis': '-z' if outward_sign < 0 else '+z',
        }

    # Left and full-plane diagnostics.
    axial_interval('LEFT_MOUTH_FULL', layout['left_mouth_z'], full_r0, full_r1,
                   -1, cfg.enable_left_near, 'left_mouth_full_plane')
    axial_interval('RIGHT_MOUTH_FULL', layout['right_mouth_z'], full_r0, full_r1,
                   +1, cfg.enable_right_near, 'right_mouth_full_plane_leakage_control')

    # The principal M1 scientific monitor: a diagnostic annular interval in the fully open right half-space.
    axial_interval('M1_OUTPUT_ANNULUS', layout['right_mouth_z'], g.output_min, g.output_max,
                   +1, cfg.enable_right_near, 'diagnostic_annular_channel_primary_M1_open_right_half_space')

    axial_interval('LEFT_NEAR', layout['left_near_z'], full_r0, full_r1,
                   -1, cfg.enable_left_near, 'left_near_external_focal_reference_plane')
    axial_interval('LEFT_FAR', layout['left_far_z'], full_r0, full_r1,
                   -1, cfg.enable_left_far and cfg.stage != 'modal', 'left_far_plane')
    axial_interval('RIGHT_NEAR', layout['right_near_z'], full_r0, full_r1,
                   +1, cfg.enable_right_near, 'right_near_external_focal_reference_plane')
    axial_interval('RIGHT_FAR', layout['right_far_z'], full_r0, full_r1,
                   +1, cfg.enable_right_far and cfg.stage != 'modal', 'right_far_plane_M2')

    if g.has_left_input_window and cfg.stage != 'modal':
        # Physical R1 aperture monitor at the actual truncated-funnel plane.
        # The source is placed slightly outside this plane; the monitor is
        # shifted slightly inside so it measures power entering the resonator.
        iz = g.left_input_z + float(cfg.aperture_offset_cells) / max(int(cfg.resolution), 1)
        axial_interval('INPUT_R1_LEFT_WINDOW', iz, 0.0, g.input_window_radius,
                       +1, True, 'physical_R1_input_window_into_resonator')

    if cfg.enable_top and cfg.stage != 'modal':
        z0, z1 = layout['top_z_min'], layout['top_z_max']
        zc, zw = 0.5 * (z0 + z1), z1 - z0
        dfts['TOP'] = sim.add_dft_fields(
            comps, fcen, 0, 1,
            center=mp.Vector3(layout['top_r'], 0, zc), size=mp.Vector3(0, 0, zw)
        )
        fluxes['TOP'] = sim.add_flux(
            fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(layout['top_r'], 0, zc), size=mp.Vector3(0, 0, zw), direction=mp.R)
        )
        meta['TOP'] = {
            'kind': 'radial', 'role': 'outer_radial_leakage',
            'r': layout['top_r'], 'z_min': z0, 'z_max': z1,
            'outward_sign': +1, 'outward_axis': '+r',
        }

    # Closed-boundary diagnostics immediately inside PML. Together these three
    # surfaces form a cylindrical closed budget (the r=0 axis has zero area).
    br = layout['boundary_r']
    bl, brr = layout['boundary_left_z'], layout['boundary_right_z']
    fluxes['BOUNDARY_LEFT'] = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(0.5 * br, 0, bl), size=mp.Vector3(br, 0, 0), direction=mp.Z)
    )
    fluxes['BOUNDARY_RIGHT'] = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(0.5 * br, 0, brr), size=mp.Vector3(br, 0, 0), direction=mp.Z)
    )
    fluxes['BOUNDARY_TOP'] = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(br, 0, 0.5 * (bl + brr)), size=mp.Vector3(0, 0, brr - bl), direction=mp.R)
    )
    meta['BOUNDARY_LEFT'] = {'kind': 'boundary_flux', 'role': 'closed_budget', 'outward_sign': -1}
    meta['BOUNDARY_RIGHT'] = {'kind': 'boundary_flux', 'role': 'closed_budget', 'outward_sign': +1}
    meta['BOUNDARY_TOP'] = {'kind': 'boundary_flux', 'role': 'closed_budget', 'outward_sign': +1}

    map_bounds = None
    if cfg.enable_internal_diagnostics and cfg.stage != 'modal':
        map_bounds = _v42_map_bounds(g, cfg, layout)
        dfts['FULL_MAP'] = sim.add_dft_fields(
            comps, fcen, 0, 1,
            center=mp.Vector3(
                0.5 * (map_bounds['r_min'] + map_bounds['r_max']), 0,
                0.5 * (map_bounds['z_min'] + map_bounds['z_max'])
            ),
            size=mp.Vector3(
                map_bounds['r_max'] - map_bounds['r_min'], 0,
                map_bounds['z_max'] - map_bounds['z_min']
            ),
        )
        meta['FULL_MAP'] = {'kind': 'map', 'role': 'complex_full_field_map', **map_bounds}

    return dfts, fluxes, meta, map_bounds

def _v42_orient_array_with_info(arr: np.ndarray, nr_expected: int, nz_expected: int,
                                component_name: str = "") -> Tuple[np.ndarray, Dict[str, object]]:
    """Orient a raw 2D MEEP DFT array to canonical (r,z) indexing.

    MEEP may return a 2D field array whose raw storage order is (z,r) for a
    cylindrical r-z DFT volume.  The scientific code must not assume raw axis
    order.  This function compares the raw shape with expected grid counts
    from the requested physical map bounds:

        expected canonical shape = (nr_expected, nz_expected).

    It returns the oriented array and an audit dictionary that is written to
    full_map_axis_orientation_audit.csv and metrics.json.
    """
    a = np.asarray(arr).squeeze().astype(complex)
    if a.ndim != 2:
        raise ValueError(f"Expected a 2D DFT array, got shape {a.shape}")

    raw0, raw1 = int(a.shape[0]), int(a.shape[1])
    direct_score = abs(raw0 - int(nr_expected)) + abs(raw1 - int(nz_expected))
    transpose_score = abs(raw1 - int(nr_expected)) + abs(raw0 - int(nz_expected))
    use_transpose = transpose_score < direct_score

    # If scores are equal, keep raw order but mark ambiguity.  This should only
    # be acceptable for nearly square maps; final publication runs should use a
    # deliberately non-square diagnostic window or inspect the audit file.
    ambiguous = (transpose_score == direct_score)
    oriented = a.T if use_transpose else a
    orientation = "raw_axes_(z,r)_transposed_to_(r,z)" if use_transpose else "raw_axes_already_(r,z)"
    if ambiguous:
        orientation += "_AMBIGUOUS_EQUAL_SCORE"

    info: Dict[str, object] = {
        "component": component_name,
        "raw_shape": f"{raw0}x{raw1}",
        "raw_shape_0": raw0,
        "raw_shape_1": raw1,
        "expected_canonical_shape": f"{int(nr_expected)}x{int(nz_expected)}",
        "expected_nr": int(nr_expected),
        "expected_nz": int(nz_expected),
        "direct_score_raw_as_rz": int(direct_score),
        "transpose_score_raw_as_zr": int(transpose_score),
        "transposed": bool(use_transpose),
        "ambiguous": bool(ambiguous),
        "orientation": orientation,
        "oriented_shape": f"{int(oriented.shape[0])}x{int(oriented.shape[1])}",
    }
    return oriented, info


def _v42_orient_array(arr: np.ndarray, nr_expected: int, nz_expected: int) -> np.ndarray:
    """Backward-compatible wrapper: return only oriented (r,z) array."""
    return _v42_orient_array_with_info(arr, nr_expected, nz_expected, "")[0]

def analyze_axial_profile_v42(sim, dft, r_min: float, r_max: float, z_plane: float,
                              outward_sign: int, label: str) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    require_meep()
    arrays = {name: dft_array(sim, dft, comp) for name, comp in (
        ('Er', mp.Er), ('Ep', mp.Ep), ('Ez', mp.Ez), ('Hr', mp.Hr), ('Hp', mp.Hp), ('Hz', mp.Hz)
    )}
    n = min(len(x) for x in arrays.values())
    if n <= 0:
        return {f'{label}_valid': False}, []
    for k in arrays:
        arrays[k] = arrays[k][:n]
    Er, Ep, Ez = arrays['Er'], arrays['Ep'], arrays['Ez']
    Hr, Hp, Hz = arrays['Hr'], arrays['Hp'], arrays['Hz']
    r = np.linspace(r_min, r_max, n)
    S_r = 0.5 * np.real(Ep * np.conj(Hz) - Ez * np.conj(Hp))
    S_phi = 0.5 * np.real(Ez * np.conj(Hr) - Er * np.conj(Hz))
    S_z = 0.5 * np.real(Er * np.conj(Hp) - Ep * np.conj(Hr))
    S_out = float(outward_sign) * S_z
    E2 = np.abs(Er)**2 + np.abs(Ep)**2 + np.abs(Ez)**2
    H2 = np.abs(Hr)**2 + np.abs(Hp)**2 + np.abs(Hz)**2
    u = 0.25 * (E2 + H2)
    # Cylindrical meridional area element is proportional to r*dr.  Earlier
    # versions omitted dr, which was harmless for percentiles on one monitor
    # but made flux proxies from monitors of different radial widths
    # quantitatively incomparable.  Trapezoidal endpoint weights fix this.
    dr = (r_max - r_min) / max(n - 1, 1)
    quadrature = np.full(n, dr, dtype=float)
    if n > 1:
        quadrature[0] *= 0.5
        quadrature[-1] *= 0.5
    area = np.maximum(r, 0.0) * quadrature
    forward_w = np.maximum(S_out, 0.0) * area
    abs_s = np.sqrt(S_r*S_r + S_phi*S_phi + S_z*S_z)
    abs_w = abs_s * area
    theta_pol = np.degrees(np.arctan2(np.abs(S_r), np.maximum(S_out, 1e-300)))
    theta_3d = np.degrees(np.arctan2(np.sqrt(S_r*S_r + S_phi*S_phi), np.maximum(S_out, 1e-300)))
    theta_pol = np.where(S_out > 0, theta_pol, 180.0)
    theta_3d = np.where(S_out > 0, theta_3d, 180.0)

    comp_power = {
        'Er': float(np.sum(np.abs(Er)**2 * forward_w)),
        'Ep': float(np.sum(np.abs(Ep)**2 * forward_w)),
        'Ez': float(np.sum(np.abs(Ez)**2 * forward_w)),
    }
    dominant = max(comp_power, key=comp_power.get)
    phase_map = {'Er': np.angle(Er), 'Ep': np.angle(Ep), 'Ez': np.angle(Ez)}
    total_forward = float(np.sum(forward_w))
    total_abs = float(np.sum(abs_w))
    signed_out = float(np.sum(S_out * area))
    metrics: Dict[str, object] = {
        f'{label}_valid': True,
        f'{label}_kind': 'axial_output_monitor',
        f'{label}_outward_axis': '-z' if outward_sign < 0 else '+z',
        f'{label}_z': z_plane,
        f'{label}_r_min': r_min,
        f'{label}_r_max': r_max,
        f'{label}_samples': n,
        f'{label}_signed_outward_flux_proxy': signed_out,
        f'{label}_positive_outward_flux_proxy': total_forward,
        f'{label}_outward_fraction_of_abs': total_forward / max(total_abs, 1e-300),
        f'{label}_backward_fraction_of_abs_z': float(np.sum(np.maximum(-S_out, 0.0)*area) / max(np.sum(np.abs(S_out)*area), 1e-300)),
        f'{label}_u_weighted_mean': float(np.sum(u*area) / max(np.sum(area), 1e-300)),
        f'{label}_u_peak': float(np.nanmax(u)),
        f'{label}_u_peak_r': float(r[int(np.nanargmax(u))]),
        f'{label}_dominant_E_component': dominant,
        f'{label}_S_r_abs_fraction': float(np.sum(np.abs(S_r)*area) / max(np.sum((np.abs(S_r)+np.abs(S_phi)+np.abs(S_z))*area), 1e-300)),
        f'{label}_S_phi_abs_fraction': float(np.sum(np.abs(S_phi)*area) / max(np.sum((np.abs(S_r)+np.abs(S_phi)+np.abs(S_z))*area), 1e-300)),
        f'{label}_S_z_abs_fraction': float(np.sum(np.abs(S_z)*area) / max(np.sum((np.abs(S_r)+np.abs(S_phi)+np.abs(S_z))*area), 1e-300)),
    }
    add_phase_front_metrics(metrics, label, r, np.angle(Ez), forward_w, 'Ez')
    add_phase_front_metrics(metrics, label, r, phase_map[dominant], forward_w, 'dominantE')
    for p in (50, 80, 90, 95, 99):
        metrics[f'{label}_theta{p}_poloidal_deg'] = weighted_percentile(theta_pol, forward_w, p/100)
        metrics[f'{label}_theta{p}_3d_deg'] = weighted_percentile(theta_3d, forward_w, p/100)
    for th in (1, 5, 10, 15, 20, 25, 30, 45):
        metrics[f'{label}_L{th}_poloidal_fraction'] = float(np.sum(forward_w[theta_pol <= th]) / max(total_forward, 1e-300))
        metrics[f'{label}_L{th}_3d_fraction'] = float(np.sum(forward_w[theta_3d <= th]) / max(total_forward, 1e-300))

    rows: List[Dict[str, object]] = []
    for i in range(n):
        row = {
            'i': i, 'r': float(r[i]), 'z': float(z_plane),
            'E2': float(E2[i]), 'H2': float(H2[i]), 'u_timeavg': float(u[i]),
            'S_r': float(S_r[i]), 'S_phi': float(S_phi[i]), 'S_z': float(S_z[i]),
            'S_outward': float(S_out[i]),
            'theta_poloidal_outward_deg': float(theta_pol[i]),
            'theta_3d_outward_deg': float(theta_3d[i]),
            'outward_weight': float(forward_w[i]),
        }
        for name, arr in arrays.items():
            row[f'{name}_real'] = float(np.real(arr[i]))
            row[f'{name}_imag'] = float(np.imag(arr[i]))
            row[f'{name}_abs'] = float(abs(arr[i]))
            row[f'{name}_phase_rad'] = float(np.angle(arr[i]))
        rows.append(row)
    return metrics, rows


def analyze_radial_profile_v42(sim, dft, r_plane: float, z_min: float, z_max: float,
                               label: str) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    require_meep()
    arrays = {name: dft_array(sim, dft, comp) for name, comp in (
        ('Er', mp.Er), ('Ep', mp.Ep), ('Ez', mp.Ez), ('Hr', mp.Hr), ('Hp', mp.Hp), ('Hz', mp.Hz)
    )}
    n = min(len(x) for x in arrays.values())
    if n <= 0:
        return {f'{label}_valid': False}, []
    for k in arrays:
        arrays[k] = arrays[k][:n]
    Er, Ep, Ez = arrays['Er'], arrays['Ep'], arrays['Ez']
    Hr, Hp, Hz = arrays['Hr'], arrays['Hp'], arrays['Hz']
    z = np.linspace(z_min, z_max, n)
    S_r = 0.5 * np.real(Ep * np.conj(Hz) - Ez * np.conj(Hp))
    S_phi = 0.5 * np.real(Ez * np.conj(Hr) - Er * np.conj(Hz))
    S_z = 0.5 * np.real(Er * np.conj(Hp) - Ep * np.conj(Hr))
    E2 = np.abs(Er)**2 + np.abs(Ep)**2 + np.abs(Ez)**2
    H2 = np.abs(Hr)**2 + np.abs(Hp)**2 + np.abs(Hz)**2
    u = 0.25 * (E2 + H2)
    # Radial cylindrical surface element is proportional to r*dz.
    dz = (z_max - z_min) / max(n - 1, 1)
    quadrature = np.full(n, dz, dtype=float)
    if n > 1:
        quadrature[0] *= 0.5
        quadrature[-1] *= 0.5
    weight_area = max(r_plane, 1e-12) * quadrature
    outward_w = np.maximum(S_r, 0.0) * weight_area
    abs_s = np.sqrt(S_r*S_r + S_phi*S_phi + S_z*S_z)
    abs_w = abs_s * weight_area
    theta_meridional = np.degrees(np.arctan2(np.abs(S_z), np.maximum(S_r, 1e-300)))
    theta_3d = np.degrees(np.arctan2(np.sqrt(S_z*S_z + S_phi*S_phi), np.maximum(S_r, 1e-300)))
    theta_meridional = np.where(S_r > 0, theta_meridional, 180.0)
    theta_3d = np.where(S_r > 0, theta_3d, 180.0)
    total_out = float(np.sum(outward_w))
    metrics: Dict[str, object] = {
        f'{label}_valid': True,
        f'{label}_kind': 'radial_output_monitor',
        f'{label}_outward_axis': '+r',
        f'{label}_r': r_plane,
        f'{label}_z_min': z_min,
        f'{label}_z_max': z_max,
        f'{label}_samples': n,
        f'{label}_signed_outward_flux_proxy': float(np.sum(S_r*weight_area)),
        f'{label}_positive_outward_flux_proxy': total_out,
        f'{label}_outward_fraction_of_abs': total_out / max(float(np.sum(abs_w)), 1e-300),
        f'{label}_u_peak': float(np.nanmax(u)),
        f'{label}_u_peak_z': float(z[int(np.nanargmax(u))]),
        f'{label}_S_r_abs_fraction': float(np.sum(np.abs(S_r)*weight_area) / max(np.sum((np.abs(S_r)+np.abs(S_phi)+np.abs(S_z))*weight_area), 1e-300)),
        f'{label}_S_phi_abs_fraction': float(np.sum(np.abs(S_phi)*weight_area) / max(np.sum((np.abs(S_r)+np.abs(S_phi)+np.abs(S_z))*weight_area), 1e-300)),
        f'{label}_S_z_abs_fraction': float(np.sum(np.abs(S_z)*weight_area) / max(np.sum((np.abs(S_r)+np.abs(S_phi)+np.abs(S_z))*weight_area), 1e-300)),
    }
    for p in (50, 80, 90, 95, 99):
        metrics[f'{label}_theta{p}_meridional_deg'] = weighted_percentile(theta_meridional, outward_w, p/100)
        metrics[f'{label}_theta{p}_3d_deg'] = weighted_percentile(theta_3d, outward_w, p/100)
    rows: List[Dict[str, object]] = []
    for i in range(n):
        row = {
            'i': i, 'r': float(r_plane), 'z': float(z[i]),
            'E2': float(E2[i]), 'H2': float(H2[i]), 'u_timeavg': float(u[i]),
            'S_r': float(S_r[i]), 'S_phi': float(S_phi[i]), 'S_z': float(S_z[i]),
            'S_outward': float(S_r[i]),
            'theta_meridional_outward_deg': float(theta_meridional[i]),
            'theta_3d_outward_deg': float(theta_3d[i]),
            'outward_weight': float(outward_w[i]),
        }
        for name, arr in arrays.items():
            row[f'{name}_real'] = float(np.real(arr[i]))
            row[f'{name}_imag'] = float(np.imag(arr[i]))
            row[f'{name}_abs'] = float(abs(arr[i]))
            row[f'{name}_phase_rad'] = float(np.angle(arr[i]))
        rows.append(row)
    return metrics, rows


def _v42_geometry_masks(g: V39OpenGeometry, cfg: UniversalRunConfig,
                        Rg: np.ndarray, Zg: np.ndarray) -> Dict[str, np.ndarray]:
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    wall_t = effective_wall_t(cfg, lam)
    stop_t = aperture_stop_t(cfg, lam)
    buffer = max(float(cfg.exclude_metal_buffer_cells), 0.0) / max(cfg.resolution, 1)

    def build(expand: float) -> np.ndarray:
        """Analytical metal mask matching make_material_function exactly."""
        st = stop_t + 2.0 * expand
        mask = np.zeros_like(Rg, dtype=bool)

        # Left physical cut edge / mouth continuation.  No metal is ever added
        # at z=-c_focus; the focal coordinate is only a geometric reference.
        if g.has_left_input_window:
            left_outer = max(g.input_window_radius, g.left_mouth_pec_outer_radius)
            mask |= ((np.abs(Zg - g.left_input_z) <= 0.5*st) &
                     (Rg >= g.input_window_radius - expand) &
                     (Rg <= left_outer + wall_t + expand))
        else:
            left_outer = g.left_mouth_pec_outer_radius
            near_left_mouth = np.abs(Zg + g.a) <= 0.5*st
            left_inner = np.full_like(Rg, float(g.R), dtype=float)
            for j in range(Zg.shape[1]):
                zz = float(Zg[0, j])
                if np.any(near_left_mouth[:, j]):
                    left_inner[:, j] = g.left_mouth_bridge_inner_radius(zz)
            mask |= (near_left_mouth &
                     (Rg >= left_inner - expand) &
                     (Rg <= left_outer + wall_t + expand))

        # Right side is fully open above the funnel edge.  The M1 annulus is
        # diagnostic only and does not create metal at output_max or beyond.
        # Keep an inner diaphragm only for an explicit output_min>R_right control.
        near_right_mouth = np.abs(Zg - g.a) <= 0.5*st
        if g.output_min > g.R_right + 0.25*wall_t:
            mask |= near_right_mouth & (Rg >= g.R_right - expand) & (Rg <= g.output_min + expand)
        return mask

    metal = build(0.0)
    metal_buffer = build(buffer)
    # Add curved wall masks row by row in z because wall_rho is scalar.
    for j in range(Zg.shape[1]):
        zz = float(Zg[0, j])
        rw = g.wall_rho(zz)
        if math.isfinite(rw):
            metal[:, j] |= (Rg[:, j] >= max(0.0, rw)) & (Rg[:, j] <= max(0.0, rw) + wall_t)
            metal_buffer[:, j] |= (Rg[:, j] >= max(0.0, rw - buffer)) & (Rg[:, j] <= max(0.0, rw) + wall_t + buffer)

    wall_rho_grid = np.full_like(Rg, np.nan, dtype=float)
    inside_horn = np.zeros_like(Rg, dtype=bool)
    for j in range(Zg.shape[1]):
        zz = float(Zg[0, j])
        rw = g.wall_rho(zz)
        if math.isfinite(rw):
            wall_rho_grid[:, j] = rw
            inside_horn[:, j] = Rg[:, j] < rw
    central_inner = (np.abs(Zg) < g.a) & (Rg < g.R)
    return {
        'metal': metal,
        'metal_or_buffer': metal_buffer,
        'vacuum_safe': ~metal_buffer,
        'wall_rho': wall_rho_grid,
        'inside_horn_envelope': inside_horn,
        'central_inner': central_inner,
    }


def _v42_zone_rows(g: V39OpenGeometry, cfg: UniversalRunConfig, Rg: np.ndarray, Zg: np.ndarray,
                   u: np.ndarray, S_r: np.ndarray, S_phi: np.ndarray, S_z: np.ndarray,
                   masks: Dict[str, np.ndarray]) -> List[Dict[str, object]]:
    safe = masks['vacuum_safe']
    # Axisymmetric volume element is proportional to r*dr*dz.  Trapezoidal
    # quadrature makes integrated energy/volume proxies stable against grid
    # refinement and map-size changes (the global 2*pi factor cancels).
    r_axis = np.asarray(Rg[:, 0], dtype=float)
    z_axis = np.asarray(Zg[0, :], dtype=float)
    dr = (r_axis[-1] - r_axis[0]) / max(len(r_axis) - 1, 1)
    dz = (z_axis[-1] - z_axis[0]) / max(len(z_axis) - 1, 1)
    qr = np.full(len(r_axis), dr, dtype=float)
    qz = np.full(len(z_axis), dz, dtype=float)
    if len(qr) > 1:
        qr[0] *= 0.5; qr[-1] *= 0.5
    if len(qz) > 1:
        qz[0] *= 0.5; qz[-1] *= 0.5
    W = np.maximum(Rg, 0.0) * qr[:, None] * qz[None, :]
    total_u = float(np.sum(u[safe] * W[safe]))
    total_v = float(np.sum(W[safe]))
    ring_outer = max(g.R, g.output_max)
    eq_r0 = g.output_min if cfg.eq_r_min is None else float(cfg.eq_r_min)
    eq_r1 = g.output_max if cfg.eq_r_max is None else float(cfg.eq_r_max)
    eq_z0 = -0.1*g.a if cfg.eq_z_min is None else float(cfg.eq_z_min)
    eq_z1 = +0.1*g.a if cfg.eq_z_max is None else float(cfg.eq_z_max)
    focal_half = max(float(cfg.focal_band_halfwidth), 1.0/max(cfg.resolution,1))

    zones: List[Tuple[str, np.ndarray, str]] = [
        ('all_safe_vacuum', safe, 'All mapped vacuum cells excluding the requested metal buffer.'),
        ('left_half', safe & (Zg < 0), 'All safe vacuum at z<0.'),
        ('right_half', safe & (Zg >= 0), 'All safe vacuum at z>=0.'),
        ('axis_zone', safe & (Rg <= cfg.axis_r_max), 'Near-axis control zone.'),
        ('central_inner_gap', safe & (np.abs(Zg) <= g.a) & (Rg <= g.R), 'Central inner gap below radius R.'),
        ('equatorial_thin_zone', safe & (Rg >= eq_r0) & (Rg <= eq_r1) & (Zg >= eq_z0) & (Zg <= eq_z1), 'User/pre-registered thin equatorial zone.'),
        ('interfocal_annular_corridor', safe & (np.abs(Zg) <= g.c_focus) & (Rg >= g.R) & (Rg <= ring_outer), 'Annular corridor between the two external focal planes.'),
        ('left_focal_band', safe & (np.abs(Zg + g.c_focus) <= focal_half) & (Rg >= g.R) & (Rg <= ring_outer), 'Vacuum band adjacent to the left focal plane.'),
        ('right_focal_band', safe & (np.abs(Zg - g.c_focus) <= focal_half) & (Rg >= g.R) & (Rg <= ring_outer), 'Vacuum band adjacent to the right focal plane.'),
        ('left_ring_channel', safe & (Zg < 0) & (Rg >= g.R) & (Rg <= ring_outer), 'Left annular air channel at and outside radius R.'),
        ('right_ring_channel', safe & (Zg >= 0) & (Rg >= g.R) & (Rg <= ring_outer), 'Right annular air channel at and outside radius R.'),
        ('left_horn_interior', safe & (Zg <= -g.a) & (Zg >= -g.L) & masks['inside_horn_envelope'], 'Vacuum inside the left funnel.'),
        ('right_horn_interior', safe & (Zg >= g.a) & (Zg <= g.right_z_end) & masks['inside_horn_envelope'], 'Vacuum inside the right funnel.'),
        ('outside_nominal_radius', safe & (Rg > ring_outer), 'Mapped external/side-leakage region beyond the outer annular radius.'),
        ('near_metal_buffer', masks['metal_or_buffer'] & (~masks['metal']), 'Vacuum-adjacent numerical buffer around PEC surfaces; excluded from safe-volume integrals.'),
    ]

    out: List[Dict[str, object]] = []
    for name, mask, description in zones:
        n = int(np.sum(mask))
        row: Dict[str, object] = {'zone': name, 'description': description, 'point_count': n}
        if n == 0:
            row.update({'valid': False})
            out.append(row)
            continue
        zone_u = float(np.sum(u[mask] * W[mask]))
        zone_v = float(np.sum(W[mask]))
        flat_u = np.where(mask, u, -np.inf)
        pi, pj = np.unravel_index(int(np.nanargmax(flat_u)), u.shape)
        sr = S_r[mask]; sp = S_phi[mask]; sz = S_z[mask]; ww = W[mask]
        sden = max(float(np.sum((np.abs(sr)+np.abs(sp)+np.abs(sz))*ww)), 1e-300)
        row.update({
            'valid': True,
            'volume_weight_fraction': zone_v / max(total_v, 1e-300),
            'energy_fraction': zone_u / max(total_u, 1e-300),
            'energy_density_enrichment': (zone_u/max(zone_v,1e-300)) / max(total_u/max(total_v,1e-300), 1e-300),
            'energy_weighted_r_mean': float(np.sum(Rg[mask]*u[mask]*ww) / max(zone_u,1e-300)),
            'energy_weighted_z_mean': float(np.sum(Zg[mask]*u[mask]*ww) / max(zone_u,1e-300)),
            'peak_r': float(Rg[pi,pj]), 'peak_z': float(Zg[pi,pj]), 'peak_u': float(u[pi,pj]),
            'S_r_signed_integral_proxy': float(np.sum(sr*ww)),
            'S_phi_signed_integral_proxy': float(np.sum(sp*ww)),
            'S_z_signed_integral_proxy': float(np.sum(sz*ww)),
            'S_r_abs_fraction': float(np.sum(np.abs(sr)*ww)/sden),
            'S_phi_abs_fraction': float(np.sum(np.abs(sp)*ww)/sden),
            'S_z_abs_fraction': float(np.sum(np.abs(sz)*ww)/sden),
            'forward_plus_z_fraction': float(np.sum(np.maximum(sz,0)*ww)/max(np.sum(np.abs(sz)*ww),1e-300)),
            'radially_outward_fraction': float(np.sum(np.maximum(sr,0)*ww)/max(np.sum(np.abs(sr)*ww),1e-300)),
        })
        out.append(row)
    return out


def save_full_map_products_v42(cdir: Path, g: V39OpenGeometry, cfg: UniversalRunConfig,
                               r: np.ndarray, z: np.ndarray, arrays: Dict[str, np.ndarray],
                               masks: Dict[str, np.ndarray]) -> None:
    if plt is None or not cfg.save_density_maps:
        return
    ensure_dir(cdir)
    Rg, Zg = np.meshgrid(r, z, indexing='ij')
    u = arrays['u_timeavg']
    safe = masks['vacuum_safe']

    def draw(data: np.ndarray, name: str, title: str, signed: bool = False, log: bool = False) -> None:
        A = np.array(data, dtype=float)
        if log:
            A = np.log10(np.maximum(A, 0.0) + 1e-30)
        finite = np.isfinite(A)
        if not np.any(finite):
            return
        fig, ax = plt.subplots(figsize=(10, 7), dpi=170)
        if signed:
            vmax = float(np.nanpercentile(np.abs(A[finite]), 99.0))
            vmax = max(vmax, 1e-30)
            im = ax.imshow(A.T, origin='lower', extent=[r.min(), r.max(), z.min(), z.max()],
                           aspect='auto', vmin=-vmax, vmax=vmax, cmap='coolwarm')
        else:
            lo = float(np.nanpercentile(A[finite], 1.0))
            hi = float(np.nanpercentile(A[finite], 99.5))
            im = ax.imshow(A.T, origin='lower', extent=[r.min(), r.max(), z.min(), z.max()],
                           aspect='auto', vmin=lo, vmax=hi)
        ax.axhline(-g.c_focus, lw=0.8, ls='--')
        ax.axhline(+g.c_focus, lw=0.8, ls='--')
        ax.axhline(-g.a, lw=0.8, ls=':')
        ax.axhline(+g.a, lw=0.8, ls=':')
        ax.axvline(g.R, lw=0.8, ls=':')
        ax.set_xlabel('r')
        ax.set_ylabel('z')
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(cdir / name)
        plt.close(fig)

    draw(u, 'MAP_01_energy_density_log.png', 'Time-averaged electromagnetic energy density log10(u)', log=True)
    draw(np.where(safe, u, np.nan), 'MAP_02_energy_density_safe_vacuum_log.png', 'Energy density excluding PEC and buffer cells', log=True)
    draw(arrays['E2'], 'MAP_03_E2_log.png', 'Electric-field intensity |E|²', log=True)
    draw(arrays['H2'], 'MAP_04_H2_log.png', 'Magnetic-field intensity |H|²', log=True)
    draw(arrays['S_r'], 'MAP_05_Sr_signed.png', 'Radial Poynting component S_r', signed=True)
    draw(arrays['S_phi'], 'MAP_06_Sphi_signed.png', 'Azimuthal Poynting component S_phi', signed=True)
    draw(arrays['S_z'], 'MAP_07_Sz_signed.png', 'Axial Poynting component S_z', signed=True)
    draw(arrays['S_abs'], 'MAP_08_Poynting_magnitude_log.png', 'Poynting-vector magnitude |S|', log=True)

    # Meridional flow map: horizontal arrow is S_r, vertical arrow is S_z.
    step_r = max(1, len(r)//45)
    step_z = max(1, len(z)//65)
    rr = Rg[::step_r, ::step_z]
    zz = Zg[::step_r, ::step_z]
    sr = arrays['S_r'][::step_r, ::step_z]
    sz = arrays['S_z'][::step_r, ::step_z]
    mag = np.sqrt(sr*sr + sz*sz)
    threshold = float(np.nanpercentile(mag[np.isfinite(mag)], 40)) if np.any(np.isfinite(mag)) else 0.0
    keep = mag >= threshold
    fig, ax = plt.subplots(figsize=(10, 7), dpi=170)
    base = np.log10(np.maximum(u, 0.0) + 1e-30)
    im = ax.imshow(base.T, origin='lower', extent=[r.min(), r.max(), z.min(), z.max()], aspect='auto')
    ax.quiver(rr[keep], zz[keep], sr[keep], sz[keep], angles='xy', scale_units='xy', scale=None, width=0.002)
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_title('Meridional Poynting flow (arrows: S_r, S_z) over log energy density')
    fig.colorbar(im, ax=ax, label='log10(u)')
    fig.tight_layout()
    fig.savefig(cdir / 'MAP_09_meridional_flow_quiver.png')
    plt.close(fig)


def analyze_full_map_v42(sim, dft, g: V39OpenGeometry, cfg: UniversalRunConfig,
                         bounds: Dict[str, float], cdir: Path) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    """Analyze the full internal/external field map with explicit axis audit.

    v42.1 critical correction
    -------------------------
    All raw 2D DFT component arrays are explicitly oriented to the canonical
    scientific convention

        array.shape == (nr, nz),  first index -> r,  second index -> z.

    The raw MEEP shape, expected shape, transpose decision, and ambiguity flag
    are saved in full_map_axis_orientation_audit.csv and in metrics.json.  This
    prevents the V41.5/V41.6 error where a raw (z,r) array was analyzed as (r,z).
    """
    require_meep()
    raw = {name: dft_array_raw(sim, dft, comp) for name, comp in (
        ("Er", mp.Er), ("Ep", mp.Ep), ("Ez", mp.Ez), ("Hr", mp.Hr), ("Hp", mp.Hp), ("Hz", mp.Hz)
    )}

    nr_expected = max(2, int(round((bounds["r_max"] - bounds["r_min"]) * cfg.resolution)) + 1)
    nz_expected = max(2, int(round((bounds["z_max"] - bounds["z_min"]) * cfg.resolution)) + 1)

    oriented: Dict[str, np.ndarray] = {}
    orientation_rows: List[Dict[str, object]] = []
    for name, arr0 in raw.items():
        arr, info = _v42_orient_array_with_info(arr0, nr_expected, nz_expected, name)
        oriented[name] = arr
        orientation_rows.append(info)
    write_csv(cdir / "full_map_axis_orientation_audit.csv", orientation_rows)

    nr = min(v.shape[0] for v in oriented.values())
    nz = min(v.shape[1] for v in oriented.values())
    for k in oriented:
        oriented[k] = oriented[k][:nr, :nz]

    r = np.linspace(bounds["r_min"], bounds["r_max"], nr)
    z = np.linspace(bounds["z_min"], bounds["z_max"], nz)
    Rg, Zg = np.meshgrid(r, z, indexing="ij")

    Er, Ep, Ez = oriented["Er"], oriented["Ep"], oriented["Ez"]
    Hr, Hp, Hz = oriented["Hr"], oriented["Hp"], oriented["Hz"]
    E2 = np.abs(Er)**2 + np.abs(Ep)**2 + np.abs(Ez)**2
    H2 = np.abs(Hr)**2 + np.abs(Hp)**2 + np.abs(Hz)**2
    u = 0.25*(E2 + H2)

    S_r = 0.5*np.real(Ep*np.conj(Hz) - Ez*np.conj(Hp))
    S_phi = 0.5*np.real(Ez*np.conj(Hr) - Er*np.conj(Hz))
    S_z = 0.5*np.real(Er*np.conj(Hp) - Ep*np.conj(Hr))
    S_abs = np.sqrt(S_r*S_r + S_phi*S_phi + S_z*S_z)

    masks = _v42_geometry_masks(g, cfg, Rg, Zg)
    zone_rows = _v42_zone_rows(g, cfg, Rg, Zg, u, S_r, S_phi, S_z, masks)
    write_csv(cdir / "internal_zone_summary.csv", zone_rows)

    safe = masks["vacuum_safe"]
    if not np.any(safe):
        return {
            "FULL_MAP_valid": False,
            "FULL_MAP_reason": "no_safe_vacuum_points_after_metal_buffer_mask",
            "FULL_MAP_expected_nr_from_r_range": int(nr_expected),
            "FULL_MAP_expected_nz_from_z_range": int(nz_expected),
        }, []

    peak_safe = np.where(safe, u, -np.inf)
    pi, pj = np.unravel_index(int(np.nanargmax(peak_safe)), u.shape)
    W = np.maximum(Rg, 0.0)

    first = orientation_rows[0] if orientation_rows else {}
    any_transposed = any(bool(x.get("transposed", False)) for x in orientation_rows)
    any_ambiguous = any(bool(x.get("ambiguous", False)) for x in orientation_rows)
    all_orientations = "; ".join(f"{x.get('component')}:{x.get('orientation')}" for x in orientation_rows)

    metrics: Dict[str, object] = {
        "FULL_MAP_valid": True,
        "FULL_MAP_axis_convention": "canonical_oriented_array_is_(r,z); first_index_r; second_index_z",
        "FULL_MAP_axis_orientation_audit_file": "full_map_axis_orientation_audit.csv",
        "FULL_MAP_expected_nr_from_r_range": int(nr_expected),
        "FULL_MAP_expected_nz_from_z_range": int(nz_expected),
        "FULL_MAP_raw_shape_Er": str(first.get("raw_shape", "")),
        "FULL_MAP_orientation_Er": str(first.get("orientation", "")),
        "FULL_MAP_transposed_Er": bool(first.get("transposed", False)),
        "FULL_MAP_any_component_transposed": bool(any_transposed),
        "FULL_MAP_any_axis_orientation_ambiguous": bool(any_ambiguous),
        "FULL_MAP_axis_orientation_all_components": all_orientations,
        "FULL_MAP_nr": int(nr), "FULL_MAP_nz": int(nz),
        "FULL_MAP_r_min": bounds["r_min"], "FULL_MAP_r_max": bounds["r_max"],
        "FULL_MAP_z_min": bounds["z_min"], "FULL_MAP_z_max": bounds["z_max"],
        "FULL_MAP_safe_energy_weighted": float(np.sum(u[safe]*W[safe])),
        "FULL_MAP_safe_volume_weighted": float(np.sum(W[safe])),
        "FULL_MAP_safe_peak_r": float(Rg[pi,pj]),
        "FULL_MAP_safe_peak_z": float(Zg[pi,pj]),
        "FULL_MAP_safe_peak_u": float(u[pi,pj]),
        "FULL_MAP_safe_peak_E2": float(E2[pi,pj]),
        "FULL_MAP_safe_peak_H2": float(H2[pi,pj]),
        "FULL_MAP_safe_peak_S_r": float(S_r[pi,pj]),
        "FULL_MAP_safe_peak_S_phi": float(S_phi[pi,pj]),
        "FULL_MAP_safe_peak_S_z": float(S_z[pi,pj]),
        "FULL_MAP_metal_buffer_cells": cfg.exclude_metal_buffer_cells,
        "FULL_MAP_metal_point_count": int(np.sum(masks["metal"])),
        "FULL_MAP_metal_or_buffer_point_count": int(np.sum(masks["metal_or_buffer"])),
        "FULL_MAP_safe_vacuum_point_count": int(np.sum(masks["vacuum_safe"])),
    }
    if any_ambiguous:
        metrics["FULL_MAP_axis_orientation_warning"] = (
            "Axis orientation score is ambiguous for at least one component. "
            "Use a non-square map window or manually inspect full_map_axis_orientation_audit.csv before scientific interpretation."
        )

    for zr in zone_rows:
        if zr.get("valid"):
            prefix = "ZONE_" + str(zr["zone"])
            for k in ("energy_fraction","volume_weight_fraction","energy_density_enrichment","peak_r","peak_z","peak_u",
                      "S_r_signed_integral_proxy","S_phi_signed_integral_proxy","S_z_signed_integral_proxy"):
                metrics[f"{prefix}_{k}"] = zr.get(k)

    arrays = {
        **oriented, "E2": E2, "H2": H2, "u_timeavg": u,
        "S_r": S_r, "S_phi": S_phi, "S_z": S_z, "S_abs": S_abs,
    }
    if cfg.save_npz:
        np.savez_compressed(
            cdir / "full_field_map_arrays.npz", r=r, z=z,
            **arrays,
            metal_mask=masks["metal"], metal_or_buffer_mask=masks["metal_or_buffer"],
            vacuum_safe_mask=masks["vacuum_safe"], wall_rho=masks["wall_rho"],
        )
    save_full_map_products_v42(cdir, g, cfg, r, z, arrays, masks)

    rows: List[Dict[str, object]] = []
    stride = 1 if cfg.save_internal_map else max(1, int(cfg.map_stride))
    ring_outer = max(g.R, g.output_max)
    for i in range(0, nr, stride):
        for j in range(0, nz, stride):
            row: Dict[str, object] = {
                "array_i_r": int(i), "array_j_z": int(j),
                "r": float(r[i]), "z": float(z[j]),
                "E2": float(E2[i,j]), "H2": float(H2[i,j]), "u_timeavg": float(u[i,j]),
                "field_intensity_proxy_E2_plus_H2": float(E2[i,j]+H2[i,j]),
                "S_r": float(S_r[i,j]), "S_phi": float(S_phi[i,j]), "S_z": float(S_z[i,j]), "S_abs": float(S_abs[i,j]),
                "metal_cell_mask": bool(masks["metal"][i,j]),
                "metal_or_buffer_mask": bool(masks["metal_or_buffer"][i,j]),
                "vacuum_safe_mask": bool(masks["vacuum_safe"][i,j]),
                "wall_rho_at_z": float(masks["wall_rho"][i,j]) if math.isfinite(float(masks["wall_rho"][i,j])) else "",
                "inside_horn_envelope": bool(masks["inside_horn_envelope"][i,j]),
                "central_inner": bool(masks["central_inner"][i,j]),
                "interfocal_corridor": bool(abs(z[j]) <= g.c_focus and g.R <= r[i] <= ring_outer),
                "left_focal_band": bool(abs(z[j]+g.c_focus) <= cfg.focal_band_halfwidth and g.R <= r[i] <= ring_outer),
                "right_focal_band": bool(abs(z[j]-g.c_focus) <= cfg.focal_band_halfwidth and g.R <= r[i] <= ring_outer),
            }
            if cfg.save_complex_fields:
                for name, arr in oriented.items():
                    row[f"{name}_real"] = float(np.real(arr[i,j]))
                    row[f"{name}_imag"] = float(np.imag(arr[i,j]))
                    row[f"{name}_abs"] = float(abs(arr[i,j]))
                    row[f"{name}_phase_rad"] = float(np.angle(arr[i,j]))
            rows.append(row)
    write_csv(cdir / "full_field_map.csv", rows)
    return metrics, rows

def make_progress_callback_v42(cdir: Path, case_label: str, interval: float):
    progress_path = cdir / 'runtime_progress.csv'
    if not progress_path.exists():
        write_csv(progress_path, [], fieldnames=['timestamp','case','meep_time','wall_elapsed_seconds'])
    started = time.time()

    def callback(sim):
        try:
            mt = float(sim.meep_time())
        except Exception:
            mt = float('nan')
        elapsed = time.time() - started
        print(
            f'[MEEP PROGRESS] case={case_label}  meep_time={mt:.6g}  wall_elapsed={elapsed:.1f}s',
            flush=True,
        )
        try:
            with progress_path.open('a', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                w.writerow([now_iso(), case_label, mt, elapsed])
        except Exception:
            pass
    return mp.at_every(max(float(interval), 0.1), callback)


def plot_geometry_v42(cdir: Path, g: V39OpenGeometry, cfg: UniversalRunConfig,
                      layout: Dict[str, float]) -> None:
    """Publication-safe full meridional 2D geometry plot.

    Physical PEC surfaces alone are solid. Open-air apertures are never
    drawn as solid segments. Diagnostic planes are thin dashed lines. Geometric focal coordinates
    are markers only and are never drawn as material.
    """
    if plt is None or not cfg.save_geometry_png:
        return

    fig, ax = plt.subplots(figsize=(12.5, 7.2), dpi=170)
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    wall_t = effective_wall_t(cfg, lam)

    # Present reflecting walls only. If R1>0 the left horn is truncated at the
    # plane where rho_wall=R1.
    zl = np.linspace(g.left_z_start, -g.a, 700)
    rl = np.array([g.wall_rho(float(x)) for x in zl])
    zr = np.linspace(g.a, g.right_z_end, 700)
    rr = np.array([g.wall_rho(float(x)) for x in zr])

    left_label = 'left PHB reflecting funnel' if g.model == 'phb' else 'left straight-generator reflecting funnel'
    right_label = 'right PHB reflecting funnel' if g.model == 'phb' else 'right straight-generator reflecting funnel'
    ax.plot(zl, +rl, lw=2.8, label=left_label)
    ax.plot(zl, -rl, lw=2.8)
    ax.plot(zr, +rr, lw=2.8, label=right_label)
    ax.plot(zr, -rr, lw=2.8)

    # Left physical opening/cut edge.
    if g.has_left_input_window:
        zin = g.left_input_z
        ax.plot([zin, zin], [-g.input_window_radius, +g.input_window_radius],
                lw=7.0, alpha=0.35, solid_capstyle='round', label='OPEN R1 input aperture')
        left_outer = max(g.input_window_radius, g.left_mouth_pec_outer_radius)
        if left_outer > g.input_window_radius:
            ax.plot([zin, zin], [g.input_window_radius, left_outer], lw=4.0,
                    label='left R1 cut-edge PEC rim')
            ax.plot([zin, zin], [-left_outer, -g.input_window_radius], lw=4.0)
        zsrc = zin - float(cfg.input_source_offset_cells) / max(int(cfg.resolution), 1)
        ax.plot([zsrc, zsrc], [-g.input_window_radius, +g.input_window_radius],
                ls=':', lw=1.8, label='R1 source plane outside aperture')
    else:
        # The interval -R..R is open.  Do not draw any solid segment there.
        # Only the physical R3 PEC continuation is solid.
        ax.annotate('INTERNAL OPEN MOUTH (air)', xy=(-g.a, 0.0),
                    xytext=(-g.a - 0.45 * g.R, 0.35 * g.R),
                    arrowprops=dict(arrowstyle='->', lw=1.1), fontsize=8.5)
        if g.left_mouth_pec_outer_radius > g.R:
            ax.plot([-g.a, -g.a], [g.R, g.left_mouth_pec_outer_radius], lw=4.0, color='saddlebrown',
                    label='PHYSICAL PEC: R3 connected left-mouth continuation')
            ax.plot([-g.a, -g.a], [-g.left_mouth_pec_outer_radius, -g.R], lw=4.0, color='saddlebrown')

    # Right side is fully open above the right-funnel edge.  Never draw the
    # measurement interval as a solid bar.  It is air, not geometry.
    ax.annotate('OPEN AIR above right funnel; NO PEC screen',
                xy=(g.a + 0.02 * g.R, g.output_max + 0.03 * g.R),
                xytext=(g.a + 0.50 * g.R, g.output_max + 0.42 * g.R),
                arrowprops=dict(arrowstyle='->', lw=1.2), fontsize=9.0, fontweight='bold')
    if g.output_min > g.R_right:
        ax.plot([g.a, g.a], [g.R_right, g.output_min], lw=4.0,
                label='right inner PEC diaphragm')
        ax.plot([g.a, g.a], [-g.output_min, -g.R_right], lw=4.0)

    # Monitors. Full planes are dashed across their full capture interval.
    monitor_specs = [
        ('LEFT_MOUTH_FULL', layout['left_mouth_z'], cfg.enable_left_near),
        ('RIGHT_MOUTH_FULL', layout['right_mouth_z'], cfg.enable_right_near),
        ('LEFT_NEAR', layout['left_near_z'], cfg.enable_left_near),
        ('LEFT_FAR', layout['left_far_z'], cfg.enable_left_far and cfg.stage != 'modal'),
        ('RIGHT_NEAR', layout['right_near_z'], cfg.enable_right_near),
        ('RIGHT_FAR/M2', layout['right_far_z'], cfg.enable_right_far and cfg.stage != 'modal'),
    ]
    for name, zz, enabled in monitor_specs:
        if enabled:
            ax.plot([zz, zz], [-layout['monitor_r_max'], layout['monitor_r_max']],
                    ls='--', lw=1.05, alpha=0.70, label=name)

    # Dedicated M1 monitor: thin cyan dashed line in air, never a solid bar.
    if cfg.enable_right_near:
        zm1 = layout['right_mouth_z']
        ax.plot([zm1, zm1], [g.output_min, g.output_max], ls='--', lw=1.8, color='tab:cyan',
                label='VIRTUAL M1 monitor in air (NO material)')
        ax.plot([zm1, zm1], [-g.output_max, -g.output_min], ls='--', lw=1.8, color='tab:cyan')

    if cfg.enable_top and cfg.stage != 'modal':
        ax.plot([layout['top_z_min'], layout['top_z_max']], [layout['top_r'], layout['top_r']],
                ls='--', lw=1.25, alpha=0.75, label='TOP +rho monitor')
        ax.plot([layout['top_z_min'], layout['top_z_max']], [-layout['top_r'], -layout['top_r']],
                ls='--', lw=1.25, alpha=0.45)

    # Focal-ring coordinates are physical only as geometric reference points.
    if g.model == 'phb':
        focus_label = 'PHB external focal-ring meridional points'
    else:
        focus_label = 'PHB focal-coordinate references (linear control has no hyperbolic foci)'
    ax.scatter([-g.c_focus, g.c_focus], [g.R, g.R], s=42, marker='x', label=focus_label)
    ax.scatter([-g.c_focus, g.c_focus], [-g.R, -g.R], s=42, marker='x')
    ax.axhline(0, lw=0.8, color='0.25')
    ax.axhline(g.R, ls=':', lw=1.0, color='0.45')
    ax.axhline(-g.R, ls=':', lw=1.0, color='0.45')
    ax.axvspan(-g.a, g.a, alpha=0.045, label='open central gap')

    x_min = min(g.left_z_start, layout['left_far_z'] if cfg.enable_left_far else g.left_z_start) - 0.10*g.R
    x_max = max(g.right_z_end, layout['right_far_z'] if cfg.enable_right_far and cfg.stage != 'modal' else g.right_z_end) + 0.10*g.R
    y_max = min(layout['boundary_r'], max(layout['monitor_r_max'], g.radial_extent_radius + 0.25*g.R))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-y_max, y_max)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.22)
    ax.set_xlabel('z: longitudinal axis / PHB focus axis')
    ax.set_ylabel('rho: meridional radius (mirrored +/- for visualization)')
    ax.set_title(
        f'V43.2 corrected full meridional geometry: model={g.model}; '
        f'output={g.output_min:g}...{g.output_max:g}; M1 is virtual; right half-space is open'
    )
    ax.legend(fontsize=7.0, loc='best')
    fig.tight_layout()
    fig.savefig(cdir / 'geometry_full_2D_meridional_corrected.png')
    fig.savefig(cdir / 'geometry.png')
    plt.close(fig)

def write_methodology_v42(outroot: Path, args, cfg: UniversalRunConfig) -> None:
    g = V39OpenGeometry(
        model='phb', a=args.a, b=args.b, R=args.R, R1=args.R1,
        R2=args.R2, R3=args.R3, R4=args.R4,
        outer_radius=args.outer_radius,
        output_r_min=args.output_r_min, output_r_max=args.output_r_max,
    )
    expected_f = cfg.a_over_lambda / max(g.a, 1e-300)
    text = f"""# PHB V43.1 corrected universal full-wave verification

This directory is produced by `PHB_CO2_v43_2_OPEN_RIGHT_CLEAN_PLOT_MEEP.py`.

## Fixed geometry semantics

- `a` is the constant mouth coordinate and half-gap. It is never split into a separate hyperbola parameter and a gap.
- `R2` is the absolute reduction of the complete right-funnel radius: `R_right = R - R2`. `R2=0` keeps a full-radius right funnel.
- `R1` is an absolute physical input-aperture radius inside the left funnel. If `R1>0`, the left funnel is truncated where its wall radius equals `R1`; old internal sources are replaced by left-side aperture injection.
- `R3` is the continuous left-mouth PEC continuation from `R` to `R+R3`. It also defines the outer edge of the primary right-side diagnostic annulus, but creates no right-side material.
- The left PHB mouth uses a grid-connected PEC bridge: on horn-side cells within the mouth slab, metal begins at the local `wall_rho(z)` and continues through `R+R3`. This prevents the singular PHB mouth slope from creating a false annular leak.
- The interval `output_min = R-R2` for `R2>0`, otherwise `R`, through `output_max = R+R3` is the primary annular measurement interval. It is not bounded by an outer PEC screen; air continues for `r>output_max`.
- `R4` is reserved and must be zero. No right-side PEC stop, support or screen is created.
- No PEC material is placed at `z=±c_focus`. External focal rings are geometric reference coordinates only.
- The linear control uses the same mouth coordinates, axial end coordinates, aperture protocol, material, sources and monitors, but replaces the hyperbolic generator by a straight line. Equal enclosed volume is not imposed.

## Corrected diagnostics

- `M1_OUTPUT_ANNULUS` is the primary output monitor and covers only the physical radial interval `{g.output_min:g} <= r <= {g.output_max:g}`.
- `RIGHT_MOUTH_FULL` is a separate full-plane leakage/control monitor and is never aliased to M1.
- `RIGHT_FAR` is the downstream M2 monitor when enabled.
- LEFT_MOUTH_FULL, LEFT_NEAR, LEFT_FAR and TOP quantify other escape channels.
- BOUNDARY_LEFT, BOUNDARY_RIGHT and BOUNDARY_TOP form the signed closed-boundary flux budget immediately inside PML.
- The full-map metal mask is generated from the same connected left-mouth bridge and curved-wall equations as the MEEP material function. It contains no fictitious focal-plane metal and no right outer stop.

## Frequency convention

For a consistent `mode_verify` run, `a_over_lambda = a * verify_frequency`.
For the current scale, `a_over_lambda={cfg.a_over_lambda:g}` corresponds to `frequency={expected_f:g}`.
A mismatch is treated as an invalid physical setup rather than silently accepted.

## Resolution policy

Low cell counts are written as warnings and do not automatically stop exploratory runs. Publication conclusions require a grid-convergence repeat. The classification uses signed forward flux through the dedicated M1 annulus and never converts backward flux into useful output by taking an absolute value.

## Current run

Geometry: a={args.a}, b={args.b}, R={args.R}, R1={args.R1}, R2={args.R2}, R3={args.R3}, R4={args.R4}
Computed output annulus: [{g.output_min}, {g.output_max}], width={g.output_width}
Stage: {cfg.stage}
Models: {','.join(cfg.models)}
Modes m: {','.join(str(x) for x in cfg.m_list)}
Source mode: {cfg.source_mode}; waveform: {cfg.source_waveform}
Resolution: {cfg.resolution}
Map domain: {cfg.map_domain}
"""
    (outroot/'METHODOLOGY_V42.md').write_text(text, encoding='utf-8')

def run_one_case_v42(args, cfg: UniversalRunConfig, model: str, m: int) -> Dict[str, object]:
    g = V39OpenGeometry(
        model=model, a=args.a, b=args.b, R=args.R, R1=args.R1, R2=args.R2, R3=args.R3, R4=args.R4,
        outer_radius=args.outer_radius, output_r_min=args.output_r_min, output_r_max=args.output_r_max
    )
    outroot = Path(cfg.outroot)
    tag = (
        f'{model}_m{m}_R1_{g.R1:.5g}_R2_{g.R2:.5g}_R3_{g.R3:.5g}_R4_{g.R4:.5g}_'
        f'out_{g.output_min:.5g}_{g.output_max:.5g}_src_{cfg.source_mode}_'
        f'aol_{cfg.a_over_lambda:.5g}_res{cfg.resolution}'
    ).replace('.','p')
    cdir = ensure_dir(outroot/tag)
    done = cdir/'metrics.json'
    if cfg.skip_existing and done.exists():
        old = json.loads(done.read_text(encoding='utf-8'))
        old['skipped_existing'] = True
        return old

    started = time.time()
    layout = compute_monitor_layout_v42(g, cfg)
    validation = validate_geometry_v42(g, cfg, layout)
    row: Dict[str, object] = {
        'timestamp': now_iso(), 'script': Path(__file__).name, 'version': V42_VERSION,
        'stage': cfg.stage, 'model': model, 'm': m,
        'a': g.a, 'b': g.b, 'R': g.R, 'R1': g.R1, 'R1_input_window_radius': g.input_window_radius,
        'R1_internal_sources_forced_off': bool(g.has_left_input_window),
        'R2': g.R2, 'R3': g.R3, 'R4': g.R4,
        'R_right': g.R_right, 'output_r_min': g.output_min, 'output_r_max': g.output_max,
        'output_width': g.output_width, 'L': g.L, 'c_focus': g.c_focus,
        'a_over_lambda': cfg.a_over_lambda, 'resolution': cfg.resolution,
        'source_mode': cfg.source_mode, 'source_waveform': cfg.source_waveform,
    }
    row.update({f'layout_{k}': v for k,v in layout.items()})
    row.update({f'validation_{k}': json.dumps(v, ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in validation.items()})
    write_json(cdir/'config.json', {'geometry': asdict(g), 'config': asdict(cfg), 'layout': layout, 'validation': validation})
    write_csv(cdir/'resolution_warnings.csv', validation['resolution_checks'])
    write_geometry_audit(cdir, g, cfg)
    plot_geometry_v42(cdir, g, cfg, layout)

    print(f'\n[V42 CASE] {tag}', flush=True)
    print(f'[GEOMETRY] a={g.a}, b={g.b}, R={g.R}, R1={g.R1}, R2={g.R2}, R3={g.R3}, R4={g.R4}', flush=True)
    if g.has_left_input_window:
        print(f'[R1 INPUT] left input window open radius={g.input_window_radius:g}; internal source modes are forced to left-input injection.', flush=True)
    print(f'[MONITORS] Lnear={layout["left_near_z"]:.6g}, Lfar={layout["left_far_z"]:.6g}, '
          f'Rnear={layout["right_near_z"]:.6g}, Rfar={layout["right_far_z"]:.6g}, top_r={layout["top_r"]:.6g}', flush=True)
    for warning in validation['warnings']:
        print('[WARNING]', warning, flush=True)
    if not validation['ok']:
        row['candidate_status'] = 'INVALID_PHYSICAL_GEOMETRY'
        row['error'] = '; '.join(validation['problems'])
        row['elapsed_seconds'] = time.time()-started
        write_json(done,row)
        return row
    if cfg.stage in ('plan','geometry'):
        row['candidate_status'] = 'GEOMETRY_ONLY'
        row['elapsed_seconds'] = time.time()-started
        write_json(done,row)
        return row

    require_meep()
    if free_gb(outroot) < cfg.min_free_gb:
        raise RuntimeError(f'Not enough free disk space: {free_gb(outroot):.2f} GB < {cfg.min_free_gb} GB')
    sim = None
    try:
        sim, layout, fcen, fwidth, source_records = make_simulation_v42(g,cfg,m)
        row.update({'fcen':fcen,'fwidth':fwidth})
        write_csv(cdir/'sources.csv', source_records)
        print(f'[MEEP START] cell rmax={layout["rmax"]:.6g}, zspan={layout["zspan"]:.6g}, '
              f'fcen={fcen:.9g}, fwidth={fwidth:.9g}', flush=True)
        if cfg.stage == 'modal':
            row.update(run_modal_v42(sim,g,cfg,cdir,fcen,fwidth,m,tag))
            row['candidate_status'] = 'MODAL_RINGDOWN_DONE' if int(row.get('modal_mode_count',0))>0 else 'MODAL_NO_MODES_FOUND'
        else:
            dfts, flux_objs, meta, map_bounds = add_full_diagnostics_v42(sim,g,cfg,layout,fcen)
            progress = make_progress_callback_v42(cdir, tag, cfg.progress_interval)
            sim.run(progress, until_after_sources=float(cfg.after_sources))

            flux_rows: List[Dict[str, object]] = []
            outward_boundary_sum = 0.0
            for name,obj in flux_objs.items():
                raw_flux = float(mp.get_fluxes(obj)[0])
                sign = int(meta.get(name,{}).get('outward_sign',1))
                outward_flux = sign*raw_flux
                row[f'flux_{name}_raw_coordinate_sign'] = raw_flux
                row[f'flux_{name}_outward_signed'] = outward_flux
                flux_rows.append({
                    'monitor': name,
                    'kind': meta.get(name,{}).get('kind',''),
                    'role': meta.get(name,{}).get('role',''),
                    'raw_coordinate_flux': raw_flux,
                    'outward_sign_multiplier': sign,
                    'outward_signed_flux': outward_flux,
                })
                if name == 'M1_OUTPUT_ANNULUS':
                    # Stable historical/publication key.  This is the signed
                    # outward MEEP flux through the actual physical annulus.
                    row['flux_M1_output_annulus_z'] = outward_flux
                    row['flux_M1_output_annulus_z_raw_coordinate_sign'] = raw_flux
                if name in ('BOUNDARY_LEFT','BOUNDARY_RIGHT','BOUNDARY_TOP'):
                    outward_boundary_sum += outward_flux
            row['closed_boundary_outward_flux_sum'] = outward_boundary_sum
            write_csv(cdir/'flux_monitor_summary.csv', flux_rows)

            for label,dft in dfts.items():
                if label == 'FULL_MAP':
                    continue
                mmeta = meta[label]
                if mmeta['kind'] == 'axial':
                    metrics, prows = analyze_axial_profile_v42(
                        sim,dft,float(mmeta['r_min']),float(mmeta['r_max']),float(mmeta['z']),
                        int(mmeta['outward_sign']),label
                    )
                else:
                    metrics, prows = analyze_radial_profile_v42(
                        sim,dft,float(mmeta['r']),float(mmeta['z_min']),float(mmeta['z_max']),label
                    )
                row.update(metrics)
                if cfg.save_profiles:
                    write_csv(cdir/f'{label}_field_profile.csv', prows)
                # Backward-compatible aliases.  Historical/publication M1 now
                # maps ONLY to the dedicated diagnostic M1 annulus, never to
                # the full right-mouth leakage plane.  M2 remains RIGHT_FAR.
                if label == 'M1_OUTPUT_ANNULUS':
                    for k,v in metrics.items():
                        if k.startswith('M1_OUTPUT_ANNULUS_'):
                            row['M1_'+k[len('M1_OUTPUT_ANNULUS_'):]] = v
                    if cfg.save_profiles:
                        write_csv(cdir/'M1_output_profile.csv', prows)
                    try:
                        subwins = parse_m1_subwindows(cfg.m1_subwindows, g)
                        for name,r0,r1 in subwins:
                            row.update(analyze_profile_rows_subwindow(prows,r0,r1,f'M1W_{name}'))
                        write_subwindow_summary(cdir/'M1_subwindows_summary.csv', row, subwins)
                        plot_subwindow_bars(cdir/'M1_subwindows_summary.csv', cdir/'M1_subwindows_summary.png', f'M1 annular subwindows: {tag}')
                    except Exception as sub_exc:
                        row['M1_subwindows_error'] = repr(sub_exc)
                if label == 'RIGHT_FAR':
                    for k,v in metrics.items():
                        if k.startswith('RIGHT_FAR_'):
                            row['M2_'+k[len('RIGHT_FAR_'):]] = v

            # Quantify how much of the full right-mouth outward flux belongs to
            # the physical M1 annulus.  Values >1 can occur only when the rest
            # of the plane contains net backward flow; both signed and positive
            # profile proxies remain available for interpretation.
            m1_pos = float(row.get('M1_positive_outward_flux_proxy', 0.0) or 0.0)
            full_pos = float(row.get('RIGHT_MOUTH_FULL_positive_outward_flux_proxy', 0.0) or 0.0)
            row['M1_fraction_of_right_mouth_positive_outward_flux_proxy'] = m1_pos / max(full_pos, 1e-300)

            if 'FULL_MAP' in dfts and map_bounds is not None:
                map_metrics, _ = analyze_full_map_v42(sim,dfts['FULL_MAP'],g,cfg,map_bounds,cdir)
                row.update(map_metrics)
            row['candidate_status'] = classify_v42(row,cfg)
    except Exception as exc:
        row['candidate_status'] = 'ERROR'
        row['error'] = repr(exc)
        row['traceback'] = traceback.format_exc()
        print('[ERROR]', repr(exc), flush=True)
        print(row['traceback'], flush=True)
        if cfg.stop_on_error:
            raise
    finally:
        try:
            if sim is not None:
                sim.reset_meep()
        except Exception:
            pass
    row['elapsed_seconds'] = time.time()-started
    write_json(done,row)
    print(f'[V42 SUMMARY] status={row.get("candidate_status")} elapsed={row["elapsed_seconds"]:.1f}s '
          f'Lmouth_full={row.get("LEFT_MOUTH_FULL_positive_outward_flux_proxy")} '
          f'Rmouth_full={row.get("RIGHT_MOUTH_FULL_positive_outward_flux_proxy")} '
          f'M1_annulus={row.get("M1_positive_outward_flux_proxy")} '
          f'Lfocal={row.get("LEFT_NEAR_positive_outward_flux_proxy")} '
          f'Rfocal={row.get("RIGHT_NEAR_positive_outward_flux_proxy")} '
          f'TOP={row.get("TOP_positive_outward_flux_proxy")}', flush=True)
    return row


def build_parser_v42() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description='PHB V43.1 corrected cylindrical full-vector FDTD verifier with annulus-only M1, full leakage planes, axis-safe maps and modal search.'
    )
    ap.add_argument('--stage', choices=['plan','geometry','modal','mode_verify','fast','confirm','full','summary'], default='plan')
    ap.add_argument('--outroot', default='/work/PHB_V43_1_OPEN_RIGHT_NO_STOP')
    ap.add_argument('--models', default='phb', help='phb, linear, both, or comma-separated phb,linear')
    ap.add_argument('--m-list', default='0,1,2')
    ap.add_argument('--a', type=float, default=1.0)
    ap.add_argument('--b', type=float, default=1.0)
    ap.add_argument('--R', type=float, default=1.0)
    ap.add_argument('--R1', type=float, default=0.0, help='Absolute radius of the physical input window in the LEFT mouth/funnel at z=-a. If R1>0, only 0<=r<=R1 is open at the left mouth and internal sources are forced off.')
    ap.add_argument('--R2', type=float, default=0.10, help='Absolute right-funnel radius reduction; R_right=R-R2.')
    ap.add_argument('--R3', type=float, default=0.0, help='Left mouth-side radial continuation/window scale: outer radius is R+R3 at z=-a; no detached focal-plane screen.')
    ap.add_argument('--R4', type=float, default=0.0, help='Reserved compatibility parameter; must be 0 in V43.1. No right PEC stop/support/screen is created.')
    ap.add_argument('--outer-radius', type=float, default=None, help='Legacy computational-extent hint only. It creates no material and does not define the M1 diagnostic interval in V43.1.')
    ap.add_argument('--output-r-min', type=float, default=None, help='Manual override only. Default is R_right=R-R2.')
    ap.add_argument('--output-r-max', type=float, default=None, help='Manual override only. Default is computed as R+R3.')
    ap.add_argument('--a-over-lambda', type=float, default=3.0)
    ap.add_argument('--min-a-over-lambda', type=float, default=1.0)
    ap.add_argument('--resolution', type=int, default=80)
    ap.add_argument('--dpml', type=float, default=None)
    ap.add_argument('--dpml-over-lambda', type=float, default=1.2)
    ap.add_argument('--wall-thickness', type=float, default=0.035)
    ap.add_argument('--aperture-stop-thickness', type=float, default=None)

    ap.add_argument('--source-mode', default='random', choices=[
        'random','coherent','coherent_volume','single','single_ring','ring',
        'symmetric','symmetric_even','antisymmetric','symmetric_odd',
        'mirrored_random','mirrored_random_even','mirrored_random_odd',
        'radial_sheet','axial_sheet','left_sheet','right_sheet',
        'left_input_sheet','left_input_ring','left_input_center','left_input_edge','left_input_random',
        'custom_csv'
    ])
    ap.add_argument('--source-waveform', choices=['gaussian','continuous'], default='gaussian')
    ap.add_argument('--source-duration', type=float, default=100.0, help='Finite duration for continuous waveform.')
    ap.add_argument('--source-components', default='Ez,Er,Ep')
    ap.add_argument('--nsrc', type=int, default=16)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--source-csv', default='')
    ap.add_argument('--source-r', type=float, default=None)
    ap.add_argument('--source-z', type=float, default=0.0)
    ap.add_argument('--source-pair-z', type=float, default=0.35)
    ap.add_argument('--source-r-min', type=float, default=None)
    ap.add_argument('--source-r-max', type=float, default=None)
    ap.add_argument('--source-z-min', type=float, default=None)
    ap.add_argument('--source-z-max', type=float, default=None)
    ap.add_argument('--source-size-r', type=float, default=0.0)
    ap.add_argument('--source-size-z', type=float, default=0.0)
    ap.add_argument('--input-source-offset-cells', type=float, default=3.0, help='For R1>0 input modes, place the source this many grid cells to the left of z=-a.')
    ap.add_argument('--fwidth-frac', type=float, default=0.18)
    ap.add_argument('--after-sources', type=float, default=300.0)

    ap.add_argument('--aperture-offset-cells', type=float, default=3.0, help='Automatic left/right near-monitor offset from z=±a.')
    ap.add_argument('--far-distance', type=float, default=8.0)
    ap.add_argument('--m2-distance', type=float, default=None, help='Legacy alias; overrides --far-distance when supplied.')
    ap.add_argument('--m2-capture-angle-deg', type=float, default=25.0)
    ap.add_argument('--far-zone-safety', type=float, default=1.0)
    ap.add_argument('--allow-near-m2', action='store_true')
    ap.add_argument('--left-near-z', type=float, default=None)
    ap.add_argument('--left-far-z', type=float, default=None)
    ap.add_argument('--right-near-z', type=float, default=None)
    ap.add_argument('--right-far-z', type=float, default=None)
    ap.add_argument('--axial-monitor-r-min', type=float, default=0.0)
    ap.add_argument('--axial-monitor-r-max', type=float, default=None)
    ap.add_argument('--top-r', type=float, default=None)
    ap.add_argument('--top-offset', type=float, default=0.25)
    ap.add_argument('--top-z-min', type=float, default=None)
    ap.add_argument('--top-z-max', type=float, default=None)
    ap.add_argument('--top-z-margin', type=float, default=0.25)
    ap.add_argument('--no-left-near', action='store_true')
    ap.add_argument('--no-left-far', action='store_true')
    ap.add_argument('--no-right-near', action='store_true')
    ap.add_argument('--no-right-far', action='store_true')
    ap.add_argument('--disable-m2', action='store_true', help='Legacy alias for --no-right-far.')
    ap.add_argument('--no-top-monitor', action='store_true')

    ap.add_argument('--no-profiles', action='store_true')
    ap.add_argument('--no-geometry-png', action='store_true')
    ap.add_argument('--no-archive', action='store_true')
    ap.add_argument('--skip-existing', action='store_true')
    ap.add_argument('--stop-on-error', action='store_true')
    ap.add_argument('--min-free-gb', type=float, default=1.0)
    ap.add_argument('--no-complex-fields', action='store_true')
    ap.add_argument('--no-npz', action='store_true')

    ap.add_argument('--narrow-theta95-deg', type=float, default=25.0)
    ap.add_argument('--useful-flux-min', type=float, default=1e-12)
    ap.add_argument('--min-output-cells', type=float, default=8.0)
    ap.add_argument('--min-wall-cells', type=float, default=4.0)
    ap.add_argument('--min-extension-cells', type=float, default=8.0)
    ap.add_argument('--min-pml-cells', type=float, default=12.0)
    ap.add_argument('--allow-underresolved', action='store_true', help='Accepted for compatibility; V42 always warns and continues.')

    ap.add_argument('--harminv-component', default='Ez,Er,Ep', help='One component or a comma-separated list, e.g. Ez,Er,Ep. Multiple components reduce nodal blind spots in modal search.')
    ap.add_argument('--harminv-points', default='')
    ap.add_argument('--modal-after-sources', type=float, default=400.0)
    ap.add_argument('--modal-fwidth-frac', type=float, default=0.70)
    ap.add_argument('--verify-frequency', type=float, default=0.0)
    ap.add_argument('--verify-fwidth-frac', type=float, default=0.025)
    ap.add_argument('--m1-subwindows', default='')

    ap.add_argument('--no-internal-diagnostics', action='store_true')
    ap.add_argument('--no-full-map-csv', action='store_true')
    ap.add_argument('--internal-r-min', type=float, default=None)
    ap.add_argument('--internal-r-max', type=float, default=None)
    ap.add_argument('--internal-z-min', type=float, default=None)
    ap.add_argument('--internal-z-max', type=float, default=None)
    ap.add_argument('--map-domain', choices=['resonator','full_cell','custom'], default='resonator')
    ap.add_argument('--map-r-margin', type=float, default=0.5)
    ap.add_argument('--map-z-margin', type=float, default=0.5)
    ap.add_argument('--map-stride', type=int, default=1)
    ap.add_argument('--exclude-metal-buffer-cells', type=float, default=2.0)
    ap.add_argument('--eq-r-min', type=float, default=None)
    ap.add_argument('--eq-r-max', type=float, default=None)
    ap.add_argument('--eq-z-min', type=float, default=None)
    ap.add_argument('--eq-z-max', type=float, default=None)
    ap.add_argument('--axis-r-max', type=float, default=0.25)
    ap.add_argument('--outer-r-min', type=float, default=None)
    ap.add_argument('--focal-band-halfwidth', type=float, default=0.05)
    ap.add_argument('--no-density-maps', action='store_true')
    ap.add_argument('--density-map-full-section', action='store_true')
    ap.add_argument('--different-model-seeds', action='store_true')

    ap.add_argument('--progress-interval', type=float, default=10.0, help='MEEP-time interval between live terminal progress messages.')
    ap.add_argument('--meep-verbosity', type=int, default=1)
    return ap


def main_v42(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser_v42().parse_args(argv)
    models = parse_models_v42(args.models)
    ms = parse_int_list(args.m_list)
    far_distance = float(args.m2_distance) if args.m2_distance is not None else float(args.far_distance)
    cfg = UniversalRunConfig(
        stage=args.stage, outroot=args.outroot, models=models, m_list=ms,
        a_over_lambda=args.a_over_lambda, min_a_over_lambda=args.min_a_over_lambda,
        resolution=args.resolution, dpml_over_lambda=args.dpml_over_lambda, dpml=args.dpml,
        wall_thickness=args.wall_thickness, aperture_stop_thickness=args.aperture_stop_thickness,
        source_mode=args.source_mode, source_waveform=args.source_waveform, source_duration=args.source_duration,
        source_components=args.source_components, nsrc=args.nsrc, seed=args.seed, source_csv=args.source_csv,
        source_r=args.source_r, source_z=args.source_z, source_pair_z=args.source_pair_z,
        source_r_min=args.source_r_min, source_r_max=args.source_r_max,
        source_z_min=args.source_z_min, source_z_max=args.source_z_max,
        source_size_r=args.source_size_r, source_size_z=args.source_size_z,
        input_source_offset_cells=args.input_source_offset_cells,
        fwidth_frac=args.fwidth_frac, after_sources=args.after_sources,
        aperture_offset_cells=args.aperture_offset_cells,
        enable_m2=not (args.disable_m2 or args.no_right_far),
        m2_distance=far_distance, m2_capture_angle_deg=args.m2_capture_angle_deg,
        far_zone_safety=args.far_zone_safety, allow_near_m2=args.allow_near_m2,
        enable_left_near=not args.no_left_near, enable_left_far=not args.no_left_far,
        enable_right_near=not args.no_right_near, enable_right_far=not (args.no_right_far or args.disable_m2),
        enable_top=not args.no_top_monitor,
        left_near_z=args.left_near_z, left_far_z=args.left_far_z,
        right_near_z=args.right_near_z, right_far_z=args.right_far_z,
        far_distance=far_distance, axial_monitor_r_min=args.axial_monitor_r_min,
        axial_monitor_r_max=args.axial_monitor_r_max, top_r=args.top_r, top_offset=args.top_offset,
        top_z_min=args.top_z_min, top_z_max=args.top_z_max, top_z_margin=args.top_z_margin,
        save_profiles=not args.no_profiles, save_geometry_png=not args.no_geometry_png,
        archive=not args.no_archive, skip_existing=args.skip_existing, stop_on_error=args.stop_on_error,
        min_free_gb=args.min_free_gb, save_complex_fields=not args.no_complex_fields, save_npz=not args.no_npz,
        narrow_theta95_deg=args.narrow_theta95_deg, useful_flux_min=args.useful_flux_min,
        min_output_cells=args.min_output_cells, min_wall_cells=args.min_wall_cells,
        min_extension_cells=args.min_extension_cells, min_pml_cells=args.min_pml_cells,
        allow_underresolved=args.allow_underresolved,
        harminv_component=args.harminv_component, harminv_points=args.harminv_points,
        modal_after_sources=args.modal_after_sources, modal_fwidth_frac=args.modal_fwidth_frac,
        verify_frequency=args.verify_frequency, verify_fwidth_frac=args.verify_fwidth_frac,
        m1_subwindows=args.m1_subwindows,
        enable_internal_diagnostics=not args.no_internal_diagnostics,
        save_internal_map=not args.no_full_map_csv,
        internal_r_min=args.internal_r_min, internal_r_max=args.internal_r_max,
        internal_z_min=args.internal_z_min, internal_z_max=args.internal_z_max,
        map_domain=args.map_domain, map_r_margin=args.map_r_margin, map_z_margin=args.map_z_margin,
        map_stride=max(1,args.map_stride), exclude_metal_buffer_cells=args.exclude_metal_buffer_cells,
        eq_r_min=args.eq_r_min, eq_r_max=args.eq_r_max, eq_z_min=args.eq_z_min, eq_z_max=args.eq_z_max,
        axis_r_max=args.axis_r_max, outer_r_min=args.outer_r_min,
        focal_band_halfwidth=args.focal_band_halfwidth,
        save_density_maps=not args.no_density_maps,
        density_map_full_section=args.density_map_full_section,
        match_model_seeds=not args.different_model_seeds,
        progress_interval=args.progress_interval, meep_verbosity=args.meep_verbosity,
    )
    outroot = ensure_dir(Path(cfg.outroot))
    write_methodology_v42(outroot,args,cfg)
    print('='*76, flush=True)
    print('PHB V43.1 MEEP VERIFIER: CONNECTED LEFT WALL + FULLY OPEN RIGHT HALF-SPACE', flush=True)
    print(f'outroot={outroot}', flush=True)
    print(f'stage={cfg.stage}; models={models}; m={ms}; source={cfg.source_mode}/{cfg.source_waveform}', flush=True)
    print('LIVE PROGRESS is enabled. Use python -u and Docker -it; another terminal may use docker logs -f.', flush=True)
    print('LOW CELL COUNTS NEVER STOP THE RUN: warnings are printed and saved.', flush=True)
    print('='*76, flush=True)

    if cfg.stage == 'summary':
        rows = collect_metrics(outroot)
        write_summary(outroot,rows)
        if cfg.archive:
            z = archive_folder(outroot)
            if z: print('Archive written:',z,flush=True)
        return 0

    for model in models:
        for m in ms:
            run_one_case_v42(args,cfg,model,m)
            write_summary(outroot,collect_metrics(outroot))
    all_rows = collect_metrics(outroot)
    write_summary(outroot,all_rows)
    if cfg.archive:
        z = archive_folder(outroot)
        if z: print('Archive written:',z,flush=True)
    bad = [r for r in all_rows if str(r.get('candidate_status','')).startswith(('ERROR','INVALID'))]
    return 1 if bad else 0


# Public entry point: external wrappers importing `main` receive the corrected
# V43.1 open-right implementation rather than the retained legacy compatibility code.
main = main_v42

if __name__ == '__main__':
    raise SystemExit(main_v42())

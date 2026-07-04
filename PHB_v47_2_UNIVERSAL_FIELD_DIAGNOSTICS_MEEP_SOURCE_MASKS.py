#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHB_v47_2_UNIVERSAL_FIELD_DIAGNOSTICS_MEEP_SOURCE_MASKS.py

Первый лёгкий MEEP/FDTD-скрипт для верификации гипотезы PHB.

Версия v47.1 UNIVERSAL FIELD DIAGNOSTICS / CENTER-SHIFTED R2 WINDOW:
  - v47.1: добавлено сохранение координат и параметров источников; добавлены метрики расстояния
    от максимумов |E|² на field-map картах до ближайшего источника для отделения source-peak
    артефактов от PHB/щелевых полевых концентраций;
  - v47.2: сохраняется material/air/slot mask около щели и в field-map областях;
  - v47.2: сохраняются source-excluded |E|²-карты и метрики пиков вне источников;
  - v47.2: добавлен M2/M1 warning, если дальний монитор собирает не только поток выбранной M1-апертуры;
  - v47: добавлены 2D DFT field-map diagnostics для внутренних карт поля,
    near-slot/pre-slot/post-slot/focal-corridor областей и Poynting-derived метрик;
  - sim.reset_meep() перенесён в finally;
  - v46: M1 и M2 лежат на общей продольной линии, параллельной F−→F+ / +z,
    с общим радиальным центром r=R+shift;
  - v46: R2 трактуется как полная абсолютная ширина окна/монитора W, не R*R2;
  - v46: выходное окно задаётся универсально как прямоугольное правое окно
    в вертикальной выходной стенке z≈+a, центрированное в точке r=R+shift:
        [R+shift-W/2, R+shift+W/2].
    shift может быть положительным (выше фокальной оси r=R) или отрицательным
    (ниже, к оси вращения). shift=0 возвращает симметричное окно v46.3.
    старый режим выреза вдоль гиперболической кривой оставлен только как --window-model horn_cut;
  - удалён мёртвый флаг save_dft_fields;
  - kappa_R2_escape_proxy теперь использует MEEP flux-monitor aperture_R2_z как основной числитель,
    а локальный DFT/Poynting flux только как fallback;
  - mp.get_fluxes(obj) сохранён намеренно: это штатный MEEP Python API;
  - добавлены правильные publication-рисунки: осесимметричный конус вокруг F−→F+ / +z;
  - добавлены theta_max_sampled/theta_at_max_density и явное поле theta_from_z_axis_deg;
  - optional M2 far monitor исключён из kappa, чтобы не двойно учитывать диагностический срез.

Базовая задача:
    гиперболическая PHB-геометрия -> low-m / low-J фокально-цилиндрический канал
    -> селективный вывод через R2 -> направленный выход энергии из щели.

Главная задача v46:
    НЕ доказывать всю лазерную физику резонатора сразу, а быстро найти или не найти
    необычные кандидаты по локальному углу выхода энергии из R2-щели.

Что принципиально исправлено относительно v41-оркестратора:
  1) Нет base64-монолитов. Вся физика и вся логика видны в этом файле.
  2) H-поля включены обязательно: локальный угол выхода считается из вектора Пойнтинга.
  3) R2 везде означает абсолютную ширину выходного окна W; это не R*R2.
  4) Главный масштаб задаётся через a/lambda; жёстко требуется a >= lambda.
  5) Главный монитор — локальный апертурный монитор M1 в исправленном R2-окне:
         r in [R+shift-W/2, R+shift+W/2], центр r=R+shift.
         v46 default: окно является правым вертикальным прямоугольным окном
         в стенке z≈+a с задаваемым сдвигом центра относительно r=R.
     Это именно монитор выхода из щели, а не дальний монитор.
  6) Для закрытых вариантов cylinder и halfring не включается дополнительный второй
     диагностический выходной монитор. Первый этап проверяет только R2-щель.
     Если M2 включён, он также центрируется на r=R и параллелен M1.
  7) Дальний монитор не используется по умолчанию. Если включить --enable-far-monitor,
     скрипт проверяет условие дальней зоны Фраунгофера и не даст молча сравнить ближнее
     поле с дальнепольной аналитикой.
  8) После завершения создаётся ZIP-архив всей папки результатов рядом с ней.
  9) Логи пишутся также в JSONL, чтобы потом строить таблицы и графики по сотням прогонов.

Запуск из WSL-папки /home/qwerty/meep-working:

  docker run --rm -it -v "$PWD":/work -e MPLBACKEND=Agg meep-working \
    python /work/PHB_v47_UNIVERSAL_FIELD_DIAGNOSTICS_MEEP.py --stage plan

Первый быстрый поиск необычного угла:

  docker run --rm -it -v "$PWD":/work -e MPLBACKEND=Agg meep-working \
    python /work/PHB_v47_UNIVERSAL_FIELD_DIAGNOSTICS_MEEP.py \
      --stage fast --outroot /work/PHB_V46_1_FIRST_FAST \
      --a 1 --b 1 --R 1 --a-over-lambda 1 \
      --R2-list 0.05,0.07,0.10,0.20,0.30 \
      --m-list 0,1,2 --phb-types open,cylinder,halfring \
      --resolution 64 --after-sources 120 --nsrc 12 --skip-existing

Более надёжный confirm-кандидат после fast:

  docker run --rm -it -v "$PWD":/work -e MPLBACKEND=Agg meep-working \
    python /work/PHB_v47_UNIVERSAL_FIELD_DIAGNOSTICS_MEEP.py \
      --stage confirm --outroot /work/PHB_V46_1_CONFIRM \
      --a 1 --b 1 --R 1 --a-over-lambda 3 \
      --R2-list 0.05,0.07,0.10,0.20,0.30 \
      --m-list 0,1,2 --phb-types open,cylinder,halfring \
      --resolution 96 --after-sources 240 --nsrc 24 --skip-existing
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
except Exception:  # plan/geometry stages can still run without meep
    mp = None


# -----------------------------
# Data classes
# -----------------------------

@dataclass(frozen=True)
class PHBGeometry:
    a: float = 1.0
    b: float = 1.0
    R: float = 1.0
    R2: float = 0.10  # absolute full output-window width W, not R*R2
    window_offset: float = 0.0  # shift relative to focal-ring axis r=R; interpretation controlled by window_offset_mode
    window_offset_mode: str = "center"  # "center": offset is center shift; "lower_edge": offset is lower-edge shift

    @property
    def c(self) -> float:
        return math.sqrt(self.a * self.a + self.b * self.b)

    @property
    def L(self) -> float:
        return self.a * math.sqrt(1.0 + (self.R / self.b) ** 2)

    @property
    def lambert_note(self) -> str:
        return "R2 is absolute output-window width W; not R*R2"

    def rho_horn(self, abs_z: float) -> float:
        if abs_z < self.a:
            return float("nan")
        val = self.R - self.b * math.sqrt(max((abs_z / self.a) ** 2 - 1.0, 0.0))
        return max(0.0, val)

    def z_for_rho(self, rho: float) -> float:
        # right horn coordinate for a given horn radius rho <= R
        return self.a * math.sqrt(1.0 + ((self.R - rho) / self.b) ** 2)

    # ------------------------------------------------------------
    # Universal center-shifted output-window / monitor geometry, v46
    # ------------------------------------------------------------
    # R2 is the absolute full width W of the selected right-side output window.
    # ``window_offset`` is now a CENTER SHIFT relative to the focal-ring axis r=R:
    #
    #         M1 interval = [R + shift - W/2, R + shift + W/2]
    #         M1 center   =  R + shift
    #
    # shift > 0 moves the aperture upward/outward; shift < 0 moves it downward
    # toward the mechanical rotation axis.  shift = 0 reproduces the v46.3
    # symmetric aperture around the external-focal-axis radius r=R.
    #
    # The cylindrical MEEP coordinate r cannot be negative; validation rejects
    # cases where the lower edge would cross r=0 instead of silently clipping.
    #
    def window_r_center(self, phb_type: str) -> float:
        if str(self.window_offset_mode).lower() in ("lower", "lower_edge", "lower-edge", "edge"):
            return self.R + self.window_offset + 0.5 * self.R2
        return self.R + self.window_offset

    def window_r_low_raw(self, phb_type: str) -> float:
        if str(self.window_offset_mode).lower() in ("lower", "lower_edge", "lower-edge", "edge"):
            return self.R + self.window_offset
        return self.window_r_center(phb_type) - 0.5 * self.R2

    def window_r_high_raw(self, phb_type: str) -> float:
        if str(self.window_offset_mode).lower() in ("lower", "lower_edge", "lower-edge", "edge"):
            return self.R + self.window_offset + self.R2
        return self.window_r_center(phb_type) + 0.5 * self.R2

    def window_r_low(self, phb_type: str) -> float:
        return self.window_r_low_raw(phb_type)

    def window_r_high(self, phb_type: str) -> float:
        return self.window_r_high_raw(phb_type)

    def window_width(self, phb_type: str) -> float:
        return max(0.0, self.window_r_high(phb_type) - self.window_r_low(phb_type))

    def z_window_horn_end(self, phb_type: str) -> float:
        # Legacy horn_cut helper only.  In v46 vertical_shift mode
        # the hyperbolic horn is not cut by the aperture model.
        return self.z_for_rho(self.R)

    # Legacy properties retained only for old helper code; in v46 legacy they
    # also return the symmetric focal-axis monitor interval.
    @property
    def slot_r_inner(self) -> float:
        return self.window_r_low("open")

    @property
    def slot_r_outer(self) -> float:
        return self.window_r_high("open")

    @property
    def slot_r_center(self) -> float:
        return self.window_r_center("open")

    @property
    def z_slot_exact_end(self) -> float:
        return self.z_window_horn_end("open")


@dataclass
class RunConfig:
    stage: str
    outroot: str
    phb_type: str
    m: int
    resolution: int
    a_over_lambda: float
    min_a_over_lambda: float
    dpml_over_lambda: float
    dpml: Optional[float]
    wall_thickness: float
    window_model: str
    window_offset: float
    window_offset_mode: str
    vertical_window_thickness: Optional[float]
    source_components: str
    source_mode: str
    nsrc: int
    seed: int
    fwidth_frac: float
    after_sources: float
    skip_existing: bool
    allow_underresolved: bool
    min_radial_slot_cells: float
    min_axial_slot_cells: float
    aperture_offset_cells: float
    enable_far_monitor: bool
    far_distance: float
    far_capture_angle_deg: float
    far_monitor_mode: str
    far_zone_safety: float
    allow_near_field_far_monitor: bool
    narrow_theta95_deg: float
    narrow_lobe5_min: float
    useful_flux_min: float
    kappa_min: float
    min_free_gb: float
    archive: bool
    stop_on_error: bool
    # v47/v46.1 full-wave field-map diagnostics.
    save_field_maps: bool
    field_map_regions: str
    field_map_components: str
    field_map_r_pad: float
    field_map_z_half: float
    field_map_corridor_half_width: float
    field_map_interior_rmax: Optional[float]
    field_map_max_points_per_region: int
    field_map_csv_stride: int
    field_map_png: bool
    field_map_npz: bool
    # v47.2 artefact-control diagnostics.
    save_material_maps: bool
    save_source_excluded_maps: bool
    source_exclusion_cells: float
    m2_extra_flux_warning_ratio: float


# -----------------------------
# Generic helpers
# -----------------------------

def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def parse_float_list(s: str) -> List[float]:
    out: List[float] = []
    for part in str(s).replace(";", ",").split(","):
        part = part.strip()
        if part:
            out.append(float(part))
    return out


def parse_int_list(s: str) -> List[int]:
    out: List[int] = []
    for part in str(s).replace(";", ",").split(","):
        part = part.strip()
        if part:
            out.append(int(part))
    return out


def parse_phb_types(s: str) -> List[str]:
    aliases = {
        "open_waist": "open",
        "openwaist": "open",
        "none": "open",
        "no_wall": "open",
        "half_ring": "halfring",
        "half-torus": "halfring",
        "half_torus": "halfring",
        "semitorus": "halfring",
        "toroid": "halfring",
        "toroidal": "halfring",
        "тороид": "halfring",
        "toroidvlnyi": "halfring",
    }
    out: List[str] = []
    for part in str(s).replace(";", ",").split(","):
        v = part.strip().lower().replace("-", "_")
        if not v:
            continue
        v = aliases.get(v, v)
        if v not in ("open", "cylinder", "halfring"):
            raise ValueError(f"Unknown PHB type: {part!r}. Use open,cylinder,halfring")
        if v not in out:
            out.append(v)
    return out


def safe_float(x, default: float = float("nan")) -> float:
    try:
        return float(x)
    except Exception:
        return default


def weighted_percentile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return float("nan")
    values = values[mask]
    weights = weights[mask]
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    c = np.cumsum(w)
    target = float(q) * c[-1]
    idx = int(np.searchsorted(c, target, side="left"))
    idx = min(max(idx, 0), len(v) - 1)
    return float(v[idx])


def weighted_circular_phase_rms(phase: np.ndarray, weights: np.ndarray) -> float:
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


def write_json(path: Path, obj: object) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def append_jsonl(path: Path, obj: object) -> None:
    ensure_dir(path.parent)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(obj, ensure_ascii=False, sort_keys=True) + "\n")


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: Optional[List[str]] = None) -> None:
    ensure_dir(path.parent)
    if fieldnames is None:
        keys: List[str] = []
        for r in rows:
            for k in r.keys():
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def free_gb(path: Path) -> float:
    usage = shutil.disk_usage(str(path if path.exists() else path.parent))
    return usage.free / (1024 ** 3)


def archive_folder(folder: Path) -> Optional[Path]:
    if not folder.exists():
        return None
    zip_base = folder.with_suffix("")
    zip_path = folder.parent / f"{folder.name}.zip"
    if zip_path.exists():
        zip_path.unlink()
    shutil.make_archive(str(zip_base), "zip", root_dir=str(folder.parent), base_dir=folder.name)
    return zip_path


# -----------------------------
# PHB geometry and validation
# -----------------------------

def wavelength_from_a(a: float, a_over_lambda: float) -> float:
    return float(a) / float(a_over_lambda)


def make_wave_report(g0: PHBGeometry, args) -> Dict[str, float | str]:
    lam = wavelength_from_a(g0.a, args.a_over_lambda)
    return {
        "lambda": lam,
        "a_over_lambda": g0.a / lam,
        "b_over_lambda": g0.b / lam,
        "R_over_lambda": g0.R / lam,
        "two_a_over_lambda": 2.0 * g0.a / lam,
        "L_over_lambda": g0.L / lam,
        "cells_per_wavelength": args.resolution * lam,
        "a_scale_status": "OK_A_GE_LAMBDA" if g0.a / lam >= args.min_a_over_lambda else "INVALID_A_LT_LAMBDA",
    }


def effective_slot_interval(g: PHBGeometry, cfg: RunConfig) -> Tuple[float, float, float, bool]:
    """Return z0,z1_effective,z1_exact,widened for the R2 output opening.

    v46 default model: ``vertical_shift``.
        The output window is a rectangular right-side aperture in the local
        vertical output wall z≈+a.  Its center is shifted relative to the
        focal-ring axis rho=R and spans r in [R+shift-W/2, R+shift+W/2].
        The z-thickness of this rectangular opening is a physical/numerical
        wall thickness, not a hyperbolic-curve widening.  Therefore
        ``slot_widened`` remains False.

    Legacy model: ``horn_cut``.
        Preserves the old v46 legacy behavior where the lower half of the window
        follows the hyperbolic horn from rho=R down to rho=R-R2/2 and may be
        numerically widened if under-resolved.  This is kept only for
        comparison and is not the default PHB working model.
    """
    model = getattr(cfg, "window_model", "vertical_shift")
    if model == "vertical_shift":
        # Thickness along z of the right output wall.  If the user does not
        # provide a physical thickness, use enough cells for FDTD to see the
        # rectangular aperture.  This is not artificial horn widening; it is
        # the explicitly modelled thickness of the output wall/aperture plane.
        min_t = float(cfg.min_axial_slot_cells) / max(float(cfg.resolution), 1.0)
        if getattr(cfg, "vertical_window_thickness", None) is not None:
            t = float(cfg.vertical_window_thickness)
        else:
            t = max(float(cfg.wall_thickness), min_t)
        t = max(t, 1.0 / max(float(cfg.resolution), 1.0))
        z0 = float(g.a) - 0.5 * t
        z1 = float(g.a) + 0.5 * t
        return z0, z1, z1, False

    # Legacy horn_cut model from v46 legacy.
    z0 = float(g.a)
    z1_exact = min(float(g.z_window_horn_end(cfg.phb_type)), float(g.L))
    dz_exact = max(0.0, z1_exact - z0)
    min_dz = float(cfg.min_axial_slot_cells) / max(float(cfg.resolution), 1.0)
    widened = False
    z1_eff = z1_exact
    if dz_exact < min_dz:
        z1_eff = min(float(g.L), z0 + min_dz)
        widened = True
    return z0, z1_eff, z1_exact, widened

def phb_central_wall_radius(z: float, g: PHBGeometry, phb_type: str) -> float:
    """Central closure radius for the corrected v46 geometry."""
    if phb_type == "cylinder":
        # User correction: cylinder wall starts from the end of the new insert,
        # i.e. from R + W/2, not from R.
        return g.R + 0.5 * g.R2
    if phb_type == "halfring":
        return g.R + math.sqrt(max(g.a * g.a - z * z, 0.0))
    return g.R


def is_right_upper_window_on_halfring(z: float, g: PHBGeometry, cfg: RunConfig, wall_t: float) -> bool:
    """True where the right halfring closure segment must be removed.

    The corrected closed-type aperture is symmetric around rho=R:
        [R-W/2, R+W/2].
    The lower half is on the hyperbolic horn.  The upper half, for halfring,
    lies on the small right part of the halfring closure where its radius
    runs from R to R+W/2.  That wall segment must be absent.
    """
    if cfg.phb_type != "halfring" or z <= 0.0 or abs(z) >= g.a:
        return False
    rwc = phb_central_wall_radius(z, g, cfg.phb_type)
    return rwc <= g.R + 0.5 * g.R2 + 0.50 * wall_t


def validate_case(g: PHBGeometry, cfg: RunConfig) -> Dict[str, object]:
    problems: List[str] = []
    warnings: List[str] = []
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    if cfg.a_over_lambda < cfg.min_a_over_lambda:
        problems.append(f"a/lambda={cfg.a_over_lambda:g} < required {cfg.min_a_over_lambda:g}")
    if g.a <= 0 or g.b <= 0 or g.R <= 0:
        problems.append("a,b,R must be positive")
    if g.R2 <= 0:
        problems.append("R2 must be positive for aperture verification")
    # v46: negative shifts are allowed physically as a downward/toward-axis displacement,
    # but the cylindrical MEEP radius cannot cross r=0.
    if g.window_r_low(cfg.phb_type) < 0.0:
        problems.append(
            f"window lower edge r_min={g.window_r_low(cfg.phb_type):g} < 0; "
            "reduce negative --window-shift or use a narrower R2"
        )
    if g.window_r_center(cfg.phb_type) < 0.0:
        problems.append(f"window center R+shift={g.window_r_center(cfg.phb_type):g} must be >= 0")
    if g.R2 > max(g.R, 1e-300):
        warnings.append(f"R2={g.R2:g} is larger than R={g.R:g}; R2 is absolute window width W")
    if cfg.phb_type == "halfring" and (g.window_r_high(cfg.phb_type) - g.R) > g.a:
        warnings.append(
            f"upper edge above focal radius is {g.window_r_high(cfg.phb_type)-g.R:g}, exceeding halfring radial span a={g.a:g}; "
            "the selected shifted band may lie outside the halfring closure region"
        )
    radial_cells = g.window_width(cfg.phb_type) * cfg.resolution
    z0, z1_eff, z1_exact, widened = effective_slot_interval(g, cfg)
    exact_axial_cells = max(0.0, z1_exact - z0) * cfg.resolution
    eff_axial_cells = max(0.0, z1_eff - z0) * cfg.resolution
    if radial_cells < cfg.min_radial_slot_cells and not cfg.allow_underresolved:
        problems.append(
            f"corrected aperture radial interval has {radial_cells:.3g} cells < {cfg.min_radial_slot_cells:g}; "
            f"increase --resolution or use larger R2"
        )
    if eff_axial_cells < cfg.min_axial_slot_cells and not cfg.allow_underresolved:
        problems.append(
            f"effective axial slot has {eff_axial_cells:.3g} cells < {cfg.min_axial_slot_cells:g}; "
            f"increase --resolution"
        )
    if widened:
        warnings.append(
            f"exact hyperbolic slot axial width is {exact_axial_cells:.3g} cells; "
            f"numerical slot widened to {eff_axial_cells:.3g} cells for fast screening"
        )
    if getattr(cfg, "window_model", "vertical_shift") == "vertical_shift":
        warnings.append(
            "vertical_shift window model: R2 is a right aperture [R+shift-W/2,R+shift+W/2] in z≈+a wall; "
            "right horn is not cut by the aperture model; axial cells describe aperture-plane thickness"
        )
    if cfg.enable_far_monitor:
        # Far-zone diagnostic.  It is a warning by default, not a hard stop,
        # because the user may intentionally use a shorter diagnostic distance.
        D = 2.0 * g.R
        fraunhofer = 2.0 * D * D / max(lam, 1e-300)
        if cfg.far_distance < cfg.far_zone_safety * fraunhofer and not cfg.allow_near_field_far_monitor:
            warnings.append(
                f"far monitor is closer than Fraunhofer estimate: distance={cfg.far_distance:g}, "
                f"recommended >= safety*2D^2/lambda = {cfg.far_zone_safety * fraunhofer:g}. "
                f"This run is still allowed; interpret far metrics as near/intermediate-field diagnostics."
            )
    return {
        "ok": not problems,
        "problems": problems,
        "warnings": warnings,
        "lambda": lam,
        "slot_z0": z0,
        "slot_z1_effective": z1_eff,
        "slot_z1_exact": z1_exact,
        "slot_widened": widened,
        "slot_exact_axial_cells": exact_axial_cells,
        "slot_effective_axial_cells": eff_axial_cells,
        "slot_radial_cells": radial_cells,
        "window_model": getattr(cfg, "window_model", "vertical_shift"),
        "vertical_window_thickness_used": max(0.0, z1_eff - z0),
        "vertical_window_center_z": 0.5 * (z0 + z1_eff),
    }


# -----------------------------
# Plot geometry
# -----------------------------

def plot_geometry(out_png: Path, g: PHBGeometry, cfg: RunConfig) -> None:
    if plt is None:
        return
    ensure_dir(out_png.parent)
    z = np.linspace(-g.L, g.L, 1200)
    rwall = np.full_like(z, np.nan, dtype=float)
    for i, zz in enumerate(z):
        az = abs(float(zz))
        if az >= g.a:
            rwall[i] = g.rho_horn(az)
        else:
            if cfg.phb_type == "open":
                rwall[i] = np.nan
            else:
                rwall[i] = phb_central_wall_radius(float(zz), g, cfg.phb_type)

    z0, z1_eff, z1_exact, widened = effective_slot_interval(g, cfg)
    z_ap = z1_eff + cfg.aperture_offset_cells / max(cfg.resolution, 1)
    r0, r1, rc = g.window_r_low(cfg.phb_type), g.window_r_high(cfg.phb_type), g.window_r_center(cfg.phb_type)

    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    ax.plot(z, rwall, lw=2, label="PHB reflecting wall / closure")
    ax.plot(z, -rwall, lw=2)
    if getattr(cfg, "window_model", "vertical_shift") == "vertical_shift":
        # v46: right rectangular aperture in the vertical output wall z≈+a.
        ax.fill_between([z0, z1_eff], [r0, r0], [r1, r1], alpha=0.28, label="right shifted vertical R2 window [R+shift-W/2,R+shift+W/2]")
        ax.fill_between([z0, z1_eff], [-r1, -r1], [-r0, -r0], alpha=0.18)
        ax.plot([g.a, g.a], [r0, r1], lw=4, alpha=0.85, label="center-shifted R2 opening: [R+shift-W/2, R+shift+W/2]")
        ax.plot([g.a, g.a], [-r1, -r0], lw=4, alpha=0.55)
    else:
        # Legacy horn_cut visualisation.
        zzs = np.linspace(z0, z1_eff, 80)
        rs = np.array([g.rho_horn(float(x)) for x in zzs])
        ax.plot(zzs, rs, lw=4, alpha=0.7, label="removed lower horn part of corrected window")
        ax.plot(zzs, -rs, lw=4, alpha=0.7)
        if cfg.phb_type in ("cylinder", "halfring"):
            ax.plot([g.a, g.a], [g.R, r1], lw=4, alpha=0.7, label="removed upper part of symmetric window")
            ax.plot([g.a, g.a], [-r1, -g.R], lw=4, alpha=0.7)
        if cfg.phb_type == "open":
            ax.plot([-g.a, -g.a], [g.R, g.R + 0.5 * g.R2], lw=4, alpha=0.7, label="left vertical add-on wall")
            ax.plot([-g.a, -g.a], [-(g.R + 0.5 * g.R2), -g.R], lw=4, alpha=0.7)

    # M1 aperture monitor: immediately after the corrected output window.
    ax.plot([z_ap, z_ap], [r0, r1], lw=4, label="M1 aperture monitor")
    ax.plot([z_ap, z_ap], [-r1, -r0], lw=4)

    # Optional M2 far monitor: parallel to M1, placed at user-selected distance.
    if cfg.enable_far_monitor:
        z_far_plot = z_ap + cfg.far_distance
        width_m1 = max(g.window_width(cfg.phb_type), 1.0 / max(cfg.resolution, 1))
        capture_half = max(0.5 * width_m1, cfg.far_distance * math.tan(math.radians(cfg.far_capture_angle_deg)))
        # v46: M2 follows the same shifted M1 window center.
        ref_r = g.window_r_center(cfg.phb_type)
        fr0 = max(0.0, ref_r - capture_half)
        fr1 = ref_r + capture_half
        ax.plot([z_far_plot, z_far_plot], [fr0, fr1], lw=3, ls="--", label="M2 far monitor")
        ax.plot([z_far_plot, z_far_plot], [-fr1, -fr0], lw=3, ls="--")

    ax.scatter([+g.c, -g.c], [g.R, g.R], marker="x", s=60, label="external focal rings in meridional section")
    ax.scatter([+g.c, -g.c], [-g.R, -g.R], marker="x", s=60)
    ax.axhline(g.R, ls="--", lw=0.8)
    ax.axhline(-g.R, ls="--", lw=0.8)
    ax.axvline(g.a, ls=":", lw=0.8)
    ax.axvline(-g.a, ls=":", lw=0.8)
    ax.set_xlabel("axial coordinate z")
    ax.set_ylabel("radius r in meridional drawing")
    ax.set_title(
        f"PHB v46 geometry: type={cfg.phb_type}, a={g.a:g}, b={g.b:g}, R={g.R:g}, R2=W={g.R2:g}, "
        f"shifted M1 center=R+shift={rc:g}, shift={g.window_offset:g}, W={g.R2:g}"
    )
    m2_note = "off"
    if cfg.enable_far_monitor:
        m2_note = f"z={z_ap + cfg.far_distance:.4g}, distance={cfg.far_distance:.4g}, mode={cfg.far_monitor_mode}"
    note = (
        f"window model={getattr(cfg, 'window_model', 'vertical_shift')}; z=[{z0:.4g},{z1_eff:.4g}], widened={widened}\n"
        f"center-shifted M1 aperture monitor: z={z_ap:.4g}, r in [{r0:.4g},{r1:.4g}], focal r=R={g.R:.4g}, center={rc:.4g}\n"
        f"M2 far monitor: {m2_note}; H-fields required for Poynting angle"
    )
    ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=8, va="bottom")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


# -----------------------------
# MEEP model
# -----------------------------

def require_meep() -> None:
    if mp is None:
        raise RuntimeError("meep is not available. Run this script inside the meep-working Docker image for FDTD stages.")


def effective_wall_t(cfg: RunConfig, lam: float) -> float:
    # Keep wall at least several grid cells and not too thin compared with wavelength.
    return max(float(cfg.wall_thickness), 2.0 / max(cfg.resolution, 1), 0.04 * lam)


def make_material_function(g: PHBGeometry, cfg: RunConfig):
    require_meep()
    air = mp.Medium(epsilon=1.0)
    metal = getattr(mp, "metal", mp.Medium(epsilon=1.0e9))
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    wall_t = effective_wall_t(cfg, lam)
    z0_slot, z1_slot, _, _ = effective_slot_interval(g, cfg)

    def matfun(pos):
        r = float(pos.x)
        z = float(pos.z)
        az = abs(z)
        if r < 0:
            return air
        if az > g.L:
            return air

        # v46 default: the output aperture is a right-side rectangular window
        # in the local vertical output wall z≈+a, with a user-selected CENTER
        # shift relative to the focal-ring axis rho=R.
        if getattr(cfg, "window_model", "vertical_shift") == "vertical_shift" and g.R2 > 0.0 and z > 0.0:
            r0_win = g.window_r_low(cfg.phb_type)
            r1_win = g.window_r_high(cfg.phb_type)
            if z0_slot <= z <= z1_slot and (r0_win - 0.50 * wall_t) <= r <= (r1_win + 0.50 * wall_t):
                return air

        # v46 open type: optional left vertical insert aligned with the selected
        # shifted window band.  This keeps the open reference geometry symmetric
        # with the selected aperture band instead of hard-coding only the upper half.
        if cfg.phb_type == "open" and z < 0.0:
            r0_win = max(0.0, g.window_r_low(cfg.phb_type))
            r1_win = g.window_r_high(cfg.phb_type)
            if abs(z + g.a) <= 0.65 * wall_t and (r0_win - wall_t) <= r <= (r1_win + wall_t):
                return metal

        # central closure: none/open, cylinder, or halfring/half-toroidal meridional closure
        if az < g.a:
            if cfg.phb_type == "open":
                return air
            rwc = phb_central_wall_radius(z, g, cfg.phb_type)

            # v46 halfring correction: remove the right upper window segment
            # where the halfring closure radius is between R and R+W/2.
            if getattr(cfg, "window_model", "vertical_shift") == "horn_cut" and is_right_upper_window_on_halfring(z, g, cfg, wall_t):
                return air

            if max(0.0, rwc) <= r <= max(0.0, rwc) + wall_t:
                return metal
            return air

        rw = g.rho_horn(az)
        if not math.isfinite(rw):
            return air

        # v46 right output: remove the lower horn part of the corrected aperture, R -> R-W/2.
        in_right_R2_slot = (
            getattr(cfg, "window_model", "vertical_shift") == "horn_cut"
            and g.R2 > 0.0
            and z > 0.0
            and z0_slot <= z <= z1_slot
            and max(0.0, rw - 0.25 * wall_t) <= r <= rw + 1.25 * wall_t
        )
        if in_right_R2_slot:
            return air

        if max(0.0, rw) <= r <= max(0.0, rw) + wall_t:
            return metal
        return air

    return matfun


def component_from_name(name: str):
    require_meep()
    table = {"Er": mp.Er, "Ep": mp.Ep, "Ez": mp.Ez, "Hr": mp.Hr, "Hp": mp.Hp, "Hz": mp.Hz}
    key = name.strip()
    if key not in table:
        raise ValueError(f"Unknown field component {name!r}; use Ez,Er,Ep")
    return table[key]


def sample_active_point(rng: np.random.Generator, g: PHBGeometry, cfg: RunConfig) -> Tuple[float, float]:
    # Rejection sampling in a conservative bounding box.
    rmax_src = max(g.window_r_high(cfg.phb_type), g.R + (g.a if cfg.phb_type == "halfring" else 0.0))
    for _ in range(20000):
        z = rng.uniform(-0.85 * g.L, 0.85 * g.L)
        az = abs(z)
        if az <= g.a:
            if cfg.phb_type == "open":
                rw = g.R
            else:
                rw = phb_central_wall_radius(z, g, cfg.phb_type)
        elif az <= g.L:
            rw = g.rho_horn(az)
        else:
            continue
        if not math.isfinite(rw) or rw <= 0:
            continue
        # avoid the axis for m>0 and avoid exact metal wall
        r_min = 0.08 * max(g.R, 1e-9) if abs(cfg.m) > 0 else 0.0
        r_max = max(r_min + 1e-9, 0.90 * min(rw, rmax_src))
        r = rng.uniform(r_min, r_max)
        return r, z
    return max(0.1 * g.R, 1.0 / max(cfg.resolution, 1)), 0.0


def make_sources(g: PHBGeometry, cfg: RunConfig):
    """Create MEEP sources and a fully reproducible source metadata table.

    v47.1 note
    ----------
    Earlier v47 field maps could show very large |E|^2 peaks, but the archive did
    not contain the source coordinates.  This made it impossible to distinguish
    a genuine PHB/edge concentration from a trivial source-near-field maximum.

    This function now returns ``source_records`` together with the MEEP sources.
    Every record is JSON/CSV-safe and includes the source type, component, center,
    size, phase and amplitude.  Downstream code also compares every field-map
    |E|^2 peak with these source records and writes nearest-source distances into
    metrics.json.
    """
    require_meep()
    rng_seed = int(cfg.seed + 1000 * cfg.m + int(10000 * g.R2))
    rng = np.random.default_rng(rng_seed)
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    fcen = 1.0 / lam
    fwidth = cfg.fwidth_frac * fcen
    comp_names = [s.strip() for s in cfg.source_components.split(",") if s.strip()]
    if not comp_names:
        comp_names = ["Ez"]
    comps = [component_from_name(s) for s in comp_names]

    source_records: List[Dict[str, object]] = []

    def add_record(
        *,
        idx: int,
        source_kind: str,
        component_name_str: str,
        center_r: float,
        center_z: float,
        size_r: float = 0.0,
        size_z: float = 0.0,
        phase_rad: Optional[float] = None,
        amplitude: complex = 1.0 + 0.0j,
    ) -> None:
        amp = complex(amplitude)
        source_records.append({
            "source_index": int(idx),
            "source_kind": str(source_kind),
            "source_mode": str(cfg.source_mode),
            "component": str(component_name_str),
            "center_r": float(center_r),
            "center_phi": 0.0,
            "center_z": float(center_z),
            "size_r": float(size_r),
            "size_phi": 0.0,
            "size_z": float(size_z),
            "r_min": float(center_r - 0.5 * size_r),
            "r_max": float(center_r + 0.5 * size_r),
            "z_min": float(center_z - 0.5 * size_z),
            "z_max": float(center_z + 0.5 * size_z),
            "phase_rad": float(phase_rad) if phase_rad is not None else "",
            "amplitude_real": float(np.real(amp)),
            "amplitude_imag": float(np.imag(amp)),
            "amplitude_abs": float(abs(amp)),
            "frequency": float(fcen),
            "fwidth": float(fwidth),
            "lambda": float(lam),
            "rng_seed": int(rng_seed),
            "m": int(cfg.m),
            "R2_abs": float(g.R2),
            "window_offset": float(g.window_offset),
            "window_offset_mode": str(g.window_offset_mode),
            "M1_r_min": float(g.window_r_low(cfg.phb_type)),
            "M1_r_max": float(g.window_r_high(cfg.phb_type)),
            "M1_r_center": float(g.window_r_center(cfg.phb_type)),
        })

    if cfg.source_mode == "coherent":
        source_rmax = max(g.window_r_high(cfg.phb_type), g.R + (g.a if cfg.phb_type == "halfring" else 0.0))
        if abs(cfg.m) > 0:
            r0 = 0.12 * source_rmax
            r1 = 0.90 * source_rmax
            center = mp.Vector3(0.5 * (r0 + r1), 0, 0)
            size = mp.Vector3(r1 - r0, 0, 1.55 * g.L)
        else:
            center = mp.Vector3(0.5 * source_rmax, 0, 0)
            size = mp.Vector3(source_rmax, 0, 1.55 * g.L)
        add_record(
            idx=0, source_kind="extended_line_volume",
            component_name_str=comp_names[0],
            center_r=float(center.x), center_z=float(center.z),
            size_r=float(size.x), size_z=float(size.z),
            phase_rad=0.0, amplitude=1.0 + 0.0j,
        )
        return [mp.Source(mp.GaussianSource(fcen, fwidth=fwidth), component=comps[0], center=center, size=size)], fcen, fwidth, source_records

    if cfg.source_mode == "single":
        r, z = sample_active_point(rng, g, cfg)
        add_record(
            idx=0, source_kind="point",
            component_name_str=comp_names[0],
            center_r=r, center_z=z,
            size_r=0.0, size_z=0.0,
            phase_rad=0.0, amplitude=1.0 + 0.0j,
        )
        return [mp.Source(mp.GaussianSource(fcen, fwidth=fwidth), component=comps[0], center=mp.Vector3(r, 0, z))], fcen, fwidth, source_records

    if cfg.source_mode != "random":
        raise ValueError("--source-mode must be random, coherent, or single")

    sources = []
    nsrc = max(1, int(cfg.nsrc))
    for i in range(nsrc):
        r, z = sample_active_point(rng, g, cfg)
        phase = float(rng.uniform(0.0, 2.0 * math.pi))
        amp = complex(math.cos(phase), math.sin(phase)) / math.sqrt(nsrc)
        comp_i = comps[i % len(comps)]
        comp_name_i = comp_names[i % len(comp_names)]
        add_record(
            idx=i, source_kind="point_random_phase",
            component_name_str=comp_name_i,
            center_r=r, center_z=z,
            size_r=0.0, size_z=0.0,
            phase_rad=phase, amplitude=amp,
        )
        sources.append(
            mp.Source(mp.GaussianSource(fcen, fwidth=fwidth), component=comp_i,
                      center=mp.Vector3(r, 0, z), amplitude=amp)
        )
    return sources, fcen, fwidth, source_records


def point_to_source_record_distance_rz(r: float, z: float, rec: Dict[str, object]) -> Tuple[float, float, float, bool]:
    """Distance from a point (r,z) to a point or extended rectangular source.

    Returns (distance, dr, dz, inside_extent).  For point sources size_r=size_z=0.
    For coherent extended sources the distance is to the rectangular source
    support in the (r,z) diagnostic section.
    """
    cr = safe_float(rec.get("center_r"))
    cz = safe_float(rec.get("center_z"))
    sr = max(0.0, safe_float(rec.get("size_r"), 0.0))
    sz = max(0.0, safe_float(rec.get("size_z"), 0.0))
    r0 = cr - 0.5 * sr
    r1 = cr + 0.5 * sr
    z0 = cz - 0.5 * sz
    z1 = cz + 0.5 * sz
    if sr > 0.0 or sz > 0.0:
        dr_out = 0.0 if r0 <= r <= r1 else min(abs(r - r0), abs(r - r1))
        dz_out = 0.0 if z0 <= z <= z1 else min(abs(z - z0), abs(z - z1))
        inside = (r0 <= r <= r1) and (z0 <= z <= z1)
        signed_dr = 0.0 if r0 <= r <= r1 else (r - cr)
        signed_dz = 0.0 if z0 <= z <= z1 else (z - cz)
        return float(math.hypot(dr_out, dz_out)), float(signed_dr), float(signed_dz), bool(inside)
    dr = float(r - cr)
    dz = float(z - cz)
    return float(math.hypot(dr, dz)), dr, dz, False


def add_source_proximity_metrics(row: Dict[str, object], source_records: List[Dict[str, object]], cfg: RunConfig) -> None:
    """Annotate field-map |E|² peaks with nearest-source distances.

    The new metrics make it possible to classify an enormous |E|² peak as:
      * likely source-near-field artefact, if it lies inside/near a source;
      * likely PHB/edge concentration, if it is far from all sources.

    All thresholds are deliberately conservative and stored in the output.
    """
    row["source_count"] = int(len(source_records))
    row["source_coordinates_available"] = bool(source_records)
    if not source_records:
        return
    cell = 1.0 / max(float(cfg.resolution), 1.0)
    near_threshold = 3.0 * cell
    row["source_near_peak_threshold_cells"] = 3.0
    row["source_near_peak_threshold_length"] = near_threshold

    # Also write global compact source coordinate lists into metrics.json.
    row["source_centers_r_json"] = json.dumps([safe_float(x.get("center_r")) for x in source_records])
    row["source_centers_z_json"] = json.dumps([safe_float(x.get("center_z")) for x in source_records])
    row["source_components_json"] = json.dumps([str(x.get("component")) for x in source_records])

    peak_prefixes: List[str] = []
    for key in row.keys():
        if key.endswith("_E2_peak_r"):
            prefix = key[:-len("_E2_peak_r")]
            if f"{prefix}_E2_peak_z" in row:
                peak_prefixes.append(prefix)

    for prefix in sorted(set(peak_prefixes)):
        pr = safe_float(row.get(f"{prefix}_E2_peak_r"))
        pz = safe_float(row.get(f"{prefix}_E2_peak_z"))
        if not (math.isfinite(pr) and math.isfinite(pz)):
            continue
        best = None
        for rec in source_records:
            dist, dr, dz, inside = point_to_source_record_distance_rz(pr, pz, rec)
            item = (dist, dr, dz, inside, rec)
            if best is None or item[0] < best[0]:
                best = item
        if best is None:
            continue
        dist, dr, dz, inside, rec = best
        idx = int(safe_float(rec.get("source_index"), -1))
        row[f"{prefix}_E2_peak_nearest_source_index"] = idx
        row[f"{prefix}_E2_peak_nearest_source_component"] = rec.get("component")
        row[f"{prefix}_E2_peak_nearest_source_kind"] = rec.get("source_kind")
        row[f"{prefix}_E2_peak_nearest_source_r"] = safe_float(rec.get("center_r"))
        row[f"{prefix}_E2_peak_nearest_source_z"] = safe_float(rec.get("center_z"))
        row[f"{prefix}_E2_peak_nearest_source_distance_rz"] = float(dist)
        row[f"{prefix}_E2_peak_nearest_source_distance_cells"] = float(dist / max(cell, 1e-300))
        row[f"{prefix}_E2_peak_nearest_source_dr"] = float(dr)
        row[f"{prefix}_E2_peak_nearest_source_dz"] = float(dz)
        row[f"{prefix}_E2_peak_inside_source_extent"] = bool(inside)
        row[f"{prefix}_E2_peak_likely_source_artifact"] = bool(inside or dist <= near_threshold)

def estimate_cell(g: PHBGeometry, cfg: RunConfig) -> Tuple[float, float, float, Optional[float]]:
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    dpml = cfg.dpml if cfg.dpml is not None else cfg.dpml_over_lambda * lam
    z0, z1_eff, _, _ = effective_slot_interval(g, cfg)
    z_ap = z1_eff + cfg.aperture_offset_cells / max(cfg.resolution, 1)
    far_z = None
    if cfg.enable_far_monitor:
        far_z = z_ap + cfg.far_distance
    r_needed = max(g.window_r_high(cfg.phb_type), g.R + (g.a if cfg.phb_type == "halfring" else 0.0)) + effective_wall_t(cfg, lam)
    if far_z is not None:
        width_m1 = max(g.window_width(cfg.phb_type), 1.0 / max(cfg.resolution, 1))
        capture_half = max(0.5 * width_m1, cfg.far_distance * math.tan(math.radians(cfg.far_capture_angle_deg)))
        # v46: M2 is centered on the same shifted radial line as M1.
        ref_r = g.window_r_center(cfg.phb_type)
        r_needed = max(r_needed, ref_r + capture_half)
    rmax = r_needed + dpml + max(0.35 * g.R, 0.8 * lam, 0.3)
    z_pos_needed = max(g.L, z_ap + 0.25 * g.R)
    if far_z is not None:
        z_pos_needed = max(z_pos_needed, far_z + 0.25 * g.R)
    z_half = z_pos_needed + dpml + max(0.35 * g.R, 0.8 * lam, 0.3)
    zspan = 2.0 * z_half
    return rmax, zspan, z_ap, far_z


def make_simulation(g: PHBGeometry, cfg: RunConfig):
    require_meep()
    rmax, zspan, z_ap, far_z = estimate_cell(g, cfg)
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    cfg.dpml = cfg.dpml if cfg.dpml is not None else cfg.dpml_over_lambda * lam
    sources, fcen, _, source_records = make_sources(g, cfg)
    matfun = make_material_function(g, cfg)
    courant = min(0.5, 1.0 / (abs(cfg.m) + 0.8)) if abs(cfg.m) > 0 else 0.5
    sim = mp.Simulation(
        cell_size=mp.Vector3(rmax, 0, zspan),
        boundary_layers=[mp.PML(cfg.dpml, direction=mp.R), mp.PML(cfg.dpml, direction=mp.Z)],
        resolution=cfg.resolution,
        dimensions=mp.CYLINDRICAL,
        m=int(cfg.m),
        sources=sources,
        material_function=matfun,
        force_complex_fields=True,
        accurate_fields_near_cylorigin=True,
        Courant=courant,
    )
    return sim, rmax, zspan, z_ap, far_z, fcen, source_records


def far_monitor_interval(g: PHBGeometry, cfg: RunConfig, rmax: float, zspan: float) -> Tuple[float, float, float]:
    """Radial interval for M2, the optional far-zone monitor.

    M2 is a z-plane monitor parallel to M1. In v46 it follows
    the same shifted M1 window center R+shift.  The width is large enough
    to cover the angular capture cone selected by --far-capture-angle-deg.
    """
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    dpml = cfg.dpml if cfg.dpml is not None else cfg.dpml_over_lambda * lam
    margin = max(0.12, 0.25 * lam)
    width_m1 = max(g.window_width(cfg.phb_type), 1.0 / max(cfg.resolution, 1))
    capture_half = max(0.5 * width_m1, cfg.far_distance * math.tan(math.radians(cfg.far_capture_angle_deg)))
    # v46: M2 follows the same shifted M1 window center, not rho=R.
    center = g.window_r_center(cfg.phb_type)
    r0 = max(0.0, center - capture_half)
    r1 = min(rmax - dpml - margin, center + capture_half)
    return r0, r1, center


def add_far_geometric_angle_metrics(profile_rows: List[Dict[str, object]], reference_r: float, distance: float, label: str) -> Dict[str, object]:
    """Add geometric divergence metrics at M2 from radial displacement.

    These metrics are different from local Poynting angle.  They estimate
    how much forward energy lies within an angular tube around the selected reference radius.
    """
    if not profile_rows or distance <= 0:
        return {}
    r = np.array([safe_float(x.get("r")) for x in profile_rows], dtype=float)
    w = np.array([max(0.0, safe_float(x.get("forward_weight"))) for x in profile_rows], dtype=float)
    theta = np.degrees(np.arctan2(np.abs(r - reference_r), max(float(distance), 1e-300)))
    fsum = float(np.sum(w))
    if fsum <= 0:
        return {f"{label}_geom_forward_flux_proxy": fsum}
    for i, row in enumerate(profile_rows):
        row[f"{label}_geom_theta_from_focal_axis_deg"] = float(theta[i])
    return {
        f"{label}_geom_reference_r": float(reference_r),
        f"{label}_geom_distance_from_M1": float(distance),
        f"{label}_geom_forward_flux_proxy": fsum,
        f"{label}_geom_theta95_deg": weighted_percentile(theta, w, 0.95),
        f"{label}_geom_theta90_deg": weighted_percentile(theta, w, 0.90),
        f"{label}_geom_theta50_deg": weighted_percentile(theta, w, 0.50),
        f"{label}_geom_lobe_fraction_5deg": float(np.sum(w[theta <= 5.0]) / max(fsum, 1e-300)),
        f"{label}_geom_lobe_fraction_10deg": float(np.sum(w[theta <= 10.0]) / max(fsum, 1e-300)),
        f"{label}_geom_lobe_fraction_20deg": float(np.sum(w[theta <= 20.0]) / max(fsum, 1e-300)),
        f"{label}_geom_density_proxy_5deg": float(np.sum(w[theta <= 5.0]) / max(math.radians(5.0), 1e-300)),
    }


# -----------------------------
# Full 2D field-map diagnostics, v47/v46.1
# -----------------------------

def parse_field_map_regions(s: str) -> List[str]:
    aliases = {
        "none": "none",
        "off": "none",
        "all": "all",
        "global": "interior",
        "whole": "interior",
        "phb": "interior",
        "slot": "near_slot",
        "slit": "near_slot",
        "aperture": "near_slot",
        "pre": "pre_slot",
        "before_slot": "pre_slot",
        "post": "post_slot",
        "after_slot": "post_slot",
        "exit": "post_slot",
        "corridor": "focal_corridor",
        "focal": "focal_corridor",
    }
    raw: List[str] = []
    for part in str(s).replace(";", ",").split(","):
        v = part.strip().lower().replace("-", "_")
        if not v:
            continue
        raw.append(aliases.get(v, v))
    if not raw or "none" in raw:
        return []
    if "all" in raw:
        return ["interior", "near_slot", "pre_slot", "post_slot", "focal_corridor"]
    allowed = ["interior", "near_slot", "pre_slot", "post_slot", "focal_corridor"]
    out: List[str] = []
    for v in raw:
        if v not in allowed:
            raise ValueError(f"Unknown field-map region {v!r}; use none,all,interior,near_slot,pre_slot,post_slot,focal_corridor")
        if v not in out:
            out.append(v)
    return out


def field_map_components_from_cfg(cfg: RunConfig):
    require_meep()
    table = {"Er": mp.Er, "Ep": mp.Ep, "Ez": mp.Ez, "Hr": mp.Hr, "Hp": mp.Hp, "Hz": mp.Hz}
    out = []
    for part in str(cfg.field_map_components).replace(";", ",").split(","):
        key = part.strip()
        if not key:
            continue
        if key not in table:
            raise ValueError(f"Unknown field-map component {key!r}; use Er,Ep,Ez,Hr,Hp,Hz")
        if table[key] not in out:
            out.append(table[key])
    # Poynting maps require all six components.  Force all components unless
    # the user intentionally provided a subset and accepts missing derived maps.
    return out or [mp.Er, mp.Ep, mp.Ez, mp.Hr, mp.Hp, mp.Hz]


def component_name(component) -> str:
    require_meep()
    table = {mp.Er: "Er", mp.Ep: "Ep", mp.Ez: "Ez", mp.Hr: "Hr", mp.Hp: "Hp", mp.Hz: "Hz"}
    return table.get(component, str(component))


def clamp_interval(lo: float, hi: float, bound_lo: float, bound_hi: float) -> Tuple[float, float]:
    lo2 = max(float(bound_lo), float(lo))
    hi2 = min(float(bound_hi), float(hi))
    if hi2 <= lo2:
        mid = 0.5 * (float(bound_lo) + float(bound_hi))
        lo2, hi2 = mid - 0.5 / 1000.0, mid + 0.5 / 1000.0
    return lo2, hi2


def estimate_region_points(r_min: float, r_max: float, z_min: float, z_max: float, resolution: int) -> int:
    nr = max(1, int(math.ceil(max(0.0, r_max - r_min) * resolution)) + 1)
    nz = max(1, int(math.ceil(max(0.0, z_max - z_min) * resolution)) + 1)
    return int(nr * nz)


def field_map_region_specs(g: PHBGeometry, cfg: RunConfig, rmax: float, zspan: float, z_ap: float, far_z: Optional[float]) -> List[Dict[str, object]]:
    """Build requested 2D DFT-map regions in the cylindrical (r,z) section.

    These field maps are diagnostic volumes/slices. They are not used to
    classify the aperture candidate; they are saved so that later scripts can
    recompute field markers, Poynting topology, phase maps, edge effects, and
    near-slot artefacts without rerunning MEEP.
    """
    regions = parse_field_map_regions(cfg.field_map_regions)
    if not cfg.save_field_maps or not regions:
        return []

    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    dpml = cfg.dpml if cfg.dpml is not None else cfg.dpml_over_lambda * lam
    margin = max(0.08, 0.25 * lam)
    safe_r0 = 0.0
    safe_r1 = max(safe_r0 + 1e-6, rmax - dpml - margin)
    safe_z0 = -0.5 * zspan + dpml + margin
    safe_z1 = +0.5 * zspan - dpml - margin

    r_low = g.window_r_low(cfg.phb_type)
    r_high = g.window_r_high(cfg.phb_type)
    rc = g.window_r_center(cfg.phb_type)
    pad_r = max(0.0, cfg.field_map_r_pad)
    z_half = max(1.0 / max(cfg.resolution, 1), cfg.field_map_z_half)
    corr_hw = max(0.5 * g.R2, cfg.field_map_corridor_half_width)

    default_interior_rmax = max(
        g.window_r_high(cfg.phb_type) + pad_r,
        g.R + (g.a if cfg.phb_type == "halfring" else 0.5 * g.R2) + effective_wall_t(cfg, lam) + pad_r,
        g.R + pad_r,
    )
    interior_rmax = cfg.field_map_interior_rmax if cfg.field_map_interior_rmax is not None else default_interior_rmax

    spec_defs: Dict[str, Tuple[float, float, float, float, str]] = {
        "interior": (0.0, interior_rmax, -g.L, min(z_ap + z_half, g.L + z_half),
                     "whole PHB/interior diagnostic section, PML excluded"),
        "near_slot": (r_low - pad_r, r_high + pad_r, g.a - z_half, g.a + z_half,
                      "near-slot box around the right shifted R2 aperture"),
        "pre_slot": (r_low - pad_r, r_high + pad_r, g.a - z_half, g.a,
                     "inside-PHB upstream side of the slot"),
        "post_slot": (r_low - pad_r, r_high + pad_r, g.a, g.a + z_half,
                      "downstream/exit side immediately after the slot"),
        "focal_corridor": (rc - corr_hw, rc + corr_hw, -g.L, min(z_ap + z_half, g.L + z_half),
                           "longitudinal corridor around the shifted output line R+shift"),
    }

    specs: List[Dict[str, object]] = []
    for name in regions:
        r0, r1, z0, z1, note = spec_defs[name]
        r0, r1 = clamp_interval(r0, r1, safe_r0, safe_r1)
        z0, z1 = clamp_interval(z0, z1, safe_z0, safe_z1)
        npoints = estimate_region_points(r0, r1, z0, z1, cfg.resolution)
        spec: Dict[str, object] = {
            "name": name,
            "r_min": float(r0), "r_max": float(r1),
            "z_min": float(z0), "z_max": float(z1),
            "r_center": float(0.5 * (r0 + r1)),
            "z_center": float(0.5 * (z0 + z1)),
            "r_size": float(r1 - r0),
            "z_size": float(z1 - z0),
            "estimated_points": int(npoints),
            "note": note,
            "skipped": False,
        }
        if npoints > int(cfg.field_map_max_points_per_region):
            spec["skipped"] = True
            spec["skip_reason"] = (
                f"estimated grid points {npoints} > --field-map-max-points-per-region "
                f"{cfg.field_map_max_points_per_region}; increase limit or reduce region/padding"
            )
        specs.append(spec)
    return specs


def add_field_map_diagnostics(sim, g: PHBGeometry, cfg: RunConfig, rmax: float, zspan: float, z_ap: float, far_z: Optional[float], fcen: float):
    require_meep()
    field_maps: Dict[str, Dict[str, object]] = {}
    plan = field_map_region_specs(g, cfg, rmax, zspan, z_ap, far_z)
    if not cfg.save_field_maps:
        return field_maps, plan
    comps = field_map_components_from_cfg(cfg)
    for spec in plan:
        if spec.get("skipped"):
            continue
        center = mp.Vector3(float(spec["r_center"]), 0, float(spec["z_center"]))
        size = mp.Vector3(float(spec["r_size"]), 0, float(spec["z_size"]))
        dft = sim.add_dft_fields(comps, fcen, 0, 1, center=center, size=size)
        field_maps[str(spec["name"])] = {"dft": dft, "components": comps, "spec": spec}
    return field_maps, plan


def get_dft_map_array(sim, dft, component, spec: Dict[str, object], resolution: int) -> np.ndarray:
    arr = np.asarray(sim.get_dft_array(dft, component, 0)).squeeze().astype(complex)
    if arr.ndim == 0:
        arr = arr.reshape((1, 1))
    if arr.ndim == 1:
        nr_est = max(1, int(round(float(spec["r_size"]) * resolution)) + 1)
        nz_est = max(1, int(round(float(spec["z_size"]) * resolution)) + 1)
        if arr.size == nr_est * nz_est:
            arr = arr.reshape((nr_est, nz_est))
        else:
            arr = arr.reshape((arr.size, 1))
    if arr.ndim > 2:
        arr = np.squeeze(arr)
        if arr.ndim > 2:
            arr = arr.reshape((arr.shape[0], int(np.prod(arr.shape[1:]))))
    # Try to orient as [r_index, z_index].
    nr_est = max(1, int(round(float(spec["r_size"]) * resolution)) + 1)
    nz_est = max(1, int(round(float(spec["z_size"]) * resolution)) + 1)
    if arr.ndim == 2 and abs(arr.shape[0] - nz_est) < abs(arr.shape[0] - nr_est) and abs(arr.shape[1] - nr_est) < abs(arr.shape[1] - nz_est):
        arr = arr.T
    return arr


def compute_field_map_arrays(arrays: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    # Missing components are filled with zero to preserve derived-map creation,
    # but raw component availability is stored in metadata.
    shape = next(iter(arrays.values())).shape
    z0 = np.zeros(shape, dtype=complex)
    Er = arrays.get("Er", z0)
    Ep = arrays.get("Ep", z0)
    Ez = arrays.get("Ez", z0)
    Hr = arrays.get("Hr", z0)
    Hp = arrays.get("Hp", z0)
    Hz = arrays.get("Hz", z0)
    E2 = 0.5 * (np.abs(Er) ** 2 + np.abs(Ep) ** 2 + np.abs(Ez) ** 2)
    H2 = 0.5 * (np.abs(Hr) ** 2 + np.abs(Hp) ** 2 + np.abs(Hz) ** 2)
    S_r = 0.5 * np.real(Ep * np.conj(Hz) - Ez * np.conj(Hp))
    S_phi = 0.5 * np.real(Ez * np.conj(Hr) - Er * np.conj(Hz))
    S_z = 0.5 * np.real(Er * np.conj(Hp) - Ep * np.conj(Hr))
    theta_poloidal = np.degrees(np.arctan2(np.abs(S_r), np.maximum(S_z, 1e-300)))
    theta_poloidal = np.where(S_z > 0, theta_poloidal, 180.0)
    theta_3d = np.degrees(np.arctan2(np.sqrt(S_r * S_r + S_phi * S_phi), np.maximum(S_z, 1e-300)))
    theta_3d = np.where(S_z > 0, theta_3d, 180.0)
    S_mag = np.sqrt(S_r * S_r + S_phi * S_phi + S_z * S_z)
    phase_Ez = np.angle(Ez)
    return {
        "E2": E2, "H2": H2, "S_r": S_r, "S_phi": S_phi, "S_z": S_z,
        "S_mag": S_mag, "theta_poloidal_from_z_deg": theta_poloidal,
        "theta_3d_from_z_deg": theta_3d, "phase_Ez_rad": phase_Ez,
    }


def approximate_phb_air_mask(RR: np.ndarray, ZZ: np.ndarray, g: PHBGeometry, cfg: RunConfig) -> np.ndarray:
    """Approximate non-metal PHB interior mask for field-map metrics.

    The raw maps are saved unmasked.  This mask is a diagnostic helper only;
    it avoids interpreting the PML/outside box as PHB interior when computing
    summary markers.
    """
    mask = np.zeros_like(RR, dtype=bool)
    for j in range(ZZ.shape[1]):
        zcol = ZZ[:, j]
        az = np.abs(zcol)
        wall = np.zeros_like(zcol, dtype=float)
        central = az < g.a
        horn = (az >= g.a) & (az <= g.L)
        if np.any(central):
            if cfg.phb_type == "open":
                wall[central] = g.R
            else:
                # vectorized halfring/cylinder central closure approximation
                if cfg.phb_type == "cylinder":
                    wall[central] = g.R + 0.5 * g.R2
                else:
                    wall[central] = g.R + np.sqrt(np.maximum(g.a * g.a - zcol[central] * zcol[central], 0.0))
        if np.any(horn):
            wall[horn] = g.R - g.b * np.sqrt(np.maximum((az[horn] / g.a) ** 2 - 1.0, 0.0))
        mask[:, j] = (wall > 0.0) & (RR[:, j] <= np.maximum(wall, 0.0))
    return mask



def diagnostic_material_code_point(r: float, z: float, g: PHBGeometry, cfg: RunConfig) -> int:
    """Classify a point in the diagnostic (r,z) section.

    Codes saved in material_mask arrays:
      0 = external air / outside diagnostic PHB body;
      1 = PHB air/interior side of the cavity;
      2 = metallic wall / closure layer proxy;
      3 = explicitly opened R2 aperture/window band.

    This reproduces the material_function logic at diagnostic resolution and is
    intended for auditing whether M1 crosses a real aperture or only a selected
    free-field stripe.  It is not used by MEEP during time stepping.
    """
    lam = wavelength_from_a(g.a, cfg.a_over_lambda)
    wall_t = effective_wall_t(cfg, lam)
    z0_slot, z1_slot, _, _ = effective_slot_interval(g, cfg)
    r = float(r); z = float(z); az = abs(z)
    if r < 0.0 or az > g.L:
        return 0

    if getattr(cfg, "window_model", "vertical_shift") == "vertical_shift" and g.R2 > 0.0 and z > 0.0:
        r0_win = g.window_r_low(cfg.phb_type)
        r1_win = g.window_r_high(cfg.phb_type)
        if z0_slot <= z <= z1_slot and (r0_win - 0.50 * wall_t) <= r <= (r1_win + 0.50 * wall_t):
            return 3

    if cfg.phb_type == "open" and z < 0.0:
        r0_win = max(0.0, g.window_r_low(cfg.phb_type))
        r1_win = g.window_r_high(cfg.phb_type)
        if abs(z + g.a) <= 0.65 * wall_t and (r0_win - wall_t) <= r <= (r1_win + wall_t):
            return 2

    if az < g.a:
        if cfg.phb_type == "open":
            return 1 if r <= g.R else 0
        rwc = phb_central_wall_radius(z, g, cfg.phb_type)
        if getattr(cfg, "window_model", "vertical_shift") == "horn_cut" and is_right_upper_window_on_halfring(z, g, cfg, wall_t):
            return 3
        if max(0.0, rwc) <= r <= max(0.0, rwc) + wall_t:
            return 2
        return 1 if r < max(0.0, rwc) else 0

    rw = g.rho_horn(az)
    if not math.isfinite(rw) or rw <= 0.0:
        return 0
    if (
        getattr(cfg, "window_model", "vertical_shift") == "horn_cut"
        and g.R2 > 0.0 and z > 0.0 and z0_slot <= z <= z1_slot
        and max(0.0, rw - 0.25 * wall_t) <= r <= rw + 1.25 * wall_t
    ):
        return 3
    if max(0.0, rw) <= r <= max(0.0, rw) + wall_t:
        return 2
    return 1 if r < max(0.0, rw) else 0


def diagnostic_material_mask_grid(r_grid: np.ndarray, z_grid: np.ndarray, g: PHBGeometry, cfg: RunConfig) -> np.ndarray:
    mask = np.zeros((len(r_grid), len(z_grid)), dtype=np.int16)
    for ii, rr in enumerate(r_grid):
        for jj, zz in enumerate(z_grid):
            mask[ii, jj] = diagnostic_material_code_point(float(rr), float(zz), g, cfg)
    return mask


def source_exclusion_mask_grid(r_grid: np.ndarray, z_grid: np.ndarray, source_records: List[Dict[str, object]], cfg: RunConfig) -> np.ndarray:
    Rg, Zg = np.meshgrid(r_grid, z_grid, indexing="ij")
    mask = np.zeros_like(Rg, dtype=bool)
    threshold = float(getattr(cfg, "source_exclusion_cells", 3.0)) / max(float(cfg.resolution), 1.0)
    for rec in source_records or []:
        cr = safe_float(rec.get("center_r")); cz = safe_float(rec.get("center_z"))
        if not (math.isfinite(cr) and math.isfinite(cz)):
            continue
        sr = max(0.0, safe_float(rec.get("size_r"), 0.0))
        sz = max(0.0, safe_float(rec.get("size_z"), 0.0))
        if sr > 0.0 or sz > 0.0:
            r0 = cr - 0.5 * sr - threshold
            r1 = cr + 0.5 * sr + threshold
            z0 = cz - 0.5 * sz - threshold
            z1 = cz + 0.5 * sz + threshold
            mask |= (Rg >= r0) & (Rg <= r1) & (Zg >= z0) & (Zg <= z1)
        else:
            mask |= ((Rg - cr) ** 2 + (Zg - cz) ** 2) <= threshold ** 2
    return mask


def source_excluded_E2_metrics(name: str, E2: np.ndarray, source_mask: np.ndarray, r_grid: np.ndarray, z_grid: np.ndarray, cfg: RunConfig) -> Dict[str, object]:
    finite = np.isfinite(E2)
    keep = finite & (~source_mask)
    out: Dict[str, object] = {
        f"fieldmap_{name}_source_exclusion_cells": float(getattr(cfg, "source_exclusion_cells", 3.0)),
        f"fieldmap_{name}_source_exclusion_fraction": float(np.mean(source_mask)) if source_mask.size else float("nan"),
        f"fieldmap_{name}_source_excluded_finite_points": int(np.sum(keep)),
    }
    if not np.any(keep):
        out[f"fieldmap_{name}_E2_peak_excl_sources"] = float("nan")
        return out
    masked_E2 = np.where(keep, E2, np.nan)
    peak_idx = np.unravel_index(int(np.nanargmax(masked_E2)), masked_E2.shape)
    mean_ex = float(np.nanmean(masked_E2))
    out.update({
        f"fieldmap_{name}_E2_peak_excl_sources": float(np.nanmax(masked_E2)),
        f"fieldmap_{name}_E2_mean_excl_sources": mean_ex,
        f"fieldmap_{name}_E2_peak_over_mean_excl_sources": float(np.nanmax(masked_E2) / max(mean_ex, 1e-300)),
        f"fieldmap_{name}_E2_peak_excl_sources_r": float(r_grid[peak_idx[0]]),
        f"fieldmap_{name}_E2_peak_excl_sources_z": float(z_grid[peak_idx[1]]),
    })
    return out


def plot_material_mask_region(cdir: Path, name: str, spec: Dict[str, object], r_grid: np.ndarray, z_grid: np.ndarray, material_mask: np.ndarray) -> None:
    if plt is None:
        return
    ensure_dir(cdir / "field_maps")
    extent = [float(spec["z_min"]), float(spec["z_max"]), float(spec["r_min"]), float(spec["r_max"])]
    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=150)
    im = ax.imshow(material_mask, origin="lower", aspect="auto", extent=extent, vmin=0, vmax=3)
    ax.set_xlabel("z")
    ax.set_ylabel("r")
    ax.set_title(f"{name}: diagnostic material mask (0 outside, 1 air, 2 metal, 3 R2 window)")
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["outside", "PHB air", "metal", "R2 window"])
    fig.tight_layout()
    fig.savefig(cdir / "field_maps" / f"fieldmap_{name}_material_mask.png")
    plt.close(fig)


def plot_source_excluded_E2_region(cdir: Path, name: str, spec: Dict[str, object], r_grid: np.ndarray, z_grid: np.ndarray, E2_source_excluded: np.ndarray) -> None:
    if plt is None:
        return
    ensure_dir(cdir / "field_maps")
    extent = [float(spec["z_min"]), float(spec["z_max"]), float(spec["r_min"]), float(spec["r_max"])]
    data = np.asarray(E2_source_excluded, dtype=float)
    finite = np.isfinite(data) & (data > 0)
    if not np.any(finite):
        return
    floor = np.nanmax(data[finite]) * 1e-8
    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=150)
    im = ax.imshow(np.log10(np.maximum(data, floor)), origin="lower", aspect="auto", extent=extent)
    ax.set_xlabel("z")
    ax.set_ylabel("r")
    ax.set_title(f"{name}: log10(|E|²) with source-near-field excluded")
    fig.colorbar(im, ax=ax, label="log10(E2), source excluded")
    fig.tight_layout()
    fig.savefig(cdir / "field_maps" / f"fieldmap_{name}_E2_log_source_excluded.png")
    plt.close(fig)


def field_map_metrics(name: str, spec: Dict[str, object], derived: Dict[str, np.ndarray], r_grid: np.ndarray, z_grid: np.ndarray, g: PHBGeometry, cfg: RunConfig) -> Dict[str, object]:
    Rg, Zg = np.meshgrid(r_grid, z_grid, indexing="ij")
    area_w = 2.0 * math.pi * np.maximum(Rg, 0.0)
    if len(r_grid) > 1:
        dr = abs(float(r_grid[1] - r_grid[0]))
    else:
        dr = max(float(spec["r_size"]), 1.0 / max(cfg.resolution, 1))
    if len(z_grid) > 1:
        dz = abs(float(z_grid[1] - z_grid[0]))
    else:
        dz = max(float(spec["z_size"]), 1.0 / max(cfg.resolution, 1))
    weights = area_w * dr * dz
    E2 = derived["E2"]
    S_r = derived["S_r"]
    S_phi = derived["S_phi"]
    S_z = derived["S_z"]
    S_mag = derived["S_mag"]
    theta_p = derived["theta_poloidal_from_z_deg"]
    theta_3d = derived["theta_3d_from_z_deg"]
    phase = derived["phase_Ez_rad"]

    mask_inside = approximate_phb_air_mask(Rg, Zg, g, cfg)
    finite = np.isfinite(E2) & np.isfinite(S_z)
    active = finite & (E2 > 0)
    forward = np.maximum(S_z, 0.0) * weights
    backward = np.maximum(-S_z, 0.0) * weights
    fsum = float(np.nansum(forward[finite]))
    bsum = float(np.nansum(backward[finite]))
    abs_z = float(np.nansum(np.abs(S_z[finite]) * weights[finite]))
    eweights = np.maximum(E2, 0.0) * weights
    esum = float(np.nansum(eweights[finite]))
    ssum = float(np.nansum(np.maximum(S_mag, 0.0)[finite] * weights[finite]))

    peak_idx = np.unravel_index(int(np.nanargmax(E2)), E2.shape)
    theta_min_idx = np.unravel_index(int(np.nanargmin(theta_p)), theta_p.shape)

    # Poynting topology proxy: curl_phi = dS_r/dz - dS_z/dr in the meridional plane.
    try:
        dSrdz = np.gradient(S_r, dz, axis=1)
        dSzdr = np.gradient(S_z, dr, axis=0)
        curl_phi = dSrdz - dSzdr
        curl_abs = np.abs(curl_phi)
        curl_mean = float(np.nanmean(curl_abs[finite]))
        curl_peak = float(np.nanmax(curl_abs[finite]))
        curl_proxy = float(np.nansum(curl_abs[finite] * weights[finite]) / max(ssum, 1e-300))
    except Exception:
        curl_mean = curl_peak = curl_proxy = float("nan")

    sr_sign = np.sign(S_r)
    sz_sign = np.sign(S_z)
    sr_zero_cross = int(np.nansum(np.abs(np.diff(sr_sign, axis=0)) > 1.0) + np.nansum(np.abs(np.diff(sr_sign, axis=1)) > 1.0))
    sz_zero_cross = int(np.nansum(np.abs(np.diff(sz_sign, axis=0)) > 1.0) + np.nansum(np.abs(np.diff(sz_sign, axis=1)) > 1.0))

    out: Dict[str, object] = {
        f"fieldmap_{name}_r_min": float(spec["r_min"]),
        f"fieldmap_{name}_r_max": float(spec["r_max"]),
        f"fieldmap_{name}_z_min": float(spec["z_min"]),
        f"fieldmap_{name}_z_max": float(spec["z_max"]),
        f"fieldmap_{name}_nr": int(E2.shape[0]),
        f"fieldmap_{name}_nz": int(E2.shape[1]),
        f"fieldmap_{name}_points": int(E2.size),
        f"fieldmap_{name}_inside_mask_fraction": float(np.mean(mask_inside)) if mask_inside.size else float("nan"),
        f"fieldmap_{name}_E2_integral_proxy": esum,
        f"fieldmap_{name}_E2_peak": float(np.nanmax(E2)),
        f"fieldmap_{name}_E2_mean": float(np.nanmean(E2[finite])) if np.any(finite) else float("nan"),
        f"fieldmap_{name}_E2_peak_over_mean": float(np.nanmax(E2) / max(np.nanmean(E2[finite]), 1e-300)) if np.any(finite) else float("nan"),
        f"fieldmap_{name}_E2_peak_r": float(r_grid[peak_idx[0]]),
        f"fieldmap_{name}_E2_peak_z": float(z_grid[peak_idx[1]]),
        f"fieldmap_{name}_S_z_forward_integral_proxy": fsum,
        f"fieldmap_{name}_S_z_backward_integral_proxy": bsum,
        f"fieldmap_{name}_S_z_backward_fraction_abs": bsum / max(abs_z, 1e-300),
        f"fieldmap_{name}_S_phi_abs_fraction_of_S": float(np.nansum(np.abs(S_phi[finite]) * weights[finite]) / max(ssum, 1e-300)) if ssum > 0 else float("nan"),
        f"fieldmap_{name}_theta_poloidal50_deg": weighted_percentile(theta_p.ravel(), forward.ravel(), 0.50),
        f"fieldmap_{name}_theta_poloidal90_deg": weighted_percentile(theta_p.ravel(), forward.ravel(), 0.90),
        f"fieldmap_{name}_theta_poloidal95_deg": weighted_percentile(theta_p.ravel(), forward.ravel(), 0.95),
        f"fieldmap_{name}_theta_3d50_deg": weighted_percentile(theta_3d.ravel(), forward.ravel(), 0.50),
        f"fieldmap_{name}_theta_3d90_deg": weighted_percentile(theta_3d.ravel(), forward.ravel(), 0.90),
        f"fieldmap_{name}_theta_3d95_deg": weighted_percentile(theta_3d.ravel(), forward.ravel(), 0.95),
        f"fieldmap_{name}_theta_poloidal_min_deg": float(np.nanmin(theta_p)),
        f"fieldmap_{name}_theta_poloidal_min_r": float(r_grid[theta_min_idx[0]]),
        f"fieldmap_{name}_theta_poloidal_min_z": float(z_grid[theta_min_idx[1]]),
        f"fieldmap_{name}_phase_Ez_rms_rad": weighted_circular_phase_rms(phase.ravel(), eweights.ravel()),
        f"fieldmap_{name}_poynting_curl_phi_abs_mean": curl_mean,
        f"fieldmap_{name}_poynting_curl_phi_abs_peak": curl_peak,
        f"fieldmap_{name}_poynting_curl_phi_proxy": curl_proxy,
        f"fieldmap_{name}_S_r_zero_crossing_count": sr_zero_cross,
        f"fieldmap_{name}_S_z_zero_crossing_count": sz_zero_cross,
    }
    if esum > 0:
        out[f"fieldmap_{name}_E2_centroid_r"] = float(np.nansum(Rg * eweights) / max(esum, 1e-300))
        out[f"fieldmap_{name}_E2_centroid_z"] = float(np.nansum(Zg * eweights) / max(esum, 1e-300))
    if fsum > 0:
        out[f"fieldmap_{name}_S_z_forward_centroid_r"] = float(np.nansum(Rg * forward) / max(fsum, 1e-300))
        out[f"fieldmap_{name}_S_z_forward_centroid_z"] = float(np.nansum(Zg * forward) / max(fsum, 1e-300))
        out[f"fieldmap_{name}_lobe_fraction_5deg"] = float(np.nansum(forward[theta_p <= 5.0]) / max(fsum, 1e-300))
        out[f"fieldmap_{name}_lobe_fraction_10deg"] = float(np.nansum(forward[theta_p <= 10.0]) / max(fsum, 1e-300))
    return out


def plot_field_map_region(cdir: Path, name: str, spec: Dict[str, object], r_grid: np.ndarray, z_grid: np.ndarray, derived: Dict[str, np.ndarray]) -> None:
    if plt is None:
        return
    ensure_dir(cdir / "field_maps")
    extent = [float(z_grid[0]), float(z_grid[-1]), float(r_grid[0]), float(r_grid[-1])]
    E2 = np.asarray(derived["E2"], dtype=float)
    theta = np.asarray(derived["theta_poloidal_from_z_deg"], dtype=float)
    Sz = np.asarray(derived["S_z"], dtype=float)
    Sr = np.asarray(derived["S_r"], dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=150)
    data = np.log10(np.maximum(E2, np.nanmax(E2) * 1e-8 if np.nanmax(E2) > 0 else 1e-30))
    im = ax.imshow(data, origin="lower", aspect="auto", extent=extent)
    ax.set_xlabel("z")
    ax.set_ylabel("r")
    ax.set_title(f"{name}: log10(|E|² proxy)")
    fig.colorbar(im, ax=ax, label="log10(E2)")
    fig.tight_layout()
    fig.savefig(cdir / "field_maps" / f"fieldmap_{name}_E2_log.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=150)
    im = ax.imshow(np.clip(theta, 0, 60), origin="lower", aspect="auto", extent=extent)
    ax.set_xlabel("z")
    ax.set_ylabel("r")
    ax.set_title(f"{name}: local poloidal Poynting angle θ from +z")
    fig.colorbar(im, ax=ax, label="θ, deg")
    fig.tight_layout()
    fig.savefig(cdir / "field_maps" / f"fieldmap_{name}_theta_poloidal.png")
    plt.close(fig)

    # Quiver map, downsampled.
    nr, nz = E2.shape
    step_r = max(1, nr // 36)
    step_z = max(1, nz // 54)
    Rg, Zg = np.meshgrid(r_grid, z_grid, indexing="ij")
    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    bg = np.log10(np.maximum(E2, np.nanmax(E2) * 1e-8 if np.nanmax(E2) > 0 else 1e-30))
    ax.imshow(bg, origin="lower", aspect="auto", extent=extent, alpha=0.55)
    ax.quiver(
        Zg[::step_r, ::step_z], Rg[::step_r, ::step_z],
        Sz[::step_r, ::step_z], Sr[::step_r, ::step_z],
        angles="xy", scale_units="xy", scale=None, width=0.0022
    )
    ax.set_xlabel("z")
    ax.set_ylabel("r")
    ax.set_title(f"{name}: meridional Poynting field (S_z,S_r)")
    fig.tight_layout()
    fig.savefig(cdir / "field_maps" / f"fieldmap_{name}_poynting_quiver.png")
    plt.close(fig)


def save_field_maps_from_dfts(sim, field_maps: Dict[str, Dict[str, object]], g: PHBGeometry, cfg: RunConfig, cdir: Path, source_records: Optional[List[Dict[str, object]]] = None) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    require_meep()
    all_metrics: Dict[str, object] = {}
    region_rows: List[Dict[str, object]] = []
    if not field_maps:
        return all_metrics, region_rows
    ensure_dir(cdir / "field_maps")
    for name, item in field_maps.items():
        spec = item["spec"]
        dft = item["dft"]
        comps = item["components"]
        arrays: Dict[str, np.ndarray] = {}
        for comp in comps:
            cname = component_name(comp)
            arrays[cname] = get_dft_map_array(sim, dft, comp, spec, cfg.resolution)
        # Align all arrays to the smallest common shape.
        nr = min(a.shape[0] for a in arrays.values())
        nz = min(a.shape[1] for a in arrays.values())
        arrays = {k: v[:nr, :nz] for k, v in arrays.items()}
        r_grid = np.linspace(float(spec["r_min"]), float(spec["r_max"]), nr)
        z_grid = np.linspace(float(spec["z_min"]), float(spec["z_max"]), nz)
        derived = compute_field_map_arrays(arrays)

        material_mask = diagnostic_material_mask_grid(r_grid, z_grid, g, cfg) if getattr(cfg, "save_material_maps", True) else None
        source_mask = None
        if getattr(cfg, "save_source_excluded_maps", True) and source_records:
            source_mask = source_exclusion_mask_grid(r_grid, z_grid, source_records, cfg)
            derived["source_exclusion_mask"] = source_mask.astype(np.uint8)
            derived["E2_source_excluded"] = np.where(source_mask, np.nan, derived["E2"])
        if material_mask is not None:
            derived["material_mask"] = material_mask.astype(np.int16)

        m = field_map_metrics(name, spec, derived, r_grid, z_grid, g, cfg)
        if source_mask is not None:
            m.update(source_excluded_E2_metrics(name, derived["E2"], source_mask, r_grid, z_grid, cfg))
        if material_mask is not None:
            total_cells = max(int(material_mask.size), 1)
            m[f"fieldmap_{name}_material_mask_air_fraction"] = float(np.sum(material_mask == 1) / total_cells)
            m[f"fieldmap_{name}_material_mask_metal_fraction"] = float(np.sum(material_mask == 2) / total_cells)
            m[f"fieldmap_{name}_material_mask_window_fraction"] = float(np.sum(material_mask == 3) / total_cells)
        all_metrics.update(m)
        row = {"fieldmap_region": name, **{k.replace(f"fieldmap_{name}_", ""): v for k, v in m.items()}}
        row["npz_file"] = f"field_maps/fieldmap_{name}.npz"
        region_rows.append(row)

        if cfg.field_map_npz:
            np.savez_compressed(
                cdir / "field_maps" / f"fieldmap_{name}.npz",
                r_grid=r_grid, z_grid=z_grid,
                metadata_json=json.dumps(spec, ensure_ascii=False),
                **arrays,
                **derived,
            )
        if cfg.field_map_csv_stride and cfg.field_map_csv_stride > 0:
            stride = max(1, int(cfg.field_map_csv_stride))
            Rg, Zg = np.meshgrid(r_grid, z_grid, indexing="ij")
            flat_rows: List[Dict[str, object]] = []
            for ii in range(0, nr, stride):
                for jj in range(0, nz, stride):
                    row_csv = {
                        "i_r": ii, "i_z": jj, "r": float(Rg[ii, jj]), "z": float(Zg[ii, jj]),
                        "E2": float(derived["E2"][ii, jj]),
                        "S_r": float(derived["S_r"][ii, jj]),
                        "S_phi": float(derived["S_phi"][ii, jj]),
                        "S_z": float(derived["S_z"][ii, jj]),
                        "theta_poloidal_from_z_deg": float(derived["theta_poloidal_from_z_deg"][ii, jj]),
                        "theta_3d_from_z_deg": float(derived["theta_3d_from_z_deg"][ii, jj]),
                        "phase_Ez_rad": float(derived["phase_Ez_rad"][ii, jj]),
                    }
                    if "material_mask" in derived:
                        row_csv["material_mask"] = int(derived["material_mask"][ii, jj])
                    if "source_exclusion_mask" in derived:
                        row_csv["source_exclusion_mask"] = int(derived["source_exclusion_mask"][ii, jj])
                        val = derived.get("E2_source_excluded", derived["E2"])[ii, jj]
                        row_csv["E2_source_excluded"] = float(val) if np.isfinite(val) else ""
                    flat_rows.append(row_csv)
            write_csv(cdir / "field_maps" / f"fieldmap_{name}_sampled_stride{stride}.csv", flat_rows)
        if cfg.field_map_png:
            plot_field_map_region(cdir, name, spec, r_grid, z_grid, derived)
            if material_mask is not None:
                plot_material_mask_region(cdir, name, spec, r_grid, z_grid, material_mask)
            if source_mask is not None and "E2_source_excluded" in derived:
                plot_source_excluded_E2_region(cdir, name, spec, r_grid, z_grid, derived["E2_source_excluded"])
    write_csv(cdir / "field_maps" / "fieldmap_region_metrics.csv", region_rows)
    return all_metrics, region_rows

def add_diagnostics(sim, g: PHBGeometry, cfg: RunConfig, rmax: float, zspan: float, z_ap: float, far_z: Optional[float], fcen: float):
    require_meep()
    comps = [mp.Er, mp.Ep, mp.Ez, mp.Hr, mp.Hp, mp.Hz]
    r_inner, r_outer, r_center = g.window_r_low(cfg.phb_type), g.window_r_high(cfg.phb_type), g.window_r_center(cfg.phb_type)
    width = max(g.window_width(cfg.phb_type), 1.0 / max(cfg.resolution, 1))
    dft_ap = sim.add_dft_fields(
        comps, fcen, 0, 1,
        center=mp.Vector3(r_center, 0, z_ap),
        size=mp.Vector3(width, 0, 0),
    )
    flux: Dict[str, object] = {}
    flux["aperture_R2_z"] = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(r_center, 0, z_ap), size=mp.Vector3(width, 0, 0), direction=mp.Z),
    )

    # Minimal escape bookkeeping. These are not used as hard angle criteria.
    dpml = float(cfg.dpml or 0.0)
    margin = max(0.12, 0.25 * wavelength_from_a(g.a, cfg.a_over_lambda))
    z_right = +0.5 * zspan - dpml - margin
    z_left = -0.5 * zspan + dpml + margin
    r_outer_boundary = rmax - dpml - margin
    if z_right > z_ap:
        flux["right_boundary_z"] = sim.add_flux(fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(0.5 * r_outer_boundary, 0, z_right), size=mp.Vector3(r_outer_boundary, 0, 0), direction=mp.Z))
    if z_left < -g.L:
        flux["left_boundary_z"] = sim.add_flux(fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(0.5 * r_outer_boundary, 0, z_left), size=mp.Vector3(r_outer_boundary, 0, 0), direction=mp.Z))
    if r_outer_boundary > g.R:
        flux["outer_boundary_r"] = sim.add_flux(fcen, 0, 1,
            mp.FluxRegion(center=mp.Vector3(r_outer_boundary, 0, 0), size=mp.Vector3(0, 0, zspan - 2.0 * dpml - 2.0 * margin), direction=mp.R))

    dft_far = None
    if cfg.enable_far_monitor and far_z is not None:
        # M2 far monitor: parallel to M1 and, by default, centered on the focal-axis radius R.
        r0, r1, _ = far_monitor_interval(g, cfg, rmax, zspan)
        if r1 > r0 and far_z < z_right:
            dft_far = sim.add_dft_fields(
                comps, fcen, 0, 1,
                center=mp.Vector3(0.5 * (r0 + r1), 0, far_z),
                size=mp.Vector3(r1 - r0, 0, 0),
            )
            flux["far_M2_z"] = sim.add_flux(
                fcen, 0, 1,
                mp.FluxRegion(center=mp.Vector3(0.5 * (r0 + r1), 0, far_z), size=mp.Vector3(r1 - r0, 0, 0), direction=mp.Z),
            )
    field_maps, field_map_plan = add_field_map_diagnostics(sim, g, cfg, rmax, zspan, z_ap, far_z, fcen)
    return dft_ap, dft_far, flux, field_maps, field_map_plan


def dft_array(sim, dft, component) -> np.ndarray:
    arr = sim.get_dft_array(dft, component, 0)
    return np.asarray(arr).squeeze().astype(complex).ravel()


def analyze_z_plane_profile(sim, dft, g: PHBGeometry, cfg: RunConfig, r_min: float, r_max: float, label: str) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    require_meep()
    Er = dft_array(sim, dft, mp.Er)
    Ep = dft_array(sim, dft, mp.Ep)
    Ez = dft_array(sim, dft, mp.Ez)
    Hr = dft_array(sim, dft, mp.Hr)
    Hp = dft_array(sim, dft, mp.Hp)
    Hz = dft_array(sim, dft, mp.Hz)
    n = min(len(Er), len(Ep), len(Ez), len(Hr), len(Hp), len(Hz))
    if n <= 0:
        raise RuntimeError(f"empty DFT profile for {label}")
    Er, Ep, Ez, Hr, Hp, Hz = [x[:n] for x in (Er, Ep, Ez, Hr, Hp, Hz)]
    dr = (r_max - r_min) / max(n, 1)
    r = r_min + (np.arange(n) + 0.5) * dr
    # Cylindrical complex Poynting components, time-averaged: 0.5 Re(E x H*)
    S_r = 0.5 * np.real(Ep * np.conj(Hz) - Ez * np.conj(Hp))
    S_phi = 0.5 * np.real(Ez * np.conj(Hr) - Er * np.conj(Hz))
    S_z = 0.5 * np.real(Er * np.conj(Hp) - Ep * np.conj(Hr))
    intensity = 0.5 * (np.abs(Er) ** 2 + np.abs(Ep) ** 2 + np.abs(Ez) ** 2)
    area_w = 2.0 * math.pi * np.maximum(r, 0.0) * abs(dr)
    forward = np.maximum(S_z, 0.0) * area_w
    backward = np.maximum(-S_z, 0.0) * area_w
    abs_z = np.abs(S_z) * area_w
    theta = np.degrees(np.arctan2(np.abs(S_r), np.maximum(S_z, 1e-300)))
    theta = np.where(S_z > 0, theta, 180.0)
    theta_3d = np.degrees(np.arctan2(np.sqrt(S_r * S_r + S_phi * S_phi), np.maximum(S_z, 1e-300)))
    theta_3d = np.where(S_z > 0, theta_3d, 180.0)
    fsum = float(np.sum(forward))
    bsum = float(np.sum(backward))
    abs_sum = float(np.sum(abs_z))
    theta95 = weighted_percentile(theta, forward, 0.95)
    theta90 = weighted_percentile(theta, forward, 0.90)
    theta50 = weighted_percentile(theta, forward, 0.50)
    theta_rms = float(math.sqrt(np.sum(forward * theta * theta) / max(np.sum(forward), 1e-300))) if fsum > 0 else float("nan")
    lobe5 = float(np.sum(forward[theta <= 5.0]) / max(fsum, 1e-300)) if fsum > 0 else float("nan")
    lobe10 = float(np.sum(forward[theta <= 10.0]) / max(fsum, 1e-300)) if fsum > 0 else float("nan")
    lobe20 = float(np.sum(forward[theta <= 20.0]) / max(fsum, 1e-300)) if fsum > 0 else float("nan")
    theta_density_5 = float(np.sum(forward[theta <= 5.0]) / max(math.radians(5.0), 1e-300)) if fsum > 0 else float("nan")
    area = math.pi * max(r_max * r_max - r_min * r_min, 0.0)
    main_lobe_area_density_5 = float(np.sum(forward[theta <= 5.0]) / max(area, 1e-300)) if fsum > 0 else float("nan")
    peak_i = int(np.nanargmax(intensity)) if len(intensity) else 0
    phase = np.angle(Ez)
    phase_rms = weighted_circular_phase_rms(phase, intensity * area_w)
    centroid_r = float(np.sum(r * forward) / max(fsum, 1e-300)) if fsum > 0 else float("nan")
    positive = forward > 0.0
    theta_max_sampled = float(np.nanmax(theta[positive])) if np.any(positive) else float("nan")
    theta_peak_forward = float(theta[int(np.nanargmax(forward))]) if fsum > 0 else float("nan")
    theta_peak_intensity = float(theta[peak_i]) if len(theta) else float("nan")

    metrics: Dict[str, object] = {
        f"{label}_points": int(n),
        f"{label}_r_min": float(r_min),
        f"{label}_r_max": float(r_max),
        f"{label}_r_center_nominal": float(0.5 * (r_min + r_max)),
        f"{label}_r_flux_centroid": centroid_r,
        f"{label}_forward_flux_proxy": fsum,
        f"{label}_backward_flux_proxy": bsum,
        f"{label}_abs_z_flux_proxy": abs_sum,
        f"{label}_backward_fraction_abs_z": bsum / max(abs_sum, 1e-300),
        f"{label}_theta95_deg": theta95,
        f"{label}_theta90_deg": theta90,
        f"{label}_theta50_deg": theta50,
        f"{label}_theta_3d95_deg": weighted_percentile(theta_3d, forward, 0.95),
        f"{label}_theta_3d90_deg": weighted_percentile(theta_3d, forward, 0.90),
        f"{label}_theta_3d50_deg": weighted_percentile(theta_3d, forward, 0.50),
        f"{label}_theta_rms_deg": theta_rms,
        f"{label}_theta_3d_rms_deg": float(math.sqrt(np.sum(forward * theta_3d * theta_3d) / max(np.sum(forward), 1e-300))) if fsum > 0 else float("nan"),
        f"{label}_azimuthal_poynting_abs_fraction": float(np.sum(np.abs(S_phi) * area_w) / max(np.sum(np.sqrt(S_r*S_r + S_phi*S_phi + S_z*S_z) * area_w), 1e-300)),
        f"{label}_theta_max_sampled_deg": theta_max_sampled,
        f"{label}_theta_at_max_density_deg": theta_peak_forward,
        f"{label}_theta_at_peak_intensity_deg": theta_peak_intensity,
        f"{label}_lobe_fraction_5deg": lobe5,
        f"{label}_lobe_fraction_10deg": lobe10,
        f"{label}_lobe_fraction_20deg": lobe20,
        f"{label}_theta_density_proxy_5deg": theta_density_5,
        f"{label}_main_lobe_area_density_5deg": main_lobe_area_density_5,
        f"{label}_intensity_peak": float(np.nanmax(intensity)),
        f"{label}_intensity_mean": float(np.nanmean(intensity)),
        f"{label}_intensity_peak_over_mean": float(np.nanmax(intensity) / max(np.nanmean(intensity), 1e-300)),
        f"{label}_phase_rms_Ez_rad": phase_rms,
        f"{label}_peak_r": float(r[peak_i]),
        f"{label}_peak_offset_from_slot_center": float(abs(r[peak_i] - 0.5 * (r_min + r_max))),
    }
    rows: List[Dict[str, object]] = []
    for i in range(n):
        rows.append({
            "i": i,
            "r": float(r[i]),
            "theta_deg": float(theta[i]),
            "theta_from_z_axis_deg": float(theta[i]),
            "S_r": float(S_r[i]),
            "S_phi": float(S_phi[i]),
            "S_z": float(S_z[i]),
            "theta_3d_from_z_axis_deg": float(theta_3d[i]),
            "poynting_azimuthal_fraction": float(abs(S_phi[i]) / max(math.sqrt(S_r[i]*S_r[i] + S_phi[i]*S_phi[i] + S_z[i]*S_z[i]), 1e-300)),
            "forward_weight": float(forward[i]),
            "intensity": float(intensity[i]),
            "phase_Ez_rad": float(phase[i]),
        })
    return metrics, rows


def classify_candidate(row: Dict[str, object], cfg: RunConfig) -> str:
    th = safe_float(row.get("aperture_theta95_deg"))
    lobe5 = safe_float(row.get("aperture_lobe_fraction_5deg"))
    flux = safe_float(row.get("aperture_forward_flux_proxy"))
    kappa = safe_float(row.get("kappa_R2_escape_proxy"))
    back = safe_float(row.get("aperture_backward_fraction_abs_z"))
    if not math.isfinite(th):
        return "INVALID_NO_APERTURE_THETA95"
    if not math.isfinite(flux) or flux <= 0:
        return "INVALID_NO_FORWARD_R2_FLUX"
    if math.isfinite(back) and back > 0.50:
        return "CHECK_STRONG_BACKWARD_COMPONENT"
    angle_ok = th <= cfg.narrow_theta95_deg
    lobe_ok = math.isfinite(lobe5) and lobe5 >= cfg.narrow_lobe5_min
    flux_ok = flux >= cfg.useful_flux_min
    kappa_ok = (not math.isfinite(kappa)) or kappa >= cfg.kappa_min
    if angle_ok and lobe_ok and flux_ok and kappa_ok:
        return "CANDIDATE_NARROW_R2_OUTPUT"
    if angle_ok and lobe_ok and not flux_ok:
        return "CHECK_NARROW_BUT_WEAK_R2_FLUX"
    if angle_ok and not lobe_ok:
        return "CHECK_LOW_THETA95_BUT_POOR_5DEG_DENSITY"
    return "FAIL_WIDE_R2_OUTPUT"


def plot_profile(out_png: Path, profile_rows: List[Dict[str, object]], metrics: Dict[str, object], title: str) -> None:
    if plt is None or not profile_rows:
        return
    ensure_dir(out_png.parent)
    r = np.array([safe_float(x.get("r")) for x in profile_rows])
    theta = np.array([safe_float(x.get("theta_deg")) for x in profile_rows])
    fw = np.array([safe_float(x.get("forward_weight")) for x in profile_rows])
    intensity = np.array([safe_float(x.get("intensity")) for x in profile_rows])

    fig, ax1 = plt.subplots(figsize=(8, 4.5), dpi=150)
    ax1.plot(r, theta, lw=1.5, label="local Poynting angle, deg")
    ax1.set_xlabel("r across R2 aperture")
    ax1.set_ylabel("local angle from +z, deg")
    ax1.grid(True, alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(r, fw / max(np.nanmax(fw), 1e-300), lw=1.2, ls="--", label="normalized forward flux density")
    ax2.plot(r, intensity / max(np.nanmax(intensity), 1e-300), lw=1.0, ls=":", label="normalized |E|²")
    ax2.set_ylabel("normalized density")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper right")
    ax1.set_title(title)
    note = (
        f"theta95={safe_float(metrics.get('aperture_theta95_deg')):.4g} deg; "
        f"lobe5={safe_float(metrics.get('aperture_lobe_fraction_5deg')):.3g}; "
        f"flux={safe_float(metrics.get('aperture_forward_flux_proxy')):.3e}; "
        f"status={metrics.get('candidate_status')}"
    )
    ax1.text(0.01, 0.01, note, transform=ax1.transAxes, fontsize=8, va="bottom")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)



def plot_publication_exit_cone(out_png: Path, profile_rows: List[Dict[str, object]], metrics: Dict[str, object], g: PHBGeometry, cfg: RunConfig, title: str) -> None:
    """Publication-safe visualization of the axisymmetric PHB output cone.

    Convention fixed for PHB papers:
      * the physical longitudinal direction is +z, i.e. F−→F+;
      * in the meridional drawing the relevant focal-ring guide is r=+R
        for the upper branch and r=-R for the mirrored lower branch;
      * the central line r=0 is only the mechanical rotation axis, not the
        focal-ring output axis;
      * theta is the half-angle relative to +z / F−→F+ and is drawn as two
        symmetric envelopes around the focal-ring guide, not as a one-sided fan.
    """
    if plt is None or not profile_rows:
        return
    ensure_dir(out_png.parent)
    z = np.linspace(-g.L, g.L, 1200)
    rwall = np.full_like(z, np.nan, dtype=float)
    for i, zz in enumerate(z):
        az = abs(float(zz))
        if az >= g.a:
            rwall[i] = g.rho_horn(az)
        else:
            if cfg.phb_type == "open":
                rwall[i] = np.nan
            else:
                rwall[i] = phb_central_wall_radius(float(zz), g, cfg.phb_type)

    z0, z1_eff, z1_exact, widened = effective_slot_interval(g, cfg)
    z_ap = safe_float(metrics.get("aperture_monitor_z"))
    if not math.isfinite(z_ap):
        z_ap = z1_eff + cfg.aperture_offset_cells / max(cfg.resolution, 1)
    r0, r1 = g.window_r_low(cfg.phb_type), g.window_r_high(cfg.phb_type)
    focal_r = g.R

    theta = np.array([safe_float(x.get("theta_from_z_axis_deg", x.get("theta_deg"))) for x in profile_rows], dtype=float)
    fw = np.maximum(np.array([safe_float(x.get("forward_weight")) for x in profile_rows], dtype=float), 0.0)
    mask = np.isfinite(theta) & np.isfinite(fw) & (fw > 0.0)
    if not np.any(mask):
        return
    theta, fw = theta[mask], fw[mask]
    wn = fw / max(float(np.nanmax(fw)), 1e-300)

    theta95 = safe_float(metrics.get("aperture_theta95_deg"))
    theta90 = safe_float(metrics.get("aperture_theta90_deg"))
    theta50 = safe_float(metrics.get("aperture_theta50_deg"))
    theta_peak = safe_float(metrics.get("aperture_theta_at_max_density_deg"))
    if not math.isfinite(theta_peak):
        theta_peak = float(theta[int(np.nanargmax(fw))])
    theta_max = safe_float(metrics.get("aperture_theta_max_sampled_deg"))
    if not math.isfinite(theta_max):
        theta_max = float(np.nanmax(theta))

    fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
    ax.plot(z, rwall, lw=2.0, label="PHB wall / closure")
    ax.plot(z, -rwall, lw=2.0)

    # R2 output window and M1 aperture monitor, mirrored for the full meridional drawing.
    if getattr(cfg, "window_model", "vertical_shift") == "vertical_shift":
        ax.fill_between([z0, z1_eff], [r0, r0], [r1, r1], alpha=0.30, label="right shifted vertical R2 window [R+shift-W/2,R+shift+W/2]")
        ax.fill_between([z0, z1_eff], [-r1, -r1], [-r0, -r0], alpha=0.18)
        ax.plot([g.a, g.a], [r0, r1], lw=4, alpha=0.75, label="center-shifted R2 opening: [R+shift-W/2, R+shift+W/2]")
        ax.plot([g.a, g.a], [-r1, -r0], lw=4, alpha=0.45)
    else:
        zzs = np.linspace(z0, z1_eff, 80)
        rs = np.array([g.rho_horn(float(x)) for x in zzs])
        ax.plot(zzs, rs, lw=4, alpha=0.70, label="R2 window on horn, effective segment")
        ax.plot(zzs, -rs, lw=4, alpha=0.70)
        if cfg.phb_type in ("cylinder", "halfring"):
            ax.plot([g.a, g.a], [g.R, r1], lw=4, alpha=0.70, label="upper part of corrected window")
            ax.plot([g.a, g.a], [-r1, -g.R], lw=4, alpha=0.70)
        elif cfg.phb_type == "open":
            ax.plot([-g.a, -g.a], [g.R, g.R + 0.5 * g.R2], lw=3, alpha=0.55, label="left add-on wall")
            ax.plot([-g.a, -g.a], [-(g.R + 0.5 * g.R2), -g.R], lw=3, alpha=0.55)
    ax.plot([z_ap, z_ap], [r0, r1], lw=4, label="M1 aperture monitor")
    ax.plot([z_ap, z_ap], [-r1, -r0], lw=4)

    # Correct PHB focal-ring axes: horizontal guides through external focal rings.
    ax.axhline(+focal_r, lw=2.0, label="upper focal-ring axis F−→F+ / +z")
    ax.axhline(-focal_r, lw=2.0, label="lower mirrored focal-ring axis")
    ax.axhline(0.0, lw=0.9, ls=":", label="mechanical rotation axis r=0")
    ax.axvline(-g.c, ls=":", lw=0.9)
    ax.axvline(+g.c, ls=":", lw=0.9)
    ax.scatter([-g.c, +g.c], [focal_r, focal_r], marker="x", s=70, label="F−, F+ on upper focal ring")
    ax.scatter([-g.c, +g.c], [-focal_r, -focal_r], marker="x", s=70)

    ray_len = max(0.8, 0.18 * (g.L + 1.0))

    # Weighted angular density inside the cone.  Draw symmetric lines around the
    # focal-ring guide, not a one-sided fan from +r only.
    for th, weight_norm in zip(theta, wn):
        rad = math.radians(float(th))
        dz = ray_len * math.cos(rad)
        dr = ray_len * math.sin(rad)
        lw = 0.4 + 1.7 * float(weight_norm)
        alpha = 0.10 + 0.55 * float(weight_norm)
        for sgn in (+1.0, -1.0):
            axis_r = sgn * focal_r
            ax.plot([z_ap, z_ap + dz], [axis_r, axis_r + dr], lw=lw, alpha=alpha)
            ax.plot([z_ap, z_ap + dz], [axis_r, axis_r - dr], lw=lw, alpha=alpha)

    # Cone half-angle envelopes for theta95 and theta_max.  If they coincide,
    # draw one guide to avoid misleading overprinted labels.
    cone_guides: List[Tuple[float, str, str]] = []
    if math.isfinite(theta95):
        cone_guides.append((theta95, "--", "θ95"))
    if math.isfinite(theta_max):
        if math.isfinite(theta95) and abs(theta_max - theta95) < 0.15:
            cone_guides[-1] = (theta95, "--", "θ95≈θmax")
        else:
            cone_guides.append((theta_max, ":", "θmax"))
    for th, ls, name in cone_guides:
        rad = math.radians(th)
        dz = 1.12 * ray_len * math.cos(rad)
        dr = 1.12 * ray_len * math.sin(rad)
        for sgn in (+1.0, -1.0):
            axis_r = sgn * focal_r
            ax.plot([z_ap, z_ap + dz], [axis_r, axis_r + dr], ls=ls, lw=2.1)
            ax.plot([z_ap, z_ap + dz], [axis_r, axis_r - dr], ls=ls, lw=2.1)
        ax.text(z_ap + dz + 0.02, focal_r + dr, f"{name}={th:.2f}°", fontsize=8)
        ax.text(z_ap + dz + 0.02, -focal_r - dr, f"{name}={th:.2f}°", fontsize=8)

    ymax_wall = np.nanmax(np.abs(rwall[np.isfinite(rwall)])) if np.any(np.isfinite(rwall)) else g.R
    ymax = max(focal_r + 1.08 * ray_len, ymax_wall + 0.20, r1 + 1.05 * ray_len)
    ax.set_xlim(-g.L - 0.1, z_ap + 1.25 * ray_len + 0.1)
    ax.set_ylim(-ymax, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("axial coordinate z: direction F−→F+ / +z")
    ax.set_ylabel("meridional radius r; upper/lower focal-ring axes shown at ±R")
    ax.set_title(title)
    note = (
        f"θ is the cone half-angle from the horizontal F−→F+ / +z direction.\n"
        f"The cone is centered on the focal-ring guides r=±R, not on a tilted one-sided ray.\n"
        f"θ50={theta50:.3g}°, θ90={theta90:.3g}°, θ95={theta95:.3g}°, "
        f"θpeak={theta_peak:.3g}°, θmax={theta_max:.3g}°."
    )
    ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=8, va="bottom")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_publication_angular_density(out_png: Path, profile_rows: List[Dict[str, object]], metrics: Dict[str, object], title: str, prefix: str = "aperture") -> None:
    """Publication-safe angular energy-density plot D(theta)."""
    if plt is None or not profile_rows:
        return
    ensure_dir(out_png.parent)
    theta = np.array([safe_float(x.get("theta_from_z_axis_deg", x.get("theta_deg"))) for x in profile_rows], dtype=float)
    fw = np.maximum(np.array([safe_float(x.get("forward_weight")) for x in profile_rows], dtype=float), 0.0)
    mask = np.isfinite(theta) & np.isfinite(fw) & (fw > 0.0)
    if not np.any(mask):
        return
    theta, fw = theta[mask], fw[mask]

    theta50 = safe_float(metrics.get(f"{prefix}_theta50_deg"))
    theta90 = safe_float(metrics.get(f"{prefix}_theta90_deg"))
    theta95 = safe_float(metrics.get(f"{prefix}_theta95_deg"))
    theta_max = safe_float(metrics.get(f"{prefix}_theta_max_sampled_deg"))
    if not math.isfinite(theta_max):
        theta_max = float(np.nanmax(theta))
    theta_peak = safe_float(metrics.get(f"{prefix}_theta_at_max_density_deg"))
    if not math.isfinite(theta_peak):
        theta_peak = float(theta[int(np.nanargmax(fw))])

    if len(theta) >= 3 and float(np.nanmax(theta)) > float(np.nanmin(theta)):
        nbins = max(6, min(16, len(theta)))
        bins = np.linspace(float(np.nanmin(theta)), float(np.nanmax(theta)), nbins + 1)
        hist, edges = np.histogram(theta, bins=bins, weights=fw)
        centers = 0.5 * (edges[:-1] + edges[1:])
        width = np.diff(edges)
    else:
        centers = theta
        hist = fw
        width = np.full_like(theta, 0.25, dtype=float)
    hist_norm = hist / max(float(np.nanmax(hist)), 1e-300)

    order = np.argsort(theta)
    theta_s = theta[order]
    cum = np.cumsum(fw[order])
    cum = cum / max(float(cum[-1]), 1e-300)

    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=160)
    ax.bar(centers, hist_norm, width=width, align="center", alpha=0.75, label="normalized angular energy density D(θ)")
    ax.plot(theta_s, cum, marker="o", ms=3, lw=1.3, label="cumulative forward energy")
    for value, ls, lab in [(theta50, ":", "θ50"), (theta90, "-.", "θ90"), (theta95, "--", "θ95"), (theta_max, "-", "θmax")]:
        if math.isfinite(value):
            ax.axvline(value, ls=ls, lw=1.2)
            ax.text(value, 1.02, f"{lab}={value:.2f}°", rotation=90, va="bottom", fontsize=8)
    if math.isfinite(theta_peak):
        ax.scatter([theta_peak], [1.0], marker="x", s=60, label=f"θpeak={theta_peak:.2f}°")
    ax.set_ylim(0, 1.13)
    ax.set_xlabel("local output angle θ from the F−→F+ / +z axis, degrees")
    ax.set_ylabel("normalized density / cumulative energy")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


# -----------------------------
# Run one case and stages
# -----------------------------

def case_name(phb_type: str, m: int, R2: float, resolution: int, aol: float) -> str:
    def f(x: float) -> str:
        return (f"{x:.6g}".replace("-", "m").replace(".", "p"))
    return f"{phb_type}_m{m}_R2_{f(R2)}_res{resolution}_aol_{f(aol)}"


def run_one_case(args, phb_type: str, m: int, R2: float) -> Dict[str, object]:
    g = PHBGeometry(a=args.a, b=args.b, R=args.R, R2=R2, window_offset=args.window_offset, window_offset_mode=args.window_offset_mode)
    cfg = RunConfig(
        stage=args.stage, outroot=args.outroot, phb_type=phb_type, m=m, resolution=args.resolution,
        a_over_lambda=args.a_over_lambda, min_a_over_lambda=args.min_a_over_lambda,
        dpml_over_lambda=args.dpml_over_lambda, dpml=args.dpml, wall_thickness=args.wall_thickness,
        window_model=args.window_model, window_offset=args.window_offset, window_offset_mode=args.window_offset_mode, vertical_window_thickness=args.vertical_window_thickness,
        source_components=args.source_components, source_mode=args.source_mode, nsrc=args.nsrc,
        seed=args.seed, fwidth_frac=args.fwidth_frac, after_sources=args.after_sources,
        skip_existing=args.skip_existing, allow_underresolved=args.allow_underresolved,
        min_radial_slot_cells=args.min_radial_slot_cells, min_axial_slot_cells=args.min_axial_slot_cells,
        aperture_offset_cells=args.aperture_offset_cells,
        enable_far_monitor=args.enable_far_monitor, far_distance=args.far_distance,
        far_capture_angle_deg=args.far_capture_angle_deg, far_monitor_mode=args.far_monitor_mode, far_zone_safety=args.far_zone_safety,
        allow_near_field_far_monitor=args.allow_near_field_far_monitor,
        narrow_theta95_deg=args.narrow_theta95_deg, narrow_lobe5_min=args.narrow_lobe5_min,
        useful_flux_min=args.useful_flux_min, kappa_min=args.kappa_min, min_free_gb=args.min_free_gb,
        archive=args.archive, stop_on_error=args.stop_on_error,
        save_field_maps=args.save_field_maps,
        field_map_regions=args.field_map_regions,
        field_map_components=args.field_map_components,
        field_map_r_pad=args.field_map_r_pad,
        field_map_z_half=args.field_map_z_half,
        field_map_corridor_half_width=args.field_map_corridor_half_width,
        field_map_interior_rmax=args.field_map_interior_rmax,
        field_map_max_points_per_region=args.field_map_max_points_per_region,
        field_map_csv_stride=args.field_map_csv_stride,
        field_map_png=args.field_map_png,
        field_map_npz=args.field_map_npz,
        save_material_maps=args.save_material_maps,
        save_source_excluded_maps=args.save_source_excluded_maps,
        source_exclusion_cells=args.source_exclusion_cells,
        m2_extra_flux_warning_ratio=args.m2_extra_flux_warning_ratio,
    )
    outroot = Path(args.outroot)
    cname = case_name(phb_type, m, R2, args.resolution, args.a_over_lambda)
    cdir = ensure_dir(outroot / "cases" / cname)
    metrics_path = cdir / "metrics.json"
    if args.skip_existing and metrics_path.exists():
        row = json.loads(metrics_path.read_text(encoding="utf-8"))
        row["skipped_existing"] = True
        return row

    event_base = {
        "time": now_iso(), "case": cname, "phb_type": phb_type, "m": m, "R2_abs": R2, "window_offset": args.window_offset, "window_offset_mode": args.window_offset_mode,
        "resolution": args.resolution, "a_over_lambda": args.a_over_lambda,
    }
    append_jsonl(outroot / "logs" / "runs.jsonl", {**event_base, "event": "start"})

    row: Dict[str, object] = {
        "case": cname,
        "phb_type": phb_type,
        "m": m,
        "R2_abs": R2,
        "a": args.a,
        "b": args.b,
        "R": args.R,
        "resolution": args.resolution,
        "a_over_lambda": args.a_over_lambda,
        "lambda": wavelength_from_a(args.a, args.a_over_lambda),
        "started_at": now_iso(),
        "window_model": args.window_model,
        "vertical_window_thickness_requested": args.vertical_window_thickness if args.vertical_window_thickness is not None else "auto",
    }
    valid = validate_case(g, cfg)
    row.update({k: v for k, v in valid.items() if k not in ("ok", "problems", "warnings")})
    # v46 audit fields: M1 is intentionally shifted above the focal-axis line.
    row["M1_r_min"] = g.window_r_low(cfg.phb_type)
    row["M1_r_max"] = g.window_r_high(cfg.phb_type)
    row["M1_r_center"] = g.window_r_center(cfg.phb_type)
    row["M1_focal_axis_r"] = g.R
    row["window_offset"] = g.window_offset
    row["window_offset_mode"] = g.window_offset_mode
    row["M1_lower_edge_on_focal_axis"] = abs(g.window_r_low(cfg.phb_type) - g.R) <= 1e-12
    row["M1_lower_edge_offset_from_focal_axis"] = g.window_r_low(cfg.phb_type) - g.R
    row["M1_center_offset_from_focal_axis"] = g.window_r_center(cfg.phb_type) - g.R
    row["M1_is_shifted_above_focal_axis"] = g.window_r_low(cfg.phb_type) >= g.R - 1e-12
    row["M1_is_symmetric_about_focal_axis"] = False
    row["validation_warnings"] = "; ".join(valid["warnings"])
    row["validation_problems"] = "; ".join(valid["problems"])

    plot_geometry(cdir / "geometry_monitor.png", g, cfg)

    if not bool(valid["ok"]):
        row["candidate_status"] = "INVALID_NUMERICS_OR_GEOMETRY"
        write_json(metrics_path, row)
        append_jsonl(outroot / "logs" / "runs.jsonl", {**event_base, "event": "invalid", "problems": valid["problems"]})
        if args.stop_on_error:
            raise RuntimeError(row["validation_problems"])
        return row

    if free_gb(outroot) < args.min_free_gb:
        row["candidate_status"] = "INVALID_LOW_DISK_SPACE"
        row["validation_problems"] = f"free disk {free_gb(outroot):.2f} GB < requested {args.min_free_gb:.2f} GB"
        write_json(metrics_path, row)
        append_jsonl(outroot / "logs" / "runs.jsonl", {**event_base, "event": "invalid_low_disk"})
        if args.stop_on_error:
            raise RuntimeError(row["validation_problems"])
        return row

    t0 = time.time()
    sim = None
    try:
        sim, rmax, zspan, z_ap, far_z, fcen, source_records = make_simulation(g, cfg)
        dft_ap, dft_far, flux, field_maps, field_map_plan = add_diagnostics(sim, g, cfg, rmax, zspan, z_ap, far_z, fcen)
        row.update({"rmax": rmax, "zspan": zspan, "aperture_monitor_z": z_ap, "far_monitor_z": far_z if far_z is not None else ""})
        row["source_count"] = int(len(source_records))
        row["source_coordinates_available"] = bool(source_records)
        row["source_coordinates_csv"] = "sources/source_coordinates.csv"
        row["source_coordinates_json"] = "sources/source_coordinates.json"
        write_csv(cdir / "sources" / "source_coordinates.csv", source_records)
        write_json(cdir / "sources" / "source_coordinates.json", source_records)
        row["field_maps_requested"] = bool(args.save_field_maps)
        row["field_map_regions_requested"] = args.field_map_regions
        row["field_map_plan_json"] = json.dumps(field_map_plan, ensure_ascii=False)
        skipped_maps = [x for x in field_map_plan if x.get("skipped")]
        if skipped_maps:
            row["field_map_skipped_regions"] = ",".join(str(x.get("name")) for x in skipped_maps)
            row["field_map_skip_reasons"] = "; ".join(str(x.get("skip_reason")) for x in skipped_maps)
        append_jsonl(outroot / "logs" / "runs.jsonl", {**event_base, "event": "meep_run_begin", "rmax": rmax, "zspan": zspan})
        sim.run(until_after_sources=float(args.after_sources))

        ap_metrics, ap_profile = analyze_z_plane_profile(sim, dft_ap, g, cfg, g.window_r_low(cfg.phb_type), g.window_r_high(cfg.phb_type), "aperture")
        row.update(ap_metrics)
        flux_vals: Dict[str, float] = {}
        for name, obj in flux.items():
            try:
                vals = mp.get_fluxes(obj)
                flux_vals[name] = float(vals[0]) if vals else float("nan")
            except Exception:
                flux_vals[name] = float("nan")
        row.update({f"flux_monitor_{k}": v for k, v in flux_vals.items()})
        # Selective output coefficient: use the MEEP aperture flux monitor as the
        # primary numerator so numerator and escape monitors have the same origin.
        # The local DFT/Poynting aperture flux is kept as a fallback and as the
        # primary source for the local angle metrics.
        ap_flux_local = max(0.0, safe_float(row.get("aperture_forward_flux_proxy")))
        ap_flux_monitor_raw = safe_float(flux_vals.get("aperture_R2_z"))
        ap_flux_monitor = max(0.0, ap_flux_monitor_raw) if math.isfinite(ap_flux_monitor_raw) else float("nan")
        if math.isfinite(ap_flux_monitor) and ap_flux_monitor > 0.0:
            ap_flux_for_kappa = ap_flux_monitor
            row["kappa_numerator_source"] = "meep_flux_monitor_aperture_R2_z"
        else:
            ap_flux_for_kappa = ap_flux_local
            row["kappa_numerator_source"] = "local_dft_poynting_fallback"
        row["aperture_flux_for_kappa"] = ap_flux_for_kappa
        # Only true escape-boundary monitors are used for this proxy.
        # Diagnostic internal slices, especially far_M2_z, must not enter the
        # denominator because they can re-measure the same aperture power.
        escape_keys = ("left_boundary_z", "right_boundary_z", "outer_boundary_r")
        escape_abs = sum(abs(flux_vals.get(k, float("nan"))) for k in escape_keys if math.isfinite(safe_float(flux_vals.get(k))))
        row["escape_abs_flux_proxy"] = escape_abs
        row["escape_flux_keys_for_kappa"] = ",".join(k for k in escape_keys if k in flux_vals)
        if "far_M2_z" in flux_vals:
            row["far_M2_z_excluded_from_kappa"] = True
        row["kappa_R2_escape_proxy"] = ap_flux_for_kappa / max(ap_flux_for_kappa + escape_abs, 1e-300)

        if dft_far is not None and far_z is not None:
            fr0, fr1, fref = far_monitor_interval(g, cfg, rmax, zspan)
            far_metrics, far_profile = analyze_z_plane_profile(sim, dft_far, g, cfg, fr0, fr1, "far_M2")
            far_geom_metrics = add_far_geometric_angle_metrics(far_profile, fref, max(far_z - z_ap, 1e-300), "far_M2")
            row.update(far_metrics)
            row.update(far_geom_metrics)
            row["far_M2_r_min"] = fr0
            row["far_M2_r_max"] = fr1
            row["far_M2_r_center"] = 0.5 * (fr0 + fr1)
            row["far_M2_reference_r"] = fref
            row["far_M2_focal_axis_r"] = g.R
            row["far_M2_shifted_window_center_r"] = g.window_r_center(cfg.phb_type)
            row["far_M2_center_offset_from_focal_axis"] = 0.5 * (fr0 + fr1) - g.R
            row["far_M2_reference_offset_from_focal_axis"] = fref - g.R
            row["far_M2_is_shifted_above_focal_axis"] = fref >= g.R
            row["far_M2_is_symmetric_about_focal_axis"] = False
            row["far_M2_mode"] = cfg.far_monitor_mode
            write_csv(cdir / "far_M2_profile.csv", far_profile)
            plot_profile(cdir / "far_M2_profile.png", far_profile, {**row, "aperture_theta95_deg": row.get("far_M2_theta95_deg"), "aperture_lobe_fraction_5deg": row.get("far_M2_lobe_fraction_5deg"), "aperture_forward_flux_proxy": row.get("far_M2_forward_flux_proxy"), "candidate_status": "far_M2"}, f"M2 far monitor profile: {cname}")
            plot_publication_angular_density(
                cdir / "far_M2_publication_angular_density.png", far_profile, row,
                f"M2 local angular energy density: {cname}", prefix="far_M2"
            )

            # v47.2: warn when M2 collects more forward flux than the selected M1 aperture.
            # This means M2 cannot be interpreted as the propagation of only the M1-selected beam.
            m2_flux_monitor = safe_float(flux_vals.get("far_M2_z"))
            m2_flux_forward = safe_float(row.get("far_M2_forward_flux_proxy"))
            m1_flux_monitor = safe_float(flux_vals.get("aperture_R2_z"))
            m1_flux_forward = safe_float(row.get("aperture_forward_flux_proxy"))
            ratio_monitor = abs(m2_flux_monitor) / max(abs(m1_flux_monitor), 1e-300) if math.isfinite(m2_flux_monitor) and math.isfinite(m1_flux_monitor) else float("nan")
            ratio_forward = m2_flux_forward / max(m1_flux_forward, 1e-300) if math.isfinite(m2_flux_forward) and math.isfinite(m1_flux_forward) else float("nan")
            row["far_M2_to_M1_flux_ratio_monitor"] = ratio_monitor
            row["far_M2_to_M1_flux_ratio_forward_proxy"] = ratio_forward
            row["far_M2_extra_flux_warning_ratio"] = float(cfg.m2_extra_flux_warning_ratio)
            row["far_M2_contains_extra_field_not_only_aperture_beam"] = bool(
                (math.isfinite(ratio_monitor) and ratio_monitor > cfg.m2_extra_flux_warning_ratio)
                or (math.isfinite(ratio_forward) and ratio_forward > cfg.m2_extra_flux_warning_ratio)
            )

        if args.save_field_maps:
            fmap_metrics, fmap_rows = save_field_maps_from_dfts(sim, field_maps, g, cfg, cdir, source_records)
            row.update(fmap_metrics)
            row["field_map_regions_saved"] = ",".join(str(x.get("fieldmap_region")) for x in fmap_rows)
            row["field_map_region_metrics_csv"] = "field_maps/fieldmap_region_metrics.csv" if fmap_rows else ""
            add_source_proximity_metrics(row, source_records, cfg)

        row["candidate_status"] = classify_candidate(row, cfg)
        row["runtime_seconds"] = time.time() - t0
        row["completed_at"] = now_iso()
        write_csv(cdir / "aperture_profile.csv", ap_profile)
        plot_profile(cdir / "aperture_profile.png", ap_profile, row, f"Aperture R2 output: {cname}")
        plot_publication_exit_cone(
            cdir / "publication_exit_cone.png", ap_profile, row, g, cfg,
            f"PHB publication exit cone: {cname}"
        )
        plot_publication_angular_density(
            cdir / "publication_angular_density.png", ap_profile, row,
            f"PHB angular energy density: {cname}", prefix="aperture"
        )
        write_json(metrics_path, row)
        append_jsonl(outroot / "logs" / "runs.jsonl", {**event_base, "event": "done", "status": row["candidate_status"], "theta95": row.get("aperture_theta95_deg")})
        return row
    except Exception as e:
        row["candidate_status"] = "ERROR_EXCEPTION"
        row["error"] = repr(e)
        row["traceback"] = traceback.format_exc()
        row["runtime_seconds"] = time.time() - t0
        write_json(metrics_path, row)
        append_jsonl(outroot / "logs" / "runs.jsonl", {**event_base, "event": "error", "error": repr(e)})
        if args.stop_on_error:
            raise
        return row
    finally:
        if sim is not None:
            try:
                sim.reset_meep()
            except Exception:
                pass


def collect_existing_metrics(outroot: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for p in sorted((outroot / "cases").glob("*/metrics.json")):
        try:
            rows.append(json.loads(p.read_text(encoding="utf-8")))
        except Exception:
            pass
    return rows


def write_summary(outroot: Path, rows: List[Dict[str, object]], args) -> None:
    # Stable field order for the important first-stage columns.
    primary = [
        "case", "candidate_status", "phb_type", "m", "R2_abs", "a", "b", "R", "lambda", "a_over_lambda", "resolution",
        "window_model", "vertical_window_thickness_used", "vertical_window_center_z",
        "aperture_theta95_deg", "aperture_theta90_deg", "aperture_theta50_deg", "aperture_theta_rms_deg",
        "aperture_theta_max_sampled_deg", "aperture_theta_at_max_density_deg", "aperture_theta_at_peak_intensity_deg",
        "aperture_lobe_fraction_5deg", "aperture_lobe_fraction_10deg", "aperture_lobe_fraction_20deg",
        "aperture_forward_flux_proxy", "aperture_main_lobe_area_density_5deg", "aperture_theta_density_proxy_5deg",
        "aperture_phase_rms_Ez_rad", "aperture_backward_fraction_abs_z",
        "far_M2_geom_theta95_deg", "far_M2_geom_lobe_fraction_5deg", "far_M2_geom_density_proxy_5deg",
        "far_M2_theta95_deg", "far_M2_theta90_deg", "far_M2_theta50_deg", "far_M2_theta_max_sampled_deg",
        "far_M2_theta_at_max_density_deg", "far_M2_lobe_fraction_5deg", "far_M2_forward_flux_proxy",
        "M1_r_min", "M1_r_max", "M1_r_center", "M1_focal_axis_r", "M1_lower_edge_on_focal_axis", "M1_center_offset_from_focal_axis", "M1_is_shifted_above_focal_axis", "M1_is_symmetric_about_focal_axis",
        "far_M2_reference_r", "far_M2_r_min", "far_M2_r_max", "far_M2_r_center", "far_M2_focal_axis_r", "far_M2_center_offset_from_focal_axis", "far_M2_reference_offset_from_focal_axis", "far_M2_is_symmetric_about_focal_axis", "far_M2_mode",
        "flux_monitor_aperture_R2_z", "flux_monitor_far_M2_z", "aperture_flux_for_kappa", "kappa_R2_escape_proxy", "kappa_numerator_source",
        "slot_widened", "slot_radial_cells", "slot_exact_axial_cells", "slot_effective_axial_cells",
        "validation_warnings", "validation_problems", "runtime_seconds",
    ]
    # append extra keys without losing them
    extra: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in primary and k not in extra:
                extra.append(k)
    write_csv(outroot / "PHB_V46_UNIVERSAL_R2_APERTURE_SUMMARY.csv", rows, primary + extra)

    good = [r for r in rows if str(r.get("candidate_status")) == "CANDIDATE_NARROW_R2_OUTPUT"]
    checks = [r for r in rows if str(r.get("candidate_status", "")).startswith("CHECK")]
    valid = [r for r in rows if not str(r.get("candidate_status", "")).startswith("INVALID") and not str(r.get("candidate_status", "")).startswith("ERROR")]
    by_theta = sorted(valid, key=lambda r: safe_float(r.get("aperture_theta95_deg")))[:20]
    by_density = sorted(valid, key=lambda r: -safe_float(r.get("aperture_main_lobe_area_density_5deg"), -1.0))[:20]

    lines: List[str] = []
    lines.append("PHB v46 first aperture verification — short Russian report")
    lines.append(f"Created: {now_iso()}")
    lines.append("")
    lines.append("Главный критерий первого этапа:")
    lines.append("  локальный угол выхода энергии из R2-щели по вектору Пойнтинга на апертурном мониторе")
    lines.append("  aperture_theta95_deg + плотность/доля энергии в ±5°, ±10°, ±20°.")
    lines.append("")
    lines.append(f"Total cases: {len(rows)}")
    lines.append(f"Narrow candidates: {len(good)}")
    lines.append(f"CHECK cases: {len(checks)}")
    lines.append("")
    lines.append("TOP by smallest aperture theta95:")
    for r in by_theta:
        lines.append(
            f"  theta95={safe_float(r.get('aperture_theta95_deg')):.6g} deg; "
            f"lobe5={safe_float(r.get('aperture_lobe_fraction_5deg')):.4g}; "
            f"density5={safe_float(r.get('aperture_main_lobe_area_density_5deg')):.4e}; "
            f"kappa={safe_float(r.get('kappa_R2_escape_proxy')):.4g}; "
            f"status={r.get('candidate_status')}; type={r.get('phb_type')}; m={r.get('m')}; R2={r.get('R2_abs')}"
        )
    lines.append("")
    lines.append("TOP by 5-degree main-lobe area density:")
    for r in by_density:
        lines.append(
            f"  density5={safe_float(r.get('aperture_main_lobe_area_density_5deg')):.4e}; "
            f"theta95={safe_float(r.get('aperture_theta95_deg')):.6g} deg; "
            f"lobe5={safe_float(r.get('aperture_lobe_fraction_5deg')):.4g}; "
            f"status={r.get('candidate_status')}; type={r.get('phb_type')}; m={r.get('m')}; R2={r.get('R2_abs')}"
        )
    lines.append("")
    lines.append("Interpretation:")
    if good:
        lines.append("  Есть кандидаты для confirm/full: гипотеза не опровергнута на первом быстром этапе.")
    elif checks:
        lines.append("  Есть CHECK-аномалии: стоит пересчитать лучшие случаи с большим resolution/a-over-lambda.")
    else:
        lines.append("  Узкий направленный выход из R2-щели на этой грубой сетке не найден.")
    lines.append("")
    lines.append("Important: this is not a final laser-resonator proof. It is only first-stage aperture-angle screening.")
    lines.append("v46.1 note: publication_exit_cone.png draws a symmetric cone around the horizontal focal-ring guides F−→F+ / +z at r=±R; theta is a half-angle from +z, not a vertical/r angle.")
    lines.append("v46.1 note: optional far_M2_z is excluded from kappa to avoid double-counting a diagnostic downstream slice.")
    lines.append("v46 note: default window_model=vertical_shift uses a right rectangular aperture in z≈+a wall with center shift relative to r=R.")
    lines.append("v46 note: kappa uses the MEEP aperture flux monitor when available; mp.get_fluxes(obj) is intentionally kept as standard MEEP API.")
    (outroot / "PHB_V46_UNIVERSAL_R2_APERTURE_REPORT_RU.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def estimate_runtime_minutes(args) -> Tuple[float, float, float]:
    """Very rough wall-clock estimate for the user's local Docker/MEEP setup.

    This is only an engineering planning note written before the run.  Actual
    time depends strongly on CPU, Docker/WSL load, geometry, PML and MEEP internals.
    """
    phb_types = parse_phb_types(args.phb_types)
    R2s = parse_float_list(args.R2_list)
    ms = parse_int_list(args.m_list)
    ncases = max(1, len(phb_types) * len(R2s) * len(ms))
    base_per_case_min = 5.0
    scale = (max(args.resolution, 1) / 64.0) ** 2
    scale *= max(args.after_sources, 1.0) / 120.0
    scale *= max(args.a_over_lambda, 0.25) / 1.0
    scale *= (max(args.nsrc, 1) / 12.0) ** 0.25
    if args.enable_far_monitor:
        scale *= 1.25
    if args.stage == "confirm":
        scale *= 1.35
    elif args.stage == "full":
        scale *= 2.0
    elif args.stage in ("plan", "geometry", "summary"):
        scale = 0.05
    center = ncases * base_per_case_min * scale
    return max(0.1, 0.45 * center), max(0.1, center), max(0.1, 2.5 * center)


def write_runtime_estimate(outroot: Path, args) -> Path:
    low, mid, high = estimate_runtime_minutes(args)
    phb_types = parse_phb_types(args.phb_types)
    R2s = parse_float_list(args.R2_list)
    ms = parse_int_list(args.m_list)
    ncases = len(phb_types) * len(R2s) * len(ms)
    lines = [
        "PHB v46 rough runtime estimate / примерная оценка времени расчёта",
        "===============================================================",
        "This file is generated before MEEP/FDTD execution and is not a guarantee.",
        "Эта оценка создаётся перед расчётом и не является гарантией времени.",
        "",
        f"stage: {args.stage}",
        f"number of cases: {ncases}",
        f"resolution: {args.resolution}",
        f"a/lambda: {args.a_over_lambda}",
        f"after_sources: {args.after_sources}",
        f"nsrc: {args.nsrc}",
        f"M2 enabled: {args.enable_far_monitor}; far_distance: {args.far_distance}",
        "",
        f"rough range: {low:.1f} - {high:.1f} minutes",
        f"central estimate: {mid:.1f} minutes",
        "",
        "Window convention:",
        "  W = R2 is absolute full aperture/monitor width, not R*R2.",
        "  shift = window_offset/window_shift is the radial CENTER shift from the focal-ring axis r=R.",
        "  interval = [R+shift-W/2, R+shift+W/2].",
        "  positive shift moves outward/upward; negative shift moves inward/downward toward r=0.",
    ]
    path = outroot / "RUNTIME_ESTIMATE_V46_RU.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_plan(args) -> None:
    outroot = ensure_dir(Path(args.outroot))
    g0 = PHBGeometry(a=args.a, b=args.b, R=args.R, R2=parse_float_list(args.R2_list)[0])
    wave = make_wave_report(g0, args)
    phb_types = parse_phb_types(args.phb_types)
    R2s = parse_float_list(args.R2_list)
    ms = parse_int_list(args.m_list)
    lines = []
    lines.append("PHB v46 first aperture verification plan")
    lines.append(f"Created: {now_iso()}")
    lines.append("")
    lines.append("Scientific target:")
    lines.append("  Verify/falsify first-stage PHB hypothesis by local Poynting angle of energy exiting the selected shifted output window M1.")
    lines.append("  v46 convention: M1 and optional M2 lie on one longitudinal line parallel to F−→F+ / +z, with common radius r=R+shift.")
    lines.append("  The script intentionally separates local aperture directionality from later far-field and resonator-level claims.")
    lines.append("")
    lines.append("Geometry:")
    lines.append(f"  a={args.a}, b={args.b}, R={args.R}")
    lines.append(f"  R2 output-window width list W={R2s}; R2 is the absolute full monitor/window width")
    lines.append(f"  window shift={args.window_offset}; interval = [R+shift-W/2, R+shift+W/2]")
    lines.append(f"  PHB types={phb_types}")
    lines.append(f"  window_model={args.window_model}; default v46 vertical_shift means right rectangular aperture centered at r=R+shift in z≈+a wall")
    if "halfring" in phb_types:
        lines.append("  Note: halfring means semicircular meridional rounded/toroid-like closure; it is a first-stage proxy, not a strict analytical torus.")
    lines.append(f"  m list={ms}")
    lines.append("")
    lines.append("Wave scale:")
    for k, v in wave.items():
        lines.append(f"  {k}: {v}")
    lines.append("")
    lines.append("Main monitors:")
    lines.append("  M1 aperture plane immediately after the corrected right output window")
    lines.append("  radial interval for M1 and default M2 is [R+shift-W/2, R+shift+W/2], center r=R+shift")
    lines.append("  H fields are always included; Poynting vector is required for the angle.")
    lines.append(f"  M2 far monitor enabled={args.enable_far_monitor}; distance from M1={args.far_distance}; mode={args.far_monitor_mode}")
    lines.append("  M2 is parallel to M1 and by default centered on the same shifted radial line r=R+shift.")
    lines.append("")
    lines.append("Main outputs:")
    lines.append("  cases/*/geometry_monitor.png")
    lines.append("  cases/*/aperture_profile.csv")
    lines.append("  cases/*/aperture_profile.png")
    lines.append("  cases/*/publication_exit_cone.png  # correct symmetric cone around focal-ring guides F−→F+ / +z")
    lines.append("  cases/*/publication_angular_density.png  # D(theta) energy-density plot")
    lines.append("  cases/*/far_M2_profile.csv, if --enable-far-monitor")
    lines.append("  cases/*/far_M2_publication_angular_density.png, if --enable-far-monitor")
    lines.append("  cases/*/far_M2_profile.png, if --enable-far-monitor")
    lines.append("  cases/*/metrics.json")
    lines.append("  PHB_V46_UNIVERSAL_R2_APERTURE_SUMMARY.csv")
    lines.append("  PHB_V46_UNIVERSAL_R2_APERTURE_REPORT_RU.txt")
    lines.append("  logs/runs.jsonl")
    lines.append("  ZIP archive beside the result folder, if --archive is on")
    (outroot / "RUN_PLAN_V46_RU.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    est_path = write_runtime_estimate(outroot, args)
    print(outroot / "RUN_PLAN_V46_RU.txt")
    print(est_path)

def run_grid(args) -> int:
    outroot = ensure_dir(Path(args.outroot))
    ensure_dir(outroot / "logs")
    write_runtime_estimate(outroot, args)
    phb_types = parse_phb_types(args.phb_types)
    R2s = parse_float_list(args.R2_list)
    ms = parse_int_list(args.m_list)
    if not R2s:
        raise ValueError("--R2-list is empty")
    if not ms:
        raise ValueError("--m-list is empty")

    # One geometry plot per type/R2 first, even before MEEP.
    for phb_type in phb_types:
        for R2 in R2s:
            g = PHBGeometry(a=args.a, b=args.b, R=args.R, R2=R2, window_offset=args.window_offset, window_offset_mode=args.window_offset_mode)
            cfg = RunConfig(
                stage=args.stage, outroot=args.outroot, phb_type=phb_type, m=ms[0], resolution=args.resolution,
                a_over_lambda=args.a_over_lambda, min_a_over_lambda=args.min_a_over_lambda,
                dpml_over_lambda=args.dpml_over_lambda, dpml=args.dpml, wall_thickness=args.wall_thickness,
                window_model=args.window_model, window_offset=args.window_offset, window_offset_mode=args.window_offset_mode, vertical_window_thickness=args.vertical_window_thickness,
                source_components=args.source_components, source_mode=args.source_mode, nsrc=args.nsrc,
                seed=args.seed, fwidth_frac=args.fwidth_frac, after_sources=args.after_sources,
                skip_existing=args.skip_existing, allow_underresolved=args.allow_underresolved,
                min_radial_slot_cells=args.min_radial_slot_cells, min_axial_slot_cells=args.min_axial_slot_cells,
                aperture_offset_cells=args.aperture_offset_cells,
                enable_far_monitor=args.enable_far_monitor, far_distance=args.far_distance,
                far_capture_angle_deg=args.far_capture_angle_deg, far_monitor_mode=args.far_monitor_mode, far_zone_safety=args.far_zone_safety,
                allow_near_field_far_monitor=args.allow_near_field_far_monitor,
                narrow_theta95_deg=args.narrow_theta95_deg, narrow_lobe5_min=args.narrow_lobe5_min,
                useful_flux_min=args.useful_flux_min, kappa_min=args.kappa_min, min_free_gb=args.min_free_gb,
                archive=args.archive, stop_on_error=args.stop_on_error,
                save_field_maps=args.save_field_maps,
                field_map_regions=args.field_map_regions,
                field_map_components=args.field_map_components,
                field_map_r_pad=args.field_map_r_pad,
                field_map_z_half=args.field_map_z_half,
                field_map_corridor_half_width=args.field_map_corridor_half_width,
                field_map_interior_rmax=args.field_map_interior_rmax,
                field_map_max_points_per_region=args.field_map_max_points_per_region,
                field_map_csv_stride=args.field_map_csv_stride,
                field_map_png=args.field_map_png,
                field_map_npz=args.field_map_npz,
                save_material_maps=args.save_material_maps,
                save_source_excluded_maps=args.save_source_excluded_maps,
                source_exclusion_cells=args.source_exclusion_cells,
                m2_extra_flux_warning_ratio=args.m2_extra_flux_warning_ratio,
            )
            plot_geometry(outroot / "geometry" / f"geometry_{phb_type}_R2_{str(R2).replace('.', 'p')}.png", g, cfg)

    rows: List[Dict[str, object]] = []
    if args.stage == "geometry":
        rows = collect_existing_metrics(outroot)
        write_summary(outroot, rows, args)
        if args.archive:
            archive_folder(outroot)
        return 0

    for phb_type in phb_types:
        for m in ms:
            for R2 in R2s:
                print(f"\n=== Running PHB v46: type={phb_type}, m={m}, R2={R2}, res={args.resolution}, a/lambda={args.a_over_lambda} ===", flush=True)
                row = run_one_case(args, phb_type, m, R2)
                rows.append(row)
                write_summary(outroot, collect_existing_metrics(outroot), args)

    all_rows = collect_existing_metrics(outroot)
    write_summary(outroot, all_rows, args)
    if args.archive:
        z = archive_folder(outroot)
        if z:
            print(f"Archive written: {z}")
    bad = [r for r in all_rows if str(r.get("candidate_status", "")).startswith("ERROR")]
    return 1 if bad else 0


# -----------------------------
# CLI
# -----------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--stage", choices=["plan", "geometry", "fast", "confirm", "full", "summary"], default="plan")
    ap.add_argument("--outroot", default="/work/PHB_V46_UNIVERSAL_SHIFTED_WINDOW")
    ap.add_argument("--phb-types", default="open,cylinder,halfring")
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--b", type=float, default=1.0)
    ap.add_argument("--R", type=float, default=1.0)
    ap.add_argument("--R2-list", default="0.05,0.07,0.10,0.20,0.30", help="absolute output-window widths W; R2 is the window size, not R*R2")
    ap.add_argument("--m-list", default="0,1,2")
    ap.add_argument("--a-over-lambda", type=float, default=1.0, help="primary scale: a/lambda; must be >= min")
    ap.add_argument("--min-a-over-lambda", type=float, default=1.0)
    ap.add_argument("--resolution", type=int, default=64)
    ap.add_argument("--dpml", type=float, default=None)
    ap.add_argument("--dpml-over-lambda", type=float, default=1.2)
    ap.add_argument("--wall-thickness", type=float, default=0.035)
    ap.add_argument("--window-model", choices=["vertical_shift", "vertical_rect", "vertical_upper", "upper_rect", "offset_upper", "horn_cut"], default="vertical_shift",
                    help="v46 default: right rectangular window [R+shift-W/2,R+shift+W/2] in vertical z≈+a wall; horn_cut keeps legacy hyperbolic-cut behavior")
    ap.add_argument("--window-offset", "--window-shift", dest="window_offset", type=float, default=0.0,
                    help="radial shift from focal-axis radius r=R; interpretation is selected by --window-offset-mode")
    ap.add_argument("--window-offset-mode", choices=["center", "lower_edge"], default="center",
                    help="center: window=[R+offset-W/2,R+offset+W/2]; lower_edge: window=[R+offset,R+offset+W], compatible with v45.5 scans")
    ap.add_argument("--vertical-window-thickness", type=float, default=None,
                    help="physical z-thickness of the vertical rectangular output window; default=max(wall_thickness,min_axial_slot_cells/resolution)")
    ap.add_argument("--source-components", default="Ez,Er,Ep")
    ap.add_argument("--source-mode", choices=["random", "coherent", "single"], default="random")
    ap.add_argument("--nsrc", type=int, default=12)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--fwidth-frac", type=float, default=0.18)
    ap.add_argument("--after-sources", type=float, default=120.0, help="fixed time after source turn-off; avoids indefinite run")
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--allow-underresolved", action="store_true")
    ap.add_argument("--min-radial-slot-cells", type=float, default=4.0)
    ap.add_argument("--min-axial-slot-cells", type=float, default=4.0)
    ap.add_argument("--aperture-offset-cells", type=float, default=2.0, help="place monitor this many grid cells after effective R2 horn segment")
    ap.add_argument("--enable-far-monitor", action="store_true", help="optional M2 far-zone monitor; works for open,cylinder,halfring")
    ap.add_argument("--far-distance", type=float, default=20.0)
    ap.add_argument("--far-capture-angle-deg", type=float, default=20.0, help="angular half-capture used to size the optional M2 far monitor")
    ap.add_argument("--far-monitor-mode", choices=["focal_axis", "same_window"], default="focal_axis", help="M2 radial reference: M2 follows shifted M1 center in v46; option retained for compatibility")
    ap.add_argument("--far-zone-safety", type=float, default=1.0)
    ap.add_argument("--allow-near-field-far-monitor", action="store_true")
    ap.add_argument("--narrow-theta95-deg", type=float, default=10.0)
    ap.add_argument("--narrow-lobe5-min", type=float, default=0.50)
    ap.add_argument("--useful-flux-min", type=float, default=1e-12)
    ap.add_argument("--kappa-min", type=float, default=0.02)
    ap.add_argument("--min-free-gb", type=float, default=2.0)

    # v47/v46.1 diagnostic field-map output.  Keep disabled by default for fast scans;
    # enable for confirm/forensic runs around a promising candidate.
    ap.add_argument("--save-field-maps", action="store_true",
                    help="save 2D complex DFT maps and derived Poynting/phase maps inside PHB and near the slot")
    ap.add_argument("--field-map-regions", default="near_slot,focal_corridor",
                    help="comma list: none,all,interior,near_slot,pre_slot,post_slot,focal_corridor")
    ap.add_argument("--field-map-components", default="Er,Ep,Ez,Hr,Hp,Hz",
                    help="DFT components to save; all six are needed for Poynting maps")
    ap.add_argument("--field-map-r-pad", type=float, default=0.06,
                    help="radial padding around R2 window for near_slot/pre_slot/post_slot maps")
    ap.add_argument("--field-map-z-half", type=float, default=0.25,
                    help="half-size in z for near-slot field maps around z≈+a")
    ap.add_argument("--field-map-corridor-half-width", type=float, default=0.08,
                    help="half-width of longitudinal field corridor around r=R+shift")
    ap.add_argument("--field-map-interior-rmax", type=float, default=None,
                    help="optional r_max for the interior full PHB field map; default follows wall span")
    ap.add_argument("--field-map-max-points-per-region", type=int, default=800000,
                    help="skip a requested 2D field-map region if estimated grid points exceed this limit")
    ap.add_argument("--field-map-csv-stride", type=int, default=8,
                    help="write downsampled CSV field map with this stride; 0 disables CSV sampling")
    ap.add_argument("--no-field-map-png", dest="field_map_png", action="store_false")
    ap.add_argument("--no-field-map-npz", dest="field_map_npz", action="store_false")
    ap.set_defaults(field_map_png=True, field_map_npz=True)

    # v47.2 additional verification output. These are on by default because they
    # are cheap compared with the MEEP run and are critical for artefact control.
    ap.add_argument("--no-material-maps", dest="save_material_maps", action="store_false",
                    help="disable diagnostic material/metal/window mask maps in field_map regions")
    ap.add_argument("--no-source-excluded-maps", dest="save_source_excluded_maps", action="store_false",
                    help="disable |E|² maps with source-near-field cells excluded")
    ap.add_argument("--source-exclusion-cells", type=float, default=3.0,
                    help="radius around point sources, in grid cells, masked in source-excluded |E|² maps")
    ap.add_argument("--m2-extra-flux-warning-ratio", type=float, default=1.25,
                    help="flag M2 as containing extra field if M2/M1 flux ratio exceeds this value")
    ap.set_defaults(save_material_maps=True, save_source_excluded_maps=True)

    ap.add_argument("--no-archive", dest="archive", action="store_false")
    ap.set_defaults(archive=True)
    ap.add_argument("--stop-on-error", action="store_true")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "window_model", "vertical_shift") in ("upper_rect", "offset_upper", "vertical_upper", "vertical_rect"):
        args.window_model = "vertical_shift"
    outroot = ensure_dir(Path(args.outroot))
    if args.stage == "plan":
        write_plan(args)
        if args.archive:
            z = archive_folder(outroot)
            if z:
                print(f"Archive written: {z}")
        return 0
    if args.stage == "summary":
        rows = collect_existing_metrics(outroot)
        write_summary(outroot, rows, args)
        if args.archive:
            z = archive_folder(outroot)
            if z:
                print(f"Archive written: {z}")
        return 0
    return run_grid(args)


if __name__ == "__main__":
    raise SystemExit(main())

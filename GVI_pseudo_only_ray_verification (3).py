#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GVI PSEUDOHYPERBOLOID-ONLY RAY-BILLIARD VERIFICATION SCRIPT
==================================================================

Author of the research programme: Vladimir I. Khaustov
Topic: Geometric Wave Engineering / second-order vertical and horizontal pseudohyperboloids

Что делает скрипт
-----------------
Единый скрипт считает лучево-бильярдную проверку для двух родственных объектов:

1) vertical pseudohyperboloid, ось вращения X:
   - фокальная зона = тонкое экваториальное кольцо радиуса R;
   - условие зоны: |x| <= a и |rho - R| <= h_abs;
   - rho = sqrt(y^2 + z^2).

2) horizontal pseudohyperboloid, ось вращения Y:
   - фокальная зона = две тонкие дисковые зоны около y=+R и y=-R;
   - радиус каждого диска = a;
   - условие зоны: ||y|-R| <= h_abs и rho <= a;
   - rho = sqrt(x^2 + z^2).

Главное исправление по замечаниям рецензента
--------------------------------------------
В старой версии h задавался как абсолютная величина. Здесь h задаётся одинаково
для обоих типов как безразмерная величина h/a:

    H_OVER_A_VALUES = [0.05, 0.10, 0.20, 0.30]
    h_abs = (h/a) * a

То есть для a=0.1 и h/a=0.20 фактическая полуширина фокальной зоны h_abs=0.02.

Также исправлена/добавлена метрика strict-плотности:

1) legacy metric:
   strict_sequence_density_per_100_bounces_pct
   = 100 * среднее число strict-контактов / n_bounces.
   Это старая метрика, сохранена только для сопоставления.

2) corrected metric:
   strict_sequence_density_per_actual_reflection_pct
   = 100 * сумма strict-контактов / сумма фактических отражающих контактов.
   Именно её следует использовать как основную научную метрику.

3) open-режимы:
   для open-геометрий добавлена отдельная метрика focal_exit_contacts_per_100_bounces_pct,
   потому что контакт с открытой фокальной апертурой является выходом, а не отражающим streak.

Ограничение модели
------------------
Это НЕ полноволновой расчёт. В модели нет длины волны, фазы, интерференции,
дифракции, добротности, мод Гельмгольца/Максвелла. Результат относится только
к геометрической маршрутизации лучей и зеркальному отражению.

Новые входные условия источников
--------------------------------
1) center_point_100:
   100 лучей из одного центрального точечного источника, равномерно по сфере;
   каждый луч по умолчанию отслеживается 1000 отражений/контактов.

2) surface_100x100:
   100 точек, равномерно по площади внутренней поверхности вращения;
   из каждой точки запускается 100 лучей равномерно по внутренней полусфере направлений;
   всего 10 000 лучей, каждый луч по умолчанию отслеживается 1000 отражений/контактов.

Как запускать
-------------
1) Обычный запуск:
       python GVI_universal_ray_verification_corrected.py

2) Быстрый тест, чтобы проверить, что всё работает:
       python GVI_universal_ray_verification_corrected.py --fast-test

3) Изменять входные условия удобно в блоке USER INPUT BLOCK ниже.

Выходные файлы
--------------
В папке OUTDIR создаются:
- all_results.csv                         полный результат;
- compact_pseudo_J_only.csv               компактная таблица только pseudo и J=low/mid/high;
- summary_by_object_j.csv                 сводка по типу объекта и J-классам;
- top_cases_corrected_metric.csv          топ случаев по исправленной метрике;
- run_metadata.json                       параметры запуска;
- figures/*.png                           базовые графики.
"""

# =============================================================================
# USER INPUT BLOCK — ОСНОВНЫЕ НАСТРОЙКИ, КОТОРЫЕ МОЖНО МЕНЯТЬ
# =============================================================================

OUTDIR = "gvi_pseudo_only_ray_output"

# Какие объекты считать: ["vertical"], ["horizontal"] или оба.
OBJECT_TYPES = ["vertical", "horizontal"]

# Только псевдогиперболоидные поверхности. Контрольные фигуры полностью исключены из расчётного sweep.
BODY_MODELS = ["pseudo"]

# Границы: open — открытая фокальная/горловинная область; closed — отражающее замыкание.
GEOMETRY_TYPES = ["open", "closed"]

# Источники:
# center_point_100 — лучи из центра;
# surface_100x100 — точки на поверхности и направления внутрь.
SOURCES = ["center_point_100", "surface_100x100"]

# h задаётся как h/a одинаково для vertical и horizontal.
# В коде h_abs = h_over_a * a.
H_OVER_A_VALUES = [0.05, 0.10, 0.20, 0.30]

# h_abs трактуется как полуширина зоны:
# vertical: |rho - R| <= h_abs;
# horizontal: ||y|-R| <= h_abs.
H_IS_HALF_WIDTH = True

# Число отражений. Для сходимости можно поставить [100, 500, 1000].
# При полном 100x100 и всех моделях 500/1000 может считаться долго.
N_BOUNCES_LIST = [1000]

# Число отрезков для поиска пересечения. Для сходимости можно поставить [64, 128].
N_SCAN_LIST = [64]

# Monte Carlo оценка объёмов тела и фокальной зоны.
# 0 отключает MC и оставляет только быстрые аналитические proxy-объёмы.
MONTE_CARLO_VOLUME_SAMPLES = 20000
MONTE_CARLO_VOLUME_SEED = 12345

# CI95: лёгкая доверительная оценка по лучам внутри каждого seed + агрегирование по seed.
# BOOTSTRAP_RESAMPLES=0 оставляет быструю normal-approx CI по ray-level значениям.
BOOTSTRAP_RESAMPLES = 0

# Дополнительные сервисные диагностики. По умолчанию выключены, чтобы полный расчёт не растягивался.
RUN_DIAGNOSTICS = False
DIAGNOSTIC_TOP_CASES = 3
DIAGNOSTIC_RAYS_PER_CASE = 5
DIAGNOSTIC_BOUNCES = 300
POINCARE_RAYS_PER_CASE = 3

# Размеры ансамблей.
CENTER_RAYS = 100
SURFACE_POINTS = 100
DIRECTIONS_PER_POINT = 100

# Источник направлений:
# "quasi"  — детерминированный golden-angle, как в старой версии;
# "random" — случайный с фиксированными seed, лучше для статистической проверки.
SAMPLING_MODE = "random"
RANDOM_SEEDS = [1]       # для seed-статистики можно [1,2,3,4,5,6,7,8,9,10]

# Пороговые J-классы. Это эвристическая, но физически полезная сегментация.
LOW_J_THRESHOLD = 0.2
MID_J_THRESHOLD = 0.6

# Строгий критерий: streak должен быть не короче этого числа.
MIN_STRICT_STREAK = 2

# 10 геометрий, использованных в прежних отчётах.
BEST10 = [
    dict(geom_index=1,  previous_label="G1_vertical",              R=5.0,  a=0.10, b=0.5),
    dict(geom_index=2,  previous_label="G2_helical_whisper",       R=5.0,  a=0.50, b=0.5),
    dict(geom_index=3,  previous_label="G3_sticky",                R=50.0, a=0.50, b=0.5),
    dict(geom_index=4,  previous_label="G4_ldos_label_only",       R=5.0,  a=0.10, b=1.0),
    dict(geom_index=5,  previous_label="G5_horn_axial_shuttle",    R=10.0, a=0.50, b=1.0),
    dict(geom_index=6,  previous_label="G6_tip_trap",              R=5.0,  a=0.05, b=1.0),
    dict(geom_index=7,  previous_label="G7_cap_cycle",             R=15.0, a=0.50, b=1.0),
    dict(geom_index=8,  previous_label="G8_periodic",              R=50.0, a=0.05, b=1.0),
    dict(geom_index=9,  previous_label="G9_directed_exit_label",   R=45.0, a=0.05, b=0.5),
    dict(geom_index=10, previous_label="G10_open_retention",       R=25.0, a=0.10, b=0.1),
]

# =============================================================================
# IMPLEMENTATION
# =============================================================================

import argparse
import json
import math
import time
import sys
import platform
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from numba import njit, prange, get_num_threads
    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False
    def njit(*args, **kwargs):
        def deco(f):
            return f
        return deco
    prange = range
    def get_num_threads():
        return 1

PI = math.pi
GOLD = 2.39996322972865332

OBJECT_ID = {"vertical": 0, "horizontal": 1}
OBJECT_NAME = {0: "vertical", 1: "horizontal"}
BODY_ID = {"pseudo": 0, "null_cylinder": 1, "null_ellipsoid": 2, "null_sphere": 3}
BODY_NAME = {0: "pseudo", 1: "null_cylinder", 2: "null_ellipsoid", 3: "null_sphere"}
GEOM_ID = {"open": 0, "closed": 1}
GEOM_NAME = {0: "open", 1: "closed"}
SOURCE_ID = {"center_point_100": 0, "surface_100x100": 1}
SOURCE_NAME = {0: "center_point_100", 1: "surface_100x100"}
J_LABELS = {0: "low-J", 1: "mid-J", 2: "high-J"}


@njit(cache=False)
def L_vertical(a, b, R):
    return a * math.sqrt(1.0 + (R / b) * (R / b))


@njit(cache=False)
def horizontal_center_radius(a, b, R):
    return a * math.sqrt(1.0 + (R / b) * (R / b))


@njit(cache=False)
def axis_extent_nb(object_id, a, b, R, geom_type, body_model):
    if body_model == 3:  # null sphere baseline, radius R, same selected axis
        return R
    if object_id == 0:
        return L_vertical(a, b, R)
    # horizontal
    if body_model == 0:
        return R + (a if geom_type == 1 else 0.0)
    return R + (a if geom_type == 1 else 0.0)


@njit(cache=False)
def reference_radius_nb(object_id, a, b, R, geom_type, body_model):
    if body_model == 3:  # null sphere, radius R
        return max(R, 1e-12)
    if object_id == 0:
        # vertical: main ring radius R; closed cap may extend to R+a, but J is normalized by R.
        return max(R, 1e-12)
    # horizontal: normalize by the maximum central radius.
    return max(horizontal_center_radius(a, b, R), 1e-12)


@njit(cache=False)
def rmax_axis_nb(object_id, s, a, b, R, geom_type, body_model):
    """Maximum transverse radius at axial coordinate s.
    object_id=0: s=x, radius=sqrt(y^2+z^2)
    object_id=1: s=y, radius=sqrt(x^2+z^2)
    """
    ax = abs(s)
    if object_id == 0:
        L = L_vertical(a, b, R)
        if body_model == 0:  # vertical pseudo
            if ax <= a:
                if geom_type == 1:
                    q = a * a - s * s
                    if q < 0.0:
                        q = 0.0
                    return R + math.sqrt(q)
                return R
            if ax <= L:
                q = (ax / a) * (ax / a) - 1.0
                if q < 0.0:
                    q = 0.0
                rr = R - b * math.sqrt(q)
                if rr < 0.0:
                    rr = 0.0
                return rr
            return -1.0
        if body_model == 1:  # null cylinder
            return R if ax <= L else -1.0
        if body_model == 3:  # null sphere, radius R
            if ax <= R:
                q = 1.0 - (s / R) * (s / R)
                if q < 0.0:
                    return -1.0
                return R * math.sqrt(q)
            return -1.0
        # null ellipsoid
        if ax <= L:
            q = 1.0 - (s / L) * (s / L)
            if q < 0.0:
                return -1.0
            return R * math.sqrt(q)
        return -1.0

    # horizontal object
    E = R + (a if geom_type == 1 else 0.0)
    r0 = horizontal_center_radius(a, b, R)
    if body_model == 0:  # horizontal pseudo
        if ax <= R:
            q = (R - ax) / b
            return a * math.sqrt(1.0 + q * q)
        if geom_type == 1 and ax <= R + a:
            q = a * a - (ax - R) * (ax - R)
            if q < 0.0:
                return -1.0
            return math.sqrt(q)
        return -1.0
    if body_model == 1:  # null cylinder, same axial extent and central radius as pseudo reference
        return r0 if ax <= E else -1.0
    if body_model == 3:  # null sphere, radius R
        if ax <= R:
            q = 1.0 - (s / R) * (s / R)
            if q < 0.0:
                return -1.0
            return R * math.sqrt(q)
        return -1.0
    # null ellipsoid
    if ax <= E:
        q = 1.0 - (s / E) * (s / E)
        if q < 0.0:
            return -1.0
        return r0 * math.sqrt(q)
    return -1.0


@njit(cache=False)
def drds_axis_nb(object_id, s, a, b, R, geom_type, body_model):
    ax = abs(s)
    sign = 1.0 if s >= 0.0 else -1.0
    if object_id == 0:
        L = L_vertical(a, b, R)
        if body_model == 0:
            if ax <= a:
                if geom_type == 1:
                    q = a * a - s * s
                    if q < 1e-18:
                        q = 1e-18
                    return -s / math.sqrt(q)
                return 0.0
            q = (s / a) * (s / a) - 1.0
            if q < 1e-18:
                q = 1e-18
            return -b * s / (a * a * math.sqrt(q))
        if body_model == 1:
            return 0.0
        if body_model == 3:
            q_s = R * R - s * s
            if q_s < 1e-18:
                q_s = 1e-18
            return -s / math.sqrt(q_s)
        q = 1.0 - (s / L) * (s / L)
        if q < 1e-18:
            q = 1e-18
        return R * (-s / (L * L)) / math.sqrt(q)

    # horizontal
    E = R + (a if geom_type == 1 else 0.0)
    r0 = horizontal_center_radius(a, b, R)
    if body_model == 0:
        if ax <= R:
            q = (R - ax) / b
            den = math.sqrt(1.0 + q * q)
            if den < 1e-18:
                return 0.0
            return -a * q * sign / (b * den)
        if geom_type == 1 and ax <= R + a:
            q = a * a - (ax - R) * (ax - R)
            if q < 1e-18:
                q = 1e-18
            return -(ax - R) * sign / math.sqrt(q)
        return 0.0
    if body_model == 1:
        return 0.0
    if body_model == 3:
        q_s = R * R - s * s
        if q_s < 1e-18:
            q_s = 1e-18
        return -s / math.sqrt(q_s)
    q = 1.0 - (s / E) * (s / E)
    if q < 1e-18:
        q = 1e-18
    return r0 * (-s / (E * E)) / math.sqrt(q)


@njit(cache=False)
def inside_nb(x, y, z, object_id, a, b, R, geom_type, body_model):
    if object_id == 0:
        s = x
        rad2 = y * y + z * z
    else:
        s = y
        rad2 = x * x + z * z
    rr = rmax_axis_nb(object_id, s, a, b, R, geom_type, body_model)
    if rr < 0.0:
        return -1.0
    return rr * rr - rad2


@njit(cache=False)
def normal_nb(x, y, z, object_id, a, b, R, geom_type, body_model):
    """Outward unit normal for surface of revolution.
    Detects cylinder end caps separately.
    """
    if object_id == 0:
        s = x
        rho = math.sqrt(y * y + z * z)
        E = axis_extent_nb(object_id, a, b, R, geom_type, body_model)
        # cylinder end caps
        if body_model == 1 and abs(abs(s) - E) < 1e-6 * max(E, 1.0):
            return (1.0 if s >= 0.0 else -1.0), 0.0, 0.0
        if rho < 1e-18:
            return (1.0 if s >= 0.0 else -1.0), 0.0, 0.0
        rp = drds_axis_nb(object_id, s, a, b, R, geom_type, body_model)
        nx = -rp
        ny = y / rho
        nz = z / rho
        nn = math.sqrt(nx * nx + ny * ny + nz * nz)
        return nx / nn, ny / nn, nz / nn

    # horizontal: axial coordinate is y
    s = y
    rho = math.sqrt(x * x + z * z)
    E = axis_extent_nb(object_id, a, b, R, geom_type, body_model)
    if body_model == 1 and abs(abs(s) - E) < 1e-6 * max(E, 1.0):
        return 0.0, (1.0 if s >= 0.0 else -1.0), 0.0
    if rho < 1e-18:
        return 0.0, (1.0 if s >= 0.0 else -1.0), 0.0
    rp = drds_axis_nb(object_id, s, a, b, R, geom_type, body_model)
    nx = x / rho
    ny = -rp
    nz = z / rho
    nn = math.sqrt(nx * nx + ny * ny + nz * nz)
    return nx / nn, ny / nn, nz / nn


@njit(cache=False)
def reflect_nb(vx, vy, vz, nx, ny, nz):
    dot = vx * nx + vy * ny + vz * nz
    vx2 = vx - 2.0 * dot * nx
    vy2 = vy - 2.0 * dot * ny
    vz2 = vz - 2.0 * dot * nz
    nn = math.sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2)
    if nn > 0.0:
        vx2 /= nn
        vy2 /= nn
        vz2 /= nn
    return vx2, vy2, vz2


@njit(cache=False)
def next_hit_nb(x0, y0, z0, vx, vy, vz, object_id, a, b, R, geom_type, body_model, n_scan):
    E = axis_extent_nb(object_id, a, b, R, geom_type, body_model)
    rr0 = reference_radius_nb(object_id, a, b, R, geom_type, body_model)
    scale = max(max(E, rr0 + a), 1.0)
    tmax = 4.0 * math.sqrt(E * E + (rr0 + a) * (rr0 + a)) + 2.0 * scale
    if tmax < 1.0:
        tmax = 1.0
    tprev = 1e-8 * scale
    fprev = inside_nb(x0 + vx * tprev, y0 + vy * tprev, z0 + vz * tprev,
                      object_id, a, b, R, geom_type, body_model)
    dt = tmax / n_scan
    for k in range(n_scan):
        t = (k + 1) * dt
        f = inside_nb(x0 + vx * t, y0 + vy * t, z0 + vz * t,
                      object_id, a, b, R, geom_type, body_model)
        if f <= 0.0 and fprev > 0.0:
            lo = tprev
            hi = t
            for _ in range(42):
                mid = 0.5 * (lo + hi)
                fm = inside_nb(x0 + vx * mid, y0 + vy * mid, z0 + vz * mid,
                               object_id, a, b, R, geom_type, body_model)
                if fm > 0.0:
                    lo = mid
                else:
                    hi = mid
            th = 0.5 * (lo + hi)
            return True, th, x0 + vx * th, y0 + vy * th, z0 + vz * th
        tprev = t
        fprev = f
    return False, 1e300, 0.0, 0.0, 0.0


@njit(cache=False)
def angular_momentum_axis_nb(x, y, z, vx, vy, vz, object_id):
    if object_id == 0:
        return y * vz - z * vy  # Lx
    return z * vx - x * vz      # Ly


@njit(cache=False)
def classify_j_nb(Jhat, low_thr, mid_thr):
    if Jhat < low_thr:
        return 0.0
    if Jhat < mid_thr:
        return 1.0
    return 2.0


@njit(cache=False)
def is_open_exit_nb(xh, yh, zh, object_id, a, b, R, geom_type, body_model):
    if geom_type != 0:
        return False
    if object_id == 0:
        # vertical: open central equatorial band, as in the original model
        return abs(xh) <= a + 1e-9
    # horizontal: open at two throat/end disks. For pseudo this is y=±R, radius a.
    rho = math.sqrt(xh * xh + zh * zh)
    if body_model == 0:
        scale = max(R, horizontal_center_radius(a, b, R), 1.0)
        return (abs(abs(yh) - R) <= max(1e-8, 1e-6 * scale)) and (rho <= a * (1.0 + 1e-6))
    # null controls: open at axial ends of the control body
    E = axis_extent_nb(object_id, a, b, R, geom_type, body_model)
    scale = max(E, reference_radius_nb(object_id, a, b, R, geom_type, body_model), 1.0)
    return abs(abs(yh) - E) <= max(1e-8, 1e-6 * scale)


@njit(cache=False)
def in_focal_zone_nb(xh, yh, zh, object_id, a, b, R, h_abs):
    if object_id == 0:
        rho = math.sqrt(yh * yh + zh * zh)
        return (abs(xh) <= a + 1e-12) and (abs(rho - R) <= h_abs)
    rho = math.sqrt(xh * xh + zh * zh)
    # Corrected horizontal definition: two disks of radius a near y=±R, not 0<=rho<=2a.
    return (abs(abs(yh) - R) <= h_abs) and (rho <= a + 1e-12)


@njit(parallel=True, cache=False)
def simulate_batch_nb(pos, vel, object_id, a, b, R, geom_type, body_model,
                      h_abs_values, n_bounces, n_scan, low_thr, mid_thr, min_streak):
    n_rays = pos.shape[0]
    nh = h_abs_values.shape[0]
    # columns: Jhat, Jclass, exit, nohit, bounces_done, L_drift_norm, speed_drift, actual_reflections, per-h blocks of 6
    out = np.zeros((n_rays, 8 + 6 * nh), dtype=np.float64)
    refR = reference_radius_nb(object_id, a, b, R, geom_type, body_model)
    E = axis_extent_nb(object_id, a, b, R, geom_type, body_model)
    scale = max(max(E, refR + a), 1.0)

    for i in prange(n_rays):
        x0 = pos[i, 0]
        y0 = pos[i, 1]
        z0 = pos[i, 2]
        vx = vel[i, 0]
        vy = vel[i, 1]
        vz = vel[i, 2]

        L0 = angular_momentum_axis_nb(x0, y0, z0, vx, vy, vz, object_id)
        Jhat = abs(L0) / max(refR, 1e-12)
        jclass = classify_j_nb(Jhat, low_thr, mid_thr)
        out[i, 0] = Jhat
        out[i, 1] = jclass

        exit_flag = 0.0
        nohit_flag = 0.0
        bounces_done = 0.0
        actual_reflections = 0.0

        old_refl = np.zeros(nh, dtype=np.float64)
        old_all = np.zeros(nh, dtype=np.float64)
        strict_seq = np.zeros(nh, dtype=np.float64)
        streak = np.zeros(nh, dtype=np.float64)
        max_streak = np.zeros(nh, dtype=np.float64)
        focal_exit = np.zeros(nh, dtype=np.float64)

        for ib in range(n_bounces):
            ok, t, xh, yh, zh = next_hit_nb(x0, y0, z0, vx, vy, vz, object_id, a, b, R,
                                             geom_type, body_model, n_scan)
            if not ok:
                for ih in range(nh):
                    if streak[ih] >= min_streak:
                        strict_seq[ih] += streak[ih]
                    streak[ih] = 0.0
                nohit_flag = 1.0
                bounces_done = ib
                break

            is_exit = is_open_exit_nb(xh, yh, zh, object_id, a, b, R, geom_type, body_model)
            for ih in range(nh):
                in_focal = in_focal_zone_nb(xh, yh, zh, object_id, a, b, R, h_abs_values[ih])
                if in_focal:
                    old_all[ih] += 1.0
                if in_focal and (not is_exit):
                    old_refl[ih] += 1.0
                    streak[ih] += 1.0
                    if streak[ih] > max_streak[ih]:
                        max_streak[ih] = streak[ih]
                else:
                    if streak[ih] >= min_streak:
                        strict_seq[ih] += streak[ih]
                    streak[ih] = 0.0
                    if in_focal and is_exit:
                        focal_exit[ih] += 1.0

            if is_exit:
                for ih in range(nh):
                    if streak[ih] >= min_streak:
                        strict_seq[ih] += streak[ih]
                    streak[ih] = 0.0
                exit_flag = 1.0
                bounces_done = ib + 1
                break

            nx, ny, nz = normal_nb(xh, yh, zh, object_id, a, b, R, geom_type, body_model)
            vx, vy, vz = reflect_nb(vx, vy, vz, nx, ny, nz)
            actual_reflections += 1.0
            tiny = 1e-8 * scale
            x0 = xh + vx * tiny
            y0 = yh + vy * tiny
            z0 = zh + vz * tiny
            bounces_done = ib + 1
        else:
            for ih in range(nh):
                if streak[ih] >= min_streak:
                    strict_seq[ih] += streak[ih]
                streak[ih] = 0.0
            bounces_done = n_bounces

        Lf = angular_momentum_axis_nb(x0, y0, z0, vx, vy, vz, object_id)
        out[i, 2] = exit_flag
        out[i, 3] = nohit_flag
        out[i, 4] = bounces_done
        out[i, 5] = abs(Lf - L0) / max(refR, 1e-12)
        out[i, 6] = abs(math.sqrt(vx * vx + vy * vy + vz * vz) - 1.0)
        out[i, 7] = actual_reflections
        base = 8
        for ih in range(nh):
            out[i, base + 6 * ih + 0] = old_refl[ih]
            out[i, base + 6 * ih + 1] = strict_seq[ih]
            out[i, base + 6 * ih + 2] = max_streak[ih]
            out[i, base + 6 * ih + 3] = 1.0 if max_streak[ih] >= min_streak else 0.0
            out[i, base + 6 * ih + 4] = old_all[ih]
            out[i, base + 6 * ih + 5] = focal_exit[ih]
    return out


# Python helper functions for sampling and volume proxies

def rmax_py(object_name, s, a, b, R, geom_type_name, body_model_name):
    return float(rmax_axis_nb(OBJECT_ID[object_name], float(s), float(a), float(b), float(R),
                              GEOM_ID[geom_type_name], BODY_ID[body_model_name]))


def drds_py(object_name, s, a, b, R, geom_type_name, body_model_name):
    return float(drds_axis_nb(OBJECT_ID[object_name], float(s), float(a), float(b), float(R),
                              GEOM_ID[geom_type_name], BODY_ID[body_model_name]))


def normal_py_from_axis(object_name, s, phi, a, b, R, geom_type_name, body_model_name):
    rr = max(0.0, rmax_py(object_name, s, a, b, R, geom_type_name, body_model_name))
    if object_name == "vertical":
        x = s
        y = rr * math.cos(phi)
        z = rr * math.sin(phi)
    else:
        x = rr * math.cos(phi)
        y = s
        z = rr * math.sin(phi)
    nx, ny, nz = normal_nb(x, y, z, OBJECT_ID[object_name], a, b, R, GEOM_ID[geom_type_name], BODY_ID[body_model_name])
    return np.array([x, y, z], dtype=np.float64), np.array([nx, ny, nz], dtype=np.float64)


def surface_point_arrays(object_name, a, b, R, geom_type_name, body_model_name, n_points=100):
    E = float(axis_extent_nb(OBJECT_ID[object_name], a, b, R, GEOM_ID[geom_type_name], BODY_ID[body_model_name]))
    eps = 1e-6 * max(E, 1.0)
    ss = np.linspace(-E + eps, E - eps, 8000)
    rs = np.array([max(0.0, rmax_py(object_name, s, a, b, R, geom_type_name, body_model_name)) for s in ss])
    rps = np.array([drds_py(object_name, s, a, b, R, geom_type_name, body_model_name) for s in ss])
    weights = rs * np.sqrt(1.0 + np.minimum(rps * rps, 1e8))
    weights[~np.isfinite(weights)] = 0.0
    weights = np.maximum(weights, 0.0)
    cdf = np.cumsum(weights)
    if cdf[-1] <= 0.0:
        sq = np.linspace(-0.9 * E, 0.9 * E, n_points)
    else:
        cdf = cdf / cdf[-1]
        qs = (np.arange(n_points) + 0.5) / n_points
        sq = np.interp(qs, cdf, ss)
    pts = np.zeros((n_points, 3), dtype=np.float64)
    norms = np.zeros((n_points, 3), dtype=np.float64)
    for i, s in enumerate(sq):
        phi = (i * GOLD) % (2 * PI)
        p, n = normal_py_from_axis(object_name, float(s), float(phi), a, b, R, geom_type_name, body_model_name)
        pts[i] = p
        norms[i] = n
    return pts, norms


def random_unit_sphere(n, rng):
    u = rng.random(n)
    v = rng.random(n)
    zc = 1.0 - 2.0 * u
    rc = np.sqrt(np.maximum(0.0, 1.0 - zc * zc))
    phi = 2.0 * PI * v
    return np.column_stack([rc * np.cos(phi), zc, rc * np.sin(phi)])


def quasi_unit_sphere(n):
    arr = np.zeros((n, 3), dtype=np.float64)
    for i in range(n):
        u = (i + 0.5) / n
        zc = 1.0 - 2.0 * u
        rc = math.sqrt(max(0.0, 1.0 - zc * zc))
        phi = (i * GOLD) % (2 * PI)
        arr[i] = [rc * math.cos(phi), zc, rc * math.sin(phi)]
    return arr


def make_initial_conditions(object_name, source_name, a, b, R, geom_type_name, body_model_name,
                            n_center, n_surface_points, n_dirs, sampling_mode="random", seed=1):
    rng = np.random.default_rng(seed)
    if source_name == "center_point_100":
        pos = np.zeros((n_center, 3), dtype=np.float64)
        if sampling_mode == "random":
            dirs = random_unit_sphere(n_center, rng)
        else:
            dirs = quasi_unit_sphere(n_center)
        if object_name == "vertical":
            # Use X as main axis: map generated z-axis component to vx.
            vel = np.column_stack([dirs[:, 1], dirs[:, 0], dirs[:, 2]]).astype(np.float64)
        else:
            # Use Y as main axis: generated dirs already have y as central column.
            vel = dirs.astype(np.float64)
        return pos, vel

    pts, norms = surface_point_arrays(object_name, a, b, R, geom_type_name, body_model_name, n_surface_points)
    n = n_surface_points * n_dirs
    pos = np.zeros((n, 3), dtype=np.float64)
    vel = np.zeros((n, 3), dtype=np.float64)
    E = float(axis_extent_nb(OBJECT_ID[object_name], a, b, R, GEOM_ID[geom_type_name], BODY_ID[body_model_name]))
    refR = float(reference_radius_nb(OBJECT_ID[object_name], a, b, R, GEOM_ID[geom_type_name], BODY_ID[body_model_name]))
    scale = max(E, refR, 1.0)
    idx = 0
    for ip in range(n_surface_points):
        p = pts[ip]
        outward = norms[ip]
        inward = -outward
        ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if abs(np.dot(ref, inward)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        t1 = np.cross(inward, ref)
        t1 = t1 / np.linalg.norm(t1)
        t2 = np.cross(inward, t1)
        t2 = t2 / np.linalg.norm(t2)
        for j in range(n_dirs):
            if sampling_mode == "random":
                cos_th = rng.random()
                phi = 2.0 * PI * rng.random()
            else:
                cos_th = (j + 0.5) / n_dirs
                phi = ((j * GOLD) + (ip * 0.371)) % (2 * PI)
            sin_th = math.sqrt(max(0.0, 1.0 - cos_th * cos_th))
            v = cos_th * inward + sin_th * (math.cos(phi) * t1 + math.sin(phi) * t2)
            v = v / np.linalg.norm(v)
            pos[idx] = p + v * (1e-7 * scale)
            vel[idx] = v
            idx += 1
    return pos, vel


def body_volume_proxy(object_name, a, b, R, geom_type_name, body_model_name, n=4000):
    E = float(axis_extent_nb(OBJECT_ID[object_name], a, b, R, GEOM_ID[geom_type_name], BODY_ID[body_model_name]))
    ss = np.linspace(-E, E, n)
    rs = np.array([max(0.0, rmax_py(object_name, s, a, b, R, geom_type_name, body_model_name)) for s in ss])
    vol = np.trapezoid(PI * rs * rs, ss)
    return float(vol)


def focal_zone_volume_proxy(object_name, a, b, R, h_abs):
    # Operational volume proxy, used only for comparison of zone scale.
    # vertical: cylinder slab height 2a and annulus R±h_abs.
    # horizontal: two disks radius a, half-thickness h_abs each.
    if object_name == "vertical":
        r1 = max(0.0, R - h_abs)
        r2 = R + h_abs
        return float((2.0 * a) * PI * (r2 * r2 - r1 * r1))
    return float(2.0 * (2.0 * h_abs) * PI * a * a)


def summarize_ray_output(ray_out, meta, h_over_a_values, h_abs_values, n_bounces):
    rows = []
    actual_ref_all = ray_out[:, 7]
    for ih, h_abs in enumerate(h_abs_values):
        base = 8 + 6 * ih
        for jlabel, jid in [("all", -1), ("low-J", 0), ("mid-J", 1), ("high-J", 2)]:
            if jid == -1:
                mask = np.ones(ray_out.shape[0], dtype=bool)
            else:
                mask = (ray_out[:, 1].astype(int) == jid)
            sub = ray_out[mask]
            row = meta.copy()
            row.update({
                "h_over_a": float(h_over_a_values[ih]),
                "h_abs": float(h_abs),
                "J_class": jlabel,
                "n_rays": int(len(sub)),
            })
            if len(sub) == 0:
                fields = [
                    "Jhat_mean", "exit_pct", "nohit_pct", "stay_pct", "mean_bounces_done",
                    "actual_reflections_mean", "actual_reflections_sum",
                    "old_reflective_density_per_100_bounces_pct",
                    "strict_sequence_density_per_100_bounces_pct",
                    "old_reflective_density_per_actual_reflection_pct",
                    "strict_sequence_density_per_actual_reflection_pct",
                    "strict_sequence_density_mean_per_ray_actual_pct",
                    "qualified_rays_pct", "focal_exit_contacts_per_100_bounces_pct",
                    "rays_with_focal_exit_pct", "old_reflective_contacts_mean",
                    "strict_sequence_contacts_mean", "max_reflective_streak_mean",
                    "max_reflective_streak_max", "rays_with_ge1_refl_contact_pct",
                    "rays_with_ge2_consecutive_pct", "rays_with_ge5_consecutive_pct",
                    "rays_with_ge10_consecutive_pct", "max_L_drift_norm", "max_speed_drift",
                    "strict_sequence_density_per_actual_reflection_ray_ci95_low",
                    "strict_sequence_density_per_actual_reflection_ray_ci95_high",
                    "strict_sequence_density_per_actual_reflection_ray_ci95_halfwidth",
                ]
                for f in fields:
                    row[f] = np.nan
                rows.append(row)
                continue

            old_refl = sub[:, base + 0]
            strict_seq = sub[:, base + 1]
            max_streak = sub[:, base + 2]
            qual = sub[:, base + 3] > 0.5
            old_all = sub[:, base + 4]
            focal_exit = sub[:, base + 5]
            actual_ref = sub[:, 7]
            actual_sum = float(np.sum(actual_ref))

            old_per_100 = float(np.mean(old_refl) / n_bounces * 100.0)
            strict_per_100 = float(np.mean(strict_seq) / n_bounces * 100.0)
            if actual_sum > 0.0:
                old_per_actual = float(np.sum(old_refl) / actual_sum * 100.0)
                strict_per_actual = float(np.sum(strict_seq) / actual_sum * 100.0)
            else:
                old_per_actual = np.nan
                strict_per_actual = np.nan
            per_ray = np.zeros_like(strict_seq)
            ok = actual_ref > 0.0
            per_ray[ok] = strict_seq[ok] / actual_ref[ok] * 100.0
            per_ray_ci_low = np.nan
            per_ray_ci_high = np.nan
            per_ray_ci_half = np.nan
            if np.any(ok):
                vals_ci = per_ray[ok]
                if len(vals_ci) > 1:
                    per_ray_ci_half = float(1.96 * np.std(vals_ci, ddof=1) / math.sqrt(len(vals_ci)))
                    per_ray_ci_low = float(np.mean(vals_ci) - per_ray_ci_half)
                    per_ray_ci_high = float(np.mean(vals_ci) + per_ray_ci_half)
                else:
                    per_ray_ci_low = float(vals_ci[0])
                    per_ray_ci_high = float(vals_ci[0])
                    per_ray_ci_half = 0.0

            row.update({
                "Jhat_mean": float(np.mean(sub[:, 0])),
                "exit_pct": float(np.mean(sub[:, 2]) * 100.0),
                "nohit_pct": float(np.mean(sub[:, 3]) * 100.0),
                "stay_pct": float(100.0 - np.mean(sub[:, 2]) * 100.0 - np.mean(sub[:, 3]) * 100.0),
                "mean_bounces_done": float(np.mean(sub[:, 4])),
                "actual_reflections_mean": float(np.mean(actual_ref)),
                "actual_reflections_sum": actual_sum,
                "old_reflective_density_per_100_bounces_pct": old_per_100,
                "strict_sequence_density_per_100_bounces_pct": strict_per_100,
                "old_reflective_density_per_actual_reflection_pct": old_per_actual,
                "strict_sequence_density_per_actual_reflection_pct": strict_per_actual,
                "strict_sequence_density_mean_per_ray_actual_pct": float(np.mean(per_ray[ok])) if np.any(ok) else np.nan,
                "qualified_rays_pct": float(np.mean(qual) * 100.0),
                "focal_exit_contacts_per_100_bounces_pct": float(np.mean(focal_exit) / n_bounces * 100.0),
                "rays_with_focal_exit_pct": float(np.mean(focal_exit > 0.0) * 100.0),
                "old_all_contacts_per_100_bounces_pct": float(np.mean(old_all) / n_bounces * 100.0),
                "old_reflective_contacts_mean": float(np.mean(old_refl)),
                "strict_sequence_contacts_mean": float(np.mean(strict_seq)),
                "max_reflective_streak_mean": float(np.mean(max_streak)),
                "max_reflective_streak_max": float(np.max(max_streak)),
                "rays_with_ge1_refl_contact_pct": float(np.mean(old_refl >= 1.0) * 100.0),
                "rays_with_ge2_consecutive_pct": float(np.mean(max_streak >= 2.0) * 100.0),
                "rays_with_ge5_consecutive_pct": float(np.mean(max_streak >= 5.0) * 100.0),
                "rays_with_ge10_consecutive_pct": float(np.mean(max_streak >= 10.0) * 100.0),
                "max_L_drift_norm": float(np.max(sub[:, 5])),
                "max_speed_drift": float(np.max(sub[:, 6])),
                "strict_sequence_density_per_actual_reflection_ray_ci95_low": per_ray_ci_low,
                "strict_sequence_density_per_actual_reflection_ray_ci95_high": per_ray_ci_high,
                "strict_sequence_density_per_actual_reflection_ray_ci95_halfwidth": per_ray_ci_half,
            })
            rows.append(row)
    return rows


def estimate_volumes_monte_carlo(object_name, a, b, R, geom_type_name, body_model_name,
                                 h_abs_values, n_samples=20000, seed=12345):
    """Monte Carlo estimate of body volume and focal-zone volume.

    This is a service diagnostic, not the primary ray metric. It estimates how large
    the operational focal zone is relative to the body for each h/a.
    """
    if n_samples <= 0:
        return np.nan, np.full(len(h_abs_values), np.nan, dtype=float)
    rng = np.random.default_rng(int(seed))
    oid = OBJECT_ID[object_name]
    gid = GEOM_ID[geom_type_name]
    bid = BODY_ID[body_model_name]
    E = float(axis_extent_nb(oid, float(a), float(b), float(R), gid, bid))
    ss_grid = np.linspace(-E, E, 3000)
    r_grid = np.array([max(0.0, rmax_axis_nb(oid, float(s), float(a), float(b), float(R), gid, bid)) for s in ss_grid])
    rbox = float(np.max(r_grid)) if len(r_grid) else float(R)
    rbox = max(rbox, 1e-9)
    s_vals = rng.uniform(-E, E, n_samples)
    u_vals = rng.uniform(-rbox, rbox, n_samples)
    v_vals = rng.uniform(-rbox, rbox, n_samples)
    inside_count = 0
    focal_counts = np.zeros(len(h_abs_values), dtype=np.int64)
    for i in range(n_samples):
        if oid == 0:
            x, y, z = s_vals[i], u_vals[i], v_vals[i]
        else:
            x, y, z = u_vals[i], s_vals[i], v_vals[i]
        if inside_nb(float(x), float(y), float(z), oid, float(a), float(b), float(R), gid, bid) > 0.0:
            inside_count += 1
            for ih, h_abs in enumerate(h_abs_values):
                if in_focal_zone_nb(float(x), float(y), float(z), oid, float(a), float(b), float(R), float(h_abs)):
                    focal_counts[ih] += 1
    vbox = (2.0 * E) * (2.0 * rbox) * (2.0 * rbox)
    body_v = vbox * inside_count / float(n_samples)
    focal_v = vbox * focal_counts.astype(float) / float(n_samples)
    return float(body_v), focal_v


def make_seed_statistics(df):
    """Aggregate per-seed rows into mean/std/CI95 across seeds."""
    metric_cols = [
        "strict_sequence_density_per_actual_reflection_pct",
        "strict_sequence_density_per_100_bounces_pct",
        "qualified_rays_pct",
        "exit_pct",
        "actual_reflections_mean",
        "max_reflective_streak_max",
        "focal_zone_volume_monte_carlo_fraction_pct",
    ]
    metric_cols = [c for c in metric_cols if c in df.columns]
    group_cols = [
        "object_type", "body_model", "geometry_type", "source", "geom_index", "previous_label",
        "R", "a", "b", "beta_b_over_a", "rho_R_over_a", "h_over_a", "h_abs",
        "n_bounces", "n_scan", "sampling_mode", "J_class"
    ]
    group_cols = [c for c in group_cols if c in df.columns]
    rows = []
    for key, g in df.groupby(group_cols, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        row = dict(zip(group_cols, key))
        row["n_seeds"] = int(g["seed"].nunique()) if "seed" in g.columns else 1
        for c in metric_cols:
            vals = pd.to_numeric(g[c], errors="coerce").dropna().to_numpy(dtype=float)
            if len(vals) == 0:
                row[c + "_mean"] = np.nan
                row[c + "_std"] = np.nan
                row[c + "_ci95_low"] = np.nan
                row[c + "_ci95_high"] = np.nan
                row[c + "_ci95_halfwidth"] = np.nan
            else:
                mean = float(np.mean(vals))
                std = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
                half = float(1.96 * std / math.sqrt(len(vals))) if len(vals) > 1 else 0.0
                row[c + "_mean"] = mean
                row[c + "_std"] = std
                row[c + "_ci95_low"] = mean - half
                row[c + "_ci95_high"] = mean + half
                row[c + "_ci95_halfwidth"] = half
        rows.append(row)
    return pd.DataFrame(rows)


def step_ray_python(state, object_name, a, b, R, geom_type_name, body_model_name, n_scan):
    x, y, z, vx, vy, vz = state
    oid = OBJECT_ID[object_name]
    gid = GEOM_ID[geom_type_name]
    bid = BODY_ID[body_model_name]
    ok, t, xh, yh, zh = next_hit_nb(float(x), float(y), float(z), float(vx), float(vy), float(vz),
                                    oid, float(a), float(b), float(R), gid, bid, int(n_scan))
    if not ok:
        return False, state, t
    nx, ny, nz = normal_nb(float(xh), float(yh), float(zh), oid, float(a), float(b), float(R), gid, bid)
    vx2, vy2, vz2 = reflect_nb(float(vx), float(vy), float(vz), nx, ny, nz)
    E = float(axis_extent_nb(oid, float(a), float(b), float(R), gid, bid))
    refR = float(reference_radius_nb(oid, float(a), float(b), float(R), gid, bid))
    tiny = 1e-8 * max(E, refR, 1.0)
    new_state = np.array([xh + vx2*tiny, yh + vy2*tiny, zh + vz2*tiny, vx2, vy2, vz2], dtype=float)
    return True, new_state, float(t)


def lyapunov_estimate_python(ray_state, object_name, a, b, R, geom_type_name, body_model_name,
                             n_bounces=300, n_scan=64, eps0=1e-8, renorm_every=5):
    """Simple finite-distance Lyapunov diagnostic for the billiard map."""
    state_a = np.array(ray_state, dtype=float).copy()
    v = state_a[3:6]
    v = v / max(np.linalg.norm(v), 1e-15)
    state_a[3:6] = v
    # deterministic small perturbation orthogonal to velocity
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, v)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    dv = ref - np.dot(ref, v) * v
    dv = dv / max(np.linalg.norm(dv), 1e-15)
    state_b = state_a.copy()
    state_b[3:6] = state_b[3:6] + eps0 * dv
    state_b[3:6] /= max(np.linalg.norm(state_b[3:6]), 1e-15)
    sum_log = 0.0
    total_t = 0.0
    n_renorm = 0
    for k in range(int(n_bounces)):
        ok_a, state_a, ta = step_ray_python(state_a, object_name, a, b, R, geom_type_name, body_model_name, n_scan)
        ok_b, state_b, tb = step_ray_python(state_b, object_name, a, b, R, geom_type_name, body_model_name, n_scan)
        if not ok_a or not ok_b:
            break
        total_t += 0.5 * (ta + tb)
        if (k + 1) % renorm_every == 0:
            diff = state_b - state_a
            dist = float(np.linalg.norm(diff))
            if dist <= 0 or not np.isfinite(dist):
                break
            sum_log += math.log(dist / eps0)
            n_renorm += 1
            state_b = state_a + diff * (eps0 / dist)
            state_b[3:6] /= max(np.linalg.norm(state_b[3:6]), 1e-15)
    if total_t <= 0.0 or n_renorm == 0:
        return np.nan
    return float(sum_log / total_t)


def poincare_section_python(ray_state, object_name, a, b, R, geom_type_name, body_model_name,
                            n_bounces=300, n_scan=64):
    """Returns a compact Poincare-like section: azimuth angle and tangential velocity at hits."""
    state = np.array(ray_state, dtype=float).copy()
    state[3:6] /= max(np.linalg.norm(state[3:6]), 1e-15)
    pts = []
    for _ in range(int(n_bounces)):
        ok, state, _t = step_ray_python(state, object_name, a, b, R, geom_type_name, body_model_name, n_scan)
        if not ok:
            break
        x, y, z, vx, vy, vz = state
        if object_name == "vertical":
            rho = math.sqrt(y*y + z*z) + 1e-15
            phi = math.atan2(z, y)
            vt = (-z/rho)*vy + (y/rho)*vz
        else:
            rho = math.sqrt(x*x + z*z) + 1e-15
            phi = math.atan2(z, x)
            vt = (-z/rho)*vx + (x/rho)*vz
        pts.append((phi, vt))
    return np.array(pts, dtype=float)


def run_optional_diagnostics(outdir, df):
    """Run Lyapunov/Poincare diagnostics for top corrected-metric cases."""
    if not RUN_DIAGNOSTICS:
        return
    metric = "strict_sequence_density_per_actual_reflection_pct"
    work = df[(df["J_class"] != "all") & (df["body_model"] == "pseudo") & (df["n_rays"] > 0)].copy()
    if work.empty or metric not in work.columns:
        return
    diag_dir = Path(outdir) / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)
    top = work.sort_values(metric, ascending=False).head(int(DIAGNOSTIC_TOP_CASES))
    rows = []
    for idx, r in top.iterrows():
        object_name = str(r.object_type)
        geom_type_name = str(r.geometry_type)
        source_name = str(r.source)
        a = float(r.a); b = float(r.b); R = float(r.R)
        n_scan = int(r.n_scan)
        seed = int(r.seed)
        pos, vel = make_initial_conditions(object_name, source_name, a, b, R, geom_type_name, "pseudo",
                                           CENTER_RAYS, SURFACE_POINTS, DIRECTIONS_PER_POINT,
                                           sampling_mode=SAMPLING_MODE, seed=seed)
        n_take = min(int(DIAGNOSTIC_RAYS_PER_CASE), pos.shape[0])
        case_label = f"case{len(rows)+1}_{object_name}_G{int(r.geom_index)}_{geom_type_name}_{source_name}"
        poincare_saved = []
        for ir in range(n_take):
            ray_state = np.concatenate([pos[ir], vel[ir]])
            lam = lyapunov_estimate_python(ray_state, object_name, a, b, R, geom_type_name, "pseudo",
                                           n_bounces=int(DIAGNOSTIC_BOUNCES), n_scan=n_scan)
            rows.append({
                "case_label": case_label,
                "ray_index": ir,
                "object_type": object_name,
                "geom_index": int(r.geom_index),
                "geometry_type": geom_type_name,
                "source": source_name,
                "R": R, "a": a, "b": b,
                "h_over_a": float(r.h_over_a),
                "J_class": str(r.J_class),
                "main_metric_pct": float(r[metric]),
                "lyapunov_estimate": lam,
            })
            if ir < int(POINCARE_RAYS_PER_CASE):
                pts = poincare_section_python(ray_state, object_name, a, b, R, geom_type_name, "pseudo",
                                              n_bounces=int(DIAGNOSTIC_BOUNCES), n_scan=n_scan)
                poincare_saved.append(pts)
        if poincare_saved:
            np.savez(diag_dir / f"poincare_{case_label}.npz", *poincare_saved)
    pd.DataFrame(rows).to_csv(diag_dir / "lyapunov_top_cases.csv", index=False)


def make_figures(outdir, df):
    outdir = Path(outdir)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    work = df[(df["J_class"] != "all") & (df["body_model"] == "pseudo") & (df["n_rays"] > 0)].copy()
    if work.empty:
        return

    metric = "strict_sequence_density_per_actual_reflection_pct"
    top = work.sort_values(metric, ascending=False).head(25)
    fig, ax = plt.subplots(figsize=(13, 6))
    labels = [
        f"{r.object_type[:1].upper()} G{int(r.geom_index)} {r.geometry_type[0]} {r.source.split('_')[0]} h/a={r.h_over_a:.2f} {r.J_class}"
        for _, r in top.iterrows()
    ]
    ax.bar(range(len(top)), top[metric].fillna(0.0))
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("strict density, % of actual reflections")
    ax.set_title("Top cases by corrected strict metric")
    fig.tight_layout()
    fig.savefig(figdir / "top_cases_corrected_metric.png", dpi=180)
    plt.close(fig)

    # Mean corrected density by object and J class
    pivot = work.pivot_table(index="J_class", columns="object_type", values=metric, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    pivot.plot(kind="bar", ax=ax)
    ax.set_ylabel("mean strict density, % of actual reflections")
    ax.set_title("Corrected metric by object and J-class")
    fig.tight_layout()
    fig.savefig(figdir / "mean_corrected_metric_by_object_J.png", dpi=180)
    plt.close(fig)

    # h/a dependence
    htab = work.pivot_table(index="h_over_a", columns="object_type", values=metric, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    htab.plot(marker="o", ax=ax)
    ax.set_xlabel("h/a")
    ax.set_ylabel("mean strict density, % of actual reflections")
    ax.set_title("Dependence on normalized focal-zone thickness h/a")
    fig.tight_layout()
    fig.savefig(figdir / "h_over_a_dependence_corrected_metric.png", dpi=180)
    plt.close(fig)


def run_universal(outdir=OUTDIR, fast_test=False):
    t0 = time.time()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    object_types = list(OBJECT_TYPES)
    body_models = list(BODY_MODELS)
    geometry_types = list(GEOMETRY_TYPES)
    sources = list(SOURCES)
    h_over_a_values = list(H_OVER_A_VALUES)
    n_bounces_list = list(N_BOUNCES_LIST)
    n_scan_list = list(N_SCAN_LIST)
    seeds = list(RANDOM_SEEDS)
    center_rays = CENTER_RAYS
    surface_points = SURFACE_POINTS
    directions_per_point = DIRECTIONS_PER_POINT

    if fast_test:
        # Быстрый тест не заменяет научный расчёт. Он нужен только для проверки запуска.
        object_types = ["vertical", "horizontal"]
        body_models = ["pseudo"]
        geometry_types = ["closed"]
        sources = ["center_point_100", "surface_100x100"]
        h_over_a_values = [0.10, 0.30]
        n_bounces_list = [30]
        n_scan_list = [32]
        seeds = [1]
        center_rays = 20
        surface_points = 20
        directions_per_point = 20

    all_rows = []
    total_jobs = (len(object_types) * len(body_models) * len(geometry_types) * len(sources) *
                  len(BEST10) * len(n_bounces_list) * len(n_scan_list) * len(seeds))
    job = 0
    print("GVI universal corrected ray verification")
    print("Numba:", HAVE_NUMBA, "threads:", get_num_threads())
    print("Total jobs:", total_jobs)
    print("Output:", outdir)

    for object_name in object_types:
        for body_model_name in body_models:
            for geom_type_name in geometry_types:
                for source_name in sources:
                    for g in BEST10:
                        a = float(g["a"])
                        b = float(g["b"])
                        R = float(g["R"])
                        h_abs_values = np.array([hv * a for hv in h_over_a_values], dtype=np.float64)
                        body_vol = body_volume_proxy(object_name, a, b, R, geom_type_name, body_model_name)
                        fvols = np.array([focal_zone_volume_proxy(object_name, a, b, R, h) for h in h_abs_values])
                        mc_body_vol, mc_fvols = estimate_volumes_monte_carlo(
                            object_name, a, b, R, geom_type_name, body_model_name, h_abs_values,
                            n_samples=int(MONTE_CARLO_VOLUME_SAMPLES),
                            seed=int(MONTE_CARLO_VOLUME_SEED + int(g["geom_index"]) * 100 + BODY_ID[body_model_name] * 10 + GEOM_ID[geom_type_name])
                        )

                        for n_bounces in n_bounces_list:
                            for n_scan in n_scan_list:
                                for seed in seeds:
                                    job += 1
                                    print(f"[{job}/{total_jobs}] {object_name} {body_model_name} {geom_type_name} "
                                          f"{source_name} G{g['geom_index']} seed={seed} "
                                          f"N={n_bounces} scan={n_scan}", flush=True)
                                    pos, vel = make_initial_conditions(
                                        object_name, source_name, a, b, R, geom_type_name, body_model_name,
                                        center_rays, surface_points, directions_per_point,
                                        sampling_mode=SAMPLING_MODE, seed=int(seed)
                                    )
                                    ray_out = simulate_batch_nb(
                                        pos, vel, OBJECT_ID[object_name], a, b, R, GEOM_ID[geom_type_name],
                                        BODY_ID[body_model_name], h_abs_values, int(n_bounces), int(n_scan),
                                        float(LOW_J_THRESHOLD), float(MID_J_THRESHOLD), int(MIN_STRICT_STREAK)
                                    )
                                    meta = dict(
                                        object_type=object_name,
                                        body_model=body_model_name,
                                        geometry_type=geom_type_name,
                                        source=source_name,
                                        geom_index=int(g["geom_index"]),
                                        previous_label=g["previous_label"],
                                        R=R,
                                        a=a,
                                        b=b,
                                        beta_b_over_a=b / a,
                                        rho_R_over_a=R / a,
                                        n_bounces=int(n_bounces),
                                        n_scan=int(n_scan),
                                        seed=int(seed),
                                        sampling_mode=SAMPLING_MODE,
                                        min_strict_streak=int(MIN_STRICT_STREAK),
                                        body_volume_proxy=body_vol,
                                    )
                                    rows = summarize_ray_output(ray_out, meta, h_over_a_values, h_abs_values, int(n_bounces))
                                    for r in rows:
                                        ih = h_over_a_values.index(r["h_over_a"])
                                        r["focal_zone_volume_proxy"] = float(fvols[ih])
                                        r["focal_zone_volume_proxy_fraction_pct"] = float(fvols[ih] / body_vol * 100.0) if body_vol > 0 else np.nan
                                        r["body_volume_monte_carlo"] = float(mc_body_vol)
                                        r["focal_zone_volume_monte_carlo"] = float(mc_fvols[ih]) if np.ndim(mc_fvols) > 0 else np.nan
                                        r["focal_zone_volume_monte_carlo_fraction_pct"] = (float(mc_fvols[ih] / mc_body_vol * 100.0) if np.isfinite(mc_body_vol) and mc_body_vol > 0 else np.nan)
                                    all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    all_path = outdir / "all_results.csv"
    df.to_csv(all_path, index=False)

    compact = df[(df["body_model"] == "pseudo") & (df["J_class"] != "all")].copy()
    compact.to_csv(outdir / "compact_pseudo_J_only.csv", index=False)

    metric = "strict_sequence_density_per_actual_reflection_pct"
    summary = df[df["J_class"] != "all"].pivot_table(
        index=["object_type", "body_model", "J_class"],
        values=[metric, "strict_sequence_density_per_100_bounces_pct", "qualified_rays_pct", "exit_pct", "actual_reflections_mean"],
        aggfunc=["mean", "max"]
    )
    summary.to_csv(outdir / "summary_by_object_j.csv")

    top = df[(df["J_class"] != "all") & (df["body_model"] == "pseudo")].sort_values(metric, ascending=False).head(100)
    top.to_csv(outdir / "top_cases_corrected_metric.csv", index=False)

    seed_stats = make_seed_statistics(df)
    seed_stats.to_csv(outdir / "seed_statistics_ci95.csv", index=False)

    metadata = {
        "created_unix_time": time.time(),
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest() if "__file__" in globals() and Path(__file__).exists() else "n/a",
        "elapsed_seconds": time.time() - t0,
        "object_types": object_types,
        "body_models": body_models,
        "geometry_types": geometry_types,
        "sources": sources,
        "h_over_a_values": h_over_a_values,
        "h_definition": "h_abs = h_over_a * a; vertical |rho-R|<=h_abs; horizontal ||y|-R|<=h_abs and rho<=a",
        "n_bounces_list": n_bounces_list,
        "n_scan_list": n_scan_list,
        "center_rays": center_rays,
        "surface_points": surface_points,
        "directions_per_point": directions_per_point,
        "sampling_mode": SAMPLING_MODE,
        "random_seeds": seeds,
        "LOW_J_THRESHOLD": LOW_J_THRESHOLD,
        "MID_J_THRESHOLD": MID_J_THRESHOLD,
        "MIN_STRICT_STREAK": MIN_STRICT_STREAK,
        "MONTE_CARLO_VOLUME_SAMPLES": MONTE_CARLO_VOLUME_SAMPLES,
        "MONTE_CARLO_VOLUME_SEED": MONTE_CARLO_VOLUME_SEED,
        "BOOTSTRAP_RESAMPLES": BOOTSTRAP_RESAMPLES,
        "RUN_DIAGNOSTICS": RUN_DIAGNOSTICS,
        "DIAGNOSTIC_TOP_CASES": DIAGNOSTIC_TOP_CASES,
        "DIAGNOSTIC_RAYS_PER_CASE": DIAGNOSTIC_RAYS_PER_CASE,
        "DIAGNOSTIC_BOUNCES": DIAGNOSTIC_BOUNCES,
        "main_metric": metric,
        "legacy_metric": "strict_sequence_density_per_100_bounces_pct",
        "warning": "Ray-billiard model only: no wavelength, phase, interference, diffraction or full-wave resonance.",
        "control_figures_removed": True
    }
    with open(outdir / "run_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    make_figures(outdir, df)
    run_optional_diagnostics(outdir, df)

    print("\nDone.")
    print("Rows:", len(df))
    print("Main metric:", metric)
    print("Saved:", all_path)
    print("Elapsed seconds:", round(time.time() - t0, 2))
    return df


def _parse_csv_list(value):
    if value is None or str(value).strip() == "":
        return None
    return [x.strip() for x in str(value).split(",") if x.strip()]


def _parse_float_list(value):
    xs = _parse_csv_list(value)
    return None if xs is None else [float(x) for x in xs]


def _parse_int_list(value):
    xs = _parse_csv_list(value)
    return None if xs is None else [int(x) for x in xs]


def main():
    parser = argparse.ArgumentParser(description="GVI pseudohyperboloid-only corrected ray-billiard verification; control figures are not included")
    parser.add_argument("--outdir", default=OUTDIR)
    parser.add_argument("--fast-test", action="store_true", help="Small quick run for checking script execution")
    parser.add_argument("--diagnostics", action="store_true", help="Run optional Lyapunov/Poincare diagnostics for top cases")
    parser.add_argument("--object-types", default=None, help="Comma list: vertical,horizontal")
    parser.add_argument("--geometry-types", default=None, help="Comma list: open,closed")
    parser.add_argument("--sources", default=None, help="Comma list: center_point_100,surface_100x100")
    parser.add_argument("--geom-indices", default=None, help="Comma list of BEST10 geometry indices, e.g. 1,6")
    parser.add_argument("--h-over-a", default=None, help="Comma list, e.g. 0.05,0.10,0.20,0.30")
    parser.add_argument("--seeds", default=None, help="Comma list, e.g. 1,2,3,4,5")
    parser.add_argument("--n-bounces", default=None, help="Comma list, e.g. 1000")
    parser.add_argument("--n-scan", default=None, help="Comma list, e.g. 64")
    parser.add_argument("--center-rays", type=int, default=None)
    parser.add_argument("--surface-points", type=int, default=None)
    parser.add_argument("--directions-per-point", type=int, default=None)
    parser.add_argument("--mc-samples", type=int, default=None, help="Monte Carlo samples for volume estimates; 0 disables MC")
    # parse_known_args avoids Jupyter/Colab -f argument crash
    args, _unknown = parser.parse_known_args()

    global RUN_DIAGNOSTICS, OBJECT_TYPES, GEOMETRY_TYPES, SOURCES, BEST10, H_OVER_A_VALUES
    global RANDOM_SEEDS, N_BOUNCES_LIST, N_SCAN_LIST, CENTER_RAYS, SURFACE_POINTS
    global DIRECTIONS_PER_POINT, MONTE_CARLO_VOLUME_SAMPLES, BODY_MODELS

    BODY_MODELS = ["pseudo"]  # hard guard: no control bodies in this script
    if args.diagnostics:
        RUN_DIAGNOSTICS = True
    parsed = _parse_csv_list(args.object_types)
    if parsed is not None:
        OBJECT_TYPES = parsed
    parsed = _parse_csv_list(args.geometry_types)
    if parsed is not None:
        GEOMETRY_TYPES = parsed
    parsed = _parse_csv_list(args.sources)
    if parsed is not None:
        SOURCES = parsed
    parsed = _parse_int_list(args.geom_indices)
    if parsed is not None:
        keep = set(parsed)
        BEST10 = [g for g in BEST10 if int(g["geom_index"]) in keep]
    parsed = _parse_float_list(args.h_over_a)
    if parsed is not None:
        H_OVER_A_VALUES = parsed
    parsed = _parse_int_list(args.seeds)
    if parsed is not None:
        RANDOM_SEEDS = parsed
    parsed = _parse_int_list(args.n_bounces)
    if parsed is not None:
        N_BOUNCES_LIST = parsed
    parsed = _parse_int_list(args.n_scan)
    if parsed is not None:
        N_SCAN_LIST = parsed
    if args.center_rays is not None:
        CENTER_RAYS = args.center_rays
    if args.surface_points is not None:
        SURFACE_POINTS = args.surface_points
    if args.directions_per_point is not None:
        DIRECTIONS_PER_POINT = args.directions_per_point
    if args.mc_samples is not None:
        MONTE_CARLO_VOLUME_SAMPLES = args.mc_samples

    run_universal(args.outdir, fast_test=args.fast_test)


if __name__ == "__main__":
    main()

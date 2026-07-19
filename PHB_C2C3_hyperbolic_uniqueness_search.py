# -*- coding: utf-8 -*-
"""
PHB_C2C3_hyperbolic_uniqueness_search.py

Purpose
-------
Search for wave signatures that could distinguish the closed second-order
vertical pseudohyperboloid (PHB) from non-hyperbolic controls, beyond the
usual annular-energy metric eta(|x|<=c, |rho-R|<=0.10R).

Why this script exists
----------------------
The previous C2/C3 metric confirmed strong high-m annular localization, but
it did not prove uniqueness of the hyperbolic generatrix: a torus-cap + linear
horn reference can give nearly the same eta in the ±10%R annulus. This script
therefore adds metrics designed to detect the expected wave imprint of the
hyperbolic focal billiard:

  1) axial inter-focal confinement:       E(|x|<=c) / E(|x|>c)
  2) horn leakage suppression:            E(|x|>c)
  3) focal-plane barrier/drop ratio:       energy just inside x=±c / energy just outside x=±c
  4) matched PHB-vs-control comparisons at similar kR, same m, same boundary proxy
  5) C3-like clusters using joint conditions: eta_100R >= threshold AND CF_axial >= threshold

Model status
------------
Reduced boundary-fitted axisymmetric scalar Helmholtz eigenproblem:
    -Delta_m U = k^2 U,
with u(x,rho,phi)=U_m(x,rho) exp(i m phi).

This is NOT full-vector Maxwell, NOT a true vector TE/TM calculation, NOT an
open-output laser calculation, and NOT ray tracing. Natural/Dirichlet boundary
conditions are scalar TE-like/TM-like proxies only.

Default target
--------------
    a = 0.3, b = 0.6, R = 3.0
    kR in [10,25], step 0.5
    m = 15

Recommended commands
--------------------
Quick smoke test:
    python PHB_C2C3_hyperbolic_uniqueness_search.py --preset quick

Target run for the current article:
    python PHB_C2C3_hyperbolic_uniqueness_search.py --preset target

Broader search over m and controls:
    python PHB_C2C3_hyperbolic_uniqueness_search.py --preset broad

Fuller but slower run:
    python PHB_C2C3_hyperbolic_uniqueness_search.py --preset full
"""
from __future__ import annotations

import argparse
import json
import math
import time
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S_MIN = 2e-3
FOCAL_BANDS = [0.025, 0.050, 0.100]
DEDUP_KR_TOL = 0.06
EPS = 1e-300

# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------

def common_length(a: float, b: float, R: float) -> float:
    return float(a * math.sqrt(1.0 + (R / b) ** 2))


def _profile_horn_by_t(absx: np.ndarray, a: float, b: float, R: float, L: float, kind: str) -> np.ndarray:
    """Horn profile for absx>a. All non-PHB controls share the same L."""
    x = np.asarray(absx, dtype=float)
    t = np.clip((x - a) / max(L - a, 1e-14), 0.0, 1.0)

    if kind == "PH":
        return np.maximum(0.0, R - b * np.sqrt(np.maximum(0.0, (x / a) ** 2 - 1.0)))
    if kind == "LINEAR":
        return np.maximum(0.0, R * (1.0 - t))
    if kind.startswith("POLY"):
        # POLY0p5, POLY1p5, POLY2, POLY3 etc.
        pstr = kind.replace("POLY", "").replace("p", ".")
        p = float(pstr)
        return np.maximum(0.0, R * (1.0 - t ** p))
    if kind == "SMOOTHSTEP":
        s = 3.0 * t**2 - 2.0 * t**3
        return np.maximum(0.0, R * (1.0 - s))
    if kind == "CIRCULAR":
        # Rounded control: same endpoints, circular-like decay.
        return np.maximum(0.0, R * np.sqrt(np.maximum(0.0, 1.0 - t**2)))
    if kind == "COSINE":
        return np.maximum(0.0, R * 0.5 * (1.0 + np.cos(np.pi * t)))
    raise ValueError(f"Unknown geometry kind: {kind}")


def wall_profile(xs: np.ndarray, kind: str, a: float, b: float, R: float, L: Optional[float] = None) -> np.ndarray:
    """Outer wall radius r_wall(x) for the mapped meridional domain."""
    if L is None:
        L = common_length(a, b, R)
    x = np.abs(np.asarray(xs, dtype=float))
    out = np.zeros_like(x)
    tor = x <= a
    # Shared half-toroidal cap: this is the deliberate closed-cavity common part.
    out[tor] = R + np.sqrt(np.maximum(0.0, a * a - x[tor] ** 2))
    horn = (x > a) & (x <= L)
    out[horn] = _profile_horn_by_t(x[horn], a, b, R, L, kind)
    return out


def wall_derivative(xs: np.ndarray, kind: str, a: float, b: float, R: float, L: Optional[float] = None) -> np.ndarray:
    if L is None:
        L = common_length(a, b, R)
    h = 1e-5 * max(L, R, a, b)
    xs = np.asarray(xs, dtype=float)
    x1 = np.clip(xs - h, -L + 1e-8, L - 1e-8)
    x2 = np.clip(xs + h, -L + 1e-8, L - 1e-8)
    r1 = wall_profile(x1, kind, a, b, R, L)
    r2 = wall_profile(x2, kind, a, b, R, L)
    d = (r2 - r1) / np.maximum(x2 - x1, 1e-15)
    return np.clip(d, -300.0, 300.0)


def geometry_labels() -> Dict[str, str]:
    return {
        "PH": "PHB hyperbolic horns",
        "LINEAR": "torus cap + linear horns",
        "POLY0p5": "torus cap + polynomial p=0.5 horns",
        "POLY1p5": "torus cap + polynomial p=1.5 horns",
        "POLY2": "torus cap + polynomial p=2 horns",
        "POLY3": "torus cap + polynomial p=3 horns",
        "SMOOTHSTEP": "torus cap + smoothstep horns",
        "CIRCULAR": "torus cap + circular-like horns",
        "COSINE": "torus cap + cosine horns",
    }

# -----------------------------------------------------------------------------
# Sparse operators
# -----------------------------------------------------------------------------

def derivative_matrix(n: int, h: float) -> sp.csr_matrix:
    rows, cols, vals = [], [], []
    for i in range(n):
        if i == 0:
            rows += [i, i]; cols += [0, 1]; vals += [-1 / h, 1 / h]
        elif i == n - 1:
            rows += [i, i]; cols += [n - 2, n - 1]; vals += [-1 / h, 1 / h]
        else:
            rows += [i, i]; cols += [i - 1, i + 1]; vals += [-0.5 / h, 0.5 / h]
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def build_problem(Nx: int, Ns: int, m: int, kind: str, bc: str, a: float, b: float, R: float) -> dict:
    L = common_length(a, b, R)
    xmin = -L + 1e-4 * max(L, R)
    xmax =  L - 1e-4 * max(L, R)
    xs = np.linspace(xmin, xmax, Nx)
    ss = np.linspace(S_MIN, 1.0, Ns)
    dx = (xmax - xmin) / (Nx - 1)
    ds = (1.0 - S_MIN) / (Ns - 1)

    Dx1 = derivative_matrix(Nx, dx)
    Ds1 = derivative_matrix(Ns, ds)
    Dx = sp.kron(Dx1, sp.eye(Ns, format="csr"), format="csr")
    Ds = sp.kron(sp.eye(Nx, format="csr"), Ds1, format="csr")

    X, S = np.meshgrid(xs, ss, indexing="ij")
    rw = wall_profile(xs, kind, a, b, R, L)[:, None]
    rpx = wall_derivative(xs, kind, a, b, R, L)[:, None]

    alpha = (S * rpx / np.maximum(rw, 1e-10)).ravel()
    W = (S * rw * rw).ravel()
    invr = (1.0 / np.maximum(rw, 1e-10)) * np.ones_like(S)
    area = dx * ds

    Gx = Dx - sp.diags(alpha, 0, format="csr") @ Ds
    Gr = sp.diags(invr.ravel(), 0, format="csr") @ Ds

    K = Gx.T @ sp.diags(W * area, 0, format="csr") @ Gx + Gr.T @ sp.diags(W * area, 0, format="csr") @ Gr
    if m > 0:
        K = K + sp.diags((m * m / np.maximum(S.ravel(), 1e-8)) * area, 0, format="csr")
    M = sp.diags(W * area, 0, format="csr")

    n = Nx * Ns
    fixed = np.zeros(n, dtype=bool)
    def idx(i: int, j: int) -> int:
        return i * Ns + j

    if bc == "natural_TE_like_scalar_proxy":
        # Regularity axis condition for m>0.
        if m > 0:
            fixed[[idx(i, 0) for i in range(Nx)]] = True
    elif bc == "dirichlet_TM_like_scalar_proxy":
        fixed[[idx(i, Ns - 1) for i in range(Nx)]] = True
        fixed[[idx(0, j) for j in range(Ns)]] = True
        fixed[[idx(Nx - 1, j) for j in range(Ns)]] = True
        if m > 0:
            fixed[[idx(i, 0) for i in range(Nx)]] = True
    else:
        raise ValueError(f"Unknown boundary proxy: {bc}")

    free = np.where(~fixed)[0]
    K = K[free][:, free]
    M = M[free][:, free]
    K = (K + K.T) * 0.5
    M = (M + M.T) * 0.5

    mdiag = M.diagonal()
    Sscale = sp.diags(1.0 / np.sqrt(np.maximum(mdiag, 1e-300)), 0, format="csr")
    A = Sscale @ K @ Sscale

    return {
        "A": A, "M": M, "Sscale": Sscale, "free": free,
        "xs": xs, "ss": ss, "X": X, "S": S, "rw": rw,
        "L": L, "Nx": Nx, "Ns": Ns, "dx": dx, "ds": ds,
        "kind": kind, "bc": bc, "a": a, "b": b, "R": R, "m": m,
    }

# -----------------------------------------------------------------------------
# Mode metrics
# -----------------------------------------------------------------------------

def _safe_ratio(num: float, den: float) -> float:
    return float(num / max(den, 1e-300))


def mode_metrics(vec_free: np.ndarray, eigval: float, prob: dict, target_kR: Optional[float] = None) -> Optional[dict]:
    Nx = prob["Nx"]; Ns = prob["Ns"]
    xs = prob["xs"]; X = prob["X"]; S = prob["S"]; rw = prob["rw"]
    a = prob["a"]; b = prob["b"]; R = prob["R"]
    kind = prob["kind"]; bc = prob["bc"]; m = prob["m"]; L = prob["L"]
    free = prob["free"]

    U = np.zeros(Nx * Ns)
    U[free] = vec_free
    U = U.reshape((Nx, Ns))
    dx = (xs[-1] - xs[0]) / (Nx - 1)
    ds = (1.0 - S_MIN) / (Ns - 1)

    RHO = S * rw
    W = S * rw * rw
    amp2 = U * U
    Etot = float(np.sum(amp2 * W) * dx * ds)
    Vtot = float(np.sum(W) * dx * ds)
    if Etot <= EPS:
        return None

    c = math.sqrt(a * a + b * b)
    k = math.sqrt(max(float(eigval), 0.0))
    kR = k * R

    mask_inter = np.abs(X) <= c
    mask_horn = np.abs(X) > c
    E_inter = float(np.sum(amp2 * W * mask_inter) * dx * ds / Etot)
    E_horn = float(np.sum(amp2 * W * mask_horn) * dx * ds / Etot)

    # Focal-plane barrier/drop.  The strip width is deliberately not one cell only;
    # it samples a small physical neighborhood on each side of |x|=c.
    strip_w = max(0.05 * c, 2.5 * dx)
    mask_inside_strip = (np.abs(X) >= (c - strip_w)) & (np.abs(X) <= c)
    mask_outside_strip = (np.abs(X) > c) & (np.abs(X) <= (c + strip_w))
    E_inside_strip = float(np.sum(amp2 * W * mask_inside_strip) * dx * ds / Etot)
    E_outside_strip = float(np.sum(amp2 * W * mask_outside_strip) * dx * ds / Etot)

    # Axial density moments.
    x_density = np.sum(amp2 * W, axis=1) * ds
    x_density_norm = x_density / max(float(np.sum(x_density) * dx), 1e-300)
    x_abs_mean = float(np.sum(np.abs(xs) * x_density_norm) * dx)
    x2_mean = float(np.sum(xs**2 * x_density_norm) * dx)
    sigma_x = math.sqrt(max(x2_mean, 0.0))

    near_wall = float(np.sum(amp2 * W * (S >= 0.9)) * dx * ds / Etot)
    center_waist = float(np.sum(amp2 * W * (np.abs(X) <= a)) * dx * ds / Etot)

    row = {
        "geometry": kind,
        "geometry_label": geometry_labels().get(kind, kind),
        "bc": bc,
        "a": a, "b": b, "R": R, "c": c, "L": L,
        "m": m,
        "kR": float(kR), "k_geom": float(k), "R_over_lambda": float(kR / (2 * math.pi)),
        "Nx": Nx, "Ns": Ns,
        "target_kR": float(target_kR) if target_kR is not None else np.nan,
        "E_interfocal_absx_le_c": E_inter,
        "E_horns_absx_gt_c": E_horn,
        "CF_axial_inter_over_horns": _safe_ratio(E_inter, E_horn),
        "horn_leakage": E_horn,
        "E_inside_focal_plane_strip": E_inside_strip,
        "E_outside_focal_plane_strip": E_outside_strip,
        "focal_plane_drop_ratio": _safe_ratio(E_inside_strip, E_outside_strip),
        "x_abs_mean": x_abs_mean,
        "sigma_x": sigma_x,
        "E_center_waist_absx_le_a": center_waist,
        "near_wall_s_ge_0p9": near_wall,
    }

    for band in FOCAL_BANDS:
        mask = (np.abs(X) <= c) & (np.abs(RHO - R) <= band * R)
        ef = float(np.sum(amp2 * W * mask) * dx * ds / Etot)
        vf = float(np.sum(W * mask) * dx * ds / max(Vtot, 1e-300))
        key = f"{int(round(1000 * band)):03d}"
        row[f"eta_{key}R"] = ef
        row[f"volfrac_{key}R"] = vf
        row[f"enrichment_{key}R"] = _safe_ratio(ef, vf)

    # Joint diagnostics.  Large values are not proofs; they rank candidates.
    row["focus_not_wall_ratio_eta100_over_wall"] = _safe_ratio(row["eta_100R"], near_wall)
    row["hyperbolic_signature_candidate_score"] = (
        row["eta_100R"]
        * math.log1p(min(row["CF_axial_inter_over_horns"], 1e6))
        * math.log1p(min(row["focal_plane_drop_ratio"], 1e6))
    )
    return row


def axial_density_from_vec(vec_free: np.ndarray, prob: dict) -> Tuple[np.ndarray, np.ndarray]:
    Nx = prob["Nx"]; Ns = prob["Ns"]; free = prob["free"]
    xs = prob["xs"]; S = prob["S"]; rw = prob["rw"]
    U = np.zeros(Nx * Ns); U[free] = vec_free; U = U.reshape((Nx, Ns))
    W = S * rw * rw
    ds = (1.0 - S_MIN) / (Ns - 1)
    dens = np.sum(U * U * W, axis=1) * ds
    area = np.trapz(dens, xs)
    if area > 0:
        dens = dens / area
    return xs, dens


def dedup_modes(rows: List[dict], tol: float = DEDUP_KR_TOL) -> List[dict]:
    rows = sorted(rows, key=lambda r: (r["kR"], -r.get("hyperbolic_signature_candidate_score", 0.0)))
    out: List[dict] = []
    for r in rows:
        if not out or abs(r["kR"] - out[-1]["kR"]) > tol:
            out.append(r)
        else:
            # Keep the better candidate by joint signature; tie-break by eta_100R.
            s1 = r.get("hyperbolic_signature_candidate_score", 0.0)
            s0 = out[-1].get("hyperbolic_signature_candidate_score", 0.0)
            if (s1, r.get("eta_100R", 0.0)) > (s0, out[-1].get("eta_100R", 0.0)):
                out[-1] = r
    return out

# -----------------------------------------------------------------------------
# Solver wrappers
# -----------------------------------------------------------------------------

def solve_shift_scan(
    Nx: int,
    Ns: int,
    a: float,
    b: float,
    R: float,
    kind: str,
    bc: str,
    m: int,
    kR_targets: np.ndarray,
    nev: int = 6,
    tol: float = 1e-6,
    maxiter: int = 2600,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    prob = build_problem(Nx, Ns, m, kind, bc, a, b, R)
    A = prob["A"]; Sscale = prob["Sscale"]
    all_rows: List[dict] = []
    envelope_rows: List[dict] = []
    k_eigs = min(nev, max(1, A.shape[0] - 2))

    for target in kR_targets:
        sigma = (float(target) / R) ** 2
        try:
            vals, y = spla.eigsh(A, k=k_eigs, sigma=sigma, which="LM", tol=tol, maxiter=maxiter)
            vecs = Sscale @ y
            order = np.argsort(vals)
            vals = vals[order]; vecs = vecs[:, order]
            local_rows = []
            for j, val in enumerate(vals):
                if val <= 1e-9:
                    continue
                row = mode_metrics(vecs[:, j], float(val), prob, target_kR=float(target))
                if row is not None:
                    local_rows.append(row)
                    all_rows.append(row)
            if local_rows:
                # Envelope records the best joint signature found near each shift.
                best = max(local_rows, key=lambda r: r["hyperbolic_signature_candidate_score"])
                er = dict(best)
                er["scan_target_kR"] = float(target)
                er["n_modes_returned"] = len(local_rows)
                envelope_rows.append(er)
            else:
                envelope_rows.append({"geometry": kind, "bc": bc, "a": a, "b": b, "R": R, "m": m, "scan_target_kR": float(target), "error": "no valid modes"})
        except Exception as e:
            envelope_rows.append({"geometry": kind, "bc": bc, "a": a, "b": b, "R": R, "m": m, "scan_target_kR": float(target), "error": str(e)[:240]})
            continue

    distinct = dedup_modes(all_rows)
    return pd.DataFrame(envelope_rows), pd.DataFrame(distinct)


def solve_single_target_with_density(
    Nx: int, Ns: int, a: float, b: float, R: float, kind: str, bc: str, m: int, kR_target: float, nev: int = 8
) -> Optional[Tuple[dict, np.ndarray, np.ndarray, dict]]:
    prob = build_problem(Nx, Ns, m, kind, bc, a, b, R)
    A = prob["A"]; Sscale = prob["Sscale"]
    sigma = (float(kR_target) / R) ** 2
    k_eigs = min(nev, max(1, A.shape[0] - 2))
    vals, y = spla.eigsh(A, k=k_eigs, sigma=sigma, which="LM", tol=1e-6, maxiter=3000)
    vecs = Sscale @ y
    order = np.argsort(vals); vals = vals[order]; vecs = vecs[:, order]
    best = None
    best_vec = None
    best_val = None
    for j, val in enumerate(vals):
        if val <= 1e-9:
            continue
        row = mode_metrics(vecs[:, j], float(val), prob, target_kR=kR_target)
        if row is None:
            continue
        if best is None or row["hyperbolic_signature_candidate_score"] > best["hyperbolic_signature_candidate_score"]:
            best = row; best_vec = vecs[:, j]; best_val = val
    if best is None:
        return None
    xs, dens = axial_density_from_vec(best_vec, prob)
    return best, xs, dens, prob

# -----------------------------------------------------------------------------
# Analysis: clusters and matched comparisons
# -----------------------------------------------------------------------------

def analyze_joint_clusters(df: pd.DataFrame, eta_thr: float, cf_thr: float, max_gap: float = 2.0) -> pd.DataFrame:
    """Find clusters of distinct modes satisfying eta and axial confinement thresholds."""
    if df.empty:
        return pd.DataFrame()
    rows = []
    for (geom, bc, m), g in df.groupby(["geometry", "bc", "m"]):
        d = g.sort_values("kR").reset_index(drop=True)
        ok = (d["eta_100R"] >= eta_thr) & (d["CF_axial_inter_over_horns"] >= cf_thr)
        idxs = list(np.where(ok.values)[0])
        if not idxs:
            continue
        clusters = []
        cur = [idxs[0]]
        for idx in idxs[1:]:
            if float(d.loc[idx, "kR"] - d.loc[cur[-1], "kR"]) <= max_gap:
                cur.append(idx)
            else:
                clusters.append(cur); cur = [idx]
        clusters.append(cur)
        for ci, cl in enumerate(clusters, 1):
            sub = d.loc[cl]
            width = float(sub["kR"].max() - sub["kR"].min()) if len(sub) > 1 else 0.0
            center = float(0.5 * (sub["kR"].max() + sub["kR"].min()))
            rows.append({
                "geometry": geom, "bc": bc, "m": int(m), "cluster_id": ci,
                "n_modes": int(len(sub)),
                "kR_min": float(sub["kR"].min()), "kR_max": float(sub["kR"].max()),
                "delta_kR": width,
                "Q_cluster_proxy": float(center / width) if width > 1e-12 else np.inf,
                "mean_eta_100R": float(sub["eta_100R"].mean()),
                "max_eta_100R": float(sub["eta_100R"].max()),
                "mean_CF_axial": float(sub["CF_axial_inter_over_horns"].mean()),
                "max_CF_axial": float(sub["CF_axial_inter_over_horns"].max()),
                "mean_drop_ratio": float(sub["focal_plane_drop_ratio"].mean()),
                "max_signature_score": float(sub["hyperbolic_signature_candidate_score"].max()),
                "kR_list": ";".join(f"{x:.4f}" for x in sub["kR"].tolist()),
            })
    return pd.DataFrame(rows).sort_values(["max_signature_score", "n_modes"], ascending=[False, False]) if rows else pd.DataFrame()


def matched_ph_vs_controls(df: pd.DataFrame, match_tol: float) -> pd.DataFrame:
    """For every PH mode, find nearest control mode at same m and bc."""
    if df.empty:
        return pd.DataFrame()
    rows = []
    ph = df[df["geometry"] == "PH"].copy()
    controls = df[df["geometry"] != "PH"].copy()
    for _, r in ph.iterrows():
        cand = controls[(controls["bc"] == r["bc"]) & (controls["m"] == r["m"])]
        if cand.empty:
            continue
        cand = cand.copy()
        cand["abs_delta_kR"] = np.abs(cand["kR"] - r["kR"])
        # one row per nearest geometry kind
        for ctrl, g in cand.groupby("geometry"):
            g2 = g.sort_values("abs_delta_kR")
            c0 = g2.iloc[0]
            row = {
                "bc": r["bc"], "m": int(r["m"]),
                "PH_kR": float(r["kR"]), "control_geometry": ctrl,
                "CTRL_kR": float(c0["kR"]), "abs_delta_kR": float(c0["abs_delta_kR"]),
                "matched_within_tol": bool(c0["abs_delta_kR"] <= match_tol),
                "PH_eta100": float(r["eta_100R"]), "CTRL_eta100": float(c0["eta_100R"]),
                "delta_eta100_pp": 100.0 * float(r["eta_100R"] - c0["eta_100R"]),
                "PH_CF": float(r["CF_axial_inter_over_horns"]), "CTRL_CF": float(c0["CF_axial_inter_over_horns"]),
                "CF_ratio_PH_over_CTRL": _safe_ratio(float(r["CF_axial_inter_over_horns"]), float(c0["CF_axial_inter_over_horns"])),
                "PH_horn_leakage": float(r["horn_leakage"]), "CTRL_horn_leakage": float(c0["horn_leakage"]),
                "horn_leakage_ratio_CTRL_over_PH": _safe_ratio(float(c0["horn_leakage"]), float(r["horn_leakage"])),
                "PH_drop_ratio": float(r["focal_plane_drop_ratio"]), "CTRL_drop_ratio": float(c0["focal_plane_drop_ratio"]),
                "drop_ratio_PH_over_CTRL": _safe_ratio(float(r["focal_plane_drop_ratio"]), float(c0["focal_plane_drop_ratio"])),
                "PH_signature_score": float(r["hyperbolic_signature_candidate_score"]),
                "CTRL_signature_score": float(c0["hyperbolic_signature_candidate_score"]),
                "signature_ratio_PH_over_CTRL": _safe_ratio(float(r["hyperbolic_signature_candidate_score"]), float(c0["hyperbolic_signature_candidate_score"])),
            }
            rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["matched_within_tol", "signature_ratio_PH_over_CTRL", "CF_ratio_PH_over_CTRL"], ascending=[False, False, False])


def summarize_uniqueness(df_modes: pd.DataFrame, df_match: pd.DataFrame, eta_thr: float, cf_ratio_thr: float, drop_ratio_thr: float) -> pd.DataFrame:
    rows = []
    if df_modes.empty:
        return pd.DataFrame()
    for (bc, m), g in df_modes.groupby(["bc", "m"]):
        ph = g[g["geometry"] == "PH"]
        ctr = g[g["geometry"] != "PH"]
        best_ph = ph.sort_values("hyperbolic_signature_candidate_score", ascending=False).head(1)
        if best_ph.empty:
            continue
        bp = best_ph.iloc[0]
        # Matched comparisons for the best PH kR.
        subm = df_match[(df_match["bc"] == bc) & (df_match["m"] == m) & (np.abs(df_match["PH_kR"] - bp["kR"]) <= DEDUP_KR_TOL + 1e-9)]
        if subm.empty:
            # fallback nearest PH kR in match table
            subm = df_match[(df_match["bc"] == bc) & (df_match["m"] == m)].copy()
            if not subm.empty:
                subm["dPH"] = np.abs(subm["PH_kR"] - bp["kR"])
                subm = subm.sort_values("dPH").head(10)
        n_controls = int(subm["control_geometry"].nunique()) if not subm.empty else 0
        n_ph_wins_cf = int((subm["CF_ratio_PH_over_CTRL"] >= cf_ratio_thr).sum()) if not subm.empty else 0
        n_ph_wins_drop = int((subm["drop_ratio_PH_over_CTRL"] >= drop_ratio_thr).sum()) if not subm.empty else 0
        n_ph_wins_eta = int((subm["PH_eta100"] - subm["CTRL_eta100"] >= 0.05).sum()) if not subm.empty else 0
        rows.append({
            "bc": bc, "m": int(m),
            "best_PH_kR": float(bp["kR"]),
            "best_PH_eta100": float(bp["eta_100R"]),
            "best_PH_CF": float(bp["CF_axial_inter_over_horns"]),
            "best_PH_drop_ratio": float(bp["focal_plane_drop_ratio"]),
            "best_PH_horn_leakage": float(bp["horn_leakage"]),
            "best_PH_signature_score": float(bp["hyperbolic_signature_candidate_score"]),
            "n_controls_compared": n_controls,
            "PH_wins_CF_ratio_count": n_ph_wins_cf,
            "PH_wins_drop_ratio_count": n_ph_wins_drop,
            "PH_wins_eta_by_5pp_count": n_ph_wins_eta,
            "uniqueness_status_hint": (
                "promising_hyperbolic_axial_signature" if (bp["eta_100R"] >= eta_thr and n_controls > 0 and n_ph_wins_cf >= max(1, n_controls//2) and n_ph_wins_drop >= max(1, n_controls//2))
                else "not_yet_unique_or_needs_other_metric"
            ),
        })
    return pd.DataFrame(rows).sort_values("best_PH_signature_score", ascending=False)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

def plot_geometry_profiles(a: float, b: float, R: float, kinds: List[str], outfile: Path) -> None:
    L = common_length(a, b, R)
    x = np.linspace(-L, L, 1100)
    c = math.sqrt(a * a + b * b)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for kind in kinds:
        r = wall_profile(x, kind, a, b, R, L)
        ax.plot(x, r, label=geometry_labels().get(kind, kind))
        ax.plot(x, -r, color=ax.lines[-1].get_color(), lw=0.8)
    ax.axvline(-c, color="k", ls=":", lw=1, label="focal planes x=±c")
    ax.axvline(c, color="k", ls=":", lw=1)
    ax.axhline(R, color="tab:green", ls="--", lw=1, label="rho=R")
    ax.axhline(-R, color="tab:green", ls="--", lw=1)
    ax.fill_between([-c, c], [R - 0.1 * R, R - 0.1 * R], [R + 0.1 * R, R + 0.1 * R], alpha=0.13, label="diagnostic ±10%R")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("signed rho")
    ax.set_title(f"Geometry controls for uniqueness test: a={a}, b={b}, R={R}")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout(); fig.savefig(outfile, dpi=180); plt.close(fig)


def plot_metric_spectra(df: pd.DataFrame, outdir: Path, bc_filter: str = "natural_TE_like_scalar_proxy") -> None:
    if df.empty:
        return
    d = df[df["bc"] == bc_filter].copy()
    if d.empty:
        d = df.copy()
    for m in sorted(d["m"].unique()):
        gm = d[d["m"] == m]
        if gm.empty:
            continue
        fig, axes = plt.subplots(3, 1, figsize=(9.5, 10.2), sharex=True)
        for kind, g in gm.groupby("geometry"):
            g = g.sort_values("kR")
            label = kind
            axes[0].plot(g["kR"], 100 * g["eta_100R"], marker="o", ms=3, lw=1, label=label)
            axes[1].semilogy(g["kR"], g["CF_axial_inter_over_horns"], marker="o", ms=3, lw=1, label=label)
            axes[2].semilogy(g["kR"], g["focal_plane_drop_ratio"], marker="o", ms=3, lw=1, label=label)
        axes[0].axhline(70, color="red", ls="--", lw=1, label="C2 eta threshold")
        axes[0].set_ylabel("eta ±10%R, %")
        axes[1].set_ylabel("CF = E(|x|≤c)/E(|x|>c)")
        axes[2].set_ylabel("focal-plane drop ratio")
        axes[2].set_xlabel("distinct eigenmode kR")
        axes[0].set_title(f"Hyperbolic uniqueness search metrics, m={m}, {bc_filter}")
        for ax in axes:
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7, ncol=2)
        fig.tight_layout(); fig.savefig(outdir / f"spectra_eta_CF_drop_m{int(m)}.png", dpi=180); plt.close(fig)


def plot_match_bars(df_match: pd.DataFrame, outdir: Path) -> None:
    if df_match.empty:
        return
    # Top matched rows where control is close in kR.
    d = df_match[df_match["matched_within_tol"]].copy()
    if d.empty:
        d = df_match.copy()
    d = d.sort_values("signature_ratio_PH_over_CTRL", ascending=False).head(20)
    if d.empty:
        return
    d["label"] = d.apply(lambda r: f"m={int(r.m)}, {r.control_geometry}, ΔkR={r.abs_delta_kR:.2f}", axis=1)
    fig, ax = plt.subplots(figsize=(10.5, 7))
    y = np.arange(len(d))
    ax.barh(y, d["signature_ratio_PH_over_CTRL"])
    ax.axvline(1.0, color="k", lw=1)
    ax.axvline(2.0, color="red", ls="--", lw=1, label="2× heuristic")
    ax.set_yticks(y); ax.set_yticklabels(d["label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("PHB / control joint signature ratio")
    ax.set_title("Matched-mode hyperbolic uniqueness score")
    ax.grid(axis="x", alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "matched_mode_signature_ratios.png", dpi=180); plt.close(fig)


def plot_axial_density_examples(args, df_modes: pd.DataFrame, df_match: pd.DataFrame, outdir: Path) -> None:
    if df_modes.empty:
        return
    ph = df_modes[df_modes["geometry"] == "PH"].copy()
    if ph.empty:
        return
    # choose the strongest PH signature among natural boundary if possible
    ph_nat = ph[ph["bc"] == "natural_TE_like_scalar_proxy"]
    if not ph_nat.empty:
        ph = ph_nat
    best_ph = ph.sort_values("hyperbolic_signature_candidate_score", ascending=False).iloc[0]
    m = int(best_ph["m"]); bc = str(best_ph["bc"]); kR = float(best_ph["kR"])
    examples = [("PH", kR)]
    if not df_match.empty:
        sub = df_match[(df_match["bc"] == bc) & (df_match["m"] == m)]
        sub = sub.iloc[(sub["PH_kR"] - kR).abs().argsort()].head(20) if not sub.empty else sub
        # choose LINEAR if available, otherwise best nearest control
        lin = sub[sub["control_geometry"] == "LINEAR"]
        if not lin.empty:
            ctrl = lin.sort_values("abs_delta_kR").iloc[0]
            examples.append(("LINEAR", float(ctrl["CTRL_kR"])))
        elif not sub.empty:
            ctrl = sub.sort_values("abs_delta_kR").iloc[0]
            examples.append((str(ctrl["control_geometry"]), float(ctrl["CTRL_kR"])))
    fig, ax = plt.subplots(figsize=(9, 5.4))
    c = math.sqrt(args.a * args.a + args.b * args.b)
    for kind, kk in examples:
        try:
            res = solve_single_target_with_density(args.Nx, args.Ns, args.a, args.b, args.R, kind, bc, m, kk, nev=max(args.nev, 8))
            if res is None:
                continue
            row, xs, dens, prob = res
            ax.plot(xs, dens, lw=1.5, label=f"{kind}: kR={row['kR']:.3f}, CF={row['CF_axial_inter_over_horns']:.2g}")
            # Also field plot for this mode
            save_field_plot(row, prob, outdir / f"field_{kind}_m{m}_kR{row['kR']:.3f}.png")
        except Exception as e:
            print(f"Warning: failed axial density example for {kind}: {e}")
    ax.axvline(-c, color="k", ls=":", lw=1, label="x=±c")
    ax.axvline(c, color="k", ls=":", lw=1)
    ax.set_xlabel("x")
    ax.set_ylabel("normalized axial energy density")
    ax.set_title("Axial energy-density profile: does PHB suppress horn leakage beyond |x|=c?")
    ax.grid(alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "axial_density_PH_vs_matched_control.png", dpi=180); plt.close(fig)


def save_field_plot(row: dict, prob: dict, outfile: Path) -> None:
    # Re-solve at row['kR'] and plot best field by joint signature for robustness.
    kind = row["geometry"]; bc = row["bc"]; m = int(row["m"]); kR = float(row["kR"])
    a = float(row["a"]); b = float(row["b"]); R = float(row["R"])
    res = solve_single_target_with_density(prob["Nx"], prob["Ns"], a, b, R, kind, bc, m, kR, nev=8)
    if res is None:
        return
    rrow, xs, dens, prob2 = res
    # Need actual field: re-run quickly inside solve_single_target? We did not return U. Simpler: plot density map from new solve is not available.
    # For a 2D visual, plot geometry + axial density as a pseudo-heat strip. This keeps the script reproducible without storing eigenvectors.
    X = prob2["X"]; S = prob2["S"]; rw = prob2["rw"]
    RHO = S * rw
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    # Pseudo-density repeated over rho for visual context; not used as scientific field map.
    D = np.tile(dens[:, None], (1, prob2["Ns"]))
    im = ax.pcolormesh(X, RHO, np.log10(D + 1e-12), shading="auto", cmap="magma")
    ax.plot(xs, rw[:, 0], color="white", lw=1.0)
    c = math.sqrt(a * a + b * b)
    ax.axvline(-c, color="cyan", ls=":", lw=1)
    ax.axvline(c, color="cyan", ls=":", lw=1)
    ax.axhline(R, color="lime", ls="--", lw=1)
    ax.axhline(R - 0.1 * R, color="lime", ls=":", lw=0.8)
    ax.axhline(R + 0.1 * R, color="lime", ls=":", lw=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x"); ax.set_ylabel("rho")
    ax.set_title(f"Axial-density visual context: {kind}, m={m}, kR≈{rrow['kR']:.3f}")
    fig.colorbar(im, ax=ax, label="log10 normalized axial density")
    fig.tight_layout(); fig.savefig(outfile, dpi=170); plt.close(fig)

# -----------------------------------------------------------------------------
# Workflow
# -----------------------------------------------------------------------------

def parse_list(s: str, cast=str) -> List:
    if isinstance(s, list):
        return s
    if not s:
        return []
    return [cast(x.strip()) for x in str(s).split(",") if x.strip()]


def apply_preset(args) -> None:
    if args.preset == "quick":
        args.Nx = 34; args.Ns = 24; args.nev = 4
        args.m_list = "15"
        args.kR_min = 14.0; args.kR_max = 18.0; args.kR_step = 1.0
        args.controls = "PH,LINEAR"
    elif args.preset == "target":
        args.Nx = 70; args.Ns = 48; args.nev = 6
        args.m_list = "15"
        args.kR_min = 10.0; args.kR_max = 25.0; args.kR_step = 0.5
        args.controls = "PH,LINEAR,POLY2,POLY3,SMOOTHSTEP,CIRCULAR"
    elif args.preset == "broad":
        args.Nx = 60; args.Ns = 42; args.nev = 5
        args.m_list = "6,10,12,15,20,25"
        args.kR_min = 10.0; args.kR_max = 28.0; args.kR_step = 1.0
        args.controls = "PH,LINEAR,POLY1p5,POLY2,POLY3,SMOOTHSTEP,CIRCULAR,COSINE"
    elif args.preset == "full":
        args.Nx = 70; args.Ns = 48; args.nev = 6
        args.m_list = "0,3,6,10,12,15,20,25"
        args.kR_min = 8.0; args.kR_max = 30.0; args.kR_step = 0.5
        args.controls = "PH,LINEAR,POLY0p5,POLY1p5,POLY2,POLY3,SMOOTHSTEP,CIRCULAR,COSINE"


def run(args) -> Path:
    apply_preset(args)
    outdir = Path(args.outdir)
    resdir = outdir / "results"
    figdir = outdir / "figures"
    resdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    kinds = parse_list(args.controls, str)
    m_values = parse_list(args.m_list, int)
    bcs = parse_list(args.bc_list, str)
    kR_targets = np.arange(args.kR_min, args.kR_max + 1e-12, args.kR_step)

    if "PH" not in kinds:
        kinds = ["PH"] + kinds

    meta = {
        "script": "PHB_C2C3_hyperbolic_uniqueness_search.py",
        "purpose": "Search for wave-level hyperbolic uniqueness beyond simple annular eta.",
        "model_status": "reduced axisymmetric scalar Helmholtz eigenproblem; not full-vector Maxwell; not 2D vector TE; not laser output",
        "target_geometry": {"a": args.a, "b": args.b, "R": args.R},
        "grid": {"Nx": args.Nx, "Ns": args.Ns},
        "m_values": m_values,
        "kR_targets": kR_targets.tolist(),
        "geometries": kinds,
        "boundary_proxies": bcs,
        "C2_eta_threshold": args.eta_threshold,
        "CF_threshold_for_cluster": args.cf_threshold,
        "match_tol_kR": args.match_tol_kR,
        "core_hypothesis": "If the hyperbolic focal billiard leaves a scalar wave imprint, PHB should show stronger inter-focal axial confinement, lower horn leakage, and a stronger focal-plane energy drop than controls at matched m and kR.",
    }
    (resdir / "run_metadata.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    plot_geometry_profiles(args.a, args.b, args.R, kinds, figdir / "geometry_profiles_controls.png")

    all_env = []
    all_distinct = []

    total_jobs = len(kinds) * len(bcs) * len(m_values)
    job = 0
    print("=" * 80)
    print("PHB hyperbolic uniqueness search")
    print("Model: reduced scalar axisymmetric Helmholtz eigenproblem")
    print(f"Target: a={args.a}, b={args.b}, R={args.R}; grid={args.Nx}x{args.Ns}")
    print(f"Geometries: {kinds}")
    print(f"m values: {m_values}")
    print(f"kR scan: {args.kR_min}..{args.kR_max}, step={args.kR_step}")
    print("=" * 80)

    for kind in kinds:
        for bc in bcs:
            for m in m_values:
                job += 1
                print(f"[{job}/{total_jobs}] geometry={kind}, bc={bc}, m={m}")
                try:
                    env, dist = solve_shift_scan(args.Nx, args.Ns, args.a, args.b, args.R, kind, bc, m, kR_targets, nev=args.nev, tol=args.eig_tol, maxiter=args.maxiter)
                    if not env.empty:
                        all_env.append(env)
                    if not dist.empty:
                        all_distinct.append(dist)
                        dist.to_csv(resdir / f"distinct_modes_{kind}_{bc}_m{m}.csv", index=False)
                except Exception as e:
                    print(f"  ERROR in {kind} {bc} m={m}: {e}")
                    err = pd.DataFrame([{"geometry": kind, "bc": bc, "m": m, "error": str(e)}])
                    all_env.append(err)

    df_env = pd.concat(all_env, ignore_index=True) if all_env else pd.DataFrame()
    df_modes = pd.concat(all_distinct, ignore_index=True) if all_distinct else pd.DataFrame()

    df_env.to_csv(resdir / "shift_window_envelope_all.csv", index=False)
    df_modes.to_csv(resdir / "distinct_modes_all_geometries.csv", index=False)

    df_clusters = analyze_joint_clusters(df_modes, args.eta_threshold, args.cf_threshold, max_gap=args.cluster_max_gap_kR)
    df_clusters.to_csv(resdir / "joint_eta_CF_modal_clusters.csv", index=False)

    df_match = matched_ph_vs_controls(df_modes, args.match_tol_kR)
    df_match.to_csv(resdir / "matched_PH_vs_controls_by_kR.csv", index=False)

    df_summary = summarize_uniqueness(df_modes, df_match, args.eta_threshold, args.cf_ratio_threshold, args.drop_ratio_threshold)
    df_summary.to_csv(resdir / "hyperbolic_uniqueness_summary.csv", index=False)

    # Plots
    plot_metric_spectra(df_modes, figdir, bc_filter="natural_TE_like_scalar_proxy")
    plot_match_bars(df_match, figdir)
    plot_axial_density_examples(args, df_modes, df_match, figdir)

    # Create top tables for quick inspection.
    if not df_modes.empty:
        top_sig = df_modes.sort_values("hyperbolic_signature_candidate_score", ascending=False).head(40)
        top_sig.to_csv(resdir / "top40_by_hyperbolic_signature_score.csv", index=False)
        top_cf = df_modes.sort_values("CF_axial_inter_over_horns", ascending=False).head(40)
        top_cf.to_csv(resdir / "top40_by_axial_confinement_CF.csv", index=False)
        top_eta = df_modes.sort_values("eta_100R", ascending=False).head(40)
        top_eta.to_csv(resdir / "top40_by_eta100.csv", index=False)

    # Human-readable summary.
    lines = []
    lines.append("PHB C2-C3 hyperbolic uniqueness search summary")
    lines.append("=" * 78)
    lines.append("Question: does the hyperbolic generatrix leave a wave signature beyond generic high-m annular localization?")
    lines.append("")
    lines.append("New diagnostic idea:")
    lines.append("  The earlier eta±10%R metric can be dominated by high-m annular/toroidal topology. The new test asks whether PHB also suppresses energy leakage beyond the external-focal planes |x|=c, as expected from the hyperbolic focal billiard mechanism.")
    lines.append("")
    lines.append("Metrics added:")
    lines.append("  CF_axial = E(|x|<=c) / E(|x|>c). Large CF means inter-focal axial confinement.")
    lines.append("  horn_leakage = E(|x|>c). Small leakage means less field penetration into the horns beyond the focal planes.")
    lines.append("  focal_plane_drop_ratio = E(c-strip <= |x| <= c) / E(c < |x| <= c+strip). Large value means a sharper energy barrier at x=±c.")
    lines.append("  matched PH/control ratios compare the same m and boundary proxy at nearest kR.")
    lines.append("")
    if not df_modes.empty:
        ph = df_modes[df_modes.geometry == "PH"].copy()
        best_ph_sig = ph.sort_values("hyperbolic_signature_candidate_score", ascending=False).head(1)
        best_ph_eta = ph.sort_values("eta_100R", ascending=False).head(1)
        if not best_ph_sig.empty:
            r = best_ph_sig.iloc[0]
            lines.append("Best PHB candidate by joint hyperbolic-signature score:")
            lines.append(f"  m={int(r.m)}, bc={r.bc}, kR={r.kR:.4f}, eta±10%R={100*r.eta_100R:.2f}%, CF={r.CF_axial_inter_over_horns:.3g}, horn leakage={100*r.horn_leakage:.4f}%, drop ratio={r.focal_plane_drop_ratio:.3g}")
        if not best_ph_eta.empty:
            r = best_ph_eta.iloc[0]
            lines.append("Best PHB candidate by usual eta±10%R:")
            lines.append(f"  m={int(r.m)}, bc={r.bc}, kR={r.kR:.4f}, eta±10%R={100*r.eta_100R:.2f}%, CF={r.CF_axial_inter_over_horns:.3g}, horn leakage={100*r.horn_leakage:.4f}%")
    if not df_match.empty:
        close = df_match[df_match.matched_within_tol]
        lines.append("")
        lines.append("Matched PHB-vs-control summary:")
        lines.append(f"  matched rows within ΔkR <= {args.match_tol_kR}: {len(close)}")
        if len(close) > 0:
            lines.append(f"  median CF ratio PH/control: {close.CF_ratio_PH_over_CTRL.median():.3g}")
            lines.append(f"  median focal-plane drop ratio PH/control: {close.drop_ratio_PH_over_CTRL.median():.3g}")
            lines.append(f"  median eta advantage PH-control: {close.delta_eta100_pp.median():.3g} percentage points")
            strong = close[(close.CF_ratio_PH_over_CTRL >= args.cf_ratio_threshold) & (close.drop_ratio_PH_over_CTRL >= args.drop_ratio_threshold)]
            lines.append(f"  rows where PHB beats controls by CF>={args.cf_ratio_threshold} and drop>={args.drop_ratio_threshold}: {len(strong)}")
            lines.append("  Top matched rows by signature ratio:")
            cols = ["m", "bc", "PH_kR", "control_geometry", "CTRL_kR", "abs_delta_kR", "PH_eta100", "CTRL_eta100", "CF_ratio_PH_over_CTRL", "drop_ratio_PH_over_CTRL", "signature_ratio_PH_over_CTRL"]
            lines.append(close.sort_values("signature_ratio_PH_over_CTRL", ascending=False)[cols].head(12).to_string(index=False))
    if not df_clusters.empty:
        lines.append("")
        lines.append("Joint eta+CF modal clusters:")
        lines.append(df_clusters.head(12).to_string(index=False))
    else:
        lines.append("")
        lines.append("No joint eta+CF modal clusters were detected at the selected thresholds. This does not refute the idea; it means the thresholds or searched m/kR range may need adjustment.")
    if not df_summary.empty:
        lines.append("")
        lines.append("Uniqueness status hints by m and boundary proxy:")
        lines.append(df_summary.to_string(index=False))
    lines.append("")
    lines.append("Interpretation rules:")
    lines.append("  Positive uniqueness evidence would require PHB to beat several controls at matched kR by axial confinement CF, horn-leakage suppression, and focal-plane drop, not merely by eta±10%R.")
    lines.append("  If PHB and controls have similar CF and drop ratios, then the closed scalar wave result still supports generic high-m annular/toroidal localization rather than unique hyperbolic wave localization.")
    lines.append("  If PHB wins strongly in CF/drop while eta remains comparable, that would be exactly the kind of wave imprint expected from the hyperbolic focal billiard.")
    lines.append("")
    lines.append(f"Runtime: {time.time() - t0:.1f} s")

    summary = "\n".join(lines)
    (resdir / "RUN_SUMMARY_hyperbolic_uniqueness_search.txt").write_text(summary, encoding="utf-8")
    print(summary)

    # README for Zenodo/checking.
    readme = f"""# PHB C2-C3 hyperbolic uniqueness search

This package tests whether the closed second-order vertical pseudohyperboloid has a wave-level signature of the hyperbolic focal billiard that is not visible in the usual annular-energy metric alone.

## Model status
Reduced axisymmetric scalar Helmholtz eigenproblem. This is not full-vector Maxwell and not a 2D vector TE/TM calculation.

## Main command
```bash
python PHB_C2C3_hyperbolic_uniqueness_search.py --preset target
```

## Important outputs
- `results/distinct_modes_all_geometries.csv` — deduplicated eigenmodes and all metrics.
- `results/matched_PH_vs_controls_by_kR.csv` — matched PHB/control comparisons at similar kR.
- `results/joint_eta_CF_modal_clusters.csv` — modal clusters satisfying eta and axial-confinement thresholds.
- `results/hyperbolic_uniqueness_summary.csv` — compact status hints.
- `results/RUN_SUMMARY_hyperbolic_uniqueness_search.txt` — human-readable summary.

## Key metrics
- `eta_100R`: energy in |x|<=c and |rho-R|<=0.10R.
- `CF_axial_inter_over_horns`: E(|x|<=c)/E(|x|>c).
- `horn_leakage`: E(|x|>c).
- `focal_plane_drop_ratio`: energy just inside x=±c divided by energy just outside x=±c.

Positive hyperbolic specificity would require PHB to outperform controls in CF and drop ratio at matched m, boundary proxy, and kR, not merely in eta_100R.
"""
    (outdir / "README_REPRODUCIBILITY.md").write_text(readme, encoding="utf-8")

    # Zip package
    package = Path(args.package)
    if package:
        with zipfile.ZipFile(package, "w", zipfile.ZIP_DEFLATED) as z:
            z.write(Path(__file__), Path(__file__).name)
            for p in outdir.rglob("*"):
                if p.is_file():
                    z.write(p, p.relative_to(outdir.parent))
        print(f"ZIP package written: {package}")
    return outdir


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="target", choices=["quick", "target", "broad", "full", "custom"], help="Run preset. Use custom to keep explicit CLI values.")
    ap.add_argument("--a", type=float, default=0.3)
    ap.add_argument("--b", type=float, default=0.6)
    ap.add_argument("--R", type=float, default=3.0)
    ap.add_argument("--m-list", default="15")
    ap.add_argument("--controls", default="PH,LINEAR,POLY2,POLY3,SMOOTHSTEP,CIRCULAR")
    ap.add_argument("--bc-list", default="natural_TE_like_scalar_proxy,dirichlet_TM_like_scalar_proxy")
    ap.add_argument("--kR-min", type=float, default=10.0)
    ap.add_argument("--kR-max", type=float, default=25.0)
    ap.add_argument("--kR-step", type=float, default=0.5)
    ap.add_argument("--Nx", type=int, default=70)
    ap.add_argument("--Ns", type=int, default=48)
    ap.add_argument("--nev", type=int, default=6)
    ap.add_argument("--eig-tol", type=float, default=1e-6)
    ap.add_argument("--maxiter", type=int, default=2600)
    ap.add_argument("--eta-threshold", type=float, default=0.70)
    ap.add_argument("--cf-threshold", type=float, default=50.0, help="CF threshold used for joint eta+CF cluster search.")
    ap.add_argument("--cluster-max-gap-kR", type=float, default=2.0)
    ap.add_argument("--match-tol-kR", type=float, default=0.30)
    ap.add_argument("--cf-ratio-threshold", type=float, default=2.0)
    ap.add_argument("--drop-ratio-threshold", type=float, default=2.0)
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parents[1] / "run_core_search"))
    ap.add_argument("--package", default=str(Path(__file__).resolve().parents[1] / "PHB_C2C3_hyperbolic_uniqueness_package.zip"))
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    run(args)


if __name__ == "__main__":
    main()

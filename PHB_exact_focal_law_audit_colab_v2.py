# -*- coding: utf-8 -*-
"""
Exact PHB focal-law audit and focus-directed ray tracing
=======================================================

Standalone Google Colab script.

This script fixes the weak point of the previous ray-tracing scripts:
it does NOT rely on a polygonal tangent to test the focal law. It uses
the exact hyperbola equation and the exact analytic normal at the
reflection point.

What is checked
---------------
1. Local law of reflection:
       angle of incidence = angle of reflection

2. Exact hyperbolic focal property in every meridional section:
       a ray directed toward one external focus is reflected toward
       the opposite external focus, in the virtual-focal sense.

3. Multi-bounce focus-directed tracing:
       a central internal source launches rays strictly toward the
       external PHB foci F± = (±c, ±R).  Each reflection is audited:
       incoming focus -> opposite outgoing focus.

4. Linear control:
       the same focus-directed rays are traced against a straight
       generator with the same endpoints.  It is not expected to pass
       the hyperbolic focal-law audit.

Scientific restriction
----------------------
This script verifies the geometrical ray law. It does not prove Q,
laser oscillation, gain, wave localization, or far-field diffraction-
limited collimation.

Outputs
-------
PHB_exact_focal_law_audit_results/
    00_parameters.json
    01_PHB_local_focal_law_audit.csv
    02_PHB_multibounce_focus_trace.csv
    03_LINEAR_multibounce_focus_trace.csv
    04_summary_metrics.csv
    05_PHB_exact_focus_traces.png
    06_LINEAR_focus_traces.png
    07_focal_law_errors.png
    README_scientific_interpretation.txt
    PHB_exact_focal_law_audit_results.zip
"""

import math
import json
import zipfile
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# PARAMETERS
# =============================================================================

@dataclass
class Params:
    # Exact PHB geometry
    a: float = 0.5
    b: float = 0.5
    R: float = 2.0
    R1: float = 0.0
    R2: float = 0.0
    R3: float = 0.10
    R4: float = 0.0

    # Ray-optical truncated tip aperture for optional output detection.
    # This does not define the focal theorem. It only defines where a traced
    # ray is considered to have left through the axial tip.
    tip_aperture: float = None

    # Central source. For the exact central test keep both jitter values zero.
    z_source: float = 0.0
    rho_source: float = 0.0
    z_source_jitter_fraction_of_a: float = 0.0
    rho_source_jitter_fraction_of_R: float = 0.0

    # Ray ensemble: rays are divided among the four signed foci.
    n_rays: int = 160
    random_seed: int = 42

    # Reflection tracing
    max_bounces: int = 10
    root_samples_per_branch: int = 120
    bisection_iters: int = 55
    eps: float = 1e-10

    # Direct local audit of the exact focal property
    n_audit_points_per_branch: int = 250

    # Output
    outdir: str = "PHB_exact_focal_law_audit_results"


P = Params()
if P.tip_aperture is None:
    P.tip_aperture = P.R3 if P.R3 > 0 else 0.05 * P.R


# =============================================================================
# BASIC GEOMETRY
# =============================================================================

def L_half(p=P):
    return p.a * math.sqrt(1.0 + (p.R / p.b)**2)


def c_focus(p=P):
    return math.sqrt(p.a*p.a + p.b*p.b)


def z_cut_for_aperture(model, p=P):
    """Positive z where wall radius equals p.tip_aperture."""
    L = L_half(p)
    ap = max(min(float(p.tip_aperture), 0.98*p.R), 1e-12)
    if model == "phb":
        q = (p.R - ap) / p.b
        return p.a * math.sqrt(1.0 + q*q)
    if model == "linear":
        return L - ap * (L - p.a) / p.R
    raise ValueError(model)


def rho_wall_phb(z, p=P):
    z = np.asarray(z, dtype=float)
    az = np.abs(z)
    rho = np.empty_like(az)
    central = az <= p.a
    rho[central] = p.R
    horn = ~central
    val = (az[horn]/p.a)**2 - 1.0
    rho[horn] = p.R - p.b*np.sqrt(np.maximum(val, 0.0))
    return np.maximum(rho, 0.0)


def rho_wall_linear(z, p=P):
    z = np.asarray(z, dtype=float)
    L = L_half(p)
    az = np.abs(z)
    rho = np.empty_like(az)
    central = az <= p.a
    rho[central] = p.R
    horn = ~central
    rho[horn] = p.R * (L - az[horn]) / max(L - p.a, 1e-15)
    return np.maximum(rho, 0.0)


def rho_wall(model, z, p=P):
    if model == "phb":
        return rho_wall_phb(z, p)
    if model == "linear":
        return rho_wall_linear(z, p)
    raise ValueError(model)


def foci_for_side(side, p=P):
    """External foci for signed meridional branch side=+1 upper, -1 lower."""
    c = c_focus(p)
    return {
        "left": np.array([-c, side*p.R], dtype=float),
        "right": np.array([+c, side*p.R], dtype=float),
    }


def all_foci(p=P):
    c = c_focus(p)
    return {
        "right_upper": np.array([+c, +p.R], dtype=float),
        "left_upper":  np.array([-c, +p.R], dtype=float),
        "right_lower": np.array([+c, -p.R], dtype=float),
        "left_lower":  np.array([-c, -p.R], dtype=float),
    }


def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v/n if n > 0 else v*np.nan


def angle_deg(u, v):
    u = unit(u); v = unit(v)
    d = float(np.clip(np.dot(u, v), -1.0, 1.0))
    return math.degrees(math.acos(d))


def exact_phb_normal(z, rho, side, p=P):
    """
    Analytic normal of signed hyperbola:
        z^2/a^2 - (rho - side*R)^2/b^2 = 1
    """
    n = np.array([2.0*z/(p.a*p.a), -2.0*(rho - side*p.R)/(p.b*p.b)], dtype=float)
    return unit(n)


def exact_phb_tangent(z, rho, side, p=P):
    n = exact_phb_normal(z, rho, side, p)
    return np.array([-n[1], n[0]], dtype=float)


def reflect_by_normal(d_in, n):
    d = unit(d_in)
    n = unit(n)
    return unit(d - 2.0*np.dot(d, n)*n)


def reflected_direction(model, z, rho, side, d_in, p=P):
    if model == "phb":
        n = exact_phb_normal(z, rho, side, p)
        return reflect_by_normal(d_in, n), n
    if model == "linear":
        # Exact normal for the straight-generator control.
        # The horn line slope is constant on each side.
        L = L_half(p)
        az = abs(z)
        if az <= p.a:
            # central cylindrical segment: tangent along z, normal radial
            t = np.array([1.0, 0.0])
        else:
            # rho = R*(L-|z|)/(L-a); derivative depends on z sign and signed side.
            # signed boundary rho_signed = side*rho_wall.
            drdz_unsigned = -p.R/(L-p.a) * (1.0 if z > 0 else -1.0)
            drdz_signed = side * drdz_unsigned
            t = unit(np.array([1.0, drdz_signed]))
        n = unit(np.array([-t[1], t[0]]))
        return reflect_by_normal(d_in, n), n
    raise ValueError(model)


def branch_intervals(model, p=P):
    zc = z_cut_for_aperture(model, p)
    # Only active horns are used for the focal-law tracing.
    # The central gap is not treated as a hyperbolic reflector.
    return [(-zc, -p.a), (p.a, zc)]


def signed_wall_value(model, z, side, p=P):
    return side * float(rho_wall(model, np.array([z]), p)[0])


# =============================================================================
# EXACT RAY-CURVE INTERSECTION
# =============================================================================

def f_branch(model, p0, d, t, side, p=P):
    z = p0[0] + t*d[0]
    rho = p0[1] + t*d[1]
    return rho - signed_wall_value(model, z, side, p)


def t_interval_for_z_range(p0, d, zmin, zmax, p=P):
    dz = d[0]
    if abs(dz) < 1e-15:
        z0 = p0[0]
        if zmin <= z0 <= zmax:
            return (p.eps, 1e6)
        return None
    t1 = (zmin - p0[0]) / dz
    t2 = (zmax - p0[0]) / dz
    lo, hi = min(t1, t2), max(t1, t2)
    lo = max(lo, p.eps)
    if hi <= lo:
        return None
    return lo, hi


def bisect_root(model, p0, d, side, lo, hi, p=P):
    flo = f_branch(model, p0, d, lo, side, p)
    fhi = f_branch(model, p0, d, hi, side, p)
    if abs(flo) < 1e-12:
        return lo
    if abs(fhi) < 1e-12:
        return hi
    if flo*fhi > 0:
        return None
    a, b = lo, hi
    fa, fb = flo, fhi
    for _ in range(p.bisection_iters):
        m = 0.5*(a+b)
        fm = f_branch(model, p0, d, m, side, p)
        if abs(fm) < 1e-13 or (b-a) < 1e-13:
            return m
        if fa*fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5*(a+b)


def find_next_intersection(model, p0, d, p=P):
    """
    Find nearest intersection with active upper/lower horn branches.
    Returns dict or None.
    """
    candidates = []
    for zmin, zmax in branch_intervals(model, p):
        interval = t_interval_for_z_range(p0, d, zmin, zmax, p)
        if interval is None:
            continue
        tlo, thi = interval
        if not np.isfinite(thi) or thi <= tlo:
            continue

        for side in (+1, -1):
            # sample the valid t interval to find sign changes
            ts = np.linspace(tlo, thi, p.root_samples_per_branch)
            vals = np.array([f_branch(model, p0, d, t, side, p) for t in ts])
            # exact close hits
            close = np.where(np.abs(vals) < 1e-9)[0]
            if len(close):
                t = float(ts[close[0]])
                if t > p.eps:
                    candidates.append((t, side))
                continue
            # sign changes
            sgn = vals[:-1] * vals[1:]
            idxs = np.where(sgn <= 0)[0]
            for idx in idxs:
                lo, hi = float(ts[idx]), float(ts[idx+1])
                root = bisect_root(model, p0, d, side, lo, hi, p)
                if root is not None and root > p.eps:
                    candidates.append((root, side))
                    break

    if not candidates:
        return None

    t, side = min(candidates, key=lambda x: x[0])
    hit = p0 + t*d
    return {"t": float(t), "side": int(side), "hit": hit}


# =============================================================================
# FOCAL-LAW AUDITS
# =============================================================================

def local_focal_law_audit(p=P):
    """
    Direct audit independent of ray-tracing:
    Sample exact PHB boundary points and check:
        ray toward Fleft reflects toward Fright;
        ray toward Fright reflects toward Fleft.
    """
    rows = []
    zc = z_cut_for_aperture("phb", p)
    # avoid the cusp/endpoints for numerical stability
    for zmin, zmax, branch_name in [(-zc, -p.a, "left_branch"), (p.a, zc, "right_branch")]:
        zs = np.linspace(zmin + 1e-8, zmax - 1e-8, p.n_audit_points_per_branch)
        for side in (+1, -1):
            F = foci_for_side(side, p)
            for z in zs:
                rho = signed_wall_value("phb", z, side, p)
                Pnt = np.array([z, rho])
                n = exact_phb_normal(z, rho, side, p)
                for incoming_focus, outgoing_focus in [("left", "right"), ("right", "left")]:
                    d_in = unit(F[incoming_focus] - Pnt)
                    d_out = reflect_by_normal(d_in, n)
                    target_dir = unit(F[outgoing_focus] - Pnt)

                    theta_in_n = math.degrees(math.acos(abs(float(np.clip(np.dot(d_in, n), -1, 1)))))
                    theta_out_n = math.degrees(math.acos(abs(float(np.clip(np.dot(d_out, n), -1, 1)))))
                    rows.append({
                        "branch": branch_name,
                        "side": side,
                        "z": z,
                        "rho": rho,
                        "incoming_focus": incoming_focus,
                        "expected_outgoing_focus": outgoing_focus,
                        "angle_incidence_to_normal_deg": theta_in_n,
                        "angle_reflection_to_normal_deg": theta_out_n,
                        "angle_law_error_deg": abs(theta_in_n - theta_out_n),
                        "focus_switch_error_deg": angle_deg(d_out, target_dir),
                    })
    return pd.DataFrame(rows)


def generate_central_focus_rays(p=P):
    """
    Central source. Rays are aimed at the four external signed foci.
    With zero jitter, there are four unique meridional directions repeated.
    """
    rng = np.random.default_rng(p.random_seed)
    F = all_foci(p)
    keys = ["right_upper", "left_upper", "right_lower", "left_lower"]

    origins = []
    dirs = []
    target_keys = []

    for i in range(p.n_rays):
        z0 = p.z_source
        rho0 = p.rho_source
        if p.z_source_jitter_fraction_of_a > 0:
            z0 += rng.uniform(-p.z_source_jitter_fraction_of_a*p.a, p.z_source_jitter_fraction_of_a*p.a)
        if p.rho_source_jitter_fraction_of_R > 0:
            rho0 += rng.uniform(-p.rho_source_jitter_fraction_of_R*p.R, p.rho_source_jitter_fraction_of_R*p.R)

        key = keys[i % len(keys)]
        target = F[key]
        p0 = np.array([z0, rho0], dtype=float)
        d = unit(target - p0)
        origins.append(p0)
        dirs.append(d)
        target_keys.append(key)

    return np.asarray(origins), np.asarray(dirs), target_keys


def trace_focus_rays(model, p=P):
    origins, dirs, target_keys = generate_central_focus_rays(p)
    zc = z_cut_for_aperture(model, p)

    ray_rows = []
    bounce_rows = []
    paths = []

    for i in range(p.n_rays):
        pcur = origins[i].copy()
        dcur = dirs[i].copy()
        path = [pcur.copy()]
        status = "max_bounces"

        for bidx in range(p.max_bounces):
            # Optional axial tip output detection.
            if dcur[0] > 1e-14:
                t_exit = (zc - pcur[0]) / dcur[0]
                if t_exit > p.eps:
                    cand = pcur + t_exit*dcur
                    if abs(cand[1]) <= p.tip_aperture:
                        path.append(cand.copy())
                        status = "right_tip_output"
                        break
            if dcur[0] < -1e-14:
                t_exit = (-zc - pcur[0]) / dcur[0]
                if t_exit > p.eps:
                    cand = pcur + t_exit*dcur
                    if abs(cand[1]) <= p.tip_aperture:
                        path.append(cand.copy())
                        status = "left_tip_output"
                        break

            inter = find_next_intersection(model, pcur, dcur, p)
            if inter is None:
                # free escape outside active horn surfaces
                path.append(pcur + 0.5*p.R*dcur)
                status = "escaped_no_intersection"
                break

            hit = inter["hit"]
            side = inter["side"]
            z, rho = float(hit[0]), float(hit[1])

            d_before = dcur.copy()
            d_after, n = reflected_direction(model, z, rho, side, d_before, p)
            path.append(hit.copy())

            # Audit reflection law and focal switch for PHB/linear.
            theta_in_n = math.degrees(math.acos(abs(float(np.clip(np.dot(d_before, n), -1, 1)))))
            theta_out_n = math.degrees(math.acos(abs(float(np.clip(np.dot(d_after, n), -1, 1)))))

            F = foci_for_side(side, p)
            # Determine closest incoming focal direction and opposite.
            err_in_left = angle_deg(d_before, F["left"] - hit)
            err_in_right = angle_deg(d_before, F["right"] - hit)
            incoming_focus = "left" if err_in_left <= err_in_right else "right"
            expected = "right" if incoming_focus == "left" else "left"
            focus_switch_error = angle_deg(d_after, F[expected] - hit)

            bounce_rows.append({
                "ray_id": i,
                "model": model.upper(),
                "bounce": bidx + 1,
                "side": side,
                "z": z,
                "rho": rho,
                "incoming_focus_closest": incoming_focus,
                "expected_outgoing_focus": expected,
                "incoming_focus_error_deg": min(err_in_left, err_in_right),
                "angle_incidence_to_normal_deg": theta_in_n,
                "angle_reflection_to_normal_deg": theta_out_n,
                "angle_law_error_deg": abs(theta_in_n - theta_out_n),
                "focus_switch_error_deg": focus_switch_error,
            })

            dcur = d_after
            pcur = hit + 100.0*p.eps*dcur
        else:
            status = "max_bounces"

        paths.append(np.vstack(path))
        ray_rows.append({
            "ray_id": i,
            "model": model.upper(),
            "initial_target_focus": target_keys[i],
            "status": status,
            "n_path_points": len(path),
            "n_bounces_recorded": sum(1 for r in bounce_rows if r["ray_id"] == i and r["model"] == model.upper()),
            "final_z": float(path[-1][0]),
            "final_rho": float(path[-1][1]),
        })

    return pd.DataFrame(ray_rows), pd.DataFrame(bounce_rows), paths


def summarise_bounces(bdf, rdf, model):
    if len(bdf) == 0:
        return {
            "model": model,
            "n_rays": len(rdf),
            "n_bounces": 0,
            "max_angle_law_error_deg": np.nan,
            "mean_angle_law_error_deg": np.nan,
            "max_focus_switch_error_deg": np.nan,
            "mean_focus_switch_error_deg": np.nan,
            "median_focus_switch_error_deg": np.nan,
            "directed_output_fraction": float(np.mean(rdf["status"].isin(["right_tip_output", "left_tip_output"]))) if len(rdf) else np.nan,
        }
    return {
        "model": model,
        "n_rays": int(len(rdf)),
        "n_bounces": int(len(bdf)),
        "max_angle_law_error_deg": float(bdf["angle_law_error_deg"].max()),
        "mean_angle_law_error_deg": float(bdf["angle_law_error_deg"].mean()),
        "max_focus_switch_error_deg": float(bdf["focus_switch_error_deg"].max()),
        "mean_focus_switch_error_deg": float(bdf["focus_switch_error_deg"].mean()),
        "median_focus_switch_error_deg": float(bdf["focus_switch_error_deg"].median()),
        "directed_output_fraction": float(np.mean(rdf["status"].isin(["right_tip_output", "left_tip_output"]))),
    }


# =============================================================================
# PLOTS
# =============================================================================

def plot_geometry(ax, model, p=P):
    zc = z_cut_for_aperture(model, p)
    z_left = np.linspace(-zc, -p.a, 500)
    z_right = np.linspace(p.a, zc, 500)
    for zs in [z_left, z_right]:
        rw = rho_wall(model, zs, p)
        ax.plot(zs, +rw, color="black", lw=1.5)
        ax.plot(zs, -rw, color="black", lw=1.5)

    # central gap indicator
    ax.axvspan(-p.a, p.a, color="0.85", alpha=0.20, label="central gap")

    # foci
    for name, pt in all_foci(p).items():
        ax.plot(pt[0], pt[1], "o", ms=4, color="crimson")
    ax.text(+c_focus(p), +p.R*1.05, "F+", color="crimson", ha="center", fontsize=9)
    ax.text(-c_focus(p), +p.R*1.05, "F−", color="crimson", ha="center", fontsize=9)

    # source
    ax.plot(p.z_source, p.rho_source, marker="*", ms=10, color="tab:purple")
    return zc


def plot_traces(model, ray_df, paths, outpath, p=P, max_paths=120):
    fig, ax = plt.subplots(figsize=(9.5, 6.2), dpi=170)
    zc = plot_geometry(ax, model, p)
    ids = np.arange(len(paths))
    if len(ids) > max_paths:
        ids = ids[:max_paths]

    for i in ids:
        path = paths[i]
        status = ray_df.loc[i, "status"]
        col = "tab:blue" if status in ("right_tip_output", "left_tip_output") else "0.60"
        alpha = 0.55 if col == "tab:blue" else 0.25
        lw = 0.75 if col == "tab:blue" else 0.45
        ax.plot(path[:,0], path[:,1], color=col, alpha=alpha, lw=lw)

    ax.set_title(f"{model.upper()}: exact-normal focus-directed ray tracing")
    ax.set_xlabel("z")
    ax.set_ylabel("signed radius ρ")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-zc - 0.2*p.R, zc + 0.2*p.R)
    ax.set_ylim(-1.18*p.R, 1.18*p.R)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def plot_errors(local_df, phb_bounce, lin_bounce, outpath):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=170)

    axes[0].hist(np.log10(np.maximum(local_df["focus_switch_error_deg"], 1e-16)), bins=80, color="tab:blue", alpha=0.8)
    axes[0].set_title("Direct PHB focal-law audit")
    axes[0].set_xlabel("log10(focus-switch error, degrees)")
    axes[0].set_ylabel("count")
    axes[0].grid(True, alpha=0.25)

    data = []
    labels = []
    if len(phb_bounce):
        data.append(np.log10(np.maximum(phb_bounce["focus_switch_error_deg"], 1e-16)))
        labels.append("PHB")
    if len(lin_bounce):
        data.append(np.log10(np.maximum(lin_bounce["focus_switch_error_deg"], 1e-16)))
        labels.append("LINEAR")
    axes[1].boxplot(data, labels=labels, showfliers=False)
    axes[1].set_title("Multi-bounce focus-switch error")
    axes[1].set_ylabel("log10(error, degrees)")
    axes[1].grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    outdir = Path(P.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    (outdir / "00_parameters.json").write_text(json.dumps(asdict(P), indent=2), encoding="utf-8")

    print("=== Exact PHB focal-law audit ===")
    print(json.dumps(asdict(P), indent=2))
    print(f"L = {L_half(P):.10f}")
    print(f"c = {c_focus(P):.10f}")
    print("External foci in signed meridional section:")
    for k, v in all_foci(P).items():
        print(f"  {k}: z={v[0]:.10f}, rho={v[1]:.10f}")
    print()

    # Direct exact audit of the hyperbolic focal law.
    local_df = local_focal_law_audit(P)
    local_df.to_csv(outdir / "01_PHB_local_focal_law_audit.csv", index=False)

    # Multi-bounce tracing.
    phb_ray, phb_bounce, phb_paths = trace_focus_rays("phb", P)
    lin_ray, lin_bounce, lin_paths = trace_focus_rays("linear", P)

    phb_ray.to_csv(outdir / "02_PHB_ray_status.csv", index=False)
    phb_bounce.to_csv(outdir / "02_PHB_multibounce_focus_trace.csv", index=False)
    lin_ray.to_csv(outdir / "03_LINEAR_ray_status.csv", index=False)
    lin_bounce.to_csv(outdir / "03_LINEAR_multibounce_focus_trace.csv", index=False)

    # Summary.
    summary = pd.DataFrame([
        {
            "model": "PHB_DIRECT_LOCAL_THEOREM_AUDIT",
            "n_rays": np.nan,
            "n_bounces": len(local_df),
            "max_angle_law_error_deg": float(local_df["angle_law_error_deg"].max()),
            "mean_angle_law_error_deg": float(local_df["angle_law_error_deg"].mean()),
            "max_focus_switch_error_deg": float(local_df["focus_switch_error_deg"].max()),
            "mean_focus_switch_error_deg": float(local_df["focus_switch_error_deg"].mean()),
            "median_focus_switch_error_deg": float(local_df["focus_switch_error_deg"].median()),
            "directed_output_fraction": np.nan,
        },
        summarise_bounces(phb_bounce, phb_ray, "PHB_MULTIBOUNCE"),
        summarise_bounces(lin_bounce, lin_ray, "LINEAR_MULTIBOUNCE"),
    ])
    summary.to_csv(outdir / "04_summary_metrics.csv", index=False)
    print(summary.to_string(index=False))

    # Plots.
    plot_traces("phb", phb_ray, phb_paths, outdir / "05_PHB_exact_focus_traces.png", P)
    plot_traces("linear", lin_ray, lin_paths, outdir / "06_LINEAR_focus_traces.png", P)
    plot_errors(local_df, phb_bounce, lin_bounce, outdir / "07_focal_law_errors.png")

    readme = f"""Scientific interpretation
=========================

This run verifies the local geometrical focal law of the PHB meridional
hyperbola using the exact analytic normal.

The tested PHB equation in each signed meridional half-plane is:
    z^2/a^2 - (rho - sR)^2/b^2 = 1,  s = +1 upper, -1 lower.

For each sampled boundary point, an incident ray directed toward one external
focus is reflected by:
    d_out = d_in - 2 (d_in · n) n,

where n is the exact analytic normal. The script then checks:
    1. equality of incidence and reflection angles;
    2. angular mismatch between the reflected ray and the direction toward
       the opposite external focus.

The direct theorem audit is the cleanest result. The multi-bounce part is a
numerical ray-dynamics test and depends on truncation, apertures and source
definition.

Important wording for a paper:
    The PHB inherits the hyperbolic meridional virtual-focal reflection law.
    This is a local one-reflection theorem. Global multi-bounce concentration
    and narrow output formation are stronger dynamical hypotheses and must be
    reported as numerical tests, not as direct consequences of the theorem.

Parameters:
    a = {P.a}
    b = {P.b}
    R = {P.R}
    L = {L_half(P):.10f}
    c = {c_focus(P):.10f}
    tip_aperture = {P.tip_aperture}
"""
    (outdir / "README_scientific_interpretation.txt").write_text(readme, encoding="utf-8")

    # Zip archive.
    zippath = outdir / "PHB_exact_focal_law_audit_results.zip"
    with zipfile.ZipFile(zippath, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for f in outdir.iterdir():
            if f.is_file() and f.name != zippath.name:
                z.write(f, arcname=f.name)

    print()
    print(f"Saved folder: {outdir.resolve()}")
    print(f"Zip archive: {zippath.resolve()}")
    print()
    print("For Google Colab manual download:")
    print("from google.colab import files")
    print(f"files.download('{zippath.as_posix()}')")

    try:
        from google.colab import files  # type: ignore
        files.download(str(zippath))
    except Exception as exc:
        print(f"Auto-download skipped outside Colab: {exc}")


if __name__ == "__main__":
    main()

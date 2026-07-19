# -*- coding: utf-8 -*-
"""
Интерактивная визуализация вертикального псевдогиперболоида 2-го порядка
с выбором центральной заглушки:

1) open      — открытый PHB, без центральной заглушки;
2) halftorus — полутороидальная заглушка;
3) cylinder  — цилиндрическая заглушка.

Важно для текущей PHB/MEEP-семантики:
    R2 = ΔR — абсолютная радиальная ширина выходного окна,
    а не доля R и не R*R2.

R1 оставлен как относительный входной радиус: R_in = R * R1.

Для Jupyter/Colab: графики выводятся прямо в ноутбуке.

Изменения v2:
    - R1 и R2 по умолчанию начинаются с 0;
    - подписи/readout ползунков уменьшены и расширена сетка управления.
"""

import numpy as np
import matplotlib.pyplot as plt

import ipywidgets as widgets
from IPython.display import display, clear_output, HTML

try:
    get_ipython().run_line_magic("matplotlib", "inline")
except Exception:
    pass


# ============================================================
# СТИЛЬ
# ============================================================
LW_WALL_2D = 0.75
LW_CAP_2D = 0.85
LW_EDGE_2D = 0.60
LW_AUX_2D = 0.40
LW_ARROW_2D = 0.90

FOCUS_MARKER_SIZE = 3.0
FOCUS_FONT_SIZE = 8

CLOSURE_OPTIONS = [
    ("1. Открытый — без заглушки", "open"),
    ("2. Полутороидальная заглушка", "halftorus"),
    ("3. Цилиндрическая заглушка", "cylinder"),
]


# ============================================================
# ГЕОМЕТРИЯ
# ============================================================
def r_horn(x, a, b, R):
    """
    Гиперболическая образующая:
        rho(x) = R - b * sqrt((|x|/a)^2 - 1)

    Область:
        |x| >= a
    """
    x = np.asarray(x, dtype=float)
    val = (np.abs(x) / a) ** 2 - 1.0
    val = np.maximum(val, 0.0)
    return R - b * np.sqrt(val)


def half_torus_outer_rho(x, a, R):
    """
    Наружная полутороидальная заглушка.
    В 2D это две наружные полуокружности:

        верхняя: rho =  R + sqrt(a^2 - x^2)
        нижняя: rho = -R - sqrt(a^2 - x^2)

    при |x| <= a.
    """
    x = np.asarray(x, dtype=float)
    val = a ** 2 - x ** 2
    val = np.maximum(val, 0.0)
    return R + np.sqrt(val)


def plot_hyperbola_foci(ax, a, b, R):
    """
    Фокусы образующих гипербол.
    Верхняя гипербола:
        (x/a)^2 - ((rho - R)/b)^2 = 1
    Нижняя гипербола:
        (x/a)^2 - ((rho + R)/b)^2 = 1
    """
    c = np.sqrt(a ** 2 + b ** 2)

    foci = [
        (-c,  R, "F1+"),
        ( c,  R, "F2+"),
        (-c, -R, "F1-"),
        ( c, -R, "F2-"),
    ]

    for xf, yf, label in foci:
        ax.plot(xf, yf, "mo", ms=FOCUS_MARKER_SIZE, zorder=10)
        dy = 0.13 * max(1.0, R / 2.0) if yf > 0 else -0.13 * max(1.0, R / 2.0)
        va = "bottom" if yf > 0 else "top"
        ax.text(
            xf,
            yf + dy,
            label,
            color="m",
            fontsize=FOCUS_FONT_SIZE,
            ha="center",
            va=va,
            fontweight="bold",
            zorder=11,
        )


def compute_geometry(a, b, R, R1, R2_abs, closure):
    """
    R1 — доля от R: R_in = R * R1.
    R2_abs — абсолютная радиальная ширина выходного окна ΔR.

    Для визуализации выходное окно строится симметрично относительно
    линии фокального кольца rho=R:
        dR_total = R2_abs
        dR_half  = R2_abs / 2

    Если dR_half выходит за область отображаемой геометрии, он временно
    обрезается только для построения, чтобы график не падал.
    """
    R_in = R * R1
    dR_total = max(float(R2_abs), 0.0)
    dR_half_raw = 0.5 * dR_total

    L = a * np.sqrt(1.0 + (R / b) ** 2)
    c = np.sqrt(a ** 2 + b ** 2)

    show_inlet = R_in > 0.0
    show_outlet = dR_total > 0.0

    # Визуальное ограничение: нельзя вырезать больше доступного радиуса horn.
    max_half = 0.999 * R
    if closure == "halftorus":
        max_half = min(max_half, 0.999 * a)

    clipped = False
    if show_outlet and dR_half_raw > max_half:
        dR_half = max_half
        clipped = True
    else:
        dR_half = dR_half_raw

    if show_inlet:
        x_edge_l = -a * np.sqrt(1.0 + ((R - R_in) / b) ** 2)
    else:
        x_edge_l = -L

    if show_outlet and dR_half > 0.0:
        rho_horn_cut = R - dR_half
        rho_cap_cut = R + dR_half

        # r_horn(x_horn_out) = R - dR_half
        x_horn_out = a * np.sqrt(1.0 + (dR_half / b) ** 2)

        # Для полутороида: R + sqrt(a^2 - x_torus_out^2) = R + dR_half
        x_torus_out = np.sqrt(max(a ** 2 - dR_half ** 2, 0.0))

        # Для цилиндра это только условная точка визуального разрыва справа.
        x_cylinder_out = a
    else:
        rho_horn_cut = R
        rho_cap_cut = R
        x_horn_out = a
        x_torus_out = a
        x_cylinder_out = a

    return {
        "closure": closure,
        "R_in": R_in,
        "dR_total": dR_total,
        "dR_half_raw": dR_half_raw,
        "dR_half": dR_half,
        "clipped": clipped,
        "L": L,
        "c": c,
        "show_inlet": show_inlet,
        "show_outlet": show_outlet,
        "x_edge_l": x_edge_l,
        "rho_horn_cut": rho_horn_cut,
        "rho_cap_cut": rho_cap_cut,
        "x_horn_out": x_horn_out,
        "x_torus_out": x_torus_out,
        "x_cylinder_out": x_cylinder_out,
    }


# ============================================================
# 2D
# ============================================================
def draw_open_window_arrows_2d(ax, a, R, L):
    """Показывает центральный открытый радиальный выход для open PHB."""
    scale = max(1.0, a, R)
    for x0 in [-0.55 * a, 0.0, 0.55 * a]:
        ax.annotate(
            "",
            xy=(x0, R + 0.70 * scale),
            xytext=(x0, R + 0.05 * scale),
            arrowprops=dict(arrowstyle="->", color="darkorange", lw=LW_ARROW_2D),
        )
        ax.annotate(
            "",
            xy=(x0, -R - 0.70 * scale),
            xytext=(x0, -R - 0.05 * scale),
            arrowprops=dict(arrowstyle="->", color="darkorange", lw=LW_ARROW_2D),
        )

    ax.text(
        0,
        R + 0.78 * scale,
        "открытое окно",
        color="darkorange",
        fontsize=8,
        ha="center",
        va="bottom",
        fontweight="bold",
    )


def draw_2d(ax, a, b, R, R1, R2_abs, closure, g):
    R_in = g["R_in"]
    dR_total = g["dR_total"]
    dR_half = g["dR_half"]
    L = g["L"]
    show_inlet = g["show_inlet"]
    show_outlet = g["show_outlet"]
    x_edge_l = g["x_edge_l"]
    rho_horn_cut = g["rho_horn_cut"]
    x_horn_out = g["x_horn_out"]
    x_torus_out = g["x_torus_out"]
    rho_cap_cut = g["rho_cap_cut"]

    # Левая гиперболическая воронка
    x_left = np.linspace(-L, -a, 900)
    rho_left = r_horn(x_left, a, b, R)
    mask_left = rho_left >= R_in - 1e-12 if show_inlet else np.ones_like(x_left, dtype=bool)

    ax.plot(x_left[mask_left], rho_left[mask_left], color="black", lw=LW_WALL_2D, label="гиперболическая стенка")
    ax.plot(x_left[mask_left], -rho_left[mask_left], color="black", lw=LW_WALL_2D)

    # Входное окно
    if show_inlet:
        ax.plot([x_edge_l, x_edge_l], [-R_in, R_in], color="royalblue", lw=LW_EDGE_2D)
        ax.annotate("", xy=(x_edge_l - 0.02 * max(1, a), 0.55 * R_in), xytext=(x_edge_l - 0.22 * L, 0.55 * R_in),
                    arrowprops=dict(arrowstyle="->", color="royalblue", lw=LW_ARROW_2D))
        ax.annotate("", xy=(x_edge_l - 0.02 * max(1, a), -0.55 * R_in), xytext=(x_edge_l - 0.22 * L, -0.55 * R_in),
                    arrowprops=dict(arrowstyle="->", color="royalblue", lw=LW_ARROW_2D))
        ax.text(x_edge_l, 0, "In-1", color="royalblue", fontsize=8, ha="right", va="center", fontweight="bold")

    # Правая гиперболическая воронка: вырез около rho=R под R2-окно
    x_right = np.linspace(a, L, 900)
    rho_right = r_horn(x_right, a, b, R)
    mask_right = rho_right <= rho_horn_cut + 1e-12 if show_outlet else np.ones_like(x_right, dtype=bool)

    ax.plot(x_right[mask_right], rho_right[mask_right], color="black", lw=LW_WALL_2D)
    ax.plot(x_right[mask_right], -rho_right[mask_right], color="black", lw=LW_WALL_2D)

    # Линии фокусов и фокусы
    ax.plot([-L, L], [R, R], color="purple", ls="--", lw=LW_AUX_2D, alpha=0.35)
    ax.plot([-L, L], [-R, -R], color="purple", ls="--", lw=LW_AUX_2D, alpha=0.35)
    plot_hyperbola_foci(ax, a, b, R)

    # Центральная область / заглушка
    if closure == "open":
        draw_open_window_arrows_2d(ax, a, R, L)
        cap_label = "открытый PHB"
    elif closure == "halftorus":
        x_tor = np.linspace(-a, a, 900)
        rho_tor = half_torus_outer_rho(x_tor, a, R)
        mask_torus = x_tor <= x_torus_out + 1e-12 if show_outlet else np.ones_like(x_tor, dtype=bool)

        ax.plot(x_tor[mask_torus], rho_tor[mask_torus], color="darkgreen", lw=LW_CAP_2D, label="полутороидальная заглушка")
        ax.plot(x_tor[mask_torus], -rho_tor[mask_torus], color="darkgreen", lw=LW_CAP_2D)
        cap_label = "полутороидальная заглушка"
    elif closure == "cylinder":
        x_cyl = np.linspace(-a, a, 600)
        ax.plot(x_cyl, np.full_like(x_cyl, R), color="darkgreen", lw=LW_CAP_2D, label="цилиндрическая заглушка")
        ax.plot(x_cyl, np.full_like(x_cyl, -R), color="darkgreen", lw=LW_CAP_2D)
        cap_label = "цилиндрическая заглушка"
    else:
        cap_label = "неизвестная заглушка"

    # Размер 2a
    ax.axvline(-a, color="blue", lw=LW_AUX_2D, ls=":", alpha=0.35)
    ax.axvline(a, color="blue", lw=LW_AUX_2D, ls=":", alpha=0.35)
    ax.annotate("", xy=(a, R + 0.18 * max(1, a)), xytext=(-a, R + 0.18 * max(1, a)),
                arrowprops=dict(arrowstyle="<->", color="gray", lw=LW_ARROW_2D))
    ax.text(0, R + 0.02 * max(1, a), "2a", ha="center", color="gray", fontsize=9, fontweight="bold")

    # Выходное окно R2 справа
    if show_outlet:
        ax.annotate("", xy=(x_horn_out + 0.25 * max(1, a), R + dR_half),
                    xytext=(x_horn_out + 0.25 * max(1, a), R - dR_half),
                    arrowprops=dict(arrowstyle="<->", color="red", lw=LW_ARROW_2D))
        ax.text(x_horn_out + 0.35 * max(1, a), R, f"R2=ΔR={dR_total:.3g}", color="red", fontsize=8,
                va="center", fontweight="bold")

        # Кромки окна, разные для разных заглушек
        ax.plot([x_horn_out], [rho_horn_cut], marker="o", color="red", ms=3)
        ax.plot([x_horn_out], [-rho_horn_cut], marker="o", color="red", ms=3)
        if closure == "halftorus":
            ax.plot([x_torus_out], [rho_cap_cut], marker="o", color="red", ms=3)
            ax.plot([x_torus_out], [-rho_cap_cut], marker="o", color="red", ms=3)
        elif closure == "cylinder":
            ax.plot([a], [R], marker="o", color="red", ms=3)
            ax.plot([a], [-R], marker="o", color="red", ms=3)

        ax.annotate("", xy=(L + 0.30 * L, R), xytext=(x_horn_out + 0.04 * L, R),
                    arrowprops=dict(arrowstyle="->", color="red", lw=LW_ARROW_2D))
        ax.annotate("", xy=(L + 0.30 * L, -R), xytext=(x_horn_out + 0.04 * L, -R),
                    arrowprops=dict(arrowstyle="->", color="red", lw=LW_ARROW_2D))

    # Заливка внутренней области
    ax.fill_between(x_left[mask_left], rho_left[mask_left], -rho_left[mask_left], color="lightsteelblue", alpha=0.08)
    ax.fill_between(x_right[mask_right], rho_right[mask_right], -rho_right[mask_right], color="lightsteelblue", alpha=0.08)
    x_mid = np.linspace(-a, a, 300)
    ax.fill_between(x_mid, R, -R, color="lightsteelblue", alpha=0.035)

    ax.set_xlabel("x")
    ax.set_ylabel("rho")
    ax.set_title(
        f"2D: {cap_label}; a={a:.2f}, b={b:.2f}, R={R:.2f}, R1={R1:.2f}, R2=ΔR={R2_abs:.2f}",
        fontsize=10,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.22)

    y_max = max(R + a, R + dR_half + a * 0.2, 0.5)
    ax.set_xlim(-L * 1.25 - 0.2 * L, L * 1.35)
    ax.set_ylim(-y_max * 1.45, y_max * 1.45)

    try:
        ax.legend(loc="lower right", fontsize=7)
    except Exception:
        pass


# ============================================================
# 3D
# ============================================================
def plot_ring_3d(ax3d, x0, rho0, color="red", lw=1.2, alpha=1.0, n=240, ls="-"):
    phi = np.linspace(0, 2 * np.pi, n)
    ax3d.plot(
        np.full_like(phi, x0),
        rho0 * np.cos(phi),
        rho0 * np.sin(phi),
        color=color,
        lw=lw,
        alpha=alpha,
        ls=ls,
    )


def draw_3d(ax3d, a, b, R, R1, R2_abs, closure, g):
    R_in = g["R_in"]
    L = g["L"]
    c = g["c"]
    show_inlet = g["show_inlet"]
    show_outlet = g["show_outlet"]
    x_edge_l = g["x_edge_l"]
    rho_horn_cut = g["rho_horn_cut"]
    x_horn_out = g["x_horn_out"]
    x_torus_out = g["x_torus_out"]
    rho_cap_cut = g["rho_cap_cut"]

    phi = np.linspace(0, 1.65 * np.pi, 90)
    phi_full = np.linspace(0, 2 * np.pi, 240)

    # Левая воронка
    x_lh = np.linspace(-L, -a, 120)
    Xl, PHIl = np.meshgrid(x_lh, phi)
    RHOl = r_horn(Xl, a, b, R)
    RHOl_draw = np.where(RHOl >= R_in - 1e-12, RHOl, np.nan) if show_inlet else RHOl
    ax3d.plot_surface(
        Xl,
        RHOl_draw * np.cos(PHIl),
        RHOl_draw * np.sin(PHIl),
        color="lightsteelblue",
        alpha=0.58,
        edgecolor="gray",
        linewidth=0.035,
        rstride=2,
        cstride=2,
    )

    # Кромка входного окна
    if show_inlet:
        plot_ring_3d(ax3d, x_edge_l, R_in, color="royalblue", lw=1.0)
        for ph in np.linspace(0, 2 * np.pi, 10, endpoint=False):
            y0 = 0.75 * R_in * np.cos(ph)
            z0 = 0.75 * R_in * np.sin(ph)
            ax3d.quiver(x_edge_l - L * 0.20, y0, z0, L * 0.16, 0, 0,
                        color="royalblue", lw=0.8, arrow_length_ratio=0.22, alpha=0.85)

    # Правая воронка
    x_rh = np.linspace(a, L, 120)
    Xr, PHIr = np.meshgrid(x_rh, phi)
    RHOr = r_horn(Xr, a, b, R)
    RHOr_draw = np.where(RHOr <= rho_horn_cut + 1e-12, RHOr, np.nan) if show_outlet else RHOr
    ax3d.plot_surface(
        Xr,
        RHOr_draw * np.cos(PHIr),
        RHOr_draw * np.sin(PHIr),
        color="lightsteelblue",
        alpha=0.58,
        edgecolor="gray",
        linewidth=0.035,
        rstride=2,
        cstride=2,
    )

    # Центральная заглушка / открытое окно
    if closure == "open":
        # Показываем два фокальных кольца и открытое центральное окно пунктиром.
        plot_ring_3d(ax3d, 0.0, R, color="darkorange", lw=1.0, alpha=0.85, ls="--")
        for ph in np.linspace(0, 2 * np.pi, 12, endpoint=False):
            y0 = R * np.cos(ph)
            z0 = R * np.sin(ph)
            # радиальный выход наружу от центрального кольца
            ax3d.quiver(0.0, y0, z0, 0, 0.35 * R * np.cos(ph), 0.35 * R * np.sin(ph),
                        color="darkorange", lw=0.8, arrow_length_ratio=0.22, alpha=0.85)
        cap_title = "открытый PHB"

    elif closure == "halftorus":
        theta = np.linspace(0, np.pi, 95)
        Theta, PhiT = np.meshgrid(theta, phi)
        Xtor = a * np.cos(Theta)
        RHOtor = R + a * np.sin(Theta)
        RHOtor_draw = np.where(Xtor <= x_torus_out + 1e-12, RHOtor, np.nan) if show_outlet else RHOtor
        ax3d.plot_surface(
            Xtor,
            RHOtor_draw * np.cos(PhiT),
            RHOtor_draw * np.sin(PhiT),
            color="mediumseagreen",
            alpha=0.62,
            edgecolor="darkgreen",
            linewidth=0.035,
            rstride=2,
            cstride=2,
        )
        cap_title = "полутороидальная заглушка"

    elif closure == "cylinder":
        x_c = np.linspace(-a, a, 100)
        Xc, PHIc = np.meshgrid(x_c, phi)
        RHOc = np.full_like(Xc, R)
        ax3d.plot_surface(
            Xc,
            RHOc * np.cos(PHIc),
            RHOc * np.sin(PHIc),
            color="mediumseagreen",
            alpha=0.62,
            edgecolor="darkgreen",
            linewidth=0.035,
            rstride=2,
            cstride=2,
        )
        cap_title = "цилиндрическая заглушка"
    else:
        cap_title = "неизвестная заглушка"

    # Фокальные кольца
    plot_ring_3d(ax3d, -c, R, color="magenta", lw=1.0, alpha=0.75, ls="--")
    plot_ring_3d(ax3d, c, R, color="magenta", lw=1.0, alpha=0.75, ls="--")

    # 3D-выходное окно отмечено красным
    if show_outlet:
        plot_ring_3d(ax3d, x_horn_out, rho_horn_cut, color="red", lw=2.0)
        if closure == "halftorus":
            plot_ring_3d(ax3d, x_torus_out, rho_cap_cut, color="red", lw=2.0)
        elif closure == "cylinder":
            plot_ring_3d(ax3d, a, R, color="red", lw=2.0)
        elif closure == "open":
            plot_ring_3d(ax3d, a, R, color="red", lw=1.5, alpha=0.9, ls="--")

        for ph in np.linspace(0, 2 * np.pi, 14, endpoint=False):
            y0 = R * np.cos(ph)
            z0 = R * np.sin(ph)
            ax3d.quiver(x_horn_out, y0, z0, L * 0.45, 0, 0,
                        color="red", lw=1.0, arrow_length_ratio=0.18, alpha=0.90)

    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.set_title(f"3D: {cap_title}; R2=ΔR={R2_abs:.2f}", fontsize=10)
    ax3d.view_init(elev=14, azim=36)

    mx = max(R + a, R + g["dR_half"], 0.5) * 1.10
    ax3d.set_xlim(-L * 1.25, L * 1.30)
    ax3d.set_ylim(-mx, mx)
    ax3d.set_zlim(-mx, mx)

    try:
        ax3d.set_box_aspect((2 * L * 1.275, 2 * mx, 2 * mx))
    except Exception:
        pass


# ============================================================
# WIDGETS
# ============================================================
# Увеличиваем ширину колонок и уменьшаем подписи/readout,
# чтобы числа возле ползунков не накладывались друг на друга.
SLIDER_STYLE = {"description_width": "42px"}
SLIDER_LAYOUT = widgets.Layout(width="260px")
DROPDOWN_STYLE = {"description_width": "70px"}
DROPDOWN_LAYOUT = widgets.Layout(width="390px")

closure_dropdown = widgets.Dropdown(
    options=CLOSURE_OPTIONS,
    value="open",
    description="Тип PHB",
    continuous_update=False,
    layout=DROPDOWN_LAYOUT,
    style=DROPDOWN_STYLE,
)

a_slider = widgets.FloatSlider(
    value=0.5,
    min=0.1,
    max=5.0,
    step=0.1,
    description="a",
    continuous_update=False,
    readout_format=".2f",
    layout=SLIDER_LAYOUT,
    style=SLIDER_STYLE,
)

b_slider = widgets.FloatSlider(
    value=1.0,
    min=0.1,
    max=5.0,
    step=0.1,
    description="b",
    continuous_update=False,
    readout_format=".2f",
    layout=SLIDER_LAYOUT,
    style=SLIDER_STYLE,
)

R_slider = widgets.FloatSlider(
    value=1.0,
    min=0.1,
    max=30.0,
    step=0.1,
    description="R",
    continuous_update=False,
    readout_format=".2f",
    layout=SLIDER_LAYOUT,
    style=SLIDER_STYLE,
)

R1_slider = widgets.FloatSlider(
    value=0.0,
    min=0.0,
    max=1.0,
    step=0.01,
    description="R1",
    continuous_update=False,
    readout_format=".2f",
    layout=SLIDER_LAYOUT,
    style=SLIDER_STYLE,
)

R2_slider = widgets.FloatSlider(
    value=0.0,
    min=0.0,
    max=2.0,
    step=0.01,
    description="R2",
    continuous_update=False,
    readout_format=".2f",
    layout=SLIDER_LAYOUT,
    style=SLIDER_STYLE,
)

out = widgets.Output()


def redraw(change=None):
    closure = str(closure_dropdown.value)
    a = float(a_slider.value)
    b = float(b_slider.value)
    R = float(R_slider.value)
    R1 = float(R1_slider.value)
    R2_abs = float(R2_slider.value)

    g = compute_geometry(a, b, R, R1, R2_abs, closure)

    with out:
        clear_output(wait=True)

        if g["clipped"]:
            display(HTML(
                f"<b style='color:#b00000'>Предупреждение:</b> "
                f"для выбранных параметров половина выходного окна R2/2 = "
                f"{g['dR_half_raw']:.3f} выходит за отображаемую область. "
                f"Для построения временно использовано R2/2 = {g['dR_half']:.3f}. "
                f"Уменьшите R2 или увеличьте a/R."
            ))

        fig = plt.figure(figsize=(15.2, 6.4))
        ax2d = fig.add_subplot(1, 2, 1)
        ax3d = fig.add_subplot(1, 2, 2, projection="3d")

        draw_2d(ax2d, a, b, R, R1, R2_abs, closure, g)
        draw_3d(ax3d, a, b, R, R1, R2_abs, closure, g)

        plt.tight_layout()
        display(fig)
        plt.close(fig)


for w in [closure_dropdown, a_slider, b_slider, R_slider, R1_slider, R2_slider]:
    w.observe(redraw, names="value")

controls = widgets.GridBox(
    children=[closure_dropdown, a_slider, b_slider, R_slider, R1_slider, R2_slider],
    layout=widgets.Layout(
        grid_template_columns="390px repeat(5, 260px)",
        grid_gap="8px 10px",
        align_items="center",
    ),
)

display(HTML(
    """
    <style>
    /* Jupyter/Colab widget label/readout compacting */
    .widget-label { font-size: 11px !important; }
    .widget-readout { font-size: 11px !important; min-width: 44px !important; }
    .jupyter-widgets .widget-label { font-size: 11px !important; }
    .jupyter-widgets .widget-readout { font-size: 11px !important; min-width: 44px !important; }
    </style>
    """
))

display(HTML(
    "<h3>Псевдогиперболоид: интерактивная форма с выбором заглушки</h3>"
    "<div>Сверху отображаются 2D и 3D. Снизу меняйте параметры "
    "<b>тип PHB, a, b, R, R1, R2</b>. "
    "<br><b>R1</b> — доля радиуса входного окна: R_in = R·R1. "
    "<br><b>R2</b> — абсолютная радиальная ширина выходного окна ΔR, "
    "а не относительная величина R2/R.</div>"
))

display(widgets.VBox([out, controls]))

redraw()

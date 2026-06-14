#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Скрипт для расчёта Moore's Modified Rayleigh test
и построения круговой диаграммы по данным из Excel.

Формат входного Excel-файла:
    1-й столбец: deg  — направление в градусах
    2-й столбец: r    — длина вектора / выраженность ориентации

Названия столбцов могут быть любыми: скрипт берёт именно первые два столбца.

Пример запуска:
    python moore_rayleigh_excel_plot_final_dots_CI_lines_poster_style_v11.py data.xlsx --title "Group 1"

Если нужно оставить только строки, где r > 0.2:
    python moore_rayleigh_excel_plot_final_dots_CI_lines_poster_style_v11.py data.xlsx --title "Group 1" --threshold 0.2

Быстрый тестовый запуск:
    python moore_rayleigh_excel_plot_final_dots_CI_lines_poster_style_v11.py data.xlsx --title "Test" --n-bootstrap 1000 --n-monte-carlo 5000

Главные изменения:
    - график рисуется не через polar-ось matplotlib, а в обычных декартовых координатах
      с равным масштабом по X и Y, поэтому стрелки не сплющиваются;
    - добавлены три границы значимости для вектора по Рейли: p = 0.05, 0.01 и 0.001;
    - индивидуальные стрелки сделаны серыми;
    - средняя стрелка сделана чёрной;
    - верхняя подпись 0 заменена на mN/gN.
"""

import argparse
import math
import random
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, FancyArrowPatch, FancyBboxPatch, Rectangle
from tkinter import Tk, filedialog


# ============================================================
# 1. Вспомогательные функции
# ============================================================

def to_float(value):
    """
    Преобразует значение из Excel в число.

    Функция понимает оба варианта записи десятичных чисел:
        0.87
        0,87
    """
    if pd.isna(value):
        return np.nan

    if isinstance(value, str):
        value = value.strip().replace(",", ".")

    return float(value)


def average_ranks(values):
    """
    Возвращает ранги для списка значений.

    Самое маленькое значение получает ранг 1.
    Самое большое значение получает ранг n.

    Если есть одинаковые значения, им присваивается средний ранг.
    """
    idx = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)

    current_rank = 1
    i = 0

    while i < len(idx):
        j = i

        while j < len(idx) and abs(values[idx[j]] - values[idx[i]]) < 1e-12:
            j += 1

        avg_rank = (current_rank + current_rank + (j - i) - 1) / 2

        for k in range(i, j):
            ranks[idx[k]] = avg_rank

        current_rank += (j - i)
        i = j

    return ranks


def angular_diff(a, b):
    """
    Вычисляет круговую разницу между двумя углами.

    Возвращает значение a - b в диапазоне [-180, 180).
    """
    return ((a - b + 180) % 360) - 180


def xy_from_bearing(deg, radius=1.0):
    """
    Переводит направление и длину вектора в координаты x, y.

    Используется биологическая / навигационная система:
        0°   = север
        90°  = восток
        180° = юг
        270° = запад

    Угол увеличивается по часовой стрелке.
    """
    theta = math.radians(deg)

    x = radius * math.sin(theta)
    y = radius * math.cos(theta)

    return x, y


def make_sector_points(ci_low, ci_high, radius=1.0, n_points=500):
    """
    Создаёт точки для заливки доверительного сектора.

    Корректно работает, если интервал пересекает 0°.
    Например: 350–20°.
    """
    if ci_low <= ci_high:
        degs = np.linspace(ci_low, ci_high, n_points)
    else:
        # Интервал пересекает 0°.
        degs = np.linspace(ci_low, ci_high + 360, n_points)

    points = [(0.0, 0.0)]

    for deg in degs:
        points.append(xy_from_bearing(deg, radius))

    points.append((0.0, 0.0))

    return points


def add_constant_arrow(
    ax,
    angle_deg,
    length,
    color="purple",
    linewidth=2.2,
    mutation_scale=16,
    alpha=1.0,
    zorder=5,
):
    """
    Рисует стрелку постоянной визуальной толщины.

    Важно:
        Это не ax.arrow() на polar-оси.
        FancyArrowPatch рисуется в обычных X/Y координатах при ax.set_aspect("equal"),
        поэтому стрелки не сплющиваются и не меняют ширину в зависимости от направления.

    В этой версии стрелки сделаны без сглаженных / округлённых углов:
        joinstyle="miter"
        capstyle="butt"
    """
    x_end, y_end = xy_from_bearing(angle_deg, length)

    arrow = FancyArrowPatch(
        (0.0, 0.0),
        (x_end, y_end),
        arrowstyle="-|>",
        mutation_scale=mutation_scale,
        linewidth=linewidth,
        color=color,
        alpha=alpha,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
        joinstyle="miter",
        capstyle="butt",
        antialiased=True,
    )

    ax.add_patch(arrow)

    return arrow


# ============================================================
# 2. Moore's Modified Rayleigh test
# ============================================================

def moore_modified_rayleigh(deg, r):
    """
    Считает Moore's Modified Rayleigh test.

    Вход:
        deg — массив направлений в градусах
        r   — массив длин индивидуальных векторов

    Алгоритм:
        1. Значения r заменяются рангами.
        2. Направления взвешиваются этими рангами.
        3. Считается результирующий вектор.
    """
    deg = np.asarray(deg, dtype=float)
    r = np.asarray(r, dtype=float)

    if len(deg) != len(r):
        raise ValueError("Массивы deg и r должны иметь одинаковую длину.")

    if len(deg) < 2:
        raise ValueError("Для расчёта нужно минимум 2 вектора.")

    ranks = np.asarray(average_ranks(r), dtype=float)

    angles_rad = np.deg2rad(deg)

    C = np.sum(ranks * np.cos(angles_rad))
    S = np.sum(ranks * np.sin(angles_rad))

    L = math.hypot(C, S)

    n = len(deg)

    alpha_deg = math.degrees(math.atan2(S, C)) % 360

    r_rank_weighted = L / np.sum(ranks)

    R_star = L / (n ** 1.5)

    return {
        "n": n,
        "alpha_deg": alpha_deg,
        "r_rank_weighted": r_rank_weighted,
        "R_star": R_star,
        "C": C,
        "S": S,
        "L": L,
        "ranks": ranks,
    }


def monte_carlo_p_value(r, observed_L, n_iter=200000, seed=123):
    """
    Оценивает p-value методом Monte Carlo.

    Нулевая гипотеза:
        направления распределены равномерно по кругу.
    """
    rng = np.random.default_rng(seed)

    ranks = np.asarray(average_ranks(r), dtype=float)
    n = len(ranks)

    count = 0

    for _ in range(n_iter):
        random_angles = rng.uniform(0, 2 * np.pi, n)

        C = np.sum(ranks * np.cos(random_angles))
        S = np.sum(ranks * np.sin(random_angles))
        L = math.hypot(C, S)

        if L >= observed_L:
            count += 1

    p_value = (count + 1) / (n_iter + 1)

    return p_value


def bootstrap_confidence_interval(deg, r, n_boot=100000, seed=123):
    """
    Считает 95% доверительный интервал для среднего направления методом bootstrap.

    Используется центрированный доверительный сектор:
        alpha ± q95(|alpha_bootstrap - alpha|)
    """
    rng = random.Random(seed)

    deg = list(map(float, deg))
    r = list(map(float, r))
    n = len(deg)

    observed = moore_modified_rayleigh(deg, r)
    alpha = observed["alpha_deg"]

    abs_diffs = []

    for _ in range(n_boot):
        indices = [rng.randrange(n) for _ in range(n)]

        deg_sample = [deg[i] for i in indices]
        r_sample = [r[i] for i in indices]

        boot = moore_modified_rayleigh(deg_sample, r_sample)

        diff = abs(angular_diff(boot["alpha_deg"], alpha))
        abs_diffs.append(diff)

    half_width = float(np.quantile(abs_diffs, 0.95))

    ci_low = (alpha - half_width) % 360
    ci_high = (alpha + half_width) % 360

    return ci_low, ci_high, half_width


# ============================================================
# 3. Загрузка Excel
# ============================================================

def load_excel_two_columns(path, threshold=None):
    """
    Загружает Excel-файл.

    Требование к файлу:
        1-й столбец = deg
        2-й столбец = r

    Названия столбцов не важны.
    """
    df = pd.read_excel(path)

    if df.shape[1] < 2:
        raise ValueError("В Excel-файле должно быть минимум два столбца: deg и r.")

    df = df.iloc[:, :2].copy()
    df.columns = ["deg", "r"]

    df["deg"] = df["deg"].apply(to_float)
    df["r"] = df["r"].apply(to_float)

    df = df.dropna(subset=["deg", "r"]).copy()

    df["deg"] = df["deg"] % 360

    if threshold is not None:
        df = df[df["r"] > threshold].copy()

    if len(df) < 2:
        raise ValueError("После фильтрации осталось меньше 2 строк. Расчёт невозможен.")

    return df


# ============================================================
# 4. Построение кругового графика
# ============================================================

def draw_polar_plot(df, stats, ci_low, ci_high, p_value, title, output_prefix):
    """
    Строит круговую диаграмму в аккуратном постерном / публикационном стиле.

    В этой версии нижний информационный блок полностью собран заново:
        - убрана вертикальная серая линия;
        - Legend и Results расположены симметрично;
        - отступы и интервалы выровнены;
        - заголовки подняты ближе к верхней границе блока;
        - нижний блок собран как единая спокойная композиция.
    """
    output_prefix = Path(output_prefix)

    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    svg_path = output_prefix.with_suffix(".svg")

    alpha = stats["alpha_deg"]
    r_group = stats["r_rank_weighted"]
    R_star = stats["R_star"]
    n = stats["n"]

    # Rayleigh-style r threshold circles
    r_05 = math.sqrt(-math.log(0.05) / n)
    r_01 = math.sqrt(-math.log(0.01) / n)
    r_001 = math.sqrt(-math.log(0.001) / n)

    fig = plt.figure(figsize=(8.8, 11.1), dpi=300, facecolor="white")

    # ------------------------------------------------------------
    # Верхний график
    # ------------------------------------------------------------
    ax = fig.add_axes([0.08, 0.405, 0.84, 0.535])

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-1.22, 1.22)
    ax.set_ylim(-1.22, 1.22)
    ax.axis("off")

    # 1. Серый сектор 95% CI
    sector_points = make_sector_points(ci_low, ci_high, radius=1.0, n_points=500)
    ci_sector = Polygon(
        sector_points,
        closed=True,
        facecolor="#c4c4c4",
        edgecolor="none",
        zorder=0,
    )
    ax.add_patch(ci_sector)

    # 2. Светло-серые осевые линии
    ax.plot([0, 0], [-1.0, 1.0], color="#d3d3d3", linewidth=1.15, zorder=1)
    ax.plot([-1.0, 1.0], [0, 0], color="#d3d3d3", linewidth=1.15, zorder=1)

    # 3. Границы значимости Rayleigh
    significance_circles = [
        (r_05, "--", 1.6, "0.05"),
        (r_01, ":", 2.0, "0.01"),
        (r_001, "-.", 2.0, "0.001"),
    ]

    for radius, linestyle, linewidth, label in significance_circles:
        circle = Circle(
            (0, 0),
            radius,
            facecolor="none",
            edgecolor="black",
            linestyle=linestyle,
            linewidth=linewidth,
            zorder=2,
        )
        ax.add_patch(circle)

        if radius < 0.98:
            label_angle = 55
            x_label, y_label = xy_from_bearing(label_angle, radius)
            ax.text(
                x_label + 0.024,
                y_label + 0.024,
                label,
                fontsize=8.5,
                ha="left",
                va="center",
                color="black",
                zorder=3,
            )

    # 4. Внешняя граница круга
    outer_circle = Circle(
        (0, 0),
        1.0,
        facecolor="none",
        edgecolor="black",
        linewidth=1.8,
        zorder=4,
    )
    ax.add_patch(outer_circle)

    # 5. Индивидуальные стрелки
    for _, row in df.iterrows():
        add_constant_arrow(
            ax=ax,
            angle_deg=row["deg"],
            length=row["r"],
            color="#6e6e6e",
            linewidth=2.0,
            mutation_scale=15,
            alpha=0.72,
            zorder=5,
        )

    # 6. Голубые точки по краю круга
    for _, row in df.iterrows():
        x_dot, y_dot = xy_from_bearing(row["deg"], 1.03)
        ax.scatter(
            x_dot,
            y_dot,
            s=110,
            facecolor="#18a8e1",
            edgecolor="black",
            linewidth=1.05,
            zorder=10,
            clip_on=False,
        )

    # 7. Средний групповой вектор
    add_constant_arrow(
        ax=ax,
        angle_deg=alpha,
        length=r_group,
        color="black",
        linewidth=3.0,
        mutation_scale=23,
        alpha=1.0,
        zorder=8,
    )

    # 8. Подписи направлений
    ax.text(0, 1.16, "mN/gN", ha="center", va="bottom", fontsize=14, fontweight="bold")
    ax.text(1.16, 0, "90", ha="left", va="center", fontsize=14, fontweight="bold")
    ax.text(0, -1.12, "180", ha="center", va="top", fontsize=14, fontweight="bold")
    ax.text(-1.16, 0, "270", ha="right", va="center", fontsize=14, fontweight="bold")

    # 9. Заголовок
    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        pad=24,
    )

    # ------------------------------------------------------------
    # Нижний единый информационный блок
    # ------------------------------------------------------------
    panel_ax = fig.add_axes([0.06, 0.045, 0.88, 0.315])
    panel_ax.set_xlim(0, 1)
    panel_ax.set_ylim(0, 1)
    panel_ax.axis("off")

    # Внешняя рамка блока.
    # Используем прямоугольник без скруглений, чтобы верхняя и нижняя
    # чёрные линии были ровными и с прямыми углами.
    main_box = Rectangle(
        (0.00, 0.03),
        1.00,
        0.94,
        linewidth=1.2,
        edgecolor="#2f2f2f",
        facecolor="#f7f7f7",
        transform=panel_ax.transAxes,
        zorder=0,
    )
    panel_ax.add_patch(main_box)

    # ------------------------------------------------------------
    # Две симметричные секции без вертикального разделителя
    # ------------------------------------------------------------
    legend_left, legend_right = 0.06, 0.46
    results_left, results_right = 0.54, 0.94

    legend_center = (legend_left + legend_right) / 2
    results_center = (results_left + results_right) / 2

    # Заголовки ближе к верхней границе блока
    panel_ax.text(
        legend_center, 0.88, "Legend",
        fontsize=12.4, fontweight="bold",
        ha="center", va="center", transform=panel_ax.transAxes,
    )
    panel_ax.text(
        results_center, 0.88, "Results",
        fontsize=12.4, fontweight="bold",
        ha="center", va="center", transform=panel_ax.transAxes,
    )

    panel_ax.text(
        legend_center, 0.81, "Graph elements",
        fontsize=8.8, color="#555555",
        ha="center", va="center", transform=panel_ax.transAxes,
    )
    panel_ax.text(
        results_center, 0.81, "Statistical summary",
        fontsize=8.8, color="#555555",
        ha="center", va="center", transform=panel_ax.transAxes,
    )

    # ------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------
    # 7 ровных строк
    legend_rows = [0.68, 0.58, 0.48, 0.38, 0.27, 0.18, 0.09]

    # Геометрия символов и текста в левой секции
    symbol_center_x = legend_left + 0.095
    sample_half_w = 0.050
    text_x = legend_left + 0.185

    # individual vector
    panel_ax.annotate(
        "",
        xy=(symbol_center_x + sample_half_w, legend_rows[0]),
        xytext=(symbol_center_x - sample_half_w, legend_rows[0]),
        xycoords=panel_ax.transAxes,
        textcoords=panel_ax.transAxes,
        arrowprops=dict(
            arrowstyle="-|>",
            color="#6e6e6e",
            lw=2.0,
            mutation_scale=14,
            shrinkA=0,
            shrinkB=0,
            joinstyle="miter",
            capstyle="butt",
        ),
    )
    panel_ax.text(
        text_x, legend_rows[0], "individual vector",
        fontsize=9.5, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # mean vector
    panel_ax.annotate(
        "",
        xy=(symbol_center_x + sample_half_w, legend_rows[1]),
        xytext=(symbol_center_x - sample_half_w, legend_rows[1]),
        xycoords=panel_ax.transAxes,
        textcoords=panel_ax.transAxes,
        arrowprops=dict(
            arrowstyle="-|>",
            color="black",
            lw=2.6,
            mutation_scale=16,
            shrinkA=0,
            shrinkB=0,
            joinstyle="miter",
            capstyle="butt",
        ),
    )
    panel_ax.text(
        text_x, legend_rows[1], "mean vector",
        fontsize=9.5, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # individual direction
    panel_ax.scatter(
        [symbol_center_x], [legend_rows[2]],
        s=64,
        facecolor="#18a8e1",
        edgecolor="black",
        linewidth=0.9,
        transform=panel_ax.transAxes,
        zorder=3,
    )
    panel_ax.text(
        text_x, legend_rows[2], "individual direction",
        fontsize=9.5, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # 95% CI sector
    sector_w = 0.065
    sector_h = 0.038
    panel_ax.add_patch(
        FancyBboxPatch(
            (symbol_center_x - sector_w / 2, legend_rows[3] - sector_h / 2),
            sector_w,
            sector_h,
            boxstyle="round,pad=0.002,rounding_size=0.005",
            linewidth=0.8,
            edgecolor="#7c7c7c",
            facecolor="#c4c4c4",
            transform=panel_ax.transAxes,
            zorder=2,
        )
    )
    panel_ax.text(
        text_x, legend_rows[3], "95% CI sector",
        fontsize=9.5, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # Rayleigh 0.05
    panel_ax.plot(
        [symbol_center_x - sample_half_w, symbol_center_x + sample_half_w],
        [legend_rows[4], legend_rows[4]],
        color="black",
        linestyle="--",
        linewidth=1.8,
        transform=panel_ax.transAxes,
        zorder=2,
    )
    panel_ax.text(
        text_x, legend_rows[4], "Rayleigh 0.05",
        fontsize=9.2, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # Rayleigh 0.01
    panel_ax.plot(
        [symbol_center_x - sample_half_w, symbol_center_x + sample_half_w],
        [legend_rows[5], legend_rows[5]],
        color="black",
        linestyle=":",
        linewidth=1.8,
        transform=panel_ax.transAxes,
        zorder=2,
    )
    panel_ax.text(
        text_x, legend_rows[5], "Rayleigh 0.01",
        fontsize=9.2, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # Rayleigh 0.001
    panel_ax.plot(
        [symbol_center_x - sample_half_w, symbol_center_x + sample_half_w],
        [legend_rows[6], legend_rows[6]],
        color="black",
        linestyle="-.",
        linewidth=1.8,
        transform=panel_ax.transAxes,
        zorder=2,
    )
    panel_ax.text(
        text_x, legend_rows[6], "Rayleigh 0.001",
        fontsize=9.2, ha="left", va="center", transform=panel_ax.transAxes,
    )

    # ------------------------------------------------------------
    # Results — компактная аккуратная таблица
    # ------------------------------------------------------------
    table_left = results_left + 0.015
    table_right = results_right - 0.015
    table_top = 0.74
    table_bottom = 0.12

    panel_ax.add_patch(
        FancyBboxPatch(
            (table_left, table_bottom),
            table_right - table_left,
            table_top - table_bottom,
            boxstyle="round,pad=0.008,rounding_size=0.015",
            linewidth=0.8,
            edgecolor="#d0d0d0",
            facecolor="white",
            transform=panel_ax.transAxes,
            zorder=1,
        )
    )

    results_rows = [
        ("n", f"{n}"),
        ("mean direction α", f"{alpha:.1f}°"),
        ("95% CI", f"{ci_low:.1f}°–{ci_high:.1f}°"),
        ("rank-weighted r", f"{r_group:.3f}"),
        ("R*", f"{R_star:.3f}"),
        ("Monte Carlo p", f"{p_value:.3g}"),
        ("Rayleigh r0.05", f"{r_05:.3f}"),
        ("Rayleigh r0.01", f"{r_01:.3f}"),
        ("Rayleigh r0.001", f"{r_001:.3f}"),
    ]

    n_rows = len(results_rows) + 1
    row_h = (table_top - table_bottom) / n_rows
    split_x = table_left + 0.63 * (table_right - table_left)

    # Шапка таблицы
    panel_ax.add_patch(
        FancyBboxPatch(
            (table_left, table_top - row_h),
            table_right - table_left,
            row_h,
            boxstyle="round,pad=0.004,rounding_size=0.010",
            linewidth=0,
            edgecolor="none",
            facecolor="#ececec",
            transform=panel_ax.transAxes,
            zorder=2,
        )
    )

    # Внутренние линии таблицы
    panel_ax.plot(
        [split_x, split_x],
        [table_bottom, table_top],
        color="#dddddd",
        linewidth=0.8,
        transform=panel_ax.transAxes,
        zorder=3,
    )

    for i in range(1, n_rows):
        y0 = table_top - i * row_h
        panel_ax.plot(
            [table_left, table_right],
            [y0, y0],
            color="#e6e6e6",
            linewidth=0.8,
            transform=panel_ax.transAxes,
            zorder=3,
        )

    # Лёгкая зебра
    for idx in range(len(results_rows)):
        if idx % 2 == 0:
            y_top = table_top - (idx + 1) * row_h
            y_bottom = table_top - (idx + 2) * row_h
            panel_ax.fill_between(
                [table_left, table_right],
                [y_bottom, y_bottom],
                [y_top, y_top],
                color="#fafafa",
                transform=panel_ax.transAxes,
                zorder=1.5,
            )

    # Заголовки таблицы
    panel_ax.text(
        (table_left + split_x) / 2,
        table_top - row_h / 2,
        "Metric",
        fontsize=9.0,
        fontweight="bold",
        ha="center",
        va="center",
        transform=panel_ax.transAxes,
        zorder=4,
    )
    panel_ax.text(
        (split_x + table_right) / 2,
        table_top - row_h / 2,
        "Value",
        fontsize=9.0,
        fontweight="bold",
        ha="center",
        va="center",
        transform=panel_ax.transAxes,
        zorder=4,
    )

    for idx, (label, value) in enumerate(results_rows, start=1):
        y_center = table_top - (idx + 0.5) * row_h

        panel_ax.text(
            table_left + 0.012,
            y_center,
            label,
            fontsize=8.5,
            ha="left",
            va="center",
            transform=panel_ax.transAxes,
            zorder=4,
        )
        panel_ax.text(
            (split_x + table_right) / 2,
            y_center,
            value,
            fontsize=8.7,
            fontweight="bold",
            ha="center",
            va="center",
            transform=panel_ax.transAxes,
            zorder=4,
        )

    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(pdf_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(svg_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    return png_path, pdf_path, svg_path


# ============================================================
# 5. Выбор Excel-файла
# ============================================================

def choose_excel_file():
    """
    Открывает окно выбора Excel-файла, если путь не был передан через командную строку.
    """
    root = Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Select Excel file with circular data",
        filetypes=[
            ("Excel files", "*.xlsx *.xls"),
            ("All files", "*.*"),
        ],
    )

    if not file_path:
        raise SystemExit("Файл не выбран.")

    return file_path


# ============================================================
# 6. Основной запуск скрипта
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Расчёт Moore's Modified Rayleigh test и построение круговой диаграммы по Excel-файлу."
    )

    parser.add_argument(
        "excel_file",
        nargs="?",
        default=None,
        help="Путь к Excel-файлу. Если не указан, откроется окно выбора файла.",
    )

    parser.add_argument(
        "--title",
        default="Moore's Modified Rayleigh test",
        help="Заголовок графика.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Фильтр: оставить только строки с r > threshold. Например: --threshold 0.2",
    )

    parser.add_argument(
        "--output",
        default=None,
        help="Префикс выходных файлов. Будут сохранены PNG, PDF и SVG. Если не указан, файлы сохраняются рядом с Excel-файлом.",
    )

    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=100000,
        help="Число bootstrap-повторов для 95%% CI.",
    )

    parser.add_argument(
        "--n-monte-carlo",
        type=int,
        default=200000,
        help="Число Monte Carlo-повторов для p-value.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Случайное зерно для воспроизводимости расчётов.",
    )

    args = parser.parse_args()

    if args.excel_file is None:
        args.excel_file = choose_excel_file()

    if args.output is None:
        input_path = Path(args.excel_file)
        args.output = str(input_path.with_name(input_path.stem + "_moore_rayleigh_result"))

    df = load_excel_two_columns(args.excel_file, threshold=args.threshold)

    deg = df["deg"].to_numpy()
    r = df["r"].to_numpy()

    stats = moore_modified_rayleigh(deg, r)

    ci_low, ci_high, ci_half_width = bootstrap_confidence_interval(
        deg,
        r,
        n_boot=args.n_bootstrap,
        seed=args.seed,
    )

    p_value = monte_carlo_p_value(
        r,
        observed_L=stats["L"],
        n_iter=args.n_monte_carlo,
        seed=args.seed,
    )

    png_path, pdf_path, svg_path = draw_polar_plot(
        df=df,
        stats=stats,
        ci_low=ci_low,
        ci_high=ci_high,
        p_value=p_value,
        title=args.title,
        output_prefix=args.output,
    )

    print("\nMoore's Modified Rayleigh test")
    print("--------------------------------")
    print(f"Файл: {args.excel_file}")
    print(f"Использовано строк: {stats['n']}")

    if args.threshold is not None:
        print(f"Фильтр: r > {args.threshold}")

    print(f"Среднее направление, alpha: {stats['alpha_deg']:.3f}°")
    print(f"Рангово-взвешенное r:       {stats['r_rank_weighted']:.6f}")
    print(f"R*_n:                       {stats['R_star']:.6f}")
    print(f"95% CI:                     {ci_low:.3f}–{ci_high:.3f}°")
    print(f"Полуширина 95% CI:          ±{ci_half_width:.3f}°")
    print(f"Monte Carlo p-value:        {p_value:.6g}")

    print("\nСохранённые файлы:")
    print(f"  {png_path}")
    print(f"  {pdf_path}")
    print(f"  {svg_path}")


if __name__ == "__main__":
    main()

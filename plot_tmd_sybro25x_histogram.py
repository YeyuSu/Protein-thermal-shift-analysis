#!/usr/bin/env python3
"""
从 tm_summary.csv 筛选 sybro2.5x、排除 12*/12_、WT 排除 4*/4_ 后，
绘制 TmD 均值 ± SD 柱状图；首柱为 IL10-standard（标签 WT_standard），
其余按真实 MUT 编号排序（mut1, mut13, …）。
阈值：WT_standard 的 TmD 均值（黑色虚线）。
色带与 analyze_thermal_shift.PALETTE_HEX 一致。
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

PALETTE_HEX = ["#92B1D9", "#C1D8E9", "#DBDDEF", "#F6C8B6", "#D4D4D4"]


def _hex_to_rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def make_gradient_cmap() -> LinearSegmentedColormap:
    colors = [_hex_to_rgb(x) for x in PALETTE_HEX]
    return LinearSegmentedColormap.from_list("thermal_shift_palette", colors, N=256)


def is_sybro25x(sample_name: str) -> bool:
    return bool(re.search(r"sybro2\.5x", sample_name, flags=re.I))


def has_12_condition(sample_name: str) -> bool:
    return "12*" in sample_name or "12_" in sample_name


def is_wt_row(sample_name: str) -> bool:
    return bool(re.search(r"IL10-WT(?:-|$)", sample_name, re.I))


def wt_has_4_condition(sample_name: str) -> bool:
    return "4*" in sample_name or "4_" in sample_name


def is_standard_row(sample_name: str) -> bool:
    return bool(re.search(r"IL10-standard", sample_name, re.I))


def parse_mut_number(sample_name: str) -> int | None:
    m = re.search(r"IL10-MUT(\d+)", sample_name, re.I)
    return int(m.group(1)) if m else None


def pick_tm_columns(df: pd.DataFrame) -> tuple[str, str]:
    """返回 (mean_col, sd_col)，优先 ROA1。"""
    for k in range(1, 7):
        mean_c = f"TmD_ROA{k}_mean"
        sd_c = f"TmD_ROA{k}_sd"
        if mean_c in df.columns and sd_c in df.columns:
            return mean_c, sd_c
    raise ValueError("未找到 TmD_ROA*_mean / *_sd 列")


def filter_summary(df: pd.DataFrame) -> pd.DataFrame:
    keep: list[pd.Series] = []
    for _, row in df.iterrows():
        sn = str(row["Sample_Name"])
        if not is_sybro25x(sn):
            continue
        if has_12_condition(sn):
            continue
        if is_wt_row(sn) and wt_has_4_condition(sn):
            continue
        keep.append(row)
    return pd.DataFrame(keep).reset_index(drop=True)


def build_plot_table(df: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """返回按作图顺序的表（含 x_label）与 standard 的 TmD 阈值。"""
    mean_c, sd_c = pick_tm_columns(df)
    std_rows = df[df["Sample_Name"].map(is_standard_row)]
    if std_rows.empty:
        raise ValueError("筛选结果中未找到 IL10-standard 行，无法设定 WT_standard 阈值。")
    if len(std_rows) > 1:
        raise ValueError(f"存在多行 standard，请检查数据：{std_rows['Sample_Name'].tolist()}")
    threshold = float(std_rows.iloc[0][mean_c])

    mut_df = df[~df["Sample_Name"].map(is_standard_row)].copy()
    mut_df["_mut"] = mut_df["Sample_Name"].map(parse_mut_number)
    if mut_df["_mut"].isna().any():
        bad = mut_df.loc[mut_df["_mut"].isna(), "Sample_Name"].tolist()
        raise ValueError(f"以下行无法解析 MUT 编号：{bad}")

    dup = mut_df["_mut"].duplicated(keep=False)
    if dup.any():
        dups = mut_df.loc[dup, ["Sample_Name", "_mut"]]
        raise ValueError(f"同一 MUT 编号出现多行，请合并或去重后再作图：\n{dups}")

    mut_df = mut_df.sort_values("_mut", kind="stable")
    std_row = std_rows.iloc[0]
    plot_rows: list[dict] = [
        {
            "x_label": "WT_standard",
            "Sample_Name": std_row["Sample_Name"],
            "mean": float(std_row[mean_c]),
            "sd": float(std_row[sd_c]),
            "mut_num": None,
        }
    ]
    for _, r in mut_df.iterrows():
        n = int(r["_mut"])
        plot_rows.append(
            {
                "x_label": f"mut{n}",
                "Sample_Name": r["Sample_Name"],
                "mean": float(r[mean_c]),
                "sd": float(r[sd_c]),
                "mut_num": n,
            }
        )
    return pd.DataFrame(plot_rows), threshold


def plot_bars(plot_df: pd.DataFrame, threshold: float, out_path: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "axes.unicode_minus": False,
        }
    )
    cmap = make_gradient_cmap()
    n = len(plot_df)
    x = np.arange(n)
    means = plot_df["mean"].to_numpy(dtype=float)
    sds = plot_df["sd"].to_numpy(dtype=float)
    colors = [cmap(i / max(n - 1, 1)) for i in range(n)]

    fig, ax = plt.subplots(figsize=(max(10.0, 0.45 * n), 5.0), dpi=150)
    ax.bar(x, means, yerr=sds, color=colors, edgecolor="white", linewidth=0.6, capsize=3, zorder=2)
    ax.axhline(threshold, color="k", linestyle="--", linewidth=1.2, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["x_label"], rotation=45, ha="right")
    ax.set_ylabel("Tm(°C)", fontfamily="Arial")
    ax.set_ylim(40, 65)
    ax.set_title("IL-10 mutants thermostability", fontsize=16, fontfamily="Arial")
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontfamily("Arial")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description="TmD 柱状图（sybro2.5x 筛选 + PALETTE 渐变）")
    p.add_argument(
        "--summary",
        type=Path,
        default=here / "output_thermal_shift_pymol" / "tm_summary.csv",
        help="tm_summary.csv 路径",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=here / "output_thermal_shift_pymol" / "figures" / "tmd_bars_sybro25x_filtered.png",
        help="输出 PNG",
    )
    args = p.parse_args()

    df = pd.read_csv(args.summary)
    filt = filter_summary(df)
    print(f"筛选后 {len(filt)} 行（来自 {args.summary}）:")
    for sn in filt["Sample_Name"]:
        print(f"  - {sn}")

    plot_df, thr = build_plot_table(filt)
    print(f"作图 {len(plot_df)} 柱，WT_standard 阈值 TmD = {thr:.4f} °C")
    plot_bars(plot_df, thr, args.out)
    print(f"已保存 {args.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
IL-10 QuantStudio xlsx 的 Protein thermal shift 后处理（对齐 PTS 说明书思路）。

阶段 1：按 ROA 拾 TmD——默认对 Fluorescence 求数值 dF/dT；可选 --tmd-source instrument_derivative
        改为对仪器导出的 Derivative 列在 ROA 内拾峰（与 QuantStudio 曲线一致）。
阶段 2：在同一 ROA 内对 F(T) 做 Boltzmann（两态 S 型）拟合得 TmB 与 B_fit（R²）。

Results 仅用于 Well→Sample Name 与过滤；不使用 Results 中 Tm 列。
可选对照：全区间仪器 Derivative 上全局两峰 Tm_inst_1/2。
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

try:
    from scipy import stats
    from scipy.optimize import curve_fit
    from scipy.signal import savgol_filter

    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False
    stats = None  # type: ignore[misc, assignment]
    curve_fit = None  # type: ignore[misc, assignment]
    savgol_filter = None  # type: ignore[misc, assignment]

HEADER_ROW = 45
MAX_ROAS = 6
PALETTE_HEX = ["#92B1D9", "#C1D8E9", "#DBDDEF", "#F6C8B6", "#D4D4D4"]

def _analysis_notes_zh(tmd_source: str) -> str:
    if tmd_source == "instrument_derivative":
        tmd_para = (
            "3. TmD：使用 Melt 表中的仪器 Derivative 列；若 --negate 则拾峰与作图用 -Derivative。"
            "在每个 ROA 内取 Derivative 的局部极大为 TmD；--sg-window 不参与仪器导数拾峰。\n"
        )
    else:
        tmd_para = (
            "3. TmD：对 Fluorescence 与 Temperature 使用 numpy.gradient 得 dF/dT；若 --negate，则拾峰与作图用 -dF/dT。"
            "若设置 --sg-window（奇数 ≥3），则先用 Savitzky–Golay 平滑荧光再求导（Boltzmann 仍用原始 F）。"
            "在每个 ROA 内取拾取用曲线的局部极大为 TmD；若无局部极大则退回为该 ROA 内全局最大点温度。\n"
        )
    return f"""\
IL-10 Protein thermal shift 自动分析说明（由 analyze_thermal_shift.py 生成）

1. 数据：Results / Melt Curve Raw Data 表头均在第 46 行。Results 仅 Well→Sample Name、Omit、对照过滤。

2. ROA：温度分析区间（Region of Analysis）。配置文件为 base-dir 下的 roa_config.json（可用 --roa-config 指定）。未提供配置时，每孔使用「该孔全部温度点」作为单个 ROA。可在 JSON 中为每块板或每个 Sample Name 配置多段 roas（至多 6 段）。

{tmd_para}
4. TmB / B_fit：在每个 ROA 内对原始 F(T) 拟合 F = fmin + (fmax-fmin)/(1+exp((tm-T)/s))；TmB=tm，B_fit=R²；每孔另给 TmB_ROAk_fit_se（tm 的渐近标准误）。需 scipy；若未安装则 TmB/B_fit/fit_se 为 NA。

5. 仪器 Derivative 全区间两峰：Tm_inst_1/2（间隔见 --peak-min-sep），作对照。

6. 汇总：按 plate + 完整 Sample Name；各 ROA 的 TmD/TmB/B_fit 及 Tm_inst 分别跨孔 mean±SD。对 TmD、TmB 的组均值另给 95% 置信区间（技术重复，t 分布，df 为非缺失孔数减 1）。

7. B_fit 为单孔 ROA 内 Boltzmann 的 R²；汇总含 B_fit 组内 min、median。每孔 TmB_ROA{{k}}_fit_se 为拟合参数 tm 的渐近标准误 sqrt(pcov[2,2])，边界解时需谨慎。

8. TmD 为导数拾峰，不定义与 Boltzmann 同义的 R²。

9. 图：上子图为 F–T；下子图为 TmD 所用信号；竖线标该组 TmD/TmB 各 ROA 均值。
"""


# 兼容旧引用；默认说明对应荧光数值导数模式
ANALYSIS_NOTES_ZH = _analysis_notes_zh("fluorescence")

DEFAULT_XLSX = [
    ("16mutants", "2026-04-28 IL10 16mutants zhengshizu.xlsx"),
    ("stability_test", "2026-04-28 IL10 stability test.xlsx"),
]


def _hex_to_rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def make_gradient_cmap() -> LinearSegmentedColormap:
    colors = [_hex_to_rgb(x) for x in PALETTE_HEX]
    return LinearSegmentedColormap.from_list("thermal_shift_palette", colors, N=256)


def is_control(sample_name: str) -> bool:
    s = str(sample_name).strip().lower()
    if not s or s == "nan":
        return True
    if s.startswith("h2o"):
        return True
    if "ntc" in s:
        return True
    if s in ("blank", "empty"):
        return True
    return False


def mutant_core(sample_name: str) -> str:
    return re.sub(r"(?i)-sybro.+$", "", str(sample_name).strip())


def safe_filename(s: str) -> str:
    x = re.sub(r"[^\w.\-]+", "_", str(s), flags=re.UNICODE)
    x = re.sub(r"_+", "_", x).strip("_")
    return x[:120] if len(x) > 120 else x


def read_results(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Results", header=HEADER_ROW)
    df["Well"] = pd.to_numeric(df["Well"], errors="coerce")
    df = df.dropna(subset=["Well"])
    df["Well"] = df["Well"].astype(int)
    return df


def read_melt(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Melt Curve Raw Data", header=HEADER_ROW)
    df["Well"] = pd.to_numeric(df["Well"], errors="coerce")
    df = df.dropna(subset=["Well"])
    df["Well"] = df["Well"].astype(int)
    return df


def filter_results(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "Omit" in out.columns:
        omit = out["Omit"].astype(str).str.lower().isin(("true", "1", "yes"))
        out = out.loc[~omit]
    out = out.loc[out["Sample Name"].notna()]
    out["Sample Name"] = out["Sample Name"].astype(str).str.strip()
    out = out.loc[~out["Sample Name"].apply(is_control)]
    return out


def load_roa_config(path: Path | None) -> dict:
    if path is None or not path.is_file():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def resolve_roas(
    plate_id: str,
    sample_name: str,
    cfg: dict,
    t_min: float,
    t_max: float,
) -> list[tuple[float, float]]:
    """返回 1..MAX_ROAS 个 (T_lo, T_hi) 闭区间内的分析区。"""
    plates = cfg.get("plates") or {}
    p = plates.get(plate_id) or {}
    samples = p.get("samples") or {}
    raw: list | None = None
    sn = str(sample_name).strip()
    if sn in samples and samples[sn].get("roas") is not None:
        raw = samples[sn]["roas"]
    elif p.get("roas") is not None:
        raw = p["roas"]
    elif (cfg.get("defaults") or {}).get("roas") is not None:
        raw = (cfg.get("defaults") or {}).get("roas")
    else:
        return [(float(t_min), float(t_max))]
    out: list[tuple[float, float]] = []
    for pair in raw[:MAX_ROAS]:
        if pair is None or len(pair) < 2:
            continue
        out.append((float(pair[0]), float(pair[1])))
    return out if out else [(float(t_min), float(t_max))]


def numerical_dfdt(temperature: np.ndarray, fluorescence: np.ndarray) -> np.ndarray:
    t = np.asarray(temperature, dtype=float)
    f = np.asarray(fluorescence, dtype=float)
    if t.size < 2:
        return np.full_like(t, np.nan)
    return np.gradient(f, t)


def smooth_fluorescence_sg(
    fluorescence: np.ndarray,
    window: int,
    polyorder: int,
) -> np.ndarray:
    """
    对荧光序列做 Savitzky–Golay 平滑后再求导，可抑制 dF/dT 毛刺。
    window 须为奇数且 >= polyorder+2；会自动裁剪到合法范围。
    """
    f = np.asarray(fluorescence, dtype=float)
    n = f.size
    if window is None or window < 3 or savgol_filter is None:
        return f
    w = int(window)
    if w % 2 == 0:
        w += 1
    w = min(w, n if n % 2 == 1 else n - 1)
    if w < 3 or w > n:
        return f
    po = int(polyorder)
    po = max(1, min(po, w - 2))
    try:
        return savgol_filter(f, w, po, mode="nearest").astype(float)
    except Exception:
        return f


def pick_peak_in_roa(
    temperature: np.ndarray,
    signal: np.ndarray,
    t_lo: float,
    t_hi: float,
) -> float:
    t = np.asarray(temperature, dtype=float)
    s = np.asarray(signal, dtype=float)
    m = (t >= t_lo) & (t <= t_hi)
    if int(np.sum(m)) < 1:
        return float("nan")
    Ti = t[m]
    Si = s[m]
    if Ti.size < 3:
        j = int(np.argmax(Si))
        return float(Ti[j])
    best_t, best_v = float("nan"), -np.inf
    for i in range(1, len(Si) - 1):
        if Si[i] >= Si[i - 1] and Si[i] >= Si[i + 1]:
            if float(Si[i]) > best_v:
                best_v = float(Si[i])
                best_t = float(Ti[i])
    if not np.isfinite(best_t):
        j = int(np.argmax(Si))
        best_t = float(Ti[j])
    return best_t


def boltzmann_4p(T: np.ndarray, fmin: float, fmax: float, tm: float, s: float) -> np.ndarray:
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return fmin + (fmax - fmin) / (1.0 + np.exp((tm - T) / s))


def fit_boltzmann_roa(
    temperature: np.ndarray,
    fluorescence: np.ndarray,
    t_lo: float,
    t_hi: float,
) -> tuple[float, float, float]:
    """返回 (TmB, R², TmB 拟合渐近标准误)。"""
    if not _HAS_SCIPY:
        return float("nan"), float("nan"), float("nan")
    t = np.asarray(temperature, dtype=float)
    f = np.asarray(fluorescence, dtype=float)
    m = (t >= t_lo) & (t <= t_hi)
    if int(np.sum(m)) < 8:
        return float("nan"), float("nan"), float("nan")
    tt = t[m]
    ff = f[m]
    lo, hi = float(np.min(ff)), float(np.max(ff))
    span = max(hi - lo, 1e-6)
    p0 = [lo + 0.05 * span, hi - 0.05 * span, float(np.median(tt)), max((tt[-1] - tt[0]) / 12.0, 0.3)]
    low = [lo - 3 * span, lo - 3 * span, float(tt.min()), 0.05]
    high = [hi + 3 * span, hi + 3 * span, float(tt.max()), max(float(tt.max() - tt.min()), 1.0)]
    if low[1] >= high[0]:
        low[1], high[0] = lo, hi + span
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            popt, pcov = curve_fit(
                boltzmann_4p,
                tt,
                ff,
                p0=p0,
                bounds=(low, high),
                maxfev=50000,
            )
        fmin, fmax, tm, s = (float(popt[0]), float(popt[1]), float(popt[2]), float(popt[3]))
        pred = boltzmann_4p(tt, fmin, fmax, tm, s)
        ss_res = float(np.sum((ff - pred) ** 2))
        ss_tot = float(np.sum((ff - np.mean(ff)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else float("nan")
        if not np.isfinite(r2):
            r2 = float("nan")
        try:
            tm_se = float(np.sqrt(float(pcov[2, 2])))
            if not np.isfinite(tm_se) or tm_se < 0:
                tm_se = float("nan")
        except Exception:
            tm_se = float("nan")
        return tm, r2, tm_se
    except Exception:
        return float("nan"), float("nan"), float("nan")


def pick_tm1_tm2_from_derivative(
    temperature: np.ndarray,
    derivative: np.ndarray,
    min_sep_c: float,
) -> tuple[float, float]:
    t = np.asarray(temperature, dtype=float)
    d = np.asarray(derivative, dtype=float)
    if t.size < 3 or d.size != t.size:
        return float("nan"), float("nan")
    interior = np.arange(1, len(d) - 1)
    is_locmax = (d[interior] > d[interior - 1]) & (d[interior] > d[interior + 1])
    idx = interior[is_locmax]
    if idx.size == 0:
        return float("nan"), float("nan")
    peak_t = t[idx]
    peak_d = d[idx]
    order = np.argsort(-peak_d)
    peak_t = peak_t[order]
    tm1 = float(peak_t[0])
    tm2 = float("nan")
    for j in range(1, len(peak_t)):
        if abs(float(peak_t[j]) - tm1) >= min_sep_c:
            tm2 = float(peak_t[j])
            break
    return tm1, tm2


def load_plate(xlsx_path: Path, plate_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw_r = read_results(xlsx_path)
    melt = read_melt(xlsx_path)
    fr = filter_results(raw_r)
    meta = fr[["Well", "Well Position", "Sample Name"]].drop_duplicates(subset=["Well"])
    melt_j = melt.merge(meta, on="Well", how="inner")
    melt_j["plate"] = plate_id
    return fr, melt_j


def compute_per_well_analysis(
    melt_all: pd.DataFrame,
    roa_cfg: dict,
    negate: bool,
    inst_min_sep: float,
    sg_window: int,
    sg_poly: int,
    tmd_source: str,
) -> tuple[pd.DataFrame, int]:
    """返回 per_well 表与 ROA 数量（列宽）。tmd_source: fluorescence | instrument_derivative。"""
    rows: list[dict] = []
    max_k = 1
    for (plate, well), g in melt_all.groupby(["plate", "Well"], sort=False):
        g = g.sort_values("Temperature")
        t = g["Temperature"].to_numpy(dtype=float)
        f = g["Fluorescence"].to_numpy(dtype=float)
        der_inst = g["Derivative"].to_numpy(dtype=float)
        sname = str(g["Sample Name"].iloc[0]).strip()
        roas = resolve_roas(plate, sname, roa_cfg, float(np.min(t)), float(np.max(t)))
        roas = roas[:MAX_ROAS]
        max_k = max(max_k, len(roas))

        if tmd_source == "instrument_derivative":
            d_pick = (-der_inst if negate else der_inst).astype(float)
        else:
            f_for_deriv = smooth_fluorescence_sg(f, sg_window, sg_poly) if sg_window >= 3 else f
            dfdt = numerical_dfdt(t, f_for_deriv)
            d_pick = (-dfdt if negate else dfdt).astype(float)

        row: dict = {"plate": plate, "Well": int(well), "Sample Name": sname}
        for i, (lo, hi) in enumerate(roas):
            k = i + 1
            row[f"TmD_ROA{k}"] = pick_peak_in_roa(t, d_pick, lo, hi)
            tmb, r2, tm_se = fit_boltzmann_roa(t, f, lo, hi)
            row[f"TmB_ROA{k}"] = tmb
            row[f"Bfit_ROA{k}"] = r2
            row[f"TmB_ROA{k}_fit_se"] = tm_se

        d_inst = -der_inst if negate else der_inst
        ti1, ti2 = pick_tm1_tm2_from_derivative(t, d_inst, inst_min_sep)
        row["Tm_inst_1"] = ti1
        row["Tm_inst_2"] = ti2
        rows.append(row)

    pw = pd.DataFrame(rows)
    for k in range(1, max_k + 1):
        for col in (f"TmD_ROA{k}", f"TmB_ROA{k}", f"Bfit_ROA{k}", f"TmB_ROA{k}_fit_se"):
            if col not in pw.columns:
                pw[col] = np.nan
    return pw, max_k


def _agg_mean_sd(series: pd.Series) -> tuple[float, float, str]:
    v = pd.to_numeric(series, errors="coerce").dropna()
    n = int(v.shape[0])
    if n == 0:
        return float("nan"), float("nan"), "NA"
    if n == 1:
        m = float(v.iloc[0])
        return m, float("nan"), f"{m:.2f} ± NA"
    m = float(v.mean())
    s = float(v.std(ddof=1))
    return m, s, f"{m:.2f} ± {s:.2f}"


def _mean_ci95(series: pd.Series) -> tuple[float, float, str]:
    """组内技术重复：均值 95% CI（t 分布）；n<2 或无数 scipy 时返回 NA。"""
    v = pd.to_numeric(series, errors="coerce").dropna()
    n = int(v.shape[0])
    if n < 2 or stats is None:
        return float("nan"), float("nan"), "NA"
    m = float(v.mean())
    s = float(v.std(ddof=1))
    se = s / np.sqrt(float(n))
    tcrit = float(stats.t.ppf(0.975, df=n - 1))
    lo = m - tcrit * se
    hi = m + tcrit * se
    return lo, hi, f"{lo:.2f}–{hi:.2f}"


def summarize_tm(per_well: pd.DataFrame, max_roas: int) -> pd.DataFrame:
    g = per_well.groupby(["plate", "Sample Name"], sort=False)
    rows = []
    for (plate, sname), sub in g:
        n = int(len(sub))
        base = {
            "plate": plate,
            "Sample_Name": sname,
            "mutant_core": mutant_core(sname),
            "n": n,
        }
        for k in range(1, max_roas + 1):
            m, s, st = _agg_mean_sd(sub[f"TmD_ROA{k}"])
            base[f"TmD_ROA{k}_mean"] = round(m, 4) if np.isfinite(m) else np.nan
            base[f"TmD_ROA{k}_sd"] = round(s, 4) if np.isfinite(s) else np.nan
            base[f"TmD_ROA{k}_mean_sd"] = st
            lo_d, hi_d, st_ci_d = _mean_ci95(sub[f"TmD_ROA{k}"])
            base[f"TmD_ROA{k}_ci95_low"] = round(lo_d, 4) if np.isfinite(lo_d) else np.nan
            base[f"TmD_ROA{k}_ci95_high"] = round(hi_d, 4) if np.isfinite(hi_d) else np.nan
            base[f"TmD_ROA{k}_mean_ci95"] = st_ci_d

            mb, sb, stb = _agg_mean_sd(sub[f"TmB_ROA{k}"])
            base[f"TmB_ROA{k}_mean"] = round(mb, 4) if np.isfinite(mb) else np.nan
            base[f"TmB_ROA{k}_sd"] = round(sb, 4) if np.isfinite(sb) else np.nan
            base[f"TmB_ROA{k}_mean_sd"] = stb
            lo_b, hi_b, st_ci_b = _mean_ci95(sub[f"TmB_ROA{k}"])
            base[f"TmB_ROA{k}_ci95_low"] = round(lo_b, 4) if np.isfinite(lo_b) else np.nan
            base[f"TmB_ROA{k}_ci95_high"] = round(hi_b, 4) if np.isfinite(hi_b) else np.nan
            base[f"TmB_ROA{k}_mean_ci95"] = st_ci_b

            br, sr, strr = _agg_mean_sd(sub[f"Bfit_ROA{k}"])
            base[f"Bfit_ROA{k}_mean"] = round(br, 4) if np.isfinite(br) else np.nan
            base[f"Bfit_ROA{k}_sd"] = round(sr, 4) if np.isfinite(sr) else np.nan
            base[f"Bfit_ROA{k}_mean_sd"] = strr
            bv = pd.to_numeric(sub[f"Bfit_ROA{k}"], errors="coerce").dropna()
            base[f"Bfit_ROA{k}_min"] = round(float(bv.min()), 4) if len(bv) else np.nan
            base[f"Bfit_ROA{k}_median"] = round(float(bv.median()), 4) if len(bv) else np.nan

            sev = pd.to_numeric(sub[f"TmB_ROA{k}_fit_se"], errors="coerce").dropna()
            base[f"TmB_ROA{k}_fit_se_mean"] = round(float(sev.mean()), 4) if len(sev) else np.nan
            base[f"TmB_ROA{k}_fit_se_median"] = round(float(sev.median()), 4) if len(sev) else np.nan
        m1, s1, t1 = _agg_mean_sd(sub["Tm_inst_1"])
        m2, s2, t2 = _agg_mean_sd(sub["Tm_inst_2"])
        base["Tm_inst_1_mean"] = round(m1, 4) if np.isfinite(m1) else np.nan
        base["Tm_inst_1_sd"] = round(s1, 4) if np.isfinite(s1) else np.nan
        base["Tm_inst_1_mean_sd"] = t1
        base["Tm_inst_2_mean"] = round(m2, 4) if np.isfinite(m2) else np.nan
        base["Tm_inst_2_sd"] = round(s2, 4) if np.isfinite(s2) else np.nan
        base["Tm_inst_2_mean_sd"] = t2
        rows.append(base)
    return pd.DataFrame(rows)


def curve_colors(n: int, cmap: LinearSegmentedColormap) -> list[tuple]:
    if n <= 0:
        return []
    if n == 1:
        return [cmap(0.5)]
    return [cmap(i / (n - 1)) for i in range(n)]


def plot_group(
    melt_j: pd.DataFrame,
    plate: str,
    sample_name: str,
    tm_d_means: list[float],
    tm_b_means: list[float],
    max_roas: int,
    out_path: Path,
    negate: bool,
    cmap: LinearSegmentedColormap,
    sg_window: int,
    sg_poly: int,
    tmd_source: str,
) -> None:
    sub = melt_j.loc[(melt_j["plate"] == plate) & (melt_j["Sample Name"] == sample_name)]
    wells = sorted(sub["Well"].unique())
    n_w = len(wells)
    colors = curve_colors(n_w, cmap)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 5.5), dpi=150, sharex=True)
    for idx, w in enumerate(wells):
        wdata = sub.loc[sub["Well"] == w].sort_values("Temperature")
        x = wdata["Temperature"].to_numpy(dtype=float)
        f = wdata["Fluorescence"].to_numpy(dtype=float)
        pos = wdata["Well Position"].iloc[0] if "Well Position" in wdata.columns else str(w)
        ax0.plot(x, f, color=colors[idx], linewidth=1.0, label=f"Well {w} ({pos})")
        if tmd_source == "instrument_derivative":
            der = wdata["Derivative"].to_numpy(dtype=float)
            y1 = -der if negate else der
            ax1.plot(x, y1, color=colors[idx], linewidth=1.0, label=f"Well {w} ({pos})")
        else:
            f_d = smooth_fluorescence_sg(f, sg_window, sg_poly) if sg_window >= 3 else f
            dfdt = numerical_dfdt(x, f_d)
            if negate:
                dfdt = -dfdt
            ax1.plot(x, dfdt, color=colors[idx], linewidth=1.0, label=f"Well {w} ({pos})")

    ax0.set_ylabel("Fluorescence (raw)")
    ax0.set_title(f"{plate}: {sample_name}")
    ax0.grid(True, alpha=0.25)
    ax0.legend(loc="best", fontsize=6)

    if tmd_source == "instrument_derivative":
        ylab = "-Derivative (instrument)" if negate else "Derivative (instrument)"
    elif sg_window >= 3:
        ylab = f"-dF/dT (SG{sg_window},p{sg_poly})" if negate else f"dF/dT (SG{sg_window},p{sg_poly})"
    else:
        ylab = "-dF/dT (numerical)" if negate else "dF/dT (numerical)"
    ax1.set_xlabel("Temperature (°C)")
    ax1.set_ylabel(ylab)
    ax1.grid(True, alpha=0.25)

    for k in range(max_roas):
        if k < len(tm_d_means) and np.isfinite(tm_d_means[k]):
            ax1.axvline(
                tm_d_means[k],
                color=plt.cm.tab10(k % 10),
                linestyle="--",
                linewidth=1,
                label=f"TmD ROA{k+1} mean={tm_d_means[k]:.2f}",
            )
        if k < len(tm_b_means) and np.isfinite(tm_b_means[k]):
            ax1.axvline(
                tm_b_means[k],
                color=plt.cm.tab10(k % 10),
                linestyle=":",
                linewidth=1.2,
                label=f"TmB ROA{k+1} mean={tm_b_means[k]:.2f}",
            )
    ax1.legend(loc="best", fontsize=6)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def run(
    base_dir: Path,
    out_dir: Path,
    negate: bool,
    inst_min_sep: float,
    roa_config_path: Path | None,
    sg_window: int,
    sg_poly: int,
    tmd_source: str,
    xlsx_specs: list[tuple[str, str]] | None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    cfg_path = roa_config_path
    if cfg_path is None:
        cand = base_dir / "roa_config.json"
        if cand.is_file():
            cfg_path = cand
    roa_cfg = load_roa_config(cfg_path)

    specs = xlsx_specs or DEFAULT_XLSX
    melt_parts: list[pd.DataFrame] = []
    for plate_id, fname in specs:
        path = base_dir / fname
        if not path.is_file():
            raise FileNotFoundError(path)
        _fr, melt_j = load_plate(path, plate_id)
        melt_parts.append(melt_j)

    melt_all = pd.concat(melt_parts, ignore_index=True)
    if sg_window >= 3 and not _HAS_SCIPY:
        print("Warning: --sg-window 需要 scipy；已忽略平滑，仍用原始 F 求导。")
        sg_window = 0
    per_well, max_roas = compute_per_well_analysis(
        melt_all, roa_cfg, negate, inst_min_sep, sg_window, sg_poly, tmd_source
    )

    summary = summarize_tm(per_well, max_roas)
    summary.to_csv(out_dir / "tm_summary.csv", index=False, encoding="utf-8-sig")
    per_well.to_csv(out_dir / "tm_per_well.csv", index=False, encoding="utf-8-sig")
    (out_dir / "analysis_notes.txt").write_text(_analysis_notes_zh(tmd_source), encoding="utf-8")

    if not _HAS_SCIPY:
        print("Warning: scipy not installed; TmB/B_fit columns are NA. pip install scipy")

    cmap = make_gradient_cmap()
    for _, row in summary.iterrows():
        plate = row["plate"]
        sname = row["Sample_Name"]
        tm_d_means = []
        tm_b_means = []
        for k in range(1, max_roas + 1):
            tm_d_means.append(float(row[f"TmD_ROA{k}_mean"]) if pd.notna(row.get(f"TmD_ROA{k}_mean")) else float("nan"))
            tm_b_means.append(float(row[f"TmB_ROA{k}_mean"]) if pd.notna(row.get(f"TmB_ROA{k}_mean")) else float("nan"))
        fn = f"melt_{plate}_{safe_filename(sname)}.png"
        plot_group(
            melt_all,
            plate,
            sname,
            tm_d_means,
            tm_b_means,
            max_roas,
            fig_dir / fn,
            negate=negate,
            cmap=cmap,
            sg_window=sg_window,
            sg_poly=sg_poly,
            tmd_source=tmd_source,
        )

    print(f"Wrote {out_dir / 'tm_summary.csv'} ({len(summary)} rows, {max_roas} ROA cols max, tmd_source={tmd_source})")
    print(f"Wrote {out_dir / 'tm_per_well.csv'} ({len(per_well)} wells)")
    print(f"Wrote {out_dir / 'analysis_notes.txt'}")
    print(f"Wrote figures to {fig_dir} ({len(summary)} PNG)")


def main() -> None:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description="PTS 风格：ROA + 数值 dF/dT 拾峰 (TmD) + Boltzmann (TmB)")
    p.add_argument("--base-dir", type=Path, default=here, help="含 xlsx 与可选 roa_config.json 的目录")
    p.add_argument("--out-dir", type=Path, default=here / "output_thermal_shift", help="输出目录")
    p.add_argument(
        "--negate",
        "--negate-derivative",
        action="store_true",
        dest="negate",
        help="对数值 dF/dT 与仪器 Derivative 取负后再拾峰/作图（主峰朝上）；--negate-derivative 为别名",
    )
    p.add_argument(
        "--peak-min-sep",
        type=float,
        default=3.0,
        metavar="DEG_C",
        help="仪器 Derivative 全区间两峰最小间隔（°C），默认 3",
    )
    p.add_argument(
        "--roa-config",
        type=Path,
        default=None,
        help="ROA JSON 路径；默认若 base-dir/roa_config.json 存在则自动读取",
    )
    p.add_argument(
        "--sg-window",
        type=int,
        default=0,
        metavar="N",
        help="Savitzky–Golay 窗口长度（奇数 ≥3）；0 表示不平滑荧光再求导。例如 15",
    )
    p.add_argument(
        "--sg-poly",
        type=int,
        default=2,
        help="SG 多项式阶数，须 < 窗口，默认 2",
    )
    p.add_argument(
        "--tmd-source",
        choices=("fluorescence", "instrument_derivative"),
        default="fluorescence",
        help="TmD 拾峰信号：荧光数值导数（默认）或仪器 Derivative 列",
    )
    p.add_argument(
        "--dual-output",
        action="store_true",
        help="连续跑两套：先写入 base-dir/output_thermal_shift_derivative（仪器 Derivative 拾 TmD），"
        "再写入 --out-dir（默认 output_thermal_shift，荧光数值导数拾 TmD）",
    )
    args = p.parse_args()
    sg_w = int(args.sg_window)
    if sg_w > 0 and sg_w % 2 == 0:
        sg_w += 1
    if args.dual_output:
        # 先写 derivative 目录，避免默认 output_thermal_shift 下文件被占用时整次失败
        run(
            args.base_dir,
            args.base_dir / "output_thermal_shift_derivative",
            args.negate,
            args.peak_min_sep,
            args.roa_config,
            sg_w,
            int(args.sg_poly),
            "instrument_derivative",
            None,
        )
        run(
            args.base_dir,
            args.out_dir,
            args.negate,
            args.peak_min_sep,
            args.roa_config,
            sg_w,
            int(args.sg_poly),
            "fluorescence",
            None,
        )
    else:
        run(
            args.base_dir,
            args.out_dir,
            args.negate,
            args.peak_min_sep,
            args.roa_config,
            sg_w,
            int(args.sg_poly),
            args.tmd_source,
            None,
        )


if __name__ == "__main__":
    main()

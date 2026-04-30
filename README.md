# Protein thermal shift analysis

**Description (for GitHub “About”):**  
PTS-style post-processing of QuantStudio / Applied Biosystems **Protein thermal shift** melt exports: **ROA**-local **TmD** from numerical **dF/dT** or instrument **Derivative**, **Boltzmann** **TmB** with **R²** and asymptotic **fit SE**, per-group **mean ± SD**, **95% CI**, and optional **dual-output** comparison; plus an optional **filtered TmD bar chart** for sybro2.5x-style summaries.

---

## Overview

This repository contains Python scripts to analyze **Melt Curve Raw Data** Excel exports (with **Results** metadata), aligned with a **PTS-like** workflow:

1. **TmD** — pick the dominant peak of a derivative signal inside each **Region of Analysis (ROA)**.  
   - Default: `numpy.gradient` of **Fluorescence** vs **Temperature** (optional Savitzky–Golay smoothing before differentiation).  
   - Alternative: instrument **`Derivative`** column (`--tmd-source instrument_derivative`).
2. **TmB** — four-parameter Boltzmann / logistic fit to **raw** fluorescence inside the same ROA; reports **R²** (`B_fit`) and **`tm` standard error** from `curve_fit` covariance.
3. **Tm_inst_1 / Tm_inst_2** — optional global two-peak reference on instrument derivative (minimum peak separation configurable).

**Results** sheet wells are used only for **Sample Name** and QC filtering; **instrument Tm from Results is not used** for these metrics.

---

## Repository layout

| Path | Role |
|------|------|
| `analyze_thermal_shift.py` | Main pipeline: read xlsx → per-well CSV → summary CSV → `analysis_notes.txt` → melt figures. |
| `plot_tmd_sybro25x_histogram.py` | Optional bar chart from `tm_summary.csv` with sybro2.5x / concentration filters (project-specific naming). |
| `roa_config.example.json` | Template for `roa_config.json` (plate- and sample-level ROA lists). |
| `docs/experiment_record_zh.md` | Full **experimental record** in Chinese: mathematics, ROA / TmD / TmB definitions, CLI, and reproducible commands (飞书-friendly). |

---

## Requirements

- Python 3.10+ recommended  
- Install: `pip install -r requirements.txt`  
- **SciPy** is required for Boltzmann fitting, SG smoothing, and **95% CI**; without it, TmB / B_fit / CI columns are missing or `NA`.

---

## Quick start

1. Place your **`.xlsx`** exports in a working folder (same folder as the scripts, or any `--base-dir`).
2. Edit **`DEFAULT_XLSX`** in `analyze_thermal_shift.py`: list of `(plate_id, filename)` tuples matching your files. The shipped example targets an **IL-10** study; replace with your plate IDs and filenames.
3. Optionally copy `roa_config.example.json` to `roa_config.json` in `--base-dir` and edit temperature windows.
4. Run from that folder:

```bash
pip install -r requirements.txt
python analyze_thermal_shift.py --base-dir . --out-dir ./output_thermal_shift
```

Useful options:

- `--tmd-source instrument_derivative` — TmD from instrument Derivative peaks.  
- `--dual-output` — run instrument-derivative TmD into `output_thermal_shift_derivative/` under `base-dir`, then fluorescence-based TmD into `--out-dir`.  
- `--negate` — flip derivative sign for peak picking / plotting.  
- `--sg-window 15 --sg-poly 2` — SG-smooth fluorescence before numerical derivative (Boltzmann still uses raw *F*).  
- `--out-dir` — always use a **new** directory if Excel might lock `tm_summary.csv`.

Optional bar chart (after you have `tm_summary.csv`):

```bash
python plot_tmd_sybro25x_histogram.py --summary ./path/to/tm_summary.csv --out ./figures/tmd_bars.png
```

Filter rules in that script are **name-based** (sybro2.5x, exclude `12*`/`12_`, WT excludes `4*`/`4_`); adjust the script if your naming differs.

---

## Outputs (under `--out-dir`)

- `tm_per_well.csv` — per well: `TmD_ROAk`, `TmB_ROAk`, `Bfit_ROAk`, `TmB_ROAk_fit_se`, `Tm_inst_*`, …  
- `tm_summary.csv` — grouped by plate + sample: means, SDs, 95% CI, B_fit min/median, etc.  
- `analysis_notes.txt` — auto-generated Chinese notes for the chosen `tmd_source`.  
- `figures/*.png` — dual-panel melt plots per summary row.

---

## Theory and parameter reference

See **[docs/experiment_record_zh.md](docs/experiment_record_zh.md)** for detailed **mathematics**, **ROA / TmD / TmB** definitions, and column glossary (same content as the lab’s Feishu-oriented experiment record).

---

## License

Scripts are shared as-is for research reuse; add a license file if you need a formal open-source terms.

## Maintainer

Repository: [github.com/YeyuSu/Protein-thermal-shift-analysis](https://github.com/YeyuSu/Protein-thermal-shift-analysis)

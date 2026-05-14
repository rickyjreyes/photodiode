#!/usr/bin/env python3
"""
compare_captures.py

Run the full spectral analysis pipeline on every .npy capture in a directory
and produce a side-by-side comparison table.

Usage
-----
python compare_captures.py --input-dir captures/20260309_010802 [other dirs...]
python compare_captures.py --files induced.npy virgin.npy disconnected.npy
python compare_captures.py --input-dir captures/ --sample-rate 20000 --nperm 500

Each file gets one row in the output table:
  label | n_peaks | exact_2_3_pairs | exact_koide_triplets | best_delta_chi2 | p_shuffle

Output
------
  comparison_results/
    comparison_table.csv
    comparison_table.txt   (human-readable)
    <label>_spectrum.png   (one per capture)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Re-use functions from the canonical analysis script.
sys.path.insert(0, str(Path(__file__).parent))
from scan_npy_logcos import (
    load_npy_waveform,
    compute_fft_spectrum,
    select_freq_window,
    detect_peaks,
    make_ratio_pairs,
    koide_triplets,
    build_log_domain,
    scan_logcos_cpu,
    shuffle_null,
    p_value_geq,
    CUPY_AVAILABLE,
)

try:
    from scan_npy_logcos import scan_logcos_gpu
    _HAS_GPU_SCAN = True
except ImportError:
    _HAS_GPU_SCAN = False

OUT_DIR = Path("comparison_results")


def analyze_one(path: Path, sample_rate_arg: float, fmin: float, fmax: float,
                peak_height_frac: float, top_n: int, kmin: float, kmax: float,
                nk: int, nperm: int, seed: int, use_gpu: bool) -> dict:
    t, y, sr = load_npy_waveform(path, sample_rate_arg)

    spec = compute_fft_spectrum(y, sr, use_cupy=use_gpu)
    spec_win = select_freq_window(spec, fmin, fmax)
    peaks = detect_peaks(spec_win, peak_height_frac, top_n)
    ratios = make_ratio_pairs(peaks)
    trips = koide_triplets(peaks)

    logdata = build_log_domain(spec_win)
    ell = logdata["ell_ln_f"].to_numpy(float)
    yy = logdata["y_detrended"].to_numpy(float)
    k_grid = np.linspace(kmin, kmax, nk)

    if use_gpu and _HAS_GPU_SCAN:
        scan_df, best = scan_logcos_gpu(ell, yy, k_grid)
    else:
        scan_df, best = scan_logcos_cpu(ell, yy, k_grid)

    p = None
    if nperm > 0:
        null_df = shuffle_null(ell, yy, k_grid, nperm, seed, use_gpu=use_gpu)
        p = p_value_geq(float(best["delta_chi2"]),
                        null_df["null_best_delta_chi2"].to_numpy(float))

    # Exact 2/3 pairs: delta_from_2_over_3 == 0.0
    exact_pairs = 0
    if len(ratios):
        exact_pairs = int((ratios["delta_from_2_over_3"] == 0.0).sum())

    # Exact Koide triplets: koide_error == 0.0
    exact_trips = 0
    if len(trips):
        exact_trips = int((trips["koide_error"] == 0.0).sum())

    return {
        "file": str(path),
        "samples": len(y),
        "sample_rate_Hz": sr,
        "duration_s": len(y) / sr,
        "n_peaks": len(peaks),
        "exact_2_3_pairs": exact_pairs,
        "exact_koide_triplets": exact_trips,
        "best_delta_chi2": round(float(best["delta_chi2"]), 4),
        "best_k": round(float(best["k"]), 4),
        "p_shuffle_scanmax": round(p, 6) if p is not None else None,
        "null_n": nperm,
        "_peaks": peaks,
        "_spec_win": spec_win,
    }


def plot_spectrum_for(label: str, spec_win: pd.DataFrame, peaks: pd.DataFrame,
                      out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(spec_win["frequency_Hz"], spec_win["amplitude"], lw=1)
    if len(peaks):
        ax.scatter(peaks["frequency_Hz"], peaks["amplitude"], s=30, zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("FFT amplitude")
    ax.set_title(f"Spectrum — {label}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{label}_spectrum.png", dpi=150)
    plt.close(fig)


def collect_files(args) -> list[Path]:
    files: list[Path] = []
    if args.files:
        for f in args.files:
            p = Path(f)
            if not p.exists():
                print(f"[warn] file not found: {p}")
            else:
                files.append(p)
    if args.input_dir:
        for d in args.input_dir:
            dp = Path(d)
            if dp.is_dir():
                files.extend(sorted(dp.rglob("*.npy")))
            elif dp.suffix == ".npy" and dp.exists():
                files.append(dp)
            else:
                print(f"[warn] not a directory or .npy file: {dp}")
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", nargs="+", help="Directory/directories of .npy files")
    ap.add_argument("--files", nargs="+", help="Explicit list of .npy files")
    ap.add_argument("--sample-rate", type=float, default=20000.0)
    ap.add_argument("--fmin", type=float, default=20.0)
    ap.add_argument("--fmax", type=float, default=1000.0)
    ap.add_argument("--peak-height-frac", type=float, default=0.05)
    ap.add_argument("--top-peaks", type=int, default=30)
    ap.add_argument("--kmin", type=float, default=0.5)
    ap.add_argument("--kmax", type=float, default=80.0)
    ap.add_argument("--nk", type=int, default=4000)
    ap.add_argument("--nperm", type=int, default=500,
                    help="Shuffle null trials per file (0 = skip)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out-dir", default="comparison_results")
    args = ap.parse_args()

    use_gpu = (not args.cpu) and CUPY_AVAILABLE
    OUT_DIR_local = Path(args.out_dir)
    OUT_DIR_local.mkdir(exist_ok=True, parents=True)

    files = collect_files(args)
    if not files:
        print("No .npy files found. Use --input-dir or --files.")
        sys.exit(1)

    print(f"[compare] {len(files)} file(s), nperm={args.nperm}, gpu={use_gpu}")
    print()

    rows = []
    for path in files:
        label = path.stem
        print(f"[analyze] {path}")
        try:
            result = analyze_one(
                path,
                sample_rate_arg=args.sample_rate,
                fmin=args.fmin, fmax=args.fmax,
                peak_height_frac=args.peak_height_frac,
                top_n=args.top_peaks,
                kmin=args.kmin, kmax=args.kmax, nk=args.nk,
                nperm=args.nperm, seed=args.seed,
                use_gpu=use_gpu,
            )
            plot_spectrum_for(label, result.pop("_spec_win"), result.pop("_peaks"),
                              OUT_DIR_local)
            result["label"] = label
            rows.append(result)
            print(f"  peaks={result['n_peaks']}  exact_pairs={result['exact_2_3_pairs']}"
                  f"  exact_triplets={result['exact_koide_triplets']}"
                  f"  delta_chi2={result['best_delta_chi2']}"
                  f"  p={result['p_shuffle_scanmax']}")
        except Exception as e:
            print(f"  [error] {e}")
            rows.append({"label": label, "file": str(path), "error": str(e)})
        print()

    display_cols = [
        "label", "file", "n_peaks", "exact_2_3_pairs", "exact_koide_triplets",
        "best_delta_chi2", "best_k", "p_shuffle_scanmax", "null_n",
        "samples", "sample_rate_Hz", "duration_s",
    ]

    df = pd.DataFrame(rows)
    df = df.reindex(columns=[c for c in display_cols if c in df.columns])

    csv_path = OUT_DIR_local / "comparison_table.csv"
    df.to_csv(csv_path, index=False)

    txt_path = OUT_DIR_local / "comparison_table.txt"
    with open(txt_path, "w") as f:
        f.write(df.to_string(index=False))
        f.write("\n")

    print("=" * 80)
    print("COMPARISON TABLE")
    print("=" * 80)
    print(df.to_string(index=False))
    print()
    print(f"Saved: {csv_path}")
    print(f"Saved: {txt_path}")
    print()
    print("Interpretation guide:")
    print("  exact_2_3_pairs       — pairs with ratio = 2/3 exactly (expect 2 in induced)")
    print("  exact_koide_triplets  — triplets with koide_error = 0 (expect 2 in induced)")
    print("  best_delta_chi2       — log-cos scan max (higher = more structure)")
    print("  p_shuffle_scanmax     — fraction of shuffled spectra that exceeded real")
    print("  Controls should score 0/0 on pairs/triplets and higher p than induced.")


if __name__ == "__main__":
    main()

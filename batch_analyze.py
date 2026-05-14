#!/usr/bin/env python3
"""
batch_analyze.py  —  Blinded batch analyzer

Processes every .npy file in a directory, assigns each a numeric ID,
and writes the full spectral results WITHOUT condition labels.

The researcher assigns condition labels (induced / virgin / sham /
disconnected / etc.) to the output CSV only after inspecting the results,
preventing unconscious parameter tuning per condition.

Workflow
--------
1. Place all captures in a single directory (any mix of conditions):

       captures_blind/
         001.npy
         002.npy
         003.npy
         ...  (name them however you like — they are sorted alphabetically)

2. Run:
       python batch_analyze.py --input-dir captures_blind/

3. Inspect:
       batch_results/results.csv   ← numeric IDs only, no labels

4. Add a "condition" column to results.csv mapping each ID to its condition
   (induced / virgin / sham / disconnected / etc.).

5. Re-read to assess whether the induced capture is distinguishable.

Output
------
  batch_results/
    results.csv             — one row per file (ID, filename, metrics)
    results_labeled.csv     — written after --label-map is provided
    <id>_spectrum.png       — spectrum per file
    summary.json
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

OUT_DIR = Path("batch_results")


def analyze_one(path: Path, sample_rate_arg: float,
                fmin: float, fmax: float,
                peak_height_frac: float, top_n: int,
                kmin: float, kmax: float, nk: int,
                nperm: int, seed: int, use_gpu: bool) -> dict:
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
        _, best = scan_logcos_gpu(ell, yy, k_grid)
    else:
        _, best = scan_logcos_cpu(ell, yy, k_grid)

    p = None
    if nperm > 0:
        null_df = shuffle_null(ell, yy, k_grid, nperm, seed, use_gpu=use_gpu)
        p = p_value_geq(float(best["delta_chi2"]),
                        null_df["null_best_delta_chi2"].to_numpy(float))

    exact_pairs = int((ratios["delta_from_2_over_3"] == 0.0).sum()) if len(ratios) else 0
    exact_trips = int((trips["koide_error"] == 0.0).sum()) if len(trips) else 0

    # Peak frequencies as comma-separated string for reference.
    peak_freqs = ",".join(f"{f:.0f}" for f in
                          peaks["frequency_Hz"].sort_values().tolist()) if len(peaks) else ""

    return {
        "n_peaks": len(peaks),
        "peak_freqs_Hz": peak_freqs,
        "exact_2_3_pairs": exact_pairs,
        "exact_koide_triplets": exact_trips,
        "best_delta_chi2": round(float(best["delta_chi2"]), 4),
        "best_k": round(float(best["k"]), 4),
        "p_shuffle_scanmax": round(p, 6) if p is not None else None,
        "null_n": nperm,
        "samples": len(y),
        "sample_rate_Hz": sr,
        "duration_s": round(len(y) / sr, 4),
        "_peaks": peaks,
        "_spec_win": spec_win,
    }


def plot_spectrum(file_id: str, spec_win: pd.DataFrame, peaks: pd.DataFrame,
                 out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(spec_win["frequency_Hz"], spec_win["amplitude"], lw=1)
    if len(peaks):
        ax.scatter(peaks["frequency_Hz"], peaks["amplitude"], s=30, zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("FFT amplitude")
    ax.set_title(f"Spectrum — ID {file_id}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{file_id}_spectrum.png", dpi=150)
    plt.close(fig)


def apply_label_map(df: pd.DataFrame, label_map_path: Path) -> pd.DataFrame:
    """
    label_map is a CSV with columns: id, condition
    Example:
      id,condition
      001,induced
      002,virgin
      003,disconnected
    """
    lm = pd.read_csv(label_map_path, dtype=str)
    lm.columns = [c.strip() for c in lm.columns]
    lm = lm.rename(columns={"id": "id"})
    merged = df.merge(lm, on="id", how="left")
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True,
                    help="Directory containing .npy files")
    ap.add_argument("--label-map",
                    help="CSV mapping file IDs to conditions (id,condition). "
                         "Produces results_labeled.csv when provided.")
    ap.add_argument("--sample-rate", type=float, default=20000.0)
    ap.add_argument("--fmin", type=float, default=20.0)
    ap.add_argument("--fmax", type=float, default=1000.0)
    ap.add_argument("--peak-height-frac", type=float, default=0.05)
    ap.add_argument("--top-peaks", type=int, default=30)
    ap.add_argument("--kmin", type=float, default=0.5)
    ap.add_argument("--kmax", type=float, default=80.0)
    ap.add_argument("--nk", type=int, default=4000)
    ap.add_argument("--nperm", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out-dir", default="batch_results")
    args = ap.parse_args()

    use_gpu = (not args.cpu) and CUPY_AVAILABLE
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    input_dir = Path(args.input_dir)
    files = sorted(input_dir.rglob("*.npy"))
    if not files:
        print(f"No .npy files found in {input_dir}")
        sys.exit(1)

    # Assign zero-padded numeric IDs so results are not connected to filenames
    # during blind review.
    id_width = max(3, len(str(len(files))))

    rows = []
    print(f"[batch] {len(files)} file(s), nperm={args.nperm}, gpu={use_gpu}")
    print()

    for i, path in enumerate(files):
        file_id = str(i + 1).zfill(id_width)
        print(f"[{file_id}] {path.name}")
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
            plot_spectrum(file_id, result.pop("_spec_win"), result.pop("_peaks"), out_dir)
            result["id"] = file_id
            result["filename"] = path.name
            rows.append(result)
            print(f"  peaks={result['n_peaks']}  pairs={result['exact_2_3_pairs']}"
                  f"  triplets={result['exact_koide_triplets']}"
                  f"  delta_chi2={result['best_delta_chi2']}"
                  f"  p={result['p_shuffle_scanmax']}")
        except Exception as e:
            print(f"  [error] {e}")
            rows.append({"id": file_id, "filename": path.name, "error": str(e)})
        print()

    display_cols = [
        "id", "filename", "n_peaks", "peak_freqs_Hz",
        "exact_2_3_pairs", "exact_koide_triplets",
        "best_delta_chi2", "best_k", "p_shuffle_scanmax", "null_n",
        "samples", "sample_rate_Hz", "duration_s",
    ]

    df = pd.DataFrame(rows)
    df = df.reindex(columns=[c for c in display_cols if c in df.columns])

    csv_path = out_dir / "results.csv"
    df.to_csv(csv_path, index=False)
    print(f"[saved] {csv_path}")

    if args.label_map:
        labeled = apply_label_map(df, Path(args.label_map))
        labeled_path = out_dir / "results_labeled.csv"
        labeled.to_csv(labeled_path, index=False)
        print(f"[saved] {labeled_path}")
        print()
        print(labeled.to_string(index=False))
    else:
        print()
        print(df.to_string(index=False))
        print()
        print("Next step: add condition labels.")
        print("Create a CSV file with columns: id, condition")
        print("Then re-run with: --label-map your_labels.csv")

    with open(out_dir / "summary.json", "w") as f:
        json.dump({
            "n_files": len(rows),
            "parameters": {
                "fmin": args.fmin, "fmax": args.fmax,
                "kmin": args.kmin, "kmax": args.kmax, "nk": args.nk,
                "nperm": args.nperm, "seed": args.seed,
                "peak_height_frac": args.peak_height_frac,
            },
            "rows": [
                {k: v for k, v in r.items() if not k.startswith("_")}
                for r in rows
            ],
        }, f, indent=2)


if __name__ == "__main__":
    main()

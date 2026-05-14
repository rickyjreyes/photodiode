#!/usr/bin/env python3
"""
persistence_scan.py

Capture waveforms at fixed time intervals after optical excitation and track
whether spectral structure persists or decays.

At each timepoint the script:
  1. Captures a waveform from the oscilloscope.
  2. Runs the full FFT + Q=2/3 + log-cos analysis.
  3. Appends one row to a running time-series CSV.
  4. Saves the raw .npy waveform with a timestamp label.

After all captures it plots the persistence curve:
  Delta chi2 and exact-pair count vs. time since excitation.

Usage
-----
# 10 captures, one every 60 seconds, starting immediately after excitation:
python persistence_scan.py --interval 60 --count 10

# More captures with GPU-accelerated scan and 500 shuffle trials:
python persistence_scan.py --interval 30 --count 20 --nperm 500

# Analyze previously saved .npy files in time order (no live capture):
python persistence_scan.py --replay-dir persistence_captures/ --sample-rate 20000

Output
------
  persistence_captures/
    t0000_<epoch>.npy
    t0001_<epoch>.npy
    ...
  persistence_results/
    timeseries.csv
    persistence_curve.png
    summary.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
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

CAP_DIR = Path("persistence_captures")
RES_DIR = Path("persistence_results")


# ============================================================
# Oscilloscope capture
# ============================================================

def capture_waveform(scope_ip: str, channel: str) -> tuple[np.ndarray, float]:
    """Capture one waveform from a Siglent SDS scope. Returns (volts, sample_rate)."""
    import pyvisa
    import re

    rm = pyvisa.ResourceManager("@py")
    scope = rm.open_resource(f"TCPIP::{scope_ip}::INSTR")
    scope.timeout = 30000

    scope.write(":STOP")

    def _num(cmd):
        resp = scope.query(cmd)
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", resp.split()[-1])
        return float(m.group(0)) if m else 0.0

    vscale = _num(f"{channel}:VDIV?")
    voffset = _num(f"{channel}:OFST?")
    sample_rate = _num(":ACQ:SRAT?")

    scope.write("WFSU SP,0,NP,0,FP,0")
    scope.write("WAV:MODE RAW")
    scope.write(f"WAV:SOURCE {channel}")

    data = scope.query_binary_values(":WAV:DATA?", datatype="B")
    scope.close()
    rm.close()

    wave = np.array(data, dtype=float)
    volts = (wave - 128.0) * (vscale / 25.0) + voffset
    volts -= volts.mean()

    return volts, sample_rate


# ============================================================
# Single-capture analysis
# ============================================================

def analyze_waveform(volts: np.ndarray, sample_rate: float,
                     fmin: float, fmax: float,
                     peak_height_frac: float, top_n: int,
                     kmin: float, kmax: float, nk: int,
                     nperm: int, seed: int, use_gpu: bool) -> dict:
    spec = compute_fft_spectrum(volts, sample_rate, use_cupy=use_gpu)
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

    return {
        "n_peaks": len(peaks),
        "exact_2_3_pairs": exact_pairs,
        "exact_koide_triplets": exact_trips,
        "best_delta_chi2": float(best["delta_chi2"]),
        "best_k": float(best["k"]),
        "p_shuffle_scanmax": p,
    }


# ============================================================
# Replay mode (analyze saved files)
# ============================================================

def replay_mode(replay_dir: Path, sample_rate_arg: float,
                fmin: float, fmax: float,
                peak_height_frac: float, top_n: int,
                kmin: float, kmax: float, nk: int,
                nperm: int, seed: int, use_gpu: bool,
                out_dir: Path) -> list[dict]:
    files = sorted(replay_dir.rglob("*.npy"))
    if not files:
        print(f"[replay] no .npy files in {replay_dir}")
        return []

    print(f"[replay] {len(files)} file(s)")
    rows = []
    t0 = None

    for i, path in enumerate(files):
        epoch = None
        # Try to extract epoch from filename (e.g. t0003_1773046194.npy)
        parts = path.stem.split("_")
        for part in parts:
            try:
                v = int(part)
                if v > 1_000_000_000:
                    epoch = v
            except ValueError:
                pass
        if epoch is None:
            epoch = int(path.stat().st_mtime)

        if t0 is None:
            t0 = epoch
        elapsed = epoch - t0

        print(f"  [{i+1}/{len(files)}] {path.name}  elapsed={elapsed}s")
        t_arr, y, sr = load_npy_waveform(path, sample_rate_arg)
        result = analyze_waveform(
            y, sr, fmin, fmax, peak_height_frac, top_n,
            kmin, kmax, nk, nperm, seed, use_gpu,
        )
        result.update({"index": i, "epoch": epoch, "elapsed_s": elapsed, "file": str(path)})
        rows.append(result)
        print(f"    peaks={result['n_peaks']}  pairs={result['exact_2_3_pairs']}"
              f"  triplets={result['exact_koide_triplets']}"
              f"  delta_chi2={result['best_delta_chi2']:.4f}"
              f"  p={result['p_shuffle_scanmax']}")

    return rows


# ============================================================
# Persistence plot
# ============================================================

def plot_persistence(rows: list[dict], out_dir: Path) -> None:
    df = pd.DataFrame(rows)
    if "elapsed_s" not in df.columns or df.empty:
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    ax1.plot(df["elapsed_s"], df["best_delta_chi2"], "o-", lw=1.5)
    ax1.set_ylabel("Best Δχ² (log-cos)")
    ax1.set_title("Spectral structure persistence after excitation")
    ax1.grid(True, alpha=0.3)

    ax2.plot(df["elapsed_s"], df["exact_2_3_pairs"], "s-", lw=1.5, label="exact 2/3 pairs")
    ax2.plot(df["elapsed_s"], df["exact_koide_triplets"], "^-", lw=1.5, label="exact Koide triplets")
    ax2.set_xlabel("Time since excitation [s]")
    ax2.set_ylabel("Count")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "persistence_curve.png", dpi=150)
    plt.close(fig)
    print(f"[plot] {out_dir / 'persistence_curve.png'}")


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser()

    # Live capture mode
    cap_group = ap.add_argument_group("Live capture")
    cap_group.add_argument("--interval", type=float, default=60.0,
                           help="Seconds between captures")
    cap_group.add_argument("--count", type=int, default=10,
                           help="Number of captures to take")
    cap_group.add_argument("--channel", default="C3")
    cap_group.add_argument("--cap-dir", default="persistence_captures")

    # Replay mode
    ap.add_argument("--replay-dir",
                    help="Analyze previously saved .npy files instead of live capture")

    # Analysis parameters (match canonical run defaults)
    ap.add_argument("--sample-rate", type=float, default=20000.0)
    ap.add_argument("--fmin", type=float, default=20.0)
    ap.add_argument("--fmax", type=float, default=1000.0)
    ap.add_argument("--peak-height-frac", type=float, default=0.05)
    ap.add_argument("--top-peaks", type=int, default=30)
    ap.add_argument("--kmin", type=float, default=0.5)
    ap.add_argument("--kmax", type=float, default=80.0)
    ap.add_argument("--nk", type=int, default=4000)
    ap.add_argument("--nperm", type=int, default=200,
                    help="Shuffle null trials per capture (0 = skip; lower for speed)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out-dir", default="persistence_results")

    args = ap.parse_args()

    use_gpu = (not args.cpu) and CUPY_AVAILABLE
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    rows: list[dict] = []

    if args.replay_dir:
        rows = replay_mode(
            Path(args.replay_dir),
            sample_rate_arg=args.sample_rate,
            fmin=args.fmin, fmax=args.fmax,
            peak_height_frac=args.peak_height_frac,
            top_n=args.top_peaks,
            kmin=args.kmin, kmax=args.kmax, nk=args.nk,
            nperm=args.nperm, seed=args.seed,
            use_gpu=use_gpu,
            out_dir=out_dir,
        )

    else:
        from dotenv import load_dotenv
        load_dotenv()
        scope_ip = os.environ.get("SCOPE_IP", "")
        if not scope_ip:
            print("[error] SCOPE_IP not set in environment or .env file.")
            sys.exit(1)

        cap_dir = Path(args.cap_dir)
        cap_dir.mkdir(exist_ok=True, parents=True)

        t_excitation = time.time()
        print(f"[start] excitation epoch = {int(t_excitation)}")
        print(f"[plan]  {args.count} captures, interval = {args.interval}s")
        print()

        for i in range(args.count):
            epoch = int(time.time())
            elapsed = epoch - int(t_excitation)
            label = f"t{i:04d}_{epoch}"
            print(f"[capture {i+1}/{args.count}] elapsed={elapsed}s  label={label}")

            volts, sr = capture_waveform(scope_ip, args.channel)
            npy_path = cap_dir / f"{label}.npy"
            np.save(npy_path, volts)
            print(f"  saved {npy_path}  samples={len(volts)}  sr={sr}")

            result = analyze_waveform(
                volts, sr,
                fmin=args.fmin, fmax=args.fmax,
                peak_height_frac=args.peak_height_frac,
                top_n=args.top_peaks,
                kmin=args.kmin, kmax=args.kmax, nk=args.nk,
                nperm=args.nperm, seed=args.seed,
                use_gpu=use_gpu,
            )
            result.update({
                "index": i,
                "epoch": epoch,
                "elapsed_s": elapsed,
                "file": str(npy_path),
            })
            rows.append(result)
            print(f"  peaks={result['n_peaks']}  pairs={result['exact_2_3_pairs']}"
                  f"  triplets={result['exact_koide_triplets']}"
                  f"  delta_chi2={result['best_delta_chi2']:.4f}"
                  f"  p={result['p_shuffle_scanmax']}")

            if i < args.count - 1:
                next_at = t_excitation + (i + 1) * args.interval
                wait = max(0.0, next_at - time.time())
                print(f"  waiting {wait:.1f}s until next capture...")
                time.sleep(wait)
            print()

    if not rows:
        print("[done] no rows collected.")
        return

    df = pd.DataFrame(rows)
    csv_path = out_dir / "timeseries.csv"
    df.to_csv(csv_path, index=False)
    print(f"[saved] {csv_path}")

    plot_persistence(rows, out_dir)

    summary = {
        "n_captures": len(rows),
        "elapsed_range_s": [
            int(df["elapsed_s"].min()),
            int(df["elapsed_s"].max()),
        ] if "elapsed_s" in df.columns else None,
        "delta_chi2_first": float(df["best_delta_chi2"].iloc[0]),
        "delta_chi2_last": float(df["best_delta_chi2"].iloc[-1]),
        "exact_pairs_first": int(df["exact_2_3_pairs"].iloc[0]),
        "exact_pairs_last": int(df["exact_2_3_pairs"].iloc[-1]),
    }
    with open(out_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[saved] {out_dir / 'summary.json'}")

    print()
    print(df[["elapsed_s", "n_peaks", "exact_2_3_pairs",
              "exact_koide_triplets", "best_delta_chi2",
              "p_shuffle_scanmax"]].to_string(index=False))


if __name__ == "__main__":
    main()

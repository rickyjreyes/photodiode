#!/usr/bin/env python3
"""
block_bootstrap.py  —  Block-bootstrap log-cos significance test

The standard null in scan_npy_logcos.py shuffles FFT amplitudes.
That destroys all phase structure but preserves the marginal amplitude
distribution. A stricter null shuffles contiguous blocks of the raw
time-domain waveform. This preserves short-range temporal correlations
(e.g., 60 Hz oscillations within a block) while destroying long-range
structure. If the log-cos result survives this null, the claim is stronger.

Two block sizes are tested:
  - SHORT:  half a 60 Hz cycle (~83 samples at 20 kHz)
  - LONG:   one full 60 Hz cycle (~333 samples at 20 kHz)

For each block size:
  1. Divide the waveform into non-overlapping blocks of that size.
  2. Shuffle the block order.
  3. Recompute FFT → log-cos scan.
  4. Record the scan-max Δχ².
  5. Repeat nperm times to build the null distribution.
  6. Report p = fraction of null trials exceeding the real scan-max.

Usage
-----
python block_bootstrap.py \\
    --input captures/20260309_010802/waveform_raw.npy \\
    --sample-rate 20000 --nperm 1000

Output
------
  block_bootstrap_results/
    null_short_block.csv
    null_long_block.csv
    null_amplitude_shuffle.csv   (standard null, for comparison)
    summary.json
    null_distributions.png
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
    build_log_domain,
    scan_logcos_cpu,
    p_value_geq,
    CUPY_AVAILABLE,
)

try:
    from scan_npy_logcos import scan_logcos_gpu
    _HAS_GPU_SCAN = True
except ImportError:
    _HAS_GPU_SCAN = False

OUT_DIR = Path("block_bootstrap_results")


def logcos_scan_max(y_waveform: np.ndarray, sample_rate: float,
                    fmin: float, fmax: float,
                    k_grid: np.ndarray, use_gpu: bool) -> float:
    """Compute the scan-max Δχ² for a given waveform."""
    spec = compute_fft_spectrum(y_waveform, sample_rate, use_cupy=False)
    spec_win = select_freq_window(spec, fmin, fmax)
    logdata = build_log_domain(spec_win)
    ell = logdata["ell_ln_f"].to_numpy(float)
    yy = logdata["y_detrended"].to_numpy(float)

    if use_gpu and _HAS_GPU_SCAN:
        from scan_npy_logcos import scan_logcos_gpu
        _, best = scan_logcos_gpu(ell, yy, k_grid)
    else:
        _, best = scan_logcos_cpu(ell, yy, k_grid)

    return float(best["delta_chi2"])


def block_shuffle(y: np.ndarray, block_size: int, rng: np.random.Generator) -> np.ndarray:
    """Shuffle non-overlapping blocks of length block_size."""
    n_blocks = len(y) // block_size
    usable = n_blocks * block_size
    blocks = y[:usable].reshape(n_blocks, block_size).copy()
    idx = rng.permutation(n_blocks)
    shuffled = blocks[idx].ravel()
    # Append any remainder unchanged (short tail, typically <1 block).
    remainder = y[usable:]
    return np.concatenate([shuffled, remainder])


def amplitude_shuffle(y: np.ndarray, sample_rate: float,
                      fmin: float, fmax: float,
                      k_grid: np.ndarray, use_gpu: bool,
                      rng: np.random.Generator) -> float:
    """Standard amplitude-shuffle null: permute FFT amplitudes, not waveform."""
    spec = compute_fft_spectrum(y, sample_rate, use_cupy=False)
    spec_win = select_freq_window(spec, fmin, fmax)
    logdata = build_log_domain(spec_win)
    ell = logdata["ell_ln_f"].to_numpy(float)
    yy = logdata["y_detrended"].to_numpy(float)

    yy_perm = rng.permutation(yy)

    if use_gpu and _HAS_GPU_SCAN:
        _, best = scan_logcos_gpu(ell, yy_perm, k_grid)
    else:
        _, best = scan_logcos_cpu(ell, yy_perm, k_grid)

    return float(best["delta_chi2"])


def run_null(null_type: str, y: np.ndarray, sample_rate: float,
             block_size: int, fmin: float, fmax: float,
             k_grid: np.ndarray, nperm: int, seed: int,
             use_gpu: bool) -> np.ndarray:
    rng = np.random.default_rng(seed)
    null_vals = np.empty(nperm, dtype=float)
    report_every = max(1, nperm // 10)

    for i in range(nperm):
        if null_type == "amplitude":
            null_vals[i] = amplitude_shuffle(y, sample_rate, fmin, fmax,
                                             k_grid, use_gpu, rng)
        else:
            y_shuf = block_shuffle(y, block_size, rng)
            null_vals[i] = logcos_scan_max(y_shuf, sample_rate, fmin, fmax,
                                           k_grid, use_gpu)

        if (i + 1) % report_every == 0:
            print(f"  [{null_type}] {i+1}/{nperm}  running max={null_vals[:i+1].max():.4f}")

    return null_vals


def plot_nulls(real_stat: float, results: dict[str, np.ndarray], out_dir: Path) -> None:
    colors = {
        "amplitude_shuffle": "#1f77b4",
        "block_short": "#ff7f0e",
        "block_long": "#2ca02c",
    }
    labels = {
        "amplitude_shuffle": "Amplitude shuffle (standard)",
        "block_short": "Block shuffle — short block",
        "block_long": "Block shuffle — long block",
    }

    fig, ax = plt.subplots(figsize=(10, 5))
    for key, vals in results.items():
        ax.hist(vals, bins=40, alpha=0.6,
                color=colors.get(key, "gray"),
                label=labels.get(key, key), density=True)

    ax.axvline(real_stat, color="red", lw=2, ls="--",
               label=f"Real scan-max Δχ² = {real_stat:.3f}")
    ax.set_xlabel("Scan-max Δχ²")
    ax.set_ylabel("Density")
    ax.set_title("Block-bootstrap null distributions vs. real scan-max")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "null_distributions.png", dpi=150)
    plt.close(fig)
    print(f"[plot] {out_dir / 'null_distributions.png'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--sample-rate", type=float, default=20000.0)
    ap.add_argument("--fmin", type=float, default=20.0)
    ap.add_argument("--fmax", type=float, default=1000.0)
    ap.add_argument("--kmin", type=float, default=0.5)
    ap.add_argument("--kmax", type=float, default=80.0)
    ap.add_argument("--nk", type=int, default=4000)
    ap.add_argument("--nperm", type=int, default=1000,
                    help="Null trials per method (total trials = 3 × nperm)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--mains-hz", type=float, default=60.0,
                    help="Mains frequency in Hz (determines block sizes)")
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out-dir", default="block_bootstrap_results")
    args = ap.parse_args()

    use_gpu = (not args.cpu) and CUPY_AVAILABLE
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    print(f"[load] {args.input}")
    t, y, sample_rate = load_npy_waveform(Path(args.input), args.sample_rate)
    print(f"[data] samples={len(y):,}  sr={sample_rate:.0f} Hz  "
          f"duration={len(y)/sample_rate:.3f}s")

    k_grid = np.linspace(args.kmin, args.kmax, args.nk)

    # Real scan-max
    print("[scan] computing real scan-max...")
    real_stat = logcos_scan_max(y, sample_rate, args.fmin, args.fmax,
                                k_grid, use_gpu)
    print(f"[scan] real Δχ² = {real_stat:.6f}")
    print()

    # Block sizes: half-cycle and full-cycle of mains frequency
    samples_per_mains_cycle = int(round(sample_rate / args.mains_hz))
    block_short = max(2, samples_per_mains_cycle // 2)   # half-cycle
    block_long = samples_per_mains_cycle                  # full cycle

    print(f"[blocks] mains={args.mains_hz} Hz  "
          f"short={block_short} samples ({1000*block_short/sample_rate:.1f} ms)  "
          f"long={block_long} samples ({1000*block_long/sample_rate:.1f} ms)")
    print()

    null_results: dict[str, np.ndarray] = {}
    p_values: dict[str, float] = {}

    configs = [
        ("amplitude_shuffle", 0),
        ("block_short", block_short),
        ("block_long", block_long),
    ]

    for null_name, bsize in configs:
        print(f"[null] {null_name}  nperm={args.nperm}")
        vals = run_null(
            null_type="amplitude" if null_name == "amplitude_shuffle" else "block",
            y=y, sample_rate=sample_rate,
            block_size=bsize,
            fmin=args.fmin, fmax=args.fmax,
            k_grid=k_grid,
            nperm=args.nperm,
            seed=args.seed,
            use_gpu=use_gpu,
        )
        null_results[null_name] = vals
        p = p_value_geq(real_stat, vals)
        p_values[null_name] = p
        print(f"  p_{null_name} = {p:.6f}  "
              f"(real={real_stat:.4f}, null_max={vals.max():.4f}, "
              f"null_mean={vals.mean():.4f})")
        pd.DataFrame({f"null_delta_chi2": vals}).to_csv(
            out_dir / f"null_{null_name}.csv", index=False)
        print()

    plot_nulls(real_stat, null_results, out_dir)

    summary = {
        "input": args.input,
        "real_scan_max_delta_chi2": real_stat,
        "fmin": args.fmin,
        "fmax": args.fmax,
        "kmin": args.kmin,
        "kmax": args.kmax,
        "nk": args.nk,
        "nperm_per_method": args.nperm,
        "mains_hz": args.mains_hz,
        "block_short_samples": block_short,
        "block_long_samples": block_long,
        "p_values": p_values,
        "interpretation": {
            "amplitude_shuffle": "Permutes detrended log-spectrum amplitudes (standard null)",
            "block_short": f"Shuffles {block_short}-sample waveform blocks "
                           f"(~half mains cycle); preserves intra-cycle correlations",
            "block_long": f"Shuffles {block_long}-sample waveform blocks "
                          f"(~full mains cycle); preserves full-cycle correlations",
        },
    }

    with open(out_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 70)
    print("BLOCK BOOTSTRAP SUMMARY")
    print("=" * 70)
    print(f"Real scan-max Δχ²   = {real_stat:.6f}")
    print()
    for null_name, p in p_values.items():
        print(f"  p_{null_name:<25} = {p:.6f}")
    print()
    print("Lower p under stricter nulls → result is not explained by")
    print("short-range temporal correlations in the waveform.")
    print(f"\nSaved: {out_dir.resolve()}")


if __name__ == "__main__":
    main()

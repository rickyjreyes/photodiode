#!/usr/bin/env python3
"""
block_bootstrap.py — corrected log-cos null/stress tests.

This version fixes the earlier interpretation problem:
- amplitude_shuffle is the PRIMARY scan-max null for the log-spectrum diagnostic.
- block/cycle shuffles are STRESS TESTS only.
- time-domain block permutation can create boundary discontinuities and spectral leakage,
  so high block-null scores are marked as inconclusive, not as evidence that the
  measured harmonic ratios disappear.

Usage:
  python block_bootstrap.py --input captures/20260309_010802/waveform_raw.npy --sample-rate 20000 --nperm 1000
"""
from __future__ import annotations

import argparse, json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from scan_npy_logcos import (
    load_npy_waveform, compute_fft_spectrum, select_freq_window,
    build_log_domain, scan_logcos_cpu, p_value_geq, CUPY_AVAILABLE,
)
try:
    from scan_npy_logcos import scan_logcos_gpu
    HAS_GPU_SCAN = True
except Exception:
    HAS_GPU_SCAN = False


def logcos_scan_max(y, sample_rate, fmin, fmax, k_grid, use_gpu):
    spec = compute_fft_spectrum(y, sample_rate, use_cupy=False)
    spec_win = select_freq_window(spec, fmin, fmax)
    logdata = build_log_domain(spec_win)
    ell = logdata["ell_ln_f"].to_numpy(float)
    yy = logdata["y_detrended"].to_numpy(float)
    if use_gpu and HAS_GPU_SCAN:
        _, best = scan_logcos_gpu(ell, yy, k_grid)
    else:
        _, best = scan_logcos_cpu(ell, yy, k_grid)
    return float(best["delta_chi2"])


def amplitude_shuffle_stat(y, sample_rate, fmin, fmax, k_grid, use_gpu, rng):
    """Primary null: shuffle detrended log-spectrum amplitudes."""
    spec = compute_fft_spectrum(y, sample_rate, use_cupy=False)
    spec_win = select_freq_window(spec, fmin, fmax)
    logdata = build_log_domain(spec_win)
    ell = logdata["ell_ln_f"].to_numpy(float)
    yy = logdata["y_detrended"].to_numpy(float)
    yy = rng.permutation(yy)
    if use_gpu and HAS_GPU_SCAN:
        _, best = scan_logcos_gpu(ell, yy, k_grid)
    else:
        _, best = scan_logcos_cpu(ell, yy, k_grid)
    return float(best["delta_chi2"])


def estimate_reference_phase(y, sample_rate, mains_hz):
    """Fit y ≈ a cos(wt) + b sin(wt). Return delta in cos(wt-delta)."""
    t = np.arange(len(y), dtype=float) / sample_rate
    w = 2.0 * np.pi * mains_hz
    X = np.column_stack([np.cos(w*t), np.sin(w*t)])
    a, b = np.linalg.lstsq(X, y, rcond=None)[0]
    return float(np.arctan2(b, a))


def cycle_sync_boundaries(n, sample_rate, mains_hz, phase_delta, half_cycle=False):
    """Boundaries at a fitted constant mains phase; avoids fixed 333-sample drift."""
    step = np.pi if half_cycle else 2.0*np.pi
    w = 2.0*np.pi*mains_hz
    duration = n / sample_rate
    m_min = math.floor((-phase_delta) / step) - 2
    m_max = math.ceil((w*duration - phase_delta) / step) + 2
    idxs = [0]
    for m in range(m_min, m_max+1):
        t = (phase_delta + m*step) / w
        idx = int(round(t * sample_rate))
        if 0 < idx < n:
            idxs.append(idx)
    idxs.append(n)
    idxs = sorted(set(idxs))
    clean = [idxs[0]]
    for x in idxs[1:]:
        if x - clean[-1] >= 8:
            clean.append(x)
    if clean[-1] != n:
        clean.append(n)
    return np.asarray(clean, dtype=int)


def cycle_sync_shuffle(y, sample_rate, mains_hz, rng, half_cycle=False):
    phase_delta = estimate_reference_phase(y, sample_rate, mains_hz)
    b = cycle_sync_boundaries(len(y), sample_rate, mains_hz, phase_delta, half_cycle)
    blocks = [y[b[i]:b[i+1]].copy() for i in range(len(b)-1)]
    lengths = np.array([len(x) for x in blocks], dtype=int)
    if len(blocks) < 3:
        return y.copy(), {"n_blocks": len(blocks)}
    order = rng.permutation(len(blocks))
    ys = np.concatenate([blocks[i] for i in order])
    meta = {
        "n_blocks": int(len(blocks)),
        "mean_block_len": float(np.mean(lengths)),
        "min_block_len": int(np.min(lengths)),
        "max_block_len": int(np.max(lengths)),
        "half_cycle": bool(half_cycle),
        "phase_delta": float(phase_delta),
    }
    return ys, meta


def hf_leakage_ratio(y0, y1, sample_rate, min_hz=1000.0):
    """High-frequency energy ratio: shuffled/original. Large values imply boundary leakage."""
    def energy(y):
        spec = compute_fft_spectrum(y, sample_rate, use_cupy=False)
        f = spec["frequency_Hz"].to_numpy(float)
        a = spec["amplitude"].to_numpy(float)
        m = f >= min_hz
        return float(np.sum(a[m]**2))
    return energy(y1) / max(energy(y0), 1e-12)


def run_amplitude_null(y, sr, fmin, fmax, k_grid, nperm, seed, use_gpu):
    rng = np.random.default_rng(seed)
    vals = np.empty(nperm, dtype=float)
    for i in range(nperm):
        vals[i] = amplitude_shuffle_stat(y, sr, fmin, fmax, k_grid, use_gpu, rng)
        if (i+1) % max(1, nperm//10) == 0:
            print(f"  [amplitude] {i+1}/{nperm} running max={vals[:i+1].max():.4f}")
    return vals


def run_cycle_stress(y, sr, fmin, fmax, k_grid, nperm, seed, use_gpu, mains_hz, half_cycle):
    rng = np.random.default_rng(seed)
    vals = np.empty(nperm, dtype=float)
    leak = np.empty(nperm, dtype=float)
    meta0 = None
    for i in range(nperm):
        ys, meta = cycle_sync_shuffle(y, sr, mains_hz, rng, half_cycle=half_cycle)
        if meta0 is None:
            meta0 = meta
        vals[i] = logcos_scan_max(ys, sr, fmin, fmax, k_grid, use_gpu)
        leak[i] = hf_leakage_ratio(y, ys, sr)
        if (i+1) % max(1, nperm//10) == 0:
            print(f"  [cycle {'half' if half_cycle else 'full'}] {i+1}/{nperm} "
                  f"running max={vals[:i+1].max():.4f} median leakage={np.median(leak[:i+1]):.3g}")
    diag = {
        "role": "stress_test_not_primary_pvalue",
        "block_meta": meta0,
        "median_high_frequency_leakage_ratio": float(np.median(leak)),
        "mean_high_frequency_leakage_ratio": float(np.mean(leak)),
        "note": "Cycle/block shuffles may alter harmonic coherence or inject leakage; high p is inconclusive, not erasure of measured ratios.",
    }
    return vals, diag


def plot_nulls(real_stat, results, out_dir):
    fig, ax = plt.subplots(figsize=(10,5))
    for name, vals in results.items():
        ax.hist(vals, bins=40, alpha=0.55, density=True, label=name)
    ax.axvline(real_stat, color="red", lw=2, ls="--", label=f"real Δχ²={real_stat:.3f}")
    ax.set_xlabel("scan-max Δχ²")
    ax.set_ylabel("density")
    ax.set_title("Log-cos nulls / stress tests")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "null_distributions.png", dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--sample-rate", type=float, default=20000.0)
    ap.add_argument("--fmin", type=float, default=20.0)
    ap.add_argument("--fmax", type=float, default=1000.0)
    ap.add_argument("--kmin", type=float, default=0.5)
    ap.add_argument("--kmax", type=float, default=80.0)
    ap.add_argument("--nk", type=int, default=4000)
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--mains-hz", type=float, default=60.0)
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out-dir", default="block_bootstrap_results")
    args = ap.parse_args()

    use_gpu = (not args.cpu) and CUPY_AVAILABLE
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    print(f"[load] {args.input}")
    _, y, sr = load_npy_waveform(Path(args.input), args.sample_rate)
    print(f"[data] samples={len(y):,} sr={sr:.0f} Hz duration={len(y)/sr:.3f}s")

    k_grid = np.linspace(args.kmin, args.kmax, args.nk)
    real_stat = logcos_scan_max(y, sr, args.fmin, args.fmax, k_grid, use_gpu)
    print(f"[scan] real Δχ² = {real_stat:.6f}\n")

    results = {}
    diagnostics = {}

    print("[null] amplitude_shuffle PRIMARY")
    results["amplitude_shuffle"] = run_amplitude_null(
        y, sr, args.fmin, args.fmax, k_grid, args.nperm, args.seed, use_gpu)
    diagnostics["amplitude_shuffle"] = {"role": "primary_scanmax_null"}

    print("\n[stress] cycle_sync_block full cycle")
    vals, diag = run_cycle_stress(
        y, sr, args.fmin, args.fmax, k_grid, args.nperm, args.seed+1, use_gpu, args.mains_hz, False)
    results["cycle_sync_block"] = vals
    diagnostics["cycle_sync_block"] = diag

    print("\n[stress] half_cycle_sync_block")
    vals, diag = run_cycle_stress(
        y, sr, args.fmin, args.fmax, k_grid, args.nperm, args.seed+2, use_gpu, args.mains_hz, True)
    results["half_cycle_sync_block"] = vals
    diagnostics["half_cycle_sync_block"] = diag

    p_values = {}
    for name, vals in results.items():
        p = p_value_geq(real_stat, vals)
        p_values[name] = float(p)
        pd.DataFrame({"null_delta_chi2": vals}).to_csv(out_dir / f"null_{name}.csv", index=False)
        print(f"[result] {name:<22} p={p:.6f} null_max={vals.max():.4f} null_mean={vals.mean():.4f}")

    plot_nulls(real_stat, results, out_dir)

    summary = {
        "input": args.input,
        "real_scan_max_delta_chi2": float(real_stat),
        "nperm_per_method": args.nperm,
        "mains_hz": args.mains_hz,
        "p_values": p_values,
        "diagnostics": diagnostics,
        "interpretation": {
            "primary_claim": "Use amplitude_shuffle as the primary log-spectrum scan-max null.",
            "stress_tests": "Cycle-sync block tests assess dependence on periodic block structure; they are not automatic rejection tests for measured harmonic ratios.",
            "important": "Exact 2/3 pairs and Koide-style triplets are direct frequency-ratio outputs and are not invalidated by block-shuffle stress p-values.",
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n" + "="*70)
    print("UPDATED NULL/STRESS SUMMARY")
    print("="*70)
    print(f"Real scan-max Δχ² = {real_stat:.6f}\n")
    for name, p in p_values.items():
        print(f"  p_{name:<22} = {p:.6f}")
    print("\nSaved:", out_dir.resolve())


if __name__ == "__main__":
    main()

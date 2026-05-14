#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scan_npy_logcos.py

Scan a raw oscilloscope waveform saved as .npy for:
  1. FFT spectral peaks
  2. 2/3, 3/2, octave, and harmonic-ratio pairs
  3. log-cos structure in frequency space:
        y(ell) = C + a cos(k ell) + b sin(k ell), ell = ln(f)
  4. optional shuffle/null test
  5. optional Koide-style triplet search among spectral peaks

Input
-----
A .npy file containing either:
  - 1D waveform voltage samples
  - 2D array with columns [time, voltage]
  - dict-like saved object with keys such as "volts", "signal", "time", "sample_rate"

For 1D waveform arrays, you MUST provide --sample-rate.

Examples
--------
python scan_npy_logcos.py --input waveform_raw.npy --sample-rate 12000 --fmax 1000

python scan_npy_logcos.py --input waveform_raw.npy --sample-rate 12000 --fmin 20 --fmax 500 --nperm 5000

Outputs
-------
outputs_npy_logcos/
  spectrum.csv
  top_peaks.csv
  ratio_pairs.csv
  koide_triplets.csv
  logcos_scan.csv
  logcos_null.csv
  summary.json
  waveform.png
  spectrum.png
  logcos_scan.png
  logcos_best_fit.png
"""

from __future__ import annotations

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.signal import find_peaks


# ============================================================
# Optional CuPy
# ============================================================

try:
    import cupy as cp
    CUPY_AVAILABLE = True
except Exception:
    cp = None
    CUPY_AVAILABLE = False


# ============================================================
# Defaults
# ============================================================

OUT_DIR = Path("outputs_npy_logcos")
OUT_DIR.mkdir(exist_ok=True, parents=True)

KOIDE_Q = 2.0 / 3.0


# ============================================================
# Loading
# ============================================================

def load_npy_waveform(path: Path, sample_rate_arg: float | None):
    obj = np.load(path, allow_pickle=True)

    sample_rate = sample_rate_arg
    t = None
    y = None

    # Handle 0-d object arrays containing dicts.
    if isinstance(obj, np.ndarray) and obj.shape == () and obj.dtype == object:
        obj = obj.item()

    if isinstance(obj, dict):
        # Common possible keys.
        time_keys = ["time", "t", "time_axis", "timestamp", "timestamps"]
        y_keys = ["volts", "voltage", "waveform", "signal", "y", "data", "samples"]
        sr_keys = ["sample_rate", "fs", "sampling_rate", "sampleRate"]

        for k in sr_keys:
            if k in obj and sample_rate is None:
                sample_rate = float(obj[k])

        for k in time_keys:
            if k in obj:
                t = np.asarray(obj[k], dtype=float)
                break

        for k in y_keys:
            if k in obj:
                y = np.asarray(obj[k], dtype=float)
                break

        if y is None:
            raise ValueError(f"Dict .npy did not contain waveform keys. Keys: {list(obj.keys())}")

    else:
        arr = np.asarray(obj)

        if arr.ndim == 1:
            y = arr.astype(float)

        elif arr.ndim == 2:
            # Assume columns [time, signal] if at least 2 columns.
            if arr.shape[1] >= 2:
                t = arr[:, 0].astype(float)
                y = arr[:, 1].astype(float)
            elif arr.shape[0] >= 2:
                t = arr[0, :].astype(float)
                y = arr[1, :].astype(float)
            else:
                raise ValueError(f"Unsupported 2D shape: {arr.shape}")

        else:
            raise ValueError(f"Unsupported .npy array shape: {arr.shape}")

    y = np.asarray(y, dtype=float).ravel()
    good = np.isfinite(y)
    y = y[good]

    if t is not None:
        t = np.asarray(t, dtype=float).ravel()
        t = t[good]
        order = np.argsort(t)
        t = t[order]
        y = y[order]

        dt = float(np.nanmedian(np.diff(t)))
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Could not infer positive dt from time column.")

        inferred_sr = 1.0 / dt
        if sample_rate is None:
            sample_rate = inferred_sr
        else:
            # Prefer explicit sample-rate if provided.
            pass

    if sample_rate is None:
        raise ValueError(
            "1D waveform has no time axis. Provide --sample-rate from the oscilloscope, "
            "for example: --sample-rate 12000"
        )

    if t is None:
        dt = 1.0 / float(sample_rate)
        t = np.arange(len(y), dtype=float) * dt

    # Remove DC.
    y = y - np.mean(y)

    return t, y, float(sample_rate)


# ============================================================
# FFT / spectrum
# ============================================================

def compute_fft_spectrum(y: np.ndarray, sample_rate: float, use_cupy: bool):
    dt = 1.0 / sample_rate

    # Window to reduce leakage.
    win = np.hanning(len(y))
    yw = y * win

    if use_cupy and CUPY_AVAILABLE:
        yg = cp.asarray(yw)
        fftg = cp.fft.rfft(yg)
        spec = cp.asnumpy(cp.abs(fftg))
        freq = cp.asnumpy(cp.fft.rfftfreq(len(yw), d=dt))
    else:
        spec = np.abs(np.fft.rfft(yw))
        freq = np.fft.rfftfreq(len(yw), d=dt)

    return pd.DataFrame({
        "frequency_Hz": freq.astype(float),
        "amplitude": spec.astype(float),
    })


def select_freq_window(spec: pd.DataFrame, fmin: float, fmax: float):
    d = spec.copy()
    d = d[np.isfinite(d["frequency_Hz"]) & np.isfinite(d["amplitude"])]
    d = d[(d["frequency_Hz"] >= fmin) & (d["frequency_Hz"] <= fmax)]
    d = d[d["frequency_Hz"] > 0]
    d = d.sort_values("frequency_Hz").reset_index(drop=True)

    if len(d) < 8:
        raise RuntimeError(f"Too few spectrum points in frequency window: {len(d)}")

    return d


def detect_peaks(spec: pd.DataFrame, peak_height_frac: float, top_n: int):
    freq = spec["frequency_Hz"].to_numpy(float)
    amp = spec["amplitude"].to_numpy(float)

    height = np.max(amp) * peak_height_frac
    peaks, props = find_peaks(amp, height=height)

    df = pd.DataFrame({
        "frequency_Hz": freq[peaks],
        "amplitude": props["peak_heights"],
    })

    df = df.sort_values("amplitude", ascending=False).head(top_n).reset_index(drop=True)
    df["rank_by_amplitude"] = np.arange(1, len(df) + 1)

    # Also add frequency-sorted rank.
    df = df.sort_values("frequency_Hz").reset_index(drop=True)
    df["rank_by_frequency"] = np.arange(1, len(df) + 1)

    if len(df):
        f0 = float(df.iloc[0]["frequency_Hz"])
        a0 = float(df.iloc[0]["amplitude"])
        df["ratio_to_lowest_peak"] = df["frequency_Hz"] / f0
        df["amplitude_ratio_to_lowest_peak"] = df["amplitude"] / a0

    return df


# ============================================================
# Ratio / Koide reports
# ============================================================

def make_ratio_pairs(peaks: pd.DataFrame):
    rows = []
    fs = peaks["frequency_Hz"].to_numpy(float)
    amps = peaks["amplitude"].to_numpy(float)

    for i, j in combinations(range(len(fs)), 2):
        f_lo, f_hi = sorted([fs[i], fs[j]])
        a_lo = amps[i] if fs[i] == f_lo else amps[j]
        a_hi = amps[j] if fs[j] == f_hi else amps[i]

        r = f_lo / f_hi

        rows.append({
            "f_low": float(f_lo),
            "f_high": float(f_hi),
            "ratio_low_over_high": float(r),
            "amp_low": float(a_lo),
            "amp_high": float(a_hi),
            "amp_ratio_low_over_high": float(a_lo / a_hi) if a_hi != 0 else np.nan,
            "delta_from_2_over_3": float(abs(r - 2.0 / 3.0)),
            "delta_from_1_over_2": float(abs(r - 0.5)),
            "delta_from_3_over_4": float(abs(r - 0.75)),
            "delta_from_3_over_2_inverse": float(abs(r - 1.0 / 1.5)),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    return out.sort_values("delta_from_2_over_3").reset_index(drop=True)


def koide_triplets(peaks: pd.DataFrame):
    """
    For ordered frequency triplets f1 < f2 < f3, compute:
        Q_low  = f1 / f2
        Q_high = f3 / (2 f2)

    Koide sideband match:
        Q_low ~= Q_high ~= 2/3

    This is the same geometry form used in the LHC winding analysis:
        (n-, n0, n+) = n0 * (Q, 1, 2Q)
    """
    rows = []
    d = peaks.sort_values("frequency_Hz").reset_index(drop=True)
    fs = d["frequency_Hz"].to_numpy(float)
    amps = d["amplitude"].to_numpy(float)

    for i, j, k in combinations(range(len(fs)), 3):
        f1, f2, f3 = fs[i], fs[j], fs[k]
        q_low = f1 / f2
        q_high = f3 / (2.0 * f2)
        q_mean = 0.5 * (q_low + q_high)
        koide_error = math.sqrt((q_low - KOIDE_Q)**2 + (q_high - KOIDE_Q)**2)

        rows.append({
            "f1": float(f1),
            "f2": float(f2),
            "f3": float(f3),
            "A1": float(amps[i]),
            "A2": float(amps[j]),
            "A3": float(amps[k]),
            "Q_low_f1_over_f2": float(q_low),
            "Q_high_f3_over_2f2": float(q_high),
            "Q_mean": float(q_mean),
            "koide_error": float(koide_error),
            "integer_120_240_360_error": float(
                math.sqrt((f1 - 120.0)**2 + (f2 - 240.0)**2 + (f3 - 360.0)**2)
            ),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    return out.sort_values(["koide_error", "integer_120_240_360_error"]).reset_index(drop=True)


# ============================================================
# Log-cos scan
# ============================================================

def build_log_domain(spec: pd.DataFrame):
    f = spec["frequency_Hz"].to_numpy(float)
    amp = spec["amplitude"].to_numpy(float)

    ell = np.log(f)

    # Work on log amplitude to avoid one huge peak dominating.
    y_raw = np.log1p(amp)

    # Remove smooth trend in ell using quadratic polynomial.
    X = np.column_stack([np.ones_like(ell), ell, ell**2])
    beta, *_ = np.linalg.lstsq(X, y_raw, rcond=None)
    trend = X @ beta
    y = y_raw - trend

    # Standardize.
    sd = np.std(y)
    if sd > 0:
        y = y / sd

    out = spec.copy()
    out["ell_ln_f"] = ell
    out["log1p_amplitude"] = y_raw
    out["trend"] = trend
    out["y_detrended"] = y
    return out


def fit_logcos_cpu(ell, y, k):
    X0 = np.ones((len(ell), 1))
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    r0 = y - X0 @ b0
    chi0 = float(np.sum(r0 * r0))

    X1 = np.column_stack([np.ones_like(ell), np.cos(k * ell), np.sin(k * ell)])
    b1, *_ = np.linalg.lstsq(X1, y, rcond=None)
    r1 = y - X1 @ b1
    chi1 = float(np.sum(r1 * r1))

    C, a, b = map(float, b1)
    A = float(math.sqrt(a*a + b*b))
    phi = float(math.atan2(-b, a))

    return {
        "k": float(k),
        "delta_chi2": float(chi0 - chi1),
        "chi2_null": float(chi0),
        "chi2_logcos": float(chi1),
        "C": C,
        "a": a,
        "b": b,
        "A": A,
        "phi": phi,
    }


def scan_logcos_cpu(ell, y, k_grid):
    rows = [fit_logcos_cpu(ell, y, k) for k in k_grid]
    scan = pd.DataFrame(rows)
    best = scan.sort_values("delta_chi2", ascending=False).iloc[0].to_dict()
    return scan, best


def scan_logcos_gpu(ell, y, k_grid):
    """
    Fast GPU least-squares scan over all k.
    """
    if not CUPY_AVAILABLE:
        return scan_logcos_cpu(ell, y, k_grid)

    ell_g = cp.asarray(ell, dtype=cp.float64)
    y_g = cp.asarray(y, dtype=cp.float64)
    k_g = cp.asarray(k_grid, dtype=cp.float64)

    n = len(ell)
    one = cp.ones((len(k_grid), n), dtype=cp.float64)
    phase = k_g[:, None] * ell_g[None, :]
    c = cp.cos(phase)
    s = cp.sin(phase)

    # Build normal matrices for each k.
    A00 = cp.sum(one * one, axis=1)
    A01 = cp.sum(one * c, axis=1)
    A02 = cp.sum(one * s, axis=1)
    A11 = cp.sum(c * c, axis=1)
    A12 = cp.sum(c * s, axis=1)
    A22 = cp.sum(s * s, axis=1)

    M = cp.stack([
        cp.stack([A00, A01, A02], axis=1),
        cp.stack([A01, A11, A12], axis=1),
        cp.stack([A02, A12, A22], axis=1),
    ], axis=1)

    rhs0 = cp.sum(y_g)
    rhs1 = cp.sum(c * y_g[None, :], axis=1)
    rhs2 = cp.sum(s * y_g[None, :], axis=1)
    rhs = cp.stack([cp.full_like(rhs1, rhs0), rhs1, rhs2], axis=1)

    beta = cp.linalg.solve(M, rhs[:, :, None])[:, :, 0]

    C = beta[:, 0]
    a = beta[:, 1]
    b = beta[:, 2]

    yhat = C[:, None] + a[:, None] * c + b[:, None] * s
    resid = y_g[None, :] - yhat
    chi1 = cp.sum(resid * resid, axis=1)

    # Null constant-only chi2.
    ymean = cp.mean(y_g)
    chi0 = cp.sum((y_g - ymean)**2)
    delta = chi0 - chi1

    Aamp = cp.sqrt(a*a + b*b)
    phi = cp.arctan2(-b, a)

    scan = pd.DataFrame({
        "k": cp.asnumpy(k_g),
        "delta_chi2": cp.asnumpy(delta),
        "chi2_null": float(cp.asnumpy(chi0)),
        "chi2_logcos": cp.asnumpy(chi1),
        "C": cp.asnumpy(C),
        "a": cp.asnumpy(a),
        "b": cp.asnumpy(b),
        "A": cp.asnumpy(Aamp),
        "phi": cp.asnumpy(phi),
    })

    best = scan.sort_values("delta_chi2", ascending=False).iloc[0].to_dict()
    return scan, best


def active_winding(k, fmin, fmax):
    delta_ell = math.log(fmax / fmin)
    return float(k * delta_ell / (2.0 * math.pi))


def shuffle_null(ell, y, k_grid, nperm, seed, use_gpu):
    rng = np.random.default_rng(seed)
    null_best = np.empty(nperm, dtype=float)
    null_k = np.empty(nperm, dtype=float)

    for i in range(nperm):
        yp = rng.permutation(y)
        if use_gpu:
            _, best = scan_logcos_gpu(ell, yp, k_grid)
        else:
            _, best = scan_logcos_cpu(ell, yp, k_grid)
        null_best[i] = best["delta_chi2"]
        null_k[i] = best["k"]

        if (i + 1) % max(1, nperm // 10) == 0:
            print(f"[null] {i+1}/{nperm}")

    return pd.DataFrame({
        "null_best_delta_chi2": null_best,
        "null_best_k": null_k,
    })


def p_value_geq(real, null_vals):
    return float((1 + np.sum(null_vals >= real)) / (1 + len(null_vals)))


# ============================================================
# Plotting
# ============================================================

def plot_waveform(t, y, out):
    plt.figure(figsize=(10, 4))
    plt.plot(t, y, lw=0.8)
    plt.xlabel("Time [s]")
    plt.ylabel("Voltage / signal [DC removed]")
    plt.title("Raw waveform")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()


def plot_spectrum(spec, peaks, out, fmax_plot):
    plt.figure(figsize=(10, 4))
    plt.plot(spec["frequency_Hz"], spec["amplitude"], lw=1)
    if len(peaks):
        plt.scatter(peaks["frequency_Hz"], peaks["amplitude"], s=30)
    plt.xlim(0, fmax_plot)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("FFT amplitude")
    plt.title("FFT spectrum")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()


def plot_logcos_scan(scan, out):
    plt.figure(figsize=(10, 4))
    plt.plot(scan["k"], scan["delta_chi2"], lw=1)
    plt.xlabel("log-frequency k")
    plt.ylabel("Delta chi2")
    plt.title("Log-cos scan")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()


def plot_best_fit(logdata, best, out):
    ell = logdata["ell_ln_f"].to_numpy(float)
    y = logdata["y_detrended"].to_numpy(float)

    k = float(best["k"])
    C = float(best["C"])
    a = float(best["a"])
    b = float(best["b"])
    yhat = C + a * np.cos(k * ell) + b * np.sin(k * ell)

    order = np.argsort(ell)

    plt.figure(figsize=(10, 4))
    plt.scatter(ell, y, s=8, label="detrended log spectrum")
    plt.plot(ell[order], yhat[order], lw=2, label=f"best k={k:.4f}")
    plt.xlabel("ell = ln(f)")
    plt.ylabel("standardized residual log amplitude")
    plt.title("Best log-cos fit")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help=".npy waveform file")
    ap.add_argument("--sample-rate", type=float, default=None, help="oscilloscope sample rate in Hz")
    ap.add_argument("--fmin", type=float, default=20.0)
    ap.add_argument("--fmax", type=float, default=1000.0)
    ap.add_argument("--plot-fmax", type=float, default=500.0)
    ap.add_argument("--peak-height-frac", type=float, default=0.05)
    ap.add_argument("--top-peaks", type=int, default=30)
    ap.add_argument("--kmin", type=float, default=0.5)
    ap.add_argument("--kmax", type=float, default=80.0)
    ap.add_argument("--nk", type=int, default=4000)
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--no-null", action="store_true")
    ap.add_argument("--cpu", action="store_true", help="force CPU even if CuPy is installed")
    args = ap.parse_args()

    use_gpu = (not args.cpu) and CUPY_AVAILABLE

    print("[load]", args.input)
    t, y, sample_rate = load_npy_waveform(Path(args.input), args.sample_rate)

    print(f"[data] samples={len(y):,}")
    print(f"[data] sample_rate={sample_rate:.6g} Hz")
    print(f"[data] duration={len(y)/sample_rate:.6g} s")
    print(f"[gpu] CuPy available={CUPY_AVAILABLE}, using_gpu={use_gpu}")

    spec = compute_fft_spectrum(y, sample_rate, use_cupy=use_gpu)
    spec_win = select_freq_window(spec, args.fmin, args.fmax)

    peaks = detect_peaks(spec_win, args.peak_height_frac, args.top_peaks)
    ratios = make_ratio_pairs(peaks)
    triplets = koide_triplets(peaks)

    logdata = build_log_domain(spec_win)
    ell = logdata["ell_ln_f"].to_numpy(float)
    yy = logdata["y_detrended"].to_numpy(float)

    k_grid = np.linspace(args.kmin, args.kmax, args.nk)

    print("[scan] log-cos")
    if use_gpu:
        scan, best = scan_logcos_gpu(ell, yy, k_grid)
    else:
        scan, best = scan_logcos_cpu(ell, yy, k_grid)

    best["active_delta_ell"] = float(math.log(args.fmax / args.fmin))
    best["active_winding_n"] = active_winding(float(best["k"]), args.fmin, args.fmax)

    null_df = pd.DataFrame()
    if not args.no_null and args.nperm > 0:
        print("[null] shuffle scanmax")
        null_df = shuffle_null(ell, yy, k_grid, args.nperm, args.seed, use_gpu=use_gpu)
        p = p_value_geq(float(best["delta_chi2"]), null_df["null_best_delta_chi2"].to_numpy(float))
        best["p_shuffle_scanmax"] = float(p)
        best["null_n"] = int(args.nperm)
    else:
        best["p_shuffle_scanmax"] = None
        best["null_n"] = 0

    OUT_DIR.mkdir(exist_ok=True, parents=True)
    spec_win.to_csv(OUT_DIR / "spectrum.csv", index=False)
    logdata.to_csv(OUT_DIR / "spectrum_log_domain.csv", index=False)
    peaks.to_csv(OUT_DIR / "top_peaks.csv", index=False)
    ratios.to_csv(OUT_DIR / "ratio_pairs.csv", index=False)
    triplets.to_csv(OUT_DIR / "koide_triplets.csv", index=False)
    scan.to_csv(OUT_DIR / "logcos_scan.csv", index=False)
    if len(null_df):
        null_df.to_csv(OUT_DIR / "logcos_null.csv", index=False)

    plot_waveform(t, y, OUT_DIR / "waveform.png")
    plot_spectrum(spec_win, peaks, OUT_DIR / "spectrum.png", args.plot_fmax)
    plot_logcos_scan(scan, OUT_DIR / "logcos_scan.png")
    plot_best_fit(logdata, best, OUT_DIR / "logcos_best_fit.png")

    summary = {
        "input": args.input,
        "samples": int(len(y)),
        "sample_rate_Hz": float(sample_rate),
        "duration_s": float(len(y) / sample_rate),
        "frequency_window_Hz": [float(args.fmin), float(args.fmax)],
        "k_scan": [float(args.kmin), float(args.kmax), int(args.nk)],
        "cupy_available": bool(CUPY_AVAILABLE),
        "using_gpu": bool(use_gpu),
        "best_logcos": best,
        "top_peaks": peaks.head(12).to_dict(orient="records"),
        "closest_2_over_3_pairs": ratios.head(12).to_dict(orient="records") if len(ratios) else [],
        "best_koide_triplets": triplets.head(12).to_dict(orient="records") if len(triplets) else [],
    }

    with open(OUT_DIR / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 80)
    print("NPY LOG-COS SCAN COMPLETE")
    print("=" * 80)

    print("\nBest log-cos:")
    print(f"  k                 = {best['k']:.6f}")
    print(f"  Delta chi2        = {best['delta_chi2']:.6f}")
    print(f"  A                 = {best['A']:.6f}")
    print(f"  phi               = {best['phi']:.6f}")
    print(f"  active Delta ell  = {best['active_delta_ell']:.6f}")
    print(f"  active winding n  = {best['active_winding_n']:.6f}")
    if best["p_shuffle_scanmax"] is not None:
        print(f"  p_shuffle_scanmax = {best['p_shuffle_scanmax']:.6g}")

    print("\nTop peaks:")
    if len(peaks):
        print(peaks.head(12).to_string(index=False))
    else:
        print("  none found; lower --peak-height-frac")

    print("\nClosest 2/3 pairs:")
    if len(ratios):
        print(ratios[["f_low", "f_high", "ratio_low_over_high", "delta_from_2_over_3"]].head(12).to_string(index=False))
    else:
        print("  not enough peaks")

    print("\nBest Koide-style triplets:")
    if len(triplets):
        print(triplets[["f1", "f2", "f3", "Q_low_f1_over_f2", "Q_high_f3_over_2f2", "Q_mean", "koide_error"]].head(12).to_string(index=False))
    else:
        print("  not enough peaks")

    print(f"\nSaved outputs to: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()

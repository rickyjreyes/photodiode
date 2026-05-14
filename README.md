# Optically Induced Harmonic State in a Silicon Photodiode

**Reproducible artifact package — prediction ledger, oscilloscope captures, FFT analysis, Q = 2/3 sideband geometry, and log-cos spectral scan with shuffle null.**

---

## Overview

This repository is a self-contained reproducibility package for a captured electrical waveform from a silicon photodiode following a brief optical excitation protocol.

The central measured fact is this: the post-excitation waveform contains a discrete 60 Hz-spaced harmonic ladder with exact Q = 2/3 ratio geometry across two independent triplets, and the log-frequency spectrum carries statistically significant cosine modulation (Δχ² = 21.139, p\_shuffle\_scanmax = 0.000999 over 1000 null trials).

All scripts are executable. All outputs are committed. The prediction ledger records what was expected before the capture.

---

## Repository Contents

```
README.md                    — this file
proof_ledger.json            — pre-declared protocol and prediction record
scan_npy_logcos.py           — main analysis script (FFT → ratio → Koide → log-cos)
siglent.py / siglent_2.py    — oscilloscope capture scripts (Siglent SDS1104X HD)
siglent_proof.py             — capture + analysis pipeline
photodiode_proof.py          — induction protocol driver

captures/
  20260309_010802/
    waveform_raw.npy         — primary artifact waveform

csv/                         — raw spectral run CSVs and combined harmonic summaries
proof_runs/                  — baseline and post-excitation FFT figures and peak CSVs
proof/                       — annotated spectral figures
outputs_npy_logcos/          — canonical analysis outputs (see below)
control/                     — control run data
pic/                         — device and setup photographs
```

### Canonical Analysis Outputs (`outputs_npy_logcos/`)

| File | Contents |
|---|---|
| `summary.json` | All scalar results: peaks, ratio pairs, Koide triplets, log-cos best fit |
| `top_peaks.csv` | Detected FFT peaks with frequency and amplitude |
| `ratio_pairs.csv` | All detected peak pairs with ratio and δ from 2/3 |
| `koide_triplets.csv` | All ordered triplets with Q\_low, Q\_high, Q\_mean, koide\_error |
| `logcos_scan.csv` | Δχ² vs. k across the full scan range |
| `logcos_null.csv` | Shuffle null distribution (1000 trials, scan-max statistic) |
| `spectrum.png` | FFT amplitude spectrum |
| `logcos_scan.png` | Δχ² scan across k |
| `logcos_best_fit.png` | Best-fit log-cos overlay on log-frequency spectrum |
| `waveform.png` | Raw time-domain waveform |

---

## Primary Artifact

**Device:** Siglent Technologies SDS1104X HD  
**Channel:** C3  
**Capture date:** 2026-03-09  
**Sample rate:** 20 000 Hz  
**Samples:** 10 000  
**Duration:** 0.5 s  
**Vertical scale:** 3.0 V/div, offset 3.0 V

The raw waveform is stored as `captures/20260309_010802/waveform_raw.npy`.

---

## Reproduce the Scan

### 1. Install dependencies

```bash
pip install numpy pandas scipy matplotlib python-dotenv
```

GPU acceleration (optional, used for the committed outputs):

```bash
pip install cupy
```

### 2. Run the log-cos scan

```bash
python scan_npy_logcos.py \
  --input captures/20260309_010802/waveform_raw.npy \
  --sample-rate 20000 \
  --fmin 20 \
  --fmax 1000
```

For a stronger null test (5000 shuffles instead of 1000):

```bash
python scan_npy_logcos.py \
  --input captures/20260309_010802/waveform_raw.npy \
  --sample-rate 20000 \
  --fmin 20 \
  --fmax 1000 \
  --nperm 5000
```

Outputs are written to `outputs_npy_logcos/`. Compare `summary.json` against the committed version to verify reproduction.

### 3. Capture a live waveform (optional)

Configure the oscilloscope IP:

```bash
cp .env.example .env
# edit .env: SCOPE_IP=<your scope address>
```

Then run:

```bash
python siglent_proof.py
```

---

## Observed Spectral Structure

The FFT of the captured waveform resolves discrete peaks at integer multiples of 60 Hz:

| Frequency (Hz) | Amplitude (arb.) | Rank by amplitude |
|---|---|---|
| 60 | 60.43 | 3 |
| 120 | 191.00 | 1 |
| 180 | 39.13 | 5 |
| 240 | 114.99 | 2 |
| 300 | 14.05 | 7 |
| 360 | 47.97 | 4 |
| 480 | 17.84 | 6 |
| 600 | 13.94 | 8 |

The ladder is 60 Hz-spaced throughout. Because 60 Hz harmonics can arise from mains or environmental pickup, controls are essential to interpret the origin of this structure (see [Controls and Falsification](#controls-and-falsification)). The ratio and log-cos analyses are reproducible measurements of the artifact regardless of source.

---

## Q = 2/3 Sideband Geometry

### Exact ratio pairs

Two peak pairs satisfy the 2/3 ratio exactly (δ = 0.000):

```
120 Hz / 180 Hz = 2/3
240 Hz / 360 Hz = 2/3
```

### Exact Koide-style triplets

For ordered triplets (f₁, f₂, f₃) the diagnostic computes:

```
Q_low  = f₁ / f₂
Q_high = f₃ / (2 f₂)
Q_mean = (Q_low + Q_high) / 2
```

Two triplets satisfy Q\_low = Q\_high = 2/3 exactly:

| Triplet (Hz) | Q\_low | Q\_high | Q\_mean | koide\_error |
|---|---|---|---|---|
| (120, 180, 240) | 2/3 | 2/3 | 2/3 | 0.0 |
| (240, 360, 480) | 2/3 | 2/3 | 2/3 | 0.0 |

Both triplets reduce to the integer ratio 2:3:4. The same sideband geometry appears at two independent frequency scales in the same capture.

This is called **Q = 2/3 sideband geometry** throughout this repository. It is a statement about the dimensionless ratio structure of the measured frequency spectrum, not a claim about particle masses or the Koide formula for leptons. The shared diagnostic form is noted because it provides a compact, exact three-frequency locking condition.

---

## Log-Cos Scan

### Method

The FFT amplitude spectrum is mapped to log-frequency space:

```
ℓ = ln(f)
```

over the window f ∈ [20, 1000] Hz. The tested model at each candidate wavenumber k is:

```
y(ℓ) = C + a cos(k ℓ) + b sin(k ℓ)
```

The scan searches k ∈ [0.5, 80] at 4000 steps, fitting by least squares at each k and recording:

```
Δχ² = χ²_null − χ²_logcos
```

where χ²\_null is the variance-normalized residual of the constant-only fit.

### Null test

The significance of the scan maximum is assessed against a shuffle null: the FFT amplitudes are permuted 1000 times, the full k-scan is repeated on each shuffled spectrum, and the scan-max Δχ² is recorded. The reported p-value is the fraction of shuffled trials whose scan-max exceeded the real scan-max.

This is a **scan-max shuffle null** — it accounts for the look-elsewhere effect across the full tested k range.

### Observed result

```
Best k             =  44.5341
Δχ²                =  21.139
χ²_null            = 491.000
χ²_logcos          = 469.861
Amplitude A        =   0.293
Null trials        =  1000
p_shuffle_scanmax  =  0.000999   (1/1001 — resolution floor of this run)
```

No shuffled spectrum exceeded the real scan statistic in this run.

The p-value of 0.000999 is the measured scan-max shuffle result for this run and this analysis window. It is not a Gaussian sigma conversion and should not be generalized beyond the tested k range and sample.

---

## Prediction Ledger

`proof_ledger.json` records the pre-declared protocol and expected observables.

The ledger was written before the canonical capture and specifies:

- the excitation protocol sequence,
- the expected transition from baseline to structured output,
- the expected presence of a discrete harmonic ladder,
- the expected persistence of structure after illumination ceases,
- the expected ratio geometry,
- the falsification conditions.

The ledger establishes that the reported structure was predicted, not reconstructed after the fact. The baseline run (epoch 1773046194) captures a virgin device prior to excitation and is included in `proof_runs/` for direct comparison against the post-excitation artifact.

---

## Controls and Falsification

The following controls are necessary to distinguish an induced device state from ordinary environmental pickup or instrument artifact. They are the direct falsification tests for the central claim.

| Control | What it tests |
|---|---|
| Virgin photodiode, excitation omitted | Does the harmonic structure appear without the protocol? |
| Sham-excited control (protocol run, diode dark or absent) | Does the excitation sequence itself drive the observation? |
| Disconnected input / terminated channel | Is the structure present in the measurement chain alone? |
| Dark enclosure, no optical path | Is ambient light driving the observation? |
| Scope channel swap | Is the signal channel-specific or instrument-wide? |
| Battery-powered or preamplifier-isolated supply | Does mains coupling in the supply account for the 60 Hz ladder? |
| Shielding and grounding variation | Does the structure track environmental EM conditions? |
| Repeat captures across devices and days | Is the effect reproducible across hardware instances and time? |
| Pre-declared parameter lock (fixed fmin, fmax, k range) | Are scan parameters chosen post-hoc to fit the result? |

The central claim is **falsified** if:

- virgin devices (excitation omitted) produce identical spectral structure under blinded comparison,
- the structure vanishes immediately after excitation ceases,
- the Q = 2/3 triplets do not recur in independent captures under identical protocol,
- disconnected-channel or terminated-input controls reproduce the full harmonic ladder.

All control runs should be processed through the same `scan_npy_logcos.py` pipeline and committed alongside the primary artifact so outputs are directly comparable.

---

## Interpretation

Three independently computable structures are present in the artifact:

1. **A 60 Hz-spaced discrete harmonic ladder** across at least 8 resolved peaks.
2. **Exact Q = 2/3 sideband geometry** in two non-overlapping triplets at different frequency scales.
3. **Statistically significant log-cos modulation** in ln(f) with p\_shuffle\_scanmax at the 1/1001 floor over 1000 trials.

The ratio geometry is scale-free: both triplets (120, 180, 240) and (240, 360, 480) encode the same 2:3:4 integer structure at different absolute frequencies. The spectral structure is therefore not a single-frequency artifact but a geometric relationship preserved across scales within the same capture.

Whether this structure reflects an optically induced device state or a reproducible environmental/pickup pattern is the open experimental question. The controls listed above are the direct path to resolving it. This repository provides the baseline artifact, the canonical analysis outputs, and the pipeline against which all control runs are to be compared.

---

## Citation / Author

**Repository:** `rickyjreyes/photodiode`  
**Capture instrument:** Siglent Technologies SDS1104X HD (firmware 5.5.6.1.1.0.2)  
**Primary capture:** 2026-03-09, epoch 1773046275

If citing this artifact package, reference the repository commit hash of the canonical run together with `summary.json` for exact scalar reproducibility.

> [Placeholder — author affiliation and formal citation to be added.]

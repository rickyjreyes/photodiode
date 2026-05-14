# Photodiode State Induction

Prediction ledger, oscilloscope captures, and reproducible spectral analysis for an optically induced long-lived harmonic state in a silicon photodiode.

## Overview

This repository contains a reproducible artifact package for testing whether a brief optical excitation protocol can induce a persistent harmonic electrical state in a standard silicon photodiode.

The central claim is not that a photodiode detects light. Photodiodes do that conventionally.

The claim tested here is stronger:

> A short optical excitation protocol may drive a silicon photodiode into a long-lived structured electrical state whose later oscilloscope readout contains discrete harmonic structure, exact ratio locks, and statistically significant log-frequency modulation.

The repository includes:

- A prediction and protocol ledger.
- Siglent oscilloscope capture scripts.
- Raw or curated waveform artifacts.
- FFT spectral outputs.
- Ratio-pair analysis.
- Koide-style `Q = 2/3` triplet analysis.
- Log-cos scans in `ell = ln(f)`.
- Shuffle/null test outputs.

## Key Result

A captured waveform artifact shows a structured harmonic ladder with peaks at:

```text
60 Hz
120 Hz
180 Hz
240 Hz
300 Hz
360 Hz
480 Hz
600 Hz
```

The artifact contains exact `2/3` frequency-ratio locks:

```text
120 / 180 = 2/3
240 / 360 = 2/3
```

It also contains exact Koide-style sideband triplets:

```text
(120, 180, 240)
(240, 360, 480)
```

For these triplets:

```text
Q_low  = f1 / f2      = 2/3
Q_high = f3 / (2*f2) = 2/3
Q_mean = 2/3
```

with Koide-style geometry error:

```text
0.0
```

The log-cos scan also detects statistically significant structure in log-frequency space:

```text
k = 44.5341
Delta chi2 = 21.139
p_shuffle_scanmax = 0.000999
```

The p-value is the floor for 1000 shuffle/null trials, meaning no shuffled spectrum exceeded the real scan statistic in that run.

## Repository Purpose

This repository is intended to make the artifact independently inspectable.

The goal is simple:

```text
data -> scripts -> outputs -> ratios -> null test
```

Anyone should be able to clone the repository, run the scan, and reproduce the reported spectral structure.

This repository does not require accepting any microscopic mechanism in advance. The primary claim is empirical:

> The waveform artifact contains persistent, structured harmonic behavior with exact `Q = 2/3` sideband geometry and significant log-cos modulation.

## Prediction Ledger

The accompanying protocol ledger records a pre-declared claim that a brief optical excitation sequence can induce a long-lived harmonic state in a silicon photodiode.

The predicted observables include:

- transition from baseline to structured output,
- a discrete harmonic ladder,
- stable frequency ratios,
- persistence after illumination ceases,
- invariance under later readout,
- falsifiability by absence of the effect in virgin devices or by immediate decay after excitation.

The ledger is included to establish that the central claim was protocol-defined rather than reconstructed after the spectral scan.

## Reproduction

### 1. Install dependencies

```bash
pip install numpy pandas scipy matplotlib python-dotenv
```

Optional GPU acceleration:

```bash
pip install cupy
```

### 2. Configure oscilloscope IP

Create a private `.env` file:

```env
SCOPE_IP=
```

Do not commit your real scope IP.

The repository includes `.env.example` as a blank template.

### 3. Run the `.npy` artifact scan

```bash
python scripts/scan_npy_logcos.py --input data/waveform_raw.npy --sample-rate 20000 --fmin 20 --fmax 1000
```

For a stronger null test:

```bash
python scripts/scan_npy_logcos.py --input data/waveform_raw.npy --sample-rate 20000 --fmin 20 --fmax 1000 --nperm 5000
```

### 4. Expected outputs

The scan writes:

```text
outputs_npy_logcos/
  spectrum.csv
  spectrum_log_domain.csv
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
```

The most important files are:

```text
summary.json
top_peaks.csv
ratio_pairs.csv
koide_triplets.csv
logcos_scan.csv
logcos_null.csv
```

## Main Analysis Steps

### FFT Spectrum

The raw waveform is transformed into a frequency-domain spectrum. Peaks are detected from the FFT amplitude spectrum.

### Ratio-Pair Scan

Detected peaks are compared pairwise to identify exact or near-exact rational locks, especially:

```text
2/3
```

The artifact contains:

```text
120 / 180 = 2/3
240 / 360 = 2/3
```

### Koide-Style Triplet Scan

For ordered triplets:

```text
f1 < f2 < f3
```

the code computes:

```text
Q_low  = f1 / f2
Q_high = f3 / (2*f2)
```

The exact sideband condition is:

```text
Q_low = Q_high = 2/3
```

The artifact contains exact triplets:

```text
120, 180, 240
240, 360, 480
```

### Log-Cos Scan

The spectrum is mapped into log-frequency space:

```text
ell = ln(f)
```

The tested model is:

```text
y(ell) = C + a*cos(k*ell) + b*sin(k*ell)
```

The scan searches over `k`, fits the best log-cos mode, and compares the real scan statistic against shuffled spectra.

Observed result:

```text
k = 44.5341
Delta chi2 = 21.139
p_shuffle_scanmax = 0.000999
```

## Interpretation

The artifact demonstrates three independent structures:

1. A discrete harmonic ladder.
2. Exact `Q = 2/3` ratio geometry.
3. Significant log-cos modulation in `ln(f)`.

The key interpretation is scale-free:

> The photodiode state does not merely preserve a single frequency. It preserves a ratio geometry.

The triplets:

```text
120 : 180 : 240
240 : 360 : 480
```

both reduce to:

```text
2 : 3 : 4
```

This means the same sideband geometry appears at different scale levels.

## Falsification

The claim is falsified or weakened if:

- virgin devices do not enter a structured state under the protocol,
- later captures show only broadband noise,
- the `2/3` and triplet geometry do not recur,
- the signal vanishes immediately after excitation,
- the same structure appears identically in controls without induction,
- the effect is fully explained by ordinary mains pickup under blinded controls.

## Controls and Next Steps

Recommended next tests:

- Run multiple virgin photodiodes.
- Include sham-excited controls.
- Include disconnected-channel controls.
- Include terminated-input oscilloscope controls.
- Repeat captures over multiple days.
- Repeat log-cos scans with fixed pre-declared parameters.
- Compare induced and non-induced devices under identical shielding and measurement conditions.

## Notes

The repository is organized as an artifact and reproduction package. The strongest use is to inspect the data, rerun the scripts, and compare the reported outputs against controls.


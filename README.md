# Optically Induced Harmonic State in a Silicon Photodiode

**Reproducible artifact package — prediction ledger, oscilloscope captures, FFT analysis, Q = 2/3 sideband geometry, and log-cos spectral scan with shuffle null.**

---

## Overview

This repository is a self-contained reproducibility package for a captured electrical waveform from a silicon photodiode following a brief optical excitation protocol.

The central measured fact is this: the post-excitation waveform contains a discrete 60 Hz-spaced harmonic ladder with exact Q = 2/3 ratio geometry across two detected triplets, and the log-frequency spectrum carries a strong cosine modulation under the canonical shuffle analysis (Δχ² = 21.139, p_shuffle_scanmax = 0.000999 over 1000 null trials).

All scripts are executable. All outputs are committed. The prediction ledger records what was expected before the capture.

The physical origin of the locked harmonic state is **unresolved**. Its 60 Hz spacing means mains-related or environmental coupling must remain in the alternative set, but the spacing alone does not identify the source. The archived runs also show that the high-amplitude harmonic lock is **condition-dependent rather than continuously present**.

The experimenter additionally reports that the harmonic state can be turned **on and off with the photodiode while the surrounding measurement setup is otherwise held fixed**. That observation points to **diode dependence** of the locked state rather than a continuously present oscilloscope/background signal alone. It is currently documented as retrospective experimental context; a dedicated labeled or blinded ON/OFF raw-control series is still required to quantify the causal separation directly from committed waveforms.

---

## Repository Contents

```
README.md                    — this file
proof_ledger.json            — pre-declared protocol and prediction record
scan_npy_logcos.py           — main analysis script (FFT → ratio → Koide → log-cos)
siglent.py / siglent_2.py    — oscilloscope capture scripts (Siglent SDS1104X HD)
siglent_proof.py             — capture + analysis pipeline
photodiode_proof.py          — induction protocol driver

docs/
  LOCKED_VS_NO_LOCK_RESULT_2026-08-31.md
                              — archived state-switch evidence and limitations
  DIODE_DEPENDENCE_CONTROL_2026-08-31.md
                              — diode ON/OFF observation and decisive replication protocol

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
| `koide_triplets.csv` | All ordered triplets with Q_low, Q_high, Q_mean, koide_error |
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
|---|---:|---:|
| 60 | 60.43 | 3 |
| 120 | 191.00 | 1 |
| 180 | 39.13 | 5 |
| 240 | 114.99 | 2 |
| 300 | 14.05 | 7 |
| 360 | 47.97 | 4 |
| 480 | 17.84 | 6 |
| 600 | 13.94 | 8 |

The ladder is 60 Hz-spaced throughout. This frequency coincidence makes 60 Hz-related electrical or environmental coupling an important control target, but **it is not a source attribution by itself**. The ratio and log-cos analyses are reproducible measurements of the artifact regardless of source.

### Archived locked vs no-lock contrast

The repository contains a later retained raw capture in which the harmonic lock is absent under the same audit pipeline:

| Run | Archived phase | RMS voltage | 60 Hz amplitude | Detected harmonics | Ladder statistic | State |
|---|---|---:|---:|---:|---:|---|
| `1773046275` | post | 0.0763359644 | 63.3030966 | 12 | 286.8951175 | **LOCKED** |
| `1773047118` | unknown | 0.0011999400 | 0.1081480 | 0 | 0 | **NO LOCK** |

The epoch identifiers are separated by **843 seconds (~14.1 minutes)**. In the later capture, the 60 Hz amplitude is approximately **585× smaller**, the detected harmonic count falls from **12 to 0**, and the ladder statistic falls from **286.9 to 0**.

This establishes a narrow but important result:

> **The high-amplitude 60 Hz-spaced harmonic lock is not continuously present across the retained captures. A later raw capture has no detected harmonic ladder under the same audit pipeline.**

This directly challenges a simplistic invariant-background model in which the full locked ladder is assumed to be continuously present in the measurement stream. It does **not** exclude intermittent or condition-dependent electrical, environmental, grounding, magnetic, or instrument coupling.

The original archive labels run `1773047118` as `phase: unknown`, and its device metadata are incomplete. Therefore it is not silently relabeled as a matched post-excitation or diode-OFF control.

See [`docs/LOCKED_VS_NO_LOCK_RESULT_2026-08-31.md`](docs/LOCKED_VS_NO_LOCK_RESULT_2026-08-31.md) and [`tables_r/statistical_audit/locked_vs_no_lock_comparison.csv`](tables_r/statistical_audit/locked_vs_no_lock_comparison.csv) for the explicit result and machine-readable values.

### Diode-dependent switching observation

The experimenter reports that the pronounced harmonic state can be turned **on and off with the photodiode**, with the surrounding oscilloscope/acquisition environment otherwise held fixed.

If reproduced in a dedicated matched control series, this would establish a stronger causal statement than the run-to-run contrast alone:

> **The probability of entering the locked harmonic state depends on the diode condition.**

That would substantially weaken a simple scope-only or continuously present ambient-background explanation. It would **not**, by itself, distinguish an unusual confined state from conventional diode physics. Ordinary alternatives would still include photodiode rectification, nonlinear junction response, diode-dependent impedance/capacitance, carrier trapping, leakage/surface-state effects, and diode-mediated mixing or amplification of an ambient 60 Hz field.

The current repository therefore treats the diode ON/OFF observation as **retrospective experimental evidence requiring a dedicated archived replication**, not as a completed blinded causal control.

See [`docs/DIODE_DEPENDENCE_CONTROL_2026-08-31.md`](docs/DIODE_DEPENDENCE_CONTROL_2026-08-31.md) for the claim boundary and the recommended ABAB/randomized protocol.

### Experimental context annotation

The experimenter reports that the photodiode session used improvised electromagnetic shielding (a makeshift Faraday enclosure), that the acquisition system was restarted, and that a subsequent repeat showed a delayed response without entering the earlier locked state. This history is recorded as **retrospective experimental context**, not as preregistered metadata. The existing archive does not independently map every element of that recollection to run `1773047118`, so the repository-supported direct claim remains the locked/no-lock state contrast above.

Improvised shielding also does not necessarily remove conducted coupling, grounding effects, or low-frequency magnetic fields. The appropriate conclusion is therefore that **mains-related coupling remains unexcluded, but is not established as the origin**.

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

Two triplets satisfy Q_low = Q_high = 2/3 exactly:

| Triplet (Hz) | Q_low | Q_high | Q_mean | koide_error |
|---|---:|---:|---:|---:|
| (120, 180, 240) | 2/3 | 2/3 | 2/3 | 0.0 |
| (240, 360, 480) | 2/3 | 2/3 | 2/3 | 0.0 |

Both triplets reduce to the integer ratio 2:3:4. They occur at two frequency scales in the same capture and share the 240 Hz peak, so they should not be treated as statistically independent observations.

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

where χ²_null is the variance-normalized residual of the constant-only fit.

### Null test

The canonical significance diagnostic shuffles the FFT amplitudes 1000 times, repeats the full k-scan on each shuffled spectrum, and records the scan-max Δχ². The reported p-value is the add-one empirical fraction of shuffled trials whose scan maximum meets or exceeds the observed value.

This is a **scan-max shuffle null**: it accounts for the look-elsewhere effect across the tested k range. It does not preserve every possible colored-noise or frequency-local correlation and therefore should not be interpreted as a complete physical source test.

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

No shuffled spectrum exceeded the real scan statistic in this canonical run.

The p-value of 0.000999 is the measured scan-max shuffle result for this run and this analysis window. It is not a Gaussian sigma conversion and should not be generalized beyond the tested statistic, range, and sample.

Structured colored-noise and harmonic-ladder nulls included in the statistical audit test different model classes. Compatibility with a synthetic 60 Hz-harmonic null means that the log-cos statistic alone cannot distinguish that model class; it does **not** establish that physical mains pickup generated the recorded waveform.

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

The ledger establishes chronology: the reported target structure was specified before the canonical capture rather than reconstructed afterward. It does not substitute for physical controls.

The baseline run (epoch 1773046194) is included in `proof_runs/` for comparison, although its archived form is not a matched raw waveform under the same normalization as the canonical `.npy` capture.

---

## Controls and Falsification

The following controls are the direct path to distinguishing a condition-dependent device state from environmental, conducted, magnetic, grounding, or instrument coupling.

| Control | What it tests |
|---|---|
| **Predeclared diode ON/OFF ABAB or randomized sequence** | Does the probability of lock track diode condition while the rest of the acquisition setup remains fixed? |
| Virgin photodiode, excitation omitted | Does the harmonic structure appear without the protocol? |
| Sham-excited control (protocol run, diode dark or absent) | Does the excitation sequence itself drive the observation? |
| Disconnected input / terminated channel | Is the structure present in the measurement chain alone? |
| Dark enclosure, no optical path | Is ambient light driving the observation? |
| Scope channel swap | Is the signal channel-specific or instrument-wide? |
| Battery-powered or preamplifier-isolated supply | Does conducted supply coupling account for the ladder? |
| Shielding and grounding variation | Does the structure track environmental EM or grounding conditions? |
| Mains-frequency monitor | Does the observed frequency/phase track the local power waveform? |
| Repeat captures across devices and days | Is the effect reproducible across hardware instances and time? |
| Pre-declared parameter lock (fixed fmin, fmax, k range) | Are scan parameters chosen post-hoc to fit the result? |

The highest-value immediate control is the diode-state series. With all other settings fixed, repeated ON/OFF or randomized conditions should be captured raw and analyzed blind using a predeclared lock classifier. The primary causal comparison is `P(lock | diode ON)` versus `P(lock | diode OFF)`.

The induced-state interpretation would be strongly challenged if matched virgin, disconnected, terminated-input, or diode-OFF controls reproducibly generate the same locked spectrum under identical acquisition conditions. Conversely, repeatable diode-specific or protocol-specific locking with flat measurement-chain controls would substantially weaken simple background explanations.

All matched control runs should be processed through the same `scan_npy_logcos.py` pipeline and committed alongside the primary artifact so outputs are directly comparable.

---

## Interpretation

Three directly or computationally reproducible structures are present in the canonical artifact:

1. **A 60 Hz-spaced discrete harmonic ladder** across multiple resolved peaks.
2. **Exact Q = 2/3 sideband geometry** in two detected 2:3:4 triplets.
3. **A strong log-cos scan maximum** in ln(f), with the canonical shuffle analysis reaching its 1/1001 Monte Carlo floor.

The ratio geometry is scale-free: the triplets (120, 180, 240) and (240, 360, 480) encode the same 2:3:4 integer structure at different absolute frequencies. Because both arise from the same 60 Hz harmonic ladder and share 240 Hz, the ratio geometry should be treated as a structural property of that ladder, not as independent statistical confirmations.

The **physical source remains unresolved**. The 60 Hz spacing makes mains-related coupling a required alternative to test; it does not justify assigning the observation to mains. Separately, the archived `1773046275` → `1773047118` comparison demonstrates that the full harmonic lock is not continuously present across retained raw captures. That state dependence is evidence against a simple invariant-background description, while still leaving intermittent or condition-dependent coupling mechanisms open.

The experimenter's report that the state can be switched on and off with the diode further motivates a **diode-dependent state / nonlinear-transduction hypothesis**. If a predeclared matched ON/OFF series reproduces that control, the result would move beyond generic 60 Hz coincidence by demonstrating that the diode condition predicts whether the locked state appears. Conventional diode nonlinearity, rectification, impedance/capacitance changes, carrier trapping, and diode-mediated environmental mixing would still need to be separated from any more unusual mechanism.

The shielding/restart/delayed-no-lock and diode ON/OFF histories are retained as retrospective experimental context and are not used to overstate what the original machine-readable provenance establishes.

The repository therefore reports what was measured, what was predicted, what reproduces computationally, what changes across runs, what the experimenter reports about controllability, and which causal questions remain open.

---

## Citation / Author

**Repository:** `rickyjreyes/photodiode`  
**Capture instrument:** Siglent Technologies SDS1104X HD (firmware 5.5.6.1.1.0.2)  
**Primary capture:** 2026-03-09, epoch 1773046275

If citing this artifact package, reference the repository commit hash of the canonical run together with `summary.json` for exact scalar reproducibility.

> [Placeholder — author affiliation and formal citation to be added.]

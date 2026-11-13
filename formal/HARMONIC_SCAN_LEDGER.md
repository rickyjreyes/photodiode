# Proof ledger: photodiode #5 (harmonic ladder / ratios / log-frequency) and #6 (scan-max inference)

Toolchain: Lean `v4.23.0`, Mathlib `v4.23.0` (`37df177aaa770670452312393d4e84aaad56e7b6`, pinned in
`lake-manifest.json`). Build: `cd formal && lake exe cache get && lake build`.
Every declaration listed below is kernel-checked and depends only on `propext`,
`Classical.choice`, `Quot.sound` (see the axiom audit in `.github/workflows/lean-protocol.yml`).

Status vocabulary: **PROVED** (kernel-checked), **REFUTED** (kernel-checked counterexample),
**CONDITIONAL** (proved under hypotheses not established for the apparatus), **OPEN**, **EMPIRICAL**
(outside Lean by design).

## #5 — harmonic ladder, ratio and log-frequency statistics

| Issue item | Source | Lean declaration(s) | Status / hypotheses |
|---|---|---|---|
| Sample-rate / bin conventions, FFT scaling | `compute_fft_spectrum` | `Spectral.hann`, `Spectral.dftRe/dftIm`, `Spectral.binFreq`, `Spectral.amplitude_smul` (pre-existing) | PROVED (exact-real model of `rfft`; NumPy agreement is a test, not a proof) |
| Frequency-index identities | `rfftfreq`, `detect_peaks` | `HarmonicLadder.binPosition_binFreq`, `onGrid_ladder`, `binPosition_rung`, `HarmonicRatios.binFreq_multiple/step` | PROVED (`N ≠ 0`, `fs ≠ 0`) |
| Harmonic ladder ratio identities | ladder `n·60 Hz` | `rung_ratio`, `rung_succ_sub`, `log_rung`, `log_gap` | PROVED (`f0 > 0`, `n > 0`) |
| Ladder vs log-periodicity | `build_log_domain` (`ell = ln f`) | `log_gap_strict_anti`, `harmonic_ladder_not_log_uniform` | PROVED: a harmonic ladder has **no** common log-period |
| Detected-triplet ratios | `koide_triplets` on `freq[peaks]` | `lowQ_binFreq`, `highQ_binFreq`, `binTriplet_exact_two_thirds_iff` | PROVED: detected ratios are bin-index ratios; exact 2/3 ⇔ `3k1 = 2k2 ∧ 3k3 = 4k2` |
| Koide-like identities, order | `koide_triplets` (sorts first) | `exact_two_thirds_implies_double_extremes` (pre-existing), `sorted_triplet_bounds`, `triplet_order_dependent` | PROVED; sorting is load-bearing (unsorted input changes the statistic) |
| Zero denominators | `Q_low`, `Q_high` | `lowQ_zero_denominator`, `highQ_zero_denominator` (pre-existing) | PROVED (Lean's `x/0 = 0` is a convention, not a physical ratio) |
| Finite precision / tolerance | `koide_error` | `approx_two_thirds_bound`: `|f3 − 2f1| ≤ 4 ε f2`; `koideErrorSquared_eq_zero_iff`, `koideErrorSquared_common_scale` | PROVED (`f2 > 0`) |
| Log-cos objective and invariances | `fit_logcos_cpu` | `logCosRSS_offset` (pre-existing), `logCosModel_shift`, `logCosRSS_shift`, `logCosRSS_unit_change`, `rot_amplitude`, `logCosRSS_scale` | PROVED: unit change = phase rotation preserving amplitude; ordinate scale law `λ²` |
| Active winding | `active_winding` | `activeWinding_scale`, `activeWinding_phase` | PROVED (`fmin, fmax > 0`) |
| Failure mode: harmonic indexing | nearest-bin bookkeeping | `round_harmonic_defect_bound`, `nearest_bin_not_harmonic`, `nearest_bin_witness_in_range` | REFUTED ("nearest bin of n·f = n × nearest bin of f" is false; witness inside the one-sided range) |
| Failure mode: thresholding | `detect_peaks` height gate | `passesHeight_gain` (PROVED), `passesHeight_offset_counterexample` (REFUTED offset invariance) | — |
| Failure mode: preprocessing gain | `log1p(amp)` | `log1p_gain_not_constant_offset` | REFUTED ("gain = constant log offset") |
| Released CSV/NPY outputs cross-check | analyzer outputs | `tests/test_harmonic_ladder_scanmax_contracts.py`, `tests/test_harmonic_ratio_contracts.py` | EMPIRICAL conformance (13 + 4 tests call the real functions) |
| Physical / WCT origin of 60 Hz-spaced peaks | — | none | Not a mathematical claim; nothing here supports an optical or WCT cause |

## #6 — scan-max permutation inference and decision gates

| Issue item | Source | Lean declaration(s) | Status / hypotheses |
|---|---|---|---|
| Exchangeability, permutation action, scan-max statistic | `shuffle_null`, `scan_logcos_cpu` | `ExchangeableModel` (pre-existing), `ScanMax.scanMax`, `le_scanMax`, `scanMax_attained`, `scanMax_mono` | PROVED |
| Rank/counting basis of a valid p-value | `p_value_geq` | `card_small_upper_tail_le`, `randomization_level`, `mc_pvalue_eq_tail`, `exchangeable_level` (pre-existing) | PROVED under the stated exchangeability / random-assignment hypothesis (CONDITIONAL for the apparatus) |
| Scan-max p-value validity | `shuffle_null` stores per-permutation scan maxima | `scanMax_randomization_level` | PROVED (scan max is a single fixed statistic) |
| Why a single-bin p-value is insufficient | best-`k` reporting | `selected_pValue_le_scan` (selected-bin p ≤ scan p, always), `minP_inflation` / `minP_exceeds_level` (null rejection prob. `3/4 > 1/2`) with `stat_marginally_valid` | PROVED / REFUTED (min-p validity) |
| Decision rule for matched ON/OFF captures, missing metadata | protocol | `Decision.classify`, `classify_pass_iff`, `missing_group_inconclusive`, `manifest_violation_inconclusive`, `too_many_missing_inconclusive`, `pass_level` (pre-existing) | PROVED (predeclared rule; INCONCLUSIVE on incomplete metadata) |
| Broken exchangeability invalidates inference | bin permutation of windowed spectra | `exchangeable_swap_symmetric`, `not_exchangeable_of_strict_order`, `mcPValue_two_strict` | REFUTED (validity without exchangeability): p ≡ 1/2 surely |
| Uncontrolled preprocessing / post-hoc selection | data-driven bin choice | `selected_pValue_le_scan`, `minP_inflation` | REFUTED (post-hoc single-bin inference) |
| Regression vs frozen scan configuration | `shuffle_null` replay | `test_shuffle_null_stores_scan_maxima_and_selection_effect` | EMPIRICAL conformance |
| Diode-exclusive causality, persistence, Q-factor | — | none | Out of scope: a valid p-value is not a causal or physical claim |

## Remaining obligations

* Exchangeability of permuted Hann-windowed spectral bins under the physical null is **not** established
  (and is generally false for correlated neighbouring bins); the device-level randomization test in
  `Decision.lean` is the primary inferential route.
* The Lean model of `rfft` is exact-real; floating-point agreement with NumPy is tested, not proved.


## Additional source-conformance checks after PR #8 initial build

`tests/test_harmonic_scan_exact_orbit.py` adds a full finite-orbit
enumeration of scan maxima over a predeclared frequency grid using the
real `scan_logcos_cpu` function, including the finite-sample
superuniformity inequality for a synthetic uniformly sampled orbit.
It also pins the two-test min-p familywise failure, a strictly ordered
nonexchangeable null, and the one-sided Nyquist convention for odd/even N.

**These added tests require a fresh local run.** The previously reported
7,222-job Lean build applies to the unchanged theorem modules only; it
does not establish empirical exchangeability or the correctness of a
device-level randomization schedule. This issue remains review-open
until source/test evidence is accepted.

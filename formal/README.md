# Photodiode analysis protocol: formal specification (Lean 4)

This directory specifies the finite-data analysis and decision procedure used to evaluate a
persistent post-excitation spectral signal, and proves its correctness properties in Lean
(`v4.23.0`, Mathlib `v4.23.0`).

> **Lean cannot prove that a physical effect exists.** What is verified is that the frozen
> protocol implements its declared mathematical tests, and that the decision rule has
> type-I error at most `α` *under explicitly stated hypotheses*. A PASS does **not** establish
> a novel optical-memory mechanism.

Build: `cd formal && lake exe cache get && lake build`.

## Modules and theorems

### `PhotodiodeProtocol/Randomization.lean`
| Theorem | Statement | Hypotheses |
|---|---|---|
| `card_small_upper_tail_le` | At most `k` reference points have upper-tail count `≤ k` (ties count toward the tail) | none |
| `randomization_level` | `P_unif(p ≤ α) ≤ α` for the exact p-value `#{b : T b ≥ T a}/|A|` | `0 ≤ α`, `A ≠ ∅`, uniform assignment on `A`, sharp null |
| `randomization_negative_level` | The rejection event is empty for `α < 0` | none |
| `pValue_pos`, `pValue_le_one` | `0 < p ≤ 1` | `a ∈ A` |
| `mc_pvalue_eq_tail` | `scan_npy_logcos.p_value_geq` formula `(1+#{null ≥ real})/(1+B)` equals the upper-tail rank of the observed value | none |
| `exchangeable_level` | `P(p ≤ α) ≤ α` when `(X_0,…,X_B)` is exchangeable | `ExchangeableModel` (weight-preserving permutation action) |
| `intersection_union_level` | Requiring both comparisons to reject keeps level `α` without multiplicity correction | as `randomization_level` |
| `card_perm_fixed_image` | `#{σ ∈ S_{n+1} : σ 0 = i} = n!` | none |

### `PhotodiodeProtocol/Spectral.lean`
This is an exact real-arithmetic model of the Hann window, `rfft`, amplitude, `rfftfreq` and band
selection, plus the declared normalized band score and persistence score.
| Theorem | Statement |
|---|---|
| `amplitude_smul` | `|Y_{cy}[k]| = |c|·|Y_y[k]|` |
| `normScore_gain_invariant` | Normalized score is invariant under any nonzero detector gain |
| `persistence_gain_invariant` | Persistence score is invariant under independent nonzero gains on baseline and late captures |
| `bins_subset_rfft` | Selected bins are valid one-sided bins |

### `PhotodiodeProtocol/Decision.lean`
| Theorem | Statement |
|---|---|
| `classify_pass_iff`, `classify_inconclusive_iff` | Exact characterisation of PASS and INCONCLUSIVE |
| `protocol_mismatch_inconclusive` | A bundle declaring any other protocol id (e.g. changed thresholds) is INCONCLUSIVE |
| `missing_group_inconclusive`, `absent_sham_inconclusive` | An empty or under-filled group (including controls) is INCONCLUSIVE |
| `manifest_violation_inconclusive` | A present capture with an incomplete manifest, wrong code revision, or timestamp on the wrong side of excitation is INCONCLUSIVE |
| `duplicate_label_inconclusive` | Non-unique blinded labels give INCONCLUSIVE |
| `too_many_missing_inconclusive` | More than `maxMissing` incomplete devices gives INCONCLUSIVE |
| `score_depends_only_on_captures` | Scores never read the assignment (blinding) |
| `refSet_eq_of_mem` | The conditional reference set is the same from every assignment in it |
| **`pass_level`** | **`P(classify = PASS) ≤ α`** under the sharp null with uniform assignment on the conditional reference set |

## Frozen protocol (data model)
* Device: unique blinded label, excitation time, optional baseline and late captures.
* Capture: `N` samples plus a manifest: acquisition timestamp, hardware id, firmware, environment,
  analysis-code revision.
* `FrozenProtocol`: protocol id, code revision, `N`, target band, reference band, `α ≥ 0`, `δ`,
  `nMin ≥ 1`, `maxMissing`. It is a parameter of `classify`, never a field of the evidence bundle.
* Missing data: an absent capture is missing (complete-case analysis). A present but invalid
  capture is a violation, and gives INCONCLUSIVE.
* Statistic: induced mean minus control mean of per-device persistence scores, computed
  separately for sham and virgin controls.
* PASS: both exact p-values `≤ α` **and** both effects `≥ δ`. Otherwise FAIL, if every gate holds.

## Explicit statistical assumptions
1. **Sharp null.** Each device's captures would be identical under any assignment.
2. **Random assignment.** The realised assignment is uniform on the conditional reference set
   (same labels off the complete set, same group counts on it). This holds for complete
   randomization with fixed group sizes when missingness is a property of the device that does
   not depend on the assignment. That implication is stated, **not** proved here.
3. The repository's existing within-waveform `shuffle_null` is valid (`exchangeable_level`) only
   if spectral bins are exchangeable under the null. Hann-windowed physical spectra generally
   violate this (`1/f` trends, leakage correlations). It tests "no frequency structure", not
   persistence versus controls, so it is **not** the primary test.

## Remaining obligations (open)
* Floating-point conformance between `scan_npy_logcos.py` and the exact model in `Spectral.lean`.
* A proof that complete randomization plus assignment-independent missingness implies
  assumption 2.
* A frozen numeric protocol instance (band edges, `α`, `δ`, `nMin`, `maxMissing`), committed
  before any replication data are collected.
* Real randomized replication experiments. The present repository data are not a randomized
  induced/sham/virgin design.

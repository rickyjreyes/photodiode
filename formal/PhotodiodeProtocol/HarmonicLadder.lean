import Mathlib
import PhotodiodeProtocol.HarmonicRatios

/-!
# Harmonic ladder, nearest-bin indexing, ordered triplets and log-frequency invariances (#5)

Sources (exact functions):

* `scan_npy_logcos.compute_fft_spectrum` — bin frequencies `k·fs/N` (`numpy.fft.rfftfreq`);
* `scan_npy_logcos.detect_peaks` — height gate `amp ≥ peak_height_frac · max(amp)`, and the
  reported peak frequencies are *bin-centre* frequencies `freq[peaks]`;
* `scan_npy_logcos.koide_triplets` — frequency-sorted triplets `f1 ≤ f2 ≤ f3` with
  `Q_low = f1/f2`, `Q_high = f3/(2 f2)`, `koide_error = sqrt((Q_low-2/3)²+(Q_high-2/3)²)`;
* `scan_npy_logcos.build_log_domain` / `fit_logcos_cpu` — `ell = ln f`,
  `y = log1p(amp)` detrended and standardized, model `C + a cos(k ell) + b sin(k ell)`;
* `scan_npy_logcos.active_winding` — `k ln(fmax/fmin) / (2π)`.

Everything below is exact real arithmetic. Nothing here asserts that a spectral peak has an
optical, electrical-mains or WCT origin; the identities only constrain what the *statistics*
can mean.

## Main results

* `rung_ratio`, `log_rung`, `log_gap`: harmonic ladder `n·f0`, its ratios and log-frequencies.
* `harmonic_ladder_not_log_uniform`: a harmonic ladder is **not** equally spaced in `ln f`
  (log gaps strictly decrease), so it has no single log-period.
* `lowQ_binFreq`, `highQ_binFreq`, `binTriplet_exact_two_thirds_iff`: because detected peak
  frequencies are bin centres, detected triplet ratios are ratios of bin indices; an exact
  two-thirds triplet on the grid is exactly `3·k1 = 2·k2` and `3·k3 = 4·k2`.
* `round_harmonic_defect_bound`, `nearest_bin_not_harmonic`: nearest-bin indexing does not
  commute with harmonic multiplication (failure mode for "index × n" harmonic bookkeeping),
  with the sharp-in-form bound `|round(n x) - n round x| ≤ (n+1)/2`.
* `sorted_triplet_bounds`, `triplet_order_dependent`: frequency sorting is load-bearing.
* `koideErrorSquared_eq_zero_iff`, `koideErrorSquared_common_scale`,
  `approx_two_thirds_bound`: exact and tolerance versions of the two-thirds criterion.
* `logCosRSS_shift`, `logCosRSS_unit_change`, `logCosRSS_scale`: change of frequency units
  is a phase rotation of `(a, b)` preserving the amplitude `sqrt(a²+b²)`; ordinate scaling
  multiplies the residual sum of squares by `λ²` (removed by the source's standardization).
* `passesHeight_gain`, `passesHeight_offset_counterexample`,
  `log1p_gain_not_constant_offset`: the relative peak-height gate is gain invariant but not
  offset invariant; `log1p` turns a detector gain change into a *non-constant* ordinate change.
-/

namespace PhotodiodeProtocol.HarmonicLadder

open PhotodiodeProtocol.Spectral PhotodiodeProtocol.HarmonicRatios

/-! ## 1. The harmonic ladder and its log-frequencies -/

/-- The `n`-th rung of a harmonic ladder with fundamental `f0` (the source's 60 Hz spacing). -/
def rung (f0 : ℝ) (n : ℕ) : ℝ := n * f0

theorem rung_ratio (f0 : ℝ) (hf0 : f0 ≠ 0) (m n : ℕ) (hn : n ≠ 0) :
    rung f0 m / rung f0 n = (m : ℝ) / n := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  unfold rung
  field_simp

theorem rung_succ_sub (f0 : ℝ) (n : ℕ) : rung f0 (n + 1) - rung f0 n = f0 := by
  unfold rung
  push_cast
  ring

theorem log_rung (f0 : ℝ) (hf0 : 0 < f0) (n : ℕ) (hn : 0 < n) :
    Real.log (rung f0 n) = Real.log n + Real.log f0 := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  unfold rung
  rw [Real.log_mul hn' hf0.ne']

/-- Consecutive rungs are separated in `ln f` by `ln((n+1)/n)`. -/
theorem log_gap (f0 : ℝ) (hf0 : 0 < f0) (n : ℕ) (hn : 0 < n) :
    Real.log (rung f0 (n + 1)) - Real.log (rung f0 n) = Real.log (((n : ℝ) + 1) / n) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  rw [log_rung f0 hf0 (n + 1) (Nat.succ_pos n), log_rung f0 hf0 n hn,
    Real.log_div (by positivity) hn'.ne']
  push_cast
  ring

/-- The log gaps of a harmonic ladder strictly decrease. -/
theorem log_gap_strict_anti (n : ℕ) (hn : 0 < n) :
    Real.log (((n : ℝ) + 2) / (n + 1)) < Real.log (((n : ℝ) + 1) / n) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  apply Real.log_lt_log (by positivity)
  rw [div_lt_div_iff₀ (by positivity) hn']
  nlinarith

/-- A harmonic ladder has no common log-period: it is not uniformly spaced in `ln f`. -/
theorem harmonic_ladder_not_log_uniform (f0 : ℝ) (hf0 : 0 < f0) :
    ¬ ∃ P : ℝ, ∀ n : ℕ, 0 < n → Real.log (rung f0 (n + 1)) - Real.log (rung f0 n) = P := by
  rintro ⟨P, hP⟩
  have h1 := hP 1 one_pos
  have h2 := hP 2 two_pos
  rw [log_gap f0 hf0 1 one_pos] at h1
  rw [log_gap f0 hf0 2 two_pos] at h2
  have h := log_gap_strict_anti 1 one_pos
  norm_num at h1 h2 h
  linarith

/-! ## 2. Detected triplets live on the FFT grid -/

/-- Detected-peak ratios are ratios of bin indices (`detect_peaks` reports `freq[peaks]`). -/
theorem lowQ_binFreq (N : ℕ) (fs : ℝ) (hN : N ≠ 0) (hfs : fs ≠ 0) (k1 k2 : ℕ) (hk2 : k2 ≠ 0) :
    lowQ (binFreq N fs k1) (binFreq N fs k2) = (k1 : ℝ) / k2 := by
  have hN' : (N : ℝ) ≠ 0 := by exact_mod_cast hN
  have hk2' : (k2 : ℝ) ≠ 0 := by exact_mod_cast hk2
  unfold lowQ binFreq
  field_simp

theorem highQ_binFreq (N : ℕ) (fs : ℝ) (hN : N ≠ 0) (hfs : fs ≠ 0) (k2 k3 : ℕ) (hk2 : k2 ≠ 0) :
    highQ (binFreq N fs k2) (binFreq N fs k3) = (k3 : ℝ) / (2 * k2) := by
  have hN' : (N : ℝ) ≠ 0 := by exact_mod_cast hN
  have hk2' : (k2 : ℝ) ≠ 0 := by exact_mod_cast hk2
  unfold highQ binFreq
  field_simp

/-- An exact two-thirds triplet of detected (bin-centre) peaks is the integer condition
`3·k1 = 2·k2 ∧ 3·k3 = 4·k2`; e.g. bins `(2,3,4)`, never the reference `(1,2,3)` ladder. -/
theorem binTriplet_exact_two_thirds_iff (N : ℕ) (fs : ℝ) (hN : N ≠ 0) (hfs : fs ≠ 0)
    (k1 k2 k3 : ℕ) (hk2 : k2 ≠ 0) :
    (lowQ (binFreq N fs k1) (binFreq N fs k2) = 2 / 3 ∧
        highQ (binFreq N fs k2) (binFreq N fs k3) = 2 / 3) ↔
      (3 * k1 = 2 * k2 ∧ 3 * k3 = 4 * k2) := by
  have hk2' : (k2 : ℝ) ≠ 0 := by exact_mod_cast hk2
  rw [lowQ_binFreq N fs hN hfs k1 k2 hk2, highQ_binFreq N fs hN hfs k2 k3 hk2]
  constructor
  · rintro ⟨h1, h3⟩
    rw [div_eq_iff hk2'] at h1
    rw [div_eq_iff (mul_ne_zero two_ne_zero hk2')] at h3
    constructor
    · have : (3 * k1 : ℝ) = 2 * k2 := by linarith
      exact_mod_cast this
    · have : (3 * k3 : ℝ) = 4 * k2 := by linarith
      exact_mod_cast this
  · rintro ⟨h1, h3⟩
    have h1' : (3 * k1 : ℝ) = 2 * k2 := by exact_mod_cast h1
    have h3' : (3 * k3 : ℝ) = 4 * k2 := by exact_mod_cast h3
    constructor
    · rw [div_eq_iff hk2']
      linarith
    · rw [div_eq_iff (mul_ne_zero two_ne_zero hk2')]
      linarith

/-! ## 3. Nearest-bin indexing does not commute with harmonic multiplication -/

/-- Fractional FFT-bin position of frequency `f`: `f·N/fs`. -/
noncomputable def binPosition (N : ℕ) (fs f : ℝ) : ℝ := f * N / fs

/-- Nearest FFT bin (round half up). -/
noncomputable def nearestBin (N : ℕ) (fs f : ℝ) : ℤ := round (binPosition N fs f)

theorem binPosition_binFreq (N : ℕ) (fs : ℝ) (hN : N ≠ 0) (hfs : fs ≠ 0) (k : ℕ) :
    binPosition N fs (binFreq N fs k) = k := by
  have hN' : (N : ℝ) ≠ 0 := by exact_mod_cast hN
  unfold binPosition binFreq
  field_simp

/-- If the fundamental is exactly on bin `b0`, the whole ladder is on the grid. -/
theorem onGrid_ladder (N : ℕ) (fs : ℝ) (b0 n : ℕ) :
    rung (binFreq N fs b0) n = binFreq N fs (n * b0) := by
  unfold rung
  rw [binFreq_multiple]

theorem binPosition_rung (N : ℕ) (fs f0 : ℝ) (n : ℕ) :
    binPosition N fs (rung f0 n) = n * binPosition N fs f0 := by
  unfold binPosition rung
  ring

/-- General defect bound: `|round(n·x) − n·round(x)| ≤ (n+1)/2`. -/
theorem round_harmonic_defect_bound (x : ℝ) (n : ℕ) :
    |((round ((n : ℝ) * x) : ℤ) : ℝ) - n * ((round x : ℤ) : ℝ)| ≤ ((n : ℝ) + 1) / 2 := by
  have h1 := abs_sub_round ((n : ℝ) * x)
  have h2 := abs_sub_round x
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  obtain ⟨a1, b1⟩ := abs_le.mp h1
  obtain ⟨a2, b2⟩ := abs_le.mp h2
  rw [abs_le]
  constructor
  · nlinarith [mul_le_mul_of_nonneg_left b2 hn]
  · nlinarith [mul_le_mul_of_nonneg_left a2 hn]

theorem round_seven_fifths : round ((7 : ℝ) / 5) = 1 := by
  rw [round_eq, Int.floor_eq_iff]
  norm_num

theorem round_twentyone_fifths : round ((21 : ℝ) / 5) = 4 := by
  rw [round_eq, Int.floor_eq_iff]
  norm_num

/-- Concrete failure mode: with `N = 70`, `fs = 3000` (Nyquist 1500 Hz, one-sided bins
`0..35`), the 60 Hz fundamental sits at bin position `1.4` (nearest bin 1) but its third
harmonic 180 Hz sits at `4.2` (nearest bin 4), not at `3 × 1 = 3`. -/
theorem nearest_bin_not_harmonic :
    nearestBin 70 3000 (rung 60 3) ≠ 3 * nearestBin 70 3000 60 := by
  have e1 : binPosition 70 3000 (rung 60 3) = (21 : ℝ) / 5 := by
    norm_num [binPosition, rung]
  have e2 : binPosition 70 3000 60 = (7 : ℝ) / 5 := by
    norm_num [binPosition]
  unfold nearestBin
  rw [e1, e2, round_seven_fifths, round_twentyone_fifths]
  norm_num

/-- Both frequencies of the witness lie inside the one-sided `rfft` range `[0, N/2]`. -/
theorem nearest_bin_witness_in_range :
    binPosition 70 3000 (rung 60 3) ≤ (70 : ℝ) / 2 ∧ binPosition 70 3000 60 ≤ (70 : ℝ) / 2 := by
  constructor <;> norm_num [binPosition, rung]

/-! ## 4. Ordered triplets and the two-thirds criterion -/

/-- For frequency-sorted positive triplets (the source sorts before `combinations`). -/
theorem sorted_triplet_bounds (f1 f2 f3 : ℝ) (h0 : 0 < f1) (h12 : f1 ≤ f2) (h23 : f2 ≤ f3) :
    0 < lowQ f1 f2 ∧ lowQ f1 f2 ≤ 1 ∧ (1 / 2 : ℝ) ≤ highQ f2 f3 := by
  have hf2 : 0 < f2 := lt_of_lt_of_le h0 h12
  refine ⟨div_pos h0 hf2, (div_le_one hf2).mpr h12, ?_⟩
  unfold highQ
  rw [le_div_iff₀ (by positivity)]
  linarith

/-- The statistics are not permutation invariant: unsorted input gives a different value. -/
theorem triplet_order_dependent :
    lowQ (240 : ℝ) 120 ≠ lowQ (120 : ℝ) 240 ∧ highQ (120 : ℝ) 360 ≠ highQ (240 : ℝ) 360 := by
  constructor <;> norm_num [lowQ, highQ]

theorem highQ_common_scale (f2 f3 c : ℝ) (hc : c ≠ 0) :
    highQ (c * f2) (c * f3) = highQ f2 f3 := by
  unfold highQ
  rw [show 2 * (c * f2) = c * (2 * f2) by ring]
  exact mul_div_mul_left _ _ hc

/-- Unit changes (e.g. Hz → kHz) do not change the Koide-style mismatch. -/
theorem koideErrorSquared_common_scale (f1 f2 f3 c : ℝ) (hc : c ≠ 0) :
    koideErrorSquared (c * f1) (c * f2) (c * f3) = koideErrorSquared f1 f2 f3 := by
  unfold koideErrorSquared
  rw [lowQ_common_scale _ _ _ hc, highQ_common_scale _ _ _ hc]

theorem koideErrorSquared_nonneg (f1 f2 f3 : ℝ) : 0 ≤ koideErrorSquared f1 f2 f3 := by
  unfold koideErrorSquared
  positivity

theorem koideErrorSquared_eq_zero_iff (f1 f2 f3 : ℝ) :
    koideErrorSquared f1 f2 f3 = 0 ↔ lowQ f1 f2 = 2 / 3 ∧ highQ f2 f3 = 2 / 3 := by
  unfold koideErrorSquared
  rw [add_eq_zero_iff_of_nonneg (sq_nonneg _) (sq_nonneg _), sq_eq_zero_iff, sq_eq_zero_iff,
    sub_eq_zero, sub_eq_zero]

/-- **Tolerance version.** If both ratios are within `ε` of `2/3`, then `f3` is within
`4 ε f2` of `2 f1`. (Finite-precision counterpart of `exact_two_thirds_implies_double_extremes`.) -/
theorem approx_two_thirds_bound (f1 f2 f3 ε : ℝ) (hf2 : 0 < f2)
    (hlo : |lowQ f1 f2 - 2 / 3| ≤ ε) (hhi : |highQ f2 f3 - 2 / 3| ≤ ε) :
    |f3 - 2 * f1| ≤ 4 * ε * f2 := by
  unfold lowQ at hlo
  unfold highQ at hhi
  have key : f3 - 2 * f1 =
      2 * f2 * (f3 / (2 * f2) - 2 / 3) - 2 * f2 * (f1 / f2 - 2 / 3) := by
    field_simp
    ring
  have h2f2 : (0 : ℝ) < 2 * f2 := by positivity
  rw [key]
  calc |2 * f2 * (f3 / (2 * f2) - 2 / 3) - 2 * f2 * (f1 / f2 - 2 / 3)|
      ≤ |2 * f2 * (f3 / (2 * f2) - 2 / 3)| + |2 * f2 * (f1 / f2 - 2 / 3)| := abs_sub _ _
    _ = 2 * f2 * |f3 / (2 * f2) - 2 / 3| + 2 * f2 * |f1 / f2 - 2 / 3| := by
        rw [abs_mul (2 * f2), abs_mul (2 * f2), abs_of_pos h2f2]
    _ ≤ 2 * f2 * ε + 2 * f2 * ε := by gcongr
    _ = 4 * ε * f2 := by ring

/-! ## 5. Invariances of the log-cos objective -/

/-- Rotation of the `(cos, sin)` coefficient pair by angle `θ`. -/
noncomputable def rotA (a b θ : ℝ) : ℝ := a * Real.cos θ + b * Real.sin θ

noncomputable def rotB (a b θ : ℝ) : ℝ := b * Real.cos θ - a * Real.sin θ

theorem logCosModel_shift (ell s k C a b : ℝ) :
    logCosModel (ell + s) k C a b =
      logCosModel ell k C (rotA a b (k * s)) (rotB a b (k * s)) := by
  unfold logCosModel rotA rotB
  rw [show k * (ell + s) = k * ell + k * s by ring, Real.cos_add, Real.sin_add]
  ring

theorem rot_amplitude (a b θ : ℝ) : rotA a b θ ^ 2 + rotB a b θ ^ 2 = a ^ 2 + b ^ 2 := by
  unfold rotA rotB
  linear_combination (a ^ 2 + b ^ 2) * Real.sin_sq_add_cos_sq θ

theorem rot_inverse (a b θ : ℝ) :
    rotA (rotA a b (-θ)) (rotB a b (-θ)) θ = a ∧ rotB (rotA a b (-θ)) (rotB a b (-θ)) θ = b := by
  unfold rotA rotB
  simp only [Real.cos_neg, Real.sin_neg]
  constructor
  · linear_combination a * Real.sin_sq_add_cos_sq θ
  · linear_combination b * Real.sin_sq_add_cos_sq θ

/-- Shifting all log-frequencies by `s` is the same objective with rotated `(a, b)`. -/
theorem logCosRSS_shift {N : ℕ} (ell y : Fin N → ℝ) (s k C a b : ℝ) :
    logCosRSS (fun i => ell i + s) y k C a b =
      logCosRSS ell y k C (rotA a b (k * s)) (rotB a b (k * s)) := by
  unfold logCosRSS
  simp only [logCosModel_shift]

/-- **Frequency-unit invariance.** Rescaling all (positive) frequencies by `c > 0` maps every
parameter of the log-cos fit to a parameter with identical residual sum of squares and the
same amplitude `a² + b²` (only the phase changes). -/
theorem logCosRSS_unit_change {N : ℕ} (f y : Fin N → ℝ) (hf : ∀ i, 0 < f i) (c : ℝ)
    (hc : 0 < c) (k C a b : ℝ) :
    ∃ a' b', logCosRSS (fun i => Real.log (c * f i)) y k C a' b' =
        logCosRSS (fun i => Real.log (f i)) y k C a b ∧ a' ^ 2 + b' ^ 2 = a ^ 2 + b ^ 2 := by
  set s := Real.log c
  have hlog : (fun i => Real.log (c * f i)) = fun i => Real.log (f i) + s := by
    funext i
    rw [Real.log_mul hc.ne' (hf i).ne']
    ring
  refine ⟨rotA a b (-(k * s)), rotB a b (-(k * s)), ?_, ?_⟩
  · rw [hlog, logCosRSS_shift]
    obtain ⟨h1, h2⟩ := rot_inverse a b (k * s)
    rw [h1, h2]
  · exact rot_amplitude a b _

/-- Ordinate scaling multiplies the objective by `λ²` (the source then divides by `std(y)`). -/
theorem logCosRSS_scale {N : ℕ} (ell y : Fin N → ℝ) (lam k C a b : ℝ) :
    logCosRSS ell (fun i => lam * y i) k (lam * C) (lam * a) (lam * b) =
      lam ^ 2 * logCosRSS ell y k C a b := by
  unfold logCosRSS logCosModel
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  ring

/-- `active_winding(k, fmin, fmax) = k ln(fmax/fmin)/(2π)`. -/
noncomputable def activeWinding (k fmin fmax : ℝ) : ℝ :=
  k * Real.log (fmax / fmin) / (2 * Real.pi)

theorem activeWinding_scale (k fmin fmax c : ℝ) (hc : c ≠ 0) :
    activeWinding k (c * fmin) (c * fmax) = activeWinding k fmin fmax := by
  unfold activeWinding
  rw [mul_div_mul_left _ _ hc]

/-- The phase `k·ell` advances by exactly `2π · active_winding` across the band. -/
theorem activeWinding_phase (k fmin fmax : ℝ) (h1 : 0 < fmin) (h2 : 0 < fmax) :
    k * Real.log fmax - k * Real.log fmin = 2 * Real.pi * activeWinding k fmin fmax := by
  unfold activeWinding
  rw [Real.log_div h2.ne' h1.ne']
  field_simp

/-! ## 6. Thresholding and preprocessing failure modes -/

/-- `detect_peaks` height gate: `frac · max(amp) ≤ amp`. -/
def passesHeight (frac maxAmp a : ℝ) : Prop := frac * maxAmp ≤ a

/-- The relative height gate is invariant under a positive detector gain. -/
theorem passesHeight_gain (frac maxAmp a c : ℝ) (hc : 0 < c) :
    passesHeight frac (c * maxAmp) (c * a) ↔ passesHeight frac maxAmp a := by
  unfold passesHeight
  rw [show frac * (c * maxAmp) = c * (frac * maxAmp) by ring]
  exact mul_le_mul_iff_of_pos_left hc

/-- …but it is **not** invariant under an additive noise-floor offset. -/
theorem passesHeight_offset_counterexample :
    ¬ passesHeight (1 / 2) 2 (9 / 10) ∧ passesHeight (1 / 2) (2 + 2) (9 / 10 + 2) := by
  unfold passesHeight
  norm_num

/-- A gain change before `log1p` is not a constant ordinate offset: it depends on amplitude,
so gain invariance of the *log-domain* pipeline is not automatic (contrast
`Spectral.normScore_gain_invariant`). -/
theorem log1p_gain_not_constant_offset :
    Real.log (1 + 2 * 1) - Real.log (1 + 1) ≠ Real.log (1 + 2 * 3) - Real.log (1 + 3) := by
  rw [← Real.log_div (by norm_num) (by norm_num), ← Real.log_div (by norm_num) (by norm_num)]
  exact ne_of_lt (Real.log_lt_log (by norm_num) (by norm_num))

end PhotodiodeProtocol.HarmonicLadder

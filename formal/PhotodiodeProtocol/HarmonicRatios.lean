import Mathlib
import PhotodiodeProtocol.Spectral

/-!
# Harmonic bin, ordered-triplet and log-cos algebra

Source: scan_npy_logcos.py:compute_fft_spectrum, make_ratio_pairs,
koide_triplets, fit_logcos_cpu and photodiode_proof.py harmonic summary.

These are exact arithmetic identities only. FFT bin frequencies are not
asserted to coincide with detected peak locations, and the separate
120/240/360 comparison is not identical to the two-thirds triplet target.
-/
namespace PhotodiodeProtocol.HarmonicRatios

open PhotodiodeProtocol.Spectral

/-- The Fourier bin frequency is k*sample_rate/N. -/
theorem binFreq_multiple (N m k : ℕ) (sampleRate : ℝ) :
    binFreq N sampleRate (m * k) = (m : ℝ) * binFreq N sampleRate k := by
  simp only [binFreq, Nat.cast_mul]
  ring

/-- Exact bin increments; valid even for N=0 in Lean's total division
semantics, but physical acquisition requires N>0 and fs>0. -/
theorem binFreq_step (N k : ℕ) (sampleRate : ℝ) :
    binFreq N sampleRate (k + 1) - binFreq N sampleRate k =
      sampleRate / N := by
  simp only [binFreq, Nat.cast_add, Nat.cast_one]
  ring

/-- f1/f2 for strictly frequency-sorted input peaks. -/
noncomputable def lowQ (f1 f2 : ℝ) : ℝ := f1 / f2

/-- f3/(2*f2) for strictly frequency-sorted input peaks. -/
noncomputable def highQ (f2 f3 : ℝ) : ℝ := f3 / (2 * f2)

/-- Source's squared error before taking the square root. -/
noncomputable def koideErrorSquared (f1 f2 f3 : ℝ) : ℝ :=
  (lowQ f1 f2 - (2 / 3 : ℝ)) ^ 2 +
    (highQ f2 f3 - (2 / 3 : ℝ)) ^ 2

theorem lowQ_zero_denominator (f : ℝ) : lowQ f 0 = 0 := by
  simp [lowQ]

theorem highQ_zero_denominator (f : ℝ) : highQ 0 f = 0 := by
  simp [highQ]

theorem lowQ_common_scale (f1 f2 c : ℝ) (hc : c ≠ 0) :
    lowQ (c * f1) (c * f2) = lowQ f1 f2 := by
  unfold lowQ
  exact mul_div_mul_left _ _ hc

/-- Exact two-thirds triplet matching forces f3 = 2*f1, not a 1:2:3
harmonic ladder. Requires a nonzero middle frequency. -/
theorem exact_two_thirds_implies_double_extremes (f1 f2 f3 : ℝ)
    (hf2 : f2 ≠ 0)
    (hlo : lowQ f1 f2 = (2 / 3 : ℝ))
    (hhi : highQ f2 f3 = (2 / 3 : ℝ)) :
    f3 = 2 * f1 := by
  have hden : (2 : ℝ) * f2 ≠ 0 := mul_ne_zero (by norm_num) hf2
  have h1 : f1 = (2 / 3 : ℝ) * f2 :=
    (div_eq_iff hf2).mp (by simpa [lowQ] using hlo)
  have h3 : f3 = (2 / 3 : ℝ) * (2 * f2) :=
    (div_eq_iff hden).mp (by simpa [highQ] using hhi)
  linarith

/-- The hard-coded (120,240,360) reference is a *different* criterion:
its exact source ratios are 1/2 and 3/4, not 2/3 and 2/3. -/
theorem reference_ladder_ratios :
    lowQ (120 : ℝ) 240 = (1 / 2 : ℝ) ∧
      highQ 240 360 = (3 / 4 : ℝ) := by
  constructor <;> norm_num [lowQ, highQ]

theorem reference_ladder_not_exact_koide :
    lowQ (120 : ℝ) 240 ≠ (2 / 3 : ℝ) := by
  norm_num [lowQ]

/-- C+a*cos(k*log f)+b*sin(k*log f) in the log-frequency variable ell. -/
noncomputable def logCosModel (ell k C a b : ℝ) : ℝ :=
  C + a * Real.cos (k * ell) + b * Real.sin (k * ell)

/-- Shifting every ordinate by d and the fitted intercept by d leaves
each least-squares residual unchanged. -/
theorem logCosResidual_offset (ell y k C a b d : ℝ) :
    (y + d) - logCosModel ell k (C + d) a b =
      y - logCosModel ell k C a b := by
  unfold logCosModel
  ring

noncomputable def logCosRSS {N : ℕ} (ell y : Fin N → ℝ)
    (k C a b : ℝ) : ℝ :=
  ∑ i : Fin N, (y i - logCosModel (ell i) k C a b) ^ 2

/-- The residual sum of squares is invariant under a common ordinate
offset when its free intercept shifts by the same amount. -/
theorem logCosRSS_offset {N : ℕ} (ell y : Fin N → ℝ)
    (k C a b d : ℝ) :
    logCosRSS ell (fun i => y i + d) k (C + d) a b =
      logCosRSS ell y k C a b := by
  simp [logCosRSS, logCosResidual_offset]

end PhotodiodeProtocol.HarmonicRatios

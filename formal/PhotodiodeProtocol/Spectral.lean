import Mathlib

/-!
# Exact finite-data spectral score

This is an exact real-arithmetic model of the per-capture spectral pipeline in
`scan_npy_logcos.compute_fft_spectrum`:

* Hann window `w[n] = ½ - ½ cos(2πn/(N-1))` (`numpy.hanning`; `w = 1` when `N = 1`),
* one-sided DFT `Y[k] = Σₙ w[n] y[n] e^{-2πikn/N}` (`numpy.fft.rfft`), bins `k = 0..⌊N/2⌋`,
* amplitude `|Y[k]|` (`numpy.abs`), and bin frequency `k·fs/N` (`numpy.fft.rfftfreq`),
* band selection `fmin ≤ f ≤ fmax` (`select_freq_window`).

On top of that sits a **declared** score: the band amplitude sum normalized by a reference-band
amplitude sum. The per-device persistence score is the late-minus-baseline difference of
normalized scores.

Proved here:
* `amplitude_smul`: `|Y_{cy}[k]| = |c| |Y_y[k]|`, so amplitudes scale exactly with detector gain.
* `normScore_gain_invariant`, `persistence_gain_invariant`: the normalized score and the
  persistence score are invariant under an independent nonzero gain on *each* capture.
* `bins_subset_rfft`: selected bins are valid one-sided bins.

This is a mathematical model of the declared transform. That the Python/NumPy implementation
agrees with it to floating-point accuracy is a separate conformance obligation, not proved here.
-/

namespace PhotodiodeProtocol.Spectral

open Finset Real

/-- `numpy.hanning(N)[n]`. -/
noncomputable def hann (N n : ℕ) : ℝ :=
  if N ≤ 1 then 1 else 1 / 2 - (1 / 2) * Real.cos (2 * π * n / (N - 1))

/-- Real part of `rfft` of the windowed signal at bin `k`. -/
noncomputable def dftRe {N : ℕ} (y : Fin N → ℝ) (k : ℕ) : ℝ :=
  ∑ n : Fin N, hann N n * y n * Real.cos (2 * π * k * n / N)

/-- Imaginary part (sign convention `e^{-iθ}`). -/
noncomputable def dftIm {N : ℕ} (y : Fin N → ℝ) (k : ℕ) : ℝ :=
  -∑ n : Fin N, hann N n * y n * Real.sin (2 * π * k * n / N)

/-- `numpy.abs(rfft(...))[k]`. -/
noncomputable def amplitude {N : ℕ} (y : Fin N → ℝ) (k : ℕ) : ℝ :=
  Real.sqrt (dftRe y k ^ 2 + dftIm y k ^ 2)

/-- `numpy.fft.rfftfreq(N, 1/fs)[k]`. -/
noncomputable def binFreq (N : ℕ) (fs : ℝ) (k : ℕ) : ℝ := k * fs / N

/-- Bins inside `[fmin, fmax]`. -/
noncomputable def bandBins (N : ℕ) (fs fmin fmax : ℝ) : Finset ℕ :=
  (range (N / 2 + 1)).filter fun k => fmin ≤ binFreq N fs k ∧ binFreq N fs k ≤ fmax

theorem bins_subset_rfft (N : ℕ) (fs fmin fmax : ℝ) :
    bandBins N fs fmin fmax ⊆ range (N / 2 + 1) := filter_subset _ _

theorem amplitude_nonneg {N : ℕ} (y : Fin N → ℝ) (k : ℕ) : 0 ≤ amplitude y k :=
  Real.sqrt_nonneg _

theorem dftRe_smul {N : ℕ} (c : ℝ) (y : Fin N → ℝ) (k : ℕ) :
    dftRe (c • y) k = c * dftRe y k := by
  unfold dftRe
  rw [mul_sum]
  refine sum_congr rfl fun n _ => ?_
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

theorem dftIm_smul {N : ℕ} (c : ℝ) (y : Fin N → ℝ) (k : ℕ) :
    dftIm (c • y) k = c * dftIm y k := by
  unfold dftIm
  rw [mul_neg, mul_sum]
  congr 1
  refine sum_congr rfl fun n _ => ?_
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

/-- Exact gain law for amplitudes. -/
theorem amplitude_smul {N : ℕ} (c : ℝ) (y : Fin N → ℝ) (k : ℕ) :
    amplitude (c • y) k = |c| * amplitude y k := by
  unfold amplitude
  rw [dftRe_smul, dftIm_smul]
  have : (c * dftRe y k) ^ 2 + (c * dftIm y k) ^ 2 = c ^ 2 * (dftRe y k ^ 2 + dftIm y k ^ 2) := by
    ring
  rw [this, Real.sqrt_mul (sq_nonneg c), Real.sqrt_sq_eq_abs]

/-- Band amplitude sum. -/
noncomputable def bandScore {N : ℕ} (y : Fin N → ℝ) (band : Finset ℕ) : ℝ :=
  ∑ k ∈ band, amplitude y k

theorem bandScore_nonneg {N : ℕ} (y : Fin N → ℝ) (band : Finset ℕ) : 0 ≤ bandScore y band :=
  sum_nonneg fun k _ => amplitude_nonneg y k

theorem bandScore_smul {N : ℕ} (c : ℝ) (y : Fin N → ℝ) (band : Finset ℕ) :
    bandScore (c • y) band = |c| * bandScore y band := by
  unfold bandScore
  rw [mul_sum]
  exact sum_congr rfl fun k _ => amplitude_smul c y k

/-- Normalized spectral score (target band / reference band). Division by zero gives `0`. -/
noncomputable def normScore {N : ℕ} (y : Fin N → ℝ) (band ref : Finset ℕ) : ℝ :=
  bandScore y band / bandScore y ref

theorem normScore_gain_invariant {N : ℕ} (c : ℝ) (hc : c ≠ 0) (y : Fin N → ℝ)
    (band ref : Finset ℕ) : normScore (c • y) band ref = normScore y band ref := by
  unfold normScore
  rw [bandScore_smul, bandScore_smul]
  have : |c| ≠ 0 := abs_ne_zero.mpr hc
  rw [mul_div_mul_left _ _ this]

/-- Primary per-device persistence score: late minus baseline normalized score. -/
noncomputable def persistence {N : ℕ} (baseline late : Fin N → ℝ) (band ref : Finset ℕ) : ℝ :=
  normScore late band ref - normScore baseline band ref

theorem persistence_gain_invariant {N : ℕ} (c₁ c₂ : ℝ) (h₁ : c₁ ≠ 0) (h₂ : c₂ ≠ 0)
    (baseline late : Fin N → ℝ) (band ref : Finset ℕ) :
    persistence (c₁ • baseline) (c₂ • late) band ref = persistence baseline late band ref := by
  unfold persistence
  rw [normScore_gain_invariant c₁ h₁, normScore_gain_invariant c₂ h₂]

end PhotodiodeProtocol.Spectral

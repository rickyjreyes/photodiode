import Mathlib

/-!
# Exact finite periodic autocorrelation: zero lag and normalization

Source: src/wct_entropy/diagnostics.py centered_autocorrelation.

The finite-grid correlation is the periodic sum of x[i] * x[shift i].
This module proves exact time-domain zero-lag and normalization identities.
DFT/Wiener–Khinchin correspondence and fftshift indexing are separate
remaining obligations, not silently asserted here.
-/
namespace EntropyProofs.DiscreteAutocorrelation

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

noncomputable def correlation (x : ι → ℝ) (shift : Equiv.Perm ι) : ℝ :=
  ∑ i : ι, x i * x (shift i)

noncomputable def zeroLag (x : ι → ℝ) : ℝ :=
  ∑ i : ι, x i ^ 2

theorem correlation_identity (x : ι → ℝ) :
    correlation x (Equiv.refl ι) = zeroLag x := by
  simp [correlation, zeroLag, pow_two]

theorem zeroLag_nonneg (x : ι → ℝ) : 0 ≤ zeroLag x := by
  unfold zeroLag
  exact Finset.sum_nonneg (fun i hi => sq_nonneg (x i))

noncomputable def normalizedCorrelation (x : ι → ℝ)
    (shift : Equiv.Perm ι) : ℝ :=
  correlation x shift / zeroLag x

theorem normalizedCorrelation_identity (x : ι → ℝ)
    (hzero : zeroLag x ≠ 0) :
    normalizedCorrelation x (Equiv.refl ι) = 1 := by
  simp [normalizedCorrelation, correlation_identity, div_self hzero]

theorem zero_field_zero_lag :
    zeroLag (fun _ : ι => (0 : ℝ)) = 0 := by
  simp [zeroLag]

theorem zeroLag_scale (x : ι → ℝ) (c : ℝ) :
    zeroLag (fun i => c * x i) = c ^ 2 * zeroLag x := by
  unfold zeroLag
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i hi
  ring

theorem correlation_scale (x : ι → ℝ) (c : ℝ)
    (shift : Equiv.Perm ι) :
    correlation (fun i => c * x i) shift = c ^ 2 * correlation x shift := by
  unfold correlation
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i hi
  ring

theorem normalizedCorrelation_scale (x : ι → ℝ) (c : ℝ)
    (shift : Equiv.Perm ι) (hc : c ≠ 0) :
    normalizedCorrelation (fun i => c * x i) shift =
      normalizedCorrelation x shift := by
  unfold normalizedCorrelation
  rw [correlation_scale, zeroLag_scale]
  exact mul_div_mul_left _ _ (pow_ne_zero _ hc)

end EntropyProofs.DiscreteAutocorrelation

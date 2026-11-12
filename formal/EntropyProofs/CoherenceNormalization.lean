import EntropyProofs.DiscreteAutocorrelationParseval

/-!
# fftshift and normalized correlation at the center of a 2×2 real grid

The zero-lag value of NumPy's ifft2(|fft2(field)|²) is at (0,0).
For an even 2×2 grid, fftshift moves that value to (1,1).
The source divides by that value, but rejects zero intensity.

No statement here identifies this diagnostic with a physical
coherence length or with an arbitrary-size-grid theorem.
-/

namespace EntropyProofs.CoherenceNormalization

open EntropyProofs.DiscreteAutocorrelationParseval

/-- The index rearrangement performed by np.fft.fftshift for the 2×2 case. -/
noncomputable def shifted (a b c d : ℝ) (i j : Fin 2) : ℝ :=
  if i = 0 then
    if j = 0 then invPower11 a b c d else invPower10 a b c d
  else
    if j = 0 then invPower01 a b c d else invPower00 a b c d

theorem shifted_center_is_energy (a b c d : ℝ) :
    shifted a b c d 1 1 = energy a b c d := by
  simp [shifted, wiener_khinchin_zero_lag]

theorem shifted_origin_is_diagonal_lag (a b c d : ℝ) :
    shifted a b c d 0 0 = 2 * (a * d + b * c) := by
  simp [shifted, wiener_khinchin_diagonal_lag]

/-- Source's zero-lag normalizer, with no clipping or regularization. -/
noncomputable def normalizedAtCenter (a b c d : ℝ) : ℝ :=
  shifted a b c d 1 1 / energy a b c d

/-- The normalized fftshift center is exactly one for a nonzero field. -/
theorem normalized_center_one (a b c d : ℝ)
    (hnz : energy a b c d ≠ 0) :
    normalizedAtCenter a b c d = 1 := by
  rw [normalizedAtCenter, shifted_center_is_energy]
  exact div_self hnz

/-- Explicit zero-field pathology: the source deliberately raises rather
than accepting Lean's total-real convention 0/0=0. -/
theorem normalized_zero_field : normalizedAtCenter 0 0 0 0 = 0 := by
  simp [normalizedAtCenter, shifted, invPower00,
    spectral00, spectral01, spectral10, spectral11,
    fft00, fft01, fft10, fft11, energy]

/-- Signed correlation at the column lag is bounded above by zero-lag
intensity. This is a scoped finite 2×2 Cauchy–Schwarz consequence. -/
theorem column_lag_le_energy (a b c d : ℝ) :
    invPower01 a b c d ≤ energy a b c d := by
  rw [wiener_khinchin_column_lag]
  dsimp [energy]
  nlinarith [sq_nonneg (a - b), sq_nonneg (c - d)]

theorem neg_energy_le_column_lag (a b c d : ℝ) :
    -energy a b c d ≤ invPower01 a b c d := by
  rw [wiener_khinchin_column_lag]
  dsimp [energy]
  nlinarith [sq_nonneg (a + b), sq_nonneg (c + d)]

end EntropyProofs.CoherenceNormalization

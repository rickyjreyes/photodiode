import Mathlib

/-! R18 (DOI 10.5281/zenodo.19705254).
A universal amplitude bound for a log-periodic C9 ansatz.
It does NOT prove a WCT effect or a fitted LHCb anomaly. -/
namespace ReleaseProofs.R18

noncomputable def logPeriodicC9Correction
    (amplitude frequency phase q2 reference : ℝ) : ℝ :=
  amplitude * Real.cos (frequency * Real.log (q2 / reference) + phase)

theorem logPeriodicC9Correction_bounded
    (amplitude frequency phase q2 reference : ℝ) :
    |logPeriodicC9Correction amplitude frequency phase q2 reference| ≤
      |amplitude| := by
  unfold logPeriodicC9Correction
  rw [abs_mul]
  calc
    |amplitude| * |Real.cos (frequency * Real.log (q2 / reference) + phase)| ≤
        |amplitude| * 1 :=
      mul_le_mul_of_nonneg_left (Real.abs_cos_le_one _) (abs_nonneg amplitude)
    _ = |amplitude| := by ring

end ReleaseProofs.R18

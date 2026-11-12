import ReleaseProofs.R18_C9

/-!
# R18: domain and phase identities for the declared C9 log-periodic ansatz.

The q² and reference inputs are strictly positive in the physical domain.
The underlying logarithm is totalized by Lean outside that domain, but no
physics is asserted there. These statements are about the declared ansatz
and do not prove an LHCb fit or model-identification claim.
-/

namespace ReleaseProofs.R18

/-- C9 baseline plus log-periodic correction. -/
noncomputable def logPeriodicC9
    (baseline amplitude frequency phase q2 reference : ℝ) : ℝ :=
  baseline + logPeriodicC9Correction amplitude frequency phase q2 reference

/-- Exact cancellation of the baseline in the correction. -/
theorem logPeriodicC9_minus_baseline
    (baseline amplitude frequency phase q2 reference : ℝ) :
    logPeriodicC9 baseline amplitude frequency phase q2 reference - baseline =
      logPeriodicC9Correction amplitude frequency phase q2 reference := by
  simp [logPeriodicC9]

/-- The correction vanishes identically when the amplitude is zero. -/
theorem logPeriodicC9_zero_amplitude
    (baseline frequency phase q2 reference : ℝ) :
    logPeriodicC9 baseline 0 frequency phase q2 reference = baseline := by
  simp [logPeriodicC9, logPeriodicC9Correction]

/-- Model-independent amplitude enclosure, without assuming the sign of A. -/
theorem logPeriodicC9_enclosed
    (baseline amplitude frequency phase q2 reference : ℝ) :
    baseline - |amplitude| ≤
      logPeriodicC9 baseline amplitude frequency phase q2 reference ∧
    logPeriodicC9 baseline amplitude frequency phase q2 reference ≤
      baseline + |amplitude| := by
  have bound := logPeriodicC9Correction_bounded amplitude frequency phase q2 reference
  rw [abs_le] at bound
  constructor <;> dsimp [logPeriodicC9] <;> linarith

/-- Phase shifts by a full turn leave the ansatz invariant. -/
theorem logPeriodicC9_phase_period
    (baseline amplitude frequency phase q2 reference : ℝ) :
    logPeriodicC9 baseline amplitude frequency (phase + 2 * Real.pi) q2 reference =
      logPeriodicC9 baseline amplitude frequency phase q2 reference := by
  unfold logPeriodicC9 logPeriodicC9Correction
  have h :
      frequency * Real.log (q2 / reference) + (phase + 2 * Real.pi) =
        (frequency * Real.log (q2 / reference) + phase) + 2 * Real.pi := by ring
  rw [h, Real.cos_add_two_pi]

/-- Common change of units for positive q² and reference scale does
not modify the ratio (or the log-periodic correction). -/
theorem logPeriodicC9_common_units
    (baseline amplitude frequency phase q2 reference unit : ℝ)
    (hunit : unit ≠ 0) (href : reference ≠ 0) :
    logPeriodicC9 baseline amplitude frequency phase
        (q2 * unit) (reference * unit) =
      logPeriodicC9 baseline amplitude frequency phase q2 reference := by
  have hratio : (q2 * unit) / (reference * unit) = q2 / reference := by
    field_simp
    <;> ring
  simp [logPeriodicC9, logPeriodicC9Correction, hratio]

/-- At zero frequency the ansatz has no q² dependence. -/
theorem logPeriodicC9_zero_frequency
    (baseline amplitude phase q2 ref q2' ref' : ℝ) :
    logPeriodicC9 baseline amplitude 0 phase q2 ref =
      logPeriodicC9 baseline amplitude 0 phase q2' ref' := by
  simp [logPeriodicC9, logPeriodicC9Correction]

end ReleaseProofs.R18

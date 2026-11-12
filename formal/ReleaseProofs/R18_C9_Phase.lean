import ReleaseProofs.R18_C9_Domain

/-!
# Release 18: logarithmic phase, reference scaling, and exact period

The ansatz is C9 = baseline + A*cos(ω*log(q²/q0²) + phase).
The q², reference and positive-frequency hypotheses are stated explicitly.
These are model-algebra results, NOT experimental inferences.
-/

namespace ReleaseProofs.R18

/-- Compensating for a change of reference by a phase shift leaves the
log-periodic argument unchanged on the positive domain. -/
theorem logPeriodicC9_reference_phase_shift
    (baseline amplitude frequency phase q2 reference newReference : ℝ)
    (hq : 0 < q2) (href : 0 < reference) (hnew : 0 < newReference) :
    logPeriodicC9 baseline amplitude frequency
        (phase + frequency * Real.log (newReference / reference))
        q2 newReference =
      logPeriodicC9 baseline amplitude frequency phase q2 reference := by
  unfold logPeriodicC9 logPeriodicC9Correction
  have hr1 : q2 / newReference > 0 := div_pos hq hnew
  have hr2 : newReference / reference > 0 := div_pos hnew href
  have hlog : Real.log (q2 / newReference) +
      Real.log (newReference / reference) = Real.log (q2 / reference) := by
    rw [← Real.log_mul (ne_of_gt hr1) (ne_of_gt hr2)]
    congr 1
    field_simp
  have harg :
      frequency * Real.log (q2 / newReference) +
          (phase + frequency * Real.log (newReference / reference)) =
        frequency * Real.log (q2 / reference) + phase := by
    linear_combination frequency * hlog
  rw [harg]

/-- On the positive q² domain, scaling q² by exp(2π/ω) advances the
oscillatory phase by precisely 2π, assuming nonzero ω. -/
theorem logPeriodicC9_log_period
    (baseline amplitude frequency phase q2 reference : ℝ)
    (hq : 0 < q2) (href : 0 < reference) (hf : frequency ≠ 0) :
    logPeriodicC9 baseline amplitude frequency phase
        (q2 * Real.exp (2 * Real.pi / frequency)) reference =
      logPeriodicC9 baseline amplitude frequency phase q2 reference := by
  unfold logPeriodicC9 logPeriodicC9Correction
  have hexp : 0 < Real.exp (2 * Real.pi / frequency) := Real.exp_pos _
  have hratio : (q2 * Real.exp (2 * Real.pi / frequency)) / reference =
      (q2 / reference) * Real.exp (2 * Real.pi / frequency) := by ring
  rw [hratio]
  have hlog : Real.log ((q2 / reference) *
        Real.exp (2 * Real.pi / frequency)) =
      Real.log (q2 / reference) + 2 * Real.pi / frequency := by
    rw [Real.log_mul (ne_of_gt (div_pos hq href)) (ne_of_gt hexp),
        Real.log_exp]
  rw [hlog]
  have harg :
      frequency * (Real.log (q2 / reference) + 2 * Real.pi / frequency) + phase =
        (frequency * Real.log (q2 / reference) + phase) + 2 * Real.pi := by
    field_simp [hf]
    ring
  rw [harg, Real.cos_add_two_pi]

end ReleaseProofs.R18

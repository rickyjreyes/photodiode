import Mathlib

/-!
# Exact pointwise phase mass conservation for the CQ Strang splitter

Source: src/wct_pde/projection_free.py, SpectralSplitStepper.step.
The nonlinear step multiplies the field by exp(i*theta(|psi|²)).
On real/imaginary components this is a 2D rotation.

This module proves arbitrary sitewise unit-circle rotations preserve the
discrete weighted L² mass *exactly*, even when the angle depends on the
field amplitude. The real Fourier linear half-step is mass preserving
under a separate Parseval/unitary FFT theorem, not proved here.
Floating-point roundoff and Strang convergence are separate questions.
-/

namespace ProjectionFree.PhaseMassConservation

def density (re im : ℝ) : ℝ := re ^ 2 + im ^ 2

def rotateRe (c s re im : ℝ) : ℝ := c * re - s * im

def rotateIm (c s re im : ℝ) : ℝ := s * re + c * im

/-- A pointwise unit-circle rotation leaves squared amplitude invariant. -/
theorem rotation_preserves_density (c s re im : ℝ)
    (hunit : c ^ 2 + s ^ 2 = 1) :
    density (rotateRe c s re im) (rotateIm c s re im) =
      density re im := by
  unfold density rotateRe rotateIm
  calc
    (c * re - s * im) ^ 2 + (s * re + c * im) ^ 2 =
        (c ^ 2 + s ^ 2) * (re ^ 2 + im ^ 2) := by ring
    _ = re ^ 2 + im ^ 2 := by rw [hunit]; ring

/-- Arbitrary nonuniform quadrature weights and finite grid indices. -/
def gridMass {n : ℕ} (weights re im : Fin n → ℝ) : ℝ :=
  ∑ i, weights i * density (re i) (im i)

/-- All local phase shifts, including amplitude-dependent phase shifts,
preserve the exact discrete mass when each phase has unit modulus. -/
theorem sitewise_phase_preserves_mass {n : ℕ}
    (weights re im c s : Fin n → ℝ)
    (hunit : ∀ i, (c i) ^ 2 + (s i) ^ 2 = 1) :
    gridMass weights
        (fun i => rotateRe (c i) (s i) (re i) (im i))
        (fun i => rotateIm (c i) (s i) (re i) (im i)) =
      gridMass weights re im := by
  unfold gridMass
  apply Finset.sum_congr rfl
  intro i hi
  rw [rotation_preserves_density (c i) (s i) (re i) (im i) (hunit i)]

/-- The trigonometric nonlinear phase is a concrete unit-modulus rotation.
Identification with the source's complex exponential is Euler's formula. -/
theorem sine_cosine_phase_preserves_mass {n : ℕ}
    (weights re im phase : Fin n → ℝ) :
    gridMass weights
        (fun i => rotateRe (Real.cos (phase i)) (Real.sin (phase i))
            (re i) (im i))
        (fun i => rotateIm (Real.cos (phase i)) (Real.sin (phase i))
            (re i) (im i)) =
      gridMass weights re im := by
  apply sitewise_phase_preserves_mass
  intro i
  nlinarith [Real.sin_sq_add_cos_sq (phase i)]

end ProjectionFree.PhaseMassConservation

import Mathlib

/-!
# The rigid-translation finite-difference obstruction

Source: wct_topology_flow.py, functions curvature, self_avoidance_energy,
total_energy, and energy_gradient.

The model follows the source's closed-chain neighbors (supplied as finite index
maps), cross-product norm, edge norm, regularized denominator, ordered-pair
self-avoidance mask, and squared-energy sum.

All statements are in *exact real arithmetic*. The implementation uses CuPy
floating-point operations and is tested separately.
-/

namespace ProteinProofs.RigidTranslationNullDerivative

abbrev Point3 := Fin 3 → ℝ
abbrev Configuration (n : ℕ) := Fin n → Point3

/-- Translating every vertex by the same displacement. -/
def translate {n : ℕ} (p : Configuration n) (t : Point3) : Configuration n :=
  fun i a => p i a + t a

/-- A directed edge, also used to construct pairwise distances. -/
def difference {n : ℕ} (p : Configuration n) (i j : Fin n) : Point3 :=
  fun a => p i a - p j a

/-- The translation drops out of *every* pairwise difference. -/
theorem difference_translate {n : ℕ}
    (p : Configuration n) (t : Point3) (i j : Fin n) :
    difference (translate p t) i j = difference p i j := by
  funext a
  simp only [difference, translate]
  ring

noncomputable def norm3 (v : Point3) : ℝ :=
  Real.sqrt (v 0 ^ 2 + v 1 ^ 2 + v 2 ^ 2)

/-- Coordinate formula matching the three-dimensional cross product. -/
def cross3 (u v : Point3) : Point3 :=
  fun a =>
    if a = 0 then u 1 * v 2 - u 2 * v 1
    else if a = 1 then u 2 * v 0 - u 0 * v 2
    else u 0 * v 1 - u 1 * v 0

/-- CuPy curvature proxy at a vertex. Prev/next are cyclic index maps in Python.
The one-trillionth regularizer is the exact-real counterpart of 1e-12. -/
noncomputable def curvature {n : ℕ} (p : Configuration n)
    (prev next : Fin n → Fin n) (i : Fin n) : ℝ :=
  let v1 := difference p i (prev i)
  let v2 := difference p (next i) i
  norm3 (cross3 v1 v2) / (norm3 v1 ^ 3 + (1 / 1000000000000 : ℝ))

noncomputable def curvatureEnergy {n : ℕ} (p : Configuration n)
    (prev next : Fin n → Fin n) : ℝ :=
  ∑ i : Fin n, curvature p prev next i ^ 2

/-- CuPy's ordered-pair distance mask: dist < eps and dist > 1e-6.
In particular diagonal pairs are excluded and pairs are counted twice. -/
noncomputable def selfAvoidTerm {n : ℕ} (p : Configuration n)
    (eps : ℝ) (i j : Fin n) : ℝ :=
  let dist := norm3 (difference p i j)
  if dist < eps ∧ dist > (1 / 1000000 : ℝ) then (eps - dist) ^ 2 else 0

noncomputable def selfAvoidEnergy {n : ℕ}
    (p : Configuration n) (eps : ℝ) : ℝ :=
  ∑ i : Fin n, ∑ j : Fin n, selfAvoidTerm p eps i j

noncomputable def totalEnergy {n : ℕ} (p : Configuration n)
    (prev next : Fin n → Fin n) (eps weight : ℝ) : ℝ :=
  curvatureEnergy p prev next + weight * selfAvoidEnergy p eps

theorem curvature_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (prev next : Fin n → Fin n) (i : Fin n) :
    curvature (translate p t) prev next i = curvature p prev next i := by
  simp [curvature, difference_translate]

theorem curvatureEnergy_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (prev next : Fin n → Fin n) :
    curvatureEnergy (translate p t) prev next = curvatureEnergy p prev next := by
  simp [curvatureEnergy, curvature_translate]

theorem selfAvoidTerm_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (eps : ℝ) (i j : Fin n) :
    selfAvoidTerm (translate p t) eps i j = selfAvoidTerm p eps i j := by
  simp [selfAvoidTerm, difference_translate]

theorem selfAvoidEnergy_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (eps : ℝ) :
    selfAvoidEnergy (translate p t) eps = selfAvoidEnergy p eps := by
  simp [selfAvoidEnergy, selfAvoidTerm_translate]

/-- Exact translation invariance of both terms in the implemented energy. -/
theorem totalEnergy_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (prev next : Fin n → Fin n) (eps weight : ℝ) :
    totalEnergy (translate p t) prev next eps weight =
      totalEnergy p prev next eps weight := by
  simp [totalEnergy, curvatureEnergy_translate, selfAvoidEnergy_translate]

def axisShift (axis : Fin 3) (h : ℝ) : Point3 :=
  fun a => if a = axis then h else 0

/-- The quantity written to *every* vertex in one column of energy_gradient.
Unlike a vertexwise variation, this perturbs the entire configuration. -/
noncomputable def centeredBulkDifference {n : ℕ} (p : Configuration n)
    (prev next : Fin n → Fin n) (eps weight : ℝ)
    (axis : Fin 3) (h : ℝ) : ℝ :=
  (totalEnergy (translate p (axisShift axis h)) prev next eps weight -
     totalEnergy (translate p (axisShift axis (-h))) prev next eps weight) / (2 * h)

/-- Even the centered quotient with h=0 has value zero in Lean's total
division convention; the physically meaningful finite-difference formula
is only intended for nonzero h. -/
theorem centeredBulkDifference_zero {n : ℕ} (p : Configuration n)
    (prev next : Fin n → Fin n) (eps weight : ℝ) (axis : Fin 3) (h : ℝ) :
    centeredBulkDifference p prev next eps weight axis h = 0 := by
  simp [centeredBulkDifference, totalEnergy_translate]

end ProteinProofs.RigidTranslationNullDerivative

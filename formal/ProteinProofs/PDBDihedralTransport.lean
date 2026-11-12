import ProteinProofs.DimensionalTransportCounterexample

/-!
# PDB pseudo-dihedral transport: exact source-coordinate obligations

Source: wct_pdb_transport.py::dihedral_angles.

The algorithm forms three consecutive edge vectors; regularizes and
normalizes them; forms two cross-product normals; normalizes those normals;
then passes two scalar invariants into atan2. We formalize the exact pair
of inputs and translation invariance, without identifying an open-chain
pseudo-dihedral phase sum with holonomy.

The existence of an atan2 value at a degenerate segment does not produce
a geometric Frenet frame there.
-/

namespace ProteinProofs.PDBDihedralTransport

open ProteinProofs.RigidTranslationNullDerivative
open ProteinProofs.DimensionalTransportCounterexample

/-- Source's regularized normalization, δ=10⁻⁸. -/
noncomputable def regularizeUnit (v : Point3) : Point3 :=
  fun a => v a / (norm3 v + (1 / 100000000 : ℝ))

/-- The raw edge normals are defined geometrically only when both
consecutive cross products are nonzero. -/
def validNormals {n : ℕ} (p : Configuration n) (i j k l : Fin n) : Prop :=
  cross3 (difference p j i) (difference p k j) ≠ 0 ∧
  cross3 (difference p k j) (difference p l k) ≠ 0

/-- Source's two real scalar arguments to numpy.arctan2(y, x). -/
noncomputable def dihedralXY {n : ℕ}
    (p : Configuration n) (i j k l : Fin n) : ℝ × ℝ :=
  let b1 := regularizeUnit (difference p j i)
  let b2 := regularizeUnit (difference p k j)
  let b3 := regularizeUnit (difference p l k)
  let n1 := regularizeUnit (cross3 b1 b2)
  let n2 := regularizeUnit (cross3 b2 b3)
  (inner3 n1 n2, inner3 (cross3 n1 b2) n2)

/-- Translation leaves each argument supplied to atan2 identical. -/
theorem dihedralXY_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (i j k l : Fin n) :
    dihedralXY (translate p t) i j k l = dihedralXY p i j k l := by
  simp [dihedralXY, difference_translate]

/-- This is valid for any deterministic atan2 or phase convention. -/
theorem dihedralPhase_translate {n : ℕ} (p : Configuration n)
    (t : Point3) (i j k l : Fin n) (phase : ℝ → ℝ → ℝ) :
    phase (dihedralXY (translate p t) i j k l).2
          (dihedralXY (translate p t) i j k l).1 =
      phase (dihedralXY p i j k l).2 (dihedralXY p i j k l).1 := by
  rw [dihedralXY_translate]

/-- A cross product of parallel equal vectors is zero. -/
theorem cross3_self (v : Point3) : cross3 v v = 0 := by
  funext a
  fin_cases a <;> simp [cross3] <;> ring

/-- For a constant chain, the source's finite regularization gives
a value, but the nondegenerate geometric-normal conditions fail. -/
theorem constant_chain_invalid {n : ℕ} (v : Point3) (i j k l : Fin n) :
    ¬ validNormals (fun _ => v) i j k l := by
  simp [validNormals, difference, cross3]

/-- Circular angular distance applied to a remainder r in [0,period). -/
def circularDistance (period r : ℝ) : ℝ :=
  min r (period - r)

/-- A symmetric circular-distance test treats δ and period−δ equally. -/
theorem circularDistance_symmetric (period delta : ℝ) :
    circularDistance period delta =
      circularDistance period (period - delta) := by
  simp [circularDistance, min_comm]

/-- For δ < tolerance < period−δ, the *remainders* δ and
period−δ lead to opposite one-sided decisions. Correspondence of these
remainders to angle modulo values is exercised against numpy below. -/
theorem one_sided_threshold_asymmetry (period delta tolerance : ℝ)
    (hsmall : delta < tolerance) (hlarge : tolerance < period - delta) :
    delta < tolerance ∧ ¬ (period - delta < tolerance) := by
  constructor
  · exact hsmall
  · linarith

end ProteinProofs.PDBDihedralTransport

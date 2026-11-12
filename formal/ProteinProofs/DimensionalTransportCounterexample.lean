import ProteinProofs.RigidTranslationNullDerivative

/-!
# Dimension embedding obstruction to the unsigned-turning closure statistic

Source: wct_dimensional_locking.py:transport_angle.

The implemented observable is the sum of arccosines of consecutive
regularized normalized tangent dot products, clipped into [-1, 1].
Appending a zero coordinate preserves each operation exactly.

This result does *not* identify tangent-angle accumulation with an ordered
frame-holonomy operator, and does not infer a physical preference for any
dimension from stochastic sampling.
-/

namespace ProteinProofs.DimensionalTransportCounterexample

open ProteinProofs.RigidTranslationNullDerivative

abbrev Point4 := Fin 4 → ℝ

/-- The standard isometric coordinate inclusion R³ ↪ R⁴. -/
def embed3 (v : Point3) : Point4 :=
  fun a =>
    if a = 0 then v 0
    else if a = 1 then v 1
    else if a = 2 then v 2
    else 0

noncomputable def norm4 (v : Point4) : ℝ :=
  Real.sqrt (v 0 ^ 2 + v 1 ^ 2 + v 2 ^ 2 + v 3 ^ 2)

def inner3 (u v : Point3) : ℝ :=
  u 0 * v 0 + u 1 * v 1 + u 2 * v 2

def inner4 (u v : Point4) : ℝ :=
  u 0 * v 0 + u 1 * v 1 + u 2 * v 2 + u 3 * v 3

theorem norm4_embed3 (v : Point3) :
    norm4 (embed3 v) = norm3 v := by
  simp [norm4, norm3, embed3]

theorem inner4_embed3 (u v : Point3) :
    inner4 (embed3 u) (embed3 v) = inner3 u v := by
  simp [inner4, inner3, embed3]

/-- Python divides a tangent by its norm plus 1e-8. -/
noncomputable def normalize3 (v : Point3) : Point3 :=
  fun a => v a / (norm3 v + (1 / 100000000 : ℝ))

noncomputable def normalize4 (v : Point4) : Point4 :=
  fun a => v a / (norm4 v + (1 / 100000000 : ℝ))

theorem normalize4_embed3 (v : Point3) :
    normalize4 (embed3 v) = embed3 (normalize3 v) := by
  funext a
  fin_cases a <;> simp [normalize4, normalize3, embed3, norm4_embed3]

/-- Exact real-valued counterpart of numpy.clip(x, -1, 1). -/
def clipDot (x : ℝ) : ℝ := max (-1) (min 1 x)

noncomputable def turn3 (u v : Point3) : ℝ :=
  Real.arccos (clipDot (inner3 (normalize3 u) (normalize3 v)))

noncomputable def turn4 (u v : Point4) : ℝ :=
  Real.arccos (clipDot (inner4 (normalize4 u) (normalize4 v)))

theorem turn_embed3 (u v : Point3) :
    turn4 (embed3 u) (embed3 v) = turn3 u v := by
  simp [turn4, turn3, normalize4_embed3, inner4_embed3]

/-- The actual unsigned turning sum over a periodic list of edges. -/
noncomputable def turningScore3 {n : ℕ}
    (edges : Fin n → Point3) (next : Fin n → Fin n) : ℝ :=
  ∑ i : Fin n, turn3 (edges i) (edges (next i))

noncomputable def turningScore4 {n : ℕ}
    (edges : Fin n → Point4) (next : Fin n → Fin n) : ℝ :=
  ∑ i : Fin n, turn4 (edges i) (edges (next i))

/-- A 3D tangent sequence and its 4D embedding have *identical* scores,
even with the source's norm regularizer and clipping convention. -/
theorem turningScore_embedding {n : ℕ}
    (edges : Fin n → Point3) (next : Fin n → Fin n) :
    turningScore4 (fun i => embed3 (edges i)) next =
      turningScore3 edges next := by
  simp [turningScore4, turningScore3, turn_embed3]

end ProteinProofs.DimensionalTransportCounterexample

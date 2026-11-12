import Mathlib

/-! R14 (DOI 10.5281/zenodo.17887562): exact Koide geometry
and conditional spin-harmonic algebra, not a prediction of particle masses. -/
namespace ReleaseProofs.R14

noncomputable def koideFromRoots (x y z : ℝ) : ℝ :=
  (x ^ 2 + y ^ 2 + z ^ 2) / (x + y + z) ^ 2

noncomputable def angleCosSquared (x y z : ℝ) : ℝ :=
  (x + y + z) ^ 2 / (3 * (x ^ 2 + y ^ 2 + z ^ 2))

theorem koide_geometric_invariant (x y z : ℝ)
    (hs : x + y + z ≠ 0) (hn : x ^ 2 + y ^ 2 + z ^ 2 ≠ 0) :
    koideFromRoots x y z * angleCosSquared x y z = (1 : ℝ) / 3 := by
  unfold koideFromRoots angleCosSquared
  field_simp [hs, hn]
  <;> ring

noncomputable def harmonicK2 (s : ℝ) : ℝ := 2 * s / (2 * s + 1)
noncomputable def harmonicQ (s : ℝ) : ℝ := (2 * s + 1) / (6 * s)

theorem spin_harmonic_invariant (s : ℝ)
    (hs : s ≠ 0) (hd : 2 * s + 1 ≠ 0) :
    harmonicQ s * harmonicK2 s = (1 : ℝ) / 3 := by
  unfold harmonicQ harmonicK2
  field_simp [hs, hd]
  <;> ring

end ReleaseProofs.R14

import ReleaseProofs.R14_Koide

/-!
# R14: Nonnegative mass-vector domain and exact geometric identity

This lifts the pre-existing R14 root-vector algebra to the paper's
three nonnegative masses. It does not identify the spin-harmonic K(s)
with a measured lepton angle; that connection is a separate model premise.

See ReleaseProofs.R14_Koide for the previously formalized conditional
spin algebra.
-/
namespace ReleaseProofs.R14PhysicalMass

open ReleaseProofs.R14

/-- The mass-based Koide ratio, defined only as a total real function.
For zero total square-root mass the denominator is zero and Lean returns
zero; this is NOT a physical definition for the zero triplet. -/
noncomputable def koideMass (m1 m2 m3 : ℝ) : ℝ :=
  (m1 + m2 + m3) /
    (Real.sqrt m1 + Real.sqrt m2 + Real.sqrt m3) ^ 2

/-- The actual normalized mass-vector geometric factor cos²(theta). -/
noncomputable def massCosSquared (m1 m2 m3 : ℝ) : ℝ :=
  (Real.sqrt m1 + Real.sqrt m2 + Real.sqrt m3) ^ 2 /
    (3 * (m1 + m2 + m3))

/-- For nonnegative masses, the squared roots reproduce the sum of masses. -/
theorem root_square_sum (m1 m2 m3 : ℝ)
    (h1 : 0 ≤ m1) (h2 : 0 ≤ m2) (h3 : 0 ≤ m3) :
    (Real.sqrt m1) ^ 2 + (Real.sqrt m2) ^ 2 +
        (Real.sqrt m3) ^ 2 = m1 + m2 + m3 := by
  rw [Real.sq_sqrt h1, Real.sq_sqrt h2, Real.sq_sqrt h3]

/-- The physical mass ratio is exactly the established R14 root-vector ratio. -/
theorem koideMass_eq_roots (m1 m2 m3 : ℝ)
    (h1 : 0 ≤ m1) (h2 : 0 ≤ m2) (h3 : 0 ≤ m3) :
    koideMass m1 m2 m3 =
      koideFromRoots (Real.sqrt m1) (Real.sqrt m2) (Real.sqrt m3) := by
  unfold koideMass koideFromRoots
  rw [root_square_sum m1 m2 m3 h1 h2 h3]

/-- The mass-vector cos² factor is the root-vector one. -/
theorem massCosSquared_eq_roots (m1 m2 m3 : ℝ)
    (h1 : 0 ≤ m1) (h2 : 0 ≤ m2) (h3 : 0 ≤ m3) :
    massCosSquared m1 m2 m3 =
      angleCosSquared (Real.sqrt m1) (Real.sqrt m2) (Real.sqrt m3) := by
  unfold massCosSquared angleCosSquared
  rw [root_square_sum m1 m2 m3 h1 h2 h3]

/-- Exact identity Q cos²(theta)=1/3 for a nonzero admissible triplet. -/
theorem koideMass_geometry (m1 m2 m3 : ℝ)
    (h1 : 0 ≤ m1) (h2 : 0 ≤ m2) (h3 : 0 ≤ m3)
    (hroot : Real.sqrt m1 + Real.sqrt m2 + Real.sqrt m3 ≠ 0)
    (hmass : m1 + m2 + m3 ≠ 0) :
    koideMass m1 m2 m3 * massCosSquared m1 m2 m3 =
      (1 : ℝ) / 3 := by
  rw [koideMass_eq_roots m1 m2 m3 h1 h2 h3,
      massCosSquared_eq_roots m1 m2 m3 h1 h2 h3]
  apply koide_geometric_invariant
  · exact hroot
  · rw [root_square_sum m1 m2 m3 h1 h2 h3]
    exact hmass

/-- Equivalent inverse-cos² identity on the same nondegenerate domain. -/
theorem koideMass_inv_cosSquared (m1 m2 m3 : ℝ)
    (h1 : 0 ≤ m1) (h2 : 0 ≤ m2) (h3 : 0 ≤ m3)
    (hroot : Real.sqrt m1 + Real.sqrt m2 + Real.sqrt m3 ≠ 0)
    (hmass : m1 + m2 + m3 ≠ 0) :
    koideMass m1 m2 m3 =
      1 / (3 * massCosSquared m1 m2 m3) := by
  unfold koideMass massCosSquared
  field_simp [hroot, hmass]

/-- The all-zero mass vector is outside the meaningful angular domain. -/
theorem koideMass_zero : koideMass 0 0 0 = 0 := by
  simp [koideMass]

theorem massCosSquared_zero : massCosSquared 0 0 0 = 0 := by
  simp [massCosSquared]

/-- Both spin functions have positive values for positive spin parameter,
but their identification with any physical mass vector is not established. -/
theorem spin_quantities_positive (s : ℝ) (hs : 0 < s) :
    0 < harmonicK2 s ∧ 0 < harmonicQ s := by
  unfold harmonicK2 harmonicQ
  constructor
  · apply div_pos (by positivity) (by linarith)
  · apply div_pos (by linarith) (by positivity)

end ReleaseProofs.R14PhysicalMass

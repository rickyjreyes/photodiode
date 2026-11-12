import ProjectionFree.PhaseMassConservation

/-!
# Mass-conserving split maps do not necessarily conserve energy

An exact two-site structural counterexample to the implication
"mass-preserving substeps imply Hamiltonian-preserving full steps".

swap is a two-site Fourier-multiplier map with eigenvalues (1,-1)
in the DFT basis. phaseFirst rotates the first complex site by 90 degrees,
leaving the second unchanged. Their palindromic composition has exact
mass conservation but does not preserve the specified coupled energy.

This is a two-site representative in exact arithmetic, not proof
of energy drift for the full production nonlinear coefficients.
-/
namespace ProjectionFree.SplitEnergyDriftLimits

open ProjectionFree.PhaseMassConservation

structure TwoSite where
  ar : ℝ
  ai : ℝ
  br : ℝ
  bi : ℝ

def mass (z : TwoSite) : ℝ :=
  density z.ar z.ai + density z.br z.bi

def swap (z : TwoSite) : TwoSite :=
  ⟨z.br, z.bi, z.ar, z.ai⟩

def phaseFirst (z : TwoSite) : TwoSite :=
  ⟨-z.ai, z.ar, z.br, z.bi⟩

def split (z : TwoSite) : TwoSite :=
  swap (phaseFirst (swap z))

theorem swap_preserves_mass (z : TwoSite) :
    mass (swap z) = mass z := by
  cases z
  simp [mass, swap]
  ring

theorem phaseFirst_preserves_mass (z : TwoSite) :
    mass (phaseFirst z) = mass z := by
  cases z
  simp [mass, phaseFirst, density]
  ring

theorem split_preserves_mass (z : TwoSite) :
    mass (split z) = mass z := by
  unfold split
  rw [swap_preserves_mass, phaseFirst_preserves_mass, swap_preserves_mass]

/-- Coupled quadratic energy, distinct from pure L2 mass. -/
def coupledEnergy (z : TwoSite) : ℝ :=
  mass z + (z.ar + z.br) ^ 2 + (z.ai + z.bi) ^ 2

def witness : TwoSite := ⟨1, 0, 1, 0⟩

theorem witness_mass : mass witness = 2 := by
  norm_num [witness, mass, density]

theorem witness_split_mass : mass (split witness) = 2 := by
  rw [split_preserves_mass, witness_mass]

theorem witness_energy_before : coupledEnergy witness = 6 := by
  norm_num [coupledEnergy, mass, density, witness]

theorem witness_energy_after : coupledEnergy (split witness) = 4 := by
  norm_num [coupledEnergy, mass, density, split, swap, phaseFirst, witness]

theorem split_not_energy_preserving :
    coupledEnergy (split witness) ≠ coupledEnergy witness := by
  rw [witness_energy_after, witness_energy_before]
  norm_num

theorem noncommuting_substeps :
    swap (phaseFirst witness) ≠ phaseFirst (swap witness) := by
  intro h
  have hc := congrArg TwoSite.ar h
  norm_num [swap, phaseFirst, witness] at hc

end ProjectionFree.SplitEnergyDriftLimits

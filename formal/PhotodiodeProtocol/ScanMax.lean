import Mathlib
import PhotodiodeProtocol.Randomization

/-!
# Scan-max permutation inference over a declared frequency family (#6)

Source: `scan_npy_logcos.scan_logcos_cpu` (maximum of `delta_chi2` over the declared `k_grid`),
`scan_npy_logcos.shuffle_null` (each permuted spectrum is re-scanned over the *whole* grid and
its best `delta_chi2` is stored) and `scan_npy_logcos.p_value_geq`.

## Main results

* `scanMax`, `le_scanMax`, `scanMax_attained`, `scanMax_mono`: the scan statistic over a finite
  nonempty frequency family; it dominates every single-frequency statistic and can only grow
  when the family grows.
* `selected_pValue_le_scan` (**selection effect**): for any reference set, the single-frequency
  p-value evaluated at the data-selected frequency is never larger than the scan-max p-value.
  Reporting the selected bin's own p-value therefore overstates significance.
* `scanMax_randomization_level`: the scan-max randomization p-value is valid at every level,
  because the scan maximum is one fixed statistic (direct corollary of `randomization_level`).
* `minP_inflation`, `stat_marginally_valid`: an exact finite counterexample in which two
  single-frequency p-values are each exactly valid, but "reject if either is `≤ 1/2`" has
  null probability `3/4 > 1/2`.
* `exchangeable_swap_symmetric`, `not_exchangeable_of_strict_order`, `mcPValue_two_strict`:
  broken exchangeability. Under any exchangeable law `P(X₀ > X₁) = P(X₁ > X₀)`; a law with
  `X₀ > X₁` almost surely is therefore not exchangeable, and for it the Monte Carlo p-value is
  identically `1/2`, so a level-`1/2` test rejects with probability `1`.

These are statements about the declared statistics and finite reference sets only. Whether
bin permutation (`shuffle_null`) is a valid null for windowed physical spectra is an empirical
modelling assumption that Lean does not certify; neighbouring Hann-windowed bins are generally
correlated, which is exactly the kind of exchangeability failure illustrated here.
-/

namespace PhotodiodeProtocol.ScanMax

open Finset PhotodiodeProtocol

variable {Φ α : Type*}

/-! ## 1. The scan statistic -/

/-- Scan maximum over the declared finite frequency family `K`. -/
noncomputable def scanMax (K : Finset Φ) (hK : K.Nonempty) (T : Φ → α → ℝ) (a : α) : ℝ :=
  K.sup' hK fun k => T k a

theorem le_scanMax (K : Finset Φ) (hK : K.Nonempty) (T : Φ → α → ℝ) (a : α) {k : Φ}
    (hk : k ∈ K) : T k a ≤ scanMax K hK T a :=
  Finset.le_sup' (fun k => T k a) hk

theorem scanMax_attained (K : Finset Φ) (hK : K.Nonempty) (T : Φ → α → ℝ) (a : α) :
    ∃ k ∈ K, T k a = scanMax K hK T a := by
  obtain ⟨k, hk, h⟩ := Finset.exists_mem_eq_sup' hK (fun k => T k a)
  exact ⟨k, hk, h.symm⟩

/-- Enlarging the frequency family cannot decrease the scan statistic. -/
theorem scanMax_mono (K L : Finset Φ) (hK : K.Nonempty) (hL : L.Nonempty) (hKL : K ⊆ L)
    (T : Φ → α → ℝ) (a : α) : scanMax K hK T a ≤ scanMax L hL T a := by
  unfold scanMax
  exact Finset.sup'_le hK _ fun k hk => Finset.le_sup' (fun k => T k a) (hKL hk)

/-! ## 2. Selection effect and validity -/

/-- **Selection effect.** If frequency `k` attains the scan maximum on the observed data `a`,
the single-frequency upper-tail p-value at `k` is at most the scan-max p-value, for every
reference set `A` (e.g. all bin permutations). -/
theorem selected_pValue_le_scan (A : Finset α) (K : Finset Φ) (hK : K.Nonempty)
    (T : Φ → α → ℝ) (a : α) (k : Φ) (hk : k ∈ K) (hsel : T k a = scanMax K hK T a) :
    pValue A (T k) a ≤ pValue A (scanMax K hK T) a := by
  unfold pValue tailCount
  apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
  exact_mod_cast card_le_card fun b hb => by
    simp only [mem_filter] at hb ⊢
    exact ⟨hb.1, by rw [← hsel]; exact le_trans hb.2 (le_scanMax K hK T b hk)⟩

/-- The scan-max randomization p-value is valid at every level `α ≥ 0`. -/
theorem scanMax_randomization_level (A : Finset α) (K : Finset Φ) (hK : K.Nonempty)
    (T : Φ → α → ℝ) (lvl : ℝ) (hlvl : 0 ≤ lvl) (hne : A.Nonempty) :
    uniformProb A (fun a => pValue A (scanMax K hK T) a ≤ lvl) ≤ lvl :=
  randomization_level A _ lvl hlvl hne

/-! ## 3. Minimum single-bin p-values are not valid -/

/-- Natural-number valued indicator statistics for two "frequencies" on a uniform four-point
assignment space `Fin 2 × Fin 2`: frequency `0` looks at the first coordinate, frequency `1`
at the second. -/
def natStat (j : Fin 2) (ω : Fin 2 × Fin 2) : ℕ :=
  if (if j = 0 then ω.1 else ω.2) = 0 then 1 else 0

noncomputable def stat (j : Fin 2) (ω : Fin 2 × Fin 2) : ℝ := natStat j ω

theorem card_univ_omega : (univ : Finset (Fin 2 × Fin 2)).card = 4 := by
  rfl

theorem tail_stat (j : Fin 2) (ω : Fin 2 × Fin 2) :
    tailCount univ (stat j) ω = (univ.filter fun b => natStat j ω ≤ natStat j b).card := by
  unfold tailCount stat
  congr 1
  exact filter_congr fun b _ => Nat.cast_le

theorem pValue_stat_le_half_iff (j : Fin 2) (ω : Fin 2 × Fin 2) :
    pValue univ (stat j) ω ≤ 1 / 2 ↔
      (univ.filter fun b => natStat j ω ≤ natStat j b).card ≤ 2 := by
  unfold pValue
  rw [tail_stat, card_univ_omega, div_le_iff₀ (by norm_num)]
  norm_num

/-- Each single-frequency p-value is exactly valid under the uniform assignment. -/
theorem stat_marginally_valid (j : Fin 2) (lvl : ℝ) (hlvl : 0 ≤ lvl) :
    uniformProb univ (fun ω => pValue univ (stat j) ω ≤ lvl) ≤ lvl :=
  randomization_level univ (stat j) lvl hlvl univ_nonempty

/-- **Multiplicity counterexample.** Rejecting when *either* single-frequency p-value is
`≤ 1/2` has null probability `3/4`, exceeding the nominal level `1/2`. -/
theorem minP_inflation :
    uniformProb univ (fun ω : Fin 2 × Fin 2 =>
      pValue univ (stat 0) ω ≤ 1 / 2 ∨ pValue univ (stat 1) ω ≤ 1 / 2) = 3 / 4 := by
  unfold uniformProb
  simp only [pValue_stat_le_half_iff]
  rw [card_univ_omega]
  have h : (univ.filter fun ω : Fin 2 × Fin 2 =>
      (univ.filter fun b => natStat 0 ω ≤ natStat 0 b).card ≤ 2 ∨
        (univ.filter fun b => natStat 1 ω ≤ natStat 1 b).card ≤ 2).card = 3 := by
    decide
  rw [h]
  norm_num

theorem minP_exceeds_level :
    (1 / 2 : ℝ) < uniformProb univ (fun ω : Fin 2 × Fin 2 =>
      pValue univ (stat 0) ω ≤ 1 / 2 ∨ pValue univ (stat 1) ω ≤ 1 / 2) := by
  rw [minP_inflation]
  norm_num

/-! ## 4. Broken exchangeability -/

/-- In any exchangeable two-coordinate model, `P(X₁ < X₀) = P(X₀ < X₁)`. -/
theorem exchangeable_swap_symmetric (M : ExchangeableModel 2) :
    M.prob (fun ω => M.X ω 1 < M.X ω 0) = M.prob (fun ω => M.X ω 0 < M.X ω 1) := by
  classical
  unfold ExchangeableModel.prob
  rw [sum_filter, sum_filter]
  rw [← (M.act (Equiv.swap 0 1)).sum_comp]
  refine sum_congr rfl fun ω _ => ?_
  rw [M.w_act, M.X_act]
  simp [Function.comp, Equiv.swap_apply_left, Equiv.swap_apply_right]

/-- Two disjoint events have total probability at most one. -/
theorem prob_disjoint_le_one {n : ℕ} (M : ExchangeableModel n) (E F : M.Ω → Prop)
    [DecidablePred E] [DecidablePred F] (hdis : ∀ ω, ¬ (E ω ∧ F ω)) :
    M.prob E + M.prob F ≤ 1 := by
  classical
  unfold ExchangeableModel.prob
  rw [sum_filter, sum_filter, ← sum_add_distrib, ← M.w_sum]
  refine sum_le_sum fun ω _ => ?_
  have hw := M.w_nonneg ω
  by_cases hE : E ω <;> by_cases hF : F ω
  · exact absurd ⟨hE, hF⟩ (hdis ω)
  · simp [hE, hF]
  · simp [hE, hF]
  · simp [hE, hF, hw]

/-- A law with `X₁ < X₀` almost surely cannot be exchangeable. -/
theorem not_exchangeable_of_strict_order (M : ExchangeableModel 2)
    (h : M.prob (fun ω => M.X ω 1 < M.X ω 0) = 1) : False := by
  have hsym := exchangeable_swap_symmetric M
  have hle := prob_disjoint_le_one M (fun ω => M.X ω 1 < M.X ω 0)
    (fun ω => M.X ω 0 < M.X ω 1) (fun ω hω => lt_asymm hω.1 hω.2)
  rw [← hsym, h] at hle
  norm_num at hle

/-- …and for such a law the source's Monte Carlo p-value (`B = 1`) is identically `1/2`:
a level-`1/2` test rejects with probability one. -/
theorem mcPValue_two_strict (x y : ℝ) (h : y < x) :
    mcPValue (fun j : Fin 2 => if j = 0 then x else y) = 1 / 2 := by
  unfold mcPValue
  have hfilter : ((univ.erase (0 : Fin 2)).filter
      fun j => (if (0 : Fin 2) = 0 then x else y) ≤ if j = 0 then x else y) = ∅ := by
    apply filter_eq_empty_iff.mpr
    intro j hj
    have hj0 : j ≠ 0 := (mem_erase.mp hj).1
    simp only [if_neg hj0, not_le]
    exact h
  rw [hfilter]
  norm_num

end PhotodiodeProtocol.ScanMax

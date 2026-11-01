import Mathlib

/-!
# Exact finite randomization tests: type-I error control with ties

Everything here is finite and exact. No asymptotics and no simulation are involved.

* `card_small_upper_tail_le`, the **rank bound**: for any finite reference set `A` and any
  statistic `T`, the number of `a ∈ A` whose upper-tail count `#{b ∈ A : T b ≥ T a}` is at
  most `k` is itself at most `k`. Ties count toward the tail, so the test is conservative.
* `randomization_level`: under a **sharp null with random assignment**, where the assignment is
  uniform on `A` and the statistic is a fixed function of it, the exact p-value
  `p(a) = #{b : T b ≥ T a} / |A|` satisfies `P(p ≤ α) ≤ α` for every real `α`.
* `mc_pvalue_eq_tail`: the repository's Monte Carlo p-value
  `(1 + #{j ≠ 0 : X_j ≥ X_0}) / (1 + B)` (`scan_npy_logcos.p_value_geq`) equals the upper-tail
  rank of `X_0` among all `B + 1` values.
* `exchangeable_level`: if the law of `(X_0, …, X_B)` is **exchangeable**, that p-value satisfies
  `P(p ≤ α) ≤ α`. The law is modelled as a finite weighted sample space with a
  weight-preserving permutation action.
* `intersection_union_level`: if PASS needs *both* control comparisons to reject, and each
  p-value is valid, then `P(PASS)` is at most the level of either one.

Assumptions are explicit hypotheses. Nothing here says exchangeability or random assignment
holds for any physical apparatus. In particular, permuting spectral bins inside one waveform
(`shuffle_null`) is valid only if bins are exchangeable under the null. Windowed physical
spectra generally violate that, so the device-level randomization test is the primary one.
-/

namespace PhotodiodeProtocol

open Finset

variable {α : Type*}

/-- Upper-tail count of `a` within `A`; ties are included. -/
noncomputable def tailCount (A : Finset α) (T : α → ℝ) (a : α) : ℕ :=
  (A.filter fun b => T a ≤ T b).card

/-- **Rank bound.** At most `k` reference points have upper-tail count `≤ k`. -/
theorem card_small_upper_tail_le (A : Finset α) (T : α → ℝ) (k : ℕ) :
    (A.filter fun a => tailCount A T a ≤ k).card ≤ k := by
  classical
  set S := A.filter fun a => tailCount A T a ≤ k with hS
  rcases S.eq_empty_or_nonempty with h | hne
  · rw [h]; exact Nat.zero_le _
  · obtain ⟨a₀, ha₀S, hmin⟩ := S.exists_min_image T hne
    have hsub : S ⊆ A.filter fun b => T a₀ ≤ T b := by
      intro a ha
      have haA : a ∈ A := (mem_filter.mp ha).1
      exact mem_filter.mpr ⟨haA, hmin a ha⟩
    calc S.card ≤ (A.filter fun b => T a₀ ≤ T b).card := card_le_card hsub
      _ = tailCount A T a₀ := rfl
      _ ≤ k := (mem_filter.mp ha₀S).2

/-- Exact randomization p-value. -/
noncomputable def pValue (A : Finset α) (T : α → ℝ) (a : α) : ℝ :=
  (tailCount A T a : ℝ) / A.card

theorem pValue_pos (A : Finset α) (T : α → ℝ) {a : α} (ha : a ∈ A) : 0 < pValue A T a := by
  unfold pValue tailCount
  have h1 : 0 < (A.filter fun b => T a ≤ T b).card :=
    card_pos.mpr ⟨a, mem_filter.mpr ⟨ha, le_rfl⟩⟩
  have h2 : 0 < A.card := card_pos.mpr ⟨a, ha⟩
  positivity

theorem pValue_le_one (A : Finset α) (T : α → ℝ) (a : α) : pValue A T a ≤ 1 := by
  unfold pValue tailCount
  rcases Nat.eq_zero_or_pos A.card with h | h
  · simp [h]
  · rw [div_le_one (by exact_mod_cast h)]
    exact_mod_cast card_filter_le _ _

/-- Probability, under the uniform assignment on `A`, of the event `E`. -/
noncomputable def uniformProb (A : Finset α) (E : α → Prop) [DecidablePred E] : ℝ :=
  ((A.filter E).card : ℝ) / A.card

/-- **Type-I error control of the exact randomization test**, at every level `α ≥ 0`. -/
theorem randomization_level (A : Finset α) (T : α → ℝ) (lvl : ℝ) (hlvl : 0 ≤ lvl)
    (hne : A.Nonempty) :
    uniformProb A (fun a => pValue A T a ≤ lvl) ≤ lvl := by
  classical
  have hN : (0 : ℝ) < A.card := by exact_mod_cast card_pos.mpr hne
  set k : ℕ := ⌊lvl * A.card⌋₊ with hk
  have hevent : A.filter (fun a => pValue A T a ≤ lvl) ⊆
      A.filter (fun a => tailCount A T a ≤ k) := by
    intro a ha
    obtain ⟨haA, hp⟩ := mem_filter.mp ha
    refine mem_filter.mpr ⟨haA, ?_⟩
    have : (tailCount A T a : ℝ) ≤ lvl * A.card := by
      unfold pValue at hp
      rwa [div_le_iff₀ hN] at hp
    exact Nat.le_floor this
  have hcard := (card_le_card hevent).trans (card_small_upper_tail_le A T k)
  unfold uniformProb
  rw [div_le_iff₀ hN]
  calc ((A.filter fun a => pValue A T a ≤ lvl).card : ℝ) ≤ k := by exact_mod_cast hcard
    _ ≤ lvl * A.card := Nat.floor_le (by positivity)

/-- For a negative level the rejection event is empty (p-values are positive). -/
theorem randomization_negative_level (A : Finset α) (T : α → ℝ) (lvl : ℝ) (hlvl : lvl < 0) :
    uniformProb A (fun a => pValue A T a ≤ lvl) = 0 := by
  classical
  have : A.filter (fun a => pValue A T a ≤ lvl) = ∅ := by
    apply filter_eq_empty_iff.mpr
    intro a ha hle
    have := pValue_pos A T ha
    linarith
  simp [uniformProb, this]

/-! ## Monte Carlo p-value of `scan_npy_logcos.p_value_geq` -/

/-- `(1 + #{j ≠ 0 : X_j ≥ X_0}) / (1 + B)`, where index `0` is the observed statistic and
indices `1..B` are the null draws. -/
noncomputable def mcPValue {B : ℕ} (X : Fin (B + 1) → ℝ) : ℝ :=
  (1 + ((univ.erase (0 : Fin (B + 1))).filter fun j => X 0 ≤ X j).card : ℝ) / (B + 1)

theorem mc_pvalue_eq_tail {B : ℕ} (X : Fin (B + 1) → ℝ) :
    mcPValue X = pValue univ X 0 := by
  classical
  unfold mcPValue pValue tailCount
  have hsplit : (univ.filter fun j => X 0 ≤ X j) =
      insert 0 ((univ.erase (0 : Fin (B + 1))).filter fun j => X 0 ≤ X j) := by
    ext j
    by_cases hj : j = 0 <;> simp [hj]
  rw [hsplit, card_insert_of_notMem (by simp), card_univ, Fintype.card_fin]
  push_cast
  ring

/-! ## Exchangeability -/

/-- A finite weighted sample space carrying an exchangeable random vector `X : Ω → (Fin n → ℝ)`.
Permutations act on outcomes, preserve the weights, and permute the coordinates of `X`. -/
structure ExchangeableModel (n : ℕ) where
  Ω : Type
  [fin : Fintype Ω]
  [dec : DecidableEq Ω]
  w : Ω → ℝ
  w_nonneg : ∀ ω, 0 ≤ w ω
  w_sum : ∑ ω, w ω = 1
  act : Equiv.Perm (Fin n) → Ω ≃ Ω
  w_act : ∀ σ ω, w (act σ ω) = w ω
  X : Ω → Fin n → ℝ
  X_act : ∀ σ ω, X (act σ ω) = X ω ∘ σ

attribute [instance] ExchangeableModel.fin ExchangeableModel.dec

/-- Probability of an event in the model. -/
noncomputable def ExchangeableModel.prob {n : ℕ} (M : ExchangeableModel n)
    (E : M.Ω → Prop) [DecidablePred E] : ℝ :=
  ∑ ω ∈ univ.filter E, M.w ω

/-- All fibres of `σ ↦ σ 0` have the same size. -/
theorem card_perm_fibre_eq {n : ℕ} (i j : Fin (n + 1)) :
    ((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter fun σ => σ 0 = i).card =
      ((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter fun σ => σ 0 = j).card := by
  classical
  refine card_bij (fun σ _ => Equiv.swap i j * σ) ?_ ?_ ?_
  · intro σ hσ
    simp only [mem_filter, mem_univ, true_and] at hσ ⊢
    simp [Equiv.Perm.mul_apply, hσ]
  · intro a _ b _ h
    exact mul_left_cancel h
  · intro τ hτ
    refine ⟨Equiv.swap i j * τ, ?_, ?_⟩
    · simp only [mem_filter, mem_univ, true_and] at hτ ⊢
      simp [Equiv.Perm.mul_apply, hτ]
    · show Equiv.swap i j * (Equiv.swap i j * τ) = τ
      rw [← mul_assoc, Equiv.swap_mul_self, one_mul]

/-- Number of permutations of `Fin (n+1)` sending `0` to a prescribed `i` is `n!`. -/
theorem card_perm_fixed_image {n : ℕ} (i : Fin (n + 1)) :
    ((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter fun σ => σ 0 = i).card =
      n.factorial := by
  classical
  have htot : (univ : Finset (Equiv.Perm (Fin (n + 1)))).card =
      ∑ j : Fin (n + 1), ((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter
        fun σ => σ 0 = j).card :=
    card_eq_sum_card_fiberwise (f := fun σ : Equiv.Perm (Fin (n + 1)) => σ 0) (t := univ)
      (by simp)
  simp only [card_perm_fibre_eq _ i, sum_const, card_univ, Fintype.card_fin,
    Fintype.card_perm, smul_eq_mul, Nat.factorial_succ] at htot
  exact (Nat.eq_of_mul_eq_mul_left (Nat.succ_pos n) htot).symm

/-- Counting lemma: for a fixed vector `v`, the permutations `σ` for which index `0` of
`v ∘ σ` has p-value `≤ lvl` number at most `lvl · (n+1)!`. -/
theorem perm_count_le {n : ℕ} (v : Fin (n + 1) → ℝ) (lvl : ℝ) (hlvl : 0 ≤ lvl) :
    (((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter
        fun σ => pValue univ (fun j => v (σ j)) 0 ≤ lvl).card : ℝ) ≤
      lvl * (n + 1).factorial := by
  classical
  have hre : ∀ σ : Equiv.Perm (Fin (n + 1)),
      pValue univ (fun j => v (σ j)) 0 = pValue univ v (σ 0) := by
    intro σ
    unfold pValue tailCount
    congr 2
    exact_mod_cast (card_bij (fun j _ => σ j) (by simp)
      (fun _ _ _ _ h => σ.injective h)
      (fun b hb => ⟨σ.symm b, by simpa using hb, by simp⟩))
  simp only [hre]
  set G := (univ : Finset (Fin (n + 1))).filter fun i => pValue univ v i ≤ lvl with hG
  have hfib : ((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter
      fun σ => pValue univ v (σ 0) ≤ lvl).card = G.card * n.factorial := by
    rw [card_eq_sum_card_fiberwise (f := fun σ : Equiv.Perm (Fin (n + 1)) => σ 0)
      (t := G) (by intro σ hσ; simpa [hG] using hσ)]
    rw [sum_const_nat (m := n.factorial)]
    intro i hi
    have : ((univ.filter fun σ : Equiv.Perm (Fin (n + 1)) => pValue univ v (σ 0) ≤ lvl).filter
        fun σ => σ 0 = i) = univ.filter fun σ : Equiv.Perm (Fin (n + 1)) => σ 0 = i := by
      ext σ
      simp only [mem_filter, mem_univ, true_and]
      constructor
      · exact fun h => h.2
      · intro h; refine ⟨?_, h⟩; rw [h]; simpa [hG] using hi
    rw [this]
    exact card_perm_fixed_image i
  have hGc : (G.card : ℝ) ≤ lvl * (n + 1) := by
    have h1 := randomization_level (univ : Finset (Fin (n + 1))) v lvl hlvl univ_nonempty
    unfold uniformProb at h1
    simp only [card_univ, Fintype.card_fin] at h1
    rw [div_le_iff₀ (by positivity)] at h1
    push_cast at h1
    exact h1
  rw [hfib, Nat.factorial_succ]
  push_cast
  have hf : (0 : ℝ) ≤ n.factorial := Nat.cast_nonneg _
  nlinarith

/-- **Level of the Monte Carlo / exchangeability p-value.** If `(X_0, …, X_n)` is
exchangeable, then `P(p(X) ≤ α) ≤ α`, where `p` is the upper-tail rank of `X_0` (equivalently
`mcPValue`, by `mc_pvalue_eq_tail`). -/
theorem exchangeable_level {n : ℕ} (M : ExchangeableModel (n + 1)) (lvl : ℝ) (hlvl : 0 ≤ lvl) :
    M.prob (fun ω => pValue univ (M.X ω) 0 ≤ lvl) ≤ lvl := by
  classical
  set f : M.Ω → ℝ := fun ω => if pValue univ (M.X ω) 0 ≤ lvl then 1 else 0 with hf
  have hprob : M.prob (fun ω => pValue univ (M.X ω) 0 ≤ lvl) = ∑ ω, M.w ω * f ω := by
    unfold ExchangeableModel.prob
    rw [sum_filter]
    refine sum_congr rfl fun ω _ => ?_
    by_cases h : pValue univ (M.X ω) 0 ≤ lvl <;> simp [hf, h]
  have hinv : ∀ σ, ∑ ω, M.w ω * f ω = ∑ ω, M.w ω * f (M.act σ ω) := by
    intro σ
    rw [← (M.act σ).sum_comp (fun ω => M.w ω * f ω)]
    simp [M.w_act]
  have havg : ((n + 1).factorial : ℝ) * ∑ ω, M.w ω * f ω =
      ∑ ω, M.w ω * ∑ σ : Equiv.Perm (Fin (n + 1)), f (M.act σ ω) := by
    calc ((n + 1).factorial : ℝ) * ∑ ω, M.w ω * f ω
        = ∑ _σ : Equiv.Perm (Fin (n + 1)), ∑ ω, M.w ω * f ω := by
          simp [sum_const, card_univ, Fintype.card_perm]
      _ = ∑ σ : Equiv.Perm (Fin (n + 1)), ∑ ω, M.w ω * f (M.act σ ω) :=
          sum_congr rfl fun σ _ => hinv σ
      _ = ∑ ω, ∑ σ : Equiv.Perm (Fin (n + 1)), M.w ω * f (M.act σ ω) := sum_comm
      _ = _ := by simp [mul_sum]
  have hcount : ∀ ω, ∑ σ : Equiv.Perm (Fin (n + 1)), f (M.act σ ω) ≤
      lvl * (n + 1).factorial := by
    intro ω
    have : ∑ σ : Equiv.Perm (Fin (n + 1)), f (M.act σ ω) =
        (((univ : Finset (Equiv.Perm (Fin (n + 1)))).filter
          fun σ => pValue univ (fun j => M.X ω (σ j)) 0 ≤ lvl).card : ℝ) := by
      simp only [hf, M.X_act]
      rw [sum_boole]
      rfl
    rw [this]
    exact perm_count_le _ lvl hlvl
  have hfac : (0 : ℝ) < (n + 1).factorial := by exact_mod_cast Nat.factorial_pos _
  have hbound : ∑ ω, M.w ω * ∑ σ : Equiv.Perm (Fin (n + 1)), f (M.act σ ω) ≤
      lvl * (n + 1).factorial := by
    calc ∑ ω, M.w ω * ∑ σ : Equiv.Perm (Fin (n + 1)), f (M.act σ ω)
        ≤ ∑ ω, M.w ω * (lvl * (n + 1).factorial) :=
          sum_le_sum fun ω _ => mul_le_mul_of_nonneg_left (hcount ω) (M.w_nonneg ω)
      _ = lvl * (n + 1).factorial := by rw [← sum_mul, M.w_sum, one_mul]
  rw [hprob]
  rw [← havg] at hbound
  nlinarith

/-- **Intersection–union test.** PASS needs both p-values to be `≤ α`. If either comparison
is valid at level `α` under the uniform assignment, so is the joint rejection. No multiplicity
correction is needed. -/
theorem intersection_union_level (A : Finset α) (T₁ T₂ : α → ℝ) (lvl : ℝ) (hlvl : 0 ≤ lvl)
    (hne : A.Nonempty) :
    uniformProb A (fun a => pValue A T₁ a ≤ lvl ∧ pValue A T₂ a ≤ lvl) ≤ lvl := by
  classical
  have hN : (0 : ℝ) < A.card := by exact_mod_cast card_pos.mpr hne
  have hsub : A.filter (fun a => pValue A T₁ a ≤ lvl ∧ pValue A T₂ a ≤ lvl) ⊆
      A.filter (fun a => pValue A T₁ a ≤ lvl) := by
    intro a ha
    simp only [mem_filter] at ha ⊢
    exact ⟨ha.1, ha.2.1⟩
  have h1 := randomization_level A T₁ lvl hlvl hne
  unfold uniformProb at h1 ⊢
  calc _ ≤ ((A.filter fun a => pValue A T₁ a ≤ lvl).card : ℝ) / A.card := by
        apply div_le_div_of_nonneg_right _ hN.le
        exact_mod_cast card_le_card hsub
    _ ≤ lvl := h1

end PhotodiodeProtocol

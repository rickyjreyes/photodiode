import Mathlib
import PhotodiodeProtocol.Randomization
import PhotodiodeProtocol.Spectral

/-!
# Frozen decision protocol: PASS / FAIL / INCONCLUSIVE

**Design.** `m` devices are indexed by `Fin m`. Each carries a unique blinded label, an optional
baseline capture and an optional late (post-excitation) capture. Each capture has a full
acquisition manifest. Treatment groups are `induced`, `sham` and `virgin`, and the assignment
`a : Fin m → Group` is revealed only at unblinding. Scores are computed from captures alone. The
score function has **no access to the assignment**, which is how blinding is modelled.

**Frozen thresholds.** All thresholds live in a `FrozenProtocol` value that is a *parameter* of
`classify`. An evidence bundle carries only a protocol identifier. A bundle declaring any other
identifier is INCONCLUSIVE (`protocol_mismatch_inconclusive`). A bundle therefore cannot move a
threshold after the fact.

**Missing data.** An absent capture is missing. A *present* capture with an incomplete or
inconsistent manifest is a protocol violation, and the result is INCONCLUSIVE. Analysis is
complete-case. More than `maxMissing` missing devices, or fewer than `nMin` complete devices in
any group (including an absent control group), gives INCONCLUSIVE.

**Test.** For complete devices `C`, the reference set `refSet` contains every relabelling that
agrees with `a` off `C` and has the same group counts on `C`. For each of the comparisons
induced−sham and induced−virgin, the statistic is the mean-difference of persistence scores and
the p-value is exact (`pValue`). PASS needs both p-values `≤ α` *and* both effects `≥ δ`.

**Main theorem `pass_level`.** Assume the **sharp null**: the captures do not depend on the
assignment. Assume **random assignment**: `a` is uniform on `refSet`. Then
`P(classify = PASS) ≤ α`. Ties are handled conservatively.

Non-claims: a PASS is a statistical statement about this design. It does not establish a novel
optical-memory mechanism, any physical interpretation, or that random assignment was physically
carried out.
-/

namespace PhotodiodeProtocol.Decision

open Finset Spectral

inductive Group | induced | sham | virgin
  deriving DecidableEq, Fintype, Repr

inductive Outcome | pass | fail | inconclusive
  deriving DecidableEq, Repr

/-- Acquisition manifest of one capture. -/
structure Manifest where
  timestamp : ℕ        -- acquisition time, epoch seconds
  hardwareId : String
  firmware : String
  environment : String
  codeRevision : String

/-- One capture: an `N`-sample waveform plus its manifest. -/
structure Capture (N : ℕ) where
  samples : Fin N → ℝ
  manifest : Manifest

/-- Per-device record, before unblinding. -/
structure Device (N : ℕ) where
  blindedLabel : String
  excitationTime : ℕ
  baseline : Option (Capture N)
  late : Option (Capture N)

/-- Thresholds and constants fixed before data collection. -/
structure FrozenProtocol where
  protocolId : String
  codeRevision : String
  N : ℕ
  band : Finset ℕ
  ref : Finset ℕ
  alpha : ℝ
  delta : ℝ
  nMin : ℕ
  maxMissing : ℕ
  alpha_nonneg : 0 ≤ alpha
  nMin_pos : 0 < nMin

/-- Evidence bundle: devices plus the protocol id it claims to follow. -/
structure Bundle (P : FrozenProtocol) (m : ℕ) where
  declaredProtocol : String
  devices : Fin m → Device P.N

variable {P : FrozenProtocol} {m : ℕ}

def manifestOK (P : FrozenProtocol) (x : Manifest) : Prop :=
  x.hardwareId ≠ "" ∧ x.firmware ≠ "" ∧ x.environment ≠ "" ∧ x.codeRevision = P.codeRevision

instance (x : Manifest) : Decidable (manifestOK P x) := by unfold manifestOK; infer_instance

/-- A present capture is valid iff its manifest is complete and its timestamp lies on the
correct side of the excitation time. -/
def captureOK (P : FrozenProtocol) (d : Device P.N) : Prop :=
  (∀ c ∈ d.baseline, manifestOK P c.manifest ∧ c.manifest.timestamp < d.excitationTime) ∧
  (∀ c ∈ d.late, manifestOK P c.manifest ∧ d.excitationTime < c.manifest.timestamp)

instance (d : Device P.N) : Decidable (captureOK P d) := by
  unfold captureOK; infer_instance

def isComplete (d : Device P.N) : Prop := d.baseline.isSome ∧ d.late.isSome

instance (d : Device P.N) : Decidable (isComplete d) := by unfold isComplete; infer_instance

/-- Complete devices. -/
def completeSet (B : Bundle P m) : Finset (Fin m) :=
  univ.filter fun i => isComplete (B.devices i)

/-- Persistence score of a device. It uses captures only, never the assignment. -/
noncomputable def score (B : Bundle P m) (i : Fin m) : ℝ :=
  match (B.devices i).baseline, (B.devices i).late with
  | some b, some l => persistence b.samples l.samples P.band P.ref
  | _, _ => 0

def groupCount (C : Finset (Fin m)) (a : Fin m → Group) (g : Group) : ℕ :=
  (C.filter fun i => a i = g).card

noncomputable def groupMean (B : Bundle P m) (C : Finset (Fin m)) (a : Fin m → Group)
    (g : Group) : ℝ :=
  (∑ i ∈ C.filter (fun i => a i = g), score B i) / groupCount C a g

/-- Mean-difference statistic `induced − g`. -/
noncomputable def stat (B : Bundle P m) (g : Group) (a : Fin m → Group) : ℝ :=
  groupMean B (completeSet B) a .induced - groupMean B (completeSet B) a g

/-- Conditional randomization reference set. -/
def refSet (B : Bundle P m) (a : Fin m → Group) : Finset (Fin m → Group) :=
  univ.filter fun a' => (∀ i ∉ completeSet B, a' i = a i) ∧
    ∀ g, groupCount (completeSet B) a' g = groupCount (completeSet B) a g

theorem self_mem_refSet (B : Bundle P m) (a : Fin m → Group) : a ∈ refSet B a := by
  simp [refSet]

/-- The reference set is the same from every point inside it. -/
theorem refSet_eq_of_mem (B : Bundle P m) {a a' : Fin m → Group} (h : a' ∈ refSet B a) :
    refSet B a' = refSet B a := by
  simp only [refSet, mem_filter, mem_univ, true_and] at h
  ext a''
  simp only [refSet, mem_filter, mem_univ, true_and]
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨fun i hi => (h1 i hi).trans (h.1 i hi), fun g => (h2 g).trans (h.2 g)⟩
  · rintro ⟨h1, h2⟩
    exact ⟨fun i hi => (h1 i hi).trans (h.1 i hi).symm, fun g => (h2 g).trans (h.2 g).symm⟩

noncomputable def pVal (B : Bundle P m) (g : Group) (a : Fin m → Group) : ℝ :=
  pValue (refSet B a) (stat B g) a

/-- Gate conditions (everything except the test itself). -/
def gatesOK (B : Bundle P m) (a : Fin m → Group) : Prop :=
  B.declaredProtocol = P.protocolId ∧
  (∀ i, captureOK P (B.devices i)) ∧
  (∀ i j, (B.devices i).blindedLabel = (B.devices j).blindedLabel → i = j) ∧
  m - (completeSet B).card ≤ P.maxMissing ∧
  ∀ g, P.nMin ≤ groupCount (completeSet B) a g

open Classical in
/-- **The frozen classification function.** -/
noncomputable def classify (B : Bundle P m) (a : Fin m → Group) : Outcome :=
  if gatesOK B a then
    if pVal B .sham a ≤ P.alpha ∧ pVal B .virgin a ≤ P.alpha ∧
        P.delta ≤ stat B .sham a ∧ P.delta ≤ stat B .virgin a then .pass else .fail
  else .inconclusive

theorem classify_pass_iff (B : Bundle P m) (a : Fin m → Group) :
    classify B a = .pass ↔ gatesOK B a ∧
      (pVal B .sham a ≤ P.alpha ∧ pVal B .virgin a ≤ P.alpha ∧
        P.delta ≤ stat B .sham a ∧ P.delta ≤ stat B .virgin a) := by
  unfold classify
  split_ifs with h1 h2 <;> simp_all

theorem classify_inconclusive_iff (B : Bundle P m) (a : Fin m → Group) :
    classify B a = .inconclusive ↔ ¬ gatesOK B a := by
  unfold classify
  split_ifs with h1 h2 <;> simp_all

/-- A bundle that declares a different protocol (e.g. changed thresholds) is INCONCLUSIVE. -/
theorem protocol_mismatch_inconclusive (B : Bundle P m) (a : Fin m → Group)
    (h : B.declaredProtocol ≠ P.protocolId) : classify B a = .inconclusive :=
  (classify_inconclusive_iff B a).mpr fun hg => h hg.1

/-- An empty or under-filled control group (or any group) can never yield PASS or FAIL. -/
theorem missing_group_inconclusive (B : Bundle P m) (a : Fin m → Group) (g : Group)
    (h : groupCount (completeSet B) a g < P.nMin) : classify B a = .inconclusive :=
  (classify_inconclusive_iff B a).mpr fun hg => absurd (hg.2.2.2.2 g) (not_le.mpr h)

theorem absent_sham_inconclusive (B : Bundle P m) (a : Fin m → Group)
    (h : groupCount (completeSet B) a .sham = 0) : classify B a = .inconclusive :=
  missing_group_inconclusive B a .sham (h ▸ P.nMin_pos)

/-- A present capture with an incomplete manifest, or the wrong code revision, or a
timestamp on the wrong side of excitation, gives INCONCLUSIVE. -/
theorem manifest_violation_inconclusive (B : Bundle P m) (a : Fin m → Group) (i : Fin m)
    (h : ¬ captureOK P (B.devices i)) : classify B a = .inconclusive :=
  (classify_inconclusive_iff B a).mpr fun hg => h (hg.2.1 i)

theorem duplicate_label_inconclusive (B : Bundle P m) (a : Fin m → Group) (i j : Fin m)
    (hij : i ≠ j) (h : (B.devices i).blindedLabel = (B.devices j).blindedLabel) :
    classify B a = .inconclusive :=
  (classify_inconclusive_iff B a).mpr fun hg => hij (hg.2.2.1 i j h)

theorem too_many_missing_inconclusive (B : Bundle P m) (a : Fin m → Group)
    (h : P.maxMissing < m - (completeSet B).card) : classify B a = .inconclusive :=
  (classify_inconclusive_iff B a).mpr fun hg => absurd hg.2.2.2.1 (not_le.mpr h)

/-- The score of a device depends only on its own captures. In particular it does not depend
on the assignment, and relabelling which group a device is in cannot change any score. -/
theorem score_depends_only_on_captures (B B' : Bundle P m) (i : Fin m)
    (hb : (B.devices i).baseline = (B'.devices i).baseline)
    (hl : (B.devices i).late = (B'.devices i).late) : score B i = score B' i := by
  unfold score; rw [hb, hl]

/-- **Type-I error control of the whole decision rule.** Under the sharp null (fixed captures)
and uniform random assignment on the conditional reference set, `P(PASS) ≤ α`. -/
theorem pass_level (B : Bundle P m) (a₀ : Fin m → Group) :
    uniformProb (refSet B a₀) (fun a => classify B a = .pass) ≤ P.alpha := by
  classical
  set A := refSet B a₀ with hA
  have hne : A.Nonempty := ⟨a₀, self_mem_refSet B a₀⟩
  have hsub : A.filter (fun a => classify B a = .pass) ⊆
      A.filter (fun a => pValue A (stat B .sham) a ≤ P.alpha) := by
    intro a ha
    obtain ⟨haA, hp⟩ := mem_filter.mp ha
    refine mem_filter.mpr ⟨haA, ?_⟩
    have := ((classify_pass_iff B a).mp hp).2.1
    unfold pVal at this
    rwa [refSet_eq_of_mem B haA] at this
  have hlvl := randomization_level A (stat B .sham) P.alpha P.alpha_nonneg hne
  have hN : (0 : ℝ) < A.card := by exact_mod_cast card_pos.mpr hne
  unfold uniformProb at hlvl ⊢
  calc ((A.filter fun a => classify B a = .pass).card : ℝ) / A.card
      ≤ ((A.filter fun a => pValue A (stat B .sham) a ≤ P.alpha).card : ℝ) / A.card := by
        apply div_le_div_of_nonneg_right _ hN.le
        exact_mod_cast card_le_card hsub
    _ ≤ P.alpha := hlvl

end PhotodiodeProtocol.Decision

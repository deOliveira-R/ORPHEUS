# P6 Adjoint Verification Spec — `normal ↔ adjoint` carve

Verification strategy for (1) the **adjoint flux φ\*** (adjoint
k-eigenvalue + reciprocity) and (2) the **adjoint-weighted
(eigenvalue-consistent) homogenization** consumer. Design-time
artefact; no production code yet.

Author: test-architect. Pre-reads cited inline by `file:line`.

---

## 0. What the carve actually changes (structural facts, verified)

The adjoint algebra is **already fully plumbed**; P6 adds exactly two
leaf transposes.

- `_AdjointOperator` (`orpheus/numerics/operator.py:615`) swaps
  `apply`↔`apply_transpose` and implements the G-adjoint
  `A* = G⁻¹ Aᵀ G` (`:676`). `OperatorSum` (`:705`), `OperatorProduct`
  (`:801`), `ScaledOperator` (`:893`) all propagate
  `CAP_APPLY_TRANSPOSE` **iff both operands carry it** (`:762`, `:841`,
  `:900`). So `A_loss.H` falls out the moment every leaf advertises
  `apply_transpose`.
- **L, C, B already advertise `CAP_APPLY_TRANSPOSE`.** `StreamingOperator`
  (`streaming.py:321` — `{CAP_APPLY, CAP_APPLY_TRANSPOSE}`, the μ→−μ
  reverse-direction matvec `Lᵀ` at `:425`); `MultiplicationOperator`
  (C, self-adjoint diagonal); `SNBoundaryOperator` (`boundary.py:141`,
  propagates from its laws). This is what makes
  `tests/sn/operators/test_g_adjoint_reciprocity.py` already GREEN for
  `(L+C−B)`.
- **S and F do NOT.** `ScatteringOperator.capabilities =
  frozenset({CAP_APPLY})` (`scattering.py:411`); its docstring (`:108`)
  states "no `apply_transpose` (would require an adjoint Pℓ
  transposition that the current ORPHEUS solver does not need)".
  `FissionOperator.capabilities = frozenset({CAP_APPLY})`
  (`fission.py:167`); its docstring (`:96`) names the exact adjoint:
  "`FᵀΨ* = νΣf · (χ · Ψ*)`, structurally distinct (it sums over the
  **adjoint** energy distribution)". **The χ↔νΣf swap IS the canonical
  adjoint-fission trap.**

So the carve = `apply_transpose` + `CAP_APPLY_TRANSPOSE` on
`ScatteringOperator` and `FissionOperator`, then an adjoint posing in
`power_iteration`'s `EigenvalueSolver` boundary
(`eigenvalue.py:179` — the adjoint row is already a "documented future
seam" at `eigenvalue.py:28`).

### How the P6 transport-adjoint gates DIFFER from the existing G-adjoint reciprocity

`tests/sn/operators/test_g_adjoint_reciprocity.py` pins the **Hilbert
adjoint** `A.H = G⁻¹AᵀG` of the within-group composite `(L+C−B)` — a
pure linear-algebra identity (`⟨Aψ,φ⟩_G = ⟨ψ,A.Hφ⟩_G`, foundation,
exact). It is **metric/transpose correctness on the operators that
already have a transpose**. It says NOTHING about:

- **S† / F†** (not in `(L+C−B)`; they are the new transposes P6 adds);
- the **physics** of the adjoint (μ-reversal + kernel transpose
  g↔g′ + emission/production swap);
- **eigenvalues** (it is a within-group fixed identity, no `F`, no `k`).

The P6 gates are a STRICT SUPERSET layered on top: they reuse the same
G-adjoint reciprocity *machinery* (the independent `_g_inner`,
`test_g_adjoint_reciprocity.py:154`) but extend the operator from
`(L+C−B)` to the full loss `A_loss = L+C−S−B` and the production `F`,
and add the eigenvalue + transport-reciprocity claim layers the
algebraic identity cannot reach. **The existing test is a foundation
dependency of P6, not a substitute for it.**

---

## 1. Claim-layer + pillar gate (vv §1.5, MANDATORY before writing)

| Gate | Claim layer | Pillar | Reference (structurally independent) | Why it can prove the layer |
|------|-------------|--------|--------------------------------------|----------------------------|
| A1 adjoint-keff = forward-keff | **eigenvalue** | closed-form (algebraic identity: same char. poly) | the FORWARD keff (itself L1-anchored vs `kinf_and_spectrum_homogeneous`) — equality is exact | k is an eigenvalue of both A and A† |
| A2 reciprocity ⟨ψ\*,q⟩=⟨ψ,Σ_d⟩ | **model** (flux-shape, fixed-source) | closed-form (the reciprocity identity) | the identity itself; both halves are independent solves | bilinear duality, solver-internal-independent |
| A3 analytic ∞-medium adjoint spectrum | **eigenvalue + flux-shape** | closed-form | LEFT eigenvector of `M=A⁻¹F` via `np.linalg.eig(Mᵀ)` | the adjoint flux of the 0-D k-problem IS the left eigenvector |
| A4 bi-orthogonality ⟨ψ\*_i,Fψ_j⟩=0 (i≠j) | **flux-shape** (intrinsic law) | closed-form | the spectral decomposition (eig of M and Mᵀ) | the defining law of a forward/adjoint eigenbasis pair |
| C1 PG-vs-Galerkin discriminator (φ\*≠φ) | **model** | closed-form (hand arithmetic) | per-region Python loop with explicit φ\* weight | the bilinear `⟨φ\*,Σφ⟩/⟨φ\*,φ⟩` |
| C2 keff-preservation under coarse-resolve | **eigenvalue** (first-order stationary) | closed-form (the fine keff) | the fine-mesh keff (L1-anchored) | perturbation theory: error O(δφ²) |
| C3 Mode-11 routing sentinel | **structural** (not a value claim) | — (call-graph) | monkeypatch counter | proves φ\* reaches the test basis |

**Pillar discipline checks (all pass):**
- A1/A2/A3/A4/C2 are eigenvalue/flux claims paired with **closed-form**,
  NOT MMS. (MMS does NOT prove eigenvalues — vv §pillars. No row here
  uses MMS, by design: the adjoint k-problem has a closed-form 0-D
  reference, so MMS is unnecessary AND would be the wrong pillar.)
- Structural independence: every chain terminates in `np.linalg.eig`
  (A1/A3/A4), the reciprocity identity (A2), or explicit hand
  arithmetic (C1) — NOT another ORPHEUS solver, NOT a procedurally-
  different fold of the same identity (vv L11).
- A1's reference is the FORWARD keff, which looks like "another ORPHEUS
  solver" — but the claim is **exact algebraic EQUALITY** (k is a root
  of the same characteristic polynomial), not "agrees to tolerance".
  The forward keff is itself anchored to `kinf_and_spectrum_homogeneous`
  (`derivations/common/eigenvalue.py:48`), so A1's trust terminates in
  `np.linalg.eig`, not in the SN solver. This is L4-agreement promoted
  to L1 by an exactness argument — declare it as such.

---

## 2. Config-blindness audit (vv §0.6 — every gate's regime activates its term)

The adjoint-specific blindnesses (compounding the standing ones):

- **Symmetric `SigS` hides the scattering kernel transpose.** S† transposes
  the group-transfer `g↔g'`. If `SigS[g_from,g_to]` is symmetric, S†=S and
  the transpose is invisible (the ERR-002 family, vv anti-#2/#6). **Every
  reciprocity + adjoint-keff config MUST use ASYMMETRIC `SigS`** (strong
  downscatter, ~zero upscatter — the physical asymmetry).
- **χ ∝ νΣf hides the fission transpose.** F† swaps emission χ ↔ production
  νΣf. If the test mixture has χ proportional to νΣf (or 1-group, where
  both are scalars), F†=F up to scale and the swap is invisible. **Configs
  MUST have χ and νΣf NON-proportional across groups** (e.g. χ=[1,0]
  fast-emission, νΣf=[ν₁,ν₂] both nonzero — the standard 2G reactor shape
  already does this; verify per-mixture).
- **1-group nulls BOTH transposes AND the eigenvalue claim** (Cardinal
  Rule; vv §1-group). 1G: SigS, χ, νΣf all scalars → S†=S, F†=F, and
  k=νΣf/Σa is flux-shape-blind. **≥2G everywhere for A1–A4.** (A2
  reciprocity could technically run 1G as a smoke test, but the
  kernel-transpose teeth require ≥2G asymmetric — so 1G adds nothing.)
- **Flat flux (homogeneous ∞-medium) is FINE for A1/A3/A4** (the 0-D
  k-problem IS flat by construction — that is the closed-form reference's
  domain). But it nulls spatial redistribution, so **A1 must ALSO run
  heterogeneous** (the eigenvalue equality must hold with spatial
  structure live, where a spatial-adjoint bug — e.g. a wrong μ-reversal
  in Lᵀ at a curvilinear pole — could break it). A3/A4 are intrinsically
  0-D (they verify the ENERGY adjoint), so flat is correct there.
- **Self-adjoint problem nulls C1.** The PG-vs-Galerkin discriminator
  needs φ\*≠φ MATERIALLY. A symmetric/near-symmetric problem has φ\*≈φ
  and the discriminator is a dud (the same trap as
  `test_homogenization.py:243`'s "fixture too flat" guard). **C1's config
  MUST be strongly non-self-adjoint** (asymmetric SigS + strong absorber
  gradient) and MUST assert the discrimination actually fired
  (`assert not np.allclose(pg, galerkin)`).

**Minimum activating configs:**
- A1/A2/A3/A4 reference mixture: **Mixture A 2G AND 4G**, asymmetric
  SigS, non-proportional χ/νΣf. 4G mandatory — Mixture A 2G's spectrum is
  coincidentally flat (`[0.707,0.707]`,
  `test_reaction_rate_functional.py:24,153`) so its forward AND adjoint
  spectra coincide → A3/A4 are flux-shape-blind at 2G. **4G φ\* is
  genuinely non-flat AND φ\*≠φ** — it carries the energy-adjoint teeth.
- A1 heterogeneous leg: a 2-region slab (fuel/absorber) + a **sphere**
  (the μ-reversal in Lᵀ is non-trivial at the pole — reuses the
  `_make_sphere` builder from `test_g_adjoint_reciprocity.py:94`).
- C1/C2 config: heterogeneous slab, vacuum→reflective tilt (reuse the
  `test_homogenization.py:243` flux-tilt fixture), asymmetric SigS so
  φ\*≠φ.

---

## 3. PHASE 1 — φ\* correctness (the hard half). Order is a dependency chain.

Gates ordered so each rests only on already-passed lower gates. φ\* is
NOT "correct" until **P1.0 → P1.5 all pass**.

### P1.0 — Leaf transpose unit gates (L0, foundation). MUST PASS FIRST.

Before any composite/eigenvalue claim, pin S† and F† **in isolation**
against hand arithmetic — the term-level defense (vv §1, enumerate every
term).

**P1.0a — F† = `νΣf · (χ·ψ*)`** (`fission.py:96` names this exactly).
- Config: 2G + 4G, χ=[1,0,…] (or non-proportional), νΣf all-nonzero,
  single cell (energy-only — the transpose is purely in the group axis).
- Reference: explicit Python loop `out_g = νΣf_g · Σ_{g'} χ_{g'} ψ*_{g'}`
  (structurally independent — NOT a re-call of F.apply).
- Tolerance: `assert_array_almost_equal_nulp(nulp≈8)` (algebraic
  identity; exact up to reduction order).
- **Mutation (proves teeth):** swap χ↔νΣf in F†'s implementation → the
  loop reference and F† diverge O(1) at 2G with non-proportional χ/νΣf →
  RED. (This is THE canonical adjoint trap; the gate that catches it is
  the headline P1.0 deliverable.) Confirm under `-O`.
- **Anti-blindness:** the mutation is INVISIBLE if χ∝νΣf — so the gate's
  fixture MUST assert `not np.allclose(χ, νΣf/νΣf.sum()*χ.sum())` (a
  positive "the swap would actually move the answer here" precondition,
  mirroring `test_reaction_rate_functional.py:153`).

**P1.0b — S† = group-transpose of the scattering kernel.** S aggregates
P0 (`SigS[g_from,g_to]`) + Pℓ (`R∘Λ∘M`) + (n,2n). S† transposes each:
- **P0†:** `Q*_g = Σ_{g'} SigS[g,g'] φ*_{g'}` (note: forward is
  `Σ_{g'} SigS[g',g] φ_{g'}` — the transpose swaps the matmul to
  `SigS @ φ*` from `SigSᵀ @ φ`). Reference: explicit loop. **Mutation:**
  use `SigS` instead of `SigSᵀ` (forget the transpose) → RED with
  asymmetric SigS, invisible with symmetric (precondition-guard the
  asymmetry).
- **Pℓ†:** the §5.6 kernel is `R∘Λ∘M` (`scattering.py:649`). Its
  transpose is `Mᵀ∘Λᵀ∘Rᵀ`. Λ (per-ℓ group matmul) transposes the group
  axis; M and R swap roles (the frame already earns `R.H` from
  `_AdjointOperator` — Phase D of the frame carve; **reuse it**). This is
  the subtle one. Reference: a dense `Gᵀ` fold of the kernel matrix on a
  tiny mesh (the structurally-independent dense-probe oracle pattern,
  `diag_p42_adjoint_oracle.py`-style, cited at
  `test_g_adjoint_reciprocity.py:50`). **Mutation:** transpose Λ but not
  the M/R role-swap (or vice versa) → RED.
- Tolerance: nulp≈16 (the einsum chain is deeper).
- Claim layer: L0 term. Foundation.

> **Carve-elegance note (for the implementer, not a gate):** S† should
> fall out of the existing operator algebra — `kernel.H` where
> `kernel = R∘Λ∘M` gives `Mᵀ∘Λᵀ∘Rᵀ` via `OperatorProduct.H`
> (`operator.py:801`) IF Λ advertises `apply_transpose`. So the ONLY new
> leaf transpose to hand-write is `LegendreMomentScattering.apply_transpose`
> (the per-ℓ group matmul transpose, `scattering.py:211`) + the P0/(n,2n)
> matmul transpose. P1.0b's dense oracle pins that the composed `S.H`
> equals the hand `Gᵀ S Gᵀ` fold — catching a wrong composition even if
> each leaf transpose is individually right (the ERR-061 frame-consistency
> lesson: every component correct, the bug is in how they compose).

### P1.1 — Per-ordinate / per-group reciprocity residual (L0, foundation).

Extend `test_g_adjoint_reciprocity.py`'s `_g_inner` machinery from
`(L+C−B)` to the **full loss `A_loss = L+C−S−B`** (now that S† exists).
- Gate: `⟨A_loss ψ, φ⟩_G = ⟨ψ, A_loss.H φ⟩_G`, random non-flat ψ,φ,
  asymmetric SigS, ≥2G. Reuse `_g_inner` (independent metric, anti-R1).
- **This is the structural keystone:** it proves S† composes correctly
  into the full loss adjoint via `OperatorSum.H`, with the independent
  G-metric (NOT the production metric — `test_g_adjoint_reciprocity.py:20`
  anti-metric-blind discipline).
- Tolerance: rel < 1e-12 (algebraic).
- **Per-group, NOT weight-summed (vv L27):** evaluate the residual
  per-group-pair so a kernel-transpose bug in ONE group channel cannot
  be masked by cancellation in the weight sum. (The existing test sums;
  P6 adds a per-group breakdown because S† is where group-channel bugs
  live.)
- **Mutation:** the P1.0b mutations re-run here → RED (confirms the
  composite path, not just the leaf, is gated). PLUS: a `−S` vs `+S` sign
  error in the loss posing → RED.

### P1.2 — Reciprocity ⟨ψ\*,q⟩ = ⟨ψ,Σ_d⟩ (L1, the FIXED-SOURCE adjoint gate).

THE adjoint correctness gate, structurally independent of solver
internals.
- Setup (mirrors a known-correct forward problem — see §AUG(b)):
  forward solve `A_loss ψ = q` for a source/detector pair on a 2G
  ASYMMETRIC heterogeneous slab; adjoint solve `A_loss† ψ* = Σ_d`.
- Gate: `⟨Σ_d, ψ⟩ == ⟨ψ*, q⟩` (volume inner product). Reciprocity is
  EXACT for the discrete system (not just continuous) when ψ* solves the
  exact discrete adjoint — so tolerance = solver `inner_tol` × a small
  safety factor (≈ `SAFETY(10) × inner_tol`, per the regression-tolerance
  idiom), NOT a loose physical tolerance.
- **Config MUST expose a wrong adjoint:**
  - **Asymmetric 2G SigS** — a S† group-transpose error breaks
    reciprocity O(1) (with symmetric SigS it would PASS even with S†=S,
    a false green).
  - **A source/detector pair where q and Σ_d live in DIFFERENT groups**
    (q in fast, Σ_d in thermal) — forces the downscatter chain through
    S†; a transpose error in the WRONG direction (upscatter) breaks it.
    This is the config that makes the kernel-transpose un-hideable.
- **Mutations (each must redden THIS gate):**
  1. S† → S (drop the group-transpose): RED.
  2. χ↔νΣf in F† — N/A here (fixed-source, no F); covered by P1.0a/P1.3.
  3. Lᵀ uses +μ instead of −μ (no ordinate reversal in streaming
     transpose): RED (the streaming direction is the spatial half of the
     adjoint).
- Claim layer: model (flux-shape, fixed-source). Pillar: closed-form
  (the reciprocity identity). **NOT eigenvalue** — fixed-source.

### P1.3 — Adjoint-keff = forward-keff (L1, the EIGENVALUE gate).

- Config: **2G AND 4G homogeneous ∞-medium** (Mixture A) — the
  closed-form domain — PLUS a **2-region heterogeneous slab AND a
  sphere** (spatial-adjoint live).
- Gate: `k_adjoint == k_forward` to **rtol≈1e-10** (exact algebraic
  equality; the only slack is the two power-iterations' convergence
  floors — both run to `inner_tol`, so 1e-10 is conservative). For the
  ∞-medium legs, BOTH equal the closed-form `kinf_homogeneous`
  (`derivations/common/eigenvalue.py:18`) — a triple equality
  (k_fwd == k_adj == k_closed) that anchors the whole gate to
  `np.linalg.eig`.
- **Mutations:**
  1. F† = F (no χ↔νΣf swap): the adjoint k SHIFTS off the forward k at
     ≥2G non-proportional χ/νΣf → RED. (At 1G or χ∝νΣf it would stay
     equal — false green — hence the config discipline.)
  2. S† = S (asymmetric SigS): adjoint k shifts → RED.
  3. The heterogeneous/sphere leg with a wrong Lᵀ μ-reversal → adjoint k
     shifts (the ∞-medium leg would NOT catch this — spatial term nulled
     by flat flux; this is why the heterogeneous + sphere legs are
     mandatory, vv §0.6).
- **Config-blindness REJECTION criterion:** if the plan ships ONLY the
  ∞-medium leg, REJECT — a spatial-adjoint bug (wrong μ-reversal) is
  invisible to flat flux. Both legs required.
- Claim layer: eigenvalue. Pillar: closed-form.

### P1.4 — Analytic ∞-medium adjoint spectrum (L1, flux-shape of φ\*).

The adjoint flux of the 0-D k-problem has a CLOSED FORM: the **left
eigenvector** of `M = A⁻¹F`, i.e. the right eigenvector of `Mᵀ = Fᵀ A⁻ᵀ`.
- Reference: extend `kinf_and_spectrum_homogeneous`
  (`derivations/common/eigenvalue.py:48`) — it already builds
  `M = solve(A,F)` and `eig(M)`. Add the LEFT eigenvector (eig of `Mᵀ`,
  i.e. `eig(M.T)` dominant, sign-normalised). This is a tiny, structurally-
  independent (`np.linalg.eig`) addition — a Branch-1 closed-form
  reference per algebra-of-record.
- Gate: the SN adjoint solve on a homogeneous reflective ∞-medium (flat
  in space) → its per-group spectrum == the closed-form left eigenvector
  (rtol≈1e-8).
- **4G mandatory** (2G Mixture A coincidentally flat → adjoint spectrum
  flat too → flux-shape-blind, `test_reaction_rate_functional.py:153`).
  At 4G the left eigenvector φ\* is genuinely non-flat AND ≠ the right
  eigenvector φ — this gate verifies φ\* has the RIGHT energy shape, not
  merely "some non-flat shape".
- **Mutation:** F†=F (χ↔νΣf) → the SN adjoint spectrum matches the RIGHT
  eigenvector (=forward φ) instead of the LEFT → diverges from the
  reference at 4G → RED. (This is the single sharpest test of F†: it pins
  the adjoint flux SHAPE, not just the eigenvalue, against a closed form.)
- Claim layer: eigenvalue + flux-shape. Pillar: closed-form.

### P1.5 — Bi-orthogonality ⟨ψ\*_i, F ψ_j⟩ = 0 for i≠j (foundation, intrinsic law).

The project standard: every math-bearing type ships a test of its
DEFINING law (user memory: test intrinsic properties). The defining law
of a forward/adjoint eigenbasis pair is bi-orthogonality w.r.t. F.
- Config: ∞-medium 4G (the 0-D problem has ng distinct eigenmodes;
  compute all via `eig(M)` for φ_j and `eig(Mᵀ)` for ψ\*_i).
- Gate: the cross-Gram `⟨ψ*_i, F φ_j⟩` is DIAGONAL (off-diagonal ≈ 0 to
  rtol≈1e-10; diagonal ≠ 0). This is a closed-form check (all from
  `np.linalg.eig`), and it pins that the SN-computed φ\* (P1.4) is the
  TRUE left-eigenbasis, not an accidental match on the dominant mode
  alone.
- **Discriminator:** a φ\* that matches the dominant left eigenvector by
  luck (e.g. a symmetric-problem coincidence) but is wrong on the
  subdominant modes FAILS the off-diagonal check. So this gate is the
  "intrinsic law" backstop to P1.4's "dominant-mode value" check.
- **Mutation:** F†=F → ⟨φ_i, Fφ_j⟩ is NOT diagonal (the right
  eigenvectors are F-orthogonal only in the self-adjoint case) → RED at
  4G asymmetric.
- Claim layer: flux-shape (intrinsic). Pillar: closed-form. Foundation
  (no theory `:label:`; it is the algebraic law).

**φ\* is "correct" iff P1.0a, P1.0b, P1.1, P1.2, P1.3, P1.4, P1.5 are all
GREEN and every named mutation reddens its named gate under `-O`.** Only
then may Phase 2 (the consumer) be built.

---

## 4. PHASE 2 — adjoint-weighted homogenization consumer (the small half).

Built ONLY after Phase 1 closes. The forward homogenization is the
**Galerkin degenerate** (φ\*=φ) already shipped
(`test_homogenization.py`, `discrete_ordinates.rst:15909`); P6 makes it
genuinely Petrov-Galerkin with φ\*≠φ.

### C1 — PG-vs-Galerkin discriminator where φ\*≠φ MATERIALLY (L0, foundation).

The φ\*-vs-φ analogue of `test_frame.py:339`
(`test_petrov_galerkin_project_differs_from_galerkin_when_test_neq_trial`)
and `test_homogenization.py:243` (φV-vs-dV). Template = those two; the
NEW axis is **test = φ\*·1_R ≠ trial = φ·1_R**.
- Config: strongly non-self-adjoint heterogeneous slab (asymmetric SigS +
  absorber gradient) so φ\* and φ differ materially across the region.
- Gate: `Σ_R = ⟨φ*, Σφ⟩_R / ⟨φ*, φ⟩_R` (the bilinear,
  `discrete_ordinates.rst:15887` `sn-homogenization-bilinear`) ==
  explicit per-region Python loop with φ\* as the test weight
  (structurally independent — NOT a re-call of the production matmul,
  `test_homogenization.py:18`). AND assert it DIFFERS from both the
  forward φ-weighted (Galerkin degenerate) AND the dV (volume) average:
  `assert not np.allclose(Σ_adjoint_weighted, Σ_forward_weighted)`.
- **Dud-guard (mandatory):** assert the discrimination fired
  (`test_homogenization.py:283` "fixture too flat" pattern, generalised:
  here "φ\* not different enough from φ"). If φ\*≈φ the PG type is
  ceremony and C1 proves nothing — the gate MUST verify φ\*≠φ first.
- Tolerance: rtol≈1e-12 (hand arithmetic). Claim: model. Pillar:
  closed-form.

### C2 — keff-preservation under coarse-resolve (L2, first-order stationary).

The bilinear-preservation invariant, currently
`sn-homogenization-bilinear` documented-but-untested
(`docs/verification/matrix.rst:1004`).
- Setup: solve fine (k_fine) → adjoint-weighted homogenize onto coarse →
  re-solve coarse (k_coarse).
- Gate: `|k_coarse − k_fine|` is **first-order stationary** — the error is
  **O(δφ²)** (perturbation theory: adjoint weighting makes k stationary
  w.r.t. flux perturbations). Concretely: refine the fine mesh in 2-3
  steps; the `|k_coarse − k_fine|` gap must SHRINK QUADRATICALLY (or at
  least faster than the forward-weighted/Galerkin homogenization's gap on
  the SAME ladder — the adjoint-weighted gap must be strictly smaller,
  which IS the whole point of eigenvalue-consistent homogenization).
- **Tolerance / claim discipline:** this is a CONVERGENCE-ORDER claim
  (the gap → 0 faster than first order), NOT a value claim — do NOT assert
  `k_coarse == k_fine` to a fixed tolerance (that would be O(h²) to a
  possibly-wrong limit, vv anti-#5). Assert the ORDER (gap_ratio shrinks)
  AND anchor the limit: k_fine is L1 (P1.3-anchored). The pairing of
  "order shrinks" + "limit is independently correct" is what makes C2 a
  real gate (vv: convergence rate necessary, NOT sufficient).
- **Discriminator vs the forward case:** run the SAME ladder with forward
  (φ=φ\*) homogenization; the adjoint-weighted keff gap must be
  SMALLER/higher-order. If they are equal, either the problem is
  accidentally self-adjoint (REJECT the config) or φ\* is not actually on
  the test side (C3 catches that).
- **Mutation:** use φ (forward) instead of φ\* as the test weight → the
  keff gap degrades to first-order (loses stationarity) → the order check
  REDS. Claim: eigenvalue (first-order stationary). Pillar: closed-form
  (k_fine reference).

### C3 — Mode-11 routing sentinel: φ\* (not φ) reaches the test basis (structural).

Template = `test_homogenization.py:310`
(`test_homogenize_routes_through_the_petrov_galerkin_frame`). The forward
sentinel counts `WeightedIndicatorBasis.analyze` to prove φ moved onto the
test side. P6's analogue proves **φ\*** (the ADJOINT) is the weight.
- Sentinel: monkeypatch-count the construction of the test basis and
  assert the weight array it receives IS φ\* (not φ). Strictly: an
  in-process autouse wrap (vv Mode-11 sharpening) on the adjoint-flux
  reader, asserting counter>0 AND that the captured weight `np.array_equal`
  the adjoint flux, NOT the forward flux.
- **Why it is NOT vacuous:** a bit-identity-preserving regression that
  silently used φ (forward) as the weight would produce a DIFFERENT number
  than C1 expects (C1 catches the value), but on a NEAR-self-adjoint
  fixture C1's value gap could be small — C3 is the structural backstop
  that fires regardless of how close φ\*≈φ. Pair C1 (value) + C3
  (structural), the exact split Mode-11 exists for.
- **Mutation:** wire φ as the test weight → C3's `array_equal(weight, φ*)`
  RED (even when the value barely moves). Claim: structural. Not a pillar
  (call-graph).

---

## 5. Ordering summary (the gate dependency DAG)

```
P1.0a (F† leaf)  ─┐
P1.0b (S† leaf)  ─┤→ P1.1 (full-loss reciprocity, composite)
                  │      │
                  │      ↓
                  │   P1.2 (fixed-source reciprocity ⟨ψ*,q⟩=⟨ψ,Σ_d⟩)  [L1]
                  │      │
                  └──────┴→ P1.3 (adjoint-k = forward-k)   [L1 eigenvalue]
                                 │
                                 ↓
                           P1.4 (∞-medium adjoint spectrum vs left eigenvector) [L1]
                                 │
                                 ↓
                           P1.5 (bi-orthogonality, intrinsic law) [foundation]
                                 │
        ════════ φ* CERTIFIED ════════
                                 │
                                 ↓
                    C1 (PG≠Galerkin, φ*≠φ) ─→ C2 (keff stationary) ─→ C3 (Mode-11 sentinel)
```

- **Before φ\* can be called correct:** P1.0 → P1.5 GREEN + every mutation
  reddens its gate under `-O`.
- **Before the consumer is built:** φ\* certified (the whole P1 chain).
- **L1 vs L4:** ALL Phase-1 value gates are **L1** (closed-form /
  algebraic-identity references terminating in `np.linalg.eig` or the
  reciprocity identity). There is **no L4 gate** in this plan — and there
  should not be one until the full L1 chain holds. (A future MC adjoint or
  a second-code adjoint cross-check would be L4, valid ONLY after P1.5 —
  vv §L4-parallel-to-ladder.)

---

## 6. Failure-mode coverage (new rows for the test-design table)

The carve introduces failure modes not yet in the `vv-principles`
failure-mode table. **File these BEFORE delivering** (self-improvement
trigger 1):

| Failure mode | Test-design row that catches it |
|--------------|--------------------------------|
| **Adjoint-fission role swap NOT applied** (F†=F, χ↔νΣf missing) | P1.0a leaf gate + P1.4 spectrum-shape gate; config MUST have χ NON-proportional to νΣf (else invisible) |
| **Scattering kernel transpose NOT applied** (S†=S) | P1.0b leaf + P1.1 per-group reciprocity + P1.2 cross-group source/detector; config MUST have ASYMMETRIC SigS |
| **Streaming μ-reversal missing in Lᵀ** (adjoint uses +μ) | P1.3 heterogeneous + sphere legs (∞-medium flat-flux leg is BLIND to it) |
| **Adjoint-weighted homogenization uses φ not φ\*** | C1 value + C3 Mode-11 sentinel (pair: value catches it on non-self-adjoint, sentinel catches it structurally regardless) |
| **Composite S.H wrong even with correct leaf transposes** (frame-consistency, ERR-061 family) | P1.0b dense `Gᵀ S Gᵀ` oracle (catches a wrong composition that per-leaf gates pass) |

→ propose as `numerical-bug-signatures` Signature 11 (adjoint role-swap
invisibility under symmetric kernels) once a real bug is caught in
implementation; until then, skill-table rows above.

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
| §4.0 `adjoint=None` degenerate pin | **foundation** (behaviour-neutral) | closed-form (bit-identity by inheritance) | 0-ULP `array_equal` vs the shipped forward output | the `None` arm IS the shipped code path |
| C1 homogenize PG discriminator (w=φ\*⊙φ ≠ φ) | **model** | closed-form (hand arithmetic) | per-region Python loop with the **bilinear** weight φ\*⊙φ | the bilinear `⟨φ\*,Σφ⟩/⟨φ\*,φ⟩` ⟹ XS field weighted by φ\*⊙φ |
| C2 homogenize keff-gap COMPARATIVE order | **eigenvalue** (first-order stationary) | closed-form (the fine keff) | the fine-mesh keff (L1-anchored) + the forward-weighted gap on the SAME ladder | perturbation theory: adjoint gap O(δφ²) < forward gap O(δφ) |
| C3 homogenize Mode-11 CAPTURE sentinel | **structural** (not a value claim) | — (call-graph) | in-process weight CAPTURE, `array_equal(·, φ\*⊙φ)` | proves the sigma-frame weight IS φ\*⊙φ (not bare φ\*, not φ) |
| C4 condense discriminator (average moves / marginalize frozen) | **model** + **structural** | closed-form (hand loop) + bit-identity | hand loop (Gate A) + `array_equal` of the marginalize channels (Gate B) | adjoint enters ONLY the `average` morphism |
| C5 condense Mode-11 CAPTURE sentinel | **structural** | — (call-graph) | in-process spectrum-weight CAPTURE | proves the average-frame spectrum weight IS φ\*_s⊙φ_s |
| Cχ χ-simplex positive control | **foundation** (intrinsic law) | closed-form (simplex constraint) | `EmissionSpectrum` validation (raises under -O) | a convex average of simplices is a simplex |

**Pillar discipline checks (all pass):**
- A1/A2/A3/A4/C2 are eigenvalue/flux claims paired with **closed-form**,
  NOT MMS. (MMS does NOT prove eigenvalues — vv §pillars. No row here
  uses MMS, by design: the adjoint k-problem has a closed-form 0-D
  reference, so MMS is unnecessary AND would be the wrong pillar.)
- Structural independence: every chain terminates in `np.linalg.eig`
  (A1/A3/A4), the reciprocity identity (A2), or explicit hand
  arithmetic (C1/C4 per-region/per-group loops with the bilinear weight)
  — NOT another ORPHEUS solver, NOT a procedurally-different fold of the
  same identity (vv L11). The homogenize/condense CAPTURE sentinels
  (C3/C5) and the §4.0 pins are structural (array identity / call-graph),
  not value claims, so they need no pillar.
- **C2 is COMPARATIVE, not a value claim.** Its reference (k_fine) IS an
  ORPHEUS solver — but C2 asserts the RATE relation `gap_adj(P) <
  gap_fwd(P)` and its faster shrinkage, NOT `k_coarse == k_fine` to a
  tolerance. k_fine is the L1-anchored LIMIT (P1.3-certified) the gaps
  shrink toward (vv anti-#5). The forward-weighted gap on the SAME ladder
  is the structurally-matched control that isolates the XS-collapse worth
  (§4 honest-scope). So C2's trust does not terminate in the SN solver's
  absolute value — it terminates in the perturbation-theory ORDER
  relation, anchored by an L1 limit.
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
- **Self-adjoint problem nulls the whole Phase-2 consumer.** The
  adjoint-weighted collapse only DIFFERS from the shipped forward when
  φ\*≠φ MATERIALLY. Reactor-physics keystone re-verified this session:
  the shipped forward `Solution.homogenize` (`orpheus/sn/solution.py:641-644`)
  is the **φ\*=1 degenerate** (sigma-frame test weight = φ, single power,
  `Σ_R = Σ_i V_i φ_i Σ_i / Σ_i V_i φ_i`); the bilinear
  `sn-homogenization-bilinear` (`docs/theory/foundations/frame.rst:999`,
  `Σ_R = ∫_R φ*Σφ dV / ∫_R φ*φ dV`) weights the XS field by the
  **elementwise product w = φ\*⊙φ** (the FIXED rule for the VECTOR
  channels). Two consequences the gates MUST encode:
  - **The bilinear is φ\*-SCALE-INVARIANT** (`∫φ*Σφ/∫φ*φ` unchanged by
    φ\*→c·φ\*), and φ\*, φ carry independent eigenvector normalizations,
    so any "φ\*≠φ" dud-guard MUST compare NORMALIZED SHAPES
    (`φ*/‖φ*‖` vs `φ/‖φ‖`), never raw arrays.
  - **The doc's "φ→φ\*" shorthand (`frame.rst:3458`) is a TRAP.** The
    correct sigma-frame weight is φ\*⊙φ, NOT bare φ\*. An implementer who
    reads "a no-op change of the test_basis (φ→φ\*)" literally writes
    `WeightedIndicatorBasis(trial, φ*)` — WRONG. C3/C5 pin the weight IS
    φ\*⊙φ by CAPTURE-and-compare (a candidate new failure mode, §6).
  Configs MUST be strongly non-self-adjoint (asymmetric SigS + absorber
  gradient) AND assert the discrimination fired (the
  `tests/sn/test_homogenization.py:284` "fixture too flat" pattern,
  generalised to the normalized-shape guard).
- **The USER RULING at P6-open (2026-07-25) — derivation-first channel
  taxonomy.** B1 ships a `derivations` module settling the collapse rule
  PER CHANNEL CLASS: **vector** w=φ\*⊙φ (fixed above); **matrix** channels
  (SigS, Sig2) per-pair `φ*_g·φ_{g'}` vs source-product — the B1
  derivation adjudicates; the factored **fission-dyad** χ rule. Gates
  below reference "**the B1-derived rule**" where the weight is not yet
  fixed (matrix source axis, χ), but make gate STRUCTURE (configs,
  references, mutations, tolerances) fully concrete NOW.
- **The condense MARGINALIZE morphism is adjoint-INVARIANT — a free
  structural discriminator.** `Mixture.condense`
  (`orpheus/data/macro_xs/mixture.py:388-408`) has TWO morphisms:
  `average` (`frame.project`, spectrum test weight — vectors + matrix
  SOURCE `g_from` axis) and `marginalize` (bare `@ table`, NO weight —
  matrix SINK `g_to` axis + χ). The adjoint enters ONLY the `average`
  morphism (no test weight on `marginalize` to adjoint-weight), so the
  condense discriminator (C4) gets a bit-identity backstop: adjoint- vs
  forward-condensed MUST differ on the average-collapsed channels AND be
  BYTE-IDENTICAL on the marginalize-collapsed channels. Needs a
  WITHIN-GROUP-VARYING spectrum (else the fine→coarse overlap weighting is
  flat and the average morphism does not discriminate).

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
- C1/C3 config (homogenize, L0): strongly non-self-adjoint heterogeneous
  slab, vacuum→reflective tilt (the `tests/sn/test_homogenization.py:254`
  flux-tilt pattern), asymmetric SigS + absorber gradient so φ\*≠φ in
  SHAPE. The module fixture (`test_homogenization.py:58`: 2G, χ=[1,0],
  νΣf=[0.24,0.48] non-proportional, SigS=[[0.60,0.10],[0.0,0.90]]
  asymmetric) already carries the SigS/χ asymmetry; C1 adds the absorber
  gradient for a materially non-flat φ\*.
- C2 config (homogenize, L2 — GREENFIELD, no existing test re-solves a
  homogenized mesh): a **16-cell** ≥2G heterogeneous asymmetric-SigS slab
  with vacuum→reflective tilt; coarse-partition LADDER P ∈ {2, 4, 8}
  coarse cells (each a contiguous union of fine cells — 16 divisible by
  all three; identity partition P=16 is exact). Re-solve
  `solve_sn(mm.materials, mm.mesh, quad, …)` (`MaterialMesh` carries
  `.materials` dict + `.mesh` coarse geometry).
- C4/C5 config (condense, L0): 4-fine-group → 2-coarse-group (reuse
  `tests/sn/test_condensation.py:_EG_FINE`/`_EG_COARSE`), asymmetric SigS,
  WITHIN-GROUP-VARYING representative spectrum, strongly non-self-adjoint
  so φ\*_spectrum differs in shape from φ_spectrum.

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

> **CORRECTION (2026-07-24, #276 A4 execution — caught by the SN daggered
> solve on first contact):** this section's original formula was WRONG.
> φ\* solves the daggered eigenproblem `Aᵀφ* = Fᵀφ*/k`, i.e. it is the
> dominant right eigenvector of **`(Aᵀ)⁻¹Fᵀ`** — the left eigenvector of
> the REVERSED product `FA⁻¹`, **NOT** of `M = A⁻¹F`. `Mᵀ = FᵀA⁻ᵀ` is
> merely SIMILAR to `(Aᵀ)⁻¹Fᵀ` (conjugation by `Aᵀ` — the same algebra as
> the #226 step-5b `A·(A⁻¹F)·A⁻¹ = FA⁻¹` Mode-12 finding), so every
> k-level check passes on both; but for the rank-1 `F` the wrong product's
> dominant eigenvector degenerates to EXACTLY `ν̂Σf` (`Fᵀx ∝ νΣf` for all
> x) — a reference carrying ZERO A-physics. Measured (Mixture A): the
> spec-formula vector ≡ ν̂Σf to all digits at 2G AND 4G; the corrected
> vector = [0.684, 0.730] (2G) / [0.470, 0.486, 0.518, 0.524] (4G), which
> the SN daggered solve reproduces digit-for-digit. The reference
> implementation + its pin (`kinf_and_adjoint_spectrum_homogeneous`,
> `tests/derivations/test_adjoint_spectrum_reference.py`) carry the
> corrected law, the defining-law pin `Aᵀφ* = Fᵀφ*/k`, and a
> ν̂Σf-degeneracy trap-catcher row. The original text below is kept for
> the record; read it with this correction.

The adjoint flux of the 0-D k-problem has a CLOSED FORM: ~~the **left
eigenvector** of `M = A⁻¹F`, i.e. the right eigenvector of `Mᵀ = Fᵀ A⁻ᵀ`~~
(corrected above: the dominant right eigenvector of `(Aᵀ)⁻¹Fᵀ`).
- Reference: extend `kinf_and_spectrum_homogeneous`
  (`derivations/common/eigenvalue.py:48`) — it already builds
  `M = solve(A,F)` and `eig(M)`. Add the adjoint eigenvector
  (`dominant_eigenpair(solve(A.T, F.T))`, sign-normalised). This is a
  tiny, structurally-independent (`np.linalg.eig`) addition — a Branch-1
  closed-form reference per algebra-of-record.
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

## 4. PHASE 2 — adjoint-weighted homogenization + condensation consumer (the small half).

> **Delta-refresh (2026-07-25, P6 open).** Re-verified against the landed
> carrier (`orpheus/sn/solution.py`) and the shipped forward frames.
> Phase 1 (§3) is GREEN and merged — φ\* is CERTIFIED (P1.0–P1.5 + the
> sphere φ\*-shape row); Phase 2 is now LIVE. Six substantive corrections
> vs the pre-carve draft:
> **(a)** C2 ladder direction REVERSED — the stationarity limit is the
> COARSE partition → the FIXED fine reference (δφ→0), NOT fine-mesh
> refinement under a fixed partition (that leaves a finite homogenization
> defect and no measurable rate); the claim is now COMPARATIVE
> (adjoint gap < forward gap, per #281 acceptance), C2 is greenfield.
> **(b)** C3 upgraded from COUNT to CAPTURE — the shipped Mode-11 sentinel
> (`tests/sn/test_homogenization.py:311-369`) only counts `.analyze`;
> C3 must CAPTURE the weight and `array_equal` it to φ\*⊙φ.
> **(c)** C1 dud-guard made NORMALIZATION-AWARE — the bilinear is
> φ\*-scale-invariant, so compare normalized shapes.
> **(d)** Condense-axis gates ADDED (C4/C5) — the ratified API covers
> `condense`; the adjoint enters ONLY the `average` morphism.
> **(e)** Degenerate pins made EXPLICIT (§4.0) — `adjoint=None` is the
> shipped path, bit-identical.
> **(f)** χ-simplex positive control ADDED (Cχ).
> Stale paths fixed: `tests/sn/test_homogenization.py`,
> `tests/numerics/test_frame.py:358`,
> `docs/theory/foundations/frame.rst:999`/`:3423`,
> `docs/theory/verification/matrix.rst`.

> **B3 implementation addendum (2026-07-26, gates landed).** Three
> corrections measured AGAINST this refresh while wiring the gates, plus
> the teeth record (probe script `b3_mutation_probes.py`, session tmp):
>
> **(g) C2 re-redesigned: SAME-MESH replacement + CONTRAST ladder.** The
> coarse-RESOLVE comparative of (a) FAILED measurably (adjoint gap 9.50e-3
> vs forward 9.41e-3 at P=2; equal at P=8): (i) the coarse-mesh DD
> discretization deviation between the two arms is NOT negligible at
> coarse P and inverts the ~1e-4 worth delta; (ii) single-material coarse
> regions (the [0,0,1,1] tiling at P=8) null the weight entirely; (iii) a
> partition ladder holds within-region heterogeneity CONSTANT under an
> alternating pattern — P is not the smallness knob. The landed C2
> (`TestC2ComparativeKeffOrder`): re-solve the FINE 16-cell mesh with the
> region-collapsed materials (the same-mesh XS-replacement worth — exactly
> what T0/T5 prove; zero discretization confounder) on a MATERIAL-CONTRAST
> ladder m1 = m0 + ε·Δ, ε ∈ {1, ½, ¼}. Measured: forward ratios 2.05/2.01
> (first order), adjoint ratios 6.08/9.24 (≥ second order), adjoint gap
> smaller on every rung. Gate: per-rung adjoint < forward; ratio_fwd < 3;
> ratio_adj > 3.
> **(h) C4 Gate B superseded by the B&G convention.** Correction (d)'s
> "the adjoint enters ONLY the average morphism; marginalize FROZEN"
> encoded the pre-literature expectation. B&G Ch. 6 (memo Source E) +
> derivation T6: the sink axis GAINS the flux-weighted-average adjoint
> carrier Ψ†_G (B&G 6.136) and χ takes the χ†/s factored rule — the sink
> morphism is NOT frozen under the adjoint. Landed C4: every channel ==
> the B&G hand rule (nested blocks), discrimination vs forward, dud-guard;
> the frozen-χ/sink expectation is replaced by the hand-rule equalities.
> C5 amended likewise: the plain-φ frame legitimately appears (the B&G
> source-axis / νΣf flux-average); the sentinel forbids BARE φ\* only.
> **(i) Teeth record (all probes RED as designed).** C2 reds on
> adjoint-DROPPED (the forward arm's own ratios ~2 < 3 violate
> ratio_adj > 3). The φ\*≡φ (φ²-weight) substitution is NOT C2's catcher —
> measured ratio 4.98 (super-first-order: the within-region fluctuation of
> ANY solution-derived weight is O(ε) on the contrast fixture); its
> committed catchers are C1's hand-rule equality + C3's capture, exactly
> the C1 mutation-row assignment. C3 reds on bare-φ\* AND forward-φ weight
> mutations (both probed via in-process `project_through_bilinear` wraps).
> C4 reds on adjoint-dropped at the `Mixture.condense` level. Cχ raises on
> the χ-FREE-group sign flip (χ₃=0: Ψ†₁ flips, the χ†₁ numerator cannot);
> blockwise-uniform and χ-supported flips SELF-CANCEL through the χ†/s
> rescale — a recorded structural sign-robustness of the factored rule.

Built ONLY after Phase 1 closes (φ\* certified). The **RATIFIED #281 P6-B2
API** (`solution.py:484-490` — `Solution.homogenize(coarse, *, adjoint:
AdjointSolution | None = None)`, SAME for `condense`) lands the
Petrov-Galerkin implementation WITH these gates — there is no
advertised-but-unwired arm in between. `adjoint=None` ⇒ today's
flux-weighted (Galerkin, φ\*=φ *reading* of the φ\*=1 degenerate)
collapse, BIT-IDENTICAL. Obtain φ\* via
`solve_sn_adjoint(materials, fine_mesh, quad, …)` (`solver.py:2348` →
`AdjointSolution`); pass it as `adjoint=` to the FORWARD solution's
`homogenize`/`condense`.

**The weight, made precise (drives every gate).** The XS field is weighted
by the bilinear product **w = φ\*⊙φ** for the VECTOR channels (fixed);
the MATRIX-channel and χ (fission-dyad) weights are settled by the **B1
derivations module** (USER RULING, §2). Gates are concrete on the vector
rule NOW and reference "the B1-derived rule" for matrix/χ — with fully
concrete configs, references, mutations, and tolerances regardless.

### 4.0 — Degenerate pins: `adjoint=None` is the SHIPPED path, bit-identical (foundation).

The `None` arm MUST keep the EXACT current code path (no reshuffled
reduction), so the pin has TWO teeth:
- **Tooth 1 — the no-arg default equals the explicit `None`.**
  `sol.homogenize(coarse)` ≡ `sol.homogenize(coarse, adjoint=None)` —
  `np.array_equal` on every channel of every coarse `Mixture` (0-ULP; ANY
  rtol here is a RED FLAG). SAME for `condense` over the `dict[int,
  Mixture]`. Uses the `tests/numerics/test_frame.py:343-354` 0-ULP
  template. Catches a `None` arm that diverges from the default.
- **Tooth 2 — no shared drift from the pre-P6 baseline.** Tooth 1 is
  blind to a refactor that routes BOTH arms through one new
  `_weighted(weight)` helper whose reduction tree differs from the
  shipped `frame.project` (both drift together, stay equal to each
  other). The GUARD is the EXISTING forward suite staying green — its
  rate gates carry structurally-independent HAND references
  (`test_homogenization.py:244-287` φV loop; the identity test :194), so
  a shared reduction-tree drift reds THEM even when Tooth 1 passes.
- **Claim:** foundation (behaviour-neutral relabel — vv Mode-12: the
  neutrality is proven by DIRECT value comparison, not a proxy). **Pillar:**
  closed-form (bit-identity by inheritance).
- **Mutations (each REDs under -O):**
  1. Make the `None` arm apply φ² (the φ\*=φ misreading, weight φ⊙φ)
     while the default stays φ → Tooth 1 `array_equal` REDs (a real value
     difference, NOT a multiply-by-1.0 no-op, which IEEE-754 leaves
     bit-exact).
  2. Reshape/reorder the shared reduction so the pre-P6 forward values
     drift → Tooth 2 (the φV hand-reference rate gate) REDs.

### C1 — homogenize PG discriminator: the vector-channel weight is φ\*⊙φ ≠ φ (L0, foundation).

The φ\*-vs-φ analogue of `tests/numerics/test_frame.py:358`
(`test_petrov_galerkin_project_differs_from_galerkin_when_test_neq_trial`)
and `tests/sn/test_homogenization.py:244-287` (the φV-vs-dV dud-guard).
Template = those two; the NEW axis is the bilinear weight **w = φ\*⊙φ** on
the test side (`sn-homogenization-bilinear`,
`docs/theory/foundations/frame.rst:999`) — NOT the forward w=φ, NOT dV's
w=1.
- **Config:** strongly non-self-adjoint heterogeneous slab (asymmetric
  SigS + absorber gradient, vacuum→reflective tilt) so φ\* and φ differ in
  SHAPE across the region. `fwd = solve_sn(...)`,
  `adj = solve_sn_adjoint(...)`.
- **Gate (VECTOR channels — concrete NOW):** for SigT (and each of SigC,
  SigL, SigF, SigP),
  `fwd.homogenize(coarse, adjoint=adj).materials[R].SigT[g]` equals the
  explicit per-region Python loop with weight w = φ\*⊙φ:
  `Σ_R = Σ_{i∈R} V_i φ*_{i,g} φ_{i,g} Σ_{i,g} / Σ_{i∈R} V_i φ*_{i,g} φ_{i,g}`
  (structurally independent — NOT a re-call of the production
  `project_through`, per `test_homogenization.py:18` / vv L11). AND assert
  it DIFFERS from BOTH the forward φ-weighted degenerate (`adjoint=None`)
  AND the dV volume average:
  `assert not np.allclose(Σ_adjoint, Σ_forward)` and `… != Σ_dV`.
- **Gate (MATRIX + χ channels — structure NOW, weight B1-derived):** the
  same production-vs-hand-loop equality, but the hand loop uses **the
  B1-derived rule** (matrix: per-pair `φ*_g·φ_{g'}` vs source-product; χ:
  the factored fission-dyad rule). Do NOT hard-code φ\*⊙φ for the matrix
  SOURCE axis before B1 confirms it — the source-axis rule may legitimately
  differ from the vector rule. Ship config/mutation/tolerance now.
- **Dud-guard (mandatory, NORMALIZATION-AWARE — correction c):** the
  bilinear `∫φ*Σφ/∫φ*φ` is INVARIANT under φ\*→c·φ\*, and φ\*, φ carry
  independent eigenvector normalizations, so compare SHAPES:
  `assert not np.allclose(φ*/‖φ*‖, φ/‖φ‖, rtol=…)` per group — φ\* must be
  a materially different importance SHAPE (else the "adjoint weighting" is
  a rescaling of the forward and C1 proves nothing). PLUS assert the
  OUTPUT discrimination fired (`Σ_adjoint ≠ Σ_forward`) — the
  `test_homogenization.py:284` "fixture too flat" pattern generalised.
- **Tolerance:** rtol 1e-12 (hand arithmetic). **Claim:** model. **Pillar:**
  closed-form (per-region hand loop).
- **Mutations (each REDs C1 under -O):**
  1. Production ignores `adjoint` (uses w=φ, the forward degenerate) →
     production ≠ φ\*⊙φ hand reference → RED.
  2. Production uses BARE φ\* (`WeightedIndicatorBasis(trial, φ*)` — the
     `frame.rst:3458` "φ→φ\*" trap) → w=φ\* ≠ φ\*⊙φ → RED.
  3. Production uses φ² (adjoint weight taken as φ, i.e. the φ\*=φ
     misreading) → ≠ φ\*⊙φ hand reference at a genuinely non-self-adjoint
     config → RED.
- **Config-blindness:** self-adjoint (φ\*≈φ in shape) → the
  normalized-shape dud-guard fires FIRST (REJECT the config), so C1 can
  never false-green on a dud.

### C2 — homogenize keff-gap COMPARATIVE order under coarse-partition refinement (L2). REDESIGNED.

The eigenvalue-consistency payoff: adjoint weighting keeps k stationary
w.r.t. the within-region flux perturbation δφ (first-order perturbation
theory), so the coarse-resolve keff gap is HIGHER-ORDER than the
forward-weighted gap. Verifies `sn-homogenization-bilinear`
(documented-untested, `docs/theory/verification/matrix.rst`).
- **Ladder direction (correction a — the pre-carve draft was WRONG).**
  δφ shrinks as the COARSE partition approaches the FIXED fine reference,
  NOT as the fine mesh refines under a fixed partition (that limit leaves
  a finite homogenization defect and NO measurable rate). So the ladder is
  a SEQUENCE OF COARSE PARTITIONS on ONE fixed fine problem: a **16-cell**
  ≥2G heterogeneous asymmetric-SigS slab, vacuum→reflective tilt; coarse
  partitions **P ∈ {2, 4, 8}** coarse cells (each a contiguous union of
  fine cells). The IDENTITY partition (P=16) is exact for BOTH weightings
  — `tests/sn/test_homogenization.py:194`
  (`test_identity_homogenization_recovers_per_cell_materials`) proves the
  XS are unchanged when coarse=fine — so the RATE of gap-shrinkage as
  P→16 is what discriminates.
- **Setup / re-solve route (GREENFIELD — no existing test re-solves a
  homogenized mesh; nail every fixture parameter):**
  `fwd = solve_sn(materials, fine16, quad)`;
  `adj = solve_sn_adjoint(materials, fine16, quad)`;
  per P: `mm_adj = fwd.homogenize(coarse_P, adjoint=adj)`,
  `mm_fwd = fwd.homogenize(coarse_P)`;
  `k_adj(P) = solve_sn(mm_adj.materials, mm_adj.mesh, quad).keff`,
  `k_fwd(P) = solve_sn(mm_fwd.materials, mm_fwd.mesh, quad).keff`;
  `k_fine = fwd.keff`.
- **Gate (COMPARATIVE — PRIMARY, per #281 acceptance):**
  1. `gap_adj(P) = |k_adj(P) − k_fine|` is STRICTLY SMALLER than
     `gap_fwd(P) = |k_fwd(P) − k_fine|` on EVERY rung P ∈ {2,4,8}.
  2. `gap_adj` shrinks FASTER: `gap_adj(8)/gap_adj(4) <
     gap_fwd(8)/gap_fwd(4)` (adjoint is higher-order in δφ).
- **Gate (absolute order — SECONDARY/optional):** a log-log fit of
  `gap_adj` vs coarse cell width `h_coarse` with slope ≳ 2 while `gap_fwd`
  ≈ 1, on the 3 rungs. Optional because the coarse-DISCRETIZATION error
  (below) contaminates the absolute slope.
- **HONEST SCOPE (correction a — what C2 does NOT prove).** The coarse
  re-solve carries TWO error sources: (i) the XS-COLLAPSE worth error
  (what adjoint weighting reduces) and (ii) the coarse-mesh spatial
  DISCRETIZATION error (diamond-difference on 2/4/8 cells vs 16) — SHARED
  by both weightings and NOT reduced by adjoint weighting. Since both arms
  solve on the SAME coarse mesh, (ii) is IDENTICAL for both, so the
  COMPARATIVE delta `gap_fwd − gap_adj` cleanly isolates (i). C2 proves
  *adjoint weighting reduces the XS-collapse contribution to the keff
  gap*, NOT *the coarse solve is exact*. This is precisely why comparative
  is primary and absolute-order is secondary.
- **Anti-#5 pairing:** k_fine is the L1-anchored limit (P1.3-certified),
  so C2 anchors the value the gaps shrink toward — "order shrinks" +
  "limit independently correct" (vv: convergence rate necessary, NOT
  sufficient — never O(h²) to a possibly-wrong limit).
- **Mutation (REDs C2 under -O):** feed φ (forward) as the test weight
  instead of φ\* → `gap_adj` degrades to `gap_fwd` → the strict-inequality
  assertion `gap_adj(P) < gap_fwd(P)` collapses to equality → RED.
- **Config-blindness:** an accidentally self-adjoint config gives
  `gap_adj ≈ gap_fwd` (no separation) → REJECT (the C1 normalized-shape
  dud-guard, run on the same fixture, catches it). **Claim:** eigenvalue
  (first-order stationary). **Pillar:** closed-form (k_fine reference).

### C3 — homogenize Mode-11 CAPTURE-and-compare: the sigma-frame weight IS φ\*⊙φ (structural).

The shipped sentinel `tests/sn/test_homogenization.py:311-369`
(`test_homogenize_routes_through_the_petrov_galerkin_frame`) only COUNTS
`WeightedIndicatorBasis.analyze` / `FrameBase.project` calls — it NEVER
captures the weight ARRAY (correction b). C3 is the CAPTURE-and-compare
UPGRADE (vv Mode-11 sharpening — the in-process wrap):
- **Sentinel:** monkeypatch `WeightedIndicatorBasis.__init__` (or the
  sigma-frame construction site, `solution.py:642-644`) to CAPTURE the
  weight array handed to the test basis, then assert `counter > 0` AND
  `np.array_equal(captured_sigma_weight, φ*⊙φ)` — the elementwise product
  in the "ij"/C flat-cell order (`solution.py:641`), NOT bare φ\* and NOT
  bare φ. Fires REGARDLESS of how close φ\*≈φ (a value gate cannot — C1's
  gap → 0 as φ\*→φ, but C3's array identity is exact).
- **Why NOT vacuous (Mode-11):** a bit-identity-preserving regression that
  silently kept the forward weight φ (adjoint dropped) produces a
  DIFFERENT number than C1's φ\*⊙φ reference expects — but on a
  near-self-adjoint fixture C1's value gap is small, so C3 is the
  STRUCTURAL backstop that reds regardless. Pair C1 (value) + C3
  (structural) — the exact split Mode-11 exists for.
- **Mutations (each REDs C3 under -O):**
  1. Wire BARE φ\* as the weight → `array_equal(weight, φ*⊙φ)` RED (the
     `frame.rst:3458` "φ→φ\*" trap — C3 is its committed catcher).
  2. Wire forward φ (adjoint dropped) → RED even when the value barely
     moves.
- **Claim:** structural. **Not a pillar** (call-graph).

### C4 — condense energy-axis discriminator: AVERAGE moves, MARGINALIZE frozen (L0).

The energy-axis analogue of C1, on `Solution.condense` /
`Mixture.condense` (`orpheus/data/macro_xs/mixture.py:313-410`). The
ratified API covers condense (`solution.py:485`). `Mixture.condense` has
TWO morphisms (`mixture.py:388-408`): **average** (`frame.project`,
spectrum test weight — vectors + matrix SOURCE `g_from` axis) and
**marginalize** (bare `@ table`, NO weight — matrix SINK `g_to` axis +
χ). The adjoint enters **ONLY the average morphism**, giving C4 a free
STRUCTURAL discriminator (correction d):
- **Config:** 4-fine-group → 2-coarse-group (reuse
  `tests/sn/test_condensation.py:_EG_FINE`/`_EG_COARSE`), asymmetric SigS,
  WITHIN-GROUP-VARYING representative spectrum (else the average morphism
  does not discriminate), strongly non-self-adjoint so φ\*_spectrum
  differs in shape from φ_spectrum.
- **Gate A (average morphism MOVES):** adjoint-condensed vectors (SigT, …)
  and matrix SOURCE axis ≠ forward-condensed, and == the hand reference
  with the **B1-derived** average weight (vector spectrum bilinear =
  φ\*_spectrum⊙φ_spectrum; matrix source axis per the B1 rule).
- **Gate B (marginalize morphism FROZEN — the sharp discriminator):** the
  MARGINALIZE-collapsed channels (χ, matrix SINK axis) are BIT-IDENTICAL
  between `adjoint=None` and `adjoint=φ*`:
  `np.array_equal(condense(adjoint=adj).chi, condense().chi)` and the sink
  axis. A regression that leaks φ\* into the χ/sink collapse REDs Gate B.
- **Dud-guard:** normalized-shape guard on the spectra
  (`φ*_spectrum/‖·‖ vs φ_spectrum/‖·‖`), as C1.
- **Tolerance:** rtol 1e-12 (Gate A hand loop); `array_equal` (Gate B).
- **Mutations (each REDs C4 under -O):**
  1. Adjoint leaks into χ / the sink axis (marginalize) → Gate B
     `array_equal` RED.
  2. Adjoint dropped from average (uses spectrum φ) → Gate A ≠ hand
     reference → RED.
- **Claim:** model (Gate A) + structural (Gate B). **Pillar:** closed-form
  (hand loop).

### C5 — condense Mode-11 CAPTURE sentinel: the average-frame spectrum weight IS φ\*_spectrum⊙φ_spectrum (structural).

The condense sibling of C3 (correction b — "add the condense sibling
sentinel"). Capture the spectrum weight handed to the average-frame's
`WeightedIndicatorBasis` (`mixture.py:385`) and assert
`np.array_equal(captured_spectrum_weight, φ*_spectrum⊙φ_spectrum)` — NOT
bare φ\*_spectrum, NOT bare φ_spectrum. Assert counter>0.
- **Mutation (REDs under -O):** wire bare φ\*_spectrum (the `frame.rst:3458`
  trap on the energy axis) or forward φ_spectrum → RED.
- **Claim:** structural. (The condense `adjoint=None` degenerate pin lives
  in §4.0.)

### Cχ — χ-simplex positive control under the B1-derived emission rule (foundation).

Constructing the adjoint-weighted `Mixture` VALIDATES χ via
`EmissionSpectrum` (`orpheus/data/emission_spectrum.py:96,101` — `raise
ValueError` on a non-simplex, so it FIRES under -O; NOT a bare `assert`,
Mode-8-safe). The gate is a POSITIVE control (correction f; vv anti-#11 —
a validation invariant needs a MUST-NOT-RAISE positive test, not only a
negative one): the χ_R produced under the B1-derived fission-dyad rule
MUST still construct a valid Mixture (no raise) AND pass
`assert_normalized`. A convex combination of simplices is a simplex, so IF
the B1 rule keeps χ_R a convex average of the fine χ_i it holds by
construction — this gate pins that the B1 rule did NOT break convexity.
- **Config:** the C1 non-self-adjoint fixture (χ=[1,0]-type fast-peaked
  fine spectra), for BOTH the homogenize χ collapse and the condense χ
  (marginalize) collapse.
- **Mutation (REDs under -O):** a B1 χ-rule that admits a NEGATIVE φ\*
  weight (importance can be signed for a non-fundamental mode) → χ_R leaves
  the simplex → `EmissionSpectrum` RAISES → the positive control REDs
  (confirms the rule guards sign).
- **Claim:** foundation (intrinsic simplex law). **Pillar:** closed-form
  (the simplex constraint).

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
                    §4.0 degenerate pins (adjoint=None ≡ shipped, bit-identical)
                                 │
                 ┌───────────────┴────────────────┐
        homogenize axis                     condense axis
                 │                                │
     C1 (w=φ*⊙φ ≠ φ, L0) ──┐           C4 (average moves / marginalize frozen, L0)
                 │          │                     │
     C2 (keff-gap comparative order, L2)   C5 (Mode-11 CAPTURE: spectrum weight = φ*_s⊙φ_s)
                 │          │                     │
     C3 (Mode-11 CAPTURE: sigma weight = φ*⊙φ) ───┤
                 └───────────────┬────────────────┘
                                 ↓
                    Cχ (χ-simplex positive control, BOTH axes)
```

- **Before φ\* can be called correct:** P1.0 → P1.5 GREEN + every mutation
  reddens its gate under `-O`.
- **Before the consumer is built:** φ\* certified (the whole P1 chain).
- **§4.0 gates the axes:** the `adjoint=None` bit-identity pins run FIRST
  (they need no φ\*); C1–C3 (homogenize) and C4–C5 (condense) are
  independent axes; Cχ spans both (the emission-rule simplex control).
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
| **Adjoint-weighted homogenization drops φ\*** (uses forward φ) | C1 value (φ\*⊙φ hand ref) + C3 Mode-11 CAPTURE (pair: value catches it on non-self-adjoint, capture catches it structurally regardless of φ\*≈φ) |
| **Adjoint weight is BARE φ\* not the bilinear product φ\*⊙φ** (the `frame.rst:3458` "φ→φ\*" doc trap) | C3/C5 CAPTURE-and-compare: `array_equal(weight, φ*⊙φ)` REDs on bare φ\*; C1 value REDs vs the φ\*⊙φ hand loop. **The most likely P6 implementation bug** — the doc language invites it |
| **Adjoint leaks into the condense MARGINALIZE morphism** (χ / matrix sink axis) | C4 Gate B: `array_equal(condense(adjoint=adj).chi, condense().chi)` — the marginalize channels MUST stay byte-frozen (adjoint touches only `average`) |
| **`adjoint=None` arm is a re-derived look-alike, not the shipped path** | §4.0 degenerate pin: 0-ULP `array_equal(homogenize(), homogenize(adjoint=None))` REDs on any FP drift (a ones-weight PG route drifts at ULP) |
| **B1 χ-rule breaks the simplex** (negative φ\* weight on a non-fundamental mode) | Cχ positive control: adjoint-weighted `Mixture` construction RAISES `ValueError` via `EmissionSpectrum` (fires under -O) |
| **Composite S.H wrong even with correct leaf transposes** (frame-consistency, ERR-061 family) | P1.0b dense `Gᵀ S Gᵀ` oracle (catches a wrong composition that per-leaf gates pass) |

→ propose as `numerical-bug-signatures` Signature 11 (adjoint role-swap
invisibility under symmetric kernels) once a real bug is caught in P6
implementation; the BARE-φ\* vs φ\*⊙φ weight confusion (the doc-trap row
above) is a strong Signature-11-sibling candidate — a plausible-wrong
adjoint-weight that passes a COUNT sentinel and a near-self-adjoint value
gate but fails the CAPTURE gate. Until a real bug lands, skill-table rows
above.

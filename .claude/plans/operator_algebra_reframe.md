# operator_algebra.rst reframe — intrinsic types + full restructure

**Status:** PLAN (awaiting user sign-off on the split architecture). Part of the SN-book
re-architecture (Task #3) / doc-corpus campaign (#231). Branch `docs/sn-doc-architecture`.

**Quality bar:** [[feedback-articulation-lossless-disassembly]] — organize = eliminate noise,
retain signal; lead with the intrinsic types; every derivation expressed (manual-leverage here —
see §Algebra-of-record). **Scope ruling (user, 2026-07-19): FULL RESTRUCTURE** — reframe + relocate
dev-history + split/retire the advanced deep-dives to their own pages.

---

## 1. Why

`foundations/operator_algebra.rst` is 9,196 lines — a whole book crammed into one foundations
page. It still calls itself a "Phase 0 stub," its Key Facts assert the **retired** algebra
`(L+C−S−F−B)ψ=q` on a stale 3-block state, and the intrinsic-type content §2.2 wants (C =
multiplication, F/S = kernels) is buried at lines 1600/2174 under "§5.6/§5.7 grand-report"
framing, wrapped in ~7,000 lines of Wave-labeled campaign-log + advanced deep-dives.

The reframe: **lead with the intrinsic-types spine + the honest algebra; split the deep-dives to
their own pages; relocate the campaign-log to a Development-history changelog.**

## 2. Code ground truth (explorer, 2026-07-19 — authoritative over the stale page)

| Sym | Class (file) | Intrinsic type |
|---|---|---|
| **C** | `MultiplicationOperator` (`transport/operators/multiplication_operator.py`) | diagonal multiplication `M[σ_t]`; `*`-homomorphism laws; `CollisionOperator` subclass retired (#261) |
| **F** | `FissionOperator` (`transport/operators/fission.py`) | rank-1 dyad `|χ⟩⟨νΣ_f|` via `outer(chi, ReactionRateFunctional(νΣ_f))`; no `solve` (singular); adjoint = dual dyad |
| **S** | `ScatteringOperator` (`transport/operators/scattering.py`) | **nonlocal kernel `R∘Λ∘M`** (Funk–Hecke) — **NOT a projection**; projection is only the analysis face `M` |
| **L** | `StreamingOperator` (`sn/operators/streaming.py`) | σ-free `Ω·∇ψ` leaf; not invertible alone |
| **B** | `SNBoundaryOperator` (`sn/operators/boundary.py`) | boundary trace reflection `R·G` |

- **Honest algebra = `A = L + C − S − B`**, F on the RHS (`Aψ=(1/k)Fψ`), 2-block state
  (`bulk⊕boundary`). Assembled at `coupled_system.py::build_within_group_system` (`A_AA = LC − S − B`).
- **`A_loss = L + C` = `InvertibleOperator`** whose `.inverse()` **is** the `SweepOperator` —
  "the algebraic foundation of the entire SN method." matvec/adjoint/sweep = 3 actions of ONE op.
- **SN extends S for anisotropy** via the Λ factor `LegendreMomentScattering` = `Σ_ℓ P_ℓ⊗Σ_{s,ℓ}`,
  retained order = `scattering_order` (0 = P0 isotropic; ≥1 = anisotropic).
- **Operator surface:** `LinearOperator[Domain,Codomain]` base (NOT forced endomorphism); NO
  capability registry — three recursive predicates `is_invertible`/`is_adjointable`/`is_assemblable`;
  verbs apply/solve/inverse/`.H` (metric-weighted); composition `+`/`−`/`@`/scalar-`*`.
- **⚠ In-code drift to fix in tandem:** `transport/full_field.py:9` module docstring still carries
  the retired `(L+C−S−F−B)ψ=q` spelling.
- **⚠ No algebra-of-record** for the operator composition (grep of `derivations/` empty). The page
  is manually-leveraged from the (rich) production docstrings. Closest symbolic source =
  `derivations/common/transport_equation.py` (continuous Boltzmann PDE only). → follow-up issue.

## 3. Blast radius (measured 2026-07-19)

47 std-labels, 78 eq-labels. Inbound: 30 `:doc:` links; labels referenced from other docs incl.
`operator-algebra`→8, `wavefront-flux-cochain`→4, `g-adjoint`/`green-operator`/`integral-kernel-category`/
`inverse-application-driver`/`matrix-inverse-operator`/`tensorial-framing`→3 each; code docstrings ref
`operator-algebra`(×3)/`eigenvalue-posing`(×1)/`tensorial-framing`(×1). **Rule: a split carries its
labels with the content** (labels are path-immune — inbound `:ref:`/`:eq:` survive the page move);
`:doc:` links to a relocated section re-point to the new page. Gate every step `-E -W`.

## 4. Section disposition (the 26 H1 sections)

**STAY — the reframed foundation spine (the intrinsic types + core algebra + posing):**
- Intro (rewrite from "Phase 0 stub") · Key Facts (reconcile → honest algebra)
- The operator surface (merge Definitions + "three-layer surface" + "apply-linear/solve-not" + "Composition algebra")
- The type partition (Operator/Kernel/Functional, locality)
- **C = the multiplication operator** (§ line 1600) · The functional category (§1818)
- **F + S = the integral-kernel category** (§2174) — F dyad, S = `R∘Λ∘M` (projection located at `M`)
- The composite `A=L+C−S−B` + the invertible sub-composite → sweep (from "Pure-L streaming")
- How methods extend (SN anisotropy via Λ; forward-link SN book)
- Tensor-product algebra (§3295) · Diagonal-operator-on-a-tagged-axis (§1523) — supporting the kernels/C
- Eigenvalue posing + power-iteration (§8528) · Trace spaces Γ∓ (§8447) — *candidates to split, see fork*

**SPLIT to own pages — CROSS-METHOD deep-dives (foundation/):**
- `operator_inverse_family.rst` ← #226 taxonomy: solver-builds-inverse (§6438) + Green operator
  (§6762) + materialising functor / matrix-inverse (§7152) + assembly axis (§7592). ~1,300 ln.
  Labels: `inverse-application-driver`, `green-operator`, `matrix-inverse-operator`, `operator-algebra-assembly-axis`.
- `operator_tensor_network.rst` ← Tensor-Network Decomposition Wave T (§3810). ~970 ln.
  Labels: `wave-t-tensor-network`, `wave-t-shape-table`, `wave-t-orchestrated-apply`.
- G-adjoint (§5305) → **own page `foundations/operator_adjoint.rst`** (REVISED 2026-07-19 from the
  frame.rst fold). Reason: the G-adjoint is the OPERATOR's composite metric adjoint
  `op.H = G⁻¹AᵀG` over `FullFieldSpace` (block-diagonal metric, singular-trace pseudo-inverse) —
  operator-algebra content. `frame.rst`'s adjoint is the `R.H` / Petrov-Galerkin **test-space**
  adjoint — a different level. Cross-link both. Label `g-adjoint` (×3 docs) travels to the new page.

**SPLIT — SN-LEANING deep-dives → INTERIM FOUNDATION PAGES (user ruling 2026-07-19: decouple
from the SN-book campaign; SN book can absorb later if genuinely SN-only):**
- Affine-typed field algebra (§4783, ~520 ln) → `foundations/field_algebra.rst`. Label `affine-typed-residual`.
- Interior face-flux cochain C¹_int (§5886, ~550 ln) → `foundations/wavefront_cochain.rst`. Labels `wavefront-flux-cochain`(×4), `wavefront-cochain-biproduct`.
- Coupled block operator ψ½ System-B (§7778, ~670 ln) → `foundations/coupled_block_operator.rst`. Label `coupled-block-operator`.

**RELOCATE to a Development-history changelog (page bottom):**
- The existing "Development history" (§9080) + the Wave-labeled rationale threaded through Key Facts
  and keeper sections. "What was tried and failed" is KEPT (Cardinal Rule 3) — distilled, not deleted.

**EVALUATE for retirement (3-search audit each):**
- When-a-moment-earns-a-type #263 (§3537) → likely fold to the type-vs-property rule / changelog.
- Boundary-conditions-as-Wave-0/1 primitives (§3644) → likely relocate to `foundations/boundary_conditions.rst`.
- Cross-solver-consumption forward-ref (§1483) → fold into the composite/consumption section.
- Orphan-tracking notes (SumOfTensorProductsOperator #260, retired Wave-T split #238) → distill to changelog/issues.

## 5. Phasing (each phase = one or more gated `-E -W` commits)

- **P1 — Spine reframe IN PLACE (no splits).** Rewrite intro + Key Facts → honest algebra; add the
  intrinsic-types lead framing (type partition → C → functional → F/S); reframe headings from
  "§5.x" to intrinsic types; fix `full_field.py:9`. Highest value, self-contained. ⏸ COMPACTION POINT.
- **P2 — Split cross-method deep-dives** (inverse-family, tensor-network, g-adjoint→frame). One
  commit per page: extract (carry labels + machine header + toctree), replace with a `:doc:` pointer,
  verify inbound refs survive `-E -W`. ⏸ COMPACTION POINT after each.
- **P3 — SN-specific deep-dives** per the fork ruling (SN book vs interim page).
- **P4 — Dev-history → changelog; retire/relocate the evaluate-set** (3-search audits).
- **P5 — Final polish; file the algebra-of-record issue; update memory router + SN "Development history".**

## 5a. STATUS

- **✅ REFRAME COMPLETE (2026-07-20). P1–P5b all landed + committed on `docs/sn-doc-architecture`.**
  `operator_algebra.rst` **9223 → 4599 ln (−50%)**. The page now LEADS with the honest algebra
  `A = L+C−S−B` (F on the RHS) + the intrinsic-type partition (Operator / Kernel / Functional, surfaced
  before the first type section); the 6 cross-method + SN-leaning deep-dives live on their own
  foundation pages; the campaign-log is a page-bottom Development-history changelog (10 rows); Key Facts
  is provenance-free. **User ruling (2026-07-20): streaming stays EARLY** — it is apply/solve's
  foundation, so relocating it late (the planned "Move A") would invert premise↔conclusion, and the
  page leads with the types regardless. Move A abandoned + Move B (surface preamble) skipped; **P5b =
  the additive `intrinsic-operator-types` preview H1 before the C block ONLY** (commit `2d2d7595`).
  Tail reorder (eigenvalue/trace above BC) DEFERRED as out-of-scope. Follow-ups: **#301** (algebra-of-
  record, FILED); `_flux_role.py:115` code message (§6); optional T2/#263 design-principles home. This
  plan is now ARCHAEOLOGY — the durable record is the doc pages + the archivist L-026..L-029 calibrations.

- **⏸ COMPACTION POINT — P5a COMPLETE, P5b (SPINE REORDER) NEXT (2026-07-19).** Tree clean, HEAD
  builds `-E -W` clean. `operator_algebra.rst` now **4573 ln**. **NEXT = P5b: the spine reorder** —
  physically lead the page with the intrinsic types. Two CORE moves: (1) MERGE the four scattered
  surface sections {`Definitions`, "apply is linear / solve is not", "three-layer operator surface",
  "Composition algebra"} into one **"The operator surface"** H1; (2) PULL the composite/sweep section
  ("Pure-L streaming + the affine collision split", currently too EARLY at ~L528, interleaved between
  Definitions and apply/solve) DOWN to after the integral-kernel (F/S) section — it is the ASSEMBLY of
  the typed operators, belongs after C/S/F are introduced. Target body order: KeyFacts → [Operator
  surface] → [type partition framing] → C (diagonal/multiplication) → Functional → F/S (integral-kernel)
  → [composite A=L+C−S−B + sweep] → supporting (tensor-product) → moment-type #263 → eigenvalue/trace →
  BC + 6 deep-dive stubs → Development-history changelog. **Dispatch a FRESH archivist** (the P3–P5a
  instance is ~900K tokens deep). Optional tail: the T2/#263 relocate-to-a-design-principles-home
  question; the #301 algebra-of-record is FILED.
- **P5b SCOPED (fresh archivist `a35ca096`, 2026-07-19) — core decision ESCALATED to the user.**
  Close-read findings: (1) **Move A's hazard** — relocating "Pure-L streaming + the affine collision
  split" (522–717) after the types creates a **pedagogical inversion**: the EARLY apply/solve section's
  argument ("invert the SUM, not the leaves") rests on the streaming affine split
  (`streaming-action-pure-l`/`-cell-balance`); moving streaming late makes apply/solve forward-ref its
  own foundation (4 `:eq:` + prose L983 "the previous subsection" → wrong). Survivable but a real
  premise↔conclusion flip → **user taste call** (keep-early recommended; the page already leads with
  C/S/F after the surface regardless). (2) Re-title "The composite A=L+C−S−B" **over-claims** — the
  section assembles only the `L+C` LOSS sub-composite (S,B are source-iteration; the sweep/resolvent
  lives in apply/solve). Use "Assembling the loss operator: the streaming leaf L and the L+C composite".
  (3) **Hidden siblings:** C = TWO H1s (`diagonal-operator` 1529 + `multiplication-operator-promotion`
  1606 — move as a unit); "Trace spaces Γ∓" (3747) is REAL content, not a stub (it is the "trace" in
  the tail). (4) **Move B** = LIGHT grouping (a short "The operator surface" preamble H1 after Key
  Facts) — reject the full H1→H2 heading-cascade merge (`-E -W` risk). (5) Type-partition framing
  already exists (Functional §5.6 suffix-law + codomain table) but appears mid-stream; optional light
  preview H1 before the C block (preview + `:ref:`, no table copy). (6) **Ambiguity C — tail reorder**
  (eigenvalue/trace above BC) is a THIRD move, NOT P5b scope → **DEFER, leave tail as-is.** Full
  proposal in the archivist transcript. Blast radius: 2 moving eq-labels path-immune; only real prose
  break = L983. **HOLD all P5b execution until the streaming-positioning ruling.**
- **P5a DONE** (Key-Facts↔changelog reconciliation; committed this session) — 10 Wave-metadata strips
  from the high-traffic Key Facts (hashes / round-labels / branches / dates → relocated to the
  changelog; the invariants + the 2 eq-labels + their vv-status + the active gotchas KEPT) + **6
  changelog rows** added (the 5 scoped milestones + a 6th for the coupled-block / walk-unification
  #280/#282 @ `6732778a`, all git-verified as HEAD ancestors). Also fixed the affine-space straggler
  `A`→`\mathbb{A}` in Key Facts (split #1's rename had not reached the summary bullet; the 7
  operator-`A` sites left intact). `-E -W` clean.
- **P4 dispositions (adjudicated + verified 2026-07-19):** **T1** cross-solver — ¶2 (live
  `_as_scipy_linop` SSOT + #2-DSA "generalise-don't-resurrect" guidance) FOLDED into Composition-
  algebra, reframe-fixed (`loss_minus_gains = (L+C)ψ − Σ g ψ = Aψ`); ¶1 forward-ref DELETED (distilled
  to one clause). **T2** moment-type #263 — KEPT in place (durable invariant; legit floor-rule/ceiling-
  theory split). **T3** BC-primitives — **DISPOSITION INVERTED by the close read: relocate+merge →
  DE-DUP** (Cardinal Rule 2). The section was near-total duplication of `boundary_conditions.rst`
  (§1797 realization-map G_α table is RICHER; §2424 descriptor-tree algebra covers the rank-N eq /
  Marshak / two-families / Wave-11-β1-β2 history in more detail; T3's `R=Σ G_α⊗A_α, A_α∈[0,1]` is just
  the scalar `Σ c_α G_α` in tensor dress — nothing lost). → 24-ln `:doc:` pointer KEEPING the bridge
  (B a first-class sibling in A=L+C−S−B); 2 labels re-anchored onto canonical content, eq-label
  `bc-rank-n-as-sum-of-products` dropped (duplicated). **T4** SumOfTensorProducts — KEPT the live-type
  H2 + rationale; #261 dropped from the deferred-list (verified LANDED @ `bbe8a51d`), #260 kept open.
  All three verified in review (T3 de-dup checked against both boundary_conditions sections).
- **NEW DISCIPLINE for every P4/P5 move (discovered P3 split 3; now archivist L-026):** a split's
  inbound `:ref:`/`:eq:` auto-repoint (path-immune), BUT **consuming PROSE that names the old host
  page** ("see `:ref:`X`` in `:doc:`operator_algebra``") silently misdirects the reader and is
  **`-W`-blind** (the link resolves; only the prose lies). Every move MUST run a whitespace-flattened
  consuming-prose sweep (the `prose_sweep2.py` pattern in `$CLAUDE_JOB_DIR/tmp`) and repoint each
  `:doc:` to the label's real home. P2's inverse-family/adjoint splits were swept retroactively here.
- **P3 split 3/3 DONE @ 7ae71f60** — `coupled_block_operator.rst` (709 ln; the ψ½ ray as System B of
  the 2×2 block-triangular operator; 10 labels). Clean contiguous relocation, NO reconciliation
  (block notation `A_AA/A_AB/A_BA/A_BB` reframe-consistent — the diagonal `A_AA` IS `L+C−S−B`);
  trace-spaces stays on the spine. source 5315→4657.
- **P3 split 2/3 DONE @ f45527d4** (+ chore `d190d47c`) — `wavefront_cochain.rst` (592 ln; the
  interior face-flux cochain C¹_int; 5 labels). Pure mechanical, no reconciliation. Boundary catch:
  the anchorless "inverse family" pointer stub (P2's) sat in the briefed range; extracted only
  3875–4426, the stub stays. Carry-forward (`V_inflow⊕V_outflow`) pre-verified legit face-level.
  source 5856→5315.
- **Consuming-prose hygiene @ 31829fe0** — tree-wide repoint of prose naming the old host page for
  all 6 moved sections: g-adjoint→`operator_adjoint` (boundary_conditions); green-operator /
  matrix-inverse-operator / inverse-application-driver → `operator_inverse_family` (sn/history,
  sn/index). STAY-spine refs (`capability-set-semantics`, `field-type-vs-property-criterion`) left.
- **P3 split 1/3 DONE @ 84dcb512** (+ memory chore `22c36de0`) — `field_algebra.rst` (581 ln; the
  affine torsor + state/displacement/residual/source role grid). 6 labels carried (path-immune);
  source 6365→5856. RECONCILIATION applied (clean Option-Y, escape hatch did NOT fire — the torsor
  rests on flux-state geometry, operator-independent): operator `A = L+C−S−B` page-wide; the sweep
  demoted to `(L+C)⁻¹` = inner kernel of `A⁻¹` (never `A⁻¹`=the sweep); SI-increment
  `M = (L+C)⁻¹(S+B)`; affine SPACE renamed off bare `A` → `\mathbb{A}` (+ head disambiguation note).
  BOUNDARY CATCH: the brief's range (3844–4383) OVERSHOT — an anchorless "composite metric adjoint"
  pointer stub (P2's link to `operator_adjoint`) sat inside it; extracted only 3844–4372, the stub
  STAYS on the spine (operator's adjoint face, not field content). Carry-forward #2
  (`V_inflow⊕V_outflow`) pre-verified legit face-level (no change). `-E -W` clean, 0 warnings.
- **P2 COMPLETE (3/3 split) @ 69abc073** — `operator_adjoint.rst` extracted (the composite
  metric adjoint `op.H = G⁻¹AᵀG` over FullFieldSpace; 581 ln; `g-adjoint` + 6 eq-labels carried;
  source 6928→6365). Own page (revised from frame.rst fold); intro carries an `.. important::`
  drawing the operator-adjoint vs frame `R.H` test-space-adjoint distinction; cross-links frame.rst.
  The 2× `(L+C−S−F−B)†` confirmed LEGIT (full-operator "intentionally unreachable" + the honest-algebra
  "never fused" affirmation) — left verbatim, matches the carry-forward. Also fixed a pre-existing
  dangling ref ("reachability table below" → the real Supersession note). `-E -W` clean.
  **Page arc: ~9223 → 6365 ln (−31%); 3 new sub-pages (inverse-family 1388 / tensor-network 1023 /
  adjoint 630). NEXT = P3 SN-leaning splits (field_algebra, wavefront_cochain, coupled_block_operator).**
- **P2 split 2/3 DONE @ 077ed7bc** — `operator_tensor_network.rst` extracted (Wave-T
  tensor-network decomposition; 973 ln; 10 labels — 5 anchors + 5 eq — carried; source 7891→6928).
  Pointer + toctree/list-table wired. `-E -W` clean. 1 L35 fix (a cross-page "above"→`:doc:`);
  the S\ :sub:`N` inline-role title mirrors the SN book H1 convention. **NEXT P2 split 3/3 =
  `operator_adjoint.rst` (g-adjoint, OWN page — see §4 revision).**
- **P2 split 1/3 DONE @ fa3afdf4** — `operator_inverse_family.rst` extracted (#226
  taxonomy: driver-applied inverse / Green / dense matrix-inverse / assembly axis; ~1,340 ln;
  6 anchors + 6 eq-labels carried; source 9223→7891). Pointer stub + toctree/list-table wired.
  `-E -W` clean. Archivist L-026: an f-string mangled `:math:`A^{-1}`` → `A^-1` in authored prose
  (`-W`-blind, valid LaTeX) — spot-check authored math after any programmatic templating.
  **NEXT P2 split 2/3 = `operator_tensor_network.rst` (Wave T).**
- **P1 DONE @ 167ac25b** — honest-algebra reconciliation: intro rewritten (intrinsic-types
  lead + `A=L+C−S−B`, F on RHS, `L+C`=invertible sub-composite→sweep, `A`=loss operator); the stale
  3-block `(L_full+C−S−F−B)ψ=q` Key Fact → honest 2-block; the `operator-fixed-source`/`-eigenvalue`
  equations → `Aψ=q` / `Aψ=(1/k)Fψ`; `full_field.py:9` docstring fixed. `-E -W` clean.
- **Reconciliation carry-forward (verify during the owning split phase, NOT stale-blanket):**
  - α-eigenvalue `(L+C−S−F−B)ψ=−αTψ` (~L8800/8849, eigenvalue-posing §) — **CORRECT** (F on LHS for
    time modes); leave.
  - G-adjoint `(L+C−S−F−B)^†` (~L5715/5729) — review at the g-adjoint→frame split (P2).
  - Cochain `V_inflow⊕V_outflow` (~L6063) — verify it is the FACE-level boundary cochain (legit), not
    the retired top-level 3-block state, at the cochain split (P3).
  - affine-field-algebra `A = L+C` naming (~L4832/4867) — clashes with the reframed `A=L+C−S−B`;
    reconcile at the field_algebra split (P3).
  - 40× `A_{\rm loss}` — a valid synonym for `A` (both = the loss operator `L+C−S−B`); no rename needed.

## 6. Follow-ups
- **Issue FILED @ #301** (2026-07-19): build an operator-algebra algebra-of-record (SymPy for
  C/S/F/L + composition) so the page becomes generator-expressed, not manually-leveraged. Carries
  the `_flux_role.py:115` minor code follow-up as a note.
- The `eigenvalue-posing` home may reconcile with the path-integral root's eigenvalue anchor (Phase H).
- **Code follow-up (out of scope this doc campaign; flagged by archivist during P3 split #1):**
  `orpheus/transport/fields/_flux_role.py:115` — `FluxRole.__add__`'s `TypeError` message reads
  "rate-density L-actions (sources / residuals)", using the pre-reframe `L` symbol. Rate-density is
  the codomain of the full loss operator `A = L+C−S−B`, not of `L` (streaming) — the message is a
  candidate for the same `A`-reframe. One-line wording fix; ride it with a future SN-message pass.

## 7. Compaction protocol
Commit → checkpoint this plan (STATUS + last commit) → /compact re-anchors from this file + git log,
never the summary alone. [[feedback-compaction-points-in-campaign-plans]]

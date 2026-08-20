# CS3-R — clear-context adversarial review of the cone carve: findings memo

Reviewer: main agent, post-compaction (2026-08-19). Charter:
`.claude/plans/space_and_kernel_binding_campaign.md` §4-R. Method: Phase 1
harsh → Phase 2 verdicts, every withdrawal with its structural reason.
Written incrementally; the completeness denominators are at the end.

## Instruments run

- Sphinx rebuilt at HEAD `f43758d8` (clean, 0 warnings) → graph fresh; the
  SessionStart hook's 2 dead references were the STALE graph (`a4d416a4`,
  mid-carve) — `[M]` `dead_references` at HEAD = **0 dead / 56 checked**.
- Smell sweeps: `twin_paths` (20) / `dead_functions` (30) / `native_place`
  (20) / `discriminations` — **no CS3-surface findings**; all rows are
  pre-existing (test mesh-builder fixture twins, `tools/research._throttle`
  family, thermal-hydraulics solver twins, derivations diagnostics). Noted,
  out of CS3-R scope (sharpening-order: file-don't-fix territory).
- Independent `qa` census (second grep vocabulary, guarding my filter
  blindness): `scratch/cs3r_census_qa.md` — folded in below when returned.

## Phase 1 findings — COMPLETENESS (torsor machinery surviving as live content)

⛔ **The step-3 claim "all production prose swept to dated past tense" is
REFUTED.** `[M]` `grep -rn "displacement" orpheus/ --include="*.py"` minus
dated-history spellings → a surviving family the carve's own sweeps missed.
The pattern in almost every case: the file WAS touched, and the sweep fixed
the pointed-at line while sibling lines in the same file survived (the
`_coefficient_role` half-fix shape, recurring).

| # | site | defect | verdict |
|---|---|---|---|
| C1 | `tests/numerics/test_si_diagnostic_trajectory.py:9-10` | two `:meth:`…Displacement.*`` roles on the DELETED type (the `…` prefix hides them from `dead_references` — instrument-invisible, concept-real; `c634919f` fixed bullet 1 of 3) | FIX inline |
| C2 | `orpheus/transport/radial_characteristic_field.py:28,101,114,120,126` | present-tense `⊖`-minted displacement in module + class docstrings, a guard comment, and two error messages (messages unpinned — `[M]` no test matches "displacement leaf") | FIX inline |
| C3 | `orpheus/transport/source_sinks/radial_characteristic_{interior,boundary}_source_sink.py:27,30` | three-role enumeration "source, flux, and displacement" as the class gate's present-tense job | FIX inline |
| C4 | `orpheus/sn/operators/boundary.py:1099` + `orpheus/sn/operators/radial_characteristic.py:486` | guard rationales name a "displacement-role member / flux-displacement composite" that is unspellable at HEAD | FIX inline |
| C5 | `orpheus/numerics/moment_layout.py:74` | "the field/displacement shape builders" — the displacement shape builders are deleted | FIX inline |
| C6 | `orpheus/numerics/coupled_system.py:181` | ⭐ present-tense-FALSE **design rationale**: "the member family's AFFINE-TORSOR signatures (flux − flux → displacement …) fit neither `Vector`'s `Self + Self → Self` spelling" — at HEAD the members fit it exactly; the file's own "collapse trigger" (converge the three member-contract concepts onto a named ravellable Protocol) is thereby UNBLOCKED | FIX prose inline + FILE the collapse |
| C7 | `orpheus/numerics/iteration.py:435` ("today's displacement-typed increment" — "today" moved), `:568`, `:766`; `orpheus/sn/solver.py:769,1150(msg),3128` | word-level vestiges; mechanism claims still true, noun stale | FIX inline |
| C8 | `orpheus/sn/acceleration/dsa.py` (12 sites incl. `apply(self, displacement:)` param) + `orpheus/derivations/discrete/sn/dsa.py` (4 prose sites) | the file's domain noun is still "displacement" while its own CS3 admission test already says "full-angular **increment**"; `[M]` no keyword callers, message pin "moment-windowed" preserved | FIX inline (rename → increment) |
| C9 | `tests/sn/mesh/test_radial_characteristic_split_leaves.py:219-222` | "Unlike the flux leaves (affine points, `⊖` mints displacements)" — present-tense dead contrast | FIX inline |
| C10 | `tests/sn/solve/test_affine_carve_bit_identity.py:1-12` | wall header narrates the #208 affine pieces as current AND claims "the SI stopping norm STAYS the flat `_l2_norm(displacement…)`" — doubly stale (type gone; the stop is the rhs-residual since #340 step 5). Prose-only fix; frozen values untouched | FIX inline |
| C11 | `orpheus/transport/spatial/scheme.py:528` | "`is_positivity_preserving` gates negative-flux diagnostics" — aspirational-as-fact (0 production readers; #390 tracks) | FIX inline (tense + cite #390) |
| C12 | five operator gates — `test_collision_operator.py:209`, `test_scattering_operator.py:461`, `test_streaming_operator.py:300`, `test_operators_apply_typed.py:333`, `test_phase_c_gates.py:258` | shared present-tense-false premise: "a general α·ψ₁ + β·ψ₂ (α+β≠1) **is illegal** on affine flux STATES"; assertion body uses the blend form `op(ψ₁ + λ(ψ₂−ψ₁))`, which an AFFINE op satisfies identically (affine maps preserve affine combinations) — only the homogeneity leg in the same test makes the pair discriminating | FIX inline: premise + re-spell to the direct law (see E1) |

**Checked and CLEAN (first-class results):**

- `orpheus/transport/__init__.py:31`, `coupled_system.py:110`,
  `angular_flux.py:121`, both RC flux leaves, `test_flux_vector_algebra.py:20`,
  `test_radial_characteristic_split_leaves.py:15` — dated past tense, correct.
- `docs/theory/methods/sn/curvilinear_one_group.rst:4166-4173` — the
  displacement row is a well-formed ⛔ RETIRED record.
- `field_algebra.rst:1215,1325` — `\ominus` only inside the dated
  overturned-design section.
- `orpheus/fuel/solver.py:612` — MECHANICAL displacement (fuel-clad contact
  gap); different domain, legitimately keeps the word.
- `tests/sn/operators/test_declared_law_is_linear.py` full header (lead 7):
  coherent — the `[M]` battery table carries its own "(ran against the
  pre-CS3 spelling; the direct rows DO red — re-measured at the carve)"
  re-scope. WITHDRAWN as a finding.
- `affine-typed-residual` section + label (lead 2, label half): the kept
  `affine-` name is deliberate, annotated, and 8-citer-load-bearing; the
  claim was never affine. WITHDRAWN.

## Phase 1 findings — EXPRESSIVENESS (leads 1–8, re-derived)

| lead | Phase-2 verdict | structural reason |
|---|---|---|
| 1 (`where_largest`/`cone_violations` twin) | PARTIAL — extract the shared concept | The selection logics differ structurally (top-k by magnitude vs predicate filter); the genuinely shared concept is the **flat-position → values-layout index-tuple decode** (2 identical lines each). Extract `_index_tuples`; do NOT merge the selections (a parameterized magnitude-map would be a boolean-ish mode knob — anti-pattern). |
| 2 (residual expressiveness) | WITHDRAWN (no code change now) | The stop's flat `_l2_norm` and the diagnostic's space-norm are two DOCUMENTED conventions; the principled successor of both is CS2's metric relocation (ruled: CS2 owns re-deriving the frozen numbers). The loop's `rhs_prev − rhs` free-identity residual is layer-correct: it never escapes the expression, so the `from_balance` role-transition doctrine (which governs exposed surfaces) does not bind it; routing through `evaluate_residual` would cost a matvec per pass for no claim gain. |
| 3 (cone vocabulary at consumers) | WITHDRAWN (docs already coherent; one tense fix = C11) | The three sites are three different CLAIM KINDS — element predicate (`cone_violations`), law admission (realizer's 𝒜=−1 refusal), positive-scalar ray normalization (`power_iteration`) — sharing ONE documented definition: the cone chapter names all three (`[M]` field_algebra.rst:113,327,412,645,696,1148). A shared code surface would be ceremony across claim kinds. #390 already tracks the flag-reader wiring. |
| 4 (`CoefficientRole` empty marker) | KEEP, no edit | The docstring already confronts the post-carve question head-on ("historically the *absence* of an affine gate was its content; today it is simply the family's shared shape"). Justification 2 (the #257 S3 multiplier product's designated home, with the 1/v member imminent) is a concrete chartered consumer — retiring now re-adds at S3. Justification 1 (isinstance taxonomy) has `[M]` 1 consumer, a test — stated as capability, not usage; honest. |
| 5 (DSA vestigial vocabulary) | CONFIRMED, widened | = C8. The whole file's noun, not just the param. |
| 6 (typed Krylov iterate) | STATE on the issue, don't build | Under V a GMRES iterate (a linear combination of fluxes) is representable as a typed flux, so `to_flat` is now an ADAPTER concern (scipy's interface) rather than an ontology necessity — the flat path's reason-to-exist narrowed. Posted to #289 (the role-erasure/FullField-generic issue). No stale prose found claiming the old necessity in iteration.py. |
| 7 (declared-law header) | WITHDRAWN | See Checked-and-CLEAN. |
| 8 (`_principal_bulk_leaf` native place) | CONFIRMED — relocate | The walk hard-codes two foreign types' anatomy (Composite's `interior`, CoupledField's `systems[0]`) in numerics' iteration module; its stated "type-agnostic by design" rationale was mid-carve scaffolding (surviving the step-2 flip) and is spent. Each owner already documents its own convention. Move: `Field.principal_bulk_leaf → self`; `Composite.principal_bulk_leaf → self.interior`; `CoupledField.principal_bulk_leaf →` first system's (the system-order convention's home). The loop degrades to one `getattr`; helper + `_NormedLeaf` retire with content migrated. |

### E1 — the five gates strengthen to the direct law (C12's constructive half)

`[M]` an affine `op(ψ) = Aψ + b` satisfies the blend identity exactly
(`op(ψ₁+λ(ψ₂−ψ₁)) = (1−λ)op(ψ₁)+λop(ψ₂)` — affine maps preserve affine
combinations), so the blend leg alone cannot red on the P3-class regression;
only paired with homogeneity is the old spelling discriminating. The direct
law `op(ψ₁+ψ₂) = op(ψ₁)+op(ψ₂)` reds on affine ALONE (b survives once left,
twice right — the declared-law file's measured result). Re-spell all five to
direct additivity, keeping homogeneity; comments state the law, with the
constrained-era spelling as one dated line.

## Filed / commented

- NEW issue (module:numerics, type:improvement): converge the three member
  vector-contract concepts (`Vector`, iteration's ravellable pair,
  `SystemField`) onto one named Protocol — the `coupled_system.py` "collapse
  trigger", whose blocker (torsor signatures unfit for `Self+Self→Self`)
  CS3 dissolved. → #391
- Comment on #289: the typed-Krylov design note (lead 6). → posted

## Completeness denominators (what was actually swept)

Stated at close — see the end of this memo after the census agent's
independent sweep is folded in.

### E1 measurement (2026-08-19, this review)

`[M]` slab `MultiplicationOperator` with affine offset `b = C.apply(ψ₁)`
(a genuine image element), seeds 51/52: blend-form residual **2.220e-16**
(blind — bit noise), direct-form residual **1.168e+00** (loaded). The five
re-spelled gates: 9/9 green post-re-spell (`test_apply_is_linear` ×
geometries, `test_apply_linearity`, `test_full_algebra_linearity`,
`test_apply_linearity_under_sweep_frame` × 2 geoms).

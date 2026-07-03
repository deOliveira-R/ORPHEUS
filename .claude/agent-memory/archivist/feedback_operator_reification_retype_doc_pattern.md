---
name: operator-reification-retype-doc-pattern
description: Doc recipe for a #226-taxonomy carve that REIFIES a duck-typed operator (confused apply/solve → matrix splitting; output-mode arg → typed composition) AND its step-3 SEQUEL (solver builds .inverse() / driver applies it; probe retirement; SupportsSeededApply family table; _seeded_inverse narrowing + #285). Category-error framing + reifying math + source-subspace honesty + tombstone + rejected-proxy-gate + deleted-class→past-tense-literal + close-the-transitional-forward-ref. Instance: #226 steps 2 (W2/W1) & 3.
metadata:
  type: feedback
---

When an operator-inverse taxonomy carve (#226 family) REIFIES a
duck-typed operator that violated the two-layer law of
`:ref:operator-algebra` ("two operators, one substrate — never two
views, one operator"), the theory-doc update follows this recipe.

**Why:** these carves recur (step 2 did W2 the G-S resolvent + W1 the
windowed path; step 3 DID delete `_MomentWindowedResolvent` and now holds
`P @ A.inverse()` directly — prediction FULFILLED, docs updated). Each
leaves a duck-typed predecessor whose doc still describes the OLD shape.
The mechanical repoint (L-002 grep gate) is the floor; the
derivation-level retyping story is the ceiling (Cardinal Rule 3).

**How to apply.** Two carve shapes, one recipe:

- **A confused apply/solve pairing → a regular matrix splitting.** The
  old object paired `apply = (L+C)ψ` with `solve = (L+C−B_lower)⁻¹` —
  inverses of DIFFERENT operators (round-trip defect O(1) = 2.667). The
  reification is the honest splitting `(L+C−B) = M − N`,
  `M = (L+C)−B_lower` reified as an `OperatorSum` whose `apply` is the
  leaf sum and `inverse()` is a genuine `SweepOperator` on M — the two
  FACES of one operator (round-trip now 5.2e-16/4.4e-16).
- **An output-mode arg that silently changes codomain → a typed
  composition.** `solve_moments` was a public method whose output-mode
  argument changed the operator's codomain. A codomain change is a
  DIFFERENT arrow (§1 two-layer law) = a COMPOSITION, never a config.
  The reification is `windowed = P ∘ A⁻¹`, `P` a block COISOMETRY
  (`analysis ∘ reconstruction = 4π·I`, NEVER `= I` — ERR-051), NOT
  invertible ⇒ the composite advertises `is_invertible=False` + no
  `CAP_SOLVE` (makes NO round-trip promise). `P @ A.inverse()` spells it;
  the fused `apply` deforests, the INHERITED `OperatorProduct.apply` is
  the retained fuller-view oracle (aggressive-retirement exception).

**The 6-move doc shape (both shapes share it):**

1. **NEW subsection** (own `.. _anchor:`) with the CATEGORY-ERROR
   framing FIRST: why the old shape was "a composition wearing a config"
   / "an operator whose faces disagree". This is the load-bearing WHY.
2. **The reifying math** as unlabelled `.. math::` (no new `:label:` —
   the gates are `@pytest.mark.foundation`, software invariants that
   NEVER carry `verifies()` and point at NO label, so a label would be a
   documented-only orphan; L-004). The splitting `(L+C−B)=M−N` + the
   row-split predicate (`SweepSchedule.lower_inflow_rows`: row in B_lower
   iff swept strictly AFTER its face's reflect), OR the coisometry
   identity + the two-arm composition.
3. **The source-subspace / domain HONESTY note.** The sweep re-derives
   the outflow-definition rows, so M⁻¹ is exact on the SOURCE SUBSPACE
   {y : y.outflow-rows=0} — which contains every production rhs — whose
   preimage is the trace-CONSISTENT states. Same property as the landed
   `(L+C).solve`; the feedback just makes it visible. The round-trip
   gate round-trips a CONSISTENT state → machine precision on bulk AND
   trace (STRONGER than the bulk-only falsifier). Record this so a
   future reader doesn't mistake the subspace restriction for a bug.
4. **Tombstone the historical design in-place** (L-007): a
   `.. note:: Retyping (#226 step N)` at the spot describing the retired
   shape (the "named methods solve vs solve_moments" impl-map block),
   flipping tense and forward-pointing to the new subsection. Do NOT
   delete the Phase-history narrative — it is accurate history.
5. **Repoint dead refs + the "what was tried and rejected".** Grep the
   WHOLE tree (L-002 — code-xrefs render plain-text, `-W` blind); the
   deleted class's blast radius spans theory + api + cross_section pages.
   Include the falsified alternative: the whole-face-OVERWRITE reflect
   dropped y_row (benign in production, O(1)-wrong as an inverse); a
   "moment-proxy RESIDUAL gate" was rejected as CATEGORY-CONFUSED — a
   coisometry composition has no `CAP_SOLVE`, so no round-trip to take a
   residual OF; its only honest correctness statements are
   representation-equivalence to the deforested oracle + a
   structurally-independent anchor (SI≡Krylov-full / closed-form k∞).
6. **Numbers + gates by name** (foundation-tagged): the round-trip /
   split-exactness / FP-invariance (Mode-9, on a DIAGONAL cubature that
   breaks the degenerate coincidence — never the isotropic-reflective
   box) + the mutation pair that replaced an unrepresentable spec
   mutation (here M-SPLIT became unrepresentable because the mask
   single-sources both `apply` and the reflect → replaced by
   M-SPLIT-DIR flip-the-split-direction + M-SPLIT-PART break-the-
   partition). Add the Development-history changelog row (reverse-chron,
   `*(in development)*` + branch name for an unmerged tree).

**Trap caught in passing:** the windowing gate predicate had drifted
(`reduced is None` was the pre-C5.4/#225 proxy; the genuine live gate is
`is_cartesian and ndim == 2`, which ALSO excludes 3-D Cartesian). The
live code explicitly names the proxy as rejected — so it was a
documented-WRONG claim (correct to the live gate, Cardinal Rule 1, flag
as beyond-brief scope). Verify the gate against the live `_maybe_window`
/ `_select_si_resolvent`, not the doc's frozen predicate.

**STEP-3 SEQUEL — the driver-consumption carve (not another
reification).** Step 3 ("solver builds `.inverse()`; SI/Krylov apply it")
is a DIFFERENT shape from steps 1–2 and gets its own moves on TOP of the
6-move recipe:

1. **A NEW `====` h1 section, not an h2 under the topic.** The
   driver-consumption model is GENERAL (SourceIteration/KEigenvalue/
   KrylovAcceleration all consume inverse operators), broader than the
   specific optimization the step-2 material sat under. Placed right after
   the step-2 h2 cluster (so the #226 step1→2→3 narrative stays
   contiguous), as a top-level h1. Anchor `inverse-application-driver`.
2. **Close the step-2 section's TRANSITIONAL forward-ref as the BRIDGE.**
   Step 2's doc ended with "...that adapter is transitional — step 3 makes
   the driver apply inverses directly and deletes it." Step 3 flips that
   paragraph to landed truth ("step 3 removed even that") and turns it into
   the forward-pointer INTO the new h1. (Same as closing a `.. todo::` —
   the promise becomes the bridge.)
3. **A DELETED class with NO successor class → past-tense double-backtick
   LITERALS, never a repointed `:class:`.** Unlike a rename (repoint to the
   new class), `_MomentWindowedResolvent` dissolved into the driver — there
   is nothing to repoint TO. So each of the 6 dead `:class:`/`:meth:` refs
   became a ``literal`` in past-tense prose ("the transitional
   ``_MomentWindowedResolvent`` adapter that once wrapped it is gone"). A
   `:class:` would render plain-text AND falsely imply existence; the
   literal correctly names history (L-002 de-role / L-007 tombstone). The
   ONE present-tense-FALSE claim ("...is therefore never even constructed"
   attributed to the live adapter) is the highest-priority fix — rewrite to
   the landed truth (`_maybe_window` returns the plain sweep off-gate; the
   driver holds the product directly).
4. **The new content the driver-consumption model needs:** (a) the
   seeded-kwarg FAMILY TABLE (`.. list-table::`: SweepOperator-curvilinear
   THREADS / SweepOperator-2D-Cartesian + InverseOperator + WindowedSweep
   accept-and-IGNORE, each with the WHY — the discriminator is "is the seed
   CONSUMED", so one class in two regimes = two rows); (b) the
   probe-retirement narrative (why it EXISTED = heterogeneous duck-typed
   `solve` sigs forced an `inspect.signature` probe; why it could GO = the
   family sig is now canonical `SupportsSeededApply`, so thread
   unconditionally); (c) CAP_SOLVE-ctor-gate → CAP_APPLY-only + the
   invertibility obligation MOVED to the `.inverse()` builder, so an
   apply-only step operator (the coisometry-factored product) is legitimate
   BY DESIGN; (d) `_seeded_inverse` as the single narrowing home + its
   SCOPED conformance (only the REACHABLE inverses carry the seeded sig,
   NOT every family `.inverse()`) → the OPEN structural-vs-convention
   decision (grep the issue #, cite it as an OPEN GitHub link).
5. **Adjacent-staleness: the driver-contract section goes HALF-stale.** The
   pre-existing "variadic driver" section said "the resolvent it must
   invert" — TRUE for Krylov (keeps forward L), now STALE for SI (applies a
   pre-inverted op). Per L-007 this is a live self-CONTRADICTION with the
   new section, so refine the ONE clause + forward-point (not a section
   rewrite); flag it as a scope note.
6. **vv-status clean by construction.** NO new eq-labels: the reifying/step
   math is UNLABELLED (matches the step-2 sibling), the gates are
   foundation/sentinel (seed spy) or verify EXISTING labels (the
   windowed×G-S pin reuses `harmonic-moment-projection` etc.). GREP the
   `verifies(...)` on the new gates first (L-004) to confirm no orphan
   target. Changelog: a SEPARATE step-3 row ABOVE step-2 (same date/branch,
   reverse-chron), naming the GATES not counts (the changelog omits counts).

Instance: #226 steps 2 & 3, branch `refactor/inverse-as-operator` —
step 2 `si-gauss-seidel-reification` (`discrete_ordinates.rst`) +
`windowing-retyped` (`operator_algebra.rst`); step 3
`inverse-application-driver` (new h1 in `operator_algebra.rst`) +
Development-history step-3 row (`discrete_ordinates.rst`). Source seed
(L-005): the rich docstrings of `iteration.py` (`SupportsSeededApply` /
`_seeded_inverse` / SourceIteration / KEigenvalue / KrylovAcceleration) /
`windowing.py` / `sweep_operator.py` / `operator.py` (InverseOperator) +
`test_seed_threading_spy.py` / `test_2d_anisotropic_windowing.py` /
`test_si_single_primitive_contract.py`.
Complements [[feedback-operator-classes-to-frame-faces-rehoming]].

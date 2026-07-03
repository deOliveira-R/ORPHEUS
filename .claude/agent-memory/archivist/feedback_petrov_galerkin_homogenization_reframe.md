---
name: petrov-galerkin-homogenization-reframe
description: Doc-recipe for re-framing a "flux-weighted homogenization" theory section from the FORWARD-ONLY "Galerkin in L²(φV)" reading to the honest PETROV-GALERKIN framing (#268 P3) — the bilinear/eigenvalue-consistency argument, measure-never-carries-discipline, two-frame design, the anchor-rename + cross-doc-ref hygiene.
metadata:
  type: feedback
---

The live successor to [[galerkin-natural-metric-reframe]] (which the #268
P3 carve REVERSED). When a section argues a flux-weighted
homogenization/condensation is "Galerkin in the natural L²(φV) (resp.
spectrum) metric", that framing is **forward-only** and must be rewritten
to **Petrov-Galerkin**. The production code lands the honest framing first
(an explicit `WeightedIndicatorBasis` test side on a `PetrovGalerkinFrame`);
the archivist brings the theory page into agreement.

**The two load-bearing rulings (the physics correction):**

1. **Homogenization & condensation are PETROV-GALERKIN, not Galerkin.**
   The "Galerkin in L²(φV)" reading folds the solution (flux) into the
   *metric* — legitimate ONLY for forward-flux, reaction-rate-only
   reduction. Eigenvalue-consistent homogenization (what reactor physics
   requires) preserves the **bilinear** ⟨φ*, Σφ⟩ (k is stationary w.r.t.
   the ADJOINT-weighted residual by 1st-order perturbation theory),
   giving Σ_R = ∫_R φ*Σφ dV / ∫_R φ*φ dV with **test=φ*·1_R ≠ trial=φ·1_R
   ⟹ M*≠R** — irreducibly PG. The forward φ*=φ is its **Galerkin
   degenerate** (a *legal* Galerkin reading because of the coincidence
   φ*=φ — but the coincidence is NOT the structure). The genuine
   adjoint-weighted PG (φ*≠φ) is a later phase (P6); set it up as the
   non-degenerate sibling and forward-reference it.
2. **The measure carries the axis + the fixed L² metric, NEVER the
   discipline.** "Folding flux into the metric" = "making the measure
   carry the discipline" = the mistake. The flux is a **test-weighting the
   solution emits**, on the test side = the frame TYPE.
   `DiscreteMeasure.weights` stays the pure geometric dV.

**Section-rewrite order (Sphinx-as-brain, full derivation + WHY):**
(a) headline sentence: homogenization = the coefficient extraction
`project = G⁻¹M` of a PG frame (trial = plain cell indicator, test =
flux-weighted indicator χ_R=φ·1_R, measure = geometric dV); (b) a
top-of-section `.. warning::` retraction tombstone (L-007) naming the OLD
"L²(φV)-Galerkin" claim, one-line WHY it's wrong, forward-ref to the
why-PG subsection; (c) trial space / test space / **cross-Gram**
G_RS=⟨χ_R,1_S⟩=δ_RS·Φ_R (diagonal by disjoint support — but the two
factors are DIFFERENT bases, contrast the SH GalerkinFrame's *symmetric*
Gram); (d) a dedicated **"Why PG not Galerkin"** subsection: the
metric-fold ⟨φ·1_R,Σ⟩_dV=⟨1_R,Σ⟩_φV is a real identity FOR φ*=φ ONLY →
the bilinear ⟨φ*,Σφ⟩ eigenvalue-consistent case → "no metric makes
test=trial" → forward = Galerkin degenerate; (e) measure-carries-axis-
never-discipline (the test weight rides the test basis, weights stays
1-D); (f) mesh-yields-basis + n-D (membership ravel_multi_index "ij" =
measure node order); (g) the G⁻¹M verb + Moore–Penrose zero-flux→0; (h)
**two-frame design** as a `.. list-table::` (Σ→reaction rate, flux test;
χ→emission rate, production test p=Σ_g νΣ_f φ_g — convex avg of simplices
→ simplex; both share trial indicator + dV; route via
`MaterialXSField.project_through`); (i) Mode-11 sentinel note (the
TEST-side reader `WeightedIndicatorBasis.analyze` is the load-bearing
count — a metric-fold regression keeps the numbers but leaves that
counter at 0).

**Hard cross-ref/label hygiene (the traps this task hit):**

- **The Mode-11 sentinel test was RENAMED in the concurrent code carve.**
  The old doc cited `test_homogenize_routes_through_the_indicator_frame`;
  the live test is `test_homogenize_routes_through_the_petrov_galerkin_frame`
  (counts `WeightedIndicatorBasis.analyze` + `FrameBase.project` now). A
  dead test-name in a doc is a plain-text staleness bug (L-002). ALWAYS
  re-inventory the cited test NAMES against the live test file
  (`grep "def test_"`), not the brief — the test file is reframed
  CONCURRENTLY with the code (sharpens L-001/L-005).
- **A NEW n-D test** (`test_homogenize_2d_rate_preservation`) was added as
  a THIRD `verifies("sn-homogenization-rate-preservation")` binding — the
  matrix auto-regen shows the verifies-count 2→3. Add an n-D bullet to the
  Verification subsection; it's not my drift (L-008).
- **Anchor rename + cross-doc ref.** The section anchor was
  `sn-homogenization-galerkin-frame` (now WRONG word). RENAME →
  `...-petrov-galerkin-frame`, update all 3 in-page `:ref:` + the **1
  cross-doc `:ref:`** in the sibling `galerkin_projection.rst` (a cross-doc
  dangling `:ref:` renders plain-text with NO -W warning, L-002). The
  sibling page had already been #268-corrected and forward-pointed here
  with "until it is rewritten" — UPDATE that prose to past-tense ("was
  rewritten") since the rewrite is now done. Confirm in the built HTML the
  cross-doc ref renders as `href=...#new-anchor` and the old anchor has
  ZERO live links.
- **Label keep/retire (L-003/L-004).** KEEP `-rate-preservation` by name +
  body (verifies-target, 3 tests, no vv-status — it's the L0 claim). KEEP
  every documented-only label whose math is unchanged. The
  `-normal-equations` label: keep name, rewrite BODY to the PG cross-Gram
  normal eqns. RETIRE the WRONG `-galerkin-equals-petrov` "two-readings-
  same-map" label; mint 3 new documented-only labels: `-test-functions`
  (χ_R=φ·1_R), `-metric-fold` (the forward-only identity, flagged), and
  `-bilinear` (the eigenvalue-consistent ⟨φ*,Σφ⟩, tagged "built P6, not
  yet a solver claim"). Each new label gets `.. vv-status: <label>
  documented`. matrix.rst auto-regen handles the +2 net (do NOT hand-edit).
- The `.. note::` "space-only, **1-D** slice" must become "**space-only**,
  dimension-agnostic" (the 1-D guard was dropped). The intro `:meth:`
  signature `takes a coarse Mesh1D` → `Mesh1D or Mesh2D`.

**Quality self-assessment (Directive 3):** Derivation-depth 5 (added the
full bilinear/perturbation-theory eigenvalue-consistency argument, the
cross-Gram with the SH-symmetric-Gram contrast, the forward=degenerate
framing — the old draft LACKED the WHY-it's-PG); Cross-refs 5 (every new
target import+hasattr-verified live; cross-doc ref confirmed as a live
href in built HTML; 3 cited test names confirmed green under -O);
Numerical-evidence 4 (the φV-vs-dV discriminator + Mode-11 TEST-side-reader
sentinel + the n-D gate are the evidence; no convergence table — exact
projection, structurally absent per the rubric note, not a deficit);
Failed-approaches 5 (the OLD L²(φV)-Galerkin claim preserved as a retraction
`.. warning::` with the WHY-it-looked-right metric-fold; the
forward-degenerate explained); Code-traceability 5; Derivation-source 5
(the frame.py / weighted_indicator_basis.py / material_xs_field.py /
solution.homogenize docstrings ARE the algebra-of-record — read them first
per L-005; they were already richly written and confirmed every ruling).
Weakest = numerical evidence, STRUCTURALLY (a forward-flux exact projection
+ a type/discipline re-frame produce no NEW flux to tabulate; the
distinction is invisible numerically on this branch — it's an architecture
distinction, so I said so rather than manufacture a table).

See [[lessons]] L-001 (verify cited test NAMES against live file, not brief
— they reframe concurrently), L-002 (cross-doc `:ref:` renders plain-text;
numerics.frame/basis NOT automodule'd → `:class:`/`:meth:` plain-text by
page convention, matched — no half-surfacing), L-003 (keep verifies-target
label, retire the wrong one), L-004 (vv-status on the 3 new representational
labels; the `-bilinear` P6-future label tagged "not yet a solver claim"),
L-007 (retraction tombstone preserves the WHY), L-008 (matrix auto-regen
+ concurrent test counts not my drift), L-010 (V&V vocab: "MMS/eigenvalue"
discipline — I wrote "k is stationary w.r.t. the adjoint-weighted residual
by 1st-order perturbation theory", the eigenvalue-consistency claim, matched
to Hébert §6 not vaguely "analytical").

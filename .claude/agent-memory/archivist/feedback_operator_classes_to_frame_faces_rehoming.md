---
name: operator-classes → frame-faces re-homing doc sweep
description: Doc-sweep recipe when a refactor retires standalone operator classes (a projection M + a reconstruction R) into the TWO FACES of one abstraction (a discrete frame) — re-home onto the abstraction, don't find-replace; add the harmonic-analysis framing it earns; KEEP documented-only eq-labels by name. NOW ALSO the canonical record of the #268 REVERSAL: discipline is a TYPE (FrameBase→PetrovGalerkinFrame→GalerkinFrame), NOT a property/role-marker; homogenization/condensation are PETROV-GALERKIN (test=φ/φ*-weighted basis), NOT "Galerkin in L²(φV)"; measure carries axis+fixed-L²-metric, never discipline. Plus: sibling-page-staleness-FLAG-not-silent-flip when prose contradicts the new ruling (repoint to discipline-neutral FrameBase, flag the substantive reversal as a separate task); M-not-Π total relabel; auto-regen-from-concurrent-tests is not my drift.
type: feedback
---

A second-generation sibling of [[feedback_decomposition_leaf_retirement_rationale]]
and [[feedback_type_retirement_concept_survives]]. Those retired a typed
field or a decomposition-into-leaves. This one retired **two standalone
operator classes** (`MomentProjection` = M = Y*W, and
`HarmonicMomentReconstruction` = R = (2ℓ+1)·S₀) into the **two operational
FACES of one generic abstraction** — the discrete
`Frame` (`orpheus.numerics.frame.Frame`): `frame.analysis` (M) and
`frame.reconstruction` (R). The (2ℓ+1) literal moved onto
`SphericalHarmonicBasis.addition_theorem_factor`. The ABCs
`ProjectionOperator`→`AnalysisOperator` were RENAMED (kept as
forward-looking discipline vocabulary alongside the new mechanism).

**SUPERSEDED RULING (2026-06-24, #268): discipline is a TYPE, not a
PROPERTY.** The entry below originally recorded "discipline becomes a
PROPERTY of the frame". The P1 discipline-type carve REVERSED that: the
shipped architecture is a Liskov hierarchy `FrameBase → PetrovGalerkinFrame
→ GalerkinFrame` (the single `Frame` class was retired). The
`GalerkinProjection`/`PetrovGalerkinProjection` marker ABCs on the operator
ROLE were retired too; `projection.py` keeps ONLY the two discipline-free
operator roles `AnalysisOperator`/`ReconstructionOperator`. AND:
**homogenization/condensation are PETROV-GALERKIN, not "Galerkin in
L²(φV)"** — the L²(φV) reading folds the solution into the metric, legit
only for forward-flux reaction-rate reduction, breaks under eigenvalue-
consistent (adjoint-weighted) homogenization. So **the measure carries the
axis + the FIXED L² metric, NEVER the discipline; the solution-weighting
(φ, φ*) is a first-class TEST BASIS on the frame TYPE.** The eigenvalue-
consistent (φ*-weighted) case is the headline future PG consumer. The
re-home-onto-the-abstraction rule (next paragraph) STILL HOLDS; only the
discipline-as-property ruling flipped.

**Rule: re-home onto the abstraction, do NOT find-replace names.** The
brief said "MomentProjection → frame.analysis". A mechanical swap would
have left the prose describing two unrelated operator classes that happen
to be spelled differently. The correct move (Cardinal Rule 3) is to
re-narrate the page so the NEW abstraction's concept is the spine: the
(R, M) pair IS one frame's reconstruction/analysis faces; the discipline
(Galerkin vs Petrov-Galerkin) is carried by the frame TYPE (test=trial vs
test≠trial), not baked into a class name and NOT a property. ADD the
framing the new abstraction earns and the old classes
lacked — here, discrete-frame theory (Christensen 2016): analysis operator
T, synthesis operator T* = S₀, frame operator S = T*T, tight frame
(the 4π-tightness IS the Π R = 4π I identity, c_V = the tightness
constant), canonical dual. That is the load-bearing value-add over a
rename.

**The brief's "X is being retired" framing meant the code carve was
LANDING CONCURRENTLY.** I read projection.py early (saw the OLD
`ProjectionOperator` + `MomentProjection` still present), wrote docs to
the TARGET vocabulary per the brief, then re-verified at the end — and
projection.py had CHANGED mid-session to `AnalysisOperator` + ABCs-only
(MomentProjection/HarmonicMomentReconstruction gone from `__all__` and the
file). The end-of-sweep `.venv/bin/python -c "import + hasattr"` liveness
check on every NEW cross-ref target is what caught this and confirmed my
docs were correct against the FINAL live tree, not the early read. LESSON
(sharpens lessons L-001): on a "being retired / migrating" brief, the code
is a MOVING TARGET during your session — the authoritative check is a
programmatic import+hasattr of every new target at the END, via the venv
(`/Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python` — system python has
no numpy), NOT the read you did at the start.

**Co-located SEPARATE retirement — the scoping discipline.** While
sweeping I found a DIFFERENT, earlier retirement whose doc-sweep was never
done: `evaluate_real_sh` / `orpheus.numerics.spherical_harmonics` (module
deleted in commit `09c4241`, P3.2) and `AngularQuadrature` /
`orpheus.sn.quadrature` (→ `Quadrature` / `orpheus.numerics.quadrature`).
~39 dead `:func:`/`:class:` across 9 files, all rendering plain-text (no
-W warning). DECISION: fix these phantoms ONLY where they co-occur in a
prose block I was already rewriting for the in-scope M/R migration (the
block must point at a REAL evaluator — `SphericalHarmonicBasis.evaluate` —
to be correct; leaving a dead `:func:` in a block I'm editing is a
Cardinal-Rule-1 bug). For the rest (blocks/files I'm NOT otherwise
touching — e.g. boundary_conditions.rst, discrete_measures.rst, and a
distant discrete_ordinates §) I FLAGGED them as a distinct follow-up
rather than scope-creeping a second retirement's full sweep into this
brief. The discriminator: "am I already rewriting this exact block for my
brief?" → fix the co-located phantom; "is this a separate section/file?"
→ flag, don't chase.

**Eq-label discipline (L-003 applied, clean outcome).** The 4 galerkin
labels (`galerkin-pair`/`-self-adjoint`/`-construction`,
`petrov-galerkin-construction`) and `pi-r-equals-4pi-i` are all
`.. vv-status: <label> documented` (structural/representational, NOT
verifies-targets). The math is IDENTICAL — only the operator NAMES change.
So: KEPT every label by name, rewrote only the BODY prose / operator-class
refs. NO label added or removed → NO matrix regen, NO new orphan target.
Confirmed by grepping `docs/verification/matrix.rst` for each label FIRST.

**No-automodule discipline (L-002 exception, load-bearing here).** Did NOT
add an `automodule` for `frame`/`basis`/`projection` — the api/numerics.rst
convention automodules ONLY `eigenvalue` + `field`; everything else renders
`:class:` as plain-text by page convention. CRITICAL extra reason:
`spherical_harmonic_basis.py` carries `.. math:: :label: real-sh-addition-theorem`
in its docstring, and that label is ALREADY defined in
`spherical_harmonics.rst` — automodule'ing the basis would DUPLICATE the
label (a real warning). The api listing uses prose `:class:` bullets (no
automodule), matching the page + dodging the collision. When a new
abstraction's docstring carries a `:label:` already owned by a theory
page, NEVER automodule it — cross-reference in prose.

**Tombstone-the-old-name in the History section (L-007).** The galerkin
page's "History — from operator classes to the discrete frame" section
KEEPS the old names as `` ``HarmonicMomentProjection`` `` literal
code-spans in an explicit old→new mapping
(``HarmonicMomentProjection`` → ``frame.analysis``). A final grep for the
retired symbols showing ONLY these literal code-spans (zero live
`:class:`/`:meth:`) is the correct end-state — the history record names
what was retired; the live cross-refs all point at the survivors.

**SIBLING-PAGE STALENESS FLAG, not silent-flip, when the substantive
prose contradicts the new ruling (P1 #268 correctness pass).** The dead
`:class:`Frame`` xref sweep (retirement-audit search #2) reached the
`discrete_ordinates.rst` HOMOGENIZATION section, whose entire prose still
argues the WRONG "Galerkin in L²(φV)" framing (a `.. _sn-homogenization-
galerkin-frame:` section with its own eq-labels + vv-status directives).
The correct xref type per the new ruling is `PetrovGalerkinFrame`, BUT
flipping the xref to PG while leaving prose that says "it sees the Galerkin
case" creates an internal contradiction. DISCIPLINE: repoint the dead
symbol to the GENERIC base `FrameBase` (accurate regardless of which
discipline the surrounding prose argues — the binding mechanism IS
FrameBase) and FLAG the homogenization-section prose reversal as a distinct
follow-up, do NOT scope-creep a whole-section correctness rewrite into a
P1 xref+single-page pass. The brief scoped the correctness rewrite to ONE
page (`galerkin_projection.rst`); the sibling page's stale derivation is
its own task (Issue #268). Generalizes the co-located-separate-retirement
discriminator below: when a dead-symbol hit lands in prose that makes a
SUBSTANTIVE claim the new ruling reverses, the minimal-correct fix is the
discipline-NEUTRAL base type, and the substantive reversal gets a flag +
forward-pointer (I added a cross-link from the rewritten page's status note
to the stale sibling section). On the rewritten page itself, the two
rejected framings (marker-ABC-on-role; discipline-as-property) BOTH get a
named "Rejected (a)/(b)" subsection in the discipline-as-type section —
preserve the WHY of each dead path (L-007), don't just delete them.

**M-not-Π relabel must be TOTAL.** The new code names the analysis face
`M` (= `frame.analysis`); the old page used `Π` throughout. Relabeled
EVERY `\Pi` → `M` (analysis operator) and `\Pi^*` → `M^*` (Hilbert
adjoint), kept `R` (reconstruction). Final `grep '\\Pi'` on the page = zero
is the consistency gate — a stray `Π` would read as a different operator
than the code's `M`. The Galerkin promise stays `M^* = R`; PG stays
`M^* ≠ R`.

**Auto-regen of matrix.rst from CONCURRENT test changes is NOT my drift.**
`matrix.rst` showed ` M` post-build — but the +4 rows were the user's
new `test_frame.py` PG-hierarchy tests (8→12 foundation), picked up by the
`-E` regen, NOT from my label edits (I preserved all 4 documented labels by
name). Per L-008, report the real regen number; do NOT revert a generated
artifact. Same for `frame.py`/`projection.py`/test files showing ` M`:
those are the user's P1 carve (the brief told me to READ them as ground
truth), NOT my edits — I touched only docstrings in OTHER files + the RST.

**Title-rename underline (L-009).** Renamed the galerkin page title
(58 code points) but its overline/underline were 60 (from the longer old
title). Docutils tolerates over-length, but I resized to 58 for
cleanliness. Size with `len(title)` in python; the em-dash `—` is 1 code
point. New sections added were all `=`-level (the page is flat = only) —
introduced no sub-section marker, so no skip-level risk.

**Quality self-assessment (Directive 3):** Derivation depth 5 (added full
frame-theory derivation — T/T*/S/tight-frame/canonical-dual — the old
operator-class pages lacked); Cross-references 5 (every new target
import+hasattr-verified live); Numerical evidence 4 (preserved the
existing Π R = 4π I residual table + Lebedev convergence; reframed as
tightness — no NEW table because no flux moves in a re-homing); Failed
approaches 4 (preserved the ERR-039 four-operator-conflation history +
the rejected strict-ΠR=I-by-dividing-weights note); Code traceability 5;
Derivation source n/a (no SymPy module — this is an operator-algebra
re-homing, the "derivation" is the frame-theory identities which are
test-pinned, not derivations/-scripted). Weakest = numerical evidence,
STRUCTURALLY (a name/abstraction re-homing produces no new flux to
tabulate) — per the rubric note, say so rather than manufacture one.

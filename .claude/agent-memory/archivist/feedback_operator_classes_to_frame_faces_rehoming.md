---
name: operator-classes → frame-faces re-homing doc sweep
description: Doc-sweep recipe when a refactor retires two STANDALONE operator classes (a projection M + a reconstruction R) into the TWO FACES of one generic abstraction (a discrete Frame) — re-home onto the abstraction, don't find-replace names; add the harmonic-analysis framing the new abstraction earns; KEEP the documented-only eq-labels by name. Instance: Frame/Basis carve (refactor/operator-inverse-algebra), MomentProjection/HarmonicMomentReconstruction → frame.analysis/frame.reconstruction. Also: the concurrent-code-carve vindication of re-verify-live, and the co-located-separate-retirement scoping discipline.
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

**Rule: re-home onto the abstraction, do NOT find-replace names.** The
brief said "MomentProjection → frame.analysis". A mechanical swap would
have left the prose describing two unrelated operator classes that happen
to be spelled differently. The correct move (Cardinal Rule 3) is to
re-narrate the page so the NEW abstraction's concept is the spine: the
(R, Π) pair IS one frame's reconstruction/analysis faces; the discipline
(Galerkin vs Petrov-Galerkin) becomes a PROPERTY of the frame (which dual
it uses: canonical-dual vs oblique-dual), not a fact baked into a class
name. ADD the framing the new abstraction earns and the old classes
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

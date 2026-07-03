# Archivist — Lessons

Read at the START of every invocation. Behavioral corrections only:
"what documentation/Sphinx/knowledge-architecture mistake did I make,
and what discipline did it teach?" The mechanical HOW lives in
`AGENT.md` ("Build-Gating & Cross-Ref Reality", "Close-Out Narrative
Arc") and in the preloaded skills (`vv-principles`, `algebra-of-record`,
`nexus-verification`). This file holds the PROCESS lessons those don't.
Campaign play-by-play (commit hashes, codenames) is retired — a lesson
here stands WITHOUT knowing what the campaign was.

The cross-cutting standard behind most of these: **a doc page is not
done when the prose reads well — it is done when every cross-ref
resolves against the LIVE tree, every claim's V&V level matches the
skill verbatim, every retired symbol leaves no dangling reference, and
the build's WARNING/ERROR/CRITICAL set is unchanged from the pre-edit
`-E` baseline.** Each lesson is one face of that.

---

## L-001 — Verify every claim against the LIVE code, NOT the quoted/docstring prose

→ **The standing directive now lives in AGENT.md Quality Checklist item 6**
(read the live source before citing any convention/shape/decision/result;
brief, docstring, and verdict memo can ALL be stale). The three war-story
faces below are kept for forensic value — they show HOW each surface lies.

The single most recurring trap. A task brief quotes "current stale
text", or a docstring describes a return shape, or a verdict memo
states a recommendation — and ALL THREE can be wrong relative to what
shipped. Three faces:

- **Brief quotes "stale text" that was already fixed.** The user's
  snapshot may pre-date a landed fix. Read the file (or `git show
  <prior-commit>:<path>`) before editing; if the current text already
  matches the desired "after", report "already fixed" rather than
  re-introducing the stale wording.
- **A docstring lies about a convention/shape.** A physics convention
  in a docstring (e.g. an index `r_i` vs `r_j`), or a return layout
  `(N, nx, 1, 1)` claimed while the array is `(N, ng, nx, ny)`, is
  ground-truth ONLY in the live code body. Read the array build / the
  consuming code; never transcribe the docstring's claim into a theory
  page.
- **A verdict memo records the RECOMMENDATION, not the OUTCOME.** An
  elegance/qa verdict's "BLOCKING NIT: drop X" may have been overruled
  — the shipped code took the alternative ("keep + strengthen test").
  Read the CODE to learn which resolution won.
- **A retirement-SHIM docstring freezes a claim that a LATER refinement
  RE-PERMITTED in the canonical forms.** A deprecation shim
  (`sn/sources.py`) carries a docstring stating "the cross-class `iso +
  aniso` dunder is RETIRED, use `from_isotropic` + within-class add" —
  and a brief built on the shim faithfully repeats it. But the shim is
  FROZEN at the commit that minted it; the canonical L2 forms it
  re-exports (`transport/source_sinks/{scalar,angular}_source_sink.py`)
  EVOLVED PAST it: the refined Issue #207 principle (recorded later)
  RE-PERMITTED the cross-class dunder via canonical subspace-containment
  injection (`iso → 1 ⊗ iso`, returning the larger type, commutative).
  Both `__add__` bodies + `__radd__` are live and wired; the module
  docstring (the algebra-of-record) says "PRESERVED". Following the
  brief's "retired" premise would have wrongly PAST-TENSED an accurate
  doc section (the "load-bearing dunder" §) — the live code MATCHED the
  doc, only the names were stale. RULE: when a brief's behavioral
  premise is sourced from a SHIM docstring, verify it against the
  CANONICAL form the shim re-exports (read both the `__add__` body AND
  the canonical module docstring), never against the shim — the shim's
  job is to redirect, and its prose is a snapshot that the canonical
  layer can overrule. (Worked: `sn.sources` retirement, branch
  `refactor/operator-inverse-algebra` — corrected the doc to the live
  subspace-containment narrative, not the shim's "retired" claim.)
- **The brief's OWN discriminator can be over-broad — it names a
  concept, but the LIVE code applies that concept at MULTIPLE layers
  with DIFFERENT types.** A "type-confinement" sweep brief (e.g.
  "confine `TimedFullField` to the driver; the operator apply contract /
  the **solve boundary** / the composite the operators speak → now
  `FullField`") gives a discriminator that is a STARTING HEURISTIC, not
  a per-ref rule. The phrase "solve boundary" named the *within-group
  operator* solve (W-C confined it → `FullField`) — but the SAME page
  also documents the *public* `solve_sn_fixed_source(external_source:
  …)` source argument, whose live signature is STILL
  `np.ndarray | TimedFullField` (the driver-iterate-compatible composite
  source). Applying the discriminator blindly would have wrongly flipped
  an entire accurate section. Resolve EACH ref by reading the LIVE
  signature of the exact symbol it describes: the within-group
  `evaluate_residual` RETURNS `FullField` (UPDATE the doc's stale
  `TimedFullField(bulk=AngularResidual,…)`); the public
  `_build_fixed_source_rhs` RETURNS `TimedFullField` (KEEP). One brief
  phrase, two live boundaries, opposite verdicts. The operator-contract
  refs ("`X.apply` acting on …", the `@overload` surface, `B.apply`
  operates on …) flip to `FullField`; the driver-iterate refs (the
  composite rhs the SI/Krylov inner BUILDS, the `to_flat`/`from_flat`
  GMRES round-trip, `TimedFullField.zeros` initial-guess, the class-gate
  arithmetic, "Layer-4 never sees … carriers") stay `TimedFullField`.
  When a gate EXERCISES an operator on `TimedFullField` via MRO (the
  G-adjoint reciprocity gate), that is the inheritance path, NOT the
  declared contract — the doc states the CONTRACT (`FullField`); add a
  one-clause "the driver's iterate reaches it via MRO" so the KEEP-side
  reader isn't confused.

- **REPRODUCING a cited numerical result during a rewrite routinely
  surfaces a flatly-WRONG pre-existing worked example — fix it (Cardinal
  Rule 1 is supreme), don't transcribe around it.** When the brief gives
  you numbers to verify (here: the n2n double-count moves k by 0.43), run
  the LIVE derivation to confirm them — and while you have the repro
  harness open, sanity-check the page's OTHER worked numbers in the same
  section. A convention-fix brief surfaced a 2-group worked example whose
  `M₂₁=2.0833`/`M₂₂=3.666` gave trace 3.875 while the page then claimed
  k=1.875 three lines later — internally self-contradictory. The live
  `M=A⁻¹F` is rank-1 `[[0.2083,1.6667],[0.2083,1.6667]]` (trace=1.875=k,
  the OTHER eigenvalue is exactly 0 because `F=χ⊗νΣf` is a rank-1 dyad).
  The page's `two-group-M` formula carried a spurious `+ν₂Σf₂/Σr₂` term
  (assumed fission emission in group 2, but `χ₂=0`). This is a teaching
  doc — a self-contradictory worked example is a Cardinal-Rule-1 defect.
  Fix the formula + the worked numbers + tie the result to the trace,
  KEEP the verifies-target labels, and FLAG the scope-expansion in the
  return (the fix went beyond the brief's n2n convention scope). The
  CORRECT final value (k=1.875) was always right — only the intermediates
  were wrong — so this is an arithmetic CORRECTION, not a falsification:
  no tombstone, just a self-consistent rewrite.
- **A deleted class's doc blast radius is the WHOLE `docs/` tree, not the
  brief's named page — grep it and repoint the API pages too.** A brief
  scoped to `docs/theory/X.rst` after a class deletion is the FLOOR; the
  retirement audit's doc search (`grep -rn "<DeletedClass>" docs/`, minus
  `_build`) is the real blast radius. Here #276 deleted `HomogeneousSolver`
  and the brief named only the theory page — but `docs/api/homogeneous.rst`
  AND `docs/api/numerics.rst` BOTH cited the dead class AND carried the
  SAME retired n2n-in-both-matrices convention (the API page is a parallel
  surface that goes stale identically). These render plain-text (no `-W`
  warning, L-002), so the grep gate is the only catch. Repoint the dead
  `:class:` to the live `:func:`, rewrite the parallel convention error,
  and remove the dead class from any "reference implementations of
  protocol P" list (the refounded solver may no longer satisfy P — here
  the direct function is NOT an `EigenvalueSolver`). FLAG (don't silently
  fix) deeper adjacent staleness in those pages outside the deletion's
  radius (here `moc.solver.MOCSolver` is a stale module path — live class
  is `moc.core.MOCSolver` — pre-existing, not #276's; reported, not fixed).
- **The brief's RATIONALE can be subtly wrong, not just its facts — verify the
  ARGUMENT against the math and CORRECT it with a qualifying note, never
  transcribe.** A brief may give a CORRECT conclusion via a FLAWED reason:
  "0-D homogeneous uses the direct engine because an iterative one would be
  dominance-ratio-fragile at the 1e-12 gate" is true IN GENERAL but FALSE for
  THIS problem — the homogeneous `F=χ⊗νΣf` is RANK-1, so `A⁻¹F` has a single
  nonzero eigenvalue (dominance ratio 0) and power iteration would converge in
  ONE step here. The honest doc KEEPS the general argument AND adds a
  `.. note::` recording the rank-1 subtlety (one-step is a fragile consequence
  of F's rank-1 structure, not a guarantee; the direct engine is exact for ANY
  F). Transcribing the reason verbatim would mint a Cardinal-Rule-1-wrong
  teaching claim. The TEST file usually already encodes the subtlety (here
  `TestGeneralizedEigenproblem`'s docstring states F is rank-1 with a single
  nonzero eigenvalue) — read it. Same family: a brief's SIGNATURE can be wrong
  (`power_iteration(A,F)` → really `power_iteration(solver: EigenvalueSolver)`;
  it sees a Protocol boundary, the dense `(A,F)` is never formed — only
  `direct_eigenvalue`/`rqi` take dense `(A,F)`). Read the live `def` first.
- **A MECHANICAL vocabulary-translation can restate a substantively-FALSE
  claim in fresh words — verify the CLAIM'S TRUTH, not just its vocabulary,
  before re-spelling it.** A retirement pass that swaps a retired term for its
  successor (`CAP_APPLY_TRANSPOSE`→`is_adjointable`, `MissingCapability`→
  `MissingAdjoint`) FEELS purely lexical and scoped — but the underlying CLAIM
  may be stale for an ORTHOGONAL reason, and a faithful vocabulary-swap then
  mints a NEW false claim in cleaner language (WORSE than the old dead-ref: it
  reads authoritative). Worked (frozenset-retirement W3): the `operator_algebra`
  G-adjoint section asserted "S/F are not adjointable → the full-loss `.H` is
  unreachable/raises". Translating `CAP_APPLY_TRANSPOSE`→`is_adjointable` I first
  restated it as "S/F are not adjointable" — but LIVE `ScatteringOperator`/
  `FissionOperator` `.is_adjointable` return **True** (they gained `apply_transpose`
  via #112/#118/#276, orthogonal to the frozenset). The whole reachability
  argument was pre-existing-stale. Cardinal Rule 1 forbids shipping (or fresh-
  minting) the false claim, so I corrected it to the verified truth (metric lives
  on the shared `full_field_space`, every leaf carries it → no composite is
  metric-blind; the within-group loss never fuses S/F/B — `_within_group_triple`
  returns `(L+C, S, B)`) + an L-007 supersession `.. note::`, and FLAGGED the
  scope-expansion. RULE: any `:class:`/CAP_* symbol you translate, grep the LIVE
  class for the property it asserts (`grep 'def is_adjointable' … | show return`)
  — the vocabulary swap is safe only after the CLAIM is re-verified. A clean
  `-W` build never catches this (plain-text refs + a false prose claim both
  build green).
- **A carve's OWN sibling changes (adjoint work, pre-inversion redesign) can
  have staled the SAME sections the carve's doc-pass touches — fix the carve's
  dead refs, but FLAG (don't silently deep-rewrite) the sibling-staleness.** The
  `discrete_ordinates` "Capability requirements" posing section carried BOTH
  step-6 dead refs (`CAP_APPLY`/`MissingCapability` — fix required) AND step-3
  pre-inversion staleness (the `inverter`-hook narrative — `SourceIteration` now
  takes a pre-inverted `A_inv`; `KrylovAcceleration.inverter`→`preconditioner`).
  I rewrote the bullets to the verified live contract (SI→`A_inv`+`TypeError`;
  `KEigenvalue`→`A.is_invertible`+`NotInvertible`) AND added a `.. note::`
  flagging the surrounding `inverter` narrative as step-3-stale + deferred. The
  dead-ref fix is in-scope; the deep behavioral rewrite is a separate task.

How to apply: before citing any convention, shape, or design decision,
confirm it against the live source this session. Cross-refs that render
plain-text carry NO warning (see L-002), so a dead/stale `:func:` is a
Cardinal-Rule-1 correctness bug `-W` will never catch — grep the symbol
across the WHOLE `docs/` tree, not just the brief's named page.

---

## L-002 — Unresolvable code-xrefs render as PLAIN TEXT with no warning; this repo is NOT nitpicky

→ **The standing directive now lives in AGENT.md Quality Checklist item 3**
(grep-gate cross-refs; `-W` is blind to a dead code-xref). The detail below —
which ref classes DO warn, and the not-member-`automodule`'d page convention —
is kept as the recall companion for when a cross-ref edit gets subtle.

`-W` does NOT catch a dead `:func:`/`:class:`/`:meth:`/`:attr:` or a
stale alias-xref — they silently render as plain text. The acceptance
gate (count-unchanged from the `-E` baseline) only proves you added no
NEW warning; it is BLIND to staleness.

- After any carve that DELETES or RENAMES a symbol, `grep -rn "<symbol>"
  docs/` and repoint every hit on correctness grounds.
- Distinguish what DOES warn: undefined `[Key]_` citations warn (even
  non-nitpicky); intra-doc dangling `:ref:` warns under `-W`; cross-doc
  dangling `:ref:` renders plain-text. A new `:ref:` to a not-yet-
  existing section MUST create the labelled section in the SAME edit.
- **FORWARD-ref corollary: a NOT-YET-BUILT deliverable's symbol is a
  code LITERAL, never a `:class:`/`:meth:` cross-ref.** When a docs pass
  documents a DEFERRED seam (a changelog row for an open issue, a
  "future hook" bullet — e.g. #280's `SweepOperator.apply_transpose`,
  the `A.H.inverse()` swap law), the planned method/class does NOT exist
  yet. A `:meth:`-ref to it renders plain-text with NO `-W` warning
  (L-002 forward-facing), so the build is BLIND — but it is a
  Cardinal-Rule-1 stale ref regardless (points at a non-existent
  symbol). VERIFY absence before choosing the spelling: `hasattr(Cls,
  "planned_method")` / `.venv/bin/python -c "import ...; hasattr(...)"`
  → if False, write ``planned.method`` as a literal (honest: names the
  deliverable without claiming it links). Same `hasattr` gate distinguishes
  a LANDED seam (flip to a live `:meth:`, per L-007's landed-seam bullet)
  from a still-deferred one (literal). This is the FORWARD twin of the
  RETIRED-symbol rule above (that one repoints dead refs; this one
  refuses to mint premature ones).
- Packages that are not member-`automodule`'d (transport.*, the
  operator/scheme/numerics leaves on several pages) render their
  `:class:`/`:meth:` as plain text BY PAGE CONVENTION. Match the page —
  do NOT half-surface 1–2 leaves by adding an `automodule` while the
  rest of the package stays plain. The `:class:` ref staying plain-text
  is NOT a regression to fix — it is the convention; repoint a dead one
  to the LIVE module path (still plain-text, but now correct) and move on.
- **`:noindex:`-automodule'd is xref-invisible too — a WHOLE package can
  be plain-text page-wide even though it IS `automodule`'d.** An api page
  that `automodule`s every module of a package with `:noindex:`
  (`docs/api/diffusion_1d.rst` does this for all of `orpheus.diffusion.*`)
  registers NO cross-reference targets, so EVERY `:class:`/`:func:` to
  that package renders plain-text everywhere, while sibling packages
  automodule'd WITHOUT `:noindex:` (transport.method, numerics.eigenvalue,
  data.macro_xs.mixture) link normally. Diagnose by HTML link-audit
  (`grep 'href="[^"]*Symbol"' built.html` — empty ⇒ plain-text) + read the
  api page's automodule options; the module appearing in Sphinx's
  "highlighting module code" list means viewcode processed it, NOT that it
  has an xref target. This is NOT a defect to fix by editing the api page
  (often out of scope / forbidden) — keep the semantically-correct
  `:class:` markup (greppable, import-verified, auto-links if `:noindex:`
  is ever lifted), and FLAG the package-wide `:noindex:` as a candidate
  infra fix. The grep/import gate (symbol EXISTS) is the real cross-ref
  check; the link is governed by the untouchable api page. (Worked: #290
  P8 — all diffusion `:class:` refs plain-text via the api page's
  `:noindex:`; all 31 cited symbols import-resolved regardless.)
- **`automodule`-readiness is a MULTI-gate test; the 0-`:label:` check is
  necessary but NOT sufficient.** A leaf with NO `.. math:: :label:` is
  safe from the *duplicate-label* collision — but automodule RENDERS the
  whole docstring under the project's strict config, so it ALSO trips on
  any of: (a) a `:pydata:` (or other non-registered) role → `ERROR:
  Unknown interpreted text role`; (b) a section-underline shorter than its
  title inside the docstring → `WARNING: Title underline too short`; (c) a
  malformed field-list / inline-literal → docutils WARNINGs; (d) **a
  member attribute whose NAME collides globally** (a class `ng`/`n`/`name`
  attr surfaced by automodule makes EVERY pre-existing bare `:attr:\`ng\``
  across OTHER pages ambiguous → `more than one target for cross-reference`
  WARNING ×N, on pages you never touched); (e) **a malformed inline-role
  in ANY ONE rendered docstring** — the classic is a closing role-backtick
  immediately followed by a word char with no whitespace/punctuation
  (``:class:`X`s``) → `WARNING: Inline interpreted text ... start-string
  without end-string`. (e) fires even on a module that passes (a)–(d)
  cleanly (0 `:label:`, 0 `:pydata:`): ONE bad docstring line in ONE method
  blocks the whole `automodule`. The `-E -W` build is the only
  way to see (a)–(e) — a plain build with cached env MASKS them (one
  session: plain build EXIT 0 while `-E` showed 4 ERRORs + a 7-page `ng`
  cascade). A docstring WITH `:label:` must be cross-referenced in prose
  instead (automodule re-registers the label → duplicate-label / "equation
  not found").
- **Signature (e) intersects the "report-don't-edit-.py" constraint — an
  (e)-blocked module is automodule-UNREADY by the SAME rule as (a)–(d):
  cross-reference it in prose (plain-text, consistent with the
  un-automodule'd family) and REPORT the one-line docstring fix that
  unblocks the autodoc.** The plain-text refs are NOT a defect — they MATCH
  the page convention (grep the built HTML: the sibling family on the SAME
  page renders plain-text too). The "surface this type" intent is still met
  for refs pointing at an ALREADY-automodule'd symbol (those link), plus a
  prose bullet-list of the new types with their theory `:ref:` + a
  `.. note::` recording the exact docstring fix and the autodoc block to
  add post-fix. Verify with TWO gates `-W` is blind to: the grep gate
  (every cited symbol exists in live code) AND the rendered-HTML link audit
  (automodule'd→`<a>` link, un-automodule'd→plain `<code>` by convention).
  (Worked: P5 condensation #274 — `energy_grid.py` was (a)–(d)-clean but
  `GroupCondensation.from_grids` had ``:class:`EnergyGrid`s`` → reverted
  the `automodule::` to a prose bullet-list+ref+fix-note, build 3→0
  warnings; `Mixture.condense` refs link (Mixture IS automodule'd),
  `EnergyGrid`/`OverlapBasis` render plain-text like their `numerics.basis`
  siblings.) **SEQUEL — P5.5 reshape (same module family, one reshape
  later):** the (e) blocker `GroupCondensation.from_grids` was DELETED by
  the reshape, so `energy_grid.py` became automodule-ready — I surfaced it
  (anchors for `EnergyGrid.as_measure`/`.as_basis`/`.overlap_to` now LIVE,
  the L-002-deferred fix CLOSED). But a NEW blocker surfaced on the SIBLING
  `wims.py` — signature **(c)**: a module docstring's plain-text
  2-space-indented numbered list ("  1. …  2. …") after an unindented para
  → `ERROR: Unexpected indentation` + `WARNING: Block quote ends without a
  blank line`. Same resolution: automodule the clean pair
  (`energy_grid` + `ornl`), cross-ref `wims` in prose + `.. note::` the
  re-flow fix. LESSON: a retirement that deletes one automodule blocker can
  UNBLOCK a module — re-test automodule-readiness on the reshaped tree
  (the deferred-fix note may now be actionable); but a sibling in the SAME
  package can carry a DIFFERENT blocker class, so `-E -W`-build EACH
  automodule you add, never assume the package is uniformly ready.
- **The scoped resolution when a cohesive cluster is automodule-UNREADY:
  automodule ONLY the clean module(s), cross-reference the rest to the
  theory page** (the `api/numerics.rst` operator-family pattern: it
  automodules `field` but cross-refs `operator` because `operator.py` has
  `:label:` docstrings). Surfacing the whole cluster (fixing `:pydata:`,
  the `ng` collision via `:noindex:`, underlines) is a SEPARATE
  architectural docs task — DEFER it, do not block the carve-doc on it.
  (Worked: the dyad-carve task — `reaction_rate_functional.py` was clean
  (0 `:pydata:`, 0 `:label:`) → automodule'd to render both functional
  classes; `fission`/`scattering`/`multiplication` used `:pydata:` + the
  `ScatteringOperator.ng` global collision → kept as theory cross-refs.)

How to apply: treat the grep gate (symbol exists / page-convention
matched) as the real cross-ref check; the warning count proves only
"added nothing new". Before adding ANY `automodule`, `-E -W` build it in
isolation and watch for the (a)–(d) signatures above — especially the
cross-page `ng`-style cascade, which fires on pages you did not edit.
(The forced-`-E` build, the three-severity grep, and the venv/worktree
facts live in AGENT.md — do not re-derive them.)

---

## L-003 — GREP the V&V matrix for an eq-label BEFORE renaming or removing it

An equation `:label:` may be a `@pytest.mark.verifies(...)` TARGET. If
it is, the test's oracle points at that exact name — renaming or
deleting it breaks the verification edge (and bumps a "no matching
equation node" line). The recurring mistake is rewriting a stale
equation's BODY and changing its LABEL in one motion, orphaning a
verifies-target.

- For a stale equation that IS a verifies-target: KEEP the label name,
  rewrite only the BODY (the claim is unchanged). Split a busy
  derivation into a label-preserving primary + NEW sub-labels for the
  decomposed steps.
- A label-rename ripples: an in-page symbol rename (`s_a → c_a`) flows
  to every `:eq:` and prose site that referenced it — a whole-page
  sweep, not a one-line edit.
- Section-label vs equation-label are different namespaces. `.. math::
  :label: X` → `id="equation-X"`, resolved by `:eq:`. A `.. _X:` anchor
  → `id="X"`, resolved by `:ref:`. When you need a section anchor and
  an equation sharing a name, suffix the section `X-section`. A
  `:label:` inside a CODE docstring is rendered by autodoc but is NOT a
  global `:eq:` target — cite the owning `:class:` and inline the math.
  COROLLARY (caught here): do NOT write `:ref:`X`` where `X` is an
  EQUATION label — it renders plain-text cross-doc / warns intra-doc.
  Point at the equation with `:eq:`X`` and describe the surrounding
  prose ("the note under the production matrix :eq:`fission-matrix`").
- **When an ALGORITHM is replaced (not just an equation rewritten), a
  retired-STEP verifies-target is usually KEPT-AND-REPOINTED to a
  conceptual survivor of the NEW algorithm — not retired.** The
  recurring trap is reflexively retiring a power-iteration step label
  (`fission-source`/`fixed-source-solve`/`keff-update`) because "the
  iteration is gone", which orphans 4–5 test edges and forces test-side
  edits. Instead, ask whether the CONCEPT survives the new method: the
  per-iteration fission source `Q_f=(χ/k)νΣf·φ` → the single dyad
  application `Fφ`; the fixed-source solve `Aφ=Q_f` → the loss-matrix
  solve `M=A⁻¹F`; the production/absorption ratio `k=prod/absn` → the
  eigenvalue extraction `k=λ_max(M)`. Each repointed equation is a REAL
  step of the direct method, the `:label:` name is preserved, and the
  reconciliation table reports "kept-and-repointed → NO test-side edit".
  A `.. note::` under each repointed label states what it "historically
  named" and what direct analogue it now carries (preserves the WHY per
  L-007 without a tombstone — the equation evolved, it wasn't falsified).
  Only RETIRE a step label when the concept genuinely has no survivor
  (here `convergence-rate`: a direct dense eig has no iteration whose
  rate a dominance ratio governs — and it was documented-only, NO test
  edge, so retiring just drops one auto-regen row). (Worked: #276
  homogeneous refounding — power iteration → direct `λ_max(A⁻¹F)`;
  4 step labels kept-and-repointed, `convergence-rate` retired.)

How to apply: grep the generated matrix (`docs/.../matrix.rst` or the
audit output) for every label you intend to touch FIRST; preserve
verifies-targets by name (kept-and-repoint when the concept survives;
retire only documented-only labels with no survivor). Also `grep
':label:'` repo-wide before MINTING a new label — duplicate labels
across files collide (a real warning), and a new label can collide with
a same-named partition predicate already living in another page.

---

## L-004 — Representational/structural eq-labels get a `.. vv-status: <label> documented` DIRECTIVE, not prose

A NEW equation that is a field-typing identity, a governing iteration,
a literature-transcribed definition, or a derivation-decomposition step
is NOT a solver claim. It must be tagged so the V&V matrix files it
under "Documented-only", not flagged as an unverified solver claim. The
recurring mistakes: (a) writing the status as prose instead of the
machine-read DIRECTIVE (a `--strict` audit then regresses); (b) leaving
a label that a NEW test's `verifies(...)` points at as an untagged
orphan.

- Structural/representational label → `.. vv-status: <label>
  documented` with a one-line rationale comment naming the bit-identity
  / foundation gate that pins the verifiable content.
- A label a NEW test `@pytest.mark.verifies(...)` points at is
  `implemented` (code + test, no eigenvalue/flux claim). If that test
  verifies a label that does not exist yet, CREATE the label in the new
  section — never leave a verifies-target orphaned.
- Pure derivation-decomposition labels (the affine-cell sub-steps, the
  facewise-separable tensor identity) sit untagged in a verified
  narrative chain by established page convention — match the page's
  siblings rather than inventing a status.

How to apply: for every eq-label you add, classify it (solver-claim /
representational / verifies-target) and apply the matching status
discipline. This is the V&V-vocabulary-curation duty (Directive 5) at
the label level — the matrix is the audit's source of truth.

---

## L-005 — A stub→rich expansion reads its sources in a fixed priority order; the docstrings are the prose SEED

Expanding a `.. todo:: Archivist` stub (the algebra-of-record handoff)
into rich narrative has a load-bearing source-reading order, and the
recurring mistake is writing from the brief alone:

1. **The close-out / verdict memo** — carries the bug ("confirmed live
   pre-fix"), the architecture-settled framing, the named tests, the
   verification numbers, the honest-interim state.
2. **The production docstrings** (the scheme/operator/class bodies, the
   `supports()` predicate, the SymPy `derive_*` docstrings) — the
   VERBATIM prose seed: the numerical-PDE statement, the contrast notes,
   the per-case table, the lit cites. These are the algebra-of-record.
3. **The test files** — the verification subsection: the named gates,
   what each asserts, the bit-identity-vs-principled-equivalence split.
4. **The SymPy module** (when present) — NEVER expand a stub without
   reading it; the narrative narrates it, does not compete with it. If
   you find an algebra error, return a DISPATCH_REQUEST for the
   method-implementer — never edit the SymPy yourself.

How to apply: read memo → docstrings → tests → SymPy, in that order,
before drafting. The honest scope (what shipped vs what is OWED to a
follow-on slice) comes from the memo's interim-state note — preserve it
verbatim; do not over-claim a wired-but-not-yet-iterating capability.

---

## L-006 — Cross-document duplicate citations are a real warning class; resolve, do not redefine

A `.. [Key]` bib entry duplicated across two standalone theory pages
emits a duplicate-citation warning. The recurring mistake is adding a
fresh `.. [Key]` block on a new page for a reference already defined
elsewhere.

- Before citing, `grep '^\.\. \[Key\]'` — if the entry exists on
  another page, cite it cross-doc (resolve), do NOT redefine.
- Where a page's existing convention already dodges the collision by
  citing a reference as PLAIN TEXT in the Literature list-table
  (because the `.. [Key]_` form would cross-doc-collide), MATCH that
  convention — add new literature as plain text too, no new bib entry.
- Pre-existing duplicate-citation warnings (cross-document cite
  collisions) are a known trade-off for standalone pages and do NOT
  need elimination during a close-out — verify the COUNT is unchanged
  pre/post, not that they are gone.
- FLAG (don't silently use) a conflated bib key when you spot one (e.g.
  a key whose title is one paper but whose cited equation content is
  another) — that is a Cardinal-Rule-1 correctness defect for the
  method-implementer/literature-researcher, surfaced from your
  cross-page vantage.

---

## L-007 — A retirement/relocation doc preserves the WHY and tombstones the wrong claim; it never deletes evidence

When an issue closes by FALSIFICATION (the approach cannot work) or a
type/decomposition is retired but its CONCEPT survives, the close-out's
archival value is HIGHER than a success story — it stops future
sessions re-attempting a dead path. The recurring mistakes are
rewriting history and deleting numerical evidence.

- **Preserve the motivation that LED to the investigation** — flip
  tenses ("is expected to" → "was expected to") but keep the logic. A
  future reader asking "why did anyone try this?" must find the answer.
- **Tombstone, don't delete.** When a new finding invalidates a
  published table/claim on the same page, add a `.. note:: **Retraction
  (date, Issue #N).**` immediately above it: (a) what the claim was, (b)
  why it's wrong (one line), (c) forward-pointer to the new section.
  Numerical VALUES stay; the INTERPRETATION gets the tombstone.
- **Retitle to the concept, KEEP the anchor.** When a type is retired
  but the concept survives in new realizations, retitle the section to
  the concept (not the dead type), KEEP the section anchor (cross-doc
  `:ref:`s auto-pick up the new title), and add a prominent succession
  note naming the new realizations. De-role dead module-path
  `:class:`/`:meth:` to literals (grep gate per L-002, not `-W`),
  past-tense the type, repoint present-tense claims to the realizations.
- The full 9-step CLOSED post-mortem arc + the PARTIAL/OPEN variant +
  the multi-issue audit-table pattern live in AGENT.md "Close-Out
  Narrative Arc" — follow it; don't re-derive it here.
- **A surgical RENAME/re-point brief routinely uncovers adjacent
  SUBSTANTIVE staleness — repoint the dead ref (Cardinal Rule 1), but
  FLAG, don't silently rewrite, the surrounding behavioral-claim
  staleness.** A "re-point `X`→`Y` after `X` was retired" pass crosses
  sections that describe a now-superseded ARCHITECTURE (an intermediate
  collapse state, a since-reversed `domain=None` design, a retired
  factory/test still cited). The dead `:class:`/:meth:` ref to `X` MUST
  be repointed to the live `Y` regardless (it renders plain-text, no
  `-W` warning — L-002 — so the grep gate is the real check). But the
  surrounding PROSE staleness (a "thin subclass" narrative #261 fully
  dissolved; a "stays at `domain=None`" claim W-D reversed; a dead
  `:meth:`SNMesh.zeros_*`` / `:file:` test path) is a SEPARATE task: a
  behavioral-claim rewrite needs its own verify-against-live pass and
  often its own issue. Repoint-in-passing is correct; rewrite-in-passing
  risks minting a NEW false live claim (worse than a dead ref to true
  history). The brief's named file/line list is the scope FLOOR (the
  full grep is the blast radius for the RENAME); substantive-narrative
  fixes beyond the rename are the scope CEILING — report them as flagged
  findings, don't smuggle them in. (Worked: #261 CollisionOperator→
  MultiplicationOperator — the §5.7 "thin subclass" §, the
  "C/S/F stay at domain=None" §, and the whole PR-TYPED-3
  IsotropicSource/PerOrdinateSource § were all substantively stale
  beyond the rename; repointed the dead refs, flagged the three §§.)
- **The dedicated behavioral-rewrite follow-up (the task L-007's prior
  bullet hands off to): when a doc describes a FUTURE/DEFERRED seam and
  a commit since LANDED it, verify the SHAPE that shipped — don't just
  flip "deferred"→"done".** The realized change can close the seam via a
  DIFFERENT mechanism than the doc envisioned, and a naive flip mints a
  new false claim. Read the commit message AND the live code; tell the
  honest "envisioned X, shipped stronger/different Y" story. (Worked,
  #271: the "Deferred typing-completion seam" § envisioned giving `S` a
  BULK `V_bulk` domain so `OperatorSum` would REJECT an `S+B` fold;
  W-D actually gave `S` the COMPOSITE `full_field_space` — same instance
  as L/C/B — so the within-group `(L+C)-S` guard VALIDATES, and the
  once-envisioned "bulk-S ≠ trace-B" rejection NO LONGER applies at all.
  A flip-only rewrite would have claimed the rejection tripwire is now
  live — FALSE. The honest rewrite said the seam closed with a different,
  stronger choice and the fold stays gone for an UNRELATED reason (the
  variadic-driver redesign).) Corollary: when a reversal touches an
  adjoint/metric mechanism, separate the conclusion that SURVIVES (here:
  "the metric applies ONCE at the op level via `_AdjointOperator` reading
  the composite domain") from the stale PREMISE that justified it (here:
  "C carries domain=None so the metric propagates from L by
  first-non-None"). Preserve the conclusion, rewrite the premise — and
  re-derive WHY the conclusion still holds (the `_AdjointOperator` reads
  the SUM's domain, never per-leaf, so a leaf now carrying the composite
  domain is no double-application risk). Read the actual wrapper/composer
  source (`operator.py`) to ground the new "why", never reason from the
  old prose.

---

## L-008 — Generated artifacts are NEVER hand-edited; materialize them, report the REAL number

Several Sphinx inputs are generated on every build by a `builder-inited`
hook (the V&V matrix via `generate_matrix`, the capability matrices via
`generate_capability_matrices`, the `docs/_generated/*.rst` includes via
`generate_rst`). The recurring mistakes are hand-editing the rendered
output and estimating a count instead of running the generator.

- A matrix/capability-table drift clears by running the generator (or
  just building) — never hand-edit the `.inc.rst`. Edit the registry-
  side metadata (`capability_rows()`) or the test that drives it.
- When a refactor changes how many rows the matrix auto-regenerates,
  REPORT the real post-regen number (it can differ from a brief's
  estimate — e.g. an auto-regen dropping 67→54 rows, not the est. −5).
  An auto-regen row delta is NOT a warning.
- In a FRESH worktree, missing generated artifacts (`.. plot::`
  FileNotFoundError on `*.h5`; `.. include::` CRITICAL on
  `_generated/`) are ENV gaps, NOT doc defects — materialize them
  (run the converter / generator, both write only gitignored dirs,
  confirm `git check-ignore` first). Do NOT "fix" the docs to route
  around a missing artifact — that corrupts correct documentation.
- Hard-coded test counts in RST go stale; `pytest --collect-only -q`
  for the current number rather than transcribing one.

---

## L-009 — Section-marker ladder is FILE-LOCAL; underline length is in CODE POINTS; reuse the file's existing same-depth marker

Recurring build-breakers from title machinery:

- The `=`/`-`/`~`/`^`/`'`/`"` underline ladder is per-FILE. Scan the
  file's first-appearance markers (grep the single-char underline rows,
  tally by depth) before picking a level, or you get "Inconsistent
  title style: skip from level N to N+2" — often at sections you did
  NOT touch (introducing a marker at a SHALLOWER depth than the file's
  existing one for that level pushes the old sections down a level).
- Underline length is measured in CODE POINTS, not bytes. An em-dash
  `—` or `÷`/`χ` is 1 code point but multiple bytes — size the
  underline with `len(title)` in Python, never `wc -c`. A global
  normalize-pass that touches PRE-EXISTING tolerated over-runs must
  RESTORE them (scope edits to your own lines).
- When adding a section at an existing depth, REUSE the file's existing
  same-depth marker character (grep for it). Introducing a different
  marker for the same logical depth is the classic "skip level" trigger.

How to apply: map the file's marker ladder first; size underlines with
`len()`; reuse existing depth markers. (The catalog of which file uses
which ladder is in AGENT.md — this lesson is the GENERAL discipline.)

---

## L-010 — V&V vocabulary in prose MUST match the skill verbatim; you are the curator (Directive 5)

You write the prose future readers quote when reasoning about
verification status. The recurring failure is paraphrasing a level
definition or over-claiming what a reference proves. The hard rules
(from `vv-principles`):

- **MMS does NOT verify eigenvalues** — it is source-driven, reaches
  flux-shape / convergence-order only. NEVER write "MMS verifies the
  eigenvalue". Eigenvalue claims need closed-form (`k_∞ = νΣ_f/Σ_a`) or
  semi-analytical references; SI≡Krylov twin agreement is
  necessary-but-NOT-sufficient (needs a structurally-independent leg).
- **L4 (code-to-code) proves zero correctness** — every L4 claim names
  its L0–L2 backing. NEVER "L4 proves correctness".
- **1-group eigenvalue is degenerate** (`k = νΣ_f/Σ_a` flux-shape
  independent) — NEVER "the 1-group test verifies the solver". (1G IS
  fine for a rate/convergence-order claim — declare the claim layer.)
- **Name the pillar** (closed-form / MMS / semi-analytical), not vaguely
  "analytical"; respect each pillar's evidence boundary.
- A **Mode-10 sub-floor term** is closed by STRUCTURAL teeth
  (producer-threads-at-machine-precision + consumed-sign-flip-moves-flux
  + a no-op control leg), NOT a tightened value band — and when no
  isolating regime exists (a boundary-trace slope below the bulk floor
  everywhere), say so: there is NO value-improvement leg to add. A
  prophylactic `.. warning::` ("do NOT write 'S9 recovers 2nd order at
  the boundary'") in the doc itself pre-empts the future over-claim.
- Bug attribution cites the failure mode (1–11) and matches
  `error_catalog.md`; a `catches(ERR-NNN)` claim is a coverage edge,
  not a topic tag.
- **"Adjoint" operator prose: distinguish the EUCLIDEAN transpose (`Aᵀ`, the
  plain group/angle matvec adjoint under the L² inner product) from the metric
  HILBERT adjoint (`A†=G⁻¹AᵀG`, the `.H` wrapper).** A campaign/commit may name
  a new capability "S†" colloquially while the method actually computes the
  Euclidean transpose `Sᵀ` (pinned by reciprocity `⟨Sψ,χ⟩=⟨ψ,Sᵀχ⟩`, NO angular
  Gram). Write the precise object (transpose), note the metric adjoint is the
  separate `.H`, and say the campaign names it "†" colloquially. When the
  FORWARD keeps a fast-path while the ADJOINT rides a different (frame) form,
  that structural asymmetry is what makes the reciprocity gate a GENUINE
  cross-check (two structurally-different representations of one operator), NOT
  a tautology — say so explicitly; it is the structural-independence argument
  at the operator-identity level.
- **The FIRST iterative member of an inverse/solver family that was previously
  all-EXACT has NO bit-id twin to inherit — its claims rest on
  structural-independence anchors ONLY, and it carries NEITHER an eigenvalue
  NOR an MMS claim.** Documenting such a member (e.g. `GreenOperator`, the
  first iterative inverse in the #226 family): state the claim layer as
  foundation / flux-shape against a structurally-independent reference
  (closed-form dense-LU + the multiple-scattering / Neumann expansion), NOT
  inheritance. An iterative *sum inverse* is neither an eigenvalue solver (so
  NO eigenvalue claim) nor source-driven (so NO MMS reference — both pillars
  are INAPPLICABLE, not merely unused). And when the parent's `solve` is
  DEFINED as `self.inverse().apply` (the `OperatorSum` contract), the
  `inverse().apply ≡ solve` equivalence that anchored the family's EXACT
  members is a TAUTOLOGY for this member — exclude it as evidence EXPLICITLY,
  so a later reader doesn't add a green tautology and mistake it for coverage.
  The name-earning distinguishing invariant (here G-Neumann: the splitting a
  generic `A⁻¹` cannot satisfy) is the load-bearing correctness evidence, not
  round-trip.

**Skill-uplift duty:** when you notice a recurring published-prose
anti-pattern, a new failure-mode signature in a close-out, or a new
pillar-evidence-boundary case that `vv-principles`/`error_catalog.md`/
`algebra-of-record` does NOT yet capture, PROPOSE the skill edit in your
return. You read across all close-outs — the skill grows when you feed
it back. (Do not duplicate skill content into this file; point to it.)

---

## L-011 — An overloaded-symbol convention sweep: inventory EVERY meaning, classify each site, replace_all only unambiguous strings

A project-wide symbol-convention change (here: `L` the invertible loss
composite → `A = L + C`, because `L` is reserved for the streaming
leaf) is NOT a blanket rename — the letter is OVERLOADED, and a naive
find-replace corrupts the legitimate uses.

- **Inventory every meaning of the letter FIRST, then classify each
  site mathematically.** `L` meant five different things across two SN
  theory pages: streaming leaf (`\Omega\cdot\nabla`, `L+C`) → KEEP;
  loss composite (`(L,S,F)` triple, `(L-S-F)\psi=q`, `L^{-1}`/`L.solve`
  the sweep, `L.apply`) → FIX to `A`; Legendre order (`L=0` isotropic,
  `L=1`) → KEEP; slab length (`L=5 cm`) → KEEP; generic operator in
  the ABSTRACT-algebra sections (protocol `apply: x\mapsto Lx`,
  `ScaledOperator (\alpha L)`, closure `L=A+0`) → KEEP + FLAG (it does
  NOT denote the composite; swapping it to a neutral symbol is a
  separate style task, out of scope for an L→A fix).
- **VERIFY the target spelling against LIVE code before committing to
  the scope** (L-001 applied to scoping). When the sweep felt alarmingly
  large (~40 `(L,S,F)` sites), reading `iteration.py` signatures
  (`SourceIteration(A_inv)`, `KEigenvalue`'s "`A` the FORWARD invertible
  loss operator", the "(A,S,F) operator triple" comment) proved the
  docs' `(L,S,F)` was STALE vs code — so the comprehensive sweep was a
  Cardinal-Rule-1 doc/code alignment, NOT over-reach. Live signatures
  resolve "is this pervasive spelling stale, or intended?".
- **replace_all ONLY on unambiguous multi-char strings** (`(L, S, F)`,
  `(L - S - F)`, `(L-S)^{-1}`); targeted context-edits for ambiguous
  bare `L`/`L^{-1}`/`L.apply`. Enumerate spacing/punctuation variants:
  `(L - S)` (spaces) MISSES `(L-S)` (no space); `(L, S, F)` MISSES
  `(L, S, F, \psi)` (trailing arg). A FINAL grep of ALL remaining
  `L`-forms — each explicitly re-classified KEEP/FIX — catches the
  variant misses (found 3 stragglers that way).
- **A convention-target symbol can COLLIDE with an existing same-letter
  use in ONE section.** Here `A` (loss operator) collided with `A` (the
  affine SPACE in the #208 torsor section) — but that section already
  used operator-`A` in `r = A\psi - q`, so it was already internally
  inconsistent with its own "connected by `L`" prose. Resolution: use
  the convention's symbol (it was already there), disambiguate AT
  INTRODUCTION by spelling the definition (`the loss operator
  A = L + C`), and FLAG the pre-existing space/operator collision for a
  future disambiguation pass — do NOT rename the other use (out of
  scope), do NOT abandon the convention.
- **Eq-label BODIES change, labels + vv-status stay (L-003).** The
  sweep crossed verifies-target equations (`operator-fixed-source`,
  `sn-streaming-reciprocity`, `streaming-inverse-direct-sum`,
  `octant-direct-sum-tensor-product`): change `L`→`A` in the math,
  keep the label name — even when the label reads "streaming-inverse"
  but the body is now `A^{-1}` (the octant block-structure IS
  streaming-induced; say so in one clause so the label name still
  makes sense). The clean `-W` build confirms every `:eq:` still
  resolves; matrix.rst auto-regens (its only diff was a test-count
  bump from a CONCURRENT agent, not my labels — L-008).

How to apply: for any symbol-convention sweep, map every meaning of the
letter, verify the target against live code, replace_all only the
unambiguous strings, targeted-edit the rest, grep-sweep to re-classify
every survivor, and flag (don't fix) same-letter collisions + the
generic-operator-symbol residue.

---

## L-012 — Merging a RE-STAGED branch's authored docs into a DIVERGED tree: programmatic verbatim splice + path translation + same-merge forward-ref reconciliation

When a long-lived feature branch's DOCS commit can't be raw-applied
(both trees rewrote the pages since the fork), the job is to integrate
the branch's authored CONTENT into today's pages, not cherry-pick the
diff. The disciplines:

- **Extract the branch's `+` block PROGRAMMATICALLY; never hand-retype
  1000+ lines.** For a large authored block carrying verifies-target
  eq-labels + math that MUST be verbatim, hand-transcription is a
  Cardinal-Rule-1 risk. Slice the diff's `+` lines (`diff[start:end]`,
  assert each `startswith('+')`, strip index [1:]), apply path
  translations as string-replaces, write to a temp file, then splice
  via a single `src.count(target)==1`-asserted replace. VERIFY the
  extract before splicing: label counts (each verifies-target ==1),
  anchor counts, first/last lines, residual-old-path==0. A machine
  splice can't introduce a transcription error; an Edit new_string of
  1000 lines can.
- **The branch was authored PRE-reorg — translate EVERY module path to
  the LIVE layout.** Inventory the paths in the diff
  (`grep -oE 'orpheus\.sn\.[A-Za-z._]*'`), map each against the reorg
  (here `orpheus.sn.geometry` → `orpheus.sn.mesh.augmented_mesh`), and
  confirm one prefix-replace covers all occurrences with **zero
  residual**. Verify each translated symbol exists in LIVE code (the
  grep gate — a dead xref renders plain-text, no `-W` warning, L-002).
  Distinguish the moved package (`orpheus.sn.geometry` → `mesh.*`) from
  a same-named-but-UNMOVED one (`orpheus.geometry.reduced_operator`,
  top-level, unchanged) — a blind replace corrupts the latter.
- **Find insertion points by CONTENT/anchor, NEVER the diff's line
  numbers.** The `-` hunk numbers are fork-relative; the tree moved
  (here a section shifted fork-4288 → today-5487). Match the exact
  unchanged target string (assert single match), and MAP the file-local
  marker ladder at the insertion point (the enclosing-section chain,
  L-009) to confirm the block's `~~~~`/`^^^^` levels are valid children
  there — the branch's markers happened to match today's ladder, but
  that is VERIFIED, not assumed.
- **A forward-reference to an issue that LANDED in the SAME merge is now
  stale — flip tense + cross-link (the L-007 landed-seam bullet, sharp
  same-merge instance).** The re-stage bundled Phase-2 content (which
  forward-references #248 as *pending*: "survives Step C … retirement
  tracked under #248") WITH #248 itself (which retired it). Integrating
  both makes the pending-tense a Cardinal-Rule-1 self-contradiction.
  Verify the SHAPE that shipped (`grep 'def __call__'` → fully gone from
  live code), preserve the WHY (why Step C deliberately scoped it out),
  flip "is tracked / survives" → "was subsequently retired under #248",
  and cross-link to the #248 note.
- **Retirement-audit: separate dead CROSS-REFS from stale LITERALS
  before acting.** The requirement's "unresolved-xref hazard" targets
  `:class:`/:meth:`/:attr:` refs (render plain-text, no warning). Grep
  specifically for the XREF FORMS of the retired symbols
  (`:(attr|meth|class):\`[^\`]*\.retired_field\``); if ALL references
  are double-backtick literals (as here for `tau_mm`), there is NO
  build-invisible hazard — the remaining work is pure Cardinal-Rule-1
  literal-claim correction. Fix the literals that DIRECTLY contradict
  the just-integrated content (clean scope: "τ single-sourced in
  spherical_streaming" → in the closure); TOMBSTONE/FLAG the ones
  describing a retired MECHANISM, ESPECIALLY when its root staleness is
  a DIFFERENT issue than the one you're integrating (here the
  `alpha_in is None` slab/curvilinear discrimination was retired by an
  earlier Issue #196 Phase G Step 2.5, surfaced during a #236 audit →
  add a correction `.. note::` grounded verbatim on the live
  `StreamingTerms` docstring, flag for a dedicated rewrite; do NOT
  fabricate the current mechanism, per L-007).

How to apply: read the fork-diff as the CONTENT source, not a
patch; splice the large block programmatically + translate paths +
verify the extract; place by anchor; then run the standing
retirement-audit (repoint xrefs, reconcile same-merge forward-refs,
fix-or-tombstone stale literals) and the `-E -W` build-count gate.

---

## Quality self-assessment rubric (Directive 3)

Rate each output 1–5 and log the weakest dimension in the return:
Derivation depth · Cross-references · Numerical evidence · Failed
approaches · Code traceability · Derivation source (from `derivations/`,
never hand-written). The recurring weak dimension on TERMINOLOGY/ROUTING
fixes is "numerical evidence" — structurally absent (no flux moves → no
convergence table), not a deficit; say so rather than manufacturing one.

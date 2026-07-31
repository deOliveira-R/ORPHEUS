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
- **In a CLAIM-CLASSIFICATION correction sweep (classify each repeated site
  A/B/C, correct only the false class), the SAME phrase can be TRUE at one site
  and FALSE at another — the discriminator is a LIVE config detail, not the
  prose.** Worked (#280 Phase 2.5b, `discrete_ordinates.rst` + `tests/sn/`): the
  recurring claim "the cylinder α-dome telescopes the seed away / was already
  exact" is a FALSE mis-attribution for a **product** quadrature (the starting
  direction coincides with the first-swept ordinate #229, so
  `c_in[m0]=(1−τ)/τ ≠ 0` is a LIVE self-coupling on the m0 diagonal; the cold
  `(L+C).solve` was seed-lagged err≈0.57 until the direct-seed fold
  `c_out→c_out−c_in`) yet genuinely TRUE for a **level-symmetric** one
  (`c_in[m0]=0` at raw τ=1 — a DEAD first-ordinate weight annihilating the seed
  at source, NOT telescoping). So a "cylinder was already exact" docstring is
  only classifiable after grepping the site's LIVE quadrature: a
  `Quadrature.level_symmetric(4)` fixture makes it TRUE (LEAVE — e.g.
  `test_282_direct_seed_fixed_point._operator`), a `product` fixture makes it
  FALSE (CORRECT). Second discriminator, SAME word, different class: "α-dome
  telescopes" splits (A) seed-absorption-of-the-SOLVE (false, level-symmetric-
  only) vs (B) weight-summed-scalar-gate-blindness (anti-pattern #8, TRUE, LEAVE
  — the `Σ_n w_n(α_{n+½}−α_{n−½})=0` identity) by WHAT it claims telescopes —
  the seed's effect on the fixed point (A) vs a per-ordinate error's effect on a
  scalar residual (B). Disambiguate by the OBJECT, not the verb. Also: the fix's
  new changelog bullet + the reframed load-bearing note reference #280 Phase
  2.5b; leave the main agent's already-corrected model sites untouched (grep the
  issue tag to find them first).

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
- **PHANTOM verifies-target (the INVERSE of an orphan): a `verifies("X")`
  marker whose `X` has NO `.. math:: :label:` anywhere under `docs/`.**
  The audit flags it (`tests/_harness/audit.py` `_phantom_verifies`) — a
  silently-dropped V&V edge. Fix = EITHER add `:label: X` to the equation
  the test verifies (only if that equation is UNLABELED) OR repoint the
  marker to the real label the equation ALREADY carries. Decide by reading
  what each test verifies + grepping for the equation's existing label:
  when the topic's equations are ALL already labeled (here the LD slab is
  `ld-ubld-d1-reduction` (operator) + `ld-ubld-slope-angular-reduction`
  (the ERR-061 reduction)), the "add a label" branch is impossible (one
  `.. math::` = one label) → REPOINT, per-test, to the accurate label
  (convergence test → the operator; diffusion-limit/frame tests → the
  reduction). BONUS: if the repoint targets were FORMER ORPHANS (labeled,
  0 tests), the repoint kills the phantom AND covers the orphans in one
  move — a net V&V-hygiene win that VALIDATES repoint over add-a-label or
  delete-the-marker. Watch for a SIBLING phantom synonym (here
  `ld-cartesian-1d`, the un-homed 1-D umbrella parallel to the real
  `ld-cartesian-2d`) co-declared on the same tests: if it has no free
  equation home (the natural one is taken by a verifies-target you can't
  rename per this lesson), it is OUT of a single-label task's scope —
  FLAG it, don't half-fix. (Worked: stencil-assembly 2b — `ld-slab`
  4-marker phantom repointed; `ld-cartesian-1d` flagged.) **SEQUEL
  (Task #55): the flagged `ld-cartesian-1d` was RESOLVED by repoint →
  `ld-ubld-d1-reduction`** (the 1-D LD operator equation ALREADY labels
  the natural home; one eq = one label ⇒ mint impossible ⇒ repoint, and
  dedup the mark that already carried it). The SAME task took the OTHER
  branch for a sibling phantom `inverse-as-operator`: **MINT**, because
  its law (`A.inverse().apply ≡ A.solve`, the #226 keystone) was stated
  in PROSE but UNLABELED. Discriminator: repoint when the law's equation
  is ALREADY labeled; mint (with `.. vv-status: <label> documented` for a
  foundation/structural law, L-004, matching a sibling documented label's
  style) when the law is stated but carries no `.. math:: :label:` yet.
  Verify the mint lands in the matrix's "Documented-only" bucket.

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

**Orphan-slice adjudication triage (the V7 backfill task, #231 #10).**
When handed a slice of EXISTING orphan labels to sort into SENTINEL /
WIRE / GAP, three discriminators do almost all the work, in order:

1. **SIBLING-CONSISTENCY is the dominant signal.** On the foundation
   operator-algebra pages every label family already has adjudicated
   members — `apply-distributes`/`solve-does-not-distribute` govern the
   whole apply-solve family; `eigen-standard-form` governs the eigen-*
   posing family; `g-adjoint-definition` the g-adjoint family;
   `wdd-forward-recurrence`/`…-three-terms` the WDD/tensor-network family;
   `green-neumann-series`/`matrix-functor-out` the inverse family. An
   orphan that is a mathematical identity / posing statement / structural
   decomposition / literature-transcribed bound in the SAME section as a
   SENTINELED sibling → SENTINEL, with the same rationale shape. This
   made ~30/31 labels SENTINEL in one batch.
2. **A doc that says its pins are `@pytest.mark.foundation` "no
   verifies() by design" GOVERNS its whole label family toward
   SENTINEL, never WIRE.** foundation tests never carry `verifies()`
   (vv taxonomy) — so a defining identity (e.g. g-adjoint reciprocity,
   pinned only by a foundation reciprocity/oracle suite) is a documented
   sentinel even though a test "exercises" it. The page's own explicit
   V&V framing ("the claim is a software/algebra invariant, anchored to
   the structurally-independent oracle") is the ruling, not your instinct
   to wire the named test. Read the test's `pytestmark`/docstring.
   **BUT `@foundation` ALONE does NOT imply SENTINEL** (correction from
   the REFERENCES-part batch, #231 #10 A5): the algebra-of-record
   reference pages (peierls / peierls_nystrom / trajectory_resolvent /
   fn_method / galerkin_spectral / singular_eigenfunction) pin their
   V_αN / V_cg / V_se SymPy-identity labels with tests marked
   `@foundation` AND `@verifies("<label>")` TOGETHER — and a
   verifies-mark on a foundation test DOES close the orphan
   (**audit-proven**: the sphere `peierls-greens-V-alpha-2` identity is
   "covered", NOT an orphan, though its ONLY gates are `@foundation`+
   `@verifies` in `test_peierls_greens_function_symbolic.py`; harness.rst's
   "foundation tests carry no verifies" is a general EXPECTATION the
   reference pages routinely break). So the real discriminator is the
   SIBLING's decorator, not the `@foundation` tag: (a) foundation-NO-verifies
   pins (operator-algebra oracle/reciprocity suites; V_se-cyl.N SymPy
   gates; the variant-α-core "one foundation test per primitive, no
   verifies" file) ⟹ SENTINEL; (b) foundation-WITH-verifies pins (the
   SymPy V_n derivation gates whose file MIXES `@foundation`+`@verifies`
   on the identity-establishing tests) ⟹ WIRE the orphan into the
   PARALLEL foundation+verifies test (e.g. the slab `T_00^slab=2E_3`
   V_α2_slab tests mirror the sphere V_α2 convention → I minted
   `@verifies("peierls-greens-slab-V-alpha-2")` on the closed-form +
   numerical + overall-pass tests, leaving the substitution-algebra tests
   alone since a 2E_3 sign flip doesn't red them).
3. **WIRE (or DEFERRED-WIRE) only when a NON-foundation test tightly
   pins THAT exact equation with sign-flip teeth.** The discriminator:
   `tests/conftest.py` records an untagged test as `level=None` (NOT
   foundation) — so a `np.testing.assert_array_equal` law-test
   (`test_sum_law`/`test_product_law`/`test_scaled_law` pinning
   `[A+B]=[A]+[B]` etc.) is wire-eligible, unlike an `@foundation`
   reciprocity gate. If that catcher lives in a test file OWNED by a
   concurrent agent (do-not-edit list), report it as DEFERRED WIRING
   with exact node ids — do NOT sentinel over it (that papers over a
   genuine gate; harness normative §vv-status forbids it). The task's
   own do-not-edit list ANTICIPATES the wire targets it names.
- **Placement is same-FILE-only (the audit enforces file, not
  position).** BEFORE-the-math `.. (vv-status rationale) …` +
  `.. vv-status: <label> documented` (blank line, then `.. math::`) is
  robust and audit-valid in ANY file: anchor the Edit on the unique
  `.. math::` + `:label: <label>` two-line string and prepend. Match a
  bullet-indented math block's indent (2-space comment, 5-space
  continuation). Self-check = a Python pass asserting every
  `vv-status: X` has a same-file `:label: X` (a typo'd sentinel is a
  hard exit-2 abort); then the permitted single end-of-run audit
  (exit 0 ⇒ no typo/misplacement) — IGNORE foreign orphans/violations
  from concurrent sibling batches.
- **De-freeze a live suite-total** ("**12 tests**:") by dropping ONLY
  the drifting total; KEEP the designed case-list breakdown ("(5) … (5)
  … (2)") — those are structural, not drifting. A foundation suite is
  not an equation-matrix row, so a `:doc:`/…/matrix`` pointer does not
  literally fit; the honest move is drop-total-keep-structure.
- **The "don't bulk-sentinel row-sum / T-matrix / escape-probability"
  caution resolves by the GATE'S MARK CLASS, not the topic.** A brief may
  flag those as WIRE-likely; check each gate individually. On the peierls
  pages the infinite-medium row-sum `Σ K_ij Σt=Σt` and the
  finite-cell-deficit are derivation-context identities whose TESTED
  realisation is a DIFFERENT label (the finite-cell
  `peierls-vacuum-bc-row-sum-gate`, already wired) ⟹ SENTINEL; the
  T-matrix rank-1 gates are `@foundation`-no-verifies ⟹ SENTINEL (matching
  the sphere sibling `peierls-specular-T-matrix`); ONLY escape-probability
  had a real l0-WITH-`verifies` gate (`TestSlabPescClosedForm`,
  "Factor-level verification — slab P_esc") ⟹ WIRE. Corollary: a test
  NAMED for a SIBLING label can pin the ORPHAN's distinguishing content —
  `test_g_prefactor_is_4_over_pi` (verifies `peierls-cyl-3d-mode-formula`,
  the P_esc twin) asserts the 4/π G-prefactor that is the ONLY thing
  distinguishing the orphan `peierls-cyl-3d-gbc-mode-formula` ⟹ add the
  orphan as a SECOND `verifies` target on that same test. Net for a
  continuous-reference-derivation slice: ~94 % SENTINEL is the CORRECT
  ratio (50/53 here), and 0 GAPs is legitimate — every label is either a
  governing/definitional/literature-identity (SENTINEL case a) or has a
  real verifies-eligible gate (WIRE); a derivation page rarely hides a
  load-bearing contract with NO test anywhere.

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

## L-013 — A fix that RETIRES a failed-approach family gets a SUCCESS-resolution chapter, not the 9-step CLOSED arc — and proportionate treatment of the big historical narrative it supersedes

When a landed fix works BY retiring a whole failed-strategy family that
a large historical arc built up (here #282 route (a) retired the
`PsiHalfAngleSeed` seed zoo threaded through a 2800-line Phase-D→F→ERR-058
saga), the doc is a **resolution chapter**, not a CLOSED post-mortem. It
still borrows the close-out arc's spine (status banner, what-was-tried-
and-failed table, numerical before/after) but its verdict is "the fix
shipped", not "the path is dead". The disciplines that made this
proportionate and correct:

- **New resolution SECTION as the saga's final chapter** (a `-` h2
  sibling of the close-outs it resolves, placed right before the next
  `=` h1), with a status banner, the structural defect (the back edge),
  the fix, the derivation, the **failed-strategy `.. list-table::`**
  (each retired class as a LITERAL + its failure mode), numerical
  evidence, and the honest-scope caveat. This is the primary deliverable.
- **ONE loud supersession banner at the ARC HEAD** (`.. attention::`
  right under the first superseded section's title) telling the reader
  every "current default / retained strategy" claim in the sections
  below is HISTORICAL, with the forward `:ref:`. This lets the whole
  historical narrative stay as legitimate history WITHOUT rewriting it —
  proportionate. Do NOT tombstone every stale sentence.
- **Targeted retraction tombstones ONLY on the bald factual REVERSALS**
  (an "X is retained (not deleted)" that route (a) deleted; an
  "Infrastructure retained" table listing the retired classes as
  production). `.. note:: **Retraction (date, Issue #N).**` per L-007 —
  the numerical/historical evidence stays; the interpretation gets the
  tombstone.
- **The prior close-out's "Open research paths" is a GOLDMINE — flip the
  one that LANDED.** If the fix implemented a predicted research path
  (here #1 "TRUE-source-driven sweep-side seed" + the full Legendre
  fold, AND its proposed "sweep quadrature order" probe was EXACTLY the
  N-sweep discriminator used), flip it "→ **LANDED as #N**", tell the
  "predicted exactly / one refinement over the prediction" story
  (L-007 landed-seam). This is a powerful principled-resolution
  validation AND closes the research-path loop on the same page.
- **Literalize the retired family's dead xrefs via `replace_all`**
  (L-002/L-011): each `:class:`X``/`:meth:`X.__call__``/`:attr:`X.a``
  is backtick-delimited and unambiguous → `` `X` `` literal naming the
  historical artifact (NOT the surviving successor — a retired STRATEGY
  class and a surviving free FUNCTION of a similar name are different
  concepts; name the class as history, point prose at the function).
  The `-W` build is BLIND to these (plain-text) — grep-gate is the check.
- **Reuse the topic's existing verifies-target labels** (`hebert-3-43x`
  were already the full-fold equations) via `:eq:`; mint new
  `documented`-status labels only for the fix's OWN structural identities
  (augmented-composite, block-triangular, the streaming-manufactures-
  anisotropy identity) — each gets `.. vv-status: <label> documented`
  (L-004) or it lands in the matrix's "Orphan equations" bucket, not
  "Documented-only". Verify the bucket in the regenerated `matrix.rst`.
- The L-012 programmatic-splice + placeholder-underline generator works
  for a large NEW authored section too (not just merging a branch diff)
  — it guarantees code-point-correct unicode underlines (L-009).

How to apply: for a "the fix worked by retiring the failed family"
doc, write the resolution chapter + ONE arc-head supersession banner +
targeted reversal tombstones + the landed-research-path flip; literalize
the dead family xrefs; reuse existing verifies-labels, tag new
structural ones `documented`. Preserve the history, make the
supersession loud, don't rewrite the saga.

---

## L-014 — Deepening a thorough resolution chapter with physics RULINGS + the current-truth vs PLANNED-design-direction split

Two distinct moves recur when a design session digs deeper physics /
architecture for a feature whose *representation* is ALREADY thoroughly
documented (here: the SN curvilinear ψ½ pole seed, whose #282 route-(a)
resolution chapter already covered the defect/fix/representation, and a
`facefield_codim1_design.md` note that dug the deeper physics + a
planned refactor).

- **Augment with the deeper WHY, cross-link the existing WHAT — never
  duplicate.** When the chapter already documents the representation
  (the composite, the block-triangular normal form, the role-quadruple),
  the augmentation's whole value is the PHYSICS beneath it: *why* the
  direct solve exists (the pole is a straight characteristic → pure 1-D
  ODE, `(1−μ²)=0` kills redistribution), *why* the presence-predicate
  looks as it does (topology of the redistribution axis: a periodic CIRCLE
  gives edge-inclusion free + spectral, an INTERVAL's open GL makes you pay
  a seed).  ⚠ **The "zero metric is structural / entirely-grazing angular
  face" ruling this bullet ORIGINALLY carried was REFUTED — see L-015.**
  The ψ½ block's Hilbert STATE metric is the SPD `G_sd = V_cell`, NOT zero;
  `(1−μ²)|_pole = 0` is an OPERATOR coefficient (through-flux M2), not the
  state metric (M3); keep the three measures distinct (M1 moment weight /
  M2 through-flux coefficient / M3 state metric) — only M2 is zero. Give each ruling
  its OWN labelled `~~~~` subsection (`documented` vv-status on any
  literature-transcribed derivation eq per L-004), and cross-link the
  existing representation eqs (`:eq:`…-block-triangular``) rather than
  re-deriving them. The "what was tried and declined" belongs here too —
  a numerically-affordable-but-architecturally-declined alternative (the
  Gauss–Lobatto pole-node study: ~1.2× penalty, but declined to keep the
  cell-centred bulk clean) is exactly the Cardinal-Rule-3 "tried and why
  not adopted" content; cite the scratch-artifact location, mark it
  uncommitted/promote-only-if-adopted.
- **A PLANNED refactor on a current-truth page gets a LOUD
  design-direction admonition + a paired "current state" subsection.**
  When the design note proposes a refactor NOT yet built (a `FaceField`
  codim-1 parent, a `face_streaming_normal` measure, mesh-derived
  presence), the current-truth page must NOT read as if it exists.
  Recipe: (1) a `.. admonition:: Design direction — … (PLANNED, not
  built)` `:class: note` opening with an explicit "not implemented — no
  X/Y/Z type exists yet" + a pointer to the design-note path; (2) unbuilt
  types are LITERALS (``FaceField``), never `:class:` refs (L-002
  forward-ref — a `:class:` to a non-existent symbol renders plain-text,
  no `-W` warning, but is a Cardinal-Rule-1 stale ref); an ASCII
  hierarchy tree in a `::` block sidesteps the issue (all literal); (3)
  a SEPARATE "Current state — the `Optional` block and the N guards"
  subsection stating what IS built, so plan and reality never blur.
  VERIFY the guard/DOF counts against LIVE code, not the design note — a
  scoped grep undercounts (a `_require_*` call site lived in a sibling
  module `operators/streaming.py`, outside the `loss_representation/`
  scope I first grepped; the note's "3+4=7" was right, my scoped grep
  said "2+4"). Verify the enclosing CLASS of a cited call site too (the
  streaming.py site was in `InvertibleOperator`, not `StreamingOperator`).
- **The same-topic forward-ref in an ADJACENT subsystem is usually
  already stale — flag with the shipped evidence, don't rewrite (L-007).**
  A "when #N lands, gate X flips" future-tense claim in a neighbouring
  section (here the assembly-mode "Cartesian-only scope" §) is a
  landed-seam once #N ships; the TEST comment confirms the flip
  (`test_assembly_mode.py`: "route (a) makes the augmented matrix EXACTLY
  block-lower-triangular … replaces, not relaxes, the RED
  characterization"). But if it's a behavioral-claim rewrite in a
  different subsystem (assembly, not the pole-seed physics), FLAG it with
  the verified shipped shape + exact lines for the main agent — a rewrite
  needs its own verify-against-live pass and may interact with a claim
  that DIDN'T change (the bulk assembler staying Cartesian-only).

How to apply: for a "deepen an already-documented feature" task, add the
physics rulings as labelled subsections that cross-link (not duplicate)
the existing representation; give any planned refactor a loud
PLANNED-design-direction admonition (literals for unbuilt types) paired
with a current-state subsection; verify every count/class against live
code; flag adjacent same-topic landed-seam staleness rather than
smuggling a behavioral rewrite into scope.

---

## L-015 — The SUCCESS-CORRECTION doc pass: a doc's own PHYSICS FRAMING was the bug, proven this session; rewrite every site to the corrected story

Distinct from L-013 (retiring a failed APPROACH) and L-014 (deepening a
resolution chapter): here a framing DOCUMENTED AS CORRECT by prior
sessions — even by THIS lessons file (L-014's since-corrected "zero metric
is structural" clause) — was PROVEN WRONG this session and the code fix
landed (SN ψ½ block metric: the retired "ghost metric" `G_sd ≡ 0` → the
SPD state metric `G_sd = V_cell`). The doc job is to make every site tell
the corrected story. The disciplines:

- **Blast radius = the refuted CONCEPT grepped tree-wide, NOT the brief's
  file list.** The brief named 6 sites; `grep "ghost metric\|G_sd = 0\|
  all-zero ghost"` across `orpheus/` + `docs/theory/` found 3 MORE
  (transport-layer field docstrings) carrying the SAME refuted claim. Fix
  them (clean one-line Cardinal-Rule-1 corrections) and FLAG the
  scope-expansion. EXCLUDE frozen archaeology (`.claude/plans/*`, other
  agents' `agent-memory/*`, sibling worktrees) — they keep the old framing
  as history.
- **RENAME (not keep) an anchor whose NAME encodes the REFUTED concept**,
  updating every referencing site — the INVERSE of L-007's
  keep-the-anchor-when-the-concept-SURVIVES. `sn-282-ghost-metric-face` →
  `sn-282-pole-state-metric` because the section now REFUTES the ghost
  framing. Safe only because all inbound refs (here 2, both cross-doc in
  the page I was already rewriting) were updated in the same pass; VERIFY
  in built HTML that the new `id=` exists, the old is gone, and each
  cross-doc `:ref:` resolved to an `<a href>` (section anchors DO resolve
  cross-doc, unlike eq-labels — L-002).
- **Preserve the retired-bug WHY in PAST TENSE (L-007), don't erase.**
  Every `G_sd = 0` / "ghost metric" that SURVIVES the sweep must read as
  history ("the retired ghost bug installed…", "the shipped `G_sd = 0` WAS
  a wrong adjoint", "a reverted ghost would ALSO leave a defect"). A final
  grep re-classifies each survivor KEEP-as-history vs FIX; the ones in the
  fix's OWN validation body (the error message explaining why `0` is
  REJECTED) are correct as-is.
- **The load-bearing content is the CATEGORY-ERROR framing.** When N
  quantities were conflated into one wrong value, the crisp correction is
  a per-quantity `.. list-table::` (M1 moment weight / M2 operator
  coefficient / M3 state metric — each "where it lives / what it governs")
  + the WHY only ONE is zero: an operator coefficient equals the state
  metric ONLY when the face's operator self-block is trivial (spatial
  trace `A_tt` = restriction map → through-flux = state metric; ψ½ `A_ss` =
  banded radial transport op → through-flux 0 ≠ state metric `V_cell`).
  Ground "trivial self-block" in MEASURED norms from the
  derivation-of-record (`A_tt` offdiag ≈2 vs `A_ss` ≈71, ratio ≈35×).
- **A refuted-framing fix that CLOSES a vv failure mode gets a POSITIVE
  reframe naming the mechanism.** The old Mode-12 gotcha ("G-recip is
  IDENTICALLY blind to the seed rows") → "G-recip CATCHES a seed-row error
  (Mode 12 CLOSED, ERR-067)". NEW closure mechanism (Directive-5 skill
  proposal, below): a Mode-12 blindness closes EITHER by gating the OBJECT
  (the skill's canonical remedy) OR by REPAIRING the functional's METRIC so
  the error class LEAVES the invariance group — available exactly when the
  metric WAS the bug (correctness fix ≡ Mode-12 closure). Cite the LIVE
  gate name (grep it — the test was RENAMED `…_is_blind_to_…` →
  `…_catches_…`) and its both-legs subtlety (control leg: unmutated recip
  holds `<1e-12`; mutated leg: flip reds `>1e-6` — without the control a
  broken baseline mimics "caught"). ERR-067 lands in `error_catalog.md`
  from a concurrent QA agent — verify it exists + your wording matches
  before citing (it did).

How to apply: grep the refuted CONCEPT tree-wide (blast radius > brief);
rename refuted-concept anchors + update all refs (verify in HTML);
preserve the bug's WHY in past tense; correct via a per-quantity
category-error table grounded in measured norms; give any closed vv
failure mode a positive reframe naming the closure mechanism; and SHARPEN
any prior lesson (L-014 here) that carried the now-refuted framing.

---

## L-016 — The EVICTION/re-homing doc-pass: a sub-object leaves a block-ON-a-carrier for its OWN coupled composite — the PHYSICS survives, so reframe the CARRIER narrowly

Distinct from L-007 (retirement) / L-013 (retiring a failed approach) /
L-015 (a refuted framing): here a correct sub-object (the ψ½ ray) is
RE-HOMED — it moves from an optional third block ON a carrier
(`FullField`) to its OWN 2-block composite (`RadialCharacteristicComposite`
= System B) coupled to System A via a `CoupledField[ψ_A, ψ_B]` pair. The
governing insight: **an eviction changes the CARRIER, not the PHYSICS**,
so the doc job is a NARROW carrier reframe, not a chapter rewrite.

- **The brief's "~N stale refs" over-counts because most of the
  "3-block" chapter is PHYSICS that survived.** The route-(a) resolution
  chapter (pole straight-characteristic, M1/M2/M3 metric, block-triangular
  walk order, source fold, R12a, circle-vs-interval) describes System B's
  physics and is UNCHANGED by the eviction. Grep for the STALE CARRIER
  FRAMING ("grows a third summand", "augmented (bulk ⊕ trace ⊕ seed)
  composite", "third (seed) block", "optional third block",
  "mixed-presence law", "N guards"), NOT the physics terms. The actual
  stale set was ~4 sites (heading + intro + one gotcha + one gate-norm
  phrase), not the brief's ~18. Reframe those; add ONE focused
  "**Where X lives — System B (the eviction)**" prose paragraph after the
  role-table that states the end-state (2-block carrier + own composite +
  `CoupledField` + retired presence-law + live-illegal-state
  unrepresentable + honest DOF + the `Solution.<member>` biconditional)
  and CROSS-LINKS the coupled 2×2 / M−N algebra to its home page — do not
  re-derive it. KEEP the total-phase-space eq-label (it is still the
  honest DOF sum; L-003) and just reframe the prose around it.
- **A RENAME (helper → new name) is a CLEAN `:func:` repoint, not a
  literalize-as-dead (L-013).** When the retired symbol SURVIVES under a
  new name (`_within_group_triple` → `build_within_group_system`, which
  also grew a record `WithinGroupSystem` carrying the named `A = M − N`
  splitting), repoint every `:func:`/literal to the LIVE name; keep ONE
  historical literal ONLY where the doc explicitly narrates the retirement
  ("the former ``_within_group_triple`` retired into this builder").
  Verify the rename is LANDED in live code first (`grep` the module: the
  old name gone, the new imported) — 13 dead `:func:` refs render
  plain-text with no `-W` warning (L-002), so the grep gate is the only
  catch. On the general/Cartesian sites, the record DEGRADES to the old
  triple (`.implicit_operator=(L+C)`, `.explicit_gains=(S,B_a)`) — say so;
  the coupled M−N grid is the CARRYING-mesh case only. (Those two fields
  were spelled `.resolvent`/`.gains` until 2026-07-28; `resolvent` was a
  misnomer — the field holds the un-inverted forward `M`, not `M⁻¹`.)
- **An IN-FLIGHT concurrent deliverable (the main agent editing prod while
  you write docs) gets post-state prose + a FLAGGED forward-dependency.**
  The LC-triplication collapse was being done in solver.py concurrently
  (`self.L` still live when I read it). The brief directed post-collapse
  prose; I wrote "the within-group L+C is spelled in ONE place — the
  builder" WITHOUT asserting solver internals (line numbers / `self.L`),
  and FLAGGED that the prose asserts the concurrent collapse for the main
  agent to reconcile at commit. Write what the brief directs, avoid the
  internals the brief forbids, flag the dependency.
- **A PLANNED-design admonition whose STRUCTURAL half landed via a
  SEPARATE commit: correct the false claims, keep the still-planned
  vision.** The `FaceField` codim-1 admonition said "no `FaceField` ABC
  exists yet" — but it LANDED (a separate C-series commit,
  `grep 'class FaceField'` + git log confirmed). Cardinal Rule 1: correct
  "PLANNED, not built" → "PARTIALLY built" (structural parent landed +
  ERR-067 metric refutation stands), note the SIBLING mechanism that
  reached the goal DIFFERENTLY (the eviction made presence
  unconstructable-by-design via System B, not via the still-unbuilt
  `PhaseSpaceCarrier`), and preserve the genuinely-unbuilt vision. Verify
  the guard COUNT is not restated (the old "7 guards" was already stale at
  17 live sites) — describe the guard FAMILY (`_require_/_refuse_/_require_leg_pair`)
  and the six-signature leaf-kwarg PROTOCOL (a `.. list-table::` of the
  apply/solve/transpose read-vs-fill contracts) instead of a call-site
  count, which sidesteps the L-014 count-drift trap.

How to apply: for an eviction/re-homing doc-pass, grep the CARRIER
framing not the physics; reframe narrowly + add one end-state paragraph
cross-linking the coupled algebra; clean-repoint the renamed builder
(verify landed, keep one historical literal); write post-state prose for
an in-flight sibling deliverable and flag it; correct a PLANNED
admonition's now-false claims while keeping its unbuilt vision.

---

## L-017 — The freed-name REMINT collision: a retired symbol's name reused for a DIFFERENT live object — disposition every mention by PASSAGE MEANING, not by name; and record a solve-leg un-weave in current-architecture passages only

The sequel to L-016 (the eviction): a refactor RETIRES a symbol AND
remints its freed name onto a **different, still-live** object. Here the
unified ψ½ leaf ``RadialCharacteristicField`` was retired (split into
``RadialCharacteristicInteriorField`` / ``…BoundaryField``) and its freed
name reminted onto System B's **composite** (the ``FullField`` mirror,
``Composite[interior, boundary]``). The name becomes a **homonym across
the remint commit** — the load-bearing discipline is to disposition each
mention by WHAT THE PASSAGE DESCRIBES, never by the name:

- **Current-architecture passage → REPOINT** to the live object (the
  composite ``~orpheus.transport.radial_characteristic_field.RadialCharacteristicField``).
- **Historical record → REWRITE as history** ("the unified single-buffer
  leaf, which then held the ``RadialCharacteristicField`` name; reminted
  onto the composite at 4e-e1b") — a literal, not an xref. A blanket
  find-replace is WRONG; the same name flips meaning at the remint.
- **Grep the FULL MODULE PATH, not the bare name.** A partial mechanical
  rename that already ran (e1b renamed the docs' ``RadialCharacteristicComposite``
  → ``RadialCharacteristicField``) leaves SOME refs correct (the composite
  path ``radial_characteristic_field.RadialCharacteristicField``) and
  OTHERS dead (the retired leaf's ``fields._bases.RadialCharacteristicField``,
  the unified space's ``…space.RadialCharacteristicSpace``, the unified
  ``…radial_characteristic_flux.RadialCharacteristicFlux``). The bare-name
  grep can't tell them apart; the module-path grep + import-gate can.
  Watch for a role FAMILY that split into ``interior ⊕ boundary`` (flux /
  source-sink / displacement / residual → 8 leaves): the unified per-role
  module paths (``…_source_sink.RadialCharacteristicSourceSink``,
  ``…_displacement.RadialCharacteristicDisplacement``) ALL die together —
  audit the whole family, not just the head symbol the brief names.

**Recording a solve-leg UN-WEAVE (inline orchestration → named resolvent).**
When a walk's welded inline orchestration is extracted to a NAMED operator
(``RadialCharacteristicOperator.solve`` / ``.solve_transpose``, the ``A_BB``
resolvent), the ENGINE ref (``carlson_inward_sweep_from_source``) stays
where it names the **march mechanism** (fine — the engine still exists and
is the single source of the march); the **orchestration** ref becomes the
operator's ``solve``. Record the un-weave in the CURRENT-architecture
passages ONLY — the Key-Facts bullet, the walk-triple SOLVE bullet, a
dedicated Cardinal-Rule-2 ``.. note::``, the six-signature protocol note,
and the Development-history changelog. The HISTORICAL saga (the Phase D–F
"one helper, two consumers" sections) that describes the OLD inline
architecture is PRESERVED per L-013 — it already carries the supersession
banner; do NOT rewrite it (proportionate). Two traps:

- **A re-aimed sentinel is a stale V&V claim (Directive 5).** The Mode-11
  wrap-sentinel was re-aimed from the engine onto the operator (class-level
  wrap of ``RadialCharacteristicOperator.solve``) + a NEW S2 "walk source
  has zero ``carlson`` references" tripwire. A doc saying "the sentinel
  confirms the solve executes ``carlson_inward_sweep_from_source``" is
  STALE — VERIFY against the live test (read the test body) and repoint to
  the operator + name the S2 tripwire. "carlson refs went N → 0" is a
  greppable fact — state the measured count.
- **A refuted-framing survivor on an adjacent HISTORICAL changelog line**
  (a "zero-metric" the ERR-067 pass corrected to SPD ``V_cell``, which the
  L-015 tree-wide sweep missed in the Development-history table): fix it
  when you edit the DEAD ref on the SAME line (Cardinal Rule 1 is supreme —
  a changelog is current-truth, not licensed to carry a refuted claim), and
  FLAG the scope-expansion. Don't retroactively inject the LATER
  architecture (System B) into an EARLIER-dated entry — keep the entry's
  date-accurate framing, only correct the dead ref + the refuted claim.

How to apply: for a freed-name remint, grep the full MODULE PATH (+ the
whole split role family), disposition each site by passage meaning
(repoint-live vs rewrite-history), record the un-weave in
current-architecture passages only (preserve the historical saga per
L-013), re-verify any re-aimed sentinel against the live test, and
fix-plus-flag refuted-framing survivors on lines you touch.

---

## L-018 — The CAPSTONE pass: documenting a COMPLETED multi-block coupled-operator architecture as one new taxonomy-culminating section

The completion (step-7) of the eviction/remint arc (L-016/L-017): the ψ½
ray, having moved to its own System B, is now one leaf of a full 2×2
**coupled block operator**, and the campaign LANDED. The doc job is a NEW
capstone `=` section documenting the whole architecture — not an
incremental reframe. Disciplines:

- **Place the capstone as the taxonomy's CULMINATION, cross-linking not
  duplicating.** The coupled block operator is the block-level
  generalization of the page's existing operator-surface taxonomy
  (three-layer surface → materialise → assemble → **N×N block grid**):
  apply→block matvec, assemble→block-offset scatter, solve→block
  substitution. Place it right after the last taxonomy sibling (the
  assembly axis), and OPEN by naming the generalization + `:ref:`-linking
  the sections it generalizes — the section EARNS its narrative place,
  it isn't a bolted-on appendix. Reconcile the pre-existing forward-ref
  ("the record bridges System B into a coupled M−N grid") by repointing
  it to the new same-page `:ref:` and softening stale framing ("bridges"
  evokes the retired bridge object).
- **A naming-dense brief on a landed architecture is an L-001 minefield —
  verify EVERY named object/line-ref/helper against the module-of-record
  before writing.** The brief conflated `A_BA` with
  `RadialCharacteristicReconstruction` (:955) — but LIVE `A_BA =
  RadialCharacteristicEmission` (:1187); Reconstruction is the Fold
  FACTOR *within* A_BA (`A_BA = Fold ∘ K ∘ integrate`). The brief's fold
  helper `fold_moments_to_starting_direction`
  (`starting_direction_space.py`) was renamed to
  `fold_moments_to_radial_characteristic` (`radial_characteristic_space.py`)
  by the campaign's OWN step-1 rename. Read the module docstring + class
  defs + import-verify each symbol before minting a cross-ref; a
  naming-dense brief's line-refs and class-names are the FIRST thing to
  go stale on a fast-moving branch.
- **Document a symbol OVERLOAD as an explicit gotcha, don't paper over
  it.** The class `RadialCharacteristicOperator` calls itself "A_BB" (the
  bare radial march μ∂_r+σ_t) AND the builder's local `A_BB = march −
  B_b` (the loss-grid self-block) — two live meanings. Faithful to the
  code-of-record: define A_BB = the bare march, spell the loss-grid (B,B)
  = A_BB − B_b explicitly ("a naming gotcha to spell out"), tie it to the
  System-A parallel (A_AA = L+C−S−B_a; both self-blocks = transport −
  gains − boundary; both boundary gains lagged in N). Reader inherits the
  overload cleanly instead of tripping on it.
- **The resolvent M and the loss grid A are DIFFERENT grids — state both
  precisely.** The brief's loose "M grid [[LC, Seeding], [None, A_BB]]"
  is the RESOLVENT M (bare-march diagonal, (B,A)=None → upper-triangular
  → direct block substitution). The LOSS grid A is [[L+C−S−B_a,
  +Seeding], [−Emission, A_BB−B_b]] ((B,A)=−Emission present). The
  splitting A = M − N puts B_b + the emission gain in N. Give both grids
  their own `.. math::`; note M(B,B)=bare-march while A(B,B)=march−B_b so
  A=M−N recovers. (B,A)=None-in-M IS the Schur/lag argument: the emission
  is the iterating scattering gain (ρ(M⁻¹N)=0.371), it belongs lagged on
  the rhs, not folded into the one-pass resolvent.
- **D5 for a campaign that RETIRES symbols: it owns its retirement
  doc-debt, but a BROAD PRE-EXISTING stale surface with NO 1:1 successor
  is FLAGGED, not rewritten-in-passing (L-007).** The campaign's OWN
  retired symbols (the ψ½ kwargs, the fused `CoupledInvertibleOperator`
  bridge, the presence guards, the `_within_group_triple`) were ALREADY
  correctly narrated as history by the incremental passes (the step-6
  "retired estate" section + dated changelog entries — literals in
  history blocks, not live xrefs) — verify each is history (grep the
  site's enclosing section title), leave it. But `transport_sweep`
  (retired by THIS campaign's step 6) had 55 refs across ~15 SN-page
  sections from many prior waves, MANY presenting it as current API, with
  NO mechanical successor (the sweep is now the resolvent's `.solve` /
  `sweep_schedule` per context — a behavioral rewrite per site).
  Rewriting ~40 current-API sites in a coupled-block docs pass is the
  exact L-007 anti-pattern. FLAG it (count + no-successor nature) as a
  dedicated sweep-entry-point-retirement pass; fix only campaign-adjacent
  current-API sites (here: zero — no transport_sweep site is in a
  coupled-block section).
- **New structural eq-labels for a landed architecture are all
  `documented` (L-004); one per load-bearing identity, grep-collision-
  check first.** Block matrix, block matvec, the fold, the M−N splitting,
  the loss grid, the free-identity residual, the block substitution —
  each representational/structural (not a solver claim), so `.. vv-status:
  <label> documented`. Grep `:label:` repo-wide before minting (all 9
  collision-free); verify in built HTML each `id="equation-<label>"`
  rendered and each in-prose `:eq:` resolved to `<a>`. Code-xrefs to
  the coupled_system/radial_characteristic classes render PLAIN-TEXT by
  page convention (not member-automodule'd — L-002; the pre-existing
  `build_within_group_system` ref already renders plain-text) — the
  import-gate is the real check, plain-text is NOT a regression.

How to apply: for a completed-campaign capstone, write ONE new
taxonomy-culminating `=` section (cross-link the siblings it generalizes,
don't duplicate); verify every named object/line-ref/helper against the
live module-of-record (naming-dense briefs go stale first); document
symbol overloads + the resolvent-vs-loss-grid distinction as explicit
gotchas; leave the campaign's own already-historical retired symbols,
FLAG the broad pre-existing no-successor surface; tag new structural
labels `documented`.

---

## L-019 — The context-dependent ENTRY-POINT retirement pass: per-site a/b/c disposition grounded in the LIVE successor, never a 1:1 rename

The EXECUTION of the L-018-flagged "dedicated retirement pass": a widely-
referenced entry point (here `transport_sweep`, 56 sites × 5 theory
pages) retires and its successor is **context-dependent** — the same
retired name maps to DIFFERENT live surfaces per site (production
resolvent `.solve` vs the scheduling layer vs a raise-guard vs an
SI-sweep twin method). A mechanical find-replace is FORBIDDEN (per-site
false-claim risk). Disposition EACH site into one of three, grounded in
what the LIVE code does THERE:

- **(a) behavioral rewrite** — the section teaches CURRENT API
  (present-tense claim, OR the whole section framing IS the retired
  entry). Two sub-shapes: an inline symbol-repoint (a present-tense
  claim swaps to the live surface — e.g. "the sweep at `X` consumes both"
  → "the within-group sweep (the resolvent `solve`) consumes both"); and
  a WHOLE-FRAMING rewrite (a section titled/built around the retired
  entry, with a stale code block — e.g. "Quadrature Dispatch" / "Typed
  input") gets its framing rewritten to the current architecture +
  cross-ref, and the stale code block **DELETED**, not symbol-swapped.
- **(b) past-tense history literal** — dev-history / changelog /
  retired-estate / diagnostic narrative. KEEP the name but as a
  double-backtick LITERAL (never a `:func:`), framed past-tense with the
  retirement citation ("the then-production ``X`` entry, since retired at
  step N (R-N.N)"). The campaign's OWN retirement passes usually already
  literalized their sites (grep: they're `` ``X`` `` not `:func:`X``) —
  LEAVE those; only literalize the dead `:func:`/`:meth:` refs OTHER
  waves left behind (they render plain-text, no `-W` warning, L-002 — so
  grep is the only catch).
- **(c) delete** — the clause carries no content once the symbol is gone
  (a "bit-identical to a direct ``X`` call" line; ``X`` as one entry in a
  list whose OTHER members survive → drop it, keep the survivors).

Disciplines that made it correct:
- **Ground EVERY successor in live code THIS session (L-001).** The
  live-grounded successor table is the load-bearing artifact — build it
  FIRST. A brief's "the successor is X" is a STARTING heuristic; the
  per-site live read is the rule (a raise-guard's kwarg was
  `moment_projection` in the doc but `moment_frame` in live code; a
  separate arg `Q_aniso` was GONE entirely — folded into the source; the
  SI-sweep twin pairs with the matvec `_apply_walk` as `…ScanWalk.sweep`,
  found by reading the walk class's methods).
- **For a big superseded-architecture SECTION prefer the L-013 arc-head
  supersession banner + past-tense over a full rewrite** (proportionate).
  A dev-history section framed "Wave-X did Y" that presents the retired
  entry with a code block: add ONE `.. note:: **Superseded (step N)**`
  naming the current architecture + `:doc:` link, past-tense the intro
  verbs, let the historical code block stand under the banner. Reserve
  the full framing-rewrite for CURRENT-architecture reference sections
  (early-doc, NOT "Wave-X"-framed).
- **Retitle a heading that NAMES the retired entry** only after grepping
  for inbound `:ref:` (autosectionlabel); size the new underline in code
  points (L-009).
- **Acceptance = the three-severity `-E -W` count UNCHANGED-from-baseline
  + a `git grep <symbol> -- docs` survivor audit** where every survivor
  is a past-tense literal (no live xref; none presenting the symbol as
  current). The skip-line DIFF (baseline vs post) proves you orphaned no
  OTHER verifies-target while resolving the in-scope ones.

How to apply: build the live-grounded successor table FIRST; disposition
each site a/b/c by what live code does there; banner-not-rewrite the big
history sections; literalize every dead `:func:` (grep-gate — `-W` is
blind); prove the skip-line diff changed only what you intended.

---

## L-020 — The retired symbol whose deletion is a COROLLARY of a design unification: the enclosing section's THESIS is stale, not just the symbol

The sharpening of L-019 for the hardest retirement case. L-019 dispositions
each SITE a/b/c and banners a "Wave-X did Y"-FRAMED history section. L-020
is the case L-019's per-site list UNDERSELLS: the brief hands you a few
"dead-role lines" but those lines are the visible symptoms of an entire
doc-SECTION whose THESIS (its design-rationale premise, presented as
CURRENT) has gone false — because the symbol you were sent to retire is one
FACET of a deeper architectural UNIFICATION that dissolved the very design
the section exists to explain.

Worked (Task #57, `transport_operator_matvec` family + `psi_bc`): the brief
listed 6 sites (discrete_ordinates 13543-13545/13576/13700-13704 +
index_convention 488/1507/1524). Reading the LIVE code (streaming.py apply →
`loss_action` → `_apply_walk`; `loss_representation/__init__.py:1458`
"matvec ≡ sweep, ONE discretization, L21") revealed the matvec deletion was
a corollary of the #206 Phase-C **matvec ≡ sweep = ONE loss-representation
walk** unification — which DISSOLVED the "two distinct discretisations
(FD-operator apply vs WDD sweep) / packed-vector-vs-structured-array layout
difference / deliberately-legacy-pending-PR-INDEX-7" design that FIVE
sections were built to teach as CURRENT. The stale unit was the section
PREMISE, not the line.

The tell (distinct from L-019's "Wave-X-framed" history): a section stated
its stale design in the PRESENT tense as a live rationale — "apply and solve
use different closures **by design**", "What **stayed** deliberately legacy",
"`apply` operates on the **packed 1-D solution vector**". A per-site literalize
leaves the false THESIS standing.

Disposition ladder for a thesis-stale section:
- **THESIS-stale reference/rationale section → SUPERSESSION BANNER at the
  section head** stating the unification + the current one-truth + the campaign
  cite, then past-tense the body verbs + literalize dead roles, PRESERVING the
  historical reasoning under the banner (the Wave-D two-closure narrative +
  ERR-026 stayed — it is still pedagogically load-bearing). Retitle
  "…**(historical)**". This is L-019's banner move promoted from
  history-framed sections to present-tense-rationale sections.
- **The ONE genuinely stale-AS-CURRENT contract → full behavioral section
  REWRITE** to the unified live contract (here the psi_bc/Q_aniso "Vector
  layouts" bullet list → the `FullField` composite: source on
  `rhs.interior.values`, boundary on typed `AngularBoundaryFlux` face views),
  with the retired triple recorded in a trailing `.. note::`. This rewrite is
  what actually KILLS the "Persistent boundary-flux dict psi_bc carrying state"
  bullet (grep it → 0).
- **A moot FUTURE-WORK section (the unification made the planned migration
  unnecessary) → retitle "…(obsoleted)" + `.. note:: Obsoleted by deletion`,
  PRESERVING its `:ref:`-target label** (an `.. _future-work:` anchor is often
  referenced from elsewhere; deleting it dangles those refs — keep the label,
  make the CONTENT truthful). Grep the label before touching it.
- **Co-literalize deleted-SIBLING roles only INSIDE clauses you're already
  rewriting** (a `:func:solution_to_angular_flux` sitting in the same sentence
  as your matvec fix — leaving a live dead-role in a line you're editing is a
  self-inflicted Rule-1 staleness). FLAG the standalone sibling-cluster
  (EquationMap/codec) for its OWN pass — don't chase it tree-wide (L-007).

Disciplines carried from L-019 unchanged: ground every successor in live code
FIRST (L-001); grep-gate every new `:func:`/`:class:`/`:meth:` xref against
live code AND (for Protocol methods) a python `hasattr` probe (L-002 — dead
roles render plain-text, `-W` blind); size retitled underlines in code points
(L-009); acceptance = the three-severity `-E -W` count UNCHANGED-from-baseline
+ a `git grep <symbol> -- docs orpheus` survivor audit where every survivor is
a past-tense literal.

How to apply: after the successor table, read each flagged line's ENCLOSING
SECTION and ask "is the section's PREMISE still true under live architecture?"
If the symbol's deletion is a corollary of a unification (grep the live code
for the "one X" / "matvec ≡ sweep" / "dissolved" fact), the premise is stale
too — banner the rationale section, rewrite the one stale-as-current contract,
obsolete-but-preserve the moot future-work section, and preserve the historical
reasoning under every banner.

---

## L-021 — The bulk-scanner staleness sweep: I am the JUDGMENT layer over a precision-over-recall Haiku pass — re-verify EVERY finding vs LIVE, reject stale-evidence, and a scanner suggestion is a STARTING POINT not the truth

A 200+-finding automated staleness sweep (Haiku scanners, "precision over recall
with command evidence") is compiled by a weaker model; the dispatch brief itself
says "you are the judgment layer". The recurring failures are (a) trusting a
scanner's suggested TARGET, (b) trusting a scanner's stale WARNING-evidence, and
(c) transcribing a fix that mints a NEW false claim. Disciplines that held across
29 files / 120→15 dead-role reduction:

- **A scanner's suggested target can name a symbol that does NOT exist.** The
  scanner conflated a retired Protocol name with the live class:
  `AngularQuadrature.spherical_harmonics` → the scanner said "retarget to
  `orpheus.numerics.quadrature.AngularQuadrature.spherical_harmonics`" but
  `AngularQuadrature` does not exist anywhere (`hasattr` False) — the live class
  is `Quadrature`. ALWAYS `hasattr`/import-verify the SUGGESTED target, not just
  confirm the OLD one is dead. The census resolver (getattr-chain longest
  importable prefix) is the cheap gate: batch-probe every candidate target before
  editing.
- **A scanner's "sphinx-build warns: X not found" evidence can be STALE — check
  the CURRENT clean build.** One finding claimed a `:mod:` bare ref "warns py:mod
  target not found"; the current `-W` build had NO such warning (bare `:mod:legacy_name`
  renders plain-text WITHOUT warning, L-002). REJECT the finding (evidence
  doesn't hold), and if the ref is a page-wide legacy-naming convention used N×,
  leave it (fixing all N is a rename beyond flagged scope) — record the rejection.
- **Reproducing a fix routinely surfaces a NEW false claim the scanner didn't
  flag — fix the CLAIM, not just the role (Cardinal Rule 1 > scope).** Two
  instances THIS pass: (i) a doc said "verified against **SymPy**'s
  `sympy.integrate(...)`" but the live test uses `scipy.integrate.quad` singular
  quadrature — corrected the pillar attribution, not just the test-name role.
  (ii) I nearly wrote "`wigner_seitz_pin_cell` default of 10 fuel + 3 clad + 7
  coolant sub-cells" transcribing a `pwr_pin_equivalent` default — but
  `wigner_seitz_pin_cell` produces region THICKNESSES only; the sub-cell COUNT is
  a `RegionMesh(n_cells=...)` choice at `Mesh1D.from_geometry`. Caught my own
  mid-edit false claim by reading the successor's live body (L-001 applied to my
  OWN draft). RULE: when a retarget crosses a numerical/structural CLAIM, read
  the successor's live def AND re-verify the surrounding claim's truth before
  re-spelling it.
- **THESIS-STALE beats symbol-stale when a whole class/design was deleted.** A
  deleted carrier class (`GeometrySpec`) staled not just its `:class:` refs but
  the SCHEMA table, the migration narrative (prospective "Phase B will split"
  when the split is DONE + the carrier since deleted), the bridge-test bullets,
  AND the rationale section — a full-section rewrite to the current
  `geometry_kind: str` + `to_geometry()` form. Similarly a scattering per-ℓ
  ladder retired for a Funk-Hecke `R∘Λ∘M` kernel staled a `:class:OperatorSum of
  per-ℓ leaves` TABLE ROW (present-tense architecture), not just the Code-Anchors
  entry. Read the ENCLOSING claim's premise against live architecture (L-020),
  not just the flagged line.
- **A `:ref:` anchor placed AFTER a section title does not resolve as a title
  target — put `.. _label:` ABOVE the title.** The only real WARNING I introduced:
  a forward `:ref:` to an anchor I'd defined between a title and its body. Move
  the anchor above the underlined title (blank line between) so `:ref:` picks up
  the section title.
- **Mechanical families are the bulk win — replace_all on unique multi-char
  strings, count-asserted.** CP/MoC test-filename drift (`test_cp_*.py`/`test_moc_*.py`
  → `tests/{cp,moc}/test_*.py`), a module rename (`test_galerkin_spectral_symbolic`
  → `test_carlvik_galerkin_symbolic`, ×8, fn names unchanged), short-path role
  normalization (`~geometry/numerics/data/sn.*`→`~orpheus.*`), numbered-dir relics
  (`07.Thermal.Hydraulics/...`→`orpheus/thermal_hydraulics/solver.py`). Use a
  count-asserting Python script (assert count==1 for uniques, log count>1 for
  known-duplicates) on the LARGE files where reading 20k lines is impractical —
  it's auditable and can't introduce a transcription error. Order matters:
  replace the LONGER stem first (`thermal_hydraulics_dae.py` before
  `thermal_hydraulics.py`) even when disjoint, to be safe.
- **The residual census MUST be fully attributed.** After the sweep, every
  remaining dead-role census hit must map to a KNOWN false-positive class
  (dataclass-field/ctor-param without default → resolver getattr-chain
  limitation; historical framing "no longer exists"; planned/future `(planned)`/
  "scheduled to promote"). Enumerate them; re-read the doc context of any you're
  unsure of; list any that DON'T fit and fix them. 15/15 residuals attributed
  this pass — none required fixing.

How to apply: for any bulk-scanner-compiled fix job, batch-probe every SUGGESTED
target (not just the old one) against live code FIRST; reject findings whose
build-warning evidence the current clean build contradicts; when a retarget
crosses a numerical/structural claim, read the successor's live body and
re-verify the claim; rewrite the enclosing THESIS when a deleted class/design
staled the whole section; and attribute every census residual to a known FP class.

---

## L-022 — The pedagogical RESTRUCTURE + cross-page THEORY-EVICTION pass: keep-anchor makes it ref-safe; the real traps are marker-depth shift, the part-boundary blank line, and the general-vs-consumer split

A page RENAME + section-tree reorg that also EVICTS large theory blocks
from a sibling page into it (here `galerkin_projection.rst` → `frame.rst`,
retitled + reorganized into a Petrov-Galerkin-first tree, absorbing the SN
page's ~2100-line homogenization + condensation GENERAL theory). Distinct
from L-012 (merging a re-staged BRANCH diff into a diverged tree): this is
a clean intra-repo relocation + pedagogical reflow. Disciplines:

- **A cross-page THEORY MOVE is fundamentally ref-safe when you KEEP the
  labels (don't rename them) — Sphinx labels, eq-labels, AND citations
  all resolve GLOBALLY cross-doc.** Moving a `.. _anchor:` + its content
  to another file keeps every `:ref:`/`:eq:` working with ZERO referrer
  edits, as long as each label stays defined exactly ONCE (move, don't
  copy). This is L-007's keep-the-anchor applied to a relocation: the
  brief's "move the label with the content and fix referrers" needs NO
  referrer fixes at all. Proof-of-pattern for citations: `[Hebert2009]_`
  was already defined in `collision_probability.rst` and used in
  `discrete_ordinates.rst` — a live cross-doc citation — so migrating
  `[WIMSD]_`/`[Rahnema2008]_` REFERENCES into a third page resolves to
  their definitions on the SN page identically (verify with the build).
  The ONLY rename that breaks refs is the DOC name (`:doc:X`): fix the
  toctree entry + every `:doc:old` (3 sites here); the page's own
  `:ref:`-label (`galerkin-projection`) is KEPT so its 6 external
  referrers need no touch.
- **The PART-BOUNDARY blank-line trap (the one real build break).** When
  you assemble a page programmatically by concatenating slices, a slice
  that ENDS in content (a migrated block sliced to its last content line,
  no trailing blank) joined directly before the next part's `.. _anchor:`
  GLUES the anchor to the preceding paragraph → the label silently fails
  to register. The symptom is **"undefined label"** (NOT "duplicate") at
  every referrer, even though `grep` shows the anchor present and at
  column 0. Fix: join parts with `\n\n` (guarantee ≥1 blank everywhere)
  or ensure every content-ending slice carries a trailing blank; then a
  triple-blank normalizer caps runs at 2. A col-0 grep for
  "`.. _x:` whose previous line is non-blank" is the pre-build catcher.
- **Marker-depth shift when content re-nests under a DEEPER parent.**
  Evicted `=`-level sections whose subsections were `-`/`~` land under a
  new `-` subsection (`§2c Applied to …`), so every migrated underline
  demotes one level (`-`→`~`, `~`→`^`), LENGTH-PRESERVING (char-for-char
  swap keeps the code-point count matching the unchanged title, L-009).
  Detect a section underline robustly as a **col-0** all-one-marker line
  (len≥4) whose previous **col-0** line is a plain title (non-blank, not
  `..`, not itself all-marker) — col-0 disambiguates it from bullet `-`,
  math, and `* -` list-table rows. Strip the evicted block's top `=`
  title (hand-write the new `-` title + KEEP its anchor); shift only the
  body.
- **The general-theory vs consumer-orchestration split is the judgment
  call — prefer MOVING general theory, keep the consuming stub lean.**
  Homog/cond theory (rate preservation, the PG-frame derivation, the
  metric-fold-vs-bilinear adjoint argument, fractional-overlap, the
  asymmetry law, AND the verification gates — they verify the GENERAL
  property even when the tests live under `tests.sn.*`) → the theory
  page. Only the SN-LAYER orchestration stays in the stub: which driver
  invokes it (`Solution.homogenize`→`MaterialMesh`→re-promote loop;
  `Solution.condense`→per-material representative spectrum→`dict[int,
  Mixture]`), plus the ONE SN-specific equation
  (`energy-condensation-representative-spectrum` — moved to the stub, not
  the theory page) + `:doc:` links to the full treatment. Split the
  evicted section's INTRO too: its general "what X is" framing →
  theory page; its "in ORPHEUS it lives as `Driver.verb` returning T"
  sentences → stub.
- **PROMOTING a buried `.. note::` to a proper subsection (the crux
  content):** preserve the argument VERBATIM (inline math copied
  exactly), convert `.. note::` prose to subsection prose, KEEP its
  `.. _anchor:` (referenced from elsewhere — here from the "unifying
  principle" §), add ONLY the ONE current-architecture design-rationale
  sentence the brief asks for ("the frame was first posed as Galerkin;
  the adjoint requirement forced the re-posing as Petrov-Galerkin" — a
  REASON, not dated process-narrative, L-010), and fix its now-intra-page
  `(:doc:otherpage)` parentheticals.
- **Clean up `(:doc:sibling)` parentheticals that became intra-page.**
  Every `:ref:X (:doc:discrete_ordinates)` where X migrated INTO this page
  is now a wrong forward-pointer (the reader is sent to the wrong page).
  `-W` does NOT catch it (the `:ref:` still resolves globally; only the
  parenthetical lies). Grep `discrete_ordinates` in the destination page
  → 0 after the pass is the gate. Distinguish from a `:ref:` to SN
  content that STAYS (`sn-scattering-adjoint`) — that keeps its cross-doc
  form.
- **Mechanics that held:** programmatic slice-and-reassemble with
  boundary ASSERTS (`F[idx]==expected_title`) fail-loud on line drift;
  compute the new strings and run ALL structural asserts on the in-memory
  result BEFORE writing either file (a failed assert then leaves the tree
  untouched — no `git checkout` recovery needed, process-discipline);
  assert label counts (each verifies-target defined once, in the
  destination; the SN-only eq NOT leaked), the new headers present, the
  retitled originals GONE, `discrete_ordinates`∉destination. The clean
  `-W` build + the regenerated `matrix.rst` mapping every moved
  verifies-target to its tests (17/8/3) is the final proof no edge
  orphaned.

How to apply: for a rename+reorg+eviction, KEEP labels (move don't copy →
ref-safe); fix only `:doc:` (toctree + `:doc:old` sites); `\n\n`-join or
trailing-blank every slice (the glued-anchor→undefined-label trap);
depth-shift migrated underlines col-0-detected + length-preserving;
move general theory + verification, keep a lean SN-orchestration stub with
the one SN-specific eq; promote a buried crux note verbatim + one
rationale sentence; scrub now-intra-page `(:doc:sibling)`; gate on
in-memory asserts before write + the clean `-W` build + the matrix.

---

## L-023 — The template-skeleton ADDITIVE front-matter pass: machine-header dropdown, Key-Facts+Overview→Synopsis fold, gotchas consolidation, essay eviction — under a HARD "don't reorder the middle" non-goal

Distinct from L-022 (a rename+reorg+cross-page eviction): here the SAME
flagship page gets the 9-section template's FRONT MATTER imposed
ADDITIVELY, with an explicit SCOPE non-goal forbidding physical reorder of
the large middle blocks (the aggressive re-level deferred to a later
phase). The job is: add §1 machine header, fold the front matter into §2
Synopsis, consolidate scattered gotchas at the tail (§8), stub §5/§6 with
automation-pending notes, and evict the narrative history essays (§9 keeps
only the changelog). Disciplines:

- **Machine header = a collapsed sphinx-design `.. dropdown::` wrapping a
  `.. code-block:: yaml`, NOT a bespoke directive.** VERIFY the intended
  ingestion directive is unregistered first (`grep add_directive` in the
  INSTALLED pkg — here `sphinxcontrib.nexus` registers `nexus-graph` /
  `verifies` / `implements`, NOT `nexus-meta`); an unknown directive fails
  `-W`. The YAML-in-code-block sidesteps EVERY RST-parse hazard (unicode
  μ/Σ/ᵀ/≡/→/∈, `#` comments, quotes, no roles/labels) because a code-block
  is literal text. Author only what the graph can't cheaply derive
  (conventions, invariants, operator glossary, retrieval aliases); keep
  entry_points/key_types/governing_equation MINIMAL (Nexus derives the
  fuller lists). Populate conventions/invariants FROM the current Key-Facts
  bullets, each fact verified against live code (L-001) — the header is
  CURRENT structured fact, so a stale convention there is a Rule-1 bug.
- **Key Facts + Overview → ONE Synopsis (front-matter FOLD, not two
  sections).** The structured convention/invariant bullets become the
  machine-header YAML (data for the machine); the prose framing becomes a
  dense, NAMED, retrieval-targeted synopsis (the primary embedding target —
  cite methods/operators/decisions BY NAME); the load-bearing clickable
  `:ref:` navigation list stays as a `.. admonition:: Conventions` block
  (the data moved to YAML, but the reader still needs the nav — keeping
  both is NOT a twin-source-of-truth violation: YAML serves the machine,
  the admonition serves the human). Fold Overview IN (don't keep a near-
  duplicate section beside the synopsis — that WOULD be a twin path).
  PRESERVE every citation usage (`[Bailey…]_` etc.) so no definition
  orphans. FLAG the two-sections→one fold for main-agent review (it is the
  one consolidation judgment call in an otherwise-additive pass).
- **`autosectionlabel` OFF ⇒ a section-heading RENAME is fully ref-safe**
  (no implicit anchor exists to break; only explicit `.. _label:` anchors
  are targets). Confirm it's off (`grep autosectionlabel docs/conf.py`) +
  grep the tree for inbound `:ref:` to the old heading's slug (here: none
  to Key-Facts/Overview) before renaming Key Facts → Synopsis.
- **Essay eviction: an intra-page `:ref:` to a deleted essay anchor IS
  `-W`-CAUGHT (unlike a cross-doc `:ref:`, which renders plain-text —
  L-002).** Grep referrers to BOTH essay anchors BEFORE deleting: essay-1's
  anchor had zero referrers (safe outright delete); essay-2's anchor had a
  single INTRA-page referrer, so deleting it would dangle → `-W` fail.
  Repoint the referrer to the DISTILLED-gotcha anchor — the distillation
  (the "why the two errors cancel for homogeneous problems" →
  homogeneous/uniform-rescale gotcha) IS the semantic successor of the
  evicted essay, so the repoint is meaning-preserving, not a band-aid.
  Paste each evicted essay's FULL text into the return for issue
  relocation (git also preserves it) — the outcome already lives in the
  changelog + ERR catalog, but the narrative goes to the originating issue.
- **HARD SCOPE non-goal vs a deliverable's own example → the non-goal
  governs; point, don't move; FLAG.** §D said "move the scattered `~`
  gotcha subsections", naming one at ~L11835 — but that one sits INSIDE the
  SCOPE section's protected range ("do NOT touch the Sweep mega-section
  internals L2534–11895, that is Phase 1d"). When a hard non-goal collides
  with a deliverable's example, the non-goal wins: LEAVE the protected-range
  item, add a `.. seealso::` pointer to its anchor from the consolidated
  section (consolidation-by-discovery), FLAG for the main agent. Extra
  confirmation here: the sweep gotcha's `.. _sn-282-gotchas:` anchor was
  referenced ONLY from within the protected range — moving it would drag
  refs across the Phase-1d boundary. The OTHER named gotcha (~L15444, under
  SNSolver, OUTSIDE the range, no anchor, no referrers) moved cleanly.
- **The three degeneracy gotchas were `**Gotcha**:` BULLETS in Key Facts,
  not a section** — they are load-bearing facts that belong in §8, so lift
  them OUT of the front matter INTO the Gotchas section (as
  consequence→manifests→catcher warning boxes), don't leave them to
  duplicate in the synopsis. Verify the named catcher tests EXIST
  (`grep def test_…` — both `test_dd_per_cell_recurrence_matches_symbolic_derivation`
  and `test_heterogeneous_absolute_keff` did; the L0 one had been RENAMED
  by an earlier fold, so cite the successor the evicted essay already
  named, L-001).
- **§5/§6 stub notes are ADDITIVE `.. note::` blocks at the section HEAD,
  no restructure.** Each names the blocking Nexus issue (flow-graph
  nexus#20 for the Implementation-map/Architecture section; label↔test
  linking for the Verification slice) + "hand-authored until then" + a live
  `:doc:` to the real surface (`../api/numerics`, `../verification/matrix`
  — verify the target exists).
- **Mechanics (L-022 reused):** one programmatic slice-and-reassemble
  script — content-anchored `.index()`/count-asserted `.replace()`,
  code-point underlines via `len(title)`, structural asserts (new anchors
  ==1, old anchors ==0, dangling-ref ==0, machine-header<synopsis<Architecture,
  moved-bullet ==1) that FAIL LOUD before any write (no partial write, no
  git-checkout recovery). Extract the moved gotcha bullets by substring
  (between markers) so no retype. Acceptance = the `-E -W` build EXIT=0 with
  0 WARNING/ERROR/CRITICAL (this branch's Phase-0/1a baseline is a CLEAN
  build, not 1) + HTML render audit (dropdown `sd-*` classes present, new
  `id=` anchors present, repointed `:ref:` renders `<a href="#…gotcha">`).

How to apply: for a template-skeleton additive pass, (1) machine header =
collapsed dropdown + YAML code-block (verify the directive is
unregistered); (2) fold Key-Facts+Overview → one named Synopsis (structured
→ YAML, prose → synopsis, nav-refs → Conventions admonition, flag the
fold); (3) evict essays after grepping BOTH anchors' referrers (intra-page
`:ref:` IS `-W`-caught → repoint to the distilled gotcha), paste full text
into the return; (4) when a hard non-goal collides with a move instruction,
leave-and-point-and-flag; (5) programmatic slice + asserts + clean `-W`
build + HTML audit.

---

## L-024 — In a structural chapter-split, the single-homed check is on the anchor DEFINITION, not raw label mentions: a router page KEEPS its bare `:ref:` back-refs to labels it exported

The `sn_split_catalog.md` STEP-5 wording ("`grep -c '<label>'` in
index.rst MUST be 0") is a *proxy* that only holds when EVERY inbound
ref to the moved label lives in a DIFFERENT file — true for ch3
(`quadrature-types`, sole inbound from MoC), FALSE the moment the
exported label is back-referenced from the router page itself. The
load-bearing check is the anchor **DEFINITION**: `grep -c '^\.\.
_<label>:' index.rst` == 0 (source no longer owns it) + `grep -rn
'^\.\. _<label>:' docs/` == 1 (new chapter owns it). BC was the first
Phase-C cut to hit this: `boundary-conditions` is back-referenced
from index.rst:171 and :16169 (was :16577, shifted −408), so raw
`grep -c` = 2 while anchor-def = 0. Both survivors are bare
`:ref:`boundary-conditions`` — after the cut they become path-immune
CROSS-doc refs (resolve globally, NO `-W` warning; L35 family). Do
NOT "recut" on a nonzero raw count — inspect WHICH lines matched; a
bare `:ref:` with no directional word and no `:doc:` page-qualifier
STAYS. The genuine recut triggers are exactly two: a surviving `..
_<label>:` DEFINITION in source, or a phantom `duplicate label` under
`-E` (L36).

Second BC-specific trap: a bystander page-qualifier on the general
foundations page — `at :ref:`bc-tensor-decompositions` (in
:doc:`/theory/methods/sn/index`)` — is a TRUE falsehood (L35 case c),
even though `-W` stays silent. The `:ref:` is path-immune, but the
adjacent `:doc:` NAMES the old home; `:doc:` targets ARE
path-sensitive, so after the cut the prose sends the reader to the
page the label just LEFT. Repoint the `:doc:` to the new chapter
(`/theory/methods/sn/boundary_conditions`). NB this creates a brief
window where the `:doc:` DANGLES (target file not yet created) — so
order the moves: fix the qualifier (STEP 1) → create the chapter file
(STEP 3) → build (STEP 6); never build in the gap.

Third trap (Verification chapter, 2625-ln / 41-label cut): a name can
be a `.. _X:` **section anchor** (std:label domain, `:ref:`-resolved)
AND a `.. math:: :label: X` **equation label** (math domain,
`:eq:`-resolved) SIMULTANEOUSLY — DIFFERENT Sphinx domains, so they
coexist with NO duplicate-label warning (a clean `-E` baseline proves
it). When the two live in different parts of the page and the split
moves ONLY one, single-homing MUST check the CORRECT namespace: a
MOVING eq-label is verified with `grep -c ':label: <name>'` (== 0 in
source, == 1 tree-wide), NOT `grep -c '^\.\. _<name>:'` — the latter
shows the same-named SECTION anchor that legitimately STAYS and would
falsely read as a recut trigger. Instance:
`sn-mms-{spherical,cylindrical}-aniso-spatial-convergence` each exist
as a section anchor in the sweep area (stays in index.rst) AND an
eq-label in the Verification block (moves to verification.rst); after
the split `:ref:`→section (index) and `:eq:`→equation (verification)
both resolve globally cross-doc. So a name owning BOTH namespaces is
TWO independent single-home checks, not one — and the `-E` build is
the collision oracle (splitting the pair across two files stays clean
BECAUSE the domains differ). Also note the END-boundary rule for a
next-section whose `.. _anchor:` sits ABOVE its header: stop the cut
BEFORE that anchor (it belongs to the STAYING section), not at
header−1 — else the split drags the next section's anchor into the
new file.

How to apply: for any verbatim chapter cut, run the L35 three-way
grep and disambiguate the moving section-anchor label from any
SUPERSTRING label it collides with (`boundary-conditions` vs the
foundations `theory-boundary-conditions`). Verify single-homing with
`grep -c '^\.\. _<label>:'` (definition), not `grep -c '<label>'`
(mentions). Report the raw count too, naming each surviving line as a
legitimate bare back-ref. Fix ONLY (a) intra-source directional prose
whose target left, (b) moved-block prose whose stay-behind target it
now mis-qualifies, (c) bystander `:doc:` page-qualifiers naming the
old home — leaving bare `:ref:`/`:eq:` (path-immune) untouched.

---

## L-025 — AUTHORING a NEW keystone foundational chapter (gather method-specific verified math → GENERALIZE to the abstract object): the within-doc symbol-collision hunt is the sharpest new gate, and the algebra-of-record is the correctness spine

Distinct task-shape from the split/eviction/close-out passes: writing a
NEW shared method-invariant chapter (`foundations/discretization.rst` —
cell balance + Step/DD/LD, derived once so SN/CP/MoC never re-derive).
The move is *gather the already-verified method-specific derivations,
GENERALIZE them to the abstract object, author generically*. Disciplines:

- **The algebra-of-record is the correctness spine — READ it, then RUN
  it, before writing a single equation.** The SN-specific SymPy
  (`orpheus.derivations.discrete.sn.balance`, 7 `derive_*` foundation
  tests) IS the source of truth for the cell balance / DD / WDD /
  curvilinear / flat-flux math (algebra-of-record skill). Generalizing
  = stripping the SN-specific specialization (face areas → 1 for the
  planar cell, keep the generic `|μ|(ψ_out−ψ_in)+Σ_t h ψ̄ = q h`) while
  quoting the SymPy identity verbatim. `pytest`-run the module (7/7
  green) so the presented equations are grounded, not transcribed
  (Cardinal Rule 1). A DIFFERENT concept in the same page can have a
  DIFFERENT algebra-of-record module (LD's 2×2 is `ld_ubld.py` +
  `_ubld.py`, NOT `balance.py`) — name each precisely in the seealso;
  don't let one `:mod:` cite stand for all the math.
- **A NEW page assembled from MULTIPLE overloading sources is the prime
  site for a WITHIN-document symbol collision — hunt it before the
  build (the build is BLIND to it).** L-011 was a cross-PAGE convention
  sweep; this is the NEW-page twin: when you gather from sources that
  each overload a letter, the SAME letter can carry two meanings in ONE
  page. Here `w` = the cell-average BLEND weight (Step w=1, DD w=½) from
  the scheme code AND `w` = the angular QUADRATURE weight in the SymPy
  `ΔA/w` geometry factor — a genuine ambiguity a teaching page must not
  ship (Cardinal Rule 1). `-W` NEVER catches it (both render fine). The
  gate is a manual re-read hunting every reused glyph across the
  gathered sources; resolve by subscripting the lower-frequency meaning
  (`w_m` for the quadrature weight) + a one-line disambiguation note at
  first collision ("`w_m` is the quadrature weight — *not* the blend
  weight `w`; they share a letter in the literature"). Grep the fix is
  complete (`grep -nE '\\Delta A\}\{w\}'` → 0 residual).
- **A designed-but-UNBUILT scheme is a code LITERAL + a traits-anticipated
  note, never a `:class:` (L-002 forward-ref rule, applied to a THIRD
  sibling).** Step is the open #158 arm — `class Step` exists only as a
  docstring EXAMPLE in `scheme.py`, not an importable class. Grep
  confirms (`grep -rn 'class Step' orpheus/` → only the docstring line).
  So write ``Step`` (literal), derive its math on its own merits (the
  w=1 case of the verified blend framework — code-grounded even without
  a Step-specific SymPy), and add a `.. note:: ORPHEUS status` naming
  the anticipated traits on the base Protocol (`is_positivity_preserving
  = True`, the w=1 upwind blend) + the issue. Honest: documents the
  concept fully without claiming a link to a symbol that isn't there.
- **A LOAD-BEARING worked example gets RUN against the live code before
  it's written (Cardinal Rule 1 / Quality item 6).** The thick-cell
  negative-flux contrast (Step +1/5 positive, DD −1/3 strongly negative,
  LD −1/19 mildly negative — "LD better positivity than DD" made
  concrete) was computed by hand THEN reproduced through the real
  `DiamondDifference.update` / `LinearDiscontinuous.update` (build a
  `StreamingTerms` + `CellVisit` + `UpstreamState` directly; `slab_streaming`
  is a mesh-level factory, not a single-cell builder — construct the
  frozen dataclass fields by hand). Exact fractions matched to machine
  precision → the example is verified, not plausible.
- **The two-axis articulation catch: a "spectrum" parameter is often ONE
  of two orthogonal axes — say so, or a fresh reader over-collapses.**
  The blend weight `w` (Step 1 → DD ½ → LD adaptive 1/(1+k)) orders the
  FACE reconstruction; the MOMENT count (Step/DD: 1, LD: 2) is a SECOND,
  orthogonal axis, and the diffusion limit lives in the moment count NOT
  the blend. "LD→Step as w→1 (thick)" is TRUE of the face blend and
  FALSE of the scheme (LD keeps the slope Step lacks). An `.. important::`
  block pinning "the blend is one axis, the moment count another"
  pre-empts the natural mis-read — the doc-side inoculation move
  (vv-principles Mode-10 doc companion generalized to any seductive
  over-collapse).
- **Zero-warning gate on a NEW standalone page ⇒ ALL references PLAIN-TEXT
  (no `.. [Key]` defs, no `[Key]_` cites).** L-006 says standalone pages
  accept duplicate-citation warnings as a trade-off — but a strict
  "0 warnings" brief forbids introducing ANY (a `[Key]_` cite is
  undefined on the new page → warns; a `.. [Key]` def collides with the
  existing def elsewhere → dup-citation warns). Resolve by citing every
  reference as PLAIN TEXT in a Literature `.. list-table::` (author,
  year, title, journal, the load-bearing eq/§ numbers inline) — zero
  citation machinery, zero warnings, AND higher articulation (the
  equation numbers sit in the prose). Mixing plain-text + `[Key]_` is
  inconsistent; go all-plain-text.

How to apply: for a NEW foundational chapter — (1) read+RUN the
algebra-of-record module(s), one per distinct concept; (2) get the
disassembly/outline right before prose (the 7-part skeleton), namespace
every section anchor AND every eq-label to the page (`<page>-` prefix,
grep-confirmed collision-free); (3) hunt within-doc symbol collisions
across the gathered sources (subscript the rarer meaning + a note);
(4) code-literal any unbuilt sibling; (5) RUN every load-bearing worked
number through the live code; (6) all-plain-text references; (7) tag
every eq-label `.. vv-status: <label> documented` with a rationale
comment naming the SymPy/foundation gate (L-004 — they land in
Documented-only, not flagged as unverified solver claims).

---

## L-026 — The split-to-new-page pattern (extract N contiguous H1 sections into a foundations page): identify by STABLE title, slice PROGRAMMATICALLY, carry every label; and the build-INVISIBLE f-string LaTeX-brace trap in the AUTHORED header

The #231 corpus campaign repeatedly splits a monolith page's advanced
deep-dives into their own pages (the operator_algebra reframe alone owes
`operator_inverse_family` / `operator_tensor_network` /
`coupled_block_operator` / `field_algebra` / `wavefront_cochain`). The
mechanical MOVE is a solved recipe; the sharp, build-invisible hazard is
in the AUTHORED header you wrap around it.

- **Locate by STABLE title, never the plan's line numbers** (they drift
  as the monolith is edited). `grep -n "<exact title>"` the four/N
  section titles + the FOLLOWING section's title (the exclusive upper
  bound). Then **prove contiguity**: `awk` the full-width `===` H1
  underline rows in the range and confirm ONLY your N titles appear —
  nothing else lives between the first and last section.
- **Inventory the traveling labels FIRST.** `awk` the range for BOTH
  `^\.\. _<label>:` (section/ref anchors) AND `:label:` (equation
  labels). ALL of them travel verbatim with the content — the plan
  usually names only the headline anchors ("and others"); the grep is
  authoritative. Labels are **path-immune**: inbound `:ref:`/`:eq:` from
  other docs resolve by NAME wherever the label lives, so they survive
  the move with zero edits on the consuming pages.
- **Slice PROGRAMMATICALLY (L-012), never hand-retype.** A Python splice
  reads the source, `block = lines[start-1:end]`, writes `header + intro
  + "".join(block)` to the new page and `prefix + pointer + suffix` back
  to the source. **Guard-assert the boundaries on the LIVE file**: first
  block line == the start `.. _<label>:`, the line after the block ==
  the next section's `.. _<label>:`, and `len(block)`. The verbatim
  block via `"".join` is transcription-safe.
- **⚠ THE TRAP — a Python f-string mangles LaTeX braces in the AUTHORED
  header, and `-W` is BLIND to it.** The moved block is safe (`"".join`,
  no interpolation); the header/intro YOU author is the risk. In an
  `f"""..."""`, ``:math:`A^{-1}` `` becomes ``A^-1`` (the f-string
  evaluates `{-1}` → the string "-1"); `\tfrac{1}{k}`, `{\rm loss}`,
  `\frac{a}{b}` corrupt identically. The mangled ``A^-1`` is **valid
  LaTeX math** — it renders (wrongly, as A⁻1 not A⁻¹), so `-W` never
  warns. This is a Cardinal-Rule-1 teaching defect on a NEW page.
  **Defense:** prefer plain concatenation (not an f-string) for
  math-bearing prose, OR escape every literal brace `{{ }}`; and ALWAYS
  grep the AUTHORED region before building — `grep -nE '\^-1|\{[^}]*\}'`
  over the header/intro lines only (the block's correct `A^{-1}` must
  not be touched) — then eyeball the rendered head. (Worked: split #1
  operator_inverse_family — the intro's two `A^{-1}` mangled to `A^-1`;
  caught by the visual head-read + a `grep 'A\^-1'`, fixed before the
  final build. The 1339-line verbatim block was untouched.) This is
  L-002 (build-blind correctness) ∩ L-012 (programmatic splice).
- **New page shape** (model `discretization.rst`): top `.. _<label>:` →
  over+under `=` title (size the bar with `len(title)` in CODE POINTS,
  L-009) → `.. contents::` `:local:` `:depth: 2` → a PROVISIONAL
  `.. dropdown:: Machine header — \`\`nexus-meta\`\` schema (PROVISIONAL)`
  `:color: muted` with a `code-block:: yaml` → a 1–2¶ intro that links
  UP to the parent (`:doc:`) and gives a semantic-TOC of the N sections
  (same-page `:ref:` to each moved anchor — guaranteed to resolve). The
  source's excised block becomes a **1-paragraph `:doc:` pointer**
  section (do NOT reuse the new page's top label). Wire the new page
  into the `index.rst` toctree (right after its sibling) AND add the
  intro `list-table` row.
- **Orphaned-HTML audit noise.** A `-E` build regenerates every LIVE
  page but does NOT garbage-collect HTML from renamed/deleted sources.
  A safety grep for `oldpage.html#<moved-label>` in the built tree WILL
  hit those orphans and look like a live stale ref. **Discriminate by
  "does the source `.rst` still exist?"** — absent ⇒ stale build
  artifact, out of scope (do not chase it); present ⇒ a real stale ref
  to fix. (Worked: `discrete_ordinates.html`/`loss_representations.html`
  orphans from prior renames carried the pre-move `operator_algebra.html#green-operator`
  href; their sources were long gone — irrelevant to the split.)
- **Split-#2 calibrations (extracting a SINGLE H1 section, vs #1's four).**
  (a) **Single-H1 → near-duplicate page/section title is ACCEPTED.** When the
  extracted block is ONE H1 section, the new page-top title and the block's own
  verbatim H1 are near-identical (page "Tensor-Network Decomposition of S_N
  Operators" over section "Tensor-Network Decomposition of SN Operators (Wave
  T)"). Do NOT "fix" it by rewriting the block's H1 — verbatim is the rule; the
  page-title + block-H1 pair is exactly split-#1's structure (page-title `===`
  followed by moved `===` sections) and builds `-E -W` clean. The two anchors
  differ (`operator-tensor-network` vs `wave-t-tensor-network`) and the titles
  aren't textually identical, so no duplicate-implicit-target warning.
  (b) **Inline `:sub:` in an over+under `=` title is a PROVEN pattern** — the SN
  book's own H1 is `Discrete Ordinates Method (S\ :sub:`N`)`. Mirror it
  (`... of S\ :sub:`N` Operators`), size the bar with `len(raw_title)` in code
  points (L-009; the role markup counts — 53 here), NO redundant `**bold**`
  (titles are already styled; the house convention uses `:sub:`/`:math:` bare).
  (c) **Directional-prose (L35) is mostly a NO-OP — fix ONLY phrases pointing at
  content that STAYS on the source page.** A grep of the block for
  above/below/earlier/later flagged 6 hits; FIVE were intra-block ("MA-Q1
  master condition above", "equations below", "coupling below" — all reference
  sibling subsections/eqs that TRAVEL with the block, so the relative direction
  is preserved on the new page) and one temporal ("later retired at commit").
  Exactly ONE was a genuine cross-page pointer ("operator-algebra types
  documented above" → the tensor-product TYPES that stay on the parent page) →
  rewrote to "documented in :doc:`/theory/foundations/operator_algebra`" (verify
  the target genuinely stays: grep the types on the source page OUTSIDE the
  moved range first). Calibration: in a self-contained deep-dive block, ~1-in-6
  directional phrases need a fix; classify each by whether its referent travels
  or stays, don't blanket-rewrite.
- **Split-#3 calibrations (g-adjoint → `operator_adjoint.rst`; the recipe's THIRD
  single-H1 run held perfectly, ZERO L35 fixes).**
  (a) **STRONGEST f-string-trap defense — author the header/intro AND the pointer
  as pure literals via the Write tool, not through ANY Python string.** The trap is
  a Python string layer corrupting LaTeX; eliminate the layer. `Write` the head to
  `/tmp/head.rst` and the pointer to `/tmp/pointer.rst` — the Write `content` is a
  pure literal (no f-string interpolation AND no backslash-escape processing, so
  `\frac`/`\nabla`/`\dagger`/`{-1}` are ALL safe — even safer than a raw
  `r"""..."""`). The Python assembler then ONLY does read+slice+concat+write; its
  sole string literals are file paths and the guard-assert boundary strings
  (`.. _g-adjoint:` + newline), NEVER math. Still run the mangle-grep (`A\^-1|G\^-1`)
  on the temp files as the gate (confirmed clean). Zero Python string layer over
  authored math = the trap cannot fire.
  (b) **A page's OWN other sections can `:ref:` the to-be-moved label — those become
  source→new-page cross-refs and are path-immune too.** Beyond EXTERNAL inbound refs,
  `grep <label>` the SOURCE page OUTSIDE the moved block: operator_algebra.rst had 5
  in-section `:ref:`g-adjoint`` (Key Facts + eigenvalue-posing region) that, post-move,
  resolve cross-page (HTML-audited: all 5 href `operator_adjoint.html#g-adjoint`, ZERO
  stale to the old location). No edits — same-file path-immunity works exactly like
  cross-file. Orphan gate: `grep -oE 'href="[^"]*oldpage.html#<label>"'` over the WHOLE
  built tree must be 0, then discriminate any hit by source-`.rst` existence.
  (c) **L35 gains two more no-fix categories: the TEMPORAL false-positive and the GONE
  referent.** Split #2 had 1 fix; #3 had ZERO. The 5 grep hits split: 3 intra-block
  ("(above)" / "following argument" / "below (≤3.6e-15)" — travel), 1 **temporal** ("an
  earlier version of this section" — a prior REVISION, not a spatial section; no fix),
  and 1 **gone referent** ("the reachability table below … see its retraction note" —
  the below/retraction-note direction still resolves to the intra-block Supersession
  `.. note::` that travels, but the "reachability table" itself no longer exists
  ANYWHERE on the page, removed in a prior edit). A gone referent is neither travels nor
  stays → OUT of L35 scope (L35 fixes only move-INDUCED breaks) → FLAG as pre-existing
  staleness (L-007), do NOT rewrite. Add temporal + gone to the per-phrase triage so
  the grep's false-positives are dispatched without a fix.
- **Split-#4/#5 calibrations (affine field algebra → `field_algebra.rst`; wavefront cochain →
  `wavefront_cochain.rst`; the brief's line-range OVERSHOT the genuine section BOTH times — the
  contiguity proof EARNED its keep, and this is now a STRUCTURAL, RECURRING pattern, not a
  one-off).**
  (a) **"Prove contiguity" must count ALL H1 `===` underlines, not just anchored ones — an
  ANCHORLESS sibling H1 inside the brief's range is INVISIBLE to the anchor-grep the brief
  author set the boundary with.** Split #4's brief "3844–4383 / up to the next anchor
  `wavefront-flux-cochain`" overshot an anchorless "The composite metric adjoint" H1; split
  #5's brief "3875–4434 / up to `coupled-block-operator`" overshot an anchorless "The inverse
  family" H1 — BOTH `:doc:` pointer stubs a PRIOR split left behind, sitting between the target
  section and the next `.. _label:`. A boundary set by "section anchor → next `.. _label:`"
  jumps clean over a stub (stubs carry no anchor). The `awk`-all-`===`-underlines proof (step
  1) CAUGHT both (2 H1s in the range, not 1). **This ALWAYS happens when the extracted section
  is immediately followed by a prior split's leftover pointer-stub H1 — endemic to a multi-split
  campaign; EXPECT it and run the count-all-H1 proof every time** (in split #5 the coordinator
  even pre-warned "confirm no anchorless sibling H1 inside the range"). The contiguity proof is
  NOT a formality, it is the gate for exactly this. (b) **Narrow the
  extraction to the GENUINE titled section when the trailing sibling's prose ties it to the
  SOURCE page.** The stub said "the adjoint face of the operator algebra **on this page**" —
  a bare directional coda (L35) that BREAKS if moved (the operator algebra stays on source),
  AND it is thematically the operator's adjoint, not the flux field algebra → LEAVE on source,
  move only `affine anchor → composite-adjoint H1`. The stable-title method (which L-026
  PRIORITIZES over line numbers) ends the section at the next H1; the brief's line-range was
  the author's anchor-scan estimate. FLAG the narrowing prominently — the coordinator set the
  wider range and reviews before commit. (c) **A three-way symbol collision folded into the
  split is a mini-L-011 done via the Edit tool BEFORE slicing (brace-safe), NOT in the
  f-string-risky head.** The reframe reserved `A` = full operator `L+C−S−B`; the block misused
  `A` for the sub-composite `L+C` (apply/solve headline + SI-increment machinery `M=A⁻¹(S+B)`,
  `Δψ=A⁻¹r̃`) AND the affine SPACE. Classify each `A`: full-op defect/solve → keep `A` (the
  residual `r=Aψ−q` was already full — reframe the headline TO match it); SI-increment
  machinery → genuinely the SWEEP `(L+C)⁻¹` (the honest spelling once `A` is the full solve —
  else Cardinal-Rule-1-wrong); affine space → `\mathbb{A}`. Apply on the SOURCE via Edit tool
  (exact-literal replacement is brace-safe) BEFORE the programmatic slice; the reconciled block
  then travels verbatim. Escape hatch (does the affine ARGUMENT require `L+C`?) did NOT fire —
  the torsor structure rests on flux-state geometry (no natural zero), independent of which
  operator connects the universes; test the ARGUMENT, not the symbol. (d) **Two guard-gotchas:**
  the mangle-grep FALSE-positives on `\mathbb{A}` + closing backtick (`mathbb\{A\}[^ ]` matches
  the `` ` `` after intact `{A}`) — the REAL gate is "any BARE `A^-1` (no braces)?" = ZERO; and a
  `.count(r"\mathbb{A}")` guard-assert counts OCCURRENCES not grep-LINES (a
  `\mathbb{A}\times V\to \mathbb{A}` line carries 2), so an `==10`-from-a-line-grep assert
  RED-flags its own miscount at 11 — a GOOD failure (asserts run before any write; content was
  right, the assert wrong). Fix the number, re-run; never loosen the assert to make it pass.
- **Split-#3/final calibration (coupled block operator → `coupled_block_operator.rst`; the NEW,
  build-INVISIBLE break class my #1/#2 closeouts MISSED — "named-source-page" consuming prose).**
  The inbound-ref audit is NOT just "does the `:ref:`/`:eq:` resolve?" — it ALWAYS does
  (path-immune, auto-repoints to wherever the label now lives). The break is in CONSUMING-FILE
  PROSE that NAMES the source page: ``see :ref:`X` in :doc:`.../operator_algebra` `` — after the
  move the `:ref:` links to the NEW page while the adjacent "in operator_algebra" sends the
  reader to the OLD page (now just a pointer stub). `-W` is BLIND (the ref resolves; only the
  prose lies) — a Cardinal-Rule-1 staleness the build never catches. **The L35 scan MUST extend
  beyond the moved block + source boundary to a whole-tree sweep of every consuming file** for
  ``:ref:`<moved-label>` … :doc:`<source-page>` `` (whitespace-FLATTENED / multi-line-aware —
  the "in" and the `:doc:` routinely wrap across two RST lines, so a line-grep misses them; use
  a `re.sub(r'\s+',' ',text)` python scan). Repoint the stale `:doc:` page-pointer to the new
  page; leave bare `:ref:`s (no page-pointer prose) and "see X **for the** …" phrasings alone.
  This session the sweep caught 2 in `sn/history.rst` (this split) AND a LEFTOVER from the
  committed split #2 (`sn/index.rst`: ``wavefront-flux-cochain` in :doc:`operator_algebra` ``) —
  fixed both, flagged the #2 leftover. LESSON: run this flattened consuming-prose sweep on EVERY
  split (the "all inbound refs resolve cross-page" HTML audit is necessary but NOT sufficient —
  it proves the LINKS work, not that the PROSE naming the old page is current). When the boundary
  is a real anchor (not an anchorless stub), the contiguity proof is clean — split #3 had exactly
  ONE H1 as the coordinator pre-scanned; still ran the count-all-H1 proof to confirm.
- **Split-#6 calibration (BISECTION: one page → TWO new children, source page DELETED — NOT the
  "extract-N-into-one-page, source-survives-as-a-`:doc:`-stub" model of #1–#5). Two structural
  facts the survivor-stub model never surfaces.** (a) **The source page-LABEL is RE-HOMED to the
  MAJORITY child, not left on a stub.** The source dissolves entirely, so its top
  `.. _<page-label>:` travels (verbatim) above the majority child's NEW title; every external
  inbound `:ref:` (audit the count — 9 sites here) auto-repoints to the child's new FILENAME with
  ZERO edits on the consumers (label path-immunity, exactly like a moved section anchor — HTML-audit
  that they land on `<majority-child>.html#<page-label>`, not the deleted file). The MINORITY child
  gets a BRAND-NEW label (absent pre-split). Single-homing is still the anchor-DEFINITION check
  (L-024): each `^\.\. _<name>:` == 1 tree-wide. (b) **Extracting a MIDDLE H1 while KEEPING its
  trailing H2 subsection → the H2 AUTO-REPARENTS to the PRECEDING H1.** When the brief keeps a
  subsection whose PARENT H1 is the one being extracted (keep "Cell-flattening invariant" H2, extract
  its parent "Cross-section convention" H1), removing the parent header re-attaches the orphaned H2 to
  the nearest preceding higher section (the previous chapter's H1 — "Derivation"). VERIFY TWO things:
  no title-level SKIP results (H1→H2 is legal and builds clean; H1→H3 would warn), AND the reparent is
  SEMANTICALLY intended (the kept subsection must belong under its new parent — the layout round-trip
  genuinely IS derivation content, even though it uses `sig_t` as the example). The brief's
  heading-LEVEL claim can be flatly wrong ("an H3 inside the Derivation chapter" for what the live
  `awk`/grep proves is an H1 SIBLING of Derivation) — L-001: the live underline rows are authoritative,
  and the programmatic guard-asserts on the EXACT header+underline strings catch a drifted boundary
  before any write. (c) **Blank-line glue at each splice follows the file's LOCAL convention** (this
  file: 1 blank between H2 siblings, 2 at H1 transitions) — after removing a middle span, DISCARD the
  bracketing structural blanks so the survivor spacing matches its new neighbours (char-identity applies
  to the moved CONTENT, L-013; the glue matches the destination, L-009). Boundary-crossing `:ref:`
  (a moved span references a label that STAYS, or vice-versa) needs NO syntax change — it silently
  becomes a cross-doc ref and resolves by name; HTML-audit it lands (here the XS page's
  `:ref:`…<sn-cell-flattening-invariant>`` → `indexing_and_layout.html#…`).

How to apply: title-locate → prove contiguity (count ALL `===` H1s) →
grep-inventory every label → programmatic guarded slice → author the
header WITHOUT an f-string over math (or escape braces) and grep it →
build `-E -W` to the unchanged baseline → HTML link-audit the inbound
refs land on the new page → **flattened consuming-prose sweep for
`<moved-label> … in :doc:<source-page>` and repoint the stale
page-pointer** → discriminate orphan artifacts by source existence.

---

## L-027 — The "relocate to page X" brief that the CLOSE READ reveals is ALREADY-fully-on-X → Cardinal-Rule-2 DE-DUPLICATE, not relocate+merge; plus the `ref.ref` caption-gotcha on an alias anchor, and reframe-consistency in a FOLD

A doc-cleanup brief can say "RELOCATE section S from page A to page B,
merging its additive parts in" — and the mandated close read of page B
reveals S is **already fully documented on B**, in equal or greater
detail. The scoping estimate (from a partial read) said "partial
overlap, ~1 additive artifact"; the CLOSE read (of every candidate
landing section) inverts it to "total duplication, ZERO additive."

- **When the content is already canonical on the target, the correct
  action is DE-DUPLICATION, not relocate+merge (Cardinal Rule 2).**
  Replace the source copy with a brief `:doc:` pointer that preserves
  the CONCEPTUAL BRIDGE (why the topic matters to the source page) and
  names the sub-topics, pointing at the canonical target sections; merge
  in NOTHING that already exists. Worked (#231 P4-T3): the "Boundary
  conditions as Wave-0/1 primitives" section on `operator_algebra.rst`
  was ~fully duplicated by `boundary_conditions.rst` — the G_α
  "primitives table" (its supposed one additive artifact) already lived
  there as the richer "SN realization map" list-table (§1794, with α=1
  fast-path columns + a bit-identity note the operator-algebra copy
  lacked); the rank-N eq, the Marshak example, the descriptor-vs-operator
  separation, AND the Wave-11/β1 predecessors were ALL already present.
  So: replace with a `:doc:` pointer, merge nothing.
- **The brief's "additive parts to MERGE IN" list is the SCOPING
  estimate, NOT ground truth — the close read overrides it. FLAG the
  inversion loudly** (the reviewer built the brief on the partial-read
  model). Report each "additive" item as "already at §X of the target,
  deduped" with the section reference, so the reviewer's close review
  can restore any specific piece. This is the T-relocate analogue of the
  split-#4 anchorless-sibling discovery: the contiguity/close-read gate
  is what catches the brief's wrong structural assumption.
- **Carry the moved section's `:ref:`-able labels as ALIAS anchors onto
  the canonical target content** (zero-inbound-ref labels still get
  carried "for outbound-ref integrity" per the brief) — BUT an EQ-label
  (`:label:`, used via `:eq:`) CANNOT be aliased onto a different eq;
  drop it if its eq is duplicated and it has zero `:eq:` refs (flag the
  drop). Std `.. _` anchors alias fine.
- **⚠ NEW `-W`-CAUGHT warning class — `ref.ref` "A title or caption not
  found".** A `.. _label:` placed before a **paragraph** (or any element
  with no title/caption) makes a BARE `:ref:`label`` FAIL under `-W`
  ("Failed to create a cross reference. A title or caption not found").
  Unlike a dead code-xref (plain-text, L-002-silent), THIS one IS gated
  by `-W`. Two fixes: (a) place the alias anchor before a TITLED or
  CAPTIONED element (a section title, or a ``.. list-table:: Caption`` /
  figure — then the caption becomes the link text), OR (b) use
  EXPLICIT-text `:ref:`link text <label>`` (resolves regardless of what
  the label precedes). An anchor before a SECTION TITLE is already safe
  (bare ref gets the title); an anchor before a paragraph/list needs (a)
  or (b). (Worked: `bc-descriptor-tree-vs-operator-tree` before §2419's
  title was fine; `bc-tensor-primitives` before a paragraph warned — I
  moved it before the captioned list-table AND made the pointer
  explicit-text.)
- **Reframe-consistency applies to a FOLD, not only a split/move.** When
  folding a paragraph into a keeper section, the SAME overloaded-symbol
  reconciliation applies: T1's `loss_minus_gains(psi) = A.apply − Σ
  g.apply` used `A` for the sub-composite `(L+C)` (gains = S,B) — the
  pre-reframe collision. Verify against LIVE code (`iteration.py:903` +
  its docstring "the matvec IS the honest `(L+C−S−B)·ψ`"), then fold
  reframe-consistently: spell the sub-composite `(L{+}C)` explicitly and
  identify the result as the full `A = L+C−S−B` applied. A fold is a
  move; a move re-exposes every reconciliation the source had.
- **A stale deferred-follow-up ISSUE-tag is distilled by git, not
  guessed.** T4's "Deferred follow-ups: #260 …, #261 (core relocation…)"
  — `git merge-base --is-ancestor <#261-commit> HEAD` proved #261 landed
  → drop it, keep the still-open #260 (singularize "follow-ups"→
  "follow-up"). Keep the DESIGN RATIONALE ("considered-and-rejected: R/M
  are rank-changing einsums, not valid tensor factors") — only the dated
  tracking tail distills.
- **⚠ Dropping a duplicated EQ-label ALSO drops its `.. vv-status:` — check the
  V&V-matrix consequence and MOVE the status to the survivor.** A de-duplicated
  concept can carry TWO `:label:`s with DIFFERENT V&V status — a `documented`
  twin and an `orphan` twin (same math). Dropping the `documented` one (because
  it is the duplicate) silently DEMOTES the concept to `orphan` (untracked), and
  `-W` is BLIND: the orphan-equation gate is a `docs/verification/matrix.rst`
  REPORT auto-regenerated by `conf.py`'s `generate_matrix` hook, NOT a
  build-breaking check. If the concept is genuinely documented-only (a
  declarative/definitional eq with no numerical result to test — e.g. a
  BC-algebra decomposition), ADD `.. vv-status: <survivor-label> documented` to
  the surviving canonical eq so the accounting is preserved through the de-dup.
  Worked (#231 P4-T3, caught by the main agent in review): dropping
  `bc-rank-n-as-sum-of-products` (documented) left only the orphan
  `bc-rank-n-tensor-decomposition` — fixed by marking the survivor `documented`.
  The retirement-audit's "handle the V&V edges" includes the vv-status
  directive, not just the `:eq:` refs.

How to apply: on a "relocate to X" brief, CLOSE-read every candidate
landing section on X FIRST; if the content is already there, de-dup by
`:doc:` pointer (merge nothing) and FLAG the inversion; carry `.. _`
aliases onto the canonical content (before a titled/captioned element,
or use explicit-text refs to dodge `ref.ref`); a fold is a move — apply
the reframe reconciliation; distill stale issue-tags by `git
merge-base`, keeping the design rationale; when dropping a duplicated
eq-label, MOVE its `.. vv-status:` to the survivor (the `-W`-blind
orphan-demotion). **Mechanical:** a full `-E`
rebuild here EXCEEDS the 120s foreground cap — use `run_in_background`
for the authoritative gate (a foreground poll-loop gets SIGTERM'd at
2 min, killing the build at the final line before "build succeeded"
prints, so the summary never lands even though zero warnings were
raised).

---

## L-028 — The Key-Facts↔changelog metadata-RELOCATION (L10 Sphinx-as-brain ≠ Sphinx-as-history): strip campaign-provenance from the high-traffic invariant section INTO the page-bottom changelog — move it, don't lose it

A theory page's **Key Facts** (highest-traffic section) accretes
campaign-provenance clauses over a long refactor — each invariant bullet
tagged "(Wave O step B.5.2, Issue #208, ``6ef5063``, 2026-06-03)". The
page-bottom **Development-history changelog** is where that provenance
belongs (L10). The task is a RELOCATION, not a rewrite: the invariant
stays in Key Facts, the campaign-metadata moves to the changelog.

- **The strip-list is EXACT and narrow: commit hashes, round/Wave/Phase
  labels ("Wave O step X", "Phase 5a", "carve P4", "S6.4(f)"), branch
  names, landing dates.** KEEP everything else in place — the invariant
  statement, the production-formula **eq-labels + their ``.. vv-status:``
  comments**, every ``:ref:``, every **active gotcha** ("F never enters
  A", "NOT deprecated", the ``(N,ng,nx,ny)`` convention), AND numerical
  evidence (an "18.3× shrink" is a datum, not campaign-metadata). Verify
  the KEEP-set survives with a grep after (eq-label + vv-status count
  unchanged; gotcha phrases present).
- **Issue-# refs (``#208``) are NOT in the strip-list — KEEP them inline
  and FLAG the decision.** They are lightweight GitHub cross-refs a
  reader values, distinct from dated/hash history; the changelog's Issue
  column carrying the same number is acceptable redundancy. (If the
  reviewer wants them stripped too, that's a one-line follow-up.)
- **"Move it, don't lose it" gates what you may strip: NEVER strip a
  DATED milestone whose provenance has no changelog home.** Map each
  bullet's provenance to a destination FIRST: (a) an EXISTING changelog
  row, (b) one of the NEW rows you add, or (c) NONE. For (a)/(b) strip
  freely. For (c) — a dated milestone missing from the changelog and NOT
  in the sanctioned new-row set (worked: the coupled-block 2026-07-12
  bullet) — KEEP its provenance inline and FLAG it (recommend a row);
  stripping it would DELETE the date/branch, violating the relocation
  principle. A round-label whose milestone IS covered elsewhere (the
  "since stencil-assembly 2b" assembly-axis note, already a changelog
  row) strips cleanly — nothing dated is lost.
- **Verify every commit hash is a HEAD ancestor before citing it in a
  new row** (``git merge-base --is-ancestor <hash> HEAD``); pull the
  ``When`` date from ``git show -s --format=%cs <hash>`` (don't trust the
  bullet's inline date — get it from git). New rows go in reverse-chron,
  matching the existing list-table's 4 columns (When / milestone / Issue
  / Where), iteration-rates omitted per the changelog's own preamble; the
  ``Where`` column carries the stripped hashes (``main (Phase 5a,
  ``hash`` / ``hash``)``). The new-row milestone text is the SAME
  invariant you kept in Key Facts, distilled — so Key Facts and the
  changelog agree, neither is lost.
- **A metadata-strip pass routinely surfaces an ORTHOGONAL staleness in
  the same bullets — FLAG, don't fix (out of scope).** Worked: a
  Key-Facts affine bullet still wrote the affine space as bare ``A``
  while the split-out ``field_algebra.rst`` deep-dive had renamed it to
  ``\mathbb{A}`` (the split-#1 reconciliation never propagated to the
  summary bullet). That is a reframe-consistency fix, not metadata —
  report it as a found-in-passing defect, leave it for a scoped pass.

How to apply: map each bullet's provenance to a changelog destination
(existing row / new row / none); strip only hashes/round-Wave-Phase-
labels/branches/dates that HAVE a home; keep invariants + eq-labels +
vv-status + gotchas + issue-refs + numerical data; keep-and-flag any
dated milestone with no changelog home; git-verify every cited hash is a
HEAD ancestor and git-source its date; grep the KEEP-set is intact;
flag orthogonal staleness found in passing.

---

## L-029 — The additive "surface the taxonomy up front" framing pass: verify the gap is REAL (a sibling taxonomy present ≠ THIS one surfaced), then PREVIEW + `:ref:` to the SSOT — never a twin table

When a scoping pass finds a governing taxonomy (a 3-way type partition,
a classification law) stated only MID-STREAM — in the 2nd of the
sections it frames, so a linear reader meets the 1st with no roadmap —
the fix is ONE short additive framing section, not a reorder. (Meta:
the scoping earned its keep by KILLING the big mechanical move — the
flagged apply/solve↔streaming dependency-INVERSION led the coordinator
to KEEP-the-early-section and ABANDON the relocation, leaving this one
high-value additive item. A load-bearing dependency flag can collapse a
multi-move plan to a single preview.)

- **Escape-hatch FIRST — verify the gap is REAL, and MATCH the specific
  taxonomy, not "is there any framing."** The trap: the upstream
  "already there" candidate (Key Facts) stated a DIFFERENT partition —
  the Representation×Role CARRIER taxonomy — which is NOT the
  Operator/Kernel/Functional CODOMAIN partition. The presence of A
  taxonomy is not the presence of THIS one; two coexist in the same
  page. Read/grep for the EXACT partition before concluding it's
  un-surfaced. If it IS already surfaced clearly upstream, FLAG that
  (decline to add ceremony) rather than force the section — the
  coordinator's explicit escape hatch, and Cardinal Rule 2 (no
  ceremony).
- **SSOT-vs-twin — the preview is a prose ROADMAP + `:ref:`, NEVER a
  copied table.** Naming the three arms in prose and `:ref:`-ing the
  canonical codomain-partition table is a pointer that cannot harmfully
  drift; COPYING the table into the preview is the Cardinal-Rule-2 twin.
  When offered "a light RELOCATION of the existing framing" as the
  elegant alternative (SSOT over preview+ref), prefer the additive
  preview when relocation would (a) touch a section you were told NOT to
  modify, or (b) orphan the framing's double-duty — here the partition
  statement is ALSO the Functional section's OWN opening paragraph, so
  relocating it guts that intro. The SSOT stays where it is; the preview
  points at it.
- **`:ref:` the SSOT SECTION anchor, verified live — not a sub-anchor
  that doesn't exist.** The canonical table lived inside the Functional
  section under an UN-anchored subsection; point `:ref:` at the section
  anchor that HOSTS it (``functional-category``), not the subsection.
  Intra-doc `:ref:` warns if dangling (L-002), so the clean `-E -W`
  build confirms resolution, and the HTML link-audit
  (``href="#functional-category"`` with the target title as link text,
  not a plain ``<code>``) confirms it rendered as a hyperlink.
- **Baseline is FRESH, not frozen.** Measure the `-E` baseline THIS
  session before crediting "count-unchanged" — the AGENT.md "1 warning
  (mesh.py ``:paramref:``)" note was STALE (true baseline 0 on
  ``docs/sn-doc-architecture``; baselines drift 9→1→0). An additive
  preview with no new `:label:`/citation and one live `:ref:` is
  provably warning-neutral; the pre/post SET-diff (both 0) proves it,
  and a pure-additive H1 needs the ``len(title)``-sized underline
  (L-009) — 28 code points here.

How to apply: verify the taxonomy is genuinely un-surfaced upstream
(match the SPECIFIC partition, not any framing; flag-and-decline if it's
already there) → author a prose-roadmap preview naming the arms +
`:ref:` to the SSOT section anchor (never copy the table) → prefer the
additive preview over relocation when relocation touches a frozen or
double-duty section → point `:ref:` at a live anchor and HTML-audit the
rendered hyperlink → measure the `-E` baseline fresh and diff the
WARNING/ERROR/CRITICAL set pre/post.

---

## L-030 — Additive `:label:` backfill on a derivation-mirror page: the skeleton is already labeled, so BARE dominates; fill only 5 recognizable gap classes

The #231 G1 batches insert descriptive `:label:` under unlabeled
`.. math::` on the literature-mirror corpus (`docs/theory/references/*`).
PURE ADDITIVE — a `:label:` line only, at option indent, immediately
after the `.. math::`; touch no content/prose/headings. Recount per
file yourself (`grep -cE "^\s*\.\. math::"` vs `grep -cE "^\s+:label:"`),
never trust the brief's count. The adjudication discipline:

- **On a page that mirrors a published derivation CHAIN, the checkpoints
  are usually ALREADY labeled** (governing eq, named kernels/Green's/
  resolvents, final boxed results, paper-numbered key eqs). So the
  unlabeled residue is dominated by TRUE intermediates that correctly
  STAY BARE — "label the skeleton, not every vertebra." The BARE classes
  (each recurred dozens of times): step-to-step algebra between two
  labeled checkpoints; substitution/change-of-variable steps; an
  immediate RESTATEMENT of an already-labeled eq (often cross-ref'd with
  `:eq:` right there — a dead giveaway); a COMPANION definition under a
  labeled eq (the `K_ij = …` under a labeled Nyström eq; the `where
  T_vol is …` operator/kernel defs under a labeled operator form; the
  `where B_LR = …` under a labeled closure); a STANDARD special-function
  integral representation (E₁, Ki₁, Ki₃ — labeled once via its own
  derivation, then restated bare); a test-gate closed-form realization;
  a cost/complexity model; a geometric parametrisation (chord length,
  ρ_max, law of cosines); a SCHEMATIC contrast (`∝`, a template with a
  text "(geometry factor)" placeholder); a numerical RESULT / sanity-
  check evaluation (`k_eff = 1.00421`, a vacuum-limit = 4); and
  deferred-investigation or falsified-heuristic components. Ratios land
  low BY DESIGN — trajectory_resolvent was 2 labeled / 31 bare,
  peierls_nystrom 8 / 88; that is the guidance applied correctly, not
  under-labeling.

- **The GENUINE gaps are a small, recognizable set — fill these:**
  (a) **Governing-eq parallel** — a page's BTE / 3-D transport eq when a
  SIBLING page already labels its own (galerkin_spectral `-bte` →
  mirror on fn_method `fn-method-bte`, peierls_nystrom
  `peierls-transport-equation-3d`); label for cross-page consistency.
  (b) **Unlabeled named-object definition** the page/corpus uses BY NAME
  but never labeled — escape probability `P_esc`, a continuum dispersion
  function `λ(ν)`, a discrete pseudo-eigenfunction, an L² inner product /
  Galerkin-orthogonality principle. (c) **Geometry-parallel gap** —
  cylinder+sphere carry `-nystrom`/`-row-sum-identity` but slab doesn't
  (fill the slab), or vice-versa. (d) **Sibling-parallel result** —
  sphere T-matrix labeled but cylinder/slab aren't; slab labels BOTH P
  and G mode formulas but cylinder only labels P (fill the cylinder G) —
  so grep finds the whole set. (e) **Paper-numbered key eq** in the
  page's ESTABLISHED `<page>-eqNN` family (`galerkin-spectral-eq3`
  joining `-eq4`; `singular-eigenfunction-eq47` joining `-eq46/-eq54`;
  `wm72-eq21d-normalization` joining `wm72-eq30/31/32`). The `-eqNN`
  domain-form-with-number IS the page's precedent, NOT the forbidden
  bare-positional `eq7` — match the family, don't invent a scheme.

- **Mechanics that bit this task (all avoidable):**
  - **zsh does NOT word-split an unquoted `$var`** (bash does). A
    `for n in $names` uniqueness loop silently ran ONCE on the whole
    concatenated string and printed a false "0 collisions". Use an
    explicit literal list (`for n in a b c …`) or `${=names}`. The
    per-name `grep -rn ":label: $n\$" docs/theory` gate (must return
    NOTHING before finalizing) is mandatory — labels are PROJECT-GLOBAL,
    a dup fails `-W`.
  - **Edit `old_string` must match LIVE bytes, not a prior render.** A
    block read earlier as `\mathrm dy` was actually `\mathrm{d}y`
    (braces); the edit missed → re-read the exact lines and retry.
  - **An aligned block (`A &= … \\ B &= …`, NO blank line) is a SINGLE
    align environment** — one `:label:` labels the whole block, SAFE.
    The "don't label multi-equation blocks" rule targets only
    BLANK-line-separated sub-equations (none appeared in this corpus).
  - **A list-nested `.. math::` (2-space indent) takes its `:label:` at
    5 spaces** (directive-indent + 3), not 3.
  - Do NOT run sphinx-build (main agent runs the gate once post-batch).
    Verify instead by recount (math count UNCHANGED ⇒ no block corrupted;
    label count up by exactly N) + the per-name uniqueness grep.

How to apply: recount yourself → learn each page's existing label
FAMILIES first (`grep -nE "^\s+:label:"`) → read each unlabeled block's
surrounding prose and default BARE unless it hits one of the 5 gap
classes → name in the page's family word-order → per-name uniqueness
grep (with a shell that actually word-splits) → apply, then recount to
confirm math unchanged / labels +N.

---

## L-031 — The docutils→bibtex citation migration: whitelist-scope the swap (auto-skips non-keys), indentation-key the block remover (preserves nested notes + footnotes), 3-signal-gate every heading removal

A corpus-wide `[Key]_` → ``:cite:`Key``` migration (sphinxcontrib-bibtex,
#231 Phase G2: 233 swaps + 78 def-blocks across 46 files) is mechanical
but has four traps a blanket regex walks straight into.

- **Scope the swap by a WHITELIST built from `refs.bib` keys (+ ruled
  consolidation aliases), NOT a blanket `\[\w+\]_` regex.** A non-key is
  simply never a replacement target, so the whitelist auto-skips every
  pseudo-site with ZERO line-number logic: `[A]_{:,j}` matrix notation,
  a ``[Foo1234]_`` syntax example, and footnote uses `[#name]_` (the `#`
  is outside `[A-Za-z0-9_]`). Literal `str.replace('[K]_', ':cite:`K`')`
  preserves surrounding punctuation exactly and can NOT touch the def
  line (`.. [K] ` ends in a space, not `]_`), so swap-order is irrelevant.
  Consolidations map alias→canonical (`[PS1982]_`→``:cite:`Pomraning…```);
  the alias key is REMOVED from refs.bib by the bib owner, so emit only
  the canonical and verify zero leaked ``:cite:`AliasKey``` after.
- **Def-block removal MUST be INDENTATION-based, not blank-delimited.** A
  citation body can carry an INTERNAL blank line — a nested `.. note::`
  admonition inside the `.. [Key]` block — so a "consume until blank"
  remover ORPHANS the note (a real dry-run hit). Consume the `.. [Key]`
  line + every following line indented deeper than the marker (internal
  blanks folded in via lookahead: a blank whose next non-blank is
  deeper-indented is internal), stopping at the first dedent-to-base.
  KEY THE REMOVER TO THE WHITELIST so footnotes (`.. [#name]`) and any
  non-citation `.. [x]` survive — a `[^\]]+` remover eats footnotes.
- **Footnotes are a DIFFERENT docutils construct — preserve both halves.**
  Always DRY-RUN-categorize every `[x]_` use AND every `.. [x]` def
  against the whitelist first, printing UNKNOWN keys + STRAY brackets;
  the dry-run is what surfaces the `.. [#name]` footnote family (3 here)
  and any typo'd key before a single byte is written.
- **Emptied-"References"-section removal needs a 3-SIGNAL cleanliness gate
  (grep BEFORE removing any heading):** (1) is `autosectionlabel` enabled
  in conf.py? — if OFF, a bare `References` heading is NOT a cross-ref
  target; (2) any inbound `:ref:` to the section's explicit `.. _anchor:`
  (grep tree-wide — here the sole `.. _bib-*:` citation anchor had zero
  referrers → safe to drop with its citation); (3) directional prose
  ("listed below", "the references section"). All three clean ⟹ REMOVE
  heading+underline+preceding-blanks (asserted script: strip trailing
  blanks, assert underline, assert heading text, pop). Referenced ⟹ keep
  + one pointer line to the bibliography page. MIXED section (docutils
  defs + a plain-text further-reading bullet list / "Internal references"
  subsection) ⟹ KEEP the heading, delete only the def blocks.
- **A NOTE describing the RETIRED docutils cross-doc citation mechanism is
  Cardinal-Rule-1 obsolete post-migration — remove it.** ("Citations
  shared across pages resolve cross-document via Sphinx's docutils
  citation index; only local citations defined here" is now false under a
  central refs.bib.) It housed the ``[Foo1234]_`` skip example, which
  correctly vanishes with it — SKIPPING a pseudo-site (don't swap) is not
  PRESERVING it forever; report the removal as a per-page decision.
- **keylabel style ⟹ the migration is INVISIBLE to readers.**
  `bibtex_default_style = 'keylabel'` renders the label AS the key, so
  ``:cite:`Hebert2009``` displays `[Hebert2009]` — character-identical to
  the retired bracket label. State this in the bibliography page's lead.
- **Python docstring-only constraint is git-diff-verifiable.** A `]_`
  operator can't appear in executable Python, so citation uses/defs live
  ONLY in docstrings — but still GATE it: `git diff <pyfiles> | grep
  '^[-+]' | grep -vE '^[-+]{3}'` and confirm every changed line is
  docstring/reference text, zero `def`/`import`/logic. Confirm the one
  math-notation skip file (`operator.py` `[A]_`) has an EMPTY diff.
- **The new bibliography page is NOT a reference SOLVER.** `.. _anchor:`
  above title, lead para (single citation home; entries in refs.bib /
  Zotero upstream; keylabel; per-page lists retired), then a BARE
  `.. bibliography::` (only cited entries render). Grep for a pre-existing
  `.. bibliography::` (a second one warns) + label collision first. Place
  it in its OWN labelled toctree subsection — dropping it into a
  reference-SOLVER toctree miscategorises it; give it a `-` subsection
  under "Pages" (size the underline in code points, L-009).

How to apply: dry-run-categorize uses+defs against the refs.bib whitelist
and report unknowns; literal-swap only whitelist keys (aliases→canonical);
indentation-key the whitelist-scoped block remover (nested notes +
footnotes survive); 3-signal-gate every heading removal (autosectionlabel
+ inbound `:ref:` + directional prose); remove stale mechanism-notes;
verify `.py` diffs are docstring-only; give the bibliography its own
toctree subsection. Full inline-output identity comes free from keylabel.

---

## L-032 — P10 `:label:` re-namespacing (a label follows its heading's ruling): the SELF-DESCRIPTION oracle, the section-vs-eq-label scope split, and the delimiter-anchored replace that survives a prefix-overlapping sibling

Deferred from a Phase-F HEADING retitle campaign: labels carry tree-wide
`:ref:` blast radius, so the anchor rename is its own pass. The rule: a
label follows its heading's design/record ruling (design-named heading →
RENAME the label to the heading's vocabulary, keep the page prefix;
record/charter heading → KEEP). The disciplines that pass taught:

- **The section's OWN self-description is the strongest oracle** — stronger
  than the label name, and it resolves the brief's genuine hedge. A `sn-282-*`
  family "may be a #282 record" per the label, but the section literally said
  "this section is the **resolution chapter** … those [other] sections are
  preserved as **the record**" — a page drawing the design/record line WITH
  ITSELF ON THE DESIGN SIDE. Combined with (a) the Phase-F map retitling its
  top heading to a design name + "P7's own worked example", (b) all subsection
  headings design-named, (c) it living on the DESIGN page while the record
  lives in the charter chapter — the verdict was RENAME, unambiguously. A
  charter page states its own status too: `curvilinear_numerics.rst` opens
  "This chapter is Part B's **campaign record**" → KEEP-all (incl. its
  issue-styled `sn-issue-196-*` anchors — an issue number inside a charter
  chapter keeps).
- **When the labels form a design FAMILY, P10 is FULL — section anchors AND
  equation labels drop the campaign token together; excluding eq labels is the
  WRONG default.** My pass-1 instinct (rename anchors, EXCLUDE eq labels because
  they carry V&V-matrix + `.. vv-status:` + `:eq:` weight, "flag the mixed
  namespace") was OVERRULED: a `sn-r12a` section anchor beside a
  `sn-282-r12a-predicate` eq label is a *two-spellings* state, and the standing
  naming-consistency rule (a family follows ONE pattern; fix off-pattern members
  in the SAME change) forbids shipping it. The correct pass unifies the whole
  family (`sn-282-*` sections + eq labels → `sn-direct-seed-*`; `wave-t-*` →
  `tensor-network-*`). The eq-label V&V weight is NOT a reason to defer — it is
  a **gate**: a documented-only eq label (silent-class grep of `orpheus/`+`tests/`
  = ZERO, so no `@pytest.mark.verifies("label")`/`catches` edge) renames
  mechanically (the `matrix.rst` auto-regens, L-008; `.. vv-status:` directives
  move with the label). ONLY a verifies-target eq label (a silent-class HIT)
  would orphan a test edge (L-003 phantom-verifies) — THAT is the one to
  flag/defer, and you don't edit `tests/` regardless (report the hit). So: run
  the silent-class grep FIRST; empty ⟹ full unification is safe and MANDATORY;
  hits ⟹ rename docs, report the test edge. (The lone brief-named eq residue
  `region-areas-pincell`→`-pin-cell` was always in scope.)
- **A section anchor can be a PREFIX of a sibling equation label**
  (`sn-282-r12a` ⊂ `sn-282-r12a-predicate`). A bare string replace corrupts
  the sibling. Replace ONLY fully-delimited forms — `.. _X:` · `` `X` `` ·
  `<X>` (the `:ref:`txt <X>`` angle form) · `.. vv-status: X ` (trailing
  space) — where the char after `X` (`:` `` ` `` `>` space) can't match a
  longer label. Do the rename with a script that reports a per-file, per-form
  count and ASSERT it against the pre-computed inbound tally; then a
  corruption grep (the sibling still present, `sn-r12a-predicate` = ZERO)
  proves the delimiters held.
- **Collision + genericness force per-label disambiguation off the clean
  token drop.** Dropping the campaign token gives `sn-282-gotchas`→`sn-gotchas`
  — but `sn-gotchas` already anchored the index page, and bare
  `sn-numerical-evidence` is one of three sibling `*-numerical-evidence`
  families. Disambiguate with a design family prefix tied to CODE vocabulary
  (`sn-direct-seed-*`, from the test file `test_282_direct_seed_*`), grep-check
  every proposed new name for 0 collisions BEFORE the script runs. And once a
  disambiguator prefix exists for TWO members, the naming-consistency rule
  pulls the WHOLE family onto it (a 7-bare/2-prefixed split is off-pattern) —
  the clean endpoint is `sn-direct-seed-*` for all 9 sections + all 5 eq labels,
  not a per-label heading-mirror (see the FULL-P10 bullet above).
- **The deterministic gate is the grep, not the build** (L-002): a cross-doc
  `:ref:` to a renamed-away anchor renders plain-text with NO `-W` warning.
  Proof-by-construction = {no OLD delimited form survives anywhere} ∧ {every
  NEW anchor exists exactly once}; confirm with a rendered-HTML `href=` audit
  on one cross-doc ref (`api/*.html` → `foo.html#new-anchor`). The clean `-W`
  build only catches intra-doc dangling. A code-string ref in `orpheus/`
  (`f"…§peierls-phase5-retreat"`) is silent-class — REPORT file:line + new
  label for the main agent (you don't edit `orpheus/`/`tests/`).
- **Prose two-spelling harmonization rides the same pass** (Surface 3): after
  renaming an anchor, event-name prose ("the Phase 5 retreat" → design name
  "the continuous-µ retreat") harmonizes design-first, KEEPING the historical
  tag where it carries record value ("… (Phase 5's terminal decision)") and
  leaving pure-provenance/file-path mentions ("Phase-5 Round-3 provenance",
  `diag_phase5_*.py`) untouched — fix ONLY where two-naming reads as two
  events. Match the canonical code-point (event-name "continuous-µ" = literal
  µ U+00B5 per the heading; the math object is `:math:`\mu``) — source it from
  the file, don't retype. A consistent issue-TAG beside an established design
  name (loss_representation's "the sweep-inverse-contract discharge (#284)")
  is NOT two-naming → WON'T-FIX with a one-line justification.

How to apply: read each candidate's heading TODAY + the Phase-F batch map +
the section's self-description; charter pages KEEP-all; design headings →
rename anchors (not eq labels — flag the mixed namespace) via a
delimiter-anchored counted script; grep-prove no OLD form survives + new
anchors exist; HTML-audit one cross-doc ref; harmonize event-name prose
design-first; report silent-class `orpheus/`/`tests/` hits.

---

## L-033 — Code-prose rebalance (#231 Phase 2): an operator file's teaching is ALREADY TWIN — expect ZERO MOVED, verify by token-invariance, and keep the CONTRACT tier by the "wrong-simplification guard" test

The pilot (P2-A, `transport/operators/scattering.py`, 73%→63% prose,
docstring 1127→721, comment 196→121) established the calibration for
rebalancing a "documentation-with-code-attached" operator file. The
classify-each-block-into-one-verdict rubric
(CONTRACT/TWIN/MOVED/HISTORY/COMMENT-cut) has a dominant outcome on
operator files that the instinct fights:

- **The operator-algebra book is EXHAUSTIVELY COMPLETE → expect ZERO
  MOVED.** On a heavily-prose operator file the reflex is "some design
  rationale must be book-worthy-but-unwritten". It almost never is. Three
  concepts that FELT unique to the code — the forward-fast-path-vs-adjoint-
  frame asymmetry, N2N-as-a-distinct-moment-operator, the foldable/σ_r
  split + the #215 trap — were EACH fully TWIN after grepping the landing
  chapters (`adjoint.rst` §sn-scattering-adjoint-source; `operator_algebra.rst`
  §integral-kernel-category / §scattering-as-tensor-product-sum;
  `slab_one_group.rst` §si-sigma-r-fold-mismatch + `loss_representation.rst`
  §loss-rep-removal-sigma). **Discriminator: grep the landing chapter for
  the concept BEFORE assuming novelty.** Budget the operator-file batches
  (fission, streaming, boundary, multiplication, isotropic_scattering) as
  TWIN-CUTTING, not MOVED-writing. If you think you found a MOVED, grep
  harder first — Cardinal Rule 3 means the theory shipped with the code.

- **CONTRACT-vs-{HISTORY,TWIN} has three sharp discriminators the pilot
  named** (these are the reusable judgment rules, not the file specifics):
  (a) **A keep-vs-retire decision on a currently-orphan symbol is CONTRACT,
  not HISTORY** — even phrased historically ("Deliberately retained W-F,
  user steered keep"). It is a live constraint: the arm is an *intentional*
  orphan kept for a named future consumer (an OPEN issue), and a naive
  retirement audit would delete it as dead. Keep the keep-decision + the
  open-issue rationale; cut only the date/steering provenance. (b) **A
  `⚠ LATENT … TRAP` / "do NOT" imperative is COMMENT-keep even when its
  EXPLANATION is TWIN** — the derivation goes to a `§`-pointer, but the
  imperative + the falsifying number (46–56 %) + the tracking-issue
  pointers (#2/#215) stay inline, because a future within-group-solver
  editor reads THIS file, not the theory page. (c) **A type-annotation-
  choice rationale that guards a plausible wrong "simplification" is
  CONTRACT** — "returns the concrete `OperatorSum`, not the bare
  `LinearOperator` erasure, so `apply_transpose` stays visible to the
  checker" prevents a modifier from silently breaking a consumer by
  "tidying" the return type. General test for the CONTRACT tier: **"would a
  competent modifier who never leaves the file do the wrong thing without
  this line?"** — if yes, it is CONTRACT regardless of how history-flavored
  the prose is.

- **HISTORY-cut only after confirming the fact is in the record.** The
  module-head Wave-D-extraction narration ("Per Cardinal Rule 2 this
  lifts… bit-identical extraction… SNSolver retains thin delegators") is
  HISTORY *because* slab_multigroup.rst 439–444/578–582 already carries it
  AND the delegators verifiably still live (solver.py:1884/1892/1921).
  Verify BOTH — the record home and the live truth — before cutting; a
  HISTORY claim that is novel-and-recorded-nowhere becomes a
  Development-history dropdown MOVE, not a cut.

- **Verification discipline unique to this task class:** (1) the edits are
  docstring/comment-ONLY, so PROVE zero code change by a **non-string/non-
  comment token comparison vs HEAD** (`tokenize`, drop COMMENT/STRING/
  NL/INDENT…) — 2397==2397 is stronger evidence than `pytest --collect-only`
  (which also passed, 6652 unchanged). The `code lines` count (484→484 via
  `ast`) corroborates. (2) **The sphinx gate is CONDITIONAL on automodule
  status — check it FIRST (`grep -rn "automodule:: <module>" docs/`).** A
  not-`automodule`'d file (P2-A scattering.py) has build-invisible docstrings →
  SKIP the multi-minute build (say why). But an `automodule`'d file (P2-G
  streaming/augmented_mesh/boundary, all in `api/discrete_ordinates.rst`)
  RENDERS its docstrings → the `-E -W` build gate is MANDATORY (capture the
  baseline BEFORE editing; acceptance = W/E/C set unchanged). `:noindex:` does
  NOT exempt it — `:noindex:` only makes the module xref-invisible (L-002); the
  docstrings still render and malformed markup still warns. Two automodule-safe
  moves when trimming: (a) NO `.. math:: :label:` in any of the three docstrings
  → cutting math blocks orphans no `:eq:` (grep-confirm the file's `:label:`
  count is 0 first); (b) KEEP section-title underlines VERBATIM (over-long is
  allowed) — trim only the prose body under a heading, never resize the
  underline, or you risk "Title underline too short". (3)
  **Pointer form = literal greppable `docs/theory/<part>/<file>.rst §<label>`**;
  `§` may point at an EQ-label (`:label:`) when no co-named section anchor
  exists (`mg-inscatter-source`, `sn-scattering-adjoint-source`) — it is a
  human marker, not a rendered role (the file isn't automodule'd), so it
  resolves via grep. Gate every label with `grep -E "^\.\. _X:|:label: X$"`
  and every file with `[ -f ]`; never invent.

- **P2-E CONFIRMATION — ZERO MOVED generalizes PAST operator files to the
  spatial-scheme file class, and the Haiku pre-inventory's MOVED column is
  ~100 % noise.** P2-E was `transport/spatial/{scheme,diamond,linear_discontinuous}.py`
  — NOT operators, different landing chapters (`foundations/discretization.rst`,
  `methods/sn/cartesian_multid.rst`, `methods/sn/loss_representation.rst`). The
  pre-inventory graded **13 MOVED**; re-adjudication overturned ALL 13 to
  TWIN/CONTRACT. Two "needs a new theory page" candidates (an advection–reaction
  interface, a reverse-mode transpose section) each already had a complete home
  — `§discretization-closures` even cross-references the exact code symbols
  (`outgoing_face_from_average`, `reaction_xs`). The ZERO-MOVED result is not
  operator-specific; it is a **Cardinal-Rule-3 consequence** (the theory shipped
  with the code), so budget ANY Phase-2 batch as TWIN-cut + CONTRACT-trim and
  grep the landing chapter before crediting a single MOVED. TWO sharp new
  discriminators the pre-inventory got wrong: (a) **a method that teaches
  d-generic / Kronecker / tensor-product structure is usually BOTH layers** —
  the LAYOUT theory is TWIN (→ `cartesian_multid.rst §spatial-moment-space`),
  but the reconstruction GOTCHAS (a d=1 trailing-axis-append; "keys on
  `face_moment_tail`, NOT a shape probe"; trace-order == inflow-order
  consistency) are CONTRACT. Same object, two layers — KEEP the contract, POINT
  the theory; never let "this is tensor-product teaching" auto-MOVE a method
  whose gotchas fail the wrong-simplification-guard test (overturned 3 LD
  methods: `_ubld_inflow`, `_ubld_outgoing_faces`, `moment_scan_closure`). (b)
  **the bit-identity operation-order discipline is a single-source across
  DOCSTRINGS** (Cardinal Rule 2, not a TWIN-cut): the canonical "explicit left
  fold, do NOT regroup" statement lives at the one helper
  (`_cartesian_streaming_diagonal`); sibling kernels that REPEAT it get trimmed
  to a pointer AT the helper, and the ⚠ do-NOT-regroup imperative stays only at
  the single source. Finally: **a contraction routinely surfaces a
  Cardinal-Rule-1 stale claim** (here a Protocol `is_affine_scannable`
  description said LD "does not qualify", false since #158 Increment B) — L-001
  applies MID-TRIM: verify the claim against LIVE code, FIX + report the
  scope-expansion, never transcribe the stale text into the trimmed form.

How to apply: for a Phase-2 operator-file batch, grep the landing chapters
for every teaching concept FIRST (expect all-TWIN); classify each block by
the "wrong-simplification guard" test for CONTRACT; keep latent-trap
imperatives + open-issue keep-decisions inline; cut TWIN/HISTORY to
greppable literal-path `§`-pointers (eq-label OK); prove zero code change
by token-invariance; skip the sphinx gate iff the file isn't automodule'd.
The lossless map is per-block `file:line | verdict | destination | id`,
written INCREMENTALLY, ending with verdict counts + before/after prose
lines + the 3–5 hardest calls that calibrate the siblings.

---

## L-034 — Code-prose rebalance, CONTRACT-DENSE file classes (#231 Phase 2, batches B + C + D + F + G): a machinery PACKAGE, a DRIVER file, an ABC/algebra-DEFINITION file, a CONTRACT-heavy OPERATOR file, a MESH/phase-space file, and a CURVILINEAR ψ½-operator pair are all contract-dense, so the honest cut is far smaller than the teaching-operator pilot (and that is CORRECT) — but the cut SURFACE differs by class

L-033 calibrated the rebalance on an OPERATOR file (`scattering.py`,
docstring −36 %). Batch B (`sn/loss_representation/{__init__,sweep_graph}.py`,
the SN sweep machinery, docstring −2.6 %) and batch C (`sn/solver.py` +
`numerics/iteration.py`, the SN driver + iteration primitives, docstring
−5 % / comment −16 %) both established that a **contract-dense file is a
DIFFERENT file-class from a teaching file** and the same rubric yields a
much smaller, honest cut. The file-class discriminator is the load-bearing
lesson — and each contract-dense class has its OWN dominant cut-surface:

- **Operator-SURFACE file** (in `transport/operators/`, `sn/operators/`): its
  prose TEACHES the operator algebra, which is 100 % TWIN in the
  operator-algebra book → aggressive TWIN-cutting (−30-40 %). **NUANCE (P2-G,
  batch G): the −30-40 % is for a TEACHING-heavy operator file (`scattering.py`
  was 73 % prose teaching the kernel). A CONTRACT-heavy operator file cuts FAR
  less** — `streaming.py` (docstring −16 %) and `boundary.py` (−4 %) carry the
  apply / solve / adjoint / reflect / split CAPABILITY contracts + the
  `_require_typed_composite` guard + the `_reflect_trace` adjoint-spelling
  ⚠-trap (a het-VACUUM-sphere-only catch) + the `reflect_rows_inplace`
  additive-not-overwrite ⚠-trap. The discriminator WITHIN the operator class is
  teaching-density vs contract-density; on the contract-heavy end the cut is
  campaign-TAGS-on-live-contracts (trim the tag, keep the contract), not
  standalone teaching. Latent-trap imperatives get an explicit `⚠` marker + the
  falsifying detail inline (L-033b), derivation → `§`-pointer. **CURVILINEAR
  sharpening (P2-F, `pole_angular_closure.py` docstring −23 % +
  `radial_characteristic.py` −15 %):** in the ERR-026/ERR-053 subsystem, KEEP
  the MATH FORMULA at point of use even when the teaching AROUND it is TWIN —
  the α-recursion index convention (`c_out=α/τ`, `c_in=(1−τ)/τ·α+α`), the
  `faces[g,m,i]` off-by-one, the τ_raw split formula, the seed-extrapolation
  `t` are each a file-local dependency a modifier depends on (a sign/index slip
  IS the historical hazard); cut the derivation, keep the formula + `§`-pointer.
  Two LATENT-TRAP imperatives were the load-bearing keeps — "do NOT tidy the τ
  arithmetic into a call to `contamination.morel_montry_weights`" (collapses the
  Leg-1 cross-check into a reference-contamination tautology, vv L11) and "do
  NOT re-implement the march at a call site" (single-source orchestration) —
  both fail the "would a file-local modifier do the wrong thing without this
  line?" test, so they stay inline even though their EXPLANATION is TWIN.
- **Sweep-MACHINERY / package file** (the `(L+C)` traversal realization, the
  DAG, the cell kernels): its prose STATES the local contract that *references*
  a book concept — "returns the FULL loss `(L+C)ψ`, Resolution A" is NOT
  teaching Resolution A (that lives at `§loss-rep-resolution-a`), it is THIS
  method's return contract. The "would a modifier who never leaves the file do
  the wrong thing without this line?" test keeps the vast majority. ZERO MOVED
  still holds (grep-confirm the book carries every concept), but the TWIN-cut
  surface is only the **module-head essays + campaign-relocation provenance +
  duplicated measured numbers**, not the method bodies.
- **DRIVER file** (the solver orchestration `sn/solver.py`, the iteration
  primitives `numerics/iteration.py`): its docstrings are the estimator /
  convergence / threat-model CONTRACT a modifier needs (kept near-whole —
  the `compute_keff` R7 role split, #291 leakage, balance identity,
  scale-bridge, the #282 lag-death certificate, the ERR-053 restart trap are
  ALL wrong-simplification guards). The dominant cut surface is therefore
  **COMMENTS, not docstrings**: standalone `#`-comment RETIREMENT TOMBSTONES
  (`_GaussSeidelResolvent`/`_MomentWindowedResolvent`/`_make_sweep_preconditioner`/
  the P1.7 `_build_rhs_*` block — git owns them; the live design is on the
  surviving function) + campaign-STATUS blocks (a 23-line "Issue #168 status"
  Phase-A/B/C/D narration annotating a 1-line default) + the HISTORY TAILS of
  SPLIT method docstrings (a "Scope"/"Verified"/"History:" section narrating a
  landed campaign under a CONTRACT algorithm). Batch C: comment −16 % (−101 ln)
  DWARFED docstring −5 % (−61 ln). Hunt the `#`-comment tombstones FIRST on a
  driver file, not the docstrings.
- **ABC / algebra-DEFINITION file** (the base-class file the book is ABOUT —
  `numerics/operator.py`, the LinearOperator ABC + composers/adjoint/inverse):
  the LEANEST cut of all (P2-D: docstring −2.4 %, comment −2.8 %, −51 ln). Its
  docstrings STATE the binding laws (closure/composition/adjoint-swap/
  homomorphism), the raise-conditions, and the typing-rationale. The trap: a
  Haiku classifier proposes MASS **MOVED** for "closure law"/"composition law"/
  "role classification"/"three-layer surface" (it did — 28). On an ABC file
  those verdicts **INVERT to CONTRACT**: the law IS the in-file contract at
  point-of-use (`OperatorSum.is_invertible`'s "ONLY the LEADING term need be
  invertible" — cut it, the next modifier "fixes" it to require both). The book
  teaches the concept's *derivation* (TWIN); the docstring states the *law the
  class obeys* (CONTRACT — never MOVED, the book already carries the concept →
  ZERO MOVED still holds). Dominant cut surface: inline **campaign-step
  provenance** (citation-vs-narration rule below) + multi-clause HISTORY
  narration stories, NOT the laws. The BATCH-SPECIAL row-8 dual-A bridge was
  ALREADY-SATISFIED here too (verify + report, per L-034's special sub-rule);
  the rebalance-read surfaced + fixed a stale `:ref:`operator-algebra-adjoint``
  → `operator-adjoint` (per the staleness-audit sub-rule).
- **MESH / phase-space file** (P2-G: `sn/mesh/augmented_mesh.py`, the `SNMesh`
  construction + property surface): like the DRIVER class, the dominant cut
  surface is COMMENTS, not docstrings (comment −26 % / −79 ln DWARFED docstring
  −0.6 % / −4 ln). The property + classmethod docstrings ARE the mesh's public
  API contract (the `bc` face-inventory-IS-BC-inventory invariant, `is_1d`'s
  ny==1-phantom gotcha, `full_field_space`'s G-adjoint metric, `dag_walk`'s
  XOR-signature) — kept near-whole. The narration lives in the `_init_core`
  CONSTRUCTION-BODY comment cluster (a 56-line Phase-C/D angular-closure-flip
  story annotating a 2-line default; the Wave-D/E 6-site migration essays; the
  deprecated-accessor tombstones) → hunt those `#`-comment clusters FIRST,
  exactly as on a driver file. The CONSTRAINT still stays even inside a
  construction comment (the CLASS-not-instance closure-bind reason, the "mesh
  provides shape only" B.5.A rule, the how-to-add-a-BC-kind recipe) — cut only
  the wave-flip STORY around it. Also surfaced a CORRECTNESS fix here's sibling:
  `streaming.apply_transpose`'s summary said "Hilbert transpose" while its body
  + the sibling `boundary.py` correctly say **Euclidean** transpose (`.H` is the
  metric Hilbert adjoint) — fixed per L-010 (a rebalance-read surfaces stale
  V&V vocabulary, per the ABC bullet's staleness-audit sub-rule).

Sharp sub-rules (machinery + driver classes):

- **Provenance trimming = citation-vs-narration, applied UNIFORMLY (internal
  consistency).** On a CONTRACT-dense file the inline provenance tags ARE the
  cut. Draw the "constraint stays / narration cuts" line at
  *citation-vs-narration*: TRIM landed campaign-STEP codes (`Wave O`/`carve
  PN`/`taxonomy §NN step N`/`spec §NN`/`Phase 2.5x`/`né _as_dense`/`O.2b`/`W-A
  collapse` — git + the theory page own them) but KEEP bare `#NNN` issue
  anchors (rubric: live-issue one-liners, the more durable traceability) and
  NAMED PATTERNS with theory anchors (`Design C`, `coding-elegance Pattern 2`).
  A bare `#280` is a citation (keep); "carried as documented twins until the
  3rd sibling fired the extraction trigger" is narration (cut). Apply to EVERY
  tag of a class or NONE — a half-stripped file violates internal consistency.
  A retired-SYSTEM lineage note ("the RUNTIME successor to the `CAP_SOLVE` tag")
  is HISTORY even when terse — `CAP_*` is fully retired, no live code references
  it, so it is pure archaeology (aggressive-retirement); the predicate law
  stands alone.

- **A hand-transposed-adjoint / reverse-scan / boundary-block comment body is
  the ALGEBRA-OF-RECORD — KEEP even though it reads like narration.** The
  cotangent routing, the seed-fold transpose, the degenerate-diagonal adjoint,
  the O.4b active-trace block, the moment-frame involution (an ERR-061-class
  diffusion-limit root cause): a modifier editing the adjoint/kernel MUST have
  these. Cutting them is the Cardinal-Rule-1 hazard the brief flags as "3
  constraint-bearing blocks misgraded HISTORY", at scale. A Haiku-style
  pre-classifier marks most of these COMMENT-cut [low-conf]; re-adjudicate EACH
  — nearly all are CONTRACT. Trimming a shape-annotation to save 1-2 lines
  strips a real invariant (the "DD/Step byte-identical" note is a *testable
  negative control*, not chatter).
- **Duplicated measured performance numbers → single-source to the canonical
  theory anchor, but keep ONE inline at the point-of-decision.** A perf basis
  repeated across N class docstrings (e.g. a Fork-B2 0.57-0.84× sweep basis) is
  TWIN with its theory home; point the *descriptive* docstrings there, but LEAVE
  the headline number in the *factory/selector* docstring where the choice is
  made. DRY (Cardinal Rule 2) without stripping numerical evidence from the code.
- **automodule + `:noindex:` makes the Sphinx gate LIVE (divergence from the
  L-033 pilot).** The pilot's `scattering.py` was NOT automodule'd → no build
  gate. A package rendered by `automodule:: … :members: :undoc-members:
  :noindex:` (here `discrete_ordinates.rst`) still RENDERS the docstrings (only
  the xref *targets* are suppressed — L-002), so a malformed docstring breaks
  `-E -W`. RUN the build gate both sides for an automodule'd file; grep-confirm
  0 `.. math:: :label:` / 0 `vv-status` first (cutting then orphans no
  `:eq:`/`verifies` target). Pointer nuance: the ratified literal
  `docs/theory/…rst §<label>` form is brief-correct, but in an automodule'd file
  a `:ref:`<label>`` role would render as a working link — flag it, don't
  unilaterally switch. **P2-F confirmed this on `radial_characteristic.py`
  (automodule'd, 0 warnings both sides) vs `pole_angular_closure.py` (NOT
  automodule'd → invisible, its 3 module-docstring `:label:` blocks cut-safe
  after the grep). TWO positive moves a RENDERED file affords that a
  non-rendered one cannot:** (1) **promote a LATENT-TRAP to a rendered
  `.. warning::` admonition** (the fission-double-apply HAZARD → a visible
  warning box, L-010 prophylactic-warning — better than an inline comment); mind
  the 3-space content indent under the directive. (2) **a section-RENAME during
  the cut stales an in-file back-ref** — renaming module-docstring headings
  ("Scope of this realization"→"Realized surface", "The ONE solve
  orchestration"→"Single source") staled two class-docstring back-references;
  grep the file for the OLD heading text after ANY rename and repoint (the
  rebalance-read staleness sub-rule, applied to in-file section refs).
- **A "fix if MISSING/drifted" BATCH-SPECIAL that turns out ALREADY-SATISFIED
  → verify against the oracle + REPORT satisfied; do NOT touch the correct
  CONTRACT block.** Batch C's SPECIAL 1 (the `notation.rst` row-8 dual-A
  bridge must survive in `iteration.py`'s module head) read as an instruction
  to edit — but the module head ALREADY stated it verbatim (posing +
  A=invertible-resolvent-operand + SN binding `A=L+C` gains `(S,B)`→`L+C−S−B`
  + fission-never-a-gain). The disciplined move: READ the oracle (`notation.rst`
  row 8), READ both ends, confirm the match, REPORT "SATISFIED, no fix" — a
  correct CONTRACT block a special protects is a KEEP, not an edit target.
  Same for a posing-drift special (SPECIAL 2, dated `(A−S−F)` posings): grep
  BOTH files, find zero, report CONFORMANT. A special is a *verification*
  obligation first, an edit obligation only on failure.
- **A rebalance READ surfaces a Cardinal-Rule-1 staleness bug — FIX it in-pass,
  flag it.** Trimming a comment means READING it, which is the only gate that
  catches a stale RAW PATH in prose (`-W` is blind to path strings — L-002; a
  `:class:`/`:func:` renders plain-text, a raw `orpheus/sn/scattering.py`
  string warns nowhere). Batch C: the source-helper comment cited
  `orpheus/sn/scattering.py`; the class lives at
  `transport/operators/scattering.py` (grep the live import). Fixed to the
  class ref (Cardinal Rule 1 supreme), folded into the tag-trim, REPORTED as a
  discrepancy. The rebalance is a free staleness-audit of every comment you open.

How to apply: FIRST classify the file — operator-surface vs machinery/package
vs driver (the folder + "does the prose teach the algebra, state a local
contract that uses it, or orchestrate/interface?"). Contract-dense (machinery
OR driver) ⟹ budget a small, surgical cut and KEEP the method bodies +
constraint comments; the cut surface = machinery's module-head essays +
provenance + duplicated numbers, driver's standalone `#`-comment tombstones +
status blocks + SPLIT-docstring HISTORY tails (hunt comments FIRST on a driver).
Run the Sphinx gate iff automodule'd (a `:noindex:` automodule STILL renders →
gate LIVE). Verify every batch-special as a check first (edit only on failure),
fix any stale raw-path you open, and REPORT the small-cut-is-correct finding
with the file-class rationale so the reviewer doesn't read −2-5 % as timidity.
Cross-links [[lessons-L33]] (the operator-file twin).

---

## L-035 — Orphan-slice adjudication (V7 backfill): the WIRE-vs-SENTINEL discriminator + the conceptual-root / foundation-coexistence corollaries, and the FAST theory-scan self-check

Adjudicating a batch of orphan eq-labels (RST `:label:`s with zero
`verifies` + no sentinel) into WIRE / SENTINEL / GAP has ONE sharp
discriminator, and it is NOT "is it definitional?" (almost everything on
a foundations page is):

- **WIRE** iff an existing test's PRIMARY assertion IS this exact equation
  against a STRUCTURALLY-INDEPENDENT reference — "would a sign/factor flip
  in the equation red this test?" YES. (`inflow_mask == flatnonzero(mu<-eps)`;
  `condensed SigT == fractional flux-weighted hand-sum`; `assert_balanced`
  on the collapsed mixture = the balance-PRESERVATION claim; `[K] ==
  np.linalg.solve(A,F)`.) Spelling MUST copy the `:label:` verbatim
  (a typo'd verifies = a matrix-flagged phantom).
- **SENTINEL** iff one of THREE structural shapes, NOT merely "it reads
  definitional": (a) a GENERAL SCHEMA / CONTINUOUS definition / LITERATURE
  identity whose CONCRETE / DISCRETE / TERMINAL instance is tested under a
  *different* label (general `M R = c_V I` → concrete `pi-r-equals-4pi-i`;
  continuous `Γ_±={Ω·n≷0}` → discrete `inflow-mask-discrete`); (b) a
  NATIVE-vs-LEGACY **bit-identity** regression (`axis_widths == legacy`,
  L-004 representational) — distinct from an independent-reference predicate
  test, which is WIRE; (c) documents code that does not exist yet
  (adjoint-weighted homogenization, blocked on an open issue).
- **GAP** only for a load-bearing COMPUTED contract with NO test anywhere.
  In a mature tree a whole 38-label slice can legitimately be 8 WIRE / 30
  SENTINEL / **0 GAP** — every "gap" turned out either tested (WIRE) or a
  definition/schema/literature identity verified downstream (SENTINEL). Do
  not manufacture a GAP to look thorough.

**Conceptual-root corollary.** A ROOT narrative page (e.g. `path_integral`,
"one object, five methods") states equations that the METHOD pages realize
and verify. ALL its orphans are SENTINEL (harness case-a: "a derivation step
whose terminal result is tested downstream") — EVEN when a formula IS tested,
because it is tested under the METHOD page's OWN label
(`path-integral-transport-correction`'s `D=1/(3Σtr)` is verified via the
diffusion page's `diffusion-coefficient`; `path-integral-generation-series`'s
`k=ρ(A⁻¹F)` via `matrix-eigenvalue`). Wiring the method-page test to the
root-page label is redundant double-labeling — SENTINEL with a rationale that
NAMES the downstream gate so a reader knows it IS tested, just elsewhere.

**Foundation-coexistence corollary.** A `:label:` backfilled in a late
label-pass (#231 Phase G) often has its test in a module-`@pytest.mark.foundation`
file whose docstring still says "software invariant — no theory :label:;
foundation carries NO verifies". That premise is STALE (the label now exists).
Resolve per-test: WIRE the ONE class that pins a COMPUTED physics formula (the
production-weighted `χ_mix` teeth-test `TestChiMixHandReference`), SENTINEL the
pure software-INVARIANT (the `Σχ=1` simplex law — the canonical foundation
case, vv-principles "foundation NEVER carries verifies"). Module-foundation +
class/method-`verifies` COEXIST and produce a real edge — the audit's
`_equation_coverage` reads `m.equations` regardless of the level tag
(`test_mixture_condense.py` is the in-tree precedent: module-foundation, class
`verifies("energy-condensation-rate-preservation")`). Add ONLY the decorator;
don't rewrite the stale docstring (scope-creep).

**FLAGGED-line-range can be stale — the doc's OWN named catcher wins.** A
stage-plan "wire to test at file:266-300" pointer had drifted: 266-300 was a
SPECTRAL cross-engine test that is Mode-12-BLIND to the mutation class (k∞
moves by *exactly* 0 under factor-swap/transpose — similarity + `eig(Mᵀ)=eig(M)`).
The RST prose itself NAMED the correct catcher (`test_K_operator_as_matrix_is_the_resolvent`,
the intrinsic `[K]==solve(A,F)` OBJECT gate). Trust the doc's named gate over
the brief's line number; verify it pins the OBJECT, not the spectrum.

**FAST self-check when N sibling batches edit concurrently.** Do NOT run the
full `python -m tests._harness.audit` (slow pytest collection, and its
tree-wide theory scan trips on sibling batches' in-progress sentinels). Call
`tests._harness.audit._scan_theory_equations(Path('docs/theory'))` DIRECTLY —
it validates every sentinel (same-file rule, spelling, `documented` set) and
resolves wired labels WITHOUT collection, in <1 s. Filter `scan.violations`
to YOUR file set; assert your new sentinels ∈ `scan.documented` and your wired
labels ∈ `scan.all_labels` \ `scan.documented`. Sentinel placement is
indent-agnostic to the parser (`line.strip()` first), but MATCH the enclosing
block's indent for RST rendering (3-space inside a list item / `.. warning::`,
2-space inside a bullet). The `.. (vv-status rationale)` prefix is parser-safe
(the regex needs `vv-status:` immediately after `.. `, and `(vv-status` fails
it). Template-B retitle: COPY the proven underline from the model page
(`collision_probability.rst`'s "Verification — what pins this chapter" =
37 `=`), never re-count by hand.

How to apply: classify each orphan by the 3-shape SENTINEL test / the
independent-reference WIRE test; treat a root narrative page as all-SENTINEL;
resolve foundation-file labels per-test (computed formula → WIRE, invariant →
SENTINEL); trust the doc's named gate over a stale line-range; self-check with
the direct theory-scan, not the full audit. Cross-links [[lessons-L03]]
(phantom/verifies-target hygiene), [[lessons-L04]] (representational →
documented), [[lessons-L10]] (Mode-12 spectral-invisibility vocabulary).

---

## L-036 — GROWING a thin "honest-stub" chapter to full at campaign close: flip the stale status, PRESERVE the landed-earlier section, RECONCILE sibling taxonomies, and deferred-wire the verifies-targets

The A6/ch15 shape (a campaign's closing docs phase): the chapter already
EXISTS as a deliberately-thin honest stub (`methods/sn/adjoint.rst` was
"deliberately thin today ... two layers landed, the third in flight") and
is ALREADY in the part toctree. The task is GROWTH, not authoring-from-
scratch, and it has a fixed anatomy:

- **The primary staleness is the stub's own "in flight / not yet landed"
  status** — flip it (L-007 tense-flip) the moment the campaign's phases
  (here A4/A5) merged. Verify the merge against git/the campaign plan's
  STATUS log, never the stub's frozen prose. The "three layers, two
  landed" framing becomes "three layers, all landed".
- **PRESERVE the already-landed section verbatim.** The thin chapter's one
  substantive section (here the #276 P3 `sn-scattering-adjoint` record)
  stays byte-for-byte — its `:label:`s are live #309 wiring-backlog items
  with exact vv-status rationales; touching them re-opens an adjudicated
  question. GROW AROUND it: new sections slot before (physics + route +
  taxonomy) and after (mechanics + carrier + verification + consumers +
  Development history), the preserved section sitting where it belongs in
  the new flow (S^T is a concrete instance of the "Euclidean transpose"
  category, so it lands right after the three-transposes taxonomy).
- **TAXONOMY RECONCILIATION is the sharpest new move.** When the charter
  asks for a NEW canonical "named landmine" section (the three transposes:
  Euclidean / Hilbert / continuous) that SUBSUMES ≥2 pre-existing sibling
  framings (loss_representation's warning contrasts {walk-Euclidean,
  μ-reversal, continuous}; the thin Key Facts named {Euclidean, Hilbert,
  walk-orientation}), do NOT contradict them — write the chapter as the
  authoritative RECONCILIATION with explicit subset relations:
  walk-orientation ⊂ Euclidean (the streaming realisation of Aᵀ),
  μ-reversal = the continuous adjoint's discrete SIGNATURE, Hilbert rides
  ON TOP of Euclidean via G. A reader who meets any sibling framing
  elsewhere lands on the same taxonomy. State the reconciliation in prose
  ("all three framings are the same taxonomy") so no future reader reads a
  contradiction.
- **Deferred-wire the verifies-target eq-labels you mint** (L-004 #3, the
  concurrent-owner case, made concrete): the certification/entries tests
  carry a comment "verifies un-linked until A6/ch15 mints the
  daggered-eigenproblem label" — they are WAITING for your `:label:`. Mint
  it UN-sentineled (a solver claim with a real L1 gate is a genuine gate —
  NEVER sentinel to paper over the transient orphan), and report DEFERRED
  WIRING with exact `test node → label` node-ids for the main agent (who
  owns the test files this phase). The `-E -W` build passes regardless
  (the audit is a SEPARATE gate the main agent runs AFTER wiring); flag
  the audit dependency loudly ("these N labels are orphans until wired").
  Definitional/literature siblings you mint alongside (the continuous
  adjoint equation) DO get `.. vv-status: <label> documented` — audit-clean
  at your build. Net: for a chapter minting a mix, some labels are
  documented-sentinels (clean now) and 1–2 are un-sentineled deferred-wires
  (clean after the main agent's marker commit).

**Two teaching-doc CORRECTNESS catches this class surfaces** (Cardinal
Rule 1 over faithful transcription):

- **Don't fuse an eigenvalue term and a fixed-source term in one
  equation.** The continuous adjoint written with BOTH a `1/k`-scaled
  fission gain AND an external source `q*` is inconsistent (the `1/k`
  belongs to the eigenvalue problem where `q*=0`). Present the eigenvalue
  form (labeled), then the fixed-source form (`q*=Σ_d`, no fission) in
  prose. A code/brief that hands you "the adjoint equation" generically
  can hide this fusion.
- **A code docstring's operator NAME can be sloppy where the theory must
  be exact.** The `KEigenvalue(A.H, (S+B).H, F.H)` spelling in the
  eigenvalue.py seam + the `solve_sn_adjoint` docstring uses "A" for the
  FIRST arg — which is the RESOLVENT `(L+C)`, NOT the loss
  `A_loss = L+C−S−B` (reading it as A_loss double-subtracts). The code is
  `KEigenvalue(resolvent.H, gain.H, F.H)`; spell `(L+C).H` in the teaching
  doc, unambiguous, and show `A_loss† = (L+C).H − (S+B).H` formed inside.

**Mode-12 live-application EXACTNESS is the vv-curator load** (Directive
5; this campaign twice caught a wrong-WHY here, so the prose is the
corpus's quoted spelling). For the daggered adjoint: `k` is EXACTLY blind
to (i) the factor-ORDER / similarity family (`eig(Mᵀ)=eig(M)` — ALL
factors transposed is a similarity), (ii) ALL vector content, and (iii)
**the G-metric itself** (`G'⁻¹AᵀG'` is metric-similar to `Aᵀ` for ANY
invertible `G'`, so the metric is a free parameter no eigenvalue gate can
EVER see — the sphere vector row is its SOLE catcher, ERR-067 family). But
leaf-transpose **DROPS** (F†→F etc.) are **NOT** k-blind — transposing ONE
factor is not a pencil similarity, k MEASURABLY moves (F†=F: 1.488→0.153
on the 4G ∞ fixture). Get the "blind-to vs not-blind-to" boundary exactly;
the corpus page and the `vv-principles` Mode-12 text must carry the same
measured spelling. Pair every `k`-row prose claim with the vector/pairing
catcher (spectrum, bi-orthogonality, duality, sphere-residual).

**Xref reality for a solver-return-type chapter.** Only the module
`automodule`'d WITHOUT `:noindex:` links (here `numerics.eigenvalue` →
`power_iteration` links); a return type in a NOT-automodule'd module
(`sn.solution` → `AdjointSolution`/`Solution`/`SolutionBase`) and a
`:noindex:`-automodule'd module (`sn.solver` → `solve_sn_adjoint`) BOTH
render plain-text BY CONVENTION (L-002) — consistent with the pre-A6
chapter's own `ScatteringOperator` refs. NOT a defect and NOT to be
"fixed" by adding an automodule: `sn.solution` carries `.. math:: :label:`
docstrings (homogenize/condense derivations), so automodule'ing it trips
duplicate-label collisions. Spot-check by (a) `-E -W` build EXIT 0
(catches intra-doc `:ref:`, all `:eq:`, `:cite:`), (b) grep the built HTML
for RAW `:class:`...`` role markup leaking (means a broken role) — none
should leak; every role renders as `<a>` link OR plain `<code>`, and (c)
confirm the ONE indexable module's refs actually link (typo-catch).

---

## L-037 — FLIPPING a "documented-future seam" to LANDED across an existing rich page: the stale-status blast radius is the WHOLE page (grep it), the wrong-rule can hide in a literature-table cell, and the wired⟹no-sentinel flip is verified against the LIVE test

The sibling of L-036 (which GROWS a thin stub chapter). Here an
EXISTING, rich foundations page carries a `documented-future seam` for
campaign X that just landed; the task is (a) flip every stale-status
claim to landed, (b) grow the seam section into the full landed
taxonomy, (c) un-sentinel the now-wired labels. Three distinctive
disciplines, none of which the brief's named file/line list fully
scoped:

- **The stale-status blast radius is the WHOLE page, not the brief's
  named 2–3 sites.** A "campaign X landed" brief names the capstone
  status block + one table row + one passage — but the SAME
  "blocked / not built / pending / documented-future seam / lands with
  P6 / in flight" claim is SCATTERED across Key Facts, the chapter
  overview bullet, section-body prose, a second table, AND a
  "V&V evidence lands with P6" closing line. Grep the page for EVERY
  future-tense/blocked token about X
  (`blocked|not built|not yet|pending|in flight|future seam|documented
  theory only|lands with P6`) and flip each on correctness grounds
  (Cardinal Rule 1) — the brief's list is the FLOOR. Worked (#281 P6):
  the brief named 3 sites; the grep surfaced SEVEN. Watch the
  "one remaining not-built discipline is X" overview bullet especially
  — when X lands, that sentence must RE-POINT to the *actually*-still-
  unbuilt sibling (here the least-squares dense-cross-Gram frame), not
  just drop X; a naive delete leaves the "one remaining" count wrong.
- **A "fix the loose (φ→φ*) phrasing" sub-task has a blast radius too —
  and the wrong rule hides in a LITERATURE-TAXONOMY TABLE CELL.** The
  brief named 2 prose sites; grep for the concept surfaced a THIRD in a
  "canonical pairs" table row. The tell for the bare-φ* trap is not the
  test-basis cell alone — it is the test/trial PAIRING: a row that lists
  `test = φ*·1_R` against an INDICATOR `trial = 1_R` silently encodes
  the bare-adjoint rule `∫φ*Σ/∫φ*` (worth-nonzeroing), NOT the bilinear
  `∫φ*Σφ/∫φ*φ`. The correct cell, matching the landed code (a capture
  gate asserts the frame weight IS the pair), is the PRODUCT weight
  `(φ*⊙φ)·1_R` against the indicator trial. Discriminator: with an
  indicator trial, the weight must be the PRODUCT; only a φ-weighted
  trial makes a bare-φ* test correct. Fix it, flag the scope-expansion.
- **The wired⟹no-sentinel flip is verified against the LIVE test, not
  the brief's assertion.** The brief said "both labels NOW carry
  verifies() (C1 stacks both; C4 stacks both) — REMOVE the sentinels."
  Per L-036's deferred-wire case a brief can say "wired" when the test
  is actually a WAITING verifies-target, so READ the live test files
  FIRST (`grep -n 'verifies\|class Test' <file>`) to confirm the
  `@pytest.mark.verifies("<label>")` decorators are really present and
  stacked — here they were (C1 stacks both labels, C2/C4 stack them).
  Then remove BOTH the `.. vv-status: <label> documented` directive AND
  rewrite (do not delete) its `.. (vv-status rationale)` comment to a
  plain `.. (Wired P6, #281 — no vv-status sentinel.) …` note naming the
  gates — the note prevents a future auditor from "helpfully" re-adding
  a sentinel to a long-documented-only label whose neighbours still
  carry them. Self-check with the FAST theory-scan (L-035): the flipped
  labels must show `label_exists=True, documented=False` with 0
  file-local violations (the label left the documented set cleanly and
  is now a covered verifies-target).
- **Grow the taxonomy by INCLUDING the generated fragment, adding
  UNLABELED supporting math, and keeping the ONE preserved verifies-
  target label byte-identical.** The five per-channel collapse rules
  come from `.. include:: ../../_generated/<name>.inc.rst` (same
  `../../_generated/` depth from `docs/theory/foundations/` as from
  `docs/theory/verification/`), NEVER hand-transcribed (L-008). The T0
  keystone / T1b angular / T4 balance / T6 carrier equations I add as
  the narrative are SUPPORTING identities — leave them UNLABELED (no new
  orphan/sentinel obligation; `git diff | grep ':label:'` must show ZERO
  net label change). The single labelled equation the section owns
  (`sn-homogenization-adjoint-weighted`, a verifies-target) is re-emitted
  BYTE-IDENTICAL inside the rewrite so git matches it as unchanged
  context — never rename a verifies-target while rewriting its prose
  (L-003).

How to apply: for a "campaign landed, modernize the page" task, (1)
grep the WHOLE page for stale-status tokens and flip each; (2) grep the
wrong-rule concept (not just the brief's 2 sites) — it hides in tables;
(3) verify verifies() markers in the LIVE test before un-sentineling,
then fast-theory-scan; (4) include the generated fragment, add
supporting math UNLABELED, keep the verifies-target label byte-identical.
Cross-links [[lessons-L36]] (the stub-growth sibling), [[lessons-L35]]
(the fast theory-scan self-check + WIRE/SENTINEL discriminator),
[[lessons-L03]] (never rename a verifies-target), [[lessons-L08]]
(generated artefacts are included, never hand-edited).

---

## L-038 — Auditing a "is the terminal docs phase done?" charter: a multi-phase campaign's LAST docs phase is often already-executed incrementally by the earlier phases' doc passes — verify by the page's OWN self-identification + build + cross-ref gate, don't infer a gap from an open plan line

A read-only "how much of Phase-N docs is already satisfied?" audit (here
the frame-projection campaign's P7 charter) has a recurring answer:
**effectively-done**, because each earlier phase (P3/P5/P6) ran an
archivist doc-pass that landed its own slice INTO the eventual capstone
page, so the "final docs phase" was executed piecewise before it was
formally reached. The audit discipline:

- **The plan's phase-line is a STALE tracking artifact, not the ground
  truth** (process-discipline: trust git/the shipped page, not a frozen
  plan claim). The driver plan still read "NEXT = P4.5 … P7 pending" while
  P4.5–P6 had all landed and the P7 page was written. NEVER infer "P7 is
  a gap" from an open plan bullet — read the shipped page.
- **A campaign's capstone page usually SELF-IDENTIFIES.** The decisive
  evidence was the page's own front-matter note titled "What shipped since
  (P3 / P5 / P7)" stating *"This page (P7) is the capstone…"*. Grep the
  candidate page's intro/Key-Facts/notes for the phase tag (`P7`, the
  issue #) — the campaign often already declared the page done in prose.
- **The plan's task-number ≠ the GitHub issue number.** The plan said
  "Tasks #46–#52 (P1–P7)"; #46–#52 are actually unrelated
  Thermal-Hydraulics/Kinetics issues — the plan used INTERNAL task
  numbering that COLLIDES with real issue numbers. The real trackers were
  the phase issues (#268/#226/#281/#275). Resolve "which issue tracks
  this?" by reading each candidate's title, never by trusting a plan's
  bare `#N`. A terminal-docs phase frequently has NO dedicated issue —
  its deliverables ride the phase issues' doc-passes.
- **Per-item verdict method for a content charter:** for each named
  deliverable, (1) locate the anchor/label (`grep -rn "<label>"
  docs/`), (2) READ the section (not the heading) and judge it against
  the articulation standard (does it carry the rejected alternatives, the
  structural WHY, the honest-scope seam?), (3) confirm the `-E -W` build
  is clean (charter's "-W clean" clause), (4) grep-gate the cross-doc
  `:ref:`/`:eq:` targets the section uses (the -W-BLIND plain-text class,
  L-002) — a Feynman-grade section with a dangling cross-doc ref is not
  actually done.
- **Distinguish a documented SEAM from a GAP.** A charter is DONE even
  when the page carries "stays a documented seam until consumer X exists"
  notes (here: the anisotropic-order Σ_{s,ℓ} moment-resolved pairing; the
  LeastSquaresFrame/GEC-rank>0 #275). A correctly-declared future-consumer
  seam (L-002 forward-ref discipline: literal not premature `:class:`) is
  the OPPOSITE of a gap — it is the honest-scope boundary the charter
  never asked to cross. Do NOT list a documented seam as owed work.
- **A charter's literal "the condensation PAGE" can be correctly
  delivered as a SECTION of a shared page** (DRY, Cardinal Rule 2): one
  frame page with space + energy as sibling PG sections beats two pages
  that duplicate the PG machinery. Flag the wording-vs-form deviation as
  INTENTIONAL, not a missing page — recommending a standalone page would
  MINT a twin-path violation.

Net verdict shape for such an audit: "effectively-done; residuals are
bookkeeping (mark the plan line ✅, no dedicated issue needed for
already-shipped work) — NOT a gap-fill pass." Cross-links [[lessons-L37]]
(the flip-a-seam-to-landed sibling), [[lessons-L02]] (the -W-blind
cross-ref grep-gate), and AGENT.md process-discipline (trust git, not the
frozen plan claim).

---

## L-039 — AUTHORING a campaign-CAPSTONE theory page (a completed feature's WHOLE story) from an algebra-of-record + plan memos + the error catalog: the narrative arc is motivation→derivation-of-record→design→discoveries→evidence→scope, and the vv-status decision for algebra-of-record SymPy-identity labels is verifies-COVERED (peierls foundation+verifies), NOT documented

Distinct from L-025 (a NEW shared-INVARIANT foundations chapter: gather
method-specific → generalize) and L-013/L-018 (a resolution/capstone
chapter of an ARC on an EXISTING page): here the whole task is a NEW
standalone page telling a COMPLETED campaign's whole story (consistent
DSA #2 — Fourier motivation, the four-step derivation, the design
decisions, the discoveries, the measured evidence, the honest scope).
The source-reading order and the label-status decision are the
load-bearing lessons.

- **Source-reading order for a capstone (extends L-005):** (1) the
  ROADMAP/plan-of-record (the phase structure, the RULINGS with dates —
  R4/R5/R6, the deviations); (2) the LITERATURE MEMO (the equations with
  paper-numbers + the errata/normalization watch-items — here Alcouffe
  (17)/(23) sign errata, the Σw=1-vs-2 map); (3) the ALGEBRA OF RECORD
  (`derivations/discrete/sn/dsa.py` — the SymPy `derive_*` functions ARE
  the equations; READ it, don't transcribe the memo's paraphrase); (4)
  the PRODUCTION code (the shipped shape — admission guards, the trace
  arm, the foldable accessors); (5) the ERROR CATALOG (the discoveries'
  full stories — ERR-070/071 parts); (6) the EVIDENCE PACK (the
  authoritative measured tables — SELECT the load-bearing ones). The
  memo is the NAVIGATION layer; the SymPy module + production code are
  the CORRECTNESS spine (algebra-of-record skill). VERIFY every DD-member
  collapse against the LIVE production body (`cell_update` returns the
  edge average; `moment1_update` the (28b) form) — the page's worked
  forms are code-grounded, not memo-transcribed.
- **The vv-status decision for an algebra-of-record page is the sharp
  new call (extends L-004 to the AUTHORING case).** When the brief says
  "mint `:label:`s on the key equations, then add `verifies()` markers",
  the derivation-identity labels (Larsen (27), (23a–f), (28), Marshak,
  the (33) synthesis) are algebra-of-record SymPy-identity gates — WIRE
  them foundation+verifies → **covered**, the peierls precedent (L-004
  case b: `test_case_method_*` carry BOTH `@foundation` AND
  `@verifies`), NOT `.. vv-status: documented`. The audit ACCEPTS
  foundation+verifies (confirmed: `test_dsa_rules` foundation gates
  verifies-cover their labels, 0 orphans). Reserve `documented` for the
  PURE-LITERATURE / STRUCTURAL labels with NO tight test (the ρ_SI=c
  motivating collapse, the 0.2247c continuum bound, the M=(I+𝒞)∘(L+C)⁻¹
  composition identity). Discriminator: does a test genuinely PIN this
  exact equation? derivation identity / object law / rate bound with a
  gate → verifies-covered; motivating/definitional literature with no
  gate → documented. **Mechanics:** the audit's `testable_labels =
  theory_labels − documented_labels`, so a `documented` label is
  EXCLUDED from the orphan gate and a `verifies`-covered label is NOT an
  orphan — either avoids orphan regression, but documented+verifies
  TOGETHER is muddy (a documented label with a test edge); prefer the
  clean split. If you wrote `documented` on a label you then decide to
  wire, REMOVE the directive (keep the rationale as a plain `..` comment
  naming the catcher) before adding `verifies()`, and re-run the audit
  to confirm 0-orphan.
- **The capstone's Key-Facts + narrative sections earn 13 labels; the
  page-prefix (`sn-dsa-*`) + a grep-collision check before writing is
  mandatory (L-003).** Author the head/intro as pure Write/Edit literals
  (no Python f-string over math — the L-026 brace trap); the f-string
  mangle grep (`A\^-1|G\^-1`) confirmed clean because I never routed math
  through a Python string.
- **Migrate a gate-file docstring's MEASURED-FACTS record to the page,
  keep its CONTRACT (enforcer NOTE f1; extends L-028/L-033).** A rate-tier
  test module's docstring carrying a "Measured design facts" section
  (D11/S2/ladder numbers) MOVES to the theory page's evidence tier; shrink
  the docstring to a greppable pointer (`docs/.../acceleration.rst`
  §`sn-dsa-rate-and-stability`) that keeps only the test CONTRACT (the
  #215-catcher mutation-matrix statement). A SIBLING test docstring that
  is already a pure contract statement + ERR narrative (the sweep-inverse
  gate) STAYS — the enforcer note says "keep its contract statement", and
  the ERR-071 story's canonical homes are the page + `error_catalog.md`.
- **The brief's paraphrased numbers/paths/targets are STARTING HEURISTICS
  — verify against the evidence + the tree (L-001 across the board):**
  (a) the brief's "Krylov 2554→21" was SI+DSA reflective in the evidence
  pack (Part A) — used the authoritative table, flagged the paraphrase;
  (b) the brief's `docs/theory/methods/diffusion/diffusion_1d.rst` path
  had no `diffusion/` subdir (real: `methods/diffusion_1d.rst`); (c) the
  brief's "point the drifted xref at field_algebra" — the drifted target
  was `operator_algebra` (only 1 passing DSA mention), and the NEW
  capstone page is the authoritative DSA home, so I pointed it THERE
  (better target, flagged the deviation); (d) the brief's "flip the
  field_algebra as_dsa_source promise" was ALREADY the landed-truth
  (L-001 already-fixed — reported, no action). FLAG each deviation in the
  return; the brief is the floor, the live tree is the rule.

---

## L-040 — RETIRING a per-X flag from the docs: the blast radius includes the TABLE COLUMNS that paraphrase it without naming it, and the flag's own JUSTIFICATION prose is usually independently FALSE

The symbol-grep (L-002) is the FLOOR of a retirement's doc radius, not the
ceiling. A `ClassVar`/field is documented in two registers: by NAME (which
greps) and by CONCEPT (which does not). Both must die.

- **Grep the symbol AND its human-readable paraphrase.** A brief listing
  N literal hits is complete only for the NAME register. Worked
  (`BoundaryTraceLaw.creates_sweep_cycle`, 2026-07-30): the brief's 7
  literal hits were exact — and MISSED 17 more cells, because the
  foundations page tabulated the flag under the header **"Sweep-cycle
  flag"** (6 per-law `True`/`False` values) and the SN page's resolution
  table carried 10 value cells under a header that DID grep. The
  paraphrase grep (`grep -rni "sweep.cycle"`) is what found them. RULE:
  after the symbol grep, grep the CONCEPT the symbol names (hyphen/space
  variants, the column header you'd write for it) — a `list-table`
  column is a documentation surface with no symbol in it.
- **Dropping a table column is a 3-part edit**: header cell, every row's
  value cell, AND the `:widths:` list (which must still match the new
  column count). Verify in the RENDERED HTML, not the source: parse the
  built page for `<col>` count + `<th class="head">` list + per-`<tr>`
  `<td>` counts. A widths/column mismatch is a real `-W` warning, but a
  silently-wrong-but-consistent table is not.
- **A column can often be REPLACED rather than deleted** — with the true,
  intrinsic property the false one was gesturing at. The `Sweep-cycle
  flag` column became `Trace-edge family` (none-the-inflow-is-data /
  same-face back-edge / opposite-face pair), read off the replacement's
  algebra of record. This keeps the pedagogical slot, makes the claim
  true, and teaches the exact distinction that killed the flag (the law
  owns its EDGE STRUCTURE; the configuration owns CYCLE-NESS). Prefer
  replace-with-the-true-invariant over delete-and-leave-a-hole.
- **The paragraph that JUSTIFIED the flag is the highest-risk prose on
  the page — re-verify it against the replacement, don't just delete the
  flag's name from it.** A per-X boolean's doc always ships a "and here
  is why X is False for these kinds" closing paragraph, and that
  paragraph inherits the flag's wrongness. Here: "Vacuum, white, albedo,
  and prescribed-inflow are all cycle-free" — FALSE for white, since
  `white|white` is cyclic for the same reason `reflective|reflective` is
  (the replacement gate's OWN docstring says so: "white on BOTH faces is
  not [acyclic]"). Read the replacement module + its test docstrings and
  rewrite the justification to the structural truth (only laws that
  supply the inflow as DATA contribute no edge at all, hence are
  unconditionally cycle-free). Cardinal Rule 1 outranks minimal-diff.
- **Option (b) — keep a retired-note section — over option (a) delete —
  whenever the retirement carries a DESIGN LESSON, and the L-007
  retitle-to-the-concept/KEEP-the-anchor move is what makes it free.**
  Retitling `The ``creates_sweep_cycle`` signal` → `Sweep cycles: a
  configuration property, not a per-law flag` kept the `bc-sweep-cycle`
  anchor, so both live `:ref:`s (one CROSS-DOC, which would have
  silently rendered plain-text if the anchor died — L-002) kept
  resolving AND auto-picked up the new title as their link text. Verify
  that payoff in the HTML (`grep 'href=".*#anchor"'` → confirm it is an
  `<a>` carrying the NEW title).
- **Structure the retired-note as an increasing-importance ladder, and
  put the un-fixable reason LAST.** Three findings: (1) zero production
  readers, (2) the attached claim was false, (3) it could not have
  worked in principle. Only (3) generalizes — so it gets the space, the
  measured truth table, and a named design rule in an `.. admonition::`
  ("a law may carry only what is intrinsic to it"). (1) and (2) are
  facts; (3) is the archaeology that stops re-invention. Sharpen (3)
  with the "one value, two different facts" tell: the flag read `True`
  for reflective meaning "can take part in a loop others close" and
  `True` for periodic meaning "closes a loop alone" — one boolean
  carrying two structurally different claims IS the proof the property
  does not live on the type.
- **Name the replacement's V&V level from the gate's marks, not from
  instinct.** The SCC gate is `@pytest.mark.foundation` with NO
  `verifies(...)` (verifies ⊥ level) — so the prose says "software/
  structural invariant of a discrete construction, not an equation
  claim", and cites the mutation-teeth test by what it proves (dropping
  the boundary edge FALSELY certifies acyclic). Never upgrade a
  foundation gate to an L-level in prose to make the section sound
  better-verified.
- **Running the mandated `-E` build REGENERATES `docs/theory/verification/
  matrix.rst`** (L-008 `builder-inited` hook) — and on a dirty branch it
  will absorb rows from OTHER uncommitted work (here +126 foundation
  tests from a sibling campaign, plus the ERR count). That is a
  legitimate by-product, NOT your edit; never revert it (it is
  generated), and REPORT it explicitly so the committer knows what the
  extra modified file is and can choose to stage it.

---

## Quality self-assessment rubric (Directive 3)

Rate each output 1–5 and log the weakest dimension in the return:
Derivation depth · Cross-references · Numerical evidence · Failed
approaches · Code traceability · Derivation source (from `derivations/`,
never hand-written). The recurring weak dimension on TERMINOLOGY/ROUTING
fixes is "numerical evidence" — structurally absent (no flux moves → no
convergence table), not a deficit; say so rather than manufacturing one.

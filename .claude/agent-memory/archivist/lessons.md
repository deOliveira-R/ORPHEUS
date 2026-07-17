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
  triple (`.resolvent=(L+C)`, `.gains=(S,B_a)`) — say so; the coupled M−N
  grid is the CARRYING-mesh case only.
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

## Quality self-assessment rubric (Directive 3)

Rate each output 1–5 and log the weakest dimension in the return:
Derivation depth · Cross-references · Numerical evidence · Failed
approaches · Code traceability · Derivation source (from `derivations/`,
never hand-written). The recurring weak dimension on TERMINOLOGY/ROUTING
fixes is "numerical evidence" — structurally absent (no flux moves → no
convergence table), not a deficit; say so rather than manufacturing one.

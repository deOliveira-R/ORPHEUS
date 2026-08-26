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

## L-041 — The DOC-ONLY "retire the false promises" pass on a subsystem under structural restoration: keep the declaration, make the CLAIM true; and prove doc-only by AST, not by eyeball

A B0 "clean before extending" phase hands you N measured docstring/prose
claims naming a consumer, capability or behaviour that does not exist. The
job is NOT deletion — it is making each claim TRUE while preserving what a
later phase needs. Disciplines that held across 18 items in 8 files
(boundary machinery, 2026-07-30):

- **Retire the CLAIM, keep the DECLARATION, and name the phase that will
  fill it.** When a brief says "do NOT delete these properties — phase B1
  populates them", the honest rewrite is three-part: state the measured
  present ("**currently unpopulated** — returns ``None`` on every law, read
  by nothing"), state HOW production reaches the same information today
  (the realizers recover `G` from the law's CLASS and `R` from
  ``law.albedo``), and name the landing ("**B1** mints the typed spec and
  populates this; the declaration is kept for that landing"). A reader then
  cannot mistake the property for live machinery *or* delete it as dead.
- **MEASURE the override lattice; never trust a "where applicable" hedge.**
  A doc saying "each subclass overrides the five universal invariants where
  applicable" is unfalsifiable prose over a lattice you can compute in ten
  lines: `ast.parse` each base body and count non-docstring statements;
  compare `getattr(Law, m).__func__ is not getattr(ABC, m).__func__` per
  (law × method). Here it produced the load-bearing numbers — 4 of 5 bases
  EMPTY, 2 of 5 overridden by NOBODY, 4 of 7 laws overriding nothing — and
  turned one hedged sentence into a per-row "*Intended* … / *Implemented*
  …" table. Also grep `raise <TypedError>` per row: a "Pinned error" column
  can name an exception **never raised anywhere in production** (ERR-040
  here), which is a second, independent hollowness.
- **A false-promise item's real scope is the CLAIM'S grep, not the brief's
  file:line.** Two of the strongest finds came from the closing gate, not
  the list: the same "downstream consumers (sensitivity adjoints) require
  the outflow trace" justification lived in a SECOND file
  (`sn/boundary/realizer.py`), and a self-contradiction the brief located
  in a docstring had a TWIN in a `#` comment 400 lines away (`B_b` …
  "present-zero bulk and trace" vs the class docstring's "the composite has
  no such slots"). Build the gate as one grep per item, run it AFTER the
  edits, and treat every extra hit as in-scope.
- **When two symbols name the SAME concept, the fix is a note that TYPES
  them, plus the collision's own tell.** A package writing `R = G_refl · α`
  (R = composite) beside an ABC writing `γ₋ψ = R G γ₊ψ + q` (R = response
  factor) will corrupt the phase that mints `R` as a type. Resolve by
  splitting every bullet into `G = …` / `R = …` and adding ONE `.. note::`
  stating the convention and naming what the composite is called instead
  (`R ∘ G`, never `R`). Then grep the concept: sibling modules using a
  DIFFERENT decomposition (`R_white = G_diff ⊗ α`, `R = Σ_α c_α G_α`) are a
  separate framing — FLAG them, don't sweep them into a scoped fix.
- **An error-message string inside `raise` is an EXECUTABLE statement —
  report it, don't edit it, under a doc-only constraint.** `apply_transpose`
  routing through a shared `_apply_faces` whose message says "``.apply``" is
  a real defect, but tests `pytest.raises(match=...)` on those strings.
  Report with file:line + the exact string + the reason it was withheld.
- **Prove doc-only with an AST gate.** Strip docstrings from every
  `Module`/`ClassDef`/`FunctionDef` body and compare `ast.dump(HEAD)` vs
  `ast.dump(now)` per file — "AST IDENTICAL" is the claim, a diff hunk is
  the exception. On a file the USER also edited this session, `ast.unparse`
  the stripped trees and `difflib` them: the printed hunk should be exactly
  their change (here one `kind` property) and nothing of yours. This is
  stronger and faster than reading a 450-line diff.
- **A brief's enumerated item can be a LATER phase's acceptance-gate text —
  leave it and say so.** The `SNBoundaryOperator` docstring's "``.H`` — the
  one channel by which the white-BC adjoint becomes available" is FALSE
  today (measured: `B.H` raises with a white face) but is quoted verbatim
  as phase B5's gate. Editing it would destroy the gate. Not on the item
  list ⟹ report as a deliberate non-edit with the owning phase.
- **Verify every issue number you cite** (`gh issue view N --json state`) —
  a docs pass that redirects a retired claim to "#183 tracks this" is
  minting a new claim, and a closed/mis-numbered issue is the same class of
  defect you are removing.
- **Baseline drift is real: this repo's `-E -W` baseline is now ZERO**
  WARNING/ERROR/CRITICAL (AGENT.md still says 1 — the `Mesh1D.from_geometry`
  `:paramref:` ERROR is gone). Re-measure the baseline every session; never
  assume the recorded number.

---

## L-042 — Auditing a corpus against a just-landed multi-commit REFACTOR: the phase-lag, the letter collision, and the retraction that INVERTS

The successor task to L-041's doc-only B0 pass: three commits land
(factor re-assignment · a new primitive · a domain narrowing) and you
must find and fix every `docs/theory/**` claim they falsified. 46 claims
adjudicated across 6 files; disciplines that were NOT already in L-041:

- **A brief's lead can be directionally inverted — settle it with
  `git show <fix>^:<path>`, not the ⚠ alone.** The lead said "the page
  says `apply_transpose` writes the `outflow_indices_for_face` slots;
  the code's ⚠ names that the WRONG spelling." Reading the ⚠ *precisely*:
  the wrong spelling is scattering over the law's own **codomain**
  (Γ₋); post-change the transpose genuinely DOES land on the Γ₊ rows,
  so the sentence is incidentally true NOW. The doc was still wrong —
  because it described an output **projection** that the PRE-change code
  never performed either (it wrote the whole face). Only the pre-commit
  body settled which of three readings was right. A ⚠ names a HAZARD,
  not what the code did; when a lead claims "the doc documents a known
  bug as the contract", read BOTH bodies.
- **A phase-N doc pass leaves phase-(N−1)'s falsifications behind — audit
  the PARAGRAPH FAMILY, not the commit's diff.** Phase B1 populated the
  two factor properties on all seven laws and never touched the theory
  page; B3.0's doc pass fixed only what B3.0 moved. Result: a correctly
  re-typed G/R section sitting three screens from "`geometry_map` and
  `response_kernel` return the ABC's `None` on **every** law and are read
  by nothing". The reader cannot tell which phase staled what, so a
  scoped-to-phase-N audit ships a self-contradicting page. Fix it and
  FLAG the scope expansion. (The `-E -W` build is blind to all of it.)
- **Replace an unfalsifiable inventory sentence with a MEASURED table**
  (the L-041 override-lattice move, now for property VALUES): one
  `python -c` over the seven laws printing `law.geometry_map` /
  `law.response_kernel` turned the false sentence into a 7-row
  ground-truth table the next reader can re-run.
- **One letter, two decompositions, two pages — type them both in ONE
  `.. warning::`.** The rank-N expansion `B = Σ G_α ⊗ A_α` (a sum over
  TERMS) collides with the affine factorisation `R G` (a factorisation of
  ONE term). Tell: a Marshak formula `R = c₁G_refl + c₂G_diff` that uses
  `R` for the whole composite AND files a Lambertian average under `G`.
  Fix = rename the colliding prose symbol (`R`→`B`), add ONE warning
  stating both decompositions and naming the composite **`R ∘ G`, never
  `R`**, then correct the mis-tiered rows. Extends L-041's same-module
  collision to the cross-page case.
- **A retraction can INVERT a claim, not just kill it — give each item its
  own `**Disposition:**`.** Three published "why this matters"
  consequences: #1 *measured not to exist* (the phantom future consumer —
  the declared-capability-no-consumer pattern), #2 **inverted** (the
  argument was "the realization is a self-adjoint idempotent projector";
  once the domain narrows, the operator is not an endomorphism, so
  idempotence is not even a well-typed question and the type tag the
  argument rejected is now the right one), #3 *right observation, wrong
  layer* (the uniformity was real; the mechanism was one layer too
  shallow). Per-item dispositions preserve the intellectual content that
  a blanket tombstone destroys — and #2's inversion is the single most
  instructive line on the page.
- **A "the gate still does X" claim is verified against the TEST BODY,
  and you must COUNT the rows.** A snapshot suite's vacuum case had been
  re-posed in the body while its class docstring still described the old
  semantics (so quoting the docstring would have re-minted the falsehood
  — L-001 in a test file). And only **3 of 7** cases were re-posed: the
  un-narrowed four still feed the full face. My first draft said "every
  case was re-posed" — a fresh falsehood, caught by reading all seven
  bodies. Also: the mixed-BC row is an `xfail(strict=True)`; document it
  as an **honest red that flips on the next phase**, never a suppression.
- **Put the Mode-12 blindness IN the prose, beside every table of
  realized operators.** Measured `|Γ₊| = |Γ₋|` on every quadrature × face
  ⇒ a shape assertion cannot distinguish `Γ₊→Γ₊` from `Γ₊→Γ₋`. The next
  reader's instinct is to "check" a typing claim by output shape; one
  sentence ("read the *declared spaces*, never the output shape") is
  worth more than the table it annotates. Same family: a three-way
  partition (inflow ⊔ outflow ⊔ **tangential**) needs "**not inflow** is
  NOT **outflow**" restated wherever a complement could be spelled.
- **`python -c` every numeric constant a doc asserts.** A page said the
  tangential band is "default `ε = 1e-12`"; measured
  `TANGENTIAL_EPS = 4·np.finfo(float64).eps ≈ 8.9e-16` — four orders out,
  never warned, pre-existing.
- **A "mitigation" sentence can be a NON-SEQUITUR wearing a hedge.**
  "The adapters carry every tangential ordinate at μ = 0 **so** ψ = 0
  there for a properly-initialised flux": μ≈0 is the *definition* of
  tangential, and it does not make ψ vanish. The honest weaker statement
  (no operator writes tangential slots, so a zero-initialised carrier
  keeps them zero) is an INITIALISATION property; the genuinely
  structural neighbour (tangential rows carry zero *metric weight*) is a
  different fact. Split them.
- **The taxonomy the brief asks you to ADD may already exist in a
  docstring the concurrent agent just wrote.** Grep the label before
  authoring. Here `.. _bc-method-realizability:` was already DEFINED in a
  package `__init__` docstring — inert only because that package is not
  `automodule`'d. The project pattern is **page owns the label, docstring
  `:ref:`s it**; use the same label name on the page so the two agree,
  and report the latent duplicate rather than inventing a second name.

---

## L-043 — The brief's "MEASURED evidence, do not re-derive" block is a CLAIM, not data: reproduce it, because a wrong attribution there propagates into the corpus as fact

L-001 says verify the brief's *facts*; L-043 is its numeric face, and it
is sharper because the brief actively tells you **not** to check. A
doc-repair brief that hands you a measured-evidence block ("bit-identical
on A and B, 1 ULP on C and D — do not re-derive") is handing you the
single most quotable content in the whole task: numbers go into a table,
and a table is what the next session cites. If the attribution is wrong,
you have laundered an error into the algebra of record with your own
credibility on it.

- **"Do not re-derive" means "don't burn a session on it", not "don't
  check".** Cost-scope it: reconstructing a retired operator body from
  the `git diff` you already read is ~10 lines. Do that. Reserve the
  deference for things that genuinely need a solver run.
- **Worked (B3.4a):** brief said white was *bit-identical* on
  `gauss_legendre(8)` AND `product(2,4)`, 1 ULP on `lebedev(17)` /
  `level_symmetric(6)`. Reproduction: `product(2,4)` is NOT bit-identical
  — and it is precisely the quadrature where the retired `> 0.0`
  classifier disagrees with `TANGENTIAL_EPS`, so it CANNOT be
  (the two operators are not computing the same functional there).
  Publishing the brief's line would have told a future reader "the
  change was FP-neutral everywhere", **exactly wrong on the one
  quadrature that motivated the phase**.
- **The deeper find, and the reason to reproduce at all: two mechanically
  DIFFERENT effects can measure the same.** Reduction-order drift
  (padding removed → numpy's pairwise tree reassociates) and a genuine
  mis-classified-row VALUE bug both read ≤ 1 ULP on an `O(1)` probe —
  because the offending weight is itself `O(ε)` (measured `cos_w` =
  7.85e-17 against a norm of 2.5651, so Δnorm = **0.0 exactly** and the
  entire discrepancy is ψ-weighted in the NUMERATOR). So the error scales
  with the **flux ratio**, is unbounded by floating point, and reaches
  6.1e-05 at a 1e12 ratio. A ULP table therefore CANNOT justify such a
  change; the justification is structural ("one classifier, not two").
  Say that in the doc — and add a `.. warning::` that the 1-ULP row is
  not evidence of equivalence, or the next reader re-derives the wrong
  conclusion from your own table.
- **While the probe harness is open, AUDIT THE WHOLE INVENTORY — the
  brief's sample is never the population.** Sweeping all production
  quadratures × all six faces (~15 lines) turned a two-row anecdote into
  a scoping law: the disagreement occurs ONLY for the `product` family
  and ONLY on `xmax`/`xmin`/`ymax` (`ymin` has the same tangential count
  and zero mis-admissions — the sign flip moves the round-off across
  zero). Two publishable consequences fall out that no one had stated: a
  **tangential-count audit is not a sufficient screen** (`lebedev` has 12
  per face and mis-admits none; `level_symmetric` has none at all), and
  the exposure is **face-asymmetric within one quadrature**, so a fixture
  exercising one face can be green while its opposite is wrong.
- **The audit routinely falsifies a claim in a file you may not edit.**
  Here it killed "that is every production quadrature but
  `gauss_legendre`" in a `.py` docstring outside the permitted set. FLAG
  it with the measured replacement sentence verbatim, so the fix is a
  paste for whoever owns the file — a flag without the corrected text is
  a ticket, not a hand-off.

How to apply: treat every number the brief marks MEASURED as
unverified-until-cheap-repro; reproduce, widen to the full inventory,
and publish only what you measured — attributing anything you could not
reproduce to its source rather than adopting it.

---

## L-044 — A retirement's doc blast radius in UN-autodoc'd modules is invisible to `-W` AND to `-n`; and the same retirement silently DEMOTES the rewired tests' claim class

Two findings from one sweep (retirement of `orpheus/sn/quadrature.py`:
4 per-family adapter classes + an `AngularQuadrature` Protocol → named
classmethod factories on the single `Quadrature`).

**(a) Nitpicky mode is NOT the missing gate.** L-002 says a dead
Python-domain role renders plain-text with no `-W` warning and suggests
`-n` as the stronger probe. MEASURED here: `-n` saw **zero** of the 22
dead `orpheus.sn.quadrature` refs — because Sphinx only nitpicks what it
RENDERS, and none of the carrying modules (`numerics/measure.py`,
`numerics/operator.py`, `numerics/quadrature/*`, `geometry/reduced_operator.py`)
is `automodule`'d, nor is any `tests/**` file ever read. So for a
docstring in an un-autodoc'd module the ONLY gate is `grep`. Before
concluding "`-n` would have caught this", check `grep -c "docstring of
<module>" <nitpick.log>` — 0 means the module is invisible at every
severity. Corollary for the fix: edits to such docstrings CANNOT move the
warning count, so "count unchanged" proves nothing about them — the
acceptance evidence is the grep inventory plus a per-hit
KEEP/FIX adjudication, and the build is only a no-regression control.
(The `-n` counts were 6493 py-domain lines / 6964 total, byte-identical
pre and post by set-diff; the normal `-E -W` build is at **0** warnings
now — the `mesh.py` `:paramref:` baseline in AGENT.md is STALE.)

**(b) The sharper find — a retirement can DEMOTE a gate's claim class
without touching one line of the test body.** Four tests named
`test_*_bit_identical_to_legacy_adapter` had their comparison target
`sed`-migrated by the retirement commit (`LebedevSphere(order)` →
`Quadrature.lebedev(order)`), keeping the local name `legacy`. But the
factory *calls* the very rule function under test, so what was a
two-implementation bit-identity gate became a value compared with itself
routed through a wrapper — it can NEVER detect the node drift its
docstring still advertised. The test stayed green, the name stayed
authoritative, and the docstring kept the stronger claim. RULE: when a
retirement rewires a comparison target, re-ask **"are the two sides still
independently produced?"** — if the survivor is the caller of the other,
the gate has silently dropped to a pass-through check and every doc/
docstring crediting it must be re-scoped (name the real pin — here the
cylindrical regression snapshots — rather than deleting the gate's
description). The tell in a diff is a variable still called `legacy`
beside a brand-new API. The tree had already self-corrected exactly ONE
of the four (`test_rules_1d.py`, with the excellent phrasing "compares a
value with itself routed through a wrapper … the drift gate is the
separate test below") — reuse an in-tree honest framing when one exists
rather than inventing vocabulary.

**Adjudication shape that worked.** Every surviving mention became a
past-tense double-backtick LITERAL ("the four SN-side wrappers this
docstring used to point at (``…LebedevSphere``) were retired into
classmethod factories on the one ``Quadrature`` type"); every live claim
got a `:meth:` at the successor. Two mechanical gates make the sweep
auditable: `grep -rn "<retired path>" … | grep -v '``'` must be EMPTY (no
surviving hit outside a literal) and a role-regex over the retired names
must be EMPTY. Both beat reading the diff. Widen the grep from the
module PATH to the bare CLASS NAMES — that surfaced 3 more dead roles
the path-grep missed, and separates the genuinely-dead from the
live-homonym (here `name="LebedevSphere"` is a **live registry
`QuadratureSpec` identifier** asserted by ~30 tests — never "fix" it).
The residue — a retired class name used as an informal FAMILY LABEL in
prose ("ProductQuadrature(2x4)", "LS-family (``LevelSymmetricSN``…)"),
~25 sites — is a separate sweep: FLAG it, don't half-do it.

How to apply: for any retirement, run the grep inventory yourself and
adjudicate hit-by-hit (history KEEPS, present-tense-false FIXES); do not
let a clean `-n` build imply the docs are clean, and re-derive what each
rewired test can still SEE before repeating its old claim.

---

## L-045 — `tools/check_docstring_xrefs.py` is the gate L-044 said did not exist; and in `tests/` a dead xref is a TRIPWIRE for a false claim about what the gate proves

L-044 concluded that for an un-`automodule`'d docstring "the ONLY gate is
grep". That is now superseded by a real gate, committed this session:
**`tools/check_docstring_xrefs.py`** resolves every *fully-qualified*
Python-domain role by **importing** it, so render coverage is irrelevant.
`.venv/bin/python tools/check_docstring_xrefs.py tests --quiet` →
`DEAD TARGETS : 0` is a hard, cheap acceptance criterion; exit 1 gates CI.
It deliberately ships an EMPTY `ALLOWLIST` — never add to it. Two design
facts to respect: (i) it skips UNQUALIFIED refs (`:meth:`Quadrature.product``)
because Sphinx resolves those against module context and flagging them
manufactures false positives — so an unqualified dead ref still needs
grep; (ii) it separates "module exists but raises on import" (a TOOLING
problem) from "genuinely absent" (a dead ref) — opposite fixes.

**Why `tests/` is the sharpest surface.** It carried 495 fully-qualified
xrefs and **nothing had ever checked one of them**: Sphinx never reads
`tests/` at any severity, so `-W` and `-n` are both structurally blind
(L-044(a)). First run: **41 dead targets across 62 sites**.

**The load-bearing finding — a dead ref in a TEST docstring is a
tripwire, not a typo.** A test docstring states what the test PINS, so
the retirement that killed the ref usually also invalidated the
surrounding CLAIM. Rate measured here: 3 of 41 dead targets sat inside a
present-tense-FALSE claim, and one was a whole-file misdescription
(`test_unified_matvec_sphere.py` still advertised a six-step
unified-vs-legacy bit-identity chain; BOTH implementations had been
deleted and the surviving file holds two σ_t/zero sanity gates — the
accurate story was already written in a CLASS docstring three screens
down). Another (`test_native_matvec.py`) had a 7-item pin list whose
item 5 asserted the **inverse** of the live gate (docstring: "face
residual zero at non-outflow ordinates"; test name:
`test_outer_face_inflow_slots_carry_the_identity`), plus one retired and
one inverted item. RULE: on every dead ref in `tests/`, read the test
BODY, not just the sentence — then REPORT the false claim explicitly
rather than quietly repointing the link, because a wrong claim about
what a gate proves is worse than a dead link (Mode-11/Mode-12 adjacent:
the docs are the only place the claim is written down).

**The four-way adjudication that worked** (per site, never blanket):
REPOINT (symbol moved — the majority, 46/62 here) · PAST-TENSE LITERAL
(the sentence is history; flip the tense and demote the role to
double-backticks, since a role PROMISES a live link) · REWRITE (the
claim is present-tense-false) · DELETE (rare; here **zero** — every
sentence carried content). A useful sub-case: a not-yet-built module
(`orpheus.derivations.registry`, "gets promoted into … per the wave3
plan") is a LITERAL, and while there also verify the cited PLAN FILE
still exists — `.claude/plans/wave3/` did not.

**The brief's own successor map is a hypothesis.** It said
`orpheus.sn.angular_flux.AngularFlux` "now lives at
`orpheus/transport/fields/angular_flux.py`". `git log --diff-filter=D`
on the old path showed the opposite: the legacy class was a *different*
object (bulk + conflated boundary buffer + history) DELETED at
`d8843ba9`, replaced by the composite `TimedFullField` whose bulk is the
same-named L2 class. Same name, different object ⇒ the 5 sites split
into past-tense literals (history) and one REWRITE (a module docstring
naming the wrong return carrier while the gate asserts
`isinstance(state, TimedFullField)`). ALWAYS run the deletion-commit
read before trusting a "just moved" claim.

**Mechanics.** A module-path rename is a mechanical `str.replace` over
`tests/**` — but ORDER the mapping longest-first (`_quadrature_recipes`
before `_quadrature`) and first prove every occurrence is inside a
docstring/comment (a live import cannot contain a dead path, so
`grep -v` the role spellings and eyeball the remainder). 43 replacements
landed that way; the ~19 judgment sites were hand-edited.

**Proving doc-only.** The reviewer's check IS the author's check: parse
HEAD and worktree, strip docstrings, compare `ast.dump`. 42 files,
0 diffs. This also proves no `@pytest.mark.verifies/catches` moved —
those are V&V registry edges. Note what it does NOT cover: comments
(absent from the AST, and legitimately editable), and **f-string
assertion messages, which ARE code** — two stale `cells_view` mentions
had to be left in `_require(...)` messages and REPORTED instead (L-041).

**Free catches worth taking.** The tool's `DECIDABLE_ROOTS` excludes
`tests`, so intra-tree `:mod:`tests.…`` refs are ungated — a 20-line
ad-hoc resolver found 5 more dead ones (a partly-executed
`test_trajectory_resolvent_*` rename family). Raw path strings with LINE
NUMBERS (`orpheus/derivations/peierls_geometry.py:2906`, ~10 sites) are
the opposite call: repointing the directory without re-verifying the
line implies a verification you did not do — FLAG, don't half-fix.

How to apply: run the tool as the acceptance gate on any tests/ or
docstring sweep; treat each dead ref as a claim to re-verify against the
test body; adjudicate four-way; prove doc-only by AST; report every false
claim you find and never fix the underlying gate yourself.

---

## L-046 — An import-resolving xref gate OVER-reports: an annotation-only class attribute is live to Sphinx and dead to `getattr` — prove the anchor before degrading a ref to make a gate green

The `orpheus/` half of L-045's sweep: **30 dead targets / 37 sites**, all
residue of module MOVES the corpus-side edit had already fixed while the
docstrings were left behind (`8cda6b73` rewrote four theory pages for the
`TransportSolver` retirement; the four *modules* naming it kept their
refs for months).

**The finding that matters: 5 of the 30 were NOT dead.**
`resolve()` imports the longest importable prefix and `getattr`s the
remainder — so a class attribute that exists only as a **class-level
annotation** reports `missing` while being a perfectly live symbol:

| target | declaration | verdict |
|---|---|---|
| `Field.UNITS` | `UNITS: ClassVar[Unit]` + `#:` comment | **renders as a LIVE `<a href="#…Field.UNITS">`** |
| `DiscreteMeasure.nodes` | `nodes: np.ndarray` | dataclass field |
| `AngularSymmetry.continuous_isotropy` | `continuous_isotropy: SubgroupOfO3` | dataclass field |
| `AngularTraceSpace.omega_dot_n` | `omega_dot_n: NDArray = field(...)` | dataclass field |
| `WithinGroupSystem.loss` | `loss: "CoupledOperator"` | dataclass field |

The decisive evidence is a **rendered-anchor grep in a FRESH build**:
`grep -o 'id="orpheus.numerics.field.Field.UNITS"' <build>/api/numerics.html`
→ present, with three inbound `href="#…"`. So autodoc DOES emit
`py:attribute` for annotation-only members; the gate's probe cannot see
them. Turning that ref into a literal to make the gate green would
DELETE a working link — the exact inversion of the gate's purpose. Leave
it, report it with the anchor as proof, and hand over the resolver patch
(after the `getattr` `AttributeError`, accept
`attribute in getattr(obj, "__annotations__", {})`). **Never edit a gate
you were not asked to edit, and never mutilate a true ref to satisfy
one.**

**The mirror-image class — genuinely unresolvable, and worth fixing.**
`napoleon_use_ivar = True` (docs/conf.py) makes numpydoc
`Attributes`/`Parameters` entries render as `:ivar:`/`:param:` FIELDS,
which mint **no** `py:attribute` target. So an instance attribute
assigned in `__init__` (`self.bc = …`, `self.pole_angular_closure = …`,
`self.eigenvalue_method = …`) is unresolvable in EVERY build no matter
how well documented. Honest spelling: a live `:class:` role on the owner
+ the attribute as a literal — "the realized ``bc`` dict on
:class:`~…DiffusionMesh`". Three such sites here. Check
`napoleon_use_ivar` before assuming an `Attributes` section creates
targets.

**Adjudication mix (25 genuine).** 12 REPOINT · 4 past-tense literal ·
9 REWRITE. The REWRITE rate is much higher than a rename sweep would
suggest because a *move* leaves a true-but-relocated symbol, whereas a
**deletion** leaves a sentence whose PREMISE died: the paragraph that
justified the deleted thing inherits its wrongness. Three shapes recurred:
(i) a **self-contradicting file** — `geometry/boundary/__init__.py`
past-tensed the dissolved realizer registry at line 173 and still
present-tensed "three stub realizers self-register at import time" at
line 466; (ii) a **completed migration still written as future work** —
`reduced_operator.py` said consumers "will migrate in Wave G" while
`SNMesh.__init__` already calls its factories; (iii) a **docstring
contradicting its own body** — `BasisSpace.solve_critical` documented a
`d=None` fallback the code deletes with a `ValueError`.

**A retirement DEMOTED a gate, again (L-044's rule, live).** The
`reduced_operator` docstring credited
`tests/geometry/test_reduced_operator.py` with hash-equality "vs the
legacy SNMesh setup methods". Those methods are gone; the test now
compares `spherical_streaming(mesh, quad)` against `sn_mesh.reduced.*`
— the value that same factory produced, through the mesh constructor —
and the two remaining `SNMesh.face_areas`/`delta_A` legs are DEPRECATED
read-throughs to the same object. It pins the WIRING, not the math. The
fix is a `.. warning::` naming what survives (the SN curvilinear
snapshots + the τ producer-equivalence floor), not a deleted claim. Same
stale claim sits in `docs/theory/foundations/structured_geometry.rst`
§"Bit-identical contract" — whose very next section says the methods
"no longer exist".

**Also worth a look every time:** a `:mod:` naming an UNBUILT plan target
(`orpheus.transport.problems`) must be a literal (L-002/L-014); a type
name invented in prose and never minted (`FlatVectorLike` — grep says it
never existed in `orpheus/`) is a dead ref with no successor, so name the
actual construct (here the duck-typed *ravellable* `to_flat` /
`from_flat` pair, deliberately unnamed so `numerics` need not import
`transport`); and a **code block whose first import raises `ImportError`**
(`docs/theory/references/trajectory_resolvent.rst:3816` still imports
`GeometrySpec`) is the hardest MUST-FIX of the tense discriminator.

How to apply: run the tool, then adjudicate EVERY hit against the live
tree before touching it — a `getattr`-based gate over-reports on
annotation-only members and under-reports on unqualified refs (L-045).
Prove a contested hit with a rendered-anchor grep in a FRESH build.
Fix the mirror class (instance attributes) properly. Expect the deletion
residue to be a stale PARAGRAPH, not a stale token.

---

## L-047 — The `docs/` half of the xref sweep: a dead ref in an `api/` page is a retired API SURFACE (rewrite the section), and `:noindex:` makes a whole page's roles plain text so "would `-n` catch it?" is moot

Closing the trilogy after L-045 (`tests/`, 41 dead / 62 sites) and
L-046 (`orpheus/`, 30 / 37). `docs/`: **20 dead targets across 24
sites in 15 pages → 0**, `-E -W` EXIT 0 with the diagnostic set
unchanged from a freshly-measured **0**-warning baseline.

**The tree-specific shape: `api/` prose dies in BLOCKS, not tokens.**
7 of the 24 sites were one page's one section — `docs/api/geometry.rst`
"Factories", listing `Zone` / `mesh1d_from_zones` / `pwr_pin_equivalent`
/ `pwr_slab_half_cell` / `homogeneous_1d` / `slab_fuel_moderator`, all
retired in one Phase-F commit (`81b083be`). Six repoints would have been
six lies: the successors are not renames but a **re-layering** (geometry
`StructuredGeometry` + `Region`, then mesh `Mesh1D.from_geometry` +
`RegionMesh`), so a section whose thesis is "the factory layer is the
recommended construction path" is stale as a THESIS. The tell that
distinguishes this from a rename sweep: the retired names' own module
docstring already carried the successor map (`factories.py` opens with
"Phase F retired … The 1-D path is now …"). **Read the surviving
module's docstring before planning repoints** — a well-retired module
tells you whether you owe N edits or one rewrite. Ratio here: `docs/`
came out **7 REPOINT · 5 mirror-class · 4 past-tense literal · 8
REWRITE**, i.e. a third rewrites, matching L-046's finding that
deletions leave stale PARAGRAPHS.

**`:noindex:` suppresses the WHOLE page's anchors — measured.**
`docs/_build/html/api/method_of_characteristics.html` contains **zero**
`id="orpheus.*"` anchors; so does `api/discrete_ordinates.html`. Every
`automodule` on those pages carries `:noindex:`, which renders docstrings
but mints no `py:` targets — so *every* py-domain role there is plain
text whether or not the symbol exists, and a live `href` to
`#orpheus.sn.solver.SNSolver` sits in the tree pointing at an anchor that
was never created. Consequence for L-002/L-044's "would `-n` have caught
it?": on this corpus the question is doubly moot — nitpicky can only
nitpick what it RENDERS, and here rendering doesn't mint targets either.
The import-checker is the only gate. (Corollary when rewriting such a
page: adding an `automodule` is still worth it for the *docstrings*, but
do not expect it to make roles link.)

**The `napoleon_use_ivar` mirror class is not a curiosity — it was 5 of
24 sites** (`SNMesh.scheme`, `SNMesh.pole_angular_closure`,
`SNSolver.inner_solver`, `KEigenvalue.eigenvalue_method`, plus
`MOCQuadrature.n_azi_2` which never existed at all). An
`__init__`-assigned attribute is unresolvable in EVERY build; adding
autodoc coverage will NOT revive it. Standard rewrite, applied five
times: a live `:class:` on the owner + the attribute as a literal, phrased
so the sentence says where the value comes from — "the ``scheme``
attribute that :class:`SNMesh` realizes in its constructor", "the
``eigenvalue_method`` constructor selector on :class:`KEigenvalue`".
Contrast the OVER-report class (L-046): the checker now falls back to
`__annotations__` across the MRO, so dataclass fields and `ClassVar`s
resolve correctly and need no defensive verification.

**Three dead refs sat on a claim that a rename had INVERTED, not just
moved.** (a) `operator_inverse_family.rst` published a 3-line code block
for `_seeded_inverse` reading
`return cast(SupportsSeededApply, cast(SupportsInverse, A).inverse())`
with prose crediting "the CALLER's `A.is_invertible` guard". The live
`seeded_inverse` (public since #276 A4, `cd000c2e`) has **no cast at
all** — two `TypeGuard` bridges — and runs its **own** `invertible()`
guard raising `NotInvertible`. The published sentence named the wrong
mechanism AND the wrong guard owner; a repoint alone would have left both
falsehoods with a working link. (b) `api/method_of_characteristics.rst`
claimed "azimuthal angles are adjusted slightly from an even distribution
to satisfy the cyclic condition" — `MOCQuadrature.create` is a plain
midpoint-even `linspace`; cyclicity is reached from the **ray-spacing**
side (`effective_ts = (t_max−t_min)/n_rays` per angle). (c) A
`.. code-block:: python` introduced by "A user constructs one with:"
opened on `from orpheus.derivations.common.geometry_spec import
GeometrySpec` (module deleted at `81b083be`) *and* used `np` with no
`import numpy`. **RUN every doc code block that a present-tense sentence
promises works** — this one now does, verified end-to-end
(`k_eff = 1.0000000000000002`), as does the new geometry construction
block. A dead import is the loudest possible dead ref and no build sees it.

**Scope discipline against an owner issue.** Two sites sat on a page an
OPEN issue (#286) already owns. Correct move was neither "defer, it's
theirs" nor "fix everything": fix the dead refs, fix the **measured**
adjacent falsehoods (`ReducedStreamingOperator` has no `tau_mm` field;
2 of 8 claimed deprecated properties survive on `SNMesh`), leave the
issue's genuine mechanism-rewrite item alone, then comment with a
measurement table and the residue's **corrected path** (the issue named
`docs/theory/discrete_ordinates.rst`, which no longer exists — the
section is now `docs/theory/methods/sn/index.rst:472`). An owner issue
whose cited paths have rotted is worth less than the five minutes it
costs to re-point them.

**A retirement DEMOTED a gate, again — and a sibling agent proved it
independently.** `structured_geometry.rst` credited
`tests/geometry/test_reduced_operator.py` with bit-identity "vs the
legacy SNMesh setup methods"; those methods are gone and `SNMesh.__init__`
now calls the factories itself, so the surviving legs compare a fresh
factory call against the value that same factory produced. Fix: past-tense
the history, `.. warning::` the demotion, name the gates that DO carry the
math (`sphere_*`/`cyl_*` regression snapshots,
`test_tau_producer_equivalence.py`, `test_alpha_closed_form.py`). A
concurrent `tests/geometry/` pass was simultaneously renaming those legs
to `test_*_is_the_factory_value` with a docstring reading "`array_equal(x,
x)` for any face-area math whatsoever" — same verdict, reached
separately. **When you suspect a gate was demoted, say so in the doc; the
convergence is evidence, and the doc is the only place the claim is
written down.**

**THE SESSION'S SHARPEST CORRECTION — I named a gate that was BLIND to
the claim I credited it with.** Having correctly demoted the
bit-identity gate, I wrote "the mathematical content is pinned
elsewhere: … `test_tau_producer_equivalence.py` + `test_alpha_closed_form.py`
for the closure coefficients themselves." The τ half was FALSE, and it
was falsified by measurement, not argument: under fully-garbaged
`spherical_streaming`/`cylindrical_streaming` factories that file passes
**5 tests in 0.03 s**. Cause — the same #236 Step C move I had just
documented two screens earlier in my OWN τ-ownership note: τ left the
reduced operator for the angular closure, so the gate compares
`morel_montry_tau_per_level` against `morel_montry_weights`, both
derived from `(μ, w)` alone. It pins a quantity `reduced_operator.py`
no longer produces. **I had written the premise and still drew the
wrong conclusion from it**, because "τ is a closure coefficient" and
"the reduced operator carries the closure coefficients" were both true
sentences on the page and I composed them without checking the
referent.

The rule this earns, and it is the doc-prose analogue of
`vv-principles`' "a `catches` marker is a COVERAGE CLAIM, not a topic
tag": **a doc sentence of the form "gates X, Y pin claim C" IS a
coverage claim, and it must be justified the same way — by a mutation
that reddens X and Y, never by topical adjacency.** Corollaries:

* **Replacing a demoted gate is the moment of maximum risk.** Having
  just proven gate A blind, the reflex is to reach for the
  nearest-sounding sibling. That reflex is exactly what produced the
  error; the sibling inherits neither A's scope nor the claim's.
* **The correct citation is PER FIELD, not per topic.** The measured
  replacement is five different files, one per array — `delta_A` →
  the closed-form `test_delta_A_magnitude` (its **sole** catcher; the
  snapshots are legitimately blind because `delta_A` has no production
  consumer); `alpha_half` → `test_per_ordinate_flat_flux_consistency[SPHERICAL]`
  (`catches` ERR-006/007); `alpha_per_level` → `test_alpha_closed_form.py`,
  **cylindrical-α only** (every fixture is `CoordSystem.CYLINDRICAL`);
  `redist_dAw` → `test_streaming_equilibrium_curvilinear.py`'s L0
  `φ = Q/(Σ_t(1−c))` identity — and NOT the flat-flux identity, which
  recomputes `ΔA/w` instead of reading the production array;
  `face_areas` → `tests/geometry/test_geometry.py` on the producer
  `compute_areas_1d`. A single "these gates cover it" sentence cannot
  be true at that granularity.
* **The SAME gate cited for TWO different claims can be right once and
  wrong once — narrow the correction, do not sweep it.** The τ gate
  appears twice on this page: at the `morel-montry-clamp` vv-status
  rationale (**correct** — it does pin τ) and in the demotion warning
  (**false** — it does not pin the reduced-operator arrays). The
  coordinator had to send a second message narrowing the first,
  because a blanket "that gate is wrong" would have destroyed a true
  citation. When told "citation of X is false", ask *false for WHICH
  claim* and grep every occurrence before editing any.
* **When a sibling agent corrects a shared claim in code you may not
  edit, MIRROR its wording rather than re-deriving prose.** The
  corrected `reduced_operator.py` docstring already carried the
  measurement and the per-field table; re-verifying each catcher
  against the live tree (5 for 5) and mirroring cost minutes and
  guarantees code and corpus say the same thing.

**Free catches en route.** A `scipy` role can die by UPSTREAM removal —
`scipy.special.sph_harm` was deprecated in 1.15 and removed in 1.17 (the
tree runs 1.17.1); successor `sph_harm_y` has a **swapped `(n, m)`**
argument order, which belongs in the fixed sentence. A `:mod:` naming a
planned promotion (`orpheus.derivations.common.chord_oracle`) needs the
literal AND a check of the plan file it cites
(`.claude/plans/trajectory_resolvent_hindsight_refactor.md` — gone; the
"scheduled" framing had to become "open opportunity, trigger = a second
consumer"). And a rename ripples past the roles: `_select_si_resolvent →
_select_si_splitting` had 1 dead role plus **3 present-tense literal
mentions** on two other pages, invisible to the checker (literals aren't
roles) — after fixing the flagged role, grep the OLD NAME tree-wide and
adjudicate every literal by tense.

**Sanity note on the by-product.** The `-E` build regenerated
`docs/theory/verification/matrix.rst` (+2 foundation tests,
`test_docstring_xrefs`), because the committed matrix predated that test
file landing in HEAD. Legitimate; report it, never revert it (L-008).

How to apply: run the checker as the gate on `docs/` too; before planning
repoints read the surviving module's docstring to learn whether you owe N
edits or one section rewrite; expect ~⅓ REWRITE and ~⅕ mirror-class; RUN
every promised code block; grep the old name for literals after fixing the
roles; and when a page has an owner issue, fix + measure + comment with
corrected paths rather than defer or annex.


## L-048 — An equation in a doc has TYPES and a SCOPE, and NO gate checks either: a type-incoherent identity and a one-instance proof stated for a whole class both build clean forever

**The task.** Correct two defects in the boundary-condition *algebra of
record* (`docs/theory/foundations/boundary_conditions.rst`, the
`bc-factor-roles` section) as step 0 of campaign step G6.3. Both defects
were in *published mathematics* — not in a cross-reference, a symbol name,
or a status claim — so every gate this agent owns was green on them:
`-E -W` at 0 warnings, `check_docstring_xrefs` at DEAD TARGETS 0, and the
V&V matrix regenerating without complaint, for as long as the claims had
shipped.

### Defect class A — the identity does not TYPE-CHECK

The rank-one theorem concluded, verbatim, **`R ∘ G = R`**. With
`G : Γ₊ → Γ₋` and the classifying `R : Γ₋ → Γ₋`, the left side is
`Γ₊ → Γ₋` and the right side is `Γ₋ → Γ₋`: they cannot be equal as
operators. The step had silently identified `Gᵀv` with `v` by treating
`v = |Ω·n|` as a *function* without tracking which half-trace it is
restricted to — **harmless until campaign phase B3.2 narrowed the SN law
onto `γ±` and made the two halves genuinely distinct spaces**, a type
abuse from that commit onward.

Note the shape: nobody introduced an error. A *code* carve (B3.2)
retroactively falsified a *math* sentence three chapters away, and the
falsification is a **type mismatch inside an equation** — a defect class
with no gate at all. The correct statement (`R ∘ G = u ⊗ (Gᵀv)`, and
`Gᵀv = v` *as a function*, so the composite is `G`-INDEPENDENT rather than
equal to `R`) preserves the theorem's whole content and every downstream
conclusion, including "the B3 correction leaves the composite unchanged".

⇒ **Read an equation's DOMAINS AND CODOMAINS, not just its symbols.** For
every displayed identity on a page you touch, name the space each side
lives in and check they agree. When a page's spaces have recently been
NARROWED (a half-trace split, a domain restriction, a role-keyed retype),
that narrowing is a licence to re-type-check every identity mentioning the
affected spaces — the old ones were written when the spaces coincided.

Bonus catch from doing this properly: the slogan's real hypothesis is that
`v` is `G`-INVARIANT, not merely that `R` is rank-one (a rank-one `R` with
`v = δ_{Ω₀}` makes `G` fully observable). Keep the memorable slogan — it is
cited by `:ref:` from 5 code/test sites — and add ONE sentence naming the
hypothesis that actually does the work.

### Defect class B — proven for ONE instance, stated for the CLASS

A paragraph headed "**The crossing is geometric**" argued — correctly —
that the specular mirror is the unique ambient isometry fixing the face,
exchanges the hemispheres, and preserves `|Ω·n|`. It then closed *"which
is why `G` and not `R` carries it"* — for **every** law. False for a law
with no isometry: a wall is not a quotient (the page's own sufficient test
says so), so nothing provides the crossing geometrically and the
**physics** does it, by integrating the outflow and re-emitting an inflow.

⇒ **On any "which is why X" generalisation, ask: for WHICH instance is the
argument given, and is the conclusion stated for the CLASS?** The tell is a
paragraph that opens on one concrete object (a mirror, a slab, a 1-group
case) and closes on a universal. Scope the existing argument to the case it
proves — keep it verbatim in substance — and *add* the missing case. Do not
rewrite the proof; it was never wrong, only over-quantified.

Corollary that made the fix stronger: the honest law ("whichever factor is
non-trivial carries the crossing") has **boundary cases worth writing out**
rather than hiding. Rank-0 laws cross vacuously; a bare scalar response on
an ANGULAR trace is non-trivial in *magnitude* and trivial in *angular
structure*, so NEITHER factor crosses — which is exactly the
already-documented realizer REFUSAL. Boundary cases that turn out to be
shipped refusals are the best possible evidence that the reformulated law
is the right one.

### The framing that DISSOLVED both — taxonomy ≠ computational factorization

Both defects had one root: `R ∘ G` is a **TAXONOMY** (does this law's
content come from geometry or from physics? — decided by multiplicativity +
the quotient test) and the page presented it as also the **computational
factorization**. As a classification `G : Γ₊→Γ₋`, `R : Γ₋→Γ₋` is coherent;
as a realized typing it is false, and `[M]` **no realized response in the
tree is an endomorphism of `Γ₋`** — the Lambertian's realization types
itself `Γ₊ → Γ₋` in its own first line.

⇒ **When a doc carries a factored / classified form of an object, state
explicitly whether it is a classification or a recipe** — and check the
declaration tier against the REALIZATION before writing either. Two
successive campaign designs were built on this conflation, each refuted by
one read of the realizing class. A dedicated short section ("the taxonomy
and the factorization are different questions") is the right shape: it
gives the distinction a citable anchor so the next design cannot re-make
the assumption silently.

### What I got RIGHT and should repeat

- **Reproduced every number before publishing it.** Re-ran the design
  probe (`R* = G₊⁻¹RᵀG₋`, `G_S` over 11 orders): `0.0` / `1.110e-16` /
  `1.110e-16` / `7.628e-01` degenerate, weighted adjoint law at exactly
  `0.0`. Verified the claimed one-line transpose
  `Rᵀ(φ) = (cos_w/norm)·Σφ` against the DENSE transpose of the operator's
  own `apply` — `max err 0.0`. Verified `is_adjointable is False`.
  Re-measured the plan's `2/2, 4/4, 49/49` `|Γ₊| = |Γ₋|` claim and
  WIDENED it (6 quadratures × every face, incl. `level_symmetric(6)`
  24/24) instead of transcribing the sample.
- **Refused to cite an ephemeral path.** The probe lived under
  `$CLAUDE_JOB_DIR/tmp/`, which no future reader can open, and
  `scratch/g6_design_measurements.md` did NOT record it. Instead I
  DESCRIBED the construction (shapes, metrics, comparison) so the table is
  reproducible from the page alone. **A path that will not exist is a
  stale raw path the moment it is written** — describe, don't cite.
- **Reported rather than edited a concurrently-edited code file.**
  `orpheus/geometry/boundary/_factors.py` was being rewritten by the main
  agent DURING my session (its mtime moved mid-task; it had already fixed
  defect B and the taxonomy framing there, but NOT the `R ∘ G = R`
  identity). Two residual sites in it were REPORTED, not touched.
- **Kept the heavily-referenced anchor.** `bc-factor-roles` has 8 inbound
  `:ref:`s across docs + code + tests; the fix went INSIDE it and the two
  new ideas got their own new anchors, so zero referrers moved.

### Two mechanical facts worth keeping

- **`⭐` and `⛔` have ZERO occurrences anywhere in `docs/`.** They are the
  agent-memory / plan / code-docstring vocabulary. `⚠` (9 sites) and `✓`
  (43) ARE corpus vocabulary. I drafted with `⭐`/`⛔` carried over from
  the plan and had to strip them — **grep a glyph in `docs/` before
  importing a marker from a plan or a docstring.**
- **A Python `SyntaxWarning` appears in the Sphinx log and is MISSED by a
  case-sensitive `grep 'WARNING:|ERROR:|CRITICAL:'`.** A non-raw docstring
  containing `\Gamma` produced `SyntaxWarning: "\G" is an invalid escape
  sequence` mid-build; it does not bump the exit code either. Add
  `SyntaxWarning` to the build-log grep — and before reporting one in a
  file the main agent is editing, re-`py_compile` the LIVE file (mine was
  already fixed a minute later).

### vv-status for an equation whose CODE does not exist yet

The new adjoint identity got `.. vv-status: bc-response-factored-adjoint
documented` — correct, because the FACTORED spelling is not built (the
Lambertian ships as ONE operator with `is_adjointable=False`; the
factorization is a later campaign step), so no `verifies` marker has a
function to point at. The rationale comment names the **precondition**
gates that DO exist (`test_the_half_trace_metric_is_strictly_positive`
pins the non-degeneracy the identity requires;
`test_the_metric_is_not_euclidean` pins that the metric is load-bearing at
all) and ends with *"when the factorization lands, WIRE a test to this
label and REMOVE this sentinel"* — so the sentinel carries its own exit
condition. Documented-only count moved 524 → 525 on regen, as expected.

⚠ **UPDATE 2026-08-04 — that exit condition FIRED six days later** (G6.3
step 3b) and the sentinel's rationale became present-tense-false in the
exact words quoted above. See **L-049** for what to do when you are not
the owner of the generated artefact the un-sentineling would move: keep
the DIRECTIVE, rewrite the RATIONALE to state that the precondition
expired + why the directive is still there + the gate that now exists,
and flag it. A sentinel that carries its own exit condition still needs
somebody to notice the condition fired — nothing in the build does.

⚠ **The `-E` regen also absorbed a row from the main agent's uncommitted
work** (`numerics/test_angular_face_trace_space`, +56 foundation tests,
8478 → 8534 total, shifting every share percentage). Legitimate
by-product; report it, never revert it (L-008).

How to apply: type-check every displayed identity against the spaces its
sides live in, especially after a narrowing carve; scope every "which is
why" to the instance actually proven and add the missing case, writing out
the boundary cases; say whether a factored form is a classification or a
recipe, and check the declaration against the realization FIRST; reproduce
and WIDEN every measured number; describe a probe instead of citing an
ephemeral path; grep a glyph in `docs/` before importing it; add
`SyntaxWarning` to the build-log grep.

---

## L-049 — A class retirement's docs blast radius is THREE tense classes, not two; and a composite `[M]` cannot certify its factors

**Task (2026-08-04):** repoint the dead xrefs left by deleting the welded
`AngularAverageOperator` (`orpheus/sn/boundary/angular.py`), replaced by
`IsotropicEmissionOperator @ PartialCurrentOperator`. 19 mentions in
`docs/` source: 12 role-bearing (the checker's list) + 7 literals the
checker cannot see. Bounded, "do not restructure".

### The THIRD tense class — a falsified PREDICTION

The brief's discriminator was the standard two: present-tense claim ⟹
REPOINT; past-tense history ⟹ DE-ROLE to a literal. Both applied and both
were right. But 4 of the 12 role sites, and the richest prose on the page,
were neither — they were **future-tense predictions written while the
replacement was still a plan**, and the landing falsified them:

* *"The type that will host it exists —* `ScalarTraceSpace` *… the
  per-face accessor and the factored `B ∘ C` spelling are* **not built
  yet**." Both halves shipped — but the host is NOT `ScalarTraceSpace`. It
  is a new `Γ`-ladder tier (`AngularTraceSpace.current_space`, a unit
  metric, per-face) and the shipped `current_space` docstring goes out of
  its way to say `ScalarTraceSpace` is a *different object* (the `(J⁺,J⁻)`
  pair for the whole boundary under the face-AREA metric — hosting one in
  the other double-counts the area weight). "Not built yet" → "built"
  would have shipped a page asserting the wrong type.
* *"that is phase **B5**, which is what makes its adjoint structurally
  available"* — landed at **G6.3 step 3b**, and by *factoring* rather than
  by the predicted `u ⊗ v` typing, which **dissolved** the transpose
  ambiguity instead of resolving it.

⭐ **A deferral contract names THREE things — the MECHANISM, the
HOST/TYPE, and the PHASE — and a landing can falsify any subset
independently.** Check each separately against the shipped code; do not
let "the seam closed" license a blanket tense-flip. (Sharpens lessons §4's
"verify the SHAPE that shipped": the shape is three fields, not one.) The
honest repair preserves the prediction and tombstones it — I kept the
sentence and added a `.. note::` naming what was predicted, that it did
not hold, and *why the shipped host is a different object* — which is
strictly more informative than the correct sentence alone would have been.

### A composite `[M]` measurement cannot certify its own FACTORIZATION

The page carried, correctly, `Rᵀ(φ) = (cos w / norm)·Σφ` verified `[M]`
bit-exactly. It built that from `Cᵀ(s) = cos w · s / norm` and
`Bᵀ(φ) = Σφ` — and the SHIPPED split is the other way round:
`Cᵀ(s) = cos_w ⊗ s` (no `/norm`) and `Bᵀ(φ) = Σφ / norm`. Composite
identical, per-factor formulas both wrong. Re-measured live: every one of
the four identities reads `0.0` on `product(2,4)`/`xmax`, so **no
measurement on the page could ever have caught it.**

⭐ The structural reason matters and belongs in the prose: the
normalisation lives in `B` because `C` produces a **current** and `B` must
produce an **intensity**, so the division is the unit-changing step —
which is what leaves `S(f)` carrying an honest `J⁺` and lets an albedo
enter as the pure scalar law `J⁻ = α J⁺`. A doc that gets the split
backwards silently deletes that argument. **When a doc factors an operator,
verify EACH factor's formula against its own live `apply_transpose`, not
the composite** — and re-verify the design-probe description too (this
page's probe also built `C` with the `/norm`; one clause fixed it).

### The vv-status sentinel whose exit condition fired, when you don't own the artefact

`.. vv-status: bc-response-factored-adjoint documented` (minted at L-048)
was correct while no production function realized the identity. The
factorization landed and the label now has a `verifies` marker on the
SHIPPED chain — so the directive is removable. But removing it
re-categorises a **generated** artefact (`theory/verification/matrix.rst`)
in a session whose remit is dead-xref repair, with the code owner mid-carve.

⭐ Resolution: **keep the DIRECTIVE, rewrite the RATIONALE.** The comment
now opens `⚠ THE SENTINEL'S PRECONDITION EXPIRED AT G6.3 step 3b AND THE
DIRECTIVE BELOW IS NOW REMOVABLE — left in place only because
un-sentineling re-categorises a GENERATED artefact and is owed the
regeneration`, quotes the superseded rationale verbatim as history, and
names the exact gate (`tests/sn/operators/test_lambertian_chain.py::
TestReciprocityAgainstTheMirrorFace::test_H_is_pointwise_the_mirror_face_kernel`).
Zero false text, zero silent V&V-category change, next session has the
whole decision. Contrast the two failure modes it avoids: flipping the
category unasked (invisible to `-W`, changes a generated table) vs.
leaving the quoted-false rationale (a future reader re-derives a dead
precondition).

### Mechanics worth keeping

* **`git status` at task start named the branch, not the one in the
  session snapshot** (`refactor/operator-strategy-layers`, not `main`),
  and by task end the main agent had COMMITTED the code side underneath
  me (`b4f0f5c9`). Re-read git before every claim about tree state.
* **The baseline `-E -W` warning was 1 and it was TRANSIENT** —
  `verification matrix regeneration failed: pytest collection failed
  (exit 2)`, caused by the main agent's half-saved test edits. `pytest
  --collect-only -q` succeeded moments later. Before attributing a
  baseline warning to the corpus, re-run the underlying tool directly.
  Final build: 0 warnings, EXIT 0.
* **`orpheus.sn.boundary.angular` is not `automodule`'d anywhere**, so
  every role on the BC page — the dead ones I removed AND the live ones I
  added — renders as plain text with no `href`. Verified in built HTML
  (`xref py py-class` spans, zero `<a>`). Matching the page convention is
  correct; adding an `automodule` for the two leaves I touched would be
  half-surfacing. The refs are still right and become links the day the
  package is surfaced.
* **A verbatim historical ERROR MESSAGE quoted in past tense
  (`` ``AngularAverageOperator.apply: psi.shape[0] = |Γ₊|, expected N`` ``)
  is correct history and already a literal** — leave it entirely. That was
  1 of the 19 sites and the only one needing no action.
* **Inline literals wrapping across two source lines render as ONE
  `<code>`** — the pre-existing table cells already did this, so
  ``` ``(IsotropicEmissionOperator(...) @ PartialCurrentOperator(...))\n
  & IdentityOperator()`` ``` is safe; `-W` catches the `:widths:`/column
  mismatch if you get a row wrong.
* **Measure the operator tree, don't infer it.** The two "MEASURED" code
  comments claim a realized `repr` shape; I walked the live tree by
  `__slots__`/`__dict__` (`realize_recursively` → `OperatorSum` →
  `ScaledOperator` → `TensorProductOperator` → `OperatorProduct(B, C)` →
  `IdentityOperator`) before writing it, and it matched
  `orpheus/geometry/boundary/white.py`'s own spelling exactly. Use the
  production module's spelling when it has one — that is the SSOT.

**Residue flagged, not fixed (owner = main agent):** 7 present-tense
`orpheus/` mentions the import-based checker cannot see (unqualified role
at `sn/boundary/angular.py:271`; "the angular primitives the realizer
**consumes**: `AngularAverageOperator`" at `sn/boundary/__init__.py:10`;
the realization-map line at `sn/boundary/realizer.py:50`; the live
`WhiteBoundary`-arm comment at `realizer.py:774`; two at
`geometry/boundary/__init__.py:241,349`) plus a claim INVERSION at
`geometry/boundary/_factors.py:1050` — `SpecularReemission.is_adjointable`
still reads *"TRUE, **unlike the Lambertian's**"* three lines after the
Lambertian's flipped to `True`.

How to apply: split a retirement's doc sites into THREE tense classes
(present-false / past-history / falsified-prediction) and check a
prediction's mechanism, host and phase separately; verify each FACTOR of a
factored operator against its own live method, never via the composite's
measurement; when un-sentineling would move a generated artefact you don't
own, keep the directive and rewrite the rationale; re-run the underlying
tool before crediting a baseline build warning.

---

## L-050 — A brief's "measured, post-carve" number can be a PRE-carve number in disguise; re-measure BOTH sides of a carve, and never let a 12-dp printout stand in for bit-exactness

The brief told me the post-carve two-channel equivalence was "to solver
tolerance, NOT bit-identically: `|φ_D0 − φ_C|_inf = 1.998e-13`,
`array_equal = False`", and asked me to document that as the current
state with a `⚠`. It was a real measurement — of the **pre**-carve tree,
between a delivery through the operator block and a delivery through the
source channel (two structurally different computations reaching one
fixed point). Post-carve the two channels are the SAME float program, so
the difference collapses to **exactly `0.0`** with `array_equal = True`
on every fixture and both inner solvers. Documenting the brief's framing
would have published the inverse of the carve's most important
consequence — and would have justified an `rtol` gate that is
structurally BLIND to the defect the carve removed (`2.9e-14` relative
sails through `rtol = 10 × inner_tol`).

Second, worse trap in the same brief: it asserted the converged inflow
trace is exact — "measured `5.000000000000` / `2.500000000000` /
`0.000000000000` at 12 dp on both fixtures" — and specified a keystone
gate as one `assert_array_equal` parameterized over BOTH inner solvers.
Measured: exact on source iteration (the sweep *writes* the seed, so `q`
arrives as a copy — `0.0` deviation, every fixture, every tolerance) and
**NOT** exact on Krylov (GMRES *solves* for the trace rows, so the
reading carries the iteration residual: 1–23 ULP at `inner_tol = 1e-13`,
and **27 580 ULP** at `1e-10`). Twelve decimal places of `2.500000000000`
cannot resolve `8e-15` at 2.5 — the printout that "proved" exactness was
blind to the effect by three orders of magnitude, and the plan never
caught it because the pre-carve Krylov leg *raised* instead of producing
a reading. The specified gate would be red on Krylov for a reason with
nothing to do with the claim.

Third: my OWN verification probe hid a failure the same way. A widened
bit-identity check used a bare `assert` inside a script I ran under
`python -O` — vv Mode 8, in my own instrument. It reported clean; the
assertion never executed.

How to apply: (1) when a brief hands you a "measured" number about a
change, ask WHICH SIDE of the change it was measured on — a pinned
worktree at the pre-carve commit (`git worktree add <dir> <hash>`) makes
both sides cheap, and re-measuring both turns one asserted number into a
before/after table that is the strongest content on the page. (2) A
venv's setuptools *editable* install hooks `sys.meta_path`, which
OUTRANKS `sys.path` — `PYTHONPATH=<worktree>` silently loads the MAIN
tree instead; strip the editable finder before importing, and PRINT
`orpheus.__file__` as the proof (my first pre-carve run measured the live
tree and looked plausible). (3) Never accept a fixed-decimal printout as
evidence of bit-exactness — assert `x == v` or print `float(x).hex()`.
(4) An exactness claim that holds on one inner solver and not another is
a per-leg gate, not a weakened gate: say `assert_array_equal` on the
exact leg and `rtol = SAFETY × inner_tol` on the iterative leg, and say
explicitly "do not relax the exact leg to match". (5) Run your own
probes WITHOUT `-O`, or raise instead of asserting.

---

## L-051 — a "measured, already-done" brief: reproduce it, and expect the
## reproduction to move the SCOPE (#341 "regular splitting" corpus sweep)

**Context.** Brief: strike the false term *regular splitting* at "seven
live doc sites"; the code sites were "already done — do NOT redo"; the
gate was `-W` EXIT 0. Delivered: 9 doc sites + 2 missed code sites + a
blocking build ERROR + a second, independent false claim.

**(a) A multi-word term takes an INFIX — grep each token, not the pair.**
The brief's grep was `"regular splitting|regular-splitting"`. The corpus
also spelled it **`regular matrix splitting`**, so two live present-tense
sites were invisible: `cartesian_multid.rst:3840` ("exactly a **regular
matrix splitting**") and `history.rst:836`. One `grep -rn "regular"` over
`docs/theory` (37 hits, ~30 seconds to triage) found both. This is the
cheap mechanic under theme 4's "grep the CONCEPT": for an
`<adjective> <noun>` term, grep the ADJECTIVE alone and triage, because
the noun phrase routinely grows a word in the middle.

**(b) "Already done in code" ≠ "gating green" — measure the `-E` baseline
before believing the brief's premise.** The `-E` baseline was **EXIT 0,
1 diagnostic: `ERROR: Malformed table`** in
`orpheus/sn/solver.py:docstring of solve_sn_fixed_source:86` — the ⚠ rate
table the SAME upstream code pass had just added, with every data row's
numerals straddling the `===`-defined column gaps (offsets 44 and 57).
So the brief's own acceptance gate (`-W` EXIT 0) could not have passed on
the tree it handed me. Fixing it was blocking AND in scope (same issue,
doc-only). Diagnose a simple table by reconstructing the column spans
from the separator with `re.finditer(r'=+', sep)` and flagging non-space
characters inside the gaps — instant, and it names the offending offsets.

**(c) ⭐⭐ A ratio is a ratio OF AN OBSERVABLE — ask which one before
citing it.** The brief handed me a spectral result and an investigation
memo whose §6 argued a conclusion from `n_GS/n_J` values. Those were
**ρ-DERIVED** (`ln ρ_J / ln ρ_GS` from an ARPACK eigen-solve of `M⁻¹N`);
every table already published in the corpus reports **SWEEP COUNTS** from
a real solve. `[M]` I re-measured five memo rows as sweep counts with a
control first: the control reproduced the published `1631 / 838 = 1.946`
**exactly**, and then **4 of 5 rows disagreed in SIGN** (memo `0.576`
"G-S wins" vs measured `1499/599 = 2.503` "G-S loses"). Both instruments
are individually sound — the memo validates ρ against a fitted residual
decay to 4 decimals. They simply measure different things: when
`ρ_GS ≈ ρ_J` to a fraction of a percent, the asymptotic ratio is wildly
sensitive while the sweep count is dominated by the transient, the
residual constant, and (per the memo's own §1.3) a frozen null-space
component the stopping test cannot see. **Rule: publish only the
observable you measured, name it in the caption, and never let a
rate-ratio and a count-ratio share a column heading.** The memo's
CONCLUSION survived (I established it independently on the other side);
its specific rows did not.

**(d) A technical term with a theorem attached is the ONLY carrier —
absence of a paraphrase is the finding, not reassurance.** I grepped
`Varga`, `comparison theorem`, `Stein-Rosenberg`, `monotone`, `M-matrix`,
`no slower than`, `never slower`, `at least as fast`, `ρ_GS`, `ρ_J` and
every `splitting` in `docs/theory`. The corpus **never once** wrote the
guarantee in prose. That is exactly why the word survived nine sites: a
reader who knows Varga supplies `ρ_GS ≤ ρ_J` silently, and a reader who
does not sees a decoration. So for a named theorem-bearing term, every
occurrence is load-bearing and none can be dismissed as decoration —
and the correction owes the corpus ONE place where the theorem is stated
and its failing hypothesis named, or the next author re-adds the word.

**(e) One home for the reason, `:ref:` from the rest.** The brief said
"point at the canonical code warning rather than restate the derivation".
Resolved by Cardinal Rule 3: a module docstring is a construction-site
note, a theory page is the brain. New H3
`sn-boundary-gs-not-regular` in `cartesian_multid.rst` (the page that
owns the boundary-G-S schedule) carries the whole derivation; the other
8 doc sites and 6 code sites `:ref:` it. Verified in built HTML: 5 of 5
cross-doc referrers resolve to real `href`s (a cross-doc `:ref:` renders
plain text with NO warning — the build cannot tell you this).

**(f) Importing algebra from a memo/code into a theory page imports its
SYMBOL COLLISIONS.** `A_a` (face area) collided with `A = L+C−S−B` (the
loss operator the whole section is about); `Σ` (transmission matrix)
collided with `Σ_t`. Resolution that satisfies the ratified
internal-consistency directive: **keep the code's spelling** and pay for
it with an explicit `.. note::` naming both overloads and their
disambiguators ("`A_a` always carries its axis subscript"; "`Σ` never
carries a `t`/`s` subscript"), stating that consistency with the
construction site outranks the local awkwardness. Do not silently rename
into the docs — that creates a code↔corpus spelling twin.

**(g) `.. (vv-status rationale)` is NOT machine-read — verified, not
assumed.** The scanner regex is `^\.\.\s+vv-status:(.*)$`
(`tests/_harness/audit.py:405`), which the `.. (vv-status rationale)`
comment does not match. So the rationale is free prose and the directive
line is the contract. Self-check without pytest:
`A._scan_theory_equations(Path("docs/theory"))` → fields
`all_labels / documented / skipped / violations`; I read **0 violations,
860 labels, 532 documented**, matching the auto-regenerated matrix
(531 → 532 from my one new label).

**(h) Scope ruling for an OPEN issue's investigation memo.** Published
(i) what the brief handed me and I re-derived (the `Σ = (2/D)1wᵀ − I`
spectrum; checked numerically at `d ∈ {2,3,4}`, plus the step-differencing
contrast `{0}^{d−1} ∪ {(D'−Σ_tV)/D'}` which makes the undamped subspace a
property of the DIAMOND closure, not of transport), and (ii) the two
counterexamples I measured myself. **Withheld** the memo's octant-order
law, its 25-pattern enumeration and its `max_a L_a > Σ_b L_b` predicate:
live findings with an unacted recommendation on an OPEN issue are the
main agent's to publish. Tombstoned, did not delete, the refuted
`ndim` reading — in all three places it appeared (theory page, a sibling
page's Key Facts, and the production docstring's user-facing
recommendation), because a half-done correction leaves a page
contradicting itself (vv anti-pattern #21's aggravator).

---

## L-052 — two dead-ref instruments DISAGREE BY DESIGN; the disagreement IS the triage

Task: the `nexus dead-references` SessionStart hook reported **21 targets /
30 sites** while `tools/check_docstring_xrefs.py` reported **0 dead across
14 914 roles**. The instinct is "one of them is broken". Neither was.

**The three-part scope story (measured 2026-08-09, `refactor/operator-strategy-layers`):**

1. **Trees.** The gate's default `roots = ["orpheus", "tests", "docs"]`
   (`tools/check_docstring_xrefs.py:199`). `examples/`, top-level
   `derivations/`, `scratch/` and `tools/` are never walked. Nexus walks
   the whole project. ⇒ 7 of 30 sites invisible to the gate for this
   reason alone. NOTE the gate's NAME understates it: it DOES read `.rst`
   whole-file (`iter_text_blocks`), so `doc:` sites in theory pages ARE
   in scope — "docstring" in the filename is wrong by omission.
2. **Roots.** `DECIDABLE_ROOTS = ("orpheus","numpy","scipy","sympy","pytest","matplotlib")`
   excludes `tests`, `tools`, `derivations` — all three ARE importable, so
   this is a coverage gap, not an impossibility. And UNQUALIFIED refs
   (`:func:`compute_G_bc``) are skipped BY DESIGN (the tool refuses to
   emulate Sphinx's module-context resolution rather than manufacture
   false positives). ⇒ 6 more targets.
3. **⭐ The semantic split, and the one that matters.** The gate resolves
   by **IMPORT**; nexus resolves by **RENDERED TARGET**. A live module
   that no `api/` page `automodule`s is *resolved* to the gate and *dead*
   to the hook. That is the ENTIRE bucket-C class (12 of 21 targets here),
   and both tools are right: the symbol exists (gate) and the role has no
   link target (hook). Neither number is a bug; **the set difference is
   the triage** — hook-minus-gate ≈ "un-surfaced but live" (issue #302),
   gate-and-hook-agree ≈ "genuinely retired or moved".

⇒ never write "all trees at 0" in a memory index from a gate run. Write
the trees, the roots, and the resolution SEMANTICS, or the claim is
present-tense-false the day someone points a second instrument at it.

**Two false-negative classes found in the gate while establishing this:**

- **PEP-420 namespace packages resolve.** `orpheus/derivations/continuous/{pn_method,
  spn_method,spectral_collocation,spectral_resolvent,escape_probability}/`
  each contain ONLY a `README.md` — no `__init__.py`, no Python.
  `importlib.import_module` succeeds anyway (`__file__ is None`, 0 members),
  so `resolve()` returns True. A `:mod:` role at such a target can NEVER
  resolve in Sphinx. Discriminator to add if the gate is ever hardened:
  `mod.__file__ is None`.
- **⭐ A role wrapped INSIDE its dotted path is invisible to BOTH tools.**
  `:func:`~orpheus.numerics.eigenvalue.\n        power_iteration`` — docutils
  collapses the newline+indent to a space, so the target becomes
  `...eigenvalue. power_iteration` and never resolves; the gate's
  `extract_target` returns `None` on any target containing whitespace, so
  it SKIPS it. Measured tree-wide: **15 such roles** (orpheus ×13,
  tests ×1, examples ×1). The discriminator is NOT "role spans lines" —
  ~180 multi-line roles are FINE because they break at the
  `display <target>` boundary. The regex that finds only the broken class:
  `\.\s*\n\s*\w` (or `\w\s*\n\s*\.\w`) inside the pre-`<` head.

**Bucket discipline that made the triage cheap.** (A) MOVED / (B) RETIRED
/ (C) NEVER-AUTODOC'D is decided by TWO probes, not one: does the symbol
import (A/B vs C), and is its module in the `automodule` set
(`grep -o "automodule:: (\S+)"` over `docs/**/*.rst`, 49 here). Anything
that imports AND is un-automodule'd is C — hands off, it is #302's.

**⚠ The reported TARGET NAME can be an artifact of a THIRD tree.** Six
"dead targets" were named `orpheus.derivations.peierls_geometry.*` — a
module deleted at `bda76faf`. No doc page contains that string. The node
exists ONLY because three `scratch/derivations/diagnostics/*.py` still
`import` from the dead path, minting an `unresolved` ast_only node; nexus
then attached the theory pages' UNQUALIFIED `:func:`compute_G_bc`` roles
to it by suffix match. The functions are alive at
`...peierls_nystrom.geometry`. ⇒ before believing a dead target's NAME,
`SELECT source,type FROM edges WHERE target=?` on `graph.db` and read the
edge TYPES: `documents` = a doc page, `references` = a docstring,
`type_uses` = a **type annotation** (i.e. real code, not prose), `calls` =
an import/call that MINTED the name.

**A `type_uses` site is a CODE bug wearing a doc-ref costume.** Three of
the 30 "sites" were `TYPE_CHECKING` imports in `examples/` and
`derivations/`, i.e. `ModuleNotFoundError` at runtime (or a pyright red)
— not prose at all. Fix the import, then ask whether the annotated body
still matches the successor type: `examples/discrete_ordinates/plotting.py`
annotated `result: SNResult`, and the live `Solution` has no `.geometry`,
no `.eg`, and a `ScalarFlux` with **no `__getitem__`** — so a bare repoint
would have replaced a dead name with a live LIE. Published the repoint
PLUS a measured `.. warning::` naming the three surviving stale accesses.

**The unit of repair is the TARGET, not the reported site.** Nexus counted
3 sites for `orpheus.sn.geometry.SNMesh`; the live tree had **13**
(12 × `derivations/diagnostics/*.py` + 1 × `examples/`), because
dead-references only reports `documents`/`references`/`type_uses` edges
and most were plain `import` statements. Fixing 3 and leaving 10 would
have been the half-correction anti-pattern. Also: nexus counts doc sites
**per PAGE**, so "compute_P_esc ×2 sites" was 9 role occurrences.

**Adjacent classes found, REPORTED not fixed** (each is its own pass):
(a) `derivations/diagnostics/` is broadly bit-rotted — after my repoint,
**15 of 39** files still carry an unresolvable `orpheus` import
(`orpheus.sn.quadrature`, `orpheus.sn.operator`, `orpheus.sn.boundary_realizer`
are all gone); say so, or the repoint reads as "these scripts run now".
(b) `docs/theory/references/trajectory_resolvent.rst` carries **31 stale
`:file:` paths** — a page rename (`peierls_greens.rst` →
`trajectory_resolvent.rst`) was applied to the TEST FILENAMES in prose,
which were never renamed. Raw paths warn NOWHERE. Fixing the 1 in my
brief and leaving 30 would have made the page contradict itself; report
the class.

**Concurrency hazard on a busy branch.** The after-measurement grew ONE
new dead target that was not mine: another agent's #340 N2b docstring
landed in `orpheus/cp/solver.py` mid-session with a wrapped role. `git
status` at session start is the only way to attribute; state the
before/after BOTH ways (raw, and net-of-not-mine). Same for
`docs/theory/verification/matrix.rst`, which my `-E` rebuild regenerated
to absorb their +1 foundation test — legitimate by-product, report it,
never revert it.

**Result.** 21 → 14 targets / 30 → 17 sites (13 + 1-not-mine); the 8
cleared are exactly buckets A and B. Gate stayed at 0 dead. `-E -W`-clean
baseline and after builds both EXIT 0 with a byte-identical (empty)
WARNING/ERROR/CRITICAL set.

**Quality self-assessment.** Derivation depth n/a · Cross-references 5 ·
Numerical evidence 4 (before/after counts, per-tree measurements, the
15-role scan — no convergence table is possible on a reference-integrity
pass) · Failed approaches 4 (recorded WHY the peierls_geometry name is an
artifact, and why 3 adjacent classes were left) · Code traceability 5 ·
Derivation source n/a.

---

## Quality self-assessment rubric (Directive 3)

Rate each output 1–5 and log the weakest dimension in the return:
Derivation depth · Cross-references · Numerical evidence · Failed
approaches · Code traceability · Derivation source (from `derivations/`,
never hand-written). The recurring weak dimension on TERMINOLOGY/ROUTING
fixes is "numerical evidence" — structurally absent (no flux moves → no
convergence table), not a deficit; say so rather than manufacturing one.

---

---

## L-053 — a per-site adjudication TABLE is an instrument too: audit its SKIP clause, its "retired" verdicts, and its `hasattr` evidence

**Context.** #346 W1: a 91-site / 64-distinct ruling table (`scratch/issue_346_w1_adjudication.md`),
built by `w1_list.py` over `graph.db` + the `check_docstring_xrefs` resolver, sorted into
Groups 0–5 with a per-site fix already decided. The brief said "do not re-derive the rulings;
apply them", and warned the instrument had mislabelled things three times. Applying it
faithfully still produced **five ruling corrections**, each from a different structural cause.

### (a) The SKIP clause is a false-NEGATIVE source, symmetric to the false positives the brief warned about

`w1_list.py` keeps a site only if the target is absent from the graph by **tail match**
(`s.endswith("."+t)`). That clause is what suppressed the brief's *known* false positives — and
it equally suppressed **alive-but-unqualified** roles whose tail happens to exist:
`:class:`Field`` on the very line I had to edit, `:class:`AngularFlux``,
`:class:`~geometry.mesh.BC`` (8 lines above three BC-member sites I *was* given), and the whole
`~geometry.*` / `~transport.*` prefix-omission family across ≥6 pages. `7d7661b2`'s own commit
message names it: **1431 relative roles across 49 `.rst` pages**, deliberately deferred.

⟹ When a work-list is generated by *"report if NOT in X"*, the deliverable owes a sentence on
what the NOT-clause hides. And the boundary decision is **file-convention**, not tidiness:
qualifying `~geometry.mesh.BC` on one line while lines 15/16/59/103 of the same file keep the
project-wide short form makes the file *less* internally consistent. I fixed `Field` (whose page
convention IS fully-qualified) and left the `~geometry.*` family, reporting the count.

### (b) "RETIRED → literal" must first ask *retired, or MOVED?*

Seven sites (`peierls_slab` ×6, `peierls_cylinder`) were listed under "Group 4 — RETIRED …
past-tense history → literals". `git log --oneline --diff-filter=D --all -- '*peierls_slab*'`
→ `bda76faf refactor(derivations): reorganize into common/discrete/continuous architecture` —
a pure `git mv`. The module is **alive** at
`orpheus.derivations.continuous.peierls_nystrom.slab`, and the page's own next line says it is
*"retained indefinitely, not retired, as an independent cross-check implementation"*. Six of the
seven sentences are PRESENT tense or IMPERATIVE ("modifications … should preserve"). A literal
would have thrown away a live link and past-tensed a live module.

⟹ Run `--diff-filter=D` on the old path **before** accepting any "retired" verdict. One command.

### (c) A dead `:attr:` can be a TRUE claim — `hasattr(Cls, x)` is the wrong probe

The table's `[M]` for Group 3 item 2 was `hasattr(SNMesh, "mesh") is False`, ruling
*"After C5.1, `SNMesh.mesh` is inbound provenance only … `None` when built from axes"* as
present-tense-FALSE and demanding a rewrite. Constructing the object says otherwise:
`SNMesh(Mesh1D(...), Quadrature.gauss_legendre(4), {0: mix}).mesh` → a `Mesh1D`. The attribute is
set on the **instance** by the base `MaterialMesh._init_data`, with no class-level annotation —
so `getattr(cls, …)` fails, autodoc mints no target, and the role renders plain text while the
*sentence* is exactly right. `augmented_mesh.py` even spells the doc's other half verbatim:
`mesh = legacy_mesh_from_axes(...) if len(axes) <= 2 else None`.

⟹ Group-3 ("prose is false") and Group-0 ("alive, unlinkable") are separated by **constructing
the object**, never by a class-level `hasattr`. The fix is a literal + a live `:class:` ref +
one clause saying *why* it cannot be a role — not a rewrite of a true paragraph.

### (d) A "role misuse → literal" ruling loses to the page's own convention six lines up

`:meth:`B.apply`` was ruled Group 5 ("names an instance, no qualification can ever make it
resolve"). Six lines above, the same page already writes
`:meth:`B.apply <orpheus.sn.operators.boundary.SNBoundaryOperator.apply>`` — `B` *is* the
`SNBoundaryOperator` and `.apply` *is* a real method. Literalising it would have been a
regression against a live, working, explicit-title link in the same paragraph block.

Conversely `:class:`ReflectiveBoundary.apply`` / `:class:`WhiteBoundary.apply``, ruled Group 2
(*repoint to `:meth:`~orpheus.geometry.boundary.reflective.ReflectiveBoundary.apply``*), resolve
**DEAD** — those classes carry `realize` / `geometry_map` / `response_kernel` / `source`, no
`apply` — and sit in a subsection headed *"What was tried and rejected"*, twenty lines under the
page's own sentence *"Calling `law.apply(psi)` raises `AttributeError` at runtime"*. Repointing
would have minted a live role at a method the page explicitly says the type system removed.

⟹ Two adjacent `X.apply` sites, opposite correct answers. **Resolve `Instance.method` by asking
what the instance's TYPE is and whether that type has the method** — never by the shape of the
target string.

### (e) Fixing one role can drop you into a page-wide SYMBOL-INDEX collision — measure, mark, refuse

The Group-4 row `_ki4_lookup` (one table cell) opened this:

| claim | site | measured |
|---|---|---|
| cylinder kernel is `Ki_4` | `collision_probability.rst` ×8 + 5 sibling pages | code ships `_ki3_mp`; `[M]` `_ki3_mp(0) = 0.7853961`, `_ki3_mp(1) = 0.2378450` = **standard** `Ki_3` (`π/4`, `0.2378450`) |
| `:eq:`ki3-def`` defines `Ki_3` | `= ∫₀^{π/2} e^{−τ/sinθ} sinθ dθ` | that is the **standard `Ki_2`** (one power of `sinθ` short) |
| `F(0) ≈ 0.4244` | geometry-comparison table | `= 4/(3π)`; matches **neither** convention |
| `self._kernel = Ki_4` | prose | live line is `self._kernel = _ki3_kernel` (`= _ki3_mp`) |
| "20 000-point lookup table" | Key Facts line 24 | retired at Phase B.4 / #94 → degree-63 Chebyshev interpolant of `e^τ Ki_3(τ)` |

So the page runs a **consistent internal index one below the standard**, under which its `Ki_4`
IS the shipped function — the *structure* of `:eq:`second-diff-cyl`` / `:eq:`self-cyl`` is right
and only the *symbol* is wrong; but `ki3-def` is separately off by a power, and `F(0)` is simply
wrong. Those three labels carry **64 / 24 / 54** `verifies()` tests
(`docs/theory/verification/matrix.rst`).

⟹ **The tests cannot see this.** They pin the code's numbers; the equation's *symbol* is outside
everything they measure — a documentation-side Mode-12: the measured functional (the test) is
invariant to the error class (the name). No gate exists, and none can be added without changing
what the equations say.
⟹ Correct move: repoint the role, **measure** the discrepancy into a `.. warning::` with an
anchor, fix only the unambiguously-wrong number (`F(0)`), and explicitly REFUSE the 142-test
re-indexing as a numerics adjudication. Renaming across ~30 corpus sites inside an xref pass
would be an unreviewable physics edit riding in a hygiene commit.

### (f) "Qualify so it resolves" is TWO claims, and they come apart — measure which one you bought

`[M]` post-build, live `href` counts for repointed targets:
`EigenvalueSolver` **43**, `Field` **30**, `numpy.array_equal` **6** — real links.
`KEigenvalue` **0**, `AngularFlux.zeros_on` **0**, `BC.vacuum` **0**, `SNMesh.axes` **0**,
`peierls_nystrom.slab` **0** — still plain text, because their modules are `:noindex:`-`autoclass`'d
(`docs/api/geometry.rst`) or not `automodule`'d at all. The gate and the graph are satisfied; the
reader is not.

⟹ Never write "these now resolve" unqualified. Write "**import-** and **graph-**resolvable;
N of M also render as links, the rest await #302 surfacing" — and check with
`grep -o 'href="[^"]*#<target>"'` in the built HTML, which is one command.

### What worked (repeat)
* **Every structural assert on the in-memory result BEFORE any write.** Four scripts aborted on a
  bad anchor/count and left the tree untouched — no `git checkout` recovery, which this tree forbids.
* **Anchor on CONTENT, never the table's line numbers.** They shifted by tens of lines after the
  first insert into `boundary_conditions.rst`; a line-keyed second pass would have corrupted it.
* **`resolve()` from the gate itself as the pre-flight probe.** Every target import-verified
  *before* it was written, so the final gate run was a confirmation, not a discovery.
* Assert failures on my OWN new text (a literal `` ``fn_method.core.x_function`` `` inside the note
  that says the module does not exist) are a good sign the assert is tight — narrow it to the ROLE
  form, don't weaken it to a substring.

---

## L-054 — a landed CARVE's docs pass: grep the CONCEPT in every SPELLING, and the brief's page count is an over-count AND an under-count at once

**Instance.** Q5.6.4 (branch `refactor/operator-strategy-layers`, `3dda18ca` + `d5067c4d`):
the cylinder's angular cell partition moved from the η-midpoint (chord) to the midpoint in ω,
and the `[½,1]` τ absorber retired. Brief asked for a label rename
(`morel-montry-clamp` → `morel-montry-closure`), a `cases` split into closure + partition,
a two-page claim repair, and a concept sweep it enumerated as **17 pages**.

### The page count was wrong in BOTH directions, from ONE cause: I grepped SPELLINGS, not the concept

The brief's list came from a grep of `clamp|absorber|[½,1]|τ_raw`.

* **Over-count — 11 of 17 pages were FALSE POSITIVES**, all from one word: `absorber` is
  ALSO a *material* ("pure absorber", "thick absorber", "cavity-absorber"), and `clamp`
  is also a GMRES `restart` clamp, an interpolant clamping to zero, and the LD weight's
  own legitimate `[½,1]`. Real M-M-clamp content lived on **6** pages
  (`foundations/structured_geometry`, `verification/sn`, `methods/sn/{curvilinear_one_group,
  curvilinear_numerics,history,curvilinear_multigroup}`). Cost of not checking: a sweep that
  "fixes" a physics term.
* **Under-count — one page the brief did NOT list carried a present-tense-false BOUND.**
  `methods/sn/angular_quadrature.rst:369` read *"the raw march-start weights satisfy
  `τ_raw ∈ [1/5, 4/5]` with the **bit-exact** reversal identity"* — both halves stale
  (now `[1/4,3/4]`, 0.5–12 ULP). It was invisible to the brief's grep because the page
  spells it **`\tau_{\rm raw}`** — a LaTeX spelling that matches neither `tau_raw` nor
  `τ_raw`. ⟹ **for any math symbol, grep at least three spellings: the ASCII identifier
  (`tau_raw`), the Unicode prose form (`τ_raw`), and the LaTeX role body
  (`tau_{\rm raw}` / `tau_{{\rm raw}` / `tau^{\rm raw}`).** I found it only because I ran a
  *residual* sweep on the NUMBER (`tfrac15|tfrac45`) after the build was green.

### A retirement propagates to a page's BOUNDS and its DIAGNOSTIC-CLASS claims, not only its symbols

Three shapes I had to fix that no symbol grep reaches:

1. **A numeric BOUND is a claim about the retired object.** `[1/5,4/5]` and "bit-exact
   reversal" were properties of the CHORD partition. Present in a theory eq
   (`morel-montry-folded-arc`, a live `verifies()` target — keep the label, rewrite the
   BODY), in a production docstring (`directional.py:717`), and on an unlisted page.
   `[M]` I re-measured all five rows from the live producer before publishing:
   `[0.292893,0.707107] / [0.259892,0.740108] / [0.252425,0.747575] / [0.250603,0.749397]
   / [0.250151,0.749849]` with 0.5/0.5/2/7/12 ULP — matched the live test docstring exactly.
2. **A section's THESIS can rest on a requirement that no longer exists.** `verification/sn.rst`
   argued the cylinder floor "structurally blocked" partly via *"No partition
   (midpoint/cumulative-weight/ordinate-interior) gives `τ_raw ∈ [½,1]` with bounded edges."*
   `[½,1]` was never a requirement — it was the ABSORBER's box. So the paragraph searched for
   a partition satisfying a condition no reference imposes, and a partition satisfying the
   REAL predicate (P3, `[0,1]`) had just shipped. Fix = `.. note:: Retraction` quoting the
   old text, keeping the surviving half (the azimuthal-duplication argument, which is about
   the MARCH not the partition) and re-deriving the cumulative-weight observation as a
   **P3** failure with the measured ladder.
3. **A "the floor is independent of X" claim can be REFUTED by the very carve that retires X.**
   `curvilinear_one_group.rst` said the residual floor "is independent of the spatial closure,
   the default, and the τ-clamp". `[M]` removing the clamp moves the floor 1.8–3.4×. The
   quadrature-scaling attribution survives; the independence clause does not. ⟹ **after a
   retirement, grep the retired thing's name inside "independent of" / "unaffected by" /
   "does not depend on" sentences** — a negative claim about it is exactly the claim the
   retirement can falsify.

### Publishing a `[M]` number the commit hands you: check its CONFIGURATION, because agreement can DEGRADE with refinement

The commit stated the closed form `τ_m = ½ + ½·cot(ω_m)·tan(Δω/4)` "verified to `1.67e-16`".
`[M]` I re-measured on `folded_product(n_mu=4, n_phi)`, max over all four levels:
`1.1e-16 / 2.2e-16 / 7.8e-16 / 7.4e-15 / 2.3e-14` at `n_φ = 4/8/16/32/64`. The single figure
is a small-`n_φ` reading; publishing it bare would have implied a machine-epsilon identity
that **degrades two orders** by `n_φ=64` — and the shipped gate knows this (`atol=1e-13`).
⟹ **an agreement number is a LADDER unless proven flat. Measure the ladder, publish the
ladder, name the gate's tolerance, and say a finer row must widen it (`vv-principles` #16).**
Same trap on a per-level ratio: the code's "(spread 0.30→1.53)" looks like a convergence
sequence and is actually the **`n_φ=16` row across one level** (`[M]` 0.59–1.41 at 8,
0.30–1.53 at 16, 0.08–1.57 at 64) — read a two-number "spread" as *one* measurement until
you have re-run it.

### A retired MODULE whose FUNCTION survives by name is a semantic trap, not a repoint

`derivations/discrete/sn/contamination.py` retired into `angular_differencing.py`, and
`morel_montry_weights` **survives by name** — but it now DELEGATES to production, so it is
**no longer an independent reference**. Four `:func:`…contamination.morel_montry_weights``
refs on `curvilinear_one_group.rst` sat inside sentences crediting it as
*"the structurally-independent reference"*. A mechanical repoint to the live path would have
produced four working links to four false COVERAGE claims (lessons §7). Correct treatment:
**de-role to a past-tense literal in the historical narrative + ONE new anchored
`.. note::` (`sn-tau-reference-migration`) that states the delegation, why it was chosen
(a reference must not be free to drift into a second definition — which is exactly how the
old module's cylinder arm went wrong, `[M]` τ off by 6.8e-2), and what each arm compares
against NOW** (hand-authored: an inline cumulative-weight expression on the sphere, the
analytic closed form on the cylinder, with the retired chord convention as its
`vv-principles` #19 negative control). One anchor, four referrers — not four rewrites.

### A BLIND diagnostic earns a `.. warning::` on the page that recommends it

`contamination_beta` is *identically zero on a σ_y-folded arc for ANY antisymmetric edge set*
— `[M]` production `+6.94e-18`, edges×0.5 `+3.47e-18`, edges **CUBED** `+1.73e-18`, random
antisymmetrised `−3.47e-18`; only breaking antisymmetry moves it (`−3.53e-03`). It certified
a convention that **diverges the solve**. The theory page said "the module computes β for any
quadrature and geometry" with no caveat. ⟹ when a carve discovers a published diagnostic is
Mode-12 blind on the shipped rule, the doc owes a `.. warning::` with the *garbage-passes*
table, and must name the instrument that DOES discriminate (here ν-closure, solve-free —
`[M]` reproduced the whole table exactly: arc/chord `1.000000`, clamped `1.016389→1.000030`,
`τ≡½` `1.164784→1.002412`).

### Two build defects my own new text introduced — both from writing prose, both -W-caught

* `*"… (*:math:`X`*), …"*: an italic run interrupted by a `:math:` role gives
  **"Inline interpreted text or phrase reference start-string without end-string"**.
  Fix: escape the seam — `(*\ :math:`X`\ *)`.
* A hand-aligned simple table whose header cell (`:math:`n_\varphi` 8→16`, 22 chars) overflows
  its `===` column (18) → **`ERROR: Malformed table`**. Fix: use `list-table` for anything with
  a role in a header cell. **Never hand-align a simple table containing a `:math:` role** —
  the role's source length is not its rendered length and the column arithmetic is invisible.

⚠ Both were caught only because I re-ran `-E -W` after writing. Two of my four builds were
**wasted** by launching before the last edit landed: ⟹ **finish EVERY edit and EVERY residual
grep, then build once.** Sequence the session as: baseline `-E -W` → all edits → all greps →
xref gate → AST doc-only proof → ONE verification build.

### What went right, keep doing

* **Baseline re-measured this session: EXIT=0, ZERO W/E/C** — so the gate was count-**equal**,
  not count-unchanged-from-a-quoted-number. Final: EXIT=0, zero W/E/C.
* **`tools/check_docstring_xrefs.py orpheus tests docs --quiet` → `DEAD TARGETS: 0`** before
  and after (11 272 → 11 274 decidable). This is what proved the retired-symbol sweep, since
  `pole_angular_closure` has no `automodule` and no build severity can see its refs.
* **AST-with-docstrings-stripped vs `HEAD`** proved both production edits doc-only
  (`reduced_operator.py`, `directional.py` → `identical = True`).
* **Rendered-href check in the built HTML** for every new `:ref:`/`:eq:` I minted
  (`sn-tau-absorber-retirement` 10 hrefs, `angular-cell-partition` 11, etc.) — cross-doc
  `:ref:` renders plain text with no warning, so this is the only proof.
* **Three orphaned July HTML files** (`theory/structured_geometry.html`,
  `theory/discrete_ordinates.html`, `theory/methods/sn/verification.html`) still carry
  `equation-morel-montry-clamp`. Discriminated by `test -f <source>.rst` → no source ⟹
  orphaned build output from a page split, NOT a live stale ref. Do not `rm -rf docs/_build`.
* **Programmatic splice with all structural asserts run BEFORE the write** caught two of my
  own mistakes (a surviving `morel-montry-clamp` on a line ABOVE my slice; a wrong expected
  ref count) with the tree left untouched.

### Scope discipline that held

`sn-tau-mm-raw` on `verification/sn.rst` is the SAME naming-honesty disease the brief flagged
(a label spelling "raw" on an equation with no raw/clamped distinction), 60 lines from the
rename. I did **not** rename it — a label rename has a V&V-matrix footprint and the brief
named exactly ONE. I fixed the equation BODY's notation (`\tau_n^{\rm raw}` → `\tau_n`),
left the label, and put the reasoning + the exact 3-site rename cost in a `.. NOTE` comment
beside its `vv-status` directive so the next session can execute it in one pass. ⟹ **when a
carve reveals a second instance of the disease it fixes, fix the CLAIM and document the
rename cost in place; do not annex the rename.**

---

## L-055 — the brief's measured number can be TRUE and its FRAMING already refuted, by a test module committed the SAME DAY

**Instance.** Q5.6.4 follow-on, branch `refactor/operator-strategy-layers`, 2026-08-11: a
citation-defect repair (Hébert-vs-BMC as the source of the weighted M-M `τ`) carried one
side-task — *"`docs/theory/verification/sn.rst` ~1186 asserts 'Positivity is never needed'.
We have just measured the HALF-ANGLE flux reaching `min ψ̂ ≈ −77` on a normalised cylinder
problem under the shipped convention. Determine whether the claim is scoped to the converged
SCALAR flux (fine) or is general (falsified)."*

### What I nearly published

I reproduced `−77.1643` exactly (`folded_product(4,32)` level 0, `exp(−6cos ω)` profile,
positive constant seed, through the production `compute_psi_half_per_level`), extended it to
all four τ conventions (`chord −229.7 / chord+absorber −23.3 / arc −77.2 / τ≡½ −24.2`), and
drafted a `.. warning::` ending in *"⚠ **Coverage gap**: `[M]` there is **no** ψ̂-positivity
gate on either arm."* Evidence: an untracked `scratch/` QA memo saying exactly that, plus my
own grep of the 15 `tests/` files mentioning `psi_half` for an assert-on-the-same-line.

### What the tree actually said

`tests/sn/sweep/curvilinear/test_psi_half_positivity.py` — **19 `foundation` rows, committed
the same day**, and it is a CHARACTERISATION module whose docstring *pre-emptively refutes the
`−77` framing*: `[M]` on a heterogeneous 2G cylinder with the **marched ψ½ seed** (the
production value path) ψ̂ is strictly **POSITIVE** — `+0.1337/+0.1286/+0.1287` at
`n_φ = 6/8/16`, i.e. 0.88/0.93/0.98 × `min ψ`. Only an **inconsistent** (zero) seed goes
negative — `−12.09/−16.35/−25.89` on the *same* converged flux — bounded by the recurrence's
worst partial amplification `A(M) = max_m Π(1−τ_k)/τ_k = 2.73/3.36/4.73`. ⟹ *the sign is a
property of the SEED's consistency, not of the scheme*, and my `−77` is an inconsistent-seed
statement. Reproduced all of it; 19 passed in 3.7 s.

⭐ **Rule: when a brief hands you a measured number to PUBLISH, grep `tests/` for a module
NAMED after the phenomenon before you write its interpretation.** A `scratch/` memo is by
construction OLDER than the tests it motivated, and a same-day characterisation module is the
likeliest home of the corrective framing. My line-based `psi_half` + assert grep missed it
(`vv-principles` #21: the subject and its assertion sit on different lines, and the module's
real evidence is in `pytest.fail` messages and the docstring, not in `assert` lines) — the
grep that found it was the read-only `tests/` sweep the brief asked for *for reporting*.

⭐ **Corollary — the correct verdict on the claim was a THIRD option the brief did not
offer.** Not "scoped to the scalar flux (fine)" and not "general (falsified)": the claim is
*general in wording, sphere-in-evidence, and substantively TRUE on the cylinder's production
path as a characterisation*. What is false is the word **"never"**. So the edit scopes the
heading (`The clamp buys no positivity on the SPHERE'S converged solve`), states the seed
taxonomy with both measured tables, keeps W1's conclusion standing (the clamp reduces the
excursion ~10× but does not remove it, and neither does `τ≡½` — the exposure is
`(1−τ)/τ`, a property of the *angular diamond family*), and points at the owning module. When
a brief offers a binary and the tree supports neither pole, say so and publish the third.

### `ref.ref` also fires for an anchor before a **BOLD RUN-IN HEADING**

`.. _label:` immediately above `**(2) Some Title.**  Prose…` gives
*"Failed to create a cross reference. A title or caption not found"* at every bare `:ref:` to
it — **5 hits across 4 files, all from one anchor**, and `-W` turns them into EXIT=1. A
run-in bold heading *looks* like a heading and is not one. Two fixes, both used here:
promote it to a real titled subsection (what I did for `sn-tau-source-of-record`), or use
explicit link text `` :ref:`the β-blindness warning <sn-tau-beta-diagnostic-blind>` `` (what
I did for an anchor that legitimately sits above a `.. warning::`). ⚠ **Do not open the new
title with `(1)` / `(2)`** — that is an enumerated-list marker in RST; use
`Correction 1 — …`.

### Build sequencing, again

Four `-E -W` builds where two would have done: baseline (0 W/E/C), a verify launched *before*
the last edit landed (wasted), a third that caught the 5 `ref.ref` warnings (earned), a
fourth that closed at EXIT=0 / 0 W/E/C. The wasted one is exactly L-054 §9's warning. The
earned one is the argument for never skipping the post-edit build even when the xref gate is
green — `tools/check_docstring_xrefs.py` reported `DEAD TARGETS: 0` while all 5 `ref.ref`
warnings were live, because it gates **Python-domain roles**, not `:ref:`.

### Doc-only proof and the generated-artefact by-product

9 production `.py` files touched, all proved DOC-ONLY by AST comparison against `HEAD` with
docstrings stripped. The `-E` build regenerated `docs/theory/verification/matrix.rst`
(9544 → 9628 tests) absorbing rows from **another agent's uncommitted `tests/sn/sweep/` work**
(`test_psi_half_positivity` +19, `test_angular_cell_partition` +56,
`test_tau_producer_equivalence` 5→14) — the L-008 by-product: never revert it, report it.

---

## L-056 — A labelled equation drifted from its own prose; a scope revocation landed mid-edit

**Task.** Repair `docs/theory/foundations/discrete_measures.rst` §"Quadrature selection
algorithm", which described the pre-2026-08-02 quadrature selector (4 stages, declared-tag
symmetry gate) while the code had shipped 5 stages with a computed `Sym(Q)`. Baseline and
final `-E -W` builds both EXIT=0 / **0** W-E-C. One tracked file edited: +415/−74.

### 1. ⭐⭐ The DEFINITION LIST is the tell that a labelled equation drifted

The commit that fixed the design (`e7d44f3c`, 2026-08-02) rewrote the geometry table, the
worked examples, the rejection messages **and the predicate quoted inside the equation's own
`.. (vv-status rationale)` comment** — and left the `.. math::` body alone. So the page
carried the corrected claim in a comment eight lines below the false claim it annotated.

The mechanical tell needs no code and no build: **the "where" list under the equation defined
`𝒟_Q` and `Sym(Q)` — two symbols absent from the equation — and omitted `G_Q`, which was in
it.** A definition list that does not match its own equation is a correction that stopped one
line short. Add this to the reading pass on ANY page with `.. math:: :label:` + a where-list.

Why nothing catches it: the label EXISTS, so every `:eq:` resolves, `-W` is silent at every
severity including `-n`, and the V&V matrix lists the label as covered — coverage is recorded
against the *label*, not against what the label says. This is `coding-standards`' "a labelled
equation is an API" clause, met in the wild.

⟹ I published the tell IN the page as a `.. admonition::`, not just in the fix. A reading
skill archived where the next reader of that equation will hit it.

### 2. ⛔ A mid-task scope REVOCATION on a file already edited: revert by re-editing, and
**publish the patch you backed out**

The brief said "correct BOTH sites", naming `registry.py:106-107`. I fixed it. A mid-task
addendum then said *"do not edit `registry.py` … confine yourself to the .rst"* — the
coordinator was editing that file concurrently.

- Reverted by **re-editing** (never `git checkout` — the tree carries uncommitted-by-policy
  state), then **proved it** with `git diff --quiet -- <path>`. Do not report "reverted"
  without that proof; an `Edit` round-trip can leave whitespace.
- The backed-out content is not lost: it goes in the RETURN, verbatim, as an apply-ready
  patch. Four `registry.py` sites were owed and the addendum named only two.
- ⚠ **`git status` then showed the file Modified again** — the coordinator's concurrent pass.
  Discriminate by grepping YOUR OWN distinctive strings (`grep -c "teaching artifact"` → 0),
  not by the porcelain flag.

### 3. ⭐ A tombstone that asserts the state of ANOTHER file is false the moment either file moves

My §"Why a registry" tombstone read *"…and the module docstring said the same thing … until
2026-08-14"*. After the revert that was FALSE (the module still said it); after the
coordinator's pass it would be false the other way. Both directions wrong from one sentence.

⟹ **Write a twin's history in the PAST tense of the CLAIM, never of the file's state.** Landed
form: *"The promise was minted twice: the module's own docstring stated it more concretely
still, as '…'. That duplication is itself the tell — a claim about a *rendering* mechanism had
never been checked against the rendering, and the module asserting it is the one that is not
rendered."* True whatever the other file does next. Same rule as §2's patch-in-the-return:
**your page may only assert what your page controls.**

### 4. ⭐ Fixing half a claim in one file is WORSE than fixing none (vv #21, met head-on)

The brief scoped me to the §starting at line 682. Three screens ABOVE it, the same page's
"Symmetry groups for quadrature invariance" section still opened *"Quadrature selection in
ORPHEUS therefore reduces to a containment check"* with the retired whole-group mapping
(slab → `SO(2)×σ_x`, sphere → `O(3)`) and closed by calling containment *"sufficient to
preserve every symmetry the geometry exhibits"*. Had I stopped at the brief's boundary the
page would contradict itself and be citable for EITHER sentence.

Repair shape, with the equation and its `vv-status: documented` label untouched (it is
`implements`-cited from `orpheus/numerics/symmetry.py:371,481`): **re-scope the equation to
what it actually is** — the order relation on the `O(3)` lattice, `H ⊆ K ⟺ ∀g∈H, g∈K`,
decided by `SubgroupOfO3.contains` — then a `.. warning::` saying *this relation is not the
selection gate*, with both reasons (a rule's symmetry is a question about NODES; the geometry
side is not one group), then the ⛔ preserving the retired sentences.

⟹ **After repairing a section, grep the WHOLE FILE for the retired predicate's spellings and
adjudicate every hit by tense.** Here: `G_Q`, `G_{\text{geom}}`, `G_{\text{quad}}`,
`four-stage`, `docstring narrates`. Final state: 6 survivors, every one inside a ⛔/history
admonition.

### 5. ⭐ A stale FORMULA can be right on a biased subset of the grid — spot-checking confirms it

Stage 2 gave the level-symmetric degree as `max(3, N−1)` and the positivity frontier as
`S_12` with `[M] −0.027 @ S_14`. Both were honest measurements **of the pre-#337 seed**
`μ₁² = 4/(N(N+2))`; #337 (`59bb38a0`, 2026-08-08) adopted the moment-matched root and moved
both. `[M]` at HEAD over even N ∈ [2,24]: degrees `3,5,7,9,11,11,15,15,17` at S2…S18 (**no
clean formula in N**), min weight at S14 `+0.01299` (positive), first refusal at **S20**.

The sharp part: `max(3, N−1)` is **right at S2, S12, S16, S18** and wrong at S4/S6/S8/S10/S14
— 4 of 9 buildable orders confirm it, **including S12, the order the retired frontier itself
made salient**. A spot-check drawn from the stale claim's own neighbourhood is biased toward
confirming it. (vv #13's congruence-class disguise, one level up: not a sampled group but a
sampled *parameter grid*.) I published this as its own ⚠ and pointed at the sweeping gate
`tests/numerics/test_advertised_degree_is_measured.py` (verified: it sweeps S2…S18).

⟹ And the fix for the drift is **not a better number** — it is a POINTER. The SSOT
(`docs/theory/methods/sn/angular_quadrature.rst` `quadrature-ls-positivity` +
`rules_sphere.py`) was already correct; `discrete_measures.rst` was the only page carrying the
stale copy as current. I replaced the numbers with a `:ref:` and said *why*: the frontier is
discovered by attempting the construction, so a second copy is exactly the thing that drifts.

### 6. Measured, not asserted — the evidence I generated before writing

Ran each before it appeared on the page: the 5 worked examples (all reproduce, incl.
`cylinder d=4 → LevelSymmetricSN(6)` fallback, 48 nodes); the stage-0/1 independence table
(4 rows, `Lebedev(5)` slab `✗/✓`, `GL(8)` cylinder `✗/✓`, `product(3,5)` cylinder `✓/✗`,
`product(3,6)` `✓/✓`); `spec.__doc__ is type(spec).__doc__` **True for all four** with one
shared `id()` and no instance `__doc__`; declared tags now honest (`σ_x`, `O_h`, `O_h`,
`D_{n_φ h}`, and `D_5h`/`D_6h` are **computed** by `_derived_product_group`, so citing them as
evidence about NODES is sound); the lattice-route failure `understated GL → lattice False /
nodes True`; `select_quadrature` has **no production consumer** (grep: def + export + 1
docstring + the test module).

### 7. Mechanics that worked

- Splice by line index with **guard asserts on live boundary strings** + structural asserts on
  the in-memory result **before** the write. First run FAILED on a bad assert (`\;` in a
  non-raw string) and the tree was untouched — the whole point.
- New content authored as a **pure literal** via `Write` to `scratch/`, so no Python string
  layer touched the LaTeX; removed after splicing (it would be a twin).
- An **enumerated list starting at `0.`** is legal: docutils sets `start="0"` and emits only
  an INFO (`report_level=1`), which Sphinx suppresses at default verbosity. Verified with a
  standalone `publish_doctree` probe before committing to the numbering that matches the code.
- ⚠ **Do NOT wrap quoted stale text in `*…*` when it contains `:math:`/`:eq:` roles** —
  docutils does not nest inline markup. Use the page's own idiom instead:
  ``⛔ X read :math:`...` until <date>``.
- Glyph check before use: `⛔` 12 files, `⚠` 17, `⭐` 10, `✓` 10, `✗` 2 in `docs/`. **The
  digest's "`⭐`/`⛔` have ZERO occurrences in docs/" is STALE** — they are corpus vocabulary
  now. Grep, don't recall.
- Two builds, not four (L-054 §9 held): baseline → all edits → all greps → xref gate → verify.

---

## L-057 — a new theorem lands ONE HOP from a page that already owns half of it; and its universal falsifies five sibling claims

**Task (2026-08-15, #344, branch `refactor/track-b-remainder`).** Document `LossKernelGauge`:
`A = L+C−S−B` is exactly singular on a `d ≥ 2` Cartesian diamond box with ≥2 reflective axis
pairs; a converged solve returns an arbitrary member of a solution manifold; a `G`-orthogonal
projector returns the canonical one. Brief named three candidate homes.

### 1. ⭐⭐ The home is decided by WHERE THE MECHANISM'S LOCAL HALF ALREADY LIVES, not by topic

The obvious homes were "boundary_conditions" (the BC is half the precondition) and "solver"
(the gauge fires at the exit). Both wrong. `cartesian_multid.rst` already carried
`.. math:: :label: dd-face-transmission-spectrum` — `Σ = (2/D)·1wᵀ − I`, its `−1` eigenvalue of
multiplicity `d−1`, the "undamped sawtooth", and the #340 measurement of that sawtooth's
signature at convergence — used there **only negatively**, as the obstruction to a Varga
comparison theorem. The kernel IS that local mode closed around a reflective loop. Filing it
anywhere else would have re-derived the transmission spectrum (a Cardinal-Rule-2 twin) and left
the existing section's reader with an incomplete story.

⟹ **before choosing a home, grep the corpus for the new result's LOCAL half.** If a page
already derives and labels the local fact, the global result is a downstream H1 on that page,
and its opening sentence should say *what the existing section stopped short of*. The module
docstring's own `:doc:` pointer named the wrong page; repointing it to the new `:ref:` was
part of the deliverable (a doc-only edit, AST-proved).

### 2. ⭐⭐ A new theorem is also a QUANTIFIER AUDIT — the doctrine it amends is asserted in N places

The result's headline ("a splitting shares a solution SET, not a POINT, when A is singular")
falsifies an unqualified universal that the corpus states **9 times** across 7 files. A
windowed regex (`fixed point` within ±2 lines of `invarian|same|shares`) found them; the
per-site adjudication was NOT uniform:

| site | verdict |
|---|---|
| `cartesian_multid` Key Facts "Only the fixed point is schedule-invariant" | present-tense FALSE → scope to **bulk** |
| `cartesian_multid` FP-invariance ¶ + *What survives* ¶ | FALSE as a universal → keep the sentence, add a `.. note::` tombstone (§3 discipline) |
| `solver.rst` Key Facts + "Two Inner Solvers" ¶ | FALSE → scope + tombstone |
| `foundations/cross_section_data` "identical to the plain (Jacobi) sweep" | FALSE → scope to bulk |
| `slab_one_group` Key Facts "(any consistent splitting … shares ψ*)" | true IN CHAPTER (d=1 kernel-free) but an explicit universal → one scope clause |
| `slab_one_group:852` (slab Krylov ≡ SI) | chapter-scoped, kernel-free → LEAVE, report |
| `foundations/boundary_conditions` C5.5 gate prose | ⭐ the gate's own FIXTURE is singular — see §3 |
| `loss_representation` scan-march gate | measurand is `scalar_flux` → name the measurand |
| `verification/sn.rst` T4 three-splittings | both configs kernel-free → say WHY it is legitimate |
| `foundations/boundary_conditions:2368/2433` (two source-DELIVERY channels) | ⛔ NOT a splitting claim — leave |

⟹ the tell that separates "leave" from "fix" is **which object the sentence quantifies over**:
a claim about *two delivery routes of one iteration* is untouched; a claim about *two splittings*
is not. And a chapter-scoped truth still needs a clause when the sentence contains an explicit
`any …` — a reader arrives at Key Facts by search, not by reading the chapter title.

### 3. ⭐⭐ Auditing a gate's PROSE, I found the gate's own fixture is in the new pathological class

`foundations/boundary_conditions.rst` credits a C5.5 "Mode-9 G-S ≡ Jacobi FP-invariance" gate on
a box described as *breaking every degenerate coincidence*: **x-reflective / y-vacuum /
z-reflective**, cells `(5,3,4)`. That is **two** reflective axis pairs. `[M]` I built it:
`dim ker A = 36` (`= n_g·(N/4)·n_y`). So G-S and Jacobi do **not** return the same trace there —
the gate survives only because what it asserts (`keff` + normalized flux *shape*) is mirror-even
and therefore blind by the new theorem. Publishing "the gate is sound BECAUSE its measurands are
mirror-even, and a future strengthening must not reach for the raw trace" is worth more than
either deleting the claim or leaving it.

⟹ **when a new result defines a pathological configuration class, MEASURE every gate fixture the
corpus names against the class predicate.** It is one `_as_sn_mesh(...)` + one attribute read,
and it converts a prose repair into a durable design constraint.

### 4. ⭐ The brief's headline rule was TRUE-on-its-fixtures and FALSE as stated

Brief: *"Excited iff the FIRST axis has an ODD cell count."* Reproduced exactly (5 excited rows
to 7 s.f., 6 inert). But two further measurements I ran unprompted:

- `dim ker A = 12` at `(2,2) (3,4) (4,4) (5,6) (6,8)` — **parity does not touch the kernel**;
- at even `n_x`, an **anisotropic** source `(1+μ_x)/W` excites it anyway: `1.756363e-02`
  (vs `6.7e-14` for the uniform isotropic source on the same mesh).

So the rule is about **excitation by a symmetric source**, never about the operator, and
"even `n_x` is safe" would have shipped as a mesh property. Published as a `.. warning::` with
both tables and the imperative *assert `dim ker == 0`, never infer kernel-freedom from a mesh
property*. This is `vv-principles` #13's congruence-class trap seen from the doc side.

### 5. ⭐ The "refuted witness gate" is the highest-value paragraph in the page

The campaign's FIRST acceptance-gate design (the `|Ω·n|⁰` full-face moment) **could not have
failed** — a face-summed moment is mirror-even and every kernel mode is mirror-odd, so the gate
was unfalsifiable by the campaign's own theorem, three sections above its own design. I
reproduced the null (`0.0` on 7 of 8 face rows, `~1e-16` on the rest, while the trace moves
6.08 %) and published the refutation beside the theorem that predicts it. A close-out's
*falsified-design* paragraph is what stops the next session re-deriving a green-forever gate.

### 6. ⛔ The changelog entry was BLOCKED by the page's own contract — report, do not fake

`history.rst`'s header states: *"Every entry below is merged to main … a new entry lands with
its merge hash or not at all."* `[M]` `git merge-base --is-ancestor f934ff57 main` → **NO**
(branch is 15 ahead). Every one of the 21 existing entries carries a merge hash. ⟹ writing the
entry now would mint exactly the class of falsehood the page's own 2026-07-24 repair retired.
Delivered the ready-to-paste row in the RETURN instead, with the reason.

⟹ **a deliverable can be blocked by a page's own stated contract; that is a finding, not a
failure to deliver.** Check a changelog's header rules before writing to it.

### 7. Measured, not asserted — everything on the page came from my own probes

Parity table (11 meshes) · Jacobi + Krylov controls · tangential vs normal currents on 2
quadratures · 8 mirror-even functionals before/after · `dim ker A` on 12 configurations
(incl. graded, mixed-BC at two different vacuum axes, LD, d=1) · closure registry
damping × ndim · T-component `G ≡ 0` on 3 rules · projector idempotence + uniform-trace
annihilation (`5.0e-18`) · build/apply cost · **the strongest single row**: gauged trace
recovers the analytic flat answer `6.09e-02 → 1.04e-13` with `‖Π(t−t_exact)‖/‖·‖ = 1.00000000`.

⚠ Two brief numbers I could NOT reproduce as stated and therefore did not publish: the plan's
`‖t_SI − t_Krylov‖ = 1.320828` (a *heterogeneous* fixture I did not have; mine reads `0.124184`
on the homogeneous one — published with MY configuration), and the memo's `41.1 ms` build
(mine: `22.0 ms` at d=3, the fused-SVD implementation). **Publish your own number with your own
configuration; never relay one whose fixture you cannot state.**

### 8. Mechanics

- Fragment authored as a **pure literal** via `Write` to `/tmp`, then spliced by a script with
  (a) underline-length + title-level-skip checks, (b) a `list-table` `:widths:`-vs-cell-count
  checker, (c) an odd-backtick scan — **all run on the fragment before any write**. The
  `:widths:` checker caught a real 3-cell row in a 2-column table.
- ⚠ `Edit` failed on a paragraph whose text I had just read: the source carried a **typographic
  apostrophe** `’` (U+2019) in *"splitting's"*. Shorten the anchor or `repr()` the live line.
- Two builds only (baseline `-E -W` → all edits → all greps → xref gate → AST doc-only proof →
  one verify build). Both EXIT=0 with **0** WARNING/ERROR/CRITICAL — the set unchanged.
- Post-build link audit is what proved the refs: **9 of 9** cross-doc `:ref:`s render as real
  `href`s, the cross-doc `:eq:` resolves, and 6 of 9 code-xrefs link. The 3 that render plain
  text (`solve_sn`, `SNMesh.loss_kernel_gauge`, `IterationHistory.*`) are the page's
  `:noindex:`-automodule convention — `[M]` **0** hrefs and **0** anchors tree-wide for those
  targets *before* my edit, so it is not a regression and half-surfacing one module is forbidden.
- V&V matrix auto-regenerated `534 → 539` documented sentinels (+5, exactly my labels); orphan
  count unchanged at 2. Never hand-edited.

---

## L-058 — a path census keyed to ONE artefact's FILENAME misses its SIBLINGS in the same directory

**Task** (2026-08-15, branch `chore/nexus-project-config`): repoint four instruction surfaces after
the Nexus graph's location moved from a hardcoded `docs/_build/html/_nexus/graph.db` to a single
declaration in `.nexus/config.toml` (`[graph]`), resolving in ORPHEUS to
`docs/_build/html/graph/graph.db`. Brief was explicit and correct about its own limit: *"I censused
with `grep -rln "graph\.db"` which only catches that exact spelling — if the concept is stated some
other way, I want to know."*

### The finding the brief's census could not reach

A build directory holds an artefact **family**, not one file. Measured after a fresh `-E` build:
`docs/_build/html/graph/` contains `graph.db`, `graph.json`, `graph.html`, `traces/` — and
`docs/_build/html/_nexus/` **does not exist at all**.

`docs/index.rst` carried, under a `Knowledge Graph` heading:

```rst
`Open interactive graph explorer <_nexus/graph.html>`_
```

`[M]` the shipped homepage rendered `href="_nexus/graph.html"` → **404**, while
`docs/_build/html/graph/graph.html` (627 053 bytes) sat un-linked. A dead link on the docs
homepage, invisible to every census keyed to `graph.db`, and invisible to every build:
**a raw relative hyperlink is checked by Sphinx at NO severity** (unlike `:doc:`/`:ref:`, which
warn). So `-W`, `-n` and `check_docstring_xrefs.py` are all silent — the last one by design, it
gates Python-domain roles.

⟹ **When a DIRECTORY moves, census the directory, not the file.** Grep the parent segment
(`_nexus/`) and each sibling extension (`graph\.(db|json|html)`), not just the one filename the
brief names. One extra alternation caught a user-facing 404 that four rounds of `graph.db` grepping
could not.

⚠ **Grep hygiene that cost two wasted rounds:** `_nexus` as a bare pattern matches `mcp__nexus__*`
— every MCP tool name in every agent/rule/skill file, 559 KB of output. Anchor it as `_nexus/`
(with the slash) and the census collapses to 9 real lines.

### The residual, stated rather than hidden

A static RST hyperlink has no mechanism to ASK where the graph lives — it cannot run
`nexus config db`. So the repair necessarily mirrors `[graph].output`, i.e. it IS a second
declaration. Rather than pretend otherwise, the fix ships an RST comment above the link naming the
coupling and the reason nothing verifies it. A second declaration that **announces itself as a
mirror** is the honest floor when the single-source mechanism is unavailable; a silent one is the
defect.

### `--db` optionality — the fix is DELETION, not substitution

`[M]` `--db` is optional on all 10 subcommands probed; `resolve_db(explicit, start)`
(`sphinxcontrib/nexus/project.py:88`) encodes **flag > `[graph].db` in the nearest
`.nexus/config.toml` > `LEGACY_DB = _nexus/graph.db`**. So in documentation the correct edit for
most examples is to **delete the flag**, not to update its value — updating it is what mints the
next stale literal. `nexus-exploring/reference.md` taught `--db <path>` on **16** command lines;
all 16 came out, replaced by one header sentence stating the precedence.

⭐ **But keep the flag where naming a file is the POINT of the example.** Per-example judgement,
not find-and-replace: `nexus analyze . --db /tmp/scratch-graph.db` is now *better* documentation
than before, because with a default in place the flag finally has a teachable meaning — a
deliberate override (a scratch graph, a second checkout, a graph you are diffing against).

### The asymmetry a config-driven CLI introduces, which no doc stated

`[M]` `--project-root` exists on `analyze`/`serve`/`config`/`file-brief`/`staleness`/`retest`/
`changes`/`rename`/`briefing`/`audit`, and **`status` REJECTS it**
(`nexus: error: unrecognized arguments: --project-root`). Since `resolve_db` walks up from
`--project-root` *or cwd*, the read-only query subcommands are **cwd-anchored**. Run `nexus status`
from a scratch directory and you get `Error: _nexus/graph.db does not exist / Run 'nexus analyze'
or 'sphinx-build' first` — a message that reads as *"the graph was never built"* when the truth is
*"you are standing in the wrong directory"*. Both skills now state it.

### My own guard asserted the wrong thing — and that is the system working

The splice guard `assert len(out) < len(src)` fired red. The content was fine; the *guard* was
wrong (a 5-line header replacing a 1-line one outweighs 16 × 12 removed characters). The file was
**untouched** — which is exactly the §5 discipline's payoff: every structural assert runs on the
in-memory result BEFORE any write, so a false red costs one re-run and never a `git checkout`.
⭐ The transferable half is `vv-principles` #4's VERIFY sharpening turned on my own instrument:
**a failed check is not a finding until you have diagnosed WHOSE failure it is.** An earlier assert
in the same script had already caught a real miscount (16 flags, not the 15 I eyeballed) — so the
instrument had a positive control before it produced its false red.

### Verified-not-assumed, and left alone

The brief predicted `AGENT.md:104` (`_build/html/ ← Build output (includes Nexus graph.db)`) was
still true and said to verify rather than assume. `[M]` `find docs/_build/html -name graph.db` →
`docs/_build/html/graph/graph.db`, i.e. still *under* `_build/html/`, and the line names no
subpath. TRUE → left. `AGENT.md:508` names `graph.db` bare with no path, and its claim (no graph in
a fresh worktree until the first `-E` build) still holds — `.nexus/config.toml` is tracked, so the
path resolves fine while the file does not exist. TRUE → left. **Two lines correctly NOT edited is
a deliverable**, and saying so is what stops the next session re-opening them.

### Gate

`-E -W --keep-going` baseline **re-measured this session**: `0` WARNING/ERROR/CRITICAL, EXIT=0.
Post-edit identical: `0`, EXIT=0. Built-HTML proof rather than source proof:
`href="graph/graph.html"` with the target present, and `grep -c '_nexus'` = **0** on both
`index.html` and `development.html`.

---

## L-059 — a machine-readable DECLARATION is a doc surface, and its blast radius is the CLAIM it displaces

**Task (2026-08-17, nexus #82).** Author `.. implements::` declarations for one theory page
(`docs/theory/methods/sn/loss_representation.rst`), because a nexus change made *declaring any
implementer of an equation stand the guessing down for that whole equation*. Plus a fix-on-sight
for a falsified docstring. Landed 28 directives over 14 equations, a new H2 recording the three
equations that deliberately declare nothing, and 19 code-docstring corrections.

### 1. ⭐⭐ Measure the instrument you are replacing — the number is the section's whole argument

The brief supplied the mechanism ("all 14004 `implements` edges are guesses"). It did not supply
what the guesses *were*. Four SQL queries against `.nexus/graph.db` turned a policy statement into
the page's load-bearing content:

| measured, pre-declaration, over the page's 14 declared equations | value |
|---|---|
| inferred `implements` edges the heuristic wrote | **397** |
| true implementers (what I then declared) | **28** |
| of those, found by the heuristic | **6** (21 % recall) |
| precision | **1.5 %** |
| guess sets for `loss-rep-LpC` vs `loss-rep-facewise-separable` | **identical, 23 for 23** |
| of `loss-rep-LpC`'s 23 guesses, how many are its 2 real implementers | **0** |

The last row is the one that explains the mechanism rather than scoring it. Both real implementers
live in `orpheus.sn.operators.streaming`; the shared token is the *package name*
`loss`/`representation`, so the guess set is exactly "the membership list of
`orpheus.sn.loss_representation`" — a set that **cannot contain them by construction**. So the
failure is not "imprecise"; it is *not a claim about the equation at all*.
⟹ **When you replace an instrument, publish its measured error, not its described one.** One query
per equation; the table writes the section.

### 2. ⭐⭐ Writing the explanation MINTED new guesses — the honest edit made the metric worse

Post-build the three UNDECLARED equations' guess counts had *risen* (23→24, 23→25, 23→24).
Cause, diffed: my own prose. Explaining why `MaterialXSField.foldable_sigma` is **not** the
implementer of `loss-rep-removal-sigma` added a `:meth:` xref to it — and `foldable_sigma` shares
the token `sigma`, so it was promptly inferred as a new implementer *of that equation*. Same for
`LossRepresentation.streaming_action` on all three.
⟹ **Citing a symbol in order to say it is NOT the implementer is enough to make the heuristic name
it as one.** Two consequences: (a) NEVER publish a live guess count — quote the frozen
pre-declaration measurement or tell the reader to re-run; (b) this is a real finding for the tool,
not a curiosity: an equation with no declaration gets *worse* every time its page is improved.

### 3. ⭐ The authoring contract inverts the usual risk: INCOMPLETE is worse than ABSENT

Because declaring stands the guessing down *per equation, not per pair*, an equation declared with
one of its two implementers shows **only the one** — the guess that used to cover the second is
switched off. So the failure mode of this doc surface is silent under-coverage produced by an act
that looks like an improvement. Discipline: for every equation, ask *what else computes this?*
before writing the first directive. Seven of the fourteen needed 2–4 directives (DD arm + LD arm;
forward + transpose; scheme door + the schedule that folds its transverse term).

### 4. ⭐⭐ The brief's site census was 6 of 18 — a windowed CONCEPT grep found the rest

Deliverable 2 named `streaming.py:546-548` plus six "may be inherited" sites. A windowed regex
(`subtract` within ±4 lines of `Resolution A|StreamingOperator.apply|operator subtract`) found
**18 sites in one file**, all asserting the same present-tense falsehood — that
`StreamingOperator.apply` *subtracts* σ. It does not: since #257 S8b its whole body is
`streaming_action(psi)` = `loss_action(0, psi)`. Leaving 12 copies of the sentence I was fixing
would have been the exact half-fix vv #21 warns about. Fixed all 18 + the brief's one, and
REPORTED the expansion.
⚠ The tell that the file already knew: at `__init__.py:371-375` the *corrected* framing
("`StreamingOperator.apply` calls this directly (#257 S8b) so L reads no σ") sits **one line above**
the stale sibling docstring at :376. A correction pass that stops at the method it was looking at
leaves the file citable for either sentence.

### 5. ⭐ The same falsehood was in the RST — including the page's own Key Facts card

The brief scoped Deliverable 2 to code. The theory page carried the identical claim in **four**
places, one of them the **Key Facts** admonition ("*Gotcha — the operator subtracts C once*"), and
another the prose wrapping `loss-rep-resolution-a` — the very equation I was adding four
`.. implements::` blocks to. Not fixing it would have put my new prose three screens from its own
contradiction, with me as the second voice.
Fix shape, all four: keep the **equation** (the identity `Lψ = (L+C)ψ − σ_t⊙ψ` is TRUE), keep the
**label** (two live `verifies()` markers point at it), correct only the *attribution*, and tombstone
the history in the past tense with a `.. note::` naming the two carves (#240 Step B removed the leaf
sum; #257 S8b removed the subtraction). Retitled the section from "the operator's only glue" to
"one action, two readings of σ" — safe, because the section carries no `.. _anchor:` and the only
tree hits were orphaned `_build/` HTML.

### 6. ⭐ A quantifier in a page table needs its own census, not the brief's

The brief's ⛔ ruling said "the only sites forming `σ_t − σ_s0` are `derivations/.../dsa.py:632`/`:1023`".
`[M]` there are **four**, and the two the brief missed are the more interesting ones:
`orpheus/sn/acceleration/dsa.py:328` (`quarter = 0.25 * (sigma_t - sigma_s0) * h  # ¼ σ̂_R h`) is
**production**, not a derivations mirror; and `orpheus/derivations/continuous/mms/sn.py:1892`
(`SigC=np.array([sigma_t - sigma_s0])`) is a **capture cross-section** — identical arithmetic, a
*material datum* rather than an operator diagonal, numerically coincident with σ_r only because the
fixture has one group. The brief's grep keyed on `sigma_r =`/`sig_r =`, i.e. on the NAME; two of
four sites never bind the name. The ruling survived and got stronger — but the count did not.
⚠ Beware the mirror trap: `sig_r` in `thermal_hydraulics/` and `kinetics/` is a **radial stress**.
A short suffix collides as badly as a one-letter symbol.

### 7. ⭐ "Implemented by nothing" is a CLASSIFICATION, and it is the durable half

Three labels declare nothing, for three *different* reasons, and naming the kinds is what makes the
section reusable: **superseded path** (`loss-rep-leaf-sum` — two independent retirements, the #240
override and the #257 σ-removal, each alone sufficient to make the route unreachable; the identity
stays true, the code that walked it is gone) · **notation** (`loss-rep-removal-sigma` — a definition
with no production caller; every site that computes the arithmetic computes a *different operator*)
· **declared tag** (`loss-rep-facewise-separable` — a `ClassVar[bool]` a human set after doing the
tensor-product argument by hand; the implementer of the *criterion* is the page).
⟹ *A statement can be true, labelled, and implemented by nothing* — and the cases where that is
correct are enumerable, which is precisely what an inference cannot know.

### 8. Gates (all re-measured this session)

* `-E -W --keep-going` baseline **0** W/E/C, EXIT=0 → post-edit **0**, EXIT=0.
* `tools/check_docstring_xrefs.py orpheus tests docs`: HEAD baseline (via `git archive` into a temp
  tree) **81 dead / 124 sites**; post-edit **81 / 124** — identical, while adding **80** xref roles.
  Measuring the baseline from `git archive HEAD` is the cheap way to get a true before/after on a
  dirty tree.
* AST doc-only proof vs `HEAD` for both edited `.py` files (docstrings stripped, `ast.dump` compared)
  — re-run after the *last* edit, not the first.
* Graph confirmation is the real acceptance test for this deliverable: query
  `edges WHERE type='implements'` per equation after the rebuild. 397 → **28**, each equation's count
  equal to the multiplicity declared. A directive that fails to resolve DOES warn
  (`target %r not found in graph — skipping`), so `-W` gates the `:by:` paths — but only the graph
  query proves the *count*.
* Pre-flighted every `:by:` path and every equation label against `graph.db` **before** writing a
  line (23 targets, 17 labels, all resolved). Cheaper than a 6-minute build per typo.

### 9. Reported, not fixed — `tests/` is not mine

Three test-module docstrings assert the retired mechanism. All three were left alone and reported:

* `tests/sn/operators/test_loss_action_convention.py:3-9,20-22` — "*the operator's `apply` applies
  the ONLY algebra glue, the Resolution-A collision subtraction*" and "*`apply` is DEFINED as
  `loss_action − σ_t·ψ`*". ⭐ Its own **function**-level docstring (`:133,:141`) is already correct
  ("the **+C glue**", "the affine relation") — the module header lags the body it introduces.
* `tests/sn/operators/test_streaming_operator.py:8-19` — a `:=` **definition** section titled
  "Resolution A — subtractive definition", plus "*L carries σ_t at constructor time*".
  `[M]` `StreamingOperator` is a dataclass with **one** field (`sn_mesh`) and no `sigma_t` attribute.
  A docstring asserting a constructor signature that does not exist is the loudest class of stale.
* `tests/sn/sweep/core/test_phase_c_gates.py:22,25,371` — names `:class:`CollisionOperator``
  (`[M]` retired at #261, importable from nowhere) and attributes the composite matvec to the leaf
  sum. Its *conclusion* `(L+C).apply(ψ) = M(ψ;σ_t)` is TRUE; only the mechanism is stale.

---

## L-060 — a SPEC has a headline and a table, and they keep different clocks; plus: a node that EXISTS can still not RESOLVE

**Task (2026-08-17, nexus #82, sibling of [[L-059]]).** Author `.. implements::` declarations for
`docs/theory/foundations/operator_algebra.rst` from an explorer-written spec, record the
no-implementer taxonomy, and repair four drift findings the exercise exposed. Landed **57**
directives over **32** equations, a new H1 recording the contract + the 8 un-declarable equations,
and repairs A–E. `-E -W` EXIT 0 / 0 W-E-C, unchanged from a freshly measured baseline of 0.

### 1. ⭐⭐ Count the spec's TABLE; its headline is a summary and summaries rot

The spec file's own headline read **"21 of 40 declarable, 19 NONE"**, and the dispatching brief
inherited it verbatim ("The 19 NONE equations…", with a kind breakdown of 4+1+1+2 = **8**, which
does not sum to 19). One `re.findall` over the table rows: **40 rows, 32 declarable, 8 NONE, 57
implementer slots over 55 distinct symbols.** The table was right and internally consistent — its
`§3` kind taxonomy names exactly the 8 — and only the headline was wrong.

⟹ **A spec is its table.** Before designing to a brief's counts, re-derive them from the artefact
in one command; a wrong headline propagates into the brief, into the section you write, and into
the return. Had I written "19 NONE" into the page it would have been a published falsehood with an
authoritative-looking provenance chain (explorer → main agent → me).

### 2. ⭐⭐ "All N node IDs resolve" ≠ "all N `:by:` targets bind" — the DIRECTIVE's resolver is narrower than the graph

The spec certified *"every dotted path existence-checked against `.nexus/graph.db` — 55 of 55 node
IDs resolve"*. True. Two of them still would not have bound.

`_node_id_for_target` (`sphinxcontrib/nexus/directives.py`) tries the literal string, then
`py:function:`, `py:method:`, `py:class:` — **and nothing else**. `Domain` / `Codomain` are
`TypeVar`s, so their nodes are `py:data:orpheus.numerics.operator.Domain`. A bare
`:by: orpheus.numerics.operator.Domain` logs *"target not found in graph — skipping"* and lands
nothing — which under `-W` is a red build, and without `-W` is a silent no-op.

Fix: the directive **accepts an already-prefixed node id** (`if target in graph: return target`),
so `:by: py:data:orpheus.numerics.operator.Domain` binds. Measured post-build: 57 directive edges,
0 skipped.

⟹ Pre-flight `:by:` paths **through the resolver's own prefix list**, not by asking "does a node
with this name exist". [[L-059]]'s graph pre-flight would have passed both. And write the reason
into the page — a future author copying the two `py:data:` lines needs to know why they differ.

### 3. ⭐⭐ A brief's "sharpest observation" is a HYPOTHESIS with a computable confusion matrix

Both the spec and the brief pressed one finding: *the page already labels its own two classes* —
every NONE row's `.. (vv-status rationale)` says "Mathematical identity" / "Definitional", while
every declarable row's names a verb, a value, or a file. The brief even pre-corrected one half of
it. It is still not a classifier, and the check is four greps:

| over the audited 40 | NONE (8) | declarable (32) |
|---|---|---|
| carries a rationale block at all | **6** | **22** |
| rationale contains *identity* | 5 of 6 | **11 of 22** |
| contains *"not a solver claim"* | 1 of 6 | **5 of 22** ← points the wrong way |
| cites a `tests/` file (the "declarable" signal) | **2 of 6**, incl. `operator-solve` | 13 of 22 |

A third of the page carries no rationale at all, and the word *identity* appears in half the
**declarable** rows. The cause is a real ambiguity, and naming it is the finding that survives:
an **identity between quantities** (`apply-solve-parallel-identity`) can have no carrier, while an
**identity between types** (`carrier-grid-operator-typing`, `harmonic-frame-is-galerkin`,
`product-solve-reroute`) is *exactly* a claim about a class declaration. Both are honestly called
identities.

⟹ Publish the measured split (`{identity, law, canonical-form} → NONE`;
`{typing-rule, definition} → look for a declaration site`) and publish the **refutation of the
keyword tell** beside it, because the next reader will otherwise re-derive the keyword heuristic
and ship it.

### 4. ⭐⭐ Before repairing a stale equation, ask whether the corpus already states it CORRECTLY

Finding A (`keff-as-integrated-rates` present-tense-false: `(n,2n)` in the numerator, no leakage)
was real and independently confirmed against `SNSolver.compute_keff`. But a 3-command corpus
census after the repair found `docs/theory/methods/sn/solver.rst` already carrying the shipped form
as `:eq:`sn-keff-update`` under `:ref:`sn-keff-estimator`` — with the divergence-telescoping
derivation, the leakage functional, and the wiring to the cross-engine gate. My repair had just
minted a **twin**.

Correct shape, applied: keep the equation (a labelled equation is an API and must not state a
falsehood) + an `.. important::` naming the SSOT, saying what *this* page's claim actually is (the
**typing** claim: both ends of the ratio are the same typed functional), and instructing that
future changes are edited there and mirrored here. Same census on finding C was *evidence the
repair was right*: `frame.rst`, `sn/slab_multigroup.rst`, `sn/cartesian_multid.rst` and the class
docstring **all already spelled it `Λ`** — `operator_algebra.rst` was the sole holdout writing `S`,
contradicting its own rejection note 900 lines up.

⟹ The census is cheap and it does two jobs: it stops the repair becoming a twin, and it tells you
whether you are fixing an outlier or inventing a convention. Run it **before** drafting, not after.
⚠ It also surfaced a symbol collision: the SSOT writes leakage `L`, this page writes
`L = Ω·∇` everywhere. Resolved by `L_{\rm leak}` **plus a note naming both spellings** — never
silently.

### 5. ⭐ My own uniqueness guards were substring bugs — twice, and labels are PREFIXES of each other

Two splice guards fired before any write (so cost nothing, per the assert-before-write discipline):

* `result.count(":label: operator-apply") == 1` — **fails**, because `:label: operator-apply` is a
  substring of `:label: operator-apply-transpose`.
* `result.count("Development history") == 1` — **fails**, because an `.. admonition:: Development
  history — G6.3 step 8.0` sits 1000 lines above the section.

⟹ A uniqueness guard over labels/titles compares **exact lines**
(`sum(1 for l in lines if l.strip() == …)`), never substrings. Eq-label families are built by
suffixing (`X`, `X-transpose`, `X-section`), so the prefix collision is the *normal* case here, not
an edge case. And diagnose a red guard as possibly the GUARD's error first — both of these were.

### 6. ⭐⭐ The plan-authoring quantifier clause applies to the prose YOU publish — I broke it twice in one section

Both caught by re-measuring my own sentences before the final build:

* *"…:eq:`scattering-carrier-grid`, :eq:`scattering-aniso-composite` and
  :eq:`reaction-rate-kinf-oracle` with three each"* reads as an enumeration of the
  multi-implementer equations. There are **15**, and **five** have three. Replaced with the count.
* *"60 on :eq:`operator-solve` alone, where every ``solve`` in the tree matches the label by name"*.
  Measured the 60 sources: only **5** are named `solve`; the rest are five `apply` methods, three
  `solve_fixed_source`, `is_invertible`, `is_adjointable`, `outer`, and whole classes. The match is
  on the **other** token. The corrected sentence is strictly more damning *and* true.

⟹ Every "each / every / all" in a doc paragraph is a universal over a set you can count in one
command. Count it. The measured version is almost always the better sentence.

### 7. The mechanism verified itself mid-session, by accident

A concurrent archivist's Sphinx build rebuilt `.nexus/graph.db` **between** my repair pass (which
carried 6 declarations) and my declaration pass. That snapshot therefore held both sides of the
comparison in one graph: the **3 declared** equations carried **0** inferred edges; the **37
undeclared** ones carried **771** (median 11, max 58). Post-pass: **57 directive, 0 inferred** on
the 32 declared; **166 inferred** remaining, every one on the 8 that cannot be declared.

⟹ When an accident hands you a clean before/after in one artefact, **say so in the page** — a
measurement that carries its own control is worth more than two dated ones. And note the corollary
the page now records: an equation that legitimately has NO implementer keeps its guesses forever,
because the stand-down is triggered by a declaration and **there is no way to declare an absence**.

### 8. Placement rule for a directive whose body RENDERS

`.. implements::` with a body emits a `<div class="docutils container">` — visible prose, no visual
marker. 57 of them dropped "immediately after the `.. math::`" would land mid-sentence wherever the
equation is followed by its own where-list. Rule adopted and stated in the splice script: **after
the math block, unless the next paragraph is a grammatical continuation of the equation's sentence
(`where …`, `so …`, `with …`, `and identically for …`) — then after that paragraph.** Encoded as a
per-label `skip` flag, previewed before writing. Every body opens `**Implemented by** …` so it
reads as an annotation rather than as body text.

### 9. Reported, not fixed

* The spec's headline count (§1) — flagged in the return, spec file not edited (not mine).
* `nexus` `_node_id_for_target` should try `py:data:` (§2) — a real tool gap; a bare TypeVar name
  is the natural thing an author writes.
* Eight further labelled equations on the page were outside the audit's scope
  (`carrier-grid-interchange-witness`, `tensor-product-axis-wise-composition`,
  `sum-of-tensor-products`, `octant-direct-sum-tensor-product`, and the four `eigen-*`).
  `[M]` all eight attract **zero** edges of either kind — unfinished, not wrong. Recorded in the
  page's own coverage subsection so the next pass finds them.

---

## L-061 — a mechanical port's WARNING COUNT is a non-representative sample of its DEFECT COUNT

**Task.** Clear 20 Sphinx warnings in `docs/theory/verification/error_catalog.rst` — a 5790-line
RST port of the 79-entry L0 error catalogue from `.claude/skills/vv-principles/error_catalog.md`,
done by a throwaway script. Brief: "It handled the bulk correctly… 20 warnings remain, and they
are genuine per-entry judgement calls." Result: 20 → 0, `EXIT=0`, 79 entries / 258 catchers intact,
and the xref gate's `error_catalog.rst` rows gone (75 → 71 dead sites tree-wide).

### The premise the measurement refuted

"The bulk is correct" was false, and one command showed it. **In RST there is no legal run of 3+
backticks outside a literal block**, so a run-length histogram is a total census:

```
RST: {1: 186, 2: 4010, 3: 678, 4: 152}      MD: {1: 4332, 2: 846, 3: 46}
```

**830 mangled delimiters on 339 lines, zero inside code blocks.** The 20 warnings were the ~2 %
of that class where the imbalance failed to cancel *within a paragraph*. Rendered HTML proved
it visible: `<code>`psi_right = fi[:, n, i, 0]``</code>`, `<code>`de8822d`</code>`.

⟹ **Before fixing warning #1 of a port, census the delimiter alphabet of the target language.**
The warning count measures where a *parser* choked; it does not measure where the *render* is wrong.

### One root cause, three surface families

The script's `` `x` `` → ``` ``x`` ``` regex was **LINE-LOCAL**, and a code span that WRAPS a line
defeats it three different ways:

| MD form | what the script did | symptom | count |
|---|---|---|---|
| `` ``x`` `` on one line | added a pair → 3–4 backtick run | mostly silent; stray backticks | 830 runs / 339 lines |
| `` `x` `` wrapping a line | converted ONE side → 1-vs-2 | **warns**, or cancels silently | 14 spans |
| `` `x` `` wrapping a line | converted NEITHER side | silent `<cite>` (italic, not code) | 16 spans |

`grep -c '<cite>' built.html` is the census for the third — **`default_role` is unset in this
project**, so every surviving single-backtick span renders *italic* instead of monospace. It is
the smoking gun for any Markdown→RST port, and it is invisible at every build severity.

### The port's own SOURCE is the oracle — it turns a 415-site blanket edit into a proof

Normalise, then check every restored literal's CONTENT against the Markdown:

```
inline literals: 2443   content not found verbatim in the .md: 5
   → 2 authored by the port itself (its new header note), 3 the `\|` class below
```

Same for prose, filtering the expected transformations (MD `#` headers → `:title:`, MD `|` tables
→ `list-table`, the replaced preamble): **3648 of 3653 MD prose lines ≥45 chars survive verbatim**;
the 5 exceptions were 1 intentional repoint, 3 artefacts of my own per-line de-markup on wrapped
`:math:` roles, and 1 correct `[x](url)` → `` `x <url>`_ `` conversion.

⟹ **A bulk delimiter edit is guarded by `src.replace('`','') == new.replace('`','')`** — proves
only backticks moved — plus an exact character-count delta and an unchanged line count. Cheap,
total, and it converts "risky mass edit" into "verified transformation". Do the write only after
the guards pass, so a failed assert leaves the tree untouched (no `git checkout` recovery needed).

### ⭐⭐ PROBE docutils, do not reason about it

I predicted the ERR-079 warning came from *emphasis containing an inline literal*. **Wrong.** One
`publish_doctree` probe with 6 one-line cases settled three entries at once:

| construct | docutils |
|---|---|
| `*"… ``lit`` …"*` (emph ⊃ literal) | **0 warnings** — and renders RAW backticks |
| `*"… **strong** …"*` (emph ⊃ strong) | **WARNS** ← the actual culprit |
| `key=``"zero"``` | **WARNS** (`=` forbidden before inline markup) |
| `key=\ ``"zero"``` | clean, identical render |
| `γ_-` in prose | `ERROR: Unknown target name: "γ"` |
| `γ\_-` | clean, identical render |
| `1. text:` then `   - a` (no blank line) | `Unexpected indentation` + `Block quote ends…` |
| `(Wave 2)\n+ the typed error` (no blank line above) | clean — `+` mid-paragraph is not a bullet |

A stub-directive/stub-role `publish_doctree` harness (register `error-entry` etc. as pass-throughs)
re-checks a 5790-line file in **under a second** vs a ~4-minute `-E` build. Build twice, iterate
in docutils.

### The Markdown discriminator for an indented block

CommonMark: **an indented code block cannot interrupt a paragraph.** So the fix differs by whether
a blank line precedes:

- **blank line before** ⟹ genuine code block ⟹ `.. code-block:: text` (mandatory when the body
  contains `*` — ERR-023's `w *= 2.0` was read as emphasis).
- **no blank line** ⟹ lazy paragraph continuation ⟹ blank lines around it → **block quote**, which
  is also what the port's other 14 indented blocks became, so it is the consistent choice.

### RST forbids inline markup after most characters — a port hits this constantly

Openers must follow whitespace or one of ``- : / ' " < ( [ {``. Markdown has no such rule, so the
port left 9 literals and **2 `:math:` roles** opening after `=`, `.`, `~`, `§`, `↔`, `*`. The
roles **were not rendering at all** — `~:math:`\mathcal{O}(h^{1.3})`` produced
`<cite>mathcal{O}(h^{1.3})</cite>` (role dead, LaTeX backslash silently eaten) and **no build at
any severity said so**. Fix is one character: `~\ :math:`…``. Tell = `<cite>` in built HTML.

### `\|` is right in prose and WRONG inside a literal

The script escaped 37 pipes; the MD had 0. In prose `\|` renders `|` (fine). Inside `` `` `` RST
does not process escapes, so it renders a **visible backslash** — measured in the HTML. Exactly 3
sites; the other 34 are harmless. Discriminate by context, don't blanket-revert.

### ⭐ Adjudicating dead `:mod:` in a HISTORICAL narrative — tense AND object survival

Four dead `orpheus.sn.spatial.*` targets, all inside ERR-026's "What Wave H Phase A/B added"
narrative. `git log --diff-filter=D` split them, and the answer is **not** uniform even though
three of four modules survived a pure `git mv`:

| site | sentence | module fate | named object fate | verdict |
|---|---|---|---|---|
| `boundary_face_flux` | "What Phase A added … Protocol (X)" | **DELETED** `3fd1302f` | Protocol retired | ``literal`` |
| `pole_angular_closure` | "What Phase B added … Protocol (X) with three strategies" | renamed `588f2429` | Protocol + 2 of 3 strategies retired (#248) | ``literal`` |
| `pole_angular_closure` | "**Documented in** A, X, and Y" | renamed | claim still true there | **repoint** |
| `diamond` | "**Citations updated in** A, B, X" | moved `5b6598f0` → `transport/spatial` | corrected BMC 2010 citation IS at `:51` | **repoint** |

⟹ **A surviving module does not license a repoint; a surviving CLAIM does.** Row 2 is the trap —
same file, same rename, opposite verdict from row 3, because the *sentence* names objects that no
longer exist there. Three corroborations, all free: the live tree's own prose spells the retired
names as ``literals`` (`sn/sweep/__init__.py:37`, `pole_angular_closure.py:93-95`); the SAME
catalogue entry already spelled the deleted path as a literal **130 lines below** the site; and a
list of three `:mod:` roles where two are live argues against making the third a literal.

### ⭐ A dead `:doc:` from a Markdown port is a PATH-FORM error, not a missing page

`:doc:`docs/theory/methods/sn/index`` — the page **exists**. MD authors write repo-root paths;
Sphinx wants a docname. Fix `/theory/methods/sn/index`, don't rewrite the prose. **Check the page
exists before concluding the reference is dead** — the brief and the warning both read as
"pointing at nothing".

### What I deliberately did NOT fix, and why

`[M]` **32 rendered `<strong>`/`<em>` elements contain raw `` `` ``** (86 raw pairs in page text
outside `<code>`) — Markdown bold/italic *containing* a code span, which RST cannot nest. Zero
warnings. Unlike the delimiter class (pure arithmetic, provably content-identical), each repair
must **choose where to break the emphasis run**, i.e. exactly the per-site judgement the brief
scoped to 20. Reported with the line list instead. Post-fix the class is *clean* `` `` `` pairs,
so a follow-up pass is mechanical.

⟹ **The scope line that held: fix what is provably content-identical, report what needs a
choice.** State the expansion loudly (830 + 30 + 3 + 2 sites vs a 20-warning brief) and give the
measurement that forced it.

---

## L-062 — a cross-reference inside HISTORY is a category error, and BOTH gates are blind to it

**Task (2026-08-18, branch `docs/err026-history-is-not-a-crossref`).** `docs/theory/verification/
error_catalog.rst` ERR-026 carried **29 python-domain roles, 20 unique, 15 of them dead**, all in
one 154-line block of Wave-E/Wave-H project archaeology. Ruling applied: **an ERR entry's body is
past-tense archaeology; a role is a present-tense claim that the symbol exists NOW at THAT path.
The two cannot be combined** — the catalogue exists *because* the code moved on, so a role inside a
historical narrative is guaranteed to go false. 29 → 13 roles, 0 unresolved; `-E -W` EXIT 0, W/E/C
count 0 = 0 baseline.

### 1. Why nothing caught it — and the ONE-LINE fix, measured

Two instruments, both silent, for two *different* reasons:

* **nexus dead-references** judged only 3 of 15 — the 10 bare roles (``:class:`SNStreamingOperator```)
  are "undecidable" and filtered out.
* **`tools/check_docstring_xrefs.py`** — my digest calls this "THE gate". `[M]` it reported the same
  **81 dead / 124 sites** before AND after. It is not that it judged them alive: **it DECLINED all
  15**, and 3 of those it could have decided.

`[M]` the mechanism, by direct call with `namespaces=()` (what an `.rst` page has — the project
carries zero `currentmodule`):

| target | role | `resolve()` | `judge()` |
|---|---|---|---|
| `orpheus.geometry.boundary.BoundaryOperator` | `class` | `(False,'missing')` | **DECLINED** |
| `orpheus.sn.geometry.SNMesh` | `class` | `(False,'missing')` | **DECLINED** |
| `tests.sn.test_snstreamingoperator.test_apply_…` | `func` | `(False,'missing')` | **DECLINED** |
| `orpheus.sn.spatial.pole_angular_closure` | **mod** | `(False,'missing')` | **DEAD** ✅ |

`judge()`'s last clause re-checks the target's HEAD *carrying the original role*:
`candidate_paths(head, namespaces, role)`. For a single-segment head like `orpheus` under a
non-`mod` role, `bare_module_guess` fires (`"." not in target and role != "mod" and not
hasattr(builtins, root)`), so the head is treated as relative → with no namespaces the candidate
tuple is `()` → `any(())` is False → DECLINED. `:mod:` is exempt from that guard, which is exactly
why only `:mod:` dead targets are ever reported on a page.

⟹ **on an `.rst` page the gate reports `:mod:` and NOTHING else.** One line fixes it — the head of a
dotted path *is* a module reference:

```python
head_role = "mod" if "." in target else role
if not any(lookup(c)[0] for c in candidate_paths(head, namespaces, head_role)):
```

`[M]` blast radius, patched COPY vs shipped, both run on a pristine `git archive HEAD` tree so
REPO_ROOT stays self-consistent: `docs/` goes **49 dead / 71 sites → 207 dead / 255 sites**. The
gate is blind to **158 dead targets across 184 sites in `docs/` alone**, every one a fully-qualified
`:class:`/`:func:`/`:meth:`/`:attr:` on a page. The patched copy flags exactly the 3 ERR-026 roles
pre-edit and zero post-edit — a positive control on the instrument (vv #17) that the count-diff
alone could never give.

⚠ **My first attempt to measure this in-process was itself broken** and read "0 dead" for BOTH arms
while a subprocess on the same tree read 49 — monkeypatching `g.judge` and calling `g.main()` twice
does not work (module-level memo/lru_cache state). Caught only because 0 contradicted a 14 I had
already measured on a subdirectory. Patch a COPY and run it as a subprocess.

### 2. The corpus's own prose is the corroborating oracle — count both spellings

Before de-roling, count how the SAME name is already spelled elsewhere. `[M]` inside the ERR-026
entry: `MorelMontryAngularSweep` **5 literals / 3 roles**, `SNMesh` **4 / 2**, `BoundaryFaceFlux`
**2 / 2**, `transport_operator_matvec` **2 / 1**, `LegacyTauSymmetricInterpolation` **2 / 3**. The
later phases (D, E, …) had already settled on literals; the roles were confined to the earlier
sections. So the entry was *already internally inconsistent* and de-roling made it consistent — that
census IS the evidence the ruling is right, and it is one command. Same trick found a sibling page
(`docs/theory/methods/sn/curvilinear_one_group.rst:2525`) already spelling the deleted test as
``a literal`` and calling it "Phase B's empirical test" — the exact phrasing to copy.

### 3. ⭐⭐ A SURVIVING CLASS does not license keeping the role — the surviving CLAIM does

The brief's LIVE table said keep `MorelMontryAngularSweep` as a role (it exists, at
`orpheus/sn/sweep/pole_angular_closure.py:1308`). I literalised all 3 anyway, because the criterion
is not *does the symbol resolve* but *does a working link mislead about what THIS sentence says* —
the same reasoning the brief itself used to literalise the live `SNMesh`. `[M]` the Phase-B site
describes the class "with starting condition ``ψ_{1/2}=0``", and the SAME entry's Phase-D section
records that as `ZeroSeed`, "Phase B's hardcoded `psi_half_left = 0`", replaced by
`psi_half_seed: … = field(default_factory=CarlsonInwardSweep)`. A link from that sentence lands on a
class whose default contradicts it.

⟹ the rule, stated in the page so it is checkable: **a name is a ``literal`` whenever the sentence
around it describes the code as it then was; a role is used only where the sentence is a
present-tense claim about something that exists now.** That single sentence adjudicates all 29 —
including why the five `:mod:` roles STAY (their sentences are "Documented in X" / "Citations
updated in X", present-tense claims I verified: `[M]` Bailey-Morel-Chang 2010 present 9× in
`reduced_operator.py`, 2× in `transport/spatial/diamond.py`, 11× in `sn/sweep/pole_angular_closure.py`).

**Nothing is lost by literalising, because the live pointers move to ONE place** where their tense
is present — a head-of-block `.. note::` that declares the convention AND says where the objects
went. That is the brief's own "live pointers belong in the status/catcher fields", realised.

### 4. Two brief classifications refuted, both by the same probe error class

* **`orpheus.derivations.continuous.sood_registry` is LIVE**, not "file missing" — it is a
  **package** (`sood_registry/__init__.py` + `la13511.py` + …), imports clean. A `.py`-only
  existence check misses a package. → kept as a role (6 live targets, not 5).
* **`SNMesh.pole_angular_closure` is LIVE**, not a dead attr — set on the INSTANCE at
  `orpheus/sn/mesh/augmented_mesh.py:399` (`self.pole_angular_closure: PoleAngularClosureBase =
  closure_cls(self)`). My own AST index missed it for the same reason a `hasattr(Cls, …)` probe does
  (L-053c). It still became a literal, but for the *sentence-tense* reason, not a dead-target reason —
  and getting the reason right is what stops the next reader "repairing" it back.
* Third, smaller: `orpheus.geometry.boundary` is a live **package**; it is the CLASS `BoundaryOperator`
  that is gone — and a live homonym exists at `orpheus/numerics/operator.py:437`
  (a `_BlockRoleMeta` marker, unrelated). Repointing would have been a false attribution (L-017).

⟹ **before calling a dotted target dead, decide WHICH segment died** — package, module, class, or
attribute. The four have different repairs and only one of them is "de-role".

### 5. The same category error one register DOWN: raw file paths

The ruling fixes roles. It does not touch the *other* present-tense claim a history block makes:
a ``tests/…/foo.py`` **path**. `[M]` in the ERR-026 entry, **14 of 14** distinct `tests/*.py` paths
no longer exist (`tests/sn/spatial/` → `tests/sn/sweep/`, `tests/sn/l1_analytical/` →
`tests/sn/verification/…`). Catalogue-wide: **40 of 100** distinct raw file paths written as
literals are gone — 31 of 72 `tests/`, 9 of 24 `orpheus/`. A raw path warns at no severity, is
invisible to the xref gate (which judges roles), and to nexus (which judges targets).

⟹ the note's third sentence is the prophylactic that matters most: ***which* tests catch ERR-NNN is
never prose — it is the `@pytest.mark.catches("ERR-NNN")` marker set**, `nexus errors` /
`context('vv:error:ERR-NNN')`. Write that once at the head of a history block and the whole class
stops being minted.

### 6. Mechanics

* **Guarded splice, all asserts before the write**: per-replacement counts; an exact
  `len(out) == len(src) + Σ n·(len(new)−len(old))` arithmetic delta; unchanged line count for the
  swap step; the final role list compared to an explicit expected list; `not re.search(r"`{3,}")`;
  `.. error-entry::` count unchanged at 79; and the decisive one — **`src[:i] == out[:k]` and
  `src[j:] == out[m:]` around the entry's own boundaries**, which proves byte-identity of the other
  78 entries in one line.
* **Roles resolve ≠ roles link.** `[M]` none of the 5 classes in my note has an `id=` anchor
  ANYWHERE in the fresh build, so all 5 render plain text — as do their 16–29 sibling sites each
  elsewhere in the corpus. Keep the role anyway: it is the corpus convention, it becomes a link the
  moment the module is surfaced, and — the real argument — **a role is machine-checked by the xref
  gate and a literal is unchecked forever.**
* **A role→literal sweep leaves ONE ragged paragraph per shrunken run.** `[M]` 25 sub-55-char lines
  in the edited region, 24 pre-existing; exactly one paragraph got 3 short lines in a row. Re-wrap
  that one (guard: `new.replace(" ","").replace("\n","") == s.replace(...)`), leave the rest — a
  line-local diff where every changed line shows exactly one swap is what makes the review cheap.
* ⚠ **I broke my own build-sequencing rule twice** (L-054): launched the verification build, then
  found a re-wrap, then found an over-claim, then found a self-inconsistency — four builds. Each
  find was correct and cheap in isolation; the fix is to run the *self-consistency* pass on new
  prose (does my own declared rule hold for every name I wrote?) BEFORE the first build, not after.
* ⚠ **Verify a successor claim against the retiring COMMIT BODY, not the successor's existence.**
  I first wrote that `solution_to_angular_flux*` / `transport_operator_matvec*` "were absorbed into
  the SN operator algebra". `[M]` `4a53737e` says the codec family "became orphan in production"
  after the bare-ndarray contract collapsed at every leaf, and `975edc51` deleted the matvec helpers
  as "without a remaining call site" — they were **retired outright, with no successor**.
  `SNStreamingOperator` really was re-layered (`400ca33d`: `SNSolver.L` → `StreamingOperator` +
  collision multiplier). Same paragraph, two different fates; "absorbed" was true of one and false
  of two.

---

## L-063 — An ONTOLOGY OVERTURN: rewriting the page whose thesis was refuted

**Task (2026-08-19, CS3 step 5, branch `refactor/cone-field-algebra`).** The code carve
had landed (4 commits): flux moved from an *affine space* `𝔸` over a difference space `V`
to the *positive cone* `K ⊂ V` of an ordered vector space. `FluxRole` and the whole
`transport/displacements/` package (8 modules, 7 leaves) were deleted. My job: make the
corpus teach the cone, with the affine era kept as dated history.

### (a) ⭐⭐ A DOC-SIDE OVERTURN IS NOT A RETIREMENT SWEEP — the unit is the ARGUMENT, and
### the load-bearing edit is re-deriving arguments whose CONCLUSION survives.

The retirement grep (32 dead Python-domain refs) was the *easy* half and it finished in one
pass. The hard half had no dead symbol in it at all: `operator_algebra.rst` carried a
**five-obstruction proof** that `Carrier[Representation, Role]` is structurally impossible,
and its obstruction **(a)** read *"the Flux role must make `flux + flux` **raise** while the
Source role must make `source + source` **succeed**"*. That premise is now false — and the
**conclusion is still true**. A retirement sweep either deletes the obstruction (destroying
a correct proof) or leaves it (a false premise under a true theorem). Neither is right.
⟹ **re-derive the argument from what survives, keep the conclusion, tombstone the example.**
`[M]` the live tree hands you the replacement: `AngularFlux.__dict__` has **no `__add__`/
`__sub__`** (MRO `AngularField → BulkField → Field → ABC`) while `AngularSourceSink` **does**
(the iso→per-ordinate containment injection), so the axis that "changes the arithmetic
interface" **inverted** — Source, not Flux. Obstruction (a) survives verbatim in force with
a different worked example plus a second leg the old text never needed (class identity *is*
units identity, and erasure would collapse `type(self) is type(other)` across every role).

⚠ The same shape recurred five times on one page: the role-axis asymmetry section, the
`(Moment, Displacement)` contrast, the fibration note, the conclusion sentence, the
vv-kind table. **Grep the retired symbol to FIND the sites; then read the enclosing ARGUMENT
to decide the edit.** A per-site symbol swap would have shipped five false premises.

### (b) ⭐⭐ RETIRE the eq-LABEL when its NAME encodes the refuted concept; KEEP it when only
### its ADJECTIVE is stale — and the discriminator is the label's BODY, not its name.

Four `:eq:`-cited labels. `[M]` 0 `@pytest.mark.verifies` markers on any (grep
`orpheus/ tests/`), 3 external `:eq:` citers, all in prose I was rewriting anyway.

| label | body states | fate |
|---|---|---|
| `affine-torsor-algebra` | the RETIRED claim (4 torsor axioms) | **RETIRED** → new `flux-vector-algebra`; 2 citers repointed |
| `affine-contraction-ratio` | ρ = ‖Δψⁱ‖/‖Δψⁱ⁻¹‖ — still TRUE, still shipped | **RENAMED** `iterate-contraction-ratio` |
| `affine-true-error` | ‖Δψ‖/(1−ρ) — still TRUE, still shipped | **RENAMED** `iterate-true-error` |
| `affine-typed-residual-eq` | r = (L+C−S−B)ψ − q — untouched by the overturn | **KEPT + annotated** |

⟹ The residual one is the interesting call. Its `affine-` prefix is a **historical artefact
of the page's former title**, not a claim: the residual role was never affine. Its *section*
anchor `affine-typed-residual` has **8 cross-doc `:ref:` citers** (boundary_conditions ×6,
coupled_block_operator ×1, self), and a cross-doc dangling `:ref:` renders plain text with
**no warning at any severity**. Renaming buys cosmetics and risks a silent break.
⟹ **KEEP, and put a `.. note::` at the anchor saying the prefix is stale and why** — so the
next reader greps `affine`, lands there, and is told in one paragraph rather than
re-litigating it. A stale NAME is not a false CLAIM; only a body can be false.

⭐ And the retired equation still had to be *shown* (the history section needs it). Solution:
display it as an **UNLABELLED `.. math::`** with one parenthetical saying why —
*"a labelled equation is an API; these lines state a retired claim, so they must not be
citable."* An unlabelled block cannot become an `:eq:` API by accident.

### (c) ⭐⭐ The sentinel bookkeeping is a SAME-FILE constraint you can check in 1 s — and
### the matrix is GENERATED, so never hand-edit it.

`tests/_harness/audit.py:405` — `.. vv-status: <label> documented` must name a `:label:`
in **the same file**, and `documented` is the only legal status. Renaming a label without
its sentinel is a hard audit violation. Run the scanner directly (sub-second, no pytest):
`from tests._harness.audit import _scan_theory_equations as scan; scan(Path("docs/theory"))`
→ `.violations` / `.documented`. `[M]` 0 violations; population **539 → 540** (4 old → 5 new;
I added `positive-cone-definition` because the cone is the page's new subject).
`docs/theory/verification/matrix.rst` regenerates at `builder-inited` from `conf.py`'s
`_GENERATORS`, so the sentinel list fixes itself — report the post-regen number, never edit.

### (d) ⭐⭐ REPRODUCE the witness, and the reproduction may REFUTE the gate's own prose.

The ruling's decisive measurement is *"DD does not preserve K, so a ψ≥0 type would refuse
production output"*. The gate `tests/sn/solve/test_cone_membership_witness.py` freezes
`min ψ = −6.399383e-01`. I reproduced it through the public entry — exact to the digit —
**and the gate's docstring is wrong about its own fixture**: it says the pair differs in
*"ONE parameter (`nx`) … half the optical cell size"*, but `_solve(nx=2, width=20)` and
`_solve(nx=4, width=40)` both have `Δx = 10`, i.e. `Δx·Σ_t = 100` **identical in both legs**.
⟹ **The argument is STRONGER than its prose** (holding the cell size fixed kills the
"different discretization scale" explanation outright), so the fix is to publish the correct
framing, not to weaken the claim. I also ran two scans nobody asked for and they turned one
frozen number into the mechanism: the **cell-SIZE** scan reproduces the textbook DD
positivity limit exactly (`Δx·Σ_t = 1` in K at `+5.8e-2`; `= 2` already out at `−8.7e-1`),
and the **cell-COUNT** scan at fixed `Δx·Σ_t = 100` shows `nx=2` is the *only* in-cone row
(`nx=3,4,5,6` → `−6.42e-1, −6.40e-1, −6.38e-1, −6.36e-1`). Two tables, ~90 s, and the page
teaches a mechanism instead of quoting a constant.
⚠ I may not edit `tests/` — so the docstring correction is REPORTED, not applied.

### (e) ⭐ Two SKILL files the brief did not scope carried the retired ontology as a POSITIVE
### precedent — and my own repair would have imported the falsehood via its cross-reference.

Brief scope was "the one skill file" (`coding-elegance` #18). `[M]` grep `.claude/skills/`:
`cross-domain-frames/reference.md` (the A.1 frame row + §192/§201 fix-suggestions) and
`numerical-bug-signatures/SKILL.md` §479/§488 also cite `FluxDisplacement` / "flux states
are an affine space" as live. And #18's *corrected* text points readers at A.1 for the frame
— i.e. **the repair cites a stale page** (`coding-standards`' "a cross-reference is a
load-bearing dependency"). ⟹ I flagged the staleness **inline at the pointer** (*"A.1's
frame is sound; its ORPHEUS worked example is NOT"*) so the repair cannot import the
falsehood, and reported both files as owed follow-ups rather than editing out of scope.

⭐ The **reversal** of an anti-pattern is more valuable than its statement, so #18 was
rewritten to lead with what survives (*"NEVER STRAND the convergence data — give it a home
on the object that knows 'previous', which is the ITERATION"*) and to carry the falsified
version verbatim beneath it, plus the checkable test the reversal yields: **(a) is there a
canonical zero?** (distinguished by the domain, not chosen) **(b) is superposition
physical?** Two yeses ⟹ vector space, one type, diagnostics on the record.

### (f) ⭐⭐ The two-sided rule the whole overturn distils to — worth quoting into any page
### that invokes "make illegal states unrepresentable".

Mint the invariant **iff** (1) every value the type admits is legal **AND** (2) every legal
value is admitted. **Half 2 is the one that gets skipped, because it is a claim about the
PRODUCERS, not about the concept.** When it fails the invariant does not prevent a bug — it
**refuses correct output**, and the pressure is then to weaken it, silence it, or route
around the type. Here half 2 failed twice independently: algebraically (K is not closed
under difference or negative scaling, and increments/errors/Krylov directions all live
outside it) and numerically (DD ships negative flux).

### (g) The mechanical residue, all measured

- **Dead Python-domain refs: 32, not the briefed 23.** The brief (and the step-3 commit
  body) counted `orpheus.transport.displacements.*` only; `orpheus.transport.fields._flux_role.*`
  is another **9**. Grep BOTH retired module paths, not the one the commit message names.
- ⚠ **`orpheus/transport/displacements/` still IMPORTS** — an untracked `__pycache__` leaves
  a PEP-420 namespace package (`__file__ is None`, 0 members), so a naive
  `importlib.import_module` probe reports it LIVE (L-052's known false negative). Probe a
  SUBMODULE, or check `__file__ is None`.
- **`tools/check_docstring_xrefs.py`: HEAD 1 dead → working tree 0.** The one it saw was the
  `:mod:` I fixed. It is BLIND to the other 31 (L-062's unlanded `head_role` bug) — my own
  import probe over **727** orpheus-rooted roles across the 8 edited pages is the real gate.
- **`fuel_behaviour.rst:303` "Displacement-Based Constraint" is a MECHANICAL displacement**
  (fuel pellet) — the overloaded-word false positive the brief's grep list contained. Triage
  by MEANING before touching.
- **Build: `-E -W --keep-going` EXIT=0, WARNING/ERROR/CRITICAL set byte-identical to the
  pre-edit `-E` baseline (both empty), 0 `SyntaxWarning`.** 11 anchors + 5 equations + 26
  live code links rendered on the new page.

### (h) ⭐ The changelog contract BLOCKED the obvious home — and the page-local one was right

`docs/theory/methods/sn/history.rst` states *"a new entry lands with its merge hash or not at
all"*, and CS3 is unmerged. `operator_algebra.rst`'s history has the escape hatch
(*"entries marked (in development) live on an unmerged feature branch"*). ⟹ I gave
`field_algebra.rst` its **own** Development history following that convention verbatim, put
a short row on `operator_algebra.rst`'s (its Role axis genuinely moved), and left
`history.rst` alone except to tombstone its 2026-06 row's *affine half* while explicitly
preserving its *typed-residual half* — one row, two halves, opposite fates.

---

## L-064 — Seeding a NEW page from a design-dialogue record: the dialectic is the deliverable, and a retrodiction table is a plan-§2 trap in the corpus

**Task.** CS1 step 5, campaign 1 ("operators born bound"), branch
`feature/cs1-energy-space`, 2026-08-20. Write `docs/theory/foundations/spaces.rst`
from scratch as `field_algebra.rst`'s sibling (that page owns the ELEMENT algebra;
this one owns the SPACES), register it, add one `automodule`, and micro-edit three
`cone_violations` sites. Source of record: `.claude/plans/cs1_energy_space_design.md`
§A/§B/§F + Appendix A (a preserved user⇄agent design dialogue, rounds 2–6).
Result: 1158-line new page, `-E -W` EXIT=0, warning set unchanged (0 both sides),
141 role targets 0 dead, `DEAD TARGETS: 0`, vv violations 0, sentinels 540 → 541.

### 1. ⭐⭐ A DIALECTICAL SEED PAGE is its own doc shape — and it is NOT the 9-step close-out arc

The close-out arc is for a CLOSED investigation whose answer is *"this cannot
work"*. A **seed page** is the opposite event: a design dialogue CONVERGED, the
first slice shipped, and the page exists so the next phase builds on the reasoning
rather than re-deriving it. The user's own instruction was the giveaway — *"make
sure we don't lose it so the archivist can write it later — it was hard to steer
until we got it out."*

The shape that worked, in order:

1. **Key Facts** — including the doctrine's *one-line discriminator tests* verbatim.
   A doctrine that cannot be applied in one sentence has not been articulated.
2. **The taxonomy** (what the new type IS — four slots, a `list-table`).
3. **The theorem** (why one instance's measure is FORCED, not chosen).
4. **The doctrine, DIALECTICALLY** — the question, then version 1 REFUTED with its
   refuting question, then version 2 REFUTED with its refuting question, then the
   standing doctrine, then the retrodictions.
5. **Fences** (what is NOT built, per phase).
6. **Development history**.

⭐ **The refutations must be typographically first-class.** I used
`.. admonition:: ⛔ The refuting question — …` / `:class: error`, one per refuted
version, titled with the QUESTION rather than with the verdict. Rationale: a
reader skimming for the answer stops at the standing doctrine; a reader who meets
only the final statement re-derives version 1 within a week, because **both refuted
versions are almost right**. Version 1 (compactness) reproduced both prior rules
and failed only on energy; version 2 was a one-word patch. The refuting question is
the transferable content — *"where does energy sit?"* and *"what is the measure of
(0,∞)?"* — not the verdict.

⭐ **Name what the doctrine does to the tension it settled.** The two prior rules
(the report's §I.9 retract vs §I.11 quotient) were BOTH right, about different
clauses; what was missing was a second FORK nobody had stated. Writing *"it does
not pick a winner"* explicitly is what stops the next reader hunting for the loser.

### 2. ⭐⭐ A RETRODICTION TABLE is `plan-authoring` §2's aspirational-row trap, moved into the corpus

I wrote *"Every entry below is a layout the tree already ships"* over six rows —
and row 6 was the **buckling member**, which is a prediction of campaign 2. Caught
by my own "count every universal you publish" habit, before the build.

The defect is exactly the plan rule: a table headed by a property of the tree reads
ENTIRELY as a survey of what IS, so one aspirational row is indistinguishable from
the observations. In a PLAN that costs a session; in the CORPUS it is worse — the
table is the doctrine's evidence, and a reader who later discovers one row is
unbuilt has grounds to discount all six.

⟹ **A published confirmation table gets a STATUS column, in the row.** Not prose
above or below it. Mine became `[M] **ships**` × 5 and
`⛔ **NOT built** — a prediction (campaign 2)`, plus a `⚠` lead-in naming the row.
And the honest re-heading is *"rows the doctrine was NOT built from"* — which is
the actual epistemic claim (a retrodiction is a prediction of something not used to
build the theory), not *"layouts the tree ships"*.

### 3. ⭐ Cross-referencing an SSOT: name the REGISTER your page owns, not just the fact

The brief flagged that `energy-condensation-counting-measure` already existed
(`frame.rst`, inside `sn-energy-condensation`) and forbade a twin. L-060's rule is
"cite the SSOT + say which claim THIS page owns". The sharpening this task added:
**the two treatments are the same fact in two REGISTERS, and naming the register is
what makes the second treatment not a twin.**

- `frame.rst` owns it in the **measure register**: `w_g = 1`, not `w_g = Δu_g`,
  derived from Hébert Eq. 3.96/3.97 (distribution vs function averaging), gated by
  rate preservation.
- `spaces.rst` owns it in the **metric register**: `G_E = I`, hence `V ≅ V*`
  isometrically along energy, hence the adjoint there is the plain transpose —
  and the fact that `EnergyAxis` now REFUSES weights at construction.

I derived it a third way (covariant group-INTEGRALS × contravariant group-AVERAGES)
as an *unlabelled* `.. math::` chain, cited the SSOT's label for the claim, and
opened with `.. important:: Single source of truth … Edited there, consumed here.`
**Net new labels on a 1158-line page: ONE** (`spaces-axis-product`, the space =
⊗ axes / shape = concatenation / metric = ⊗ of factor measures identity), sentinelled
`documented` with a rationale naming the foundation battery.

### 4. ⚠ Two agreeing sources can both be wrong about a DATE — git is the arbiter

The brief said the byte gate held *"dated 2026-08-21"*, and
`tests/sn/architecture/test_monomorphic_leaves.py:668` independently says
*"CS1 step 3b (2026-08-21)"*. `[M]` `git log --date=short` puts every CS1 commit
(`1afff47b` … `6da1b23c`) on **2026-08-20**, and the session date was 2026-08-20 —
i.e. both surfaces carry a FUTURE date. Two independent agreeing sources felt like
corroboration; they were one mis-dating copied forward. Publish the git date.

### 5. ⭐ The "only producer" claim, and the collision the grep that proves it exposes

To publish *"every harmonic-moment space is still legacy"* I needed the closure
argument, not the universal: `of_axes` is the only ROOT producer of an `axes`
record (`*` and `dual` merely THREAD one, so both need an axis-built ancestor).
Publishing the derivation instead of the universal is strictly better — it stays
true as the tree grows.

⭐ The grep that established it (`grep -rn "axes=" orpheus/`) surfaced a live
gotcha worth its own note: **`mm.axes` is the GEOMETRIC tuple (`(AxisMesh,)`) and
`mm.bulk_space.axes` is the SPACE-FACTOR tuple (`(EnergyAxis, Axis)`) — the same
attribute NAME on one object, neither derived from the other.** Production already
imports around it (`from orpheus.numerics.axis import Axis as SpaceFactorAxis`).
Three senses of "axis" live in this corpus (space-factor / geometric / symmetry);
the page opens with a note enumerating them and pointing at the rename issue (#393).

### 6. The `V` collision, and the sibling page that had already solved it

A NEW page assembled from multiple sources is the prime site for a within-document
symbol collision (L-011/L-034) — and here it was `V` (the function space, the
page's subject) vs `V` (cell volume, the weight that makes clause 1 work). They
meet in EVERY clause-1 sentence. `field_algebra.rst` had already ruled on it
(`V_cell` written out in full), so the fix was to **adopt the sibling's ruling and
say so**, not to invent one — 5 sites swept programmatically with a `count == 1`
assert per substitution, plus a note citing the sibling for the same reason.

### 7. Numerical evidence a doc page can carry when nothing converges

A terminology/architecture page has no convergence table (the routine "weakest
dimension"). Three measurements carried real weight here and cost ~90 s:

- **The identity bridge, demonstrated**: quotient point (`volumes=[1.0]`) and a
  one-cell slab (`volumes=[2.0]`) mint `energy(2,)*spatial(1,)#<digest A|B>` — same
  shape `(2,1)`, different digests, `==` is `False`. ⭐ Published the name FORM, not
  the hex: a digest in prose is a stale-number risk with no SSOT but the code.
- **The quotient weight is CONSUMED, not decorative**: same flux `φ=(1,1)`, same
  mixture, production rate `0.225` (V=1) vs `0.450` (V=2) through
  `IntegratedReactionRate.evaluate → mesh.volume_measure`. That turns a doctrinal
  claim ("the pairing consumes it") into an arithmetic one.
- **The axis laws**, each run: ones→None canonicalization (`==` and `hash` agree),
  `-0.0`≡`+0.0` bytes, non-finite refused, signed weights legal, `synthetic !=
  from_grid`, `EnergyAxis != Axis`, weighted `EnergyAxis` refused.

### 8. Findings reported, not fixed (scope discipline held)

- `docs/theory/foundations/infinite_medium.rst:1153` and its `:1243-1244` code
  block still teach `basis_shape=(ng, 1)` and a bare `from_solver_data(mat_xs=…)`
  as the current homogeneous spelling. `[M]` I RAN both the doc's four lines and
  production's: **both execute and give `k_inf = 1.8750000000000009`** — the
  keyword survives as an optional override, so this is *stale description*, not a
  broken block. That measurement is what turned a "fix on sight" impulse into a
  correct FLAG: a behavioral rewrite of another chapter's worked example is exactly
  the in-passing rewrite my own rule forbids.
- `orpheus.numerics.space` and `orpheus.transport.mesh.material_mesh` are
  `automodule`'d NOWHERE, and `docs/api/homogeneous.rst` uses `:noindex:` — so
  `FunctionSpace.of_axes`, `has_coordinate_cone`, `MaterialMesh.bulk_space` and
  `solve_homogeneous_infinite` render PLAIN TEXT (measured: 0 `href`s). Page
  convention, not a defect; surfacing them is its own task.

### 9. Gate sequencing (what it cost)

I built **five** times where two would do — every extra build bought by an edit
made after launching a build. The self-consistency pass (universals, symbol
collisions, aspirational rows) must be run to EXHAUSTION before the first
verification build, not interleaved with it. What DID work: a single python
self-check script asserting short-underlines / ladder order / per-table column
consistency / widths-sum / label+anchor uniqueness / role import-resolution /
`:eq:`+`:ref:`+`:doc:` corpus resolution — re-runnable in ~2 s, and it caught every
structural defect before any build.

`docs/theory/verification/matrix.rst` regenerated on the `-E` build and shows the
CS1 campaign's uncommitted TEST work (9868 → 9920 collected; new rows for
`numerics/test_axis`, `numerics/test_space_of_axes`,
`homogeneous/test_operator_spaces`, `homogeneous/test_byte_stability`;
`architecture/test_monomorphic_leaves` 102 → 98, corroborating the four deleted
strict-xfails) alongside my own `540 → 541` sentinel row. Legitimate by-product —
report it, never revert it.

---

## L-065 — Resolving an N-WAY CONTRADICTION: when the corpus states one object three
## incompatible ways, the disagreement IS the diagnosis (a hidden parameter is unnamed)

**Task, 2026-08-23.** Record F-0 of `.claude/plans/frame_square_recarve.md` — a landed
metric repair — across `foundations/spherical_harmonics.rst`, `foundations/frame.rst`,
`conventions/normalization.rst`, `foundations/operator_algebra.rst`,
`verification/error_catalog.rst`. Branch `feature/cs1-energy-space`. Code/tests already
landed and OFF LIMITS (read-only). 6 files, +1246/−81. `-E -W` EXIT=0, warning set
unchanged (0 ↔ 0); vv violations 0, sentinels 541 → 545; `DEAD TARGETS 0`; my own
import probe over the 16 roles on added lines: 0 dead.

### 1. ⭐⭐ The brief handed me THREE published statements of "the adjoint of M". All
### three were internally consistent. That is not three bugs — it is ONE missing parameter.

The corpus said, in three places:

| site | claim |
|---|---|
| `frame.rst` eq `galerkin-strict-adjoint-vs-reconstruction` | "the strict adjoint is the NAKED synthesis (no factor)" |
| `spherical_harmonics.rst` eq `hilbert-adjoint-equals-metric-times-S0` | `Π* = g_C·S₀`, `g_C = 4π/(2ℓ+1)` |
| `normalization.rst` prefactor-ledger row | "The adjoint Π*: carries **Nothing** — the naked reconstruction" |

Worse, `frame.rst`'s **equation** and the **prose four lines below it** disagreed with
each other *inside one admonition* (naked vs `g_C·S₀`) and had shipped that way for
months, warning-free.

⟹ **The reflex "find which one is right and fix the other two" is WRONG here.** An
adjoint is defined by a PAIR of inner products; every one of the three was the correct
adjoint under a *different* coefficient metric (Euclidean / continuum Gram / Parseval
inverse), and none of them named its metric. The repair is therefore **not** N local
corrections — it is *naming the parameter*, once, in a table at the point of definition:

    | Coefficient metric | Where it lives | The Π* it induces |
    | Euclidean          | the bare-transpose reading | S₀ |
    | continuum g_C      | `SphericalHarmonicSpace.from_L` | g_C·S₀   ⛔ pre-F-0 |
    | Parseval G⁻¹       | `FrameBase.basis_space`         | S₀∘G⁻¹ = R/W  shipped |

Each of the three sites then becomes a POINTER into one row, and none can rot
independently again.

⭐ **The generalisable tell, and it costs nothing to look for: two published statements
of the same object that disagree, where each is defended by a correct-looking argument,
means a parameter both arguments quietly fixed differently.** Do not adjudicate; find the
parameter. (Same shape as `vv-principles` #24(b): when a ranking is explained by a
mechanism nobody was debating, the debate was mis-framed.)

### 2. ⭐⭐ A DESIGN PROBE goes stale against the repair it motivated — and it does so
### SILENTLY, because it still runs and still prints plausible numbers.

The plan cited `scratch/probe_f1_parseval.py` for its headline `Parseval ratio 118.7`.
The probe reads `G_stored = frame.test_space.inner_product_weights` — which pre-repair
was the continuum metric and **post-repair IS the Parseval metric**. Run today it prints
`ratio = 1.000` in the row labelled *stored*. It has not broken; it has silently changed
what it measures.

⟹ Three consequences, all mandatory:
- **Never cite a pre-repair probe path as the reproducer of a post-repair page.** Publish
  the **CONSTRUCTION** instead — I wrote the exact 6-line recipe (build the frame, draw
  `default_rng(1234)` unmasked, synthesize, analyse, read five residuals off five named
  attributes) so the table regenerates from the page with no file dependency. (L-048's
  "describe a probe, never cite an ephemeral path", sharpened: the reason is not only that
  `scratch/` is untracked — it is that the probe's SEMANTICS moved.)
- **Re-measure every number against the LIVE tree with your OWN probe.** Mine
  (`probe_f0_doc.py`, ~60 lines) reproduced the theorem/Parseval/closure table for all
  6 sphere families and refuted two inherited figures (below).
- **Report the probe's staleness upward** — the main agent owns `scratch/`.

### 3. ⭐⭐ A SEED-DEPENDENT number must be published as its BOUND, never its value —
### and then find the exact quantity hiding behind it.

Plan: `Parseval ratio 118.7`. Mine, same rule (LS4, L=1): **81.4**. Both correct; the
ratio is a *moment-energy-weighted average* of the per-ℓ factors `(4π/(2ℓ+1))²`, so it
moves with the coefficient draw. Publishing either bare number invites a future session
to "fail to reproduce" a true result.

⟹ Publish the **draw-independent** statement: *it lies between the extreme factors
PRESENT AT THAT L* — `[17.5, 157.9]` at L=1, `[6.3, 157.9]` at L=2 — *and can therefore
never be 1*. ⚠ note the quantifier: I first wrote `[6.3, 157.9]` for a sentence covering
both L values, which is false at L=1 (the ℓ=2 factor does not exist there). A bound is a
universal; `plan-authoring` §2 applies to it.

⭐ **Then look one level down for the number that IS exact.** The *ratio of the two
adjoints* on a single-ℓ unit input is `(4π/(2ℓ+1))²` — measured to `≤2.8e-16` relative at
every ℓ, seed-free, and strictly more useful than the average it produces. A
draw-dependent aggregate almost always has an exact per-mode parent; find it and publish
that.

### 4. ⭐⭐ `.. no-implementation::` has a class the taxonomy was missing:
### AN IDENTITY BETWEEN TWO QUANTITIES THAT ARE EACH COMPUTED.

L-059/L-060 give the classes `{identity, law, canonical-form} → NONE` and
`{typing-rule, definition} → look for a carrier`. This task produced a sharper case, and
it is the one most likely to be mis-declared because it *looks* declarable:

`φ = Mψ = Gc`. Both sides ship — `Mψ` is `_FrameAnalysis.apply`, `Gc` is
`FrameBase.discrete_gram` — and the **identity itself is evaluated nowhere**. Same for
`d_ℓ·G_ℓ = W`: both factors ship, their product is never formed (that is the POINT — the
identity is what lets the kernel carry one `1/W` scalar instead of a per-ℓ table).

⟹ **Declaring either side asserts that one of them IS the identity.** Use
`no-implementation :kind: identity` and say in the block *which* symbol computes *which
side*, plus what the suite measures instead (here: the identity's CONSEQUENCE, the
isometry). `[M]` this stood 17 + 16 name-token guesses down to 0 on two labels; a third
(`galerkin-strict-adjoint-vs-reconstruction`, a *contrast* between two separately-declared
faces) went 2 → 0 and one of its two guesses was `solve_sn_adjoint`, an SN solver entry
point that never touches a spherical-harmonic face.

⭐ The mirror, same session: `sh-space-metric` had NO declaration and 3 genuine
implementers (`metric_per_ell` → `_padded_metric_tensor` → `from_L`). The re-derivation
of a sibling equation's `implements::` set is where you find them: `metric_per_ell` LEFT
`hilbert-adjoint-…` (it is the continuum Gram, no longer that equation's factor) and had
to LAND somewhere — an equation losing an implementer is a prompt to ask which equation
gains it.

### 5. The `documented` sentinel marks the KIND, not the coverage — sibling consistency
### decides, and the matrix legitimately lists a label in TWO places.

`hilbert-adjoint-equals-metric-times-S0` is verifies-covered (`[M]` 2 → **9** tests after
F-0) *and* carries `.. vv-status: … documented`, so the generated matrix lists it under
both "verified by N tests" and "Documented-only equations". That reads like a bug and is
not: on this corpus `documented` marks *representational / face-distinction / literature*
KIND, and its siblings on the same page (`scattering-spectral-theorem`,
`galerkin-strict-adjoint-vs-reconstruction`, `moment-projection-transpose-T`) all do the
same. **Do not "clean it up"** — you would be re-categorising a whole convention from one
label, and it moves a generated artefact. Keep the directive; write the RATIONALE comment
if it is missing (this one had none — I added one naming all four gates).

### 6. ⭐ A three-way SYMBOL COLLISION, and the resolution that keeps every established
### spelling: rename only the one with no constituency.

The frame page already used `W` for **the coefficient space** (`R : W → V`) and
`⟨·,·⟩_W` for **the quadrature-weighted nodal metric**; the code and
`normalization.rst`'s ledger use `W = Σ w_n` for **the scalar total weight**. My
derivation needed a fourth: the weight **matrix**.

⟹ L-051's rule (keep the code's spelling, pay with a `.. note::`) resolves this cleanly
once you notice the four are not symmetric: three have constituencies (a page convention,
a page convention, the code + the ledger), the matrix has none. So write the matrix as
`\mathrm{diag}(w)` — never `W` — and open the section with a `.. warning::` naming all
three survivors and stating the rule. Cost: one admonition; benefit: every equation in a
650-line section is unambiguous.

### 7. ⭐⭐ The reusable close-out shape for a LATENT defect: "why nothing caught it" is
### THREE independent shields, and the third one is the dangerous sentence.

I wrote it three times (frame.rst, spherical_harmonics.rst, error_catalog.rst) and it is
the load-bearing pedagogy of the whole record:

1. **Consistency is not correctness.** The defining adjoint identity
   `⟨Mψ,c⟩_g = ⟨ψ,M*c⟩_W` held at the round-off floor (`[M]` `9.5e-16` at L=1, **exactly
   `0.0`** at L=2) — because `.H` is *built from* the stored metric. It is true for
   **every** SPD metric and therefore carries **zero information about which one is
   installed**. The instrument that CAN fail compares the metric to something defined
   without it (here Parseval: the field's own norm).
2. **Composed chains are immune** — interior metrics cancel, so the production kernel
   never reads a face's `.H`, and the 0-ULP canary is green by construction.
3. **Only end-of-chain adjoints are exposed, and there were none.** `[M]`
   `grep -rn "analysis\.H\|reconstruction\.H" orpheus/` → exactly one hit, a docstring.

⛔ **Shield 3 is where a close-out goes wrong.** "No consumer exists" is not safety, it is
**latency** — write it that way, with the clock: *the defect becomes live with the first
adjoint consumer, which is why the metric had to be right before those land*. A page that
reports shield 3 as reassurance teaches the next session to defer.

### 8. Extending an ERR entry vs minting a new one — the decision, and the one thing to
### check first.

F-0 is the THIRD chapter of ERR-039 (Wave 0 → Phase 1 → F-0), all "metric / transpose /
adjoint conflation on the same operator pair", each one level deeper: *wrong operator* →
*right operator, unasked metric* → *right Gram, WRONG SIDE*. Extending was right, and the
decisive check is not narrative tidiness — it is that **the landed gates already carry
`catches("ERR-039")`**, so a new number would silently orphan them and I cannot edit
`tests/`. ⟹ **Read the catching tests' markers BEFORE choosing the ERR number**; the
marker set is the constraint, the narrative is not.

Also: mark the superseded chapter IN PLACE (`⛔ superseded 2026-08-23; see the F-0
chapter below`) on the *bullet* that states the retired formula, not only in the new
chapter — a reader who lands on the Phase-1 list must not read it as current.

### 9. Corrections that the WIDENED sweep produced (the brief named 4 of 6 sites)

The brief's grep list found `spherical_harmonics.rst`, `frame.rst` ~2700, and
`normalization.rst:161`. Running the sweep myself added:
- **`frame.rst` "Numerical evidence"** — item 2 still read `M* = g_C S_0`, ~1100 lines
  from the note the brief pointed at (L-044's "audit the PARAGRAPH FAMILY, not the diff").
- **`frame.rst` Schur bullet** — `g_C` as "the SO(3) Plancherel weight": TRUE (it is the
  continuum Gram) but now one reciprocal away from the frame's metric; qualified rather
  than changed.
- **`operator_algebra.rst:3295`** — "the addition-theorem `R`, not the W-weighted adjoint
  of `M`": the *negation* survived the repair while its vocabulary died. Post-F-0 `R = W·M*`
  exactly, so the sentence is improved by stating the relation instead of denying one.
- **`normalization.rst`'s ledger row was the THIRD contradiction** and the brief did not
  know it existed — it was found by grepping `W-weighted`, not `g_C`.
⭐ And the payoff nobody asked for: that page's own "unification the canon misses"
(`(2ℓ+1)/W`, W = 4π sphere / 2 slab) **IS** the Parseval metric `G⁻¹ = (2ℓ+1)/W`. A sweep
for staleness turned into the strongest single piece of corroboration on the page —
*always read what the stale site was TRYING to say.*

### 10. Findings reported, NOT fixed (code/tests are off limits)

- **`orpheus/numerics/frame.py:116-119`** (the `_DISCRETE_GRAM_DIAGONALITY_RTOL` docstring):
  says the slab live off-diagonals "sit at ~0.5 of the Cauchy–Schwarz scale
  `√(G_jj G_kk)`". `[M]` relative to the C–S scale they are **0.9347**; **0.5774** is
  relative to the largest DIAGONAL. Verdict unaffected (threshold `1e-10`), but the
  *stated normalisation* is wrong and the same wording is copied into two
  `tests/numerics/test_frame.py` docstrings. ⟹ when a docstring quotes a ratio, check
  WHICH denominator — two plausible ones differ by 1.6× here.
- **`spherical_harmonic_space.py`** class docstring's `inner_product_weights` parameter
  still says "row ℓ holds 4π/(2ℓ+1)". `from_L`'s docstring WAS updated by the F-0 commit;
  the class-level parameter description was not — and the frame-dressed object IS a
  `SphericalHarmonicSpace` (built by `dataclasses.replace`), so it is present-tense-false
  for the majority instance. A half-done docstring sweep, exactly `vv-principles` #21.
- `scratch/probe_f1_parseval*.py` no longer reproduce their own headline (see §2).
- The `check_docstring_xrefs.py` `.rst` blind spot (L-062's one-line `head_role` fix) is
  **still unlanded** — re-confirmed: it gates only `:mod:` on `.rst` pages, so my own
  import probe was the acceptance evidence for the 16 roles I added.

### Quality self-assessment

| dimension | score | note |
|---|---|---|
| Derivation depth | 5 | φ = Gc from three re-associated products; the metric SOLVED for, not asserted; both general adjoint sandwiches collapsed step by step; the SH-specific collapse isolated as its own labelled identity |
| Cross-references | 5 | 16 roles added, 0 dead by import probe; 5 new anchors, every cross-doc `:ref:` verified as a rendered `href` |
| Numerical evidence | 5 | 7-frame × 5-residual table + the slab refusal table + the indicator instance + the pre/post ratio table — every figure re-measured this session, with the construction published |
| Failed approaches | 5 | the three-metric table IS the failed-approach record; ⛔ pre-F-0 equation preserved unlabelled; the 3-shield "why nothing caught it"; the refuted LS₈ 24 % claim |
| Code traceability | 5 | 17 declared `implements::` across 4 labels + 3 `no-implementation` blocks, all pre-flighted against `graph.db` and verified in the rebuilt graph |
| Derivation source | 4 | derived from the LIVE code + my own probe (no `derivations/` script exists for frame algebra; the plan's probes are pre-repair scratch — flagged) |

## L-066 — A FACTORY-TIER retirement is a THREE-TENSE sweep, and the disposition is decided by the SENTENCE, not by the symbol

**Task (2026-08-24, branch `feature/cs1-energy-space`, campaign 1 CS4b S5).** CS4b S5
retired the mesh-keyed sugar tier on every transport field leaf — `from_mesh(values, mesh)`,
`zeros_on(mesh)`, `from_ndarray(arr, mesh)` deleted from `AngularField` / `ScalarField` /
`FaceField` and every concrete leaf, `MomentField.from_ndarray` too, plus the
`spatial_moments=` int on those factories. Replacement is SPACE-primary: `Leaf(values=…,
space=…)` and `Leaf.zeros(space)` on the carrier's cached mints (`angular_bulk_space`,
`bulk_space`, `angular_trace`, `scalar_trace`, the two ψ½ spaces), with a NEW mint
`SNMesh.angular_trial_space` replacing the retired int. Composites went space-keyed
(`FullField.zeros(..., space=)`, `RadialCharacteristicField.flux_zeros/source_zeros`,
presence-gated). Brief supplied a ~30-hit grep list + three disposition rules.
Gates: `-E -W` EXIT=0, W/E/C/SyntaxWarning set unchanged (**0 ↔ 0**); `DEAD TARGETS 0`;
my own import probe over the 37 qualified roles on added lines = **0 dead**; vv violations
0, sentinels 545 unchanged; 6 doc files, +360/−106.

### The load-bearing finding: three tenses, one symbol, three different repairs

A factory retirement scatters the SAME token across three grammatical registers, and the
symbol grep cannot tell them apart. Sorting the 30 hits by TENSE gave a clean 3-way split
that turned out to be the whole job:

| register | example (verbatim) | repair |
|---|---|---|
| **live guidance** — prose/tables/code telling the reader how to build a field TODAY | *"``from_mesh_laws`` returns exactly ``zeros_on``"*; the ladder table's bottom rung; a runnable `code-block` | re-word to the space-primary spelling **with the right carrier mint**; a wrong mint is a fresh falsehood, so measure it |
| **history** — a landed change's own narrative, an ERR-NNN post-mortem, a "before X, callers hand-rolled…" | *"a hand-rolled ``from_mesh(trace.values.copy())`` beside a Pattern-4 factory"* | **prose STAYS** (past tense is correct history); a `:meth:` role pointing at the deleted target is DOWNGRADED to a ``literal`` keeping the exact old name |
| **landed-but-written-as-future** — a step record whose "NEXT sub-step" already shipped | *"When they land, the only change at the factory call sites is passing ``spatial_moments=…``"* | re-tense to past + a dated `.. note::` saying where the capability lives now |

⭐ **The third register is the one a symbol grep is worst at and the one that costs the
most**, because it reads as a *plan*, not as a claim, so nobody re-checks it. `[M]` in
`cartesian_multid.rst` an entire "Construct-general, select-narrow — what this step does and
does NOT do" subsection was present tense about a state two campaigns old: *"No production
field selects it. The ``spatial_moments`` factory parameter defaults to ``1`` at EVERY call
site and is NOT auto-read from ``mesh.scheme.spatial_basis_per_axis``"* — while
`solver.py:920` and `streaming.py:999` both pass exactly that expression, and the leaf route
had since moved to `angular_trial_space` entirely. Repair shape that worked: **re-tense the
bullets in place** (adding *(as of S3)* to the one that is a dated observation), flip the
section TITLE's verb too ("does and does NOT do" → "did and did NOT do" — free, since the
section carried no `.. _anchor:`), and append ONE dated `.. note::` carrying the live
mechanism + the measured shapes + the named survivors. Do NOT delete the bullets: they are
the reason the capability was built default-OFF, which is the durable content.

### Two-step ownership history is worth its own paragraph (the brief's disposition (c))

The `indexing_and_layout.rst` allocator section had ALREADY been corrected once — a 2026-08-10
`.. note:: **Correction (Issue #346)**` recording that the allocator moved OFF `SNMesh` and
ONTO the leaf. S5 moved it again (off the mesh KEY onto the space KEY). ⟹ the honest shape is
a **second `.. note::` beside the first**, not a rewrite of the first: *"#346 moved the owner;
S5 moved the key; the leaf is still the owner"*. Both notes carry a one-command `[M]`
(`[n for n in dir(SNMesh) if "zero" in n.lower()] == []` / `hasattr(AngularField,"zeros_on")
is False`). A reader landing on either correction now sees the whole arc.

⭐ **And the trap inside that section: it QUOTED a production docstring that no longer
exists.** The old text read *"The production docstrings say so in as many words — «the uniform
leaf-side allocator … replaces the retired `SNMesh.zeros_*` mesh-side factories»"*. `[M]`
`grep -rn "uniform leaf-side" orpheus/` → **0 hits**; `Field.zeros`'s docstring now says the
opposite-keyed thing. A quotation is a claim about a FILE, invisible to every gate (no role,
no label, no path). ⟹ **grep the quoted STRING, not just the symbol, whenever a doc quotes
code.** Same class as the raw-file-path defect of L-062, one register up.

### Measure the replacement before you write it (the mint is a choice, and one is wrong)

The re-word target is not mechanical: `zeros_on(mesh)` maps to **five** different mints
depending on the leaf family, and `angular_bulk_space` vs `angular_trial_space` is a real
choice at every angular site. I built one probe and published its table:

`[M]` vacuum slab, `N = 4` (GL), `ng = 2`, `nx = 4`:
`angular_bulk_space (4,2,4)` · `bulk_space (2,4)` · `angular_trace (16,)` ·
`full_field_space` blocks `(4,2,4)+(16,)`; **DD: `angular_trial_space is angular_bulk_space`
→ True** (same cached instance). Same slab under `LinearDiscontinuous`: trial reads
`(4,2,4,2)`, axes `angular(4)·energy(2)·spatial(4)·spatial_moment(2)`, the moment factor
MODAL with measure `(1, 1/3)`. Sphere (4 cells, GL(4)): ψ½ blocks `(16,)` cells + `(4,)`
corner; the same `flux_zeros` call on the slab RAISES the R12a diagnosis.

⭐ The `is`-identity at DD is the single most useful sentence in the rewrite — it is what
makes "read the trial mint" safe advice at *every* scheme width, so the doc's worked example
can carry one line instead of a branch.

### The doc code block was ALREADY broken, in a way unrelated to the retirement

`verification/sn.rst`'s composite-source example built `TimedFullField(bulk=…, boundary=…)`.
`[M]` the dataclass fields are `interior` / `boundary` — `bulk=` has not existed for
campaigns. Found only by RUNNING the block (lesson-2 rule), which I had to do anyway to pick
the mint. Fixed in the same edit and pinned the result in prose (`[M] max φ = 1.8265` on the
stated fixture) so the next reader can falsify it. ⟹ **a retirement sweep that touches a code
block owes that block a RUN**, and the run finds defects the sweep was not looking for.

### Housekeeping observed, not caused

* `docs/theory/verification/matrix.rst` is regenerated by the `builder-inited` hook, so the
  FIRST `-E` build of a session materialises whatever drift the landed code/test commits left
  (here `10063 → 10067` collected, `test_meshless_construction 5 → 8`). It shows up in
  `git diff` and is **not** a hand edit — report it, never revert it, never edit it.
* My added-role probe's one "DEAD" hit was an **unqualified** role
  (`:meth:`BulkField._compose_spatial_moments``) on a line I only re-tensed. Unqualified roles
  resolve by Sphinx module context and the gate skips them by design — a false positive of the
  probe, not a defect. Check the diff (`git diff | grep '^[+-].*symbol'`) before believing the
  probe on an unqualified target.

---

## L-067 — Documenting a MEASURED machinery: the brief's numbers were a sample, the record's `[M]` carried a confound, and the gate that certifies me is blind to its own class

**Task.** CS4b S7 docs half (branch `feature/cs1-energy-space`, 2026-08-24): write the
flagship "axis collapse pair" section on `spaces.rst`, repair `infinite_medium.rst` for the
EE-1 typed-rate landing + the K2 pose split, fix the `field_algebra` mesh-identity rows, the
frame-square `Π* = R` contradiction, and the changelogs. Six pages, +1407/−104. Final
`-E -W` EXIT=0 with the WARNING/ERROR/CRITICAL/SyntaxWarning **set unchanged (0 ↔ 0)** from a
freshly measured baseline; vv violations 0, sentinels 545 → 549; `dead_references` 1 (a
pre-existing PRODUCTION docstring, not mine).

### 1. ⭐⭐ "Bit-exact" is a property of the DRAW, and a gate's seed is a sample of size ONE

The brief and the gate docstrings both said `R ∘ E = id` is `[M]` **BIT-EXACT**. It is not a
property of the operators. `R∘E` computes `Σ_n w_n·(φ/Σw)`; whether the round-off floor is
*zero* depends on how those products re-associate for the particular numbers involved.

`[M]` on the gate's own synthetic fixture (`w = [0.3,0.7,0.5,0.5]`, `Σw = 2.0`):

| row | `np.array_equal` holds on |
|---|---|
| `R∘E = id` (G6.1) | **1156 of 2000** seeds — fails on 844, worst rel `1.480e-16` (~1 ULP) |
| `P = E∘R` idempotent (G6.2) | **143 of 200** seeds — fails on 57 |
| both, on the shipped SN carrier (GL4 slab) | **200 of 200** |

Both gate rows are green because the seeds they hard-code (`1` and `2`) land in the exact
set. **Change the seed and they red.** The SN carrier is robustly exact because `Σw = 2`
exactly *and* the symmetric GL weights re-associate cleanly — which is what licenses
`np.array_equal` on the production-facing rows and does NOT license it on the synthetic one.
⟹ Publish a **bound over ≥200 draws** with the norm and the seed family written out
(`max‖a−b‖_∞/‖b‖_∞` over `default_rng(1000+k)`), never a single reading; and say WHICH
fixture is robustly exact and why. Reported the seed-fragility upward — I do not edit `tests/`.

⚠ The dual, same session: the **tightness** rows (minted kernels vs the literal frame's face
contents) ARE robustly bit-exact — `[M]` 200/200 on all three correspondences — because both
sides evaluate the same reduction in the same order. So "bit-exact" is sometimes a property
of the construction and sometimes of the draw, and only the measurement separates them.
Saying which is the whole content of the row.

### 2. ⭐⭐ A COINCIDENCE claim needs its family, and this one is false exactly where it matters

Brief + production docstring + gate docstring all carry: *"the gram einsum is bit-identical to
`weights.sum()` on 8 of 8 probed fixtures (n ∈ {2,4,5,6,16,64} incl. GL64's inexact Σw)"*.
`[M]` mine, two weight families:

| n | `leggauss` weights | `linspace(0.1,1.3,n)` |
|---|---|---|
| 2, 4, 5, 6 | identical | identical |
| **16, 64** | **NOT identical** (`2.0000000000000004` vs `2.0`) | **NOT identical** |

The Gram is `einsum("n,nj,nk->jk", w, T, T)`; `ndarray.sum` is a pairwise reduction. And on
the **shipped SN quadratures** the split lands at **GL8**: divisor `1.9999999999999998` vs
`quad.weights.sum() = 2.0`, so `AngularSourceSink.from_isotropic` differs from a hand-written
`Q/Σw` by `2.0e-16` relative in production. ⟹ The structural claim (*the divisor IS the
frame's `discrete_gram[0,0]`*) is exact by construction and is what a gate must pin; the
coincidence with `weights.sum()` is a fixture accident and must never be relied on. Published
the ladder as a table with its two columns side by side.

### 3. ⭐⭐ A design record's `[M]` can carry a CONFOUND — two settings moved together, one got the credit

The record read: *"sphere GL L=1 (DIAGONAL Gram): `face.H(e₀φ) == E(φ)` to 2.2e-16; Slab L=2
(DENSE Gram): un-physical"* — attributing the split to **geometry**. It cannot be geometry:
the angular frame is built from `sn.quad`, which knows nothing about the spatial coordinate
system. Running the crossed cell:

| fixture | Gram max off-diag | `face.H(e₀φ)` vs `E φ` | `reconstruct(e₀φ)/W` vs `E φ` |
|---|---|---|---|
| slab L=1 | 5.6e-17 | 5.6e-17 | **0.0** `array_equal` |
| sphere L=1 | 5.6e-17 | 1.1e-16 | **0.0** |
| slab L=2 | 1.155 | 16.17 | **0.0** |
| sphere L=2 | 1.155 | **16.17** | **0.0** |

⟹ the discriminator is the **Gram's diagonality**, i.e. **L**, and geometry is inert. A 1-D
polar rule has no azimuthal nodes, so the m≠0 modes are not orthogonal under it — `[M]`
`gauss_legendre(8)` at L=2 reads the SAME `1.155` / `16.17`, so refining the order does not
fix it. The clean, **metric-free** statement (`E = reconstruct(e₀·)/W`) is bit-exact in all
four. ⟹ **When a record's two arms differ in more than one setting, run the crossed cell
before publishing either as the cause.** The correction made the section stronger: it is the
reason the collapse pair is minted from an *indicator* frame instead of lifted out of the
harmonic one — it must keep working where the harmonic metric does not exist.

⚠ And: the committed probe `scratch/probe_s6_q5_dissolution.py` is the SLAB L=2 arm — the arm
the record itself calls un-physical — so run as committed it prints `1.617e+01` while a
**production docstring** cites that path for `2.2e-16`. A scratch probe cited from `orpheus/`
is a raw-path claim about a file that no instrument checks (L-062, one register up).

### 4. ⭐⭐ The xref gate's `head_role` bug is ROLE-scoped, not `.rst`-scoped — my own memory was too narrow

L-053/L-062 recorded this as *"on an `.rst` PAGE that gate reports `:mod:` and nothing else"*.
`[M]` it is worse and simpler: `judge(target, ns, role)` re-checks the target's HEAD **carrying
the original role**, and `candidate_paths("orpheus", ns, "meth")` returns
`('<namespace>.orpheus',)` — which never resolves — so every **dead** fully-qualified
`orpheus.*` target under a non-`mod` role is DECLINED, in `.py` docstrings exactly as in
`.rst`. Live ones return ALIVE earlier and are unaffected, which is why the gate looks healthy.

- `judge("orpheus.numerics.space.FunctionSpace.definitely_not_here", role="meth")` → **DECLINED**
- `judge("orpheus.numerics.does_not_exist", role="mod")` → **DEAD**

⟹ **`DEAD TARGETS: 0` certifies `:mod:` targets and nothing else.** The one-line fix
(`head_role = "mod" if "." in target else role`) applied to a COPY, run as a subprocess from
inside the repo (it resolves paths against `REPO_ROOT`, so a `/tmp` copy scans 0 files), read
**1 dead target / 2 sites** where the stock gate read 0 — one of them my own new xref, one a
pre-existing production docstring. ⟹ the acceptance evidence for a page is still YOUR OWN
import probe; the gate is a `:mod:` check.

### 5. ⭐ Two independently-vocabularied instruments agreeing IS the acceptance evidence

nexus `dead_references` (resolves by RENDERED target) and the patched gate (resolves by
IMPORT) returned **exactly the same single finding** — `FaceField.from_face_arrays` at
`face_layout.py:355`. Neither alone would have been persuasive: the stock gate said 0, and
nexus's set-difference with the gate is normally noisy (L-052). Convergence from two different
resolution mechanisms is what makes a one-line report actionable.

### 6. ⭐ A brief can name the wrong CLASS for a method — and the same error is already in production

The brief said *"`FaceField.from_face_arrays` is the typed entry"*. `[M]` `hasattr(FaceField,
"from_face_arrays")` is **False**; it lives on `BoundaryField`. I wrote the brief's spelling
into a changelog row and my own selfcheck caught it — and the SAME wrong class is in
`orpheus/numerics/face_layout.py:363`, which is where the brief's author almost certainly read
it. ⟹ a brief's symbol claim inherits the tree's own errors; `hasattr` every method-on-class
before minting a role.

### 7. Repair shapes worth reusing

- **A self-contradicting Key Facts block, 12 lines apart.** `frame.rst` promised
  `GalerkinFrame ⟹ Π* = R` and, twelve lines below in the same admonition, `M* = R/W`. The
  post-F-0 truth is neither: Galerkin fixes *which basis* the adjoint re-synthesises on (the
  trial one, `M* = S₀∘G⁻¹` — a **canonical** dual), and the metric stays. Fix at BOTH poles
  (`Π* ∝ R` in the diagram + a ⚠ clause naming ERR-039/051 and the indicator counter-example),
  never at one.
- **A "single-sourced through X" claim is two claims.** `operator_algebra.rst` said the
  `iso + aniso` dunder is single-sourced through `from_isotropic`. `[M]` the dunder's body is
  `self.values[None] + other.values` — the **plain** broadcast — while `from_isotropic` applies
  `1/Σw`. They differ by exactly the axis's total weight, i.e. they are the two arrows the whole
  section I was writing exists to keep apart. The repair writes the ⚠ *and* points at the new
  section, so the falsehood becomes the worked example.
- **A retired guard tier leaves a stale REASON attached to a surviving FACT.** `verification/sn.rst`
  said the composite is re-homed "because `TimedFullField` algebra enforces mesh identity". The
  re-home still happens; the reason is now space CONTENT. Keep the instruction, replace the
  reason, and say what changed — a reader who trusts the old reason will "optimise away" the
  re-home for a twin carrier.
- **The production helper's own docstring lied the same way.** `_require_typed_composite`'s
  docstring says *"(2) `field.interior.mesh` is the operator's SAME `sn_mesh` instance"* over a
  body that compares `field.interior.space != field.interior.space_on(sn_mesh)`. Reported.
- **A colliding bare step number.** `spaces.rst` said the V/V* condensation morphisms are
  "scheduled for S7" — a plan-internal number that collides with CS4b's own step S7, which
  landed that day and built none of it. Disambiguated at BOTH sites (`plan-authoring` §9b in the
  corpus).
- **A fence row that fell.** "Only the scalar bulk is axis-built; every other space is legacy"
  → `[M]` the angular bulk and the scheme-widened trial space are axis-built too and report
  `has_coordinate_cone is True`; what is still legacy is the composite and the flat traces.
  Re-title the fence to what is actually still fenced.

### 8. The changelog routing, again

`methods/sn/history.rst` contracts *"a new entry lands with its merge hash or not at all"*, so
an unmerged branch is BLOCKED there — while `spaces.rst`, `field_algebra.rst` and
`operator_algebra.rst` each carry the *(in development)* escape hatch. ⟹ route the entry to
the page whose SUBJECT moved and which permits the hatch; report the SN row ready-to-paste.

---

## L-068 — Discharging a merge-hash contract: the blast radius is the BRANCH NAME, not the blocked page

**Task (2026-08-24).** `feature/cs1-energy-space` ff-merged to `main` at `55bb47b9`
(90 commits, 264 files). My held contract — `methods/sn/history.rst` contracts
*"a new entry lands with its merge hash or not at all"* — was finally
dischargeable, and the dispatch was scoped to "write the row(s)".
Landed `68d265ef` on `docs/sn-history-campaign1-landing`.

### 1. ⭐⭐ The merge event falsifies the SIBLING pages the same instant

L-067 taught the *routing* rule: when `history.rst` blocks you, route the entry to
the page that carries the `*(in development)*` hatch (`spaces.rst`,
`field_algebra.rst`, `operator_algebra.rst`). What it did not say is that the
hatch is a **debt**, and the merge is what calls it in. The moment the branch
merges, every `*(in development)* branch ``<name>``` cell is present-tense-FALSE —
and nothing points at them, because the dispatch names only the blocked page.

⟹ **On discharging the contract, grep the BRANCH NAME across `docs/` first.**
`[M]` `grep -rn "cs1-energy-space" docs --include="*.rst"` → exactly **3 cells**
(2 on `spaces.rst`, 1 on `field_algebra.rst`), all replaced with
`merged @ ``55bb47b9`` —`. One command, one minute, and it is the difference
between a corpus that agrees with git and one that says three campaigns are still
in flight. ⚠ Also grep `"in development)\*"` corpus-wide — that catches a hatch
whose branch was named differently, and it found the standing *explanatory
sentence* on `operator_algebra.rst` (a convention note with no instances, which
correctly STAYS).

### 2. ⭐ A DATE in a prose history block is a git question, and it drifts by one

`frame.rst:4744` read `**2026-08-24 — step F-1, the mint**`. `[M]`
`git log --format="%h %ad" --date=iso -1 3dfea889` → **2026-08-23 16:19:54**.
Its F-0 sibling four paragraphs up was right; the S6.0b block below was right.
One block, one day off — written from "which session was I in", not from git.
This is lessons §1's *"Git is the arbiter for dates"* (L-064) applied to a prose
changelog rather than a claim, and the tell is free: **when you cite a commit's
date into a NEW row, you have already looked it up — diff it against every
existing prose block naming the same commit.**

### 3. Row granularity: group by THESIS, and let the page's own precedent settle it

The merge carried five architecturally distinct milestones. One consolidated row
or five? The page answers itself: `[M]` the #280 coupled-block campaign holds
**six** rows across 2026-07-05…07-12, every one stamped `(merged @ ``3f0b8c74``)`.
So per-milestone rows sharing one merge hash IS the convention, and the `Where`
format is `` `<step hash>` (merged @ ``<merge hash>``) ``.

⟹ **Group by the THESIS that moved, never by the plan's phase labels.** The
campaign's own step boundaries (S4 / S4-amendment / F-0 / F-1) cut across
subjects: S4 is field-layer, the S4-amendment is operator-layer, and they landed
in one session. The rows I wrote are *"a field is an element of a space"*,
*"the frame owns its metric and mints bound faces"*, *"an operator is not an
operator without its two spaces"*, *"S/F/C are kernels"*, *"the space layer gains
axes"* — five theses, not five phases.

⚠ And **strip the plan-internal tokens on the way in.** My first draft carried
*"the standing R2 hazard"* and *"a ``§6b`` call-site set"* — a plan's risk label
and a rules-file section number, both meaningless in the corpus and both colliding
with live campaigns' own numbering (L-067's bare-step-number rule). Rewrote to
*"the hazard the monomorphic-leaves suite had catalogued"* and *"a call-site set
that is complete by symbol grep"*. The corpus says what the thing IS.

### 4. ⭐⭐ A row whose MECHANISM a later row in the same merge overturned gets an in-place ⛔, not a rewrite

The S4-amendment's step A1 (`6e04a749`, 2026-08-22) bound `HarmonicFrame` itself
to its two field spaces at construction. F-1 (`3dfea889`, 2026-08-23) **reverted
it** — a frame is a shared FACTORY, so the binding belongs on the faces it mints
(`[M]` live: `HarmonicFrame.__init__(self, basis, measure)`). Both landed in the
same merge, one day apart.

Writing A1 in the present tense would have shipped a falsehood; deleting it would
have destroyed the reason the correction happened. So the 2026-08-22 row carries
its (d) clause plus, in place:

> ⛔ **Superseded the next day by F-1** (the row above): the binding belongs on the
> *faces* a frame mints, not on the shared factory. The amendment's *demand* stands
> unchanged — it is what made the misplacement visible in the first place.

That last sentence is the load-bearing one: it says which HALF survived. This is
`plan-authoring` §3 (edit the refuted premise in place) landing in the corpus, and
in a reverse-chronological table *"the row above"* is a correct pointer.

### 5. The gate stack for a changelog-only edit — and why the standard gate is not enough

Baseline and post, both forced `-E`: EXIT=0, W/E/C/SyntaxWarning **0 ↔ 0**;
`check_docstring_xrefs.py` `DEAD TARGETS 0`; nexus `dead_references`
`total_dead 0`; vv-status `violations 0 / sentinels 549` (unchanged — the rows add
no `:label:`); `verification/matrix.rst` regenerated byte-identical.

But per L-062/L-067 that gate is ROLE-scoped-blind, so the acceptance evidence was
**my own import probe over the added lines**: parse every
`:role:`target`` out of `git diff --unified=0 -- docs/` (flatten whitespace FIRST —
roles wrap), strip the `<display <target>>` form and the `~`, then walk
`importlib` + `hasattr`. `[M]` 36 distinct qualified roles, **0 dead** — and it
caught two before the build: `orpheus.data.mixture.Mixture` (lives at
`orpheus.data.macro_xs.mixture`) and `orpheus.numerics.operator.AdjointOperator`
(private `_AdjointOperator`; rewrote the prose to "the generic adjoint wrapper"
rather than xref a private class).

⭐ **And the render check the build cannot do**: slice the built HTML between the
new rows' first and last distinctive phrases, strip tags, unescape, and count
**visible backticks** and **surviving `:role:` spellings**. `[M]` 0 and 0 in my
rows. This is the only instrument that proves the markup MEANT what it said —
see §6 for what it found in the rest of the page.

### 6. ⭐⭐ RST cannot nest inline markup — and the census is exhaustive, so measure it with the ISSUE'S OWN instrument

The same render check reported **84 visible backticks on the page** — none mine,
all pre-existing, all one mechanism: `**bold naming ``a symbol``**`. RST forbids
nested inline markup, so the inner delimiters render literally. Silent at every
severity.

`#379` already owned this, scoped to the error catalogue at `[M]` 32 runs. Rather
than file a sibling, I re-ran **#379's own grep** corpus-wide:
`<(strong|em)>[^<]*``[^<]*</(strong|em)>` → **125 runs across 25 live pages**,
with the catalogue's 32 reproducing exactly (a control that the instruments agree)
at 26 % of the total. Posted as a widening comment with a suggested retitle.

Three things worth carrying:

- ⚠ **Exclude `_build` pages whose `.rst` no longer exists.** `[M]` 12 orphaned
  pages carry a further **76 runs** no source edit can reach — counting them
  inflates the figure by 60 %. (Also exclude `_modules/`: a literal backtick in a
  viewcode source listing is correct output.) Same trap my memory already warns
  about for stale-ref greps, here in a *measurement*.
- ⭐ **The strictly worse sibling ranks above it**: the same RST rule (inline
  markup may not open after `. * ~ § ↔ =`) kills **104 interpreted-text ROLES
  across 28 live pages** — they survive into rendered prose as their own source
  spelling, with the LaTeX backslash eaten. `[M]` the commonest survivor is
  `` :math:`mu` `` (10 sites), then `` :math:`tau` `` (6). A visible backtick is
  ugly; a dead `:math:` is a **missing equation**.
- ⭐ **Both greps are a CENSUS, not a sample**, and say so when publishing: RST
  admits no role spelling in rendered prose and no stray delimiter outside a
  literal block, so a hit is a defect by construction. That is the L-061 argument
  (a warning count is a non-representative sample of a fidelity-loss class) reused
  as the *justification* for the number rather than as a caution about it.

### 7. Numbers re-derived rather than relayed

- The byte gate: the plan says "D5 8/8". Ran it — `tests/homogeneous/test_byte_stability.py`,
  **8 passed**, and its own fixture docstring says *"exhaustive over what the tree
  ships"*, which is the phrase the row now uses.
- The GL8 correction: the plan's banner said the probe's ladder "skipped GL8". My
  draft embellished it to *"a ladder of eight fixtures that broke every arithmetic
  pattern"* — a property of the ORIGINAL probe I could not verify (only the gate's
  CURRENT list is knowable). Cut back to exactly what the gate docstring asserts.
  Same reflex as lessons §1's *"count your own universals"*, applied to an
  adjective.
- CS1's slot claim: my draft read *"the slot the kernel-binding phase then tightens
  to MANDATORY"*. `[M]` CS4a K2b made **F**'s space mandatory, not S's
  (`ScatteringOperator.__init__(..., space: FunctionSpace | None = None)` still
  ships). Scoped to `MANDATORY on :math:`F``.

---

## L-069 — The RENDER is the only gate for inline markup, and a LITERAL is not a role

**Task (2026-08-26).** Archive a three-probe investigation (literature sweep +
SymPy derivation + an original asymptotic derivation) into
`docs/theory/methods/sn/curvilinear_one_group.rst` and
`docs/theory/foundations/discretization.rst`: the tensor-product factorization
of the curvilinear angular-redistribution operator, the τ-arity theorem, the
Padé positivity ladder, the seed cone risk, and a refutation of Morel–Montry's
own 1984 summary rule. Mid-task the code carve I had been told not to name
LANDED, adding three stale-reference repairs.

### 1. ⭐⭐ A double-backtick LITERAL renders a backslash VERBATIM

In a `list-table` of measured values I wrote ``` ``1.4\times10^{-6}`` ``` for two
cells. A literal does exactly what a literal promises: the built page carried
the characters `1.4\times10^{-6}` in prose. `-W` EXIT=0, `-n` blind,
`check_docstring_xrefs.py` blind (it gates TARGETS, not whether a span parsed),
nexus `dead_references` blind. **The only instrument that saw it was slicing
the built HTML and counting raw TeX outside the MathJax spans.**

⟹ **a number in scientific notation is `:math:`1.4\times10^{-6}``, never a
literal.** The discriminator: does the cell contain a backslash? If yes it is
math, not code.

### 2. ⭐⭐ `**``value``**` is the commonest way to mint the nested-markup defect

RST cannot nest inline markup (L-068), and the shape I reached for *fourteen
times in one session* is a bold-wrapped literal in a numeric table cell —
`- **``-0.200000``**` — because I wanted the negative rows to stand out. Every
one rendered as ``` ``-0.200000`` ``` with four visible backticks. Plus one
`**`[M]` this is what ships**` (bold wrapping a `<cite>` span), two more.

⟹ **in a table cell a literal already carries its own visual weight; NEVER wrap
it in `**`.** Emphasis goes in the surrounding prose, which is where the reader
needs the interpretation anyway. Guard: `assert "**``" not in text` before the
write — one line, catches the whole family.

### 3. ⭐ A bare `:ref:` to a section whose TITLE contains `:math:` leaks raw TeX

`The dome closes — :math:`\alpha_{M+1/2} = 0` as an admission contract` is a
section title. A bare `:ref:` to it pulls the title as link text and the math
arrives as the literal characters `\alpha_{M+1/2} = 0`. Pre-existing on **4**
sites of that page — so it is page behaviour, not something I introduced — but
my two new sites made it five, and the fix is free: explicit link text
`` :ref:`the dome-closure contract <sn-alpha-dome-closes>` ``.

⟹ before adding a bare `:ref:`, look at the TARGET'S TITLE. Math in a title ⟹
explicit link text. (Same reflex as the admonition-anchor rule, different
cause: that one WARNS, this one is silent.)

### 4. ⭐⭐ Build the render checker carefully — its own regex is a false-negative source

Two instrument bugs, both of which made the checker useless in opposite
directions:

- **Sphinx emits display math as `<div class="math notranslate nohighlight"
  id="equation-X">`**, so a `<div class="math[^"]*">` strip misses EVERY
  numbered equation and the checker reports ~1000 false raw-TeX hits (it is
  reading correct MathJax source). Fix: `<div class="math[^"]*"[^>]*>`. Same
  for the inline `<span>`.
- **The `<head>` MathJax macro configuration is raw TeX in the page**
  (`"Sigt": ["\\Sigma_{\\mathrm{t},#1}", ...]`), so a whole-page scan always
  reports hits. Slice by `<section id="...">` to the id of the NEXT section.

⟹ **and the source-side regex alternative does not work at all.** A
`\*\*[^*]*``[^*]*\*\*` scan over my new blocks returned **26 suspects, 0 real**:
`**A** … **B**` matches as one run whenever no `*` sits between them, so every
adjacent pair of bold runs is a false positive. The render check found 3 real
classes with 0 false positives. **The rendered page is the instrument; the
source is not.**

### 5. ⭐⭐ Expand the series yourself — "monotone and positive" can still be INCONSISTENT

The source memo offered a positivity/accuracy trade inside the lumped-LD family
and named `(λ,ν) = (0,½)` as *"genuinely monotone at the cost of dropping to
**first** order"*, with transmission `2/((1+τ)(2+τ))`. I re-derived the family
from scratch (nodal DG cell, one free parameter per row) and the transmission
reproduces exactly — and the order label is wrong. `a'(0) = −3/2`, not `−1`:

| cells over a fixed `Σ_t X/|μ| = 1` | 10 | 100 | 1000 | 10000 |
|---|---|---|---|---|
| `(0,½)` | 0.2367 | 0.2245 | 0.2233 | **0.2231** |

It converges cleanly — to `e^{−3/2} = 0.223130`, not `e^{−1} = 0.367879`. It is
`vv-principles` #5 in its purest form: **a correct rate to the wrong limit**,
and both of the properties the memo checked (sign-preservation, `A⁻¹ ≥ 0`) are
perfectly true of it. Consistency is a THIRD property neither test sees.

⭐ The correction was cheap and produced a better object: solving
`a'(0) = −1` symbolically gives `ν = 1 − λ` (a ONE-parameter family, not two),
monotonicity gives `λ ≤ 0`, and the nearest monotone consistent member is
`(0,1)` with `a = 1/(1+τ_opt/2)²` — strictly positive, `A⁻¹ ≥ 0`, genuinely
first order. **A refuted memo claim replaced by a derived one is the best
possible outcome of "verify every number you cite".**

⟹ when a memo states an ORDER, expand the series. One `sp.series(a - exp(-t))`.

### 6. ⭐ Read the class docstring of the object you are theorising about

The carve landed mid-task and `AngularRedistribution`'s own docstring **already
states the tensor-product factorization** and cites the same memo. Two
consequences: (a) my chapter is the theory home for a structure the code
asserts, not a twin — say so; (b) **align to the code's exact spelling**
(`R_spatial ⊗ A_angular(τ, α, w)`), because internal consistency between code
and corpus outranks brevity (L-051). I had drafted `R ⊗ A_ang`.

### 7. ⭐ A `.. vv-status:` sentinel works INDENTED — check before relocating one

`tests/_harness/audit.py`'s `sentinel_re.match(stripped)` matches the STRIPPED
line, and the only rule is same-FILE. I nearly moved one out of a `.. warning::`
block on the belief it needed column 0. Read the scanner (30 s) instead of
reasoning about it. `[M]` 15/15 new labels registered, 0 violations,
documented 549 → 564.

### 8. ⭐ A retirement's stale REASON outlives its stale NAME

Site 1 of the carve's blast radius read *"``alpha_half`` … stay on the geometry
side — they are genuinely geometric"*. The NAMES were the greppable half; the
load-bearing half was the **reason**, which the factorization refutes (the dome
is the ANGULAR factor, a function of `(quadrature, coord)` alone). The same
false reason appeared a second time 1200 lines away in a development-history
item, where the names were correctly past-tensed and the reason was not.
⟹ after fixing a retired name, read the sentence that JUSTIFIES it.

### 9. ⭐ Reproduce a claim from the SHIPPED function, with the shapes it wants

`affine_scan_coefficients` takes `V` at `(N, nx)`, not `(nx,)` — my first two
attempts died on `V[:, None, :]`. Fed correctly, DD and LD reproduce the Padé
ladder to `1.1e-16` / `1.2e-16` over six optical depths, which converts "the
closed forms are the shipped scheme's" from an assertion into a measured bound.
Same for the seed: `carlson_inward_sweep_from_source` on 8 cells at
`Σ_t Δr = 3` returns `+0.4, −0.08, +0.016, …` — ratio `−0.2 = (2−3)/(2+3)`
exactly, i.e. the shipped seed march sign-alternates, measured on production.

### 10. Verified first-hand against the rendered scans (not the memos)

- **Adams–Martin 1992 App. A, printed p. 160** — read the page: (A.1a)/(A.1b)
  carry `+r_kΔr_k`, `−Δr_k²/6` / `+Δr_k²/6`, `−r_kΔr_k/3`. Two minus signs on
  the `ψ^x`-coupled entries; magnitudes match the Gram exactly. The sibling
  removal block `σ_tk[V_kψ + W_kψ^x]` / `σ_tk[W_kψ + X_kψ^x]` is symmetric on
  the same page — the typo argument is visible without leaving the page.
- **Hill 1975 ONETRAN, printed pp. 9–11** — Eq. (30) plain angular diamond
  *pointwise in r*; Eq. (32) applies it to the two-point spatial AVERAGE; (35a)
  shows the redistribution as `(α/w)[ΔA_i; z_5]⊗[1,1]`, manifestly rank-1;
  (36)–(38) the starting direction. The rank contradiction is real.
- **MWS 1996, printed p. 452** — Eqs. (74)/(75)/(76) and the verdict quote,
  verbatim; and **they name the Padé degrees themselves**, so the whole ladder
  framing is literature-backed rather than ours.
- **Palmer–Adams 1993 = UCRL-JC-111847** (the code docstring says
  UCRL-ID-114256, which is Palmer's *thesis* — reported, not edited).
  Their LD verdict is quoted as PREVIOUS work (their ref [5], Palmer–Adams
  1991), a nuance worth preserving.

### 11. My own derivations that replaced relayed numbers

- flat-flux row-1 identity: sphere `A_+ + A_- − 2V/h = 4πh²/3 = R_10` **exactly**
  (symbolic); cylinder both `= 0`, so the gate reads `0 = 0` there — the
  "run it on the SPHERE" rule is a theorem, not a measurement.
- `R_01/R_00 = h/(3(r_-+r_+)) ≤ 1/3` with equality iff `r_- = 0`, so `R` is SPD
  on every admissible cell and `det R ≠ 0` — which is what makes "β = 0 is
  necessary as well as sufficient" true.
- the 2×1 rectangular column `[ΔA ; ΔA·h/(6r_c)] = [ΔA ; 4πh²/3]`, matching
  ONETRAN's own `[ΔA_i ; z_5]`.
- `β⁻/β⁺` half-range split at the M-M τ: `+0.101808 / −0.101808` (N=4) …
  `+0.124610 / −0.124610` (N=32), sum ~1e-17 — reproduces the memo to every
  digit and is the evidence for "β = 0 is a GLOBAL identity across μ = 0".
- `β_e` sign flip `+9.107e-01 (N=2) → −1.111e-01 (N=4)`, and the `|μ_s+1|`
  equivalents `0.1132 / 0.0161 / 0.0053 / 0.0015 / 0.0004`; and
  `morel_montry_beta = 1.5 × β` bit-for-bit, so the shipped instrument IS the
  object the seed analysis needs.

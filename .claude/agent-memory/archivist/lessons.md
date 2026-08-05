# Archivist — Lessons (hot digest)

Read FIRST, every dispatch. One rule per lesson, imperative, with the failure→correction core that
earned it. War stories, evidence and `file:line` detail live in **`lessons_archive.md`** — open a
`→ L-0NN` section on demand. Mechanical HOW is NOT repeated here: build-gating, venv/worktree facts
and the 9-step close-out arc in `AGENT.md`; V&V vocabulary in `vv-principles`; Branch-1/Branch-2 in
`algebra-of-record`.

**THE SPINE — a page is DONE when** every cross-ref resolves against the LIVE tree · every claim
was verified against live code THIS session · every claim's V&V level matches the skill verbatim ·
every retired symbol leaves no present-tense-false mention · the build's WARNING/ERROR/CRITICAL
**set** is unchanged from a freshly-measured `-E` baseline. Every rule below is one face of that.

---

## 1. Ground truth is the LIVE tree — every other surface lies eventually

**Meta-rule: the brief is the FLOOR; live code is the rule.** Brief, docstring, verdict memo,
retirement shim, scanner finding, plan line and "MEASURED" block are point-in-time snapshots.
Verify, then write, then FLAG every scope-expansion the verification forced.

- **Read the live `def`/body before citing any convention, shape, signature or design decision.**
  Seen: a docstring lying about an index convention and a return layout; a verdict memo recording
  the RECOMMENDATION while the code shipped the alternative; a brief naming args the live Protocol
  never takes. → L-001
- **A retirement SHIM's docstring is frozen at the commit that minted it — verify against the
  CANONICAL form it re-exports.** A shim called a cross-class dunder "retired"; the canonical
  modules had since RE-PERMITTED it, so the brief would have past-tensed an accurate section.
  → L-001
- **A brief's discriminator is a heuristic, not a per-ref rule — one phrase can be TRUE at one site
  and FALSE at another.** Resolve EACH ref by the exact symbol's live signature / the site's live
  fixture, never by the phrase. → L-001
- **A brief's RATIONALE can be wrong while its conclusion is right, and a MECHANICAL vocabulary
  swap can restate a FALSE claim in fresh, authoritative words.** Re-verify the CLAIM before
  re-spelling it — grep the live class for the property the sentence asserts; `-W` catches neither.
  → L-001
- **Reproduce every number you cite, and sanity-check its neighbours while the harness is open** —
  one worked example's intermediates contradicted its own result three lines later. → L-001
- **A COMPOSITE's measured identity cannot certify its FACTORS — verify each factor against its own
  live method.** A page's `Rᵀ = (cos w/norm)·Σφ` was `[M]` bit-exact while BOTH per-factor formulas
  it was built from were wrong (the shipped split puts `1/norm` in `B`, not `C`, because `C` makes a
  *current* and `B` an *intensity*); every measurement on the page reads `0.0` either way, so no
  gate could see it. Re-verify the design-probe's description too. → L-049
- **A "MEASURED, do not re-derive" block is a CLAIM** — that means "don't burn a session", not
  "don't check". A bit-identity attribution was wrong on exactly the configuration that motivated
  the change; widen the repro to the WHOLE inventory, since the brief's sample is never the
  population. Two mechanically different effects can also measure the same (reduction-order drift
  vs a real value bug both read ≤1 ULP when the offending weight is `O(ε)`), so a ULP table cannot
  justify such a change — give the structural reason and `.. warning::` that the ULP row is NOT
  evidence of equivalence. → L-043
- **You are the judgment layer over any bulk scanner.** Import-verify the SUGGESTED target (one
  named a class existing nowhere); reject findings whose evidence the current clean build
  contradicts; when a retarget crosses a numerical/structural claim read the successor's live body;
  attribute every residual to a known false-positive class. → L-021
- **A naming-dense brief on a fast-moving branch goes stale FIRST** — import-verify every class,
  helper and line-ref before minting a cross-ref. → L-018, L-039
- **A plan's phase-line and internal task numbers are stale tracking artifacts.** Never infer
  "phase N is a gap" from an open plan bullet — read the shipped page; never trust a plan's bare
  `#N` (internal numbering COLLIDES with real issues); verify every issue you cite. → L-038, L-041
- **Verify a "the gate still does X" claim against the TEST BODY, and COUNT the rows.** A class
  docstring described semantics its body had abandoned, and only 3 of 7 cases were re-posed — my
  draft's "every case" was a fresh falsehood. Likewise `python -c` every numeric constant a doc
  asserts (one was four orders off, unwarned, for months). → L-042
- **An EQUATION has TYPES and a SCOPE, and NO gate checks either — read the domains/codomains, and
  ask which instance the proof covers.** A published `R∘G = R` could not type-check
  (`Γ₊→Γ₋` vs `Γ₋→Γ₋`); nobody introduced it — a CODE carve narrowing the spaces retroactively
  falsified a math sentence three chapters away. Separately, "the crossing is geometric … which is
  why G carries it" was proven for the mirror and stated for EVERY law. A narrowing carve is a
  licence to re-type-check every identity naming the affected spaces; a "which is why X" closing a
  one-instance argument is a licence to re-quantify it. Fix by SCOPING the proof and ADDING the
  missing case — never rewrite a proof that was only over-quantified — and write out the boundary
  cases (they turned out to be a shipped realizer REFUSAL, the best evidence the new law is right).
  Root cause of both: a factored form presented as BOTH a classification and a computational
  recipe; say which it is, and check the declaration tier against the REALIZATION first. → L-048
- **Describe a probe, never cite an ephemeral path.** A `$CLAUDE_JOB_DIR/tmp/` script no reader can
  open is a stale raw path the moment it is written (and `scratch/` is untracked). Publish the
  construction — shapes, metrics, comparison — so the table regenerates from the page. Reproduce AND
  WIDEN every measured number a plan hands you (a 3-sample `|Γ₊|=|Γ₋|` claim became 6 quadratures ×
  every face). → L-048

---

## 2. The build is BLIND to most doc-correctness defects — grep is the gate

**Meta-rule: `-W` proves only "I added no NEW warning". The acceptance evidence for a correctness
sweep is a grep inventory with a per-hit KEEP/FIX adjudication.**

- **Unresolvable `:func:`/`:class:`/`:meth:`/`:attr:` render as PLAIN TEXT with no warning.** After
  any carve that deletes or renames a symbol, `grep -rn "<symbol>" docs/` and repoint every hit.
  → L-002
- **`-n` (nitpicky) is NOT the missing gate.** MEASURED: `-n` saw ZERO of 22 dead refs, because
  Sphinx only nitpicks what it RENDERS and the carrying modules were not `automodule`'d
  (`tests/**` is never read). Edits to such docstrings cannot move the warning count, so "count
  unchanged" proves nothing. → L-044
- **`tools/check_docstring_xrefs.py` IS the gate — run it, don't grep blind.** It resolves every
  FULLY-QUALIFIED role by IMPORTING it, so render coverage is irrelevant;
  `… <tree> --quiet` → `DEAD TARGETS : 0` is the acceptance criterion. Never touch its empty
  ALLOWLIST. It skips UNQUALIFIED refs by design (Sphinx resolves those by module context), and it
  is blind to LITERALS — so after fixing a renamed symbol's roles, grep the OLD NAME tree-wide and
  adjudicate every ``literal`` by tense (`_select_si_resolvent`: 1 dead role + 3 live-prose
  literals on two other pages). All three trees are now at 0: `tests/` 41 dead/62 sites,
  `orpheus/` 30/37, `docs/` 20/24 in 15 pages. → L-045, L-046, L-047
- **Beyond AGENT.md's warn-list, two more DO warn:** a `:widths:`/column mismatch, and `ref.ref`
  "*A title or caption not found*" — a bare `:ref:` to an anchor sitting before a PARAGRAPH (fix:
  anchor a titled/captioned element, or use `` :ref:`text <label>` ``). Raw path strings in prose
  warn NOWHERE. → L-002, L-027, L-040
- **Grep `SyntaxWarning` in the build log too — a case-sensitive `WARNING:|ERROR:|CRITICAL:` MISSES
  it** and it does not bump the exit code. A non-raw docstring with `\Gamma` emits
  `SyntaxWarning: "\G" is an invalid escape sequence` mid-build. Before reporting one in a file
  another agent is editing, `py_compile` the LIVE file — mine was fixed a minute later. → L-048
- **`⭐`/`⛔` have ZERO occurrences in `docs/`** (they are plan / agent-memory / docstring
  vocabulary); `⚠` and `✓` ARE corpus vocabulary. Grep a glyph in `docs/` before importing a marker
  from a plan. → L-048
- **A not-yet-built symbol is a code LITERAL, never a `:class:`/`:meth:`.** Gate with `hasattr`;
  the same probe flips a LANDED seam to a live role. → L-002, L-014, L-025
- **Plain-text refs are often the page CONVENTION, not a defect** — un-`automodule`'d packages, and
  `:noindex:`-automodule'd ones, are plain-text page-wide. MEASURED: `api/method_of_characteristics.html`
  and `api/discrete_ordinates.html` carry **zero** `id="orpheus.*"` anchors, so `:noindex:` renders
  docstrings but mints NO targets and leaves live `href`s pointing at anchors that never existed.
  Adding an `automodule` there is still worth it for the DOCSTRINGS — just don't expect roles to
  link. Match the page, repoint dead refs to the LIVE path, never half-surface 1–2 leaves.
  → L-002, L-034, L-047
- **`automodule`-readiness is MULTI-gate; "0 `:label:`" is necessary, not sufficient** — it also
  trips on an unregistered role, a short docstring underline, a malformed field-list, a member-name
  collision cascading onto pages you never touched, and a closing role-backtick followed by a word
  char. `-E -W`-build EACH in isolation; if a cluster is unready, automodule only the clean module,
  prose-ref the rest, REPORT the unblocking fix. → L-002
- **Labels are PATH-IMMUNE; `:doc:` is PATH-SENSITIVE.** A moved label needs zero referrer edits —
  the break is CONSUMING PROSE naming the old page (`` see :ref:`X` in :doc:`.../oldpage` ``): the
  link goes to the new page while the prose sends the reader to the old one. Sweep the tree with a
  whitespace-FLATTENED scan (the `:doc:` routinely wraps). → L-024, L-026
- **Any doc-cleanup pass is a free staleness audit** — reading a line to trim it is the only gate
  catching a stale RAW PATH or a stale V&V word. → L-034
- **That gate also OVER-reports — its `getattr` probe is blind to an annotation-only class
  attribute (`ClassVar`, dataclass field), which autodoc DOES publish.** 5 of 30 `orpheus/` hits
  were live; `Field.UNITS` renders a real `href` in a FRESH build. Prove a contested hit with a
  rendered-anchor grep, then LEAVE it and report — never mutilate a true ref to green a gate,
  never edit a gate you weren't asked to edit. Mirror class, genuinely unresolvable and worth
  fixing: with `napoleon_use_ivar = True` an `Attributes`/`Parameters` entry mints NO target, so
  an `__init__`-assigned attribute needs a live `:class:` + a literal — **5 of 24 `docs/` sites**,
  and autodoc coverage will NEVER revive them. Phrase the replacement so the sentence says where
  the value comes from ("the ``scheme`` attribute that :class:`SNMesh` realizes in its
  constructor"). → L-046, L-047
- **In `docs/api/`, dead refs cluster by SECTION: the unit of repair is the retired API SURFACE.**
  7 of 24 sites were ONE section listing 6 factories retired in one commit; the successors were a
  re-LAYERING, not renames, so 6 repoints would have been 6 lies. Read the surviving module's own
  docstring FIRST — a well-retired module states its successor map and tells you whether you owe N
  edits or one rewrite. Expect ~⅓ REWRITE on a deletion-driven sweep. → L-047
- **RUN every doc code block a present-tense sentence promises works.** One opened on an import of
  a module deleted months earlier AND used `np` with no `import numpy`; no build sees either. A
  dead import is the loudest possible dead ref. → L-047
- **A `scipy`/3rd-party role can die by UPSTREAM removal** — `scipy.special.sph_harm` was removed
  in 1.17; the successor `sph_harm_y` has a SWAPPED `(n, m)` order, which belongs in the fixed
  sentence, not just the target. → L-047

---

## 3. A `:label:` is a V&V edge — grep the matrix before touching it

- **NEVER rename or delete a label a `@pytest.mark.verifies(...)` targets.** For a stale equation
  that IS a verifies-target, keep the label and rewrite only the BODY. Run the silent-class grep of
  `orpheus/`+`tests/` FIRST: empty ⟹ safe to rename; a hit ⟹ report the test edge (you don't edit
  `tests/`). → L-003, L-032
- **When an ALGORITHM is replaced, a retired-STEP label is usually KEPT-AND-REPOINTED to a
  conceptual survivor** (reflexively retiring iteration-step labels orphans test edges). Ask whether
  the CONCEPT survives, `.. note::` what it historically named, and retire only a documented-only
  label with no survivor. → L-003
- **PHANTOM verifies (marker whose label exists nowhere): repoint if the equation already carries a
  label, MINT only if the law is prose-stated but unlabeled** (one `.. math::` = one label). → L-003
- **Classify every label you add.** Structural / representational / literature-transcribed ⟹ the
  machine-read DIRECTIVE `.. vv-status: <label> documented` + a rationale comment naming the gate
  (prose status does NOT count — a `--strict` audit regresses). A NEW test's verifies-target ⟹
  leave it un-sentineled; never sentinel to paper over a transient orphan. → L-004, L-036
- **Algebra-of-record SymPy-identity labels are verifies-COVERED, not documented** — foundation +
  verifies COEXIST and produce a real edge. Reserve `documented` for motivating/definitional
  literature with no tight gate; the two together is muddy. → L-039, L-035
- **Orphan adjudication.** WIRE iff an existing test's PRIMARY assertion IS this equation against a
  structurally-independent reference ("would a sign flip red it?"). SENTINEL for exactly three
  shapes: a general/continuous identity whose concrete instance is tested under a DIFFERENT label;
  a native-vs-legacy bit-identity regression; code that does not exist yet. GAP only for a
  load-bearing computed contract with no test anywhere — never manufacture one to look thorough (a
  38-label slice legitimately came out 0 GAP). Sibling-consistency dominates, and a ROOT narrative
  page's orphans are ALL sentinel (its formulas are tested under the METHOD pages' own labels —
  name that downstream gate in the rationale). → L-035, L-004
- **Un-sentineling is verified against the LIVE test, not the brief** (a brief says "wired" when the
  marker is still WAITING). After removing the directive, rewrite — don't delete — its rationale
  comment to a plain note naming the gates, so a future auditor doesn't re-add a sentinel. → L-037
- **When the exit condition FIRES but you don't own the generated artefact un-sentineling would move:
  keep the DIRECTIVE, rewrite the RATIONALE.** Open it `⚠ PRECONDITION EXPIRED … REMOVABLE`, quote
  the superseded text verbatim as history, name the exact gate that now exists. Avoids both silently
  re-categorising a generated table and leaving quoted-false text. A sentinel carrying its own exit
  condition still needs somebody to NOTICE it fired — no build does. → L-049, L-048
- **Dropping a duplicated eq-label ALSO drops its `.. vv-status:`, silently DEMOTING the concept to
  orphan** — and `-W` is blind (the orphan gate is a generated REPORT, not a build check). Move the
  status to the survivor. → L-027
- **Backfilling labels on a derivation-mirror page: BARE dominates** (the checkpoints are already
  labeled; the residue is true intermediates). Fill only the recognizable gap classes: a governing
  eq parallel to a sibling page · an unlabeled object the corpus uses BY NAME · a geometry/sibling
  parallel gap · a paper-numbered eq in the page's established family. 2-of-31 is correct. → L-030
- **Section-label and equation-label are DIFFERENT namespaces** that coexist under one name with no
  warning — a name owning both is TWO independent single-home checks. Verify with
  `grep -c '^\.\. _X:'` or `grep -c ':label: X'`, never a raw mention count. → L-024, L-003

---

## 4. Retirement & staleness: three greps, and the unit is the THESIS

**Meta-rule: grep the SYMBOL, the full MODULE PATH, and the CONCEPT'S human paraphrase. Then ask of
each hit's ENCLOSING SECTION: "is the PREMISE still true?"**

- **A DELETION (unlike a MOVE) leaves a stale PARAGRAPH, not a stale token** — a move leaves a
  true-but-relocated symbol, a deletion leaves a sentence whose premise died. Three shapes recur:
  a file that past-tenses a retirement in one section and PRESENT-tenses it in another; a LANDED
  migration still written as future work ("consumers will migrate in Wave G"); a docstring
  contradicting its own body (a documented `None` fallback the code `raise`s on). → L-046
- **Preserve the WHY; tombstone, don't delete.** Flip tenses, keep the logic. When a finding
  invalidates a published table, add `.. note:: **Retraction (date, Issue #N).**` above it — values
  stay, the INTERPRETATION gets the tombstone. → L-007
- **Retitle to the CONCEPT and KEEP the anchor when the concept survives; RENAME the anchor when its
  name encodes a REFUTED concept** (updating every inbound ref in the same pass, verified in built
  HTML). Keeping the anchor is what makes a retired-note section free — cross-doc `:ref:`s keep
  resolving and auto-pick up the new title. → L-007, L-015, L-040
- **A retirement that REMINTS the freed name onto a different live object makes the name a homonym
  across one commit — disposition each mention by what the PASSAGE DESCRIBES, never by the name.**
  Blanket find-replace is wrong; grep the full module path (the bare name cannot tell live from
  dead) and audit the whole split role FAMILY, not the head symbol. → L-017
- **Per-site ladder for a widely-referenced retired entry point:** (a) behavioral rewrite where the
  section teaches CURRENT API — including DELETING a stale code block rather than symbol-swapping
  it; (b) past-tense double-backtick LITERAL in history/changelog narrative; (c) delete where the
  clause carries no content. Build the LIVE-grounded successor table FIRST — the successor is
  context-dependent and a 1:1 rename is forbidden. → L-019
- **When the deletion is a COROLLARY of a design unification, the SECTION'S THESIS is stale, not
  just the line.** The tell is a stale design stated in the PRESENT tense as live rationale ("by
  design", "what stayed deliberately legacy"). Banner the rationale section preserving its
  reasoning, fully rewrite only the one genuinely-stale-as-current contract, and retitle a moot
  future-work section "(obsoleted)" while KEEPING its `:ref:` anchor. → L-020, L-013
- **Grep the CONCEPT, not only the symbol — a `list-table` COLUMN is a documentation surface with no
  symbol in it** (a brief's 7 literal hits missed 17 cells under a paraphrased header). Dropping a
  column is a 3-part edit (header · every value cell · `:widths:`) verified in RENDERED HTML; prefer
  REPLACING it with the true intrinsic property. And the paragraph that JUSTIFIED the retired flag
  inherits the flag's wrongness — re-verify it. → L-040
- **In `tests/`, a dead xref is a TRIPWIRE for a false CLAIM, not a typo** — a test docstring says
  what the test PINS, so the retirement that killed the ref usually invalidated the sentence too.
  Read the test BODY, then REPORT the false claim (never quietly repoint, never fix the gate). Seen:
  a module docstring advertising a unified-vs-legacy bit-identity chain when BOTH implementations
  were deleted; a pin list whose item asserted the INVERSE of the live gate. Adjudicate four ways —
  REPOINT (majority) · PAST-TENSE LITERAL (a role PROMISES a live link) · REWRITE · DELETE (rare;
  0 of 62 here). A not-yet-built module is a LITERAL — and check its cited PLAN FILE exists too.
  → L-045
- **The brief's successor map is a HYPOTHESIS — run `git log --diff-filter=D` on the old path.**
  "X now lives at Y" hid a DELETED legacy class whose replacement merely reuses the name; same name,
  different object splits the sites into history-literals and a rewrite. → L-045
- **A dead ref can sit on a claim a rename INVERTED, not merely moved** — a published code block
  showed a `cast`-based helper whose live body has no cast and owns the guard the prose credited to
  the CALLER; a repoint alone leaves two falsehoods with a working link. Read the live body, then
  re-state the mechanism. → L-047
- **When a page has an OPEN owner issue: fix the dead refs + the MEASURED adjacent falsehoods,
  leave the issue's genuine rewrite item, and comment with a measurement table AND the residue's
  CORRECTED path** (#286 named a page that no longer exists). Neither "defer, it's theirs" nor
  "annex it". → L-047
- **A retirement can DEMOTE a gate's claim class without touching the test body.** When a rewire
  points a comparison at the successor, re-ask "are the two sides still INDEPENDENTLY produced?" —
  if the survivor CALLS the other, the gate became a pass-through check and every doc crediting it
  must be re-scoped (name the real pin). The tell in a diff is a variable still called `legacy`
  beside a brand-new API. → L-044
- **Doc-only "retire the false promise": keep the DECLARATION, make the CLAIM true** — state the
  measured present, how production reaches that information today, and the phase that fills it.
  → L-041
- **Replace an unfalsifiable inventory sentence with a MEASURED table** — "each subclass overrides
  these where applicable" is prose over a lattice computable in ten lines. → L-041, L-042
- **Give each retracted consequence its own `**Disposition:**` — a retraction can INVERT a claim,
  not just kill it** (one conclusion flipped to the opposite type-tag once the domain narrowed).
  → L-042
- **A phase-N doc pass leaves phase-(N−1)'s falsifications behind — audit the PARAGRAPH FAMILY, not
  the commit's diff.** A correctly re-typed section three screens from a sentence the PREVIOUS phase
  falsified ships a self-contradicting page. → L-042
- **FLAG, don't silently rewrite, adjacent SUBSTANTIVE staleness.** Repoint-in-passing is correct;
  behavioral-rewrite-in-passing risks minting a NEW false live claim (worse than a dead ref to true
  history). Exception: fix a refuted-claim survivor on a line you are ALREADY editing. → L-007,
  L-014, L-018
- **When a doc describes a DEFERRED seam that since LANDED, verify the SHAPE that shipped — don't
  just flip "deferred"→"done".** The change routinely closes the seam by a DIFFERENT mechanism;
  separate the surviving conclusion from the stale PREMISE, and re-derive why it still holds.
  → L-007
- **A retirement leaves THREE tense classes, not two — the third is a falsified PREDICTION.** Beyond
  present-false (repoint) and past-history (de-role to a literal), future-tense prose written while
  the replacement was a plan names **a MECHANISM, a HOST/TYPE and a PHASE**, and a landing can
  falsify any subset independently — check each against shipped code. Seen: "the type that WILL host
  it exists — `ScalarTraceSpace`" (shipped host was a NEW ladder tier the live docstring explicitly
  distinguishes from it) and "that is phase B5" (landed at G6.3, by factoring not by the predicted
  `u⊗v` typing). Preserve the prediction and `.. note::` WHY it didn't hold — more informative than
  the corrected sentence alone. → L-049
- **A stale-status blast radius is the WHOLE page.** Grep every future-tense/blocked token
  (`blocked|not built|not yet|pending|in flight|future seam|lands with`) — a brief naming 3 sites
  had 7. A "the one remaining unbuilt X" sentence must RE-POINT to the still-unbuilt sibling when X
  lands, not just drop X. → L-037

---

## 5. Page surgery: slice programmatically, assert before writing

- **Never hand-retype a large block.** Read → slice → `"".join` → write, with guard-asserts on the
  LIVE file's boundary strings and lengths, and ALL structural asserts run on the in-memory result
  BEFORE any write (a failed assert then leaves the tree untouched — no `git checkout` recovery
  needed). A machine splice cannot mis-transcribe. → L-012, L-022, L-023, L-026
- **⚠ An f-string mangles LaTeX braces in the header YOU author, into valid-but-wrong LaTeX `-W`
  never sees** (`A^{-1}` → `A^-1`). Strongest defense: author heads/intros/pointers as pure literals
  via the Write tool so no Python string layer touches math; then grep for bare `^-1`. → L-026
- **Locate by STABLE TITLE, never by the brief's line numbers — and prove contiguity by counting ALL
  H1 underlines in the range, not just anchored ones.** An ANCHORLESS sibling H1 (typically a prior
  split's leftover `:doc:` pointer stub) is invisible to the anchor-grep the brief's author used, so
  the range overshoots. Endemic to multi-split campaigns. → L-026
- **A cross-page theory MOVE is ref-safe if you KEEP the labels** (move, don't copy — defined exactly
  once). Only `:doc:` needs fixing (toctree + every `:doc:old`), plus now-intra-page
  `(:doc:sibling)` parentheticals, which lie without warning. → L-022, L-026
- **Splice mechanics that broke builds:** a slice ending in content, joined directly before the next
  `.. _anchor:`, GLUES the anchor to the preceding paragraph — the label silently fails to register
  and referrers report "undefined label" (NOT "duplicate") though grep shows it at column 0; join
  parts with `\n\n`. Re-nesting under a deeper parent demotes every migrated underline one level,
  LENGTH-PRESERVING (detect an underline as a col-0 all-one-marker line whose previous col-0 line is
  a plain title). Removing a middle H1 while keeping its trailing H2 AUTO-REPARENTS the H2 — verify
  no title-level SKIP and that the new parent is intended. → L-022, L-026
- **When the brief says "relocate to page X" and the CLOSE READ shows the content is already
  canonical on X, the action is DE-DUPLICATION (Cardinal Rule 2), not relocate+merge.** Replace with
  a `:doc:` pointer preserving the conceptual bridge, merge nothing, FLAG the inversion (the brief
  was built on a partial read). A FOLD is a MOVE — it re-exposes every symbol reconciliation the
  source had. → L-027
- **Prefer an additive prose ROADMAP + `:ref:` to the SSOT over copying a table (the twin) or
  relocating a double-duty section** — and first verify the gap is REAL (the presence of *a*
  taxonomy is not the presence of *this* one). → L-029
- **Metadata relocation, not deletion:** strip campaign provenance (hashes, phase labels, dates)
  from a high-traffic invariant section INTO the changelog, KEEPING invariants, eq-labels +
  vv-status, active gotchas, issue-`#N` refs and numerical data. Map each item to a destination
  FIRST — a dated milestone with NO changelog home keeps its provenance inline, flagged. → L-028
- **An overloaded-symbol sweep: inventory every MEANING of the letter first, classify each site
  mathematically, `replace_all` only unambiguous multi-char strings (enumerate spacing variants),
  targeted-edit the rest, then re-classify EVERY survivor in a final grep.** Flag same-letter
  collisions rather than renaming out of scope. A NEW page assembled from multiple sources is the
  prime site for a WITHIN-document collision the build cannot see — hunt every reused glyph and
  subscript the rarer meaning. A section RENAME can also stale an in-file back-reference: grep the
  file for the OLD heading text after a retitle. → L-011, L-025, L-034

---

## 6. Match the doc SHAPE to the event class

- **Stub → rich narrative: read memo → production docstrings → tests → SymPy, in that order.** The
  docstrings are the VERBATIM prose seed; the memo carries the honest interim scope — preserve it,
  don't over-claim. Never expand a stub without reading the SymPy; on an algebra error
  DISPATCH_REQUEST the method-implementer, never edit it yourself. → L-005
- **Campaign capstone (a completed feature's whole story): roadmap → literature memo → algebra of
  record → production code → error catalog → evidence pack.** The memo NAVIGATES; the SymPy module
  and production code are the CORRECTNESS spine. Arc: motivation → derivation-of-record → design →
  discoveries → evidence → honest scope. → L-039
- **A fix that works BY RETIRING a failed-approach family gets a SUCCESS-RESOLUTION chapter, not the
  9-step CLOSED arc**, and treats the superseded saga PROPORTIONATELY: ONE loud `.. attention::`
  supersession banner at the arc head + targeted tombstones on the bald factual REVERSALS only.
  Don't tombstone every stale sentence; don't rewrite the history. Flip any prior close-out's "open
  research path" that LANDED. → L-013
- **Deepening an already-documented feature: add the WHY, cross-link the WHAT, never duplicate.** A
  PLANNED refactor on a current-truth page gets a loud `PLANNED, not built` admonition (literals for
  unbuilt types) PAIRED with a "current state" subsection so plan and reality never blur. Verify
  every count against live code — a scoped grep undercounts; prefer describing the guard FAMILY over
  a call-site count. → L-014
- **An EVICTION changes the CARRIER, not the PHYSICS** — grep the stale carrier framing, not the
  physics terms (the brief over-counts; most of the chapter survived). Reframe narrowly plus ONE
  end-state paragraph cross-linking the new algebra. → L-016
- **A completed architecture earns ONE new taxonomy-CULMINATING section, opening by naming the
  generalization and `:ref:`-linking what it generalizes** — never a bolted-on appendix. Document
  symbol OVERLOADS as explicit gotchas. → L-018
- **A NEW foundational chapter: READ then RUN the algebra-of-record module(s) — one per distinct
  concept — before writing one equation.** Generalize by stripping the method-specific
  specialization while quoting the identity verbatim; RUN every load-bearing worked number through
  live code so the example is verified, not plausible. → L-025
- **Growing a thin honest-stub chapter at campaign close:** flip its own stale "in flight" status
  (verified against git, not its frozen prose), PRESERVE the already-landed section byte-for-byte
  and grow AROUND it, and RECONCILE sibling taxonomies explicitly with subset relations rather than
  contradicting them. Report DEFERRED WIRING with exact `test node → label` ids when tests you may
  not edit await labels you mint. → L-036
- **"Is the terminal docs phase done?" almost always answers effectively-DONE** — each earlier
  phase's doc pass landed its slice into the eventual capstone. Verify by the page's own
  SELF-IDENTIFICATION plus the build plus a cross-ref grep-gate. A documented SEAM is the OPPOSITE
  of a gap; and a charter's literal "the X page" can be correctly delivered as a SECTION of a shared
  page (a standalone page would MINT a twin path). → L-038
- **Merging a re-staged branch's docs into a diverged tree:** read the fork-diff as a CONTENT
  source, not a patch; splice programmatically; translate EVERY module path to the live layout (a
  moved package vs a same-named unmoved one) with zero residual; place by anchor, never by the
  diff's line numbers; reconcile forward-refs landing in the SAME merge. → L-012

---

## 7. V&V vocabulary — you are the curator (Directive 5)

You write the prose future readers QUOTE about verification status. Match `vv-principles` verbatim;
never paraphrase a level definition. → L-010

- **Never** "MMS verifies the eigenvalue" (source-driven: flux-shape / convergence-order only) ·
  never "L4 proves correctness" (name its L0–L2 backing) · never "the 1-group test verifies the
  solver". NAME the pillar (closed-form / MMS / semi-analytical), not vaguely "analytical". → L-010
- **Never upgrade a `@pytest.mark.foundation` gate to an L-level in prose to make a section sound
  better-verified** — read the marks and say "software/structural invariant of a discrete
  construction, not an equation claim". → L-040
- **A doc sentence "gates X, Y pin claim C" IS a coverage claim** — the prose analogue of
  `vv-principles`' "a `catches` marker is a COVERAGE CLAIM, not a topic tag". Justify it by a
  MUTATION that reddens X and Y, never by topical adjacency. Cite **per field**, not per topic
  (5 arrays needed 5 different files; one had a SOLE catcher, another was cylindrical-only).
  Highest-risk moment is REPLACING a gate you just demoted — the nearest-sounding sibling
  inherits neither scope. I credited a τ gate for reduced-operator arrays it passes in 0.03 s
  under fully-garbaged factories, two screens after writing the note explaining τ had LEFT that
  operator. → L-047
- **The SAME gate cited for TWO claims can be right once and wrong once — narrow, never sweep.**
  On "citation of X is false", ask *false for WHICH claim* and grep every occurrence before
  editing any; a blanket fix destroys the true citation. → L-047
- **Distinguish the EUCLIDEAN transpose `Aᵀ` from the metric HILBERT adjoint `A† = G⁻¹AᵀG`.** A
  campaign may colloquially call the former "†"; write the precise object. A docstring summary
  saying "Hilbert transpose" over a body computing the Euclidean one is a real defect. → L-010,
  L-034
- **A Mode-10 sub-floor term is closed by STRUCTURAL teeth, not a tightened value band** — and when
  no isolating regime exists, SAY there is no value-improvement leg to add. Pair the honest-scope
  note with a prophylactic `.. warning::`: the test pins the math, the warning pins the LANGUAGE.
  → L-010
- **Get a Mode-12 blindness boundary EXACTLY right** — a `k`-row is blind to the
  factor-order/similarity family, to all vector content, and to the metric itself, but NOT to a
  single-leaf transpose drop. Pair every `k`-claim with its vector/pairing catcher; a METRIC-REPAIR
  closure needs its CONTROL leg described, or a still-broken baseline mimics "caught". → L-036, L-015
- **The FIRST iterative member of a previously all-exact family has no bit-id twin to inherit** —
  claim foundation / flux-shape against a structurally-independent reference, and EXCLUDE the
  round-trip tautology explicitly so a reader doesn't mistake it for coverage. Related: never fuse
  an eigenvalue and a fixed-source term in one equation. → L-010, L-036
- **Skill-uplift duty:** propose the `vv-principles` / `error_catalog.md` / `algebra-of-record` edit
  in your return whenever you meet a published-prose anti-pattern or evidence-boundary case the
  skill doesn't capture. The skill grows when you feed it back. → L-010

---

## 8. Code-prose rebalance (docstring/comment trimming)

- **Expect ZERO MOVED.** Cardinal Rule 3 means the theory shipped WITH the code, so a concept that
  FEELS unique to the file is almost always already TWIN in the landing chapter. Grep the landing
  chapter before crediting one MOVED; a pre-classifier's MOVED column is ~100 % noise. → L-033
- **The CONTRACT test: "would a competent modifier who never leaves this file do the wrong thing
  without this line?"** If yes it is CONTRACT however history-flavored — including a keep-vs-retire
  decision on an intentional orphan, a ⚠ latent-trap imperative (keep the imperative + the
  falsifying number inline, derivation to a `§`-pointer), and a type-annotation rationale guarding a
  plausible wrong "simplification". → L-033
- **FILE-CLASS sets the size and SURFACE of the honest cut.** Teaching-heavy operator ⟹ aggressive
  TWIN-cut. Contract-heavy operator / machinery / ABC ⟹ small cut; the surface is module-head
  essays, campaign provenance and duplicated numbers, NOT method bodies. Driver / mesh ⟹ hunt
  standalone `#`-COMMENT tombstones and campaign-status blocks FIRST (comments dwarf docstrings). A
  −2 to −5 % cut is CORRECT — report the file-class rationale so it isn't read as timidity. → L-034
- **Provenance trimming = citation-vs-narration, applied uniformly** (trim landed campaign-STEP
  codes; KEEP bare `#NNN` anchors and named patterns; half-stripping violates internal
  consistency). But a hand-transposed-adjoint / reverse-scan comment body IS the algebra of record
  — KEEP it though it reads like narration. → L-034
- **A batch "special" is a VERIFICATION obligation first, an edit obligation only on failure** —
  read the oracle, read both ends, report SATISFIED, never touch a correct CONTRACT block. → L-034
- **Prove the edit is doc-only by AST/token comparison vs HEAD** (`tokenize` dropping
  COMMENT/STRING, or `ast.dump` with docstrings stripped), not by reading the diff. The AST check
  also proves no `verifies`/`catches` marker moved — but it is BLIND to comments (fine, they are
  editable) and an **f-string assertion message is CODE**: leave it and REPORT. → L-041, L-045 Run the Sphinx
  gate **iff** the file is `automodule`'d — `:noindex:` does NOT exempt it. A RENDERED file affords
  two extra moves: promote a latent trap to a real `.. warning::`, and repoint in-file back-refs
  after a heading rename. → L-033, L-034, L-041

---

## 9. Gates, generated artefacts, tooling

- **Generated artefacts are NEVER hand-edited** (V&V matrix, capability tables,
  `_generated/*.inc.rst`) — fix the registry-side metadata and report the REAL post-regen number. A
  `-E` build on a dirty branch absorbs rows from OTHER uncommitted work: a legitimate by-product —
  never revert it, REPORT it. In a fresh worktree a missing generated artefact is an ENV gap, not a
  doc defect: materialize it, never route the docs around it. Never transcribe a hard-coded test
  count. And orphaned built HTML from a renamed source looks like a live stale ref — discriminate by
  "does the source `.rst` still exist?". → L-008, L-040, L-026
- **RE-MEASURE the `-E` baseline every session; never assume a recorded number** (it has drifted
  9 → 1 → 0). Diff the WARNING/ERROR/CRITICAL SET pre/post, not just the count. A full `-E` rebuild
  can exceed the 120 s foreground cap — background it, or the poll-loop is SIGTERM'd at the final
  line. → L-029, L-041, L-027
- **Title markers (AGENT.md owns the rule): prefer COPYING a proven underline from a model page
  over re-counting code points.** → L-009, L-035
- **Citations: `grep '^\.\. \[Key\]'` before citing or defining** — resolve cross-doc, never
  redefine; match a page's plain-text convention; on a strict-zero-warning NEW page go ALL
  plain-text (a Literature `list-table` with equation numbers inline — higher articulation, zero
  machinery). Pre-existing duplicate-citation warnings are a known trade-off: verify the count is
  unchanged. FLAG a conflated bib key. → L-006, L-025
- **A corpus-wide mechanical migration is dry-run-first and WHITELIST-scoped** (which auto-skips
  every pseudo-site); key any block remover to INDENTATION too, or it eats footnotes. → L-031
- **Self-check the V&V scan directly, not via the full audit** — the theory-equation scanner runs in
  <1 s, avoids pytest collection, and doesn't trip on sibling batches' in-progress sentinels.
  → L-035
- **zsh does NOT word-split an unquoted `$var`** — a uniqueness loop ran once on the concatenated
  string and printed a false "0 collisions". An `Edit` `old_string` must match LIVE bytes. → L-030
- **An error-message string inside `raise` is EXECUTABLE — report it, don't edit it, under a doc-only
  constraint** (tests `pytest.raises` match on those strings). Same for a brief item that is a LATER
  phase's acceptance-gate text: leave it, name the owning phase. → L-041

---

## Quality self-assessment (Directive 3)

Rate 1–5 and log the weakest dimension: Derivation depth · Cross-references · Numerical evidence ·
Failed approaches · Code traceability · Derivation source. On TERMINOLOGY / ROUTING / retirement
passes the weak dimension is routinely "numerical evidence" — structurally ABSENT (no flux moves ⟹
no convergence table), not a deficit. Say so, don't manufacture one.

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

- **⭐ TWO AGREEING SOURCES CAN BOTH BE WRONG — corroboration is not independence.** A brief said
  a gate held *"dated 2026-08-21"* and an in-tree comment independently agreed; `[M]`
  `git log --date=short` puts every commit of that step on **2026-08-20** (a FUTURE date on both
  surfaces — one mis-dating copied forward). Git is the arbiter for dates and merge status; two
  prose surfaces agreeing is one surface. → L-064
- **⭐ Publish the CLOSURE ARGUMENT, not the universal it implies.** Instead of *"every harmonic
  space is legacy"*, write *"`of_axes` is the only ROOT producer of an `axes` record (`*` and
  `dual` merely THREAD one), therefore …"* — the derivation stays true as the tree grows, and the
  grep that establishes it is free reconnaissance (mine surfaced `mm.axes` = the GEOMETRIC tuple
  vs `mm.bulk_space.axes` = the SPACE-FACTOR tuple: one attribute NAME, one object, neither
  derived from the other — a publishable gotcha). → L-064
- **Read the live `def`/body before citing any convention, shape, signature or design decision.**
  Seen: a docstring lying about an index convention and a return layout; a verdict memo recording
  the RECOMMENDATION while the code shipped the alternative; a brief naming args the live Protocol
  never takes. → L-001
- **⭐⭐ A SPEC is its TABLE — re-derive the brief's counts from the artefact in one command.** A
  spec file's headline read "21 of 40 declarable, 19 NONE" over a table that is `[M]` **32/8**;
  the brief inherited "19" verbatim (its own kind breakdown summed to 8). A headline is a summary
  and summaries rot while the table stays right. → L-060
- **⭐⭐ A brief's "sharpest observation" is a HYPOTHESIS with a computable confusion matrix.** A
  briefed *"the page already labels its own two classes"* keyword-tell died on four greps: the
  word *identity* sits in 5 of 6 un-implementable rationales AND **11 of 22** implementable ones;
  *"not a solver claim"* points the WRONG way (1 of 6 vs 5 of 22); a third of the page carries no
  rationale at all. Publish the measured split AND the refutation, or the next reader re-derives
  the heuristic and ships it. The surviving distinction was real: an identity between
  **quantities** has no carrier, an identity between **types** is a claim about a class
  declaration. → L-060
- **⭐⭐ Every "each / every / all" in prose YOU publish is a universal you can count in one
  command — count it.** I broke this twice in one new section and caught both by re-measuring:
  "…and X, Y, Z with three each" enumerated 3 of 5 (and read as all 15 multi-implementer rows);
  "where every ``solve`` in the tree matches the label" was `[M]` **5 of 60**. The measured
  sentence was strictly better both times. (plan-authoring §2, applied to the corpus.) → L-060
- **⭐ Verify a SUCCESSOR claim against the RETIRING COMMIT'S BODY, not against the successor's
  existence.** A live class with the right shape is not evidence the dead one became it. `[M]` I
  wrote that three symbols "were absorbed into the operator algebra"; the commit bodies say ONE was
  re-layered (`SNSolver.L` → two leaves) and the other two were **retired outright** — one "became
  orphan in production", the other "without a remaining call site". One paragraph, two fates, and
  "absorbed" was false for two of three. `git log --diff-filter=D` / `-S <symbol>` then read `%b`.
  → L-062
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
- **Ask WHICH SIDE of a carve a brief's "measured" number came from — and re-measure BOTH.** A
  briefed "post-carve: agrees to solver tolerance, `1.998e-13`, `array_equal=False`" was a PRE-carve
  number (two structurally different deliveries reaching one fixed point); post-carve the two
  channels are the same float program and collapse to **`0.0`**, `array_equal=True`. Publishing the
  brief would have inverted the carve's headline AND justified an `rtol` gate blind to the very
  defect (`2.9e-14` sails through `10×inner_tol`). A pinned pre-carve worktree makes both sides
  cheap and turns one number into a before/after table. ⚠ the venv's **editable install hooks
  `sys.meta_path`, which OUTRANKS `sys.path`** — `PYTHONPATH=<worktree>` silently loads the MAIN
  tree; strip the editable finder and PRINT `orpheus.__file__` as proof. → L-050
- **Never accept a fixed-decimal printout as evidence of bit-exactness.** "`2.500000000000` at 12 dp"
  cannot resolve `8e-15` at 2.5. Measured: the converged inflow trace is exact on SI (the sweep
  *writes* the seed) and NOT on Krylov (GMRES *solves* the trace rows — 1–23 ULP at `tol=1e-13`,
  27 580 ULP at `1e-10`, i.e. the iteration residual, not FP noise). An exactness claim true on one
  inner solver is a **per-leg** gate (`array_equal` on the exact leg, `rtol=SAFETY×inner_tol` on the
  iterative one) — say "do not relax the exact leg to match". Assert `x == v` or print
  `float(x).hex()`. And run your OWN probes without `-O`: a bare `assert` in my widened check was
  stripped (vv Mode 8, in my own instrument). → L-050
- **⭐⭐ Before you PUBLISH a number a brief hands you, grep `tests/` for a module NAMED after
  the phenomenon — and read its docstring, not its `assert` lines.** A `scratch/` memo is by
  construction OLDER than the tests it motivated. A briefed "we just measured `min ψ̂ ≈ −77`
  under the shipped convention" was reproducible to the digit AND its framing was already
  refuted by a 19-row `foundation` module committed the SAME DAY: on the production
  (marched-seed) path ψ̂ is strictly POSITIVE (`+0.134/+0.129/+0.129`, 0.88–0.98 × `min ψ`);
  `−77` is an INCONSISTENT-seed statement. I had drafted "⚠ coverage gap: no ψ̂ gate on either
  arm" — false. My line-based grep missed the module because its evidence lives in the
  docstring and in `pytest.fail` messages, not in `assert` lines (vv #21). → L-055
- **⭐⭐ Reproduce the witness — the reproduction can REFUTE the GATE'S OWN PROSE while
  STRENGTHENING its claim.** A cone-violation gate's docstring said its two legs differ in
  *"ONE parameter … half the optical cell size"*; `[M]` both legs have `Δx·Σ_t = 100`
  **identically** (`nx=2,width=20` vs `nx=4,width=40`). The argument is stronger than its
  prose — holding cell size fixed kills the rival explanation outright — so publish the
  CORRECT framing, never weaken the claim, and REPORT the docstring (you don't edit `tests/`).
  ⭐ Then run the two scans nobody asked for: a cell-SIZE scan reproduced the textbook DD
  positivity limit exactly (`Δx·Σ_t = 1` in K, `= 2` already out) and a cell-COUNT scan showed
  the benign row is the only one — ~90 s turned one frozen constant into a mechanism. → L-063
- **⭐⭐ "Make illegal states unrepresentable" is TWO-sided, and half 2 is the one skipped.**
  Mint the invariant iff (1) every admitted value is legal AND (2) every legal value is
  admitted. Half 2 is a claim about the PRODUCERS, not the concept; when it fails the
  invariant does not prevent a bug, it **refuses correct output**. Quote this into any page
  invoking the pattern. → L-063
- **⭐ When a brief offers a BINARY verdict and the tree supports neither pole, publish the
  third.** "Scoped to the scalar flux (fine)" vs "general (falsified)" both missed it: the
  claim was general-in-wording, sphere-in-evidence, and substantively true on the cylinder's
  production path *as a characterisation*. Only the word "never" was false. Scope the heading,
  publish the seed taxonomy with both measured tables, keep the original conclusion standing
  where it survives, and point at the owning gate. → L-055
- **⭐⭐ A DESIGN PROBE goes stale against the repair it motivated, SILENTLY — it still
  runs and still prints plausible numbers.** `[M]` the plan's probe read
  `frame.test_space.inner_product_weights` as "the stored metric"; post-repair that IS
  the repaired metric, so the row labelled *stored* now prints `1.000` and the headline
  `118.7` is unreproducible from the file the plan cites. ⟹ never cite a pre-repair probe
  path as a post-repair page's reproducer — publish the **CONSTRUCTION** (build the object,
  name the seed, name the five attributes the residuals come off) so the table regenerates
  from the page; re-measure with your OWN probe; report the staleness upward. → L-065
- **⭐⭐ A SEED-DEPENDENT number is published as its BOUND, never its value — then find the
  exact per-mode parent behind it.** `[M]` plan `118.7` vs mine `81.4`, both correct: the
  Parseval ratio is a moment-energy-weighted average of the per-ℓ factors. Publish the
  draw-independent claim (*lies between the extreme factors PRESENT AT THAT L —
  `[17.5,157.9]` at L=1, `[6.3,157.9]` at L=2 — so never 1*), and note a bound is a
  universal (plan-§2: I first wrote the L=2 range over a sentence covering both). ⭐ Then
  look one level down: the RATIO OF THE TWO ADJOINTS on a single-ℓ unit input is exactly
  `(4π/(2ℓ+1))²`, `[M]` to `≤2.8e-16` — seed-free and strictly more useful. → L-065
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
- **⭐⭐ Publish YOUR number with YOUR configuration; never relay one whose fixture you cannot
  state.** Two brief/plan numbers did not reproduce because they were measured on fixtures I did
  not have (a heterogeneous SI-vs-Krylov gap; a pre-fusion build time). Re-measuring gave
  different, correct values on a fixture I could name — and the re-measurement also refuted the
  brief's headline *as stated*: "excited iff `n_x` is ODD" is true of the symmetric fixtures and
  false as a mesh property (`[M]` `dim ker A = 12` at even AND odd `n_x`; an anisotropic source
  excites the even one at `1.76e-2`). Publish the scoped rule + the counter-row. → L-057
- **⭐⭐ A ratio is a ratio OF AN OBSERVABLE — name it before citing it.** A memo's `n_GS/n_J` was
  ρ-DERIVED (`ln ρ_J/ln ρ_GS` from an eigen-solve) while every published table reported SWEEP
  COUNTS. My control reproduced the published `1631/838` **exactly**, then **4 of 5** memo rows
  disagreed in **SIGN** (`0.576` "wins" vs measured `2.503` "loses") — two individually sound
  instruments, different observables, near-degenerate rates. Publish only the observable you
  measured; never let a rate-ratio and a count-ratio share a column heading. → L-051
- **"Already done in code" ≠ "gating green".** A brief's "don't redo the code" premise said nothing
  about the build: the `-E` baseline carried an `ERROR: Malformed table` in the very docstring that
  pass had just edited, so the brief's own `-W` gate could not have passed. Fixing it was blocking
  AND in scope. Diagnose a simple table by rebuilding column spans from the `===` separator
  (`re.finditer(r'=+', sep)`) and flagging non-space chars in the gaps. → L-051
- **Importing algebra from code/memo imports its SYMBOL COLLISIONS** (`A_a` face-area vs the loss
  operator `A`; `Σ` transmission vs `Σ_t`). KEEP the code's spelling — internal consistency outranks
  local awkwardness — and pay with an explicit `.. note::` naming each overload AND its
  disambiguator. Never silently rename into the docs: that mints a code↔corpus twin. → L-051

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
  literals on two other pages). → L-045, L-046, L-047
- **⭐⭐ …but ON AN `.rst` PAGE that gate reports `:mod:` and NOTHING ELSE — never read an unmoved
  `DEAD TARGETS` as "my page is clean".** `[M]` I fixed 15 dead roles and the count sat at
  **81/124 both sides**: `judge()`'s head-check re-checks the target's head *carrying the original
  role*, so a single-segment head (`orpheus`) trips `bare_module_guess` under any non-`mod` role,
  and with a page's empty namespace the candidate tuple is `()` → DECLINED. One line fixes it
  (`head_role = "mod" if "." in target else role`); `[M]` on a pristine `git archive HEAD` tree that
  takes `docs/` from **49 dead/71 sites → 207/255** — the gate is blind to 158 dead
  `:class:`/`:func:`/`:meth:`/`:attr:` targets in `docs/` alone. Until it lands, the acceptance
  evidence for a page is YOUR OWN import probe over its roles. ⚠ Measure such a patch as a COPY run
  as a SUBPROCESS — monkeypatching `judge` and calling `main()` twice in-process read `0` for both
  arms while a subprocess read 49. → L-062
- **⭐⭐ A per-site adjudication TABLE is an instrument — audit its SKIP clause, its "retired"
  verdicts, and its `hasattr` evidence.** Applying a 91-site ruling table faithfully still yielded
  FIVE corrections: (a) its *keep-if-absent-from-graph* filter hid ~1400 alive-but-unqualified
  roles, incl. one on a line I had to edit — a NOT-clause is a false-NEGATIVE source; (b) 7
  "RETIRED" sites were a pure `git mv` (`--diff-filter=D` on the old path answers this in one
  command) and 6 of them read PRESENT/IMPERATIVE, so a literal would have killed a live link;
  (c) `hasattr(Cls,"mesh") is False` ruled a TRUE paragraph false — the base class sets it on the
  INSTANCE, so **construct the object**, never probe the class; (d) two adjacent `X.apply` sites
  had OPPOSITE right answers (one qualifies — the page already links it 6 lines up; one is dead —
  the type has no `apply`), so resolve `Instance.method` by the instance's TYPE, never by the
  target string's shape; (e) one table cell opened a page-wide SYMBOL-INDEX collision. → L-053
- **⭐ "Qualify so it resolves" is TWO claims that come apart — say which one you bought.** `[M]`
  post-build hrefs: `EigenvalueSolver` 43, `Field` 30, `numpy.array_equal` 6 — real links;
  `KEigenvalue`/`SNMesh.axes`/`BC.vacuum`/`zeros_on`/`peierls_nystrom.slab` **0**, still plain
  text (`:noindex:` autoclass, or no automodule). Import-/graph-resolvable ≠ rendered link; check
  with `grep -o 'href="[^"]*#<target>"'` in the built HTML. → L-053
- **⭐⭐ A page can run a symbol convention ONE INDEX BELOW the code's, and no test can see it.**
  CP's `Ki_4` IS the shipped standard-`Ki_3` (`[M]` `_ki3_mp(0)=0.7853961=π/4`), its `ki3-def` is
  the standard `Ki_2`, its `F(0)=0.4244=4/(3π)` matches neither — across 3 labels carrying
  64/24/54 `verifies()` tests that pin the code's NUMBERS and are structurally blind to the
  equation's NAME (a doc-side Mode 12). Repoint the role, MEASURE the collision into an anchored
  `.. warning::`, fix only the unambiguously-wrong number, and REFUSE the ~30-site re-indexing as
  a numerics adjudication — a physics rename must not ride inside an xref pass. → L-053
- **⭐ NEVER say "all trees at 0" — say the TREES, the ROOTS and the SEMANTICS.** The gate walks
  only `orpheus tests docs` (NOT `examples/`, top-level `derivations/`, `scratch/`, `tools/`) and
  judges only `DECIDABLE_ROOTS` (which omits importable `tests`/`tools`/`derivations`); its NAME
  also understates it — it DOES read whole `.rst`, so `doc:` sites ARE covered. Decisive: the gate
  resolves by **IMPORT**, the `nexus dead-references` hook by **RENDERED TARGET**, so a live
  un-`automodule`'d module is *resolved* to one and *dead* to the other. Both right — the SET
  DIFFERENCE is the triage (hook-only ⇒ un-surfaced-but-live = #302; both ⇒ really moved/retired).
  Measured 21/30 vs 0/14 914 on the same tree. → L-052
- **Two gate false-negatives to know:** a **PEP-420 namespace package** (a dir with only a
  `README.md`) IMPORTS fine (`__file__ is None`, 0 members) so `:mod:` refs at it read "resolved"
  though Sphinx can never link them; and a role **wrapped INSIDE its dotted path**
  (`:func:`~a.b.\n  c``) is skipped by the gate (whitespace ⇒ `extract_target` → `None`) AND
  renders plain text. 15 such roles tree-wide; the discriminator is `\.\s*\n\s*\w` in the pre-`<`
  head — ~180 multi-line roles that break at the `display <target>` seam are FINE. → L-052
- **Before believing a dead target's NAME, read its graph EDGES**
  (`SELECT source,type FROM edges WHERE target=?` on `graph.db`). A name can be an artifact minted
  by a THIRD tree: six `orpheus.derivations.peierls_geometry.*` targets existed only because
  `scratch/` scripts still import the deleted path, and nexus suffix-matched the theory pages'
  unqualified roles onto it. Edge type decodes the site class: `documents`=page ·
  `references`=docstring · `type_uses`=**a type annotation, i.e. a CODE bug** · `calls`=the import
  that minted it. And nexus counts doc sites **per PAGE** (2 "sites" was 9 roles), while the unit
  of repair is the TARGET tree-wide (3 reported sites, 13 real). → L-052
- **⭐ Two defects YOUR OWN new prose introduces, both `-W`-caught, both mechanical:** (a) an
  italic run interrupted by a role — `*"… (*:math:`X`*) …"*` — gives *"Inline interpreted text …
  start-string without end-string"*; escape the seam: `(*\ :math:`X`\ *)`. (b) **NEVER hand-align a
  simple `===` table containing a `:math:` role** — the role's SOURCE length is not its rendered
  length, so the column arithmetic is invisible and you get `ERROR: Malformed table`; use
  `list-table`. → L-054
- **Beyond AGENT.md's warn-list, two more DO warn:** a `:widths:`/column mismatch, and `ref.ref`
  "*A title or caption not found*" — a bare `:ref:` to an anchor sitting before a PARAGRAPH,
  **including a BOLD RUN-IN HEADING** (`**(2) Title.**  Prose…` looks like a heading and is
  not one; ONE such anchor cost 5 warnings across 4 files and EXIT=1). Fix: promote to a real
  titled subsection — but never open the title with `(1)`/`(2)`, an RST enumerated-list
  marker; use `Correction 1 — …` — or use explicit text `` :ref:`text <label>` `` when the
  anchor legitimately sits above an admonition. ⚠ **`check_docstring_xrefs.py` is BLIND to
  this class** — it gates Python-domain roles, and reported `DEAD TARGETS: 0` while all 5
  were live; only the `-E -W` build sees them. Raw path strings in prose warn NOWHERE.
  → L-002, L-027, L-040, L-055
- **Grep `SyntaxWarning` in the build log too — a case-sensitive `WARNING:|ERROR:|CRITICAL:` MISSES
  it** and it does not bump the exit code. A non-raw docstring with `\Gamma` emits
  `SyntaxWarning: "\G" is an invalid escape sequence` mid-build. Before reporting one in a file
  another agent is editing, `py_compile` the LIVE file — mine was fixed a minute later. → L-048
- **Grep a glyph in `docs/` before importing a marker from a plan — and re-grep, because this
  answer MOVED.** `[M]` 2026-08-14: `⛔` 12 files, `⚠` 17, `⭐` 10, `✓` 10, `✗` 2. All are corpus
  vocabulary now; the old "⭐/⛔ are zero in docs/" reading is retired. → L-048, L-056
- **⭐⭐ A mechanical PORT's warning count is a non-representative sample of its defect count —
  census the target language's delimiter alphabet before fixing warning #1.** RST has no legal
  run of 3+ backticks outside a literal block, so a run-length histogram is a TOTAL census: a
  briefed "handled the bulk correctly, 20 warnings left" was `[M]` **830 mangled delimiters on
  339 lines** + 30 broken spans; the warnings were the ~2 % where the imbalance failed to cancel
  inside a paragraph. One root cause — a LINE-LOCAL `` `x` ``→``` ``x`` ``` regex meeting a code
  span that WRAPS a line — with three surfaces: over-added pair (silent stray backticks),
  one-side-only (warns), neither-side (silent `<cite>`). **The port's SOURCE is the oracle and
  makes the mass edit a PROOF**: post-fix 2438/2443 literal contents and 3648/3653 prose lines
  appear VERBATIM in the `.md`, every exception explained. → L-061
- **⭐⭐ `<cite>` in the built HTML is the smoking gun of a Markdown port** (`default_role` unset
  ⟹ a surviving single-backtick span renders ITALIC, not monospace) — and its sibling: **RST
  forbids inline markup after most chars** (`= . ~ § ↔ *`), which Markdown does not, so a port
  leaves roles that DO NOT RENDER AT ALL. `[M]` two `:math:` roles opening after `~` produced
  `<cite>mathcal{O}(h^{1.3})</cite>` — role dead, LaTeX backslash eaten, **silent at every build
  severity**. Fix is one char (`~\ :math:`…``). Census with `grep -c '<cite>'`. And `\|` added by
  a port is RIGHT in prose, WRONG inside a literal (renders a visible backslash). → L-061
- **⭐ A dead `:doc:` from a Markdown port is usually a PATH-FORM error, not a missing page** —
  MD authors write repo-root paths (`docs/theory/…/index`), Sphinx wants a docname
  (`/theory/…/index`). Check the page EXISTS before rewriting the prose; the warning and the
  brief both read as "points at nothing". → L-061
- **⭐⭐ The WHERE-LIST is the tell that a labelled equation drifted from its own prose.** A
  definition list defining symbols ABSENT from the equation it annotates, and omitting one that is
  IN it, is a correction that stopped one line short. Needs no code and no build — and no build can
  help: the label EXISTS, so every `:eq:` resolves, `-W`/`-n` are silent, and the V&V matrix reports
  the label covered because coverage is keyed to the LABEL, not to what it says. Seen: a 2026-08-02
  fix rewrote the geometry table, the worked examples, the rejection messages AND the predicate
  quoted in the equation's own vv-status comment — leaving the `.. math::` body stating the retired
  claim 8 lines above its own correction. Publish the tell IN the page. → L-056
- **⭐⭐ A `.. implements::`/`.. verifies::` DECLARATION is a doc surface whose failure mode is
  INCOMPLETENESS, and declaring is the only thing that stands the guessing down.** Nexus infers
  code↔equation links from shared name TOKENS; declaring ONE implementer switches inference off
  for the WHOLE equation, so an equation declared with 1 of its 2 implementers shows only 1 —
  under-coverage produced by an act that looks like an improvement. Ask *what else computes this?*
  before the first directive (7 of 14 needed 2–4: DD arm + LD arm; forward + transpose; the scheme
  door + the schedule that folds its term). Pre-flight every `:by:` path and label against
  `graph.db` before writing (a bad `:by:` DOES warn, so `-W` gates paths — but only a
  post-rebuild `edges WHERE type='implements'` query proves the COUNT). MEASURE the instrument you
  displace: `[M]` 397 guesses vs 28 true implementers on one page, **21 % recall, 1.5 % precision**,
  and the guess sets for two unrelated equations **identical 23-for-23** — because the matched token
  was the PACKAGE name, so the set is a module membership list that *cannot contain* implementers
  living elsewhere (`loss-rep-LpC`: **0 of 23**). → L-059
- **⭐⭐ Writing the explanation MINTS new guesses — citing a symbol to say it is NOT the
  implementer makes the heuristic name it as one.** `[M]` adding two `:meth:` xrefs while
  explaining three undeclared equations raised their guess counts 23→24/25/24. ⟹ NEVER publish a
  live guess count (quote the frozen pre-declaration measurement or say "re-run"), and know that an
  undeclared equation gets WORSE every time its page is improved. → L-059
- **⭐⭐ "All N node IDs resolve" ≠ "all N `:by:` targets bind" — the DIRECTIVE's resolver is
  NARROWER than the graph.** `_node_id_for_target` tries the literal string, then
  `py:function:`/`py:method:`/`py:class:` and **nothing else**, so a `TypeVar` (a `py:data:` node)
  fails a bare dotted name while a graph pre-flight says "exists". Fix: pass the already-prefixed
  node id (`:by: py:data:pkg.mod.Domain`) and write the reason into the page. Pre-flight through
  the RESOLVER's prefix list, not by node existence. → L-060
- **⭐⭐ An equation with NO implementer keeps its guesses forever — you cannot declare an
  absence.** `[M]` post-pass on one page: 57 directive / **0** inferred on the 32 declared; **166**
  inferred remaining, every one on the 8 that cannot be declared (60 on `operator-solve` alone).
  That residue is the ceiling on what authoring can retire, and it is the argument for a
  machine-readable KIND. → L-060
- **⭐⭐ `no-implementation` has a class that LOOKS declarable: an IDENTITY BETWEEN TWO
  QUANTITIES THAT ARE EACH COMPUTED.** `φ = Mψ = Gc` — `Mψ` is the analysis face, `Gc` is
  `discrete_gram`, and the identity is evaluated nowhere; same for `d_ℓ·G_ℓ = W`, whose
  two factors ship and whose product is never formed (that IS the point — it lets the
  kernel carry one `1/W` scalar). Declaring either side asserts that one of them IS the
  identity. Use `:kind: identity`, say which symbol computes which SIDE, and name what the
  suite measures instead (the CONSEQUENCE). `[M]` 17+16 guesses → 0 on two labels, and a
  *contrast* label went 2 → 0, one of its guesses being an SN solver entry point that never
  touches the faces it names. ⭐ Mirror: re-deriving one equation's `implements::` set is
  where you find another's — an implementer that LEAVES one equation has to LAND somewhere
  (`metric_per_ell` left the adjoint, gained `sh-space-metric`, 3 implementers, 0 declared
  before). → L-065
- **⭐ The `documented` sentinel marks the KIND, not the coverage — a label can honestly sit
  in BOTH matrix lists.** `[M]` `hilbert-adjoint-…` is verified by 9 tests AND sentineled;
  so are its page siblings. On this corpus `documented` = representational /
  face-distinction / literature KIND. Do NOT "clean up" the redundancy from one label — that
  re-categorises a convention and moves a generated artefact. Keep the directive; ADD the
  rationale comment if missing. → L-065
- **⭐ "Implemented by nothing" is a CLASSIFICATION worth a section — and the classes are
  enumerable, which is what an inference cannot know.** `{identity, law, canonical-form}` → NONE;
  `{typing-rule, definition}` → look for a declaration site (a typing rule CAN have a materialized
  carrier — a class, a Protocol parameter list, typed methods; an identity cannot). Three further
  classes seen: **superseded path** (the
  identity stays TRUE, the code that walked it is gone — check for INDEPENDENT retirements, there
  were two) · **notation** (a definition whose arithmetic IS computed elsewhere, for a *different
  operator* — declaring it is a false attribution at `confidence=1.0`) · **declared tag** (a
  `ClassVar[bool]` a human set after doing the math by hand; the implementer of the *criterion* is
  the page). → L-059
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

- **⭐⭐ An ONTOLOGY OVERTURN is not a retirement sweep — grep the retired symbol to FIND the
  sites, then read the enclosing ARGUMENT to decide the edit.** The dead-ref half finishes in
  one pass; the load-bearing half has no dead symbol in it. A five-obstruction proof that
  `Carrier[Rep, Role]` is impossible rested on *"the Flux role must make `flux + flux`
  raise"* — premise now false, **conclusion still true**. Deleting it destroys a correct
  proof; leaving it ships a false premise. ⟹ re-derive from what survives, keep the
  conclusion, tombstone the example — and the live tree hands you the replacement: `[M]`
  `AngularFlux` now defines NO `__add__`/`__sub__` while `AngularSourceSink` does, so the
  "changes the arithmetic interface" axis **inverted**. Same shape hit 5 sites on one page.
  → L-063
- **⭐⭐ An eq-LABEL is RETIRED when its NAME encodes the refuted concept and RENAMED when
  only its ADJECTIVE is stale — the discriminator is the label's BODY.** 4 labels, 4 fates:
  body states the retired claim ⟹ retire + repoint every `:eq:` citer; body still true ⟹
  rename to a live name; body untouched by the overturn ⟹ **KEEP + `.. note::` at the anchor
  saying the prefix is a historical artefact** (that one had **8 cross-doc `:ref:` citers**,
  and a cross-doc dangling `:ref:` renders plain text at every severity — renaming buys
  cosmetics and risks a silent break). A stale NAME is not a false CLAIM. ⭐ The retired
  equation still gets SHOWN in the history section — as an **UNLABELLED `.. math::`** with
  one line saying why, so it cannot become an `:eq:` API by accident. → L-063
- **⭐ Check the vv-status bookkeeping in 1 s, and never hand-edit the matrix.** A sentinel
  must name a `:label:` in the SAME file (`tests/_harness/audit.py`), so a rename without its
  sentinel is a hard violation: `from tests._harness.audit import _scan_theory_equations`
  → `.violations` / `.documented`, sub-second, no pytest. `matrix.rst` regenerates at
  `builder-inited`; report the post-regen count (`[M]` 539 → 540). → L-063
- **⭐⭐ Two SKILL files outside a brief's scope can carry the retired ontology — and your own
  repair can IMPORT the falsehood through its cross-reference.** #18's corrected text points
  at `cross-domain-frames` A.1, whose worked example is the retired type. ⟹ flag the
  staleness **inline at the pointer** ("A.1's frame is sound; its example is NOT"), report the
  out-of-scope files, never silently edit them. ⭐ And a REVERSED anti-pattern leads with what
  survives, carries the falsified version verbatim beneath, and ships the checkable test the
  reversal yields. → L-063

- **⭐⭐ When the corpus states ONE object N incompatible ways and each is internally
  consistent, that is not N bugs — a hidden PARAMETER is unnamed.** `[M]` three published
  `Π*` (naked `S₀` / `g_C·S₀` / `S₀∘G⁻¹`), plus one admonition whose EQUATION and PROSE
  disagreed with each other, warning-free for months. All three are the correct adjoint
  under a different coefficient metric; none named its metric. ⟹ do NOT adjudicate — name
  the parameter ONCE in a `list-table` at the point of definition (metric | where it lives
  | the adjoint it induces), then make every site a POINTER into one row, so none can rot
  independently again. The tell is free: two defended statements of the same object that
  disagree. → L-065
- **⭐⭐ The reusable close-out shape for a LATENT defect is THREE shields, and shield 3 is
  the dangerous sentence.** (1) *Consistency is not correctness* — the defining identity
  held at the round-off floor because `.H` is BUILT FROM the stored metric, so it is true
  for every SPD metric and carries ZERO information about which is installed; the
  instrument that can fail compares the metric to something defined without it. (2)
  *Composed chains are immune* — interior metrics cancel. (3) *No end-of-chain consumer
  existed* (`[M]` one grep hit, a docstring). ⛔ Write (3) as **latency, with the clock**
  ("live with the first adjoint consumer, which is why the metric had to be right before
  those land") — reported as reassurance it teaches the next session to defer. → L-065
- **⭐ Extending an ERR entry vs minting a new number: the LANDED MARKERS decide, not the
  narrative.** F-0 became ERR-039's third chapter because the shipped gates already carry
  `catches("ERR-039")` and I cannot edit `tests/` — a new number would silently orphan
  them. Read the catching tests' markers BEFORE choosing. And mark the superseded chapter
  IN PLACE, on the bullet stating the retired formula, not only in the new chapter. → L-065
- **⭐ A three-way SYMBOL COLLISION: rename only the one with NO constituency.** `W` was the
  coefficient space (page convention), the quadrature metric subscript (page convention),
  and the scalar total weight (code + ledger); my derivation needed a fourth, the weight
  MATRIX — the only one with no constituency. Write it `\mathrm{diag}(w)`, keep the other
  three, and open the section with a `.. warning::` naming all three survivors. → L-065, L-051
- **⭐⭐ Before repairing a stale equation, census the CORPUS for a page that already states it
  right — the census does two jobs and both are load-bearing.** (a) It stops the repair minting a
  TWIN: my `keff-as-integrated-rates` fix restated a formula whose SSOT already ships as
  `:eq:`sn-keff-update`` under `:ref:`sn-keff-estimator`` *with* the derivation and the gate — the
  correct shape is keep the equation (a label is an API and must not be false) **+** an
  `.. important::` naming the SSOT and saying which claim THIS page owns, plus "edited there,
  mirrored here". (b) It tells you whether you are fixing an outlier or inventing a convention: on
  the `S → Λ` repair, `[M]` **3 sibling pages + the class docstring already wrote `Λ`** and the
  edited page was the sole holdout — that census IS the evidence the repair is right. Run it
  BEFORE drafting. ⚠ It also surfaces symbol collisions across pages (SSOT writes leakage `L`,
  my page writes `L = Ω·∇`): resolve with a local subscript **plus a note naming both
  spellings** — never silently. → L-060

- **⭐⭐ A brief's SITE CENSUS is a sample; run the windowed CONCEPT grep yourself before fixing
  one.** `[M]` a brief naming 1 + 6 "may be inherited" sites: a regex for the predicate within ±4
  lines of the subject found **18 in one file**, all present-tense-false the same way. Leaving 12 is
  the exact half-fix vv #21 warns about. ⚠ The tell that the file already knew better: the
  CORRECTED framing sat **one line above** a stale sibling docstring. And the same falsehood was in
  the RST too — four places, one of them the page's own **Key Facts** card and one the prose
  wrapping the very equation I was declaring against. Fix shape: keep the equation (the identity is
  TRUE), keep the label (live `verifies()`), correct only the ATTRIBUTION, tombstone the mechanism
  in past tense naming the carve. Retitling is safe iff the section carries no `.. _anchor:` — and
  tree hits in `_build/` are orphaned HTML, not references. → L-059
- **⭐ A ⛔ ruling's quantifier needs YOUR census — the brief greps the NAME, and the sharpest sites
  never bind it.** `[M]` "the only two sites forming `σ_t − σ_s0`" was **four**: one in PRODUCTION
  (inline, `# ¼ σ̂_R h`) and one an MMS **capture** cross-section — identical arithmetic, a material
  datum, coincident only at 1 group. The ruling survived and got stronger; the count did not.
  ⚠ Mirror trap: `sig_r` in `thermal_hydraulics/`+`kinetics/` is a RADIAL STRESS — a short suffix
  collides as badly as a one-letter symbol. → L-059
- **⭐⭐ A math symbol has THREE spellings — ASCII id (`tau_raw`), Unicode prose (`τ_raw`), LaTeX
  role body (`tau_{\rm raw}` / `tau^{\rm raw}`) — and a brief's page list is built from ONE.** So
  it is over- AND under-counted at once: 11 of a briefed 17 pages were false positives from one
  overloaded word (`absorber` is also a MATERIAL; `clamp` is also a GMRES restart), while an
  UNLISTED page carried a present-tense-false bound spelled only in LaTeX. Also grep the NUMBER
  (a stale `[1/5,4/5]` found the page no symbol grep reached). → L-054
- **⭐⭐ When a DIRECTORY moves, census the DIRECTORY and every artefact SIBLING — a grep keyed to
  one filename is blind to the rest of the family.** A `graph.db` census (4 rounds, brief-supplied)
  missed `docs/index.rst`'s `<_nexus/graph.html>` — a **404 on the docs homepage**, since a fresh
  `-E` build leaves no `_nexus/` at all while `graph/graph.html` (627 KB) sits un-linked. Grep the
  parent segment plus each extension (`_nexus/`, `graph\.(db|json|html)`); anchor the SLASH or
  `_nexus` matches every `mcp__nexus__*` tool name (559 KB of noise). ⚠ **A raw relative hyperlink
  is checked at NO severity** — unlike `:doc:`/`:ref:`. And when a static link genuinely cannot ASK
  (no CLI from RST), ship the mirror WITH an RST comment naming the coupling: a second declaration
  that announces itself beats a silent one. → L-058
- **⭐ A now-optional flag is fixed by DELETION, not by updating its value** — updating mints the
  next stale literal. `--db` resolves via `.nexus/config.toml`, so 16 `--db <path>` lines came out
  of one reference for one precedence sentence. KEEP it where naming a file is the example's POINT
  (a scratch/override graph), which makes the flag *better* documented than when it was mandatory.
  Report the resolution ASYMMETRY too: `status` REJECTS `--project-root`, so read-only subcommands
  are cwd-anchored and their "does not exist" error reads as *"never built"*. → L-058
- **⭐⭐ A surviving MODULE does not license a repoint; a surviving CLAIM does — adjudicate a dead
  `:mod:` by SENTENCE TENSE *and* by whether the NAMED OBJECT survived.** `[M]` 4 dead
  `orpheus.sn.spatial.*` in one historical entry: 3 of 4 modules survived a pure `git mv`, yet the
  verdicts split 2 literal / 2 repoint. The trap is two sites in the SAME file with the SAME rename
  and opposite answers — "**Documented in** X" (claim still true there ⟹ repoint) vs "what Phase B
  added … Protocol (X) with three strategies" (Protocol + 2 of 3 strategies since retired ⟹
  literal). Three free corroborations: `git log --diff-filter=D` on the old path; the LIVE tree's
  own prose (it spells the retired names as ``literals``); and the same page already spelling the
  deleted path as a literal 130 lines below. A list where 2 of 3 `:mod:` are live argues against
  literalising the third. → L-061
- **⭐⭐ A CROSS-REFERENCE INSIDE HISTORY IS A CATEGORY ERROR — a body that records the code *as it
  then was* cannot carry a role, which claims the symbol exists NOW at THAT path.** The rule that
  adjudicates every site, and belongs IN the page as a head-of-block `.. note::`: *a name is a
  ``literal`` whenever the sentence around it describes the code as it then was; a role is used only
  where the sentence is a present-tense claim about something that exists now.* So a SURVIVING CLASS
  does not license keeping the role — the surviving CLAIM does: I literalised a live class because
  its sentence stated a `ψ_{1/2}=0` default the SAME entry's later section records as replaced.
  Nothing is lost — the live pointers move to ONE forward-orientation paragraph where their tense is
  present. Corroborate before editing by counting BOTH spellings of each name in the same file
  (`[M]` 5 literals/3 roles, 4/2, 2/2 — the page had already settled on literals and was
  self-inconsistent). → L-062
- **⭐ Before calling a dotted target dead, decide WHICH SEGMENT died** — package / module / class /
  attribute have different repairs and only one is "de-role". `[M]` two brief classifications
  refuted by the same probe error: `…continuous.sood_registry` is a live **package** (a `.py`-only
  check misses it) and `SNMesh.pole_angular_closure` is a live **instance** attribute
  (`self.x = …` in `__init__`; L-053c again). And a dead CLASS under a live PACKAGE may have a live
  HOMONYM elsewhere — repointing is then a false attribution. → L-062
- **⭐⭐ A raw FILE PATH in a literal is the same category error one register down, and no
  instrument sees it** (roles are gated, paths are not). `[M]` in one ERR entry **14 of 14**
  `tests/*.py` paths no longer existed; catalogue-wide **40 of 100**. Fix the class, not the
  paths: state once, at the head, that *which* tests catch ERR-NNN is never prose — it is the
  `@pytest.mark.catches("ERR-NNN")` marker set, read with `nexus errors` /
  `context('vv:error:ERR-NNN')`. → L-062
- **A retirement propagates to BOUNDS and to NEGATIVE claims, which no symbol grep sees.** A
  numeric bound is a claim about the retired object (re-measure it from the live producer). And
  grep the retired name inside `independent of|unaffected by|does not depend on` — a negative
  claim about X is exactly what retiring X can falsify (`[M]` "the floor is independent of the
  τ-clamp" → removing it moved the floor 1.8–3.4×). → L-054
- **A retired MODULE whose FUNCTION survives by NAME is a semantic trap, not a repoint** — if the
  survivor now DELEGATES to production it is no longer an *independent reference*, so a mechanical
  repoint yields working links to false COVERAGE claims. Fix: past-tense literals in the history +
  ONE anchored `.. note::` stating the delegation, WHY, and what each arm compares against now.
  → L-054
- **When a carve proves a published diagnostic is Mode-12 BLIND on the shipped rule, the page owes
  a `.. warning::` with the garbage-passes table AND the instrument that DOES discriminate.**
  → L-054

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
- **⭐⭐ Fixing HALF a claim in one file is worse than fixing none — after repairing a section,
  grep the WHOLE FILE for the retired predicate's spellings and adjudicate every hit by tense.**
  A brief scopes you to a section; the same falsehood routinely survives three screens above it,
  and a self-contradicting page is citable for EITHER sentence (vv #21's aggravator). Seen: the
  repaired selection §, with the upstream § still opening "quadrature selection therefore reduces
  to a containment check". Re-scope the survivor to what the equation ACTUALLY is (an order
  relation), add a `.. warning::` saying it is NOT the gate, keep the label + `vv-status`
  untouched (it was `implements`-cited from production). Grep list = every spelling of the retired
  symbols plus the stage COUNT. → L-056
- **⭐ A tombstone may only assert what YOUR page controls.** "…and the module docstring said the
  same until <date>" is false the moment the other file is fixed — or reverted. Write a twin's
  history in the past tense of the CLAIM ("the promise was minted twice"), never of the file's
  state. → L-056
- **⭐ A stale FORMULA can be right on a biased subset of the grid, so a spot-check CONFIRMS it.**
  `max(3, N−1)` matched the measured level-symmetric degree at S2/S12/S16/S18 and missed at
  S4/S6/S8/S10/S14 — and S12 is the order the stale frontier itself made salient. Cite the
  SWEEPING gate, and fix a drifted number with a `:ref:` to the SSOT, not a fresher copy: the
  producer discovers it by construction, so any second copy is the thing that drifts. → L-056
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

- **⛔ A mid-task scope REVOCATION on a file you already edited: revert by RE-EDITING, prove it
  with `git diff --quiet -- <path>`, and publish the backed-out patch in your RETURN.** Never
  `git checkout`. Afterwards the file may show Modified again — that is the concurrent editor, so
  discriminate by grepping YOUR OWN distinctive strings, not by the porcelain flag. The addendum
  named 2 of the 4 sites I had found, so the return is the only place the other 2 survive.
  → L-056
- **An enumerated list starting at `0.` is legal** (docutils sets `start="0"`, INFO-level only,
  suppressed at Sphinx's default verbosity) — use it when the code numbers its stages 0..N, so a
  runtime message names the paragraph that explains it. Probe with `publish_doctree` first.
  ⚠ And never wrap quoted stale text in `*…*` when it contains `:math:`/`:eq:` roles — docutils
  does not nest inline markup; use the page's own ``⛔ X read Y until <date>`` idiom. → L-056
- **Never hand-retype a large block.** Read → slice → `"".join` → write, with guard-asserts on the
  LIVE file's boundary strings and lengths, and ALL structural asserts run on the in-memory result
  BEFORE any write (a failed assert then leaves the tree untouched — no `git checkout` recovery
  needed). A machine splice cannot mis-transcribe. ⚠ **A red guard may be the GUARD's error —
  diagnose WHOSE failure it is before touching content** (`vv` #4's VERIFY sharpening, turned on
  your own instrument): a `len(out) < len(src)` guard fired on a splice that legitimately GREW,
  while an earlier assert in the same script had already caught a real miscount — that positive
  control is what made the false red cheap. → L-012, L-022, L-023, L-026, L-058
- **⭐⭐ PROBE docutils, never reason about it — and stand up a stub-directive harness so you
  iterate in 1 s, not 4 min.** I predicted a warning came from *emphasis ⊃ inline literal*: WRONG
  (that is silent and renders raw backticks); *emphasis ⊃ **strong*** is what warns. One
  `publish_doctree` call with 6 one-liners settled three entries at once (`key=``x``` warns,
  `key=\ ``x``` is clean; `γ_-` errors, `γ\_-` is clean; a nested list needs a blank line, a `+`
  mid-paragraph does not). Register `error-entry`/roles as pass-throughs and re-check a 5790-line
  file in under a second. **Markdown discriminator for an indented block:** blank line before ⟹
  a real code block (⟹ `.. code-block:: text`, mandatory if the body holds a `*`); no blank line
  ⟹ a lazy paragraph continuation (⟹ blank lines around it → block quote). → L-061
- **⭐⭐ When the edit is confined to ONE section of a big file, the strongest guard is BOUNDARY
  BYTE-IDENTITY:** `src[:i] == out[:k]` and `src[j:] == out[m:]`, with `i/j` and `k/m` the section's
  own delimiters. Two lines prove the other 78 entries of a 5800-line catalogue are untouched — far
  stronger than any per-edit count, and it makes a 22-site sweep defensible. Pair it with an exact
  `len(out) == len(src) + Σ n·(len(new)−len(old))` arithmetic delta per replacement table. → L-062
- **⭐ Run a SELF-CONSISTENCY pass on prose YOU authored before the first build, not after.** New
  prose that DECLARES a rule must obey it: I built four times because I kept finding my own note
  violating its own convention (one live class left a literal), over-claiming a successor, and
  leaving one paragraph ragged. Ask of every name in your new text: *which branch of my stated rule
  does this take?* — then build once. → L-054, L-062
- **⭐ Guard a bulk DELIMITER edit with `src.replace('`','') == new.replace('`','')`** — proves
  only backticks moved — plus an exact char-count delta and unchanged line count, all asserted
  BEFORE the write. That is what makes a 415-site blanket edit defensible rather than reckless.
  → L-061
- **⭐ A uniqueness guard over LABELS or TITLES compares EXACT LINES, never substrings.** Two of
  mine fired (before any write, so free): `count(":label: operator-apply")==1` fails because it is
  a substring of `:label: operator-apply-transpose`, and `count("Development history")==1` fails
  because an `.. admonition:: Development history — …` sits 1000 lines above the section.
  Eq-label families are BUILT by suffixing (`X`, `X-transpose`, `X-section`), so the prefix
  collision is the normal case here. Use `sum(1 for l in lines if l.strip() == …)`. → L-060
- **⭐ A directive whose BODY renders needs a placement rule, or 50 of them land mid-sentence.**
  `.. implements::`/`.. verifies::` with a body emit a plain `<div class="docutils container">` —
  visible prose, no marker. Rule: **after the `.. math::` block, unless the next paragraph is a
  grammatical continuation of the equation's sentence** (`where …`, `so …`, `with …`, `and
  identically for …`) — then after that paragraph. Encode it as a per-label `skip` flag and
  PREVIEW before writing; open every body `**Implemented by** …` so it reads as an annotation.
  → L-060
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

- **⭐⭐ A DIALECTICAL SEED PAGE (a design dialogue CONVERGED, first slice shipped) is NOT the
  9-step close-out arc** — that arc is for a CLOSED "cannot work". Order: Key Facts *carrying
  the doctrine's one-line discriminator tests verbatim* → the taxonomy → the theorem → **the
  doctrine dialectically** (question → v1 REFUTED → v2 REFUTED → standing → retrodictions) →
  fences per phase → dev history. ⭐ Give each refutation its own `.. admonition:: ⛔` **titled
  with the REFUTING QUESTION, not the verdict** — both refuted versions are *almost* right, so
  the question is the transferable content and a reader who meets only the final statement
  re-derives v1 within a week. ⭐ And say explicitly what the doctrine does to the tension it
  settled ("it does not pick a winner — both prior rules are right about different clauses");
  otherwise the next reader hunts for the loser. → L-064
- **⭐⭐ A RETRODICTION / confirmation table is `plan-authoring` §2's aspirational-row trap moved
  into the CORPUS — and it costs more here.** A table headed by a property of the tree reads
  ENTIRELY as a survey of what IS, so one unbuilt row is indistinguishable from the observations
  and, once found, discredits every row. ⟹ **STATUS column, in the row** (`[M] ships` vs
  `⛔ NOT built — a prediction`), never prose above or below; and head it *"rows the doctrine was
  NOT built from"* (the real epistemic claim) rather than *"layouts the tree ships"*. Caught in
  my own draft by counting my own universal. → L-064
- **⭐ Citing an SSOT: name the REGISTER your page owns, not just the fact** — the same fact in a
  different register is not a twin. `frame.rst` owns the counting measure in the MEASURE register
  (`w_g = 1` vs `Δu_g`, Hébert-derived, rate-preservation-gated); my page owned it in the METRIC
  register (`G_E = I` ⟹ `V ≅ V*` ⟹ adjoint = plain transpose ⟹ construction refuses weights).
  Derive it a third way as UNLABELLED math, cite the SSOT's label, open with
  `.. important:: … Edited there, consumed here.` **Net new labels on a 1158-line page: ONE.**
  → L-064
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
- **⭐⭐ A NEW THEOREM: home it where its LOCAL half is already derived, then audit the
  UNIVERSAL it amends tree-wide.** The obvious homes (the BC page, the solver page) were both
  wrong — the page that already labelled the local `−1` face mode owned half the mechanism, so
  the global result is a downstream H1 there and anywhere else mints a twin. Then: the headline
  ("a splitting shares a solution SET, not a POINT, when `A` is singular") falsified an
  unqualified claim asserted **9× across 7 files**; a windowed regex finds them and the
  adjudication is NOT uniform — scope-to-**bulk** where the measurand is bulk, `.. note::`
  tombstone where the sentence stays, one clause where a chapter-scoped truth still says "any",
  and LEAVE where the sentence quantifies over something else (two source-DELIVERY routes of one
  iteration is not a splitting claim). ⭐ Auditing the prose, MEASURE every gate fixture the
  corpus names against the new pathological predicate — one of them WAS singular
  (`[M]` `dim ker A = 36`), and saying *why the gate survives anyway* (its measurands are
  mirror-even) beats deleting or keeping the claim. ⛔ And the changelog entry can be BLOCKED
  by the page's own contract ("merged to main, with its hash, or not at all") — report the
  ready-to-paste row, never fake the hash. → L-057
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
- **⭐ An ONTOLOGY-OVERTURN changelog goes on the page whose THESIS moved.** `history.rst`
  contracts "a new entry lands with its merge hash or not at all", so an unmerged carve is
  BLOCKED there; `operator_algebra.rst`'s history has the `*(in development)* <branch>`
  escape hatch. ⟹ give the rewritten page its OWN Development history following that
  convention verbatim, a short row on the sibling whose axis genuinely moved, and on the
  blocked page tombstone only the falsified HALF of its row (one row, two halves, opposite
  fates). → L-063
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
- **⭐ A retirement's dead-ref count is under-reported by the commit that made it** — a brief
  and a commit body both said "23 dead refs" (one retired package); the second retired module
  path added **9** more. Grep every retired path, not the one the message names. ⚠ And a
  deleted package can still IMPORT: an untracked `__pycache__` leaves a PEP-420 namespace
  package (`__file__ is None`, 0 members) that a naive `import_module` probe calls LIVE —
  probe a SUBMODULE. `[M]` the xref gate saw 1 of 32 (L-062's unlanded `head_role` bug); my
  own import probe over 727 roles across 8 edited pages is the real gate. → L-063
- **⭐ Measure the xref-gate baseline from `git archive HEAD` into a temp tree** — the cheap way to
  get a TRUE before/after on a dirty working tree (`git archive HEAD orpheus tests docs tools | tar
  -x -C <tmp>`, then run that tree's own copy of the gate on it). `[M]` 81 dead / 124 sites both
  sides while adding 80 xref roles. Its file-count will differ (untracked files); the DEAD number is
  the gate. → L-059
- **RE-MEASURE the `-E` baseline every session; never assume a recorded number** (it has drifted
  9 → 1 → 0; measured 0 again 2026-08-11). Diff the WARNING/ERROR/CRITICAL SET pre/post, not just
  the count. A full `-E` rebuild can exceed the 120 s foreground cap — background it, or the
  poll-loop is SIGTERM'd at the final line. → L-029, L-041, L-027
- **⭐ SEQUENCE the session so you build TWICE, not four times:** baseline `-E -W` → *all* edits →
  *all* residual greps → xref gate → AST doc-only proof → ONE verification build. Two of four
  builds were wasted by launching before the last edit landed (a residual grep always finds one
  more site). ⚠ **Re-broken on a NEW page: FIVE builds**, every extra one bought by an edit made
  after launching — the self-consistency pass (universals, symbol collisions, aspirational rows)
  must run to EXHAUSTION *before* the first verification build, never interleaved with it.
  → L-054, L-064
- **⭐ ONE re-runnable python self-check beats the build for structure, at ~2 s:** short-underline
  detection + ladder-order (first-appearance) + per-table column consistency across EVERY row +
  `:widths:` sum + label/anchor uniqueness (EXACT-line compare) + role import-resolution +
  `:eq:`/`:ref:`/`:doc:` resolution against the whole `docs/` corpus. It caught every structural
  defect on a 1158-line new page before any build ran. → L-064
- **An agreement `[M]` number handed to you is a LADDER unless proven flat** — a briefed "verified
  to 1.67e-16" was a small-fixture reading that degrades to `2.3e-14` two refinements later, and
  the shipped gate already knew (`atol=1e-13`). Measure the ladder, publish the ladder, name the
  gate's tolerance, say a finer row must widen it (`vv-principles` #16). Likewise read a
  two-number "spread 0.30→1.53" as ONE row, not a sequence, until re-run. → L-054
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

# QA Lessons — hot digest

Read every dispatch. **Behavioral rules only**: one imperative, the check that
makes it decidable, and its `→ L-0NN` archive pointer.

- **War stories, evidence, `file:line`, measured tables** live in
  `lessons_archive.md` (`## L-0NN`, L-001..L-106). Open only the `L-0NN` a rule
  points at; never read it whole.
- **Doctrine is NOT restated here.** The preloaded skills own it: `vv-principles`
  (#1–#38, Modes 7–12, bit-identity, 1-group degeneracy, the `catches`
  directive), `numerical-bug-signatures` (Sig 1–11, H1–H5),
  `instrument-doctrine` (X1–X4 and their procedures), `retirement-audit` (items
  1–25), `coding-elegance`; the qa definition (`qa/AGENT.md`) owns make-it-RED,
  the behaviour-neutral retype, coverage-by-what-ran, the promised-vs-shipped
  read, the prose audit (§4) and the verdict stamps (§5). A `[skill: …]` tag
  means the principle is there and what follows is only the ORPHEUS mechanic.
  §I says where the retired rules went.
- New lesson ⟹ append `L-0NN` to the archive, land a 2–5 line rule here. Sharpen
  in place; never let this file grow narrative.

---

## A. Making a gate RED — mutation mechanics

`vv-principles` #17 (a)–(j) owns the arm taxonomy, #18 the algebraic-class rule,
Mode 8 the fires-but-cannot-fail classes; "a zero is a claim about the
instrument first" is the definition's §3. Below: the ORPHEUS mechanics, and the
arms that taxonomy does not name.

**A3. Revert a mutation by RE-EDITING.** [rule: `process-discipline`, "Mutation-testing
an uncommitted file"] Untracked files make `git diff` empty, so the revert proof
is gate-green-again + zero mutation markers; a `-p <module>` plugin needs
`PYTHONPATH`. When no file edit is licensed at all, hand the function a COPIED
input (a mutated entry dict): control leg silent, mutated leg fires, nothing to
restore. → L-039, L-043, L-052, L-086

**A4. Your OWN mutation needs a bite check.** [skill: Mode-8 METHOD WARNING] A
capability REFUSAL is TWO-part — adding `apply_transpose` does not lift
`is_adjointable`/`is_invertible` (predicates defaulting False). A **0-call
counter is a FINDING** (path unreachable), not an inert mutation. → L-061, L-062

**A5. When `simplify` is pathologically slow on the MUTATED expression, call the
`derive_*` builder with concrete Rationals and read the residual** instead of
waiting on pytest. Seconds, and decisive. → L-029

**A6. Cripple a GENERATOR, not a value, for the sharpest coverage verdict** —
replacing `O_h`'s 48 ops with its 8 diagonal sign-flips (= `D_2h`) left a
182-test suite green. → L-062

**A7. ONE mutation direction is almost never enough — enumerate the leaks.** (a)
*capability default-OFF*: the factory AUTO-SELECTS the wider shape or appends a
PHANTOM length-1 axis (control: `not hasattr(space,"factors")`); (b) *`xfail`→live
flip*: red against the re-introduced bug AND the EMULATED PRE-change behaviour —
only the latter rules out a gate already green at HEAD; (c) *polymorphic hook*:
override returns the base type, AND override DROPPED (base `replace()` keeps
state, so only the empty-state tooth reds, and only if the test ADVANCES first). → L-032, L-038, L-041

**A10. Run the mutation over the WHOLE module tree in BOTH arms; the symmetric
difference turns "does an external pin exist?" from an argument into a LIST.**
`[M]` old-τ vs HEAD over `tests/gates/sn`: 7 red only at HEAD, **32 red only under
old-τ** — the 32 named the analytic pins I had just concluded did not exist.
Bite-check first: the target gates must FLIP, with a non-zero call count. → L-069

**A13-r. Split READS from WRITES before reading a field's test-hit count as
coverage** (5 of 5 hits on `StreamingTerms.mu_start` were constructor kwargs);
and a zero-reader field naming a REAL contract is usually RESPELLED, not dead —
ask how production answers that question today before wiring it. [skill: #17(e)] → L-075

**A14-r. ORPHEUS traffic-census mechanics.** [skill: #29 owns the claim] Patch
the singledispatch **registry**, or wrap `cls.__dict__["apply"]` through the
descriptor protocol; check `apply is _apply_impl` by identity; attribute per SITE
with `__init__` + `extract_stack`; keep instances alive so `id()` cannot recycle.
Controls are a LADDER — instrument → installation marker → **per-ARM activation**
(without it, 8 of 23 zeros are unreadable) → headline bit-identity. Bound the
workload with a STATIC reference census so "measures its workload only" names the
real residual. → L-076, L-073

**A15. A refusal DUPLICATED one frame down leaves the outer guard witnessless:
run BOTH arms — swallow the inner raise, then let it propagate.** `[M]` deleting
a normaliser check reddened **0 of 670** because `Quotient.induced_action`
refuses the same motions; without the swallow it reddened **16**, so the outer
guard's real job is converting a raise into a `False`, which its gate does not
assert. A 0-red arm is *guard with a twin* or *provable optimisation* — decide by
algebra, not red count. Per-ARM pays: one body → 1 red, its five arms → 7/27/1/1/1. → L-077

**A18. Scope a mutation to ONE PHASE by rewriting the RECORD that phase reads,
never by flagging a shared verb.** A global mutation on a verb the driver also
calls moves the CONVERGED answer instead of the reconstruction and reports false
coverage; wrapping the producer and rewriting the record after it returns leaves
the driver honest (`[M]` 8 of 14 rows red, exactly the 4 arms with a live `B`).
⭐ The wrapper is a free CENSUS: 161 inner solves / **0** scheduled splittings
named, in one number, an arm the new gate module never reaches. → L-080

**A22. A generator asserts source→target; ask who asserts target→source.** `[M]`
deleting a manifest entry left its generated always-on rule (≈708 tokens) on disk
with `--check` reporting `0 problems, 0 drifted`. For a blind spot with an EMPTY
present population: call it LATENT, give the denominator, and prove the class
with a SYNTHETIC positive control rather than a tree witness. → L-082

**A23. Retirement notes: the residue `retirement-audit` item 19 does not carry.**
[skill: item 19 owns "a note is a carrier claim; count the source's rules; read
its closing paragraph"] A ⛔ PROTECTION justified by "loaded at every session
start" survives into a COLD archive when a hot file is split: grep a split
archive for self-references to the hot file's LOADING. → L-083

**A24. For a PROSE census the positive control is a phrase PASTED from screen and
`\s+`-normalised** — a remembered phrase carries neither the source's line
wrapping nor its `**` markup, and both failed silently in one session (one nearly
published "no carrier exists"). Before an `L<n>`/`#<n>` count means anything,
discriminate the namespaces: a V&V level is not a lesson, "<agent> lessons L9"
names THAT agent's memory, "`<skill>` L11" may name nothing at all. → L-083

**A25. A rebinding loop compares against the ORIGINAL captured before the loop,
never the defining module's live attribute** — rebinding the definer first
turns every later `is` test false, and the activation counter still reads
non-zero from the definer's own internal callers. Attribute activations PER
CALLER MODULE, or the count certifies the wrong site. → L-095

**A26. A source-transform mutant exec'd into a COPY of the module namespace reads
module globals from the copy, so a test's `monkeypatch.setattr(module, name, spy)`
never reaches it** — a recording-spy row reads RED ("never called") for a harness
reason. Exec into `module.__dict__` itself, or read every spy row's message before
counting it (`[M]` 2026-10-12: G1b red under `ratio2` and `no_middle`, the only
red of `no_middle`). Not `vv-principles` #17's parent/child process mechanism. → L-102

---

## B. Where a gate is structurally blind (ORPHEUS shapes)

**B-DISPATCH. A claim about what a DISPATCH RECEIVES is measured by
DISPATCHING** — frontmatter, a plan ruling and a `keep − omit` token probe all
read TRUE on a false one (the measured facts are now the `workflows` rule,
invariant 4). A DIFFERENCE instrument cannot name what sits in BOTH arms, and a
probe with 0 tools carries no role block, so role-block growth cannot move it:
price each part with a probe that CARRIES it. ⭐ A `skills:` preload is a claim
about a TURN, not a dispatch — a skill generated mid-session arrived only at
turn 2 — so a first-turn absence is **unadjudicated**, never "broken". → L-085, L-089

**B1. Instrument a CALL COUNT first, to find which twin the gate actually runs,
and mutate the SHARED source rather than the dead-for-this-path method.**
[skill: Mode 11; X1 activation count] Sweep and matvec share only precomputed
coefficients, so three apply-path mutations gave call-count 0 and identical error
ladders. ORPHEUS traffic facts: SI sweeps never touch `loss_action`, which runs
only under `inner_solver="krylov"` (1600 / 0 kernel calls on an MMS solve); and a
"fires under quadrature Q" claim is a 3-line probe
(`count_nonzero(|mu_x|<1e-15)` → zero at every LS order). → L-016, L-018, L-021, L-033, L-036, L-059

**B5. A fixture SYMMETRIC in the axis under test cannot see that axis.** Three
ORPHEUS shapes: a SQUARE `nx==ny` mesh hides axis ORDERING (and the algebra-law
suite is swap-invariant anyway, so the catcher is a broadcast oracle at `nx≠ny`);
a UNIFORM fixture makes a per-cell and a global-mean check indistinguishable, so
one fixture must VARY along the non-reduced axis; two SAME-AXIS faces make
`|Ω·n|` bit-identical and annihilate the packing gate's only knob-reader (`[M]`
0/10 red, `changed=False` every call; a y-face moves it 0.963). → L-030, L-040, L-065

**B9. A "the matrix says the operator is healthy" argument must cite a
certificate for the EXACT gated BC** — a sibling-BC certificate plus "same
mechanism" is inference. Build the fully-coupled matrix: `ρ_prod > ρ_matrix` ⟹
splitting/wall lag (honest); `ρ_matrix ≈ 1` ⟹ real consistency failure. → L-059

**B10. The headline category gate is usually the WEAK one.** `runtime_checkable`
checks member PRESENCE only — `isinstance` flips True on a monkeypatched attr and
stays False under a realistic PARTIAL leak; the direct `not hasattr(...)`
negatives are the defense. Credit them, not the headline. → L-039, L-042, L-047

**B11. The MMS refinement ladder is BLIND to the diffusion limit — probe
`σ_t·h ≫ 1` on a COARSE mesh**, where users actually run (refinement drives
`σ_t·h` thin, where flat-source is fine). Probe against DD with an ε-scaled
diffusive material; a reflective `c≈1` probe is a TRAP, both schemes reading
~82 % wrong from non-convergence. → L-017

**B12. A CONVERGENCE flag at an eigenvalue entry is the OUTER fact only.**
`solve_sn`/`solve_sn_adjoint` warn on `max_outer`/`keff_tol`; a within-group
solve hitting `max_inner` never reaches the warning. Ask WHICH LOOP the flag
belongs to and wrap the within-group exit check for the other (`[R]` that check
is `_check_convergence_claim` since `6e04f1f0`, formerly
`_certify_within_group_exit`). Companion holes: a suppressed warning, an
`xfail`-absorbed one, `-m "not slow"` deselection (#36). → L-067, L-053

**B13-r. Gate a published COMMAND STRING with
`_pytest.config.parse_warning_filter(s, escape=False)`.** [skill: Mode 8(8)]
`[M]` `-W error::ConvergenceWarning` at 4 doc sites does NOT parse (an undotted
category resolves against `builtins`; pytest exits ERROR, 0 collected) while the
file's own "it is escalatable" test passes, installing the filter
programmatically. → L-067

**B14. Coverage by what ran: the ORPHEUS ladder and its blind rung.** [definition
§3, "Coverage by what ran"; `nexus-verification` owns `run=` and the capture
config] `[M]` 0 of 21 claiming tests reach `Quadrature.ordinate_permutation`
statically, **7** do at runtime, because every production call site is
annotation-mediated: `nexus callers` → 0 and `dead-functions` FLAGS a live method.
Join recipe: `[json] show_contexts = True`, joined to
`sphinxcontrib.nexus.runtime.build_node_index` spans. ⚠ It measures
CO-EXECUTION, not co-constraint — a candidate list to mutation-verify. The ladder
is CLAIMED 21 → EXERCISED 7 → ASSERTED ≤2 → MUTATION-VERIFIED 0, and **no edge
quality separates rungs 2 and 3.** → L-070

**B16. When a MANIFEST becomes DISCOVERY, audit by asking "for which dropped
check did the input stay REPRESENTABLE?" — not "which checks are gone".** Two
survived that filter and were deleted anyway: an empty-front-matter skill (fatal
for the skill, `--check` green) and "target must lie under the harness dir"
(demoted to a per-IMPLEMENTATION test, so a second `Harness` escapes — the
Protocol seam is what makes that insufficient). Riders: a function parsing the
same language twice inherits the first parse's error contract at the second site
(an unguarded `yaml.safe_load` crashed `sphinx-build`); an AST import census
reading only `ImportFrom.module` misses `from . import X`; a hand-typed tuple as
a gate's population excludes every future member; and a rename can mint a summary
line whose arithmetic is exact and whose LABEL names a population the tool cannot
see (`always-on ≈20341` omitted 7479 hand-maintained tokens, 27 %). → L-087

**B17. A quadrature whose resolution is a caller's COUNT (grading layers, points
per piece) is attacked at the optical scale of the extreme PIECE — the thinnest
panel (Sigma h -> 0) and the thickest region (Sigma r >> 100) — under a law whose
functional reads that rule (vacuum, alpha < 1), never only under closed-body
conservation, which each line satisfies on its own.** check: one mpmath closed
form (P_esc, P_ss, white chain) at tau 1e-6 and 1e2 beside the gates' range.
→ L-098

**B18. A normalisation (a gauge, a ScaleGauge section) is attacked with a state its
functional annihilates ANALYTICALLY — a higher mode of a closed homogeneous body, a
null-space mode of a rank-deficient production — never only the fundamental; an
`== 0.0` refusal is blind to the rounding zero.** check: |n(ψ)| against
ε·n(|ψ|) on mode 2 of a closed body. → L-099

**B19. A grading law added at the site where a defect was measured (the rim) is
attacked at every site its MECHANISM reaches (each interior tangency, each point
on a radius, each region below a thinner one), and a non-converging plateau
labelled "[R] rounding" is split before it is believed.** check: one probe per
site of the mechanism; and the rule's sum of the exact integrand at its exact
nodes against the same sum at the nodes the object REALISES from its stored
coordinate. → L-100

**B20. A row moved behind a memo is blind to in-process mutation in the warm run, so a
bypass battery proves its teeth only under bypass.** [skill: #17's closing paragraph owns
"mutate a copied tree, or `bypass()`"] check: for each arm whose reds fell ONLY on re-posed
rows, re-run it WITHOUT bypass on a store target over the default tier; a red names the
surviving in-process pin, a green is a lost pin; and a `catches` marker's row is re-read the
same way. → L-103

**B21. When a carve makes the algebra NORMALISE (flatten, collect, collapse), census every fact
STAMPED on an operator instance (`op.block_role = …`, `A.system_role = …`) and read it again after
`-op`, `2·op` and `op + other`, at every amplitude.** A normaliser mints fresh objects and the stamp
does not travel; the reading is amplitude-dependent when `1·X is X` keeps it and `α·X` does not. → L-104

**B22. A tool gated only on a synthetic fixture TREE is blind to every input the REAL repository
supplies that the fixture lives outside of (the recorder, the tool's own code): they are not
first-party there.** check: one micro-capture of 5-10 real tests (about 30 s) read through the
tool before crediting any fixture gate; compare against the spec's named expectations. → L-106

---

## C. Reference contamination & structural independence

**C1. The circularity test is: does the bug live on BOTH sides?** [skill: #6, #7,
X4] Two ORPHEUS tells that a corroboration is only PROCEDURAL: the two sides ride
one **antiderivative identity** (a discrete recursion summing `f` "confirms" a
claim about `F=∫f` only because `F'=f` — true whatever the claim is); and the
corroborating gate's own docstring cites, as ITS reference, the claim being
corroborated. → L-029, L-068

**C2. For a WEIGHTED value-pin, independence is not enough:** the hand-reference
must carry EVERY weight factor AND the fixture must make a factor-BLIND formula
give a different answer. Hand-compute the blind number, then mutate production
blind and confirm only that gate reds. → L-046

**C3-r. When a carve leaves the composite byte-identical, the anchor for a
re-baselined `.npy` is `composite − collision`** — never "whatever the changed
leaf emits". [skill: § Bit-identity, Sig-10: never old-vs-new ULP] → L-034, L-049

**C4. Two independent implementations IS independence;
producer-vs-its-own-projector is not.** Compare production's emitted slot against
a TEST-SIDE `leggauss` reference, then separately pin the two projectors'
agreement. Same family: a brief's "0 ULP" after an API migration is a CLAIM until
you recompute the OLD contraction on a structurally-independent table; and a
two-paths oracle's analytical anchor is often TRANSITIVE and in another file —
confirm that file is green before crediting analytical grounding. → L-028, L-044, L-051, L-052

**C6. For an adjoint, the independent reference is a DENSE matrix built by LOOPS,
transposed directly, composed with metrics by hand.** Re-derive the inner-product
identity first, then prove `(A⁻¹)ᵀ=(Aᵀ)⁻¹`. An ASYMMETRIC transpose pair is
CORRECT when each transpose mirrors its OWN forward. Where the metric commutes
with `Aᵀ` (a spatially-diagonal operator, Mode 12's `[G,Aᵀ]=0`) no `.H` gate can
witness a metric flip; the honest witness there is a construction refusal. → L-052, L-060, L-073 §2

**C8-r. Prove a helper's independence MECHANICALLY:**
`dis.Bytecode(f).codeobj.co_names` — the forbidden names usually appear only in
the DOCSTRING. The INPUT-RESOLUTION axis closes silently: `[M]` an axis-letter
x↔y swap left the "genuinely independent routes" file 15/15 GREEN and reddened 78
siblings. [skill: #22 owns the two axes] → L-064

---

## D. Re-baseline & bit-identity integrity

**D2. Run the MASKING-CHECK on any loosened gate or regenerated baseline.**
Loosened → re-run the untouched arms and confirm they STILL hard-fail ≫ the
bound. Regenerated → OLD-snapshot-vs-NEW-code must hard-fail (load-bearing) AND
NEW-snapshot-vs-OLD-code must hard-fail (live gate). → L-022, L-028

**D3. Characterize drift from the BINARY** (`git show <c>~1:x.npy` vs
`git show <c>:x.npy`, then ULP-diff): live-code vs regenerated-snapshot is
necessarily 0 ULP and characterizes nothing. → L-022

**D5. A HARD nULP floor and a STRICT bit-identity floor are different gates —
verify WHICH invocation ran.** Strict = `-W error::DriftWarning` layered on top;
`tests/gates/sn/regression/conftest.py` downgrades it for its own directory and does
NOT leak to siblings (measured — assume neither way). Prove a strict floor live
by perturbing the baseline 1 ULP (`np.nextafter`). → L-014, L-015

**D6. Settle a byte-identity dispute with the IEEE micro-fact and
`git status --short '**/*.npy'` — NEVER a docstring.** `0.5*(a+b) ==
0.5*a+0.5*b` bit-for-bit for all doubles (so a `w=½` affine closure IS
byte-identical to DD, contra its own docstring); `2*X/D ≠ 2*X*(1/D)` at 1 ULP
(the real re-baseline trigger); an einsum spectator lift `fc->gc` ⇒ `fc...->gc...`
is `array_equal` at rank-2. → L-020, L-028, L-032

**D7. When a carve preserves the COMPOSITE and not the leaf, prove byte-identity
on the composite DIRECTLY** (both emitted against a read-only baseline worktree);
a brief's "≤16 ULP" can understate leaf drift ~7×. Say which object is pinned. → L-049

**D9. For an ADDITIVE-only change, grep for ANY importer of the new module
(excluding its own tests): empty ⟹ it cannot perturb a pre-existing outcome** —
stronger evidence than re-running the baseline reds. → L-042

**D10. Prove a `singledispatch` alias rename via
`Cls.__dict__['apply'] is Cls.__dict__['_apply_impl']`**: `Cls.apply is
Cls._apply_impl` is False (a fresh descriptor per access), a red herring. → L-050

**D11. Do NOT trust "byte-identical EXCEPT one LATENT collision" — PROBE it.**
Compute BOTH branches on every call and `array_equal` them across the FULL gate
suite (the plugin reassigns the symbol in EVERY importing module; attribute by
`item.nodeid`; read under `-s`): 48 divergences at 70 % ⟹ REACHED. The two-paths
gate that found them is Mode-11 blind (shared callee), and "latent via the public
entry" can be TRUE while "latent everywhere" is FALSE. → L-035

**D13. MOVING a method to a sibling object can convert a SELF-consistency into a
cross-object coupling guarded only by EXTENT.** Ask what array the old owner's
callers relied on: if the new owner carries a COPY cross-checked by shape/length
only, a same-size-different-values pair is now ACCEPTED where it used to RAISE
(`to_local` moved operator→space; the gather reads `op.indices`, the remap
`space.ordinate_indices`), and a round-trip gate harmless while ONE array existed
becomes the gap the moment there are two. → L-065

**D14-r. Before judging a re-baseline's LEGITIMACY, `git log` the snapshot's own
directory for a commit that ALREADY made the decision** — the reds may be its
REMAINDER, and the question is then completeness. `[M]` `39b46a31`'s universal
"all 23 snapshots … the only two that changed" was scoped to ONE directory while
7 further references had moved. [skill: #25 owns the per-mechanism null-check]
⚠ Sig-10's sibling-pass discriminator is VOID when the changed code is
single-geometry (the siblings never ran it, so their green carries no
information): bisect with two revision-fingerprinted worktrees instead. → L-069

**D15. A carve that DISSOLVES branches of HEAD (a `len(radii) == 1` fast path, a
solid-only shortcut, a per-kind arm) owes its old-vs-new probe one fixture inside
each dissolved branch's domain; a probe on the general fixture certifies the
general branch only.** check: grep HEAD's body for `if len(`, `== 0.0`, `kind ==`
guards the new body lacks, and run HEAD vs new once per guard. → L-101

## E. Markers, levels, and the ORPHEUS audit surface

**E3. An audit-MISSING `catches` has FOUR outcomes — grep the production RAISE
SITE first.** (1) genuine catcher → tag it, mutation-verified; (2) the catalog's
L0 test was RETIRED and the marker did not migrate → re-tag the successor; (3)
the typed error is exported but NEVER raised → dead scaffolding, NO CATCHER, do
not invent a marker; (4) `assert_X` delegates to a WEAKER sibling → NO CATCHER. → L-054

**E1. Level conflation is SILENT in three shapes.** `@foundation` stacked with
`@verifies("<physics-eq>")` is recorded twice, so Nexus credits a physics
equation with a foundation test's parametrizations (tell: a `documented` equation
whose ONLY coverage is a foundation test); `foundation` under a module
`pytestmark=l1` emits `conflicting V&V level markers` and the intended level is
DROPPED; a self-generated regression baseline wearing `l1` is conflation — fix
the marker, not the file, and its `_load_or_skip` should HARD-FAIL, not skip. → L-007, L-058, L-061

**E6. Audit mechanics.** Orphan triage order D→B→A→C (class D = an existing test
needing only the label; ~25 % of orphans). `matrix.rst` LAGS a label rename —
re-run `python -O -m tests._harness.audit --gaps` for the live spelling.
`vv-status` rationale comments use (parens), never [brackets] (docutils reads
those as citations). → L-002, L-003, L-004

**E7. NEVER read a Nexus V&V number as a coverage claim — the surface is a SEARCH
relation wearing a PROOF relation's name.** `[M]` all 2748 `tests` edges are
`test→equation` (there is **no** test→code edge); `verified` is set iff
`len(tests)>0` with no confidence floor, so **351 of 692** "verified" equations
have no declared test, and a CP test "verifies" an SN cell-flatten identity
through the shared token `"cell"`. Ask of any status: *what predicate sets it,
and what is its weakest admissible evidence?* [skill: `nexus-verification` owns
the inferred-`implements` and `claims_*` warnings] → L-070

**E8. The graph's AST marker surface is PARTIAL.** [definition §3: a marker
census reads the resolved manifest] `[M]` `foundation` (1515 usages / 308 files)
and `regression` have no AST node attribute at all; only
`verifies`/`vv_level`/`catches`/`slow` are lifted. `catches` is an ATTRIBUTE, not
an edge, and no `ERR-NNN` node exists, so its claims cannot be joined to the
catalogue by traversal; no `.npy` snapshot is a node, so a frozen reference
cannot even be NAMED. The resolved manifest is `runtime_markers` (a
`--collect-only` capture). → L-070

**E20. A production change that makes the SUT converge FASTER re-baselines every
catcher whose defect needs the slow path to bite** (a restart cap, a budget, a
stall). [skill: § "Log every caught bug" owns decay by fixture, tolerance or
budget; this is decay by a change to the SUT] check: for each ERR the changed
solver's catchers name, re-drop the defect on the new SUT and read the red set. → L-097

---

## E′. Auditing a rewrite, a distillation or a relocation

The definition's §4 owns the method (imperative, `check:`, `tell:` per clause,
at the cited target, per modal; a text- or total-keyed diff reads clean). Also
recount PER BUCKET, never the total, and ask whether the clause's AUDIENCE moved
(E18). The instances:

**E19. Resolve the three parts SEPARATELY.** Consolidation leaves the imperative
at the citing site and moves its `check:` to the cited home, so both ends look
populated and an identity-keyed recount reads `[M]` 323 of 323, 0 lost, while one
`check:` survives at NEITHER site. Four siblings: a `**check:**` re-laid-out
INSIDE the last bullet narrows a nine-class check to one with zero text deleted;
"deleted because X holds it" and "deleted and pointed at X" are identical in a
diff, so audit the decision table per ROW; a count CREATED by a move inherits the
replaced text's arithmetic — recount it from the MOVED text; a HEADING outlives
its body, invisibly, when every citation uses the § PREFIX. → L-089

**E11+E12. A compression's defects are MODAL, not numeric** (`[M]` over 30
restored items every number was exact and all four defects were modal) — and it
also ADDS lede-universals. When one sentence is dropped it is the one saying WHY,
which is usually the CHECK: ask of every compressed clause *what experiment would
I run?* — if the text does not say, it is recognition-only. → L-084

**A20. The brief's denominator may assume a structure the SOURCE does not have**
(a brief asked for "every bold-tagged clause" of an original that had none).
State the predicate you actually counted, report that the briefed one was
unanswerable, and re-pose the answerable question — here *"did the source's
prescribed mechanical CHECK survive in executable form?"*, which localised all
four losses. → L-081

**E13. An index line and its `[body]`/`[case]` target are repaired SEPARATELY** —
a `[REMEDIED]` on the index ROW does not reach the BODY; grep the target for
every mechanism the index line names. → L-084

**E14. A rule written from a defect the SAME commit repairs ships without its
measurement** — the writer can still see the defect, the tree no longer can. The
commit adding the rule adds its evidence entry, marked `[REMEDIED @<hash>]`. → L-084

**E18. RELOCATING a clause out of the always-on tree changes its AUDIENCE, not
its home** — a text-diff audit returns clean (`[M]` 82 of 82 clauses) while a
clause addressed to a brief's READER lands where only its AUTHOR looks. check,
per clause and never per section: grep the always-on tree for its most
distinctive token, with a positive control proving the filter works, then OPEN
the named carrier and read its NEIGHBOURS — a carrier that lacked the clause has
been operating without it, and `[M]` one held the opposite instruction. A named
carrier is a claim (X3): true for four siblings, false for the fifth. → L-088

---

## F. Claim-scope — the claim is broader than the evidence

**F2-r. Exercised ≠ constrained: never collapse the middle state** (nulled /
exercised-but-unconstrained / verified). [skill: Mode 10 owns the three states
and the structural pair] ORPHEUS residue: calibrate the consumption tolerance
LIVE — a deterministic SI re-solve floors at 0.0. → L-026, L-037, L-038

**F4. A cited mutation MAGNITUDE for a metric-adjoint SOLVE must be the
full-solve value — RUN it, never the angular-collapsed 0-D proxy** (metric
conjugation of a MUTATED operator is not spectrum-preserving). A never-asserted
cited number is still a plausible-substitution error. → L-058

**F5. For a "no missed site" dedup claim the PLAN is the scope authority, not the
closeout:** a residual hit is a defect only if (a) a direct reconstruction, (b)
not transitively routed one level deeper, AND (c) in declared scope. → L-025

**F8. Check what the test HELPER tolerates before crediting an enforcement
claim** — a `squeeze_density` helper made the suite agnostic to `keepdims`, so
the bit-identity claim held only up to a squeeze. → L-042

**F9. "Matvec twin verified" is KERNEL-level; end-to-end Krylov≡SI is a separate
claim.** A loud `NotImplementedError` on the deferred half is the CORRECT
interim — but say so, and do not let a spec's wording credit the un-shipped half. → L-031, L-033

**F10. Every brief-named symbol or file is a CLAIM — confirm with `find`/grep
before editing** (two phantoms in one brief). Byte-compile no-test generator
SCRIPTS after a rewire: a broken import there is a breakage no test run surfaces. → L-051

**F11. ACCEPT a floor-CHARACTER gate; do not demand a floor-REMOVAL gate.**
[definition §5: an over-demand is a review defect] When a fix cleans a RATE but
leaves a floor, the honest claim is a falsifiable scaling pin
(`err(S32) < err(S16)/2`): a closure-BUG floor is quadrature-independent
(ratio ≈ 1 → fails). → L-009

**F12. A retired/tombstoned claim is a CONJUNCTION — enumerate its legs; one
routinely survives its siblings' death.** Legs come per SUCCESSOR (a tombstone
naming N gates may name one PER LEG: a dropped `codomain` binding reddened only
the periodic gate, the split gate's `a.codomain is b.codomain` being
`None is None`-satisfiable) and per PARTITION CLASS (SN faces are THREE-way, so
"residual zero at non-outflow" is inflow ⊔ tangential). Mutation-check per leg,
and check the fixtures can EXPRESS the survivor at all. → L-063, L-064

**F14. A plan's ⏳PENDING decisive row: use the plan's own anchors as the probe's
positive control.** [definition §5 owns "run the affordable measurement"] The
plan usually states what X must reproduce — that IS the control. Cache the
expensive reference to disk so the arm sweep is cheap; grep
`Solution.__dataclass_fields__` before assuming an attribute name (`keff`, not
`k_eff`). → L-068

**F16. An inherited blast-radius number counts a NAME, not a TYPE — re-measure it
with an in-process wrap before it sizes (or DEFERS) the work.** `[M]` "~87 reads"
was a grep for `.converged`, of which **72** belonged to an unrelated family
sharing the attribute name; wrapping the class's `__init__` + `__getattribute__`
gave **33 constructions / 0 without the field / 2 reads** — 43× over, in the
direction that defers a zero-churn fix. Route: Nexus for PRODUCERS (an attribute
read is not an edge — `degree: 1` is the graph being right), a dynamic wrap for
READERS, grep only to enumerate candidates. Pair a dynamic `0` with a static
no-other-path proof (no `**` splat / `asdict` / `replace`) or it is "not
observed", not "none". → L-066

**F20. In a multi-assembly review, read every RIVAL's self-attacks as a checklist
against your target, then push one level past the argument each answers** — a
self-attack marks the SEAM, not the depth, and the prepared defence is the tell
that the author stopped there. → L-072

**F23. A CALIBRATION of your cross-check is itself a ratio and owes the
share-a-population test; an uncalibrated corroboration beats a mis-calibrated
refutation.** `[M]` calibrating a second route on the elastic channel gave
`ΔTR/ΔP1 = 0.60`, which would have "shown" the direct route 2.6× too small — but
**327 of 421** groups run a negative corrected diagonal in the elastic leg
against **6** in the (n,2n) leg, so the factor is not transferable. ⟹ report a
corroboration at its measured accuracy class ("sign decisive; cannot adjudicate a
factor of 2"), and when a second route's convention risk exceeds the claim's,
DON'T run it — a wrong reproduction of yours impeaches a correct result. → L-078

**F24. A ladder's last STEP is not the coarse rung's ERROR.** A monotone step
at refinement ratio r and order p understates the error by 1/(1 - r^-p) (4/3 at
r = 2, p = 2); only an alternating, shrinking sequence brackets its limit
within the last step. check: before a "reference bound <= T/10" floor is
believed, convert each monotone step with the MEASURED p and re-test the floor;
ask whether the argument used to deny a sibling a bound applies here too. → L-091

**F25. When a change makes a dead arm REACHABLE, census every method the new
key unlocks, not only the one the change routes.** check: list the methods that
dispatch on the newly reachable key (`grep` the kind string), and for each ask
which gate CONSTRUCTS it through the new door; a sibling with 0 positive gates
is a live defect habitat (`[M]` `sphere_mr` fixed-source: `n_groups` 1 on 2G). → L-092

**F26. After a guard is tightened, the red loop is not the fixture census:
withdrawn and `-m "not slow"`-deselected tests never run, so their
pre-strictness declarations survive.** check: grep the old declaration's
spelling over every file that calls the guarded entry, withdrawn modules
included, and diff the hits against the relabelled set (`[M]` ERR-094: one
`white|white` CP slab left in a withdrawn, slow Peierls test). Extends
`retirement-audit` census member (8). → L-093

**F27. A guard in `__post_init__` and a retired dataclass field are both
blind to unpickling: a value persisted before the carve loads as the new type
with the old field in `__dict__`.** check: pickle the retired spelling in the
pre-carve tree (`git archive <base>` into tmp, run from there) and load it in
the carved tree; a silent load is the finding, a `__setstate__` that refuses
is the fix. Extends `retirement-audit` (no persisted-state clause) and
`coding-elegance` Pattern 4's `replace` corollary. → L-094

---

## G. Doc / prose correctness (Cardinal Rule 3 findings, not V&V)

**G6. DERIVATIVE staleness — a correction's own TODO note outlives the
correction.** A sentence of the form *"file X is stale and owes a dated fix"* is
a claim ABOUT ANOTHER FILE, and nothing in X's repair prompts anyone to retire
it; it then instructs readers to distrust a file that is now correct. ⟹ after
fixing X, **grep for pointers AT X**, not only inside it. Same shape for a
retirement TOMBSTONE naming its own successor test: `[M]` one named a test that
exists NOWHERE in the tree, renamed by the same carve. → L-071

**G1. Campaign-narration staleness: the FIX bar is "provably lies about CURRENT
code", VERIFIED before ruling** — grep the named symbol or wiring tree-wide,
`gh issue view N`. Default = KEEP; a stale line inside a RUNTIME STRING is
behavioral and a "failure here HALTs Phase X" banner is a record, both KEEP.
⭐ A CAMPAIGN-STEP NAME in a forward-looking claim is a self-expiring token: when
the step lands, grep the step's own name minus the retrospective forms
(`since|at|\(|—`) — `[M]` that found the 2 survivors among 33 hits in one
command, both claiming the step retired something it deliberately did not, one of
them 146 lines from the same file's corrected twin (#21's aggravator). Also: a
brief declaring "the known baseline reds" declares the reds of the batteries IT
ran — widen the scope and reconcile against the PARENT commit before attributing. → L-055, L-065

**G3. Reviewing a results-compilation page:** a count DE-FREEZE is CERTIFIABLE
(a live `--collect-only` proves the old literal lied); a doc RETITLE can beat the
test's own stale name and docstring — verify against the live `assert` body; a
run-book delegating detail to a config file may point at a contradicting note. → L-057

**G4. A test's own prose is the least reliable thing in the file.** A "frozen /
bit-identical to the pre-carve path" docstring stales SILENTLY when its `.npy` is
regenerated and the test file is untouched — on any regen, grep consumers for
"frozen". A cited issue number can be wrong (trust git archaeology). A prose "the
ERR-NNN class" citation is a nit; the same string inside `catches()` is a defect. → L-020, L-028, L-030, L-034

**G-NX. After a retirement, grep `:by:` in `docs/` for every deleted name: a `.. implements::` / `.. verifies::`
target that is gone is a `nexus.directive` WARNING, hence a red `-W` CI, not stale prose to defer to the archivist.**
Re-point it at the successor in the same change. → L-105a

---

## H. Mechanics, environment, probe hygiene

**Situational tripwires** — one line each; the procedure is at the archive
pointer, open it when you are in that situation.

- *Judging a `tests/` Mode-8 hypothesis:* settle it in 2 min with a synthetic
  control + a falsified COPY of a real file, both modes, then PIVOT — the premise
  is usually REFUTED (0/676 inert) and the real surface is `orpheus/` plus
  NON-COLLECTED helpers; the productive pivot is an AST census of *what the
  asserts assert* (only ~29 % pinned a VALUE). → L-006, L-010
- *Judging a value claim:* replicate the test's OWN solve helper — a naive
  `solve_sn_fixed_source(...)` defaults to vacuum, so a divergent hand-replication
  usually means YOU dropped the BC. → L-011
- *Baselining:* a READ-ONLY worktree + `PYTHONPATH`, and VERIFY it took — the
  editable `.venv` otherwise resolves to the MAIN tree (`orpheus.__file__`;
  `code-search` carries the `python -c` CWD variant); a worktree pyright count
  needs the main `.venv` symlinked into its root. → L-041, L-045, L-049
- *pyright deltas:* apples-to-apples only after line-stripping and per-file
  reconciliation (a `(file, rule, msg)` diff gives FALSE net-new when a
  type-RENDERING string shifts); and removing a `NoReturn`-poisoned return
  UNMASKS every latent error downstream, so net-new ≠ per-file delta. → L-027, L-039, L-050
- *Hunting slow/timing-out tests:* batch into runs that COMPLETE — a SIGTERM'd
  run writes no junit-xml and loses the `-rfE` reasons; mark slow PARAMS with
  `pytest.param(..., marks=...)`, not the function. → L-005
- *Running a battery beside a reaper:* a reaper keyed on a process-name pattern
  owns every process with that pattern. Run nothing that spawns the pattern
  while a battery reaps, or reap by process group. check: a red whose child
  stderr is EMPTY is a kill until shown otherwise; re-run it uncontended. → L-096

**H11. `full_output=True` does NOT make a scipy status readable — `disp=False`
is the load-bearing half:** with `disp` defaulted True a non-converged
`brentq`/`root_scalar` RAISES instead of returning `converged=False`, so the
False leg is an unreachable branch wearing an honest name. → L-066

**H12. The SUBJECT of your review can move while you review it.** [definition
§5 owns the HEAD and `git diff | shasum` stamps, the end-of-review predicate
re-run and the separate remediated set] The DOCUMENT under review moves too,
and a commit landing on it mid-review changes the verdict's subject before the
end stamp is read: re-`wc -l` and `git log -1` it before writing the verdict
(`[M]` a memo grew 721→879 lines mid-dispatch, adding the strongest defence of
the link I was refuting). In a SHARED tree diff only your OWN touched
files, and trust git over the harness's session-start snapshot of the branch. → L-068, L-071

**H13. Never pipe a BACKGROUND command through `grep`** — the task file then
holds only the filtered output and the evidence cannot be re-extracted. Always
`--color=no`, redirect FULL output to a file, filter afterwards. [rule:
`code-search` carries the coloured-`^FAILED` false all-green] → L-069

**H14+H15. THREE greps, three answers — run ≥2 and reconcile numerically.**
[rule: `code-search` owns ugrep, `--exclude-dir=_build` and `git grep` as
tracked-only] ⛔ The wrapper does NOT honour `.gitignore`: `[M]` same query over
`docs/` — wrapper **514**, `command grep` (BSD) **793**, `git grep` **11**. Use
`git grep` as the source-truth filter and `command grep` as the ignore-blind upper
bound; `git grep -- .claude` sweeps `plans/` and `agent-memory/` too, so restrict
both counts to the SAME trees. On a FILE-COUNT claim also exclude `__pycache__`
and `.nexus`, and prove the excluded tree untracked (`git check-ignore`,
`git ls-files <tree> | wc -l`): the inflation always argues for the conclusion
its author already reached (`[M]` 529 files, 503 in stale gitignored
`docs/_build/html_*` trees, true radius 26). → L-071, L-090

**H16. The session scratchpad is SHARED by every agent of the session — a copy or
probe at a generic path (`cp/`, `p1.py`) may be another agent's.** `mkdir -p` and
an overwrite both succeed silently, and a stale tree under the same name serves
imports the extraction did not write. ⟹ every copy and probe lives under one
fresh directory named for the review (`qa_<unit>_<hash>/`, `test -e` first), and
a background battery's default 30-min limit is raised (`timeout`) or the battery
split, then the mutated file is diffed against `git show <hash>:<path>` after a
kill. `[M]` 2026-10-24: `p1.py` printed another agent's probe output; the battery
was killed at its limit with `reference_body.py` left mutated. → L-105b

---

## I. Retired INTO the skills and rules — point, don't restate

These digest rules retired into a skill or rule clause; the ledger (rule →
clause → archive ID, with the ORPHEUS figures) is the `## Ledger` section at the
end of `lessons_archive.md`. Before writing a new rule, check it is not one of
these. `vv-principles`: A8 (#18), A12 (#17d), A13 (#17e), A16 (#34 + X4), A17
(#17g), A19 (#13, fifth disguise), A21 (#17i, #17a), the A14 claim (#29), B4
(#36), D4 (§ Bit-identity + Sig-10), E2 (§ "Log every caught bug"), E4 (Mode
8(4)), F1 (#12), F3 (Mode 12), F15 (#24c), F19 (#30), F21 and F22 (#35), H1
(Mode 8(1)), H10 (Mode 8(6)), plus the optional-metadata guard (#28), the
congruence-class ladder (#13), the adjudicating instrument (#24), Mode 11 and
test count ≠ coverage (#3, #4, H1–H5). `instrument-doctrine`: A9 (X2 census
(7)), B8 (X1 activation count), the production ∪ TEST construction census (X2
census (1)). `retirement-audit`: D12 (B.7, A.2), the A23 core (item 19).
`plan-authoring` §10 and `instrument-doctrine` X1 full-success reading: F18 (DESIGNED-RED). The qa definition: F6 (§3 "the
promised against the shipped"), make-it-RED and the behaviour-neutral retype
(§3).

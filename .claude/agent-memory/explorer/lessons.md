# Explorer — Lessons

Behavioral corrections only: "what mistake did I make exploring, and what
did I learn that improved my behaviour?" The HOW of each Nexus tool lives
in the preloaded skills (`nexus-exploring`, `nexus-guide`, and the wider
`nexus-impact` / `nexus-debugging` / `nexus-verification` / `nexus-refactoring`
family) — point there, never duplicate the workflow here. Per-campaign
`file:line` maps are archaeology; they go stale in days. A lesson earns its
place only if it changes how the NEXT exploration is run.

The cross-cutting spine: **an exploration answer is not "I found the symbol"
— it is "I found EVERY consumer the next action will touch, I verified the
premise against the current tree (not the issue text, not a frozen memory),
and I separated the durable subsystem-shape from the line numbers that will
drift."** This spine is now codified as standing directives in `AGENT.md`
Operating Principles 4–7 (blast-radius, premise-verification, git-merge-status,
durable-vs-line). L-001/L-002/L-003/L-005 below are RETAINED for forensic
value (the war-story behind each directive); the directive itself, not the
incident, governs behaviour. L-004 and L-006 stay lesson-only — they fire on
narrower question shapes (carve-verdict / probe-collapse), not every task.

---

## L-001 -- A retirement/rename blast radius = graph callers AND text grep AND direct constructors AND doc nodes

→ **Now AGENT.md Operating Principle 4** (the four-search discipline). War-story kept below for the specific misses each search catches.

`mcp__nexus__callers` / `impact` find the *graph* consumers, but a retirement
audit that stops there under-scopes — and under-scoping a retirement forces a
mid-session re-plan (the documented ~2× cost behind the proactive-explorer
trigger). The consumers a single `callers` query misses:

- **Property-reached leaves.** A method only reachable by reading a
  `cached_property` and calling `.apply` shows `callers() == 0` while still
  being live through the property. Audit the property's readers, not just the
  method's callers.
- **Bypass-trick / class-name consumers.** A test that uses an orphan via its
  CLASS NAME for a side purpose (a validation-bypass) is invisible to a
  method-level `callers`. A repo-wide grep for the `_ClassName` surfaces it.
- **Direct constructors of a guarded type.** A guard-at-the-data-source has
  blast radius = EVERY direct `Foo(...)` call, not just the factory path.
- **Doc nodes that will dangle.** A retired symbol referenced from a theory
  page leaves a broken `:ref:`. `graph_query` for the doc→symbol edge (or the
  symbol name in `docs/`) catches it so the archivist hand-off is complete.

How to apply: for ANY retirement/rename audit, run BOTH graph (`callers`/
`impact`) AND a text grep of the symbol/class name AND a constructor audit (if
a guarded type) AND a doc-node scan. Four searches, not one. (Reinforces the
proactive-explorer-before-retirement trigger; sibling to method-implementer
L-004.)

**Sharpening (two graph-blinding patterns endemic to the operator algebra).**
Two constructs make `callers()` systematically lie in this codebase, and BOTH
appeared in the W-F scope audit — when an `apply`-dispatch arm's liveness is the
question, grep is not a cross-check, it is the *primary* evidence:
(a) **runtime-aliased dispatch** — `apply = _apply_impl` with `@singledispatchmethod`
arms means the graph attributes calls to the alias, not to the per-type arm; a
`@_apply_impl.register` leaf shows near-zero callers while a whole solver feeds it.
(b) **Protocol-typed receivers** — a call `solver.compute_fission_source(...)` where
`solver: EigenvalueSolver` (a Protocol) is unresolvable to the concrete `SNSolver`
method, so the concrete method reads `callers==0` even though `power_iteration`
drives it every outer step. The liveness of a dispatch arm is decided by the
ACTUAL input TYPE at the production call site — trace `power_iteration →
solver.compute_X → op.apply(<what type?>)` by READING the chain, and let grep, not
the graph, enumerate the `op.apply` sites.

---

## L-002 -- The issue text is a stale premise; verify it against the current tree FIRST

→ **Now AGENT.md Operating Principle 5.** War-story kept below for the worked examples (diamond-coefficients / 2-D-matvec premises already landed).

Repeatedly, an audit's first deliverable was "the premise the issue describes
is STALE — that work already landed." An issue body is written at one moment
and the natural trigger for its work (a related carve) often lands it early
under a *different* campaign. Examples of the same shape: a "lift the inline
diamond coefficients" issue was already folded onto the scheme; a "2-D matvec
recomputes inline" concern was resolved by an earlier phase.

How to apply: before mapping HOW to do an issue, spend one query confirming it
still NEEDS doing. Grep the named symbol / read the current body of the named
function. If the premise is stale, the deliverable flips from "implement" to
"CLOSE-VERIFY (regression-pin + issue hygiene)" — say so up front. This pairs
with the git-authoritative discipline (L-005): code state, not the issue's
prose, is ground truth.

---

## L-003 -- Separate the DURABLE subsystem-shape from the line numbers that will drift

→ **Now AGENT.md Operating Principle 7.** War-story kept below for the home-placement detail (durable → AGENT.md durable-shape section; transient → topic file flagged with the HEAD it was current at).

Every audit I wrote mixed two things with opposite shelf-lives: the durable
STRUCTURE (what couples to what, which seam is polymorphic, which path is
canonical) and the perishable `file:line` map. The structure survives years of
churn; the line map is wrong within a sprint. A memory that fronts the line map
reads as authoritative long after it has rotted, and a future session trusts a
dead address.

How to apply: lead every finding with the durable claim ("the within-group
operator is the variadic `(L+C, S, B)`; the sweep reads `ψ.boundary.inflow` and
does NOT re-apply `R·G` internally"), then mark line numbers as
re-derive-via-Nexus, never as the headline. The durable subsystem-shape belongs
in `AGENT.md` (its "SN operator-algebra subsystem — durable shape" section is
the canonical home); transient maps belong in a topic file flagged
"line numbers current at HEAD X, re-derive if drifted" — and are deletable once
the campaign merges.

---

## L-004 -- A clean carve verdict names BOTH the retire case and the keep-as-anchor case, with the discriminator

The strongest audit verdicts were not "retire it" or "keep it" but "RETIRE-eligible
BY <discriminator>, AND here is the defensible documented-KEEP, AND the call is
the user's because it turns on a judgment the explorer can surface but not make."
The discriminator that decides it is the `coding-standards` aggressive-retirement
rule's own test: same-math-available-via-the-surviving-helper ⟹ retire (genuine
redundancy); genuine-independent-consumer-need (even a future one a named typed
leaf would serve) ⟹ keep-as-anchor is defensible. The honest counter-weight is
Cardinal Rule 2 cutting both ways: a clean typed surface a future preconditioner/
DSA would consume is an architectural asset, not just dead weight.

How to apply: when asked "does X earn its keep?", deliver the dependency surface
+ the retire-with-rewire map + the keep-as-anchor counter-weight + the explicit
discriminator, and hand the value judgment to the user. Do not pre-decide a
retirement that turns on "will a future open issue consume this."

---

## L-005 -- Git is authoritative for merge-status; a memory's "in-flight / NOT pushed" freezes mid-flight

→ **Now AGENT.md Operating Principle 6** (and the always-on `process-discipline.md` rule). War-story kept below for the SN-campaign pattern that motivated it.

Memory notes captured a campaign as "uncommitted on branch X / NOT pushed," but
nearly every SN campaign merged in a later session — the note froze the moment it
was written. A future dispatch that trusts the frozen "in-flight" wastes effort
re-deriving landed work or, worse, treats merged code as still-pending.

How to apply: NEVER trust a memory's merge-status. Reconcile every "resume X /
in-flight X" against `git merge-base --is-ancestor <hash> HEAD` (or
`... <branch> main`) before acting. Active-state in MEMORY.md should say only
what git confirms; when in doubt, the answer is "check git." (Now an always-on
rule: `.claude/rules/process-discipline.md` §"Trust git for merge-status".)

---

## L-007 -- On a branch under ACTIVE edit, re-run the census immediately before reporting

During the F2 cast-family recon (2026-07-02, `refactor/pyright-burndown`), a cast
site moved 1532 → 1552 BETWEEN two of my greps: the main session was editing the
same files concurrently (uncommitted C3 carve in flight), and part of my brief
(the scattering `apply_transpose` item) was being FIXED while I explored. This is
intra-session drift — a different failure shape from L-002's stale-issue drift.

How to apply: when `git status` shows uncommitted edits in the subsystem being
audited, (1) re-run the position census as the LAST step before writing the
report, (2) diff the uncommitted hunks against the brief's items — an item may
already be mid-fix, flipping that deliverable to "confirm the in-flight fix +
report the alternative," and (3) timestamp reported line numbers as "at final
read; tree moving."

---

## L-008 -- zsh: an unquoted separator starting with `=` aborts the WHOLE compound command

`echo ===` (or any unquoted word starting with `=`) triggers zsh's `=cmd`
expansion; the lookup fails ("== not found") and the ENTIRE command line is
aborted — including greps sequenced after the echo, silently costing a
round-trip. Quote separators (`printf 'NAME\n'` or `echo "---"`), never bare
`===`, when batching multiple searches into one Bash call.

---

## L-009 -- A dataclass-FIELD rename audit is a grep problem, not a graph problem — and the field name may be a substring of an English word

Two independent findings from the `WithinGroupSystem.resolvent`/`.gains` rename
audit, both of which change how the NEXT field-rename audit is run:

**(a) Nexus does not model dataclass fields as nodes.** `context` on
`py:class:…WithinGroupSystem` returns class-level edges (doc pages, implemented
equations, referencing functions) but the only `py:attribute:` node was
`…WithinGroupSystem.loss` (degree 2) — `.resolvent` and `.gains` had NO nodes at
all, despite ~75 consumer lines. The graph surface contributed **0 of the 75**
sites. This is L-001's "graph alone under-scopes" at its extreme: for a FIELD
(as opposed to a function/class/method) rename, **text-grep is the primary
evidence and the graph is at best a way to find the owning class's doc pages.**
Don't spend a round-trip on `impact`/`callers` for a field; spend it on
`grep -rn "\.<field>\b|<field>=|<field>:"` plus a `replace(obj, <field>=` /
`getattr` / `asdict` sweep for dynamic access.

**(b) Check whether the OLD token is a substring of a common word before
proposing any replace strategy.** `gains` is a substring of **`against`**
(a-**gains**-t) — 679 occurrences in `orpheus/`+`tests/` `.py` alone. A bare
`sed s/gains/…/g` or an unanchored `replace_all` corrupts every one. The same
class of trap: `loss` ⊂ `lossless`, `space` ⊂ `namespace`, `role` ⊂ `payroll`,
`gain` ⊂ `bargain`. It also poisons the CENSUS, not just the edit — my first
grep of the owning file reported 21 `gains` hits where 2 were `against` in
prose.

How to apply: for any rename audit, (1) skip the graph for fields and go
straight to anchored greps; (2) run one `grep -c "<newword-containing-old>"`
sanity probe — or simply grep the old token with `-w`/`\b` anchors AND without,
and report the delta as a hazard line in the deliverable. Report the anchoring
requirement as an explicit instruction to the implementer, since a rename is
usually executed by a different agent than the one that audited it.

---

## L-010 -- "Complement" ≠ "the named sibling": a two-way selector split by a signed predicate has a THIRD bucket, and whether it is populated is DATA-dependent

Asked to verify a measured claim that projector `M` "is the projector onto the outflow
subspace" (i.e. `M == P_out`), the measurement reproduced on a slab and **failed on the
production cylinder**. `M` is `I − P_in` *by construction*; `I − P_in == P_out` only when
`inflow ∪ outflow` exhausts the index set. The trace's selectors are `< -eps` / `> +eps`
over a signed projection, so ordinates with `|·| <= eps` (tangential/grazing) fall in
NEITHER — and the CYL production quadrature has 4 of 8 there (Lebedev: always). Rank 18
vs rank 6 on the same face. The original claim was true only on the geometry it was
measured on.

How to apply: when an audit hands you "X is the complement of Y, hence X == Z", (1) find
the PREDICATE that defines the split and check whether it is exhaustive (a strict `<`/`>`
pair with an epsilon never is), and (2) re-run the measurement across the **production
data** that reaches it — enumerate the real quadratures / meshes / grids, not the one the
first fixture used. The slab is the degenerate case for nearly every SN index question;
CYL (`Quadrature.product`) and Lebedev are the discriminating ones. Same family as L-006
(split the probe KINDS before proposing the collapse), one level up: split the SET
partition before accepting an identity between two spellings.

---

## L-006 -- A "shape probe" is not always a missing predicate — split boolean-presence from integer-width before proposing a typed swap

Asked to collapse N value-based `arr.shape[-1] > 1`-style probes into one typed
predicate, the load-bearing finding was that the probes split into two KINDS with
opposite fates: Kind-A pure-presence ("does this axis exist?", boolean → swap to
the typed predicate) vs Kind-B width/count ("the actual `2^d`, needed for buffer
ALLOCATION → these are honest counts, KEEP them). Proposing to delete the width
derivations would have broken allocation. A second constraint that governs such
work: a typed factor may live on the FIELD, but the inner-walk sites that do the
probing often see only a bare ndarray + `mesh.scheme` — so the "clean predicate
swap" is really a small-plumbing change, not a one-line rename.

How to apply: before recommending "collapse these probes into one predicate,"
classify each probe boolean-presence vs integer-width, and check whether the
probe site even HAS the typed object in scope. Report the verdict as
"(B) small plumbing," not "(A) clean swap," when the factor isn't reachable at
the site.

---

## L-012 -- On a "blast radius ahead of a carve" brief, run `git diff --stat` as the FIRST tool call — the carve may already be underway, and that flips the deliverable's shape

During the B3.4c periodic/face-name audit (2026-08-01, `refactor/operator-strategy-layers`)
the main session executed the carve **while I was auditing it**: 8 production
files went clean → modified (+511/−76) across the dispatch, including the
"primitive to be minted" (already minted), all five named transcription sites
(already repointed), and the `apply_transpose` defect I had just written up as
the top risk (already fixed). I discovered this by accident — an import probe
caught `__all__` listing three names the module did not yet define, a transient
mid-write state.

This is L-007 (intra-session drift) one step earlier and with a bigger
consequence. L-007 says *re-run the census before reporting*; L-012 says **run
the diff BEFORE the census**, because on a pre-carve brief the premise "X is not
yet built" is the thing most likely to be hours stale, and when it is, the
deliverable is no longer a blast radius — it is a **done-vs-remaining
reconciliation**, which is a different document.

**Sharpening (2026-08-01, #326 half-range map).** `git status --short` at OPEN is
not enough to call a file "tracked at HEAD". Three test modules I cited as
"exists, tracked" were **untracked** — the main session created them mid-dispatch,
and a `?? ` line is easy to miss in a 38-line status. **Run
`git ls-files --error-unmatch <path>` on every file you are about to describe as
landed**, and re-run `git diff --name-only` at close (a helper went clean → `+118`
under me). Then tag each reference *(untracked, in-flight)* vs *(at HEAD)* in the
deliverable — the distinction changes what the reader may build on.

**Sharpening 2 (2026-08-03, non-SN geometry census).** The drift can invalidate a
**VERDICT**, not just a line number. I wrote the gap "`roots_of_unity` has zero
production consumers" — true at the opening HEAD, **false by the closing one**:
the main session landed `rules_circle.periodic_trapezoid` mid-dispatch, and that
new rule turned out to BE the thing MoC hand-rolls (its upper half, measured
identical to 5e-16). The deliverable flipped from "extract a missing primitive"
to "consume the rule that now exists" — a different recommendation. What caught
it was **re-running the census greps verbatim at close and following one
surprising hit**, not reading the diff. So: at close, re-run the *searches whose
EMPTINESS is a finding* (a "zero consumers" / "zero hits" claim is the most
drift-fragile kind of claim there is), and when `git log <open>..<close>` shows
commits in the neighbourhood you audited, read the NEW code — it may be your
answer.

How to apply: for any brief phrased "ahead of a surgical carve / before we change
X", open with `git status --short` + `git diff --stat` (and re-run both at the
end). If the carve is underway: (1) keep the audit sections as taken — they are
the record of what the carve was walking into — but (2) add a terminal
reconciliation section that verifies each finding **by runtime probe against the
final tree, not by reading the diff** (probes caught that `SpatialWrap.is_adjointable`
flipped `False→True` while `permutes_ordinates` correctly did NOT — a distinction
the diff hunk alone rendered ambiguous), and (3) lead the report with the
still-open items, since the closed ones are now archaeology. The highest-value
finding in such an audit is the item the carve did NOT reach.

---

## L-011 -- A docstring that DELEGATES ("the sweep handles it", "whoever orchestrates it") is the highest-yield falsity shape — grep the named MECHANISM, never the named symbol

Investigating the periodic-BC claim "*the SN sweep handles the spatial wrap via
its own face-pair indexing*", the mechanism (`face-pair indexing`) did not exist
anywhere in `orpheus/sn/` — no face→face map, no `partner_face`, nothing. Three
more falsities sat in the same file family, all the same shape: prose that hands
a responsibility to an unnamed other layer. Delegation claims are unfalsifiable
by the reader (the work is "over there"), so they survive refactors that would
have deleted an ordinary wrong sentence.

Two techniques that made it cheap:

- **Grep the MECHANISM NOUN, not the symbol.** `grep "PeriodicBoundary"` returns
  40 live hits and looks healthy; `grep "face-pair"` returns 2 — both of them the
  claim itself, restated. The noun the claim invents is the fastest disproof.
- **The SIBLING METHOD that REFUSES the same thing is a free oracle.** The
  diffusion realizer raises `BoundaryError` on the identical law with the exact
  structural reason ("*couples a face to its OPPOSITE face … no slot for
  cross-face coupling*"), while SN silently realizes an identity. When one method
  in a polymorphic family refuses and another accepts, the accepting one is the
  suspect — read its acceptance, not its docstring.

Also: the claim may be **half-stale**. The brief attributed it to two files; one
had been rewritten already and now carried a DIFFERENT unsatisfied claim. Always
re-locate the quoted prose before judging it (Operating Principle 5 applied to
prose, not just issues), and report where it actually lives.

How to apply: for any "is this doc claim true?" task, (1) extract the invented
noun and grep THAT; (2) check whether a sibling implementation refuses the same
input and why; (3) check the strict-`xfail` set — a `pytest.mark.xfail(strict=True)`
row naming the gap is the project's own admission that the claim is false, and is
better evidence than any prose.

---

## L-013 -- For "what breaks if this numeric primitive changes?", SWAP IT AND RUN — a grep-classification of exact assertions is guesswork

Auditing #325 (trig-evaluated → algebraically-generated quadrature nodes), a grep
for `assert_array_equal` across the 50 consuming test files returned ~200 hits.
Classifying those by reading would have been slow AND wrong. Instead: a **pytest
plugin that swaps the primitive at `pytest_configure`** and a run of the consuming
surface. 3024 tests, ~9 min, answer = **exactly 1 failure + 1 DriftWarning**. The
audit's central number went from "~200 candidates to triage" to "2 items", measured.

Three sub-lessons that generalize:

- **The dangerous class is a FROZEN right-hand side, not an exact comparison.**
  `assert_array_equal(route_A, route_B)` computed in the SAME process from the
  same inputs is *immune* to an input perturbation — both sides move together.
  Only a comparison against a stored `.npz`/`.npy`, a hardcoded literal, or a
  hash can move. Classify by "is the RHS frozen?", never by "is the comparison
  exact?". Nearly all of the ~200 hits were route-equivalence and immune.
- **Patch every re-export, not just the definition.** `from X import f` binds a
  name; patching `X.f` misses `directional.f`, `registry.f`, and both
  `__init__.f`, plus any dataclass field that CAPTURED the function object (the
  registry's `QuadratureSpec.factory` needed `object.__setattr__`). Patch the
  module list AND the captured references, then print a confirmation line from
  `pytest_configure` so a silent no-op swap is impossible.
- **The consuming-file list from grep is INCOMPLETE — run the sibling suites
  too.** `test_dd_regression.py` never spells `.product(`; it reaches it through
  `_generate_snapshots.CASES`. That file held the ONLY moving snapshot. A
  second batch over the whole owning directories (`tests/sn/regression`,
  `tests/moc`, …) is what found it.

Two more findings from the same audit, both re-usable question shapes:

- **A guard test's FIXTURE ENUMERATION is where vacuity hides.** A test asserting
  "no shipped quadrature has a cosine in the round-off band" built its list from
  `gauss_legendre` + `lebedev` only — excluding `product`, the one family that
  violates it. The assertion was strong and the *sample* was empty of violators
  (vv Mode 7). When asked "does this guard bite?", read the parametrize/fixture
  LIST first, then the assertion.
- **Check whether the "new" hazard already exists on the sibling.** #325's ties
  looked like a new reproducibility hazard — until measuring `level_symmetric`
  (already algebraic, already exact) showed 18–216 ties per rule TODAY, with
  `np.argsort` kinds already disagreeing at LS6+. Exact symmetry CREATES ties;
  the already-exact family is the free oracle for "is this consequence new or
  pre-existing?" (Same oracle move as L-011's sibling-that-refuses.)

---

## L-014 -- To adjudicate an ALGORITHM against the literature, read the source's DERIVATION and its INDEX-DOMAIN sentence — the equation alone is permutation- and domain-agnostic

Asked whether the cylindrical α-recursion's ordinate ordering is "correct or merely a
convention" (#326), the theory page and the code agreed with each other and both
matched the published *equation* — so an equation-level check said "fine". The answer
was in two places the equation is not:

- **The source's DERIVATION names what the quantity IS.** Hébert doesn't just state
  `α_{q+1/2} = α_{q−1/2} + 𝒲μ`; two lines earlier he *defines* `α ≡ 𝒲_p·η_{q+1/2}` —
  the tangential cosine at a real boundary. That single definition converts "which
  ordering is conventional?" into "which ordering reproduces the closed form?", i.e.
  from a taste question into a decidable one. **Always read the paragraph that
  PRODUCES the equation, not the equation.**
- **The load-bearing sentence is PROSE that bounds the index range.** "Each axial
  level contains 2ℱ(p) base points in interval `0 ≤ ω ≤ π`" and "the weights are
  normalized on each level to sum to `2√(1−ξ²)`" are the two sentences that decide
  the whole issue (the level is a HALF range; ORPHEUS spans the full circle). Neither
  is an equation; a grep for `\alpha` finds neither. Grep the OCR sidecar for the
  *domain* words — `interval`, `octant`, `normalized to`, `range`, `for m = 1 … M` —
  right after you find the equation.

Two corollaries that generalize past this issue:

- **A recursion is only as meaningful as its enumeration.** Any cumulative recursion
  (`x_{k+1} = x_k + f_k`) telescopes under EVERY permutation, so its closure test, its
  sum rule, and its "step" identity are all permutation-invariant — structurally blind
  to the very ordering they appear to certify (`vv-principles` Mode 12). When a brief
  asks "does this gate adjudicate X?", check whether the gate's quantity is a
  telescoping sum before reading its assertion.
- **Find the closed form first; it is a cheaper oracle than an MMS.** The whole
  question collapsed to a pointwise identity (`α == −W·ξ` at the boundary, exact via a
  Dirichlet kernel) — a millisecond quadrature-only check, versus the proposed
  "run the L1 MMS suite under each candidate ordering". And the MMS turned out to be
  Mode-7 blind anyway (both ansatzes lived inside the symmetric sector). **When a plan
  proposes an expensive reference to settle a discretization choice, spend one query
  asking whether the coefficient has a closed form.**

How to apply: for any "is this convention or is it determined?" brief, (1) pull the
source's derivation paragraph, not its equation; (2) grep the sidecar for the index-domain
prose; (3) classify each candidate gate as telescoping-invariant / fixture-restricted /
frozen-RHS before believing it bites; (4) hunt the closed form before endorsing a
reference solve.

---

## L-015 -- To test a "this DOF is redundant, fold it" hypothesis, enumerate the FUNCTIONALS, not the algorithm — and first ask whether the fold is of the ALGORITHM or of the STATE

Asked to map a half-range azimuthal level (#326) and to try to REFUTE the framing
"the redundancy IS the bug", every attempted refutation aimed at the marching
algorithm (does alpha still close? does the redistribution term survive doubled
weights? do the specular BCs still pair?) came back CLEAN. The one real break was
somewhere the algorithm never looks: a **functional of the state**. A fold with
doubled weights reproduces every *even*-parity spherical-harmonic moment to 5e-16
and turns the *odd*-parity one from its structural `-1.3e-16` into `+2.94`.

The reusable move, in two parts:

- **The sweep is parity-blind; the analysis face is not.** A quotient by a
  symmetry group G is exact on the G-invariant sector and meaningless on its
  complement. So enumerate every INTEGRAL the code takes of the state — moments,
  currents, leakage, inner products, the adjoint metric — and classify each
  integrand's parity under G. Even-parity functionals survive the fold with
  reweighting; odd-parity ones are *out of the space*, not "inaccurate". That
  reframing is what turns a defect into a typed obligation (restrict the trial
  space / a Petrov-Galerkin analysis face) instead of a patch.
- **Split "fold the algorithm" from "fold the state" BEFORE scoping.** They read
  as the same proposal and have opposite blast radii. Folding only the MARCH and
  lifting the partners back into the full state buffer makes the symmetry hold by
  construction, kills the ordering ambiguity, AND makes the functional break
  vanish (the integrals see the lifted, exactly-symmetric state) — at the price of
  no memory saving. Folding the STATE is the memory win and pulls in the trial-
  space restriction, the partition-by-sign consumers, and every `n_dof`/Krylov
  resize. The brief usually means the first; the headline number ("2x fewer
  unknowns") advertises the second. Name the split in the deliverable.

A third, cheaper corollary from the same audit: when two candidate constructions
both satisfy the criterion the issue argues from (here, both half-range rules
reproduce the alpha closed form with the SAME constant), **the criterion does not
discriminate them — go find the predicate that does.** It was a one-line structural
predicate elsewhere in the tree (`0 < tau_raw[0] < 1`) that flipped an entire
solver route on for one candidate and not the other. Sweep the STRUCTURAL
PREDICATES the codebase already keys behaviour on, and evaluate each candidate
against them; that is where the real cost difference lives.

---

## L-016 -- A stored NUMERIC PROPERTY (an exactness/order/rank tag) is a claim: sweep it, and first establish what the SYMMETRY gives for free

Surveying the quadrature landscape (2026-08-01), the single highest-value finding
came from ~40 lines of probe: brute-force sweep the monomials and ask "what is the
LARGEST degree at which EVERY monomial is exact?", per rule, per parameter. Result:
`level_symmetric_sn` tags `degree_of_exactness = N-1` and measures **3 for every
N** — a 12-degree over-claim at S16, with a live consumer (the registry selector
returns it when you ask for degree 15). Two moves made it decisive:

- **Establish the FREE baseline before crediting the construction.** I built a
  RANDOM `O_h` orbit with equal weights summing to 4π and measured it: degree-3
  exact, fails at 4 — identical to the real rule. So the rule's entire measured
  exactness is a consequence of the `invariance_group` tag it already carries, and
  the level construction contributes nothing. Without that control I would have
  reported "degree 3, not N−1" instead of "the number is not just wrong, it is
  redundant with another field". **For any "does this property hold?" audit, first
  measure what a structurally-trivial object with the same declared symmetry gives.**
- **Read the WEIGHTS, not the nodes, when a moment claim fails.** `n_distinct_weights == 1`
  for every N immediately named the root cause (equal-weight, not the cited
  Carlson-Lathrop moment-matched construction). A node-level diff would not have.

Why no test caught it (the reusable diagnosis): every test asserted the **tag**
(`assert m.degree_of_exactness == sn_order - 1`) and every *property* test stopped
at degree 2 — inside the sector the symmetry makes free (vv Mode 7,
fixture-restricted). **A tag-pinning assert is not a property test; when a brief
says "audit claim X", grep the tests for `== <the tag>` and treat every such line
as evidence that the property is UNTESTED.**

Also from the same survey, two cheap discriminators worth reusing:
- **When a docstring says a `min()` is "conservative" over two incommensurable
  units, look for the unstated MAPPING before calling it a bug.** The product
  rule's `min(2n_mu-1, n_phi-1)` measured CORRECT and SHARP on both branches,
  because for `x^a y^b z^c` on `S²` the max azimuthal frequency is exactly `a+b`.
  The defect was the missing mapping, not the arithmetic — a different (and much
  cheaper) deliverable than "the formula is wrong".
- **A structural FLAG's docstring and its registry ENTRIES can disagree; measure
  to decide which is wrong.** `half_range_clean`'s attribute docstring said
  "Lebedev and level-symmetric are not"; the entries said LS=True, Lebedev=False.
  Measured `w(z>0)/w_tot`: LS/product exactly 0.5, Lebedev 0.33-0.43 (equator
  nodes). The ENTRIES were right. Never assume the prose is the ground truth just
  because it is longer.

---

## L-017 -- Before counting a retirement's blast radius, check whether a NON-target sibling shares the target's name; and never accept a test's self-description of what it pins

Auditing four Gauss-rule retirements, a bare `grep gauss_legendre tests/` returned
~570 lines. The true blast radius of the target was **2 files**. The other **450
lines** were `Quadrature.gauss_legendre(...)` — a *classmethod* in a sibling module
of the SAME package, spelled identically to the module-level function being retired,
and emphatically not retiring. Reporting the unanchored number would have inflated
the scope ~200x and buried the one finding that mattered.

This is L-009's substring trap one level up: not `gains` ⊂ `against` (a lexical
accident) but a genuine **namespace collision inside one package** — the factory
classmethod named after the rule it wraps. It is the NORM in a
`Facade.rule()` → `rules_x.rule_on_y()` layering, not an exception.

- **How to apply:** the FIRST action of a retirement audit is a two-number probe —
  `grep -c '<name>'` vs `grep -c '<anchored form>'`. Report the delta as a
  named hazard, and hand the implementer the anchored pattern explicitly (they
  are usually a different agent). Anchor on the call shape (`[^.]name(`) or on
  import lines, not on the bare token.

**Second half — a test's docstring is not evidence of what the test pins.** The
audited `test_..._bit_identical_to_legacy_adapter` self-declared as "the
**load-bearing contract** for the refactor: if the nodes drift even at the last
bit, the regression snapshots will silently mis-compare", and used
`np.array_equal` for emphasis. Reading the RHS's call chain showed it was
`Quadrature.gauss_legendre(n)` — which *calls the LHS function*. Same process,
same source: pure route-equivalence, immune to the exact drift it advertised
(L-013's frozen-RHS rule, but disguised by a confident docstring and by the two
sides having DIFFERENT SPELLINGS). The real characterization surface was a set of
`.npz`/`.npy` snapshots that never name the symbol.

- **How to apply:** for every test you classify as CHARACTERIZATION, resolve the
  right-hand side's call chain one hop and ask "does this move when the SUT
  moves?" Two differently-spelled call routes that converge on one implementation
  read as independent and are not. Then go find the real frozen baselines by
  `find tests -name '*.npz' -o -name '*.npy'` — a grep of the symbol will never
  surface them.

A third, cheap corroboration from the same audit: a **captured function object in
a dataclass field** (`QuadratureSpec(factory=gauss_legendre_on_mu)`) is a live
production consumer with ZERO graph edges — `callers()` reported it nowhere. L-013
already flags this for *patching*; it applies identically to *auditing*.

---

## L-018 -- Before scoping a change to a STATIC TABLE, measure which rows are still consulted; and when one tag routes through two dispatch branches, the discriminating fixture is the one that FAILS

Mapping the blast radius of parameterizing `SubgroupOfO3.Z2` by its mirror plane,
the two findings that reshaped the answer both came from measurement, not reading:

- **Half the table was dead code.** `_NAMED_LATTICE` looks like 5 load-bearing
  `Z2` edges. But `_contains` decides **finite × finite by computed matrix
  containment** and only falls back to the table when one side is *continuous* —
  so `(Trivial,Z2)`, `(Z2,Oh)`, `(Z2,Ih)` are never read. The change needed **2**
  relations, not 5. Generalizes: **a hand-written lookup table that sits behind a
  computed fast path has dead rows; enumerate the table and ask, per row, "which
  branch answers this?" before costing a change.** (Corollary for a NEW tag type:
  both `_contains` and `_check_invariance_1d` end in a bare `return False`, so an
  unhandled tag gets a *wrong-but-silent* answer — measured `O3.contains(Mirror('x'))
  = False`. Check the fallthrough, not just the arms.)
- **A tag with two dispatch branches is invisible on the fixture both accept.**
  `Z2` means "σ_z" on 3-D nodes and plane-free "`x → −x`" on 1-D nodes. Every
  shipped gate uses Gauss-Legendre / Lebedev — sets closed under **all three**
  coordinate mirrors — so σ_x and σ_z agree and the overload cannot be seen. The
  discriminator was a set that *breaks* the symmetry: embed an asymmetric μ two
  ways and the answers split (`True` vs `False`), exposing a false certification.
  **When a brief asks "is tag X consistent across its consumers?", build the input
  that FAILS the property and check whether the two routes still agree** — an
  input that passes proves nothing about which question was asked. (Sibling to
  L-016's free-baseline control and vv Mode 7.) The cheapest version: hunt a
  SHIPPED datum that already discriminates — `product(4,3)` is closed under σ_z
  and not σ_x, so no synthetic fixture was needed at all.

---

## L-020 -- The BRIEF's own timeline is a claim; and a prior audit's "structurally cannot express X" EXPIRES the moment a substrate lands

Re-auditing a 1-day-old boundary-layer audit (2026-08-03, `refactor/operator-strategy-layers`),
the two findings that reshaped the deliverable were both about **inherited claims
that had a timestamp**, and neither was in the list of things I was asked to check:

- **The dispatch brief said "what landed AFTER that audit: B3.4a/b/c".** Measured:
  B3.4a/b/c landed **2026-08-01**, the audit was written **2026-08-03 00:43**, and
  the audit *cites B3.4a and B3.4c by name in its own body*. The brief's framing
  was backwards, and had I trusted it I would have gone hunting for changes that
  could not exist, and mis-attributed `SpatialWrap`'s `is_adjointable` flip to a
  post-audit event. **One `git log -1 --format="%h %ad" --date=iso <hash>` per
  named commit, plus the target document's mtime (`stat -f "%Sm"` on macOS —
  `ls --time-style` is GNU-only), settles it in one call.** Operating Principle 5
  says verify the ISSUE's premise; this extends it: **verify the BRIEF's own
  timeline the same way.** A brief is written from the dispatcher's memory, which
  freezes exactly like an issue body does.
- **The audit's negative capability claim was half-stale in 6 hours.** It wrote
  "`symmetry.py` **CANNOT express** the periodic translation — `SubgroupOfO3` is
  origin-fixing; explicitly out of reach today." True at write time. Six hours
  later a *different* campaign step (G3) landed `RigidMotion`, an E(d) element
  carrying a translation part with a `translation_by` constructor. The tag layer
  still can't name it — so the claim is now true of `SubgroupOfO3` and false of
  the substrate, which is a materially different recommendation.

The reusable rule: **an "X cannot be expressed / X does not exist / X has zero
consumers" verdict is the most perishable kind of finding, and the thing that
expires it is usually a sibling campaign, not a change to the audited code.** So
when re-auditing, don't only re-run the emptiness greps (L-012 Sharpening 2) —
**read the NEW module any intervening commit added**, and ask the negative claim
against it directly. Here that was ~10 lines of probe (`RigidMotion.translation_by`,
`on_points` vs `on_directions`, and `_orbit_closure` fed the translation) and it
turned "out of reach" into a precise three-tier answer: the ELEMENT exists, the
TAG cannot name it, and the certifier correctly REJECTS it because a translation
is the identity on directions.

Corollary worth reusing: **when an audit reports N spellings of one concept, check
whether they differ in DOMAIN/CODOMAIN before endorsing a unification.** The
"four vocabularies of σ_a" turned out to be four different categories (deck
transformation `Γ₊→Γ₋` / constitutive kernel `Γ₋→Γ₋` / a curried factory /
a subgroup tag), and three of the four are a deliberate, documented,
sweep-schedule-load-bearing split. The genuine duplication was somewhere the
audit never looked: **three live spellings of the axis-letter→index table, beside
a docstring asserting it has "ONE home"**. Count the tiers first; the real twin is
usually the boring shared primitive, not the named types.

---

## L-019 -- Hunting a hidden TRANSFORMATION: read the chart-defining ASSIGNMENT and COUNT the partition's parts — the matrix and the docstring are the two places it isn't

Auditing "where does the angular layer rotate/reflect without naming a group
element?", the two highest-value findings were invisible to every grep I would
naturally run (`rotat`, `mirror`, `np.eye(3)`, `wigner` — all ran, all missed
them), and both were found by a different move:

- **A coordinate CONVENTION applied by variable choice is a group element with
  no matrix.** `cos_theta = mu_x` (one assignment, in a basis module) makes every
  `Y_ℓ^m` in the project the textbook harmonic composed with the 120° rotation
  about `(1,1,1)` — a real `O_h` element, in the checker's own `_octahedral_ops()`,
  and *not expressible as a tag* (`Cn(3)` is about z and measurably excludes it).
  Nothing can test it: there is no matrix for an invariance check, no adjoint,
  and a rename breaks nothing. **How to apply:** for any "find the hidden
  transformation" brief, grep the *chart-defining assignments* — `cos_theta =`,
  `= arctan2(`, `polar axis`, `_ = nodes[:, k]` — and reconstruct the implied
  frame matrix by hand, then ask `_group_elements(tag)` whether the machinery can
  NAME it. Constructing the 3×3 and testing membership is ~10 lines and decides it.
- **A partition predicate's LABEL SET is an orbit-type stratification — count the
  parts before believing the name.** `Quadrature.octants` is documented as the
  8-way sign decomposition; measured it returns **26** parts on `lebedev(17)`
  (8 chambers + 18 walls) and 2 on a slab. The `0`-component labels are exactly
  `Fix(σ_a)` — the singular set the same package computes EXACTLY elsewhere. One
  `len(...)` per shipped rule turned "ad-hoc sign classification?" into a table.

Two corollaries that generalize past this audit:

- **The tolerance-family census is a cheap, high-yield side product.** Three
  epsilons for one question (`1e-15`, `8.88e-16`, `1e-14`), all provably idle
  (measured min genuine `|cos|` = `1.57e-1`), and the one comment defending the
  first points at `_DEGENERATE_ABS_MU_THRESHOLD` — a symbol that exists **nowhere**.
  That is L-011's delegation-shaped falsity in a `#:` comment rather than a
  docstring: whenever a constant is justified by "keep in lockstep with X", grep X.
- **Measuring a docstring's claim on the DEGENERATE input finds the bug the test
  fixture can't.** "For slab GL1D only the `m=0` harmonics are non-zero" is false
  at `ℓ≥2` (measured ~0.83 in the `m>0` slots; a 4.4× reconstruction difference),
  because the slab's `(μ,0,0)` embedding makes `(cos φ, sin φ) = (0,0)` — not a
  point of `S¹` — and the on-axis guard never fires. The only `P≥2` test in the
  tree uses a 3-D Lebedev rule, where the chart is fine. **When a docstring says
  "for the degenerate/1-D case, X vanishes", evaluate X on that case.** It costs
  one probe and it is exactly where the fixtures aren't.

---

## L-023 -- "N spellings of one concept" is a SYMPTOM: find the primitive that DISCARDED the information. And a brief's named EXEMPLAR consumer is a claim to verify

Auditing "convergence has three spellings in `sn/solver.py`" (2026-08-08), three
moves reshaped the deliverable, and all three generalize to any
"one concept, many transcriptions" brief:

- **Count the spellings, then go UP one hop and ask who threw the answer away.**
  There were FIVE, not three (the brief's own count was low: a 4th differed only
  by `<=` vs `<`, a 5th was the Krylov arm). But the finding that changed the
  scope was that `power_iteration` *knows* whether it broke out of its loop or
  exhausted `max_iter` and **returns a 3-tuple with no flag**. Four production
  solvers consume that tuple; one re-derives the fact correctly, one hardcodes
  `True`, three don't try. **The N transcriptions are not N independent bugs —
  they are N callers re-deriving a fact the callee already had.** So: for any
  "why is this spelled N ways?" brief, read the CALLEE's `return` statement
  before mapping the callers. If the callee computes and drops the quantity,
  the fix is one hop up and the deliverable flips from "unify the N sites" to
  "the primitive's return type is wrong". (Corollary that decided the
  local-vs-shared scope: the eigenvalue path records only a COUNT, never the
  inner residuals — so a derived flag there is structurally outer-only until
  new plumbing lands. That constraint is invisible from the construction site.)
- **A default value can be the mechanism behind the hardcoded lie.** The type's
  own `converged: bool = True`, plus a delegate returning `True` when the
  history is `None`, means "forgot the kwarg" and "no diagnostics" both read as
  success. Grep the *defaults* of any boolean claim field before writing up its
  hardcoded sites — the literal `=True` at one call site is usually the type's
  posture made visible, and the type-level fix reds the unit tests that PIN the
  default (find those; they are the contract change made visible, and they get
  rewritten, not deleted).
- **The brief's named EXEMPLAR is a claim, like its timeline (L-020).** The
  brief said "`test_dsa_rate.py` caps `max_inner=50` on purpose" — measured, its
  helper defaults to **4000** and the `=50` sites are HEADROOM on tests asserting
  landing in 2–3 iterations. Had I inherited it, the opt-out population would
  have been sized wrong in the direction that decides the default policy. One
  `sed` of the helper's signature settled it. Same move as verifying an issue's
  premise, applied to the one datum a brief states with most confidence.

Cheapest high-leverage probe of the whole dispatch: **grep the pytest config for
`filterwarnings`**. Zero occurrences ⟹ a `warn`-by-default policy costs *zero*
test churn, which resolves the raise-vs-warn question that the rest of the audit
could only frame. When a brief asks "should this raise or warn?", one grep of the
consumer's escalation config often decides it outright.

---

## L-022 -- On a "what REMAINS of issue X" recon, the highest-yield staleness class is the campaign's own MID-FLIGHT PROSE

Mapping the #325/#326 remainders (2026-08-08) after a campaign had landed the
fixes in steps, the four falsified-prose finds all had the same shape: an
in-code note the campaign ITSELF wrote at an intermediate step, describing the
then-true remainder — falsified by a LATER step of the same campaign. Found in
one pass: a "See Also" saying the consumer is "today still welded in as
`linspace`" (the repoint landed later); a `.. caution::` saying "the remaining
half is the checker's own C_n/σ_v operators, which still [evaluate trig]" (the
checker was repointed to the shared generator two steps later); a docstring
citing a gate by a test name that never existed; a theory page still saying
"reflection partners … cached at construction" (retired). L-020 says claims
expire; this names the SEARCH STRATEGY for remainder recon: don't only re-run
emptiness greps — grep the campaign's own interim honesty vocabulary ("today
still", "remaining half", "until then", "not yet", the caution/note directives
in the touched modules) and re-verify each hit against the tree. The notes
written to be honest mid-flight are precisely the ones nobody returns to.

Also from the same dispatch, two cheap reusable probes: (a) a "did the
acceptance gate land?" check must read the gate's FIXTURE ENUMERATION, not its
assertion — the eps-gap gate existed and still enumerated only GL+Lebedev,
skipping the family (`product`) that motivated the issue (L-013's finding,
still unfixed 7 days later — report it as the remainder, citing the earlier
find); (b) when a brief hands you six "the campaign landed X @ hash" claims,
verifying them costs one read each and all six held — but the SAME session's
plan had already existence-checked its own NEXT pointer (per plan-authoring
§1), which is why; trust rises with the tree's own hygiene, never with the
brief's confidence.

---

## L-021 -- A brief's TYPE table is a claim about MATERIALIZED objects; and an all-green "does it break?" run may have measured INERT, not SAFE

Mapping the G6.3 boundary-operator binding sites (2026-08-04,
`refactor/operator-strategy-layers`), the two findings that reshaped the deliverable
were both about the *shape of the question*, not about the code:

- **The brief handed me a four-row typing table** (`γ₊ : Γ→Γ₊`, `G : Γ₊→Γ₋`,
  `R : Γ₋→Γ₋`, …) and asked "where is each constructed?". The table was correct
  MATHEMATICS and had **no code counterpart**: `law.geometry_map` returns a
  `SelfPairedDeck`/`SpatialWrap` and `law.response_kernel` a `SpecularReemission` —
  *descriptors*, not `LinearOperator`s. The realizer emits ONE operator per law with
  `G` and `R` already collapsed, so two of the four rows had **zero** construction
  sites and the "endomorphism of `Γ₋`" row was un-bindable in principle. Had I filled
  the table by finding the nearest plausible object per row, I'd have shipped a map
  that told the implementer to bind a `PermutationOperator` as `G` — which is the
  SAME body two different laws reach while declaring the mirror lives in different
  tiers, i.e. structurally unfillable.
  **How to apply:** when a brief (or a theory page) names an arrow `A --f--> B`, the
  FIRST query is *"is `f` a materialized object?"* — one grep of its accessor's return
  type. Only then map construction sites. This is Operating Principle 5 applied to a
  TYPING premise: OP5 catches a premise that went stale; this catches one that was
  never true in code, only on paper. The deliverable's headline flips from "here are
  the sites" to "two of your four rows do not exist, and here is why they can't".

- **I did the L-013 move (install the change, run the suite) and it came back
  ALL GREEN — 4 941 bindings, ~5 100 tests, zero new failures. That result was
  nearly a false reassurance.** The composition G6.3 types is spelled as three raw
  `.apply` calls, never `@`/`+`, so **no composability gate exists on that path at
  all**. Green did not mean "the binding is right"; it meant "the binding is inert".
  Those are opposite messages to an implementer: the first says ship it, the second
  says ship it AND schedule the step that makes it bite.
  **How to apply:** after any all-green "what breaks?" measurement, spend one query on
  *"could a gate have fired?"* — find the gate (here: `OperatorSum`/`OperatorProduct`
  `__init__`), and check whether the audited path routes through it. Report
  `inert` vs `verified` explicitly. A measurement that CANNOT fail is Mode-7 blind at
  the whole-suite scale, and an audit that reports it as "safe" has laundered a
  no-op into evidence.

Two cheaper corollaries from the same dispatch, both reusable:

- **A predecessor survey from the SAME campaign expires like any other audit.** Two of
  the G6.0 survey's Tier-1 "pure omission, both types exist" calls were written against
  the pre-B3.4a operator shapes and were wrong by G6.1 (`AngularAverageOperator` is
  `Γ₊→Γ₋`, not angular→scalar; the boundary `TensorProductOperator`'s space is NOT
  `TensorProductSpace(a.domain, b.domain)` because the face space ALREADY is the
  product). L-020 says a sibling campaign expires a claim; this says **an earlier PHASE
  of your own campaign does too** — re-measure the predecessor's per-item calls, don't
  inherit them.
- **The cheapest cross-face threading is usually already computed and thrown away.**
  The periodic partner face — the one genuinely hard-looking argument in the whole map —
  is *returned* by the guard the realizer already calls and *discarded* at the call
  site. Before writing "X would have to be threaded", read the helper's return
  statement.

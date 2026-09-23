# Test Architect — Lessons (hot digest)

Read whole at the START of every dispatch. ONE imperative per entry, stripped to
the correction. Every measured number, war story and `file:line` lives in
`lessons_archive.md` (sections L1–L86, war stories; L87–L88, cross-campaign
lookup tables for families 6 and 8) — open ONLY the section a pointer names.

**Cite, never restate.** The failure-mode taxonomy (AI modes 1–6, test-design
modes 7–12), the three pillars and anti-patterns #1–#36 are `vv-principles`; the
four evidence statements are `instrument-doctrine`; the retirement census is
`retirement-audit`; the denominator clauses are `plan-authoring` §2; the `-O`
scope and the guard-is-debt rule are `coding-standards`; the canonical
invocation, the tolerance contract and the math-type law rule are `vv-testing`;
the reference inventory and the XS mixtures are `AGENT.md`; the per-carve
recipes are `MEMORY.md` §3. An entry exists HERE only because no rule and no
skill carries it.

**THE SPINE** (standing: `AGENT.md` §0.5 / §0.6 / §1.5). A plan is done not when
the tests pass but when, for EVERY gate: (a) a named mutation reddens it under
`python -O`; (b) the reference is structurally INDEPENDENT of the SUT; (c) the
regime ACTIVATES the term the bug lives in.

**Maintenance.** A new lesson adds its RULE here with a `→ LNN` and its war
story as a NEW archive section — never a campaign-named block appended here. An
entry recording a GAP gains a LANDING note the day the gap closes: `L3`'s "no SN
MMS exercises `q.boundary ≠ 0`" was true when written, the §4.6 family closed
it, and the stale entry then generated a whole phase brief for work already
done (`L40a`).

---

## 0. Eight meta-lessons — most entries below are an instance of one

**M5 — THE TREE MOVES UNDER YOU.** The checkpoint protocol (`git status
--porcelain` and `git log` at the START and the END; a pristine copy before the
first mutation, `diff -q` after) is `process-discipline` ("Trust git",
"Mutation-testing an uncommitted file"), cited, not restated. What this tree
adds: `ls`/`shasum` the target before any pre-carve claim and bracket every
measurement with the mutated thing's state, because a concurrent write can
DEMOTE your measurement to a value-compared-with-itself; another agent's
deliberate mutation is indistinguishable from a production bug (re-read with
`inspect.getsource` before reporting one); a carve can land twice, three times,
or mid-sentence, and then every gate becomes MEASURED rather than predicted.
→ `L28`, `L33`, `L43a`, `L44h`, `L72a`, `L75`, `L84a`

## 1. Gates that cannot red — the shapes `vv` Mode 8 / anti-#17 do not carry

- **⛔⛔ RETIRING A DEFERRAL: never read the reason string — RUN the row.** The
  fix an exclusion CITES is usually not the fix that healed it. Method: a `-p`
  plugin no-op'ing `pytest.xfail` that ASSERTS its own installation, one run per
  row; attribute the cure afterwards, and say UNDISCRIMINATED rather than pick
  between two candidates. Riders: a reason's named BUDGET may not be a live
  knob; an UNCONDITIONAL stub whose body is only the `pytest.xfail` call must be
  SUPPLIED a body with exactly one failable statement; the healed row's doc
  claims go present-tense-FALSE. Sharpens `vv` Mode 8(9). → `L45`
- **⛔ Refusal gates: the message IS the gate.** `pytest.raises` without `match=`
  legs KEYED to the argument that triggers it is teeth-less, and a blanket "the
  message names both completions" pins a FALSE reason on the row whose defect is
  different. A `-> T | None` refusal verb collapses N guards into one value:
  isolate from the INPUT side (an input passing all guards but one) and mint
  DISJOINT fragments, asserting the disjointness once. When the fix is "raise
  the project's TYPED error instead of the builtin", `pytest.raises(<builtin>)`
  is green BEFORE and after (`BoundaryError(ValueError)`) — name the SUBCLASS
  and require the pre-fix state to be RED. → `L31`, `L35f`, `L36f`, `L75e`
- **⛔ A new guard wired AFTER an existing one: the earlier guard's inputs are a
  DISCRIMINATION row, not a negative row** — assert, on that same input, old
  fragment PRESENT + new fragment ABSENT, which pins the WIRING ORDER. When a
  new law and an existing clause share one `__post_init__`, the old clause's
  historical witness violates BOTH, so each owes a DISCRIMINATING input and the
  ORDER owes its own row. → `L43c`, `L73g`
- **⛔ A guard's REACHABILITY is decided by whatever runs before it.** Enumerate
  the guards ABOVE a refusal and the PRODUCERS of its predicate before
  re-wording it; a refusal no shipped input reaches is RETIRED, not re-worded.
  Put a typed guard at the TOP of the promoting classmethod — the caller's own
  attribute reads raise a frame earlier, so the error TYPE is the placement
  gate. The audit's fifth search: a retirement can ORPHAN a guard by removing
  its only public ROUTE while guard and witness both survive. → `L60f`, `L79d`,
  `L81c`
- **⛔⛔ `warnings.warn(stacklevel=N)` is a claim about EVERY call site's DEPTH
  and NO message gate can see it**; worse, the obvious gate is blind to half the
  class ("the attributed file is outside the package" reds for `→2`, stays green
  for `→4`). Ship TWO legs — `not is_relative_to(pkg_root)` over every entry,
  PLUS `w.lineno == inspect.currentframe().f_lineno + 1` above a DIRECT call —
  and never gate the literal `stacklevel == 3`. The structural fix beats the
  gate: hoist every emission into the public entry. → `L46b`, `L46c`
- **⛔ A TYPE-ANNOTATION widen has NO runtime witness** ("hand it the wider type,
  assert it constructs" is green before and after). Gate it with
  `tests/test_pyright_ratchet.py`, and SAY which done-when items are grep
  OBLIGATIONS rather than gates, or they read as covered. → `L59d`
- **⛔ Constructing a break-exactly-ONE-invariant mutant is a design problem** —
  `np.roll(arange(N),1)` breaks measure AND sign AND involution, and no ODD
  cycle can isolate an involution. Carry the fixture PER ROW. → `L31`
- **⭐⭐ THE ROUTE GATE, for "consumer X now reads from OWNER A, not B" — no
  value or identity gate states it.** Pose, SWAP the old owner's object for a
  decoy, require the answer UNMOVED (or MOVED, for a re-point). Three traps:
  mutating ONE consumed surface certifies one route; a DRIVER that re-poses
  internally measures nothing; a surviving cache MASKS the swap, so the gate
  needs an ACTIVATION leg. The DECOY must clear the PRODUCTION ADMISSION GUARDS
  of the arm the gate lives on, not merely discriminate — print its
  discriminating array first. → `L64a`, `L65c`, `L66b`
- **⭐ Ship the arms designed to go GREEN, and read the ones that do not.** A
  DECLARED PARTIAL NULL arm is the only instrument that can state a flagship
  gate's own Mode-12 blindness. A DECLARED-BLIND arm that REDDENS is a finding
  about the gate SET — read its red SET, which can PARTITION the suite as no
  green arm could. Rows NO arm reds may be DECLARED reference-side controls: say
  so, because their green is the LICENCE to read the neighbours as coverage, not
  coverage. → `L72f`, `L73i`, `L75a`, `L78f`, `L81g`, `L83h`
- **⛔ A declared blindness must name the RIGHT symbol and stay inside the
  problem's CONVERGENT regime**, or the fixture's convergence guard fires first
  and every red reads *"did not fully converge"*. Two more: a defect with TWO
  ENDS needs TWO arms, and the arm you write covers the end you were looking at;
  and a solver entry's STRATEGY DEFAULT decides which production branch a whole
  module ever poses — `inspect.signature` the entry and enumerate the defaults
  before claiming branch coverage. → `L78l`, `L78m`, `L78o`
- **⛔ When a "remove the step" mutation is NULL, try "corrupt the step"** —
  idempotence at a fixed point makes REMOVAL invisible and CORRUPTION loud. A
  step measured inert is a NUMBER in the commit message, never a gate. And a
  "complementary pair" of guards can be complementary in ONE variable and both
  silent in another: draw both complements and check they COVER. → `L78b`,
  `L78e`
- **⛔ Two defect classes with the IDENTICAL law residual need a STRUCTURAL leg**
  — ship the partition as two separately-messaged assertions. And a `(1+ε)`
  VALUE mutation cannot red an `is` row; a battery predicting it predicts the
  impossible. → `L84j`, `L84l`
- **⛔ The identity family: three legs that are false reds waiting.** (a) A
  done-when spelled `is` on a VALUE type — `__post_init__` normalises, so the
  constructor returns a FRESH instance; gate a value-merge with
  `==`/`hash`/`name`/`repr`/container-dedup and say `is` is NOT asserted. (b)
  `hash(a) != hash(b)` is NOT a legal "these differ" leg (a frozen dataclass
  hashes the FIELD TUPLE, not the class) — assert separation through the
  CONTAINER, and run a constant-`__hash__` arm before ANY identity carve: its
  red set is the re-pose list. (c) A value-MERGE collapses roster denominators
  silently — assert `len(set(roster)) == len(roster) == N` and grep the roster's
  docstring for the sentence the merge falsified. → `L70a`, `L72b`, `L72g`,
  `L79e`
- **⛔ Price every guard at EACH TIER, and ask which arm produces the headline.**
  A guard load-bearing at its API tier can be INERT where its end-to-end test
  lives (a later stage refuses first) — write the inert tier into the docstring.
  A design's HEADLINE consequence can be its least falsifiable one: ask whether
  its arm READS THE DATA. And before writing a control that separates two
  groups, ask whether their INVARIANT RINGS differ — coinciding rings make the
  control unwritable at every fixture, and the honest deliverable is the
  measured inertness as a NAMED blindness row. → `L70b`, `L71c`, `L71f`
- **⛔ REPAIRING a decayed gate is a different design problem from writing one.**
  Re-pose onto a REGIME-INDEPENDENT mechanism, never drive the fixture back into
  the regime; check reachability BEFORE trying to reach it; never compute the
  reference with the routine that ESTABLISHES it. The acceptance measurement is
  the PER-GUARD red table with every residual miss categorically out of scope
  and said so; re-run the AUDITOR's own harness, never a re-implementation.
  → `L28`
- **⛔ A DIMENSION count (SVD nullity, rank) is blind to a kernel ROTATION.**
  When `A = A_RR ⊕ 0_K` with `A_RR` nonsingular, a defect that fills `A[R,K]`
  leaves the nullity EXACTLY unchanged (the kernel becomes `e_t − A_RR⁻¹A_RK e_t`),
  so a law stated as `dim ker A = … + |K|` stays green. Gate the OBJECT: `A e_t == 0`
  and `(Ax)[K] == 0` bitwise, one leg per side (forward rows, transpose rows), each
  with its own tooth. And a reciprocity gate on a singular metric sees only
  `A[range, ker]`: writes into kernel rows and kernel→kernel maps are its
  stabiliser. → `L90`

- **⛔ An "adjoint == closed form" gate whose reference is built from the SUT's
  own transpose is invariant under every change of that transpose that commutes
  with the metric** (a scaled extension: `r.H = 2ι` vs a reference `2ι`). Pair it
  with a RECIPROCITY leg `⟨r x, v⟩ = ⟨x, r.H v⟩` against `apply`, from raw
  weights. The two legs have complementary stabilisers (reciprocity is blind to
  anything on the metric's kernel); only a singular member separates a Euclidean
  `.H` from the metric one when the law says they coincide. → `L91`

## 2. Harness discipline — the instrument lies before the code does

`vv` anti-#17's nine checks and `instrument-doctrine` X1 are the rule. Below:
only the ORPHEUS mechanisms they do not name.

- **⛔ A red set entirely INSIDE the new gate class is not automatically `vv`
  #17(e)'s "mirror, not a gate"** — discriminate by asking whether the consumers
  EXIST and are blind FOR A STATED REASON, told apart by the PRE-CARVE consumer
  census, never by the count. → `L73l`
- **⛔ Attribute an out-of-scope red by AUDITING THE DIFF for arithmetic**, not
  by re-running before/after (which this tree forbids): `git diff -U0 orpheus/ |
  <strip comments> | grep -E '=|def |raise |return '` is the complete
  added/removed CODE line list. And a pre-existing red must be CHARACTERISED,
  not counted — "1 pre-existing red" yields nothing; "the GL case, 1 ULP, on the
  one rule the change does not touch" yields a `-k` filter and a do-not-absorb
  instruction. → `L42g`, `L44l`
- **⚠ Two shell mechanisms that read as "the battery found nothing" rather than
  "the battery ran nothing".** zsh does NOT word-split an unquoted parameter
  expansion, so `pytest $SCOPE` passes ONE argument and collects 0 — use
  `${=SCOPE}`, an array, or a driver taking `"$@"`, and print the collected
  count. And a `nohup … &` chained AFTER an `until` loop in one tool call never
  launches when the call is killed at its timeout, while `pgrep -f <script>`
  then matches the dead shell's own heredoc and prints RUNNING — launch long
  jobs with `run_in_background: true` and ONE command. → `L63h`, `L68f`, `L80h`
- **⚠ After adding a field to a type, grep the tests for REFLECTION walkers**
  (`vars(`, `asdict`, `fields(`) — a walker over arbitrary objects sweeps the
  new field's arrays into an unrelated count and reddens for the wrong reason.
  → `L65h`
- **⚠ `ast.col_offset` counts UTF-8 BYTES; `ast.get_source_segment` is the only
  safe reader — and quantify such a hazard against the CORRECT implementation,
  never a proxy** (the proxy said 128 of 128 spans corrupt; the honest
  instrument read 0 of 128). A latent hazard with 0 witnesses is a RIDER, not a
  defect. → `L74g`
## 3. Config blindness — the ORPHEUS fixture-fact inventory

Generic rule: `AGENT.md` §0.6, `vv` §H2 / anti-#3 / anti-#4 / Mode 7. Below: the
builders and corpora that SILENTLY null a channel. Each is a measurement with a
shelf life — check it against a concrete row before trusting a green.

- **`make_mixture` nulls TWO channels**: `sig_2` defaults to all-zero and there
  is NO `sig_l` parameter (it hardcodes `SigL = zeros(ng)`), so any (n,2n) or
  (n,α) term is identically nulled and a "balanced" fixture built through it is
  IMBALANCED by exactly `sig_l`. Build `Mixture(...)` DIRECTLY. → `L1`
- **Every `xs_library` mixture (A/B/C/D × 1g/2g/4g) ships `Sig2 = 0` and
  `len(SigS) == 2`.** So any (n,2n) leg there is vacuous (manufacture it —
  `Mixture.Sig2` is a Legendre STACK), and `min(L, len(SigS)−1)` maps every
  request ≥ 1 to 1, so a row spelled "P2 vs P3" asserts `1 != 1`. The only
  non-vacuous abstract-library truncation pair is `(0, ≥1)` (a real 476 pcm
  discrimination); the 421-group library is the other option. Assert the premise
  IN the row. → `L78d`, `L82a`, `L86g`
- **Only mixture A is fissile, and EVERY fissile 0-D mixture is SUPERCRITICAL**
  (`k_inf = 1.5 / 1.875 / 1.4878`); B/C/D give `rank(A⁻¹F) = 0`. Hence a
  `solve_sn` eigen snapshot on a moderator mixture is `k = 0/abs → nan`, a
  silent dead test (reformulate as fixed-source against `φ = (diagΣ_t −
  Σ_s0ᵀ)⁻¹Q`), and a "subcritical infinite medium" row is unreachable as shipped
  — manufacture it with `replace(mix, SigP=0.4·SigP, SigF=0.4·SigF)` → `k = 0.75`
  EXACTLY (k is linear in the scale because `F` is rank-1), with 0.6 flipping
  the sign as the refusal control. → `L7`, `L84d`
- **`A − F` is NEVER positive-stable on the SN composite** (8 eigenvalues at
  exactly −1.0, the TRACE rows, at every k) while the BULK block is. A predicate
  spelled "positive-stable" refuses every SN problem; the honest one is
  `ρ(A⁻¹F) = k_eff < 1`. → `L84e`
- **The SN operator fixtures carry `placeholder_materials`** (SigS / χ / νΣf all
  zero) ⟹ `F` is the ZERO operator and its reciprocity row is `0 == 0`. Every
  `tests/sn/architecture/_config` mesh is NON-FISSILE (`solve_sn` raises
  *"leakage scale bridge is degenerate"*) and reads `|N2N·x| = 0`. "Reuse the
  existing fixture" is a hypothesis: measure it, and record WHY each neighbour
  was rejected. → `L26`, `L83`, `L84j`
- **A cylindrical `SNMesh` admits only CARRYING quadratures** (15 of 15
  `folded_product`; 0 of 20 `product`/`lebedev`/`level_symmetric`). So a
  folded-vs-UNFOLDED equivalence gate is UNWRITABLE there, every admissible rule
  IS the σ_y quotient, and the σ_y half of the cylinder's declared symmetry is
  Mode-12 blind BY CONSTRUCTION — the honest gate asserts the structural
  commitment (the solver refuses every unfolded rule). → `L68e`, `L75c`
- **A 1-D rule carries `μ_y = μ_z = 0` on every ordinate** (ERR-080), so any
  σ_y/σ_z evenness leg there reads `0.0` with IDENTITY permutations: owe every
  such leg a `perm moves k/N` VACUITY guard. Read a fold's quadrant census off
  the ORDINATE cosines, never the orbit barycentres (post-#434-R4 those are
  `P_H p`, so the mirror column is exactly zero and `np.sign` reports a
  plausible, flattering, wrong 0 of 4). → `L75c`
- **The SN regression corpus pins NO full-solve `angular_flux`** — the 14 DD
  cases pin `keff` + `scalar_flux`; `2d_octant_equivalence_*` pins angular
  per-SWEEP; the only full-solve angular wall is #448's 32 anchors.
  `material_xs_field()` is a FRESH MINT per call; `geometry_cache_for` fires
  1182× per 1-D eigen solve (intern-absorbed). → `L83k`
- **`-W error::tests.sn.regression._regression_assert.DriftWarning` is a 1-ULP
  wall but NOT absolute on this tree** (19 passed plain, 9 failed / 10 passed
  escalated). Used absolutely it is 9 false reds; used as a DELTA (the drift SET
  and each case's ULP count unchanged) it is exact and free in ~1.6 s. The
  bit-exact cases are all four fixed-source cases (both P1-aniso and the
  windowed 2-D) plus `cyl_2g_3reg_folded_4x6` — the anisotropic and windowed
  paths are the strongest free anchors in the tree. Verify both that the `-W`
  string PARSES and that it bites. → `L58c`, `L77d`, `L83k`
- **⭐ Re-run the activation question PER PHASE and PER CLAIM — the answer can
  INVERT inside one campaign**, and the two halves of one step can have DISJOINT
  activating configs (per-cell scheme dispatch: slab 80 / curvilinear 0; closure
  dispatch: slab 0 / curvilinear thousands), so neither geometry family alone is
  an acceptance set. → `L64f`
- **⛔ A branch added to DODGE a rank/carrier hazard CREATES the congruence
  blindness** (the new path runs on one carrier kind only), and a second MINT
  SITE hides on the branch where the producer does not exist — a Pattern-2 twin
  one label typo apart. Parametrize the gate over the BRANCH and make the other
  arm assert the OPPOSITE claim. → `L65d`, `L66d`
- **⛔ Before promoting an observed regularity to an assertion, run it on every
  channel the same code path serves** — the (n,2n) and elastic Legendre moments
  decay monotonically, thermal does NOT. → `L76d`
- **⛔ A pseudo-inverse round trip is `P_range(G)`, not `id`, and a corpus can
  dodge the null space entirely** — all four SN ledger fixtures carry 0
  tangential (`|Ω·n| = 0`) trace slots while a legal 2-D `product(4,4)` mesh has
  32 of 64, so a naive `== id` gate is blind on the whole corpus AND a false red
  in production. → `L67c`
- **⛔ Before pricing a change to a diagnostic, check whether the tree's pin of
  that diagnostic EXCLUDES the arm the change lives on** — one SI-trajectory pin
  refuses a WINDOWED fixture by construction while the moment metric is read
  ONLY on the windowed arm. → `L80c`
- **⛔ A functional that silently accepts the WRONG SHAPE is more dangerous than
  one that raises, and a 1-group fixture cannot tell them apart** — one consumer
  REFUSED a stored `(ng,)` flux loudly while its sibling ACCEPTED it and
  returned a different number, and at 1g the two spellings agree exactly. ⟹ a
  design handing one datum to two consumers must probe BOTH, and **≥2G is
  required for a SHAPE reason**, independently of the 1-group eigenvalue
  degeneracy. → `L86d`
- **⛔ "Re-point the space" is not plumbing when today's space has NO metric** —
  a `Field.l2` moved 41 %, not ULP. Check `space.inner_product_weights is None`
  before believing any re-point is neutral. → `L62c`
- **⭐ Subcritical SN slabs are cheap and the anchors already own one** — 2g
  fuel|moderator GL-8 4+4: `L=2.0` refl|vac → `k = 0.435195214`; `L=4.0` →
  `0.907457573` (`1/(1−k) = 10.8`, the strong discriminator); `L=8.0` refl|refl
  → `1.374233987` (the supercritical refusal witness). Three solves 1.17 s;
  three dense 160×160 pencil assemblies 0.23 s. → `L84f`

## 4. Reference, claim layer, and the proactive refutation

- **⛔⛔ A brief saying a relocated computation uses "plain / flat / simple"
  arithmetic has named TWO conventions — enumerate the candidate spellings and
  MEASURE the spread before writing the pin; the spread IS the pin's
  discriminating power and its tolerance.** The dangerous rival is the CONVENIENT
  one (it deletes a helper). ⚠ A "harmless" rival can be phase-scoped, so make
  the later phase's legitimate red a blocking ruling, not a surprise. → `L58a`
- **⭐ A "the fall-through lands on X after we delete Y" claim is a RUNNABLE
  experiment, pre-carve** — invoke the base implementation the carve will expose
  on the inputs the guard must still refuse. ⚠ Assert the DISCRIMINATOR in-test:
  `a.space == b.space` can be True (`FunctionSpace.__eq__` is `(name, shape)`)
  while `a.space is b.space` is False, so a row omitting that precondition
  degrades silently. → `L58d`
- **⛔⛔ "Solve for X on interval I" is a CONSTRUCTION only if the root is UNIQUE
  there — run a sign-change scan and report the root COUNT before accepting the
  design** (a briefed root-find had TWO roots at 4 of 9 orders, so `brentq`
  raises before it starts). Two roots ⟹ the plan owes a SELECTION RULE with a
  two-legged gate: (a) the shipped value IS the selected root, (b) the discarded
  root is exhibited and measured bad — (b) is what makes the rule a reason
  instead of a coincidence, and the only attributable leg. → `L42a`
- **⭐ When a plan's blocker is "is this small number real or is it arithmetic?",
  answer with arbitrary precision instead of escalating** — one mpmath probe at
  50/60 dps closed a blocking user ruling. Keep the reasoning ✅ ANSWERED, not
  deleted (`plan-authoring` §3), and pin the MARGIN VALUE beside the sign so a
  conditioning regression reds first. → `L42e`
- **⛔⛔ AN ADJUDICATING INSTRUMENT (one that RANKS designs) is a different object
  from a gate**; its four checks are now `vv` anti-#24, promoted from here. What
  anti-#24 does NOT carry: **a CONTINUOUS homotopy beats "rank agreement"** —
  require MONOTONICITY in `w` along `blend(w) = (1−w)A + wB`, five solves,
  falsified by one non-monotone triple; STRATIFY the ensemble (NEAR = a 2 %
  jitter · MID = the real rival · FAR = garbage) and require the threshold on
  **NEAR∪MID alone**, since a ρ over all three is dominated by the garbage split;
  the ensemble MUST contain the pair inside the stabiliser you fear. Declare each
  instrument CONSTRAINT / RANKER / DIAGNOSTIC in its own docstring — the
  graveyard died of silent promotion, and 4 of 6 dead instruments died at the
  `<1 s` solve-free pre-flight. → `L48a`, `L48b`, `L48f`
- **⛔ A re-pose can INVERT a migration gate's SENSITIVITY partition — the
  inherited `[M]` characterisation dies by being FIXED, not refuted.** The old
  anti-claim arm becomes a must-RED arm and a brand-new must-stay-GREEN arm
  appears (the *un-wiring* proof) that could not be stated before. Run both at
  BOTH HEADs and put the 2×2 in the docstring. → `L61c`
- **⭐ Before excluding a field from an identity key on DOCTRINAL grounds, check
  whether the exclusion is also MANDATORY** — the stronger, more durable gate.
  Including a `Quadrature`/`DiscreteMeasure` makes `__eq__` RAISE and `hash`
  RAISE, not merely disagree; pin the REASON with `pytest.raises` legs on the
  generator TYPES. → `L65b`
- **⭐⭐ Measure an ADJOINT/metric objection on the RANGE OF THE PRODUCER, not on
  `randn(space.shape)` — a claim about inputs the producer cannot emit is not a
  claim.** A recorded "moves 10 of 33 rows" was 5 of 33 on random draws and 0 of
  33 on a covariant moment `φ = Mψ`, because on a folded rule the σ-odd harmonics
  are identically zero at every node so `G⁺` projects them out. For any `.H`
  claim on a producer's CODOMAIN, the fixture is `producer(x)`. → `L80a`
- **⭐ A registry field asserting a PHYSICS claim can be gated at the SOLVER tier
  for ~1 s, and that is the most a gate can say** — solve a deliberately
  asymmetric fixed source and compare ψ at ordinate n with ψ at the ordinate g
  maps it to, one positive leg per element IN the group and one NEGATIVE leg per
  element outside it. ⭐ And do NOT gate a table relation you cannot DERIVE:
  record it as an observation WITH its denominator, nowhere as an assertion.
  → `L75c`, `L75`
- **⭐ A declared ONE-SIDEDNESS needs its own battery arm — GREEN on the blind
  gate and RED on its two-sided partner.** The physics bound `|Σ_ℓ| ≤ Σ_0`
  catches an inflation and is blind to a deflation; the two-sided catcher is a
  RATIO-INVARIANCE row, and the arm proving the pair is `scale**ℓ`. ⚠ Threshold
  `1 + 1e-9`, never tighter. → `L76c`

## 5. Tolerance is a claim — choose it per law, from measurement

- **State the law in the direction that IS a float theorem, and normalise a
  residual that scales with its input.** `on_points − on_directions == t` is NOT
  exact (`fl(a+t) − a ≠ t`); `on_points == on_directions + t` IS bit-exact
  6000/6000 because it recomputes the same expression — and is the stronger
  assertion. An ABSOLUTE `atol` for a residual scaling `O(ops × ‖t‖ × eps)` reds
  on correct code for large draws and is too loose for small: normalise by
  `max(1, ‖desired‖_∞)`. → `L35j`
- **⛔ A DIFFERENCE of two members of an affine operator family is NOT
  bit-identical — the MATRIX form is.** `at(σ+τ) − at(σ) == −τ·M` differs on 144
  of 200 draws (catastrophic cancellation, not an identity) while
  `at(σ).as_matrix() == A − σ·F` is `array_equal`. Same family: `Σ(Aψ) − Σ(Sψ)`
  is NOT `Σ((A−S)ψ)`, so a re-signature folding two operands into one changes an
  estimator's arithmetic — and it can hold EXACTLY on a synthetic fixture whose
  hand values are representable, so the committed row survives the carve for a
  reason that does not generalise. → `L84b`, `L85f`
- **⛔⛔ For a ROOT-FIND gate the tolerance is `noise / slope`, and the slope can
  collapse 4 orders across ONE parameter family** — so a single rtol is a false
  red at one end and a dead gate at the other (1.0 ULP at S4 → 40 653 ULP at
  S18). Derive it (`Δx ≈ evaluation noise / |f'|`), tabulate PER ROW ×10
  decade-rounded, put the arbitrary-precision value in the literal, and STATE
  what the floor leaves ungated. → `L42b`
- **⛔ "Bit-identical at the degenerate fixture" is usually 1 ULP — asserting
  `array_equal` on it reds your OWN control** (`np.cos(np.pi/2) = 6.12e-17`, not
  0). Assert "15 orders below the signal", never "the bits match"; and such a
  blindness CONTROL still earns teeth — it must red when one convention moves and
  stay green when production BECOMES the other one. → `L47g`
- **Gating a MORE-accurate implementation against the less-accurate one it
  replaces gates it against the error it exists to remove.** State the criterion
  against the arbitrary-precision value, and honestly ("within 0.57 ulp
  everywhere", never "always closer"). → `L34c`
- **⭐ An instrument's ORDER is a measured choice, and MORE CAN BE WORSE** — an
  orbit-circle trapezoid mean is exact for `n ≥ 3`, so its residual GROWS from
  2.2e-16 at n=8 to 2.6e-14 at n=1024. Ship the small `n` and SAY in the
  docstring that raising it degrades the gate, or a later session "strengthens"
  it into a false red. → `L73h`
- **⛔ A conditioning bound is the amplification AVERAGED over the integrand, not
  the argument's maximum** — `exp(-tau)` condition `max tau` is wrong when the
  argument's own rounding is non-uniform (a chord `sqrt(R^2-h^2)` at the rim) and
  the integrand concentrates there: compute the weighted mean by mpmath and write
  the law from it. And a probe's sampled radii/regime are a SAMPLE: re-measure the
  edge (`r -> R`, thick cell) before inheriting its "spectral at n=64". → `L92`
## 6. Carve archetypes — where the load-bearing gate lives, by carve shape

Reference material, not a per-dispatch rule: moved to `lessons_archive.md`
§L87 on 2026-09-21 (60+ shapes indexed by their bold name — relocation,
un-weld, type-collapse, admission/refusal, rename-a-field-and-move-its-
meaning, kernel/datum-mint, ...). Skim §L87's bold names when a dispatch's
carve matches one; open only that row.

## 7. Snapshots, generators, and exactness

- **When a completion supersedes a retired spelling, INHERIT a frozen artefact
  generated by the SIBLING law rather than regenerating** — it predates every line
  under test, so re-baseline criterion 2 holds by construction. → `L31`
- **Making a SYMMETRY EXACT is NOT a ≤1-ULP change: exactness manufactures TIES.**
  Values moved 1.06e-14; the downstream `argsort` ordering changed in 36 of 36
  configurations and the end-to-end flux by 1.008 % — twelve orders apart, under
  one justification sentence. Checklist: (1) grep consumers for
  `argsort`/`argmin`/`sort`/`unique`/`set(`/dict-keying ON the quantity made
  exact; (2) is the sort key INJECTIVE? (non-injective ⟹ the ordering was never
  determined by the physics — a latent defect to be RULED on); (3) does the
  ambiguity converge away? (4) SPLIT the commits, the ordering ruling first,
  alone. And **"the level is sorted by η" is the wrong functional** — sortedness
  is invariant under permuting equal elements. Gate the full INDEX TUPLE against
  an independently-constructed one, plus a `kind=`-invariance row (quicksort /
  stable / heapsort / mergesort must agree bit-identically), the operational proof
  the key is injective. → `L34`
- **A nearest-neighbour partner search is not automatically a bug — measure the
  MARGIN first** (separation 5.0e-3 vs a 1e-16 perturbation REFUTED the claim);
  keep the search rather than minting an index-formula twin path. And `ref[ref] ==
  id` passes on a residual-0.94 garbage map: any self-inverse pairing satisfies an
  involution law, and the only functional outside that stabiliser is the RESIDUAL.
  → `L34b`, `L33`
- **⛔⛔ A principled re-read can sit INSIDE the hard tolerance band and OUTSIDE
  bit identity — name the DRIFT TRIPWIRE's new population, not just the gate's
  headroom.** Deriving the eigen `scalar_flux` as `∫ψ dΩ` was `array_equal` 0 of
  16 cases against a band with 135× headroom (every pin holds), while the
  escalated `-W error::DriftWarning` run goes 2 → ~28 reds and STAYS there. A
  bit-identity claim in this tree faces TWO gates; say which one you mean. → `L86f`
- **⭐ When a carve only RELOCATES where a measure is stored, simulate it with
  `dataclasses.replace` pre-carve and claim `array_equal` — but MEASURE it.**
  `dataclasses.replace` round-trips a frozen canonicalizing dataclass
  bit-identically (bytes, read-only flag, `eq`/`hash`, kw-only subclass fields,
  idempotent canonicalization), so a `replace`-based field upgrade is
  `coding-elegance` Pattern 4∩2 safe and needs no hand-written constructor.
  → `L65h`, `L80d`
- **⛔ A derived space NAME is a future landmine when axis identity is
  per-SUBCLASS** — pin axis CONTENT and relative identity, never the name literal.
  ⭐ And an identity derived from a RENDERING is only as injective as the
  renderer: ndarray `repr` TRUNCATES with `...`, so two distinct long weight
  vectors give the SAME name. Derive from `.tobytes()` through a digest, state the
  float caveats (`-0.0` vs `+0.0`; `nan`), and give the injectivity gate at least
  one pair whose SHAPES are identical. → `L59e`, `L62f`
- **⛔ "Prove it does not allocate" cannot be gated by asking a densifier to
  `MemoryError`** — a 550 GB outer product is OOM-KILLED (exit 137), which fails
  the RUN, not the test. Size the gate for SEPARATION (asserted on reachable
  `ndarray.nbytes`), add an EXACT structural leg, and add a BEHAVIOURAL leg —
  "never densify" implemented by DROPPING the metric passes the first two.
  → `L59c`

## 8. Verifying a pure-math PRIMITIVE (a group / algebra type)

`vv-testing` § "A type that embodies a mathematical concept ships the test of its
defining laws" is the obligation. Reference material for the fourteen further
gate shapes (homomorphism laws, involution/affine blindness, the G-fixed
centroid, bijectivity vs match window, the coset-search inverse-direction
theorem, the stabiliser-maximality gate, ...) moved to `lessons_archive.md`
§L88 on 2026-09-21 — open it when the dispatch is a group/algebra-type gate.

## 9. Pointers

- **Two instrument vocabularies, do not mix them.** THEOREM / REFERENCE / RECORD
  (§4) is a TEST's claim kind. CONSTRAINT / RANKER / DIAGNOSTIC —
  `docs/development/evidence/lessons.md` § "L51 instrument type" — is for a
  DESIGN-RANKING instrument, and is the wrong vocabulary for a test (every
  collected pytest test is a CONSTRAINT by construction, so that partition has one
  non-empty cell).
- **V&V tagging idioms** (`foundation` must not carry `verifies`; module
  `pytestmark`; `slow` on params not functions) → `feedback_vv_tagging.md` (`L9`).
  **Cross-method agreement infra** (reuse the registry schema; `max(tol_a,
  tol_b)`; tag L1 not L4; truth values MUST trace to primary citations) →
  `feedback_cross_method_protocol.md` (`L8`). **Per-carve RECIPES** → `MEMORY.md`
  §3. **Measured scope costs and the serial-host gate** → `MEMORY.md` §2 and the
  main memory's `reference_test_execution_env`.

# Rule, skill and definition candidates from the memory distillation of 2026-11-16

Tracked by #718. Each section is one owner's `candidates.md`, verbatim. A candidate stays in its owner's digest until ruled.

---

# Owner: qa

## qa — PROMOTE and UPLIFT candidates (2026-11-16)

No PROMOTE row. The two lessons qa applies on essentially every task (stamp the
review's subject; make the gate red) are already the definition's §5 and §3.
The candidates below generalise across agents, so each goes to a rule or a
skill and not to the qa definition. Each one stays in the digest until the
orchestrator applies it; once it lands, the digest rule retires into §I. The
carrier census is
`/private/tmp/claude-501/-Users-rodrigo-git-nuclear-ORPHEUS/54d619d4-06ec-4f9f-b3c1-bc030bf79384/scratchpad/qa_distill/census.py`
(session scratchpad, may be gone), run over 20 carrier files (11 rules, 8
preloaded skills, the qa definition). Its positive control was `ugrep`: 5 hits,
in `code-search`, `plan-authoring` and `instrument-doctrine`.

### U1. Sig-10's sibling-pass discriminator is void for a single-geometry change (digest D14-r)

- **Target:** `.claude/skills/numerical-bug-signatures/SKILL.md`, Signature 10,
  "Discriminator", first bullet.
- **Proposed sentence:** "The sibling-pass test is VOID when the changed code
  path is itself single-geometry: the siblings never execute it, so their green
  carries no information; bisect between two revision-fingerprinted worktrees
  instead (print the module's `__file__` and probe a post-change symbol)."
- **Evidence:** L-069 (Task 51 CYL snapshot reds), cross-refs paragraph: "Sig-10
  (whose sibling-pass discriminator is VOID when the changed code is
  single-geometry — SLB/SPH green carries no information)". The census token
  `single-geometry` has 1 hit in the skill, the symptom line ("a single-geometry
  huge-ULP red: check the siblings"). That line prescribes the check the caveat
  voids. The digest's §I table carried this row as "retired into Sig-10", which
  was false.

### U2. A `.. implements::` / `.. verifies::` `:by:` target is a retirement surface (digest G-NX)

- **Target:** `docs/development/skills/retirement-audit.md`, item A.2, in the
  surface-to-instrument list.
- **Proposed sentence:** "A `.. implements::` or `.. verifies::` directive's
  `:by:` target naming a deleted symbol: a `nexus.directive` warning, so a `-W`
  build reddens; grep `:by:` in `docs/` for every retired name and re-point it
  in the same change."
- **Evidence:** L-105a (W2.3 S4, `0f5e1c5d`). The census tokens `:by:` and
  `nexus\.directive` have 0 hits each over the 20 carriers.

### U3. A persisted pickle is a retirement surface no constructor guard sees (digest F27)

- **Target:** `docs/development/skills/retirement-audit.md`, section D, as a new
  item, or as a rider on item 18.
- **Proposed sentence:** "A retired dataclass field and a `__post_init__` guard
  are both blind to unpickling: a value persisted before the carve loads as the
  new type with the old field in `__dict__`. check: pickle the old spelling in
  the pre-carve tree (`git archive <base>`) and load it in the carved tree; a
  `__setstate__` that refuses is the fix."
- **Evidence:** L-094 (W3 reflective cleanup, `65f1dc7a`). The census token
  `unpickl` has 0 hits. The digest rule itself says "Extends `retirement-audit`
  (no persisted-state clause)".

### U4. The session scratchpad is shared by every agent of the session (digest H16)

- **Target:** `docs/development/rules/code-search.md`, the probe-hygiene
  `check:` bullets (beside "a probe is written to a FILE under a QUOTED
  heredoc").
- **Proposed sentence:** "The session scratchpad is shared by every agent of the
  session: a copy or probe at a generic path may be another agent's, and
  `mkdir -p` and an overwrite both succeed silently. Put every copy and probe
  under one fresh directory named for the task (`test -e` first)."
- **Evidence:** L-105b (WF2.3 S5, `db2a3871`): `p1.py` printed another agent's
  probe output. The census token `scratchpad` has 0 hits. The mechanism belongs
  to the harness, not to qa, so every agent that probes is exposed.

### U5. The `grep` wrapper does not honour `.gitignore` (digest H14+H15)

- **Target:** `docs/development/rules/code-search.md`, the `docs/_build` bullet.
- **Proposed sentence:** "The ugrep wrapper does not honour `.gitignore`, so it
  and `command grep` both count ignored trees: `[M]` one query over `docs/` read
  514 (the wrapper), 793 (BSD grep) and 11 (`git grep`). For a file-count claim,
  prove an excluded tree untracked with `git ls-files <tree> | wc -l`."
- **Evidence:** L-071 and L-090. A previous qa digest claimed the opposite (that
  the wrapper honours `.gitignore`), and that claim would have certified an
  inflated count. The census pattern `gitignore|_build` hits `code-search` 3
  times, all on `_build`. None of the 3 says what the wrapper does with
  `.gitignore`.

### Considered and not proposed

- **H11 (`disp=False` makes a scipy status readable)** goes nowhere. The nearest
  carrier is Signature 8 of `numerical-bug-signatures` (a discarded `info`
  flag). There the defect is a flag that is dropped; here it is a raise that
  pre-empts the flag. Different mechanism, one recorded instance: below the
  uplift bar.
- **B-DISPATCH** goes nowhere. Its measured facts already reached the
  `workflows` rule (invariant 4). What remains (measure a dispatch by
  dispatching) is a qa review method.

---

# Owner: elegance-enforcer

## elegance-enforcer — PROMOTE and UPLIFT candidates, 2026-11-16

Each candidate goes to the user; none is applied by this pass. The digest keeps
every source entry whole until a candidate lands.

### P1. PROMOTE to the definition, §4 (the instruments), one sentence

> A design with no code yet is graded the same way: build a stand-in from
> today's objects, call today's doors on it, and count the producers and readers
> at runtime, because the design's own prose is the claim under review.

Evidence: seven digest entries apply this one move to a W5 design or an
unwritten slice, each with its own `[M]`: L-025 (a `-p` plugin counted 167 of
584 constructions imbalanced), L-027 (every placement prototyped as one numpy
chain), L-028 third part (the "bit-identical" claim prototyped on the corner
fixture), L-052 (the spec's identities instantiated on the kept value), L-053
(the builder's live signature diffed against the new type), L-054 (the plan's
Protocol written in scratch and run under pyright), L-058 (a stand-in operator
read through `.H` and `.dual()`). The definition's §4 says "verified against
the live tree" and names instruments only for an artefact that exists (L-058
records this gap in its own words). Census: `prototype` and `today's objects`
read 0 in the 15 carriers. W5 is one of my two dispatch shapes, so the move
applies on roughly every design review. The seven entries stay as their
specific probes.

### U1. UPLIFT to the `code-search` rule (or `coding-standards`), cross-agent

> `pyproject.toml`'s `[tool.pyright]` excludes `scratch/`, so a pyright probe
> file written under `scratch/` reads a clean 0 on anything, `x: int = "s"`
> included: put the probe outside the excluded trees or pass a `-p` config
> whose own `exclude` admits it, with one known error as the positive control.

Evidence: L-055 `[M]` 2026-10-25 (an in-tree probe printed 0 on everything);
`pyproject.toml` lines 141–145 state the exclusion. Census: `excluded` hits 3
carriers, none on pyright (plan-authoring EXCLUSION-IS-A-PREDICATE, a census
list; process-discipline, qa's report exclusions; nexus-elegance, a tool
column). Any agent that runs a typing probe (qa, test-architect, me) meets it;
it is the X1 silent-zero class the `code-search` rule already collects.

### U2. UPLIFT to `retirement-audit` item 9 (a labelled equation is an API)

> The `:by:` target of an `.. implements::` directive is the same API: after a
> body moves between methods and the old method survives as a shell, the
> directive still resolves, no build or `dead_references` warning fires, and the
> graph keeps the equation on the shell; grep `:by:.*<Class>\.<old>`.

Evidence: L-043 `[M]` 2026-10-10, #592. Census: `:by:` and `implements::` read 0
in the 15 carriers. Item 9 covers `:eq:` citers only. Applies to the archivist
and qa on any move.

### U3. UPLIFT to `retirement-audit` item 2 or the census list ("Members spelled without the symbol")

> An AST `Name`/`Attribute` census is blind to a quoted forward-reference
> annotation (`def f(m: "T[OpT]")` is a string constant): before any
> zero-consumer verdict also grep `["']T\b` and `T\[`.

Evidence: L-028 `[M]` 2026-09-28 (a Pattern-6 trim drafted on a false 0).
Census: `quoted annotation` hits retirement-audit item 2 once, which routes a
quoted annotation to `dead_references` for a DEAD reference; it does not warn
that a consumer census misses one. Applies to any agent's consumer census.

### U4. UPLIFT to `retirement-audit` "The mirror — landing a deferred capability stales its DEFERRAL CONTRACT"

> Two milder residues: a blocker that cleared while the deferral stays true
> (the stated reason is stale, SHOULD-FIX), and an in-diff comment describing a
> plan the same diff's later commit executed.

Evidence: L-016b and archive L-016 (2nd heading). Census: `reason.stale`,
`blocker that cleared` read 0. If it lands, L-016b in the digest retires
entirely into the mirror section.

---

# Owner: test-architect

## test-architect — PROMOTE and UPLIFT candidates (2026-11-16)

Each candidate's digest entry stays in the proposed `lessons.md` until the
orchestrator applies it; once landed, that entry retires with a pointer to its
carrier. Carrier census: `census.py` and `census_out.txt` in this directory
(18 carrier files, positive controls `ugrep` 5 hits and `rests_on` 2 hits).

### P1. PROMOTE to the definition: read the tree's state before claiming anything about it (digest 0.1, M5)

- **Target:** `docs/development/agents/test-architect.md`, §3 "Each row",
  beside "Construct and measure first".
- **Proposed sentence:** "Before any claim about today's tree, record the state
  of what you will read (`git status --porcelain`, `shasum` of the target files)
  at the start and the end of the dispatch; when the orchestrator may be writing
  the carve concurrently, measure from a detached worktree of HEAD."
- **Evidence:** M5 cites eight archive sections (L28, L33, L43a, L44h, L72a,
  L75, L84a, L129), each a measurement demoted or misread because the tree
  moved under it; L129 is a battery that measured a half-written carve. It
  applies to every pre-carve dispatch, which is the role's W3 phase.
  `process-discipline` carries the pristine copy and the start/end `git status`
  ("Trust git", "Mutation-testing an uncommitted file"); `shasum`: 0 hits and
  `detached worktree`: 0 hits over the 18 carriers.

### U1. A tie broken by `set` iteration is a draw of `PYTHONHASHSEED` (digest N.hashseed, from C674.1)

- **Target:** `docs/development/skills/instrument-doctrine.md`, X2,
  "Draw-stable statistics".
- **Proposed sentence:** "A tie broken by iterating a `set` or a `dict` of
  strings is a draw of `PYTHONHASHSEED`: break every tie on a total key and run
  the producer twice, comparing outputs byte for byte, before quoting a table."
- **Evidence:** sub-cone r2 (2026-10-25): `static_closure.py` chose a test's
  home with `max(maxima, key=share)` over a `set`, so tied homes moved between
  processes; after keying on `(share, name)`, 6 of 6 outputs were
  byte-identical (`scratch/verification_tree/subcones_r2/`). The clause lists a
  seed, a wall-clock timing and a hand-typed prediction as draws;
  `PYTHONHASHSEED`: 0 hits, `set iteration`: 0 hits over the carriers. Any agent
  writing a census or a table is exposed.

### U2. `pytest -k <word>` also matches the module path (digest 2.9)

- **Target:** `docs/development/rules/code-search.md`, the check "for a probe
  that reads another program's output or runs from another tree".
- **Proposed sentence:** "`pytest -k <word>` matches the module path as well as
  the test name, so a word that appears in the file name selects every row of
  that file: select by the test-name prefix and print the collected count."
- **Evidence:** L106 (#405 test-runtime gates repair): a battery scoped by `-k`
  ran the rows it meant to exclude. `-k `: 0 hits over the carriers. qa and the
  numerics-investigator run scoped batteries too.

### U3. A scratch script named after a stdlib module shadows it (digest N.stdlib, from C790.2)

- **Target:** `docs/development/rules/code-search.md`, the same check as U2.
- **Proposed sentence:** "Never name a probe script after a standard-library
  module (`select.py`, `types.py`, `json.py`): every script run from that
  directory imports it from `sys.path[0]`, and the failure surfaces far from the
  cause (`import subprocess` dies inside `selectors`)."
- **Evidence:** 2026-11-09, the reading pass: `second/select.py` broke
  `second/check_group.py`. `stdlib`: 0 hits; `sys.path[0]`: 1 hit, the
  `python -c` CWD clause of `code-search`, a different mechanism.

### U4. A probe script imports the main checkout through the editable install (digest 2.13)

- **Target:** `docs/development/rules/code-search.md`, the clause "`python -c`
  puts the CWD at `sys.path[0]` ... (run a probe SCRIPT from outside the
  repository and print `module.__file__` first)".
- **Proposed sentence:** "Running the script from outside the repository is not
  enough for a worktree: the editable install resolves `orpheus` to the main
  checkout, so set `PYTHONPATH` to the worktree; and a pytest run of a draft
  gate under the main repository's `scratch/` uses the main `pyproject.toml`'s
  `pythonpath`, whatever `PYTHONPATH` says."
- **Evidence:** L99 (#405 P1 step 5), L118 (numerics untangling S7), L119
  (wave 2 S2b): three dispatches met it separately. `editable`: 0 hits;
  `module.__file__`: 1 hit, the clause this sentence refines.

### Considered and not proposed

- **pyright in a worktree (digest 2.18):** generalises, but the main memory's
  `reference_pyright_lsp_rooting.md` is the pyright home and the definition's
  Return clause names `npx pyright`; the right carrier is undecided, so it stays
  memory.
- **The workflow-worktree Bash guard (digest N.wfguard):** a harness fact that
  may change with the guard; its home would be `docs/development/harness.md`,
  and one dispatch measured it. Not yet general.

---

# Owner: archivist

## archivist — PROMOTE and UPLIFT candidates, 2026-11-16

No PROMOTE: the definition (rewritten 2026-11) already absorbed the lessons I apply on every task
(the live-tree check, the gate and its redirect, the silent-role list, programmatic editing, the
end-of-pass negatives). The candidates below generalise beyond the archivist; each stays in my
digest until its target carries it.

### U1 — `code-search`: a wrapper's exit code is not the program's

- **Target:** `.claude/rules/code-search.md`, the bullet "check, for a probe that reads another
  program's output or runs from another tree" (source page `docs/development/rules/code-search.md`).
- **Proposed sentence:** "A wrapper reports its own exit, not the program's: `cmd | tee log` returns
  tee's status and `nohup cmd &` the shell's, so a failing gate reads as a pass; redirect to a file
  and read `$?`."
- **Evidence:** archive L-089 and L-114 (a failing `sphinx -W | tee` read as a pass; a background
  `nohup … &` reporting the shell). The definition carries only the archivist's Sphinx half
  ("Read Sphinx's own exit code with a redirect, never through `tee`"); every agent that runs pytest
  or a gate through a pipe has the same exposure, and `tee` appears in no rule in this sense
  (census: `tee` 8 hits; 7 are inside other words — steering, teeth, guaranteed — and the eighth
is the definition's own clause).

### U2 — `code-search`: a probe script named after a stdlib module shadows it

- **Target:** `.claude/rules/code-search.md`, beside the `sys.path[0]` clause.
- **Proposed sentence:** "A probe script is `sys.path[0]`, so one named after a stdlib module
  (`numbers.py`, `types.py`) shadows it and `import numpy` dies with a circular-import error: give
  every probe an agent prefix."
- **Evidence:** archive L-118 (2026-10-02). `shadows`: 0 hits over the 19 carriers; the
  `sys.path[0]` clause in `code-search` covers the CWD-beats-`PYTHONPATH` mechanism only. The
  failure is agent-independent.

### U3 — `instrument-doctrine` skill, X2 STRUCTURAL-DENOMINATOR: the attribute-existence ladder

- **Target:** `.claude/skills/instrument-doctrine/SKILL.md` (source
  `docs/development/skills/instrument-doctrine.md`), X2, the STRUCTURAL-DENOMINATOR shape.
- **Proposed sentence:** "Whether an object HAS a field is not `hasattr`: it is False for a
  dataclass field with no class default and for an `__init__`-assigned attribute; climb `hasattr` →
  `dataclasses.fields`/`__annotations__` across the MRO → `self.x =` across the MRO → construct
  the object and print it."
- **Evidence:** archive L-093, L-076, L-053, L-105 (four separate false negatives on doc claims of
  the form "X has no field Y"). `dataclasses.fields`: 1 hit, X2's "enumerate against the TYPE",
  which names the instrument but not its false negatives. qa and explorer make the same membership
  claims; conservative, since the always-on block is over its target and this is a skill.

### U4 — the page template (#498): the peer template for a campaign landing into one page

- **Target:** issue #498 ("the page template gets a home in docs/development", OPEN), which the
  definition names as the template's future home.
- **Proposed content:** the digest §6 shape: lead (issue, one line placing the step against its
  predecessor; commits in History) · Key Facts, 5–7 bullets, one honest-scope · labelled equations
  each with its `vv-status` · one canonical at-a-glance list-table · the rationale as a
  rejected-alternative catalogue · a numerical-evidence list-table · an honest-scope `.. warning::`
  framed as the attacker's absence · cross-refs to the predecessors.
- **Evidence:** archive L-125 (verbatim original, which put the commit chain in the lead; the
  `documentation` rule's tell "a commit list" moves it to History). A comment on #498 rather than a
  rule edit; the orchestrator decides.

---

# Owner: method-implementer

## method-implementer — PROMOTE and UPLIFT candidates, 2026-11-16

No PROMOTE: no lesson here applies on essentially every task; each fires on one carve shape.

### UPLIFT 1 — three crosswalk axes the skill does not list (from L-006)

Target: `docs/development/skills/coding-elegance.md` (the source page; `SKILL.md` is generated), section "Convention crosswalk", the sentence listing the seven axes.

Proposed sentence, appended after the seven axes: "Three more seams bit in shipped carves: ×V versus ÷V reductions (which quantity carries θ, which carries V); a mass-weighted right-hand side (`M·S`) against a source subtracted raw (the average row agrees by coincidence when `M_00 = 1`); and raw `.values` arrays against the typed field space whose shape gate the carve must widen."

Evidence: L-006's four instances, each the load-bearing bug of a shipped carve (#240's LD work); census tokens `mass-weighted` 0, `field-space` 0, `Schur` 0 across the 17 carriers. The scan seam stays in memory (narrower, LD-specific). Generalises across agents: the test-architect and the main agent write crosswalks for W3 carves. All four seams are outside the seven axes; the scan seam is LD-specific and stays in memory. If accepted, L-006 shrinks to the scan seam and the d=1 link-test check.

### UPLIFT 2 — the one-substitution design-fork check (from L-015)

Target: `docs/development/skills/vv-principles.md`, Mode 12, its check sentence.

Proposed sentence: "For a design fork, compute `candidate_A − candidate_B` symbolically and substitute the intended gate's fixture; if it vanishes, the gate cannot rank the fork, and the plan names a fixture that does not annihilate it."

Evidence: two instances in one derivation (`ld_curvilinear_shape.md` rulings 2 and 5: the cylinder's per-moment flat-flux gate reads `0 = 0`; the starting-direction seed's gap is proportional to `Σ_t ψ_in − q0`, zero at the flat-flux fixed point, and to the slope source moment, zero in every ORPHEUS fixture). Mode 12 states the design-time stabiliser enumeration but not this recipe (`subs`-style recipe absent; `DESIGN time` 1 hit, the principle). Applies to the test-architect at least as much as to this agent. If accepted, L-015 retires as COVERED.

---

# Owner: numerics-investigator

## numerics-investigator — PROMOTE and UPLIFT candidates, 2026-11-16

No PROMOTE candidate: the 2026-09-21/22 pass already moved every lesson I apply on essentially
every dispatch into my definition (§2 the question kinds, §3 the cascade steps, §4 the
instruments), and this pass found the spine (M1, M2, M5) carried there in full.

Two UPLIFT candidates, both test-design corrections that generalise to test-architect and qa. The
skill is generated: the target file is the source page `docs/development/skills/vv-principles.md`.

### U1 — L8 (restored): the causally-irrelevant-region mutation probe → `vv-principles` anti-pattern #4

**Proposed sentence** (a `check:` appended to #4, "NEVER accept homogeneous-only verification"):
"**check:** a masked O(1) defect (a source/sink index swap, a wrong material lookup) surfaces only on
a spatially VARYING field: mutate a region the answer must not depend on and assert invariance, and
classify every sibling assembly module as same-bug or different-structure-correct, since the
correct sibling is the structural ground."

**Evidence.** ERR-063 (fission χ read at the SINK node where the SOURCE node belongs): invisible on
uniform material fields, settled by a mutation of a causally-irrelevant region plus a sibling solver
of a different kernel family that already source-indexed χ. Census: `causally-irrelevant` and
`irrelevant region` 0 hits in the 11 rules, 7 preloaded skills and my definition. The 2026-09-22
retirement claimed `vv-principles` #6 carried the whole entry; #6 carries only the one-source half.

### U2 — L22·4: a reference that re-associates production's own expression → `vv-principles` § Bit-identity

**Proposed sentence:** "**NEVER** gate production against a reference that computes the same
expression in a different floating-point association (`src + (A + B)` against `(src + A) + B`): the
pass is decided by which slots happen to round alike, so re-associate the reference to match or
demote the assertion to a tolerance."

**Evidence.** `[M]` the stale-reference triage (cold: `_archive/krylov_composite_restart_and_stale_references.md`):
68 of 80 slots bit-identical, 12 differ by at most 1 ULP, and the passing sibling passed by index
luck. Census: `association` 1 hit in `vv-principles`, at #31 ("prove the float re-association
exact"), which is about sweeping seeds before claiming bit-exactness, not about the reference's own
association.

---

# Owner: cross-domain-attacker

## cross-domain-attacker — PROMOTE and UPLIFT candidates (2026-11-16, HEAD `0d68fc28`)

No PROMOTE row. The lesson applied on essentially every dispatch (measure the
premise first) is already the definition's step 0; the naming procedure is step 5;
the two standing bars carry the old M2 and M3.

The four UPLIFT rows target `.claude/skills/cross-domain-frames/reference.md`. Each
was marked "PROPOSED" in the digest's Part 3 ledger at its second independent
sighting, as the definition's Return paragraph prescribes, and none has landed:
the census (`scratchpad/cda/census.py`, 15 carrier files, positive control `ugrep`
3 + 1 hits) reads `multiplier` 0 hits in `reference.md`, `marginali` 0 everywhere,
and the `collapse` / `DEC` / `exterior calculus` hits in `reference.md` are Smell
##16 text and the A.1 differential-geometry row, read at the clause. Each row is
written to `reference.md`'s own Growth Protocol (one-sentence trigger, at least
three levers, a concrete payoff). If any was declined before, the digest row
should say so; I found no record of a decision (NEEDS).

### U1. Part C smell: collapse-morphism-blind (digest Part 3 row 3, D4)

- **Proposed smell (Part C, structural, #18):** "A marginalization treated as a
  weight = 1 average. A conserved-quantity collapse is fixed by WHAT it conserves:
  a reaction rate is averaged (`G⁻¹M`), a probability or a mass is marginalized
  (`M`); a 'weight = 1 degenerate of project' divides by the bin count and breaks
  `Σχ = 1`. First test: `Σχ = 1` after the collapse, and order non-commutativity on
  a multi-axis channel."
- **Evidence:** sighting 1, XS energy coarsening (L-010,
  `xs_coarsening_collapse_marginalize_vs_average.md`); sighting 2, a different
  problem class, the MC tally: the raw collision tally is the marginal `M` of the
  empirical measure and `flux_per_lethargy = tally / np.abs(xs.du)` is the average
  `G⁻¹M` (`orpheus/mc/solver.py:578` raw tally `tally[ig] += w / sig_t`, `:769` per-lethargy, read 2026-11-16).

### U2. Part C smell: a unary predicate that is a product (digest Part 3 row 8)

- **Proposed smell (Part C, semantic, #19):** "A property stated of one object
  that belongs to a product (parameter × scheme × chart × reduction; mode ×
  functional × symmetry group). TELL: 'the outer extent is non-affine', 'DD vs LD
  splits on the sweep strategy', 'the fundamental mode admits a positive gauge',
  with no scheme, chart or group named. First test: evaluate the predicate per
  factor value (the second difference per scheme in the parameter's natural chart;
  the gauge predicate over all poles on an asymmetric and a symmetric fixture)."
- **Evidence:** four sightings, three classes: transverse-coupling order and
  affinity in a width (SN, `d5_trait_and_mms_frames.md`,
  `shift_ontology_taxonomy_frames.md`); the diffusion P1 face closure, a non-SN
  scheme, rational in the slab dilation where SN-DD is an exact pencil in 1/λ;
  mode admissibility, `[M]` 2026-09-28: a positive section on 20 of 20 and 240 of
  240 modes on the asymmetric body, on the reflection-even modes only on the
  symmetric body (D19, `fundamental_mode_gauge_attack.md`).

### U3. Part A.3 row: multiplier-algebra embedding

- **Proposed row:** Frame **Multiplier algebra / multiplication-operator
  embedding**. Trigger: scalar "coefficients" attached to operator terms that are
  really fields per region, group or cell. Levers: the embedding `f ↦ M_f` of
  `L^∞` into the operator algebra; constant sections as the scalar sub-case; the
  commutant of the multiplication algebra (an operator commuting with every `M_f`
  is itself a multiplier; lever added in this pass, not from a sighting); a one-parameter family whose parameter multiplies a
  diagonal is a coefficient family, not a new summand. Payoff example: XS fields
  promoted from per-term scalars to multipliers; the question space's cell
  coordinates are the multiplier directions; the curvilinear LD fibre metric is an
  SPD(2)-valued multiplier and the redistribution term the multiplier `(d−1)/r` in
  the Gauss-nodal basis.
- **Evidence:** three sightings in three problem classes (carrier typing,
  `coefficient_field_promotion_frames.md`; question posing,
  `posing_ontology_clean_attack_frames.md`; frame metric,
  `frame_bridge_pullback_riesz_fibre_metric.md`). `vv-testing` already lists the multiplier
  algebra's laws as gates (`M_f M_g = M_{fg}`, `M_1 = I`); the row adds the trigger, not the laws.

### U4. Part A.1 row: DEC grading with a per-axis Hodge star

- **Proposed row:** Frame **Discrete exterior calculus, per-axis degree**.
  Trigger: a per-axis "density vs integrated" or "degree" label on a
  tensor-product carrier. Levers: form degree per axis written as the pair (native
  degree, stored degree); the per-axis Hodge star as the storage change; the
  restriction (pullback) to a face as the boundary trace; the pairing of a
  k-cochain with an (n−k)-cochain (lever added in this pass, not from a sighting). Breaks: an upwind closure is a contraction, not
  a star; no complex on harmonic or collocated axes. Payoff example: the SN face
  flux is a 1-cochain and the boundary trace its restriction, so storage locus and
  variance stop being one integer.
- **Evidence:** two sightings in two classes (face-flux storage locus,
  `field_role_typing_faceflux_frames.md`; carrier pairing/variance 2026-11-12,
  `true_ends_chart_variance_frames.md`).

---

# Owner: explorer

## explorer — candidates, 2026-11-16

Two candidates, both for the user. Neither is applied by this pass.

### P1. PROMOTE, from L-046, into the definition §3

Target: `docs/development/agents/explorer.md`, section 3, the sentence that
closes an audit ("close it by re-running verbatim every search whose
emptiness is a finding").

Proposed sentence, appended to that one: "and by `git status --short` on the
audited scope, because a re-run search sees changed content and not a new
file."

Evidence. The definition takes `git status --short` at OPEN only
(`status --short`: 1 hit in the 15 carrier files, in that opening clause).
The lesson's case: a coverage map was built while a peer wrote the successor
gates; line numbers shifted about 40 lines in one file, and a new untracked
gate file appeared after the census, which no re-run of an earlier search
could report. The retired meta-lesson M-2 held the sibling check
(`git ls-files --error-unmatch` on every file called landed: 0 hits in the
15 carriers), and the definition did not take it over when M-2 retired on
2026-09-22. The check applies to every audit on a tree that a peer is
editing, which under parallel dispatch is most of them. If the user declines,
L-046 stays in the digest as proposed and nothing is lost.

### U1. UPLIFT, from L-048, into the `code-search` rule

Target: `docs/development/rules/code-search.md`, the check that already says
of `git grep -E` that "`\b` is the thing that does not exist (POSIX ERE has
no word boundary)".

Proposed sentence, appended to that check: "The same holds for `\w`: `git
grep -E 'class \w*Space'` prints nothing over a tree where `class
[A-Za-z_]*Space` prints the classes; `git grep -P` honours both."

Evidence, `[M]` 2026-11-16 at `0d68fc28`, three commands over `orpheus/`:
`git grep -cE 'class \w*Space' -- orpheus` printed nothing;
`git grep -cE 'class [A-Za-z_]*Space' -- orpheus` and
`git grep -cP 'class \w*Space' -- orpheus` each printed the same files
(`orpheus/derivations/continuous/characteristic/system.py` first). The
original case (2026-10-30): a class-name census spelled
`class \w*(Covector|Dual\w*)` read 0 while `DualSpace` existed. The rule's
existing clause covers `\b` only (`git grep`: 3 hits in `code-search.md`,
none naming `\w`). The failure is the rule's own class, a confident empty
answer, and it is not specific to this agent. For this agent the uplift
does not retire L-048: the rule's body does not reach an `omitClaudeMd`
agent, so L-048 stays in the digest marked `[rule body]` unless the
`brief:` line of `code-search` is extended, which this candidate does not
ask for (the line is already long).

---

# Owner: literature-researcher

## literature-researcher — candidates for the user, 2026-11-16

Two candidates. Both are for the user to rule on; neither is applied by this
proposal. No lesson is proposed for an always-on rule.

### U1 — UPLIFT: correct two sentences of the `research` skill

**Target.** `.claude/skills/research/SKILL.md`, section "Tier 0 — the local
OCR'd library", first paragraph.

**The sentences as they stand.** "`scratch/literature/` is the user's
actively-maintained PDF folder (ALL Nuclear Science & Engineering volumes are
local). Every PDF has a Mistral-OCR **sidecar** at
`scratch/literature_ocr/<stem>.md` plus a raw cache `<stem>.mocr.json`."

**Why both are false as written** (`[M]` 2026-11-16, Python `pathlib` globs
from the repository root):

- `scratch/literature/` holds 130 PDFs. The complete Nuclear Science and
  Engineering run is not there: it is a directory of per-volume zip bundles
  outside the repository (760 `Vol_*.zip` files), described in the memory file
  `user_nse_volume_archive.md`.
- `scratch/literature_ocr/` holds 77 sidecars (75 with a `.mocr.json` cache,
  2 named `*.textlayer.md`). Matched by file stem, 74 of the 130 PDFs have a
  sidecar and 56 have none (control: the Reed-Lathrop 1970 PDF, whose sidecar
  was read in this pass, is among the 74).

**The cost of the sentences.** The skill's first search step is a search of
the sidecars. A reader who believes both sentences reports an NSE paper as
"not local" after a search that could see neither the zip archive nor the 56
PDFs without a sidecar. This is the false negative that the retired lesson
L-012 records twice (2026-08-12 and 2026-08-25). The definition already says
"A PDF without a sidecar is OCRed once", so the definition and the skill
disagree.

**Proposed replacement.** "`scratch/literature/` is the user's
actively-maintained PDF folder. The complete Nuclear Science & Engineering run
is kept outside the repository as per-volume zip bundles; search their file
listings before calling an NSE paper not local (the agent's memory file
`user_nse_volume_archive.md` gives the path and the search). A PDF has a
Mistral-OCR **sidecar** at `scratch/literature_ocr/<stem>.md`, with a raw
cache `<stem>.mocr.json`, only once `tools/ocr_literature.py` has been run on
it; list the PDFs without a sidecar before trusting a sidecar search that
finds nothing."

### P1 — PROMOTE: one sentence for the definition, section 2 "Provenance"

**Proposed sentence.** "A section or equation number inside a correct citation
is a separate claim: before opening a page, confirm on the sidecar that the
equation falls in the section named and that the section treats the geometry
named."

**Evidence.** Lesson L-003, sharpening of 2026-08-11: a brief cited Hébert
§3.9.4 and Eqs. 3.418–3.439 for the cylindrical closure; the book and chapter
were right, the section is the sphere's, and the cylinder is §3.9.3. Two
searches of the sidecar (section headings, then `\tag{<number>}`) showed it
before any page was read. No carrier states the check: the token `subsection`
has 0 hits in the 13 carrier files, and `tag\{` has 1 hit, in the skill's
description of the sidecar format (census by Python `re`; positive control
`ugrep`, 4 hits).

**Why it is a definition sentence and not a lesson.** It applies whenever a
brief hands over a section together with an equation number, which most
extraction briefs do. If the user declines, the check stays as L-003 in the
proposed `lessons.md`, where it already is.

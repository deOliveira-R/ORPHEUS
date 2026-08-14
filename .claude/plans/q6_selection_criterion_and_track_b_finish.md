# Q6 + finish track B

> **Plan of record, 2026-08-14.** Both threads reconciled against HEAD
> `0f5ca91c` by `explorer`. User rulings recorded inline. *(This file previously
> held the #364/doc-fix plan — that work is COMPLETE and is overwritten here.)*

---

## ⏹ BOTH THREADS COMPLETE — 2026-08-14. Q6-A…Q6-E all landed.

**Re-anchor from `git log`, not from this file's older sections.** Branch
`refactor/track-b-remainder`, **10 commits** ahead of `main`, ⚠ not merged —
verify with `git merge-base --is-ancestor <hash> HEAD`.

| commit | what |
|---|---|
| `9ba35a8f` `05d35a93` `65efbb52` `66942848` | Thread B (see the older compaction point below) |
| `e17116f9` | **Q6-A** — 9 present-tense falsehoods in the quadrature layer |
| `2cb3c960` | **Q6-B(i)** — `expected_node_count` retired; stage 4 reads `measure.n_points` |
| `37f3f880` | **Q6-B(ii)** — the Lebedev order set is ASKED of scipy |
| `626dc855` | **Q6-C** — the V conjunct checked WHOLE (reference + degree) |
| `b02dd536` | **Q6-D + Q6-E** — registry frozen; inverter asks for the degree |
| `83aa324e` | the theory page + `coding-standards`' cross-reference clause |

✅ **Wide gate GREEN** — `[M]` 2026-08-14, `python -O -m pytest -m "not slow"`,
SERIAL: **9425 passed, 0 failed**, 19 skipped, 227 deselected, 76 xfailed,
**1:03:29**. Counts reconcile exactly with `tests._harness.audit`'s 9747
collected (`9425 + 19 + 227 + 76 = 9747`). The 133 warnings are the #340
convergence machinery firing in tests that deliberately starve a solve.
⟹ **the branch is ready to `--ff-only` merge.**

**Filed, not built** — each with a measured body: **#369** (the three-spelling
`ReferenceMeasure` identity), **#370** (`folded_product`'s two registration
gaps), **#371** (`spec.parameters` is a dead field), **#372** (the 27–32 s
selector cost).

### The three findings most worth carrying
1. ⭐⭐ **Q6-B's blocker DISSOLVED rather than needing a patch.** The second
   Lebedev transcription existed *only* to serve `expected_node_count`, which
   computed a number the selector already held (built at stage 2, formula
   consulted 45 lines later at stage 4). Retiring the twin left one copy to
   replace. "Clean before extending" was not stylistic ordering — it is what
   made the fix a one-copy change.
2. ⭐⭐ **A guard is not a gate.** Three arms were added to the V conjunct; `[M]`
   two bit immediately and the **`claim is None` arm was inert** — disabling it
   left all 54 gates green. Caught only by mutating each arm separately. The
   same discipline surfaced that 2 of 4 rows in the argmin gate are vacuous
   (singleton survivor sets on slab/sphere).
3. ⭐ **The tabulate-vs-ask ruling was in the file, applied HALF-WAY.**
   `_ls_sn_invert` discovered the *frontier* by attempting construction
   (vindicated when #337 moved it with zero edits) and still *tabulated* the
   **degree**. That is why the `N-1` prose went stale AND why the family was
   priced up to 2× too high. Both halves now ask.

---

## ⏸ COMPACTION POINT — 2026-08-14. **Thread B is COMPLETE. Thread A (Q6) is NOT STARTED.**
⛔ **SUPERSEDED** by the block above — Q6 is now complete. Kept for its Thread-B
detail, which is still accurate.

**Re-anchor from THIS FILE + `git log`, never from a conversation summary.**

### Where the tree is
Branch **`refactor/track-b-remainder`**, 4 commits ahead of `main`, **not
merged, not pushed** — ⚠ a point-in-time claim; verify with
`git merge-base --is-ancestor <hash> HEAD` before trusting it
(`process-discipline`).

| commit | what |
|---|---|
| `9ba35a8f` | `docs:` six present-tense-false claims (the #350 ruling's byproducts) |
| `05d35a93` | `docs(plans):` four campaign plans reconciled + three `plan-authoring` clauses |
| `65efbb52` | `refactor(numerics):` **B3.3** — `IncomingOrdinateMaskTensor` retired |
| `66942848` | `test(sn):` **B3.5** — RG-3b → B4; battery homed; its N9 leg repaired |

⚠ **The wide gate is OWED before any merge**: `python -O -m pytest -m "not slow"`,
SERIAL, budget **≥90 min**, via `Popen(start_new_session=True)` + a persistent
`Monitor` writing to a LOG (a harness background task dies at ~30 min). Tiers run
so far, all green: `tests/numerics` 2306 · `tests/sn/operators`+`tests/geometry`
+phase-C 1973 · `tests/geometry` 834 · Sphinx `-W` clean.

### What Thread B delivered
B3.3 and B3.5 landed; **B4–B7 were reconciled and REPORTED, not built**, per the
user's ruling. Filed: **#367** (B4 — the generic `R ∘ G` dispatch; the plan's
"`OperatorProduct` used by no boundary operator" sub-claim is now FALSE),
**#368** (B5's remainder — retire `ReflectiveBoundary.albedo`; the rest of B5 was
dissolved by G6.3), and a reconciliation comment on **#189** (B6 live; B7
scope-split, its MoC/MC/CP rows out of scope under the priority law).

### The two findings worth carrying
1. ⭐ **B3.3's audit predicted 6 test rewires; 5 were already covered** on
   `TraceRestrictionOperator`. Checking first saved writing five duplicates —
   *and* avoided pinning the MASK's message fragments (`"out of range"`,
   `"duplicates"`) where the successor raises `"must lie in"` / `"must be
   unique"`, which would have passed for the wrong reason.
2. ⭐⭐ **Homing the mutation battery exposed that its N9 leg had been silently
   INERT** since G6.5 moved `to_local` onto the space — a bare class-attribute
   assignment *creates* a missing attribute rather than failing. `[M]` N9 read
   *1167 passed, identical to CONTROL* (→ "the gates are blind") and now reads
   **42 failed**. Guarded with a `raise`, not an `assert`, because `-O` strips
   asserts.

### ▶ Where work resumes
**Thread A — Q6**, designed below and untouched in the tree. Start at **Q6.0**,
the exactness-reference axis, which is the PREREQUISITE that makes Q6.4 safe.
Nothing in `orpheus/numerics/quadrature/` has been modified.

⚠ Before designing to any Q6 step, re-read its **⛔ premise corrections** —
`invariance_group` is not a field, item 4 is unsafe as written (`[M]` 0.696),
and the two `ParameterSpec` arms want opposite things.

---

## Context — why these two, now

User direction: **tackle Q6, and finish track B.** Both are campaign
remainders, and the reason for deferring #364 applies in reverse — a large
reshaping campaign is imminent, and these are the *unfinished* work that should
not be carried into it.

⚠ **Naming hazard.** "Q6" and the B-substeps are **internal** task numbers.
`[M]` GitHub **#30 is "Homogeneous: Sensitivity of k_inf to boron
concentration"** — unrelated to Q6. Never resolve one of these from a bare `#N`
(`plan-authoring` §6).

⛔ **Do not read `.claude/worktrees/nexus-workspace-wiring/`.** It is a separate
checkout carrying pre-`e7d44f3c` copies of `registry.py` (*with* the retired
`invariance_group`) and ~60 stale duplicate mask hits. It is how the stale field
list reached the Q6 task text, and it will inflate any blast-radius grep.

---

# Thread B — finish track B

**User ruling:** land **B3.3** and the **B3.5 remainders**; then **reconcile
B4–B7 and report** — do not design or build them.

## ⭐ Premise corrections (all `[M]`, and they change the work)

1. **B3.5's gate half is ALREADY SHIPPED.** RG-1…RG-5 landed at **B3.2
   `7f02de15`** in `tests/sn/operators/test_sn_boundary_operator.py` (RG table
   verbatim in the module docstring `:30-70`; implementations `:185-845`). Task
   #15's subject is wrong. **Two remainders only** — see B3.5 below.
2. ⭐ **RG-3b belongs to B4, not here.** It is 4 live `xfail(strict=True)` at
   `:745`. The crosswalk (`:418-421`) parked it in "B3.4, the one place the
   guard belongs — the phase that restructures around `R ∘ G`". **B3.4a/b/c
   landed and did NOT restructure around `R ∘ G`; that is B4.** ⟹ under the
   user's ruling it **stays xfail**, and this campaign's job is to *re-attribute
   it in place*, not to build it.
3. **B3.3's "false docstring" clause is already discharged** at `b11a2ce3`: the
   "projection onto the outflow subspace" claim was replaced by a `.. warning::`
   refuting it, plus a `.. deprecated:: B3.2` block reading *"This class has no
   production consumer and retires in B3.3."* Nothing to correct — it deletes
   with the class.
4. **The suite is 13 tests, not 10.** `[M]` `grep -c '^def test'` = 13. The task
   *and* `docs/theory/foundations/boundary_conditions.rst:1434* both say ten.

## B3.3 — retire `IncomingOrdinateMaskTensor`

**Goal.** The zero map `Γ₊ → Γ₋` is the only spelling of vacuum's realization,
and no live text claims otherwise.

`[M]` **Premise confirmed twice (open and close): zero production construction
sites** — `grep 'IncomingOrdinateMaskTensor(' orpheus/` returns only the `class`
statement. Replacement path: `orpheus/sn/boundary/realizer.py:908-912`
(`_outflow_restriction` → `_narrowed_zero_operator`), probed to return
`ZeroOperator` on both space hooks.

### Three-search audit — results, not a to-do

- ⭐ **Graph finds what grep cannot: retiring the class deletes a doc→code
  `implements` edge on `math:equation:bc-extraction-loss-operator` (eq. 10).**
  That is a V&V provenance link. No grep and no build surfaces it — re-point or
  retire the edge deliberately.
- ⛔ **`orpheus.numerics.operator` is NOT Sphinx-rendered.** `[M]`
  `objects.inv` has **zero** entries for it, and the 8 mask xrefs render with
  **no `<a>` wrapper — they are already dangling plain text today**. ⟹ the
  retirement emits **zero warnings at any severity, `-n` included**. Grep is the
  only gate; an unchanged warning count proves nothing.
- **No marker migration owed.** `[M]` the suite carries 13 × `@pytest.mark.l0`
  and nothing else — no `catches(...)`, no `verifies(...)`; absent from
  `error_catalog.md` and from `tests/_harness/`.

### MUST-FIX, split by tense (16 past-tense sites STAY untouched)

⭐ **Three sites are present-tense-FALSE *today*** — false since B3.2, 13 days,
missed by `b11a2ce3`'s own doc-repair pass. These are bugs independent of the
retirement:
- `orpheus/geometry/boundary/vacuum.py:28-32` — "The SN realisation **is** an
  `IncomingOrdinateMaskTensor` that zeroes the per-face inflow ordinates only".
- `orpheus/geometry/boundary/_realizer.py:14-16` — same claim, in the module
  docstring that *teaches the three-layer split*.
- `tests/sn/sweep/core/test_phase_c_gates.py:551-553` — same claim, in a
  `@pytest.mark.foundation` gate's docstring.

Then the retirement proper: `operator.py:179` (`__all__`), `:2779-2900` (class +
3 messages), `:2568` (a row in **`TraceRestrictionOperator`'s own** docstring
inventory), `:2801-2812` (the deprecation promise);
`orpheus/geometry/boundary/_bound_compat.py:37-38` and `__init__.py:573-574`
(present-tense inventories); `boundary_conditions.rst:1431-1434` (future-tense
promise **and** the wrong count); and one paraphrase-only site,
`tests/sn/operators/test_snmesh_realizer_wiring.py:101` (a stale section banner).

### Test migration — 13 classified

- **DELETE (6):** rows 1, 2, 5, 6, 7 — masking semantics, "non-inflow
  untouched" (⚠ built on the **two-way** assumption RG-2 refuted), degenerate
  parameters, an in-process route-equivalence "pin", and a generic `M @ I` law.
- **REWIRE to `TraceRestrictionOperator` (6):** rows 8 (axis-parameterised
  action), 9/10/11 (guards — ⭐ the successor's triple is **strictly stronger**:
  range + unique + sorted; pin the *short* fragments `"out of range"`,
  `"duplicates"`, `"1-D"`), 12 (non-aliasing), and **13** — the
  adjointable-not-invertible quadrant the task names.
- ⚠ **Check before claiming lost coverage:**
  `tests/numerics/test_trace_restriction_operator.py` already ships the defining
  laws (`γι = I`, `ιγ = P` idempotent+symmetric, `ι = γᵀ` materialised,
  three-way partition), so rows 3, 4, 6 are **largely already covered there**.
- **The second live construction site**
  (`tests/numerics/test_operator_capability_predicates.py:27,65,189`) — move the
  `_LEAVES` representative and the `_CONTRACT_ROWS` row onto
  `TraceRestrictionOperator`. `[M]` it is not yet in that file, so the slot is
  free. ⚠ Per `coding-standards`, **verify the row still MOVES** — its teeth come
  from `assert_inverse_adjoint_contract`, and `STRUCTURAL_ABSENT` must still be
  the right verdict.

## B3.5 — the two real remainders

1. **Re-attribute RG-3b in place.** Its 4 `xfail(strict=True)` reasons still say
   "delete this marker when B3.4 lands the guard". Correct them to name **B4**,
   with the reason (B3.4 did not restructure around `R ∘ G`). No code change.
2. **Give the hand-written mutation battery a discoverable home.**
   ✅ **RULED 2026-08-14: move it now as a tracked plugin; file the CI meta-test
   refactor as its own issue.** Move `scratch/b3_2_mutations.py` to a
   discoverable tracked home carrying a README that names the invocation and all
   10 hazards, plus a marker so it is findable by search rather than by memory.
   This closes the evaporation risk *immediately* and preserves the mechanism
   that is **measured to work**.

   Rationale the ruling encodes: the in-process monkeypatch is not an accident of
   convenience — it is what makes the battery *correct* (it never touches a file,
   so it cannot trip the never-`git checkout` rule). Converting it to run under
   pytest is a separate, riskier job, and ⛔ it must **not** be implemented as
   pytest-inside-pytest: `vv-principles` #17 records that monkeypatching the
   parent while a child re-imports a clean module reads **GREEN for every
   mutation** — a mutation harness failing in the safe-looking direction. The
   honest refactor lifts each RG assertion into a callable helper; that is the
   filed issue's content, not this campaign's.

### The battery, measured
`scratch/b3_2_mutations.py` — 263 lines, **TRACKED** at `7f02de15`, an
**in-process pytest plugin** (monkeypatches `SNBoundaryOperator._reflect_trace`
/ the realizer / `to_local`; touches no file on disk, deliberately, per
`process-discipline`'s never-`git checkout` rule). **10 mutations, not 31** —
N1–N5, N7–N9, M1, M2 — each a plausible transcription of a real B3.2 hazard,
plus a **positive control** (`_bite_fingerprint`, `:217-263`) so a silent no-op
patch is impossible.

`[M]` **it still runs and still bites at HEAD:**

| leg | result |
|---|---|
| CONTROL | 1167 passed, 1 skipped, 5 xfailed — **38.94 s** |
| `ORPHEUS_B32=N1` | **66 failed**, 1101 passed — 138.74 s |

⚠ A full 10-leg sweep is **~20–25 min**: each mutated leg costs ~140 s because a
broken boundary makes downstream fixed-source solves iterate to their caps. And
it runs *all* of `tests/sn/operators` (1174 tests), so its reds are not scoped
to the RG gates. **Both facts argue for scoping the battery to the RG gates when
it is promoted.**

### ⭐⭐ Why this item is real work, not bookkeeping — the failure has already recurred
`[M]` `tests/sn/operators/__pycache__/` contains
`conftest_mutate_kernel.cpython-314-pytest-9.0.2.pyo` and
`conftest_mutate_pr...pyo`, but **the sources do not exist and were NEVER
tracked** (`git log --all -- 'tests/sn/operators/conftest_mutate*'` → empty).
Two hand-written mutation plugins were written into `tests/`, used, and lost,
leaving only bytecode.

⚠ **`tests/_mutation/` is NOT the home** — it holds **cosmic-ray** configs
(tool-driven, module-scoped AST mutation, scored as a percentage). And its own
README recipe step 5 is `git checkout -- <module>` — **exactly the practice
`b3_2_mutations.py` exists to avoid**. `[M]` **no hand-written battery anywhere
in the repo is wired as a real pytest test**; all four live in `scratch/`.

## B4–B7 — reconciliation report (deliverable, not design)

| phase | verdict | evidence |
|---|---|---|
| **B4** — *the realizer composes `R ∘ G` generically; the descriptor algebra gains `∘`* | 🔴 **LIVE, not dissolved** — the natural next substantive step | 8-arm `isinstance` ladder at `realizer.py:791,831,836,911,936,955,1049,1117`; `_composition.py` has **only** `LawScaled`+`LawSum`, no `∘`; `as_operator` still deliberately absent (`_factors.py:242,252`). ⛔ **One sub-claim now FALSE**: "`OperatorProduct` used by no boundary operator" — `realizer.py:755` is `emission @ contraction` (G6.3 step 3). Owns 2 live B2.2 defects: the **unscaled** specular corner swap (`boundary.py:1029-1037`) and `compute_keff`'s leakage predicate (`solver.py:1603`). ⚠ Its stated "bit-identical" gate is the **§15.12 blind shape** — `[M]` a shared body left a 60-passed equivalence suite green while the operator was wrong on 3 of 5 quadratures; the catchers are the independent hand-built anchors, which must not later be deleted as redundant. |
| **B5** — *white typed as rank-one, adjointable* | ⚠ **~2/3 DISSOLVED** | `AngularAverageOperator` **RETIRED** (`b4f0f5c9`) — **factored**, not re-typed, into `PartialCurrentOperator @ IsotropicEmissionOperator`. C-2's flip **DONE** at G6.3 step 3, whose own note says it *"absorbed the B5 step this comment was waiting on"*. `[M]` white realizes adjointable today. **Remainder: retire reflective's `albedo` parameter** (3 stale forward refs at `_factors.py:130,1197`, `boundary_conditions.rst:627`); the composite `B.H` gate is really a **B6/#189 dependency**. ⚠ Re-derive B5's scope from the tree — its section describes an operator that no longer exists. |
| **B6** — *every law reachable from a tag; rank-1 via the walker* | 🔴 **LIVE** | `[M]` `SNMesh.BOUNDARY_OPERATOR_REGISTRY == ['reflective','vacuum']` (2 of 7 laws); diffusion has 4. #189/#244 OPEN. **Standing user ruling: rank-N is WIRED, not retired — do not re-litigate.** |
| **B7** — *theory re-derivation + sibling realizers* | ⚠ **SCOPE-SPLIT** | #179/#180/#181 (MoC/MC/CP realizers) are **out of scope** under the priority law (CP/MoC/MC are later targeted campaigns). In-scope half is the theory-page re-derivation. #219 OPEN. |
| **#335** | **genuinely later** | its own body: *"today's meshes are structured…"* — a capability extension conditioned on unstructured meshes. B4's `R ∘ G` is its prerequisite algebra, not its delivery. |
| **#330 / #331** | #330 **partly landed, and B4 is a #330 increment**; #331 **separate** | G6.3 step 3 is tagged `(#330)`; `[M]` #331 is real — `SNBoundaryOperator` raises on an `AngularBoundaryDisplacement` composite. B4 must not worsen it. |

### Plan-vs-tree drifts to fix while here
- ⛔ **`SpatialWrap` is GONE** (`e13313a8`) — both plans still instruct in terms
  of `SpatialWrap.permutes_ordinates`. The predicate now lives on
  `PairedDeck`/`SelfPairedDeck` and is **derived from the motion**.
- ⛔ The cited tripwire `test_reemission_closure.py:1134-1141` no longer holds
  the claim it is cited for.
- ⛔ `tests/geometry/test_law_composition.py:64-69` says a gate is worthless
  *"until B3.4b lands"* — **B3.4b and B3.4c both landed 2026-08-01**, so the gate
  is now stateable and the prose is present-tense-false.
- ⛔ Both plans' cold-pickup anchors: the review's `:126-133` ("B2 measured, not
  started") is contradicted by its own `:1039` ("B2 — RESULT ✅"). **Late wins.**
  The review's `:6` still says `NEXT = #325, then B4`; **#325 CLOSED 2026-08-13**
  — the third time that anchor has lied forward, exactly as its own `:8-10`
  warns.

---

# Thread A — Q6

**User ruling (2026-08-14a):** *reference axis first, then the full four items.*
**User ruling (2026-08-14b), superseding the sequencing only:** *stop and
re-plan Q6 against the findings below, then resume from the corrected plan.*
The four items and their order are unchanged; what follows is the reconciliation
that ruling asked for.

---

## ⛔⛔ RE-PLAN 2026-08-14 — six premises refuted, all `[M]` at HEAD `66942848`

Everything in this block was measured after the plan was written. Per
`plan-authoring` §3 the original text stays below with its refutation beside
it; read this block FIRST, because it is the errata and it governs.

### R-1 ⛔ "the selector's V stage compares `deg(Q) ≥ d`" — there IS no such comparison
`registry.py:900-910`. Stage 2 is a parameter **inverter**, not a comparator:
`params = spec.degree_of_exactness_for(target_degree)`, and the ONLY V-stage
rejection is `params is None`. No `>=` on a degree exists anywhere in
`select_quadrature`. The measure — which carries the whole `ExactnessClaim` —
does not exist until `spec.build(params)` at `:912`, i.e. **after** V and
**before** stage 0 at `:915`.
⟹ a reference check is a **NEW STAGE in the `:912-915` gap**, not a conjunct
added to an existing predicate. The plan's phrasing ("the selector becomes
reference-aware") is an outcome and survives; its implied mechanism does not.

### R-2 ⛔ `OrthogonalSystem` CANNOT discriminate the case item 4 is dangerous for
This was the natural fix and it is refuted by measurement. `[M]`:

| | `LEGENDRE.gauss(4)` | `CHEBYSHEV_T.gauss(4)` |
|---|---|---|
| `system` | ALGEBRAIC | ALGEBRAIC ← **same** |
| `reference.support` | `[-1,1]` | `[-1,1]` ← **same** |
| `degree` | 7 | 7 ← **same** |
| `reference.name` | `legendre` | `chebyshev_t` ← the ONLY difference |
| `Σ w·x⁶` | `0.285714` ✅ | **`0.981748`** (exact `2/7`; err **0.696**) |

`exactness.py:23-35` names two failure modes. `OrthogonalSystem` catches the
**second** ("different SPACES, same number"). Item 4's hazard is the **first**
("different MEASURES, same number"). Only the reference itself separates them.
⟹ the new stage must compare the **reference**, and nothing weaker.

### R-3 ⛔⛔ "Lebesgue on `[-1,1]`" has THREE incomparable spellings — Cardinal Rule 2
`[M]` all three denote the same measure; all three compare `!=`:

| spelling | `.name` | `.support` |
|---|---|---|
| `LEGENDRE` | `legendre` | `[-1,1]` |
| `LEGENDRE.on(-1,1)` | `legendre_on[-1.0,1.0]` | `[-1.0,1.0]` |
| `UniformMeasure(...)` (built by `equispaced`, `measure.py:1516`) | `uniform([-1.0,1.0])` | `[-1.0,1.0]` |

`.on()` does **not** canonicalise its identity case, and the support strings
differ by float formatting alone. Meanwhile `exactness.py:161-164` promises
`name` is *"canonical mathematical identity — what equality compares."*
`[M]` that promise holds **within** each class (`GeneratingMeasure.__eq__`
compares `(name, support)`; `UniformMeasure.__eq__` compares
`(support, orthogonal_system)` and its `name` is a derived property that does
**not** participate) and is **silent across** them — cross-class dataclass
`__eq__` returns `False` before it looks at any field.
⟹ **This is a real architectural defect and it is NOT Q6's to fix.** It is
shared-numerics surgery (`exactness.py` + `generating_measure.py` +
`measure.py` + `test_exactness.py`'s 21 gates). File it; do not absorb it.
⚠ But Q6.0 must be **designed so the defect cannot bite it** — see the
corrected Q6.0 below.

### R-4 ⛔ Q6.0 is INERT if it lands alone — its only witness is item 4's own entry
`[M]` on today's 4-rule registry, domain (stage 0) already implies reference on
all 16 (geometry × rule) rows: every `S²` rule carries `UNIFORM_ON_SPHERE` and
the only interval rule is GL. So a reference stage added now can never reject
anything — a gate that cannot fail, wearing an authoritative name
(`coding-standards`, the tautological-gate hazard).
The witness that makes it bite is **Gauss-Chebyshev**, which is item 4's own
registration. ⟹ **Q6.0 and the Gauss-family half of Q6.4 are ONE step.** The
plan's step boundary cut across the thing that makes the step verifiable —
`plan-authoring` §6b, whose founding case was a signature's call sites; this is
the same defect with a *witness* instead of a call site.

### R-5 ⛔ `test_every_registry_rule_declares_a_symmetry_it_actually_has` has lost its example
`tests/numerics/test_symmetry.py:966-1010` (⚠ **not** in `test_registry.py`, as
the plan's §"three named gates" implies). Its docstring justifies being
*"deliberately WEAKER than equality with the computed maximal group"* by citing
*"one shipped rule genuinely is [true-but-not-maximal]: `gauss_legendre_on_mu`
declares `SO(2)`"*. Both halves are false at HEAD:
- `[M]` `gauss_legendre_on_mu(8).invariance_group` is **`Mirror('x')`**, not `SO(2)`;
- `[M]` walking `maximal_invariance_groups` over all four registry rules:
  **0 of 4** are true-but-not-maximal (`GL → sigma_x ∈ {sigma_x,sigma_y,sigma_z}`;
  Lebedev `Oh`; LS_N `Oh`; Product `D_6h` — each declares one of its own maximal groups).

⚠ **The gate should NOT be strengthened to equality.** A declaration is
*permitted* to be non-maximal, and a gate must assert the contract, not the
accident that no current rule exercises the slack (`vv-principles` #16). What
is broken is the **evidence in the docstring**, and the repair is to re-justify
the weakness from the contract instead of from a vanished example.

### R-6 ⛔ "each `QuadratureSpec` docstring narrates the trade-off" describes something that CANNOT exist
`registry.py:106-107` and `discrete_measures.rst:974-975`. `[M]` `QuadratureSpec`
is a dataclass, so all four *instances* share the class docstring
(`spec.__doc__ is type(spec).__doc__` for all four); the entries carry inline
`#` comments only. And the module is not `automodule`-rendered, so nothing
"becomes a Sphinx table row" either.
⭐ **This is Q6.2's real motivation**, and a far firmer one than the plan's
abstract "key, identity, label are three things": the registry *promised*
per-entry narration and the dataclass-instance design **cannot deliver it**.
The label is what discharges a promise the tree already made and broke.

### Also found, being repaired in parallel (not Q6's content)
The theory page's whole "Quadrature selection algorithm" section
(`discrete_measures.rst`, from ~`:682`) describes the **pre-2026-08-02** design:
it says "four-stage filter" (the code has five), its stage 1 is the retired
*declared-tag* form `G_geom ⊆ G_Q`, **stage 0 (domain) is absent entirely**,
and the labelled equation `quadrature-selection-criterion` carries the retired
3-conjunct body while the "where" list beneath it defines `𝒟_Q` and `Sym(Q)` —
symbols the equation does not contain. The same label is declared a second time
at `registry.py:61` with the *correct* 4-conjunct body, silently, because that
module is not rendered. → dispatched to `archivist`.

### R-7 ⛔ Item 3's Lebedev numbers are wrong, and the real gap is 12 orders wider
`[M]` scipy 1.17.1 in this venv serves Lebedev to order **131** (32 orders), not
59. `_LEBEDEV_ORDERS` (`registry.py:242-244`) covers **18 of 32**, topping at 47.
Its own comment says it was *"probed at registry-construction time on the SciPy
installed in this dev container"* — `[M]` that probe is `60f9fb29`, **2026-05-06**,
never re-run.
`[M]` the live consequence: `_lebedev_invert(53) is None` while
`Quadrature.lebedev(53)` serves 974 nodes summing to `4π`. **The selector refuses
a degree the tree can deliver.** That is a correctness defect, not a tidiness one.
⚠ **A blocker the plan does not carry:** `[M]` `_lebedev_node_count({'order': 53})`
raises `KeyError`. `_LEBEDEV_NODE_COUNTS` (`:424-428`) is a **second independent
transcription** of the same order set, so refreshing one without the other makes
the selector crash at stage 4. They must be single-sourced in one change.
`[M]` **3 tests** pin the stale copy, not the 1 the plan names:
`test_lebedev_invert_too_high_returns_none` (both assertions flip),
`test_no_rule_fits_raises_with_log` (its only-Lebedev registry stops raising),
and `test_lebedev_invert`'s parametrize (survives — all rows ≤ 47).
`[M]` **9 transcription sites tree-wide, 5 in `orpheus/`, already diverged 3 ways**
— incl. `angular_trace_space.py:82` claiming "Lebedev orders 3..53", which
`registry.py` asserts do not exist.

### R-8 ⛔ `ParameterSpec` is not "hardening a channel" — the field it targets is DEAD
`[M]` `select_quadrature` **never reads `spec.parameters`**. It calls
`degree_of_exactness_for(target)` and passes the returned dict straight into
`spec.build(...)` → `factory(**params)`. The field's only consumers are its own
docstring and a Sphinx table that (per R-6) does not render.
⟹ the honest framing is **"give a dead field its first job"**, which has a
different risk profile and a different done-when than "add validators to an
existing channel". And the two arms the plan wanted to reconcile do not need one:
Lebedev's fix is R-7 (derive the set), and level-symmetric's is *no change at all*
(R-9).

### R-9 ✅ The level-symmetric counter-ruling is VINDICATED — by a real event
`[M]` `_ls_sn_invert` discovers the bound by ATTEMPTING `level_symmetric_sn(n_min)`
(`:474-481`), with no literal. `59bb38a0` (#337) moved the frontier S12→S18 and
**did not touch `registry.py` at all**; the inverter tracked it correctly
(`[M]` `_ls_sn_invert(17) = {'sn_order': 18}`). ⟹ **do not tabulate this bound.**
⚠ But its *prose* became the exact failure it warned against — see R-10.

### R-10 ⛔ `registry.py`'s `_ls_sn_invert` docstring is stale THREE ways, and the irony is exact
`[M]` even N from 2 to 24: **S_18 is the last that builds; S_20 is the smallest
that raises.** So:
- `:445-451` — "no POSITIVE solution above `S_12`, `[M]` min weight `-0.027` at
  `S_14`" is **FALSE**. It is the pre-#337 convention seed's frontier;
  `rules_sphere.py:1072-1079` already carries the correction, dated 2026-08-13,
  **for its own file only**.
- `:436` — "conservative `deg = N - 1`" is **FALSE**: `[M]` the realized degrees
  are S2→3, S4→5, S6→7, S8→9, S10→11, S12→11, S14→15, S16→15, S18→17. It
  under-claims by 2 at S4/S6/S8/S14 and is right only at S2/S12/S16/S18 — so a
  spot-check on the wrong order confirms it.
- ⭐ **The irony, and it is the lesson:** the paragraph at `:453-459` warning that
  *"a literal `12` here would be a second copy of a frontier that lives in the
  solve, and the two would drift"* sits **eight lines below a `12` that has
  drifted**. The design's CODE half was vindicated; its PROSE half became the
  failure it predicted. Anyone citing this ruling must cite `:453-459` **and**
  flag `:445-451` as the counter-example above it.
- `[M]` a **third** site carries the same stale pair: `discrete_measures.rst:715-719`.

### R-11 ⛔ `registry.py:300-303` "so the selector stays cheap" — mechanism true, consequence dead
The node-counters really are pure. But `[M]` `select_quadrature` instantiates
**every** surviving candidate at `:912` (stages 0/1 need the nodes), and
`_ls_sn_invert` instantiates too. Measured wall-clock:

| geometry | deg 21 | deg 31 |
|---|---|---|
| slab | 0.787 s | **26.985 s** |
| sphere | — | **32.419 s** |
| cylinder | — | 29.040 s |
| cartesian2d | — | 28.987 s |

⭐ **A slab pays 27 s to construct a level-symmetric rule it can never select**,
because stage 2 (V) runs before stage 0 (domain) — and `registry.py:81-89` argues
that ordering is a consequence of the theorem, so it is not casually reversible.
⚠ The file **contradicts itself**: `:457-459` says *"Selection is not a hot path
(`[M]` the selector has no production consumers at all today), so the honest check
is affordable"* — which is TRUE — while `:300-303` claims cheapness. Keep
`:457-459`; retire the claim at `:300-303`.

### R-12 ⛔⛔ Item 4's `folded_product` half is BLOCKED — and not by anything `ParameterSpec` reaches
`[M]` `folded_product(4,8)` returns a measure with **`exactness=None`,
`degree_of_exactness=None`, `invariance_group=None`**, on support
**`'S^2/sigma_y'`**. It is the only shipped family carrying no exactness claim at
all. Registering it therefore fails **two** stages structurally:
- **stage 2** has no degree to invert;
- **stage 0** compares `measure.support` against `AngularSymmetry.support`, which
  `[M]` only ever yields `[-1,1]` or `S^2` (`registry.py:648-659`) — no geometry
  can equal a quotient support.

⟹ **`folded_product` registration is not a Q6 step.** It needs an exactness claim
for the quotient (a theorem, not a tag) and a domain notion that admits quotients.
File it.
⚠ Related, and it removes item 3's fourth-parameter-kind motivation: `[M]` `shift`
has **zero independently-selectable call sites** — `product` welds
`NODE_ALIGNED` (`rules_product.py:684`), `folded_product` welds `STAGGERED`
(`directional.py:746-751`), and `folded_product`'s docstring argues the welding is
*deliberate* (offset-selection and fold-intent must be co-selected or the
incoherent combination becomes spellable). There is nothing for a validator to
validate.

### R-13 ⛔ Item 2's precedent does not exist; the real one is one package over
- `[M]` **`SubgroupOfO3` is not a registering constructor** — no `create`, no
  `registry`, no `from_name`; its 7 named members are `ClassVar` singletons
  assigned at module scope (`symmetry.py:1384-1390`). And it is **not frozen and
  not a dataclass**: `[M]` `g = SubgroupOfO3.Mirror('z'); g._tag = 'MUTATED'`
  **succeeds**.
- `[M]` **"no precedent exists" is FALSE.** `orpheus/geometry/boundary/albedo.py:33-34`
  is `@dataclass(frozen=True) class AlbedoBoundary(BoundaryTraceLaw, key="albedo")`
  — a frozen dataclass with parameters, mounted on `RegistryMixin` with a key.
  That is exactly the shape the plan wanted, already shipped.
- `[M]` **`RegistryMixin` is not contradicted "by design."** Its claim says
  *"rather than maintaining a parallel **dict**"*; `quadrature_registry` is a
  **list**, and it matches `LOSS_REPRESENTATIONS`
  (`sn/loss_representation/__init__.py:2645`) — the tree's *other* established
  idiom, serving the *other* job. The split is by JOB: `RegistryMixin` =
  **name-addressed lookup**; ordered sequence = **rank-and-select**.
  `select_quadrature` is rank-and-select and has **no key parameter**.
  `[M]` `RegistryMixin.create()` has **zero production call sites** tree-wide; it
  functions as a completeness ledger for audits, not as dispatch.

### R-14 ⛔ Item 2's worked example is INVERTED
`[M]` `QuadratureSpec.name` is **already** a label, and correctly declared as one
(`:275-277` "Human-readable identifier"); all six production uses are f-string
interpolation and **zero** are lookups. So the missing concept is not the *label*
— it is a **key** (there is none) and a mathematical **identity** (there is none).
⭐ **The genuine evidence for the third concept is `SubgroupOfO3` itself**: it is
the one type that already split identity (`_tag`) from label (`.name`) — and it
had to invent a **third** string to key on, `repr()`. Its own comment
(`symmetry.py:454-459`) says the invariant is enforced by that comment rather than
by a type. That is the plan's thesis, stated by the code.
⚠ Real latent bug worth fixing regardless: two specs sharing a `name` would make
`dict(log.rejected)` silently drop a row — an idiom the theory page **teaches** at
`discrete_measures.rst:957`. `[M]` `quadrature_registry` is also the tree's only
**mutable list** registry; every sibling is a `tuple`.

---

## ⛔ Premise corrections *(as written before the re-plan; R-1…R-14 above govern)*
1. **`invariance_group` is NOT a field** — retired at `e7d44f3c`, absence
   gated. Its retirement **binds this design**: a rule's symmetry is
   *parameter-dependent* (`product_mu_phi`'s is `D_{n_φ h}`), so a parameter-free
   field cannot state it truthfully; the field that tried shipped `SO(2)`, false
   for every finite point set on `S²` (`registry.py:261-273`). ⟹ **nothing
   parameter-dependent may be stored on the spec.**
2. ⭐⭐ **Item 4 is UNSAFE as written.** `[M]` the most natural
   `GaussChebyshev1D` entry the current shape admits **passes all five selector
   stages**; list position alone picks the winner, and when it wins
   `Σw = 3.14159…` (not 2) and `Σ w·x⁶ = 0.9817` against exact `2/7 = 0.2857` —
   **off by 0.696**. Cause: `degree_of_exactness_for: Callable[[int], …]` takes a
   bare degree with **no notion of which measure it is against**. Latent only
   because `select_quadrature` has **zero production callers**.
3. **The twin-source premise holds, and the two arms want OPPOSITE things.**
   Lebedev's admissible set has **5 copies, already diverged** — `[M]` scipy
   serves to order **59**, `_LEBEDEV_ORDERS` stops at **47**, so
   `_lebedev_invert(53) is None` for a degree scipy can serve, and
   `test_lebedev_invert_too_high_returns_none` **pins the stale copy**. But
   level-symmetric carries a dated ruling *against* tabulating its bound
   (`registry.py:452-459`) — **vindicated** when the frontier moved S12→S18 at
   #337 and the inverter followed with no literal moved. ⚠ `_ls_sn_invert` costs
   `[M]` **26.9 s** to refuse at degree 31.

## ⭐ The spine: the family axis IS the exactness-reference axis

*A family is exactly "what the degree is measured against."* Golub–Welsch-on-`W`
has reference `W`; the circle rule `UNIFORM_ON_CIRCLE`; the tabulated sphere rule
`UNIFORM_ON_SPHERE`; a product **derives** from its factors. The tree already
types this one layer down and the registry reads none of it:
`ExactnessClaim.reference` (`exactness.py:148+`), `OrthogonalSystem` (`:116-147`
— `ALGEBRAIC`/`TRIGONOMETRIC`/`SPHERICAL_HARMONIC`), and
`GeneratingMeasure.orthogonal_system` (`generating_measure.py:260-281`) which
proves `ALGEBRAIC` **as a theorem of the construction, not a choice**.

⟹ **derive the family; do not mint a parallel taxonomy.** One addition closes the
T12b gap that stopped one layer short *and* makes item 4 safe.

## ✅ CORRECTED ORDERED STEPS — 2026-08-14, derived from R-1…R-14

⚠ The original six steps are kept verbatim below this block (§3). Where a step
survived it is named; where it dissolved, its residue is named. **Three of the
four chartered items changed shape**, so the list below is what the campaign
actually is.

**Goal of the whole thread (unchanged, and it is an OUTCOME):** *a quadrature's
degree is never readable without the measure it is a degree against, and the
registry cannot hand back a rule that is exact against a different measure than
the caller asked for.*

### Q6-A — the falsified-prose sweep ▸ do FIRST, no design risk
**Goal.** No live text in the quadrature layer states a frontier, a cost, or a
mechanism that the tree contradicts.
**Why first:** every item is a present-tense falsehood *today* (the standing rule:
a falsified doc is a BUG, fix on sight); none depends on any design decision; and
two of them (R-10's `S_12`, R-7's `47`) are numbers a later step would otherwise
copy forward. Highest value-to-risk ratio in the thread.
Sites, all `[M]`: `registry.py:445-451` (S_12/−0.027), `:436` (`deg = N-1`),
`:300-303` (cheapness — retire; `:457-459` already states the truth),
`:106-107` (per-entry docstrings that cannot exist),
`discrete_measures.rst:715-719` (same S_12 pair),
`test_registry.py:185` ("N=12 — the LAST order the family can serve", sitting
directly above the test that serves S14/S16),
`test_symmetry.py:976-983` (the vanished true-but-not-maximal example — re-justify
from the CONTRACT, do **not** strengthen the gate: `vv-principles` #16),
`numerics/registry.py:31-49` (`RegistryMixin`'s usage example).
▸ `discrete_measures.rst`'s selection section is dispatched to `archivist`
(incl. its own copy of the `S_12` pair at `:715-719`).
**Done when:** each site states what the tree does, with the superseded text kept
past-tense and dated; `sphinx-build -W` clean.

#### ✅ Q6-A LANDED (uncommitted) — plus two findings and one self-refutation

⛔ **`angular_trace_space.py:82` is NOT a falsehood — struck from this list.** The
recon flagged "Lebedev orders 3..53" as contradicting the registry's 47. It does
not: the sentence describes the coverage of an **empirical eps probe**
(`eps_probe.py`), not the selector's admissible set, and `[M]` `lebedev_sphere(53)`
really does build (974 nodes, `Σw = 4π`). Only `_lebedev_invert` refuses 53. ⟹ a
grep hit that shares a *number* with a real defect is not a site;
`coding-standards` already requires triaging concept-grep hits **by meaning**, and
this is that clause earning its place a second time.

⭐ **`RegistryMixin`'s example was worse than a dangling name.** `[M]`
`orpheus.numerics.operator.BoundaryOperator` **does exist** — as a
`_BlockRoleMeta` class with no `registry` and no `create` — so a reader who greps
the example's root *finds a class* and is confirmed in the wrong model. Same
mechanism as `plan-authoring` §1's one-letter-symbol clause (the check certifies
the error), here on a class name. Re-pointed at the real root `BoundaryTraceLaw`
with real members, so the example is now checkable rather than invented.

⭐⭐ **The `deg = N-1` staleness is not narration — it is a live 2× cost defect**
→ **Q6-E**, below. `[M]` feeding each realized degree back through `_ls_sn_invert`
returns an order **2 higher than the one that achieves it at 6 of 9 buildable
orders** (target degree 5 → `S_6`/48 nodes, when `S_4`/24 nodes has realized
degree 5). Correctness is unaffected — `N-1` is a true lower bound — but stage 4
ranks by node count, so level-symmetric is priced at up to 2× and can lose a
tie-break it should win.

`[M]` **164 passed** (`test_registry` + `test_symmetry` + `test_boundary_trace_law`),
**45 passed** on `test_registry` after the second `registry.py` pass, and
`sphinx-build -W --keep-going` **succeeded with zero warnings**
(`tools/check_docstring_xrefs.py` → `DEAD TARGETS: 0`).

⭐ **The archivist found an EIGHTH falsehood three screens ABOVE its brief** —
`discrete_measures.rst:555` opened *"Quadrature selection in ORPHEUS therefore
reduces to a containment check"* with the retired whole-group mapping (slab →
`SO(2)×σ_x`, sphere → `O(3)`) and closed by calling containment *"sufficient to
preserve every symmetry the geometry exhibits"*. Repairing only the named section
would have left the page citable for either sentence. ⟹ **when a section is
falsified by a design change, grep the page for the change's OTHER framings before
declaring the repair done** — a correction pass scoped to the section a reviewer
noticed inherits that reviewer's scope, not the defect's.

⭐ **And it handed back two survivors in MY file that my own sweep missed**, both
"four-stage": `registry.py:229` (the `See Also`, now a precise
`:ref:`quadrature-selection-algorithm``) and `:888` (`select_quadrature`'s own
docstring, *"The four-stage filter — G, V, structural, minimum-points"* — the
count AND the enumeration wrong, domain missing). I had grepped the *numbers*
(`S_12`, `47`) and the *concepts* I knew were stale; "four-stage" was a third
vocabulary I never searched. Plus `:1065`'s tie-break comment claiming registry
order "puts the most specialised rule first" — `[M]` nothing establishes that, and
no test exercises the tie-break at all (no shipped pair produces an equal node
count), so the policy is unstated AND unpinned. All three fixed.

### Q6-E — the inverter asks the family for the DEGREE, as it already does for the FRONTIER
**Goal.** The order the selector returns is the cheapest that meets the target.
⭐ The ruling is already in the file and applied **half-way**: `_ls_sn_invert`
discovers the *frontier* by attempting the construction (R-9, vindicated) and
still *tabulates* the **degree** as `N-1`. Reading `degree_of_exactness` off the
build it already performs, then stepping down while the target is still met,
applies one principle to both — and is why the `N-1` prose could go stale at all.
⚠ **Behaviour change with pinned tests**: 3 of 4 `test_ls_sn_invert` rows flip
(target 3 → `S_2`, target 4 → `S_4`, target 11 → `S_10`). Re-pose them as
"cheapest order meeting the target", not as "the `N-1` inversion".
⚠ Costs one extra build per step-down; R-11 already measures the selector at
27–32 s, so quantify before landing.

### Q6-B — Lebedev's admissible set is DERIVED, not transcribed ▸ a real refusal bug
**Goal.** The selector serves every Lebedev order the installed scipy serves.
**Why it is not a `ParameterSpec`:** R-8 — there is no channel to harden. The fix
is to compute the order set (and the node counts) from the source of truth instead
of freezing a 2026-05-06 probe.
⚠ **Single-source `_LEBEDEV_ORDERS` and `_LEBEDEV_NODE_COUNTS` in ONE change** or
the selector crashes at stage 4 (R-7's `KeyError`).
⚠ Then re-pose the 3 pinning tests (R-7). ⭐ Two of them exist to prove a
*refusal*; a refreshed ceiling must leave them able to prove one — pick a target
above 131 rather than deleting them, or the change trades a bug for lost coverage
(`coding-standards`, the single-sourcing demotion clause).
**Done when:** `[M]` `_lebedev_invert(53)` returns an order scipy serves; the two
tables cannot disagree by construction; `test_registry.py` green.
⛔ **Do NOT tabulate level-symmetric's bound** — R-9, vindicated by a real event.

#### ✅ Q6-B LANDED — `2cb3c960` (clean) + `37f3f880` (extend)

⭐ **The blocker DISSOLVED instead of needing a patch.** `expected_node_count`
was computing a number the selector already held: the measure is built at stage 2
(`:998`, because stages 0/1 ask its nodes) and the formula was consulted **45
lines later** at stage 4. `_LEBEDEV_NODE_COUNTS` existed *only* to serve that
formula — so retiring the twin removed one of the two Lebedev transcriptions, and
left exactly one copy to replace. "Clean before extending" was not a stylistic
ordering here; it is what made the fix a one-copy change.
`[M]` the formula and the measure agreed on **all 25** shipped configurations —
which is what makes a twin dangerous, not safe. The family that breaks them
already exists: `folded_product` quotients by a mirror, so `n_mu·n_phi`
over-counts it **2×**.

⭐ **Asking scipy is FREE, so there was no reason to tabulate at all.** `[M]` a
miss costs ~0 s (scipy raises *before* constructing); the largest hit — order 131,
5810 nodes — costs **1.85 ms**; the full sweep costs **16.6 ms**, `@cache`d, and
sweeping to 1001 measures identical to sweeping to 201. The refusal bug is gone:
`_lebedev_invert(48)` → order 53, `(100)` → 101, `(132)` → `None`.
⭐ `_LEBEDEV_SEARCH_CEILING` is a **search window, not a claim** — a claim can be
wrong, a window can only be too small, and the function **refuses** rather than
truncating when it binds.

⭐⭐ **The two instruments catch DISJOINT halves, `[M]` measured, and the
docstring says so** — so neither is retired later as redundant. Ceiling on a
*served* order (47) ⟹ the **guard** fires, the oracle gate is never reached.
Ceiling in a *gap* (49, 51) ⟹ the guard is **silent** (the last order found is 47,
not the ceiling) and only the **advertised-list oracle** catches the truncation.
The oracle compares two independent channels of scipy's own truth — discovery by
*constructing* vs scipy's *advertised* list in its `NotImplementedError` — which
is what lets the two refusal gates derive their thresholds from `_lebedev_orders`
without becoming self-referential: they check the *inverter* against the set, it
checks the *set* against scipy.

⚠ All three pinning tests **re-posed, none deleted** — two exist to prove a
REFUSAL, and both were pinning the defect (`48 is None`, `100 is None`). Their
thresholds are now derived, so they keep proving a refusal when scipy's table
grows. `[M]` **291 passed**; `npx pyright` 0 errors.

### Q6-C — the reference stage, landing WITH its witness ▸ the thread's core
**Goal.** A rule exact against the *wrong measure* cannot be returned. (The 0.696
defect made unspellable.)
**Shape** (R-1: a NEW stage in the `:912-915` gap, not a conjunct):
`AngularSymmetry` gains a **derived property** `reference` — `SO2 → LEGENDRE`,
`Trivial → UNIFORM_ON_SPHERE`, else `NotImplementedError` — mirroring the existing
`support` property one-for-one. ⭐ A *property*, never a field: `[M]` a new field
reddens 5 exact-equality constructor pins; a property reddens none.
Then the stage compares `measure.exactness.reference` against it.
⚠ **R-3 is the hazard, and the design must neutralise it rather than fix it.**
`==` on references is sound here only because the registry's reference vocabulary
happens to be `{LEGENDRE, UNIFORM_ON_SPHERE}`. ⟹ **close the vocabulary and gate
the closure**: one test asserting every registered rule's reference is one of the
canonical constants, so the day a rule arrives spelled the third way it reddens
**at registration** rather than silently mis-selecting. That converts an accident
into an enforced invariant.
⭐ **Lands WITH the Gauss-family registration** (R-4) — Gauss-Chebyshev is the
only witness that makes the stage bite; alone it is inert on all 16 rows.
⚠ Its rejection string must satisfy `test_log_rejected_list_carries_reasons`,
which requires every reason to name a stage from a fixed marker list.
**Done when:** `[M]` a `GaussChebyshev1D` entry is REJECTED for slab/sphere with a
stage-naming reason, GL still wins, and an in-process mutation removing the stage
turns that red green.

### Q6-D — the registry becomes a tuple; duplicate keys cannot ship ▸ small
**Goal.** Two rules cannot share a name, and the registry cannot be mutated.
`[M]` R-14: a shared `name` makes `dict(log.rejected)` silently drop a row — an
idiom `discrete_measures.rst:957` **teaches**. And `quadrature_registry` is the
tree's only mutable-**list** registry; every sibling is a `tuple`.
⛔ **Do NOT mount on `RegistryMixin`** (R-13): it is the *name-addressed-lookup*
idiom, `select_quadrature` is *rank-and-select* and has no key parameter, and
`create()` is production-dead tree-wide. `LOSS_REPRESENTATIONS` is the correct
sibling.

### ▸ FILED, NOT BUILT — each has a measured body, none belongs to Q6
1. **The three-spelling reference identity** (R-3) — `exactness.py` +
   `generating_measure.py` + `measure.py`; makes `exactness.py:163`'s promise true
   across realizations. Shared-numerics surgery; natural opener for the reshaping
   campaign.
2. **`folded_product` cannot be registered** (R-12) — needs a quotient exactness
   claim (a theorem) and a domain notion admitting quotients. Two structural gaps.
3. **`ParameterSpec`, re-scoped** (R-8/R-12) — the plan's version is refuted, but
   "should `spec.parameters` have a job at all, or be retired?" is a live question
   with a real answer either way.
4. **The selector costs 27–32 s at degree 31 on every geometry** (R-11) — benign
   today (zero production callers) and a blocker the moment one appears.

---

## ⛔ Ordered steps *(as written before the re-plan; the corrected list above governs)*
- **Q6.0 — the reference axis (PREREQUISITE for Q6.4).** ⚠ SURVIVES as **Q6-C**,
  fused with the Gauss registration. `QuadratureSpec` gains
  the exactness reference; `degree_of_exactness_for` and the selector become
  reference-aware so a degree is always *against a named measure*. Signature
  changes are free — zero production callers.
- **Q6.1 — `ParameterSpec` as a VALIDATOR channel**, modelled on
  `ClosureRecipe.validators: tuple[Callable, ...]`
  (`derivations/continuous/peierls_nystrom/geometry.py:6066-6082`). This is what
  reconciles the two arms: **Lebedev's domain DERIVED from scipy** (5 copies →
  1, closing a real refusal bug), **level-symmetric keeps its
  discovered-by-construction bound** plus a cheap obviously-out-of-range
  pre-filter — never a second copy of the frontier. ⚠ Must express the **fourth**
  parameter kind the task missed: `periodic_trapezoid(n, *, shift: Fraction)`,
  `NODE_ALIGNED`/`STAGGERED`, no default, quotient `ℚ/ℤ`
  (`rules_circle.py:178-271`) — and whether a parameter is **not independently
  selectable**, since `folded_product` *co-selects* it
  (`directional.py:706-712`). ⚠ `ClosureRecipe`'s own `name` is duplicated (dict
  key + `.name`) — cite it as the **negative** example.
- **Q6.2 — `label` is a THIRD concept.** `name` is deliberately a *mathematical
  identity*: `_jacobi_name(a,b)` canonicalises so `jacobi(0,0) == LEGENDRE`
  (`generating_measure.py:570-585`). Key, identity, label are three things.
- **Q6.3 — the registering constructor.** ⚠ **No precedent exists.** Follow the
  house idiom for a closed sum type *with parameterized members*:
  **`SubgroupOfO3`** (`numerics/symmetry.py:346`) — private
  `_NamedSubgroup(Enum)`, separate frozen dataclasses for parameterized members,
  a public frozen wrapper giving a uniform surface
  (`SubgroupOfO3.OctahedralOh` beside `SubgroupOfO3.Dnh(6)`). House idiom is
  `Enum`, **not** `StrEnum` (1 use) and **not** `Literal`. ⚠ Reconcile
  `RegistryMixin`'s claim that *"every extensible family should mount on it
  rather than maintaining a parallel dict"* (`numerics/registry.py:19-22`), which
  `quadrature_registry` contradicts by design.
- **Q6.4 — register what is missing**, now safe: the Gauss families, plus
  ⭐ **`Quadrature.folded_product`** — a **fifth shipped angular family with no
  registry entry** the task never mentions (the σ_y-quotient the cylindrical fold
  uses).
- **Q6.5 — the falsified prose** (bugs): `registry.py:445-465` states the LS
  frontier as `S_12` (tree does `S_18`); `test_symmetry.py:978-983` justifies
  itself with the retired `SO(2)` tag; `discrete_measures.rst` claims *"each
  `QuadratureSpec` docstring narrates the trade-off"* — false twice;
  `registry.py:300-303`'s cheapness claim; `test_registry.py:185`'s comment.

⚠ **The quadrature package is deliberately NOT `automodule`-rendered**
(`docs/api/discrete_ordinates.rst:76-78` — to avoid duplicate-label collisions
with the theory pages), so a `:class:` xref to a new type **silently renders as
plain text**. The remedy is *not* adding an automodule.

---

# Verification

**Baselines to establish before touching anything** (`[M]` where already known):
- `tests/numerics/test_registry.py` + the symmetry row → **46 passed, 2.40 s**.
- `tests/sn/operators` CONTROL leg → **1167 passed, 1 skipped, 5 xfailed, 38.94 s**.
- `tests/numerics/test_incoming_ordinate_mask_tensor.py` → 13 tests.

**Per step:**
- **B3.3**: the 3 present-tense-false sites fixed first (independent bugs); then
  delete + rewire; `tests/numerics` + `tests/sn/operators` + `tests/geometry`
  green; grep proves no live reference survives while all **16 past-tense sites
  remain**; the `implements` edge on eq. 10 deliberately re-pointed or retired.
  ⚠ **Sphinx cannot gate this** — grep is the instrument.
- **B3.5**: RG-3b's 4 xfail reasons name B4; the battery's CONTROL leg green
  from its new home and **at least one mutated leg still reddens its RG gate**
  (a promotion that loses the teeth is worse than none).
- **Q6**: the three named gates re-posed, not deleted —
  `test_registry_population` (the independent literal pin),
  `test_registry_specs_are_well_formed`,
  `test_every_registry_rule_declares_a_symmetry_it_actually_has` (RED three ways
  on registration: Gauss families are *deliberately* `invariance_group=None`;
  `equispaced` is degree 1 so cannot reach degree 5; the closing keyset). Plus a
  new gate that **`GaussChebyshev1D` can no longer be returned for an unweighted
  query** — the 0.696 defect made unspellable.
- **Wide gate before merge**: `python -O -m pytest -m "not slow"`, **SERIAL**,
  budget **≥90 min**, via `Popen(start_new_session=True)` + a persistent
  `Monitor` writing to a LOG (a harness background task dies at ~30 min).
- **Sphinx** `-W` clean after every doc touch.

**Sequencing.** Track B first (a half-done retirement is the state that breeds
"which path is canonical?"), Q6 second. Within B: the 3 doc bugs → B3.3 →
B3.5's two remainders → the B4–B7 report.

# The streaming path's objects say what they are, live where they belong, and carry only what they use

> ## ▶ RESUME STATE — ⏸ COMPACTION POINT, rewritten 2026-08-28 after P4.1a/b/c
>
> **This file is the resume surface. Trust it over any summary.**
> ⛔ It replaces a post-P3 header. A summary is quoting dead text if it says
> any of: "▶ NEXT = P2" · P2 moves α · P4 has a four-way fork about where the α
> type lives · `face_areas`/`delta_A` are fields · `streaming_terms` has three
> chart arms · P4.3 comes before P4.4 · `AngularMeasure` "dies" at P4.2.
> Every one of those was true once and is refuted below.
>
> ### Where the tree is — `[M]` reconcile with git, never with this line
>
> `main` == `origin/main`, **no open branch for this campaign**, working tree
> clean. ✅ Phase B, P1, P2, P3a/b and **P4.1a + P4.1b + P4.1c** are all merged
> and pushed. Reconcile with
> `git merge-base --is-ancestor 7c5a8fb3 main`, never with this sentence.
>
> ⚠ This repo has **no CI** (`.github/workflows/` absent). "main is always
> green" is enforced entirely by the local gate; an empty `gh run list` is not
> a pass.
>
> ⚠ **History was REWRITTEN 2026-08-27 and force-pushed.** A `git add -A` swept
> 212 untracked `scratch/` files (74 670 lines) into the P4.1a commits; the user
> ruled a rewrite, so `5ff113c5..main` was replayed without them and **every
> hash in that range changed**. `[M]` verified before the push: 16 subjects
> identical, every non-`scratch` path bit-identical to the commit the gate ran
> on (so the gate carried over unre-run), tracked `scratch/` back to its
> `5ff113c5` state, and all 745 working files restored to disk byte-for-byte
> from a `cp -a` backup — **`filter-branch`'s checkout DELETES files it
> un-tracks.** ⟹ if a summary or memory quotes a hash in `5ff113c5..main` that
> `git cat-file -e` cannot resolve, it predates the rewrite; re-find it by
> commit SUBJECT. `[M]` 2026-08-28: all **20** hashes this plan cites resolve
> and are ancestors of `main`. The rule landed in `process-discipline.md`
> ("never `git add -A` here").
>
> ### What has LANDED
>
> | commit | phase | item |
> |---|---|---|
> | `cc01dd27` | Phase B | the angular closure takes its two tensor factors |
> | …`1c93b14f` | P1 | nothing in the streaming path is held that nothing reads (7 commits) |
> | `dcd6a9f6` | **P2** | `sn/sweep/pole_angular_closure.py` → `sn/angular/closure.py`, byte-identical |
> | `75571c4d` | **P3a** | `redistribution_gram` → `redistribution_pairing` |
> | `4a3f2390` | **P3b** | `GeometryCoefficients` → `StreamingCoefficientCache` |
> | `8ffddfb9` | — | the post-P3 compaction point |
> | `af801d76` | — | the α-defect diagnostic measured the normalization, not the defect (+ its 8-case gate) |
> | `d48f4bf4` | **P4.1a** | the chart lives on the mesh, not on a copy beside it |
> | `7f08d1d7` | **P4.1b** | the slab is the sphere's zero-curvature case — 3 arms → 1 body |
> | `7c5a8fb3` | **P4.1c** | the `SNMesh` deprecation shims retire |
> | `16501ca0` | **P4.4** | the connection coefficients come home — 4 symbols to `sn/mesh/`, tests follow |
> | `5940deba` | **P4.2** | the angular factor comes home — 6 symbols to `sn/angular/`, and the L0 ladder takes α as a keyword |
>
> `[M]` exit gates, full fast set, 13 trees all `rc=0` each time:
> **P2 9815/0** (+1 vs P1 — the new `sn/angular/__init__.py` adds one case to
> the import-linter's `rglob` parametrize), **P3 9815/0** (delta 0 on all four
> axes), **P4.1a 9823/0** (+8 = the `test_alpha_defect_normalization` gate that
> rode in on the same landing, `[M]` collects exactly 8; delta 0 elsewhere),
> **P4.1b+c 9823/0 — delta 0 on every tree AND every axis**, and
> ⭐ **P4.4 9824/0 — `[M]` 13 trees `rc=0`, 64 min, and the +1 was PREDICTED
> before the run.** Delta 0 on skipped/deselected/xfailed, and localized to
> exactly three trees: `geometry` **−45**, `sn` **+45** — the migrated tests,
> conserved at TREE granularity, which is an independent check of the
> 18+45=63 collection count — and `root+harness` **+1**, the import-linter
> `rglob`ping one new module. Nine other trees +0. pyright 0, `sphinx -W` 0
> and `dead_references` 0 throughout.
>
> ### ▶ NEXT — **P4.3.** (P4.4 ✅ `16501ca0`, P4.2 ✅ `5940deba`.) Read §4ter first.
>
> ⛔⛔ **The order was corrected TWICE. `P4.2 → P4.4` (written 2026-08-28
> in this very header) is itself REFUTED — P4.2 cannot precede P4.4 either.**
> `[M]` 2026-08-28 by AST: the three `*_streaming` factories **stay in
> `geometry/`** and CALL `angular_redistribution` at **runtime**
> (`reduced_operator.py:1198/:1262/:1360`), and `ReducedStreamingOperator.angular`
> is annotated `AngularRedistribution` (`:628`) — so moving the α cluster to
> `sn/angular/` while they stay creates `geometry → sn`, forbidden with no
> tolerance (`geometry` is INPUT; the `TYPE_CHECKING` carve is `L1|L2` only).
> ⭐ **Injected and run** (the §6d discipline, not reasoning): the gate goes red
> **two ways** — `test_no_forbidden_imports[…reduced_operator.py]` *and* a hard
> circular import killing `import orpheus.geometry`, `orpheus.numerics` and
> `orpheus.sn.solver` in a fresh interpreter. Same shape as P2's α half.
> `[M]` the reverse was checked: **P4.4 does not depend on P4.2** — after P4.4
> the factories sit in `sn/` and name `angular_redistribution` / `AngularMeasure`
> / `StreamingTerms` back in `geometry/`, which is `sn → geometry`, **legal**;
> and the geometry residual (`AngularMeasure`, `StreamingTerms`, the 5 α
> symbols) names only `CoordSystem`, `Mesh1D`, numpy — no L2/L3 at all.
> ⟹ **reorder, do not fuse**, exactly as ruled for P4.3/P4.4 one paragraph down.
>
> ⭐ **The lesson, and it generalises past this campaign: a re-home's forbidden
> edge can be INVISIBLE TO AN IMPORT GREP, because the edge does not exist yet
> — the split CREATES it.** This header's own P4.2 entry ran a thorough §6d
> check and found the `derivations → sn` hazard, because that one is a literal
> `import` line. It missed `geometry → sn`, which is four *intra-file* call
> sites — precisely where an import-based audit does not look. ⟹ **when the
> mover shares a FILE with its callers, the file is the caller: enumerate
> intra-file references by AST, not only imports.** (`plan-authoring` §6d.)
>
> ⛔ **The earlier correction, which still stands.** The plan's own sequence was 4.2 → 4.3 →
> 4.4, and `[M]` P4.3-before-P4.4 creates `geometry → transport` (**0** today
> against **16** the other way) because `StreamingTerms`' only constructor
> stayed behind in `geometry/` at a **runtime** call site
> (`geometry/reduced_operator.py:856`, pre-P4.4) — a declared `FORBIDDEN_EDGES`
> violation with a red waiting for it. ✅ **REMEDIED by P4.4 `16501ca0`**: `[M]`
> that constructor is now `orpheus/sn/mesh/reduced_operator.py:531`, so the
> ordering constraint is DISCHARGED and P4.3 is unblocked. `[M]` the reverse dependency was then checked too and is
> clean (`sn → geometry` is legal, `geometry` being INPUT), so the two
> **reorder rather than fuse**. Full record, both directions, in **P4.3's ⛔
> block** — including the two scope changes it forces on P4.4 (the module
> DELETION moves to P4.3; `geometry/__init__.py`'s re-export must be **deleted,
> not re-pointed**, with 11 consumer files).
>
> **P4.2 — the angular factor comes home.** `[M]` re-verified 2026-08-27/28:
> - the movers are **5 symbols and only 3 contain "alpha"** — `alpha_dome`,
>   `_assert_alpha_dome_closes`, `_ALPHA_CLOSURE_ATOL`, **`AngularRedistribution`**
>   and **`angular_redistribution`**. A grep-driven reading finds 3 of 5.
> - `[M]` the destination `orpheus/sn/angular/redistribution.py` **does not
>   exist yet** (existence-checked); `sn/angular/closure.py` does.
> - ⛔ the hazard is real and enforced: `derivations/discrete/sn/angular_differencing.py:163-164`
>   imports `alpha_dome` at module scope, and `[M]` `FORBIDDEN_EDGES` has 15
>   keys including `"derivations": L2|L3`, with the `TYPE_CHECKING` tolerance
>   scoped to `L1|L2` only — so **L0 is not covered even for typing.** The fix
>   is that file's own ruled precedent: **the L0 ladder accepts α as a
>   keyword**; its delegating local `alpha_dome` retires.
> - ⛔ **there is a WRONG shortcut, and it is one line.** `[M]` the linter ships
>   a `WHITELIST` whose existing entries are all `derivations → L2/L3`
>   exemptions, so adding a row lands the move green. **Do not** — a whitelist
>   records a violation where the keyword fix removes one.
> - ⚠ `AngularMeasure` does **not** move here and does **not** "die" here (see
>   the corrections below).
>
> ### Corrections that supersede older text in this file
>
> | ⛔ the older text says | ✅ the tree says, `[M]` |
> |---|---|
> | P4 has a **four-way fork** about where the α type lives so `geometry` may name it | superseded — the module holds no geometry and **dissolves**, so the edge never arises. §4ter. |
> | the spatial-field retirement is **bit-identical** | ✅ **REMEDIED by P4.1b.** It was true for `coord` (3/3 charts) and false for `face_areas` (2/3). Both are now **derived accessors, not fields**, so the claim has no subject left. |
> | `streaming_terms` dispatches on three chart arms | ✅ **REMEDIED by P4.1b** — one body plus a 2-way resolution of what `direction_idx` MEANS (global ordinate vs within-level index). That surviving `if` is P4.7's, not a chart dispatch. |
> | P4.2 retires `AngularMeasure` — "the 3 factories read 0 of its 6 members" | ⛔ **SUSPENDED into CS5.** That argued CONSUMERS; the open question is TYPE DESIGN (should `DiscreteMeasure` branch into Spatial/Angular/Energy?). `[M]` it has **5 use sites, all inside `reduced_operator.py`**, zero elsewhere ⟹ it travels with the factories at P4.4 and strands nothing. CS5 rules whether it should EXIST; P4.4 rules where it LIVES. |
> | `AngularMeasure` "travels with the factories at P4.4" | ⛔ **REFUTED 2026-08-28.** The measurement (5 in-module use sites) was right; the inference was not. `[M]` `AngularRedistribution` NAMES `AngularMeasure` (`:1082`), so 2 of those 5 sites are in the α cluster. Sending it to `sn/` with the factories while α stays would be `geometry → sn` — forbidden. ⟹ it **stays in `geometry/` across P4.4** (factories read it as `sn → geometry`, legal) and **travels with α at P4.2**. Lands in `sn/` either way, so CS5 is not prejudged. |
| `delta_A` moves to `sn/` as a field | ✅ **RULED 2026-08-27 (user): it DISSOLVES into R.** ΔA is R's rank-1 realization, not a second object. P4.1b is its first step — `delta_A` is already a derivation over `mesh.areas`, so what P4.4 moves is a derivation, not an array. |
> | `SNMesh.face_areas` / `.delta_A` are deprecated read-throughs | ✅ **REMEDIED by P4.1c** — retired. `[M]` 11 readers, 0 in `orpheus/`; every consumer was a test and 4 were the tests verifying the shims. |
>
> ### Still BLOCKED, with a named blocker (not a defer)
>
> Two of §4ter's forks — **the producer's home** and the **`ChartConnection`
> name** — resolve as consequences of **CS5**
> (`space_and_kernel_binding_campaign.md` §5.5, *an axis can name the generator
> that made it*), whose NODAL half is the prerequisite. Everything else is
> unblocked: **P4.2, P4.4, P4.3, P4.5, P4.6, P4.7 and P4b.**
> ⚠ **P3 is NOT finished** — `ReducedStreamingOperator` → `ChartConnection`
> (P3c) is deliberately sequenced **after P4**, so the name describes what
> remains once `streaming_terms` is out of it.
>
> ### ⭐ The rules P1–P4.1 paid for — read before the next rename or re-home
>
> - `plan-authoring` **§6d** (P2): **a RE-HOME must check the import edge in
>   BOTH directions — and look for an import-linter FIRST.** The done-when that
>   failed was *true*: it asked whether the MOVER depends on its old home, when
>   the question is whether the OLD HOME depends on the mover.
> - `plan-authoring` **§6d, 2026-08-27** (P4.3): **a plan row written in a
>   RULE'S OWN VOCABULARY reads as that rule having been applied.** P4.3 said
>   *"no new edge **either way**"* — §6d's phrase, verbatim, with no `[M]`
>   beside it. The check had not run. ⟹ cite the rule's MEASUREMENT, not its
>   phrasing; rule-shaped language reads as compliance and never as doubt.
> - `plan-authoring` **§2, 2026-08-27**: **parse structured source with AST,
>   never a line window or a regex, when the answer is a MEMBERSHIP question.**
>   Two filters in a row said `derivations` was not a `FORBIDDEN_EDGES` key — a
>   `sed` window one line short, then a regex terminated by a `}` nested inside
>   `L3_PACKAGES - {"sn"}`. Both wrong, both toward "no constraint". A wrong
>   filter is not a coin flip when you are filtering for a **prohibition**.
> - `plan-authoring` **§1** (P3b): **a proposed NAME's ABSENCE can be evidence
>   AGAINST it.** `SweepCoefficientCache` was free in code *and* git history
>   **because it was rejected** (`lessons.md` L15). Grep the PROSE corpus.
> - `coding-standards`, **2026-08-27** (P4.1a): **retiring a duplicate promotes
>   whatever KEPT THE COPIES EQUAL from redundant to load-bearing** — usually a
>   production guard with no test. `[M]` the three factory chart guards had
>   **zero witnesses** tree-wide; the tests that asserted the retired copy were
>   tautologies and became the guards' witnesses.
> - `process-discipline`, **2026-08-27**: **never `git add -A` here**, and a
>   history rewrite DELETES what it un-tracks — `cp -a` first.
> - `.claude/lessons.md` **L61**: an unvalidated filter and a clean tree print
>   the same thing; ⛔ `grep` here is `ugrep`. Run completeness checks in
>   Python **with a positive control**.
>
> ### ⭐⭐ The finding that outlives every phase
>
> **A lesson's own worked example is not exempt from the lesson.** `lessons.md`
> L15 is *"a cache shape that mixes immutability strata hurts twice"* and holds
> `GeometryCoefficients` up as the RIGHT shape — while `[M]` it mixes **three**
> (7 fields bit-identical across 3 meshes, 4 mesh-bound, 2 traversal-only).
> Nothing prompts that re-check, because the exhibit is what taught you the
> rule. → §4bis, and **P4b**.
>
> ### The rulings that bind (✅ user; full text + rejected alternatives at §9)
>
> `ChartConnection` · `sn/angular/` · the chart sweep SPLITS · P4 mints the
> pairing, the moment measure is **P7** · the cache is
> **`StreamingCoefficientCache`** · `delta_A` **dissolves into R** · P4.1b
> **collapses the arm**.
>
> ⭐⭐ **The governing ruling:** *"scheduled to the end of this plan, **not
> defer and forget**."* Deferred work earns a PHASE NUMBER or a named external
> blocker. P4b, P6 and P7 exist because of it.
>
> ### ⚠ `scratch/` is UNTRACKED — and a `git clean` destroys it
>
> `[M]` **745 files on disk, 142 tracked.** It holds the four cited audit
> memos, every gate log (`_p2`/`_p3`/`_p41a`/`_p41bc_fast_gate.log`), the gate
> drivers and ~190 probes. ⛔ Do **not** fix this with `git add -A` — that is
> what caused the 2026-08-27 rewrite. If any of it should be tracked, stage it
> by path, deliberately, in its own commit.
>
> ### Running the gate
>
> `.venv/bin/python -O -m pytest -p no:randomly -m "not slow" -q` — SERIAL,
> `[M]` **~62 min** wall (P4.1a 20:52→21:56, P4.1b+c 22:33→23:34). Two trees
> are the whole cost: `[M]` **derivations ~37 min, sn ~15 min**; the other
> eleven total ~9 min. Use the per-tree driver (copy
> `scratch/_p41bc_gate_driver.sh`, change its `LOG=` line) launched DETACHED
> via `Popen(start_new_session=True)` — harness background tasks die at
> ~30–90 min. ⛔ Never pipe it through `tail`; watch the log with a `Monitor`.
> ⚠ `-m "not slow"` is part of the number: it deselects 227, including the
> known-inherited `cyl_2g_3reg_folded_4x8_dd_n40` red (#404, dup of #397,
> `[M]` bit-identical pre-branch).
>
> **The acceptance baseline for the next step is `9825 passed / 0 failed,
> 22 skipped / 227 deselected / 70 xfailed`** — state the delta against it and
> account for every unit. ⚠ Each phase that adds a MODULE under `orpheus/` adds
> **+1** collected test, because the import-linter `rglob`s every one. Measured
> three times now: P2 `sn/angular/__init__.py` +1, P4.4
> `sn/mesh/reduced_operator.py` +1 (9823→9824, layer gate 342→343), P4.2
> `sn/angular/redistribution.py` +1 (layer gate 343→**344**, and the full gate
> confirmed **9825/0, 13 trees `rc=0`, 62 min** — delta 0 on the other three
> axes, localized to `geometry` −18 / `sn` +18 / `root+harness` +1, everything
> else +0; ⭐ `derivations` **1637, Δ 0**, so retiring the L0 wrapper and
> re-signing three graders cost that tree nothing). ⟹ **P4.3 DELETES
> a module and moves `StreamingTerms` into an existing one, so expect −1 →
> `9824`** — the first phase of this arc where the count goes DOWN, and a
> reading of `9825` there means the file did not actually die.

**Status.** ⛔ ~~Proposed 2026-08-26, un-ruled.~~ **IN EXECUTION** — Phase B,
P1, P2 and P3a/P3b have landed (see the table above); all four §9 forks and
three further design questions are ✅ RULED. Successor work to the un-weld arc's
Phase B (`27576937`, `6859ca05`); governed by
`.claude/plans/posing_filtration_charter.md` (the posing filtration, R1–R23)
and `.claude/rules/plan-authoring.md`.

**What this rests on** — three measured artefacts, all self-contained. Read
them before designing to any phase here; this plan CITES their numbers rather
than copying them (§9 of plan-authoring):

| memo | what it settles |
|---|---|
| `scratch/ld_curvilinear_shape_derivation.md` | the 1-D spherical LD moment balances; the three matrices `M` / `G` / `R` |
| `scratch/tau_under_ld_dip_analysis.md` | the tensor factorization; τ's arity theorem; the seed cone risk |
| `scratch/adjoint_gram_ownership_audit.md` | what shapes the adjoint; what `R` IS; the ownership table |
| `scratch/specialization_audit.md` | coordinate/scheme specialization; the weld list; both smuggling directions |

---

## 0. The three findings that reorganize everything

Everything below is a consequence of these. If one is refuted, re-derive the
phase that rests on it rather than patching the step.

**F1 — the redistribution operator FACTORS, and the factorization is the
ownership map.** `[M]` `𝓡 = R_spatial ⊗ A_angular(τ, α, w)` (tau memo §3).
`R`'s **measure comes from the chart**; its **basis comes from the scheme**;
its **angular partner is the measure's weights**. No single existing object
owns it, which is why it has been homeless and kept getting fused onto
whichever neighbour was convenient.

**F2 — `R` is the Galerkin matrix of `∇·ê_r`, and `R₀₀ = ΔA` is the
divergence theorem.** `[M]` adjoint audit A3:
`R_kj = ∫ b_k b_j (∇·ê_r) dV`, `∇·ê_r = (d−1)/r`, `∫∇·ê_r dV = ∮ê_r·n̂ dA`.
⭐ The consequence that matters for naming and homes: the *"geometry-dependent
factor-of-2 in the α normalization"* this tree has carried as a per-chart quirk
is an **artefact of spelling the measure `r^(d−2)`**. Written as `∇·ê_r` the
constant is `(d−1)` and **nothing per-chart remains** — so a great deal of
apparent sphere-vs-cylinder specialization is a spelling, not a fact.

**F3 — the closure STRADDLES TWO CACHE STRATA, and that is the weld nobody
named.** `[R]` from the audits, newly assembled here: the closure's `τ`, `α`,
`c_in`, `c_out` and the angular partition are **mesh-free** (quadrature ×
chart only — `[M]` bit-identical across 5 meshes incl. a graded one, Phase B's
probe), while its `_dAw_per_level` is **mesh-dependent** (`ΔA ⊗ 1/w`). Today
one object holds both, so a re-mesh rebuilds algebra that cannot have changed.
⟹ **§5c's ruled function/evaluated-table split applies to the ANGULAR closure
verbatim.** The charter already ruled it for the SPATIAL closure; the angular
one has the identical structure and the mesh-binding hid it.

---

## 1. Renames — every one of these names is currently false

⚠ Renames are `git mv` + import sweep + the THREE-search retirement audit
(code, tests, **docs** — a Python-domain xref renders as plain text with no
warning at any severity). Each row's evidence is the reason, not decoration.

| current | proposed | why — with evidence |
|---|---|---|
| `ReducedStreamingOperator` | ✅ **`ChartConnection`** *(ruled §9.1; `StreamingCoefficients` rejected — it lends `L`'s own word to a non-operator, reproducing the very defect being retired)* | `[M]` **0 of 13** operator-surface members (`apply`/`domain`/`codomain`/`H`/`inverse`/`solve`/`apply_transpose`/`__matmul__`/`__add__`/`is_adjointable`/`block_role`/`system_role`); `SweepOperator` has 12 of 13. It maps nothing to anything. ⛔ And the name is **TAKEN** — the real reduced streaming operator is `L`, in the SN algebra, with genuine spaces. "Operator" is this codebase's most load-bearing word (the S4 amendment: *an operator is not an operator without its two spaces*); a struct wearing it teaches every reader the wrong thing. |
| `redistribution_gram` | `redistribution_pairing` | `[M]` its own `(n_mom, n_thread)` axes admit the **rectangular** ONETRAN case (Hill 1975 Eq. 32, the angular index closed on the cell average only). A rectangular object is a **pairing, not a Gram** — the word over-claims off the diagonal of its own design space. ⚠ **This name is MINE, from `6859ca05`.** |
| `GeometryCoefficients` | ✅ **`StreamingCoefficientCache`** — ⛔ **BOTH plan candidates were REFUTED on measurement, 2026-08-26; see §4bis** | `[M]` **0 of ~~15~~ 13 fields** (⚠ the 15 is a **pre-P1** count — P1 @ `ebe5d22f` retired `mu_start`; §3's REMEDIED-fact clause) are un-permuted chart data: 4 are Morel–Montry constants, 3 angular, 2 pure traversal, and the closest 3 are **chain-ordered** (i.e. already traversal artefacts). The name promises geometry and delivers a permuted sweep cache. |
| module `sn/sweep/pole_angular_closure.py` | `…/angular/closure.py` — ⚠ **executed in P2, not P3** (the re-home and the module rename are one `git mv`) | Two lies in one path. **"pole"** names the special case (the pole cell) for a family that closes the *whole* angular axis — `IdentityAngularClosure` never sees a pole. **"sweep"** is traversal (see §2). R15 already ruled the family `AngularClosure`. |
| `SNMesh.pole_angular_closure` (the ATTRIBUTE + the constructor kwarg) | `angular_closure` | ⭐ **Found during P2, 2026-08-26 — the module carried TWO occupants of the name and only one moved.** P2 renamed the module `pole_angular_closure` → `sn.angular.closure` on the grounds that *"pole"* names a special case for a family closing the whole angular axis (`IdentityAngularClosure` never sees a pole). The identical lie sits on the mesh attribute, and it is the SURVIVING spelling: `[M]` **62** of the tree's `pole_angular_closure` hits after the move are this attribute, across `sn/mesh`, `sn/loss_representation`, `sn/operators`, `sn/sweep/cache`, `transport/`, and 10 test modules. ⚠ It is also why P2's import sweep had to match QUALIFIED forms only — a bare-token rewrite would have corrupted every one of the 62. Left deliberately: renaming a public constructor kwarg is not bit-identical. |
| `SNMesh.curvature: str \| None` | **retire** — use `coord` | `[M]` a **stringly-typed duplicate** of the `CoordSystem` enum (`"cylindrical"` / `"spherical"` / `None`), with `is_cartesian` defined as `curvature is None`. Read at two matvec entries through a **defaulted `getattr`** (`vv-principles` #28's temporal twin). ⚠ And `augmented_mesh.py:124` still documents it as storing `α_n/r_i` — a numeric quantity it has never held. |
| the chart's four spellings | one | `[M]` `coord` 228 hits / `curvature` 68 / `geometry_kind` 105 across `orpheus/`, with 8 / 3 / several case-spellings respectively. One concept, four vocabularies — a grep for any one returns a confident partial answer. |

⚠ **Do NOT rename `ReducedStreamingOperator` in isolation** — §3 shows it is a
bundle of three things at three stages, and a rename would name a bundle that
is about to stop existing. Sequenced accordingly in §7.

---

## 2. Re-homes — what the filtration says, versus where things sit

**The governing test**, from the charter's own ontology: an object's home is the
stage at which **all** its inputs exist. A module that must invent a Protocol to
avoid importing one of its own inputs is in the wrong stage.

| object | today | belongs | evidence |
|---|---|---|---|
| `AngularRedistribution`, `angular_redistribution()`, `alpha_dome`, `_assert_alpha_dome_closes` (+ `_ALPHA_CLOSURE_ATOL`, which no symbol grep for the four returns) | `orpheus/geometry/reduced_operator.py` | ✅ **`sn/angular/`** (ruled §9.2; `transport/angular/` rejected — 0 candidate consumers) — ⛔ **but in P4, NOT P2**: `[M]` 2026-08-26 the move creates a hard `geometry → sn` import cycle while the three `*_streaming` factories still call `angular_redistribution`; reproduction at P2's ⛔ block. The destination is unchanged; only the phase moved. | `[M]` `alpha_dome(cosines, weights)` takes **no geometry argument at all**; the `1/r` lives in `ΔA` (the spatial factor), not in α. Chart-dependence is a *selection*, not spatial data. `[M]` the whole object is buildable from `(quad, coord)` with **no mesh** — and is already called that way on the multi-D path (`augmented_mesh.py:417`). |
| the `AngularMeasure` Protocol | same file | **dissolves** (in **P4**, with α) | Its own docstring gives the reason it exists: *"the geometry layer needs no import from the quadrature package at all."* ~~`[M]` trace which members survive moving α out — **4 of 6 go with it**, and the other 2 exist only because `StreamingTerms` bundles angular data into a geometry packet (one of those two is dead, §3).~~ ⛔ **RE-MEASURED 2026-08-26 — it is stronger than this row claimed, and the row's arithmetic was wrong.** The question is not which members *survive* but which the SURVIVING CODE READS, and `[M]` the three `*_streaming` factories that stay behind read **0 of 6** — they only forward the object. All four real reads (`mu_x`, `weights`, `mu_z`, `level_indices`) sit inside `angular_redistribution`'s body and travel with it; `N` and `eta` are read **nowhere in the file at all**. ⟹ the Protocol does not shrink to ≤2, it **dies**. ⛔ And its stated reason for existing is already false: `[M]` `geometry/boundary/` imports `Quadrature` from `orpheus.numerics.quadrature` under `TYPE_CHECKING` in **5** modules — so the replacement ships in its own package, five times over. ⟹ **a boundary Protocol whose surviving side reads none of it is the shadow of a misplacement, not a boundary.** |
| the `AngularClosure` family | `sn/sweep/` | `sn/angular/` | `sn/sweep/` is the **traversal** package (`scan.py`, `cache.py`, `pairing.py`). The closure is a *discretization* object the sweep consumes. §5c's hard guard says traversal is contraband in the scheme; this is the exact mirror. And it breaks R15's symmetry — `SpatialClosure` lives at `transport/spatial/`, its sibling at `sn/sweep/`. |
| `R` / the pairing | a property on the geometry object | minted by the **scheme** (§5) | `[M]` today it lives on the geometry object and **never reads the scheme** — correct only because DD's basis is the constant. F1 says the basis is the scheme's. |
| `StreamingTerms` | `geometry/reduced_operator.py` | stays in geometry — but **shed its angular fields** | `[M]` it is **not pure geometry**: `mu`, `abs_mu`, `mu_start` and the `/w` inside `delta_A_over_w` are all quadrature reads. Those four are exactly what forces the 6-member Protocol. |
| `StreamingCoefficientCache` (was `GeometryCoefficients`) / `CollisionCache` | `sn/sweep/` | **stay** — ✅ the name was REMEDIED in P3b 2026-08-26; ⛔ but §4bis then measured the CONTENTS to weld three strata, so "home right" is now the only surviving half of this verdict | `[M]` chain-ordered ⟹ genuinely traversal-time. This is the audit's D2 verdict: mostly NO smuggling, but the names lie. |

---

## 3. Un-welds — what is held versus what is used

| weld | what is held | what is used | verdict |
|---|---|---|---|
| `requires_upstream_angular_state`, `angular_marching_axis` | two fields on `ChartConnection` | `[M]` **0 production readers**; 12 test assertions. Made RED rather than grepped: flipping both on **997** constructed operators over 4 test trees gives **6 failures — exactly the 6 assertions that name the fields**. | ⛔ **RETIRE, not wire.** The concept is respelled twice already (`upstream_state.angular_upstream is None`, `SNMesh.is_cartesian`). |
| the three-link dead chain — ✅ **REMEDIED by P1 @ `ebe5d22f`** (both downstream links retired; ⚠ the third name below is now spelled `StreamingCoefficientCache`, P3b) | `AngularRedistribution.mu_start_per_level` → `StreamingTerms.mu_start` → `GeometryCoefficients.mu_start` → **∅** | `[M]` the terminal has **zero readers** of any kind (dynamic access checked). So `StreamingTerms.mu_start`'s only production consumer is *the write into it* — while its docstring claims `MorelMontryAngularSweep` consumes it, which reads **the owner** instead. | Retire the two downstream links; the owner stays. |
| `ChartConnection._quadrature` | a second reference to the measure | `[M]` `_quadrature IS angular.quadrature` → **True**. A twin reference to one object. | Retire; read through the angular factor. |
| `μ_start`'s edge-extrapolation branch | a code path + the `argsort` hazard **recorded in #361 (CLOSED 2026-08-13, `dde93b64`)** — ⚠ §9(b): a RECORD, not a work item; its own body is titled *"Why this was left alone rather than fixed inline"*, so `gh issue view 361` reading CLOSED is correct and means "decision documented", **not** "hazard removed" | `[M]` **structurally unreachable** on the probed fixtures (perturbation `0.0`, control `4.04e-01`): its only consumer runs on NON-carrying levels, and every level carries. ⚠ `[M]` on 2 fixtures — **a sample.** | ⛔ **BOTH THE VERDICT AND THIS ROW'S DECISION RULE ARE REFUTED — measured 2026-08-26**, see below. **KEEP** the branch. |

> ⛔⛔ **REFUTED 2026-08-26 — and the refutation is of the METHOD, not the
> number.** The row instructed: *"evaluate `consumes_independent_seed` over
> every rule `assert_carrying_quadrature` admits. All-carrying ⟹ #361 is a
> **retirement**, not a repair."* Run exactly as written it returns
> **0 of 88** — true, and it licenses the wrong conclusion.
>
> `[M]` **`assert_carrying_quadrature` has exactly ONE call site**
> (`augmented_mesh.py:347`), inside `case CoordSystem.CYLINDRICAL`. The
> `SPHERICAL` arm (`:352-356`) calls **no admission gate at all**, and
> admission is by *structure, not provenance* (the gate's own docstring). So
> the prescribed denominator was **cylinder-only** while the branch's live
> witness is a **sphere** rule: a hand-built **Gauss–Lobatto** polar rule
> builds a production `SNMesh(SPHERICAL)` and **fires the branch** at
> `n = 6, 8, 9, 10, 11, 17` — `[M]` **6 of 11 orders**.
>
> ⭐ **And the naive check is TAUTOLOGICAL**, which the census caught and
> reported rather than banking: the gate and `_carrying_levels` read the
> *same producer* on the same `(quad, coord)`, so "an admitted rule with a
> non-carrying level" is a contradiction in terms. Its `0` carries **zero**
> information. What decided the question was execution: positive control
> (forcing `_carrying_levels := frozenset()` moves state `rel 4.246e-01` /
> `3.555e-01`) plus a gate-bypass separation — **carrying 105 built / 0
> fired; non-carrying 13 built / 13 fired.**
>
> ⭐⭐ **The #361 hazard is ANNIHILATED where it is reachable, by theorem.**
> `_edge_seed_stencil:1696` computes
> `t = (μ_start − μ[m0]) / (μ[m1] − μ[m0])`, and `on_edge_node` *means*
> `μ_start == μ[m0]` — so the numerator is exactly `0` and the ambiguous
> `m1` is multiplied by `0.0`. `[M]` over 75 reachable non-carrying levels:
> tie-break live in **0 of 75**, `t == 0` exactly in **74 of 75**;
> perturbing slot `m1` by `1e3` on the production-admitted Lobatto witness
> moves the seed by **`0.000000e+00`** (control on `m0`: `1.0e+03`).
>
> ⟹ **the disposition is a THIRD state the row did not offer** — neither a
> live repair nor a clean retirement. **Keep** the branch (retiring it
> removes the only seed path a `μ=−1`-noded sphere rule has, and that rule
> class is live project work — `rules_sphere.py:876` names "a future
> Gauss–Lobatto sphere" in `MarchStart.on_edge_node`'s own docstring); write
> the measured fraction into its docstring; and land any `argsort` change
> **with a Lobatto witness row in the same step** (§6c — otherwise the change
> ships with no gate that can redden).
>
> ⚠ Two riders. The 5 Lobatto *refusals* are **round-off-gated**
> (`tau_raw[N-1] = 1 + O(1e-15)` against an exact `[0,1]` guard), not
> structural — so the 6-of-11 is a floor, not a property. And `M` is closed
> over **shipped factories** only: Gauss–Radau (same `on_edge_node` witness,
> unprobed), Gauss–Chebyshev, `quotient(g≠Mirror('y'))` and arbitrary
> hand-built `Quadrature(measure=…)` are named holes, the last genuinely
> unclosable.
| `ChartConnection` itself | three stages in one object | pure chart (`face_areas`, `delta_A` — need **no** quadrature); the head-stage angular pairing (extracted at Phase B); and `streaming_terms`' evaluated per-`(cell, ordinate)` view (**O-3's layer**) | Split three ways. The Phase B carve already thinned it toward this; the two dead flags are the residue. |
| `redistribution_pairing` allocates per access | — | `[M]` `g1 is not g2` → **True**: `delta_A[:, None, None]` re-allocates on every call. Harmless today (one call per mesh) but it is a *cache-shaped* object with no cache. ⚠ **MINE, from `6859ca05`.** | Compute once at the mint (§5). |
| `#407` — DD's blend inverse in scheme-neutral code | `2.0 = 1/w_DD` at `cell_balance.py:248/:343` | `[M]` `cell_balance_terms` has ONE production consumer (`diamond.py:194`); `linear_discontinuous.py` never imports either helper; and the matvec **declares** the specificity the producer does not (`loss_representation:3091`: *"Curvilinear (DD-only) … with DD's diamond march inlined"*). | The whole module is a DD body in a scheme-neutral name and home. Carve with O-3. |
| `#407`'s second instance | three more `2.0 = 1/w_DD` at `psi_half_angle_seed.py:180-185` | Latent, guarded only by `_require_slab` **in a different package, on the other tensor factor**. | ⭐ Sharper than the first: that march is the **angular** factor of a product `pairing.py` declares orthogonal. |

---

## 4. Strata and laziness — F3 made operational

The tree already has the right pattern and it is the model to extend: the
two-strata cache (`sn/sweep/cache.py`) separates *geometry × quadrature*
(survives Σ_t rebinds) from *× Σ_t* (rebuilt once per epoch). The question this
plan answers is **which stratum each object in §1–§3 belongs to.**

| stratum | invalidated by | members |
|---|---|---|
| **S0 — mesh-free angular algebra** | a new quadrature or chart ONLY | `α` dome, the angular cell partition, `τ`, `c_in`, `c_out`, `μ_start`, the carrying-level set. `[M]` bit-identical across 5 meshes incl. graded. |
| **S1 — chart × basis** | a re-mesh, or a scheme change | `face_areas`, `delta_A`, the **pairing** `R`, the closure's `dAw = R ⊗ 1/w`, the moment measure `M/V` |
| **S2 — × Σ_t** | a cross-section rebind | `CollisionCache` (already correct) |
| **S3 — traversal** | a sweep-order change | chain ordering, leg decomposition, the pole-mirror permutation, inflow/outflow partition |

⟹ **The laziness rule this yields:** an object is built at the stage of its
*latest* input and cached against its *earliest* invalidator. Today the angular
closure is built once per **mesh** while `[M]` its S0 half cannot change under
a re-mesh — that is the F3 weld, and it is why the split is worth doing rather
than merely tidy.

⚠ **What NOT to make lazy.** `CollisionCache`'s eagerness is load-bearing —
`_build_count` instruments a cardinal invariance gate (*exactly once per Σ_t
epoch*). Do not convert a gate-instrumented eager build into a lazy one; the
gate measures the build.

---

## 4bis. ⛔ The cache class welds THREE strata — measured, and it refuted both proposed names

`[M]` 2026-08-26, during P3b. Three meshes against **one** quadrature (uniform
`nx=6`, uniform `nx=20`, GRADED `nx=6`), comparing every field of the class §1
proposed to rename:

| stratum | fields | reading |
|---|---|---|
| **S0 — mesh-free** | `abs_mu`, `c_in`, `c_out`, `tau_inv`, `mm_a_in_coeff`, `is_degenerate`, `level_ordinates` | **bit-identical on all three** — 7 of 13 |
| **S1 — chart × basis** | `A_down`, `A_total`, `dA_w`, `V` | **differ on every re-mesh** — 4 of 13 |
| **S3 — traversal** | `chain_idx`, `chain_idx_inv` | identical between uniform and **graded** at equal `nx` ⟹ turn on ordinate sign and cell COUNT, not edge positions — 2 of 13 |

⟹ **The object is not one stratum, so no single-stratum name is true of it**, and
§1's own first candidate is the worst available: `ChainScanCoefficients` names
the **S3** half — **2 of 13** — where the `Geometry` it replaces at least names
4. A name asserting a stratum would *bless* the weld.

⛔⛔ **And §1's SECOND candidate is worse still: `SweepCoefficientCache` is the
name of a REFUTED DESIGN.** `.claude/lessons.md` **L15** records
`SweepCoefficientCache (N, nx, ng)` — a single cache mixing geometry and Σ_t
fields — as *"the wrong shape"*, the very monolith whose split produced this
class and `CollisionCache`. `[M]` no class by that name ever existed
(`git log -S` empty): it was a proposal, rejected. Landing it would make an
**always-session-loaded** lesson read as condemning the shipped design.
⟹ promoted to `plan-authoring` §1: *a proposed name's ABSENCE from the tree can
be evidence AGAINST it — grep the PROSE corpus, not just the code.*

✅ **Ruled `StreamingCoefficientCache`** (user, 2026-08-26). It claims only what
is true of all three strata, and it realizes the algebra **L15 itself states**:
*"L (streaming + curvature) lives in [this class]; C joins via
`1/(g_streaming + Σ_t·g_volume)` to form `CollisionCache`."* This is **L**'s
coefficient cache; its sibling is **L + C**'s — so the pair now names its
lifetime discriminator consistently, which `Geometry`/`Collision` did not.

⚠⚠ **The finding that outlives the name: L15's own lesson applies to L15's own
worked example.** L15 is *"cache shape that mixes immutability strata hurts
twice"* and it holds this class up as the RIGHT shape — while the measurement
above shows it mixes three. **A lesson's exhibit is not exempt from the lesson**,
and nothing prompts the re-check, because the exhibit is what taught you the
rule.

### ▶ The split earns a phase — **P4b**, not an issue-and-forget

Per §9.3's governing ruling. The whole object rebuilds on any re-mesh including
7 fields proven unable to change. **Done when:** the S0 half is cached against
`Quadrature × CoordSystem` alone and survives a re-mesh by identity
(`cache_a.s0 is cache_b.s0` for two meshes over one quadrature); the S1 half is
mesh-bound; `chain_idx` moves to the traversal layer. ⚠ Do **not** convert
`CollisionCache`'s eager build to lazy — §4 records that its `_build_count`
instruments a cardinal gate.

---

## 4ter. ⛔ The module in `geometry/` holds no geometry — and that re-poses P4

**Measured 2026-08-27**, in a session that set out only to price P4's α blocker.
The blocker dissolved because its premise was wrong: the question was never
*"where should this geometry type live"*.

### The measurement

`[M]` **Intra-`geometry/` consumers, by AST over geometry files OTHER than the
definer, code references only (a `__init__` re-export is not a consumer):**

| symbol | definer | other geometry files using it |
|---|---|---|
| `CoordSystem` | `coord.py` | 4 |
| `Mesh1D` | `mesh.py` | 2 |
| `StructuredGeometry` | `structured_geometry.py` | 1 |
| `RigidMotion` | `transformation.py` | 1 |
| **`StreamingTerms`** | `reduced_operator.py` | **0** |
| **`ReducedStreamingOperator`** | `reduced_operator.py` | **0** |
| **`AngularRedistribution`** | `reduced_operator.py` | **0** |

Every genuine geometry primitive has 1–4 intra-package consumers. Everything in
`reduced_operator.py` has **zero** — it is an island inside its own package.
⚠ This measurement DISCRIMINATES only because the control row exists; "0 intra-
package consumers" would otherwise read as normal for an input layer.

`[M]` **What the module's fields actually are** (`spherical_streaming` /
`cylindrical_streaming` bodies, 4 and 5 executable statements):

```python
face_areas = mesh.areas                        # a VERBATIM COPY of a Mesh1D property
delta_A    = face_areas[1:] - face_areas[:-1]  # np.diff of that same array
coord      = mesh.coord                        # the guard reads mesh.coord two lines above
mesh       = mesh                              # …and it holds the source anyway
```

`[M]` `delta_A == ∫(∇·ê_r) dV` per cell to **8.9e-16** (sphere, `_spherical_mesh()`,
`gauss_legendre(8)`) — the contracted connection coefficient.

`[M]` **Who wants which**, ten filter families over 334 `.py` files under
`orpheus/`, every filter validated against a literal of its own shape BEFORE its
zero was read, tree control = the two known `delta_A` sites:

| datum | sn | geometry | transport | cp | moc | mc | diffusion | numerics |
|---|---|---|---|---|---|---|---|---|
| `delta_A` | 9 | 24 | **0** | **0** | **0** | **0** | **0** | **0** |
| `face_areas` | 6 | 23 | 0 | 0 | 0 | 0 | **2** | **8** |
| `.areas` | 4 | 3 | **5** | **2** | 0 | 0 | **4** | 0 |

⟹ `face_areas` is **shared geometry** with four independent non-SN consumers —
`cp/solver.py:412` (`S_cell` in the escape/entrance reciprocity),
`numerics/spaces/scalar_trace_space.py` (the trace space is *built from* a face-area
inventory — the Stokes restriction), `diffusion/`, `transport/`. `delta_A` has
**none**, and its absence is STRUCTURAL, not incidental: a conservative
finite-volume method needs `A_{i±1/2}` *separately* (`A₊J₊ − A₋J₋`), while
curvilinear SN needs their *difference*, because `ΔA` multiplies an ANGULAR
difference `[α ψ₊ − α ψ₋]`. The difference is formed only where one geometric
factor multiplies an angular difference.

⭐ **The discriminator is already in the code**: `StreamingTerms` carries BOTH
forms side by side — `face_area_inner`/`face_area_outer` (separate, spatial
balance) and `delta_A_over_w` (the difference, angular redistribution) — and the
Cartesian arm sets `delta_A_over_w = 0.0` with the comment *"no curvature
redistribution"*.

⛔ **`StreamingTerms`' docstring claims "Per-(cell, direction) purely geometric
inputs". It carries `mu`, `abs_mu`, and `ΔA` divided by a QUADRATURE WEIGHT**
(`:781`, `:812` — the divisor is `_weight_of` → `angular.quadrature.weights[n]`).
A per-direction quantity over a quadrature weight is not geometry.

### The ruling — the layer test, and what it forces

> **A datum belongs to the layer that can define it without naming a method.
> Everything else is posing, and posing belongs to the method head that poses it.**

`[M]` exactly one datum in the module survives that test — `face_areas` — and it
is a copy of something already single-sourced in `geometry/coord.py`
(`compute_areas_1d`: CARTESIAN → `ones_like(edges)`, CYLINDRICAL → `2πr`,
SPHERICAL → `4πr²`). ⟹ **the module dissolves.**

| what | goes to | layer | why |
|---|---|---|---|
| `face_areas` | **retired** — read `mesh.areas` | INPUT | already single-sourced; 4 non-SN consumers read it there |
| `coord` | **retired** — read `mesh.coord` | — | duplicate, and what lets `transport` reach THROUGH the bundle |
| `StreamingTerms` | `transport/spatial/` | **L2** | its 2 runtime importers are `cell_balance.py:79`, `scheme.py:71` |
| the five α symbols | `sn/angular/redistribution.py` | **L3** | α is chart × measure at the intersection; only a collocated angular derivative forms it |
| `AngularMeasure` | **retired** | — | factories read 0 of its 6 members; its stated reason is already false (`geometry/boundary/` TYPE_CHECKING-imports `Quadrature` in 5 modules) |
| `delta_A`, `redistribution_pairing`, `streaming_terms()` | `sn/` | **L3** | `np.diff(mesh.areas)`, derived where its only consumer lives |

⭐ The forbidden edge never arises — not routed around, not whitelisted, not
`TYPE_CHECKING`-ed. `geometry` stops NAMING an SN type because it stops HOLDING
one.

### The step order, and why it is forced

**P4.1a — `coord` retires.** ✅ **LANDED `d48f4bf4`** 2026-08-27 (exit gate: 13 trees `rc=0`,
9823 passed / 0 failed, delta 0 on skipped/deselected/xfailed). `[M]` duplicate on **3/3**
charts (`op.coord is op.mesh.coord`). Bit-identical. **This is the load-bearing
prerequisite**: it is what severs `transport`'s reach into the bundle
(`transport/radial_characteristic_field.py:361-362` reads `mesh.reduced.coord`
→ `mesh.coord`). Without it, P4.4 breaks an L2 consumer.

⭐ **Execution upgraded the 3/3 measurement to a THEOREM, and that is worth
carrying:** each of the three factories *validates* `mesh.coord` against the
literal it then stores (`slab_streaming` raises unless
`mesh.coord is CARTESIAN`, and so on), so `op.coord is op.mesh.coord` held **by
construction** — no reachable state could have made them disagree. `[M]` all
three guards confirmed by AST. The fixture measurement was true and weaker than
the fact.

⭐⭐ **And the finding the step paid for, which is a `coding-standards` MIRROR
instance nothing prompted: the retirement PROMOTED those three guards from
redundant to load-bearing, and `[M]` they had ZERO witnesses tree-wide**
(`grep "requires .* mesh"` → only the 3 production lines). Before P4.1a a broken
guard was survivable — the stored literal and the mesh said the same thing two
ways. After it, the guard is the *sole* reason `op.mesh.coord` is the operator's
chart. Discharged in the same commit: `TestProperties`' three chart tests were
tautologies (a stored literal read back) and are now the guards' witnesses in
`vv-principles` #11 form — one positive leg, two negative legs each, matching on
the production message.

`[M]` **the actual blast radius, by AST + a validated string-form filter** (not
the 2 sites the plan named): **13 production reads across 7 files** —
`transport/radial_characteristic_field.py` (1, the only L2 site),
`sn/loss_representation/__init__.py` (2), `sn/solver.py` (2),
`sn/coupled_system.py` (1), `sn/sweep/cache.py` (1),
`sn/mesh/augmented_mesh.py` (4), plus the 4 in-module `self.coord` reads inside
`streaming_terms()`. **Every consumer had an `SNMesh` in scope**, so not one
site needed `reduced.mesh.coord`: the reach-through died at all seven files, not
only at `transport`. Four dead preamble lines (`reduced = sn_mesh.reduced` +
its narrowing assert) died with the reads they narrowed, at `solver.py` ×2 and
`transport` ×1. `sn/sweep/cache.py`'s assert **stays** — it is an admission
contract on `reduced is None`, not narrowing (⚠ and it is a bare `assert`, so
`-O` strips it — P4.5 owns that predicate, see below).

⚠ **The §8 question this step DOES raise, and its answer.** Re-pointing
`self.reduced.coord` → `self.coord` is not a pure re-spelling: the old
expression **raises** when `reduced is None` (d≥2 Cartesian) and the new one
returns `CARTESIAN`. So the edit could silently turn a latent crash into a
proceed at every site — §8's "an enabler still has its own blast radius", in the
direction nothing complains about. `[M]` **all 7 sites are provably unreachable
with `reduced is None`**, and the proof is from the guards, not from a fixture:

| site | why `reduced` cannot be `None` there |
|---|---|
| `augmented_mesh.radial_characteristic_levels` | `if self.is_cartesian: return` fires first |
| `augmented_mesh.dag_walk` | `raise ValueError` on `reduced is None`, above the read |
| `augmented_mesh.dag_walk_cell_indices` | same `raise` |
| `augmented_mesh._representative_ordinate` | private; `[M]` its only 2 callers are the two above |
| `loss_representation._OneDimScanWalk._run` | a **1-D-only** frame, and every 1-D mesh has a `reduced` (d≥2 Cartesian is `_OctantWalk`) |
| `…._run_transpose` | same frame |
| `transport.source_from_angular` | raises above on `vals.ndim != 3` — carrying ⇒ 1-D curvilinear |

⟹ no reachable input distinguishes old from new. ⚠ Note the two bare `assert
… is not None` lines that were doing this narrowing are **stripped by `-O`**, so
they were never the reason; the enclosing `raise`/branch is.

⚠ **Test migration, per `coding-standards`.** `test_snmesh_consumes_reduced.py`
stated its invariant 1 as *"a `ReducedStreamingOperator` **with the matching
`CoordSystem`**"* — a proxy for "the right factory ran" that P4.1a makes
tautological (`reduced.mesh is sn.mesh`, so it would restate `sn.coord`). It now
asserts `sn.reduced.mesh is sn.mesh` — the wiring claim the file's own header
says it is for, strictly stronger than a chart tag (which agrees for *any* two
meshes on the same chart) and chart-independent, so P4.1b does not move it.

**P4.1b — the spatial fields take their neutral element.** ⛔ **NOT
bit-identical — this plan asserted it was, and the check refuted it.**
`[M]` `face_areas` is a duplicate on **2 of 3** charts: on the slab the field is
`None` while `mesh.areas` is `[1,1,1,1,1,1]`, because `compute_areas_1d` returns
a real unit cross-section for CARTESIAN and `slab_streaming` never asks for it.

⭐ **What this step actually is: finishing P1's un-weld on the SPATIAL half.**
P1 converted the ANGULAR fields from `None` to the neutral element on
2026-08-26 — its gate calls that *"a STRICTLY STRONGER claim … the neutral
element says the dome IS identically zero"* and *"the per-coordinate `Optional`
union died with it (Pattern 4)"*. The spatial fields were left as `None`, and
`tests/geometry/test_reduced_operator.py:343` records the asymmetry two lines
away (*"The SPATIAL chart is still absent for a slab"*) without closing it.

`[M]` **blast radius, by the §8 grep** (`is None` / `is not None` / `getattr`
default over `orpheus/` + `tests/`, control = 9 branches in the definer): **14
branches** —
* **1** production control-flow branch: `reduced_operator.py:694`,
  `if self.delta_A is None: return np.zeros((nx,1,1))`;
* **6** type-narrowing asserts (legitimate per `coding-standards` — `-O` strips
  them and the downstream failure is an immediate `AttributeError`);
* **2** test pins that go RED: `test_reduced_operator.py:343-344`;
* 5 curvilinear asserts, unaffected.

`[M]` **the special case dissolves bit-exactly**: feeding the general body
`face_areas=mesh.areas, delta_A=np.diff(mesh.areas)` on the slab gives
`array_equal=True` against the `None` branch (both `(5,1,1)`, `max|·| = 0.0`).
⟹ the branch is a value wearing a conditional (Pattern 4).

⭐⭐ **RE-SCOPED 2026-08-27, measured during P4.1a — P4.1b is not a two-field
un-weld, it is a Pattern-2 arm collapse, and THAT is the prize.** The plan
framed this step as Pattern 4 (an `Optional` dies). Reading
`streaming_terms()`'s three arms at HEAD shows the Pattern-2 half nobody had
written down: **the CARTESIAN arm differs from the SPHERICAL arm only in three
hardcoded literals**, and those literals are hand-transcriptions of what the
spherical body computes.

| slot | slab arm | sphere arm | slab value once populated |
|---|---|---|---|
| `face_area_inner` | **`1.0`** literal | `face_areas[i]` | `areas[i] = 1.0` ✓ |
| `face_area_outer` | **`1.0`** literal | `face_areas[i+1]` | `areas[i+1] = 1.0` ✓ |
| `delta_A_over_w` | **`0.0`** literal | `delta_A[i]/w` | `0.0/w = 0.0` ✓ |

every other slot (`chord_length`, `mu`, `volume`, `abs_mu`) is already spelled
identically, and both arms read `direction_idx` as the global ordinate.

`[M]` `scratch/_p41b_arm_collapse_probe.py`, 2026-08-27 — populate the slab's
two spatial fields, force the packet down the SPHERICAL body, and compare the
whole `StreamingTerms` tuple over the 5×8 (cell × ordinate) grid: **40 of 40
bit-identical, 0 mismatches**, with a positive control (a perturbed `delta_A`
IS detected). ⟹ **the CARTESIAN arm is provably redundant.**

⭐ And the file already knows the principle — it states it 25 lines below the
code that violates it, about the α-dome: *"until 2026-08-12 each spelled the
recursion out in its own body (sphere and cylinder, identical arithmetic — the
sphere IS the cylinder's single-level case). That twin path is …"*. The same
sentence is true one level up: **the slab IS the sphere's zero-curvature case**,
and sphere IS cylinder's single-level case, so the three-way `if` reduces to one
chart-dependent step — *resolve the global ordinate and the radial direction
cosine* (`mu_x[direction_idx]` vs `eta[level_indices[p][direction_idx]]`) — over
one shared body.

⟹ **Why this belongs to P4.1b rather than to a later step:** populating the
fields without collapsing the arm leaves a branch that is dead by construction
(`coding-standards` retire-as-you-go), and the collapse is what makes P4.4's
"move the residue to `sn/`" a move of **one** body instead of three. It also
makes P4.7's retirement of `StreamingTerms.delta_A_over_w` a one-place edit.

✅ **RULED 2026-08-27 (user): collapse the arm.**

#### How far the collapse actually goes — measured, and it is not "all the way"

⛔ **A first design read said the three arms collapse to ONE body with NO chart
dispatch. That was wrong, and the check that caught it is the §2 quantifier
check.** The reasoning was sound from the *definitions*: `[M]`
`Quadrature.eta` and `Quadrature.mu_x` are **the same accessor** (both
`axis_cosines(0)`; `eta`'s own docstring says *"Same data as `mu_x` (column
0)"*), and `[M]` `level_indices` returns `[arange(N)]` **by construction**
whenever `level_structure is None`. If both held universally, then
`level_indices[p or 0][direction_idx]` would reduce to `direction_idx` on
slab and sphere and one expression would serve all three charts.

`[M]` over the shipped registry — **40 rules constructed** (GL 2–18,
Gauss–Lobatto 2–18, `level_symmetric` 2–18, `product`/`folded_product` over
`n_mu ∈ {2,4,6} × n_φ ∈ {4,6,8,16}`), positive control included: **9 are
single-level, and 1 of those 9 has a PERMUTED level list.**
`level_symmetric(2)` is `N = 8`, one level, `level_indices[0] =
[2,3,0,1,6,7,4,5]` — it carries a `LevelStructure`, so the `arange` shortcut
does not apply to it. ⟹ **`level_indices[0][direction_idx]` is NOT
interchangeable with `direction_idx`**, and unifying through it would silently
re-index any slab/sphere posed on that rule (nothing in the factories forbids
it — they accept "any `AngularMeasure`-shaped object").

⭐ **What the surviving `if` is actually about, and why keeping it is the
finding rather than the compromise.** `direction_idx` means **the global
ordinate** for slab/sphere and **the within-level azimuthal index** for
cylinder. That is one parameter name carrying two contracts — the real defect,
and it is *not* a chart-dispatched arithmetic. So the honest end state for
P4.1b is **3 arms → 1 body + a 2-way index resolution**:

```python
if self.mesh.coord is CoordSystem.CYLINDRICAL:
    if mu_level_idx is None:
        raise ValueError(...)                    # the admission stays
    n = int(quad.level_indices[mu_level_idx][direction_idx])
else:
    n = direction_idx                            # already the global ordinate
radial_cosine = float(quad.mu_x[n])              # == quad.eta[n], one accessor
return StreamingTerms(...)                       # ONE body, all three charts
```

⟹ the packet arithmetic is single-sourced (the Pattern-2 win), and what remains
is a named statement about *what the caller's index means* — which P4.7 owns,
since it is the same producer-contract question as the `delta_A_over_w`
retirement.

#### ⭐ The MEANS changes too — they were never fields

The plan's proposed means was *"populate the two fields on the slab"*. Executing
it that way would **create a Pattern-2 triplicate**: `[M]` the sphere and
cylinder factory bodies are already character-identical —

```
face_areas = mesh.areas
delta_A    = face_areas[1:] - face_areas[:-1]
```

— so "populate the slab too" means writing those two lines a **third** time, in
the step whose whole point is that the values were always derivable. That is
`coding-standards`' clean-before-extend read backwards.

⟹ **`face_areas` and `delta_A` stop being fields and become derived accessors
on the mesh** — `face_areas` a plain `@property` returning `self.mesh.areas`
(`[M]` an eagerly-computed frozen attribute on `Mesh1D`, so this is one
attribute hop, no recompute, and the same object the 4 external
`.reduced.face_areas` readers already get), and `delta_A` a
`functools.cached_property` over `np.diff` (⚠ **must** be cached — `streaming_terms`
is per-(cell, direction), so a plain property would recompute an `nx`-element
diff in the hot path).

This serves the ruled GOAL by a better means (`plan-authoring` §1 — the goal is
*"P1's un-weld finished on the spatial half; the per-coordinate `Optional`
union dies"*; "populate the fields" was the proposed mechanism):

* the `Optional` dies, and so does `delta_A_column`'s `if self.delta_A is None`
  branch — **Pattern 4**, same as P1's angular half;
* **no factory changes at all**, and all three bodies collapse to one shape
  (guard + `angular_redistribution(measure, CHART)`), which is what P4.4 wants;
* it is strictly *closer* to this plan's own ruled end state — the §4ter table
  says `face_areas` is **retired, read `mesh.areas`** — while keeping the 4
  external readers and `SNMesh`'s deprecated accessor working, so the re-point
  stays P4.4's;
* ⭐ and it removes a latent defect nobody filed: `@dataclass` generates
  `__eq__` over the fields, so `op1 == op2` on two ndarray fields raises
  *"truth value of an array is ambiguous"*. Dropping them from the field list
  makes the comparison well-defined.

⚠ `[M]` the only 3 construction sites are the factories themselves, and there
are **0** `dataclasses.replace(...)` calls naming either field, so removing them
from the constructor signature has no other consumer.

#### The §8 check on P4.1b's OWN blast radius — the collapse widens a data dependence

The collapse is value-identical (40/40), but it is not dependence-identical, and
§8 asks what a step CHANGES rather than what it computes. The retired CARTESIAN
arm returned the literal `delta_A_over_w = 0.0` and **read no quadrature weight
at all**; the shared body computes `delta_A[i] / self._weight_of(ordinate)`. So
a slab packet now depends on `quad.weights[n]` where it did not — and `0.0 / w`
is `0.0` for any finite non-zero `w`, but `NaN` for `w = 0`.

`[M]` 2026-08-27 over the 40-rule shipped registry, positive control included:
**0 rules carry a zero or non-finite weight**, and the smallest `|w|` anywhere
is `1.750e-04` (`level_symmetric(18)`). ⟹ safe on everything that ships.

⚠ And the reason this needs no guard and no code comment: a zero weight is
already pathological for the curvilinear charts, where the numerator is *not*
zero and `ΔA/0` diverges. The slab is the **more** robust case, not a newly
fragile one — the collapse extends an existing precondition to a chart where it
happens to be satisfied trivially. Recording the measurement, not a hazard.

#### The coverage check `coding-standards` demands before this lands — run, and it inverts

Single-sourcing a duplicate obliges the *"hunt for an EXTERNAL hand-written pin
before concluding you lost coverage"* step, **and** its rider: check the pin
**moves under the old value**, or it is blind.

`[M]` 2026-08-27 — the pin exists, twice, and neither is blind:

* `tests/geometry/test_geometry.py:113` pins the producer directly against a
  literal written independently of everything here —
  `compute_areas_1d(CARTESIAN, [0, 1, 3]) == [1.0, 1.0, 1.0]`;
* `tests/geometry/test_reduced_operator.py:425-427` **and**
  `tests/sn/sweep/core/test_discretization_scheme_protocol.py:403-406` (a
  different tree) each pin the slab packet: `face_area_inner == 1.0`,
  `face_area_outer == 1.0`, `delta_A_over_w == 0.0`. Under the collapse those
  values come from `mesh.areas`, so a wrong `compute_areas_1d` reddens them —
  they move.

⭐⭐ **And the two packet pins are PROMOTED by this step — the MIRROR clause
again, one phase after P4.1a hit it.** Today they are tautologies: the CARTESIAN
arm *returns the literal `1.0`* and the test *asserts `1.0`*. After the collapse
they assert that the areas path flows end-to-end through the unified body. Same
assertions, no line of either test body changing, claim class strengthened.

⟹ **re-describe both docstrings/comments in the same commit.** Their current
comment — *"Neutral curvature: slab carries the values that make the unified
cell-balance algebra collapse to the slab form"* — describes a hardcoded arm
that will no longer exist.

⚠ Note the symmetry with P4.1a, because it is the reason to look both times:
P4.1a's promotion landed on a **guard**, which had *no* witness and needed one
written; P4.1b's lands on **gates**, which *are* the witnesses and need only an
honest description. Same clause, opposite deliverable.

⚠ **Naming item for P4.7's family.** The unified body needs a chart-neutral
spelling of `axis_cosines(0)`: `mu_x` is slab-flavoured, `eta` is
cylinder-flavoured, and `[M]` they are one accessor with two names. `mu_x`'s
own docstring already concedes it (*"the column index, not the name, is the
actual semantic"*). Bind a well-named local (`radial_cosine`) for now;
⛔ do **not** reach for `axis_cosines(0)` directly — it is not on the
`AngularMeasure` Protocol, so using it would narrow what the factories accept,
and that Protocol's fate is suspended into CS5.

**P4.1c — the `SNMesh` deprecation shims retire.** ✅ **UNBLOCKED NOW** — it
depends on nothing in P4 and can land beside P4.1b.

`SNMesh.face_areas` and `SNMesh.delta_A` are `@property` shims that emit
`DeprecationWarning` and route to `self.reduced`. They date from Wave E Round 2
(#164) — many merge cycles ago, against `coding-standards`' *"deprecation
aliases live for **one merge cycle only**"*.

`[M]` 2026-08-27, by AST over `orpheus/` + `tests/`: **11 reads, 0 of them in
`orpheus/`.** Every consumer is a test, and `[M]` the tests are
`test_sphere_deprecated_properties_warn` /
`…_route_to_reduced` / `…_cylinder_…` — i.e. **the only thing keeping the shims
alive is the suite that verifies the shims.** That is the closed loop
`coding-standards` describes: a shim with no consumer, wearing coverage.

⟹ delete both properties and the 4 tests that exist only to exercise them
(API-smoke, per the test-migration rule — there is no behavioural contract to
rewire, since `reduced.face_areas` is the thing they route to and it is pinned
directly elsewhere). ⚠ Grep the **warning message strings** as well as the
symbol — the shortest distinctive fragment is `is deprecated; use SNMesh` — and
re-read `test_snmesh_consumes_reduced.py`'s header, whose invariant **2** and
**3** are *about* these accessors and retire with them.

**P4.2 — the angular factor comes home.** The five α symbols to
`sn/angular/redistribution.py`; `SNMesh` builds the redistribution and hands both
tensor factors to the closure — which it **already does** on the d≥2 Cartesian
path at `augmented_mesh.py:412`. The 1-D curvilinear path is the outlier being
brought into line, not a new pattern.
⚠ `derivations/discrete/sn/angular_differencing.py:164` imports `alpha_dome` at
module scope and `derivations → sn` is FORBIDDEN with no tolerance. The fix is
the precedent that file's own ⛔ banner set for `tau`/`edges` on 2026-08-12:
**the L0 ladder accepts α as a keyword** (*"the defect was L0 reaching UP for
it"*). Its delegating local `alpha_dome` retires; nothing re-implements it.

✅ **LANDED `5940deba`** 2026-08-28. Six symbols to
`orpheus/sn/angular/redistribution.py`; `geometry/reduced_operator.py` is down
to `StreamingTerms` alone and `[M]` now imports **nothing but `dataclasses`** —
a file in `geometry/` that names no geometry, which is the whole finding the
arc was built on. P4.3 deletes it.

⭐ **The plan said FIVE movers; it is SIX.** `AngularMeasure` had to travel:
`[M]` `AngularRedistribution` annotates a field with it, so leaving it behind
is `geometry` naming an `sn` type — the same forbidden edge in miniature. This
was caught by the pre-execution memo, not during execution.

⭐ **The L0 fix, and the thing worth carrying: `w` was only ever a means to α.**
`[M]` in all three graders `w` fed **nothing** but the `alpha_dome` call, so
injecting α made it redundant — it left `morel_montry_beta`'s signature
entirely, and `alpha_defect_beta` lost `quad`/`geometry` too (it is
`β = α − (1−e²)`, a function of α and the edges alone). ⟹ *the keyword fix is
not additive; it reveals which parameters were only ever a route to the thing
you are now handed.* Leaving them would have made dead parameters read as
load-bearing.

⭐⭐ **A retired gate was PROMOTED rather than deleted, and its claim class
changed.** `test_the_derivations_alpha_dome_IS_the_production_one` pinned that
the L0 wrapper *delegated* — a comparison no input could fail, and its own
docstring said so. With the wrapper retired the twin is **unspellable** (L0
cannot import the production recursion at all), so the gate became an
**anti-twin** gate: the module must define no `alpha_dome`, and all three
graders must take α as a **required keyword-only** parameter. `[M]`
mutation-checked in-process — 3 arms (twin re-added / default added / made
positional) all redden, positive control passes. `coding-standards`' mirror
clause, applied deliberately for once instead of caught late.

⚠ `[M]` the production α is **bit-identical** to what the L0 module computed
internally, on both charts (GL(8) sphere, folded(4,8) cylinder) — checked
BEFORE editing the tests, so no test value moved.

⭐ `dead_references` caught a dead `:attr:` **I introduced myself** in the new
docstring (`StreamingCoefficientCache` is in `sn/sweep/cache.py`, not
`sn/mesh/reduced_operator.py`). Same instrument, same session, now catching the
author rather than the inheritance — 1 → 0 dead / 52 checked.

✅ **The blocker record (kept, now past tense) — P4.4 landed at `16501ca0`, which is what unblocked this.**
`[M]` the three `*_streaming` factories are in `orpheus/sn/mesh/reduced_operator.py`;
their call to `angular_redistribution` is `sn → geometry` today and becomes
`sn → sn` the moment α moves. The ordering block below is kept as the record of
why it could not go first (`plan-authoring` §3), and is now PAST TENSE.

⚠ **What P4.2 still owes, unchanged:** `AngularMeasure` travels WITH α (it is
`AngularRedistribution`'s field type); the `derivations → sn` hazard bites here
for the first time and the fix is the **keyword**, never the `WHITELIST` row;
and the destination `orpheus/sn/angular/redistribution.py` still does not exist.
⚠ `geometry/__init__.py` re-exports `AngularMeasure` and `alpha_dome` — both
must be **DELETED**, not re-pointed.

⛔⛔ **ORDERING REFUTED 2026-08-28 (the record) — P4.2 RUNS AFTER P4.4.** The reconciliation
below is about the `derivations → sn` hazard and it still holds verbatim. What
it never asked is whether the move is landable *at this point in the order*, and
`[M]` it is not. The α cluster's callers are **inside its own file**:

| caller (stays in `geometry/` until P4.4) | names | kind |
|---|---|---|
| `slab_streaming` `:1198` | `angular_redistribution(...)` | **runtime call** |
| `spherical_streaming` `:1262` | `angular_redistribution(...)` | **runtime call** |
| `cylindrical_streaming` `:1360` | `angular_redistribution(...)` | **runtime call** |
| `ReducedStreamingOperator.angular` `:628` | `AngularRedistribution` | field annotation |

⟹ moving the α cluster out first makes `geometry/` import `sn/` at module
runtime. `[M]` `FORBIDDEN_EDGES["geometry"] = L2|L3` ∌ tolerance (`geometry` is
INPUT; the `TYPE_CHECKING` carve is `L1|L2`). **Injected and run** rather than
reasoned: the gate reddens on
`test_no_forbidden_imports[…reduced_operator.py]` **and** on a circular import
that kills `import orpheus.geometry` / `orpheus.numerics` / `orpheus.sn.solver`
in a fresh interpreter. File restored byte-identical (`diff -q` + `git diff
--quiet`).

⭐ **Why the thorough check missed it:** §6d says *enumerate the mover's
callers*. Three of these four are `angular_redistribution(...)` calls in the
same file — there is no import line to grep, because **the split is what
creates the edge**. The audit that found the `derivations` hazard found it
precisely because that one IS an import. ⟹ intra-file references are the
caller set too; enumerate by AST.

⟹ **After P4.4 this is clean**: the factories are in `sn/`, so the call becomes
`sn → sn`; `AngularMeasure` travels with the α cluster here (it is
`AngularRedistribution`'s field type at `:1082`, not only the factories'
parameter type); and the `derivations → sn` hazard below bites for the first
time, which is where the keyword fix lands.

✅ **RECONCILED against the tree 2026-08-27 — the claim holds, verbatim.**
`[M]` by AST: `FORBIDDEN_EDGES` has **15** keys and `"derivations"` is one of
them (`tests/test_layer_imports.py:71`), mapping to `L2_PACKAGES |
L3_PACKAGES` — so `sn` is forbidden. `[M]` the TYPE_CHECKING tolerance is
`if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES)` and
`L0_PACKAGES = {"derivations"}`, so **derivations is not covered even for
typing** — the same shape as the `geometry` finding that killed P2's α half,
one layer down. `[M]` the import is live at `:163-164` and the module still
re-exports a delegating `alpha_dome` at `:231` whose body is
`return _production_alpha_dome(mu, w)`.

⛔ **The move's inventory is 5 symbols, and they are not all "α" symbols** —
`[M]` top-level in `reduced_operator.py`: `alpha_dome`,
`_ALPHA_CLOSURE_ATOL`, `_assert_alpha_dome_closes`, **`AngularRedistribution`**
and **`angular_redistribution`**. The last two are the factor and its producer;
a summary that reads "the five α symbols" and greps `alpha` finds **3 of 5**.

⭐ **There IS a second, WORSE option, and it is worth naming so nobody
rediscovers it as a shortcut.** `[M]` the linter ships a `WHITELIST` whose
existing entries are all `derivations → L2/L3` exemptions
(`cases/diffusion.py → diffusion`, `mms/moc.py → moc`,
`sood_registry/builders.py → cp`, `mms/sn.py → transport`), so "add a whitelist
row" is an established idiom and would make the move land green in one line.
⟹ **Do not.** The whitelist records a transitional violation; the keyword fix
removes the violation. The file's own banner already ruled the principle for
`tau`/`edges` — *the defect was L0 reaching UP* — and taking the exemption would
re-open, for α, exactly what that banner closed for τ.

⚠ `[M]` the d≥2 precedent the paragraph above cites is real and now sits at
`augmented_mesh.py:408-419` (not `:412`): `if self.reduced is not None:` reads
both factors off the bundle, `else:` builds them from `(quad, coord)` alone —
with the comment *"That they are buildable from ``(quad, coord)`` alone is the
un-weld's own point: the closure's operands were never mesh facts."*

**P4.3 — `StreamingTerms` to L2.** ⛔ ~~`transport → geometry` stays legal (it
reads `mesh.areas`); no new edge either way.~~ **REFUTED 2026-08-27 — the
"either way" half is false, and it is §6d's own failure mode.** Fix the
*"purely geometric"* claim while the file is open.

`[M]` by AST: **`geometry → transport` = 0** today and **`transport → geometry`
= 16**, so the move runs against a 16:0 gradient — the same shape as the 24:0
that made P2's α half unlandable. And the mover's producer **stays behind**:
`ReducedStreamingOperator.streaming_terms` lives in `geometry/` and
`[M]` **CONSTRUCTS** `StreamingTerms(...)` at `reduced_operator.py:856`, inside
the function body — a **runtime** use, so the `TYPE_CHECKING` tolerance could
not save it even if `geometry` were covered by it (it is `INPUT_PACKAGES`, so it
is not). `FORBIDDEN_EDGES["geometry"] = L2_PACKAGES | L3_PACKAGES` and
`transport` is L2 ⟹ **a declared violation with a red waiting for it.**

✅ **REMEDIED 2026-08-28 by P4.4 (`16501ca0`) — this blocker is DISCHARGED, and
the paragraph above is now PAST TENSE.** `[M]` `StreamingTerms`' only
constructor is `orpheus/sn/mesh/reduced_operator.py:531`; it is in `sn/`, and
`sn → transport` is the established direction (`[M]` 171 edges). The consumers
re-point cleanly: `transport/spatial/scheme.py:92` (RUNTIME) and
`cell_balance.py:79` (TYPE_CHECKING) become `transport → transport`;
`sn/mesh/augmented_mesh.py:32` and `sn/mesh/reduced_operator.py:207` become
`sn → transport`, legal.
⚠ **`geometry/__init__.py:30` still re-exports `StreamingTerms`, and it must be
DELETED, not re-pointed** — `geometry → transport` is forbidden, the identical
trap P4.4 hit with its own four names.

⟹ **ORDERING CONSTRAINT: P4.3 cannot precede P4.4.** P4.4 re-homes
`streaming_terms()` into `sn/`; once the producer is there the import is not
merely legal but the established direction (`[M]` `sn → transport` = **171**).

#### ✅ And the REVERSE dependency, checked 2026-08-28 — reordering is clean

⚠ The clause above asserted a one-way dependency after measuring only one
direction, which is the very defect §6d exists to name. Checked properly:
**P4.4 does NOT depend on P4.3**, so the two reorder rather than needing to
fuse.

`[M]` `FORBIDDEN_EDGES["sn"] = L3_PACKAGES - {"sn"}` and `geometry` is
`INPUT_PACKAGES` ⟹ **`sn → geometry` is LEGAL.** So the interval state — the
producer in `sn/`, `StreamingTerms` still in `geometry/` — is clean in every
direction:

| interval edge | verdict |
|---|---|
| `sn` → `geometry` (the producer constructs `StreamingTerms`) | legal — `geometry` is INPUT |
| `transport` → `geometry` (`cell_balance.py:79`, `scheme.py:92`) | unchanged from today (`[M]` 16 such edges) |
| `geometry` → anything L2/L3 | none — the residual file needs only `Mesh1D` / numpy |

⛔ **But P4.4's SCOPE changes twice, and both are easy to miss:**

**(a) The module DELETION moves out of P4.4 and into P4.3.** P4.4's row says
*"the residue to `sn/`, **and the module is deleted**"* — it cannot be, because
`StreamingTerms` is still defined there and `[M]` `transport` imports it at
**runtime** from that exact module. The file survives P4.4 holding
`StreamingTerms`, and dies when P4.3 moves it out.

**(b) `geometry/__init__.py`'s re-export must be DELETED, not re-pointed** —
and this is a §6d item for **P4.4 itself**. `[M]` `geometry/__init__.py:30`
re-exports 7 movers; re-pointing that line at `sn/` would be `geometry → sn`,
**forbidden**. `[M]` its consumers are **11 files** — 3 in `orpheus/sn/`
(`augmented_mesh.py:32`, `radial_characteristic.py:142`,
`angular/closure.py:176`) and 8 under `tests/` — every one of which re-points
to the new home in the same step.
⚠ `[M]` the linter walks `ORPHEUS_ROOT = .../orpheus` only, so the 3
`tests/transport/spatial/*` files that import `slab_streaming` from
`orpheus.geometry` are **not gated** and will not go red when they start
importing from `sn/`. Legal, and worth a look: a *transport* test reaching into
`sn/` for a fixture is the layering smell the linter cannot see.

⭐ **And the thing that is NOT a blocker, contrary to how it reads:**
`AngularMeasure`'s retirement is suspended into CS5, which sounds like it
strands the file. `[M]` it has **5 use sites, all inside `reduced_operator.py`
itself** (the three factory signatures plus `angular_redistribution`), and
**zero** users anywhere else. So it travels with the factories and does not
keep the module alive. CS5 decides whether it should EXIST; P4.4 decides where
it LIVES — separable, and moving it prejudges nothing.

⟹ **Ruling: reorder, do not fuse.** §6b says fuse when a step boundary cuts a
call-site set; here it does not — the interval compiles and violates no edge.
Two smaller gated steps beat one large one, and P4.4 is already moving ~8
symbols and re-pointing 11 files.

⚠ **Why this survived being written, and it is worth more than the fix:** the
row says *"no new edge **either way**"* — which is **§6d's own vocabulary**. It
reads as a §6d check that was run and passed. It was not run; the phrase was
borrowed. See the surprise-log row of the same date.

**P4.4 — the residue to `sn/`, and the module is deleted.** Safe only after
P4.1a. `[M]` 0 intra-geometry consumers, so nothing stays behind.

✅ **LANDED `16501ca0`** 2026-08-28. ⛔ The row's second half is void as
written: *"the module is deleted"* moved to **P4.3** (the ⛔ block above says
why — `StreamingTerms` is still defined there and `transport` imports it at
runtime), and *"nothing stays behind"* is false at this phase — `AngularMeasure`,
`StreamingTerms` and the five α symbols all stay, by design.

**What moved** (4 symbols → `orpheus/sn/mesh/reduced_operator.py`, ✅ user-ruled
destination): `ReducedStreamingOperator` + `slab_streaming` /
`spherical_streaming` / `cylindrical_streaming`. ⭐ **Beside `SNMesh`, which
`[M]` is their ONLY constructor and holder** (`augmented_mesh.py:333/:340/:357`);
`sn/mesh/__init__.py`'s own docstring already described this object ("the
precomputed streaming stencil, the α / geometry-factor / Morel–Montry weights").
⚠ **The filename is deliberately NOT renamed** — `ReducedStreamingOperator` is
not a streaming operator, but §9.1's name fork is OPEN, so P3c renames class and
module together rather than this step minting a name it does not own
(`plan-authoring` §1).

**§6b call-site set — complete by TWO filters, and the second one found what the
first could not.**

| surface | method | count |
|---|---|---|
| production imports | AST | **3** — `geometry/__init__.py` (re-export **DELETED**), `sn/mesh/augmented_mesh.py`, `sn/operators/radial_characteristic.py` (`TYPE_CHECKING`) |
| test imports | AST | **13** — 8 direct, 5 via the `orpheus.geometry` façade |
| `docs/*.rst` refs | regex + positive control | **28** in 5 files, incl. the V&V `:by:` directive at `index.rst:625` |
| **`.py` docstring xrefs** | **`dead_references`** | **9 in 5 files** ⭐ |

⭐⭐ **The docstring row is the finding.** The import audit was complete and the
residual check returned 0; the full fast suite was green and `sphinx -W` built
clean — and **5 targets / 9 sites** still named `orpheus.geometry.reduced_operator.<mover>`
inside `:class:`/`:func:`/`:attr:` roles in `orpheus/sn/angular/closure.py`,
`transport/spatial/{scheme,cell_balance}.py` and two test modules. Nothing else
in the toolchain can see them: the module is **not** `automodule`'d, so Sphinx
renders none of it and `-W` is silent at every severity (`coding-standards`'
silent class). `mcp__nexus__dead_references` reported exactly 9; the sweep fixed
exactly 9; re-run reads **0 dead / 52 checked / 52 rescued**. ⟹ *a symbol has
more than one spelling SURFACE, and "the import grep is complete" is a claim
about one of them.*

**Test migration** (`coding-standards`: retiring means the tests follow the SUT).
`[M]` by AST the file split cleanly — **29 of 36** test functions name only the
movers, **7** (`TestAlphaDomeClosureContract`) only the residue, and all three
mesh fixtures belong to the movers. So
`tests/geometry/test_reduced_operator.py` → `tests/sn/mesh/test_reduced_operator.py`
for the 29, with the 7 staying (they travel again to `tests/sn/angular/` at
P4.2). `[M]` **63 collected → 18 + 45 = 63**, conserved — and confirmed
independently by `docs/theory/verification/matrix.rst`, which the build
regenerates and which I did not write.

**Two docs claims this move REFUTED, corrected in place (§3):**
* `structured_geometry.rst` — *"lifts the math into a **geometry-layer**
  primitive rather than a solver-side one"*. The single-sourcing half is
  Cardinal Rule 2 and stands; the **layer** half is precisely what P4.4
  disproved. Kept, banner-dated, with the layer test that replaced it.
* the same page's file **enumeration** — re-measured at **15** files, and `[M]`
  `transport/spatial/linear_discontinuous.py` was missing from the published
  list **already at HEAD** (it named `StreamingTerms` there too). Found only by
  re-running the page's own predicate rather than reading it — `plan-authoring`
  §2, a published enumeration owes its own re-run.

**Gates.** layer-import **342 → 343** (+1, predicted before running: the linter
`rglob`s every module in `orpheus/` and there is one more — the same mechanism
as P2's `sn/angular/__init__.py`); affected trees 1669 passed / 0 failed; the
pre-migration sweep/transport/operators run 2342 passed / 0 failed; pyright
**0**; `sphinx -W` **0**; `dead_references` **0**. Fresh-interpreter import
smoke passes with `orpheus.geometry` FIRST — the order that breaks under a
cycle, and the one a naive smoke test gets wrong.

**P4.5 — the hollow Cartesian object dies.** `[M]` d≥2 Cartesian gets
`reduced = None`; the slab gets an object whose every meaningful field is
`None`/neutral, minted *"so `sn_mesh.reduced` is always populated"* — a promise
the d≥2 arm already breaks. `reduced is None` today conflates **"no chain scan"**
(`ndim ≥ 2`) with **"no curvature"** (`coord is CARTESIAN`), and `SNMesh` answers
both directly. `[M]` **12 of 36** code reads of `.reduced` are `None`-guards
paying for the conflation. The object should exist **iff there is a radial
reduction**.

⭐ **Two sites P4.1a surfaced that belong to THIS step — recorded here so they
earn a phase rather than an issue-and-forget.** Both are the same defect: a
predicate spelled defensively because the object it interrogates does not say
what it is.

* `sn/sweep/cache.py:269` — `assert sn_mesh.reduced is not None, "…requires a
  ReducedStreamingOperator (1-D Cartesian / spherical / cylindrical). 2-D
  Cartesian wavefront uses anti-diagonal scheduling, not the chain scan."`
  This is an **admission contract**, not type-narrowing, and it is a bare
  `assert` — so `[M]` the canonical `-O` runner **strips it** and the cache
  accepts a 2-D Cartesian mesh silently (`coding-standards`: a bare `assert` in
  `orpheus/` is not a contract). P4.1a left it untouched deliberately: its
  subject is exactly the `reduced is None` conflation this step re-poses, so
  converting it now would mean writing the guard twice. ⟹ **P4.5 owes it a real
  `raise` keyed on the honest predicate** ("does this mesh have a chain scan?"),
  not on `reduced is None`.
* `sn/solver.py:2316-2324` — `_is_curvilinear` reads
  `getattr(mesh, "coord", None)` then `getattr(coord, "name", str(coord))`
  *after* an `isinstance(mesh, Mesh1D)` narrowing. `coord` is a required
  dataclass field on `Mesh1D`, so **both defaults are unreachable**, and the
  function then compares an upper-cased *string* against `("SPHERICAL",
  "CYLINDRICAL")` rather than comparing the enum. It is the stringly-typed
  spelling of `mesh.coord is not CoordSystem.CARTESIAN`, which `SNMesh` already
  answers as `is_cartesian`. ⚠ Not a P4.1a site (it reads the *mesh*, which is
  the successor), and not a bug today — but it is a third spelling of this
  step's own predicate, so it retires with the other two.

**P4.6 — the moment mass consumes the MEASURE, not a chart tag.**
✅ RULED 2026-08-27 (user): *"it is very important to correct the flaw"*, in the
move, after compaction.

`[M]` shipped `M` does **not consume the measure at all**. `LinearDiscontinuous.
moment_mass_diagonal` builds it from `mass_1d(np.ones(()), self.theta)` — the
1-D mass factor at **UNIT WIDTH** — Kronecker'd over `ndim`. That is exact on a
slab (where `M/V` is width-independent, which is *why* the shortcut works) and
wrong on a curved chart. `_assert_moment_mass_is_expressible` keeps it honest by
**refusing** curvilinear multi-moment rather than by computing.

⭐ **This is the SAME fix as P4 item 1** (retire the TRANSITIONAL
`coord: CoordSystem` tag), not a second one: the tag is a **proxy for the
measure**. Its own docstring already says so — the predicate is
`is not CARTESIAN`, i.e. *"is the volume element constant across the cell?"*,
a property of the MEASURE's behaviour. Once `M` is handed the measure, the tag
has nothing left to stand for and both retire together.

⚠ **The refusal does not vanish — it moves to an honest predicate**, and the
three arms are not equally affected:

* **single-moment (DD), any chart** — `∫b₀b₀ dV / V = V/V = 1`. Unaffected on
  every chart, by construction.
* **multi-moment on Cartesian** — `[1, θ]`. Already correct.
* **multi-moment on curvilinear** — the true `M/V` is **non-diagonal AND
  cell-dependent** (`[[1, 0.5], [0.5, 0.4]]` at a spherical pole cell), and
  `FunctionSpace`'s metric is Hadamard **by construction**, so the VALUE is
  currently unspellable — ORPHEUS **#409**.

⟹ the guard's predicate changes from *"is the chart Cartesian?"* (a proxy) to
*"can this metric be expressed?"* (the question). That is strictly better and is
a **step toward #409 rather than a step blocked on it**: the honest refusal names
the machinery gap instead of a chart, so the day #409 lands the guard retires by
itself. ⚠ Do NOT read this ruling as "compute the curvilinear multi-moment mass"
— that value needs #409's non-Hadamard metric, and #158 to give it a consumer.

**P4.7 — naming, the fusion, and ⭐ the gap both of them point at.**
✅ RULED 2026-08-27 (user).

**(a) The fusion.** *"A fusion is accepted to be cached as a performance
optimization at the place that will assemble a hot path … we're not writing
algebra to conform to performance optimization and welding things."* `[M]` the
two fusion sites are NOT alike:
* `StreamingCoefficientCache.dA_w` — an `(N, nx)` array built **once** at solver
  init. That IS pre-operating at the hot-path assembly point. **KEEP**, renamed.
* `StreamingTerms.delta_A_over_w` — built **per (cell, direction)** inside the
  loop, so it buys no performance; its inputs are recomputed as often.
  ⟹ **RETIRE.**
⛔ And a fourth present-tense-false docstring: `_weight_of` says *"Callers form
`ΔA/w` … rather than reading a stored product that fused a geometric with a
quadrature quantity (Pattern 2 — neither side owned the fusion)"* while the
packet stores exactly that product. The fusion moved one frame out; it was never
eliminated.

**(b) Names — one greppable family.** *"this is not Fortran punch-cards anymore
… It's more important for the name to be self-explanatory."*

| current | `[M]` what it IS | new |
|---|---|---|
| `dA_w` | | `delta_A_over_w` |
| `V` | | `volume` |
| `A_down` | `v.face_area_downstream` — the **sweep-direction** downstream face | `face_area_downstream` |
| `A_total` | `face_area_inner + face_area_outer` (`cache.py:329-330`) — the **sum** | `face_area_total` |
| `face_area_inner`/`_outer` | the two faces by **chart position** | unchanged |

⭐ `A_down` needs **no coinage — the tree already has the name.**
`CellVisit.face_area_downstream` is used in `cell_balance.py` and `diamond.py`,
and `cache.py:324` literally reads `v.face_area_downstream` then stores it under
a shorter one. A Pattern-2 unification, not a rename. ⟹ the family greps as
`face_area_`: `_inner`, `_outer`, `_downstream`, `_total`.
⚠ They are **three distinct concepts, not two** — chart-frame position,
sweep-frame direction, and their sum.

**(c) ⭐⭐ The gap the producer and `AngularMeasure` both point at.**
User, on the producer holding mesh + quadrature: *"it makes it seem like it
needs the spatial axis of space and the angular axis of space … if so, then it's
just space and R, not a stage-2 generator."*

`[M]` **Half true, and the failing half is the finding.**

| the producer reads | does a space carry it? |
|---|---|
| `volumes`, `widths`, `weights` | ✅ `Axis.weights` |
| `mu_x`, `eta`, `mu_z`, `level_indices` | ⛔ **0 hits** in BOTH `numerics/space.py` and `numerics/axis.py` |

Because `Axis` is a **lossy projection of `DiscreteMeasure`**:

```
DiscreteMeasure : nodes, weights, support, invariance_group, exactness
Axis            : label, shape,          weights,            kind
```

— it keeps the weights and **drops the nodes**, and the direction cosines ARE
the nodes. `[M]` the exact loss site is `augmented_mesh.py:1170-1175`, which
builds the angular axis from `quad.N` + `quad.weights` alone and declares
`kind=BasisKind.NODAL` — **a nodal basis carrying no nodes.**

⟹ **the producer IS `(space, R)`; the space is currently too thin to say so.**
Give the angular axis its measure and the producer becomes an operator bound to
one space — which is this campaign's own root ruling ("an operator is not an
operator without its two spaces").

⛔ **This SUSPENDS `AngularMeasure`'s retirement** (P4 item 3). That ruling
argued CONSUMERS (*"the 3 factories read 0 of its 6 members"*); the question
here is TYPE DESIGN — should `DiscreteMeasure` branch into
`SpatialMeasure` / `AngularMeasure` / `EnergyMeasure`, so an axis can carry a
real measure? A consumer count cannot refute a type design. `AngularMeasure` is
at `orpheus/geometry/reduced_operator.py:304-360` (a structural `Protocol`:
`mu_x`, `weights`, `N`, `eta`, `mu_z`, `level_indices`); its implementer
`Quadrature` **wraps** a `DiscreteMeasure`, and the family already branches
(`BundleMeasure`, `DiscreteMeasurePartition`, `GeneratingMeasure`,
`ReferenceMeasure`, `UniformMeasure`, `ProductMeasure`). **Investigate before
retiring.**

**P4.8 — ⭐⭐ an axis is a FORGETFUL MAP and should carry its generator.**
User hypothesis 2026-08-27: *"The space axes are forgetful maps. They get a
Discrete Measure and get the data they need from the Measure, but they don't own
the measure. They should hold an accessor to the Measure that generated them as
provenance data."* `[M]` **It holds — with one refinement the code already
supplies.**

**(a) The forgetting is real, and located.** `Axis` is
`@dataclass(frozen=True, eq=False)` over `label, shape, weights, kind`, against
`DiscreteMeasure`'s `nodes, weights, support, invariance_group, exactness`. It
**keeps the weights and drops the nodes** — and the direction cosines ARE the
nodes. Loss site: `augmented_mesh.py:1170-1175`, which builds the angular axis
from `quad.N` + `quad.weights` alone while declaring `kind=BasisKind.NODAL`.
`[M]` `mu_x`/`eta`/`mu_z`/`level_indices` → **0 hits** in BOTH
`numerics/space.py` and `numerics/axis.py`.

**(b) ⭐ But a MODAL axis is generated by a BASIS, not a measure** — so the rule
generalises rather than excepting. `[M]` `FrameBase`'s fields are literally
**`['basis', 'measure']`**, with `analysis` (nodal→modal), `reconstruction`
(modal→nodal), `discrete_gram`, `basis_space` and `measure_space`. And the one
shipped MODAL axis confirms it: `scheme.py:1517` sets
`weights=self.moment_mass_diagonal(...)` — i.e. **`diag(Gram)` of its
tensor-Legendre basis.**

| axis kind | generator | its `weights` are |
|---|---|---|
| **NODAL** | a **measure** (nodes + weights) | the quadrature weights |
| **MODAL** | a **basis** (functions) | `diag(Gram)` of that basis |

⟹ **the general statement: every axis is a forgetful map from its GENERATOR;
`kind` already names WHICH KIND of generator, and the accessor would name the
generator itself.** That is verbatim what `Axis`'s own docstring claims its
identity is — *"the identity is **what kind of generator produced this factor**,
not a bag of fields"* — a claim the fields cannot currently honour. The frame is
the object holding BOTH sides, so a frame-minted pair is the two-sided case.

**(c) ⛔ CORRECTION — energy is NODAL, not MODAL over an indicator basis.** A
session recollection had it the other way. `[M]` `EnergyAxis.__post_init__`
(`axis.py:217-362`) **refuses** anything else, verbatim:

> *"an EnergyAxis is NODAL by construction: groups are the CELLS of a 1-D energy
> mesh (components are per-group integrals, not spectral coefficients)"*

and separately refuses weights at all: *"the multigroup convention makes the
energy measure COUNTING as a theorem (group integrals × group averages pair
without widths)."*

⭐⭐ **The misunderstanding is worth naming precisely, because the ENUM'S NAMES
INVITE IT.** An indicator expansion IS "modal" in the classical FEM sense —
coefficients, not point values. But `BasisKind` does not track that distinction.
Its own docstring gives the predicate it *does* track:

> `NODAL` — *"components are point/cell VALUES (**indicator-like basis**): the
> factor carries a **coordinate cone**, so per-component positivity is a
> meaningful predicate."* `MODAL` — *"expansion COEFFICIENTS (spectral basis):
> no coordinate cone; a positive function may have negative coefficients."*

⟹ **`BasisKind` answers *"does the coefficient cone coincide with the function
cone?"*, not *"are these point values or coefficients?"*** A per-group integral
of a non-negative function is non-negative, so indicators are cone-preserving
and classify NODAL. The enum exists to answer `has_coordinate_cone` — this
campaign's own cone ruling — and its labels are classical-FEM-flavoured for a
predicate that is about the cone. ⚠ A naming item, not a defect: the values are
correct under the stated predicate.

**(d) Design consequences.**
* `[M]` **No import cycle**: `numerics/axis.py` imports `measure` **0**;
  `numerics/measure.py` imports `axis` **0**. New edge, no reverse — §6d clean.
* **Exclude provenance from `_identity_key`.** Two axes with identical
  `label/shape/weights/kind` are the same axis whatever instance produced them;
  provenance is not identity. This also sidesteps ORPHEUS **#403**
  (content-equal `DiscreteMeasure` instances make `frame ==` RAISE).
* `[M]` **3 of 4 shipped axis sites are NODAL** — including both spatial ones —
  so this is not an angular-only fix. A mesh is a discrete measure in principle
  (nodes at cell centres, weights the volumes).

**(e) ⟹ ⛔ THE WORK HAS A HOME, AND P4'S REMAINDER IS BLOCKED ON IT — a named
blocker, NOT a defer** (user's standing ruling: deferred work earns a phase
number or a named external blocker). The machinery is **space-layer
architecture**, not names/homes/welds, so by this plan's own §5b precedent it
belongs to the space campaign: **`space_and_kernel_binding_campaign.md` §5.5,
Phase CS5 — "an axis can name the generator that made it"**, chartered
2026-08-27. It carries the three machinery items, the sequencing (NODAL half
first), the done-when, and `AngularMeasure`'s suspended fate.

⟹ **Blocked on CS5's NODAL half:** P4's producer home, and the `ChartConnection`
name. ⟹ **NOT blocked:** P4.1a/1b, P4.2, P4.3, P4.4, P4.5, P4.6, P4.7 and P4b —
all of which can land first.

**This is the ROOT of the remaining cluster, and smaller than what depends
on it.** Give NODAL axes their measure → the angular axis can answer `mu_x`,
`eta`, `mu_z`, `level_indices` → the producer binds to `(space, R)` and becomes
an operator with its two spaces → **which is this campaign's root ruling**, and
the `ChartConnection` name resolves as a consequence rather than a decision.
`AngularMeasure`'s fate follows too: not a retirement but *which measure type* —
the `DiscreteMeasure` → `Spatial/Angular/Energy` branching question.

### ⛔ The four options this supersedes — kept per §3, with their refutations

The P4 fork (§7) offered four. All four were answers to the wrong question, and
three are independently refuted:

* **(a) `ChartConnection` leaves `geometry/`** — reframed, not refuted. It does
  leave; but as *posing*, and only after P4.1a, and `StreamingTerms` splits off
  to L2 rather than travelling with it.
* **(b) the α type to `numerics/`** — ⛔ REFUTED **on the concept**, though it
  passes the layer contract cleanly (`numerics → geometry` is permitted and
  already exists 2 RT + 1 TC; import-probed under 5 entry orders, no cycle).
  `alpha_dome` is not a shareable primitive: **one consumer, identity morphism**,
  and the one other family that would want this — DG/FE-in-angle — **cannot call
  it**, needing `f` at its own edges and in the interior. The honest shared
  object is `f = r·dx/ds`, a different thing. `(ndarray, ndarray)` looks
  universal, but *the arrays are a quadrature* — the object every non-consumer
  lacks by construction (`lessons.md` L60).
* **(c) `ChartConnection` stops NAMING the angular type** — ⛔ REFUTED. `[M]` the
  class reads the field's interface at **5 sites** (`:658` `weights`, `:707`/`:725`
  `mu_x`, `:755` `level_indices`, `:757` `eta`). ⭐ But every one dereferences
  `.quadrature` first, and **none touches `alpha_per_level`, `mu_start_per_level`,
  `coord` or `n_levels`** — geometry never reads the dome at all. That is the
  weld: the class needs the MEASURE and is handed the REDISTRIBUTION.
* **(d) a `WHITELIST` entry** — ⛔ DEAD TWICE. A whitelist silences the linter; it
  cannot silence `test_entry_point_imports_in_a_fresh_interpreter`, which
  cold-imports `orpheus.geometry` in a subprocess. And `[M]` the whitelist has no
  staleness gate — 1 of its 4 entries already covers **0** imports (ORPHEUS #411).

### Open forks — RULE BEFORE EXECUTING

1. ⛔ **§9.1's `ChartConnection` name rests on a residue this decomposition
   changes.** It was ruled because *"the residue is `face_areas` + `delta_A` …
   connection data, not metric data"*. Under P4.1a/1b `face_areas` is retired and
   `delta_A` is derived, so the residue is a **producer of sweep coefficients**,
   not connection data. Re-rule the name. ⚠ `plan-authoring` §1: grep the PROSE
   corpus (`lessons.md`, `plans/`, `agent-memory/`, closed issues) before adopting
   any candidate — the P3b precedent is that a free name can be free *because it
   was rejected*.
2. ✅ **PARTLY RULED 2026-08-27 (user) — `redistribution_pairing` (R) goes to
   `transport/` (L2).** The reason is its AUDIENCE, established by separating R
   from `M` by consumer:

   | | consumes | audience |
   |---|---|---|
   | `M = ∫ b_k b_j dV` | basis + volume measure | broad — `[M]` **8** references, ALL in `transport/`, **0** in `sn/`; plus LS-MoC, FE/DG diffusion, nodal, CP-with-linear-source in principle |
   | `R = ∫ b_k b_j (∇·ê_r) dV` | basis + volume measure **+ the connection** | `{SN actual, Pn potential}` |

   They differ by **exactly one ingredient**, and *neither consumes anything
   angular* — no `μ`, `w`, `τ` or `α`. That is what makes them the spatial
   factor.

   ⭐ **R is NOT SN-specific, and this plan said it was.** `[R]` reasoned (there
   is no `pn/` package): a method consumes R iff it has a `1/r` term multiplying
   the unknown AND projects it onto a spatial moment basis. The curvilinear **Pn**
   moment equations carry `c_ℓ·φ_{ℓ±1}/r`; projecting onto `b_k` with
   `φ = Σ_j c_j b_j` gives `Σ_j c_j ∫ b_k b_j (1/r) dV` — that IS R. SPn and
   FE/DG diffusion do NOT: their `1/r` lives inside `∇²`, so the weak form yields
   the **stiffness** matrix `∫∇b_k·∇b_j dV`, a different bilinear form. MoC/CP/MC
   have no `1/r` term at all.

   ⭐⭐ **So α and R split, and this plan had treated them as one placement
   question.** Pn needs **no α** — its angular derivative folds into the Legendre
   recurrence, so no angular differencing coefficient is ever formed — but it
   does need R. ⟹ **α is SN-specific; R is specific to methods that discretize
   an angular variable in a ROTATING FRAME, by collocation (SN) or expansion
   (Pn).**

   ⚠ **STILL OPEN: where the PRODUCER goes.** `[M]` `streaming_terms()` reads
   `angular.quadrature.{mu_x, eta, level_indices, weights}` (angular),
   `mesh.{volumes, widths}` + `face_areas` (geometry), `delta_A` (the
   connection) and `coord`. **It consumes BOTH factors — it is the fusion
   site**, so the R ruling does not place it. Its inputs are all L2/L1/INPUT, so
   `sn/` is legal; it produces an SN sweep packet, which argues it is SN's. The
   packet TYPE `StreamingTerms` must be ≤ L2 regardless (`transport/spatial/`
   imports it at runtime and `transport → sn` is forbidden).
3. ✅ **RULED 2026-08-27 (user) — `delta_A` goes to `sn/`, NOT to `coord.py`.**
   Verbatim reason, and it is the durable half: *"`coord.py` should have geometry
   inside, method agnostic, serves all methods equally. `delta_A` in `coord.py`
   means a concept from curvilinear SN leaked into geometry."*

   ⭐ **This supersedes the weaker argument this plan first offered, and the
   difference is load-bearing.** The proposal reasoned from CONSUMER COUNT — "it
   has exactly one, and the project defers a shared primitive until two"
   (`coding-elegance` Pattern 6). That argument **flips** the moment a second
   SN-ish consumer appears. The ruling reasons from `coord.py`'s **CONTRACT**:
   method-agnosticism is what the module is *for*, so an SN-specific concept is a
   leak at ANY consumer count. It does not flip. When a placement question has
   both an anticipation argument and a contract argument, the contract argument
   is the one to record — the other is a guess about the future wearing a rule.

   `[M]` 2026-08-27, and it upgrades the ruling from a preference to an invariant
   the module already holds: `coord.py`'s entire public surface is `CoordSystem`,
   `compute_volumes_1d`, `compute_areas_1d`, `compute_volumes_2d` — and a regex
   census over seven method vocabularies (`sn`/`S_N`/discrete-ordinate, `moc`/
   characteristic, collision-probability, monte-carlo, `diffusion`, `sweep`,
   ordinate/angular/quadrature), positive control `CoordSystem` present, returns
   **0 hits in every one**. `delta_A` would have been the first violation.

   ⟹ **the standard, stated so a later phase can check it:** nothing enters
   `orpheus/geometry/coord.py` that names or presupposes a solution method.
   A cheap tell is the census above — if it stops reading 0/7, something leaked.

   ✅ **RULED 2026-08-27 (user): `delta_A` DISSOLVES into R.** ⛔ ~~Not recorded
   as ruled: it is a CONSEQUENCE derived here, not a decision the user stated.~~

   R's rank-1 entry **is** `ΔA`, so `delta_A` and R are not two objects;
   `delta_A` is R's rank-1 realization. If R lived in `transport/` and
   `delta_A` in `sn/`, one quantity would span two layers and R's L2 producer
   would need an L3 import (**forbidden**). So `delta_A` does not travel to
   `sn/` as a named field at all — **R's producer derives the connection
   integral from `mesh.areas` itself**, which honours this fork's contract
   argument (nothing SN-flavoured enters `coord.py`) while removing the layer
   conflict.

   ⭐ **This SUPERSEDES this fork's DESTINATION while keeping its REASON.** The
   ruling above ("`delta_A` goes to `sn/`") answered *"where does the field
   live?"*; the answer is now *"it stops being a field."* The contract argument
   that produced it — `coord.py` is method-agnostic, so an SN-flavoured concept
   is a leak at any consumer count — is untouched and still binds.

   ⟹ **P4.1b is the first step of this dissolution**, not merely adjacent to
   it: making `delta_A` a `cached_property` over `np.diff(mesh.areas)` is
   already *"derive the connection integral from `mesh.areas`"* — it stops
   being stored data one phase early, so what P4.4 later moves is a derivation,
   not an array.
4. **Freezing.** `redistribution_pairing`'s docstring defers `frozen=True` to
   "the re-home, where the class definition is being touched anyway" — because
   it synthesises `__hash__` and the fields hold ndarrays. `eq=False`? per-field
   `compare=False`? Comes due in P4.4.

---

## 5. Where `L` and the ray-characteristic operators MEET — the scheme's mint

**The question the user posed:** if two objects share something, and the shared
thing is discretization-related, that is where the shared concept should be
constructed.

`[M]` the sharers and the shared: `StreamingOperator` (`L`) and the
radial-characteristic family (`A_BB`, `Emission`, `Reconstruction`, `Seeding`)
**both consume the angular closure** — `L` through `cell_contribution` /
`angular_adjoint` in the matvec, the RC family at
`radial_characteristic.py:886/:894/:961`. They also both need the pairing,
transitively, because the closure carries it.

⟹ **The meeting point is the scheme's minted package** (R14/R19: the scheme is
a stage-2 generator on the Frame pattern), and the shape follows from F1 —
the scheme owns the basis, the chart owns the measure, the quadrature owns the
angular side:

```
package = scheme.mint(chart)              # S1: needs the mesh's chart
    .spatial_closure                       #   the cell relation (DD / LD)
    .redistribution_pairing                #   R  — (nx, n_mom, n_thread)
    .moment_measure                        #   M/V for the space's moment axis  ⚠ §6

angular_algebra = angular_redistribution(quad, coord)      # S0: mesh-free
angular_closure = device(angular_algebra,                  # device ∈ {WDD/M-M,
                         package.redistribution_pairing)   #   plain diamond, angular-LD}
```

and **both** `L` and the RC operators receive the minted `angular_closure` —
neither builds it, and neither reaches for a mesh to get it.

⭐ **Why the device is a separate axis from the scheme** (and why the scheme
must not simply *be* the angular closure's factory): `[M]` the rank
contradiction. Adams–Martin close the angular index **per spatial moment**
(square `R`); ONETRAN closes it on the cell **average only** (rank-1 column).
Those are two different *angular devices* over the same spatial scheme. Folding
the device into the scheme would make that choice unspellable — the very error
the `(n_mom, n_thread)` axes exist to prevent.

**Where the closure's own split lands (F3):** `angular_algebra` is S0 and
shareable across meshes; the `device(…)` binding is S1 and rebuilds with the
chart. That is the function/evaluated-table split, on the angular axis.

---

## 5b. ⭐ The three factors — make the algebra the implementation

**Added 2026-08-26 after a user question this plan had NOT anticipated.** The
adjoint audit's closing observation was *"signature of a product of one diagonal
factor and two triangular ones"*, and the question is whether those three become
first-class. This section says yes, states what must be verified first, and is
deliberately scoped as a SUCCESSOR campaign, not a phase of §7.

### The measured state — they are not explicit

`[M]` 2026-08-26:

* the **composition machinery already exists**: `OperatorProduct`
  (`numerics/operator.py:1650`) composes and is a `LinearOperator`, so `.H`
  propagates through it (`:983`); `PermutationOperator` (`:2458`),
  `ScaledOperator`, `PointwiseOperator`, `TensorProductOperator` and
  `_AdjointOperator` all ship;
* there is **no triangular-factor concept anywhere** — `grep` for
  `Triangular` / `forward_substitution` over `orpheus/` returns **one** hit,
  and it is `OperatorProduct`'s unrelated line;
* and the SN loss matvec is a **fused body** in `sn/loss_representation/`, not
  a composition — so today `apply` and `apply_transpose` are **two hand-written
  implementations of one relation**.

### Why this is the highest-value item the audits surfaced

`[M]` `apply_transpose ≡ apply.T` to `3.05e-16` / `1.99e-16` / `0.0`
(sphere / cylinder / slab). That gate exists **because they are twins** — it is
Pattern 2's signature, and the project rule is to make the twin unspellable
rather than to test that the two copies agree.

If the factorization is real and declared:

```
L  =  D  @  T_ang  @  T_spatial
L.H  =  T_spatial.H  @  T_ang.H  @  D          # D diagonal ⇒ self-adjoint
```

then **the adjoint stops being written and starts being derived.** The audit's
own inventory says exactly what each factor's `.H` must do, and it is small:
`[M]` every coefficient enters the adjoint **at the same value**, and the
reversal needs only **two order reversals and one block swap**. Those are
`T_ang.H` (the reversed march), `T_spatial.H` (the reversed scan), and the
trace-block swap. ⟹ `angular_adjoint` and the hand-written reverse scan RETIRE,
and the seam where a future angular device (ONETRAN's rank-1 column,
angular-LD) plugs in becomes ONE factor instead of two hand-mirrored bodies.

This is the *"build the machinery; realize the operator algebra ALWAYS"* ruling
applied to the last welded operation in the streaming path: a welded, un-named
composition is a failure to realize the algebra.

### ⛔ P0 — VERIFY THE FACTORIZATION BEFORE DESIGNING TO IT

⚠ **The three-factor claim is `[R]`, not `[M]`.** It was *inferred* from the
symmetry characters of `∂L/∂k` (`0` skew / `2` symmetric / `1.414`
triangular-like), which is evidence of triangular STRUCTURE and **not** a proof
that `L` equals a product of exactly three named factors. Designing to an
unverified factorization is how a campaign inherits a false premise.

**Done when** — all three answered, on the sphere (⛔ never the cylinder: `[M]`
`R/ΔA` is bit-exactly `diag(1, 1/3)` there and cannot discriminate):

1. `‖D·T_ang·T_spatial − L‖ / ‖L‖` at machine precision on ≥2 charts, with the
   factors constructed independently of `L`'s assembly;
2. **the factor COUNT is confirmed, not assumed** — the audit says "two order
   reversals AND one block swap", so the honest arity may be four
   (`D · P_mirror · T_ang · T_spatial`) once the pole-mirror permutation and
   the inflow/outflow trace partition are placed. Report the count you MEASURE;
3. each factor's triangularity is checked **in its own index order** (`T_ang`
   in the ordinate index at fixed cell; `T_spatial` in chain order), because a
   matrix is triangular only with respect to an ordering, and the ordering is
   itself a traversal artefact (§4, S3).

A refutation here is the most valuable outcome, not the least — it would mean
the adjoint is genuinely irreducible to a composition and the twin must stay,
which is a thing worth knowing before anyone tries.

### The constraint that shapes the design — algebra ≠ realization

⚠ The sweep is the **hot loop**, and `CumprodScan` (Blelloch) exists to make it
fast. Applying three factors sequentially through typed carriers would allocate
three times and destroy that fusion. **Do not read this section as "apply three
operators."**

The resolution is already this tree's own vocabulary, twice over: `WindowedSweep`
is *"the FUSED typed product `P @ A.inverse()`"* (a declared composition with a
fused evaluation), and §5c's ruled **function / evaluated-table** split says the
same thing. So:

* **declare** the factors — they are the algebra, and they carry `.H`;
* **realize** the composite's `apply` through the existing fused scan;
* and **gate** that the two agree.

⭐ That gate is the `3e-16` measurement, **repurposed**: today it compares two
hand-written twins (and so can only ever say "the copies match"); afterwards it
compares the DERIVED adjoint against the FUSED realization — an
algebra-vs-implementation equivalence, which is a strictly stronger claim from
the same number.

### Scope

This is **not** a phase of §7 — §7 is names, homes and welds, all bit-identical
or nearly so, and this is an architecture change to the operator algebra. It
belongs with the operator-strategy campaign / O-3's `StreamingOperator` binding
(`posing_filtration_charter.md` §5c), and it should follow §7 so the factors are
declared in objects that already have honest names and homes.

⚠ **Dependency, stated because it is easy to miss:** `T_ang`'s content is the
angular closure's S0 algebra (§4), so F3's function/evaluated-table split is a
PREREQUISITE, not a parallel. Splitting the closure after declaring `T_ang`
would cut the same seam twice.

---

## 6. The one CORRECTNESS item — not cosmetics

⛔ **A curvilinear cell is being given the slab's moment metric.** `[M]` an LD
spherical `SNMesh` **constructs today** and installs `spatial_moment weights =
[1, 0.3333]` (the slab mass) while the true `M/V` at the pole cell is
`[[1, 0.5], [0.5, 0.4]]`. It does not bite only because `loss_representation`
refuses one layer out (`IncompatibleRepresentation`, #158).

At `n_mom > 1` there are **three different matrices on one index** — `R/ΔA`,
the true `M/V`, and the shipped `diag(1, θ)` — and **only the third is a
metric.**

⚠ The hard part is not the value, it is that the true `M/V` is **cell-dependent
and non-diagonal**, which a per-axis `Axis(weights=…)` cannot express. `[M]` a
Gram–Schmidt basis diagonalises `M` and does **not** diagonalise `R`, so a
basis change cannot make both diagonal.

⟹ **This plan does not fix it** (it needs #158's cell solve to have a
consumer). It must **refuse honestly** instead: the metric installer should
reject a curvilinear multi-moment space rather than install the slab mass. §7
P1 carries that guard.

#### ⭐ Where the slab mass ENTERS — `[M]` located 2026-08-26

The installer is `SpatialScheme.moment_axis(ndim)` →
`Axis(…, weights=self.moment_mass_diagonal(ndim))`
(`transport/spatial/scheme.py:1375-1416`), and LD's override
(`linear_discontinuous.py:623-636`) builds it as

```python
base = np.diagonal(mass_1d(np.ones(()), self.theta))   # [1, θ]
```

⭐⭐ **`mass_1d` evaluated at UNIT WIDTH is the Cartesian-only shortcut, and
it is exactly where the slab assumption enters.** For a slab, `M/V` really is
width-independent — which is *why* the shortcut is valid and why nothing has
ever caught it. On a curvilinear chart `M/V` depends on `r_i`, so unit width
is not a normalisation, it is a different matrix. The docstring's claim
*"This IS `diag(M)/V`"* is true on exactly one chart family.

#### ✅ The §6c witness — CONSTRUCTED, not argued (`[M]` 2026-08-26)

§6c says a gate must name an input that exists in the tree *today* and that
it rejects. Built by hand rather than reasoned about:

```
slab      LD SNMesh: CONSTRUCTS. moment weights = [1.  0.33333333]
sphere    LD SNMesh: CONSTRUCTS. moment weights = [1.  0.33333333]
cylinder  LD SNMesh: CONSTRUCTS. moment weights = [1.  0.33333333]
```

⭐ **Wider than §6 stated** — §6 named only the sphere; the **cylinder**
constructs too, and all three charts install a **bit-identical** slab mass.
The witness is therefore
`SNMesh(Mesh1D(coord=SPHERICAL), gauss_legendre(4), …, scheme=LinearDiscontinuous())`
and its cylindrical sibling, and the gate must name one of them.

✅ **And the blast radius is much smaller than "a mesh stops constructing".**
`[M]` `moment_axis` is reached from a **`@cached_property`**
(`augmented_mesh.py:1227`, the trial space), **not** from `__init__`. So the
refusal does not make an LD curvilinear `SNMesh` unconstructible — it fires
when a consumer actually reaches for the moment metric, i.e. exactly where the
wrong value would otherwise be installed. The mesh is not the wrong object;
*asking it for a slab mass on a curved chart* is the wrong question.

⛔ **And the guard §7 P1 promises cannot be written where the plan assumed.**
`[M]` **neither** `moment_axis(self, ndim: int)` **nor**
`moment_mass_diagonal(self, ndim: int)` receives a chart — so the installer
cannot refuse a curvilinear space *even in principle*: it is never told which
one it is on. ⟹ P1's item is **not** a guard insertion; it is a signature
change (the installer must be told the chart) plus the refusal, i.e. a §6b
call-site change. Sized accordingly in P1.
⛔ **CORRECTED 2026-08-26 — the first count here was mine and was wrong two
ways.** It read *"8 references total, 5 in `orpheus/` … and 5 in `tests/`"*:
`5 + 5 = 10 ≠ 8`, so it did not even add up, and both halves were built from
a partial grep that mixed defs and prose into a "call-site" count.

`[M]` re-derived by grepping the CALL form (`\.moment_axis\(` /
`\.moment_mass_diagonal\(`), which is the predicate that matters:
**13 call sites — 3 production, 10 test.**
Production: `transport/fields/_bases.py:266`, `transport/spatial/scheme.py:1414`
(internal, inside `moment_axis`), `sn/mesh/augmented_mesh.py:1227`.
Test: 10, concentrated in `tests/sn/mesh/test_angular_bulk_space.py` (8) and
`tests/numerics/test_spatial_moment_field_space.py` (2).

✅ **And the signature question is settled by a measurement, not a taste
call:** `[M]` **both** production callers already hold a mesh that carries the
chart — `MaterialMesh.coord` (`material_mesh.py:230`) and `SNMesh.coord`. So
passing the chart costs the callers nothing; it only stops the callee
*assuming* it. `moment_mass_diagonal` computes `∫ b_k b_j dV`, which needs the
basis (the scheme's) AND the measure (the chart's) — taking only `ndim` is
precisely the missing half.

⭐ **Two expressiveness failures, not one** — worth stating separately because
they need different machinery, and #409 must cover both:
1. the true metric is **non-diagonal within the moment axis** ⟹ a matrix
   metric (#409);
2. it is **cell-dependent**, i.e. it varies along a *different* axis ⟹ it is
   not a per-axis weight at all, but a coupling between the cell axis and the
   moment axis (block-diagonal in cell, dense in moment).
   `Axis(weights=…)` fails **twice over**. #409's operator-valued framing does
   cover both — the metric is an operator on the product space — but a
   "matrix metric on the moment axis" framing would cover only (1).

⛔⛔ **And the instrument that cannot adjudicate it:** `[M]` reciprocity
(`A† ≡ G⁻¹AᵀG`) is an **identity for every invertible `G`** — `1.4e-16` under
Euclidean, random, and `(V·w)³`; mismatch control `8.22`. A reciprocity gate
proves *loadedness and consistency*, **never choice of metric**. What pins `G`
is the physical functionals (scalar flux, reaction rate). Do not cite a
reciprocity gate as evidence that a metric is right — here or anywhere.

---

## 7. Phasing — bit-identical first, behaviour last

Ordering principle: every phase that can be bit-identical **is**, and says so in
its done-when, so a value move is never hidden inside a move-and-rename.

### P1 — the dead stops being carried *(bit-identical; no renames, no moves)*

> ## ✅ P1 CODE COMPLETE — 2026-08-26, branch `refactor/p1-carry-only-what-is-read`
>
> 6 commits, 27 files, +439/−183. `npx pyright orpheus/` **0 errors**.
> ⏳ The full canonical gate (`-m "not slow"`, serial, whole tree) is the
> exit check — read `scratch/_p1_exit_gate.log`'s summary before merging.
>
> | commit | item |
> |---|---|
> | `37d6d1af` | 1–2 · `is_cartesian` → `coord`; both dead flags retired |
> | `fc8b1a0a` | 3 · the `_quadrature` twin, dissolving an `Optional` |
> | `ebe5d22f` | 4 · the three-link `mu_start` chain |
> | `983b36f9` | 5–7 · pairing identity; two lying docstrings |
> | `bc1fb804` | 8 · `curvature` retired |
> | `500de1b4` | 9 · the moment mass is told its chart |
>
> ⭐ **What the phase actually bought, beyond deletions.** Three of the nine
> items dissolved an `Optional` or a guard rather than merely removing a
> field: `_quadrature`'s retirement made a `None` state unrepresentable and
> took its narrowing `assert` with it; `curvature`'s removed a *runtime
> re-validation of an enum's own domain*; item 9 replaced a silent wrong
> value with a refusal that names both of its blockers. Pattern 4 three
> times over — the un-weld is what made each one spellable.
>
> ⛔ **The one defect P1 introduced, and its rule.** Item 8's residual grep
> (`\.curvature\b|curvature\s*[:=]`) could not see
> `getattr(sn_mesh, "curvature", None)` — a name inside a STRING — so a test
> branched on a `None` default that *was* the Cartesian case, and every
> curvilinear mesh silently took the slab arm. Two reds, caught only because
> the wrong branch's assertion happened to be falsifiable. ⚠ The aggravator:
> the same string-form check had been run for `mu_start` **one item
> earlier** and come back clean. Now a standing clause in
> `coding-standards`' retirement audit.

**Goal.** Nothing in the streaming path is held that nothing reads.
**Means.** Retire the two flags (with their 12 test assertions **migrated**, not
deleted — the successor spellings exist); retire the three-link dead chain's two
downstream links; retire `_quadrature`; cache the pairing instead of
re-allocating. ✅ **`μ_start`'s branch is MEASURED and the answer is KEEP**
(§3's row, refuted in place 2026-08-26) — so P1's item is no longer
"retire or repair" but: write the measured fraction into
`_edge_seed_stencil`'s docstring, and **retire the stale decision rule in
the `#361` comment at `:1683-1697`**, which still says *"this path may be
dead on the cylinder arm, in which case the answer is retirement"* — the
same too-narrow denominator this plan inherited from it.
⭐ Also fix the falsified docstring the census found: `precompute_psi_state`
(`pole_angular_closure.py:1707`) advertises `psi_view` as
`(N, ng, nx, 1) canonical`; `[M]` the production caller
(`loss_representation/__init__.py:3589`) passes **3-D** `(N, ng, nx)` and a
4-D array **raises** `too many values to unpack` at `:1752`. A falsified doc
is a bug (deferred only to the end of the in-flight merge gate, so the gate
covers one tree).
Add §6's honest refusal for a curvilinear
multi-moment metric — ⚠ **re-sized 2026-08-26**: `[M]` neither
`moment_axis(ndim)` nor `moment_mass_diagonal(ndim)` receives a chart, so this
is a **signature change + refusal** (a §6b call-site change over 8 references),
NOT the guard insertion this line originally implied. See §6's
"Where the slab mass ENTERS".
**Done when.** `grep` finds no reader for each retired name; the µ_start
measurement is recorded with its denominator; the metric guard REFUSES a
constructible input (⚠ §6c — name that input in the gate, or the guard is not
gated); full fast set bit-identical.

#### P1's §6b call-site set — `[M]` audited 2026-08-26, complete BY SYMBOL

⚠ "Complete **by symbol grep**" is the method, stated per §2's FILTER clause.
The 2026-08-24 surprise (a duck-typed kwarg surrogate + an attribute read that
no symbol grep returns) says the consumers' test run is part of this
enumeration, not a check performed after it.

| target | `orpheus/` | `tests/` | `docs/` sources |
|---|---|---|---|
| `requires_upstream_angular_state` | 7 | 6 | 4 — `methods/sn/index.rst:876,884`; `foundations/structured_geometry.rst:733,741` |
| `angular_marching_axis` | 6 | 6 | 2 — `foundations/structured_geometry.rst:734,742` |
| `_quadrature` (the FIELD) | **8** | **2** | 2 (both past-tense, about a *different* object) |

⚠ `index.rst:884` is a Python-domain `:attr:` xref — it renders as plain text
with **no warning at any severity**. Grep is the only gate on it.

⭐ **The audit changed this step: `_quadrature` is NOT dead.** The plan's §3
row called it "a twin reference to one object", which is true, and reads as
*delete it*. `[M]` it has **five live production reads** —
`reduced_operator.py:667`, `:686` (`mu_x[direction_idx]`), `:717`
(`level_indices`), `:719` (`eta[global_n]`), plus the narrowing assert at
`:650` — three writes (`:1065/:1137/:1243`), one declaration (`:565`), and
**two test reads** (`sn/sweep/core/_c_surrogate.py:163,167`, one of which is
an error-message string naming the attribute — grep its SHORTEST distinctive
fragment, not your own rewording).
⟹ the retirement is a **re-point of 7 reads** onto `self.angular.quadrature`,
not a deletion.

⭐⭐ **And that re-point is worth more than tidiness — it dissolves an
`Optional`.** `_quadrature: AngularMeasure | None` is optional; `AngularRedistribution.quadrature`
is **not**. So routing through the angular factor makes the `None` state
unrepresentable and **retires the narrowing `assert` at `:650`** with it —
Pattern 4, and one fewer bare `assert` in `orpheus/` (which `python -O`
strips anyway, so it was never a guard). This is the un-weld paying out: the
twin reference existed *because* the angular factor had no owner.

⚠ **Three false-positive families in the `_quadrature` grep — triage by
MEANING before treating any hit as a site** (raw substring: 16 code / 10 test,
i.e. **2× and 5× over-counted**):
1. every `*_quadrature` symbol (`chord_quadrature`, `select_quadrature`,
   `angular_quadrature`, …) — use `(^|[^a-zA-Z0-9])_quadrature\b`;
2. a **retired BC-shim** attribute of the same name (#176 / C176.1) —
   `test_snmesh_realizer_wiring.py:372,373,433`, `test_bound_compat.py:258,266`,
   and both `docs/` hits. Past-tense history about a *different object*: it
   **STAYS**;
3. a test-local helper *function* `_quadrature(coord)` —
   `test_kinf_homogeneous_tolerance.py:105,159,208`. Unrelated.

### P2 — the angular CLOSURE comes home *(pure `git mv` + imports; bit-identical)*

⛔⛔ **THIS PHASE WAS RE-SCOPED 2026-08-26. Its α half was UNLANDABLE as
chartered — measured, with a reproduction. ✅ Re-scope RULED by the user the
same day (option 1 of 3).** The original text is kept verbatim below because the
refutation is worth more than either half alone (§3).

> **~~Goal.~~** ~~The angular half of the factorization lives in an angular
> package, and the geometry layer stops declaring a Protocol to avoid its own
> input.~~
> **~~Means.~~** ~~`sn/angular/redistribution.py` ← `AngularRedistribution`,
> `angular_redistribution`, `alpha_dome`, `_assert_alpha_dome_closes`.~~
> `sn/angular/closure.py` ← the `AngularClosure` family (out of `sn/sweep/`).
> ~~Shed `StreamingTerms`' angular fields; the `AngularMeasure` Protocol shrinks
> or dies.~~
> **~~Done when.~~** ~~`alpha_dome` has no geometry import;~~ the Protocol's
> member count is ≤2 or the file is gone; `sn/sweep/` contains only traversal;
> `dead_references` 0; bit-identical.

#### ⛔ Why the α half cannot move yet — `[M]` 2026-08-26, reproduced

The four α symbols have **callers that stay behind**: the three `*_streaming`
factories call `angular_redistribution(...)` at `reduced_operator.py:1105`,
`:1172`, `:1275` — **module-runtime calls**, not annotations. Moving the callee
to `sn/` therefore forces a module-scope `geometry → sn` import, and that is a
hard cycle. Not argued — **injected and run**:

```
orpheus/geometry/__init__.py:30      → .reduced_operator
  → orpheus.sn.angular._probe        → orpheus/sn/__init__.py:13 → .solver
    → sn/mesh/augmented_mesh.py:32   → back into the HALF-BUILT reduced_operator
ImportError: cannot import name 'cylindrical_streaming' from partially
initialized module 'orpheus.geometry.reduced_operator' (circular import)
```

`import orpheus.geometry` — the canonical entry — dies. ⚠ And it is
**order-dependent**: importing `orpheus.sn` FIRST happens to work, so a
naive smoke test can pass while the package façade is broken. Note the cycle
does **not** run through `geometry.coord`; it runs through `orpheus/sn/__init__.py`
eagerly importing the solver, so *any* module-scope `geometry → orpheus.sn.*`
import breaks it, whatever the mover needs.

`[M]` the edge does not exist today and this would have created it: over the
whole `orpheus/` tree, by AST — **`geometry → sn` = 0** import statements,
**`sn → geometry` = 24**. The established direction is the opposite one.

⭐⭐ **And that 0 is not history — it is ENFORCED. `[M]` found during P2's own
execution, after the cycle probe.** `tests/test_layer_imports.py` declares
`FORBIDDEN_EDGES["geometry"] = L2_PACKAGES | L3_PACKAGES` and gates it with
`test_no_forbidden_imports` (`@pytest.mark.foundation`, parametrized over every
module in `orpheus/`). So the α move was never a latent hazard to be discovered
— it was a **declared layer violation with a foundation red waiting for it**,
and one grep for `FORBIDDEN_EDGES` would have answered the whole question
before any probe. ⚠ P4 inherits a harder constraint than "avoid the cycle":
`[M]` the linter's `TYPE_CHECKING` tolerance is
`if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES)`, and `geometry` is
`INPUT_PACKAGES` — so **geometry may not import `sn` even for typing**. See P4
for the fork that follows.

⭐ **The done-when named the WRONG DIRECTION of a symmetric relation, and that
is why the plan did not see it.** *"`alpha_dome` has no geometry import"* is
**true** and irrelevant: it asks whether the MOVER depends on its old home. The
failing direction is whether the OLD HOME depends on the mover. One clause, one
direction, reads as complete. ⟹ promoted to `plan-authoring` §6d.

#### ✅ P2's re-scoped charter

**Goal.** The angular closure stops living in the traversal package, and
`sn/sweep/` contains only traversal.
**Means.** `git mv orpheus/sn/sweep/pole_angular_closure.py
orpheus/sn/angular/closure.py` + the import sweep. `[M]` **28** import
statements (26 `from …pole_angular_closure import …` + 2
`from …sweep import pole_angular_closure`), of which **5 are production**
(`loss_representation/__init__.py` ×2 function-local, `sn/mesh/augmented_mesh.py:60`,
`sn/sweep/__init__.py:28` re-export, `sn/sweep/pairing.py:50`). Complete by AST,
so multi-line `from X import (…)` blocks are inside the count — a line-grep
returns 26 and misses two.
⚠ This subsumes §1's `module sn/sweep/pole_angular_closure.py → …/angular/closure.py`
row. **P3 must not redo it.**
**Done when.** `sn/sweep/` contains only `scan.py`, `cache.py`, `pairing.py`,
`psi_half_angle_seed.py`; `sphinx -W` 0; **`grep -rn "sn\.sweep\.pole_angular_closure\|sn/sweep/pole_angular_closure" docs/ orpheus/ tests/`
returns only dated past-tense history**; bit-identical.

⛔ **`dead_references 0` was in this done-when and is REMOVED — it is a
designed-green tell (`plan-authoring` §10, third shape), caught by running it
before designing to it.** `[M]` 2026-08-26 at HEAD: `total_dead: 0,
total_checked: 52, rescued: 52` — and it will read 0 *after* the move too,
whether or not a single xref is repaired, because it judges a narrow decidable
population of 52 while this move touches ~195 references. A tell that cannot
change carries no information (`vv-principles` #19); keeping it would have
certified the sweep instead of the tree.

⭐ **What DOES gate this move, `[M]` and it is exactly three lines.** Neither
`orpheus.sn.sweep.pole_angular_closure` nor `orpheus.geometry.reduced_operator`
is `automodule`'d — `[M]` **0 of 49** autodoc directives in `docs/` name either,
and the built HTML resolves **0** py-domain links into them against positive
controls of 226 (`sn.operators.boundary`) and 90 (`geometry.mesh`) out of 2957
build-wide. So Sphinx never renders these docstrings and `-n` gives **no
differential signal** (≈142 warnings before, ≈142 after). The one real tripwire
is Nexus's own directive resolver: the three `.. implements:: mm-weights`
declarations with `:by: orpheus.sn.sweep.pole_angular_closure.*` at
`docs/theory/methods/sn/curvilinear_one_group.rst:657/:666/:669` emit
`logger.warning(… target not found in graph — skipping)`, fatal under `-W`.
⚠ **Without `-W` that warning is non-fatal and the directive SKIPS** — equation
`mm-weights` silently loses all three `implemented_by` edges and its V&V
coverage degrades with no red anywhere. ⛔ And note the asymmetry for P4: the
four α symbols have **zero** `:by:` declarations, so **the redistribution half
of this arc has no build gate at all** — there, grep is the whole instrument.

⚠ **The scale, `[M]` measured before the move so it cannot be mistaken for
damage this step caused:** 142 Python-domain roles in `docs/**.rst` name a
moving symbol (122 `pole_angular_closure`, 98 of them in `curvilinear_one_group.rst`),
and **they already render as plain text today**. The move does not break them —
it makes them *wrong* as well as unlinked, which is worse, because a
fully-qualified confident path is what a future session copies into an import.
⭐ Precedent, directly on point: `orpheus/sn/spatial` was renamed at task #54 and
`[M]` **47** references survive in `docs/**.rst` — but **0 of the 47 are
Python-domain roles**. The roles were migrated; only literals inside dated
error-catalog history were left. That is the standard to match.

⛔ **The TEST tree does NOT move with the module, and the reason is a coupling
worth knowing before any future test re-home: the DIRECTORY IS THE CAPABILITY.**
`tests/sn/**/conftest.py` × **12** call `stamp_capability_marker(items,
__file__, "<cap>")`, which stamps `@pytest.mark.cap(...)` on every test
collected at or below it — *"The directory IS the capability (single source of
truth)"*, in the conftest's own words. So `git mv`-ing a test file out of
`tests/sn/sweep/curvilinear/` **silently strips its `cap("sweep_curvilinear")`
marker** — no error, no red, just a test that has left the taxonomy.
`[M]` the capability names are hard-listed in `pyproject.toml:97` (**10** nodes),
and the `sentinel` marker on the next line declares *"one sentinel per capability
node"* — so a `tests/sn/angular/` tree is not a directory, it is an **11th
taxonomy node** plus its sentinel. That is a real decision and it is not a
bit-identical move.
⟹ P2 leaves every test file where it is. Only `tests/sn/sweep/curvilinear/test_pole_angular_closure.py`
now carries a filename naming a module that moved; ⚠ it is **not** the only
misplaced one (`test_angular_cell_partition`, `test_march_start_structure`,
`test_tau_arc_wellposedness`, `test_angular_beta_identity`,
`test_tau_producer_equivalence`, `test_compute_psi_half_per_level` all test the
closure from `tests/sn/sweep/`), so moving the one would split the closure's
tests across two trees for no gain. The whole set is one item, and its home is
`.claude/plans/test_architecture_redesign.md` — where the capability DAG lives.

**Out of P2, with its home named** (the §9.3 ruling — a phase number, never an
issue-and-forget): the four α symbols + `_ALPHA_CLOSURE_ATOL` re-home in **P4**,
because §6b's unit of work is the call-site set and α's call sites ARE the three
factories P4 exists to dissolve. `psi_half_angle_seed.py` stays in `sn/sweep/`
for now — `[M]` its five functions are radial-characteristic marches
(`carlson_inward_sweep_from_source`, `radial_characteristic_*`) consumed by
`sn/operators/radial_characteristic.py`; it is neither traversal nor angular
closure, and finding it a home is **not** this plan's business.

### P3 — the names stop lying *(mechanical; bit-identical)*
**Goal.** Every §1 row reads true.
**Means.** The §1 table, one rename per commit, each with the three-search audit.
`ChartConnection` LAST — after P4 has taken `streaming_terms` out of it, so the
name describes what remains.
**Done when.** Each old name has zero live references (past-tense history stays);
`sphinx -W` clean.
⛔ `dead_references 0` is **NOT** a done-when here either — same designed-green
tell P2 removed from its own: `[M]` it judges 52 decidable references and read
`total_dead: 0` unchanged across a 195-reference move. Use the per-symbol
residual grep instead (three filters: qualified form, **string form**,
**attribute access** — see the pre-audit below).

#### §6b pre-audit — `[M]` measured 2026-08-26 at P2's close-out

Each row is the full call-site set, three filters, **nothing clipped through
`head`** (`plan-authoring` §2's VIEWPORT clause — a clipped listing silently
turns a population into a sample). Order the small ones first; `ChartConnection`
is last by ruling.

| symbol | code+test lines / files | string-form | docs `.rst` (roles / total) |
|---|---|---|---|
| `redistribution_gram` → `redistribution_pairing` | ⛔ **9 was the SYMBOL; the CONCEPT is 70** — see below | 0 | 0 / 0 |
| `GeometryCoefficients` → ~~`ChainScanCoefficients`~~ ✅ **`StreamingCoefficientCache`** (both plan candidates refuted — §4bis) | **35** / 7 → `[M]` **45** incl. docs | 2 | 4 / 10 |
| `ReducedStreamingOperator` → `ChartConnection` | 36 orpheus + 15 tests | 3 | — / 20 |

⛔⛔ **The `redistribution_gram` row's `9` was measured on the SYMBOL, and the
work is the CONCEPT — `[M]` 2026-08-26, 7× larger.** The identifier appears 9
times; the *word* `gram` names the same object **70 times** across 7 files: the
guard `_require_single_moment_gram`, the `gram` parameter in 3 signatures, the
`cls(angular, gram)` constructor contract, the `self._gram` field, two message
strings, and a test class `TestGramContract`. Renaming the property alone would
leave the concept spelled **two ways** — a `redistribution_pairing` returning
something every consumer still calls a `gram` — which is strictly worse than
either doing it fully or not at all (Pattern 2, and
[[feedback-naming-consistency-greppable]]). ⟹ done as a concept rename: **60
lines rewritten, 6 deliberately KEPT.**

⭐ **The 6 keeps are the point of the rename, not an exception to it.** A Gram is
⟨bᵢ, bⱼ⟩ over ONE basis — square, symmetric, PSD. This object is
`R_kj = ∫ b_k^scheme · b_j^thread · r dr`, pairing **two different** bases, so it
is rectangular in general. Kept as *Gram*, correctly: `reduced_operator.py:611`
("closing per spatial moment gives a **square** Gram" — that case genuinely is
one), `augmented_mesh.py:1151/1198/1247/1271` (the field-space Hilbert metric
`G_bulk = V·w_n ⊗ diag(1,θ,…)`), and `linear_discontinuous.py:280`. ⚠ The word
appears **493** times repo-wide, nearly all of it legitimate
(`numerics/frame.py`, `numerics/basis/*`) — this is `coding-standards`' concept-grep
dual hazard, so every hit was triaged by MEANING before any was called a site.

⛔ **And a blanket replace would have corrupted `closure.py:1445`** — the word
`pro`**`gram`**`s`. `[M]` substring match = 71, word-bounded = 70. One line, and
it is production prose.

⚠ **All string-form hits are benign and are listed so nobody re-hunts them**:
`GeometryCoefficients` (pre-rename) — a forward return annotation
(`cache.py:205`) and `__all__` (`:506`), both carried through by P3b; `ReducedStreamingOperator` — two `__all__` entries
(`geometry/__init__.py:55`, `reduced_operator.py:1281`) and a forward annotation
(`radial_characteristic.py:155`). None is a `getattr(x, "name", default)`, which
is the form that fails *silently in the default's direction* — that is what the
filter was run to rule out, and it is ruled out.
⭐ The filter was validated against a positive control before being trusted
(`MorelMontryAngularSweep` → 2 known hits). ⚠ Two earlier attempts at this same
census returned a **false `0`** from a zsh quoting error that made the pattern
unparseable — a broken filter and a clean tree are indistinguishable in the
output, which is precisely §2's FILTER clause. Re-run it in Python, not zsh.

#### `SNMesh.curvature`'s audit — `[M]` 2026-08-26, and it is NOT a rename

⭐ **ORDER IT FIRST, and the reason is a dependency the plan did not state.**
`is_cartesian` is the successor spelling P1 migrates the two dead flags ONTO
— and `[M]` `augmented_mesh.py:521` implements it as `return self.curvature
is None`. So the migration target is built on the field P3 retires. Retiring
`curvature` before re-basing `is_cartesian` breaks **12 production sites**
(`solver.py:894,1105`; `dsa.py:203`; `loss_representation` ×6;
`windowing.py`). §6b's shape, with a *definition* in place of a call site.

✅ **But the fix is small, because `is_cartesian` is the CONTRACT and
`curvature` is its implementation detail.** Re-base the one-line body onto
`coord` and all 12 consumers — plus both duck-typed `SimpleNamespace`
surrogates (`test_si_gate_dispatch.py:68`,
`test_unified_sweep_dispatch.py:354`) — never notice, because the surrogates
stub the **property**, not the field. ⚠ That is luck worth naming: those two
files are the *same family* as the 2026-08-24 kwarg-surrogate surprise; had
they stubbed `curvature=` they would have been invisible to every grep here.

**The real consumer set** (prose hits excluded — `[M]` "curvature" appears in
~50 more lines of docstring where it is the *physics word*, not the field):

| site | what | becomes |
|---|---|---|
| `augmented_mesh.py:348,355,1839` | the 3 writes | deleted with the field |
| `augmented_mesh.py:521` | `is_cartesian` body | `coord` — **do this first** |
| `augmented_mesh.py:879` | `if self.curvature is None` | `coord` |
| `loss_representation:712, 2727` | 2 repr / error strings | `coord` |
| `loss_representation:3129-3134, 3169, 3228, 3324` | matvec entry **1** | see below |
| `loss_representation:3479-3481, 3610, 3684, 3764` | matvec entry **2** | see below |
| 5 test reads | `test_radial_characteristic_slot_coordination:94`, `test_native_matvec:392`, `test_axis_native_construction:214`, `test_snmesh_consumes_reduced:192`, `test_sweep_regression:225` | migrate |

⭐⭐ **The two matvec entries are where the value is, and they are TWINS.**
Each opens with the identical three lines —

```python
curvature_raw = getattr(sn_mesh, "curvature", None)
curvature = curvature_raw if curvature_raw is not None else "cartesian"
if curvature not in ("spherical", "cylindrical", "cartesian"):
    raise ValueError(f"Unknown curvature: {curvature!r}")
```

— then dispatch on `curvature == "cartesian"` / `!= "cartesian"`. That is
**stringly-typed dispatch re-derived from a value that is already an enum**,
guarded by a **runtime re-validation of that enum's own domain**, written
**twice** (Pattern 2). Routing to `coord: CoordSystem` makes the `ValueError`
branch unrepresentable (Pattern 4) and collapses the twin. ⟹ this row is a
correctness-adjacent carve wearing a rename's clothes, and it should be
costed as such rather than batched with the cosmetic rows.

### P4 — the pairing is minted where its two halves meet *(behavioural at the seam)*
**Goal.** `R` is produced by the object that owns its basis, from the object that
owns its measure — and `L` and the RC family receive one minted closure.
**Means.** §5's `scheme.mint(chart)`; the angular device becomes the separate
axis it is; the closure splits S0 algebra from S1 binding (F3).

⭐⭐ **AND P2's α half lands HERE — inherited 2026-08-26, with its blocker
discharged by this phase's own work.** Not a punt: §6b's unit of work is the
call-site set, and α's three production call sites *are* the `*_streaming`
factories this phase dissolves. Moving α before that is the measured cycle at
P2's ⛔ block; moving it after costs nothing, because the caller is gone.

The sequence inside P4, and it only works in this order:

1. The factories stop MANUFACTURING the angular factor and are HANDED it —
   `augmented_mesh.py` (already in `sn/`, `[M]` the only production caller, 3
   sites at `:326/:333/:350`) builds the `AngularRedistribution` and passes it.
2. With no caller left in `geometry/`, `alpha_dome`, `_assert_alpha_dome_closes`,
   `_ALPHA_CLOSURE_ATOL`, `AngularRedistribution` and `angular_redistribution`
   move to `sn/angular/redistribution.py`. No runtime `geometry → sn` edge is
   ever created, so the P2 cycle cannot arise.

   ⛔ **BUT STEP 1 IS NOT INDEPENDENTLY LANDABLE, and this text said it was —
   caught 2026-08-26 by re-reading it against the linter, hours after writing
   it.** Killing the runtime CALL does not kill the reference: whatever holds
   the angular factor — the factory's parameter, or (if the factories collapse
   into it, see below) `ChartConnection`'s own field — still needs the **NAME**
   for its annotation. And `geometry` may not import `sn` for typing either.
   ⟹ **the sharp statement of P4's constraint, which supersedes the step
   ordering above:** *`ChartConnection`'s `angular` field needs a type, and that
   type cannot live in `sn/` while `ChartConnection` lives in `geometry/`.* Steps
   1 and 2 are therefore **one step**, and it cannot be taken until the fork
   below is ruled — which is §6b's rule (the unit of work is the call-site set)
   applied to a *type* reference rather than a call.
   ⚠ Do not read this as "step 1 first, then decide". Decide first.

   ⛔⛔ **AND THE CONSTRAINT IS HARDER THAN THE CYCLE — `[M]` 2026-08-26, found
   during P2 execution.** The tree ships an **import-linter**,
   `tests/test_layer_imports.py`, and `geometry → sn` is a **declared forbidden
   edge**, not an accident of history:
   `FORBIDDEN_EDGES["geometry"] = L2_PACKAGES | L3_PACKAGES`, gated by
   `test_no_forbidden_imports` (`@pytest.mark.foundation`, parametrized over
   every module). So the P2 move would not merely have failed at import — it
   would have gone **red on a foundation gate**, which is the better outcome and
   the one the plan should have anticipated by reading the linter.

   ⛔ **REFUTED, same day, by its own author.** This block first read
   *"cycle never arises: geometry annotates via `TYPE_CHECKING`"*. That is
   **false**, and it was written without reading the tolerance it invoked.
   `[M]` `_check_module` applies it as
   `if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES): continue` — the
   tolerance covers **`numerics` and `transport` only**, and `geometry` is in
   `INPUT_PACKAGES`. ⟹ **`geometry` may not import `sn` under `TYPE_CHECKING`
   either.** A `TYPE_CHECKING` annotation is not an escape hatch here.

   ⛔⛔ **SUPERSEDED 2026-08-27 — the fork was an answer to the wrong question.
   The four options below STAY (§3), with their refutations, at §4ter's
   "The four options this supersedes". Read §4ter FIRST; it re-poses P4.**
   `[M]` `geometry/reduced_operator.py` holds no geometry: its `face_areas` is a
   verbatim copy of `mesh.areas`, `coord` a copy of `mesh.coord`, `delta_A` a
   `np.diff` of the former with **0** non-SN consumers, and the module has **0**
   intra-`geometry/` consumers against 1–4 for every other geometry primitive.
   The module dissolves; the forbidden edge never arises.

   ⟹ **the ORIGINAL fork text, kept because its refutations are the record:**
   (a) `ChartConnection` leaves `geometry/` for `transport/` (L2 — which *does*
       carry the tolerance) or for `sn/`; (b) the α type lands somewhere
   `geometry` may import — `numerics/` (L1) is arguable, since
   `alpha_dome(cosines, weights)` is pure angular-measure math with no chart
   and no neutron; (c) `ChartConnection` stops naming the type at all, because
   after step 1 it only *stores* what it was handed; (d) a `WHITELIST` entry
   with a `RETIRE_IN_…` trigger — the tree's own idiom for a transitional
   violation (`[M]` 4 live entries), but an exemption is not a design.
   ⚠ Note (b) sits in tension with §9.2's fork-2 ruling (`sn/angular/` over
   `transport/angular/`) — that ruling weighed *consumers*, and the layer
   contract is a *different* constraint it did not see. Re-open it explicitly
   at P4 rather than treating either as settled.

   ⭐ And read `tests/test_layer_imports.py`'s own module docstring first: it
   already documents this exact failure mechanism for the `numerics ↔ geometry`
   back-edge — *"a partially-initialised package can serve
   `from orpheus.numerics.measure import X` (a SUBMODULE import) but NOT
   `from orpheus.numerics import X` (an ATTRIBUTE lookup on a package whose body
   has not reached that line yet)"*. That is precisely what the P2 probe hit
   (`from orpheus.geometry.reduced_operator import cylindrical_streaming`, an
   attribute lookup on a half-built module). The knowledge was authored, gated
   and inert — nobody read it.
3. `AngularMeasure` **dies outright** — it does not "shrink to ≤2". `[M]` the
   three surviving factories read **0 of its 6** members (`mu_x`, `weights`,
   `N`, `eta`, `mu_z`, `level_indices`); they only forward the object. Its four
   real readers are all inside `angular_redistribution`'s body and travel with
   it. ⛔ And its docstring's stated reason for existing — *"the geometry layer
   needs no import from the quadrature package at all"* — is already false:
   `[M]` `geometry/boundary/` imports `Quadrature` from
   `orpheus.numerics.quadrature` under `TYPE_CHECKING` in **5** modules
   (`_base.py:63`, `white.py:20`, `_specular.py:94`, `albedo.py:27`,
   `reflective.py:25`), plus a runtime function-local at `_realizer.py:297`. So
   the replacement for the Protocol is a `TYPE_CHECKING` import, with five
   precedents in its own package.
4. ⚠ `[M]` `reduced_operator.py:1229` reads `getattr(angular_measure,
   "level_structure", None)` — a member the Protocol never declared. It is a
   string-form read that no symbol grep for the six members returns
   (`coding-standards`, the clause P1 paid for). Retiring the Protocol must
   account for it, not just for the declared six.

⭐ **What the factories BECOME is the elegance question, and it is a
Pattern-2 finding, not a re-home.** Once handed their angular factor, the three
differ *only* in the enum they validate and forward:
`if mesh.coord is not X: raise` then `ChartConnection(coord=X, mesh=mesh,
angular=angular)`. That is one function written three times — and the honest
collapse is the dataclass's own constructor with a `mesh.coord == angular.coord`
consistency check in `__post_init__`, at which point the three factories are a
**retirement**, not a move. Decide it here; do not pre-commit to it in a plan
row.

⭐⭐ **AND the moment mass comes with it — they are ONE bilinear form.**
*(added 2026-08-26 on a user challenge to P1 item 9's signature.)* F2 gives

```
M_kj = ∫ b_k b_j          dV      # the moment mass  — Galerkin matrix of 1
R_kj = ∫ b_k b_j (∇·ê_r)  dV      # the pairing      — Galerkin matrix of ∇·ê_r
```

— the **same construction against different integrands**. So
`moment_mass_diagonal` and `redistribution_pairing` are two products of one
form, and the mint produces **both**. This is not a new abstraction; it is
noticing that the object §5 already mints has a second output.

⛔ **P1 item 9 left a TRANSITIONAL parameter that P4 must retire, and it must
not be inherited as a decision.** `moment_mass_diagonal(ndim, coord)` and
`moment_axis(ndim, coord)` take the chart as a **tag** — and a tag is the
wrong thing to take, because what the mass needs is the *measure*.
Discriminating on a tag to recover a measure is the missing-type smell
`nexus discriminations` hunts.
⭐ The tell that it is a proxy: the guard's predicate is not three-way but
`is not CARTESIAN` — *"is the volume element constant across the cell?"* — a
property of the **measure's behaviour**, not the chart's identity.
⚠ It shipped because the guard is a correctness repair (`[M]` the slab's mass
was installed on sphere AND cylinder, bit-identical, silently) and could not
exist while the producer was never told its chart. One tag-dispatch was the
smallest thing that makes the wrong value unspellable *today*.
`[M]` nothing else could supply it: `ndim` does not determine the chart (1-D
is slab/cylinder/sphere), and the schemes carry **no** chart state by design
(`DiamondDifference` none, `LinearDiscontinuous` only `theta`) — giving a
scheme a chart field would re-weld two axes this campaign separated.

**Done when.** The scheme is the only producer of `R`; two meshes over one
`(quad, coord)` share one `angular_algebra` object by identity; `L` and the RC
operators receive the closure and neither constructs one; ⚠ a gate pins that the
rank-1 (ONETRAN) pairing is *expressible* — §6c: it must be constructible, not
merely refused; **and `grep -n "coord: CoordSystem" transport/spatial/scheme.py`
is empty** — the mint has replaced the tag, and the moment mass is evaluated
against a measure rather than branched on a chart.
**Plus P2's inherited half:** `grep -rn "alpha_dome\|AngularRedistribution"
orpheus/geometry/` is empty; `grep -c "AngularMeasure" orpheus/` is 0;
`pytest tests/test_layer_imports.py` green **with no new WHITELIST entry**
(a whitelisted violation is a deferral, not a done-when); and
`.venv/bin/python -c "import orpheus.geometry"` succeeds **in a fresh
interpreter with nothing else imported first** — ⚠ that qualifier is the whole
gate, since `[M]` importing `orpheus.sn` first masks the cycle entirely.

### P5 — `ChartConnection`'s three stages separate *(rides O-3)*
**Deferred on a DEPENDENCY, not a punt.** Its third piece is `streaming_terms`'
evaluated per-`(cell, ordinate)` view, which is exactly the layer O-3 redesigns.
Splitting it here would churn that seam twice.
**Unblocked when:** O-3 lands its per-`(cell, ordinate)` layer. That is the
whole gate — checkable, and named here so this row cannot become a
defer-and-forget (the ruling at §9.3 applies to every "later" item in this
plan, and P5 is the one whose blocker is genuinely external).

### P6 — the chart speaks ONE vocabulary *(mechanical; bit-identical)*
**Goal.** One concept, one spelling. A grep for the chart returns the whole set
rather than a confident partial answer.
**Why it is a phase and not an issue.** §9.3: scheduled, not deferred. It is
last because it is a ~400-hit sweep across files P1–P4 never open — landing it
inside a phase whose done-when is bit-identity would dominate that phase's diff
and dilute its evidence. Its own phase, its own diff.
**Means.** Collapse `geometry_kind` (105) → `coord` (228). `curvature` (68) is
already gone at P1/P3. ⚠ Three-search audit per §1's warning, and grep the
**case-spellings** (`[M]` 8 / 3 / several) — a Python-domain xref renders as
plain text with no warning at any severity.
**Done when.** `grep -rn "geometry_kind" orpheus/` is empty; `coord` is the only
live spelling; `dead_references` 0; `sphinx -W` clean; full fast set
bit-identical.

### P7 — the space metric admits a NON-DIAGONAL `G` *(behavioural; the terminal phase)*
**Goal.** A function space can carry a metric that is not a Hadamard weight, so
a cell-dependent non-diagonal `M/V` is *expressible* — and §6's honest refusal
can become an honest value.

⭐ **This is NOT "choose `G` for curvilinear LD", and the difference is what
makes it schedulable.** Choosing `G` needs #158's cell solve to have a
consumer (§6) and is genuinely blocked. What is **not** blocked is the
machinery, and `[M]` 2026-08-26 the machinery is absent by construction:

| surface | `orpheus/numerics/space.py` | shape |
|---|---|---|
| `inner_product` (dense arm) | `float(np.sum(w * x * y))` `:542` | Hadamard |
| `inner_product` (axes arm) | `np.sum(_apply_axes_weights(x) * y)` `:538` | per-axis diagonal |
| `apply_metric` | docstring spells it `G⊙x` `:570` | Hadamard **by contract** |
| `apply_inverse_metric` | `(1/G)⊙x` `:592` | Hadamard |
| the one extension point | composite spaces override `apply_metric` per BLOCK `:577` | each block still diagonal |

⟹ there is **no** path for a non-diagonal `G`, at any level. `[M]` the
`__post_init__` guard at `:214` calling `inner_product_weights` "dense" means
*a dense array of diagonal weights*, not a matrix — which is exactly the sort
of near-miss vocabulary that makes an absent capability read as present.

⭐⭐ **Cardinal Rule 2 — TWO campaigns demand this one primitive**, which is
what promotes it from this plan's loose end to a shared foundation:

| campaign | the non-diagonal `G` it needs | what it does today |
|---|---|---|
| CS4b frame square, F-0 | the **slab** frame's live discrete Gram, off-diag ~1.15 | records `parseval unavailable` as a limitation; the Parseval gate SKIPS slab |
| this plan, §6 | curvilinear multi-moment `M/V` = `[[1, 0.5], [0.5, 0.4]]`, cell-dependent | installs the **slab mass** `[1, 0.3333]`; refused one layer out (#158) |

`[M]` recorded independently in three places as a "CS4c matrix-metric debt"
(`frame_square_recarve.md:122,325`; `collapse_pair_frame_induction.md:103`;
`orpheus-operator-machinery-report-v2.md:634`) — three records, one missing
capability, no owner. This phase gives it one.

**Means** (proposed, NOT verified — design at the phase, per §1's discipline):
the metric becomes a thing that is *applied* rather than a thing that is
*multiplied* — i.e. the Hadamard weight is the diagonal SPECIAL CASE of an
operator-valued metric, and `apply_metric`/`apply_inverse_metric` are its two
faces. That is the same shape CS4c already named ("the leg-as-operator
machinery"), so the two campaigns should converge here rather than each
minting one. ⚠ Check CS4c's state before designing — if it has landed legs,
this phase CONSUMES them and shrinks to the curvilinear installer.

**Done when.** A space carrying a non-diagonal `G` satisfies
`⟨x, y⟩ = yᵀGx` and `apply_inverse_metric ∘ apply_metric == id` on `G`'s
range; the CS4b slab arm's Parseval gate STOPS skipping (⚠ §6c — that
un-skip IS this phase's witness, and it must be named in the gate); §6's
curvilinear refusal is re-derived — either it becomes a value, or it states
which of the two blockers still stands. Euclidean and diagonal paths
bit-identical.

---

## 8. Explicitly NOT in scope

- **#158's cell solve.** `n_mom > 1` stays refused (`6859ca05`'s guard).
- **#407's carve.** Both instances are recorded; the fix belongs with O-3, where
  `cell_balance.py` reorganises under the scheme family anyway.
- **Re-baselining the `cyl_2g_3reg_folded_4x8_dd_n40` red.** `[M]` triply
  confirmed pre-existing and bit-identical across three independent carves
  (#404); it is `@pytest.mark.slow` and deselected by the canonical gate.
- **Choosing `G`** — the VALUE. §6 says what CANNOT adjudicate it; choosing it
  is its own campaign, pinned by physical functionals, and it needs #158's cell
  solve to have a consumer.
  ✅ **AMENDED 2026-08-26 (user ruling, §9.4) — the row was too wide.** It read
  as excluding the *machinery* along with the *value*, and those have different
  blockers: `[M]` the value is blocked on #158, the machinery is blocked on
  nothing and has **two** waiting consumers. The machinery is now **P7**, a
  scheduled terminal phase. Only the VALUE stays out of scope.
  > ⚠ This is plan-authoring §3's REMEDIED-fact case: the row was true when
  > written and is now half-false, and nothing would have prompted its edit
  > because being true is exactly why nobody re-reads it.

---

## 9. Forks — ✅ ALL FOUR RULED 2026-08-26 (user)

The fork text stays (plan-authoring §3: a decided premise is edited in place,
never dropped — the alternative's REASON is what stops it being re-litigated).

1. ✅ **`ChartConnection`** — the first names what it is (connection
   coefficients on a chart, the differential-geometry term the α-dome already
   cites); the second is more familiar and less right.
   ⛔ **`StreamingCoefficients` REJECTED**, and the reason is structural, not
   taste: *"streaming" is what `L` DOES*. Lending `L`'s vocabulary to a
   non-operator struct is the exact defect `ReducedStreamingOperator` is being
   retired for, so the alternative reproduces the bug it is meant to fix.
   ⭐ And the positive evidence is **F2**: after P1's retirements and P4's
   extraction the residue is `face_areas` + `delta_A`, and
   `ΔA = ∮ê_r·n̂ dA = ∫∇·ê_r dV` with `∇·ê_r = (d−1)/r` — the *contracted
   connection coefficient*. The residue is connection data, not metric data,
   so the name is literally right rather than merely evocative.
   ⛔ **PREMISE MOVED 2026-08-27 — re-rule before using this name.** The ruling
   rests on *"the residue is `face_areas` + `delta_A`"*. §4ter retires
   `face_areas` (it is `mesh.areas`) and derives `delta_A`, so the residue is a
   **producer of sweep coefficients**, not connection data. The name may still
   be right; the ARGUMENT for it no longer holds as written. ⚠ §1: grep the
   prose corpus before adopting any replacement.
2. ✅ **`sn/angular/`.** `[M]` today the closure family has
   **zero** non-`sn` importers, so `sn/` is honest. But R19 requires the *scheme*
   family to serve diffusion; if the angular family is ever expected to serve a
   second method, `transport/` is the forward-compatible home. My read: the
   angular-cell closure is genuinely SN's (diffusion has no angular axis to
   close; MoC marches space; CP is integral) ⟹ `sn/angular/`.
   ⛔ **`transport/angular/` REJECTED** — it asserts cross-method generality
   with **0** candidate consumers. A later move is a pure `git mv`.
3. ✅ **SPLIT, and the remainder is SCHEDULED — not deferred.** Retire
   `SNMesh.curvature` inside P1/P3 (it is a stringly-typed *duplicate* of the
   `CoordSystem` enum read through a defaulted `getattr` — a live defect, and
   already load-bearing for P1's flag retirements). The `geometry_kind` ↔
   `coord` *synonym* unification is **P6**, a terminal phase of THIS plan.
   `[M]` 228 + 68 + 105 hits is a large mechanical sweep that touches
   files this plan otherwise never opens — which is why it gets its own phase
   and its own diff, **not** why it gets postponed out of the plan.
   > ⭐ The user's ruling, and it governs #4 too: *"we'll do option 1 to not
   > dilute the diff, but the rest of the work should be scheduled to the end
   > of this plan, **not defer and forget**."*
4. ✅ **Mint the pairing (+ spatial closure) at P4; the moment measure is
   P7, a terminal phase of THIS plan.** Not `.moment_measure` at P4 — that
   member requires choosing `G`, and §6 measured that no gate available here
   can adjudicate the choice. But per the ruling above it is **scheduled, not
   deferred**: see **P7**, whose scope is materially different from "choose
   `G`" and is unblocked (§7).

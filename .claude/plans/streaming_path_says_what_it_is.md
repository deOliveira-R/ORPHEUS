# The streaming path's objects say what they are, live where they belong, and carry only what they use

> ## ▶ RESUME STATE — rewritten at the 2026-08-26 post-P3 compaction point
>
> **This file is the resume surface. Trust it over any summary.**
> ⛔ This replaces a post-P1 header. If a summary still says "▶ NEXT = P2", or
> describes P2 as moving α, or sizes a rename in single-digit lines, it is
> quoting dead text: P2 was re-scoped, P3 landed, and BOTH plan rename
> candidates for the cache class were refuted. See §4bis and §7 P2's ⛔ block.
>
> ### Where the tree is — `[M]` reconcile with git, never with this
>
> `main` @ **`8ffddfb9`** == `origin/main`. ✅ **P2 AND P3 are both MERGED and
> PUSHED; no open branch for this campaign** (`refactor/p2-angular-comes-home`
> and `refactor/p3-names-stop-lying` both deleted). `[M]` P3's exit gate ran
> BEFORE the merge: **9815 passed / 0 failed**, 13 trees all `rc=0`, and
> **delta 0 on every one of the four axes** against the P2 gate — two renames,
> nothing moved anywhere.
> ⚠ Reconcile with git, never with this line
> (`git merge-base --is-ancestor 4a3f2390 main`).
> ⚠ This repo has **no CI** (`.github/workflows/` absent). "main is always
> green" is enforced entirely by the local gate; an empty `gh run list` is not
> a pass.
>
> ### What has LANDED
>
> | commit | phase | item |
> |---|---|---|
> | `cc01dd27` | Phase B | the angular closure takes its two tensor factors |
> | …`1c93b14f` | P1 | nothing in the streaming path is held that nothing reads (7 commits) |
> | `dcd6a9f6` | **P2** | `sn/sweep/pole_angular_closure.py` → `sn/angular/closure.py`, byte-identical |
> | `75571c4d` | **P3a** | `redistribution_gram` → `redistribution_pairing` (concept: 60 lines, 7 files) |
> | `4a3f2390` | **P3b** | `GeometryCoefficients` → `StreamingCoefficientCache` (45 sites, 12 files) |
> | `8ffddfb9` | — | this compaction point (P2+P3 merged to main here) |
>
> `[M]` P2's exit gate: **9815 passed / 0 failed**, 13 trees all `rc=0`; against
> P1's baseline **+1** passed — the new `sn/angular/__init__.py` adds one case to
> the import-linter's `rglob` parametrize — and **0 on every other axis**
> (skipped 22, deselected 227, xfailed 70). `[M]` P3: `tests/sn` +
> `tests/transport` **3763 passed / 0 failed**, then the full gate **9815/0**,
> **delta 0 on all four axes** vs P2. pyright 0 and `sphinx -W` 0 throughout.
>
> ### ▶ NEXT — P4, RE-POSED 2026-08-27. **Read §4ter before anything else.**
>
> ⛔ P4's four-way fork ("where does the α type live so `geometry` may name it?")
> was an answer to the WRONG QUESTION and is superseded. `[M]`
> `geometry/reduced_operator.py` **holds no geometry** — `face_areas` is a
> verbatim copy of `mesh.areas`, `coord` a copy of `mesh.coord`, `delta_A` a
> `np.diff` of the former with **0** non-SN consumers, and the module has **0**
> intra-`geometry/` consumers against 1–4 for every other geometry primitive.
> ⟹ the module **dissolves**, and the forbidden edge never arises rather than
> being routed around. Full record, the five ordered steps, the four superseded
> options with their refutations, and **three OPEN forks** (a fourth — where
> `delta_A` lives — was ruled 2026-08-27: `sn/`, on a CONTRACT argument, not a
> consumer count):
> **§4ter**.
>
> ⚠ Two things §4ter corrects that a summary will otherwise carry forward:
> **(a)** this plan asserted the field retirement was bit-identical — `[M]` it is
> for `coord` (3/3 charts) and is **NOT** for `face_areas` (2/3; the slab holds
> `None` where `mesh.areas` is `ones`). That step is *finishing P1's un-weld on
> the spatial half*, not a duplication retirement. **(b)** §9.1's
> `ChartConnection` name rests on a residue §4ter changes — re-rule it.
>
> The mint (below, item 1) and P4b (item 4) are unchanged by any of this.
>
> Full text at §7 P4 **and §4bis**. ⚠ Read both; P4 has absorbed work from two
> other phases.
>
> 1. **The mint** — `scheme.mint(chart)`; `R` produced by the object that owns
>    its basis. Retires P1's TRANSITIONAL `coord: CoordSystem` tag. Tell:
>    `grep -n "coord: CoordSystem" transport/spatial/scheme.py` → empty.
> 2. **The moment mass comes with it** — `M` and `R` are the Galerkin matrices
>    of `1` and `∇·ê_r`: one bilinear form, two products, one mint.
> 3. **P2's α half** (inherited) — the four α symbols + `_ALPHA_CLOSURE_ATOL`
>    move to `sn/angular/redistribution.py`; `AngularMeasure` **dies** (`[M]`
>    the 3 surviving factories read **0 of its 6** members).
> 4. **P4b** (new, §4bis) — split the cache's three strata.
>
> ⛔⛔ **P4's hard constraint is ENFORCED, not advisory.**
> `tests/test_layer_imports.py` declares
> `FORBIDDEN_EDGES["geometry"] = L2_PACKAGES | L3_PACKAGES`, gated on every
> module by a `@pytest.mark.foundation` parametrized test. `[M]` its
> `TYPE_CHECKING` tolerance is `if is_tc and src_pkg in (L1_PACKAGES |
> L2_PACKAGES)` — `geometry` is `INPUT_PACKAGES`, so **geometry may not import
> `sn` even for typing.** ⟹ the sharp form: *`ChartConnection`'s `angular` field
> needs a type, and that type cannot live in `sn/` while `ChartConnection` lives
> in `geometry/`.* Four options at §7 P4 — **rule the fork BEFORE editing**;
> its steps 1 and 2 are one step, not two.
>
> ⚠ **P3 is NOT finished.** `ReducedStreamingOperator` → `ChartConnection`
> (P3c) is sequenced deliberately **after P4**, so the name describes what
> remains once `streaming_terms` is out of it.
>
> ### ⭐ The three rules P1–P3 paid for — read before the next rename or re-home
>
> - `plan-authoring` **§6d** (P2): **a RE-HOME must check the import edge in
>   BOTH directions — and look for an import-linter FIRST.** The done-when that
>   failed was *true*: it asked whether the MOVER depends on its old home, when
>   the question is whether the OLD HOME depends on the mover.
> - `plan-authoring` **§1** (P3b): **a proposed NAME's ABSENCE can be evidence
>   AGAINST it.** `SweepCoefficientCache` was free in code *and* in git history
>   **because it was rejected** — `lessons.md` L15 records it as "the wrong
>   shape". Grep the PROSE corpus (lessons, plans, agent-memory, issues).
> - `.claude/lessons.md` **L61**: **an unvalidated filter and a clean tree print
>   the same thing.** Six false negatives in one session, two mechanisms — zsh
>   eating quotes/backticks from double-quoted patterns, and ⛔ **`grep` here is
>   `ugrep`**, where an anchor inside an alternation group returns 0 with **no
>   error at all** (now in `nexus-tools.md`). ⟹ run every COMPLETENESS check in
>   Python, with a positive control.
>
> ### ⭐⭐ The finding that outlives all three phases
>
> **A lesson's own worked example is not exempt from the lesson.** `lessons.md`
> L15 is *"cache shape that mixes immutability strata hurts twice"* and holds
> `GeometryCoefficients` up as the RIGHT shape — while `[M]` (3 meshes, one
> quadrature) it mixes **three**: 7 fields bit-identical across meshes, 4
> mesh-bound, 2 traversal-only. Nothing prompts that re-check, because the
> exhibit is what taught you the rule. → §4bis, and P4b.
>
> ### The rulings that bind (✅ user; full text + rejected alternatives at §9)
>
> `ChartConnection` · `sn/angular/` · the chart sweep SPLITS · P4 mints the
> pairing, the moment measure is **P7** · P2 re-scoped to the closure move only
> · the cache is **`StreamingCoefficientCache`**.
>
> ⭐⭐ **The governing ruling:** *"scheduled to the end of this plan, **not
> defer and forget**."* Deferred work earns a PHASE NUMBER or a named external
> blocker. P4b, P6 and P7 exist because of it.
>
> ### ⚠ `scratch/` is UNTRACKED — a `git clean` destroys it
>
> The four audit memos below, the gate logs (`_p2_fast_gate.log`,
> `_p3_fast_gate.log`) and ~190 probes.
>
> ### Running the gate
>
> `.venv/bin/python -O -m pytest -p no:randomly -m "not slow" -q` — SERIAL,
> **~63 min**. Use the per-tree driver (`scratch/_p3_gate_driver.sh`; copy it,
> change its `LOG=` line) launched DETACHED via
> `Popen(start_new_session=True)` — harness background tasks die at ~30–90 min.
> ⛔ Never pipe it through `tail`. ⚠ `-m "not slow"` is part of the number: it
> deselects 227, including the known-inherited `cyl_2g_3reg_folded_4x8_dd_n40`
> red (#404, dup of #397, `[M]` bit-identical pre-branch).

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

**P4.1a — `coord` retires.** `[M]` duplicate on **3/3** charts
(`op.coord is op.mesh.coord`). Bit-identical. **This is the load-bearing
prerequisite**: it is what severs `transport`'s reach into the bundle
(`transport/radial_characteristic_field.py:361-362` reads `mesh.reduced.coord`
→ `mesh.coord`). Without it, P4.4 breaks an L2 consumer.

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

**P4.3 — `StreamingTerms` to L2.** `transport → geometry` stays legal (it reads
`mesh.areas`); no new edge either way. Fix the *"purely geometric"* claim while
the file is open.

**P4.4 — the residue to `sn/`, and the module is deleted.** Safe only after
P4.1a. `[M]` 0 intra-geometry consumers, so nothing stays behind.

**P4.5 — the hollow Cartesian object dies.** `[M]` d≥2 Cartesian gets
`reduced = None`; the slab gets an object whose every meaningful field is
`None`/neutral, minted *"so `sn_mesh.reduced` is always populated"* — a promise
the d≥2 arm already breaks. `reduced is None` today conflates **"no chain scan"**
(`ndim ≥ 2`) with **"no curvature"** (`coord is CARTESIAN`), and `SNMesh` answers
both directly. `[M]` **12 of 36** code reads of `.reduced` are `None`-guards
paying for the conflation. The object should exist **iff there is a radial
reduction**.

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
2. **Where the residue lands inside `sn/`** — `sn/mesh/` (beside its only
   constructor) or `sn/sweep/`. Depends on whether §5's `scheme.mint(chart)`
   absorbs `redistribution_pairing`, which would leave much less behind. Decide
   WITH the mint, not before it.
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

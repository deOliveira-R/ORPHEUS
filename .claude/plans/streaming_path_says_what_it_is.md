# The streaming path's objects say what they are, live where they belong, and carry only what they use

> ## ▶ RESUME STATE — ⏸ COMPACTION POINT, rewritten 2026-08-28 after P4.3 (P4.4 + P4.2 landed earlier the same day)
>
> **This file is the resume surface. Trust it over any summary.**
> ⛔ A summary is quoting dead text if it says any of: "▶ NEXT = P2" · P2 moves
> α · P4 has a four-way fork about where the α type lives ·
> `face_areas`/`delta_A` are fields · `streaming_terms` has three chart arms ·
> **"P4.3 before P4.4"** · **"P4.2 before P4.4"** · `AngularMeasure` "dies" at
> P4.2 · **the α movers are five** · **`transport/spatial/` should move to
> `sn/spatial/`** · the connection operator or the α cluster is still in
> `geometry/` · **"▶ NEXT = P4.3"** · **`geometry/reduced_operator.py` still
> exists** · `StreamingTerms` is "purely geometric" or lives in `geometry/` ·
> P4.9's edit list is five small edits · **"▶ NEXT = P4.9" or "P4.9a is
> un-landed"** · the Morel–Montry twin still sits in `diamond.py` ·
> `cell_balance_terms` exists · `CellVisit` carries `tau`/`c_in`/`c_out` ·
> the 59 %/204-ULP figure quoted as a point value (it is one draw's
> reading — the seed-stable discriminator is max |Δ| = 1.776e-15) ·
> **"P4.9b is un-landed" or "▶ NEXT = P4.9b"** · `StreamingOperator(sn_mesh)`
> takes only the mesh · **"SNMesh sheds scheme + pole_angular_closure at
> P4.9b"** (⛔ REVISED BY RULING — the hub KEEPS both; consumers flow through
> the operator) · `pole_angular_closure` / `PoleAngularClosureBase` exist
> (renamed `angular_closure` / `AngularClosureBase`; only the retired-Protocol
> HISTORY and the `sn-pole-angular-closure-protocol` anchor keep the old
> spelling) · the walk reads `mesh.scheme` at apply time · the scan cache is
> solver-built / mesh-memoized (`_geom_cache` is retired; the strategy-layer
> intern `geometry_cache_for` is the one home; `_coll_cache` +
> `_pole_mirror_cache` deliberately remain) · the 24.65 % figure quoted
> without its fixture (the operator count is fixture-dependent, 6–10 vs
> 38–43/solve; the COUNT gate is the instrument) ·
> **`L = D·T_ang·T_spatial` stated as a fact, or §5b ⛔ P0 as
> unblocked-but-unrun** (⛔ P0 ANSWERED 2026-08-29: the product form is
> REFUTED with a measured witness; the measured algebra is the SUM
> `L = D − E_sp·S − E_ang·A_ang` — `scratch/p0_product_form_measurement.md`)
> · the cylinder called "non-carrying" as a live class (slab-only since
> Q5.6.3 — the admitted fold is always carrying) · **P4.5 as un-landed, the
> slab's `reduced` as "hollow", or `reduced is None` as conflating two
> facts** (✅ P4.5 landed `260ddc64`; P4.1b made the slab object
> load-bearing; `[M]` presence ⟺ `is_1d` and the guards now SAY so) ·
> `_is_curvilinear` exists · **the moment-mass family takes `(ndim, coord)`**
> or its guard re-interprets a chart tag (✅ P4.6 `2ec73b80`: it consumes
> `tuple[Axis1D, ...]` and asks `has_constant_volume_element`) ·
> **`StreamingTerms` carries `delta_A_over_w`/`chord_length`/`mu`**, the
> cache spells `dA_w`/`A_down`/`A_total`/`V`, or "P4.7 became a
> consequence" (✅ P4.7 `3456dd37`: the packet is 4 fields; the family
> speaks the long names; the `[R]` was refuted — two deliberate readers).
> Every
> one was true once, or was proposed and refuted, below.
>
> ⭐ **Three phases landed 2026-08-28 — P4.4, P4.2, P4.3 — and P4.9 was
> chartered, then SHARPENED (the direction-supplier ruling), then
> PRE-MEASURED.** If a summary does not mention P4.9, it predates the charter;
> if it calls P4.9 five small edits, it predates
> `scratch/p4_9_design_measured.md`.
>
> ⚠ **Line numbers re-measured at `fe0f43ba` (post-P4.7, 2026-08-29).**
> `streaming.py` — `StreamingOperator` `:174`, **`.pose` `:297`**
> (unshifted); `closure.py` — **`AngularClosureBase` `:242`**,
> **`march_psi_half_step` `:1372`** (the owner; +2 vs the P4.9b reading);
> `loss_representation/__init__.py` — **`geometry_cache_for` `:479`**
> (the strategy-layer intern), `_OctantWalk` `:1014`, `default_for`
> `:2776`, `_OneDimScanWalk` `:2909` (all unshifted); `cache.py` —
> `StreamingCoefficientCache` `:126`, **the fields RENAMED at P4.7**
> (`face_area_downstream`/`face_area_total`/`delta_A_over_w`/`volume`,
> `:240-243` — any note spelling `dA_w`/`A_down`/`A_total`/`V` predates
> P4.7), `from_mesh_and_quad(sn_mesh, angular_closure)` `:259`, the
> `is_1d`-keyed admission `raise` `:286` (re-keyed at P4.5);
> `coupled_system.py` — `build_streaming_collision` `:423`;
> `transport/spatial/scheme.py` — `StreamingTerms` `:108` / `CellVisit`
> `:304` UNSHIFTED; `augmented_mesh.py` — `scheme` stored `:268`,
> **`angular_closure`** bound `:422` (the attr's new name; both DELIBERATELY
> retained — the hub ruling). Live sections were refreshed 2026-08-29; DATED
> records (incl. §5b's pre-carve measurements and P4.9a's own landmark
> paragraph if quoted anywhere) keep the numbers true at their date.
>
> ### Where the tree is — `[M]` reconcile with git, never with this line
>
> `main` == `origin/main`, **no open branch for this campaign**, working tree
> clean. ✅ Phase B, P1, P2, P3a/b, P4.1a/b/c, P4.4, P4.2 and **P4.3** are all
> merged and pushed. Reconcile with
> `git merge-base --is-ancestor da507e3d main`, never with this sentence.
> `[M]` 2026-08-28 (post-P4.3): **all hashes this file cites resolve AND are
> ancestors of `main`** (re-verified at this compaction point). ⚠ `refactor/operator-inverse-algebra` exists locally and is **already
> merged** — a stale local branch, not open work.
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
> | `da507e3d` | **P4.3** | `StreamingTerms` to L2 beside its contract; `geometry/reduced_operator.py` DELETED; docs sweep `590c12d0` |
> | `7a0f434c` | **P4.9a** | the closure owns its march; DD spatial-only; `cell_balance_terms` dead; the visit family purely spatial; the handing; the degenerate frozen corpus; docs `ca852c44` (12 commits, `cb65c4cc`…) |
> | `3456dd37` | **P4.7** | the packet sheds ΔA/w+chord+mu (fusion: closure owns / cache interns / packet dead); the scan family renamed concept-complete; 4th+5th false docstrings fixed |
> | `2ec73b80` | **P4.6** | the moment mass consumes the AXES — `has_constant_volume_element` per axis; the chart tag leaves the family; `ndim` → `len(axes)`; mixed-tuple witness |
> | `260ddc64` | **P4.5** | the predicates say what they ask (`is_1d` re-keys; the annotation tells the truth; `_is_curvilinear` dead) — chain-scan reading; the title's premise dissolved at P4.1b |
> | `a60e5c0f` | **P4.9b** | the operator is POSED with its two closures (`.pose`; no defaults, no guards); the walk consumes the HANDED pair (keystone: post-pose hub swaps INERT); the strategy-layer intern (count-gated); the pole misnomer dead (`angular_closure`/`AngularClosureBase`); docs `9c3eb60a` (16 commits `b253732f`…`a60e5c0f` + tail `5a1591d2`, incl. archivist) |
> | `83c1ccc8` | **P4b** | the closure block has ONE home (closure caches five per-ordinate arrays read-only; cache 13→8 fields, builder takes no closure); dead `level_ordinates` retired; both strata refuse writes; the §4bis done-when REFUTED by the opener and re-ruled (3 commits `3ebb45d9`…`83c1ccc8`) |
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
> and `dead_references` 0 throughout. **P4.5 9847/0 — `[M]` 13 trees
> rc=0, delta 0 on EVERY axis vs the P4.9b baseline** (22 sk / 227 des /
> 70 xf; sn 3291, derivations 1637, numerics 2445, geometry 727, root 409
> all baseline-exact) — the chain-scan reading's "bit-identical by
> construction" measured true; pyright orpheus/ 0, tests/ −7 net.
> **P4.6 9848/0 — `[M]` 13 trees rc=0, delta exactly +1 = the new
> mixed-axes witness, PREDICTED before the run** (22 sk / 227 des / 70 xf
> unchanged; matrix.rst regenerated 10166 → 10167 by the same +1);
> pyright orpheus/ 0, tests/ at baseline; sphinx -W 0; dead_references
> 0/52. ⟹ the acceptance baseline for the next step is **9848**.
> **P4.7 9846/0 — `[M]` 13 trees rc=0, delta exactly −2 = the two
> retired-field smoke tests (their resolution claims live on in the
> `abs_mu` siblings), PREDICTED before the run** (22 sk / 227 des / 70 xf
> unchanged; matrix.rst 10167 → 10165 by the same −2); family battery
> 1240/0; pyright orpheus/ 0, tests/ baseline-exact 1386; sphinx -W 0;
> dead_references 0/52. ⟹ **the acceptance baseline is now 9846.**
> **P4b 9847/0 — `[M]` 13 trees rc=0, delta exactly +1 = the new
> read-only refusal witness, PREDICTED before the run** (22 sk / 227 des /
> 70 xf unchanged; the two gate re-poses are 1:1 so they move no count);
> sweep+operators pre-gate 2163/0; pyright orpheus/ 0, touched test files
> net −2; sphinx -W 0; dead_references 0/52. ⟹ **the acceptance baseline
> is now 9847.**
> **CS5 9888/0 at `b0bfc06c` — `[M]` 13 trees rc=0, delta exactly +41 =
> the new gate suite, PREDICTED before the run** (22 sk / 227 des / 70 xf
> unchanged), **then +8 measured on the post-fix re-run** (the roster's
> fifth factory `folded_product`, an archivist census finding: numerics
> 2484 → 2492, sn/angular unchanged) ⟹ **the acceptance baseline is now
> 9896** (= 9888 + 8, both legs measured at `f8c69117`). Teeth: 6-arm
> in-process mutation battery, every gate reddened by its named mutation
> (A→6, B→4, C→24, D→4, E→1, F→1), positive control 79/79 then 87/87;
> pyright orpheus/ 0 tree-wide; sphinx -E 0; dead_references 0/52.
>
> ### ▶ NEXT — **the campaign's next act is a SEQUENCING ruling** (P4.9 is
> COMPLETE: P4.9a ✅ `7a0f434c`, P4.9b ✅ merged 2026-08-29; **§5b ⛔ P0
> ✅✅ ANSWERED 2026-08-29** — ruled first at the post-P4.9b sequencing round
> and run the same day, probe-only). ⭐⭐ **P0's verdict: the 3-factor
> product form is REFUTED** — a theorem with a measured witness (the
> residual of `D·T_ang·T_sp` is EXACTLY `+G_ang·D⁻¹·G_sp`, ≤ 3.8e-16), and
> order-free (`L·T_spatial⁻¹` has 42–54 % of its mass outside the same-cell
> algebra, so no 3-factor form with that `T_spatial` exists in ANY order).
> **The measured algebra is a SUM of resolvents**
> `L = D − E_sp·S − E_ang·A_ang` (`[M]` rel ≤ 3.3e-16, 12 of 12 configs),
> each accumulator EXACTLY triangular in its OWN traversal order, and on
> the augmented space `[ψ̄; f; h]` `L` is the SCHUR COMPLEMENT of a sparse
> `Ã`. Answers in §5b's ✅✅ ANSWERED block; record
> `scratch/p0_product_form_measurement.md` (read its VERDICT block first).
> ⟹ **§5b's build is RE-POSED, not dead**: the adjoint still DERIVES
> (reverse each traversal independently — `[M]` which the shipped
> `loss_action_transpose` already implements), but the algebra a §5b
> campaign would declare is the sum-of-resolvents, not a product — a design
> fork for the user at the Campaign-2/O-3 boundary. ⭐ **The §7-remainder order was RULED
> 2026-08-29 (user): P4.5 → P4.6 → P4.7**, and **P4.5 is ✅ LANDED
> `260ddc64`** under the CHAIN-SCAN reading (the title's premise dissolved
> at P4.1b — see the section's banner; ground:
> `scratch/p4_5_ground_remeasure.md`), **P4.6 is ✅ LANDED `2ec73b80`**
> (ground `scratch/p4_6_ground_measure.md`), and **P4.7 is ✅ LANDED
> `3456dd37`** (ground `scratch/p4_7_ground_measure.md` — its opening
> re-measure REFUTED the `[R]` consequence note; the packet is four honest
> fields; the scan family speaks one name per concept). ⟹ **the ruled
> §7-remainder order is COMPLETE.** ✅ **RULED 2026-08-29 (user, the
> post-§7-remainder sequencing round): ▶ NEXT = P4b (the cache strata
> carve), with a context COMPACTION before execution.** P4b's opener
> owes the same shelf-life re-measure every §7-remainder phase paid:
> §4bis's strata census (`[M]` 2026-08-26: 7 fields bit-identical across
> 3 meshes / 4 mesh-bound / 2 traversal-only) predates P4.9b's intern,
> P4.5's `is_1d` re-key AND P4.7's renames + the ΔA/w
> formed-from-factors build — re-derive the per-field strata on the
> RENAMED cache before designing, and keep §4bis's ⚠ (do NOT convert
> `CollisionCache`'s eager build to lazy; its `_build_count` instruments
> a cardinal gate — now the P4.9b COUNT gate). P0's shape census
> corroborates the split ((N,) = T_ang's algebra, (N,nx) = D/T_spatial).
> ✅ **P4b is LANDED 2026-08-29** (`3ebb45d9` + `55d5783f` + `83c1ccc8`,
> ground `scratch/p4b_ground_measure.md` — the opener REFUTED the §4bis
> done-when and the phase was re-ruled *"Un-hoist to owners"*; full
> record + residue at §4bis's ✅/⛔ banners). ✅ **RULED 2026-08-29
> (user, post-P4b sequencing round): ▶ NEXT = CS5, with a context
> COMPACTION before execution.** CS5's HOME is
> `.claude/plans/space_and_kernel_binding_campaign.md` **§5.5** ("an
> axis can name the generator that made it") — resume THERE; this plan
> holds the tail CS5 gates (the P4-mint remainder + P3c's
> `ChartConnection` re-name). CS5's opener owes the campaign-pattern
> shelf-life re-measure (5 of 5 openers have corrected their own
> section): §5.5's `[M]` 2026-08-27 findings — the both-sides
> `generator↔space` census, the Axis loss site, the `mu_x`/`level_indices`
> 0-hit claims — predate P4.2/P4.3/P4.4/P4.9a/P4.9b/P4.5–P4.7/P4b, all
> of which moved SN mesh/space surfaces. `[M]` anchor drift at
> `01088fb5`: the angular-axis loss site is now
> `augmented_mesh.py:1175-1179` (§5.5 cites `:1170-1175`; same code);
> `Basis` ABC `numerics/basis/base.py:117` holds; `Axis`
> `numerics/axis.py:103`; `FrameBase` `numerics/frame.py:131`.
> Behind CS5: P6 (terminal sweep by §9.3's ruling) · the re-posed §5b
> build (Campaign-2/O-3 boundary fork). P5 rides O-3; P7 = #409. Read
> `scratch/p4_9b_design.md` §§7–10 before designing anything that
> touches the posing surfaces.
> ✅ **CS5's NODAL HALF IS LANDED 2026-08-29** (`4e7b8977` machinery +
> mint sites + gates; `b0bfc06c` the Protocol declares `level_structure`;
> `cb3cd15b` the archivist-census gate repairs; `f8c69117` the doctrine
> into `spaces.rst` §spaces-axis-generator + `discrete_measures.rst`) —
> full record + refutation banners at the §5.5 home; ground
> `scratch/cs5_ground_measure.md` (the 6th consecutive opener to correct
> its own section: [M] a DiscreteMeasure-typed accessor answers 3 of the
> done-when's 4 names — the angular generator is the QUADRATURE, measure
> ⊕ level fibration); triage `scratch/cs5_reach_past_triage.md`;
> verification plan `scratch/cs5_verification_plan.md`.
> ⟹ **P4's remainder is UNBLOCKED — the next act is a SEQUENCING ruling
> (user):** the P4-remainder phase (the producer binds to `(space, R)`
> reading the angular axis's generator; the §4ter WELD — the
> `AngularRedistribution.quadrature` courier, [M] read set
> `reduced_operator.py:448,517,528` + `closure.py:1659,2214` + 2 test
> sites — dissolves WITH the binding, and the `ChartConnection` name
> resolves as a consequence; §6c: gates **G5 (refusal) + G7 (route
> keystone)** land WITH that phase — specs ready at the verification
> plan §3.5/§3.7, the decoy mechanism pre-validated by mutation F)
> **versus P6 versus the re-posed §5b build.** ⚠ CS5's MODAL half
> (§5.5 items 2+3: the #409 non-diagonal metric + the polynomial
> `Basis`) stays open at the §5.5 home and does not block the
> P4-remainder.
>
> ⭐ **Before designing P4.9, read `scratch/p4_9_design_measured.md`** — the
> §6b pre-measurement (taken 2026-08-28 while P4.3's gate ran) re-sizes the
> charter's edit list: rows 1–3 (kill the twin, retire `cell_balance_terms`,
> shed the protocol's angular members) are the small coherent step; row 4 is
> `[M]` a **~150-call-site migration** (`StreamingOperator(sn_mesh)`, 1
> production + ~148 test sites in ~40 files) and row 5 touches `[M]` **~60
> production reads** of `SNMesh.scheme`/`pole_angular_closure` in 11 files
> (incl. an L2 reach-through at `transport/radial_characteristic_field.py:360`).
> ✅ The split was **RULED 2026-08-28 (user, design round)**: **P4.9a** =
> rows 1–3 + the `mm_a_in_coeff` handing + the `is`-identity gate; **P4.9b**
> = rows 4–5, its migration lever RULED a **production posing-head factory**.
> Full ruling block in the charter section.
>
> ⭐ **P4.9 was chartered 2026-08-28 (user)** — *each axis is closed by its
> own closure, and the operator composes them* — and its **P4.9a half is
> LANDED** (`7a0f434c`): #407 closed, the Morel–Montry twin at `diamond.py`
> is DEAD, the visit family is purely spatial, the closure mints its scan
> constants. **P4.9b** (the 4-arg ctor + `SNMesh` shedding) is the
> **prerequisite for §5b's ⛔ P0** (P0 must apply each factor separately).
> Before designing P4.9b: read the charter section's ruled blocks + landing
> audit, `scratch/p4_9_design_measured.md` rows 4–5 (the ~150-ctor-site +
> ~60-read migration), and §5b's FACTOR INVENTORY (its structural content
> feeds P0; its "realized today" line numbers are the PRE-carve tree's).
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
> | the α cluster is **five** symbols | ⛔ **SIX.** `[M]` `AngularRedistribution` annotates a field with `AngularMeasure`, so it had to travel too — leaving it behind is `geometry` naming an `sn` type. A grep for `alpha` finds **3 of 6**. |
> | `StreamingTerms` / the connection operator / the α cluster live in `geometry/` | ✅ **All moved.** `ReducedStreamingOperator` + the 3 factories → `sn/mesh/reduced_operator.py` (P4.4 `16501ca0`); the 6 α symbols → `sn/angular/redistribution.py` (P4.2 `5940deba`). `[M]` `geometry/reduced_operator.py` now holds **only `StreamingTerms`** and imports **nothing but `dataclasses`**. ✅ P4.3 `da507e3d` moved the packet to `transport/spatial/scheme.py:107` and DELETED the module. |
> | the L0 ladder imports `alpha_dome` from production | ✅ **REMEDIED by P4.2.** It ACCEPTS α as a keyword, like `tau`/`edges` since 2026-08-12. `[M]` `w` fed nothing but that call in all three graders, so it left `morel_montry_beta`'s signature and `alpha_defect_beta` lost `quad`/`geometry` too. ⛔ NOT the `WHITELIST` row. |
> | `transport/spatial/` is SN-only, so it should be `sn/spatial/` | ⛔ **PROPOSED AND WITHDRAWN 2026-08-28.** The import census measured MATURITY, not architecture — SN is the deliberate vanguard. `[R]` every method solving a cell-local balance forms a closure (FD diffusion IS box/DD, FE/DG IS LD's basis, LS-MoC needs it for the source); only MC is exempt. `transport/spatial/` is correctly placed. |
> | the spatial scheme legitimately closes the angular axis | ⛔ **NO — that is the twin.** `[M]` `DiamondDifference.update` evaluates the Morel–Montry relation inline, duplicating `closure.py:1327-1330`; both are LIVE on a data branch and **nothing gates them**. Forced by the layer (L2 cannot call L3). **P4.9** removes the obligation. ✅ **REMOVED at P4.9a `8a78be1d`** — the twin is dead, the is-identity gate + M8 anti-twin gate stand where nothing gated them. |
> | `AngularMeasure` "travels with the factories at P4.4" | ⛔ **REFUTED 2026-08-28.** The measurement (5 in-module use sites) was right; the inference was not. `[M]` `AngularRedistribution` NAMES `AngularMeasure` (`:1082`), so 2 of those 5 sites are in the α cluster. Sending it to `sn/` with the factories while α stays would be `geometry → sn` — forbidden. ⟹ it **stays in `geometry/` across P4.4** (factories read it as `sn → geometry`, legal) and **travels with α at P4.2**. Lands in `sn/` either way, so CS5 is not prejudged. |
> | `delta_A` moves to `sn/` as a field | ✅ **RULED 2026-08-27 (user): it DISSOLVES into R.** ΔA is R's rank-1 realization, not a second object. P4.1b is its first step — `delta_A` is already a derivation over `mesh.areas`, so what P4.4 moves is a derivation, not an array. |
> | `SNMesh.face_areas` / `.delta_A` are deprecated read-throughs | ✅ **REMEDIED by P4.1c** — retired. `[M]` 11 readers, 0 in `orpheus/`; every consumer was a test and 4 were the tests verifying the shims. |
>
> ### Still BLOCKED, with a named blocker (not a defer)
>
> Two of §4ter's forks — **the producer's home** and the **`ChartConnection`
> name** — resolve as consequences of **CS5**
> (`space_and_kernel_binding_campaign.md` §5.5, *an axis can name the generator
> that made it*), whose NODAL half is the prerequisite. Everything else is
> unblocked: **P4.5, P4.6, P4.7 and P4b** (P4.2, P4.4 and P4.3 ✅ landed).
> ⚠ **P3 is NOT finished** — `ReducedStreamingOperator` → `ChartConnection`
> (P3c) is deliberately sequenced **after P4**, so the name describes what
> remains once `streaming_terms` is out of it.
>
> ### ⭐ The rules P4.4 + P4.2 + P4.3 + P4.9a + P4.9b paid for — the ones that generalise
>
> **+2 from the P0 → P4.5 → P4.6 → P4.7 stretch (2026-08-29):**
> - **A phase OPENS with the shelf-life re-measure of its own section, as a
>   dispatched ground memo — and three of three openers were paid.** P4.5's
>   re-measure found the section's title premise dissolved by P4.1b and its
>   "12 of 36" stale (the fork the design round then ruled); P4.6's found the
>   operand home the section never named (`Axis1D.coord`) and the hub already
>   binding `axes`; P4.7's REFUTED the plan's own `[R]` consequence note and
>   re-counted the migrations (23 → 27 + ~57 uncounted). Cost of each memo:
>   one explorer dispatch during the previous phase's gate; value: zero
>   designs built to stale numbers.
> - **A full-gate reading is only evidence when its DELTA was predicted
>   before the run** (the §10 discipline applied to the suite count): P4.5
>   predicted 0-on-every-axis, P4.6 predicted +1, P4.7 predicted −2 — each
>   confirmed exactly, so each gate certifies the change rather than merely
>   ending green. An unexplained delta of even ±1 is a finding, not noise.
>
> - **P4.9b — a hazard CLAIM about runtime behaviour is run, not reasoned**
>   (F12): my "the wrong-family closure is silent, plausible-wrong k" was
>   refuted by EXECUTING the doctored state — it raises on every geometry
>   where it could matter. §6d's inject-and-run, generalised from import
>   cycles to any claimed failure mode; the false sentence was one commit
>   from living in a ctor docstring as licence for the guard the ruling
>   forbids.
> - **P4.9b — a no-guard ruling needs its enforcement MEASURED**: mutate
>   the by-construction path (M5: `.pose` mints instead of reads) and
>   record that ONLY structural legs redden while every value gate stays
>   green — that measurement, not the argument, is what the ruling rests
>   on going forward (`p4_9b_verification_plan.md` §9).
> - **P4.9b — a perf consequence is a COUNT, never a percentage**: the
>   operator-count behind F2's 24.65 % proved fixture-dependent (6–10 vs
>   38–43/solve between two honest measurers); the shipped instrument
>   pins builds==1 exactly. §4's configuration clause at instrument-choice
>   scale.
> - **P4.9b — a §6b census over literal names must declare that predicate**
>   (the surprise-log row): registry-variable calls, the changed class's
>   own internal call, and monkeypatch surrogates are member spellings no
>   name-keyed AST census can see; the red loop is their enumerator and
>   must be budgeted as such.

>
> - **P4.9a — an acceptance ARTIFACT must be qualified by ACTIVATION, not
>   name** (`plan-authoring` §10's new log row): the charter's canary
>   provably could not execute the carved path (`n_phi ≡ 2 (mod 4)` only;
>   the whole frozen corpus at 4/8), so its bit-identity was unconditional.
>   Ask *"if the carve broke, would this artifact move?"* — and the anchor
>   that qualifies must land PRE-carve on unmodified production.
> - **P4.9a — a stochastic measurement's configuration includes its DRAW**
>   (`plan-authoring` §4's draw sharpening): the 59 %/204-ULP two-forms
>   figure was one seed's reading (46–51 % / 1e2–1e5 ULP across 200); only
>   `max |Δ| = 1.776e-15` is stable, and the weld gate pins THAT — a nulp
>   band would have certified a seed.
> - **P4.9a — the phase-proof shape for any twin-kill: the owner-mutation
>   red set must be a strict SUPERSET of the pre-carve twin's catcher set**
>   (13 → 50 here), measured with the same battery before and after —
>   else the carve was cosmetic.
>
> - **P4.3 — a `:mod:` role in a `#` COMMENT is invisible to every instrument
>   except grep.** Sphinx never renders comments, `dead_references` reads
>   docstrings, `-W` is silent — so P4.4's clean docstring sweep left
>   `augmented_mesh`'s Wave-B comment saying the math "lives in geometry",
>   and only P4.3's concept-grep caught it. ⟹ the spelling-surface list is:
>   imports, docstrings, `.rst`, strings, **and comments**.
> - **P4.3 — a HOME is ratified per-field, by ownership ontology, not by any
>   location argument.** The user refused both offered file options until
>   every field was classified by what it IS (mesh's / quadrature's / the
>   weld's); the home then followed from *which contract declares the type*
>   (`CellVisit.streaming_terms` ⟹ beside `DiscretizationSchemeBase`). The
>   vanguard/contract lesson, extended from packages to types — and it also
>   produced the direction-supplier ruling as a by-product.
>
> - `plan-authoring` **§6d, INTRA-FILE clause** (P4.4): **a re-home's forbidden
>   edge can be INVISIBLE TO AN IMPORT GREP, because the split CREATES it.**
>   §6d says *enumerate the mover's callers* and every reading reaches for
>   imports — but a caller in the mover's OWN FILE imports nothing. That is what
>   made "P4.2 before P4.4" wrong after "P4.3 before P4.4" was already corrected.
> - `coding-standards` (P4.4): **a COMPLETE import audit is still a PARTIAL
>   audit — a symbol has more than one spelling SURFACE.** AST-enumerated
>   imports + a positive-control residual check reading **0** + green suite +
>   pyright 0 + `sphinx -W` clean, and `dead_references` still found **9 live
>   refs in DOCSTRINGS**. The filter was validated and pointed at the wrong
>   corpus. ⟹ **run `dead_references` before calling any re-home done** (P4.2 ran
>   it proactively and it caught one the author had just written).
> - `plan-authoring` surprise log (P4.3 design): **an import census in a codebase
>   with a deliberate VANGUARD measures MATURITY, not architecture.** I measured
>   `transport/spatial` SN-only against a real control and proposed moving it;
>   the user refuted the METHOD — SN is deliberately first, so "no importer"
>   means "not built out yet". ⟹ a placement argument must be a **CONTRACT**
>   argument; a consumer census may only corroborate one.
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
> re-signing three graders cost that tree nothing). ⟹ ✅ **P4.3 CONFIRMED
> the −1: `9824 / 0, 22 sk / 227 des / 70 xf`, localized to root+harness
> (layer gate 344 → 343), all other trees and axes +0** — the first phase of
> this arc where the count went DOWN, and it went down because the file died.
> ✅ **P4.9b: `9847 / 0, 22 sk / 227 des / 70 xf` — `[M]` 13 trees `rc=0`
2026-08-29, ~62 min; delta vs P4.9a exactly +11 = the phase's new gates
(3 witnesses + 4 keystone rows + 2 read-set rows + 2 cache gates); per-tree:
sn 3291 (+11), every other tree identical (numerics 2445, transport 566,
geometry 727, derivations 1637, root 409, …).**
✅ P4.9a: `9836 / 0, 22 sk / 227 des / 70 xf` — `[M]` 13 trees `rc=0`,
> 64 min, delta exactly +12 passed and 0 on the other three axes,
> localized to transport +1 (the anti-twin gate) / sn +11 (2 is-identity +
> 3 minted-constants + 2 cache-handing gates, 1 snapshot + 1 walk-baseline
> + 2 CYL_DEG affine rows), every other tree +0 on every axis. No new
> production module ⟹ the layer gate stayed at 343.**
> **The acceptance baseline for the next step is `9836 passed / 0 failed,
> 22 skipped / 227 deselected / 70 xfailed`.**

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
| `#407` — DD's blend inverse in scheme-neutral code | `2.0 = 1/w_DD` at `cell_balance.py:246/:341` (post-P4.3 numbering) | `[M]` `cell_balance_terms` has ONE production consumer (`diamond.py:194`); `linear_discontinuous.py` never imports either helper; and the matvec **declares** the specificity the producer does not (`loss_representation:3091`: *"Curvilinear (DD-only) … with DD's diamond march inlined"*). | The whole module is a DD body in a scheme-neutral name and home. Carve with O-3. |
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

> ✅ **P4b (2026-08-29) resolved the CACHE-side face of this weld** — the
> ownership carve, not a lifetime carve: the closure block's twin storage is
> gone (one home per datum), but the closure itself is STILL minted per mesh.
> `[M]` (`scratch/p4b_ground_measure.md` §B.3) the closure is not
> (quad, coord)-pure as an object — it holds `_gram` `(nx,1,1)` and
> `_dAw_per_level` `(nx, M_p)`, consumed by its own kernels — so no
> (quad, coord)-lifetime owner exists on this path today, and hoisting one
> out requires splitting the closure's own two strata. **That carve rides
> O-3** (where the closure/factory family splits anyway; Campaign 2
> territory), recorded here so the F3 sentence above cannot be read as an
> unaddressed defect of the cache.

⚠ **What NOT to make lazy.** `CollisionCache`'s eagerness is load-bearing —
`_build_count` instruments a cardinal invariance gate (*exactly once per Σ_t
epoch*). Do not convert a gate-instrumented eager build into a lazy one; the
gate measures the build.

---

## 4bis. ⛔ The cache class welds THREE strata — measured, and it refuted both proposed names

> ✅ **REMEDIED by P4b, 2026-08-29 (`3ebb45d9` + `55d5783f` + `83c1ccc8`).**
> The weld below is HISTORY: the closure block (4 fields) has one home on the
> closure, `level_ordinates` (production-dead, `[M]` zero readers of any
> kind) is retired, the cache is the 8-field chain-ordered spatial table, and
> both strata's arrays refuse writes. The phase's opening re-measure
> (`scratch/p4b_ground_measure.md` — the 4th consecutive opener to correct
> its own section) REPRODUCED the S0-seven verdict on all three charts and
> REFINED it: the S1 quartet is chart-dependent (on the SLAB only `volume`
> is edge-sensitive; the other three are neutral constants behaving as S3),
> a denominator the census below never stated.

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

> ⛔ **The done-when above was REFUTED by the phase's own ground measure
> (2026-08-29, `scratch/p4b_ground_measure.md` §A.0.3/§A.4/§F5) and the
> phase was re-ruled** (user, P4b design round: *"Un-hoist to owners"* over
> the literal split):
> - *"`Quadrature × CoordSystem` **alone**"* cannot be an honest key —
>   `[M]` 4 of the 7 S0 fields are CLOSURE-minted, the closure class is
>   user-overridable per mesh, and two shipped gates (the intern's identity
>   validation; the keystone's activation leg) require a doctored closure
>   handed at pose to get its own values.
> - The identity criterion had no owner to live on — `[M]` no
>   (quad, coord)-lifetime object exists on this path (closure,
>   `AngularRedistribution`, reduced operator are all per-mesh instances),
>   and the shareable win is O(N) trivial allocations, unmeasured — cache
>   machinery against the algebra-eager/performance-lazy rule.
> - *"`chain_idx` moves to the traversal layer"* died with the same
>   measure: the chain-ordered tensors are STORED BY that permutation —
>   table and permutation are one strategy artifact (the P4.9b binding
>   principle's fused table), and no consumer asks for them apart.
>
> ✅ **LANDED 2026-08-29 as the ownership carve** (`3ebb45d9` retire dead
> `level_ordinates` + sweep-docstring truth; `55d5783f` the shed — closure
> caches all five per-ordinate arrays read-only, cache 12→8 fields, builder
> takes no closure, walks read their HANDED closure, `from_geometry` handed
> the closure, fidelity gate re-posed as the accessor-stability witness;
> `83c1ccc8` both builders' arrays refuse writes + witness). **Done-when as
> landed:** grep `geom.c_in|c_out|tau_inv|mm_a_in_coeff` → `[M]` 0 hits
> tree-wide; the closure's five per-ordinate accessors identity-stable +
> read-only (witnessed); `level_ordinates` unspellable; `dead_references`
> 0/52; bit-identical (sweep+operators trees 2163/0 pre-gate; full gate
> `[M]` **9847/0, 13 trees rc=0, delta exactly the predicted +1** — the
> 4th consecutive exact pre-run prediction).
> ⚠ Residue, recorded not deferred: `geometry_cache_for`'s closure-identity
> validation is now VESTIGIAL cache-side (a doctored closure's rebuild is
> bit-identical — the walk reads the doctored closure directly), retained
> as the intern's declared key; it dissolves with the strategy layer at
> Campaign 2 (Q1's own rationale). The closure's own (quad,coord)×mesh
> weld → §4's O-3 rider note.

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

✅ **ALL SIX ROWS EXECUTED as ruled, 2026-08-27/28** — `face_areas`/`coord`
retired (P4.1a/b), the α symbols + `AngularMeasure` to
`sn/angular/redistribution.py` (P4.2 `5940deba`), `StreamingTerms` to
`transport/spatial/scheme.py:107` (P4.3 `da507e3d`), `delta_A` a derivation
on the operator (P4.1b), and **the module is deleted**. `[M]` 2026-08-28:
14 files name the primitive family, zero under `orpheus/geometry/`.

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
`orpheus/sn/angular/redistribution.py`; `geometry/reduced_operator.py` was
left holding `StreamingTerms` alone, `[M]` importing **nothing but
`dataclasses`** — a file in `geometry/` that named no geometry, which is the
whole finding the arc was built on. ✅ P4.3 `da507e3d` deleted it.

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
reads `mesh.areas`); no new edge either way.~~

⛔⛔ **The DESTINATION's reason was wrong, and the destination is still right —
corrected 2026-08-28 after a user challenge.** A pre-execution memo argued
`transport/spatial/scheme.py` from **where the consumers are**, and separately
proposed moving the whole `transport/spatial/` subpackage to `sn/spatial/` on an
import census (`[M]` 10 consumer files, all `sn/`, 0 intra-transport, 0 other
methods, against a control showing every genuinely shared transport subpackage
has 2–3 consumer packages plus intra-transport use).

**Both inferences are void.** Verbatim (user): *"the discretization scheme was
born in Sn first and moved to transport, because we considered it useful for all
transport methods that need to build a spatial discretization. Measuring no
import today is a signal that the other transport methods are in infancy stage
… we're using Sn specifically to bring the entire machinery to state of the art
before using this machinery in the other methods. So don't assume by import."*

⟹ **an import census in a codebase with a deliberate VANGUARD measures maturity,
not architecture**, and it always argues for narrowing what was deliberately
built wide. (→ `plan-authoring` surprise log, 2026-08-28.)

`[R]` **The first-principles argument, which is the one that decides it.** A
cell balance is ONE equation in (ψ_in, ψ_out, ψ̄) — under-determined — so every
method solving a cell-local balance must posit how ψ varies inside the cell:
FD diffusion **is** the box/DD scheme; FE/DG diffusion **is** LD's basis; nodal
is a transverse-integrated expansion; LS-MoC needs it for the **source**;
CP-with-linear-source likewise. Only MC is exempt. ⟹ `transport/spatial/` is
correctly placed and the `sn/spatial/` proposal is **withdrawn**.

⭐ And `StreamingTerms` field-by-field: `[M]` exactly **one** of its seven fields
is SN-specific — `delta_A_over_w`, whose divisor is a **quadrature weight**.
`chord_length` / `face_area_inner` / `face_area_outer` / `volume` are any
finite-volume method's; `mu` / `abs_mu` belong to any **directional** method.
And P4.7 retires the one, so the packet becomes method-agnostic by a step
already scheduled. ⟹ the destination stands on the **contract** argument (fork
3's own ruling: a consumer-count argument flips, a contract argument does not).

⛔ **TWO FIELDS ARE PRODUCTION-DEAD — a P1 remainder P1 missed.** `[M]`
2026-08-28, two independent filters (AST scoped to the 13 packet-naming files,
and a direct grep over all of `orpheus/`):

| field | production readers | test readers | constructed at |
|---|---|---|---|
| `chord_length` | **0** | 15 | `sn/mesh/reduced_operator.py:534` (was `:532` pre-P4.3 docstring edit) |
| `mu` | **0** | 8 | `:535` |

Production reads `abs_mu` and **never** `mu`. P1's charter was *"nothing in the
streaming path is held that nothing reads"*, so these two are exactly its
subject and it did not reach them. ⟹ **their retirement is scheduled into P4.7**
(the phase that already owns `StreamingTerms` field retirement); it is NOT
P4.3's — P4.3 moves the packet, it does not change its shape. ⚠ 23 test readers
migrate with them (`coding-standards`: retiring means test migration).
 **REFUTED 2026-08-27 — the
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

✅ **P4.3 LANDED `da507e3d` (+ docs sweep `590c12d0`) 2026-08-28.** Exit gate:
**9824 passed / 0 failed, 22 skipped / 227 deselected / 70 xfailed** — the
predicted **−1**, delta **0 on every other axis**, localized to root+harness
(410 → 409: the import-linter's `rglob` losing the deleted module, layer gate
344 → 343), all 12 other trees +0; 13 trees `rc=0`, 62 min wall. pyright **0**;
`sphinx -W` **0**; `dead_references` **0 dead / 52 checked** — ⭐ the first
phase of this arc where the docstring-surface sweep was complete BEFORE the
instrument ran (P4.4 and P4.2 each had 9 found). Fresh-interpreter smoke with
`orpheus.geometry` imported FIRST. The verification matrix (build-regenerated,
not hand-written) independently confirms 10144 → 10143.

**The home was ratified through an OWNERSHIP evaluation, not a location
argument** (user, at execution: *"evaluate ownership of the terms and propose
their home based on ontology reasoning"*). `[M]` field-by-field, AST
receiver-attributed (`scratch/p4_3_design_measured.md` §"Ownership by
ontology"): the metric triple (`face_area_inner`/`_outer`, `volume`) is the
MESH's, single-sourced there; `abs_mu` is the QUADRATURE's, consumed as the
closure's direction parameter — both members read it, and
direction-parameterized ≠ angular-closure-aware; `delta_A_over_w` is the weld,
owned by neither (P4.7); `chord_length`/`mu` production-dead (P4.7). ⟹ the
packet is the **evaluation point** of a spatial closure — (cell metric) ×
(direction coordinate) — and the TYPE is the family contract's argument
vocabulary, consumed by BOTH members only through
`CellVisit.streaming_terms`. A type declared by a contract lives with the
contract: `scheme.py` (`StreamingTerms` now at `:107`), module-level sibling
of `DiscretizationSchemeBase` — the "base both DD and LD subclass from" the
user hypothesized; it exists, and both do. ⚠ LD's non-read of the face areas
is MATURITY (it refuses curvilinear; slab A ≡ 1), not ontology — the vanguard
rule forbids reading that 0 as "LD doesn't need them".

**Stale AT HEAD, found by this phase's sweep and fixed in passing** (§3 tense
cases, each): the `structured_geometry.rst` enumeration was pre-P4.2 — re-run
at its own predicate: **14 files, 3 definers, and for the first time ZERO
under `orpheus/geometry/`**, the arc's done-when now a property of the
enumeration; `augmented_mesh`'s Wave-B **comment** still said the math "lives
in geometry" — ⭐ a `:mod:` role in a `#` comment is invisible to Sphinx AND
`dead_references` (which reads docstrings); only grep sees that surface;
index.rst's "geometry-layer … geometry-only"; two "purely geometric"
parentheticals (index.rst Step-C note, history.rst changelog);
`test_diamond`'s fixture docstring naming α fields the packet lost at #236
Step C; "namedtuple" for a frozen dataclass.

**Recorded seam, NOT built:** when diffusion's box scheme joins the family
(R19), it consumes the metric triple with no ordinate — the {CellMetric} ×
directional-extension carve is O-3's decomposition-by-destination, cheap at
the family site. And the NAME `StreamingTerms` is §9.1's own defect
("streaming is what `L` does") — a fork for P4.7-time, when the shape
settles. ⚠ `CellBalanceTerms` is TAKEN (`cell_balance.py`).

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

**P4.5 — the hollow Cartesian object dies.**

✅✅ **EXECUTED 2026-08-29 under the CHAIN-SCAN reading @ `260ddc64`**
(user-ruled at the design round; ground re-measured FIRST —
`scratch/p4_5_ground_remeasure.md`, whose per-claim verdict table governs
every number below). ⛔ **This section's TITLE names a premise P4.1b
dissolved one day after the section's numbers were taken**: the slab object
is no longer hollow — P4.1b ("the slab is not a special case") made it the
zero-curvature case of the ONE chain-scan carrier (real widths/volumes,
neutral angular element), and `[M]` at HEAD `reduced is not None ⟺ is_1d`
bit-exactly, so the charter sentence ("exists iff there is a chain scan")
was ALREADY TRUE with zero existence change. The ruled repair is the
honest-predicate package, bit-identical: the three admission raises
(`dag_walk`, `dag_walk_cell_indices`, the scan-cache builder) + the
scan-path dispatch key on `is_1d` (the QUESTION) instead of
`reduced is None` (the realization); the annotation becomes
`ReducedStreamingOperator | None` (the d≥2 type-ignore dropped; the
doubly-false mint comment replaced by the presence contract);
`_is_curvilinear` retired; 35 test-side Optional narrowings ride along
(tests/ nets −7 pyright errors, orpheus/ stays 0). The CURVATURE reading
(the slab loses the object) was REJECTED — it would re-introduce the
special case P4.1b deleted.

`[M]` d≥2 Cartesian gets
`reduced = None`; the slab gets an object whose every meaningful field is
`None`/neutral, minted *"so `sn_mesh.reduced` is always populated"* — a promise
the d≥2 arm already breaks. `reduced is None` today conflates **"no chain scan"**
(`ndim ≥ 2`) with **"no curvature"** (`coord is CARTESIAN`), and `SNMesh` answers
both directly. `[M]` **12 of 36** code reads of `.reduced` are `None`-guards
paying for the conflation. (⛔ STALE at execution — the 12 reproduces at the
plan-date tree as 12 of 31 production Loads; the 36 reproduces under no
predicate the re-measure ran; at HEAD the honest fraction was 11 of 17 —
memo §6.) The object should exist **iff there is a radial
reduction**. (⟹ under the chain-scan reading this sentence was ALREADY TRUE
at HEAD — see the banner above.)

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
  not on `reduced is None`. (✅ DISCHARGED in two steps: the raise landed at
  P4.9b step 2c with its 2-D witness; the `is_1d` re-key landed at P4.5
  `260ddc64`.)
* `sn/solver.py:2316-2324` — `_is_curvilinear` reads
  `getattr(mesh, "coord", None)` then `getattr(coord, "name", str(coord))`
  *after* an `isinstance(mesh, Mesh1D)` narrowing. `coord` is a required
  dataclass field on `Mesh1D`, so **both defaults are unreachable**, and the
  function then compares an upper-cased *string* against `("SPHERICAL",
  "CYLINDRICAL")` rather than comparing the enum. It is the stringly-typed
  spelling of `mesh.coord is not CoordSystem.CARTESIAN`, which `SNMesh` already
  answers as `is_cartesian`. ⚠ Not a P4.1a site (it reads the *mesh*, which is
  the successor), and not a bug today — but it is a third spelling of this
  step's own predicate, so it retires with the other two. (✅ RETIRED at P4.5
  `260ddc64` — and the re-measure sharpened it: `[M]` ZERO callers since birth
  2026-05-09, so the retirement was pure deletion.)

**P4.6 — the moment mass consumes the MEASURE, not a chart tag.**
✅ RULED 2026-08-27 (user): *"it is very important to correct the flaw"*, in the
move, after compaction.

✅✅ **EXECUTED 2026-08-29 @ `2ec73b80`** (operand RULED at the design round:
**per-axis**; ground `scratch/p4_6_ground_measure.md` — every claim below
CONFIRMED at HEAD, the tag surface exactly 4 definitions with 2 plumbing-free
production callers). Realization: the family consumes `tuple[Axis1D, ...]`
and asks each axis the NEW protocol capability
`has_constant_volume_element` (AxisMesh True / RadialAxisMesh False) — the
measure's behaviour read polymorphically at the granularity the Kronecker'd
mass is built at; `ndim` died into `len(axes)`; whole-mesh `CoordSystem`
left the signatures and both files' runtime imports; the guard names the
offending AXIS kinds and a mixed (z, r) tuple is admissible + witnessed
(the case the whole-mesh enum structurally could not pose). ⚠ Exactly as
ruled, the REFUSAL stays — the curvilinear multi-moment VALUE still needs
#409 + #158; what changed is who is asked. The transitional-tag admonition
carries its resolution note; the full M+R mint remains CS5's P4-mint.
P4-item-1's tag scope is thereby DISCHARGED for this family (`[M]` the
explorer's census: no other member anywhere).

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

✅✅ **EXECUTED 2026-08-29 @ `3456dd37`** (ground
`scratch/p4_7_ground_measure.md`, taken FIRST — and it REFUTED the P4.9
ordering block's `[R]` "P4.7 became a consequence": `delta_A_over_w` had
TWO production readers, the cache builder AND the walk's degenerate
cylinder arm, the second by P4.9a Q2's own bit-identity ruling). What
landed: (1) the packet sheds `delta_A_over_w`/`chord_length`/`mu` — the
first was one of **THREE** stored spellings of the ΔA/w fusion (the
2026-08-27 ruling's 2-row table predates P4.9a/b minting the closure's
`_dAw_per_level`): the CLOSURE owns the fusion (algebra), the CACHE
interns the strategy-side copy formed from its two factors at build,
the packet's per-(cell,direction) rebuild is DEAD; both production
readers re-sourced bit-identically from the factors; ~90 test readers
migrated (⛔ the section's "23 test readers" undercounted — `[M]` `mu`
had 12, not 8; and the ΔA/w reader population was ~57 lines the section
never counted). (2) The rename family landed CONCEPT-complete per the
2026-08-26 lesson — cache fields AND the scan-family kw-only params AND
the two string-spelled consumers together (310 word-bounded + surgical
`V`; `numerics.vector.V` excluded): `delta_A_over_w` /
`face_area_downstream` (the predicted Pattern-2 unification — CellVisit
already owned the name) / `face_area_total` / `volume`. (3) The 4th and
a newly-found 5th present-tense-false docstring rewritten
(`_weight_of`'s "neither side owned the fusion"; `slab_streaming`'s
"stay ``None``", false since P4.1a/b); the theory pages' walk-idiom
examples updated; the error catalog keeps its historical spellings.
Section (c) stays a RECORD — P4.8 routes the axis-measure gap to CS5.

**(a) The fusion.** *"A fusion is accepted to be cached as a performance
optimization at the place that will assemble a hot path … we're not writing
algebra to conform to performance optimization and welding things."* `[M]` the
two fusion sites are NOT alike: (⛔ CORRECTED at execution — by then the
tree held **three** spellings, P4.9a/b having minted the closure's own
`_dAw_per_level`; the honest 3-row statement is closure = OWNER (algebra),
cache = strategy INTERN (this row's KEEP-renamed), packet = the dead
rebuild (this row's RETIRE). The two rows below are the 2026-08-27
population.)
* `StreamingCoefficientCache.dA_w` — an `(N, nx)` array built **once** at solver
  init. That IS pre-operating at the hot-path assembly point. **KEEP**, renamed.
⛔ **P4.7 also retires the two PRODUCTION-DEAD fields** `[M]` 2026-08-28:
`StreamingTerms.chord_length` (**0** production readers, 15 test) and
`StreamingTerms.mu` (**0** production readers, 8 test — production reads
`abs_mu` and never `mu`). Both are constructed at
`sn/mesh/reduced_operator.py:534-535` and read by nobody in `orpheus/`. They are
a P1 remainder — P1's charter was *"nothing is held that nothing reads"* and it
did not reach them. Full measurement + the two-filter method: P4.3's section.
⚠ 23 test readers migrate. (⛔ `[M]` at execution: 27 for chord/mu —
`mu` had 12, not 8 — plus ~57 ΔA/w reader lines never counted here; all
migrated at `3456dd37`.)

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
the nodes. `[M]` the exact loss site is `augmented_mesh.py:1171-1176`, which
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
nodes. Loss site: `augmented_mesh.py:1171-1176`, which builds the angular axis
from `quad.N` + `quad.weights` alone while declaring `kind=BasisKind.NODAL`.
`[M]` `mu_x`/`eta`/`mu_z`/`level_indices` → **0 hits** in BOTH
`numerics/space.py` and `numerics/axis.py`.

**(b) ⭐ But a MODAL axis is generated by a BASIS, not a measure** — so the rule
generalises rather than excepting. `[M]` `FrameBase`'s fields are literally
**`['basis', 'measure']`**, with `analysis` (nodal→modal), `reconstruction`
(modal→nodal), `discrete_gram`, `basis_space` and `measure_space`. And the one
shipped MODAL axis confirms it: `scheme.py:1716` sets
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

⛔ **P0 ANSWERED 2026-08-29 — the product form above is REFUTED; the
section's GOAL survives on a different algebra.** The display
`L = D @ T_ang @ T_spatial` fails at 62–136 % relative on every curvilinear
chart, and the failure is a THEOREM with a measured witness: the residual is
exactly `+G_ang·D⁻¹·G_sp`, the second-order angular×spatial cross term a
product mints and the assembled SUM does not contain (`[M]` identity verified
≤ 3.8e-16). What IS true (`[M]` 12 of 12 configs, record
`scratch/p0_product_form_measurement.md`): `L = D − E_sp·S − E_ang·A_ang` at
≤ 3.3e-16 — three SUMMANDS; each accumulator EXACTLY triangular in its own
traversal order but a RESOLVENT of a sparse shift (`S = 2(I+N_sp)⁻¹N_sp`
bit-exactly — a dense lower triangle), and on the augmented space
`[ψ̄; f; h]` `L` is the SCHUR COMPLEMENT of a sparse `Ã` (≤ 2.0e-16). The
payoff SURVIVES re-posed: `Lᵀ` = reverse each traversal independently —
which is precisely what the shipped `loss_action_transpose` does — so a §5b
campaign declares the sum-of-resolvents algebra (each term carrying `.H`),
not a product of triangular factors.

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

✅✅ **ANSWERED 2026-08-29** (`numerics-investigator` probe, 9 scripts
`scratch/p0_probe_*.py` + `p0_lib.py`; memo
`scratch/p0_product_form_measurement.md` — read its VERDICT block first; the
subject was tied back to the `apply_transpose ≡ apply.T` gate's own dense
harness, and 6 positive controls all bite):

1. **Form/order — REFUTED, all orders.** 9 candidate forms (5 3-factor
   orders + 4 diagonal-absorbed 2-factor forms) rejected at 42–136 %
   relative on 10 curvilinear configs; the slab control passes exactly
   (both cross terms ≡ 0.000000e+00, `E_ang ≡ 0` there). The measured truth
   is the SUM `L = D − E_sp·S − E_ang·A_ang`, rel 3.4e-17…3.3e-16, `[M]`
   12 of 12 configs (sphere GL 4–14 / cylinder fp(2..4, 6..10) / slab;
   N odd at 9 and 15, nx ∈ {4,5,7}; non-uniform widths, heterogeneous Σ_t;
   ng=1 — the cache is ng-free by construction, `cache.py:183`).
2. **Arity — 3 SUMMANDS, not 3 (or 4) product factors.** `P_mirror` is NOT
   a factor: the pole continuation is `[M]` ONE edge of `N_sp` per chain,
   absorbed by chain CONCATENATION. On the augmented space `[ψ̄; f; h]` the
   operator is sparse (1.4–5.7 %) and `L` is its SCHUR COMPLEMENT (`[M]`
   ≤ 2.0e-16, 6 of 6); `Ã` itself block-LU-factors (≤ 9.7e-18), but that
   factorization's Schur block either omits `L` or IS `L` — circular as a
   construction of `L`.
3. **Triangularity — YES, exactly, each in its OWN order** (`[M]`
   max|strict-upper| ≡ 0.000000e+00, 12 of 12; angular order =
   (cell, level, WITHIN-position) — not global ordinate index; spatial
   order = (chain class, chain position), the leg-derived chain matching
   the cache's `chain_idx` on every non-degenerate ordinate). But the
   factors are RESOLVENTS of sparse shifts — `S = 2(I+N_sp)⁻¹N_sp`
   BIT-exactly — i.e. dense lower triangles with sub-diagonal reach =
   chain length − 1, NOT the nearest-neighbour bands the inventory table
   reads as (see the ⛔ C-1 correction at the inventory). ⚠ Denominator:
   GL-sphere + folded-cylinder rules. On a Gauss-Lobatto sphere rule
   `edge_extrapolated_seed` FIRES (the tree's own census,
   `closure.py:1805-1830`, #361) and triangularity there is UNMEASURED —
   **#415** carries the witness item.

⚠ The paragraph above guessed that a refutation would mean *"the adjoint is
genuinely irreducible … and the twin must stay"* — HALF wrong, in the good
direction: the product is refuted AND the adjoint still derives, from the
sum (per-summand traversal reversal, which `loss_action_transpose` already
implements). Ride-alongs of this record: `test_282`'s "non-carrying control"
docstrings corrected (its cylinder is a second CARRYING case since
`384d62e4`), the corpus-wide sweep of the same staleness, and **#415**.

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

### ✅ THE MECHANISM — the operator is the composition site (user, 2026-08-28)

**Ruled shape.** `StreamingOperator` is constructed from **a domain, a codomain,
a spatial closure and an angular closure**. Discretization schemes are
**factories returning a spatial closure**; the angular scheme returns an
**angular closure**; those four arguments fully determine a well-posed operator
with a derived adjoint.

⭐⭐ **This is not speculative — the SAME unweld already SHIPPED for the
vectorized path, and only the per-cell path is the holdout.** `[M]`
`transport/spatial/cell_balance.py` carries **two** balance helpers:

| helper | signature | closure-blind? | consumer |
|---|---|---|---|
| `cell_balance_for_streaming` | `(*, abs_mu, A_downstream, A_total, total_xs, volume, psi_face_in, **angular_denom_term**, **angular_numer_upstream**)` | ✅ **YES — injected** | the vectorized matvec |
| `cell_balance_terms` | `(st, A_downstream, total_xs, upstream_state, *, **c_in**, **c_out**)` | ⛔ no — names M-M, reads `upstream_state.angular_upstream` | `[M]` **one** production consumer: `diamond.py:194` (#407) |

The blind one's own docstring records the move: *"Phase 2.11 pushes those names
into the `PoleAngularClosureBase` strategy: the matvec body calls
`closure.cell_contribution(...)` to obtain `(angular_denom_term,
angular_numer_upstream)` … `cell_balance_for_streaming` is now **geometry-blind
by interface** (no M-M names)."* ⟹ the mechanism is **precedented in the same
file**; what remains is to finish it on the per-cell arm.

#### Why the layer works — and it is the ONLY layer that does

`orpheus/sn/operators/streaming.py` is **L3**. `L3 → L2` (`transport/spatial`)
is legal; `L3 → L3` (`sn/angular`) is the same package; and **L2 never names
L3**. The composition happens at the only layer permitted to see both factors,
so the twin disappears because the OBLIGATION disappears — not because the
relation was injected into L2.

#### `[M]` the ingredients already exist — they are bound to the WRONG object

| ingredient | today | under the mechanism |
|---|---|---|
| spatial closure | `SNMesh.scheme` (`augmented_mesh.py:268`) | a constructor argument of the operator |
| angular closure | `SNMesh.pole_angular_closure` (`:422`, built from a class) | a constructor argument of the operator |
| domain / codomain | derived **properties** of `StreamingOperator` | constructor arguments |
| the operator itself | `[M]` `@dataclass class StreamingOperator` with **exactly ONE field, `sn_mesh`** | `(domain, codomain, spatial_closure, angular_closure)` |

⟹ a *mesh* is carrying the *method's* two operators. That is the same
"bound at construction" ruling the campaign already made for `Kernel × Frame →
Operator`, unapplied to the streaming path.

#### `[M]` both closures ALREADY ship their own adjoint — so the adjoint derives

| factor | its adjoint, today | status |
|---|---|---|
| `T_spatial` | `DiamondDifference.streaming_cell_transpose` (`:529`), `residual_kernel_batch_transpose` (`:485`), `reflect_scan_coefficients_transpose` (`:767`) — with `has_transpose_kernel` **DERIVED from the registration**, never declared (`scheme.py:1109`, #310 ruling 2) | ✅ ships |
| `T_ang` | `angular_adjoint` — a polymorphic family (`closure.py:557` base, `:1915` M-M, `:2128` Identity) | ✅ ships |
| `D` | real diagonal ⟹ self-adjoint | ✅ free |

⟹ `L.H = T_spatial.H @ T_ang.H @ D` is a **composition of adjoints that already
exist**. §5b's payoff — *"the adjoint stops being written and starts being
derived"* — is one binding away, not one implementation away.

⛔ **P0 (2026-08-29): the composition rule two lines up is REFUTED as
written** — the product it transposes does not equal `L`. What survives,
strengthened: both adjoints DO ship, and the measured derivation is
per-SUMMAND traversal reversal (`Lᵀ = D − Sᵀ·E_sp − A_angᵀ·E_ang`), which
the shipped `loss_action_transpose` already implements. See the ✅✅ ANSWERED
block below.

#### What the unweld concretely deletes

> ✅ Items 1–3 REMEDIED by P4.9a (`8a78be1d` + `737f8b32`); item 4 is
> P4.9b's. (2026-08-28; the owner's formula body now lives at
> `march_psi_half_step`, `closure.py:1370`.)

1. `DiamondDifference.update`'s `# ── Angular closure (Morel-Montry) ──` block —
   DD returns spatial-only, and **the twin of `closure.py:1327-1330` is gone**.
2. `cell_balance_terms` retires onto `cell_balance_for_streaming` (already
   imported at `diamond.py:87`). `[M]` one production consumer ⟹ #407's
   *"the whole module is a DD body in a scheme-neutral name and home"* closes
   with it, and its hard-coded `2.0 = 1/w_DD` goes with the blend weight the
   spatial closure supplies.
3. The L2 protocol sheds its angular members — `CellVisit.{tau, c_in, c_out}`,
   `UpstreamState.angular_upstream`, `CellResult.outgoing_angular_state` —
   leaving `transport/spatial/scheme.py` **purely spatial, which is what it
   claims to be**.
4. `SNMesh` sheds `scheme` and `pole_angular_closure`.

⭐ **And the coupling resolves itself.** `[M]` `closure.cell_contribution(...)`
already returns BOTH halves of `ΔA/w` — `angular_denom_term = (ΔA/w)·c_out` into
`D`, and `angular_numer_upstream = (ΔA/w)·c_in·ψᵃ` into `T_ang`. So the angular
closure already owns the coupling's USE, and the operator supplies `ΔA` as data.
That is fork 3's *"`delta_A` dissolves into R"* with a concrete consumer: R's
producer hands `ΔA` to the operator, the operator hands it to the closure.

⚠ **`abs_mu` stays with the SPATIAL closure, and that is correct** — a spatial
closure for a *directional* method is legitimately parameterized by direction
(`a = 2|μ|A_total/denom − 1`). Direction-parameterized ≠ angular-closure-aware.
✅ SHARPENED 2026-08-28 (user): *a parameter of the closure's EVALUATION,
supplied by the OPERATOR, stored by neither* — see P4.9's direction-supplier
ruling; μ reaches the scheme only inside assembled face coefficients.

⚠ **P0 still binds.** This mechanism is about WHERE the factors are composed and
WHO owns each adjoint. It does not establish the product form, its order, or the
arity — those remain the ⛔ P0 measurement above, and the mechanism is what makes
that measurement cheap to run (each factor becomes separately applicable).

### The FACTOR INVENTORY — what goes into each term, measured 2026-08-28

⚠ Read with the ⛔ P0 above: the **product form and its order are `[R]`**. What
follows is `[M]` — the assembled equation the tree actually computes, and which
coefficient of it belongs to which candidate factor. P0 is what turns the second
column from an assignment into a factorization.

#### The assembled equation `[M]` (verbatim from the two producers)

`cell_balance_for_streaming` (`transport/spatial/cell_balance.py:121`) and
`DiamondDifference.affine_scan_coefficients` (`diamond.py:572`) agree on:

```
denom[n,g,i] = 2|μ_n|·A_down[n,i]      ← streaming face        (spatial)
             + (ΔA[n,i]/w_n)·c_out[n]  ← curvature redistrib.  (COUPLING)
             + Σ_t[g,i]·V[n,i]         ← collision             (local)

numer_up[g,n] = |μ_n|·A_total[n]·ψˢ_in      ← upstream CELL     (spatial)
              + (ΔA/w)·c_in[n]·ψᵃ_in        ← previous ORDINATE (angular)

ψ̄ = (source + numer_up) / denom
a[n,g,i] = 2|μ_n|·A_total[i] / denom[n,g,i] − 1
```

closed by **spatial** `ψˢ_out = 2ψ̄ − ψˢ_in` (DD, `w = ½`) and **angular**
`ψᵃ_out = (ψ̄ − (1−τ)ψᵃ_in)/τ` (Morel–Montry), with `[M]` (`closure.py:1592`)

```
c_out = α_out/τ           c_in = (1−τ)/τ·α_out + α_in
```

⟹ `c_in`, `c_out` are pure `(α, τ)`: **F1's `A_angular(τ, α, w)` exactly.**

#### Which coefficient belongs to which factor

⛔ **CORRECTED by P0 (2026-08-29): the "its coefficients" column names
DIAGONAL scalings, not off-diagonals.** `[M]` `|μ|·A_total` and `(ΔA/w)·c_in`
are the diagonal coefficient matrices `E_sp` / `E_ang` that MULTIPLY the
accumulators; the accumulators' own off-diagonals are the DD alternating ±2
chain and the M-M `τ⁻¹` / `(1−τ)τ⁻¹` recurrence, and each accumulator is a
RESOLVENT (a dense lower triangle, sub-diagonal reach = chain length − 1),
not a nearest-neighbour band. Reading this column as nearest-neighbour
triangular factors is what made the product form look plausible. The table
STAYS per §3 — its coefficient-to-mechanism assignment is exactly what P0
consumed; only that reading is corrected.

| factor | what it is | triangular in | its coefficients | realized today by | its adjoint |
|---|---|---|---|---|---|
| **`D`** | the cell-local denominator — **not** "the collision term": all three mechanisms put their DIAGONAL part here | nothing (diagonal in `(n,g,i)`) | `2\|μ\|A_down` + `(ΔA/w)c_out` + `Σ_t V` | `cell_balance_for_streaming`; cached as `inverse_denom` | `D.H = D` — real diagonal, self-adjoint. This is what makes the reversal cheap |
| **`T_spatial`** | the chain scan: cell `i` ← its upstream neighbour, closed by DD | **chain order** (per ordinate) | off-diagonal `\|μ\|·A_total`; recurrence `ψ_out = a·ψ_in + b`, `a = 2\|μ\|A_total/denom − 1` | `affine_scan_coefficients` → `StreamingCoefficientCache` → `CumprodScan`/`ScanMarch` (`cumprod_a` is the *schedule* transform, deliberately not in the coefficient method) | the **reversed scan** — same coefficients, reversed order. ⭐ `[M]` its index already ships: `StreamingCoefficientCache.chain_idx_inv`, *"inverse permutation"* |
| **`T_ang`** | the ordinate march: `m` ← `m−1` within a μ-level, closed by M-M | **ordinate index within a level** | off-diagonal `(ΔA/w)·c_in`; recurrence `ψᵃ_{m+½} = (ψ̄_m − (1−τ_m)ψᵃ_{m−½})/τ_m` | `closure.compute_psi_half_per_level` / `precompute_psi_state` (batch), `cell_contribution` (per-cell); `tau_inv`, `mm_a_in_coeff`, `c_in` in the cache. ⛔ **twinned inline in `diamond.py`** (✅ twin REMEDIED by P4.9a `8a78be1d`) — see the layer constraint below | the **reversed march** — today hand-written as `angular_adjoint`, which is what §5b retires |

#### ⭐ `ΔA/w` is NOT a fourth factor — it is the scalar `D` and `T_ang` SHARE

`[M]` it appears exactly twice: `·c_out` inside **`D`'s diagonal**, and
`·c_in·ψᵃ_in` as **`T_ang`'s off-diagonal**. Its two halves are
**F2's `R₀₀ = ΔA`** (the spatial rank-1, by the divergence theorem) and
**`1/w`** (the angular weight).

⟹ `ΔA/w` *is* the spatial ⊗ angular coupling — which is why §4ter calls
`streaming_terms()` **"the fusion site"**, why fork 3 rules *"`delta_A`
dissolves into R"*, and why F3 names `_dAw_per_level` as **"the weld nobody
named"**. It is one quantity wearing two factors' clothes, and that is the
whole reason it has been homeless.

#### ⭐⭐ F3's two strata are visible in the CACHE'S OWN SHAPES

`[M]` `StreamingCoefficientCache` (`sn/sweep/cache.py:233-244`):

| shape | fields | stratum | factor |
|---|---|---|---|
| `(N,)` | `abs_mu`, `c_in`, `c_out`, `tau_inv`, `mm_a_in_coeff` | **mesh-free** — quadrature × chart only | `T_ang`'s algebra (+ `abs_mu`) |
| `(N, nx)` | `A_down`, `A_total`, `dA_w`, `V`, `chain_idx`, `chain_idx_inv` | **mesh-bound** | `D` and `T_spatial` |

That is F3 stated in the type system rather than in prose: the `(N,)` fields
cannot change under a re-mesh and the `(N, nx)` ones must. ⭐ And it is the same
split **P4b** was scheduled for from the L15 direction (7 fields bit-identical
across 3 meshes, 4 mesh-bound, 2 traversal-only) — two independent routes to one
carve, which is corroboration worth having.

#### ⛔ What is therefore NOT yet established, restated against this inventory

1. **The FORM.** The equation above is a **sum** — a diagonal plus two upstream
   couplings. That it equals a **product** `D @ T_ang @ T_spatial` is the P0
   claim, and the assignment table does not prove it.
2. **The ORDER.** Nothing measured here fixes `D @ T_ang @ T_spatial` over any
   other ordering; the sweep applies angular-outer / spatial-inner, which
   *suggests* one, and suggestion is not the check.
3. **The ARITY** may be 4 — the audit says *"two order reversals AND one block
   swap"*, and the pole-mirror permutation and the inflow/outflow trace
   partition are unplaced in the table above.
4. **Triangularity is order-relative**, so each factor is checked in its OWN
   index order — `T_ang` in the ordinate index at fixed cell, `T_spatial` in
   chain order.
5. ⛔ **Verify on the SPHERE.** `[M]` `R/ΔA` is bit-exactly `diag(1, 1/3)` on the
   cylinder and cannot discriminate.

### ⛔ The LAYER constraint on `T_ang` — measured 2026-08-28, and it is not in the section above

> ✅ **REMEDIED by P4.9a `8a78be1d` (2026-08-28), the same day it was
> measured.** The twin this section diagnoses is DEAD: DD returns
> spatial-only, the walk applies the owner's march
> (`closure.march_psi_half_step` via `advance_psi_half`), and the M8
> anti-twin gate makes the relation unspellable in `orpheus/transport/`.
> Of the two designs offered at the end, the second (remove the angular
> obligation from the spatial protocol) was ruled and landed. The section
> STAYS as the design record; its line numbers are the pre-carve tree's.

`T_ang` cannot simply be declared and consumed, because **its only per-cell
consumer may not name it.** `[M]` `DiamondDifference.update`
(`orpheus/transport/spatial/diamond.py`) evaluates the Morel--Montry angular
closure INLINE, under a comment that says so:

```python
# ── Angular closure (Morel-Montry) ──
psi_angle_out = (psi_avg - (1.0 - tau) * upstream_state.angular_upstream) / tau
```

and `[M]` that is a **Pattern-2 twin** of the owner,
`orpheus/sn/angular/closure.py:1327-1330`:

```python
psi_half[:, m + 1, :] = (
    psi_level[:, m, :] - (1.0 - tau_m) * psi_half[:, m, :]
) / tau_m
```

⭐ **Both spellings are LIVE, inside ONE class, selected by a data branch.**
`[M]` `_OneDimScanWalk._run` (`loss_representation/__init__.py`) routes normal
ordinates through `closure.precompute_psi_state(...)` (the owner) and
**degenerate cylindrical-axis ordinates** — `if geom.is_degenerate[global_n]:`,
its own comment *"slow per-cell path"* — through `scheme.update(...)`, i.e.
diamond's copy (`:4291`, consumed at `:4298`). That is the ERR-026 shape: one
relation, two implementations, a fix to one silently missing the other.

`[M]` **nothing gates them against each other** — 8 test files name
`precompute_psi_state`, 3 name `outgoing_angular_state`, and the intersection
is **empty** (both counts are positive controls for the negative).

⭐⭐ **And the duplication is FORCED, not sloppy.** `diamond.py` is **L2**
(`transport`), `closure.py` is **L3** (`sn`), and `transport → sn` is a declared
`FORBIDDEN_EDGES` violation. DD *cannot* call the owner, so it re-spells the
relation. ⟹ **declaring `T_ang` does not by itself remove the twin** — its L2
consumer still may not reference it. The design must therefore choose:

* **inject the closure** — L2 is HANDED the relation (a callable/strategy on the
  visit), staying blind to its provenance. This is the same *"hand it the
  closure, don't go shopping for one"* principle the L0 ladder took for
  `tau`/`edges`/`alpha` (P4.2), one layer up; **or**
* **remove the angular obligation from the spatial protocol** — the scheme
  closes the SPATIAL axis and returns the cell average; the SN walk applies the
  angular closure. `UpstreamState.angular_upstream` and
  `CellResult.outgoing_angular_state` leave `transport/spatial/scheme.py`. This
  makes the illegal state unrepresentable rather than merely well-injected, and
  it is what `𝓡 = R_spatial ⊗ A_angular` says the two axes are.

⚠ **This is the ANGULAR sibling of `#407`**, which records the same defect on
the SPATIAL side (`2.0 = 1/w_DD` hard-coded in scheme-neutral `cell_balance.py`).
#407's own diagnosis applies verbatim here: *"the containment is being done by a
capability fence, not by the algebra"* — `[M]` `LinearDiscontinuous` REFUSES
curvilinear (`linear_discontinuous.py:148-153`), so DD is the sole carrier of an
angular obligation the protocol declares for every member. Carve both with O-3.

⛔ *"Carve both with O-3" is SUPERSEDED 2026-08-28 — **P4.9** was chartered
(user-ruled, the same day) to own exactly this carve: it closes #407 and kills
the twin above. O-3 keeps the scheme-family mint design (charter §5c), which
P4.9 does not touch. §8's #407 row carries the matching amendment.*

⭐ Fair to the file, so nobody re-flags it: DD's consumption of `tau`, `c_in`,
`c_out`, `abs_mu` and `dA_w` as DATA is **correct** and the file insists on it
(*"DD must NOT rebuild them from `st.alpha_*` / `st.tau_mm`"*). The leak is
narrower and sharper than "diamond.py knows about angle": DD is closure-blind
about the CONSTANTS and closure-aware about the RELATION. And
`_reflection_coeffs`' `w` is DD's **blend weight** (½), not a quadrature
weight — a false alarm.

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
(`transport/spatial/scheme.py:1677-1718`), and LD's override
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
(`augmented_mesh.py:1228`, the trial space), **not** from `__init__`. So the
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
Production: `transport/fields/_bases.py:266`, `transport/spatial/scheme.py:1716`
(internal, inside `moment_axis`), `sn/mesh/augmented_mesh.py:1228`.
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
⟹ P2 leaves every test file where it is. Only `tests/sn/sweep/curvilinear/test_angular_closure.py`
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

### P4.9 — each axis is closed by its own closure, and the operator composes them *(behavioural; closes #407)*

**Goal.** A spatial discretization scheme closes the **spatial** axis and
nothing else. The angular closure closes the **angular** axis. Neither knows the
other exists, and the object that owns both — the operator — is the only thing
that composes them. ⟹ the Morel–Montry relation has exactly ONE spelling in the
tree, and `transport/spatial/` is method-agnostic in fact and not only in name.

**Proposed means** (✅ user-ruled 2026-08-28; the *shape* is ruled, the edit
list below is proposed and unverified in detail):
`StreamingOperator` is constructed from **(domain, codomain, spatial closure,
angular closure)**. Discretization schemes become **factories returning a
spatial closure**; the angular scheme returns an **angular closure**; those four
arguments fully determine a well-posed operator whose adjoint is derived.

⭐⭐ **SHARPENED 2026-08-28 (user, ratified after a feasibility check): the
scheme is direction-agnostic in SIGNATURE, not only in state — the OPERATOR
supplies μ.** Ruled on the Frame analogy (the stage-2 generator hands its
minted product to the operator; the operator owns the evaluation point).
`[M]` three facts already in its favor: DD carries **no fields** and LD only
`theta` — the schemes are direction-STATELESS today, on every path; the
vectorized path already injects direction from the SN side
(`affine_scan_coefficients(abs_mu=geom.abs_mu, …)`, called by the cache
builder); and in the assembled equation μ enters ONLY through the products
`|μ|·A_down` / `|μ|·A_total` (LD's `g = |μ|/h` is `(|μ|A_down)/V` on the
slab) — so the agnostic closure consumes **assembled face coefficients**,
and a diffusion box scheme fills the same slots with face conductances
`D·A/Δx`, no direction anywhere. ⟹ the scheme family becomes agnostic of
**directionality itself**, not merely of which method — this RULES the
direction of O-3's decomposition-by-destination for the μ-bearing slots.
⚠ The boundary, so nobody over-reads it: the closure FUNCTION is
direction-free; the evaluated TABLE cannot be (τ depends on |μ|) — the
evaluation relocates to the operator's side (§5c's two-layer split), it
does not disappear. Sweep-direction resolution (sign μ →
`face_area_downstream`) and degenerate-ordinate signalling (`|μ|A ≈ 0`
data) stay operator-side by construction. ⚠ Until CS5's NODAL half lands,
the operator reads μ from the quadrature; after it, from its own domain
space's angular axis — which closes the campaign's root ruling.

⭐⭐ **This is a COMPLETION, not a new design — `[M]` the same unweld already
shipped on the vectorized path at Phase 2.11 (#197 PR-TYPED-6.5).**
`cell_balance_for_streaming` is already **closure-blind by interface**: it takes
`angular_denom_term` / `angular_numer_upstream`, injected from
`closure.cell_contribution(...)`. Its twin `cell_balance_terms` still names
`c_in`/`c_out` and has `[M]` **one** production consumer (`diamond.py:194`). The
per-cell arm is the sole holdout, and it is the arm carrying the duplicate.

**Why the layer permits it, and only here.** `sn/operators/streaming.py` is
**L3**: `L3 → L2` (`transport/spatial`) is legal, `L3 → L3` (`sn/angular`) is
the same package, and **L2 never names L3**. The twin dies because the
*obligation* dies — strictly better than injecting the relation into L2.

**The edit list** (proposed; each row needs its own §6b enumeration at design
time):

1. `DiamondDifference.update` — delete the `# ── Angular closure (Morel-Montry) ──`
   block; return spatial-only. Route through `cell_balance_for_streaming`
   (already imported at `diamond.py:87`) instead of `cell_balance_terms`.
2. `cell_balance_terms` **retires** onto the blind helper. `[M]` one production
   consumer ⟹ **#407 closes with it**, and its hard-coded `2.0 = 1/w_DD` goes
   with the blend weight the spatial closure now supplies.
3. `transport/spatial/scheme.py`'s protocol sheds its angular members —
   `CellVisit.{tau, c_in, c_out}`, `UpstreamState.angular_upstream`,
   `CellResult.outgoing_angular_state`.
4. `StreamingOperator` gains `(domain, codomain, spatial_closure,
   angular_closure)`. `[M]` today it is a `@dataclass` with **exactly one
   field** (`sn_mesh`) and derives `domain`/`codomain` as properties.
5. `SNMesh` sheds `scheme` (`augmented_mesh.py:268`) and `pole_angular_closure`
   (`:422`) — a *mesh* is carrying the *method's* two operators.

✅ **RULED 2026-08-28 (user, design round) — the phase SPLITS, and two forks
with it:**

* **P4.9a** = rows 1–3 **plus the `mm_a_in_coeff` handing** (next bullet) plus
  the `is`-identity gate. Closes **#407**, kills the twin. **P4.9b** = rows 4–5
  (the ~150-ctor-site + ~60-read migration; sizing in
  `scratch/p4_9_design_measured.md`). §6b at the boundary, verified: the two
  sub-phases change **disjoint interfaces** — rows 1–3 change the scheme
  protocol's members (consumers: diamond, cell_balance, linear_discontinuous,
  the walk, 12 test files — all inside P4.9a's own enumeration); rows 4–5
  change the operator ctor + mesh fields (~150 + ~60 sites). No signature's
  call-site set is cut by the boundary. Interval state: the mesh still BUILDS
  both closures (`augmented_mesh.py:268`/`:422` untouched by P4.9a); the scheme
  just no longer APPLIES the angular one — the walk does, reaching the same
  object the vectorized branch reads, which is exactly what turns the
  `is`-identity gate green at P4.9a. The done-when below splits: bullets
  1/2/4/5 → P4.9a; bullet 3 (ctor unconstructable without closures) → P4.9b.
* **`mm_a_in_coeff` folds into P4.9a.** `[M]` `sn/sweep/cache.py:377` derives
  `(1.0 - tau) / tau` itself — the cache spelling a fragment of the M-M
  relation (the memo's "same smell one notch down"). The closure exposes its
  march constants; the cache is handed them. That spelling sits outside the
  done-when's `transport/` grep scope, so the handing gets its own pin (the
  cache's `(N,)` closure-algebra fields bit-identical before/after).
* **P4.9b's migration lever is a PRODUCTION posing-head factory** (the
  Frame/stage-2 pattern): the posing head mints both closures and hands them
  to `StreamingOperator`; tests consume the same production site — a
  tests-only builder would twin the assembly (Pattern 2). Binds P4.9b's
  design round; details ruled there.

✅ **RULED 2026-08-28 (user, execution round — after the two pre-measurements
`scratch/p4_9a_enumeration.md` + `scratch/p4_9a_verification_plan.md`):**

* **Q1 — the degenerate branch routes through the OWNER's form, via a closure
  march-step method** (Form A, the closure's own operation order; its batch
  march delegates to the same body). `[M]` the expanded scan form differs from the
  owner's by max |Δ| = 1.776e-15 on real τ (the seed-stable discriminator;
  %/ULP figures are draw-dependent — memo §F2 band correction) and breaks
  `keff` bit-equality on 3 of 4 configs; and the `is`-identity gate is unbuildable without a closure METHOD
  call — two independent reasons, one ruling. The fast path keeps the minted
  scan-form constants (bit-identity of every cylindrical solve), welded to the
  march by an `array_equal` field gate (`[M]` the realistic drift,
  `tau_inv − 1.0`, is 1–2 ULP — any tolerance is a non-catcher).
* **Q2 — CellVisit sheds all three** (`tau`, `c_in`, `c_out`); `update` /
  `residual` take the two ASSEMBLED angular contributions
  (`angular_denom_term`, `angular_numer_upstream`) — `residual`'s existing
  internal idiom and `cell_contribution`'s return type. The walk assembles
  from the visit's own `st.delta_A_over_w` (NOT the cache's separately-produced
  `dA_w` twin) to preserve bit-identity; readership of the packet field moves
  from `cell_balance_terms` to the walk — **P4.7's enumeration must include
  it**. The mesh's closure-stamp in `_make_cell_visit` dies with the fields.
* **Row 3b rides P4.9a**: `affine_scan_coefficients(dA_w=, c_out=)` re-poses
  onto an assembled `angular_denom_term=` (caller assembles; DD's in-L2
  `dA_w·c_out` at `diamond.py:676` dies). LD's scan guard re-keys on the
  assembled term; the per-cell guard re-keys on
  **`face_area_inner != face_area_outer`** — the geometric curvature truth,
  mesh-free and **P4.7-proof** (`delta_A_over_w`, the architect's candidate
  signal, is a field P4.7 retires; the assembled term is non-universal at dome
  ends where `c_out = 0`).
* ⭐ **DIRECTION, binding P4.9b (user, same ruling):** the fast path should
  stop having *"its own totally special implementation"* — the fused scan
  table (the rearranged+welded form: scan-normal coefficients + the
  `ΔA/w·c_in` spatial⊗angular fusion) is the **OPERATOR's artifact**, minted
  near the hot path from the two closures the operator holds. P4.9a lands the
  first half (the closure mints the constants; the spelling is welded by
  gate); P4.9b's operator assembles the table and the walk consumes operator
  artifacts plus one named kernel.

✅ **P4.9a LANDED 2026-08-28 — 12 commits `cb65c4cc`…`7a0f434c` (docs
`ca852c44`), merged with the full fast gate green (`9836/0`, 13 trees
`rc=0`, 64 min, delta exactly the +12 new gates).** The done-when audit is
under each bullet below; the battery table is the tracked memo §9; #407
was closed (by the user, mid-phase) and carries the completion comment.
P4.9b inherits: rows 4–5, the posing-head-factory lever, the
operator-minted scan-table direction, and the is-identity gate's control
leg (the fast path stays call-free) as a standing constraint.

✅✅ **P4.9b DESIGN ROUND COMPLETE 2026-08-28 (user, two ask rounds) — the
full record is `scratch/p4_9b_design.md` §§7–8; the rulings, compact:**

* ⛔ **Row 5 is REVISED BY RULING — the mesh sheds NOTHING.** `SNMesh` is
  recognized as a misnomer: it is the **save-state / data hub**, and it
  KEEPS `scheme` (stage-2 generator: space induction nodal/modal, and
  cross-consumer consistency — DSA must read the SAME generator) and the
  bound `pole_angular_closure` (shared machinery; one instance ever).
  The 67 consumer reads dissolve by consumers flowing through the
  OPERATOR, not by the mesh losing fields — no space-side re-pose, no
  `__eq__` change, no `_bases.py:257` re-key. (The rename-the-misnomer
  question is filed as its own issue, out of scope.)
* **The ctor**: `StreamingOperator(sn_mesh, scheme, pole_angular_closure)`
  — three required fields, **no defaults** ("the scheme is an active
  choice"), **no guards** (the user's no-guard position was attacked four
  ways and HELD — pose-path unspellable; the raw ctor is the declared
  expert seam. ⛔ CORRECTED per verification plan F12 `[M]`: the
  wrong-FAMILY arm RAISES at first sweep (typed sphere / untyped
  cylinder), it is NOT silent — only the cross-mesh-smuggling arm is
  silent. Both documented in the ctor docstring, not guarded: a guard
  would forbid the seam's own use case, doctored diagnostic probes). `spatial_closure` property = the identity
  extraction seam until O-3 splits the closure/factory family; the
  `scheme` field IS §5c item (d)'s provenance accessor — **item (d)
  discharged**.
* **The posing surface**: classmethod **`StreamingOperator.pose(sn_mesh)`**
  — the INTERMEDIATE while the migration runs — reads the hub's two
  objects and passes them. `build_streaming_collision` routes through it;
  ⭐ solve entries UNCHANGED (`scheme=` keeps entering the hub at
  `_as_sn_mesh`; the hub ctor stays the active-choice site and its DD
  default SURVIVES). `[M]` 136 ctor sites migrate (1 production + 135 in
  40 test files — recount `scratch/p4_9b_row45_recount.md`).
* ⭐⭐ **The algebra/performance principle (user, verbatim in the memo):
  the algebra stays unwelded and expressible as long as possible;
  performance welds resolve LAZILY, as close to the solution STRATEGY as
  possible.** ⟹ the operator owns/exposes the ALGEBRA (two closures +
  their minted per-ordinate constants); the fused scan table is the
  STRATEGY's artifact, lazily resolved from the operator's objects — the
  eager `SNSolver.__init__:1434` build retires. ⛔ CORRECTED per
  verification plan F2 `[M]`: the memo's LIFETIME is a ruled question
  (Q1) — the operator is built 6–10×/solve, so a per-operator memo costs
  up to 24.65 % of a slab solve; the perf gate is a COUNT (today 1).
  (Memory: `feedback_algebra_eager_performance_lazy.md`.)
* **The end state (recorded, NOT this phase)**: cross-method
  `(domain, codomain, spatial-discretization[, angular-discretization])`
  ctor with `sn_mesh` out — rides O-3/CS5. The mesh field + derived
  spaces are the declared transitional weld.
* ✅ **THIRD ask round (2026-08-28) — Q1 + Q5 ruled by criterion, picks
  delegated** (full synthesis: `scratch/p4_9b_design.md` §9): the table's
  interning mechanism lives IN the `loss_representation` layer (WeakKey on
  the hub, closure-pair-validated; count-1 gated; the layer is
  retirement-bound at Campaign 2, so the interim machinery retires with
  it — the user's survives-the-lazy-strategy criterion). The slot
  vocabulary is the SYMMETRIC pair **`spatial_closure` +
  `angular_closure`** on operator AND representation (today
  `spatial_closure` receives the hub's scheme instance; extraction =
  identity until O-3); the hub keeps `scheme` (generator, provenance);
  ⭐ **the pole misnomer dies at step 3**: `pole_angular_closure` →
  `angular_closure` (91 hits) + `PoleAngularClosureBase` →
  `AngularClosureBase` (119 hits), member names untouched. Adopted from
  the verification plan: representation `.pose` classmethod (58 sites),
  route gate self-drives, selection predicates operator-side, F5's
  21 space-side reads STAY hub-side, both permanent gates, F10's
  diagnostic-probe carve-out. ⛔ F12/F2 corrections above stand.
* **Steps**: (1) ctor + `.pose` + the 136-site migration, bit-identical
  by construction (same objects flow); (2) the operator feeds the walk —
  `default_for` + representation take the closure pair, 43+17+5 reads
  re-plumb, the 2 L2 reach-throughs are HANDED their objects, the
  strategy's lazy table, `cache.py:273` assert→raise ride-along; (3)
  docs (~127 prose mentions) + owner-mutation battery + `dead_references`
  + banners + landing record. Test-architect verification plan BEFORE
  step 1 (the proactive carve trigger).
  ✅ **EXECUTED 2026-08-28/29 with three recorded refinements** (the full
  trail is the commit sequence `b253732f`…`5a1591d2` + the memos): step 1
  staged as witness-red → pose-FORWARDER → mechanical migration →
  atomic flip, so every commit stayed green and the §6b unit collapsed
  to a point; step 2's §6b census gained three member SPELLINGS a
  name-keyed AST census cannot see (calls through a VARIABLE
  `rep_cls(sn)`, the base class's own internal `supports` call, a
  monkeypatch surrogate lambda — all caught by red loop; surprise-log
  row below); the L2 `_bases.py:257` read SURVIVES by the hub ruling
  (F9 voided the recount's red predictions), and
  `radial_characteristic_field.source_from_angular` stays as the one
  posing-time residual (#414). One deviation from the architect's
  numbers, recorded in the memo §10: the F2 operator count is
  fixture-dependent (6–10 vs 38–43/solve) — the COUNT gate, not a
  percentage, is the instrument.

**Done when** (checkable):

* `[M]` `grep -c "1.0 - tau"` over `orpheus/transport/` returns **0** — the
  Morel–Montry relation is spelled once, in `sn/angular/closure.py`.
  ✅ **MET at `8a78be1d`** and made PERMANENT by the M8 anti-twin gate
  (`tests/transport/spatial/test_no_angular_closure_twin.py` — text-surface
  spellings + AST-surface owner names, per-pattern positive controls).
  ⚠ **SCOPED 2026-08-28** (`scratch/p4_9a_verification_plan.md` §F2): this
  makes `transport/` clean and single-homes the relation's DERIVATION in the
  closure — which after P4.9a owns BOTH its representations (the Form-A march
  + the minted scan-form constants, welded by gate). The walk's fast path
  consumes the minted constants (restructuring it onto operator-minted
  artifacts is the P4.9b direction above); the transpose's hand-written
  adjoints are **§5b's named target**, untouched here. "Exactly ONE spelling
  in the tree" is the arc's end-state, not this phase's.
* `transport/spatial/scheme.py` contains no `tau` / `c_in` / `c_out` /
  `angular_*` member on any protocol type. ✅ **MET at `737f8b32`** (row 3 +
  row 3b: the scan API's `c_out=`/`dA_w=` kwargs also re-posed onto the
  assembled term), pinned by the structural gate's `dataclasses.fields`
  legs.
* `StreamingOperator(...)` cannot be constructed without both closures — the
  illegal state is unrepresentable, not merely undesirable.
  ✅ **MET at P4.9b** (three required fields, no defaults; the arity witness
  pins it with `match="missing 2 required positional"`; every un-migrated
  spelling is a loud collection-time `TypeError`). The companion safety
  argument is the pose-identity gate (no ctor guards, by ruling — the
  four-attack record is `scratch/p4_9b_design.md` §8, its M5 enforcement
  measurement `scratch/p4_9b_verification_plan.md` §9).
* `[M]` the `_OneDimScanWalk` degenerate-ordinate branch and the vectorized
  branch reach the **same** closure object — a `is`-identity assertion, which is
  the gate the twin never had (`[M]` today the intersection of tests naming
  `precompute_psi_state` and `outgoing_angular_state` is **empty**).
  ✅ **MET** — landed RED at `066982cf` (strict xfail, verified failing on
  leg 1), flipped at `8a78be1d`; the matvec's degenerate arm is the second
  per-cell leg and fp(4,8) is the zero-calls control (the fast path stays
  on minted constants — Q1's perf half, now P4.9b's standing constraint).
* ~~the aniso curvilinear canary is bit-identical~~ ⛔ **REFUTED 2026-08-28
  `[M]`** (`scratch/p4_9a_verification_plan.md` §F1): the canary CANNOT
  witness this carve — `DiamondDifference.update` executes only on cylinders
  with `n_phi ≡ 2 (mod 4)` (counting spy: **0** calls on slab, sphere, and
  every `n_phi ∈ {4, 8}` config in full eigenvalue solves; the canary uses
  8/16), so its bit-identity is unconditional. **Replaced by the BUILT
  degenerate set**: the `CYL_DEG` (`folded_product(4,6)`) affine-carve
  baseline row — landed PRE-carve on unmodified production — plus a fp(4,6)
  eigenvalue `array_equal` pin, the degenerate sibling of the frozen matvec
  baseline, and the `[M]` 13 existing twin-catcher rows staying green. Slab +
  sphere stay as CONTROLS (they carry zero information about the carve; the
  plan says so, so nobody cites them as evidence). Plus `sphinx -W` 0;
  pyright 0.

⛔ **§6c — this phase's gate must land with the case it catches.** The
`is`-identity gate above is the witness: before the unweld it CANNOT pass (there
are two objects), which is exactly the red that proves it has teeth. Record that
reading in the gate's docstring.

**Ordering.**
* **After P4.3** — both touch `transport/spatial/scheme.py`; P4.3 is the smaller
  diff and its module deletion should not be entangled with a protocol change.
* **Before §5b's ⛔ P0** — P0 must apply each factor separately to measure the
  product form, and this phase is what makes them separately applicable. Running
  P0 first would mean measuring a factorization whose factors are still welded.
* `[R]` **it probably makes P4.7 a consequence rather than a task**: once the
  per-cell arm takes the coupling from `closure.cell_contribution(...)`,
  `StreamingTerms.delta_A_over_w`'s readers reduce to the cache builder, which
  can read `ΔA` from the connection directly. **Not measured** — verify before
  reordering P4.7. (⛔ ADJUDICATED 2026-08-29, `scratch/p4_7_ground_measure.md`:
  REFUTED — two production readers, the second the degenerate cylinder arm
  kept deliberately by P4.9a Q2. P4.7 executed as a task, `3456dd37`.)

**What it buys beyond the deletion.** `[M]` both factors already ship their own
adjoint — `DiamondDifference.streaming_cell_transpose` (with
`has_transpose_kernel` **derived from the registration**, `scheme.py:1109`) and
the polymorphic `angular_adjoint` family (`closure.py:557/:1915/:2128`). So
after this phase `L.H = T_spatial.H @ T_ang.H @ D` composes adjoints that
already exist, and §5b's *"the adjoint stops being written and starts being
derived"* is one binding away. ⛔ P0 (2026-08-29) refuted the PRODUCT reading
of this sentence — see §5b's ✅✅ ANSWERED block: the derived adjoint is
per-summand traversal reversal, not `T_spatial.H @ T_ang.H @ D`; the
"adjoints already ship" half stands.

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
  ✅ **AMENDED 2026-08-28 — P4.9 now owns the carve and closes #407** (chartered
  the same day, user-ruled). The row predates the charter: it was written when no
  phase of this plan could hold the fix. What stays with O-3 is the scheme-family
  mint design (charter §5c) and `cell_balance.py`'s broader reorganisation — not
  #407. ⚠ §3's REMEDIED-fact case, and the second amended row in this section.
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

# Phase A — Reference blast-radius map for the `docs/theory/` restructure

**Anchor:** HEAD `49135ab6`, branch `docs/sn-doc-architecture`, 2026-07-15.
**Scope:** whole repo. **Excluded:** `.claude/worktrees/nexus-workspace-wiring/`
(a *separate git checkout*, git-excluded via `.git/info/exclude:11`; it has its
own `docs/` and its own state — it is NOT part of this restructure). `docs/_build/`
excluded. `scratch/` reported but untracked.

**Nothing was moved, edited, or created** except this file.

---

## 0. Executive framing — which reference classes actually break

Restructuring is a **path** change. Sorting every "reference" by whether it is
path-bearing and whether a build can catch it is the whole game:

| Class | Path-bearing? | `sphinx-build -W` catches a break? | Count | Verdict |
|---|---|---|---|---|
| `:ref:` role | **NO** — resolves by global label name | yes (undefined label) | 1076 in `docs/` | **IMMUNE to `git mv`** |
| `:eq:` / `:numref:` | **NO** — global label | yes | — | **IMMUNE** |
| `:doc:` role (rendered) | **YES** | yes | 10 live in `.py` + 188 in `.rst` | **BREAKS, loudly** |
| toctree entry | **YES** | yes (`toctree contains reference to nonexisting document`) | 62 entries / 14 directives | **BREAKS, loudly** |
| `.. include::` | **YES** (relative) | yes (include file not found) | 10 | **BREAKS, loudly** |
| `:doc:`/`:ref:` in **non-autodoc'd** `.py` | YES | **NO — inert text** | 71 `:doc:` + 54 `:ref:` | **BREAKS, SILENTLY** |
| `:file:` role | YES | **NO — renders as literal text** | 156 | **BREAKS, SILENTLY** |
| Raw path string in prose/comment | YES | **NO** | 837 | **BREAKS, SILENTLY** |
| **Python `Path` constants** (`ROOT / "docs" / "theory"`) | YES | **NO** — and a `docs/theory/` grep **cannot see them** | 4 (+1 test) | **BREAKS, SILENTLY + RED TESTS** |

The last row is a class the brief did not anticipate and is the **#1 hazard**
(§7). The `:ref:` immunity (row 1) is the single biggest piece of good news:
**1076 of the ~1150 total cross-references need no work at all.**

---

## 1. Every `:doc:` role usage in the repo

**Totals.** 269 genuine `:doc:` role usages in Sphinx-relevant sources
(`.rst` + `.py`), worktree excluded:

| Source bucket | Count | Sphinx renders it? |
|---|---:|---|
| `docs/**/*.rst` | 188 | **YES** |
| `orpheus/**/*.py` (docstrings) | 54 | only if the module is `automodule`'d — **10 live / 44 inert** |
| `tests/**/*.py` (docstrings) | 27 | **NO** — nothing under `tests/` is ever `automodule`'d → all inert |
| **TOTAL (real roles)** | **269** | |

Plus **39** occurrences of the *text* `:doc:` inside `.claude/**/*.md`
(36 of them in `agent-memory`) — these are prose *about* the role in planning /
memory documents, not Sphinx input. Not build-relevant; listed in §4 treatment.

**Absolute vs relative split** (this decides the fix strategy):

| Form | Count | Behaviour under the move |
|---|---:|---|
| **ABS** `:doc:`/theory/x`` | 125 | source-dir-relative → breaks iff **the target** moves |
| **REL** `:doc:`x`` | 144 | resolved against **the source page's dir** → breaks if **either** the source *or* the target moves |

**REL is the dangerous form**: 132 of the 144 REL roles live *inside*
`docs/theory/` and are bare sibling names (`:doc:`operator_algebra``). The moment
`discrete_ordinates.rst` and `operator_algebra.rst` land in *different*
subdirectories, every one of those 132 needs a prefix. Converting them all to
**ABS** (`/theory/methods/operator_algebra`) makes them depend only on the target
— strongly recommended.

### 1a. ABS `/theory/*` `:doc:` roles by source bucket (break when the target page moves)

| Source | Count |
|---|---:|
| `orpheus/**` | 46 (**38 of them inert/silent** — see §1d) |
| `docs/**` | 22 |
| `tests/**` | 9 (all inert) |
| **total** | **77** |

### 1b. Inbound `:doc:` count per target page — the move-map priority list

Pages with the most inbound `:doc:` (each is a fix site when that page moves):

| Target page | inbound REL | inbound ABS | total |
|---|---:|---:|---:|
| `operator_algebra` | 28 | 3 | **31** |
| `discrete_ordinates` | 27 | 21 | **48** |
| `collision_probability` | 21 | 2 | **23** |
| `peierls_nystrom` | 7 | 15 | **22** |
| `diffusion_1d` | 2 | 9 | **11** |
| `boundary_conditions` | 2 | 9 | **11** |
| `reference_solvers` | 1 | 7 | **8** |
| `frame` | 7 | 0 | **7** |
| `trajectory_resolvent` | 6 | 0 | **6** |
| `loss_representations` | 5 | 2 | **7** |
| `verification/reference_solutions` | 2 | 13 | **15** |
| `development` | 0 | 11 | **11** |
| `verification/index` | 0 | 5 | **5** |
| `verification/matrix` | 4 | 3 | **7** |
| `sood_registry` | 4 | 0 | 4 |
| `structured_geometry` | 3 | 2 | 5 |
| `singular_eigenfunction` | 3 | 2 | 5 |
| `discrete_measures` | 3 | 0 | 3 |
| `homogeneous` | 1 | 2 | 3 |
| `monte_carlo` | 2 | 0 | 2 |
| `fn_method` | 2 | 1 | 3 |
| `api/numerics` | 3 | 0 | 3 |

### 1c. Full `docs/**` `:doc:` list (188)

*(exact `file:line` → target; `ABS` = leading `/`)*

**Cross-directory (`docs/api`, `docs/architecture`, `docs/development`)**

| file:line | kind | target |
|---|---|---|
| `docs/api/data.rst:90` | REL | `derivations` |
| `docs/api/data.rst:91` | REL | `../theory/verification` |
| `docs/api/diffusion_1d.rst:9` | ABS | `/theory/diffusion_1d` |
| `docs/api/monte_carlo.rst:59` | REL | `homogeneous` |
| `docs/api/transport.rst:24` | ABS | `/theory/boundary_conditions` |
| `docs/architecture/layering.rst:109` | ABS | `/verification/index` |
| `docs/architecture/layering.rst:232` | ABS | `/verification/index` |
| `docs/development.rst:258` | REL | `testing/architecture` |
| `docs/development.rst:259` | REL | `verification/index` |

> `docs/api/data.rst:91` (`../theory/verification`) and `docs/api/monte_carlo.rst:59`
> (`homogeneous` → resolves to `api/homogeneous`) are **REL crossing a directory
> boundary** — the most fragile form in the corpus.

**`docs/theory/` → `docs/theory/` siblings (bare REL; the 132-strong bulk)**

`boundary_conditions.rst`: 126→`operator_algebra`; 2079, 2214, 2612, 3503→`discrete_ordinates`
`collision_probability.rst`: 1738→`discrete_ordinates`; 3306, 3354, 3409, 4051, 4326→`peierls_nystrom`
`cross_section_data.rst`: 207→`sood_registry`
`discrete_ordinates.rst`: 88, 247, 6333→`loss_representations`; 108, 153, 7315, 7461, 7480, 11071, 13246, 14085, 14779, 14785, 15003, 18354, 18406, 18498, 18624, 18882, 18933, 18985, 19025, 19063→`operator_algebra`; 2331, 2507, 12053→`structured_geometry`; 5588, 6130, 6258→`discrete_measures`; 7317, 9814→`boundary_conditions`; 15655→`method_of_characteristics`; 18121, 18138, 18145, 18172, 18177→`frame`; 18822→`homogeneous`; 95, 152→`../api/numerics`; 15491→`../verification/matrix`
`fn_method.rst`: 71, 2450→`trajectory_resolvent`; 110, 2453→`peierls_nystrom`; 179, 2451→`singular_eigenfunction`; 571→`reference_solvers`; 715→`sood_registry`
`frame.rst`: 2706, 2762, 2804→`operator_algebra`
`homogeneous.rst`: 1544→`verification`
`index_convention.rst`: 1314→`operator_algebra`; 1597→`../api/discrete_ordinates`
`loss_representations.rst`: 18, 154, 301, 2820, 2865→`discrete_ordinates`; 27, 754, 2869→`operator_algebra`
`method_of_characteristics.rst`: 806→`discrete_ordinates`; 823→`../api/numerics`
`monte_carlo.rst`: 294→`discrete_ordinates`
`operator_algebra.rst`: 2143, 4444, 5015, 5098, 5103, 7010, 7529, 7625, 8001, 8036, 8670, 8684, 10036→`discrete_ordinates`; 2247→`frame`; 9759→`collision_probability`; 9839, 10612→`loss_representations`; 9840, 9991→`diffusion_1d`
`peierls_nystrom.rst`: 1214, 5440→`trajectory_resolvent`; 1450, 1492, 1646, 2102, 2129, 2168, 5485, 5560, 5633, 5897, 5899, 5913, 6005, 6066, 6098, 6104, 6158, 6720, 8470, 8474→`collision_probability`; 8095, 8165→`monte_carlo`; 5642→`../verification/reference_solutions`
`singular_eigenfunction.rst`: 1952, 2076→`sood_registry`; 2072→`fn_method`; 2074→`trajectory_resolvent`
`sood_registry.rst`: 914→`fn_method`; 916→`singular_eigenfunction`; 918→`trajectory_resolvent`
`spherical_harmonics.rst`: 213→`frame`
`structured_geometry.rst`: 288, 342→`discrete_ordinates`

**`docs/theory/` ABS forms**

| file:line | target |
|---|---|
| `collision_probability.rst:718, 2104` | `/verification/reference_solutions` |
| `collision_probability.rst:2103, 4337` | `/theory/peierls_nystrom` |
| `collision_probability.rst:4338` | `/theory/discrete_ordinates` |
| `cross_section_data.rst:647` | `/verification/index` |
| `diffusion_1d.rst:179, 431` | `/theory/boundary_conditions` |
| `diffusion_1d.rst:189, 514` | `/theory/operator_algebra` |
| `diffusion_1d.rst:324, 369` | `/theory/homogeneous` |
| `diffusion_1d.rst:1193` | `/verification/reference_solutions` |
| `discrete_ordinates.rst:5387` | `/theory/reference_solvers` |
| `discrete_ordinates.rst:6487, 7095, 8312` | `/development` |
| `discrete_ordinates.rst:6488` | `/verification/matrix` |
| `discrete_ordinates.rst:6632` | `/theory/boundary_conditions` |
| `galerkin_spectral.rst:709` | `/theory/reference_solvers` |
| `operator_algebra.rst:2025, 4318, 4412, 6352, 11138` | `/theory/diffusion_1d` |
| `peierls_nystrom.rst:344` | `/api/derivations` |
| `reference_solvers.rst:202` | `/verification/matrix` |
| `reference_solvers.rst:344` | `/verification/index` |
| `reference_solvers.rst:346` | `/verification/reference_solutions` |
| `sood_registry.rst:694` | `/testing/architecture` |
| `trajectory_resolvent.rst:2351, 2633, 2848, 2885, 2978, 3251, 3274, 3286` | `/development` |
| `verification.rst:509` | `/theory/diffusion_1d` |

**`docs/verification/`**

| file:line | kind | target |
|---|---|---|
| `docs/verification/index.rst:25` | REL | `reference_solutions` |
| `docs/verification/index.rst:29, 53` | REL | `matrix` |
| `docs/verification/index.rst:32, 133` | ABS | `/theory/reference_solvers` |
| `docs/verification/index.rst:54` | REL | `reference_solutions` |
| `docs/verification/reference_solutions.rst:152, 410` | REL | `matrix` |

### 1d. `orpheus/**` `:doc:` roles (54) — **LIVE vs INERT**

Only **45 modules** are ever `automodule`'d anywhere in `docs/`. A `:doc:` in a
module outside that set is **never rendered** — so `-W` cannot catch it, and it
rots exactly like a raw path string.

**LIVE (10 `:doc:` — a target move breaks the `-W` build):**

| file:line | target |
|---|---|
| `orpheus/derivations/common/quadrature.py:45` | `/theory/peierls_nystrom` |
| `orpheus/derivations/common/quadrature.py:300` | `/theory/peierls_nystrom` |
| `orpheus/derivations/continuous/analytical/homogeneous.py:35` | `/verification/reference_solutions` |
| `orpheus/derivations/continuous/cases/diffusion.py:48` | `/theory/diffusion_1d` |
| `orpheus/derivations/continuous/cases/diffusion.py:50` | `/verification/reference_solutions` |
| `orpheus/derivations/continuous/flat_source_cp/cylinder.py:9` | `/theory/peierls_nystrom` |
| `orpheus/derivations/continuous/flat_source_cp/slab.py:6` | `/theory/peierls_nystrom` |
| `orpheus/sn/loss_representation/__init__.py:101` | `/theory/loss_representations` |
| `orpheus/sn/loss_representation/__init__.py:120` | `/theory/loss_representations` |
| `orpheus/transport/method.py:105` | `/theory/boundary_conditions` |

**INERT (44 — silent staleness; 38 of them ABS `/theory/*`):**

`orpheus/derivations/__init__.py:32`→`/verification/reference_solutions`;
`common/continuous_reference.py:30`→`/verification/reference_solutions`;
`common/solution_types.py:87`→`/skills/algebra-of-record` **(already broken)**;
`common/solution_types.py:92`→`/skills/vv-principles` **(already broken)**;
`continuous/cases/sn.py:23, 62, 733`→`/theory/discrete_ordinates`;
`continuous/cases/sn.py:64`→`/verification/reference_solutions`;
`continuous/flat_source_cp/geometry.py:4`, `sphere.py:8`→`/theory/peierls_nystrom`;
`continuous/fn_method/moment_space.py:62`→`/theory/fn_method`;
`continuous/galerkin_spectral/basis_space.py:140`→`/theory/galerkin_spectral`;
`continuous/galerkin_spectral/basis_space.py:144`→`/theory/transport_solver_protocol` **(already broken — page does not exist)**;
`continuous/galerkin_spectral/basis_space.py:417`→`/skills/algebra-of-record` **(already broken)**;
`continuous/galerkin_spectral/basis_space.py:440`→`/theory/reference_solvers`;
`continuous/mms/moc.py:85`→`/theory/method_of_characteristics`;
`continuous/mms/sn.py:51, 646`→`/theory/discrete_ordinates`;
`continuous/peierls_nystrom/cylinder.py:24`, `geometry.py:11, 113, 257`, `reference.py:6`, `sphere.py:24`→`/theory/peierls_nystrom`;
`continuous/peierls_nystrom/cylinder.py:27`, `sphere.py:25`→`/theory/collision_probability`;
`continuous/singular_eigenfunction/__init__.py:105`, `spectrum.py:129`→`/theory/singular_eigenfunction`;
`continuous/singular_eigenfunction/spectrum.py:132, 482`→`/theory/reference_solvers`;
`orpheus/geometry/boundary/__init__.py:6, 232, 371`, `_realizer.py:69`→`/theory/boundary_conditions`;
`orpheus/geometry/boundary/__init__.py:378`→`/theory/discrete_ordinates`;
`orpheus/geometry/boundary/__init__.py:385`→`/theory/operator_algebra`;
`orpheus/geometry/reduced_operator.py:72, 143`→`/theory/structured_geometry`;
`orpheus/geometry/reduced_operator.py:74`→`/theory/discrete_ordinates`;
`orpheus/sn/sweep/pole_angular_closure.py:160`→`/theory/discrete_ordinates`;
`orpheus/transport/spatial/diamond.py:74`, `linear_discontinuous.py:167`, `scheme.py:153, 730`→`/theory/discrete_ordinates`

### 1e. `tests/**` `:doc:` roles (27) — **all INERT**

Nothing under `tests/` is `automodule`'d (`nexus_extra_source_dirs = ['tests']`
feeds *Nexus*, not autodoc). All 27 are human-facing text that will rot silently:

`tests/cross_method/__init__.py:23`→`/testing/cross_method`;
`tests/cross_method/protocol.py:24, 38, 50, 57, 230, 344`→`/skills/vv-principles` **(broken)**;
`tests/cross_method/protocol.py:35`→`/testing/cross_method`;
`tests/cross_method/protocol.py:39`→`/skills/algebra-of-record` **(broken)**;
`tests/cross_method/test_eigenvalue.py:51, 220, 503`→`/skills/vv-principles` **(broken)**;
`tests/derivations/test_kernels.py:19`→`/verification/reference_solutions`;
`tests/derivations/test_quadrature.py:409`→`/theory/peierls_nystrom`;
`tests/derivations/test_sn_mms_anisotropic_symbolic.py:40`→`/theory/discrete_ordinates`;
`tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py:47`→`/theory/discrete_ordinates`;
`tests/diffusion/test_continuous_reference.py:40`→`/theory/diffusion_1d`;
`tests/diffusion/test_continuous_reference.py:41`→`/verification/reference_solutions`;
`tests/homogeneous/test_continuous_reference.py:28`→`/verification/reference_solutions`;
`tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py:53, 143`→`/.claude/skills/vv-principles/SKILL.md` **(broken — not a docname at all)**;
`tests/sn/sweep/core/test_phase_c_gates.py:4`→`/theory/discrete_ordinates`;
`tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py:48`→`docs/theory/discrete_ordinates` **(broken — REL form with `docs/` prefix; would never resolve)**;
`tests/sn/verification/mms/test_mms_2d.py:13`, `test_mms_aniso.py:15`, `test_mms_curvilinear.py:33`, `test_mms_heterogeneous.py:20`→`/theory/discrete_ordinates`

### 1f. Already-broken `:doc:` targets (17) — pre-existing, all INERT

| Target | Sites | Why broken |
|---|---:|---|
| `/skills/vv-principles` | 10 | `docs/skills/` does not exist — skills are `.md` outside the Sphinx tree |
| `/skills/algebra-of-record` | 3 | same |
| `/.claude/skills/vv-principles/SKILL.md` | 2 | not a docname |
| `/theory/transport_solver_protocol` | 1 | page never existed / was never created |
| `docs/theory/discrete_ordinates` (REL in a docstring) | 1 | REL from a non-docs source has no defined base |

All 17 sit in **non-autodoc'd** modules → `-W` has never flagged them. This is
direct proof that the inert-role class rots undetected.

---

## 2. Every toctree entry — the ownership tree

**14 toctree directives. 62 entries. Zero `:glob:`. Zero `:hidden:`. Zero
`:numbered:`.** (Verified by full-corpus grep — no glob hazard exists.)

| # | Owner file | line | options | entries |
|---|---|---:|---|---|
| 1 | `docs/index.rst` | 4 | `:maxdepth: 2`, `:caption: Theory & Derivations` | `theory/index` |
| 2 | `docs/index.rst` | 10 | `:maxdepth: 2`, `:caption: Architecture` | `architecture/index` |
| 3 | `docs/index.rst` | 16 | `:maxdepth: 2`, `:caption: Testing & Verification` | `testing/index`, `verification/index` |
| 4 | `docs/index.rst` | 23 | `:maxdepth: 2`, `:caption: Development` | `development` |
| 5 | `docs/index.rst` | 29 | `:maxdepth: 2`, `:caption: API Reference` | `api/numerics`, `api/transport`, `api/data`, `api/geometry`, `api/homogeneous`, `api/discrete_ordinates`, `api/collision_probability`, `api/method_of_characteristics`, `api/monte_carlo`, `api/diffusion_1d`, `api/fuel_behaviour`, `api/thermal_hydraulics`, `api/reactor_kinetics`, `api/derivations` |
| 6 | `docs/theory/index.rst` | 69 | `:maxdepth: 1` (**no caption**) | `boundary_conditions`, `cross_section_data`, `discrete_measures`, `frame`, `glossary`, `homogeneous`, `index_convention`, `operator_algebra`, `spherical_harmonics`, `structured_geometry`, `verification` |
| 7 | `docs/theory/index.rst` | 85 | `:maxdepth: 1`, `:caption: Discrete (production) solver theory` | `transport_methods`, `diffusion_1d`, `fuel_behaviour`, `thermal_hydraulics`, `reactor_kinetics` |
| 8 | `docs/theory/index.rst` | 96 | `:maxdepth: 1`, `:caption: Continuous (reference) solver theory` | `reference_solvers` |
| 9 | `docs/theory/transport_methods.rst` | 20 | `:maxdepth: 2` | `collision_probability`, `discrete_ordinates`, `loss_representations`, `method_of_characteristics`, `monte_carlo` |
| 10 | `docs/theory/reference_solvers.rst` | 307 | `:maxdepth: 1` | `peierls`, `peierls_nystrom`, `trajectory_resolvent`, `fn_method`, `singular_eigenfunction`, `galerkin_spectral`, `sood_registry` |
| 11 | `docs/theory/reference_solvers.rst` | 327 | `:maxdepth: 1` | `spectral_resolvent`, `pn_method`, `spn_method`, `escape_probability`, `bn_method`, `galerkin_sn_hybrid`, `spectral_collocation` |
| 12 | `docs/testing/index.rst` | 16 | `:maxdepth: 2` | `architecture`, `cross_method` |
| 13 | `docs/verification/index.rst` | 157 | `:maxdepth: 2` | `reference_solutions`, `matrix` |
| 14 | `docs/architecture/index.rst` | 10 | `:maxdepth: 2` | `layering` |

**Every toctree entry is a bare relative name.** All 62 are resolved against the
owning file's directory. Consequence: moving `transport_methods.rst` into
`methods/` silently re-bases all 5 of its children to `methods/` too — which is
*probably what you want*, but it means the toctree entries do **not** need editing
if hub and children move together. Conversely, moving a **child without its hub**
requires editing the hub's toctree. **Design the move map hub-first.**

**Toctree entries pointing at a non-existent doc: zero.** The toctree graph is
currently perfectly consistent.

---

## 3. `:ref:` roles and label definitions

### 3a. Totals

| Metric | Value |
|---|---:|
| `:ref:` usages in `docs/**/*.rst` | **1076** |
| `:ref:` usages in `orpheus/**` + `tests/**` `.py` | **73** (19 live / 54 inert) |
| **Defined labels** (`.. _name:`) in `docs/` | **522** |
| Duplicate label names | **0** |
| Hyperlink targets (`.. _name: url`) | **0** |
| Unresolved `:ref:` **inside `docs/`** | **3 — all Sphinx built-ins** |

**The 3 "unresolved" are `genindex` / `modindex` / `search`** at
`docs/index.rst:58/59/60` — Sphinx built-in targets, **not broken**.

> **The corpus has ZERO genuinely-broken `:ref:` inside `docs/`.** That is a
> clean baseline and it is worth protecting. And because `:ref:` resolves by
> **global label name, not path**, **all 1076 are immune to `git mv`.**

### 3b. Broken `:ref:` **outside** `docs/` (pre-existing)

| file:line | target | Status |
|---|---|---|
| `orpheus/numerics/operator.py:2994` | `operator-algebra-adjoint` | **BROKEN** — no such label. Near-misses: `bc-extraction-g-adjoint` (`operator_algebra.rst:6409`), `eager-adjoint-behavior-change` (`operator_algebra.rst:1337`). Module **not** autodoc'd → inert → never warned. |
| `orpheus/sn/boundary/realizer.py:137` | `bc-tensor-decomposition` | **BROKEN** — near-miss `bc-tensor-primitives` (`operator_algebra.rst:3641`). A rename victim. Module not autodoc'd → inert. |
| `tests/sn/solve/test_krylov_curvilinear_precond_safety.py:10` | `lessons-l19` | inert (agent-memory anchor, not a doc label) |
| `tests/sn/operators/test_invertible_operator.py:645` | `lessons-l18` | inert |
| `tests/sn/verification/analytical/test_phase_c_crosscheck.py:415, 424` | `vv-principles`, `algebra-of-record` | inert (skill names) |
| `tests/sn/sweep/core/test_sweep_graph.py:7, 9, 11, 13, 16` | `assert_upwind_orientation`, `assert_topologically_sorted`, `assert_cell_coverage`, `assert_face_pairing_consistent`, `apply_invariants` | inert (these are *function* names — should be `:func:`) |
| `tests/sn/sweep/curvilinear/test_coupled_pole_mu_level_invariant.py:63` | *(empty)* | inert; malformed `` :ref:`` `` |
| `tests/geometry/_generate_bc_equivalence_snapshots.py:25` | `vv-principles` | inert |
| `tests/cp/test_cylinder_pss.py:22` | `peierls-cyl-Pss-homogeneous` | inert |

The 19 **live** (autodoc-rendered) `:ref:` all resolve correctly today:
`orpheus/cp/solver.py:29`, `data/energy_grid.py:5,109`,
`derivations/common/kernels.py:389`, `common/quadrature.py:12,299`,
`common/quadrature_recipes.py:14,138`, `common/xs_library.py:334`,
`fuel/solver.py:16`, `homogeneous/solver.py:32`, `kinetics/solver.py:25`,
`mc/solver.py:19`, `moc/solver.py:12`, `sn/loss_representation/__init__.py:3252`,
`sn/operators/boundary.py:47`, `sn/operators/streaming.py:98`,
`sn/solver.py:30`, `thermal_hydraulics/solver.py:17`.

### 3c. Label census grouped by DEFINING PAGE (the mis-filing oracle)

| Defining page | #labels | prefix histogram |
|---|---:|---|
| `theory/discrete_ordinates.rst` | 109 | `sn`:66 `ld`:13 `sweep`:7 `balance`:3 `bc`:3 `phase`:2 `si`:2 `theory`:1 `quadrature`:1 |
| `theory/operator_algebra.rst` | 73 | `bc`:14 `sn`:14 `wave`:5 `carrier`:4 `design`:3 `operator`:2 `affine`:2 `wavefront`:2 `green`:2 |
| `theory/fn_method.rst` | 37 | `fn`:36 `theory`:1 |
| `theory/peierls_nystrom.rst` | 36 | `peierls`:15 `theory`:11 `section`:5 `cp`:2 `app`:2 `part`:1 |
| `theory/boundary_conditions.rst` | 35 | `bc`:25 `sn`:8 `theory`:1 `affine`:1 |
| `theory/singular_eigenfunction.rst` | 34 | `theory`:33 `spectrum`:1 |
| `theory/frame.rst` | 22 | `sn`:12 `frame`:8 `galerkin`:1 `petrov`:1 |
| `theory/trajectory_resolvent.rst` | 20 | `peierls`:16 `theory`:2 `billiards`:1 `trajectory`:1 |
| `theory/loss_representations.rst` | 19 | `loss`:19 |
| `theory/collision_probability.rst` | 14 | `peierls`:3 `why`:2 + 9 singletons |
| `theory/diffusion_1d.rst` | 13 | `diffusion`:12 `theory`:1 |
| `theory/homogeneous.rst` | 13 | 9 distinct prefixes, all ≤3 |
| `theory/galerkin_spectral.rst` | 11 | `carlvik`:8 `theory`:2 `galerkin`:1 |
| `theory/sood_registry.rst` | 11 | `sood`:10 `theory`:1 |
| `theory/index_convention.rst` | 10 | `sn`:7 `theory`:3 |
| `verification/reference_solutions.rst` | 9 | `verification`:4 + 5 singletons |
| `theory/cross_section_data.rst` | 8 | `emission`:3 + 5 singletons |
| `theory/verification.rst` | 6 | 6 singletons |
| `theory/reference_solvers.rst` | 4 | `reference`:2 `theory`:1 `orbit`:1 |
| `verification/index.rst` | 4 | `verification`:4 |
| `theory/method_of_characteristics.rst` | 3 | `theory`:1 `flat`:1 `moc`:1 |
| `theory/monte_carlo.rst` | 3 | `mc`:2 `theory`:1 |
| `theory/thermal_hydraulics.rst` | 3 | `theory`:1 `bdf`:1 `th`:1 |
| `theory/discrete_measures.rst` | 2 | `discrete`:2 |
| `theory/index.rst` | 2 | `theory`:2 |
| `theory/peierls.rst` | 2 | `theory`:1 `peierls`:1 |
| `theory/spectral_resolvent.rst` | 2 | `theory`:1 `spectral`:1 |
| `theory/spherical_harmonics.rst` | 2 | `spherical`:2 |
| `testing/architecture.rst` | 2 | `vv`:2 |
| `theory/{bn_method, escape_probability, fuel_behaviour, galerkin_sn_hybrid, glossary, pn_method, reactor_kinetics, spectral_collocation, spn_method, structured_geometry, transport_methods}.rst` | 1 each | `theory`:1 (page anchor only) |
| `architecture/layering.rst` | 1 | `architecture`:1 |
| `api/numerics.rst` | 1 | `field`:1 |
| **TOTAL** | **522** | |

### 3d. Global prefix census (first token)

`sn-`:107 · `theory-`:79 · `bc-`:42 · `fn-`:36 · `peierls-`:35 · `loss-`:19 ·
`diffusion-`:13 · `ld-`:13 · `sood-`:10 · `frame-`:8 · `carlvik-`:8 ·
`verification-`:8 · `sweep-`:7 · `wave-`:5 · `section-`:5 · `emission-`:4 ·
`carrier-`:4 · `reference-`:4 · `affine-`/`cp-`/`balance-`/`direct-`/`operator-`/`design-`/`vv-`:3 each · 20 more at ≤2.

### 3e. **Mis-filing oracle — prefixes SPLIT across pages**

| Prefix | total | split across |
|---|---:|---|
| `sn-*` | 107 | `discrete_ordinates`=66, **`operator_algebra`=14**, **`frame`=12**, **`boundary_conditions`=8**, **`index_convention`=7** |
| `bc-*` | 42 | `boundary_conditions`=25, **`operator_algebra`=14**, `discrete_ordinates`=3 |
| `peierls-*` | 35 | `trajectory_resolvent`=16, `peierls_nystrom`=15, `collision_probability`=3, `peierls`=1 |
| `theory-*` | 79 | **32 pages** — but 33 of the 79 are on `singular_eigenfunction` alone, 11 on `peierls_nystrom` |
| `diffusion-*` | 13 | `diffusion_1d`=12, `verification`=1 |
| `verification-*` | 8 | `verification/index`=4, `verification/reference_solutions`=4 |
| `emission-*` | 4 | `cross_section_data`=3, `operator_algebra`=1 |
| `reference-*` | 4 | `reference_solvers`=2, `verification`=1, `verification/reference_solutions`=1 |

**Readings for the filing decision:**

1. **`theory-*` has two incompatible meanings.** For 30 pages it is the
   *page-top anchor* (`theory-monte-carlo` = "the Monte Carlo page"). But
   `singular_eigenfunction.rst` uses `theory-` as its **section** prefix for 33
   labels, and `peierls_nystrom.rst` for 11. Those two pages are off-pattern
   (cf. the naming-consistency preference for one word-order per family). A
   `theory-` label on a page whose own anchor is *also* `theory-…` is ambiguous.
   **Not a blocker for the move — but this is the one prefix that cannot serve as
   a filing oracle**, because it means "any theory page".
2. **`bc-*` genuinely spans two pages by design.** `boundary_conditions.rst`
   holds the *law/realizer* layers (25); `operator_algebra.rst` holds the
   **`bc-extraction-*` family** (14: `bc-extraction`, `-design-corrections`,
   `-variadic-driver`, `-two-routes`, `-reflect-trace`, `-scope`,
   `-scope-future`, `-2d-si-krylov-twin`, `-numerical-evidence`,
   `-operator-output-typing`, `-operator-output-o2`, `-g-adjoint`, plus
   `bc-tensor-primitives`, `bc-descriptor-tree-vs-operator-tree`). That is
   **correct filing** — BC-extraction *is* an operator-algebra concept (the `B`
   sibling), not a boundary-law concept. If `boundary_conditions.rst` and
   `operator_algebra.rst` land in **different** subdirectories, this 14-label
   family is the one to sanity-check.
3. **`sn-*` scattered over 5 pages is the real signal.** 41 of the 107 `sn-`
   labels live *outside* `discrete_ordinates.rst`. If the target layout puts
   `discrete_ordinates` under `methods/` and `operator_algebra`/`frame`/
   `index_convention` under `foundations/` or `conventions/`, then **38% of the
   `sn-` namespace crosses the new branch boundary.** That is either (a) fine —
   `sn-` is a *topic* tag, not a *location* tag — or (b) a hint that some of those
   14 `operator_algebra` `sn-*` labels are mis-filed. **Worth one review pass;
   not a move blocker** (refs are path-immune).
4. **`peierls-*` split 16/15 across `trajectory_resolvent`/`peierls_nystrom` is
   deliberate** — the documented three-page Peierls architecture. Keep those two
   pages in the same subdirectory to preserve the reading order.

---

## 4. Raw path strings that name a docs page — **the silent class**

**837 total matches** (worktree excluded): 760 name a **specific page**, 77 name
only a **directory** (`docs/theory/` — survives, since the directory still
exists).

### 4a. Bucket census

| Bucket | specific-page | dir-only | already-STALE |
|---|---:|---:|---:|
| `CLAUDE.md` + root | 2 | 2 | 0 |
| `.claude/rules/` | 2 | 1 | 0 |
| `.claude/hooks/` | 1 | 0 | 0 |
| `.claude/agents/` | 9 | 4 | 0 |
| `.claude/skills/` | 20 | 8 | **2** |
| `docs/` (prose in `.rst`) | 7 | 7 | 0 |
| `orpheus/` | 16 | 0 | **1** |
| `tests/` | 34 | 9 | **4** |
| `tools/` | 6 | 5 | 0 |
| **ACTIONABLE SUBTOTAL** | **97** | **36** | **7** |
| `.claude/plans/` | 160 | 11 | 23 |
| `.claude/agent-memory/` | 337 | 29 | 16 |
| `.claude/scratch/` (tracked) | 166 | 1 | 6 |
| **ARCHAEOLOGY SUBTOTAL** | **663** | **41** | **45** |
| **TOTAL** | **760** | **77** | **52** |

**Recommendation on scope:** fix the **97 actionable** rows. `.claude/plans/`,
`.claude/agent-memory/`, `.claude/scratch/` (663 rows) are **archaeology** —
frozen point-in-time records whose whole purpose is to describe the tree *as it
was*. Rewriting history-notes to point at today's paths makes them lie about
their own moment. **Leave them.** (This matches the project's own
process-discipline stance that a campaign note's terminal state is "merged @
hash", after which it is archaeology.)

### 4b. **THE PROOF THAT THIS HAZARD IS REAL: 52 raw paths already point at pages that do not exist**

`git log --diff-filter=R -- docs/` shows the rename history, and each rename left
orphans behind:

| Rename | commit | raw refs left dangling |
|---|---|---:|
| `peierls_greens.rst` → `trajectory_resolvent.rst` | `d7fa25bf` | **16** |
| `carlvik_galerkin.rst` → `galerkin_spectral.rst` | `d7fa25bf` | 2 |
| `_peierls_capability_matrix.inc.rst` → `_peierls_nystrom_capability_matrix.inc.rst` | `e2f2e829` | 3 |
| `peierls_unified.rst` → `peierls_nystrom.rst` | `1a1c7c85` | — |
| **`galerkin_projection.rst` → `frame.rst`** | **`3de597a3` — THIS BRANCH, 2 days ago** | **2** |

Plus never-created / long-dead pages still named in prose: `sn_dim_agnostic.rst`
(4), `sn_verification_matrix.rst` (4), `geometry.rst` (3), `sn_adjoint.rst` (3),
`sn_operator_algebra.rst` (3), `case_method.rst` (3), `diffusion.rst` (2),
`transport_solver_protocol.rst` (1), `greens_functions.rst` (1),
`peierls_nystrom_advanced.rst` (1), `registry.rst` (1),
`verification/agreement.rst` (1), `theory/capabilities.rst` (1).

**Every one of these was left by a `git mv` exactly like the one being planned.**
The corpus's own history is the evidence that a build gate never catches this.

### 4c. The 7 already-stale rows in ACTIONABLE buckets (fix opportunistically)

| file:line | names | reality |
|---|---|---|
| `.claude/skills/algebra-of-record/SKILL.md:731` | `docs/theory/peierls_greens.rst` | → `trajectory_resolvent.rst` |
| `.claude/skills/subagent-handoff-protocol/SKILL.md:511` | `docs/theory/peierls_greens.rst` | → `trajectory_resolvent.rst` |
| `orpheus/transport/mesh/axis.py:27` | `docs/theory/sn_dim_agnostic.rst` | never created ("lands in C8") |
| `tests/sn/regression/_generate_snapshots.py:4` | `docs/testing/sn_verification_matrix.rst` | never created |
| `tests/sn/regression/test_dd_regression.py:41` | `docs/testing/sn_verification_matrix.rst` | never created |
| `tests/sn/verification/analytical/README.md:9` | `docs/testing/sn_verification_matrix.rst` | never created |
| `tests/sn/verification/analytical/test_kinf_homogeneous.py:12` | `docs/testing/sn_verification_matrix.rst` | never created |

### 4d. ACTIONABLE specific-page raw paths — full list (97)

**`CLAUDE.md` (2)**
- `CLAUDE.md:49` — ``` `docs/theory/` contain theory, derivations… ``` *(dir-only — survives)*
- `CLAUDE.md:150` — ``` See `docs/development.rst` for the full workflow ``` *(survives unless `development.rst` moves)*

**`.claude/rules/` (2 specific + 1 dir)**
- `.claude/rules/vv-testing.md:43` — `docs/testing/architecture.rst` **← moves if `docs/testing/` consolidates**
- `.claude/rules/vv-testing.md:50` — `` in `docs/theory/` `` *(dir-only)*
- `.claude/rules/vv-testing.md:57` — `docs/verification/matrix.rst` **← moves if `docs/verification/` consolidates**

**`.claude/hooks/` (1)**
- `.claude/hooks/session-start.txt:38` — `Read docs/development.rst`

**`.claude/agents/` (9)** — **note `explorer/AGENT.md` names 5 theory pages in a literal path block**
- `.claude/agents/archivist/AGENT.md:110` — `docs/conf.py`
- `.claude/agents/archivist/AGENT.md:114` — `docs/conf.py`
- `.claude/agents/archivist/AGENT.md:468` — `api/numerics.rst`
- `.claude/agents/explorer/AGENT.md:107` — `docs/theory/operator_algebra.rst`
- `.claude/agents/explorer/AGENT.md:162` — `docs/theory/discrete_ordinates.rst`
- `.claude/agents/explorer/AGENT.md:163` — `docs/theory/collision_probability.rst`
- `.claude/agents/explorer/AGENT.md:164` — `docs/theory/homogeneous.rst`
- `.claude/agents/explorer/AGENT.md:165` — `docs/theory/method_of_characteristics.rst`
- `.claude/agents/explorer/AGENT.md:166` — `docs/theory/monte_carlo.rst`
- *(also `.claude/agents/elegance-enforcer/AGENT.md:42`, `method-implementer/AGENT.md:156,183`, `archivist/AGENT.md:112` — dir-only `docs/theory/`)*

**`.claude/skills/` (20)**
- `algebra-of-record/SKILL.md:731` — `docs/theory/peierls_greens.rst` **STALE**
- `algebra-of-record/SKILL.md:780` — `docs/conf.py`
- `doc-issue-relocation/SKILL.md:113` — `docs/theory/peierls_nystrom.rst`
- `doc-issue-relocation/SKILL.md:118` — `docs/conf.py`
- `doc-issue-relocation/SKILL.md:143` — `docs/theory/peierls_nystrom.rst`
- `doc-issue-relocation/SKILL.md:256` — `verification/matrix.rst`
- `doc-issue-relocation/SKILL.md:272` — `docs/theory/peierls_nystrom.rst` *(+ explicit line range 5410–5476)*
- `nexus-cli/SKILL.md:58` — `docs/conf.py`
- `nexus-migration/SKILL.md:29` — `theory/collision_probability.rst`
- `numerical-bug-signatures/SKILL.md:296` — `docs/theory/fn_method.rst`
- `numerical-bug-signatures/SKILL.md:601` — `docs/theory/diffusion_1d.rst`
- `subagent-handoff-protocol/SKILL.md:511` — `docs/theory/peierls_greens.rst` **STALE**
- `vv-principles/SKILL.md:128` — `docs/theory/homogeneous.rst`
- `vv-principles/error_catalog.md:176` — `docs/theory/discrete_ordinates.rst`
- `vv-principles/error_catalog.md:1054` — `docs/theory/diffusion_1d.rst`
- `vv-principles/error_catalog.md:2109` — `docs/theory/peierls_nystrom.rst`
- `vv-principles/error_catalog.md:4477` — `docs/theory/discrete_ordinates.rst`
- `vv-principles/error_catalog.md:4575` — `docs/theory/peierls_nystrom.rst`
- `vv-principles/scripts/err006_convergence_to_wrong_limit.md:48` — `docs/theory/discrete_ordinates.rst`
- `vv-principles/scripts/issue226_spectral_invisibility.md:57` — `docs/theory/homogeneous.rst`

**`docs/` prose (7)**
- `docs/conf.py:90` — `docs/verification/matrix.rst` *(comment)*
- `docs/theory/peierls.rst:719` — `:file:`docs/theory/_peierls_nystrom_capability_matrix.inc.rst``
- `docs/theory/peierls_nystrom.rst:6286` — `See docs/theory/peierls_nystrom.rst` *(self-reference in prose)*
- `docs/theory/peierls_nystrom.rst:6397` — `See docs/theory/peierls_nystrom.rst §12`
- `docs/theory/verification.rst:436` — ``` ``docs/testing/architecture.rst`` ```
- `docs/verification/matrix.rst:836` — ``` ``docs/testing/architecture.rst`` ``` **(auto-generated — fix the generator, not the file)**
- `docs/verification/matrix.rst:1220` — ``` ``docs/testing/architecture.rst`` ``` **(auto-generated)**

**`orpheus/` (16)**
- `orpheus/data/micro_xs/gendf.py:346` — `docs/theory/cross_section_data.rst`
- `orpheus/derivations/README.md:131` — `docs/theory/verification.rst`
- `orpheus/derivations/README.md:134` — `docs/verification/reference_solutions.rst`
- `orpheus/derivations/README.md:137` — `docs/api/derivations.rst`
- `orpheus/derivations/common/transport_equation.py:80` — `docs/theory/verification.rst`
- `orpheus/derivations/continuous/cases/diffusion.py:133` — `docs/theory/diffusion_1d.rst`
- `orpheus/derivations/continuous/cases/diffusion.py:167` — `docs/theory/diffusion_1d.rst`
- `orpheus/derivations/continuous/cases/sn.py:160` — `docs/theory/diffusion_1d.rst`
- `orpheus/derivations/continuous/mms/sn.py:329` — `docs/theory/discrete_ordinates.rst`
- `orpheus/derivations/continuous/mms/sn.py:1403` — `docs/theory/discrete_ordinates.rst`
- `orpheus/derivations/continuous/mms/sn.py:2761` — `docs/theory/discrete_ordinates.rst`
- `orpheus/derivations/continuous/peierls_nystrom/geometry.py:5507` — `docs/theory/peierls_nystrom.rst` **(inside an f-string — a user-visible runtime warning message)**
- `orpheus/derivations/continuous/peierls_nystrom/origins/cylinder_g_bc_3d.py:100` — `docs/theory/peierls_nystrom.rst:peierls-cyl-Gbc-3d-final`
- `orpheus/derivations/discrete/moc/equations.py:4` — `docs/theory/method_of_characteristics.rst`
- `orpheus/derivations/discrete/sn/balance.py:9` — `docs/theory/discrete_ordinates.rst`
- `orpheus/transport/mesh/axis.py:27` — `docs/theory/sn_dim_agnostic.rst` **STALE**

**`tests/` (34)**
- `tests/_harness/__init__.py:3`, `tests/_harness/registry.py:35`, `tests/conftest.py:14`, `tests/conftest.py:78`, `tests/data/test_cross_section_data.py:22`, `tests/geometry/test_geometry.py:49` — `docs/testing/architecture.rst`
- `tests/cp/test_cylinder.py:22`, `tests/cp/test_properties.py:30` — `docs/theory/collision_probability.rst`
- `tests/cp/test_peierls_rank_n_protocol.py:37, 52`, `tests/derivations/test_peierls_rank2_bc.py:475`, `tests/derivations/test_peierls_rank_n_bc.py:10, 39`, `tests/derivations/test_peierls_rank_n_class_b_mr_mg.py:515` — `docs/theory/peierls_nystrom.rst` *(the last is inside an **f-string assert message**)*
- `tests/data/test_cross_section_data.py:24`, `tests/homogeneous/test_homogeneous.py:90` — `docs/theory/homogeneous.rst`
- `tests/data/test_energy_grid.py:10`, `tests/data/test_gendf_canonical_order.py:11`, `tests/mc/test_gaps.py:594` — `docs/theory/cross_section_data.rst`
- `tests/data/test_mixture.py:3` — `docs/theory/homogeneous.rst:622` **(carries a LINE NUMBER — already fragile)**
- `tests/data/test_mixture_transport_xs.py:14` — `docs/theory/diffusion_1d.rst`
- `tests/derivations/test_sn_mms_nonvacuum_symbolic.py:40`, `tests/sn/operators/test_invertible_operator.py:687`, `tests/sn/solve/test_krylov_curvilinear_precond_safety.py:64`, `tests/sn/verification/analytical/test_mms_prescribed_inflow.py:43`, `tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py:26`, `tests/sn/verification/mms/test_mms.py:11`, `tests/sn/verification/mms/test_space_angle_separability.py:23`, `tests/transport/spatial/test_ld_slope_frame.py:31` — `docs/theory/discrete_ordinates.rst`
- `tests/numerics/test_spherical_harmonic_space.py:10` — `docs/theory/spherical_harmonics.rst`
- `tests/sn/regression/_generate_snapshots.py:4`, `tests/sn/regression/test_dd_regression.py:41`, `tests/sn/verification/analytical/README.md:9`, `tests/sn/verification/analytical/test_kinf_homogeneous.py:12` — `docs/testing/sn_verification_matrix.rst` **STALE ×4**

**`tools/` (6)**
- `tools/verification/generate_matrix.py:1` — `docs/verification/matrix.rst` *(module docstring)*
- `tools/verification/generate_matrix.py:21` — `docs/testing/architecture.rst`
- `tools/verification/generate_matrix.py:24` — `docs/conf.py`
- `tools/verification/generate_matrix.py:32` — `docs/verification/matrix.rst`
- `tools/verification/generate_matrix.py:221` — `docs/testing/architecture.rst` **(emitted INTO the generated matrix.rst)**
- `tools/verification/generate_matrix.py:263` — `docs/testing/architecture.rst` **(emitted INTO the generated matrix.rst)**

**`pyproject.toml` (1)**
- `pyproject.toml:88` — `"verifies(*labels): Sphinx equation labels (:label: in docs/theory/*.rst) that this test exercises"` *(glob, dir-level — survives; but see §7 for the scanner)*

### 4e. `:file:` roles (156) — a second silent class

`:file:` renders as literal text and produces **no link and no warning**. 156
occurrences, concentrated in `trajectory_resolvent.rst` (77),
`discrete_ordinates.rst` (31), `peierls_nystrom.rst` (30),
`index_convention.rst` (5), `architecture/layering.rst` (4). Almost all point at
`tests/…`, `scratch/…`, `.claude/…` — **not** at docs pages, so the restructure
barely touches them. The one that does: `docs/theory/peierls.rst:719` →
`:file:`docs/theory/_peierls_nystrom_capability_matrix.inc.rst`` (see §5/§7).

---

## 5. `.. include::` directives and media

### 5a. `.. include::` — 10 total, **all relative, all in `docs/theory/`**

| Source | line | include target | breaks if… |
|---|---:|---|---|
| `docs/theory/verification.rst` | 293 | `../_generated/verification_table.rst` | **`verification.rst` changes depth** |
| `docs/theory/verification.rst` | 314 | `../_generated/homogeneous_derivation.rst` | ″ |
| `docs/theory/verification.rst` | 344 | `../_generated/sn_derivation.rst` | ″ |
| `docs/theory/verification.rst` | 363 | `../_generated/cp_slab_derivation.rst` | ″ |
| `docs/theory/verification.rst` | 382 | `../_generated/cp_cylinder_derivation.rst` | ″ |
| `docs/theory/verification.rst` | 461 | `../_generated/moc_derivation.rst` | ″ |
| `docs/theory/verification.rst` | 479 | `../_generated/mc_derivation.rst` | ″ |
| `docs/theory/verification.rst` | 546 | `../_generated/diffusion_derivation.rst` | ″ |
| `docs/theory/fn_method.rst` | 65 | `_fn_method_capability_matrix.inc.rst` | **`fn_method.rst` moves** (§7) |
| `docs/theory/peierls_nystrom.rst` | 223 | `_peierls_nystrom_capability_matrix.inc.rst` | **`peierls_nystrom.rst` moves** (§7) |

**`verification.rst` is the single most move-sensitive page in the corpus** — it
carries 8 relative `../_generated/` includes. Moving it to
`docs/theory/foundations/verification.rst` requires all 8 → `../../_generated/`.
Moving it to `docs/verification/verification.rst` requires all 8 → `../_generated/`
(**unchanged, correct by accident** — do not let that coincidence hide the issue).

### 5b. `literalinclude` / `figure` / `image` / `csv-table :file:`

**ZERO.** No `.. literalinclude::`, `.. figure::`, or `.. image::` directive
exists anywhere in `docs/`. No media paths to remap. (The `:file:` hits in §4e are
the *inline role*, not the directive.)

### 5c. `.. plot::` (matplotlib plot_directive)

`matplotlib.sphinxext.plot_directive` is loaded in `conf.py:27`. Any `.. plot::`
with an external script path would be move-sensitive — none found with a path
argument.

---

## 6. `docs/conf.py` — every path-valued setting

| line | setting | value | move impact |
|---:|---|---|---|
| 8–9 | `project_root` / `sys.path.insert` | `Path(conf.py).parent.parent` | none (repo root) |
| 37 | `bibtex_bibfiles` | `['refs.bib']` | none (`docs/refs.bib` stays) |
| **47** | **`templates_path`** | **`['_templates']`** | **`docs/_templates` DOES NOT EXIST.** Sphinx warns for `html_static_path`/`html_extra_path`/`logo`/`favicon` but **not** for `templates_path` → silent. Pre-existing, harmless, but it is a config path that does not resolve. |
| **48** | **`exclude_patterns`** | **`['_build', '_generated', 'Thumbs.db', '.DS_Store']`** | **SEE 6a — this is effectively a NO-OP** |
| 53 | `html_static_path` | `['_static']` | resolves (`sortable.css`, `sortable.js`); unaffected |
| 54–55 | `html_js_files` / `html_css_files` | `['sortable.js']` / `['sortable.css']` | relative to `_static`; unaffected |
| **61** | **`nexus_extra_source_dirs`** | **`['tests']`** | repo-relative — unaffected by a `docs/` move. **This is why `tests/` docstrings are Nexus nodes but NOT autodoc pages** (§1e). |
| **71** | **`nexus_source_exclude_patterns`** | **`['student_resources/*']`** | source-tree pattern, not a docs pattern — unaffected |
| 93–111 | `_regenerate_verification_matrix` | `subprocess … tools.verification.generate_matrix`, `cwd=project_root` | **fires on `builder-inited` — writes `docs/verification/matrix.rst`** (§7) |
| 124–137 | `_regenerate_capability_matrices` | `subprocess … tools.verification.generate_capability_matrices`, `cwd=project_root` | **fires on `builder-inited` — writes `docs/theory/_*.inc.rst`** (§7) |
| 140–142 | `setup(app)` | connects both hooks to `builder-inited` | both run **before** Sphinx reads sources |

### 6a. **`exclude_patterns` is a no-op for `_generated`** — verified against Sphinx 9.1.0

```
_translate_pattern('_build')     -> '_build$'
_translate_pattern('_generated') -> '_generated$'

Matcher(exclude_patterns)('_generated/sn_derivation.rst')  ->  False   # NOT excluded
Matcher(exclude_patterns)('_generated')                    ->  True
Matcher(exclude_patterns)('_build/html/x.rst')             ->  False   # NOT excluded
```

`'_generated'` compiles to the regex `_generated$`, which matches only the bare
docname `_generated` — **not** `_generated/sn_derivation.rst`. So **all 8
`docs/_generated/*.rst` files ARE read by Sphinx as first-class source
documents.** They escape the "document isn't included in any toctree" warning
**only** because `docs/theory/verification.rst` `.. include::`s all 8 of them.

(`'_build'` is likewise a no-op, but harmless: Sphinx auto-excludes the outdir
when it lives inside srcdir.)

**Implication:** the include chain in §5a is not just a convenience — it is the
**sole** thing keeping 8 files out of the orphan-warning list under `-W`. If the
restructure moves `verification.rst` and any include is dropped or mis-pathed,
you get *both* an include error *and* up to 8 orphan warnings.

**Correct fix (independent of the move):** `exclude_patterns` should be
`['_build', '_generated/*', 'Thumbs.db', '.DS_Store', 'theory/_*.inc.rst']` —
or better, `**/_generated/**`. Worth doing in the same commit; it makes the
`_generated` and `.inc.rst` handling robust rather than accidental.

### 6b. The orphan-suppression mechanism (Sphinx 9.1.0 `BuildEnvironment.check_consistency`)

```python
included = set().union(*self.included.values())
for docname in sorted(self.all_docs):
    if docname not in self.files_to_rebuild:
        if docname == self.config.root_doc:   continue
        if docname in included:               continue      # <-- THE SUPPRESSOR
        if 'orphan' in self.metadata[docname]: continue
        logger.warning("document isn't included in any toctree", ...)
```

`self.included` is populated by `.. include::`. **This is load-bearing for 10
files** (8 `_generated`, 2 `.inc.rst`). None of the 10 carries an `:orphan:`
metadata field — verified.

---

## 7. **NEW CLASS — Python `Path` constants (the #1 hazard)**

The brief asked for raw path *strings*. There is a fourth class that a
`grep 'docs/theory/'` **structurally cannot find**, because the path is built
from **segments**:

```python
DOCS_THEORY_DIR = REPO_ROOT / "docs" / "theory"      # contains no "docs/theory/"
```

Exhaustive search for `"docs"` as a Python path segment (whole repo, worktree
excluded) — **5 hits, 4 of them live couplings:**

| file:line | constant | consequence of the move |
|---|---|---|
| **`tools/verification/generate_capability_matrices.py:83`** | `DOCS_THEORY_DIR = REPO_ROOT / "docs" / "theory"` | **writes `.inc.rst` to the FLAT path forever** |
| **`tools/verification/generate_capability_matrices.py:84`** | `MATRIX_GLOB = "_*_capability_matrix.inc.rst"` | used with **non-recursive `.glob()`** at lines 301 & 331 |
| **`tools/verification/generate_matrix.py:44`** | `DEFAULT_OUT = REPO_ROOT / "docs" / "verification" / "matrix.rst"` | **breaks if `docs/verification/` consolidates** |
| **`tests/derivations/test_capability_matrices.py:29`** | `DOCS_THEORY_DIR = REPO_ROOT / "docs" / "theory"` | **a `@pytest.mark.foundation` test goes RED** |
| `orpheus/derivations/generate_rst.py:17` | `OUTPUT_DIR = … / "docs" / "_generated"` | safe — `docs/_generated/` is not being moved |

### 7a. The capability-matrix triple-coupling — walk the failure

Suppose `fn_method.rst` → `docs/theory/references/fn_method.rst`:

1. `conf.py:141` fires `_regenerate_capability_matrices` on **`builder-inited`**,
   i.e. *before Sphinx reads any source*.
2. The generator writes to `DOCS_THEORY_DIR / "_fn_method_capability_matrix.inc.rst"`
   = **`docs/theory/_fn_method_capability_matrix.inc.rst`** — the OLD flat path.
   It has no idea the page moved.
3. `docs/theory/references/fn_method.rst:65` says
   `.. include:: _fn_method_capability_matrix.inc.rst`, which resolves relative to
   **`docs/theory/references/`** → **file not found → Sphinx include error → `-W` FAILS.**
4. Meanwhile the flat `docs/theory/_fn_method_capability_matrix.inc.rst` is now
   included by **nothing** → drops out of `env.included` → **"document isn't
   included in any toctree" → a SECOND `-W` failure.**
5. `_wipe_existing_matrices()` (line 301) does
   `DOCS_THEORY_DIR.glob(MATRIX_GLOB)` — **non-recursive**. If you hand-move the
   `.inc.rst` into `references/`, the wipe can no longer see it, and the
   generator **re-creates the flat copy** → two divergent files, one stale forever.
6. `tests/derivations/test_capability_matrices.py` (marked
   `@pytest.mark.foundation`) runs the generator as a **subprocess** with
   `--check` and asserts exit 0; `--check` compares rendered paths against
   on-disk and flags "stale on-disk matrices that do not correspond to any
   discovered package" → **RED**.

So one `git mv` of `fn_method.rst` produces: 2 Sphinx `-W` errors + 1 red
foundation test + a silently duplicated generated file. **None of it is visible
to a `:doc:`/`:ref:`/toctree audit.**

**Fix:** the generator must learn the page's new home. Either (a) keep both
`.inc.rst` files at `docs/theory/` and make the moved pages include
`../_fn_method_capability_matrix.inc.rst`, or (b) parameterise
`DOCS_THEORY_DIR` per package and switch the wipe to `rglob`. **(b) is the
principled one** — the generator currently hard-codes a layout assumption that the
restructure invalidates. Whichever is chosen, `tests/derivations/test_capability_matrices.py:29`
must move in lockstep.

### 7b. The label scanner is SAFE — `rglob`

Good news, verified:

| file:line | call | verdict |
|---|---|---|
| `tests/_harness/audit.py:339` | `theory_dir.rglob("*.rst")` | **recursive → survives subdirectories** |
| `tests/_harness/audit.py:373` | `docs_dir.rglob("*.rst")` | **recursive → survives** |
| `tests/_harness/audit.py:426` | `default=Path("docs/theory")` | argparse default; `docs/theory/` still exists |
| `tests/_harness/audit.py:445` | `_scan_all_doc_labels(args.theory_dir.parent)` → `docs/` | **recursive → survives** |

The orphan-equation gate and the `verifies(...)` label census therefore **survive
the restructure untouched**. (They will additionally start seeing the 2
`.inc.rst` files — they already do today, since `rglob` matches `*.inc.rst`.)

---

## 8. The three questions

### Q1 — Does `docs/theory/index.rst` reference labels defined on OTHER pages?

**YES — 2 of the 3.** The `list-table` at `docs/theory/index.rst:12–42` names
three branch labels:

| `:ref:` site | label | **Defined at** | Same page? |
|---|---|---|---|
| `docs/theory/index.rst:20` | `theory-infrastructure` | **`docs/theory/index.rst:60`** | **YES — local**, it is the heading of the page's own "Infrastructure" section |
| `docs/theory/index.rst:26` | `theory-transport-methods` | **`docs/theory/transport_methods.rst:1`** | **NO — another page** |
| `docs/theory/index.rst:33` | `theory-reference-solvers` | **`docs/theory/reference_solvers.rst:1`** | **NO — another page** |

(There is no `theory-foundations` label anywhere; the first branch is called
`theory-infrastructure`. `docs/theory/index.rst` defines exactly 2 labels:
`theory-index` at line 1, `theory-infrastructure` at line 60.)

**The architectural finding this exposes — and it matters for the target layout:**
the index presents a symmetric 3-branch taxonomy, but the branches are **not
structurally symmetric**:

- Branch 1 "Infrastructure" — a **local section on `index.rst`** owning a bare,
  caption-less toctree of 11 pages. **It has no hub page of its own.**
- Branch 2 "Transport methods" — a **real hub page** (`transport_methods.rst`,
  955 bytes) owning 5 children.
- Branch 3 "Reference solvers" — a **real hub page** (`reference_solvers.rst`,
  13.8 KB) owning 14 children across 2 toctrees.

So a `foundations/` (or `conventions/`) subdirectory has **no existing hub to
carry into it** — you must either mint `foundations/index.rst` (and move the
`theory-infrastructure` label onto it, which is free: `:ref:` is path-immune) or
leave branch 1's toctree on `theory/index.rst` pointing into the subdirectory.
**The former is the symmetric choice** and makes the tree self-describing.

Other inbound refs to these labels (all path-immune, no edits needed):
`reference_solvers.rst:342` → `theory-transport-methods`;
`discrete_measures.rst:789` → `theory-reference-solvers`.

### Q2 — The current toctree tree (as it exists today)

```
docs/index.rst                                    [root_doc]
│
├─[toctree :maxdepth:2 :caption: Theory & Derivations]
│  └── theory/index.rst  ......................... (.. _theory-index:)
│      │
│      ├─[toctree :maxdepth:1]  (NO caption — branch "Infrastructure",
│      │                         section .. _theory-infrastructure: is LOCAL at L60)
│      │  ├── theory/boundary_conditions.rst  ......... 165 KB   35 labels
│      │  ├── theory/cross_section_data.rst  ..........  53 KB    8 labels
│      │  ├── theory/discrete_measures.rst  ...........  35 KB    2 labels
│      │  ├── theory/frame.rst  ....................... 180 KB   22 labels
│      │  ├── theory/glossary.rst  ....................   5 KB    1 label
│      │  ├── theory/homogeneous.rst  .................  68 KB   13 labels
│      │  ├── theory/index_convention.rst  ............  72 KB   10 labels
│      │  ├── theory/operator_algebra.rst  ............ 577 KB   73 labels
│      │  ├── theory/spherical_harmonics.rst  .........  24 KB    2 labels
│      │  ├── theory/structured_geometry.rst  .........  21 KB    1 label
│      │  └── theory/verification.rst  ................  25 KB    6 labels
│      │        └── .. include:: ../_generated/*.rst  × 8   <-- MOVE-CRITICAL
│      │
│      ├─[toctree :maxdepth:1 :caption: Discrete (production) solver theory]
│      │  ├── theory/transport_methods.rst  ...........  955 B    1 label   [HUB]
│      │  │   └─[toctree :maxdepth:2]
│      │  │      ├── theory/collision_probability.rst  . 174 KB   14 labels
│      │  │      ├── theory/discrete_ordinates.rst  .... 942 KB  109 labels
│      │  │      ├── theory/loss_representations.rst  .. 146 KB   19 labels
│      │  │      ├── theory/method_of_characteristics.rst 52 KB    3 labels
│      │  │      └── theory/monte_carlo.rst  ...........  47 KB    3 labels
│      │  ├── theory/diffusion_1d.rst  ................  64 KB   13 labels
│      │  ├── theory/fuel_behaviour.rst  ..............  19 KB    1 label
│      │  ├── theory/thermal_hydraulics.rst  ..........  28 KB    3 labels
│      │  └── theory/reactor_kinetics.rst  ............  13 KB    1 label
│      │
│      └─[toctree :maxdepth:1 :caption: Continuous (reference) solver theory]
│         └── theory/reference_solvers.rst  ...........  14 KB    4 labels   [HUB]
│             ├─[toctree :maxdepth:1]
│             │  ├── theory/peierls.rst  ..............   34 KB    2 labels
│             │  ├── theory/peierls_nystrom.rst  ......  401 KB   36 labels
│             │  │     └── .. include:: _peierls_nystrom_capability_matrix.inc.rst
│             │  ├── theory/trajectory_resolvent.rst  .  273 KB   20 labels
│             │  ├── theory/fn_method.rst  ............  109 KB   37 labels
│             │  │     └── .. include:: _fn_method_capability_matrix.inc.rst
│             │  ├── theory/singular_eigenfunction.rst    88 KB   34 labels
│             │  ├── theory/galerkin_spectral.rst  ....   35 KB   11 labels
│             │  └── theory/sood_registry.rst  ........   37 KB   11 labels
│             └─[toctree :maxdepth:1]   (reserved / stub slots)
│                ├── theory/spectral_resolvent.rst  ...    4 KB    2 labels
│                ├── theory/pn_method.rst  ............    2 KB    1 label
│                ├── theory/spn_method.rst  ...........    2 KB    1 label
│                ├── theory/escape_probability.rst  ...    3 KB    1 label
│                ├── theory/bn_method.rst  ............    2 KB    1 label
│                ├── theory/galerkin_sn_hybrid.rst  ...    2 KB    1 label
│                └── theory/spectral_collocation.rst  .    2 KB    1 label
│
├─[toctree :maxdepth:2 :caption: Architecture]
│  └── architecture/index.rst
│      └─[toctree :maxdepth:2]
│         └── architecture/layering.rst  ..............  1 label
│
├─[toctree :maxdepth:2 :caption: Testing & Verification]
│  ├── testing/index.rst
│  │   └─[toctree :maxdepth:2]
│  │      ├── testing/architecture.rst  ...............  2 labels
│  │      └── testing/cross_method.rst  ...............  0 labels
│  └── verification/index.rst  ........................  4 labels
│      └─[toctree :maxdepth:2]
│         ├── verification/reference_solutions.rst  ...  9 labels
│         └── verification/matrix.rst  ................  0 labels  [AUTO-GENERATED]
│
├─[toctree :maxdepth:2 :caption: Development]
│  └── development.rst
│
└─[toctree :maxdepth:2 :caption: API Reference]
   ├── api/numerics.rst  .............................  1 label
   ├── api/transport.rst
   ├── api/data.rst
   ├── api/geometry.rst
   ├── api/homogeneous.rst
   ├── api/discrete_ordinates.rst
   ├── api/collision_probability.rst
   ├── api/method_of_characteristics.rst
   ├── api/monte_carlo.rst
   ├── api/diffusion_1d.rst
   ├── api/fuel_behaviour.rst
   ├── api/thermal_hydraulics.rst
   ├── api/reactor_kinetics.rst
   └── api/derivations.rst

NOT IN ANY TOCTREE (read by Sphinx; orphan-warning suppressed by .. include::):
   ├── theory/_fn_method_capability_matrix.inc.rst        [GENERATED, TRACKED]
   ├── theory/_peierls_nystrom_capability_matrix.inc.rst  [GENERATED, TRACKED]
   └── _generated/*.rst  × 8   [GENERATED, GITIGNORED; exclude_patterns is a NO-OP on them]
```

**Tree facts:** 63 documents read after `exclude_patterns`; 61 reachable from
`docs/index.rst`; 2 orphans. 37 theory pages + 2 `.inc.rst` = the 39 files in
`docs/theory/`.

### Q3 — Orphans / glob-only reachability

**Orphans (built but not in any toctree): exactly 2.**

| file | Why it does not warn today |
|---|---|
| `docs/theory/_fn_method_capability_matrix.inc.rst` | `.. include::`d by `fn_method.rst:65` → lands in `env.included` → warning suppressed |
| `docs/theory/_peierls_nystrom_capability_matrix.inc.rst` | `.. include::`d by `peierls_nystrom.rst:223` → suppressed |

Neither carries an `:orphan:` metadata field; neither is matched by
`exclude_patterns`. Their non-orphan status is **100% contingent on the include
resolving** — see §7a.

*(`docs/_generated/*.rst` — 8 files — are in the same category: read as source,
not in any toctree, suppressed only by `verification.rst`'s includes. They are
gitignored build artefacts, so they don't appear in a `git ls-files` orphan
audit, but Sphinx sees them.)*

**Reachable only via `:glob:`: ZERO.** There is **no `:glob:` option anywhere in
`docs/`** (nor `:hidden:`, nor `:numbered:`). Confirmed by full-corpus grep. The
glob re-scoping hazard named in the brief **does not exist in this corpus** —
and introducing one during the restructure would be a regression, because a
`:glob:` under a new subdirectory silently absorbs anything later dropped there.

---

## 9. Hazard ranking — what actually bites

| # | Hazard | Detected by `-W`? | Severity |
|---|---|---|---|
| **1** | **Capability-matrix triple-coupling** (§7a): hardcoded `DOCS_THEORY_DIR`, non-recursive wipe glob, relative `.. include::`, and a `@pytest.mark.foundation` test pinning the path. Fires on `builder-inited`, *before* Sphinx reads anything. **Invisible to every grep the brief specified.** | 2 errors, but only *after* you've already broken it | **CRITICAL** |
| **2** | **`verification.rst`'s 8 relative `../_generated/` includes** (§5a) + `exclude_patterns` being a no-op (§6a). Moving this one page can produce an include error *and* 8 orphan warnings. | yes | **HIGH** |
| **3** | **132 bare-sibling REL `:doc:` inside `docs/theory/`** (§1). Any two pages landing in different subdirs breaks every REL between them. Fix: convert to ABS. | yes | **HIGH (volume)** |
| **4** | **`tools/verification/generate_matrix.py:44`** hardcodes `docs/verification/matrix.rst`; lines 221/263 **emit `docs/testing/architecture.rst` INTO the generated page**. If `docs/verification/`+`docs/testing/` consolidate, fix the *generator*, not `matrix.rst`. | no (writes fine to a dead path) | **HIGH** |
| **5** | **97 actionable raw path strings + 44 inert `orpheus/` `:doc:` + 27 inert `tests/` `:doc:` + 54 inert `:ref:`** = **~222 silent references.** No build catches any of them. §4b proves 52 are *already* dead from prior renames. | **NO** | **HIGH (silent)** |
| 6 | `.claude/agents/explorer/AGENT.md:162–166` — the explorer agent's own "Sphinx Documentation" block names 5 theory pages by path. A stale map here **mis-steers every future exploration dispatch**. | no | MEDIUM |
| 7 | Pre-existing broken refs: `operator-algebra-adjoint`, `bc-tensor-decomposition` (§3b); 17 broken `:doc:` (§1f). All inert. Free to fix in the same pass. | no | LOW (but free) |
| 8 | `templates_path=['_templates']` → non-existent dir (§6). Silent, harmless. | no | COSMETIC |

### 9a. Recommended commit shape

1. **Pre-flight (no move):** fix `exclude_patterns` → `['_build', '_generated/*', 'Thumbs.db', '.DS_Store']`; land it; confirm `-W` still green. This de-risks #2 by making `_generated` handling deliberate instead of accidental.
2. **Convert all 132 REL `:doc:` in `docs/theory/` → ABS `/theory/…`.** Bit-identical rendering, zero move yet, `-W` green. Now the move only has to touch *targets*, not sources.
3. **Parameterise the capability-matrix generator** (`DOCS_THEORY_DIR` per package + `rglob` wipe) and move `tests/derivations/test_capability_matrices.py:29` in lockstep. Green before the move.
4. **`git mv` hub-first** (hub + its toctree children together — §2), fixing: toctree entries only where hub/child separate; the 8 `../_generated/` includes; the 2 `.inc.rst` includes; all ABS `:doc:` targets.
5. **Same commit:** the 97 actionable raw paths + the ~125 inert `.py` roles + `explorer/AGENT.md:162–166`.
6. **Leave** `.claude/plans/`, `.claude/agent-memory/`, `.claude/scratch/` (663 rows) — archaeology.
7. **Gate:** `sphinx-build -W docs docs/_build/html` exit 0 **AND** `pytest tests/derivations/test_capability_matrices.py tests/test_vv_harness_audit.py` green **AND** a re-run of §4's stale-path detector showing 52 → (52 − 7 fixed) and no NEW stale entries.

> **A `-W` clean build is necessary but nowhere near sufficient here.** It cannot
> see: 837 raw path strings, 156 `:file:` roles, ~125 inert `.py` roles, or 4
> Python `Path` constants. The stale-path detector (script preserved at
> `/private/tmp/claude-501/-Users-rodrigo-git-nuclear-ORPHEUS/e9c6341d-23b8-475b-b00d-7258a15b0a4a/scratchpad/rawpaths3.py`)
> is the only gate that covers the silent class — consider promoting it to
> `tools/` as a permanent guard, since §4b shows this failure recurs on **every**
> doc rename this project has ever done.

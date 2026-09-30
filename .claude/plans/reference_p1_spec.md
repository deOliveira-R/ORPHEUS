# #405 step 2, phase P1: the specification and the question (verification specification)

Workflow W3 (surgical carve), phase P1 of `.claude/plans/reference_cache.md` ("P1, the carve order", with the rulings of "P1 opened", "P1, second exchange" and "P1, third exchange"). Written by the test-architect on 2026-09-25 against `refactor/reference-specification` at `26a4f95c` (tracked `orpheus/` and `tests/` clean; `git status --porcelain -- orpheus tests` empty). The main agent writes the code, step by step with the user; this file is written before any code and specifies, per step, the gates that land with it.

Status: a specification only. No production or test code has been written. Every probe named below is a re-runnable file in `scratch/reference_architecture/p1probe/` (untracked), run from the repository root with `.venv/bin/python -O`. The canonical test invocation is `.venv/bin/python -O -m pytest`. Markers: `[M]` measured (the probe or command is named), `[R]` reasoned, `[HYPOTHESIS]`, `[REFUTED 2026-09-25]`.

## 0. What changed from the plan (read first)

Rulings received during this specification (the orchestrator, 2026-09-25, reported to the user):

1. **The carve order is swapped at steps 6 and 7, and step 6 is renamed.** Step 6 is now *the phase-space functions* (`Symbolic`, `RegionwiseConstant`, the anisotropy predicate, `srepr` storage), named by structure, not by role. Step 7 is the question, whose adjoint `FixedSource` holds a detector of that function type. The specification's fields name the roles (`source`, `detector`). Reason: the adjoint fixed-source question needs a detector, which is a mesh-free function on phase space, the same kind of object as the source (§1.6, §1.7).
2. **The equal-width edge rule is R2, one body.** A spacing rule is affine in the measure coordinate `T(r)`: equal width uses `T(r) = r` in every coordinate system, equal volume uses `T(r) = r`, `r²`, `r³`. Edges are `T⁻¹(T(r_k) + f_j · (T(r_{k+1}) − T(r_k)))` with `f = linspace(0, 1, n+1)`, and both interval ends are the geometry's breakpoints exactly. On a slab the two rules are one code path, so #495 is fixed at its root. The rival R1 (`np.linspace(r_k, r_{k+1}, n+1)`) is refused because it moves the SN regression corpus (§3, §7).
3. **Hollow geometries are admitted.** The geometry is the interval of breakpoints `(r_0, …, r_R)` with `r_0 ≥ 0` for cylinder and sphere (any `r_0` for a slab). Its boundary is the topological boundary of that region in its coordinate system, so the number of boundary laws is derived: 2 for a slab; 1 for a cylinder or sphere with `r_0 = 0` (the centre is an interior point); 2 otherwise. A law at the centre is refused, and a hollow curvilinear geometry without an inner law is refused. `origin=` retires.
4. **`None` retires as a boundary declaration at step 3.** The geometry requires one law per boundary point. Each of today's `None` sites migrates to the law its consumer applied, and a second capture proves it (§3.2). The four consumer defaults (§0 item 8) are a weld and retire.
5. **#511 (filed by the orchestrator): SN silently ignores a declared inner law on a hollow curvilinear mesh.** At step 3 SN refuses any inner law other than `reflective`, with a `SCOPE-BOUNDARY[guard]` naming the missing machinery (an inner-surface trace on `RadialAxisMesh`); the curvilinear SN build is deferred by the user.
5b. **#513 (filed by the orchestrator): CP, MoC and MC read only some of a geometry's declared laws.** Ruling: every method that reads only some of the declared laws refuses a geometry whose other laws it would drop, a `SCOPE-BOUNDARY[guard]` recorded on #513 (CP and MC; the CP slab's mirror asymmetry is a comment there), #514 (MoC) and #511 (SN's inner surface); what a method realises on the unread face decides which law it admits there, and a site migrates to that law so it stays bit-identical (S3.12 to S3.14). #512 (a CP slab with a vacuum right law gives k = 1.9137, above k∞ = 1.875) is out of P1's scope.

Premise corrections measured while writing this specification:

6. **"Equal width: equal widths" is false bitwise for every realisable edge formula.** `[M]` `partition_laws.py`: the differences of `np.linspace` edges are not all equal in 772 of 852 (interval, n) cases, and the equal-volume sum of measures differs from the interval's closed-form measure by up to 6 ULP. The laws of §1.3 are therefore stated on the stored MEASURES (exact per rule) and on the realised widths within a derived ULP bound, never as bitwise equality of edge differences.
7. **"Every cell width ≤ h" is false for the natural count rule on realised edges.** `[M]` `maxwidth.py`: with `n = ceil(L/h)`, some realised width exceeds `h` in 14 of 114 399 (L, h) pairs, among them `L = 1, h = 0.1` (a width `0.10000000000000009`). The count law of `CellsByMaxWidth` is stated on the rule's nominal width `fl(L/n)`, and the realised width is bounded by `h + 2·ulp(r_{k+1})` (measured worst 1.38 ULP, `misc_laws.out`).
8. **`None` means four different laws today** `[M]` (grep, `orpheus/`): reflective in `SNProblem` and `resolve_boundary_conditions` (`transport/method.py`, whose docstring's "uniform across methods" is false); the entry's `boundary_condition` argument (default `"vacuum"`) in `solve_sn_fixed_source`, `solve_sn_adjoint_fixed_source` and `solve_sn_multiplying_source` when both mesh laws are `None` (`sn/solver.py:120-140`); `BC("white")` in CP (`cp/solver.py:228`); reflective in MoC (`moc/geometry.py:308`); `BC("periodic")` in MC (`mc/solver.py:184`).
9. **The Mixture "dtype" identity leg of the brief is vacuous.** `[REFUTED 2026-09-25]` "a change of the arrays' dtype changes the digest": `Mixture.__post_init__` stores every dense field as float64 and every CSR block with int32 indices, so a float32 or int64-index input constructs a mixture EQUAL to the float64 one `[M]` (`identity.py`: `Mixture == with int64 CSR indices: True`; float32 `SigT` stored as float64). The explorer's caveat on CSR index dtype is refuted by the same line. The fact that survives: dtype is canonicalised at construction, so the digest reads canonical bytes and needs no dtype leg; a leg asserting that a float32-ROUNDED value (a value change) moves the digest is the honest replacement, and it is vacuous on `xs_library` mixture A (its `SigT` is float32-representable; `[M]` equal after rounding), so it runs on a manufactured mixture.
10. **Signed zero is two conventions today** `[M]` (`identity.py`): `Mixture` equality is by bytes, so `SigL = −0.0` and `SigL = +0.0` are UNEQUAL mixtures, while `BC` compares its params with float `==`, so `BC("p", {"a": −0.0}) == BC("p", {"a": 0.0})`. The shared digest encoder needs one convention (NEEDS 1).
11. **`Mesh1D == Mesh1D` raises today** `[M]`: `ValueError: The truth value of an array ... is ambiguous` (even at N = 1), and `hash(mesh)` raises `TypeError`. The step-3 law "two rules giving the same cells give equal meshes" is new, and today's `==` is its first red.
12. **The geometry package imports the mesh module at two sites that the move must cut** `[M]` (`move_census.py`, AST): `orpheus/geometry/__init__.py:29` (`BC, Mesh1D, Mesh2D, RegionMesh`) and `orpheus/geometry/structured_geometry.py:115` (`BC`). `orpheus/geometry/factories.py` also constructs `Mesh1D` and `Mesh2D` (`pwr_pin_2d`), so the whole `factories` module moves into `orpheus/mesh/`, not only `_subdivide_zone`; otherwise the new layer row is violated.
13. **The layer linter cannot see a relative import** `[M]` (`misc_laws.py`): `_imports_of` returns `("transport", False)` for `from ..transport import x` and drops it at the `startswith("orpheus.")` filter. Today 0 relative imports cross a top-level package in `orpheus/` `[M]` (`move_census.py`), so resolving them changes no current reading; the new `mesh` package would otherwise be able to import `transport` unseen.
14. **Hollow curvilinear sites exist in the tree** `[M]`: 13 of the 263 direct sites whose edges evaluate statically start at `r_0 ∈ [0.01, 0.1]` on a cylinder or sphere, 5 of them pass `bc_left = BC("reflective")`; 1 slab starts at 1.0. The 0.01 cylinders are documented as avoiding the pole (`test_affine_carve_baseline.py:143`), so their honest geometry is hollow with an inner reflective law, which is also what SN computes (#511).

## 1. The ladder, per step

Each step lands as one commit, green under `.venv/bin/python -O -m pytest` on what it touches. For each gate: its claim kind (THEOREM: a law true for every admissible input; REFERENCE: a structurally independent route; RECORD: what the code printed on a given day), its first red (the input existing in the tree when it lands that the gate rejects, `plan-authoring` §6c; for a new type, its defining refusals), its mutation witness, its level, and the rungs it rests on. Every new, improved or reused gate carries `@pytest.mark.rests_on(...)` with the node ids named here. Unless stated otherwise a gate is `foundation` (a software or mathematical invariant with no theory-page label, so no `verifies`).

A refusal row is a `pytest.raises(<typed error>, match=<fragment>)` whose fragment is keyed to the argument that triggers it; the fragments of one constructor's refusals are disjoint, and one row asserts their disjointness (test-architect lessons §1, "the message IS the gate").

### 1.1 Step 1: the mesh module, a pure move

What moves: `Mesh1D`, `Mesh2D`, `RegionMesh` (until step 3), the whole of `geometry/factories.py` (`_subdivide_zone` and `pwr_pin_2d`, §0 item 12), and `transport/mesh/axis.py` into `orpheus/mesh/`; `BC` into `orpheus/geometry/boundary/`. No re-export shim.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S1.1 | Linter rows, in `tests/gates/test_layer_imports.py`: a new layer `MESH_PACKAGES = {"mesh"}`; `FORBIDDEN_EDGES["mesh"] = L2 ∪ L3`; `"mesh"` added to the forbidden set of `geometry`, `data` and `numerics`; `"mesh"` joins the set `test_input_layer_imports_numerics_only_by_submodule` ranges over. The existing parametrised `test_no_forbidden_imports` then covers every file of the new package | THEOREM (the layer contract) | the two geometry-to-mesh import sites of §0 item 12, re-spelled as a shim `from orpheus.mesh... import Mesh1D` in `geometry/__init__.py`: the linter must red on that file. `[R]` (the file does not exist before the move; the check is S1.2's leg b, which is measurable today) | delete `"mesh"` from `FORBIDDEN_EDGES["geometry"]` → S1.2 leg b reds | foundation; none |
| S1.2 | Linter self-test (new): the linter exposes a pure `_check_source(rel_path, source) -> list[str]`, and six legs feed it synthetic sources: (a) `mesh/x.py` importing `orpheus.transport.fields` → a violation; (b) `geometry/__init__.py` importing `orpheus.mesh.mesh1d` → a violation; (c) `data/x.py` importing `orpheus.mesh` → a violation; (d) `numerics/x.py` importing `orpheus.mesh` → a violation; (e) `mesh/x.py` with `from ..transport import fields` → a violation; (f) positive legs, no violation: `mesh/x.py` importing `orpheus.geometry.boundary` and `orpheus.numerics.measure`, `transport/x.py` importing `orpheus.mesh`, `derivations/x.py` importing `orpheus.mesh` | THEOREM | leg (e) is red on today's linter `[M]` (`misc_laws.py`: the relative import is dropped); legs (a)-(d) are red until S1.1's rows exist | one arm per leg: drop the relative-import resolution → (e) reds; drop each row of S1.1 → its leg reds; a forbid-everything linter → leg (f) reds | foundation; rests on nothing |
| S1.3 | Cold-interpreter imports: `orpheus.mesh`, each moved module, `orpheus.geometry.boundary` join the parametrisation of `test_entry_point_imports_in_a_fresh_interpreter` | THEOREM (no import cycle) | today `import orpheus.mesh` fails (the package does not exist); the load-bearing red is a cycle the move could introduce (`BC` now lives in `geometry.boundary`, which `mesh` imports, and whose lazy imports were written against the old layout, `mesh.py:73`) | re-introduce a module-level `from orpheus.mesh import Mesh1D` inside `geometry/boundary/__init__.py` → the entry reds in a fresh interpreter | foundation; none |
| S1.4 | No shim (new, `tests/gates/mesh/test_module_layout.py`): `importlib.util.find_spec` is `None` for `orpheus.geometry.mesh`, `orpheus.geometry.factories` and `orpheus.transport.mesh.axis`; `orpheus.geometry` has none of the attributes `Mesh1D`, `Mesh2D`, `RegionMesh`, `pwr_pin_2d`; `orpheus.transport.mesh` has none of `AxisMesh`, `RadialAxisMesh`, `Axis1D`, `AxisCoord`; `BC.__module__` starts with `orpheus.geometry.boundary`, and if the `orpheus.geometry` package root still exports `BC` (a name of the same package, not a shim across the move; the main agent's choice), `orpheus.geometry.BC is orpheus.geometry.boundary.BC` (one definition) | THEOREM (the layout) | every leg is red today `[M]` by construction (the three modules exist; `move_census.out`) | add a shim module `orpheus/geometry/mesh.py` doing `from orpheus.mesh.mesh1d import *` → the first leg reds | foundation; none |
| S1.5 | Behaviour-neutral proof, not a pytest gate, recorded in the commit body: (a) `git diff -M --stat` shows each moved file as a rename, and the added/removed CODE-line list (`git diff -U0 -M orpheus/ \| <strip comments and docstrings> \| grep -E '=\|def \|raise \|return '`) contains only import statements; (b) the SN regression drift SET is unchanged: `-W error::tests.gates.sn.regression._regression_assert.DriftWarning` on `tests/gates/sn/regression/` gives the same failed-node set and per-case ULP counts before and after (about 1.6 s; the delta, never the absolute count, lessons §3); (c) the full canonical suite, serial, gives the same passed / skipped / xfailed / xpassed counts as the pre-move baseline (§5); (d) `sphinx-build -W -n` with the `-n` warning set diffed before and after over the pages that change (`docs/api/geometry.rst` holds `autoclass:: orpheus.geometry.mesh.BC`, `automodule:: orpheus.geometry.mesh`, `automodule:: orpheus.geometry.factories`); (e) `dead_references` (§4) | RECORD | the pre-move baseline counts, captured at `26a4f95c` before the first edit | (b): a one-ULP perturbation of `_subdivide_zone`'s volume (a monkeypatch in the regression conftest) → the drift set grows | not a gate |

Import census the move must re-point `[M]` (`move_census.py`, AST over `git ls-files '*.py'`): statements importing a moving MODULE, orpheus 25 (`geometry.mesh` 12, `transport.mesh.axis` 11, `geometry.factories` 2), tests 68 (47 + 21), examples 1; statements importing a moved NAME from the `orpheus.geometry` package root, orpheus 16, tests 279 (252 files), examples 7. String constants naming a moving path: 46, all docstrings (43 in `orpheus/`, 3 in `tests/`), which is `dead_references`'s surface. Documentation: 104 lines in 20 `.rst` files under `docs/` (build tree excluded), and 47 lines in 22 files under `.claude/` `[M]` (`git grep -c`). The full suite's collection is the completeness check for the code imports: a stale import is a collection error, and `pytest --collect-only -q tests/gates` must report the same item count before and after with 0 errors (8.75 s, P0's measurement).

### 1.2 Step 2: the geometry value

`StructuredGeometry` carries `coord: CoordSystem` (the `"SLB"/"CYL"/"SPH"` tag retires), its breakpoints `(r_0, …, r_R)`, one material id per interval, and one boundary law per boundary point (§0 item 3). `from_thicknesses(coord, thicknesses, mat_ids, boundaries, r_0=0.0)` serves the registries that speak in thicknesses.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S2.1 | Breakpoint laws, refusals keyed and disjoint: fewer than 2 breakpoints; not strictly increasing (one leg equal, one leg decreasing); a non-finite entry (`nan`, `inf`); `r_0 < 0` on a cylinder or sphere; a material-id count other than `R`; a non-int material id. Positive legs: a slab with `r_0 = −1.0`; a sphere with `r_0 = 0.1`; adjacent intervals with the same material (the 27 irregular sites need it, census §(2)) | THEOREM (defining refusals of the value) | defining refusals of a new constructor signature; the re-posed rows are `TestRegion::test_zero_thickness_rejected` and `test_negative_thickness_rejected` (a thickness ≤ 0 is a non-increasing breakpoint pair) | drop each clause → exactly its leg reds (a per-arm table, 7 arms) | foundation; none |
| S2.2 | The boundary is derived: the admitted law count is 2 for a slab, 1 for a curvilinear geometry with `r_0 = 0`, 2 for `r_0 > 0`, over the six cells coord × {`r_0 = 0`, `r_0 > 0`}; refusals keyed and disjoint: a law at the centre (curvilinear, `r_0 = 0`, two laws), a hollow curvilinear geometry with one law, a slab with one law, a `None` in the laws (`TypeError`, §0 item 4). The laws are indexed by boundary point (inner, outer), never by a coordinate-specific name | THEOREM | hollow geometries are inexpressible today (13 hollow direct sites, §0 item 14); `StructuredGeometry("SPH", …, bcs=(BC.vacuum, BC.vacuum))` is refused today for a different reason (it requires 1), so the centre-law leg is a DISCRIMINATION row: on that input the new message names the centre and the old "requires 1 BC" fragment is absent | replace the derived count by the old table `{SLB: 2, CYL: 1, SPH: 1}` → the two hollow legs red | foundation; S2.1 |
| S2.3 | `from_thicknesses` is the left fold of today's accumulation, bitwise: `breakpoints == tuple(itertools.accumulate(thicknesses, initial=r_0))`, so the registries (`wigner_seitz_pin_cell`, `pwr_slab_half_cell`, `la13511`) keep their bits; `np.cumsum` equals that fold in 20 000 of 20 000 random cases `[M]` (`misc_laws.py`) | THEOREM (a float identity, both sides the same sequential sum) | defining (new constructor) | a pairwise or `math.fsum`-based cumulative sum → reds on a 5-term input (a fold with a different association) | foundation; S2.1 |
| S2.4 | Breakpoints are stored, never re-derived: `StructuredGeometry(breakpoints=E, …).breakpoints == E` bitwise for `E = (0, .17, .45, .62, 1.0)` | THEOREM | `[M]` (`misc_laws.py`): today's route re-adds thicknesses and gives `(0.0, 0.17, 0.45000000000000007, 0.6200000000000001, 1.0)`, the census's 1-ULP site `test_g_adjoint_reciprocity.py:244` | store `accumulate(diff(E))` → red | foundation; S2.1 |
| S2.5 | The kind tag retires: a `str` coordinate is refused with `TypeError` keyed on the argument; `coord` is a `CoordSystem` member; `_GEOMETRY_TO_COORD` and `_GEOMETRY_TO_N_ENDPOINTS` no longer exist (an `AttributeError` leg on the module) | THEOREM | today `StructuredGeometry("SLB", …)` constructs `[M]` (every census row) | accept `"SLB"` through a map → the str leg reds | foundation; none |

The 13 closed kind vocabularies in `orpheus/` (explorer report §3) are not migrated here: derivations' 287 kind strings move with their families in P4 (the plan's ruling). Step 2 retires only V2 (`_GEOMETRY_TO_COORD`) and the test-side V13 (`tests/gates/cp/test_properties.py:47`).

### 1.3 Step 3: the discretisation and the one constructor

Types: `Partition` (per interval, the cell edges and the cell measures), the spacing rules `EqualWidth` and `EqualVolume` (the names are the main agent's; §0 item 2 fixes their bodies), `CellsByCount(counts, spacing)` with `.uniform_width(n)` and `.uniform_volume(n)`, `CellsByMaxWidth(widths, spacing)`, `2 * d`, and `Mesh1D(geometry, discretization)`. `precomputed_volumes`, `RegionMesh`, `from_geometry`, `bc_left`/`bc_right` and `None` as a boundary declaration retire. `Mesh2D` keeps its constructor.

The measure `m(a, b)` of an interval `[a, b]` is `b − a` (slab), `π(b² − a²)` (cylinder), `(4/3)π(b³ − a³)` (sphere).

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S3.1 | Partition nests in the geometry: in every interval the first edge is `r_k` and the last is `r_{k+1}`, BITWISE; edges strictly increasing; measures > 0; one cell count ≥ 1 per interval. An explicit `Partition` is checked against the geometry: refusals keyed and disjoint for an end edge 1 ULP off its breakpoint (`np.nextafter`), an interval-count mismatch, a non-increasing interior edge | THEOREM | defining refusals (new type). `[M]` (`partition_laws.py` (a)): today's `_subdivide_zone` already lands both ends exactly in 852 of 852 cases per coordinate, so the law costs the equal-volume sites nothing | drop the end pinning in the R2 body → the equal-width rows red at the non-dyadic lengths | foundation; S2.1 |
| S3.2 | The sum law: per interval, `\|Σ measures − m(r_k, r_{k+1})\| ≤ (2 + ⌈log₂ n⌉) · ulp(m)`, for both rules and all three coordinates, over `n ∈ {1, 2, 3, 5, 8, 13, 64, 100, 1000, 4096}` and a seeded population of intervals, hollow ones included | THEOREM with a derived tolerance: a stored measure `fl(m/n)` carries ½ ULP of `m/n`, `n` of them carry at most about 1 ULP of `m`, and the pairwise sum adds at most `⌈log₂ n⌉` half-ULPs; measured worst 6 ULP over 3 × 852 cases and 4 ULP over 3 × 4000 `[M]` (`partition_laws.out` (d), `misc_laws.out`) | none in today's rules (a defining law) | a measure that drops the interval's inner radius, `V(0, r_{k+1})/n`, → red on every interval beyond the first (the second leg of today's ERR-020 row, which this law generalises) | foundation; S3.1 |
| S3.3 | Equal volume (ERR-020 preserved): per interval the stored measures are bitwise identical and each is `fl(m/n)`. REUSE, re-spelled onto `CellsByCount.uniform_volume`: `test_equal_volume_cylindrical_invariant`, `test_equal_volume_spherical_invariant`, `test_equal_volume_multi_region_invariant[SLB\|CYL\|SPH]`, `test_multi_region_cylinder_equal_volume` (all `catches("ERR-020")`, carried with them, `retirement-audit` C.13) and `test_equal_volume_edges_bound_the_volumes` (the radius-formula witness) | THEOREM | reused; their first red is ERR-020 itself | measures from `compute_volumes_1d(coord, edges)` (the ERR-020 defect) → the catchers red on all three coordinates: on the three-region fixture (radii 0, 0.5, 1.5, 2.0; cells 5, 7, 11) the re-derived volumes take 3/3/2 distinct values per region on a slab, 4/3/3 on a cylinder and 5/6/5 on a sphere `[M]` (one-line probe on `_subdivide_zone`, 2026-09-25) | foundation; S3.1, S3.2 |
| S3.4 | Equal width: the edges are the R2 body, `r_k + f_j (r_{k+1} − r_k)`, bitwise; every realised width is within `2 · ulp(r_{k+1})` of `fl((r_{k+1} − r_k)/n)`; the stored measure is `fl((r_{k+1} − r_k)/n)` on a slab and the shell measure between the realised edges on a cylinder or sphere | THEOREM | `[M]` (`issue495.py`): today's `"uniform"` stores `diff(edges)` on a slab, 3, 5, 3 and 5 distinct volumes at n = 5, 7, 11, 9 on `[0, 3]` | store `diff(edges)` as the slab measure (today's rule) → the measure leg reds at n = 5 | foundation; S3.1 |
| S3.5 | **The #495 law** (§2): on a slab, `uniform_width(n)` and `uniform_volume(n)` give EQUAL partitions (edges and measures bitwise) and equal meshes, for n ∈ {1, …, 64} ∪ {100, 127, 255, 1000} on `[0, 3]`, `[0, 2.872]`, `[1.1, 1.8]` and a 3-interval slab; on a cylinder and a sphere they differ at every n ≥ 2 (the negative leg, which proves the comparison can fail) | THEOREM | `[M]` `issue495.out` (§2) | give `EqualWidth` its own Cartesian body `np.linspace(r_k, r_{k+1}, n+1)` (R1) → the edge leg reds in 576 of 852 cases `[M]` (`partition_laws.out` (b)) | foundation; S3.3, S3.4 |
| S3.6 | Refinement: `(2 * d).partition(g)` refines `d.partition(g)`: per interval the fine edges at even indices equal the coarse edges BITWISE, for `CellsByCount` under both spacings and for `CellsByMaxWidth`; `2 * d` doubles the counts `d` realises on `g` | THEOREM | `[M]` (`partition_laws.out` (e)): 0 of 840 cases per (coordinate, rule) violate it for today's formulas; exact by argument `[R]`: `fl(1/(2n)) = fl(1/n)/2`, so `2k · fl(1/(2n)) = k · fl(1/n)`, and the affine map then sees the same fraction. For `CellsByMaxWidth`, `2 * d` spelled as "halve the width" is `[REFUTED 2026-09-25]` as a refinement: `ceil(2L/h)` is odd whenever the fractional part of `L/h` exceeds ½, and an odd count does not nest | spell `2 * CellsByMaxWidth(h)` as `CellsByMaxWidth(h/2)` → the `CellsByMaxWidth` leg reds at `L = 1, h = 0.3` (counts 4 and 7) | foundation; S3.1 |
| S3.7 | `CellsByMaxWidth` count law: per interval, `n = min{n : w(n) ≤ h}` where `w(n)` is the rule's nominal widest cell: `fl(L/n)` for equal width, the first cell's width from the equal-volume body for equal volume (the widest, and monotone in n: `[M]` 0 exceptions over 2 × 1995 cases, `maxwidth.out`); minimality leg `w(n − 1) > h`; realised-width leg `max width ≤ h + 2 · ulp(r_{k+1})` (measured worst 1.38 ULP over 132 220 pairs, `misc_laws.out`); discriminating rows `(L, h) = (1.0, fl(1/3)) → 3` and `(3.0, fl(1/3)) → 9`. Refusals: `h ≤ 0`, a width count other than the interval count | THEOREM | defining (new type). The naive spelling is the witness below | `n = ceil(L/h)` checked on realised widths → the realised-width leg reds at `L = 1, h = 0.1` (width `0.10000000000000009`, §0 item 7); the exact-rational count (`ceil(Fraction(L)/Fraction(h))`) → the discriminating rows red (it returns 4 and 10, `misc_laws.out`) | foundation; S3.4 |
| S3.8 | The one constructor: `Mesh1D(geometry, discretization)` stores the geometry and the partition, never the rule. (a) `Mesh1D(g, CellsByCount.uniform_width(n)) == Mesh1D(g, <the Partition it produced>)` and the hashes agree; (b) `mesh.mat_ids` is each interval's material broadcast over its cells; (c) `mesh.volumes` IS the partition's stored measures (`np.shares_memory` or `is`), never recomputed from the edges; (d) `mesh.boundaries is geometry.boundaries`; (e) `widths`, `centers`, `areas` as today. Refusals: an interval-count mismatch between the discretisation and the geometry; `CellsByCount(counts)` without a spacing (`TypeError`: no default spacing, the ruling of 2026-09-25); a count ≤ 0 or non-int (the re-posed `TestRegionMesh` rows); `CellsByCount.uniform` does not exist (`AttributeError`) | THEOREM | `[M]` today `Mesh1D == Mesh1D` raises `ValueError` (§0 item 11); today `RegionMesh(n).method == "equal-volume"` is silently applied (the default the ruling removes) | recompute `volumes` from the edges in `Mesh1D` (the pre-ERR-020 path) → leg (c) and S3.3 red; store the rule instead of the partition → leg (a) reds | foundation; S3.1, S2.2 |
| S3.9 | SN refuses an inner law other than `reflective` on a hollow curvilinear geometry: `SCOPE-BOUNDARY[guard]` naming the machinery (an inner-surface trace on `RadialAxisMesh`) and #511; `pytest.raises` keyed on the fragment naming #511, over the inner laws `vacuum`, `white` and `albedo(0.5)` on a cylinder and a sphere; positive leg: an inner `BC.reflective` solves and its flux equals the pre-carve capture of the same fixture bitwise | THEOREM (a declared boundary of the realised machinery) | `[M]` `hollow_inner_law.out`: today an inner `vacuum` and an inner `reflective` give bit-identical fluxes on a cylinder and a sphere with `r_0 = 0.5` (the declared law is dropped); the slab control moves | remove the guard → the vacuum legs solve silently → red | foundation; S2.2 |
| S3.10 | `None` is not a boundary declaration: the geometry refuses `None` (S2.2's leg); `Mesh1D` has no `bc_left`, `bc_right` or `precomputed_volumes` field (`AttributeError` legs); `boundary_condition` is not a parameter of `solve_sn_fixed_source`, `solve_sn_adjoint_fixed_source`, `solve_sn_multiplying_source` or `solve_sn` (`inspect.signature`) | THEOREM (the retirement's witness, `retirement-audit` D.18) | every leg is red today `[M]` (the fields and parameters exist; `sn/solver.py:157, 2829, 3246, 3542`) | re-add a `boundary_condition=None` pass-through → the signature leg reds | foundation; none |
| S3.11 | The re-posed constructor refusals of `tests/gates/geometry/test_geometry.py` (the census's class G, 5 sites at 218, 223, 228, 233, 382). Non-monotone and equal edges → S3.1's explicit-partition refusal (RE-POSED); too few edges → S2.1's breakpoint-count refusal (RE-POSED); a wrong `mat_ids` length → DEAD (the material map is per interval and a per-cell array can no longer be passed: delete, never repair by passing the new argument); an invalid `bc_left` type → S2.2 and the existing `test_bcs_must_be_a_tag_or_a_typed_law` (DEMOTED onto the geometry) | each as its new home | each is its own re-pose | as the home gate | foundation |
| S3.12 | CP refuses a slab whose two laws differ (`SCOPE-BOUNDARY[guard]`, #513): the only CP slab declaration that means what CP computes is left law = right law, because CP reads only `bc_right` (`cp/solver.py:228`); refusal keyed on #513 for (left, right) ∈ {(vacuum, white), (white, vacuum)}; positive leg (white, white) solves. Every CP slab site migrates with left := right, bit-identical because the left law is not read | THEOREM (a declared boundary of the realised machinery) | `[M]` (`cp_slab_left_law.out`): left `white`, `vacuum` and `None` give k = 1.8749980808246423 to all 16 digits (dropped); the right law is read (control). Not specified: which law CP realises on the left face. `[M]` (`cp_mirror.out`): a fuel\|moderator slab with BOTH faces white is not mirror-symmetric (k 1.2183664822172493 against 1.218306071627994, relative 5.0e-5, `keff_tol = 1e-13`), and neither k equals the mirrored double slab (1.2188137), so the left face is neither the right law nor a mirror: a CP defect reported to the orchestrator, out of P1's scope | remove the guard → the (vacuum, white) leg solves → red | foundation; S2.2 |
| S3.13 | MoC refuses any geometry other than a solid cylinder (`CYLINDRICAL`, `r_0 = 0`), whose single outer law it reads (registry: `reflective`), `SCOPE-BOUNDARY[guard]` on #514: its `Mesh1D` is a Wigner-Seitz proxy for a 2-D square pin cell (`moc/geometry.py:253`) | THEOREM (a declared boundary) | `[M]` (`partial_law_readers.out`): a CARTESIAN `Mesh1D` with the cylinder's edges constructs and gives exactly the cylinder's k (0.7197934814031874: the coordinate is not read); a hollow cylinder `r_0 = 0.1` constructs and gives 0.71516, because the region areas read `edges[0]` while the tracks treat region 0 as a disk (`moc/geometry.py:272-287`), an inconsistent realisation | remove the guard → both legs construct → red | foundation; S2.2 |
| S3.14 | MC refuses a geometry unless every law is `periodic` (its registry's only kind), `SCOPE-BOUNDARY[guard]` on #513 | THEOREM (a declared boundary) | `[M]` (`partial_law_readers.out`): `MCMesh` constructs a slab with a `vacuum` left law and a `periodic` right law (the left law is dropped) | remove the guard → the vacuum-left leg constructs → red | foundation; S2.2 |

The equal-volume `2 * d` and the equal-width `2 * d` legs of S3.6 use one population: intervals `(0, 1)`, `(0, 3)`, `(0, 2.872)`, `(0.9, 1.1)`, `(1.1, 1.8)`, `(0, 0.41)`, `(0.01, 2)`, `(2, 7)`, and n up to 500, the probe's own.

### 1.3a Step 3 split into 3a, 3b and 3c (added 2026-09-29)

⛔ The 3a column of the placement table and the round-1 battery table below describe the retired `Partition` API; read "Round 2" at the end of this section first.

Written by the test-architect after the user's rulings of 2026-09-29 (`.claude/plans/reference_cache.md`, "P1 step 3 opened": ruling 1 splits step 3 into three sub-commits, ruling 3 adds the named-face geometry constructors, ruling 4 moves the discretisation digest to step 5). The step-3a production was read in the main tree (uncommitted, on `refactor/reference-specification` at `08fe3e72`): `orpheus/mesh/partition.py`, `CoordSystem.measure_exponent`/`measure_constant`/`interval_measure` in `orpheus/geometry/coord.py`, and the constructors `slab`, `cylinder`, `sphere`, `uniform_boundary` and `from_homogeneous` on `StructuredGeometry`.

**Premises measured or corrected while writing the 3a gates.**

1. `[M]` The strict xfail §2 placed at step 1 never landed: `git grep 495` over `tests/gates/mesh/`, `test_structured_geometry.py` and `test_geometry.py` returns nothing. The #495 law is therefore born green in 3a, on partitions, with no mark to remove; its mesh-level leg lands in 3b.
2. `[M]` (scratchpad probe, 21 035 random (interval, n) cases per coordinate) The R2 body without end pinning lands its last edge off the breakpoint in 168 cylinder and 21 sphere cases (0 on a slab). `_subdivide_zone` does not pin, so a `from_geometry` curvilinear mesh can in principle differ from its partition at an interval's last edge; capture 1 found 0 such sites in the suite (§3.4, all `from_geometry` rows `ev`), so `[R]` 3b's capture should show none.
3. `[M]` `k * rule` nests bit for bit only when `k` is a power of two: fine[::k] differs from coarse in 0 of 7176 cases for k = 2, 4, 8, in 3909 for k = 3 and 5300 for k = 5 (worst 2 ULP). Ruled 2026-09-29 (the orchestrator, option (a)): a factor that is not a power of two is refused, ValueError "a refinement factor is a power of two".
4. `[M]` The nominal widest cell of the equal-volume body on a slab must be `fl(L/n)`, not the first realised cell: on `[0, 3]` with `h = 0.3` the realised first cell `fl(fl(1/10)·3) = 0.30000000000000004` exceeds `h` (1 of 54 144 slab pairs), which would give 11 cells for equal volume and 10 for equal width. The landed `Spacing.widest_cell` takes the `p = 1` branch for both rules on a slab, so the one-body law holds; the row `TestIssue495::test_max_width_is_one_body_on_a_slab[0-3-h0.3]` pins it.
5. `[REMEDIED 2026-09-29]` `uniform_boundary` decided hollowness a second time, on the RAW first breakpoint (`== 0.0`), while the geometry decides it on the PARSED float. Witness: `Fraction(1, 10**400)` (a `Real` whose float is 0.0) was refused. The orchestrator moved both decisions onto `_is_hollow`/`_boundary_points` over parsed breakpoints; the witness row is `TestUniformBoundary::test_hollowness_is_the_geometry_s_decision`.
6. `[M]` The spelling of the interval measure (`c·(b^d − a^d)` against `c·b^d − c·a^d`) is invisible on most intervals: the two agree on 16 of 16 intervals of the law population on a sphere and 15 of 16 on a cylinder. `(0.3, 0.7)` is added to the population as the sphere's witness, so the bitwise `fl(m/n)` pin sees a re-association.
7. `[M]` The spec's "3 distinct volumes at n = 5 on `[0, 3]`" is a property of today's `np.linspace` edges; the R2 edges' differences take 2 distinct values there. The discriminating row asserts "not all equal".

**Where every row of §1.3 lands.** Gate files: `tests/gates/mesh/test_partition.py` (195 rows) and `tests/gates/geometry/test_named_face_constructors.py` (58 rows); `[M]` 253 passed under `.venv/bin/python -O -m pytest` in 1.7 s, and `npx pyright` reports 0 errors on both files. Every row is `foundation` and carries `rests_on` edges into S2.1/S2.2 (`test_structured_geometry.py`) and within the two files (21 distinct targets, all collected node ids, `[M]`).

| row | 3a (landed with this spec section) | 3b | 3c / 5 |
|---|---|---|---|
| S3.1 nesting | `TestPartitionValue` (12 keyed refusals and their disjointness row, equality bitwise, copy and read-only, the flat views stored once, `partition(g) is p`, the unhashable RECORD row), `TestNesting` (both spacings × three coordinates × 17 intervals × 10 counts, plus a three-interval body; the pinning's own witness row) | — | the unhashable row reds at step 5 by design and is re-posed as the eq/hash contract |
| S3.2 sum law | `TestSumLaw` (6 rows, 170 cases each) | — | — |
| S3.3 equal volume | `TestEqualVolume`: `fl(m/n)` broadcast, per interval, bitwise (`catches("ERR-020")`, earned: see the battery), the R2 edges, and the keystone against `_subdivide_zone` | the keystone row RETIRES with `_subdivide_zone`; the six `from_geometry` rows of `test_structured_geometry.py` (`test_multi_region_cylinder_equal_volume`, `test_equal_volume_cylindrical_invariant`, `test_equal_volume_spherical_invariant`, `test_equal_volume_multi_region_invariant[×3]`, all `catches("ERR-020")`, and `test_equal_volume_edges_bound_the_volumes`) are re-spelled onto `Mesh1D(g, CellsByCount.uniform_volume(n))` and keep their markers (`retirement-audit` C.13); the ERR-020 arm is re-run on them before the markers are believed | — |
| S3.4 equal width | `TestEqualWidth` (R2 edges; widths within `2 ulp(b)` of `fl(L/n)`; slab measure `fl(L/n)`, curvilinear shells between the realised edges; the slab discriminator) | — | — |
| S3.5 #495 | `TestIssue495`: partitions equal on a slab for n ∈ 1..64 ∪ {100, 127, 255, 1000} over four slab bodies; n = 5, 7, 9, 11 by name; the curvilinear negative leg; the max-width one-body law | the mesh leg `Mesh1D(g, uniform_width(n)) == Mesh1D(g, uniform_volume(n))` on a slab (needs S3.8's equality), and the from_geometry `"uniform"` rows (`test_single_region_slab_uniform`, `test_multi_region_slab`, `test_the_first_breakpoint_is_the_origin`) re-spelled onto `uniform_width` | — |
| S3.6 refinement | `TestRefinement` (k ∈ {2, 4, 8} × both spacings × three coordinates × 8 intervals × 70 counts; measures additive under refinement where they are `fl(m/n)`; max-width `2 * d` is the doubled realised count; the refuted halving exhibited) and `TestRefinementRefusals` (the factor, including 3 and 6; `2 * Partition`, `2 * CellEdges`) | — | — |
| S3.7 max width | `TestMaxWidthCount` (equal width: the least `n` with `fl(L/n) ≤ h`, checked by brute force over every smaller count; equal volume on a cylinder or sphere: the widest realised cell within `2 ulp(b)` of `h` and no smaller count fits; the discriminating rows (1, fl(1/3)) → 3, (3, fl(1/3)) → 9, (1, 0.1) → 10, (3, 0.3) → 10; one width per interval) and `TestCellsByMaxWidthRefusals` | — | — |
| rules' refusals | `TestCellsByCountRefusals` (no default spacing; no `CellsByCount.uniform`; a string spacing refused), `TestCellsByMaxWidthRefusals`, `TestCellEdges` | the re-posed `TestRegionMesh` rows of `test_structured_geometry.py` (S3.11: count ≤ 0 and non-int → `TestCellsByCountRefusals`; the unknown `method=` row is DEAD, deleted with `RegionMesh`) | — |
| posing seed 2 | `TestEveryCellInExactlyOneInterval` (15 rows: five rules × three coordinates; the owner decided by grouping and by containment, independently) | the mesh leg: `mesh.mat_ids` is each interval's material broadcast over its cells (S3.8(b)) | the discretisation digest: step 5 (ruling 4) |
| S3.8 the one constructor | — | all of (a) to (e) and its refusals; (a)'s hash leg moves to step 5 (ruling 4: the geometry is unhashable until `BC` is); (c) is `mesh.volumes is mesh.partition.all_measures` and `mesh.edges is mesh.partition.all_edges` (the orchestrator's note: both are stored fields of `Partition`) | the hash leg: step 5 |
| S3.9, S3.12–S3.14 | LANDED in step 2 | — | — |
| S3.10 retirements | — | the `Mesh1D` field legs (`bc_left`, `bc_right`, `precomputed_volumes` raise `AttributeError`) | the `None` refusal on the mesh path and the `boundary_condition=` signature leg, proven by capture 2 |
| S3.11 re-posed refusals | — | all five sites (spec §1.3) | — |
| S4.1 `from_homogeneous` | `TestFromHomogeneous`, `TestFromHomogeneousRefusals` | — | — |
| S4.2 k∞ legs | — | re-spelled through `from_homogeneous(w, BC.reflective)` with the new mesh constructor | — |
| ruling 3 constructors | `TestSlab`, `TestRadial`, `TestUniformBoundary` and their refusal classes: each constructor equals the bare constructor field by field over `dataclasses.fields(StructuredGeometry)`, laws by identity; each refusal is the bare constructor's own message (compared as strings, the route witness of one boundary check) | the test-site migration consumes them (N1's `_bcs_for` helpers retire onto `uniform_boundary`) | — |

**The step-3a mutation battery, run 2026-09-29** `[M]`. Instrument: `scratch/reference_architecture/p1step3/battery3a/battery3a.py`, a `-p` plugin installed at `pytest_configure` in the pytest process; the arm is chosen by `ORPHEUS_ARM`; each arm raises `Uninstallable` unless the mutant's answer on a probe input differs from the honest one (probes compare array BYTES: an ndarray `repr` truncates to 8 digits and read two arms as uninstallable before the fix), and prints `installed=1 rebinds=N`. Driver: `run_battery.sh <arm>...`; logs in `battery_out/`. Scope: the two new files (253 rows); `test_structured_geometry.py` is excluded because its mesh rows still build through `from_geometry` and `_subdivide_zone`, which no arm touches. Each run takes about 2 s. Every arm reddened its target row; 0 collection errors in every run.

| arm | mutation (installed in-process) | target rows | red | verdict |
|---|---|---|---|---|
| control-err020 (positive control) | every measure re-derived from the realised edges (`compute_volumes_1d`), ERR-020 itself | S3.3 measure rows | 19 | red on all three coordinates: `test_measures_are_the_interval_measure_over_n[×3]`, `test_multi_interval_measures_are_per_interval[×3]`, keystone ×3, the slab equal-width measure row and its discriminator, the four named #495 rows, measure additivity ×3, the CellEdges slab row. DECLARED BLIND: `test_equal_width_is_equal_volume_on_a_slab` stays green, because the mutant corrupts both rules identically (the equality law's stabiliser); the named rows see it through the distinct-measure count |
| unpinned-ends | the spacing body without `edges[0], edges[-1] = a, b` | the pinning witness row | 2 | red (cylinder, sphere). DECLARED BLIND: the population rows of `TestNesting` stay green, the population holding 0 witnesses |
| measure-from-origin | equal measures `m(0, b)/n` | S3.2, S3.3 | 14 | red wherever `a > 0` |
| interval-measure-association | `c·b^d − c·a^d` | S3.3 bitwise, keystone | 4 | red on cylinder and sphere through the added `(0.3, 0.7)` witness (premise 6) |
| r1-linspace-width | equal width placed by `np.linspace(a, b, n+1)` | S3.4 edges, S3.5 | 11 | red |
| slab-width-measure-is-diff | equal-width measures from the edges (today's `"uniform"`) | S3.4 slab measure, S3.5 | 15 | red |
| max-width-halved | `2 * CellsByMaxWidth(h)` as `CellsByMaxWidth(h/2)` | S3.6 max-width rows | 7 | red |
| count-ceil-realised | `n = ceil(L/h)`, raised until the realised widths fit | S3.7 rows `(1, 0.1)` | 18 | red |
| count-exact-rational | `n = ceil(Fraction(L)/Fraction(h))` | S3.7 rows `fl(1/3)` | 15 | red |
| widest-is-realised-first-cell | the nominal widest cell read off the realised first cell for every spacing | the one-body max-width row `(0, 3), h = 0.3` | 7 | red |
| any-factor | the power-of-two clause removed | the factor 3 and 6 refusals | 7 | red |
| eq-ignores-measures | `Partition.__eq__` compares edges only | the equality `measure` row | 2 | red |
| no-nesting-check | `Partition.partition` returns `self` unchecked | interval-count and end-edge refusals | 5 | red |
| all-edges-duplicates-breakpoints | `all_edges = concatenate(edges)` | seed 2 | 16 | red on all 15 seed-2 rows and `test_derived_views` |
| cell-edges-slab-formula | `CellEdges` measures `diff(edges)` in every coordinate | the shell rows | 4 | red on cylinder and sphere |
| radial-swapped-laws | `cylinder`/`sphere` store `(outer, inner)` | hollow equals-bare | 2 | red |
| uniform-boundary-raw-hollowness | the pre-fix `uniform_boundary` | the Fraction witness | 2 | red (and green on the fixed tree) |
| homogeneous-one-face | `from_homogeneous` puts vacuum on the right face | S4.1 identity legs | 12 | red |
| homogeneous-no-own-check | `from_homogeneous` without its width check | S4.1 keyed refusals | 4 | red: the rows pin the verb's own message, not the bare constructor's |

**The ERR-020 markers.** The control arm re-drops ERR-020 into the production rule, and the two rows it reddens on all three coordinates, `TestEqualVolume::test_measures_are_the_interval_measure_over_n` and `TestEqualVolume::test_multi_interval_measures_are_per_interval`, now carry `catches("ERR-020")`. The keystone reddened too and carries none, since it retires at 3b. At 3b, after the six `from_geometry` catchers are re-spelled, the same arm is re-run against them; a re-spelled row that stays green loses its marker.

**3b owes, beyond the table:** the capture comparison of §3 over the migrated tree; the retirement audit of §4's step-3 row minus the 3c items; and a per-arm re-run of this battery (the gates of 3a must stay red under their arms once `Mesh1D` consumes the partition).

#### Round 2 (2026-09-29, after the review of `2c62f9f2`): the interval rules and the one measure

This block SUPERSEDES the 3a column of the placement table above and the round-1 battery table, both of which describe the `Partition` API that retired (`reference_cache.md`, "P1 step 3a, after review": `Partition`/`PartitionRule` retire; rules are interval rules `rule.cells(geometry, (a, b)) -> (edges, measures)`; the measure has one correctly-rounded definition `coord.measure(edges) = c·diff(T(edges))`; whole-geometry realisation moves to the `Mesher` in 3b). The round-1 text is kept as history.

**What moved.**
- `[REFUTED 2026-09-29]` "S3.1 nesting across intervals and posing seed 2 are 3a's": there is no free-standing multi-interval object in 3a any more. S3.1 keeps its PER-INTERVAL legs in 3a; the cross-interval legs (every breakpoint is a cell edge, every cell in exactly one interval, the flat views) are 3b's, on `Mesher(g).partition(rule).mesh`, together with the mesh-level #495 leg.
- `[REFUTED 2026-09-29]` the keystone "equal volume is bitwise `_subdivide_zone`" as a law (ruling 5: numpy's power replaces Python's scalar `**`). The row is deleted; its successor is the ERR-020 invariant (equal shares bitwise identical, each `m(a, b)/n` with `m` the one definition) plus the one-definition ROUTE row.
- `Partition`'s value rows (equality, hash RECORD, flat views) are deleted with the type; `CellEdges` is now the value type with rows of its own (value equality and hash, `-0.0` canonicalisation, copy and read-only).
- The refinement factor law is on `Refined` (the type), gated through the constructor and through the operator; `Refined` of `Refined` flattens by the product.

**The 3a rows, round 2** (`[M]` 255 passed under `.venv/bin/python -O -m pytest` in 1.2 s, `npx pyright` 0 errors; 21 `rests_on` targets, all collected node ids). `tests/gates/mesh/test_partition.py` (193 rows):

| class | rows | law |
|---|---|---|
| `TestTheOneMeasure` | 14 | the measure is `c·diff(edges**d)` on arrays (bitwise, 3 coordinates × 17 intervals × 3 counts); ROUTE: a decoy `CoordSystem.measure` moves `compute_volumes_1d`, `geometry.measure`, equal shares and realised shells by exactly its factor; `compute_volumes_1d` is a one-line delegate (RECORD, retires in 3b); the measure coordinate per system and its inverse; `MeasureCoordinate` refuses exponents 0, 4, −1, 1.5; `geometry.intervals` and `geometry.measure` |
| `TestNesting` | 11 | S3.1 per interval (both spacings × 3 coordinates × 17 intervals × 10 counts); the pinning's own witness row; returned arrays are fresh |
| `TestSumLaw` | 6 | S3.2 |
| `TestEqualVolume` | 9 | S3.3: shares bitwise `m(a, b)/n` (`catches("ERR-020")`); each interval its own share on `(0, .5, 1.5, 2)` meshed 5/7/11 (`catches("ERR-020")`); the R2 edges |
| `TestEqualWidth` | 7 | S3.4 |
| `TestIssue495` | 14 | S3.5 per interval: equal cells on a slab over 5 intervals (one with a negative origin) × 68 counts; n = 5, 7, 9, 11 named; the curvilinear negative leg; `CellsByMaxWidth` one body on a slab |
| `TestRefinement`, `TestRefinementRefusals` | 36 + 26 | S3.6: nesting for k ∈ {2, 4, 8}; equal shares additive; `2 * (4 * r) == 4 * (2 * r) == Refined(r, 8)`; the `4 * (2 * d)` row (qa F3); the refined max-width rule; the refuted halving; factor refusals (0, −2, 3, 6: power of two; 2.0, `True`: int) through the operator AND the constructor; `Refined(CellEdges)` and `2 * CellEdges` refused |
| `TestMaxWidthCount` | 28 | S3.7 as in round 1, plus qa F2: a width below the interval's float resolution and an estimate above 2**31 cells are refused at once (5 rows) |
| `TestRuleRefusals` | 28 | `CellsByCount`, `CellsByMaxWidth` (a tuple count or width is refused: the per-interval choice is the Mesher's), `CellEdges` (qa F5: strings and bools refused), one disjointness row each |
| `TestCellEdges` | 9 | the geometry's measures between the given edges; curvilinear equal width is its own edges, the slab is not; value equality and hash; qa F4 `-0.0`; copy and read-only |
| two module rows | 5 | every rule is an `IntervalRule`; the spacings are values |

`tests/gates/geometry/test_named_face_constructors.py` (62 rows): round 1's rows, plus elegance F8 (`uniform_boundary("SLB", (−1, 1), …)` gives the bare constructor's keyed retirement refusal, message compared as a string), qa F4 (`slab((−0.0, 1.0), …)` equals the slab from `0.0`; `from_homogeneous(−0.0, …)` is refused as non-positive).

**3b owes, updated:** S3.1's cross-interval legs and posing seed 2 on the mesh; the mesh-level #495 leg; S3.8 on `Mesh1D(geometry, edges, volumes)` including the construction law of ruling 2 (each stored volume against `geometry.measure` of its cell within the derived bound; every breakpoint a cell edge) with its wrong-coordinate control; the `Mesher` (`partition(rule)` and `partition((r0, r1, …))`, `.mesh`, `refine(2)`, and `refine` nesting the current mesh); the six ERR-020 catchers of `test_structured_geometry.py` re-spelled through the Mesher and re-adjudicated by the control arm.

**The battery, round 2** `[M]` 2026-09-29. Instrument `scratch/reference_architecture/p1step3/battery3a/battery3a.py` (rewritten for the interval rules), driver `run_battery.sh`, logs `battery_out_round2/` (round 1's logs stay in `battery_out/`). Scope: the two files, 255 rows. 24 arms, every one installed (`installed=1`, a byte-compared precondition) and red on its target, 0 collection errors. Two arms first hung or would have: `count-ceil-realised` walked into the 2**31-cell rows because the mutant also dropped the refusals; the count arms now call the honest `count` first (keeping the refusals) and mutate only the count.

| arm | red | target read red |
|---|---|---|
| control-err020 (shares re-derived from the edges) | 16 | both ERR-020 rows on 3 coordinates, the slab width rows, the named #495 rows, share additivity. DECLARED BLIND: the #495 equality rows (both rules mutated alike) |
| unpinned-ends | 2 | the pinning witness rows (the population rows hold no witness) |
| measure-from-origin | 11 | shares and the sum law where `a > 0` |
| measure-association (`c·T(b) − c·T(a)`) | 8 | the one-measure rows, shares, shells, `CellEdges` (cylinder, sphere) |
| second-measure-spelling (`compute_volumes_1d` given its own body) | 2 | the ROUTE row and the delegate RECORD row |
| r1-linspace-width | 12 | S3.4 edges, S3.5 |
| slab-width-measure-is-diff | 16 | S3.4 slab measure, S3.5 |
| max-width-halved | 7 | the refined max-width rows, the halving row |
| count-ceil-realised | 17 | `(1, 0.1)` and the equal-width count rows |
| count-exact-rational | 15 | the `fl(1/3)` rows |
| no-resolution-guard (the pre-fix `ZeroDivisionError`) | 2 | the two below-resolution rows |
| widest-is-realised-first-cell | 7 | the one-body max-width row `(0, 3), h = 0.3` |
| any-factor | 9 | factors 3 and 6, both routes |
| factor-sum (qa F3's arm) | 9 | the `4 * (2 * d)` rows and the flattening row (0 red in round 1) |
| negative-zero-kept | 4 | the three `-0.0` rows and `CellEdges` (0 red in round 1) |
| cell-edges-np-array-parse | 5 | the string and bool entries (qa F5) |
| cell-edges-no-end-check | 2 | the wrong-ends refusal and its disjointness row |
| cell-edges-slab-formula | 6 | the geometry-measure rows (cylinder, sphere) |
| cell-edges-identity-eq | 2 | value equality and `-0.0` |
| measure-coordinate-any-exponent | 4 | the exponent refusals |
| radial-swapped-laws | 2 | hollow equals-bare |
| uniform-boundary-raw-hollowness | 2 | the Fraction witness |
| uniform-boundary-coord-late (elegance F8's arm) | 1 | the `(−1, 1)` string-coordinate row; the `(0, 1)` row stays green, the input on which the late parse was already right |
| homogeneous-one-face | 13 | S4.1 identity legs |

#### Step 3b's gates and the migration of `tests/gates/mesh` and `tests/gates/geometry` (2026-09-29)

The constructor ruled after round 2 is `Mesh1D(coord, edges, volumes, mat_ids, face_laws)` (option A: the mesh holds no geometry). The `Mesher` is the one call site that lifts a geometry: `Mesher(g).partition(rule | rules).mesh`, `.refine(k)`. S3.8's legs (a) "equal to the mesh of its partition" and (d) "`mesh.boundaries is geometry.boundaries`" are re-posed accordingly: (a) is bitwise equality of the value, and (d) is the lift, where each face law IS the geometry's law object.

New gate files:
- `tests/gates/mesh/test_mesh1d.py`, 46 rows:
  - S3.8 and S3.11: 19 keyed construction refusals and their disjointness row. The S3.11 rows are the old `test_geometry.py::TestMesh1D` refusals, re-posed; the material-id count is a law again under option A.
  - The volume law: a ±4-ulp volume is admitted, a 64-ulp one refused, on both the first and the last cell, and the wrong-coordinate control is refused in four pairs. The stored volume is not recomputed.
  - The value: bitwise equality over each field; the unhashable RECORD row for step 5; frozen and read-only; the inputs copied; the derived views; `with_distinct_cell_ids`.
  - S3.10: `bc_left`, `bc_right`, `precomputed_volumes`, `from_geometry`, `RegionMesh` and `_subdivide_zone` are gone; the `__init__` field list is exactly the five fields.
- `tests/gates/mesh/test_mesher.py`, 93 rows:
  - The lift, 5 rules (one of them a per-interval mixed tuple) × 3 coordinates × solid and hollow bodies: the mesh is the concatenation of each interval's `rule.cells`, with the materials repeated per interval and the laws taken by identity.
  - Posing seed 2 and S3.1 across intervals: every breakpoint is an edge at the count offsets, bitwise, and every cell lies in exactly one interval (by containment), the one whose material it carries.
  - Broadcast and per-interval rules give one mesh; the session chains.
  - `refine` nests, composes by the product, is refused for `CellEdges` and for a factor that is not a power of two.
  - Six keyed session refusals.
  - The #495 mesh leg: slab meshes are equal over 4 bodies × 68 counts, and the curvilinear negative leg holds.
  - A population control row.

**The migration** (`migration_brief.md`).
- Migrated in place: `test_mesh.py`, `test_bound_compat.py`, `test_boundary_factor_consumers.py`, `test_geometry.py`, `test_structured_geometry.py`, `test_hollow_inner_law.py`, `test_dropped_laws_are_refused.py`, `test_module_layout.py`.
- The ERR-020 catchers of `test_structured_geometry.py` keep their node ids and markers, and so do the rest of the `TestMesh1DFromGeometry` rows (the class was renamed `TestMeshingAGeometry` by the archivist's docs pass of 2026-09-29, with the error catalogue's and the CP page's citations updated in the same pass).
- Deleted or re-posed (brief rule 5):
  - `TestRegionMesh` (5 rows) is deleted. The count refusals are `test_partition.py::TestRuleRefusals`; the `method=` row is dead with `RegionMesh`.
  - `test_geometry.py`: the 4 `TestMesh1D` refusal rows are re-posed in `test_mesh1d.py`; `test_mesh1d_bc_defaults_none` and `test_mesh1d_backward_compat` are deleted, since their subject was the retired `None` default.
  - `test_dropped_laws_are_refused.py`: the two `[undeclared]` parameters are deleted, because `None` is refused at the mesh.
  - `test_hollow_inner_law.py::test_an_undeclared_inner_law_is_the_reflective_one` is re-posed as the `None` refusal.
- The capture gate `[M]`:
  - mesh tree: 364 passed; meshes `a` 29; 6 `pre_unmatched`, all in the deleted or re-posed `None` rows; 1019 `post_unmatched`, all in the two new gate files (991 in `test_mesher.py`, 28 in `test_mesh1d.py`); 4 law mismatches, the same four `None` rows.
  - geometry tree: 821 passed; `a` 39; `c` 5 (`test_mesh.py`, `linspace(0, 5, 6)` against the R2 edges, 1 ulp per edge, volumes unchanged); 3 `pre_unmatched`, all in the two deleted `None`-default rows; 0 law mismatches.

**The step-3b battery** `[M]` 2026-09-29. Instrument: `scratch/reference_architecture/p1step3/battery3a/battery3b.py` with `run_battery3b.sh`; logs in `battery_out_3b/`. Scope: `tests/gates/mesh`, plus `test_structured_geometry.py`, `test_geometry.py` and `test_mesh.py` (500 rows, `-m "not slow"`). 13 arms, every one installed with a byte-compared precondition and red on its target; 0 collection errors.

| arm | red | target |
|---|---|---|
| control-err020 (equal shares re-derived from the edges, through the mesher) | 27 | the six re-spelled ERR-020 catchers of `test_structured_geometry.py` (all red: their markers are re-adjudicated), the partition-level catchers, and the #495 named rows |
| lift-materials-by-position | 66 | the lift and seed-2 rows |
| lift-default-laws | 53 | the lift's law-identity leg |
| volume-law-off | 10 | the volume-law rows and the two volume refusals |
| volume-law-tight (0 ulp) | 99 | the band's admitted legs, and every mesher mesh whose equal shares are not the realised shells |
| volumes-recomputed | 34 | `test_the_volumes_are_stored_not_recomputed` and the ERR-020 catchers |
| eq-ignores-face-laws | 1 | the equality row |
| refine-forgets-its-factor | 12 | `test_refine_nests` (the second refine is the product) |
| mesh-before-partition-is-none | 2 | the session refusal and its disjointness row |
| slab-face-count-everywhere | 49 | the face-count refusals and every solid curvilinear build |
| no-negative-radius-check | 2 | the negative-radius refusal |
| retired-field-back (`Mesh1D.bc_left`) | 1 | the retirement row |
| slab-width-measure-is-diff | 20 | the #495 mesh leg (4 rows) and the partition-level slab rows |

#### Step 3b, review round (2026-09-29): the derived volume band, positivity, the span check, one boundary-law check, the axis adapter

The fixes of `review3b_qa.md` and `review3b_elegance.md` landed in production. The gates changed as follows.

- **The volume band** is `_volume_ulps(coord) = 2p + 5`: 7 on a slab, 9 on a cylinder, 11 on a sphere. The derivation is in its docstring; the flat `_VOLUME_ULPS = 8` is gone.
  - `TestTheVolumeLaw::test_the_band` asserts the formula per coordinate. It tests two edges of the band: `band − 1` ulp is admitted and `band + 1` refused, on the first and the last cell. Adding `k` ulp rounds by at most half an ulp, so the two edges are unambiguous.
  - New `test_the_worst_legal_output_is_admitted`: the worst measured equal-volume outputs are admitted, and each must stay within 1 ulp of its recorded gap. `[M]` 2026-09-29, macOS libm, 4000 random intervals per coordinate (scales 1e-3 to 1e3, n < 3000):

    | coordinate | interval | n | worst gap | band |
    |---|---|---|---|---|
    | sphere (the elegance case) | (0.06675753124937545, 0.9507586186918889) | 2400 | 7.504 ulp | 11 |
    | sphere | (43.0446014117623, 81.72750086937991) | 2677 | 7.460 ulp | 11 |
    | cylinder | (5.8405877264467225, 23.15910405796302) | 997 | 6.623 ulp | 9 |
    | slab | (7.13e-05, 0.00278) | 2827 | 3.082 ulp | 7 |

  - The claim kind is a THEOREM with a derived tolerance; the measured worst is the RECORD.
  - The slab's derivation is conservative by a factor of 2.3.
- **qa F4, positivity.** Two refusal rows sit on a cell one ulp wide, `[1, nextafter(1)]`, with volumes `0.0` and `−1e-300`. Both lie inside the absolute band of the measure, so positivity is the law that refuses them (fragment "a cell volume is positive").
- **Elegance Q3, the span check.** `TestEveryBreakpointIsAnEdge` uses two stub interval rules that pass the protocol door: their edges are shifted by one ulp at the right end and at the left end. They are refused on 3 coordinates with a keyed rule index, and there is a new session-refusal row.
- **One boundary-law check.** `TestOneBoundaryLawCheck` has two rows. The identity row asserts that `Mesh1D`'s `parse_boundary_laws` is the geometry's. The route row rebinds a decoy in every binding, needing at least 2, and both constructors must raise it.
- **The axis adapter.** New `tests/gates/mesh/test_axis_adapter_laws.py` (8 rows):
  - a hollow `RadialAxisMesh` gets a reflective inner face and its outer law by identity (SCOPE-BOUNDARY #511); the solid control gets no inner face;
  - an undeclared `AxisMesh` or radial outer law becomes reflective (ELEGANCE-DEBT #405), with a declared control.

**Battery, re-run** `[M]` (`run_battery3b.sh`, logs in `battery_out_3b_round2/`). 19 arms over 524 rows; every one installed and red on its target, 0 errors. New and re-posed arms:

| arm | red | target |
|---|---|---|
| volume-law-off | 10 | the band rows and the two volume refusals |
| volume-law-flat-8 (7 ulp everywhere) | 4 | the cylinder and sphere band rows and both sphere worst-case rows; the slab stays green, as it must, since 7 is the slab's value |
| volume-law-tight (0 ulp) | 103 | every mesh whose shares are not the shells |
| no-positivity-law | 3 | the two thin-cell rows and their disjointness row |
| no-span-check | 11 | the six stub rows, their session row and disjointness rows |
| second-boundary-check (a weaker copy in the mesh module) | 11 | the three count rows, the `None` rows and both route rows |
| adapter-inner-is-outer | 3 | the hollow inner-face row |
| adapter-undeclared-is-vacuum | 3 | the three undeclared rows |

The control arm now reds 31 rows, including the six ERR-020 catchers. The other arms are as in the table above; their counts moved only by the rows added.

### 1.4 Step 4: `from_homogeneous(width, boundary)`

The infinite medium realised by the test as a finite slab: one interval `[0, width]`, one material, both boundary points carrying `boundary`.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S4.1 | Construction law: one interval, a slab, `breakpoints == (0.0, width)`, both laws equal to `boundary`; refusals: `width ≤ 0`, a non-finite width, a `None` boundary | THEOREM | defining | pass the boundary to one side only → red | foundation; S2.1, S2.2 |
| S4.2 | IMPROVED, not new: the 1-D legs of `tests/gates/sn/solve/test_d3_admission.py::test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective[2g\|4g]` (SN) and `tests/gates/diffusion/test_solver.py::TestInfiniteMedium::test_reflective_slab_reproduces_k_infinity[2g\|3g]` (diffusion) are re-spelled through `from_homogeneous(w, BC.reflective)`; their assertions and tolerances are unchanged. The SN leg is a `None` site today (it relies on `SNProblem`'s reflective default) and migrates to the explicit law, which is what makes it a `from_homogeneous` consumer | REFERENCE (k∞ from the matrix eigenvalue, which never touches the sweep) | none: a re-spelling. Its reading must not move (§3) | the SN leg with `from_homogeneous(w, BC.vacuum)` → k drops below k∞ by far more than the tolerance (the arm the re-spelling must stay sensitive to) | L1 (SN, as today) and L2 (diffusion, as today); rests on the homogeneous k∞ gates `tests/gates/derivations/test_adjoint_spectrum_reference.py` (the dense pencil) |

`from_homogeneous` is a slab verb by its name's content (a width). A curvilinear "infinite medium" (a sphere with a white law) is not an infinite medium in continuous transport; it is not built.

### 1.5 Step 5: content identity

One encoder, reused from `Axis._structural_bytes` (`orpheus/numerics/axis.py:274`: type-tagged, length-prefixed, no `hash()`, no dict order) and `hashlib.blake2b`. The digest covers `Mixture` (over `_identity_key`), `Materials` (equality follows it; `eq=False` retires), `BC` (made hashable: `params` frozen), every registered boundary-law type, and the geometry.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S5.1 | Seed stability: a subprocess under `PYTHONHASHSEED=1` and another under `2` print the digests of one `Mixture` (mixture A, 4G), one `Materials` with two keys, `BC("albedo", {"albedo": 0.3})`, one hollow sphere geometry with an albedo inner law; the two outputs are equal. Control leg in the same row: `hash(mixture)` differs between the two runs (proves the harness sees salting) | THEOREM | defining for the new digest; the control leg is `[M]` (explorer, `probe_mixture_hash.py`: 8190618255204151145 vs 4571189084716701260) | a digest built from `hash()` → red | foundation; none |
| S5.2 | Equal content gives equal digest, `==` and `hash` (the eq/hash contract), over pairs built independently: distinct array objects; `BC` params in two insertion orders; `{"albedo": 1}` and `{"albedo": 1.0}` (equal under `==` today `[M]`, so their digests must agree); a CSR block with int64 indices (`[M]` canonicalised to int32, equal); `Materials({0: m})` and `Materials.of({0: m})` and `Materials({0: replace(m)})` | THEOREM | `[M]` (`identity.py`): `Materials({0: m}) == Materials({0: m})` is False today; `hash(BC.vacuum)` and `hash(StructuredGeometry(…))` raise `TypeError` | a digest over `repr(params)` → the int/float leg reds; a digest over the dict's insertion order → the order leg reds | foundation; S5.1 |
| S5.3 | Any field change gives a different digest, quantified over `dataclasses.fields(T)` for each type (the population is the type, never a hand list: a field added later reds this gate until it is perturbed): each dense field of `Mixture` (`nextafter` on one entry), each Legendre block of `SigS` and `Sig2`, `chi`, `eg` (`None` against a grid); `Materials` (a relabelled key `{0: m}` against `{1: m}`, an added key, a changed value); `BC` (kind, each param value, an added param); every law type in the boundary registry (`boundary/__init__.py:175`, 7 keys: one default instance per key, each field perturbed); the geometry (coord, each breakpoint, each material id, each law). One leg on a manufactured mixture perturbs `SigT` by float32 rounding (a value change; mixture A's `SigT` is float32-representable, `[M]`, so that leg is vacuous there) | THEOREM | the field enumeration is red against any encoder that omits a field; positive control on today's `Mixture._identity_key`: dropping `eg` from it → the `eg` leg reds | drop one field from the encoder → exactly that field's leg reds (a per-field arm table) | foundation; S5.1 |
| S5.4 | Signed zero and NaN, per the convention ruled (NEEDS 1). Recommended: the shared encoder canonicalises `−0.0` to `+0.0` and refuses NaN (a cross section or law parameter that is NaN is not a value). Legs: `SigL = −0.0` against `+0.0` digest-equal and `==`; a NaN entry refused with a keyed `ValueError`; `BC("p", {"a": nan})` refused | THEOREM | `[M]` (`identity.py`): today `Mixture` with `SigL = −0.0` is UNEQUAL to the `+0.0` one (bytes), while `BC` treats them as equal | an encoder over raw bytes → the signed-zero leg reds | foundation; S5.2 |
| S5.5 | A law with no content identity refuses the digest: a `PrescribedInflow` whose `_source` is a function object with no content (the test-side `_ManufacturedFaceInflow`, `tests/gates/sn/verification/analytical/test_mms_declared_inflow.py:74`) raises a typed error keyed on the law; `ConstantInflowSource` and `NoSource` (frozen dataclasses) digest | THEOREM (defining refusal) | the refusal is defining; the input exists in the tree today (the test-side manufactured inflow) | a digest falling back to `id()` → the refusal leg reds (it returns a value) | foundation; S5.3 |
| S5.6 | The shared encoder is the only route (X4): a decoy rebinding of the encoder moves the digests of `Mixture`, `Materials`, `BC`, a law and the geometry, all five (a ROUTE gate, lessons §1) | THEOREM | a second encoder written for one type → that type's leg stays unmoved → red | as the first-red column | foundation; S5.3 |
| S5.7 | RECORD fingerprint: the digests of three canonical objects (mixture A 2G, a 3-interval hollow sphere geometry, a `Materials` of two mixtures) are pinned as hex literals, with the message "the digest bytes moved: every cache key of P3 is invalidated; re-pin only with the reason" | RECORD (designed to red on any encoder change, `vv-principles`, the producer fingerprint owed to frozen-byte consumers) | defining | any encoder edit → red | foundation; S5.3 |

Before the carve, a constant-`__hash__` arm (lessons §1, the identity family): run the suite subset that constructs `Materials` and `MaterialMesh` with `Materials.__hash__` returning 0, and read its red set as the list of consumers that key on `Materials` identity. `[M]` (grep, `orpheus/`): the in-process caches are `sn/loss_representation/__init__.py:493-494`, keyed by content tuples and by `SNProblem`; `MaterialMesh._contractibility_key` already folds each mixture's `_identity_key`, so `[R]` content equality of `Materials` changes no key there. Tests that assert identity, which content equality must not break: `tests/gates/diffusion/test_augmented_mesh.py:142`, `tests/gates/transport/test_material_mesh.py:157` and `tests/gates/sn/primitives/test_snmesh_materials_pr_typed_0.py:101` assert `is`, which stays true (the same object is passed through).

### 1.6 Step 6: the phase-space functions

`Symbolic`: a function `q(r, μ, φ)` per energy group, a tuple of SymPy expressions stored through `srepr`, with a declared density convention (per steradian). `RegionwiseConstant`: one value per (region, group), an angle-integrated rate (cm⁻³ s⁻¹), which is `Symbolic`'s isotropic, piecewise-constant special case. The projections onto a method's unknowns are P4.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S6.1 | `srepr` round trip and seed stability: `Symbolic.from_srepr(s.srepr) == s` and the digests agree, over a population: a polynomial in `r, μ`; `sin(πr)·exp(−μ) + cos(φ)√(1−μ²)·r`; a `Piecewise` with a 30-digit `Float` and a `Rational`; a 12-term sum. Symbol assumptions survive the round trip. The digests printed under `PYTHONHASHSEED` 0 and 1 in two subprocesses are equal | THEOREM | `[M]` (`srepr_probe.py`): 4 of 4 round-trip, assumptions kept, identical output under seeds 0, 1, 2, 3; so the law holds for SymPy, and the first red is a digest over `hash(expr)` or over `str(expr)` (the witness) | a digest over `hash(expr)` → the seed leg reds | foundation; S5.1 |
| S6.2 | The phase-space coordinates are a declared vocabulary: an expression whose free symbols are not a subset of `{r, μ, φ}` (the symbols the type owns, with their assumptions) is refused, keyed on the stray symbol's name | THEOREM (defining refusal) | the input exists: an MMS-style expression with free parameters `a0 … a11` (`srepr_probe.py`, "many") | drop the check → the stray-symbol leg constructs → red | foundation; none |
| S6.3 | The anisotropy predicate: isotropic iff, per group, `simplify(∂q/∂μ) == 0` and `simplify(∂q/∂φ) == 0`; an undecided case counts as anisotropic (the refusal-safe side), and the docstring says so. Positive control `1 + μ` → anisotropic; negative controls `2` and `r² + 1` → isotropic; trap rows `sin²φ + cos²φ` and `μ² + (1 − μ²)cos²φ + (1 − μ²)sin²φ` (Ω·Ω) → isotropic | THEOREM | `[M]` (`srepr_probe.out`): a free-symbols predicate calls both trap rows anisotropic; SymPy folds `μ² + (1 − μ²)` to 1 on construction, so that one is not a trap | the free-symbols predicate → the two trap rows red | foundation; none |
| S6.4 | `RegionwiseConstant` is the same function as the `Symbolic` it lowers to: sampled at points strictly inside each region, in each group, the `Symbolic`'s value is the table value divided by 4π, exactly for `Rational` tables and to 1 ULP for floats; `is_isotropic` is computed by S6.3's predicate on the lowered function (a route leg: a decoy predicate moves `RegionwiseConstant`'s answer), never a hard-coded `True`; the region convention at a breakpoint (which region owns `r_k`) is stated and one row pins it | THEOREM | defining | lower with `Q` instead of `Q/(4π)` → the value leg reds (mode 3, a missing factor) | foundation; S6.3, S6.5 |
| S6.5 | The density convention, fixed at the definition site (`coding-elegance` Pattern 7): for a `RegionwiseConstant` with value `Q` in a region, `∫_{−1}^{1} ∫_0^{2π} q dφ dμ == Q` exactly (SymPy, `Rational` `Q`), and `∫_0^{2π} q dφ == Q/2` (the per-unit-μ density on the 1-D measure, W = 2, `docs/theory/foundations/normalization.rst:36`) | THEOREM (SymPy, Branch 1) | defining. The probe's evidence that the convention drifts today: `solve_sn_fixed_source`'s docstring equation applies `/W` twice (explorer, source report, NEEDS) | lower with `Q/2` (the slab W) → the sphere-integral leg reds; lower with `Q` → both red | L0 (SymPy); none |
| S6.6 | The group axis: the number of expressions (or of table columns) is the group count; the specification checks it against the materials (S8.1 leg g) | THEOREM | defining | — | foundation |

No positivity law: a manufactured source is signed (`q = Lψ − Sψ` changes sign), so the cone law belongs to the consumers that need it (a Monte Carlo sampler in P4), which refuse there `[R]`. The boundary source `InflowSourceSpec` is the Γ₋ twin of this type (explorer, specification report §3); whether it moves onto `Symbolic` is not ruled and is not specified here.

### 1.7 Step 7: the question

The closed sum type `Eigen(k)`, `Eigen(c)`, `FixedSource(σ)`, `CriticalParameter`, each with forward and adjoint. The adjoint `FixedSource` holds a detector, a step-6 function. No `Eigen(α)`.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S7.1 | Closed: the cases are exactly the four named, asserted by one exhaustive `match` helper whose `case _` raises, and by the union's `typing.get_args`; an eigen question over `ALPHA_MAP` (or the name `"alpha"`) is refused, keyed | THEOREM | `[M]` `ALPHA_MAP` exists at L1 (`numerics/posing.py:83`), so an `Eigen(spectral_map)` constructor admits it unless refused | admit any `SpectralMap` → the alpha leg constructs → red | foundation; none |
| S7.2 | One vocabulary: `Eigen(k)`'s spectral map IS `K_MAP` (`is`, not a copy, X4). `Eigen(c)` names no spectral map of the pencil `(A, F)`: the c-eigenvalue splits the collision operator differently and is a different pencil (explorer, specification report §2); what it names is NEEDS 3 | THEOREM | `[R]` a local copy of the map (`SpectralMap(…, name="k")`) → the `is` leg reds | as the first-red column | foundation; none |
| S7.3 | The adjoint's arity mirrors L1: `Eigen` and `CriticalParameter` adjoints are nullary (`EigenPosing.H()` is nullary, `k† = k`); `FixedSource`'s adjoint takes one detector of the step-6 function type (`SourcePosing.H(detector)` is unary). Legs: the adjoint `FixedSource` without a detector → `TypeError`; with a bare `ndarray` detector → refused (an extensional, mesh-bound array); a forward `FixedSource` given a detector → refused; `inspect.signature` arity parity with the L1 `.H` per case; involution `q.adjoint().adjoint() == q` for the nullary cases | THEOREM | `[R]` a boolean `adjoint=True` field (the explorer's warning) makes the detector-less adjoint spellable; the arity leg is its red | a boolean flag → the detector-less leg constructs → red | foundation; S6.1 |
| S7.4 | The default: a question derived from a source is `FixedSource(σ = 1)` (F11: the physical multiplying-source question), and from no source `Eigen(k)`; the explicit and the derived spellings are equal and digest-equal | THEOREM | defining | default σ = 0 → red | foundation; S7.1 |
| S7.5 | Content identity: every case digests (the digest names `K_MAP` by its `name`, never through its lambdas, which have no content); `Eigen(k) ≠ Eigen(c)`, `FixedSource(1) ≠ FixedSource(0)`, forward ≠ adjoint, and two detectors differing in one coefficient give different digests | THEOREM | `[M]` a `SpectralMap` holds three lambdas (`posing.py:61`), so a digest that walks the map's fields cannot be written; the gate's red is a digest over `repr` of the map | digest the map's `repr` (it contains the lambdas' addresses) → the seed leg of S5.1's style reds | foundation; S5.1, S6.1 |

`CriticalParameter`'s parameter (which degree of freedom of the geometry varies: an outer breakpoint, a scale of all breakpoints, and for what target) is not ruled; its gates wait for NEEDS 3.

### 1.8 Step 8: the specification

`Specification(materials, geometry | None, question, source | None)`, with its digest. The role fields name the pairing: `source` enters the right-hand side; the adjoint question's detector is paired with the flux.

| id | gate | kind | first red | mutation witness | level, rests on |
|---|---|---|---|---|---|
| S8.1 | Composition refusals, keyed and disjoint: (a) no geometry and more than one material (the infinite medium has one); (b) a geometry material id missing from the materials; (c) materials with different group counts; (d) a `FixedSource` without a source; (e) an `Eigen` or `CriticalParameter` with a source (the multiplying-source question is `FixedSource(1)`); (f) an adjoint `FixedSource` with a source (its right-hand side is the detector); (g) a `RegionwiseConstant` whose region count differs from the geometry's intervals, or whose group count differs from the materials'; (h) a `CriticalParameter` without a geometry. Rule (b)-(c) is the rule `MaterialMesh._validate_materials` applies; it lives once, in the input layer, and `MaterialMesh` calls it (Pattern 2; `MaterialMesh` is in transport and cannot be imported by the specification) | THEOREM (defining refusals) | defining; leg (b) mirrors `MaterialMesh`'s existing refusal, whose message fragment is pinned by today's gates (grep the shortest distinctive fragment before re-wording it, `retirement-audit` B.11) | drop each clause → its leg reds | foundation; S5.2, S6.6, S7.1 |
| S8.2 | Digest: equal content → equal digest; each field perturbed (a material, a breakpoint, the question, one source coefficient) → different; seed-stable; the derived-default specification equals the explicit one (S7.4); a specification holding a `RegionwiseConstant` and one holding the `Symbolic` it lowers to on the same geometry are equal (one function, S6.4) | THEOREM | defining | drop the source from the digest → the source leg reds | foundation; S5.3, S6.4, S7.5 |
| S8.3 | The layer contract: in a fresh interpreter, `import <the specification module>` followed by building one specification leaves no module of `orpheus.transport` or of any L3 package in `sys.modules` | THEOREM (the plan's ban: derivations never import `MaterialMesh`) | a specification module importing `MaterialMesh` for its material check (the tempting spelling of S8.1 (b)) | that import → red | foundation; S1.1 |
| S8.4 | The assembly: `MaterialMesh(spec.geometry.mesh(d), spec.materials)` builds for a 2-interval, 2-group specification, and for `geometry=None` the test's lift `from_homogeneous(w, BC.reflective)` builds an `SNProblem` whose k equals the dense pencil k∞ within 10 × `keff_tol` (read from the configuration that drove the solve) | REFERENCE (the dense pencil, assembled in the test from the raw `Mixture` arrays, independent of every SN operator) | none (an integration rung) | an assembly that drops the second material → k moves | L1; S4.2, S8.1 |
| S8.5 | RECORD fingerprint of one canonical specification's digest, as S5.7 | RECORD | defining | any encoder edit | foundation; S8.2 |

**Role confusion, and what can catch it.** The source and the detector are one function type by ruling, so the type cannot tell them apart. Three confusions are spellable; two are refusals. (i) A forward `FixedSource` given a detector: refused by the question's arity (S7.3). (ii) An adjoint `FixedSource` given a source: refused by S8.1 (f). (iii) The right function in the wrong field (the forward source passed as the detector): no P1 gate can see it, because the two are values of one type with no units. Its catcher is a solver-level VALUE gate owed by P4: the response `⟨R, A⁻¹q⟩` compared with an independent reference computed with the roles as intended, on a fixture where `A` is not self-adjoint (streaming with vacuum boundaries, anisotropic scattering, upscatter) and `q ≠ R` with different supports, so that the swapped response `⟨q, A⁻¹R⟩` differs from it. The adjoint identity `⟨R, A⁻¹q⟩ = ⟨A^{−†}R, q⟩` itself cannot see the swap: it holds for either assignment of the two functions, so the swap is inside its stabiliser (Mode 12), and on a self-adjoint fixture the swapped response is even equal.

### 1.9 Proving each gate can fail: the batteries

One battery per step, a `-p` plugin installed at `pytest_configure` that rebinds its target in every `sys.modules` binding (re-installed per test, because `importlib.reload` undoes a rebinding, lessons §2), refuses to run when a rebinding count is 0 or when a mutant's answer equals the honest one on a probe input (a distinct `Uninstallable`), and prints its rebind count in its result line. Each arm's verdict is a row of a per-arm table (red set against the TARGET row, never a count), and each battery has a positive control: step 1, delete the whole `"mesh"` row (S1.1, S1.2 (a)-(d) red); step 2, the old `{SLB: 2, CYL: 1, SPH: 1}` table (S2.2 red); step 3, the ERR-020 mutant (S3.3 and S3.2 red); step 5, an encoder that drops every field but the type tag (S5.3 red on every leg); step 6, the free-symbols predicate (S6.3); step 7, a boolean adjoint flag (S7.3); step 8, a specification with no composition check (S8.1). Each arm names the gate it must red in the "mutation witness" column above; an arm that reds nothing is reported as such with the four causes ruled out (not installed, did not bite, a twin predicate guards the path, the fixture annihilates the degree of freedom). Scope: the new and re-spelled gate files of the step, plus `tests/gates/geometry/` and `tests/gates/sn/regression/`, which the positive controls red.

## 2. The #495 gate

`tests/gates/mesh/test_partition.py::test_equal_width_is_equal_volume_on_a_slab` (S3.5). Its first red is measured on today's spelling `[M]` (`issue495.py`, slab `[0, 3]`, one region):

| n | `"uniform"` distinct volumes | `"equal-volume"` distinct volumes | edges that differ | volumes that differ |
|---|---|---|---|---|
| 4 | 1 | 1 | 0 | 0 |
| 5 | 3 | 1 | 4 | 3 |
| 7 | 5 | 1 | 1 | 5 |
| 9 | 5 | 1 | 1 | 7 |
| 11 | 3 | 1 | 4 | 3 |
| 16 | 1 | 1 | 0 | 0 |

The dyadic counts are blind (n = 4, 16, and today's `test_single_region_slab_uniform` at n = 4 on a width 4). The gate therefore includes n = 5, 7, 9, 11 by name, not only a range.

How it reddens and flips in the tree. It lands at step 1 as a strict xfail spelled on today's API, `Mesh1D.from_geometry(g, region_meshes=(RegionMesh(n, "uniform"),))` against `RegionMesh(n, "equal-volume")`, with `reason="#495: 'uniform' stores diff(edges) on a slab"` and a companion row asserting the mark's strictness by introspection. At step 3 the migration re-spells it onto `uniform_width(n)` and `uniform_volume(n)` and removes the mark in the same commit. The flip is not a no-op (Mode 8(10)): the assertion reads the stored measures and the edges, both of which the production change determines. The commit body of step 3 records the strict xfail's XPASS under the new code, run once before the mark is deleted.

## 3. The re-baselines of step 3: the migration protocol

The census's rebuild probe covers 255 of the 450 direct sites (those whose edges evaluate statically) and none of the ~102 `from_geometry` sites. The protocol that proves each site bit-identical or explains it is a runtime capture, because meshes are built in fixtures, helpers, parametrised loops and at collection time, where no static reading reaches.

### 3.1 Capture 1: every `Mesh1D` the suite builds

`scratch/reference_architecture/p1probe/mesh_capture.py` is the prototype: a `-p` plugin installed at `pytest_configure` that wraps `Mesh1D.__post_init__` and `Mesh1D.from_geometry`, asserts its own installation (`installed=1` in its result line), and records per construction the test node id (or `<collection>`), the construction kind (direct, `from_geometry`, direct on a global `linspace` with materials decided by position), the coordinate system, the cell count, and today's edges and volumes. Its positive control (`test_capture_control.py`) reads the #495 slab as `eV` under R1 (volumes move, edges not), the irregular literal as `ev` (nothing moves) and the hollow split sphere as `EV`.

For step 3 the capture runs twice, at the pre-carve commit and after the migration, over the full canonical suite, keyed by (test node id, ordinal of the construction within that test). The comparison is bitwise on edges, volumes, material ids, coordinate system and resolved boundary laws. Each differing row must fall in exactly one class below with a magnitude inside its bound; a row outside every class stops the migration.

| class | sites | what moves | bound | why it is correct |
|---|---|---|---|---|
| (a) bit-identical | the majority (§3.4) | nothing | 0 | — |
| (b) the #495 fix | Cartesian equal-width meshes whose `diff(edges)` is not constant | volumes, to `fl(L/n)` | each measure within `2·(ulp(r_k)+ulp(r_{k+1}))` of the pre-carve one (§3.4) | S3.4: the stored measure is the rule's exact measure (the ERR-020 precedent) |
| (c) R1 to R2 edge rounding | direct `np.linspace` sites whose `linspace` differs from the affine body | edges, and the volumes that follow them on a cylinder or sphere | ≤ 2 ULP per edge | §0 item 2 (ruled) |
| (d) materials decided by position | 5 sites (census §(2)) | edges, per class (c), once the interfaces become breakpoints | as (c) | the honest geometry has a breakpoint at each material interface |
| (e) hollow curvilinear | 13 statically evaluated sites (§0 item 14) and any at run time | nothing (an explicit inner `BC.reflective`) | 0 | #511: SN computes the reflective cavity today |
| (f) `None` boundary declarations | every site with a `None` law | nothing in the solution; the law becomes explicit | 0 in capture 2 | §3.2 |
| (g) the 58 curvilinear uniform sites | migrate to `uniform_width`, never `uniform_volume` | as (c) at most | ≤ 2 ULP; a larger move is a mis-migration to `uniform_volume` | the census's rebuild: `"equal-volume"` differs on 58 of 190 |
| (i) CP, MoC and MC sites | every site whose mesh reaches `solve_cp`, `solve_moc` or `MCMesh` | nothing: a CP slab declares left := right, a MoC geometry is a solid cylinder with its outer law, an MC geometry declares `periodic` everywhere | 0 | S3.12 to S3.14: the unread laws are not realised, so declaring them changes no computation |
| (h) the one-ULP irregular site | `test_g_adjoint_reciprocity.py:244` | nothing (the literal breakpoints are stored) | 0 | S2.4 |

### 3.2 Capture 2: the resolved boundary law per face

A second plugin records, per test and per face, the law each consumer resolves: at `resolve_boundary_conditions` (SN, diffusion), `_apply_default_bcs` (the three SN fixed-source entries), `cp/solver.py:228`, `moc/geometry.py:308` and `mc/solver.py:184`. Before the carve it records which of the four defaults each `None` resolved to; after the carve the same face must carry the same law, now declared. This is the proof of class (f), and the only instrument that can see a `None` site migrated to the wrong law, since capture 1 compares no boundary. Prototype: `scratch/reference_architecture/p1probe/law_capture.py` (installs at 5 sites: `resolve_boundary_conditions` in every module binding it, `_apply_default_bcs`, `CPMesh._resolve_bc`, `MOCMesh._resolve_bc`; MC's resolution is inline in `MCMesh.__init__` (`mc/solver.py:172-184`) and needs an `__init__` wrap, not yet written). Its positive control `[M]` (`test_law_capture_control.py`, `law_control.jsonl`): one `Mesh1D` with both laws `None` resolves to `BC('reflective')` on both faces in `SNProblem`, to an injected `BC('vacuum')` in `solve_sn_fixed_source`, and to `BC('white')` in `solve_cp`, three laws from one declaration. The CP and MoC rows also record the ignored `bc_left`: CP drops a declared left law on a slab `[M]` (`cp_slab_left_law.py`: k = 1.8749980808246423 with a white, vacuum or `None` left law; the right law is read), the #511 pattern in another method, reported to the orchestrator. `[R]` The population is large: of the 435 direct test sites, `bc_left` is given at 274 and `bc_right` at 304 (census; lower bounds, 5 sites splat `**kwargs`).

### 3.3 Which gates read a moved mesh, and how each is re-derived

A frozen reference whose producer builds a mesh in class (b), (c) or (d) changes its reading. Such a reference is RE-DERIVED, never re-pinned and never given a wider tolerance:

1. regenerate it with its own generator at the post-carve commit;
2. check that the drift is floating-point non-associativity of the size the moved mesh explains: for a direct evaluation, `reduction depth × ULP` relative to the moved cell measures; for an iterative one, `iterations × condition × ULP` (`vv-principles`, the three criteria);
3. check that the reference's independent anchor, when it has one, still holds at its own tolerance (for example the SN regression corpus's closed-form rows);
4. record the before/after `max|a − b| / max|b|` beside the ULP count in the commit body (a nulp count alone is uninterpretable).

The frozen artefacts in `tests/` (`find tests -name '*.npz' -o -name '*.npy'`, `[M]`): `sn/regression/snapshots` 25, `sn/_data/finalize_reconstruction_448` 32, `sn/_data/bc_extraction_baseline` 21, `sn/_data/affine_carve_baseline` 12, `geometry/snapshots` 7, `sn/_data/affine_carve_converged` 6, `sn/_data/bc_extraction_2d_baseline` 3, `numerics/data` 2, `sn/_fixtures/wave_t_t3` 1, `sn/_fixtures/wave_t_t4` 1, and 1 more under `sn/`. Their producers' meshes:

- `tests/gates/sn/regression/_generate_snapshots.py`: `from_geometry` with the default equal-volume on slab, sphere and cylinder. Under R2 the `from_geometry` cases keep their affine edges and `fl(L/n)` measures, and the curvilinear cases keep `_subdivide_zone` verbatim: 14 of 15 constructions bit-identical `[M]` (§3.4); `slab_fixed_source_dd_n20` is a direct `linspace` slab and moves by class (b) under either rule, so it is re-derived.
- `sn/_fixtures/wave_t_t4/_capture_pre_t4_snapshots.py`: `Mesh1D(np.linspace(0, 4, nx + 1), …)` direct (class (b) or (c), depending on `nx`); the 2-D producers (`wave_t_t3`, 2-D octant snapshots, `_test_helpers.py:253, 457`) build `Mesh2D`, which keeps its constructor and does not move.
- `sn/_data/affine_carve_baseline` (`test_affine_carve_baseline.py:150-190`): a slab and a sphere on `linspace(0, 2)`, and the hollow cylinder `linspace(0.01, 2)`, all direct (classes (b), (c), (e)).

The per-artefact verdict (moves or not, and by how much) is read from capture 1 joined to the test that consumes each artefact; §3.4 records it.

### 3.4 Capture results

`[M]` Capture 1 over `tests/gates` at `26a4f95c`, serial, `.venv/bin/python -O -m pytest -p mesh_capture tests/gates` (2026-09-25, 20:25 to 23:28): 1 failed, 12 302 passed, 315 skipped, 56 xfailed in 10 796 s (2 h 59 min, the plugin's cost included); `installed=1`, 3493 `Mesh1D` constructions over 2944 tests (1 at collection). Summary: `summarise_capture.py` over `cap_full.jsonl`.

| kind, coordinate | constructions | R2 verdict (`e`/`E` edges bitwise same/moved, `v`/`V` volumes) |
|---|---|---|
| `from_geometry`, slab / cylinder / sphere | 269 / 169 / 169 | all `ev` (bit-identical) |
| direct, slab | 1697 | 1170 `ev`, 317 `eV` (class b), 205 `Ev` (class c, edges only), 5 `EV` |
| direct, cylinder | 463 | 373 `ev`, 90 `EV` (class c) |
| direct, sphere | 672 | 605 `ev`, 67 `EV` (class c) |
| direct on a global `linspace`, materials by position (class d) | 22 / 23 / 9 | all `EV` |

Totals under R2: 738 of 3493 constructions move, in 514 tests (136 test files); maximum relative edge move 3.05e-16 (one ULP), maximum relative volume move 5.69e-14. Under R1: 381 move, in 317 tests; the same maximum volume move. The volume bound of classes (b) to (d) is therefore not "≤ 2 ULP of the cell measure": a one-ULP move of an edge `x` moves a cell of width `w` by about `ulp(x)/w` relative (the 5.69e-14 is a 640-cell slab, `test_heterogeneous_transport.py::test_sn_2region_reflective_case`, `test_keff_slab.py::test_heterogeneous_absolute_keff`). The class bound in §3.1 is corrected to: every edge within 2 ULP, and every measure within `2·(ulp(r_k) + ulp(r_{k+1}))·|∂m/∂r|` of the pre-carve one. `[REFUTED 2026-09-25]` "the SN regression corpus is bit-identical under R2" (§3.3): 14 of its 15 constructions are, and `test_dd_regression[slab_fixed_source_dd_n20]` builds its slab directly on a `linspace` whose differences are not constant, so the #495 fix moves its volumes by up to 1.39e-15 relative under R1 and R2 alike (class b). That snapshot is re-derived at step 3 by §3.3's four steps; the other 14 cases stay a free bit-level anchor.

The failure: `test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck[cyl_2g_3reg_folded_4x8_dd_n40]` (per-group max `|Δφ_norm|` 0.1268 against a tolerance of 0.12). The plugin only reads `Mesh1D`'s output after production's `__post_init__` has run, so `[R]` it is not caused by the capture; `[M]` a plain re-run of the node without the plugin (`.venv/bin/python -O -m pytest <node>`, 580.7 s) FAILS the same way, so it is a pre-existing red at `26a4f95c`, not the capture's. The pre-carve baseline of step 1 must carry it characterised (this node, the 0.1268 against 0.12 reading), never counted, and `process-discipline` owes it a fix or an issue before feature work proceeds.

## 4. The retirement audit, per step

Each retirement runs `retirement-audit` A.1 to A.3 (graph callers, a text grep over code, tests, `docs/` without `_build/`, and `.claude/`, and direct constructors), then `dead_references` after the fix, and `sphinx-build -W -n` as a set difference over the touched pages. Every grep starts with a positive control on a known member, and every completeness claim is re-run in Python (`re` over `pathlib.rglob`), because `grep` here is ugrep (`code-search`).

| step | retires | the three searches must find | and beyond the symbol |
|---|---|---|---|
| 1 | the module paths `orpheus.geometry.mesh`, `orpheus.geometry.factories`, `orpheus.transport.mesh.axis` | 0 import statements naming them (S1.4 and the collection count); 0 string constants outside past-tense history; `dead_references` 0 new dead targets | the 46 docstring constants and 104 `.rst` lines (§1.1), the three `autoclass`/`automodule` directives in `docs/api/geometry.rst`, the 47 `.claude/` lines sorted by tense (B.4, E.19); the whitelist comment and the cycle essay in `test_layer_imports.py:142-160`, which name `orpheus.geometry.mesh` in the present tense |
| 2 | the kind tag (`"SLB"`, `"CYL"`, `"SPH"`) on `StructuredGeometry`; `_GEOMETRY_TO_COORD`; `_GEOMETRY_TO_N_ENDPOINTS`; `Region.outer_thickness_cm` | 0 `StructuredGeometry("SLB"\|"CYL"\|"SPH"` constructions (the spec census: SLB 53, SPH 28, CYL 14, `HEX` 1, `sph` 1 at construction); the string inside `getattr` and in `match` statements (B.6); test-side `_COORD_TO_GEOMETRY_TAG` (`tests/gates/cp/test_properties.py:47`) | the tag as prose ("geometry kind", the docstring tables of `structured_geometry.py`); derivations' own vocabularies V4-V11 stay until P4 (ruled) and are NOT in this search's zero |
| 3 | `RegionMesh`, `Mesh1D.from_geometry`, `precomputed_volumes`, `Mesh1D.bc_left`/`bc_right`, `origin=`, `None` as a boundary declaration, SN's `boundary_condition=` parameters, the CP/MoC/MC `or BC(...)` defaults, `_apply_default_bcs` | `RegionMesh` 114 calls in 46 files plus the aliased `_RM` (`cp/solver.py:933`, `moc/solver.py:144`), which a last-segment filter misses (explorer); `from_geometry` 102 text sites; `CollisionCache.from_geometry` is a homonym and stays; `replace(…, bc_left=…)` 5 sites (`field_census.out`: `sn/solver.py:140`, `geometry/mesh.py:465`, `test_mms_declared_inflow.py:167, 376, 432`); `getattr(…, "edges")` 1 (`tests/gates/sn/acceleration/test_dsa_low_order.py:67`); attribute loads of `bc_left`/`bc_right`: orpheus 9, tests 17 (receiver unresolved, a ceiling) | the false docstring of `resolve_boundary_conditions` ("uniform across methods"); `error_catalog.rst:788` (ERR-020's fix paragraph names the retired `mesh1d_from_zones` and `precomputed_volumes` in the present tense); `docs/theory/foundations/structured_geometry.rst` (8 lines); the stale "#290 P5 owns the typed refusal" comment in `geometry/boundary/_source.py` (#290 is closed); `Mesh1D.with_distinct_cell_ids` (a `replace(self, mat_ids=…)`, which becomes a geometry with one interval per cell) |
| 4 | none | — | — |
| 5 | `Materials(eq=False)` and its docstring's "until content identity joins" | the `eq=False` and the sentence; `MaterialMesh._law_key` (`material_mesh.py:94`), the workaround for the unhashable `BC`, which becomes the law's own hash (retire, and re-point its callers) | every `hash(...)`-keyed persistence: the Sood cache's `_canonicalize` (`sood_registry/cache.py:181`) refuses a `Mixture` today and keys on `hash`-free JSON; it is superseded in P3, not here |
| 6-8 | none | — | — |

## 5. Costs

`[M]` measured on the host `.venv`, serial.

- Every probe of this file runs in under 10 s; `maxwidth.py` in about 6 s (114 399 pairs); the SymPy integrals of S6.5 and the derivative checks of S6.3 in under 1 s each `[R]` from the probe's run.
- A cold `import orpheus.geometry`: 0.09-0.15 s over three runs; each subprocess gate (S1.3's new entries, S5.1, S6.1, S7.5, S8.3) costs about 1 s including the interpreter.
- The SN regression drift-set delta (S1.5 (b)): about 1.6 s (lessons §3).
- `pytest --collect-only -q tests/gates`: 8.75 s (P0).
- The slow items are the full-suite runs: step 1's count identity and step 3's two captures (each run before and after the carve). §3.4 records the wall time of one full capture run. Timing is never a gate.

No gate of this file is marked `slow`: every new gate runs in under 2 s. S4.2 and S8.4 solve small problems (2 or 3 cells, 2 or 4 groups), under 1 s each as today.

## 6. The ladder, as a graph

| rung | gates | rests on |
|---|---|---|
| foundation: the layout | S1.1-S1.4 | nothing |
| foundation: the geometry value | S2.1 → S2.2, S2.3, S2.4, S2.5 | S2.1 |
| foundation: the partition laws | S3.1 → S3.2, S3.4, S3.6 → S3.3, S3.7 → S3.5 | S2.1 |
| edge: the #495 law | S3.5 (equal width = equal volume on a slab; the foundations are S3.3 and S3.4) | S3.3, S3.4 |
| the constructor | S3.8, S3.10, S3.11 | S3.1, S2.2 |
| scope boundaries | S3.9 (SN inner surface, #511), S3.12 (CP slab, #513), S3.13 (MoC, #514), S3.14 (MC, #513) | S2.2 |
| interior, reused | S4.2 (k∞ on a reflective slab, SN and diffusion) | the dense-pencil k∞ gates |
| foundation: identity | S5.1 → S5.2 → S5.3 → S5.4, S5.5, S5.6, S5.7 | S5.1 |
| foundation: phase-space functions | S6.1, S6.2, S6.3 → S6.4 ← S6.5 | S5.1 |
| foundation: the question | S7.1 → S7.3, S7.4, S7.5; S7.2 | S6.1 |
| the specification | S8.1 → S8.2 → S8.5; S8.3; S8.4 | S5.3, S6.4, S7.5, S1.1, S4.2 |

A cycle: none. A test on no rung: none. The migration captures (§3) are instruments of the carve, not gates, and sit on no rung.

## 7. Refuted or rejected candidates

- **R1, `np.linspace` as the equal-width body:** rejected by ruling (§0 item 2). The structural reason: it gives equal width and equal volume two bodies on a slab, so #495 is fixed by a special case rather than by construction, and it moves the edges of every slab `from_geometry` mesh, the SN regression corpus included.
- **"Equal widths" as a bitwise law:** `[REFUTED 2026-09-25]` for this question (what can a partition gate assert bitwise?). The fact: realised edges cannot give equal differences in general (772 of 852 cases for `linspace`); the stored measures can, which is where the law lives.
- **`ceil(L/h)` on realised widths as `CellsByMaxWidth`'s count:** refuted (§0 item 7): it violates "every width ≤ h" at `L = 1, h = 0.1`.
- **The exact-rational count `ceil(Fraction(L)/Fraction(h))`:** rejected, `[R]`: it returns 4 cells for `L = 1, h = fl(1/3)` because the float `h` is below 1/3, which no user writing `h = 1/3` means. The nominal width `fl(L/n)` is correctly rounded and monotone in n, so the least admissible n is well defined in the arithmetic the rule uses.
- **`2 * CellsByMaxWidth(h)` as `CellsByMaxWidth(h/2)`:** refuted as a refinement (S3.6): an odd count does not nest.
- **A dtype leg for the `Mixture` digest:** refuted (§0 item 9): dtype is canonicalised at construction.
- **CSR index dtype as an identity hazard:** refuted (§0 item 9): `_read_only_csr` canonicalises it; the fact that survives is that the key reads canonical bytes.
- **A free-symbols anisotropy predicate:** refuted (S6.3): it calls `sin²φ + cos²φ` anisotropic.
- **A self-adjoint reciprocity gate as the role-confusion catcher:** refuted (§1.8): the swap is inside its stabiliser.
- **Keeping `pwr_pin_2d` in `geometry/factories.py`:** refuted (§0 item 12): it constructs `Mesh2D`, an import the new row forbids.

## NEEDS

1. **Signed zero and NaN in the shared digest encoder (S5.4).** Today `Mixture` separates `−0.0` from `+0.0` and `BC` does not. Recommendation: canonicalise `−0.0` to `+0.0` and refuse NaN, in the one encoder.
2. **Where `RegionwiseConstant`'s regions come from (S6.4, S8.2).** Either it holds region indices and becomes a function of `r` only when paired with a geometry (then S6.4 needs a geometry argument and the "one function" law lives at S8.2), or it holds its own breakpoints (then S8.1 needs an agreement refusal against the geometry's, a second copy of one datum). Recommendation: region indices, lowered by the specification.
3. **What `Eigen(c)` names and what `CriticalParameter` varies (S7.2, §1.7).** No c spectral map exists on `(A, F)`; `CriticalParameter` needs its degree of freedom (an outer breakpoint, a uniform scale of the breakpoints) and its target.
4. **The region convention at a breakpoint (S6.4).** Which region owns `r_k` in a lowered `Symbolic`: measure-zero for every projection, but point evaluation at a cell edge (the SN starting direction at `r = 0`, a face value) reads it.
5. **Whether the step-3 captures run over the full suite twice (before and after), as specified, or over a named subset.** The full suite is the honest census (§3.1); a subset needs its own completeness argument.

## Orchestrator's rulings on the NEEDS (2026-09-25, reported to the user)

1. **Signed zero and NaN in the shared encoder:** −0.0 is canonicalised to +0.0 before digesting; NaN is refused at digest time with a message naming the field. So `Mixture` and `BC` agree.
2. **`RegionwiseConstant`'s regions:** they are the geometry's own region indices. The specification checks that their count matches the geometry, and lowers them to per-cell values through the mesh's nesting of cells in intervals.
3. **What `Eigen(c)` names and what `CriticalParameter` varies:** `CriticalParameter` varies the outer extent, in cm (the user's ruling). The stored eigenvalue question is WAITING for the posing-sequence review (`.claude/plans/posing_sequence.md`, W5 dispatched 2026-09-25). Step 7's gates are respecified after that ruling; steps 1 to 6 do not depend on it.
4. **Breakpoint ownership:** none is needed. Cells nest inside intervals, so lowering a function integrates cell by cell and never asks which region owns a breakpoint. Evaluating a function exactly at a breakpoint is refused.
5. **Step 3's captures:** they run over the full suite, before and after, as a detached driver with a log (the canonical gate takes about 90 minutes or more, so it is not run in the foreground).

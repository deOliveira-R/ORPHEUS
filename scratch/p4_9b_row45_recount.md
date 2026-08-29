# P4.9b — rows 4–5 RECOUNT (post-P4.9a), 2026-08-28

Re-verification of `scratch/p4_9_design_measured.md` rows 4–5, mandated by that
memo's own banner (plan-authoring §7). Measured at HEAD `10314dfa` (post-P4.9a
compaction; P4.9a landed `ca852c44`..`7a0f434c`).

Method: Python `ast` + `re` over `pathlib.rglob` (AST authoritative for
membership; regex for the name/string surface). Every enumeration ran WHOLE —
no `| head`. Every filter states its predicate and positive control inline.

⚠ **Brief-path correction:** the dispatch brief located the production ctor at
`orpheus/numerics/coupled_system.py:441`. `[M]` that file has **0**
`StreamingOperator` hits; the site is `orpheus/sn/coupled_system.py:441`
(`return StreamingOperator(sn_mesh) + MultiplicationOperator(`) — same line
number, different package. The prior memo said bare `coupled_system.py:441`;
the brief resolved the package wrong.

---
## Deliverable 1 — Row 4 recount: `StreamingOperator` construction sites

**Predicate (ctor):** `ast.Call` whose `func` is `Name(id="StreamingOperator")`
or `Attribute(attr="StreamingOperator")`, over every `*.py` under `orpheus/`
and `tests/` (`pathlib.rglob`, run whole). **Positive control:** asserted
`("orpheus/sn/coupled_system.py", 441) in ctor` — passed. **Negative
controls:** `\bStreamingOperator\b` asserted NOT to match
`ReducedStreamingOperator` (no word boundary inside CamelCase) and asserted TO
match `isinstance(x, StreamingOperator)`.

### Counts (prior → now)

| quantity | prior (pre-P4.9a) | now (`10314dfa`) |
|---|---|---|
| ctor calls TOTAL | ~150 | **136** |
| production ctor | 1 (`coupled_system.py:441`) | **1** — `orpheus/sn/coupled_system.py:441` |
| test ctor | ~148 in ~40 files | **135 in 40 files** |

Delta ≈ −14 test ctor calls: P4.9a landed between the two measurements
(retired the M-M twin + reworked the twin-equivalence gates; the file count is
unchanged at 40, so the shrinkage is within-file, not whole-file retirement).

### Per-file ctor counts — top 10 test files (of 40)

| ctor calls | file |
|---|---|
| 25 | `tests/sn/operators/test_streaming_collision_operator.py` |
| 19 | `tests/sn/operators/test_streaming_operator.py` |
| 10 | `tests/sn/operators/test_ld_adjoint_deferral.py` |
| 9 | `tests/sn/operators/test_removal_form_matvec_sweep.py` |
| 7 | `tests/sn/sweep/core/test_phase_c_gates.py` |
| 5 | `tests/sn/operators/test_bc_extraction_2d.py` |
| 5 | `tests/sn/operators/test_capability_survival.py` |
| 3 | `tests/sn/_test_helpers.py` |
| 3 | `tests/sn/operators/test_apply_full_field_codomain.py` |
| 3 | `tests/sn/operators/test_g_adjoint_reciprocity.py` |

(7 more files at 3, 6 at 2, 20 at 1 — full list reproducible by the script's
predicate above. ⚠ `tests/sn/_test_helpers.py` carries 3 ctor sites — the
prior memo's "test-side builder is probably the migration's lever" already has
a natural home.)

### Non-call references (the §6b members a symbol-ctor grep misses)

| kind | count | sites |
|---|---|---|
| class definition | 1 | `orpheus/sn/operators/streaming.py:172` |
| subclass bases | 0 | — |
| imports | 46 (2 production, 44 tests), **0 aliased** | production: `sn/coupled_system.py:184` (module scope), `sn/loss_representation/__init__.py:191` (**TYPE_CHECKING block** — type-only, erased at runtime) |
| isinstance guard | 1 production | `streaming.py:576` — `StreamingCollisionOperator.__init__` REFUSES a non-`StreamingOperator` first arg (`raise TypeError`). Ctor-shape change does not touch it, but it is the class-name consumer that makes duck-typed surrogates fail loudly. |
| generic subscript | 1 production | `streaming.py:450` — `OperatorSum["FullField","FullField","StreamingOperator","MultiplicationOperator"]` (string forward-ref type parameter) |
| string form (`getattr` hazard) | 6 total | `streaming.py:115` (`__all__`), `:450/:573/:628` (forward-ref annotations), `tests/.../test_streaming_collision_operator.py:262` (`pytest.raises(..., match="StreamingOperator")`). **No `getattr`/`hasattr` string dispatch.** |
| docstring/comment/prose | ~127 word-occurrences | raw regex total 313 across 61 files vs 186 AST-visible — remainder is prose (docstrings citing `StreamingOperator(sn_mesh)` etc.). These stale when the ctor changes; sweep at migration time. |

**No aliased imports and no string dispatch ⟹ the AST ctor count is the
complete call-spelling surface**; the extra §6b members are the isinstance
guard, the generic subscript, the forward-ref strings, and ~127 prose mentions.

---
## Deliverable 2 — Row 5 recount: `.scheme` / `.pole_angular_closure` reads, classified

**Predicates.** Attribute reads: regex `([\w.]+)\.scheme\b` and
`([\w.]+)\.pole_angular_closure\b` over every line of every `*.py` (receiver
captured for triage). Positive controls: `sn_mesh.scheme.spatial_basis_per_axis`
→ recv `sn_mesh`; `type(self.sn_mesh.scheme).has_transpose_kernel` → recv
`self.sn_mesh` — both passed. Negative controls: `.scheme_name` and `.schemes`
do NOT match (word boundary) — passed. String form: `grep -rnE
"['\"](scheme|pole_angular_closure)['\"]"` (the `getattr` hazard). Enclosing
defs resolved by AST span walk; code-vs-prose by reading each hit line.

**⚠ Filter hazard measured, worse than briefed:** `.scheme` is not only a
common word — it is a MODULE (`orpheus/transport/spatial/scheme.py`), so **74
of 160** raw orpheus/ hits have receiver `orpheus.transport.spatial`: import
lines and Sphinx `:meth:`/`:class:` xrefs to the module path. Not mesh reads;
excluded. A further 8 mesh-receiver hits are prose/comments (`__init__.py:236`,
`pairing.py:78`×2, `closure.py:579`, `_bases.py:235`, `augmented_mesh.py:266/:283`,
+1). ⚠ One crude-filter correction: `augmented_mesh.py:1235` (starts with `*` —
starred expression) is CODE, initially mis-tagged prose; corrected below.

### Headline totals (code reads off a mesh-typed receiver, orpheus/)

| bucket | count | prior (pre-P4.9a) |
|---|---|---|
| **consumer reads** (outside SNMesh + sweep_graph own-field) | **67** = 56 `.scheme` + 10 `.pole_angular_closure` + 1 string-form `getattr` | "~60 in 11 files" |
| sweep_graph OWN-field `self.scheme` (downstream fill; hazard (a)) | **6** — `sweep_graph.py:905,:921,:1022,:1023,:1108,:1109` (`_CellSolve.cell`, `_CellResidual.cell`, `_CellResidualTranspose.cell`) | 6 (unchanged) |
| SNMesh-INTERNAL (definition side, dissolves with the class change) | **7** — `:268` (fill `self.scheme`), `:422` (fill `self.pole_angular_closure`), `:587`×2 (`is_same_phase_space`), `:1230`, `:1235` (`angular_trial_space`), `:1363` (`boundary_face_layout`) | n/a |
| tests/ (total only; predicate: mesh-like receiver = `*mesh`/`sn_mesh`/`self.mesh`/`probe`/`sn2`) | 93 raw = 37 module-path + ~56 mesh reads (69 `.scheme` + 24 `.pole` raw) | not counted |

### ⛔ NEW FIND — a consumer INVISIBLE to the attribute regex (string form)

`orpheus/transport/fields/_bases.py:257` — `scheme = getattr(mesh, "scheme", None)`
inside the moment-tail widening derivation. Today `None` means "bare transport
carrier, scheme binds at augmentation" and RAISES TypeError. After row 5 an
SNMesh also has no `.scheme`, so **every moment-tailed (LD) production field
hits this raise** — loud, but the message becomes false, and no symbol grep
finds the site. Class (iv) REACH-THROUGH, duck-typed. This is the
`getattr`-string audit member (`coding-standards` clause) the prior memo missed.

### Classification: file → count per class

(i) CACHE/TABLE BUILD · (ii) PER-CELL/APPLY DISPATCH · (iii) SOLVER
ENTRY/ASSEMBLY · (iv) REACH-THROUGH · (v) OTHER

| file | (i) | (ii) | (iii) | (iv) | total |
|---|---|---|---|---|---|
| `sn/loss_representation/__init__.py` | 1 | 39 | 5 | 0 | 45 |
| `sn/solver.py` | 3 | 0 | 5 | 0 | 8 |
| `sn/operators/loss_kernel_gauge.py` | 0 | 0 | 4 | 0 | 4 |
| `sn/operators/streaming.py` | 0 | 2 | 0 | 0 | 2 |
| `sn/operators/radial_characteristic.py` | 0 | 2 | 0 | 0 | 2 |
| `sn/sweep/cache.py` | 1 | 0 | 0 | 0 | 1 |
| `sn/acceleration/dsa.py` | 0 | 0 | 1 | 0 | 1 |
| `sn/loss_representation/assembly.py` | 0 | 0 | 1 | 0 | 1 |
| `transport/radial_characteristic_field.py` | 0 | 0 | 0 | 1 | 1 |
| `transport/fields/_bases.py` (string form) | 0 | 0 | 0 | 1 | 1 |
| `derivations/continuous/mms/sn.py` | 0 | 0 | 1 | 0 | 1 |
| **totals** | **5** | **43** | **17** | **2** | **67** |

### Class (iv) — every hit (the reads that must be HANDED the closure)

- `transport/radial_characteristic_field.py:360` — `RadialCharacteristicField.source_from_angular`:
  `level_indices = mesh.pole_angular_closure.level_indices`. The L2 module the
  prior memo flagged — STILL PRESENT, unchanged.
- `transport/fields/_bases.py:257` — the `getattr(mesh, "scheme", None)` above
  (string form; L2 module).

### Class (ii) — every hit (43)

`sn/loss_representation/__init__.py` (39):
- `_LossRepresentation` helpers (⚠ construction/apply boundary — re-plumb =
  hand the representation the closure pair at construction): `:556`
  (`_n_face_moments`), `:573` (`_moment_frame_signs`), `:589` (`_spatial_moment_tail`)
- `_OctantWalk.loss_action` `:1135`; `.loss_action_transpose` `:1250`
- **cell-callable field FILLS** — `scheme=self.mesh.scheme` kwargs constructing
  `_CellSolve`/`_CellResidual`/`_CellResidualTranspose` (the SOURCE of
  sweep_graph's 6 own-field reads; hazard (a)'s upstream): `:1635`, `:1644`
  (`MovingFrontierWindow._sweep_interior`), `:1723` (`._loss_action_interior`),
  `:1807` (`._loss_action_transpose_interior`), `:2023`
  (`FullFieldWavefront._sweep_interior`), `:2094`, `:2184`
- `ScanMarch._sweep_interior` `:2405`; `._loss_action_interior` `:2502`;
  `._loss_action_transpose_interior` `:2621`
- `_OneDimScanWalk`: `._dag_legs` `:2868` (pole); `._degenerate_positions`
  `:2943` (pole); `._apply_walk` `:3147` (pole), `:3159`, `:3217`, `:3221`,
  `:3259` (**`sn_mesh.scheme.residual_kernel_batch(` — the matvec's per-cell
  kernel dispatch**); `.loss_action_transpose` `:3496`, `:3505`, `:3512` (pole),
  `:3513`, `:3539`, `:3541`, `:3543`, `:3554`; `._run` `:3921`, `:4040`
  (`source_emission`), `:4057` (`cell_average`), `:4212` (pole), `:4339`,
  `:4359`; `._run_transpose` `:4441`, `:4629` (pole)
- `_sweep_scheduled` `:4863` (apply-path layout read)

`sn/operators/streaming.py` (2) — dissolve into the operator's OWN closure
fields at row 4: `:277` (`StreamingOperator.is_adjointable` —
`type(self.sn_mesh.scheme).has_transpose_kernel`), `:992`
(`StreamingCollisionOperator._solve_timed_full_field` — `spatial_basis_per_axis`)

`sn/operators/radial_characteristic.py` (2): `:877`
(`RadialCharacteristicSeeding.apply` — `closure = mesh.pole_angular_closure`),
`:951` (`.apply_transpose`)

### Classes (i)/(iii) — representative lines

(i) CACHE/TABLE BUILD (5): `solver.py:1433` (`is_affine_scannable` guard on the
scan-cache build), `:1440` + `:1486` (`..., sig_t_1d, sn_mesh.scheme,` into
`CollisionCache.from_geometry` at `SNSolver.__init__` / `rebind_cross_sections`),
`loss_representation/__init__.py:3860` (`CollisionCache.from_geometry(geom,
sig_t_1d, self.mesh.scheme)` in `_ensure_coll_cache`), `cache.py:384`
(`closure = sn_mesh.pole_angular_closure` in
`StreamingCoefficientCache.from_mesh_and_quad` — post-P4.9a the closure MINTS
the scan constants; the cache reads the closure, no longer `.scheme`).

(iii) SOLVER ENTRY/ASSEMBLY (17): `loss_representation` strategy selection —
`:240/:244` (`_curvilinear_capability`), `:1380` (`CumprodScan.supports`),
`:2268/:2275` (`ScanMarch.supports`); `solver.py:474`
(`_residual_is_expressible`), `:920` (`_windowed_cold_start`,
`spatial_moments=` thread), `:3400/:3496/:3525` (rhs/source-lift/averaging
layout); `loss_kernel_gauge.py:375-419` (`gauge_freedom` — face transmission
spectrum + three message reads); `dsa.py:214` (`DSALowOrderSystem.from_sn_mesh`
— ⚠ `getattr(sn_mesh.scheme, "key", None)`, a stringly scheme-tag read);
`assembly.py:316` (`assemble_ordinate_blocks`); `derivations mms/sn.py:1601`
(`prescribed_inflow` layout — L3 tree).

### Prior-vs-now deltas (P4.9a landed between)

- `loss_representation/__init__.py`: prior "~45 `.scheme` + 6 pole" → now
  **39 + 6 = 45** (P4.9a shrank the scheme reads ~6).
- `cache.py`: prior counted as a `.scheme` reader → now **0 `.scheme`, 1 pole
  read** (`2340db08` — the closure mints its own scan constants; the builder
  now asks the CLOSURE).
- Everything else file-for-file unchanged: solver 8, gauge 4, streaming 2,
  radial_characteristic 2, dsa 1, assembly 1, mms 1, transport RC field 1.
- NEW members the prior audit missed: `_bases.py:257` (string form) and the
  explicit SNMesh-internal/own-field sub-buckets.
- Test-side surface exists but was not previously counted: ~56 mesh reads
  (+2 SNMesh-in-docstring refs; `test_unified_sweep_dispatch.py:691` pins
  "`SNMesh.scheme` defaults to DiamondDifference" — a gate that goes red at
  row 5 by design).

---
## Deliverable 3 — where the two closures are BUILT and CHOSEN today

**Construction surfaces (3, all funneling into `_init_core`)** — each takes the
same two selection parameters:

- `SNMesh.__init__(mesh, quadrature, materials, scheme=None, pole_angular_closure=None)` — `augmented_mesh.py:206`
- `SNMesh.from_axes(..., scheme=None, pole_angular_closure=None)` — `:687`
- `SNMesh.from_material_mesh(..., scheme=None, pole_angular_closure=None)` — `:764`

**Parameter shapes differ deliberately:**
- `scheme: DiscretizationSchemeBase | None` — an **INSTANCE**. Fill at
  `:268-269`: `self.scheme = scheme if scheme is not None else DiamondDifference()`.
- `pole_angular_closure: type[PoleAngularClosureBase] | None` — a **CLASS**
  (comment `:283-289`: a closure binds at construction and the mesh doesn't
  exist yet when the caller assembles ctor args). Stashed as
  `_user_supplied_closure`, instantiation DEFERRED past the coord dispatch.

**Closure instantiation (`:404-425`)** — post-Phase-B contract:
```
closure_cls = self._user_supplied_closure or default_angular_closure_class(self.coord)
angular, pairing = (self.reduced.angular, self.reduced.redistribution_pairing)   # 1-D reduced
                 | (angular_redistribution(self.quad, self.coord), zeros)        # multi-D Cartesian
self.pole_angular_closure = closure_cls(angular, pairing)
```
⭐ The closure takes the TWO TENSOR FACTORS, **not the mesh** — the block's own
comment: "the closure's operands were never mesh facts", buildable from
`(quad, coord)` (+ reduced data on 1-D). This is the P4.9b-relevant fact: the
posing head can mint the closure without a finished mesh.

**Default map** (`sn/angular/closure.py:2260-2277`, `default_angular_closure_class`):
CARTESIAN → `IdentityAngularClosure`; SPHERICAL, CYLINDRICAL →
`MorelMontryAngularSweep`. ⚠ DOC DRIFT: that function's docstring (`:2263-2265`)
still says "the factory dispatch (instantiation with ``sn_mesh``) is the
caller's job" — the pre-Phase-B contract; the real contract is
`cls(angular, pairing)`. Present-tense-false; fix at P4.9b.

**The concrete families (whole-tree census, `... , key=` subclass grep):**
- `DiscretizationSchemeBase` (RegistryMixin, ABC, `transport/spatial/scheme.py:894`):
  `DiamondDifference` (`diamond.py:165`, key `diamond_difference`),
  `LinearDiscontinuous` (`linear_discontinuous.py:243`, key `linear_discontinuous`).
  (`Step` at `scheme.py:918` is a DOCSTRING example, not code.)
- `PoleAngularClosureBase` (RegistryMixin, ABC, `sn/angular/closure.py:240`):
  `MorelMontryAngularSweep` (`:1504`, key `morel_montry_angular_sweep`),
  `IdentityAngularClosure` (`:2137`, key `identity_angular_closure`).
  (`MyClosure` at `:268` is a docstring example.)
- Registry string-key route: `.create(key)` has **0 production call sites**
  (grep `\.create(` × scheme/closure — docstring mentions only). Selection is
  by INSTANCE/CLASS object, never by tag.

**Direct constructions of the four concrete classes in `orpheus/`** (AST census;
predicate: `ast.Call` with func name in the 4-name set; positive control:
augmented_mesh's default fill — passed):
- **Exactly 1 name-keyed ctor call:** `augmented_mesh.py:269` `DiamondDifference()`
  (the default fill).
- The closure instantiation at `:422` is a call through the VARIABLE
  `closure_cls(...)` — the second construction site, invisible to a name-keyed
  census by construction (stated per the filter-honesty rule).
- `LinearDiscontinuous`: **0** production ctor calls — LD enters only when a
  caller passes `scheme=LinearDiscontinuous()` (tests/derivations do).
- `MorelMontryAngularSweep` / `IdentityAngularClosure`: **0** direct production
  instantiations outside the `closure_cls` route.
- Non-ctor production refs (4): closure.py's own default-map returns (2), and ⚠
  `loss_representation/__init__.py:4213` + `:4630` —
  `if not isinstance(closure, MorelMontryAngularSweep)` inside
  `_OneDimScanWalk._run` / `._run_transpose`: the walk still TYPE-DISCRIMINATES
  on the closure class (the degenerate/identity bailout). A residual dispatch
  site P4.9b's "operator holds the closure" design should either dissolve
  (polymorphic method) or knowingly keep.

⟹ **Where the posing choice enters today:** at SNMesh construction, through the
two kwargs on any of the 3 surfaces; the mesh is the composition site. No other
production code chooses a scheme or closure.

---

## Deliverable 4 — the cache/table build path

**`StreamingCoefficientCache`** (`sn/sweep/cache.py:125`, frozen slots
dataclass — Stratum 1, geometry+quadrature only). Builder:
`from_mesh_and_quad(cls, sn_mesh)` at `:258`.

**What the builder reads off the mesh** (`:258-411`): `sn_mesh.quad` (N, `mu_x`),
`.nx`, `.coord`, `.reduced` (⚠ `:273` — a bare `assert sn_mesh.reduced is not
None`: type-narrowing flavor, but stripped under the canonical `-O` runner),
`sn_mesh.dag_walk(ordinate_idx=...)` (iterates per-cell `StreamingTerms`,
unpacked once into chain-ordered tensors `chain_idx/A_down/A_total/dA_w/V`),
`quad.level_indices` (cylindrical arm), and — the one row-5 attr —
`sn_mesh.pole_angular_closure` (`:384`), from which it copies the four
closure-MINTED scan constants: `c_out_per_ordinate`, `c_in_per_ordinate`,
`tau_inv_per_ordinate`, `march_a_in_coeff_per_ordinate`.
**It does NOT read `mesh.scheme`** — post-P4.9a (`2340db08` family) the closure
mints its own constants (the `mm_a_in_coeff = (1-τ)/τ` derivation the prior
memo flagged at `cache.py:377` is GONE from the builder; the comment block at
`:370-383` records exactly this handing).

**`CollisionCache`** (`:418` — Stratum 2, σ_t-epoch). Builder
`from_geometry(geom: StreamingCoefficientCache, sig_t, scheme)` at `:479`:
takes the SCHEME as an **explicit parameter** (callers read it off the mesh —
the three class-(i) sites), and delegates the `(a, 1/denom)` recurrence
coefficients to `scheme.affine_scan_coefficients` (cache keeps storage +
`cumprod_a` scan-schedule transform; scheme owns the math).

**Construction sites, production (predicate: grep `StreamingCoefficientCache` /
`from_mesh_and_quad` / `from_geometry` over `orpheus/`, run whole):**

| cache | site | context | passes |
|---|---|---|---|
| Streaming | `sn/solver.py:1434` | `SNSolver.__init__`, guarded `:1433` `if sn_mesh.reduced is not None and sn_mesh.scheme.is_affine_scannable` | `sn_mesh` |
| Streaming | `sn/loss_representation/__init__.py:3791` | `_OneDimScanWalk._ensure_geom_cache` — lazy fallback, ⚠ memoized ONTO the mesh: `self.mesh._geom_cache = cache` (`getattr(self.mesh, "_geom_cache", None)` probe at `:3789`) | `self.mesh` |
| Collision | `sn/solver.py:1440` | `SNSolver.__init__` | `self.geom_cache, sig_t_1d, sn_mesh.scheme` |
| Collision | `sn/solver.py:1486` | `SNSolver.rebind_cross_sections` (σ rebind keeps Stratum 1) | `self.geom_cache, sig_t_1d, self.sn_mesh.scheme` |
| Collision | `sn/loss_representation/__init__.py:3860` | `_OneDimScanWalk._ensure_coll_cache` — lazy fallback for solver-bypassing callers (docstring: invariant is build-in-`SNSolver.__init__`, consume-per-sweep) | `geom, sig_t_1d, self.mesh.scheme` |

**Data flow (mesh → cache → walk).** The posing choice lands on the mesh at
`SNMesh._init_core` (`augmented_mesh.py:268` scheme, `:422` closure). At solver
posing, `SNSolver.__init__` (`solver.py:1433-1440`) builds Stratum 1 once per
solver lifetime (`StreamingCoefficientCache.from_mesh_and_quad(sn_mesh)`, which
walks `sn_mesh.dag_walk` and copies the closure's minted constants at
`cache.py:384`) and Stratum 2 per σ-epoch (`CollisionCache.from_geometry(geom,
sig_t, sn_mesh.scheme)`, delegating to `scheme.affine_scan_coefficients`).
Solver-bypassing callers get the same two objects lazily from
`_ensure_geom_cache` / `_ensure_coll_cache`
(`loss_representation/__init__.py:3787-3868`). The walk then consumes them as
explicit parameters — `_OneDimScanWalk._run(Q, sig_t, boundary_flux, geom,
coll)` (`:4399`) and `._run_transpose(..., geom, coll)` — whose scan fast path
is the closure's own representation `geom.tau_inv·ψ_avg − geom.mm_a_in_coeff·ψ_in`.
⟹ for P4.9b: the cache is ALREADY closure-fed (row 5's `.pole_angular_closure`
read at `cache.py:384` is a hand-over site, class (i)); the scheme reaches
Stratum 2 only through caller parameters, so re-plumbing rows 4-5 moves three
`(i)`-class reads (solver ×2, loss_rep ×1) plus the builder's one closure read.

---

## Compact conclusions for the P4.9b design

1. Row 4 is **smaller than briefed**: 136 ctor sites (not ~150), still 1
   production + 40 test files. No aliases, no string dispatch — the migration
   surface is exactly the AST list + 1 isinstance guard + 1 generic subscript
   + prose.
2. Row 5 consumer surface: **67 reads** = 5 (i) + 43 (ii) + 17 (iii) + 2 (iv),
   plus 6 sweep-graph own-field reads (filled at the 7 `scheme=` kwarg sites
   inside the walks) and 7 SNMesh-internal sites.
3. ⛔ The two class-(iv) reach-throughs are BOTH in `transport/` (L2):
   `radial_characteristic_field.py:360` (briefed) and the NEW
   `fields/_bases.py:257` `getattr(mesh, "scheme", None)` — invisible to
   attribute grep, breaks loudly-but-misleadingly post-row-5.
4. The closure ctor contract is already mesh-free (`cls(angular, pairing)`) —
   the posing head can own both choices without mesh-order gymnastics; the
   walk's two `isinstance(closure, MorelMontryAngularSweep)` sites are the
   residual type-discrimination to dissolve or keep knowingly.
5. Tests pin the current defaults: `test_unified_sweep_dispatch.py:691`
   ("SNMesh.scheme defaults to DiamondDifference") and two ctor-kwarg helpers
   (`test_angular_bulk_space.py:65`, `test_unified_sweep_dispatch.py:258`) go
   red at row 5 by design.

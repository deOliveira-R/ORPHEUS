# P4.9b — verification plan (PRE-carve)

**Authored** 2026-08-28 by `test-architect`, at `10314dfa` (post-P4.9a,
compaction point). **Scope** charter rows 4–5 as REVISED by the §7/§8 rulings
in `scratch/p4_9b_design.md`: the operator gains
`(sn_mesh, scheme, pole_angular_closure)`, `.pose` is the migration lever, the
hub sheds NOTHING, and the fused scan table becomes the strategy's lazily
resolved artifact.

Precedent for shape and rigour: `scratch/p4_9a_verification_plan.md`.

Every `[M]` below carries the probe that produced it. Probes live in
`/private/tmp/claude-501/-Users-rodrigo-git-nuclear-ORPHEUS/292c0918-53a4-478b-b616-9156e051dc68/scratchpad/`
(`probe_counts*.py`, `probe_cachecost.py`, `probe_identity*.py`,
`probe_readset.py`, `probe_swap*.py`, `probe_two_schemes.py`, `probe_poison.py`,
`probe_seam*.py`, `spy_p49b.py` — the activation census, and `mut_p49b.py` —
the M1 battery plugin; both assert their own installation).
All in-process monkeypatch / subclass; **no production file was edited on disk
at any point**. Nexus MCP was available and used (`provenance_chain`,
`dead_references`).

---

## 0. HEADLINE — twelve findings that change the brief

### F1 ⭐⭐ The step-2 §6c witness EXISTS, is BUILDABLE, and is RED TODAY — and it is the plan's keystone

The brief asks for four §6c witnesses, none of which is red-first for step 2
(the ctor-arity witness is trivially red; the pose-identity gate is green the
moment `.pose` is written; the lazy-table and memo-retirement gates are
post-hoc). The phase's actual claim — *"the walk's method-flavoured needs come
from the OPERATOR, not from the mesh"* — is a **route** claim (`vv` #26), and a
route claim needs a route instrument.

**The instrument: pose the operator, then SWAP the hub's bound objects for
deliberately mutant ones, and drive the solve through the ALREADY-POSED
operator.** Pre-carve the walk re-reads the hub at apply time, so the answer
MOVES. Post-carve the operator holds the pre-swap objects, so the answer must
be `array_equal`.

`[M]` (`probe_swap2.py`, 2-G heterogeneous, nx=8, `(L + C).solve(rhs)` on an
operator posed **before** the swap, mesh-attr memos dropped so the mutant is
actually consulted):

| config | hub object swapped | mutated surface (×1.05) | max \|Δ\| | rel | `array_equal` |
|---|---|---|---:|---:|---|
| SLAB `gauss_legendre(8)` | `mesh.scheme` | `source_emission` + `residual_kernel_batch` | 7.637211e-02 | **5.000e-02** | **False** |
| CYL `folded_product(4, 6)` | `mesh.pole_angular_closure` | `advance_psi_half` + `c_out_per_ordinate` | 1.458240e-01 | **4.596e-02** | **False** |
| CYL `folded_product(4, 8)` | same | same | 1.057619e-01 | **5.313e-02** | **False** |
| SPHERE `gauss_legendre(8)` | same | same | 1.246737e-01 | **1.196e-01** | **False** |

So the gate has a 5–12 % signal on every geometry, costs four sweeps, and its
red-today reading is measured rather than argued. **Spec in §2.1.**

⚠ **Three traps travel with it, all measured, all of which make it silently
green for the wrong reason:**

* ⛔ **`cell_contribution` alone is an INSUFFICIENT mutation, not a blind
  gate.** `[M]` mutating only `MorelMontryAngularSweep.cell_contribution`
  reads `array_equal = True` on all three curvilinear rows — because P4.9a's
  Q1 ruling split the surfaces: `.solve` consumes `advance_psi_half` plus the
  **minted scan constants**, while `cell_contribution` is the MATVEC's
  per-cell arm. A gate mutating one surface certifies one route.
  (`lessons` L49c: read a null arm as an insufficient mutation first.)
* ⛔ **`tests/sn/_test_helpers.sweep_once` is the WRONG driver.** `[M]` it
  builds `StreamingOperator(sn_mesh)` internally at `:814` — i.e. it re-poses
  **after** the swap, so post-carve it would still read the mutant and the
  gate would stay red for a reason that has nothing to do with the carve. The
  gate must drive `(L + C).solve` on the operator it posed itself. (Either
  that, or `sweep_once` grows an `operator=` parameter — a production-shaped
  change; see Q3.)
* ⛔ **The mesh-attr memos MASK the swap.** `[M]` without dropping
  `_geom_cache`/`_coll_cache` the cached table survives the swap and the gate
  passes *because of the cache*, not because of the operator. Post-carve the
  memo moves, but the gate must still carry an **activation leg** proving the
  mutant object is consulted at all on the pre-swap route (else it is `X == X`).

### F2 ⛔⛔ "Numpy assembly, measured-cheap" is REFUTED — the lazy table costs up to 24.65 % of a slab solve

The design memo (§5 rebind/perf note) says the operator-held table "re-mints
per `build_within_group_system` call (per solve) — numpy assembly,
measured-cheap; if contested, time it at execution." Timed now.

`[M]` (`probe_counts2.py` + `probe_cachecost.py`) — `StreamingOperator` is
constructed **6** times per slab eigenvalue solve and **10** on the sphere,
**independent of `nx`, `ng` and inner solver**; `default_for` fires exactly
once per operator (the `cached_property`), so a per-operator memo means one
table build per operator:

| config | scan-cache build (min of 5) | whole solve | 8 builds as % of solve |
|---|---:|---:|---:|
| CARTESIAN `GL8`, nx=8 | 0.22 ms | 169.6 ms | 1.05 % |
| **CARTESIAN `GL16`, nx=200** | **8.78 ms** | **284.8 ms** | **24.65 %** |
| CYLINDRICAL `fp(4,6)`, nx=8 | 0.36 ms | 668.0 ms | 0.43 % |
| CYLINDRICAL `fp(8,16)`, nx=200 | 39.10 ms | 10 246 ms | 3.05 % |
| SPHERICAL `GL16`, nx=200 | 8.81 ms | 1 471.9 ms | 4.79 % |

At the ruled 6 builds/slab-solve that is `6 × 8.78 = 52.7 ms` on a `313.6 ms`
solve — **16.8 %**. `[M]` today the count is **exactly 1** per solve on every
config (the eager `SNSolver.__init__` build + the mesh memo).

⟹ **The plan's perf gate is a COUNT, not a wall clock** (`lessons` L24/L25):
pin `StreamingCoefficientCache.from_mesh_and_quad` calls per solve. Today `1`.
Post-carve, whatever number the design chooses, it must be **written into the
gate**, so a memo-scoping regression is one red away instead of a silent 17 %.
**This is also Q1 (blocking): the memo's LIFETIME must span the operators
built within one solve, and a `cached_property` on a per-solve-rebuilt operator
does not.**

### F3 ⛔ Step 2's §6b call-site set is 79, not "~28"

The design memo sizes the representation signature change as
"`default_for(mesh)` → 1 production call at `streaming.py:407`, ~28 test calls
in 6 files". Measured by AST over `orpheus/` + `tests/` (positive control: the
production site asserted present):

| surface | count | files |
|---|---:|---|
| `default_for(...)` | **22** = 1 production + **21 tests** | 5 (`test_unified_sweep_dispatch` 10, `test_mms_ld_2d` 5, `test_mms_ld_slab` 4, `test_streaming_operator_decomposition` 2) |
| **direct representation ctors** (`CumprodScan` / `ScanMarch` / `MovingFrontierWindow` / `FullFieldWavefront`) | **58**, **0 production** | **11 test files** (`test_multi_d_reverse_walk` 24, `test_ld_adjoint_deferral` 8, `test_scan_march_equivalence` 6, …) |

`[M]` **every one of the 58 passes exactly one positional argument** (`sn`,
`sn_mesh`, `ld_mesh`, `_build_mesh("cart2d")`, …) and `_LossRepresentation
.__init__(self, mesh)` is the shared base signature for all four classes. So
if the representation gains two REQUIRED fields, all 58 break at collection.

⟹ Either (a) the 58 migrate too (making step 2's §6b set 79 sites in 14 test
files), or (b) the representations grow their own `.pose`-style classmethod, or
(c) the two closure fields are keyword-with-mesh-read defaults — which
re-instates the read the phase exists to remove and is therefore the wrong
answer. **Q2, blocking.**

### F4 ⛔ There are THREE mesh-attr memos, not one — and the REBIND contract's only witness dies with the solver slots

The design memo names `_geom_cache` (loss_rep:3791). `[M]` grep over
`orpheus/` + `tests/`:

| memo slot | stamped at | read at |
|---|---|---|
| `_geom_cache` | `solver.py:1444`, `loss_rep:3792` | `loss_rep:3789` |
| **`_coll_cache`** | `solver.py:1445`, **`solver.py:1488` (`rebind_cross_sections`)**, `loss_rep:3861` | `loss_rep:3855` |
| **`_pole_mirror_cache`** | `loss_rep:_ensure_pole_mirror` | same |

* `_coll_cache` is the σ-stratum and `rebind_cross_sections` **re-stamps it on
  the mesh** — that is what keeps a depletion/thermal-feedback rebind from
  serving a stale σ. `[M]` `_ensure_coll_cache` reads the memo and **never
  validates σ**: it is fresh only because the rebind re-stamps.
* `_pole_mirror_cache` is a third memo of the same idiom, derived from the
  quadrature's mirror deck. Untouched by the design memo.
* ⛔ **`tests/sn/sweep/core/test_cache.py:295-340`
  (`test_geometry_coefficients_invariance_under_sigma_t_change`) is the ONLY
  witness of the two-stratum rebind contract, and it asserts on
  `solver.geom_cache` / `solver.coll_cache` — the solver-held slots the design
  retires.** When those slots go, the test loses its subject (`lessons` L61h's
  DIE flavour). It must be re-posed at the new tier **in the same commit**, or
  the rebind contract ships unwitnessed.
* `[M]` one further consumer: `tests/sn/operators/test_loss_transpose_solve.py:401`
  calls `w._ensure_geom_cache()` directly.

### F5 ⛔ The 67 reads do NOT all re-plumb — 21 of them are SPACE-side facts the ruling puts hub-side

The design memo's step 2 says "the 43 (ii) + 17 (iii) + 5 (i) reads
re-plumb". The §7.1 ruling says the hub keeps the scheme **because** it induces
the space and supplies cross-consumers (DSA). Those two statements are in
tension, and the attribute census resolves it. `[M]` by attribute over
`orpheus/` (mesh-typed receivers, comments/prose excluded):

| attribute | count | class | verdict under the ruling |
|---|---:|---|---|
| `scheme` (bare — passed/bound) | 34 | mixed | per-site |
| **`scheme.spatial_basis_per_axis`** | **15** | space/layout | **STAYS hub-side** — the induced space fact |
| `pole_angular_closure` (bare) | 9 | mixed | per-site |
| **`scheme.is_multi_moment`** | **6** | space/layout | **STAYS hub-side** |
| `scheme.is_affine_scannable` | 4 | strategy selection | operator or hub — Q4 |
| `pole_angular_closure.level_indices` | 3 | per-cell | operator |
| `scheme.residual_kernel_batch` | 2 | per-cell kernel | **operator** |
| `scheme.source_emission` | 2 | per-cell kernel | **operator** |
| `scheme.cell_average` | 2 | per-cell kernel | **operator** |
| `scheme.residual_kernel_batch_transpose` | 1 | per-cell kernel | **operator** |
| `scheme.cell_kernel_batch` | 1 | per-cell kernel | **operator** |
| `scheme.transverse_coupling_is_facewise` | 1 | per-cell layout | operator |
| `scheme.supports_curvilinear` | 1 | strategy selection | operator or hub — Q4 |
| `scheme.face_transmission_spectrum` | 1 | diagnostic (`loss_kernel_gauge`) | **STAYS hub-side** |

`spatial_basis_per_axis` alone is **15 sites across 6 files**
(`solver.py` 5, `loss_representation` 6, `augmented_mesh` 1, `streaming.py` 1,
`derivations/mms/sn.py` 1, `_bases.py` 1 — the last a docstring). Re-plumbing
those through the operator inverts the ruling: a *layout* is a property of the
space the hub induces, not of the operator posed on it.

⟹ **Step 2's done-when must name the PARTITION, not the count.** A row reading
"the 67 reads re-plumb" is designed-red under the phase's own non-goals
(`plan-authoring` §10, third shape). §2.2 gives the executable form.

### F6 ⛔ The walk TYPE-DISCRIMINATES on the closure, and it refuses a proxy

`[M]` (`probe_readset.py`) wrapping `sn.pole_angular_closure` in a transparent
recording proxy makes the curvilinear sweep raise
`TypeError: _OneDimScanWalk curvilinear scan requires the Morel-Montry closure
(…)` — the residual `isinstance(closure, MorelMontryAngularSweep)` dispatch at
`loss_representation/__init__.py:4213` / `:4630` that the recount flagged
(Deliverable 3, note 4).

⟹ Any closure-side route/identity instrument must be a **SUBCLASS** of
`MorelMontryAngularSweep`, never a proxy or a duck-typed stand-in. F1's
`MutantMM` is built that way and passes. Record it so nobody re-derives the
failure. (It is also a live argument for dissolving that `isinstance`, but that
is O-3's item, not this phase's.)

### F7 ⭐ Step-1 activation is UNIVERSAL — the opposite of P4.9a; step-2's cache re-pose is NOT

`[M]` (`spy_p49b.py`, an activation census over `tests/sn/regression` +
`test_affine_carve_baseline.py`, 27 tests; the plugin `raise`s unless it binds
≥ 6 symbols — it bound 8):

| counted | total | tests with ZERO |
|---|---:|---:|
| `StreamingOperator.__init__` | 90 | **1 of 27** (`walk_matvec_baseline[cart2d_2g]`) |
| `default_for` | 90 | 1 of 27 (identical set — one per operator) |
| **`StreamingCoefficientCache.from_mesh_and_quad`** | **15** | **12 of 27** |
| `CollisionCache.from_geometry` | 15 | 12 of 27 |
| DD kernel batch | 407 670 | (2-D dominated) |
| `MorelMontry.cell_contribution` | 587 416 | (curvilinear) |

⟹ **The frozen corpus DOES pin step 1** (26 of 27 artifacts construct an
operator) — no P4.9a-style "build the bit-identity set" problem here. But
**the scan-cache re-pose is activated by only 15 of 27**: the 12 blind ones are
the `matvec_bulk_unmoved[*]` rows, the `walk_matvec_baseline[*]` rows and the
two 2-D DD regressions (which run the wavefront, not the chain scan). Cite the
15, never the 27, as the cache re-pose's anchor.

⭐ And a second activation asymmetry inside step 2, measured
(`probe_identity3.py`, full eigenvalue solves):

| config | per-cell **scheme** dispatch | per-cell **closure** dispatch |
|---|---|---|
| SLAB `GL8` | `residual_kernel_batch` **80** | **0** |
| CYL `fp(4,6)` | **0** | `cell_contribution` **14 496** |
| CYL `fp(4,8)` | **0** | `cell_contribution` **14 048** |
| SPHERE `GL8` | **0** | `cell_contribution` **3 192** |

⟹ **The two halves of step 2 have DISJOINT activating configs.** A curvilinear
fixture cannot witness the scheme re-plumb and a slab fixture cannot witness
the closure re-plumb. Every step-2 gate must carry both.

### F8 The `__eq__` change is real but small, and NOTHING hashes the operator

`[M]` `StreamingOperator` is `@dataclass(eq=True, frozen=False)` ⟹
**`__hash__` is `None`** — the operator is unhashable, so no dict/set/`lru_cache`
can be keyed on it today (any such site would already raise). And
`SNMesh.__eq__` is `object.__eq__`, so `StreamingOperator.__eq__` today is
**mesh IDENTITY**.

After the carve `__eq__` compares three fields. `[M]` `DiamondDifference.__eq__`
is dataclass VALUE equality (and it is hashable) while
`MorelMontryAngularSweep.__eq__` is object identity. So two operators over the
same mesh but with independently minted closures become UNEQUAL where they were
equal. `[M]` reachable only through the raw expert-seam ctor — under `.pose` the
instances are the hub's, so equality is unchanged.

`[M]` an AST sweep of every `Compare` node in the 61 `StreamingOperator`-touching
files finds **no `==`/`!=` between two operator instances** anywhere, and
`repr(op)` appears only inside `tests/_harness/predicates.py` failure MESSAGES
(never asserted). ⟹ the `__eq__`/`__repr__` change is inert for the shipped
suite; state it in the ctor docstring rather than gating it.

### F9 ⛔ The recount memo's row-5 red-by-design predictions are VOID

`scratch/p4_9b_row45_recount.md` Deliverable 2, note 5 says
`test_unified_sweep_dispatch.py:691` ("SNMesh.scheme defaults to
DiamondDifference") plus two ctor-kwarg helpers "go red at row 5 by design".
Row 5 was **revised by ruling** — the hub keeps both fields and its
`DiamondDifference()` default. `[M]` `TestDefaultDiscretizationScheme` at
`test_unified_sweep_dispatch.py:689-699` is untouched by the ruled design.
⟹ **Zero tests are predicted-red by design in P4.9b.** Any red in the step-1
scope is a defect, which makes the red-loop attribution free.

### F10 ⚠ A second `DiamondDifference` instance exists in every solve — and it is legitimate

`[M]` (`probe_two_schemes.py`) instrumenting `DiamondDifference.__init__` over a
slab eigenvalue solve: **two** instances are born.

1. `augmented_mesh.py:269` — the hub's default fill.
2. `transport/spatial/scheme.py:566` (`_face_transmission_spectrum`), reached
   from `loss_kernel_gauge.py:375 gauge_freedom` — a `@cache`d, **type-keyed**
   probe that deliberately constructs `scheme_type()` because the verdict "is a
   property of the closure and does not move with the probe cell".

`[M]` `cell_kernel_batch` dispatches on that throwaway instance (2 calls) while
`residual_kernel_batch` (80) and `affine_scan_coefficients` (1) dispatch on the
hub's.

⟹ Any "every consumer reaches ONE object" identity gate must carve this out
BY NAME, or it is a false red; and any instance-level mutation arm must know
that the diagnostic probe is invisible to it. (Also: `_face_transmission_spectrum`
is `@cache`d process-globally on `(type, ndim)` — a hazard for anything that
later keys on the scheme INSTANCE.)

### F11 ⚠ Naming: the brief's `op.angular_closure` does not exist in the ruled design

The brief's §6c(b) asks for `op.angular_closure is sn_mesh.pole_angular_closure`.
`[M]` the ruled fields (`p4_9b_design.md` §8 FINAL) are
`(sn_mesh, scheme, pole_angular_closure)` with `spatial_closure` as a *property*
over `scheme`. `[M]` `spatial_closure`, `angular_closure` and `.pose(` are all
**free tree-wide** (4 / 1 / 1 hits, all inside the campaign plan itself) — so
either spelling is available.

The asymmetry is the finding: `op.scheme` + `op.spatial_closure` (property) +
`op.pole_angular_closure` + *no* `op.angular_closure` reads as three
vocabularies for two slots, against
`feedback_naming_consistency_greppable` (one word-order pattern, role token
always present, so grep finds the whole set). ⟹ **Q5**: either mint the
symmetric `angular_closure` property in the same commit, or drop
`spatial_closure` until O-3 actually splits the family. Do not ship one of two.

### F12 ⛔ The ruling's own justification for expert-seam arm 3 is REFUTED — the doctored state RAISES, it is not silent

`[M]` `IdentityAngularClosure` on a curvilinear mesh CONSTRUCTS and then raises
at the first sweep — typed on the sphere (`TypeError: … requires the
Morel-Montry closure`), untyped on the cylinder (`IndexError: tuple index out
of range`) — and is bit-identically inert on the slab (rel deviation
`0.0000e+00`, because it IS the Cartesian default). The design memo's *"silent,
plausible-wrong k"* is not what the tree does; the same `isinstance` residual
that refuses a proxy (F6) refuses the wrong family.

**The ruling survives; the sentence must not be transcribed into the ctor
docstring.** Full table, the characterization test that freezes the ruling, and
the step-1→step-2 interval hazard: **§2.7**.

---

## 1. Step 1 — the operator's contract

### 1.0 What could still drift, given "bit-identical by construction"

The claim is that `.pose` flows the SAME instances, so nothing moves. `[M]` the
enumerated ways that could be false, each with its verdict:

| mechanism | verdict | evidence |
|---|---|---|
| `is_adjointable` re-points `type(self.sn_mesh.scheme)` → `type(self.scheme)` | inert under `.pose` (same object) | F8 |
| `loss_representation` cached_property re-selects | inert — `default_for` sees the same objects | F7 (`default_for` count == operator count) |
| dataclass `__eq__` widens from 1 field to 3 | inert for the suite | F8 (no operator `==` anywhere) |
| dataclass `__repr__` widens | inert — repr only in failure messages | F8 |
| `__hash__` | already `None`; unchanged | F8 |
| `streaming.py:992` `spatial_basis_per_axis` dissolves into own field | inert under `.pose` | F5 |
| a migrated site passes the two closures POSITIONALLY in the wrong order | **unspellable via `.pose`**; spellable at the raw ctor | §1.2 gate |
| `.pose` reads the wrong hub attribute (a SWAP) | catchable, and LOUD: `[M]` neither closure class carries `has_transpose_kernel` / `is_affine_scannable` / `spatial_basis_per_axis` / `is_multi_moment`, so a swap raises `AttributeError` at the first `is_adjointable` read — but only if that property is read, which the forward SI path does not do | §1.2 gate is the reliable catcher |

⟹ **Step 1's numerical acceptance is `array_equal` on the whole frozen corpus,
and it is genuinely activated** (F7: 26 of 27 artifacts construct an operator).
This is the P4.9a §F1 question asked and answered POSITIVELY — state it that
way in the close-out, because the previous phase's answer was the opposite and a
reader will assume it carried over.

### 1.1 §6c witness (a) — the ctor is unconstructable without both closures

**File** `tests/sn/operators/test_streaming_operator.py` (the ctor's own home,
19 existing ctor sites). `@pytest.mark.foundation`.

```
StreamingOperator(sn_mesh)                    -> TypeError
StreamingOperator(sn_mesh, sn_mesh.scheme)    -> TypeError
StreamingOperator(sn_mesh, sn_mesh.scheme, sn_mesh.pole_angular_closure)  -> constructs
```

⛔ **`match=` is mandatory on both negative legs, AND the fragment must be the
ARITY message.** A bare `pytest.raises(TypeError)` in this area is satisfied by
`StreamingCollisionOperator.__init__`'s own two `TypeError` guards
(`streaming.py:576`, `:583`) — `vv` Mode 8's self-satisfied-`raises` class one
step over. ⚠ And `[M]` the obvious fragment is **already taken and already
ambiguous**: `tests/sn/operators/test_streaming_collision_operator.py:262`
reads `pytest.raises(TypeError, match="StreamingOperator")` against the
`isinstance` refusal — and Python's own arity error is
`StreamingOperator.__init__() missing 2 required positional arguments: 'scheme'
and 'pole_angular_closure'`, which contains that same fragment. Match the
argument NAMES (`missing 2 required positional argument`), which pins the ARITY
rather than merely "it raises", and cannot be satisfied by either isinstance
guard.

**Red-today reading to commit in the docstring:** pre-carve
`StreamingOperator(sn_mesh)` constructs happily — `[M]` it is what all 136
sites spell today.

### 1.2 §6c witness (b) — the pose-identity gate (this REPLACES ctor guards)

The ruling is explicit that the raw ctor carries **no guards** and that
production consistency is *by construction* through `.pose`. That makes the
pose path the whole safety argument, so it needs its own gate.

**File** `tests/sn/operators/test_streaming_operator.py`, `@pytest.mark.foundation`.

Legs, all four:

1. **IDENTITY, spatial** — `L.scheme is sn.scheme`.
2. **IDENTITY, angular** — `L.pole_angular_closure is sn.pole_angular_closure`.
3. **NON-VACUITY** — a SECOND mesh over the same `(mesh, quad, materials)` has
   `sn2.pole_angular_closure is not sn.pole_angular_closure` (`[M]` true —
   the P4.9a is-identity gate already relies on it at `:144`), and
   `StreamingOperator.pose(sn2).pole_angular_closure is not sn.pole_angular_closure`.
   Without this leg an `is`-assertion could pass against a singleton.
4. **NO-SWAP** — `isinstance(L.scheme, DiscretizationSchemeBase)` and
   `isinstance(L.pole_angular_closure, PoleAngularClosureBase)`, asserted as a
   PAIR. `[M]` a swapped `.pose` reds legs 1 and 2 anyway; leg 4 makes the
   failure MESSAGE say which slot, and it survives a future `.pose` that
   returns copies.

**The mutation that reddens it:** `.pose` mints fresh objects
(`cls(sn_mesh, DiamondDifference(), default_angular_closure_class(sn_mesh.coord)(...))`)
instead of reading the hub. That is the realistic drift — it type-checks, it
solves correctly, and it silently breaks the "one instance ever" invariant the
DSA-consistency ruling rests on. Legs 1+2 red; nothing else in the tree does.

⭐ **Do NOT fold this into the P4.9a is-identity gate** — see §1.3.

### 1.3 The fate of `test_angular_closure_is_single_object.py`

`[M]` the gate asserts (`:132`, `:139`, `:149`) that every per-cell consumer
reaches **`sn.pole_angular_closure`** — the MESH's instance — and its control
leg (`:151-158`) asserts the fp(4,8) fast path fires the march recorder
**zero** times.

**Verdict: it SURVIVES unchanged in meaning, gains ONE leg, and its control
still binds.** Under the ruling the hub keeps the instance and `.pose` passes
it, so `sn.pole_angular_closure` remains the right identity target. Concretely:

* `:77` `L = StreamingOperator(sn_mesh)` → `StreamingOperator.pose(sn_mesh)` —
  one of the 135 migrated sites, no meaning change.
* **New leg to add, at the top of `test_every_per_cell_consumer_reaches_the_mesh_closure`:**
  `assert StreamingOperator.pose(sn).pole_angular_closure is sn.pole_angular_closure`
  — closing the loop between the operator's field and the identity the rest of
  the test asserts. Without it, post-carve the gate proves the consumers reach
  the HUB's object while the operator could hold a different one.
* **Control leg re-verification (mandatory, `[M]` reasoning + measurement):**
  the fp(4,8) zero-calls control patches `advance_psi_half` only. `[M]` at
  fp(4,8) `cell_contribution` fires **14 048** times per eigenvalue solve while
  the *march* recorder must stay at 0 — the control is about the march surface
  alone. Step 2 mints the scan table from the closure's **constants**
  (`c_out_per_ordinate` etc., attributes/`cached_property`), not from
  `advance_psi_half`, so the control survives. ⛔ **It must be re-run after the
  step-2 table re-pose and the number recorded** — a lazy table that mints by
  calling the march would silently break the standing constraint, and this is
  the only gate that would see it.
* ⚠ **Do not weaken the gate to `L.pole_angular_closure`.** Asserting the
  consumers reach *the operator's* field is strictly weaker post-carve: with the
  operator's field and the hub's being the same object it is the same claim, and
  the moment they diverge (the expert seam) the hub-keyed assertion is the one
  that still says something.

### 1.4 The migration — enumeration procedure and red-loop order

**Enumeration (§6b-complete, re-runnable).** AST, never `grep | head`
(`plan-authoring` §2's VIEWPORT clause):

```
predicate: ast.Call whose func is Name(id="StreamingOperator")
           or Attribute(attr="StreamingOperator"), over
           pathlib.Path(root).rglob("*.py") for root in {orpheus, tests}
positive control: ("orpheus/sn/coupled_system.py", 441) must be in the result
negative controls: \bStreamingOperator\b must NOT match "ReducedStreamingOperator";
                   must match "isinstance(x, StreamingOperator)"
```

`[M]` at `10314dfa`: **136 calls in 41 files = 1 production + 135 in 40 test
files**, all 40 under `tests/sn/`:

| dir | files |
|---|---:|
| `tests/sn/operators` | 24 |
| `tests/sn/solve` | 3 |
| `tests/sn/sweep/core` | 3 |
| `tests/sn/sweep/cartesian_2d` | 2 |
| `tests/sn/sweep/curvilinear` | 2 |
| `tests/sn/{_fixtures/wave_t_t4, ., architecture, regression, sweep, verification/mms}` | 1 each |

**The non-call §6b members and what happens to each** (recount Deliverable 1,
re-verified):

| member | site | disposition |
|---|---|---|
| `isinstance` refusal | `streaming.py:576` | **untouched** — a class-name check, blind to arity. But it is what makes a duck-typed surrogate fail loudly, so any test double for `L` must be a real `StreamingOperator` (`coding-standards`' algebra-routing clause). |
| generic subscript | `streaming.py:450` `OperatorSum["FullField","FullField","StreamingOperator","MultiplicationOperator"]` | **untouched** — a string forward-ref type parameter |
| production imports | `sn/coupled_system.py:184` (module scope), `sn/loss_representation/__init__.py:191` (`TYPE_CHECKING`) | untouched |
| `__all__` / forward-ref strings | `streaming.py:115`, `:450`, `:573`, `:628` | untouched |
| `pytest.raises(match="StreamingOperator")` | `test_streaming_collision_operator.py:262` | re-read it: `[M]` it matches the `isinstance` refusal's message, not the ctor's |
| prose | `[M]` **149** literal `StreamingOperator(sn_mesh|mesh|sn)` spellings across 54 files, minus 136 AST calls ⟹ ~13 prose; ~127 bare-name mentions | step 3 archivist sweep |

**Red-loop order** (cheapest-first, and chosen so a break is attributable):

1. `orpheus/sn/coupled_system.py:441` + the `.pose` classmethod + the ctor —
   then `tests/sn/operators/test_streaming_operator.py` (19 sites) alone,
   `[M]` the ctor's own file.
2. `tests/sn/operators` (24 files, 82 of the 135 sites) — `[M]` **1240 passed
   / 1 deselected / 5 xfailed / 74.0 s**.
3. `tests/sn/sweep/{core,curvilinear,cartesian_2d}` (7 files) — inside the
   `tests/sn/sweep` baseline, `[M]` **911 passed / 1 skipped / 8 deselected /
   31 xfailed / 282.2 s**.
4. The remaining 9 singleton files.
5. The whole 40-file set as the step gate — `[M]` **663 passed / 1 deselected
   / 19 xfailed / 56.2 s**.

⚠ Every un-migrated site is a **collection-time `TypeError`**, not a runtime
one, so a partial migration reports `ERROR`, not `FAILED`. Run with
`--continue-on-collection-errors` and count `^ERROR` separately (`vv` #17's
third pipeline failure; `lessons` L47d).

---

## 2. Step 2 — the operator feeds the walk and the strategy owns the table

### 2.1 ⭐⭐ Keystone: the hub-swap route gate (§6c witness, RED today)

**File** `tests/sn/operators/test_operator_feeds_the_walk.py` (new).
`@pytest.mark.foundation` — a software-invariant/route claim, no theory
`:label:`, so **no `verifies()`** (`feedback_vv_tagging`).

**Contract.** After `L = StreamingOperator.pose(sn)`, the answer produced by
driving `(L + C).solve(rhs)` is **independent of any later mutation of
`sn.scheme` / `sn.pole_angular_closure`** — because `L` holds the objects.

**Structure, per row:**

```
sn  = <fixture>
L   = StreamingOperator.pose(sn)
_   = L.loss_representation            # realize the selection now
base = drive(sn, L)                    # (L + C).solve(rhs) -> ndarray

# ACTIVATION leg (mandatory, run FIRST):
#   a freshly posed operator over the MUTANT hub must MOVE, proving the
#   mutant surface is consulted on this fixture at all.
assert not array_equal(drive(sn_mut, StreamingOperator.pose(sn_mut)), base)

swap sn's bound object for the mutant; drop _geom_cache/_coll_cache/_pole_mirror_cache
assert np.array_equal(drive(sn, L), base)      # <- the claim
```

**Rows** (both halves — F7's disjoint-activation finding):

| row | fixture | hub object swapped | mutant surface | `[M]` pre-carve rel |
|---|---|---|---|---:|
| `scheme_slab` | SLAB `GL8`, nx=8, 2G het | `sn.scheme` | `MutantDD(DiamondDifference)` overriding `source_emission` + `residual_kernel_batch` ×1.05 | **5.000e-02** |
| `closure_cyl_deg` | CYL `fp(4,6)`, nx=8, 2G het | `sn.pole_angular_closure` | `MutantMM(MorelMontryAngularSweep)` overriding `advance_psi_half` + `c_out_per_ordinate` ×1.05 | **4.596e-02** |
| `closure_cyl` | CYL `fp(4,8)` | same | same | **5.313e-02** |
| `closure_sphere` | SPHERE `GL8` | same | same | **1.196e-01** |

**Non-negotiable construction details, each of which is a measured trap
(F1):**

* The mutant closure MUST be a **subclass** of `MorelMontryAngularSweep`
  (F6 — the walk's `isinstance` refuses anything else), built as
  `MutantMM(sn.reduced.angular, sn.reduced.redistribution_pairing)`.
* `drive()` MUST NOT call `tests/sn/_test_helpers.sweep_once` — `[M]` it
  re-poses internally at `:814`. Build `(L + C).solve(TimedFullField(...))`
  in the test body from the operator the test posed.
* BOTH closure surfaces must be mutated. `[M]` `cell_contribution` alone reads
  `array_equal = True` on every curvilinear row.
* The three mesh memos must be dropped before the post-swap drive, or the
  cached table masks the swap. Post-carve the memo has moved; keep the drop
  anyway and say why (it must not be the reason the gate is green).

**Docstring text to commit (the §6c red-today record):**

```
⛔ RED BEFORE P4.9b step 2, by construction — and that is this gate's whole
point.  Pre-carve the representation reads ``self.mesh.scheme`` /
``self.mesh.pole_angular_closure`` at APPLY time (43 class-(ii) reads,
`scratch/p4_9b_row45_recount.md`), so an operator posed BEFORE the swap
still marches with the mutant.  [M] 2026-08-28 at 10314dfa, rel deviation
5.000e-02 (slab / scheme), 4.596e-02 (cyl fp(4,6) / closure), 5.313e-02
(cyl fp(4,8)), 1.196e-01 (sphere).  After the carve the operator holds the
pre-swap objects and every row must read np.array_equal.
⚠ ``cell_contribution`` alone is an INSUFFICIENT mutation: the .solve route
consumes ``advance_psi_half`` plus the minted scan constants (P4.9a Q1), so a
gate mutating only the matvec's per-cell arm reads array_equal=True and
certifies nothing.
```

### 2.2 The read-set gate — making F5's partition executable

The route gate (§2.1) pins the VALUE. This pins the CONTRACT: which hub
attributes the walk is still allowed to read.

**Instrument.** Wrap `sn.scheme` / `sn.pole_angular_closure` in a recording
descriptor after the operator is posed, run one sweep and one matvec, and
assert the recorded attribute set is a subset of the declared allowlist.

`[M]` the pre-carve read-sets (`probe_readset.py`, warmed then recorded):

| fixture | op | `mesh.scheme` attributes read | `mesh.pole_angular_closure` attributes read |
|---|---|---|---|
| SLAB `GL8` | sweep | `is_affine_scannable` ×2, `spatial_basis_per_axis` ×2, `source_emission` ×2, `cell_average` ×2 | — |
| SLAB `GL8` | matvec | `is_affine_scannable` ×2, `spatial_basis_per_axis` ×3, `is_multi_moment` ×2, `residual_kernel_batch` ×16 | `precompute_psi_state` ×1, `level_indices` ×2 |
| SPHERE `GL8` | sweep | `supports_curvilinear` ×3, `is_affine_scannable` ×2, `spatial_basis_per_axis` ×1 | `level_indices` ×2, `precompute_psi_state` ×1, `cell_contribution` ×8 |
| CYL `fp(4,6)` | sweep | `supports_curvilinear` ×3, `is_affine_scannable` ×2, `spatial_basis_per_axis` ×1 | `level_indices` ×2, `precompute_psi_state` ×1, `cell_contribution` ×32 |

⚠ **My own probe carried a vacuous arm and I am recording it rather than the
number**: the SPHERE/CYL *matvec* rows read `{}` because `do_matvec` guarded on
`if seed is None:` and `het_operands` returns a seed on curvilinear meshes — the
apply never ran. Those two cells are **not evidence**; the gate must assert its
own activation (a non-zero read count pre-carve) so the same silence cannot
recur inside the committed test (`vv` #17).

**Post-carve allowlist to gate** (the F5 partition, as a checkable predicate):

```
ALLOWED_HUB_SCHEME_READS  = {"spatial_basis_per_axis", "is_multi_moment"}
ALLOWED_HUB_CLOSURE_READS = set()
```

Everything else — `source_emission`, `cell_average`, `residual_kernel_batch`,
`precompute_psi_state`, `level_indices`, `cell_contribution` — must come from
the operator. `is_affine_scannable` / `supports_curvilinear` are **Q4**
(strategy selection: operator-side or hub-side?), and the gate's allowlist is
where that ruling becomes executable.

⭐ This gate is what stops step 2 from being over-scoped: it states the
partition as a set the tree can check, rather than as a count of reads that the
phase's own non-goals make unreachable (`plan-authoring` §10, third shape).

### 2.3 §6c witness (c) — the lazy table, and the COUNT gate F2 demands

**Two gates, and they answer different questions.**

**(c-i) EQUIVALENCE, in the transition commit only.** Build both — the retired
eager `StreamingCoefficientCache.from_mesh_and_quad(sn_mesh)` and the
strategy's lazily resolved artifact — and assert **`np.array_equal` field by
field**, over the 4-config grid `{SLAB GL8, SPHERE GL8, CYL fp(4,6), CYL
fp(4,8)}`.

⛔ `array_equal`, never `allclose`. `[M]` the P4.9a precedent measured the
realistic drift for exactly this kind of move (the cleaner algebraic spelling)
at **1–2 ULP**, so any tolerance ≥ 1e-15 is a non-catcher
(`p4_9a_verification_plan.md` §5.2). The fields are a gather over `dag_walk`
plus a copy of the closure's minted constants — reduction depth 0 ⟹
`array_equal` is the honest assertion, not a strict one.

⚠ **Then DELETE the eager build in the same commit.** A surviving eager builder
next to a lazy one is a Pattern-2 twin, and this gate would go tautological the
moment the lazy path is implemented by calling the eager one
(`coding-standards`' single-sourcing clause). The equivalence gate's honest
lifetime is one commit; what survives is (c-ii) plus the frozen corpus.

**(c-ii) COUNT — permanent.** A counting spy on
`StreamingCoefficientCache.from_mesh_and_quad`, asserting an EXACT number of
builds per solve.

`[M]` today: **exactly 1** per solve on every config (spy over full eigenvalue
solves, `probe_counts.py`), while `StreamingOperator` is constructed **6**
(slab, every nx/ng/inner-solver) to **10** (sphere) times. F2 prices the
difference at up to **24.65 %** of a slab solve.

The gate asserts the ruled number and names it in the message. If Q1 rules "one
build per solve", the gate is `== 1` and it is also the memo-lifetime gate. If
Q1 rules "one per operator", the gate is `== n_operators` and the docstring
must carry F2's table so the cost is a recorded decision rather than an
accident.

⚠ Scope: a wall-clock assertion is forbidden (flaky proxy, `lessons` L24/L25).
The count is exact, free, and it is the only instrument that can see a
memo-scoping regression.

### 2.4 §6c witness (d) — the memo retirements, all THREE of them

`[M]` F4: `_geom_cache`, `_coll_cache`, `_pole_mirror_cache`.

| leg | assertion | note |
|---|---|---|
| d-1 | after a full solve, `not hasattr(sn_mesh, "_geom_cache")` | the design's named target |
| d-2 | after a full solve, `not hasattr(sn_mesh, "_coll_cache")` | ⛔ **the design memo does not name it**; it is stamped from THREE sites incl. `rebind_cross_sections` |
| d-3 | `_pole_mirror_cache` — **ruling needed (Q6)**: same idiom, untouched by the design. Either retire it with the others or say in the plan why it stays | |
| d-4 | **the rebind contract, re-posed** — see below | |

⛔ **d-4 is the load-bearing one.** `[M]`
`tests/sn/sweep/core/test_cache.py:295-340` is the ONLY witness that a σ_t
rebind keeps Stratum 1 and rebuilds Stratum 2, and it asserts on
`solver.geom_cache` / `solver.coll_cache` — the slots being retired. Retiring
them without re-posing it **loses the contract silently** (`lessons` L61h's
DIE flavour: the gate's whole subject becomes unconstructible).

Re-posed form, at the surviving tier:

```
solver = SNSolver(sn_mesh=sn)
run one solve                                   # resolve the lazy artifacts
geom_before = <the strategy's resolved Stratum-1 artifact>
solver.rebind_cross_sections(sig_t * 2.0)
run one solve
assert <Stratum-1 artifact> is geom_before      # geometry survives the rebind
assert <Stratum-2 artifact> is not coll_before  # sigma stratum rebuilt
assert the post-rebind answer reflects the NEW sigma   # the staleness leg
```

⭐ The third assertion is **net-new teeth**. `[M]` `_ensure_coll_cache` reads
the memo and never validates σ — today it is fresh only because
`rebind_cross_sections` re-stamps the mesh slot at `solver.py:1488`. Once the
memo moves, "the rebind is seen" stops being structural and becomes something
that must be asserted. Write it in the commit that moves the memo, not later
(`coding-standards`: the retirement creates the exposure).

### 2.5 The `default_for` / representation signature — §6b (F3)

`[M]` 22 `default_for` calls (1 production + 21 tests in 4 files) **and** 58
direct representation ctors (0 production, 11 test files), all passing one
positional argument. Whatever Q2 rules, the enumeration procedure is:

```
predicate: ast.Call with func name in {"default_for"} ∪ {c.__name__ for c in LOSS_REPRESENTATIONS}
           over rglob("*.py") under {orpheus, tests}
positive control: orpheus/sn/operators/streaming.py must appear for default_for
```

⚠ `LOSS_REPRESENTATIONS` is a registry — derive the name set from it at census
time, never from a hand-written list, or a fifth representation is invisible.

### 2.6 The two L2 reach-throughs, and the ride-along

| site | today | disposition | gate |
|---|---|---|---|
| `transport/radial_characteristic_field.py:360` | `mesh.pole_angular_closure.level_indices` | HANDED the closure | covered by §2.2's allowlist (`level_indices` must leave the closure read-set) |
| `transport/fields/_bases.py:257` | `getattr(mesh, "scheme", None)` | ⭐ **SURVIVES** — the hub keeps `scheme`, so the `None` arm keeps meaning "a bare transport carrier". **No re-key, no message change** (the recount predicted one under the old row 5; that prediction is void with F9) | assert the raise still fires on a scheme-less carrier, with its CURRENT message — grep the shortest distinctive fragment first |
| `sn/sweep/cache.py:273` | bare `assert sn_mesh.reduced is not None` | ride-along → typed `raise` | see below |

**The `cache.py:273` ride-along needs its own witness, and `[M]` it has never
had one.** Under the canonical `-O` runner the assert is a NO-OP (`vv` Mode 8;
`coding-standards`' bare-assert clause) — it is a **domain/admission contract**
("the chain scan needs a `ReducedStreamingOperator`"), not type-narrowing, so it
must become a `raise`. `[M]` the guard's own message fragment
`"StreamingCoefficientCache requires"` returns **0 hits across `tests/`** ⟹
nothing asserts this refusal anywhere in the tree, so the conversion's teeth are
**NET-NEW coverage, not a migration** (`lessons` §1, the retiring-guard clause:
grep before crediting a mechanism swap as behaviour-identical).

Gate: a 2-D Cartesian mesh (`sn_mesh.reduced is None`) fed to the builder must
raise, with `match=` on the shortest distinctive fragment of the NEW message,
**under `python -O`**. Prove the conversion four lines the honest way: run the
guard's own arithmetic on the bad input under plain `python` and under `-O`;
if the second returns instead of raising, the contract is still inert.

### 2.7 The expert seam — what documents it, and ⛔ the ruling's stated reason is REFUTED

The ruling is *no operator-side guard*: the raw dataclass ctor is a **declared
expert seam** whose two residual arms (design memo §8, attacks 3 and 4) are
"recorded in the ctor docstring as what the seam does not check — articulation,
not enforcement". The ruling stands. **Its stated justification for arm 3 does
not**, and the docstring must not transcribe it.

`[M]` attack 3 (a wrong-FAMILY closure — `IdentityAngularClosure` on a
curvilinear mesh), built through the hub kwarg and driven through one sweep
(`probe_seam.py`, `probe_seam2.py`, 4 of 4 geometries):

| geometry | constructs? | first sweep |
|---|---|---|
| CYL `folded_product(4, 6)` | ✅ | ⛔ **`IndexError: tuple index out of range`** — loud, but untyped and undiagnostic |
| CYL `folded_product(4, 8)` | ✅ | ⛔ same `IndexError` |
| SPHERE `gauss_legendre(8)` | ✅ | ⛔ **`TypeError: _OneDimScanWalk curvilinear scan requires the Morel-Montry closure (…)`** — loud and typed |
| SLAB `gauss_legendre(8)` | ✅ | ✅ runs, **max rel deviation exactly `0.0000e+00`** — `IdentityAngularClosure` IS the Cartesian default, so this is not a hazard at all |

⟹ The memo's characterisation — *"τ≡1, c≡0 then close the pole with no
redistribution: **silent, plausible-wrong k**"* — is **not what the tree does**.
The same `isinstance(closure, MorelMontryAngularSweep)` residual that refuses a
proxy (F6) refuses the wrong family. Arm 3 is loud on every geometry where it
could matter and inert where it could not.

⭐ Two consequences:

1. **The ctor docstring must record the MEASURED behaviour**, not the memo's
   prose. A hazard note that says "silent" about something that raises teaches
   the next reader the opposite of the truth, and it is the kind of sentence
   that survives into a future audit as licence to add the guard the ruling
   forbids.
2. **The cylinder's `IndexError` is a worse diagnosis than the sphere's typed
   refusal** — same defect, two qualities of message, because the cylinder path
   indexes before the guard runs. Out of scope for P4.9b (the ruling forbids a
   new guard), worth a one-line issue.

**The characterization test that freezes the ruling** — recommended, cheap,
and it is the artefact that stops a well-meaning future guard:

```
tests/sn/operators/test_streaming_operator.py
@pytest.mark.foundation
def test_the_raw_ctor_is_a_declared_expert_seam():
    """The raw ctor accepts a doctored (mesh, scheme, closure) triple BY
    DESIGN — production consistency is by construction through .pose, and a
    guard here would forbid the seam's own use case (a diagnostic probe that
    deliberately switches redistribution off).  Ruled 2026-08-28.
    [M] the doctored state is not silent: it raises at the first sweep on
    every curvilinear geometry (TypeError on sphere, IndexError on cylinder)
    and is bit-identically inert on the slab.
    """
    L = StreamingOperator(cyl_mesh, cyl_mesh.scheme, IdentityAngularClosure(...))
    assert L.pole_angular_closure is not cyl_mesh.pole_angular_closure
```

⚠ It asserts **constructibility**, which is the ruling; it must NOT assert the
raise (that is the walk's `isinstance`, a different owner, and O-3 may dissolve
it). One positive leg, no negative leg — deliberately, and say why in the
docstring, because "a contract-validation test with no negative leg" is
otherwise a review flag (`vv` #11).

⛔ **Interval hazard, step 1 → step 2.** Between the two steps the operator
HOLDS the closures and the walk still READS the mesh, so a raw-ctor operator
carrying a different closure is **silently ignored** — arm 3's behaviour in the
interval is neither the memo's "silent wrong answer" nor the measured raise, but
"no effect at all". Land the characterization test at step 2, not step 1, or it
pins the interval rather than the contract.

---

## 3. The bit-identity set and its ACTIVATION

The P4.9a §F1 question, asked of every named artifact: *does it EXECUTE the
changed path?*

| # | artifact | pins | activates step 1? | activates step 2? | `[M]` |
|---|---|---|---|---|---|
| A1 | `TestAffineCarveSweepBaseline[SLB\|SPH\|CYL\|CYL_DEG]` | raw `sweep_once` angular+scalar, `array_equal` | ✅ all 4 (`op_ctor` ≥ 1) | ✅ all 4 build the scan cache (`geom` = 1) | spy census |
| A2 | `TestAffineCarveMatvecBaseline[SLB\|SPH\|CYL\|CYL_DEG]` | matvec bulk | ✅ all 4 | ⛔ **NO — `geom` = 0** on all four | spy census |
| A3 | `test_dd_regression[*]`, **14** solve-level snapshots | keff + scalar flux | ✅ **14 of 14** | ✅ the **11** 1-D rows (`geom` = 1); ⛔ the **three** 2-D rows are wavefront, `geom` = 0 | spy census |
| A4 | `test_walk_matvec_baseline[slab\|sphere\|cyl\|cyl_deg\|cart2d]` | frozen matvec | ✅ 4 of 5 (**`cart2d_2g` constructs NO operator**) | ⛔ `geom` = 0 on all 5 | spy census |
| A5 | the 40-file ctor scope | end-to-end behaviour through the migration | ✅ by construction | partial | `[M]` 663 passed / 56.2 s |
| C1 | `cart2d_2g` walk baseline | **CONTROL** — it constructs no operator, so its bit-identity carries ZERO information about step 1 | ⛔ | ⛔ | spy census |

⟹ **Step 1's anchor is A1+A2+A3+A4 — `[M]` 26 of the 27 collected artifacts
construct an operator.** Step 2's cache re-pose anchor is **A1 (4) + the eleven
1-D rows of A3 = 15 artifacts, not 27** (`[M]` the census counts `geom` = 15
builds over 27 tests; 12 tests build zero). Cite that denominator; A2, A4 and C1
are blind to the cache and must not be counted (`vv` anti-#20).

⚠ Step 2's *read re-plumb* anchor is different again, and split by geometry
(F7): `[M]` the SLAB rows are the only ones that execute the scheme's per-cell
kernel dispatch (`residual_kernel_batch` 80 calls; curvilinear = **0**), and
the CURVILINEAR rows are the only ones that execute the closure's
(`cell_contribution` 3 192 – 14 496; slab = **0**). Neither family alone is an
acceptance set.

⭐ **Nothing new needs to be BUILT.** This is the one place P4.9b is easier
than P4.9a: the corpus already covers step 1 universally, and P4.9a's own
`CYL_DEG` additions cover the degenerate arm. The plan's job is to name the
right subset per claim, not to manufacture a fixture.

---

## 4. The mutation battery

**Harness discipline.** Every arm is an **in-process monkeypatch/subclass
installed by a pytest plugin** (`-p mut_p49b`) — no production file is edited
on disk, so a `SIGTERM` at the harness timeout cannot leave the tree mutated.
Strictly stronger than copy-aside + `diff -q` (`lessons` L63h). The plugin
**asserts its own installation** (`raise RuntimeError` unless it binds every
symbol it names) and prints a `sessionfinish` census, so a missing instrument
is visible per row rather than reading as a clean zero.

⚠ `--continue-on-collection-errors`, and count `^ERROR` separately from
`^FAILED`. A mutation that makes production raise kills collection and pytest
reports `FAILED = 0` — the flattering direction.

⛔ **Patch at the DEFINING class, resolved through the MRO — four of the nine
mutated surfaces do NOT live on the concrete class.** `[M]` measured while
building this battery (my own plugin's installation assertion is what caught it,
`vv` #17):

| surface | defined on | kind |
|---|---|---|
| `residual_kernel_batch`, `cell_kernel_batch`, `affine_scan_coefficients`, `residual_kernel_batch_transpose` | `DiamondDifference` | function |
| **`source_emission`, `cell_average`** | **`DiscretizationSchemeBase`** | **staticmethod** |
| `cell_contribution`, `precompute_psi_state` | `MorelMontryAngularSweep` | function |
| **`advance_psi_half`** | **`PoleAngularClosureBase`** | function |
| **`c_out_per_ordinate`** | **`PoleAngularClosureBase`** | property |
| `level_indices` | not on the MRO — an instance attribute | — |

A battery that patches only `DiamondDifference` / `MorelMontryAngularSweep`
binds **5 of 9** and reports a clean, confident, partial zero. Resolve
`next(c for c in cls.__mro__ if name in c.__dict__)`, handle `staticmethod` and
`property` separately from plain functions, and `raise` unless every named
symbol binds.

⭐ **And the design consequence, which is NOT about the battery:
`source_emission` and `cell_average` are STATICMETHODS on the base**, so
`mesh.scheme.source_emission` and `op.spatial_closure.source_emission` resolve
to the *same function object* regardless of the instance. Re-plumbing those two
reads is therefore **value-inert by construction** — no value gate can ever see
it, at any tolerance, on any fixture (`vv` Mode 12 at the dispatch). Their only
possible witness is the structural read-set gate (§2.2). Say so in §2.2's
docstring, or a future audit will count them as covered by the route gate's
numbers.

⚠ **Scope every arm to the files that can redden**, budgeted off the MUTATED
cost (a garbage mutation destroys convergence and blows a baseline-sized
timeout — `[M]` this battery's first attempt at the 40-file + regression scope
overran a 10-minute budget against a 56 s baseline). `[M]` scope costs at
`10314dfa`:

| scope | tests | wall |
|---|---:|---:|
| the 40 ctor-site files | 663 p / 1 d / 19 xf | **56.2 s** |
| `tests/sn/operators` | 1240 p / 1 d / 5 xf | **74.0 s** |
| `tests/sn/sweep` | 911 p / 1 s / 8 d / 31 xf | **282.2 s** |
| `tests/transport` | 566 p / 1 s | **22.3 s** |
| `tests/sn/{solve,regression,architecture}` | 327 p / 2 d / 15 xf | **311.9 s** |

### 4.1 The battery

| arm | mutation | scope | expected | `[M]` pre-carve |
|---|---|---|---|---|
| **M0 control** | `DiamondDifference.residual_kernel_batch` ×1.05 (class-level) | operators + sweep-core | many reds | run at step 0; a small red set means the harness is dead, not that the tree is untested |
| **M1 superset — the phase's proof** | class-level mutation of the OWNER surfaces the walk consumes, split into `m1_scheme` (`residual_kernel_batch`, `source_emission`, `cell_average`) and `m1_closure` (`advance_psi_half`, `cell_contribution`) | frozen corpus, then operators + sweep | ⛔ post-carve red set must be a **SUPERSET** of the pre-carve red set, **per arm** | **`[M]` TAKEN — see §8** |
| **M2 route, scheme** | swap `sn.scheme` for `MutantDD` AFTER pose | §2.1's gate | must red PRE-carve, GREEN post | `[M]` rel **5.000e-02**, `array_equal` False |
| **M3 route, closure** | swap `sn.pole_angular_closure` for `MutantMM` AFTER pose | §2.1's gate | must red PRE-carve, GREEN post | `[M]` rel **4.596e-02 / 5.313e-02 / 1.196e-01** |
| **M3b insufficient-mutation control** | swap the closure, mutate **`cell_contribution` only** | §2.1's gate | ⛔ **must stay GREEN pre-carve** — the arm exists to prove M3's surface choice was necessary | `[M]` `array_equal` **True** on all 3 curvilinear rows |
| **M4 pose reads the wrong hub attribute** | `.pose` swaps its two arguments | §1.2 | reds legs 1+2 of the pose-identity gate | post-carve |
| **M5 pose mints instead of reads** | `.pose` constructs fresh `DiamondDifference()` / `closure_cls(...)` | §1.2 + `test_angular_closure_is_single_object` | reds the identity legs; ⚠ **must NOT move any value gate** — that is the point (it is a silent invariant break) | post-carve |
| **M6 lazy table drift** | the lazy artifact computes `mm_a_in_coeff` as `tau_inv - 1.0` | §2.3(c-i) + A1 | reds `array_equal`, NOT `allclose` | `[M]` the P4.9a measurement: **1–2 ULP** |
| **M7 memo scoping** | force the lazy artifact to rebuild on every access | §2.3(c-ii) | reds the COUNT gate and nothing else | `[M]` the count is **1** today; the regression is worth up to **24.65 %** |
| **M8 rebind staleness** | `rebind_cross_sections` no longer invalidates the Stratum-2 artifact | §2.4 d-4 | reds the "post-rebind answer reflects the new σ" leg | ⛔ `[M]` **no such leg exists today** — the teeth are net-new |
| **M9 `cache.py:273` guard** | delete the converted `raise` | §2.6 | reds the 2-D negative test | ⛔ `[M]` the bare `assert` is inert under `-O` today ⟹ pre-carve the gate cannot exist |
| **M10 diagnostic carve-out control** | mutate `_face_transmission_spectrum`'s probe instance only | §1.2 / any identity gate | ⛔ must NOT red an identity gate | `[M]` F10 — a legitimate second `DiamondDifference` exists in every solve |

### 4.2 The proof obligation

**M1 is the phase's whole point, and it must be read by SET, not by count**
(`vv` #17's identity clause). Post-carve the owner mutation's red set must
CONTAIN the pre-carve set. A same-sized but disjoint set means the re-plumb
reached a *different* object — exactly the failure mode `.pose` exists to
prevent, and the only one no value gate can see.

⚠ **Take the pre-carve M1 red set at step 0, before any edit.** It cannot be
reconstructed afterwards.

⚠ **M5 is the arm that carries the ruling.** The user's no-guard position rests
on "the hub creates and passes, so disagreement is unspellable". `.pose` minting
fresh objects makes it spellable again while every value gate stays green — so
the identity legs are the *entire* enforcement of the phase's safety argument.
Record M5's red set explicitly.

---

## 5. Gate order and costs

`[M]` all timings at `10314dfa`, `.venv/bin/python -O -m pytest -p no:randomly
-q -m "not slow"`, SERIAL. Baseline for the whole fast set:
`9836 / 0, 22 sk / 227 des / 70 xf`, ~64 min (plan header).

**Step 0 — baselines, before any edit.** (a) Re-take the five scope numbers in
§4's table on the branch tip. (b) **The M1 pre-carve red sets are already taken
— §8.** Widen them to `tests/sn/operators` + `tests/sn/sweep` if the battery is
to claim anything outside the frozen corpus. (c)
`[M]` `mcp__nexus__dead_references` = **0 dead / 52 checked** — the step-3
denominator. (d) `git status --porcelain` (this tree carries uncommitted
state; a concurrent writer invalidates every number — `lessons` §2).

**Step 1a — the ctor + `.pose` + the production route.**
`streaming.py` ctor gains two required fields; `.pose` classmethod;
`coupled_system.py:441` routes through it; `is_adjointable` and
`streaming.py:992` dissolve onto own fields.
*Gates:* §1.1 arity witness (with `match=`), §1.2 pose-identity gate.
*Cost:* `tests/sn/operators/test_streaming_operator.py` alone first, then
`tests/sn/operators` **74 s**.

**Step 1b — the 135-site migration.** Order per §1.4. Land in ONE commit with
1a: `[M]` the production call site and the test call sites are one §6b set, and
a partial migration is a collection kill, not a test failure.
*Gates:* the 40-file scope **56.2 s** → `663 passed / 19 xfailed`; then the
frozen corpus (A1–A4) `array_equal`; then `tests/sn/sweep` **282 s**.

**Step 2a — `default_for` + the representation fields (Q2's answer).**
*Gates:* §2.5's census re-run to zero; `tests/sn/sweep` + `tests/sn/operators`.

**Step 2b — the 43 (ii) + the class-(iv) reads re-plumb onto the operator's
fields.**
*Gates:* ⭐ **§2.1's route gate flips RED → GREEN — this is the step's
acceptance**; §2.2's read-set allowlist; §2.7's expert-seam characterization
test (land it HERE, not at step 1 — the interval hazard); A1 + the ten 1-D A3
rows `array_equal`; `tests/sn/regression` (inside the **311.9 s**
solve+regression+architecture scope).

**Step 2c — the lazy table + the three memo retirements + `cache.py:273`.**
*Gates:* §2.3(c-i) equivalence (this commit only), §2.3(c-ii) COUNT, §2.4's
four legs including the re-posed rebind contract, §2.6's `-O` raise witness.
⛔ Delete the eager build in the SAME commit as (c-i), or (c-i) goes
tautological.

**Step 3 — battery + docs + audit.** M0–M10 with the table recorded, not a
boolean. Then: the ~127 prose mentions + the ~13 `StreamingOperator(sn_mesh)`
literals; the present-tense-false `default_angular_closure_class` docstring
(`[M]` verified rather than relayed — `orpheus/sn/angular/closure.py:2261-2263`
reads *"The factory dispatch (instantiation with ``sn_mesh``) is the caller's
job"*, while the shipped contract is `closure_cls(angular, pairing)`; the
recount's claim reproduces, `plan-authoring` §4's VERIFY move);
`mcp__nexus__dead_references` back to **0 dead / 52 checked**; `sphinx -W`;
`pyright`.

**Full fast gate** before merge (~64 min; budget ≥ 90 min per the reference
memory).

---

## 6. Open rulings — BLOCKING

**Q1 ⛔⛔ (blocking, F2). What is the lazily-resolved table's memo LIFETIME?**
`[M]` the operator is constructed **6–10 times per solve**, so a
`cached_property` on the operator means 6–10 table builds where today there is
**1** — up to **24.65 %** of a slab solve at `GL16`/nx=200. The ruling
("performance welds resolve lazily, as close to the strategy as possible") is
satisfied by any of: (a) keep a mesh/geometry-keyed memo and make it lazy only
(retires the eager build, keeps the count at 1); (b) memoise on the closure
pair; (c) accept the rebuild and record the cost. Whichever is chosen, §2.3's
COUNT gate must assert the number. **A recommendation is not mine to make —
but the plan cannot gate the artifact until the number is ruled.**

**Q2 ⛔ (blocking, F3). Do the 58 direct representation constructions migrate?**
Adding two required fields to `_LossRepresentation.__init__` breaks all 58 (11
test files) at collection. Options: migrate them (step 2's §6b set becomes 79
sites / 14 files); give the representations their own `.pose`-style
classmethod; or keyword-with-mesh-read defaults (⛔ re-instates the read the
phase removes).

**Q3 (blocking for §2.1's spec). Does `sweep_once` grow an `operator=`
parameter?** `[M]` it re-poses internally, so it cannot drive any route gate.
Either the gate builds `(L + C).solve` itself (available today, no production
change) or the helper takes the operator. The first is recommended for the gate;
the second is a broader question about whether test helpers should re-pose at
all.

**Q4 (blocking for §2.2's allowlist). Are `is_affine_scannable` /
`supports_curvilinear` operator-side or hub-side reads?** They are *strategy
selection*, consumed by `default_for`/`supports` before the walk exists. If
`default_for` receives the closure pair (the design's step 2), they become
operator-side; if selection stays a mesh query, they stay hub-side. The
read-set gate's allowlist is exactly this ruling, so it cannot be written
without it.

**Q5 (F11). Naming: `spatial_closure` + `angular_closure`, or neither?**
`[M]` all three names are free tree-wide. Shipping `spatial_closure` alone
leaves three vocabularies for two slots.

**Q6 (F4). Does `_pole_mirror_cache` retire with the other two memos?**
Same idiom, same tier, unnamed by the design.

**Q7. Where does the M8 anti-twin-style permanent gate live for this phase?**
P4.9a made its done-when permanent as an AST/text gate. P4.9b's analogue is
"no production module outside `sn/mesh/` reads `mesh.scheme.<per-cell
attribute>`" — i.e. §2.2's allowlist promoted from a runtime read-set gate to a
static AST gate. Recommend both: the AST gate is cheap and total, the runtime
gate catches the dynamic route.

---

## 7. What this plan does NOT cover

* The end-state cross-method ctor `(domain, codomain, spatial, angular)` with
  `sn_mesh` out — recorded, rides O-3/CS5.
* Dissolving `isinstance(closure, MorelMontryAngularSweep)`
  (`loss_rep:4213`/`:4630`). F6 shows it constrains how the closure can be
  instrumented; retiring it is O-3's.
* The `SNMesh`-the-misnomer rename — filed as its own issue by ruling.
* `#392` (`face_areas` shim), `#398` (ctor admissibility), `#402` (LD
  moment-tailed residual) — adjacent, must not conflict; `[M]` F5's partition
  leaves `moment_axis` / `spatial_basis_per_axis` reachable exactly as today,
  which is #402's precondition.
* Whether the space-side induction should read an induced PRODUCT rather than
  the generator — the design's dead D3, kept dead by the ruling.

---

## 8. `[M]` PRE-CARVE battery record — the M1 denominator, taken 2026-08-28 at `10314dfa`

Plugin: in-process monkeypatch (`mut_p49b.py`, session scratchpad),
MRO-resolved, asserting its own installation (`RuntimeError` unless every named
symbol binds; the printed `BOUND n` is grepped into each result line).
Scope: `tests/sn/regression` + `tests/sn/sweep/core/test_affine_carve_baseline.py`
= **27 collected**, `-m "not slow"`, `--continue-on-collection-errors`,
`^FAILED` and `^ERROR` counted separately (zero `ERROR` in every arm).

| arm | mutated (×1.05) | bound | wall | reds |
|---|---|---:|---:|---:|
| baseline | — | 0 | 48 s | **0** |
| **`m1_scheme`** | `DiamondDifference.residual_kernel_batch` + `DiscretizationSchemeBase.{source_emission, cell_average}` | 3 | 122 s | **20** |
| **`m1_closure`** | `PoleAngularClosureBase.advance_psi_half` + `MorelMontryAngularSweep.cell_contribution` | 2 | 60 s | **16** |
| union | | | | **26 of 27** |

**The red-set PARTITION — an independent corroboration of F7's disjoint
activation:**

| bucket | n | rows |
|---|---:|---|
| **scheme-only** | **10** | `test_dd_regression[{slab_2g_3reg, slab_2g_homogeneous, slab_2g_p1_aniso, slab_fixed_source, 2d_2g_LS4_dd_8x4_het_si, 2d_2g_p1_aniso_dd_8x4_het_si}]`, `test_walk_matvec_baseline[{slab_2g, cart2d_2g}]`, `TestAffineCarveMatvecBaseline[SLB]`, `TestAffineCarveSweepBaseline[SLB]` |
| **closure-only** | **6** | `test_walk_matvec_baseline[{cyl_2g, cyl_deg_2g, sphere_2g}]`, `TestAffineCarveMatvecBaseline[{CYL, CYL_DEG, SPH}]` |
| **both** | **10** | the curvilinear solve-level regressions + the curvilinear sweep baselines |
| **never red** | **1** | `test_dd_regression[2d_1g_LS4_dd_15x15]` |

⚠ **The survivor is my ARM's composition, not a blind gate.** `[M]` that row is
2-D wavefront and its per-cell surface is `cell_kernel_batch`
(**406 564** calls, activation census F7), which `m1_scheme` does **not**
mutate. Adding `DiamondDifference.cell_kernel_batch` to the arm is the fix; I
am recording the gap rather than banking a 26/27 as though it were a coverage
measurement (`vv` #17's granularity trap).

**How to read this post-carve.** For each arm separately, the post-carve red set
must CONTAIN the pre-carve one. A same-sized but disjoint set means the
re-plumb reached a *different* object — the one failure mode `.pose` exists to
prevent and the only one no value gate can see (`vv` #17's identity clause).
⛔ Do NOT compare the union: the two arms partition the corpus by geometry, so
a union comparison would hide a scheme-side regression behind closure-side reds.

---

## 9. `[M]` POST-CARVE battery record — run 2026-08-28 at `d14dd545` (main agent)

Same plugin (`mut_p49b.py`), same scope as §8 (27 collected), plus the
step-3 arms in `mut_p49b_step3.py` (both in the session scratchpad;
in-process, self-asserting, `--continue-on-collection-errors`, zero `ERROR`
in every arm).

| arm | bound | reds | verdict |
|---|---:|---:|---|
| baseline (none) | 0 | **0 / 27** | clean tree |
| `m1_scheme` | 3 | **20** | ⭐ **EQUALITY with §8's set** — same 20 rows (scheme-only 10 ∪ both 10); superset obligation met as identity, per arm |
| `m1_closure` | 2 | **16** | ⭐ **EQUALITY with §8's set** — same 16 rows (closure-only 6 ∪ both 10) |
| M2/M3/M3b | — | — | LIVE in the keystone gate: pre-carve reds verified at its strict-xfail landing (`700f00e4`); post-carve all four rows green with activation legs moving |
| M4 (pose swaps args) | 1 | **1** | exactly the pose-identity gate |
| M5 (pose MINTS — the ruling-carrier) | 1 | **5** | pose-identity + is-single-object **leg 0** (the §1.3 loop-closer — it caught it) + the keystone's three closure-row ACTIVATION legs. ⭐ Every VALUE assertion stayed green — the invariant break is silent to values, structural legs are its ONLY catchers, which is the no-guard ruling's entire enforcement, measured |
| M7 (intern bypassed) | 1 | **2** | the COUNT gate + the invariance witness — the two memo instruments, nothing else |
| M8 (rebind stale) | 1 | **1** | exactly the re-posed rebind witness — the net-new staleness teeth bite |
| M9 (guard deleted) | 1 | **1** | exactly the −O raise witness |
| M10 (diagnostic probe) | — | N/A | no identity gate references the type-keyed `_face_transmission_spectrum` probe instance (F10's carve-out is by construction; the gauge suite is green in every wide gate) |
| M6 (lazy-table respelling) | — | VACUOUS BY DESIGN | the strategy's table COPIES the closure's P4.9a-minted constants — there is no second derivation to respell; the P4.9a weld gate (`TestMintedScanConstants`) owns that surface |

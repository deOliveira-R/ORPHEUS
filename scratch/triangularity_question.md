# Is `L + C` ALWAYS triangular in ORPHEUS? — evidence dossier

Branch `refactor/operator-naming-honesty` @ `8d996f53`. Line numbers current at
this HEAD; the *structural* claims are durable, the `file:line` map is not.

**Short answer up front:** No — and ORPHEUS already knows it, says so in the
theory corpus, and certifies it per case. `docs/theory/foundations/path_integral.rst:1078`
states the position in one sentence:

> **Acyclicity is a property of the (mesh, closure, boundary) triple**, and
> ORPHEUS certifies it *per case* rather than assuming it.

What makes `(L+C)` *appear* always-triangular in the code is an **architectural
choice plus a scope restriction**, not a theorem about the operator: the one
cycle source that exists (`B`) was structurally evicted from `(L+C)` into a
sibling operator, and the one boundary law ORPHEUS cannot handle acyclically
(periodic) is refused at mesh construction.

---

## 1. The forward/inverse split as it stands

### Where triangularity is ASSERTED (prose only)

- `orpheus/sn/operators/sweep_operator.py:86` — "The inverse operator
  :math:`A^{-1}` of a **schedule-triangular** forward ``A``".
- `orpheus/sn/operators/sweep_operator.py:56-57` — the ctor-guard rationale,
  verbatim:

  > its ctor guard is the ``SweepInvertible`` TYPE itself
  > (schedule-triangularity is what makes the forward sweep-invertible —
  > **no value check needed**).

- `orpheus/sn/operators/sweep_operator.py:77-80` — `SweepInvertible =
  StreamingCollisionOperator | ScheduledInvertibleOperator`, "one wrapper for the
  whole schedule-triangular family".

### Where it is CHECKED at runtime — nowhere

`StreamingCollisionOperator.__init__` (`orpheus/sn/operators/streaming.py:564-603`)
runs exactly four guards, none about ordering:

1. `isinstance(streaming, StreamingOperator)` — `:569`
2. `isinstance(diagonal, MultiplicationOperator)` — `:574`
3. mesh-identity `streaming.sn_mesh is diagonal.coefficient.mesh` — `:585`
4. `np.all(diagonal.coefficient.values > 0)` — `:593`. This is
   *non-singularity of the per-cell WDD divide* ("The sweep would emit NaN at
   those cells"), NOT a sweep-ordering check.

`is_invertible` is a hard-coded constant (`streaming.py:608-615`):

```python
@property
def is_invertible(self) -> bool:
    # (L+C) is sweep-invertible: the WDD forward-substitution sweep IS its
    # inverse operator ... This is the SOLE invertible OperatorSum
    return True
```

`SweepOperator` has **no `__init__` at all** — the back-half comes from
`InverseWrapMixin`, and the "guard" is the static type parameter
`InverseWrapMixin["SweepInvertible"]` (`sweep_operator.py:83-85`), which is a
`TYPE_CHECKING`-only `Union` (`:69-80`) — i.e. erased at runtime.

> **Verdict for Q1: triangularity is asserted by CONSTRUCTION and by STATIC
> TYPE. There is no runtime guard anywhere that a sweep ordering exists.**

Two near-misses worth naming, because they look like the guard and are not:

- **`default_for`** (`orpheus/sn/loss_representation/__init__.py:2647-2690`) raises
  `IncompatibleRepresentation` — but on **(scheme × geometry)** pairing (e.g.
  Linear-Discontinuous on a sphere), not on ordering. It also fires **lazily**,
  on the first read of the `cached_property`
  `StreamingOperator.loss_representation` (`streaming.py:384-404`), so a
  `(L+C)` object can be constructed successfully and only fail at first use.
- **`CoupledOperator._triangular_orientation()`**
  (`orpheus/numerics/coupled_system.py:850-888`) IS a real structural check —
  it inspects the block grid, returns `"upper"|"lower"|None`, and
  `CoupledSubstitutionOperator.__init__` (`:1248`) raises when the grid is not
  block-triangular. But it operates on the **2×2 block grid** (the ψ½ system),
  not on the cell DAG inside `(L+C)`.

---

## 2. The sweep dependency graph — index arithmetic, not a topological sort

`orpheus/sn/loss_representation/sweep_graph.py`.

**Verdict: closed-form index arithmetic derived from the octant sign pattern.
Nothing traverses a graph; nothing detects a cycle; acyclicity is assumed.**

`SweepDependencyGraph.from_cartesian` (`sweep_graph.py:512-610`) is the whole
construction. The load-bearing lines:

```python
face_in  = tuple(0 if s >= 0 else 1 for s in signs)      # :572
face_out = tuple(1 if s >= 0 else 0 for s in signs)      # :573
axis_map = [np.arange(n) if s >= 0 else np.arange(n)[::-1]
            for n, s in zip(shape, signs)]               # :577-580
local     = np.indices(shape).reshape(d, -1)             # :586
level_of  = local.sum(axis=0)                            # :587
n_levels  = sum(n - 1 for n in shape) + 1                # :588
```

`signs` is the octant label, i.e. `sign(mu_axis)` per axis (projected in-plane
at the single site `sweep_schedule._octant_sweep`, `sweep_schedule.py:217-238`).
Cells on level `k` are those with `Σ_a local_a == k` — the **anti-hyperplane of
the index lattice**. The class docstring is explicit that the absence of a
topological sort is deliberate (`sweep_graph.py:381-386`):

> The closed-form precompute on a Cartesian grid lives entirely inside
> `from_cartesian` and never appears in the sweep loop. This is structural, not
> hand-rolled — the "library version" (a generic topological-sort over an
> explicit DAG) would be over-engineering for a regular pattern that collapses
> to ~5 lines of arange + diagonal extraction.

So "DAG" is a *name for the ordering*, not a walked data structure. `walk_full`
(`:721-764`) and `walk_windowed` (`:766-814`) iterate `self.levels` in storage
order. No visited-set, no back-edge test, no `networkx`.

**The only rejection in the builder** is the degenerate all-zero octant
(`sweep_graph.py:554-559`) — an ordinate with no in-plane streaming has no
ordering, and the walk short-circuits it with `Q/Σ_t`. That is a *degenerate*
case, not a *cyclic* one.

### The declared-but-unimplemented cycle invariant

`sweep_graph.py:43-49` cites the Grand Report §15A.2 invariant set the module's
tests "pin": `assert_upwind_orientation`, `assert_topologically_sorted`,
`assert_face_pairing_consistent`, `assert_boundary_trace_classification`,
**`assert_cycles_are_declared`**.

Repo-wide grep: **`assert_cycles_are_declared` occurs ONLY in that docstring.**
`tests/sn/sweep/core/test_sweep_graph.py:7-19` lists the other four (plus
`assert_cell_coverage`, `apply_invariants`) and silently drops the cycles one;
`grep -n cycle` on that test file returns nothing.

> **The cycle invariant is named in the design and never realized.** Consistent
> with §3 — on a structured Cartesian lattice there is nothing for it to catch.

### The 1-D path is a SCAN, not a wavefront

`CumprodScan` (`loss_representation/__init__.py:1330`) with executor
`_OneDimScanWalk` (`:2769`) evaluates the first-order linear recurrence in
**closed form** (Blelloch 1990 §1.5 — `cumprod`/`cumsum`; cited at `:2706`), not
by walking anything. `_WalkLeg` (`:2724-2749`) factorizes the 1-D node set into
(μ-level × direction) chains via the sign trichotomy `_MU_DIRECTION_EPS = 1e-15`
(`:2721`). Ordering again from the sign of μ, arithmetic, never a traversal.

(This is the durable **scan-vs-wavefront** distinction: 1-D and curvilinear do
NOT ride `SweepDependencyGraph`. Any "unify the sweeps" proposal must respect it.)

---

## 3. Dimensionality

`LOSS_REPRESENTATIONS` (`loss_representation/__init__.py:2639-2644`) is the
selection policy in priority order; `default_for` returns the first whose
`supports` admits the mesh:

| Representation | Applies to | Ordering mechanism |
|---|---|---|
| `CumprodScan` (`:1330`) | every 1-D mesh (slab, sphere, cylinder) | parallel-prefix scan over μ-sign legs |
| `ScanMarch` (`:2166`) | multi-D Cartesian — production default (S6.9 Fork-B2) | row-march |
| `MovingFrontierWindow` (`:1508`) | d=2 selectable peer | rolling (d−1)-frontier anti-hyperplane wavefront |
| `FullFieldWavefront` (`:1788`) | **any-d Cartesian** — "the never-stuck any-d spine" | full-cochain anti-hyperplane wavefront |

`default_for`'s docstring answers d≥3 directly (`:2669-2671`):

> a d≥3 Cartesian mesh — axis-native via ``SNMesh.from_axes`` since C5.5 (#225)
> — → ``FullFieldWavefront``, the never-stuck any-d spine.

`from_cartesian` is written dimension-generic — "``d = len(shape) ∈ {1, 2, 3}``"
(`sweep_graph.py:521`), d=3 levels are "a 2-D lattice slice of the ``i+j+k = ℓ``
hyperplane" (`:537`), and `_build_frontier_plan` carries an explicit `d ≥ 3`
fancy-simplex-index branch (`:653-662`). **d≥3 is BUILT, not seam-only** — #227
is performance/scale work (and the "future direction" note at `sweep_graph.py:61-66`
proposes `spsolve_triangular` for very large 3-D), not existence.

**Is the per-octant ordering provably acyclic on a structured Cartesian mesh in
d dims — and does the code rely on that?** Yes to both. The theorem is stated in
the corpus (`path_integral.rst:1070-1077`):

> For a fixed ordinate on a structured Cartesian mesh with a Cartesian closure,
> the cell-dependency digraph is acyclic by a lattice product-order argument:
> `sign(Ω_x)·i + sign(Ω_y)·j` is a strict potential that every dependency edge
> increases. **Note what the proof uses: the mesh's lattice order — not the
> characteristics. It is a mesh theorem, and it is exactly as strong as its
> hypotheses.**

The code relies on it by **constructing the order from that very potential**
(`level_of = local.sum(axis=0)`, `sweep_graph.py:587`) rather than discovering
it — exactly the "generate from sign(mu_axis) per axis" hypothesis in the
question. The check that the assembled *matrix* actually respects the order is
empirical (§6), not structural.

---

## 4. Where non-triangularity ACTUALLY enters

### The architectural headline: `B` is evicted, so `(L+C)` never contains the cycle

`creates_sweep_cycle: ClassVar[bool]` is declared on the ABC
(`orpheus/geometry/boundary/_base.py:146`, default `False`) and set `True` on
`ReflectiveBoundary` (`reflective.py:102`) and `PeriodicBoundary`
(`periodic.py:58`); `False` on vacuum/white/albedo/prescribed (pinned by
`tests/geometry/test_bc_universal_invariants.py:387-443`).

**But the flag is WRITE-ONLY.** Every occurrence in `orpheus/` is a declaration
or a docstring; the only *attribute reads* in the entire repo are the eight
assertions in `tests/geometry/test_bc_universal_invariants.py` and
`tests/geometry/test_boundary_trace_law.py`. `docs/theory/foundations/boundary_conditions.rst:1842-1846`
describes it as "the structural signal that lets **the §15A.2 sweep-cycle
detector** identify these BCs without inspecting the realization" — that detector
does not exist. The flag cannot gate anything, because ORPHEUS solved the problem
structurally instead:

The bare sweep reads `ψ.boundary.inflow` as GIVEN and writes `ψ.boundary.outflow`;
it does NOT re-apply `R·G` internally (`streaming.py:51-70`, and the "KEYSTONE
DELETED" comment at `loss_representation/__init__.py:3269-3275`). The reflective
coupling lives in the sibling `−B` (`SNBoundaryOperator`), and the outer SI/Krylov
loop drives `ψ.inflow − B·ψ.outflow → 0`. So the cycle is real, it is just **not
inside the operator whose name is at stake**.

Partial recovery is available and *named*: `SweepSchedule.gauss_seidel`
(`sweep_schedule.py:111-168`) + `SNBoundaryOperator.split`
(`orpheus/sn/operators/boundary.py:431-448`) partition `B = B_lower + B_upper`,
fold `B_lower` into the reified `M = (L+C−B_lower)`
(`ScheduledInvertibleOperator`), and lag `B_upper` — "the **cyclic back-edges**
plus every row of a never-reflected face" (`boundary.py:439-440`). A
both-faces-reflective axis is explicitly called a **2-cycle** giving only partial
one-pass G-S (`sweep_schedule.py:22-25`; `cartesian_multid.rst:3796-3800`).
`_reflective_faces` (`sweep_schedule.py:257-268`) excludes white — "couples ALL
ordinates on the face ⟹ the octant-order G-S degenerates to Jacobi".

**And the one law ORPHEUS cannot handle acyclically is refused at the door.**
`SNMesh.BOUNDARY_OPERATOR_REGISTRY` (`orpheus/sn/mesh/augmented_mesh.py:171-187`)
admits exactly **two** kinds:

```python
BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]] = {
    "vacuum": VacuumInflow,
    "reflective": ReflectiveBoundary,
}
```

with the comment:

> The 5 other kinds the realizer handles today (``white``, ``periodic``,
> ``albedo``, ``prescribed_inflow``, ``mixed``) are NOT registered here —
> **adding them requires SN-sweep-side wiring (sweep cycles for periodic, etc.)**

`resolve_boundary_conditions` raises `ValueError(... does not support boundary
condition '<kind>' on face '<face>')` for anything unregistered
(`orpheus/transport/method.py:289-298`). So `BC("periodic")` on an `SNMesh` is a
**hard construction-time failure**. This is the single strongest *enforced*
guard in the whole story — and note it is a guard on the MESH, not on `(L+C)`.

### (a) Curvilinear angular redistribution — triangular only in the AUGMENTED sense

Yes: the sphere is triangular **only after the space × angle × ψ½-seed space is
augmented**, and that is exactly what the `[[LC, Seeding], [None, march]]` grid
statement means.

- The Morel–Montry μ-derivative couples ordinates within an angular level. The
  **cylinder** is fine as-is — non-carrying per R12a; `test_assembly_mode.py:40-43`
  records "The cylinder control shows exact triangularity ... **NOT** α-dome
  'telescoping the seed away', corrected in #280 Phase 2.5b".
- The **sphere** was genuinely non-triangular: `test_assembly_mode.py:521-524`,
  "The #280 Phase 2.5d fix retired the lagged Morel–Montry ψ½ pole seed (**a
  walk-order BACK edge that made the spherical operator non-triangular in sweep
  order**)". The falsification is on record with a number —
  `path_integral.rst:1086-1092`: "a cold residual of `5.2 × 10^5` in the
  operator-level probes, with the seed iteration diverging geometrically".
- Route (a) (#282) made ψ½ a **first-class STATE block**, so the augmented matrix
  is `[[A_ss, 0], [A_bs, A_bb]]` — triangular in
  `[seed⁻ march, seed⁺ march, ordinates↑μ]` order
  (`test_assembly_mode.py:660-673`).
- **So yes — the `[[LC, Seeding], [None, march]]` statement IS the
  augmented-triangularity claim.** `orpheus/sn/coupled_system.py:87-94`:
  "`M = [[L+C, +Seeding], [None, march]]` — since step 5 an **HONEST
  upper-triangular** `CoupledOperator` grid ... `solve` is the numerics block
  back-substitution (`ψ_B = march⁻¹ q_B` — System B first, exactly the
  curvilinear sweep order — then `ψ_A = LC⁻¹(q_A − Seeding·ψ_B)`)".
- Inside System B, `A_BB` is itself triangular *in radius*:
  `orpheus/sn/operators/radial_characteristic.py:44-49` — "Because the recurrence
  is triangular in radius (**the #284 forward-substitution certificate**), the
  two-leg Carlson march ... is the EXACT direct inverse `A_BB⁻¹` — no iteration."
- **This is the case where the composite guard is real** (§1): the grid's
  triangularity IS checked, by `CoupledOperator._triangular_orientation()`.

### (b) LD / DG within-cell moment coupling — BLOCK-triangular at cell granularity

Confirmed, and pinned. `test_g2_ld_block_triangular_and_lapack_solve_equals_sweep`
(`tests/sn/sweep/test_assembly_mode.py:436-498`):

> the assembled block is **BLOCK-lower-triangular in walk order** (exact — the
> **moment blocks are dense within a cell**, so the structural zero lives ABOVE
> the cell blocks)

The mask is `kron(triu(1,1), ones(cm×cm))` with `cm = 2**d`
(`_block_upper_mask`, `:385-389`), and the DOF permutation is cell-major /
moment-minor (`_dof_order`, `:392-395`; same convention documented on
`ordinate_walk_order`, `orpheus/sn/loss_representation/assembly.py:168-173`).
Note the LD arm's solve leg uses **`lu_solve`, not `solve_triangular`** (`:490`)
— an honest admission that the dense diagonal blocks need a real factorization.

### (c) The r=0 pole mirror — a specular reflection kept INSIDE the walk

The certification #300 cites. Two places:

**The doctrine** (`docs/theory/foundations/path_integral.rst:1093-1109`):

> **Reflection does not automatically force extraction.** The folklore "a
> reflective boundary creates a cycle, so the reflected coupling must be
> extracted from the sweep" is *false as stated*: ORPHEUS keeps the curvilinear
> `r = 0` pole mirror — a specular reflective coupling — **inside the walk** as a
> forward edge, because the sweep order visits `μ < 0` before `μ > 0` and the
> mirror feeds information only downstream (`orpheus/sn/loss_representation`,
> certified lower-triangular). A *single* reflecting face is acyclic; a cycle
> needs a *closed loop* — e.g. both faces of a slab reflecting. The honest
> extraction criterion is therefore not a boolean on the boundary type but an
> **SCC decomposition** of the (face, ordinate) dependency digraph ... (The S_N
> grand report *proposes* names for these components — `SweepStrongComponent`,
> `ReflectiveSweepCycle` — as a design direction; **ORPHEUS does not yet ship
> them**.)

**The code** (`orpheus/sn/loss_representation/__init__.py:3278-3289`) — the
inward sweep runs first, then:

```python
if curvature != "cartesian":
    # Carlson coupled-pole spatial seed (ERR-058 a, Issue #195):
    # at r = 0 the outward characteristic is the CONTINUATION of
    # the inward one — ψ(0, +μ) = ψ(0, −μ) — so the +1 sweep's
    # pole-face seed is the −1 sweep's pole-face outflow at the
    # mirror ordinate (already computed above: data, propagated
    # from the outer boundary, lower-triangular).
    pole_face_seed = outflow_at_inner.T[quad.reflection_index("x")]
```

The parenthetical "(already computed above: data, propagated from the outer
boundary, **lower-triangular**)" is the certification-in-comment. The *executable*
certification is `test_282_augmented_walk_order_is_triangular`
(`test_assembly_mode.py:650-701`) on `[SPHERICAL, CYLINDRICAL]` with
`bc_left=BC("reflective")` — i.e. the pole-mirror fixture — asserting
`triu(permuted, k=1) == 0` exactly.

Note #300's line reference (`__init__.py:2840-2846,3654-3660`) has **drifted**;
the live site is `:3278-3289` (forward) with the adjoint mirror at `:3639-3641`.

### Any other source? No.

Sweeping the space of what could break triangularity in ORPHEUS today:

| Candidate | Status |
|---|---|
| reflective `B` | real cycle source — **evicted** to sibling `−B`; partial in-sweep fold via G-S schedule |
| periodic `B` | real cycle source — **refused at `SNMesh` construction** |
| white / albedo `B` | rank-1 per face; no octant order to respect (excluded from G-S, degenerates to Jacobi) |
| curvilinear μ-redistribution (sphere) | was a back edge; **fixed by augmentation** (#282 route (a)) |
| curvilinear (cylinder) | never was — non-carrying (R12a) |
| LD/DG moment coupling | dense *within* a cell → **block**-triangular, not a cycle |
| r=0 pole mirror | forward edge, certified lower-triangular |
| unstructured mesh | **does not exist** (§5) |
| grazing / pure-z ordinates | no ordering needed — `Q/Σ_t` short-circuit (`sweep_graph.py:554-559`) |

---

## 5. Unstructured meshes

**ORPHEUS has none.** `orpheus/geometry/mesh.py` defines exactly `RegionMesh`
(`:134`), `Mesh1D` (`:211`), `Mesh2D` (`:542`); `orpheus/geometry/structured_geometry.py`
defines `Region` (`:142`) and `StructuredGeometry` (`:211`). There is no tet, no
polyhedral, no connectivity-array mesh anywhere in the package. The canonical
mesh surface is axis-native (`SNMesh.from_axes`, `augmented_mesh.py:607`) — a
tensor-product lattice by construction.

**The seam is documented — as prose, not as a type or a `NotImplementedError`.**

- `docs/theory/methods/sn/placement.rst:127-130`:

  > The third [price] is the sweep schedule itself: the triangularity that makes
  > the solve cheap is a *theorem about the mesh*
  > (:ref:`path-integral-method-map`), and **unstructured or cyclically-coupled
  > configurations must first be decomposed before they sweep.**

- `docs/theory/methods/sn/cartesian_multid.rst:4049-4053` names the literature
  (`:cite:`Pautz2002``; Adams & Larsen 2002 §VI on parallel sweeps) in the
  context of the cyclic back-edge / partial one-pass G-S discussion.
  `docs/refs.bib:531-533` = Pautz, "An Algorithm for Parallel S_n Sweeps on
  Unstructured Meshes" — **the Pautz/Plimpton problem is cited but not solved.**
- `docs/theory/methods/sn/acceleration.rst:922` and the DSA literature memos
  (`.claude/plans/dsa_literature_memo.md:868,936`;
  `.claude/plans/phase_i_survey_adams_larsen_2002.md:340-341,730`) record
  unstructured-grid acceleration as an *open problem in the canon*.
- `.claude/plans/neutron_transport_grand_report_v3.md:2148-2160` is the design
  sketch: "If acyclic, a topological order gives a sweep. If reflective or
  periodic boundaries create cycles, the graph decomposes into strongly connected
  components" — proposing `SweepStrongComponent` / `CyclicSweepBlock` /
  `ReflectiveSweepCycle`. `path_integral.rst:1109` states flatly: **"ORPHEUS does
  not yet ship them."** Same at `.claude/plans/neutron_transport_grand_report_v3.md:1254`
  (`space.assert_sweep_graph_acyclic_for_direction(omega)` — proposed, not built).
- `.claude/plans/sn_operator_realization_and_posing_plan.md:779-780, 851-852,
  1170-1171, 1283` is the most recent (uncommitted, on this branch) planning
  thread: "acyclic trace dependencies can be folded into a schedule; cyclic trace
  dependencies form SCCs that require a coupled solve or lag." And at `:854`:
  "The existing `creates_sweep_cycle` boolean is useful as a **conservative
  hint**, but ..." — i.e. the plan itself already treats the flag as insufficient.

**Gap:** there is no code artifact — no Protocol, no `MeshTopology` seam, no
raising stub — marking where an unstructured mesh would enter. The acknowledgment
is 100% documentation + plans.

---

## 6. The certificate (#284) — what it actually certifies

**#284 itself is not a triangularity issue.** `gh issue view 284` (CLOSED,
`module:sn`, `type:improvement`) is titled *"decide the sweep-inverse contract on
non-source rhs (outflow-row projection) before the first M⁻¹A consumer"*. Its
content: every sweep-realized inverse is exact only on the **source subspace**
`{y : y.outflow-rows = 0}`, because the sweep's shed re-derives the
outflow-definition rows; a residual `r = b − Ax` carries non-zero outflow rows and
was being silently projected. Three options were posed (documented projection /
raise-guard / consistent lift).

It was closed by the **consistent lift** — the ERR-071 one-site restore now
visible in `streaming.py:952-979` (forward) and `:1081-1102` (transpose), with the
motivation recorded inline: the dropped term "made the preconditioner
`M = (I + 𝒞)∘(L+C)⁻¹` SINGULAR on the outflow-trace subspace (measured
‖Mq‖/‖q‖ = 1e-15 — GMRES stalled at an O(1) true residual)".

**The name "forward-substitution certificate" attached to #284 because the
DISCHARGE is a triangularity gate.** The tests live in
`tests/sn/sweep/test_assembly_mode.py`:

| Test | Line | Level | What it certifies |
|---|---|---|---|
| `test_g2_walk_order_triangularity_is_exact` | `:187-205` | `l0` | `triu(PᵀMP, 1) == 0` **exactly** (zero tolerance) for every ordinate × group, on `slab` + `cartesian_2d` |
| `test_g2_lapack_forward_substitution_equals_sweep` | `:211-238` | `l2` | `scipy.linalg.solve_triangular(PᵀMP, Pᵀq) ≡ (L+C).solve(q)` — LAPACK `dtrtrs` is a structurally-independent forward substitution |
| `test_g2_ld_block_triangular_and_lapack_solve_equals_sweep` | `:438-498` | `l2` | LD arm: block-lower-triangular + `lu_solve` ≡ sweep |
| `test_282_augmented_walk_order_is_triangular` | `:655-701` | `foundation` | sphere + cylinder: the **augmented** (ψ½ ⊕ bulk) matrix is exactly block-lower-triangular |
| `test_282_teeth_coupling_direction_swap_reds` | `:705-754` | `foundation` | mutation test — inject a back edge, confirm `triu != 0` |
| `test_teeth_shared_kernel_sign_flip_moves_all_three_modes` | `:288-356` | `foundation` | one-source teeth: sweep, matvec, assembly move together |
| `test_curvilinear_refuses_the_cartesian_walk` | `:502-516` | `foundation` | the Cartesian assembler refuses a sphere (`AttributeError, "Cartesian-only"`) |

Module docstring `:8-16`, verbatim:

> **G2** ... LAPACK ``dtrtrs`` is a structurally-INDEPENDENT forward substitution
> vs the ORPHEUS sweep recurrence, which EARNS this gate its L2 cross-check
> status and **discharges #284's sweep-inverse contract question at the object
> level**: on the source subspace (bulk source, zero trace — every production rhs
> today), the sweep IS forward substitution on the assembled
> walk-order-triangular matrix. The triangularity leg ``triu(PᵀMP, 1) == 0`` is a
> structural EXACT zero (**the one 0-tolerance assertion in the family**).

The certificate covers **slab + 2-D Cartesian (DD and LD) + sphere + cylinder
(augmented)**. It does **not** cover d≥3 (no d=3 fixture), and it does **not**
cover a reflective-both-faces slab with `B` folded in (the SCC case).

---

## 7. One more subtlety: WHICH `(L+C)` is triangular?

`(L+C)` in ORPHEUS is endomorphic on the **composite** `V_bulk ⊕ V_trace`
(`streaming.py:278-310`), not on the bulk alone. Three distinct claims live
under one word:

- **(A) The per-(ordinate, group) BULK block** is cell-(block-)triangular under
  the walk order — the theorem of `path_integral.rst:1070-1077`, certified by G2.
  Also exactly block-diagonal over (ordinate, group) on Cartesian —
  `test_g1_assembled_matvec_equals_apply` asserts "leaked off-block" == 0
  (`test_assembly_mode.py:173-178`).
- **(B) The composite (bulk ⊕ trace)** is triangular only because `L_full`'s
  trace–trace `(s,s)` block carries the **identity** — the inflow identity row
  `I·ψ.inflow` plus the outflow defect row. `path_integral.rst:1051-1058`:
  "`L_full`'s unit-lower-triangular trace block carrying the identity that makes
  the whole within-group operator triangular (**load-bearing twice over, Issue
  #298**)". Issue #298 (CLOSED) was precisely a doc defect writing that block as
  `0`; its title says it: "*and that identity is what makes the walk triangular*".
- **(C) The augmented (bulk ⊕ trace ⊕ ψ½)** grid, curvilinear only — §4(a).

A name that asserts triangularity is asserting (A)+(B) jointly, on the *current*
admitted mesh/BC scope.

---

## Verdict

**`L + C` is not always triangular. It is triangular in ORPHEUS today because of
a theorem with hypotheses, plus an architectural eviction, plus a scope
restriction — and only the last of those is enforced by a guard.**

The precise conditions, and how each is held:

| # | Condition for `(L+C)` to admit an exact one-pass forward substitution | How it is held today |
|---|---|---|
| 1 | The spatial mesh is a **structured tensor-product lattice** (`Mesh1D`/`Mesh2D`/`from_axes` axes) | **True by current scope — no guard, no seam.** No unstructured mesh type exists to violate it. Enforced only by the absence of an alternative. |
| 2 | The per-cell closure is **upwind** so dependency edges increase the lattice potential `Σ_a sign(Ω_a)·i_a` | **Structural / by construction.** The order is *generated from* the potential (`level_of = local.sum(axis=0)`), never verified against the matrix at runtime. Verified per-case at *test* time only (G2). |
| 3 | The boundary law contributes **no closed loop** in the (face, ordinate) digraph | **Held by ARCHITECTURE, not by a check**: `B` is a sibling operator, structurally outside `(L+C)`. The `creates_sweep_cycle` flag exists but is **read by nothing in `orpheus/`**. |
| 4 | If a cycle-creating law is present, it is one ORPHEUS can evict or schedule | **ENFORCED — the one real guard.** `SNMesh.BOUNDARY_OPERATOR_REGISTRY` admits only `{vacuum, reflective}`; `resolve_boundary_conditions` raises `ValueError` for `periodic`/`white`/`albedo`/`mixed`/`prescribed`. This is a **mesh-construction** guard, not an operator guard. |
| 5 | On curvilinear, the ψ½ starting-direction flux is a **first-class state block**, so the triangularity is in the AUGMENTED space | **Enforced for the composite** by `CoupledOperator._triangular_orientation()` + `CoupledSubstitutionOperator.__init__`'s raise — the only runtime triangularity check in the codebase. Bare `(L+C)` on a carrying mesh is the ray-DECOUPLED `(A,A)` diagonal block; the joint object is the grid. |
| 6 | σ > 0 everywhere | **ENFORCED** at `StreamingCollisionOperator.__init__` (`streaming.py:593`) — but this is invertibility of the *diagonal*, not triangularity. |
| 7 | Multi-moment closures are triangular only at **CELL-block** granularity (dense `2^d × 2^d` diagonal blocks) | **True by construction**, certified by the LD G2 gate — which uses `lu_solve`, not `solve_triangular`, i.e. the diagonal blocks are genuinely dense. |

**Enforced by types/guards:** #4 (mesh BC registry — hard `ValueError`), #5
(coupled-grid orientation — hard raise), #6 (σ>0 — hard `ValueError`).

**Merely true-by-current-scope, with no guard anywhere:** #1 (structured mesh),
#2 (upwind closure ⇒ acyclic), #3 (`B` evicted rather than checked). These are
the three that a future unstructured mesh, a non-upwind closure, or a
periodic/lattice BC would break — and all three would break it **silently at the
operator level**, surfacing only as a wrong answer or as a red assembly gate.

The corpus already states the honest position, and it is worth quoting as the
one-sentence answer to the mission question
(`docs/theory/foundations/path_integral.rst:1078-1092`):

> Acyclicity is a property of the **(mesh, closure, boundary) triple**, and
> ORPHEUS certifies it *per case* rather than assuming it ... **A theorem about a
> triple can be false for a new triple; certify, don't assume.**

And the certificate has fired in anger: #282's back edge produced a cold residual
of `5.2 × 10^5`.

### Gaps found (not asked for, but they bear on the naming decision)

1. `assert_cycles_are_declared` is named in `sweep_graph.py:49` and implemented
   nowhere.
2. `creates_sweep_cycle` is a `ClassVar` that **no production code reads** — its
   only attribute-read consumers are two test files asserting its values. The
   uncommitted plan on this branch already calls it "a conservative hint"
   (`.claude/plans/sn_operator_realization_and_posing_plan.md:854`).
5. **Doc-internal tension on exactly this question.**
   `docs/theory/foundations/boundary_conditions.rst:1832-1840` states flatly that
   a reflective face makes "the sweep DAG acquire a cycle";
   `docs/theory/foundations/path_integral.rst:1093-1103` calls that exact
   proposition "*false as stated*" and demands a **closed loop** (SCC), citing
   the r=0 pole mirror as the counterexample. Both pages are current. Given the
   project's ratified internal-consistency prime directive, the BC page's
   sentence should be qualified to "a reflective face contributes a trace
   back-edge; a *cycle* requires a closed loop" — this is a live archivist item,
   and it is the same conflation the naming decision has to avoid.
3. The unstructured seam has **zero code presence** — no Protocol, no raising
   stub. It exists only in `docs/theory/methods/sn/placement.rst:127-130`, the
   grand report, and the DSA literature memos.
4. `SweepOperator`'s "ctor guard is the `SweepInvertible` TYPE itself" is a
   `TYPE_CHECKING`-only `Union` (`sweep_operator.py:69-80`) — erased at runtime,
   so the claimed guard is static-only.

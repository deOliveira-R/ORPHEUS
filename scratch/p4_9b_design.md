# P4.9b — design memo (2026-08-28, post-P4.9a `7a0f434c`)

The operator is posed with its two closures; the mesh returns to being a
phase space. Charter: `streaming_path_says_what_it_is.md` §P4.9 rows 4–5.
Inherited rulings binding this design (all user, 2026-08-28 unless noted):

- **The 4-arg shape**: `StreamingOperator` constructed from *(domain,
  codomain, spatial closure, angular closure)*; schemes are factories
  returning a `SpatialClosure`, the angular scheme returns an
  `AngularClosure` (charter §5c, R15's `<Axis-role>Closure` family).
- **The migration lever is a PRODUCTION posing-head factory** (Frame/stage-2
  pattern) — tests consume the same production site; a tests-only builder
  would twin the assembly.
- **The operator-minted scan table**: the fused table is the OPERATOR's
  artifact, minted near the hot path from the two closures it holds; the
  walk consumes operator artifacts plus one named kernel.
- **The hot loop stays call-free** (fp(4,8) fires the march recorder ZERO
  times — the is-identity gate's control leg).
- **Item (d), explicitly this round's question** (§5c): whether the operator
  keeps a provenance accessor to its generator — *"StreamingOperator should
  have all information it can leverage to be tested and diagnosed."*
- **One-mint doctrine** (§5c): the scheme induces on BOTH sides (space and
  operator), and *the two inductions must be minted together, at one site* —
  ERR-026's shape is what happens when they are not.

## 1. The measured ground (`[M]` 2026-08-28 at `10314dfa`, mine)

- `StreamingOperator` today: `@dataclass`, ONE field `sn_mesh`;
  `domain`/`codomain` derived (`sn_mesh.full_field_space`, :282/:311);
  `loss_representation = default_for(self.sn_mesh)` cached_property (:407);
  `is_adjointable` reads `type(self.sn_mesh.scheme).has_transpose_kernel`
  (:277). The ONE production ctor site is `build_streaming_collision`
  (`orpheus/sn/coupled_system.py:441`) — ⚠ the memo's
  "`orpheus/numerics/coupled_system.py:441`" was a path slip; the numerics
  module is the generic layer, the SN builder family lives in
  `orpheus/sn/coupled_system.py`.
- **The posing family already exists**: `build_streaming_collision(sn_mesh,
  mat_xs)` ("THE one LC spelling"), `build_within_group_system` →
  `WithinGroupSystem` (the posed record every within-group solve consumes),
  `_adjoint_posing_parts` (the daggered twin). Solve entries take `scheme=`
  and thread it INTO THE MESH via `_as_sn_mesh` (`solver.py:170/:172`,
  consumed at `:3150`/`:3716`).
- **Method choice enters at the mesh ctor today**: `_init_core` takes
  `scheme=` (instance, default `DiamondDifference()`, stored :268) and
  `pole_angular_closure=` (a CLASS, coord-dispatched default, bound
  `closure_cls(angular, pairing)` at :422 from `self.reduced` /
  neutral elements — the closure itself is mesh-free since Phase B).
- `[M]` kwarg triage by AST (validated against `solver.py`'s
  `SNMesh(scheme=)` as positive control): tests pass `scheme=` to **SNMesh
  directly at 36 sites in 15 files** (13 in `test_mms_ld_2d.py`), to local
  mesh helpers (`_slab` 20, `_mesh` 2, `_curvilinear_mesh` 2 — the edit is
  in each helper body), and to `solve_sn_fixed_source` at **19 sites**
  (STABLE — the entry keeps the kwarg, re-threading it to the posing head).
  `pole_angular_closure=` into SNMesh: 4 test sites (2 `_test_helpers`,
  2 `test_phase_c_gates`).
- `default_for(mesh)`: production calls it ONCE (the operator, :407);
  tests call it ~28× in 6 files (selection pins + representation access).
  Its `supports(mesh)` admission reads the scheme off the mesh (the
  curvilinear-LD rejection).

## 2. ⭐ THE FINDING THAT RE-SIZES ROW 5 — the space-side induction is real
## and already in the tree

Row 5 ("SNMesh sheds `scheme` + `pole_angular_closure`") is NOT a
two-field deletion. `[M]` the scheme co-determines the mesh's SPACES —
§5c's space-side induction, realized:

| site | reads | induces |
|---|---|---|
| `augmented_mesh.py:1230-1235` | `scheme.is_multi_moment`, `scheme.moment_axis(ndim, coord)` | `angular_trial_space` (LD's φ̂ moment tail) |
| `:1285` | (via `angular_trial_space`) | `full_field_space` — the operator's domain/codomain |
| `:1363` | `face_moment_count(scheme.spatial_basis_per_axis, ndim)` | the TRACE layout's face moments |
| `:587` | `type(self.scheme) is type(other.scheme)` | `SNMesh.__eq__` — mesh identity includes the method |
| `transport/fields/_bases.py:257` | `getattr(mesh, "scheme", None)` → `scheme.moment_axis(...)` | a field's own moment-tail width re-mint (L2; raises LOUD if absent — message must move) |
| `solver.py:992` | `sn_mesh.scheme.spatial_basis_per_axis` | (solver-side moment plumbing) |

⟹ shedding the generator obliges a decision about where the INDUCED SPACE
FACTS live. Three options in §4/F3 below; the recommended one keeps the
mesh as the carrier of the *induced* space-side products (stage-2
forgetting: retain the induced parts, forget the generator) and makes the
posing head the ONE site that receives the generator and distributes both
inductions.

## 3. Row 4/5 recount (`[M]` explorer 2026-08-28 at `10314dfa`,
## `scratch/p4_9b_row45_recount.md` — the detail lives there)

- **Row 4: 136 ctor sites** (prior ~150; the −14 is P4.9a's test
  shrinkage): **1 production + 135 across 40 test files** (top:
  `test_streaming_collision_operator.py` 25, `test_streaming_operator.py`
  19, `test_ld_adjoint_deferral.py` 10). AST census with positive AND
  negative controls; **zero aliased imports, zero getattr-string
  dispatch** — the call-spelling surface is closed. Extra §6b members: the
  `isinstance` refusal in `StreamingCollisionOperator.__init__`
  (`streaming.py:576`, unaffected by the ctor change), the
  `OperatorSum[...]` string subscript (`:450`), 2 production imports,
  ~127 prose mentions (the step-4 archivist sweep's population).
- **Row 5: 67 production reads** (prior ~60) = **5 (i) cache-build +
  43 (ii) per-cell dispatch + 17 (iii) posing-time + 2 (iv)
  reach-through**, plus 6 `sweep_graph` OWN-field reads (downstream —
  filled at 7 `scheme=self.mesh.scheme` kwarg sites inside the walks,
  which are themselves in class (ii)) and 7 SNMesh-internal sites (the
  §2 space-side table). Class (ii) is dominated by `_OneDimScanWalk`
  (39 in `loss_representation/__init__.py`, incl. the matvec kernel
  dispatch `sn_mesh.scheme.residual_kernel_batch(` at `:3259`);
  `streaming.py`'s 2 dissolve into the operator's own fields at row 4;
  `radial_characteristic.py:877/:951` read the closure at APPLY time
  (that operator gets HANDED the closure at posing). Class (iv), both
  L2: `transport/radial_characteristic_field.py:360`
  (`.level_indices` off the mesh closure) and `transport/fields/_bases.py:257`
  (the string-form `getattr` — found only by the string-form scan; ⚠ 74
  of 160 raw attribute hits were module-path references, the filter
  hazard measured).
- **Closure provenance**: three mesh-ctor surfaces (`__init__:206`,
  `from_axes:687`, `from_material_mesh:764`) funnel to `_init_core`;
  scheme default `DiamondDifference()` at `:268` is the ONLY name-keyed
  concrete ctor call in production; the closure default map is
  `closure.py:2260` (`default_angular_closure_class` — ⚠ its docstring
  still describes the retired `cls(sn_mesh)` contract: present-tense-false,
  fix on sight in step 3). `closure_cls(angular, pairing)` is mesh-free
  (Phase B) — the head minting from an existing mesh reads
  `sn_mesh.reduced.{angular, redistribution_pairing}` / the Cartesian
  neutral elements, exactly the `_init_core` body. Registry
  `.create(key)`: 0 production callers. Residual `isinstance(closure,
  MorelMontryAngularSweep)` dispatch in the walk (`:4213`/`:4630`) —
  pre-existing, out of scope, survives the re-plumb reading the handed
  closure.
- **Cache path**: `StreamingCoefficientCache.from_mesh_and_quad(sn_mesh)`
  (`cache.py:258`) is built at `solver.py:1434` (`SNSolver.__init__`) or
  LAZILY at `loss_representation/__init__.py:3791`, **memoized onto the
  mesh as `_geom_cache`** — two build sites, one mesh-attr memo slot; the
  operator-minted design single-sources this (the operator's cached
  artifact replaces the mesh memo). It reads ONLY `pole_angular_closure`
  (`:384`, the P4.9a mint reads — `2340db08` git-verified).
  `CollisionCache.from_geometry(geom, sig_t, scheme)` already takes the
  scheme EXPLICITLY at its 3 sites (`solver.py:1440/:1486`,
  `loss_rep:3860`) — the σ-stratum is already handed, corroborating the
  handing design. The walk consumes `geom`/`coll` as explicit `_run`
  parameters. ⚠ Ride-along fix at step 2: `cache.py:273` is a bare
  `assert` (stripped under `-O`) — convert to a `raise` while editing
  the file (coding-standards).

## 4. The design (proposed — every fork marked, none ruled yet)

### D1 — the ctor (fork F1)

`StreamingOperator` fields become:

```python
@dataclass
class StreamingOperator(LinearOperator["FullField"]):
    sn_mesh: "SNMesh"                                # the geometric substrate
    spatial_closure: "DiscretizationSchemeBase"      # required, no default
    angular_closure: "PoleAngularClosureBase"        # required, no default
```

- Done-when bullet 3 holds: `StreamingOperator(sn_mesh)` is a `TypeError`
  at every un-migrated site — loud at collection, never silent.
- `domain`/`codomain` stay DERIVED. ⚠ This deviates from the ruled literal
  *(domain, codomain, spatial, angular)* — the deviation and its reason
  are F1, for ruling. Reason: (a) the representation/walk is mesh-welded
  until O-3 (`default_for`, volumes, bc, connection — the operator cannot
  function from spaces alone today); (b) spaces do not carry the mesh, so
  passing `(domain, codomain)` ALONGSIDE the mesh makes the mismatch
  spellable — the Pattern-4 inversion; (c) `ScatteringOperator` (the §5c
  precedent) passes data + `space:` — and its space field exists because
  it has NO substrate object to derive from; L has one.
  When O-3/CS5 un-welds the representation, the mesh field recedes and
  the literal 4-tuple becomes reachable.
- `is_adjointable` reads `type(self.spatial_closure).has_transpose_kernel`;
  `loss_representation` becomes `default_for(self.sn_mesh,
  self.spatial_closure, self.angular_closure)`.

### D2 — the posing head (fork F2)

A free function in the existing posing family, `orpheus/sn/coupled_system.py`:

```python
def build_streaming_operator(
    sn_mesh, *, scheme=None, pole_angular_closure=None,
) -> StreamingOperator:
```

- Mints defaults exactly as `SNMesh._init_core` does today (DD;
  coord-dispatched closure class; binds `closure_cls(angular, pairing)`
  from `sn_mesh.reduced` / neutral elements) — the mint MOVES, verbatim.
- GUARDS the one-mint agreement (see D3): the scheme's space-side products
  must match the mesh's retained ones; mismatch raises with both named.
- `build_streaming_collision` gains the two kwargs and routes through it;
  `build_within_group_system` and the solve entries thread `scheme=` here
  instead of into `_as_sn_mesh`.
- Tests migrate `StreamingOperator(sn_mesh)` →
  `build_streaming_operator(sn_mesh)` — same production site as production
  (the ruled lever). Naming follows the family (`build_*`), so the set
  stays greppable.

### D3 — the space-side resolution (fork F3, the load-bearing one)

**Recommended: the mesh retains the INDUCED space facts and forgets the
generator; the head is the one generator site.**

- `SNMesh._init_core` keeps accepting `scheme=` but stores only the
  space-side products — `self.spatial_basis` (the scheme's basis
  descriptor: `is_multi_moment`, `moment_axis(ndim, coord)`,
  `spatial_basis_per_axis`) — and does NOT retain the instance. The five
  space-minting reads re-point at the retained product. `__eq__` compares
  the retained product (semantic gain: two meshes with identical phase
  spaces and different methods are EQUAL — the method is the operator's).
- `pole_angular_closure=` on the mesh RETIRES outright (its 4 test sites
  migrate to the posing head; the mesh binds no closure — :422 dies).
- The posing head receives the generator and distributes: space-side to
  the mesh it builds (production entries) or VALIDATES against a
  pre-built mesh (the guard, for tests that assemble separately); the
  operator-side closure to `StreamingOperator`.
- Alternatives (not recommended): (b) mesh ctor takes the product
  explicitly (`spatial_basis=...`) — purer, but migrates all 36+helper
  sites for no behavioural gain and makes the common DD case ceremonial;
  (c) full space-first re-plumb (fields minted off the posed space, mesh
  scheme-free) — CS5-sized, not P4.9b.

### D4 — the representation + cache re-plumb (fork F4: scope)

- `LossRepresentation` gains the two closure fields; `default_for(mesh,
  spatial_closure, angular_closure)` selects with the HANDED scheme; the
  walk's internal `self.mesh.scheme` reads (13) → `self.spatial_closure`;
  `sn_mesh.pole_angular_closure` (:3147) → `self.angular_closure`;
  `sweep_graph`'s own `scheme` field is filled from the representation.
- The scan cache: `[M]` every field is σ-free (geometry × quadrature ×
  closure algebra) ⟹ mintable by L alone. The ruled direction: the
  OPERATOR owns the artifact — a cached property on `StreamingOperator`
  minted from `(sn_mesh, spatial_closure, angular_closure)`; the walk
  consumes `op`'s artifact. [Exact current build path: explorer §4 —
  integrate before sizing the step.]

### D5 — item (d) discharges structurally

With D1 the operator holds the mesh, both closures, and (D4) the minted
table — every input it computes from is a readable field. The "provenance
accessor to its generator" needs no NEW surface: TODAY the spatial closure
IS the scheme instance (the closure/factory family split is O-3's retained
task), so `op.spatial_closure` is the generator access, honestly named
after the role it plays HERE. Proposal: record item (d) as discharged by
construction, no accessor minted; O-3 revisits when the family splits.

### D6 — hazard sites (each gets its own §6b row at execution)

- `transport/fields/_bases.py:257` — the `getattr(mesh, "scheme", None)`
  arm re-keys on the retained product; its raise-message re-words (grep
  the shortest distinctive fragment for pins).
- `transport/radial_characteristic_field.py:360` reach-through — [explorer
  classifies; expected: re-plumb to a handed closure].
- `SNMesh.__eq__` semantics change (type(scheme) → retained product) —
  enumerate __eq__/__hash__ consumers (caching keyed on mesh identity!).
- The `-k`/marker-registry surfaces naming `scheme` on the mesh
  (docstrings, theory pages) — the archivist sweep at landing.

## 5. Step decomposition (draft; §6b-staged so ONE closure instance ever exists)

⚠ **The interval hazard that fixes the staging**: the is-identity gate
(`test_angular_closure_is_single_object`) asserts every per-cell consumer
reaches `sn.pole_angular_closure` — the MESH's instance. Any interval in
which the posing head MINTS a closure while the mesh still BINDS one holds
two instances, and the gate (or bit-identity, or both) breaks mid-phase.
⟹ the head's defaults are transitional in step 1 and flip in step 3:

- **Step 1 — the ctor + the posing head; the head READS, it does not yet
  mint.** `StreamingOperator` gains the two required closure fields;
  `build_streaming_operator(sn_mesh, *, scheme=None, pole_angular_closure=None)`
  defaults by READING the mesh's existing bound objects
  (`sn_mesh.scheme`, `sn_mesh.pole_angular_closure`) — no new mint, one
  instance ever, the is-identity gate passes UNREWIRED (the operator's
  closure IS the mesh's object). `build_streaming_collision` + solve-entry
  threading + ALL ctor-site migrations land here (the §6b unit of the
  signature change; every un-migrated site is a loud `TypeError` at
  collection). Transitional read-off-mesh default is declared scaffolding,
  dies in step 3 (one merge cycle — aggressive-retirement compliant).
  Gate: the ctor-unconstructability witness (2-arg call raises) + full
  fast set bit-identical (same instances flow).
- **Step 2 — the representation, walk, and cache consume the OPERATOR's
  fields.** `default_for(mesh, spatial, angular)`; the representation
  gains the two fields; its 13 internal `self.mesh.scheme` reads +
  `:3147` closure read re-plumb; `sweep_graph`'s own `scheme` field fills
  from the representation; the scan cache mints as the OPERATOR's
  artifact (σ-free — Stratum-1's own docstring already says geometry-only;
  the solver's held cache re-plumbs to the operator's). Same instances,
  bit-identity pins + the minted-field `array_equal` gates + frozen corpus.
  The is-identity gate re-points its identity target at
  `op.angular_closure` (== `sn.pole_angular_closure` during the interval —
  both asserts hold, no flag day).
- **Step 3 — the mesh sheds; the head's defaults flip to MINTING.** The
  `_init_core` mint body (scheme default at :268, closure binding at :422)
  MOVES verbatim into the head; `pole_angular_closure=` retires from the
  mesh ctor (4 test sites migrate); `scheme=` becomes extract-and-forget
  (D3: retained induced product `spatial_basis_per_axis` + minted
  `moment_axis`); the five space-side reads re-point at the product;
  `__eq__`'s last conjunct compares the product; `_bases.py:257` re-keys
  (message re-word — grep the pin fragment first); `solver.py:992` reads
  the mesh product. §6c witnesses: `hasattr(sn_mesh, "pole_angular_closure")`
  is False; a DD-mesh and an LD-mesh compare UNEQUAL via the product while
  two DD-built meshes stay equal; the `_bases` raise still fires on a
  scheme-less carrier with the new message.
- **Step 4 — docs + audit.** Archivist sweep (theory pages naming the
  mesh as method-carrier), `dead_references` before done, the battery
  (owner-mutation superset per the M1 shape), the landing record.

Rebind/perf note (D4): today the cache lives on the SOLVER and "survives
rebinds; Stratum 1 is geometry-only" (`solver.py:1447`). ⛔ **"numpy
assembly, measured-cheap" was NOT measured and is REFUTED (verification
plan F2, `[M]` timed): the operator is built 6×/solve on the slab
(10× sphere), so a per-operator memo costs up to 24.65 % of a slab solve
(GL16/nx=200); today's build count is exactly 1.** The memo LIFETIME is
Q1 (blocking) — the perf gate is a COUNT pinned to the ruled number, not
a wall clock. The rebind contract itself is untouched (L is σ-free; C is
rebuilt fresh per solve already), but its ONLY witness asserts on the
retiring solver slots and must be re-posed in the same commit (F4).

## 5b. Adjacent open issues (checked 2026-08-28 — not in scope, must not conflict)

- **#402** (LD moment-tailed residual mint, un-skip the exit certificate):
  reads the same space-side induction D3 re-poses — the retained-product
  design must leave `moment_axis`/`spatial_basis_per_axis` as reachable as
  the scheme made them (it does: the mesh carries the induced product
  directly, which is STRICTLY easier to consume than `mesh.scheme.…`).
- **#398** (SNMesh ctor bypasses the quadrature admissibility relation):
  a ctor-surface repair adjacent to D3's ctor edit — do not fix here, but
  do not make it harder: the D3 edit keeps the ctor funnel at `_init_core`.
- **#392** (`face_areas` shim retirement, 5 test files): same file
  (`augmented_mesh.py`), independent — a candidate ride-along at step 3
  ONLY if the user wants it; otherwise untouched.

## 6. Forks for the user

- **F1** ctor spelling: substrate + 2 closures (recommended) vs literal
  4-arg (unreachable until O-3) — rule the deviation.
- **F2** posing surface: `build_streaming_operator` free function in the
  family (recommended) vs classmethod `StreamingOperator.pose`.
- **F3** space-side: mesh retains induced product, head owns the
  generator, guard for separate assembly (recommended) vs explicit
  product kwarg vs defer-to-CS5.
- **F4** operator-minted cache: in P4.9b (ruled direction; recommended)
  vs a P4.9c follow-on if sizing demands.
- **F5** item (d): discharge-by-construction (recommended) vs mint an
  explicit accessor now.


## 7. ⛔⛔ DESIGN REVISED at the fork round (user, 2026-08-28) — §§4–6 above
## are superseded where they conflict; they STAY per plan-authoring §3

The user's steer, verbatim anchors in brackets:

1. **SNMesh is the SAVE-STATE / DATA HUB, and it KEEPS the scheme.**
   ["SNMesh is more of a save state that we can dump… It would be keeping
   the discretization scheme in SNMesh not because it needs the
   discretization scheme, but because it stores shared machinery. For
   example, DSA discretization scheme needs to be consistent with the Sn
   solver, so SNMesh can provide the same spatial discretization scheme to
   both."] ⟹ D3 (extract-and-forget on the mesh) is DEAD in all variants;
   the space-side induction stays hub-side, reading the retained GENERATOR
   ["It also induces the space as nodal or modal during space construction
   (and in particular the spatial axis construction)"]. No `__eq__` change,
   no `_bases.py:257` re-key, no step-3 mesh shed. The charter's row 5
   ("SNMesh sheds scheme + pole_angular_closure") is REVISED BY RULING —
   the plan row gets its ⛔ banner at landing.
2. **The stage-2 extraction happens at the OPERATOR.** ["Discretization
   Scheme should be a stage 2 generator. It gets passed as an argument to
   the StreamingOperator (and the DSA), which extracts the spatial
   closure, and leaves just an accessor for provenance."] ⟹ the ctor
   REQUIRES the scheme (["it shouldn't have a default spatial
   discretization scheme, because this is an active choice"]); the
   extraction is identity until O-3 splits the closure/factory family
   (spelled as a `spatial_closure` property over the stored scheme, the
   ready seam); the scheme field IS item (d)'s provenance accessor —
   discharged by the user's own words.
3. **`.pose` is the intermediate classmethod; the end state is the
   cross-method literal ctor.** ["StreamingOperator should have eventually
   a constructor that expects the arguments it needs. Maybe .pose is an
   intermediate class method to construct it while we do the entire
   migration to pose StreamingOperator as a proper operator. But the
   streaming operator of every transport method will need a domain,
   codomain, some sort of way to discretize space, and some will need a
   way to discretize angle. SNMesh is probably unnecessary."] ⟹ F2's
   free-function option is dead; the migration lever is
   `StreamingOperator.pose(sn_mesh)`; the arc's exit (recorded, not
   P4.9b) removes `sn_mesh` from the ctor.
4. **The operator supplies the walk's method-flavored needs; traversal
   stays solver/cost-side.** ["the walk is not a StreamingOperator concern
   but a solver strategy concern… but maybe the right way is for the
   streaming operator to do the walk, not the SNMesh. The StreamingOperator
   is a FullOperator… figure out in the StreamingOperator how to do the
   walk?"] ⟹ step 2's mechanism stands with sharpened framing: the
   representation is SELECTED as strategy (cost-side seam survives) and
   FED by the operator (closures, minted table, kernels — answer-side);
   the mesh stays the geometry supplier, which the hub frame makes
   legitimate.
5. **Q4 ruled**: the operator-minted scan table lands in P4.9b step 2
   (option selected). ⚠ The accompanying note arrived TRUNCATED — "the
   ideal would be if." — completion requested before execution; do not
   guess its content.

### The revised step decomposition (shrunk — the mesh-shed step is gone)

- **Step 1 — the operator's contract.** Ctor requires
  `(sn_mesh, scheme, pole_angular_closure)` with IDENTITY GUARDS
  (`scheme is sn_mesh.scheme`, `closure is sn_mesh.pole_angular_closure`)
  — the args are the END-STATE contract, the mesh field + guards the
  DECLARED transitional weld (retired when O-3/CS5 remove the mesh);
  one closure instance ever, the is-identity gate untouched.
  `.pose(sn_mesh)` reads the hub's two objects and passes them.
  Migrate the 136 ctor sites → `.pose`; `build_streaming_collision`
  routes through it. ⭐ Solve entries DON'T re-thread — `scheme=` keeps
  flowing into the hub (`_as_sn_mesh`), which pose reads: the hub stays
  the active-choice site. `is_adjointable`/`spatial_basis_per_axis`
  reads dissolve into own fields. Bit-identical by construction.
- **Step 2 — the operator feeds the walk and owns the table.**
  `default_for` + representation ctor receive the closure pair from the
  operator; the 43 (ii) + 17 (iii) + 5 (i) reads re-plumb; the two
  reach-throughs are HANDED their objects; `op.scan_cache` becomes the
  minted artifact (Q4 ruling); the solver slot + `mesh._geom_cache` memo
  retire; `cache.py:273` assert→raise ride-along.
- **Step 3 — docs + battery + audit** (was step 4): archivist sweep
  (~127 prose mentions), `dead_references`, the owner-mutation superset
  battery, charter/plan REVISED banners for old row 5, landing record.

### Still-open (the next ask round)

- The truncated Q4 note ("the ideal would be if.") — complete it.
- The ANGULAR half of the hub: the scheme's hub residency is ruled; the
  bound `pole_angular_closure` staying hub-side (shared machinery: walk,
  seeding operator, cache all consume ONE instance) is the symmetric
  INFERENCE — confirm.
- "it shouldn't have a default" — read as the OPERATOR (ctor + pose pass
  explicit/hub objects). The HUB's own `DiamondDifference()` default at
  `_init_core:268` survives as the active-choice site's default — or
  should the hub default die too (blast radius: every default-DD SNMesh
  construction, hundreds of sites)? Confirm.
- The identity guards (strongest consistency: the operator's objects ARE
  the hub's) — confirm, or prefer looser space-side agreement guards.
- `SNMesh`-the-misnomer: file the rename as an issue (save-state/hub
  name), out of P4.9b's scope — confirm.


## 8. Second ask round — all four resolved (user, 2026-08-28)

1. **The truncated note, completed — a BINDING principle**: *"correctness
   concerns separated from performance concerns. So if a particular
   solution method goes through a scan in a certain way or welds terms for
   performance reasons, then this should be done as close to the solution
   strategy as possible. The algebra should be unweld and highly
   expressible for as long as possible, and performance optimizations
   should be lazily resolved."* ⟹ REFINES the table's home: the OPERATOR
   owns/exposes the ALGEBRA (the two closures + their minted per-ordinate
   constants); the fused/scan-normal TABLE is the STRATEGY's artifact,
   resolved LAZILY (first need) from the operator's objects. The eager
   build at `SNSolver.__init__:1434` and the mesh-attr memo
   `_geom_cache` (loss_rep:3791) both retire; the memo slot moves to the
   representation (cached per operator). "Minted near the hot path" from
   the earlier ruling reads correctly as strategy-side.
2. **Angular half**: hub keeps the bound instance (as today, :422); pose
   reads it; the operator requires it as a field. One instance ever.
3. **DD default**: the hub's `DiamondDifference()` default survives — the
   hub ctor is the active-choice site. The OPERATOR never defaults.
4. **Guards — the user's position tested and CONFIRMED; no operator-side
   guard.** The exercise (tried to beat "hub creates + passes ⟹
   unspellable"):
   - *Attack 1, the pose path*: unbeatable — pose reads the hub's own
     objects; no disagreement constructible.
   - *Attack 2, raw ctor with a fresh spatial scheme* (LD op on DD hub):
     spellable but LOUD — the operator's domain derives from the hub's
     spaces (no moment tail) while LD kernels index the tail → shape
     error at first apply. (Caveat noted: DD-op-on-LD-hub could in
     principle numpy-broadcast a trailing axis silently; requires
     deliberately crossing the seam.)
   - *Attack 3, raw ctor with a wrong-FAMILY closure* (`[M]`
     `IdentityAngularClosure.__init__` checks only
     `_require_single_moment_pairing` — a non-neutral curvilinear dome is
     ACCEPTED at construction). ⛔ **My "silent, plausible-wrong k"
     consequence claim was REFUTED by the test-architect (F12, `[M]` run,
     not reasoned): the doctored state RAISES at the first sweep — typed
     on the sphere (`TypeError: … requires the Morel-Montry closure`),
     untyped on the cylinder (`IndexError`), bit-inert on the slab (it IS
     the Cartesian default). The walk's residual `isinstance` dispatch
     (:4213/:4630) refuses the wrong family.** The no-guard verdict
     SURVIVES on the remaining grounds (the seam's legitimate doctored-
     probe use; pose-path unspellability); the refuted sentence must not
     reach the ctor docstring.
   - *Attack 4, cross-mesh closure smuggling* (mesh-A's closure into
     mesh-B's operator at equal N): silent wrong pairing; needs two
     meshes and deliberate crossing; guarding it needs the closure to
     remember its mesh — re-welding what Phase B un-welded.
   ⟹ VERDICT: the position holds. Production consistency is by
   construction (pose); the raw dataclass ctor is the DECLARED expert
   seam, its two residual arms (3, 4) recorded in the ctor docstring as
   what the seam does not check — articulation, not enforcement.

### FINAL P4.9b design (all ruled)

- `StreamingOperator(sn_mesh, scheme, pole_angular_closure)` — three
  required fields, no defaults, no guards; `spatial_closure` property =
  the identity extraction seam (O-3 splits the family later); the
  `scheme` field IS item (d)'s provenance accessor.
- `StreamingOperator.pose(sn_mesh)` classmethod — the intermediate
  posing surface reading the hub's two objects; `build_streaming_collision`
  routes through it; solve entries UNCHANGED (scheme= keeps entering the
  hub at `_as_sn_mesh`); 135 test ctor sites migrate to `.pose`.
- The hub (SNMesh) keeps `scheme` (generator; space induction;
  cross-consumer consistency incl. DSA) and the bound closure — NO mesh
  shedding; charter row 5 gets its ⛔ REVISED banner.
- Step 2: the representation receives the closure pair from the operator
  (`default_for(mesh, spatial, angular)` + fields); 43 (ii) + 17 (iii) +
  5 (i) reads re-plumb; the 2 reach-throughs are handed their objects;
  the strategy lazily resolves the scan table from operator algebra;
  `cache.py:273` assert→raise ride-along.
- Step 3: docs (~127 prose mentions), battery (owner-mutation superset),
  `dead_references`, plan/charter banners, landing record.
- End state (recorded, NOT this phase): cross-method
  `(domain, codomain, spatial-discretization[, angular-discretization])`
  ctor, `sn_mesh` out — rides O-3/CS5.


## 9. Third ask round — Q1 and Q5 ruled by criterion (user, 2026-08-28);
## picks delegated to the main agent within the stated bounds

**Q1 — the table memo's lifetime.** User's criterion, verbatim core: *"we
want the operator to be, well, the operator… if something changes during
solve, that is a solve problem… The important aspect is to separate the
algebra from the computation. This probably can only be realized when the
solver strategy itself is lazy… **Pick the thing that would be best at
surviving the change to a lazy solution strategy**"* — and the hint that
the `loss_representation` is itself RETIREMENT-BOUND (Campaign 2 rebuilds
the consumer side: GeneralizedEigenPencil, resolvent, partitioning).
⟹ **PICK: the interning mechanism lives IN the loss_representation layer**
— a module-level `WeakKeyDictionary` keyed on the hub, VALIDATED against
the (spatial, angular) closure-pair identity (rebuild on mismatch). Why it
best survives the transition:
- it sits in the strategy layer (= "as close to the solution strategy as
  possible", the ruled principle), so when the lazy solution strategy
  materializes at Campaign 2 the mechanism retires WITH the layer it
  serves — no machinery stranded on operator or hub;
- the operator stays pure algebra (nothing parked on it); the hub stays
  save-state (the `_geom_cache` ATTR-stamping retires — computation stops
  accumulating on the hub);
- count-1 per solve: all 6–10 operators of one solve share the hub's pair
  → one key → one build (the F2 COUNT gate pins **1**);
- the closure-pair validation dissolves F1's "a surviving cache masks the
  swap" trap: a doctored pair handed at pose gets its OWN build.
The eager `SNSolver.__init__` build + its `geom_cache`/`coll_cache`-slot
surface retire; F4's rebind-contract witness re-poses at the new tier IN
THE SAME COMMIT. `_coll_cache` / `_pole_mirror_cache`: untouched this
phase (Campaign 2 / O-3 territory).

**Q5 — naming.** User: *"make it symmetrical and easy to reason about…
Either it's spatial closure or it's spatial scheme. And it's either
angular closure or angular scheme… A cylinder does not have a pole… the
use of 'pole' seems like a misnomer. The important thing is that
something is spatial, and something is angular… symmetrical and
expressive (and grep works well)."*
⟹ **PICK: X = closure — the operator's and representation's two slots are
`spatial_closure` + `angular_closure`**, one symmetric greppable pair.
Today `spatial_closure` receives the hub's scheme INSTANCE (the
extraction is identity until O-3 splits generator/product — one docstring
sentence at the field); the hub keeps `scheme` (the GENERATOR — provenance
lives on the hub, reachable through the operator's transitional
`sn_mesh`; item (d) discharged with zero aliases). **The pole misnomer
dies**: hub attr `pole_angular_closure` → `angular_closure` (`[M]` 91
hits: 40 orpheus / 30 tests / 21 docs) and the class base
`PoleAngularClosureBase` → `AngularClosureBase` (`[M]` 119 hits incl. 66
docs) — both at step 3, riding the archivist sweep; MEMBER names
(`MorelMontryAngularSweep`, `IdentityAngularClosure`) untouched (only the
family-defining spellings move — the naming-consistency rule's
"swap only the colliding member"). Pole-flavored PROSE (~147 hits) is the
archivist's triage population, meaning-read per hit (a sphere's polar cap
is legitimately a pole; the CLOSURE vocabulary is not).

### Adopted from the verification plan without further ruling
Q2 = representation `.pose`-style classmethod (58 sites migrate to it);
Q3 = the route gate drives `(L+C).solve` itself; Q4 = selection
predicates operator-side (they consume the handed pair); Q6 =
`_pole_mirror_cache` untouched; Q7 = both permanent gates (AST + runtime
read-set); F5's partition = step 2's done-when (21 space-side reads STAY
hub-side); F9 = zero predicted-red; F10's second DiamondDifference
(the type-keyed diagnostic probe) carved out of identity gates BY NAME.

## 10. ⛔ Two corrections from the archivist's independent re-measurement
## (2026-08-29, its docs pass — §4's verify move, working again)

1. **`SNMesh.__eq__` does NOT compare the scheme type — §2's table row is
   WRONG on attribution.** `[M]` (archivist + verification plan F8):
   `SNMesh.__eq__` is `object.__eq__` (identity); the
   `type(self.scheme) is type(other.scheme)` conjunct I quoted at ":587"
   belongs to **`is_same_phase_space`**. Zero design consequence (the
   ruled design left the mesh untouched and no step keyed on `__eq__`),
   but the claim was published in this memo's §2 and the D3 option text —
   both now read with this banner.
2. **The F2 operator count is FIXTURE-DEPENDENT and my relays kept one
   fixture's number.** Architect at `10314dfa`: 6–10 operators/solve;
   archivist on its own fixtures: **38–43**/solve (the count scales with
   outer iterations — "independent of nx/ng/inner solver" was the
   architect's measured scope, not a universal). Consequence restated
   honestly: a per-operator memo costs 16.8–24.65 % on the architect's
   fixtures and **up to +68 %** on the archivist's. The RULED design
   made the number moot for cost (`[M]` builds == 1 on all four configs,
   both measurers) — the COUNT gate, not the percentage, is the
   instrument.

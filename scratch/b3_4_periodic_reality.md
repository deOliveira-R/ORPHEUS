# B3.4 — Periodic BC: what the code actually does

Branch `refactor/operator-strategy-layers`, HEAD `498f94b0`. Line numbers current
at that HEAD — re-derive via Nexus/grep if the tree has moved. Method: read the
bodies; docstrings are claims to be checked (L33), never evidence.

**Headline: the claim is FALSE. There is no face-pair / wrap-around indexing
anywhere in the SN sweep, schedule, DAG, or trace space. The "sweep handles it"
sentence describes a mechanism that does not exist and never has.**

---

## 0. Where the claim actually lives (premise correction)

The brief attributes the sentence

> "the SN sweep handles the spatial wrap via its own face-pair indexing and the
> BC operator only needs to pass the angular trace through unchanged"

to BOTH `orpheus/geometry/boundary/periodic.py` and `orpheus/numerics/operator.py`.
**That is half-stale.** A repo-wide grep for `face-pair` finds it in exactly two
production places:

- `orpheus/numerics/operator.py:2622-2624` — `PeriodicWrapOperator` docstring,
  verbatim: "*this matches the `PeriodicBoundary` semantics, where the SN sweep
  handles the spatial wrap via its own face-pair indexing and the BC operator
  only needs to pass the angular trace through unchanged.*"
- `orpheus/sn/boundary/realizer.py:466-468` — a code comment restating it:
  "*The current `PeriodicWrapOperator` body is identity-with-copy (the SN sweep
  handles the spatial wrap via its face-pair indexing)*".

`orpheus/geometry/boundary/periodic.py` **no longer carries that sentence** — it
was rewritten at `a0fd17b4` (boundary B1). It now makes a DIFFERENT, weaker, and
also-unsatisfied claim (periodic.py:35-44):

> "Realising periodicity at the *sweep* level requires coupling the two faces'
> boundary-flux buffers -- which is a sweep-orchestration concern not modelled by
> `apply` alone. … The two-face plumbing is handled by whoever instantiates
> `PeriodicBoundary` and orchestrates the sweep."

So the tree now carries TWO mutually-inconsistent prose claims: `operator.py`
says the sweep DOES the wrap; `periodic.py` says an un-named "whoever … orchestrates
the sweep" must do it. Neither corresponds to code that exists (§1).

**Third falsity, same file:** `operator.py:2636-2639` xrefs
`orpheus.geometry.boundary.periodic.PeriodicBoundary.apply` as the "legacy
aliasing-safety contract". **`PeriodicBoundary` has no `apply`** — it is a pure
descriptor (`periodic.py:42-44`; `_base.py:112-122` "No `apply`"; the strip landed
at `da414ebc`). This is a dangling Python-domain xref: exactly the class the `-W`
Sphinx gate does NOT catch.

---

## 1. Q1 — Any face-pair / wrap-around indexing in the SN sweep or trace machinery?

**NO. No such mechanism exists.** Every candidate home was read; all negative.

### 1a. Searches run
- `grep -rn "partner\|wrap\|opposite"` over `orpheus/sn/`,
  `orpheus/numerics/spaces/`, `orpheus/geometry/`. Every `partner` hit in
  `orpheus/sn/solver.py` (1569, 1578, 1615, 2145) is the **specular/angular**
  partner (the ψ½ mirror ordinate) or the `#291` scale partner — none is a
  spatial face partner. Every `wrap` hit outside the descriptors is Python
  object-wrapping.
- `grep` for `partner_face` / `opposite_face` / a literal `{"xmin": "xmax"}` map:
  **zero production hits.** The only `partner_face` mention is `_factors.py:385`
  saying the field was deliberately NOT added.

### 1b. `orpheus/sn/loss_representation/sweep_schedule.py` — no face→face map
- `_outgoing_faces(label)` (241-254) derives a face set from the octant's own
  sign per axis. One-directional; no partner.
- `_reflective_faces(sn_mesh)` (257-281) selects faces whose law
  `permutes_ordinates`. `SpatialWrap.permutes_ordinates` is **`False`**
  (`_factors.py:405-408`, comment: "*Spatial, not angular: ordinate n at face A
  feeds ordinate n at face B*"). So a periodic face is **excluded** from the
  reflective set, never enters any group's `reflect_faces`, and every one of its
  inflow rows falls into `B_upper` (lagged) — stated at `sweep_schedule.py:187-191`.
  That is the schedule correctly *declining* to realize a law it cannot express;
  it is not the schedule implementing a wrap.

### 1c. `orpheus/sn/loss_representation/sweep_graph.py` — the DAG never wraps
The named invariant `assert_face_pairing_consistent` is **interior** cell pairing
within an octant: "*the outgoing face of cell (i, j) matches the incoming face of
(i + sx, j)*" (`tests/sn/sweep/core/test_sweep_graph.py:13-15`). The rolling
frontier (`_FrontierPlan` / `_MovingFrontier`, sweep_graph.py:172-361) advances
along the octant direction and terminates at the domain edge; there is no
modular/wrap index anywhere in it. This is the interior 1-cochain `C¹_int`, not
the domain trace.

### 1d. `orpheus/numerics/spaces/angular_trace_space.py` — strictly per-face
`inflow_indices_for_face` (434-441), `outflow_indices_for_face` (443-450),
`_face_restrictions` (454-482), `outflow_restriction` (484-497),
`inflow_restriction` (499-506) — every one takes a single `face: str` and reads
one row of `omega_dot_n`. Nothing maps a face to another face. `_checked_face`
(508-511) only validates existence.

### 1e. `orpheus/sn/mesh/augmented_mesh.py` — no wrap at BC install
`realize_boundary_law` (397-442) builds an `SNMethodSpace.for_face(...)` for ONE
face and hands it to `SNBoundaryRealizer().realize`. The method space carries
that face's `inflow_indices` only; there is no channel for a partner face.

### 1f. The ONLY place in the repo that models periodic correctly is a derivation
`orpheus/derivations/discrete/sn/sweep_acyclicity.py:361-366` does the real thing:
```python
elif kind == "periodic":
    other = "right" if face == "left" else "left"
    tgt = TraceSlot(other, s.ordinate, "inflow")
```
That is a graph-theoretic acyclicity derivation over an abstract slab trace
digraph — **not the production sweep**, and nothing in `orpheus/sn/` consumes it.
It is evidence the concept is understood, and evidence it was never built.

### 1g. Contrast — the diffusion realizer REFUSES periodic, loudly
`orpheus/diffusion/boundary_realizer.py:295-315` dispatches on
`isinstance(law.geometry_map, SpatialWrap)` and raises `BoundaryError`:
"*its geometry map is a SPATIAL wrap, which couples a face to its OPPOSITE face
(a trace-block permutation) — the one geometry P1 cannot integrate away into a
per-face albedo*". Pinned by `tests/diffusion/test_boundary_realizer.py:301-305`
and `:394-395`.

**SN, facing the identical structural obstruction (a block-diagonal-over-faces
`B`), silently realizes an identity instead.** That asymmetry is the finding: the
diffusion side already wrote down the exact reason SN cannot do this either.

---

## 2. Q2 — Can `PeriodicBoundary` be reached from a user-facing SN configuration?

**Not from a `BC` tag. Refused, loudly, and the refusal is pinned by two tests.**

### 2a. The tag path (the only user-facing route)
`SNMesh.__init__` → `_init_core` → `self.bc = resolve_boundary_conditions(self)`
(`augmented_mesh.py:363-365`)
→ `orpheus/transport/method.py:218-262` `resolve_boundary_conditions`, which per
face calls `_law_from_tag(method, tag, label)`
→ `orpheus/transport/method.py:289-298`:
```python
law_cls = method.BOUNDARY_OPERATOR_REGISTRY.get(tag.kind)
if law_cls is None:
    raise ValueError(f"{type(method).__name__} does not support boundary condition ...")
```
The SN vocabulary is defined at **`orpheus/sn/mesh/augmented_mesh.py:175-178`**:
```python
BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]] = {
    "vacuum": VacuumInflow,
    "reflective": ReflectiveBoundary,
}
```
The comment right below (`augmented_mesh.py:186-195`) is honest and correct:
"*The 4 other kinds `SNBoundaryRealizer` dispatches today (`white`, `periodic`,
`albedo`, `prescribed_inflow`) are NOT registered here — so they are declarable
only by constructing the law directly, never from a `BC(...)` tag; admitting them
requires SN-sweep-side wiring (sweep cycles for periodic, etc.) and is issue #189.*"

Confirmed by test:
- `tests/sn/operators/test_snmesh_realizer_wiring.py:456-460` —
  `assert set(SNMesh.BOUNDARY_OPERATOR_REGISTRY) == {"vacuum", "reflective"}`
- `tests/sn/operators/test_snmesh_realizer_wiring.py:463-472` and
  `tests/sn/operators/test_boundary_conditions.py:97-103` — `BC("periodic")`
  raises `ValueError` matching `"'reflective'.*'vacuum'"`.

### 2b. Other routes (constructible, but nothing installs it)
1. **`SNBoundaryRealizer().realize(PeriodicBoundary(), method_space)`** — WORKS.
   `realizer.py:464-473` returns `stamp_boundary_role(PeriodicWrapOperator() & IdentityOperator())`.
   This is the route every existing test takes (§5).
   `assert_realizable` (`realizer.py:279-282`) does not object.
2. **`SNMesh.realize_boundary_law(law, face)`** (`augmented_mesh.py:397-442`) —
   a PUBLIC method that takes a **law object**, bypassing the tag registry
   entirely. `sn_mesh.realize_boundary_law(PeriodicBoundary(), "xmin")` returns a
   valid `_BoundBoundaryOperator`. It does not install it.
3. **`sn_mesh.bc` is a plain mutable instance dict** (`augmented_mesh.py:363`), so
   `sn_mesh.bc["xmin"] = sn_mesh.realize_boundary_law(PeriodicBoundary(), "xmin")`
   would install a periodic law post-construction. **Nothing in `orpheus/` or
   `tests/` does this** (grepped `sn.bc[...] =` / `sn_mesh.bc[...] =`: no
   assignment sites outside construction).
4. `realize_recursively(tree, ms, SNBoundaryRealizer())` — the rank-N law-tree
   walker (`orpheus/geometry/boundary/_realizer.py`) would also realize a
   `PeriodicBoundary` leaf. Same conclusion: realizable, never installed.

**Verdict: `PeriodicBoundary` is realizable but structurally unreachable in a
solve.** No production or test path ever puts one in `sn_mesh.bc`.

### 2c. What WOULD happen if one were installed (the latent defect)
Since it is one dict assignment away, it is worth being precise. In
`_reflect_trace` (`boundary.py:381-388`) the face action is
`ι₋ ∘ law ∘ γ₊`:
- `gamma_out.apply(face_in)` gathers the face's **own** |Γ₊| outflow rows;
- `PeriodicWrapOperator & IdentityOperator` passes them through (a copy);
- `gamma_in.apply_transpose(image)` scatters those values into the **same face's**
  |Γ₋| inflow rows, positionally.

Result: inflow row *i* of Γ₋ receives outflow row *i* of Γ₊, **on the same face**.
That is neither periodicity nor a specular mirror — it is a positional
self-face pairing with no physical meaning. And **it does not raise**, because
`|Γ₊| == |Γ₋|` on every quadrature × face pair in the tree (measured and
documented at `tests/sn/operators/test_b3_domain_narrowing.py:304-308`). Silent
wrong physics, gated only by the registry admission in §2a.

---

## 3. Q3 — What `SNBoundaryOperator` structurally assumes about per-face independence

The module docstring states it (`orpheus/sn/operators/boundary.py:45`): "*`B` is
therefore block-diagonal over faces — it never mixes faces.*" It is not a comment;
it is load-bearing in **eight** places. Every one breaks if a law must read a
different face.

1. **The loop** — `_reflect_trace`, 375-428:
   `for face, law in face_laws.items():` → `face_in = boundary.face_view(face)`
   (376) → `out_boundary.face_view(face)[...] = …` (386 / 397 / 426). The law's
   only input is `law.apply(gamma_out.apply(face_in))` (384): **one face's
   restricted array and nothing else.** There is no channel for partner data —
   `LinearOperator.apply(x)` is 1-arg by contract (`_base.py:112-135`).
2. **The output buffer** — `out_boundary = AngularBoundarySourceSink.zeros_on(mesh)`
   (347), written with whole-slot **assignment** `[...] =` (386, 426). Two faces
   writing into each other's slots would need accumulate semantics plus an
   ordering rule; neither exists.
3. **The `faces` subset argument** (290, 364-372). Exactness is justified ONLY by
   block-diagonality — comment at 351-354: "*`B` is block-diagonal over faces, so
   the subset action is the EXACT restriction (no cross-face coupling is
   dropped).*" With a face-coupling law, a subset excluding the partner silently
   drops a real coupling.
4. **The `rows` subset argument** (291, 373-374, 396-399, 421-425). Same argument
   one level finer — 360-362: "*Row-granular restriction is exact for the same
   reason the face restriction is: the projected action writes each target row
   independently.*"
5. **The transpose path** (400-428). `Bᵀ = ι₊ ∘ lawᵀ ∘ γ₋` **per face**, with
   `gamma_in` (420/425) and `gamma_out` (426) both taken from `trace.…(face)` for
   the SAME `face`. The transpose of an `f→f'` block is an `f'→f` block; this code
   cannot express it — it would scatter the transposed image onto the wrong face.
6. **`is_adjointable`** (256-272) — the per-face **intersection**
   `all(law.is_adjointable for law in laws)`, valid only for a direct sum.
   Note it queries the **realized operator**, which for periodic is the identity
   body → `True`. The **law factor disagrees**: `SpatialWrap.is_adjointable` is
   `False` (`_factors.py:410-417`) with the comment "*`False` reports an
   IMPLEMENTATION gap: the realized operator is currently an angular identity with
   the spatial pushforward unbuilt (#183). B3.4 builds it and this flips, WITH its
   gate.*" See divergence D2 in §6.
7. **`split`** (554-588). `upper_rows` is `setdiff1d` against **that face's own**
   inflow indices, looped `for face in self._face_laws`; the exactness claim
   (565-566) is per-face by construction.
8. **The in-place reflects** — `reflect_inflow_inplace` (512-552, `for face in
   selected:` → whole-face `[inflow] =` at 550-552) and
   `SNMaskedBoundaryOperator.reflect_rows_inplace` (920-960, `[rows] +=` at
   958-960). The latter is **additive and in place**: with a face-coupling `B`,
   face A's update would be read by face B's update inside the same loop — a
   read-write hazard the block-diagonal framing exists to rule out.

**Summary of what breaks:** the law has no way to *receive* partner data (1);
the write discipline has no way to *combine* two faces' contributions (2, 8); the
G-S `faces=`/`rows=` restrictions stop being exact and silently drop the
cross-face edge (3, 4); `apply_transpose` becomes structurally wrong rather than
merely absent (5); and the composite predicate rule loses its justification (6).

---

## 4. Q4 — Who calls with a `rows` / subset argument, and what is the subset?

### 4a. `rows=` (per-face inflow ORDINATE row subsets) — 2 call sites, 1 producer
| Site | Call | Subset |
|---|---|---|
| `boundary.py:918` | `SNMaskedBoundaryOperator.apply` → `inner._apply_faces(psi, "apply", rows=self.rows)` | `self.rows` — the half's per-face inflow rows |
| `boundary.py:954-956` | `SNMaskedBoundaryOperator.reflect_rows_inplace` → `inner._reflect_trace(bf, "apply", faces=tuple(selected), rows=selected)` | `self.rows ∩ faces` |

Producer: **`SNBoundaryOperator.split(schedule)`** (`boundary.py:554-588`).
- `lower_rows = schedule.lower_inflow_rows(self.sn_mesh)`
  (`sweep_schedule.py:170-214`) — per face `f`, the intersection of `f`'s inflow
  ordinates with the ordinates of every octant group swept **strictly after**
  `f`'s reflect group. Those rows read the FRESH current iterate ⟹ `B_lower`.
- `upper_rows` = per-face `setdiff1d(inflow_indices_for_face(f), lower_rows[f])` —
  the cyclic back-edges plus **every row of a never-reflected face (vacuum, white,
  albedo, periodic)** (`boundary.py:562-563`, `sweep_schedule.py:187-191`).

Consumers of the split:
- `orpheus/sn/solver.py:772-773`: `parts = B.split(SweepSchedule.gauss_seidel(sn_mesh))`
  → `return LC - parts.lower, (S, parts.upper)` — i.e. `M = (L+C) − B_lower`,
  gains `= (S, B_upper)`. Gated to **multi-D Cartesian, seedless** meshes
  (solver.py:755, 765-771).
- `orpheus/sn/operators/streaming.py:709-737` — `__sub__` dispatch on
  `SNMaskedBoundaryOperator` producing the `ScheduledInvertibleOperator`.
- `orpheus/sn/operators/scheduled_invertible.py:259-261` — passes
  `reflect=self.lower.reflect_rows_inplace` into the sweep's forward
  substitution.
- `orpheus/sn/loss_representation/__init__.py:4904-4910` — the walk calls
  `reflect(boundary_flux, group.reflect_faces)` after each octant group.

### 4b. `faces=` (face subsets, no row restriction) — 2 call sites
| Site | Call | Subset |
|---|---|---|
| `solver.py:2002` | `SNBoundaryOperator(sn_mesh).reflect_inflow_inplace(boundary_flux, faces=faces)` | `None` (whole trace) for the reconstruction sweep + direct SI seed; a face subset for the G-S inter-group reflect (documented solver.py:1979-1990) |
| `loss_representation/__init__.py:4910` | `reflect(boundary_flux, group.reflect_faces)` | `OctantSweepGroup.reflect_faces` — the reflective faces whose LAST outflowing octant group is this one (`sweep_schedule.py:150-164`) |

**Every one of these subsets is justified by block-diagonality over faces.** A
face-coupling law invalidates all four.

---

## 5. Q5 — Is there any test that exercises `PeriodicBoundary` end-to-end in SN?

**No. Zero solves. Every test is realization-only or descriptor-only.** Full
inventory (`grep -rn "PeriodicBoundary\|PeriodicWrapOperator\|SpatialWrap"
tests/ examples/`):

### (a) Genuine SN solves with periodic — **NONE**
No test constructs an `SNMesh` carrying a periodic law and calls any solver.
Structurally impossible via the tag path (§2a) and nothing takes the
`realize_boundary_law` / `bc[...] =` route.

### (b) Realization-only (SN)
- `tests/sn/operators/test_sn_boundary_realizer.py:474-495` — periodic realizes to
  `TensorProductOperator(PeriodicWrapOperator, IdentityOperator)`; asserts the type.
- `tests/sn/operators/test_operator_block_role.py:152, 158-170` — the realized op
  advertises `BlockRole.BOUNDARY`.
- `tests/sn/operators/test_b3_domain_narrowing.py:283, 288-334` — the Γ₊→Γ₋
  domain gate. **Periodic is in the `deferred=True` set ⟹ `xfail(strict=True)`.**
  The xfail reason (`:255-268`) is measured and blunt: "*albedo and periodic
  silently accept it and echo Γ₊ back — i.e. Γ₊ → Γ₊, the wrong codomain,
  invisible to a shape check because |Γ₊| == |Γ₋| on every fixture. Deciding
  their Γ₊ → Γ₋ action is B3.4's job (#183, #189).*"

### (c) Descriptor / factor tests (geometry)
- `tests/geometry/test_boundary_factors.py:83` — `(PeriodicBoundary(), SpatialWrap,
  ScalarResponse, 1.0, "periodic")`: the two factors are the declared types.
- `tests/geometry/test_boundary.py:469-480` — "*PeriodicBoundary is identity on
  the angular axis (smoke test)*"; `:590`, `:790` — registry key.
- `tests/geometry/test_bc_universal_invariants.py:656` — universal invariants.
- `tests/geometry/test_bc_equivalence_snapshot.py:402-410` +
  `_generate_bc_equivalence_snapshots.py:218-221` — a bit-exactness snapshot of
  the realized op, described in the generator as "*Identity …*". **This snapshot
  pins the IDENTITY body as a baseline** — i.e. it will go RED when B3.4 builds
  the real wrap. Flag it for the implementer.

### (d) Pure operator tests (numerics)
- `tests/numerics/test_periodic_wrap_operator.py` (whole file) — copy semantics,
  self-adjointness, shape passthrough. Line 68 repeats the stale claim: "*The
  legacy `PeriodicBoundary.apply` body is `psi_out.copy()`*" — that method no
  longer exists (§0).
- `tests/numerics/test_operator_capability_predicates.py:67, 192` —
  `("PeriodicWrap", PeriodicWrapOperator(), False, True, STRUCTURAL_ABSENT)`:
  not invertible, adjointable.

### (e) Refusal tests (diffusion — the useful precedent)
- `tests/diffusion/test_boundary_realizer.py:301-305` and `:394-395` —
  `pytest.raises(BoundaryError, match="OPPOSITE face")` for a bare law and for a
  `LawSum` tree.
- `tests/geometry/test_bc_errors.py:127, 148` — same refusal contract.

---

## 6. Docstring-vs-code divergences (each is a claim the code does not support)

**D1 — the claim under investigation. FALSE.**
Claim (`operator.py:2622-2624`): "*the SN sweep handles the spatial wrap via its
own face-pair indexing and the BC operator only needs to pass the angular trace
through unchanged.*"
Code: no face-pair indexing exists anywhere in `orpheus/sn/` (§1). The second
clause is true only in the vacuous sense that the operator IS an identity — but
its stated justification is false, so the identity is unbacked, not "sufficient".
Restated at `realizer.py:466-468`.

**D2 — periodic "advertises `apply_transpose`". Contradicted by its own factor.**
Claim (`boundary.py:220-221`): "*`apply_transpose` is advertised iff EVERY
per-face law advertises it — reflective (involution), vacuum, periodic and albedo
do*".
Code: `SpatialWrap.is_adjointable` is **`False`** (`_factors.py:410-417`) and says
so explicitly ("*`False` reports an IMPLEMENTATION gap*"). The composite predicate
reads the **realized operator** (identity → `True`), not the law, so the two
disagree. Same divergence at `boundary.py:259` and `boundary.py:679`.

**D3 — `PeriodicBoundary.apply` xref. Dangling.**
Claim (`operator.py:2636-2639`): xrefs
`orpheus.geometry.boundary.periodic.PeriodicBoundary.apply`.
Code: no such attribute — descriptors are not callable (`_base.py:112-122`).
Also repeated in prose at `tests/numerics/test_periodic_wrap_operator.py:12, 68`.

**D4 — "the two-face plumbing is handled by whoever … orchestrates the sweep."**
Claim (`periodic.py:42-44`).
Code: nobody does. There is no such orchestrator, in production or tests. This
sentence is a placeholder that reads as a delegation.

**D5 — `_reflect_corner`'s law list is right for the wrong stated reason.**
`boundary.py:731-736` puts periodic in the "loud-deferred" set. That IS the
current behaviour and it is correct — but the stated basis
(`_has_ruled_corner_action`, boundary.py:118-121: "*a spatial wrap needs the
partner face's ray*") is the same partner-face concept the module elsewhere claims
the sweep provides. The two statements cannot both be true. `_reflect_corner` is
the honest one.

---

## 7. Answers, condensed

1. **Face-pair / wrap indexing in the SN sweep?** **No such mechanism exists.**
   Not in `sweep_schedule.py`, not in `sweep_graph.py` (its "face pairing" is
   interior cell-to-cell), not in `angular_trace_space.py` (strictly per-face),
   not in `augmented_mesh.py`. The only correct model of periodic in the repo is
   a derivation script (`derivations/discrete/sn/sweep_acyclicity.py:361-366`)
   that production does not consume.
2. **Reachable from a user-facing SN config?** **No.** `SNMesh.BOUNDARY_OPERATOR_REGISTRY`
   = `{"vacuum", "reflective"}` (`augmented_mesh.py:175-178`), enforced in
   `_law_from_tag` (`transport/method.py:289-298`), pinned by two tests.
   Constructible by direct realization (`realizer.py:464-473`) and installable in
   principle via the public `realize_boundary_law` + the mutable `bc` dict, but
   **no code does it**. If it were installed it would silently produce a
   meaningless same-face positional pairing (§2c).
3. **Per-face independence assumptions in `SNBoundaryOperator`:** eight
   structural sites (§3) — the loop, the assignment-write buffer, `faces=`,
   `rows=`, the transpose, `is_adjointable`, `split`, and both in-place reflects.
   A law needing another face's data has **no channel to receive it** and would
   break the exactness justification of every subset restriction.
4. **`rows`/subset callers:** `SNMaskedBoundaryOperator.apply` (918) and
   `.reflect_rows_inplace` (954), both produced by `SNBoundaryOperator.split`
   from `SweepSchedule.lower_inflow_rows`; consumed by the multi-D Cartesian
   Gauss-Seidel path (`solver.py:772`, `scheduled_invertible.py:260`,
   `loss_representation/__init__.py:4910`). Face-only subsets: `solver.py:2002`
   and the walk's `group.reflect_faces`.
5. **End-to-end SN periodic test:** **none.** 11 test sites: realization-only (3
   SN), descriptor/factor (5 geometry, one of which — the bit-exactness snapshot
   at `test_bc_equivalence_snapshot.py:402` — pins the IDENTITY body and will go
   RED when the real wrap lands), pure operator (2 numerics), refusal (diffusion).
   `test_b3_domain_narrowing.py` already carries periodic as
   `xfail(strict=True)` — **B3.4's todo list is already spelled as a strict
   xfail**, so building the wrap flips that row green automatically.

# The boundary-condition realization chain — a complete level-by-level map

**Scope**: from a user declaring `BC("reflective")` down to the array op on a face,
for BOTH SN and diffusion.
**Tree state**: branch `refactor/operator-strategy-layers`, HEAD `73627b71`
(2026-07-30). All line numbers read from the working tree at that HEAD.
**Nexus caveat**: the graph DB was built at `8d996f53`; three commits since
(`ddc1ee10`, `1fd15f64`, `73627b71`) touched `orpheus/geometry/boundary/`.
Every structural claim below was verified by *reading the file*; Nexus was used
only as a lead generator. Divergences are flagged inline.

---

## 1. The level-by-level chain

### Level 0 — the tag type: `BC`

`orpheus/geometry/mesh.py:44-69`

```python
@dataclass(frozen=True)
class BC:
    kind: str
    params: dict[str, float] = field(default_factory=dict)
```

A **solver-agnostic string tag** plus a numeric parameter bag. The geometry
module makes no claim about what a `kind` means — semantics are resolved per
method (docstring `mesh.py:46-52`; note that docstring still says
`BC_REGISTRY`, the pre-#290-P7b name — the live name is
`BOUNDARY_OPERATOR_REGISTRY`. **Doc drift, cosmetic.**)

`BC` also carries an unrelated escape hatch `BC.to_alpha()`
(`orpheus/geometry/mesh.py:71-120`) that maps `vacuum→0.0`, `reflective→1.0`,
`partial→params["albedo"]` for the *trajectory-resolvent / Birkhoff–Sinai*
reference solvers. This is a **second, parallel tag→physics map** that does NOT
go through the law/realizer chain at all. Note `"partial"` is a tag that exists
here and nowhere in either method's registry.

### Level 1 — declaration: the axis endpoints

`orpheus/transport/mesh/axis.py`

- `Axis1D` Protocol — `axis.py:154-216`. Its BC surface is
  `bc -> dict[str, BC | None]` (`axis.py:203-204`) and
  `with_uniform_bc(bc) -> Axis1D` (`axis.py:206-216`).
- `AxisMesh` (Cartesian) — `axis.py:223-283`. Fields `bc_low`, `bc_high`
  (`axis.py:247-248`); endpoints `("min", "max")` (`axis.py:274-275`);
  `bc` property returns `{label_low: bc_low, label_high: bc_high}`
  (`axis.py:277-279`).
- `RadialAxisMesh` (solid sphere / cylinder) — `axis.py:286-357`. ONE endpoint
  `("outer",)` (`axis.py:343-345`), field `bc_outer`. **The pole at r=0 is
  deliberately NOT an endpoint** (rationale `axis.py:160-182`): it carries an
  angular-closure regularity condition, not a trace law, so a "BC that is not
  a BC" is unrepresentable.
- `FaceLabel` — `axis.py:88-147`. The `(axis_index, endpoint)` structural key.
  `face_name` (`axis.py:118-147`) renders `"{AXIS_NAMES[axis_index]}{suffix}"`
  through `_ENDPOINT_SUFFIX = {"min":"min","max":"max","outer":"max"}`
  (`axis.py:85`) — so a solid radial axis's `outer` renders as `xmax`.
- `face_labels(axes)` — `axis.py:369-381`. The canonical face inventory:
  `(axis ascending, endpoint in axis.endpoints order)`.
- `axes_from_legacy_mesh(mesh)` — `axis.py:493-580`. The `Mesh1D`/`Mesh2D`
  inbound adapter: `bc_left/bc_right` → `bc_low/bc_high`;
  spherical/cylindrical `bc_right` → `bc_outer` and `bc_left` is **dropped**
  (`axis.py:508-513` — the centreline is implicit-reflective, and the pole is
  not an endpoint); 2-D `bc_xmin/bc_xmax/bc_ymin/bc_ymax` → two `AxisMesh`.

So the user-facing declaration surfaces are: (a) `Mesh1D(..., bc_left=BC("vacuum"))`
etc., which get adapted; (b) `AxisMesh(edges=..., bc_low=BC(...), bc_high=BC(...))`
directly; (c) `axis.with_uniform_bc(BC(...))` — used by `solve_sn`'s
`boundary_condition` parameter.

### Level 2 — tag → law: the per-method admission registry + `RegistryMixin`

Two independent registries exist and must not be confused:

**(2a) `BoundaryTraceLaw.registry`** — the *self-registration* dict of
ALL law classes, `orpheus/geometry/boundary/_base.py:142`:

```python
registry: ClassVar[dict[str, type["BoundaryTraceLaw"]]] = {}

@classmethod
def _registry_base(cls) -> type:
    return BoundaryTraceLaw
```

Populated by `RegistryMixin.__init_subclass__` via the class-creation kwarg,
e.g. `class VacuumInflow(BoundaryTraceLaw, key="vacuum")`. Each law also gets a
`.key` attribute — read back at `augmented_mesh.py:427` to tag the SN shim.

**(2b) `<Method>.BOUNDARY_OPERATOR_REGISTRY`** — the per-method *admission
table*, a `ClassVar[dict[str, type[BoundaryTraceLaw]]]` declared on the
Protocol at `orpheus/transport/method.py:190`. This is what actually gates
which tags a method accepts:

| | SN (`orpheus/sn/mesh/augmented_mesh.py:171-174`) | Diffusion (`orpheus/diffusion/augmented_mesh.py:158-163`) |
|---|---|---|
| `"vacuum"` | `VacuumInflow` | `VacuumInflow` |
| `"reflective"` | `ReflectiveBoundary` | `ReflectiveBoundary` |
| `"albedo"` | — | `AlbedoBoundary` |
| `"zero_flux"` | — | `ZeroFluxBoundary` |

**SN admits only TWO tags.** `white`, `periodic`, `albedo`, `prescribed_inflow`
are realizable by `SNBoundaryRealizer` but NOT reachable from a `BC` tag —
explicitly noted at `augmented_mesh.py:182-187` ("adding them requires
SN-sweep-side wiring"). Diffusion deliberately omits `"white"`
(`diffusion/augmented_mesh.py:154-157`: at P1 white coincides with reflective).

### Level 3 — the shared face loop

`resolve_boundary_conditions(method)` — `orpheus/transport/method.py:218-262`.
Body (`method.py:254-262`), verbatim:

```python
default = BC("reflective")
resolved: dict[str, OpT] = {}
for label in face_labels(method.axes):
    tag = method.axes[label.axis_index].bc[label.endpoint] or default
    law = _law_from_tag(method, tag, label)
    resolved[label.face_name] = method.realize_boundary_law(
        law, label.face_name,
    )
return resolved
```

It **owns**: the face loop over `face_labels(method.axes)`; the
`BC("reflective")` default for an unset endpoint (the infinite-lattice /
eigenvalue convention); the tag → law parse.
It **delegates**: everything method-specific, through
`method.realize_boundary_law(law, face_name)`.

`_law_from_tag` — `method.py:265-311`. Three arms:
- unknown tag → `ValueError` naming the supported set (`method.py:289-298`);
- `ReflectiveBoundary` → `ReflectiveBoundary(axis=AXIS_NAMES[label.axis_index], albedo=1.0)`
  (`method.py:299-302`) — the reflection axis derives from the FACE, so a
  z-face gets the right permutation at any dimension;
- `AlbedoBoundary` → `AlbedoBoundary(albedo=float(tag.params["albedo"]))`, with
  a loud refusal for a parameter-less `BC("albedo")` (`method.py:303-310`);
- everything else → `law_cls()` (`method.py:311`).

**Note the hard-coded albedo=1.0 on the reflective arm**: a `BC("reflective")`
can never produce a partially-reflecting law through this path. The
`ReflectiveBoundary.albedo` field is only reachable by constructing the law
object directly.

`TransportMethod` Protocol itself — `method.py:145-215`; the two conformers are
`SNMesh` and `DiffusionMesh`, structurally (nobody imports the Protocol to
conform; `method.py:60-73`).

### Level 4 — the per-method hook `realize_boundary_law`

**SN** — `orpheus/sn/mesh/augmented_mesh.py:385-427`. Body:

```python
method_space = SNMethodSpace.for_face(
    mesh=self.mesh,
    quadrature=self.quad,
    face=face,
    trace=self._trace,
)
realized = SNBoundaryRealizer().realize(law, method_space)
return _BoundBoundaryOperator(realized, kind=law.key)
```

Called from the construction body at `augmented_mesh.py:351-353`
(`self.bc = resolve_boundary_conditions(self)`), immediately after the
`AngularTraceSpace` is built at `augmented_mesh.py:336-339`.

**Diffusion** — `orpheus/diffusion/augmented_mesh.py:294-315`. Body:

```python
method_space = DiffusionMethodSpace.for_face(mesh=self, face=face)
return DiffusionBoundaryRealizer().realize(law, method_space)
```

Called at `diffusion/augmented_mesh.py:247-249`, after the `ScalarTraceSpace`
is built at `diffusion/augmented_mesh.py:224-237`. Note the diffusion arm does
NOT wrap the result in any shim — `OpT` is the bare `LinearOperator`.

### Level 5 — method-space construction

**`SNMethodSpace`** — `orpheus/sn/mesh/method_space.py:71-175`. Frozen
dataclass, five fields (`method_space.py:92-96`):
`quadrature`, `face`, `inflow_indices`, `mesh`, `trace`.
`for_face` (`method_space.py:111-158`) populates
`inflow_indices = trace.inflow_indices_for_face(face)` (`method_space.py:150-151`)
— **the only field the realizer actually consumes beyond `quadrature`**.
`mesh` is explicitly documented as dead metadata: "nothing in the realizer
chain reads it" (`method_space.py:134-137`).
`minimal(quadrature)` (`method_space.py:98-109`) gives a quadrature-only space;
vacuum realization on it raises.

**`DiffusionMethodSpace`** — `orpheus/diffusion/method_space.py` (read below).

### Level 6 — the realizer `isinstance` chains

**`SNBoundaryRealizer.realize`** — `orpheus/sn/boundary/realizer.py:147-363`.
Pre-dispatch, in order:
1. `ZeroFluxBoundary` → **refused** with `BoundaryError`
   (`realizer.py:164-180`): "a negative angular inflow has no transport
   realization (ψ ≥ 0 ⟹ J± ≥ 0)".
2. `quad = method_space.quadrature` (`realizer.py:182`).
3. `law.assert_realizable(quad, inflow_indices=method_space.inflow_indices)`
   (`realizer.py:194-197`) — the §16A.12 invariant aggregate, fired once.

Then the seven-arm dispatch (full enumeration in §2 below).

**`DiffusionBoundaryRealizer.realize`** —
`orpheus/diffusion/boundary_realizer.py:176-195`. It is NOT a monolithic
isinstance chain: it is a two-stage composition —

```python
albedo = self._partial_current_albedo(law)
return stamp_boundary_role(_albedo_operator(albedo))
```

`_partial_current_albedo` (`boundary_realizer.py:197-249`) is the isinstance
chain law → scalar 𝒜; `_albedo_operator` (`boundary_realizer.py:142-157`) is
the structure-keyed scalar → operator collapse. **This is already a
factorization** — see §3.

### Level 7 — the realized operator classes: what `apply` actually computes

All in `orpheus/numerics/operator.py` unless noted.

**`PermutationOperator`** — `operator.py:2171-2281`. Body (`operator.py:2251-2255`):

```python
def apply(self, x: np.ndarray) -> np.ndarray:
    return np.take(x, self.perm, axis=self.axis)

def apply_transpose(self, x: np.ndarray) -> np.ndarray:
    return np.take(x, self.inverse_perm, axis=self.axis)
```

An integer gather along `axis`. `inverse_perm = np.argsort(perm)` at
construction; `is_involution` = `perm[perm] == arange(n)` (`operator.py:2247-2249`).

**`IncomingOrdinateMaskTensor`** — `operator.py:2284-2375`. Body
(`operator.py:2355-2365`):

```python
def apply(self, x: np.ndarray) -> np.ndarray:
    out = np.asarray(x).copy()
    if self.inflow_indices.size == 0:
        return out
    if self.axis == 0:
        out[self.inflow_indices] = 0.0
    else:
        idx: list = [slice(None)] * out.ndim
        idx[self.axis] = self.inflow_indices
        out[tuple(idx)] = 0.0
    return out
```

A copy-then-zero-the-inflow-rows projector. **Idempotent, self-adjoint,
rank-deficient**; distinct from `ZeroOperator` because the OUTFLOW trace
survives (`operator.py:2301-2312`).

**`PeriodicWrapOperator`** — `operator.py:2378-2417`. Body
(`operator.py:2407-2413`):

```python
def apply(self, x: np.ndarray) -> np.ndarray:
    return np.asarray(x).copy()

def apply_transpose(self, x: np.ndarray) -> np.ndarray:
    return np.asarray(x).copy()
```

**It is the identity-with-copy.** The docstring is explicit
(`operator.py:2382-2394`): the SN sweep does the spatial wrap through its own
face-pair indexing; the type is a placeholder "reserved for a future
spatial-pushforward extension … See follow-up issue". So today the periodic BC
operator computes NOTHING.

**`IdentityOperator`** — `operator.py:1834-1859`. `apply` returns `x` (no copy).
**`ZeroOperator`** — `operator.py:1862-1950`. `apply` returns `0.0 * x`
(endomorphism default) or `codomain_zero(x)` when a codomain zero-builder was
supplied.
**`ScaledOperator`** — `operator.py:1700-1831`. `apply` =
`self.scalar * self.op.apply(x)` (`operator.py:1758-1759`). Refuses `scalar == 0`
(`operator.py:1729-1737` — "use ZeroOperator explicitly"). It PROPAGATES
`block_role` and `system_role` from its operand (`operator.py:1740-1743`) —
which is why `stamp_boundary_role(α * base)` works even though the stamp lands
on the outer wrapper.

**`TensorProductOperator`** — `operator.py:2420-2558`. `apply` is a
LEFT-TO-RIGHT FOLD over the factors (`operator.py:2503-2507`):

```python
def apply(self, x: np.ndarray) -> np.ndarray:
    out = x
    for op in self.ops:
        out = op.apply(out)
    return out
```

Built by the `&` operator via `_build` which FLATTENS nested TPs
(`operator.py:2491-2501`). Factors must act on disjoint axes and broadcast on
the rest, so order is irrelevant. Note: because the fold is sequential apply,
`X & IdentityOperator()` is bit-identical to bare `X` — which is exactly the
"bit-identity preserved" claim the SN realizer's Wave-T comments make.

**`AngularAverageOperator`** — `orpheus/sn/boundary/angular.py:36-225`.
**The load-bearing body** (`angular.py:194-225`), verbatim:

```python
def apply(self, psi: np.ndarray) -> np.ndarray:
    if psi.shape[0] != self.n_ordinates:
        raise ValueError(
            f"AngularAverageOperator.apply: psi.shape[0] = "
            f"{psi.shape[0]}, expected {self.n_ordinates}"
        )
    # Cosine-weighted scalar average over the outgoing hemisphere.
    psi_avg = (
        self._cos_w.reshape((-1,) + (1,) * (psi.ndim - 1)) * psi
    ).sum(axis=0) / self._norm
    # Broadcast the (..., ) average back over every ordinate.
    return np.broadcast_to(
        psi_avg[None, ...], psi.shape,
    ).copy()
```

with, from `from_quadrature` (`angular.py:179-192`):

```python
weights = quadrature.weights
outgoing_mask = (outward_sign * mu_n) > 0.0
cos_w = weights * (outward_sign * mu_n)
cos_w = np.where(outgoing_mask, cos_w, 0.0)
norm = float(cos_w.sum())
```

**`IncomingSourceOperator`** — `orpheus/sn/boundary/angular.py:228-347`. Body
(`angular.py:323-347`):

```python
def apply(self, psi_out: np.ndarray) -> np.ndarray:
    q = self.source.evaluate(tuple(int(s) for s in psi_out.shape))
    if self._inflow_mask is None:
        return q
    if q.shape[0] != self._inflow_mask.shape[0]:
        raise ValueError(...)
    return q * self._inflow_mask.reshape(
        (-1,) + (1,) * (q.ndim - 1)
    )
```

Rank-0: the input is used ONLY for its shape. `_inflow_mask` is a dense 0/1
float vector built at construction (`angular.py:318-320`).

### Level 8 — `stamp_boundary_role` and the role attachment

`stamp_boundary_role` — `orpheus/geometry/boundary/_realizer.py:109-129`:

```python
def stamp_boundary_role(op: "LinearOperator") -> "LinearOperator":
    op.block_role = BlockRole.BOUNDARY
    return op
```

An **instance** stamp (not a class attribute) — deliberate: the realized op is
a generic numerics primitive that plays the BOUNDARY role only in this context
(`_realizer.py:114-122`).

- `BlockRole` — `operator.py:208-231`: `BULK` (`A_bb` only) / `FULL`
  (off-diagonal `A_bs`/`A_sb`; only `L`) / `BOUNDARY` (`A_ss` only).
- `SystemRole` — `operator.py:234-280`: `A` / `B` / `COUPLED` — the ORTHOGONAL
  curvilinear two-system axis. **Boundary realization never touches
  `system_role`**; it stays `None` on every realized BC leaf, except where a
  composer propagates one (`ScaledOperator.__init__`, `operator.py:1743`).
- The join lattices: `_join_block_roles` (`operator.py:283-314`) and
  `_join_system_roles` (`operator.py:317-335`) — "same role stays, anything
  different becomes FULL/COUPLED, `None` propagates". This is what makes
  `(L + C − S − F − B)` carry `FULL` by construction.
- `isinstance` markers `BulkOperator`/`FullOperator`/`BoundaryOperator` via
  `_BlockRoleMeta.__instancecheck__` reading `op.block_role`
  (`operator.py:338-387`).

**Who is stamped, per realizer:**

| law | SN stamped? | diffusion stamped? |
|---|---|---|
| vacuum | YES (`sn/boundary/realizer.py:256`) | YES (`boundary_realizer.py:195`) |
| reflective | YES (`realizer.py:287`, `:288`) | YES |
| white | YES (`realizer.py:304`, `:305`) | YES (P1-coincident with reflective) |
| albedo | YES (`realizer.py:309`, `:311`, `:319`) | YES |
| periodic | YES (`realizer.py:332`) | REFUSED |
| zero_flux | REFUSED | YES |
| **prescribed_inflow** | **NO** (`realizer.py:344`, `:349` — no stamp) | REFUSED |

The prescribed-inflow exclusion is deliberate and documented in three places
(`_realizer.py:123-126`, `operator.py:384-386`,
`diffusion/boundary_realizer.py:82-88`): it is the boundary *source*
`q.boundary`, not a linear boundary operator `B`.

### The SN `kind` shim — `_BoundBoundaryOperator`

`orpheus/geometry/boundary/_bound_compat.py:77+`. Wraps the realized operator
with a `kind` string so `sn_mesh.bc["xmin"] == "vacuum"` comparisons work
(`sn/mesh/augmented_mesh.py:160-163`, `:427`). Strict 1-arg passthrough since
#186/C-B3.4. Diffusion has NO such shim — `DiffusionMesh.bc` values are bare
`LinearOperator`s.

### Level 9 — the composite: the whole-trace `B`

**`SNBoundaryOperator`** — `orpheus/sn/operators/boundary.py:123-475`.
`block_role = BlockRole.BOUNDARY` is a CLASS attribute here (`boundary.py:175`),
not an instance stamp. `_face_laws` (`boundary.py:180-194`) is a `@property`
that re-reads `sn_mesh.bc` on every access:

```python
return {
    face: self.sn_mesh.bc[face]
    for face in self.sn_mesh.angular_trace.layout.faces
}
```

**The load-bearing walk is `_reflect_trace`** (`boundary.py:221-320`). The
per-face core (`boundary.py:294-320`):

```python
for face, law in face_laws.items():
    face_in = boundary.face_view(face)
    sel = (
        rows[face] if rows is not None
        else trace.inflow_indices_for_face(face)
    )
    if method == "apply":
        full = law.apply(face_in)
        out_boundary.face_view(face)[sel] = full[sel]
    else:
        # ``(P_sel ∘ law)ᵀ = lawᵀ ∘ P_sel``: mask the INPUT to the
        # forward's codomain rows, write the FULL transpose image.
        if not adjointable(law):
            raise MissingAdjoint(...)
        masked = np.zeros_like(face_in)
        masked[sel] = face_in[sel]
        out_boundary.face_view(face)[...] = law.apply_transpose(masked)
```

**This is architecturally the most important line in the whole chain.** The
realized law is applied to the WHOLE face slot (`face_in` — all N ordinates,
inflow AND outflow), and the composite THEN projects onto the inflow rows:
`B_face = P_inflow ∘ law`. The docstring says so explicitly
(`boundary.py:231-242`): "The realized per-face law … is a *full-face*
operator, so a non-zero outflow emission would corrupt the outflow-definition
residual". **The `G_in` half of the hypothesized factorization already exists —
it just lives in the composite, not in the law.**

The ⚠ note at `boundary.py:244-251` records a real bug class: output-projecting
`lawᵀ` onto the OUTFLOW rows extracts a law's DIAGONAL block, giving a spurious
`+1` outflow diagonal for the vacuum mask — caught only by the A2a
grid-reciprocity arm on the het-VACUUM sphere (permutation laws are
bit-identical under either spelling).

Other members: `_apply_faces` lifts onto a zero-bulk `FullField`
(`boundary.py:322-363`); `apply` (`:365-367`); `reflect_into_inflow`
(`:369-395`) — the trace-only entry the SI loop uses; `reflect_inflow_inplace`
(`:397-429`) — the mutating façade; `split(schedule)` (`:431-465`) →
`B_lower`/`B_upper` as `SNMaskedBoundaryOperator`s for the Gauss-Seidel regular
splitting; `apply_transpose` (`:467-475`); `is_adjointable` (`:196-205`) is the
per-face INTERSECTION — **white makes the composite non-adjointable** because
`AngularAverageOperator` advertises no transpose.

Siblings: `RadialCharacteristicBoundaryOperator` (`B_b`, System B's ray-corner
boundary, `boundary.py:478+`) and `SNMaskedBoundaryOperator` (`boundary.py:715+`).

**`DiffusionBoundaryOperator`** — `orpheus/diffusion/operators.py:544-668`
(**NOT `:407` as the brief has it — `:407` is `LeakageOperator.apply` on the
current tree; the file moved under it**). `face_laws` is a plain dict snapshot
taken in `__init__` (`operators.py:580`), unlike SN's live property. The walk
(`operators.py:610-618`):

```python
out_boundary = ScalarBoundarySourceSink.zeros_on(self.mesh)
for face, law in self.face_laws.items():
    out_boundary.face_view(face)[ScalarTraceSpace.INFLOW_ROW] = (
        law.apply(trace.outflow_view(face))
    )
return FullField(
    interior=ScalarSourceSink.zeros_on(self.mesh),
    boundary=out_boundary,
)
```

Here the row selection is STRUCTURAL, not quadrature-derived: read
`OUTFLOW_ROW`, write `INFLOW_ROW`. It also has an assembly mode
(`operators.py:626-668`) that extracts each face's `(ng, ng)` block through
`law.as_matrix(basis_shape=(ng,))` — ng probes of `law.apply` — never re-reading
the albedo scalar.

### Level 10 — the actual numerics, two concrete cases

#### (A) S4 slab, reflective left face (`BC("reflective")` on `xmin`)

Measured on this tree with `Quadrature.gauss_legendre(4)`:

```
mu_x           = [-0.861136, -0.339981, +0.339981, +0.861136]
weights        = [ 0.347855,  0.652145,  0.652145,  0.347855]
reflection_index("x") = [3, 2, 1, 0]
omega_dot_n["xmin"]   = [+0.861136, +0.339981, -0.339981, -0.861136]
inflow  at xmin = [2, 3]      outflow at xmin = [0, 1]
```

Chain: `AxisMesh(bc_low=BC("reflective"))` → `face_labels` yields
`FaceLabel(0,"min")` → `face_name == "xmin"` → `_law_from_tag` builds
`ReflectiveBoundary(axis="x", albedo=1.0)` → `SNMesh.realize_boundary_law`
builds `SNMethodSpace.for_face(...)` with `inflow_indices = [2,3]` →
`SNBoundaryRealizer.realize` takes the `albedo == 1.0` fast path
(`sn/boundary/realizer.py:282-287`) and returns

```python
stamp_boundary_role(PermutationOperator([3,2,1,0], axis=0) & IdentityOperator())
```

wrapped as `_BoundBoundaryOperator(..., kind="reflective")`.

**The array op.** With a face slot `psi_face` of shape `(4, ng)`:

1. `TensorProductOperator.apply` folds: `IdentityOperator.apply(np.take(psi_face, [3,2,1,0], axis=0))`
   → `full = psi_face[[3,2,1,0]]`, i.e. `full[n] = psi_face[3-n]`.
2. `SNBoundaryOperator._reflect_trace` writes only `sel = [2,3]`:
   `out[2] = full[2] = psi_face[1]`, `out[3] = full[3] = psi_face[0]`.
   `out[0] = out[1] = 0`.

So the net numerical action is **two scatter-writes**: inflow ordinate 2
(µ=−0.34) receives outflow ordinate 1 (µ=+0.34); inflow ordinate 3 (µ=−0.86)
receives outflow ordinate 0… wait — careful: `reflection_index` is
`[3,2,1,0]` so `full[3] = psi_face[0]` (µ=−0.861 ← µ=−0.861?). No:
`psi_face[0]` is the ordinate with µ=−0.861136 which is OUTFLOW at xmin
(`omega_dot_n[0] = +0.861` > 0). Consistent: at `xmin` the outward normal is
−x̂, so µ<0 ordinates are outgoing. Inflow rows `[2,3]` are µ>0, receiving from
`[1,0]` = µ<0 — the specular partner. The **half of the computed permutation
that lands on outflow rows is computed and then discarded.**

#### (B) Diffusion 1-D, albedo face (`BC("albedo", {"albedo": 0.35})`)

Chain: `_law_from_tag` → `AlbedoBoundary(albedo=0.35)` (`method.py:303-310`) →
`DiffusionMesh.realize_boundary_law` → `DiffusionBoundaryRealizer.realize`:
`_partial_current_albedo` returns `0.35` (`boundary_realizer.py:215-216`),
`_albedo_operator(0.35)` returns `0.35 * IdentityOperator()`
(`boundary_realizer.py:157`) = `ScaledOperator(0.35, IdentityOperator())`,
stamped BOUNDARY.

**The array op.** With a face slot of shape `(2, ng)`
(rows = `OUTFLOW_ROW`, `INFLOW_ROW`):

`DiffusionBoundaryOperator.apply` does
`out[INFLOW_ROW] = 0.35 * trace.outflow_view(face)`, i.e. a scalar multiply of
a `(ng,)` vector: `J⁻_g = 0.35 · J⁺_g` for every group `g`.
The assembled form is `ng` diagonal entries `0.35` at
`(j_minus[g], j_plus[g])` in the global sparse matrix.

---

## 2. Per-law table: law → realizer arm → operator → array op

`BoundaryTraceLaw.registry` (verified at runtime on this tree) holds SEVEN
laws: `albedo`, `periodic`, `prescribed_inflow`, `reflective`, `vacuum`,
`white`, `zero_flux`.

### 2a. SN (`SNBoundaryRealizer.realize`, `orpheus/sn/boundary/realizer.py:147-363`)

| Law | Arm | Construction expression (verbatim) | The array operation |
|---|---|---|---|
| `ZeroFluxBoundary` | `:164-180` | — | **REFUSED** `BoundaryError`: 𝒜=−1 has no non-negative angular realization |
| `VacuumInflow` | `:199-263` | `stamp_boundary_role(IncomingOrdinateMaskTensor(inflow_indices=ms.inflow_indices, n_ordinates=quad.N, axis=0) & IdentityOperator())` | `out = x.copy(); out[inflow_indices] = 0.0` — a projector onto the outflow subspace |
| `ReflectiveBoundary(α=1)` | `:265-287` | `stamp_boundary_role(PermutationOperator(perm, axis=0) & IdentityOperator())` | `np.take(x, perm, axis=0)` — an integer gather |
| `ReflectiveBoundary(α≠1)` | `:288` | `stamp_boundary_role(float(law.albedo) * base)` | `α * np.take(x, perm, axis=0)` |
| `WhiteBoundary(α=1)` | `:290-304` | `stamp_boundary_role(AngularAverageOperator.from_quadrature(quad, law.axis, law.outward_sign) & IdentityOperator())` | `broadcast(Σ_n cos_w[n]·x[n] / Σ cos_w)` — **a rank-1 contract-then-broadcast** |
| `WhiteBoundary(α≠1)` | `:305` | `stamp_boundary_role(float(law.albedo) * base)` | `α ·` the above |
| `AlbedoBoundary(0.0)` | `:307-309` | `stamp_boundary_role(ZeroOperator())` | `0.0 * x` |
| `AlbedoBoundary(1.0)` | `:310-311` | `stamp_boundary_role(IdentityOperator())` | `x` (no copy) |
| `AlbedoBoundary(α)` | `:312-321` | `stamp_boundary_role(float(law.albedo) * (IdentityOperator() & IdentityOperator()))` | `α * x` |
| `PeriodicBoundary` | `:323-332` | `stamp_boundary_role(PeriodicWrapOperator() & IdentityOperator())` | `np.asarray(x).copy()` — **the identity** |
| `PrescribedInflow` | `:334-349` | `IncomingSourceOperator(law.source, inflow_indices=…, n_ordinates=…)` — **NOT stamped** | `source.evaluate(shape) * inflow_mask` — input ignored |
| anything else | `:351-363` | — | `BoundaryError` naming the six available cases |

### 2b. Diffusion (`DiffusionBoundaryRealizer`, `orpheus/diffusion/boundary_realizer.py`)

Two-stage: `_partial_current_albedo(law) -> float` (`:197-249`) then
`_albedo_operator(albedo)` (`:142-157`), then `stamp_boundary_role` (`:195`).

| Law | 𝒜 (arm) | Operator | The array operation |
|---|---|---|---|
| `VacuumInflow` | `0.0` (`:209-210`, Marshak `J⁻=0`) | `ZeroOperator()` | `0.0 * J⁺` |
| `ZeroFluxBoundary` | `-1.0` (`:211-212`) | `float(-1.0) * IdentityOperator()` = `ScaledOperator(-1, I)` | `-1.0 * J⁺` |
| `ReflectiveBoundary(α)` | `float(law.albedo)` (`:213-214`) | `I` (α=1) / `α·I` | `α * J⁺` |
| `WhiteBoundary(α)` | `float(law.albedo)` (`:213-214`) | same as reflective — **P1-coincident** | `α * J⁺` |
| `AlbedoBoundary(α)` | `float(law.albedo)` (`:215-216`) | `Zero`/`I`/`α·I` | `α * J⁺` |
| `PeriodicBoundary` | — (`:218-227`) | **REFUSED** | needs the opposite-face wrap (a trace-block permutation) |
| `PrescribedInflow` | — (`:228-238`) | **REFUSED** | `J⁻=q` is the affine source, not an operator |

Note the α=0 collapse asymmetry: `_albedo_operator(0.0)` returns `ZeroOperator`
because `ScaledOperator` REFUSES a zero scalar (`operator.py:1729-1737`).

### 2c. The two REFUSAL asymmetries (they are exactly complementary)

- `zero_flux`: diffusion-only. SN refuses (negative angular inflow).
- `prescribed_inflow`: SN realizes as a NON-stamped affine source; diffusion
  refuses pending the P5 fixed-source wiring.
- `periodic`: SN realizes (as a no-op); diffusion refuses.
- `white`: SN realizes a genuine angular average; diffusion collapses it onto
  reflective — and deliberately does not ADMIT the tag.

---

## 3. The factorability analysis

**Notation.** I use PIPELINE order throughout, because the hypothesis is stated
that way: `ψ_out --[G_out]--> (reduced) --[R]--> (reduced) --[G_in]--> ψ_in`.
As operator composition that is `B = G_in ∘ R ∘ G_out`. (Note the brief writes
`B = G_out ∘ R ∘ G_in`; the historical source writes `γ₋ψ = R G γ₊ψ`, i.e.
`B = R ∘ G` with `G` applied first — see §4.)

### (a) `AngularAverageOperator` — ALREADY a contract-then-broadcast rank-one form

**Answer: yes, internally it is exactly a rank-one dyad — but written inline,
not composed from `RankOneOperator`.**

The body (`orpheus/sn/boundary/angular.py:218-225`) is two statements:

```python
        # Cosine-weighted scalar average over the outgoing hemisphere.
        psi_avg = (
            self._cos_w.reshape((-1,) + (1,) * (psi.ndim - 1)) * psi
        ).sum(axis=0) / self._norm
        # Broadcast the (..., ) average back over every ordinate.
        return np.broadcast_to(
            psi_avg[None, ...], psi.shape,
        ).copy()
```

Statement 1 IS the contraction `⟨w, ψ⟩` with `w = cos_w / norm`; statement 2 IS
the broadcast `𝟙 · ⟨w, ψ⟩`. The dyad is

```
AngularAverage  =  |𝟙⟩⟨ cos_w / norm |     with   cos_w[n] = 1{Ω_n·n̂ > 0} · w_n · |Ω_n·n̂|
```

Note `cos_w` is **already zeroed on the non-outgoing hemisphere**
(`angular.py:181-184`): `cos_w = np.where(outgoing_mask, weights*(sign*mu_n), 0.0)`.
**So `G_out` (the outgoing-hemisphere restriction) is FOLDED INTO the covector.**

I verified numerically on this tree, S8 Gauss-Legendre, `psi` shape `(8,3,5)`:

```
outer(np.ones((8,1,1)), InnerProductFunctional(cos_w.reshape(8,1,1)/norm, axis=0))
   vs  AngularAverageOperator.from_quadrature(q, "x", +1)
max abs diff = 2.22e-16     bit identical: False
```

Not bit-identical only because of division ORDER (the operator divides the
SUM by `norm`; the dyad divides the WEIGHT). It is the same rank-1 map to
round-off.

**Three-factor reading of the white BC** (the hypothesis's worked example, and
it holds):

| factor | what it is | where it lives today |
|---|---|---|
| `G_out` | restrict to the outgoing hemisphere, contract with the partial-current covector `w·\|Ω·n̂\|` | **folded into `cos_w`** (`angular.py:181-184`) |
| `R` | `1/norm` (conservation normalization) × `α` (the albedo) | **`1/norm` folded into `apply`; `α` correctly hoisted out into `ScaledOperator` by the realizer (`realizer.py:305`)** |
| `G_in` | isotropic rebroadcast over ordinates, then restrict to the inflow set | **broadcast folded into `apply`; the inflow restriction lives in `SNBoundaryOperator._reflect_trace`** (`sn/operators/boundary.py:302`) |

So all three factors EXIST — spread across three files, one of them the
composite operator, none of them a first-class object.

### (b) Every other realized operator, factored

| Law | `G_out` | `R` | `G_in` | Factorable as a tensor network? |
|---|---|---|---|---|
| **reflective** | — (the permutation is not a contraction; it is a bijection on the whole face) | `α·I` (diagonal, hoisted to `ScaledOperator` when α≠1) | `P_inflow` (in `_reflect_trace`) | **Degenerate.** `Π` is full-rank, so there is no rank reduction to contract to. The honest form is `B = P_in · (αI) · Π` — and R already IS separated. The sandwich buys nothing here. |
| **vacuum** | `Mask` zeroes the inflow rows, PRESERVING the outflow | `0` implicitly | `P_inflow` | **Pathological — see below.** |
| **albedo** | `I` | `α` | `P_inflow` | **Perfect degenerate fit**: `G_out = G_in = I`, all content in `R`. This is literally the grand report's `IdentityBoundaryMap` pattern. |
| **periodic** | *the spatial wrap — but it is NOT in the operator* | `1` | `P_inflow` | **The factor is MISSING.** `PeriodicWrapOperator.apply` is `np.asarray(x).copy()` (`operator.py:2407-2410`); the actual geometry (opposite-face pairing) lives in the SN sweep's face-pair indexing. The type is a named placeholder for a `G` that was never built. |
| **prescribed inflow** | — (rank-0, input ignored) | `0` | `_inflow_mask` **inside the operator** (`angular.py:345-347`) | Affine: `B=0`, `q = mask ⊙ source.evaluate(shape)`. Its `G_in` is a SECOND, independent copy of the inflow mask that `_reflect_trace` also applies. |
| **zero-flux (diffusion)** | read `OUTFLOW_ROW` | `−1` | write `INFLOW_ROW` | **Already exactly three factors** — see below. |
| **diffusion, all laws** | `trace.outflow_view(face)` | `𝒜·I` | `out.face_view(face)[INFLOW_ROW] = …` | **The tensor network already exists**, degenerately: the trace is 1 DOF per (face, group), so `G` is a structural row selection, and the realizer's ONLY job is to produce the scalar `𝒜`. `DiffusionBoundaryRealizer` is literally `law → 𝒜 → operator` in two named stages (`boundary_realizer.py:194-195`). |

**Two findings that fall out of doing this exercise:**

**(i) The vacuum arm's composite action is the ZERO MAP, and the whole
outflow-preserving design is discarded.** `IncomingOrdinateMaskTensor` was
written to zero ONLY the inflow ordinates so "the outflow trace survives …
downstream consumers (sensitivity adjoints, post-BC field readers) require [it]"
(`operator.py:2301-2306`; `sn/boundary/realizer.py:23-28` calls this the §16A.5
"semantic correction"). But its one production consumer,
`SNBoundaryOperator._reflect_trace`, writes only `out[inflow] = full[inflow]` —
and `full[inflow]` is exactly what the mask just set to zero. Verified on this
tree:

```
face_in      = [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]   inflow=[2,3]
law.apply    = [[1,2,3],[4,5,6],[0,0,0],[0,0,0]]
out[inflow]  = full[inflow]  ->  [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
IS THE ZERO MAP: True
```

So the SN vacuum BC computes a full face copy, zeroes half of it, and then
copies the zeroed half into a zero buffer. Under the composite it is
`ZeroOperator` — which is exactly what the DIFFUSION realizer produces for the
same law (`boundary_realizer.py:209-210`). The mask's extra structure is
unreachable through `B`. (It is not *wrong* — the outflow-preserving semantics
matter if a consumer ever applies the law to a face slot directly — but no such
consumer exists today, and `_reflect_trace`'s ⚠ note at
`sn/operators/boundary.py:244-251` records that the transpose spelling that
WOULD have exposed the difference was a bug.)

**(ii) `P_inflow` — the `G_in` factor — is implemented THREE times, in three
different places, in three different spellings.**
1. `SNBoundaryOperator._reflect_trace`: `out[sel] = full[sel]` with
   `sel = trace.inflow_indices_for_face(face)` (`sn/operators/boundary.py:296-302`).
2. `IncomingOrdinateMaskTensor`: a `[indices] = 0.0` scatter (the COMPLEMENT
   projector) — `operator.py:2355-2365`.
3. `IncomingSourceOperator._inflow_mask`: a dense 0/1 float multiply
   (`angular.py:318-320`, `:345-347`).
Three encodings of one geometric object. That is the strongest single piece of
evidence FOR the hypothesis: `G_in` is a real, load-bearing concept that has no
type, so it got re-derived at every site that needed it.

### (c) The primitives a tensor network would compose from — all present, none used at the boundary

| Primitive | `operator.py` | Surface | Boundary use? |
|---|---|---|---|
| `RankOneOperator` | `:2898-3056` | `apply = recon * functional.evaluate(x)`; `apply_transpose` = the dual dyad `\|w⟩⟨v\|`; `is_adjointable` iff the row is an `InnerProductFunctional`; **structurally non-invertible** | **NO.** Only production consumer is `FissionOperator.kernel` (`transport/operators/fission.py:311`) |
| `outer(recon, functional)` | `:3059-3083` | the readable constructor | NO |
| `InnerProductFunctional` | `numerics/functional.py:182-252` | `evaluate = (w*x).sum(axis, keepdims=True)` | NO |
| `TensorProductOperator` | `:2420-2558` | `apply` = sequential per-factor fold; `&` builds+flattens; adjoint/inverse distribute factor-wise | **YES but decoratively** — 5 sites in `sn/boundary/realizer.py`, EVERY one of the form `X & IdentityOperator()` (or `I & I`). The second factor is always a no-op; the comments say so ("the TP fold reduces to the bare mask's `apply`", `realizer.py:252-255`). It buys type-visibility, not structure. |
| `SumOfTensorProductsOperator` | `:2561-2652` | `apply` = summand-wise sum; `assert_separable()` | **ZERO production consumers anywhere in `orpheus/`.** Only its definition, the `numerics/__init__` re-export, and a docstring mention. |
| `PermutationOperator` | `:2171-2281` | integer gather, algebra-closed inverse | YES (reflective) |
| `DiagonalOperator` | `:2655+` | `np.expand_dims(c, broadcast_axes) * x`; self-adjoint; invertible iff no zero entry | NO (boundary uses `ScaledOperator` + `IdentityOperator` instead of a diagonal) |
| `OperatorProduct` | `:1504+` | the `@` composition composer | NO — **no boundary operator is ever built by composition.** |

**The decisive observation for (c):** the boundary layer already imports and
uses `TensorProductOperator` — the ⊗ composer — five times, and NEVER uses
`OperatorProduct` (∘) or `RankOneOperator` (the contraction⊗broadcast). Those
are precisely the two composers a `G_in ∘ R ∘ G_out` network needs. The
machinery is complete and sitting unused.

### Verdict on the hypothesis

**Supported, with one important correction.**

- **Confirmed:** `geometry_map` and `response_kernel` are declared on
  `BoundaryTraceLaw` (`_base.py:148-164`) and are `None` on **all seven**
  concrete laws (verified at runtime). Nothing populates them; nothing reads
  them. The only references anywhere are the declarations, two docstrings, and
  three tests that assert the defaults are `None`
  (`tests/geometry/test_boundary_trace_law.py:95-101`). The realizers build the
  operator monolithically via `isinstance`.
- **Confirmed:** the white BC IS a rank-one contract-then-broadcast, and it is
  written inline rather than composed.
- **Confirmed:** the three factors genuinely exist in the data flow — but they
  are distributed across the law, the realizer, the operator body, and the
  composite, with `G_in` triplicated.
- **Correction:** the historical intent was a **TWO**-factor `R ∘ G`, not a
  three-factor sandwich (§4). The third factor (`G_in`, the inflow
  restriction) emerged later, from the composite's `A_ss` block structure, and
  was never folded back into the law algebra.
- **Counter-evidence worth weighing:** for `reflective` and `albedo` — the only
  two laws SN actually ADMITS from a `BC` tag — the factorization is degenerate
  (`G` is a bijection or the identity, `R` is a scalar already hoisted into
  `ScaledOperator`). The network would pay off for `white`, `periodic`
  (genuinely missing its `G`), and any future anisotropic-albedo / Marshak /
  energy-redistributing law — none of which is tag-reachable today.

---

## 4. The historical intent — what `R` and `G` were meant to be, in the source's own words

### 4.1 The primary source: Grand Report v3 §16A

`.claude/plans/neutron_transport_grand_report_v3.md`. This is where
`geometry_map` / `response_kernel` come from, and it is unambiguous.

**§16A.1, the vocabulary list** (`:2742-2753`) — note these are named as
SEPARATE first-class types:

```text
BoundaryGeometryMap       # geometry-induced map on trace points
BoundaryResponseKernel    # physical response on mapped traces
BoundarySource            # prescribed incoming source
BoundaryTraceLaw          # affine law gamma_- = R G gamma_+ + q
BoundaryRealizer          # method-specific realization
MethodBoundaryOperator    # concrete operator/constraint/kernel/tensor for one method
```

**§16A.2, the intended skeleton** (`:2780-2802`) — verbatim:

```python
class BoundaryGeometryMap(Protocol):
    def map_outgoing_to_incoming_geometry(self, outgoing_trace_point): ...
    def as_tensor(self, trace_space): ...

class BoundaryResponseKernel(Protocol):
    def apply(self, mapped_outgoing_trace): ...
    def as_tensor(self, trace_space): ...

class BoundarySource(Protocol):
    def evaluate(self, incoming_trace_point): ...

@dataclass(frozen=True, slots=True)
class BoundaryTraceLaw:
    geometry_map: BoundaryGeometryMap
    response: BoundaryResponseKernel
    source: BoundarySource

    def as_operator(self, trace_space):
        G = self.geometry_map.as_operator(trace_space)
        R = self.response.as_operator(trace_space)
        q = self.source.as_field(trace_space.incoming)
        return AffineBoundaryOperator(linear=R @ G, source=q)
```

**And vacuum, spelled as a composition of the three** (`:2806-2815`):

```python
@dataclass(frozen=True, slots=True)
class ZeroInflowTrace(BoundaryTraceLaw):
    def __init__(self):
        super().__init__(
            geometry_map=IdentityBoundaryMap(),
            response=ZeroBoundaryResponse(),
            source=ZeroBoundarySource(),
        )
```

**§16A.5, the SN realizer sketch** (`:2946-2960`) — verbatim:

```python
@dataclass(frozen=True, slots=True)
class DiscreteOrdinatesBoundaryRealizer:
    def realize(self, law: BoundaryTraceLaw, space: DiscreteOrdinatesPhaseSpace):
        trace = DiscreteOrdinatesBoundaryTrace(space)

        if isinstance(law.response, ZeroBoundaryResponse):
            return ZeroInflowOrdinateConstraint(
                incoming_mask=trace.incoming_mask()
            )

        G = trace.discretize_geometry_map(law.geometry_map)
        R = trace.discretize_response(law.response)
        q = trace.discretize_source(law.source)
        return SNBoundaryOperator(linear=R @ G, source=q)
```

**This is the decisive citation.** The intended realizer dispatches on
`isinstance(law.response, ZeroBoundaryResponse)` — on the RESPONSE OBJECT, not
on the law class — and otherwise *discretizes each declared factor separately*
and composes `R @ G`. The shipped `SNBoundaryRealizer` does the opposite: it
dispatches on the LAW CLASS and builds the whole operator per arm.

**§16A.10, "Boundary tensor networks by method"** (`:3142-3166`) — verbatim:

```text
The tensor-network factorization should separate geometry and response:

outgoing trace
    -> BoundaryGeometryTensor
    -> BoundaryResponseTensor
    -> incoming trace

For `SN`:

BoundaryGeometryTensor        -> sparse face/direction permutation or mask
BoundaryResponseTensor        -> angular/energy albedo matrix
IncomingOrdinateMaskTensor    -> zero-inflow enforcement
```

**§16.2, the albedo tensor product** (`:2591-2614`) — verbatim:

```text
B = G_patch ⊗ K_omega ⊗ K_g
```
```
- `G_patch` maps boundary geometry
- `K_omega` redistributes angle
- `K_g` redistributes energy
```
```python
B = BoundaryTensorNetwork([
    PatchIncidenceFactor(patch),
    BoundaryGeometryMap(patch),
    AlbedoAngularKernel(...),
    AlbedoEnergyKernel(...),
])
```

followed by the stated payoff, verbatim:

> This makes it possible to swap geometry without changing response, or swap
> response without changing geometry.

**§16.3, the intended composition algebra** (`:2616-2638`) — verbatim:

```python
bc = BoundaryResponseKernel(...) @ BoundaryGeometryMap(...)   # layered responses
bc = PatchMap & AngularAlbedo & EnergyAlbedo                  # separable response
```

So the vision has **BOTH** composers: `@` (∘, sequential `R` after `G`) AND `&`
(⊗, per-axis separability). The shipped boundary layer uses ONLY `&`.

### 4.2 The theory pages

**`docs/theory/foundations/boundary_conditions.rst:183-227`** — the affine form
as the page of record. `G` is described as:

> :math:`G : \Gamma_+ \to \Gamma_+` is the **geometric map** — a
> measure-preserving permutation, pushforward, spatial wrap-around, or
> hemispheric cosine-weighted average. It carries pure geometry (it changes
> nothing about the physical interaction at the boundary; it just relabels the
> angular fluxes that meet there).

and `R`:

> :math:`R : \Gamma_+ \to \Gamma_-` is the **response kernel** — a scalar
> amplitude in :math:`[0, 1]` for the standard sub-Markov BCs (albedo, white,
> partial-current) or a full angular kernel in general weak-form BCs
> (deferred; …).

Note `G : Γ₊ → Γ₊` and `R : Γ₊ → Γ₋`. **So in the page's own typing, `R` is the
factor that crosses from outflow to inflow** — the inflow restriction is `R`'s
codomain, not a third factor. That is the two-factor reading.

**`boundary_conditions.rst:245-256`** states the intended payoff of the split:

> The split lets cross-method realizers introspect the law's geometric and
> response components separately — the SN realizer dispatches on the law's
> class today, but a future weak-form realizer might dispatch on the geometry /
> response / source independently.

**This sentence is the codebase admitting, in the theory page of record, that
the properties are declared-but-unused and naming the intended future
consumer.** It is the single best piece of "what was the original vision"
evidence in the docs.

**`docs/theory/methods/sn/boundary_conditions.rst:154-218`** — the
`bc-tensor-decomposition` equation (referenced from `sn/boundary/realizer.py:137`
and `geometry/boundary/vacuum.py:25`):

```
ψ_in(Ω) = (R ψ_out)(Ω) = Σ_α (G_α ψ_out)(Ω) · A_α
```

with a per-law table of `G_α` / `A_α` — attributed to Lewis & Miller (1984)
§3.4. **This table is prose only.** Its rows are exactly the values that
`geometry_map` / `response_kernel` were supposed to carry:

| Class | `G_α` | `A_α` |
|---|---|---|
| `VacuumInflow` | `0` (no operator) | 0 |
| `ReflectiveBoundary` | permutation under reflection axis | albedo |
| `WhiteBoundary` | cosine-weighted hemispheric average | albedo |
| `PeriodicBoundary` | spatial pushforward (caller-supplied) | 1 |
| `AlbedoBoundary` | identity in angle | albedo |
| `PrescribedInflow` | 0 | 0 |

**`docs/theory/foundations/operator_tensor_network.rst`** — the Wave-T
retrospective, and the most honest document in the set. `:73-83`:

> What landed is **not** the uniform :math:`A = \sum_k A_x^{(k)} \otimes
> A_\omega^{(k)} \otimes A_g^{(k)}` aspiration of Grand Report v3 §15-§16A.
> The shipped state is **five algebraically distinct shapes**, chosen per
> operator by what the underlying physics actually couples.

`:114-118`:

> **Zero production consumers** of :class:`SumOfTensorProductsOperator`. …
> Only T.1 (BC realizers) and T.2 (fission rank-1) cleanly admit the clean
> tensor-product factorisation.

and the BC row of the shape catalogue (`:172-185`) says the shipped form is
`§16A.10 B = G_patch ⊗ K_omega ⊗ K_g` **"with two factors degenerate to
`IdentityOperator`"**.

### 4.3 The git record — the near-miss

The decisive archaeological fact: **`RankOneOperator` was invented in the SAME
campaign, one step AFTER the BC realizers were lifted, and was not applied to
the white BC.**

- `fa13e78e` (2026-05-30 16:06) — *"Wave T step T.1 — lift remaining BC
  realizers to TensorProductOperator"*. Its own table records white as
  `AngularAverageOperator.from_quadrature(.)` → `… & I (TP base)`. The
  operator body was not touched.
- `0b2848be` (2026-05-30 16:24, **18 minutes later**) — *"Wave T step T.2 —
  fission as rank-1 TensorProductOperator"*. Message, verbatim:

  > The §15 form is mathematically a rank-1 outer product on the group axis
  > (``F[g,g'] = χ[g]·νΣ_f[g']``), not a disjoint-axes tensor product — both χ
  > and νΣ_f touch the group axis (χ expands, νΣ_f contracts). Introduces a
  > new L1 primitive `RankOneOperator(left, right, axis)` for the
  > ``|left⟩⟨right|`` structure…

  T.2 diagnosed *contraction-then-broadcast on a shared axis* as needing a new
  primitive — and minted it. T.1 had just walked past the identical structure
  in `AngularAverageOperator` (contraction-then-broadcast on the ORDINATE
  axis) and wrapped it in a decorative `& I` instead.

Other relevant history: `f71a32cb` (Issue #186) made `BoundaryTraceLaw` a pure
descriptor and removed `apply`; `fda1a59e` (Wave 3 / C3.1) is where
`geometry_map` / `response_kernel` first appeared — as `None`-returning
properties. `git log -S "response_kernel" -- orpheus/` returns **four commits,
none of which populates them**.

### 4.4 The live campaign is about to build the same shape — elsewhere

`.claude/plans/operator_strategy_realization_campaign.md` (the plan of record
for the branch this tree is on) lists, in its dependency-readiness table
(`:317`):

| capability | state | evidence |
|---|---|---|
| Realizer seam precedent | **LANDED** | `BoundaryTraceLaw → SNBoundaryRealizer → fixed-signature operator` |

and phase P2b (`:486-490`) proposes:

> **P2b — `SNReactionRealizer`**, on the `SNBoundaryRealizer` precedent. Each
> method builds its own monomorphic realized operator as a composition
> `j_bulk ∘ E₀ ∘ F_G ∘ M₀ ∘ π_bulk`. **Composition is the algebra of record and
> the equivalence oracle; a fused path may follow, pinned bit-identical.**

That target shape — `inject ∘ response ∘ project` — **is structurally the same
sandwich the hypothesis proposes for the boundary**, being planned for the
reaction operators, citing the boundary realizer as the precedent. The boundary
realizer is being cited as the model for a composition it does not itself
perform.

---

## 5. What surprised me

Ordered by how much I think it should change a decision.

### S1. The SN vacuum BC's composite action is literally the zero map

Covered in §3(b)(i). The `IncomingOrdinateMaskTensor`'s entire reason for
existing — the §16A.5 "semantic correction" that made vacuum preserve the
outflow trace instead of zeroing everything (`sn/boundary/realizer.py:19-28`;
`operator.py:2301-2306`) — is annulled by the one production consumer, which
keeps only `full[inflow]`. Under `B`, SN vacuum ≡ `ZeroOperator`, which is
exactly what diffusion produces. Verified numerically on this tree.
The stated beneficiary ("sensitivity adjoints … need the outflow trace
preserved") does not exist as a caller.

### S2. `G_in` (the inflow restriction) is implemented three times, three ways

`_reflect_trace`'s `out[sel] = full[sel]`; `IncomingOrdinateMaskTensor`'s
complement scatter; `IncomingSourceOperator._inflow_mask`'s dense 0/1 multiply.
Three encodings of one geometric object, none of them named. This is the
`Pattern 2 / single-source-of-truth` smell in its textbook form, and it is the
strongest concrete argument that the missing factor is real rather than
theoretical.

### S3. The `kind` STRING round-trips through the whole three-layer chain and is
### read back out for dispatch in four production sites

`BC.kind` (str) → `BOUNDARY_OPERATOR_REGISTRY` lookup → law class → `law.key`
(str) → `_BoundBoundaryOperator(kind=…)` → and then:

- `orpheus/sn/acceleration/dsa.py:213-215` — `str(sn_mesh.bc["xmin"].kind)`,
  checked against `_SUPPORTED_BC` to gate the DSA low-order rows.
- `orpheus/sn/loss_representation/sweep_schedule.py:265-268` —
  `if sn_mesh.bc[face] == "reflective"` (via `_BoundBoundaryOperator.__eq__`
  against a string, `_bound_compat.py:172-175`) to pick the specular faces the
  Gauss-Seidel schedule can forward-substitute through.
- `orpheus/sn/operators/boundary.py:556-559` and `:610-612` —
  `getattr(self.sn_mesh.bc["xmax"], "kind", None) in _RULED_CORNER_KINDS`, both
  for `is_adjointable` and for the ray-corner action.
- `orpheus/sn/solver.py:1362-1364` — `if op.kind == "vacuum"` to find the
  leakage faces.

The whole point of the three-layer split was to replace tag dispatch with type
dispatch, and the #290 P7b registry dissolution deleted a string-keyed registry
for exactly this reason (`transport/method.py:74-92`). Four call sites still
discriminate on the string, on the far side of the realization. Three of them
(`dsa`, `sweep_schedule`, ray-corner) are asking questions the type system
could answer — "does this law couple ordinates?", "is this law a specular
permutation?", "does this face leak?" — which are exactly `geometry_map` /
`response_kernel` questions. **This is the latent consumer the theory page
predicted, already present, answering itself with strings.**

### S4. `RankOneOperator` was minted 18 minutes after the BC lift, for the same structure

§4.3. The T.2 commit message diagnoses contraction-then-broadcast on a shared
axis as needing a dedicated primitive and creates one; T.1, immediately prior,
had wrapped the identical structure (`AngularAverageOperator`) in `& I` and
moved on. Fission got `outer(χ, ⟨νΣf|)`; white got a decorative tensor product.

### S5. `SumOfTensorProductsOperator` has zero production consumers, and the docs know

`operator_tensor_network.rst:114-118` states it outright. The type, its
`assert_separable` contract, and its `__init__` re-export exist for a form the
physics rejected. (It is not dead-by-accident — it is a documented aspiration.
But under the aggressive-retirement rule it is a live candidate, and the
`ddc1ee10` precedent from two commits ago on this very branch is the template.)

### S6. Two commits ago, this branch retired a `BoundaryTraceLaw` ClassVar for
### precisely the reason `geometry_map` / `response_kernel` would be retired

`ddc1ee10 refactor(geometry): retire creates_sweep_cycle — a per-law flag cannot
carry a configuration property`, whose message opens: *"The ClassVar had ZERO
production readers."* The same sentence is true of `geometry_map` and
`response_kernel`.

**But the discriminator differs, and it matters.** `creates_sweep_cycle` was
retired because it *could not have worked in principle* — "whether a face's
back-edge closes a cycle depends on the WHOLE FACE CONFIGURATION … which a
boolean on the boundary KIND cannot express." No such objection applies to
`geometry_map` / `response_kernel`: they CAN carry their fact (each law's `G`
and `R` are per-law facts, tabulated in prose in
`docs/theory/methods/sn/boundary_conditions.rst:187-218`), and there is a live
plan (P2b) building the analogous composition for reaction operators. So the
honest verdict here is the `coding-standards` L4 shape: **either populate them
this campaign or retire them** — the one thing that is not defensible is
leaving them declared, documented as first-class, and empty for another cycle.

### S7. Two theory-page claims are factually false against the tree

1. `docs/theory/foundations/boundary_conditions.rst:3167-3170` — "Every
   concrete BC … carrying only its parameters …, its :attr:`kind` tag, its
   :attr:`geometry_map` / :attr:`response_kernel` / :attr:`source` **property
   overrides**, and the relevant :meth:`assert_*` invariants." **No concrete
   law overrides `geometry_map` or `response_kernel`** — verified at runtime,
   all seven return `None`. Only `source` is genuinely overridden (by
   `PrescribedInflow`).
2. `orpheus/geometry/boundary/_base.py:89-91` — "Concrete subclasses
   (`VacuumInflow`, `ReflectiveBoundary`, `WhiteBoundary`, `AlbedoBoundary`,
   `PeriodicBoundary`, `PrescribedInflow`) **populate these**." They do not.

Minor, but they are the reason a reader would believe the tensor network exists.

### S8. `PeriodicWrapOperator` is a named type whose `apply` is `x.copy()`

`operator.py:2407-2410`. Its geometry lives in the sweep's face-pair indexing,
i.e. outside the operator algebra entirely. It is the one law where the
`geometry_map` factor is not merely un-extracted but genuinely **absent** — the
type is a reserved name with a "See follow-up issue" pointer
(`operator.py:2390-2394`). It is also SN-realizable-but-not-admitted
(`BOUNDARY_OPERATOR_REGISTRY` omits `"periodic"`), so nothing can reach it from
a `BC` tag today. Diffusion refuses it outright.

### S9. `BC.to_alpha()` is a second, parallel tag → physics map

`orpheus/geometry/mesh.py:71-120`. Maps `vacuum→0`, `reflective→1`,
`partial→params["albedo"]` for the trajectory-resolvent reference solvers,
bypassing the law/realizer chain entirely — and it recognizes a tag
(`"partial"`) that neither method's registry contains, while raising
`NotImplementedError` for `"white"`. Two tag vocabularies, one `BC` type.

### S10. SN's admitted-tag set is two, and the interesting laws are unreachable

`SNMesh.BOUNDARY_OPERATOR_REGISTRY` = `{"vacuum", "reflective"}`. `white`,
`albedo`, `periodic`, `prescribed_inflow` are fully realized by
`SNBoundaryRealizer` — with tests, invariants, and a documented tensor
decomposition — but no `BC` tag reaches them (`augmented_mesh.py:182-187`).
So the two laws that would most benefit from the factorization (white =
rank-one; periodic = missing `G`) are, today, only reachable by constructing a
law object and calling the realizer by hand. That materially lowers the urgency
of the refactor — and materially raises the risk that the ONE tag-reachable
non-trivial law (`reflective`, whose factorization is degenerate) is used to
judge whether the network is worth building.

---

## Appendix — drift notes

* **Brief vs tree.** `DiffusionBoundaryOperator` is at
  `orpheus/diffusion/operators.py:544`, not `:407` (`:407` is
  `LeakageOperator.apply`). `SNMethodSpace` is at
  `orpheus/sn/mesh/method_space.py:71` (class), `:72` is its docstring.
  `resolve_boundary_conditions` is at `orpheus/transport/method.py:218` —
  matches. `SNBoundaryRealizer.realize` at `:147` — matches.
  `DiffusionBoundaryRealizer.realize` at `:176` — matches.
  `DiffusionMethodSpace` at `orpheus/diffusion/method_space.py:56` (class
  statement), `:57` is the docstring.
* **Nexus vs tree.** The graph DB is at `8d996f53`; `ddc1ee10` /
  `1fd15f64` / `73627b71` have landed since. `ddc1ee10` deleted
  `BoundaryTraceLaw.creates_sweep_cycle` and `1fd15f64` made
  `BoundaryRealizer` generic (`Protocol[MethodSpaceT_contra]`) and tightened
  `law: Any → BoundaryTraceLaw`. Nexus was used only for orientation; every
  claim above was read from the working tree or executed against it.
* **The `_BoundBoundaryOperator` shim** (`orpheus/geometry/boundary/_bound_compat.py:77-178`)
  is a transparent 1-arg passthrough that forwards `is_invertible`,
  `is_adjointable`, `inverse()`, `block_role`, `apply`, `apply_transpose`, and
  implements `__eq__` against a raw string. It is what makes S3's
  `sn_mesh.bc[face] == "reflective"` comparisons work.
* **Realizer output verification** (executed against this tree, S4 GL,
  `face="xmin"`, `inflow_indices=[2,3]`):

  ```
  law                  α      -> realized type            factors                                        block_role  is_adjointable
  VacuumInflow         None      TensorProductOperator    [IncomingOrdinateMaskTensor, IdentityOperator]  BOUNDARY    True
  ReflectiveBoundary   1.0       TensorProductOperator    [PermutationOperator, IdentityOperator]         BOUNDARY    True
  ReflectiveBoundary   0.5       ScaledOperator           -                                               BOUNDARY    True
  WhiteBoundary        1.0       TensorProductOperator    [AngularAverageOperator, IdentityOperator]      BOUNDARY    False
  WhiteBoundary        0.5       ScaledOperator           -                                               BOUNDARY    False
  AlbedoBoundary       0.0       ZeroOperator             -                                               BOUNDARY    True
  AlbedoBoundary       1.0       IdentityOperator         -                                               BOUNDARY    True
  AlbedoBoundary       0.3       ScaledOperator           -                                               BOUNDARY    True
  PeriodicBoundary     None      TensorProductOperator    [PeriodicWrapOperator, IdentityOperator]        BOUNDARY    True
  ```

  Every `TensorProductOperator` here is 2-factor with the second factor an
  `IdentityOperator`. White is the sole non-adjointable law, which makes
  `SNBoundaryOperator.is_adjointable` False for any mesh carrying a white face.
* **Law-registry verification** (runtime): all seven registered laws
  (`albedo`, `periodic`, `prescribed_inflow`, `reflective`, `vacuum`, `white`,
  `zero_flux`) return `geometry_map is None` and `response_kernel is None`.

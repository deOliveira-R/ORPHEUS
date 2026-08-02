# Quadrant review — the SN COMPOSITE boundary machinery

**Scope:** `orpheus/sn/operators/boundary.py` (825 lines) and everything it
composes with. Branch `refactor/operator-strategy-layers` @ `73627b71`.
**Evidence key** as in `.claude/plans/boundary_machinery_review.md`:
**[M]** measured · **[R]** read+quoted · **[G]** exhaustive grep · **[U]** unverified.

Given (not re-derived): `_reflect_trace` writes only
`sel = trace.inflow_indices_for_face(face)` (`:302`), so `B(psi)`'s outflow rows
are zero for every law.

---

## 1. Class inventory — `orpheus/sn/operators/boundary.py`

**[R]** Four public names, `__all__` at `:85-90`. There are **no other classes**
in the file **[G]** (`grep -n "^class"` returns exactly the four below).

| class | line | carrier (domain = codomain) | role | `apply` signature |
|---|---|---|---|---|
| `SNBoundaryOperator` | `:123` | `FullField` (bulk ⊕ trace); `domain`/`codomain` = `sn_mesh.full_field_space` (`:215`, `:219`) | `B_a` — System A's trace `A_ss` block | `apply(self, psi: "FullField") -> "FullField"` (`:365`) |
| `RadialCharacteristicBoundaryOperator` | `:478` | `RadialCharacteristicField`; spaces = `sn_mesh.radial_characteristic_field_space` (`:567`, `:571`) | `B_b` — System B's ψ½ ray-corner `A_ss` block | `apply(self, ray: "RadialCharacteristicField") -> "RadialCharacteristicField"` (`:676`) |
| `SNMaskedBoundaryOperator` | `:715` | `FullField` — declared generically: `LinearOperator["FullField", "FullField"]` | one half of the schedule split (`B_lower` / `B_upper`) | `apply(self, psi: "FullField") -> "FullField"` (`:766`) |
| `BoundarySplit` | `:821` | — (`NamedTuple`) | the named `(lower, upper)` pair | n/a |

### 1.1 Monomorphism — confirmed

**[G]** `grep -n "singledispatch\|@overload" orpheus/sn/operators/boundary.py`
→ **zero hits**. Every `apply` / `apply_transpose` in the file is a plain
`def` with one positional carrier argument. Confirms the review's §1 claim
(the realized operators are monomorphic) for all three operator classes.

**[R]** The carrier separation is *structural, not a dispatch arm*: `B_a` and
`B_b` are different classes with **different `domain`/`codomain` spaces**, and
`B_b`'s constructor refuses to exist on a seedless mesh (`:539-548`):

```python
if sn_mesh.radial_characteristic_field_space is None:
    raise ValueError(
        "RadialCharacteristicBoundaryOperator: the mesh carries no "
        "ψ½ ray (radial_characteristic_field_space is None) — ...")
```

### 1.2 Class-attribute inventory (what each advertises)

**[R]**

| | `block_role` | `system_role` | `is_adjointable` | `apply_transpose` | `inverse` / `solve` |
|---|---|---|---|---|---|
| `SNBoundaryOperator` | `BOUNDARY` (`:175`) | *(unset → base default)* | per-face **intersection** (`:196-205`) | yes (`:467`) | none (base `is_invertible=False`) |
| `RadialCharacteristicBoundaryOperator` | `BOUNDARY` (`:534`) | `SystemRole.B` (`:537`) | **per-leaf kind test** (`:550-559`) | yes (`:680`) | none |
| `SNMaskedBoundaryOperator` | `BOUNDARY` (`:738`) | *(unset)* | *(not overridden → base)* | **absent by design** (`:732-735`) | none |

**[R]** `SNMaskedBoundaryOperator` docstring `:732-735`: *"A masked half is NOT
invertible and does not advertise a transpose (`B_lowerᵀ` masks input rows, not
output rows — mint it when the adjoint-inverse carve #280 produces a consumer),
so it is apply-only and the two-axis contract holds by the base defaults."*

---

## 2. `_reflect_corner` — the ruled/loud-deferred corner story

`RadialCharacteristicBoundaryOperator._reflect_corner` (`:573-628`).

### 2.1 What it is

**[R]** `:576-583` — *"The ``A_ss`` CORNER action on System B's boundary member
(R13, 2.5d). The (R, μ = ∓1) corner pair closes the ray boundary on a
seed-carrying mesh: the inward seed leg's r = R inflow is BC data, and for a
specular-reflective outer face the reflected partner of the outward ray μ = +1 is
EXACTLY the inward one μ = −1 (its own mirror — an off-quadrature ray, so the
per-face law OPERATOR cannot act on it; the specular fact is applied directly)."*

The body is 10 lines of real work (`:619-628`): allocate a zero
`RadialCharacteristicBoundarySourceSink`, and for `kind == "reflective"` copy
`corner(level,+1) → corner(level,−1)` per carried level (transpose: the mirror
image, `corner(level,−1) → corner(level,+1)`). For `kind == "vacuum"` the
all-zero `out` falls through — **no branch at all**, only a comment (`:620-621`).

### 2.2 When it fires

**[G]** Three call sites, all in this file: `_apply_faces` (`:673`) — i.e. both
`apply` (`:678`) and `apply_transpose` (`:689`) — and `reflect_corner_inplace`
(`:708`).
**[G]** `reflect_corner_inplace` has exactly ONE production caller:
`orpheus/sn/solver.py:1990`, inside `_reflect_boundary_inplace`. `B_b` itself is
constructed at `orpheus/sn/coupled_system.py:499` (the (B,B) grid slot) and at
`solver.py:1990`.
**[M]** The class is unconstructable on a seedless mesh (ctor guard `:540-547`) —
so it fires only on a sphere (1-D curvilinear seed-carrying mesh).

### 2.3 What it raises and why — `_RULED_CORNER_KINDS`

**[R]** `:100` — `_RULED_CORNER_KINDS = frozenset({"vacuum", "reflective"})`,
with the module-level rationale at `:92-99`:

> *"The outer-face law kinds with a RULED ψ½ corner action (RULING P1's ray
> carrier). Single source for BOTH `RadialCharacteristicBoundaryOperator.is_adjointable`
> (a ruled kind's corner map is Euclidean-adjointable) AND
> `RadialCharacteristicBoundaryOperator._reflect_corner` (an unruled kind is
> loud-deferred). white / albedo / periodic are absent — their off-quadrature
> μ = ±1 re-emission is a design decision not yet ruled (2.5d plan-of-record);
> add a kind here AND its branch in `_reflect_corner` when it is ruled."*

**[R]** The raise (`:612-618`) is a `NotImplementedError`:

```python
if kind not in _RULED_CORNER_KINDS:  # single source with is_adjointable
    raise NotImplementedError(
        f"RadialCharacteristicBoundaryOperator: the outer-face law kind "
        f"{kind!r} has no ruled corner action yet (white / albedo / "
        f"periodic at the off-quadrature μ = ±1 ray — loud-deferred, "
        f"2.5d plan-of-record).")
```

**[R]** The mathematical reason, quoted from the docstring's dispatch note
(`:592-602`):

> *"Law dispatch — on the declared law KIND (the #186 shim's `kind` tag, the
> same registry key `sn_mesh.bc[face] == "reflective"` comparisons read), NOT on
> the realized operator's composition tree (**which is an (ordinate ⊗ group)
> operator over the QUADRATURE rows and structurally cannot act on the
> off-quadrature μ = ±1 ray**)"* … *"anything else (white / albedo / periodic /
> prescribed) — **loud-deferred** (`NotImplementedError`) per the 2.5d
> plan-of-record (e.g. **white re-emission at the off-quadrature ray needs the
> `|Ω·n|`-weighted outflow average for μ = −1, not yet ruled**)."*

So: the ψ½ starting-direction ray sits at μ = ±1, which is **not a quadrature
ordinate**. The realized per-face law is a permutation / mask / average indexed
by ordinate row, so it has no row to act on for that ray. Reflective and vacuum
are ruled because their corner action is a *fact* independent of quadrature
(specular mirror of μ=+1 is exactly μ=−1; vacuum emits nothing). White's would
require re-deriving the hemispherical average at an off-quadrature direction.

**[M] Reachability of the deferral.** `SNMesh.BOUNDARY_OPERATOR_REGISTRY`
(`orpheus/sn/mesh/augmented_mesh.py:171-174`) admits **only** `{"vacuum",
"reflective"}` — identical to `_RULED_CORNER_KINDS`. **⟹ the
`NotImplementedError` at `:613` is UNREACHABLE from any `BC(...)` declaration
today.** (This also raises the parent review's §6 [U] to **[R]**: the registry
literally is `{"vacuum": VacuumInflow, "reflective": ReflectiveBoundary}`, with
`:181-187` stating the 5 other kinds *"are NOT registered here — adding them
requires SN-sweep-side wiring"*.)

**[G]** No test constructs the unruled arm either (`grep -rn "_reflect_corner"
tests/` → one hit, `tests/sn/operators/test_psi_half_coupling.py:165`, a
docstring noting the vacuum floor does not exercise the *reflective* arm).
**[G]** The reflective sphere IS production-reachable and exercised:
`tests/sn/solve/test_fixed_source_g1.py:74-77` builds `CoordSystem.SPHERICAL`
with `bc_right=BC("reflective")`.

---

## 3. `split` / `B_lower` / `B_upper` — the schedule split

### 3.1 What produces them

**[R]** `SNBoundaryOperator.split(schedule)` (`:431-465`) returns a
`BoundarySplit` NamedTuple of two `SNMaskedBoundaryOperator`s built from the row
partition:

```python
lower_rows = schedule.lower_inflow_rows(self.sn_mesh)
upper_rows = {face: np.setdiff1d(trace.inflow_indices_for_face(face),
                                 lower_rows.get(face, np.empty(0, dtype=np.intp)))
              for face in self._face_laws}
```

**[R]** The row law lives in `SweepSchedule.lower_inflow_rows`
(`orpheus/sn/loss_representation/sweep_schedule.py:169-214`):

> *"face `f` is reflected exactly once, after its LAST outflowing octant group
> (`reflect_faces`), at which point every outflow feeding `f` is complete. An
> inflow ordinate row `(f, m')` therefore reads the FRESH current-iterate
> reflection iff its octant group is swept strictly AFTER `f`'s reflect group —
> those rows are `B_lower` (realized in-sweep by the forward substitution); rows
> swept at-or-before the reflect keep the lagged seed — `B_upper`."*

### 3.2 Who consumes them

**[G]** The only `B.split(...)` production call: `orpheus/sn/solver.py:772`
(`parts = B.split(SweepSchedule.gauss_seidel(sn_mesh))`), guarded at `:765-768`
by an `isinstance(B, SNBoundaryOperator)` check (multi-D Cartesian ⟹ seedless ⟹
plain `B`). The halves then feed `ScheduledInvertibleOperator`
(`orpheus/sn/operators/scheduled_invertible.py`) as `M = (L+C) − B_lower` and
the gain pair `(S, B_upper)`.
**[R]** `scheduled_invertible.py:115-118` re-checks the type:
`if not isinstance(lower, SNMaskedBoundaryOperator): raise ...` — a constructor
guard naming `SNBoundaryOperator.split` as the only legal producer.
**[R]** `StreamingOperator.__sub__` special-cases the masked half
(`orpheus/sn/operators/streaming.py:709`, `:735-737`) so `(L+C) − B_lower` is
*"the canonical spelling — the dispatch lives"* there.

### 3.3 The exactness invariant — quoted

**[R]** `boundary.py:441-446`:

> *"The partition is exact: **the specular map has no octant-diagonal**, and the
> two row sets are complementary within each face's inflow by construction here.
> Returns a named pair so the two construction sites cannot be swapped silently:
> `M = (L + C) - parts.lower` and `gains = (S, parts.upper)`."*

**[R]** The same invariant restated in the schedule (`sweep_schedule.py:184-186`):
*"The specular map flips one direction-cosine sign, so a row and its source
always sit in DIFFERENT octants — `B` has no octant-diagonal and
`B = B_lower + B_upper` is an exact partition."*

**[R]** The row-granularity claim inside `_reflect_trace` (`:275-281`):

> *"`rows` (#226 step 2) restricts WITHIN a face: per face, only the given
> ordinate rows of the codomain projection are emitted (a subset of the inflow
> rows — the schedule-split `B_lower`/`B_upper` halves of `split`). A face absent
> from `rows` emits nothing. **Row-granular restriction is exact for the same
> reason the face restriction is: the projected action writes each target row
> independently.**"*

That is the load-bearing invariant: because the composite's forward action is
`out[sel] = law.apply(face_in)[sel]` — a *scatter into a zero buffer indexed by
target row* — restricting `sel` restricts the operator exactly, with no
cross-row coupling dropped. It holds for ANY law, including white (whose input
support is the whole face) — the restriction is on the CODOMAIN, and the
codomain rows are independent.

**[M]** Verified numerically: `_reflect_trace` is a scatter whose image support
is exactly `sel` — measured `nnz(B|trace) = 8` on `slab_seedless` (32 trace dof)
and `216` on `cart2d_seedless` (864 trace dof), all confined to inflow rows of
reflective faces (see §4.2).

---

## 4. `is_adjointable` and the adjoint path

### 4.1 The two different rules

**[R]** `SNBoundaryOperator.is_adjointable` (`:196-205`) — the per-face
**intersection**:

```python
laws = self._face_laws.values()
return bool(laws) and all(law.is_adjointable for law in laws)
```

**[R]** `RadialCharacteristicBoundaryOperator.is_adjointable` (`:550-559`) — a
per-leaf **kind membership test**, NOT an intersection (`B_b` has one face):
`return kind in _RULED_CORNER_KINDS`, with `kind = getattr(self.sn_mesh.bc["xmax"], "kind", None)`.

### 4.2 Which laws are adjointable — MEASURED

**[M]** Every law the SN realizer produces, realized standalone on a GL-S8
`xmax` method space, reporting `op.is_adjointable`:

| law | realized type | `is_adjointable` | `block_role` |
|---|---|---|---|
| `VacuumInflow` | `TensorProductOperator` | **True** | BOUNDARY |
| `ReflectiveBoundary(α=1)` | `TensorProductOperator` | **True** | BOUNDARY |
| `ReflectiveBoundary(α=0.5)` | `ScaledOperator` | **True** | BOUNDARY |
| `WhiteBoundary` | `TensorProductOperator` | **False** | BOUNDARY |
| `AlbedoBoundary(0.0)` | `ZeroOperator` | **True** | BOUNDARY |
| `AlbedoBoundary(1.0)` | `IdentityOperator` | **True** | BOUNDARY |
| `AlbedoBoundary(0.4)` | `ScaledOperator` | **True** | BOUNDARY |
| `PeriodicBoundary` | `TensorProductOperator` | **True** | BOUNDARY |
| `PrescribedInflow(NoSource)` | `IncomingSourceOperator` | **False** | **None** |
| `ZeroFluxBoundary` | — | raises `BoundaryError` | — |

**[M]** On the three production fixtures, the composite is adjointable in every
case (`slab_seedless` reflective\|vacuum → True; `cart2d_seedless` 4 faces →
True; `sphere_carrying` (xmax vacuum) → True).

> **⟹ `SNBoundaryOperator.is_adjointable` can only return `False` for a law the
> SN `BOUNDARY_OPERATOR_REGISTRY` does not admit** (white / prescribed). The
> intersection rule is a real rule computing a value that, on every declarable
> mesh, is `True`. See §8.

### 4.3 What happens on `.H` — the promise that does NOT hold

**[R]** `LinearOperator.adjoint()` (`orpheus/numerics/operator.py:847-854`)
gates `.H` **eagerly on the very same predicate**:

```python
if not adjointable(self):
    raise MissingAdjoint(f"{type(self).__name__} is not adjointable — .H/.adjoint() ...")
return _AdjointOperator(self)
```

**[R]** `SNBoundaryOperator`'s class docstring claims the opposite (`:151-154`,
`:161-164`):

> *"``B`` acts on the composite as the ``A_ss`` block (zero bulk; non-zero only
> on the trace block, where the cosine-weighted ``|Ω·n|·w`` partial-current
> metric lives). That block metric is what makes the Hilbert adjoint ``B.H`` the
> physically correct partial-current adjoint — **the one channel by which the
> white-BC adjoint becomes available**."*
> … *"**white does NOT** [advertise `apply_transpose`] (it is self-adjoint only
> under the ``|Ω·n|·w`` metric, so **its adjoint routes through ``B.H``** on the
> weighted trace space at O.2)."*

**[M] That channel is closed.** With a realized `WhiteBoundary` installed at
`slab_seedless`'s `xmin` (built through `SNBoundaryRealizer` exactly as the
realizer would, wrapped in `_BoundBoundaryOperator(kind="white")`):

```
per-face: {'xmin': ('white', False), 'xmax': ('vacuum', True)}
B.is_adjointable = False
B.H  RAISED MissingAdjoint: SNBoundaryOperator is not adjointable — .H/.adjoint() requires
             is_adjointable=True (a working apply_transpose on every constituent).
_reflect_trace(..., "apply_transpose") RAISED MissingAdjoint:
             face 'xmin' law _BoundBoundaryOperator has no Euclidean transpose
forward with white: |B psi|_inf = 1.0     (apply works fine)
```

`.H` is gated on `is_adjointable`; `is_adjointable` is False exactly when white
is present; therefore **`B.H` is unreachable in precisely the case the docstring
says it exists to serve.** `_AdjointOperator.apply` would in any case route back
through `inner.apply_transpose` → the same per-face `MissingAdjoint` raise
(`:310-316`). This is a **promise the code cannot deliver** (§8, item 1).

### 4.4 Is the metric channel doing anything today?

**[M]** The trace metric IS populated — `full_field_space.apply_metric` on a
unit trace returns non-trivial per-ordinate values (`slab_seedless`: 4 distinct
values `{0.0665, 0.0972, 0.1649, 0.1772}` = `|μ_n|·w_n`; `cart2d_seedless`: 2
distinct). The bulk block of `B.H`'s output is exactly zero, as `BlockRole.BOUNDARY`
requires.

**[M]** But for every registry-reachable law the metric is **inert**:
`‖B.H(x) − Bᵀ(x)‖_∞ = 2.8e-17` (slab, `‖Bᵀx‖_∞ = 2.33`) and `4.4e-16`
(cart2d, `‖Bᵀx‖_∞ = 3.11`) on random trace data. Structurally: `B = P_in ∘ Π`
with `Π` the specular permutation; `G = diag(|μ|w)` is invariant under the
mirror, so `G⁻¹ Π G = Π` and the metric cancels; the vacuum face is the zero
map. So `.H` and `ᵀ` coincide for everything declarable.

### 4.5 The `(P_sel ∘ law)ᵀ = lawᵀ ∘ P_sel` identity — is it correct?

**[R]** The claim, `boundary.py:239-242` and `:304-306`:

> *"the forward action is projected onto the inflow rows: `B_face = P_inflow ∘
> law`, and the Euclidean transpose is `B_faceᵀ = lawᵀ ∘ P_inflow` — mask the
> INPUT to the forward's codomain rows, write the full `lawᵀ` image."*

**YES — the identity is correct, and the implementation realizes it exactly.**
Three things have to hold, and all three do:

1. **The algebra.** `(P·L)ᵀ = Lᵀ·Pᵀ`. `P_sel` is the 0/1 *diagonal* selector on
   the face's ordinate axis (`out[sel] = ·[sel]` written into a zero buffer), so
   `P_selᵀ = P_sel` and `(P_sel·L)ᵀ = Lᵀ·P_sel`. **[R]** Domain and codomain of
   a face block are the SAME index set (that face's `(N, …)` slot), so `P_sel`
   is an endomorphic orthogonal projector and the transpose is well-posed.
2. **The code.** `:317-319` — `masked = np.zeros_like(face_in); masked[sel] =
   face_in[sel]` is literally `P_sel·face_in`; then
   `out_boundary.face_view(face)[...] = law.apply_transpose(masked)` writes the
   FULL image of `lawᵀ`. That is `lawᵀ ∘ P_sel`, with no output projection —
   correct, and deliberately asymmetric with the forward (which projects the
   OUTPUT and not the input).
3. **`law.apply_transpose` really is `lawᵀ`.** **[R]**
   `PermutationOperator.apply_transpose` gathers with `inverse_perm =
   np.argsort(perm)` (`operator.py:2254-2255`) — the genuine transpose, not a
   re-application of `perm` (which would be right only for an involution).
   `IncomingOrdinateMaskTensor.apply_transpose` returns `self.apply` and is
   honestly self-adjoint (`operator.py:2367-2373`).

**[M] Verified numerically, not just argued.** Materializing the trace-block
matrices column-by-column (`_reflect_trace` applied to every unit trace vector):

| fixture | trace dof | `‖F‖∞` | `‖T‖∞` | **`‖T − Fᵀ‖∞`** | nnz(F) |
|---|---|---|---|---|---|
| `slab_seedless` | 32 | 1 | 1 | **0** (exact) | 8 |
| `cart2d_seedless` | 864 | 1 | 1 | **0** (exact) | 216 |

`apply_transpose` is the **bit-exact** Euclidean transpose of `apply`.

**Caveat worth naming (not a bug today).** The identity is correct for a
*fixed* `sel`. `_reflect_trace` computes `sel` from `rows` when `rows is not
None` — but the transpose branch is only ever reached with `rows=None`
(**[G]** `SNMaskedBoundaryOperator` calls `_apply_faces(..., "apply", rows=...)`
only; `apply_transpose` never passes `rows`). If a future consumer transposes a
masked half, `sel` would correctly be the mask rows — the code would still be
right; but `SNMaskedBoundaryOperator` deliberately does not advertise a
transpose (`:732-735`), so the composition is currently unspellable.

**The ⚠ warning at `:244-251` is also accurate.** Output-projecting `lawᵀ` onto
the outflow rows would give, for vacuum (`law = P_out`, self-adjoint),
`P_out·P_out = P_out ≠ 0` while the forward `P_in·P_out = 0` — a spurious `+1`
outflow diagonal. For a pure permutation law both spellings agree bitwise (the
map is off-diagonal), which is exactly why *"every reflective-fixture gate stayed
green over the wrong one"*.

---

## 5. `SNMaskedBoundaryOperator` — why it exists, who builds it

### 5.1 Why

**[R]** It is the row-restricted half of the boundary Gauss-Seidel regular
splitting. Class docstring `:716-730`:

> *"One half of the schedule split `B = B_lower + B_upper` — the whole-trace
> `SNBoundaryOperator` restricted to a per-face set of inflow ordinate ROWS
> (#226 §17 W2). … Which rows belong to which half is SCHEDULE-order semantics
> (`SweepSchedule.lower_inflow_rows`), so the instance carries its `schedule` —
> the reified `M = (L+C−B_lower)` reads the walk order off its lower operand
> rather than pairing a foreign schedule with a mismatched mask."*

### 5.2 Who constructs it

**[G]** Exactly one construction site: `SNBoundaryOperator.split` (`:462-465`).
**[G]** One production caller of `split`: `orpheus/sn/solver.py:772`
(`_within_group_splitting`, the `inner_schedule == "gauss_seidel"` and
`is_cartesian and not is_1d` arm) → `return LC - parts.lower, (S, parts.upper)`.
**[R]** `orpheus/sn/operators/streaming.py:735-741`: `StreamingCollisionOperator.__sub__`
`isinstance`-dispatches on `SNMaskedBoundaryOperator` to build the
`ScheduledInvertibleOperator`. **[R]** `scheduled_invertible.py:115-120`
re-guards the type in its constructor.

### 5.3 How it differs from `SNBoundaryOperator`

| | `SNBoundaryOperator` | `SNMaskedBoundaryOperator` |
|---|---|---|
| owns the mesh | yes (`self.sn_mesh`) | delegates (`self.inner.sn_mesh`, `:754-756`) |
| owns `domain`/`codomain` | yes | delegates (`:758-764`) |
| forward | `_apply_faces(psi, "apply")` — all faces, all inflow rows | `_apply_faces(psi, "apply", rows=self.rows)` — the same body, restricted |
| `apply_transpose` | present (`:467`) | **absent by design** |
| `is_adjointable` | intersection rule → **True** on every production mesh | **not overridden ⟹ base `False`** |
| extra verb | `reflect_inflow_inplace` (whole-face ASSIGNMENT) | `reflect_rows_inplace` (row-restricted **ADDITIVE** `+=`) |

**[R]** The additive-vs-assign distinction is the load-bearing one (`:773-795`):

> *"solving `M z = y` on a strictly-lower inflow row reads
> `z_in = y_row + (B z)_row` — the buffer already holds the seed `y_row`, so
> ACCUMULATING the fresh reflection completes the inhomogeneous row exactly. …
> ⚠ Additive, NOT whole-face overwrite: the dissolved resolvent's OVERWRITE
> dropped `y_row` — benign in production (zero on a reflective face) but
> O(1)-wrong as an inverse."*

**[R]** `reflect_rows_inplace` is passed as the `reflect=` callback into the
scheduled walk: `scheduled_invertible.py:260` (`reflect=self.lower.reflect_rows_inplace`).

### 5.4 MEASURED — the split is exact, and it costs the adjoint

**[M]** `cart2d_seedless` + `SweepSchedule.gauss_seidel`:

- `lower.rows` = `{xmin: [4,5,6,7,12,…,23] (12), ymin: [6,7,14,15,22,23] (6)}` → 18 rows
- `upper.rows` = `{xmin: [] , xmax: (12), ymin: (6), ymax: (12)}` → 30 rows
- per face: `union(lower, upper) == inflow_indices` **True** and
  `intersect(lower, upper) == ∅` **True**, for all four faces.
- `‖B − (B_lower + B_upper)‖∞ = 0` on the materialized 864×864 trace matrix;
  `nnz(B)=216 = nnz(lower)=168 + nnz(upper)=48`. **The partition is exact,
  measured, not just argued.**
- **[M]** 24 of `B_upper`'s 30 rows sit on the two VACUUM faces and contribute
  exactly zero (nnz 48 comes only from `ymin`'s 6 reflective rows).

**[M] The adjoint cost:** because `SNMaskedBoundaryOperator` does not override
`is_adjointable`, the whole G-S composite loses its adjoint —

```
M = (L+C) - B_lower  -> ScheduledInvertibleOperator  is_adjointable=False  is_invertible=True
(M - B_upper).is_adjointable = False        <- the G-S within-group A
((L+C) - B).is_adjointable   = True         <- the Jacobi within-group A
```

So the boundary-G-S path is structurally **adjoint-less**, while the Jacobi path
is adjointable. Worth flagging against the #276 adjoint work: an adjoint SN
solve on a 2-D Cartesian mesh must take the Jacobi arm, or `B_lowerᵀ`/`B_upperᵀ`
must be minted (the docstring's own "mint it when #280 produces a consumer").

---

## 6. `RadialCharacteristicBoundaryOperator` — genuinely a different carrier

**Confirmed: a separate class for a separate carrier, NOT a dispatch arm.**
**[R]/[M]** Four independent pieces of evidence:

1. **Different spaces.** `domain = codomain = sn_mesh.radial_characteristic_field_space`
   (`:562-571`), vs `B_a`'s `full_field_space`. They are not even comparable
   carriers — the `CoupledOperator` grid type-checks the placement
   (`:563-566` comment).
2. **Different `system_role`.** `system_role = SystemRole.B` (`:537`); `B_a`
   carries none at class level and is stamped `SystemRole.A` on the composed
   `A_AA` at `coupled_system.py:469`.
3. **Structurally unconstructable where the other lives** — the ctor raises on
   a seedless mesh (`:539-548`). **[M]** `slab_seedless` / `cart2d_seedless`
   have `radial_characteristic_field_space is None`.
4. **Different grid slots.** `orpheus/sn/coupled_system.py`: `B_a` at `:465`
   feeds `A_AA = LC - S - B_a` (the (A,A) slot); `B_b` at `:499` feeds
   `A_BB = march - B_b` (the (B,B) slot) — `march` being
   `RadialCharacteristicOperator`. **[R]** The loss grid is
   `[[A_AA, A_AB], [-emission, A_BB]]` (`:515-517`).

**What it does.** It closes the ψ½ ray's `r = R` boundary by the `(R, μ = ∓1)`
corner pair: forward copies `corner(level, +1) → corner(level, −1)` for a
reflective outer face (the specular partner of the outward μ=+1 ray IS the
inward μ=−1 ray); the transpose copies the mirror image; vacuum emits nothing.
It emits a `RadialCharacteristicField` with a **zero-source interior**
(`:671-674`) — System B's own composite, no bulk/trace padding (`:634-640`).
Its `reflect_corner_inplace` (`:691-712`) is the in-place façade called once per
solve from `orpheus/sn/solver.py:1990`, as the ray sibling of `B_a`'s
`reflect_inflow_inplace` at `:1981` — *"one reflect per system (RULING P1)"*.

**Adjoint metric.** **[R]** `:513-523`: `B_b` advertises the **Euclidean**
transpose only and **no per-leaf `.H`**, because *"The ray corner gauge is
symmetric (`g₊ = g₋ = V(R)` — both corners at r = R), so
`B_b.H = G_sd⁻¹ B_bᵀ G_sd = B_bᵀ`."* That is a genuine, stated reason — unlike
`B_a`'s metric claim in §4.3.

---

## 7. The `G_in` triplication — RESOLVED, and the prior claim needs correcting

The prior review asserts the inflow injection is *"implemented three times, three
ways."* **[M] Half right.** All three are built from the SAME index set
(`trace.inflow_indices_for_face(face)` / `method_space.inflow_indices`), but they
are **not three spellings of one map** — they are **two spellings of `P_in` plus
one spelling of its complement `I − P_in`.**

Materialized as 24×24 matrices on `cart2d_seedless`'s `xmin` face
(inflow = `[4,5,6,7,12,13,14,15,20,21,22,23]`):

| # | site | code | matrix | measured |
|---|---|---|---|---|
| (a) | `boundary.py:302` | `out[sel] = full[sel]` into a `zeros_on` buffer | diag = **1 on inflow**, 0 elsewhere | — |
| (b) | `operator.py:2355-2365` | `out = x.copy(); out[inflow_indices] = 0.0` | diag = **0 on inflow**, 1 elsewhere | `‖b − (I − a)‖∞ = 0` |
| (c) | `angular.py:345-347` | `q * self._inflow_mask.reshape(...)` where `mask[indices]=1.0` | diag = **1 on inflow**, 0 elsewhere | `‖a − c‖∞ = 0` |

**[M]** `‖a − c‖∞ = 0` — (a) and (c) are the **same linear map**, bit-identical,
realized once as an index scatter-write and once as a dense diagonal multiply.
**[M]** `‖b − (I − a)‖∞ = 0` and `a @ b = 0` — (b) is the **orthogonal
complement**, and the two compose to the zero map (which is exactly why
`B(ψ)|vacuum ≡ 0`).
**[M]** All three are diagonal (`offdiag nnz = 0` for each) — pure projectors,
no ordinate mixing.

**[R]** The three are named as if unrelated: (b) is a typed operator class
called `IncomingOrdinateMaskTensor` whose docstring says it is *"projection onto
the outflow subspace"* (`operator.py:2306`) — i.e. **the class named "Incoming…"
IS the outgoing projector**; (c) is a private attribute `_inflow_mask` on
`IncomingSourceOperator`; (a) is an anonymous slice-write with no name at all.

**Architectural reading (the load-bearing part):** the concept — *the face's
inflow indicator, and the projector `P_in` it induces* — is a **single object
materialized three times in three representations, one of which is the
complement, none of which is typed as `P_in`.** The only typed member of the
family is the complement `P_out`, named "Incoming". A `P_in`/`P_out` pair
(one type, a `.complement`) would make `P_in ∘ P_out = 0` a *provable* property
instead of an emergent measurement, and would delete two of the three spellings.

---

## 8. Promises this file does not deliver — fresh-eyes inventory

Eight findings, all raised to [M]/[R]/[G]. Ordered by severity.

### 8.1 `B.H` as the white-BC adjoint channel — **the promise is structurally impossible**

Covered in §4.3. **[R]** `:152-154` and `:161-164` claim white's adjoint *"routes
through `B.H`"*. **[M]** `.H` is gated on `is_adjointable`, which is `False`
exactly when a white face is present ⟹ `MissingAdjoint` at construction. **[M]**
And for every law that CAN reach `.H`, the metric is inert (`‖B.H − Bᵀ‖∞ ≲ 4e-16`).
The stated architectural justification for `B`'s `FullFieldSpace` domain (*"That
block metric is what makes the Hilbert adjoint `B.H` the physically correct
partial-current adjoint"*) currently buys nothing measurable.

### 8.2 `_zero_bulk_source` — "single source for BOTH `B_a` and `B_b`" is false

**[R]** `:112-114`: *"Single source (Cardinal Rule 2) for the zero bulk both
`B_a` (`SNBoundaryOperator`) and `B_b` (`RadialCharacteristicBoundaryOperator`)
emit."*
**[G]** `grep -rn "_zero_bulk_source" orpheus/ tests/` → **two hits, both in this
file**: the `def` at `:103` and one call at `:361` (inside
`SNBoundaryOperator._apply_faces`). **[R]** `B_b._apply_faces` builds
`RadialCharacteristicInteriorSourceSink.zeros_on(mesh)` (`:672`) instead — a
different carrier entirely, since the B.2b re-type. So the Cardinal-Rule-2
"single source for both" is a one-consumer helper whose second consumer left.

### 8.3 The class docstring says `B_a` emits a ray block; `_apply_faces` says it does not

**[R]** A direct self-contradiction inside one file:
- `:140-141`: *"`B_a` emits a **present-zero** ray block so `B_a + B_b` sums
  bit-identically; the ray corner is entirely `B_b`'s."*
- `:336-337`: *"B.2d: System B is its own composite, so `B_a` **neither reads nor
  pads a ray block**."*

**[M]** The code agrees with `:336-337`: on `sphere_carrying`, `B_a.apply(ψ)`
returns a `FullField` whose dataclass fields are exactly `['interior',
'boundary']` — no ray member exists to be present-zero.
**[M]** And the sum is not merely unnecessary, it is **unspellable**:
```
B_a + B_b  ->  IncompatibleOperatorComposition: OperatorSum requires equal domains;
               got FullFieldSpace('full_field', shape=(112,))
               and FullFieldSpace('radial_characteristic', shape=(28,)).
```
`B = B_a + B_b` is a *conceptual* direct sum over the coupled grid
(`coupled_system.py:99` — `N = [[S + B_a, ∅], [+Emission, B_b]]`), never an
`OperatorSum`. The module docstring (`:25-30`) and `:136`, `:140-141` all read as
if the `+` were real.

### 8.4 `B_b.domain`'s comment cites a RETIRED object

**[R]** `:565-566`: *"the FullField-summed production `B = B_a + B_b` rides the
transient adapter, which declares `full_field_space`."*
**[G]/[R]** The transient adapter is gone — `tests/sn/operators/test_psi_half_coupling.py:2126-2128`
says so explicitly: *"G-b3.3 (the B.2b adapter byte-identity + delegation gate)
RETIRED at B.2d d1 with its referent: **the transient FullField-gain adapters are
gone** (the driver consumes the blocks natively through the gain grid)."*
**[G]** No adapter class exists under `orpheus/sn/`.

### 8.5 `reflect_inflow_inplace` names a consumer that passes a DIFFERENT method

**[R]** `:408-411`: *"the `_sweep_scheduled` inter-group reflect passes THIS bound
method (#226 step 2 — the reified `M = (L+C−B_lower)` supplies
`boundary.reflect_inflow_inplace`)"*.
**[R]** The reified `M` passes the OTHER verb:
`scheduled_invertible.py:260` — `reflect=self.lower.reflect_rows_inplace`
(`SNMaskedBoundaryOperator`'s ADDITIVE row update).
**[G]** `reflect_inflow_inplace` has exactly ONE caller repo-wide:
`orpheus/sn/solver.py:1981`, inside `_reflect_outflow_into_inflow` — which IS the
docstring's *second* named consumer, so half the sentence is right.
This matters beyond hygiene: `reflect_rows_inplace`'s own docstring (`:788-790`)
warns that whole-face OVERWRITE *"dropped `y_row` … O(1)-wrong as an inverse"* —
i.e. the two verbs are **not** interchangeable, and `:408-411` names the wrong one
at the exact site where the distinction is load-bearing.

### 8.6 A stringly-typed `method` tag, unvalidated, discriminated in two classes

**[R]** Both `_reflect_trace` (`:300-319`) and `_reflect_corner` (`:622-627`)
dispatch `if method == "apply": … else: <transpose>` on a bare `method: str`.
There is no `elif method == "apply_transpose"` and no validation.
**[M]** Consequence — an unknown tag silently performs the TRANSPOSE:

```
B._reflect_trace(bf, "banana")  == apply_transpose result :  True
B._reflect_trace(bf, "")        == apply_transpose result :  True
B_b._reflect_corner(seed, "banana") -> corner(+1)=3.0  (the transpose image)
```

The callers are all internal literals today (**[G]** 9 call sites, all in this
file, all `"apply"` / `"apply_transpose"`), so it does not bite now — but this is
the `discriminations` smell (one tag, two classes, no type) sitting inside the
subsystem the parent review is holding up as the typed-architecture template.

### 8.7 The `is_adjointable` intersection rule is gated on a condition never met

**[M]** §4.2: every law the SN `BOUNDARY_OPERATOR_REGISTRY` admits (`vacuum`,
`reflective`) is adjointable. The only non-adjointable realizations are
`WhiteBoundary` and `PrescribedInflow`, **neither of which is registered**
(`augmented_mesh.py:171-174` + the `:181-187` note). So
`SNBoundaryOperator.is_adjointable` returns `True` for every mesh that can be
declared, and the per-face `MissingAdjoint` raise at `:310-316` — described in
its own comment as *"unreachable-in-practice"* — is unreachable in fact too.
Same for `_reflect_corner`'s `NotImplementedError` (`:613`, §2.3).
**Not necessarily wrong** — these are honest guards for a registry that is
*expected* to grow. But it means three of this file's capability-gates are
currently theatre, and the review's "declares capability ahead of consumption"
pattern (parent §4) extends to them.
Minor sibling: `is_adjointable`'s `bool(laws)` guard (`:205`) protects against an
empty face set; **[M]** all three fixtures have 1/2/4 faces and the SN mesh
builds faces from `face_labels`, so the empty case looks unreachable too.

### 8.8 Smaller items

- **[M] Dead compute in `B_upper.apply`.** `_reflect_trace` filters faces by
  `f in rows` (`:293`) but does NOT skip an EMPTY row array; it still calls
  `law.apply(face_in)` and then scatters through an empty index
  (`:301-302`). Measured: `parts.upper.rows['xmin']` has size 0 yet `xmin`
  survives the filter, so every `B_upper.apply` computes and discards the xmin
  specular permutation. `reflect_rows_inplace` DOES skip empties (`:797-802`) —
  the two paths disagree.
- **[R] The class docstring describes the UN-projected action.** `:129-132`:
  *"`B_a.apply(ψ)` returns a `FullField` with zero bulk and, on each face,
  `bc[<face>].apply(ψ.boundary.face_view(<face>))`"* — omitting the inflow
  projection that `_reflect_trace` exists to apply. Since that projection is
  precisely what zeroes the outflow rows (the parent review's §4.1 finding), the
  class-level description is materially incomplete at the exact point of
  confusion.
- **[R]/[M] A wrong reason attached to a right conclusion.** `:202-203`:
  *"is_invertible inherits base False (**a BC reflection map is not
  invertible**)."* **[M]** The reflective per-face law reports
  `is_invertible=True` (it is a `PermutationOperator`). `B` is non-invertible
  because `P_in ∘ law` is rank-deficient, not because reflection is.
- **[R] `apply_transpose` failures report `apply`.** `_apply_faces`'s two guards
  hardcode the method name: `"SNBoundaryOperator.apply: input field and operator
  must share the same SNMesh instance…"` (`:345`) fires for `apply_transpose`
  too. `B_b`'s sibling does it correctly (`context=f"…{method}"`, `:658`).
- **[G] `BoundarySplit` is exported but never imported.** It is in `__all__`
  (`:86`) yet no module imports the name — consumers use `parts.lower` /
  `parts.upper` off the return value.
- **[R] Friend-class coupling.** `SNMaskedBoundaryOperator` reaches into
  `self.inner._apply_faces` (`:768`) and `self.inner._reflect_trace` (`:804`) —
  two private methods of a different class. Deliberate (it IS a restriction of
  `B`), but it means the "masked half" is not expressible through
  `SNBoundaryOperator`'s public surface.

---

## 9. What is genuinely sound here (so the plan does not over-correct)

- **[M]** `apply_transpose` is the **bit-exact** Euclidean transpose of `apply`
  (`‖T − Fᵀ‖∞ = 0` on both fixtures), and its `(P_sel∘law)ᵀ = lawᵀ∘P_sel`
  reasoning is mathematically correct, correctly implemented, and gated by
  tests (`tests/sn/operators/test_psi_half_coupling.py:644-660` pins
  `dense(B_aᵀ) ≡ dense(B_a)ᵀ`; `test_g_adjoint_reciprocity.py` and
  `test_capability_survival.py` exercise the predicate).
- **[M]** The `B = B_lower + B_upper` partition is **exactly** complementary and
  exhaustive, measured on the materialized 864×864 trace matrix.
- **[M]** The three operator classes are genuinely monomorphic; `B_b` is a real
  separate carrier, not a dispatch arm (four independent confirmations, §6).
- **[R]** `B_b`'s "Euclidean = Hilbert because the ray corner gauge is symmetric"
  (`:513-523`) is a *stated, checkable* reason — the model the `B_a` metric claim
  in §8.1 should be held to.

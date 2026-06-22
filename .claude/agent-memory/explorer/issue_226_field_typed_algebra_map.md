---
name: issue-226-field-typed-algebra-map
description: Architecture-campaign map for the FIELD-TYPED operator algebra — operator census by true (Vin,Vout), the BoundaryMomentField gap, P3.6 DirectSumSpace/TimedFullField→Field promotion status, the scipy flat-ndarray serialization boundary, and the retype/consumer blast radius. Distinct from issue-226-operator-generics-map (the narrower Generic[V] pyright change).
metadata:
  type: project
---

# Field-typed operator-algebra campaign map (read-only)

Point-in-time `file:line` are current at **2026-06-19** on `main`. The SHAPE
(buckets, gap analysis, serialization mechanics, blast-radius dimensions) is the
durable part; line numbers drift — re-derive via Nexus `context`/`query`.

> **Scope note / naming hazard.** The task framed this as "#226". GH **#226 is
> actually** *"pyright: ~691 real type findings…"* (OPEN, `module:tests`), and
> the sibling memo [[issue-226-operator-generics-map]] covers the narrow
> `Generic[V]` typing change. **THIS memo is the broader ARCHITECTURE campaign**:
> push the operator algebra to speak typed FIELDS (not ndarrays), close the
> `BoundaryMomentField` gap, land the P3.6 `DirectSumSpace`/`TimedFullField→Field`
> promotion, and demote flat-ndarray to an explicit serialization adapter. If a
> GH issue is filed for this campaign it is NOT #226 — file a fresh one.

## 0. Headline verdicts (for the campaign plan)

1. **`DirectSumSpace` DOES NOT EXIST.** It is named ONLY in docstrings/plans
   (`timed_full_field.py:46-51,98,149`, `space.py:42-59` "future direction",
   Grand Report v3 §5.3) and one test docstring (`test_spherical_harmonic_space.py:416`).
   There is NO class, NO file. P3.6 = *write it + promote `TimedFullField` onto it*.
   (⚠ The "Phase 3.6" in `derivations/continuous/mms/sn.py` + `discrete_ordinates.rst`
   is a DIFFERENT phase numbering — the anisotropic-curvilinear-MMS phase — UNRELATED
   to the TimedFullField promotion. Do not conflate.)

2. **The flat ndarray is a SERIALIZATION boundary, not an algebra vector space.**
   `as_scipy_linop` (`operator.py:1661`) has **ZERO production callers** (only 5 sites
   in `tests/numerics/test_operator.py`). The SN Krylov path builds its OWN
   `spla.LinearOperator` inline (`iteration.py:755`) whose `A_matvec`
   (`iteration.py:743-753`) ravels field→flat→`apply`→flat via the `to_flat`/`from_flat`
   ravellable protocol. So operators NEVER actually consume the flat vector in
   production — scipy does, at one closure. A "typed serialization adapter" formalizes
   THIS one site, not `as_scipy_linop`.

3. **The genuinely-flat primitives are mostly UNUSED-in-production or
   axis=0-on-the-ordinate-trace.** Of the §6.c primitives, only
   `PermutationOperator`/`IncomingOrdinateMaskTensor`/`PeriodicWrapOperator`
   (boundary realizer) and `RankOneOperator`/`TensorProductOperator`/`IdentityOperator`
   (fission kernel) have production instantiation. ALL operate on a single tagged numpy
   axis (the ordinate axis `axis=0`, or the group axis) of a field's `.values` ndarray —
   they are **field-AXIS operators wearing a flat costume**, expressible as a
   `DiagonalOperator`-style "act-on-this-axis-of-this-space" once spaces carry axis tags.
   None is irreducibly flat in a way a typed-field algebra couldn't host. (§6 details.)

4. **The `BoundaryMomentField` gap is real and small.** `HarmonicMomentField`
   (`FluxRole, MomentField(BulkField)`) lets `TimedFullField` carry bulk moments TODAY,
   but there is NO `BoundaryMomentField` under `BoundaryField`. To window the boundary
   trace into moments (or to type a moment-valued boundary source) you mint one leaf
   mirroring `BoundaryFlux`'s machinery. The trace-space side already broadcasts over
   trailing axes (`trace_space.py` metric/selectors), so storage is the only real lever.
   (§2.)

---

## 1. Operator census by TRUE (Vin, Vout)

Every class implementing the operator interface (subclass of `LinearOperatorMixin`
or `OperatorSum`, or structural `LinearOperator`). The Protocol pins
`apply(x: np.ndarray) -> np.ndarray` (`operator.py:317`); SN/projection leaves
OVERRIDE with field signatures → the 7 `reportIncompatibleMethodOverride`.

### (a) Container-endomorphism TODAY (`TimedFullField → TimedFullField`)

These consume + emit the composite carrier. Their `apply` body reads
`psi.bulk`/`psi.boundary` and returns a fresh `TimedFullField`. The INNER bulk
type can change role (flux→source) at runtime, but the CONTAINER type is fixed.

| Class | file:line | apply: in→out (inner bulk/boundary) | role | block_role |
|---|---|---|---|---|
| `StreamingOperator` (L) | `sn/operator.py:212` | TFF(AngularFlux,BoundaryFlux) → TFF(AngularSourceSink, lpc.boundary) — `apply` 358-424; `apply_transpose` 426-468 | FULL | `BlockRole.FULL` (317) |
| `CollisionOperator` (C) | `sn/operator.py:511` | TFF → TFF(AngularSourceSink, BoundarySourceSink.zeros) — `apply` 592-618; `solve` 620-649 (→AngularFlux/BoundaryFlux); `apply_transpose`=apply 651-657 | BULK | `BlockRole.BULK` (585) |
| `InvertibleOperator` (L+C) | `sn/operator.py:681` (extends `OperatorSum`) | OVERRIDES apply/apply_transpose to the loss-rep walk; `.solve` = WDD sweep | FULL (derived) | from operands |
| `SNBoundaryOperator` (B) | `sn/boundary_operator.py:72` | TFF(any-bulk,BoundaryFlux) → TFF(AngularSourceSink.zeros, BoundarySourceSink) — `_apply_faces` 230-262, `apply` 264-266, `apply_transpose` 296+ | BOUNDARY | `BlockRole.BOUNDARY` (115) |
| `ScatteringOperator` (S) | `sn/scattering.py:297` | `singledispatchmethod` apply (962): TFF→TFF, AngularFlux→AngularFlux, ScalarFlux→ScalarSourceSink, ndarray→ndarray | BULK | `BlockRole.BULK` (343) |
| `FissionOperator` (F) | `sn/fission.py:103` | `singledispatchmethod` apply (222): TFF→TFF, AngularFlux→AngularFlux, ScalarFlux→ScalarSourceSink, ndarray→ndarray | BULK | `BlockRole.BULK` (134) |

⭐ **`S` and `F` are the imperfect-fit arms**: `singledispatch` apply means their TRUE
`(Vin,Vout)` is a UNION (TFF, AngularFlux, ScalarFlux, ndarray) — they straddle (a)
container-endo AND (b) under-typed field-map AND the legacy ndarray path in one method.
A field-typed algebra must pick: keep `singledispatch` (multi-arm) or split into typed
siblings. `apply_transpose` is NOT advertised on S/F (`capabilities={CAP_APPLY}`) — so the
full `(L+C−S−F−B).H` is intentionally unreachable (`operator.py:344-349` capability lattice).

### (b) Under-typed field map that is REALLY container/field→field (ndarray costume)

`apply(x: np.ndarray) -> np.ndarray` but `domain()`/`codomain()` return real
`FunctionSpace`s — a field→field map serialized as ndarray. With `TimedFullField`
as carrier these become container-ENDOMORPHIC (TFF→TFF, inner field differs at runtime).

| Class | file:line | TRUE inner map | domain→codomain |
|---|---|---|---|
| `MomentProjection` (Π) | `numerics/projection.py:241` | AngularFlux-values `(N,…)` → HarmonicMomentField-values `(L+1,2L+1,…)`; `apply` 424-432 = ONE einsum on axis 0; `apply_transpose` 434-482 | `angular_ordinate (N,)` w/w_n (404-422) → `SphericalHarmonicSpace (L+1,2L+1)` w/g_C (380-401) |
| `HarmonicMomentReconstruction` (R) | `numerics/projection.py:486` | HarmonicMomentField-values → AngularFlux-values; `apply` 589-598 = ONE einsum | (carries no domain/codomain props; CAP_APPLY only) |
| `ProjectionOperator`/`GalerkinProjection`/`PetrovGalerkinProjection`/`ReconstructionOperator` | `numerics/projection.py:108/143/184/208` | ABCs — no impl | — |
| `AngularAverageOperator` (G_diff, white BC) | `sn/angular_operator.py:37` | angular `(N,…)`-axis-0 Lambertian average → broadcast; CAP_APPLY only (not self-adjoint unweighted) | — |
| `IncomingSourceOperator` | `sn/angular_operator.py:230` | per-face inflow source builder (boundary trace) | — |
| `LegendreMomentScattering` | `sn/scattering.py:174` | per-ℓ `R·Λ·M`-style scattering kernel; CAP_APPLY | — |
| `ReducedStreamingOperator` / `StreamingTerms` | `geometry/reduced_operator.py:352/183` | curvilinear reduced-streaming τ/α carrier (1-D); structural `LinearOperator` | — |

⭐ **`MomentProjection`/`R` are the canonical "field→field map wearing an ndarray
costume"** the task names. Once `TimedFullField` can carry a `HarmonicMomentField` bulk
(it can — §2) AND a moment boundary (it can't yet — `BoundaryMomentField` gap), Π/R
become container-endomorphic TFF→TFF maps whose inner bulk type changes
(AngularFlux↔HarmonicMomentField). This is the seam where the windowed-moment SI
(Phase 5a/5c, `_MomentWindowedResolvent`) already lives at the ndarray level.

### (c) Genuine sub-field flat-vector PRIMITIVES (numerics/operator.py) — the KEY question

For EACH: what axis of what field it operates on, who uses it, irreducibly-flat verdict.

| Primitive | file:line | Operates on | Production user | Irreducibly flat? |
|---|---|---|---|---|
| `PermutationOperator` | 961 | ONE tagged numpy `axis` (gather `np.take`); reflective uses **axis=0 = ordinate axis** of a per-face trace slot | `boundary_realizer.py:199` (reflective `B`), folded `& IdentityOperator` into a 2-factor `TensorProductOperator`; `quad.reflection_index(axis)` is an involution | NO — "permute the ordinate axis of the trace field". A typed `PermuteAxis(space, axis)` hosts it. |
| `IncomingOrdinateMaskTensor` | 1065 | tagged `axis` (zeroes inflow ordinate indices); realizer uses **axis=0** | `boundary_realizer.py:176` (vacuum `B`), folded `& IdentityOperator` | NO — "zero inflow ordinates on the trace's ordinate axis". |
| `PeriodicWrapOperator` | 1155 | identity-copy (angular passthrough today) | `boundary_realizer.py:252` (periodic `B`), folded `& IdentityOperator` | NO — currently a typed-no-op; reserved for spatial pushforward (follow-up). |
| `DiagonalOperator` | 1417 | tagged `axis` broadcast-multiply (the §9 `AngularWeightMatrix` W) | **NO SN-production instantiation** — tests/derivations + `from_measure`; it IS the angular-quadrature-weight primitive | NO — "multiply axis n by w_n". The metric `apply_metric` already does this inline on spaces. |
| `RankOneOperator` | 1535 | tagged `axis` outer-product `|ℓ⟩⟨r|`; fission uses **axis=0 = group axis** | `fission.py:218` (`F.kernel` = `RankOneOperator(χ,νΣf,axis=0) & IdentityOperator`) | NO — "rank-1 on the group axis (χ⊗νΣf)". |
| `TensorProductOperator` | 1203 | per-axis sequential apply of factors | `fission.py:219`, every realized `B` (`& IdentityOperator` fold), `RankOneOperator & I` | NO — it IS the typed per-axis composition; pairs with `TensorProductSpace`. |
| `SumOfTensorProductsOperator` | 1324 | sum of TP summands (§15 streaming/scattering canonical form) | structural type; no live production sum yet | NO — typed separable-sum. |
| `IdentityOperator` | 877 | echoes input (any type) | ubiquitous (TP folds, `ZeroOperator` companions) | type-agnostic by design. |
| `ZeroOperator` | 899 | `codomain_zero(x)` builds the codomain's zero (FLUX→SOURCE honesty) | within-group zero-fission slot (`operator.py:932-940`) | type-agnostic; `codomain_zero` is the pre-#208 typed-output hook. |

⭐ **KEY OPEN-QUESTION VERDICT: NONE of the (c) primitives is irreducibly flat.**
Every one acts on a SINGLE named tensor axis (ordinate / group) of a field's `.values`,
and is already wrapped `& IdentityOperator` into a `TensorProductOperator` so the OTHER
axes broadcast. The flat-ndarray signature is a CONVENIENCE (numpy broadcasting handles
the untyped axes), not a necessity. The campaign target: give these primitives an
"act-on-axis-`a`-of-`FunctionSpace`-`V`" typed face (the `TensorProductSpace.find_factor`
machinery already lets a composed space recover which axis is the ordinate/group factor).
The honest constraint: they operate on `.values` ndarrays INSIDE a field, so they need
the field to expose its axis layout — which `BulkField`/`TraceSpace` already do via shape.

### (d) scipy serialization adapter

| Symbol | file:line | Role |
|---|---|---|
| `as_scipy_linop` | `operator.py:1661` | wraps `op.apply`/`apply_transpose` into `spla.LinearOperator`. **ZERO prod callers** (tests only). |
| inline `spla.LinearOperator((n,n), matvec=A_matvec)` | `iteration.py:755` (Krylov), `760-767` (precond M) | the LIVE serialization site — see §4. |
| `_AdjointOperator` | `operator.py:580` | not serialization — the metric-correct G-adjoint wrapper `A†=G⁻¹AᵀG` (reads domain/codomain `apply_metric`/`apply_inverse_metric`). |

---

## 2. Field hierarchy + the `BoundaryMomentField` gap

### MRO (from `numerics/field.py` + `transport/fields/_bases.py`)

```
Field (ABC, L1)                         numerics/field.py:143
  values:NDArray, space:FunctionSpace; UNITS:ClassVar[Unit]; frozen,eq=False,kw_only
  dunders: __add__/__sub__/__neg__/__mul__/__rmul__/__truediv__ (262-280)
  _check_partner (235): class-identity + space-equality gate
  _from_balance (285): the ONLY sanctioned class-transition (→ residual leaves)
  zeros (220), copy (368), linf/l2/inner_product (347-366)
 ├─ BulkField (ABC, L2)                  _bases.py:97   + mesh:SNMesh, ng, _phase_space_shape (abstract),
 │   │                                                   _compose_spatial_moments (157), spatial_moments_per_axis (231)
 │   ├─ AngularField (ABC)               _bases.py:267  + N, (N,ng,*spatial), _SPACE_NAME, from_mesh/zeros_on
 │   │   ├─ AngularFlux                  transport/fields/angular_flux.py   (FluxRole)
 │   │   └─ AngularSourceSink            transport/source_sinks/…           (SourceRole)
 │   │   └─ AngularResidual              transport/residuals/…              (ResidualRole)
 │   ├─ ScalarField (ABC)                _bases.py:372  + (ng,*spatial), _SPACE_NAME
 │   │   ├─ ScalarFlux / ScalarSourceSink / ScalarResidual
 │   └─ MomentField (ABC)                _bases.py:460  family marker (shape leaf-specific)
 │       └─ HarmonicMomentField          harmonic_moment_field.py:113 (FluxRole, MomentField)
 │                                        shape (L+1,2L+1,ng,*spatial[,sm**d]); UNITS=SCALAR_FLUX_UNITS (180; see [[harmonic-moment-field-units-convention]])
 └─ BoundaryField (ABC, L2)              _bases.py:481  + mesh, TraceSpace contract (space.layout),
     │                                                   layout/face_view/face_views, from_mesh/zeros_on/from_face_arrays
     ├─ BoundaryFlux                     transport/fields/boundary_flux.py   (FluxRole)
     ├─ BoundarySourceSink               transport/source_sinks/…            (SourceRole)
     └─ BoundaryResidual                 transport/residuals/…               (ResidualRole)
```

Role mix-ins (`FluxRole`/`SourceRole`/`ResidualRole`) carry `UNITS` + role identity; the
storage×role×locus grid is documented at `_bases.py:12-34`.

### The grid + the gap

```
            Flux              Source/Sink            Residual          Displacement
Angular     AngularFlux       AngularSourceSink      AngularResidual   (FluxDisplacement, #208)
Scalar      ScalarFlux        ScalarSourceSink       ScalarResidual
Moment      HarmonicMomentF.  —                      —                  ← bulk moments EXIST
Boundary    BoundaryFlux      BoundarySourceSink     BoundaryResidual
BdyMoment   ✗ MISSING         ✗                      ✗                  ← THE GAP
```

**`BoundaryMomentField` does not exist.** Every existing `BoundaryField` subclass:
`BoundaryFlux`, `BoundarySourceSink`, `BoundaryResidual` (3 leaves, all scalar-per-face).

### What a `BoundaryMomentField` must mirror from `BoundaryFlux`/`BoundaryField`

(`BoundaryField` base provides MOST of it — a moment leaf inherits, like
`HarmonicMomentField` inherits from `MomentField(BulkField)`):

- **Storage**: single flat backing `(layout.total_size,)`; `__post_init__`
  (`_bases.py:504-526`) asserts `space is TraceSpace` carrying a `FaceLayout` +
  `values.shape == (total_size,)`. A moment leaf needs the layout to allocate a moment
  axis per face slot — TODAY `TraceSpace`/`FaceLayout` are scheme-blind/scalar-per-face
  (see [[issue-251-trace-space-widening]] §"trace is purely SCALAR-per-face today";
  the single lever = `geometry.py:boundary_face_layout` appends a face-tail to the slot
  shape; `trace_space.py` UNCHANGED, metric+selectors broadcast over trailing axes).
- **mesh-trace binding**: `mesh:SNMesh` + cross-mesh `_check_partner` (530-549) +
  layout-identity check.
- **dunders**: inherited from `Field` (same-class/same-space), + the `BoundaryField`
  mesh/layout guard. A moment leaf would add an `L`-match like
  `HarmonicMomentField._check_partner` (harmonic_moment_field.py:213).
- **factories**: `zeros_on` (623), `from_mesh` (591), `from_face_arrays` (642) — all
  `cls`-based, so they construct the subclass for free; `zeros_on` is what
  `TimedFullField.zeros(boundary=…)` calls.
- **`UNITS`**: per [[harmonic-moment-field-units-convention]] a stored SH moment carries
  SCALAR-flux units (no-prefactor SH). A boundary-moment flux leaf would likely mirror
  `BoundaryFlux.UNITS` (angular-flux units on the trace) OR scalar, depending on whether
  it's a per-face angular-moment or a transverse-spatial-moment carrier — DECIDE per
  use-case (windowed-trace moments vs #251 Leg-B transverse face-slope are different).

**Where it lives**: `orpheus/transport/fields/boundary_moment_field.py` (sibling of
`harmonic_moment_field.py`), `BoundaryMomentField(FluxRole, BoundaryField)` —
OR, if it carries the #251 transverse-spatial-moment tail rather than an angular-moment
axis, it may just be `BoundaryFlux` widened (no new class; see [[issue-251-trace-space-widening]]).
⭐ **Open design fork**: is the missing leaf a (i) per-face ANGULAR-moment trace (the
boundary twin of `HarmonicMomentField`, for windowing the trace into SH moments) or
(ii) a per-face TRANSVERSE-SPATIAL-moment trace (#251 Leg B, the boundary face-slope)?
These are DIFFERENT axes. The campaign should name which it needs.

---

## 3. P3.6 — `DirectSumSpace` / `Field`-promotion status

### What exists vs deferred

- **`DirectSumSpace`: does NOT exist** (verdict §0.1). The only direct-sum machinery that
  ships is `FunctionSpace.__mul__` (tensor product, `space.py:255`) → `TensorProductSpace`
  (334) and `DualSpace` (474). The direct-sum dunder `S + T` is listed as FUTURE
  (`space.py:54-59` "S + T — direct sum; NOT shipped in 9.6").
- **`FullFieldSpace`** (the bulk⊕trace composite for the G-adjoint) DOES exist and is the
  closest prior art — `SNMesh.full_field_space` (referenced `operator.py:351,356`); it
  overrides `apply_metric`/`apply_inverse_metric` to apply a PER-BLOCK metric (volume on
  bulk, `|Ω·n|·w` on trace). It is the metric-bearing direct-sum space for `L.H`/`B.H`.
  ⭐ **`FullFieldSpace` is the seed of `DirectSumSpace`** — P3.6 likely generalizes it
  (find it via `query("FullFieldSpace")` / `query("full_field_space")`).

### What `TimedFullField` currently provides vs what `Field` ABC demands

`TimedFullField` (`timed_full_field.py:123`, `@dataclass(frozen=True,kw_only=True)`,
**NOT** a `Field` subclass) HAS:

| TimedFullField member | line | `Field` ABC counterpart |
|---|---|---|
| `bulk: BulkField`, `boundary: BoundaryField`, `_history`, `history_depth` | 176-179 | `Field` wants `values:NDArray` + `space:FunctionSpace` — TFF has NEITHER (it has two sub-fields + history) |
| `zeros(bulk:type, boundary:type, mesh)` | 183 | `Field.zeros(space, **fields)` (field.py:220) — different signature |
| `_check_partner` (container-level) | 260 | `Field._check_partner` (class+space) — TFF's is class-only (delegates member checks to leaves) |
| `__add__`/`__sub__`/`__neg__`/`__mul__`/`__rmul__`/`__truediv__` | 282-325 | matches `Field` dunders (propagate to members) ✓ |
| `advance`/`at_lag`/`history_length` | 329-427 | NOT on `Field` — history is TFF-specific (the "Timed") |
| `to_flat`/`from_flat`/`copy` | 439-519 | `copy` ✓; `to_flat`/`from_flat` are the ravellable protocol (NOT `Field` members) |
| (missing) `values`, `space`, `UNITS`, `linf`/`l2`/`inner_product` | — | `Field` requires `values`/`space`/`UNITS`; `l2`/`inner_product` need a space metric |

**What P3.6 must do to promote**: give `TimedFullField` a `space: DirectSumSpace` whose
`.shape` = bulk⊕boundary flat size, a `values` view = `to_flat()` (or a `DirectSumStorage`
backing both sub-arrays), a `UNITS` (composite — likely raises / is the bulk's), and route
`l2`/`inner_product` through the direct-sum metric (already prototyped in `FullFieldSpace`).
The docstring (`timed_full_field.py:50-52`) promises the public API stays the same →
non-breaking. The deferral reason (`:48`): `[[feedback_unify_after_two_instances]]` — wait
for the SECOND direct-sum use case (kinetics `flux ⊕ precursors`) before generalizing.

### Plan/issue scoping P3.6

- `timed_full_field.py:95` references `.claude/plans/depth_b_field_on_function_space.md`
  §3.8 (Container architecture) — the authoritative design doc.
- `docs/architecture/layering.rst:77,283` mention P3.2–P3.6 structural steps and a
  layering relationship that "dissolves under P3.6".
- No OPEN GH issue is dedicated to `DirectSumSpace`/P3.6 (grep found none) — file one if
  the campaign commits to it.

---

## 4. The scipy serialization mechanics (the typed-adapter target)

The ONE live serialization site is the Krylov `A_matvec` closure (NOT `as_scipy_linop`).

### Field → flat → field round-trip (the ravellable protocol)

- **Protocol detection** (`iteration.py:178-188`): `_is_ravellable(x)` = has `to_flat()`
  instance method + `from_flat(flat, template)` classmethod. Canonical instance =
  `TimedFullField`. Bare ndarrays fall through to `np.ravel`/`reshape`.
- **`_ravel`/`_unravel_like`/`_zeros_like`/`_l2_norm`** (`iteration.py:191-223`): the four
  protocol helpers. `numerics` MUST NOT import `transport` (L1↛L2) — so the protocol is
  DUCK-TYPED, not an ABC import. This decoupling is load-bearing.
- **`TimedFullField.to_flat`** (`timed_full_field.py:439-458`):
  `np.concatenate([bulk.values.ravel(), boundary.values])`. Boundary is already 1-D
  (`(layout.total_size,)`); bulk is reshaped. History NOT packed.
- **`TimedFullField.from_flat(flat, template)`** (`460-503`): splits at
  `n_bulk = template.bulk.values.size`, reshapes the bulk slice to
  `template.bulk.values.shape`, `dataclasses.replace(template.bulk, values=…)` (+ boundary),
  empty history. **`template` is the structural source of truth** (shapes/spaces/meshes).

### Where shape `n` comes from + the matvec

`iteration.py:695-755` (`KrylovAcceleration.solve`):
- `b = _ravel(q_ext)` (727); `n = b.size` (728).
- `solution_template` (737-741): the FLUX space (from `initial_guess` when ravellable,
  else `q_ext`) — B.5.2 fix: ψ lives in the SOLUTION space, not q_ext's source space.
- `A_matvec(psi_flat)` (743-753): `psi = _unravel_like(solution_template, psi_flat)`;
  `out = self.L.apply(psi)`; `for g in gains: out = out - g.apply(psi)`; `return _ravel(out)`.
  ⭐ This IS the honest `(L+C − S − B)·ψ` matvec — composed at the FIELD level (typed
  `apply`/`−`), serialized ONLY at the `_ravel` boundary.
- `A_scipy = spla.LinearOperator((n,n), matvec=A_matvec)` (755); preconditioner M (757-769).

### What a "typed serialization adapter" needs

A `as_scipy_linop_typed(op, template)` (or a `FieldSerializingLinearOperator`) that:
1. probes `n = template.to_flat().size` (no manual shape arg — kills the `as_scipy_linop`
   shape-guessing footgun at `operator.py:1679-1683`).
2. wraps `matvec = lambda flat: op.apply(template.from_flat(flat, template)).to_flat()`.
3. wraps `rmatvec` from `apply_transpose` when advertised (mirroring `as_scipy_linop:1706`).
4. lets the SI/Krylov drivers DROP the inline `_ravel`/`_unravel_like` closure — operators
   then NEVER touch flat ndarrays; the adapter owns the boundary. This is the concrete
   "operators never carry ndarray" deliverable. (Note the `Generic[V]` change
   [[issue-226-operator-generics-map]] is a PREREQUISITE-or-parallel typing move.)

---

## 5. Retype / consumer surface (blast radius)

### Annotations that change (`LinearOperator`/`LinearOperatorMixin`/`OperatorSum`/`q_ext: np.ndarray`)

- **`q_ext: np.ndarray`** — only 2 PRODUCTION sites, both in `iteration.py`:
  `SourceIteration.solve` (450) and `KrylovAcceleration.solve` (697). Both ALREADY accept
  typed flux via the ravellable protocol at runtime; the annotation is stale. Retype to a
  `Vector`/`Ravellable` protocol (not `Field` — `TimedFullField` isn't a `Field`; same
  bound conclusion as [[issue-226-operator-generics-map]]).
- **`LinearOperator`/`LinearOperatorMixin`/`OperatorSum` annotations**: the full surface
  is the §1 census + the composer family (`operator.py:660/746/819/877/899/...`) + the
  iteration drivers' `L: LinearOperator`, `*gains: LinearOperator` (`iteration.py:347/354,
  658/659`), the estimators (`_default_keff_estimator` 271, `_default_production_estimator`
  250). See [[issue-226-operator-generics-map]] §"consumer scope" for the exhaustive
  `Generic[V]` retype list — that memo IS the annotation-blast-radius map; THIS campaign
  consumes it.

### Callers of the under-typed projection/reconstruction operators

`MomentProjection`/`HarmonicMomentReconstruction` production consumers (grep §findings):
`sn/scattering.py` (the `R·Λ·M` aniso source), `sn/solver.py`, `sn/operator.py`,
`sn/loss_representation.py` (the windowed-moment path), `numerics/spaces/
spherical_harmonic_space.py`, `numerics/basis/spherical_harmonic_basis.py`,
`transport/fields/harmonic_moment_field.py`. The windowed-moment resolvent
(`_MomentWindowedResolvent`, Phase 5a/5c) is where Π/R run at the ndarray level INSIDE
the SI loop — the retype's highest-value target (it's the moment SI that a
`HarmonicMomentField`-carrying `TimedFullField` would type end-to-end).

---

## 6. Related issues / work — done vs open (what this campaign depends on / subsumes)

| # | title | state | relationship |
|---|---|---|---|
| **#208** | Operator-algebra completeness: Bulk/Full/BoundaryOperator + Source/Residual as FullField (Wave O) | **CLOSED** | LANDED the block-role typing, `FullFieldSpace`, G-adjoint, typed Source/Residual leaves, the `FluxDisplacement` affine gate, box-7 `evaluate_residual` (`solver.py:225`). This campaign BUILDS ON #208 — the typed field leaves + composite carrier are its substrate. See [[project-wave-o-operator-algebra]] (the ⭐ block has the full #208 record). |
| **#2** | SN: Diffusion Synthetic Acceleration (DSA) | OPEN | The consistent-DSA low-order correction CONSUMES the typed residual `evaluate_residual` (`solver.py:225-271`); in-algebra `A_diff = L+C−S`. A field-typed algebra is the substrate DSA wants. |
| **#200** | SN: block-inverse preconditioner for Krylov on the typed AngularFlux algebra | OPEN | Directly downstream — block-inverse precond on the typed algebra; the `M_scipy` precond closure (`iteration.py:760-769`) is where it lands. |
| **#226** | pyright: ~691 real type findings… | OPEN | The `Generic[V]` typing change ([[issue-226-operator-generics-map]]) is the TYPE-LEVEL prerequisite/parallel to this campaign's RUNTIME field-typing. NOT the same work — #226 makes the annotations sound; this campaign makes the algebra speak fields. |
| **#251** | TraceSpace widening (Leg B of #247) | (filed; see [[issue-251-trace-space-widening]]) | The trace-moment-axis storage lever — directly relevant to the `BoundaryMomentField` gap (§2). The transverse-spatial-moment trace is one of the two `BoundaryMomentField` forks. |
| **TransportState** | (renamed `TimedFullField` 2026-05-28) | done | `TransportState` IS `TimedFullField` — the rename happened; [[project-transport-state-container]] is the original decision memo. No separate `TransportState` exists. |

### Anti-overlap with [[issue-226-operator-generics-map]]

That memo = the NARROW `Generic[V]` parameterization (TypeVar bound = UNBOUNDED, no mixed
ndarray/field `OperatorSum` in prod, pure typing change, singledispatch S/F are the
imperfect arms). THIS memo = the BROADER campaign (BoundaryMomentField mint,
DirectSumSpace/P3.6 promotion, the flat-ndarray serialization-adapter demotion,
irreducible-flat verdict on the primitives). They share the S/F-singledispatch and
"no-Field-bound" findings; they do NOT share the gap/promotion/serialization analysis.
```

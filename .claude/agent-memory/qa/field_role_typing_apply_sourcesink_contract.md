---
name: field-role-typing-apply-sourcesink-contract
description: Re-checkable role-typing contract for SN operators — .apply bulk=AngularSourceSink, .solve bulk=AngularFlux, boundary=BoundaryFlux/BoundarySourceSink; plus the A2D-1 source-hash-pin update procedure. Re-verify on any apply/solve/SourceSink/operator-output change.
metadata:
  type: feedback
---

The #208 / field-role-typing work (MERGED, origin/main `8c2f355`→`9504146`)
makes operator output **roles** part of the type, not just `.values`. This
is a durable verification contract: re-check it whenever an operator's
`.apply`/`.solve`/matvec output type, or a `*SourceSink`/residual type, is
touched.

**The role contract (the substitution-error tripwire):**

- `.apply` / matvec output bulk → `AngularSourceSink`
  (`orpheus.transport.source_sinks.angular_source_sink`) — it is a
  rate-density source/sink, NOT a flux.
- `.solve` output bulk → `AngularFlux` (it IS the flux ψ).
- boundary fields → `BoundaryFlux` (`transport/fields/boundary_flux.py`)
  or `BoundarySourceSink` (`transport/source_sinks/boundary_source_sink.py`)
  by the same apply-vs-solve logic. (The earlier "boundary retype deferred"
  note is RESOLVED — #208 landed the boundary + residual retyping.)
- typed equation residual → `AngularResidual` / `BoundaryResidual`
  (`transport/residuals/`), minted via `from_balance`.

**Why (elegance Pattern 3 + make-illegal-states-unrepresentable):** a
source-density and a flux carry different physical roles even at identical
`.values`; the type label makes apply/solve role-confusion *unspellable*.
`AngularSourceSink.from_mesh` and `AngularFlux.from_mesh` produce
bit-identical `.values` — the retype is **behavior-neutral**, a type-label
change, NOT a numerical change. So the QA gate is: assert the role TYPE, and
separately assert `.values` are unchanged (bit-identity inheritance).

**Test-migration rule (substitution-error signature to guard):** flip
`isinstance(<apply-out>.bulk, AngularFlux)` → `AngularSourceSink` ONLY for
operator `.apply`/matvec outputs. NEVER flip a `.solve` output, nor a
genuine flux input (`TimedFullField.zeros(bulk=AngularFlux, ...)`,
`AngularFlux.from_mesh(...)`, `integrate_angular()` results). A test class
that mixes both (e.g. `tests/sn/operators/test_collision_operator.py`
`TestApply` vs `TestSolve`) needs per-assertion judgment — `C.solve` shape
test stays `AngularFlux`. Mis-flipping a `.solve` assertion to SourceSink is
the exact failure this contract exists to catch.

**A2D-1 source-hash pin** —
`tests/sn/operators/test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin`
pins `sha256(inspect.getsource(StreamingOperator._apply_2d_cartesian))`. ANY
source edit to that method (even a from_mesh retype or an added import line)
trips it. The pin is intentionally brittle — a tripwire, not a soft
constraint. **Update procedure:** (1) READ the method and confirm
behavior-neutrality FIRST; (2) recompute the hash with the test's own
method; (3) replace `EXPECTED_SHA256`; (4) ADD a dated history-log comment
entry — do NOT weaken or delete the pin.

`tests/sn/operators` is a fast tier (~6s under `-O -n 6`); its pre-existing
conflicting-marker warnings (`foundation`+`l1` / `foundation`+`l0`) are
noise, not failures.

**Affine-gate test-migration playbook (`FluxRole` leaves — `AngularFlux`,
`ScalarFlux`, `HarmonicMomentField`, `BoundaryFlux`; mixin
`transport/fields/_flux_role.py`).** When a `FluxRole` leaf migration breaks
stale `+`/`-` tests, the new contract is fixed and re-checkable:

- `flux + flux` → `TypeError` matching `"affine_combination"` (the message
  names the legal alternative). Canonical mirror:
  `tests/transport/fields/test_angular_flux.py::TestAlgebra::test_flux_add_flux_forbidden_torsor_allowed`.
- `flux − flux` → mints the sibling displacement (`_DISPLACEMENT_CLS`:
  `AngularDisplacement` / `ScalarDisplacement` / `MomentDisplacement` /
  `BoundaryDisplacement` in `transport/displacements/`), NOT another flux;
  values = `self.values − other.values`. Assert the displacement type AND
  `not isinstance(c, <FluxLeaf>)`.
- Torsor update `flux + (flux − flux) → flux` is the legal "recombine".
- Partner-compat checks (cross-mesh "mesh-bound", cross-L "equal space",
  cross-class "same-class") moved OFF `+` ONTO `__sub__` — the flux+flux
  gate fires FIRST on same-class `+`, so retarget those assertions to `−`.
- Operator **linearity** on an affine domain CANNOT be spelled
  `A.apply(α·u + β·v)` (the combined input is flux+flux). Reformulate at
  the `.values` array level: set `combined.bulk.values[...] = α·u.values +
  β·v.values`, apply, compare to `α·A.apply(u).values + β·A.apply(v).values`
  (`TimedFullField.__add__` delegates to `bulk + bulk`, so it inherits the
  gate). Keep a `pytest.raises(TypeError)` on the illegal spelling so the
  gate stays pinned. This preserves the /W-projection + bulk-boundary
  convention-drift catch.
- Decomposition completeness (`iso + aniso == self`) is array-level:
  `iso.values + aniso.values == self.values` (disjoint slices), with a
  `pytest.raises(TypeError)` pinning that the flux+flux spelling is illegal.

**Realized-BC tensor-product drill (`tests/geometry/test_bound_compat.py`).**
Post-Wave-T.1/D-B+1 a realized boundary op's `.inner` is a
`TensorProductOperator` (`orpheus.numerics.operator`, factors at `.ops`),
NOT the bare angular primitive. The §16A.10 `B = G_patch ⊗ K_omega ⊗ K_g`
form puts the angular operator at `.inner.ops[0]` (vacuum →
`IncomingOrdinateMaskTensor`, reflective →`PermutationOperator`; albedo==1.0
takes the bare-permutation TP fast path) and `IdentityOperator` at `ops[1]`.
Strongest assertion: `isinstance(.inner, TensorProductOperator)` AND
`isinstance(.inner.ops[0], <AngularPrimitive>)` AND the `== "vacuum"/
"reflective"` kind tag. Bit-identical to the legacy single-axis form
(IdentityOperator fold is a no-op) — structure changed, semantics did not.

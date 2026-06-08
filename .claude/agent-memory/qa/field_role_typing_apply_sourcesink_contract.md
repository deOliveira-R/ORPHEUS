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

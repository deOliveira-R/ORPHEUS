---
name: SN reshape Issue 9.6 closeout
description: Phase B of curried-swimming-prism plan. FunctionSpace + RegistryMixin + operator dunders + BoundaryOperator structural lift + CellUpdateBase + Mesh.volume_measure
type: project
---

# Issue 9.6 closeout — operator-algebra dunders + auto-registration

**Branch**: `refactor/sn-operator-algebra`
**Plan**: `/home/vscode/.claude/plans/curried-swimming-prism.md`
**Date**: 2026-05-10

## What shipped (B1-B9)

- **B1**: `orpheus/numerics/space.py` (NEW, 283 LOC) — `FunctionSpace`
  primitive with `(name, shape)`-keyed identity; `inner_product` /
  `norm` (weighted L² with weights metadata, falls back to Euclidean);
  factory functions `angular_flux_space`, `scalar_flux_space`,
  `boundary_trace_space`. Module docstring anticipates Grand Report
  v3 §5.3 hierarchy (MeshFunctionSpace, TraceSpace, etc.).
- **B2**: `orpheus/numerics/registry.py` (NEW, 174 LOC) —
  `RegistryMixin` per Grand Report v3 §4. `key=` class-creation kwarg;
  `__init_subclass__` auto-registration; `create()` factory with
  available-keys-listed `KeyError`. **Critical fix**: tolerates
  `@dataclass(slots=True)` replacement pattern via qualname+module
  identity check (the slots decorator creates a new class object so
  the registered entry is the original; we update on second
  init-subclass call when key is None and existing matches qualname).
- **B3**: `orpheus/numerics/operator.py` (+200 LOC) —
  `IncompatibleOperatorComposition(ValueError)`; `domain`/`range`
  properties on Protocol + Mixin (default None for backward compat);
  `__call__` (alias for `apply`, `*args, **kwargs` so 2-arg BC apply
  composes); `__pow__(n: int)` (n=0 → IdentityOperator, n>=1 →
  repeated `@`); `adjoint()` + `.H` property (weight-aware Hilbert
  adjoint via `_AdjointOperator` wrapper); composition compatibility
  checks on `OperatorSum`/`OperatorProduct`/`ScaledOperator` (skip when
  either operand has None for domain/range — backward-compat);
  `__repr__` showing class/domain/range/caps; OperatorSum/
  ScaledOperator forward `*args, **kwargs` so BC composition works.
- **B4**: `orpheus/numerics/measure.py` (+69 LOC) — `integrate`
  accepts `np.ndarray` overload (load-bearing for B5 wiring);
  `__call__`/`__iter__`/`__len__`/`__getitem__`/`__repr__`.
- **B5**: `orpheus/geometry/mesh.py` (+58 LOC) — `Mesh1D.volume_measure`
  and `Mesh2D.volume_measure` properties. `space="spatial_R1"` /
  `"spatial_R2"` string tags.
- **B6**: `orpheus/geometry/boundary.py` (+~165 LOC, ~−100 LOC for
  retired apply_to_incoming) — `BoundaryOperator` lifted from
  Protocol to concrete ABC inheriting `(LinearOperatorMixin,
  RegistryMixin, ABC)`; concrete subtypes self-register via `key=`
  ("vacuum", "reflective", "white", "periodic", "albedo", "mixed").
  `apply_to_incoming` DROPPED entirely; `apply(psi_out, quadrature)`
  is canonical. `SpecularBoundaryOperator.apply_transpose` shipped
  (involution permutation, capability bumped to {apply, apply_transpose}).
  All 28 call sites in `sn/sweep.py` and `sn/operator.py` migrated.
- **B7**: `orpheus/sn/spatial/cell_update.py` (+70 LOC) —
  `CellUpdateBase(RegistryMixin, ABC)` ships alongside the existing
  `CellUpdate` Protocol (Protocol stays for structural typing).
- **B8**: `orpheus/sn/spatial/diamond.py` (+2 LOC) — `DiamondDifference`
  inherits `CellUpdateBase` with `key="diamond_difference"`.
- **B9**: NO sites migrated. After inspecting all candidate sites in
  `cp/solver.py`, `sn/operator.py`, `sn/solver.py`:
  - `cp/solver.py:350,365,368` — per-cell access (`V[i]` /
    `V[:, np.newaxis] * flux`). NOT integration.
  - `cp/solver.py:500,503,550,554` — element-wise multiplication
    `volumes * Q_g`. NOT integration.
  - `cp/solver.py:589,611,617` — `V = self.volumes` followed by
    `np.sum(... * V, axis=...)` over (cells, groups). Sums over
    BOTH axes; volume_measure's 1D-cell-axis design is awkward fit.
  - `cp/solver.py:624,628` — per-cell `volumes[k]`. NOT integration.
  - `cp/solver.py:679-681` — subset-sum filtered by mat_id. NOT
    full-mesh integration.
  - `cp/solver.py:687-695,784,792` — pass-through to internal
    helpers; not integration consumers.
  - `sn/operator.py:905,917` — `sn_mesh.volumes` passed as data to
    inner helpers (per-cell access in
    transport_operator_matvec_*). NOT integration.
  - `sn/solver.py:233,297,302` — same shape (cells, groups);
    `np.sum(... * volume[..., None])` sums over both axes.
    Awkward fit.

  The affordance ships shipped + tested but no production site is
  cleanly amenable. The wiring step is honest: shipping
  volume_measure as the natural integration measure, with array +
  callable overloads tested, leaves the production sites for a
  follow-up campaign that re-architects them around the measure
  primitive (likely Wave 1 amendment).

## Test additions

- `tests/numerics/test_space.py` (NEW, 17 tests, 158 LOC) —
  equality/hash on (name, shape); inner product Euclidean +
  weighted; norm consistency; factory functions.
- `tests/numerics/test_registry_mixin.py` (NEW, 10 tests, 153 LOC)
  — `__init_subclass__` registration; `create()` factory;
  duplicate detection; registry-root isolation. Named
  `test_registry_mixin.py` because `test_registry.py` already
  exists for the quadrature-rule selection registry.
- `tests/numerics/test_operator.py` (+~210 LOC, 17 new tests) —
  __call__, __pow__, .H aliases, Hilbert adjoint identity (Euclidean
  + weight-aware), IncompatibleOperatorComposition, __repr__.
- `tests/numerics/test_measure.py` (+~95 LOC, 9 new tests) —
  __call__, array overload, __iter__, __len__, __getitem__, __repr__.
- `tests/geometry/test_mesh.py` (NEW, 10 tests, 168 LOC) — Mesh1D
  and Mesh2D volume_measure with constant + arbitrary value tests.
- `tests/geometry/test_boundary.py` (+~115 LOC, 8 new tests) —
  registry membership, create() factory, specular transpose
  reciprocity, self-inverse identity, OperatorSum-of-BCs matches
  MixedBoundaryOperator baseline.
- `tests/sn/spatial/test_diamond.py` (+~25 LOC, 2 new tests) —
  CellUpdateBase registry membership.

## Verification gate results

- `pytest -m regression -q` → **11/11 PASSED bit-identical**
  (763.23s ≈ 12:43).
- `pytest tests/numerics tests/geometry tests/sn/spatial -q` →
  529 passed.
- `pytest tests/numerics tests/geometry -q` → 500 passed in 0.94s.
- `sphinx-build -W docs docs/_build/html` → exit 0.
- `python -m tests._harness.audit` → 23 orphans, 36/38 ERR coverage
  (matches Phase A baseline).

## Critical findings

1. **`@dataclass(slots=True)` + `__init_subclass__` interaction**:
   The slots decorator creates a NEW class object replacing the
   original in the module namespace. Registry initially registers
   the pre-slots class (correct at registration time); but
   `from module import Foo` returns the post-slots class. The fix:
   second `__init_subclass__` call with `key=None` triggers a
   silent registry update if (a) the class has inherited
   `cls.key`, (b) the existing entry has the same qualname/module,
   (c) existing is not cls. This preserves duplicate-key detection
   for actual duplicates while tolerating slots replacements.

2. **OperatorSum/ScaledOperator must forward *args, **kwargs**:
   BoundaryOperator's `apply(psi_out, quadrature)` takes 2 positional
   args. Existing `OperatorSum.apply(self, x)` only forwarded `x`,
   so `(0.7 * specular + 0.3 * white).apply(psi, quad)` would crash.
   Fix: forward `*extra, **kwextra` in OperatorSum.apply,
   apply_transpose, and ScaledOperator.apply, apply_transpose.
   OperatorProduct deliberately NOT changed (products of BCs don't
   make sense; ``A @ B`` requires A.apply to consume B.apply's
   output, so `(specular @ specular).apply(psi)` would receive the
   inner permutation result and try to call `.apply(psi_out, quad)`
   with only 1 arg — correctly fails).

3. **Plan deviation (B6 SpecularBoundaryOperator transpose)**: Plan
   said "Verify the index layout matches by reading the current
   apply body — the slicing pattern must mirror exactly." Current
   apply uses `psi_out[ref]` (basic indexer), not `psi_out[..., ref, :]`
   (positional indexing). Transpose mirrors with `psi_in[perm]`.
   The reciprocity test (`<B(psi_out), phi_in> == <psi_out,
   B^T(phi_in)>`) and self-inverse test
   (`B(B(x)) == albedo^2 * x`) both pass at 1e-13 — confirms the
   layout is correct as written.

## Out of scope (deferred)

- B9 wiring at production sites (no clean integration touchpoints —
  re-architect into a follow-up campaign).
- WhiteBoundaryOperator / PeriodicBoundaryOperator /
  AlbedoBoundaryOperator / MixedBoundaryOperator transposes.
- CellUpdateBase consumption from the sweep dispatcher (Wave D
  task).
- v3 Space hierarchy (MeshFunctionSpace, TraceSpace, etc.) — file
  is structured to admit them additively, but Phase B ships only
  `FunctionSpace` itself.

## Ergonomics notes for future agents

- `orpheus/numerics/registry.py:RegistryMixin._registry_base` MUST
  be overridden in each registry root. Default returns `cls`, so
  forgetting the override means subclasses register under whatever
  `cls` they were defined as — usually the wrong scope.
- `Mesh1D.volume_measure` does a local import of
  `orpheus.numerics.measure` to avoid circular imports.
  `orpheus.numerics.measure` is imported by lots of code; adding
  `from orpheus.geometry.mesh import ...` at the top of measure.py
  would create a cycle.
- The `_AdjointOperator` wrapper's `apply_transpose` raises
  `NotImplementedError` — taking the transpose of an adjoint is a
  weight-swap-and-divide that no current consumer needs. Add later
  when sensitivity-analysis pipelines demand it.

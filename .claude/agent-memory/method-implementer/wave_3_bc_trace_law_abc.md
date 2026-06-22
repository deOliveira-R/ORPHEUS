---
name: Wave 3 BoundaryTraceLaw ABC + 8 named errors + BoundarySource
description: BC refactor Wave 3 closeout 2026-05-11. Pure architectural ABC + Protocol + named-error catalog. boundary.py promoted to boundary/__init__.py via git mv; new submodules _base/_errors/_source. Registry isolation between BoundaryTraceLaw and legacy BoundaryOperator verified.
type: project
---

# Wave 3 — BoundaryTraceLaw ABC + BoundarySource + 8 named errors

`feature/bc-trace-law-abc` 2026-05-11. Wave 3 of the boundary-operator
refactor (`transient-giggling-cake.md`). Pure architecture: ABC +
Protocol + typed-error catalog. No constructive math (no Branch-1 /
Branch-2 split — algebra-of-record discipline doesn't apply to this
wave; bifurcation reappears in Wave 7 when concrete BCs ship the
invariant overrides).

## Path-conflict resolution

The plan's text contradicted its critical-files list (text said
`boundary.py` "stays unchanged"; list put new files under
`boundary/_base.py`). Resolved: `boundary.py` → `boundary/__init__.py`
via `git mv` (preserves history), new submodules added alongside.
Every legacy import site (`from orpheus.geometry.boundary import X`)
resolves unchanged because Python serves the package
`__init__.py` for that path. Wave 4 will split legacy concretes out
of `__init__.py` into per-BC submodules.

## Deliverables (manifest)

1. `orpheus/geometry/boundary.py` → `orpheus/geometry/boundary/__init__.py`
   via `git mv` (582 lines, content unchanged) — re-exports for the
   new Wave 3 names appended at the bottom.
2. NEW `orpheus/geometry/boundary/_errors.py` — 1 base
   (`BoundaryError(ValueError)`) + 8 typed subclasses, each carrying
   a `law: str = ""` field and the inherited `ValueError` MRO.
3. NEW `orpheus/geometry/boundary/_source.py` — `BoundarySource`
   `@runtime_checkable` Protocol + `NoSource` sentinel
   (`@dataclass(frozen=True)`) + `ConstantInflowSource(value: float)`.
4. NEW `orpheus/geometry/boundary/_base.py` — `BoundaryTraceLaw(LinearOperatorMixin, RegistryMixin, ABC)`
   with `registry: ClassVar[dict]` independent of legacy
   `BoundaryOperator.registry`; `creates_sweep_cycle: ClassVar[bool] = False`;
   three first-class properties (`geometry_map`, `response_kernel`,
   `source` defaulting to `None / None / NoSource()`); five `assert_*`
   universal invariants as no-op defaults; abstract `apply(psi_out, *args, **kwargs)`;
   `realize(method_space)` raises `NotImplementedError` (Wave 5 wires
   the realiser).
5. NEW `tests/geometry/test_bc_errors.py` — 11 foundation tests
   (one per error class + base + 2 default-kwarg contracts).
6. NEW `tests/geometry/test_boundary_trace_law.py` — 17 foundation
   tests covering ABC non-instantiability, stub-concrete construction,
   `__call__` delegation, default-property contract, no-op `assert_*`,
   `realize` Wave 5 deferral, registry self-registration, and
   **registry isolation between `BoundaryTraceLaw` and legacy
   `BoundaryOperator`** (load-bearing for Wave 7's rename).
7. UPDATED `.claude/skills/vv-principles/error_catalog.md` — 8 new
   entries ERR-040 through ERR-047 (one per typed error). Each entry
   pins the failure mode (per the 6-AI-mode taxonomy), the mechanism,
   how it hides, the Wave-7 catching test that will fire
   `@pytest.mark.catches("ERR-NNN")`, and the lesson.

## Test counts

- Wave 3 NEW: 28 tests (11 in `test_bc_errors.py` + 17 in
  `test_boundary_trace_law.py`); all `@pytest.mark.foundation`.
- Regression `tests/geometry/` + `tests/numerics/`: **652/652 pass**
  in 1.05s (no regression from Wave 2 baseline of 624).
- Import-site smoke (`test_boundary.py` + `test_registry_mixin.py` +
  `test_angular_average_operator.py` + `test_snstreamingoperator.py`):
  **84/84 pass** in 1.1s.
- `tests/sn/` full directory (non-slow): exit code 0.

## Architecture decisions

- **Why two independent registries?** Wave 7 will rename concrete BCs
  from `VacuumBoundaryOperator` (legacy) to `VacuumInflow` (new); both
  inherit `RegistryMixin`. Without distinct `_registry_base()` returns,
  the same string key (`"vacuum"`) on both hierarchies would collide
  at import via the `KeyError` guard in `RegistryMixin.__init_subclass__`.
  The Wave-3 test `test_legacy_and_new_registries_are_disjoint` pins
  this — without it, Wave 7's first new concrete crashes at import.
- **Why `apply(psi_out, *args, **kwargs)` and not `apply(psi)`?**
  Transition-period contract. During Waves 7-9 the new concretes will
  delegate to the legacy `BoundaryOperator.apply` body which takes
  `(psi_out, quadrature)`; the abstract signature must accept both
  forms. Wave 10 collapses to 1-arg `(psi)` once the realiser captures
  the quadrature at construction.
- **Why `@runtime_checkable` Protocol for `BoundarySource`?**
  Sources need to be type-checked by realisers (Wave 5) without
  forcing subclassing — custom user sources (spatially-varying beam
  injection, time-dependent inflow) only need the `evaluate` method.
  Duck typing is the right contract; the runtime_checkable decorator
  makes `isinstance(s, BoundarySource)` meaningful.
- **Why `BoundaryError(ValueError)` and not `BoundaryError(Exception)`?**
  Backward compat. The legacy boundary layer raises `ValueError`
  in numerous spots; downstream code likely catches it. Making
  `BoundaryError` a `ValueError` keeps every `except ValueError`
  consumer working while letting new code match by typed subclass.
  Verified by `test_boundary_error_base_class_contract`.
- **Why no `MAIN_AGENT_AUGMENTATION` response?** The main agent's
  two concerns were both satisfied by the original brief: the
  `*args, **kwargs` stub signature is in `_StubLaw.apply` and the
  registry-disjointness contract has a dedicated test.

## Architectural deviations from brief

None of substance. The brief's `__all__` snippet used Python `*globals().get()`
which is unreliable inside a module; I rewrote as
`__all__ = [*__all__, "BoundaryTraceLaw", ...]` after the legacy
`__all__` definition. Same end result; cleaner.

## ERR-NNN allocations

| ERR     | Error class                                          | Failure mode    |
| ------- | ---------------------------------------------------- | --------------- |
| ERR-040 | IncomingOutgoingTraceClassificationError             | #5 (index)      |
| ERR-041 | VacuumAppliedToOutgoingTraceError                    | #6 (convention) |
| ERR-042 | BoundaryGeometryMapNotMeasurePreservingError         | #5 + #6         |
| ERR-043 | BoundaryResponseNotPositiveError                     | #1 (sign)       |
| ERR-044 | ReflectionNotInvolutiveError                         | #5 (index)      |
| ERR-045 | ReflectionDidNotMapInflowToOutflowError              | #5 (index)      |
| ERR-046 | SubmarkovViolationError                              | #4 (factor)     |
| ERR-047 | BoundarySourceNotOnIncomingTraceError                | #6 (convention) |

Each entry's "Catching test" line names the Wave-7 concrete BC test
that will fire the invariant. The `@pytest.mark.catches("ERR-NNN")`
decorators are NOT added in Wave 3 — they ship in Wave 7 alongside
the concrete BCs.

## Open issues / follow-ups

- **Wave 4**: split legacy concretes (`VacuumBoundaryOperator`,
  `SpecularBoundaryOperator`, etc.) out of `__init__.py` into per-BC
  submodules. Legacy import paths must continue to resolve.
- **Wave 5**: ship `BoundaryRealizer` Protocol + per-method realisers.
  `BoundaryTraceLaw.realize` will defer to a registry; remove the
  Wave-3 `NotImplementedError`.
- **Wave 7**: concrete `BoundaryTraceLaw` subclasses with invariant
  overrides that fire the 8 typed errors. Add `@pytest.mark.catches("ERR-NNN")`
  decorators on the L0 tests then.

## Commit

ONE atomic commit including the `git mv`:
`feat(geometry): BoundaryTraceLaw ABC + BoundarySource + 8 named errors (Wave 3 / C3.1-C3.3)`.

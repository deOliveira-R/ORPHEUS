# Test Suite — Improvement Tracker

Central registry for test infrastructure improvements.
See `02.Discrete.Ordinates/IMPROVEMENTS.md` for solver improvements.

## Tracking Number Format

`TS-YYYYMMDD-NNN` (TS = Test Suite)

---

## OPEN — Not Yet Implemented

### TS-20260405-001 — Reorganize tests/ by model

**Priority**: Medium | **Effort**: Moderate (rename-heavy)

Restructure `tests/` from flat layout to model-based folders:

```
tests/
├── sn/
│   ├── test_cylindrical.py
│   ├── test_spherical.py
│   ├── test_quadrature.py
│   ├── test_solver_components.py
│   └── test_cartesian.py
├── cp/
│   └── test_collision_probability.py
├── diffusion/
│   └── test_diffusion_1d.py
└── conftest.py
```

Requires updating all `pytest` commands in documentation, CLAUDE.md,
agent instructions, and any CI configuration.

### TS-20260405-002 — Add `@pytest.mark.l0` marker to all L0 tests

**Priority**: Medium | **Effort**: Small

Register an `l0` marker in `pyproject.toml` and tag all term-level
verification tests across the suite.  This enables:
- `pytest -m l0` to run only L0 tests
- Publication table generation from test docstrings
- Coverage analysis: which terms have L0 tests, which don't

Requires: identify all L0 tests across `test_sn_*.py` files and tag
their docstrings with `L0-SN-NNN` problem numbers.

### TS-20260405-003 — Sphinx L0 verification page

**Priority**: Medium | **Effort**: Moderate

Create `docs/theory/verification_l0.rst` with sections per solver
(SN, CP, etc.).  Each section lists L0 problems with:
- Problem ID (L0-SN-NNN)
- Physical setting
- Term being verified
- Expected value (analytical/hand-calc)
- Test function that verifies it

Generated from test docstrings + the error catalog.

### TS-20260403-001 — Richardson reference caching

**Priority**: Medium | **Effort**: Small

SN and MOC heterogeneous verification cases recompute Richardson
extrapolation references on every test run (~15 min for all 6 SN cases).
These should be computed once and cached (e.g., as JSON or in a
`_cached_references.py` generated module).  The cache is invalidated
when the solver code or XS library changes.

### TS-20260403-002 — MMS (Method of Manufactured Solutions) tests

**Priority**: Medium | **Effort**: Moderate

Fixed-source verification for SN, MOC, and diffusion.  Tests the
transport operator independently of the eigenvalue iteration.
See DV-20260403-007.

### TS-20260405-004 — L0 test coverage for CP solver

**Priority**: Low | **Effort**: Moderate

No `l0_cp.py` or CP-specific L0 tests exist yet.  When the CP
solver gets curvilinear improvements, L0 tests should be written
following the SN pattern.

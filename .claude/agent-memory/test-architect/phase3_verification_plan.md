---
name: phase3-verification-plan
description: Verification plan for Phase 3 of the moment-space-and-layering refactor (layered re-architecture). Covers the P3.1 import-linter test specification with FORBIDDEN_EDGES + tolerances + skeleton + expected initial failures; per-step verification gates (P3.2 numerics reorg, P3.3 transport migration + AngularFlux split, P3.4 problem/solver split, P3.5 range→codomain rename, P3.6 kinetics restructure); AngularFlux design-note CC.4 resolution; bit-identity vs principled-equivalence judgments; risk register with concurrence note on the REVISED sequencing.
metadata:
  type: project
---

# Phase 3 — Verification Plan (Moment-Space & Layering)

**Branch**: `refactor/moment-space-and-layering` (worktree:
`.claude/worktrees/moment-space-and-layering/`)
**Plan**: `.claude/plans/moment_space_and_layering_plan.md`
**Phase 1 closed**: commits 1ab6233..ea02ab5 (11 commits).
**Phase 3 sequence (REVISED)**: P3.1 → P3.5 → P3.0 → P3.2 → P3.3 →
P3.4 → P3.6.

## CRITICAL — claim layers and pillar gating (per `vv-principles`)

Phase 3 is a **package-layout refactor**. It moves modules, renames
attributes, and splits classes. It touches NO equation, NO algorithm,
NO discretization. The claim layers therefore are:

- **Foundation claims** — software invariants on the import graph
  (P3.1 import-linter), the attribute rename (P3.5), and the dunder
  shape of `transport/` types. Pillar: software contract; no theory-page
  label. Tag: `@pytest.mark.foundation`.
- **Bit-identity invariance claims** — for every mechanical move,
  the existing L0/L1/L2 verification suite stays GREEN, with results
  bit-identical to the pre-move baseline OR drifting only within the
  FP-non-associativity bound per `vv-principles` §"Bit-identity vs
  principled-equivalence". Pillar: inherited (frozen snapshots +
  L1 MMS gates + existing analytical references). No new pillar work
  is required by Phase 3 itself.
- **NO new eigenvalue, flux-shape, or convergence-order claims** are
  made by Phase 3. The existing references in L1 gates remain the
  truth-of-record; Phase 3's contract is "do not break them."

Cardinal-rule check (`vv-principles` §"1-group degeneracy"): the
existing L1 gates already contain ≥2G heterogeneous mesh-refined cases
(`tests/sn/test_mms_aniso.py`, `tests/sn/l1_analytical/...`). Phase 3
inherits this coverage; the gates checklist below pins these as STRICT
green at every commit.

---

## SECTION 1 — The P3.1 import-linter test specification

### 1.1 Test file path

`tests/test_layer_imports.py` (top-level — the test is cross-cutting;
it walks the whole `orpheus/` tree).

### 1.2 Test mechanism

Pytest test that:

1. Walks `orpheus/**/*.py` excluding `__pycache__`.
2. Parses each file's import statements via Python's `ast` module
   (NOT regex — regex misses `from ... import ...` inside `if
   TYPE_CHECKING:` blocks, multi-line imports, and conditional
   imports).
3. For each module's parent package, looks up the `FORBIDDEN_EDGES`
   target set.
4. Reports every `from orpheus.X` or `import orpheus.X` where `X` is
   in the forbidden set.
5. Parametrised one test per module — a failing test names the
   offender in the test ID.

The test is `@pytest.mark.foundation` (per `vv-principles` §"V&V
level taxonomy") — a software contract, no theory-page `:label:`.

### 1.3 Layer assignment dictionary — EVERY existing package

Reading the package contents:

| Package | Layer | Rationale |
|---|---|---|
| `numerics/` | **L1** | Per plan §P3.0 — math primitives; knows no neutrons. |
| `geometry/` | **(input)** | Per plan §P3.0 — geometry input. |
| `data/` | **(input)** | Per plan §P3.0 — nuclear-data input. |
| `transport/` | **L2** | Per plan §P3.0 — transport vocabulary, method-agnostic (created by P3.3). |
| `sn/` | **L3** | Per plan §P3.0 — one method. |
| `moc/` | **L3** | Per plan §P3.0. |
| `cp/` | **L3** | Per plan §P3.0. |
| `mc/` | **L3** | Per plan §P3.0. |
| `diffusion/` | **L3** | Per plan §P3.0. |
| `kinetics/` | **L3** (transitional) | Plan §P3.6 dissolves this package; while it exists it is L3-peer (an ODE solver). Imports `data/` only (verified). After P3.6, the contents migrate to `transport/problems/` (L2) + `numerics/solvers/` (L1); the directory empties or retires. The linter MUST tolerate `kinetics/` as L3 during the transition window and the assignment retires automatically once the package is empty. |
| `pn/` | **L3** (planned) | Plan §P3.0 names it as L3; not yet a real package. The linter declares the assignment so the package can land into a green-linter codebase without a re-spec. |
| `fuel/` | **L3 (peer)** | `fuel/solver.py` imports only `data/`. It is a thermophysics solver that depends on `data/`. Per the plan §P3.0 layer rule, the lowest layer whose vocabulary suffices is L3 — it is a domain solver (cf. `diffusion/`). Treat as a sibling of `diffusion/`. |
| `thermal_hydraulics/` | **L3 (peer)** | `thermal_hydraulics/solver.py` imports only `data/`. Same argument as `fuel/` — domain solver, L3 peer of `diffusion/`. |
| `homogeneous/` | **L3 (peer)** | `homogeneous/solver.py` imports `data/` AND `orpheus.numerics.eigenvalue.power_iteration` — a legitimate L3 → L1 import. The package is a 0-D reactor-physics reduced-order solver, L3 peer of `diffusion/`. (CC.8 retires `power_iteration`; `homogeneous/` rewires to `KEigenvalue` in P3.4 alongside `sn/`.) |
| `derivations/` | **L0 reference** | `derivations/` ships analytical / SymPy / mpmath REFERENCE solvers (`algebra-of-record` Branch 1A/1B/1C). Per `algebra-of-record` §"Structural independence", references must NOT share upstream primitives with production. However, `derivations/continuous/cases/diffusion.py:152` does `from orpheus.diffusion.solver import ...` INSIDE a function body — this is an L0-USES-L3 PATTERN (the reference USES the production solver as a black-box at a higher-level cross-check), not L3-USES-L0. The package classifies as **L0 (below L1)** with **explicit exemptions** for the L0-uses-L3 inline-imports inside functions (which are documented as benchmark cross-checks, not algebra sharing). See §1.5 exemptions. |
| `plotting.py` | **L4 (driver)** | Single-file plotting utilities — orchestration / display. L4. Treated as a special-case rule (no forbidden-edge set; consumes everything). |

### 1.4 The `FORBIDDEN_EDGES` dictionary

```python
# Layer-package assignment (frozenset for hashing in tuple keys)
L1_PACKAGES: frozenset[str] = frozenset({"numerics"})
INPUT_PACKAGES: frozenset[str] = frozenset({"geometry", "data"})
L2_PACKAGES: frozenset[str] = frozenset({"transport"})
L3_PACKAGES: frozenset[str] = frozenset({
    "sn", "pn", "moc", "cp", "mc", "diffusion",
    "kinetics",            # transitional; retires under P3.6
    "fuel", "thermal_hydraulics", "homogeneous",  # L3 peers
})
L0_PACKAGES: frozenset[str] = frozenset({"derivations"})

# FORBIDDEN_EDGES: importer package → forbidden target packages.
FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
    # L1 (numerics) imports nothing above itself.
    "numerics": L2_PACKAGES | L3_PACKAGES,

    # Input layers (geometry, data) — neither depends on transport
    # nor on any method; geometry may use numerics primitives.
    "geometry": L2_PACKAGES | L3_PACKAGES,
    "data":     L2_PACKAGES | L3_PACKAGES,

    # L2 (transport) imports L1 and inputs only.
    "transport": L3_PACKAGES,

    # L3 methods cannot import sibling L3 packages.
    "sn":                  L3_PACKAGES - {"sn"},
    "pn":                  L3_PACKAGES - {"pn"},
    "moc":                 L3_PACKAGES - {"moc"},
    "cp":                  L3_PACKAGES - {"cp"},
    "mc":                  L3_PACKAGES - {"mc"},
    "diffusion":           L3_PACKAGES - {"diffusion"},
    "kinetics":            L3_PACKAGES - {"kinetics"},
    "fuel":                L3_PACKAGES - {"fuel"},
    "thermal_hydraulics":  L3_PACKAGES - {"thermal_hydraulics"},
    "homogeneous":         L3_PACKAGES - {"homogeneous"},

    # L0 (derivations) sits BELOW L1 in the import hierarchy.
    # References must not depend on transport / methods unless via
    # the WHITELIST below (explicit benchmark uses inside function
    # bodies — NOT module-level imports).
    "derivations": L2_PACKAGES | L3_PACKAGES,
}
```

### 1.5 Tolerances (the two whitelists)

**(a) `TYPE_CHECKING` guard exemption.** Imports inside an `if
TYPE_CHECKING:` block are skipped — they exist only for type
checkers and produce no runtime edge. The parser walks `ast.If`
nodes; when the condition is `ast.Name(id="TYPE_CHECKING")`, the body
is parsed with `in_type_checking=True` and forbidden-edge violations
are suppressed for any combination of (source layer, target layer)
where the source is L1 or L2.

This implements plan §P3.1: "TYPE_CHECKING imports of L3-only types
inside L1/L2 modules are allowed (string annotations don't create
runtime edges)."

**(b) Explicit `WHITELIST` set.** A `frozenset[tuple[str, str]]` of
`(module_relative_path, target_top_level_package)` pairs that the
linter MUST pass even though the FORBIDDEN_EDGES dict would reject
them.

Initial whitelist (until each is retired in a later Phase 3 step):

```python
WHITELIST: frozenset[tuple[str, str]] = frozenset({
    # CC.8 — power_iteration legacy shim. Retires in P3.4.
    ("numerics/eigenvalue.py", "sn"),   # the function imports SN
                                        # types in its docstring
                                        # `:class:` references; the
                                        # ACTUAL module body has no
                                        # such import — this entry is
                                        # documentary placeholder
                                        # until verified at P3.1
                                        # implementation time.

    # derivations L0-uses-L3 black-box benchmark imports (inline,
    # inside function bodies — NOT module-level). These are
    # benchmark cross-checks (Branch-1-uses-production-as-black-box),
    # not algebra-sharing. Retire when each derivation moves to a
    # method-side test or to an external benchmark harness.
    ("derivations/continuous/cases/diffusion.py", "diffusion"),
    ("derivations/continuous/mms/moc.py", "moc"),
    ("derivations/continuous/sood_registry/builders.py", "cp"),
})
```

The whitelist entries are commented with their retirement trigger so a
future contributor can clean them out without reverse-engineering the
intent.

### 1.6 Expected initial failures

Direct file grep against the current codebase shows these violations.
The P3.1 commit lands the linter with these as `xfail` OR with
`WHITELIST` entries marking each as "to retire in P3.x"; the
subsequent Phase 3 commits flip the `xfail` to PASS one module at a
time. The first principled course of action is the **`WHITELIST`
path with `# RETIRE_IN_PX.X` comments**, so the actual failure set is
empty at P3.1 landing time AND the migration to-do list is visible
in the whitelist.

| Expected initial failure | Phase 3 step that retires it | Mechanism |
|---|---|---|
| `numerics/__init__.py:3` imports `from .eigenvalue import EigenvalueSolver, power_iteration` — `eigenvalue.py` itself has a docstring `:class:` reference to `orpheus.numerics.iteration.KEigenvalue` (L1-to-L1, not forbidden) but the module's CC.8 deprecation chain references SN symbols by string. NOT a real graph edge; not a violation. **No whitelist entry needed.** | — | — |
| `derivations/continuous/cases/diffusion.py:152` does `from orpheus.diffusion.solver import ...` INSIDE a function body. AST walk sees this as a module-level edge from `derivations` to `diffusion` (L0 to L3) — forbidden under the L0 rule. | P3.6 close-out (or sooner if the case file moves to a method-side test). | `WHITELIST` entry until then. |
| `derivations/continuous/mms/moc.py:96-97` does `from orpheus.moc.geometry import MOCMesh; from orpheus.moc.quadrature import MOCQuadrature` at module level — L0 to L3, forbidden. | Same as above — these MMS sources should be either method-side or import only `orpheus.geometry` primitives (verified `MOCMesh` extends `Mesh1D`). | `WHITELIST` until refactored; flagged for follow-up after Phase 3. |
| `derivations/continuous/sood_registry/builders.py:62` does `from orpheus.cp.solver import CPParams` inside a `TYPE_CHECKING` block — already in the type-checking-tolerance. **Likely no violation reported.** Confirm at P3.1 implementation. | — | TYPE_CHECKING exemption (auto-detect). |
| `numerics/spaces/`, `numerics/basis/` are NEW packages from P1.1/P1.2 — verified clean (no L2/L3 imports). | — | — |
| `numerics/iteration.py` is verified L1-clean per CC.7 — no `from orpheus.sn` at module level; cross-package flow is purely duck-typed via `_is_ravellable`. | — | — |
| `sn/`, `moc/`, `cp/`, `mc/`, `diffusion/` — sibling-L3 imports: NONE confirmed (grep `from orpheus.sn` inside `moc/` etc. returns empty across L3 packages — confirmed against current branch). | — | — |
| `cp/solver.py:43-45` imports from `orpheus.derivations.common.kernels`, `orpheus.derivations.common.quadrature`, `orpheus.derivations.continuous.flat_source_cp.geometry`. This is **L3 → L0** — REVERSE direction from L0 → L3. Per the layer rule, references sit BELOW L1; if `derivations/` is L0 (below L1), then `L3 → L0` is in fact **L3 reading from BELOW the math layer**, which violates "imports flow from more-knowledge to less-knowledge" (L0 has less knowledge than L3, so the import is allowed). **No violation.** This is an L0-reference being used by an L3 method as a (verified, structurally-independent) source of truth — the canonical Branch-1-as-reference pattern. The linter MUST permit `L3 → L0`. | — | — |

**Net expected initial failure count after the whitelist is applied:
0** (with 3 whitelist entries marked for follow-up retirement). The
P3.1 commit therefore lands the linter GREEN, and every subsequent
P3.x commit must keep it green. If a P3.x commit introduces a new
violation, the linter catches it before merge.

### 1.7 Skeleton implementation

```python
"""Import-linter: enforce the L0 / L1 / (input) / L2 / L3 layer
contract per plan §P3.0.

Layers (per plan §P3.0):
  L0  (derivations/)   — Branch-1 references (SymPy / mpmath);
                         BELOW L1 in the import hierarchy.
  L1  (numerics/)      — math primitives; knows no neutrons.
  (input) (geometry/, data/) — geometry + nuclear data.
  L2  (transport/)     — transport vocabulary; method-agnostic.
  L3  (sn/, pn/, moc/, cp/, mc/, diffusion/, kinetics/,
       fuel/, thermal_hydraulics/, homogeneous/) — one method's
       machinery; method-specific.
  L4  (plotting.py)    — orchestration; consumes everything.

Imports flow only L3 -> L2 -> L1 (and L3 -> L0 references, and any
layer -> input). Forbidden edges raise the parametrised test.

Tolerances:
  - TYPE_CHECKING imports of L3-only types inside L1/L2 modules are
    permitted (string annotations don't create runtime edges).
  - WHITELIST entries cover transitional exemptions retired by
    later Phase 3 steps; each whitelist entry carries a comment
    naming its retirement trigger.
"""
from __future__ import annotations

import ast
import pathlib
from collections.abc import Iterable

import pytest

ORPHEUS_ROOT = pathlib.Path(__file__).parent.parent / "orpheus"

# Layer assignment.
L0_PACKAGES: frozenset[str] = frozenset({"derivations"})
L1_PACKAGES: frozenset[str] = frozenset({"numerics"})
INPUT_PACKAGES: frozenset[str] = frozenset({"geometry", "data"})
L2_PACKAGES: frozenset[str] = frozenset({"transport"})
L3_PACKAGES: frozenset[str] = frozenset({
    "sn", "pn", "moc", "cp", "mc", "diffusion",
    "kinetics",
    "fuel", "thermal_hydraulics", "homogeneous",
})

FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
    "numerics":             L2_PACKAGES | L3_PACKAGES,
    "geometry":             L2_PACKAGES | L3_PACKAGES,
    "data":                 L2_PACKAGES | L3_PACKAGES,
    "transport":            L3_PACKAGES,
    "sn":                   L3_PACKAGES - {"sn"},
    "pn":                   L3_PACKAGES - {"pn"},
    "moc":                  L3_PACKAGES - {"moc"},
    "cp":                   L3_PACKAGES - {"cp"},
    "mc":                   L3_PACKAGES - {"mc"},
    "diffusion":            L3_PACKAGES - {"diffusion"},
    "kinetics":             L3_PACKAGES - {"kinetics"},
    "fuel":                 L3_PACKAGES - {"fuel"},
    "thermal_hydraulics":   L3_PACKAGES - {"thermal_hydraulics"},
    "homogeneous":          L3_PACKAGES - {"homogeneous"},
    "derivations":          L2_PACKAGES | L3_PACKAGES,
}

WHITELIST: frozenset[tuple[str, str]] = frozenset({
    # RETIRE_IN_P3_FOLLOWUP — inline-import benchmark cross-check
    ("derivations/continuous/cases/diffusion.py", "diffusion"),
    # RETIRE_IN_P3_FOLLOWUP — MMS source uses MOCMesh; move to test side
    ("derivations/continuous/mms/moc.py", "moc"),
    # RETIRE_IN_P3_FOLLOWUP — sood_registry CPParams import
    ("derivations/continuous/sood_registry/builders.py", "cp"),
})


def _iter_python_modules(root: pathlib.Path) -> Iterable[pathlib.Path]:
    for p in root.rglob("*.py"):
        if "__pycache__" in p.parts:
            continue
        yield p


def _top_level_package(rel_path: pathlib.Path) -> str:
    return rel_path.parts[0]


def _imports_of(module_path: pathlib.Path) -> list[tuple[str, bool]]:
    """Parse imports, marking TYPE_CHECKING-guarded ones."""
    src = module_path.read_text()
    tree = ast.parse(src, filename=str(module_path))
    results: list[tuple[str, bool]] = []

    def _visit(node: ast.AST, in_tc: bool) -> None:
        if isinstance(node, ast.If) and isinstance(node.test, ast.Name) \
                and node.test.id == "TYPE_CHECKING":
            for child in node.body:
                _visit(child, True)
            for child in node.orelse:
                _visit(child, in_tc)
            return
        if isinstance(node, ast.ImportFrom) and node.module:
            results.append((node.module, in_tc))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                results.append((alias.name, in_tc))
        for child in ast.iter_child_nodes(node):
            _visit(child, in_tc)

    _visit(tree, False)
    return results


def _check_module(module_path: pathlib.Path) -> list[str]:
    rel = module_path.relative_to(ORPHEUS_ROOT)
    src_pkg = _top_level_package(rel)
    if src_pkg not in FORBIDDEN_EDGES:
        return []
    forbidden = FORBIDDEN_EDGES[src_pkg]
    rel_str = str(rel).replace("\\", "/")
    violations: list[str] = []
    for module_name, is_tc in _imports_of(module_path):
        if not module_name.startswith("orpheus."):
            continue
        tgt_pkg = module_name.split(".")[1]
        if tgt_pkg not in forbidden:
            continue
        # TYPE_CHECKING tolerance: L1/L2 importing L3 names for typing.
        if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES):
            continue
        if (rel_str, tgt_pkg) in WHITELIST:
            continue
        violations.append(
            f"{rel} imports {module_name} "
            f"(forbidden: {src_pkg} → {tgt_pkg})"
        )
    return violations


_ALL_MODULES = sorted(_iter_python_modules(ORPHEUS_ROOT))


@pytest.mark.foundation
@pytest.mark.parametrize("module_path", _ALL_MODULES, ids=str)
def test_no_forbidden_imports(module_path: pathlib.Path) -> None:
    violations = _check_module(module_path)
    assert not violations, "\n".join(violations)
```

(76 LoC inside the test file, under the 80-LoC budget. The `ast`
walker handles the TYPE_CHECKING tolerance natively; the WHITELIST
documents the migration to-do list.)

---

## SECTION 2 — Cross-cutting Phase 3 verification gates

For each step in the REVISED sequencing (P3.1 → P3.5 → P3.0 → P3.2 →
P3.3 → P3.4 → P3.6), the gates below MUST pass. Every step is
behaviour-preserving by design; only P3.4 + P3.6 introduce new
production paths.

### 2.1 P3.1 — Import-linter test (FIRST)

**Tests that MUST stay green at every commit**:

- The full pytest suite — `pytest -q` — confirms no test imports
  anything the linter rejects. The linter is parametrised per
  module so a failure isolates the offender by name.

**Tests at risk**:

- None. The linter is additive — it never alters production behaviour.

**New V&V test BEFORE this step**: this STEP IS the new V&V test.
Per `lessons-L7`-style "tests before optimization" — the linter
lands FIRST so every subsequent P3.x is gated by it. The expected
initial failure list (§1.6) is the migration to-do.

**Bit-identity judgment**: N/A — no production code changes.

### 2.2 P3.5 — `range` → `codomain` rename (BEFORE P3.2/P3.3)

**Tests that MUST stay green at every commit**:

- `tests/numerics/test_operator.py` — pins the `_AdjointOperator`
  machinery that reads `inner.range` at lines 511-516 of
  `operator.py`. The rename touches this read site directly.
- `tests/numerics/test_projection_operators.py` — `MomentProjection`
  currently exposes BOTH `.codomain` and `.range` (transitional dual);
  the rename retires `.range`.
- `tests/numerics/test_spherical_harmonic_space.py` — the new test
  at line 422 (`test_moment_projection_codomain_is_spherical_harmonic_space`)
  reads `.codomain` directly; no change needed.
- `tests/numerics/test_space.py` — `FunctionSpace` is the upstream
  type carrying the attribute; the rename propagates to it.
- `tests/sn/test_scattering_operator.py` — `ScatteringOperator` reads
  operator `.range` / `.codomain` in its adjoint machinery; pin.
- All Krylov / SI tests that exercise `_AdjointOperator.apply`
  indirectly (i.e. `tests/numerics/test_iteration_angular_flux.py`).

**Tests at risk**:

- ANY test fixture that reads `op.range` directly. Grep target list
  for the rename audit (per `lessons-L20`): every occurrence of
  `\.range\b` in `tests/` and `orpheus/` that resolves to an
  operator attribute (NOT a Python builtin call `range(...)`). The
  audit MUST distinguish the two — recommend Nexus `rename` if
  available (graph-aware) or `grep -nE '\.range\b'` followed by
  per-hit triage.
- `@pytest.mark.verifies('...')` decorators are NOT affected (they
  carry equation labels, not attribute names).
- Test fixtures that construct mock operators with a hand-coded
  `.range` attribute MUST update to `.codomain` (e.g. dataclass
  test doubles in `tests/numerics/test_operator.py`).

**New V&V test BEFORE this step**: NO. The existing tests already
pin the contract. The rename is mechanical — exhaustive grep + the
new linter catches any miss.

**Bit-identity judgment**: STRICT bit-identical everywhere. A pure
attribute rename produces no numerical change; the FP reduction tree
is unchanged. Per `vv-principles` §"Bit-identity vs principled-
equivalence", criterion 1 (principled) is N/A because the
implementation is not changed — only the API surface. ANY snapshot
drift post-P3.5 is a bug, not a refactor artefact.

### 2.3 P3.0 — Documentation of the criterion + layer table

**Tests that MUST stay green at every commit**:

- The Sphinx build — `sphinx-build -W docs docs/_build/html` — MUST
  build clean. P3.0 adds or extends `docs/architecture/layering.rst`.
- `python -m tests._harness.audit` — V&V matrix unaffected.

**Tests at risk**: None.

**Bit-identity judgment**: N/A — documentation only.

### 2.4 P3.2 — `numerics/` internal reorganization

**Tests that MUST stay green at every commit**:

- `tests/numerics/` ENTIRE suite — STRICT bit-identical. The
  reorganisation moves modules; imports follow. Any test that
  imports `from orpheus.numerics.space import FunctionSpace` MUST
  continue working (either via direct rewire or a back-compat
  re-export shim in `numerics/__init__.py`).
- `tests/numerics/test_spherical_harmonic_space.py` — moved to
  `numerics/spaces/`; the test file path is unchanged (it lives in
  `tests/`); the production-side import is now
  `from orpheus.numerics.spaces.spherical_harmonic_space import ...`.
- `tests/sn/regression/test_dd_regression.py` — STRICT bit-identical.
- All L1 MMS gates (`tests/sn/test_mms_aniso.py`, curvilinear analog).

**Tests at risk**:

- Tests that import `from orpheus.numerics.space import FunctionSpace`
  (i.e. anywhere the old path is used). The P3.2 commit MUST either
  (a) leave a back-compat shim in `numerics/space.py` that re-exports
  the new location, retired in a subsequent commit; OR (b) update
  every test import in lockstep. Option (a) is the lessons-L20
  preferred ordering (retirement = test migration, but DEPRECATION
  shim allowed for one merge cycle per `feedback_aggressive_retirement`).

**New V&V test BEFORE this step**: NO. The reorganisation is
mechanical. The import-linter is the new V&V test (P3.1 already
landed).

**Bit-identity judgment**: STRICT bit-identical. Module relocation
does not touch FP reduction trees. The plan's "spherical_harmonics.py
shim deletion + four sn-quadrature delegator rewire" is the ONE risk
point: the delegators MUST consume the same `evaluate_real_sh`
function (now via `SphericalHarmonicBasis.evaluate`), bit-identical.
Verify by re-running `tests/numerics/test_spherical_harmonics.py`
post-rewire — it carries L0/L1 tests against fixed expected values
(the SH evaluator itself is unchanged).

### 2.5 P3.3 — Introduce `transport/`; migrate fields + sources

**Tests that MUST stay green at every commit**:

- `tests/sn/test_typed_fields.py` — pins `ScalarFlux`, `AngularFlux`,
  `HarmonicMomentField` dunder algebra. The L2 base must satisfy
  these tests after migration; SN-specific behaviour (e.g.
  `from_flat_with_traces`) tested on the L3 adapter.
- `tests/sn/test_angular_flux_with_boundary.py` — currently pins
  `AngularFlux.from_flat_with_traces` (lines 248, 258, 274, 281, 290,
  304). After CC.4 resolution (see §3 below), this test moves to
  `tests/sn/test_angular_flux_b1pp_adapter.py` OR keeps its current
  path with the test now exercising the L3 adapter. STRICT
  bit-identical.
- `tests/sn/test_harmonic_moment_field.py` — pins
  `HarmonicMomentField` algebra. After P3.3 the class lives in
  `transport/fields/`; tests follow.
- `tests/numerics/test_iteration_angular_flux.py` (line 35) imports
  `AngularFlux` and exercises the cross-package `_ravellable`
  Protocol. After the split, this test consumes the L2 base
  `AngularFlux` (the duck-typed Protocol the Krylov inner loop reads
  has no eq-map dependency). STRICT bit-identical.
- `tests/sn/test_scattering_operator.py` — STRICT bit-identical.
- `tests/sn/regression/test_dd_regression.py` — STRICT bit-identical.
- All L1 MMS gates.

**Tests at risk**:

- Any test importing the migrated symbols (`HarmonicMomentField`,
  `ScalarFlux`, `AngularFlux`, `IsotropicSource`, `PerOrdinateSource`)
  must update the import path. Recommend a one-commit back-compat
  re-export shim in `sn/__init__.py` deleted in the next commit per
  the deprecation-shim rule.

**New V&V test BEFORE this step**:

- **YES** — write a new test `tests/transport/fields/test_angular_flux_base_algebra.py`
  that exercises ONLY the L2 algebra (storage, dunders,
  `integrate_angular`) of the new `transport/fields/angular_flux.py`
  base class. This is the pin that says "the L2 base is functionally
  complete without the SN adapter." Per `lessons-L7`-style "tests
  before optimization" — write the L2-only test BEFORE moving the
  class, so the test is the contract.
- **YES** — write a new test `tests/transport/fields/test_harmonic_moment_field_base.py`
  for the same reason on `HarmonicMomentField`. The plan-§P3.3
  forced dependency cleanup introduces a `SpatialGroupMesh` Protocol;
  the test pins the Protocol contract (the L2 base reads only
  `mesh.ng`, `mesh.nx`, `mesh.ny`).

**Bit-identity judgment**: STRICT bit-identical. The migration is a
package rename + a class split; numerical behaviour is unchanged.
The single risk: `AngularFlux.__post_init__` currently does
`from .boundary_flux import BoundaryFlux` (line 135) — a sibling
import inside SN. After the split, the L2 base may need a different
boundary-allocation strategy (see §3 below).

### 2.6 P3.4 — Problem / Solver split (greenfield per CC.7)

**Tests that MUST stay green at every commit**:

- `tests/numerics/test_iteration.py` — current `SourceIteration` +
  `KEigenvalue` tests. After P3.4 the classes live in
  `numerics/solvers/`; tests follow path.
- `tests/numerics/test_iteration_angular_flux.py` — STRICT
  bit-identical; the Krylov path is unchanged.
- `tests/sn/test_scattering_operator.py` — STRICT bit-identical.
- `tests/sn/regression/test_dd_regression.py` — STRICT bit-identical
  AT THE NUMERICAL LEVEL; per §A.3 of the Phase 1 plan, criterion 3
  (FP-non-associativity bound) applies. The new
  `CriticalityProblem + PowerIteration` wiring goes through
  `KEigenvalue` (renamed `PowerIteration`); the call site is
  reconstructed but the FP reduction tree should be identical
  (renames only).
- L1 MMS gates — STRICT bit-identical.
- `tests/homogeneous/` (if it exists) — `homogeneous/solver.py:26`
  imports `power_iteration`; CC.8 retires this and rewires to
  `KEigenvalue`. Bit-identity AT FP-non-assoc bound (per
  `vv-principles`).

**Tests at risk**:

- Tests pinning `power_iteration` directly (the legacy function).
  Per CC.8, these MIGRATE to the `CriticalityProblem +
  PowerIteration` API in the same commit that retires
  `power_iteration` (per `feedback_retirement_means_test_migration`).
- Tests pinning `KEigenvalue` by name — the rename to
  `PowerIteration` (class) is the bigger surface. Recommend back-compat
  re-export `KEigenvalue = PowerIteration` retired in the next commit.

**New V&V test BEFORE this step**:

- **YES** — write `tests/transport/problems/test_criticality_problem.py`
  that pins the declarative `(A_loss, F)` triple → `PowerIteration`
  → keff path against a structurally-independent analytical
  reference (homogeneous infinite medium k_inf = νΣ_f/Σ_a, multi-
  group, ≥2G per cardinal rule). This is a NEW L1 verification
  test — it pins the new wiring before the SN solver rewires.
- **YES** — write `tests/transport/problems/test_fixed_source_problem.py`
  pinning the `(L, S, F, q) + SourceIteration` path against a
  similar reference. Per `vv-principles` §"L1 without L0 = compensating
  errors" — these tests pin the new construction site against an
  analytical limit, not against the OLD wiring.

**Bit-identity judgment**: **PRINCIPLED-EQUIVALENCE acceptable** per
`vv-principles` §"Bit-identity vs principled-equivalence". The new
`CriticalityProblem + PowerIteration` composition routes through the
same primitives as the old `KEigenvalue.solve` — bit-identical
expected on slabs (well-conditioned). Curvilinear cases may drift
within the existing `rtol=5e-6` iteration floor. The three criteria
apply:
1. **Principled**: each new intermediate (`Problem`, `Solver`,
   `PowerIteration.solve`) is a named, inspectable object. SATISFIED.
2. **Structurally-independent reference**: the new
   `tests/transport/problems/test_criticality_problem.py` pins keff
   against k_inf for homogeneous reflective — analytical limit. The
   curvilinear / heterogeneous snapshots inherit the existing L1 MMS
   gates as references. SATISFIED.
3. **Drift dimensionally explainable**: rename-only reconstruction
   of the call site; reduction depth unchanged; drift ≤ `outer_iters
   × ULP`. Per snapshot bound (`rtol=1e-12` slab, `5e-6` curvilinear).

If drift exceeds bound: the rewire is wrong (NOT a tolerance issue).
DO NOT loosen snapshot tolerances.

### 2.7 P3.6 — `kinetics/` restructure

**Tests that MUST stay green at every commit**:

- Any existing `tests/kinetics/` suite — STRICT bit-identical. The
  point-kinetics solver IS PHYSICS; numerical drift must respect the
  same bit-identity-or-FP-bound rule.
- L1 / L0 / L2 gates above — unaffected (kinetics doesn't intersect
  the transport pipeline directly).
- `tests/transport/problems/test_initial_value_problem.py` — NEW
  (see below).

**Tests at risk**:

- `kinetics/solver.py` will dissolve; tests pinning the old API
  MIGRATE to (a) the new `transport/problems/initial_value.py`
  + (b) `numerics/solvers/time_stepping.py` per
  `feedback_retirement_means_test_migration`. The migration is a
  per-test rewire — recommend the same one-commit-shim pattern.

**New V&V test BEFORE this step**:

- **YES** — write `tests/numerics/solvers/test_time_stepping.py`
  pinning each time-stepper primitive (BDF1, BDF2, Crank-Nicolson)
  against analytical references (decay-equation closed form for
  BDF1; oscillator for higher orders). This is L1 verification;
  pillar: closed-form (Branch 1A).
- **YES** — write `tests/transport/problems/test_initial_value_problem.py`
  pinning the declarative `InitialValueProblem + TimeStepper`
  composition against the same analytical references.

**Bit-identity judgment**: **PRINCIPLED-EQUIVALENCE** required. The
restructure changes the call wiring of every kinetics consumer; the
new path's FP reduction may differ. Per `vv-principles`:
1. **Principled**: `Problem + Solver` decomposition exposes named
   intermediates (the time-step result IS a flux + precursor field,
   not an unnamed array). SATISFIED.
2. **Structurally-independent reference**: NEW L1 tests above —
   decay-equation / oscillator analytical limits. SATISFIED if the
   tests are written FIRST.
3. **Drift**: kinetics solvers iterate; drift ≤ `iter_count × ULP`.
   Per-snapshot bound (if a kinetics snapshot suite exists).

The risk is that NO L1 kinetics references exist today and the
existing tests are mostly L4 (code-to-code self-checks). The
new V&V tests above are MANDATORY before the restructure.

---

## SECTION 3 — AngularFlux design-note verification (CC.4 resolution)

### 3.1 The split

Per plan §P3.3 "AngularFlux design note":

- **L2 base** (`transport/fields/angular_flux.py`) — pure algebra:
  storage (`values`, `mesh`, `boundary`, `history_depth`, `_history`),
  dunders (`__add__`, `__mul__`, `__rmul__`, `__neg__`, `<<`),
  reductions (`integrate_angular`), iteration history
  (`stash`, `__call__`), shape invariants (`__post_init__`).
- **L3 SN adapter** (`sn/angular_flux_b1pp_adapter.py`) — exposes
  the `from_flat_with_traces`, `to_flat_with_traces` factory +
  consumer, which import from `sn/operator.py`
  (`build_equation_map_with_traces`,
  `solution_to_angular_flux_with_traces`, `pack_with_traces`). These
  are SN B1''-eq-map primitives.

### 3.2 Tests that pin the L2 algebra

Identified by reading the worktree:

- `tests/sn/test_typed_fields.py::TestAngularFlux` (lines 153-205).
  Pins shape invariant, arithmetic, `integrate_angular`,
  `at_ordinate` — all pure L2 algebra. After the split, these tests
  EITHER (a) move to `tests/transport/fields/test_angular_flux.py`
  to consume the L2 base directly, OR (b) stay in `tests/sn/` but
  import from `orpheus.transport.fields.angular_flux` (the cross-
  package import is now legal — L3 tests using L2 base is fine).
- `tests/numerics/test_iteration_angular_flux.py::TestRavellableProtocol`
  (lines 63-128). Pins the duck-typed `_ravellable` Protocol
  (`__ravel__`, `_unravel_like`, `zeros_like`, `_l2_norm`). Per CC.7
  the Krylov / SI primitives consume this Protocol — the L2 base
  MUST implement it. After the split, this test imports the L2
  base; STRICT bit-identical.

### 3.3 Tests that pin the L3 B1''-eq-map machinery

- `tests/sn/test_angular_flux_with_boundary.py` lines 243-310 — pins
  `from_flat_with_traces` round-trip (lines 248, 258, 274), face
  decoder (line 281), error-message contract (line 304). All
  exercise the SN B1''-eq-map.
- The same file probably tests `to_flat_with_traces` round-trip at
  line 285-295 (visible in grep — `flat_out = psi.to_flat_with_traces()`).

These tests move to `tests/sn/test_angular_flux_b1pp_adapter.py` (a
new file), or stay in their current location with imports updated.
Recommend the SAME file path renamed only if the symbol set changes
(it does: the import is now `from orpheus.sn.angular_flux_b1pp_adapter
import AngularFlux` OR the adapter exposes a `classmethod` on the
L2 base via post-import injection).

### 3.4 Does the split break any test?

- **`AngularFlux.__post_init__` boundary auto-allocation**
  (`angular_flux.py:131-138`) currently does
  `from .boundary_flux import BoundaryFlux` — a sibling import
  inside SN. After the split, the L2 base **needs the same boundary
  allocation but BoundaryFlux is L3** (it depends on `SNMesh`-shaped
  faces). **Two options**:
  - **Option A (preferred)**: Move `BoundaryFlux` to
    `transport/fields/boundary_flux.py` alongside `AngularFlux`.
    Verify `BoundaryFlux` reads only `mesh.{ng, nx, ny}` (the same
    `SpatialGroupMesh` Protocol). If yes, this is a clean L2 promote.
  - **Option B**: Defer boundary allocation to the adapter. The L2
    base `__post_init__` accepts `boundary=None` and stores it; the
    L3 adapter's factories ensure a BoundaryFlux is always supplied
    at construction. This BREAKS the `coding-elegance` Pattern 4
    invariant ("illegal states unrepresentable — every AngularFlux
    has a non-None boundary by the time anything reads it") if any
    L2 consumer instantiates AngularFlux without going through the
    adapter.

  **Recommendation**: Option A. Verify by Read on
  `sn/boundary_flux.py` (lines for the L3 vs L2 imports) AS PART OF
  P3.3 PRE-MIGRATION RECON. If `BoundaryFlux` imports only
  `geometry/` / `numerics/` / `transport/`, it promotes cleanly.

- **`AngularFlux.from_flat_with_traces`** depends on
  `mesh.curvature` (line 384) — this attribute is on `SNMesh`.
  After the split, the L3 adapter still owns this branch (the L2
  base never sees `mesh.curvature`). No break.

- **`AngularFlux` mesh field type annotation** — currently
  `mesh: "SNMesh"` (string annotation). After the split, the L2
  base uses `mesh: "SpatialGroupMesh"` (the new Protocol). The
  string annotation is non-load-bearing at runtime; structural
  typing works through the Protocol. NO break.

### 3.5 Required pre-P3.3 reconnaissance

Dispatch the **explorer** sub-agent (proactive trigger per
`subagent-handoff-protocol` — "operator-algebra carve crossing
subsystem boundaries"; this is a typed-field carve crossing
`SN/transport` boundary) to:

1. Audit `sn/boundary_flux.py` — does it import anything other than
   `numerics/`, `geometry/`, `transport/` (post-AngularFlux-move)?
2. Audit every consumer of `AngularFlux.from_flat_with_traces` and
   `to_flat_with_traces`. Confirm they are all in `sn/` or
   `tests/sn/`.
3. Audit every consumer of the L2-only algebra (`integrate_angular`,
   `<<`, `__call__`, dunders) — confirm whether any cross-package
   tests rely on the L3 method surface.

The audit precedes P3.3 implementation. If finding #1 shows that
`BoundaryFlux` carries SN-specific structure (e.g. ordinate
permutations, reflective-BC face pairing tables), Option B becomes
mandatory and the L2 base accepts `boundary: Optional[...]`.

---

## SECTION 4 — Risk register

Phase 3 steps ranked by failure impact + likelihood. The risk is
expressed as "what could go wrong AND what fingerprint catches it."

| Rank | Step | Specific failure mode | Detection / catching test |
|---|---|---|---|
| **1 (highest)** | **P3.4 — Problem/Solver split + `power_iteration` retirement** | The `homogeneous/solver.py:26` import of `power_iteration` REWIRES to `KEigenvalue` (new name `PowerIteration`); the wiring of the outer eigenvalue loop changes call-site reconstruction (per memo §C2: `sn/solver.py:499, 539` and `:618, 659, 1452, 1484` all become Problem-construction sites). Risk: a convention drift (Failure Mode #6) in the new declarative `(L, S, F)` triple vs the imperative legacy wiring — e.g. fission operator F at `homogeneous` is constructed differently because the old `power_iteration` accepted a raw callable, the new `CriticalityProblem` accepts a typed `LinearOperator`. The bug fingerprint: 1G works (k = νΣ_f/Σ_a degeneracy), 2G+heterogeneous diverges or shifts by O(1). | (a) New L1 verification at `tests/transport/problems/test_criticality_problem.py` pinning k against k_inf for ≥2G homogeneous reflective; (b) `tests/sn/regression/test_dd_regression.py` slab-2G snapshots STRICT bit-identical (rtol=1e-12); (c) `tests/homogeneous/` keff against analytical limit. |
| **2** | **P3.6 — `kinetics/` restructure** | Time-stepping primitives (BDF1, BDF2, Crank-Nicolson) extracted and called via new `InitialValueProblem + TimeStepper` composition. Risk: any of (sign flip on the time-derivative; missing precursor source coupling; wrong staggering of fast/thermal flux) — all Modes #1, #3, #4. Existing kinetics tests are likely L4 (code-to-code), so an L0/L1 reference must be NEW. | NEW L1 verification at `tests/numerics/solvers/test_time_stepping.py` (decay equation, oscillator) and `tests/transport/problems/test_initial_value_problem.py`. The new tests pin against analytical decay / oscillator solutions BEFORE the restructure. Existing kinetics tests rewire and remain green. |
| **3** | **P3.3 — `AngularFlux` split + BoundaryFlux co-promote** | If Option B is chosen (boundary deferred to adapter), an L2 consumer that instantiates `AngularFlux` without going through the SN adapter gets a half-constructed object — `boundary=None` at L2 level, no auto-allocate. Risk: silent NaN propagation through subsequent reductions; Mode #2/#3 (variable swap / missing factor). | (a) NEW L2 test `tests/transport/fields/test_angular_flux_base_algebra.py` with explicit `boundary=...` and explicit `boundary=None` assertion — make the contract explicit. (b) `tests/numerics/test_iteration_angular_flux.py::test_zeros_like_returns_angular_flux` (line 115) MUST verify `zero_psi.boundary` is non-None — this is the existing L0 contract that prevents the silent-NaN failure mode. (c) Option A (promote BoundaryFlux) sidesteps the entire risk and is the recommended path. |
| **4** | **P3.5 — `range` → `codomain` rename** | A reference to `op.range` survives the rename and silently shadows Python's builtin `range(...)` — likely produces `TypeError` at runtime, not a numerical bug. But: a fixture that constructs a mock operator with `.range` and the production reads `.codomain` produces an `AttributeError` at use site, possibly only under a rarely-exercised path. | (a) The P3.1 linter does NOT catch this (it's an attribute, not an import). (b) Comprehensive `grep -nE '\.range\b'` audit (per `lessons-L20`) gating the rename commit. (c) Full pytest suite run AFTER the rename — the `AttributeError` surface in any test that exercises the operator algebra. (d) Static analysis (mypy / pyright in strict mode) flags the rename gaps if type stubs are kept current. |
| **5** | **P3.2 — `numerics/` reorganisation** | The `spherical_harmonics.py` shim retirement + four sn-quadrature delegators rewire. Risk: a delegator imports from the shim and breaks when the shim deletes; OR the new `SphericalHarmonicBasis.evaluate` API has a different signature than `evaluate_real_sh` and the rewire silently does the wrong thing. | (a) `tests/numerics/test_spherical_harmonics.py` STRICT bit-identical (it tests the SH evaluator against fixed expected values). (b) `tests/sn/test_quadrature.py::TestProductQuadrature` — pins the per-quadrature delegators (`tests/numerics/test_quadrature_directional.py:75, 393` IS a likely consumer; verify at P3.2 time). (c) Run `pytest -q tests/numerics/test_spherical_harmonics.py tests/sn/` after the shim deletion. |
| **6** | **P3.1 — Import-linter** | Risk: a TYPE_CHECKING tolerance bug allows a real import edge to slip through; OR the WHITELIST accumulates stale entries that mask real violations. | The linter parametrises per-module so failures are isolated. Coverage gap can be audited via `pytest -q tests/test_layer_imports.py --collect-only` to confirm every `orpheus/**/*.py` is in the test set. Stale whitelist entries: each carries a `RETIRE_IN_PX.X` comment; a follow-up audit (out of Phase 3 scope) verifies whether each is still load-bearing. |
| **7 (lowest)** | **P3.0 — Documentation** | Risk: stale docs after package moves. Mitigation: the Sphinx build catches dead `:mod:` references via `-W` warnings. | `sphinx-build -W docs docs/_build/html` is the gate. |

### 4.1 Concurrence with the REVISED sequencing

**I CONCUR with the REVISED sequencing** in plan §"Sequencing within
Phase 3 (REVISED post-Phase-1 per QA learnings)":

> P3.1 → P3.5 → P3.0 → P3.2 → P3.3 → P3.4 → P3.6.

Three reasons:

1. **P3.1 first** is structurally mandatory — the linter is the
   safety net for every move. Without it, every P3.x commit could
   introduce an unnoticed forbidden edge.
2. **P3.5 second** is the load-bearing INSIGHT from Phase 1's
   transitional dual-name (`MomentProjection.range` + `.codomain`).
   Doing the rename EARLY means every subsequent move lands against
   the canonical attribute name. Doing it LATE (as the original
   pre-Phase-1 sequencing had it) means P3.2 + P3.3 + P3.4 all
   carry the dual-name overhead and must be re-rebased over the
   rename. QA's recommendation is correct.
3. **P3.6 last** is correct — kinetics depends on the new
   `numerics/solvers/time_stepping.py` (extracted in P3.6 itself
   from kinetics/) AND on `transport/problems/initial_value.py`
   (created in P3.4 conceptually but landed here). It cannot precede
   P3.4 + P3.2 + P3.3 because its target homes are created by those
   steps.

**Minor concurrence note**: the proactive trigger per
`subagent-handoff-protocol` ("operator-algebra carve crossing
subsystem boundaries") fires at P3.3 (AngularFlux split crossing
SN ↔ transport boundary) AND P3.4 (problem/solver split crossing
sn/solver.py call-site convention boundary). Both should dispatch
the explorer sub-agent first for the consumer-inventory pass before
implementation begins. The P3.3 inventory specifically must answer
the `BoundaryFlux` Option-A-vs-B question identified in §3.4.

---

## Pointers

- Plan: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/.claude/plans/moment_space_and_layering_plan.md`
- Phase 1 verification plan: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/.claude/agent-memory/test-architect/phase1_verification_plan.md`
- Phase 1 QA review: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/.claude/agent-memory/qa/phase1_moment_space_review.md`
- Lessons load-bearing for Phase 3: L11 (structural-independence-via-elimination-of-FP), L17 (convention crosswalk before carve), L18 (Pattern 7 at the producer), L20 (retirement requires dependency audit)
- Skills: `vv-principles`, `coding-elegance` Pattern 7, `algebra-of-record` (for the L0/derivations classification), `subagent-handoff-protocol` (for proactive explorer dispatches at P3.3 / P3.4)
- The P3.1 test file path: `tests/test_layer_imports.py`
- The new V&V test paths (created BEFORE each step that needs them):
  - P3.3: `tests/transport/fields/test_angular_flux_base_algebra.py`, `tests/transport/fields/test_harmonic_moment_field_base.py`
  - P3.4: `tests/transport/problems/test_criticality_problem.py`, `tests/transport/problems/test_fixed_source_problem.py`
  - P3.6: `tests/numerics/solvers/test_time_stepping.py`, `tests/transport/problems/test_initial_value_problem.py`

## Self-improvement entries

### New failure mode → skill update (NONE for this plan)

Phase 3 covers existing failure modes (modes 1, 2, 3, 6 in
`vv-principles` §"6 AI failure modes" plus the cross-cutting
layering-violation failure mode that is structurally caught by the
import-linter, not by a runtime check). No new failure mode
introduced.

### Plan-rejection counter-examples (none yet)

If qa or the user rejects any row in this plan, log a one-paragraph
counter-example under
`.claude/agent-memory/test-architect/feedback_*.md`.

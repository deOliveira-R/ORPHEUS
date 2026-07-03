---
name: cross-layer-relocation-carve-verification
description: Verifying a "relocate a module DOWN a layer + route a hardcoded dispatch through a registry" carve. Two structural hazards the naive plan misses (registration timing masked by process-global state; the layer-imports gate's TYPE_CHECKING exemption not covering INPUT packages) + the fresh-process gate that catches the first.
metadata:
  type: project
---

When the carve is **relocate a method-agnostic unit from a method
package (L3 `sn/`) DOWN to a shared lower layer (`geometry/`) AND
replace its hardcoded leaf-dispatch import with a runtime registry
lookup**, the value/behavior gate is the easy part. The two
load-bearing hazards are STRUCTURAL and the naive plan misses both.
First worked: `realize_recursively` `sn/boundary_realize.py →
geometry/boundary/` (branch `refactor/operator-inverse-algebra`,
2026-06-26; spec `.claude/plans/realize_recursively_move_spec.md`).

## Hazard 1 — registry-routing introduces a registration-TIMING
regression that process-global class state HIDES in-suite.

The pre-move walker hard-imports the concrete realizer
(`from orpheus.sn.boundary_realizer import SNBoundaryRealizer`) inside
its body — that import GUARANTEES the `@Registry.register("SN")`
decorator has fired by leaf-dispatch time. After the move, routing
through `Registry.get(method)` does NOT import the realizer, and
importing the NEW lower-layer home pulls in ZERO of the method
package (it must — that is the whole point). So in a **fresh process**
`get("SN")` raises `RegistryError(Available: [])`. The registry is
**process-global class state**, so any earlier test that imported the
method package MASKS the miss — a green suite is NOT evidence. PROVE
the hazard in a fresh interpreter (`subprocess.run([sys.executable,
"-O", "-c", script])`); a same-process assertion is **vacuous for the
timing claim** (Mode-11-adjacent: the assert fires but never
exercises the empty-registry condition). The resolution contract:
"the walker dispatches `get(method)`; importing the method's realizer
module is the CALLER's responsibility (the per-method subpackage
self-registers on import)." Production callers live in the method
package and do import it → the fresh-process gate proves they satisfy
it. Pin the miss loudly with a negative gate (subprocess
`del Registry._registry["SN"]` → the walker call raises `RegistryError`,
NOT an `AttributeError` deep in the sweep). Mutation-teeth: the same
`del` in-subprocess reddens the positive fresh-process gate.

## Hazard 2 — the VERBATIM move trips the layer-imports gate it is
meant to SATISFY, via a TYPE_CHECKING import.

`tests/test_layer_imports.py` parses runtime AST imports AND records
`TYPE_CHECKING` `ImportFrom` (it visits the TC block). Its TC-tolerance
(`:148`) exempts ONLY source packages in `L1|L2` (`numerics`/
`transport`). `geometry`/`data` are **INPUT** packages — **NOT
exempted**. So a relocated file's `if TYPE_CHECKING: from orpheus.sn.
method_space import SNMethodSpace` (which was legal in `sn/`) becomes a
**flagged `geometry→sn` forbidden edge**. The gate going RED on the new
file is the CORRECT early signal that the typing import was not removed.
Resolution: drop ALL `orpheus.sn…` imports from the new home including
the TC one; retype the method-space parameter as `Any` (matching the
registry Protocol's `realize(self, law, method_space: Any)`). This is
also the *more correct* coupling — a now-method-agnostic unit should
not be statically bound to one method's space type. Verify by replaying
the gate's AST visitor on the current file BEFORE the move, so the
implementer knows the TC import is a blocker (not a surprise at gate
time).

## The "layer inversion" framing is usually DOCUMENTATION-ONLY at
runtime — refute it at design time (lessons L10).

Before crediting "removes a geometry→sn layer inversion", PROVE the
runtime edge exists: `python -c "import sys, orpheus.geometry;
print(sorted(m for m in sys.modules if m.startswith('orpheus.sn')))"`.
For `realize_recursively` it was **empty** — every `from orpheus.sn…`
in `geometry/` was a docstring `:class:`/`:func:` cross-ref or a lazy
body import in a method NO production caller invokes. The brief's cited
line numbers were all inside docstring `.. code-block::` examples. The
ONLY runtime importer of the moved unit was a TEST
(`test_law_composition.py:50`). So the real payoff is (a) remove the
one TEST runtime edge, (b) co-locate the unit with its Protocol +
registry, (c) drop docstring/TC coupling — NOT remove a production
cycle. State the payoff accurately in the commit; a "fix the runtime
inversion" claim the `sys.modules` diff cannot corroborate is an
over-claim a fresh session inherits.

## The behavior gate for a pure relocation

Bit-identity is genuinely expected (operator-assembly code is
byte-relocated). The existing composed-path wall
(`test_law_composition.py`: `OperatorSum`-of-`ScaledOperator`
structure + an L1 pointwise-numpy `0.3·a+0.7·b` reference that does
NOT call the walker = structurally independent) is sufficient — AUDIT
that it exercises the COMPOSED path (not just leaves) and rewire its
import to the successor (behavioral test → rewire, NOT delete).
Optional stronger anchor: a same-process `old≡new` `assert_array_equal`
across both import paths (unavailable on a hard no-shim move — then the
`nulp` pointwise reference + structural assertions ARE the evidence).

## Migration-table reminders specific to this carve shape

Three-search retirement audit (graph callers + grep code/tests/docs +
direct constructors): the moved module's docstring `:func:`/`:mod:`
cross-refs across the SIBLING package files (and any RUNTIME ERROR
STRING that names the old module path — user-facing) must be
grep-updated; Sphinx `-W` catches a ref that no longer RESOLVES but
NOT one that resolves by coincidence, so the explicit grep is required
alongside `-W`. Watch for a STALE test-path pointer in a docstring
(here `tests/sn/test_boundary_realize.py` no longer exists; the tests
are in `tests/geometry/test_law_composition.py`). When a concurrent
cleanup touches the SAME `__init__.py` seam, flag the one test file
that imports across BOTH tendrils and run it only after both land.

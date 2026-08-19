# Verification spec — relocate `realize_recursively` `sn/ → geometry/boundary/`

**Author:** test-architect (proactive pre-carve pass)
**Branch:** `refactor/operator-inverse-algebra`
**Env:** Host `.venv/bin/python` (3.14); canonical `python -O -m pytest`, SERIAL
(`-p no:xdist --timeout=300 -p no:cacheprovider`).
**Status:** READ-ONLY investigation complete. No production edits made.
**Mode:** spec only — the main agent owns implementation + the concurrent `SNMethodSpace`
re-export repoint (NOT covered here).

---

## TL;DR — what changed about the plan's premises (refute-at-design-time, AGENT.md L10)

Two premises in the brief are **measurably false at HEAD**, and one **new hazard** the
brief did not name dominates the whole carve. The plan below is built on the verified
state, not the brief's stated state.

1. **REFUTED — "geometry→sn layer inversion (runtime)" does not exist.**
   `import orpheus.geometry` / `import orpheus.geometry.boundary` pulls in **ZERO**
   `orpheus.sn` modules (proven: `sys.modules` diff is empty). Every `from orpheus.sn…`
   occurrence in `orpheus/geometry/` is either a **docstring `:class:`/`:func:` cross-ref**
   or a **lazy import inside a method body that no production caller invokes**
   (`BoundaryTraceLaw.realize` at `_base.py:282` is documented "no current production
   caller routes through this hook"). The brief's cited lines
   (`_base.py:129`, `__init__.py:167/180`, `__init__.py:168 imports SNMethodSpace`) are
   **all inside docstring `.. code-block:: python` examples** — NOT executable imports.
   The layer-imports gate (`tests/test_layer_imports.py`) parses **runtime AST imports**
   and therefore reports **no `geometry→sn` violation today**. So the move's headline
   payoff is a *documentation/optics* cleanup (and a real removal of the one **test**-side
   runtime edge `tests/geometry/test_law_composition.py:50`), not the removal of a
   production runtime cycle. Ship it for the right reason; do not claim a runtime-inversion
   fix the graph cannot corroborate.

2. **NEW HAZARD (dominant) — registry-routing introduces a registration-TIMING
   regression that is INVISIBLE under normal test runs.** Today the walker hard-imports
   `SNBoundaryRealizer` inside its body, which **guarantees the `@register("SN")`
   decorator has fired** at the leaf-dispatch moment. After the move, routing through
   `BoundaryRealizerRegistry.get("SN")` **does NOT import the SN realizer** — and
   `import orpheus.geometry.boundary` leaves the registry **empty** (`method_names() == ()`,
   proven in a fresh process). A fresh-process `get("SN")` raises
   `BoundaryRealizerRegistryError("…Available: [].")`. Because the registry is
   **process-global class state**, any test that imported `orpheus.sn` earlier in the
   session **masks the miss** — so a green suite is NOT evidence the routing is safe. The
   §2 gate below MUST run in a **subprocess / fresh interpreter** to be valid.

3. **NEW HAZARD (silent gate break) — the verbatim move TRIPS the layer-imports gate it
   is meant to satisfy.** `boundary_realize.py:119` carries a `TYPE_CHECKING` import
   `from orpheus.sn.method_space import SNMethodSpace`. The gate's `TYPE_CHECKING`
   tolerance (`test_layer_imports.py:148`) exempts **only** source packages in `L1|L2`
   (`numerics`/`transport`). `geometry` is an **INPUT** package — **not exempted** — so the
   relocated file's `TYPE_CHECKING` sn import would be **flagged as a `geometry→sn`
   forbidden edge** (proven by replaying the gate's AST visitor on the current file). The
   move MUST drop ALL `orpheus.sn…` imports from the new home, including the
   `TYPE_CHECKING` one (see §3.B / §4 — retype `method_space` as `Any`, matching the
   `BoundaryRealizer` Protocol the walker now dispatches through; this is also the *more
   correct* coupling for a now-method-agnostic walker).

---

## Verified pre-state (captured this session, for the implementer to diff against)

| Gate | Pre-state |
|------|-----------|
| `tests/geometry/test_law_composition.py` (the wall) | **18 passed** (`-O` serial) |
| `tests/geometry/test_boundary.py` + `test_bc_equivalence_snapshot.py` | **32 passed** |
| `tests/test_layer_imports.py` | **291 passed** (no geometry→sn edge today) |
| `import orpheus.geometry.boundary` → sn modules pulled | **NONE** |
| `BoundaryRealizerRegistry.method_names()` after `import orpheus.geometry.boundary` only | **`()`** (empty) |
| `…method_names()` after `import orpheus.sn.boundary_realizer` | **`('SN',)`**, `get("SN") → SNBoundaryRealizer` |
| fresh-process `get("SN")` without importing sn | **raises `BoundaryRealizerRegistryError`** |
| The ONLY runtime `realize_recursively` importer | `tests/geometry/test_law_composition.py:50` (a TEST; no production importer) |

Baselines the move must not regress: **the 7-and-only-7 pre-existing reds** (full serial
suite), **CLI pyright `orpheus/` = 412**, Sphinx **`-W`** exit 0.

---

## The walker's real dependency surface (what the new home is allowed to import)

| Import in `boundary_realize.py` | Layer edge | Legal in `geometry/boundary/`? |
|---|---|---|
| `from orpheus.geometry.boundary import BoundaryTraceLaw, LawScaled, LawSum` | intra-geometry | YES (becomes a relative `from . import …` / `from ._composition`/`._base`) |
| `from orpheus.numerics.operator import LinearOperator, OperatorSum, ScaledOperator` | geometry → numerics (input → L1) | YES (allowed direction) |
| `TYPE_CHECKING: from orpheus.sn.method_space import SNMethodSpace` | geometry → sn (input → L3) | **NO — must be REMOVED** (§3.B/§4) |
| (body) `from orpheus.sn.boundary_realizer import SNBoundaryRealizer` | geometry → sn (input → L3) | **NO — replaced by registry `get(method)`** (the carve's whole point) |

After the carve, the walker imports **only** geometry-internal + numerics symbols → the new
home is layer-legal **iff both sn imports are gone**.

---

## §1. Behavior-preservation gate (the core — bit-identity of the composed operator)

**Claim layer:** flux-shape (operator-output) equivalence; pillar = closed-form /
structurally-independent pointwise reference. **NOT** an eigenvalue claim — no solver runs.

**What it must pin:** `realize_recursively(law, ms)` produces the **same composed
operator** before/after the move, for COMPOSED laws (`0.3*Reflective + 0.7*White`, nested
sums/scales `0.5*(0.3*spec + 0.7*white)`), not just leaves.

**The wall already exists and is sufficient — `tests/geometry/test_law_composition.py`.**
Audited against the §1 requirements:

- **Exercises the COMPOSED path (not just leaves):** YES.
  `test_realize_recursively_lawsum_returns_operator_sum` (structure of `OperatorSum` of two
  `ScaledOperator`), `…_apply_matches_pointwise_weighted_sum` (L1 — value of
  `0.3·spec + 0.7·white` vs the structurally-independent pointwise numpy
  `0.3*spec_realised.apply + 0.7*white_realised.apply`, `nulp=4`),
  `…_walks_nested_depth_first` + `…_nested_apply_matches_distributive_form` (the
  `0.5*(0.3*spec+0.7*white)` nested case). The leaf + type-guard cases round it out
  (`…_leaf_dispatches_to_sn_realizer`, `…_raises_type_error_on_unknown_node`,
  `…_raises_on_ndarray_leaf`).
- **Reference is structurally independent of the SUT:** YES. The apply-equivalence
  reference is the explicit pointwise weighted sum using only numpy `+`/`*` (above the
  trusted-library line per `algebra-of-record`); it does NOT call `realize_recursively`.
  This is a genuine cross-check, not a self-referential snapshot.
- **Regime activates the term:** the composed `OperatorSum`/`ScaledOperator` assembly IS
  the only logic the move touches (the leaf body just swaps its dispatch source). A swap of
  `OperatorSum(a,b)` ↔ wrong order, or a dropped `ScaledOperator(scalar, …)` wrap, reddens
  `…_lawsum_returns_operator_sum` / `…_lawscaled_wraps_in_scaled_operator` /
  `…_apply_matches_pointwise_weighted_sum`. Confirmed live: the suite is **18 passed**
  pre-move.

**Import-path migration (load-bearing for this gate):** line 50
`from orpheus.sn.boundary_realize import realize_recursively` **MOVES** to
`from orpheus.geometry.boundary import realize_recursively` (assuming the walker is
re-exported from the boundary package `__init__`, which it should be — it is the package's
sole descriptor→operator transformer). See §4 table. This file is the **primary
test-migration target**; it is a behavioral (correctness-contract) test → **rewire to the
successor import**, do NOT delete.

**Gate-teeth proof (do at implementation time, AGENT.md §0.5):** the suite is already green
pre-move; after the move, mutate `OperatorSum(a_op, b_op)` → `OperatorSum(b_op, a_op)` in the
relocated walker **in-process (monkeypatch — NEVER `git checkout` the uncommitted file)** and
confirm `…_apply_matches_pointwise_weighted_sum` reddens under `-O`. (The reference is
order-independent for a 2-term sum value but the structural assertions
`op.a.scalar == 0.3` / `op.b.scalar == 0.7` red on the swap — confirm at least one reddens.)

**STRONGER bit-identity option (recommended, cheap):** because this is a pure relocation
with **zero numerical change**, add a **same-process old≡new** equivalence assertion to the
migration commit (or run it once as a throwaway): capture `realize_recursively` from BOTH
`orpheus.sn.boundary_realize` (if a one-cycle shim is kept) AND
`orpheus.geometry.boundary`, and `np.testing.assert_array_equal` their `.apply(psi)` on the
Marshak tree for a fixed-seed `psi`. If the implementer does a **hard move** (no shim — the
preferred aggressive-retirement path), this cross-import check is unavailable; the pointwise
`nulp=4` reference + the structural assertions are then the bit-identity evidence (they
already pin the exact composed structure). Bit-identity here is genuinely expected because
the operator-assembly code is byte-for-byte relocated.

---

## §2. Registry-routing equivalence gate (the carve's defining change)

**Claim:** routing leaf dispatch through `BoundaryRealizerRegistry.get(method)` instead of
the hardcoded `SNBoundaryRealizer()` resolves to the **same SN realizer**, AND the SN
realizer is **actually registered at the moment the walker needs it**.

**This is the gate the brief asked for AND the home of the §TL;DR-2 hazard.** Two distinct
sub-claims — BOTH required:

**§2a. The lookup returns the SN realizer for the SN key (value equivalence).**
A test asserting `BoundaryRealizerRegistry.get("SN") is SNBoundaryRealizer` already exists:
`tests/sn/operators/test_sn_boundary_realizer.py::TestRegistryLookup::test_get_sn_returns_sn_realizer`
(L1). It is sufficient for the *value* claim. ALSO add, to the walker's own suite, an
end-to-end equivalence: `realize_recursively(spec, ms).apply(psi)` ==
`SNBoundaryRealizer().realize(spec, ms).apply(psi)` — which
`…_leaf_dispatches_to_sn_realizer` already does (it compares walker-op vs
direct-`SNBoundaryRealizer()`-op, `assert_array_equal`). Post-move that test silently
upgrades to "walker-via-registry == direct-realizer" — exactly the §2a claim. **No new test
needed for §2a**; confirm `…_leaf_dispatches_to_sn_realizer` still passes post-move.

**§2b. (THE HAZARD) The registry is populated when the walker runs — a fresh-process
registration-timing gate.** A registry **miss** is a silent regression (`get("SN")` raises
→ every composed-BC realization dies). The existing tests CANNOT catch it because they run
in a process where `orpheus.sn` is already imported (process-global class state masks the
miss). **MUST add a new gate that runs in a SUBPROCESS:**

```python
# tests/geometry/test_law_composition.py  (or a sibling)
import subprocess, sys, textwrap

@pytest.mark.l1
def test_walker_resolves_sn_realizer_in_fresh_process() -> None:
    """Registry-routing leaf dispatch must NOT depend on an ambient
    prior `import orpheus.sn`. Runs in a fresh interpreter so the
    process-global BoundaryRealizerRegistry starts empty — the exact
    condition the registration-timing regression hides under in-suite.
    """
    script = textwrap.dedent('''
        from orpheus.geometry.boundary import (
            realize_recursively, ReflectiveBoundary, WhiteBoundary,
        )
        # Build the method space WITHOUT a top-level `import orpheus.sn`
        # at the script head — exercise whatever import the walker itself
        # guarantees. (The ms construction will legitimately import sn;
        # the POINT is the walker must not silently rely on an ambient one
        # that a real non-test caller might not have made.)
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.sn.boundary_realizer import SNMethodSpace
        ms = SNMethodSpace.minimal(Quadrature.gauss_legendre(8))
        law = 0.3 * ReflectiveBoundary(axis="x") + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
        op = realize_recursively(law, ms)
        import numpy as np
        op.apply(np.zeros((8, 3)))   # must NOT raise BoundaryRealizerRegistryError
        print("OK")
    ''')
    r = subprocess.run([sys.executable, "-O", "-c", script],
                       capture_output=True, text=True)
    assert r.returncode == 0, f"fresh-process walker failed:\n{r.stderr}"
    assert "OK" in r.stdout
```

**Design note on the resolution the implementer must choose (and this gate enforces):**
the walker has two correct ways to guarantee the registry is populated, and the gate above
passes for BOTH — pick one explicitly and document it:

- **(Preferred) Keep a body-level import of the SN realizer module purely for its
  registration side-effect** is **NOT** an option in the new home — that re-creates the
  `geometry→sn` edge the move deletes (and the layer gate would flag it). So the geometry
  walker MUST rely on the **caller** having imported the method's realizer before calling
  `get(method)`. The fresh-process gate proves the *real* production callers
  (`SNMesh._resolve_bcs` and friends, which live in `sn/` and import the realizer) satisfy
  this. **The contract becomes: "the walker dispatches `get(method)`; the method's realizer
  module is the caller's responsibility to have imported (the per-method subpackage
  self-registers on import)."** This matches the existing registry docstring
  (`get()`: "typically meaning the consumer didn't import the realizer module yet").
- The `BoundaryRealizerRegistryError` miss message is already excellent (names the missing
  key + `Available: […]`) — **add a negative gate** asserting that a fresh-process
  `realize_recursively(leaf, ms_for_unregistered_method)` raises
  `BoundaryRealizerRegistryError` (positive+negative pairing, `vv` anti-#11). Since SN is
  the only functional realizer, simulate the miss by clearing the entry in a subprocess
  (`del BoundaryRealizerRegistry._registry["SN"]` then call) and asserting the raise — this
  pins the failure mode loudly rather than as an `AttributeError` deep in the sweep.

**Gate-teeth proof:** in a subprocess, `del BoundaryRealizerRegistry._registry["SN"]` before
the walker call → the positive gate (`…_in_fresh_process`) MUST red with
`BoundaryRealizerRegistryError`; the negative gate MUST stay green. This is the mutation
that proves §2b's teeth bite. (Do it in-subprocess; the in-suite registry must not be
mutated — it would corrupt later tests.)

---

## §3. Layer-inversion assertion (the move's structural deliverable)

**Claim:** after the move, `orpheus/geometry/` does NOT import `orpheus/sn/` for the walker
(neither runtime nor `TYPE_CHECKING`).

**§3a. The existing gate IS the right gate — `tests/test_layer_imports.py`.**
- It covers `geometry → sn`: `FORBIDDEN_EDGES["geometry"] = L2_PACKAGES | L3_PACKAGES`, and
  `sn ∈ L3`. The `@pytest.mark.foundation` parametrized `test_no_forbidden_imports` walks
  every `orpheus/**/*.py` (so it will include the **new**
  `orpheus/geometry/boundary/boundary_realize.py` automatically once the file lands).
- **It parses runtime AST imports**, so it is the authoritative check (text-grep is a
  secondary confirmation). It already passes (**291 passed**) pre-move because today the
  walker lives in `sn/` (a legal `sn → geometry`/`sn → numerics` direction).
- **CRITICAL (the §TL;DR-3 hazard):** if the walker is moved **verbatim**, the gate will
  **FAIL** on the relocated file because of its `TYPE_CHECKING: from orpheus.sn.method_space
  import SNMethodSpace`. The gate's TC-tolerance (`:148`) does NOT exempt `geometry`
  (INPUT). **So this gate going RED on the new file is the *expected* signal that the
  `SNMethodSpace` import was not removed** — it is doing its job. The move is complete only
  when `test_no_forbidden_imports` is GREEN for the new file, which forces §3b/§4.

**§3b. Text-grep confirmation (secondary, cheap, run at implementation time):**
```bash
# After the move, NO sn import (runtime or TYPE_CHECKING) may remain in the walker's new home:
grep -rEn "^[[:space:]]*(from|import) orpheus\.sn" orpheus/geometry/boundary/boundary_realize.py
# Expect: NO matches.  (Docstring :class:/:func: cross-refs to orpheus.sn.* are fine — they are
# not imports; but UPDATE them to point at the new location per Cardinal Rule 3 / Sphinx -W.)

# And the whole geometry tree must still pull in zero sn modules at import:
python -c "import sys, orpheus.geometry; print(sorted(m for m in sys.modules if m.startswith('orpheus.sn')) or 'NONE')"
# Expect: NONE.
```

**§3c. Awareness flag (the concurrent `SNMethodSpace` repoint — NOT in this plan's scope).**
The main agent is concurrently repointing
`geometry/boundary/__init__.py:168` (the docstring `from orpheus.sn.boundary_realizer import
SNMethodSpace` → `sn.method_space`). **Both tendrils touch the same `__init__.py`.** Tests
that exercise the combined import surface of that `__init__`:
- `tests/geometry/test_law_composition.py` imports `from orpheus.geometry.boundary import
  (AlbedoBoundary, LawScaled, LawSum, ReflectiveBoundary, VacuumInflow, WhiteBoundary)` AND
  (today) `from orpheus.sn.boundary_realize import realize_recursively`. **After BOTH
  changes**, this file's imports become: the geometry names + `realize_recursively` now from
  `orpheus.geometry.boundary`, and `SNMethodSpace` still from `orpheus.sn.boundary_realizer`
  (re-export) or repointed to `orpheus.sn.method_space`. **Flag:** the implementer should
  land the two changes so that this single file's import block is internally consistent —
  run `test_law_composition.py` after BOTH land, not just after the walker move, since it is
  the one file touching both tendrils. (The `SNMethodSpace` re-export from
  `boundary_realizer` currently still works — confirmed `__all__ = ["SNBoundaryRealizer",
  "SNMethodSpace"]` at `boundary_realizer.py:105` — so a stale import will not break until
  the re-export is actually deleted.)

---

## §4. Test-migration table (retire-means-migrate; coding-standards floor)

Per the retirement-audit-is-three-searches rule (graph callers + text-grep code/tests/docs
+ direct constructors), the migration surface is:

| File / line | Imports today | Classification | Action |
|---|---|---|---|
| `tests/geometry/test_law_composition.py:50` | `from orpheus.sn.boundary_realize import realize_recursively` | **behavioral** (the 16-test correctness wall) | **REWIRE** → `from orpheus.geometry.boundary import realize_recursively`. Keep ALL tests. This is the primary wall. |
| `tests/geometry/test_law_composition.py:51` | `from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace` | behavioral (uses realizer + method space as the reference) | **STAYS** (these symbols remain in `sn/`; the walker leaves, the realizer + method space do not). Re-confirm `SNMethodSpace` still importable from `boundary_realizer` (re-export) OR repoint to `sn.method_space` if the main agent retires the re-export — coordinate. |
| `orpheus/sn/boundary_realize.py` (whole module) | — | the relocated unit | **MOVE** body to `orpheus/geometry/boundary/boundary_realize.py`; **DELETE** the now-empty `sn/` file (no shim, aggressive-retirement default — there is no production importer to break, only the one test, which is rewired). Re-export `realize_recursively` from `geometry/boundary/__init__.py`. |
| `orpheus/sn/boundary_realize.py:119` (`TYPE_CHECKING` `SNMethodSpace` import) | `from orpheus.sn.method_space import SNMethodSpace` | typing | **DROP** in the new home. Retype the `method_space` parameter as `Any` (matching the `BoundaryRealizer` Protocol's `realize(self, law, method_space: Any)`). Principled: the now-method-agnostic walker should not be statically bound to the SN method-space type. (§TL;DR-3.) |
| `orpheus/sn/boundary_realize.py:200` (body) `from orpheus.sn.boundary_realizer import SNBoundaryRealizer` | leaf dispatch | the carve | **REPLACE** with `BoundaryRealizerRegistry.get(<method>)()` (method key resolution — see §2 design note; the method key must come from `method_space` or be a walker parameter). Remove the sn import. |
| `tests/sn/operators/test_sn_boundary_realizer.py:429-470` (docstring + `TestRegistryLookup`) | references `tests/sn/test_boundary_realize.py` (STALE — that file no longer exists; the walker tests are in `tests/geometry/test_law_composition.py`) AND `BoundaryRealizerRegistry.get("SN")` | API-stale docstring + a kept registry test | **UPDATE** the stale `tests/sn/test_boundary_realize.py` docstring pointer → `tests/geometry/test_law_composition.py`. `TestRegistryLookup` STAYS (it is the §2a value gate). |
| `orpheus/geometry/boundary/_base.py`, `__init__.py`, `_composition.py` (docstring `:func:`/`:mod:` refs to `orpheus.sn.boundary_realize.realize_recursively`) | docstring cross-refs | docs | **UPDATE** every `:func:`/`:mod:` ref `orpheus.sn.boundary_realize…` → `orpheus.geometry.boundary…`. An unresolved Sphinx xref renders as plain text with **no `-W` warning** (the retirement-audit doc-blast-radius rule) — so grep these explicitly; the `-W` build will NOT catch a stale `:func:` that resolves-by-coincidence, but WILL catch one that no longer resolves. Run the grep AND the `-W` build. |
| `orpheus/sn/boundary_realizer.py:53,278` (docstring + error-msg `orpheus.sn.boundary_realize.realize_recursively`) | docstring + a runtime error string | docs + 1 string | **UPDATE** both to the new location. The `:278` one is in a *raised error message* (user-facing) — update it so the guidance points at the real module. |

**No API-smoke-only tests** (symbol-exists) were found for the walker — every test calls it
behaviorally → all are rewired, none deleted. **No characterization tests** (the walker has
no documented-limitation delta) — nothing to file under `tests/*/characterization/`.

---

## §5. Gate sequence (run in this order; the move is DONE only when ALL pass)

All `python -O`, SERIAL (`-p no:xdist --timeout=300 -p no:cacheprovider`).

1. **Layer-imports (the structural deliverable + the §TL;DR-3 tripwire):**
   `python -O -m pytest tests/test_layer_imports.py -p no:xdist -p no:cacheprovider -q`
   → must be **291 passed** (the new file included, GREEN ⇒ no `geometry→sn` import
   remains; RED on the new file ⇒ the `SNMethodSpace`/realizer import was not removed).
   PLUS the grep confirmations in §3b.

2. **Registry-routing fresh-process gate (the §TL;DR-2 hazard):**
   `python -O -m pytest tests/geometry/test_law_composition.py -k "fresh_process" -p no:xdist -p no:cacheprovider -q`
   → the NEW subprocess gate(s) green; mutation-confirm teeth
   (subprocess `del _registry["SN"]` → positive reds, negative stays green).

3. **Behavior-preservation wall (the core §1):**
   `python -O -m pytest tests/geometry/test_law_composition.py -p no:xdist -p no:cacheprovider -q`
   → **18 passed → 18 (+ any new fresh-process gates) passed**. Run AFTER both the walker
   move AND the concurrent `SNMethodSpace` repoint land (§3c — the one file touching both
   tendrils).

4. **Named adjacent wall:**
   `python -O -m pytest tests/geometry/test_boundary.py tests/geometry/test_bc_equivalence_snapshot.py tests/sn/operators/test_sn_boundary_realizer.py -p no:xdist -p no:cacheprovider -q`
   → **32 + (sn realizer suite) passed**, unchanged. `test_sn_boundary_realizer.py` confirms
   the SN realizer + `TestRegistryLookup` still green after the move.

5. **CLI pyright (no regression, AGENT.md L10 — assert the EXACT count, never trust it):**
   `npx pyright orpheus/` → must be **≤ 412** (the move retypes `method_space: Any`, drops
   2 sn imports; expect Δ ≤ 0). If it rises, the `Any` retype or a docstring-driven stub
   resolution regressed — investigate, do NOT `# type: ignore`.

6. **Sphinx `-W` (Cardinal Rule 3 — the doc xref blast radius):**
   `sphinx-build -W docs docs/_build/html` → exit 0. Catches a `:func:`/`:mod:` ref that no
   longer resolves after the module moved. (Does NOT catch a coincidentally-resolving stale
   ref — that is why §4's explicit grep of `orpheus.sn.boundary_realize` cross-refs is also
   required.) Nexus graph rebuilds on this build; the `geometry→sn` doc edges should drop.

7. **Full serial regression (the bit-identity / no-collateral floor):**
   `python -O -m pytest -m "not slow" -p no:xdist --timeout=300 -p no:cacheprovider`
   → the **7-and-only-7 pre-existing reds** (`reference_test_execution_env`: ~2452 passed,
   ~9 min; xdist is UNSTABLE for `tests/sn`+`tests/numerics` — SERIAL is mandatory). Any
   8th red is collateral from the move.

---

## Hazard register (carried forward for the implementer)

- **H1 (dominant) — registration timing:** the geometry walker MUST NOT import the SN
  realizer module (that re-creates the deleted edge). It relies on the **caller** having
  imported the method's realizer. Production SN callers live in `sn/` and do import it; the
  §2b fresh-process gate proves it. A non-test caller that calls `realize_recursively` for a
  method whose realizer module was never imported gets a loud `BoundaryRealizerRegistryError`
  (good) — pin that with the negative gate.
- **H2 — the verbatim move trips the layer gate** via the `TYPE_CHECKING` `SNMethodSpace`
  import (geometry is not TC-exempt). Drop it; retype `method_space: Any`. The layer gate
  going RED on the new file is the correct early signal, not a surprise.
- **H3 — process-global registry masks H1 in-suite.** Any §2 validation that does not use a
  fresh subprocess is **vacuous for the timing claim** (Mode-11-adjacent: the assertion
  fires but never exercises the empty-registry condition). The subprocess is non-negotiable.
- **H4 — stale Sphinx xrefs render silently.** `:func:`/`:mod:` refs to
  `orpheus.sn.boundary_realize…` across `_base.py`/`__init__.py`/`_composition.py`/
  `boundary_realizer.py` (incl. the `:278` runtime error string) must be grep-updated; `-W`
  alone will not catch a coincidentally-resolving one.
- **H5 — the "layer inversion" framing is documentation-only at runtime.** The move's real,
  verifiable payoff: (a) removes the one **test** runtime edge `geometry→sn`
  (`test_law_composition.py:50`), (b) co-locates the walker with its `BoundaryRealizer`
  Protocol + registry (architectural cohesion), (c) drops the docstring/`TYPE_CHECKING` sn
  coupling. It does NOT remove a production runtime cycle (there is none). State the payoff
  accurately in the commit so a fresh session does not over-claim.
- **H6 — "unify after two instances" tension (awareness, not a gate).** The walker's own
  docstring (`boundary_realize.py:84-92`) says the cross-method registry generalisation is
  *deferred until a second functional realizer ships* — and today SN is still the only
  functional realizer (MoC/MC/CP/diffusion are `NotImplementedError` stubs). This move
  performs that generalisation **early** (registry routing with one live realizer). That is
  a **design call for the main agent/user**, not a verification objection — but the spec
  flags it because the deferral was a documented, deliberate choice. The verification is
  identical either way; only note that the registry-routing branch is exercised by exactly
  one live method until a second lands (so §2's negative/miss gate is the only thing
  exercising the "method not registered" arm).

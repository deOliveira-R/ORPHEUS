# W-F verification spec — overload realignment + dead-arm retirement

**Task:** P4.5 W-F (#76), branch `refactor/operator-inverse-algebra`.
**Author:** test-architect (proactive pass for an operator-algebra carve
crossing the typed↔packed boundary).
**Status:** framework + ONE pre-written live-arm sentinel (GREEN now,
mutation-confirmed). The DEAD/LIVE arm split is the **explorer's**
deliverable — this spec CONSUMES it ("the arm the explorer marks LIVE" /
"each arm the explorer marks DEAD"), it does not pre-guess it.

---

## 0. The two changes have DIFFERENT proof obligations — never conflate

W-F bundles two edits with **orthogonal** verification surfaces. The
single most important thing this spec asserts:

| Edit | Nature | What proves it | What does NOT prove it |
|------|--------|----------------|------------------------|
| (a) `@overload` realignment (FullField-first) | **compile-time only** — `@overload` has ZERO runtime effect (the stubs live in `if TYPE_CHECKING:`; `apply = _apply_impl` at runtime) | **CLI-pyright oracle** (`npx --no-install pyright orpheus/`, baseline **412**) | ANY runtime pytest run. A green `-O` suite says NOTHING about overload correctness. |
| (b) dead-arm retirement | **runtime** — changes the dispatch table | the regression wall (below) + the live-arm sentinel | CLI-pyright (a dead arm removed cleanly leaves pyright unmoved or improved) |

The carve must NOT mistake a green runtime run for overload-correctness,
NOR a clean pyright for arm-retirement-correctness. They are checked by
**different tools**.

---

## Canonical invocation + baseline (every runtime gate below)

```
.venv/bin/python -O -m pytest <paths> -p no:xdist --timeout=300 -p no:cacheprovider -q -rfE
```

- **`-O` is mandatory** (production strips `assert`; vv Mode 8). Every
  gate named here uses `np.testing.*` / `require` (a `pytest.fail` call)
  / explicit `raise` — never a bare `assert` — so the teeth bite under
  `-O`. (Confirmed for the sentinel; the named legacy gates already
  comply.)
- **`-p no:xdist`** — xdist is UNSTABLE on `tests/sn` + `tests/numerics`
  in the Host `.venv` (worker crashes); run SERIAL.
- **Baseline = the 7-and-only-7 pre-existing reds.** Any W-F runtime gate
  run is "clean" iff the ONLY failures are these 7:
  - 3× SPH-vacuum ULP `test_bc_extraction_matvec`
  - 2× `ymin` mu_y
  - 2× curvilinear sphere snapshots

  A NEW red (an 8th) is a W-F regression. Zero new reds is the bar.

---

## 1. The @overload realignment — CLI-pyright oracle, NOT runtime

**Sites:** `scattering.py` `@overload` stubs ~1452–1462; `fission.py`
~507–513. Both live inside `if TYPE_CHECKING:` and end with the
`def apply(self, x: Any, /) -> Any: ...` shadow; runtime is
`apply = _apply_impl` in the `else:` branch. **No runtime byte changes.**

**Gate (the ONLY proof of (a)):**

```
npx --no-install pyright orpheus/        # CLI is the ORACLE (pyright 1.1.410 present)
```

- **Pass iff the error count does NOT regress above 412.** A drop
  (412 → ≤412) is acceptable and expected if the realignment resolves the
  spurious `reportOverlappingOverload`. A rise is a W-F regression — the
  realignment introduced a real overload conflict.
- **LSP `<new-diagnostics>` are #226 noise — IGNORE them.** The
  in-editor language server mis-resolves these overloads (the
  `reportOverlappingOverload` the realignment targets is itself an LSP
  artifact per the W-F caveat in `[[project-frame-projection-machinery]]`).
  Only the CLI count is authoritative.
- **Mechanically confirm (a) is runtime-inert** before trusting the
  pyright delta: `git diff` the two files and verify every changed line
  is inside the `if TYPE_CHECKING:` block. If a non-stub line moved, (a)
  has become a runtime change and item-2's wall applies to it too.

**Anti-pattern flagged:** crediting a green `python -O` run as evidence
the overloads are correct. It is not — see §0.

---

## 2. The dead-arm retirement — the runtime regression wall

The retirement IS the runtime change. After each arm is removed, the
following must run with **only the 7-and-only-7 reds** (zero new):

**Tier A — the operator-local canaries (run FIRST, seconds-fast):**

```
tests/sn/operators/test_fission_kernel_crosscheck.py
tests/sn/operators/test_scattering_operator.py
```

- `test_fission_kernel_crosscheck.py` — the **0-ULP fission canary** (B.1
  hand-loop correctness, B.2 production-rate equivalence, B.4
  route-through-functional Mode-11 sentinel) **PLUS the new W-F live-arm
  sentinel** (§3). The B.1 row `test_apply_matches_hand_derived_emission`
  feeds `phi` as **bare ndarray** → it directly EXERCISES the fission
  ndarray arm; if W-F retires that arm, B.1 reds. (This is a guard
  *against* mis-classifying the fission ndarray arm as dead.)
- `test_scattering_operator.py` — the scattering **apply-dispatch** suite
  (`TestApplySemantics`, `TestProtocolCompliance`,
  `TestAnisotropicScatteringExtraction`, `TestProducerSideNormalisation`).
  Exercises the typed `AngularFlux` / `ScalarFlux` arms — confirms the
  RETAINED typed arms still dispatch after the base-dispatcher edit.

**Tier B — the broader runtime wall (run after Tier A is clean):**

```
tests/sn/solve/          # solve + eigenvalue paths (the K-loop integration)
tests/numerics/          # 684 numerics tests
tests/sn/                # SN L0/L1 verification (incl. the 88 verification gates)
```

Map to the brief's named counts: solve+eigenvalue (110), numerics (684),
verification L0/L1 (88). The K-loop lives in `tests/sn/solve/` +
`tests/numerics/` (the `power_iteration` integration tests) — these are
the rows that would red if the LIVE arm vanished and the sentinel somehow
did not catch it (defense in depth).

**Structural note for the carve author — the operators are ASYMMETRIC:**

- **`fission.py`** registers a real `np.ndarray` arm (`fission.py:484`,
  `def _(self, phi_arr: np.ndarray)` → `self.kernel.apply(phi_arr)`).
  This arm is **K-loop-LIVE** (proven §3). Its retirement (if the
  explorer ever marks it dead) would break the K-loop — do NOT retire it.
- **`scattering.py`** registers **NO `np.ndarray` arm at all** — only
  `FullField`, `ScalarFlux`, `AngularFlux`, `HarmonicMomentFlux`, plus the
  `@singledispatchmethod` base that raises `TypeError`. So scattering's
  "dead ndarray dispatch" is about the **base-dispatcher fallback** and/or
  the `@overload` shadow, not a registered arm. The explorer's
  classification must say WHICH (a registered arm vs the base fallback)
  for each operator; this spec's per-arm proof method (§4) applies to
  whatever it names.

---

## 3. The LIVE-arm protection — the crux (PRE-WRITTEN, GREEN, mutation-confirmed)

**Claim layer:** this is a **coverage/behavior claim** ("the K-eigenvalue
loop EXECUTES the fission bare-ndarray arm today"), NOT a value claim. It
is **1-group-irrelevant** (it asserts a call happened, not an
eigenvalue), so the fixture's group count is chosen for the *other*
gates' needs (≥3G fissile heterogeneous, reused from the canary), and the
Cardinal Rule (no 1G eigenvalue claim) is not engaged — declared per
`vv-principles §1.5`.

**Pre-written gate:**

```
tests/sn/operators/test_fission_kernel_crosscheck.py
  ::TestFissionNdarrayArmIsKEigenvalueLive
  ::test_keff_solve_executes_fission_ndarray_arm        [@pytest.mark.sentinel]
```

**What it does (the Mode-11 "wrap the internal call" recipe):**

1. Reaches the **registered ndarray implementation leaf** via the
   descriptor handle
   `FissionOperator.__dict__["_apply_impl"].dispatcher.registry[np.ndarray]`
   — NOT the outer dispatcher (wrapping the dispatcher defeats
   type-routing; see the gotcha below).
2. Re-registers a counting wrapper via `descriptor.register(np.ndarray, …)`
   (the registry is a read-only `mappingproxy` — `register` is the ONLY
   supported mutation; `monkeypatch.setitem` raises `TypeError`).
3. Drives a **real keff solve** (`power_iteration(solver_4g, max_iter=60)`).
4. Restores the original leaf in `finally` (manual revert — the registry
   is class-global; there is no monkeypatch hook for it; NEVER leave the
   live dispatch table mutated).
5. `require(calls["n"] > 0, …)` — fires under `-O`.

**Why it is load-bearing documentation:** the function-call sentinel both
GUARDS against deleting the live arm AND documents *which* arm is
load-bearing (the fission `np.ndarray` arm, hit once per outer iteration
via `SNSolver.compute_fission_source` → `self.fission_op.apply(flux) / keff`,
`sn/solver.py:1055`, where `flux` is the bare ndarray from
`power_iteration`).

**Verification status (done in this pass):**

- **GREEN now** under canonical `-O` (4.2 s): the ndarray arm fires 6×
  for a 6-outer-iteration 4g solve.
- **Mutation-RED confirmed two ways, in-process, NO git checkout:**
  - *coarse* — re-register the arm to raise (simulating deletion): the
    K-loop raises `TypeError` → sentinel does NOT pass. ✓
  - *subtle (the real Mode-11 test)* — patch
    `SNSolver.compute_fission_source` to wrap the bare ndarray into a typed
    `ScalarFlux` (the K-loop **routes around** the ndarray arm and
    converges fine, no TypeError): the sentinel REDs **via its own
    counter** (`require(calls["n"] > 0)` → the exact Mode-11 message). ✓
    This proves the gate is a genuine wrap of the production reader, not
    merely catching an incidental exception — a route-around path cannot
    fake the wrap.

**Gotcha banked for the carve author (and the explorer):** you CANNOT
instrument a `singledispatchmethod` by replacing the outer callable
(`FissionOperator.apply = wrapper` or `_apply_impl = wrapper`) — the
wrapper's `__class__` becomes the dispatch key and the input falls to the
base `TypeError` arm. This is the Mode-11 hazard in miniature: the
"obvious" instrumentation silently measures nothing. Always wrap the
**registered leaf** in `dispatcher.registry[T]`, reached through the
class-`__dict__` descriptor.

---

## 4. Per-DEAD-arm proof-of-death — method per arm the explorer marks DEAD

For **each arm the explorer marks DEAD**, the carve uses ONE of two
methods (state which in the carve plan, per arm):

**Method D1 — delete-and-run-green (preferred when cheap).** Remove the
arm; run the Tier-A + Tier-B wall (§2). If the result stays at the
7-and-only-7 baseline, the arm carried no live behavior → death proven.
This is the strongest proof when the wall genuinely covers the arm's
input type (the wall feeds that carrier somewhere).

**Method D2 — zero-production-calls sentinel (when D1 is ambiguous).**
When the wall might not feed the arm's carrier (so a green wall is
*absence of evidence*, not *evidence of absence*), pin death with a
**reverse Mode-11 sentinel**: wrap the arm's registered leaf with a
counter (same handle as §3), run the **full Tier-B wall**, and
`require(calls["n"] == 0, …)` — asserting NO production path reaches it.
A non-zero count means the arm is actually LIVE and the explorer's
classification is wrong → STOP, re-classify, do not retire.

**Decision rule the carve must apply per arm:**

- Arm's input carrier IS fed by some Tier-B test → **D1** (delete, green
  wall is sufficient).
- Arm's input carrier is NOT obviously fed by the wall → **D2 first**
  (prove zero production calls), then D1.

For the `scattering.py` base-dispatcher fallback specifically (no
registered ndarray arm): the "retirement" is whatever the explorer names
(e.g. the `cast`/fallback that handles a bare carrier). If it is the
`TypeError`-raising base `_apply_impl`, that body is the dispatch
**error-reporter** — retiring it would be a behavior change (callers with
a wrong type would get a different error). Treat any edit to the base
`_apply_impl` as D2-gated (prove no production path relies on its current
message), NOT a free delete.

---

## 5. The minimal carve-plan checklist (fuse with the explorer's split)

1. **(a) overload realign** → `git diff` confirms TYPE_CHECKING-only →
   `npx --no-install pyright orpheus/` ≤ 412. (No runtime gate.)
2. **Pre-flight** → run Tier-A; confirm the **live-arm sentinel is GREEN**
   and the baseline is the 7-and-only-7 (so a later red is attributable).
3. **For the arm(s) the explorer marks LIVE** → KEEP. The §3 sentinel is
   the standing guard.
4. **For each arm the explorer marks DEAD** → apply D1 or D2 (§4) per the
   decision rule → run the Tier-A + Tier-B wall → 7-and-only-7.
5. **Final** → full `-O` serial wall (Tier-A + Tier-B) at 7-and-only-7 +
   `pyright orpheus/` ≤ 412. Both tools must pass — they prove different
   halves (§0).

---

## Scope guard

This spec covers W-F ONLY (the overload realign + dead-arm retirement at
the two named sites + the K-loop live-arm protection). It does NOT extend
to W-H (`_DISPLACEMENT_CLS`), the `transport/operators/` relocation, or
any other P4.5 wave. The pre-written sentinel is the only file touched in
this pass; it is additive (a new `foundation`/`sentinel` class in the
existing fission canary) and green on `main`'s current tree.

## Open API note (not a W-F blocker — flagged for the backlog)

The SN **eigenvalue path leaves `history.n_inner = None`** (the rate
measurand for any future SI-convergence-rate claim is unsurfaced on the
K-loop). W-F does not touch this; noted because the live-arm sentinel
drives `power_iteration` and a reviewer may ask "why not assert on
n_inner?" — because it is None on this path. See
`[[si_convergence_rate_verification]]` (the still-open gap). Out of W-F
scope.

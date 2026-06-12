# QA Lessons

Behavioral corrections only. AGENT.md has the V&V hierarchy,
anti-patterns, and error catalog format -- never duplicate here.

---

## L-001 -- Test count is not coverage

20 passing tests (homogeneous exact, conservation, balance,
non-negativity) missed a fundamental 2-term bug in cylindrical DD.
Signature: keff diverging under mesh refinement (1.15 -> 0.90 -> 0.52).

**Rule**: When reviewing "all tests pass" for any solver, first ask:
"Is there a heterogeneous mesh-refinement convergence test?"
If not, the suite proves nothing about the transport operator.

---

## L-002 -- Orphan equation triage class-D before class-B

When closing orphan equations on a Sphinx theory page, scan the
**existing** test suite for evidence the equation is already verified
under a different label BEFORE writing a new test or marking it
documented. The Peierls/case-method/V_α suites have a long history of
adding equations to a theory page and adding tests in a test file
without wiring the `@pytest.mark.verifies("label")` connection. The
audit interprets these as orphans, but the test exists — adding the
label is a 1-line fix that closes the orphan with no new verification
work.

**Rule**: For a new orphan, the search order is D -> B -> A -> C:
1. Class D (existing test, just needs the label) — most common in
   ORPHEUS, scales to 25 percent or more of orphans on the trajectory-
   resolvent / peierls-Nystrom pages.
2. Class B (definitional / derivation step / governing equation) —
   bulk of the rest; mark `.. vv-status: <label> documented` with a
   rationale comment per Cardinal Rule 3.
3. Class A (write a real test) — only when no test covers and the
   equation is a verifiable claim.
4. Class C (stale, remove) — rare; the test suite usually catches
   stale labels via the audit's drop-into-orphan signal.

---

## L-003 -- The matrix.rst orphan list can lag the live RST

`docs/verification/matrix.rst` is auto-generated from
`tools/verification/generate_matrix.py` on every Sphinx build. Until
the build runs after a label rename, the matrix snapshot lists the OLD
labels. ALWAYS re-run `python -m tests._harness.audit --gaps` to get
the live orphan list before classifying — do not trust the matrix.rst
snapshot for label spelling. Observed during the 2026-05-03 78-orphan
sweep with `case-method-eqXX` -> `singular-eigenfunction-eqXX` rename:
the snapshot showed `case-method`, the audit showed `singular-
eigenfunction`. Mass-applying labels from the snapshot would have
created 5 self-inflicted orphans.

---

## L-005 -- Locating slow/timeout tests in tests/derivations

The whole `tests/derivations` suite CANNOT be run in one bounded
process to find the `Timeout (>60.0s)` tests: with `--timeout=60
--timeout-method=signal` the per-test 60s stalls accumulate past any
sane `gtimeout` wall (even `-n 6` xdist hit the 600s cap mid-run and
junit-xml is NOT written when the process is SIGTERM-killed, so the
`-rfE` reason summary is lost). Working method:

1. Split into batches that each COMPLETE: the slow tests cluster
   entirely in `test_peierls_*` files. `test_fn*`/`*la13511*` (13
   unique files) and all 32 non-peierls/non-fn files are fast (0
   timeouts) — clear them as 2 group runs first.
2. For the peierls group, the 4 heavy files (`test_peierls_reference`,
   `_nystrom_verification`, `_convergence`, `_specular_bc`) plus
   `test_peierls_greens_function_cylinder_mr` hold ALL the timeouts;
   run each suspect file ALONE with `-n 8 -q -rfE --tb=no` so the
   per-test `FAILED ... - Failed: Timeout (>60.0s)` reason line is
   captured (only a COMPLETED run writes the short-summary reasons).
3. The other ~46 peierls files run clean as one group.

The 2026-06 sweep found exactly 20 timeout tests, all genuine
`Timeout (>60.0s)` (zero real assertion/exception failures): 5 in
`test_peierls_reference`, 11 in `_specular_bc`, 2 in
`_nystrom_verification`, 1 in `_convergence` (`cp_slab_1eg_2rg`
param), 1 in `_greens_function_cylinder_mr`.

**Param-level precision**: when only SOME params of a parametrized
test time out (sphere passes, cylinder+slab stall), mark just the
slow params with `pytest.param(..., marks=pytest.mark.slow)` rather
than the whole function/method — keeps the fast params in the
default lane. Verify with `--collect-only -m "not slow"`.

---

## L-006 -- Mode-8 (-O strip) classification: rewriter boundary + testpaths

Two facts gate whether a bare `assert` is a real `-O` false-green:

1. **Assertion-rewriter boundary.** pytest rewrites bare asserts ONLY
   in (a) collected test modules and (b) registered conftest/plugins.
   Asserts in `orpheus/` production modules (incl. `orpheus/derivations/`)
   are NEVER rewritten, so under `-O` they are inert NO-OPs. (`np.testing.*`
   / `pytest.fail` are function calls -> fire under `-O` regardless.)
2. **`testpaths = ["tests"]`.** The canonical suite collects ONLY `tests/`.
   `test_*` wrappers that live INSIDE an `orpheus/derivations/*.py` module
   (e.g. `balance.py:test_cartesian_1d`) are NOT collected -> their
   internal asserts run only on a manual `pytest orpheus/.../balance.py`
   (the docstring usage), never by `python -O -m pytest`. Those are class-D
   (dead w.r.t. the suite), NOT class-A.

**Classification recipe** for a bare assert in `orpheus/`:
- Nexus `callers` on the function node. Filter to callers whose id starts
  `tests.` (collected) vs `orpheus.derivations.` (in-module wrapper, dead).
- If a COLLECTED test calls it: read what the test independently asserts on
  the RETURN VALUE. If the test cross-checks the same property against a
  structurally-independent path (e.g. SymPy coeffs vs the DD sweep), the
  internal assert is class-B redundant (the H4 self-reference trap does NOT
  bite because producer and consumer are independent). If the test only
  consumes the return and the assert is the sole correctness gate -> class-A.
- `isinstance(...)` after a `case`/`if coord`-branch, and `x is not None`
  on an Optional the contract guarantees -> class-C type-narrowing (strip =
  designed; downstream AttributeError if ever violated).
- `assert row == n_unk` before `np.linalg.solve` (matrix-assembly row count),
  `if __debug__:` blocks, `assert <closed-form sanity>` guarding `[0]`
  indexing of a `sp.solve` result whose REAL verification is a returned
  `pass_*` boolean -> class-C.
- Import-time `validate_all()` on a HARDCODED data table (xs_library.py:307,
  `np.allclose(sig_t, sig_c+sig_f+sig_s.sum)`) with NO independent
  collected-test coverage and a constructor (`Mixture`) that does NOT
  re-validate -> class-A genuine false-green: the only consistency gate on
  the canonical XS library, silently inert under `-O`.

**The #228 audit** (56 sites): 1 class-A (`xs_library.py:307`), rest C/D.
The original premise (test_keff_2d bare asserts inert) was REFUTED -- those
ARE collected -> rewritten -> fire under `-O`.

---

## L-004 -- vv-status rationale comments must NOT use [brackets]

The `:vv-status: documented` directive lives in the same RST file as
the labelled equation, conventionally as a top-level RST comment
(`.. vv-status: <label> documented`). When attaching a rationale
comment block with `..    [category] description` formatting,
docutils parses each [xxx] as a citation reference, producing
"duplicate citation" warnings under `-W`. Use (parens) instead of
[brackets] in rationale comments, and prefer a single-line `..`
comment over a multi-line `..  / .. / ..` continuation block.

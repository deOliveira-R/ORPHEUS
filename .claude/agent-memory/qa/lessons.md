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

## L-007 -- foundation + verifies(...) is silent level conflation

The harness (`tests/_harness/registry.py`) accepts `@pytest.mark.foundation`
stacked with `@pytest.mark.verifies("<eq>")`: `_existing_level` resolves
`level="foundation"` while `_collect_str_marker_args(item,"verifies")`
SEPARATELY records the equation, so Nexus writes a `tests` edge and the
audit credits the physics equation with the foundation test's parametrizations
(observed: a documented eq showed "6 test(s)" from one 6-param foundation test
that never touches a non-flat ψ). The registry docstring forbids it verbatim
("Foundation tests never carry a `verifies(...)` marker"), but NOTHING enforces
it — collection is silent.

**Rule**: a `foundation` test verifies a SOFTWARE invariant (data-structure /
factory / reflection-index contract); it MUST NOT carry `verifies(<physics-eq>)`.
The tell is a `documented`/representational-identity equation whose ONLY
"coverage" is a foundation test. Check the theory page's `.. vv-status:`
rationale: if it names the REAL verifiers (MMS operator-admission gate,
adjoint bit-identity) and the foundation test is not among them, the marker is
a misleading edge — drop `verifies(...)`, keep `foundation`, reference the eq
in prose via `:ref:` only. Fix is 1 line (delete the marker); re-run
`python -O -m tests._harness.audit` to confirm the eq drops out of the
coverage attribution. (Caught 2026-06 on
`test_coupled_pole_mu_level_invariant.py`, eq `sn-err-058-coupled-pole-continuity`.)

---

## L-008 -- False-xfail under a stale index: verify the FAILURE REASON, not just xfail status

A `xfail(strict=True)` test is "satisfied" the moment it fails for ANY
reason -- the strict gate only checks pass/fail, never WHY. A stale array
index (e.g. `(ng,nx)` flux read as `values[:,0]` -> length-1 cell-0 slice
broadcast against a length-nx reference -> garbage L2 ~14 that DIVERGES)
makes the test a FALSE xfail: green suite, but failing for a reason
unrelated to its documented xfail reason (#229 floor). The W1 review
(2026-06-13) caught the FIXED version: `values[:,0]`->`values[0,:]`.

**Rule**: when a diff touches an array index inside a strict-xfail test,
re-run it with `pytest --runxfail` to surface the REAL failure and confirm
it matches the xfail `reason=` (here sphere orders [1.995,1.999,1.407] +
err[-1]=1.4e-3 = the genuine #229 fine-mesh floor; cylinder orders ~0
err~1.9e-2 = structural floor). Then re-run WITHOUT --runxfail (`-rxX`) to
confirm it stays XFAIL (no strict-XPASS suite break). Same bug class as
the rank-d-carve fallout (`Q_numerical[:,0,:,0]->[:,0,:]`); `(ng,nx)`
1-group MMS is ALWAYS `[0,:]` for the radial profile, `[:,0]` is cell-0.

## L-010 -- Mode-8 (-O) does NOT apply to bare asserts in collected tests/ modules

A bare `assert` in a test file UNDER `tests/` IS rewritten by pytest's
assertion rewriter at COLLECTION time and FIRES under `python -O` — the
interpreter `-O` flag strips asserts in NON-rewritten modules only
(production `orpheus/`, see L-006). The `PytestConfigWarning: assertions not
in test modules or plugins will be ignored` is about non-test-module asserts;
it does NOT mean the collected test's bare asserts are inert. Definitive
probe (run once if unsure): a collected `test_x(): assert 1==2` under
`-O` FAILS. So a W2/W3/W4-style file whose load-bearing gates are all bare
`assert np.all(orders > 1.9)` is SAFE under the canonical `-O` invocation.
Do NOT raise a Mode-8 flag for bare asserts in `tests/` — reserve Mode-8 for
bare asserts in `orpheus/` production paths (L-006 recipe).

## L-011 -- Replicate the test's OWN solve helper before judging a "value" claim

A prescribed-inflow / non-vacuum-BC MMS value-claim depends on `q.boundary`
being wired. A naive `solve_sn_fixed_source(...)` (vacuum default) reproduces
the WRONG number (measured 50% outer-cell error vs the test's 0.26%) — which
looks like the test is lying but is actually the test's own catcher firing
("vacuum inflow misses the A(R)>0 surface term"). ALWAYS import and call the
test module's internal `_solve(case, nc)` to reproduce a value claim; a
divergent hand-replication usually means YOU dropped the BC, not the test.
(Caught 2026-06-13 W3 prescribed_inflow review: test's `_solve` gave max-rel
2.637e-3 at nx=160, matching the docstring exactly.)

## L-012 -- "BC X is load-bearing because k=k_inf" holds ONLY for homogeneous

A directional-eigenvalue test's "vacuum BC is load-bearing — a reflective
sphere has k=k_inf, flux-shape independent, so P1 can't change k" reasoning
is TRUE only for a HOMOGENEOUS reflective medium. For a HETEROGENEOUS
reflective finite sphere the flux is NON-flat, so P1 DOES change k via
spectral/spatial coupling (measured Δ=2.4e-2 reflective vs 1.4e-2 vacuum —
reflective Δ was LARGER, not zero). The k=k_inf control fires cleanly only
on a homogeneous fissile sphere (measured reflective Δ=1.5e-12, machine
zero). When a docstring justifies a BC choice via k_inf, check the config:
if it is heterogeneous, the justification is wrong even though the test's
assertions (on the vacuum config) still hold — a Cardinal-Rule-3 (WHY must be
right) doc-correctness flag, not a test-validity failure.

## L-009 -- "floor scales with quadrature" is a falsifiable gate, not a tautology

When a fix is claimed to clean a convergence RATE but NOT remove a floor,
the honest gate (vv anti-pattern #5/#17) pins the floor as a verified
quadrature-SCALING quantity rather than asserting "floor removed". The W1
`test_w1_aniso_sphere_floor_scales_with_quadrature` asserts
`err(S32,nx=160) < err(S16,nx=160)/2.0`. This is FALSIFIABLE: a fixed
closure-bug floor would be quadrature-independent -> ratio ~1.0 -> gate
FAILS. Only a #229-style interpolation floor (the half-angle thread is
interpolated, scales with angular order) passes (measured ratio 3.42).
A "floor scales" gate is a CLAIM about floor character, distinct from
both the rate gate and a (false) removal gate -- accept it.

## L-004 -- vv-status rationale comments must NOT use [brackets]

The `:vv-status: documented` directive lives in the same RST file as
the labelled equation, conventionally as a top-level RST comment
(`.. vv-status: <label> documented`). When attaching a rationale
comment block with `..    [category] description` formatting,
docutils parses each [xxx] as a citation reference, producing
"duplicate citation" warnings under `-W`. Use (parens) instead of
[brackets] in rationale comments, and prefer a single-line `..`
comment over a multi-line `..  / .. / ..` continuation block.

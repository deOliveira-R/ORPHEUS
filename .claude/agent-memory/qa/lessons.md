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

## L-013 -- Verbatim-relocation claims: prove by NORMALIZED ast-diff, not by re-running gates alone

A "this code moved verbatim with deterministic substitutions X→Y" claim is
provable MECHANICALLY: extract the OLD body + NEW body, apply the claimed
substitutions to OLD (`.replace("self.sn_mesh","self.mesh")` etc.), strip
docstrings/imports/comments/blanks, and `difflib.unified_diff`. A TRUE
verbatim relocation reduces to ONLY: the signature line + the deliberately-
added fork (e.g. an `emit_angular`-guarded block) + the return-shape change.
Any other line in the diff is an unaudited edit hiding inside a "pure move"
claim. This is FAST (seconds) and far stronger than re-running a regression
gate, because the gate may admit drift (see L-014) or have a coverage hole
(see L-015). Used 2026-06-14 on #206 Phase C: the 1-D `_compute_LpC`→
`_apply_walk` move + `_compute_LpC_transpose`→`loss_action_transpose` move
both reduced to signature-only diffs (transpose: ZERO body lines changed).

## L-014 -- A regression gate's HARD floor and its STRICT (bit-identity) floor are different gates

A `kind="direct"` regression assert (`assert_array_almost_equal_nulp(
nulp=reduction_depth)`) HARD-tolerates up to `reduction_depth` ULP. "Bit-
identical (0 ULP)" for a pure-refactor PR is enforced ONLY by the
`-W error::DriftWarning` escalation layered ON TOP. So "the gate passes"
≠ "bit-identical" — verify WHICH invocation ran. Prove the strict floor is
LIVE (not a false gate): perturb the committed baseline `.npy` by 1 ULP
(`np.nextafter`), run with `-W error::DriftWarning` → MUST fail; without it
→ passes with a DriftWarning. Restore via `git checkout --` (np.save appends
`.npy` to a manual `.bak` name — don't hand-roll the backup).

## L-015 -- conftest filterwarnings overrides are SESSION-GLOBAL but do NOT cross-leak to sibling dirs (verify, don't assume)

`tests/sn/regression/conftest.py::pytest_configure` does
`config.addinivalue_line("filterwarnings","always::DriftWarning")`, which
makes `-W error::DriftWarning` INERT for that directory's own iterative DD
snapshots (they emit 100s–10000s ULP drift but never fail under `-W error`).
The fear: this leaks to a sibling gate (`tests/sn/sweep/core/`) co-collected
in the same session → false green. EMPIRICALLY DISPROVEN 2026-06-14: with a
1-ULP-perturbed `sweep/core` baseline AND `tests/sn/regression/` co-
collected, the `sweep/core` A-NEW gate STILL FAILED under `-W error::
DriftWarning` (per-item filterwarnings precedence: the `-W` CLI filter beats
the conftest `addinivalue_line` for items OUTSIDE regression/). So a "the
override leaks" worry is testable, not theoretical — perturb-and-run before
flagging it as a blind gate.

## L-016 -- "branch fires under quadrature Q" claims need a degeneracy probe, not faith

#206 claim 5 asserted the cylinder pure-azimuthal degenerate-ordinate branch
(`|mu_x|<1e-15`, `A_downstream=0`) "fires under the A-NEW matvec[CYL] leg".
FALSE: `Quadrature.level_symmetric(sn_order=2..8)` ALL have `min|mu_x| ≥
0.22` — ZERO degenerate ordinates. The branch is dead code under standard
LS cubature and is exercised by NO current test. Probe in 3 lines
(`np.count_nonzero(np.abs(q.mu_x)<1e-15)`) before accepting a branch-
coverage claim. (Not a regression — the branch was a verbatim relocation,
proven by L-013 — but the EVIDENCE for it is vacuous; flag the coverage gap.)

## L-017 -- Diffusion-limit silent-error: probe the THICK-CELL regime, not just refined mesh

A spatial scheme advertised "diffusion-limit-consistent" but shipped with a
FLAT source (slope source Q̂=0, e.g. LD Increment A #158) is O(h²) AND exact
on linears AND passes a sin-ansatz MMS -- yet SILENTLY loses the diffusion
limit on optically THICK cells. The MMS ladder hides it because every
refinement drives σ_t·h → thin where flat-source LD is fine; the failure
lives at σ_t·h ≫ 1 (coarse mesh on a diffusive medium), exactly where a
practical user runs.

**Probe recipe** (the discriminating oracle is DD, which IS interior-diffusion-
consistent via WDD): fixed coarse mesh, vacuum BC, eps-scaled diffusive
material (σ_t=1/eps, σ_a=eps, c→1, Q=eps), compare DD_mid vs scheme_mid as h
refines. Measured #158-A: DD holds ~0.950 at every refinement; flat-source LD
gives 0.401 at σ_t·h=100 (~58% deficit) and only recovers (0.884) at σ_t·h=12.5.
DO NOT trust a REFLECTIVE-BC infinite-medium probe at c≈1 -- it needs 1e5+
inner iters and both DD and LD read 81.9% wrong from non-convergence (a probe
artifact, NOT physics). Vacuum thick-cell head-to-head vs DD is the clean cut.

**The flag**: this is a SILENT wrong-answer exposure when (a) the docstring
headline claims "diffusion-limit-consistent" / "all four diffusion limits"
(true of full LD, FALSE of the flat-source cut that shipped), (b) the flat-Q̂
restriction is buried in code comments only ("flat (Q̂=0)", "Increment C"),
(c) NO user-facing warning and NO xfail/tripwire guards the interim. A
deferred-to-increment-C limitation needs EITHER a forward xfail tripwire
(strict=False, flips to xpass when C lands) OR a docstring user-warning that
the diffusion limit requires the moment source -- a staging note in a plan
file is NOT enough when the public entry (solve_sn_fixed_source cell_update=)
accepts the scheme NOW.

## L-018 -- "matvec path tested" needs an instrumented call-count, not a round-trip

A batched residual_kernel_batch round-trip test (residual(cell_kernel_batch(q))≈0)
is a SELF-consistency check (both arms share _kernel_terms) -- it is NOT the
L14 leg-2 (matvec correct) or leg-3 (matvec≡sweep). To know whether the
forward matvec is even EXERCISED end-to-end, monkeypatch-count the two kernels
during the solve: #158-A LD MMS solve = 1600 cell_kernel_batch (sweep/solve)
+ 0 residual_kernel_batch (matvec). SI sweeps never touch loss_action. The
matvec runs ONLY under inner_solver='krylov'. Probe it: LD-via-Krylov gave the
SAME flux as LD-via-SI to 4.1e-14 (matvec IS correct + matvec≡sweep holds) --
but NO committed test drives it, so the claim "matvec works" is true-but-
unverified (NEEDS-EVIDENCE: a 1-line inner_solver='krylov' MMS sibling closes it).

## L-019 -- A stress-ansatz flagged by the test-architect is a binding contract, not advice

When the test-architect memo's GATE spec mandates an angularly-non-trivial,
mixed-scale (k=1,3), heterogeneous-2G, a0>0-non-vanishing-at-boundary stress
ansatz AND the shipped MMS test instead uses build_1d_slab_mms_case() (the
canonical sin(πx/L), 1G, homogeneous -- the EXACT Mode-7 simplification bias
the spec said to override), that is a gate-downgrade, not a stylistic choice.
The sin ansatz: vanishes at both faces (BC handling untested), isotropic-in-μ
(no per-ordinate spatial variation in the moments), 1G (flux-shape degenerate
per H1), homogeneous (nulls redistribution per H2). It cannot stress LD's
slope-moment closure. The per-cell linear-exactness oracle (gate 1) IS
structurally independent and non-tautological (sign-flip breaks it by 1.88 vs
1e-12 tol -- verified), so LD is not unverified -- but the L1 MMS leg ran on
the weak ansatz. Cross-check the shipped test's case-builder against the
test-architect's ansatz spec; a mismatch is a flag even when all tests pass.

## L-020 -- "w=½ generic ops are byte-identical to DD's factored form" is TRUE (verify the IEEE micro-fact, not the docstring)

A coefficient-model refactor that replaces DD's factored closures with generic
affine ops parameterized by w=½ CAN be genuinely byte-identical, because mult-
by-0.5 is an exact power-of-2 scaling: `0.5*(a+b) == 0.5*a + 0.5*b` bit-for-bit
for ALL doubles (verified 2M random pairs, 0 differ) -- the single rounding in
`a+b` equals the single rounding in `0.5a+0.5b` (each summand exactly halved,
exponent-shifted). Likewise `QV*inv/0.5 == 2.0*QV*inv` (0 differ). So
`cell_average=(1-w)in+w*out` and `source_emission=QV*inv/w` at w=½ reproduce
DD's `0.5*(in+out)` / `2*QV*inv` EXACTLY. (#158 Inc B, 2026-06-14.)

**The trap**: the production docstring (`affine_closure.py`) CLAIMED the
opposite -- "principled-equivalent, not byte-identical, ~1 ULP, DD snapshots
re-baseline". That is STALE/WRONG: 0 `.npy` snapshots changed in the working
tree, and the sha256 byte-identity gate (`test_affine_carve_bit_identity.py`,
`si_slab_2g_het`) stayed GREEN. A Cardinal-Rule-3 doc-correctness flag, NOT a
numerics flag. ALWAYS resolve a "byte-identical?" dispute by (a) `git status
--short '**/*.npy'` (re-baseline tell) + (b) the sha256/array_equal gate +
(c) the IEEE micro-fact at the python prompt -- never by the docstring's claim.

**Liveness (L-014 applied to sha256)**: prove a sha256 gate is LIVE before
trusting its green -- monkeypatch the touched op to inject `np.nextafter(out,
inf)` (+1 ULP) and confirm the hash flips. Verified the slab psi-sha flips
under a 1-ULP perturbation of `source_emission`, so the green is real.

## L-021 -- Increment B closed the L-018 matvec-coverage gap for LD

L-018 flagged (Inc A) that LD's matvec was correct-but-UNVERIFIED (no committed
test drove `residual_kernel_batch`; SI sweeps only touch the solve kernel).
Inc B added `test_sn_1d_slab_ld_mms_krylov_matches_si` (inner_solver='krylov'),
which an instrumented call-count proves drives `residual_kernel_batch=640`,
`cell_kernel_batch=0` -- the matvec path is now genuinely exercised end-to-end
AND pinned ≡ the SI sweep. When re-reviewing a follow-up increment, re-check
whether it closes a prior increment's NEEDS-EVIDENCE item (the call-count probe
is the verification, not the round-trip self-consistency test).

## L-022 -- re-baseline masking-check: re-run the CONVERTED gate, prove the pre-existing red STILL hard-fails ≫ nulp; characterize drift via git-show OLD-vs-NEW .npy

When a commit converts a STRICT `assert_array_equal` snapshot gate to a nULP
`assert_regression(kind="direct")` AND deliberately leaves some baselines
untouched (claiming the untouched ones carry a SEPARATE pre-existing structural
red), the masking failure mode is: the looser gate silently SWALLOWS the real
red. The decisive check is NOT "the suite is green" -- it is re-running the
converted gate on the LEFT-UNTOUCHED arms and confirming they STILL HARD-FAIL
at a magnitude ≫ the nulp bound (#240 SPH bulk: ~1e15 ULP vs nulp=nx=5; the
conversion did not mask them). A re-baseline that silenced a real red would show
those arms flipping green.
Characterize the drift PRINCIPLEDLY (criterion c) by diffing the OLD vs NEW
snapshot bytes directly: `git show <commit>~1:path.npy > /tmp/old.npy` +
`git show <commit>:path.npy > /tmp/new.npy`, then ULP-diff. This shows the
EXACT re-baseline scope (#240: only 3 SLAB `_apply_bulk` keys changed; ALL
`_apply_boundary`, ALL `_solve_*`, ALL curvilinear, ALL 2-D keys byte-identical)
and proves the boundary-byte-identical claim from the binary itself, not prose.
The LIVE-code-vs-REGENERATED-snapshot ULP is necessarily 0 (snapshot was
regenerated at HEAD) -- it does NOT characterize the drift; only OLD-vs-NEW does.
A near-zero cancellation value can show a large ULP count (#240 seed=2: 64 ULP
at |val|=0.024, absΔ=2.22e-16, every other element exactly 1 ULP) -- inspect the
worst element's magnitude before calling it an algorithmic change; large-ULP at
small-magnitude is a ULP-metric artifact, not a non-associativity bound
violation. Criterion 2 (structural-independence) is the load-bearing one: run
the multi-group analytical k∞ recovery (`test_si_carve_recovers_analytical_kinf`
2eg/4eg) + LD MMS O(h²) -- old-vs-new ULP proximity is necessary-not-sufficient.
Cross-ref [[lessons-L020]] (git .npy status + ULP/sha256 over docstring),
[[lessons-L014]] (HARD nulp floor vs STRICT DriftWarning floor; the streaming
boundary gate is `assert_regression` not strict, but is exactly as strict as
`assert_array_equal` under the canonical `-W error::DriftWarning` -- prove it by
running the snapshot class under that flag: #240 = 18 passed / 0 escalations).

## L-023 -- convention-relocation re-baseline: scan the WHOLE tree for the OLD-convention literal, not just the diff's touched tests

When a kernel's input contract changes convention (#240: `SNMesh.streaming` /
`s_axes` went from pre-scaled `2|μ|/Δ` to RAW `g=|μ|/Δ`, with the scheme now
applying the diamond `2`), EVERY test that hand-feeds the kernel the OLD literal
is now passing physically-wrong input. The diff author re-baselined the
convention-encoding tests in the SAME directory as the code change
(`tests/sn/sweep/core/`) but MISSED 3 sites in a SIBLING dir
(`tests/sn/spatial/test_linear_discontinuous.py:272/303/340`, all
`s_axes=(2.0*mu/h,)`). 2 of the 3 broke (the geometry-cross-checks
`test_group1_equals_group2_flat` + `test_group3_equals_group2_scan_flat`: one
arm feeds the stale literal to `cell_kernel_batch`, the other derives `g` from
`abs_mu`/`V` correctly → divergence). The 3rd (`test_batched_round_trip`)
SURVIVED because it is a self-consistency round-trip (both arms share the stale
`s_axes`; residual at solved ψ̄ vanishes regardless of convention — the L-018
trap: a round-trip does NOT pin a convention).

**Recipe**: after confirming the touched-file gates pass, `grep -rn` the WHOLE
test tree for the OLD-convention literal (`2.0 *mu/ *h`, `2.0 *np.abs.*/widths`,
`s_axes=.*2\.0\*`, `streaming\(.\).*2\.0\*`). For each hit, classify: a
cross-check against a geometry-derived value WILL break (genuine missed
re-baseline → main goes red); a self-consistency round-trip survives but is now
feeding wrong physics (latent stale-convention test → fix for intent). Prove
the kernel is CORRECT (not buggy) by feeding the NEW literal locally and
confirming the cross-check passes — isolates "stale test input" from "code bug".
Prove the breaks are NEW (not pre-existing) by stash-pop: the 2 broke ONLY with
the diff applied; on clean HEAD they pass. (Caught 2026-06-15 #240 Step A.)

## L-024 -- affine-in-σ "value-correct leaf sum" carve: prove teeth bite by DISABLING the override

#240 Step B: `InvertibleOperator.apply` overrides the inherited
`OperatorSum.apply` (leaf sum `L.apply+C.apply`) to single-source σ from C via
`loss_action(self.sigma)`. The matvec is AFFINE in σ in the FORWARD direction
(`M(σ)ψ = streaming_action(ψ) + σ·ψ`), so the leaf sum is value-EQUAL to the
override to ≤2 ULP — a value-correct-by-coincidence twin source, NOT a bug (no
wrong value ever shipped). Verification consequences I confirmed:

1. **Teeth gate MUST be `array_equal` (0 ULP), not allclose.** Only bit-identity
   discriminates leak-vs-override (both are value-equal). PROVE the teeth bite
   by DISABLING the override (rename `apply`→`_DISABLED_apply` so
   `OperatorSum.apply` leaf-sum takes over) → all 7 teeth (4 fwd + 3 transpose)
   FAIL at exactly the predicted ULP (max 1.42e-14 / 7.99e-15 rel = ≤2 ULP).
   Restore byte-exact (sha256 match). This is the L-007/L-022 marker-removal
   masking-check applied to a strict-xfail→pass flip.
2. **NOT a `catches(ERR-NNN)`** — no wrong value shipped → `foundation` gate,
   `verifies(...)` only. The carve says so explicitly and is right.
3. **Migration loud-fail**: a missed caller passing an operator where the σ
   ARRAY is expected → `AttributeError` on `sig_t.shape[0]` / `[None]` (the
   dataclass operator has no `.shape`/`__getitem__`) — NOT a silent wrong-shaped
   array. Confirm by grepping the operator class for `shape`/`__getitem__`.
4. **Structurally-independent ref for a re-baselined Krylov golden = the SI
   golden for the SAME config.** SI rides `solve` (no apply override); Krylov
   rides matvec (override). They agree 3.9e-12 → the NEW Krylov value is CORRECT,
   not merely close-to-OLD (vv criterion 2). The SI golden stays UNCHANGED (apply-
   only blast radius) — that invariance IS the cross-check.
5. **seed0 46-ULP flag = large-ULP@small-mag artifact** (maxabs 3.55e-15 ≡ the
   CYL matvec order; rel ~7e-15). Masking-check (L-022): OLD CYL `.npy` HARD-FAILS
   under NEW code (`46 ULP ≫ nulp=reduction_depth=5`) — proves the re-baseline
   load-bearing; SPH untouched red STILL hard-fails (~1e15 ULP, #195/#209) —
   proves the gate not globally loosened.
6. **WATCH (fragility, not a blocker)**: an `array_equal` slab-apply value-pin is
   SEED-DEPENDENT (seed=7 → 0 ULP, seed=0 → 1 ULP via `TimedFullField.__add__`).
   Passes but brittle; acceptable because the teeth gate owns the structural
   distinction and 1 ULP is FP noise. (Reviewed 2026-06-15 #240 Step B.)

## L-025 -- "no missed site" for a dedup-carve: grep WHOLE tree + cross-ref the PLAN, not the closeout

A "route N inlined duplicates through one op" carve's missed-site check is NOT
satisfied by grepping the diff's touched files. Grep the WHOLE module subtree
(`orpheus/sn/`) for the OLD reconstruction literal (`= 2.0*psi... -`, the LD
`psi... + .../d2` form), then classify EACH residual hit as routed /
deferred-by-design / MISSED. Two failure modes a closeout memo can hide:
1. A deferral can belong to a THIRD category the closeout's bucket doesn't name.
   #240 D1 closeout listed only "scan-recurrence" deferrals (β-source `2ψ̄` at
   `loss_representation.py:1435`) but the genuinely-remaining direct DD `2ψ̄−in`
   at `loss_representation.py:2117` (`_OneDimScanWalk._sweep_direction`,
   curvilinear matvec) was NOT in that bucket. It is correctly deferred — but the
   authority is the PLAN (`issue_240_phase2_step_d_homing.md` scoped D1 to
   *Cartesian* inlined `2ψ̄−in`; the curvilinear-angular-fused thread is the NEXT
   campaign), NOT the closeout. ALWAYS read the plan's D-step scope line + `git
   blame` the deferral comment (here `fde76ac5`, pre-D1 → genuine standing
   rationale, not a paper-over).
2. A "missed" site can already be routed ONE LEVEL DEEP: the Cartesian arm of
   the same `_sweep_direction` (line 2056-2083) calls `residual_kernel_batch`,
   which D1 routed — so its reconstruction IS routed transitively. Check the
   call graph, not just the literal.

Verdict rule: a residual OLD-literal hit is a DEFECT only if (a) it is a direct
`ψ_out=reconstruct(ψ̄,in)` (not a scan-recurrence β-source), (b) not transitively
routed, AND (c) in the carve's DECLARED scope per the plan. Fail any of the
three → documented follow-up, not a blocker. (#240 D1, 2026-06-16: exactly ONE
residual direct DD recon tree-wide = line 2117, in-(c)-fail = deferred OK.)

## L-027 -- prove a routing-predicate fix's negative-test teeth by REVERTING ONLY production

A "close-the-misroute" change (narrow a strategy `supports()` predicate so a
mesh stops selecting the wrong sweep rep) ships negative tests that assert the
misroute is GONE. Anti-pattern #11 demands the negative test could have FAILED
against the buggy code; for a routing predicate the cheapest proof is:
`git stash push -- <production files only>` (leave the NEW tests in the working
tree), then run the new tests against the reverted-to-pre-fix production. The
negative tests MUST go red AND the red message must NAME the original bug (not
just `AttributeError`); the strategy-free trait probes correctly go red with
`AttributeError: no attribute '<trait>'` (the trait did not exist pre-fix) —
that is the EXPECTED shape, not a flaw. `git stash pop` to restore. (#240
D5-0, 2026-06-16: all 7 new tests RED pre-fix; `test_2d_ld_sweep_raises_not_
silently_dd` red = "did NOT raise — silent return = ran inline DD" = the LIVE
silent-DD hole proven, not asserted.)

Bit-identity claim for a routing-only change: the load-bearing gate is
DESELECT-the-new-tests → pre-existing count UNCHANGED (the predicate touches no
computed flux). A directory-scoped strict gate's TOTAL legitimately grows by
+N_new (new tests live in the gated dir); the invariant is "no PRE-EXISTING
test's value moved", verified by deselection, NOT "total count unchanged"
(#240 D5-0: full 512/1/4, deselect-7 → 505/1/4 = the real proof).

Pyright adjudication for a docstring/ClassVar-only diff: prove ZERO net-new
diagnostics by capturing pyright on the touched files at pre-fix AND post-fix
(stash the production), line-number-STRIP (`sed -E 's/:[0-9]+:[0-9]+//'`), and
diff. Identical-modulo-line-shift = the inserted docstring block shifted every
pre-existing diagnostic by exactly its line count (e.g. +14 from a 14-line
ClassVar block); the diagnostic SET is unchanged → all pre-existing, no
regression. A file that contributes ZERO CLI diagnostics standalone means the
user's "import-unresolved" / "not accessed" items are IDE/cross-tree config
artifacts, not blockers.

## L-004 -- vv-status rationale comments must NOT use [brackets]

The `:vv-status: documented` directive lives in the same RST file as
the labelled equation, conventionally as a top-level RST comment
(`.. vv-status: <label> documented`). When attaching a rationale
comment block with `..    [category] description` formatting,
docutils parses each [xxx] as a citation reference, producing
"duplicate citation" warnings under `-W`. Use (parens) instead of
[brackets] in rationale comments, and prefer a single-line `..`
comment over a multi-line `..  / .. / ..` continuation block.

## L-028 -- the ÷D vs ×(1/D) re-baseline: a "byte-identical" coefficient-model premise fails where the consumer still divides; verify the ORDERING before trusting a snapshot regen

When a fold routes a leftover-inline path onto an established coefficient model,
the fold is byte-identical ONLY where the consumer already uses the
`×inverse_denom` reciprocal form. A consumer still on `÷denom` DIVISION (a
leftover inline path) re-baselines ~1 ULP when it joins the model — division and
reciprocal-then-multiply are NOT IEEE-bit-identical (`2*X/D ≠ 2*X*(1/D)` at 1
ULP; verified `cartesian_scan_coefficients` reproduces the OLD inline alpha/beta
to max abs diff 2.22e-16 = exactly 1 ULP, NOT array_equal). #240 D5a: the spec's
"2-D SOLVE stays byte-identical" premise was WRONG — the pre-D5a
`ScanMarch._sweep_interior` was the ONE remaining `÷D_row` path (the 1-D
CumprodScan already rode `×inv`), so the SOLVE re-baselined too (BOTH
`si_2d_p1_aniso_het` AND `krylov_2d_p1_aniso_het` golden sha moved; the slab/1-D
sha UNCHANGED = the true negative control). The method-implementer's
load-bearing finding (caught a wrong spec premise) was CORRECT — confirm it by
direct algebraic inspection: replicate the OLD inline `alpha = 2*sx2/D - 1`,
`beta = 2*(Q+sy2*psi_y)/D` and the NEW `a, inv, w, (c_y,) =
cartesian_scan_coefficients(...)` + `source_emission(Q + c_y*psi_y, inv, w)` at
controlled input → max rel diff ~2e-16 (`c_y == 2*g_y` byte-exact, `w==0.5`).

**Snapshot-regen ORDERING masking-check** (don't let a re-baseline launder an
EARLIER untracked drift): the load-bearing claim is "pre-fold LIVE matvec ==
FROZEN snapshot at 0 ULP" — i.e. the frozen IS the correct pre-fold reference.
PROVE it: `git stash push -- <PRODUCTION + the regen'd snapshot>` (KEEP the new
test arms) → run the new arms vs the OLD code + OLD frozen snapshot under
`-W error::DriftWarning` → MUST pass at 0 ULP (proves frozen ≡ pre-fold live).
Then OLD-code + NEW(regen) snapshot → MUST hard-fail (proves the regen is real +
the gate is live). #240 D5a: 3 cart2d arms passed strict pre-fold, hard-failed
"6 ULP (max 256)" with the swap = ordering holds. OLD-vs-NEW snapshot ULP-diff:
relΔ ~1.2-3.6e-16 = 1 ULP of the O(1) field (maxabs 1.8-3.6e-15 @ |val|~17-75);
the 256-ULP metric is the L-022 large-ULP@small-mag near-zero-cancellation
artifact. Boundary trace + the `_LpC_` key + ALL non-cart2d (slab/curvilinear)
keys stayed byte-identical (0 ULP) — blast radius = the 2-D row-march ONLY.

**The STRICT-FROZEN docstring stales silently on a re-baseline** (CR3 / L-020):
`test_bc_extraction_2d.py::test_vacuum_bulk_bit_identical` uses
`np.testing.assert_array_equal` (STRICT) with a docstring claiming "must not move
a single bit" / "must stay frozen" / "E0-T1 proved bit-identical to the pre-carve
path". D5a regenerated its `.npy` baselines (relΔ ~1.5e-16) but did NOT touch the
test file → the "must stay frozen" WHY is now FALSE (the output is no longer
bit-identical to the pre-carve path; it is strict-against-the-POST-D5a value).
Gate functions correctly; rationale prose lags. Flag as a doc-correctness nit,
not a blocker. ALWAYS check: a silently-regenerated strict `.npy` baseline whose
CONSUMER test file is untouched → grep the consumer docstring for "frozen" /
"bit-identical to pre-carve" / "must not move a single bit" and flag the stale WHY.

**The two-paths oracle's analytical anchor is TRANSITIVE, not direct** (L14): the
D5a.1 `test_scan_march_equivalence` asserts `ScanMarch.sweep ≡ FullFieldWavefront`
via `assert_allclose` = TWO-PATHS-AGREE. The analytical `k_inf`/`φ=Q/Σ_t` ground
is reached SEPARATELY (`test_keff_2d::TestHomogeneousExact::test_homogeneous_exact`
pins the ScanMarch DEFAULT path to `νΣ_f/Σ_a` ≥2G; the G6 closed-form anchor in
`test_scan_march_end_to_end` runs WINDOW-forced, NOT ScanMarch). A closeout that
says "the oracle pins analytical k_inf=1.875" is LOOSE — the oracle is
transitively pinned; the direct anchor lives in a different file. Confirm the
anchor file is GREEN before crediting the oracle with analytical grounding.

## L-029 -- a re-encoded closed form is NOT L11-circular when compared against an INDEPENDENTLY-ASSEMBLED primitive (the d=1-reduction-to-production oracle)

L-026 said a token-for-token copy of a production formula is circular as a VALUE
check (a sign-flip in prod propagates into the test, stays green). The DISTINCT
case here (#240 D5b-S1, the SymPy UBLD): the test's RIGHT side re-encodes the
production `_schur_terms` S/eff_source/.slope (test lines 346-351, verbatim of
`linear_discontinuous.py:332-335,258-259`) — BUT the LEFT side is the symbolic
primitive's d=1 reduction obtained by `A⁻¹R` of a SEPARATELY-built Kronecker
matrix (`assemble_ubld([h],[mu],...)` → `LUsolve`), NOT a re-statement of
`_schur_terms`. So `diff_psi_bar==0` proves "the production Schur scalar EQUALS
the independently-assembled 2×2 solve" — one side is genuinely structurally
independent. A sign-flip in `_schur_terms` would NOT propagate into the Kronecker
assembly → the oracle WOULD catch it. The circularity test is: **does the bug
live on BOTH sides of the diff?** Re-encoded-formula-vs-re-encoded-formula =
circular (L-026); re-encoded-formula-vs-independent-construction = legitimate.

The genuinely-independent anchor in the same gate is `test_d1_symbolic_primitive_
matches_production_update`: evaluates the symbolic d=1 ψ̄/ψ_out at concrete 2-group
het numbers and asserts `LinearDiscontinuous().update` (the LIVE running algebra,
not a copied formula) reproduces them ≤1e-12 — closes the loop to production.

Sub-claims to ALSO check on a "the production X-view equals the d=1 reduction"
oracle: `diff_face` (`psi_out` vs `downstream_face_trace`) is BOTH sides the SAME
solve vector → a trace-operator-consistency check, NOT a structural-indep value
check (fine, it's a foundation closure-consistency claim, not a value claim — do
not credit it as independence). The ÷V `_kernel_terms` and ×V `affine_scan_
coefficients` views BOTH reduce to the same independently-assembled LEFT (verified
by transcribing prod lines 443-453 / 564-571 myself → diff 0); `a_source_indep`
= `Qbar not in a.free_symbols` is a real structural property of the transmission.

**Fast mutation-probe for a symbolically-SLOW gate**: when a foundation gate's
`sp.simplify(diff)` is pathologically slow on the MUTATED (garbage) expression
(#240 D5b: d2 exact-on-bilinear `simplify` of the |μ_axis|-dropped residual hit
the 400s pytest timeout — a `simplify` perf artifact, NOT gate evidence), DON'T
wait on the full pytest. Apply the mutation, then call the `derive_*` builder with
the params as CONCRETE RATIONALS from the start (no symbolic LUsolve blowup) and
read `diff` at concrete numbers: #240 mutated d2 residual = [0.596,-0.396,0.179,
-0.226] (manifestly non-zero on all 4 moments) → `is_zero_matrix` False →
`_require_zero_matrix`→`pytest.fail` → test FAILS. That is DECISIVE and seconds-
fast; the slow `simplify` only confirms the same non-zero. The d1 tests staying
GREEN under the mutation (proven: `test_d1_symbolic_primitive_matches_production_
update` passed 1.45s under mutation — d1 routes inflow inline via
`mu*fin_trace_weight()*psi_in`, never through `assemble_inflow_axis`) IS the
"d=1 oracle is blind to a per-axis factor" evidence (ERR-060 H2 multi-D analog).
ALWAYS revert the mutation (both return branches here) + re-run green (6 passed)
before closing. (Reviewed 2026-06-16 #240 D5b-S1 Branch 1 — VERDICT: all 6 claims
SUPPORTED.)

## L-031 -- #240 D5b-S2 (d≥2 UBLD LD wiring): a self-consistency round-trip + an A==A pin can BOTH be blind to the bug their docstring/marker claims; the genuine catcher is the continuous-PDE exact-on-bilinear oracle; and "matvec twin verified" ≠ end-to-end Krylov

D5b-S2 wires the d≥2 bilinear LD kernel onto the wavefront (`cell_kernel_batch`/
`residual_kernel_batch` route `len(s_axes)≥2` → `_ubld_system`+`per_cell_solve`;
`_ubld_outgoing_faces` sums the o_a=0,1 Kronecker blocks; faces in `2^{d-1}`
transverse order). DD/Step held byte-identical (513-strict gate UNCHANGED ==
S1 baseline; `tail=() if n_face_moments==1` → no length-1 axis appended; gate is
DIMENSION `len(s_axes)>1` AND trait `spatial_basis_per_axis>1`). Moment-ordering
out-face↔inflow consistency VERIFIED by hand (`_ubld_outgoing_faces` == manual
Kronecker trace; inflow consumer accepts same `2^{d-1}` object per axis, x↔y
non-square). d=1 closed form `==` dense `per_cell_solve` reduction (L29 anchor holds).

1. **MISATTRIBUTED `catches(ERR-NNN)` — the decisive finding.** The NEW
   `test_d2_assembled_matrices_match_symbolic` (entry-wise A/M/G/F_out numpy↔SymPy
   pin) carries `@catches("ERR-060")` but is BLIND to ERR-060: ERR-060 was the
   dropped `|μ_axis|` factor in `assemble_inflow_axis` (the INFLOW assembly), and
   the A==A pin checks `assemble_ubld` (the CELL matrices — A/M/G/F_out contain NO
   inflow factor). MUTATION-PROVEN: drop `|μ_axis|` in `_ubld.py:254` (`return out`
   instead of `mu_axis[...,None]*out`) → of the 3 `catches(ERR-060)` tests, ONLY
   the numpy `test_d2_exact_on_bilinear` fails; the A==A pin PASSES (matrix count
   2→3 inflated by a non-catcher). The pin IS a legit Mode-3 structural pin (a
   dropped streaming factor in G IS caught — verified by a hand-built buggy-G), but
   it does NOT catch the specific bug its marker claims. FIX = drop the
   `catches("ERR-060")` marker (keep `foundation` + the docstring's "catches a
   dropped/mis-scaled factor in the Kronecker assembly" claim, which is true);
   follow-up nit, NOT a blocker (the genuine catcher is correctly marked). RULE:
   for any `catches(ERR-NNN)` on a NEW test, MUTATE the exact documented bug and
   confirm THIS test (not just SOME test) goes red — a marker is a coverage CLAIM,
   L-007 applied to `catches`.
2. **Round-trip (D5b.1) is self-consistent → blind to inflow bugs (L-018/L-023
   reconfirmed in d≥2).** `test_residual_zero_at_solved_cell_avg_2d` PASSED under
   the |μ_axis| mutation (solve+apply share `_ubld_inflow` → the dropped factor
   cancels). It correctly feeds the FULL `psi_bar=psi_avg` (4 moments, not partial)
   + NON-flat per-axis `psi_in` (slope active) + ng∈{1,2} het → it pins solve≡apply
   SAME-system + the matvec-twin face reconstruction (both axes, x↔y-asymmetric
   replay: residual ~2e-16, out_x==rout_x/out_y==rout_y, matvec non-trivial at a
   different probe). But the VALUE-correctness of the multi-D kernel rests on
   `test_d2_exact_on_bilinear` (the L11/L26-clean continuous-PDE oracle: source
   `Q=Ω·∇ψ+Σ_tψ` from a known bilinear ψ via SymPy, asserts solved moments == exact
   Legendre projections, cross-moment xy active d=0.9 — a kernel sign-flip would NOT
   propagate into the SymPy projections) + the MMS smoke. Sound: the indep oracle IS
   committed and IS the genuine ERR-060 catcher.
3. **"Matvec twin verified" is KERNEL-LEVEL only; end-to-end Krylov RAISES
   (deferred, honest).** Brief/spec D5b.4 wanted "Krylov≡SI on the 2-D LD path",
   but `_CellResidual.cell` RAISES `NotImplementedError` for d≥2 LD (the matvec walk
   needs the `2^d`-moment spatial iterate = S3). PROVEN: 2-D LD SI solve works
   end-to-end (finite flux); `inner_solver="krylov"` RAISES loudly (deferred). This
   is the CORRECT interim per L-017 (loud raise, NOT silent-wrong) — the kernel-level
   `residual_kernel_batch` matvec IS verified (round-trip + asymmetric replay), and
   the raise blocks an accidental wrong Krylov answer. L14 leg-3 (matvec≡sweep
   END-TO-END) is genuinely deferred; the kernel twin is the strongest claim S2 can
   make. NOT a blocker (loud-fail interim); the spec's D5b.4 "Krylov≡SI" wording
   over-reaches S2's actual scope — the SHIPPED scope (kernel twin + raise) is honest
   and the test file documents it.
4. **Smoke MMS honestly scoped (Brief Q2 = SUPPORTED).** `test_ld_2d_converges_
   second_order_smoke` is `@l1` NO `@verifies`, checks BOTH rate (`orders[-1]>1.85,
   all>1.7`) AND a value band (`1e-8<err[-1]<1e-2`) → NOT a Mode-5 rate-only false
   green. The absence of `verifies("ld-cartesian-2d")` is CORRECT: the isotropic
   sin·sin ansatz (`build_2d_cartesian_heterogeneous_mms_case`) under-stresses the
   bilinear slope (Mode-7) + S2 threads only the average source moment (Q̂=0); the
   real flux-shape claim (strengthened μ-non-trivial ansatz + Q̂≠0 moment source +
   non-vanishing boundary) is deferred to S4 — a SOUND deferral, NOT a hole, because
   the kernel value-correctness is independently pinned by the exact-on-bilinear
   oracle (the slope IS exercised there, cross-moment active). `ld-cartesian-2d`
   correctly NOT minted (S4/D6). DD≠LD routing-flip (D5b.5) is 1G but legit (a
   discrimination contract, not an eigenvalue claim; DD≠LD is structural).
5. **Pyright (Brief Q6): exactly ONE net-new diagnostic** (pre/post stash-diff,
   rule-level JSON, 51→52). `linear_discontinuous.py:556` `_ubld_inflow` returns
   `np.ndarray|None` (`R=None` seeded, accumulated in a `for a in range(d)` loop
   pyright can't prove non-empty). REAL but benign type nit (d≥1 ⇒ loop always runs
   ⇒ never None at return; SI solve confirmed live). Pattern-3 nit, follow-up: seed
   `R` with the first term. The other 51 (`AngularSourceSink`/`PoleAngularClosure`/
   `sweep_graph` None-subscript) are ALL pre-existing/rooting-noise. Mode-8 (Q7):
   clean — the 2 prod bare asserts in touched files are PRE-EXISTING (scheme.py:499
   is in a DOCSTRING; loss_rep:2238 pre-existing type-narrow); new tests use
   np.testing/pytest.fail; test-file bare asserts are in `tests/` (rewritten, L-010).
   VERDICT 2026-06-16: SUPPORTED-WITH-CONCERNS — all numerics sound; 1 marker
   misattribution (follow-up) + 1 pyright nit (follow-up) + spec-vs-shipped scope
   wording on D5b.4 (honest as shipped). NO false-green, NO blocker.

## L-032 -- "construct-general-only" capability-addition: the byte-id gate needs teeth on BOTH the auto-select AND the phantom-length-1-axis mistake

#240 D5b-S3-A0 minted a typed `SpatialMomentSpace` factor + an OPTIONAL
`spatial_moments: int = 1` param on the flux/source field-space factories
(`AngularField`/`ScalarField`/`HarmonicMomentField`), DEFAULT-OFF. The load-
bearing claim = byte-identical capability addition (DD/Step/LD all unchanged in
a live solve, no production field carries the axis yet). VERDICT: SUPPORTED (all
7 brief questions). The transferable verification pattern for a
construct-general / select-narrow capability addition:

1. **Two DISTINCT teeth-proofs, not one.** A "capability default-OFF" gate can
   leak two ways: (a) the factory AUTO-SELECTS the wider shape (auto-reads
   `mesh.scheme.spatial_basis_per_axis` → LD silently widens), (b) the gate
   appends a PHANTOM length-1 axis at default (re-associates a downstream
   reduction even though "nothing widened"). PROVE BOTH bite:
   - (a) MUTATION: force the helper to auto-read the scheme → ONLY the `[ld]`
     byte-id arms red (DD scheme reads 1, stays green) = the gate discriminates
     auto-select. (#240: `test_*_default_byte_identical_all_schemes[ld]` ×2 red.)
   - (b) MUTATION: make the "append iff >1" policy append at n==1 too
     (`return (n,)` instead of `() if n==1 else (n,)`) → byte-id arms red for
     ALL schemes (DD+LD) = the gate discriminates a phantom axis. (#240: 9 red.)
   The negative-control assertion that catches (b) is `not hasattr(field.space,
   "factors")` — a default field must be a BARE `FunctionSpace`, not a length-1
   `TensorProductSpace`. Independently confirm DD `(24,2,3,4)` == LD `(24,2,3,4)`
   at default DESPITE LD's `spatial_basis_per_axis==2` (the construct-general
   proof: the scheme SAYS 2 but the factory does not read it).

2. **The "append iff >1" policy MUST be single-sourced** (Pattern 7). #240's
   `spatial_moment_tail` delegates to `_ubld.face_moment_tail` (`() if n==1 else
   (n,)`); the cell-tail and face-cochain-tail can never disagree. Verify the
   delegation by reading both + the `AVERAGE_MOMENT=0` constant the space's
   `average_moment_index` surfaces (NOT a re-spelled `0`).

3. **Einsum "spectator-broadcast" lift (`fc->gc` ⇒ `fc...->gc...`) is provably
   byte-identical at rank-2-exact** AND correct as `Σ⊗I_spatial`: at the python
   prompt, `np.array_equal(einsum('fg,fc->gc',...), einsum('fg,fc...->gc...',...))`
   on rank-2 input = True (the `...` matches nothing → no axis); on rank-3 input,
   `einsum('fg,fc...->gc...')` == per-moment-independent stack (each spatial
   moment scattered independently). The IEEE micro-fact resolves the byte-id
   dispute, NOT the docstring (L-020). (#240: all 3 `material_xs_field.py`
   einsum lifts — apply_p0/apply_n2n/legendre_moments — array_equal at DD/Step.)

4. **The strict gate baseline is 513 (not 562).** The S3 crosswalk's "562/2skip"
   was a STALE-PLAN figure; the live S2/S3 baseline is 513P/1skip/4xf under
   `-W error::DriftWarning` (matches L-031). Re-confirmed; no golden moved
   (`git status --short '**/*.npy'` = empty). When a closeout and a plan
   disagree on a count, RUN the gate — the closeout's 513 was right.

5. **Adding a dataclass FIELD (not just a space factor) to a Flux leaf ripples
   to its Displacement sibling.** `FluxRole._mint_displacement` (`φ⊖φ`) copies
   EVERY init field → `MomentDisplacement` needed the same `spatial_moments`
   field or `φ⊖φ` raises TypeError. Verify the affine round-trip (`φ⊖φ→disp`,
   `φ+disp→φ` array_equal) at BOTH default AND widened. (#240: both exact.)

6. **Pyright net-new = 0** proven apples-to-apples (L-027/L-030): run the SAME
   5 touched files PRE (stash prod + hold the new untracked file) and POST,
   path+line-strip, `comm -23`. All 8 errors pre-existing (`DualSpace.of`
   return, `MaterialXSField` Optional, `from_face_arrays` Optional layout,
   `_check_partner` `other.L` on object); the NEW file alone = 0 errors; the new
   `find_factor` RAISES KeyError (not returns None) → type-clean. The brief's
   worry about "find_factor-returns-object at space.py:521" was a brief
   mis-attribution: :521 is `DualSpace.of`, pre-existing. (VERDICT 2026-06-16:
   SUPPORTED, no blocker, no follow-up.)

## L-033 -- #240 D5b-S3 (LD diffusion-limit completion, ERR-061) VERDICT SUPPORTED-WITH-CONCERNS: the structural-independence ground is GENUINE (L11 clean) but the d≥2 matvec ships a LATENT pure_z twin-path crash no committed gate drives

D5b-S3 = the unified all-d LD moment matvec + ERR-061 frame fix (slope ψ̂_n
stored in per-ordinate SWEEP frame, summed by consumer as GLOBAL → backward
ordinates CANCEL the forward slope → φ̂ 6× under-driven → diffusion limit lost;
fixed by `octant_moment_frame_signs` = ∏_a sign_a^{o_a} involution via `_reframe`).

1. **The headline correctness claim is REAL (L11 clean, NOT L4).** The
   from-scratch LM-1989 solver (`_independent_ld_slab` in
   `tests/sn/spatial/test_ld_slope_frame.py`) is GENUINELY structurally
   independent: hand-built cell 2×2 `[[σh+μ,μ],[-μ/θ,σh+μ/θ]]`, hand SI, NO
   ORPHEUS kernel. Verified live: sweep-frame=1.4717 (== ORPHEUS pre-fix
   bit-for-bit), global-frame=2.3080 (rel 2.3% vs ANALYTICAL diffusion 2.362).
   Anchor is the closed-form diffusion VALUE (`@foundation`
   `test_independent_ld_global_frame_recovers_diffusion`), NOT "LD≈DD". The
   production LD value 2.30798 == the independent solver's global-frame value 4dp
   AND matches analytical diffusion (2.3%) BETTER than it matches DD (4.1%) → the
   chain prod-LD≡indep-LD≡analytical closes. The frame primitive verified at the
   prompt: ∏ closed form exact, DD→None, genuine involution (s·s=1 ∀octant), d=2
   x̂y flips iff ODD # axes reverse.
2. **`catches("ERR-061")` markers MUTATION-VERIFIED.** Neuter `_reframe`
   (`return arr` — single mutation, `loss_representation.py:147` imports the SAME
   helper so ALL paths hit) → the 3 prod-path catchers (slope-frame consistency,
   thick-diffusion 1G, thick-diffusion 2G) ALL go RED with the ERR-061 mechanism
   in the err_msg; the `@foundation` independent ground STAYS GREEN by design (no
   `catches`, doesn't ride `_reframe`). Correct per anti-pattern #11 / L-007.
   (catalog entry lists only the 2 thick tests; the slope-frame test ALSO carries
   `catches("ERR-061")` — minor catalog omission, not wrong.)
3. **Mode-7-at-primitive resolved** (the brief's load-bearing concern): the S1
   `assemble_ubld` exact-on-bilinear oracle nulled the diffusion slope; the NEW
   thick-cell tripwire (σ_t·h=10, c=0.99, COARSE nx=4 — NOT refined mesh, the
   L-017 thick-cell probe) + the 2G-het Mode-6 companion (asym SigS
   [[30,9.6],[0,39.6]], both groups recover, g0→g1=9.6/g1→g0=0 → transpose-
   sensitive) genuinely exercise the regime. ERR-061 catalog entry COMPLETE
   (mechanism + fix + #1+#6 classification + how-it-hid + lesson + bug-signature
   x-link).
4. **DD/Step byte-identity = TRUE** (negative control): GATE 4 513/1/4 under
   `-W error::DriftWarning` IDENTICAL pre/post; `git status --short '**/*.npy'`
   EMPTY (no golden moved). `face_moment_tail(1)==()` + `octant_moment_frame_signs(_,1)==None`
   → DD never grows the moment axis. The 7 spatial+operators reds are GENUINELY
   pre-existing (git-stash: at clean HEAD the same 7 + 2 fix-dependent tests fail;
   with diff exactly 7; the 7 = sphere matvec ×5 + 2-D mu_y BC ×2, none touch the
   slab scan / moment frame).

⚠ **THE CONCERN (the d≥2 matvec twin-path gap — a real follow-up).** "matvec≡sweep
for BOTH d=1 and d≥2" is NOT fully verified. d=1 has committed gates (scan≡DAG
+ Krylov≡SI). d≥2 has NO committed end-to-end Krylov≡SI gate (grep krylov in the
d≥2 MMS files = EMPTY; the closeout's "rel 4.99e-11" was an UNCOMMITTED manual
smoke). The d≥2 raise WAS retired (`_CellResidual.cell` comment "d≥2 raise is
RETIRED") so the path is live — but RUNNING it on the MMS case quad (N=110, **2
pure-z ordinates**) CRASHES: `loss_representation.py:742` matvec `pure_z` does
`LpC[oct_idx]=sigma*probe[oct_idx]` with NO moment-axis broadcast guard → ValueError
`(2,6,6) vs (1,2,6,6,4)`. The SWEEP `pure_z` (line 654-655) HAS the guard
(`if q.ndim>sig.ndim+1: sig=sig[...,None]`); the matvec twin does NOT — a Pattern-2
twin-path asymmetry. Latent because the d≥2 verification uses SI (smoke) + sweep-vs-
sweep two-paths on level_symmetric (ZERO pure-z); the matvec `pure_z`+moments combo
is untested. LOUD crash (ValueError), NOT silent-wrong → no false-green ships, but the
closeout's "2-D Krylov works end-to-end" is only true for no-pure-z quads. Fix = port
the sweep `pure_z` moment-broadcast guard to the matvec `pure_z` + add a committed 2-D
LD Krylov≡SI gate on a quad WITH pure-z (the L-018/L-021 "matvec needs a committed
call-count/end-to-end gate, not a round-trip" lesson, recurring a THIRD time).

OTHER follow-ups (non-blocking): (a) 4 net-new pyright nits (apples-to-apples
stash PRE 110 / POST 113 + comm: `D1ClosedForm` un-imported but used in a
`moment_scan_closure` return annotation — runtime-safe via `from __future__`;
`scheme.moment_scan_closure` LD-only method on a base-typed handle behind the
`is_moment` gate; `Q.spatial_moments_per_axis` narrowness; 1 reportReturnType) —
all Pattern-3 type-narrowness debt; (b) `verifies("ld-cartesian-1d","ld-slab")`
labels have NO `:label:` math block in docs/theory (pre-existing — already at HEAD,
audit tracks them with 6/4 tests, exit 0; `ld-cartesian-2d` correctly deferred per
L-031); minted in S4/archivist. Mode-8 clean (2 prod bare asserts in touched files
PRE-EXISTING: loss_rep:2376 type-narrow, solver:866 in `if __debug__:`). VERDICT
2026-06-17: SUPPORTED-WITH-CONCERNS. Numerics + the headline fix SOUND, byte-id real,
markers honest. The d≥2 matvec pure_z crash is a real latent defect (loud, narrow)
+ the missing d≥2 Krylov gate is a genuine coverage hole — both follow-ups, NOT
commit blockers (the SHIPPED scope is SI-verified d≥2 + matvec-verified d=1; no
false-green). No false-green found anywhere.

## L-034 -- #240 D5b-S4 (2-D LD stress MMS) VERDICT SUPPORTED-WITH-CONCERNS: the honest-scope claim must be SHARPENED -- the scattering channel EXERCISES the slope-source code path but the MMS is empirically BLIND to a sign error there

D5b-S4 = a vv-Mode-7 strengthened 2-D Cartesian LD MMS `ψ=[A+μ_x·B+μ_y·C]/W`
verifying the multi-D bilinear UBLD slope rows landed in S3. VERDICT SUPPORTED-WITH-
CONCERNS; numerics SOUND, no false-green, no blocker.

1. **L11 structural independence GENUINE.** The SymPy source is from the
   CONTINUOUS PDE (no `_LDCellTerms`/`_schur_terms`/`_ubld`). The FD-residual
   cross-check IS a genuine 2nd structural path PROVEN: corrupt `Q_closed`'s
   `μ_x·∂_xA` streaming sign → FD residual 0.047 ≫ 1e-7 tol (FD uses numpy
   central-diff of ψ, RHS embeds SymPy `diff` → a diff sign error breaks the
   equality). Branch2==Branch1 source ≤1e-13. Single-source `_LD2D_STRESS_COEFFS`
   `(num,den)` pairs (Rational∥float) = amplitudes can't drift.
2. **⭐ THE HONEST-SCOPE FINDING (sharpen point 2).** The closeout/docs say the
   slope-SOURCE half is "DEFERRED" because the EXTERNAL Q̂ is zeroed
   (`_lift_external_source_to_moments` → `lifted[...,AVERAGE_MOMENT]`, slopes=0,
   confirmed). BUT the SCATTERING source `Σ_s·φ̂` IS a genuine `(N,ng,nx,2^d)`
   moment source consumed through the SAME `Q_cells` slot as an external Q̂ would
   be (loss_rep:2814-2825 lifts BOTH into the same `Q_per_ord`→`QV_per_ord`; the
   slope-row sign code path `_reframe`+UBLD is source-AGNOSTIC). INSTRUMENTED the
   solve: iterate scalar-flux moments fed to `apply_p0_in_scatter` carry NON-ZERO
   slope rows (avg=1.31, x-slope=0.257, y-slope=0.129, xy=0.067), scattered
   `Σ_s⊗I_spatial` (`fg,fc...->gc...`) → the slope-source rows ARE populated +
   consumed. SO the slope-source CODE PATH is exercised. BUT — DECISIVE MUTATIONS
   on the slope-source rows (`_CellSolve.cell` Q_cells[...,1:]): SIGN-FLIP → order
   stays 1.97, finest in-band → NOT caught; ZERO the rows → order 1.99 NOT caught;
   ×3 magnitude → order 2.02 NOT caught. The scattering-slope source is an
   O(h)-small DG-internal forcing (slopes ~5× < average, c≤1.0) whose sign/
   magnitude affects the converged flux ABOVE O(h²) → absorbed in the floor. So
   the PRECISE honest claim is NEITHER of the brief's two options: it is
   "slope-source code path EXERCISED via scattering but the MMS is BLIND to a
   slope-source SIGN error (a sign flip is not caught) — genuinely UNVERIFIED for
   the sign convention; external-Q̂ plumbing also deferred." The docs' "DEFERRED"
   is substantively CORRECT (sign unverified) but the parenthetical "the only
   moment-valued source consumed is Σ_s·φ̂" UNDER-states that this consumed source
   does NOT verify the sign → CR3 doc-sharpen (follow-up, not blocker): the note
   should say the scattering channel exercises-but-does-not-constrain the
   slope-source sign.
3. **VALUE band REAL + tight (Mode-5 not rate-only).** Reproduced: errs
   [1.42e-2,3.54e-3,8.81e-4], orders [2.00,2.01], maxrelerr 1.78e-2→1.2e-3 (4×/
   halving), flux range matches ref. Band (1e-9,1e-2): upper 1e-2 is BELOW the
   coarsest error 1.42e-2 → a wrong-limit/non-converged flux WOULD exit the band.
   Genuine value gate.
4. **Mutation conclusion SOUND + I isolated what the closeout couldn't.** The
   slope-UNKNOWN half: sign-flip `_GRAD_1D[1,0]:-2→+2` → NaN (caught); a FINITE
   x↔y-symmetric 10% error (`_GRAD_1D` ×0.9) → order −0.06, finest 0.072 ≫ band
   (CAUGHT, NOT divergent — the missing subtle-finite discriminator). So the
   strengthening is non-vacuous + load-bearing for the slope-UNKNOWN half (catches
   both catastrophic AND finite). The strengthening's SPECIFIC x↔y-asymmetry value
   targets the slope-SOURCE same-sign trap — which (per finding 2) this MMS cannot
   reach AT ALL → the x↔y strengthening is defensive-correct-but-currently-
   untestable (ship per spec; its payoff arrives with the moment-source increment).
5. **Gate/marker integrity.** `ld-cartesian-2d` minted (1 unique `:label:`),
   audit `ld-cartesian-2d → 4 tests` exit 0, all verifies targets resolve
   (transport-cartesian-2d/multigroup/mg-balance exist). Quadrature exactness
   CONFIRMED: LS S4 `<μ_x²>=<μ_y²>=1/3`, `<μ_xμ_y>=<μ_x>=<μ_x³>=0`, ZERO pure-z →
   `φ=A` exact. Mode-8 clean (0 bare asserts new files + 0 in S4 prod additions).
   Mode-7 declaration present. ⚠ L-007 NIT: `test_v_ld2d_stress_substitution_
   identity` stacks `@foundation @verifies("ld-cartesian-2d")` = the conflation
   L-007 warns of — BUT (a) established convention (anisotropic_symbolic.py:69-70/
   147-148 do the same for algebra-of-record substitution gates), (b) the label is
   NOT solely foundation-backed (3 genuine L1 verifiers) → the L-007 tell ("ONLY
   coverage is foundation") does NOT bite → minor consistency nit, not blocker.
6. **Modes 1-6 defenses ALL present.** 2G, het (176 unique cell materials,
   spatially-varying σ_t), STRICTLY-asymmetric downscatter (SigS[0→1]=0.233,
   [1→0]=0.000 — pure downscatter, transpose-sensitive), non-square mesh 16×11 +
   domain 1.3×0.9.

VERDICT 2026-06-17 SUPPORTED-WITH-CONCERNS: 1 CR3 doc-sharpen (the scattering
channel exercises-but-doesn't-constrain the slope-source sign — the honest note
under-states this) + 1 L-007 marker nit (convention, non-biting) — both
follow-ups. The shipped scope (slope-UNKNOWN sign verified + average-moment
boundary + matvec twin + two-paths, mutation-verified non-vacuous) is honest. No
false-green, no blocker.

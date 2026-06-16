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

## L-004 -- vv-status rationale comments must NOT use [brackets]

The `:vv-status: documented` directive lives in the same RST file as
the labelled equation, conventionally as a top-level RST comment
(`.. vv-status: <label> documented`). When attaching a rationale
comment block with `..    [category] description` formatting,
docutils parses each [xxx] as a citation reference, producing
"duplicate citation" warnings under `-W`. Use (parens) instead of
[brackets] in rationale comments, and prefer a single-line `..`
comment over a multi-line `..  / .. / ..` continuation block.

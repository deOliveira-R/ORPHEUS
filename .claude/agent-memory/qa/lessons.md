# QA Lessons

Behavioral corrections only. AGENT.md has the V&V hierarchy,
anti-patterns, and error catalog format -- never duplicate here.

> **Promoted to AGENT.md (2026-06-22):** the standing stance "a green
> gate is evidence of nothing until you have made it RED — mutation-verify
> every gate's teeth under `-O`, in-process, revert by re-editing" is now
> Enforcement #11 in `.claude/agents/qa/AGENT.md`. The lessons below keep
> the per-incident *mechanics* (which mutation point, which sentinel, which
> revert proof) — those stay here as recalled technique. The *rule* lives in
> the definition. Recurring instances: L-007, L-014, L-020, L-024, L-027,
> L-031, L-033, L-036, L-039, L-040, L-042, L-045..L-050.

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

## L-004 -- vv-status rationale comments must NOT use [brackets]

The `:vv-status: documented` directive lives in the same RST file as
the labelled equation, conventionally as a top-level RST comment
(`.. vv-status: <label> documented`). When attaching a rationale
comment block with `..    [category] description` formatting,
docutils parses each [xxx] as a citation reference, producing
"duplicate citation" warnings under `-W`. Use (parens) instead of
[brackets] in rationale comments, and prefer a single-line `..`
comment over a multi-line `..  / .. / ..` continuation block.

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

## L-016 -- "branch fires under quadrature Q" claims need a degeneracy probe, not faith

#206 claim 5 asserted the cylinder pure-azimuthal degenerate-ordinate branch
(`|mu_x|<1e-15`, `A_downstream=0`) "fires under the A-NEW matvec[CYL] leg".
FALSE: `Quadrature.level_symmetric(sn_order=2..8)` ALL have `min|mu_x| ≥
0.22` — ZERO degenerate ordinates. The branch is dead code under standard
LS cubature and is exercised by NO current test. Probe in 3 lines
(`np.count_nonzero(np.abs(q.mu_x)<1e-15)`) before accepting a branch-
coverage claim. (Not a regression — the branch was a verbatim relocation,
proven by L-013 — but the EVIDENCE for it is vacuous; flag the coverage gap.)

---

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

---

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

---

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

---

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

---

## L-021 -- Increment B closed the L-018 matvec-coverage gap for LD

L-018 flagged (Inc A) that LD's matvec was correct-but-UNVERIFIED (no committed
test drove `residual_kernel_batch`; SI sweeps only touch the solve kernel).
Inc B added `test_sn_1d_slab_ld_mms_krylov_matches_si` (inner_solver='krylov'),
which an instrumented call-count proves drives `residual_kernel_batch=640`,
`cell_kernel_batch=0` -- the matvec path is now genuinely exercised end-to-end
AND pinned ≡ the SI sweep. When re-reviewing a follow-up increment, re-check
whether it closes a prior increment's NEEDS-EVIDENCE item (the call-count probe
is the verification, not the round-trip self-consistency test).

---

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

---

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

---

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

---

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

---

## L-026 -- "scattering exercises the slope-source path" does NOT mean the MMS constrains its SIGN (Mode-10 honest-scope)

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

---

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

---

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

---

## L-029 -- a re-encoded closed form is NOT L11-circular when compared against an INDEPENDENTLY-ASSEMBLED primitive (the d=1-reduction-to-production oracle)

The circularity principle: a token-for-token copy of a production formula is
circular as a VALUE check (a sign-flip in prod propagates into the test, stays
green). The DISTINCT case here (the SymPy UBLD): the test's RIGHT side re-encodes the
production `_schur_terms` S/eff_source/.slope (test lines 346-351, verbatim of
`linear_discontinuous.py:332-335,258-259`) — BUT the LEFT side is the symbolic
primitive's d=1 reduction obtained by `A⁻¹R` of a SEPARATELY-built Kronecker
matrix (`assemble_ubld([h],[mu],...)` → `LUsolve`), NOT a re-statement of
`_schur_terms`. So `diff_psi_bar==0` proves "the production Schur scalar EQUALS
the independently-assembled 2×2 solve" — one side is genuinely structurally
independent. A sign-flip in `_schur_terms` would NOT propagate into the Kronecker
assembly → the oracle WOULD catch it. The circularity test is: **does the bug
live on BOTH sides of the diff?** Re-encoded-formula-vs-re-encoded-formula =
circular; re-encoded-formula-vs-independent-construction = legitimate.

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

---

## L-030 -- Intrinsic-property gates: a PER-CELL invariant tested only with a SPATIALLY-UNIFORM fixture is untested for its per-cell-ness

`SpectrumField.__post_init__` enforces the simplex PER CELL (`values.sum(axis=0)` then
`np.allclose(col_sums,1)`), but every test fixture (`_chi` helper) builds a per-cell-UNIFORM
χ — so "per-cell sum==1" and "global mean==1" are INDISTINGUISHABLE in the suite. The code
is correct (probed: a χ summing to 1 in cell0 / 1.2 in cell1 IS rejected) but the per-cell
distinction is uncovered. **Rule**: when a validator's invariant is keyed on a specific AXIS
(per-cell, per-group, per-ordinate), at least one negative fixture MUST VARY along the
non-reduced axis so a global-collapse mis-implementation (`values.sum()` vs `.sum(axis=0)`)
would be caught. Uniform fixtures are blind to axis-collapse bugs — the spatial analogue of
the 1-group degeneracy (H1).

**Validator Mode-11 recipe for a `__post_init__`/`replace`-revalidating type** (verified live
on #257): (a) the `+`-routes-through-revalidation claim — a leaf that does NOT mix in `FluxRole`
inherits `Field.__add__` → `replace()` → re-runs `__post_init__`; PROVE by tracing
`__post_init__` call-count during `χ+χ` (fires ≥1, then raises simplex). (b) a negative test
that pins ONE branch in isolation (neg-entry with col_sum==1.0 exactly) genuinely isolates it
ONLY IF the OTHER branch (sum) would PASS that fixture AND runs SECOND — read `__post_init__`
ordering (neg-check before sum-check) + confirm the raised message names the intended branch,
not the masking one. (c) `mix`/convex-blend re-validation rides the same `replace` path.

**Doc-attribution drift (CR3, no phantom edge)**: `spectrum_field.py:32` frames a χ-drift /
depletion-feedback bug as "the ERR-039 normalization-bug class" — but ERR-039 is the moment-
projection `apply_transpose` (2ℓ+1) addition-theorem factor (nothing to do with χ or depletion).
`units.py:100` "ERR-039 normalization class" for a missing-`/4π` is a looser-but-OK framing
(ERR-039 IS an angular-norm factor). No `catches("ERR-039")` marker exists on these tests, so
no phantom coverage edge is written (contrast L-007) — pure stale-doc, flag don't block.
RULE: a prose "the ERR-NNN class" citation in a NEW docstring still warrants a 30-sec catalog
read; mis-citation in prose is a nit, but the same string in a `catches()` marker is a defect.

---

## L-031 -- a self-consistency round-trip + an A==A pin can BOTH be blind to the bug their docstring/marker claims; the genuine catcher is the continuous-PDE exact-on-bilinear oracle; and "matvec twin verified" ≠ end-to-end Krylov

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
   `test_d2_exact_on_bilinear` (the L11/L-029-clean continuous-PDE oracle: source
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

---

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

6. **Pyright net-new = 0** proven apples-to-apples (L-027): run the SAME
   5 touched files PRE (stash prod + hold the new untracked file) and POST,
   path+line-strip, `comm -23`. All 8 errors pre-existing (`DualSpace.of`
   return, `MaterialXSField` Optional, `from_face_arrays` Optional layout,
   `_check_partner` `other.L` on object); the NEW file alone = 0 errors; the new
   `find_factor` RAISES KeyError (not returns None) → type-clean. The brief's
   worry about "find_factor-returns-object at space.py:521" was a brief
   mis-attribution: :521 is `DualSpace.of`, pre-existing. (VERDICT 2026-06-16:
   SUPPORTED, no blocker, no follow-up.)

---

## L-033 -- a GENUINE structural-independence ground (L11 clean) can still ship a LATENT twin-path crash no committed gate drives (the d≥2 matvec pure_z gap)

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

---

## L-034 -- Stale-snapshot triage: a HUGE-ULP bit-identity red is "live correct, frozen stale" until you find the apply-changing commit that did NOT re-capture

A frozen-snapshot bit-identity gate failing with a catastrophic metric
("not equal to 5 ULP, max is 8.8e15"; or `assert_allclose` 100% mismatched
few-%) on ONE geometry arm while sibling arms (slab/cylinder) PASS is the
fingerprint of **a stale snapshot left by an unrelated correctness fix**, NOT a
live solver bug. The live apply is usually the MORE-correct value; the frozen
reference is stale. Triage procedure (do NOT modify anything — produce verdict):

1. **Confirm the asymmetry** — run the failing arm AND a sibling arm (other
   geometry/seed). Sibling green + this red = geometry-scoped change, not a
   broad regression. (Here SLB matvec passes @1 ULP DriftWarning; cart2d+cyl
   streaming arms pass; only SPH fails → sphere-only-by-construction fix.)
2. **Blob-hash the fixture across refs** — `git rev-parse <ref>:<fixture>`.
   If the snapshot blob is byte-IDENTICAL since the last refresh commit but the
   code moved, the snapshot is the stale side. (SPH `.npy` `501fd29` unchanged
   since the ERR-058 refresh `798372f`; the `.npz` later changed at #240 but
   only its CARTESIAN arrays — sphere arrays rode the stale `798372f` values.)
3. **Find the diverging apply commit** — `git log <refresh>..<base> -- <prod
   apply files>` and grep the messages for the geometry + math-term keyword
   (curvilinear/sphere/clamp/seed/closure/weight). The culprit is a commit that
   (a) changed the geometry's apply value, single-sourced so the matvec
   inherits it, AND (b) touched NONE of the failing suites' fixtures (staleness
   slips in SILENTLY because the commit's own `-O` sweep didn't run them). Here:
   `b2d8a6d` "unclamp spherical Morel–Montry weight (Bailey Eq. 43)" — dropped a
   spurious `[½,1]` τ-clamp in `spherical_streaming`, regenerated only its OWN
   targeted snapshot, left these 2 suites stale.
4. **Verify the call path** — Nexus `callers(spherical_streaming)` →
   `SNMesh._init_core` → `StreamingOperator.apply`: BOTH the matvec test
   (`_LpC_apply`) and the streaming-operator test (`L.apply`) consume the same
   producer, so both inherit the change.
5. **Does the unmerged sibling branch fix it?** Git-archaeology beats a worktree
   run (editable `.venv` resolves to MAIN tree; worktree creation may be denied
   anyway). `merge-base --is-ancestor <fix> <branch>` (does it contain the fix?)
   + blob-hash the fixture on the branch vs base. If fixtures are byte-identical
   AND the branch changes the apply path FURTHER, the branch does NOT fix it —
   it inherits the same stale snapshot and moves the live value further away.
   (#236 forks off the clamp fix, leaves both fixtures byte-id, reworks
   `pole_angular_closure.py` +942 → would still fail, possibly differently.)

**Verdict mapping**: stale snapshot + NO open issue owning it = (B)
NEW/UNTRACKED — recommend a `tests(sn)`+`type:bug` re-baseline issue, NOT a
correctness regression. The fix = re-capture on main, validated against the
STRUCTURALLY-INDEPENDENT grounds (the matvec file's `Q/Σ_t` L0 streaming-
equilibrium row + the curvilinear L1 MMS/closed-form k_∞ the streaming class
cites) per vv bit-id criterion 2 — NOT just old-vs-new ULP. ⚠ A test docstring's
issue attribution can be STALE/wrong: this SPH red's docstring cited "#195/#209"
(both CLOSED, both different mechanisms — ERR-058 MMS-rate + cylindrical-pole
NaN); the REAL cause `b2d8a6d` is Refs #229. Trust the git-archaeology over the
docstring's cited issue number.

---

## L-035 -- "byte-identical EXCEPT a LATENT collision" claim: instrument BOTH branches on the gate suite + classify each divergent site as reached/correct

When a refactor claims "byte-identical to all paths EXCEPT one latent S4-style
collision," do NOT trust the latency claim — PROBE it. Patch the touched primitive
to compute BOTH the OLD and NEW branch on every call and assert array_equal,
running the FULL gate suite (a pytest plugin `pytest_configure` reassigns
`module._symbol` in EVERY importing module — `sg._reframe` AND `lr._reframe`;
attribute via `pytest_runtest_setup` → `item.nodeid`; pytest captures stderr so
read the probe under `-s`). #246's `_reframe` keyed on `is_moment_valued` (typed
origin) vs OLD `arr.shape[-1] != frame_signs.shape[0]` (size probe) — claimed the
d=2 `2^d==4` collision was "latent in production." The probe found **48
divergences, NEW≠OLD by 70%** (NOT FP, NOT zero), all at `_CellSolve.cell`'s
`Q_cells` reframe, in exactly TWO tests: `test_ld_2d_two_paths_ffw_equals_mfw` +
its `_stress_` sibling. So the collision is REACHED — not latent — by the
LOW-LEVEL `MovingFrontierWindow/FullFieldWavefront.sweep(Q_flat, ...)` API with a
flat source whose anti-diagonal level has exactly `2^d` cells (`n_diag==4`).

The decisive correctness call (which branch is RIGHT) needs a STRUCTURALLY-INDEPENDENT
reference, NOT the gate (the two-paths gate is Mode-11 BLIND: both FFW+MFW legs
share `_CellSolve`, both corrupted identically under OLD → agree to 0.0 EITHER way,
so it cannot distinguish OLD-wrong from NEW-correct). The independent reference:
moment-LIFT the same flat source onto slot 0 (`face_moment_tail`, slopes=0) and
sweep — a flat source and its zero-slope lift are the SAME physics. Result: NEW
flat-sweep == moment-lifted-sweep BYTE-IDENTICAL (0.0); OLD flat-sweep ≠ lifted by
70% (OLD was inconsistent with its OWN lift — the size-probe mis-classified the
4-cell anti-diagonal as a 4-slot moment vector and applied a spurious `[1,-1,-1,1]`
involution, scrambling cells). ⟹ NEW is CORRECT, OLD was a REACHED (test-path)
silent error.

⚠ The "production" qualifier is load-bearing: `solve_sn_fixed_source` ALWAYS
moment-lifts via `_lift_external_source_to_moments` (source reaches `_CellSolve`
as ndim+1 → rank test True → never the flat collision), so production USERS never
hit it — only the low-level `.sweep()` test API does. So "latent in production
(via the public solver)" is TRUE; "latent everywhere" is FALSE. Verdict =
SUPPORTED (the fix is correct and strictly better than OLD), but flag the framing:
it FIXES a reached test-path silent error, it is not purely prophylactic. The
rank discriminator `Q.ndim > sigt.ndim + 1` is genuinely S4-safe (a rank cannot
collide the way a trailing-size can) AND correctly classifies BOTH entry points
(flat `(N_oct,ng,n_diag)` ndim-3 → False; moment-lifted `(N_oct,ng,n_diag,2^d)`
ndim-4 → True), single-source-shared with `_moment_broadcast_sigma:515`. Gate-1
(`test_reframe_moment_intent.py`) is the genuine unit catcher (mutation-verified:
emulate OLD `_reframe` → Row-1 `out==arr` FAILS, sign-flips `[0,-1,2,-3,...]`);
Gate-3 `is_multi_moment` mutation-verified (const-True reddens DD-P2, const-False
reddens LD-P1). DD byte-id confirmed (regression suite: 0 divergences, all
short-circuit on `frame_signs is None`; the 13 within-tol DriftWarnings PRE-EXIST,
not escalated under `-W error`). 1-D scan + LD-kernel suites: 0 divergences (those
sites only ever reframe genuine moment arrays where size-probe≡intent agree).

---

## L-036 -- "MMS covers the retired term" claim: mutation-verify, but the deleted code may be a SWEEP↔MATVEC twin (mutate the SHARED coefficient source, not the deleted-apply method)

A retirement that deletes "verification machinery" (here the `M_spatial`/`M_angular_redist`
separately-applicable operator-leaf split + `loss_action_decomposed` + the `emit_angular`
arm of `_apply_walk`) defends "no correctness lost" by naming a surviving MMS. The
deleted decomposition tests pinned only STRUCTURE — `TestT4c…` asserts
`(m_full − m_ang) + m_ang == m_full`, a **TAUTOLOGY by the subtraction construction**
(`m_spat = m_cell − m_ang_cell` is literally the deleted code); `TestT4b…` are
`isinstance`/`direction_sign`/`capabilities`/`cached_property` type pins +
`L == M_spatial − σ_t·ψ` (self-referential, both sides one walk); `TestT5…` is
`from_geometry == from_geometry`. None pins an INDEPENDENT correctness invariant of
the production `m_full`. The deleted `m_ang` emission wrote ONLY a separate buffer no
production code read — the fused `m_full = (denom·ψ̄ − numer_upstream)/V` already
carries the redistribution (it lives INSIDE `denom`/`numer_upstream` via the closure's
`(ΔA/w)·c_out` / `(ΔA/w)·c_in·ψ`), so the production output is byte-identical (strict DD
regression passed at the SAME documented 6920-ULP baseline, 0 `.npy`/`.npz` regenerated,
607-del/48-ins).

⭐ THE TRAP when mutation-verifying the surviving MMS: the deleted code was in the
MATVEC/apply path (`_apply_walk`), but the MMS runs `solve_sn_fixed_source` = the
SWEEP/solve path. These are genuine TWINS that share only the precomputed coefficients.
My first 3 mutations (of `MorelMontryAngularSweep.cell_contribution`, then
`_redistribution_for_level`/`__call__`, then `cell_balance_terms`) ALL showed
call-count 0 on the sphere solve and byte-identical error ladders (GREEN-BLIND =
patching dead-for-this-path code). The sphere sweep is the `ScanMarch`/`CumprodScan`
1-D scan in `loss_representation.py:3106` (`ang_contrib = dA_w·c_in·ψ` into the source
`b` + `c_out` baked into `inverse_denom`); the cylinder routes through
`diamond.update`→`cell_balance_terms`. THE RIGHT MUTATION POINT = the SINGLE SHARED
source `GeometryCoefficients.from_mesh_and_quad` (`c_out=α_out/τ`, `c_in=(1−τ)/τ·α_out+α_in`,
sweep_cache.py:309-310) — `dataclasses.replace(gc, c_in=f·gc.c_in, c_out=f·gc.c_out)`.
Confirm the factory call-count > 0 + identity-reimpl (f=1) reproduces baseline FIRST,
THEN mutate. Result: c_in/c_out sign-flip, ×3, even ×1.5 → BOTH sphere AND cylinder MMS
go NaN or land orders of magnitude outside the gate bands (1e-8..5e-3 / 1e-3..5e-2).
The redistribution is an O(1) term in the curvilinear cell balance — the MMS STRONGLY
constrains it end-to-end. Verdict: SUPPORTED, no correctness coverage lost.

PROCEDURE for "deleted machinery covered by surviving test X": (1) read each deleted
test → classify STRUCTURE-pin vs CORRECTNESS-pin (tautology-by-construction = structure);
(2) confirm the production value is byte-identical (strict gate at documented baseline +
0 snapshot regen + diff is deletion-dominated); (3) mutation-verify X catches the term —
but FIRST instrument a call-count to find which of the SWEEP/MATVEC twins X exercises and
mutate the SHARED coefficient source, not the deleted-apply method (else GREEN-BLIND on
dead-for-this-path code mis-reads as "X is blind"). Mode-8 caveat: the MMS uses bare
`assert` under `-O`, but pytest's rewriter fires asserts in `tests/` modules (L-010) —
proved live by breaking a band to red. Baseline reds (5 stale SPH snapshots #250 + 2 mu_y
#232) are pre-existing: cylinder snapshot siblings walk the SAME modified curvilinear
`_apply_walk` and PASS; only SPH fails at ~1e15 ULP (stale-snapshot signature, L-034),
not the 1-ULP FP drift a real regression would show.

---

## L-037 -- Mode-10 closeout verification recipe (the activated-but-unconstrained slope source)

#247 Leg A closed the slope-SOURCE half of the LM-1989 trap for 2-D Cartesian LD
(the vv Mode-10 gap: a term genuinely CONSUMED yet a sign flip leaves the
converged flux at O(h²) + ~1.4×, sub-floor). VERDICT SUPPORTED, NO ERR (the
slope source was UNVERIFIED, not WRONG — the production lift correctly zeroed an
honest default q̂=0). Reusable recipe for adjudicating a Mode-10 closeout:

1. **The teeth are NOT the converged-flux value-band** (the §0 trap — the slope
   error is O(h²)-small). Demand TWO O(1) structural teeth instead: (a) the
   PRODUCER threads the projection through at machine precision
   (`array_equal(lifted, Qm)` — the production-change proof; under the bug, the
   dropped slope is O(1) e.g. 0.179); (b) a CONSUMED source-row sign flip moves
   the converged answer ≫ solver tol (the consumption proof).

2. **Prove the teeth bite by re-introducing the EXACT bug in-process** (throwaway
   conftest plugin monkeypatching the producer to re-zero slopes — NO production
   edit, L28). The sign-mutation gate's red message is the tell:
   `|Δφ|/|φ|=0.000e+00 ≤ tol — the slope row is NOT consumed` (flipping a
   re-zeroed row is a no-op). Confirm the 3 new gates RED.

3. **Prove the Mode-10 ASYMMETRY**: run the EXISTING flat scalar gate under the
   SAME buggy producer → it MUST stay GREEN (it feeds a flat source → slope row
   already zero → blind to the slope sign). GREEN-flat + RED-moment IS the gap.

4. **Calibrate the consumption-tol live**: a deterministic SI solve has noise
   floor EXACTLY 0.0 on an identical re-solve (measure it). Smallest signal
   (xy slot) ~5.8e-5 clears 1e-8 by ~5.8e3× → no false-green. The tol is
   defensible iff (noise ≪ tol ≪ smallest signal).

5. **L11 check the projector source**: `leggauss` + numpy + hand-laid algebra
   ONLY. Typed source CONTAINERS (`AngularSourceSink`/`TimedFullField`) imported
   to FEED the solve do NOT contaminate the reference — the LD cell op / the lift
   must not be called. The foundation sub-gate's reference is hand-derived
   polynomial coefficients, not a production echo.

6. **Confirm no latent CONSUMPTION bug** (else mint ERR): the now-consumed slope
   path must have no sign/magnitude error. The architecture is the proof when
   the producer change RIDES an EXISTING consumer path (external + scattering
   moment vectors SUMMED into ONE global-frame array → shared rank-gated
   involution reframe `octant_moment_frame_signs`, shared mass M=diag(h,θh),
   shared Kronecker order). Dispatch explorer to trace reframe/mass/Kronecker;
   if no separate external-vs-scattering branch (no extra/missing flip, no
   transpose), the consumed path is correct.

WATCH (non-blocking doc-nits this review surfaced): a "single-source shared by
fixed-source AND eigenvalue" docstring can be STALE (grep the lift's callers —
#247's lift has ONE prod caller; the eq path wraps its sweep OUTPUT, doesn't
call the lift); a d=2 cell-unknown prose label can transpose x↔y vs the
canonical [bar,y,x,xy] (axis0=x OUTER) — prose only, slots come from
moment_layout, no code path. Both CR3/stale-doc (L-020/L-028), no ERR.

---

## L-038 -- prove an xfail→live FLIP is non-vacuous TWO ways + Mode-10 with NO dominant regime (boundary transverse face-slope)

#251 widened the 2-D Cartesian LD boundary trace to CARRY the `2^{d-1}` transverse
face-slope (boundary twin of #247 Leg A's bulk slope). The boundary slot grows a
trailing moment axis at ONE lever (`geometry.boundary_face_layout` appends
`face_moment_tail(per_axis**(ndim-1))`); `_inflow_to_moments` rank-discriminates
(`is_moment_valued_by_flat_rank(face, mesh.ndim+1)`) → scalar arm seeds slot-0 only,
moment arm passes through; 4 outflow capture-collapse sites dropped so the outflow
moments land in the now-moment-shaped slot. VERDICT SUPPORTED, NO ERR, NO blocker.

1. **An xfail→live FLIP needs TWO red-proofs, not one.** A gate the closeout says
   "was xfail-strict, now passes via production" is only non-vacuous if it (a) goes
   RED against a re-introduced post-change bug AND (b) goes RED against the EMULATED
   PRE-change behavior. The re-zero mutation (`f[...,1:]=0` in the NEW moment arm)
   proves the consumed slope is constrained; but EMULATING the old unconditional
   zero-fill (treat the `(...,2)` moment face as scalar → spurious `(...,2,2)` axis +
   slot-1 zero) is what proves the gate genuinely REQUIRES the #251 change (threading
   gate red "did not RECOGNISE the moment-resolved inflow"; width-reject red "DID NOT
   RAISE"). Only (b) rules out a gate that was already green at HEAD. Do BOTH via a
   throwaway conftest plugin under `-O` (L28: no prod edit, no git stash).
2. **Mode-11 closure for a public-solve gate = INSTRUMENT the rewired arm + count.**
   When a closeout says "the surrogate monkeypatch was dropped and the gate re-targeted
   onto the public API," confirm the rewired production line is on the LIVE call graph:
   monkeypatch the method to COUNT which arm fires during the public solve. #251:
   `_inflow_to_moments` fired 344×, moment-resolved arm 688× (0 scalar/identity) on the
   public `+slope` `solve_sn_fixed_source` → Mode-11 CLOSED (the gate drives production,
   not a recompute-on-both-sides surrogate). The consumption gate's RED under the re-zero
   mutation is THROUGH that public path (|Δφ|/|φ|=0.000e+00 — flipping a zeroed slope is
   a no-op = the exact Mode-10 signature).
3. **Mode-10 with NO O(1)-dominant regime → structural teeth are the COMPLETE
   resolution (no value-improvement leg).** A boundary-trace slope is sub-floor for ANY
   value claim, not just the sign (probed: seeding the REAL slope makes near-bdy A-err
   WORSE 2.131e-2→2.163e-2, flipped is BETTER). So "improves-on-flat" is UNACHIEVABLE
   and dropping it is HONEST, not hiding a problem (keeping it would falsely RED a
   correctly-consumed slope). The companion-gate half of the Mode-10 recipe (isolate
   the term so its error is O(1)) is UNAVAILABLE — no fixed-source problem makes a
   boundary-trace slope the dominant forcing. Positive verification = TWO O(1)
   structural signals ONLY: machine-precision threading (`array_equal` slot-1, leggauss
   reference = L11, NON-circular b/c prod's arm is a pure pass-through) + consumed-flip
   ≫ TOL (4.101e-3, triple-agrees across my re-derivation / test-architect surrogate /
   public-path, above the deterministic 0.0 noise floor). This is a NEW Mode-10
   sub-case (neither #240 D5b-S4 nor #247 Leg A had a term with no dominant regime —
   both could improve-on-flat) → warrants the test-architect's one-line vv Mode-10 row
   addition.
4. **Reflective storage pass-through ≠ reflective SIGN.** A trace-widening's reflective
   path has TWO concerns: storage (the `PermutationOperator(axis=0)` broadcasts the new
   moment axis — verify NO corruption by seeding a random moment-shaped trace, running
   `_reflect_trace`, and checking slot-1 follows slot-0's permutation: #251 = 0
   corruption over 12 matched ordinates) and SIGN (the transverse-slope sign under a
   normal-flip reflection — UNVERIFIED b/c the vacuum-BC MMS nulls the reflective
   coupling, H2). Storage-correct is shippable; the SIGN is a Mode-1 trap the vacuum
   gates CANNOT see → MUST be a follow-up (#252, filed OPEN with correct labels), NOT a
   blocker. Confirm the follow-up issue actually exists (`gh issue view`) before
   crediting "filed as #NNN."
5. **Producer-rank carve: a widened slot needs the EXISTING SCALAR producers audited,
   not just the new moment one.** When a carve widens a trace/field slot, the existing
   scalar MMS callers feed the SAME widened slot → the producer (`prescribed_inflow`)
   must accept BOTH ranks (scalar→seed slot-0; moment→write full slot). The explorer's
   "1 real producer edit" under-scoped this by 1 (the scalar-onto-moment arm). Same
   class as Leg A's field-space layer: a rigid scalar contract above a widened slot
   needs a typed-union relaxation, not just an indexing fix.

⚠ Minor scope-note (NOT a defect): spec D6 said "DD rejects any moment trace" but
the impl early-returns IDENTITY at `n==1` (DD never receives a moment inflow —
`face_moment_tail(1)==()` makes the DD trace scalar-only), so the shipped reject-gate
only tests LD wrong-width. Correct by construction; flag for the spec author only.

---

## L-039 -- runtime_checkable Protocol gate: prove teeth by DROPPING/adding a member, not by a vacuous isinstance pass

Validating a `@runtime_checkable` structural Protocol as a REAL coverage gate
(vs a vacuous isinstance pass). VERDICT REAL-GATE. The recipe, consolidated
from two reviews (a `Vector` carrier Protocol and a `TransportState` state
Protocol):

(1) **Mode-8 is cleared by the mutation, not by inspection** -- bare asserts
in a COLLECTED test module ARE rewritten by pytest and DO fire under `-O`
(the `PytestConfigWarning` refers to NON-test modules; reaffirms L-010).
Proof: a Protocol mutation made the asserts raise a real `AssertionError`
under `-O`; `pytest.fail` likewise raises `Failed` under `-O` (a print after
it never runs).

(2) **Mutate the PRODUCTION Protocol, not the test object.** Two complementary
moves:
   - **Drop a required member** -- DROP only `__rmul__` -> ONLY
     `test_no_scalar_mul_rejected` reds while `test_string_is_not_vector` STAYS
     GREEN (str still lacks `__sub__`); the asymmetry proves each negative test
     OWNS a specific dunder. Reduce the Protocol to `__add__`-only -> ALL
     negatives red, positives stay green.
   - **Build an in-memory MUTANT Protocol and monkeypatch it in** -- a
     drop-all `class M(Vector, Protocol): ...` flips `np.ndarray` AND the leaf
     type to `isinstance==True`; run the REAL test fn with the production
     name monkeypatched to the mutant (patch BOTH the defining module's name
     AND the test module's import binding) -> the discriminating negative
     (`test_ndarray_..._not_a_transport_state`) fires with the exact message.
   The RED message names the wrongly-accepted object. Revert by RE-EDITING
   (untracked `??` files make `git diff` empty -> vacuous; the real revert
   proof is gate-green-again + grep zero MUTATION markers + dunder/line-count
   match).

(3) **`scalar * vector` is `__rmul__` NOT `__mul__` -- prove with the
Python micro-fact, not the docstring** (L-020 discipline): `0.0*OnlyRmul()`
fires `__rmul__`; `0.0*OnlyMul()` raises TypeError (`float.__mul__` returns
NotImplemented -> Python falls back to RHS `__rmul__`). So a carrier with
`__mul__`-but-no-`__rmul__` genuinely breaks inside `ScaledOperator`/
`ZeroOperator` (both do `scalar*op`). WATCH: ndarray*ndarray elementwise sites
(`DiagonalOperator`, `RankOne`) are NOT scalar*vector and don't bear on the
contract.

(4) **Coverage honesty: "every leaf satisfies the Protocol" is ONE shared
base, not N independent leaves** -- AngularFlux/ScalarFlux/BoundaryFlux/
HarmonicMomentField all inherit the dunders from the `Field` base. The
genuinely-independent positives are the FAMILIES (np.ndarray native /
Field-base / delegating subclass), not the leaf instances; don't over-credit
3 leaf cases as 3 proofs. A docstring that names a specific carrier as covered
when no test exercises it is a documented-but-untested gap (L-007-flavored) --
recommend +1 line, flag don't block.

(5) **runtime_checkable + @property data members**: `isinstance` checks
PRESENCE of all members (`__protocol_attrs__`; `hasattr`-style, ignores
property-vs-attr and the property's return type). Rule out "all-True/all-False
by accident" with a Partial duck (missing one member -> False, so each member
is individually load-bearing) + a complete Duck (all members, NO inheritance
-> True, structural). An `_is_a(candidate: object, protocol)` helper whose body
is literally `return isinstance(...)` does NOT mask: the `: object` annotation
only launders the STATIC type so pyright skips its unsafe-overlap warning on a
concrete literal -- which is the asserted FACT, not a hazard.

(6) **pyright "= baseline, NO offset" rooting** (the trap: a real +N masked by
a coincidental -N). Airtight no-checkout proof = THREE facts together: (a)
full-tree total EXACTLY the stated baseline; (b) the SUT file isolation 0/0;
(c) the SEAM file (the one a reverted-risky-part touched) has an EMPTY `git
diff --stat` + isolation 0/0. The masked-offset hazard REQUIRES a nonzero diff
on the seam or a new error somewhere -> empty diff + unchanged total rules it
out. Always demand the seam file's diff be empty when a closeout says "I
reverted the risky part."

---

## L-040 -- algebra-law-suite and broadcast-oracle have DISJOINT coverage; the law-suite is swap-INVARIANT (demand a separate nx≠ny oracle for the variable-swap mode)

A multiplier-algebra law-suite (M_1=I, M_0=ZeroOp, linearity, self-adjoint, spectrum→CAP_SOLVE,
homomorphism) on a DiagonalOperator broadcast engine **cannot catch a variable-swap mode #2
(axis-ordering) bug** — linearity `M[af+bg]=aM[f]+bM[g]` and homomorphism `M[f]M[g]=M[fg]` are
ALGEBRAICALLY swap-invariant (probed: a CONSISTENT group/spatial transpose applied to all three
operands preserves both laws → allclose stays True). The axis-ordering bug is caught ONLY by the
**broadcast oracle** (`engine.apply ≡ sigma[None]*psi`) in the discriminating regime: a 2-D carrier
`(N_ord, ng, nx, ny)` with **nx≠ny** makes a spatial-axis transpose `(ng,ny,nx)` RAISE on broadcast
(`(1,2,3,5)` vs `(12,2,5,3)`), whereas a SQUARE mesh silently agrees in shape (no discrimination).
So the two test families are NOT redundant — law-suite = intrinsic-property gate (verify the laws
hold, mutation-proven: nonlinearity reds linearity, additive-offset reds homomorphism, non-unit
scale reds M[1]=I), oracle = the variable-swap catcher. **Rule**: do NOT credit an algebra-law-suite
with axis-ordering coverage; demand a SEPARATE nx≠ny broadcast oracle for the variable-swap mode.
The ≥2G-asym-het requirement on linearity/homomorphism (anti-pattern #3/#4) is about NOT NULLING
the group/spatial structure (so the laws are exercised on real coupling), NOT about catching swaps.

**CAP_SOLVE behavioral-strengthening review (anti-pattern #11, BOTH-tested)**: the promotion adds
an honest spectrum gate (CAP_SOLVE iff min|f|>0) where the legacy CollisionOperator advertised it
always → silent IEEE NaN on σ=0. POSITIVE+NEGATIVE both present and teeth-proven: emulating the
legacy always-on bug (monkeypatch engine to force CAP_SOLVE) reds `test_spectrum_cap_solve...` with
the exact `-O`-firing message. Audit-confirmed safe: 3 prod sites all use σ_t (bounded away from 0
via S2 `total_cross_section_field`); `CollisionOperator.solve` has ZERO prod callers (the WDD sweep
is `InvertibleOperator.solve`, which has its OWN stricter construction-time `σ>0` check at
operator.py:784); the σ_r removal-fold path that COULD go ≤0 is issue #200/#215 — documented, NOT a
live code path. So nothing relied on the old always-on CAP_SOLVE.

**Mode-11 (gate-reaches-new-code) on a promotion**: the resolvent gates (kinf_homogeneous,
si_carve) CONSTRUCT the promoted C (`C_init` 458/22) + fire its `__post_init__` spectrum gate +
read `C.sigma` (19900/1705 — σ threaded into the WDD sweep), but `MultiplicationOperator.apply/solve`
are NEVER called (0/0 — `InvertibleOperator` OVERRIDES apply via loss-rep and solve is the sweep).
So the resolvent gates cover C's CONSTRUCTION+σ-threading; the apply/solve ARITHMETIC is covered by
the NEW broadcast oracle. Honest, complete, non-overlapping — but state it explicitly (don't claim
the kinf gate exercises `M.apply`). **Field-promotion is a label, 0 ULP**: `CrossSectionField.from_mesh(arr).values IS arr` (same object), so the S2 σ_t→CrossSectionField rewire is a pure retype;
broadcast oracle confirmed exactly 0 ULP (both forms are `expand_dims` on axis 0, reduction_depth=1).
Mode-8 clean (0 bare asserts in new test+prod; all gates via `_require`/`pytest.fail`/`np.testing`).

---

## L-041 -- cofree base-extraction "bit-identical carrier" claim: prove the polymorphic `_recombine` hook bites TWO mutation ways + the dedicated unit test must be byte-UNCHANGED vs HEAD

S4.5 extracts a TIMELESS `FullField` base (the 6 vector dunders + `to_flat`/`from_flat` + `copy`/`zeros` + validation, lifted ONCE) out of `TimedFullField`; the timed subclass keeps `_history`/`history_depth`/`advance`/`at_lag` and OVERRIDES `_recombine` (returns TimedFullField, empty history, preserved depth). The load-bearing claim is BIT-IDENTICAL `TimedFullField` behavior (pure carrier extraction, no math).

**Verification recipe for a base-extraction "behavior-unchanged" claim:**
1. **The dedicated unit test must be byte-UNCHANGED vs HEAD** — `git diff <HEAD> -- tests/.../test_timed_full_field.py` MUST be empty. A passing-but-EDITED test is weaker evidence (the edit could have relaxed a path). #257 S4.5: `test_timed_full_field.py` diff empty, 38 pass under `-O`. The NEW `test_full_field.py` is a SUPERSET (adds the recombine teeth + discriminating membership), not a replacement of coverage.
2. **The polymorphic-hook teeth bite TWO mutation ways** (do not stop at one):
   - (a) override returns the BARE base type → "type preserved" tooth reds (3 algebra tests: `type(out) is TimedFullField` fails "got FullField").
   - (b) override DROPPED → base `replace(self,...)` runs → KEEPS history (replace copies the class AND the field) → only the EMPTY-history tooth reds. This is the realistic "forgot to override" mutation and is caught ONLY because the empty-history test ADVANCES first (builds real history_length==1 as a precondition) — without that precondition the `out._history == ()` assertion is TAUTOLOGICAL (zeros()-input already has empty history). The advance-first precondition is the load-bearing non-tautology.
3. **`from_flat` made generic (`template: T -> T`, routes through `template._recombine`)** is pinned by AttributeError-teeth: mutate it to return a bare `FullField` → `from_flat_drops_history`/`_preserves_history_depth`/`iteration_protocol_detection` red with `'FullField' has no attribute history_depth/history_length` (the timed-only attrs are the discriminator).
4. **Discriminating type-check became NOMINAL** (was runtime-checkable Protocol `TransportState`, now concrete `@dataclass FullField` isinstance): confirm at runtime `isinstance(ndarray, FullField) is False` AND `issubclass(TimedFullField, FullField) is True` AND `isinstance(FullField, type) is True` (a real class, not Protocol-only) — anti-#11 positive+negative+timeless/timed all present, Mode-8-safe via `_require`/`pytest.fail` (0 bare asserts).
5. **type:ignore accounting** — count in BOTH files vs the HEAD original: net-new must be 0. S4.5: HEAD `timed_full_field.py` had 2 (`zeros_on` `[attr-defined]`); they MOVED to the new `full_field.py` (still 2), `timed_full_field.py` now 0. Net 0. The main-agent post-fixes (generic `from_flat` removing a `[override]` ignore; `zeros` delegating to base de-duping `zeros_on` ignores back to 2) are ignore-REDUCING, not -adding.

**Baseline-red triage**: the 7 regression reds (#250 SPHERE ×5 huge-ULP ~1e15 stale-snapshot per L-034/L-036 + #232 mu_y ×2 ValueError) are geometry-scoped — SLB/CYL/Cartesian arms in the SAME files (all using the `TimedFullField.zeros(...)` public API) PASS. A carrier-refactor break would fail ALL arms (type/attr error), not just SPH/mu_y (geometry-math/quadrature). Closeout's baseline-worktree (`93aa016` + symlinked `.venv`) independently saw the same 7. pyright EXACT `2295 errors, 19 warnings`.

**Pyright baseline-comparison gotcha (from closeout, worth keeping)**: a git worktree pyright count is ONLY comparable with the MAIN `.venv` symlinked into the worktree root (pyrightconfig `venv: .venv`); a worktree without it analyzes a different file set → bogus count (2922 vs 2295).

---

## L-042 -- additive-only SUT: prove the "baseline-invariant" claim by grepping for ANY importer (empty ⇒ cannot perturb), and the variable-swap teeth are the contraction AXIS not the operand order

S5 = additive only: 2 new untracked SUT files (`orpheus/numerics/functional.py`
runtime_checkable `Functional` Protocol; `orpheus/transport/production_rate_functional.py`)
+ 1 additive `numerics/__init__.py` export. NO edit to operator.py/fission.py/solver.py
(`git diff --stat HEAD` empty for all three). **Additive proof for "baseline-7 invariant"**:
`grep` the tracked tree for any importer of the SUT EXCLUDING the 4 new test files →
EMPTY ⇒ no pre-existing consumer ⇒ S5 cannot perturb any pre-existing outcome (stronger
than re-running the reds). transport+numerics dirs 915 passed/1 skipped; fission 18 passed.

**Bit-id premise must be checked, not assumed**: `CrossSectionField.from_mesh(nsf,sn).values`
is `array_equal` to raw nsf (the producer doesn't transform), so SUT
`(nu_sigma_f.values*phi).sum(axis=0,keepdims=True)` is 0-ULP `array_equal` to the legacy
`RankOneOperator.apply` `inner` line 1776 (`(right*x).sum(axis=axis,keepdims=True)`,
right=νΣf, axis=0). Correctness rides a GENUINELY structurally-independent ref
(`hand_derived_production_density` = explicit nested-Python double-loop, no numpy reduction,
no ORPHEUS algebra) — L11 clean. B.2 RankOne-equivalence is correctly DEMARCATED as de-risk
not correctness.

**Mode-2 framing trap (test author got it RIGHT)**: a literal νΣf↔φ swap is VALUE-INVARIANT
(pointwise product commutes — verified `array_equal`). The genuine Mode-2 hazard is the
CONTRACTION AXIS, discriminated by the nx≠ny (5≠3) mesh by SHAPE. Don't accept "swaps the
operands" as the swap teeth; the axis is the teeth.

**Mutation results (L12)**: `axis=0→1` reds exactly 6 production-rate gates (matches
test-architect claim); category gates stay green (assert no number). Shape-preserving Mode-3
`×1.5` magnitude factor reds 5 gates incl. dedicated `test_density_unweighted_by_cell_volume`
— real measure-fold teeth. **CAVEAT (novel)**: `keepdims=True→False` does NOT redden — the
`squeeze_density` helper collapses a leading length-1 axis, making the WHOLE suite agnostic to
keepdims (by design per helper docstring; even B.2 bit-id passes squeezed). So the user's
"confirm keepdims reddens" is the ONE unsupported sub-claim: the suite tolerates either output
rank, it does not ENFORCE keepdims. Not a correctness bug (value+bit-id hold either way) — but
the bit-identity-with-legacy claim is enforced only up to a squeeze, not literal rank.

**Category teeth (Claim 2) — the headline gate is the WEAK one**: runtime_checkable
`LinearOperator` checks 4 members (`apply,capabilities,domain,codomain`); `isinstance(func,
LinearOperator)` only goes True on a FULL 4-member leak. The realistic leak (a Functional that
grows `apply`+`capabilities`) leaves the headline `isinstance NOT LinearOperator` BLIND
(domain/codomain missing → still False) but reds the 3 surface gates `lacks_apply`/
`lacks_capabilities`/`disjoint`. Pytest-level leaky-probe mutation: all 3 surface gates red,
headline blind on partial leak. The 3 surface gates ARE the defense (test-architect flagged
this correctly). Mode-8 clean: 0 bare `assert` in any S5 test file; all route through
`require`(pytest.fail)/`np.testing.*` (fires under -O). pyright EXACT `2307 errors, 19 warnings`
(=user's b404ae1 baseline; plan's 2295 was stale-worktree); both SUT files individually 0/0/0
→ no masked offset. 1 skip = by-design (no estimator wrapper shipped).

---

## L-043 -- Mode-11 sentinel for a NEW PRIVATE adapter = install it as a PYTEST PLUGIN that wraps the internal call (prove scipy FIRED it, not just BUILT it)

Bit-identical refactor: 2 inline scipy-LinearOperator closures (`A_matvec`/`M_matvec`)
in `KrylovAcceleration.solve` → 1 named carrier `loss_minus_gains(psi)` + 2
`_as_scipy_linop(carrier, template, n)` calls; retired public `as_scipy_linop` (0 callers)
+ orphaned `spla` import + 5 tests + 3 doc xrefs.

**CLAIM-1 byte-id (Mode-2 A/M template-swap):** prove the non-swap TWO ways — (a) read the
2 call sites (A→`solution_template`, M→`q_ext`), (b) RUNTIME binding sentinel: monkeypatch
`_as_scipy_linop` + `KrylovAcceleration.solve` (stash `q_ext`/`solution_template` ids in
solve, compare template id in the adapter) → reported `carrier=loss_minus_gains
bind=solution_template` / `carrier=<lambda> bind=q_ext`. A swap would have inverted these.
`loss_minus_gains` reduction order char-identical to old `A_matvec` (L.apply first, then
`for g in self.gains: out=out-g.apply`); `(n,n)`/`dtype=float` preserved.

**CLAIM-2 Mode-11 sentinel for a NEW PRIVATE adapter — sharper than L-031/L-033 in-process
probes:** install the sentinel as a PYTEST PLUGIN (`-p <module>`, module must be on
PYTHONPATH — `-p /tmp/x` fails "No module named", copy to cwd + `PYTHONPATH=$(pwd)`), patch in
`pytest_configure`, restore + summarize in `pytest_unconfigure`. Wrap `linop._matvec` (the
internal scipy calls, NOT `.matvec`) with a counter to prove scipy FIRED it, not just BUILT
it. Tag A vs M by `carrier.__name__` (`loss_minus_gains` vs the precond lambda). Result on
identity-precond[slab] (non-None `lambda q:q` → M built): A built=2/fired=160, M
built=2/fired=161, both on TimedFullField → M-template wiring exercised on the REAL typed path.
This is the gold-standard Mode-11 evidence WITHOUT mutating any tracked production file (L28).

**CLAIM-3 retirement:** word-boundary grep (`[^_]as_scipy_linop`) = 0 hits in orpheus/tests/docs
(the `_as_scipy_linop` private hits are noise). 5 deleted tests pinned a now-gone
`LinearOperator`-taking public adapter; its only unique assertion (`MissingCapability` on
missing `CAP_APPLY`) maps to a behavior the NEW bare-callable adapter does NOT have — the
equivalent guard moved UP to `KrylovAcceleration.__init__:422` (composition-time, STRONGER),
covered by 14 `MissingCapability` refs in test_iteration.py + 3 surviving `NoApplyOperator`
negatives in test_operator.py. No non-redundant coverage lost.

**pyright net-new=0 PROOF without mutating tree:** the 3 `reportCallIssue` at the new line 228
(`spla.LinearOperator((n,n),matvec=...,dtype=float)`) are the scipy-stub false-positive that
existed at HEAD across `grep spla.LinearOperator HEAD:iteration.py` = 2 sites (756,766) + a 3rd
in the deleted public adapter → refactor CONSOLIDATED 3 sites → 1, which is why total dropped
2307→2297 (−10). `# type: ignore` delta −1 (removed `op.apply_transpose`), 0 added. Gates:
138 Krylov/round-trip pass (-O); broad regression 7 reds = EXACTLY #250 SPH×5 (huge-ULP ~1e15
while SLB sibling 1-ULP DriftWarning-pass = L-034 stale-snap) + #232 mu_y×2, all in
tests/sn/operators/ with 0 refs to the changed code (orthogonal). No ERR (no bug caught).

---

## L-044 -- "producer-now-emits-X, test-helper-decoupled" re-baseline integrity: prove the helper rebuilds the flat baseline TEST-SIDE or it silently inherits the new emission

S9 makes `SN2DCartesianLDStressMMSCase.prescribed_inflow` EMIT the moment-resolved
face slot (slot-0 transverse cell AVERAGE, slot-1 bare transverse P1 slope) via a
case-owned leggauss-only `_project_inflow_to_face_moments`, gated on
`face_moment_count>1` (DD/Step byte-identical). NO new field type, NO value gate
(slope is sub-floor for converged flux — vv Mode-10 companion-unavailable, 3rd recurrence).

**Re-baseline-integrity recipe for a "production-now-emits-X, test-helper-decoupled"
change** (the HIGH-PRIORITY trap): when a producer (MMS `prescribed_inflow`) gains a
new emission AND a test helper's "flat baseline" branch USED to route through that
producer, the helper MUST rebuild the flat baseline TEST-SIDE or it silently inherits
the new emission → toggle collapses. PROVE the decoupling kept teeth by probing the
helper's 4 legs directly: `None==zero` byte-id (slope-free baseline + no-op control
has teeth), `|mom−None|`/`|flip−None|`/`|mom−flip|` all ≫tol (slope consumed, sign
matters), `None≠mom` (toggle not vacuous). #257 S9: None==zero byte-id, |mom−flip|≈2.19e-2.

**Sign-mutation gate teeth proof (cheap, in-process, no prod edit):** monkeypatch the
SLOPE SOURCE (`_face_transverse_buffers`) to zero the slope, re-run the mom-vs-flip
comparison → `|mom−flip|/|φ|` goes 4.10e-3 (healthy, ~5.6 orders >1e-8 `_CONSUMPTION_TOL`)
→ 0.000e+00 (bug) → gate reds. Confirms the consumed-flip is genuine, not a tautology.

**Mode-11 producer-stamp NOT circular** (L-029 applied): GATE-B compares production
`case.prescribed_inflow` slot-1 vs `_face_transverse_buffers` (test-side leggauss),
NOT vs `case._project_inflow_to_face_moments` — two INDEPENDENT leggauss impls; GATE-C
separately pins their agreement (array_equal, maxdiff 0.0). A sign error in the prod
projector would NOT propagate into the test ref → GATE-B reds. Sentinel-instrument
`LR._LossRepresentation._inflow_to_moments` (the genuine prod consumer, reached 688×/solve):
flat/zero slot1==0, mom/flip slot1=1.9e-2, `|mom−flat|`phi_sum>1e-3 (consumed), zero==flat
byte-id. Mode-11 closed: producer IS exercised, not a surrogate.

**Verdict-pin teeth** = `improves`(mom<flat) check + `|mom−flat|/flat ≤ 0.30` band.
At bc_scale=20× (strongest amplification), mom monotonically WORSE (improves all False,
rel max 0.205 < 0.30, orders [1.7,2.4]) → sub-floor wall fundamental. Pin reds if slope
ever becomes above-floor. Coherent-promise gate teeth = flat first-cell-row order ≥1.85
(measured 1.99/2.00/2.00 — average alone delivers O(h²) at boundary, no asterisk).

**DD byte-identity proven 3 ways:** (a) `np.array_equal(prod_DD.values, pre-S9 face_coords
build)==True` (1344,); (b) GATE D strict `-W error::DriftWarning` 520/1/4 = baseline, NO
DriftWarning fired; (c) no LD-stress consumer in tests/sn/sweep/core or solve (grep) → no
value/snapshot pin could shift. Gates: G1 35pass / G2 590pass,1skip,4xfail / GATE-D 520/1/4
/ pyright 2282 = baseline 0 net-new. Mode-8 clean (0 bare assert in new file or prod).
NO blocker, NO false-green, NO ERR.

---

## L-045 -- a "behavior-neutral field-zeroing" claim is valid ONLY for the ONE fission/emission contract it was proven against; re-prove inertness for EVERY consumer (ERR-063)

SUT = `EmissionSpectrum(np.ndarray)` value-object + `Mixture/Isotope.__post_init__`
simplex/null χ guard (keyed `is_fissile = bool(np.any(SigF>0))`) + a "behavior-neutral"
precursor zeroing non-fissile χ on shared `xs_library` regions B/C/D. The TYPE + guard +
intrinsic-property gates are SOUND; the precursor is NOT behavior-neutral → BLOCK.

**What's SOUND (mutation-/byte-verified):** (1) intrinsic gates vv#11 BOTH legs, hand-laid
L11 refs; negativity clause INDEPENDENT of sum (mutate prod: drop `>=0` → ONLY
`test_negative_entry_raises_even_when_sum_is_one` reds; drop sum → 2 sum legs red, negativity
green; relax `assert_null` to atol=1e-6 → `test_any_nonzero_raises` reds, pinning STRICT
exact-zero). (2) Mode-8 clean (0 bare asserts; all `_require`/`pytest.fail`/`np.testing`/
`pytest.raises`; 28+13 pass under -O incl real-GENDF). (3) SN-path behavior-neutrality REAL:
direct re-solve of het 3-region DD (fuel A + non-fissile mod B) with mod.chi=[1,0] vs [0,0]
→ keff `1.2298233055738448` BYTE-identical + flux array_equal (max abs diff 0.0); confirmed
SigP≡0 on B/C/D so SN `FissionOperator` χ·(νΣf·φ)≡0. (4) is_fissile/SigP seam (item 5):
explorer audit + GENDF MF6/MT18-co-located-with-MF3/MT18 → NO real-data production path has
nonzero-χ-∧-zero-SigF; only the synthetic billiard fixture (reads SigP/chi, never SigF —
SigF=nu_sf injection inert, confirmed billiard.py:1031-1032 + tree-grep). DD reg 13pass / TA
full 107pass,2xfail (matches closeout).

**THE BLOCK (ERR-063):** "zeroing non-fissile χ is inert" assumed the SN/`compute_macro_xs`
contract (χ gated by SAME region's νΣf). FALSE for `solve_peierls_mg`: its MG fission op
`B[i,ge,j,gs] += K[i,j]·chi[i,ge]·nu_sf[j,gs]` weights SOURCE-region νΣf by SINK-region χ, so
χ on non-fissile region B (the emission spectrum of fission BORN in A but emitted INTO B) is
LOAD-BEARING. Direct probe: region-B χ [1,0]→[0,0] moves peierls k_eff `1.0985→0.5563` (1G/2R)
/ `1.1008→0.3856` (2G/2R) — O(1), not ULP. 7 L1 tests in
`tests/derivations/test_peierls_rank_n_class_b_mr_mg.py` (cylinder/sphere hebert overshoot +
recovers_kinf[2G_2R]+RICH + mark_floor[cyl/sph]) FAIL under S10a, PASS at clean HEAD (proven
via `git worktree add c6e21c0` + PYTHONPATH=worktree: 4+4 passed). Only RICH is @slow; other 6
plain @l1. Closeout MISSED it — it relied on "0 EmissionSpectrum reds" (counts only guard
ValueErrors, blind to silent accuracy regressions) + "DD snapshots didn't move" (DD = SN-only,
the consumer where χ IS inert) + never ran the 494s peierls suite. Test authors had ALREADY
flagged this χ-dependence (commit 76b11e8, Issue #132).

**RULE (new):** a "behavior-neutral field-zeroing" claim is only valid for the ONE
fission/emission contract it was proven against. When the field is a SHARED source feeding
consumers with DIFFERENT contracts (same-region χ·νΣf vs sink-region-χ × source-νΣf), re-prove
inertness for EVERY consumer with a DIRECT old-vs-new VALUE comparison (O(1) move = not neutral),
NOT a fast proxy ("snapshots didn't move" / "no guard errors"). Run the slow accuracy-band suites
that consume the edited field. L20 shared-source hazard + H5 (test count ≠ coverage). Recommended
fix: do NOT zero the shared library χ — decouple peierls cases' χ from the guarded library, OR
key the guard on production not SigF, OR restrict the guard off placeholder library regions.
Worktree-baseline recipe (clean-HEAD confirm): `git worktree add /tmp/x <HEAD>` + run with
`PYTHONPATH=/tmp/x` (editable .venv else imports MAIN tree — verify `orpheus.__file__`).

---

## L-046 -- for a WEIGHTED value-pin, L11-independence is not enough: the hand-ref must carry EVERY weight factor AND the fixture must make a factor-BLIND formula give a different answer

SUT = NEW `production_weighted_chi(isotopes,sigF,aDen,fissile_indices)` helper
(`χ_mix = weights @ fissile_spectra`, `w_i = aDen_i·Σ_g ν̄σf_i / Σ_j(…)`) replacing the
first-fissile χ shortcut in `compute_macro_xs`. The S10b consumer that produces the
multi-fissile χ_mix the [[L-045]] S10a guard validates for free (gate-1 interlock). ALL 7
review points SUPPORTED; clean, no blocker.

**⭐ THE #1 SCRUTINY — hand-ref structural-independence for a WEIGHTED value-pin.** When the
only value-pin is a hand-laid convex average, L11-independence is NOT enough: the hand-ref
must independently carry EVERY weight factor (here `aDen` AND the `Σ_g ν̄σf` production sum),
AND the fixture must make a factor-BLIND formula give a DIFFERENT answer or the factor is
untested (a vacuous pin). PROVE it two ways: (1) the fixture discriminates — gate 2 uses
unequal `aden=[2,1]`, so aDen-aware w=[0.4545,0.5455] vs aDen-blind w=[0.2941,0.7059] → χ
differs by 0.128 ≫ atol=1e-12 (compute by hand, don't eyeball); (2) MUTATE the production to
be aDen-BLIND (drop the `aDen[i]` factor) → gate 2 reds, gate-1 simplex + gate-3 byte-id STAY
green (blind formula is still a convex average of simplices = a simplex; single-fissile
unaffected). Hand-ref here is genuinely independent: explicit scalar `p_i = aden[i]*nubar_i*sigf_i`
term-by-term (single-nonzero-group fact), NOT `weights @ fissile_spectra` re-spelled. The
aDen-blind variant IS the "shares the code's weight-derivation / forgets aDen the same way"
failure the brief warned of — proven defeated.

**Other teeth (all mutation-verified live, not trusted from closeout — L12):** legacy
first-fissile shortcut → gate 2 + real-UO2 smoke red, gate1/gate3 green; unweighted mean →
gate 2 red (sole catcher); non-convex `2·(weights@spectra)` → gate-1 S10a `assert_normalized`
interlock fires (7 red) = the interlock is LIVE. Single-fissile collapse is EXACT byte-id
(w=[1], max abs diff 0.0). Mode-7 honest-scope ("flat-flux representative, NOT flux-exactness")
declared in helper docstring + gate file. Mode-8 clean (8/8 under -O, all `_require`/`np.testing`).

**Byte-identity scoping (re-baseline list EMPTY, independently confirmed):** DD regression
13pass (within-tol DriftWarnings pre-existing FP-noise); DD path never touches `compute_macro_xs`
(grep empty — builds Mixture via `xs_library.make_mixture`). Multi-fissile `compute_macro_xs`
callers = `uo2_fuel`/`pwr_like_mix` (`fissile_indices=[0,1]`); the ONLY pytest-collected
consumer is `test_solver_components.py::test_profile_421g` = a TIMING test (prints ms, asserts
NO k_eff/flux); `pwr_like_mix`/other `uo2_fuel` refs all in `examples/` (NOT in
`testpaths=["tests"]`). So no committed test pins a converged value off a multi-fissile mixture
→ the χ-value change rests entirely on gate 2's hand-ref. pyright net-new = 0 (mixture.py 3
errors WITH==WITHOUT change via stash; all 3 pre-existing `SigP/Sig2/SigT = sum(...)` int-noise,
#226; full project 2353==baseline). The closeout's `weights @ fissile_spectra` deviation (vs
brief's `sum(generator)` which costs +1 `reportReturnType`) verified principled — a convex
average IS a matvec.

---

## L-047 -- a runtime_checkable category Protocol's member-presence loophole is REAL; the direct `not hasattr(...)` negative gates are the defense (+ the partial-coverage `kernel ≠ full apply` caveat)

S6 is ADDITIVE + bit-identical: a §5.6 `IntegralKernelOperator` Protocol
(`orpheus/transport/integral_kernel_operator.py`, `@runtime_checkable`,
sole member `kernel`), a `FissionOperator.production_rate` property (S5
`ProductionRateFunctional` over νΣf), a `ScatteringOperator.kernel`
property (`OperatorProduct(R, OperatorProduct(Λ, M))`, `skip_l0=True`).
All 5 claims SUPPORTED; 2 caveats, no blocker, no false-green.

**Claim 1 (category teeth) -- `runtime_checkable` member-presence loophole
is REAL but the direct-attr gates close it.** Monkeypatch `ident.kernel =
"fake"` on `IdentityOperator()` → `isinstance(ident, IKO)` flips True
(the documented S5 loophole; isinstance only checks PRESENCE). The 3
negative gates that assert `not hasattr(..., "kernel")` directly
(`*_lacks_kernel`, `*_lacks_kernel_and_apply`) are the defense-in-depth
that does NOT depend on the Protocol machinery. Discriminator
(`IdentityOperator` IS a LinearOperator but NOT an IKO) proves a strict
refinement, not a `LinearOperator` alias. 20/20 green under -O.

**Claim 2 (fission) -- B.1 L11-clean, B.2 Mode-11-live.** B.1 reference
`hand_derived_fission_emission` = explicit Python double-loop, shares NO
numpy reduction with production (role-swap sensitive, verified
max-rel-diff 10×). B.2 reads `op.production_rate` OFF the live operator;
MUTATION-VERIFY (point production_rate at `total_cross_section_field`
instead of νΣf) reds BOTH B.2 gates @100% mismatch. `evaluate` =
`(nu_sigma_f.values*phi).sum(axis=0, keepdims=True)` = the RankOneOperator
`inner` line byte-for-byte → bit-id is structural. ⭐ ASYMMETRY:
`fission.kernel` IS the FULL F (production reads it at fission.py:454/471);
`scattering.kernel` is the aniso ℓ≥1 part ONLY.

**Claim 3 (scattering) -- Mode-11-live + the skip_l0 blind-spot CAVEAT.**
`S.kernel.apply == _aniso_source_from_moment_values(M·ψ)` @ 0 ULP, reads
live `S.kernel`. MUTATION: drop R → 2 gates red (value + shape moment-
tensor); `skip_l0 True→False` → value gate red (subtle flag IS load-
bearing). ⚠ CAVEAT: `S.kernel.apply` (L2≈0.98) is ~5% of the full
scattering source (L2≈21.9) -- ONLY the ℓ≥1 aniso redistribution, pre-1/W;
P0 in-scatter + n2n are NOT in it. DOCUMENTED honestly (module + property
docstrings, "genuinely-nonlocal-in-angle part", "P0/n2n are LOCAL/separate
components"; test docstring says "pre-1/W") but NO production reader and NO
test POSITIVELY asserts `kernel != full apply`. A future consumer mistaking
`kernel` for full S silently loses ~95% of the source. Recommend a 1-line
gate `require(not allclose(kernel.apply(ψ.values), S.apply(ψ).values))` to
pin the partial-ness; minor follow-up, not a blocker (the only current
consumer is the cross-check, which knows the semantics).

**Claim 4 (matvec arms byte-id) -- 17/17 green.** TestAnisoMomentSourcePath,
TestProtocolCompliance, TestP0AlgebraicIdentities,
TestRankOneTensorProductKernel, TestBitIdenticalToLegacyInlinedMath all
green; aniso MMS `test_curvilinear_aniso_scattering_p1.py` 2/2 green (the L1
physics reference for scattering).

**Claim 5 (pyright) -- CONFIRMED 0 net-new on production.** CLI `npx pyright`
= 2311 errors / 19 warnings (= user's number; S5 base 2307 + 4 from the new
test skeletons). fission 8, scattering 22 — proved 0 net-new by stash-tracked
+ hide-untracked-module → TRUE baseline (S6 reverted) = 30 = 8+22 unchanged.
The 3 `cast(LinearOperator, ...)` in scattering.kernel are a legit PEP-484
bridge (MomentProjection/LegendreMomentScattering/HarmonicMomentReconstruction
all carry `apply` at runtime; composition green @ 0 ULP) for the
unparametrised-LinearOperatorMixin generic gap (#226). 0 `# type: ignore`
added (the one match is inside an explanatory comment).

**Claim 6 (baseline reds = 7) -- CONFIRMED.** 5 #250 SPHERE stale-snapshot
reds (test_streaming_operator.py TestT4cPreT4RegressionSnapshotCurvilinear
×2 + test_bc_extraction_matvec.py SPH ×3, max ULP ~8.77e15 = L-034 stale-snap
signature) + 2 #232 mu_y. S6 touches none of streaming/matvec snapshot code
(additive only) → reds pre-exist.

**Mode-8 -- the test-architect's flag on `TestRankOneTensorProductKernel`
bare asserts is WRONG (re-confirms L-010).** The 4 NEW S6 test/helper files
have ZERO bare asserts (all route through `require`=pytest.fail / np.testing.*).
The existing `TestRankOneTensorProductKernel` (lines 365-411) DOES use bare
`assert isinstance/is/==` — BUT these FIRE under -O because pytest's assertion
REWRITER rewrites asserts in collected tests/ modules at import time, BEFORE
-O would strip them. PROVEN twice: (a) broke the kernel → bare-assert
`isinstance` red under -O with AssertionError; (b) `assert 1==2` probe in a
tests/_tmp_probe module FAILS under -O. So `TestRankOneTensorProductKernel`
is NOT a Mode-8 gap; the test-architect's "S6 should fix it" is unnecessary
(Mode-8 is a concern ONLY for bare asserts in `orpheus/` production, not in
collected `tests/`).

---

## L-048 -- "behavioral-neutral codomain re-point" = TWO legs: identical failing-test IDs vs a read-only baseline worktree AND the dropped field has ZERO production `.advance(` callers

#257 S8a: SN operator matvec leaves (`StreamingOperator/InvertibleOperator/
MultiplicationOperator/SNBoundaryOperator.apply`, the `S`/`F` TimedFullField-input
arm) re-typed to EMIT the timeless `FullField` instead of history-bearing
`TimedFullField` (cofree-comonad finding: an operator is a base arrow
`FullField→FullField`; only the iteration DRIVER carries the comonad). Claimed
value-neutral / bit-identical. VERDICT SUPPORTED — recipe:

**CLAIM-1 value-neutrality = TWO independent legs.** (a) Reconcile the
baseline-red set against a READ-ONLY `git worktree add -d HEAD~ /tmp/x` checkout
(L28 — never mutate the working tree): run the EXACT same `-O` gate on both;
S8a-tree and baseline must produce IDENTICAL failing test IDs (here 7: #250
SPHERE×5 + #232 mu_y×2; pass-count delta = +14/+1xfail = exactly the new C5
file, nothing else). ZERO non-baseline reds. (b) Prove the dropped `_history`
is genuinely unused in steady-state: grep + Nexus `context` for ALL production
`.advance(` callers — if ZERO `calls` edges (only docstrings/prose), the history
shift-register is test-only and dropping `_history=()` cannot perturb any
converged value. Here confirmed 0 production callers.

**The reattach mechanism (don't take on faith).** The driver re-attaches the
timed type via `TimedFullField.__add__`'s `_recombine` hook — CONFIRM the timed
operand is on the LEFT of the `+` (`rhs = q_ext + g.apply(psi)`, q_ext timed →
`rhs.__add__(FullField)` → `self._recombine` resolves to `TimedFullField._recombine`).
The reverse order (`timeless + timed`) would resolve to the BASE hook and yield
timeless — a silent history-drop. Also confirm the resolvent `L.solve` STILL
returns TimedFullField (re-mints the iterate each step). Krylov path is
`FullField−FullField→FullField` throughout (unravels to scipy flat, reconstructs
from solution_template — never relies on `__add__`).

**CLAIM-2 scope (no math drift).** For a "type-surface only" production diff,
PROVE it by filtering the diff: every `^[+-]` line must be a type annotation
(`"TimedFullField"`→`"FullField"`), a docstring/comment, or whitespace — ZERO
numerical expressions (`sigma`, `values`, `einsum`, `out_bulk`, the `(L+C)−C`
arithmetic `lpc.bulk.values - sigma_t[None]*psi.bulk.values`). Confirm the
INPUT dispatch (`@apply.register def _(self, psi: TimedFullField)`) is UNTOUCHED
(only the return annotation + the output construction `TimedFullField(...)→
FullField(...)` dropping `_history=()`/`history_depth=` changed). A drift into
the NEXT sub-stage's behavioral change (here S8b pure-L) would falsify neutrality.

**CLAIM-3 teeth + Mode-11 (the matvec leaf has ZERO graph callers — reached only
via OperatorSum/driver).** The codomain gate (C5) MUST call `L.apply`/`C.apply`/
`F.apply` DIRECTLY (not solve-only — solve routes through the sweep/loss-rep and
never touches the matvec emit path). Mutation-verify teeth by REVERTING the
re-point on ONE leaf (make `StreamingOperator.apply` emit `TimedFullField` again,
QA-MUTATION-SENTINEL comment) → the C5a `type(out) is FullField` checks go RED
across all geometries with a precise diagnostic (`got TimedFullField`); revert +
confirm green + `grep QA-MUTATION-SENTINEL` = 0 residue. The legacy snapshot
gates (`TestT4b/c`) DO reach the matvec leaf (`L.apply(state)` directly) — verify
slab/cyl arms reproduce frozen bulk (`assert_regression kind=direct`) + STRICT
0-ULP boundary; SPHERE arms are #250 stale-snap (O(1) value diff, L-034), red on
BOTH trees = pre-existing.

**CLAIM-4 re-pointed B-tests.** The ~41 re-pointed tests are clean type-surface
updates: `isinstance(out, TimedFullField)` → `isinstance(out, FullField)` +
`not isinstance(out, TimedFullField)`, DROPPING the now-meaningless
`out.history_depth==depth`/`out._history==()` assertions (they test the EXACT
attribute S8a removes), while PRESERVING all value assertions (`.bulk.mesh is`,
`isinstance(.bulk, AngularSourceSink)`, `.values.shape`, boundary `==0.0`). The
roundtrip tests (`test_removal_form_matvec_sweep`, `test_invertible_operator`
solve∘apply=id) re-wrap the timeless `op.apply` output into a TimedFullField
before feeding `solve` — mirrors the driver; byte-identical source `.bulk`/
`.boundary`. NIT (cosmetic, not blocker): a few function NAMES still say
`..._timed_full_field` despite now asserting timeless. Gate counts: pyright
2297/19 (=baseline, 0 net-new); regression 7 baseline reds only; L1/MMS 40pass/
2xfail (converged limit unmoved); C5 14pass/1xfail (sphere-krylov #200 xfail).

---

## L-049 -- when a value-moving carve preserves the COMPOSITE not the leaf, prove byte-id on the composite directly; pin a re-baselined `.npy` to a STRUCTURALLY-INDEPENDENT reference, not "whatever the leaf emits"

VERDICT SUPPORTED. The value-moving core of the streaming carve: production's
within-group matvec uses the COMPOSITE `(L+C).apply`=`InvertibleOperator.apply`
(rides `loss_action(σ_t)` UNCHANGED) — so the value-preservation story rests on
composite byte-identity, NOT the standalone pure-L leaf.

**CLAIM-1 (composite byte-id) is the load-bearing one — prove it directly.**
Emit `(L+C).apply(ψ)` on live vs a read-only `git worktree add -d 9316321`
baseline, per geometry (slab/sphere/cyl) ≥2G het, `PYTHONPATH=$PWD` to OVERRIDE
the editable .venv (confirm baseline ran baseline code via `inspect.getsource` +
dataclass `fields` — `sigma_t` present/absent is the tell). Result: 0 ULP,
absdiff=0.000e+00 all 3 geoms. The pure-L LEAF drifts (CART 32 / SPH 12 / CYL
117 canonical numpy nulp, boundary STRICT 0-ULP exactly) — "≤16 ULP" in a brief
can UNDERSTATE the leaf drift (CYL 117); it's still genuine FP-reassoc
(large-mag → moderate-ULP, rel ~1e-15) and the leaf has ZERO graph callers
(Nexus `callers` total:0), so inconsequential. Test bound `_BULK_NULP=256`
covers it; T4b snapshot DriftWarnings up to 192 ULP all within 256.

**CLAIM-4 re-baselined .npy (the headline laundering risk) — pin to the
STRUCTURALLY-INDEPENDENT composite, not "whatever pure-L emits".** The 3
`bc_extraction_2d` `.npy` were re-captured. The decisive check is NOT
"committed==pure-L" (circular) but "committed == `(L+C).apply.bulk − σ_t·ψ`"
(= the BYTE-IDENTICAL composite minus the collision diagonal) — measured ≤64
ULP (rel ~1e-16). Because the composite didn't move a single bit (CLAIM 1) AND
pure-L = composite − collision to ULP, the frozen `.npy` IS the genuine pure-L
value. MASKING-CHECK: OLD baseline `.npy` ≠ NEW `.npy` (absdiff 7.1e-15) ⟹ the
re-baseline was LOAD-BEARING (the strict gate would trip on the un-rebaselined
snap), not cosmetic. Three re-pointed test files (`TestSubtractiveDefinition`
→`array_equal`→`assert_array_almost_equal_nulp(256)` + boundary STRICT;
`test_apply_equals…`→`test_pure_L_plus_C_recovers_loss_action_het` with
composite==loss_action byte-exact + affine ULP; `TestResolutionADifferent…`
→`TestPureLIsLossActionAtZeroSigma` array_equal vs `loss_action(0)`) — all
structurally grounded.

**CLAIM-3 (C1 σ-freedom teeth + Mode-11).** Mutation-verify BOTH: (a) sentinel
on `loss_rep.streaming_action` FIRES (hits=1) when `L.apply` runs → C1 reaches
the rewired matvec leaf (Mode-11 — the leaf has zero callers, sweep routes
around it); (b) a σ-re-reading stub (`loss_action(σ_t)` not `loss_action(0)`)
makes `L.apply` σ-dependent → outputs DIFFER by O(1) maxdiff ~12 → C1's
`array_equal` reddens. The shipped teeth test asserts the leaking stub differs;
both confirmed all 3 geoms.

**Gates.** pyright full-tree 2297/19 = baseline (the AUTHORITATIVE oracle; the
baseline-WORKTREE pyright is mis-rooted (L-027 cross-tree-config artifact) — `numpy` unresolved → renders
the SAME #226 family with `Unknown` types, so the message-multiset diff is
noise; confirm instead that NO live diagnostic references `streaming_action`/
`_zero_sigma_for`/the new apply). 0 net-new `# type:ignore`. Broad regression
`-O`: 7 reds, ALL reconciled PRE-EXISTING on the baseline worktree (5 SPHERE
#250 = 3 `test_bc_extraction_matvec[*-SPH]` ~1e15 ULP L-034-stale + 2
`test_streaming_operator` T4c sphere; 2 mu_y #232). ⚠ the spec's route-around
`-k "not (sphere_1g/2g_apply)"` named only 2 of the 5 SPHERE reds — the 3
`test_bc_extraction_matvec` SPH are also #250-family (run WITHOUT `-k` +
reconcile all 7 vs baseline = stronger). CLAIM-5: `scattering.py`/`fission.py`/
`orpheus/transport/` byte-untouched (empty diff); `loss_action` body unedited
(only `streaming_action`/`_transpose`/`_zero_sigma_for` ADDED);
`InvertibleOperator.apply`/`.solve` NOT in diff. NO blocker, NO false-green.

---

## L-050 -- a singledispatch alias rename (`apply`→`_apply_impl`, `else: apply=_apply_impl`): prove runtime bit-id via `Cls.__dict__['apply'] is Cls.__dict__['_apply_impl']`, and a removed `NoReturn`-poisoned return UNMASKS every latent downstream pyright error

The change: `FissionOperator`/`ScatteringOperator.apply` dispatch on input CARRIER
type → DISTINCT output carrier (heteromorphic, not endomorphism). S8c renamed the
`@singledispatchmethod` dispatcher `apply`→`_apply_impl` (base `-> "Any"`, was no
annotation ⇒ pyright inferred `NoReturn`), kept all `.register` arms at natural
indent, and added `if TYPE_CHECKING: @overload def apply(...)->Carrier; def apply(self,x:Any,/)->Any else: apply=_apply_impl`.

**Runtime bit-identity is BY CONSTRUCTION + the alias-identity proof.** `apply` IS
`_apply_impl` at runtime — prove via `Cls.__dict__['apply'] is Cls.__dict__['_apply_impl']`
→ True (do NOT use `Cls.apply is Cls._apply_impl` — the singledispatchmethod
DESCRIPTOR returns a fresh `_singledispatchmethod_get` wrapper each class-attr access
→ False, a red herring). The `TimedFullField` arm's `self.apply(...)` still routes the
SAME dispatcher. Confirmed empirically: 111 operator-suite + C6 PASS, 77 Section-D
MMS backstops PASS (2 pre-existing xfails #195/#252) — bit-identity-sensitive
(convergence rates + ERR-026 catches would break on any runtime change).

**Mode-11 (rewired-path reached).** C6 `test_c6_apply_dispatch_parity` calls the
PUBLIC `apply`; alias-identity guarantees reach, but PROVE it: register a sentinel
arm on `Cls.__dict__['apply'].dispatcher` (3.14 API: `.dispatcher.dispatch(Carrier)`
to grab orig, `.dispatcher.register(Carrier, wrapper)`) then call `F.apply(phi)` →
sentinel fires count=1 + returns the right type. Mode-8 OK (`pytest.fail` not bare
assert; passed under `-O`).

**Static C6 gate has TEETH (mutation-verified).** `_c6_static_typing_pins` (no
`test_` prefix → never collected; pyright-only) carries `assert_type(F.apply(phi),
ScalarSourceSink)` per carrier. Mutate one overload (`ScalarFlux→ScalarSourceSink`
⇒ `→AngularSourceSink`) → `npx pyright <testfile>` reddens EXACTLY that
`assert_type` line (`reportAssertTypeFailure`) → revert (L28: edit-revert, NOT git
stash; verify exact via grep+sha256). Clean file shows only 3 pre-existing
`BC.reflective` `reportAttributeAccessIssue` (enum-stub quirk, not S8c).

**The −15 net pyright = −19 disappeared + 4 net-new (RECONCILE EXACTLY).** Method:
back up S8c files (sha256), `git show HEAD:<f> > <f>` to restore the clean
pre-change baseline (S8c uncommitted ⇒ HEAD IS baseline; NOT git stash per L28),
full `npx pyright --outputjson` baseline + current, diff on `(file,rule,msg)` key
(robust to line shifts) THEN confirm with PER-FILE before/after counts. ⚠ the
`(file,rule,msg)` global key gives FALSE net-new when a message's TYPE-RENDERING
text shifts at the same logical error (#257 S8c:
`test_krylov_curvilinear_precond_safety.py` L174 `gains` arg showed as both −1 and
+1 — SAME error, per-file count 4==4 = net ZERO; the `LinearOperator[V@Krylov…]`
render changed). The REAL net-new = +3, ALL in the standalone capture SCRIPT
`tests/sn/_fixtures/wave_t_t3/_capture_pre_t3_snapshots.py` (L191 `aniso.values`
×2 + L204 `np.savez allow_pickle`).

**Root cause = NoReturn→unreachable SUPPRESSION lift (PRE-EXISTING LATENT, not a
regression).** Baseline `apply` (no annotation) inferred `NoReturn`; pyright treats
statements AFTER a `Never`-returning call as UNREACHABLE → suppresses ALL
downstream diagnostics. So line-175 `out=p1_op.apply(psi)` poisoned the whole
`main()` body below it → the LATENT errors at L191/L204 were hidden. S8c's `Any`
base makes the body reachable again → the pre-existing under-typing surfaces.
CLASSIFY: pre-existing-latent, ZERO runtime defect — `build_aniso_source` declares
`np.ndarray | AngularSourceSink | None` but at runtime (scattering_order=1, non-None
psi) returns `AngularSourceSink` which HAS `.values`; `np.savez allow_pickle` is a
numpy-stub `**kwargs` quirk. File is a one-shot capture script (leading `_`, in
`_fixtures/`, `def main()`+`__main__`, no `test_`, NOT pytest-collected). NOT a
blocker; optional follow-up = tighten `build_aniso_source` return or annotate the
script. RULE: removing a `NoReturn`-poisoned dispatcher return UNMASKS every
pre-existing latent error in code downstream of the FIRST poisoned call — expect
net-new ≠ (per-file delta in the two changed files); reconcile globally + classify
each unmasked error as latent-vs-regression by checking it's an under-typed
accessor with a concrete correct runtime type, NOT a real defect.

**cast/ignore hygiene.** 0 new `# type:ignore` (scattering's lone grep hit @645 is
PROSE in the S6 docstring saying "NOT a type:ignore"). 3 honest casts: 1 production
(`scattering.py:1237` `cast("AngularFlux|HarmonicMomentField", psi.bulk)` — runtime
`psi.bulk` is `AngularFlux`, both union members dispatch to `AngularSourceSink`,
verified live) + 2 test sites (under-typed `integrate_angular()→object` / `state.bulk
→BulkField`). NO blocker, NO false-green, NO ERR.

---

## L-051 -- a mechanical API-migration rewire (deleted class → new face) is bit-id-VERIFIABLE not bit-id-ASSUMED: recompute the OLD einsum on a structurally-independent table; and brief-named symbols/files are CLAIMS (two phantoms here)

Task: `MomentProjection`/`HarmonicMomentReconstruction` (orpheus/numerics/projection.py)
DELETED; rewire tests to `quad.angular_frame(L).analysis.apply` / `.reconstruction.apply`
(or `op.frame.analysis` where a `ScatteringOperator` is in scope). The new faces delegate
to `SphericalHarmonicBasis.analyze`/`reconstruct`; brief CLAIMED 0-ULP bit-id.

**Don't ASSUME the brief's bit-id claim — PROVE it before trusting the unchanged asserts.**
The 3 `test_scattering_operator.py` sites carry `np.testing.assert_array_equal` against a
FROZEN `.npz` snapshot captured by the OLD `MomentProjection`. For those to pass unchanged,
`op.frame.analysis.apply == old M.apply` must be BYTE-identical, not just close. Verified
in-process (10-line script): `frame.table == quad.spherical_harmonics(L)` (np.array_equal),
`frame.analysis.apply(psi) == np.einsum("n,nlm,n...->lm...", w, Y, psi)`, and
`frame.reconstruction.apply(c) == np.einsum("nlm,l,lm...->n...", Y, 2l+1, c)` — all True.
The recomputed einsum is the STRUCTURALLY-INDEPENDENT reference (hand-written contraction,
not the production face) — this is the bit-id leg, NOT old-vs-new ULP. Frozen snapshot needs
NO regen: the value is byte-identical (L-049 inverse — here byte-id IS preserved, so the
re-baseline question doesn't arise).

**Two PHANTOMS in the brief — confirm every named symbol/file before editing.** (a)
`SphericalHarmonicBasis.from_quadrature(quad, L).values` (in the snapshot generator's
defensive block) does NOT exist — the brief flagged it; replaced with `quad.angular_frame
(L).analysis.apply`. (b) brief item 6 said "find `_s6_stub_plugin.py` (grep for it)" — it
does NOT exist anywhere (`find` + text-grep for the name both empty); the user's original
grep had conflated it with the `from ...projection import` matches. NOT every brief-named
artifact is real — a `find`/grep confirmation is one cheap call and prevents a fabricated edit.

**Scope discipline on a concurrent-edit task.** User was editing 4 sibling test files
(test_frame/test_spherical_harmonic_space/test_scattering_kernel_crosscheck/
test_projection_operators) — touched NONE. After rewiring, byte-compiled both no-test
generator SCRIPTS (`py_compile`) and exercised the rewired `_capture_legendre_moments`
helper end-to-end (shape `(L+1,2L+1,ng,nx,ny)`) since scripts aren't pytest-collected — a
broken script import is a latent breakage no test run would surface. 100 tests PASS under
`-O`. Pure mechanical rewire, no claim-pushback → no vv-principles/ERR update fires.

---

## L-052 -- a Hilbert-adjoint-via-metric-composition (`A.H = G_dom⁻¹·Aᵀ·G_cod`) is VERIFIABLE by a dense-matrix transpose + the DEFINING inner-product law; the "weight-free transpose" choice is provable, not faith

Frame carve Phase D added `R.H` (reconstruction Hilbert adjoint) to the discrete `Frame`
via a new `reconstruct_transpose` (`einsum("nlm,l,n...->lm...")`, weight-free) + a capability
flip so the generic `_AdjointOperator` composes `R.H = G_basis⁻¹·Rᵀ·G_measure`. VERDICT
SUPPORTED, math correct, all gates have teeth. The reusable adversarial recipe for a
normal↔adjoint operator-algebra change:

1. **Re-derive the adjoint identity by hand FIRST, then numerically prove it.** `R[n,(ℓm)]
   =(2ℓ+1)Yₗᵐ`; the matrix transpose `Rᵀ[(ℓm),n]=(2ℓ+1)Yₗᵐ` is weight-free BY DEFINITION
   (a representation transpose carries NO metric). Compose with metrics: `R* = g_C⁻¹·Rᵀ·w`
   `= ((2ℓ+1)/4π)·(2ℓ+1)·Σ_n w_n Yₗᵐ v_n = (2ℓ+1)²/4π·Σ w_n Y v`. The "weight-free
   reconstruct_transpose" choice is CORRECT, not a missing-factor bug: it is ASYMMETRIC with
   `analyze_transpose` (which DOES carry `w_n`) precisely because each transpose mirrors its
   OWN forward — `analyze` bakes `w_n` in, `reconstruct` does not. A spurious `w_n` in
   `reconstruct_transpose` would give `R*` a `w_n²` and break `⟨Rc,v⟩_W = ⟨c,R*v⟩_{g_C}`.
2. **The structurally-independent reference is a DENSE matrix built by LOOPS, transposed
   directly, composed with metrics by hand** — shares zero code with production's fused
   einsums. `R^T` agreed at 0 ULP, `R.H` at 0 ULP, the closed-form `(2ℓ+1)²/4π·Σw_n Y v`
   target at ~3e-15 (FP non-assoc on the reduction). The DEFINING law `⟨Rc,v⟩_W=⟨c,R*v⟩_{g_C}`
   (Riesz) is the strongest pin — calls `R.apply`/`R.H` on both sides but asserts their
   ALGEBRAIC consistency, NOT circular (it's the adjoint definition itself).
3. **Teeth: 4 mutations, all RED via in-process monkeypatch under `-O`.** drop `(2ℓ+1)` /
   bake a per-node factor / reverse the factor array (`[::-1]`) → all 3 new tests RED;
   wrong GRAM POWER (`metric_per_ell` squared, build a FRESH frame so it flows into the
   space) → both reconstruction-adjoint tests + the analysis `R.H` test RED. Restore → green.
4. **Mode-11 cleared by a sentinel that COUNTS entries into the rewired readers:**
   `R.H.apply(v)` calls `_FrameReconstruction.apply_transpose` ×1 → `basis.reconstruct_transpose`
   ×1 (the new path IS on the gate's call graph); reverting the capability to `{CAP_APPLY}`
   makes `R.H` raise `MissingCapability` (the cap flip is the load-bearing enabler). The
   capability-assert test has teeth: under pytest it REDs with "Extra items in the right set:
   'apply_transpose'" when the face is reverted to APPLY-only (inject the class-attr mutation
   via a `PYTHONPATH=/tmp` pytest-plugin `pytest_configure`, NOT a standalone script).
5. **⚠ METHODOLOGY TRAP I hit (L-010 self-application): a bare `assert` in MY OWN `python -O`
   probe SCRIPT is STRIPPED** — my throwaway `assert rec.capabilities == ...` printed "PASSED"
   while the values were visibly unequal (Mode-8 in the PROBE, not the test). A capability/value
   teeth-check MUST run through PYTEST (assertion-rewriter active in `tests/`) or use
   `np.testing.assert_*` / explicit `if x!=y: raise` in the script — NEVER a bare `assert`
   under `-O`. The test itself was fine; my probe lied.
6. **The bit-faithful-reference choice (rtol=1e-14 per-term fold) is PRINCIPLED, provable.**
   The test's `einsum("nlm,n->lm", Y*f, v)` (2-operand, pre-scaled table) is BIT-IDENTICAL to
   production's 3-operand `einsum("nlm,l,n...->lm...")` (0 ULP, `array_equal` True — rtol=1e-14
   is generous). The REJECTED post-scaled form `f·(S0ᵀv)` drifts 112 ULP (docstring said ~2;
   direction right, magnitude under-stated) because the Σ then runs at ×f-larger magnitude. A
   third independent dense matmul also agrees 0 ULP → the value is right, not just close to one
   reference. Choosing the per-term fold over post-scaling is a bit-FAITHFULNESS choice, NOT a
   tolerance relaxation. Cross-ref [[lessons-L011]] (per-term-fold = structural independence),
   [[lessons-L051]] (recompute the einsum on a structurally-independent table for bit-id).

## L-053 -- migration-review of a `.solve -> .inverse().apply` rewire: the keystone pins the WRAPPER, the migration rewires the LOOP; and the strong catcher may be slow-deselected

(Taxonomy re-evaluation 2026-07-01, branch refactor/inverse-as-operator @ 69ed531 -- F4/W5/W2
legs; distilled by the qa leg, persisted by the main agent.)

1. **Wrapper-surface vs loop-surface.** A keystone gate pinning `inverse().apply(b) == solve(b)`
   covers only the single-call WRAPPER delegation -- bit-identical BY CONSTRUCTION (both sides
   share the same, possibly-buggy, sweep), hence robust to / orthogonal to sweep-correctness
   bugs (#282). It does NOT cover the ITERATION-LOOP per-iterate seed threading -- the surface
   a `.solve -> .inverse().apply(rhs, initial_guess=psi_prev)` migration actually rewires.
2. **Test the loop surface by simulating the exact regression under the ACTUAL run config.**
   Patch the seed-threading helper to drop `initial_guess` in-process and run the canonical
   `-m "not slow"` invocation -- do NOT trust a plan's "test_X must red" claim: the strong
   end-to-end catcher here was `@pytest.mark.slow` and thus DESELECTED (a slow-marker sibling
   of Mode-8: a gate that cannot fire under the run config). The het-2G sphere seed-drop
   (|dk| = 3.46e-2) reddened only a fragile 1G monotone margin under `-m "not slow"`.
3. **A seed-insensitive geometry makes its seeded VALUE gate vacuous for seed-drop detection**
   (cylinder telescopes -- the seed cancels identically). The Mode-11 path-spy on the
   `initial_guess` threading is the load-bearing guarantee there, not any value assertion.
4. **A "fold-config identity" inverse realized as a SCHEDULE** (`B_lower` = octant-group edge
   set, not an algebra operator) **has NO exact round-trip pin** -- no forward `.apply` exists
   to invert against, and the converged-SI-equivalence fallback is Mode-9-DEGENERATE (the
   fixed point is splitting-invariant by construction: it cannot distinguish the fold, or even
   G-S from Jacobi). Minting the forward matvec of the SAME restricted operator (reify
   `M = (L+C-B_lower)`) is the only way to make the round-trip exact.

---

## L-054 -- triaging audit-MISSING `catches("ERR-NNN")`: grep the production RAISE SITE / invariant-enforcement FIRST -- three outcomes, and only one is "add the marker"

(Metadata-only marker-patch task 2026-07-11, branch refactor/sn-walk-unification; 9 MISSING
ERRs -> 5 tagged, 4 reported NO CATCHER.)

The audit lists an ERR as MISSING when no test carries its `catches` marker. That is NOT
automatically "add a marker" -- before tagging OR reporting a gap, grep whether the cataloged
bug is even *reachable*: find the production `raise <ErrorClass>` site (or the enforcement of
the invariant the ERR broke). Four things it discriminates:

1. **Genuine catcher present -> tag it, mutation-verified.** The exact bug re-introduced reds
   THIS test (L-007). Verified empirically for ERR-020 (bit-identity `np.all(vol==vol[0])`;
   edge-derived `**3` round-trip reds it), ERR-031 (positional arg-swap -> the swapped radii
   `[2.0,0.1]` trip the strictly-increasing guard -> `ValueError`), ERR-040 (tangential
   ordinate admitted to a selector -> `test_axis_aligned_ordinates_excluded_from_both_selectors`
   reds). Probe via a `/tmp` script with EXPLICIT boolean prints, never a bare `assert` (L-052).
2. **Catalog's "L0 test" names a RETIRED test class -> the marker didn't migrate.** ERR-020's
   entry named `TestZoneSubdivision::test_equal_volume_*` -- retired in Phase F, MOVED to
   `test_structured_geometry.py` (`test_equal_volume_{cyl,sph}_invariant`), and the marker was
   left behind in a now-STALE comment claiming the decorators "stay". Retirement-means-test-
   migration (L-022 family) applies to MARKERS too: re-tag the successor asserting the SAME
   invariant. This is the usual cause of a MISSING whose catalog "L0 test" still reads plausible.
3. **Typed error defined+exported but NEVER raised = dead scaffolding -> genuine unbuilt-invariant
   gap, report NO CATCHER.** ERR-041/045/047: the error classes ship + export + have a catalog
   entry ("TYPE SHIPPED Wave 3 / catching test ships Wave 7"), but `grep -rn "raise <Class>"
   orpheus/` is EMPTY and the `assert_*` invariant is a no-op default (no concrete override).
   The promissory "Wave 7 catching test" was never built; nothing can red on a recurrence.
   Do NOT invent a marker -- the MISSING is truthful.
4. **`assert_X` delegates to a WEAKER sibling invariant -> it catches the sibling's bug, NOT its
   own.** ERR-042 (`assert_geometry_map_measure_preserving`) `self.assert_is_involutive(quad)`
   and ASSUMES weight-symmetry "by construction" -- so it reds only on non-involution (ERR-044),
   never on the weight-measure drift ERR-042 documents; no Q4.x quadrature test checks
   `weights[ref]==weights` either. Tagging its test `catches("ERR-042")` is the exact L-007
   blind-marker trap (reds on a DIFFERENT class). Report NO CATCHER; the method's name over-
   claims its body (coding-elegance #20).

The pushback rationale for outcomes 3-4 is already covered by the vv-principles `catches`-marker
directive (mutation-verify the EXACT bug reds THIS test) -- no new anti-pattern; this is its
audit-triage application. Marker-only edits: confirm green under canonical `-O` AND that the
`git status` dirt in do-not-touch files is PRE-EXISTING (grep your own ERR numbers out of their
diffs) -- a shared working tree makes every dirty file look like yours.

---

## L-055 -- adjudicating campaign-narration staleness (#304 class): FIX bar = "provably lies about CURRENT code", verified against grep+`gh` BEFORE ruling; two over-fix guards

A doc-hygiene pass over campaign narration (`Phase X` / `Wave Y` / `(#N step` tags in test
comments/docstrings) is a Cardinal-Rule-3 correctness task, NOT a numerical review -- so it does
NOT touch the vv-principles anti-patterns (those are about verification EVIDENCE). The
three-way rule: KEEP-current (open issue's genuinely pending work), KEEP-provenance (truthful
backward attribution / retirement records / plan-pointers -- keep even if the plan file is
gone), FIX-stale (narration that LIES about the present). Conservative default = KEEP.

**The FIX bar is "provably lies about CURRENT code" -- and "provably" means VERIFIED, not
inferred from the tag.** Before ruling a forward-looking deferral ("future X", "blocked on
Phase C", "Phase 2 will land", "once #N lands") stale, verify the *future* against reality:
(1) `grep` the named symbol/closure/wiring/workaround tree-wide -- the strongest FIX signal is a
"future" thing the code now SHIPS + WIRES (a callsite proves it, e.g. `_build_white_hebert_op`
calls `compute_P_ss_cylinder` at a specific line) OR a "pending" workaround that is now 0-hits
(`cast(LinearOperator` gone tree-wide); (2) `gh issue view N --json state` -- an OPEN issue whose
*named sub-phase derivation* landed while a DIFFERENT residue stays open (#112: "Phase A/C
derivation landed, 3-D-normalization/rank-N residue open") means the tag's specific claim can be
stale even though the issue is open. Landed => rewrite to present-tense truth (name the shipped
fn + the still-open residue). Ambiguous / partial-landing (a phase that shipped INFRASTRUCTURE
but the usable capability still raises `NotImplementedError` pending an open follow-on; a bare
campaign phase with no issue number and no confirmable landing) => KEEP + report as orphan-TODO,
do NOT guess.

**Two guards against over-fixing** (both bit me as tempting-but-wrong FIX targets):
1. **A stale line inside a RUNTIME STRING is behavioral -> KEEP even when its text is stale.** A
   `description="... Wave 8 will switch ..."` dataclass field, an f-string diagnostic, an assert
   message -- these are data the code may write to a snapshot / test-ID / error, so editing them
   is a behavioral change. The "never touch runtime strings" constraint PROTECTS you here: the
   one genuinely-stale line in tests/geometry (`_generate_bc_equivalence_snapshots.py:159`) was a
   `description=` field -> untouchable despite the module docstring itself confirming Wave 8
   landed.
2. **A load-bearing-gate "failure here HALTs Phase X" banner is a characterization RECORD, not a
   pending-work lie -> KEEP after Phase X lands.** It states the test's structural importance
   ("this is THE exactness gate; its failure invalidates the whole closure"), which is durable;
   the `Phase X` is provenance, usually paired with a plan-pointer. Rewriting it churns truthful
   history into design vocabulary (the task explicitly forbids that).

Mechanics: comment/docstring-only, ZERO behavioral change; verify `pytest --collect-only` clean
after; a shared working tree means `git diff --stat` shows OTHER agents' production-file edits --
diff ONLY your own touched files to confirm your edits are comment-only (cross-ref L-054's
shared-tree note). (#304 surface-2, 2026-07-22: 277 hits in scope, 3 files FIXed -- a
future-closure-now-shipped, a workaround-now-0-hits, a Phase-2-constructors-now-shipped; the rest
KEEP.)

---

## L-056 -- reviewing a skill->Sphinx DISTILLATION: verify code-anchored specifics against CODE, never against the skill source (the source carries stale specifics that propagate, and the build gate is blind to them)

When a Sphinx page is authored as a faithful distillation of a `.claude/skills/*` doctrine
(e.g. `verification/principles.rst` from `vv-principles/{SKILL,reference}.md`), the DOCTRINE
(ladder/pillars/claim-layers/modes/anti-patterns) is almost always faithful -- reading it against
the preloaded skill confirms mechanism/instance/defense with no inversion (verified the whole
modes-7..12 highest-risk block clean this way in one pass). The yield is entirely in the
**code-anchored SPECIFICS the skill states but the build never checks**: module paths, an
"evaluated in mpmath" vs `scipy.optimize.brentq`-double impl detail, a test-count, an ERR
war-story's composition. Two structural reasons the build gate misses these:

1. **Python-domain roles are not `-W`-gated.** A `:mod:`/`:class:`/`:func:` pointing at a
   NON-EXISTENT target renders as plain text with NO warning unless the build runs `-n` (nitpicky)
   -- so an "exit 0 -W" gate is ZERO evidence a `:mod:` resolves. Verify every code-pointer by
   filesystem/`find`, AND grep the WHOLE corpus for the canonical spelling: the OUTLIER count is
   the bug (caught `orpheus.derivations.continuous.peierls` used 1x on the reviewed page vs
   `...peierls_nystrom` 240x everywhere else -- the bare form is a dead module).
2. **The skill source is not the code and is not build-gated, so its stale specifics propagate
   verbatim into the corpus.** The reviewed page inherited BOTH a dead module path AND a wrong
   impl detail from `reference.md`/`SKILL.md` -- both said the same wrong thing, so cross-checking
   the page against the skill would have PASSED them. "Code outranks doc" means CODE, not the
   trusted skill twin: confirm `mpmath`-vs-`scipy`, `brentq`-vs-`findroot`, `dtype=float`-vs-`mp.mpf`
   by reading the derivation, and read the CONSUMING test's docstring (it often states the truth
   the doctrine page fumbled -- here "double-precision transfer-matrix back-substitution").

A "worked example" whose stated purpose is "every coordinate is TRUE for this case" is a MUST-FIX
magnet: check each coordinate independently (claim-layer / ladder / pillar / tier / operator-form)
-- the classification can be right (semi-analytical, T2, `operator_form=="diffusion"` all held)
while ONE parenthetical mechanism token is false. None of these pushbacks is a vv-principles
anti-pattern (they are doc-accuracy, not evidence-reasoning) -> no skill anti-pattern addition;
this is the distillation-review application of Cardinal-Rule-3 + the retirement-audit "grep docs,
Python-roles silent under `-W`" rule. (#10 stage V5 principles.rst review, 2026-07-23: 2 MUST-FIX
[dead `:mod:` peierls path; "evaluated in mpmath" on a brentq/double reference] + 1 SHOULD-FIX
[carried "twenty one-group tests" where 3 sources say "20 passing tests", a mix]; the ~40-claim
doctrine body otherwise faithful.)

---

## L-057 -- reviewing a results-COMPILATION page (de-freeze + evidence-map + run-book): the count-de-freeze is CERTIFIABLE by live collect, a retitle can beat the test's OWN stale docstring, and a run-book that cites config for "operational detail" can point at a note that CONTRADICTS its headline

#231 task-#10 stage V6 = authoring `docs/theory/verification/summary.rst` (the V&V-part
results compilation) + de-freezing 4 frozen test-counts across the per-method chapters.
VERDICT PASS; the ~50-claim page was faithful end-to-end. Three reusable techniques:

1. **A count-DE-FREEZE is certifiable, not taken on faith.** When a diff removes a frozen
   "N tests across M files" and replaces it with "the auto-generated matrix carries the live
   counts", PROVE the old number was stale by `pytest <dir> --collect-only -q | tail -1`:
   CP 106/6→**154/11**, MC 55→**57**, MOC 102→**104** — all three old literals genuinely
   lied, so the de-freeze is warranted (not an invented rationale). The NEW page states NO
   count (counts-de-freeze doctrine), so a brief's own count claim ("48 rows") is NOT a
   doc-truth defect even when the live table is **47** — but report the delta so the parent
   doesn't propagate the wrong number. Structural counts that SURVIVE (the 27-case CP grid =
   3×3×3, the 4/3/4/1/5 cross-method list lengths) ARE verifiable: `.venv/bin/python -c` the
   list lengths + `ADAPTERS_BY_NAME` (6) at runtime, don't eyeball.
2. **A doc-RETITLE can be MORE accurate than the test's own name/docstring — verify against
   the live ASSERTION body.** V6 retitled the SN property "Flux symmetry"→"Flux flatness";
   `tests/sn/primitives/test_properties.py::test_flux_symmetry` is NAMED "symmetry" and its
   docstring says "must be symmetric about the center", but the LIVE assertion is
   `assert_allclose(flux, flux[0], rtol=1e-6)` ("homogeneous slab flux is exactly flat"). The
   retitle correctly describes the assertion, not the stale name. So a retitle-faithfulness
   check reads the `assert`, NEVER the test name or its (possibly-stale) docstring — same
   "code outranks doc" as L-056, applied to test-name vs test-body. (Diffusion vacuum→Marshak
   was the same shape: the doc states the ASSERTED framing `J⁻=0 @ 1e-12·scale AND
   boundary-cell flux>0`, matching the body, not the old "flux is small" scaffold.)
3. **A run-book that cites config for "operational detail" can point at a note that
   contradicts its own headline.** V6's run-book calls `python -O -m pytest -m "not slow"`
   "the pre-merge gate: the full tree ... single-process" and says the `[test]` extra's
   pyproject notes "carry the operational detail" — but those notes say "The SN suite OOMs
   when run as ONE single-process invocation ... **NO whole-tree single-process run**"
   (per-tier is the memory-safe default). A whole-REPO single-process run executes the whole
   SN tree in one process = exactly what the cited note warns against. RECONCILABLE (the
   pyproject note is inner-loop SN-iteration memory advice; the pre-merge gate genuinely IS
   the full-tree SERIAL run — `reference_test_execution_env` memory: completes 6391/0 in
   ~52 min, xdist UNSTABLE so serial is canonical), and the xdist "within-tier" statement IS
   faithful to pyproject — so it is a NIT (surface the pre-merge-gate vs inner-loop-per-tier
   distinction), not a falsehood. The lesson: when a run-book delegates to a config file for
   detail, READ that file and check the delegated-to text doesn't read as contradicting the
   delegating headline. Not a vv anti-pattern (doc internal-consistency, #231 prime directive)
   → qa-lessons only, no skill edit. Everything else — 6 evidence-map anchors resolving to
   CLAIMED content (not just resolving), the 8-case SN MMS ladder→files, the Mode-12
   homogeneous K=A⁻¹F matrix-object gate, the `compute_kinf_*`-vs-`kinf_homogeneous` footnote,
   the 10 matrix.rst headings, Peierls `precision_digits=30`, the `sentinel` marker "run
   WITHOUT -O", `generate_rst` runnable as `-m` + `reference_values` pkgutil auto-discovery —
   verified faithful. (#231 #10 V6, 2026-07-23.)

---

## L-058 -- a "k is designed-green / functional-blind to the mutation class" (Mode-12) claim is VERIFIABLE by running the mutation — a leaf-transpose-DROP is NOT similar to forward, so its k SHIFTS

#276 A4 phase-gate: the adjoint φ* certification. The test docstrings + the standalone
NOTE claimed "F†=F leaves k EXACTLY equal / k is designed-green on the entire adjoint
mutation class (eig(A†)=eig(A))" and used that to motivate the P1.3-k-vs-P1.4-spectrum
split. **The claim is FALSE and self-contradicted by a passing sibling test.** Verified two
ways:

1. **Direct closed-form check (10-line numpy):** under F†→F on ∞-medium 4G, the daggered
   operator becomes `(Aᵀ)⁻¹F` (F NOT transposed), whose dominant eig = `χ·A⁻¹νΣf = 0.153`
   vs forward `νΣf·A⁻¹χ = 1.488` — k SHIFTS 1.488→0.153 (|Δk|=1.33). The char poly is
   `det(A−Fᵀ/k) ≠ det(A−F/k)` for asymmetric A / χ∦νΣf. Only the CORRECT adjoint
   `(Aᵀ,Fᵀ)` is similar to forward `(A,F)`; a leaf-DROP mutation `(Aᵀ,F)` is a NON-transpose
   operator and its k is unconstrained.
2. **The self-contradiction tell:** `TestP13Mutations.test_fission_role_swap_shifts_k`
   applies the SAME F†=F on `_infinite_medium("4g")` and asserts `|Δk|>1e-6` — and PASSES.
   So the note's "the ∞ row stays GREEN under F†=F, which is why the k-tooth had to ride a
   shifted-k regime" is refuted by the very ∞-medium k-tooth it describes.

**The CORRECT Mode-12 framing** (the fix): `k_adj==k_fwd` is automatic FOR THE CORRECTLY-
BUILT adjoint (eig(A†)=eig(A)), so it confirms the eigenVALUE but NOT the adjoint FLUX SHAPE
(a right-k/wrong-ψ* solver — forward φ, or the νΣf degeneracy — passes the k-legs). THAT is
why the vector gates (spectrum, biorthogonality) are needed. The leaf-transpose-DROP
mutations DO shift k, so the k-legs carry real teeth too; the blind spot is the eigenVECTOR
identity, not "the machinery." **Behavioral rule: NEVER accept a "this functional is blind
to this mutation class" narrative by inspection — RUN the mutation (or the closed-form eig).
A k-blindness claim is right ONLY for mutations that keep the operator a valid transpose;
leaf-drops break the transpose and shift k.** This is the mirror of the #226 step-5b
OVERCLAIM (Mode-12): step-5b claimed a gate CATCHES when it's blind; this UNDERCLAIMS
(claims blind when it has teeth) — same defense (run it, don't narrate it). No coverage gap
(the genuine flux-shape blind spot IS closed by P1.4/P1.5); the finding is a wrong WHY =
should-fix (CR3), not a blocker.

**Skill refinement flagged (read-only task, so surfaced not applied):** the `vv-principles`
Mode-12 "Live application" text carries the same overstatement verbatim — "#276 A4 … 'k*
matches k' carries ZERO mutation coverage on the adjoint machinery." Sharpen "ZERO mutation
coverage on the adjoint machinery" → "cannot confirm the adjoint eigenVECTOR/flux-shape; the
leaf-transpose-drop mutations still shift k." When a review pushes back on the skill's OWN
example, flag the skill edit as a finding under a read-only constraint (main agent applies) —
honors both the read-only task instruction and the self-improvement trigger.

**Other durable micro-findings this review:** (a) an angle-flat-blindness justification is
checkable — `Σχ/W==Σwχ/W when angle-flat` holds IFF W==N; false for lebedev (N=110, W=4π:
8.75 vs 1.0) — but a conservative angle-varying `require` + an in-test wrong-spelling
discriminator make it harmless (doc-precision note, not a teeth failure). (b) `foundation`
stacked on a method under a module `pytestmark=l1` → `PytestUnknownMarkWarning: conflicting
V&V level markers … using 'l1'`; the intended foundation level is silently dropped — a
level-marker hygiene nit distinct from the L-007 foundation+verifies conflation. (c) a
k-only geometry leg with no closed-form (the coupled sphere daggered posing) is honest-scope:
it rides upstream-verified transpose machinery (#280/#310) + k-equality; flag the missing
vector-shape check as a scope boundary, don't credit k as flux-shape evidence.

**A6/ch15 RE-REVIEW SHARPENING (2026-07-25): 0.153 is the 0-D PROXY, not the SN-solve
k-tooth's value (0.171) — the metric carries the angular weight.** The `1.488→0.153
(|Δk|=1.33)` in point 1 above is the 0-D char-poly `eig((Aᵀ)⁻¹F)` (angular-COLLAPSED). The
ACTUAL `TestP13Mutations.test_fission_role_swap_shifts_k` SN daggered SOLVE gives
**k_mut = 0.171 (|Δk|=1.317)**, NOT 0.153 — reproduced 4× via the test module's own helpers,
stable + converged. WHY they differ: `.H`'s metric `G = V·w_n` carries the ANGULAR weight, so
the mutated (non-transpose) fission `F.H_mut = G⁻¹FG` is angularly non-trivial (the mutated
adjoint mode is 21% non-flat across ordinates); the 0-D reduction that yields 0.153 collapses
that angular structure. The QUALITATIVE claim (leaf-drop → k moves O(1)) is right; only the
MAGNITUDE 0.153 is the wrong (0-D) number for the SN k-tooth. **A6/ch15 propagated the 0-D
proxy everywhere it describes the SN-solve k-tooth:** `docs/theory/methods/sn/adjoint.rst:941`,
`docs/theory/verification/sn.rst:4893`, cert-test docstrings `:327`/`:359`, AND vv-principles
`SKILL.md:135` all cite `1.488→0.153` for the SN-solve mutation → should be `0.171` (MUST-FIX;
the number is NEVER asserted — the test only checks `|Δk|>1e-6` — so no gate fails, but the
cited magnitude is a plausible-substitution error: a real eig of a related operator).
**Rule (sharpens the L-058 verify-by-running rule): a cited mutation-MAGNITUDE for a
METRIC-adjoint SOLVE must be the full-solve value (RUN it), never the angular-collapsed 0-D
char-poly proxy — the metric conjugation on a MUTATED (non-transpose) operator is NOT
spectrum-preserving (only the CORRECT full-dagger is similar to forward), so `0-D ≠ SN-solve`
whenever the metric carries a reduced axis (`w_n`).**

---

## L-059 -- an accelerator's PRODUCTION rate can sit ABOVE the operator's spectral radius (a splitting/wall lag); certify the operator by building the FULLY-COUPLED matrix, don't credit "matrix says healthy" from a sibling-BC certificate

#2 DSA 3c rate tier. The reflective-BC gate (D12) split the rate claim: D11's Fourier
bound (ρ ≤ 0.2247c) runs VACUUM (production ρ_est ≈ 0.18 ≈ matrix ρ), while the fully-
reflective regime uses a LOOSER one-sided band (ρ ≤ 0.35) because production measures
0.28-0.31 — attributed to a "Jacobi wall lag" (the production splitting reads the iterate's
outgoing trace one iteration late), NOT a rate bug. The question: honest split or paper-over?

**The decisive check is to BUILD the fully-coupled matrix (error-iteration) operator and
read its ρ directly** — the operator ρ is the floor the splitting can approach; if the
matrix ρ ≪ production ρ_est, the gap IS a splitting artifact (fixable: wall ordering / trace
relaxation), not an operator/consistency bug. Confirmed here:
- production refl/VAC ρ_est = 0.279 vs its matrix certificate 0.154 (D2 report) — the lag is
  real and visible on a config that HAS a committed certificate.
- I built the refl/REFL matrix (both walls resolved IN-sweep by iterating the boundary
  partner fluxes to convergence, then composed with the production low-order via the test's
  own `_t_dsa`): ρ = 0.19-0.22 — HEALTHY. Production refl/refl 0.28-0.31 is the lag,
  confirmed. Split is HONEST.

**The evidence-completeness flag (IMPROVE, not blocker):** the committed instruments
(`_wd_sweep_matrix`, D2 Part C, rate-report Part D) only certify refl/**VAC** and vac/vac —
the `_wd_sweep_matrix` hardcodes a vacuum right wall (its `bc[0]=="reflective"` branch is
DEAD in every test; L-016). So the docstring's "the matrix certificates say the operator is
healthy" for the refl/**REFL** regime that D12 actually gates rested on INFERENCE (mechanism
+ the refl/vac certificate) until the reviewer supplied the refl/refl matrix. When a rate-
split's honesty argument cites "the matrix certifies the operator", verify the certificate
covers the EXACT BC/config the runtime gate exercises — a sibling-BC certificate + "same
mechanism" is inference, and the fully-coupled matrix for the gated config is cheap to build.

**Behavioral rule:** to adjudicate "elevated production rate = splitting artifact vs rate
bug", build the fully-coupled operator matrix for the GATED config and compare ρ_matrix to
ρ_production. ρ_production > ρ_matrix ⟹ splitting lag (honest, characterize + file the
improvement); ρ_matrix ≈ 1 ⟹ genuine consistency failure. Never accept "the matrix says
healthy" when the committed matrix is a different BC than the gate runs. Companion to the
numerical-bug-signatures Sig-8 (unconverged-inner-solve masquerade) and Sig-9 (ρ-blind stop)
— a THIRD "the rate looks off" mechanism: a within-iteration boundary lag elevating the
splitting's ρ above the operator's.

---

## L-060 -- a transpose/adjoint RECIPROCITY gate is Mode-12 blind to a SYMMETRIC completion-drop; a symmetric-completion (E_out diagonal) inverse-fix is VERIFIABLE by dense (Aᵀ)⁻¹=(A⁻¹)ᵀ, and mutation-tested BOTH ways

#2 ERR-071 root fix (the honest full-space composite sweep inverse). The forward outflow-row
is the defect `streamed − ψ_out`, so `(L+C)⁻¹` emits `ψ_out = streamed − rhs_out` via a
post-march restore `buf[out_rows] −= seed[out_rows]`. The transpose half claims: `E_out` (the
restore) is a diagonal partial identity ⟹ symmetric ⟹ `(Aᵀ)⁻¹ = S_oldᵀ − E_out` = the SAME
one-site restore on the SAME forward-sense outflow selector in `solve_transpose`.

**The symmetry argument is directly VERIFIABLE (don't trust the prose).** Build the composite
flatten/unflatten, then dense `A` (apply on unit cols), `A⁻¹` (solve on cols), `(Aᵀ)⁻¹`
(solve_transpose on cols); check `A·A⁻¹=I`, `Aᵀ·(Aᵀ)⁻¹=I`, and — the crux — `(A⁻¹)ᵀ=(Aᵀ)⁻¹`.
Measured 1e-16 on slab/ld_slab/cyl_product INCLUDING the cyl free-DOF subspace (where
`A·A⁻¹` shows `1.0` on the 8 μ_r≈0 rows — A genuinely rank-deficient there, the honest
free-DOF pair, NOT a bug). This proof is independent of the bilinear reciprocity gate.

**⭐ THE Mode-12 FINDING — a reciprocity gate pins the transpose RELATIONSHIP, not
correctness.** `⟨A.solve q,p⟩=⟨q,A.solve_transpose p⟩` (the g3 gate, the named transpose
catcher) is satisfied by ANY genuine transpose pair `(S,Sᵀ)` — so mutation-test it BOTH ways:
- **MUT-T** (fix forward, undo ONLY the transpose completion: additively add `E_out` back →
  `S_oldᵀ`): g3 REDS at O(1) (measured 3.9–7.4%) — the asymmetric catch. ✓
- **MUT-BOTH** (empty the outflow selector → the true pre-fix HEAD state, BOTH halves dropped):
  g3 stays GREEN (5.7e-17) — `(S_old, S_oldᵀ)` IS a transpose pair, so reciprocity is blind.
The forward/symmetric half is caught ONLY by the one-sided identity gate
`test_sweep_inverse_identity.py` (`(L+C)∘(L+C)⁻¹≡I`, reds at 1.8 under the symmetric drop).
So the two gates are NON-REDUNDANT partners; neither alone covers ERR-071. The catalog +
streaming.py note attribute each half to its gate correctly (transpose→g3, forward→identity)
— NO overclaim, but flag for maintainers: never delete the identity gate on "g3 covers it".
**Behavioral rule:** when a fix touches BOTH `solve` and `solve_transpose` symmetrically,
a reciprocity gate is Mode-12 blind to the symmetric regression; require a one-sided
`A∘A⁻¹=I` companion and mutation-test MUT-T AND MUT-BOTH separately.

**Package A (P1-DSA d₁ arm) — clean, SUPPORTED.** (28b) `f₁=−(D/h)Δf₀+a·d₁@ρ=0`:
`moment1_update` bit-exact vs an independent (23f) recompute (max|Δ|=0); `_dh`=D/h,
`_a_coef`=a are the SAME arrays that build a_low/g_map (single-source, so transitively pinned
by the entry-for-entry reference-builder gate `test_dsa_low_order.py`). The (28b) COMBINATION
is not pinned entry-for-entry (no reference builder for the updates) but is end-to-end
CONSTRAINED by the S2-exactness anchor (angular space = span{1,μ}, so one correction must
annihilate the ℓ=1 gain): sign-flip / drop-a·d1 / 3× all blow n from 2 to 49 (Mode-10
mutation-verified — the term is constrained, not merely exercised). Anti-mint confirmed:
`angular_frame(1).table[:,1,1]==mu_x` bit-exact (a CALLED frame row, not a `w·μ` twin).
P0-forced tooth: healthy n=2, forced n=33 (large margin). so=0 is a pure rename of the P0
path (verified vs HEAD). Trace arm ℓ=0 by theorem (derivation reflecting row f₁=0; vacuum
inert). Cross-refs [[lessons-L058]] (Mode-12 verify-by-running), [[lessons-L024]] (prove teeth
by disabling), [[lessons-L016]] (product-quad needed to exercise the μ_r≈0 free-DOF branch).

# Archive — the grown-Krylov-composite restart truncation (#282) and the stale-frozen-reference triage

Moved COLD from `lessons.md` 2026-09-21 (verbatim, lines 394-439 then 749-800).
Digest successors: L15 and L22. Neither had an external record file — this archive IS the record.

## L15: To rule PRINCIPLED-vs-REGRESSION on a SEED / angular-CLOSURE re-pose, sweep ANGULAR order N at a FIXED fine mesh — NOT h at fixed N. And a carve that GROWS a Krylov composite must resize `restart` from the composite `to_flat`, not the bulk formula

Two durable points from the #282 route-(a) sphere ψ½ re-pose (OLD `edge_extrapolated_seed`
→ NEW `carlson_inward_sweep_from_source` direct march), diagnosed 2026-07-04 on
`refactor/sn-walk-unification` (the carve LANDED mid-investigation as commit `a29ab2d`
"2.5d d3"; pre-carve OLD = `5170f20` dormant-seed — re-run OLD via a worktree at `5170f20`).

**(a) The N-sweep is the discriminator; the h-sweep at fixed N MISLEADS.** A seed IS an
angular closure. Two seed treatments give DIFFERENT keff(h→0) at a FIXED angular order (their
low-N angular truncation differs) yet the SAME (h→0, N→∞) answer. Measured (het 1g fuel|mod
reflective sphere, GL8): NEW keff(h→0)=0.73825, OLD=0.73654 — a 1.7e-3 gap that does NOT
close under h-refinement (the task's "do they share keff(h→0)?" test FALSELY reads regression).
But sweeping N at fixed n=80: NEW & OLD AGREE to 1.5e-6 at N=32, 2.7e-6 at N=64 (both→~0.7368);
dd_regression 2g/3reg agrees to 8e-8 at N=64. ⇒ PRINCIPLED (the seed changes O(N) truncation,
not the converged value); the frozen N=8 snapshots (`sphere_2g_3reg_dd_n40`) are legitimate
§16.D re-baselines. This is L14/Mode-7 realized: MMS (≤linear-in-μ) is seed-blind, so MMS
O(h²) does NOT certify the seed — the N-sweep does. **Honest nuance:** at the test's N=8 the
OLD edge-extrap seed is actually CLOSER to the N-converged truth (0.7365 vs truth 0.7368) than
NEW (0.7382 overshoots); route-(a)'s justification is STRUCTURAL (an honest direct inverse —
cold residual 5e5→1e-11 for #200/#280), NOT angular accuracy. Sub-quadratic keff h-rate is the
PRE-EXISTING pole-cell O(h^1.4) (homogeneous-vacuum-sphere control isolates it, no interface)
propagated outward along characteristics — a single O(h) pole cell does NOT give O(h²) global;
shared by OLD, documented WONTFIX ([[curvilinear-tau-clamp-vs-pole-floor]]), not carve-introduced.
The n=5→10 near-coincidence (Δ8e-7) is a REAL coarse-mesh feature (persists at keff_tol=1e-12,
not iteration noise) that trips a fragile `diff_2<diff_1` ladder → robustify to n∈[10,20,40].

**(b) Krylov restart-truncation re-triggered by a grown composite (ERR-053 family).** Route-(a)
grew the within-group Krylov state from 2-block (bulk⊕trace) to 3-block (bulk⊕trace⊕
starting_direction seed), but `n_dof`→scipy-gmres `restart` is still the BULK formula
`N·ng·prod(spatial_shape)` (`solver.py:1511` eig, `:2599` fixed-src). The raveled composite
`to_flat` is LARGER (n=10: bulk 160 < composite 210 = +42 seed +8 trace), so restarted
GMRES(160) STAGNATES on the 210-dim augmented system (info=300, residual plateau, 868 s; keff
best-effort eventually right via the outer loop but WRONG=0.865 under a bounded outer cap).
Forcing restart=210 → info=0, keff=SI to 3e-10, 45× faster. So the seed block is NOT
intrinsically zero-metric-weight-unreducible (the task's hypothesis) — it just pushes the
ravel past the bulk-sized restart. NOT issue #200 (that's the IDENTITY precond, separate);
this is restart-sizing, route-(a)-introduced. OLD Krylov + the c=0.5 fixed-source path
converge clean (fit within one restart cycle); the stall needs poor conditioning (moderator
c=0.95 reflective eig) + the grown composite. Fix LANDED in the SAME `a29ab2d` d3 commit —
`n_dof=int(initial_guess.to_flat().size)` at both sites (eig + fixed-src); verified end-to-end
on HEAD (restart 210, info=0, k_SI≡k_Krylov 4.7e-11, 3.4 s). **General rule: any
operator-algebra carve that adds a block to a Krylov composite MUST re-derive restart/n_dof
from the composite dimension — grep every `restart=`/`n_dof=` against the ravel, not the bulk.**
Diagnostics `derivations/diagnostics/diag_282_{krylov_restart_truncation,sphere_repose_convergence}.py`;
probes `/Users/rodrigo/.claude/jobs/84fd66f8/tmp/probe_0{2,3,7,8}_*.py`.


## L22: A frozen reference is stale by MAGNITUDE CLASS, and the nulp count tells you neither

Triaging 9 "bit-identity" reds (2026-08-12, task 51: 5 Cartesian LS snapshot rows, 1 sphere
DD gate, 3 affine-carve sha256 arms) — all 9 were stale references, none a regression, in
two magnitude classes that the failure messages could not distinguish.

1. ⭐⭐ **A nulp count is uninterpretable in BOTH directions — always report
   `max|a−b| / max|b|` FIRST.** The received warning is "huge nulp near zero is nothing".
   The dual bites just as hard: `[M]` `1.04e+15` nulp here meant the values differ by
   **8 %** (`rel 4.3e-02 … 8.9e-02`, 216/216 elements), while the sphere's
   `array_equal → False` on two arrays that PRINT identically was **1 ULP**
   (`rel 2.06e-16`). One measurement re-sorted the whole investigation and took 5 minutes.
2. ⭐⭐ **The quadrature-family split IS the discriminator, and the passing sibling is the
   evidence.** `[M]` 5 of 5 failing rows used `level_symmetric`; 1 of 1 passing used
   `lebedev`, bit-identical at 0/126 elements — through the *same* sweep, dispatch and
   reflect helper. A shared-machinery bug cannot be family-selective, so that single green
   row refuted "the sweep regressed" instantly. Ask what the green rows have in common
   before bisecting anything.
3. ⭐ **A gross move in the SCALAR flux refutes an ordering hypothesis on sight.**
   `φ = Σ_n w_n ψ_n` is permutation-invariant to `N × ULP`; if φ moved 2.8–6.2 %, the rule's
   VALUES moved, not its order. (And the ordering hypothesis was right about scale and
   mechanism elsewhere: `[M]` the real 1-ULP mover permuted nothing — imposing a rule's
   declared symmetry makes derived ordinates bit-copies of the seed octant, changing the
   last bit of 16 of 24 nodes with the order untouched.)
4. ⭐⭐ **A rule-tier TOLERANCE pin can never warn about a consumer-tier BIT-IDENTITY pin —
   the gap is structural, not an oversight.** `[M]` `gauss_legendre` is gated against
   numpy's `leggauss` at `< 8 * ulp` (correct: neither construction is "the" answer). A
   3-ULP change is INSIDE that contract and OUTSIDE every downstream `array_equal` /
   `sha256` consumer, simultaneously. The only instrument that closes it is a **byte-level
   fingerprint beside the tolerance pin**, whose sole job is "the bytes moved — re-baseline
   the consumers".
5. ⭐⭐ **A declared re-baseline's blast radius is the set of FROZEN REFERENCES, not the set
   of `.npz` FILES.** The causing commit said "those baselines are re-captured in the
   following commit, not silenced" and that commit re-captured 22 snapshots — but missed a
   `sha256` hex string living in a `.py` module and a hand-written arithmetic expression
   living in a test body. Neither is a file a regeneration script can see. Enumerate frozen
   references by KIND (stored array / digest literal / in-test formula) before declaring a
   re-baseline done.
6. ⭐ **A gate whose reference is a different FP ASSOCIATION of production's own expression
   is a coin flip, not a contract.** `[M]` production computes `src + (A + B)` (the helper
   returns `A+B`); the test writes `src + A + B` = `(src + A) + B`. Over the full
   (5 cells × 8 ordinates × 2 groups) grid: **68 of 80 slots bit-identical, 12 differ, max
   ULP = 1**. The passing sibling passes by ordinate-index luck. Re-baselining the number
   leaves the fragility; the reference must be re-associated or the assertion demoted.
7. ⚠ **Two probe-harness self-inflicted failures, both in the flattering direction.**
   (a) `grep -cE "^FAILED"` reported `nfailed=0` at 9 of 9 commits including two already
   measured RED — pytest's ANSI codes precede `FAILED` (`vv-principles` #17, third class).
   Use `--color=no`. (b) `python -c` prepends **CWD** to `sys.path[0]`, AHEAD of
   `PYTHONPATH`, so a worktree probe silently imported the MAIN tree and printed HEAD's
   values for every commit. Run probe **script files** located outside the repo, and make
   every probe print `module.__file__`.


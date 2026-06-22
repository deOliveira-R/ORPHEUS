# Numerics Investigator — Lessons

## L1: Run the diagnostic cascade in order — no skipping

Six wrong hypotheses were wasted on cylindrical DD divergence by guessing
before isolating. Steps 3-5 (fixed-source, component isolation,
per-ordinate analysis) identify the root cause directly. The cascade
order is not optional.

## L2: Curvilinear redistribution is the prime suspect

In cylindrical/spherical DD, alpha recursion and geometry factors
(signs, delta-A/w scaling) are where bugs hide. Diverging-with-refinement
keff in curvilinear geometry means the balance equation is wrong. After
confirming spatial streaming works alone (step 4, alpha=0), go to
per-ordinate flat-flux consistency (step 5).

## L3: Rank-N closure improvements must show signed-error stability across ≥2 quadrature schemes

A single quadrature can fake apparent rank-N convergence. Any claim of
structural improvement at the optically-thick regime (σ_t·R ≥ 10)
requires at least two independent quadratures (e.g. RICH = (4, 8, 64)
vs RICH+panels = (5, 8, 64), or Gauss-Legendre vs tanh-sinh) showing
the **same signed error** — not just the same magnitude — before the
result is taken seriously. The pathology is "L17/L19 quadrature
crossing": the rank-N error and the F.4 baseline error cross zero at
different points along the quadrature-refinement axis, so a single
quadrature snapshot picks one sign for rank-N and the opposite for
F.4 and reports a false structural win.

Falsified Direction-C / Direction-Q / Direction-N attempts (Issues
#121 closed, #122 closed, #123 open at devcontainer budget) all
looked promising on a single RICH-grade quadrature and collapsed
under the two-quadrature stability check. Use
`assert_rank_n_structural_win` in
`tests/cp/test_peierls_rank_n_protocol.py` as the gate before
shipping any new closure.

## L4: Convergence-rate fingerprints discriminate failure modes

A solver stalled at 1–10 % vs reference is not enough information.
Run an n-doubling sweep and check three discriminating products
of (error × n^p):

| Empirical pattern              | Failure mode                            |
| ------------------------------ | --------------------------------------- |
| `err·sqrt(n) ≈ const`          | Schneider C^(0,1) endpoint singularity  |
|                                | in the *solution* (graded mesh repairs) |
| `err·n/log(n) ≈ const`         | Log-singular kernel diagonal truncation |
|                                | in the discrete operator (product-      |
|                                | Nyström repairs)                        |
| `err·n² ≈ const`               | Smooth-integrand quadrature undersam-   |
|                                | pling (just bump n)                     |
| `err·n^p, 2 < p < 4`           | Simpson on a piecewise-quadratic        |
|                                | trial; rate limited by solution         |
|                                | regularity                              |
| `err`-ratio not monotone in n  | Eigenvalue-iteration tolerance, not     |
|                                | the discretisation; tighten the         |
|                                | iteration test                          |

ERR-036 (Wave 2-A Path A.i, 2026-05-03) was misclassified by the
literature memo as a Schneider endpoint-singularity problem
(`err·sqrt(n)`); the empirical fingerprint was `err·n/log(n)`,
pointing to log-singular kernel diagonal truncation. The
methodology recommendation (Atkinson product-Nyström) was right
but the *justification* was wrong — graded mesh would have been
secondary, not primary. Always run the three-way fingerprint
attribution before committing to a fix; a wrong attribution
chooses the wrong fix even when both fixes happen to be
beneficial.

## L5: Read the paper's stated approximation level before assuming a code bug

When a small numerical disagreement (1-10%) with a published reference
appears, **read the paper's own caveats explicitly before assuming the
gap is a code bug**. ERR-038 (Wave 2 Front-3, 2026-05-03) was investigated
as a "singular limit at R→1 needing multi-day asymptotic analysis" when
in fact Atalay 1997 explicitly states (p.236, p.246) that the published
Tables 2-5 are **first-order Fredholm approximations** with degraded
precision at small slab thicknesses. The 5% R=0.99 gap is the paper's
own precision floor.

The discriminator that resolves "code bug vs paper floor":

| Observation                                  | Likely cause                       |
| -------------------------------------------- | ---------------------------------- |
| Uniform offset across all parameters         | Code bug (constant, sign, factor)  |
| Scales with a physical parameter (1/d, 1/τ)  | Paper floor (omitted higher-order) |
| Scales with a numerical parameter (n, dps)   | Quadrature or precision bound      |
| Insensitive to all numerical knobs but exact | Paper floor or fixed approximation |
| at one structural limit                      | level                              |

Always include a "moderate-parameter consistency check": at a parameter
value where the paper's stated approximation is tight, the solver MUST
agree to machine precision. This is the structurally-independent ground
that supports the verdict "paper floor, not code bug" — without it, the
case stays open (per V&V principles §1).

Also caught the X-function-precision false-positive: an upstream
Signature 7 (1.2% rel-diff between legacy mpmath maxdegree=14 and
tanh-substituted mpmath at ν=0.99) does NOT automatically propagate to
every downstream quantity at the same magnitude. Sensitivity analysis
(X→K_j→2d_crit: 1.2% → 0.3% → <0.1%) is mandatory before claiming the
upstream fix solves the downstream gap. Plain magnitude propagation is
NOT the default — many integrals are insensitive to upstream-spline
small shifts when the cancellation structure absorbs them.

## L6: A curvilinear matvec is verified only against a NON-FLAT per-ordinate hand reference

> **→ The standing rule is now in AGENT.md Step 5 (per-ordinate
> analysis).** Kept here for forensic value: the specific bug
> codenames and the twin-path mechanism.

Flat ψ makes every redistribution/routing bug in a curvilinear SN matvec
vanish: index permutations become the identity (all level-internal values
equal), and cell-average-vs-cell-centre mismatches → 0 in the continuous
limit. Krylov-on-apply then converges to a self-consistent fixed point of
the buggy operator — physically acceptable, still linear+conservative — so
nothing looks wrong. Two distinct cylinder matvec bugs (ERR-049 bool-mask
scatter ordering; the decoder "analytical extension" O(h) twin-path
divergence) BOTH hid for months behind flat-ψ-only coverage; see
[[cyl-matvec-twin-path-signatures]] (#197/#206). The gate: every curvilinear
matvec needs a per-ordinate hand-reference test on NON-FLAT ψ (the sphere
had `test_apply_face_fluxes_match_sweep_recurrence_spherical`; cylinder
lacked the analog until ERR-049). This is L2 sharpened into a test-design
rule, and is vv-principles H2/H3 (homogeneous + conservation are both
degenerate to redistribution/routing bugs) applied to the matvec twin of
the sweep. The architectural cure is Pattern 2: ONE `SNCellOperator`
consumed by both sweep and matvec so they cannot drift (#206).

## L7: A direction-dependent per-ordinate moment must be lifted to the global frame BEFORE the angular reduction

ERR-061 (#240 D5b-S3, 2026-06-17). An LD slope moment ψ̂_n produced in the
per-ordinate SWEEP frame (downstream-positive) and summed by the angular
reduction φ̂=Σ_n w_n ψ̂_n (which feeds the isotropic scattering slope source
Σ_s·φ̂) MUST be sign-corrected to the GLOBAL frame for backward ordinates first
— else forward and backward slopes CANCEL (φ̂ ~6× under-driven) and the scheme
loses the diffusion limit (LD nx=4: 38.9% off DD). The fix is a per-octant
moment-frame involution ∏_a (octant_sign_a)^{o_a} applied at the producer-
consumer seam (source/probe global→sweep IN, moment/residual sweep→global OUT;
the outgoing FACE stays sweep-frame, it propagates the wavefront).

Three durable methodology points:
1. **The angular reduction is the discriminator.** A per-ordinate convention bug
   is invisible until a quantity is summed across ordinates of OPPOSITE sweep
   direction. This is H2/H3 (flat flux nulls the slope; round-trip/conservation
   are telescoping-degenerate to the frame error) applied to the moment iterate.
2. **Matvec self-consistency (SI≡Krylov, round-trip≈0) is NEVER sufficient for a
   moment-iterate fold.** It proves the operator is internally consistent, not
   that its fixed point is physically correct (vv §5). Gate the converged VALUE
   against a structurally-independent reference, not the round-trip — the brief
   named this trap and it held.
3. **When all components are individually correct but the FP is wrong, build a
   STRUCTURALLY-INDEPENDENT from-scratch kernel.** It reproducing the wrong value
   bit-for-bit localizes the bug to the SHARED math (the 2×2 + the frame
   convention); the independent fix confirms the fix class without touching
   production. When the literature says the scheme IS consistent but your
   faithful implementation isn't, the bug is a CONVENTION between two correct
   pieces (the sweep↔global frame at the seam), not in either piece.
See [[issue-240-d5b-s3-diffusion-limit]].

## L8: The project's own theory page can be the contaminated reference — trace to literature + an independent kernel

ERR-063 (#257 S10a, Peierls MG fission χ sink/source swap, 2026-06-21). The
suspect line read χ at the SINK node `i` (`B[i·ng+ge,j·ng+gs] += K·chi[i,ge]·νΣf[j,gs]`)
where the SOURCE node `j` belongs (χ must share the source index with νΣf — the
fission emission density is a single LOCAL birth quantity at the fission point).
The trap: `docs/theory/peierls_nystrom.rst` Eq. `peierls-mg-operator` AND the
`solve_peierls_mg` docstring math BOTH documented the SAME wrong `χ_g(r_i)` — so
the theory page is NOT a structurally-independent ground (vv §6 reference
contamination). Confirming the code against its own doc would have confirmed the
bug. Two structurally-independent grounds settled it: (1) Hébert 2009 Eq. 3.57/3.58
(textbook integral-transport, χ shares the source argument; via literature-researcher);
(2) the sibling Variant-α `trajectory_resolvent` solver — a DIFFERENT kernel family —
already source-indexes χ (`chi_nodes=chi[region_at_node]`, local `chi_nodes·F_r`).

Durable methodology:
1. **A masked O(1) bug needs a spatially-VARYING field to surface.** χ_i≡χ_j on
   HEAD's library (every region same χ) hid the swap; the swap is visible ONLY when
   the field varies across regions. This is H2 (homogeneity nulls the term) applied
   to a material field, not a flux field. The discriminating probe MUTATES a
   non-fissile region's χ and asserts k_eff invariance — a region with νΣf=0 emits
   zero fission neutrons, so its χ is causally meaningless; only sink-indexing leaks
   it (measured 4.7% k_eff move under the bug, 0 under the fix).
2. **Byte-identity on the original data is the "no re-baseline" proof.** Stash-vs-fix
   on the uniform-χ HEAD library gave Δ=0.000e+00 — the fix changes nothing where χ
   is uniform, so the Hébert Class-B benchmark pins need no re-baseline. Always prove
   a "fixes-the-bug-without-moving-the-references" claim with a DIRECT old-vs-new
   value comparison on uniform-field data, not by inspection.
3. **Twin-path sweep (Cardinal Rule 2) extends to sibling MODULES, not just sites.**
   The fission-assembly pattern lived at 3 sites in `peierls_nystrom` (geometry.py +
   slab.py interior + slab.py white-BC) — all fixed. The sibling `trajectory_resolvent`
   module uses a DIFFERENT assembly (per-node local product, not a K·χ·νΣf matrix) and
   was already correct — which made it the independent witness. Grep ALL assembly
   modules, classify each as same-bug / different-structure-correct; the correct
   sibling is often the structural ground you need.
See [[issue_257_s10a_peierls_chi_source_swap]] (if logged) + ERR-063.

## L9: A reported "hang / non-convergence" may be a FIXTURE cost — bound the solver directly FIRST

A "SN heterogeneous-keff non-convergence hang" (#212) was NOT a solver bug: the
solver converged in ~0.3 s, but the test's `continuous_get(...)` fixture eagerly
walked the whole reference registry into O(minutes) of adaptive-mpmath Peierls
solves. A timeout is NOT a convergence failure until the solver is bounded
independently of its data fixtures.

How to apply: when a timeout appears in a test that consumes a registry-walk /
builder fixture (`continuous_get`, `*_cases()`, any auto-discovery), FIRST bound
the solver directly — tiny `max_outer`/`max_inner`, or bypass the fixture by
calling the producing module's builder. If the solver converges fast, the cost is
in the FIXTURE; do NOT bisect the solver (Step-3 fixed-source isolation already
exonerates it). Reproduce the slow path on a fast proxy
(`timeout 60 python -c "continuous_get('<name>')"`) and bisect on THAT, not the
test. Generalises: separate fixture-cost from solver-cost before attributing a
timeout to numerics. See [[sn-keff-hang-was-eager-registry]].

## L10: "Error grows with refinement" + a discarded library info-flag = an unconverged inner solve, not a discretization bug

A Krylov-inner `solve_sn` gave keff 25% wrong on a refined homogeneous sphere
(error GREW with mesh: 4e-10 @ n=5 → 5e-1 @ n=20), the classic
discretization-inconsistency fingerprint (cascade Step 1). But the cause was an
under-converged GMRES whose scipy `info > 0` flag was discarded (`solution, _info =`)
plus a hardcoded `restart=min(50, …)` that truncated the Krylov subspace below the
problem size (328). Tightening `inner_tol` did NOT move the wrong fixed point —
the discriminator that rules out tolerance and points at subspace truncation.

Two durable points:
1. **A discarded convergence flag turns "diverges with refinement" into a false
   discretization-bug signature.** Before concluding the operator is inconsistent,
   confirm the inner linear solve actually converged — re-run with the library's
   info flag checked, and sweep `restart` (info flips 200→0 and keff snaps to exact
   at restart ≈ problem size). A hardcoded subspace cap is the ERR-004 magic-constant
   anti-pattern (Pattern 7): derive `restart` from the operator's domain size.
2. **An L1 anchor can pass by NUMERICAL COINCIDENCE at one mesh size.** The anchor
   used n=10, where restart=50 GMRES happens to project onto the true eigenmode; the
   bug only declared itself at n=20. A single-mesh anchor is not coverage (vv H5) —
   always pair the analytical anchor with a mesh-refinement leg.
See [[krylov-restart-truncation-bug]].

## L11: For a ρ-honest stopping/diagnostic, measure the residual r=Aψ−q — NOT the iterate increment ‖Δψ‖

> **→ The standing rule is now in AGENT.md Step 3 (fixed-source
> diagnostic).** Kept here for forensic value: the FluxDisplacement /
> AngularResidual typed-diagnostic catalogue and the cross-family
> ρ/(1−ρ) amplification links.

The SI increment `Δψ = ψ⁽ⁱ⁾−ψ⁽ⁱ⁻¹⁾` UNDERSTATES the true error by `1/(1−ρ)`:
`‖Δψ‖ ≈ ρ‖Δψ_prev‖`, true error `= Δψ/(1−ρ)`. At c=0.99 (ρ≈0.99) a "converged at
tol" iterate is ~100·tol from the fixed point. The equation residual `r = Aψ−q`
(rate-density units) is ρ-HONEST — it measures distance from the solved equation,
independent of the iteration scheme, and does not shrink artificially as ρ→1.

How to apply: when a solver REPORTS converged but the answer is wrong by a
c-dependent factor that GROWS as c→1, suspect a ρ-blind stopping test before
suspecting the operator. The `FluxDisplacement.contraction_ratio` (>1 diverges /
≈1 stalled / <1 healthy) turns ‖Δψ‖ honest; the typed `AngularResidual.balance_map`
is the per-ordinate flat-flux residual (vv H3 — exposes per-ordinate balance failure
that telescoping conservation hides). This is the SAME ρ/(1−ρ) amplification as the
curvilinear-MG inner-tol drift and the L10 restart-truncation family. See
[[issue-208-flux-displacement-residual-typing-debug-value]].

## L12: An OFFLINE-isolated error is only THE floor after an end-to-end swap + a silent control — and AMPLIFY is the sharpest disproof

A term whose error is real and proven in isolation (the τ-clamp thread error
1e-3→1e-15; an LD boundary slope moment) can be a SECOND-order contributor masked
by a larger end-to-end error, or be genuinely consumed yet sub-floor so it never
moves the converged value. An offline magnitude alone never settles "X is THE
floor" or "X improves the answer."

The three durable probes:
1. **End-to-end swap + a term-SILENT control.** Run the full solver with the term
   toggled AND a control configuration where the term is identically zero (an
   isotropic MMS for the redistribution term; a flat/zeroed input for a slope) —
   the control byte-identity pins the asymmetry and proves the toggle is non-vacuous.
2. **AMPLIFY the suspect term.** Growing the term (bc_scale, a larger slope) is the
   sharpest disproof of "it improves the answer": if amplifying makes the converged
   value strictly WORSE, the improves-on-flat hypothesis is refuted cleanly. A
   correctly-consumed sub-floor term can legitimately make the value slightly worse.
3. **First-cell / boundary-row ORDER is the mechanism oracle, not the global L2.**
   If the boundary cell is already O(h²) with the term off, the term cannot repair an
   O(h) deficiency that isn't there. The global volume-weighted L2 DILUTES a localized
   defect (√V ~ h^1.5 at the pole) — measure the order at the cell, not in the norm.

This is L7's "matvec self-consistency ≠ correct fixed point" generalized to any
offline-isolated quantity, and pairs with vv Mode-10 (activated-but-unconstrained):
the resolution is structural teeth + a no-op control, never a tightened value band.
See [[curvilinear-tau-clamp-vs-pole-floor]], [[issue-257-s9-ld-boundary-slope-optical-verdict]].

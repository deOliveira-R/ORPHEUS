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

## L13: "Diverges/crashes for LD (or any spectator trailing axis)" = a greedy `(Ellipsis, *idx)` fancy-index that mis-targets axes

When a moment-tensor / einsum path is bit-identical for the scalar (single-moment)
case but DIVERGES or raises `IndexError` once a trailing spatial-moment (LD's 2^d φ̂
axis) — or any spectator broadcast axis — is present, the prime suspect is a per-cell
fancy index built as `cells = (Ellipsis, *idx)`. `Ellipsis` is GREEDY from the front:
it absorbs the leading (m, g, …) axes, so the spatial cell indices `*idx` land on the
LAST k axes — which under a trailing φ̂ axis are `(…, last-spatial, φ̂)` instead of the
two spatial axes. Symptom is geometry-dependent: rectangular grid (nx≠ny) → `IndexError:
index N out of bounds for axis … with size M` (the spatial index over-runs the smaller
wrong axis); square grid + asymmetric material map → SILENT wrong value (cells scattered
to wrong positions, e.g. 43% rel error). The fix is to pin the leading axes EXPLICITLY:
`cells = (slice(None), slice(None), *idx)` so `*idx` targets the spatial axes and the
trailing spectator broadcasts. Bit-identical when no trailing axis (Ellipsis ≡
(slice,slice) then), so non-LD tests are BLIND — the gate must use an LD flux
(`spatial_moments=2`). This is vv failure-mode #2 (variable/axis swap) gated by a
spectator-axis presence; the catalog signature is "frame/moment path OK for scalar,
breaks for LD". Isolation = monkeypatch the verb back to `(Ellipsis,*idx)` (NEVER git
checkout) and confirm it reddens on a RECTANGULAR LD flux. ORPHEUS instance: #276 A2
commit `0b3275d` fixed all four `MaterialXSField` moment-scatter verbs
(`apply_legendre_scattering_moments`{,`_transpose`}, `apply_n2n_moments`{,`_transpose`});
diagnostic `derivations/diagnostics/diag_276_full_scatter_kernel_ld_trailing_axis.py`.

## L14: The curvilinear `(L+C).solve` seed-lag is QUADRATURE-dependent, not geometry-uniform — slab direct; cyl is DEAD (level-symmetric) or LAGGED-but-FOLDABLE (product); sphere WAS lagged (fixed by route (a))

> **CORRECTED 2026-07-05 (#280 Phase 2.5b):** the original "direct inverse for
> slab+CYLINDER" headline was OVER-GENERALIZED from a level-symmetric probe. The
> cylinder is quadrature-dependent — see [[curvilinear-inverse-seed-taxonomy]] Verdict
> + the #280 §. **product-cyl cold err = 0.575** (NOT 0); the "α-dome telescopes the
> seed away" mechanism was a MIS-ATTRIBUTION (real cause: LS has `c_in[m0]=0` dead seed;
> product has `c_in[m0]≠0` live via the #229 clamp). The product lag is RETIRABLE by a
> **pure-diagonal fold** (κ=dA_w[m0]·c_in[m0] into the m0 cell diagonal — POC single-pass
> = M⁻¹ at 5e-16); fixed point is BIT-IDENTICAL (keff/MMS/matvec gates don't move). The
> whole MMS ladder is still blind to it (Mode-7, ≤linear-in-μ = seed's exact regime).

`(L+C).solve` (the WDD sweep) is seed-independent + machine-precision ONLY where the
angular-redistribution seed cancels. Measured (removal-form `InvertibleOperator`, ≥2G,
random non-flat ψ): SLAB seedΔ=0.0 / residual 8e-16; CYLINDER **quadrature-dependent**
(level-symmetric seedΔ=0.0 / residual 7e-16 via a DEAD seed `c_in[m0]=0`; **product cold
err 0.575** — a LIVE foldable self-coupling, NOT telescoping); SPHERE seedΔ(X1,X2)=4.6e-2 /
`‖Aψ−b‖∞/‖b‖∞ = 5e5` (WAS lagged, FIXED by route (a) `a29ab2d`, L15). So a curvilinear
inverse is NOT uniformly a `SweepOperator`. The SOLE lagged element is the M-M half-angle starting
seed ψ_{1/2} per level, read from `initial_guess.bulk.values` (`_initial_guess_values` →
`closure.psi_half_seed(psi_level,ctx)` in `loss_representation/__init__.py:3162/3197`;
`None`→zero). Default seed is `AngularEdgeExtrapolation` (NOT `CarlsonInwardSweep` —
ERR-058 superseded it), exact on flat AND linear-in-μ. **It creates a LOCAL CYCLE**: the
seed reads the two most-inward ordinate cell-AVERAGES, and ordinate-0's redistribution
consumes the seed — so the sweep breaks the cycle by lagging (everything else — pole
r=0 continuity capture, the M-M recurrence — is feed-forward within the sweep). Seed-dep
is therefore a FORMULATION choice, not intrinsic: the μ=±1 equation is CLOSED ((1−μ²)=0,
`carlson_inward_sweep_from_source` already implements it) → a direct sphere sweep needs
only to resolve that cycle (explicit ψ(·,μ=−1) state / per-level block-solve / source-
driven seed) = exactly issue #200's block-inverse face preconditioner. **V&V punchline
(Mode-7):** every curvilinear MMS ansatz is ≤ linear-in-μ (isotropic `sin(πr/R)`, or
`(A(r)+B(r)μ)/W`) — precisely the seed's EXACT regime — so SI converges O(h²) on the
whole ladder and the seed-lag instability is INVISIBLE. A genuinely higher-order-in-μ
field (a plain uniform-source sphere) makes the seed-iteration DIVERGE (→NaN under SI at
every c∈[0,0.99]); production dodges it by shipping GMRES with the IDENTITY precond
(`_within_group_krylov`, `solver.py:332`, #200) and by keff being shape-independent. This
is L6 ("curvilinear needs a NON-FLAT per-ordinate reference") realized end-to-end.
Diagnostics: `/Users/rodrigo/.claude/jobs/84fd66f8/tmp/diag_curvilinear_seed_sensitivity.py`
(+ `diag_sphere_fixedpoint_consistency.py`). See [[curvilinear-inverse-seed-taxonomy]].

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

## L16: To compare two ANGULAR QUADRATURES' accuracy in a curvilinear SN scheme where they differ in the POLE/SEED treatment — build a standalone scheme-faithful driver gated to production, reference = fine-N + a cross-quadrature CONTAMINATION GUARD (MMS is blind), validate the new pole handling with the per-ordinate flat-flux residual

A DESIGN study (Gauss-Legendre vs Gauss-Lobatto for spherical SN, 2026-07-06): would GLob
(nodes AT μ=±1) make the M-M starting-direction seed ψ½ an ordinary weighted ordinate, erasing
the seed-block type machinery, and at what accuracy cost? Method that worked:
1. **Standalone scheme-faithful driver, NOT a production point-swap.** When the quadrature sets
   the angular DISCRETIZATION (M-M τ/α + pole route), not just moment integration, a point-swap
   hits SINGULAR coefficients — GLob's μ=−1 node lands ON the lower angular edge → raw
   **τ_0=0** → the recurrence ψ_{3/2}=(…)/τ_0 divides by zero. So μ=−1 MUST be a straight
   characteristic (Carlson DD march, (1−μ²)=0), NOT run through the recurrence (caveat 3);
   production's UNCLAMPED sphere τ + its `starting_direction_levels` predicate (τ_raw,0∈(0,1))
   both break on GLob ⇒ wiring it in is real surgery. Reimplement the EXACT production M-M
   weighted-diamond (verify every coeff file:line) parametrized by (μ,w,pole_mode).
2. **GATE the driver bit-faithful to production GL** on NON-FLAT vacuum bare spheres (flat/
   homogeneous-reflective is H2-degenerate). Hit rel 1e-11 keff / 1e-10 flux vs `solve_sn`,
   THEN swap only quadrature+pole — so a driver bug can't masquerade as a GLob effect.
3. **Reference = fine-N GL + the cross-quadrature CONTAMINATION GUARD** (compute fine-N GLob
   too, confirm GL_∞==GLob_∞ to ~1e-6): makes the reference quadrature-family-UNBIASED (vv-6).
   Report BOTH error-vs-ref AND the reference-free matched-N |GL−GLob| diff.
4. **MMS is BLIND** (L14/vv Mode-7): every curvilinear MMS ansatz is ≤linear-in-μ = the seed's
   exact regime, certifying neither the seed nor the quadrature accuracy. Anchor correctness
   with a closed form (k_inf=1.875 hand-derived) + φ=Q/Σt streaming equilibrium instead.
5. **Per-ordinate flat-flux residual** (vv-H3/L6) validates the NEW pole handling per ordinate
   (angle-integrated φ is degenerate): all pole modes gave max |ψ_n−C| ~1e-15.
Verdict: GLob tracks GL at a bounded ~1.2× error penalty at resolved N (N≥8, N>L), regime/c/
anisotropy-insensitive; affordable for the architectural win. Full recipe + numbers in
[[glob-vs-gl-spherical-quadrature-study]]; artefacts `scratch/experimental/glob_sphere_study/`
+ `derivations/diagnostics/diag_glob_0{1..5}_*.py` (33 tests green).

## L17: A STATE DOF's Hilbert metric is NOT its angular-integration weight — and when apply_transpose is the EXACT Euclidean transpose (T=Aᵀ), a block metric is GAUGE-FREE (any SPD), the determining equation only forbids DEGENERACY

Deriving the SN curvilinear starting-direction (ψ½) block metric `G_sd` for the augmented
composite `A` on `bulk⊕trace⊕seed` (#282/#280 2.5d, Mode-12 closure, 2026-07-06). The
"ghost metric" `G_sd≡0` was justified as the angular through-flux coefficient
`(1−µ²)|_{µ=±1}=0`. That reasoning is WRONG: it confuses the angular-INTEGRATION weight of
the µ=±1 ray (correctly zero — which is why the seed does NOT appear in the scalar-flux
reduction Σ_n w_n ψ_n) with the STATE metric of a discrete DOF. ψ½ is NOT a quadrature node;
it is a first-class radial state field with a nonzero self-block `A_ss` (‖·‖=4.0) and a
nonzero seed→bulk coupling `A_bs` (‖·‖=6.0, the M-M recurrence). Its Hilbert metric is fixed
by its ROLE in the operator algebra, not an angular weight. The radial-field VOLUME makes it
nonzero; the pole angular weight is a red herring.

**The determining equation `Aᵀ G = G A†`.** The linchpin: is the IMPLEMENTED `A†` the
metric-daggered `G⁻¹AᵀG`, or an independent transpose kernel? Read `_AdjointOperator.apply`
(operator.py:1146): `A.H = G⁺·apply_transpose·G` — apply_transpose is INDEPENDENT (#280
`_seed_rows_transpose` + reverse walk). MEASURED (dense unit-vector probe of `A.apply` /
`A.apply_transpose` / `A.H.apply` over `to_flat`, sphere GL-S4 2G nx=4): **T = apply_transpose
== Aᵀ EXACTLY** (‖T−Aᵀ‖=3.6e-16, incl. `T_sb==A_bsᵀ`). ⟹ `A.H=G⁺AᵀG` is the honest
metric-adjoint for ANY invertible G ⟹ **reciprocity `⟨Aψ,φ⟩_G=⟨ψ,A.Hφ⟩_G` is GAUGE-FREE**:
holds for EVERY SPD `G_sd` (V_cell 1e-16 / identity 6e-17 / V·w 1e-15 random-seed defect).
The determining equation pins `G_sd` ONLY up to SPD (in the production diagonal-metric
architecture: any strictly-positive diagonal). Gauge is PHYSICAL: `A.H` is block-upper-
triangular with seed at the TOP, so its bulk⊕trace rows are BITWISE gauge-invariant (Δ=0.0
exact across identity/V_cell/10·V_cell) — only the internal φ†_seed moves; no observable
reads it. Recommended fixing = **V_cell** (radial volume, matches bulk `G_bulk=V·w_n`; the
angular w is the sole gauge d.o.f., no canonical value for a single ray).

**Three durable methodology points:**
1. **`G_block=0` is the ONE forbidden value, and it's WORSE than "blind."** Measured: with
   `G_sd=0` a random (nonzero-seed) reciprocity probe BREAKS at 1.3e-2 on the production path
   — the shipped `A.H` is a WRONG adjoint the instant the seed carries data (the `A_bs`
   coupling is unmatched, `A.H` severs the seed: `H[seed,:]=H[:,seed]=0`). It looks correct
   ONLY because the gate feeds a present-but-ZERO seed. A zero-weight block in a Hilbert
   metric is not a conservative default; it silently corrupts the adjoint off the zero-probe
   regime (vv Mode-12 sharpened).
2. **Closing a Mode-12 invariant-functional blindness needs TWO changes, not one.** (a) install
   the non-degenerate metric AND (b) feed NONZERO block data in the gate — a zero-block probe
   can't activate the block's rows (A_ss·0=0, g_s·0=0) even with a perfect metric. The
   existing `test_mode12_..._blind_to_a_seed_row_flip` POSITIVELY PINS the blindness with a
   zero seed; after the fix it must INVERT to assert the flip REDS.
3. **The dense-probe recipe for any augmented-composite block metric.** Assemble A (forward),
   T (apply_transpose), G_b/G_t (production `apply_metric`) as dense matrices via unit-vector
   probing of `FullField.to_flat`; check `T==Aᵀ` (numpy transpose = structurally-independent
   ground, NOT the operator's own machinery). If T==Aᵀ ⟹ gauge-free, need SPD; if T carries a
   weight ⟹ pinned by `g_s[i]=g_b[coupled j]`. Faithfulness `G⁺TG==A.H` (2.8e-14) ties the
   dense reconstruction to production. Forward stays bit-identical under the install (metric
   read only by A.H + inner_product, #208 trace-metric precedent). Diagnostics
   `derivations/diagnostics/diag_gsd_0{1,2,3}_*.py` (17 green). See
   [[starting_direction_metric_gauge_derivation]].

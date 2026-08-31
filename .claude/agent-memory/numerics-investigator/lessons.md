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

**L10b (2026-08-08) — the flag need not be DISCARDED; RECORDED-BUT-UNREAD gives the
identical false fingerprint, and the tell is TOLERANCE-INSENSITIVITY.** A gate red at
3.29e-10 was bit-identical across `inner_tol` ∈ {1e-9,1e-11,1e-13,1e-15} — read as "a
real discretization floor, since the exit ignores the tolerance". It was the opposite:
the running residual at the cap was 1.185e-09, i.e. ABOVE every tested tol including
1e-9, so all four runs hit the SAME `max_inner` and returned the same iterate. **A
tolerance sweep can only discriminate if at least one tol is LOOSER than the residual
the capped run reaches — otherwise the sweep is a vacuous knob** (the L18
verify-the-knob-is-reachable warning, in convergence clothing). FIRST read
`history.converged` / `n_inner` and check them against `max_*`; a plateau at exactly
`max_iter−1` is the whole diagnosis. Sharpenings: (a) the flag was honest and typed —
the hole is the CONTRACT, a best-effort and a certified answer returning as the same
type from the same call, so a certificate that no-ops on the unclaimed exit ("best-
effort returns stay legal") leaves the class wide open; (b) the *shape* discriminator
that kills "per-ordinate discretization bias" dead is a **max_iter sweep of the error
MAP** — a bias is fixed, a mode decays geometrically at constant ρ with its shape
preserved (measured ρ=0.9852±0.0008 over four intervals, shape cosine ≥0.993, the
suspicious "8 of 24 ordinates" concentration constant at 100 % throughout); (c) audit
the blast radius with an in-process pytest plugin wrapping the solver entries and
printing every `converged=False` at `sessionfinish` — it found the red gate's SIBLING
riding the same truncated exit while GREEN on a looser rtol (a latent false green), in
one run.

**L10c — an all-reflective (zero-leakage) pure absorber is SI-HARD, and the cost
EXPLODES with dimension.** `[M]` d=1 32 → d=2 258 → **d=3 1631** sweeps at
`inner_tol=1e-13`; one vacuum face collapses d=3 to 208. Mechanism: with `Σ_s=0` and no
leakage the only coupling is the reflective boundary and the only damping is absorption,
and the slow mode is the **DD face sawtooth** — faces alternating `ψ*(1±δ)` (measured
ratios 1.074414/0.925586, summing to exactly 2) so the cell AVERAGE is exactly `ψ*`.
Having zero cell average, the mode is invisible to the collision term `Σ_t V ψ_c`; it is
damped only through the inter-axis balance mismatch. Consequence, measured and useful as
a budget law: **`Σ_t · n_inner` is invariant** (0.4→3093, 0.8→1631, 1.6→850, 3.2→437;
products flat to 13 %). So a FIXED `max_inner` default cannot serve d=1 and d=3 at once
— derive it (Pattern 7). Corollary worth checking before trusting a d=3 rate default:
**boundary-G-S is 2.5× FASTER than Jacobi at d=2 (258 vs 648) and 1.95× SLOWER at d=3
(1631 vs 838)** — a rate claim measured at one dimension can INVERT at the next, so a
`is_cartesian and not is_1d` schedule gate that was tuned at d=2 is unverified at d=3.
Also: `ψ` flat does NOT imply the TRACE is flat — the converged trace kept an 11.26 %
sawtooth while the exact-uniform state also had residual 1.06e-15, i.e. the all-
reflective 3-D DD operator has a near-null space in the trace block (L6/H2 sharpened:
"homogeneous nulls redistribution" has a twin — "flat cell averages null the face mode").

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

## L18: To adjudicate a LABELING/ORDERING degeneracy in a discrete-ordinate scheme, the instrument is the operator's own SYMMETRY GROUP — MMS is exactly blind, and the answer is usually "no ordering is right, the closure is broken"

#326 (2026-08-01, cylindrical per-level ordinate order, `rules_product.py:139`
`argsort(mu_x)` where η is 2-to-1 over φ). Five durable points:

1. **An MMS whose ansatz AND source depend only on the INVARIANTS of the degenerate class
   is EXACTLY blind** (Mode 7 by design + Mode 12: the relabeling is in the measured
   functional's stabiliser). Both ORPHEUS cylindrical ansatzes are functions of (η, ξ²) —
   exactly the two things the azimuthal mirror pair shares — so the pair carries identical
   ψ_ref, Q, w. Measured: two tie-breaks agree to 3e-12 / 9e-15 on the production ladders,
   at every mesh, at every quadrature order. **Declare the ansatz's invariants BEFORE
   trusting it as an adjudicator**; a companion ansatz that leaves the symmetric sector
   (here ξ-ODD `(A+B ξ)/W`) SEES it (20.6%) yet still cannot ADJUDICATE — both orderings
   converge to the SAME angular floor from opposite sides at spatial order ~0.
2. **A within-tie permutation cannot move a cumulative-sum coefficient.** α is a partial
   sum of `w_m η_m` and `w` is constant within a product level ⟹ α and τ are BIT-identical
   across tie-breaks; only WHICH ORDINATE sits at each position moves. So the whole effect
   is a labeling, and the symmetry-defect MAGNITUDE is ordering-invariant.
3. **The leak path is any coupling that maps ACROSS the degenerate class.** The tie-break
   need not COMMUTE with it. Here the r=0 pole seed `pole_outflow[reflection_index("x")[n]]`
   sends ordinate n to the −η class where the tie-break made an independent choice —
   measured non-commuting for 24/64 ordinates — which is why even ξ-EVEN data moves 2.6e-2.
   Grep every cross-class index map (`reflection_index`, mirror/partner tables) before
   concluding a relabeling is inert.
4. **The adjudicating criterion is a SYMMETRY the continuous AND semi-discrete problems
   both have** — no reference solver, no MMS, structurally independent by construction.
   1-D cylindrical: ψ is EVEN in ξ; the product rule is closed under ξ→−ξ with equal
   weights. Production breaks it at 1.19e-1 (30% local; LS4 3.08e-1), flat in n_mu, falling
   in n_phi ⟹ it IS the #229 azimuthal floor, seen WITHOUT a reference. Verdict: no ordering
   is correct — the M-M 1-D η-march on a level with duplicate η is, plus the [½,1] clamp
   that turns the resulting ZERO-WIDTH angular cell into an arbitrary {1, ½} τ split.
5. **The constructive exit is to fold to the FUNDAMENTAL DOMAIN.** A degeneracy in a sort
   key usually means the discretization carries a redundant symmetric copy. On ω∈[0,π]
   (the independent half) η is strictly MONOTONE: no ties, ordering UNIQUE, and every
   competing criterion coincides — α is simultaneously a non-negative dome AND exactly
   `2 w_gl κ ξ(ω_{m+½})` (the user's closed form, κ=Δω/(2sin(Δω/2))=1+Δω²/24, verified
   3e-16 via the Dirichlet kernel), and the ξ-mirror holds by construction.

**Method warning (cost me a probe):** verify the varying knob is REACHABLE. My first
tie-break probe was vacuous — production's trig-evaluated nodes split the "ties" by 1 ULP,
so lexsort/stable/quicksort all agree and the tie-break is not a free variable until the
nodes are algebraically exact (#325). A control leg asserting "the two settings really
DO differ" caught it. This is vv Mode-8's SIGNATURE-tautological class in numeric clothing:
the knob existed, the value it should have varied was decided upstream by rounding noise.
See [[cylindrical-level-ordering-symmetry-adjudication]].

**[M] Promotion sharpening (2026-08-01), two measured corrections worth reusing:**
(a) **`xfail` swallows FIXTURE SETUP ERRORS too**, so "move the solve into a fixture so an
incidental failure surfaces as an ERROR" is FALSE — measured: a raising stub for the solve
left all three strict-xfail rows reporting `3 xfailed` in 0.32 s. The working structure is
the xfail row asserting ONLY its documented inequality PLUS a **reddenable un-xfailed
SIBLING** consuming the same fixtures (a well-formedness/activation band); break the solve
and the sibling ERRORs loudly. Pair it with the positive control — simulate the fix and
confirm `XPASS(strict)` — which is what proves the row measures what the repair changes.
(b) **A diagnostic that RE-IMPLEMENTS the production kernel must be rewired to CALL it at
promotion** (vv Mode 11). My α diagnostic ran a local copy of the recursion; the promoted
gate drives `cylindrical_streaming`, and the proof is a mutation applied through the test's
OWN import binding (flip the production α sign → 20 L1 rows red, the 15 foundation
derivation rows stay green — which also validates the L1-vs-foundation marker split).

---

## L19: An iterative-solver RATE question is a spectrum question — eigen-solve `M⁻¹N`, never re-time the solver; and a "G-S beats Jacobi" claim needs its comparison THEOREM checked, not its fixtures re-run

#341 (2026-08-09, boundary Gauss-Seidel vs Jacobi inverting between d=2 and d=3 SN).
Five transferable points, in order of leverage:

1. **Build the iteration matrix, don't re-time.** Any `x ← A_inv.apply(Σ gᵢ.apply(x))`
   driver IS a linear operator: wrap it over the composite's `to_flat`/`from_flat` as a
   `scipy.sparse.linalg.LinearOperator` and run `eigs(which="LM")`. Cost ≈ a few hundred
   sweeps for the whole spectrum vs thousands to converge one case, and it is **immune to
   the stopping test** — so it cannot be contaminated by a ρ-blind stop (L11) or a
   `max_inner` truncation (L10). Positive control: it must reproduce the ρ *fitted* from a
   real residual history (measured 0.98552 vs 0.985348 and 0.97541 vs 0.975014, 4 decimals).
   Ship two cheap linearity checks with it (`G(2x)=2G(x)`, `G(0)=0`) — they also prove the
   `initial_guess=` seed is inert on the geometry under test.
2. ⭐ **Before hunting for why a splitting comparison inverted, ask whether the theorem that
   forbids it still applies.** Varga's comparison (`ρ_GS ≤ ρ_J` whenever `A = M−N` are
   *regular* splittings with `N_GS ≤ N_J`) makes the inversion IMPOSSIBLE for a
   non-negative iteration matrix — so an observed inversion is first of all evidence that
   the operator is **not** non-negative, and the productive question is *which term is
   negative and why*. Here: the multi-D diamond closure's face transmission is
   `Σ = (2/D)·1wᵀ − I` (`w_a = 2|μ_a|A_a`, `D = Σ_tV + Σw_b`) — one damped eigenvalue
   `1 − 2Σ_tV/D` plus **`d−1` eigenvalues exactly `−1`**: an *undamped* zero-cell-average
   face sawtooth (`ψ_c = 0`, so `Σ_tVψ_c` cannot see it) whose dimension grows with `ndim`.
   Step differencing would give the same `d−1` modes eigenvalue **0**. Any all-reflective
   DD rate pathology should be read through that spectrum first.
3. **A per-axis SIGN is usually a gauge — check before blaming signs.** On the octant
   hypercube `Q_d` (specular reflection flips one cosine ⟹ the coupling graph IS `Q_d`),
   flipping the sign of one axis's gain is a diagonal similarity, so both `ρ_J` and `ρ_GS`
   are invariant (measured identical over all 8 sign patterns). Sign *indefiniteness* is
   necessary to void the theorem; the sign *pattern* explains nothing. The same 2^d-scalar
   model — full hypercube, exact fold, exact ordering — **never inverts**, which localises
   the mechanism to the intra-octant `d×d` block the model discards. A model that fails to
   reproduce the effect is worth as much as one that does: it deletes a whole hypothesis class.
4. **Enumerate the discrete design space instead of sampling it.** The G-S fold here is
   fully described by `L_a` = the constant-sign suffix run of the octant sweep order
   (derived, then measured: `Σ L_a` implicit rows out of `d·2^d`). All `8!` orders collapse
   to **25** patterns; measuring all 25 gave an exact separating law
   (`LOSES ⟺ max_a L_a > Σ_{b≠a} L_b`, 25/25) and revealed a **2.5× spread** in the rate
   with the shipped order 24th of 25. When a knob's reachable values are finite, sweep them
   all — a sampled sweep would have produced a fitted story instead of a law.
5. **Ask whether the two arms are racing the SAME mode before calling a change a "flip".**
   Extract the dominant eigenvector and report where its mass sits (per face, per ordinate
   class, spatial sign pattern). Measured: at d=2 both splittings race the same y-face
   sawtooth; at d=3 G-S races the z faces and Jacobi the y faces — so it was two different
   comparisons, and the deep fold had merely unmasked (and degraded) a different survivor.

**Verdict discipline that generalises:** a production default must not branch on a variable
you have only *correlated*. `ndim` was falsified on both sides by direct measurement (3 d=2
fixtures where G-S loses, 3 d=3 where it wins) — optical thickness, mesh, aspect ratio,
quadrature order and `c` all move the sign at fixed `ndim`, and near-critical `c ≥ 0.99`
removes the effect entirely. Also worth the habit: when a docstring names a *theorem*
("the regular splitting", "`ρ_GS ≈ ρ_J²`"), that word is a checkable claim, and here both
were measurably false. Full record: `scratch/issue_341_gs_jacobi_mechanism.md`.

---

## L20: A RESIDUAL cannot gate an EIGENVALUE contract — measure the TRANSFER GAIN before proposing any residual threshold, and never approximate a signed adjoint projection

#340 N5 (2026-08-10): does an outer certificate `defect = ‖Aψ − Fφ(ψ)/k‖/‖Fφ/k‖` (the
production within-group `_certify_within_group_exit` lifted one level) discriminate a
*corrupting* truncated inner from a *benign* one? **REFUTED**, on 38 solves / 8 geometries /
3 mixtures. Five transferable points.

1. ⭐⭐ **The one number that decides any "can statistic X gate contract Y?" question is the
   TRANSFER GAIN `|Δy| / X`, measured across configurations — compute it FIRST, before any
   threshold hunt.** A threshold on `X` bounds `|Δy|` only through that gain, so if the gain
   is unbounded no constant exists and the tuning exercise is void. `[M]` `|Δk|/defect` spans
   `1.152e-05 … 1.340` = **1.16e+05×** ⟹ populations overlap **634×**; a zero-false-alarm
   threshold misses **15 of 16** corrupting cases, and the whole trade-off curve is
   unusable (100 % sensitivity costs a **59 %** false-alarm rate). This is vv **Mode 12 read
   in the MIRROR**: the usual failure is a gate BLIND to the error class; here the gate is
   SIGHTED on a class the contract is blind to. `‖r‖` is up to **99.995 %** reflective-trace
   rows (`bnd_frac`), and a reflective inflow-trace defect in a zero-leakage system carries
   no net current ⟹ `k = production/absorption` cannot see it by conservation. Two
   functionals, different invariance groups, neither containing the other ⟹ simultaneously
   OVER- and UNDER-sensitive (a truncated row read a *lower* defect than the fully-converged
   one). The cure is to project onto the functional the contract reads: the angle+volume
   integrated per-group rate defect `R_g = Σ_n w_n Σ_i V_i r` cut the overlap 634× → **4.64×**
   (14/16 caught at 2/22 FA) — good enough to REPORT as a number, still not a GATE.
2. ⭐ **A SIGNED projection against an approximate weight is WORSE than no weight** — it
   manufactures accidental near-cancellations, i.e. false NEGATIVES. First-order perturbation
   theory (`δk/k = ⟨ψ†,r⟩/⟨ψ†,Fφ/k⟩`) is the correct statistic, but with a spatially-FLAT 0-D
   adjoint it degraded the overlap 4.64× → **128.95×** and the gain spread to **2.27e+05×**
   (one CORRUPT row collapsed 46×, gain 20.6). The weight was *verified* (`|k_pencil − k_inf|
   = 0.00e+00`, the hand-built `A⁻¹F` pencil as its own positive control) — it was correct
   for the WRONG PROBLEM. So: either pay for the real adjoint or use the unsigned norm; never
   the cheap signed shortcut.
3. **Answer the NULL case before the discrimination question, and separate the two levels'
   slack with a two-legged tolerance sweep.** The certificate's "pass" value was **3.47e-07 =
   3.47 × keff_tol** — not machine-zero. Sweeping `inner_tol` 1e-08→1e-14 at fixed outer moved
   it **0.2 %**; sweeping the outer at fixed inner moved it **6 decades** (→3.79e-15). So it
   was the OUTER's own increment-stop slack (L11), not a floor, and not the inner's. Without
   both legs the 3.47e-07 reads as a structural floor and the whole study is mis-anchored.
4. **When lifting a production certificate one level, the CONSTANT does not come with it —
   and copying `record.binding_criterion.tolerance` silently picks the LOOSER criterion.** The
   eigenvalue outer's binding criterion was `dphi` (tol `flux_tol`) in *every* solve measured,
   never `dk`: `SAFETY × keff_tol = 1e-6` catches 8/16, `SAFETY × flux_tol = 1e-5` catches
   **2/16**. A residual bar scaled by an INCREMENT tolerance is a category error twice over.
5. ⭐ **Gate every verification fixture on the mixture's own consistency
   `σ_t == σ_c + σ_f + Σ_to SigS[0][g,:]` — an inconsistent mixture makes two legitimate
   references DISAGREE with no bug in either.** The brief's benign pole ("keff correct to
   2.5e-11 vs `k_inf`") did not reproduce: `|k − k_inf| = 6.9e-02` (30 %). Cause: the fixture
   wrote `sig_s` in `[to,from]` (its own `# 0 -> 1` comment says so, and `σ_c = σ_t −
   s.sum(axis=0)` is the correct removal for that) while `make_mixture` reads
   `SigS[g_from,g_to]` ⟹ `σ_t` off by **±0.12**. In a zero-leakage medium the transport
   balance (removal `σ_t`) gives `0.23076923076923` and production/absorption gives
   `0.30000000000000`; the SN reports the second, `solve_homogeneous_infinite` the first. Two
   further poisons rode along: `φ₁ ≡ 0` (effectively **1-group**, vv anti-#3) and `c = 0.9`
   giving `σ_c = −0.14`. One character (`sig_s=s.T`) repairs everything — and the repaired
   fixture DOES exhibit the intended benign pole (`|Δk| = 1.10e-11`, 4/4 inners truncated at
   200/200, `ρ ≈ 0.985`). **Never trust a brief's reference value on a hand-built mixture
   until the consistency identity is printed.** Full record:
   `scratch/n5_outer_certificate_measurement.md`.

## L21: An ANGULAR-consistency claim is separated by `h → 0`, not by the physical parameter it is named after — and a "sweep the regime" design can self-destruct

#319 / #235 flux-dip discriminator, 2026-08-12, 251 solves all `converged` (record:
`scratch/q68_flux_dip_discriminator.md`). Six transferable points.

1. ⭐⭐ **When a scheme claims consistency in the limit of ANOTHER discretisation
   (angular consistency "in the diffusion limit"), the axis that separates it from a
   rival is the OTHER mesh going to zero — because the claim is exact only there.**
   `[M]` sweeping optical thickness at fixed cells-per-mfp separates the shipped
   Morel–Montry τ from plain diamond **not at all**: the defect is constant to FOUR
   figures over `Σ_t·R = 5…50` for BOTH (fitted decay rate `0.000`), on sphere and
   cylinder. Refining `h` at fixed physics separates them **without bound**: the good
   scheme's defect → 0 at exactly first order, the rival's SATURATES, ratio
   `3.2× → 204×` over 2 → 64 cells/mfp. Ask "which limit is the claim exact in?" before
   choosing the sweep axis.
2. ⭐⭐ **A regime sweep at fixed `c` self-destructs: `Σ_a·R = Σ_t·R(1−c)` grows with it,
   the interior becomes a source plateau, the current at the origin dies, and every
   scheme agrees for a reason unrelated to the question.** `[M]` at `c=0.99,
   Σ_t·R=100` three τ schemes agreed to **3 significant figures**; at 300 the metric
   was `2e-10` for all. Reading that as "equal ⟹ hypothesis refuted" is the trap. Use
   the ε-scaling `Σ_t=1/ε, Σ_a=ε, Q=ε` (holds `Σ_a·R = O(1)`) and **carry a
   fixture-liveness column** (the smooth profile's own variation over the first few
   mfp) that declares when the fixture stopped posing the question.
3. ⭐ **Build a λ-CONTINUUM through the two candidates, not an A/B.** `τ(λ)=λτ_A+(1−λ)τ_B`
   turns "A beats B" into "is A the MINIMISER?", and `λ_opt(h)` is then a falsifiable
   curve. `[M]` sphere: `λ_opt → 0.993 / 1.001` (two instruments) as `h→0` ⟹ the shipped
   τ is the optimum of the family, and any apparent optimum below 1 on a coarse mesh is
   spatial. Cylinder: `λ_opt → 0.73`, and the two instruments DISAGREE ⟹ no optimum
   claimed — the disagreement is itself the finding (two consistency conditions that
   coincide on the sphere decouple there = a missing angular DOF).
4. ⭐ **A theory scalar can be τ-loaded or τ-blind depending on which EDGES you feed it,
   and the blind version is the natural one to write.** M&M's `β` (Eq. 6a) built from the
   STANDARD weight-partition edges is τ-blind *by construction* (that substitution IS
   their β=0 proof); built from the edges the CLOSURE implies
   (`μ̃_{m+½}=(μ_m−(1−τ_m)μ̃_{m−½})/τ_m`, `μ̃_½=−1`) it is solve-free, exactly zero for the
   shipped sphere τ at every order, and the measured anomaly is LINEAR in it. It also
   catches `τ→1−τ` (the Mode-12 reflection the membership/fold-box/reversal gates are
   exactly blind to). ⛔ But it is **identically zero for BOTH schemes on a folded
   cylinder at every `n_φ`** — a spherical invariant does not transfer to a geometry
   whose angular derivative is in a different variable.
5. ⭐ **A literature diagnostic transfers between geometries only in its LEVEL-LOCAL
   form.** M&M's effective starting cosine `(ψ_s−φ)/(3J)` reads a `+2.76` artefact on a
   cylinder because the on-axis flux is azimuth-independent but genuinely POLAR-angle
   dependent; rebuilt from the level's own zeroth/first azimuthal moments it reduces to
   the published formula on the sphere bit-for-bit and gives sane cylinder values.
   ⚠ It is an S2/S4-class instrument — it presumes ψ affine in the level's angle, so at
   S8/S16 the genuine curvature dominates and it reports a fixed bias, not a defect.
6. **The benefit of a low-order-consistency fix DECAYS WITH ANGULAR ORDER and can
   invert.** `[M]` sphere `14× (S2) → 1.4× (S4) → 0.9× (S8, A worse) → 0.9× (S16)`,
   tracking `β(B)` falling 5 orders; cylinder `5.4× (n_φ=8) → 2.3× (n_φ=16)`. So an
   accuracy comparison run only at high N will report the principled scheme as a
   regression — correctly, and for a reason that is not a bug.

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

---

## L23: A discrete operator's SINGULARITY is a two-question object — measure `dim ker`
## against the count of the "benign" rows, and REFUSE the either/or the brief hands you

#344 (2026-08-14, `A = L + C − S − B` on an all-reflective Cartesian box). The brief
offered two exclusive readings — "benign tangential-slot bookkeeping" vs "real trace
underdetermination". `[M]` **both are true and additive**, and on the fixture the whole
campaign had measured (`level_symmetric`) the benign one contributes **exactly zero**.
Six transferable points.

1. ⭐⭐ **When a claim is "the null space IS class X", the decisive number is
   `dim ker` MINUS `|X|`, and you must measure `|X|` — not assume it exists.** `[M]`
   `dim ker A = 12` (d=2) / `138` (d=3) against **0** tangential `(face, ordinate)`
   pairs, because a **level-symmetric rule places every cosine on a shell `|μ| ≥ μ₁ > 0`
   and CANNOT produce `Ω·n = 0`**. The brief (and the issue) asserted the opposite family
   property. One line — `min |omega_dot_n|` — refuted the framing before any solve.
   `dim ker = |X| + R` held **exactly** on 9 (geometry × quadrature) rows, which is what
   turned an either/or into a decomposition: `product(4,4)` is `R = 0`, `level_symmetric`
   is `|X| = 0`, `lebedev(11)` and `product(8,8)` carry both.
2. ⭐⭐ **A dense SVD through the production builders is CHEAP and settles rank questions
   that ARPACK only bounds.** The prior record said "3 unit modes at d=2, ≥6 at d=3" from
   `eigs(G, k=12)` — a **lower bound**, and the true answers are 12 and 138. `[M]` unit-vector
   probing of the composite `to_flat()`: 1248 dof → 2.0 s build; 7392 dof → 30 s build +
   124 s SVD (437 MB). Report the **singular-value GAP** (`σ[-13]/σ[-12] = 9.5e+12`) so the
   rank threshold is visibly not arbitrary — a threshold anywhere in `[1e-13, 1e-2]` gave
   the same rank.
3. ⭐⭐ **Two blindness mechanisms that look alike are told apart by the METRIC, and that
   difference decides the remedy.** A tangential slot carries `G = |Ω·n|·w_n` **exactly
   `0.000000e+00`** ⟹ *no* G-weighted functional can ever see it ⟹ typing it away
   (Pattern 4) is the only fix. The real-underdetermination rows carry `G ≥ 1.83e-01` and
   the null shift measures **`3.97e-02` relative** in the G-norm ⟹ a gate CAN exist; typing
   cannot remove a rank deficiency. Always ask "is this class in `ker G`, or merely in
   `ker` of the functionals I happen to gate with?" — `vv-principles` #18 covers only the
   first, and the second is the commoner case.
4. ⭐ **A residual-based stop and a conservation projection are blind to `ker A` BY
   CONSTRUCTION — so the only informative half of that measurement is the POSITIVE
   CONTROL.** `A(ψ + αv) − q ≡ Aψ − q` is a theorem, not a finding. `[M]` both functionals
   sat at `~1e-16` under an 11.26 % trace shift while a NON-null perturbation of the **same
   flat 2-norm** moved them to `3.40e-01` and `1.16e-02` (15 and 14 orders). Budget the
   probe for the control; the "unmoved" column is free.
5. ⭐ **A converged solver's deviation from the analytic answer is in `ker A` EXACTLY —
   test it with `‖Aδ‖/(‖A‖‖δ‖)`, no null basis needed.** Any fixed point of `ψ ← M⁻¹(q+Nψ)`
   satisfies `Aψ = q`, so the difference of two fixed points is a null vector. `[M]`
   `3.97e-14` on a solve reporting `converged=True` at 1614 sweeps — which also refutes
   "it is just an undecayed `ρ = 0.985` mode". And **identify the recorded scalar before
   trusting it**: the memo's `1.1258e-01` was `max|ψ/want − 1|` over ALL ordinates on the
   face; the printed ordinate-0 row is `7.44e-02`. Both reproduce; only one is the quoted
   number.
6. ⭐⭐ **Fit the counting law, WRITE THE PREDICTION DOWN, then test it off-sample — and
   swap the SCHEME to get the mechanism.** `[M]` d=2: `ng·N/4`, mesh- and
   scattering-independent; d=3: `ng·(N/8)·(2Σnᵢ − 1)`, **3 of 3 on held-out points
   including a change of quadrature order**. Preconditions measured, not assumed:
   `d ≥ 2` (d=1 is `0` for both BCs and both families) and **≥ 2 reflective axis pairs**.
   The mechanism was closed by one substitution: **`LinearDiscontinuous` on the identical
   box is NON-singular** ⟹ it is the DD closure's `ψ_out = 2ψ̄ − ψ_in` involution, whose
   zero-cell-average eigenspace has eigenvalue `−1` and is undamped by `Σ_t V ψ_c`
   (L19-2's algebra, now confirmed end-to-end). Blast radius worth naming: `A` excludes
   `F`, and the eigenvalue entry DEFAULTS to all-reflective, so **every `d ≥ 2` Cartesian
   DD k-eigenvalue solve runs a singular within-group operator** (`cond = ∞`; a direct LU
   is rank-deficient; `A[trace,trace]` alone has `dim ker = 168` of `672`, far worse than
   `A`'s `12`). Full record: `scratch/issue_344_null_space_structure.md`; gate:
   `derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py` (10 green, 58 s).

### L23 addendum (2026-08-14) — the DISPOSITION half: CORRECTNESS or DETERMINISM?

Three more transferable points, from settling whether the #344 singularity is a bug.

7. ⭐⭐ **To decide "is a solver-selection defect a CORRECTNESS bug?", change the
   SPLITTING — not the mesh.** A splitting cannot change the equation, so anything
   that moves under `inner_schedule` is the solver's, not the operator's. `[M]` same
   operator, same source, same ZERO cold start: boundary-G-S returns a trace
   `8.97e-02` (d=2) / `1.126e-01` (d=3) from the closed form, **Jacobi returns
   `6.4e-13` / `6.8e-13`** — 5 of 5 fixtures. So the frozen component is
   `P₁ψ_exact` with `P₁` the *splitting's* OBLIQUE spectral projector, NOT a
   property of `ker A` (which is splitting-invariant). ⟹ **the standard guarantee
   "a splitting changes the rate, never the fixed point" is VOID whenever the fixed
   point is a MANIFOLD** — different splittings select different members, and any
   schedule-invariance / DSA-FP gate on the TRACE will legitimately red with no bug
   present (`vv-principles` Mode 9, sharpened).
8. ⭐⭐ **A refinement ladder can be PARITY-SPLIT, and the `vv #13` break-the-
   congruence-class rule is what finds it.** `[M]` on cells `(n,n)`: the deviation is
   **identically zero at even `n`** (`1e-12`, and `‖Ad‖/‖d‖ = O(1)` ⟹ ordinary
   iteration residual, not a null component) and `6.2e-02` at odd `n`. A 4/8/16/32
   ladder reports "nothing to see". Two method points: **(a)** carry `‖Ad‖/‖d‖`
   beside the error so you can tell a frozen null component (`~1e-11`) from leftover
   residual (`~1`) — the error column alone cannot; **(b)** run the ladder INSIDE
   one parity class, then the law was exact: `err·n = 0.311671` to **8 s.f.** over
   `n = 5…31`. It CONVERGES, at O(h) ⟹ not "the wrong limit". `[M]` only `n_x`
   parity matters (11/11; `(3,4)` deviates, `(4,3)` does not) — the x-major octant
   order is the suspect.
9. ⭐ **Before recommending a gauge, measure whether the TRUE answer IS the canonical
   representative.** `[M]` `‖P_G ψ_exact‖_G/‖ψ_exact‖_G = 1.1e-15` ⟹ the exact
   solution is the minimum-`‖·‖_G` member of the manifold, so projecting the returned
   iterate off `ker A` is an **EXACT fix** (`8.97e-02 → 5.8e-13`), not a convention.
   And when enumerating what a null direction is invisible to, **enumerate the moment
   LADDER, do not reason about it**: I predicted `J⁺ ≠ 0` cancelling in the net, then
   predicted the spatial sum annihilated it — **both wrong**. `[M]` every linear trace
   functional whose angular weight is a function of `|Ω·n|` ALONE is annihilated
   per-face-CELL at `~1.6e-15` (`φ±`, `J±`, `|Ω·n|^p` for p = 0..3); what SEES it is
   the raw per-ordinate value (75 %), the QUADRATIC G-norm (43 %), and an
   angularly-resolving detector (5.8e-03) — whose adjoint problem is then
   **INCONSISTENT** (`‖P_null Σ_d‖/‖Σ_d‖ = 5.0e-02`).

### L23 addendum II (2026-08-14) — the COHERENCE half: is the splitting a splitting?

Settling "is boundary-G-S a splitting of A, or is it reflecting an inconsistent
trace (ERR-056's failure)?". Three more transferable points.

10. ⭐⭐ **To tell a GAUGE FREEDOM from an INCOHERENT solver, REMOVE THE KERNEL and
    re-run — and remove it ≥3 structurally-different ways.** A pure-trace `ker A`
    forces "boundary moves, bulk does not", and so does an incoherent schedule seen
    from a distance; the asymmetry is NOT the discriminator. `[M]` on four
    independent `dim ker A = 0` configs (vacuum pair / mixed `xmin`-refl+`xmax`-vac
    / **LD on the ALL-reflective box** / d=3 one reflective pair) × 2 source types,
    both schedules agree: trace `≤ 1.7e-12`, bulk `≤ 2.2e-13`. ⟹ COHERENT. ⚠ The
    obvious control can be a NON-control: an "even-`n_x` box, where your parity
    finding says the kernel is absent" has `dim ker A = 12` — what is absent is the
    kernel's EXCITATION by that source, not the kernel. Assert `dim ker == 0` inside
    the control, do not infer it from a deviation being zero.
11. ⭐⭐ **`‖M M⁻¹ − I‖` over the FULL space is the wrong instrument for an
    iteration's inverse — probe the RHS SUBSPACE the driver actually supplies.**
    `[M]` boundary-G-S reads `‖M M⁻¹−I‖ = 3.3e-01` / `‖M⁻¹M−I‖ = 1.9e+00` (Jacobi
    `3e-16`) — which reads as incoherence and contradicts the measured
    `‖Aψ*−q‖/‖q‖ = 8e-14`. Resolve by BLOCK-DECOMPOSING the defect: it is
    **exactly** the (inflow-row, outflow-column) block, and the driver's
    `r = q + Sψ + B_upper ψ` has **exactly `0.000e+00`** outflow-trace content, so
    on that subspace the inverse is exact (`1.6e-15`) and the fixed-point identity
    holds (`4.9e-13`). A reified forward-substitution "inverse" is a **SUBSPACE
    inverse** by construction — fine for SI, a live hazard for a Krylov
    PRECONDITIONER, which feeds it arbitrary vectors. When two of your own
    measurements contradict, one is the wrong instrument: find which directions it
    probes that production never supplies.
12. ⚠ **My own probe was wrong first, twice, both flattering toward the alarming
    verdict**: it densified `M⁻¹` with `initial_guess=x` (a seed VARYING with the
    probe vector) and "checked linearity" on the DENSE MATRIX — linear by
    construction, so the control could never fail (`vv-principles` #17). Corrected
    controls: `initial_guess ∈ {0, random, b}` all bit-identical, operator linearity
    `0.00e+00`. **And the positive control is what makes the whole verdict mean
    anything**: the ERR-056 mutation (reflect after the FIRST outflowing octant
    group, not the LAST) drives the same comparison to trace `1.0000e+00` **and bulk
    `0.39…0.80`** on kernel-free configs — twelve orders of dynamic range, and it
    shows that an incoherent schedule hits the BULK too, which is the direct answer
    to "should a schedule change touch bulk and boundary equally?".

---

## L24: A discrete operator's KERNEL is usually a CLOSED-FORM problem — the tell is that the measured counting law is independent of a physical parameter, and the route is "substitute the degenerate branch of the scheme's own closure back in and see what CANCELS"

#344 (2026-08-14), `ker A` for `A = L+C−S−B` on an all-reflective Cartesian DD box.
Prior sessions had two FITTED counting laws (`ng·N/4` at d=2, `ng·(N/8)(2Σn−1)` at
d=3, off-sample 3/3) and no basis; the brief asked for a basis or a no-closed-form
verdict. A basis exists, in `0.05 s` where a dense SVD is `23 s` at half the size.
Record: `scratch/issue_344_kernel_basis.md`. Six transferable points.

1. ⭐⭐ **A counting law that does NOT depend on a parameter the operator plainly
   contains is telling you the governing equation is COMBINATORIAL — go derive it,
   do not fit it.** `dim ker` was measured mesh-independent at d=2, `c`-independent,
   cross-section-independent, exactly `∝ ng`. All four survive the substitution:
   set the degenerate branch (`ψ_c = 0`, read off the measured `1.1e-28` bulk share
   of the null projector), which turns DD's `ψ_out = 2ψ_c − ψ_in` into the
   involution `ψ_out = −ψ_in`, so every face field is the sawtooth
   `ψ_a(k,i_⊥) = (−1)^k φ_a(i_⊥)`; substitute into the balance and **every**
   cross-section, mesh width, weight and area cancels, leaving
   `Σ_a s_a Y_a(s_{≠a}; i_{≠a}) = 0` — *"a sum of functions, each blind to one
   coordinate and one sign, vanishes identically"*. Both laws then drop out as
   theorems. **The parameter-independence WAS the derivation hint.**
2. ⭐⭐ **Sign CHARACTERS diagonalise a specular-BC constraint system.** A specular
   BC says a quantity is blind to one SIGN; the balance says it is blind to one
   COORDINATE. Expanding in `χ_T(s) = ∏_{b∈T} s_b` splits the coupled system into
   one INDEPENDENT additive-separable (ANOVA) equation per character subset `U`,
   with `dim = κ(U)·∏_{c∉U} n_c`, `κ(U) = Σ_{a∈U}∏_{b≠a}n_b − ∏_U n_b + ∏_U(n_b−1)`
   (`κ(pair) = 1`, `κ(triple) = Σn − 1`). Basis = **pair generators**: pick two axes,
   a character, and an index tuple on the rest. ⟹ read a **SUM over axes** in a
   counting law as "the modes live on PLANES (one free coordinate)", never as a fit.
   The orbit count `N/2^d` is the number of ordinate orbits under the reflection
   group — NOT ordinates per octant (that reading is off by `2^{d−1}`, since at d=2
   the `±μ_z` ordinates are in different orbits).
3. ⭐ **Where an SVD is unaffordable, the span check is a PRODUCTION-GENERATED
   kernel vector — with a round-off NEGATIVE control.** `[M]`
   `‖(I−P)(ψ_GS − ψ_exact)‖/‖·‖ = 9.8e-13` at d=3 `ndof=7392` and **`1.000000`**
   in-span at three odd-`n_x` d=2 meshes. The control: on the Jacobi arm the
   deviation is `9e-13` (pure round-off) and reads **`1.000e+00` OUT of span** —
   proving the projector captures kernel content and is not a universal absorber.
   Without that leg, "everything I test is in the span" is unfalsifiable.
4. ⚠ **Two mechanisms can share a PARITY fingerprint, and only a kernel-CONTENT
   measurement separates them.** I hypothesised the known even-`n_x` split was
   DETECTOR blindness (the mode's transverse profile is `(−1)^{i_⊥}`, so a uniform
   detector is exactly blind at even cell counts — `[M]` **`0.000000e+00` on every
   functional including the odd controls** at `(4,4)`). True, and NOT the cause:
   `[M]` at even `n_x` the deviation is `2e-13…5e-13` with only `15–31 %` in
   `ker A` ⟹ the mode is **ABSENT**, not hidden. Measure `‖P d‖/‖d‖`, not `‖d‖`.
5. ⭐ **A blindness LIST measured on one quadrature is a sample, not a population.**
   The prior "every `|Ω·n|^p`, p=0..3, is blind" was measured on `level_symmetric`,
   where the tangential component `T = 0`. `[M]` on `lebedev(11)` the **`p = 0`**
   moment (a plain face-averaged scalar flux) reads the T modes at **`2.99e-02`**
   while `p ≥ 1` stays `<1e-17`. Honest condition: **mirror-EVEN in angle AND ≥ 1
   power of `|Ω·n|`** (a CURRENT-type functional). Corollary theorem worth reusing:
   every mode carries a non-trivial sign character on every axis it touches, so
   **any mirror-even angular weight annihilates the kernel exactly** — which is also
   why `ψ_exact ⊥_G ker A`. And matching matters: a `sign(μ_xμ_y)` weight is BLIND
   where `sign(μ_x)` sees 4.4e-2 — "angularly resolved" ≠ "sighted".
6. ⭐ **A DENSE basis is not the shippable form of a STRUCTURED nullspace.** Disjoint
   supports per (orbit, group) ⟹ `BᵀGB` is block-diagonal ⟹ `17.6 GiB → 154 MiB`
   at `(12,12,12)` S8 ng=4, apply `12 ms`. ⛔ And `ker G ∩ ker A ≠ 0` is a real
   hazard: the tangential slots have `G` **bit-zero**, so `BᵀGB` is SINGULAR and a
   `sqrt(G)`-QR gives `0/0` — there is NO minimum-norm gauge for them. Project on
   the G-positive component only. (Cost verdict: setup `0.04 s`, apply `0.094 ms`
   against a `7.9 s` solve; the basis never reads a cross-section, so it is built
   once per PHASE SPACE and cached — fissile vs absorber gave a bit-identical
   `2.799e-16` residual.)

---

## L25: An ARITY / "does this contract widen?" question is a THEOREM question — ask what COMMUTES; and a positivity question is a DIFFERENT question the asymptotic expansion is structurally blind to

Curvilinear LD × Morel–Montry τ (2026-08-25/26). Record:
`scratch/tau_under_ld_dip_analysis.md`; probes `scratch/probe_tau_ld_0{1..6}_*.py`
(**106 SymPy checks, 0 failures**). Seven transferable points.

1. ⭐⭐ **"Does accessor `f(x)` have to become `f(x, y)`?" is answered by asking what
   the defining conditions COMMUTE with — not by an expansion.** A **scalar convex
   combination commutes with every linear map**, so if a quantity is defined by (i)
   membership of a convex set and (ii) exactness of a scalar blend, then *both*
   conditions are the SAME scalar statement in every component of every linear
   representation, and the widened form is an **overdetermined system whose every
   row returns the same value**. `[M]` solved independently at `n_mom = 4`: all
   rows give `τ = (μ_n−e_−)/(e_+−e_−)`. Hypotheses to CHECK, all three cheap: the
   scalar is independent of the widened index; the projection is linear; the set is
   convex. Cost: minutes, against a multi-day expansion.
2. ⭐⭐ **A SIGNATURE is not an invariance proof, and a spatial-argument-free
   producer is the classic tell.** "β cannot accept a mesh, therefore β is
   mesh-independent" is `vv-principles` Mode 8's SIGNATURE-tautological class: the
   signature is a consequence of the DERIVATION'S SCOPE (M-M 1984 and BMC 2010 both
   hold space *continuous* and say so), not evidence about the joint problem. The
   conclusion was true; the argument carried zero information. **Always ask what
   analysis produced the object, not what its parameters are.**
3. ⭐⭐ **The asymptotic expansion and the POSITIVE CONE answer different questions,
   and the expansion is structurally blind to the cone.** A sign-alternating
   cell-to-cell mode is EXCLUDED BY THE ANSATZ, so no order of the expansion sees
   it — which is why Palmer–Adams carry "limits to a *stable* diffusion equation"
   as a SEPARATE acceptance criterion. `[M]` transmission sign ladder, all derived
   in one probe: **DD flips at `τ_opt=2` (Padé(1,1)), bare LD at `τ_opt=3`
   (Padé(1,2)), lumped LD NEVER (Padé(0,2), discriminant `1−2<0`)** — and
   `a_LD·τ_opt → −2`, i.e. sign-flipped and decaying only as `1/τ_opt`. ⟹ when a
   scheme is "verified in the diffusion limit", ask which of the two it was.
4. ⭐ **Where a redistribution/coupling operator is a TENSOR PRODUCT `R ⊗ A(θ)`
   with disjoint index sets, EVERY functional of it factors** and the θ-scalar is
   untouched by any change to `R`. `[M]` verified at `n_mom = 1,2,3,4`; `[M]` the
   scalar's free symbols contained no spatial symbol at all — **grep the free
   symbols, it is a one-line proof.** Corollary that pays: the leading-order
   discrete diffusion equation then carries no θ, so a SPATIAL defect cannot be
   repaired by θ and vice versa — **two orthogonal failure modes**, and conflating
   them is the whole confusion.
5. ⭐⭐ **When one knob is exact-by-construction, price the OTHER knob's error in
   the first knob's units — the ratio is usually the finding.** `[M]` with τ at
   Morel–Montry (β ≡ 0), a starting-cosine error of only **1.6 % at S4, 0.05 % at
   S32** reproduces the ENTIRE contamination τ exists to remove, because
   `β_eff(μ_s) = S + (μ_s+1)·S_e` is EXACTLY affine and `|S_e/S_diamond|` GROWS
   with N (24× → 333×). ⟹ the celebrated knob was never the risk; its neighbour
   was. Look for the affine law before assuming two knobs are independent.
6. ⚠ **A published RULE OF THUMB carries the order it was derived at.** M-M's
   summary "the dip is eliminated with any spatial scheme as long as the starting
   flux is not seriously UNDERestimated" is `dβ/dμ_s = S_e > 0`, and `[M]` **`S_e`
   flips sign between N=2 (their own test case, `+9.1e-1`) and N≥4 (`−1.1e-1 …
   −1.7e-5`)** — the safe DIRECTION inverts. Honest companion: `|S_e|` falls 5
   orders over the same range, so the direction inverts while the stakes collapse.
   Same family as L21-6 (a low-order-consistency benefit decaying with angular
   order). **Evaluate the quantity, never the heuristic.**
7. ⚠ **"Lumping" names ≥3 different operations and they disagree.** `[M]`
   Legendre-diagonal lumping destroys `R₀₁` and breaks the per-moment-row
   flat-flux identity; **nodal ROW-SUM lumping (which is what Palmer–Adams's FL
   actually is) PRESERVES `R₀₁` exactly** — only `R₁₁` moves, `ΔA/3 → ΔA` — so the
   identity survives untouched. And row-sum lumping the GRADIENT gives transmission
   **identically 0** (a degenerate scheme), which is why MWS lump everything except
   the gradient. ⟹ a ⛔ banner saying "X may not be lumped" must name the BASIS and
   the MATRIX; mine (inherited) condemned the wrong operation. The real freedom:
   the infinite-medium identity pins only `rowsum(L)`, leaving one parameter per
   row, and a MONOTONE (`A⁻¹ ≥ 0`) member exists inside it at the cost of dropping
   to first order — **the accuracy/positivity trade is a choice of that parameter,
   not a property of "lumping".**

## L26: An OWNERSHIP question ("who owns this quantity?") is answered by TWO measurable discriminators — does it move the METRIC, and is `∂A/∂k` DIAGONAL — never by what the quantity is called

Written 2026-08-26, the adjoint/Gram ownership audit of the 1-D curvilinear SN
streaming path (memo `scratch/adjoint_gram_ownership_audit.md`, probes
`scratch/probe_adjoint_gram_0{1..6}_*.py`).

1. ⛔ **"Does X influence the adjoint?" is the WRONG question when the adjoint is
   a reverse-mode VJP** — `[M]` `apply_transpose ≡ apply.T` to `3.05e-16`, so the
   answer is "everything" and it carries no information. The informative question
   is the **symmetry character of the INCREMENT `∂A/∂k`**, measured by perturbing
   one entry and taking the dense difference. ⟹ **build the difference operator,
   not the operator.**
2. ⭐⭐ **Calibrate the symmetry ratio before reading it.**
   `‖A+Aᵀ‖_F/‖A‖_F = √(2 + 2·tr(A²)/‖A‖²_F)`: **0 = skew, 2 = symmetric,
   √2 = `tr(A²)=0` (no entry pairs with its transpose partner, i.e.
   triangular-like)**. `[M]` controls: skew `0.000000`, symmetric `2.000000`,
   strictly-lower-triangular `1.414214`. Without the calibration a √2 reading
   looks like "very asymmetric" when it actually means *triangular*, which is a
   structural statement about a MARCH, not a magnitude.
3. ⭐⭐ **Two derived constants from the SAME two numbers on ADJACENT lines can sit
   on opposite sides of the self-adjointness split.** `[M]` `c_out = α_out/τ`
   contributes an **exactly diagonal** block (`offdiag = 0.000e+00`, sphere AND
   cylinder — self-adjoint in EVERY metric); `c_in = (1−τ)/τ·α_out + α_in`
   contributes an increment with `tr = 0` and **0 of its 8 nonzeros pairing with
   their transpose**. ⟹ a design that bags them as "the closure constants" fuses
   a reaction-like scalar with a transport-like coupling. **Diagonal ⟹ order-free
   ⟹ may live anywhere; non-diagonal ⟹ welded to the traversal that reads it.**
4. ⭐⭐ **A `⟨Aψ,φ⟩_G = ⟨ψ,A†φ⟩_G` reciprocity gate CANNOT adjudicate the choice of
   `G`.** With `A† ≡ G⁻¹AᵀG` it is an identity for EVERY invertible `G`: `[M]`
   `1.4e-16` under Euclidean `G′=1`, random `G′~U(1,2)` and adversarial
   `G′=(V·w)³`; the mismatch control (A.H on `G`, pairing on `G′`) reads
   `8.22e+00`. So the gate + its wrong-metric control prove **consistency and
   loadedness, never choice** — Mode 12 with the whole invertible-diagonal group
   as the stabiliser. ⟹ **the adjoint does not pin the metric; the physical
   FUNCTIONALS do** (scalar flux `Σw_nψ_n`, reaction rate `∫Σφ dV`), which is
   exactly why the metric belongs to the SPACE and the operator is its consumer.
5. ⭐ **One symbol, four roles, four owners — and the separating experiment is a
   GLOBAL RESCALE through the production constructors, not instance surgery.**
   `[M]` rebuilding at `w → 3.7w` (cylinder/slab, the charts that admit it) scales
   `G` **exactly** by `3.7` and leaves `L` **bit-identical** (`2.06e-16` / `0.0`),
   because `α ∝ w` ⟹ the redistribution sees only the ratio `α/w`. The sphere
   **REFUSES** the same rescale (`τ_raw ∉ [0,1]`) because its angular-cell
   partition is the *cumulative weight* and needs `Σw = 2`. And a SINGLE-ordinate
   `w` perturbation is refused on any curvilinear mesh by the pole mirror's
   weight-preservation contract. ⟹ metric / scale-free ratio / absolute mesh width
   / admission precondition. **Ask which one a consumer means before moving it.**
6. ⭐⭐ **A "Gram" in an OPERATOR is usually the mass matrix under a DIFFERENT
   measure, and naming the measure settles the ownership.** `[M]` (SymPy, exact,
   sphere+cylinder) the curvilinear redistribution Gram is
   `R_kj = ∫ b_k b_j (∇·ê_r) dV` with `∇·ê_r = (d−1)/r`, so `R₀₀ = ΔA` is the
   **divergence theorem**, not a per-chart normalization (the "factor-of-2-absorbed
   α normalization" is an artefact of writing the measure as `r^{d−2}`). It is SPD
   ⟹ a genuine inner product, but `W = M⁻¹R` is **not** `λI` ⟹ not `M` rescaled;
   `W` is the Galerkin compression of the Radon–Nikodym derivative. ⟹ **measure
   from the CHART, basis from the SCHEME ⟹ the home is the (chart × scheme) pair**,
   and it is a *coefficient*, never a metric — `[M]` perturbing it leaves `G`
   bit-identical. ⚠ And when the two axes differ (`n_mom ≠ n_thread`, the
   ONETRAN/Hill column) the object is **rectangular** ⟹ a PAIRING, not a Gram: the
   word over-claims off the diagonal of its own design space.
7. ⛔ **Do NOT justify a transport metric as "the one that makes streaming
   skew-adjoint".** `[M]` `‖ĜL+(ĜL)ᵀ‖/‖ĜL‖ → 1.4172 ≈ √2`, flat over `nx = 4…64`,
   **same number on the slab**: a face-ELIMINATED marching operator is triangular,
   not skew. Skew-adjointness is a property of the (cell ⊕ face) saddle system;
   DD's `ψ_out = 2ψ̄ − ψ_in` substitutes the interior faces out and destroys it.
8. ⚠ **Two structural ZEROES that a naive inventory would call "not used".** `[M]`
   `μ_start` is inert on every shipped 1-D curvilinear fixture (its only consumer
   runs on NON-carrying levels, and every level carries: `{0}` sphere,
   `{0,1,2,3}` `folded_product(4,8)`) — control `4.04e-01` when a level is forced
   non-carrying. And `c_in[m=0]` is inert because it multiplies the ψ½ SEED, which
   the ray-decoupled `(A,A)` block feeds with **zeros** — control: the neighbouring
   index `m=1` moves at `3.08e-04`. ⟹ **each zero needs its OWN control, and the
   cheapest one is the adjacent array entry.** (The instrument itself bit first:
   `sn.__dict__["areas"] = …` is a silent no-op because `areas` is a `property` —
   write the backing field `sn._areas`.)
9. ⭐ **The cylinder is blind to the Gram question — a THIRD member of the family.**
   `[M]` `R/ΔA` is bit-exactly `diag(1, 1/3)` on the cylinder (= the shipped
   moment-axis metric, so it reads as "not its own object") and has off-diagonal
   `h/(6r_c)` on the sphere. Sharper still: **the cylinder's MASS Gram `M/V` IS the
   sphere's REDISTRIBUTION Gram `R/ΔA`, exactly** — the two objects an ownership
   argument is trying to separate are the same matrix one geometry over. Joins the
   known cylinder blindnesses (the slope row reading `0 == 0`; σ_y-folded `β`).
   ⟹ **any curvilinear Gram/measure claim must be witnessed on the SPHERE.**

---

## L27: A "noise mode" reading of a small singular value is a HYPOTHESIS — kill it with the closed form and the RAW TABLE; and a numerical THRESHOLD has two edges, one of which its justifying instrument is structurally blind to

`_DENSE_METRIC_RCOND` re-derivation (2026-08-31, ERR-080/#429). A pinned
`pinv` cutoff was justified by *"one `~1e-16` **noise** mode … `1e-12` sits ~4
orders above the noise floor"*. Every quoted number reproduced; the reading was
inverted. Six transferable points, in order of leverage.

1. ⭐⭐ **Discriminate "round-off residue of an EXACT dependency" from "small
   real mode" by SOLVING for the dependency, then applying the candidate null
   vector to the RAW TABLE — not to the Gram.** Forming `G = AᵀWA` can itself
   manufacture rank loss, so `‖Gv‖ ≈ 0` is weak; `‖Av‖ ≈ 0` says the *columns*
   are dependent as functions on the node set and settles it. `[M]` here the
   closed form fell out in two lines (`Y_2^{+2} = (√3/2)(1−μ²)`,
   `1−μ² = (2/3)(Y_00−Y_20)` ⟹ `v = (1,0,−1,0,−√3)/√5`), matched the SVD null
   vector to `2.2e-16`, and gave `‖Av‖∞ = 1.4e-16`. ⟹ **a "noise floor" is a
   claim about a number that should have a NAME; if it has one, it is not
   noise.**
2. ⭐ **A threshold has TWO edges and the instrument that justifies one is
   usually MONOTONE in it, hence blind to the other** (vv #24(d), zero-set).
   Parseval here is `ratio = 1 − (cᵀG_dropped c)/(cᵀGc)`: *lowering* rcond only
   shrinks `G_dropped`, so a flat `1.000000000` across `[1e-15, 1e-2]` is a
   true reading that carries **zero** information about the lower half. Before
   citing a flat scan as justification, ask which direction the statistic can
   even move in.
3. ⭐ **Then go looking for the guard that already forecloses the other side —
   it is usually there and undocumented.** `[M]` every rcond below `8.696754e-17`
   is REFUSED at construction by the *pair-consistency* guard
   (`_DENSE_METRIC_PENROSE_RTOL`, `|G⁺GG⁺−G⁺| = 7.9e+16` at `1e-18`), so the
   corrupt band is **unreachable, not merely distant**, and `‖G⁺‖₂` is
   bit-constant across the whole legal window. Bisect the boundary; do not
   reason about it. (The same guard is *blind* upward — `G⁺GG⁺ = G⁺` holds for
   a truncated pinv too.)
4. ⚠ **Read which DECOMPOSITION the library actually cuts on, and re-measure
   the residue several ways.** `pinv(hermitian=True)` cuts on `eigh`; the
   shipped comment quoted `svd`. `[M]` the SAME matrix's largest residue reads
   `8.70e-17` (eigh) / `9.71e-18` (eigvalsh, svd-hermitian) / `2.27e-17` (svd)
   — **9.0× spread** — and the bisected refusal boundary equals the *eigh*
   figure to 7 s.f. A round-off number has no stable value; that, not distance,
   is the real argument for a wide margin.
5. ⭐⭐ **A threshold validated on ONE flagship fixture is a claim about that
   fixture's spectral GAP. Census the shipped grid.** `[M]` at the shipped
   rcond: **31 of 105** slab `(order, L)` rows affected (20 raise, 11 breach the
   Parseval gate's own `rtol`), min affected `L = 3`, including the DEFAULT
   `gauss_legendre(16).angular_frame(4)`; **0 of 196** 3-D rows. And the
   mechanism was the opposite of the expected one: **the LIVE modes descend to
   meet the pin**, they are not met by a rising kernel — a 1-D chart mints
   `~L²/2` phantom columns against a fixed node count, the odd-`m` ones carry a
   non-polynomial `√(1−μ²)`, and the spectrum fills in continuously. Where there
   is no gap, **no cutoff is right** and the repair is upstream.
6. ⚠ **Method, twice-bitten: my own instruments died silently and a positive
   control caught both.** (a) A docutils RST check reported "no warnings" on
   deliberately broken input — an unflushed `warning_stream`. (b) A spliced-file
   import raised inside `dataclasses` because the module was not registered in
   `sys.modules` — **my harness, not the candidate**; reporting it as a defect
   would have impeached a correct block (plan-authoring §4's VERIFY clause).
   Also useful: sizing the blast radius of a doc claim needs the *rendering*
   question first — `[M]` `metric.py` has **no `automodule`**, so `-W` cannot
   gate any of this prose at any severity.

Full record + 17 probes: `scratch/rcond_rederivation.md`,
`scratch/probe_rcond_{01..17}_*.py`, block `scratch/rcond_comment_block.txt`.

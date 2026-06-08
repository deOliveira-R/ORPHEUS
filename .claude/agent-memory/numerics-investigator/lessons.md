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

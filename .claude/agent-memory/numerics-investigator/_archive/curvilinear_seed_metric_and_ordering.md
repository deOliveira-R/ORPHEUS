# Archive — the curvilinear seed taxonomy, the GLob quadrature study, the psi-half block metric, and the cylindrical level ordering

Moved COLD from `lessons.md` 2026-09-21 (verbatim: L14 lines 354-393, L16 lines 440-468, L17 lines 469-519, L18 lines 520-579).
Digest successors: L14 / L16 / L17 / L18 (imperatives only). Records: `curvilinear_inverse_seed_taxonomy.md`,
`glob_vs_gl_spherical_quadrature_study.md`, `starting_direction_metric_gauge_derivation.md`,
`cylindrical_level_ordering_symmetry_adjudication.md` (all four TRACKED agent-memory topic files).

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


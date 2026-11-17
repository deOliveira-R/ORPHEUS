# Numerics Investigator — Lessons

Behavioral digest: **what mistake was made, and what correction changed how I work.** Measured war
stories are cold under `_archive/` and in the records each entry names. Nothing here restates a
rule, a preloaded skill or my definition; where one carries a correction, the entry is a stub
naming the clause.

**L-numbers and their `·n` sub-points are STABLE identifiers, never renumbered or reused.**
`numerical-bug-signatures` cites "lesson L10" (Signature 8) and "lesson L11" (Signature 9). An
entry or sub-point that leaves the digest keeps its number as a stub here and its full text in
`_archive/retired_digest_entries.md`.

## Retired IDs — where each correction now lives

The three spine meta-lessons (M1, M2, M5) and the 2026-09-22 retirees are carried by my
definition (`AGENT.md`, the generated block) or a preloaded skill; full text verbatim in
`_archive/retired_digest_entries.md`.

- **M1** (ask what KIND of question first: spectrum, theorem, closed form, transfer gain,
  symmetry) → definition §2 "Name the question before the cascade"; difference operator → §4.
- **M2** (a claim names a LIMIT; sweep the axis it is exact in) → definition §3 step 7 "Refine
  along the right axis"; exact-regime enumeration → `vv-principles` Mode 7.
- **M3** (converged-but-wrong is usually the solver) → definition §3 step 1 "Characterise" and §2
  "A hang or a timeout".
- **M4** (a degenerate fixture annihilates a bug class) → `vv-principles` Mode 7 ("declares which
  terms its ansatz activates and which it nulls") and Mode 12 ("found by degenerating the fixture").
- **M5** (my own instrument was wrong first; budget every probe for its control) → definition §3
  ("every probe budgets for its positive control") and `instrument-doctrine` X1.
- **M6** (two contradicting measurements: one instrument is wrong for the question) → definition
  §4 "Two of your measurements contradict".
- **M7** (a published claim carries its analysis's scope) → definition §3 step 1 ("read the
  paper's own approximation level first"); the instances stay in L4, L21·5–6, L25·5–6.
- **M8** (the reference is the first suspect) → definition §3 step 1 ("A reference that calls the
  kernel under test is contaminated") and `vv-principles` #6; instances stay in L3, L16·3, L22.
- **L1** (run the cascade in order) → `probe-cascade` § Anti-patterns "Don't skip probes".
- **L2** (curvilinear redistribution is the prime suspect) → `numerical-bug-signatures` Signature 1.
- **L6** (a curvilinear matvec is verified only on a NON-FLAT per-ordinate reference) →
  definition §3 step 6; record `cyl_matvec_twin_path_signatures.md`.
- **L9** (a "hang" may be fixture cost; bound the solver first) → definition §2 "A hang or a
  timeout"; record `sn_keff_hang_was_eager_registry.md`. Not carried: **if the bounded solver is
  fast, reproduce the slow path on a fast proxy and bisect THAT, not the test.**
- **L11** (a ρ-honest stop measures `r = Aψ − q`, not `‖Δψ‖`) → definition §3 step 3 and
  `numerical-bug-signatures` Signature 9. Where the diagnostics live now: `IterationRecord`
  (`increment_norms`, `contraction_ratios`, `true_error_estimate()`) in
  `orpheus/numerics/convergence.py`; the per-ordinate balance map is `AngularResidual` in
  `orpheus/transport/residuals/angular_residual.py`.
- **L12** (measure a localized defect's ORDER at the cell) → `vv-principles` #5's check.
- **L13** (greedy `(Ellipsis, *idx)` on a trailing axis) → `numerical-bug-signatures` Signature 11.

---

## L3: A rank-N closure claim needs the SAME SIGNED error across ≥2 quadratures

A single quadrature can fake rank-N convergence: the rank-N error and the F.4 baseline error
cross zero at *different* points along the quadrature-refinement axis, so one snapshot picks one
sign for rank-N and the opposite for F.4 and reports a false structural win. **Magnitude
agreement is not enough — the SIGN must agree.**

- Gate before shipping any closure: `assert_rank_n_structural_win(...)` in
  `tests/gates/cp/test_peierls_rank_n_protocol.py` (it refuses fewer than two quadratures).
- Protocol tracked by **#123, OPEN** (re-verified 2026-11-16). The three falsified closure
  directions with their structural reasons: `_archive/direction_c_pca_rich_adaptive.md`,
  `_archive/direction_q_lambert_marshak_derivation.md`, `frame_5_qmc_quadrature.md`.
- `vv-principles` #13 is the general form (a finite sample is not the population); here the sign
  crossing is what makes it bite.

## L4: Convergence-rate fingerprints discriminate failure modes

A solver stalled at 1–10 % is not enough information. Run an n-doubling sweep and read
`error × n^p`:

| Empirical pattern | Failure mode |
| --- | --- |
| `err·sqrt(n) ≈ const` | Schneider C^(0,1) endpoint singularity in the *solution* (graded mesh repairs) |
| `err·n/log(n) ≈ const` | log-singular kernel diagonal truncation in the *operator* (product-Nyström repairs) |
| `err·n² ≈ const` | smooth-integrand quadrature undersampling — just bump n |
| `err·n^p, 2 < p < 4` | Simpson on a piecewise-quadratic trial; rate limited by solution regularity |
| `err`-ratio not monotone in n | eigenvalue-iteration tolerance, not the discretisation |

**The lesson is the attribution, not the table.** ERR-036 was mis-attributed by the literature
memo to a Schneider endpoint singularity when the fingerprint was `err·n/log(n)`; the
recommended fix happened to be right and the *justification* wrong. **A wrong attribution picks
the wrong fix even when both fixes help.** (`numerical-bug-signatures` Signature 6 carries the
first two rows and the fix family.)

## L5: Read the paper's own stated approximation level before assuming a code bug

ERR-038: a 5 % gap was opened as "a singular limit needing multi-day asymptotic analysis", and
Atalay 1997 states on p.236/p.246 that its Tables 2–5 are first-order Fredholm approximations
with degraded precision at small thicknesses. It was the paper's own floor. My definition §3
step 1 carries the uniform-offset / physical-parameter / numerical-parameter table; the row it
does not carry: **insensitive to every numerical knob and exact at one structural limit ⟹ a
fixed approximation level.**

- **Always carry a moderate-parameter consistency check**: where the paper's approximation is
  tight the solver MUST agree to machine precision, or the verdict has no independent ground.
- **An upstream precision defect does NOT propagate at constant magnitude** (`[M]`
  1.2 % → 0.3 % → <0.1 % down one chain): run the sensitivity chain before claiming an upstream
  fix closes a downstream gap. Record: `atalay_r099_paper_floor_2026_05_03.md`.

## L7: Lift a direction-dependent per-ordinate moment to the GLOBAL frame before the angular reduction

ERR-061. A moment produced in the per-ordinate SWEEP frame and summed by `φ̂ = Σ_n w_n ψ̂_n` must
be sign-corrected to the GLOBAL frame for backward ordinates first — else forward and backward
slopes CANCEL and the scheme loses the diffusion limit. The fix is a per-octant moment-frame
involution at the producer–consumer seam.

- **The angular reduction is the discriminator**: a per-ordinate convention bug is invisible
  until a quantity is summed across ordinates of OPPOSITE sweep direction (a same-direction
  reduction is a degenerate fixture for it).
- ·2 Build a structurally-independent kernel from scratch; reproducing the wrong value bit for
  bit localises the bug to the shared math → definition §3 step 4.
- ·3 A literature-consistent scheme with a faithful-but-inconsistent implementation means a
  CONVENTION between two correct pieces → definition §3 step 4 ("a convention at a
  producer–consumer seam"). Record: `issue_240_d5b_s3_diffusion_limit.md`.

## L8: The project's own theory page can be the contaminated reference

ERR-063. The "a doc plus its own docstring are ONE source" correction is `vv-principles` #6's
check, and the sibling-of-a-different-kernel-family ground is in the same check. Two probes that
no clause carries (restored 2026-11-16; the 2026-09-22 retirement claimed full coverage):

- ·1 Two independent grounds (Hébert 2009 Eq. 3.57/3.58; a sibling solver of a DIFFERENT kernel
  family) → `vv-principles` #6's check.
- **A masked O(1) bug needs a spatially VARYING field to surface**; the discriminating probe
  MUTATES a causally-irrelevant region's value and asserts the answer is invariant.
- **Grep every sibling ASSEMBLY module and classify each same-bug /
  different-structure-correct** — the correct sibling is often the structural ground you need.

## L10: "Error grows with refinement" + an unread library info-flag = an unconverged inner solve

The symptom, mechanism, discriminators and probe are `numerical-bug-signatures` Signature 8.
Two sharpenings it does not carry: **derive `restart` from the operator's domain size** (a
hardcoded cap is the ERR-004 magic-constant anti-pattern, `coding-elegance` Pattern 7); and
**an L1 anchor can pass by NUMERICAL COINCIDENCE at one mesh size**, so **pair every analytical
anchor with a mesh-refinement leg**. Record: `krylov_restart_truncation_bug.md`.

**L10b — recorded-but-unread gives the same false fingerprint as a discarded flag.** The
tolerance-insensitivity tell and "read the convergence record against the budget FIRST" are
Signature 8's ⚠ bullet and definition §3 step 1. What they do not carry:

- ·1 (a tolerance sweep discriminates only if one tol is looser than the capped residual; read
  the record against the budget first) → Signature 8 ⚠ and definition §3 step 1.
- The hole is the CONTRACT, not the flag: a best-effort and a certified answer returning as the
  same type from the same call leaves the class open.
- **The shape discriminator that kills "per-ordinate discretization bias" dead is a `max_iter`
  sweep of the error MAP** — a bias is fixed; an undecayed mode decays geometrically at constant
  ρ with its shape preserved.
- **Audit the blast radius with an in-process pytest plugin** wrapping the solver entries and
  printing every `converged=False` at `sessionfinish`; in one run it found the red gate's SIBLING
  riding the same truncated exit while GREEN on a looser rtol — a latent false green.

**L10c — an all-reflective zero-leakage pure absorber is SI-HARD, and the cost EXPLODES with
dimension** (`[M]` d=1 32 → d=2 258 → d=3 1631 sweeps; one vacuum face collapses d=3 to 208).

- **Budget law: `Σ_t · n_inner` is invariant**, so a FIXED `max_inner` default cannot serve d=1
  and d=3 at once — derive it (`coding-elegance` Pattern 7).
- **A rate claim measured at one dimension can INVERT at the next**: boundary-G-S is 2.5× faster
  than Jacobi at d=2 and 1.95× slower at d=3, so a schedule gate tuned at d=2 is unverified at
  d=3 (L19).
- **ψ flat does NOT imply the TRACE is flat** — "homogeneous nulls redistribution" has a twin:
  *flat cell averages null the face mode*. The mechanism is L19·2 / L23 / L24.

## L14: A curvilinear `(L+C).solve` is NOT uniformly a SweepOperator — the verdict is per (geometry × quadrature)

Cold: `_archive/curvilinear_seed_metric_and_ordering.md`; record
`curvilinear_inverse_seed_taxonomy.md`.

- ⚠ **The behavioural lesson is the CORRECTION of my own headline**: it read "direct inverse for
  slab + CYLINDER", OVER-GENERALIZED from a single `level_symmetric` probe, and the mechanism I
  published was a MIS-ATTRIBUTION. **A per-family verdict needs one probe per family member**
  (`vv-principles` #13), and a mechanism story that fits one member is the cheapest thing to get
  wrong.
- The SOLE lagged element is the M-M half-angle starting seed per level, and it creates a LOCAL
  CYCLE the sweep breaks by lagging. **Seed-dependence is a FORMULATION choice, not intrinsic** —
  the `μ = ±1` equation is CLOSED.
- **The MMS ladder is exactly blind to all of it**: every curvilinear MMS ansatz is ≤ linear in μ,
  the seed's exact regime, so SI converges O(h²) over the whole ladder while the seed-iteration on
  a genuinely higher-order-in-μ field DIVERGES to NaN.

## L15: A SEED/closure re-pose is ruled by sweeping ANGULAR N; a grown Krylov composite must resize `restart`

`[LANDED a29ab2d]`. Cold: `_archive/krylov_composite_restart_and_stale_references.md`.

- **(a)** The N-sweep discriminator → definition §3 step 7; "never dress a structural win as an
  accuracy win" → definition, Return. (The founding case, where the OLD seed was closer at low N,
  is in the cold file.)
- **(b) ⭐ Any carve that ADDS a block to a Krylov composite MUST re-derive `restart`/`n_dof` from
  the composite `to_flat`, not the bulk formula** (ERR-053 family): `[M]` 160 → 210 with `restart`
  left at 160 stagnated GMRES (info=300, 868 s, keff wrong under a bounded outer cap), and
  `restart=210` gave info=0 and 45× faster. **Grep every `restart=` / `n_dof=` against the ravel.**

## L16: To compare two angular QUADRATURES that differ in the pole/seed treatment

Cold: `_archive/curvilinear_seed_metric_and_ordering.md`; record
`glob_vs_gl_spherical_quadrature_study.md`. The recipe, in order:

1. **A standalone scheme-faithful driver, NOT a production point-swap** — when the quadrature sets
   the angular DISCRETIZATION a point-swap hits SINGULAR coefficients (a node AT `μ = −1` makes
   `τ_0 = 0` and the recurrence divides by zero, so `μ = ±1` MUST be a straight characteristic).
   Reimplement the EXACT production closure, verifying every coefficient at `file:line`.
2. Gate the driver bit-faithful to production on a non-degenerate fixture → definition §4 "A probe
   driver". (Here: NON-FLAT vacuum spheres; flat or homogeneous-reflective is degenerate.)
3. **Reference = fine-N in BOTH families plus a cross-family CONTAMINATION GUARD**, and report
   both error-vs-reference and the reference-free matched-N difference.
4. **MMS is blind** (L14): anchor with a closed form (a hand-derived `k_inf`) and `φ = Q/Σ_t`
   streaming equilibrium instead.
5. **Validate the new pole handling with the per-ordinate flat-flux residual** — the
   angle-integrated φ is degenerate to it.

## L17: A STATE DOF's Hilbert metric is not its angular-integration weight

Cold: `_archive/curvilinear_seed_metric_and_ordering.md`; record
`starting_direction_metric_gauge_derivation.md`. The ψ½ "ghost metric" `G_sd ≡ 0` was justified as
the angular through-flux coefficient `(1−μ²)|_{μ=±1} = 0`, confusing the angular-INTEGRATION
weight of the `μ = ±1` ray (correctly zero) with the STATE metric of a discrete DOF. **A metric is
fixed by the DOF's ROLE in the operator algebra, never by an angular weight.**

- **`G_block = 0` is the ONE forbidden value, and it is WORSE than blind**: a nonzero-seed
  reciprocity probe BREAKS on the production path because `A.H` severs the seed. It looked correct
  only because the gate fed a present-but-ZERO seed. **A zero-weight block is not a conservative
  default.**
- **Closing a Mode-12 blindness needs TWO changes**: install the non-degenerate metric AND feed
  NONZERO block data in the gate — and a test that POSITIVELY PINS the blindness must INVERT
  to assert the flip reds.
- **The dense-probe check** (materialising by unit-vector probing is definition §4): compare `T`
  against **numpy's** `Aᵀ` (a structurally-independent ground, not the operator's own machinery).
  `T == Aᵀ` ⟹ gauge-free up to SPD and reciprocity cannot choose it (the general theorem is
  L26·4); if `T` carries a weight, the metric is PINNED.

## L18: Adjudicate a LABELING/ORDERING degeneracy with the operator's own SYMMETRY GROUP

#326 `[CLOSED; remediation LANDED dde93b64 — a level's order is the FIBER's.]` Cold:
`_archive/curvilinear_seed_metric_and_ordering.md`. Carried by the definition: declare the
ansatz's invariants first and adjudicate with a symmetry both problems have (§2 "Ordering or
labelling"); the `xfail` row beside an un-`xfail`ed sibling, and rewiring a re-implementing
driver to call production at promotion (§4 "A probe driver"). What it does not carry:

- **A companion ansatz outside the symmetric sector SEES the defect and still cannot ADJUDICATE.**
- **A within-tie permutation cannot move a cumulative-sum coefficient**, so the whole effect is a
  LABELING and the defect magnitude is ordering-invariant.
- **The leak path is any coupling that maps ACROSS the degenerate class, and it need not COMMUTE
  with the tie-break — grep every cross-class index map before calling a relabeling inert.**
- Verdict here — **no ordering is correct, the CLOSURE is broken.**
- **The constructive exit is to fold to the FUNDAMENTAL DOMAIN**: a degeneracy in a sort key
  usually means the discretization carries a redundant symmetric copy, and on the independent half
  the key is monotone and every competing criterion coincides.
- At promotion, mutate through the test's OWN import binding (`vv-principles` Mode 11).

---

## L19: An iterative-solver RATE question is a SPECTRUM question

#341 `[CLOSED; the octant-order lever is #343, OPEN]`. Cold:
`_archive/rate_spectrum_certificate_and_angular_axis.md`; record `issue_341_boundary_gs_rate.md`.
"Build the iteration matrix; never re-time" is definition §2; "a production default branches only
on a variable varied causally" is definition, Return (here `ndim` was falsified on both sides, and
thickness, mesh, aspect ratio, quadrature order and `c` all move the sign at fixed `ndim`).

1. **The recipe when a dense build is too big:** any `x ← A_inv.apply(Σ gᵢ.apply(x))` driver IS a
   linear operator; wrap it over the composite's `to_flat`/`from_flat` as a
   `scipy.sparse.linalg.LinearOperator` and `eigs(which="LM")`. It is **immune to the stopping
   test**, so no ρ-blind stop or `max_inner` truncation (L10) can contaminate it. Controls:
   reproduce the ρ FITTED from a real residual history, plus `G(2x) = 2G(x)` and `G(0) = 0`.
2. ⭐ **Before hunting for why a splitting comparison inverted, ask whether the THEOREM that forbids
   it still applies.** Varga's comparison makes the inversion IMPOSSIBLE for a non-negative
   iteration matrix, so an observed inversion is evidence the operator is NOT non-negative, and the
   productive question is *which term is negative and why*. Here the multi-D diamond face
   transmission carries **`d−1` eigenvalues exactly `−1`** — an undamped zero-cell-average face
   sawtooth invisible to `Σ_t V ψ_c`, growing with `ndim`. **Read any all-reflective DD rate
   pathology through that spectrum first.**
3. **A per-axis SIGN is usually a gauge** (a diagonal similarity leaves both rates invariant): sign
   *indefiniteness* voids the theorem, the sign *pattern* explains nothing. And **a model that fails
   to reproduce the effect is worth as much as one that does — it deletes a whole hypothesis class.**
4. ⭐ **Enumerate a finite design space instead of sampling it**: all `8!` octant orders collapse to
   **25** patterns, and measuring all 25 gave an exact separating law
   (`LOSES ⟺ max_a L_a > Σ_{b≠a} L_b`, 25/25) plus a 2.5× rate spread with the shipped order 24th
   of 25. A sampled sweep would have produced a fitted story instead of a law.
5. **Ask whether the two arms are racing the SAME mode before calling a change a "flip"** — extract
   the dominant eigenvector and report where its mass sits; at d=2 both raced one face, at d=3 they
   raced different faces, so it was two different comparisons.

## L20: A RESIDUAL cannot gate an EIGENVALUE contract — measure the TRANSFER GAIN first

#340 N5, **REFUTED** on 38 solves / 8 geometries / 3 mixtures. Cold:
`_archive/rate_spectrum_certificate_and_angular_axis.md`.

1. The transfer gain `|Δy| / X`, measured before any threshold → definition §2 and
   `instrument-doctrine` X1 "Transfer gain". What they do not say: this is `vv-principles` Mode 12
   read in the MIRROR — the gate is SIGHTED on a class the CONTRACT is blind to. **Cure: project
   onto the functional the contract reads.**
2. ⭐ **A SIGNED projection against an APPROXIMATE weight is worse than no weight** — it manufactures
   accidental near-cancellations, i.e. false NEGATIVES. The cheap flat adjoint was itself *verified*,
   and correct for the WRONG PROBLEM. Pay for the real adjoint or use the unsigned norm.
3. **Answer the NULL case before the discrimination question, with a TWO-LEGGED tolerance sweep** —
   the outer leg moved the "pass" value 6 decades and the inner 0.2 %, so it was the OUTER's
   increment-stop slack, not a floor; one leg alone mis-anchors the whole study.
4. **When lifting a production certificate one level, the CONSTANT does not come with it — and
   copying `record.binding_criterion.tolerance` silently picks the LOOSER criterion.** A residual bar
   scaled by an INCREMENT tolerance is a category error twice over.
5. The mixture consistency identity `σ_t = σ_c + σ_f + Σ_to SigS[0][g,:]` → definition §3 (every
   dispatch) and `vv-testing` (only for work under `tests/`). Founding case: a `[to, from]`
   `sig_s` read as `[from, to]` made a brief's "benign" reference 30 % off.

## L21: An ANGULAR-consistency claim is separated by `h → 0`, not by the parameter it is named after

#319 / #235 (both OPEN for later phases, re-verified 2026-11-16), 251 solves all `converged`.
`[Phase-0 gates LANDED a3121cfe]`. Cold: `_archive/rate_spectrum_certificate_and_angular_axis.md`;
record `issue_319_flux_dip_discriminator.md`.

1. The axis choice → definition §3 step 7. The number that founded it: sweeping optical thickness
   at fixed cells-per-mfp separated the two τ schemes **not at all** (fitted decay rate 0.000 for
   BOTH); refining `h` at fixed physics separated them without bound (3.2× → 204×).
2. **A regime sweep at fixed `c` self-destructs** — `Σ_a·R = Σ_t·R(1−c)` grows with it and every
   scheme agrees for a reason unrelated to the question; use the ε-scaling `Σ_t = 1/ε, Σ_a = ε,
   Q = ε`. (`vv-principles` #24(e) carries the regime check; the half it does not is ⟹ **carry a
   fixture-LIVENESS column** that declares when the fixture stopped posing the question.)
3. ⭐ **Build a λ-CONTINUUM through the two candidates, not an A/B** — `τ(λ)=λτ_A+(1−λ)τ_B` turns "A
   beats B" into "is A the MINIMISER?", and `λ_opt(h)` becomes a falsifiable curve. Sphere
   `λ_opt → 1` on two instruments ⟹ the shipped τ is the family optimum; cylinder `λ_opt → 0.73`
   with the two instruments DISAGREEING ⟹ no optimum claimed — **the disagreement is the finding.**
4. ⭐ **A theory scalar can be τ-loaded or τ-blind depending on which EDGES you feed it, and the
   blind version is the natural one to write** (from the standard weight-partition edges M&M's β is
   τ-blind *by construction* — that substitution IS their β=0 proof). ⛔ And it is identically zero
   for BOTH schemes on a folded cylinder: **a spherical invariant does not transfer to a geometry
   whose angular derivative is in a different variable.** (Main-lessons L45 is the same β annihilated
   by a σ_y fold: **feed any symmetry-suspect analytic instrument deliberate garbage in the varied
   slot.**)
5. ⭐ **A literature diagnostic transfers between geometries only in its LEVEL-LOCAL form** — rebuilt
   from the level's own azimuthal moments it reproduces the published formula on the sphere
   bit-for-bit and gives sane cylinder values, where the global form reads a `+2.76` artefact. ⚠ It
   is an S2/S4-class instrument, so at S8/S16 genuine curvature dominates and it reports a bias.
6. **The benefit of a low-order-consistency fix DECAYS with angular order and can invert** (`[M]` 14×
   at S2 → 0.9× at S8, the principled scheme measurably *worse*). **An accuracy comparison run only
   at high N will report the principled scheme as a regression — correctly, and for a reason that is
   not a bug.** Same family as L25·5.

## L22: A frozen reference is stale by MAGNITUDE CLASS, and the nulp count tells you neither

Triaging 9 "bit-identity" reds — all 9 stale references, none a regression, in two magnitude classes
the failure messages could not distinguish. Cold:
`_archive/krylov_composite_restart_and_stale_references.md`. (`numerical-bug-signatures`
Signature 10 carries the sibling-pass discriminator and the re-baseline rule; `vv-principles` #25
the re-baseline radius by reference KIND.)

1. Report `max|a−b| / max|b|` before the nulp count → `vv-principles` § Bit-identity ("print
   `max|a−b| / max|b|` BEFORE reading the nulp count", with this case's `[M]` 1.04e+15 nulp = 8 %
   and the 1-ULP `array_equal → False`).
2. **A gross move in the SCALAR flux refutes an ordering hypothesis on sight** — `φ = Σ_n w_n ψ_n` is
   permutation-invariant to `N × ULP`, so a percent-level move means the rule's VALUES moved.
3. A tolerance pin cannot warn a bit-identity consumer; the byte-level fingerprint →
   `vv-principles` § Bit-identity ("a producer's TOLERANCE pin can never warn its consumers'
   BIT-IDENTITY pins").
4. ⭐ **A gate whose reference is a different FP ASSOCIATION of production's own expression is a coin
   flip, not a contract** — `[M]` `src + (A + B)` against `(src + A) + B`: 68 of 80 slots
   bit-identical, 12 differ, max ULP 1, and the passing sibling passes by index luck.
   **Re-associate the reference or demote the assertion.**
5. The two probe-harness failures (pytest ANSI codes before `FAILED`; `python -c` importing the main
   tree) → `code-search` ("pass `--color=no`"; "`python -c` puts the CWD at `sys.path[0]`").

---

## L23: A discrete operator's SINGULARITY is a two-question object — refuse the either/or

#344 `[CLOSED; LANDED f934ff57/b51bc802; characterization gates 1a2be025.]` Cold, with every
measurement: `_archive/issue_344_singularity_kernel_and_gauge.md`.

1. ⭐⭐ **When a claim is "the null space IS class X", the decisive number is `dim ker` MINUS `|X|` —
   and you must MEASURE `|X|`, not assume it exists.** One line (`min|Ω·n|`) refuted the brief's
   framing before any solve, because a level-symmetric rule cannot produce `Ω·n = 0` at all; and
   `dim ker = |X| + R` then held exactly on 9 rows, which turned an either/or into a DECOMPOSITION.
2. Dense SVD through the production builders, reported with its singular-value gap → definition §4
   "Materialise the operator". (Here ARPACK's `eigs(k=12)` bound "3 at d=2, ≥6 at d=3" was really
   12 and 138.)
3. ⭐⭐ **Two blindness mechanisms that look alike are told apart by the METRIC, and that difference
   decides the REMEDY**: a slot whose `G` is exactly zero can never be seen by any G-weighted
   functional ⟹ typing it away (`coding-elegance` Pattern 4) is the only fix; a rank deficiency with
   `G > 0` CAN be gated, and typing cannot remove it. **Always ask "is this class in `ker G`, or
   merely in `ker` of the functionals I happen to gate with?"** — `vv-principles` #18 covers only the
   first, and the second is the commoner case.
4. ⭐ **A residual stop and a conservation projection are blind to `ker A` BY CONSTRUCTION, so the
   only informative half of that measurement is the POSITIVE CONTROL** (`A(ψ+αv) − q ≡ Aψ − q` is a
   theorem, not a finding). Budget the probe for the control; the "unmoved" column is free.
5. ⭐ **A converged solver's deviation from the analytic answer is in `ker A` EXACTLY — test it with
   `‖Aδ‖/(‖A‖‖δ‖)`, no null basis needed.** And **identify the recorded scalar before trusting it**:
   the memo's headline was a max over all ordinates while the printed row was one ordinate.
6. ⭐⭐ **Fit the counting law, WRITE THE PREDICTION DOWN, test it OFF-SAMPLE, then SWAP THE SCHEME to
   get the mechanism.** `[M]` 3 of 3 off-sample including a change of quadrature order, preconditions
   measured and not assumed; the mechanism closed with one substitution — **`LinearDiscontinuous` on
   the identical box is NON-singular** ⟹ it is DD's `ψ_out = 2ψ̄ − ψ_in` involution, undamped by
   `Σ_t V ψ_c` (L19·2 end-to-end). **Blast radius worth naming:** the eigenvalue entry DEFAULTS to
   all-reflective, so every `d ≥ 2` Cartesian DD k-solve ran a singular within-group operator.
7. ⭐ **Measure whether the TRUE answer IS the canonical representative before recommending a gauge**
   — the exact solution is the minimum-`‖·‖_G` member, so projecting the iterate off `ker A` is an
   EXACT fix (8.97e-02 → 5.8e-13), not a convention. And **when enumerating what a null direction is
   invisible to, enumerate the moment LADDER; do not reason about it** — I predicted two
   cancellations in a row and `[M]` both were wrong.
8. ⭐ **Carry `‖Ad‖/‖d‖` BESIDE the error, or you cannot tell a frozen null component from leftover
   residual.** `[M]` the deviation is identically zero at even `n` and 6.2e-02 at odd `n`, so a
   4/8/16/32 ladder reports "nothing to see" (`vv-principles` #13's break-the-congruence-class rule);
   inside ONE parity class the law was exact (`err·n` constant to 8 s.f.), i.e. O(h) and **not** the
   wrong limit.
9. ⭐⭐ **To tell a GAUGE FREEDOM from an INCOHERENT solver, REMOVE THE KERNEL and re-run — and remove
   it ≥3 structurally-different ways**, because "boundary moves, bulk does not" is produced by BOTH.
   ⚠ **The obvious control can be a NON-control**: an even-`n_x` box, where the parity finding says
   the mode is absent, still has `dim ker A = 12` — what is absent is the kernel's EXCITATION by that
   source. **Assert `dim ker == 0` INSIDE the control; never infer it from a deviation being zero.**
   (My own probe was wrong first, twice, both times toward the alarming verdict.) The fixed-point
   manifold is `vv-principles` Mode 9's false-RED premise.
10. Full-space `‖MM⁻¹−I‖` versus the RHS subspace the driver supplies → definition §4 ("re-measure on
    the subspace the driver actually feeds"). What it does not say: **a reified forward-substitution
    "inverse" is a SUBSPACE inverse by construction** — fine for SI, a live hazard for a Krylov
    PRECONDITIONER. (`[M]` 3.3e-01 full-space against `‖Aψ*−q‖/‖q‖ = 8e-14`.)

## L24: A discrete operator's KERNEL is usually a CLOSED-FORM problem

#344's basis, found in `0.05 s` where a dense SVD is `23 s` at half the size. Cold:
`_archive/issue_344_singularity_kernel_and_gauge.md`.

1. A counting law independent of a parameter the operator contains is combinatorial; derive it →
   definition §2 "Kernel or counting". The derivation here: setting the degenerate branch
   (`ψ_c = 0`) turns DD's closure into an involution, so every face field is a sawtooth; in the
   balance every cross-section, mesh width, weight and area cancels, leaving *a sum of functions,
   each blind to one coordinate and one sign, vanishing identically*.
2. ⭐⭐ **Sign CHARACTERS diagonalise a specular-BC constraint system** (a specular BC says a quantity
   is blind to one SIGN, the balance that it is blind to one COORDINATE; expanding in the sign
   characters splits the system into one additive-separable equation per character subset). ⟹ **read a
   SUM over axes in a counting law as "the modes live on PLANES (one free coordinate)", never as a
   fit.** ⚠ An orbit count is the number of ordinate ORBITS under the reflection group, NOT ordinates
   per octant — that reading is off by `2^{d−1}`.
3. ⭐ **Where an SVD is unaffordable, the span check is a PRODUCTION-GENERATED kernel vector — with a
   round-off NEGATIVE control.** The control (the other arm's pure-round-off deviation) must read
   fully OUT of span, proving the projector is not a universal absorber; without that leg, "everything
   I test is in the span" is unfalsifiable.
4. ⚠ **Two mechanisms can share a PARITY fingerprint, and only a kernel-CONTENT measurement separates
   them.** My detector-blindness hypothesis was true (a uniform detector is exactly blind to a
   `(−1)^{i_⊥}` profile at even cell counts) and NOT the cause: at even `n_x` only 15–31 % of the
   deviation is in `ker A` ⟹ the mode is **ABSENT, not hidden**. **Measure `‖P d‖/‖d‖`, not `‖d‖`.**
5. ⭐ **A blindness LIST measured on ONE quadrature is a sample, not a population.** `[M]` the prior
   "every `|Ω·n|^p`, p = 0..3, is blind" held only where the tangential component is 0; on
   `lebedev(11)` the **`p = 0`** moment reads the modes at 2.99e-02. Honest condition: **mirror-EVEN in
   angle AND ≥ 1 power of `|Ω·n|`** (a CURRENT-type functional). And matching matters — a
   `sign(μ_xμ_y)` weight is BLIND where `sign(μ_x)` sees 4.4e-2: **"angularly resolved" ≠ "sighted"**.
6. ⭐ **A DENSE basis is not the shippable form of a STRUCTURED nullspace** (disjoint supports per orbit
   and group make the Gram block-diagonal: 17.6 GiB → 154 MiB, apply 12 ms). ⛔ And
   `ker G ∩ ker A ≠ 0` is a real hazard: a bit-zero `G` makes `BᵀGB` SINGULAR and a `sqrt(G)`-QR give
   `0/0`, so **there is NO minimum-norm gauge for those directions** — project on the G-positive
   component only.

## L25: An ARITY / "does this contract widen?" question is a THEOREM question

Curvilinear LD × Morel–Montry τ; 106 SymPy checks, 0 failures. Cold:
`_archive/curvilinear_tau_ld_and_gram_ownership.md`.

1. "Ask what the defining conditions commute with" → definition §2 "Contract or arity". The
   argument it does not spell: **a scalar convex combination commutes with every linear map**, so a
   quantity defined by membership of a convex set plus exactness of a scalar blend gives the SAME
   scalar statement in every component of every linear representation, and the widened form is an
   **overdetermined system whose every row returns the same value**. Three cheap hypotheses to
   check: the scalar is independent of the widened index; the projection is linear; the set is
   convex. Minutes, against a multi-day expansion.
2. ⭐⭐ **The asymptotic expansion and the POSITIVE CONE answer different questions, and the expansion
   is structurally BLIND to the cone** — a sign-alternating cell-to-cell mode is EXCLUDED BY THE
   ANSATZ, which is why Palmer–Adams carry "limits to a *stable* diffusion equation" as a SEPARATE
   criterion. `[M]` the transmission sign ladder: DD flips at `τ_opt = 2`, bare LD at 3, lumped LD
   never. ⟹ **when a scheme is "verified in the diffusion limit", ask WHICH of the two it was.**
3. ⭐ **Where a coupling operator is a TENSOR PRODUCT with disjoint index sets, EVERY functional of it
   FACTORS — grep the free symbols, it is a one-line proof.** Corollary that pays: the leading-order
   discrete diffusion equation then carries no angular parameter, so **a SPATIAL defect cannot be
   repaired by the angular knob and vice versa** — two orthogonal failure modes, and conflating them
   is the whole confusion.
4. ⭐⭐ **When one knob is exact by construction, price the OTHER knob's error in the first knob's units
   — the ratio is usually the finding.** `[M]` with τ exact, a starting-cosine error of 1.6 % at S4
   reproduces the ENTIRE contamination τ exists to remove, because the contamination is EXACTLY affine
   in the starting cosine and its coefficient GROWS with N (24× → 333×). **The celebrated knob was
   never the risk; its neighbour was.**
5. ⚠ **A published RULE OF THUMB carries the order it was derived at.** M-M's "the dip is eliminated as
   long as the starting flux is not seriously UNDERestimated" is a claim about a derivative that `[M]`
   **flips sign between N=2 (their own test case) and N ≥ 4** — the safe DIRECTION inverts, while the
   magnitude falls 5 orders, so the stakes collapse as it flips. **Evaluate the quantity, never the
   heuristic** (same family as L21·6).
6. ⚠ **"Lumping" names ≥3 different operations and they disagree.** `[M]` Legendre-diagonal lumping
   breaks the per-moment-row flat-flux identity; **nodal ROW-SUM lumping (what Palmer–Adams's FL
   actually is) preserves it exactly**; row-sum lumping the GRADIENT gives transmission identically 0,
   a degenerate scheme. ⟹ **a ⛔ banner saying "X may not be lumped" must name the BASIS and the
   MATRIX** — mine, inherited, condemned the wrong operation. The real freedom: the infinite-medium
   identity pins only the row sums, leaving one parameter per row, so **the accuracy/positivity trade
   is a CHOICE OF THAT PARAMETER, not a property of "lumping".**

## L26: An OWNERSHIP question is answered by TWO measurable discriminators

The adjoint/Gram ownership audit of the 1-D curvilinear SN streaming path. Cold:
`_archive/curvilinear_tau_ld_and_gram_ownership.md`.

1. Measure the increment's structure; build the difference operator → definition §2 "Ownership" and
   §4. The trap it does not name: ⛔ **"does X influence the adjoint?" is the WRONG question when the
   adjoint is a reverse-mode VJP** (`apply_transpose ≡ apply.T`, so the answer is "everything").
2. ⭐⭐ **Calibrate the symmetry ratio before reading it**: `‖A+Aᵀ‖_F/‖A‖_F` is **0 for skew, 2 for
   symmetric, √2 when no entry pairs with its transpose partner (triangular-like)**. Without the
   calibration a √2 reading looks like "very asymmetric" when it means *triangular* — a structural
   statement about a MARCH.
3. ⭐⭐ **Two derived constants from the SAME two numbers on ADJACENT lines can sit on opposite sides of
   the self-adjointness split** (one contributes an exactly DIAGONAL block, self-adjoint in every
   metric; the other an increment with zero trace and no transpose pairing), so a design that bags them
   as "the closure constants" fuses a reaction-like scalar with a transport-like coupling. **Diagonal ⟹
   order-free ⟹ may live anywhere; non-diagonal ⟹ welded to the traversal that reads it.**
4. ⭐⭐ **A `⟨Aψ,φ⟩_G = ⟨ψ,A†φ⟩_G` reciprocity gate CANNOT adjudicate the choice of `G`** — with
   `A† ≡ G⁻¹AᵀG` it is an identity for EVERY invertible `G` (verified under Euclidean, random and
   adversarial metrics, with a mismatch control that reads 8.22), so the gate plus its wrong-metric
   control prove **consistency and loadedness, never CHOICE** (Mode 12 with the whole
   invertible-diagonal group as the stabiliser). ⟹ **the adjoint does not pin the metric; the physical
   FUNCTIONALS do**, which is why the metric belongs to the SPACE and the operator is its consumer.
   (L17 is the instance this generalises.)
5. ⭐ **One symbol, four roles, four owners — and the separating experiment is a GLOBAL RESCALE through
   the production constructors, not instance surgery.** `[M]` rebuilding at `w → 3.7w` scales the
   metric exactly and leaves the streaming operator bit-identical (the redistribution sees only the
   ratio `α/w`); the sphere REFUSES the same rescale because its angular-cell partition is the
   cumulative weight; a single-ordinate perturbation is refused by the pole mirror's
   weight-preservation contract. ⟹ metric / scale-free ratio / absolute mesh width / admission
   precondition. **Ask which one a consumer means before moving it.**
6. ⭐⭐ **A "Gram" inside an OPERATOR is usually the mass matrix under a DIFFERENT measure, and naming
   the measure settles the ownership.** `[M]` (SymPy, exact) the curvilinear redistribution Gram is
   `R_kj = ∫ b_k b_j (∇·ê_r) dV`, so `R₀₀ = ΔA` is the **divergence theorem**, not a per-chart
   normalization; it is SPD (a genuine inner product) but `M⁻¹R ≠ λI`, so it is not `M` rescaled. ⟹
   **measure from the CHART, basis from the SCHEME ⟹ the home is the (chart × scheme) pair**, and it is
   a *coefficient*, never a metric. ⚠ When the two axes differ the object is **rectangular** ⟹ a
   PAIRING, not a Gram.
7. ⛔ **Do NOT justify a transport metric as "the one that makes streaming skew-adjoint".** `[M]` the
   ratio is √2, flat over `nx = 4…64` and the **same number on the slab**: a face-ELIMINATED marching
   operator is TRIANGULAR, not skew. Skew-adjointness belongs to the (cell ⊕ face) saddle system, and
   DD substitutes the interior faces out and destroys it.
8. ⚠ **Each structural ZERO needs its OWN control, and the cheapest one is the ADJACENT array entry.**
   `[M]` two knobs read inert on every shipped fixture for two DIFFERENT reasons (one has no carrying
   consumer; one multiplies a seed the ray-decoupled block feeds with zeros), and each was proved live
   by its own control. A naive inventory would call both "not used". (And the instrument bit first:
   assigning into `__dict__` is a silent no-op for a `property`.)
9. ⭐ **The cylinder is BLIND to the Gram question — a THIRD member of the family.** `[M]` the cylinder's
   ratio is bit-exactly the shipped moment-axis metric (so it reads as "not its own object") while the
   sphere's carries an off-diagonal; sharper, **the cylinder's MASS Gram IS the sphere's REDISTRIBUTION
   Gram, exactly** — the two objects the ownership argument is trying to separate are the same matrix
   one geometry over. ⟹ **any curvilinear Gram/measure claim must be witnessed on the SPHERE.**

## L27: A "noise mode" reading of a small singular value is a HYPOTHESIS

`_DENSE_METRIC_RCOND` (`orpheus/numerics/metric.py`) re-derivation (ERR-080). A pinned `pinv` cutoff
was justified by *"one ~1e-16 noise mode … 1e-12 sits ~4 orders above the noise floor"*. Every quoted
number reproduced; the reading was inverted. Cold: `_archive/rcond_threshold_rederivation.md`.

1. ⭐⭐ **Discriminate "round-off residue of an EXACT dependency" from "small real mode" by SOLVING for
   the dependency, then applying the candidate null vector to the RAW TABLE — not to the Gram.**
   Forming `G = AᵀWA` can itself manufacture rank loss, so `‖Gv‖ ≈ 0` is weak; `‖Av‖ ≈ 0` says the
   COLUMNS are dependent as functions on the node set and settles it. `[M]` the closed form fell out in
   two lines and matched the SVD null vector to 2.2e-16. ⟹ **a "noise floor" is a claim about a number
   that should have a NAME; if it has one, it is not noise.**
2. A threshold has two edges and its justifying statistic is monotone → `vv-principles` #24(f).
   (Here the Parseval statistic read a flat `1.000000000` across `[1e-15, 1e-2]`.)
3. Find the guard that forecloses the other side and BISECT its refusal boundary → `vv-principles`
   #24(f). (Here every rcond below 8.696754e-17 is REFUSED at construction by a pair-consistency
   guard: the corrupt band is unreachable.)
4. ⚠ **Read which DECOMPOSITION the library actually cuts on, and re-measure the residue several ways.**
   `pinv(hermitian=True)` cuts on `eigh` while the shipped comment quoted `svd`, and `[M]` the SAME
   matrix's largest residue reads 8.70e-17 / 9.71e-18 / 2.27e-17 across the three routines — a **9.0×
   spread**, with the bisected refusal boundary matching the *eigh* figure to 7 s.f. **A round-off
   number has no stable value; that, not distance, is the real argument for a wide margin.**
5. ⭐⭐ **A threshold validated on ONE flagship fixture is a claim about THAT fixture's spectral GAP.
   Census the shipped grid.** `[M]` **31 of 105** slab `(order, L)` rows affected — 20 raise, 11 breach
   the gate's own rtol, min affected `L = 3`, including the DEFAULT
   `gauss_legendre(16).angular_frame(4)` — against **0 of 196** 3-D rows. And the mechanism was the
   opposite of the expected one: **the LIVE modes descend to meet the pin**, they are not met by a
   rising kernel. **Where there is no gap, no cutoff is right and the repair is upstream.**
6. ⚠ **A doc claim's blast radius needs the RENDERING question first**: `[M]` the module carries no
   `automodule`, so `-W` cannot gate any of that prose at any severity.

## L28: Platform drift is proven by INVERTING the frozen snapshots

Not covered by `numerical-bug-signatures` Signature 10 (stale snapshot) or `vv-principles` #31
(single-draw bit-exact): those name the classes, not the instrument for a platform change with no
old environment to run. `[M]` 2026-10-01, macOS 27.0.1 Accelerate.

1. **The old primitive output is recoverable from the snapshots it fed.** A coordinate ULP search
   over the few free numbers (GL-4: x1, x2, w1, w2) against the frozen walk-matvec `.npy` found a
   unique rule reproducing 496 of 496 elements; substituting it alone turned 20 of 26 reds green.
   That is a bit-exact positive control, which "it is probably LAPACK" never is.
2. ±3-ULP jitter over ≥10 seeds to sort gates by whether noise flips them → `vv-principles` #38.
   (Here the noise-verdict gates were two FP association orders, `(a·y)/a == y`, and a tolerance
   inside its own noise band; re-capturing does not cure them.)

Record: `scratch/platform_drift/memo.md` (untracked).

## L29: A matched-accuracy target read off a GAUGED observable can be one scalar's error, not the reading's

Not covered by `instrument-doctrine` X2 (it asks for both legs' populations, not for the error's
decomposition). `[M]` 2026-10-11: the old reference's phi_g(1.2) erred 1.19e-3 and 1.16e-3 in its two
groups, the SAME sign and size, while the gauge-free ratio phi_1/phi_0 erred 3.1e-5: the "accuracy to
match" was the old solve's production normaliser, not its transport. A rule matching 1.2e-3 can be 40x
worse in the ratio. **How to apply:** before fixing a matched-accuracy target on a gauged quantity
(eigen flux, any normalised reading), print the per-component errors and a gauge-free combination
beside them; if the components share one factor, state both criteria and which rule each selects
(`scratch/characteristic_architecture/p3/step1/errors.log`, untracked).

## L30: Retiring a mesh that resolved two scales re-exposes the second scale's convergence parameter

Not covered by definition §3 step 7 (it names which axis separates a claim, not that an incumbent's
mesh can be doing a SECOND parameter's work for free) nor by `vv-principles` "bit-identity vs
principled-equivalence" (which asks for an independent reference, not for the error's owner). `[M]`
2026-10-12, product integration along the characteristic reference's lines: production at 4 Gauss
points per piece read 3e-16 along lines only because its exponential mesh (260 pieces per slot) also
over-resolved the geometry; on the geometric pieces alone the same 4 nodes left 1.5e-6 in a block
(case 1) and 4e-5 (case 2), and matched accuracy needed n = 8 to 12. **How to apply:** before calling a
mesh-removing scheme "matched at the nominal parameter", hold everything else fixed (same lines, same
groups, same emission) and tabulate the integration error alone against a finer incumbent at n, 1.5n,
2n, 3n; and check the reference is not inside the scheme's exact class (case 2's mpmath reference had a
source constant per region, which interpolation reproduces at any n: green at n = 4 for the wrong
reason). Record: `scratch/characteristic_architecture/p3/product/` (untracked).

## L31: A left kernel is a property of the ROW scaling — derive it from the rows the unknowns enter, not from the right kernel's metric

2026-10-31: the brief's probe tested ker Aᵀ against span(Φ), span(GΦ), span(G⁻¹Φ) (all 1.00) — the
natural guesses when the metric is believed to be the bridge. The answer was neither: writing out the
trace ROWS (unit-normalised reflection and telescoped-DD transmission) and solving the column equations
gave ker Aᵀ = span(JΦ), J = −(Ω·n̂)A_f, a SIGNED weight no metric carries; the left generators are ±1
patterns. No clause covers it: `instrument-doctrine` X4 and definition §4 say how to measure, not that
ker Aᵀ moves by D⁻¹ under a row rescaling D, so a metric-image guess is a category error before it is a
wrong number. Check: before guessing a left kernel, print the rows that the right kernel's DOFs enter
and solve lᵀA = 0 column by column; completeness then follows from rank–nullity on the square A, not
from the fit. Record: `scratch/membership_covector/u0/left_basis/MEMO.md` (untracked).

## L32: A variance hidden by a unit metric is seen by a NON-UNIFORM MERGE's coarse EQUATION, never by a pairing or an adjoint

2026-11-12: on energy G_E = I, so every adjoint-vs-transpose number and every pairing is blind to E
against E\*. And a restriction test is blind too: a swapped assignment (ψ† summed, q partner-averaged)
keeps both pairings exact. The instrument is the coarse balance under a merge with UNEQUAL group sizes:
the summed chart solves A_c to 1e-15, while the mean chart misses by 1.19 (slab) and 10.3 (0-D) through
S and F (n_G′/n_G). A uniform merge is exactly blind (1e-15), so it is the control leg. The gap:
`vv-principles` Mode 12 tells you to find space-side stabilisers by degenerating the fixture, but does
not name refinement or merge as the probe for an INDEX that a unit metric hides. Second trap: a
partner-weighted average zero-fills where its weight vanishes, so a localised source contaminates every
"which weight" comparison; re-measure on a full-support source. Record:
`scratch/true_ends/consult/four_space_probe.md` (untracked).

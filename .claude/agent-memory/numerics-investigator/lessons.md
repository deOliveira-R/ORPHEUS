# Numerics Investigator — Lessons

Behavioral digest: **what mistake was made, and what correction changed how I work.** Every
measured war story is COLD — six files under `_archive/` (moved verbatim 2026-09-21) plus the
records each entry names. Nothing here restates a rule or a preloaded skill; where a correction
has been uplifted, the entry cites the clause and stops.

**The L-numbers are STABLE identifiers, never renumbered** — `numerical-bug-signatures`
SKILL.md cites "lessons L10/L9" and "lesson L11" by number, and
the #282 probe `diag_282_sphere_repose_convergence.py:23` cited "L14/L15" (the probe is retired at `f36572c8` (R19); recover with `git show f36572c8^:<old path>`). A retired
entry keeps its number as a stub naming where its correction lives.

## The spine — three meta-lessons

Read these three every dispatch. Open an `L` entry when its instance is in front of you; open
`_archive/` only when a number needs checking.

- **M1 — Ask what KIND of question this is before measuring; the answer is usually a THEOREM, a
  SPECTRUM or a CLOSED FORM, not a re-run.** A rate question is a spectrum question (L19·1);
  "does this contract widen?" is *what commutes* (L25·1); a counting or kernel law is
  combinatorial, so derive it (L24·1); an ownership question is two discriminators (L26·1);
  "can statistic X gate contract Y?" is one number, the TRANSFER GAIN (L20·1); a labeling
  degeneracy is the operator's own symmetry group (L18·4). **Build the difference operator, not
  the operator** (L26·1). → now in `AGENT.md`, "Before the cascade" (2026-09-21).
- **M2 — A claim names a LIMIT; the separating axis is the OTHER discretisation going to zero,
  and your fixture is probably inside the claim's exact regime.** Angular consistency separates
  by `h→0`, not by the parameter it is named after (L21·1); a SEED is an angular closure, so
  sweep `N` at fixed fine mesh (L15a); every curvilinear MMS ansatz is ≤ linear-in-μ = the
  seed's exact regime, so the whole ladder is seed-blind (L14, L16·4). Enumerating the field
  family a closure is EXACT on, and gating outside it, is `vv-principles` Mode 7.
- **M5 — My own instrument has been wrong FIRST, repeatedly, always in the flattering direction.
  Budget every probe for its control.** The varying knob was unreachable (L18); a grep of pytest
  output read zero failures through ANSI codes and `python -c` imported the main tree instead of
  the worktree (L22·5); a probe "checked linearity" on a matrix that was linear by construction
  (L23·9); a docutils check reported no warnings on deliberately broken input (L27·6); a
  monkeypatch onto a `property` was a silent no-op (L26·8). This is `instrument-doctrine` X1
  applied to *analytic* instruments, and **each structural zero needs its OWN control.**

---

_Retired 2026-09-22 by the agent-definitions audit, each restated by a clause or now carried by the definition: M3, M4, M6, M7, M8, L2, L8, L9, L12, L13. Numbers are not reused._

## L1: Run the diagnostic cascade in order — no skipping

DUPLICATE: `AGENT.md` § Diagnostic Cascade ("Execute in order … Do NOT skip steps") and
`probe-cascade` § Anti-patterns ("Don't skip probes"). Founding cost kept as the teeth: **six
wrong hypotheses** were spent on cylindrical DD divergence by guessing before isolating.

## L3: A rank-N closure claim needs the SAME SIGNED error across ≥2 quadratures

A single quadrature can fake rank-N convergence: the rank-N error and the F.4 baseline error
cross zero at *different* points along the quadrature-refinement axis, so one snapshot picks one
sign for rank-N and the opposite for F.4 and reports a false structural win. **Magnitude
agreement is not enough — the SIGN must agree.**

- Gate before shipping any closure: `assert_rank_n_structural_win(...)` in
  `tests/cp/test_peierls_rank_n_protocol.py` (it refuses fewer than two quadratures).
- Protocol tracked by **#123, OPEN** (verified 2026-09-21). The three falsified closure
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
with degraded precision at small thicknesses. It was the paper's own floor.

| Observation | Likely cause |
| --- | --- |
| uniform offset across all parameters | code bug (constant, sign, factor) |
| scales with a PHYSICAL parameter (1/d, 1/τ) | paper floor (omitted higher order) |
| scales with a NUMERICAL parameter (n, dps) | quadrature or precision bound |
| insensitive to every numerical knob, exact at one structural limit | fixed approximation level |

- **Always carry a moderate-parameter consistency check**: where the paper's approximation is
  tight the solver MUST agree to machine precision, or the verdict has no independent ground.
- **An upstream precision defect does NOT propagate at constant magnitude** (`[M]`
  1.2 % → 0.3 % → <0.1 % down one chain): run the sensitivity chain before claiming an upstream
  fix closes a downstream gap. Record: `atalay_r099_paper_floor_2026_05_03.md`.

## L6: A curvilinear matvec is verified only against a NON-FLAT per-ordinate reference

DUPLICATE: `AGENT.md` Step 5 carries this as a standing rule and `vv-principles` Mode 7 the
enumerate-the-exact-regime check. Forensics kept: TWO cylinder matvec bugs (ERR-049 bool-mask
scatter ordering; the decoder "analytical extension" O(h) twin path) hid for months behind
flat-ψ-only coverage — `cyl_matvec_twin_path_signatures.md` (#197/#206). The architectural cure
is `coding-elegance` Pattern 2: ONE `SNCellOperator` consumed by sweep and matvec.

## L7: Lift a direction-dependent per-ordinate moment to the GLOBAL frame before the angular reduction

ERR-061. A moment produced in the per-ordinate SWEEP frame and summed by `φ̂ = Σ_n w_n ψ̂_n` must
be sign-corrected to the GLOBAL frame for backward ordinates first — else forward and backward
slopes CANCEL and the scheme loses the diffusion limit. The fix is a per-octant moment-frame
involution at the producer–consumer seam.

- **The angular reduction is the discriminator**: a per-ordinate convention bug is invisible
  until a quantity is summed across ordinates of OPPOSITE sweep direction (M4).
- **When every component is individually correct but the fixed point is wrong, build a
  structurally-independent from-scratch kernel** — reproducing the wrong value bit-for-bit
  localises the bug to the SHARED math.
- **When the literature says the scheme IS consistent and your faithful implementation is not,
  the bug is a CONVENTION between two correct pieces.** Record:
  `issue_240_d5b_s3_diffusion_limit.md`.

## L10: "Error grows with refinement" + an unread library info-flag = an unconverged inner solve

DUPLICATE of `numerical-bug-signatures` Signature 8 (symptom, mechanism, three discriminators,
probe). Two sharpenings it does not carry: **derive `restart` from the operator's domain size**
(a hardcoded cap is the ERR-004 magic-constant anti-pattern, `coding-elegance` Pattern 7); and
**an L1 anchor can pass by NUMERICAL COINCIDENCE at one mesh size**, so **pair every analytical
anchor with a mesh-refinement leg**. Record: `krylov_restart_truncation_bug.md`.

**L10b — the flag need not be DISCARDED: recorded-but-unread gives the identical false
fingerprint, and the tell is TOLERANCE-INSENSITIVITY.** A red gate was bit-identical across
`inner_tol ∈ {1e-9 … 1e-15}`, which I read as a real discretization floor; the residual at the
cap was ABOVE every tested tol, so all four runs hit the same `max_inner`.

- ⭐ **A tolerance sweep discriminates only if at least one tol is LOOSER than the residual the
  capped run reaches**; otherwise it is a vacuous knob (M5 in convergence clothing). **Read
  `history.converged` / `n_inner` against `max_*` FIRST** — a plateau at exactly `max_iter − 1`
  is the whole diagnosis (`numerical-bug-signatures` Signature 8 carries this since
  2026-09-21).
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
  d=3 (→ L19).
- **ψ flat does NOT imply the TRACE is flat** — "homogeneous nulls redistribution" has a twin:
  *flat cell averages null the face mode*. The mechanism is L19·2 / L23 / L24.

## L11: For a ρ-honest stop or diagnostic, measure the residual `r = Aψ − q`, not `‖Δψ‖`

DUPLICATE: `AGENT.md` Step 3 carries the standing rule and `numerical-bug-signatures` Signature 9
the full symptom, mechanism and discriminator. **STALE API, landing named:** `FluxDisplacement`
was RETIRED 2026-08-19 (`coding-elegance` #18 records the reversal); the diagnostics now live on
the iteration record — `IterationRecord.increment_norms` / `.contraction_ratios` /
`.true_error_estimate()` in `orpheus/numerics/convergence.py`, with `AngularResidual` in
`orpheus/transport/residuals/angular_residual.py` carrying the per-ordinate balance map.
Catalogue: `issue_208_flux_displacement_residual_typing_debug_value.md`.

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
- **The MMS ladder is exactly blind to all of it** (M2): SI converges O(h²) over the whole ladder
  while the seed-iteration on a genuinely higher-order-in-μ field DIVERGES to NaN.

## L15: A SEED/closure re-pose is ruled PRINCIPLED-vs-REGRESSION by sweeping ANGULAR N, not `h`; and a grown Krylov composite must resize `restart`

`[LANDED a29ab2d]` (#282 route (a)). Cold:
`_archive/krylov_composite_restart_and_stale_references.md`.

- **(a) The N-sweep is the discriminator; the h-sweep at fixed N MISLEADS.** A seed IS an angular
  closure, so two treatments give different `keff(h→0)` at fixed angular order and the SAME
  `(h→0, N→∞)` answer — a gap that never closes under `h` reads FALSELY as a regression.
- **Honest nuance:** at the test's low N the OLD seed was *closer* to the N-converged truth. The
  justification was STRUCTURAL (an honest direct inverse, cold residual 5e5 → 1e-11), not
  accuracy. **Do not dress a structural win as an accuracy win.**
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
2. **GATE the driver bit-faithful to production** on NON-FLAT vacuum spheres (flat or
   homogeneous-reflective is degenerate, M4) BEFORE swapping anything — else a driver bug
   masquerades as the effect under study.
3. **Reference = fine-N in BOTH families plus a cross-family CONTAMINATION GUARD**, and report
   both error-vs-reference and the reference-free matched-N difference.
4. **MMS is blind** (M2): anchor with a closed form (a hand-derived `k_inf`) and `φ = Q/Σ_t`
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
  NONZERO block data in the gate (M4) — and a test that POSITIVELY PINS the blindness must INVERT
  to assert the flip reds.
- **The dense-probe recipe**: assemble `A`, `apply_transpose` and the production metrics as dense
  matrices by unit-vector probing `to_flat`, and check `T == Aᵀ` against **numpy's** transpose (a
  structurally-independent ground, not the operator's own machinery). `T == Aᵀ` ⟹ gauge-free up to
  SPD and reciprocity cannot choose it (the general theorem is L26·4); if `T` carries a weight, the
  metric is PINNED.

## L18: Adjudicate a LABELING/ORDERING degeneracy with the operator's own SYMMETRY GROUP

#326 `[CLOSED; remediation LANDED dde93b64 — a level's order is the FIBER's.]` Cold:
`_archive/curvilinear_seed_metric_and_ordering.md`.

- **An MMS whose ansatz AND source depend only on the INVARIANTS of the degenerate class is
  EXACTLY blind** — declare the ansatz's invariants BEFORE trusting it as an adjudicator; a
  companion ansatz outside the symmetric sector SEES the defect and still cannot ADJUDICATE.
- **A within-tie permutation cannot move a cumulative-sum coefficient**, so the whole effect is a
  LABELING and the defect magnitude is ordering-invariant.
- **The leak path is any coupling that maps ACROSS the degenerate class, and it need not COMMUTE
  with the tie-break — grep every cross-class index map before calling a relabeling inert.**
- ⭐ **The adjudicating criterion is a SYMMETRY the continuous AND semi-discrete problems both
  have**: no reference solver, no MMS, structurally independent by construction. Verdict here —
  **no ordering is correct, the CLOSURE is broken.**
- **The constructive exit is to fold to the FUNDAMENTAL DOMAIN**: a degeneracy in a sort key
  usually means the discretization carries a redundant symmetric copy, and on the independent half
  the key is monotone and every competing criterion coincides.
- Two promotion-time corrections: **`xfail` swallows FIXTURE SETUP ERRORS too**, so pair the xfail
  row with a **reddenable un-xfailed SIBLING** on the same fixtures (sharpening `vv-principles`
  Mode 8(4)); and a diagnostic that RE-IMPLEMENTS the production kernel must be rewired to CALL it
  at promotion, mutating through the test's OWN import binding (Mode 11).

---

## L19: An iterative-solver RATE question is a SPECTRUM question

#341 `[CLOSED; docstring repairs LANDED adc887d6; the octant-order lever is #343, OPEN]`. Cold:
`_archive/rate_spectrum_certificate_and_angular_axis.md`; record `issue_341_boundary_gs_rate.md`.

1. ⭐⭐ **Build the iteration matrix; do not re-time.** Any `x ← A_inv.apply(Σ gᵢ.apply(x))` driver
   IS a linear operator: wrap it over the composite's `to_flat`/`from_flat` as a
   `scipy.sparse.linalg.LinearOperator` and `eigs(which="LM")`. A few hundred sweeps buys the whole
   spectrum, and it is **immune to the stopping test**, so no ρ-blind stop (L11) or `max_inner`
   truncation (L10) can contaminate it. Controls (`instrument-doctrine` X1): reproduce the ρ FITTED
   from a real residual history, plus `G(2x) = 2G(x)` and `G(0) = 0`.
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

**Verdict discipline that generalises: a production default must never branch on a variable you
have only CORRELATED.** `ndim` was falsified on both sides by direct measurement, and thickness,
mesh, aspect ratio, quadrature order and `c` all move the sign at fixed `ndim`.

## L20: A RESIDUAL cannot gate an EIGENVALUE contract — measure the TRANSFER GAIN first

#340 N5, **REFUTED** on 38 solves / 8 geometries / 3 mixtures. Cold:
`_archive/rate_spectrum_certificate_and_angular_axis.md`.

1. ⭐⭐ **The one number that decides any "can statistic X gate contract Y?" question is the TRANSFER
   GAIN `|Δy| / X`, measured across configurations — compute it FIRST, before any threshold hunt.**
   (The `instrument-doctrine` skill's X1 carries it since 2026-09-21.)
   A threshold on `X` bounds `|Δy|` only through that gain, so an unbounded gain means no constant
   exists and the tuning is void (`[M]` gain spread 1.16e+05×; a zero-false-alarm threshold missed
   15 of 16 corrupting cases). This is `vv-principles` Mode 12 read in the MIRROR — the gate is
   SIGHTED on a class the CONTRACT is blind to. **Cure: project onto the functional the contract
   reads.**
2. ⭐ **A SIGNED projection against an APPROXIMATE weight is worse than no weight** — it manufactures
   accidental near-cancellations, i.e. false NEGATIVES. The cheap flat adjoint was itself *verified*,
   and correct for the WRONG PROBLEM. Pay for the real adjoint or use the unsigned norm.
3. **Answer the NULL case before the discrimination question, with a TWO-LEGGED tolerance sweep** —
   the outer leg moved the "pass" value 6 decades and the inner 0.2 %, so it was the OUTER's
   increment-stop slack (L11), not a floor; one leg alone mis-anchors the whole study.
4. **When lifting a production certificate one level, the CONSTANT does not come with it — and
   copying `record.binding_criterion.tolerance` silently picks the LOOSER criterion.** A residual bar
   scaled by an INCREMENT tolerance is a category error twice over.
5. ⭐ **Gate every verification fixture on the mixture's own consistency identity
   `σ_t == σ_c + σ_f + Σ_to SigS[0][g,:]` — an inconsistent mixture makes two legitimate references
   DISAGREE with no bug in either.** `[M]` the brief's "benign" reference was 30 % off because the
   fixture wrote `sig_s` as `[to, from]` while the constructor reads `[from, to]`, so the two solvers
   each honestly reported a different balance. **Never trust a brief's reference value on a
   hand-built mixture until the consistency identity is printed.** (`vv-testing` carries it since
   2026-09-21.)

## L21: An ANGULAR-consistency claim is separated by `h → 0`, not by the parameter it is named after

#319 / #235, 251 solves all `converged`. `[Phase-0 gates LANDED a3121cfe; #319 and #235 remain OPEN
for later phases.]` Cold: `_archive/rate_spectrum_certificate_and_angular_axis.md`; record
`issue_319_flux_dip_discriminator.md`.

1. ⭐⭐ **When a scheme claims consistency in the limit of ANOTHER discretisation, the axis that
   separates it from a rival is the OTHER mesh going to zero — because the claim is exact only
   there.** `[M]` sweeping optical thickness at fixed cells-per-mfp separated the two τ schemes **not
   at all** (fitted decay rate 0.000 for BOTH); refining `h` at fixed physics separated them without
   bound (3.2× → 204×). **Ask which limit the claim is exact in before choosing the sweep axis.**
2. **A regime sweep at fixed `c` self-destructs** — `Σ_a·R = Σ_t·R(1−c)` grows with it and every
   scheme agrees for a reason unrelated to the question; use the ε-scaling `Σ_t = 1/ε, Σ_a = ε,
   Q = ε`. (`vv-principles` #24(e) now carries the regime check; the half it does not is ⟹ **carry a
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
   not a bug.** Same family as L25·5 (M7).

## L22: A frozen reference is stale by MAGNITUDE CLASS, and the nulp count tells you neither

Triaging 9 "bit-identity" reds — all 9 stale references, none a regression, in two magnitude classes
the failure messages could not distinguish. Cold:
`_archive/krylov_composite_restart_and_stale_references.md`. (`numerical-bug-signatures`
Signature 10 carries the sibling-pass discriminator and the
re-baseline-against-an-independent-reference rule; `vv-principles` #25 carries "a re-baseline's
radius is the set of frozen REFERENCES by KIND — stored array, digest literal, in-test formula — not
the set of `.npz` files".)

1. ⭐⭐ **A nulp count is uninterpretable in BOTH directions — report `max|a−b| / max|b|` FIRST.**
   (`vv-principles` bit-identity since 2026-09-21, with point 3.) The
   received warning is "huge nulp near zero is nothing"; the dual bites as hard. `[M]` `1.04e+15`
   nulp meant the values differ by **8 %** (216/216 elements), while an `array_equal → False` on two
   arrays that PRINT identically was **1 ULP**. One measurement re-sorted the whole investigation.
2. **A gross move in the SCALAR flux refutes an ordering hypothesis on sight** — `φ = Σ_n w_n ψ_n` is
   permutation-invariant to `N × ULP`, so a percent-level move means the rule's VALUES moved.
3. ⭐⭐ **A rule-tier TOLERANCE pin can never warn about a consumer-tier BIT-IDENTITY pin — the gap is
   structural, not an oversight.** `gauss_legendre` is gated against numpy's `leggauss` at `< 8 ulp`
   (correct: neither construction is "the" answer), so a 3-ULP change is INSIDE that contract and
   OUTSIDE every downstream `array_equal`/`sha256` consumer simultaneously. **The only instrument
   that closes it is a byte-level fingerprint beside the tolerance pin**, whose sole job is "the
   bytes moved — re-baseline the consumers".
4. ⭐ **A gate whose reference is a different FP ASSOCIATION of production's own expression is a coin
   flip, not a contract** — `[M]` `src + (A + B)` against `(src + A) + B`: 68 of 80 slots
   bit-identical, 12 differ, max ULP 1, and the passing sibling passes by index luck.
   **Re-associate the reference or demote the assertion.**
5. ⚠ **Two probe-harness self-inflicted failures, both flattering** (M5): `grep -cE "^FAILED"` read
   `0` at nine commits including two already measured RED, because pytest's ANSI codes precede
   `FAILED` (use `--color=no`); and `python -c` prepends **CWD** to `sys.path[0]`, AHEAD of
   `PYTHONPATH`, so a worktree probe silently imported the MAIN tree and printed HEAD's values for
   every commit (run probe **script files** outside the repo, and make every probe print
   `module.__file__`; `code-search` carries both since 2026-09-21).

---

## L23: A discrete operator's SINGULARITY is a two-question object — refuse the either/or

#344 `[CLOSED; the closed-form kernel and the exit gauge LANDED f934ff57/b51bc802; the
characterization gates LANDED 1a2be025.]` Cold, with every measurement:
`_archive/issue_344_singularity_kernel_and_gauge.md`.

1. ⭐⭐ **When a claim is "the null space IS class X", the decisive number is `dim ker` MINUS `|X|` —
   and you must MEASURE `|X|`, not assume it exists.** One line (`min|Ω·n|`) refuted the brief's
   framing before any solve, because a level-symmetric rule cannot produce `Ω·n = 0` at all; and
   `dim ker = |X| + R` then held exactly on 9 rows, which turned an either/or into a DECOMPOSITION.
2. ⭐⭐ **A dense SVD through the PRODUCTION builders is cheap and settles rank questions that ARPACK
   only BOUNDS** (the prior record's "3 at d=2, ≥6 at d=3" was an `eigs(k=12)` lower bound; the truth
   was 12 and 138). **Report the singular-value GAP** so the rank threshold is visibly not arbitrary.
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
   (My own probe was wrong first, twice, both times toward the alarming verdict — M5.) That a
   splitting's fixed point is a MANIFOLD when `A` is singular, so a schedule-invariance gate
   legitimately reds with no bug present, is now `vv-principles` Mode 9's false-RED premise.
10. ⭐⭐ **`‖MM⁻¹−I‖` over the FULL space is the wrong instrument for an iteration's inverse — probe
    the RHS SUBSPACE the driver actually supplies.** `[M]` boundary-G-S reads 3.3e-01, which looks
    like incoherence and contradicts a measured `‖Aψ*−q‖/‖q‖ = 8e-14`; the defect sits exactly in the
    (inflow-row, outflow-column) block, where the driver supplies exactly zero content. **A reified
    forward-substitution "inverse" is a SUBSPACE inverse by construction** — fine for SI, a live
    hazard for a Krylov PRECONDITIONER (M6).

## L24: A discrete operator's KERNEL is usually a CLOSED-FORM problem

#344's basis, found in `0.05 s` where a dense SVD is `23 s` at half the size. Cold:
`_archive/issue_344_singularity_kernel_and_gauge.md`.

1. ⭐⭐ **A counting law that does NOT depend on a parameter the operator plainly contains is telling
   you the governing equation is COMBINATORIAL — go derive it, do not fit it.** Setting the degenerate
   branch (`ψ_c = 0`) turns DD's closure into an involution, so every face field is a sawtooth;
   substitute into the balance and **every** cross-section, mesh width, weight and area cancels,
   leaving *a sum of functions, each blind to one coordinate and one sign, vanishing identically*.
   Both fitted laws then drop out as theorems. **The parameter-independence WAS the derivation hint.**
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

1. ⭐⭐ **"Does accessor `f(x)` have to become `f(x, y)`?" is answered by asking what the defining
   conditions COMMUTE with — not by an expansion.** A **scalar convex combination commutes with every
   linear map**, so a quantity defined by membership of a convex set plus exactness of a scalar blend
   has both conditions as the SAME scalar statement in every component of every linear representation,
   and the widened form is an **overdetermined system whose every row returns the same value**. Three
   cheap hypotheses to check: the scalar is independent of the widened index; the projection is
   linear; the set is convex. Minutes, against a multi-day expansion.
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
   heuristic** (M7; same family as L21·6).
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

1. ⛔ **"Does X influence the adjoint?" is the WRONG question when the adjoint is a reverse-mode VJP** —
   `apply_transpose ≡ apply.T`, so the answer is "everything" and it carries no information. The
   informative question is the **symmetry character of the INCREMENT `∂A/∂k`**, measured by perturbing
   one entry and taking the dense difference. ⟹ **build the DIFFERENCE operator, not the operator.**
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

`_DENSE_METRIC_RCOND` re-derivation (ERR-080 / #429). A pinned `pinv` cutoff was justified by *"one
~1e-16 noise mode … 1e-12 sits ~4 orders above the noise floor"*. Every quoted number reproduced; the
reading was inverted. Cold: `_archive/rcond_threshold_rederivation.md`.

1. ⭐⭐ **Discriminate "round-off residue of an EXACT dependency" from "small real mode" by SOLVING for
   the dependency, then applying the candidate null vector to the RAW TABLE — not to the Gram.**
   Forming `G = AᵀWA` can itself manufacture rank loss, so `‖Gv‖ ≈ 0` is weak; `‖Av‖ ≈ 0` says the
   COLUMNS are dependent as functions on the node set and settles it. `[M]` the closed form fell out in
   two lines and matched the SVD null vector to 2.2e-16. ⟹ **a "noise floor" is a claim about a number
   that should have a NAME; if it has one, it is not noise.**
2. ⭐ **A threshold has TWO edges, and the instrument that justifies one is usually MONOTONE in it, hence
   blind to the other** (`vv-principles` #24(d), zero-set). `[M]` the Parseval statistic read a flat
   `1.000000000` across `[1e-15, 1e-2]` — a TRUE reading carrying **zero** information about the lower
   half. **Before citing a flat scan as justification, ask which direction the statistic can even move
   in.** (`vv-principles` #24(f) since 2026-09-21.)
3. ⭐ **Then go looking for the guard that already forecloses the other side — it is usually there and
   undocumented.** `[M]` every rcond below 8.696754e-17 is REFUSED at construction by a
   pair-consistency guard, so the corrupt band is **unreachable, not merely distant**. **Bisect the
   boundary; do not reason about it.**
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
   `automodule`, so `-W` cannot gate any of that prose at any severity. (My two instrument deaths that
   session — an unflushed docutils `warning_stream`, and a spliced-file import failing inside
   `dataclasses` because the module was not in `sys.modules` — are M5.)

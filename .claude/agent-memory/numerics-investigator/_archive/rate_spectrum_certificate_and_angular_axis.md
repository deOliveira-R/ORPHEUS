# Archive — the RATE/spectrum question, the outer certificate refutation, and the angular-axis discriminator

Moved COLD from `lessons.md` 2026-09-21 (verbatim: L19 lines 582-635, L20 lines 638-695, L21 lines 696-748).
Digest successors: L19 / L20 / L21 (imperatives only). Records: `issue_341_boundary_gs_rate.md`,
`scratch/issue_341_gs_jacobi_mechanism.md`, `scratch/n5_outer_certificate_measurement.md`,
`issue_319_flux_dip_discriminator.md`, `scratch/q68_flux_dip_discriminator.md` (the three scratch files are TRACKED).

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


# Archive — #344: the loss operator's singularity, its kernel in closed form, and the gauge

Moved COLD from `lessons.md` 2026-09-21 (verbatim, lines 803-1012 of the pre-distillation file).
Digest successor: `lessons.md` L23 / L24 (imperatives only). Records: `scratch/issue_344_null_space_structure.md`,
`scratch/issue_344_kernel_basis.md` (both UNTRACKED in the working tree — this archive copy is the durable one).

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
   `derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py` (10 green, 58 s). [gone: `git show a1c90aac^:derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py`]

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


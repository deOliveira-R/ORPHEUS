# Coupled block operator — the ψ½ Ray-Characteristic system realized in operator algebra

**Campaign plan — 2026-07-07.** Posing user-confirmed. **Supersedes** the FaceField carve's
remaining phases (the guard-retirement C3, the `StartingDirection→RadialCharacteristic`
rename C4, and the docs C5 all fold in here) and **subsumes** the coupled-system machinery.
Substrate: the **C2 FaceField carve LANDED** (`c637407` C2a / `4081c0d` C2b / `be5e7f8` C2c
on branch `refactor/sn-walk-unification`; the ψ½ field is `StartingDirectionField(FaceField
[tuple[int,int,str]])`, sibling of the spatial `BoundaryField`, `G_sd = V_cell` metric fix).

**Binding lens:** [[feedback-build-the-machinery-operator-algebra]] — the goal is the
MACHINERY (LEAN but LARGE), realize the operator-algebra formalism ALWAYS (the four
questions: *what is it in operator algebra / how is it posed / what invariant tests it / how
is it solved*). Risk/effort set the SEQUENCE (small steps so a bug stays findable), never the
WHETHER. The welded un-named `A_BA` is a defect in the machinery (a coupling block never
posed in the algebra), not an optional cleanup.

---

## The object (posing — confirmed 2026-07-07)

The within-group augmented SN system is a **2×2 coupled block operator** over two systems
(the composite biproduct re-partitioned — `Mat₂(Mat₂(𝒞)) ≅ Mat₄(𝒞)`; NOT a new type):

```
[ A_AA   A_AB ] [ transport ]      A_AA = L + C − S − B     (System A operator)
[ A_BA   A_BB ] [ ray       ]      A_BB = RayOp             (System B operator)
```

- **System A (transport):** state = `FullField` (`BulkField` bulk angular flux ⊕ spatial
  `BoundaryField` as its BC). `A_AA = L+C−S−B`.
- **System B (ray):** state = `RadialCharacteristicFlux` (a codim-1 `FaceField` — the ψ½
  cells at each radial cell) carrying **TWO boundary conditions**: r=R Dirichlet (the seed /
  given inflow data) **and** r=0 pole-reflection (a two-point BVP). `A_BB = RayOp` (the banded
  radial straight-characteristic march μ∂_r+σ_t, Hébert 3.434–3.435, + the two BCs).
- **`A_AB`** (ray→transport): the ray seeds the bulk angular recurrence — **posed as an operator**.
- **`A_BA`** (transport→ray): the bulk emission folded onto the ray source — the **un-named
  Schur elimination** currently hand-rolled into S/F — **posed as an operator** (the fold =
  named machinery, single-sourced via `fold_moments_to_starting_direction`).

**"Deciding where the block goes IS operator algebra":** once all four are posed as
operators, whether `A_BA` sits in the resolvent (folded), the explicit block operator, or a
DSA preconditioner is a *composition choice* the machinery supports. That flexibility is the
realization. The block-SOLVE stays spectral-character-keyed (the existing `InvertibleOperator`
/ `SourceIteration` / `MatrixInverseOperator` family — resolvent block-triangular-direct +
scattering iterative), which is already correctly realized.

## Specialist-verified facts (the measured algebra — 3 dispatches, 2026-07-07)

- **Within `(L+C)` the ray↔bulk coupling is block-triangular, solved DIRECTLY** (ρ=0, one
  pass): measured `A_sb = 0` EXACT (seed←bulk within L+C), `A_bs = 7.5` (`A_AB` M-M feed),
  `A_ss = 5.0` (`A_BB` march); `solve(apply(ψ))−ψ = 3.5e-16` (#284 forward-substitution).
- **Only the scattering coupling iterates:** `A_BA` lives ENTIRELY in the lagged gain S/F —
  measured `S_sb = 0.183` (bulk→ray source), `S_bs = 0` EXACT (ψ½ zero moment weight),
  `ρ(M⁻¹N) = 0.371 < c` (outer SI, Adams–Larsen). So the API reflects BOTH a block-triangular
  DIRECT solve (ray/bulk, in the resolvent) AND a block-ITERATIVE solve (scattering).
- **Extract is principled-equivalent, not bit-identical** (5.5e-16, different reduction tree);
  a NAIVE extraction returns O(1) garbage (the sweep treats inflow/seed rows as given data —
  the row-contract MUST be preserved). Principled-equivalence + the invariant test is the bar
  ([[feedback-principled-over-bit-identical]]).
- **The corrections that stand:** System B is a **two-point BVP** (r=R + r=0-pole → 2 BC data);
  the block-solve is **spectral-keyed** (reuse existing solvers); the ρ-honest **block residual**
  `r = Aψ − q` (via `evaluate_residual`, `solver.py:231`) is the convergence metric, NOT
  `‖Δψ‖`; **DSA** is a future CONSUMER (its R/P operators are a separate feature — not a
  risk-deferral); **RQI** (KKT/indefinite) and **multiphysics** (nonlinear/Newton) live in
  DIFFERENT machinery — out of scope.
- **Pre-carve baselines banked:** `derivations/diagnostics/diag_coupled_0{1,2}_*.py` (the
  block-triangular certificate + welded=exact-inverse) — promote to
  `tests/sn/operators/test_psi_half_coupling.py`.
- **Stale doc to fix:** the `"dormant until d3"` comments (`scattering.py:1562`,
  `fission.py:508`) — d3 LANDED @ `a29ab2d`; the A_BA arms are LIVE.

## The sequence — small de-risking steps, each a COMPLETE operator-algebra realization + its invariant

Each step answers the four questions and lands its invariant test (test-architect designs the
suite). Principled-equivalent, verified — not necessarily bit-identical.

1. **Pose System B.** `RadialCharacteristicFlux` (the `StartingDirection→RadialCharacteristic`
   rename lands here) + its two BCs (r=R Dirichlet, r=0 pole-reflection) + `RayOp = A_BB`.
   *Invariant:* the radial two-point BVP is well-posed; the straight-char march IS its solve
   (the r=R + r=0 data close it). *Corner→BC:* the corner is already separately-keyed;
   `_reflect_starting_direction` (`boundary.py:300`) is already its BC arm.
2. **Name `A_BA`.** Un-weld the Schur fold out of S (`scattering.py:1571`) and F
   (`fission.py:513`) into a named coupling operator (reusing `fold_moments_to_starting_
   direction`); fix the stale `"dormant"` comments. *Invariant:* it reproduces the S/F seed
   arms (single-source fold; `A_BA` = Fold ∘ bulk-emission). Low-risk (lagged gain OUTSIDE the
   resolvent — cannot touch #284).
3. **Name `A_AB`.** Pose the ray→bulk seed injection (`precompute_psi_state` reading
   `cells(p,−1)`) as an operator. *Invariant:* matches the in-sweep injection; the resolvent
   realizes it as forward-substitution; the adjoint `seed_cells_bar` is `A_AB.H`.
4. **Assemble `CoupledOperator [[A_AA,A_AB],[A_BA,A_BB]]`** + widen the role lattice
   (`SystemRole {A,B,COUPLED}` cloning the set-union `_join_block_roles`; System B first-class
   — it is a LIVE 3rd block absent from today's `BlockRole`). *Invariant:* block matvec ≡ the
   current augmented apply; block `.H` ≡ the current adjoint (the `(A.H).inverse()=(A.inverse
   ()).H` swap law extends). **Highest-risk step** if it touches the interleaved (L+C) walk —
   gate hard; may WRAP the sweep as the resolvent's block-triangular solve of the posed
   operators (naming, not re-implementing) as the first landing, extraction second.
5. **The solve.** Pose the resolvent (block-triangular DIRECT) + scattering (ITERATIVE) via
   the spectral-keyed solvers; the ρ-honest block residual as the stopping test. *Invariant:*
   coupled solve ≡ current solve (principled-equivalent, 5.5e-16 vs the dense-LU oracle).
6. **Retire the guards + `Optional` + the two-channel collapse.** Presence = the mesh-
   enumerated block-list / the operator's arity (structural, mismatch unconstructable); retire
   the 7 `_require/_refuse_starting_direction` guards + `FullField.Optional[starting_direction]`;
   retire the `sweep()` `starting_direction_source/_flux` kwargs (the seed rides only the
   composite — user ruling 2026-07-06). *Invariant:* a carrying-mesh composite missing the ray
   block is unconstructable; a non-carrying one carrying it is unconstructable.
7. **Verify + docs + close-out.** Promote `diag_coupled_0{1,2}`; the operator-algebra theory
   page (the 2×2 coupled system, the four blocks, the seed-as-BC); full SN+numerics+transport
   wall; pyright ≤ 1; sphinx -W. Comment #280/#282; the DSA seam (#2) named as a consumer.

## ⏸ Compaction points

Any step boundary is valid; recommended defaults ⏸ after **step 3** (System B + both couplings
named — the algebra is posed, the hard assembly remains) and ⏸ after **step 5** (the coupled
operator + solve landed — only retirement + docs remain). Before each: commit everything +
append a `## CHECKPOINT` block (step statuses done@hash / next, deviations, rulings) + tell the
user it is safe to `/compact`. Re-anchor from this plan + `git log` (trust git).

## Standing constraints

- **Surgical mode** — main agent writes directly, user steers; NO `method-implementer`;
  `test-architect` (proactive — this is the numerics↔transport↔sn triple-boundary operator-
  algebra carve), `explorer`, `qa`, `elegance-enforcer`, `archivist` encouraged.
- **Tests** — canonical `.venv/bin/python -O -m pytest -p no:xdist --timeout=300
  -p no:cacheprovider` SERIAL (xdist unstable for tests/sn+tests/numerics).
- **pyright** — oracle is `npx pyright` / the ratchet CLI (baseline `transport:1` = accepted
  #288); streamed `<new-diagnostics>` = documented #226 LSP artifact, IGNORE.
- **Pushes HELD**; `.claude/agent-memory/*` + `scratch/*` + `.claude/skills/*` NOT the main
  agent's to commit; never `git checkout`/`restore` on uncommitted files (monkeypatch to
  revert mutations); commits stamp **Opus 4.8**; sphinx -W must stay clean.

## Out of scope / deferred AS CONSUMERS (not risk-deferrals)

- **DSA (#2)** — SN⊗diffusion via R/P; a CONSUMER of this coupled-block machinery once its
  restriction/prolongation operators exist (a separate feature). The seam is named at step 7.
- **RQI** — KKT/indefinite saddle-point (`SaddlePointOperator`, trigger #294 mixed-form) —
  different machinery.
- **Multiphysics** — nonlinear (Newton/JFNK) — different machinery.

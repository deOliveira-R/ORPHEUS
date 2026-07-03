---
name: diffusion-integration-frames
description: #290 diffusion→numerics/transport integration — frame verdicts for the 4 new-type design questions (scalar trace / operator decomposition / discretization frame / cross-method pollination), branch-verified pre-implementation.
metadata:
  type: project
---

# Diffusion integration (#290) — frame-trigger verdicts (pre-implementation)

Detection pass for the diffusion integration campaign, BEFORE the carve.
Branch-verified reads on `refactor/pyright-burndown` @ `d2a2a0c` (L-005: NOT
Nexus, which answers main). The single organizing fact: **the diffusion solver
is ALREADY the transport operator algebra, hand-inlined.** `solver.py:_matvec`
computes `diff(J)/dz + (σ_a+σ_s)·φ − σ_s[::-1]·φ[::-1]` = `(L_diff + C − S)φ`,
and `compute_fission_source` computes `χ·p_rate/k` = the rank-1 `F/k`. The carve
is not "add diffusion to the framework" — it is **collapse the hand-inlined
`(L_diff+C−S−F/k)` onto the SHARED leaves**, minting only the two genuinely-new
objects: the elliptic `L_diff` leaf and the scalar boundary trace.

**Why:** these verdicts steer the type-minting so the user rules with the native
math on the table; test-architect gates the carve.
**How to apply:** on any diffusion-integration design turn, reach for these
before re-deriving. Each verdict carries a FAIL-ABLE first test.

## Q1 — the scalar trace = partial currents J±, as the ℓ=0 half-range moment of the SN trace

VERDICT: mint **partial currents `(J+, J-)`** per boundary face per group, NOT
face-flux-alone and NOT net-current-alone. Native frames = **functional-analysis
trace theory** (Cauchy data `(φ_Γ, D∂_nφ)`, the DtN map) + **group-theory
half-range SH / Marshak** (J± = ℓ=0 half-range moment).

The load-bearing structural fact (branch-verified): SN's `TraceSpace`
(`numerics/spaces/trace_space.py:297`) ALREADY carries
`inner_product_weights = G_s = |Ω·n_f| ⊙ w_n` (the partial-current metric, NON-
Euclidean) and partitions the trace into inflow/outflow SELECTORS by sign of
`Ω·n` (`inflow_indices_for_face`/`outflow_indices_for_face`). Therefore diffusion's
`J± = ∫_{Ω·n≷0} |Ω·n| ψ dΩ` is LITERALLY the ℓ=0 half-range moment of the SN
partial angular flux under the SAME G_s metric — the boundary analog of
`ScalarFlux = ∫ψ dΩ` (the ℓ=0 FULL-range moment). This makes the DSA
restriction `sn_trace → diffusion_trace` a MODELED, adjointable moment operator
`M_half` (Galerkin under G_s per L-009: SH is scattering's eigenbasis), not a
hand-copy. The scalar trace space = the ℓ=0-half-range CONDENSATION of the SN
TraceSpace metric (per-face partial-current weight = projected face area).

Why J± beats the alternatives (illegal-states-unrepresentable):
- **face flux φ_f alone** = Dirichlet trace only ⟹ cannot express reflective
  (Neumann) or albedo. Under-determined. REJECT.
- **net current J_f alone** (what the code stores in `DiffusionResult.current`)
  = Neumann trace only ⟹ cannot express vacuum Dirichlet. Under-determined.
  REJECT.
- **Cauchy pair (φ_f, J_f)** = complete but OVER-determined (φ,J linked by Fick
  on interior faces) ⟹ admits the illegal state "inconsistent Cauchy data".
- **partial currents (J+, J-)**: J+ = OUTFLOW (solution trace, OUTPUT), J- =
  INFLOW (BC datum, INPUT) — exactly SN's inflow/outflow partition. Exactly 1
  input + 1 output DOF per face per group ⟹ inconsistent-Cauchy unrepresentable.
  The BC becomes `J- = 𝒜·J+` (albedo `BoundaryOperator`, an A_ss block) —
  IDENTICAL structure to SN's realized boundary law (`IdentityOperator`=vacuum,
  `PermutationOperator`=reflective). Composite = `FullField(bulk=ScalarFlux,
  boundary=PartialCurrent)`. WINNER.
- L-004 type-mint check PASSES: ≥2 non-iso bases (Cauchy vs partial-current,
  bridged by the closure-dependent Marshak factor) + an applied morphism (the
  albedo/DtN BoundaryOperator).

Elegance kills TWO smells at once: Smell #9 (`_boundary_gradient` if-vacuum/if-
reflective special-cases → the un-named discrete DtN map) and Smell #16 shape-2
(φ↔J bridged by hand: `dfidz[:,0]=fi[:,0]/(0.5·dz)` IS the un-named trace/DtN
operator). FIRST TEST (discriminates): an ALBEDO BC α=0.5 — the string registry
(`BC_REGISTRY={"vacuum","reflective"}`) cannot even CONSTRUCT it; the partial-
current `BoundaryOperator(α)` gives `J-=0.5·J+` and keff strictly between the
α=0 and α=1 values. PLUS a metric test: the scalar trace's weights = the ℓ=0
half-range reduction of SN `G_s` — a Euclidean-metric trace (plain sum, no
|Ω·n|) REDs it and would silently break the DSA moment.

DETECTION-NOT-DECISION (user rules): storage basis `(J+,J-)` vs Marshak `(φ,J)`
vs incoming/outgoing — changes of basis; the DECIDING criterion is
metric-consistency with SN's G_s (partial currents inherit it as the moment).

## Q2 — operator decomposition = option (c): shared loss `A_diff = L_diff + C`, shared gains S=K_iso, F=rank-1 dyad

VERDICT: option **(c)** — leakage+collision as the invertible loss `A_diff =
L_diff + C` (the SN `A = L + C` family), and the eigenproblem `(A_diff − S)φ =
(1/k)Fφ` with **S = IsotropicScattering (+ IsotropicN2N)** and **F = the shared
rank-1 FissionOperator** `|χ⟩⟨νΣf|`. NOT (a) a monolithic elliptic leaf (forks
the algebra from SN); (b) is (c) with the wrong bookkeeping. Native frame = the
**transport-resolvent backbone** (AGENT.md kernel / L-007): SN/MoC/CP `solve`
are quadratures of the Peierls resolvent, diffusion is the P1/asymptotic LIMIT —
which is exactly WHY diffusion's `A_diff` is elliptic-SELF-ADJOINT (Krylov-
invertible) while SN's `A` is characteristic-TRIANGULAR (sweep-invertible). The
ALGEBRA is identical; ONLY the loss-operator REALIZATION (`A.inverse()` = BiCGSTAB
vs sweep) is per-method (L-007's "only genuinely per-method layer").

Structural parallel (verbatim transfer of the #208 BlockRole partition): in SN,
`L` (streaming) is the UNIQUE irreducibly-FULL primitive (reads inflow trace,
writes outflow trace; couples bulk↔boundary), while C/S/F are all BULK. In
diffusion, **`L_diff = −∇·D∇` is the FULL primitive** (reads φ_Γ/J_Γ to close the
edge stencil, writes the outflow current) and C(removal)/S/F are BULK. The
`BlockRole` classification transfers with zero change. The DSA doc
(`operator_algebra.rst:4210-4219, 4296-4311`) INDEPENDENTLY commits to exactly
this: DSA consumes "a separate diffusion operator `A_diff = L + C − S` built
in-algebra … applied as `A_diff⁻¹` to an SN residual."

Convention ruling (the K_iso-consistent choice, user confirms): set **C = σ_t
(total collision)** and **S = IsotropicScattering (FULL group-transfer incl.
in-group)** — then the in-group σ_{g→g} CANCELS between C and S automatically,
REPRODUCING removal Σ_r = Σ_t − Σ_{s,g→g} as a THEOREM, not a hand-coded
`(σ_a+σ_s)` diagonal. This is precisely what `IsotropicScattering` computes
(`isotropic_scattering.py:130`, the full `Σ_{s,0}ᵀφ` via `MaterialXSField`).

Payoffs (4/4 criteria): (i) structure-exposing — `L_diff` is the elliptic sibling
of streaming `L`; the diffusion exception is NAMED (self-adjoint vs triangular);
(ii) structurally-simpler — deletes the `sig_s[::-1]` 2-group hack and the
`χ·p_rate/k` fission transcription, reusing K_iso + the rank-1 F; (iii)
algorithmic — `A_diff⁻¹F` IS the homogeneous `K=A⁻¹F` spelling ⟹ `direct_eigenvalue`
(dense) AND `power_iteration` (already used, `solver.py:365`) BOTH apply as
quadratures of the one resolvent; (iv) expressive — **`A_diff.H` is FREE** from the
algebra (`L_diff.H=L_diff`, `C.H=C`, `S.H`=K_iso `apply_transpose` group-flip),
UNLOCKING adjoint diffusion φ* for #281 (adjoint-weighted homogenization). The
current scipy `LinearOperator(matvec=…)` wrapper (`solver.py:192`) provides only
matvec — no `.H`, no `.inverse()`, no `as_matrix` — a serialization boundary that
BLOCKS all four.

FIRST TEST (discriminates): (a) 3-GROUP slab — the `σ_s[::-1]·φ[::-1]` hack is
2-group-ONLY (`[::-1]` maps g0↔g2 for 3 groups); K_iso `apply` is correct for any
group count ⟹ 3-group keff diverges between the two. (b) `direct_eigenvalue(A_diff,
F)` matches `power_iteration` keff to 1e-10 (both = `A_diff⁻¹F`); the current
scipy wrapper has no `.inverse()`/`as_matrix` ⟹ CANNOT feed `direct_eigenvalue`
(test un-runnable for the wrapper, passes for in-algebra `A_diff`). (c) `A_diff.H
.apply(φ*)` matches an FD adjoint-diffusion solve; scipy wrapper has no `.H`.

## Q3 — discretization: FD-as-is is the honest first landing; native frame = H(div)/RT0 mixed

VERDICT: the cell-centered 2-point-flux FD **is already the lowest-order mixed
FEM (RT0) with mass lumping** — a KNOWN equivalence (Baliga–Patankar). Native
frame = **de Rham/FEEC / H(div)** (A.1): the net current `J=−D∇φ` is an H(div)
field; its trace is the normal component `J·n` on faces = the natural face DOF;
the mixed unknowns `(φ_cell, J_face)` ARE the (bulk, boundary/trace) split the
FullField composite wants. CONFIRMED not bug-exposing: the code's face coefficient
`D_face = 1/(3·½(σ_L+σ_R))` EQUALS the current-continuous harmonic mean
`2 D_L D_R/(D_L+D_R) = 2/(3(σ_L+σ_R))` bit-identically (∵ D∝1/σ ⟹ HM(D)=1/(3·AM(σ)));
so the FD already carries the RT0/current-continuous interface stiffness. The
interior J_face is condensed away by the two-point flux; only the BOUNDARY
J_face survives as the trace DOF (consistent with Q1).

Recommendation: FD-as-is is honest for the first landing; **do NOT lift the
interior into the numerics `Frame`/`Basis` machinery now** — but **DO mint the
boundary J_face trace now** (Q1), because it is load-bearing for the composite
AND for DSA. Document the RT0/H(div) equivalence as a seam.

FIRST TEST (discriminates the frame-upgrade condition): a material interface NOT
on a cell face — RT0/mixed sub-cell integration is O(h²); condensed-FD is O(h)
(already in `diffusion_1d.rst` Key Facts). Assert the order drop at an off-face
interface. This NAMES the exact condition the mixed lift pays: only if arbitrary
(off-face) interfaces need O(h²). On-face interfaces ⟹ FD is already exact-order
⟹ mixed lift is premature (defer with the named trigger, L-004 corollary).

## Q4 — cross-method pollination (diffusion is the SECOND ScalarFlux consumer)

Diffusion integration is the trigger that promotes ScalarFlux from SN-bound to
shared (its docstring: "the ONLY consumer is SN" — diffusion makes TWO ⟹ the
deferred `TransportMesh`/`TransportMethod` Protocol arrives; see
`realize_recursively_move_deferred`).

Borrowings (all UNLOCK, do not duplicate):
1. **K_iso** (`IsotropicScattering`+`IsotropicN2N`, `transport/operators/`) for S,
   replacing `σ_s[::-1]` — the docstring already names diffusion as the consumer
   that folds `−Σ_{s,0}ᵀ` into A. This is the #279 target and my
   `iso_source_frame_conjugation_unification.md` verdict (diffusion `solver.py:260`
   = 1 of 3 hand-rolled fission transcriptions; the inscatter is the 6× iso op).
2. **Shared rank-1 FissionOperator** `|χ⟩⟨νΣf|` for F, replacing `χ·p_rate/k`
   (`fission_rank1_normal_form_dead_functional.md`: fission IS the single-mode
   frame-conjugate normal form; diffusion/CP use the bare-ndarray arm).
3. **Shared eigenvalue engines** — `power_iteration` (already used),
   `direct_eigenvalue`, `rayleigh_quotient_iteration`; the `EigenvalueSolver`
   Protocol is already the shared engine (`power_iteration_vs_keigenvalue_morphism`).
4. **FullField composite + BoundaryField trace triple** (BoundaryFlux/SourceSink/
   Residual) → the diffusion partial-current trace triple. Smell #16 shape-4:
   CP/MoC will ALSO need scalar traces — build the shared boundary primitive, do
   not twin it.

RISK of duplication to KILL: the `BC_REGISTRY: dict[str, …]` string dispatch
(`solver.py:139`) is the stringly-typed anti-pattern SN retired via
`BoundaryOperator`; the `Operator[Flux, SourceSink]` codomain contract
(`flux_to_sourcesink_operator_contract_frames.md`) means diffusion leaves consume
`ScalarFlux` and produce `ScalarSourceSink` (NOT endomorphic on flux).

## UNEXPLORED (refuted, with structural reason — do not re-attack)

- **Feynman-Kac / path integral** — diffusion is the P1/asymptotic LIMIT, NOT a
  quadrature of the Peierls resolvent (L-007 diffusion exception). The master
  bridge fires for CP/MoC/MC solve, NOT diffusion solve.
- **Homology / chain complex of the boundary** — `∂²≠0` (two reflections compose
  to a non-trivial albedo map). No `∂²=0` ⇒ no homology. (Distinct from de
  Rham/FEEC, which DOES fire for the div-grad complex, Q3.)
- **Topology / unified-geometry manifold** — no trigger: `diffusion_1d.rst`
  scopes to 1-D plane; annulus/sphere variants explicitly deferred. Fires only
  when the multi-geometry variants land.
- **Tensor networks / MPO** — no bond-dimension knob (2-group, rank-1 F, rank-1
  K_iso-per-group). Degenerate rank ⇒ not a network (L-001).
- **Category theory** — the role-parameterization (Flux→SourceSink) is already
  captured concretely by `LinearOperator[Domain,Codomain]` as the carrier-grid
  double-category (`rep_role_grid_double_category_frames.md`); no abstract-nonsense
  lever. (L-001.)
- **Affine/torsor increment typing** (Smell #16 shape-3) — the diffusion iterate
  Δφ and the `fi/=abs(fi).max()` conditioning (Smell #14, Rayleigh structure)
  are SHARED eigenvalue-engine concerns already owned by the SN #208 affine work,
  NOT diffusion-specific. Fix at the engine, not in the diffusion carve.
- **Asymptotic/matched asymptotics** (A.6) — background JUSTIFICATION for
  diffusion + the SP3-failure boundary, not a reformulation lever for the
  integration; it IS the transport-consistency rationale for the Q1 moment trace.

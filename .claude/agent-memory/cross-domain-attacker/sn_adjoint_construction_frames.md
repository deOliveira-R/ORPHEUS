---
name: sn-adjoint-construction-frames
description: Frame attack on the SN adjoint transport solve (P6 prerequisite — adjoint-weighted homogenization). VERDICT — the SN adjoint is NOT a new layer; it is a DAGGERED POSING ROW over the SAME resolvent backbone (L23b confirmed in code). Construction route RANKED — build the EXACT DISCRETE transpose A_loss† (S†, F†, the already-shipped (L+C)† + B†), NOT a μ-reversal-reuse-forward-sweep, BECAUSE adjoint-keff=forward-keff is exact ONLY for the discrete transpose, and the discrete-transpose reverse-walk Lᵀ is ALREADY in the code (loss_action_transpose). F† IS literally a free dyad-swap (outer(νΣf-recon, χ-functional) — the RankOneOperator docstring already promises it). S† IS free by Funk-Hecke self-adjointness (Λ real-diagonal ⟹ S=R∘Λ∘M is Euclidean-symmetric up to the SH frame's two metrics). Cross-method seam: the daggered-triple posing is method-agnostic (lives at numerics/transport); the per-method work is ONLY the leaf transpose (already done for SN's L; CP=matrix.T, diffusion=self-adjoint).
metadata:
  type: project
---

# SN adjoint construction — foreign-frame attack (P6 prerequisite)

Grounded on main @ `68ceb9a`, branch-verified by direct Read (L-005). The brief's
five structural questions resolve to ONE native frame + a ranked construction route.

## The native frame (three confirmed, all already-named in the campaign)

1. **Daggered posing row over a shared resolvent backbone** (L-007 + eigenvalue-posing
   L23b — CONFIRMED in code). `eigenvalue.py:27-29` already documents the adjoint row
   `A_loss†ψ†=λM†ψ†` as a "future seam" at the SAME posing layer as K and α. The
   adjoint is NOT a new layer: `power_iteration` (the `fix(step)` combinator) never
   learns forward-vs-adjoint; only the role-operators dagger. The forward k-eigenvalue
   spectrum is REAL (Krein-Rutman, `eigenvalue.py:10-14`) ⟹ forward & adjoint share
   the dominant eigenvalue ⟹ λ unchanged, only the operators transpose. This is the
   `(L+C,S,B,F) → ((L+C)†,S†,B†,F†)` posing the brief asks about in Q3 — and it is
   the ONLY new structure needed.

2. **Lewis-Miller backward semigroup (`Ω→−Ω` = path reversal)** — the transport-
   resolvent backbone's adjoint face (AGENT.md kernel). The continuous SN adjoint is
   the forward operator with direction reversed + kernels transposed. BUT the code's
   `apply_transpose` does NOT realise this by μ-reversal-forward-sweep; it realises it
   by an explicit **discrete reverse-substitution walk** (`loss_representation.py:2602
   loss_action_transpose`: `for i in reversed(cells)` + boundary block swap
   outflow-cotangent→inflow-cotangent + `closure.angular_adjoint` for the nested
   angular factor). The streaming-apply-transpose frame's THREE gotchas
   ([[streaming_apply_transpose_frame]]) are ALL already handled in this body.

3. **Dagger inverse biproduct category** ([[issue_208_operator_algebra_frames]]) — the
   `†=G⁻¹AᵀG` metric-adjoint with the `|Ω·n|·w` boundary metric. `.H` (`operator.py:543
   adjoint` / `_AdjointOperator`) ALREADY wraps `apply_transpose` in the G-conjugation.
   `(L+C).H` and `(L+C−B).H` are ALREADY reachable (both transpose-closed); only
   `S.H`/`F.H` are missing because S/F advertise `{CAP_APPLY}` only — by DESIGN
   (`streaming.py:360-369`: the capability lattice makes `A_loss.H` deliberately
   unreachable so it fails LOUD, never silently Euclidean).

## The construction route — RANKED (the deliverable)

**WINNER: build the EXACT DISCRETE transpose `A_loss†` (add `S†`, `F†`; the
`(L+C)†` + `B†` already ship). NOT μ-reversal-reuse-forward-sweep. NOT
Krylov-on-apply_transpose for the inner.**

Why discrete-transpose over μ-reversal (Q1+Q2 decision):
- The adjoint-keff = forward-keff identity is **EXACT only for the discrete adjoint**
  (`A†` of the DISCRETIZED `A`), because the power method's Rayleigh quotient closes
  exactly only when `⟨φ†, Aφ⟩` uses the true matrix transpose. "transpose-then-
  discretize" (μ-reversal physical adjoint) and "discretize-then-transpose" (the code's
  reverse walk) DIFFER at: the WDD diamond weight (the `1/τ` closure is not μ-symmetric
  on curvilinear), the boundary closure, and the curvilinear angular-redistribution
  `(1−µ²)/r ∂_µ` (which the discrete reverse walk transposes via `closure.angular_adjoint`,
  a SECOND triangular factor — a μ-reversal of the forward sweep would leave it
  UN-transposed = silently wrong, gotcha #1).
- The discrete transpose is **already implemented** for the hard leaf (`(L+C)†` via
  `loss_action_transpose`, pinned by `test_g_adjoint_reciprocity` slab/sphere/cylinder
  + its L11 wrong-trace-metric negative control). μ-reversal would be a SECOND,
  divergent path to the same operator (Smell #16 shape-4 — a third hand-rolled path).

Why NOT Krylov-on-apply_transpose for the INNER resolvent: the adjoint inner solve
`(L+C)†⁻¹` should be the **reverse-direction WDD sweep** (back-substitution on the
reversed DAG), the adjoint twin of the forward sweep, NOT GMRES on `(L+C).H.apply`.
The reverse sweep is the analytic inverse of `(L+C)†` exactly as the forward sweep is
of `(L+C)` (transport-resolvent backbone: adjoint solve = backward semigroup). BUT —
the code currently ships `(L+C).apply_transpose` (the matvec) but NOT
`(L+C).solve_transpose` (the reverse-sweep SOLVE): `operator.py:650-651` defers
`(A.solve).H = (A.H).solve`. **This is the ONE genuine new sweep primitive the adjoint
needs** — a `solve_transpose` = reverse-DAG forward-substitution, the inverse twin of
the reverse-walk matvec already present. Everything else is dagger-functor algebra.

## F† IS a free dyad-swap — CONFIRMED (Q4 verdict)

`RankOneOperator` docstring (`operator.py:1752-1756`) ALREADY promises it verbatim:
"`apply_transpose` is the dual dyad `|w⟩⟨v|` (swap the column and the row) ... when
the adjoint transport machinery lands it falls out by swapping the two factors."
`F = outer(χ, ReactionRateFunctional(νΣf))`; `F† = outer(νΣf_as_vector,
InnerProductFunctional(χ, axis=0))`. The emission↔production role swap (χ↔νΣf) IS the
column↔row swap of the dyad. NO new operator code — it is the SAME `outer` primitive
with arguments swapped. ONE subtlety: `χ` becomes the FUNCTIONAL weight (a co-vector)
and `νΣf` becomes the reconstruction COLUMN — both are already per-cell `(ng,*spatial)`
arrays, so the swap is a constructor re-wiring, not new math. (`χ` as a functional is a
generic `InnerProductFunctional`, NOT a `ReactionRateFunctional` — χ is a spectrum, not
a cross-section; the typed `ReactionRateFunctional` carries `νΣf` on the FORWARD side.)

## S† IS free by Funk-Hecke self-adjointness (Q4 verdict, the deep one)

L-009 (eigenbasis ownership): `S = R∘Λ∘M` is the SPECTRAL THEOREM `UΣU*` with Λ the
REAL-diagonal Legendre moments (`scattering.py:510 frame`). A real-symmetric-spectrum
operator's transpose is `Sᵀ = Mᵀ∘Λ∘Rᵀ` — and since the SH frame's analysis/reconstruction
are an adjoint PAIR (the `GalerkinFrame` M*=R up to the two metrics, the landed P1.3
`R.H` work), `S†` is structurally `S` with M↔R swapped, which for a self-adjoint zonal
kernel collapses back toward `S` itself (up to the `(2ℓ+1)`/`4π` Gram factors that the
frame already carries). The energy transfer `Σ_{s,ℓ}(g'→g)` transposes to `(g→g')` —
the ONLY genuinely asymmetric piece, a `Λ` matrix transpose per ℓ. So `S†` = the same
`R∘Λᵀ∘M` with Λ's group-transfer transposed: a transposed `Mixture`/kernel input
(brief's Q4 hypothesis CONFIRMED — "no new operator code, transposed kernel data"),
NOT a new `apply_transpose` arm. The angular part is self-adjoint (Funk-Hecke); only
the group-transfer matrix transposes.

## Cross-method seam (Q5) — the daggered-triple posing is method-agnostic

The `(L+C,S,B,F)→dagger` posing row lives at the SAME layer as K/α (numerics/transport
posing, method-agnostic per L23b/eigenvalue-posing). The per-method work is ONLY the
leaf transpose: SN's `(L+C)†` (reverse sweep — DONE for matvec, needs `solve_transpose`),
CP's `A†=Aᵀ` (dense matrix transpose — trivial), diffusion's `A†=A` (self-adjoint,
elliptic — `A†=A` exactly, the transport-resolvent EXCEPTION L-007). So a
method-agnostic `AdjointPosing` that daggers the triple + a per-method
`leaf.apply_transpose`/`solve_transpose` is the right factoring. TRIGGER for building
the shared abstraction: the 2nd method's adjoint (CP or diffusion) — until then SN's
adjoint is built directly, but the posing row is written method-agnostically so #2
drops in (unify-after-two; the adjoint row is ALREADY anticipated in `eigenvalue.py`).

## First tests that DISCRIMINATE (L-002)

1. **Discrete-vs-continuous adjoint keff (the route decider).** Build `A_loss†` two
   ways: (a) the discrete reverse-walk transpose; (b) a μ-reversal forward sweep on the
   antipodal-permuted ordinate set. Run adjoint power iteration both ways on a
   curvilinear (sphere N≥4) heterogeneous problem. Assert `k_adjoint_discrete ==
   k_forward` to ~1e-12 (EXACT — same matrix, transposed) but `k_adjoint_mu_reversal`
   DIFFERS at the WDD-closure-bias level (~1e-3 on curvilinear). The μ-reversal route
   PASSES on a uniform slab (where the closure IS μ-symmetric) and FAILS on curvilinear
   — that divergence is the discriminator. A test that only ran on a slab could not
   tell the routes apart (cannot-fail = rejected).
2. **F† dyad-swap bit-identity.** Build a dense `F` matrix, take `F.T`. Assert
   `outer(νΣf, InnerProductFunctional(χ, axis=0)).apply(y) == F.T @ y` via
   `array_equal` (0 ULP — same `(w*x).sum` primitive, swapped operands). A "new F†
   operator that recomputes" would only `allclose`, not `array_equal`.
3. **`(L+C).H.H == (L+C)` involution** (dagger axiom f††=f) — free the moment
   `solve_transpose` lands. A wrong boundary-block swap in the reverse sweep fails the
   round-trip.
4. **G-metric load-bearing on S† (negative control).** Build `S†` with the WRONG
   angular metric (Euclidean, dropping the `(2ℓ+1)/4π` Plancherel weight). The
   Euclidean reciprocity `⟨Sφ,φ†⟩=⟨φ,S†φ†⟩` PASSES; the G-metric one FAILS. Proves the
   SH frame's metric is load-bearing for the adjoint (the [[issue_208_operator_algebra_frames]]
   metric-blindness smell, now on the adjoint path).

## REFUTED / low-signal for THIS question

- **μ-reversal-reuse-the-forward-sweep as the PRIMARY route** — refuted: it gives the
  CONTINUOUS (physical) adjoint, which makes adjoint-keff ≠ forward-keff on a
  discretized problem (the identity is exact only for the discrete transpose), AND it
  is a second divergent path to `Lᵀ` (Smell #16 shape-4) when the discrete reverse walk
  already exists. It is a legitimate VERIFICATION ORACLE (continuous-limit cross-check),
  not the production path — the coding-standards "fuller-view oracle" exception applies
  IF kept, but only as a uniform-mesh principled-equivalence check.
- **Krylov-on-apply_transpose for the inner adjoint solve** — refuted as PRIMARY: the
  reverse WDD sweep IS the analytic inverse of `(L+C)†` (backward semigroup); GMRES on
  the transpose matvec is the fallback for when no reverse-sweep exists (CP/diffusion
  dense), not for SN where the triangular structure gives the sweep for free.
- **A new `AdjointOperator` per-leaf class hierarchy** — refuted: `.H`/`_AdjointOperator`
  + the dagger functor (`(A+B).H=A.H+B.H`) already give the adjoint algebra for free;
  the leaves need `apply_transpose`/`solve_transpose` METHODS, not adjoint CLASSES
  (L-004 type-theatrics: the adjoint is a wrapper, the metric is on the space, the
  transpose is a method — no new type earns its place).

---
name: phase4-o2b-4-6-mms-ansatz-frame
description: Durable frame for designing an anisotropic SN MMS ansatz — (A(r)+μB(r))/W IS the native truncated-P1 (P0+P1) Legendre element; a LINEAR angular closure is FULLY probed by an input non-constant in μ, so linear-in-μ fully (not partially) activates the curvilinear redistribution (1−μ²)/r ∂ψ/∂μ — no P2 needed. Pole-regularity hazard: the redistribution carries a 1/r, so a curvilinear B MUST vanish at the centre (B(0)=0); a slab B has no such constraint.
metadata:
  type: project
---

# Anisotropic SN MMS ansatz (A+μB)/W — native-frame + hazard analysis

The reusable frame content from the Phase 4 / O.2b 4.6 non-vacuum-BC MMS design. The
case-specific term tables landed; what survives is the cross-domain detection.

## Frame: (A+μB)/W IS the native P1 Legendre element (not an ad-hoc ansatz)
`(A(x)+μB(x))/W` is the truncated P0+P1 Legendre form: `P0=1`, `P1(μ)=μ`, so `A` is the
P0-moment shape and `μB` the P1-moment shape. The SN closure's native angular basis (the
Carlson seed folds the moment-source via `P_ℓ(−1)=(−1)^ℓ`) IS the Legendre-moment basis.
Linear-in-μ = P1-truncated = the LOWEST non-trivial moment with a non-zero `∂_μ`. There
is NO frame mismatch — linear-in-μ is the native P1 element. DETECTION: when an MMS
ansatz is proposed in a discrete-ordinates angular variable, check whether its
μ-dependence is a Legendre-moment truncation; if it is, it is native and the question
"do I need a richer angular form?" reduces to "which operator terms does each moment
activate?".

## Frame: a LINEAR operator is fully probed by a non-constant input (the key payoff)
The curvilinear redistribution operator `(1−μ²)/r · ∂ψ/∂μ` is LINEAR in ψ (the discrete
closure — M-M half-angle recurrence + Carlson seed — is `is_linear=True`). With ψ linear
in μ, `∂ψ/∂μ = B(r)/W` is a non-zero CONSTANT in μ; multiplied by the μ-dependent dome
`(1−μ²)` it produces a genuinely μ²-structured redistribution term `(1−μ²)B/r`. A LINEAR
operator is FULLY probed by an input that is non-constant in its argument — so linear-in-μ
exercises the FULL redistribution map (half-angle recurrence + μ²-second-moment coupling),
not a partial slice. A quadratic-in-μ (P2) term adds NO new structural coverage of a
linear operator — it only changes WHICH point in the operator's already-fully-probed range
you land on. DETECTION HEURISTIC (general): to fully exercise a linear operator with an
MMS, you need an input non-constant in the operator's active variable — NOT a
high-polynomial-degree input. Enrich the ansatz degree ONLY to test a QUADRATURE exactness
(e.g. `Σw_nμ_n²`), never to "more fully" probe a linear map.

## HAZARD: pole-regularity of the redistribution term (1/r) — curvilinear-only
On the sphere the redistribution is `(1−μ²)B(r)/r`. As `r→0` with `B(0)≠0` this is
`(1−μ²)B(0)/r → ∞` — a 1/r singularity in the manufactured source. So on a CURVILINEAR
geometry `B(0)=0` is a HARD regularity constraint (only the OUTER-boundary value `B(R)` is
free). The slab has NO such constraint (no 1/r). Concrete fix: an `(r/R)`-prefactored B,
`B(r)=(r/R)·[b0+b1·cos(k_r r)]`, gives `B(0)=0` (pole-regular) while `B(R)≠0` (lights the
prescribed inflow). The P0 shape A has no 1/r companion, so `A(0)≠0` is fine. GENERAL
DETECTION: any MMS that lights an angular-redistribution / curvature term must check the
measure-induced singularity of that term at the coordinate origin — a slab-derived ansatz
silently drops BOTH the redistribution term itself (Cartesian has no `∂_μ` operator) AND
its pole-regularity constraint. This is the single largest correctness hazard when lifting
a Cartesian MMS to curvilinear.

## The measure also enters the ERROR NORM, not just the source
The volume/area measure (`4πr²dr`) does NOT enter the manufactured continuous source (it
is per-unit-ψ; the discrete consumer applies the measure via the cell-balance V_i). But
the L2 convergence norm MUST be volume-weighted — a Cartesian (unweighted) norm mis-measures
the order of accuracy. Check the convergence norm carries the geometry's measure.

## Architectural note (Cardinal Rule 2)
A boundary-non-vanishing variant of an existing proven ansatz is a BOUNDARY-TRACE
PARAMETERIZATION of that case, NOT a new dataclass: generalize A,B to admit
non-vanishing-at-the-outer-face coefficients, wire the non-vacuum prescribed inflow, keep
the curvilinear `B(0)=0` invariant. Reusing the SymPy Branch-1 identity machinery with the
generalized A,B re-proves the source for free.

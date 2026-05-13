---
name: morel-1989-si-vs-apply-equivalence
description: Morel 1989 (NSE 101) — verbatim extraction of his explicit statement on when the standard SN SI sweep is and is not equivalent to the discrete-ordinates linear system, with the lower-triangular / full-matrix distinction that maps directly onto ORPHEUS curvilinear M-M angular coupling.
metadata:
  type: project
---

# Morel 1989 — SI sweep vs. discrete-ordinates linear system

## One-line answer

**Yes, Morel addresses this — explicitly and in plain text.** His
exact criterion is: the standard SN iterative SI sweep is
*equivalent to the linear system* (i.e., reaches the system's
solution) **if and only if** the angular leakage matrix is
**lower triangular**. For a **full** angular leakage matrix the
SI sweep "cannot be used" as-is; it can only be used by
"placing the off-diagonal components on the right-hand side of
the equation and iterat[ing] on" them — i.e., as an outer
Jacobi-like fixed-point loop, not a direct sweep — and that
fixed-point loop is conditionally stable, with spectral radius
approaching unity in "certain limits."

This is the canonical literature statement of the structural
non-equivalence the numerics-investigator hypothesised at
`.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_replan_verdict.md`.

## The verbatim load-bearing quote

Source: Morel, J. E. (1989). "A Hybrid Collocation-Galerkin-SN
Method for Solving the Boltzmann Transport Equation," *Nuclear
Science and Engineering* **101**, 72–87. DOI:
[10.13182/NSE89-4](https://doi.org/10.13182/NSE89-4).
Local PDF:
`scratch/literature/Morel(1989)A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann Transport Equation.pdf`.

**Section III, page 75**, in the paragraph immediately following
the display of the collocation, Galerkin and SN angular-leakage
matrices `L_C`, `L_G`, `L_{S2}` (Eqs. 12, 13, 14 of the paper):

> "The spatial leakage matrices for the collocation and *S_n*
> methods are identical. This is always so and is a direct
> result of making the interpolation and collocation cosines
> identical. **The Galerkin matrix is full rather than diagonal.**
> Although the Galerkin approximation can produce a diagonal
> matrix, it very rarely does so. **If the matrix is full, the
> standard *S_n* iteration sweep cannot be used unless the
> off-diagonal components are placed on the right side of the
> equation and iterated on.** We have performed a simple Fourier
> analysis that indicates that this iteration process is both
> stable and rapidly convergent for this particular example
> case, but for the general case it is not clear that the
> associated increase in cost will be justified by an increase
> in accuracy. Therefore, we conclude that *S_n*-collocation
> approximation should be used for this term." (emphasis added)

And then in the immediately following paragraph (still p. 75),
about the **angular** leakage matrix:

> "The angular leakage matrix for the collocation method is not
> valid because the collocation method cannot be applied with a
> discontinuous trial space. Even if it is applied with a valid
> trial space, the collocation method often yields an angular
> leakage matrix that is not conservative. For this example
> case, the angular leakage matrix is conservative if and only
> if the elements in each column sum to zero. Note that both
> the *S_n* and Galerkin matrices have this property. **The
> Galerkin matrix is full whereas the *S_n* matrix is lower
> triangular. The standard *S_n* iteration sweep cannot be used
> with a full matrix unless the upper triangular components are
> placed on the right side of the equation and iterated on.** A
> simple Fourier analysis shows that this iteration process is
> stable for the example Galerkin matrix, but the spectral
> radius of the iteration matrix approaches unity in certain
> limits. **This means that arbitrarily slow convergence rates
> could be obtained.** Because there are significant
> disadvantages associated with both the collocation and
> Galerkin methods, we conclude that the standard *S_n* method
> should be used for the angular leakage term." (emphasis added)

## How this reads onto ORPHEUS

The numerics-investigator's hypothesis at
`issue_196_phase_g_step2_replan_verdict.md` is now textbook-confirmed:

| Morel's distinction                  | ORPHEUS analog                              |
| ------------------------------------ | ------------------------------------------- |
| `L` = angular leakage matrix         | per-cell angular coupling block (the M-M `α_{m±1/2}` ladder) |
| `L` lower-triangular                 | causal angular coupling (one direction's outgoing flux feeds the next, like the inner-to-outer μ-sweep order in `orpheus.sn.solver.sweep`) |
| `L` full                             | non-causal coupling — the curvilinear M-M cross-derivative term effectively couples *all* discrete μ values within a cell because the `α_{m±1/2}` recurrence produces a closed loop when reflective BC closes the μ-domain |
| "standard SI sweep cannot be used unless off-diagonals on RHS and iterated on" | what the SI sweep in `orpheus.sn.solver` ACTUALLY does: it advances one μ at a time, treating the other μ's coupling terms as a previous-iterate source — i.e., it's the Jacobi-on-the-full-matrix that Morel describes, not a direct solve |
| "stable but spectral radius approaches unity in certain limits → arbitrarily slow convergence rates" | the 681 outer iterations to 9.55e-15 residual reported in the verdict — *fixed point reached* but *not the system's solution* because… |
| (the implicit corollary)             | …the converged SI fixed point satisfies the equation `L_lower ψ = RHS − L_upper ψ_prev`, which only equals the full-system `(L_lower + L_upper) ψ = source` IF the two iterates `ψ` and `ψ_prev` coincide — true at convergence in scalar theory, but the M-M recurrence within a cell breaks this because the within-cell `α_{m+1/2}` and `α_{m-1/2}` couple to the **same iterate**, while across iterations the residual coupling is asymmetric. The result: SI converges, but to a different fixed point than the full discrete-ordinates linear system. |

The 22% L0 discrepancy between SI sweep and Krylov apply-matvec
on the homogeneous reflective sphere is exactly the symptom of
this asymmetry: both methods are stable, both converge to
machine precision, but they are solving structurally different
fixed-point problems on a curvilinear M-M angular coupling
matrix that is **full** (not lower-triangular).

## What Morel proposes as the remedy

Morel's specific remedy in this paper is **not** to fix the SI
sweep — he explicitly recommends against the Galerkin / full-matrix
approach for the angular leakage term and falls back to standard
SN for that term (last sentence of the page-75 passage). His
"hybrid collocation-Galerkin-SN" of the title applies the
Galerkin treatment ONLY to scattering and the spatial leakage,
**keeping the angular leakage discretization standard SN
(lower-triangular)**. From p. 78, Sec. IV opener:

> "Our hybrid method differs from previous hybrid methods in
> that a combination of finite element and *S_n* treatments is
> applied to the angular variables." […] "Based on our
> comparison of the collocation, *S_n*, and Galerkin methods,
> we define our hybrid method as follows. The collocation method
> is formally used for the spatial leakage term, the *S_n*
> method is used for the angular leakage term, and the Galerkin
> method is formally used for the removal term and for the
> inscatter term."

So Morel's response to "full angular leakage matrix breaks the
SI sweep" is to *avoid the full angular leakage matrix entirely*
by sticking with the standard SN (lower-triangular) form. He
sidesteps rather than solves the equivalence question.

## What Morel does NOT address

- He does **not** quantify the converged-SI vs. full-system
  discrepancy (no "22% L0 error" type number). His Fourier
  analysis is on the **iteration matrix's spectral radius**
  (convergence speed), not on the **fixed-point displacement**
  (correctness of the converged fixed point).
- He does **not** treat the M-M curvilinear redistribution as
  the example. His pages 73–76 use a spherical-geometry
  example with a half-range double-Gauss S_2 (μ_1 = −0.5,
  μ_2 = +0.5) trial space — comparable but not the n=40 GL-8
  ladder.
- He does **not** discuss reflective BC closure (he uses a
  generic spherical example without specifying the boundary
  closure). The reflective-BC closure of the μ-domain in the
  ORPHEUS test is what makes the upper-triangular residual
  load-bearing — for vacuum BC the upper-triangular part is
  truncated and the SI/system discrepancy would be much smaller.

## Caveats

- Morel's example (p. 74, Eqs. 9 + 9a/b) uses a **half-range
  step trial space** (`B_1(μ) = 1 for μ<0`, `B_2(μ) = 1 for μ≥0`).
  This is a much coarser angular discretization than ORPHEUS's
  n=40 GL-8. The full-vs-lower-triangular distinction he draws
  is *qualitative* — it tells you which matrix structure SI
  can handle, not how big the discrepancy will be for any given
  quadrature.
- The recommendation "place upper-triangular components on the
  RHS and iterate" is exactly the SI-with-old-iterate scheme.
  Morel calls it "stable" for his example case but flags
  "spectral radius approaches unity in certain limits." His
  qualifier "the general case it is not clear that the
  associated increase in cost will be justified" is the
  literature-level warning that the SI iteration may stagnate
  or converge to a contaminated fixed point — neither of which
  invalidate his Fourier analysis but both of which are exactly
  what ORPHEUS Phase G is observing.
- Morel writes "iteration matrix" in the modern Krylov sense
  (one Jacobi step on the linear system), NOT in the
  source-iteration / Neumann-series sense of scattering-source
  iteration. Don't conflate the two when porting his statements.

## Adjacent references — verified

- **Bailey, Morel, Chang (2010), NSE 165, 149–169** — same
  research lineage but a different concern. This paper does an
  *asymptotic diffusion-limit* analysis of the converged SN
  scheme; it implicitly assumes the SI sweep delivers the
  discrete-ordinates system solution (so a single converged
  ψ_m field is the object of study). It does NOT revisit the
  SI/system equivalence question that Morel 1989 raises.
- **Morel and Montry (1984), Trans. Theor. Stat. Phys. 13** —
  the flux-dip / weighted-diamond derivation for curvilinear
  geometry. Cited from Bailey-Morel-Chang 2010 (Ref. 2) as
  the Galerkin-based diffusion-analysis that motivated the
  weighted-diamond coefficient choice. Not in
  `scratch/literature/` based on a quick check — would be a
  next pull if the user wants the canonical α-weight derivation.
- **Pomraning (1989), "The Transport Equation in General
  Geometry"** — gives the *continuous* curvilinear transport
  operator. Useful for setting up the continuous-to-discrete
  correspondence but does not discuss the SI sweep at all.
- **Lewis & Miller (1984) Ch. 4** — standard textbook
  derivation of the curvilinear SN sweep with the M-M α
  recurrence. Lewis & Miller present the sweep as if it
  delivers the system solution, **without** flagging the
  Morel-1989 lower-triangularity caveat. Worth cross-checking
  if Morel 1989 is the only published source pinning this
  exactly — if so, the ORPHEUS Sphinx page should cite
  Morel 1989 explicitly rather than Lewis & Miller.

## What to do next if this is insufficient

If the verdict needs more than Morel 1989's qualitative statement
(e.g., a quantitative bound on the SI/system fixed-point
displacement as a function of `n_quad` and `c`), the next pulls
are:

1. **Morel and Montry 1984** (Trans. Theor. Stat. Phys.) — the
   pre-Bailey curvilinear weighted-diamond paper. May contain
   the α-recurrence stability analysis Morel calls "simple
   Fourier analysis" in the 1989 paper without showing.
2. **Larsen and Morel (2010) Springer chapter** ("Advances in
   discrete-ordinates methodology"). Modern review; would
   explicitly address whether the SI sweep on curvilinear
   M-M is now considered settled or still has the Morel-1989
   asterisk.
3. **Adams and Larsen (2002), Prog. Nucl. Energy** ("Fast
   iterative methods for discrete-ordinates particle transport
   calculations"). Comprehensive review of SI vs. Krylov; if
   anyone has tabulated the fixed-point displacement for
   curvilinear M-M, it would be there.

None of these three are in `scratch/literature/` at the time of
this memo. If the user wants them next, a fresh
`literature-researcher` dispatch should pull them.

## ORPHEUS notation translation table

| Morel 1989 symbol               | ORPHEUS analog                                                   |
| ------------------------------- | ---------------------------------------------------------------- |
| `L` (Eq. 11)                    | `orpheus.sn.solver`'s implicit per-cell linear operator         |
| `L_C` (Eq. 12) — collocation    | spatial advection block, the part that lives in `orpheus.sn.sweep` (this is the lower-triangular causal piece) |
| `L_G` (Eq. 13) — Galerkin       | the full Galerkin redistribution block — equivalent to what would arise if M-M were re-derived via Galerkin angular projection rather than collocation at the quadrature points |
| `L_{S2}` (Eq. 14) — SN          | the standard M-M α-ladder used in ORPHEUS today |
| `ψ_m`                           | `psi[..., m]` indexed by `mu_index` |
| "angular leakage term"          | the `α_{m+1/2} ψ_{m+1/2} − α_{m-1/2} ψ_{m-1/2}` redistribution in `orpheus.sn.curvilinear` |
| "spatial leakage term"          | the `μ ∂/∂r (r² ψ)` divergence, handled in `orpheus.sn.sweep` |
| "removal term"                  | `σ_t ψ` |
| "inscatter term"                | scattering source `Σs φ` |
| "standard SN iteration sweep"   | `orpheus.sn.solver.sweep_si` (the legacy SI loop) |
| "off-diagonal components placed on RHS and iterated on" | the implicit lagging of cross-μ coupling in the SI fixed-point form `apply_matvec(ψ_prev) → ψ_next` |
| Morel's Eq. 12 lower-triangular `L_{S2}` | what `apply_matvec` constructs as its *true* linear operator (the Krylov solve sees the full operator, not just the lower-triangular sweep order) |

The mapping makes precise what the Phase G verdict claims: the SI
sweep implements the **lower-triangular** projection of the M-M
angular coupling matrix and lags the upper-triangular residual
to the previous iterate; the apply-matvec presents the **full**
M-M angular coupling matrix to Krylov in a single operator
evaluation. The two fixed points differ by exactly the
non-trivial action of the upper-triangular residual at the
discrete-iterate level, which Morel's "spectral radius
approaches unity" warning captures qualitatively.

See related: [[lessons]] (no relevant entries yet), and the
numerics-investigator verdict at
`.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_replan_verdict.md`.

---
name: streaming-apply-transpose-frame
description: Durable native frame for the transpose of a transport sweep operator Lᵀ — triangular-transpose (sweep=(I−N)⁻¹, matvec=(I−N), Lᵀ=(I−Nᵀ)=reverse DAG walk with per-cell coupling transposed + face roles swapped) + Lewis-Miller adjoint-ordinate (−Ω octant gives the reversed traversal FREE). THREE load-bearing gotchas a naive spatial-only reverse silently violates: (1) angular recurrence must transpose too, (2) pole/regularity seed transposes to an outer-face term, (3) trace block swaps inflow↔outflow.
metadata:
  type: project
---

# StreamingOperator.apply_transpose — native frame for Lᵀ

`apply_transpose` returns the plain Euclidean `Lᵀ`; the G-adjoint `A†=G⁻¹AᵀG` (with
`G_s=|Ω·n|·w_n`) is applied AROUND it by the adjoint wrapper. The transpose is NOT new
code — it is a THIRD application of the existing per-cell streaming algebra. (Wave O
Phase 4 / O.2b.)

## Two confirmed frames

1. **Triangular factorization (NLA).** The forward DAG-walk sweep = forward substitution
   on the unit-triangular `L=(I−N)` in cell-visit order (propagation
   `ψ_face_in = 2·ψ_cell − ψ_face_in`). Sweep solve = `(I−N)⁻¹`; matvec = `(I−N)`;
   transpose-matvec = `(I−Nᵀ)` = REVERSE cell-visit order with the per-cell coupling
   transposed AND the down/up face roles swapped. The `denom·ψ_cell` diagonal is
   SELF-transpose; only the off-diagonal upstream propagation reverses. Hits
   structure-exposing + expressive + algorithmic-advantage (matrix-free, same cost).

2. **Lewis–Miller adjoint-ordinate (time reversal).** `Lᵀ` = FORWARD streaming of the
   reflected-ordinate (`Ω→−Ω`) problem. The DAG-walk index generator is sign-parametric
   (forward `range(nx)` for μ≥0, reversed `range(nx-1,-1,-1)` for μ<0) — the −Ω octant
   ALREADY emits the reversed traversal. So `apply_transpose` is a thin SIGN RE-ROUTE of
   the forward sweep, NOT new traversal code.

## THREE GOTCHAS (load-bearing — the durable detection content)

A naive "reverse the spatial walk only" silently produces the WRONG operator. The
transpose of a transport sweep has THREE coupled pieces, not one:

1. **Angular transpose is the subtlety.** The curvilinear angular closure (Morel-Montry
   half-angle recurrence) is a SECOND triangular factor NESTED in the ordinate index,
   NOT symmetric. On transpose the angular recurrence REVERSES in ordinate index
   (`p←p'` ⇒ `p'←p`); the in/out angular coefficients SWAP roles;
   `upstream_per_ordinate ↔ downstream_per_ordinate`. Reversing only the spatial walk
   leaves the angular block UN-transposed → silently wrong operator. Needs a transpose
   variant of the cell angular contribution.

2. **Pole / regularity seed transposes to an OUTER-face term.** The pole seed is an
   `r=0` REGULARITY condition (NOT a BC). The Carlson seed reads outer inflow (the
   most-inward ordinate's outer value). Its transpose flows the inner-cell seed BACK to
   the outer-face adjoint term. "Seed is just 0 at r=0 on transpose" DROPS this coupling.
   Verify on a sphere N=2: the adjoint of the regularity seed MUST appear in the
   outer-face slot of `Lᵀ·y`.

3. **Trace block swaps INFLOW↔OUTFLOW.** Forward reads inflow to seed, writes the
   outflow residual + inflow identity. `Lᵀ` MIRRORS: reads OUTFLOW to seed, writes
   INFLOW. This is what makes `(L+C−B).H` compose — `B` is `A_ss:V_outflow→V_inflow`,
   self-adjoint under `|Ω·n|·w`. Swap the outflow/inflow trace selectors exactly as the
   boundary operator's own transpose already does.

## Cross-domain detection heuristic (durable, generalizable beyond SN)
When transposing ANY sweep/forward-substitution operator that nests a second recurrence
(here: angular-in-ordinate inside spatial-in-cell), the transpose must reverse BOTH
recurrences and swap the role of every directional read/write (inflow↔outflow,
seed↔residual, in-coefficient↔out-coefficient). A spatial-only reverse is the classic
silent bug. The clean construction is never a hand-rolled backward loop — it is the SAME
per-cell primitive re-applied with reversed traversal + transposed local coupling. The
isolating diagnostic is a small DENSE build: form `L`, take `L.T`, assert the
reverse-walk recurrence reproduces it to ~1e-12 on the SPATIAL traversal alone, THEN
layer the angular transpose.

## Adjoint-for-free theorem (cross-poll, Krylov/eigenvalue)
Once `Lᵀ` lands, `(L+C−B).H = L.H+C.H−B.H` is a THEOREM (biproduct + dagger functor,
[[issue_208_operator_algebra_frames]]): `C.H=C` (self-adjoint collision), `B.H` exists,
the adjoint wrapper distributes G over the sum. MoC borrowing: the adjoint-face recurrence
= MoC backward-ray transmission transpose (lower→upper triangular).

## Smell #16 — fires BEFORE the code is written (confirmed this problem class)
A hand-rolled adjoint sweep would be a THIRD path to the per-cell operator already shared
by the DD residual + the vectorized sweep — the twin-path smell fires before any code
exists. Native fix: `Lᵀ` is a re-application of the shared streaming primitive, not a
twin. (One of the seven Smell-16 sightings now promoted to AGENT.md.)

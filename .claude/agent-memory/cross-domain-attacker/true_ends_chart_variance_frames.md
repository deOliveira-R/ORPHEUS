---
name: true-ends-chart-variance-frames
description: #707 W5 attack (2026-11-11) — condensed storage makes .H = dual map in the density chart; chart must live on the SPACE; side is per-AXIS variance; frame dual law M^T = ♭RD exact where the Hilbert square fails; Λ breaks the head kernel condition on under-resolved rules
metadata:
  type: project
---

Verdicts (memo `scratch/true_ends/review/xda.md`; probes `scratch/true_ends/review/xda/`):

- `A = ♯A_int` stored ⇒ `G⁺AᵀG` IS the dual map in the density chart. This is exact for every symmetric G. The involution is exact only in the INTEGRATED chart; in the density chart it is `PAP`. So the 100× option is to store V* integrated.
- The Riesz legs are built from SPACES at the adjoint's construction ⇒ a chart that is a member property is unrealisable; key it on `DualSpace`.
- Mixed variance: the moment head is covariant (M writes C*), cells are contravariant, and the group axis is integrated but hidden by `G = I`. `TensorProductSpace` already accepts `X*⊗Y`.
- The frame dual law `Mᵀ = W·R·diag(1/(2ℓ+1))` holds to 1e-16 on 7/7 rows, including the under-resolved ones. ⇒ The S_N scattering adjoint "transpose Σ_ℓ" follows in one line, because D commutes with Λ.
- The production Λ fails the kernel condition on LS4 L=5 (9e-2) and on product(4,4) L=4 (0.16), because the alias kernel mixes ℓ. In the same rows R does not kill ker G_head (0.14–0.35): the head serves two variances.

**Why:** this is the first sighting of "a binary side is the all-axes special case of per-axis variance". Promote it at a second sighting in a different problem class.

**How to apply:** before accepting any "space side / role" design over a tensor product, tabulate the variance per axis and run `‖R(I−P)‖` and `‖(I−P)M‖` on an under-resolved rule. [[frame-bridge-pullback-riesz-fibre-metric]] (D24) is the sibling.

**2026-11-12 follow-up (memo `scratch/true_ends/review/variance_xda.md`):** "density vs integrated" per axis welds TWO data: TYPE (Radon–Nikodym: measure pushes forward / function pulls back — read by CONDENSATION, `mixture.py` condense: g_to & χ @T, vectors & g_from averaged, ψ† pair-weighted) and INDEX (upper/lower wrt λ_a — read by the PAIRING). On energy λ = counting ⇒ index label has 0 consumers (`E == E*` False, metric None both); ψ† ≠ ψ in TYPE on EVERY axis. Energy counting = no Haar measure on energy (dE vs du both consumed as VIEWS). Data side = the method's BALANCE (test) axes, per equation not per operator (S local in space flips; integrates angle, doesn't). A2 blocker [M]: `V.dual()` / R8 product lose `.axis()` (TypeError) ⇒ index belongs in the axes record (`Axis.dual()`), not DualSpace factors.
**How to apply:** on any "variance / side / density" brief, write the 4-row table (TYPE × chart → INDEX, condensation rule) before accepting a binary flag; test homogenisation (flux vs σ both upper in space).

**2026-11-12 second attack (memo `scratch/true_ends/review/form_degree_xda.md`): "form degree is the native label".** Verdict: RENAMING. Degree read consistently = the PAIR (native degree, stored degree), each in {0, top}, = (TYPE, INDEX), with ⋆ the chart; one degree per axis = stored where a ⋆ exists, native on energy (no ⋆) ⇒ fails homogenisation (`sn/solution.py:501–530`: flux V, σ Vφ, χ production weights, all stored 0). Gains, all blind in A2: (a) energy — ψ ∈ E*, ψ† ∈ E reconciles R8 "spatial only" with "A⊗(E⊗X)*" (they name different reference spaces); SELF-CORRECTED my prior "INDEX same on energy" (no ⋆ ⇒ index has no content apart from type); (b) R4/#705/#706 = ⋆ degeneracy/indefiniteness; (c) curvilinear α recursion IS d on the μ-complex. Streaming "metric-free because d" FALSE [M]: collision C=diag(Vσ) has the same canonical-transpose property; L = d∘K_Ω, K_Ω upwind contraction is not a ⋆. Breaks: harmonic head (no complex), J± (swap, not negate), LD = FEEC Galerkin Hodge not DEC.
**How to apply:** on a "degree / grading" brief, write (native, stored) per axis first; a single integer per axis is a weld.

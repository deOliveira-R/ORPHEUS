---
name: frame-bridge-pullback-riesz-fibre-metric
description: 2026-10-18 attack on numerics_untangling (a)/(b)/(c) — frame = sampling map + canonical dual (pullback metric); Riesz datum has two homes; LD curvilinear = trivial bundle with fibre metric g(tau); Gauss r^{d-2} nodes diagonalise both LD forms
metadata:
  type: project
---

Memo: `scratch/numerics_untangling/frame_bridge/attack.md` (probes in `attack_probes/`).

- (a) Frame = sampling map Φ = ι*|_U + canonical dual; kind by rank: iso 3/33, isometric embedding 30/33 on the fork grid. Bridged modal metric = PULLBACK Φᵀ W Φ ⇒ the modal axis is downstream of the frame. 7 of 10 production frames are PG coarse-graining (nodal↔nodal), refuting "frame only when one factor is carried both ways".
- (b) The fork R-6.2c-1 compared against a VARIANCE-wrong arm (D24). Real discrete-vs-continuum choice: 1 grid row (folded_product(2,4) L=2, iso = Morel GQ). Frame reads Riesz datum twice (Basis.reconstruct's measure-free (2ℓ+1) vs basis_space's G_d⁺): M.H=R/W misses 0.365 there. Coherent bridges: GQ1 / GQ2-extension / standard SN (Strang crime).
- (c) Discrete base ⇒ trivial bundle; metric G = ⊕ V_c ĝ(τ_c), τ = r_in/h; Kronecker iff constant volume element. Guard's "value unknown" refuted by R-6.2c-1. Lagrange at 2-pt Gauss of r^{max(d−2,0)} diagonalises M and R exactly (≤1.3e-13): same consistent scheme, NOT lumping (Palmer–Adams FL = row-sum on Lobatto nodes, a different operator).

**Why:** the user ruled the code is never the oracle (ruling 14); these readings feed the numerics concept order.
**How to apply:** before any metric/frame/LD-curvilinear brief, re-ground the row counts against the tree; status of rulings lives in `.claude/plans/numerics_untangling.md`.

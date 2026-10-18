---
name: axis-collapse-unit-and-gspace-invariance
description: 2026-10-18 numerics-untangling F3/F6 attack — collapse pair = <unit,.>_G ⊗ id + Moore-Penrose section (axis owns unit, space owns lift, frame = oracle); invariance = relabelling ∘ pushforward ∘ support.action (G-space on the MANIFOLD; AngularMeasure refuted both ways)
metadata:
  type: project
---

Memo: `scratch/numerics_untangling/ontology_f3_f6.md` (probes beside it).

F3: `_collapse_pair` builds a rank-one frame and reads back exactly its inputs (weights array_equal; Σw 0–6 ULP off np.sum via the Gram einsum). 6/6 prod collapses are "angular" on `of_axes` products no frame minted. Universal form: unit = coordinates of the constant 1 on the axis; R = ⟨unit,·⟩_G ⊗ id; E = R⁺ (ν-free); modal Legendre exact (0.0), the docstring's refused covector G·1 off 9.1 %. Axis owns unit (F1's "induced data by name"), space owns the lift (positioned-form admission is a product clause). Agrees with #525 item 9, not with its "frame object" ownership.

F6: kernel ≡ match∘pushforward∘induced_action, 75/75 finite rows (45 True), Hungarian matcher, perturbation control. Zero-padding manufactures verdicts (EnergyGroups/IndexSet O2_x True; chart μ-rule SO2_x True, SO2_y False though it is a z-product factor). Prod-entered closures 1315/1315 on the sphere family (gates.numerics census). The action belongs to the support (G-space); measure owes only relabelling-equality (#403).

**Why:** the user's F3/F6 questions in `.claude/plans/numerics_untangling.md` rulings 6, 9.
**How to apply:** open before any brief on axis collapse, marginal, section, invariance, quotient, `Manifold.quotient` or a measure specialisation; re-ground the plan's rulings first.

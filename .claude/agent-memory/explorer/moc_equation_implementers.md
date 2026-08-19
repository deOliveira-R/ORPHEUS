---
name: moc-equation-implementers
description: The MoC theory page's equation→code truth — the kernel lives in exactly TWO places, the page carries NO vv-status rationale comments, and two measured defect leads (dual t_s^eff denominator, MMS trig-chart directions).
metadata:
  type: project
---

Measured 2026-08-18 (`HEAD` 9ba35a8f → re-verified 58e46c6f) while declaring
`implements` edges for 11 equations on
`docs/theory/methods/method_of_characteristics.rst`. Findings written to
`.claude/inventories/declare_moc.md`. Related: [[nexus82-operator-algebra-implementers]],
[[nexus82-loss-representation-implementers]].

**Why:** the V&V ledger needs declared equation→code links; 216 coverage claims
ride on these 11 equations, and the graph had **0** `implements` edges on all of
them (the name-token heuristic matches nothing — MoC labels share no token with
any symbol name).

**How to apply:** read this before any MoC equation/coverage dispatch; the
three durable claims below survive line-number churn.

## 1. The MoC kernel lives in exactly TWO places — declare BOTH, always

- `orpheus.moc.core.MOCSolver.solve_fixed_source` — production sweep.
- `orpheus.derivations.continuous.mms.moc.mms_sweep` — an **independent**
  re-implementation of the same segment kernel with a manufactured per-segment
  source, imported by `tests/moc/test_mms.py`. It is not a wrapper; it duplicates
  τ, Δψ, the `4π ω_a ω_p t_s sinθ` weight and the flux reconstruction.
- Third, symbolic: `orpheus.derivations.discrete.moc.equations.derive_bar_psi`
  (SymPy, integrates + asserts the ψ̄ identity — the statement of record for
  `bar-psi`, `delta-psi`, `attenuation`).
- ⚠ Its sibling `derive_scalar_flux_weight` is **`print`-only** — no SymPy, no
  value. It derives `boyd-eq-45` and implements nothing; a claim routed through
  it can never be corroborated by a coverage run.
- ⭐ Consequence: a "who implements MoC equation X" answer that names only
  `orpheus/moc/` is **half an answer**, and the missing half is what the L1 MMS
  gate executes.

## 2. This page carries NO `(vv-status rationale)` comments — the method's step 1 is page-dependent

`[M]` `grep -n "vv-status" <the MoC page>` → **5** hits, all bare
`.. vv-status: <label> documented` status stamps, and **none** on any of the 11
equations. The rationale form exists on `methods/index.rst` and
`methods/collision_probability.rst`. ⟹ before briefing (or believing) a fanout
method that leads with "read the authored rationale", grep the page for it; on
this page the answer comes from the equation + prose + the **claiming test
module's own comment block** (`tests/moc/test_*.py` pytestmark comments name the
gate and often the mechanism — that is where MoC's authored knowledge sits).

## 3. Two measured defect leads (independent of the declarations)

- **`t_s^eff` has two denominators.** `MOCMesh._generate_tracks` places rays with
  `(t_max−t_min)/n_rays` (the documented `effective-spacing`) but STORES
  `(t_max−t_min)/len(tracks_kept)` — the value the `boyd-eq-45` weight reads.
  They differ whenever a degenerate ray is dropped. `[M]` bit-identical on
  **16 of 16** azimuthal indices at the shipped default (standard PWR
  Wigner-Seitz cell, `n_azi=16`, `n_polar=3`, `ray_spacing=0.05`, reflective) —
  latent, not active.
- **`mms_sweep` evaluates trig on the angle chart** (`np.cos(quad.phi[a])`)
  while `MOCMesh._generate_tracks` reads the exact #325 points
  (`quad.cos_phi[a]`, docstring: *"never trig on the chart"*). `[M]` bitwise
  different on **7/15/25 of 16/32/64** components at `n_azi = 8/16/32`
  (max |Δ| 3.3e-16), and the chart **loses** the exact mirror identity
  `cos_phi[n−1−k] == −cos_phi[k]` (True on points, False on chart, all 3 orders).
  ⟹ the MMS reference sweeps along directions that are not exactly the ones its
  own tracks were traced with. Looks like a residual #325 site.

## 4. The scoping ruling that recurs

The same formula appears under a DIFFERENT chart on another theory page, and
declaring it under this page's label is wrong granularity. Two instances measured
here: the flat-source `(q/Σ_t)(1−e^{−τ})` at **12** sites under
`derivations/continuous/trajectory_resolvent/origins/specular/greens_function*.py`
(own page: `theory/references/trajectory_resolvent.rst`), and the ray-circle
quadratic in polar/chord form (`ρ = −r cos ω ± √Δ`) at ~40 sites under
`peierls_nystrom/` + `chord_oracle.py`. ⚠ I over-claimed exhaustiveness TWICE by
grepping `orpheus/` with `| grep -v derivations` — the excluded half was not
empty. Never filter `derivations/` out of an exhaustiveness grep; filter the hits
by MEANING afterwards.

---
name: psi-half-seed-angular-trace-frames
description: #282/#280 R10 ψ½(μ=−1) starting-direction augmentation — the seed is the ANGULAR-INFLOW TRACE of the (r,μ) phase-space rectangle; verdict on its type + metric + sphere/cylinder keying
metadata:
  type: project
---

# ψ½ starting-direction seed = the angular-inflow trace of phase space (#282 R10)

Attack on the AUGMENTED-STATE framing for the M-M starting-direction flux
ψ(·,μ=−1) (branch `refactor/sn-walk-unification`, roadmap Phase 2.5d, ruling
R10). Branch-verified reads (L-005): `orpheus/sn/spatial/pole_angular_closure.py`,
`psi_half_angle_seed.py`, `orpheus/geometry/reduced_operator.py`,
`orpheus/transport/full_field.py`, `orpheus/numerics/spaces/full_field_space.py`.

**Why:** design crux the execution task routes to detection — where does ψ½ live,
what type, what metric, sphere-vs-cylinder uniformity.

**How to apply:** the durable verdict for the 2.5d carrier augmentation + the
frame refutations. Fire on any "where does this boundary/seed DOF live" question
in curvilinear transport.

## The central structural verdict

The phase space is the (r,μ) rectangle. Its boundary splits ∂spatial (r=R edge,
all μ) ⊔ ∂angular (μ=±1 edges, all r). The existing `FullField = bulk ⊕ ∂spatial`
carries only ONE trace because Cartesian SN has no angular advection. Curvilinear
SN's `(1−μ²)/r ∂_μ` term IS a first-order ADVECTION IN μ with inflow at μ=−1 (a
SWEEP in angle). ψ½ is the ANGULAR-INFLOW TRACE — the μ=−1 edge — genuinely DUAL
to ∂spatial in the rectangle (F2's rectangle-edge duality is REAL). So the third
block is NOT a plain bolt-on (F1) — it is the native second half of the phase-space
boundary trace.

## Metric verdict (the load-bearing correction)

ψ½ carries ZERO angular-quadrature metric. The angular advective flux (1−μ²)
VANISHES at μ=−1 — the SAME structural fact as α_{1/2}=0 (the Bailey dome closing
at the endpoint) and F3's "w=0 kills the scattering sum". So the seed is a metric
GHOST, invisible to every angular-MOMENT functional (φ*₀, k, reaction rates) —
reciprocity ⟨Aψ,φ⟩_G legitimately does not see it (F2's zero-metric guess RIGHT;
reuse the tangential-trace-slot zero+pseudo-inverse pattern in
`FullFieldSpace`/`AngularTraceSpace`). F2's own guessed nonzero dual weight (the
redistribution α at the boundary) SELF-REFUTES: α_{1/2}=0.

DO NOT fabricate a nonzero VOLUME metric to "fix" the G-adjoint coupling. The
seed's radial dynamics (its own inward advection, Hébert 3.432–3.435) live in the
OPERATOR block A_seed, NOT in the metric G. A volume metric would over-weight a
zero-measure μ-point and make reciprocity see a non-physical quantity. The seed's
coupling to the bulk is preserved EUCLIDEAN-ly by the roadmap's already-specified
`apply_transpose` (metric applied only at the outer G-wrap: `A.inverse().H.apply(b)
= G⁺·apply_transpose(G·b)`). CONSISTENT: the Euclidean core carries the coupling to
the metric-visible bulk/trace adjoint; only the seed's OWN adjoint component is
ghost-zeroed by G⁺ — harmless to moment consumers, LATENT wrong-answer for a future
pointwise-μ=−1 adjoint consumer (adjoint detector at the pole direction / full
angular-adjoint reconstruction). This is the sharp risk: UNLIKE the tangential trace
slots (identically zero in output), the seed is metric-invisible-YET-ACTIVE (nonzero
in forward output) — the `FullFieldSpace` docstring's exact-pseudo-inverse
justification ("tangential slots identically zero in every matvec output") does NOT
transfer; it holds for moment consumers only.

## Sphere-vs-cylinder uniformity ruling

Seed is a GENUINE independent DOF ONLY on SPHERE (μ=−1 ∉ GL interior nodes). On
CYLINDER the level's inflow edge η=−sinθ COINCIDES with the first azimuthal ordinate
(η₀=−sinθ ⇒ τ_raw=0, clamped — the reason the cylinder clamp exists, #229) ⇒ the
seed IS ψ₀ ⇒ a DEAD DOF if carried (the task's "cylinder telescopes, dead state").
Cartesian: no angular advection, no seed. So the carrier augmentation is
CURVATURE+QUADRATURE-keyed: present iff `mu_start ∉ level μ-nodes` (equiv. τ_raw≠0).
Empty seed block on cylinder/Cartesian — same pattern as the pole being "not a face".

## Frame dispositions

- **F1 (third block)** — CONTAINER correct, silent on identity. Completed by F2.
- **F2 (angular-boundary-trace)** — GEOMETRY correct (rectangle-edge duality, ∂angular
  sibling of ∂spatial; the OUTWARD μ=+1 direction IS derived outflow — its r=0 value
  is the pole reflection ψ(0,+1)=ψ(0,−1) of the inward seed, NOT independent state;
  the seed's own inflow BC is the corner (R,−1) = ∂spatial's μ=−1 slice → the two
  traces are corner-coupled). METRIC-symmetry WRONG: the angular trace weight is 0,
  not a partial current |μ|w. Keep F2's geometry, reject its metric symmetry.
- **F3 (zero-weight ordinate)** — REJECTED. Mislocates a BOUNDARY object as an INTERIOR
  ordinate ⇒ breaks every (N,…) shape contract, node-count + reflection_index pairing,
  level_indices. The w=0-auto-kills-scattering win is a shallow coincidence. Sole
  durable contribution: confirms the seed's angular weight is 0 (⇒ Frame-A ghost metric).
- **F4 (descriptor/DAE)** — CONFIRMING lens. Seed = index-0 DIFFERENTIAL (not algebraic)
  variable (own ∂/∂r), block-triangular ⇒ forward-substitutable, no DAE solver. Latent
  warning: a 0-D/infinite-medium limit that drops ∂/∂r makes it index-1 (algebraic
  constraint) needing different handling.

## Recommendation

F1⊕F2 synthesis: third `FullField` block = the angular-inflow trace; ZERO-metric ghost
(reuse tangential-trace pseudo-inverse); Euclidean coupling via `apply_transpose`;
sphere-only (curvature+quadrature-keyed). Makes A block-lower-triangular [seed; μ-inc
bulk] (turns #282's lagged back-edge into a forward edge), 2.5b reverse-scan = clean DAG
reversal, 2.5c = the existing G-wrap mechanism with no new metric math. Alternative if
carrier-shape variability bites: SCHUR-eliminate the seed (its scattering-coupling folds
into S, source-part into q; (L+C) triangular, S denser; no seed DOF) — Option C.

## New elegance smell (candidate, 1st sighting)

"metric-invisible-yet-active DOF": a zero-quadrature-weight state that is nonzero in the
forward output and couples to metric-visible DOFs. The zero-metric pseudo-inverse is
honest for the intended (moment) functionals but its usual justification ("identically
zero in output") does NOT hold — the DOF is active, only weightless. Distinct from the
tangential-trace zero-in-zero-out slot. Held for a 2nd sighting before Part C promotion.

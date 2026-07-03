# Diffusion scalar-trace convention crosswalk (#290 P2 — L17, written BEFORE code)

Every producer/consumer pair below must match column-for-column. The Bridge row names
the ONE site that owns each convention (Pattern 7 — normalise at the definition site).

## The (J⁺, J⁻) ↔ (φ_b, J) dictionary (P1 half-range moments)

| Quantity | Definition | Inverse |
|---|---|---|
| J⁺ (outflow partial current) | `φ_b/4 + J/2` | `φ_b = 2(J⁺ + J⁻)` |
| J⁻ (inflow partial current)  | `φ_b/4 − J/2` | `J = J⁺ − J⁻` |

Owner of the derived quantities (`net_current`, `p1_boundary_scalar_flux`):
**`ScalarBoundaryFlux` methods** — no consumer re-derives the dictionary.

## Storage layout

| Convention | Value | Owner (single site) |
|---|---|---|
| Face slot shape | `(2, ng, *face_spatial)` — component axis 0, group axis 1 | `ScalarTraceSpace.for_faces` builds the slots |
| Component rows | row 0 = J⁺ (outflow), row 1 = J⁻ (inflow) | `ScalarTraceSpace.OUTFLOW_ROW / INFLOW_ROW` ClassVars |
| Flat packing | faces in `face_labels(axes)` order (axis asc, endpoint order); each slot C-raveled | `FaceLayout` (existing, reused unchanged) |
| Face names | `"{axis}{min|max}"`; radial `"outer"` renders `max`; pole is NOT a face | `FaceLabel.face_name` (existing) |

## Orientation / sign

All trace quantities are **face-local against the OUTWARD normal n̂_f**: J⁺ ≥ 0 leaves
the domain, J⁻ ≥ 0 enters, net `J = J⁺ − J⁻` > 0 = leakage out — at EVERY face
(xmin included). Conversion to a global axis-vector current (`J_z`: `+J_face` at xmax,
`−J_face` at xmin) happens ONLY at consumers that need the vector form (the P4
`LeakageOperator` assembly); the trace never stores axis-signed values.

## BC family (the A_ss albedo operator, P3 realizer)

`J⁻ = 𝒜 · J⁺` per face per group:

| Law | 𝒜 | Notes |
|---|---|---|
| vacuum | 0 | = Marshak zero-incoming (user ruling 3: vacuum MEANS J⁻=0) |
| reflective | 1 | J = 0 |
| white | 1 | coincides with reflective at P1 (isotropic return preserves current) — document |
| albedo(α) | α | physical range [0, 1] |
| zero-flux Dirichlet | −1 | φ_b = 0; the mathematical idealization (NEW law, honestly named); breaks J± ≥ 0 (positivity is a property of PHYSICAL laws 𝒜 ∈ [0,1], not a type invariant) |

## Metric

| Space | Weight | Meaning |
|---|---|---|
| SN `AngularTraceSpace` | `G_s = |Ω·n̂_f| ⊙ w_n` per ordinate | partial-current pairing of angular DENSITIES |
| `ScalarTraceSpace` | boundary face AREA per slot (slab 1, cyl 2πR, sph 4πR²) | surface measure for already-angle-INTEGRATED currents |

Consistency claim (DSA #2 seam): the scalar trace is the ℓ=0 half-range moment of the
SN trace under G_s; the M_half restriction operator owns that reduction — the scalar
metric does NOT try to encode the angular weights (they are already integrated out).

## Shared with bulk (unchanged project-wide conventions)

- Group ordering: index 0 = fastest (descending energy); downscatter upper-triangular.
- Units: J± `1/(cm²·s)` — same dimensional family as ScalarFlux (areal rate, eV-free).
- Composite flat: `FullField.to_flat` = bulk C-ravel ⊕ trace flat buffer (existing);
  scalar composite dim = `ng·prod(spatial) + Σ_faces 2·ng·prod(face_spatial)`.
- Mesh identity: `bulk.mesh is boundary.mesh` (the existing composite invariant is the
  scalar/angular mixing guard for free — an SNMesh bulk can never pair a MaterialMesh trace).

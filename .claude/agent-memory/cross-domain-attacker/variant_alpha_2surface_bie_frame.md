---
name: Variant α 2-surface BIE / monodromy frame
description: Match between Variant α extension to slab/annulus/hollow-sphere and boundary-to-boundary scattering operator S with T = (I − S)^{-1}. Promotes existing scalar T to a named operator resolvent.
type: project
---

Frame match: Variant α (Sanchez Eq. A1 closed-period bouncing
characteristic Green's function) → boundary-to-boundary
scattering operator S on surface phase-space, with full
Green's = F + B·(I − S)^{-1}.

**Why:** the existing sphere code expresses T(μ_surf) = 1/(1 −
α·exp(−Σ_t·L_period)) as a scalar geometric-series shortcut.
That shortcut is the rank-1 case of the resolvent of a
multiplication operator on surface phase-space. The 2-surface
generalization (slab-2BC, annulus, hollow-sphere) is the rank-2
block case of the SAME resolvent — not three new derivations.

**How to apply:** before extending Variant α from sphere to
cylinder/slab/annulus/hollow-sphere:

1. Refactor sphere code to name `S = α · R_chord` and `T = (I − S)^{-1}`.
2. Rewrite V_α1/2/3 SymPy proofs as resolvent-of-S identities
   (V_α1 = direct resolvent expansion; V_α2 = Schur complement /
   matrix-element identity for T_00; V_α3 = S = 0 trivial reduction).
3. Cylinder extension = same rank-1 S, only chord algebra changes
   (ρ_max(r, θ), r'(r, ρ, θ)) — share assembly with MOC.
4. Slab/annulus/hollow-sphere extension = rank-2 block S; vacuum
   on a face = zeroing one row of S; mixed BCs are diagonal entries
   of S.
5. Predict per-geometry pathology by ess_range(M_period) ∋ 1
   (Phase 5 A.3 multiplication-operator caveat — sphere α=1 grazing
   ray μ→0 is unbounded; slab is immune; cylinder/annulus inherit
   the pathology iff a vanishing-chord locus exists at the
   boundary).

**Cross-references:**
- scripts/validated_hilbert_schmidt_separable.md — same lever
  family (resolvent / multiplication operator), Phase 5 precedent.
- scripts/validated_unified_geometry.md — geometry-as-parameter
  pattern; the BIE frame extends it from static geometry to
  BC-absorbed Green's function.
- scripts/validated_bc_tensor_network.md — BC-as-tensor pattern;
  the boundary-to-boundary S IS the BC tensor network's contraction
  along surface legs.

**Promotion candidate:** when this refactor ships and the rank-2
S identity for slab is verified in SymPy, promote to
scripts/validated_boundary_scattering_resolvent.md and add a
trigger row to A.7 ("boundary-to-boundary scattering operator —
trigger: BC-absorbed Green's function with surface-only resolvent
structure").

**Rejected candidates (low-signal for this problem class):**
- Floquet/Bloch — bouncing is path-length-periodic, not
  lattice-periodic.
- Symbolic dynamics / ergodic billiards — 1D radial billiards
  are integrable (impact parameter conserved); would fire for
  irregular 2D cells.
- Symplectic geometry — no field, no drift, no long-trajectory
  preservation issue; trivial.
- Method of images — equivalent to the resolvent expansion of S
  for slab; not an independent frame, just a sanity check.

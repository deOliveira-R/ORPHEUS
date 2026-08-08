---
name: harmonic-frame-folded-quadrature-plumbing
description: Durable shape of the SN anisotropic-scattering harmonic-frame layer (no computed Gram anywhere on the kernel path; three shape-contract tiers; sigma_y parity rule per SH slot) — mapped for the Q5.6 folded-quadrature campaign
metadata:
  type: project
---

# Harmonic frame + folded-quadrature plumbing — durable shape

Mapped 2026-08-07 on `refactor/operator-strategy-layers` @ `f4029bb4` (Q5.6
campaign, σ_y-folded cylindrical quadrature). Line numbers re-derive via Nexus;
the SHAPE below is durable.

**Normalization architecture (the load-bearing fact).** The scattering kernel
`S_aniso = (1/W)·R∘Λ∘M` has **NO computed Gram and no pseudo-inverse anywhere**:
`M` = raw `einsum("n,nlm,n...->lm...", w, Y, ψ)` (no division;
`SphericalHarmonicBasis.analyze`), `R` carries the STORED CONTINUUM dual factor
`(2l+1) = 4π·g_C⁻¹` (`addition_theorem_factor`), and the scalar `1/W` sits at the
producer boundary (`build_aniso_source` / the windowed arm), NOT inside the
kernel. ⟹ `R∘M ≈ projector` holds ONLY by assumed quadrature-exactness of the
discrete Gram — never checked at frame construction. On a folded rule the
odd-even cross-Gram is non-zero (measured: flat flux → slot [1,2] reads +6.49 =
Σw·μ_y ≈ ∫|ξ| ≈ 2π), which IS the Q5.6 defect. `FrameBase.gram`/`project`
(row-sum probe) are the homogenise/condense verbs only — not on the scattering
path.

**Three shape-contract tiers** (decides sub-basis design):
1. **Table-derived (adapts):** the einsum outputs size from the table's
   `(N, L+1, 2L+1)`; the 2-D moment-emit buffer derives
   `n_l, n_m = moment_frame.table.shape[1:3]`.
2. **Scalar-L-derived (hardcodes rectangular):** `MomentField._phase_space_shape`
   = `(L+1, 2L+1, …)`; `SphericalHarmonicSpace.__post_init__` VALIDATES
   `shape == (L+1, 2L+1)`; `HarmonicMomentFlux.truncate`.
3. **Fixed-slot indexing:** DSA reads `angular_frame(1).table[:, 1, 1]` (=μ_x=η,
   σ_y-EVEN); `scalar_flux` reads `[0,0]`; `l_block`/Λ slice `[l, :2l+1]` (Λ is
   m-BLIND — m is an einsum spectator). ⟹ a fewer-function FLAT basis breaks
   tiers 2+3 and every "nlm" einsum; a rectangular even-sub-basis with ZEROED
   odd table columns flows through everything.

**σ_y parity rule (basis chart: polar axis μ_x, cos φ_b = μ_y/sinθ, sin φ_b =
μ_z/sinθ; ξ = μ_y, folded reps on ξ>0):** slot `[l, l+m]` (m>0, cos branch) odd
iff m odd; slot `[l, l−|m|]` (m<0, sin branch) odd iff |m| even; m=0 even.
Measured L=3: odd = {[1,2]}, {[2,0],[2,3]}, {[3,1],[3,4],[3,6]}. ⚠ The
quadrature-chart phrase "sin-φ branch" does NOT map to the basis's sin slots —
the basis azimuth is measured FROM μ_y.

**Cylinder-path consumer set is minimal:** production `angular_frame`/
`spherical_harmonics` consumers = `ScatteringOperator.frame`/`.Y` + DSA (slab-only
admission, #282-blocked on curvilinear). Windowing/moment-emit is 2-D-Cartesian-
only (`_maybe_window`: `is_cartesian and ndim == 2`; the 1-D scan RAISES on
`moment_frame`); the 1-D cylinder SI iterate is full `AngularFlux` — moments are
transient einsum intermediates. Adjoint rides the SAME `S.frame` object
(`full_scatter_kernel`) — restricting the frame at source covers it. Diffusion:
zero SH consumers.

**Test-coverage traps (route-equivalence, L-013/L-017 class):** the only
cylinder-P1 gate (`test_cylindrical_p1_source_matches_hand_reference`, L0) shares
the quadrature between SUT and hand-ref — on a folded rule BOTH compute the same
garbage ξ-odd moment and stay green. NO converged P1-scattering cylinder solve
test exists (P1 solve/keff tests are slab+sphere+2-D-Cartesian; the "#229
xfail-strict ladders" premise is STALE — xfails removed 2026-06-13, and those
ladders are the aniso-ANSATZ case with ISOTROPIC scattering,
`build_cylindrical_anisotropic_mms_case` = `Quadrature.product`, path (I) ≠ the
P1-scattering path (II)). The ξ-mirror adjudicator is
`test_azimuthal_mirror_symmetry.py` (P0, full product rule).

---
name: MMS anisotropic curvilinear closeout
description: feat/derivations/mms-anisotropic-curvilinear branch closeout — Phase 1 of SN test-architecture campaign defending vv-principles failure mode #7
type: project
---

`feat/derivations/mms-anisotropic-curvilinear` 2026-05-05. **Phase 1 of SN test-architecture campaign**: companion anisotropic curvilinear MMS factories that activate the angular-redistribution operator missed by the existing isotropic ansatz (vv-principles failure mode #7).

**Branch 1 (SymPy)** — `derive_spherical_anisotropic_mms()` and `derive_cylindrical_anisotropic_mms()` in `orpheus/derivations/continuous/mms/sn.py`. Substitute ψ_n = (A + B·ζ)/W (ζ=μ for sphere, ζ=η for cylinder; A=sin(πr/R), B=(r/R)(1-r/R)cos(πr/R)) into the continuous SN operator and prove `simplify(LHS-RHS) = 0` algebraically.

**Branch 2 (numpy)** — `SNSphericalAnisotropicMMSCase` / `SNCylindricalAnisotropicMMSCase` evaluate the closed-form Q^ext_n(r) per ordinate via vectorised numpy. Mesh + materials methods mirror the isotropic siblings exactly so consumer L1 tests can swap ansatz with one factory call.

**L1 cross-check (the gate)** — Branch-2 numpy ↔ Branch-1 SymPy `lambdify` agreement: max|Q1-Q2| = **2.220e-16** (sphere), **2.220e-16** (cylinder) at machine precision. The cylindrical L1 test additionally verifies the (η,ξ) ↔ (θ,φ) convention round-trip via `arctan2(ξ,η) = φ_az` and `arcsin(sqrt(η²+ξ²)) = θ`.

**Key math findings**:

- Spherical Q^ext_n: `μ A' + μ² B' + (1-μ²) B/r + (Σ_t-Σ_s) A + Σ_t μ B`. The `(1-μ²) B/r` term IS the angular redistribution analytic — verified by Bailey 2009 Eq. 7-10.
- Cylindrical Q^ext_n: `η A' + η² B' + ξ² B/r + (Σ_t-Σ_s) A + Σ_t η B`. The `ξ² B/r` term is the cylindrical analog. Derivation handles `ξ = ξ(η,φ)` and `ψ = ψ(r,η(θ,φ))` correctly: `∂(ξψ)/∂φ = (∂ξ/∂φ)ψ + ξ(∂ψ/∂φ) = η·(A+Bη) + ξ·(-Bξ)` then `-(1/r)∂(ξψ)/∂φ` cancels the streaming `ηA/r + η²B/r` and leaves `ξ²B/r`.
- Choice of B(r): `(r/R)(1-r/R)cos(πr/R)`. Vanishes at r∈{0,R} via the `r(R-r)` envelope. The cosine factor is non-trivial and produces extrema at r=R/3 (cos(π/3)=1/2), avoiding the trivial r=R/2 zero of cos(π·1/2). The non-zero-redistribution test gate exercises r=R/3 specifically.

**Bug-revealing potential**: the isotropic curvilinear MMS test in `tests/sn/test_mms_curvilinear.py` is currently FAILING on main with orders ≈ 0 instead of O(h²) (errors flatlining at 0.107-0.111 across n_cells=10,20,40,80,160). This is exactly ERR-026 (curvilinear sweep WDD wrong fixed point). The Phase-0 main-agent task to write the consumer `tests/sn/_l1/test_mms_spherical_anisotropic_dd_convergence_O_h2.py` will determine whether ERR-026 also fails the anisotropic case (expected: yes, possibly with different order signature so the pair narrows down where in the sweep the bug lives).

**4 new equation labels** in `docs/theory/discrete_ordinates.rst` §"Curvilinear anisotropic MMS — angular redistribution probe":

- `sn-mms-spherical-aniso-psi`
- `sn-mms-spherical-aniso-qext`
- `sn-mms-cylindrical-aniso-psi`
- `sn-mms-cylindrical-aniso-qext`

V&V matrix auto-regenerated; each label shows 1 test (foundation-tagged via `@pytest.mark.verifies`).

**12 foundation tests** at `tests/derivations/test_sn_mms_anisotropic_symbolic.py`:

- 4 spherical: substitution identity, redistribution non-vanish, BC vanish, overall pass.
- 4 cylindrical: same structure.
- 2 L1 cross-check (sphere, cylinder) pinning Branch1↔Branch2 at 1e-13 atol.
- 2 quadrature-moment (sphere, cylinder) pinning φ_discrete = A.

Total runtime: 3.6 s. Sphinx -W clean for the new section (pre-existing warnings on `paramref`/`/skills/...` not introduced).

**Atomic commits** (2):
1. `feat(mms)`: Branch 1 SymPy + Branch 2 numpy + 12 foundation tests (1087 lines added).
2. `docs(theory)`: §"Curvilinear anisotropic MMS" + 4 equation labels (172 lines added; matrix.rst auto-regenerated +10 lines).

**Out of scope (explicit)**: consumer L1 convergence test (`tests/sn/_l1/test_mms_spherical_anisotropic_dd_convergence_O_h2.py`) is Phase-0 main-agent work; not modifying the existing isotropic factories; not touching the SN reshape (Issues 1-18 of `.claude/plans/sn_reshape.md`); not 2D Cartesian / P1 aniso MMS (those exist).

**Frame finding (cross-domain-frames)**: the `(1-μ²)B/r` ↔ `ξ²B/r` mapping is structurally a **connection coefficient** (Christoffel-symbol analog for the angular phase-space) — both terms come from how the radial direction cosine rotates with position in the curved geometry. The candidate frame `candidate_cylindrical_connection.md` in `cross-domain-frames` SKILL is now empirically supported by this mapping; the spherical `(1-μ²)/r ∂_μ` and cylindrical `-(1/r)∂_φ(ξ·)` are the same operator viewed in two coordinate charts on SO(3). Worth elevating from candidate to validated if the user wants to pull the abstraction up to a shared "angular-redistribution kernel" primitive.

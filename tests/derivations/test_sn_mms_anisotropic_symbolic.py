r"""Foundation-tagged symbolic tests for the anisotropic curvilinear
MMS sources (Phase 3.6 — defends against vv-principles failure mode #7).

The Branch-1 SymPy module (in
:mod:`orpheus.derivations.continuous.mms.sn`) substitutes the chosen
anisotropic ansatz

.. math::

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + B(r)\,\zeta_n\bigr),

(:math:`\zeta = \mu` for sphere, :math:`\zeta = \eta` for cylinder)
into the continuous transport operator and verifies algebraically
that the closed-form :math:`Q^{\rm ext}_n(r)` returned by the
factory's ``external_source`` method is the unique source that
satisfies the equation.

A second tier of tests pins the **L1 cross-check gate**: the
numerical source built by Branch 2 (the dataclass's ``external_source``
method, vectorised numpy) must agree with the SymPy-evaluated
:math:`Q^{\rm ext}_n(r)` at cell centres on a sample mesh to floating-
point precision (~ :math:`10^{-13}`). Disagreement points to a
copying error between the symbolic and numerical branches.

Why this matters
----------------

The existing isotropic MMS tests (in :file:`tests/sn/test_mms_curvilinear.py`)
use ansatz :math:`\psi_n(r) = A(r)/W`. The angular-redistribution
operator — the curvilinear sweep's hardest math, where ERR-026 lives —
gives **identically zero** under that ansatz. The anisotropic ansatz
defended here makes the redistribution term non-trivially active so
the L1 convergence test can mathematically observe a redistribution
bug class.

References
----------

- ``.claude/skills/vv-principles/SKILL.md`` (failure mode #7).
- :doc:`/theory/methods/sn/verification` — anisotropic curvilinear MMS
  section (labels ``sn-mms-spherical-aniso-psi``, ``sn-mms-spherical-aniso-qext``,
  ``sn-mms-cylindrical-aniso-psi``, ``sn-mms-cylindrical-aniso-qext``).
- [BaileyMorelChang2010]_ for the spherical and cylindrical
  angular redistribution operator structure.
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.continuous.mms.sn import (
    SNCylindricalAnisotropicMMSCase,
    SNSphericalAnisotropicMMSCase,
    _cylindrical_anisotropic_symbolic,
    _spherical_anisotropic_symbolic,
    build_cylindrical_anisotropic_mms_case,
    build_spherical_anisotropic_mms_case,
    derive_cylindrical_anisotropic_mms,
    derive_spherical_anisotropic_mms,
)


# ═══════════════════════════════════════════════════════════════════════
# V_sph-aniso — spherical anisotropic MMS source identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-spherical-aniso-psi",
    "sn-mms-spherical-aniso-qext",
)
def test_v_sph_aniso_substitution_identity():
    r"""V_sph-aniso — substituting the ansatz into the spherical SN
    operator yields the closed-form :math:`Q^{\rm ext}_n` exactly
    (zero residual under :func:`sympy.simplify`).

    This is the load-bearing claim: the closed form is the unique
    manufactured source consistent with the chosen ansatz, the
    spherical streaming term, the angular redistribution term, and
    the isotropic-scattering source.
    """
    result = derive_spherical_anisotropic_mms()
    assert result["pass"], (
        f"V_sph-aniso failed: substituted Q differs from closed form by "
        f"{result['diff']}"
    )


@pytest.mark.foundation
def test_v_sph_aniso_redistribution_term_is_nonzero():
    r"""V_sph-aniso — the redistribution-coupled portion of the
    manufactured source :math:`(1-\mu^2)\,B(r)/r` is identically
    non-zero on the open interval :math:`r \in (0, R)` for any
    :math:`\mu \in (-1, 1)`.

    This is what distinguishes the anisotropic from the isotropic
    ansatz (mode #7 defence): the term that fully cancels for the
    isotropic ansatz must NOT cancel here.
    """
    r_sym, mu_sym, R, _, _, _, _, B, _, _, _ = _spherical_anisotropic_symbolic()
    redistribution = (1 - mu_sym**2) * B / r_sym
    # Avoid r=R/2: at r=R/2, cos(π r/R) = cos(π/2) = 0, so B and the
    # redistribution term vanish for trivial cosine-zero reasons unrelated
    # to mode #7. Use r=R/3 where cos(π/3) = 1/2 ≠ 0.
    val_interior = redistribution.subs({r_sym: R / 3, mu_sym: sp.Rational(0)})
    val_simplified = sp.simplify(val_interior)
    assert val_simplified != 0, (
        f"redistribution term unexpectedly cancels at r=R/2, μ=0: "
        f"{val_simplified}"
    )


@pytest.mark.foundation
def test_v_sph_aniso_ansatz_vanishes_at_boundaries():
    r"""V_sph-aniso — every angular flux :math:`\psi_n(r)` vanishes
    at :math:`r=0` (symmetry) and :math:`r=R` (vacuum).

    Both :math:`A(r)` and :math:`B(r)` must vanish at the endpoints
    so the BCs hold uniformly in :math:`\mu_n`. If only :math:`A`
    vanished, the BC would fail for non-zero :math:`\mu_n`.
    """
    r_sym, mu_sym, R, _, _, _, A, B, _, _, _ = _spherical_anisotropic_symbolic()
    assert sp.simplify(A.subs(r_sym, 0)) == 0
    assert sp.simplify(A.subs(r_sym, R)) == 0
    assert sp.simplify(B.subs(r_sym, 0)) == 0
    assert sp.simplify(B.subs(r_sym, R)) == 0


@pytest.mark.foundation
def test_v_sph_aniso_overall_pass():
    r"""V_sph-aniso — final pass marker (mirrors slab/sphere/cylinder
    Variant α convention)."""
    res = derive_spherical_anisotropic_mms()
    assert res["pass"]
    assert "Q_substituted" in res
    assert "Q_closed" in res


# ═══════════════════════════════════════════════════════════════════════
# V_cyl-aniso — cylindrical anisotropic MMS source identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies(
    "transport-cylindrical",
    "sn-mms-cylindrical-aniso-psi",
    "sn-mms-cylindrical-aniso-qext",
)
def test_v_cyl_aniso_substitution_identity():
    r"""V_cyl-aniso — substituting the η-linear ansatz into the
    cylindrical SN operator yields the closed-form :math:`Q^{\rm ext}_n`
    exactly (zero residual under :func:`sympy.simplify`)."""
    result = derive_cylindrical_anisotropic_mms()
    assert result["pass"], (
        f"V_cyl-aniso failed: substituted Q differs from closed form by "
        f"{result['diff']}"
    )


@pytest.mark.foundation
def test_v_cyl_aniso_redistribution_term_is_nonzero():
    r"""V_cyl-aniso — the azimuthal-redistribution term
    :math:`\xi_n^2\,B(r)/r` is identically non-zero on
    :math:`r \in (0, R)` for any :math:`(\theta, \varphi)` with
    :math:`\sin\theta \sin\varphi \neq 0`.

    The cylindrical analog of mode #7 defence: the isotropic ansatz
    forces this term to zero (because :math:`B \equiv 0`), so
    activating it is the entire point of the anisotropic case.
    """
    r_sym, theta_sym, phi_sym, R, _, _, _, _, B, _, xi, _, _, _ = (
        _cylindrical_anisotropic_symbolic()
    )
    term = xi**2 * B / r_sym
    val = term.subs({
        r_sym: R / 3,           # avoid r=R/2 where cos(π r/R) vanishes
        theta_sym: sp.pi / 2,   # equator ordinate (sin θ = 1)
        phi_sym: sp.pi / 4,     # 45° azimuth (sin φ ≠ 0)
    })
    assert sp.simplify(val) != 0, (
        f"cylindrical redistribution unexpectedly cancels: {sp.simplify(val)}"
    )


@pytest.mark.foundation
def test_v_cyl_aniso_ansatz_vanishes_at_boundaries():
    r"""V_cyl-aniso — :math:`\psi_n(r)` vanishes at :math:`r=0`
    and :math:`r=R` for every ordinate."""
    r_sym, _, _, R, _, _, _, A, B, _, _, _, _, _ = (
        _cylindrical_anisotropic_symbolic()
    )
    assert sp.simplify(A.subs(r_sym, 0)) == 0
    assert sp.simplify(A.subs(r_sym, R)) == 0
    assert sp.simplify(B.subs(r_sym, 0)) == 0
    assert sp.simplify(B.subs(r_sym, R)) == 0


@pytest.mark.foundation
def test_v_cyl_aniso_overall_pass():
    res = derive_cylindrical_anisotropic_mms()
    assert res["pass"]


# ═══════════════════════════════════════════════════════════════════════
# L1 cross-check — Branch 2 numerical Q^ext == Branch 1 SymPy-evaluated
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_spherical_aniso_numerical_qext_matches_sympy():
    r"""L1 cross-check (Branch 1 ↔ Branch 2): the Branch-2 numpy
    ``external_source`` agrees with the SymPy-evaluated
    :math:`Q^{\rm ext}_n(r_i)` to machine precision on a sample
    mesh.

    This is the gate that catches a copy-error between the symbolic
    derivation and the numerical implementation. The two branches
    are structurally independent above the trusted-library line:
    Branch 1 evaluates SymPy expressions via ``lambdify``; Branch 2
    uses hand-written numpy. Neither branch shares an in-house
    primitive with the other.
    """
    case: SNSphericalAnisotropicMMSCase = build_spherical_anisotropic_mms_case(
        sigma_t=1.0, sigma_s=0.5, radius=5.0, n_ordinates=8,
    )
    n_cells = 12
    mesh = case.build_mesh(n_cells=n_cells)
    Q_numerical = case.external_source(mesh)         # (N, ng=1, nx)

    # Branch 1: lambdify the SymPy closed form and evaluate at the
    # same (mu_n, r_i) grid as Branch 2.
    (r_sym, mu_sym, R_sym, Sigma_t_sym, Sigma_s_sym, _,
     _A, _B, _psi, _phi, Q_closed) = _spherical_anisotropic_symbolic()
    Q_lam = sp.lambdify(
        (r_sym, mu_sym, R_sym, Sigma_t_sym, Sigma_s_sym),
        Q_closed,
        modules="numpy",
    )
    r_centers = mesh.centers
    mu_n = case.quadrature.mu_x

    # Evaluate on the (N, nx) grid. r_centers is (nx,), mu_n is (N,).
    rr = r_centers[None, :].repeat(len(mu_n), axis=0)
    mm = mu_n[:, None].repeat(len(r_centers), axis=1)
    Q_sympy_grid = Q_lam(
        rr, mm, case.radius, case.sigma_t, case.sigma_s,
    )    # (N, nx)

    # R-1 Step 4 A1 — Branch 2's ``external_source`` emits per-ordinate
    # **density** (the raw SymPy closed form divided by :math:`\sum_n w_n`
    # at the producer boundary, Pattern 7).  The mathematical equivalence
    # check is therefore against the SymPy form scaled by ``1/sum_w``.
    sum_w = float(case.quadrature.weights.sum())
    np.testing.assert_allclose(
        Q_numerical[:, 0, :], Q_sympy_grid / sum_w,
        atol=1e-13, rtol=1e-13,
    )


@pytest.mark.foundation
def test_cylindrical_aniso_numerical_qext_matches_sympy():
    r"""L1 cross-check (Branch 1 ↔ Branch 2) for the cylindrical case.

    Branch 1 uses the SymPy expression in :math:`(r, \theta, \varphi)`;
    Branch 2 uses :math:`(\eta, \xi)` extracted from
    :class:`ProductQuadrature`. The agreement is the strongest check
    that the :math:`(\eta, \xi) \leftrightarrow (\theta, \varphi)`
    convention is consistent between the two branches.
    """
    case: SNCylindricalAnisotropicMMSCase = build_cylindrical_anisotropic_mms_case(
        sigma_t=1.0, sigma_s=0.5, radius=5.0, n_mu=4, n_phi=8,
    )
    n_cells = 10
    mesh = case.build_mesh(n_cells=n_cells)
    Q_numerical = case.external_source(mesh)

    (r_sym, theta_sym, phi_sym, R_sym, Sigma_t_sym, Sigma_s_sym,
     _, _A, _B, eta_sym, xi_sym, _psi, _phi, Q_closed) = (
        _cylindrical_anisotropic_symbolic()
    )
    Q_lam = sp.lambdify(
        (r_sym, theta_sym, phi_sym, R_sym, Sigma_t_sym, Sigma_s_sym),
        Q_closed,
        modules="numpy",
    )

    # ProductQuadrature gives (eta, xi) = (η_n, ξ_n) directly. We
    # invert to (theta_n, phi_n) on the unit sphere: η = sinθ cosφ,
    # ξ = sinθ sinφ → tan φ = ξ/η, sinθ = sqrt(η² + ξ²).
    eta = case.quadrature.eta
    xi = case.quadrature.xi
    sin_theta = np.sqrt(eta**2 + xi**2)
    # phi_az: arctan2 gives the right quadrant in [-π, π].
    phi_az = np.arctan2(xi, eta)
    theta = np.arcsin(np.clip(sin_theta, -1.0, 1.0))

    r_centers = mesh.centers
    rr = r_centers[None, :].repeat(len(eta), axis=0)
    tt = theta[:, None].repeat(len(r_centers), axis=1)
    pp = phi_az[:, None].repeat(len(r_centers), axis=1)

    Q_sympy_grid = Q_lam(
        rr, tt, pp, case.radius, case.sigma_t, case.sigma_s,
    )    # (N, nx)

    # R-1 Step 4 A1 — Branch 2 emits per-ord density (raw / sum_w).
    sum_w = float(case.quadrature.weights.sum())
    np.testing.assert_allclose(
        Q_numerical[:, 0, :], Q_sympy_grid / sum_w,
        atol=1e-13, rtol=1e-13,
    )


# ═══════════════════════════════════════════════════════════════════════
# Sanity invariants — vanishing scalar-flux moments + BC consistency
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_spherical_aniso_phi_equals_A_under_quadrature():
    r"""For symmetric Gauss-Legendre the discrete scalar flux
    :math:`\sum_n w_n \psi_n / W` reduces to :math:`A(r)` because
    :math:`\sum_n w_n \mu_n = 0`. Confirms ``phi_exact`` is the
    correct reference flux for solver convergence tests."""
    case = build_spherical_anisotropic_mms_case(n_ordinates=16)
    r_sample = np.linspace(0.5, 4.5, 7)
    psi = (
        case.A(r_sample)[None, :]
        + case.B(r_sample)[None, :] * case.quadrature.mu_x[:, None]
    )                                              # (N, nx) (without 1/W)
    weights = case.quadrature.weights
    W = float(weights.sum())
    phi_discrete = (weights[:, None] * psi).sum(axis=0) / W
    np.testing.assert_allclose(phi_discrete, case.phi_exact(r_sample), atol=1e-14)


@pytest.mark.foundation
def test_cylindrical_aniso_phi_equals_A_under_quadrature():
    r"""ProductQuadrature analog: :math:`\sum_n w_n \eta_n = 0` so
    :math:`\phi(r) = A(r)`."""
    case = build_cylindrical_anisotropic_mms_case(n_mu=4, n_phi=8)
    r_sample = np.linspace(0.5, 4.5, 7)
    psi = (
        case.A(r_sample)[None, :]
        + case.B(r_sample)[None, :] * case.quadrature.mu_x[:, None]
    )
    weights = case.quadrature.weights
    W = float(weights.sum())
    phi_discrete = (weights[:, None] * psi).sum(axis=0) / W
    np.testing.assert_allclose(phi_discrete, case.phi_exact(r_sample), atol=1e-14)

r"""Foundation-tagged symbolic tests for the Phase 4 / O.2b 4.6
NON-VACUUM prescribed-inflow MMS sources (slab + sphere).

The Branch-1 SymPy module (in
:mod:`orpheus.derivations.continuous.mms.sn`) substitutes the
non-vacuum P1 ansatz

.. math::

   \psi_n(x) = \frac{1}{W}\bigl(A(x) + \mu_n B(x)\bigr)

— with :math:`A` chosen NON-zero at the boundary (``a0>0``) so the
prescribed-inflow trace :math:`\gamma_-\psi \neq 0` — into the
continuous transport operator and verifies algebraically that the
closed-form :math:`Q^{\rm ext}_n` returned by the factory's
``external_source`` is the unique source that satisfies the equation.

The single structural delta over the existing (vacuum) MMS catalog is
the non-vanishing boundary trace; the angular form is the proven P1
element (see ``tests/derivations/test_sn_mms_anisotropic_symbolic.py``).

Tiers:

1. **Substitution identity** (slab + sphere): SymPy ``simplify`` of the
   substituted operator minus the closed form is exactly zero.
2. **Regression pin (decision A)**: the EXISTING
   :func:`derive_spherical_anisotropic_mms` still passes after the
   parameterization of :func:`_spherical_anisotropic_symbolic` — the
   no-arg default reproduces the Phase 3.6 vacuum shapes byte-for-byte.
3. **Boundary-trace non-vanishing** (the 4.6 novelty): :math:`A` (and
   hence :math:`\gamma_-\psi`) is NON-zero at the outer face; HAZARD H1
   :math:`B(0)=0` for the sphere (pole regularity), :math:`B(R)\neq 0`.
4. **L1 cross-check (Branch 1 ↔ Branch 2)**: the slab numpy
   ``external_source`` is bit-equal to the lambdified SymPy closed form
   on a sample mesh (~:math:`10^{-13}`).

See:
- ``.claude/skills/vv-principles/SKILL.md`` (failure mode #7).
- ``.claude/skills/algebra-of-record/SKILL.md`` (Branch 1 / Branch 2).
- ``docs/theory/methods/sn/verification.rst`` labels
  ``sn-mms-nonvacuum-psi``/``-qext`` (slab),
  ``sn-mms-nonvacuum-sph-psi``/``-qext`` (sphere).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.continuous.mms.sn import (
    _nonvacuum_slab_symbolic,
    _nonvacuum_spherical_AB,
    _spherical_anisotropic_symbolic,
    build_slab_nonvacuum_mms_case,
    build_sphere_nonvacuum_mms_case,
    derive_nonvacuum_slab_mms,
    derive_nonvacuum_spherical_mms,
    derive_spherical_anisotropic_mms,
)


# ═══════════════════════════════════════════════════════════════════════
# V_nonvac-slab — non-vacuum slab MMS source identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies(
    "transport-cartesian",
    "sn-mms-nonvacuum-psi",
    "sn-mms-nonvacuum-qext",
)
def test_v_nonvac_slab_substitution_identity():
    r"""V_nonvac-slab — substituting the non-vacuum ansatz into the slab
    SN operator yields the closed-form :math:`Q^{\rm ext}_n` exactly
    (zero residual under :func:`sympy.simplify`)."""
    result = derive_nonvacuum_slab_mms()
    assert result["pass"], (
        f"V_nonvac-slab failed: substituted Q differs from closed form "
        f"by {result['diff']}"
    )


@pytest.mark.foundation
def test_v_nonvac_slab_ansatz_nonvanishing_at_faces():
    r"""V_nonvac-slab — :math:`A` (and hence :math:`\gamma_-\psi`) is
    NON-zero at both faces. THIS is the 4.6 novelty: the existing slab
    MMS vanishes at the faces; this one does not, lighting up the
    prescribed-inflow ``q.boundary`` path."""
    x, mu, k, a0, a1, b0, _, _, _, A, B, _, _, _ = _nonvacuum_slab_symbolic()
    L = sp.Symbol("L", positive=True)
    # A(0) = a0 (must be > 0). a0 is a free symbol here; the factory
    # bakes a0 = 1/2 > 0. The symbolic claim is A(0) = a0 (not 0).
    assert sp.simplify(A.subs(x, 0) - a0) == 0
    # B(0) = b0 != 0 (slab has no pole — non-vanishing B is fine).
    assert sp.simplify(B.subs(x, 0) - b0) == 0
    # The shipped factory coefficients are non-vacuum:
    case = build_slab_nonvacuum_mms_case()
    assert case.a0 > 0.0
    assert float(case.A(0.0)) > 0.0
    assert float(case.A(case.slab_length)) != 0.0


# ═══════════════════════════════════════════════════════════════════════
# V_nonvac-sph — non-vacuum spherical MMS source identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-nonvacuum-sph-psi",
    "sn-mms-nonvacuum-sph-qext",
)
def test_v_nonvac_sph_substitution_identity():
    r"""V_nonvac-sph — substituting the non-vacuum ansatz into the
    spherical SN operator (REUSING :func:`_spherical_anisotropic_symbolic`
    with the 4.6 shapes) yields the closed-form :math:`Q^{\rm ext}_n`
    exactly. The redistribution term :math:`(1-\mu^2)B/r` is live."""
    result = derive_nonvacuum_spherical_mms()
    assert result["pass"], (
        f"V_nonvac-sph failed: substituted Q differs from closed form "
        f"by {result['diff']}"
    )


@pytest.mark.foundation
def test_v_nonvac_sph_pole_regularity_and_nonvacuum():
    r"""V_nonvac-sph — HAZARD H1: :math:`B(0)=0` (pole regularity — the
    redistribution :math:`(1-\mu^2)B/r` is finite at :math:`r=0`), AND
    the boundary trace is NON-vacuum: :math:`A(R)\neq 0`, :math:`B(R)\neq
    0`, :math:`A(0)=a_0>0` finite."""
    A, B = _nonvacuum_spherical_AB()
    r, R, k = sp.symbols("r R k", positive=True, real=True)
    # B(0) = 0 (the (r/R) prefactor) — pole regularity.
    assert sp.simplify(B.subs(r, 0)) == 0
    # A(0) = a0 = 1/2 > 0 (finite at pole; A has no 1/r companion).
    assert sp.simplify(A.subs(r, 0)) == sp.Rational(1, 2)
    # B(R) != 0 for generic k (non-vacuum outer inflow).
    B_at_R = sp.simplify(B.subs(r, R))
    assert B_at_R != 0
    # Concretely, at k*R = pi/2: B(R) = 3/10, A(R) = 3/4 — both non-zero.
    subs_kR = {k: sp.pi / (2 * R)}
    assert sp.simplify(B.subs(r, R).subs(subs_kR)) == sp.Rational(3, 10)
    assert sp.simplify(A.subs(r, R).subs(subs_kR)) == sp.Rational(3, 4)


@pytest.mark.foundation
def test_v_nonvac_sph_redistribution_term_is_nonzero():
    r"""V_nonvac-sph — the redistribution-coupled portion
    :math:`(1-\mu^2)B(r)/r` is non-zero on :math:`(0, R)` for
    :math:`\mu \neq \pm 1`. This is the ERR-026 path the slab nulls;
    the non-vacuum sphere must activate it (Mode-7 companion)."""
    A, B = _nonvacuum_spherical_AB()
    r, R, k = sp.symbols("r R k", positive=True, real=True)
    mu = sp.Symbol("mu", real=True)
    redistribution = (1 - mu**2) * B / r
    # Evaluate at an interior point, μ=0, generic k: must be non-zero.
    val = redistribution.subs({r: R / 3, mu: 0})
    # B(R/3)/(R/3) = (1/3)(b0 + b1 cos(kR/3)) / (R/3) = (b0+b1 cos)/R
    assert sp.simplify(val) != 0


# ═══════════════════════════════════════════════════════════════════════
# Regression pin (decision A) — the EXISTING vacuum derivation still passes
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_existing_spherical_aniso_still_passes_after_parameterization():
    r"""Decision A regression pin — parameterizing
    :func:`_spherical_anisotropic_symbolic` with optional ``A=None,
    B=None`` MUST leave the existing no-arg vacuum derivation
    byte-unchanged. The default branch reproduces the Phase 3.6 shapes
    :math:`A=\sin(\pi r/R)`, :math:`B=(r/R)(1-r/R)\cos(\pi r/R)`."""
    res = derive_spherical_anisotropic_mms()
    assert res["pass"], f"existing sphere-aniso regressed: {res['diff']}"

    # And the no-arg default reproduces the Phase 3.6 shapes exactly.
    r, mu, R, _, _, _, A, B, _, _, _ = _spherical_anisotropic_symbolic()
    assert sp.simplify(A - sp.sin(sp.pi * r / R)) == 0
    assert sp.simplify(B - (r / R) * (1 - r / R) * sp.cos(sp.pi * r / R)) == 0


# ═══════════════════════════════════════════════════════════════════════
# L1 cross-check — Branch 2 numerical Q^ext == Branch 1 SymPy-evaluated
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_slab_nonvacuum_numerical_qext_matches_sympy():
    r"""L1 cross-check (Branch 1 ↔ Branch 2): the slab numpy
    ``external_source`` agrees with the SymPy closed form to machine
    precision on a sample mesh. Catches a copy-error between the
    symbolic derivation and the numerical implementation. The two
    branches are structurally independent above the trusted-library line
    (Branch 1 = ``lambdify``; Branch 2 = hand-written numpy)."""
    case = build_slab_nonvacuum_mms_case(
        sigma_t=1.0, sigma_s=0.5, slab_length=5.0, n_ordinates=8,
    )
    mesh = case.build_mesh(n_cells=12)
    Q_numerical = case.external_source(mesh)         # (N, ng=1, nx)

    # Branch 1: lambdify the SymPy closed form over (x, mu) + parameters.
    # The symbolic Q_closed uses Σ_t, Σ_s directly (1g, c=1); the case's
    # 1g amplitude c0=1, so A_0(x) = a0 + a1 sin(kx) matches the symbolic A.
    (x_s, mu_s, k_s, a0_s, a1_s, b0_s, St_s, Ss_s, _W,
     _A, _B, _psi, _phi, Q_closed) = _nonvacuum_slab_symbolic()
    Q_lam = sp.lambdify(
        (x_s, mu_s, k_s, a0_s, a1_s, b0_s, St_s, Ss_s),
        Q_closed, modules="numpy",
    )
    x_centers = mesh.centers
    mu_n = case.quadrature.mu_x
    xx = x_centers[None, :].repeat(len(mu_n), axis=0)
    mm = mu_n[:, None].repeat(len(x_centers), axis=1)
    Q_sympy_grid = Q_lam(
        xx, mm, case.k, case.a0, case.a1, case.b0,
        float(case.sigma_t_g[0]), float(case.sigma_s_matrix[0, 0]),
    )    # (N, nx)

    # Branch 2 emits per-ord density (raw SymPy / sum_w at the producer).
    # external_source layout is (N, ng, nx); slice g=0 → (N, nx).
    sum_w = float(case.quadrature.weights.sum())
    np.testing.assert_allclose(
        Q_numerical[:, 0, :], Q_sympy_grid / sum_w,
        atol=1e-13, rtol=1e-13,
    )


@pytest.mark.foundation
def test_sphere_nonvacuum_numerical_qext_matches_sympy():
    r"""L1 cross-check (Branch 1 ↔ Branch 2) for the sphere. The
    Branch-2 numpy ``external_source`` agrees with the lambdified SymPy
    closed form (the SAME a0,a1,b0,b1 and ``kR=π/2`` baked into
    :func:`_nonvacuum_spherical_AB`) to machine precision. Confirms the
    :math:`(r/R)[b_0+b_1\cos]` :math:`B(r)` is transcribed identically in
    both branches (the redistribution term :math:`(1-\mu^2)B/r` is the
    delicate one)."""
    case = build_sphere_nonvacuum_mms_case(
        sigma_t=1.0, sigma_s=0.5, radius=5.0, n_ordinates=8,
    )
    mesh = case.build_mesh(n_cells=12)
    Q_numerical = case.external_source(mesh)         # (N, 1, nx, 1)

    A_nv, B_nv = _nonvacuum_spherical_AB()
    r_s, mu_s, R_s, St_s, Ss_s, _W, _A, _B, _psi, _phi, Q_closed = (
        _spherical_anisotropic_symbolic(A=A_nv, B=B_nv)
    )
    k_s = sp.Symbol("k", positive=True, real=True)
    Q_lam = sp.lambdify(
        (r_s, mu_s, R_s, St_s, Ss_s, k_s), Q_closed, modules="numpy",
    )
    r_centers = mesh.centers
    mu_n = case.quadrature.mu_x
    rr = r_centers[None, :].repeat(len(mu_n), axis=0)
    mm = mu_n[:, None].repeat(len(r_centers), axis=1)
    Q_sympy_grid = Q_lam(
        rr, mm, case.radius, case.sigma_t, case.sigma_s, case.k,
    )    # (N, nx)

    sum_w = float(case.quadrature.weights.sum())
    np.testing.assert_allclose(
        Q_numerical[:, 0, :], Q_sympy_grid / sum_w,
        atol=1e-13, rtol=1e-13,
    )


@pytest.mark.foundation
def test_slab_nonvacuum_phi_equals_A_under_quadrature():
    r"""Discrete :math:`\sum_n w_n \psi_n / W = A(x)` (since
    :math:`\sum_n w_n \mu_n = 0`). Confirms ``phi_exact`` is the right
    reference flux for the convergence rows."""
    case = build_slab_nonvacuum_mms_case(n_ordinates=16)
    x_sample = np.linspace(0.5, 4.5, 7)
    psi = (
        case.A(x_sample)[None, :]
        + case.B(x_sample)[None, :] * case.quadrature.mu_x[:, None]
    )                                              # (N, nx) (without 1/W)
    weights = case.quadrature.weights
    W = float(weights.sum())
    phi_discrete = (weights[:, None] * psi).sum(axis=0) / W
    np.testing.assert_allclose(
        phi_discrete, case.phi_exact(x_sample), atol=1e-14,
    )

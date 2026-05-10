r"""L1 flat-flux algebraic identity for the PoleAngularClosure family.

This test pins the **structural invariant** that the three Phase B
strategies preserve on **angularly-flat** :math:`\psi`: the per-ordinate
redistribution outputs collapse to the same algebraic form, regardless
of which strategy is used.

The Bailey :math:`\Delta A/w` flat-flux invariant (the per-ordinate
streaming + redistribution = 0 cancellation that
:class:`~orpheus.sn.spatial.pole_angular_closure.BaileyFlatFluxRedist`
explicitly enforces) is the special-case algebraic identity all three
strategies satisfy on flat :math:`\psi` — for both the canonical
M-M weighted DD recurrence
(:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`)
and the pre-Phase-B :math:`\tau`-symmetric form
(:class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`),
the algebra collapses to the same answer when the input is angle-flat.

Why this is L1 (not foundation)
-------------------------------

This is **not** a software-contract claim — it is a transport-equation
identity, the discrete analogue of

.. math::

   \int_{-1}^{1} d\mu\, \partial_\mu\bigl[(1-\mu^2)\,\psi(r,\mu)\bigr]
   = 0 \quad \text{for } \psi \text{ constant in } \mu.

The three strategies discretize this integral with different choices
of :math:`\psi_{n\pm 1/2}` evaluator; on flat :math:`\psi` all three
choices reduce to the cell-centre value :math:`\psi_n`, so the
discretized integral is the same algebraic expression
:math:`(\Delta A/w)\,(\alpha_{n+1/2} - \alpha_{n-1/2})\,\psi_n / V_i`
which by the :math:`\alpha`-recurrence equals
:math:`-\mu_n\,\Delta A_i\,\psi_n / V_i` — the per-ordinate streaming
cancellation Bailey 2009 enforces.

The test is the **consistency check** that the Phase B Protocol is a
correct generalization of the Phase A operator, not a different
operator: the family agrees on the pre-Phase-B regime (flat ψ) and
diverges only on the regime the pre-Phase-B form was inaccurate for
(angularly-varying ψ).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.reduced_operator import (
    cylindrical_streaming,
    spherical_streaming,
)
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.spatial.pole_angular_closure import (
    BaileyFlatFluxRedist,
    LegacyTauSymmetricInterpolation,
    MorelMontryAngularSweep,
)

pytestmark = pytest.mark.l1


# ═══════════════════════════════════════════════════════════════════════
# Spherical flat-flux identity
# ═══════════════════════════════════════════════════════════════════════


def _build_spherical_inputs(nx: int = 5, N: int = 8, ng: int = 2):
    """Build the Phase B closure inputs for a sphere mesh."""
    quad = GaussLegendre1D.create(n_ordinates=N)
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    red = spherical_streaming(mesh, quad)
    V = mesh.volumes
    return quad, red, V, ng


def test_spherical_flat_flux_legacy_matches_bailey_collapse():
    r"""On angularly-flat :math:`\psi`, the
    :class:`LegacyTauSymmetricInterpolation` strategy reduces to the
    :class:`BaileyFlatFluxRedist` algebraic identity bit-for-bit.

    Both strategies' :math:`\psi_{n\pm 1/2}` evaluations collapse to
    the cell-centre value when :math:`\psi` is angle-flat
    (:math:`\tau\,\psi + (1-\tau)\,\psi = \psi`), so the
    redistribution is the algebraic identity
    :math:`(\Delta A/w)\,(\alpha_{n+1/2}-\alpha_{n-1/2})\,\psi/V`.

    This is the **load-bearing** consistency check on the Phase B
    default — :class:`LegacyTauSymmetricInterpolation` IS the Phase B
    default, and it MUST agree with :class:`BaileyFlatFluxRedist` on
    flat ψ for the per-ordinate flat-flux invariant to hold (which
    Phase A's flat-flux test in
    :file:`tests/sn/test_snstreamingoperator.py` and the ERR-026
    evidence test in
    :file:`tests/sn/test_sweep_operator_inconsistency.py` depend on).
    """
    quad, red, V, ng = _build_spherical_inputs()
    N, nx = quad.N, len(V)
    psi_flat = 3.7 * np.ones((ng, N, nx))

    bff_result = BaileyFlatFluxRedist()(
        psi_flat, red.alpha_half, red.redist_dAw, red.tau_mm, V,
    )
    legacy_result = LegacyTauSymmetricInterpolation()(
        psi_flat, red.alpha_half, red.redist_dAw, red.tau_mm, V,
    )
    np.testing.assert_allclose(
        bff_result, legacy_result, rtol=1e-13,
        err_msg="LegacyTauSymmetricInterpolation should reduce to the "
                "Bailey flat-flux algebraic identity on angularly-flat "
                "ψ.  Disagreement here would mean the legacy form is "
                "inconsistent with the Bailey collapse on flat ψ.",
    )


def test_spherical_flat_flux_morel_montry_integrated_identity():
    r"""On angularly-flat :math:`\psi`, the canonical Hébert
    :class:`MorelMontryAngularSweep` recurrence does **NOT** collapse
    to the Bailey form per-ordinate (the half-angle face fluxes
    oscillate :math:`0, 2c, 0, 2c, \ldots` under the DD recurrence
    seeded at :math:`\phi_{1/2,i} = 0`), but the **angular-integrated**
    invariant

    .. math::

       \sum_n \mathcal{W}_n \cdot R_{n,i,g}^{\rm MMS}
       = \sum_n \mathcal{W}_n \cdot R_{n,i,g}^{\rm BFF}
       = 0 \quad \text{(by α-telescoping)}

    holds bit-for-bit.  This is the consistency check: the MMS
    canonical form is a correct generalization of the BFF form on
    the angular-integrated layer (the layer that physics cares about),
    but a different operator at per-ordinate granularity.
    """
    quad, red, V, ng = _build_spherical_inputs()
    N, nx = quad.N, len(V)
    psi_flat = 3.7 * np.ones((ng, N, nx))

    mms_result = MorelMontryAngularSweep()(
        psi_flat, red.alpha_half, red.redist_dAw, red.tau_mm, V,
    )
    bff_result = BaileyFlatFluxRedist()(
        psi_flat, red.alpha_half, red.redist_dAw, red.tau_mm, V,
    )
    weights = quad.weights
    mms_integrated = np.sum(
        weights[None, :, None] * mms_result, axis=1,
    )
    bff_integrated = np.sum(
        weights[None, :, None] * bff_result, axis=1,
    )
    # Both should be zero (α-telescoping).
    np.testing.assert_allclose(
        mms_integrated, np.zeros_like(mms_integrated), atol=1e-12,
        err_msg="MMS canonical form should preserve the α-telescoping "
                "angular-integrated flat-flux invariant.",
    )
    np.testing.assert_allclose(
        bff_integrated, np.zeros_like(bff_integrated), atol=1e-12,
        err_msg="BFF should preserve the α-telescoping angular-"
                "integrated flat-flux invariant.",
    )
    # And they should agree with each other on the integrated layer.
    np.testing.assert_allclose(
        mms_integrated, bff_integrated, atol=1e-12,
        err_msg="MMS and BFF should agree on the angular-integrated "
                "flat-flux invariant.",
    )


def test_spherical_flat_flux_per_ordinate_streaming_cancellation():
    r"""On flat :math:`\psi`, the BFF redistribution per ordinate
    equals :math:`-\mu_n\,\Delta A_i\,\psi/V_i` — the algebraic form
    that cancels the streaming term per-ordinate (the Bailey
    :math:`\Delta A/w` invariant).
    """
    quad, red, V, ng = _build_spherical_inputs(nx=5, N=8, ng=1)
    N, nx = quad.N, len(V)
    psi_const = 1.0
    psi_flat = psi_const * np.ones((ng, N, nx))

    redist = BaileyFlatFluxRedist()(
        psi_flat, red.alpha_half, red.redist_dAw, red.tau_mm, V,
    )
    # Build the analytic expectation: -μ_n · ΔA_i · ψ / V_i.
    expected = np.zeros((ng, N, nx))
    delta_A = red.delta_A
    for n in range(N):
        for i in range(nx):
            expected[0, n, i] = -quad.mu_x[n] * delta_A[i] * psi_const / V[i]
    np.testing.assert_allclose(redist, expected, rtol=1e-13)


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical flat-flux identity (per-level)
# ═══════════════════════════════════════════════════════════════════════


def _build_cylindrical_inputs(nx: int = 5, ng: int = 2):
    """Build Phase B closure inputs for a cylinder mesh + level structure."""
    quad = ProductQuadrature.create(n_mu=4, n_phi=8)
    mesh = Mesh1D(
        edges=np.linspace(0.01, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    red = cylindrical_streaming(mesh, quad)
    V = mesh.volumes
    return quad, red, V, ng


def test_cylindrical_flat_flux_legacy_matches_bailey():
    """Cylindrical analogue of
    :func:`test_spherical_flat_flux_legacy_matches_bailey_collapse`."""
    quad, red, V, ng = _build_cylindrical_inputs()
    N, nx = quad.N, len(V)
    psi_flat = 2.7 * np.ones((ng, N, nx))

    bff = BaileyFlatFluxRedist()(
        psi_flat,
        red.alpha_per_level,
        red.redist_dAw_per_level,
        red.tau_mm_per_level,
        V,
        level_indices=quad.level_indices,
    )
    legacy = LegacyTauSymmetricInterpolation()(
        psi_flat,
        red.alpha_per_level,
        red.redist_dAw_per_level,
        red.tau_mm_per_level,
        V,
        level_indices=quad.level_indices,
    )
    np.testing.assert_allclose(
        bff, legacy, rtol=1e-13,
        err_msg="Cylindrical LegacyTauSymmetricInterpolation should "
                "reduce to the Bailey flat-flux algebraic identity "
                "on flat ψ.",
    )


def test_cylindrical_flat_flux_morel_montry_integrated_identity():
    """Cylindrical analogue of the integrated-identity invariant —
    each :math:`\\mu`-level satisfies its own α-telescoping invariant
    when summed against the level's azimuthal weights."""
    quad, red, V, ng = _build_cylindrical_inputs()
    N, nx = quad.N, len(V)
    psi_flat = 2.7 * np.ones((ng, N, nx))

    mms = MorelMontryAngularSweep()(
        psi_flat,
        red.alpha_per_level,
        red.redist_dAw_per_level,
        red.tau_mm_per_level,
        V,
        level_indices=quad.level_indices,
    )
    bff = BaileyFlatFluxRedist()(
        psi_flat,
        red.alpha_per_level,
        red.redist_dAw_per_level,
        red.tau_mm_per_level,
        V,
        level_indices=quad.level_indices,
    )
    # Per-level weighted sum should vanish for each strategy.
    weights = quad.weights
    for level_idx in quad.level_indices:
        # Slice the per-level weights and contributions.
        w_lvl = weights[level_idx]
        mms_lvl_int = np.sum(
            w_lvl[None, :, None] * mms[:, level_idx, :], axis=1,
        )
        bff_lvl_int = np.sum(
            w_lvl[None, :, None] * bff[:, level_idx, :], axis=1,
        )
        np.testing.assert_allclose(
            mms_lvl_int, np.zeros_like(mms_lvl_int), atol=1e-12,
            err_msg="Cylindrical MMS per-level α-telescoping violated.",
        )
        np.testing.assert_allclose(
            bff_lvl_int, np.zeros_like(bff_lvl_int), atol=1e-12,
            err_msg="Cylindrical BFF per-level α-telescoping violated.",
        )

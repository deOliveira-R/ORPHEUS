"""Foundation + L1 tests for ``orpheus.numerics.quadrature.rules_product``.

The product rule :math:`\\mu_{\\text{GL}} \\times \\phi_{\\text{equispaced}}`
is the cylindrical-SN workhorse. Tests pin shape / metadata, weight
sum, bit-identical match to the legacy
:class:`~orpheus.sn.quadrature.ProductQuadrature` adapter, and
polynomial-exactness in the polar factor.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import SPACE_SPHERE, DiscreteMeasure
from orpheus.numerics.quadrature import product_mu_phi
from orpheus.numerics.quadrature.rules_sphere import LevelStructure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Foundation
# ---------------------------------------------------------------------------


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("n_mu", "n_phi"), [(4, 4), (8, 8), (8, 16), (16, 8), (12, 24)]
)
def test_product_returns_measure_and_structure(n_mu: int, n_phi: int) -> None:
    """Shape, metadata, accompanying :class:`LevelStructure`."""
    m, s = product_mu_phi(n_mu, n_phi)
    assert isinstance(m, DiscreteMeasure)
    assert m.nodes.ndim == 2
    assert m.nodes.shape == (n_mu * n_phi, 3)
    assert m.weights.shape == (n_mu * n_phi,)
    assert m.support == SPACE_SPHERE
    # D_{n_phi h}, not SO(2): the phi grid is a FINITE set, so no
    # continuous rotation preserves it. The rule carried SO2 until
    # 2026-08-02 — a claim no finite point set on S^2 can satisfy.
    assert m.invariance_group == SubgroupOfO3.Dnh(n_phi)
    assert m.degree_of_exactness == min(2 * n_mu - 1, n_phi - 1)

    assert isinstance(s, LevelStructure)
    assert s.n_levels == n_mu
    assert len(s.level_indices) == n_mu


@pytest.mark.foundation
@pytest.mark.parametrize(("n_mu", "n_phi"), [(4, 4), (8, 8), (16, 12)])
def test_product_weight_sum_is_4pi(n_mu: int, n_phi: int) -> None:
    """Weight sum :math:`= 4\\pi`."""
    m, _ = product_mu_phi(n_mu, n_phi)
    assert m.weights.sum() == pytest.approx(4.0 * np.pi, abs=1e-12)


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("n_mu", "n_phi"), [(4, 4), (8, 8), (8, 16), (16, 8), (12, 24)]
)
def test_product_bit_identical_to_legacy_adapter(n_mu: int, n_phi: int) -> None:
    """Bit-identical match against
    :class:`~orpheus.sn.quadrature.ProductQuadrature.create`.

    Pins the cylindrical regression snapshots
    (``cyl_*_product_*.npz``)."""
    m, s = product_mu_phi(n_mu, n_phi)
    legacy = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
    nodes = m.nodes
    assert np.array_equal(legacy.mu_x, nodes[:, 0])
    assert np.array_equal(legacy.mu_y, nodes[:, 1])
    assert np.array_equal(legacy.mu_z, nodes[:, 2])
    assert np.array_equal(legacy.weights, m.weights)
    assert legacy.n_levels == s.n_levels
    assert np.array_equal(legacy.level_mu, s.level_mu)
    for li_legacy, li_new in zip(legacy.level_indices, s.level_indices):
        assert np.array_equal(li_legacy, li_new)


@pytest.mark.foundation
def test_product_rejects_zero() -> None:
    """``n_mu < 1`` or ``n_phi < 1`` must raise."""
    with pytest.raises(ValueError):
        product_mu_phi(0, 8)
    with pytest.raises(ValueError):
        product_mu_phi(8, 0)


# ---------------------------------------------------------------------------
# L1 — polynomial-exactness in the polar factor
# ---------------------------------------------------------------------------


@pytest.mark.l1
@pytest.mark.parametrize("n_mu", [4, 8, 16])
def test_product_polar_polynomial_exactness(n_mu: int) -> None:
    """The polar (Gauss-Legendre) factor of the product rule
    integrates :math:`\\mu^k` exactly for :math:`k \\le 2 n_\\mu - 1`,
    after multiplying out the trivial azimuthal integral
    :math:`\\int_0^{2\\pi} d\\phi = 2\\pi`.

    Closed-form: :math:`\\int_{S^2} \\mu_z^k \\, d\\Omega = 2\\pi \\cdot
    \\int_{-1}^{1} \\mu^k \\, d\\mu` which is :math:`0` for odd ``k``
    and :math:`4\\pi/(k+1)` for even ``k``.
    """
    n_phi = 8  # azimuthal exactness comfortable
    m, _ = product_mu_phi(n_mu, n_phi)
    for k in range(2 * n_mu):
        if k % 2 == 1:
            exact = 0.0
        else:
            exact = 4.0 * np.pi / (k + 1)
        approx = m.integrate(lambda x, k=k: x[:, 2] ** k)
        assert approx == pytest.approx(exact, abs=1e-12), (
            f"product({n_mu},{n_phi}) fails at polar degree {k}: "
            f"{approx} vs {exact}"
        )

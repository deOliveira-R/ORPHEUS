"""Foundation gates for :attr:`Mixture.scattering_ratio`.

The Case–Zweifel ``c = (Σ_s + νΣ_f)/Σ_t`` is the canonical scalar
input to 1G analytical solvers (F_N method, trajectory_resolvent
specular kernels). Promoting it from an ad-hoc helper in
``tests/cross_method/adapters.py`` to a :class:`Mixture` property
eliminates the ``_extract_c`` helper and makes the per-group ratio
available to any consumer of the production XS payload.

These tests verify the formula across the regimes the property is
expected to handle:

* 1G isotropic, super-critical (c > 1) — Sood/LA-13511 column.
* 1G pure-scattering, sub-critical (c < 1).
* 2G multi-group with downscatter — per-group diagnostic.

Foundation tier: this is a software invariant on a derived property,
not a theory-equation claim. No ``verifies(...)`` decorator.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture

pytestmark = [pytest.mark.foundation]


def test_scattering_ratio_one_group_isotropic():
    r"""1G isotropic ``c = (Σ_s + νΣ_f)/Σ_t`` matches hand calc.

    Builds a 1G ``Mixture`` with ``Σ_t = 1.0``, ``Σ_s = 0.4``,
    ``νΣ_f = 0.9`` (super-critical c=1.3 — Sood Ua-1-0 column).
    ``mix.scattering_ratio[0]`` must equal ``1.3`` exactly (the
    formula is closed under double-precision for these values).
    """
    SigT = np.array([1.0])
    SigC = np.array([0.1])
    SigL = np.zeros(1)
    SigF = np.array([0.5])
    SigP = np.array([0.9])  # νσ_f
    SigS = csr_matrix(np.array([[0.4]]))
    # Synthetic abstract XS — no physical energy grid (Phase E).
    mix = Mixture(
        SigC=SigC, SigL=SigL, SigF=SigF, SigP=SigP, SigT=SigT,
        SigS=[SigS], Sig2=csr_matrix((1, 1)),
        chi=np.array([1.0]),
    )
    c = mix.scattering_ratio
    assert c.shape == (1,)
    assert float(c[0]) == 1.3, (
        f"expected c=1.3 exactly; got {float(c[0])!r}"
    )


def test_scattering_ratio_pure_scatter_below_unity():
    r"""Pure scattering (``νΣ_f = 0``) gives ``c = Σ_s/Σ_t < 1``."""
    SigT = np.array([2.0])
    SigC = np.array([0.5])
    SigL = np.zeros(1)
    SigF = np.zeros(1)
    SigP = np.zeros(1)
    SigS = csr_matrix(np.array([[1.5]]))
    # Synthetic abstract XS — no physical energy grid (Phase E).
    mix = Mixture(
        SigC=SigC, SigL=SigL, SigF=SigF, SigP=SigP, SigT=SigT,
        SigS=[SigS], Sig2=csr_matrix((1, 1)),
        chi=np.array([0.0]),
    )
    assert float(mix.scattering_ratio[0]) == 0.75


def test_scattering_ratio_two_group_per_group():
    r"""Multi-group ``scattering_ratio[g]`` is the per-group c, using the
    full P0 row sum (in + out) for ``Σ_{s,g}``.

    2G mixture with diagonal scatter 0.30/0.90 and off-diagonal 0.10
    (downscatter). Group 0 row sum is 0.40, group 1 row sum is 0.90.
    """
    SigC = np.array([0.20, 0.10])
    SigL = np.zeros(2)
    SigF = np.array([0.10, 0.20])
    SigP = np.array([0.25, 0.50])
    SigS_dense = np.array([[0.30, 0.10], [0.00, 0.90]])
    SigS = csr_matrix(SigS_dense)
    SigT = (
        SigC + SigL + SigF + np.array(SigS.sum(axis=1)).ravel()
    )
    # Synthetic abstract XS — no physical energy grid (Phase E).
    mix = Mixture(
        SigC=SigC, SigL=SigL, SigF=SigF, SigP=SigP, SigT=SigT,
        SigS=[SigS], Sig2=csr_matrix((2, 2)),
        chi=np.array([1.0, 0.0]),
    )
    c = mix.scattering_ratio
    expected_g0 = (0.40 + 0.25) / SigT[0]
    expected_g1 = (0.90 + 0.50) / SigT[1]
    np.testing.assert_allclose(c, [expected_g0, expected_g1], rtol=0, atol=0)

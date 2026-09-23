"""Gates for :attr:`Mixture.transport_xs` and :attr:`Mixture.diffusion_coefficient`.

The diffusion data seam (campaign #290, phase P1): the diffusion model
consumes ``D = 1/(3 Σ_tr)`` with ``Σ_tr = Σ_t − rowsum(Σ_s1)`` (the
outflow transport approximation), DERIVED from the ``Mixture`` channels
instead of hand-entered (the legacy ``TwoGroupXS.transport``).

Structure:

* foundation — the arithmetic of the derivation (P1 row sum, the exact
  P0-only limit, the CORE1D.m encoding bridge, the positivity guard).
  Software invariants on derived properties; no ``verifies(...)``.
* L1 — ``D = 1/(3 Σ_tr)`` against the theory equation
  ``:label: diffusion-coefficient`` (``docs/theory/methods/diffusion_1d.rst``),
  hand-computed on a mixture with a genuine P1 moment.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import get_mixture, get_xs


@pytest.mark.foundation
def test_transport_xs_is_total_minus_p1_outscatter_row_sum():
    r"""Region A 2G (carries a genuine P1 moment, :math:`\bar\mu = 0.05`):
    ``transport_xs`` equals the hand-computed ``sig_t − sig_s1.sum(axis=1)``
    bit-identically (same subtraction, same operands)."""
    xs = get_xs("A", "2g")
    mix = get_mixture("A", "2g")
    expected = xs["sig_t"] - xs["sig_s1"].sum(axis=1)
    np.testing.assert_array_equal(mix.transport_xs, expected)
    # The correction is genuinely active: Σ_tr < Σ_t in every group.
    if not np.all(mix.transport_xs < mix.SigT):
        pytest.fail("P1 out-scatter correction did not reduce Σ_t")


@pytest.mark.foundation
def test_transport_xs_p0_only_equals_total_exactly():
    r"""A P0-only mixture (``len(SigS) == 1``) has :math:`\Sigma_{tr} =
    \Sigma_t` EXACTLY — the isotropic limit is the identity, not an
    approximation, so the claim is bit-identity."""
    xs = get_xs("A", "2g")
    mix_p0 = Mixture(
        SigC=xs["sig_c"], SigL=np.zeros(2), SigF=xs["sig_f"],
        SigP=xs["nu"] * xs["sig_f"], SigT=xs["sig_t"],
        SigS=[csr_matrix(xs["sig_s"])], Sig2=[csr_matrix((2, 2))],
        chi=xs["chi"],
    )
    np.testing.assert_array_equal(mix_p0.transport_xs, mix_p0.SigT)


@pytest.mark.l1
@pytest.mark.verifies("diffusion-coefficient")
def test_diffusion_coefficient_matches_definition():
    r"""``D_g = 1/(3 Σ_tr,g)`` (:eq:`diffusion-coefficient`) hand-computed
    on region A 2G, whose P1 moment makes :math:`\Sigma_{tr} \ne \Sigma_t`
    (the gate would be vacuous on P0-only data)."""
    xs = get_xs("A", "2g")
    mix = get_mixture("A", "2g")
    sig_tr = xs["sig_t"] - xs["sig_s1"].sum(axis=1)
    np.testing.assert_array_equal(mix.diffusion_coefficient, 1.0 / (3.0 * sig_tr))
    # Dimensional sanity: D is a positive length [cm].
    if not np.all(mix.diffusion_coefficient > 0.0):
        pytest.fail("diffusion coefficient must be a positive length")


@pytest.mark.foundation
def test_diffusion_coefficient_legacy_core1d_bridge():
    """The bit-identical encoding bridge (#290 user ruling 4, part 1).

    The legacy island's ``TwoGroupXS.transport`` IS the transport XS
    (MATLAB CORE1D.m hands Σ_tr directly). Encoding that library as a
    ``Mixture`` with ``SigT := transport`` and NO P1 moment must
    reproduce the legacy ``D = 1/(3·transport)`` bit-identically, so the
    analytic references do not move. Constants pinned literally (fuel,
    ``orpheus/diffusion/solver.py::_default_xs`` at authoring) so this
    bridge outlives the island's retirement (P6).
    """
    transport = np.array([0.2181, 0.7850])  # CORE1D.m fuel Σ_tr
    mix = Mixture(
        SigC=np.array([0.0096, 0.0959]),    # legacy "absorption" as capture
        SigL=np.zeros(2), SigF=np.zeros(2), SigP=np.zeros(2),
        SigT=transport,
        SigS=[csr_matrix((2, 2))], Sig2=[csr_matrix((2, 2))],
        chi=np.zeros(2),
    )
    np.testing.assert_array_equal(mix.transport_xs, transport)
    np.testing.assert_array_equal(
        mix.diffusion_coefficient, 1.0 / (3.0 * transport)
    )


@pytest.mark.foundation
def test_diffusion_coefficient_rejects_nonpositive_transport():
    """Positive + negative pair for the Σ_tr > 0 guard: a P1 out-scatter
    row sum exceeding Σ_t must raise (a negative D is never silently
    produced); the physical sibling in the other tests must not raise."""
    mix = Mixture(
        SigC=np.array([0.1]), SigL=np.zeros(1), SigF=np.zeros(1),
        SigP=np.zeros(1), SigT=np.array([1.0]),
        SigS=[csr_matrix(np.array([[0.9]])), csr_matrix(np.array([[1.2]]))],
        Sig2=[csr_matrix((1, 1))], chi=np.zeros(1),
    )
    np.testing.assert_array_less(mix.transport_xs, 0.0)
    with pytest.raises(ValueError, match="transport_xs"):
        _ = mix.diffusion_coefficient

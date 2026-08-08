r"""The ξ-mirror pairing controls — and the TOMBSTONE of the #326 defect gates.

WHAT SURVIVES HERE
------------------
1-D cylindrical geometry is invariant under the mirror through the plane
spanned by ``z_hat`` and ``r_hat``: ``(eta, xi, mu_z) -> (eta, -xi, mu_z)``.
The two ``foundation`` rows below pin the PAIRING facts that every consumer
of ``ordinate_permutation`` (the coupled-pole seed, the BC realizer, the
MMS parity gate) relies on: the σ_y pairing is the ξ-mirror involution with
η/μ_z held and weights equal, and on a slab GL rule it is the identity.
These are quadrature-level claims — the full-product rules remain
constructible as rules; only cylindrical SNMesh admission refuses them.

TOMBSTONE — the #326 defect gates (retired at the Q5.6.3 admission flip)
------------------------------------------------------------------------
Until Q5.6.3 this file carried the ISSUE #326 defect record: the
cylindrical Morel–Montry η-march broke the geometry's ξ-mirror symmetry —
per-ordinate, ``max_n |psi_n - psi_{mirror_y(n)}| / max|psi|`` measured
2026-08-01 on isotropic-source fixed-source solves::

    product(n_mu=4, n_phi=8), nx=20, homogeneous    ->  1.19e-1
    product(n_mu=4, n_phi=8), nx=20, heterogeneous  ->  5.14e-2
    level_symmetric(4)                              ->  3.08e-1

MECHANISM (kept because the number alone does not teach): on a FULL-circle
level, two mirror partners share ``eta`` bit-exactly, so the η-midpoint
edge collapses onto the node — the lower partner got ``tau_raw = 1`` and
the upper ``tau_raw = 0``, and the structural clamp turned that into
``{1, 1/2}``: two ordinates the geometry says are identical received
DIFFERENT angular weights.  No re-ordering fixes it (the defect magnitude
was measured ordering-INVARIANT to 1e-16 across three tie-breaks); it is a
CLOSURE defect of the full-circle configuration, and the constructive exit
the original module named — "the half-range level (Hebert §3.9.3), under
which only one member of each pair exists and the symmetry holds by
construction" — is exactly ``Quadrature.folded_product``, which Q5.6.3
made the ONLY admitted cylindrical family.  On the folded arc the η-tie is
UNSPELLABLE (η is injective per level), so the defect's configurations
refuse at SNMesh construction and the three ``xfail(strict=True)`` defect
rows resolved by REFUSAL, not repair.

Retired with their subjects (all built P(4,8)/P(2,n_φ)/LS4 cylinder
meshes, unconstructible since the flip): the three ``_XFAIL_326`` defect
rows and their un-xfailed well-formedness sibling;
``test_mirror_defect_is_an_angular_not_a_spatial_defect``;
``test_mirror_defect_scales_with_the_AZIMUTHAL_order_only``;
``test_cylinder_pole_map_and_axis_crossing_differ_by_exactly_the_xi_mirror``;
``test_tie_break_permutation_does_not_commute_with_the_pole_map``.

SUCCESSOR GATES (the live coverage of the symmetry claim):

* the admission refusals themselves —
  ``tests/sn/mesh/test_cylindrical_quadrature_admission.py`` (the
  full-product and LS families refuse, with the per-family R12a reason);
* the MMS σ_y-parity gate —
  ``tests/sn/verification/mms/test_mms_ordering_blindness.py`` (the folded
  solution's parity verified on the fold's PARENT rule + the bit-exact
  restriction leg);
* the #22/#326 heterogeneous-cylinder FOLDED regression lands at Q5.6.4
  with the honest-τ change (issue #326 closes at Q5.6.5).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import SelfPairedDeck
from orpheus.numerics.quadrature import Quadrature


def _mirror_pairing(quad: Quadrature, axis: str) -> np.ndarray:
    """The ordinate pairing :math:`\\sigma_{axis}` induces — read from
    production's own source (``ordinate_permutation``, the same derivation
    the coupled-pole seed and the BC realizer consume). The control legs
    below verify its ξ-mirror facts INDEPENDENTLY (η/μ_z held, ξ negated,
    weights equal)."""
    pi = quad.ordinate_permutation(
        SelfPairedDeck.mirror(axis=axis, dimension=3).motion
    )
    if pi is None:
        raise AssertionError(f"no {axis}-mirror ordinate pairing on this rule")
    return pi.indices


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [4, 8, 16])
def test_y_mirror_pairing_is_the_xi_mirror_involution(n_phi):
    """The σ_y pairing really pairs ``(eta, xi)`` with ``(eta, -xi)``.

    Quadrature-level control (no SNMesh — full-product rules remain
    constructible as rules).  Consumers: the coupled-pole seed's pairing
    derivation and the MMS parity gate's fold-to-parent lookup.
    """
    quad = Quadrature.product(n_mu=4, n_phi=n_phi)
    partner = _mirror_pairing(quad, "y")
    np.testing.assert_array_equal(partner[partner], np.arange(quad.N))
    np.testing.assert_allclose(quad.eta[partner], quad.eta, atol=1e-15)
    np.testing.assert_allclose(quad.mu_z[partner], quad.mu_z, atol=1e-15)
    np.testing.assert_allclose(quad.xi[partner], -quad.xi, atol=1e-15)
    np.testing.assert_array_equal(quad.weights[partner], quad.weights)


@pytest.mark.foundation
def test_slab_quadrature_has_a_trivial_xi_mirror_control():
    """Control leg: on a 1-D slab GL rule the ``xi`` mirror is the identity.

    Pins that any ξ-mirror structure is specific to azimuthal-bearing
    rules, not an artefact of how the pairing is derived.
    """
    quad = Quadrature.gauss_legendre(8)
    np.testing.assert_array_equal(
        _mirror_pairing(quad, "y"), np.arange(quad.N)
    )

"""MOVED from tests/gates/numerics/test_interval_legendre_basis.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import numpy as np
import pytest
import sympy as sp
from itertools import product

from orpheus.numerics.basis import IntervalLegendreBasis, spatial_moment_basis


pytestmark = [pytest.mark.foundation]


@pytest.mark.rests_on("tests/gates/numerics/test_interval_legendre_basis.py::test_the_endpoint_values_are_plus_minus_one_to_the_degree")


def test_the_ubld_face_terms_are_the_hand_tables() -> None:
    r"""The UBLD assembler's face terms, read off the basis, are the LD hand tables.

    At d = 1, unit width, unit :math:`|\mu|`: the downstream-face matrix is
    ``outer(B(+1), B(+1)) = [[1, 1], [1, 1]]`` and the inflow of a unit
    upstream face is ``B(-1) = [1, -1]`` (Larsen and Morel 1989, the LD face
    traces; typed here). Bit-exact: products of ±1. Since the moment-basis
    carve the assembler reads these from ``endpoint_values``, so this row is
    the production pin; the endpoint row above is the basis-side law.
    """
    from orpheus.transport.spatial._ubld import assemble_inflow_axis, assemble_ubld

    one = [np.ones(())]
    basis = spatial_moment_basis(2, 1)
    np.testing.assert_array_equal(assemble_ubld(one, one, np.zeros(()), basis)["F_out"], np.array([[1.0, 1.0], [1.0, 1.0]]))
    np.testing.assert_array_equal(assemble_inflow_axis(one, one, 0, np.ones(1), basis), np.array([1.0, -1.0]))


@pytest.mark.rests_on("tests/gates/numerics/test_interval_legendre_basis.py::test_the_endpoint_values_are_plus_minus_one_to_the_degree")


def test_the_streaming_table_and_the_face_traces_satisfy_integration_by_parts() -> None:
    r"""``G + Gᵀ = −(B(+1)B(+1)ᵀ − B(−1)B(−1)ᵀ)`` for the hand-tabulated UBLD streaming table ``G``.

    ``_GRAD_1D`` is :math:`-\int_{-1}^{1} \partial_\xi B_i\, B_j \,\mathrm d\xi`
    for :math:`B = (P_0, P_1)`, typed by hand (the interval basis holds no
    derivative pairing); integration by parts ties it to the endpoint values
    the face terms read. Bit-exact: small integers. Tooth: a sign or a slot of
    the hand table changed.
    """
    from orpheus.transport.spatial._ubld import _GRAD_1D

    b = IntervalLegendreBasis(2)
    plus, minus = b.endpoint_values(1), b.endpoint_values(-1)
    np.testing.assert_array_equal(_GRAD_1D + _GRAD_1D.T, -(np.outer(plus, plus) - np.outer(minus, minus)))
    # and one entry against SymPy, so the identity cannot hold for a wrong antisymmetric part
    xi = sp.symbols("xi")
    g10 = -sp.integrate(sp.diff(xi, xi) * 1, (xi, -1, 1))
    assert float(g10) == _GRAD_1D[1, 0]




def _hand_octant_signs(octant: tuple[int, ...]) -> np.ndarray:
    r"""The sign of slot ``[o_0, ..., o_{d-1}]`` is :math:`\prod_a \mathrm{octant}[a]^{o_a}`, decoded by hand, axis 0 outermost."""
    d = len(octant)
    signs = np.ones(2 ** d)
    for flat in range(2 ** d):
        for a in range(d):
            if (flat >> (d - 1 - a)) & 1 and octant[a] < 0:
                signs[flat] = -signs[flat]
    return signs


@pytest.mark.rests_on("tests/gates/numerics/test_interval_legendre_basis.py::test_the_reflection_signs_are_the_mirrors_action")
@pytest.mark.parametrize("d", (1, 2, 3))
def test_the_octant_frame_map_is_the_reflection_signs(d: int) -> None:
    r"""``octant_moment_frame_signs(octant, 2)`` equals the basis's ``reflection_signs`` and the hand Kronecker decode.

    REFERENCE leg: since S3 the octant map reads ``reflection_signs``; it equals
    the hand decode it replaced (:func:`_hand_octant_signs`, axis 0 outermost) for
    all :math:`2^d` octants. Split from numerics'
    ``test_the_reflection_signs_are_the_mirrors_action_and_the_octant_frame_map``
    (layer DAG Phase 2 refresh): the map is transport's. Tooth: the decode with
    axis 0 innermost (d >= 2), ``SITE_ARM=kron_decode_144``.
    """
    from orpheus.transport.spatial._ubld import octant_moment_frame_signs

    basis = spatial_moment_basis(2, d)
    for octant in product((1, -1), repeat=d):
        np.testing.assert_array_equal(octant_moment_frame_signs(octant, 2), _hand_octant_signs(octant))
        np.testing.assert_array_equal(octant_moment_frame_signs(octant, 2), basis.reflection_signs(tuple(s < 0 for s in octant)))

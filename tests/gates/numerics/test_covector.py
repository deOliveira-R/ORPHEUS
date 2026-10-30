r"""The co-vector: a member of the dual V*, minted only as ``V.flat(w) = G w`` and paired canonically with V.

The space-membership unit U2 (plan ``.claude/plans/space_membership_covector_gauge.md``, spec
``scratch/membership_covector/spec/spec.md`` §U2.2). The functional ⟨w, ·⟩_V of a density w is the element
♭_V w = G w of V*, and its value on a member x of V is the metric-free sum Σ c_i x_i. The metric is read once, at the
mint. A :class:`~orpheus.numerics.space.Covector` is a :class:`~orpheus.numerics.space.Member` of V*, so the raising
leg ♯ : V* → V admits it and refuses a member of V.

The laws, one row each (``vv-testing``, "a type that embodies a mathematical concept ships the test of its defining
laws"): the mint lands in V*; the pairing is the hand integral (an oracle that uses no production primitive); a member
of another space is refused; ♯ ∘ ♭ is the identity on a non-degenerate metric; linearity in the carrier and in the
density; the number one is the space's unit, which is not ``ones`` on a frame-dressed head; the raw constructor and the
specialised one give the same object (the user's ruling R1).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
from orpheus.numerics.field import Field
from orpheus.numerics.operator import RieszRaiseOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import Covector, FunctionSpace
from orpheus.numerics.units import SCALAR_FLUX_UNITS

pytestmark = pytest.mark.l0

_EPS = np.finfo(float).eps
_E = EnergyAxis.synthetic(2)
_DV = (1.0, 2.0, 3.0)
V = FunctionSpace.of_axes(_E, Axis("spatial", (3,), weights=np.array(_DV), origin=Declared(BasisKind.NODAL)))
#: Same shape as V, another spatial measure: no shape check tells the two apart.
W = FunctionSpace.of_axes(_E, Axis("spatial", (3,), weights=np.array([3.0, 2.0, 1.0]), origin=Declared(BasisKind.NODAL)))
SIGMA = np.array([[0.5, 1.0, 2.0], [3.0, 0.25, 4.0]])
PHI = np.arange(1.0, 7.0).reshape(2, 3)
#: By hand, independent of every production primitive: sum_g sum_i sigma_gi * phi_gi * dV_i.
HAND = float(sum(SIGMA[g, i] * PHI[g, i] * _DV[i] for g in range(2) for i in range(3)))


@dataclass(frozen=True, eq=False, kw_only=True)
class _Member(Field):
    """A minimal concrete member of V, built in numerics."""

    UNITS = SCALAR_FLUX_UNITS


def _phi(values=PHI, space=V) -> _Member:
    return _Member(values=np.asarray(values, dtype=float), space=space)


class TestTheCovectorIsAMemberOfTheDual:
    def test_flat_mints_a_member_of_V_dual(self):
        c = V.flat(SIGMA)
        assert isinstance(c, Covector)
        assert c.space == V.dual() and c.space != V

    def test_the_raw_and_the_specialised_constructor_give_one_object(self):
        """R1: the raw constructor, handed what the mint computes, is the same value as ``V.flat``."""
        minted = V.flat(SIGMA)
        raw = Covector(values=SIGMA * np.array(_DV), space=V.dual())
        assert raw.space == minted.space
        np.testing.assert_array_equal(raw.values, minted.values)

    def test_the_canonical_pairing_is_the_hand_integral(self):
        """Six products summed: the tolerance is six rounding units."""
        np.testing.assert_allclose(V.flat(SIGMA).evaluate(_phi()), HAND, rtol=6 * _EPS)

    def test_a_member_of_another_space_is_REFUSED(self):
        """A co-vector minted on W pairs only with members of W: the formal tier, on every run."""
        with pytest.raises(ValueError, match="not a member of"):
            W.flat(SIGMA).evaluate(_phi())

    def test_the_density_is_a_member_of_V(self):
        with pytest.raises(ValueError, match="not a member of"):
            V.flat(_phi(SIGMA, W))

    def test_linearity_in_the_carrier_and_in_the_density(self):
        a, b = 0.75, -2.5
        f1, f2 = _phi(), _phi(PHI[::-1].copy())
        lhs = V.flat(SIGMA).evaluate(a * f1 + b * f2)
        rhs = a * V.flat(SIGMA).evaluate(f1) + b * V.flat(SIGMA).evaluate(f2)
        np.testing.assert_allclose(lhs, rhs, rtol=1e-14)
        s2 = SIGMA[::-1].copy()
        np.testing.assert_allclose(
            V.flat(a * SIGMA + b * s2).evaluate(_phi()),
            a * V.flat(SIGMA).evaluate(_phi()) + b * V.flat(s2).evaluate(_phi()),
            rtol=1e-14,
        )


class TestTheRaisingLegActsOnTheDual:
    """♯ : V* → V takes a member of V*, not of V (the U1 qa review, finding F1; it waited for a typed element of V*)."""

    def test_sharp_of_flat_is_the_identity(self):
        """♯ ∘ ♭ = G⁺ G = I on a non-degenerate metric. Two roundings per entry (a multiply and a divide)."""
        back = np.asarray(RieszRaiseOperator(V).apply(V.flat(SIGMA)))
        np.testing.assert_allclose(back, SIGMA, rtol=2 * _EPS)

    def test_a_member_of_V_is_REFUSED_by_sharp(self):
        with pytest.raises(ValueError, match="not a member of"):
            RieszRaiseOperator(V).apply(_phi())


class TestTheIntegrationCovectorReadsTheUnit:
    def test_on_an_all_nodal_space_the_unit_is_ones(self):
        np.testing.assert_array_equal(V.unit, np.ones(V.shape))
        np.testing.assert_allclose(V.flat(V.unit).evaluate(_phi()), sum(PHI[g, i] * _DV[i] for g in range(2) for i in range(3)), rtol=6 * _EPS)

    def test_an_axes_less_space_has_the_unit_ones(self):
        np.testing.assert_array_equal(FunctionSpace("g", (2, 1)).unit, np.ones((2, 1)))

    def test_on_a_frame_dressed_head_the_unit_is_NOT_ones(self):
        """On a dressed head the constant function's coordinates are its moments Φᵀw (ruling 37), so ``flat(unit)``
        integrates and ``flat(ones)`` does not. Oracle: the space's own pairing ⟨u, c⟩, a different route (it does not
        build a co-vector)."""
        head = frame_faces.basis_space(frame_faces.angular_frame(Quadrature.level_symmetric(sn_order=4), 1))
        u = head.unit
        assert not np.array_equal(u, np.ones(head.shape)), "premise: the dressed head's unit is not all ones"
        c = np.random.default_rng(3).standard_normal(head.shape)
        integral = head.inner_product(u, c)
        np.testing.assert_allclose(head.flat(u).evaluate(c), integral, rtol=1e-13)
        assert not np.isclose(head.flat(np.ones(head.shape)).evaluate(c), integral)

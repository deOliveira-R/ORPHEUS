r"""The moment-head AXIS — one generic coefficient axis since S7 (ruling 81),
:meth:`Axis.of_basis <orpheus.numerics.axis.Axis.of_basis>` (CS4c step 6 item
6.2c-ii, 2026-09-08, for the memo's gates P2 and P3; until S7 two family
subclasses, ``HarmonicAxis`` and ``LegendreAxis``).

A moment head is an axis-built space with ONE MODAL axis whose measure is the
head's metric and whose origin is the object that minted it —
``FromBasis(basis)`` for the continuum head
(:meth:`FunctionSpace.of_basis <orpheus.numerics.space.FunctionSpace.of_basis>`),
``FromFrame(frame)`` for the Parseval-dressed one (``basis_space(frame)``). Two
things this file pins:

* **P2 — the domain is part of the Legendre head's IDENTITY.** `[M]` (the
  6.2c verification round, hazard H-10) the heads about ``x`` and about ``z``
  carry ``array_equal`` measures and one shape, so an identity of ``(label,
  shape, kind, weights)`` alone COLLAPSES two physically different spaces (the
  tree carries two poles). Until S7 a ``spent_axis`` field on the Legendre
  family axis separated them; since S7 the head axis keeps its domain
  :math:`S^2/O(2)_a` as its support (ruling 71). The positive leg (same axis
  ⟹ equal) and the negative leg (different pole ⟹ unequal) are asserted with
  the measures ``array_equal`` on both sides.
* **P3 — the axis's slots.** ``MODAL``; the layout is the basis's; the measure
  IS the head's metric (the continuum Gram on the basis-minted head, the
  Parseval inverse on the frame-dressed one); the origin holds the basis or the
  frame.

Foundation mark: software invariants of the axis doctrine; no physics claim.
"""

from __future__ import annotations

from orpheus.numerics.space import FunctionSpace

import numpy as np
import pytest

from orpheus.numerics.axis import Axis
from orpheus.numerics.basis.base import TruncatedBasis
from orpheus.numerics.coefficient_layout import HarmonicLayout, LegendreLayout
from orpheus.numerics.basis.legendre_basis import LegendreBasis
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.frame import GalerkinFrame
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Declared, FromBasis, FromFrame, Modal

pytestmark = pytest.mark.foundation


def _head_axis(space) -> Axis:
    assert space.axes is not None and len(space.axes) == 1, "a moment head is a single-axis space"
    return space.axes[0]


# ═════════════════════════════════════════════════════════════════════════
# P2 — the domain separates the Legendre heads
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("L", [0, 1, 3])
def test_the_domain_is_part_of_the_legendre_heads_identity(L: int) -> None:
    """Positive leg: two mints about the same axis are one axis and one space.
    Negative leg: two different poles are two axes and two spaces — with the
    measures ``array_equal`` on both sides (the separation is the domain's,
    not the measure's)."""
    x1, x2, z = FunctionSpace.of_basis(LegendreBasis(L=L, axis="x")), FunctionSpace.of_basis(LegendreBasis(L=L, axis="x")), FunctionSpace.of_basis(LegendreBasis(L=L, axis="z"))
    ax1, ax2, az = _head_axis(x1), _head_axis(x2), _head_axis(z)
    assert ax1.layout == az.layout == LegendreLayout(L)
    assert ax1.support == LegendreBasis(L=L, axis="x").domain and az.support == LegendreBasis(L=L, axis="z").domain
    assert ax1 == ax2 and hash(ax1) == hash(ax2) and x1 == x2
    # the measures agree bit-for-bit across the two poles — the family's Gram
    # does not know the axis — so the inequality below is the spent axis's alone
    assert ax1.weights is not None and az.weights is not None
    np.testing.assert_array_equal(ax1.weights, az.weights)
    assert ax1.shape == az.shape and ax1.kind == az.kind
    assert ax1 != az and x1 != z and hash(x1) != hash(z)
    assert len({x1, x2, z}) == 2


def test_a_family_generic_axis_would_collapse_the_two_poles() -> None:
    """The counterfactual (`[M]` hazard H-10): the SAME slots on a declared
    :class:`Axis` — label, shape, kind, weights, with no domain and no layout —
    compare EQUAL across the two poles, and a space built on such axes cannot
    tell ``x`` from ``z``. The head axis's domain is what separates them."""
    ax_x, ax_z = _head_axis(FunctionSpace.of_basis(LegendreBasis(L=1, axis="x"))), _head_axis(FunctionSpace.of_basis(LegendreBasis(L=1, axis="z")))
    plain_x = Axis(ax_x.label, ax_x.shape, ax_x.weights, origin=Declared(ax_x.kind))
    plain_z = Axis(ax_z.label, ax_z.shape, ax_z.weights, origin=Declared(ax_z.kind))
    assert plain_x == plain_z, "the generic identity is blind to the spent axis — that is the collapse"
    assert ax_x != ax_z
    # and a head axis is never equal to a declared axis wearing its slots
    assert ax_x != plain_x and plain_x != ax_x


# ═════════════════════════════════════════════════════════════════════════
# P3 — the axis's slots: MODAL, the layout, the head's measure, the origin
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("L", [0, 2])
def test_the_basis_minted_head_axis_carries_the_continuum_gram_and_the_basis(L: int) -> None:
    """``of_basis`` mints ONE MODAL axis: the SH head's measure is the padded
    ``4π/(2ℓ+1)`` (zero outside ``|m| ≤ ℓ``), the Legendre head's the flat
    ``4π/(2ℓ+1)``; the layout is the basis's and the origin holds the basis."""
    sh = FunctionSpace.of_basis(SphericalHarmonicBasis(L=L))
    ax = _head_axis(sh)
    assert ax.layout == HarmonicLayout(L) and isinstance(ax.kind, Modal)
    assert ax.shape == (L + 1, 2 * L + 1)
    w = ax.weights
    assert w is not None
    per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    for ell in range(L + 1):
        np.testing.assert_allclose(w[ell, : 2 * ell + 1], per_ell[ell], rtol=1e-15)
        np.testing.assert_array_equal(w[ell, 2 * ell + 1 :], 0.0)
    assert isinstance(ax.origin, FromBasis)
    basis = ax.origin.basis
    assert isinstance(basis, SphericalHarmonicBasis) and isinstance(basis, TruncatedBasis) and basis.L == L
    assert ax.support == basis.domain

    leg = FunctionSpace.of_basis(LegendreBasis(L=L, axis="z"))
    lax = _head_axis(leg)
    assert lax.layout == LegendreLayout(L) and isinstance(lax.kind, Modal) and lax.shape == (L + 1,)
    assert lax.weights is not None
    np.testing.assert_allclose(lax.weights, per_ell, rtol=1e-15)
    assert isinstance(lax.origin, FromBasis)
    lbasis = lax.origin.basis
    assert isinstance(lbasis, LegendreBasis) and lbasis.axis == "z" and lbasis.L == L
    # the space's metric verbs read the axis's measure — the head has no other source
    assert sh.inner_product_weights is None and sh.metric is None
    assert leg.inner_product_weights is None and leg.metric is None


@pytest.mark.parametrize("label", ["gauss_legendre(8)", "lebedev(11)", "level_symmetric(4)"])
def test_the_frame_dressed_head_axis_carries_the_parseval_measure_and_the_frame(label: str) -> None:
    """``basis_space(frame)`` re-weights the head axis with the Parseval inverse
    of the discrete Gram's diagonal (zero on dead slots), keeps the layout, and
    makes the FRAME the origin — the object that can re-dress the head at
    another order."""
    rules = {
        "gauss_legendre(8)": lambda: Quadrature.gauss_legendre(8),
        "lebedev(11)": lambda: Quadrature.lebedev(11),
        "level_symmetric(4)": lambda: Quadrature.level_symmetric(4),
    }
    frame = frame_faces.angular_frame(rules[label](), 2)
    dressed = frame_faces.basis_space(frame)
    ax = _head_axis(dressed)
    cont = _head_axis(FunctionSpace.of_basis(frame.basis))
    assert ax.layout == cont.layout and isinstance(ax.kind, Modal) and ax.shape == cont.shape
    assert isinstance(ax.origin, FromFrame) and ax.origin.frame is frame
    assert ax.weights is not None and cont.weights is not None
    diag = np.diagonal(frame.discrete_gram).reshape(ax.shape)
    live = diag > 0.0
    np.testing.assert_allclose(ax.weights[live], 1.0 / diag[live], rtol=1e-15)
    np.testing.assert_array_equal(ax.weights[~live], 0.0)
    assert not np.array_equal(ax.weights, cont.weights)
    assert ax != cont and dressed != FunctionSpace.of_basis(frame.basis)


def test_the_harmonic_head_axis_is_never_a_declared_axis_wearing_its_slots() -> None:
    """A declared ``Axis`` with identical slots is a different axis: it has no
    domain and no layout (ruling 71 and 81), so a hand-built head can never
    pass for the family's (until S7 the family was a subclass)."""
    ax = _head_axis(FunctionSpace.of_basis(SphericalHarmonicBasis(L=1)))
    plain = Axis(ax.label, ax.shape, ax.weights, origin=Declared(ax.kind))
    assert ax != plain and plain != ax
    assert hash(ax) != hash(plain) or ax != plain

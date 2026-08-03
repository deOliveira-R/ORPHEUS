r"""The defining laws of :mod:`orpheus.numerics.exactness`.

An exactness claim is a mathematical object, so it is tested against
the laws that define it rather than against the way callers happen to
use it:

* ``combined_with`` is a **meet** — commutative, associative,
  idempotent, and order-reversing in the degree;
* it is **partial**: claims about different reference measures have no
  common refinement, and the honest answer is ``None`` rather than a
  smaller number;
* :attr:`ExactnessClaim.system` is **read** from the reference, so the
  claim and the space it quantifies over cannot drift apart;
* a degree is meaningless below ``0``.

Plus the structural fact that makes the design work:
:class:`~orpheus.numerics.generating_measure.GeneratingMeasure`
satisfies :class:`ReferenceMeasure` **without being told to**.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from orpheus.numerics.exactness import (
    ExactnessClaim,
    OrthogonalSystem,
    ReferenceMeasure,
)
from orpheus.numerics.generating_measure import (
    CHEBYSHEV_T,
    CHEBYSHEV_U,
    HERMITE,
    LEGENDRE,
    jacobi,
)
from orpheus.numerics.measure import SPACE_CIRCLE, Space

pytestmark = [pytest.mark.foundation]


@dataclass(frozen=True)
class _FourierRef:
    """A minimal non-generating reference — the circle's Fourier system.

    Deliberately NOT a ``GeneratingMeasure``: it has no three-term
    recurrence, which is the whole reason the claim is typed against
    the protocol. Its presence here proves the protocol admits a second,
    structurally different realization (the project's ">= 2 realizations"
    criterion), rather than being a one-instance abstraction.
    """

    name: str = "uniform_on_circle"
    support: Space = SPACE_CIRCLE
    orthogonal_system: OrthogonalSystem = OrthogonalSystem.TRIGONOMETRIC


# ── the protocol is satisfied structurally ──────────────────────────────


def test_generating_measure_satisfies_reference_measure() -> None:
    """A Gauss generator IS a reference measure — by having the
    attributes, not by inheriting."""
    for gm in (LEGENDRE, CHEBYSHEV_T, CHEBYSHEV_U, HERMITE, jacobi(1.5, 2.25)):
        assert isinstance(gm, ReferenceMeasure)
        assert gm.orthogonal_system is OrthogonalSystem.ALGEBRAIC


def test_a_second_structurally_different_reference_exists() -> None:
    """The protocol has >= 2 non-isomorphic realizations.

    Without this the abstraction would be a one-instance ceremony; the
    Fourier reference has no recurrence and a different system, which is
    exactly the case ``GeneratingMeasure`` cannot represent.
    """
    ref = _FourierRef()
    assert isinstance(ref, ReferenceMeasure)
    assert ref.orthogonal_system is not LEGENDRE.orthogonal_system


# ── the claim reads its system; it cannot drift ─────────────────────────


def test_system_is_read_from_the_reference() -> None:
    """``system`` is not stored, so it cannot disagree with the reference."""
    assert ExactnessClaim(LEGENDRE, 7).system is OrthogonalSystem.ALGEBRAIC
    assert (
        ExactnessClaim(_FourierRef(), 7).system
        is OrthogonalSystem.TRIGONOMETRIC
    )


def test_the_two_failure_modes_are_now_distinguishable() -> None:
    r"""The header's two measured bugs become different objects.

    Same degree, different *measure* (GL vs GC — the 0.696 miss), and
    same degree, different *space* (algebraic vs trigonometric — the
    ``min()``-gives-1 bug). As bare integers all three were ``7``.
    """
    gl = ExactnessClaim(LEGENDRE, 7)
    gc = ExactnessClaim(CHEBYSHEV_T, 7)
    trig = ExactnessClaim(_FourierRef(), 7)

    assert gl.degree == gc.degree == trig.degree == 7      # the old view
    assert gl != gc and gl != trig and gc != trig          # the new one
    assert gl.system is gc.system                          # same space...
    assert gl.reference != gc.reference                    # ...different measure
    assert gl.system is not trig.system                    # different space


# ── ``combined_with`` is a partial meet ─────────────────────────────────


@pytest.mark.parametrize("p,q", [(0, 0), (0, 5), (3, 7), (7, 3), (9, 9)])
def test_combined_with_takes_the_minimum_degree(p: int, q: int) -> None:
    a, b = ExactnessClaim(LEGENDRE, p), ExactnessClaim(LEGENDRE, q)
    combined = a.combined_with(b)
    assert combined is not None
    assert combined.degree == min(p, q)
    assert combined.reference == LEGENDRE


def test_combined_with_is_commutative() -> None:
    a, b = ExactnessClaim(LEGENDRE, 3), ExactnessClaim(LEGENDRE, 7)
    assert a.combined_with(b) == b.combined_with(a)


def test_combined_with_is_idempotent() -> None:
    a = ExactnessClaim(LEGENDRE, 5)
    assert a.combined_with(a) == a


def test_combined_with_is_associative() -> None:
    a = ExactnessClaim(LEGENDRE, 9)
    b = ExactnessClaim(LEGENDRE, 4)
    c = ExactnessClaim(LEGENDRE, 6)
    left = a.combined_with(b)
    right = b.combined_with(c)
    assert left is not None and right is not None
    assert left.combined_with(c) == a.combined_with(right)


def test_combining_never_strengthens_a_claim() -> None:
    """The meet is a lower bound on both — an exactness composition can
    only ever weaken, which is the direction that keeps it sound."""
    a, b = ExactnessClaim(LEGENDRE, 3), ExactnessClaim(LEGENDRE, 11)
    combined = a.combined_with(b)
    assert combined is not None
    assert combined.degree <= a.degree and combined.degree <= b.degree


# ── partiality: the condition that is load-bearing, not defensive ───────


def test_different_references_have_NO_common_claim() -> None:
    r"""`[M]` ``gauss_legendre(4) * gauss_chebyshev(4)`` advertised degree
    7 while integrating the constant 1 to 6.2832 instead of 4.

    Each factor is exact against its OWN weight; the product of two
    different weights is a measure neither ever claimed anything about.
    ``None`` is the honest answer — a smaller number would still be a
    claim about nothing.
    """
    gl = ExactnessClaim(LEGENDRE, 7)
    gc = ExactnessClaim(CHEBYSHEV_T, 7)
    assert gl.combined_with(gc) is None
    assert gc.combined_with(gl) is None


def test_different_SPACES_have_no_common_claim_either() -> None:
    """A polar (algebraic) and an azimuthal (trigonometric) claim do not
    combine by ``min`` — their composition is a theorem about the target
    space, and this method correctly declines to guess it."""
    polar = ExactnessClaim(LEGENDRE, 7)
    azimuthal = ExactnessClaim(_FourierRef(), 7)
    assert polar.combined_with(azimuthal) is None


def test_canonically_equal_references_DO_combine() -> None:
    """``jacobi(0, 0)`` IS Legendre, so their claims compose.

    Equality is on the mathematical identity, not the code path that
    built it — so a rule constructed either way carries a claim that
    composes with the other.
    """
    combined = ExactnessClaim(LEGENDRE, 7).combined_with(
        ExactnessClaim(jacobi(0.0, 0.0), 5)
    )
    assert combined is not None
    assert combined.degree == 5


# ── the degree's own invariant ──────────────────────────────────────────


@pytest.mark.parametrize("bad", [-1, -7])
def test_a_negative_degree_is_not_a_claim(bad: int) -> None:
    """A rule with no exact subspace carries ``None``, never a negative
    degree — the illegal state is unrepresentable rather than merely
    unusual."""
    with pytest.raises(ValueError, match="degree >= 0"):
        ExactnessClaim(LEGENDRE, bad)


def test_degree_zero_is_a_legitimate_claim() -> None:
    """Positive control for the guard above: exact for the constant only
    is a real, weak claim — the guard must not reject it."""
    claim = ExactnessClaim(LEGENDRE, 0)
    assert claim.degree == 0


def test_claim_is_frozen() -> None:
    claim = ExactnessClaim(LEGENDRE, 7)
    with pytest.raises(Exception):
        claim.degree = 9  # type: ignore[misc]


def test_str_names_all_three_parts() -> None:
    """The rendering must carry degree, system AND reference — reading
    any one alone is what this type exists to prevent."""
    rendered = str(ExactnessClaim(CHEBYSHEV_T, 7))
    assert "7" in rendered
    assert "algebraic" in rendered
    assert "chebyshev_t" in rendered

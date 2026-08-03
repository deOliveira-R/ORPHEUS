r"""The defining laws of :mod:`orpheus.numerics.exactness`.

An exactness claim is a mathematical object, so it is tested against
the laws that define it rather than against the way callers happen to
use it:

* ``tensor_with`` meets the **degree** and multiplies the
  **reference** — so it is order-reversing in the degree but is NOT
  idempotent, because ``a ⊗ a`` lands on the square;
* it is **partial**: factors whose orthogonal SYSTEMS differ have no
  product claim, because a mixed tensor system's degree is a theorem
  about the target space rather than a minimum;
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
    ProductMeasure,
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


# ── ``tensor_with``: a meet on the DEGREE, a product on the REFERENCE ──
#
# NOT a meet as a whole, and the distinction is worth stating because an
# earlier draft of these tests assumed it was. Only the degree meets; the
# reference ACCUMULATES into a ``ProductMeasure``. So ``tensor_with`` is
# neither idempotent (a ⊗ a is a rule on a SQUARE, not on the interval)
# nor associative as an object (the references nest differently), and
# asserting either would be asserting a law the operation does not have.


@pytest.mark.parametrize("p,q", [(0, 0), (0, 5), (3, 7), (7, 3), (9, 9)])
def test_tensor_degree_is_the_minimum(p: int, q: int) -> None:
    """Separable integrands of degree ``d`` per axis need BOTH factors
    exact to ``d``, so the surviving degree is the smaller."""
    a, b = ExactnessClaim(LEGENDRE, p), ExactnessClaim(LEGENDRE, q)
    product = a.tensor_with(b)
    assert product is not None
    assert product.degree == min(p, q)


def test_tensor_reference_is_the_PRODUCT_not_a_factor() -> None:
    """The load-bearing property: a product claims the product measure.

    Keeping a factor's reference would assert exactness against a measure
    the product is not exact against — the error the direct sum still
    demonstrates by carrying no claim at all.
    """
    product = ExactnessClaim(LEGENDRE, 3).tensor_with(
        ExactnessClaim(CHEBYSHEV_T, 7)
    )
    assert product is not None
    assert isinstance(product.reference, ProductMeasure)
    assert product.reference != LEGENDRE
    assert product.reference != CHEBYSHEV_T
    assert product.reference.name == "legendre ⊗ chebyshev_t"
    assert product.reference.support == "[-1,1] × [-1,1]"


def test_tensor_is_NOT_idempotent() -> None:
    """``a ⊗ a`` is a rule on the SQUARE, not the original rule.

    Pins that the operation is not being treated as a meet: a meet would
    satisfy ``a ∧ a == a``, and asserting that here would be asserting a
    law this operation does not have.
    """
    a = ExactnessClaim(LEGENDRE, 5)
    squared = a.tensor_with(a)
    assert squared is not None
    assert squared != a
    assert squared.degree == a.degree
    assert squared.reference.support == "[-1,1] × [-1,1]"


def test_tensor_degree_never_strengthens() -> None:
    """A composition can only ever weaken the degree — the direction that
    keeps the claim sound."""
    a, b = ExactnessClaim(LEGENDRE, 3), ExactnessClaim(LEGENDRE, 11)
    product = a.tensor_with(b)
    assert product is not None
    assert product.degree <= a.degree and product.degree <= b.degree


# ── partiality: which pairs have NO product claim, and why ──────────────


def test_different_references_DO_compose_when_the_system_agrees() -> None:
    r"""⚠ Inverted 2026-08-02, and the inversion is the finding.

    This previously asserted that Legendre ⊗ Chebyshev has **no** claim,
    citing "it integrates the constant 1 to 6.2832 where the answer is
    4". `[M]` Re-measured: :math:`2\pi = 6.2832` **is** the mass of
    ``legendre ⊗ chebyshev_t``, and the rule IS exact to degree 7 per
    axis against that measure (4.16e-13, tight at degree 8). The old
    expectation of 4 was the *Lebesgue* product — the wrong reference.

    The refusal was a workaround for having no way to NAME the product
    measure, not a law. With the name available the correct claim is
    representable, which is the whole point of the carve.
    """
    gl = ExactnessClaim(LEGENDRE, 7)
    gc = ExactnessClaim(CHEBYSHEV_T, 7)
    product = gl.tensor_with(gc)
    assert product is not None
    assert product.degree == 7
    assert product.reference.name == "legendre ⊗ chebyshev_t"


def test_different_SYSTEMS_have_no_product_claim() -> None:
    """A polar (algebraic) and an azimuthal (trigonometric) claim do NOT
    compose by ``min``.

    This is the case that survives as partial, and it is the one that
    matters: their product spans a mixed tensor system whose degree comes
    from a theorem about the target space (which monomials on
    :math:`S^2` factor into which polar and azimuthal degrees). Guessing
    a minimum here is exactly the ``min()``-gives-1 bug.
    """
    polar = ExactnessClaim(LEGENDRE, 7)
    azimuthal = ExactnessClaim(_FourierRef(), 7)
    assert polar.tensor_with(azimuthal) is None
    assert azimuthal.tensor_with(polar) is None


def test_a_mixed_product_reference_is_unconstructible() -> None:
    """The guard is on the type, not only on the composing method — so a
    mixed product cannot be built by going around ``tensor_with``."""
    with pytest.raises(ValueError, match="share an orthogonal system"):
        ProductMeasure(factors=(LEGENDRE, _FourierRef()))


def test_canonically_equal_references_compose_as_one_measure() -> None:
    """``jacobi(0, 0)`` IS Legendre, so the product names it once per
    factor rather than treating the two constructions as different."""
    product = ExactnessClaim(LEGENDRE, 7).tensor_with(
        ExactnessClaim(jacobi(0.0, 0.0), 5)
    )
    assert product is not None
    assert product.degree == 5
    assert product.reference.name == "legendre ⊗ legendre"


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

r"""The two zeros after the split: ``ZeroMorphism``'s bound ends, ``ZeroOperator``'s echo.

The zero map is the one operator whose ACTION cannot reveal which spaces it
connects — every input goes to zero, so no probe distinguishes
:math:`0 : \mathcal D \to \mathcal C` from :math:`0` on any other pair. The
S4-amendment (2026-08-22) split the two roles the old single class straddled:

* :class:`~orpheus.numerics.operator.ZeroMorphism` — the zero MAP between two
  DECLARED spaces, which therefore DEMANDS both at construction and derives
  its zeros from the bound shapes (the per-site ``codomain_zero`` /
  ``transpose_zero`` closures retired with their consumers);
* :class:`~orpheus.numerics.operator.ZeroOperator` — the ENDOMORPHIC zero,
  multiplication by ``0``: a stateless :class:`PointwiseOperator` member whose
  ``0.0 * x`` echo IS the zero of the operand's own space.

⭐ **Why the unequal-ends rows exist — the ``vv`` Mode-12 argument, inherited.**
These rows are the home of a claim that used to live on the retired
``IncomingSourceOperator`` (retired at **P3**, 2026-08-05):

    ``|Γ₊| == |Γ₋|`` on every quadrature × face pair in the tree, so a fixture
    where the two agree cannot tell "emits the CODOMAIN's zero" from "echoes
    the input's shape" — the error class sits inside the shape functional's
    invariance group. Unequal sizes are the only way to see it.

That hazard did NOT retire with the closures — it is now what
:class:`ZeroMorphism`'s bound-shape derivation exists to prevent, and
:func:`~orpheus.sn.boundary.realizer._narrowed_zero_operator`'s own docstring
still says relying on the endomorphic echo would be "wrong in principle and
merely lucky in practice". So the discrimination is made where unequal ends
are CONSTRUCTIBLE, which is here: the pair below is deliberately ``3 ≠ 7``.

V&V tags
========

``@pytest.mark.foundation`` — these are software-invariant claims (which array
shape comes out, which space object is reported), not numerical contracts. The
references are closed-form: the shapes are fixed by the bound spaces.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    PointwiseOperator,
    ZeroMorphism,
    ZeroOperator,
)
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation

#: Deliberately UNEQUAL, and deliberately different from any probe's leading
#: axis: 3 rows out of the forward, 7 out of the transpose, probed at 7 and 3.
#: (The inherited Mode-12 discrimination — see the module docstring.)
_N_CODOMAIN = 3
_N_DOMAIN = 7


def _space(name: str, n: int) -> FunctionSpace:
    """A full-shape space — ``(name, shape)`` IS its identity. The trailing
    ``2`` mirrors the trace seam's group axis, so the bound-shape derivation
    is exercised on a genuinely multi-axis shape."""
    return FunctionSpace(name, (n, 2))


def _bound_zero() -> ZeroMorphism:
    return ZeroMorphism(
        domain=_space("gamma_plus", _N_DOMAIN),
        codomain=_space("gamma_minus", _N_CODOMAIN),
    )


class TestTheForwardEmitsTheCodomainsZero:
    r"""``0 : \mathcal D \to \mathcal C`` produces an element of ``\mathcal C``."""

    def test_the_row_count_is_the_codomains_not_the_inputs(self) -> None:
        """⭐ The inherited Mode-12 row: ``3`` out for a ``7``-row probe."""
        op = _bound_zero()
        out = op.apply(np.ones((_N_DOMAIN, 2)))
        assert out.shape == (_N_CODOMAIN, 2), (
            f"emitted {out.shape} for a {_N_DOMAIN}-row probe with a "
            f"{_N_CODOMAIN}-row codomain — the forward is echoing its INPUT's "
            f"shape, which is the endomorphic `0.0 * x` default leaking "
            f"through a genuine map between spaces."
        )
        np.testing.assert_array_equal(out, 0.0)

    def test_the_transpose_lands_in_the_domain(self) -> None:
        """And the mirror: ``7`` out for a ``3``-row probe."""
        op = _bound_zero()
        out = op.apply_transpose(np.ones((_N_CODOMAIN, 2)))
        assert out.shape == (_N_DOMAIN, 2), (
            f"the transpose emitted {out.shape} — it must land in the DOMAIN "
            f"({_N_DOMAIN} rows). With the two ends swapped this row is the "
            f"only thing that can tell, since the values are zero either way."
        )
        np.testing.assert_array_equal(out, 0.0)

    def test_the_two_directions_are_not_interchangeable(self) -> None:
        """The negative leg: the forward and transpose shapes must DIFFER.

        Without this, a future edit that bound both ends to the same space
        would leave every row above green while destroying the distinction
        they exist to make.
        """
        op = _bound_zero()
        fwd = op.apply(np.ones((_N_DOMAIN, 2))).shape
        bwd = op.apply_transpose(np.ones((_N_CODOMAIN, 2))).shape
        assert fwd != bwd, (
            "this fixture's two ends have equal extent, so it cannot "
            "discriminate direction — pick unequal bound shapes"
        )


class TestTheTwoSpacesAreDemanded:
    """The S4-amendment's sharpest instance of the base demand: the zero map
    cannot be probed into a pair, so the pair is REQUIRED at construction."""

    def test_the_bound_spaces_are_reported_by_identity(self) -> None:
        domain = _space("gamma_plus", _N_DOMAIN)
        codomain = _space("gamma_minus", _N_CODOMAIN)
        op = ZeroMorphism(domain=domain, codomain=codomain)
        # ``is``, not ``==``: FunctionSpace equality is (name, shape), so
        # identity is the strictly stronger claim — and it is what the
        # composability check actually compares.
        assert op.domain is domain
        assert op.codomain is codomain

    def test_the_ends_are_not_swapped(self) -> None:
        r"""⭐ The only row a swap reddens.

        A zero map with its two ends exchanged computes the same values, has
        the same predicates, and composes against the same shapes whenever the
        extents agree. Naming which space is which END is the whole content of
        the binding, so the gate has to say it — and it has to say it with
        DIFFERENT extents, or a swap is not even shape-visible.
        """
        domain = _space("gamma_plus", _N_DOMAIN)
        codomain = _space("gamma_minus", _N_CODOMAIN)
        op = ZeroMorphism(domain=domain, codomain=codomain)
        assert op.domain is not codomain
        assert op.codomain is not domain
        # …and the reported extents match the direction the actions emit in.
        assert op.domain.shape[0] == op.apply_transpose(
            np.ones((_N_CODOMAIN, 2)),
        ).shape[0]
        assert op.codomain.shape[0] == op.apply(
            np.ones((_N_DOMAIN, 2)),
        ).shape[0]

    def test_a_missing_space_is_unconstructable(self) -> None:
        """The demand's teeth: no half-bound zero map exists. (The old
        optional binding — 'the composability check skips a None end by
        design' — is exactly the silent state the amendment retired for
        this class; the endomorphic role moved to ZeroOperator.)"""
        with pytest.raises(TypeError):
            ZeroMorphism(domain=_space("d", _N_DOMAIN))  # type: ignore[call-arg]
        with pytest.raises(TypeError):
            ZeroMorphism(codomain=_space("c", _N_CODOMAIN))  # type: ignore[call-arg]
        with pytest.raises(TypeError, match="demands BOTH spaces"):
            ZeroMorphism(domain=None, codomain=None)  # type: ignore[arg-type]

    def test_the_adjoint_is_the_swapped_zero_map(self) -> None:
        """``0^* : C → D`` — algebra-closed and metric-free (a zero is a
        zero in any inner product), so no sandwich wrapper is built."""
        op = _bound_zero()
        adj = op.H
        assert isinstance(adj, ZeroMorphism)
        assert adj.domain is op.codomain
        assert adj.codomain is op.domain
        np.testing.assert_array_equal(
            adj.apply(np.ones((_N_CODOMAIN, 2))), np.zeros((_N_DOMAIN, 2)),
        )


class TestTheEndomorphicZeroIsPointwise:
    """The other half of the split: the stateless ×0 multiplier."""

    def test_the_endomorphic_echo(self) -> None:
        """``0.0 * x`` IS the zero of the operand's own space — correct
        precisely because a pointwise member is an endomorphism there."""
        op = ZeroOperator()
        probe = np.ones((5, 2))
        np.testing.assert_array_equal(op.apply(probe), np.zeros_like(probe))
        np.testing.assert_array_equal(
            op.apply_transpose(probe), np.zeros_like(probe),
        )

    def test_membership_and_the_family_laws(self) -> None:
        """A PointwiseOperator member: no stored pair (None BY LAW,
        discriminated by type) and a metric-free self adjoint."""
        op = ZeroOperator()
        assert isinstance(op, PointwiseOperator)
        assert op.domain is None
        assert op.codomain is None
        assert op.H is op


class TestTheStructuralAxesAreUnchanged:
    """The split must not disturb the two predicate axes, on either class."""

    def test_adjointable_and_singular(self) -> None:
        for op in (_bound_zero(), ZeroOperator()):
            assert op.is_adjointable, "0ᵀ = 0"
            assert not op.is_invertible, "the zero map is singular"
            assert not hasattr(op, "inverse"), (
                f"{type(op).__name__} declares no inverse() at all — a "
                f"raising stub would be the harmful-stub anti-pattern "
                f"numerics.operator is built against"
            )

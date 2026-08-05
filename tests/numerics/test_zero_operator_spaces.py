r"""``ZeroOperator``'s two ends: the codomain's zero, and the two named spaces.

The zero map is the one operator whose ACTION cannot reveal which spaces it
connects — every input goes to zero, so no probe distinguishes
:math:`0 : \mathcal D \to \mathcal C` from :math:`0` on any other pair. Two
independent claims follow, and this module is where each is reachable:

1. **The forward emits the zero of the CODOMAIN, the transpose the zero of the
   DOMAIN.** Wrong-direction is invisible whenever the two ends have equal
   extent, which is the case on *every* reachable SN face — see below.
2. **The two spaces are NAMED** (G6.3 step 6, #330), so a zero slot composes
   under the ``OperatorSum`` / ``OperatorProduct`` composability check instead
   of short-circuiting it with ``None``.

⭐ **Why this module exists — the ``vv`` Mode-12 argument, inherited.** These
rows are the new home of a claim that used to live on the retired
``IncomingSourceOperator``
(``tests/geometry/test_bc_universal_invariants.py::TestIncomingSourceOperator::
test_apply_fills_the_codomain_not_the_input_shape``, retired at **P3**,
2026-08-05). That test's argument was:

    ``|Γ₊| == |Γ₋|`` on every quadrature × face pair in the tree, so a fixture
    where the two agree cannot tell "emits the CODOMAIN's zero" from "echoes
    the input's shape" — the error class sits inside the shape functional's
    invariance group. Unequal sizes are the only way to see it.

That hazard did NOT retire with the operator. It is precisely what
:class:`~orpheus.numerics.operator.ZeroOperator`'s two space hooks exist to
prevent, and
:func:`~orpheus.sn.boundary.realizer._narrowed_zero_operator`'s own docstring
says relying on the endomorphic ``0.0 * x`` echo would be "wrong in principle
and merely lucky in practice: ``|Γ₊| == |Γ₋|`` on every reachable fixture, so
the shapes coincide — an accident, not a contract."

So the discrimination has to be made where unequal ends are CONSTRUCTIBLE,
which is here: ``ZeroOperator`` is hand-buildable with whatever pair of hooks
you hand it, and the pair below is deliberately ``3 ≠ 7``.

V&V tags
========

``@pytest.mark.foundation`` — these are software-invariant claims (which array
shape comes out, which space object is reported), not numerical contracts. The
references are closed-form: the shapes are fixed by the hooks' definitions.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import ZeroOperator
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation

#: Deliberately UNEQUAL, and deliberately different from any probe's leading
#: axis: 3 rows out of the forward, 7 out of the transpose, probed at 7 and 3.
#: Equal sizes here would make every row below vacuous (see the module
#: docstring's Mode-12 argument).
_N_CODOMAIN = 3
_N_DOMAIN = 7


def _zero_rows(n_rows: int):
    """A space hook emitting ``n_rows`` leading entries, trailing axes intact.

    Mirrors ``orpheus.sn.boundary.realizer._zero_rows``. Spelled locally rather
    than imported so the gate does not depend on the SN boundary tier to state
    a numerics-tier claim.

    ``reportArgumentType`` at the call sites below: ``ZeroOperator`` is generic
    over ``Vector`` and MEASURED (2026-07-31) ``np.ndarray`` does not satisfy
    that protocol statically. The same annotation sits on the production hook in
    ``orpheus/sn/boundary/realizer.py``; the gap is upstream, not here.
    """
    def build(x):
        arr = np.asarray(x)
        return np.zeros((n_rows,) + arr.shape[1:], dtype=arr.dtype)

    return build


def _space(name: str, n: int) -> FunctionSpace:
    """A bare space of the given extent — ``(name, shape)`` IS its identity."""
    return FunctionSpace(name, (n,))


class TestTheForwardEmitsTheCodomainsZero:
    r"""``0 : \mathcal D \to \mathcal C`` produces an element of ``\mathcal C``."""

    def test_the_row_count_is_the_codomains_not_the_inputs(self) -> None:
        """⭐ The inherited Mode-12 row: ``3`` out for a ``7``-row probe."""
        op = ZeroOperator(
            codomain_zero=_zero_rows(_N_CODOMAIN),  # type: ignore[reportArgumentType]
            transpose_zero=_zero_rows(_N_DOMAIN),  # type: ignore[reportArgumentType]
        )
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
        op = ZeroOperator(
            codomain_zero=_zero_rows(_N_CODOMAIN),  # type: ignore[reportArgumentType]
            transpose_zero=_zero_rows(_N_DOMAIN),  # type: ignore[reportArgumentType]
        )
        out = op.apply_transpose(np.ones((_N_CODOMAIN, 2)))
        assert out.shape == (_N_DOMAIN, 2), (
            f"the transpose emitted {out.shape} — it must land in the DOMAIN "
            f"({_N_DOMAIN} rows). With the two ends swapped this row is the "
            f"only thing that can tell, since the values are zero either way."
        )
        np.testing.assert_array_equal(out, 0.0)

    def test_the_two_directions_are_not_interchangeable(self) -> None:
        """The negative leg: the forward and transpose shapes must DIFFER.

        Without this, a future edit that wired both hooks to the same size
        would leave every row above green while destroying the distinction
        they exist to make.
        """
        op = ZeroOperator(
            codomain_zero=_zero_rows(_N_CODOMAIN),  # type: ignore[reportArgumentType]
            transpose_zero=_zero_rows(_N_DOMAIN),  # type: ignore[reportArgumentType]
        )
        fwd = op.apply(np.ones((_N_DOMAIN, 2))).shape
        bwd = op.apply_transpose(np.ones((_N_CODOMAIN, 2))).shape
        assert fwd != bwd, (
            "this fixture's two ends have equal extent, so it cannot "
            "discriminate direction — pick unequal hook sizes"
        )

    def test_the_endomorphic_default_still_echoes(self) -> None:
        """Positive control for the hooks: WITHOUT them, the echo is correct.

        A zero operator with no hooks is an endomorphism (domain == codomain),
        where ``0.0 * x`` IS the zero of the codomain. Pinned so the rows above
        read as "the hooks change the behaviour", not as "the default was
        broken" — and so a change that made the hooks mandatory reddens here.
        """
        op = ZeroOperator()
        probe = np.ones((5, 2))
        np.testing.assert_array_equal(op.apply(probe), np.zeros_like(probe))
        np.testing.assert_array_equal(
            op.apply_transpose(probe), np.zeros_like(probe),
        )


class TestTheTwoSpacesAreNamed:
    """G6.3 step 6: a zero slot reports its ends, so compositions can check."""

    def test_the_bound_spaces_are_reported_by_identity(self) -> None:
        domain = _space("gamma_plus", _N_DOMAIN)
        codomain = _space("gamma_minus", _N_CODOMAIN)
        op = ZeroOperator(
            codomain_zero=_zero_rows(_N_CODOMAIN),  # type: ignore[reportArgumentType]
            transpose_zero=_zero_rows(_N_DOMAIN),  # type: ignore[reportArgumentType]
            domain=domain,
            codomain=codomain,
        )
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
        op = ZeroOperator(
            codomain_zero=_zero_rows(_N_CODOMAIN),  # type: ignore[reportArgumentType]
            transpose_zero=_zero_rows(_N_DOMAIN),  # type: ignore[reportArgumentType]
            domain=domain,
            codomain=codomain,
        )
        assert op.domain is not codomain
        assert op.codomain is not domain
        # …and the reported extents match the direction the hooks emit in.
        assert op.domain.shape[0] == op.apply_transpose(
            np.ones((_N_CODOMAIN, 2)),
        ).shape[0]
        assert op.codomain.shape[0] == op.apply(
            np.ones((_N_DOMAIN, 2)),
        ).shape[0]

    def test_an_unbound_zero_reports_none(self) -> None:
        """The pre-#208 default survives: binding stays OPTIONAL.

        Most zero slots in the tree are unbound (the within-group zero fission
        slot among them), and the composability check skips a ``None`` end by
        design. A change that made the spaces required would break those
        callers, so the default is asserted rather than assumed.
        """
        op = ZeroOperator(codomain_zero=_zero_rows(_N_CODOMAIN))
        assert op.domain is None
        assert op.codomain is None


class TestTheStructuralAxesAreUnchanged:
    """The binding must not disturb the two predicate axes."""

    def test_adjointable_and_singular(self) -> None:
        op = ZeroOperator(
            codomain_zero=_zero_rows(_N_CODOMAIN),  # type: ignore[reportArgumentType]
            transpose_zero=_zero_rows(_N_DOMAIN),  # type: ignore[reportArgumentType]
            domain=_space("d", _N_DOMAIN),
            codomain=_space("c", _N_CODOMAIN),
        )
        assert op.is_adjointable, "0ᵀ = 0"
        assert not op.is_invertible, "the zero map is singular"
        assert not hasattr(op, "inverse"), (
            "ZeroOperator declares no inverse() at all — a raising stub would "
            "be the harmful-stub anti-pattern numerics.operator is built "
            "against"
        )

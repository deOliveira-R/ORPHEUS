r"""One parse at the domain: coordinate-space membership, on every run (U1b, plan ``space_membership_covector_gauge.md``).

The parse has two legs, two different facts:

* the SHAPE leg, "x is an element of the coordinate space R^shape": the primitive
  :func:`~orpheus.numerics.coordinate_space.coordinates` (exact) and
  :func:`~orpheus.numerics.coordinate_space.leading_coordinates` (an action :math:`A \otimes I_P`: the leading axes are
  R^shape, a trailing payload rides through). It needs no space and sits at the bottom of the concept order;
* the IDENTITY leg, "which space": a typed member compared with :class:`~orpheus.numerics.space.FunctionSpace` by
  content identity.

There are no ``__debug__`` rows: the check is the parse, not a guard, and every row runs under the canonical
``python -O -m pytest``.
"""
from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

from orpheus.numerics import coordinate_space as cs

pytestmark = pytest.mark.l0

X = np.arange(1.0, 7.0).reshape(2, 3)

#: The two fragments, disjoint (asserted below). Both are the vocabulary the tree already pins
#: (`retirement-audit` B.11: prefer the established vocabulary).
EXACT = "exact shape"
LEADING = "leading shape"


class TestThePrimitive:
    r"""Membership in R^shape is a shape equality (exact), or an equality of the leading axes (an action)."""

    def test_the_fragments_are_disjoint(self):
        assert EXACT not in LEADING and LEADING not in EXACT

    def test_the_module_sits_below_every_concept(self):
        """The concept order gains a first row ``("coordinate_space",)``; the module imports nothing from ``orpheus``,
        so measure, basis and metric can call it."""
        from tests.gates.numerics import test_concept_order as order

        assert order.RANK["coordinate_space"] < min(order.RANK[m] for m in ("metric", "basis", "measure", "space"))
        tree = ast.parse(inspect.getsource(cs))
        imported = [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
        imported += [n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
        assert not [m for m in imported if m.split(".")[0] == "orpheus"], imported

    def test_the_parse_has_no_debug_tier(self):
        """The check survives ``-O`` because nothing gates it: no ``__debug__`` name and no ``assert`` in the module."""
        tree = ast.parse(inspect.getsource(cs))
        assert not any(isinstance(n, ast.Name) and n.id == "__debug__" for n in ast.walk(tree))
        assert not any(isinstance(n, ast.Assert) for n in ast.walk(tree))

    def test_an_exact_member_is_returned_as_itself(self):
        """Identity, not equality: the parse returns the operand object, so the multiply and reduction trees downstream
        are the ones they were."""
        assert cs.coordinates(X, (2, 3), owner="t") is X
        assert cs.leading_coordinates(X, (2,), owner="t") is X

    @pytest.mark.parametrize("shape", [(2, 1), (1, 3), (3,), (6,), (2, 3, 1), (1, 2, 3)], ids=lambda s: "x".join(map(str, s)))
    def test_every_non_member_is_refused_by_the_exact_parse(self, shape):
        """Each shape is in numpy's broadcast closure of (2, 3) or a reshape of it: the inputs ERR-122 (outer product),
        the implicit pullback x⊗1 and the non-member-out classes are made of."""
        with pytest.raises(cs.NotInCoordinateSpace, match=EXACT):
            cs.coordinates(np.ones(shape), (2, 3), owner="t")

    def test_the_leading_parse_admits_a_payload_and_refuses_a_wrong_lead(self):
        assert cs.leading_coordinates(np.ones((2, 3, 4)), (2, 3), owner="t").shape == (2, 3, 4)  # A ⊗ I_P, payload (4,)
        for bad in [(2, 1, 4), (1, 3), (3, 2), (2,)]:
            with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
                cs.leading_coordinates(np.ones(bad), (2, 3), owner="t")

    def test_a_wrong_shape_is_a_value_error(self):
        """A wrong shape is a wrong value of the right type (the user's ruling, 2026-10-30): not a ``TypeError``."""
        assert issubclass(cs.NotInCoordinateSpace, ValueError)
        assert not issubclass(cs.NotInCoordinateSpace, TypeError)

    def test_the_message_names_the_owner(self):
        """``retirement-audit`` B.21: a shared helper moves the raise one frame out; the owner tag keeps "this verb
        refused" from widening to "something refused"."""
        with pytest.raises(cs.NotInCoordinateSpace, match="DiscreteMeasure.integrate"):
            cs.coordinates(np.ones(2), (3,), owner="DiscreteMeasure.integrate")

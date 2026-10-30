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
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
from orpheus.numerics.basis.legendre_basis import LegendreBasis
from orpheus.numerics.frame import GalerkinFrame
from orpheus.numerics.functional import InnerProductFunctional
from orpheus.numerics.metric import DenseMetric, DiagonalMetric, FactoredMetric
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace, Member

pytestmark = pytest.mark.l0

_E = EnergyAxis.synthetic(2)
V = FunctionSpace.of_axes(_E, Axis("spatial", (3,), weights=np.array([1.0, 2.0, 3.0]), origin=Declared(BasisKind.NODAL)))
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


class TestTheSpaceVerbsParseOnEveryRun:
    """The identity leg plus the primitive. First red before U1b: each wrong call RETURNED under ``-O`` (the explorer's
    ``probe_O.log``, on production spaces; here on the synthetic ``V`` so the numerics layer can host them)."""

    def test_the_honest_values(self):
        assert V.inner_product(X, X) == 210.0
        np.testing.assert_array_equal(V.apply_metric(X), X * np.array([1.0, 2.0, 3.0]))

    @pytest.mark.parametrize("verb", ["inner_product", "norm", "flat", "integral"])
    def test_a_scalar_verb_refuses_a_broadcastable_non_member(self, verb):
        y = X[:, :1]  # (2, 1): the implicit-pullback class
        call = {
            "inner_product": lambda: V.inner_product(X, y),
            "norm": lambda: V.norm(y),
            "flat": lambda: V.flat(y),
            "integral": lambda: V.integral.evaluate(y),
        }[verb]
        with pytest.raises(cs.NotInCoordinateSpace, match=EXACT):
            call()

    @pytest.mark.parametrize("verb", ["apply_metric", "apply_inverse_metric"])
    @pytest.mark.parametrize("bad", [(2, 1), (1, 3)], ids=["pullback", "non-member-out"])
    def test_an_action_refuses_a_wrong_lead(self, verb, bad):
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            getattr(V, verb)(np.ones(bad))

    def test_an_action_still_carries_a_payload(self):
        """D2: the action A ⊗ I_P stands (the Riesz-leg and ``.H`` payload callers stay)."""
        out = V.apply_metric(np.ones((2, 3, 4)))
        np.testing.assert_array_equal(out, np.broadcast_to(np.array([1.0, 2.0, 3.0])[None, :, None], (2, 3, 4)))

    def test_the_space_verbs_have_no_debug_tier(self):
        """X1's execution evidence: no ``__debug__`` name and no ``assert`` in the space's two membership verbs."""
        import textwrap

        for verb in (FunctionSpace.coordinates_of, FunctionSpace.leading_coordinates_of):
            tree = ast.parse(textwrap.dedent(inspect.getsource(verb)))
            assert not any(isinstance(n, ast.Name) and n.id == "__debug__" for n in ast.walk(tree))
            assert not any(isinstance(n, ast.Assert) for n in ast.walk(tree))

    def test_a_member_is_built_through_the_same_primitive(self):
        """``Member.__post_init__``'s shape check is the primitive (one spelling): its message carries the new fragment
        AND the two the tree pins on it (``test_field.py``, ``test_reaction_rate_on_the_pose.py``,
        ``test_angular_boundary_flux.py``, ``test_wall_condition_algebra.py``)."""
        with pytest.raises(cs.NotInCoordinateSpace, match=EXACT) as err:
            Member(values=X[:, :1], space=V)
        assert "does not match space.shape" in str(err.value) and "values.shape" in str(err.value)


# ── objects below the space level parse their own domain (U1b step 3) ──────────────────────────────────────

_W = np.array([1.0, 2.0, 3.0])
_Q = Quadrature.gauss_legendre(n_ordinates=4)
_B = LegendreBasis(2)
_FR = GalerkinFrame(_B, _Q.measure)


class TestBelowTheSpaceEachObjectParsesItsDomain:
    """Each object that owns its domain shape parses it with the one primitive. First red before U1b: every wrong call
    RETURNED under ``-O`` (`[M]` ``scratch/membership_covector/u1b/probes/below_space_O.out``, the wrong answer's shape
    and sum recorded per row)."""

    def test_a_diagonal_block_parses_its_block(self):
        # [M] before U1b: (2, 1) against block (3,) at start 1 returned a (2, 3) array, sum 12 (the pullback)
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            DiagonalMetric(_W).apply_block(np.ones((2, 1)), start=1, block_shape=(3,), inverse=False)

    def test_a_dense_block_parses_its_block(self):
        # [M] before U1b: (4,) against a dim-2 block returned (4,), sum 8
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            DenseMetric(np.array([[2.0, 0.5], [0.5, 1.0]])).apply_block(np.ones(4), start=0, block_shape=(2,), inverse=False)

    def test_a_factored_metric_parses_its_concatenated_blocks(self):
        """A Euclidean (``None``) block applies nothing, so no block alone sees its extent: the product parses."""
        # [M] before U1b: apply on (2, 1) returned (2, 3); pairing with y (2, 1) returned 12
        fm = FactoredMetric((((2,), None), ((3,), DiagonalMetric(_W))))
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            fm.apply(np.ones((2, 1)))
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            FactoredMetric((((2,), DiagonalMetric(np.ones(2))), ((3,), None))).apply(np.ones((2, 1)))

    def test_a_legacy_diagonal_metric_is_checked_against_its_space_once(self):
        """The legacy ``inner_product_weights`` may be broadcast-shaped: ``(1, 4, 1)`` and the leading-aligned ``(3,)``,
        ``(3, 1)``, ``(1, 4)`` on a ``(3, 4, 2)`` space are legal (`[M]` each pairs ones with ones to 24.0). So the
        weights do not own a domain shape: the space parses each element, and the weights are checked against the
        space ONCE, when it is built. Before U1b, weights ``(2,)`` built and failed only at first use."""
        with pytest.raises(ValueError, match="cannot serve a space"):
            FunctionSpace(name="p", shape=(3, 4, 2), inner_product_weights=np.ones(2))
        for legal in [(1, 4, 1), (3,), (3, 1), (1, 4)]:
            space = FunctionSpace(name="p", shape=(3, 4, 2), inner_product_weights=np.ones(legal))
            assert space.inner_product(np.ones((3, 4, 2)), np.ones((3, 4, 2))) == 24.0

    def test_a_measure_parses_its_nodes(self):
        # [M] before U1b: (1, 4) on a 4-node rule returned a (4,) array, sum 8: a non-member out
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            _Q.measure.integrate(np.ones((1, 4)))
        np.testing.assert_array_equal(_Q.measure.integrate(np.ones((4, 2))), np.full(2, _Q.measure.weights.sum()))

    def test_an_inner_product_functional_parses_its_weight(self):
        # [M] before U1b: the (2,) carrier against the (2, 1) weight returned (1, 2), sum 0.3 (ERR-122's class)
        with pytest.raises(cs.NotInCoordinateSpace, match=EXACT):
            InnerProductFunctional(np.array([[0.1], [0.2]]), axis=0).evaluate(np.array([0.6, 0.4]))

    @pytest.mark.parametrize("verb", ["reconstruct", "analyze_transpose", "synthesize"])
    def test_a_basis_parses_its_coefficients_at_the_table_modes(self, verb):
        # [M] before U1b: (1, 5) coefficients for 3 modes returned (4, 5): one coefficient read for three
        args = (_FR.table, _Q.measure.weights) if verb == "analyze_transpose" else (_FR.table,)
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            getattr(_B, verb)(np.ones((1, 5)), *args)

    @pytest.mark.parametrize("verb", ["reconstruct_transpose", "analyze"])
    def test_a_basis_parses_its_values_at_the_table_nodes(self, verb):
        args = (_FR.table, _Q.measure.weights) if verb == "analyze" else (_FR.table,)
        with pytest.raises(cs.NotInCoordinateSpace, match=LEADING):
            getattr(_B, verb)(np.ones((1, 5)), *args)

    def test_every_basis_family_inherits_the_one_parse(self):
        """The five contractions are concrete on ``Basis`` and final; a family implements only the hooks
        (``_reconstruct`` …), so a new family cannot skip the parse."""
        from orpheus.numerics.basis.base import Basis

        def families(cls):
            for sub in cls.__subclasses__():
                yield sub
                yield from families(sub)

        verbs = ("synthesize", "analyze", "analyze_transpose", "reconstruct", "reconstruct_transpose")
        found = list(families(Basis))
        assert len(found) >= 6, found
        for cls in found:
            overridden = [v for v in verbs if v in vars(cls)]
            assert not overridden, (cls.__name__, overridden)

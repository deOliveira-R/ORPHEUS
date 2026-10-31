r"""D1: a plain binding admits a typed member of its end, by content; it refuses a member of another space of equal shape.

The user's ruling D1, 2026-10-30, plan
``space_membership_covector_gauge.md`` U1b). ``admit_array`` becomes the plain binding's call of its end's one parse
(``FunctionSpace.coordinates_of``): the typed leg compares the member's space by content identity, the bare leg is the
coordinate-space primitive; only the composite refusal (``BulkLift``) and the type refusal stay ``admit_array``'s own.

The return for a typed operand is the same bare array as for a bare one (the orchestrator, 2026-10-30): typed in,
bare out is a KNOWN RESIDUE owned by #633 (the role typing of arrows), asserted here only as ``array_equal``.

First reds before D1 (`[M]` at ``6aa26d00``):

* the typed member of the end is REFUSED with ``TypeError: ... the plain binding acts on bare arrays — pass
  ScalarFlux.values`` (both ``IsotropicFission`` and ``MultiplicationOperator`` on the 2g homogeneous pose);
* so the equal-shape foreign member is refused for the WRONG reason (the "pass .values" message), which the
  discrimination row below reds on.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.homogeneous.solver import HomogeneousProblem
from orpheus.numerics.axis import Axis
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.fields.scalar_flux import ScalarFlux

pytestmark = pytest.mark.l0

_OLD = "the plain binding acts on bare"        # the D1-retired refusal: must be ABSENT from every message below


def _pose():
    p = HomogeneousProblem(get_mixture("A", "2g"))
    energy, spatial = p.space.axes
    # Equal shape (2, 1), another spatial measure: the two spaces only content identity tells apart.
    foreign = FunctionSpace.of_axes(energy, Axis("spatial", (1,), weights=np.array([2.0]), origin=spatial.origin))
    assert foreign != p.space and foreign.shape == p.space.shape       # the discriminator, asserted in-test
    return p, foreign


@pytest.mark.parametrize("binding", ["production", "collision"])
class TestAPlainBindingAdmitsATypedMemberOfItsEnd:
    def test_the_typed_member_reads_the_bare_arithmetic(self, binding):
        """R1's pair (raw and formal build one object): ``array_equal``, the same multiply."""
        p, _ = _pose()
        col = np.array([[0.6], [0.4]])
        op = getattr(p, binding)
        np.testing.assert_array_equal(np.asarray(op.apply(ScalarFlux(values=col, space=p.space))), op.apply(col))

    def test_a_member_of_an_equal_shape_space_is_refused_by_content(self, binding):
        """DISCRIMINATION row (the wiring order): the new reason present, the retired one absent."""
        p, foreign = _pose()
        with pytest.raises((TypeError, ValueError)) as err:
            getattr(p, binding).apply(ScalarFlux(values=np.array([[0.6], [0.4]]), space=foreign))
        assert "is not a member of" in str(err.value)
        assert _OLD not in str(err.value)

    def test_a_bare_non_member_is_refused_by_the_one_parse(self, binding):
        """The ``(2,)`` array numpy would broadcast against ``(2, 1)`` (ERR-122's shape) is refused at the exact
        shape, a ``ValueError`` (the user's ruling: ``ValueError`` only)."""
        p, _ = _pose()
        with pytest.raises(ValueError, match="exact shape"):
            getattr(p, binding).apply(np.array([0.6, 0.4]))


def test_the_composite_refusal_stays_admit_arrays_own():
    """D1 keeps the BulkLift refusal (the composite action is the lift's). Reuse: ``test_bulk_lift.py``
    ``test_admit_array_refusals_name_the_lift_and_the_values`` loses its ``.values`` leg (D1 flips it) and keeps
    the ``BulkLift`` and type legs; this row is its positive companion on the homogeneous pose."""
    from orpheus.transport.full_field import FullField  # noqa: F401  (the composite type the refusal names)
    p, _ = _pose()
    with pytest.raises(TypeError, match="IsotropicFission"):
        p.production.apply(object())

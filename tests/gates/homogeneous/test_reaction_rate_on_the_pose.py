r"""The reaction rate on the homogeneous 0-D pose: where the 0.225 defect lived (ERR-122), and ⟨1, φ⟩.

The space-membership unit U2c (plan ``.claude/plans/space_membership_covector_gauge.md``). The defect: the 2g Mixture-A
production rate of a flat ``(2,)`` flux on the ``(2, 1)`` pose read 0.225 for the exact 0.095 (an outer product
contracted; the layer-DAG Phase 2 spec §5c). Both paths refuse the input on every run: the formal path at the scalar
flux's construction, the raw path at the space's parse (U1b; until then it was refused under ``__debug__`` only).
These rows moved here from ``tests/gates/transport/test_integrated_reaction_rate.py``, whose layer may not import
``orpheus.homogeneous`` (``tests/gates/test_test_layers.py``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.homogeneous.solver import HomogeneousProblem
from orpheus.transport.fields.scalar_flux import ScalarFlux

pytestmark = pytest.mark.foundation

class TestTheReactionRateOnThePose:
    def test_a_flat_carrier_on_the_homogeneous_column_space_is_REFUSED(self):
        """The defect this campaign closes (`[M]` the layer-DAG Phase 2 spec §5c): the 2g Mixture-A production rate of
        a flat ``(2,)`` flux on the ``(2, 1)`` pose read 0.225 for the exact 0.095 (an outer product contracted). On the
        raw path the shape is refused by the space's parse on every run; the formal path is the next row."""
        problem = HomogeneousProblem(get_mixture("A", "2g"))
        with pytest.raises(ValueError, match="exact shape"):
            problem.production_rate.evaluate(np.array([0.6, 0.4]))

    @pytest.mark.catches("ERR-122")
    def test_a_typed_flat_carrier_cannot_be_built(self):
        """On every run, ``-O`` included: a ``(2,)`` scalar flux on the ``(2, 1)`` pose is refused at its own construction,
        so the formal path cannot spell the 0.225 input."""
        problem = HomogeneousProblem(get_mixture("A", "2g"))
        with pytest.raises(ValueError, match="does not match space.shape"):
            ScalarFlux(values=np.array([0.6, 0.4]), space=problem.space)

    def test_the_homogeneous_total_flux_is_the_integral_bit_for_bit(self):
        """⟨1, φ⟩ on the counting 0-D pose: the integration co-vector and the space's own pairing with the unit are one
        arithmetic (the metric applied to the unit, then one sum), so ``==``, not a tolerance."""
        space = HomogeneousProblem(get_mixture("A", "2g")).space
        phi = np.array([[0.7], [0.3]])
        assert space.integral.evaluate(ScalarFlux(values=phi, space=space)) == space.inner_product(space.unit, phi)

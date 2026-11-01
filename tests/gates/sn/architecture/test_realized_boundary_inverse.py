"""MOVED from tests/gates/numerics/test_inverse_universal.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import numpy as np
import pytest
from orpheus.numerics.operator import NotInvertible


pytestmark = pytest.mark.foundation






















def test_realized_boundary_forwards_inverse():
    """The NINTH advertiser: a realized boundary operator forwards ``.inverse()`` (§12.5).

    Wave 2 S2b re-pose: the compatibility shim (``_BoundBoundaryOperator``) is retired and the typed realized operators of
    ``orpheus.sn.boundary.realized_operators`` are the advertisers. The mirror face (an ordinate permutation, invertible) forwards
    its inner permutation's inverse; the vacuum face (a zero map, a projection) propagates the inner operator's own refusal.
    """
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.boundary.realized_operators import ScaledOrdinatePermutation, ZeroMap
    from tests.gates.sn._test_helpers import realized_slab_faces

    quad = Quadrature.gauss_legendre(n_ordinates=8)
    vacuum, mirror = realized_slab_faces(quad, "vacuum", "reflective")
    if not isinstance(mirror, ScaledOrdinatePermutation) or not isinstance(vacuum, ZeroMap):
        pytest.fail(f"expected a mirror permutation and a zero map, got {type(mirror).__name__} and {type(vacuum).__name__}")
    x = np.random.default_rng(7).standard_normal(mirror.inner.codomain.shape)
    np.testing.assert_array_equal(mirror.inverse().apply(x), mirror.inner.inverse().apply(x))
    # carve P4: the shim's solve forward retired with the surface (the promised deletion), and the typed operators never had one:
    # solving through a realized boundary IS .inverse().apply.
    assert not hasattr(ScaledOrdinatePermutation, "solve")
    with pytest.raises(NotInvertible):
        vacuum.inverse()



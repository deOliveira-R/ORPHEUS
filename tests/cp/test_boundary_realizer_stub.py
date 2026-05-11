"""Stub-realizer tests for CP (Wave 5 / C5.4).

See :mod:`tests.moc.test_boundary_realizer_stub` for invariant rationale.
"""

import pytest

from orpheus.geometry.boundary import BoundaryRealizerRegistry
from orpheus.cp.boundary_realizer import CPBoundaryRealizer


@pytest.mark.l0
class TestStubRegistered:
    def test_method_name_attribute(self):
        assert CPBoundaryRealizer.method_name == "CP"

    def test_registered_in_registry(self):
        assert BoundaryRealizerRegistry.get("CP") is CPBoundaryRealizer

    def test_realize_raises_NotImplementedError(self):
        realizer = CPBoundaryRealizer()
        with pytest.raises(NotImplementedError, match="not yet implemented"):
            realizer.realize(law=None, method_space=None)

"""Stub-realizer tests for MC (Wave 5 / C5.4).

See :mod:`tests.moc.test_boundary_realizer_stub` for invariant rationale.
"""

import pytest

from orpheus.geometry.boundary import BoundaryRealizerRegistry
from orpheus.mc.boundary_realizer import MCBoundaryRealizer


@pytest.mark.l0
class TestStubRegistered:
    def test_method_name_attribute(self):
        assert MCBoundaryRealizer.method_name == "MC"

    def test_registered_in_registry(self):
        assert BoundaryRealizerRegistry.get("MC") is MCBoundaryRealizer

    def test_realize_raises_NotImplementedError(self):
        realizer = MCBoundaryRealizer()
        with pytest.raises(NotImplementedError, match="not yet implemented"):
            realizer.realize(law=None, method_space=None)

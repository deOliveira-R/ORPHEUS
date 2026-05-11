"""Stub-realizer tests for diffusion (Wave 5 / C5.4).

See :mod:`tests.moc.test_boundary_realizer_stub` for invariant rationale.
"""

import pytest

from orpheus.geometry.boundary import BoundaryRealizerRegistry
from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer


@pytest.mark.l0
class TestStubRegistered:
    def test_method_name_attribute(self):
        assert DiffusionBoundaryRealizer.method_name == "diffusion"

    def test_registered_in_registry(self):
        assert (
            BoundaryRealizerRegistry.get("diffusion")
            is DiffusionBoundaryRealizer
        )

    def test_realize_raises_NotImplementedError(self):
        realizer = DiffusionBoundaryRealizer()
        with pytest.raises(NotImplementedError, match="not yet implemented"):
            realizer.realize(law=None, method_space=None)

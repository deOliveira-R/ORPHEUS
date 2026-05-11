"""Stub-realizer tests for MoC (Wave 5 / C5.4).

Pins three architectural invariants of the stub realizer:
1. ``method_name`` attribute matches the registry key.
2. The stub self-registers under ``"MoC"`` at import time
   (verified by :meth:`BoundaryRealizerRegistry.get`).
3. :meth:`realize` raises :class:`NotImplementedError` with the
   "not yet implemented" sentinel string so the architecture is
   wired end-to-end but the math is deferred.
"""

import pytest

from orpheus.geometry.boundary import BoundaryRealizerRegistry
from orpheus.moc.boundary_realizer import MoCBoundaryRealizer


@pytest.mark.l0
class TestStubRegistered:
    def test_method_name_attribute(self):
        assert MoCBoundaryRealizer.method_name == "MoC"

    def test_registered_in_registry(self):
        assert BoundaryRealizerRegistry.get("MoC") is MoCBoundaryRealizer

    def test_realize_raises_NotImplementedError(self):
        realizer = MoCBoundaryRealizer()
        with pytest.raises(NotImplementedError, match="not yet implemented"):
            realizer.realize(law=None, method_space=None)

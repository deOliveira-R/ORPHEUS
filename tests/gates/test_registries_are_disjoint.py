"""MOVED from tests/gates/numerics/test_registry_mixin.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import pytest
from orpheus.numerics.registry import RegistryMixin












@pytest.mark.foundation
def test_registries_disjoint_across_codebase():
    """AngularClosureBase and DiscretizationSchemeBase registries don't share keys
    just because they both use RegistryMixin. (The third root, ``BoundaryTraceLaw``, retired in WF2.4.)"""
    from orpheus.sn.angular.closure import AngularClosureBase
    from orpheus.transport.spatial.scheme import DiscretizationSchemeBase

    # The two registries are distinct dicts.
    assert AngularClosureBase.registry is not DiscretizationSchemeBase.registry
    assert AngularClosureBase.registry and DiscretizationSchemeBase.registry
    # And they hold their own subclasses, not each other's.
    for cls in AngularClosureBase.registry.values():
        assert cls not in DiscretizationSchemeBase.registry.values()



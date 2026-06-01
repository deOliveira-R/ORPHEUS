"""Capability marker for tests under ``tests/sn/primitives/``.

The directory IS the capability (single source of truth). Every test
collected directly in this directory is stamped ``cap("primitives")``;
see ``stamp_capability_marker`` for the rationale.
"""
from tests.sn._test_helpers import stamp_capability_marker


def pytest_collection_modifyitems(items):
    stamp_capability_marker(items, __file__, "primitives")

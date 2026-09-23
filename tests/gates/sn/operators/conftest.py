"""Capability marker for tests under ``tests/gates/sn/operators/``.

The directory IS the capability (single source of truth). Every test
collected directly in this directory is stamped ``cap("operators")``;
see ``stamp_capability_marker`` for the rationale.
"""
from tests.gates.sn._test_helpers import stamp_capability_marker


def pytest_collection_modifyitems(items):
    stamp_capability_marker(items, __file__, "operators")

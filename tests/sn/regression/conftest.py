"""Capability marker for the frozen-snapshot regression gate.

``tests/sn/regression/`` is the canonical snapshot + generator store
(referenced sn-root-relative by several sweep tests), so it stays in
place rather than moving under solve/. Its one TEST (test_dd_regression)
is a solve-tier drift gate; stamp it ``cap("solve")``.
"""
from tests.sn._test_helpers import stamp_capability_marker


def pytest_collection_modifyitems(items):
    stamp_capability_marker(items, __file__, "solve")

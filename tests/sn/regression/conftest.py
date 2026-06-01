"""Capability marker + DriftWarning registration for the regression gate.

``tests/sn/regression/`` is the canonical snapshot + generator store
(referenced sn-root-relative by several sweep tests), so it stays in
place rather than moving under solve/. Its one TEST (test_dd_regression)
is a solve-tier drift gate; stamp it ``cap("solve")``.

``DriftWarning`` is the bit-identity tripwire (see
:mod:`._regression_assert`). We register it as ``always`` so pytest
surfaces it in the run-end warnings summary (deduped, visible) rather
than letting it be swallowed by the default once-per-location filter.
The strict bit-identity gate for a "pure refactor, zero numerical
change" PR is::

    pytest tests/sn/regression/ -W error::tests.sn.regression._regression_assert.DriftWarning

which escalates ANY sub-tolerance bit drift to a hard failure.
"""
import warnings

from tests.sn._test_helpers import stamp_capability_marker
from tests.sn.regression._regression_assert import DriftWarning


def pytest_configure(config):
    # Surface every DriftWarning in the summary (not just the first per
    # location) so a refactor that moves several cases is fully reported.
    warnings.simplefilter("always", DriftWarning)
    config.addinivalue_line(
        "filterwarnings", "always::tests.sn.regression._regression_assert.DriftWarning"
    )


def pytest_collection_modifyitems(items):
    stamp_capability_marker(items, __file__, "solve")

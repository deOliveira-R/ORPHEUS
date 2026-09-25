"""The test-side marks of the withdrawals in force — one per withdrawal.

A withdrawal (:class:`orpheus.derivations.common.withdrawal.Withdrawal`) is
declared once, next to the generators it locks. Its test-side marker is
minted HERE from that same value, so every test that consumes a withdrawn
generator carries ``@PEIERLS_NYSTROM_WITHDRAWN`` (or
``pytestmark = PEIERLS_NYSTROM_WITHDRAWN``) and the reason sentence is
spelled once, not at every site. ``tests/conftest.py`` resolves and parses
the marker like any other ``@pytest.mark.withdrawn(reason, issue=N)``.
"""

from __future__ import annotations

import pytest

from orpheus.derivations.continuous.peierls_nystrom import PEIERLS_NYSTROM_WITHDRAWAL

__all__ = ["PEIERLS_NYSTROM_WITHDRAWN"]

#: The #506 withdrawal of the Peierls Nyström solver half, as a test mark.
PEIERLS_NYSTROM_WITHDRAWN = pytest.mark.withdrawn(
    PEIERLS_NYSTROM_WITHDRAWAL.reason, issue=PEIERLS_NYSTROM_WITHDRAWAL.issue
)

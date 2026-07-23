"""ORPHEUS V&V test harness: markers, registry, shared helpers, audit.

See ``docs/theory/verification/harness.rst`` for the full contributor
guide.

Public API
----------
TestMetadata / TEST_REGISTRY
    In-memory V&V metadata populated by the conftest collection hook.
    One entry per pytest item keyed by ``item.nodeid``. Consumed by the
    audit CLI (``python -m tests._harness.audit``) and, via a Sphinx
    generator, by the ``docs/theory/verification/`` matrix page.

Tests tag levels with raw ``pytest.mark.*`` decorators (``l0``..``l3``
/ ``foundation`` + ``verifies`` / ``catches``) — the one convention the
whole tree uses. A ``verify.lN(...)``/``vv_cases(...)`` sugar layer
existed here until 2026-07 but was retired unused (zero consumers).
"""

from __future__ import annotations

from tests._harness.registry import TEST_REGISTRY, TestMetadata

__all__ = ["TEST_REGISTRY", "TestMetadata"]

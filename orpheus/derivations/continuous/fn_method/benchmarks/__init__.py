"""Deprecation shim — re-exports from :mod:`sood_registry`.

This package moved to
:mod:`orpheus.derivations.continuous.sood_registry`. New code should
import from there directly. Old import paths still work via
:mod:`.la13511` (which emits a DeprecationWarning on import).
"""
from __future__ import annotations

from .la13511 import (  # noqa: F401
    ALL_FIRST_SLICE,
    LA13511_CASES,
    La13511Case,
    La13511Truth,
    MeshTemplate,
    PU_2_0_IN,
    PUA_1_0_IN,
    UA_1_0_CY_STUB,
    UA_1_0_SL_STUB,
    UA_1_0_SP_STUB,
)

__all__ = [
    "La13511Case",
    "La13511Truth",
    "MeshTemplate",
    "PUA_1_0_IN",
    "PU_2_0_IN",
    "UA_1_0_SL_STUB",
    "UA_1_0_CY_STUB",
    "UA_1_0_SP_STUB",
    "ALL_FIRST_SLICE",
    "LA13511_CASES",
]

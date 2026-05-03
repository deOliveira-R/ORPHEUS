"""Deprecation shim — re-exports from :mod:`sood_registry.la13511`.

The Sood-family benchmark catalogue moved out of
``fn_method/benchmarks/`` and into the method-agnostic
:mod:`orpheus.derivations.continuous.sood_registry` package. The new
location keeps method-specific code (semi-analytical reference solvers
in ``fn_method/``, Variant α in ``trajectory_resolvent/``) separate
from method-agnostic case configurations (XS, geometry, truth).

Historical importers (e.g.
``from orpheus.derivations.continuous.fn_method.benchmarks.la13511
import PUA_1_0_IN``) still work — they receive a deprecation warning
on import. New code should import from
:mod:`orpheus.derivations.continuous.sood_registry` directly.
"""
from __future__ import annotations

import warnings

warnings.warn(
    "orpheus.derivations.continuous.fn_method.benchmarks.la13511 has "
    "moved to orpheus.derivations.continuous.sood_registry.la13511 "
    "(method-agnostic case registry). The old import path is a "
    "deprecation shim and will be removed in a future release.",
    DeprecationWarning,
    stacklevel=2,
)

from orpheus.derivations.continuous.sood_registry.la13511 import (  # noqa: E402, F401
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

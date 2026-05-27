r"""DEPRECATED re-export shim — TraceSpace moved to
:mod:`orpheus.numerics.spaces.trace_space`.

This shim exists for ONE merge cycle to give in-flight callers a window
to update their imports (per ``feedback_aggressive_retirement``). It will
be deleted in Depth B step D-K. New code MUST import from the canonical
:mod:`orpheus.numerics.spaces.trace_space` location.

Move rationale: Depth B step D-C aligns the L1 spaces layout with the
plan's §3.1 layer assignments — every :class:`FunctionSpace` subclass
lives under ``orpheus/numerics/spaces/``. See
``.claude/plans/depth_b_field_on_function_space.md`` §6 step D-C.
"""

from __future__ import annotations

import warnings as _warnings

from orpheus.numerics.spaces.trace_space import (  # noqa: F401
    InflowTraceSpace,
    OutflowTraceSpace,
    TraceSpace,
)

_warnings.warn(
    "orpheus.numerics.trace_space is deprecated; "
    "import from orpheus.numerics.spaces.trace_space instead "
    "(shim retires in Depth B step D-K).",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["InflowTraceSpace", "OutflowTraceSpace", "TraceSpace"]

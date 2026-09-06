r"""ONE spelling of the zero-trace composite operand the angular gains admit.

CS4c step 5 (*each binding acts through the body its ends select*): every
angular gain — :class:`~orpheus.transport.operators.scattering.ScatteringOperator`,
:class:`~orpheus.transport.operators.n2n.N2NOperator`,
:class:`~orpheus.transport.operators.fission.FissionOperator` — is bound on the
composite ``bulk ⊕ trace`` :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
and admits exactly the :class:`~orpheus.transport.full_field.FullField` whose
interior rides its bound end's interior
(:func:`~orpheus.transport.operators.lift.admit_composite`). A bare bulk field
handed to ``apply``/``apply_transpose`` is a ``TypeError`` naming the operator.

A gate whose subject is the BULK action therefore rides its bulk field in a
composite whose trace is ZERO — which is also exactly what the lift emits back
(a collision gain is volumetric: extension-by-zero on the trace,
:func:`~orpheus.transport.operators.lift.lift_bulk_action`). That wrap is one
concept, so it is spelled once here instead of at each of the ~40 call sites the
step-5 re-key touched (``coding-elegance`` Pattern 2). The trace VALUES are
irrelevant to a volumetric gain's output — zero is the honest choice because it
is what the operator emits, so ``interior`` in / ``interior`` out is a
round-trip through the composite that adds no content.

**The blocks are read off the OPERATOR's own bound end**, never off a mesh
handed in beside it: the operand's space is the binding's business (SPACE
FIRST), and an LD binding's interior carries a trailing spatial-moment axis the
mesh's ``angular_bulk_space`` does not. That also makes these helpers work
unchanged for a moment-domain sibling
(:meth:`~orpheus.transport.operators.angular_lift.AngularLift.on_moment_domain`),
whose domain interior is the moment composite.

⚠ These helpers are for the operand's *shape*, never for its *claim*: a row
that means to gate the trace block must build its own non-zero trace and assert
on ``out.boundary``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import AngularSourceSink

if TYPE_CHECKING:  # pragma: no cover - typing only
    from orpheus.numerics.field import Field
    from orpheus.numerics.space import FunctionSpace

__all__ = [
    "bulk_apply",
    "bulk_apply_transpose",
    "transpose_values",
    "zero_trace_composite",
]


def _trace_space(space: "FunctionSpace", owner: str) -> "FunctionSpace":
    trace = getattr(space, "trace_space", None)
    if trace is None:
        raise TypeError(
            f"{owner}: {space!r} carries no trace block — the composite "
            f"operand needs one to be built.",
        )
    return trace


def zero_trace_composite(bulk: "Field", trace_space: "FunctionSpace") -> FullField:
    """``bulk`` in the composite whose trace is the zero angular flux."""
    return FullField(
        interior=bulk,  # type: ignore[arg-type]
        boundary=AngularBoundaryFlux.zeros(trace_space),
    )


def bulk_apply(op, bulk: "Field") -> "Field":
    """``op.apply`` on the bulk block — the composite wrap, then ``.interior``."""
    trace = _trace_space(op.domain, f"{type(op).__name__}.domain")
    return op.apply(zero_trace_composite(bulk, trace)).interior


def bulk_apply_transpose(op, bulk: "Field") -> "Field":
    """``op.apply_transpose`` on the bulk cotangent — wrap, then ``.interior``."""
    trace = _trace_space(op.codomain, f"{type(op).__name__}.codomain")
    return op.apply_transpose(zero_trace_composite(bulk, trace)).interior


def transpose_values(op, values) -> np.ndarray:
    r"""``opᵀ`` on a RAW cotangent array — the bucket-(c) wrap, values in / out.

    The cotangent of a composite-bound gain is a composite too: an
    :class:`~orpheus.transport.source_sinks.AngularSourceSink` on the operator's
    own CODOMAIN interior plus the zero trace. Returns the emitted interior's
    raw values, so a caller that only wants the arithmetic keeps its arithmetic.
    """
    interior = AngularSourceSink(
        values=np.asarray(values), space=op.codomain.interior_space,
    )
    return np.asarray(bulk_apply_transpose(op, interior).values)

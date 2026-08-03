r"""L2 transport vocabulary — method-agnostic primitives.

This package is the L2 layer in the ORPHEUS layer contract (per plan
§P3.0 and Depth B plan §3.1). It carries transport concepts that are
shared across SN, CP, MoC, diffusion, Pn — fields, sources, problems,
solver protocols — without binding to any one method's discretisation
machinery.

Layer rules (enforced by ``tests/test_layer_imports.py``):

* ``transport`` imports from ``numerics`` (L1), ``geometry`` /
  ``data`` (input), and standard library / numpy / scipy. It does
  NOT import from any L3 package (``sn``, ``cp``, ``moc``,
  ``diffusion``, ``mc``, ``kinetics``, ``pn``, ``homogeneous``,
  ``fuel``, ``thermal_hydraulics``).
* L3 packages may import freely from ``transport`` — the L2 types
  are the canonical contracts the methods implement.

Today's contents:

* :mod:`orpheus.transport.fields` — typed flux fields, each inheriting
  from :class:`~orpheus.numerics.field.Field`: :class:`ScalarFlux`
  (Depth B step D-D), HarmonicMomentFlux (D-E), AngularBoundaryFlux
  (D-G), AngularFlux (D-H).
* :mod:`orpheus.transport.source_sinks` — ScalarSourceSink /
  AngularSourceSink and the boundary / moment role leaves (D-F).

Future contents (deferred to later parent-plan steps; NOT built — the
names below are plan targets, not importable modules):

* ``orpheus.transport.problems`` — Problem ABCs (P3.4): Criticality,
  FixedSource, AlphaEigen, InitialValue.

References
----------

* ``.claude/plans/depth_b_field_on_function_space.md`` §3.1
  (layer assignments).
* ``.claude/plans/moment_space_and_layering_plan.md`` §P3.0
  (layer criterion).
"""

from __future__ import annotations

from orpheus.transport.fields import (
    CrossSectionField,
    HarmonicMomentFlux,
    ScalarFlux,
)
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.integral_kernel_operator import IntegralKernelOperator
from orpheus.transport.timed_full_field import TimedFullField

__all__ = [
    "CrossSectionField",
    "FullField",
    "HarmonicMomentFlux",
    "IntegralKernelOperator",
    "ScalarFlux",
    "TimedFullField",
]

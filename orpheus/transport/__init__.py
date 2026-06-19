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

* :mod:`orpheus.transport.fields` — typed flux / source / residual
  fields, each inheriting from :class:`~orpheus.numerics.field.Field`.
  Currently houses :class:`ScalarFlux` (Depth B step D-D). Later
  steps add HarmonicMomentField (D-E), Sources (D-F), BoundaryFlux
  (D-G), AngularFlux (D-H).

Future contents (deferred to later Depth B / parent-plan steps):

* :mod:`orpheus.transport.source_sinks` — ScalarSourceSink / AngularSourceSink
  (D-F).
* :mod:`orpheus.transport.problems` — Problem ABCs (P3.4): Criticality,
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
    HarmonicMomentField,
    ScalarFlux,
)
from orpheus.transport.state import TransportState

__all__ = [
    "CrossSectionField",
    "HarmonicMomentField",
    "ScalarFlux",
    "TransportState",
]

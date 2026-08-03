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

Today's contents — the COMPLETE shipped set. Each role below is that
member's own one-line self-description, so this list is checkable
against the tree rather than a curated impression of it. (It listed 2
of 13 until 2026-08-03; an inventory that silently goes stale reads as
"this is all there is" and hides whole subsystems from a fresh reader.)

Carriers and their algebra:

* :mod:`orpheus.transport.fields` — typed flux fields, each inheriting
  from :class:`~orpheus.numerics.field.Field`.
* :mod:`orpheus.transport.source_sinks` — typed source/sink fields.
* :mod:`orpheus.transport.residuals` — typed transport residuals.
* :mod:`orpheus.transport.displacements` — typed flux displacements: the
  difference vector space :math:`V` the fields form a torsor over.
* :mod:`orpheus.transport.full_field` — the composite carrier, a generic
  ``interior ⊕ boundary`` block field.
* :mod:`orpheus.transport.timed_full_field` — the timed (history-bearing)
  full field.
* :mod:`orpheus.transport.radial_characteristic_field` — System B as an
  independent ``Composite[interior ⊕ boundary]``.

Operators, functionals and frames:

* :mod:`orpheus.transport.operators` — the §5.6 cross-method operator
  algebra (collision, fission, scattering). Import its members from the
  package, not from the defining submodule.
* :mod:`orpheus.transport.reaction_rate_functional` — the §5.6 functional
  :math:`\langle \Sigma_x, \cdot\rangle`.
* :mod:`orpheus.transport.frames` — carrier-typed specializations of the
  numerics frames.

Method-facing structure:

* :mod:`orpheus.transport.mesh` — the method-agnostic geometry/materials
  cluster.
* :mod:`orpheus.transport.method` — ``TransportMethod``, the Protocol over
  the method-mesh layer.
* :mod:`orpheus.transport.spatial` — method-generic spatial discretization
  schemes (the per-cell closure layer).

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

r"""The q½ starting-direction source block — the *source* role leaf of
:class:`~orpheus.transport.fields._bases.StartingDirectionField`.

The typed per-level source/sink block driving the starting-direction
equation (Hébert §3.9.4, Eq. 3.432: closed :math:`\mu = \mu_{\rm start}`
streaming :math:`+\;\sigma\psi_{1/2} = q_{1/2}`) — the third block of an
augmented SOURCE composite (a ``q_ext``, an operator ``.apply`` output)
on seed-carrying meshes (#282 route (a); presence per R12a).

Cells vs corner semantics (units note)
======================================

* The **cells** blocks ``(ng, nx)`` carry the within-group emission
  :math:`\bar q_{1/2}` folded to the starting direction (ruling R14:
  the full :math:`(-1)^\ell` Legendre fold) — a volumetric rate
  density, :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`
  (``1/(cm³·s·sr)``), exactly like the bulk
  :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`.
* The **corner** slots ``(ng,)`` are trace-like: the inflow corner
  carries the prescribed r = R inflow *datum* (an angular-flux value —
  the affine-BC ``q`` of the seed's identity row, mirroring the
  boundary trace's "on the trace a source does NOT pick up the
  volumetric cm⁻³" convention), and the outflow corner is the defect
  row's source slot. The declared :attr:`UNITS` is the dominant
  cells-block identity; the corner deviation is this documented note
  (the same per-slot exception the angular trace documents). UNITS is
  metadata — the arithmetic gate is class identity.

All storage, validation, algebra, ``cells``/``corner`` views, and the
``zeros_on`` / ``from_mesh`` factories are inherited from
:class:`StartingDirectionField`. Plain vector-space algebra (source
sums are closed) — no role mixin, like every SourceSink leaf.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import StartingDirectionField

__all__ = ["StartingDirectionSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class StartingDirectionSourceSink(StartingDirectionField):
    r"""L2 starting-direction source/sink — the q½ block of an augmented source composite.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)``.
    space : StartingDirectionSpace
        The R12a-keyed space (canonically
        ``mesh.starting_direction_space``).
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).

    Notes
    -----
    A thin role leaf — the class identity is what keeps seed source,
    flux, and displacement arithmetic from silently mixing (all three
    share the SAME ``mesh.starting_direction_space``, so it is the
    class gate, not the space gate, that rejects cross-role sums).
    """

    #: Dimensional identity (View-G): the cells blocks carry the folded
    #: volumetric emission rate ``1/(cm³·s·sr)`` —
    #: :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`, shared with
    #: ``AngularSourceSink``. The corner slots deviate (trace-like flux
    #: datum) — see the module docstring's units note. Metadata, not the
    #: arithmetic gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS

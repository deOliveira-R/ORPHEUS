r"""Isotropic source field :math:`Q(\vec r, g)`.

The L2 typed wrapper for an isotropic external + secondary source
density. Same storage shape as
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux` (``(ng, nx,
ny)``), but a structurally distinct physical quantity (units of
reaction-rate density / source density, not flux density).

Migration status (Depth B step D-F → B.1 → B.2 → B.5)
======================================================

Moved (D-F) from ``orpheus.sn.sources.IsotropicSource`` to here;
re-parented (B.1) onto
:class:`~orpheus.transport.fields._bases.ScalarField` (the storage-base
dedup); renamed ``IsotropicSource`` → ``ScalarSource`` (B.2, hard
rename, no shim) to complete the ``{Angular, Scalar} × {Flux, Source,
Residual}`` role grid; then ``ScalarSource`` → ``ScalarSourceSink``
(B.5, source/sink role rename — the leaf holds both production *sources*
and operator-loss *sinks*; see :mod:`orpheus.transport.fields`). The cross-class arithmetic
with :class:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink`
is **PRESERVED** — the refined Issue #207 principle (recorded in the
conversation 2026-05-27) permits cross-class dunders when there is a
canonical subspace-containment relation between the operands' spaces.
:class:`ScalarSourceSink` lives in the subspace of
:class:`AngularSourceSink` where every ordinate carries the same
value; the canonical injection ``iso → 1 ⊗ iso`` (broadcast across
the Ω axis) is the structural map the dunder applies during the
operation.

The math reads cleanly as

.. math::

    Q_{\text{total}} = Q_{\text{iso}} + Q_{\text{aniso}}

with the sum written once, in the class of the LARGER space
(:class:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink`
declares ``CONTAINS = (ScalarSourceSink,)``, and ``BulkField.__add__``
defers ``iso + aniso`` to its reflected add):

.. code-block:: python

    # Within-class: returns ScalarSourceSink (Field's inherited algebra).
    Q_total = iso_a + iso_b

    # Cross-class with containment: returns the LARGER type
    # (AngularSourceSink), broadcast inside the dunder.
    Q_total = iso + aniso          # → AngularSourceSink
    Q_total = aniso + iso          # commutative; same result

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

Reaction-rate density per cell per group:
:math:`[1/(\mathrm{cm^3 \cdot s})]`
(:data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`; eV-free, binned-
energy convention) for the isotropic form, vs
:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr})]` (``ANGULAR_RATE_UNITS``) for
the per-ordinate form. The dimensionality differs by ``1/sr``; physically
the iso → per-ordinate broadcast also requires a ``/(4π sr)``
normalisation. The dunder does NOT apply that conversion (the
broadcast is the pure structural injection); callers that need
Pattern 7 producer-side normalisation apply ``/sum_w`` explicitly
BEFORE the dunder (e.g. ``(iso / sum_w) + aniso``) — or use the
named factory :meth:`AngularSourceSink.from_isotropic` which bakes
``/sum_w`` in. The choice between caller-applied ``/sum_w`` and the
named factory is a Pattern 7 discipline concern, not a dunder concern.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, cast

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.units import SCALAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import ScalarField
from orpheus.transport.fields.scalar_flux import ScalarFlux



__all__ = ["ScalarSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarSourceSink(ScalarField, flux=ScalarFlux):
    r"""Isotropic source field :math:`Q(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(ng, nx, ny)`` in the principled layout
        (Issue #196 PR-INDEX-7).
    space : FunctionSpace
        The function space — the carrier's cached ``mesh.bulk_space``
        (CS4b S5).

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`
    (dunders ``+``, ``-``, unary ``-``, scalar ``*``, scalar ``/``,
    plus diagnostics ``linf``, ``l2``, ``inner_product``, ``copy``).
    Field's Layer 1 class-identity gate rejects cross-class arithmetic
    with :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
    (same shape, different physical kind), with
    :class:`AngularSourceSink` (different shape AND different physical
    kind), and with any other Field subclass. Same-class arithmetic is
    closed.

    The cross-class ``iso + aniso`` dunder of the pre-D-F class is
    RETIRED — see module docstring for the named-composition migration.
    """

    #: Dimensional identity (View-G, B.4): scalar rate density
    #: ``1/(cm³·s)`` — :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`,
    #: shared with ``ScalarResidual`` (same units, different role → the
    #: gate is class identity). Metadata, not the gate.
    UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS

    # ── Algebra: the containment sum ``iso + aniso`` is answered by
    # :class:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink`,
    # the class of the containing space (its ``CONTAINS`` names this class;
    # ``BulkField.__add__`` defers to it). This class names no angular type.

    # ── Conversion — RETIRED (CS4b S4) ───────────────────────────────
    #
    # ``as_per_ordinate`` (broadcast this iso source to per-ordinate,
    # WITHOUT the /Σw normalisation) retired with the field's ``mesh``
    # binding: the broadcast target needs the ANGULAR axis — carrier
    # knowledge a scalar field's space cannot carry — and the verb had
    # zero production callers. The surviving spellings of the same
    # injection: the containment dunder (``iso + AngularSourceSink`` —
    # the pure structural broadcast, above) and
    # :meth:`AngularSourceSink.from_isotropic` (the Pattern-7
    # producer-normalised factory, which takes the carrier as an
    # argument).

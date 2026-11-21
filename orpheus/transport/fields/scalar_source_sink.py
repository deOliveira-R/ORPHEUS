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
and operator-loss *sinks*; see :mod:`orpheus.transport.fields`). The algebra is
:class:`~orpheus.numerics.field.Field`'s, same-class and closed. A scalar
source enters the per-ordinate transport equation through
:meth:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink.from_isotropic`
(the angular section: divide by the mass, then scatter across the ordinates)
and then a within-class add; ``ScalarSourceSink + AngularSourceSink`` is
refused by ``Field._check_partner`` (the units differ by ``1/sr``).

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

Reaction-rate density per cell per group:
:math:`[1/(\mathrm{cm^3 \cdot s})]`
(:data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`; eV-free, binned-
energy convention) for the isotropic form, vs
:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr})]` (``ANGULAR_RATE_UNITS``) for
the per-ordinate form. The dimensionality differs by ``1/sr``; physically
the iso → per-ordinate broadcast also requires a ``/(4π sr)``
normalisation, which is the angular section's divisor, the correctly
rounded mass of the stored weights; it is applied by
:meth:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink.from_isotropic`.
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
    """

    #: Dimensional identity (View-G, B.4): scalar rate density
    #: ``1/(cm³·s)`` — :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`,
    #: shared with ``ScalarResidual`` (same units, different role → the
    #: gate is class identity). Metadata, not the gate.
    UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS

    # ── Conversion — RETIRED (CS4b S4) ───────────────────────────────
    #
    # ``as_per_ordinate`` (broadcast this iso source to per-ordinate,
    # WITHOUT the /Σw normalisation) retired with the field's ``mesh``
    # binding: the broadcast target needs the ANGULAR axis — carrier
    # knowledge a scalar field's space cannot carry — and the verb had
    # zero production callers. The surviving spelling of the same
    # injection: :meth:`AngularSourceSink.from_isotropic` (the Pattern-7
    # producer-normalised factory, which takes the carrier as an
    # argument).

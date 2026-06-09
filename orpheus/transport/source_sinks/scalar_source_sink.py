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
and operator-loss *sinks*; see :mod:`orpheus.transport.source_sinks`). The cross-class arithmetic
with :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`
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

with the dunder using ``isinstance`` dispatch (or
:func:`functools.singledispatchmethod` — the formal version) to
recognize the containment direction:

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

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.units import SCALAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import ScalarField

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh
    from orpheus.transport.source_sinks.angular_source_sink import AngularSourceSink


__all__ = ["ScalarSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class ScalarSourceSink(ScalarField):
    r"""Isotropic source field :math:`Q(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(ng, nx, ny)`` in the principled layout
        (Issue #196 PR-INDEX-7).
    space : FunctionSpace
        The function space. Must have ``shape == (mesh.ng, mesh.nx,
        mesh.ny)``. Use :meth:`from_mesh` to derive automatically.
    mesh : SNMesh
        The SN phase-space carrier.

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

    #: The :class:`FunctionSpace` name for this leaf (preserves the
    #: pre-B.1 space identity). All mesh/shape/algebra/factory machinery
    #: is inherited from :class:`ScalarField` / :class:`BulkField`.
    _SPACE_NAME: ClassVar[str] = "scalar_source_sink"

    #: Dimensional identity (View-G, B.4): scalar rate density
    #: ``1/(cm³·s)`` — :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS`,
    #: shared with ``ScalarResidual`` (same units, different role → the
    #: gate is class identity). Metadata, not the gate.
    UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS

    # ── Algebra extensions (over Field) ──────────────────────────────

    def __add__(self, other):
        r"""Add a partner source — dispatches by partner type.

        * :class:`ScalarSourceSink` partner → within-class add via
          Field's inherited algebra.
        * :class:`AngularSourceSink` partner → canonical subspace-
          containment injection: broadcast ``self.values`` across the
          Ω axis (the structural map ``iso → 1 ⊗ iso``) and add into
          the per-ordinate operand's space. Returns
          :class:`AngularSourceSink`. Commutative — see
          :meth:`AngularSourceSink.__add__`.

        See module docstring for the principled justification (refined
        Issue #207: cross-class dunders permitted via canonical
        subspace containment).
        """
        from orpheus.transport.source_sinks.angular_source_sink import (
            AngularSourceSink,
        )
        if isinstance(other, AngularSourceSink):
            if self.mesh is not other.mesh:
                raise ValueError(
                    "ScalarSourceSink + AngularSourceSink across "
                    "distinct SNMesh instances is forbidden — both "
                    "fields are mesh-bound."
                )
            # Canonical injection: broadcast iso → per-ordinate, add.
            return replace(
                other, values=self.values[None] + other.values,
            )
        # Within-class (and any other-type → NotImplemented) via Field.
        return super().__add__(other)

    def __radd__(self, other):
        r"""Reflected add — Python falls here when
        ``AngularSourceSink.__add__(ScalarSourceSink)`` returns
        NotImplemented. Delegates to :meth:`__add__` (commutative).
        """
        from orpheus.transport.source_sinks.angular_source_sink import (
            AngularSourceSink,
        )
        if isinstance(other, AngularSourceSink):
            return self.__add__(other)
        return NotImplemented

    # ── Conversion (named cross-class composition) ───────────────────

    def as_per_ordinate(self) -> "AngularSourceSink":
        r"""Broadcast this isotropic source to a per-ordinate source.

        Returns a :class:`AngularSourceSink` whose every ordinate slice
        equals ``self.values``. Uses ``np.broadcast_to(...).copy()`` so
        the resulting array is writable and owns its data.

        Use this when the producer-side normalisation (``/sum_w``,
        Pattern 7) has ALREADY been applied to ``self.values`` and the
        consumer just needs the per-ordinate broadcast. When the
        normalisation has NOT been applied yet, use
        :meth:`AngularSourceSink.from_isotropic` instead.
        """
        from orpheus.transport.source_sinks.angular_source_sink import (
            AngularSourceSink,
        )
        N = self.mesh.quad.N
        target_shape = (N, *self.values.shape)
        per_ord_values = np.broadcast_to(
            self.values[None], target_shape,
        ).copy()
        return AngularSourceSink.from_mesh(per_ord_values, self.mesh)

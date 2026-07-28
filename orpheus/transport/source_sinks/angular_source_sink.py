r"""Per-ordinate source field :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`.

The L2 typed wrapper for an anisotropic / per-ordinate source. Same
storage shape as :class:`~orpheus.transport.fields.angular_flux.AngularFlux` (``(N,
ng, nx, ny)``), but a structurally distinct physical quantity (source
density, not flux density). Carries the :math:`P_\ell \ge 1` Galerkin
reconstruction contribution AND any MMS-style per-ordinate external
source.

Migration status (Depth B step D-F → B.1 → B.2 → B.5)
======================================================

Moved (D-F) from ``orpheus.sn.sources.PerOrdinateSource`` to here;
re-parented (B.1) onto
:class:`~orpheus.transport.fields._bases.AngularField` (the storage-base
dedup — mesh / shape-check / ``from_mesh`` / ``_check_partner`` now live
there); renamed ``PerOrdinateSource`` → ``AngularSource`` (B.2, hard
rename, no shim) to complete the ``{Angular, Scalar} × {Flux, Source,
Residual}`` role grid; then ``AngularSource`` → ``AngularSourceSink``
(B.5, source/sink role rename — the leaf holds both production *sources*
and operator-loss *sinks*; see :mod:`orpheus.transport.source_sinks`). The cross-class arithmetic
with :class:`~orpheus.transport.source_sinks.scalar_source_sink.ScalarSourceSink`
is **PRESERVED** — the refined Issue #207 principle (recorded in the
conversation 2026-05-27) permits cross-class dunders when there is a
canonical subspace-containment relation. :class:`ScalarSourceSink`
lives in the subspace of :class:`AngularSourceSink` where every
ordinate carries the same value; the canonical injection
``iso → 1 ⊗ iso`` (broadcast across the Ω axis) is applied inside
the dunder. The result type is the LARGER (containing) space:

.. code-block:: python

    Q_total = aniso + iso          # → AngularSourceSink (commutative)
    Q_total = iso + aniso          # → AngularSourceSink (same result)

Two named factories complement the dunder for Pattern 7
producer-side normalisation (``/sum_w`` baked in vs not):

* :meth:`from_isotropic(values, mesh)` — broadcast + ``/sum_w`` (the
  Pattern 7 entry point).
* :meth:`ScalarSourceSink.as_per_ordinate` — broadcast WITHOUT
  ``/sum_w`` (when caller has already done it).

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

Per-ordinate reaction-rate density:
:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`; eV-free, binned-
energy convention), differing from :class:`ScalarSourceSink`'s units by a
``1/sr`` factor. The
cross-class dunder applies the structural injection (broadcast)
without ``/(4π sr)`` normalisation; the Pattern 7 ``/sum_w``
discipline is the caller's concern (apply before the dunder) or
goes through :meth:`from_isotropic` (named factory bakes it in).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import AngularField

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh


__all__ = ["AngularSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularSourceSink(AngularField):
    r"""Per-ordinate source field :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    space : FunctionSpace
        The function space. Must have ``shape == (mesh.quad.N,
        mesh.ng, *mesh.spatial_shape)``. Use :meth:`from_mesh` to derive
        automatically.
    mesh : SNMesh
        The SN phase-space carrier.

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`.
    Same-class arithmetic is closed. Cross-class arithmetic — including
    with :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
    (same shape, different physical kind) — is rejected by Field's
    Layer 1 class-identity gate.
    """

    #: The :class:`FunctionSpace` name for this leaf (preserves the
    #: pre-B.1 space identity). All mesh/shape/algebra/factory machinery
    #: is inherited from :class:`AngularField` / :class:`BulkField`.
    _SPACE_NAME: ClassVar[str] = "angular_source_sink"

    #: Dimensional identity (View-G, B.4): per-ordinate rate density
    #: ``1/(cm³·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`,
    #: shared with ``AngularResidual`` (same units, different role → the
    #: gate is class identity). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS

    # ── Algebra extensions (over Field) ──────────────────────────────

    def __add__(self, other):
        r"""Add a partner source — dispatches by partner type.

        * :class:`AngularSourceSink` partner → within-class add via
          Field's inherited algebra.
        * :class:`ScalarSourceSink` partner → canonical subspace-
          containment injection: broadcast the iso operand across the
          Ω axis (the structural map ``iso → 1 ⊗ iso``) and add into
          this per-ordinate operand's space. Returns
          :class:`AngularSourceSink`. Commutative — delegates to
          :meth:`ScalarSourceSink.__add__` for the broadcast-and-add
          (single source of truth for the iso injection logic).

        See module docstring for the principled justification (refined
        Issue #207: cross-class dunders permitted via canonical
        subspace containment). DELIBERATELY untyped — see
        :meth:`ScalarSourceSink.__add__` (the containment exception is
        statically unspellable against Field's ``(T, T) -> T``; #288).
        """
        from orpheus.transport.source_sinks.scalar_source_sink import ScalarSourceSink
        if isinstance(other, ScalarSourceSink):
            # Delegate to ScalarSourceSink.__add__ — Pattern 2 (single
            # source of truth for the broadcast-and-add).
            return other + self
        # Within-class (and any other-type → NotImplemented) via Field.
        return super().__add__(other)

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_isotropic(
        cls, iso_values: NDArray, mesh: "SNMesh",
    ) -> "AngularSourceSink":
        r"""Project an iso scalar source :math:`Q(\vec r, g)` to per-ordinate.

        Producer-side Pattern 7 normalisation: applies the ``/sum_w``
        projection and broadcasts across the :math:`N` ordinates.
        Returns a per-ordinate source whose every ordinate slice equals
        ``iso_values / sum_w``.

        This is the canonical entry for **external, user-supplied** iso
        scalar sources (fixed-source problems) that must enter the
        per-ordinate SN transport equation
        :math:`(\Omega\cdot\nabla + \sigma_t)\psi_n = Q/W + q_{\rm aniso,n}`.
        It is also the canonical migration path from the retired
        cross-class :meth:`ScalarSourceSink.__add__(AngularSourceSink)`
        dunder (see module docstring).

        Parameters
        ----------
        iso_values : NDArray
            Iso scalar source, shape ``(ng, nx, ny)``.
        mesh : SNMesh
            Phase-space carrier.

        Returns
        -------
        AngularSourceSink
            Per-ordinate source with ``Q/sum_w`` broadcast across N.

        See also
        --------
        :meth:`ScalarSourceSink.as_per_ordinate` — broadcast WITHOUT
            the ``/sum_w`` normalisation (use when caller has already
            divided by sum_w).
        """
        expected = (mesh.ng, *mesh.spatial_shape)
        if iso_values.shape != expected:
            raise ValueError(
                f"AngularSourceSink.from_isotropic expects iso shape "
                f"(ng, *spatial) = {expected}; got {iso_values.shape}"
            )
        sum_w = float(mesh.quad.weights.sum())
        N = mesh.quad.N
        per_ord_values = np.broadcast_to(
            (iso_values / sum_w)[None],
            (N, *expected),
        ).copy()
        return cls.from_mesh(per_ord_values, mesh)

    # ── Selectors ────────────────────────────────────────────────────

    def at_ordinate(self, n: int) -> NDArray:
        r"""Return the per-ordinate slice ``values[n]``, shape ``(ng, nx, ny)``."""
        return self.values[n]

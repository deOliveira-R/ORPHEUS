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
dedup — the shape-check / space derivation / ``_check_partner`` now
live there); renamed ``PerOrdinateSource`` → ``AngularSource`` (B.2, hard
rename, no shim) to complete the ``{Angular, Scalar} × {Flux, Source,
Residual}`` role grid; then ``AngularSource`` → ``AngularSourceSink``
(B.5, source/sink role rename — the leaf holds both production *sources*
and operator-loss *sinks*; see :mod:`orpheus.transport.fields`). The algebra is
:class:`~orpheus.numerics.field.Field`'s, same-class and closed: no leaf
overrides ``__add__``. A scalar source enters an angular one through
:meth:`from_isotropic` (the angular section: divide by the mass, then
scatter across the ordinates) and then a within-class add;
``ScalarSourceSink + AngularSourceSink`` is refused by
``Field._check_partner`` (the units differ by ``1/sr``).

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

Per-ordinate reaction-rate density:
:math:`[1/(\mathrm{cm^3 \cdot s \cdot sr})]`
(:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`; eV-free, binned-
energy convention), differing from :class:`~orpheus.transport.fields.scalar_source_sink.ScalarSourceSink`'s units by a
``1/sr`` factor. The Pattern 7 ``/sum_w``
discipline (the divisor is the angular section's mass, the correctly rounded
sum of the stored weights) is applied by :meth:`from_isotropic`.

The cross-class dunder of this class and ``BulkField.CONTAINS`` were retired
on 2026-11-20 (#708, #288).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import AngularField
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.numerics.axis_marginal import AxisSectionOperator

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.fields.scalar_source_sink import ScalarSourceSink
    from orpheus.transport.mesh.material_mesh import MaterialMesh


__all__ = ["AngularSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class AngularSourceSink(AngularField, flux=AngularFlux):
    r"""Per-ordinate source field :math:`Q^{\rm aniso}(\vec r, \hat\Omega_n, g)`.

    Parameters
    ----------
    values : NDArray
        Field values of shape ``(N, ng, nx, ny)`` in the principled
        layout (Issue #196 PR-INDEX-5/7).
    space : FunctionSpace
        The function space — the carrier's cached mint
        (``mesh.angular_bulk_space``, or ``mesh.angular_trial_space``
        when the scheme's within-cell moment tail is carried; CS4b S5).

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`.
    Same-class arithmetic is closed. Cross-class arithmetic — including
    with :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
    (same shape, different physical kind) — is rejected by Field's
    Layer 1 class-identity gate.
    """

    #: Dimensional identity (View-G, B.4): per-ordinate rate density
    #: ``1/(cm³·s·sr)`` — :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`,
    #: shared with ``AngularResidual`` (same units, different role → the
    #: gate is class identity). Metadata, not the gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS

    # ── Construction factories ───────────────────────────────────────

    @classmethod
    def from_isotropic(
        cls, iso_values: NDArray, problem: "MaterialMesh",
    ) -> "AngularSourceSink":
        r"""Project an iso scalar source :math:`Q(\vec r, g)` to per-ordinate.

        Producer-side Pattern 7 normalisation: applies the ``/sum_w``
        projection and broadcasts across the :math:`N` ordinates —
        realized (CS4b S6.2) by the angular axis's section
        :meth:`AxisSectionOperator.of
        <orpheus.numerics.axis_marginal.AxisSectionOperator.of>`, whose
        divisor is the correctly rounded mass. Returns a
        per-ordinate source whose every ordinate slice equals
        ``iso_values / sum_w``.

        This is the canonical entry for **external, user-supplied** iso
        scalar sources (fixed-source problems) that must enter the
        per-ordinate SN transport equation
        :math:`(\Omega\cdot\nabla + \sigma_t)\psi_n = Q/W + q_{\rm aniso,n}`.

        Parameters
        ----------
        iso_values : NDArray
            Iso scalar source, shape ``(ng, nx, ny)``.
        problem : MaterialMesh
            The angular phase-space carrier (an SNProblem); must own an
            angular bulk space.

        Returns
        -------
        AngularSourceSink
            Per-ordinate source with ``Q/sum_w`` broadcast across N.

        Raises
        ------
        ValueError
            If ``iso_values`` is not shaped ``(ng, *spatial)``.
        TypeError
            If the carrier owns no angular bulk space (not an SN carrier).
        """
        expected = (problem.ng, *problem.spatial_shape)
        if iso_values.shape != expected:
            raise ValueError(
                f"AngularSourceSink.from_isotropic expects iso shape "
                f"(ng, *spatial) = {expected}; got {iso_values.shape}"
            )
        # The section's divisor is the correctly rounded mass (math.fsum of the stored weights), the one definition every reader of 1/W divides by (#708).
        return cls.from_scalar(iso_values, cls._space_for_problem(problem))

    @classmethod
    def from_scalar(cls, iso: "NDArray | ScalarSourceSink", into: "FunctionSpace") -> "AngularSourceSink":
        r"""The one injection of an iso scalar source into the per-ordinate
        class: the angular section of ``into`` (divide by the mass, then
        scatter across the ordinates), wrapped in this class on ``into``.
        :meth:`from_isotropic` and the lift's producer-side combine both
        call it, so the injection is spelled once."""
        return cls(values=AxisSectionOperator.of(into, "angular").apply(iso), space=into)

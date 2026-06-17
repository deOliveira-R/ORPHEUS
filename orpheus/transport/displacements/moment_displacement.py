r"""Harmonic-moment flux displacement :math:`\Delta\phi_\ell^m(\vec r, g)`.

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.harmonic_moment_field.HarmonicMomentField` —
the increment between two moment-field iterates. This is the SI iterate
displacement on the **windowed 2-D Cartesian** path, where the held iterate's
``bulk`` is a :class:`HarmonicMomentField` (Phase 5a): ``moments^{(i)} ⊖
moments^{(i-1)}``. Same storage shape ``(L+1, 2L+1, ng, nx, ny)``, units, and
:class:`~orpheus.numerics.space.TensorProductSpace` as the flux; a DISTINCT
class. See :class:`~orpheus.transport.displacements._displacement.Displacement`.

Carries ``L`` (like its flux sibling — the moment shape is leaf-specific) and
mirrors :meth:`HarmonicMomentField._phase_space_shape`. Born ONLY from
``HarmonicMomentField ⊖ HarmonicMomentField`` (the :class:`FluxRole` field-copy
mint shares the flux's tensor-product space).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.spaces.spatial_moment_space import spatial_moment_tail
from orpheus.numerics.units import SCALAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import MomentField


__all__ = ["MomentDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class MomentDisplacement(Displacement, MomentField):
    r"""L2 harmonic-moment flux displacement on ``(L+1, 2L+1, ng, *spatial)``."""

    L: int

    #: Optional within-cell spatial-moment basis size per axis (#240
    #: D5b-S3-A0) — mirrors :attr:`HarmonicMomentField.spatial_moments` so
    #: the displacement minted by ``moments ⊖ moments`` (the
    #: :class:`~orpheus.transport.fields._flux_role.FluxRole` field-copy
    #: mint, which copies EVERY init field) lives in the flux's own
    #: (possibly widened) tensor-product space. Default ``1`` →
    #: byte-identical.
    spatial_moments: int = 1

    #: Same units as :class:`HarmonicMomentField` (a moment is angle-integrated,
    #: :math:`1/(\mathrm{cm^2 \cdot s})`).
    UNITS: ClassVar[Unit] = SCALAR_FLUX_UNITS

    def _phase_space_shape(self) -> tuple[int, ...]:
        r"""The ``(L+1, 2L+1, ng, *spatial[, spatial_moments**ndim])`` shape.

        Mirrors :meth:`HarmonicMomentField._phase_space_shape` (incl. the
        optional spatial-moment tail #240 D5b-S3-A0) so the shared
        :meth:`BulkField.__post_init__` validator accepts the displacement on
        the flux's tensor-product space.  Spatial tail rank-``d`` via
        ``mesh.spatial_shape`` (no phantom ``ny``); the "append iff > 1"
        gate is single-sourced from
        :func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`.
        """
        n_moments = self.spatial_moments ** self.mesh.ndim
        return (
            self.L + 1, 2 * self.L + 1,
            self.mesh.ng, *self.mesh.spatial_shape,
            *spatial_moment_tail(n_moments),
        )

r"""The homogeneous body a one-material reference generator solves on.

The continuous reference generators of the singular-eigenfunction, F_N,
Galerkin-spectral and trajectory-resolvent families (``Spectrum``,
``MomentSpace``, ``BasisSpace``, ``Billiard``) solve on ONE material
filling a slab of width :math:`L`, or a cylinder or sphere of radius
:math:`R` centred at the origin. Each reads that body off a
:class:`~orpheus.geometry.StructuredGeometry`, and a geometry can say
more than such a body can: intervals of different materials, or
a hollow cylinder or sphere (:math:`r_0 > 0`), whose width
:math:`r_R - r_0` is a shell thickness and not a radius. Read as a
body, the first would be solved as its first material alone and the
second as a solid of the wrong radius, both silently.

:func:`homogeneous_body` is the one place that reading happens: it
returns the body when the geometry is one, and refuses otherwise,
naming the generator that asked. A slab has no centre, so its
:math:`r_0` is a translation and every one-material slab is a body.
"""
from __future__ import annotations

from dataclasses import dataclass

from orpheus.geometry import CoordSystem, StructuredGeometry


@dataclass(frozen=True)
class HomogeneousBody:
    r"""One material filling a slab of width ``extent_cm`` or a solid
    cylinder or sphere of radius ``extent_cm``.

    Attributes
    ----------
    coord : CoordSystem
        The body's shape.
    extent_cm : float
        The slab width, or the radius of the cylinder or sphere (cm).
    mat_id : int
        The id of the one material that fills it.
    """

    coord: CoordSystem
    extent_cm: float
    mat_id: int


def homogeneous_body(
    geometry: StructuredGeometry, *, owner: str,
) -> HomogeneousBody:
    """The homogeneous body ``geometry`` describes, or a keyed refusal.

    Parameters
    ----------
    geometry : StructuredGeometry
        The geometry a reference generator was handed.
    owner : str
        The generator's name, for the refusal message.

    Raises
    ------
    ValueError
        If the geometry holds more than one material, or is a hollow
        cylinder or sphere. Several intervals of one material are one
        body (its interior breakpoints are not material interfaces).
    """
    materials = set(geometry.mat_ids)
    if len(materials) != 1:
        raise ValueError(
            f"{owner} solves one material filling the whole body; the "
            f"geometry holds the materials {sorted(materials)} "
            f"(breakpoints {geometry.breakpoints})."
        )
    if geometry.is_hollow:
        raise ValueError(
            f"{owner} solves a solid {geometry.coord.name.lower()} body "
            f"centred at r = 0; the geometry is hollow (r_0 = "
            f"{geometry.breakpoints[0]!r})."
        )
    return HomogeneousBody(
        coord=geometry.coord,
        extent_cm=geometry.domain_extent_cm,
        mat_id=geometry.mat_ids[0],
    )


__all__ = ["HomogeneousBody", "homogeneous_body"]

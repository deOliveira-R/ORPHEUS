r"""The body a continuous reference generator solves on, read once from a geometry.

The continuous reference generators of the singular-eigenfunction, F_N,
Galerkin-spectral and trajectory-resolvent families (``Spectrum``,
``MomentSpace``, ``BasisSpace``, ``Billiard``) solve a handful of body
shapes, the ones the reference literature (Sood, Forster & Parsons 2003,
and the papers it collects) states its benchmarks on. A
:class:`~orpheus.geometry.StructuredGeometry` can say more than any one
generator solves, so each generator must read which shape it was handed and
refuse the shapes it does not serve.

:func:`reference_body` is the one place that reading happens. It is TOTAL:
every geometry is exactly one :data:`ReferenceBody`, and the reading knows
no solver. Adjacent intervals of one material are first merged into one
material RUN (an interior breakpoint between them is not an interface).
Then, with :math:`n` the number of runs:

========================  ==============================================
shape                     the geometry
========================  ==============================================
:class:`HomogeneousBody`  :math:`n = 1`, a slab or a solid cylinder or
                          sphere
:class:`HollowBody`       :math:`n = 1`, a hollow cylinder or sphere
:class:`ReflectedSlab`    a slab with :math:`n = 3` whose outer runs share
                          one material and one width
:class:`LayeredBody`      every other geometry of :math:`n \ge 2` runs
========================  ==============================================

Each generator then matches on the shape and refuses what it does not
serve through :func:`refuse_unserved`, the one door, which names the
solver the generator lacks. Which generator serves which shape is recorded
on the structured-geometry theory page.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, NoReturn

from orpheus.geometry import BC, CoordSystem, StructuredGeometry

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BoundaryTraceLaw


@dataclass(frozen=True)
class HomogeneousBody:
    r"""One material filling a slab of width ``extent_cm``, or a solid
    cylinder or sphere of radius ``extent_cm``."""

    coord: CoordSystem
    extent_cm: float
    mat_id: int

    @property
    def mat_ids(self) -> tuple[int, ...]:
        """The material of each run: one."""
        return (self.mat_id,)


@dataclass(frozen=True)
class HollowBody:
    r"""One material filling a hollow cylinder (an annulus) or a hollow
    sphere between ``inner_radius_cm`` and ``outer_radius_cm``."""

    coord: CoordSystem
    inner_radius_cm: float
    outer_radius_cm: float
    mat_id: int

    @property
    def mat_ids(self) -> tuple[int, ...]:
        """The material of each run: one."""
        return (self.mat_id,)


@dataclass(frozen=True)
class ReflectedSlab:
    r"""A core slab of width ``core_width_cm`` with the same reflector,
    ``reflector_width_cm`` thick, on both faces (Sood's symmetric
    two-media slabs)."""

    core_width_cm: float
    reflector_width_cm: float
    core_mat_id: int
    reflector_mat_id: int

    @property
    def coord(self) -> CoordSystem:
        """A reflected slab is Cartesian."""
        return CoordSystem.CARTESIAN

    @property
    def mat_ids(self) -> tuple[int, ...]:
        """The material of each run, left to right: reflector, core, reflector."""
        return (self.reflector_mat_id, self.core_mat_id, self.reflector_mat_id)


@dataclass(frozen=True)
class LayeredBody:
    r"""Two or more material runs in any other arrangement: a layered
    cylinder or sphere (solid or hollow), or a slab that is not a symmetric
    reflected slab.

    Attributes
    ----------
    coord : CoordSystem
        The body's shape.
    breakpoints : tuple of float
        The run boundaries, :math:`r_0 < \dots < r_n`, a subset of the
        geometry's breakpoints.
    mat_ids : tuple of int
        One material per run; adjacent runs differ.
    is_hollow : bool
        The geometry's :attr:`~orpheus.geometry.StructuredGeometry.is_hollow`,
        read at classification (the one place the centre's role is decided).
    """

    coord: CoordSystem
    breakpoints: tuple[float, ...]
    mat_ids: tuple[int, ...]
    is_hollow: bool


ReferenceBody = HomogeneousBody | HollowBody | ReflectedSlab | LayeredBody


def _material_runs(geometry: StructuredGeometry) -> tuple[tuple[float, ...], tuple[int, ...]]:
    """The run boundaries and the material of each run (adjacent equal
    materials merged)."""
    breakpoints = [geometry.breakpoints[0]]
    mat_ids: list[int] = []
    for (_, outer), mat_id in zip(itertools.pairwise(geometry.breakpoints), geometry.mat_ids):
        if mat_ids and mat_ids[-1] == mat_id:
            breakpoints[-1] = outer
        else:
            mat_ids.append(mat_id)
            breakpoints.append(outer)
    return tuple(breakpoints), tuple(mat_ids)


def _equal_widths(breakpoints: tuple[float, ...]) -> bool:
    r"""The two outer runs of a three-run slab have one width, up to the
    rounding of the stored breakpoints.

    Each breakpoint is a stored float, and a registry that states its layers
    as thicknesses builds them by a sequential sum (``from_thicknesses``),
    so the right reflector's width :math:`r_3 - r_2` carries the rounding
    of :math:`r_2` and :math:`r_3`, at most one unit in the last place of
    the largest magnitude each; the difference of two widths carries at most
    four. Anything larger is a real asymmetry.
    """
    r_0, r_1, r_2, r_3 = breakpoints
    scale = max(abs(r_0), abs(r_3))
    return math.isclose(r_1 - r_0, r_3 - r_2, rel_tol=0.0, abs_tol=4.0 * math.ulp(scale))


def reference_body(geometry: StructuredGeometry) -> ReferenceBody:
    """The body shape ``geometry`` describes (total: every geometry is one)."""
    breakpoints, mat_ids = _material_runs(geometry)
    coord = geometry.coord
    if len(mat_ids) == 1:
        if geometry.is_hollow:
            return HollowBody(
                coord=coord,
                inner_radius_cm=breakpoints[0],
                outer_radius_cm=breakpoints[1],
                mat_id=mat_ids[0],
            )
        return HomogeneousBody(
            coord=coord, extent_cm=breakpoints[1] - breakpoints[0], mat_id=mat_ids[0],
        )
    if (
        coord is CoordSystem.CARTESIAN
        and len(mat_ids) == 3
        and mat_ids[0] == mat_ids[2]
        and _equal_widths(breakpoints)
    ):
        return ReflectedSlab(
            core_width_cm=breakpoints[2] - breakpoints[1],
            reflector_width_cm=breakpoints[1] - breakpoints[0],
            core_mat_id=mat_ids[1],
            reflector_mat_id=mat_ids[0],
        )
    return LayeredBody(
        coord=coord, breakpoints=breakpoints, mat_ids=mat_ids,
        is_hollow=geometry.is_hollow,
    )


def refuse_unserved(what: str, *, owner: str, missing: str) -> NoReturn:
    """Refuse a configuration ``owner`` does not serve, naming what it lacks.

    **SCOPE-BOUNDARY[guard]** — machinery: the reference solver ``missing`` names, per owner (#536).
    ruling: the user, 2026-09-29, P1 step 2b of ``.claude/plans/reference_cache.md``.
    revisit: a solver for the refused configuration lands in the owner's family (#536 lists Sood's cases).

    The one door every generator refuses through, so the refusals share
    one spelling and one ledger entry. ``what`` describes the refused
    configuration (a body shape, :func:`describe`, or a boundary law);
    ``missing`` names the solver the owner would need for it.
    """
    raise NotImplementedError(f"{owner} does not solve a {what}: {missing} (#536).")


def describe(body: ReferenceBody) -> str:
    """The body's shape in words, for a refusal."""
    kind = body.coord.name.lower()
    match body:
        case HomogeneousBody():
            return f"homogeneous {kind} body"
        case HollowBody():
            return f"hollow {kind} body"
        case ReflectedSlab():
            return "symmetric reflected slab"
        case LayeredBody():
            hollow = "hollow " if body.is_hollow else ""
            return f"{hollow}layered {kind} body of {len(body.mat_ids)} material runs {body.mat_ids}"


def specular_albedo(law: BC | BoundaryTraceLaw, *, owner: str) -> float:
    r"""The specular albedo :math:`\alpha \in [0, 1]` a boundary law declares.

    The continuous references parametrise a boundary by one specular albedo:
    :math:`\alpha = 0` is vacuum, :math:`\alpha = 1` perfect mirror
    reflection, and a value between is partial specular reflection. The laws
    that say this are:

    * vacuum (``BC.vacuum``, :class:`~orpheus.geometry.boundary.VacuumInflow`): 0;
    * mirror reflection (``BC.reflective``: 1;
      :class:`~orpheus.geometry.boundary.ReflectiveBoundary`: its albedo);
    * partial specular reflection (``BC("partial", {"albedo": a})``: ``a``;
      :class:`~orpheus.geometry.boundary.AlbedoBoundary` whose re-emission is
      :class:`~orpheus.geometry.boundary.SpecularReturn`: its albedo).

    Any other law returns neutrons in another angular shape (white,
    isotropic, an unstated re-emission, a prescribed inflow, periodic), which
    a specular albedo cannot express, and is refused through the door.
    """
    from orpheus.geometry.boundary import (
        AlbedoBoundary,
        ReflectiveBoundary,
        SpecularReturn,
        VacuumInflow,
    )

    match law:
        case BC(kind="vacuum") | VacuumInflow():
            return 0.0
        case BC(kind="reflective"):
            return 1.0
        case BC(kind="partial"):
            return law.to_alpha()
        case ReflectiveBoundary():
            return float(law.albedo)
        case AlbedoBoundary(reemission=SpecularReturn()):
            return float(law.albedo)
    refuse_unserved(
        f"boundary law {law!r}", owner=owner,
        missing="the continuous references parametrise a boundary by a specular "
        "albedo, and this law returns neutrons in another angular shape",
    )


def specular_albedos(geometry: StructuredGeometry, *, owner: str) -> tuple[float, ...]:
    """The specular albedo at each boundary point of ``geometry``, inner first."""
    return tuple(specular_albedo(law, owner=owner) for law in geometry.boundaries)


def require_vacuum(geometry: StructuredGeometry, *, owner: str, missing: str) -> None:
    """Refuse a geometry any of whose boundary points is not vacuum."""
    albedos = specular_albedos(geometry, owner=owner)
    if any(albedo != 0.0 for albedo in albedos):
        refuse_unserved(
            f"body with reflecting boundaries (specular albedos {albedos})",
            owner=owner, missing=missing,
        )


__all__ = [
    "HollowBody",
    "HomogeneousBody",
    "LayeredBody",
    "ReferenceBody",
    "ReflectedSlab",
    "describe",
    "reference_body",
    "refuse_unserved",
    "require_vacuum",
    "specular_albedo",
    "specular_albedos",
]

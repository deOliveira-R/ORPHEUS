r"""Method-agnostic geometry recipe for cross-method case construction.

:class:`MeshTemplate` is a recipe — a small, immutable description of
``(geometry_kind, critical_dimension_{mfp,cm}, n_groups, mat_id,
bc_left, bc_right)`` — that BOTH production solvers AND continuous
reference solvers can consume at the case-description boundary:

* **Production solvers** (``solve_cp``, ``solve_sn``, etc.) consume
  the discretised :class:`~orpheus.geometry.mesh.Mesh1D` produced by
  :meth:`MeshTemplate.build`.
* **Continuous reference solvers** (``solve_fn_*``,
  ``solve_greens_function_*``, etc.) consume the descriptor scalars
  directly via :attr:`MeshTemplate.domain_extent_cm`,
  :attr:`MeshTemplate.critical_dimension_cm`,
  :attr:`MeshTemplate.bc_left`, :attr:`MeshTemplate.bc_right` —
  with no need to discretise into cells.

The class lived inside
``orpheus.derivations.continuous.sood_registry.la13511`` until the
2026-05-03 R0.5 promotion (geometry-handling unification audit at
``.claude/scratch/geometry_handling_unification_audit.md`` §"Q3.b
Promote-and-unify"). The audit identified that:

1. The descriptor is already method-agnostic — its attributes
   describe the SHAPE of the domain, not any particular solver.
2. Both production-side (sood_registry → ``Mesh1D`` → ``solve_cp``)
   and reference-side (cross-method adapter → ``solve_fn_slab``)
   call sites benefit from a single canonical recipe.
3. The promotion is independent of the trajectory_resolvent
   R1/R2/R3 chord-oracle / power-iterate refactors.

The legacy import path
``from orpheus.derivations.continuous.sood_registry.la13511 import
MeshTemplate`` is preserved as a re-export for backward
compatibility — existing callers do not need to change.

References
----------

* :doc:`/skills/algebra-of-record` §"Branch 1 / Branch 2 separation"
  — :class:`MeshTemplate` lives at the math layer, NOT inside any
  specific solver's algebra-of-record.
* ``.claude/scratch/geometry_handling_unification_audit.md`` —
  the audit memo that identified the promotion opportunity.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.factories import Zone, mesh1d_from_zones
from orpheus.geometry.mesh import BC, Mesh1D


_GEOMETRY_TO_COORD: dict[str, CoordSystem] = {
    "slab": CoordSystem.CARTESIAN,
    "sphere": CoordSystem.SPHERICAL,
    "cylinder": CoordSystem.CYLINDRICAL,
}


@dataclass(frozen=True)
class MeshTemplate:
    """Geometry + critical-dimension recipe for a benchmark case.

    Method-agnostic: stores the shape of the domain (slab /
    sphere / cylinder / infinite / ISLC), the critical dimension as
    published, the group count, and the boundary conditions. The
    :meth:`build` method constructs a concrete :class:`Mesh1D` at the
    requested refinement using
    :func:`orpheus.geometry.factories.mesh1d_from_zones`.

    For ``geometry == "infinite"`` no mesh exists: :meth:`build`
    raises :class:`ValueError`. Infinite-medium consumers should use
    :func:`...sood_registry.builders.build_materials` and ignore the
    mesh entirely.

    Parameters
    ----------
    geometry : str
        One of ``"infinite"``, ``"slab"``, ``"sphere"``, ``"cylinder"``,
        ``"ISLC"``.
    critical_dimension_mfp : float | None
        Critical radius (sphere/cylinder) or half-thickness (slab) in
        mean free paths. ``None`` for infinite-medium cases.
    critical_dimension_cm : float | None
        Same as above but in cm.
    n_groups : int
        Number of energy groups (1, 2, 3, or 6).
    mat_id : int
        Material identifier used in the constructed mesh. The matching
        :class:`~orpheus.data.macro_xs.mixture.Mixture` lives in
        ``case.materials[mat_id]``.
    bc_left : BC
        Inner / left boundary condition. For sphere/cylinder default
        is :attr:`BC.reflective` (centreline); for slab default is
        :attr:`BC.vacuum`.
    bc_right : BC
        Outer / right boundary condition. Default
        :attr:`BC.vacuum`.
    """

    geometry: str
    critical_dimension_mfp: float | None
    critical_dimension_cm: float | None
    n_groups: int
    mat_id: int = 0
    bc_left: BC = field(default_factory=lambda: BC.vacuum)
    bc_right: BC = field(default_factory=lambda: BC.vacuum)

    def __post_init__(self) -> None:
        if self.geometry not in {"infinite", "slab", "sphere", "cylinder", "ISLC"}:
            raise ValueError(f"Unknown geometry {self.geometry!r}")
        if self.geometry == "infinite":
            if self.critical_dimension_cm is not None:
                raise ValueError(
                    "infinite geometry must have critical_dimension_cm=None"
                )
        else:
            if self.critical_dimension_cm is None:
                raise ValueError(
                    f"geometry={self.geometry!r} requires critical_dimension_cm"
                )

    @property
    def domain_extent_cm(self) -> float:
        """Total mesh extent in cm — what :meth:`build` actually constructs.

        Conventions:

        * **Slab**: ``2 * critical_dimension_cm`` (full symmetric slab
          ``[0, 2a]`` with vacuum BCs at both ends — F_N convention).
        * **Sphere / cylinder**: ``critical_dimension_cm`` (radius;
          mesh is ``[0, R]`` with reflective BC at the centre and
          vacuum at the outer surface unless overridden).
        * **Infinite / ISLC**: undefined, raises.
        """
        if self.geometry == "infinite":
            raise ValueError("infinite geometry has no domain_extent")
        if self.geometry == "ISLC":
            raise NotImplementedError("ISLC domain_extent not implemented")
        assert self.critical_dimension_cm is not None  # narrowed above
        if self.geometry == "slab":
            return 2.0 * float(self.critical_dimension_cm)
        return float(self.critical_dimension_cm)

    def build(self, n_cells: int = 64) -> Mesh1D:
        r"""Construct a :class:`Mesh1D` at the published critical dimension.

        Conventions (see :attr:`domain_extent_cm`):

        * **Slab** (``geometry == "slab"``): builds the FULL symmetric
          slab ``[0, 2a]`` where ``a = critical_dimension_cm``. This is
          the F_N method's natural domain — :math:`a` is the
          half-thickness, the published critical configuration is the
          full slab. Default BCs are vacuum at both ends.
        * **Sphere / cylinder**: builds ``[0, R]`` where
          ``R = critical_dimension_cm``. Default BCs: reflective at
          ``r = 0`` (centreline / axis), vacuum at the outer surface.

        Parameters
        ----------
        n_cells : int
            Number of equal-volume sub-cells.

        Returns
        -------
        Mesh1D
            A 1-D mesh with ``n_cells`` cells, homogeneous material
            ``mat_id``, and BCs from this template.

        Raises
        ------
        ValueError
            If ``geometry == "infinite"``.
        NotImplementedError
            If ``geometry == "ISLC"``.
        """
        if self.geometry == "infinite":
            raise ValueError(
                "infinite-medium cases have no mesh; consume materials directly"
            )
        if self.geometry == "ISLC":
            raise NotImplementedError(
                "ISLC mesh construction is not implemented in the first slice"
            )
        coord = _GEOMETRY_TO_COORD[self.geometry]
        zones = [Zone(
            outer_edge=self.domain_extent_cm,
            mat_id=self.mat_id,
            n_cells=n_cells,
        )]
        mesh = mesh1d_from_zones(zones, coord=coord)
        # Re-stamp with BC fields (mesh1d_from_zones doesn't set BCs).
        return Mesh1D(
            edges=mesh.edges,
            mat_ids=mesh.mat_ids,
            coord=mesh.coord,
            precomputed_volumes=mesh.precomputed_volumes,
            bc_left=self.bc_left,
            bc_right=self.bc_right,
        )


__all__ = ["MeshTemplate"]

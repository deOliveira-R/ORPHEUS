r"""Discrete-mesh adapters wrapping production CP/SN solvers.

This module is the **discrete-adapter scaffold** required by the
test-architect spec § "Discrete-solver adapter scaffold (#4a)". It
exposes :class:`CPMeshAdapter` and :class:`SNMeshAdapter` — thin
facades that wrap the production
:func:`~orpheus.cp.solver.solve_cp` and
:func:`~orpheus.sn.solver.solve_sn` calls behind the
:class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol surface.

What this scaffold ships
------------------------

* :class:`CPMeshAdapter` and :class:`SNMeshAdapter` accept the
  production-protocol input pair (``materials: dict[int, Mixture]``,
  ``geometry_spec: GeometrySpec``) plus a
  :class:`~orpheus.derivations.common.discretization_spec.DiscretizationSpec`
  for the cell count + angular quadrature order + chord-quadrature
  order.
* They construct on demand and call the production solver behind a
  :meth:`solve_critical` method that returns a
  :class:`~orpheus.derivations.common.solution_types.CriticalSolution`.
* They satisfy the Protocol via the same shape contract as
  ``Billiard`` / ``MomentSpace`` — ``materials``, ``geometry_spec``,
  ``method_name``, ``solve_critical``.

What this scaffold does NOT ship
--------------------------------

* **No L4 cross-check against the continuous reference.** That gate
  is deferred to a follow-up dispatch (see the test placeholder
  ``test_cross_method_discrete_vs_reference_TODO`` in
  :mod:`tests.cross_method.test_discrete_adapters_smoke`). The
  scaffold's promise is "construct without crash + solve_critical
  returns a non-zero eigenvalue" — NOT correctness.
* **No fixed-source path.**
* **No multi-region geometry.** The scaffold consumes the spec's
  single ``mat_id`` field; multi-region would need a richer spec.

Structural-independence guarantee
---------------------------------

Per :doc:`/skills/algebra-of-record` § "Structural independence
applies above the trusted-library line", the discrete adapters
share NO in-house primitive with the continuous-reference adapters
above the trusted-library line:

* Continuous (``Billiard``, ``MomentSpace``) — Variant α resolvent
  / F_N Galerkin projection — both consume scipy/numpy primitives
  but no shared ORPHEUS module.
* Discrete (CP, SN) — production solvers in :mod:`orpheus.cp` /
  :mod:`orpheus.sn` — consume their own scipy/numpy primitives via
  separate code paths.

References
----------

* :doc:`/theory/transport_solver_protocol` § "Discrete adapter
  scaffold".
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.cp.solver import CPParams, solve_cp
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.discretization_spec import DiscretizationSpec
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.factories import Zone, mesh1d_from_zones
from orpheus.geometry.mesh import BC, Mesh1D
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn


# Geometry → CoordSystem mapping.
_COORD_FOR_GEOMETRY = {
    "slab": CoordSystem.CARTESIAN,
    "sphere": CoordSystem.SPHERICAL,
    "cylinder": CoordSystem.CYLINDRICAL,
}


@dataclass(frozen=True)
class CPMeshAdapter:
    r"""TransportSolver-conforming wrapper for the production CP solver.

    Wraps :func:`~orpheus.cp.solver.solve_cp`. The
    :meth:`solve_critical` method builds a :class:`Mesh1D` at the
    requested cell count, calls ``solve_cp`` with the requested
    chord-integration quadrature, and re-packs the
    :class:`~orpheus.cp.solver.CPResult` into a
    :class:`CriticalSolution` matching the Protocol shape.

    The discretization spec carries:

    * ``n_cells`` — production mesh refinement (passed to
      :meth:`GeometrySpec.build`).
    * ``n_chord_quad`` — chord integration order
      (``CPParams.n_quad_y``). Slab CP is closed-form in :math:`E_3`
      and ignores this field.
    * ``n_angular`` — ignored by CP (CP integrates analytically over
      angle).
    """

    materials: dict[int, Mixture]
    geometry_spec: GeometrySpec
    discretization: DiscretizationSpec
    method_name: str = "cp"

    @classmethod
    def from_problem(
        cls,
        *,
        materials: dict[int, Mixture],
        geometry_spec: GeometrySpec,
        discretization: DiscretizationSpec | None = None,
    ) -> "CPMeshAdapter":
        r"""Construct a CP discrete adapter from production-protocol inputs."""
        if geometry_spec.geometry not in _COORD_FOR_GEOMETRY:
            raise ValueError(
                f"CPMeshAdapter cannot build on geometry "
                f"{geometry_spec.geometry!r}; supported: "
                f"{sorted(_COORD_FOR_GEOMETRY)}"
            )
        return cls(
            materials=dict(materials),
            geometry_spec=geometry_spec,
            discretization=(discretization or DiscretizationSpec()),
        )

    def solve_critical(self) -> CriticalSolution:
        r"""Run :func:`solve_cp` on this configuration."""
        mesh = _build_mesh_with_default_bcs(
            self.geometry_spec, n_cells=self.discretization.n_cells,
            outer_bc=BC("white"),
        )
        params = CPParams(n_quad_y=self.discretization.n_chord_quad)
        result = solve_cp(self.materials, mesh, params)

        return CriticalSolution(
            eigenvalue=float(result.keff),
            eigenvalue_kind="k_eff",
            parameter_value=float(self.geometry_spec.domain_extent_cm),
            parameter_kind="domain_extent_cm",
            converged=True,
            metadata={
                "method": "cp",
                "n_cells": self.discretization.n_cells,
                "n_chord_quad": self.discretization.n_chord_quad,
                "raw_result": result,
            },
        )


@dataclass(frozen=True)
class SNMeshAdapter:
    r"""TransportSolver-conforming wrapper for the production SN solver.

    Wraps :func:`~orpheus.sn.solver.solve_sn`. Currently 1-D only
    (Gauss-Legendre quadrature); 2-D / Lebedev support is future
    work.
    """

    materials: dict[int, Mixture]
    geometry_spec: GeometrySpec
    discretization: DiscretizationSpec
    method_name: str = "sn"

    @classmethod
    def from_problem(
        cls,
        *,
        materials: dict[int, Mixture],
        geometry_spec: GeometrySpec,
        discretization: DiscretizationSpec | None = None,
    ) -> "SNMeshAdapter":
        r"""Construct an SN discrete adapter from production-protocol inputs."""
        if geometry_spec.geometry not in _COORD_FOR_GEOMETRY:
            raise ValueError(
                f"SNMeshAdapter cannot build on geometry "
                f"{geometry_spec.geometry!r}; supported: "
                f"{sorted(_COORD_FOR_GEOMETRY)}"
            )
        spec = discretization or DiscretizationSpec()
        if spec.n_angular % 2 != 0:
            raise ValueError(
                f"SNMeshAdapter requires even n_angular; got "
                f"{spec.n_angular}"
            )
        return cls(
            materials=dict(materials),
            geometry_spec=geometry_spec,
            discretization=spec,
        )

    def solve_critical(self) -> CriticalSolution:
        r"""Run :func:`solve_sn` on this configuration."""
        mesh = _build_mesh_with_default_bcs(
            self.geometry_spec, n_cells=self.discretization.n_cells,
            outer_bc=BC.reflective,
        )
        quadrature = GaussLegendre1D.create(
            n_ordinates=self.discretization.n_angular
        )
        result = solve_sn(self.materials, mesh, quadrature)

        return CriticalSolution(
            eigenvalue=float(result.keff),
            eigenvalue_kind="k_eff",
            parameter_value=float(self.geometry_spec.domain_extent_cm),
            parameter_kind="domain_extent_cm",
            converged=True,
            metadata={
                "method": "sn",
                "n_cells": self.discretization.n_cells,
                "n_angular": self.discretization.n_angular,
                "raw_result": result,
            },
        )


def _build_mesh_with_default_bcs(
    spec: GeometrySpec,
    *,
    n_cells: int,
    outer_bc: BC,
) -> Mesh1D:
    r"""Build a Mesh1D from a GeometrySpec with default discrete BCs.

    GeometrySpec defaults to vacuum-vacuum; production CP/SN want
    reflective centreline (sphere/cylinder) or reflective-reflective
    (closed-loop k_inf-style smoke tests).
    """
    coord = _COORD_FOR_GEOMETRY[spec.geometry]
    extent_cm = float(spec.domain_extent_cm)
    zones = [
        Zone(outer_edge=extent_cm, mat_id=spec.mat_id, n_cells=n_cells),
    ]
    mesh = mesh1d_from_zones(zones, coord=coord)

    is_user_default = (
        spec.bc_left == BC.vacuum and spec.bc_right == BC.vacuum
    )
    if is_user_default:
        if spec.geometry in ("sphere", "cylinder"):
            bc_left = BC.reflective
            bc_right = outer_bc
        else:
            bc_left = BC.vacuum
            bc_right = BC.vacuum
    else:
        bc_left = spec.bc_left
        bc_right = spec.bc_right

    return Mesh1D(
        edges=mesh.edges,
        mat_ids=mesh.mat_ids,
        coord=mesh.coord,
        precomputed_volumes=mesh.precomputed_volumes,
        bc_left=bc_left,
        bc_right=bc_right,
    )


__all__ = [
    "CPMeshAdapter",
    "SNMeshAdapter",
]

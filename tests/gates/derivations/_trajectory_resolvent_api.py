r"""One spelling of the trajectory-resolvent reference's API for the R7b2.2–R7b2.7 gates (#405 P2 step 7b.2.2).

Every production name the gates touch is resolved HERE, once (L100): if the
build lands a different argument name or module path, this file is the one
edit. The declared API (spec §1.7b.2, "The trajectory-resolvent references as
``ReferenceSolution``s"; the orchestrator's brief of 2026-10-03):

* ``orpheus.derivations.continuous.trajectory_resolvent.reference``:
  ``trajectory_resolvent_reference(specification, quadrature, *, max_iter,
  tol, initial_k) -> ReferenceSolution`` and ``TrajectoryResolventDerivation``;
* the multi-region chord oracles take ``at=`` (the evaluation radii,
  defaulting to the nodes) separately from the spline's knots;
* ``Symbolic.steps(r_range)``: the step locations of a weight inside
  ``r_range``, one definition with production's scope edge.
"""
from __future__ import annotations

import importlib
from typing import Any

import numpy as np

from orpheus.geometry import CoordSystem

MODULE = "orpheus.derivations.continuous.trajectory_resolvent.reference"
FACTORY = "trajectory_resolvent_reference"
DERIVATION = "TrajectoryResolventDerivation"

#: The gates' resolutions, keyed as ``Billiard``'s quadrature is: the sphere's rows take seconds, the cylinder's
#: 10–125 s (``[M]`` 2026-10-03, most of it the reading's flux integrals), so the cylinder's are ``slow``.
QUADRATURE = {
    CoordSystem.SPHERICAL: {"n_r": 8, "n_mu": 8, "n_traj_quad": 16},
    CoordSystem.CYLINDRICAL: {"n_r": 8, "n_mu_axial": 4, "n_phi_az": 8, "n_traj_quad": 16},
}
TOL = 1e-10
MAX_ITER = 500
INITIAL_K = 1.0


def module() -> Any:
    return importlib.import_module(MODULE)


def reference(specification: Any, coord: CoordSystem | None = None, **overrides: Any) -> Any:
    """The factory at the gates' resolution (``coord`` picks the quadrature; default the specification's)."""
    c = coord if coord is not None else specification.geometry.coord
    kwargs = {"max_iter": MAX_ITER, "tol": TOL, "initial_k": INITIAL_K} | overrides
    quadrature = kwargs.pop("quadrature", QUADRATURE[c])
    return getattr(module(), FACTORY)(specification, quadrature, **kwargs)


def derivation_class() -> type:
    return getattr(module(), DERIVATION)


def apply_at(oracle: Any, source_profile: np.ndarray, at: np.ndarray | None, n_traj_quad: int = 16) -> np.ndarray:
    """The oracle's transport of ``source_profile`` (on its knots) evaluated at ``at`` (``None``: the default)."""
    if at is None:
        return oracle.apply_operator(source_profile, 0.0, n_traj_quad=n_traj_quad)
    return oracle.apply_operator(source_profile, 0.0, n_traj_quad=n_traj_quad, at=at)


def steps(weight: Any, low: float, high: float) -> tuple[float, ...]:
    """The step locations of a ``Symbolic`` weight inside ``[low, high]``, sorted."""
    return tuple(sorted(float(x) for x in weight.steps((low, high))))


#: The two multi-region solvers, as the factory reaches them (the lazy-solve spy patches these attributes).
SOLVERS = (
    ("orpheus.derivations.continuous.trajectory_resolvent.greens_function", "solve_greens_function_sphere_mr"),
    ("orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder", "solve_greens_function_cylinder_mr"),
)


# ── the review round of 7b.2.2 (qa F1–F3, elegance B1) ─────────────────────

#: The fragment of ``Symbolic.steps``' scope-edge refusal (a step on a non-polynomial argument).
SCOPE_EDGE = "polynomial in r"
#: Any refusal of an unlocated non-smooth construct (the allow-list's message, whatever its final wording).
UNLOCATED_REFUSAL = r"(?i)polynomial in r|not located|allow|scope boundary"


def answer(ref: Any) -> Any:
    """The derivation's converged answer (solves on first access)."""
    return ref.derivation.answer


def emission_density(ref: Any) -> np.ndarray:
    """``(G, n_r)``: q_g/(4π) on the solve's knots."""
    return np.asarray(answer(ref).emission_density)


def oracle(ref: Any, group: int) -> Any:
    """Group ``group``'s chord oracle on the solve's knots and the solve's own angular rule."""
    return answer(ref).rays[group].oracle


def with_reading(ref: Any, *, angular_points_per_piece: int, ray_points_per_segment: int) -> Any:
    """The same reference with a cheaper READING (the solve's quadrature unchanged; a fresh, unsolved derivation)."""
    import dataclasses

    from orpheus.reference.solution import ReferenceSolution

    reading = module().ReadingQuadrature(
        angular_points_per_piece=angular_points_per_piece, ray_points_per_segment=ray_points_per_segment,
    )
    return ReferenceSolution(ref.specification, dataclasses.replace(ref.derivation, quadrature=reading), None)


def derivation(specification: Any, coord: CoordSystem | None = None) -> Any:
    """The derivation's PUBLIC constructor, called directly (no factory), at the gates' parameters."""
    c = coord if coord is not None else specification.geometry.coord
    return derivation_class()(specification, QUADRATURE[c], max_iter=MAX_ITER, tol=TOL, initial_k=INITIAL_K)

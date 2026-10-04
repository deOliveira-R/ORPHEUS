r"""The trajectory-resolvent sphere and cylinder as reference solutions (#405 P2 step 7b.2.2).

What it is
----------
A :class:`~orpheus.reference.solution.ReferenceSolution` answering the k
question on a layered solid sphere or cylinder with a specular outer law,
through the multi-region trajectory resolvent (Variant α,
:class:`~.billiard.Billiard`). Its family derives no bound on its distance to
the exact answer (#566; the cylinder also #516), so it carries no
certificate and every reading is
:class:`~orpheus.reference.reading.Uncertified` (the user's ruling of
2026-10-03, step 7b).

The reading is the natural extension
------------------------------------
The solve converges to an eigenpair on its radial nodes. Its STATE is the
source its final iterate was transported from: the isotropic emission
density per steradian at the nodes,
:func:`~.greens_function.emission_density`,

.. math::

   \frac{q_g(r_i)}{4\pi} = \frac{1}{4\pi}\Big[\sum_{g'} \Sigma_{s,g'\to g}(r_i)\,\phi_{g'}(r_i)
       + \chi_g(r_i)\,\frac{\sum_{g'} \nu\Sigma_{f,g'}(r_i)\,\phi_{g'}(r_i)}{k}\Big],

kept by the solve with the fission rate it normalised that iterate by
(``last_emission_density``, ``last_fission_rate``), and never the nodal
:math:`\phi` itself (the user's ruling 1 of 2026-10-03: one reference, one
answer). The scalar flux anywhere is the transport integral of that density,
angle-integrated, in the solve's gauge,

.. math::

   \phi_g^{\rm ext}(r) = \frac{1}{F}\int_{4\pi} \big(K_\alpha\,q_g/4\pi\big)(r, \Omega)\,d\Omega,

where :math:`K_\alpha` is the solver's own operator: the chord oracle's
body, called with the evaluation radii (``at=``) apart from the density's
knots (ONE transport; a second ray integrator written for the reading would
be a twin path, X4). The density between knots is the oracle's per-region
cubic spline (ERR-090). At the knots, under the solve's own angular rule,
the extension is the solve's final angular flux, bit for bit.

Ordinary float quadrature (the user's ruling 2 of 2026-10-03)
-------------------------------------------------------------
The orders are one value, :class:`ReadingQuadrature`.

* **The angular integral** is split at every direction whose ray is tangent
  to a knot sphere or an interface below :math:`r` (a tangency is a
  square-root singularity of the integrand in the angle). On the sphere:
  :math:`\mu = \pm\sqrt{1 - (\rho/r)^2}` and 0, Gauss–Legendre in
  :math:`\mu` per piece. On the cylinder: the azimuths :math:`\varphi` with
  :math:`r\,|\sin\varphi| = \rho` and the multiples of :math:`\pi/2`,
  Gauss–Legendre per piece; the axial cosine is integrated as
  :math:`\mu = \cos\theta`, Gauss–Legendre in :math:`\theta`, because the
  3-D lift makes the integrand smooth in :math:`\sin\theta` and so singular
  at :math:`\mu = \pm 1` (a square root in :math:`\mu`). Each square-root
  singularity sits at a piece end, where Gauss–Legendre converges only
  algebraically: `[M]` (qa of 7b.2.2) 16 → 32 → 64 points per piece move a
  sphere point value by 2.9e-7, 3.4e-8, 9e-11, a cylinder one by 1.25e-6,
  6.1e-8, 1.9e-9.
* **The radial integral** of a flux integral, in the geometry's measure
  (:func:`~orpheus.geometry.coord.compute_areas_1d`: :math:`4\pi r^2`,
  :math:`2\pi r` per unit height), is Gauss–Legendre per piece, split at the
  geometry's breakpoints, the knots (the extension's own singular radii) and
  the weight's steps
  (:meth:`~orpheus.numerics.mesh_free_function.Symbolic.steps`, one
  definition with production).

Assumptions, and where they fail
--------------------------------
* Isotropic scattering, no (n,2n): the solver reads :math:`\Sigma_{s,0}`
  alone and no (n,2n) matrix, so a material carrying a higher scattering
  moment or an (n,2n) reaction is refused at construction (answering it
  would answer another problem).
* A solid layered sphere or cylinder (``Billiard``'s multi-region arms),
  specular outer law: every other body is refused (a scope boundary,
  :class:`TrajectoryResolventDerivation`).
* The emission density between knots is a cubic spline per region: its
  error is the solver's discretisation error, unbounded today (#566).
* An eigen answer's representative is gauged by the solve, so its readings
  carry meaning only as ratios; a uniform scaling of the emission density is
  invisible by design.
"""

from __future__ import annotations

import functools
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace
from functools import cached_property
from typing import Any, assert_never

import numpy as np

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.derivations.continuous.trajectory_resolvent.billiard import Billiard
from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    MultiRegionCylinderChordOracle,
    MultiRegionSphereChordOracle,
)
from orpheus.geometry.coord import compute_areas_1d
from orpheus.numerics.content import FrozenMapping
from orpheus.numerics.mesh_free_function import MeshFreeFunction, RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, Linear, PointValue
from orpheus.numerics.question import Eigen
from orpheus.reference.reading import Uncertified
from orpheus.reference.solution import ReferenceSolution
from orpheus.specification.specification import GeometrySpecification

__all__ = ["ReadingQuadrature", "TrajectoryResolventDerivation", "trajectory_resolvent_reference"]

RadialWeight = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class ReadingQuadrature:
    """The orders of the reading's float quadrature (the solve keeps its own ``n_traj_quad``).

    ``ray_points_per_segment``: Gauss–Legendre points along each ray segment;
    ``angular_points_per_piece`` and ``radial_points_per_piece``: points per
    angular and radial piece (the module docstring).
    """

    ray_points_per_segment: int = 64
    angular_points_per_piece: int = 16
    radial_points_per_piece: int = 8


def _gauss_legendre_pieces(breaks: np.ndarray, points_per_piece: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss–Legendre nodes and weights on each piece between consecutive ``breaks``."""
    x, w = np.polynomial.legendre.leggauss(points_per_piece)
    lows, highs = breaks[:-1], breaks[1:]
    half = 0.5 * (highs - lows)
    nodes = (lows[:, None] + half[:, None] * (x[None, :] + 1.0)).ravel()
    weights = (half[:, None] * w[None, :]).ravel()
    return nodes, weights


def _tangency_radii(oracle: MultiRegionSphereChordOracle | MultiRegionCylinderChordOracle, r: float) -> np.ndarray:
    """The knot spheres and interfaces below ``r``: the radii a ray from ``r`` can graze."""
    singular = np.concatenate([oracle.r_nodes, oracle.radii[:-1]])
    return np.unique(singular[(singular > 0.0) & (singular < r)])


@dataclass(frozen=True)
class _SphereRays:
    r"""One group's transport on the sphere, and its direction fibre: :math:`\mu \in [-1, 1]`, :math:`d\Omega = 2\pi\,d\mu`."""

    oracle: MultiRegionSphereChordOracle

    @classmethod
    def per_group(cls, raw: Any, billiard: Billiard) -> tuple[_SphereRays, ...]:
        """One per group, on the solve's knots and angular rule."""
        radii = np.asarray(billiard.geometry_payload["radii"], dtype=float)
        sigma_t = np.asarray(billiard.xs_payload["sigma_t"], dtype=float)
        return tuple(
            cls(MultiRegionSphereChordOracle(
                r_nodes=raw.r_nodes, mu_nodes=raw.mu_nodes, R=float(radii[-1]), radii=radii,
                sigma_t_per_region=sigma_t[:, g], alpha=billiard.alpha_payload["alpha"],
                region_at_node=raw.region_at_node,
            ))
            for g in range(sigma_t.shape[1])
        )

    def scalar_flux(self, density: np.ndarray, r: float, quadrature: ReadingQuadrature) -> float:
        """:math:`2\\pi\\int \\psi(r, \\mu)\\,d\\mu`, the :math:`\\mu` axis split at every tangency."""
        grazing = np.sqrt(1.0 - (_tangency_radii(self.oracle, r) / r) ** 2) if r > 0.0 else np.empty(0)
        breaks = np.unique(np.concatenate([[-1.0, 0.0, 1.0], grazing, -grazing]))
        mu, w_mu = _gauss_legendre_pieces(breaks, quadrature.angular_points_per_piece)
        rays = replace(self.oracle, mu_nodes=mu)
        psi = rays.apply_operator(density, 0.0, n_traj_quad=quadrature.ray_points_per_segment, at=np.array([r]))
        return float(2.0 * np.pi * (psi[0] @ w_mu))


@dataclass(frozen=True)
class _CylinderRays:
    r"""One group's transport on the cylinder, and its direction fibre: :math:`\mu = \cos\theta` and the azimuth :math:`\varphi`."""

    oracle: MultiRegionCylinderChordOracle

    @classmethod
    def per_group(cls, raw: Any, billiard: Billiard) -> tuple[_CylinderRays, ...]:
        """One per group, on the solve's knots and angular rule."""
        radii = np.asarray(billiard.geometry_payload["radii"], dtype=float)
        sigma_t = np.asarray(billiard.xs_payload["sigma_t"], dtype=float)
        return tuple(
            cls(MultiRegionCylinderChordOracle(
                r_nodes=raw.r_nodes, mu_axial_nodes=raw.mu_axial_nodes, phi_az_nodes=raw.phi_az_nodes,
                R=float(radii[-1]), radii=radii, sigma_t_per_region=sigma_t[:, g],
                alpha=billiard.alpha_payload["alpha"], region_at_node=raw.region_at_node,
            ))
            for g in range(sigma_t.shape[1])
        )

    def scalar_flux(self, density: np.ndarray, r: float, quadrature: ReadingQuadrature) -> float:
        """:math:`\\int_0^\\pi\\!\\sin\\theta\\,d\\theta\\int_0^{2\\pi}\\!\\psi\\,d\\varphi`, :math:`\\varphi` split at every tangency."""
        grazing = np.arcsin(_tangency_radii(self.oracle, r) / r) if r > 0.0 else np.empty(0)
        quarter_turns = 0.5 * np.pi * np.arange(5)
        breaks = np.unique(np.concatenate([quarter_turns, grazing, np.pi - grazing, np.pi + grazing, 2.0 * np.pi - grazing]))
        azimuth, w_azimuth = _gauss_legendre_pieces(breaks, quadrature.angular_points_per_piece)
        theta, w_theta = _gauss_legendre_pieces(np.array([0.0, np.pi]), quadrature.angular_points_per_piece)
        rays = replace(self.oracle, mu_axial_nodes=np.cos(theta), phi_az_nodes=azimuth)
        psi = rays.apply_operator(density, 0.0, n_traj_quad=quadrature.ray_points_per_segment, at=np.array([r]))
        return float((w_theta * np.sin(theta)) @ psi[0] @ w_azimuth)


_Rays = _SphereRays | _CylinderRays

#: The bodies whose reading is built: Billiard's two multi-region arms, each with its rays class.
_RAYS: dict[str, type[_SphereRays] | type[_CylinderRays]] = {
    "sphere_mr": _SphereRays,
    "cylinder_mr": _CylinderRays,
}


@dataclass(frozen=True)
class _ConvergedAnswer:
    r"""The solve's answer: its eigenvalue, the source its final iterate came from, and that iterate's normaliser.

    ``emission_density[g, i]`` is :math:`q_g(r_i)/(4\pi)`, the source the
    final iterate was transported from; ``fission_rate`` is the total fission
    rate the solve divided that iterate by. ``rays[g]`` is group g's
    transport on the solve's knots and angular rule.
    """

    k: float
    emission_density: np.ndarray
    fission_rate: float
    rays: tuple[_Rays, ...]


@dataclass(frozen=True)
class _NotConverged:
    """A solve that exhausted its budget: there is no answer to read."""

    iterations: int

    def refusal(self) -> RuntimeError:
        return RuntimeError(
            f"the trajectory-resolvent solve did not converge in {self.iterations} iterations, so it has no answer to read"
        )


@dataclass(frozen=True, eq=False)
class TrajectoryResolventDerivation:
    """The multi-region trajectory resolvent's natural extension, solved on the first evaluation and read uncertified.

    Built from the specification it answers, the solve's ``Billiard``
    quadrature (``n_r``, ``n_mu`` or ``n_mu_axial`` and ``n_phi_az``,
    ``n_traj_quad``), the power iteration's ``max_iter``, ``tol`` and
    ``initial_k`` (``None``: the solver's default) and the reading's
    :class:`ReadingQuadrature`. Construction is the one door: it refuses an
    infinite medium, another question, a material with anisotropic
    scattering or an (n,2n) reaction, every body ``Billiard`` refuses, and
    every body it serves outside its multi-region arms:

    SCOPE-BOUNDARY[guard] machinery: a rays class per further Billiard arm (homogeneous and hollow sphere, cylinder, annulus; a slab's signed-μ fibre), each over its chord oracle carved with ``at=``.
    ruling: the orchestrator, #405 P2 step 7b.2.2.
    revisit: when P4 certifies another trajectory-resolvent family (#566), or a migrated row needs one.

    Nothing is solved here. The solve runs once, on the first
    :meth:`evaluate`, and is held by a ``cached_property`` (a declared pure
    function of the fields, which writes the instance's ``__dict__``); a solve
    that did not converge is held too, as its refusal, so it is not re-run.
    """

    specification: GeometrySpecification
    solver_quadrature: Mapping[str, int]
    max_iter: int | None
    tol: float | None
    initial_k: float | None
    quadrature: ReadingQuadrature = ReadingQuadrature()
    billiard: Billiard = field(init=False, repr=False)
    rays: type[_SphereRays] | type[_CylinderRays] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        specification = self.specification
        if not isinstance(specification, GeometrySpecification):
            raise TypeError(
                f"trajectory_resolvent_reference answers a GeometrySpecification (rays need a finite geometry), "
                f"got a {type(specification).__name__}: an infinite medium has no geometry"
            )
        k_question = Eigen(CellCoefficient.every(Channel.FISSION_EMISSION).resolve(specification.materials))
        if specification.question != k_question:
            raise ValueError(
                f"trajectory_resolvent_reference answers the fundamental k-eigenvalue question at the physical "
                f"point, got {specification.question!r}"
            )
        object.__setattr__(self, "solver_quadrature", FrozenMapping(self.solver_quadrature.items()))
        billiard = Billiard(specification.geometry, dict(specification.materials.items()), dict(self.solver_quadrature))
        rays = _RAYS.get(billiard.geometry_kind)
        if rays is None:
            raise NotImplementedError(
                f"trajectory_resolvent_reference reads a layered solid sphere or cylinder; the reading of a "
                f"{billiard.geometry_kind} body is not built (a scope boundary: its rays class over its chord "
                f"oracle carved with at=)"
            )
        object.__setattr__(self, "billiard", billiard)
        object.__setattr__(self, "rays", rays)

    @cached_property
    def _solve(self) -> _ConvergedAnswer | _NotConverged:
        """The solve, once, converged or not."""
        solution = self.billiard.solve_critical(max_iter=self.max_iter, tol=self.tol, initial_k=self.initial_k)
        raw = solution.metadata["raw_result"]
        if not solution.converged:
            return _NotConverged(int(raw.iterations))
        return _ConvergedAnswer(
            float(raw.k_eff), raw.last_emission_density, float(raw.last_fission_rate), self.rays.per_group(raw, self.billiard),
        )

    @property
    def answer(self) -> _ConvergedAnswer:
        """The converged solve; a solve that did not converge has no answer to read."""
        match self._solve:
            case _ConvergedAnswer() as answer:
                return answer
            case _NotConverged() as exhausted:
                raise exhausted.refusal()
            case unreachable:
                assert_never(unreachable)

    @cached_property
    def scalar_flux(self) -> Callable[[int, float], float]:
        r""":math:`\phi_g^{\rm ext}(r)` as ``scalar_flux(g, r)``, memoised per ``(g, r)`` for the derivation's lifetime."""
        answer = self.answer

        @functools.cache
        def extension(group: int, r: float) -> float:
            transported = answer.rays[group].scalar_flux(answer.emission_density[group], r, self.quadrature)
            return transported / answer.fission_rate

        return extension

    def evaluate(self, observable: Eigenvalue | Linear) -> Uncertified:
        """The observable's value, uncertified: the family derives no bound (#566)."""
        match observable:
            case Eigenvalue():
                return Uncertified(self.answer.k)
            case PointValue(position=position, group=group):
                return Uncertified(self.scalar_flux(group, position))
            case FluxIntegral(weight=weight):
                return Uncertified(self._flux_integral(weight))
            case _:
                assert_never(observable)

    def _flux_integral(self, weight: MeshFreeFunction) -> float:
        r""":math:`\sum_g \int w_g\,\phi_g^{\rm ext}\,dV` by Gauss–Legendre per radial piece, in the geometry's measure."""
        geometry = self.billiard.geometry
        breakpoints = np.asarray(geometry.breakpoints, dtype=float)
        radial_weights, steps = _radial_weight(weight, breakpoints)
        knots = self.answer.rays[0].oracle.r_nodes
        breaks = np.unique(np.concatenate([breakpoints, knots, steps]))
        breaks = breaks[(breaks >= breakpoints[0]) & (breaks <= breakpoints[-1])]
        radii, w_radial = _gauss_legendre_pieces(breaks, self.quadrature.radial_points_per_piece)
        volume_element = compute_areas_1d(geometry.coord, radii) * w_radial  # dV at each node, cm³ (per cm of height)
        total = 0.0
        for group, radial_weight in enumerate(radial_weights):
            if radial_weight is None:
                continue
            w_g = radial_weight(radii)
            supported = np.flatnonzero(w_g != 0.0)  # the scalar flux is read only where the weight is not zero
            flux = np.array([self.scalar_flux(group, r) for r in radii[supported].tolist()])
            total += float(np.sum(w_g[supported] * flux * volume_element[supported]))
        return total


def _radial_weight(weight: MeshFreeFunction, breakpoints: np.ndarray) -> tuple[list[RadialWeight | None], tuple[float, ...]]:
    """Each group's weight as a function of the radius (``None`` where it is zero), and its steps on the geometry."""
    match weight:
        case RegionwiseConstant(values=values):
            def region_value(column: np.ndarray) -> RadialWeight:
                return lambda r: column[np.searchsorted(breakpoints[1:-1], r, side="right")]

            return [region_value(values[:, g]) if np.any(values[:, g]) else None for g in range(weight.n_groups)], ()
        case Symbolic():
            import sympy

            scalar = weight.without(Symbolic.mu, Symbolic.phi)
            steps = scalar.steps((float(breakpoints[0]), float(breakpoints[-1])))

            def evaluated(expression: sympy.Expr) -> RadialWeight:
                function = sympy.lambdify(Symbolic.r, expression, "numpy")
                return lambda r: np.broadcast_to(np.asarray(function(r), dtype=float), r.shape)

            return [None if e == 0 else evaluated(e) for e in scalar.expressions], steps
        case _:
            assert_never(weight)


def trajectory_resolvent_reference(
    specification: GeometrySpecification,
    quadrature: Mapping[str, int],
    *,
    max_iter: int | None = None,
    tol: float | None = None,
    initial_k: float | None = None,
) -> ReferenceSolution:
    r"""The trajectory resolvent's :class:`~orpheus.reference.solution.ReferenceSolution`, uncertified and lazy.

    It answers the fundamental k-eigenvalue question at the physical point
    (:math:`\mathrm{Eigen}` along every fission emission) on a layered solid
    sphere or cylinder; ``quadrature`` is the solve's (``Billiard``'s keys).
    The refusals are the derivation's (:class:`TrajectoryResolventDerivation`),
    all before any solve: nothing is solved here.
    """
    derivation = TrajectoryResolventDerivation(specification, quadrature, max_iter, tol, initial_k)
    return ReferenceSolution(specification, derivation, None)

r"""Mathematical billiards — phase-space dynamics with boundary reflection.

This module exposes :class:`Billiard`, the math-rich frame that names
the structure shared by every Variant α reference solver in
:mod:`~orpheus.derivations.continuous.trajectory_resolvent`. The
implementation is a thin facade over the geometry-specific
``solve_greens_function_*`` entry points; what it adds is a single
named *concept* (the billiard system) that covers all six geometries
(slab, sphere, cylinder, slab-asymmetric, hollow sphere, annulus) at
both orbit-space classes (one-surface compact and two-surface).

Mathematical billiards in one paragraph
---------------------------------------

A *billiard* in the Birkhoff–Sinai sense is a dynamical system on a
domain :math:`M \subset \mathbb R^d` (the "billiard table") with a
specified reflection law on :math:`\partial M`. A particle moves
freely between collisions and reflects at the boundary. The flow is
deterministic, time-reversible, volume-preserving (Liouville), and
possibly ergodic. When :math:`M` is bounded and the reflection is
specular (angle of incidence equals angle of reflection), the
billiard is *closed*; when the reflection has albedo :math:`\alpha
\in [0, 1)`, the billiard is *absorbing* and orbits decay
geometrically with rate :math:`\alpha`. Linear (neutron) transport in
a homogeneous geometry IS a mathematical billiard with the
*streaming* flow on the interior and a *partial-albedo* reflection
law on :math:`\partial M`. The Variant α resolvent
:math:`T = (I - S)^{-1}` is the Birkhoff transfer operator's
resolvent for that billiard, summed in closed form via the Neumann
series :math:`\sum_n (\alpha\, e^{-\tau})^n`.

Architectural role (Phase D)
----------------------------

Per the architectural reset, ``Billiard`` is a **reference solution
generator** that consumes a :class:`StructuredGeometry` directly via
its ``__init__``. It is *not* a production solver and does not share
a Protocol with the discrete CP/SN/MOC family — those consume
``(materials, mesh, params)`` via canonical free functions
(``solve_cp`` / ``solve_sn`` / ``solve_moc``). The architectural split
is now full: reference generators build the truth values; production
solvers chew through 1000-group industrial data. Each pillar consumes
its own geometry/mesh layer cleanly.

Why "Billiard"
--------------

The class name encodes a piece of mathematics every reader benefits
from learning: that **the rank-1 / rank-2 closure formula in
:mod:`.variant_alpha_core` is not a transport-physics special trick
but the analytic resolvent of a Birkhoff transfer operator on a
billiard table**. The same operator algebra appears in classical
mechanics (gas-of-particles in a box), in semiclassical quantum
chaos (Gutzwiller trace formula on hyperbolic billiards), and in
ergodic theory (Sinai's dispersing billiards). The transport-physics
specialization we ship here is one application of a much wider
mathematical theory; naming the class :class:`Billiard` makes that
visible at every call site.

References
----------

- Birkhoff, G.D. (1927). *Dynamical Systems*, Ch. 6 (billiard
  flows). Reissued by AMS Colloquium Publications, 1966.
- Sinai, Ya.G. (1970). "Dynamical systems with elastic reflections
  and their applications." *Russ. Math. Surveys* **25**, 137–189.
  DOI: 10.1070/RM1970v025n02ABEH003794.
- Bunimovich, L.A. (1979). "On the ergodic properties of nowhere
  dispersing billiards." *Commun. Math. Phys.* **65**, 295–312.
- Chernov, N. and Markarian, R. (2006). *Chaotic Billiards*, AMS
  Mathematical Surveys and Monographs **127**.
- Sanchez, R. (2002). "Treatment of Boundary Conditions in
  Trajectory-Based Deterministic Transport Methods." *Nucl. Sci.
  Eng.* **140**, 23–50. The transport-physics specialization of
  billiard theory we implement here. The Variant α resolvent
  :math:`T = (I - S)^{-1}` is Sanchez's Eq. (A4) in our notation.
- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* **14**.
  DOI: 10.1080/00411458608210456. The :math:`\omega_1 = 0`
  (isotropic-scattering) precursor.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry.structured_geometry import StructuredGeometry

# Use the SHARED cross-method result types — these define the contract
# every math-heart class (MomentSpace for fn_method, Billiard for
# trajectory_resolvent) populates. Method-specific rich results live
# under ``metadata["raw_result"]`` (the legacy result dataclass) and
# ``metadata["psi"]``, ``metadata["phi"]``, ``metadata["iterations"]``,
# ``metadata["mesh"]`` for direct field access.
from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)


# Re-export so consumers can do
#   from orpheus.derivations.continuous.trajectory_resolvent import (
#       Billiard, CriticalSolution, FluxSolution,
#   )
# even though the canonical home is ``derivations.common.solution_types``.
__all__ = ["Billiard", "CriticalSolution", "FluxSolution"]


# ─────────────────────────────────────────────────────────────────────
# Result-shape contract
# ─────────────────────────────────────────────────────────────────────
#
# Per the schema in ``derivations.common.solution_types``, the unified
# CriticalSolution is small (eigenvalue, eigenvalue_kind,
# parameter_value, parameter_kind, converged, metadata). Method-specific
# rich data — the angular flux ``psi``, scalar flux ``phi``, iteration
# count, and mesh metadata — lives in ``metadata`` under the canonical
# keys:
#
# - ``metadata["raw_result"]`` — the legacy per-geometry result
#   dataclass (e.g. ``GreensFunctionResult``); preserves bit-equal
#   access to every original field.
# - ``metadata["psi"]`` / ``metadata["phi"]`` — direct convenience
#   handles to the angular and scalar flux arrays.
# - ``metadata["iterations"]`` — power-iteration step count.
# - ``metadata["mesh"]`` — dict of grid arrays (``r_nodes``, ``x_nodes``,
#   ``mu_nodes``, ``mu_axial_nodes``, ``phi_az_nodes``,
#   ``region_at_node``).
# - ``metadata["n_groups"]`` — group count (1 for 1G; G for MG / MR).
# - ``metadata["geometry_kind"]`` — the Billiard's geometry tag, for
#   downstream dispatch.
#
# For Variant α power iteration, ``parameter_kind = "fixed_geometry"``:
# the geometry is fixed by the caller (we don't root-find a critical
# parameter), so the parameter_value field carries the geometry's
# characteristic length (``R`` for sphere/cylinder, ``L`` for slab,
# ``R_out`` for hollow geometries). Cross-method comparators that
# need the eigenvalue at this configuration read it from
# ``solution.eigenvalue`` directly.


# ─────────────────────────────────────────────────────────────────────
# Billiard system
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Billiard:
    r"""A mathematical billiard — phase-space dynamics with boundary
    reflection.

    A *billiard* in the Birkhoff–Sinai sense is a dynamical system
    consisting of:

    1. **A domain :math:`M`** (the "billiard table") — here the
       spatial geometry (slab, sphere, cylinder, hollow sphere,
       annulus, asymmetric slab; classified by the orbit-space M/G
       structure).
    2. **A boundary reflection law** at :math:`\partial M` — here
       the BC parametrized by :math:`\alpha`:

       - specular (:math:`\alpha = 1`) ↔ pure billiard, ergodic at
         high density;
       - vacuum (:math:`\alpha = 0`) ↔ open billiard, particles
         escape;
       - partial-albedo (:math:`\alpha \in (0, 1)`) ↔ scattering
         billiard (energy / particle loss at boundary).
    3. **Free flight inside the domain**, governed by the streaming
       operator :math:`\Omega \cdot \nabla + \Sigma_t` along
       characteristic curves (rays / geodesics in the appropriate
       metric).
    4. **A transfer operator** :math:`S` that maps boundary-flux
       moments to boundary-flux moments under one full bounce
       period, and its resolvent :math:`T = (I - S)^{-1}`.

    Construction (Phase D)
    ----------------------

    Construct directly with a :class:`StructuredGeometry` plus the
    materials dict::

        from orpheus.geometry.structured_geometry import (
            Region, StructuredGeometry,
        )
        from orpheus.geometry.mesh import BC
        from orpheus.derivations.continuous.trajectory_resolvent import (
            Billiard,
        )

        geom = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=5.0),),
            bcs=(BC.reflective,),
        )
        b = Billiard(
            geometry=geom,
            materials={0: pu_mixture},
            alpha=1.0,
            quadrature={"n_r": 24, "n_mu": 24, "n_traj_quad": 64},
        )
        sol = b.solve_critical()

    Geometry tag mapping
    --------------------

    The class accepts uppercase :class:`StructuredGeometry` tags
    (``"SLB"``, ``"SPH"``, ``"CYL"``) and dispatches internally to
    the same per-geometry ``solve_greens_function_*`` entry points.
    The asymmetric-slab branch is selected when ``alpha`` is a dict
    carrying ``alpha_left`` / ``alpha_right`` keys.

    Attributes
    ----------
    geometry : StructuredGeometry
        The pure-geometry layer object describing the billiard table.
    materials : dict[int, Mixture]
        Production-protocol cross-section payload, keyed by material
        ID. Today's :class:`Billiard` consumes only the active
        material at key ``0`` for single-region geometries.
    alpha : float | dict[str, float]
        Boundary reflectivity. Float for symmetric / one-surface
        billiards; dict for asymmetric two-surface billiards
        (``{"alpha_left": ..., "alpha_right": ...}``).
    quadrature : dict[str, int]
        Quadrature orders. Standard keys: ``n_r`` / ``n_x``,
        ``n_mu`` / ``n_mu_axial``, ``n_phi_az``, ``n_traj_quad``.
        Defaults are inherited from each geometry's solver if the
        key is absent.
    geometry_kind : str
        One of ``"sphere"``, ``"cylinder"``, ``"slab"``,
        ``"slab_asymmetric"``. Auto-derived in :meth:`__post_init__`
        from the :class:`StructuredGeometry` tag and the alpha-dict
        shape; selects which underlying solver is dispatched.
    closure_rank : int
        ``1`` for one-endpoint orbit spaces, ``2`` for two-endpoint
        orbit spaces. Auto-resolved from :attr:`geometry_kind`.

    Examples
    --------
    Closed homogeneous sphere with k_eff = k_inf::

        geom = StructuredGeometry(
            geometry="SPH",
            regions=(Region(mat_id=0, outer_thickness_cm=5.0),),
            bcs=(BC.reflective,),
        )
        b = Billiard(geometry=geom, materials={0: pu_mixture}, alpha=1.0)
        sol = b.solve_critical()
    """

    geometry: StructuredGeometry
    materials: dict[int, Mixture]
    alpha: float | dict[str, float] = 1.0
    quadrature: dict[str, int] = field(default_factory=dict)

    # Derived fields populated in __post_init__. Marked as init=False so
    # they're owned by the dataclass but not user-supplied.
    geometry_kind: str = field(init=False)
    closure_rank: int = field(init=False)
    xs_payload: dict[str, Any] = field(init=False)
    geometry_payload: dict[str, Any] = field(init=False)
    alpha_payload: dict[str, float] = field(init=False)

    def __post_init__(self) -> None:
        if not _is_mixture_dict(self.materials):
            raise ValueError(
                "Billiard.materials must be a dict[int, Mixture]; got keys "
                f"{list((self.materials or {}).keys())!r}."
            )

        geometry_kind = _infer_geometry_kind(self.geometry, self.alpha)
        geometry_payload = _geometry_payload_for_solver(
            geometry_kind, self.geometry
        )
        xs_payload = _mixture_to_solver_xs_payload(
            self.materials, geometry_kind
        )

        # Normalize alpha into the per-geometry payload.
        if isinstance(self.alpha, dict):
            alpha_payload: dict[str, float] = dict(self.alpha)
        elif geometry_kind in ("hollow_sphere", "annulus"):
            alpha_payload = {
                "alpha_in": float(self.alpha),
                "alpha_out": float(self.alpha),
            }
        else:
            alpha_payload = {"alpha": float(self.alpha)}

        # Resolve closure rank from orbit-space class.
        closure_rank = 2 if geometry_kind == "slab_asymmetric" else 1

        # Frozen dataclass: use object.__setattr__ to populate derived
        # fields.
        object.__setattr__(self, "geometry_kind", geometry_kind)
        object.__setattr__(self, "closure_rank", closure_rank)
        object.__setattr__(self, "xs_payload", xs_payload)
        object.__setattr__(self, "geometry_payload", geometry_payload)
        object.__setattr__(self, "alpha_payload", alpha_payload)

    def with_alpha(self, alpha: float | dict[str, float]) -> "Billiard":
        """Return a copy with a different boundary reflectivity."""
        return Billiard(
            geometry=self.geometry,
            materials=self.materials,
            alpha=alpha,
            quadrature=self.quadrature,
        )

    # ─────────────────────────────────────────────────────────────────
    # solve methods — thin facades over the existing solvers
    # ─────────────────────────────────────────────────────────────────

    def solve_critical(
        self,
        *,
        max_iter: int | None = None,
        tol: float | None = None,
        initial_psi: np.ndarray | None = None,
        initial_k: float | None = None,
    ) -> CriticalSolution:
        r"""Solve the *k*-eigenproblem on this billiard.

        Dispatches to the geometry-specific
        ``solve_greens_function_*`` entry point and re-packs the
        result into the SHARED cross-method
        :class:`~orpheus.derivations.common.solution_types.CriticalSolution`.

        Parameters
        ----------
        max_iter : int, optional
            Maximum power-iteration count. Defaults to the
            underlying solver's default (200 for 1G, 300 for MG,
            500 for MR).
        tol : float, optional
            Relative-:math:`k` convergence tolerance. Defaults to
            ``1e-10`` for 1G, ``1e-9`` for MG / MR.
        initial_psi : ndarray, optional
            Initial angular flux iterate. Shape must match the
            geometry's phase-space discretization.
        initial_k : float, optional
            Initial :math:`k_{\rm eff}` estimate. Defaults to
            ``k_inf`` from the homogenized infinite-medium balance.

        Returns
        -------
        CriticalSolution
        """
        return _dispatch_critical(
            self,
            max_iter=max_iter,
            tol=tol,
            initial_psi=initial_psi,
            initial_k=initial_k,
        )

    def solve_fixed_source(
        self,
        external_source: np.ndarray,
        *,
        max_iter: int | None = None,
        tol: float | None = None,
    ) -> FluxSolution:
        r"""Solve the fixed-source problem on this billiard.

        Currently only supported for :attr:`geometry_kind` ==
        ``"sphere_mr"`` (the Garcia 2021 multi-region sphere
        benchmark family). Other geometries raise
        :class:`NotImplementedError`.
        """
        return _dispatch_fixed_source(
            self,
            external_source=external_source,
            max_iter=max_iter,
            tol=tol,
        )


# ─────────────────────────────────────────────────────────────────────
# Internal dispatchers — bridge from Billiard to the per-geometry solvers
# ─────────────────────────────────────────────────────────────────────


def _filter_kwargs(d: dict[str, Any], allowed: tuple[str, ...]) -> dict[str, Any]:
    """Return only the allowed keys (drops unspecified defaults)."""
    return {k: d[k] for k in allowed if k in d}


def _dispatch_critical(
    billiard: Billiard,
    *,
    max_iter: int | None,
    tol: float | None,
    initial_psi: np.ndarray | None,
    initial_k: float | None,
) -> CriticalSolution:
    """Route to the appropriate ``solve_greens_function_*`` entry."""
    # Lazy imports — keep Billiard importable from doc-build contexts
    # that may not yet have all geometry modules wired.
    from orpheus.derivations.continuous.trajectory_resolvent import (
        greens_function as gf_sphere,
        greens_function_cylinder as gf_cyl,
        greens_function_slab as gf_slab,
        greens_function_slab_asymmetric as gf_slab_asym,
        greens_function_hollow_sphere as gf_hollow,
        greens_function_annulus as gf_annulus,
    )

    g = billiard.geometry_kind
    mats = billiard.xs_payload
    geom = billiard.geometry_payload
    a = billiard.alpha_payload
    q = billiard.quadrature

    # Build common iteration kwargs (skip None — let the underlying
    # solver use its own defaults).
    iter_kwargs: dict[str, Any] = {}
    if max_iter is not None:
        iter_kwargs["max_iter"] = max_iter
    if tol is not None:
        iter_kwargs["tol"] = tol
    if initial_psi is not None:
        iter_kwargs["initial_psi"] = initial_psi

    # Determine multi-group from materials shape.
    sigma_t = np.asarray(mats["sigma_t"])
    is_mg = sigma_t.ndim >= 1 and sigma_t.size > 1

    # Resolve the configuration parameter the eigenvalue is associated
    # with. Variant α is "fixed_geometry" — the user picks the geometry
    # and we report k_eff at that fixed config (no critical root-find).
    parameter_kind = "fixed_geometry"
    if g in ("sphere", "cylinder"):
        parameter_value = float(geom["R"])
    elif g in ("slab", "slab_asymmetric"):
        parameter_value = float(geom["L"])
    elif g in ("hollow_sphere", "annulus"):
        parameter_value = float(geom["R_out"])
    elif g == "sphere_mr":
        parameter_value = float(np.asarray(geom["radii"])[-1])
    else:
        parameter_value = 0.0
    wrap_extra = dict(
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=g,
    )

    if g == "sphere":
        if is_mg:
            res = gf_sphere.solve_greens_function_sphere_mg(
                R=geom["R"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_sphere.solve_greens_function_sphere(
                R=geom["R"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_1g(res, n_groups=1, **wrap_extra)

    if g == "sphere_mr":
        res = gf_sphere.solve_greens_function_sphere_mr(
            radii=geom["radii"],
            sigma_t=mats["sigma_t"],
            sigma_s=mats["sigma_s"],
            nu_sigma_f=mats["nu_sigma_f"],
            chi=mats.get("chi"),
            alpha=a["alpha"],
            **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
            **iter_kwargs,
            **({"initial_k": initial_k} if initial_k is not None else {}),
        )
        # Multi-region always carries (G, n_r, n_mu) shape — psi_g.
        n_groups = int(np.atleast_2d(np.asarray(mats["sigma_t"])).shape[-1]
                       if np.asarray(mats["sigma_t"]).ndim > 1
                       else 1)
        return _wrap_mr(res, n_groups=n_groups, **wrap_extra)

    if g == "cylinder":
        if is_mg:
            res = gf_cyl.solve_greens_function_cylinder_mg(
                R=geom["R"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha=a["alpha"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_cyl_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_cyl.solve_greens_function_cylinder(
                R=geom["R"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha=a["alpha"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
            )
            return _wrap_cyl_1g(res, n_groups=1, **wrap_extra)

    if g == "slab":
        if is_mg:
            res = gf_slab.solve_greens_function_slab_mg(
                L=geom["L"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_slab_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_slab.solve_greens_function_slab(
                L=geom["L"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha=a["alpha"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_slab_1g(res, n_groups=1, **wrap_extra)

    if g == "slab_asymmetric":
        if is_mg:
            res = gf_slab_asym.solve_greens_function_slab_asymmetric_mg(
                L=geom["L"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha_left=a["alpha_left"],
                alpha_right=a["alpha_right"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_slab_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_slab_asym.solve_greens_function_slab_asymmetric(
                L=geom["L"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha_left=a["alpha_left"],
                alpha_right=a["alpha_right"],
                **_filter_kwargs(q, ("n_x", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_slab_1g(res, n_groups=1, **wrap_extra)

    if g == "hollow_sphere":
        if is_mg:
            res = gf_hollow.solve_greens_function_hollow_sphere_mg(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_hollow.solve_greens_function_hollow_sphere(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
                **iter_kwargs,
            )
            return _wrap_1g(res, n_groups=1, **wrap_extra)

    if g == "annulus":
        if is_mg:
            res = gf_annulus.solve_greens_function_annulus_mg(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=mats["sigma_t"],
                sigma_s=mats["sigma_s"],
                nu_sigma_f=mats["nu_sigma_f"],
                chi=mats.get("chi"),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
                **({"initial_k": initial_k} if initial_k is not None else {}),
            )
            return _wrap_cyl_mg(res, n_groups=int(sigma_t.size), **wrap_extra)
        else:
            res = gf_annulus.solve_greens_function_annulus(
                R_in=geom["R_in"],
                R_out=geom["R_out"],
                sigma_t=float(sigma_t),
                sigma_s=float(np.asarray(mats["sigma_s"])),
                nu_sigma_f=float(np.asarray(mats["nu_sigma_f"])),
                alpha_in=a["alpha_in"],
                alpha_out=a["alpha_out"],
                **_filter_kwargs(q, (
                    "n_r", "n_mu_axial", "n_phi_az", "n_traj_quad",
                )),
                **iter_kwargs,
            )
            return _wrap_cyl_1g(res, n_groups=1, **wrap_extra)

    raise NotImplementedError(
        f"solve_critical not implemented for geometry_kind={g!r}"
    )


def _dispatch_fixed_source(
    billiard: Billiard,
    *,
    external_source: np.ndarray,
    max_iter: int | None,
    tol: float | None,
) -> FluxSolution:
    """Route fixed-source solves to the appropriate solver."""
    from orpheus.derivations.continuous.trajectory_resolvent import (
        greens_function as gf_sphere,
    )

    g = billiard.geometry_kind
    mats = billiard.xs_payload
    geom = billiard.geometry_payload
    a = billiard.alpha_payload
    q = billiard.quadrature

    iter_kwargs: dict[str, Any] = {}
    if max_iter is not None:
        iter_kwargs["max_iter"] = max_iter
    if tol is not None:
        iter_kwargs["tol"] = tol

    if g == "sphere_mr":
        res = gf_sphere.solve_greens_function_sphere_mr_fixed_source(
            radii=geom["radii"],
            sigma_t=mats["sigma_t"],
            sigma_s=mats["sigma_s"],
            external_source=external_source,
            alpha=a["alpha"],
            **_filter_kwargs(q, ("n_r", "n_mu", "n_traj_quad")),
            **iter_kwargs,
        )
        n_groups = int(
            np.asarray(external_source).reshape(-1, 1).shape[1]
            if np.asarray(external_source).ndim > 1
            else 1
        )
        # Collapse phi_g per group to a representative scalar flux for
        # the shared FluxSolution `scalar_flux` field.
        phi_for_shared = (
            res.phi_g if res.phi_g.ndim == 1 else res.phi_g[0]
        )
        metadata = {
            "raw_result": res,
            "psi": res.psi_g,
            "phi": res.phi_g,
            "phi_g": res.phi_g,
            "iterations": int(res.iterations),
            "converged": bool(res.converged),
            "mesh": {
                "r_nodes": res.r_nodes,
                "mu_nodes": res.mu_nodes,
                "region_at_node": res.region_at_node,
            },
            "n_groups": n_groups,
            "geometry_kind": g,
        }
        return FluxSolution(
            spatial_nodes=res.r_nodes,
            scalar_flux=phi_for_shared,
            angular_flux=res.psi_g,
            angular_nodes=res.mu_nodes,
            eigenvalue=0.0,
            eigenvalue_kind="none",
            spatial_units="cm",
            metadata=metadata,
        )

    raise NotImplementedError(
        f"solve_fixed_source not yet implemented for "
        f"geometry_kind={g!r}. Currently only 'sphere_mr' is supported."
    )


# ─────────────────────────────────────────────────────────────────────
# Result-shape adaptors — populate the SHARED CriticalSolution schema
# ─────────────────────────────────────────────────────────────────────


def _build_critical(
    res: Any,
    *,
    psi: np.ndarray,
    phi: np.ndarray,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
    mesh: dict[str, Any],
) -> CriticalSolution:
    """Construct the shared CriticalSolution from a legacy result."""
    metadata = {
        "raw_result": res,
        "psi": psi,
        "phi": phi,
        "iterations": int(res.iterations),
        "mesh": mesh,
        "n_groups": int(n_groups),
        "geometry_kind": geometry_kind,
    }
    return CriticalSolution(
        eigenvalue=float(res.k_eff),
        eigenvalue_kind="k_eff",
        parameter_value=float(parameter_value),
        parameter_kind=parameter_kind,
        converged=bool(res.converged),
        metadata=metadata,
    )


def _wrap_1g(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap a 1G geometry-specific result with (r_nodes, mu_nodes)."""
    return _build_critical(
        res,
        psi=res.psi,
        phi=res.phi,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


def _wrap_mg(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap an MG geometry-specific result."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


def _wrap_mr(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap a multi-region sphere result with region_at_node."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_nodes": res.mu_nodes,
            "region_at_node": res.region_at_node,
        },
    )


def _wrap_cyl_1g(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap cylinder / annulus 1G result (3D angular phase-space)."""
    return _build_critical(
        res,
        psi=res.psi,
        phi=res.phi,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_axial_nodes": res.mu_axial_nodes,
            "phi_az_nodes": res.phi_az_nodes,
        },
    )


def _wrap_cyl_mg(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap cylinder / annulus MG result."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "r_nodes": res.r_nodes,
            "mu_axial_nodes": res.mu_axial_nodes,
            "phi_az_nodes": res.phi_az_nodes,
        },
    )


def _wrap_slab_1g(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap slab / slab_asymmetric 1G result."""
    return _build_critical(
        res,
        psi=res.psi,
        phi=res.phi,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "x_nodes": res.x_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


def _wrap_slab_mg(
    res: Any,
    *,
    n_groups: int,
    parameter_value: float,
    parameter_kind: str,
    geometry_kind: str,
) -> CriticalSolution:
    """Wrap slab / slab_asymmetric MG result."""
    return _build_critical(
        res,
        psi=res.psi_g,
        phi=res.phi_g,
        n_groups=n_groups,
        parameter_value=parameter_value,
        parameter_kind=parameter_kind,
        geometry_kind=geometry_kind,
        mesh={
            "x_nodes": res.x_nodes,
            "mu_nodes": res.mu_nodes,
        },
    )


# ─────────────────────────────────────────────────────────────────────
# Internal helpers — translate StructuredGeometry + Mixture inputs
# into the solver-facing payload shapes
# ─────────────────────────────────────────────────────────────────────


def _is_mixture_dict(materials: Any) -> bool:
    """Return True iff *materials* looks like ``dict[int, Mixture]``."""
    if not isinstance(materials, dict) or not materials:
        return False
    first_key = next(iter(materials))
    if not isinstance(first_key, int):
        return False
    return isinstance(materials[first_key], Mixture)


# Map StructuredGeometry uppercase tag → Billiard's internal lowercase
# geometry_kind. The asymmetric-slab branch is selected by alpha-dict
# shape, not by the StructuredGeometry tag (slab is slab — asymmetry is
# a BC concept, not a coordinate concept).
_TAG_TO_KIND: dict[str, str] = {
    "SLB": "slab",
    "SPH": "sphere",
    "CYL": "cylinder",
}


def _infer_geometry_kind(
    geom: StructuredGeometry,
    alpha: float | dict[str, float],
) -> str:
    """Infer Billiard's internal geometry_kind from a StructuredGeometry.

    Returned values are one of: ``"slab"``, ``"slab_asymmetric"``,
    ``"sphere"``, ``"cylinder"``. The asymmetric-slab branch is
    selected when *alpha* is a dict carrying ``alpha_left`` /
    ``alpha_right`` keys.
    """
    tag = geom.geometry
    if tag not in _TAG_TO_KIND:
        raise ValueError(
            f"Billiard cannot construct on StructuredGeometry tag "
            f"{tag!r}; supported: {sorted(_TAG_TO_KIND)}."
        )
    kind = _TAG_TO_KIND[tag]
    if kind == "slab" and isinstance(alpha, dict):
        keys = set(alpha.keys())
        if "alpha_left" in keys or "alpha_right" in keys:
            return "slab_asymmetric"
    return kind


def _geometry_payload_for_solver(
    geometry_kind: str,
    geom: StructuredGeometry,
) -> dict[str, Any]:
    """Build the per-geometry payload the dispatchers consume.

    Maps the :class:`StructuredGeometry` extent onto the
    geometry-specific kwargs of the underlying
    ``solve_greens_function_*`` entry points.

    Slab convention reminder
    ------------------------
    :attr:`StructuredGeometry.domain_extent_cm` is the FULL slab width
    (sum of region thicknesses). The underlying
    :func:`solve_greens_function_slab` consumes ``L`` = full slab
    width, so this is a direct pass-through (no halving / doubling).

    For sphere / cylinder, ``domain_extent_cm`` is the outer radius.
    """
    extent = geom.domain_extent_cm
    if geometry_kind in ("sphere", "cylinder"):
        return {"R": float(extent)}
    if geometry_kind in ("slab", "slab_asymmetric"):
        return {"L": float(extent)}
    raise ValueError(
        f"_geometry_payload_for_solver: unsupported {geometry_kind!r}"
    )


def _mixture_to_solver_xs_payload(
    materials: dict[int, Mixture],
    geometry_kind: str,
) -> dict[str, Any]:
    """Translate a Mixture-keyed dict to the solver-facing XS payload.

    The ``solve_greens_function_*`` entry points consume raw numpy
    arrays (``sigma_t``, ``sigma_s``, ``nu_sigma_f``, ``chi``); this
    helper extracts those from the production-protocol
    :class:`Mixture` at ``materials[0]``.
    """
    if 0 not in materials:
        raise ValueError(
            f"Billiard: materials must contain key 0 "
            f"(the active mat_id). Got keys {sorted(materials.keys())}."
        )
    mix = materials[0]
    sig_t = np.asarray(mix.SigT)
    sig_s_p0 = np.asarray(mix.SigS[0].todense())
    nu_sigma_f = np.asarray(mix.SigP)
    chi = np.asarray(mix.chi)
    if sig_t.size == 1:
        return {
            "sigma_t": float(sig_t[0]),
            "sigma_s": float(sig_s_p0[0, 0]),
            "nu_sigma_f": float(nu_sigma_f[0]),
        }
    return {
        "sigma_t": sig_t,
        "sigma_s": sig_s_p0,
        "nu_sigma_f": nu_sigma_f,
        "chi": chi,
    }

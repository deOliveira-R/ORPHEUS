r"""Solver adapters for the cross-method test protocol.

Each adapter wraps a continuous-reference solver to the
:class:`~tests.cross_method.protocol.SolverAdapter` shape. The
adapter:

* reads cross sections directly off ``case.registry_case.materials``
  / ``case.materials`` (the production-protocol Mixture API);
* selects internal numerical parameters (n_modes for fn_method,
  n_r/n_mu/n_traj_quad for trajectory_resolvent) based on the
  requested case tolerance;
* performs unit conversions (mfp ↔ cm, half-thickness ↔ full
  slab);
* returns a :class:`ScalarResult` with the right ``tag``.

Phase D
-------

The pre-Phase-D ``mixture_to_fn_arrays`` extractor was retired as
part of the architectural reset; adapters now read
``mixture.SigT`` / ``SigS`` / ``SigP`` directly (the same pattern
the math-heart classes Billiard / MomentSpace / Spectrum / BasisSpace
already use after their direct-__init__ migration).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .protocol import CrossMethodCase, ScalarResult


# ═══════════════════════════════════════════════════════════════════
# F_N method adapters (fn_method package)
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class FNSlabAdapter:
    r"""Adapter for :func:`...fn_method.slab.solve_fn_slab_bare_critical`.

    Reports the F_N method's predicted critical half-thickness in
    mean-free paths (the ``a_critical_mfp`` tag). Internally selects
    ``n_modes = 10`` by default — Grandjean-Siewert Table XI shows
    F_10 reaches ~5e-6 absolute on the slab half-thickness across
    the c-sweep, well below typical case tolerances.
    """

    name: str = "fn_slab"
    method: str = "fn_method"
    geometry: str = "slab"
    n_modes: int = 10

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.fn_method.slab import (
            solve_fn_slab_bare_critical,
        )

        c = float(case.registry_case.materials[0].scattering_ratio[0])
        res = solve_fn_slab_bare_critical(c=c, n_modes=self.n_modes)
        return ScalarResult(
            tag="a_critical_mfp",
            value=float(res.a_critical_mfp),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                "determinant_residual": complex(res.determinant_residual),
                "nu0": float(res.nu0),
                "c": c,
            },
        )


@dataclass(frozen=True)
class FNSphereAdapter:
    r"""Adapter for :func:`...fn_method.sphere.solve_fn_sphere_bare_critical`.

    Reports ``R_critical_mfp``. Sphere F_N at ``n_modes = 10``
    reaches ~5e-8 absolute against Sood truth — exquisitely tight.
    """

    name: str = "fn_sphere"
    method: str = "fn_method"
    geometry: str = "sphere-1d"
    n_modes: int = 10

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.fn_method.sphere import (
            solve_fn_sphere_bare_critical,
        )

        c = float(case.registry_case.materials[0].scattering_ratio[0])
        res = solve_fn_sphere_bare_critical(c=c, n_modes=self.n_modes)
        return ScalarResult(
            tag="R_critical_mfp",
            value=float(res.R_critical_mfp),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                "determinant_residual": complex(res.determinant_residual),
                "c": c,
            },
        )


@dataclass(frozen=True)
class FNReflectedSlabAdapter:
    r"""Adapter for :func:`...fn_method.slab.solve_fn_slab_reflected_critical`.

    Reflected-slab F_N (Neshat-Maiorino 1980). Returns
    ``tau_critical_mfp`` — the core half-thickness at criticality
    given the reflector configuration.

    Each case must carry inline ``materials`` (mat_id 0 = core,
    mat_id 1 = reflector) and a ``structured_geometry`` with three
    regions in left-to-right order ``(reflector, core, reflector)``.
    The Case–Zweifel ``c`` for each material is read off
    :attr:`Mixture.scattering_ratio`; the reflector half-thickness in
    cm is read off ``structured_geometry.regions[0].outer_thickness_cm``
    and converted to mfp via the reflector mixture's :math:`\Sigma_t`
    (group 0).

    There is currently no trajectory_resolvent counterpart for
    reflected slab — this adapter has no agreement partner. That
    one-sided coverage is intentional; see
    ``.claude/scratch/cross_method_test_protocol_assessment.md``
    §"Out of scope".
    """

    name: str = "fn_reflected_slab"
    method: str = "fn_method"
    geometry: str = "reflected-slab"
    n_modes: int = 7

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.fn_method.slab import (
            solve_fn_slab_reflected_critical,
        )

        # Read c values off the inline materials and the reflector
        # half-thickness off structured_geometry.regions[0]. The expected
        # region layout is (reflector, core, reflector); we cross-
        # check that to fail loudly on accidental ordering bugs.
        if case.materials is None or 0 not in case.materials or 1 not in case.materials:
            raise ValueError(
                f"FNReflectedSlabAdapter: case {case.case_id!r} must "
                f"carry inline materials with mat_id=0 (core) and "
                f"mat_id=1 (reflector); got "
                f"{None if case.materials is None else sorted(case.materials)}"
            )
        geom = _structured_geometry_for(case)
        regions = geom.regions
        if len(regions) != 3:
            raise ValueError(
                f"FNReflectedSlabAdapter: case {case.case_id!r} must "
                f"carry a 3-region StructuredGeometry "
                f"(reflector, core, reflector); got "
                f"{len(regions)} regions"
            )
        if (regions[0].mat_id, regions[1].mat_id, regions[2].mat_id) != (1, 0, 1):
            raise ValueError(
                f"FNReflectedSlabAdapter: case {case.case_id!r} region "
                f"layout must be (reflector=1, core=0, reflector=1); "
                f"got mat_ids "
                f"({regions[0].mat_id}, {regions[1].mat_id}, "
                f"{regions[2].mat_id})"
            )
        c_core = float(case.materials[0].scattering_ratio[0])
        c_reflector = float(case.materials[1].scattering_ratio[0])
        # Convert the reflector cm thickness to mfp via the
        # reflector mixture's Σ_t. Under the unit-σ_t convention used
        # for Sood / NM 1980 reflected-slab cases this is the
        # identity, but we compute it explicitly to remove the
        # unit-equivalence assumption from the call site.
        sigma_t_reflector = float(case.materials[1].SigT[0])
        reflector_half_thickness_mfp = (
            float(regions[0].outer_thickness_cm) * sigma_t_reflector
        )

        res = solve_fn_slab_reflected_critical(
            c_core=c_core,
            c_reflector=c_reflector,
            reflector_half_thickness=reflector_half_thickness_mfp,
            n_modes=self.n_modes,
        )
        return ScalarResult(
            tag="tau_critical_mfp",
            value=float(res.tau_critical_mfp),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                "c_core": c_core,
                "c_reflector": c_reflector,
                "reflector_half_thickness_mfp": reflector_half_thickness_mfp,
                "converged": bool(res.converged),
            },
        )


# ═══════════════════════════════════════════════════════════════════
# trajectory_resolvent adapters (Variant α package)
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class TrajectoryResolventSlabAdapter:
    r"""Adapter for :func:`...trajectory_resolvent.greens_function_slab.solve_greens_function_slab`.

    Trajectory_resolvent solves the k-eigenvalue problem on the slab at a
    given full-width ``L``; the cross-method gate evaluates ``k_eff``
    at the **independently-known critical half-thickness from a
    different reference** (typically F_N's ``a_critical_mfp``). The
    adapter therefore reports ``k_eff`` (which should be 1.0 at the
    truth thickness).

    The continuous-albedo ``alpha`` is derived from
    ``structured_geometry.bcs[-1]`` (slab cases use symmetric BCs by
    convention) via :meth:`BC.to_alpha`; bare-critical slab registry
    cases are vacuum-on-vacuum (``α = 0``), closed slab is reflective-
    on-reflective (``α = 1``).

    Default quadrature: ``(n_x, n_mu, n_traj_quad) = (48, 128, 96)``.
    Slab vacuum has a near-cusp at μ=0 that needs ~128 angular nodes
    to resolve to ~1e-5.
    """

    name: str = "trajectory_resolvent_slab"
    method: str = "trajectory_resolvent"
    geometry: str = "slab"
    n_x: int = 48
    n_mu: int = 128
    n_traj_quad: int = 96
    max_iter: int = 500
    tol: float = 1e-9

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab import (
            solve_greens_function_slab,
        )

        sigma_t, sigma_s, nu_sigma_f = _extract_1g_xs(case)
        # Trajectory_resolvent slab takes FULL width L (not the half-
        # thickness). StructuredGeometry.domain_extent_cm IS the full
        # slab width in cm — no truth-vs-cm re-derivation needed.
        L_full_cm = _slab_L_full_cm(case)
        alpha = _outer_bc_for(case).to_alpha()

        res = solve_greens_function_slab(
            L=L_full_cm,
            sigma_t=sigma_t,
            sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha=alpha,
            n_x=self.n_x,
            n_mu=self.n_mu,
            n_traj_quad=self.n_traj_quad,
            max_iter=self.max_iter,
            tol=self.tol,
        )
        return ScalarResult(
            tag="k_eff",
            value=float(res.k_eff),
            solver_name=self.name,
            metadata={
                "n_x": self.n_x,
                "n_mu": self.n_mu,
                "n_traj_quad": self.n_traj_quad,
                "iterations": int(res.iterations),
                "converged": bool(res.converged),
                "L_full_cm": L_full_cm,
                "alpha": alpha,
            },
        )


@dataclass(frozen=True)
class TrajectoryResolventSphereAdapter:
    r"""Adapter for :func:`...trajectory_resolvent.greens_function.solve_greens_function_sphere`
    for bare-critical sphere cases.

    Reports ``k_eff`` at the **independently-known critical radius**
    (typically F_N's ``R_critical_mfp``). At ``α = 0`` and the truth
    radius, ``k_eff`` should be 1.0.

    The continuous-albedo ``alpha`` is derived from
    ``structured_geometry.bcs[-1]`` (the outer-surface BC) via
    :meth:`BC.to_alpha`. The inner BC at ``r = 0`` is the natural
    centreline reflective and is not parametrically relevant to the
    trajectory_resolvent operator.

    Closed-sphere ``α = 1`` cases (``k_eff = k_inf`` exactly) use
    :class:`TrajectoryResolventSphereClosedAdapter` instead — the
    parameter sets and convergence behaviour are different enough
    that two adapters keep the protocol clean.
    """

    name: str = "trajectory_resolvent_sphere"
    method: str = "trajectory_resolvent"
    geometry: str = "sphere-1d"
    n_r: int = 32
    n_mu: int = 32
    n_traj_quad: int = 64
    max_iter: int = 400
    tol: float = 1e-10

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
            solve_greens_function_sphere,
        )

        sigma_t, sigma_s, nu_sigma_f = _extract_1g_xs(case)
        # Read the radius in cm directly off the case's
        # StructuredGeometry. Sphere convention:
        # ``StructuredGeometry.domain_extent_cm`` IS R_cm.
        R_cm = _sphere_R_cm(case)
        alpha = _outer_bc_for(case).to_alpha()

        res = solve_greens_function_sphere(
            R=R_cm,
            sigma_t=sigma_t,
            sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha=alpha,
            n_r=self.n_r,
            n_mu=self.n_mu,
            n_traj_quad=self.n_traj_quad,
            max_iter=self.max_iter,
            tol=self.tol,
        )
        return ScalarResult(
            tag="k_eff",
            value=float(res.k_eff),
            solver_name=self.name,
            metadata={
                "n_r": self.n_r,
                "n_mu": self.n_mu,
                "n_traj_quad": self.n_traj_quad,
                "iterations": int(res.iterations),
                "converged": bool(res.converged),
                "R_cm": R_cm,
                "alpha": alpha,
            },
        )


@dataclass(frozen=True)
class TrajectoryResolventSphereClosedAdapter:
    r"""Adapter for closed-sphere (``α = 1``) k_inf cases.

    The closed sphere with perfect specular BC has rank-1 isotropic
    eigenmode and ``k_eff = k_inf = νΣ_f / Σ_a`` to machine
    precision (V_α1 algebraic identity). Useful as a multi-group
    cross-method gate where the bare-critical pillar is missing.

    Geometry, XS, and radius come from the case's inline
    ``materials`` + ``structured_geometry`` (the registry-less path).
    The continuous-albedo ``alpha`` is derived from
    ``structured_geometry.bcs[-1]`` via :meth:`BC.to_alpha`; closed
    sphere is :attr:`BC.reflective` on the outer surface, giving
    ``α = 1.0``.
    """

    name: str = "trajectory_resolvent_sphere_closed"
    method: str = "trajectory_resolvent"
    geometry: str = "closed-sphere-1d"
    n_r: int = 12
    n_mu: int = 12
    n_traj_quad: int = 24
    max_iter: int = 50
    tol: float = 1e-12

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
            solve_greens_function_sphere,
        )

        # Closed-sphere cases use the inline-materials path
        # (registry_case is None; materials + structured_geometry are set).
        sigma_t, sigma_s, nu_sigma_f = _extract_1g_xs_inline(case)
        R_cm = _sphere_R_cm(case)
        alpha = _outer_bc_for(case).to_alpha()
        res = solve_greens_function_sphere(
            R=R_cm,
            sigma_t=sigma_t,
            sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha=alpha,
            n_r=self.n_r,
            n_mu=self.n_mu,
            n_traj_quad=self.n_traj_quad,
            max_iter=self.max_iter,
            tol=self.tol,
        )
        return ScalarResult(
            tag="k_inf",
            value=float(res.k_eff),
            solver_name=self.name,
            metadata={
                "n_r": self.n_r,
                "n_mu": self.n_mu,
                "n_traj_quad": self.n_traj_quad,
                "iterations": int(res.iterations),
                "converged": bool(res.converged),
                "sigma_t": sigma_t,
                "sigma_s": sigma_s,
                "nu_sigma_f": nu_sigma_f,
                "R_cm": R_cm,
                "alpha": alpha,
            },
        )


# ═══════════════════════════════════════════════════════════════════
# Helpers — XS / parameter extraction from CrossMethodCase
# ═══════════════════════════════════════════════════════════════════


def _extract_1g_xs(case: CrossMethodCase) -> tuple[float, float, float]:
    r"""Extract :math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f)` for a 1G case
    from a registry-backed case.

    Pulls from ``case.registry_case.materials[0]`` via
    :func:`mixture_to_fn_arrays`. Raises if the case is multi-group
    (1G adapters can't consume those) or if the case carries no
    registry case (use :func:`_extract_1g_xs_inline` for that path).
    """
    if case.registry_case is None:
        raise ValueError(
            f"CrossMethodCase {case.case_id!r} has registry_case=None; "
            f"the registry-backed XS extractor cannot serve this case. "
            f"Use _extract_1g_xs_inline for inline-materials cases."
        )
    return _xs_from_materials_dict(
        case.registry_case.materials, case.case_id
    )


def _extract_1g_xs_inline(case: CrossMethodCase) -> tuple[float, float, float]:
    r"""Extract :math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f)` for a 1G case
    from inline ``case.materials``.

    Used by adapters whose case carries inline materials + geometry_spec
    (the no-registry path — closed-sphere k_inf, MMS, custom
    configurations).
    """
    if case.materials is None:
        raise ValueError(
            f"CrossMethodCase {case.case_id!r} has materials=None; "
            f"_extract_1g_xs_inline requires inline materials. Use "
            f"_extract_1g_xs for registry-backed cases."
        )
    return _xs_from_materials_dict(case.materials, case.case_id)


def _xs_from_materials_dict(
    materials: dict, case_id: str,
) -> tuple[float, float, float]:
    """Common backend: pull 1G ``(σ_t, σ_s, νσ_f)`` from a materials dict.

    Reads directly off ``Mixture.SigT`` / ``SigS[0]`` / ``SigP``
    (the production-protocol surface).
    """
    primary = materials[0]
    sigma_t_arr = np.asarray(primary.SigT, dtype=float)
    sigma_s_arr = primary.SigS[0].toarray().astype(float)
    nu_sigma_f_arr = np.asarray(primary.SigP, dtype=float)
    if sigma_t_arr.shape[0] != 1:
        raise ValueError(
            f"_xs_from_materials_dict: case {case_id!r} is "
            f"{sigma_t_arr.shape[0]}G; expected 1G"
        )
    return (
        float(sigma_t_arr[0]),
        float(sigma_s_arr[0, 0]),
        float(nu_sigma_f_arr[0]),
    )


def _structured_geometry_for(case: CrossMethodCase):
    r"""Resolve the :class:`StructuredGeometry` for a case.

    Reads from ``case.structured_geometry`` (inline / override path)
    or builds it via ``case.registry_case.to_geometry()`` (registry
    path), whichever is populated. The override path takes precedence
    when both are present (cross-method agreement gates substitute a
    predicted critical dimension by setting an inline structured
    geometry).
    """
    if case.structured_geometry is not None:
        return case.structured_geometry
    if case.registry_case is not None and hasattr(
        case.registry_case, "to_geometry"
    ):
        return case.registry_case.to_geometry()
    raise ValueError(
        f"_structured_geometry_for: case {case.case_id!r} has neither "
        f"inline structured_geometry nor a registry_case carrying one."
    )


def _sphere_R_cm(case: CrossMethodCase) -> float:
    r"""Return the sphere radius in cm from the case's StructuredGeometry.

    For SPH geometry, ``domain_extent_cm`` IS R_cm (single-region
    radius, sum trivially equals the radius).
    """
    geom = _structured_geometry_for(case)
    if geom.geometry != "SPH":
        raise ValueError(
            f"_sphere_R_cm: case {case.case_id!r} structured geometry "
            f"is {geom.geometry!r}, expected 'SPH'"
        )
    return float(geom.domain_extent_cm)


def _slab_L_full_cm(case: CrossMethodCase) -> float:
    r"""Return the slab full width in cm from the case's StructuredGeometry.

    Slab convention: ``StructuredGeometry.domain_extent_cm`` IS the
    FULL slab width :math:`[0, L]`, which is exactly what
    :func:`solve_greens_function_slab` expects as its ``L`` argument.
    """
    geom = _structured_geometry_for(case)
    if geom.geometry != "SLB":
        raise ValueError(
            f"_slab_L_full_cm: case {case.case_id!r} structured "
            f"geometry is {geom.geometry!r}, expected 'SLB'"
        )
    return float(geom.domain_extent_cm)


def _outer_bc_for(case: CrossMethodCase):
    """Return the outer-surface BC for a case.

    For ``SLB`` geometry the slab is symmetric vacuum-vacuum (or
    closed reflective-reflective); both endpoints share the same kind
    in the cases this protocol covers, so we return ``bcs[1]``
    (the right end). For ``SPH`` / ``CYL`` geometry the single endpoint
    in :attr:`StructuredGeometry.bcs` IS the outer-surface BC.
    """
    geom = _structured_geometry_for(case)
    return geom.bcs[-1]


# ═══════════════════════════════════════════════════════════════════
# Adapter registry — used by tests and (future) agreement-matrix renderer
# ═══════════════════════════════════════════════════════════════════


ADAPTERS_BY_NAME: dict[str, object] = {
    "fn_slab": FNSlabAdapter(),
    "fn_sphere": FNSphereAdapter(),
    "fn_reflected_slab": FNReflectedSlabAdapter(),
    "trajectory_resolvent_slab": TrajectoryResolventSlabAdapter(),
    "trajectory_resolvent_sphere": TrajectoryResolventSphereAdapter(),
    "trajectory_resolvent_sphere_closed": TrajectoryResolventSphereClosedAdapter(),
}
"""All registered adapters. New adapters MUST register here so the
agreement-matrix renderer (and future cross-method audit tools)
can discover them.
"""

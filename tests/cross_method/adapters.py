r"""Solver adapters for the cross-method test protocol.

Each adapter wraps a continuous-reference solver to the
:class:`~tests.cross_method.protocol.SolverAdapter` shape. The
adapter:

* extracts the right cross sections from
  ``case.registry_case.materials`` via
  :func:`...sood_registry.extractors.mixture_to_fn_arrays`;
* selects internal numerical parameters (n_modes for fn_method,
  n_r/n_mu/n_traj_quad for trajectory_resolvent) based on the
  requested case tolerance — looser tolerance → faster solve;
* performs unit conversions (mfp ↔ cm, half-thickness ↔ full
  slab);
* returns a :class:`ScalarResult` with the right ``tag``.

Adding a new solver
-------------------

Subclass the appropriate base or write a fresh dataclass with
``name``, ``method``, ``geometry`` fields and a ``solve`` method
matching the :class:`SolverAdapter` Protocol. Register the new
adapter in :data:`ADAPTERS_BY_NAME` so the agreement-matrix
renderer can find it. The adapter classes are intentionally
plain — no inheritance — so the protocol stays Protocol-typed
rather than ABC-typed.

References
----------

* :doc:`/skills/algebra-of-record` — structural-independence
  ladder; adapters trust ``numpy``/``scipy`` (above the trusted-
  library line) but not in-house primitives.
* :doc:`/skills/vv-principles` — quadrature-floor anti-patterns
  (don't pick tolerances tighter than the floor).
"""
from __future__ import annotations

from dataclasses import dataclass

from orpheus.derivations.continuous.sood_registry.extractors import (
    mixture_to_fn_arrays,
)

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

        c = _extract_c(case)
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

        c = _extract_c(case)
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
    given the reflector configuration. Each case must carry
    ``reflector_half_thickness_mfp`` and ``c_reflector`` in its
    metadata; this adapter reads them off the case's notes /
    attributes.

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

        # Reflected-slab cases carry their non-XS parameters in
        # case.registry_case attribute is None; the case is fully
        # described inline. Pull from a documented attribute path.
        params = _reflected_slab_params(case)
        res = solve_fn_slab_reflected_critical(
            c_core=params["c_core"],
            c_reflector=params["c_reflector"],
            reflector_half_thickness=params["reflector_half_thickness_mfp"],
            n_modes=self.n_modes,
        )
        return ScalarResult(
            tag="tau_critical_mfp",
            value=float(res.tau_critical_mfp),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                **params,
                "converged": bool(res.converged),
            },
        )


# ═══════════════════════════════════════════════════════════════════
# trajectory_resolvent adapters (Variant α package)
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class TrajectoryResolventSlabAdapter:
    r"""Adapter for :func:`...trajectory_resolvent.greens_function_slab.solve_greens_function_slab`
    at vacuum BC (``alpha = 0``).

    Trajectory_resolvent solves the k-eigenvalue problem on the slab at a
    given full-width ``L``; the cross-method gate evaluates ``k_eff``
    at the **independently-known critical half-thickness from a
    different reference** (typically F_N's ``a_critical_mfp``). The
    adapter therefore reports ``k_eff`` (which should be 1.0 at the
    truth thickness).

    To use this adapter the case MUST be parametrised on
    ``a_critical_mfp`` truth (or carry the F_N predicted thickness
    as a separate input). This is the canonical pattern in
    ``test_fn_la13511_slab_xverif.py``.

    Default quadrature: ``(n_x, n_mu, n_traj_quad) = (48, 128, 96)``.
    Slab vacuum has a near-cusp at μ=0 that needs ~128 angular nodes
    to resolve to ~1e-5.
    """

    name: str = "trajectory_resolvent_slab"
    method: str = "trajectory_resolvent"
    geometry: str = "slab"
    alpha: float = 0.0
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
        # Trajectory_resolvent slab takes FULL width L (not the half-thickness).
        # The caller supplies a_truth (mfp) on the case via
        # case.notes / case.tolerances mapping; resolve via helper.
        a_truth_mfp = _slab_a_truth_mfp(case)
        L_full_cm = 2.0 * a_truth_mfp / sigma_t

        res = solve_greens_function_slab(
            L=L_full_cm,
            sigma_t=sigma_t,
            sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha=self.alpha,
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
                "a_truth_mfp": a_truth_mfp,
            },
        )


@dataclass(frozen=True)
class TrajectoryResolventSphereAdapter:
    r"""Adapter for :func:`...trajectory_resolvent.greens_function.solve_greens_function_sphere`
    at vacuum BC (``alpha = 0``) for bare-critical sphere cases.

    Reports ``k_eff`` at the **independently-known critical radius**
    (typically F_N's ``R_critical_mfp``). At α=0 and the truth
    radius, ``k_eff`` should be 1.0.

    Closed-sphere α=1 cases (``k_eff = k_inf`` exactly) use
    :class:`TrajectoryResolventSphereClosedAdapter` instead — the
    parameter sets and convergence behaviour are different enough
    that two adapters keep the protocol clean.
    """

    name: str = "trajectory_resolvent_sphere"
    method: str = "trajectory_resolvent"
    geometry: str = "sphere-1d"
    alpha: float = 0.0
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
        R_truth_mfp = _sphere_R_truth_mfp(case)
        R_cm = R_truth_mfp / sigma_t

        res = solve_greens_function_sphere(
            R=R_cm,
            sigma_t=sigma_t,
            sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha=self.alpha,
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
                "alpha": self.alpha,
            },
        )


@dataclass(frozen=True)
class TrajectoryResolventSphereClosedAdapter:
    r"""Adapter for closed-sphere (α=1) k_inf cases.

    The closed sphere with perfect specular BC has rank-1 isotropic
    eigenmode and ``k_eff = k_inf = νΣ_f / Σ_a`` to machine
    precision (V_α1 algebraic identity). Useful as a multi-group
    cross-method gate where the bare-critical pillar is missing.

    Geometry XS and radius come from the case ``notes`` field
    (``sigma_t=A.AA sigma_s=B.BB nu_sigma_f=C.CC R_cm=D.DD``).
    Defaults to fuel-A-like fixture.
    """

    name: str = "trajectory_resolvent_sphere_closed"
    method: str = "trajectory_resolvent"
    geometry: str = "closed-sphere-1d"
    alpha: float = 1.0
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
        # (registry_case is None; materials + mesh_template are set).
        sigma_t, sigma_s, nu_sigma_f = _extract_1g_xs_inline(case)
        R_cm = case.mesh_template.critical_dimension_cm  # type: ignore[union-attr]
        if R_cm is None:
            raise ValueError(
                f"Closed-sphere case {case.case_id!r} has no "
                f"critical_dimension_cm in its mesh_template"
            )
        res = solve_greens_function_sphere(
            R=float(R_cm),
            sigma_t=sigma_t,
            sigma_s=sigma_s,
            nu_sigma_f=nu_sigma_f,
            alpha=self.alpha,
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
                "R_cm": float(R_cm),
                "alpha": self.alpha,
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

    Used by adapters whose case carries inline materials + mesh_template
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
    """Common backend: pull 1G ``(σ_t, σ_s, νσ_f)`` from a materials dict."""
    primary = materials[0]
    sigma_t_arr, sigma_s_arr, nu_sigma_f_arr, _chi = mixture_to_fn_arrays(
        primary
    )
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


def _extract_c(case: CrossMethodCase) -> float:
    r"""Compute :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t` for a 1G case.

    F_N method takes ``c`` as its primary input parameter; this is
    the canonical extraction for 1G isotropic-scattering cases.
    """
    sigma_t, sigma_s, nu_sigma_f = _extract_1g_xs(case)
    return (sigma_s + nu_sigma_f) / sigma_t


def _slab_a_truth_mfp(case: CrossMethodCase) -> float:
    r"""Return the case's truth half-thickness in mean-free paths.

    For bare-critical slab cases the truth tag must be
    ``"a_critical_mfp"`` and the value is the half-thickness in
    mfp.
    """
    if case.truth_tag != "a_critical_mfp":
        raise ValueError(
            f"_slab_a_truth_mfp: case {case.case_id!r} truth_tag is "
            f"{case.truth_tag!r}, expected 'a_critical_mfp'"
        )
    return float(case.truth_value)


def _sphere_R_truth_mfp(case: CrossMethodCase) -> float:
    r"""Return the case's truth radius in mean-free paths."""
    if case.truth_tag != "R_critical_mfp":
        raise ValueError(
            f"_sphere_R_truth_mfp: case {case.case_id!r} truth_tag is "
            f"{case.truth_tag!r}, expected 'R_critical_mfp'"
        )
    return float(case.truth_value)


def _parse_notes_kv(case: CrossMethodCase) -> dict[str, str]:
    """Parse ``key=value`` tokens from ``case.notes``.

    Convention: notes contains space-separated tokens of the form
    ``"key=value"``. Tokens without ``=`` are skipped. Used by
    inline-parametrised adapters (reflected slab, closed sphere)
    to read their parameters off the case without needing a
    registry entry.
    """
    return dict(
        token.split("=", 1)
        for token in case.notes.split()
        if "=" in token
    )


def _reflected_slab_params(case: CrossMethodCase) -> dict[str, float]:
    r"""Extract reflected-slab inline parameters from
    ``CrossMethodCase.notes``.

    NM 1980 / Sood Table 10 cases don't have ORPHEUS registry
    entries — they're parametrised by ``(c_core, c_reflector,
    reflector_half_thickness_mfp)`` directly. The required tokens
    are listed in case.notes as ``c_core=X.XX c_reflector=Y.YY
    reflector_half_thickness_mfp=Z.ZZ``.
    """
    tokens = _parse_notes_kv(case)
    required = {
        "c_core",
        "c_reflector",
        "reflector_half_thickness_mfp",
    }
    missing = required - set(tokens)
    if missing:
        raise ValueError(
            f"_reflected_slab_params: case {case.case_id!r} missing "
            f"keys in notes: {missing}. Got: {dict(tokens)}"
        )
    return {k: float(tokens[k]) for k in required}


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

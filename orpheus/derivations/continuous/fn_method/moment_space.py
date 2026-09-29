r"""The F_N moment space — math-heart class for the F_N method.

Phase D
-------

``MomentSpace`` is a **reference solution generator**. It consumes a
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
directly via its ``__init__`` (no mesh, no cell counts). Infinite-
medium :math:`k_\infty` is a separate static method
:meth:`MomentSpace.solve_kinf` that takes a :class:`Mixture` only —
no geometry needed. The architectural split between continuous
references and discrete-mesh production solvers is now full.

The single class :class:`MomentSpace` encapsulates everything the
F_N method needs to commit to before computing a solution:

* **Geometry** — slab vs sphere (cylinder is intentionally out of
  pillar; see :mod:`...singular_eigenfunction.cylinder`),
* **Materials** — the cross sections specifying :math:`c =
  (\Sigma_s + \nu\Sigma_f)/\Sigma_t`,
* **Boundary conditions** — bare-critical (vacuum) is the
  shipping configuration; reflected slabs use the parallel
  :class:`...slab.reflected.SlabReflectedFNResult` path,
* **F_N order** :math:`N` — the truncation of the moment basis,
* **Flux-reconstruction strategy** — none, Atkinson product-Nyström
  (the recommended ERR-036-fixed path), or the legacy plain
  Gauss–Legendre path (kept for diagnostic comparison).

The class is intentionally a thin computational specialisation of
:class:`StructuredGeometry` for the F_N method — it is the *F_N
moment space*, not the *F_N solver*. The solve methods are thin
wrappers around the existing function-level API in
:mod:`...slab.one_group`, :mod:`...sphere.one_group`,
:mod:`...slab.flux_reconstruction`, and
:mod:`...sphere.flux_reconstruction`.

Slab convention reminder
------------------------

:attr:`StructuredGeometry.domain_extent_cm` is the FULL slab width
(production-natural convention). F_N's natural half-thickness
:math:`a = L / 2` is computed inside
:meth:`MomentSpace.solve_critical` from the geometry's full extent.

References
----------

* Siewert, Benoist (1979) "The F_N Method in Neutron-Transport
  Theory. Part I: Theory and Applications", *Nuclear Science and
  Engineering* **69**, 156–169.
* Grandjean, Siewert (1979) "The F_N Method in Neutron-Transport
  Theory. Part II: Applications and Numerical Results", *Nuclear
  Science and Engineering* **69**, 161–168.
* Siewert, Thomas (1986) "Neutron Transport Calculations in Bare
  and Reflected Spheres of Multiplying Material", *Nuclear Science
  and Engineering* **94**, 264–270.
* Kaper, Lindeman, Leaf (1974) "Benchmark Values for the Slab
  and Sphere Criticality Problem in One-Group Neutron Transport
  Theory", *Nuclear Science and Engineering* **54**, 94–98.
* Atkinson, K.E. (1997) *The Numerical Solution of Integral
  Equations of the Second Kind*. Cambridge University Press.
* :doc:`/theory/references/fn_method` — the canonical theory exposition.
"""
from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.solution_types import (
    CriticalSolution,
    FluxSolution,
)
from orpheus.derivations.common.reference_body import (
    HomogeneousBody,
    ReflectedSlab,
    describe,
    reference_body,
    refuse_unserved,
    require_vacuum,
)
from orpheus.geometry import CoordSystem, StructuredGeometry


FluxReconstructionStrategy = Literal["none", "atkinson_nystrom", "legacy_gl"]


def _require_one_group_equal_total(core: Mixture, reflector: Mixture) -> None:
    """The reflected-slab F_N of Neshat and Maiorino (1980) is one-group and
    puts both media on one mean-free-path scale, so the core and the
    reflector must be one-group with one total cross section (Sood's
    two-media slabs are stated this way: Tables 2, 9 and 13)."""
    sigma_core = np.asarray(core.SigT, dtype=float)
    sigma_reflector = np.asarray(reflector.SigT, dtype=float)
    if sigma_core.shape != (1,) or sigma_reflector.shape != (1,):
        refuse_unserved(
            f"reflected slab of {sigma_core.shape[0]} and "
            f"{sigma_reflector.shape[0]} groups", owner="MomentSpace",
            missing="a multi-group reflected-slab F_N; only the one-group "
            "reflected slab (Neshat and Maiorino 1980) is built",
        )
    if not math.isclose(sigma_core[0], sigma_reflector[0], rel_tol=1e-12, abs_tol=0.0):
        refuse_unserved(
            f"reflected slab whose core and reflector total cross sections "
            f"differ ({sigma_core[0]!r}, {sigma_reflector[0]!r})",
            owner="MomentSpace",
            missing="a reflected-slab F_N for two mean-free-path scales; "
            "Neshat and Maiorino's needs the core and the reflector to "
            "share one total cross section",
        )


@dataclass(frozen=True)
class MomentSpace:
    r"""The F_N moment space — Galerkin half-range projection of the angular flux.

    The F_N method (Siewert–Benoist 1979) attacks the one-speed
    transport equation by **Galerkin projection of the angular flux
    :math:`\psi(r, \mu)` onto a finite-dimensional moment basis**, then
    collocation on a discrete set of :math:`\mu` values to close the
    system. ``MomentSpace`` is the data class that owns the basis
    (the F_N order :math:`N`), the geometry (slab vs sphere — the two
    differ by exactly one sign in the boundary attenuation term), the
    materials (encoded as :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`
    for the 1G isotropic case), and the boundary condition
    (bare-critical vacuum is the shipping configuration).

    Construction (Phase D)
    ----------------------

    Direct construction with a :class:`StructuredGeometry` (slab or
    sphere) and a ``materials: dict[int, Mixture]`` payload::

        from orpheus.geometry import BC, CoordSystem, StructuredGeometry

        geom = StructuredGeometry(
            coord=CoordSystem.CARTESIAN,
            breakpoints=(0.0, L_full),
            mat_ids=(0,),
            boundaries=(BC.vacuum, BC.vacuum),
        )
        ms = MomentSpace(geometry=geom, materials={0: mix})
        sol = ms.solve_critical()

    For infinite-medium :math:`k_\infty`, use
    :meth:`MomentSpace.solve_kinf` — no geometry needed::

        kinf = MomentSpace.solve_kinf(mix)

    Parameters
    ----------
    geometry : :class:`StructuredGeometry`
        Pure-geometry layer object: one material filling a slab
        (``CoordSystem.CARTESIAN``) or a solid sphere
        (``CoordSystem.SPHERICAL``), read by
        :func:`~orpheus.derivations.common.reference_body.reference_body`, or a
        symmetric reflected slab of one group whose core and reflector share
        the total cross section (Neshat and Maiorino 1980).
        A cylinder is out of pillar (Westfall–Metcalf 1972 — see
        :mod:`...singular_eigenfunction.cylinder`).
    materials : dict[int, Mixture]
        Production-protocol materials, keyed by material ID.
    fn_order : int
        F_N truncation order :math:`N`. Defaults to 9.
    flux_reconstruction : :data:`FluxReconstructionStrategy`
        Which flux-reconstruction path to use. Default
        ``"atkinson_nystrom"``.
    """

    geometry: StructuredGeometry
    materials: dict[int, Mixture]
    fn_order: int = 9
    flux_reconstruction: FluxReconstructionStrategy = "atkinson_nystrom"

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        """Validate the F_N method's structural preconditions.

        The F_N method as shipped applies to **slab** (Siewert–Benoist
        1979) and **sphere** (Siewert–Thomas 1986). It is **explicitly
        out of pillar** for cylinder geometry — Westfall–Metcalf 1972
        documents that the Mitsis-style Wiener–Hopf reduction is
        non-convergent for the bare cylinder. Cylinder critical
        dimensions ship via :mod:`...singular_eigenfunction.cylinder`.
        """
        body = self._body
        if body.coord is CoordSystem.CYLINDRICAL:
            refuse_unserved(
                describe(body), owner="MomentSpace",
                missing="an F_N cylinder, which is out of pillar "
                "(Westfall-Metcalf 1972; see singular_eigenfunction.cylinder, "
                "and #170). For infinite-medium k_inf, use "
                "MomentSpace.solve_kinf(mixture), with no geometry",
            )
        require_vacuum(
            self.geometry, owner="MomentSpace",
            missing="the F_N solvers here are bare or core-reflected slabs and "
            "spheres with vacuum outer faces",
        )
        if self.fn_order < 0:
            raise ValueError(f"fn_order must be ≥ 0, got {self.fn_order}")
        for mat_id in self._body_mat_ids:
            if mat_id not in self.materials:
                raise ValueError(
                    f"materials dict missing mat_id={mat_id} "
                    f"required by geometry; got keys "
                    f"{sorted(self.materials.keys())}"
                )
        if isinstance(body, ReflectedSlab):
            _require_one_group_equal_total(
                self.materials[body.core_mat_id], self.materials[body.reflector_mat_id],
            )
        if self.flux_reconstruction not in {"none", "atkinson_nystrom", "legacy_gl"}:
            raise ValueError(
                f"flux_reconstruction must be one of "
                f"{{none, atkinson_nystrom, legacy_gl}}, got "
                f"{self.flux_reconstruction!r}"
            )

    # ------------------------------------------------------------------
    # Material accessors — read directly from Mixture (no extractor)
    # ------------------------------------------------------------------

    @functools.cached_property
    def _body(self) -> HomogeneousBody | ReflectedSlab:
        """The body this generator solves on: a homogeneous slab or sphere, or
        a symmetric reflected slab; any other shape is refused
        (:func:`~orpheus.derivations.common.reference_body.refuse_unserved`)."""
        match reference_body(self.geometry):
            case HomogeneousBody() | ReflectedSlab() as body:
                return body
            case other:
                refuse_unserved(
                    describe(other), owner="MomentSpace",
                    missing="an F_N solver for a layered, hollow or asymmetric "
                    "configuration (only the symmetric reflected slab of "
                    "Neshat and Maiorino 1980 is built)",
                )

    @property
    def _body_mat_ids(self) -> tuple[int, ...]:
        """Every material the body holds."""
        match self._body:
            case HomogeneousBody(mat_id=mat_id):
                return (mat_id,)
            case ReflectedSlab(core_mat_id=core, reflector_mat_id=reflector):
                return (core, reflector)

    @property
    def _mat_id(self) -> int:
        """The multiplying material: the body's, or a reflected slab's core."""
        match self._body:
            case HomogeneousBody(mat_id=mat_id):
                return mat_id
            case ReflectedSlab(core_mat_id=core):
                return core

    @property
    def _mixture(self) -> Mixture:
        """The active :class:`Mixture` for this geometry's mat_id."""
        return self.materials[self._mat_id]

    @property
    def c(self) -> float:
        r"""Mean number of secondaries per collision, :math:`c`, of the
        multiplying material (the body's; a reflected slab's core).

        For 1G isotropic-scattering problems
        :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`.

        Raises
        ------
        ValueError
            If the active mixture has more than one energy group.
        """
        mixture = self._mixture
        sig_t = np.asarray(mixture.SigT, dtype=float)
        if sig_t.shape[0] != 1:
            raise ValueError(
                f"MomentSpace.c requires a 1G mixture (got "
                f"{sig_t.shape[0]}G). Multi-group problems should call "
                f"solve_kinf(mixture) directly."
            )
        return float(mixture.scattering_ratio[0])

    @property
    def n_groups(self) -> int:
        """Number of energy groups of the multiplying material (the body's; a
        reflected slab's core, which shares it with the reflector)."""
        return int(np.asarray(self._mixture.SigT).shape[0])

    # ------------------------------------------------------------------
    # solve_critical — the canonical critical-configuration call
    # ------------------------------------------------------------------

    def solve_critical(
        self,
        *,
        n_bracket: Optional[int] = None,
        bisect_tol: float = 1e-12,
        max_bisect: int = 80,
    ) -> CriticalSolution:
        r"""Solve the critical configuration for the F_N moment space.

        For a **bare-critical slab** this returns the critical
        half-thickness :math:`a` in mean free paths. F_N's natural
        unit is the half-thickness, computed internally from the
        geometry's full extent as :math:`a = L_{\rm cm} / 2 \cdot \Sigma_t`.

        For a **bare-critical sphere** this returns the critical
        radius :math:`R_c` in mean free paths, computed as
        :math:`R_c = R_{\rm cm} \cdot \Sigma_t`.

        Parameters
        ----------
        n_bracket : int | None, default None
            Initial bracket-scan resolution.
        bisect_tol : float, default 1e-12
            Bisection tolerance on the configuration parameter.
        max_bisect : int, default 80
            Maximum bisection iterations.

        Returns
        -------
        :class:`CriticalSolution`
        """
        if isinstance(self._body, ReflectedSlab):
            if (n_bracket, bisect_tol, max_bisect) != (None, 1e-12, 80):
                raise ValueError(
                    "MomentSpace.solve_critical: n_bracket, bisect_tol and "
                    "max_bisect control the bare slab's and sphere's bisection; "
                    "the reflected slab's F_N solver (Neshat and Maiorino 1980) "
                    "takes none of them."
                )
            return self._solve_critical_reflected_slab(self._body)
        coord = self.geometry.coord
        mixture = self._mixture
        sig_t = np.asarray(mixture.SigT, dtype=float)
        n_groups = sig_t.shape[0]

        if n_groups != 1:
            raise NotImplementedError(
                "MomentSpace.solve_critical for slab/sphere is currently "
                "1G-only — the multi-group F_N spatial extension "
                "(Siewert-Thomas 1986 2G machinery) is not yet wired "
                "through this class facade. For multi-group "
                "infinite-medium k_inf, use "
                "MomentSpace.solve_kinf(mixture)."
            )

        c = float(mixture.scattering_ratio[0])
        if c <= 1.0:
            raise ValueError(
                f"F_N bare-critical {coord.name.lower()} requires c > 1 "
                f"(multiplying medium); got c={c} from mixture "
                f"(sigma_s + nu_sigma_f) / sigma_t, sigma_t = {sig_t[0]}."
            )

        match coord:
            case CoordSystem.CARTESIAN:
                return self._solve_critical_slab(c, n_bracket, bisect_tol, max_bisect)
            case CoordSystem.SPHERICAL:
                return self._solve_critical_sphere(c, n_bracket, bisect_tol, max_bisect)

        raise NotImplementedError(  # pragma: no cover (validated above)
            f"MomentSpace.solve_critical: unhandled coordinate system {coord!r}"
        )

    def _solve_critical_reflected_slab(self, body: ReflectedSlab) -> CriticalSolution:
        r"""The critical core half-thickness of a symmetric reflected slab.

        Neshat and Maiorino (1980): given the core's and the reflector's
        :math:`c` and the reflector thickness :math:`\Delta` in mean free
        paths, the F_N method returns the core half-thickness :math:`\tau`
        at which the slab is critical. As for the bare slab, the geometry's
        core width is the question's unknown, not an input: only the
        reflector's width is read. Both media share one :math:`\Sigma_t`
        (checked at construction), which converts the width to mean free
        paths.
        """
        from .slab.reflected import solve_fn_slab_reflected_critical

        core = self.materials[body.core_mat_id]
        reflector = self.materials[body.reflector_mat_id]
        sigma_t = float(np.asarray(core.SigT, dtype=float)[0])
        res = solve_fn_slab_reflected_critical(
            c_core=float(core.scattering_ratio[0]),
            c_reflector=float(reflector.scattering_ratio[0]),
            reflector_half_thickness=body.reflector_width_cm * sigma_t,
            n_modes=self.fn_order,
        )
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.tau_critical_mfp),
            parameter_kind="core_half_thickness_mfp",
            converged=bool(res.converged),
            metadata={
                "n_groups": 1,
                "method": "solve_fn_slab_reflected_critical",
                "n_modes": self.fn_order,
                "c_core": float(res.c_core),
                "c_reflector": float(res.c_reflector),
                "reflector_half_thickness_mfp": body.reflector_width_cm * sigma_t,
                "raw_result": res,
            },
        )

    def _solve_critical_slab(
        self,
        c: float,
        n_bracket: Optional[int],
        bisect_tol: float,
        max_bisect: int,
    ) -> CriticalSolution:
        from .slab.one_group import solve_fn_slab_bare_critical

        # Pass-through kwargs only when explicitly set; otherwise the
        # underlying solver's defaults take effect (preserves bit-equality
        # for callers that don't pass overrides).
        kwargs: dict = {"c": c, "n_modes": self.fn_order,
                        "bisect_tol": bisect_tol, "max_bisect": max_bisect}
        if n_bracket is not None:
            kwargs["n_bracket"] = n_bracket
        res = solve_fn_slab_bare_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.a_critical_mfp),
            parameter_kind="half_thickness_mfp",
            # READ from the bisection (#340).  This said ``True`` with the
            # comment "bisection always converges given bracket" until
            # 2026-08-09 — present-tense false, since ``max_bisect`` is
            # caller-tunable and forwarded right above.
            converged=bool(res.converged),
            metadata={
                "n_groups": 1,
                "method": "solve_fn_slab_bare_critical",
                "n_modes": int(res.N),
                "c": float(res.c),
                "nu0": float(res.nu0),
                "determinant_residual": complex(res.determinant_residual),
                "n_bracket_points": int(res.n_bracket_points),
                "raw_result": res,
            },
        )

    def _solve_critical_sphere(
        self,
        c: float,
        n_bracket: Optional[int],
        bisect_tol: float,
        max_bisect: int,
    ) -> CriticalSolution:
        from .sphere.one_group import solve_fn_sphere_bare_critical

        kwargs: dict = {"c": c, "n_modes": self.fn_order,
                        "bisect_tol": bisect_tol, "max_bisect": max_bisect}
        if n_bracket is not None:
            kwargs["n_bracket"] = n_bracket
        res = solve_fn_sphere_bare_critical(**kwargs)
        return CriticalSolution(
            eigenvalue=1.0,
            eigenvalue_kind="k_eff",
            parameter_value=float(res.R_critical_mfp),
            parameter_kind="radius_mfp",
            # scipy's own OptimizeResult.success, threaded through (#340).
            # It was read and branched on inside the solver, then dropped.
            converged=bool(res.converged),
            metadata={
                "n_groups": 1,
                "method": "solve_fn_sphere_bare_critical",
                "n_modes": int(res.N),
                "c": float(res.c),
                "nu0": float(res.nu0),
                "determinant_residual": complex(res.determinant_residual),
                "n_bracket_points": int(res.n_bracket_points),
                "raw_result": res,
            },
        )

    # ------------------------------------------------------------------
    # solve_kinf — infinite-medium k_inf (no geometry needed)
    # ------------------------------------------------------------------

    @staticmethod
    def solve_kinf(mixture: Mixture) -> CriticalSolution:
        r"""Compute infinite-medium :math:`k_\infty` for *mixture*.

        Replaces the pre-Phase-D ``MomentSpace.from_problem`` path
        (which used the now-retired ``GeometrySpec`` carrier with
        ``geometry="infinite"``) with a direct mixture-only entry
        point. No geometry needed — :math:`k_\infty` is a material
        property of the medium.

        Dispatches to :func:`...multi_group.k_inf.compute_kinf_*`
        based on the mixture's group count (1G, 2G, multi-group).

        Parameters
        ----------
        mixture : :class:`Mixture`
            Production-protocol mixture. Must have :math:`\Sigma_t`,
            scattering matrices :math:`\Sigma_s`, fission production
            :math:`\nu\Sigma_f`, fission spectrum :math:`\chi`.

        Returns
        -------
        :class:`CriticalSolution`
            With ``eigenvalue_kind="k_inf"``,
            ``parameter_kind="k_inf_only"``, ``parameter_value=0.0``,
            and method-specific diagnostics in :attr:`metadata`.
        """
        from .multi_group.k_inf import (
            compute_kinf_1g,
            compute_kinf_2g_general,
            compute_kinf_mg,
        )

        sig_t = np.asarray(mixture.SigT, dtype=float)
        sig_s_p0 = mixture.SigS[0].toarray().astype(float)
        nu_sig_f = np.asarray(mixture.SigP, dtype=float)
        chi = np.asarray(mixture.chi, dtype=float)
        n_groups = sig_t.shape[0]

        if n_groups == 1:
            k_inf = compute_kinf_1g(
                float(sig_t[0]),
                float(sig_s_p0[0, 0]),
                float(nu_sig_f[0]),
            )
            method = "compute_kinf_1g"
        elif n_groups == 2:
            # Use the general (Eq 28) formula — it reduces to Eq 29 when
            # there's no upscatter, so it's the single source of truth.
            k_inf = compute_kinf_2g_general(sig_t, sig_s_p0, nu_sig_f, chi)
            method = "compute_kinf_2g_general"
        else:
            k_inf = compute_kinf_mg(sig_t, sig_s_p0, nu_sig_f, chi)
            method = "compute_kinf_mg"

        return CriticalSolution(
            eigenvalue=k_inf,
            eigenvalue_kind="k_inf",
            parameter_value=0.0,
            parameter_kind="k_inf_only",
            converged=True,
            metadata={
                "n_groups": n_groups,
                "method": method,
                "raw_result": k_inf,
            },
        )

    # ------------------------------------------------------------------
    # Flux reconstruction
    # ------------------------------------------------------------------

    def reconstruct_flux(
        self,
        *,
        n_panels: int = 256,
        z_eval: Optional[np.ndarray] = None,
    ) -> FluxSolution:
        r"""Reconstruct the scalar flux from F_N moments.

        For a bare-critical slab or sphere, the F_N moments
        :math:`(a_0, \ldots, a_N)` returned by
        :meth:`solve_critical` define the boundary angular flux
        :math:`\psi(\pm a, \mu)` as a polynomial in :math:`\mu`.
        Reconstructing the interior scalar flux :math:`\phi(z)`
        requires evaluating the Peierls integral with a log-singular
        kernel :math:`E_1(\tau)`. The choice of quadrature for the
        singular kernel sets the achievable accuracy:

        * ``"atkinson_nystrom"`` — Atkinson 1972/1997 §4.2 product-
          Simpson rule. **Recommended; fixes ERR-036.**

        * ``"legacy_gl"`` — plain Gauss–Legendre on the (z, μ)
          tensor product. Diagnostic only — saturates at 1–7\,\%.

        * ``"none"`` — :meth:`reconstruct_flux` raises.

        Parameters
        ----------
        n_panels : int, default 256
            Number of Simpson panels (atkinson_nystrom) or GL nodes
            (legacy_gl).
        z_eval : np.ndarray | None
            Spatial nodes (in mfp) at which to evaluate.
        """
        if isinstance(self._body, ReflectedSlab):
            refuse_unserved(
                "flux reconstruction on a symmetric reflected slab",
                owner="MomentSpace",
                missing="the two-media Peierls integral is not built; only the "
                "critical core half-thickness is",
            )
        if self.flux_reconstruction == "none":
            raise ValueError(
                "MomentSpace was constructed with flux_reconstruction='none'; "
                "reconstruct_flux is not callable."
            )

        # Solve criticality first to get the F_N expansion coefficients.
        critical = self.solve_critical()
        raw = critical.metadata.get("raw_result")
        if raw is None:
            raise RuntimeError(
                "MomentSpace.reconstruct_flux: solve_critical did not "
                "return a raw_result; cannot reconstruct without F_N "
                "expansion coefficients."
            )

        match self.geometry.coord:
            case CoordSystem.CARTESIAN:
                return self._reconstruct_slab(raw, n_panels, z_eval, critical)
            case CoordSystem.SPHERICAL:
                return self._reconstruct_sphere(raw, n_panels, z_eval, critical)
        raise NotImplementedError(  # pragma: no cover
            f"MomentSpace.reconstruct_flux: unhandled coordinate system "
            f"{self.geometry.coord!r}."
        )

    def _reconstruct_slab(
        self,
        raw,
        n_panels: int,
        z_eval: Optional[np.ndarray],
        critical: CriticalSolution,
    ) -> FluxSolution:
        from .slab.flux_reconstruction import (
            slab_scalar_flux_fn_projection,
            slab_scalar_flux_fn_projection_atkinson,
        )

        a = float(raw.a_critical_mfp)
        if z_eval is None:
            z_eval = np.linspace(-a, a, n_panels)
        z_eval_arr = np.asarray(z_eval, dtype=float)

        if self.flux_reconstruction == "atkinson_nystrom":
            phi = slab_scalar_flux_fn_projection_atkinson(
                raw, z_eval_arr, n_panels=n_panels
            )
        elif self.flux_reconstruction == "legacy_gl":
            phi = slab_scalar_flux_fn_projection(
                raw, z_eval_arr, n_quad_z=n_panels
            )
        else:
            raise ValueError(  # pragma: no cover
                f"unhandled flux_reconstruction strategy "
                f"{self.flux_reconstruction!r}"
            )

        return FluxSolution(
            spatial_nodes=z_eval_arr,
            scalar_flux=np.asarray(phi, dtype=float),
            angular_flux=None,
            angular_nodes=None,
            eigenvalue=critical.eigenvalue,
            eigenvalue_kind=critical.eigenvalue_kind,
            spatial_units="mfp",
            metadata={
                "geometry": "slab",
                "flux_reconstruction": self.flux_reconstruction,
                "n_panels": n_panels,
                "n_modes": self.fn_order,
                "a_critical_mfp": a,
                "c": critical.metadata.get("c"),
            },
        )

    def _reconstruct_sphere(
        self,
        raw,
        n_panels: int,
        z_eval: Optional[np.ndarray],
        critical: CriticalSolution,
    ) -> FluxSolution:
        from .sphere.flux_reconstruction import (
            solve_kll_sphere_continuum_coefficient,
            sphere_scalar_flux_kll,
        )

        R = float(raw.R_critical_mfp)
        if z_eval is None:
            z_eval = np.linspace(0.0, R, n_panels)
        z_eval_arr = np.asarray(z_eval, dtype=float)

        # Sphere flux reconstruction uses the KLL Fredholm path
        # (Path B in the fn_method theory page) — Path A.i for sphere
        # would require its own Atkinson product-Nyström treatment.
        c = float(critical.metadata.get("c", self.c))
        kll = solve_kll_sphere_continuum_coefficient(
            R, c, n_nodes=n_panels
        )
        phi = sphere_scalar_flux_kll(kll, z_eval_arr)

        return FluxSolution(
            spatial_nodes=z_eval_arr,
            scalar_flux=np.asarray(phi, dtype=float),
            angular_flux=None,
            angular_nodes=None,
            eigenvalue=critical.eigenvalue,
            eigenvalue_kind=critical.eigenvalue_kind,
            spatial_units="mfp",
            metadata={
                "geometry": "sphere",
                "flux_reconstruction": "kll_fredholm",
                "n_panels": n_panels,
                "n_modes": self.fn_order,
                "R_critical_mfp": R,
                "c": critical.metadata.get("c"),
            },
        )


__all__ = [
    "FluxReconstructionStrategy",
    "MomentSpace",
]

r"""Peierls integral equation reference for spherical CP verification.

Spherical specialisation of the unified polar-form Peierls Nyström
infrastructure in :mod:`orpheus.derivations.continuous.peierls.geometry`. This
module is a **registry-only façade** as of Issue #138 (2026-04-29),
mirror of :mod:`peierls_cylinder`: it owns the
:func:`~orpheus.derivations.continuous.peierls.sphere._build_peierls_sphere_case`
and
:func:`~orpheus.derivations.continuous.peierls.sphere._build_peierls_sphere_hollow_f4_case`
continuous-reference constructors and the
:data:`~orpheus.derivations.continuous.peierls.sphere.GEOMETRY` singleton
binding the canonical
:data:`~orpheus.derivations.continuous.peierls.geometry.SPHERE_1D`. The
``solve_peierls_sphere_{1g,mg}`` wrappers and the shape-specific
``PeierlsSphereSolution`` dataclass were retired in commit 99a05ab
— see :ref:`theory-peierls-api-posture` "Retired wrappers
(Issue #138)" for the migration recipe. Everything else — volume-kernel
assembly, Lagrange basis, angular/radial composite quadrature,
white-BC closures, eigenvalue power iteration — lives in
:mod:`~orpheus.derivations.continuous.peierls.geometry` and dispatches through
:class:`~orpheus.derivations.continuous.peierls.geometry.CurvilinearGeometry`
with ``kind = "sphere-1d"``.

See :doc:`/theory/peierls_nystrom` for the end-to-end derivation of
the unified structure and :doc:`/theory/collision_probability` for
the sphere-specific narrative (3-D point kernel, :math:`\sin\theta`
angular measure, rank-1 white-BC limitations paralleling the
cylinder).

The spherical Peierls equation in observer-centred polar coords:

.. math::

   \Sigma_t(r)\,\varphi(r)
     \;=\; \frac{\Sigma_t(r)}{2}\!
       \int_{0}^{\pi}\!\sin\theta\,\mathrm d\theta\!
       \int_{0}^{\rho_{\max}(r,\theta)}\!\!
         e^{-\tau(r,\rho,\theta)}\,
         q\bigl(r'(r,\rho,\theta)\bigr)\,\mathrm d\rho
     + S_{\rm bc}(r).

The prefactor :math:`1/2` absorbs the :math:`1/(4\pi)` of the 3-D
Green's function and a factor of :math:`2\pi` from trivial azimuthal
integration (the source field is radially symmetric, so only the
polar angle :math:`\theta` matters):
:math:`1/(4\pi) \cdot 2\pi = 1/2`.

The :math:`\sin\theta` weight comes from the spherical solid-angle
element :math:`\mathrm d\Omega = \sin\theta\,\mathrm d\theta\,
\mathrm d\phi` — no :math:`\pm` folding is needed since
:math:`\sin\theta \ge 0` on :math:`[0, \pi]` and the integrand
already covers the full hemisphere of directions seen from the
observer.

.. note::

   The ray-geometry primitives :math:`\rho_{\max}(r,\theta)` and
   :math:`r'(r,\rho,\theta)` are IDENTICAL to the cylinder case —
   a 1-D radial domain bounded by a spherical shell of radius
   :math:`R` has the same chord algebra regardless of whether the
   surrounding field is 2-D-symmetric (cylinder) or 3-D-symmetric
   (sphere). The only geometry-specific ingredients are the kernel
   (:math:`e^{-\tau}` vs :math:`\mathrm{Ki}_1`) and the angular
   weight (:math:`\sin\theta` vs constant).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from . import geometry as _pg
from ...common.continuous_reference import (
    ContinuousReferenceSolution,
    ProblemSpec,
    Provenance,
)
from ...common.xs_library import LAYOUTS, get_mixture, get_xs


# ═══════════════════════════════════════════════════════════════════════
# Sphere geometry singleton (binds the unified infrastructure)
# ═══════════════════════════════════════════════════════════════════════

GEOMETRY = _pg.SPHERE_1D


# ═══════════════════════════════════════════════════════════════════════
# ContinuousReferenceSolution builder
# ═══════════════════════════════════════════════════════════════════════

_MAT_IDS_SPH = {1: [2]}


def _build_peierls_sphere_case(
    ng_key: str,
    n_regions: int,
    n_panels_per_region: int = 3,
    p_order: int = 5,
    n_theta: int = 20,
    n_rho: int = 20,
    n_phi: int = 20,
    precision_digits: int = 20,
) -> ContinuousReferenceSolution:
    """Build a Peierls-sphere reference matching a cp_sph1D case."""
    if ng_key != "1g" or n_regions != 1:
        raise NotImplementedError(
            f"peierls_sphere continuous reference currently supports "
            f"1G 1-region only; got ng_key={ng_key!r}, n_regions={n_regions}"
        )

    from ..flat_source_cp.sphere import _RADII

    layout = LAYOUTS[n_regions]
    ng = int(ng_key[0])
    radii = np.array(_RADII[n_regions], dtype=float)

    xs_list = [get_xs(region, ng_key) for region in layout]
    sig_t = np.array([xs["sig_t"][0] for xs in xs_list])
    sig_s = np.array([xs["sig_s"][0, 0] for xs in xs_list])
    nu_sig_f = np.array([(xs["nu"] * xs["sig_f"])[0] for xs in xs_list])

    sol = _pg.solve_peierls_1g(
        GEOMETRY, radii, sig_t, sig_s, nu_sig_f,
        boundary="white",
        n_panels_per_region=n_panels_per_region,
        p_order=p_order,
        n_angular=n_theta, n_rho=n_rho, n_surf_quad=n_phi,
        dps=precision_digits,
    )

    r_nodes = sol.r_nodes
    _, r_wts, _ = _pg.composite_gl_r(
        radii, n_panels_per_region, p_order, dps=precision_digits,
    )
    phi = sol.phi_values[:, 0]
    integral = GEOMETRY.shell_volume_integral(r_nodes, r_wts, phi)
    if abs(integral) > 1e-30:
        phi_normed = phi / integral
        sol = replace(sol, phi_values=phi_normed[:, np.newaxis])

    def phi_fn(x: np.ndarray, g: int = 0) -> np.ndarray:
        return sol.phi(x, g)

    mat_ids = _MAT_IDS_SPH[n_regions]
    materials = {
        mat_ids[i]: get_mixture(region, ng_key)
        for i, region in enumerate(layout)
    }

    return ContinuousReferenceSolution(
        name=f"peierls_sph1D_{ng}eg_{n_regions}rg",
        problem=ProblemSpec(
            materials=materials,
            geometry_type="sphere-1d",
            geometry_params={
                "radius": float(radii[-1]),
                "radii": radii.tolist(),
                "mat_ids": mat_ids,
            },
            boundary_conditions={"outer": "white"},
            is_eigenvalue=True,
            n_groups=ng,
        ),
        operator_form="integral-peierls",
        phi=phi_fn,
        k_eff=sol.k_eff,
        provenance=Provenance(
            citation=(
                "Case & Zweifel 1967 (bare-sphere critical-R); "
                "Hebert 2020 Ch. 3 §3.5 (curvilinear Peierls)"
            ),
            derivation_notes=(
                f"Polar (θ, ρ) Nyström via the unified "
                f"CurvilinearGeometry(kind='sphere-1d'). "
                f"{n_panels_per_region} panels × {p_order} GL points on "
                f"[0, R], n_θ = {n_theta}, n_ρ = {n_rho}. Exponential "
                f"kernel (no dimensional reduction needed — the 3-D "
                f"point kernel is already in 3-D), sin θ angular weight. "
                f"White BC via rank-1 Schur closure (radial symmetry "
                f"collapses the general N_θ block to a single scalar "
                f"J⁻ = J⁺)."
            ),
            sympy_expression=None,
            precision_digits=precision_digits,
        ),
        equation_labels=(
            "peierls-unified",
        ),
        vv_level="L1",
        description=(
            f"{ng}G {n_regions}-region spherical Peierls "
            f"(exp polar Nyström via unified geometry, rank-1 white BC)"
        ),
        tolerance="O(h²)",
    )


def _build_peierls_sphere_hollow_f4_case(
    r0_over_R: float,
    ng_key: str = "1g",
    n_panels_per_region: int = 3,
    p_order: int = 5,
    n_theta: int = 24,
    n_rho: int = 24,
    n_phi: int = 24,
    precision_digits: int = 20,
) -> ContinuousReferenceSolution:
    r"""Build a hollow-sphere Peierls F.4 continuous reference.

    Uses Stamm'ler IV Eq. 34 = Hébert 2009 Eq. 3.323 (see
    :math:numref:`hebert-3-323`) — scalar rank-2 per-face closure
    assembled through
    :func:`~orpheus.derivations.continuous.peierls.geometry.compute_hollow_sph_transmission`
    (bare :math:`\exp(-\tau)` kernel with explicit :math:`\theta`
    integration, :math:`W_{oi} = (R/r_0)^2\,W_{io}` reciprocity —
    distinct from the cylinder's first-power form).

    Residuals vs :math:`k_\infty` at default quadrature and XS:

    - :math:`r_0 / R = 0.1` → 0.4 % err
    - :math:`r_0 / R = 0.2` → 1.2 % err
    - :math:`r_0 / R = 0.3` → 3.3 % err

    The sphere's 3-10× tighter residuals relative to the cylinder at
    the same :math:`r_0/R` (cylinder: 1.4 %, 5.4 %, 13 %) reflect
    the sphere's higher SO(3) symmetry capturing more angular
    structure at the scalar mode. Rank-N per-face refinement was
    falsified as a path to improve beyond F.4 — see
    :ref:`peierls-rank-n-per-face-closeout` and research-log L21.

    The reference is forward-looking: no production CP solver case
    exists at ``inner_radius > 0`` today. A future hollow
    ``cp_sph1D_hollow_*`` consumer can load this reference via
    :func:`orpheus.derivations.reference_values.continuous_get`.
    """
    from ..flat_source_cp.sphere import _RADII

    n_regions = 1
    ng = int(ng_key[0])
    R_out = float(_RADII[n_regions][-1])
    r0 = float(r0_over_R) * R_out
    if not (0.0 < r0 < R_out):
        raise ValueError(
            f"r0_over_R must be in (0, 1); got {r0_over_R}"
        )

    layout = LAYOUTS[n_regions]
    xs_list = [get_xs(region, ng_key) for region in layout]
    # (n_regions, ng) per-region per-group arrays (n_regions == 1 here).
    sig_t = np.stack([np.asarray(xs["sig_t"], dtype=float) for xs in xs_list])
    sig_s = np.stack([np.asarray(xs["sig_s"], dtype=float) for xs in xs_list])
    nu_sig_f = np.stack([
        np.asarray(xs["nu"] * xs["sig_f"], dtype=float) for xs in xs_list
    ])
    chi = np.stack([np.asarray(xs["chi"], dtype=float) for xs in xs_list])

    radii = np.array([R_out])

    geometry = _pg.CurvilinearGeometry(
        kind="sphere-1d", inner_radius=r0,
    )
    sol = _pg.solve_peierls_mg(
        geometry, radii, sig_t, sig_s, nu_sig_f, chi,
        boundary="white_f4",
        n_panels_per_region=n_panels_per_region,
        p_order=p_order,
        n_angular=n_theta, n_rho=n_rho, n_surf_quad=n_phi,
        dps=precision_digits,
    )

    r_nodes = sol.r_nodes
    _, r_wts, _ = _pg.composite_gl_r(
        radii, n_panels_per_region, p_order, dps=precision_digits,
        inner_radius=r0,
    )
    # Normalise each group's flux to unit shell-volume integral.
    phi_normed = np.empty_like(sol.phi_values)
    for g in range(ng):
        phi_g = sol.phi_values[:, g]
        integral_g = GEOMETRY.shell_volume_integral(r_nodes, r_wts, phi_g)
        phi_normed[:, g] = (
            phi_g / integral_g if abs(integral_g) > 1e-30 else phi_g
        )
    sol = replace(sol, phi_values=phi_normed)

    def phi_fn(x: np.ndarray, g: int = 0) -> np.ndarray:
        return sol.phi(x, g)

    mat_ids = _MAT_IDS_SPH[n_regions]
    materials = {
        mat_ids[i]: get_mixture(region, ng_key)
        for i, region in enumerate(layout)
    }

    r0_tag = f"{int(round(r0_over_R * 100)):02d}"
    return ContinuousReferenceSolution(
        name=f"peierls_sph1D_hollow_{ng}eg_{n_regions}rg_r0_{r0_tag}",
        problem=ProblemSpec(
            materials=materials,
            geometry_type="sphere-1d",
            geometry_params={
                "radius": R_out,
                "inner_radius": r0,
                "radii": radii.tolist(),
                "mat_ids": mat_ids,
            },
            boundary_conditions={"outer": "white_rank2"},
            is_eigenvalue=True,
            n_groups=ng,
        ),
        operator_form="integral-peierls",
        phi=phi_fn,
        k_eff=sol.k_eff,
        provenance=Provenance(
            citation=(
                "Stamm'ler & Abbate 1983 Ch. IV Eq. 34; "
                "Hébert 2009 Ch. 3 §3.8.4 Eq. 3.323"
            ),
            derivation_notes=(
                f"F.4 scalar rank-2 per-face closure on hollow sphere "
                f"(r_0/R = {r0_over_R:.2f}). Polar (θ, ρ) Nyström via "
                f"CurvilinearGeometry(kind='sphere-1d', "
                f"inner_radius={r0:.6f}). Transmission matrix built by "
                f"compute_hollow_sph_transmission (bare exp(-τ) kernel, "
                f"explicit θ integration), R_eff = (I - W)^(-1). "
                f"{n_panels_per_region} panels × {p_order} GL points on "
                f"[r_0, R], n_θ = {n_theta}, n_ρ = {n_rho}, "
                f"n_surf_quad = {n_phi}. White BC."
            ),
            sympy_expression=None,
            precision_digits=precision_digits,
        ),
        equation_labels=(
            "hebert-3-323",
            "peierls-unified",
        ),
        vv_level="L1",
        description=(
            f"{ng}G {n_regions}-region hollow spherical Peierls "
            f"(F.4 rank-2 per-face, exp kernel, r_0/R = {r0_over_R:.2f})"
        ),
        tolerance=f"O(h²) + scalar-mode residual ~{_F4_SPH_TOL[r0_over_R]}",
    )


# Measured baseline at default quadrature (see test_hollow_sph_rank2_beats_rank1_mark).
_F4_SPH_TOL = {0.1: "0.4 %", 0.2: "1.2 %", 0.3: "3.3 %"}



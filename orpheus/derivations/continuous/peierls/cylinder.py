r"""Peierls integral equation reference for cylindrical CP verification.

Cylindrical specialisation of the unified polar-form Peierls Nyström
infrastructure in :mod:`orpheus.derivations.continuous.peierls.geometry`. This
module is a **registry-only façade** as of Issue #138 (2026-04-29):
it owns the
:func:`~orpheus.derivations.continuous.peierls.cylinder._build_peierls_cylinder_case`
and
:func:`~orpheus.derivations.continuous.peierls.cylinder._build_peierls_cylinder_hollow_f4_case`
continuous-reference constructors and the
:data:`~orpheus.derivations.continuous.peierls.cylinder.GEOMETRY` singleton
binding the canonical
:data:`~orpheus.derivations.continuous.peierls.geometry.CYLINDER_1D`. The
``solve_peierls_cylinder_{1g,mg}`` wrappers and the shape-specific
``PeierlsCylinderSolution`` dataclass were retired in commit
99a05ab — see :ref:`theory-peierls-api-posture` "Retired wrappers
(Issue #138)" for the migration recipe. Everything else — volume-kernel
assembly, Lagrange basis, angular/radial composite quadrature,
white-BC closures, eigenvalue power iteration — lives in
:mod:`~orpheus.derivations.continuous.peierls.geometry` and dispatches through
:class:`~orpheus.derivations.continuous.peierls.geometry.CurvilinearGeometry`
with ``kind = "cylinder-1d"``.

See :doc:`/theory/peierls_nystrom` for the end-to-end derivation of
the unified structure. The cylinder-specific narrative (polar form,
Jacobian cancellation, rank-1 white-BC limitations) remains in
:doc:`/theory/collision_probability`.

The cylindrical Peierls equation in observer-centred polar coords:

.. math::

   \Sigma_t(r)\,\varphi(r)
     \;=\; \frac{\Sigma_t(r)}{\pi}\!
       \int_{0}^{\pi}\!\mathrm d\beta\!
       \int_{0}^{\rho_{\max}(r,\beta)}\!\!
         \mathrm{Ki}_1\!\bigl(\tau(r,\rho,\beta)\bigr)\,
         q\bigl(r'(r,\rho,\beta)\bigr)\,\mathrm d\rho
     + S_{\rm bc}(r).

The prefactor :math:`1/\pi` absorbs the :math:`1/(2\pi)` of the 2-D
Green's function (from integrating the 3-D point kernel over the
infinite z-axis) and a factor of 2 from folding
:math:`\beta \in [0, 2\pi]` to :math:`[0, \pi]` via the radially-
symmetric :math:`\pm\beta` reflection.
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
# Cylinder geometry singleton (binds the unified infrastructure)
# ═══════════════════════════════════════════════════════════════════════

# Production cylinder path: ``CYLINDER_1D`` with the natural
# :math:`\mathrm{Ki}_1` kernel. The cylinder-polar variant (explicit
# out-of-plane :math:`\varphi`-quadrature, formerly the production
# path) was archived 2026-04-19 as mathematically equivalent — see
# :file:`derivations/archive/peierls_cylinder_polar_assembly.py`.
GEOMETRY = _pg.CYLINDER_1D


# ═══════════════════════════════════════════════════════════════════════
# ContinuousReferenceSolution builder
# ═══════════════════════════════════════════════════════════════════════

_MAT_IDS_CYL = {1: [2]}


def _build_peierls_cylinder_case(
    ng_key: str,
    n_regions: int,
    n_panels_per_region: int = 3,
    p_order: int = 5,
    n_beta: int = 20,
    n_rho: int = 20,
    n_phi: int = 20,
    precision_digits: int = 20,
) -> ContinuousReferenceSolution:
    """Build a Peierls-cylinder reference matching a cp_cyl1D case."""
    if ng_key != "1g" or n_regions != 1:
        raise NotImplementedError(
            f"peierls_cylinder continuous reference currently supports "
            f"1G 1-region only; got ng_key={ng_key!r}, n_regions={n_regions}"
        )

    from ..flat_source_cp.cylinder import _RADII

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
        n_angular=n_beta, n_rho=n_rho, n_surf_quad=n_phi,
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

    mat_ids = _MAT_IDS_CYL[n_regions]
    materials = {
        mat_ids[i]: get_mixture(region, ng_key)
        for i, region in enumerate(layout)
    }

    return ContinuousReferenceSolution(
        name=f"peierls_cyl1D_{ng}eg_{n_regions}rg",
        problem=ProblemSpec(
            materials=materials,
            geometry_type="cylinder-1d",
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
                "Sanchez & McCormick 1982 (NSE 80) §IV.A; "
                "Hebert 2020 Ch. 3 §3.5"
            ),
            derivation_notes=(
                f"Polar (β, ρ) Nyström via the unified "
                f"CurvilinearGeometry(kind='cylinder-1d'). "
                f"{n_panels_per_region} panels × {p_order} GL points on "
                f"[0, R], n_β = {n_beta}, n_ρ = {n_rho}. Ki₁ kernel, "
                f"no singular Jacobian. White BC via rank-1 Schur "
                f"closure (radial symmetry collapses the general "
                f"N_β block to a single scalar J⁻ = J⁺)."
            ),
            sympy_expression=None,
            precision_digits=precision_digits,
        ),
        equation_labels=(
            "peierls-cylinder-equation",
            "peierls-cylinder-polar",
            "peierls-cylinder-ray-optical-depth",
            "peierls-unified",
        ),
        vv_level="L1",
        description=(
            f"{ng}G {n_regions}-region cylindrical Peierls "
            f"(Ki₁ polar Nyström via unified geometry, rank-1 white BC)"
        ),
        tolerance="O(h²)",
    )


def _build_peierls_cylinder_hollow_f4_case(
    r0_over_R: float,
    ng_key: str = "1g",
    n_panels_per_region: int = 3,
    p_order: int = 5,
    n_beta: int = 24,
    n_rho: int = 24,
    n_phi: int = 24,
    precision_digits: int = 20,
) -> ContinuousReferenceSolution:
    r"""Build a hollow-cylinder Peierls F.4 continuous reference.

    Uses Stamm'ler IV Eq. 34 = Hébert 2009 Eq. 3.323 (see
    :math:numref:`hebert-3-323`) — scalar rank-2 per-face closure
    assembled through
    :func:`~orpheus.derivations.continuous.peierls.geometry.compute_hollow_cyl_transmission`
    (Ki\ :sub:`3` Bickley fold, :math:`W_{oi} = (R/r_0)\,W_{io}`
    reciprocity — distinct from the sphere's squared form).

    **1G residuals** vs :math:`k_\infty` at default quadrature and
    region-B XS:

    - :math:`r_0 / R = 0.1` → 1.4 % err
    - :math:`r_0 / R = 0.2` → 5.4 % err
    - :math:`r_0 / R = 0.3` → 13 % err

    The absolute residuals are looser than sphere F.4 (0.4 %, 1.2 %,
    3.3 % at the same cavity ratios) because the cylinder's
    out-of-plane :math:`\theta` fold into Ki\ :sub:`3` averages away
    angular information the sphere's explicit 3-D integration
    retains. Rank-N per-face refinement was falsified as a path to
    improve beyond F.4 — see
    :ref:`peierls-rank-n-per-face-closeout` and research-log L21.

    The reference is forward-looking: no production CP solver case
    exists at ``inner_radius > 0`` today. A future hollow
    ``cp_cyl1D_hollow_*`` consumer can load this reference via
    :func:`orpheus.derivations.reference_values.continuous_get`.

    Parameters
    ----------
    r0_over_R
        Cavity ratio. Must be in ``(0, 1)``.
    ng_key
        XS-library group key (``"1g"``, ``"2g"``, ...). Added in
        Issue #104 (2026-04-24); the MG path routes through the
        canonical :func:`peierls_geometry.solve_peierls_mg` with
        :class:`CurvilinearGeometry` bound for the hollow case.
        ``"1g"`` preserves the legacy residuals listed above;
        ``"2g"`` runs a fresh benchmark per commit 2 of the Issue
        #104 plan.
    """
    from ..flat_source_cp.cylinder import _RADII

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
    # (n_regions, ng) per-region per-group arrays. The single-region
    # hollow cell has n_regions=1, so these are shape (1, ng).
    sig_t = np.stack([np.asarray(xs["sig_t"], dtype=float) for xs in xs_list])
    sig_s = np.stack([np.asarray(xs["sig_s"], dtype=float) for xs in xs_list])
    nu_sig_f = np.stack([
        np.asarray(xs["nu"] * xs["sig_f"], dtype=float) for xs in xs_list
    ])
    chi = np.stack([np.asarray(xs["chi"], dtype=float) for xs in xs_list])

    radii = np.array([R_out])

    geometry = _pg.CurvilinearGeometry(
        kind="cylinder-1d", inner_radius=r0,
    )
    sol = _pg.solve_peierls_mg(
        geometry, radii, sig_t, sig_s, nu_sig_f, chi,
        boundary="white_f4",
        n_panels_per_region=n_panels_per_region,
        p_order=p_order,
        n_angular=n_beta, n_rho=n_rho, n_surf_quad=n_phi,
        dps=precision_digits,
    )

    r_nodes = sol.r_nodes
    _, r_wts, _ = _pg.composite_gl_r(
        radii, n_panels_per_region, p_order, dps=precision_digits,
        inner_radius=r0,
    )
    # Normalise each group's flux to unit shell-volume integral.
    # Preserves the legacy 1G behaviour (group-0 integral = 1) and
    # extends it per-group for MG.
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

    mat_ids = _MAT_IDS_CYL[n_regions]
    materials = {
        mat_ids[i]: get_mixture(region, ng_key)
        for i, region in enumerate(layout)
    }

    r0_tag = f"{int(round(r0_over_R * 100)):02d}"
    return ContinuousReferenceSolution(
        name=f"peierls_cyl1D_hollow_{ng}eg_{n_regions}rg_r0_{r0_tag}",
        problem=ProblemSpec(
            materials=materials,
            geometry_type="cylinder-1d",
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
                f"F.4 scalar rank-2 per-face closure on hollow cylinder "
                f"(r_0/R = {r0_over_R:.2f}). Polar (β, ρ) Nyström via "
                f"CurvilinearGeometry(kind='cylinder-1d', "
                f"inner_radius={r0:.6f}). Transmission matrix built by "
                f"compute_hollow_cyl_transmission (Ki_3 Bickley fold), "
                f"R_eff = (I - W)^(-1). "
                f"{n_panels_per_region} panels × {p_order} GL points on "
                f"[r_0, R], n_β = {n_beta}, n_ρ = {n_rho}, "
                f"n_surf_quad = {n_phi}. White BC."
            ),
            sympy_expression=None,
            precision_digits=precision_digits,
        ),
        equation_labels=(
            "hebert-3-323",
            "peierls-cylinder-equation",
            "peierls-cylinder-polar",
            "peierls-unified",
        ),
        vv_level="L1",
        description=(
            f"{ng}G {n_regions}-region hollow cylindrical Peierls "
            f"(F.4 rank-2 per-face, Ki_3 fold, r_0/R = {r0_over_R:.2f})"
        ),
        tolerance=f"O(h²) + scalar-mode residual ~{_F4_CYL_TOL[r0_over_R]}",
    )


# Measured baseline at default quadrature (see test_hollow_cyl_rank2_beats_rank1_mark).
_F4_CYL_TOL = {0.1: "1.4 %", 0.2: "5.4 %", 0.3: "13 %"}



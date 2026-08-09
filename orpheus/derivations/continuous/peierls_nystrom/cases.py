r"""Unified continuous-reference registry for Peierls Nyström solvers,
organized by **orbit-space M/G class** instead of shape.

See :file:`.claude/plans/topology-based-consolidation.md` and Sphinx
§\ ``orbit-space-m-g-classification`` (the structural signature),
§\ ``theory-peierls-capabilities``, and §\ ``theory-peierls-naming``.

The two orbit-space classes are determined by the number of physical
boundary endpoints of the 1-D orbit space M/G (where M is the 3-D
problem domain and G is the symmetry group acting by isometries):

- **Class A — two-surface** (F.4 applies). M/G is a 1-D interval with
  two physical BC endpoints. Members: slab (two parallel faces, M/G
  = R³ modulo R²-translation), hollow annular cylinder (inner +
  outer ring, M/G = R³ modulo R-translation × SO(2) with inner-radius
  cut), hollow sphere (inner + outer shell, M/G = R³ modulo SO(3)
  with inner-radius cut). Shared closure class: Stamm'ler IV Eq. 34
  = Hébert 2009 Eq. 3.323 (scalar rank-2 per-face).
- **Class B — one-surface compact** (rank-1 Mark only). M/G is a 1-D
  interval with one physical BC endpoint (the inner :math:`r=0` is a
  coordinate singularity, not a physical surface). Members: solid
  cylinder (M/G = R³ modulo R-translation × SO(2)), solid sphere
  (M/G = R³ modulo SO(3)). F.4 structurally collapses.

This module is the canonical entry point for continuous-reference
registration via :func:`cases`. The per-geometry modules
(:mod:`~orpheus.derivations.continuous.peierls_nystrom.slab`,
:mod:`~orpheus.derivations.continuous.peierls_nystrom.cylinder`,
:mod:`~orpheus.derivations.continuous.peierls_nystrom.sphere`) retain their
``_build_*_case`` constructor functions — this module calls them
directly. Their ``continuous_cases()`` hooks return empty lists to
avoid double-registration; the registry-builder's auto-discovery
walks every module and this module is the single source for Peierls
continuous references.

Slab note (2026-04-24): slab has two independent verification paths:

1. **Native E₁ Nyström** (:mod:`~orpheus.derivations.continuous.peierls_nystrom.slab`) —
   classical singularity-subtraction + product-integration, multi-
   group via a block-Toeplitz assembly. Retained as an independent
   cross-check.
2. **Unified curvilinear** (:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_mg`
   with :data:`~orpheus.derivations.continuous.peierls_nystrom.geometry.SLAB_POLAR_1D`)
   — adaptive ``mpmath.quad`` with forced :math:`\mu = 0` breakpoint,
   machine precision by construction (see Phase G — Sphinx
   §theory-peierls-slab-polar).

Both paths are multi-group capable as of Issue #104 (2026-04-24).
Phase G.5 routing activation (Issue #130, 2026-04-24): after
`Issue #131 <https://github.com/deOliveira-R/ORPHEUS/issues/131>`_
diagnosed the 1.5 % discrepancy to closed-form-avoidance in the
multi-region slab branches of ``compute_P_esc_{outer,inner}`` and
``compute_G_bc_{outer,inner}`` (finite-N GL over the µ-integral
when the integral has a closed form
:math:`\tfrac12\,E_2(\tau_{\rm total})`), the fix was applied and
the unified path now matches native E₁ **bit-exactly** on the
shipped ``peierls_slab_2eg_2rg`` fixture
(``rel_diff = 5.4e-16`` at ``n_panels_per_region=2, p_order=3,
dps=20``). The ``_SLAB_VIA_UNIFIED`` flag now **defaults to True**;
set ``ORPHEUS_SLAB_VIA_E1=1`` to force the native path for
bisection. Both paths remain exercised by the test suite: the
unified path is the shipped registry route, and the native path is
exercised by every gate that calls
:func:`~orpheus.derivations.continuous.peierls_nystrom.slab.solve_peierls_eigenvalue`
directly — :mod:`tests.derivations.test_peierls_convergence` (L0
self-convergence under panel refinement),
:mod:`tests.derivations.test_peierls_multigroup` (including the
diagnostic test
:class:`tests.derivations.test_peierls_multigroup.TestSlabViaUnifiedDiscrepancyDiagnostic`,
now at ``rel_diff < 1e-10`` bound), and
:mod:`tests.derivations.test_peierls_greens_function_slab_solver`.
"""
from __future__ import annotations

import os as _os
from collections.abc import Callable

from ...common.continuous_reference import ContinuousReferenceSolution

# Issue #130 Phase G.5 routing switch. Defaults to True (unified
# path) as of 2026-04-24 — see module docstring for the benchmark
# that unblocked activation. ``ORPHEUS_SLAB_VIA_E1=1`` overrides to
# the native E₁ Nyström for bisection / testing.
_SLAB_VIA_UNIFIED: bool = (
    _os.environ.get("ORPHEUS_SLAB_VIA_E1", "0") != "1"
)


# ---------------------------------------------------------------------
# Class A — two-surface (F.4 applies)
# ---------------------------------------------------------------------


def build_two_surface_case(
    shape: str,
    ng_key: str = "1g",
    n_regions: int = 1,
    *,
    inner_radius: float | None = None,
) -> ContinuousReferenceSolution:
    r"""Build a Class-A (two-surface) continuous reference.

    Class A members share the F.4 scalar rank-2 per-face closure
    (:math:numref:`hebert-3-323`). Dispatch on ``shape``:

    - ``"slab"`` — calls
      :func:`orpheus.derivations.continuous.peierls_nystrom.slab._build_peierls_slab_case`.
      ``inner_radius`` is ignored (slab has two parallel faces at
      :math:`x=0` and :math:`x=L`, not a cavity).
    - ``"cylinder-1d"`` — requires ``inner_radius > 0``; calls
      :func:`orpheus.derivations.continuous.peierls_nystrom.cylinder._build_peierls_cylinder_hollow_f4_case`.
    - ``"sphere-1d"`` — requires ``inner_radius > 0``; calls
      :func:`orpheus.derivations.continuous.peierls_nystrom.sphere._build_peierls_sphere_hollow_f4_case`.

    Parameters
    ----------
    shape
        ``"slab"``, ``"cylinder-1d"``, or ``"sphere-1d"``.
    ng_key
        XS-library group key (``"1g"``, ``"2g"``, ``"4g"``). The
        hollow curvilinear cases currently support ``"1g"`` only;
        multi-group lift is Issue #104.
    n_regions
        Number of radial regions. Hollow curvilinear cases support
        1-region only today (single annular shell). Slab supports
        1/2/4.
    inner_radius
        Cavity radius :math:`r_0` for curvilinear hollow cases.
        **Required** for ``"cylinder-1d"`` / ``"sphere-1d"``; must
        be strictly between 0 and the outer radius. Ignored for slab.

    Raises
    ------
    ValueError
        For curvilinear shapes when ``inner_radius`` is missing or
        not in ``(0, R_outer)``. Use
        :func:`build_one_surface_compact_case` for solid geometry.
    """
    if shape == "slab":
        if _SLAB_VIA_UNIFIED:
            return _build_peierls_slab_case_via_unified(ng_key, n_regions)
        from .slab import _build_peierls_slab_case
        return _build_peierls_slab_case(
            ng_key, n_regions,
        )
    if shape == "cylinder-1d":
        if inner_radius is None:
            raise ValueError(
                "cylinder-1d is a Class A (two-surface) case only "
                "when inner_radius > 0. Use build_one_surface_compact_case "
                "for solid cylinder."
            )
        from .cylinder import _build_peierls_cylinder_hollow_f4_case
        # The hollow-f4 builder takes r_0_over_R (unitless), not
        # absolute inner_radius. The single shipped outer radius for
        # 1g 1-region is R=1 (from cp_cylinder._RADII[1][-1]) so the
        # two are numerically equal when R=1 — but be explicit.
        from ..flat_source_cp.cylinder import _RADII as _CYL_RADII
        R_out = float(_CYL_RADII[n_regions][-1])
        return _build_peierls_cylinder_hollow_f4_case(
            r0_over_R=float(inner_radius) / R_out,
            ng_key=ng_key,
        )
    if shape == "sphere-1d":
        if inner_radius is None:
            raise ValueError(
                "sphere-1d is a Class A (two-surface) case only "
                "when inner_radius > 0. Use build_one_surface_compact_case "
                "for solid sphere."
            )
        from .sphere import _build_peierls_sphere_hollow_f4_case
        from ..flat_source_cp.sphere import _RADII as _SPH_RADII
        R_out = float(_SPH_RADII[n_regions][-1])
        return _build_peierls_sphere_hollow_f4_case(
            r0_over_R=float(inner_radius) / R_out,
            ng_key=ng_key,
        )
    raise ValueError(
        f"build_two_surface_case: unknown shape {shape!r}; "
        f"expected 'slab', 'cylinder-1d', or 'sphere-1d'"
    )


# ---------------------------------------------------------------------
# Phase G.5 — slab routing through the unified adaptive-mpmath path
# (Issue #130). Default-off; enabled by ``_SLAB_VIA_UNIFIED``.
# ---------------------------------------------------------------------


def _build_peierls_slab_case_via_unified(
    ng_key: str,
    n_regions: int,
    n_panels_per_region: int = 16,
    p_order: int = 6,
    precision_digits: int = 30,
) -> ContinuousReferenceSolution:
    r"""Build the slab continuous reference through the unified
    :func:`peierls_geometry.solve_peierls_mg` path
    (:data:`peierls_geometry.SLAB_POLAR_1D` + adaptive ``mpmath.quad``
    with forced :math:`\mu = 0` breakpoint).

    Mirrors :func:`orpheus.derivations.continuous.peierls_nystrom.slab._build_peierls_slab_case`
    on inputs and on the ``ContinuousReferenceSolution`` output
    schema. The difference is the K-matrix assembly route:

    - Native path: classical :math:`E_1` Nyström with singularity
      subtraction + product integration (fast, O(h²) convergence).
    - Unified path: adaptive ``mpmath.quad`` on observer-centred
      polar coords, one adaptive double-quad per K element
      (verification-primitive precision, O(N²) adaptive cost).

    **Not the default** for the shipped reference as of Issue #130
    (2026-04-24). Benchmark at modest quadrature shows a ~1.5 %
    rel_diff on the ``peierls_slab_2eg_2rg`` fixture — too large for
    a shipped L1 reference. Enable via ``ORPHEUS_SLAB_VIA_UNIFIED=1``
    for bisection / testing; default routing stays on the native
    path until the discrepancy is resolved.
    """
    import numpy as _np

    from ...common.xs_library import LAYOUTS, get_mixture, get_xs
    from ...common.continuous_reference import ProblemSpec, Provenance
    from ..flat_source_cp.slab import _THICKNESSES
    from .geometry import SLAB_POLAR_1D, solve_peierls_mg
    from .slab import _MAT_IDS

    layout = LAYOUTS[n_regions]
    ng = int(ng_key[0])
    thicknesses = _THICKNESSES[n_regions]

    xs_list = [get_xs(region, ng_key) for region in layout]

    # (n_regions, ng) per-region per-group arrays.
    sig_t = _np.stack(
        [_np.asarray(xs["sig_t"], dtype=float) for xs in xs_list]
    )
    sig_s = _np.stack(
        [_np.asarray(xs["sig_s"], dtype=float) for xs in xs_list]
    )
    nu_sig_f = _np.stack(
        [_np.asarray(xs["nu"] * xs["sig_f"], dtype=float) for xs in xs_list]
    )
    chi = _np.stack(
        [_np.asarray(xs["chi"], dtype=float) for xs in xs_list]
    )

    # Cumulative outer-radius boundaries (slab thicknesses → radii).
    radii = _np.cumsum(_np.asarray(thicknesses, dtype=float))

    sol = solve_peierls_mg(
        SLAB_POLAR_1D, radii,
        sig_t=sig_t, sig_s=sig_s, nu_sig_f=nu_sig_f, chi=chi,
        boundary="white_f4",
        n_panels_per_region=n_panels_per_region,
        p_order=p_order,
        dps=precision_digits,
        tol=10 ** -(precision_digits - 5),
    )

    def phi_fn(x: _np.ndarray, g: int = 0) -> _np.ndarray:
        return sol.phi(x, g)

    mat_ids = _MAT_IDS[n_regions]
    materials = {
        mat_ids[i]: get_mixture(region, ng_key)
        for i, region in enumerate(layout)
    }

    return ContinuousReferenceSolution(
        name=f"peierls_slab_{ng}eg_{n_regions}rg",
        problem=ProblemSpec(
            materials=materials,
            geometry_type="slab",
            geometry_params={
                "length": sum(thicknesses),
                "thicknesses": thicknesses,
                "mat_ids": mat_ids,
            },
            boundary_conditions={"left": "white", "right": "white"},
            is_eigenvalue=True,
            n_groups=ng,
        ),
        operator_form="integral-peierls",
        phi=phi_fn,
        k_eff=sol.k_eff,
        provenance=Provenance(
            citation=(
                "Case & Zweifel 1967 Ch. 4; "
                "Kress 2014 (Nyström for Fredholm); "
                "Phase G §theory-peierls-slab-polar"
            ),
            derivation_notes=(
                f"Unified-path slab Peierls via solve_peierls_mg "
                f"(SLAB_POLAR_1D, adaptive mpmath.quad with forced "
                f"µ=0 breakpoint). {n_panels_per_region} panels × "
                f"{p_order} GL points per region, white_f4 rank-2 "
                f"per-face F.4 closure. Issue #130 routing path."
            ),
            sympy_expression=None,
            precision_digits=precision_digits,
        ),
        equation_labels=("peierls-equation", "peierls-unified"),
        vv_level="L1",
        description=(
            f"{ng}G {n_regions}-region slab Peierls "
            f"(unified adaptive mpmath.quad, white_f4 BC)"
        ),
        tolerance="verification-primitive (see Sphinx §theory-peierls-multigroup)",
    )


# ---------------------------------------------------------------------
# Class B — one-surface compact (rank-1 Mark only)
# ---------------------------------------------------------------------


def build_one_surface_compact_case(
    shape: str,
    ng_key: str = "1g",
    n_regions: int = 1,
) -> ContinuousReferenceSolution:
    r"""Build a Class-B (one-surface compact) continuous reference.

    Class B members (solid cylinder, solid sphere) ship only the
    rank-1 Mark closure. F.4 collapses to rank-1 Mark on solid
    geometry (no second-face coupling).

    **No references are registered in Class B today.** The rank-1
    Mark floor (21 % err at :math:`R = 1` MFP for cylinder per
    Issue #103) is too loose to serve as an L1 reference for the
    ``cp_{cyl,sph}1D_*`` solver tests. Lifting the floor requires
    Issue #103 (rank-N DP\ :sub:`N` on the single outer face) or
    Issue #101 (chord-based Ki₁ analytical).

    This function exists for future use when one of those lands.
    Until then it unconditionally raises ``NotImplementedError``
    with an explanatory message.
    """
    raise NotImplementedError(
        f"build_one_surface_compact_case({shape!r}) is not yet "
        f"populated. Class B (solid cylinder / solid sphere) has "
        f"no shipped continuous references because the rank-1 Mark "
        f"floor is too loose (21 % err at R=1 MFP per Issue #103). "
        f"Resolution requires rank-N DP_N (Issue #103) or chord-"
        f"based Ki_1 analytical (Issue #101)."
    )


# ---------------------------------------------------------------------
# Registry entry point — auto-discovered by reference_values.py
# ---------------------------------------------------------------------


def _class_a_cases() -> list[ContinuousReferenceSolution]:
    """Class A — two-surface cases. Slab + hollow cylinder/sphere F.4.

    Multi-group hollow cyl/sph references were added in Issue #104
    (2026-04-24) once the unified :func:`peierls_geometry.solve_peierls_mg`
    path landed. Each ``r_0/R`` sweep entry now ships a 1G and 2G
    variant — the 1G residuals against :math:`k_\\infty` are
    reference-stable (1.4 % / 5.4 % / 13 % cyl; 0.4 % / 1.2 % / 3.3 %
    sph); the 2G variants inherit the same F.4 scalar rank-2 per-face
    closure applied group-wise.
    """
    refs: list[ContinuousReferenceSolution] = []
    # Slab: 2G 2-region (current shipped default — native E₁ Nyström
    # path per peierls_cases module docstring).
    refs.append(build_two_surface_case("slab", "2g", 2))
    # Hollow cylinder F.4 at r_0/R ∈ {0.1, 0.2, 0.3}, 1G and 2G variants.
    for r0 in (0.1, 0.2, 0.3):
        refs.append(build_two_surface_case(
            "cylinder-1d", "1g", 1, inner_radius=r0,
        ))
        refs.append(build_two_surface_case(
            "cylinder-1d", "2g", 1, inner_radius=r0,
        ))
    # Hollow sphere F.4 at r_0/R ∈ {0.1, 0.2, 0.3}, 1G and 2G variants.
    for r0 in (0.1, 0.2, 0.3):
        refs.append(build_two_surface_case(
            "sphere-1d", "1g", 1, inner_radius=r0,
        ))
        refs.append(build_two_surface_case(
            "sphere-1d", "2g", 1, inner_radius=r0,
        ))
    return refs


def _class_b_cases() -> list[ContinuousReferenceSolution]:
    """Class B — one-surface compact cases. Empty today."""
    return []


def continuous_cases() -> list[ContinuousReferenceSolution]:
    r"""All Peierls continuous references across both topology classes.

    Registered by auto-discovery in
    :func:`orpheus.derivations.reference_values._build_continuous_registry`
    via the standard ``continuous_cases()`` contract.
    """
    return _class_a_cases() + _class_b_cases()


# Alias for readers who want the topology-explicit name.
cases = continuous_cases


def continuous_case_builders() -> dict[str, Callable[[], ContinuousReferenceSolution]]:
    r"""Lazy ``{name: thunk}`` map for the Class-A Peierls references (Issue #212).

    The registry walker
    (:func:`orpheus.derivations.reference_values._build_continuous_registry`)
    **prefers this contract over** :func:`continuous_cases` when present: it
    records the reference *names* cheaply and defers each O(minutes)
    adaptive-``mpmath`` eigenvalue solve until
    :func:`orpheus.derivations.reference_values.continuous_get` actually
    requests that name. Fetching an unrelated reference (e.g. an SN MMS case)
    therefore no longer pays the Peierls build cost — the root cause of the
    apparent "hang" diagnosed in Issue #212.

    The thunks delegate to :func:`build_two_surface_case`, so each *built*
    reference's :attr:`~ContinuousReferenceSolution.name` is authoritative.
    The keys here are derived by the **same** formula the inner builders stamp
    onto that attribute (``peierls_slab_{ng}eg_{nr}rg`` and
    ``peierls_{cyl1D,sph1D}_hollow_{ng}eg_1rg_r0_{r0_tag}`` with
    ``r0_tag = round(r0/R_out · 100)``); the equivalence
    ``set(continuous_case_builders()) == {c.name for c in continuous_cases()}``
    is pinned by ``tests/derivations/test_continuous_registry_lazy.py`` so the
    cheap keys can never silently drift from the built names.

    Mirrors the :func:`_class_a_cases` loop exactly — any new shipped Class-A
    reference must be added in both places (the drift-guard test enforces it).
    """
    from functools import partial

    from ..flat_source_cp.cylinder import _RADII as _CYL_RADII
    from ..flat_source_cp.sphere import _RADII as _SPH_RADII

    builders: dict[str, Callable[[], ContinuousReferenceSolution]] = {}

    # Slab — 2G 2-region (the single shipped slab default).
    builders["peierls_slab_2eg_2rg"] = partial(build_two_surface_case, "slab", "2g", 2)

    # Hollow cylinder / sphere F.4 at r_0/R ∈ {0.1, 0.2, 0.3}, 1G + 2G.
    cyl_R_out = float(_CYL_RADII[1][-1])
    sph_R_out = float(_SPH_RADII[1][-1])
    for r0 in (0.1, 0.2, 0.3):
        cyl_tag = f"{int(round((r0 / cyl_R_out) * 100)):02d}"
        sph_tag = f"{int(round((r0 / sph_R_out) * 100)):02d}"
        for ng_key, ng in (("1g", 1), ("2g", 2)):
            builders[f"peierls_cyl1D_hollow_{ng}eg_1rg_r0_{cyl_tag}"] = partial(
                build_two_surface_case, "cylinder-1d", ng_key, 1, inner_radius=r0,
            )
            builders[f"peierls_sph1D_hollow_{ng}eg_1rg_r0_{sph_tag}"] = partial(
                build_two_surface_case, "sphere-1d", ng_key, 1, inner_radius=r0,
            )
    return builders


# ---------------------------------------------------------------------
# Metadata-only enumeration (no eigenvalue solves). Source of truth
# for the Sphinx §theory-peierls-capabilities matrix. Keep this list
# synchronised with the ``_class_{a,b}_cases()`` loops above — any
# new shipped reference must appear in both. See
# ``tools/verification/generate_capability_matrices.py`` for the
# consumer (the meta-generator that auto-discovers ``cases.py`` across
# every ``orpheus.derivations.continuous`` package).
# ---------------------------------------------------------------------


def capability_rows() -> list[dict[str, object]]:
    """Static metadata for every shipped Peierls continuous reference.

    Returns one dict per registered reference. The schema contract is
    documented in
    :mod:`tools.verification.generate_capability_matrices` (required:
    ``name``, ``geometry``, ``n_groups``, ``n_regions``, ``bc``,
    ``status``; optional auto-detected: ``r0_over_R``, ``closure``,
    ``accuracy``, ``scattering_order``, ``multiplying``,
    ``orbit_space_class``).

    This function does **not** call any eigenvalue solver. It is safe
    to invoke at Sphinx build time without paying the O(minutes) cost
    of :func:`continuous_cases`. The capability-matrix infrastructure
    test
    :func:`tests.derivations.test_capability_matrices.test_check_mode_exits_zero_when_in_sync`
    pins the rendered include file against this registry.
    """
    # Lazy imports of per-shape tolerance tables so that this module
    # is still importable from doc-build contexts that may not have
    # every optional dep wired up.
    from .cylinder import _F4_CYL_TOL
    from .sphere import _F4_SPH_TOL

    f4_label = r":math:`{\rm F.4}` (Stamm'ler Eq. 34)"
    rank2_label = r"white rank-2 per-face (E\ :sub:`2`/E\ :sub:`3`)"
    # Slab rank-2 closure imposes WHITE BC at the outer face (see
    # ``slab/__init__.py``); F.4 imposes the curvilinear analogue
    # closing the cavity face.
    bc_white = "white (rank-2 per-face)"
    bc_f4 = "vacuum + F.4 cavity closure"
    shipped_status = "shipped (registry-anchored)"

    rows: list[dict[str, object]] = []

    # Class A — slab (single shipped entry: 2G, 2-region, native E₁
    # path OR unified-adaptive depending on ``_SLAB_VIA_UNIFIED``; the
    # closure class and matrix column are identical either way).
    rows.append({
        "name": "``peierls_slab_2eg_2rg``",
        "geometry": "slab",
        "n_groups": 2,
        "n_regions": 2,
        "r0_over_R": None,
        "closure": rank2_label,
        "accuracy": "O(h²), Wigner-Seitz exact",
        "orbit_space_class": "A",
        "bc": bc_white,
        "status": shipped_status,
    })

    # Class A — hollow cylinder F.4 at r_0/R ∈ {0.1, 0.2, 0.3} × {1G, 2G}.
    for r0 in (0.1, 0.2, 0.3):
        r0_tag = f"{int(round(r0 * 100)):02d}"
        tol_1g = _F4_CYL_TOL[r0]
        rows.append({
            "name": f"``peierls_cyl1D_hollow_1eg_1rg_r0_{r0_tag}``",
            "geometry": "cylinder-1d",
            "n_groups": 1,
            "n_regions": 1,
            "r0_over_R": r0,
            "closure": f4_label,
            "accuracy": f"~{tol_1g} structural (scalar mode)",
            "orbit_space_class": "A",
            "bc": bc_f4,
            "status": shipped_status,
        })
        rows.append({
            "name": f"``peierls_cyl1D_hollow_2eg_1rg_r0_{r0_tag}``",
            "geometry": "cylinder-1d",
            "n_groups": 2,
            "n_regions": 1,
            "r0_over_R": r0,
            "closure": f4_label,
            "accuracy": (
                f"2G builds, finite k_eff (``TestMG2GHollowRegistration``); "
                f"k_eff vs ``cp_cylinder`` analytical not yet gated; "
                f"structural residual expected ~{tol_1g} (group-local closure, "
                f"unverified) — Issue #104 AC"
            ),
            "orbit_space_class": "A",
            "bc": bc_f4,
            "status": "shipped (k_eff gate pending — Issue #104)",
        })

    # Class A — hollow sphere F.4 at r_0/R ∈ {0.1, 0.2, 0.3} × {1G, 2G}.
    for r0 in (0.1, 0.2, 0.3):
        r0_tag = f"{int(round(r0 * 100)):02d}"
        tol_1g = _F4_SPH_TOL[r0]
        rows.append({
            "name": f"``peierls_sph1D_hollow_1eg_1rg_r0_{r0_tag}``",
            "geometry": "sphere-1d",
            "n_groups": 1,
            "n_regions": 1,
            "r0_over_R": r0,
            "closure": f4_label,
            "accuracy": f"~{tol_1g} structural (scalar mode)",
            "orbit_space_class": "A",
            "bc": bc_f4,
            "status": shipped_status,
        })
        rows.append({
            "name": f"``peierls_sph1D_hollow_2eg_1rg_r0_{r0_tag}``",
            "geometry": "sphere-1d",
            "n_groups": 2,
            "n_regions": 1,
            "r0_over_R": r0,
            "closure": f4_label,
            "accuracy": (
                f"2G builds, finite k_eff (``TestMG2GHollowRegistration``); "
                f"k_eff vs ``cp_sphere`` analytical not yet gated; "
                f"structural residual expected ~{tol_1g} (group-local closure, "
                f"unverified) — Issue #104 AC"
            ),
            "orbit_space_class": "A",
            "bc": bc_f4,
            "status": "shipped (k_eff gate pending — Issue #104)",
        })

    # Class B — one-surface compact. No shipped references today
    # (rank-1 Mark floor is too loose; see Issues #101 / #103).

    return rows

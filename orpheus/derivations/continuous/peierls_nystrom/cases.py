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

from .naming import ShippedReference
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
    r0_over_R: float | None = None,
) -> ContinuousReferenceSolution:
    r"""Build a Class-A (two-surface) continuous reference.

    Class A members share the F.4 scalar rank-2 per-face closure
    (:math:numref:`hebert-3-323`). Dispatch on ``shape``:

    - ``"slab"`` — calls
      :func:`orpheus.derivations.continuous.peierls_nystrom.slab._build_peierls_slab_case`.
      ``r0_over_R`` is ignored (slab has two parallel faces at
      :math:`x=0` and :math:`x=L`, not a cavity).
    - ``"cylinder-1d"`` — requires ``r0_over_R``; calls
      :func:`orpheus.derivations.continuous.peierls_nystrom.cylinder._build_peierls_cylinder_hollow_f4_case`.
    - ``"sphere-1d"`` — requires ``r0_over_R``; calls
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
    r0_over_R
        The **dimensionless** cavity ratio :math:`r_0/R` for curvilinear
        hollow cases — **required** for ``"cylinder-1d"`` /
        ``"sphere-1d"``, ignored for slab. Forwarded to the hollow-F.4
        builders unchanged, because a ratio is what they take.

        .. note:: **This was ``inner_radius``, an absolute length, until
           2026-08-09 (#345).** Every caller and every docstring in the
           family described the sweep as :math:`r_0/R \in \{0.1, 0.2,
           0.3\}` — a ratio — while the parameter declared a length and
           the body divided by ``R_out`` to recover the ratio. The two
           coincided only at :math:`R = 1`, which every entry of the
           shipped ``_RADII`` tables happens to be, so the mismatch was
           inert and invisible. Naming the ratio makes the intended
           quantity the one the signature accepts, and removes the
           division that was silently un-doing a unit error.

    Raises
    ------
    ValueError
        For curvilinear shapes when ``r0_over_R`` is missing, or is not a
        ratio in ``(0, 1)``. Use :func:`build_one_surface_compact_case`
        for solid geometry.
    """
    if shape == "slab":
        if _SLAB_VIA_UNIFIED:
            return _build_peierls_slab_case_via_unified(ng_key, n_regions)
        from .slab import _build_peierls_slab_case
        return _build_peierls_slab_case(
            ng_key, n_regions,
        )
    if shape in ("cylinder-1d", "sphere-1d"):
        if r0_over_R is None:
            raise ValueError(
                f"{shape} is a Class A (two-surface) case only when it has "
                f"a cavity: pass r0_over_R (the dimensionless ratio r_0/R). "
                f"Use build_one_surface_compact_case for the solid shape."
            )
        # Refuse AT THE BOUNDARY, before the O(minutes) mpmath solve.
        # ``reference_name`` enforces the same range, but it is called by the
        # stamp at the very END of the build — so without this the caller pays
        # a full solve on a geometry whose cavity is outside its own shell and
        # then gets an error about *naming*. Same law, right layer (Pattern 4).
        if not 0.0 < r0_over_R < 1.0:
            raise ValueError(
                f"{shape}: r0_over_R must be a RATIO in (0, 1), got "
                f"{r0_over_R!r}. The cavity lies strictly inside the outer "
                f"surface; a value at or above 1 is almost always an absolute "
                f"radius passed where the dimensionless r_0/R belongs (the "
                f"#345 unit error — see the parameter's note above)."
            )
        if shape == "cylinder-1d":
            from .cylinder import _build_peierls_cylinder_hollow_f4_case
            return _build_peierls_cylinder_hollow_f4_case(
                r0_over_R=float(r0_over_R), ng_key=ng_key,
            )
        from .sphere import _build_peierls_sphere_hollow_f4_case
        return _build_peierls_sphere_hollow_f4_case(
            r0_over_R=float(r0_over_R), ng_key=ng_key,
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


#: The shipped Class-A references, enumerated ONCE (#345).
#:
#: Three consumers used to re-spell this grid independently —
#: :func:`_class_a_cases` (eager build), :func:`continuous_case_builders`
#: (lazy keys) and :func:`capability_rows` (the published matrix) — under a
#: comment asking humans to keep them synchronised. They had already drifted:
#: the matrix computed its ``r_0`` tag as ``round(r0 * 100)`` where the other
#: two used ``round(r0/R_out * 100)``, agreeing only because every shipped
#: outer radius is ``1.0``. One grid makes that disagreement unspellable.
#:
#: The order here is the published row order of
#: ``_peierls_nystrom_capability_matrix.inc.rst`` — slab, then cylinder then
#: sphere, each ascending in :math:`r_0/R` with 1G before 2G. Reordering this
#: tuple reorders the rendered matrix, which the generator's ``--check`` mode
#: will catch.
#:
#: Each entry carries IDENTITY only. Per-row prose (accuracy class, closure
#: label, production status) legitimately differs per consumer and stays with
#: the consumer.
SHIPPED_CLASS_A: tuple[ShippedReference, ...] = (
    # Slab: 2G 2-region (current shipped default — native E₁ Nyström
    # path per the module docstring).
    ShippedReference("slab", n_groups=2, n_regions=2),
    # Hollow cylinder then sphere, F.4 at r_0/R ∈ {0.1, 0.2, 0.3}, 1G + 2G.
    *(
        ShippedReference(shape, n_groups=ng, n_regions=1, r0_over_R=r0)
        for shape in ("cylinder-1d", "sphere-1d")
        for r0 in (0.1, 0.2, 0.3)
        for ng in (1, 2)
    ),
)


def _build(case: ShippedReference) -> ContinuousReferenceSolution:
    """Materialise one grid entry. The single construction call."""
    return build_two_surface_case(
        case.shape, case.ng_key, case.n_regions, r0_over_R=case.r0_over_R,
    )


def _class_a_cases() -> list[ContinuousReferenceSolution]:
    """Class A — two-surface cases. Slab + hollow cylinder/sphere F.4.

    Enumerated by :data:`SHIPPED_CLASS_A`.

    Multi-group hollow cyl/sph references were added in Issue #104
    (2026-04-24) once the unified :func:`peierls_geometry.solve_peierls_mg`
    path landed. Each ``r_0/R`` sweep entry ships a 1G and 2G variant —
    the 1G residuals against :math:`k_\\infty` are reference-stable
    (1.4 % / 5.4 % / 13 % cyl; 0.4 % / 1.2 % / 3.3 % sph); the 2G variants
    inherit the same F.4 scalar rank-2 per-face closure applied group-wise.
    """
    return [_build(case) for case in SHIPPED_CLASS_A]


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
    Both the key here and that stamp come from
    :func:`~orpheus.derivations.continuous.peierls_nystrom.naming.reference_name`
    (#345) — before that they were two hand-written formulas, and the
    equivalence ``set(continuous_case_builders()) ==
    {c.name for c in continuous_cases()}`` pinned by
    ``tests/derivations/test_continuous_registry_lazy.py`` was carrying the
    risk. That test now guards the *enumeration* (does every grid entry build,
    and does the built object stamp the grid's name) rather than a race
    between two spellings.

    Enumerated by :data:`SHIPPED_CLASS_A` — the same grid :func:`_class_a_cases`
    and :func:`capability_rows` walk, so a new shipped reference is added in
    exactly one place.
    """
    from functools import partial

    return {
        case.name: partial(_build, case) for case in SHIPPED_CLASS_A
    }


# ---------------------------------------------------------------------
# Metadata-only enumeration (no eigenvalue solves) — the PROSE half of
# the Sphinx §theory-peierls-capabilities matrix. Its IDENTITY half is
# ``SHIPPED_CLASS_A`` above, walked in order, so a new shipped reference
# is added in exactly ONE place. (Until #345 this carried a second
# hand-written copy of the grid under a comment asking humans to keep
# the two synchronised; they had already drifted.) See
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
    tol_by_shape = {"cylinder-1d": _F4_CYL_TOL, "sphere-1d": _F4_SPH_TOL}
    cp_ref_by_shape = {"cylinder-1d": "cp_cylinder", "sphere-1d": "cp_sphere"}

    for case in SHIPPED_CLASS_A:
        row: dict[str, object] = {
            # IDENTITY — from the grid, so the published name is the
            # registry name by construction, not by agreement.
            "name": f"``{case.name}``",
            "geometry": case.shape,
            "n_groups": case.n_groups,
            "n_regions": case.n_regions,
            "r0_over_R": case.r0_over_R,
            "orbit_space_class": "A",
        }
        if case.r0_over_R is None:
            # Slab — rank-2 white closure at both parallel faces.
            row |= {
                "closure": rank2_label,
                "accuracy": "O(h²), Wigner-Seitz exact",
                "bc": bc_white,
                "status": shipped_status,
            }
        else:
            # Hollow curvilinear — F.4 closes the cavity face.
            tol_1g = tol_by_shape[case.shape][case.r0_over_R]
            row |= {"closure": f4_label, "bc": bc_f4}
            if case.n_groups == 1:
                row |= {
                    "accuracy": f"~{tol_1g} structural (scalar mode)",
                    "status": shipped_status,
                }
            else:
                row |= {
                    "accuracy": (
                        f"{case.n_groups}G builds, finite k_eff "
                        f"(``TestMG2GHollowRegistration``); k_eff vs "
                        f"``{cp_ref_by_shape[case.shape]}`` analytical not yet "
                        f"gated; structural residual expected ~{tol_1g} "
                        f"(group-local closure, unverified) — Issue #104 AC"
                    ),
                    "status": "shipped (k_eff gate pending — Issue #104)",
                }
        rows.append(row)

    # Class B — one-surface compact. No shipped references today
    # (rank-1 Mark floor is too loose; see Issues #101 / #103).

    return rows

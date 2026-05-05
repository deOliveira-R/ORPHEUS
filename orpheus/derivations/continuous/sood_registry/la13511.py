r"""Sood/Forster/Parsons LA-13511 (1999) benchmark case catalogue.

Production-protocol-aligned port of the legacy
``orpheus.derivations.continuous.fn_method.benchmarks.la13511`` module.
Each case now carries:

* **`materials: dict[int, Mixture]`** — keyed by integer material ID,
  exactly what :func:`orpheus.cp.solver.solve_cp` and
  :func:`orpheus.sn.solver.solve_sn` consume. Built via
  :func:`orpheus.derivations.common.xs_library.make_mixture` from the
  raw Sood XS components (ν, Σ_f, Σ_c, Σ_s, χ).
* **`geometry_spec: GeometrySpec`** — encodes geometry kind +
  critical-dimension; ``geometry_spec.build(n_cells)`` produces a
  :class:`orpheus.geometry.mesh.Mesh1D` at the desired refinement.
* **`truth: La13511Truth`** — all published reference values bundled
  in one struct. Different cases populate different subsets.

The legacy attributes (``case.sigma_t``, ``case.sigma_s`` etc.) are
retained as **read-only properties** that delegate to the materials
dict via :func:`.extractors.mixture_to_fn_arrays`. This keeps the F_N
test suite working unmodified while migration consumers flip to the
new shape one at a time.

ORPHEUS convention vs Sood convention
-------------------------------------

Sood et al. number energy groups :math:`g=N` (fast) → :math:`g=1`
(slow), the reverse of typical nuclear-engineering convention. This
module does the conversion at construction time; consumers see XS
arrays in **ORPHEUS convention** (``g=0`` fast, ``g=N-1`` slow).

The scattering matrix uses ORPHEUS's ``[from, to]`` convention:
``sigma_s[g, h]`` is :math:`\Sigma_{s, g \to h}`. The
``Sigma_g^{rem}`` removal cross sections used by Sood Eq 26-27 are
therefore ``sigma_t[g] - sigma_s[g, g]``.

Provenance
----------

* All XS values: LA-13511 Tables 1-67.
* All :math:`k_\infty`, :math:`r_c`, flux ratios: LA-13511 Section
  VII (1G), Section VII.B (2G), with primary references cited per case.

Python value precision: published values transcribed verbatim. Where
LA-13511 reports e.g. "k_inf = 2.612903" with 6 published digits, the
catalogue stores 2.612903 exactly as a Python float. Tolerance for
verification is 1e-5 (i.e., 5 significant figures match), so the
6th-digit rounding of the published value is not load-bearing.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry.mesh import BC


# ═══════════════════════════════════════════════════════════════════
# Geometry spec + truth dataclasses
# ═══════════════════════════════════════════════════════════════════
#
# :class:`GeometrySpec` was promoted to
# :mod:`orpheus.derivations.common.geometry_spec` on 2026-05-03 (R0.5)
# per the geometry-handling unification audit and renamed from
# ``MeshTemplate`` (originally named for its discrete-mesh-builder role)
# to ``GeometrySpec`` (descriptive of what it IS, not what it produces)
# on 2026-05-03. Imported directly from the canonical location below.


@dataclass(frozen=True)
class La13511Truth:
    """Tabulated reference values for a Sood LA-13511 case.

    Different cases populate different subsets. ``None`` means the
    published reference does not tabulate that quantity.

    Parameters
    ----------
    k_eff_or_kinf : float
        Reference :math:`k_{\\rm eff}` (finite cases — usually 1.0,
        critical) or :math:`k_\\infty` (infinite cases).
    flux_ratios : Mapping[float, float] | None
        For 1G cases with published flux table: dict mapping
        ``r/r_c`` to :math:`\\phi(r)/\\phi(0)`.
    flux_ratio_groupwise : Mapping[int, float] | None
        For 2G infinite cases: dict mapping ORPHEUS group index to
        :math:`\\phi_g/\\phi_0` (ratio relative to ORPHEUS group 0 =
        fast). For ``PU-2-0-IN``, this is
        ``{0: 1.0, 1: phi_slow/phi_fast}``.
    angular_flux_at_surface : Mapping[float, Mapping[float, float]] | None
        Reserved for future cases that publish surface angular flux
        :math:`\\psi(\\mu, r=R)`. ``None`` for first-slice cases.
    critical_dimension_mfp : float | None
        Published critical dimension in mean free paths.

        For ``"slab"``: the half-thickness :math:`a` (F_N convention).
        For ``"sphere"`` / ``"cylinder"``: the radius :math:`R`.
        For ``"infinite"``: ``None``.

        Use case: registry-truth value (the published critical
        configuration). Multiply by :math:`1 / \\Sigma_t` to convert to
        cm; this is what :meth:`La13511Case.to_geometry` does internally.
        Living on Truth and not on geometry mirrors the architectural
        fact that this is a truth claim ("at this size, the configuration
        is critical"), not a geometric description.
    extrapolated_endpoint_mfp : float | None
        Published extrapolated endpoint :math:`z_0` in mean free paths.

        Standard transport-theory value used for diffusion-theory
        boundary conditions. Optional metadata; not all cases publish it.
    """

    k_eff_or_kinf: float
    flux_ratios: Mapping[float, float] | None = None
    flux_ratio_groupwise: Mapping[int, float] | None = None
    angular_flux_at_surface: Mapping[float, Mapping[float, float]] | None = None
    critical_dimension_mfp: float | None = None
    extrapolated_endpoint_mfp: float | None = None


# ═══════════════════════════════════════════════════════════════════
# Provenance — citation / publication metadata
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class Provenance:
    """Citation / publication metadata for a registry case.

    Separates *where the case came from* from *what the case is*.
    Pure metadata; not used in any solve. Mirrors (with structure) the
    legacy flat fields ``La13511Case.sood_table`` /
    ``La13511Case.primary_reference`` / ``La13511Case.notes``.

    Parameters
    ----------
    paper_id : str
        Short identifier for the source paper (e.g. ``"LA-13511"``,
        ``"Atalay-1997"``, ``"NM-1980"``).
    paper_table : int | str
        Table number (or compound table/case label) within the paper
        where this case is tabulated. Integer for simple cases (e.g.
        ``8``); string for compound labels (e.g. ``"Table 1 / Case 4"``).
    primary_reference : str
        Full bibliographic citation for the *primary* peer-reviewed
        source — the paper Sood (or the case's own author) cites as
        the source of the reference values.
    notes : str
        Any extra remarks (typo flags, convention drift, conversion
        subtleties, known issues, cross-checks). Default empty string.
    """

    paper_id: str
    paper_table: int | str
    primary_reference: str
    notes: str = ""


# ═══════════════════════════════════════════════════════════════════
# Main case dataclass
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class La13511Case:
    """A single LA-13511 benchmark configuration.

    Production-protocol form: cross sections live in a
    :class:`Mixture` (the same object production solvers consume),
    geometry lives in a :class:`GeometrySpec` (which builds a
    :class:`Mesh1D` on demand), reference values live in a
    :class:`La13511Truth`.

    Parameters
    ----------
    case_id : str
        Sood naming convention identifier
        (``<Material>-<Groups>-<Scattering>-<Geometry>``).
    problem_number : int
        Sequential problem number within LA-13511.
    description : str
        One-line human-readable description.
    materials : dict[int, Mixture]
        Macroscopic cross sections keyed by material ID. Single-region
        cases use ``{0: Mixture(...)}``; multi-region cases (none in
        the first slice) add more keys.
    geometry_spec : GeometrySpec
        Geometry specification + critical dimension. Use
        :meth:`GeometrySpec.build` to obtain a concrete mesh.
    scattering_order : int
        Legendre order of the scattering kernel (0 = isotropic,
        1 = P_1, 2 = P_2). All first-slice cases are isotropic
        (``= 0``).
    truth : La13511Truth
        All published reference values for this case.
    sood_table : int
        LA-13511 table number where this case is tabulated.
    primary_reference : str
        The peer-reviewed paper Sood cites as the source of the
        reference values.
    notes : str
        Any extra remarks (typo flags, conversion subtleties, etc).
    provenance : Provenance | None
        Structured citation metadata mirroring ``sood_table`` /
        ``primary_reference`` / ``notes``. Optional today (Phase B
        adds the field; Phase C/D migrate consumers; Phase F removes
        the flat fields). When None, no structured provenance has
        been populated yet.

    Notes
    -----
    The legacy F_N consumer interface (``case.sigma_t``,
    ``case.sigma_s``, ``case.nu_sigma_f``, ``case.chi``,
    ``case.geometry``, ``case.n_groups``,
    ``case.critical_dimension_mfp``, ``case.critical_dimension_cm``,
    ``case.k_eff_or_kinf``, ``case.flux_ratios``,
    ``case.flux_ratio_groupwise``) is exposed as **read-only
    properties** that delegate to ``materials[0]`` /
    ``geometry_spec`` / ``truth``. This keeps the F_N test suite
    working unchanged.
    """

    case_id: str
    problem_number: int
    description: str
    materials: dict[int, Mixture]
    geometry_spec: GeometrySpec
    scattering_order: int
    truth: La13511Truth
    sood_table: int
    primary_reference: str
    notes: str = ""
    provenance: Provenance | None = None

    # ── Legacy compatibility properties ────────────────────────────

    @property
    def n_groups(self) -> int:
        """Number of energy groups."""
        return self.geometry_spec.n_groups

    @property
    def geometry(self) -> str:
        """Geometry kind (``"infinite"``, ``"slab"``, ...)."""
        return self.geometry_spec.geometry

    @property
    def critical_dimension_mfp(self) -> float | None:
        """Critical dimension in mean free paths (or None for infinite)."""
        return self.geometry_spec.critical_dimension_mfp

    @property
    def critical_dimension_cm(self) -> float | None:
        """Critical dimension in cm (or None for infinite)."""
        return self.geometry_spec.critical_dimension_cm

    @property
    def k_eff_or_kinf(self) -> float:
        """Reference :math:`k_{\\rm eff}` (finite) or :math:`k_\\infty` (infinite)."""
        return self.truth.k_eff_or_kinf

    @property
    def flux_ratios(self) -> Mapping[float, float] | None:
        """Published flux ratio table, or None."""
        return self.truth.flux_ratios

    @property
    def flux_ratio_groupwise(self) -> Mapping[int, float] | None:
        """Published 2G groupwise flux ratio, or None."""
        return self.truth.flux_ratio_groupwise

    @property
    def sigma_t(self) -> np.ndarray:
        """Total XS in ORPHEUS convention, shape ``(n_groups,)``.

        Convenience accessor for legacy F_N consumers; pulls from
        ``self.materials[0]``.
        """
        from .extractors import mixture_to_fn_arrays
        sigma_t, _, _, _ = mixture_to_fn_arrays(self._primary_mixture)
        return sigma_t

    @property
    def sigma_s(self) -> np.ndarray:
        """P_0 scattering matrix in ``[from, to]`` convention, shape ``(n_groups, n_groups)``."""
        from .extractors import mixture_to_fn_arrays
        _, sigma_s, _, _ = mixture_to_fn_arrays(self._primary_mixture)
        return sigma_s

    @property
    def nu_sigma_f(self) -> np.ndarray:
        """Production XS :math:`\\nu \\Sigma_f`, shape ``(n_groups,)``."""
        from .extractors import mixture_to_fn_arrays
        _, _, nu_sigma_f, _ = mixture_to_fn_arrays(self._primary_mixture)
        return nu_sigma_f

    @property
    def chi(self) -> np.ndarray:
        """Fission spectrum, shape ``(n_groups,)``."""
        from .extractors import mixture_to_fn_arrays
        _, _, _, chi = mixture_to_fn_arrays(self._primary_mixture)
        return chi

    @property
    def _primary_mixture(self) -> Mixture:
        """The single mixture for first-slice (homogeneous) cases."""
        if len(self.materials) != 1:
            raise ValueError(
                f"Legacy property accessor expects single-region case; "
                f"case {self.case_id} has {len(self.materials)} materials. "
                f"Use case.materials[mat_id] explicitly."
            )
        return next(iter(self.materials.values()))


# ═══════════════════════════════════════════════════════════════════
# Helper: build a single-material 1G mixture from raw Sood XS
# ═══════════════════════════════════════════════════════════════════


def _mix_1g_isotropic(
    sigma_t: float,
    sigma_c: float,
    sigma_f: float,
    nu: float,
    sigma_s_self: float,
) -> Mixture:
    """Build a 1G isotropic Mixture from Sood's raw per-isotope XS components.

    Sood publishes :math:`(\\nu, \\Sigma_f, \\Sigma_c, \\Sigma_s, \\Sigma_t)`
    separately. Production solvers want the production XS
    :math:`\\nu \\Sigma_f` and the absorption XS components separately
    — :func:`make_mixture` accepts exactly that decomposition.
    """
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sigma_c]),
        sig_f=np.array([sigma_f]),
        nu=np.array([nu]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s_self]]),
    )


# ═══════════════════════════════════════════════════════════════════
# Case 1 — PUa-1-0-IN (Sood problem 1): 1G infinite medium, Pu-239 (a)
# ═══════════════════════════════════════════════════════════════════
#
# LA-13511 Table 2 (Pu-239 (a) cross sections, 1G isotropic):
#   ν = 3.24,  Σ_f = 0.0816,  Σ_c = 0.019584,  Σ_s = 0.225216,
#   Σ_t = 0.32640,  c = (Σ_s + νΣ_f)/Σ_t = 1.50.
# Reference: k_inf = 2.612903 (LA-13511 Eq 20).

PUA_1_0_IN = La13511Case(
    case_id="PUa-1-0-IN",
    problem_number=1,
    description="Pu-239 (a) bare infinite medium, 1G isotropic",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=3.24,
        sigma_s_self=0.225216,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.612903,
        flux_ratios=None,
    ),
    sood_table=2,
    primary_reference="LA-13511 Eq 20 (closed form)",
    notes=(
        "1G infinite-medium k_inf reduces to nu·Sigma_f/Sigma_a; the "
        "'c' factor in Eq 20 cancels algebraically (verified in "
        "fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_20_simplifies_to_eq_19)."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Case 5 — PU-2-0-IN (Sood problem 44): 2G infinite medium, Pu
# ═══════════════════════════════════════════════════════════════════
#
# LA-13511 Tables 30-31 (Pu 2G isotropic, no upscatter):
#
#   Sood convention (g=2 fast, g=1 slow):
#     Fast (g=2):  ν_2=3.10,  Σ_2f=0.0936,  Σ_2c=0.00480,
#                  Σ_22s=0.0792, Σ_12s=0.0432 (from g=2 to g=1),
#                  Σ_2 = 0.2208,  χ_2 = 0.575.
#     Slow (g=1):  ν_1=2.93,  Σ_1f=0.08544, Σ_1c=0.0144,
#                  Σ_11s=0.23616, Σ_21s=0.0 (no upscatter),
#                  Σ_1 = 0.3360,  χ_1 = 0.425.
#
# In ORPHEUS convention (g=0 fast, g=1 slow):
#     sigma_t = [0.2208, 0.3360]
#     sigma_s[g=0=fast, g=0=fast] = Σ_22s = 0.0792
#     sigma_s[g=0=fast, g=1=slow] = (from fast to slow) = Σ_12s = 0.0432
#     sigma_s[g=1=slow, g=0=fast] = (from slow to fast) = Σ_21s = 0.0
#     sigma_s[g=1=slow, g=1=slow] = Σ_11s = 0.23616
#     sigma_f = [0.0936, 0.08544]
#     sigma_c = [0.00480, 0.0144]
#     nu = [3.10, 2.93]
#     chi = [0.575, 0.425]
#
# Reference: k_inf = 2.683767, φ_2/φ_1 = 0.675229 (Sood Eq 28-29 + Eq 32).

PU_2_0_IN = La13511Case(
    case_id="PU-2-0-IN",
    problem_number=44,
    description="Pu bare infinite medium, 2G isotropic, no upscatter",
    materials={0: make_mixture(
        sig_t=np.array([0.2208, 0.3360]),
        sig_c=np.array([0.00480, 0.0144]),
        sig_f=np.array([0.0936, 0.08544]),
        nu=np.array([3.10, 2.93]),
        chi=np.array([0.575, 0.425]),
        sig_s=np.array([
            [0.0792,  0.0432],   # from g=0 (fast):  → fast self, → slow downscatter
            [0.0,     0.23616],  # from g=1 (slow):  → fast upscatter (none), → slow self
        ]),
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.683767,
        flux_ratios=None,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 0.675229},  # slow/fast = 1/Sood ratio
    ),
    sood_table=30,
    primary_reference="LA-13511 Eq 28-29 (k_inf) + Eq 32 (flux ratio)",
    notes=(
        "Sood Eq 28 has a typo: the chi_1 and chi_2 numerators have "
        "the wrong Sigma_g^rem factor — Eq 28 as printed reduces to "
        "2.862 (not 2.684) when Sigma_21s = 0, while Eq 29 (printed "
        "separately) reduces correctly. The SymPy derivation in "
        "fn_method.origins.k_inf_derivations.derive_kinf_2g_general_from_matrix "
        "exposes the correct general form by computing det(M)=0 from "
        "Eq 25 directly, and verifies the published Eq 29 against it."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Cases 2-4 — bare critical sphere/slab/cylinder of U-235 (a), 1G
# ═══════════════════════════════════════════════════════════════════
#
# Common XS (LA-13511 Tables 4-6, U-235 (a) 1G isotropic):
#   ν = 2.70,  Σ_f = 0.06528,  Σ_c = 0.013056,  Σ_s = 0.248064,
#   Σ_t = 0.32640,  c = (Σ_s + νΣ_f)/Σ_t = 1.30.
# All three cases are critical (k_eff = 1.0 by construction).

_UA_1G_KW = dict(
    sigma_t=0.32640,
    sigma_c=0.013056,
    sigma_f=0.06528,
    nu=2.70,
    sigma_s_self=0.248064,
)


# Case 2 — Ua-1-0-SL (Sood problem 12): 1G bare slab, U-235 (a)

UA_1_0_SL_STUB = La13511Case(
    case_id="Ua-1-0-SL",
    problem_number=12,
    description="U-235 (a) bare slab, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=2.872934,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.25: 0.9669506,
            0.50: 0.8686259,
            0.75: 0.7055218,
            1.00: 0.4461912,
        },
    ),
    sood_table=4,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94",
    notes=(
        "Slab F_N solver shipped at ≤5e-6 absolute on a_c (see "
        "fn_method.slab.solve_fn_slab_bare_critical). Critical "
        "dimension is the half-thickness; slab full width is "
        "2*critical_dimension_cm = 5.745868 cm."
    ),
)


# Case 3 — Ua-1-0-CY (Sood problem 13): 1G bare cylinder, U-235 (a)

UA_1_0_CY_STUB = La13511Case(
    case_id="Ua-1-0-CY",
    problem_number=13,
    description="U-235 (a) bare infinite cylinder, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_spec=GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=1.72500292,
        critical_dimension_cm=5.284935,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
    ),
    sood_table=5,
    primary_reference="Westfall-Metcalf 1973 NSE 52, 1",
    notes=(
        "WM-72 singular-eigenfunction cylinder solver shipped at ~1% "
        "relative accuracy (single-cell product integration on the "
        "log-singular kernel diagonal; see "
        "orpheus.derivations.continuous.singular_eigenfunction.cylinder). "
        "Variant α cylinder cross-check at 8.5e-6 holds the strict "
        "1e-5 anchor. WM-72 prototype provides the second, "
        "structurally-independent cross-check anchor (different "
        "mathematical pillar than Variant α / Bickley-Naylor)."
    ),
)


# Case 4 — Ua-1-0-SP (Sood problem 14): 1G bare sphere, U-235 (a)

UA_1_0_SP_STUB = La13511Case(
    case_id="Ua-1-0-SP",
    problem_number=14,
    description="U-235 (a) bare sphere, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.4248249802,
        critical_dimension_cm=7.428998,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            # Kaper-Lindeman-Leaf 1974 Table VII at c=1.30 — the same
            # XS as Sood Ua-1-0-SP (U-235 (a) 1G isotropic, c=1.30).
            0.25: 0.93244907,
            0.50: 0.74553332,
            0.75: 0.48095413,
            1.00: 0.17177706,
        },
    ),
    sood_table=6,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94 (Table VII)",
    notes=(
        "Sphere F_N solver shipped at ≤1e-7 absolute on R_c (see "
        "fn_method.sphere.solve_fn_sphere_bare_critical). Used as the "
        "structurally-independent L1 reference for Variant α sphere. "
        "Flux ratios populated from KLL Table VII c=1.30 row — same "
        "XS as this case (cross-check via "
        "fn_method.sphere.flux_reconstruction)."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Phase B3 — wide enumeration of currently-implementable cases
# ═══════════════════════════════════════════════════════════════════
#
# Cases below are added by Phase B3 of the wide-enumeration sweep
# (closeout: ``.claude/agent-memory/method-implementer/
# sood_registry_wide_enumeration_phase_b3.md``). They cover the subset
# of LA-13511's 75 problems that the existing ``fn_method`` machinery
# (k_inf 1G/2G/MG, slab F_N, sphere F_N) can solve TODAY.
#
# Coverage matrix:
#   * 1G k_inf infinite (PU/U-235/UD2O/URR variants):  9 cases
#   * 2G k_inf infinite (PU/U/UAL/URR/UD2O):           7 cases
#   * 3G k_inf infinite (URR-3-0-IN):                  1 case
#   * 6G k_inf infinite (URR-6-0-IN):                  1 case
#   * Slab F_N bare-critical (1G isotropic):           3 new cases
#   * Sphere F_N bare-critical (1G isotropic):         3 new cases
#   * Cylinder bare-critical (1G isotropic):           2 new STUBS
#                                                      (B1 dispatch will activate)
#   * 2G bare-critical slab/sphere:                    7 STUBS
#                                                      (machinery deferred —
#                                                      needs Siewert-Thomas
#                                                      1986 2G F_N)
# ═══════════════════════════════════════════════════════════════════


# ───────────────────────────────────────────────────────────────────
# Helper: 2-group ORPHEUS-ordered Mixture from Sood-side fast/slow XS
# ───────────────────────────────────────────────────────────────────


def _mix_2g_isotropic(
    *,
    sigma_t_fast: float,
    sigma_t_slow: float,
    sigma_c_fast: float,
    sigma_c_slow: float,
    sigma_f_fast: float,
    sigma_f_slow: float,
    nu_fast: float,
    nu_slow: float,
    chi_fast: float,
    chi_slow: float,
    sigma_22s: float,   # Sood Σ_22s = self-scatter in fast (g=2)
    sigma_11s: float,   # Sood Σ_11s = self-scatter in slow (g=1)
    sigma_12s: float,   # Sood Σ_12s = scatter from g=2 (fast) to g=1 (slow)
    sigma_21s: float,   # Sood Σ_21s = scatter from g=1 (slow) to g=2 (fast) — upscatter
) -> Mixture:
    """Build a 2G isotropic Mixture from Sood-convention raw XS.

    Sood numbers groups g=2 fast → g=1 slow. ORPHEUS uses g=0 fast,
    g=N-1 slow. This helper keeps the call site readable in
    Sood-convention names while emitting an ORPHEUS-ordered Mixture.

    The scattering matrix in ORPHEUS ``[from, to]`` convention:
    ``sigma_s[0, 0] = Σ_22s`` (fast self), ``sigma_s[1, 1] = Σ_11s``
    (slow self), ``sigma_s[0, 1] = Σ_12s`` (downscatter fast→slow),
    ``sigma_s[1, 0] = Σ_21s`` (upscatter slow→fast — zero for most
    cases, nonzero for URRb/URRc/URRd).
    """
    return make_mixture(
        sig_t=np.array([sigma_t_fast, sigma_t_slow]),
        sig_c=np.array([sigma_c_fast, sigma_c_slow]),
        sig_f=np.array([sigma_f_fast, sigma_f_slow]),
        nu=np.array([nu_fast, nu_slow]),
        chi=np.array([chi_fast, chi_slow]),
        sig_s=np.array([
            [sigma_22s, sigma_12s],   # from fast: self, downscatter
            [sigma_21s, sigma_11s],   # from slow: upscatter, self
        ]),
    )


# ═══════════════════════════════════════════════════════════════════
# 1G k_inf infinite-medium cases (LA-13511 Tables 5, 12, 16, 20)
# ═══════════════════════════════════════════════════════════════════
#
# Pure rational algebra in ν, Σ_f, Σ_s, Σ_t (Sood Eq 19/20). All
# verifiable to machine precision; tolerance is set by Sood's published
# precision (≤ 1e-5 absolute on a 6-7 digit truth value).

# Case 5 — PUb-1-0-IN: Pu-239 (b), c=1.40
PUB_1_0_IN = La13511Case(
    case_id="PUb-1-0-IN",
    problem_number=5,
    description="Pu-239 (b) bare infinite medium, 1G isotropic, c=1.40",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=2.84,
        sigma_s_self=0.225216,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.290323),
    sood_table=5,
    primary_reference="LA-13511 Eq 19 / Table 5",
    notes="Same Σ_t / Σ_s as PUa, only ν changes (3.24 → 2.84). c=1.40.",
)


# Case 11 — Ua-1-0-IN: U-235 (a), c=1.30
UA_1_0_IN = La13511Case(
    case_id="Ua-1-0-IN",
    problem_number=11,
    description="U-235 (a) bare infinite medium, 1G isotropic, c=1.30",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.25),
    sood_table=12,
    primary_reference="LA-13511 Eq 19 / Table 12",
    notes="Sood publishes 'k_inf = 2.25' (3 digits printed but algebraically exact).",
)


# Case 15 — Ub-1-0-IN: U-235 (b), c=1.3194202
UB_1_0_IN = La13511Case(
    case_id="Ub-1-0-IN",
    problem_number=15,
    description="U-235 (b) bare infinite medium, 1G isotropic, c=1.3194202",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.013056,
        sigma_f=0.065280,
        nu=2.797101,
        sigma_s_self=0.248064,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.330917),
    sood_table=12,
    primary_reference="LA-13511 Eq 19 / Table 12",
    notes="Cross-section variant (b): same Σ_t/Σ_s as Ua, ν tuned to give c=1.3194202.",
)


# Case 17 — Uc-1-0-IN: U-235 (c), c=1.3014616
UC_1_0_IN = La13511Case(
    case_id="Uc-1-0-IN",
    problem_number=17,
    description="U-235 (c) bare infinite medium, 1G isotropic, c=1.3014616",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.013056,
        sigma_f=0.065280,
        nu=2.707308,
        sigma_s_self=0.248064,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.256083),
    sood_table=12,
    primary_reference="LA-13511 Eq 19 / Table 12",
    notes="Cross-section variant (c): same Σ_t/Σ_s as Ua, ν=2.707308.",
)


# Case 19 — Ud-1-0-IN: U-235 (d), c=1.2958396
UD_1_0_IN = La13511Case(
    case_id="Ud-1-0-IN",
    problem_number=19,
    description="U-235 (d) bare infinite medium, 1G isotropic, c=1.2958396",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.013056,
        sigma_f=0.065280,
        nu=2.679198,
        sigma_s_self=0.248064,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.232667),
    sood_table=12,
    primary_reference="LA-13511 Eq 19 / Table 12",
    notes="Cross-section variant (d): same Σ_t/Σ_s as Ua, ν=2.679198.",
)


# Case 21 — UD2O-1-0-IN: U-D2O reactor, c=1.02
UD2O_1_0_IN = La13511Case(
    case_id="UD2O-1-0-IN",
    problem_number=21,
    description="U-D2O reactor bare infinite medium, 1G isotropic, c=1.02",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.70,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.133333),
    sood_table=16,
    primary_reference="LA-13511 Eq 19 / Table 16",
    notes="Heavy-water-moderated low-enrichment U; lowest c in the bare 1G family (1.02).",
)


# Case 29 — Ue-1-0-IN: U-235 reactor with Fe/Na surrounds (infinite-medium k_inf only — no Fe/Na in this case)
UE_1_0_IN = La13511Case(
    case_id="Ue-1-0-IN",
    problem_number=29,
    description="U-235 reactor (e) bare infinite medium, 1G isotropic, c=1.230",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.407407,
        sigma_c=0.01013756,
        sigma_f=0.06922744,
        nu=2.50,
        sigma_s_self=0.328042,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.1806667),
    sood_table=20,
    primary_reference="LA-13511 Eq 19 / Table 20",
    notes=(
        "U-235 (e) cross sections used in the Ue-Fe-Na multi-region "
        "case; the infinite-medium variant uses the U-235 (e) XS alone."
    ),
)


# Case 31 — PU-1-1-IN: Pu-239 with linearly anisotropic scattering — k_inf is identical to the isotropic case
PU_1_1_IN = La13511Case(
    case_id="PU-1-1-IN",
    problem_number=31,
    description="Pu-239 (a/b) bare infinite medium, 1G P_1 anisotropic, c=1.40",
    materials={0: _mix_1g_isotropic(
        sigma_t=1.0,
        sigma_c=0.0,
        sigma_f=0.266667,
        nu=2.5,
        sigma_s_self=0.733333,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,  # k_inf does not depend on anisotropy; isotropic XS suffice
    truth=La13511Truth(k_eff_or_kinf=2.5),
    sood_table=24,
    primary_reference="LA-13511 Eq 19 / Table 24",
    notes=(
        "Sood: 'The anisotropic scattering cross sections do not change "
        "k_inf' (LA-13511 p. 22). Catalogued with scattering_order=0 + "
        "Σ_s = Σ_s0 since the P_1 moment is a no-op for infinite-medium "
        "k_inf. Anisotropic data lives in the slab cases 32-35."
    ),
)


# Case 38, 40, 42 — UD2O-{a,b,c}-1-1-IN: U-D2O P_1 anisotropic — k_inf depends only on isotropic XS
# Each (a,b,c) has slightly different ν tuned to give specified c values 1.0308381 / 1.0341086 / 1.01964.
UD2OA_1_1_IN = La13511Case(
    case_id="UD2Oa-1-1-IN",
    problem_number=38,
    description="U-D2O (a) bare infinite medium, 1G P_1 anisotropic, c=1.0308381",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.808381,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.205587),
    sood_table=28,
    primary_reference="LA-13511 Eq 19 / Table 28",
    notes="P_1 anisotropy doesn't change k_inf; ν tuned for c=1.0308381.",
)

UD2OB_1_1_IN = La13511Case(
    case_id="UD2Ob-1-1-IN",
    problem_number=40,
    description="U-D2O (b) bare infinite medium, 1G P_1 anisotropic, c=1.0341086",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.841086,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.227391),
    sood_table=28,
    primary_reference="LA-13511 Eq 19 / Table 28",
    notes="P_1 anisotropy doesn't change k_inf; ν tuned for c=1.0341086.",
)

UD2OC_1_1_IN = La13511Case(
    case_id="UD2Oc-1-1-IN",
    problem_number=42,
    description="U-D2O (c) bare infinite medium, 1G P_1 anisotropic, c=1.01964",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.6964,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.130933),
    sood_table=28,
    primary_reference="LA-13511 Eq 19 / Table 28",
    notes=(
        "U-D2O (c) has *negative* P_1 scattering moment Σ_s1 = -0.27850447 "
        "(backward-peaked); k_inf still depends only on Σ_s0. Slab "
        "cases (39/41/43) inherit anisotropy."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# 2G k_inf infinite-medium cases (LA-13511 Tables 30-31, 33-34, 36-37,
# 39-40, 43-44, 46-47, 49-50)
# ═══════════════════════════════════════════════════════════════════


# Case 47 — U-2-0-IN: U-235, 2G isotropic, no upscatter
U_2_0_IN = La13511Case(
    case_id="U-2-0-IN",
    problem_number=47,
    description="U-235 bare infinite medium, 2G isotropic, no upscatter",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.2160, sigma_t_slow=0.3456,
        sigma_c_fast=0.00384, sigma_c_slow=0.01344,
        sigma_f_fast=0.06192, sigma_f_slow=0.06912,
        nu_fast=2.70, nu_slow=2.50,
        chi_fast=0.575, chi_slow=0.425,
        sigma_22s=0.078240, sigma_11s=0.26304,
        sigma_12s=0.0720, sigma_21s=0.0,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.216349,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 0.474967},
    ),
    sood_table=33,
    primary_reference="LA-13511 Eq 29 (k_inf) + Eq 32 (flux ratio) / Tables 33-34",
    notes="Sood publishes φ_2/φ_1 (fast/slow) = 0.474967.",
)


# Case 50 — UAL-2-0-IN: Uranium-Aluminum-Water assembly, 2G isotropic, no upscatter
# Note: Σ_2f = 0.0 in fast group — fission only in slow group. χ_fast = 1.0, χ_slow = 0.0.
UAL_2_0_IN = La13511Case(
    case_id="UAL-2-0-IN",
    problem_number=50,
    description="U-Al-Water assembly bare infinite medium, 2G isotropic",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.26817, sigma_t_slow=1.27698,
        sigma_c_fast=0.000222, sigma_c_slow=0.00314363958,
        sigma_f_fast=0.0, sigma_f_slow=0.06070636042,
        nu_fast=0.0, nu_slow=2.83,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.247516, sigma_11s=1.21313,
        sigma_12s=0.020432, sigma_21s=0.0,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.661745,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 3.1250},
    ),
    sood_table=36,
    primary_reference="LA-13511 Eq 29 + Eq 32 / Tables 36-37",
    notes=(
        "Slow-only fission (ν_2 = 0, χ_1 = 0). Sood publishes "
        "φ_2/φ_1 = 3.1250 — i.e. fast/slow > 1 because slow group is "
        "very absorbing (Σ_1 = 1.27698 mostly self-scatter)."
    ),
)


# Case 53 — URRa-2-0-IN: 93%-enriched U research reactor (a), 2G isotropic, no upscatter
URRA_2_0_IN = La13511Case(
    case_id="URRa-2-0-IN",
    problem_number=53,
    description="URR (a) — 93% enriched U bare infinite medium, 2G isotropic",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.65696, sigma_t_slow=2.52025,
        sigma_c_fast=0.0010046, sigma_c_slow=0.025788,
        sigma_f_fast=0.0010484, sigma_f_slow=0.050632,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.62568, sigma_11s=2.44383,
        sigma_12s=0.029227, sigma_21s=0.0,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.631452,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 2.614706},
    ),
    sood_table=39,
    primary_reference="LA-13511 Eq 29 + Eq 32 / Tables 39-40",
    notes="93%-enriched bare research reactor. χ_1=0 (all fission from fast group).",
)


# Case 56 — URRb-2-0-IN: research reactor (b) WITH thermal upscatter (Σ_21s = 0.000767)
URRB_2_0_IN = La13511Case(
    case_id="URRb-2-0-IN",
    problem_number=56,
    description="URR (b) bare infinite medium, 2G isotropic, *with* thermal upscatter",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.88721, sigma_t_slow=2.9727,
        sigma_c_fast=0.001104, sigma_c_slow=0.024069,
        sigma_f_fast=0.000836, sigma_f_slow=0.029564,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.83892, sigma_11s=2.9183,
        sigma_12s=0.04635, sigma_21s=0.000767,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.365821,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 1.173679},
    ),
    sood_table=43,
    primary_reference="LA-13511 Eq 28 (general — with upscatter) / Tables 43-44",
    notes=(
        "Has thermal upscatter (Σ_21s = 0.000767). MUST use the general "
        "Eq-28 formula (compute_kinf_2g_general or compute_kinf_mg), NOT "
        "the no-upscatter Eq-29 specialisation."
    ),
)


# Case 57 — URRc-2-0-IN: research reactor (c) WITH thermal upscatter (Σ_21s = 0.00116)
URRC_2_0_IN = La13511Case(
    case_id="URRc-2-0-IN",
    problem_number=57,
    description="URR (c) bare infinite medium, 2G isotropic, *with* thermal upscatter",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.88655, sigma_t_slow=2.9628,
        sigma_c_fast=0.001472, sigma_c_slow=0.029244,
        sigma_f_fast=0.001648, sigma_f_slow=0.057296,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.83807, sigma_11s=2.8751,
        sigma_12s=0.04536, sigma_21s=0.00116,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.633380,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 1.933422},
    ),
    sood_table=43,
    primary_reference="LA-13511 Eq 28 (general — with upscatter) / Tables 43-44",
    notes="Larger Σ_1f / Σ_21s than URRb; same upscatter structure.",
)


# Case 62 — URRd-2-0-IN: ISLC base material, 2G isotropic, no upscatter, ν=1.004 (slightly unphysical per Sood)
URRD_2_0_IN = La13511Case(
    case_id="URRd-2-0-IN",
    problem_number=62,
    description="URR (d) bare infinite medium, 2G isotropic — ISLC base material",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.650917, sigma_t_slow=2.13800,
        sigma_c_fast=0.0019662, sigma_c_slow=0.023496,
        sigma_f_fast=0.61475, sigma_f_slow=0.045704,
        nu_fast=1.004, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.0, sigma_11s=2.06880,
        sigma_12s=0.0342008, sigma_21s=0.0,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.034970,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 2.023344},
    ),
    sood_table=46,
    primary_reference="LA-13511 Eq 29 + Eq 32 / Tables 46-47",
    notes=(
        "ISLC (Infinite Slab Lattice Cell) base XS. Sood uses ν_fast=1.004 "
        "'to stress code verification' (LA-13511 p. 25 — i.e. unphysical "
        "but algebraically valid). Σ_22s = 0 (no fast self-scatter)."
    ),
)


# Case 67 — UD2O-2-0-IN: U-D2O reactor, 2G isotropic, no upscatter — k_inf is just barely critical (1.000196)
UD2O_2_0_IN = La13511Case(
    case_id="UD2O-2-0-IN",
    problem_number=67,
    description="U-D2O reactor bare infinite medium, 2G isotropic",
    materials={0: _mix_2g_isotropic(
        sigma_t_fast=0.33588, sigma_t_slow=0.54628,
        sigma_c_fast=0.008708, sigma_c_slow=0.02518,
        sigma_f_fast=0.002817, sigma_f_slow=0.097,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.31980, sigma_11s=0.42410,
        sigma_12s=0.004555, sigma_21s=0.0,
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.000196,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 26.823271},
    ),
    sood_table=49,
    primary_reference="LA-13511 Eq 29 + Eq 32 / Tables 49-50",
    notes=(
        "Heavy-water reactor; k_inf = 1.000196 is essentially at the "
        "infinite-medium critical threshold. φ_fast/φ_slow = 26.82 "
        "(very slow-flux dominated due to D2O moderation)."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# 3G k_inf infinite-medium case (LA-13511 Tables 59-61)
# ═══════════════════════════════════════════════════════════════════


# Case 74 — URR-3-0-IN: 3-group URR, no upscatter
# Sood ordering: g3=fast, g2=mid, g1=slow. ORPHEUS: g=0 fast, g=1 mid, g=2 slow.
URR_3_0_IN = La13511Case(
    case_id="URR-3-0-IN",
    problem_number=74,
    description="URR bare infinite medium, 3G isotropic, no upscatter",
    materials={0: make_mixture(
        sig_t=np.array([0.240, 0.975, 3.10]),
        sig_c=np.array([0.006, 0.040, 0.20]),
        sig_f=np.array([0.006, 0.060, 0.90]),
        nu=np.array([3.0, 2.5, 2.0]),
        chi=np.array([0.96, 0.04, 0.0]),
        # ORPHEUS [from, to] convention.
        # Sood Σ_{i,j,s} means scatter from g=j (Sood) to g=i (Sood).
        # ORPHEUS index: 0=fast=Sood 3, 1=mid=Sood 2, 2=slow=Sood 1.
        sig_s=np.array([
            [0.024, 0.171, 0.033],   # from fast: Σ_33s, Σ_23s, Σ_13s
            [0.0,   0.60,  0.275],   # from mid:  Σ_32s=0, Σ_22s, Σ_12s
            [0.0,   0.0,   2.0  ],   # from slow: no upscatter, self
        ]),
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=3,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.60,
        flux_ratio_groupwise={0: 1.0, 1: 0.480, 2: 0.150},
    ),
    sood_table=59,
    primary_reference="LA-13511 Eq 59 / Tables 59-61 / Forster's worked example",
    notes=(
        "Sood Tables 59/60/61 are constructed by Forster (Ref. 38) so "
        "that f_23 = 4 and f_13 = 15 give k_inf = 1.60 exactly with "
        "φ_2/φ_3 = 0.480 and φ_1/φ_3 = 0.150 (Eqs 60-65). All match "
        "to machine precision."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# 6G k_inf infinite-medium case (LA-13511 Tables 62-67)
# ═══════════════════════════════════════════════════════════════════
#
# URR-6-0-IN is built from 2 coupled 3-group blocks (groups 6,5,4
# mirror groups 1,2,3). Decoupled in scattering, coupled only via χ.
# Sood guarantees same k_inf and flux ratios as URR-3-0-IN.

URR_6_0_IN = La13511Case(
    case_id="URR-6-0-IN",
    problem_number=75,
    description="URR bare infinite medium, 6G isotropic, *with* thermal upscatter",
    materials={0: make_mixture(
        sig_t=np.array([0.240, 0.975, 3.10, 3.10, 0.975, 0.240]),
        sig_c=np.array([0.006, 0.040, 0.20, 0.20, 0.040, 0.006]),
        sig_f=np.array([0.006, 0.060, 0.90, 0.90, 0.060, 0.006]),
        nu=np.array([3.0, 2.5, 2.0, 2.0, 2.5, 3.0]),
        chi=np.array([0.48, 0.02, 0.0, 0.0, 0.02, 0.48]),
        # ORPHEUS [from, to].  Sood: g=6 fast → g=1 slow.
        # ORPHEUS-index ↔ Sood-index map: 0↔6, 1↔5, 2↔4, 3↔3, 4↔2, 5↔1.
        # Scattering structure: top 3 groups (Sood 6,5,4 = ORPHEUS 0,1,2)
        # downscatter only; bottom 3 (Sood 3,2,1 = ORPHEUS 3,4,5) upscatter
        # only; the two sets DECOUPLED in scattering (only χ links them).
        sig_s=np.array([
            [0.024, 0.171, 0.033, 0.0,   0.0,   0.0  ],   # from Sood 6 (fast)
            [0.0,   0.60,  0.275, 0.0,   0.0,   0.0  ],   # from Sood 5
            [0.0,   0.0,   2.0,   0.0,   0.0,   0.0  ],   # from Sood 4 (self only)
            [0.0,   0.0,   0.0,   2.0,   0.0,   0.0  ],   # from Sood 3 (self only)
            [0.0,   0.0,   0.0,   0.275, 0.60,  0.0  ],   # from Sood 2 (up to Sood 3, self)
            [0.0,   0.0,   0.0,   0.033, 0.171, 0.024],   # from Sood 1 (up to Sood 3,2, self)
        ]),
    )},
    geometry_spec=GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=6,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.60,
        flux_ratio_groupwise={
            0: 1.0,    # Sood g6 fast (ORPHEUS 0)
            1: 0.480,  # Sood g5 mid  (= φ_5/φ_6)
            2: 0.150,  # Sood g4 slow (= φ_4/φ_6)
            3: 0.150,  # Sood g3 slow (mirror of 4)
            4: 0.480,  # Sood g2 mid  (mirror of 5)
            5: 1.0,    # Sood g1 fast (mirror of 6)
        },
    ),
    sood_table=62,
    primary_reference="LA-13511 Tables 62-67 / O'Dell (Ref. 39) private comm.",
    notes=(
        "6G == 2 coupled URR-3-0-IN blocks. Same k_inf=1.60. Has "
        "thermal-upscatter pattern in the bottom 3 groups (Σ_21s=0.171, "
        "Σ_31s=0.033, Σ_32s=0.275). compute_kinf_mg is the only Branch-2 "
        "entry that handles this case — Eq-29 specialisation cannot."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# 1G slab F_N bare-critical cases (Sood Tables 6, 7, 13, 17)
# ═══════════════════════════════════════════════════════════════════


# Case 2 — PUa-1-0-SL: Pu-239 (a), c=1.50 slab
PUA_1_0_SL = La13511Case(
    case_id="PUa-1-0-SL",
    problem_number=2,
    description="Pu-239 (a) bare slab, 1G isotropic, c=1.50",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=3.24,
        sigma_s_self=0.225216,
    )},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.605055,
        critical_dimension_cm=1.853722,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=6,
    primary_reference="Lathrop-Leonard 1965 NSE 22, 115 (Ref. 9)",
    notes=(
        "F_N solver at N=12 reaches err ≤ 2e-6 vs Sood truth (well "
        "within the 1e-5 tolerance). Highest c in the bare 1G slab "
        "family (c=1.50)."
    ),
)


# Case 6 — PUb-1-0-SL: Pu-239 (b), c=1.40 slab
PUB_1_0_SL = La13511Case(
    case_id="PUb-1-0-SL",
    problem_number=6,
    description="Pu-239 (b) bare slab, 1G isotropic, c=1.40",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=2.84,
        sigma_s_self=0.225216,
    )},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.73660355,
        critical_dimension_cm=2.256751,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.25: 0.9701734,
            0.50: 0.8810540,
            0.75: 0.7318131,
            1.00: 0.4902592,
        },
    ),
    sood_table=7,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94 (Ref. 26)",
    notes="Slab F_N at N=12 reaches err ≤ 3e-6 on a_c.",
)


# Case 22 — UD2O-1-0-SL: U-D2O reactor, c=1.02 slab
UD2O_1_0_SL = La13511Case(
    case_id="UD2O-1-0-SL",
    problem_number=22,
    description="U-D2O reactor bare slab, 1G isotropic, c=1.02",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.70,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=5.6655054562,
        critical_dimension_cm=10.371065,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.25: 0.93945236,
            0.50: 0.76504084,
            0.75: 0.49690627,
            1.00: 0.13893858,
        },
    ),
    sood_table=17,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94 (Ref. 26)",
    notes=(
        "Lowest c in the bare 1G slab family. Slab F_N at N=12 reaches "
        "err ≤ 2e-6 on a_c. NOTE: F_N at N≥14 fails for low c "
        "(determinant scan loses bracket); use N=12 for this case."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# 1G sphere F_N bare-critical cases (Sood Tables 7, 13, 17)
# ═══════════════════════════════════════════════════════════════════


# Case 8 — PUb-1-0-SP: Pu-239 (b), c=1.40 sphere
PUB_1_0_SP = La13511Case(
    case_id="PUb-1-0-SP",
    problem_number=8,
    description="Pu-239 (b) bare sphere, 1G isotropic, c=1.40",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=2.84,
        sigma_s_self=0.225216,
    )},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=1.9853434324,
        critical_dimension_cm=6.082547,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.25: 0.93538006,
            0.50: 0.75575352,
            0.75: 0.49884364,
            1.00: 0.19222603,
        },
    ),
    sood_table=7,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94 (Ref. 26)",
    notes="Sphere F_N at N=10 reaches err ≤ 5e-8 on R_c (well within 1e-5).",
)


# Case 24 — UD2O-1-0-SP: U-D2O, c=1.02 sphere
UD2O_1_0_SP = La13511Case(
    case_id="UD2O-1-0-SP",
    problem_number=24,
    description="U-D2O reactor bare sphere, 1G isotropic, c=1.02",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.70,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=12.0275320980,
        critical_dimension_cm=22.017156,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.25: 0.91063756,
            0.50: 0.67099621,
            0.75: 0.35561622,
            1.00: 0.04678614,
        },
    ),
    sood_table=17,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94 (Ref. 26)",
    notes="Sphere F_N at N=10 reaches err ≤ 4e-8 on R_c.",
)


# ═══════════════════════════════════════════════════════════════════
# 1G cylinder bare-critical cases — STUBS (B1 dispatch will activate)
# ═══════════════════════════════════════════════════════════════════
#
# Truth values from Sood Tables 7, 17. These will be activated when
# the B1 cylinder solver (Westfall-Metcalf 1973 singular eigenfunction
# expansion) is shipped. NO solver tests are added in Phase B3 for
# these cases — only the registry entries.


# Case 7 — PUb-1-0-CY: Pu-239 (b), c=1.40 cylinder
PUB_1_0_CY_STUB = La13511Case(
    case_id="PUb-1-0-CY",
    problem_number=7,
    description="Pu-239 (b) bare cylinder, 1G isotropic, c=1.40 — STUB",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=2.84,
        sigma_s_self=0.225216,
    )},
    geometry_spec=GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=1.396979,
        critical_dimension_cm=4.279960,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.50: 0.8093,
            1.00: 0.2926,
        },
    ),
    sood_table=7,
    primary_reference="Westfall 1983 Trans. ANS 44, 281 / Westfall-Metcalf 1972 (Refs. 27,28)",
    notes=(
        "STUB: solver activated by B1 dispatch (Westfall-Metcalf 1973 "
        "cylinder F_N). Sood publishes flux ratios only at r/r_c = 0.5 "
        "and 1.0 (Table 8) to 4 digits. Truth values verified."
    ),
)


# Case 23 — UD2O-1-0-CY: U-D2O reactor, c=1.02 cylinder
UD2O_1_0_CY_STUB = La13511Case(
    case_id="UD2O-1-0-CY",
    problem_number=23,
    description="U-D2O reactor bare cylinder, 1G isotropic, c=1.02 — STUB",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.70,
        sigma_s_self=0.464338,
    )},
    geometry_spec=GeometrySpec(
        geometry="cylinder",
        critical_dimension_mfp=9.043255,
        critical_dimension_cm=16.554249,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=17,
    primary_reference="Westfall-Metcalf 1972/1973 (Refs. 27,28)",
    notes=(
        "STUB: solver activated by B1 dispatch. No flux ratios published "
        "for this case in Sood's tables."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# 2G bare-critical STUBS (machinery deferred — needs Siewert-Thomas
# 1986 2G F_N or equivalent)
# ═══════════════════════════════════════════════════════════════════
#
# These cases have published truth values in LA-13511 but require 2G
# F_N machinery that is NOT yet implemented in ORPHEUS. Registered as
# stubs so future implementations can pick them up directly.


def _mix_pu_2g_for_finite_geometry() -> Mixture:
    """Pu-239 2G XS shared by PU-2-0-SL and PU-2-0-SP (Sood Tables 30-31).

    Same XS used in PU-2-0-IN — the difference between the three cases
    is only the geometry (infinite vs slab vs sphere).
    """
    return _mix_2g_isotropic(
        sigma_t_fast=0.2208, sigma_t_slow=0.3360,
        sigma_c_fast=0.00480, sigma_c_slow=0.0144,
        sigma_f_fast=0.0936, sigma_f_slow=0.08544,
        nu_fast=3.10, nu_slow=2.93,
        chi_fast=0.575, chi_slow=0.425,
        sigma_22s=0.0792, sigma_11s=0.23616,
        sigma_12s=0.0432, sigma_21s=0.0,
    )


PU_2_0_SL_STUB = La13511Case(
    case_id="PU-2-0-SL",
    problem_number=45,
    description="Pu-239 bare slab, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_pu_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.396469,
        critical_dimension_cm=1.795602,
        n_groups=2,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=32,
    primary_reference="Siewert-Thomas 1986 NSE 94, 264 / Forster 1970 thesis (Refs. 8, 35, 36)",
    notes="STUB: needs Siewert-Thomas 1986 2G F_N slab machinery (not yet implemented).",
)


PU_2_0_SP_STUB = La13511Case(
    case_id="PU-2-0-SP",
    problem_number=46,
    description="Pu-239 bare sphere, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_pu_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=1.15513,
        critical_dimension_cm=5.231567,
        n_groups=2,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=32,
    primary_reference="Siewert-Thomas 1986 NSE 94, 264 (Ref. 8)",
    notes=(
        "STUB: needs Siewert-Thomas 1986 2G F_N sphere machinery. The "
        "slab and sphere F_N share the geometry-sign abstraction in "
        "fn_method.core; extending to 2G requires the matrix dispersion "
        "law (Λ matrix; Case eigenvalues are 2x2 matrix roots not "
        "scalars). High priority follow-on after B1 cylinder lands."
    ),
)


def _mix_u_2g_for_finite_geometry() -> Mixture:
    """U-235 2G XS shared by U-2-0-SL and U-2-0-SP (Sood Tables 33-34)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.2160, sigma_t_slow=0.3456,
        sigma_c_fast=0.00384, sigma_c_slow=0.01344,
        sigma_f_fast=0.06192, sigma_f_slow=0.06912,
        nu_fast=2.70, nu_slow=2.50,
        chi_fast=0.575, chi_slow=0.425,
        sigma_22s=0.078240, sigma_11s=0.26304,
        sigma_12s=0.0720, sigma_21s=0.0,
    )


U_2_0_SL_STUB = La13511Case(
    case_id="U-2-0-SL",
    problem_number=48,
    description="U-235 bare slab, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_u_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.649377,
        critical_dimension_cm=3.006375,
        n_groups=2,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=35,
    primary_reference="Siewert-Thomas 1986 / Forster 1970 thesis (Refs. 8, 35, 36)",
    notes="STUB: needs 2G F_N slab machinery.",
)


U_2_0_SP_STUB = La13511Case(
    case_id="U-2-0-SP",
    problem_number=49,
    description="U-235 bare sphere, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_u_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=1.70844,
        critical_dimension_cm=7.909444,
        n_groups=2,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=35,
    primary_reference="Siewert-Thomas 1986 NSE 94, 264 (Ref. 8)",
    notes="STUB: needs 2G F_N sphere machinery.",
)


def _mix_ual_2g_for_finite_geometry() -> Mixture:
    """U-Al-Water 2G XS shared by UAL-2-0-SL and UAL-2-0-SP (Tables 36-37)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.26817, sigma_t_slow=1.27698,
        sigma_c_fast=0.000222, sigma_c_slow=0.00314363958,
        sigma_f_fast=0.0, sigma_f_slow=0.06070636042,
        nu_fast=0.0, nu_slow=2.83,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.247516, sigma_11s=1.21313,
        sigma_12s=0.020432, sigma_21s=0.0,
    )


UAL_2_0_SL_STUB = La13511Case(
    case_id="UAL-2-0-SL",
    problem_number=51,
    description="U-Al-Water bare slab, 2G isotropic — STUB",
    materials={0: _mix_ual_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=2.09994,
        critical_dimension_cm=7.830630,
        n_groups=2,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=38,
    primary_reference="Siewert-Thomas 1986 / Forster 1970 thesis (Refs. 8, 35, 36)",
    notes="STUB: needs 2G F_N slab machinery.",
)


UAL_2_0_SP_STUB = La13511Case(
    case_id="UAL-2-0-SP",
    problem_number=52,
    description="U-Al-Water bare sphere, 2G isotropic — STUB",
    materials={0: _mix_ual_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=4.73786,
        critical_dimension_cm=17.66738,
        n_groups=2,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=38,
    primary_reference="Siewert-Thomas 1986 NSE 94, 264 (Ref. 8)",
    notes="STUB: needs 2G F_N sphere machinery.",
)


def _mix_urra_2g_for_finite_geometry() -> Mixture:
    """URRa 2G XS shared by URRa-2-0-SL and URRa-2-0-SP (Tables 39-40)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.65696, sigma_t_slow=2.52025,
        sigma_c_fast=0.0010046, sigma_c_slow=0.025788,
        sigma_f_fast=0.0010484, sigma_f_slow=0.050632,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.62568, sigma_11s=2.44383,
        sigma_12s=0.029227, sigma_21s=0.0,
    )


URRA_2_0_SL_STUB = La13511Case(
    case_id="URRa-2-0-SL",
    problem_number=54,
    description="URR (a) bare slab, 2G isotropic — STUB",
    materials={0: _mix_urra_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=4.97112,
        critical_dimension_cm=7.566853,
        n_groups=2,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            # Sood Table 42 — fast group, normalised to fast at center
            0.241394: 0.943363,
            0.502905: 0.761973,
            0.744300: 0.504012,
            1.0: 0.147598,
        },
    ),
    sood_table=41,
    primary_reference="Siewert-Thomas 1986 / Forster 1970 / Stewart 1974 (Refs. 8, 35, 36)",
    notes=(
        "STUB: needs 2G F_N slab machinery. Sood Table 42 gives 2G "
        "flux ratios at four spatial points (fast + slow). "
        "flux_ratios stored here is the FAST group; the slow-group "
        "ratio at the same points is (0.340124, 0.273056, 0.173845, "
        "0.0212324)."
    ),
)


URRA_2_0_SP_STUB = La13511Case(
    case_id="URRa-2-0-SP",
    problem_number=55,
    description="URR (a) bare sphere, 2G isotropic — STUB",
    materials={0: _mix_urra_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=10.5441,
        critical_dimension_cm=16.049836,
        n_groups=2,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=41,
    primary_reference="Siewert-Thomas 1986 NSE 94, 264 (Ref. 8)",
    notes="STUB: needs 2G F_N sphere machinery.",
)


def _mix_ud2o_2g_for_finite_geometry() -> Mixture:
    """U-D2O 2G XS shared by UD2O-2-0-SL and UD2O-2-0-SP (Tables 49-50)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.33588, sigma_t_slow=0.54628,
        sigma_c_fast=0.008708, sigma_c_slow=0.02518,
        sigma_f_fast=0.002817, sigma_f_slow=0.097,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.31980, sigma_11s=0.42410,
        sigma_12s=0.004555, sigma_21s=0.0,
    )


UD2O_2_0_SL_STUB = La13511Case(
    case_id="UD2O-2-0-SL",
    problem_number=68,
    description="U-D2O reactor bare slab, 2G isotropic — STUB",
    materials={0: _mix_ud2o_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=284.367,
        critical_dimension_cm=846.632726,
        n_groups=2,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=51,
    primary_reference="Siewert-Thomas 1986 / Forster 1970 / Stewart 1974 (Refs. 8, 35, 36)",
    notes=(
        "STUB: needs 2G F_N slab machinery. Critical dimension is "
        "VERY LARGE (284 mfp — barely-supercritical heavy-water "
        "reactor); high N_F may be needed."
    ),
)


UD2O_2_0_SP_STUB = La13511Case(
    case_id="UD2O-2-0-SP",
    problem_number=69,
    description="U-D2O reactor bare sphere, 2G isotropic — STUB",
    materials={0: _mix_ud2o_2g_for_finite_geometry()},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=569.430,
        critical_dimension_cm=1695.337621,
        n_groups=2,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=51,
    primary_reference="Siewert-Thomas 1986 NSE 94, 264 (Ref. 8)",
    notes="STUB: needs 2G F_N sphere machinery. Critical R ~ 1695 cm — heavy-water reactor.",
)


# ═══════════════════════════════════════════════════════════════════
# Public registry — name → case lookup
# ═══════════════════════════════════════════════════════════════════


ALL_FIRST_SLICE: tuple[La13511Case, ...] = (
    PUA_1_0_IN,
    PU_2_0_IN,
    UA_1_0_SL_STUB,
    UA_1_0_CY_STUB,
    UA_1_0_SP_STUB,
)
"""Phase A first-slice cases (5). Cases 1+5 (k_inf) ship with full
Branch-2 solvers; cases 2-4 are bare-critical slab/cylinder/sphere —
slab and sphere F_N shipped, cylinder F_N pending Westfall-Metcalf."""


WIDE_SLICE_KINF: tuple[La13511Case, ...] = (
    # 1G k_inf (the 9 1G isotropic + anisotropic infinite-medium cases —
    # PUa already in FIRST_SLICE, the rest are new):
    PUB_1_0_IN, UA_1_0_IN, UB_1_0_IN, UC_1_0_IN, UD_1_0_IN,
    UD2O_1_0_IN, UE_1_0_IN, PU_1_1_IN,
    UD2OA_1_1_IN, UD2OB_1_1_IN, UD2OC_1_1_IN,
    # 2G k_inf (PU already in FIRST_SLICE; these 7 are new):
    U_2_0_IN, UAL_2_0_IN, URRA_2_0_IN, URRB_2_0_IN, URRC_2_0_IN,
    URRD_2_0_IN, UD2O_2_0_IN,
    # 3G k_inf:
    URR_3_0_IN,
    # 6G k_inf:
    URR_6_0_IN,
)
"""Phase B3 wide-slice k_inf cases (19). All solved by ``compute_kinf_*``;
no spatial discretisation involved."""


WIDE_SLICE_BARE_CRITICAL_1G: tuple[La13511Case, ...] = (
    # Slab F_N (Ua already in FIRST_SLICE; these 3 are new):
    PUA_1_0_SL, PUB_1_0_SL, UD2O_1_0_SL,
    # Sphere F_N (Ua already in FIRST_SLICE; these 2 are new):
    PUB_1_0_SP, UD2O_1_0_SP,
)
"""Phase B3 wide-slice 1G bare-critical cases activated by the
existing slab/sphere F_N solvers (5)."""


# ═══════════════════════════════════════════════════════════════════
# Wave 2-C — 1G P_1 anisotropic bare-critical slab + sphere cases
# ═══════════════════════════════════════════════════════════════════
#
# Bare-critical slab/sphere with linearly anisotropic scattering.
# Verified by :mod:`...galerkin_spectral`. CRITICAL convention: Sood's
# Σ_s1 is the scattering-only anisotropy moment; Dahl-Sjostrand 1979
# uses μ̄ = mean cosine of all secondaries (scattering + fission,
# fission isotropic). Conversion: μ̄_eff = Σ_s1/(c·Σ_t).


def _mix_1g_anisotropic(
    sigma_t: float, sigma_c: float, sigma_f: float, nu: float,
    sigma_s_self: float, sigma_s1_self: float,
) -> Mixture:
    """Build a 1G P_1 anisotropic Mixture from raw Sood XS components."""
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sigma_c]),
        sig_f=np.array([sigma_f]),
        nu=np.array([nu]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s_self]]),
        sig_s1=np.array([[sigma_s1_self]]),
    )


PUA_1_1_SL = La13511Case(
    case_id="PUa-1-1-SL",
    problem_number=32,
    description="Pu-239 (a) bare slab, 1G P_1 anisotropic (forward), c=1.40",
    materials={0: _mix_1g_anisotropic(
        sigma_t=1.0, sigma_c=0.0, sigma_f=0.266667, nu=2.5,
        sigma_s_self=0.733333, sigma_s1_self=0.20,
    )},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.77032,
        critical_dimension_cm=0.77032,
        n_groups=1,
        bc_left=BC.vacuum, bc_right=BC.vacuum,
    ),
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=25,
    primary_reference="Sood Table 25 / problem 32 / Sanchez 1976 (Ref. 30)",
    notes="Σ_s1=0.20. Carlvik-Galerkin uses μ̄_eff = 0.20/1.40 = 0.142857.",
)

PUB_1_1_SL = La13511Case(
    case_id="PUb-1-1-SL",
    problem_number=34,
    description="Pu-239 (b) bare slab, 1G P_1 anisotropic (strong forward), c=1.40",
    materials={0: _mix_1g_anisotropic(
        sigma_t=1.0, sigma_c=0.0, sigma_f=0.266667, nu=2.5,
        sigma_s_self=0.733333, sigma_s1_self=0.333333,
    )},
    geometry_spec=GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.79606,
        critical_dimension_cm=0.79606,
        n_groups=1,
        bc_left=BC.vacuum, bc_right=BC.vacuum,
    ),
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=25,
    primary_reference="Sood Table 25 / problem 34 / Sanchez 1976 (Ref. 30)",
    notes="Σ_s1=0.333333 (negative scattering for μ near -1). μ̄_eff = 0.238095.",
)

UD2OA_1_1_SP = La13511Case(
    case_id="UD2Oa-1-1-SP",
    problem_number=39,
    description="U-D2O (a) bare sphere, 1G P_1 anisotropic, c=1.0308381",
    materials={0: _mix_1g_anisotropic(
        sigma_t=0.54628, sigma_c=0.027314, sigma_f=0.054628, nu=1.808381,
        sigma_s_self=0.464338, sigma_s1_self=0.056312624,
    )},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=10.0,
        critical_dimension_cm=18.30563081,
        n_groups=1,
        bc_left=BC.reflective, bc_right=BC.vacuum,
    ),
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=29,
    primary_reference="Sood Table 29 / problem 39 / Mitsis 1963 (Ref. 15)",
    notes="μ̄_eff = 0.10 — matches DS Table I row d=20, μ̄=0.10.",
)

UD2OB_1_1_SP = La13511Case(
    case_id="UD2Ob-1-1-SP",
    problem_number=41,
    description="U-D2O (b) bare sphere, 1G P_1 anisotropic, c=1.0341086",
    materials={0: _mix_1g_anisotropic(
        sigma_t=0.54628, sigma_c=0.027314, sigma_f=0.054628, nu=1.841086,
        sigma_s_self=0.464338, sigma_s1_self=0.112982569,
    )},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=10.0,
        critical_dimension_cm=18.30563081,
        n_groups=1,
        bc_left=BC.reflective, bc_right=BC.vacuum,
    ),
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=29,
    primary_reference="Sood Table 29 / problem 41 / Mitsis 1963 (Ref. 15)",
    notes="μ̄_eff = 0.20 — matches DS Table I row d=20, μ̄=0.20.",
)

UD2OC_1_1_SP = La13511Case(
    case_id="UD2Oc-1-1-SP",
    problem_number=43,
    description="U-D2O (c) bare sphere, 1G P_1 anisotropic (back-peaked!), c=1.01964",
    materials={0: _mix_1g_anisotropic(
        sigma_t=0.54628, sigma_c=0.027314, sigma_f=0.054628, nu=1.6964,
        sigma_s_self=0.464338, sigma_s1_self=-0.27850447,
    )},
    geometry_spec=GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=10.0,
        critical_dimension_cm=18.30563081,
        n_groups=1,
        bc_left=BC.reflective, bc_right=BC.vacuum,
    ),
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0),
    sood_table=29,
    primary_reference="Sood Table 29 / problem 43 / Boffi-Molinari-Spiga 1977 (Ref. 16)",
    notes="μ̄_eff = -0.50 (back-peaked). Outside Dahl-Sjostrand table coverage.",
)


WIDE_SLICE_BARE_CRITICAL_1G_P1: tuple[La13511Case, ...] = (
    PUA_1_1_SL, PUB_1_1_SL,
    UD2OA_1_1_SP, UD2OB_1_1_SP, UD2OC_1_1_SP,
)
"""Wave 2-C P_1 anisotropic bare-critical cases (5)."""


WIDE_SLICE_STUBS: tuple[La13511Case, ...] = (
    # Cylinder stubs (Ua already in FIRST_SLICE; these 2 are new):
    PUB_1_0_CY_STUB, UD2O_1_0_CY_STUB,
    # 2G bare-critical stubs (need Siewert-Thomas 1986 2G F_N):
    PU_2_0_SL_STUB, PU_2_0_SP_STUB,
    U_2_0_SL_STUB, U_2_0_SP_STUB,
    UAL_2_0_SL_STUB, UAL_2_0_SP_STUB,
    URRA_2_0_SL_STUB, URRA_2_0_SP_STUB,
    UD2O_2_0_SL_STUB, UD2O_2_0_SP_STUB,
)
"""Phase B3 wide-slice STUBS — registered but no solver tests added.
Cylinder cases activated by B1 dispatch; 2G bare-critical cases need
Siewert-Thomas 1986 2G F_N machinery."""


_ALL_CASES: tuple[La13511Case, ...] = (
    *ALL_FIRST_SLICE,
    *WIDE_SLICE_KINF,
    *WIDE_SLICE_BARE_CRITICAL_1G,
    *WIDE_SLICE_BARE_CRITICAL_1G_P1,
    *WIDE_SLICE_STUBS,
)


LA13511_CASES: dict[str, La13511Case] = {
    case.case_id: case for case in _ALL_CASES
}
"""Name → case mapping for ergonomic test access:

>>> from orpheus.derivations.continuous.sood_registry import LA13511_CASES
>>> case = LA13511_CASES["PUa-1-0-IN"]
>>> mixture = case.materials[0]

Coverage: Phase A first slice (5 cases) + Phase B3 wide enumeration
(36 cases) = 41 cases total. Of these:

* 19 ``k_inf`` cases (1G/2G/3G/6G) verified by ``compute_kinf_*``.
* 7 1G bare-critical cases (slab + sphere) verified by F_N solvers.
* 2 cylinder STUBS (B1 dispatch will activate).
* 12 2G bare-critical STUBS (Siewert-Thomas 1986 2G F_N deferred).
* 1 multi-region slab case (UA_1_0_SL_STUB has flux ratios used by
  existing F_N test_slab_xverif).
"""

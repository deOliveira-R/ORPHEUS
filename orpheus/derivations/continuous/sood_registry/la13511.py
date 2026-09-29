r"""Sood, Forster & Parsons (2003) benchmark case catalogue.

The cases are cited to the journal edition, Sood, Forster & Parsons,
"Analytical benchmark test set for criticality code verification",
*Prog. Nucl. Energy* 42(1), 55-106 (2003). The 1999 report LA-13511
is the same test set under other table, equation and reference numbers,
and for the UAL and UD2O two-group sets with other values; the module
and class names keep the report's number.

Production-protocol-aligned port of the legacy
``orpheus.derivations.continuous.fn_method.benchmarks.la13511`` module.
Each case now carries:

* **`materials: dict[int, Mixture]`** — keyed by integer material ID,
  exactly what :func:`orpheus.cp.solver.solve_cp` and
  :func:`orpheus.sn.solver.solve_sn` consume. Built via
  :func:`orpheus.derivations.common.xs_library.make_mixture` from the
  raw Sood XS components (ν, Σ_f, Σ_c, Σ_s, χ).
* **`geometry_kind: str`** — one of ``"slab"``, ``"sphere"``,
  ``"cylinder"`` (finite cases) or ``"infinite"`` (k_inf cases).
  :meth:`La13511Case.to_geometry` materialises a
  :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
  on demand (raises for the infinite kind). Discrete consumers
  build a mesh via :meth:`Mesh1D.from_geometry`; reference solvers
  consume the geometry directly.
* **`truth: La13511Truth`** — all published reference values bundled
  in one struct. Different cases populate different subsets,
  including the published ``critical_dimension_mfp`` (the registry-
  truth artifact lifted off the legacy ``GeometrySpec`` in Phase B).

ORPHEUS convention vs Sood convention
-------------------------------------

Sood et al. number energy groups :math:`g=N` (fast) → :math:`g=1`
(slow), the reverse of typical nuclear-engineering convention. This
module does the conversion at construction time; consumers see XS
arrays in **ORPHEUS convention** (``g=0`` fast, ``g=N-1`` slow).

The scattering matrix uses ORPHEUS's ``[from, to]`` convention:
``sigma_s[g, h]`` is :math:`\Sigma_{s, g \to h}`. The
``Sigma_g^{rem}`` removal cross sections used by Sood Eqs (A.9)-(A.10) are
therefore ``sigma_t[g] - sigma_s[g, g]``.

Provenance
----------

* All XS values: Sood 2003 Tables 2-64.
* All :math:`k_\infty`, :math:`r_c`, flux ratios: Sood 2003 Sections 4
  (1G), 5 (2G), 6 (3G) and 7 (6G), with primary references cited per
  case.

Python value precision: published values transcribed verbatim. Where
Sood reports e.g. "k_inf = 2.612903" with 6 published digits, the
catalogue stores 2.612903 exactly as a Python float. Tolerance for
verification is 1e-5 (i.e., 5 significant figures match), so the
6th-digit rounding of the published value is not load-bearing.
"""
from __future__ import annotations

import numpy as np

from orpheus.data.citation import Citation
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import make_mixture

from .case import La13511Case, La13511Truth


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


def _mix_pua_1g() -> Mixture:
    """The PUa one-group set, shared by every case that uses it."""
    return _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=3.24,
        sigma_s_self=0.225216,
    )


def _mix_pub_1g() -> Mixture:
    """The PUb one-group set, shared by every case that uses it."""
    return _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=2.84,
        sigma_s_self=0.225216,
    )


def _mix_ud2o_1g() -> Mixture:
    """The UD2O one-group set, shared by every case that uses it."""
    return _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.70,
        sigma_s_self=0.464338,
    )


# ═══════════════════════════════════════════════════════════════════
# Case 1 — PUa-1-0-IN (Sood problem 1): 1G infinite medium, Pu-239 (a)
# ═══════════════════════════════════════════════════════════════════
#
# Sood 2003 Table 2 (Pu-239 (a) cross sections, 1G isotropic):
#   ν = 3.24,  Σ_f = 0.0816,  Σ_c = 0.019584,  Σ_s = 0.225216,
#   Σ_t = 0.32640,  c = (Σ_s + νΣ_f)/Σ_t = 1.50.
# Reference: k_inf = 2.612903 (Sood 2003 Eq (A.3)).

PUA_1_0_IN = La13511Case(
    case_id="PUa-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 1"),
    description="Pu-239 (a) bare infinite medium, 1G isotropic",
    materials={0: _mix_pua_1g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.612903,
        sources=(Citation("SoodForsterParsons2003", "p. 69"),),
        flux_ratios=None,
    ),
    notes="1G infinite-medium k_inf reduces to nu·Sigma_f/Sigma_a; the 'c' factor in Eq (A.3) cancels algebraically (verified in fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_a3_simplifies_to_eq_a2).",
)


# ═══════════════════════════════════════════════════════════════════
# Case 5 — PU-2-0-IN (Sood problem 44): 2G infinite medium, Pu
# ═══════════════════════════════════════════════════════════════════
#
# Sood 2003 Tables 27-28 (Pu 2G isotropic, no upscatter):
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
# Reference: k_inf = 2.683767, φ_2/φ_1 = 0.675229 (Sood Eqs (A.11)-(A.12) + Eq (A.15)).

PU_2_0_IN = La13511Case(
    case_id="PU-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 44"),
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
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.683767,
        sources=(Citation("SoodForsterParsons2003", "p. 81"),),
        flux_ratios=None,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 0.675229},  # slow/fast = 1/Sood ratio
    ),
    notes="Sood Eq (A.11) (the general two-group k_inf, with upscatter) reduces to Eq (A.12) term by term when Sigma_21s = 0, and as printed gives the published k_inf = 2.683767 (checked against the 1999 page image, 2026-09-29; the 2003 body is the same equation). An earlier note here called Eq (A.11) a typo that reduces to 2.862; that was a mis-transcription of the equation (a swapped Sigma_g^rem factor), not an error in the report. The SymPy derivation in fn_method.origins.k_inf_derivations.derive_kinf_2g_general_from_matrix computes det(M)=0 from Eq (A.8) directly and verifies Eq (A.12) against it.",
)


# ═══════════════════════════════════════════════════════════════════
# Cases 2-4 — bare critical sphere/slab/cylinder of U-235 (a), 1G
# ═══════════════════════════════════════════════════════════════════
#
# Common XS (Sood 2003 Table 9, U-235 (a) 1G isotropic):
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
    problem=Citation("SoodForsterParsons2003", "problem 12"),
    description="U-235 (a) bare slab, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 10"), Citation("SoodForsterParsons2003", "Table 11"), Citation("KaperLindemanLeaf1974")),
        flux_ratios={
            0.25: 0.9669506,
            0.50: 0.8686259,
            0.75: 0.7055218,
            1.00: 0.4461912,
        },
        critical_dimension_mfp=0.93772556,
    ),
    notes="Slab F_N solver shipped at ≤5e-6 absolute on a_c (see fn_method.slab.solve_fn_slab_bare_critical). Critical dimension is the half-thickness; slab full width is 2*critical_dimension_cm = 5.745868 cm.",
)


# Case 3 — Ua-1-0-CY (Sood problem 13): 1G bare cylinder, U-235 (a)

UA_1_0_CY_STUB = La13511Case(
    case_id="Ua-1-0-CY",
    problem=Citation("SoodForsterParsons2003", "problem 13"),
    description="U-235 (a) bare infinite cylinder, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_kind="cylinder",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 10"), Citation("Westfall1983"), Citation("WestfallMetcalf1972")),
        critical_dimension_mfp=1.72500292,
    ),
    notes="WM-72 singular-eigenfunction cylinder solver shipped at ~1% relative accuracy (single-cell product integration on the log-singular kernel diagonal; see orpheus.derivations.continuous.singular_eigenfunction.cylinder). Variant α cylinder cross-check at 8.5e-6 holds the strict 1e-5 anchor. WM-72 prototype provides the second, structurally-independent cross-check anchor (different mathematical pillar than Variant α / Bickley-Naylor).",
)


# Case 4 — Ua-1-0-SP (Sood problem 14): 1G bare sphere, U-235 (a)

UA_1_0_SP_STUB = La13511Case(
    case_id="Ua-1-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 14"),
    description="U-235 (a) bare sphere, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 10"), Citation("KaperLindemanLeaf1974"), Citation("KaperLindemanLeaf1974", "Table VII")),
        flux_ratios={
            # Kaper-Lindeman-Leaf 1974 Table VII at c=1.30 — the same
            # XS as Sood Ua-1-0-SP (U-235 (a) 1G isotropic, c=1.30).
            0.25: 0.93244907,
            0.50: 0.74553332,
            0.75: 0.48095413,
            1.00: 0.17177706,
        },
        critical_dimension_mfp=2.4248249802,
    ),
    notes="Sphere F_N solver shipped at ≤1e-7 absolute on R_c (see fn_method.sphere.solve_fn_sphere_bare_critical). Used as the structurally-independent L1 reference for Variant α sphere. Flux ratios populated from KLL Table VII c=1.30 row — same XS as this case (cross-check via fn_method.sphere.flux_reconstruction).",
)


# ═══════════════════════════════════════════════════════════════════
# Phase B3 — wide enumeration of currently-implementable cases
# ═══════════════════════════════════════════════════════════════════
#
# Cases below are added by Phase B3 of the wide-enumeration sweep
# (closeout: ``.claude/agent-memory/method-implementer/
# sood_registry_wide_enumeration_phase_b3.md``). They cover the subset
# of Sood's 75 problems that the existing ``fn_method`` machinery
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


def _mix_pu_2g() -> Mixture:
    """Pu-239 2G XS shared by the infinite-medium and finite cases, PU-2-0-SL and PU-2-0-SP (Sood 2003 Tables 27-28).

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


def _mix_u_2g() -> Mixture:
    """U-235 2G XS shared by the infinite-medium and finite cases, U-2-0-SL and U-2-0-SP (Sood 2003 Tables 30-31)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.2160, sigma_t_slow=0.3456,
        sigma_c_fast=0.00384, sigma_c_slow=0.01344,
        sigma_f_fast=0.06192, sigma_f_slow=0.06912,
        nu_fast=2.70, nu_slow=2.50,
        chi_fast=0.575, chi_slow=0.425,
        sigma_22s=0.078240, sigma_11s=0.26304,
        sigma_12s=0.0720, sigma_21s=0.0,
    )


def _mix_ual_2g() -> Mixture:
    """U-Al-Water 2G XS shared by the infinite-medium and finite cases, UAL-2-0-SL and UAL-2-0-SP (Sood 2003 Tables 33-34)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.268165, sigma_t_slow=1.276976,
        sigma_c_fast=0.000217, sigma_c_slow=0.003143,
        sigma_f_fast=0.0, sigma_f_slow=0.060706,
        nu_fast=0.0, nu_slow=2.830023,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.247516, sigma_11s=1.213127,
        sigma_12s=0.020432, sigma_21s=0.0,
    )


def _mix_urra_2g() -> Mixture:
    """URRa 2G XS shared by the infinite-medium and finite cases, URRa-2-0-SL and URRa-2-0-SP (Sood 2003 Tables 36-37)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.65696, sigma_t_slow=2.52025,
        sigma_c_fast=0.0010046, sigma_c_slow=0.025788,
        sigma_f_fast=0.0010484, sigma_f_slow=0.050632,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.62568, sigma_11s=2.44383,
        sigma_12s=0.029227, sigma_21s=0.0,
    )


def _mix_ud2o_2g() -> Mixture:
    """U-D2O 2G XS shared by the infinite-medium and finite cases, UD2O-2-0-SL and UD2O-2-0-SP (Sood 2003 Tables 46-47)."""
    return _mix_2g_isotropic(
        sigma_t_fast=0.33588, sigma_t_slow=0.54628,
        sigma_c_fast=0.0087078, sigma_c_slow=0.02518,
        sigma_f_fast=0.002817, sigma_f_slow=0.097,
        nu_fast=2.50, nu_slow=2.50,
        chi_fast=1.0, chi_slow=0.0,
        sigma_22s=0.31980, sigma_11s=0.42410,
        sigma_12s=0.0045552, sigma_21s=0.0,
    )


# ═══════════════════════════════════════════════════════════════════
# 1G k_inf infinite-medium cases (Sood 2003 Tables 2, 9, 13, 17)
# ═══════════════════════════════════════════════════════════════════
#
# Pure rational algebra in ν, Σ_f, Σ_s, Σ_t (Sood 2003 Eqs (A.2)/(A.3)). All
# verifiable to machine precision; tolerance is set by Sood's published
# precision (≤ 1e-5 absolute on a 6-7 digit truth value).

# Case 5 — PUb-1-0-IN: Pu-239 (b), c=1.40
PUB_1_0_IN = La13511Case(
    case_id="PUb-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 5"),
    description="Pu-239 (b) bare infinite medium, 1G isotropic, c=1.40",
    materials={0: _mix_pub_1g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.290323, sources=(Citation("SoodForsterParsons2003", "p. 69"),)),
    notes="Same Σ_t / Σ_s as PUa, only ν changes (3.24 → 2.84). c=1.40.",
)


# Case 11 — Ua-1-0-IN: U-235 (a), c=1.30
UA_1_0_IN = La13511Case(
    case_id="Ua-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 11"),
    description="U-235 (a) bare infinite medium, 1G isotropic, c=1.30",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.25, sources=(Citation("SoodForsterParsons2003", "p. 71"),)),
    notes="Sood publishes 'k_inf = 2.25' (3 digits printed but algebraically exact).",
)


# Case 15 — Ub-1-0-IN: U-235 (b), c=1.3194202
UB_1_0_IN = La13511Case(
    case_id="Ub-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 15"),
    description="U-235 (b) bare infinite medium, 1G isotropic, c=1.3194202",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.013056,
        sigma_f=0.065280,
        nu=2.797101,
        sigma_s_self=0.248064,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.330917, sources=(Citation("SoodForsterParsons2003", "p. 71"),)),
    notes="Cross-section variant (b): same Σ_t/Σ_s as Ua, ν tuned to give c=1.3194202.",
)


# Case 17 — Uc-1-0-IN: U-235 (c), c=1.3014616
UC_1_0_IN = La13511Case(
    case_id="Uc-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 17"),
    description="U-235 (c) bare infinite medium, 1G isotropic, c=1.3014616",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.013056,
        sigma_f=0.065280,
        nu=2.707308,
        sigma_s_self=0.248064,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.256083, sources=(Citation("SoodForsterParsons2003", "p. 71"),)),
    notes="Cross-section variant (c): same Σ_t/Σ_s as Ua, ν=2.707308.",
)


# Case 19 — Ud-1-0-IN: U-235 (d), c=1.2958396
UD_1_0_IN = La13511Case(
    case_id="Ud-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 19"),
    description="U-235 (d) bare infinite medium, 1G isotropic, c=1.2958396",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.013056,
        sigma_f=0.065280,
        nu=2.679198,
        sigma_s_self=0.248064,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.232667, sources=(Citation("SoodForsterParsons2003", "p. 71"),)),
    notes="Cross-section variant (d): same Σ_t/Σ_s as Ua, ν=2.679198.",
)


# Case 21 — UD2O-1-0-IN: U-D2O reactor, c=1.02
UD2O_1_0_IN = La13511Case(
    case_id="UD2O-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 21"),
    description="U-D2O reactor bare infinite medium, 1G isotropic, c=1.02",
    materials={0: _mix_ud2o_1g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.133333, sources=(Citation("SoodForsterParsons2003", "p. 73"),)),
    notes="Heavy-water-moderated low-enrichment U; lowest c in the bare 1G family (1.02).",
)


# Case 29 — Ue-1-0-IN: U-235 reactor with Fe/Na surrounds (infinite-medium k_inf only — no Fe/Na in this case)
UE_1_0_IN = La13511Case(
    case_id="Ue-1-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 29"),
    description="U-235 reactor (e) bare infinite medium, 1G isotropic, c=1.230",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.407407,
        sigma_c=0.01013756,
        sigma_f=0.06922744,
        nu=2.50,
        sigma_s_self=0.328042,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=2.1806667, sources=(Citation("SoodForsterParsons2003", "p. 74"),)),
    notes="U-235 (e) cross sections used in the Ue-Fe-Na multi-region case; the infinite-medium variant uses the U-235 (e) XS alone.",
)


# Case 31 — PU-1-1-IN: Pu-239 with linearly anisotropic scattering — k_inf is identical to the isotropic case
PU_1_1_IN = La13511Case(
    case_id="PU-1-1-IN",
    problem=Citation("SoodForsterParsons2003", "problem 31"),
    description="Pu-239 (a/b) bare infinite medium, 1G P_1 anisotropic, c=1.40",
    materials={0: _mix_1g_isotropic(
        sigma_t=1.0,
        sigma_c=0.0,
        sigma_f=0.266667,
        nu=2.5,
        sigma_s_self=0.733333,
    )},
    geometry_kind="infinite",
    scattering_order=0,  # k_inf does not depend on anisotropy; isotropic XS suffice
    truth=La13511Truth(k_eff_or_kinf=2.5, sources=(Citation("SoodForsterParsons2003", "p. 76"),)),
    notes="Sood: 'The anisotropic scattering cross sections do not change k_inf' (Sood 2003 p. 76, section 4.2.1). Catalogued with scattering_order=0 + Σ_s = Σ_s0 since the P_1 moment is a no-op for infinite-medium k_inf. Anisotropic data lives in the slab cases 32-35.",
)


# Case 38, 40, 42 — UD2O-{a,b,c}-1-1-IN: U-D2O P_1 anisotropic — k_inf depends only on isotropic XS
# Each (a,b,c) has slightly different ν tuned to give specified c values 1.0308381 / 1.0341086 / 1.01964.
UD2OA_1_1_IN = La13511Case(
    case_id="UD2Oa-1-1-IN",
    problem=Citation("SoodForsterParsons2003", "problem 38"),
    description="U-D2O (a) bare infinite medium, 1G P_1 anisotropic, c=1.0308381",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.808381,
        sigma_s_self=0.464338,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.205587, sources=(Citation("SoodForsterParsons2003", "p. 77"),)),
    notes="P_1 anisotropy doesn't change k_inf; ν tuned for c=1.0308381.",
)

UD2OB_1_1_IN = La13511Case(
    case_id="UD2Ob-1-1-IN",
    problem=Citation("SoodForsterParsons2003", "problem 40"),
    description="U-D2O (b) bare infinite medium, 1G P_1 anisotropic, c=1.0341086",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.841086,
        sigma_s_self=0.464338,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.227391, sources=(Citation("SoodForsterParsons2003", "p. 77"),)),
    notes="P_1 anisotropy doesn't change k_inf; ν tuned for c=1.0341086.",
)

UD2OC_1_1_IN = La13511Case(
    case_id="UD2Oc-1-1-IN",
    problem=Citation("SoodForsterParsons2003", "problem 42"),
    description="U-D2O (c) bare infinite medium, 1G P_1 anisotropic, c=1.01964",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.54628,
        sigma_c=0.027314,
        sigma_f=0.054628,
        nu=1.6964,
        sigma_s_self=0.464338,
    )},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.130933, sources=(Citation("SoodForsterParsons2003", "p. 77"),)),
    notes="U-D2O (c) has *negative* P_1 scattering moment Σ_s1 = -0.27850447 (backward-peaked); k_inf still depends only on Σ_s0. Slab cases (39/41/43) inherit anisotropy.",
)


# ═══════════════════════════════════════════════════════════════════
# 2G k_inf infinite-medium cases (Sood 2003 Tables 27-28, 30-31, 33-34,
# 36-37, 40-41, 43-44, 46-47)
# ═══════════════════════════════════════════════════════════════════


# Case 47 — U-2-0-IN: U-235, 2G isotropic, no upscatter
U_2_0_IN = La13511Case(
    case_id="U-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 47"),
    description="U-235 bare infinite medium, 2G isotropic, no upscatter",
    materials={0: _mix_u_2g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.216349,
        sources=(Citation("SoodForsterParsons2003", "p. 82"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 0.474967},
    ),
    notes="Sood publishes φ_2/φ_1 (fast/slow) = 0.474967.",
)


# Case 50 — UAL-2-0-IN: Uranium-Aluminum-Water assembly, 2G isotropic, no upscatter
# Note: Σ_2f = 0.0 in fast group — fission only in slow group. χ_fast = 1.0, χ_slow = 0.0.
UAL_2_0_IN = La13511Case(
    case_id="UAL-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 50"),
    description="U-Al-Water assembly bare infinite medium, 2G isotropic",
    materials={0: _mix_ual_2g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.662437,
        sources=(Citation("SoodForsterParsons2003", "p. 83"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 3.124951},
    ),
    notes="Slow-only fission (ν_2 = 0, χ_1 = 0). Sood publishes φ_2/φ_1 = 3.124951 (2003 edition; the primary source, Siewert and Thomas 1986, gives 3.125) — i.e. fast/slow > 1 because slow group is very absorbing (Σ_1 = 1.276976 mostly self-scatter).",
)


# Case 53 — URRa-2-0-IN: 93%-enriched U research reactor (a), 2G isotropic, no upscatter
URRA_2_0_IN = La13511Case(
    case_id="URRa-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 53"),
    description="URR (a) — 93% enriched U bare infinite medium, 2G isotropic",
    materials={0: _mix_urra_2g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.631452,
        sources=(Citation("SoodForsterParsons2003", "p. 84"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 2.614706},
    ),
    notes="93%-enriched bare research reactor. χ_1=0 (all fission from fast group).",
)


# Case 56 — URRb-2-0-IN: research reactor (b) WITH thermal upscatter (Σ_21s = 0.000767)
URRB_2_0_IN = La13511Case(
    case_id="URRb-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 56"),
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
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.365821,
        sources=(Citation("SoodForsterParsons2003", "p. 85"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 1.173679},
    ),
    notes="Has thermal upscatter (Σ_21s = 0.000767). MUST use the general Eq (A.11) formula (compute_kinf_2g_general or compute_kinf_mg), NOT the no-upscatter Eq (A.12) specialisation.",
)


# Case 57 — URRc-2-0-IN: research reactor (c) WITH thermal upscatter (Σ_21s = 0.00116)
URRC_2_0_IN = La13511Case(
    case_id="URRc-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 57"),
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
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.633380,
        sources=(Citation("SoodForsterParsons2003", "p. 85"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 1.933422},
    ),
    notes="Larger Σ_1f / Σ_21s than URRb; same upscatter structure.",
)


# Case 62 — URRd-2-0-IN: ISLC base material, 2G isotropic, no upscatter, ν=1.004 (slightly unphysical per Sood)
URRD_2_0_IN = La13511Case(
    case_id="URRd-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 62"),
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
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.034970,
        sources=(Citation("SoodForsterParsons2003", "p. 87"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 2.023344},
    ),
    notes="ISLC (Infinite Slab Lattice Cell) base XS. Sood uses ν_fast=1.004 'to stress code verification' (Sood 2003 p. 79 — i.e. unphysical but algebraically valid). Σ_22s = 0 (no fast self-scatter).",
)


# Case 67 — UD2O-2-0-IN: U-D2O reactor, 2G isotropic, no upscatter — k_inf is just barely critical (1.000221)
UD2O_2_0_IN = La13511Case(
    case_id="UD2O-2-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 67"),
    description="U-D2O reactor bare infinite medium, 2G isotropic",
    materials={0: _mix_ud2o_2g()},
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.000221,
        sources=(Citation("SoodForsterParsons2003", "p. 88"),),
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 26.822093},
    ),
    notes="Heavy-water reactor; k_inf = 1.000221 is essentially at the infinite-medium critical threshold. φ_fast/φ_slow = 26.822 (very slow-flux dominated due to D2O moderation).",
)


# ═══════════════════════════════════════════════════════════════════
# 3G k_inf infinite-medium case (Sood 2003 Tables 56-58)
# ═══════════════════════════════════════════════════════════════════


# Case 74 — URR-3-0-IN: 3-group URR, no upscatter
# Sood ordering: g3=fast, g2=mid, g1=slow. ORPHEUS: g=0 fast, g=1 mid, g=2 slow.
URR_3_0_IN = La13511Case(
    case_id="URR-3-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 74"),
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
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.60,
        sources=(Citation("SoodForsterParsons2003", "p. 91"), Citation("ODell1998QA3Group")),
        flux_ratio_groupwise={0: 1.0, 1: 0.480, 2: 0.150},
    ),
    notes="Sood Tables 56/57/58 take their cross sections from O'Dell (Ref. [47]), chosen so that f_23 = 4 and f_13 = 15 give k_inf = 1.60 exactly with φ_2/φ_3 = 0.480 and φ_1/φ_3 = 0.150 (Eqs (A.43)-(A.48)). All match to machine precision.",
)


# ═══════════════════════════════════════════════════════════════════
# 6G k_inf infinite-medium case (Sood 2003 Tables 59-64)
# ═══════════════════════════════════════════════════════════════════
#
# URR-6-0-IN is built from 2 coupled 3-group blocks (groups 6,5,4
# mirror groups 1,2,3). Decoupled in scattering, coupled only via χ.
# Sood guarantees same k_inf and flux ratios as URR-3-0-IN.

URR_6_0_IN = La13511Case(
    case_id="URR-6-0-IN",
    problem=Citation("SoodForsterParsons2003", "problem 75"),
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
    geometry_kind="infinite",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.60,
        sources=(Citation("SoodForsterParsons2003", "p. 92"), Citation("ODell1998QAUpscatter")),
        flux_ratio_groupwise={
            0: 1.0,    # Sood g6 fast (ORPHEUS 0)
            1: 0.480,  # Sood g5 mid  (= φ_5/φ_6)
            2: 0.150,  # Sood g4 slow (= φ_4/φ_6)
            3: 0.150,  # Sood g3 slow (mirror of 4)
            4: 0.480,  # Sood g2 mid  (mirror of 5)
            5: 1.0,    # Sood g1 fast (mirror of 6)
        },
    ),
    notes="6G == 2 coupled URR-3-0-IN blocks. Same k_inf=1.60. Has thermal-upscatter pattern in the bottom 3 groups (Σ_21s=0.171, Σ_31s=0.033, Σ_32s=0.275). compute_kinf_mg is the only Branch-2 entry that handles this case — Eq (A.12) specialisation cannot.",
)


# ═══════════════════════════════════════════════════════════════════
# 1G slab F_N bare-critical cases (Sood 2003 Tables 3, 4, 10, 14)
# ═══════════════════════════════════════════════════════════════════


# Case 2 — PUa-1-0-SL: Pu-239 (a), c=1.50 slab
PUA_1_0_SL = La13511Case(
    case_id="PUa-1-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 2"),
    description="Pu-239 (a) bare slab, 1G isotropic, c=1.50",
    materials={0: _mix_pua_1g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 3"), Citation("KornreichGanapol1997")), critical_dimension_mfp=0.605055),
    notes="F_N solver at N=12 reaches err ≤ 2e-6 vs Sood truth (well within the 1e-5 tolerance). Highest c in the bare 1G slab family (c=1.50).",
)


# Case 6 — PUb-1-0-SL: Pu-239 (b), c=1.40 slab
PUB_1_0_SL = La13511Case(
    case_id="PUb-1-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 6"),
    description="Pu-239 (b) bare slab, 1G isotropic, c=1.40",
    materials={0: _mix_pub_1g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 4"), Citation("SoodForsterParsons2003", "Table 5"), Citation("KaperLindemanLeaf1974")),
        flux_ratios={
            0.25: 0.9701734,
            0.50: 0.8810540,
            0.75: 0.7318131,
            1.00: 0.4902592,
        },
        critical_dimension_mfp=0.73660355,
    ),
    notes="Slab F_N at N=12 reaches err ≤ 3e-6 on a_c.",
)


# Case 22 — UD2O-1-0-SL: U-D2O reactor, c=1.02 slab
UD2O_1_0_SL = La13511Case(
    case_id="UD2O-1-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 22"),
    description="U-D2O reactor bare slab, 1G isotropic, c=1.02",
    materials={0: _mix_ud2o_1g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 14"), Citation("SoodForsterParsons2003", "Table 15"), Citation("KaperLindemanLeaf1974")),
        flux_ratios={
            0.25: 0.93945236,
            0.50: 0.76504084,
            0.75: 0.49690627,
            1.00: 0.13893858,
        },
        critical_dimension_mfp=5.6655054562,
    ),
    notes="Lowest c in the bare 1G slab family. Slab F_N at N=12 reaches err ≤ 2e-6 on a_c. NOTE: F_N at N≥14 fails for low c (determinant scan loses bracket); use N=12 for this case.",
)


# ═══════════════════════════════════════════════════════════════════
# 1G sphere F_N bare-critical cases (Sood 2003 Tables 4, 10, 14)
# ═══════════════════════════════════════════════════════════════════


# Case 8 — PUb-1-0-SP: Pu-239 (b), c=1.40 sphere
PUB_1_0_SP = La13511Case(
    case_id="PUb-1-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 8"),
    description="Pu-239 (b) bare sphere, 1G isotropic, c=1.40",
    materials={0: _mix_pub_1g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 4"), Citation("SoodForsterParsons2003", "Table 5"), Citation("KaperLindemanLeaf1974")),
        flux_ratios={
            0.25: 0.93538006,
            0.50: 0.75575352,
            0.75: 0.49884364,
            1.00: 0.19222603,
        },
        critical_dimension_mfp=1.9853434324,
    ),
    notes="Sphere F_N at N=10 reaches err ≤ 5e-8 on R_c (well within 1e-5).",
)


# Case 24 — UD2O-1-0-SP: U-D2O, c=1.02 sphere
UD2O_1_0_SP = La13511Case(
    case_id="UD2O-1-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 24"),
    description="U-D2O reactor bare sphere, 1G isotropic, c=1.02",
    materials={0: _mix_ud2o_1g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 14"), Citation("SoodForsterParsons2003", "Table 15"), Citation("KaperLindemanLeaf1974")),
        flux_ratios={
            0.25: 0.91063756,
            0.50: 0.67099621,
            0.75: 0.35561622,
            1.00: 0.04678614,
        },
        critical_dimension_mfp=12.027532098,
    ),
    notes="Sphere F_N at N=10 reaches err ≤ 4e-8 on R_c.",
)


# ═══════════════════════════════════════════════════════════════════
# 1G cylinder bare-critical cases — STUBS (B1 dispatch will activate)
# ═══════════════════════════════════════════════════════════════════
#
# Truth values from Sood 2003 Tables 4, 14. These will be activated when
# the B1 cylinder solver (Westfall-Metcalf 1973 singular eigenfunction
# expansion) is shipped. NO solver tests are added in Phase B3 for
# these cases — only the registry entries.


# Case 7 — PUb-1-0-CY: Pu-239 (b), c=1.40 cylinder
PUB_1_0_CY_STUB = La13511Case(
    case_id="PUb-1-0-CY",
    problem=Citation("SoodForsterParsons2003", "problem 7"),
    description="Pu-239 (b) bare cylinder, 1G isotropic, c=1.40 — STUB",
    materials={0: _mix_pub_1g()},
    geometry_kind="cylinder",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 4"), Citation("SoodForsterParsons2003", "Table 5"), Citation("Westfall1983"), Citation("WestfallMetcalf1972")),
        flux_ratios={
            0.50: 0.8093,
            1.00: 0.2926,
        },
        critical_dimension_mfp=1.396979,
    ),
    notes="STUB: solver activated by B1 dispatch (Westfall-Metcalf 1973 cylinder F_N). Sood publishes flux ratios only at r/r_c = 0.5 and 1.0 (Table 5) to 4 digits. Truth values verified.",
)


# Case 23 — UD2O-1-0-CY: U-D2O reactor, c=1.02 cylinder
UD2O_1_0_CY_STUB = La13511Case(
    case_id="UD2O-1-0-CY",
    problem=Citation("SoodForsterParsons2003", "problem 23"),
    description="U-D2O reactor bare cylinder, 1G isotropic, c=1.02 — STUB",
    materials={0: _mix_ud2o_1g()},
    geometry_kind="cylinder",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 14"), Citation("Westfall1983"), Citation("WestfallMetcalf1972")), critical_dimension_mfp=9.043255),
    notes="STUB: solver activated by B1 dispatch. No flux ratios published for this case in Sood's tables.",
)


# ═══════════════════════════════════════════════════════════════════
# 2G bare-critical STUBS (machinery deferred — needs Siewert-Thomas
# 1986 2G F_N or equivalent)
# ═══════════════════════════════════════════════════════════════════
#
# These cases have published truth values in Sood 2003 but require 2G
# F_N machinery that is NOT yet implemented in ORPHEUS. Registered as
# stubs so future implementations can pick them up directly.




PU_2_0_SL_STUB = La13511Case(
    case_id="PU-2-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 45"),
    description="Pu-239 bare slab, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_pu_2g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 29"), Citation("SiewertThomas1986"), Citation("ForsterMetcalf1969"), Citation("Forster1970")), critical_dimension_mfp=0.396469),
    notes="STUB: needs Siewert-Thomas 1986 2G F_N slab machinery (not yet implemented).",
)


PU_2_0_SP_STUB = La13511Case(
    case_id="PU-2-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 46"),
    description="Pu-239 bare sphere, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_pu_2g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 29"), Citation("SiewertThomas1986"), Citation("ForsterMetcalf1969"), Citation("Forster1970")), critical_dimension_mfp=1.15513),
    notes="STUB: needs Siewert-Thomas 1986 2G F_N sphere machinery. The slab and sphere F_N share the geometry-sign abstraction in fn_method.core; extending to 2G requires the matrix dispersion law (Λ matrix; Case eigenvalues are 2x2 matrix roots not scalars). High priority follow-on after B1 cylinder lands.",
)




U_2_0_SL_STUB = La13511Case(
    case_id="U-2-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 48"),
    description="U-235 bare slab, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_u_2g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 32"), Citation("SiewertThomas1986"), Citation("ForsterMetcalf1969"), Citation("Forster1970")), critical_dimension_mfp=0.649377),
    notes="STUB: needs 2G F_N slab machinery.",
)


U_2_0_SP_STUB = La13511Case(
    case_id="U-2-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 49"),
    description="U-235 bare sphere, 2G isotropic, no upscatter — STUB",
    materials={0: _mix_u_2g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 32"), Citation("SiewertThomas1986")), critical_dimension_mfp=1.70844),
    notes="STUB: needs 2G F_N sphere machinery.",
)




UAL_2_0_SL_STUB = La13511Case(
    case_id="UAL-2-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 51"),
    description="U-Al-Water bare slab, 2G isotropic — STUB",
    materials={0: _mix_ual_2g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 35"), Citation("SiewertThomas1986"), Citation("ForsterMetcalf1969"), Citation("Forster1970")), critical_dimension_mfp=2.09994),
    notes="STUB: needs 2G F_N slab machinery.",
)


UAL_2_0_SP_STUB = La13511Case(
    case_id="UAL-2-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 52"),
    description="U-Al-Water bare sphere, 2G isotropic — STUB",
    materials={0: _mix_ual_2g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 35"), Citation("SiewertThomas1986")), critical_dimension_mfp=4.73786),
    notes="STUB: needs 2G F_N sphere machinery.",
)




URRA_2_0_SL_STUB = La13511Case(
    case_id="URRa-2-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 54"),
    description="URR (a) bare slab, 2G isotropic — STUB",
    materials={0: _mix_urra_2g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        sources=(Citation("SoodForsterParsons2003", "Table 38"), Citation("SoodForsterParsons2003", "Table 39"), Citation("SiewertThomas1986"), Citation("ForsterMetcalf1969"), Citation("Forster1970")),
        flux_ratios={
            # Sood 2003 Table 39 — fast group, normalised to fast at center
            0.241394: 0.943363,
            0.502905: 0.761973,
            0.744300: 0.504012,
            1.0: 0.147598,
        },
        critical_dimension_mfp=4.97112,
    ),
    notes="STUB: needs 2G F_N slab machinery. Sood Table 39 gives 2G flux ratios at four spatial points (fast + slow). flux_ratios stored here is the FAST group; the slow-group ratio at the same points is (0.340124, 0.273056, 0.173845, 0.0212324).",
)


URRA_2_0_SP_STUB = La13511Case(
    case_id="URRa-2-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 55"),
    description="URR (a) bare sphere, 2G isotropic — STUB",
    materials={0: _mix_urra_2g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 38"), Citation("SiewertThomas1986")), critical_dimension_mfp=10.5441),
    notes="STUB: needs 2G F_N sphere machinery.",
)




UD2O_2_0_SL_STUB = La13511Case(
    case_id="UD2O-2-0-SL",
    problem=Citation("SoodForsterParsons2003", "problem 68"),
    description="U-D2O reactor bare slab, 2G isotropic — STUB",
    materials={0: _mix_ud2o_2g()},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 48"), Citation("SiewertThomas1986"), Citation("ForsterMetcalf1969"), Citation("Forster1970")), critical_dimension_mfp=284.367),
    notes="STUB: needs 2G F_N slab machinery. Critical dimension is VERY LARGE (284 mfp — barely-supercritical heavy-water reactor); high N_F may be needed.",
)


UD2O_2_0_SP_STUB = La13511Case(
    case_id="UD2O-2-0-SP",
    problem=Citation("SoodForsterParsons2003", "problem 69"),
    description="U-D2O reactor bare sphere, 2G isotropic — STUB",
    materials={0: _mix_ud2o_2g()},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 48"), Citation("SiewertThomas1986")), critical_dimension_mfp=569.43),
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
    problem=Citation("SoodForsterParsons2003", "problem 32"),
    description="Pu-239 (a) bare slab, 1G P_1 anisotropic (forward), c=1.40",
    materials={0: _mix_1g_anisotropic(
        sigma_t=1.0, sigma_c=0.0, sigma_f=0.266667, nu=2.5,
        sigma_s_self=0.733333, sigma_s1_self=0.20,
    )},
    geometry_kind="slab",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 22"), Citation("LathropLeonard1965")), critical_dimension_mfp=0.77032),
    notes="Σ_s1=0.20. Carlvik-Galerkin uses μ̄_eff = 0.20/1.40 = 0.142857.",
)

PUB_1_1_SL = La13511Case(
    case_id="PUb-1-1-SL",
    problem=Citation("SoodForsterParsons2003", "problem 34"),
    description="Pu-239 (b) bare slab, 1G P_1 anisotropic (strong forward), c=1.40",
    materials={0: _mix_1g_anisotropic(
        sigma_t=1.0, sigma_c=0.0, sigma_f=0.266667, nu=2.5,
        sigma_s_self=0.733333, sigma_s1_self=0.333333,
    )},
    geometry_kind="slab",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 22"), Citation("LathropLeonard1965")), critical_dimension_mfp=0.79606),
    notes="Σ_s1=0.333333 (negative scattering for μ near -1). μ̄_eff = 0.238095.",
)

UD2OA_1_1_SP = La13511Case(
    case_id="UD2Oa-1-1-SP",
    problem=Citation("SoodForsterParsons2003", "problem 39"),
    description="U-D2O (a) bare sphere, 1G P_1 anisotropic, c=1.0308381",
    materials={0: _mix_1g_anisotropic(
        sigma_t=0.54628, sigma_c=0.027314, sigma_f=0.054628, nu=1.808381,
        sigma_s_self=0.464338, sigma_s1_self=0.056312624,
    )},
    geometry_kind="sphere",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 26"), Citation("DahlSjostrand1979")), critical_dimension_mfp=10.0),
    notes="μ̄_eff = 0.10 — matches DS Table I row d=20, μ̄=0.10.",
)

UD2OB_1_1_SP = La13511Case(
    case_id="UD2Ob-1-1-SP",
    problem=Citation("SoodForsterParsons2003", "problem 41"),
    description="U-D2O (b) bare sphere, 1G P_1 anisotropic, c=1.0341086",
    materials={0: _mix_1g_anisotropic(
        sigma_t=0.54628, sigma_c=0.027314, sigma_f=0.054628, nu=1.841086,
        sigma_s_self=0.464338, sigma_s1_self=0.112982569,
    )},
    geometry_kind="sphere",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 26"), Citation("DahlSjostrand1979")), critical_dimension_mfp=10.0),
    notes="μ̄_eff = 0.20 — matches DS Table I row d=20, μ̄=0.20.",
)

UD2OC_1_1_SP = La13511Case(
    case_id="UD2Oc-1-1-SP",
    problem=Citation("SoodForsterParsons2003", "problem 43"),
    description="U-D2O (c) bare sphere, 1G P_1 anisotropic (back-peaked!), c=1.01964",
    materials={0: _mix_1g_anisotropic(
        sigma_t=0.54628, sigma_c=0.027314, sigma_f=0.054628, nu=1.6964,
        sigma_s_self=0.464338, sigma_s1_self=-0.27850447,
    )},
    geometry_kind="sphere",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("SoodForsterParsons2003", "Table 26"), Citation("SahniDahlSjostrand1995")), critical_dimension_mfp=10.0),
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

Coverage: 47 cases, the Phase A first slice (5) and the Phase B3 wide
enumeration (20 k_inf, 5 one-group bare slabs and spheres, 5 of them with
P_1 scattering, 12 stubs). By kind:

* 22 infinite-medium ``k_inf`` cases (12 one-group, 8 two-group, one
  three-group, one six-group).
* 12 one-group bare slabs and spheres.
* 3 one-group bare cylinders.
* 10 two-group bare slabs and spheres, which wait for the two-group
  F_N solver (Siewert-Thomas 1986).
"""

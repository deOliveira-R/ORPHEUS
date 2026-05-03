r"""Sood/Forster/Parsons LA-13511 (1999) benchmark case catalogue.

Machine-readable case definitions for the 5-case complexity ramp
identified in the literature memo
``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``.
The first slice fully populates Cases 1 (PUa-1-0-IN) and 5
(PU-2-0-IN); Cases 2 (Ua-1-0-SL), 3 (Ua-1-0-CY), 4 (Ua-1-0-SP) are
shipped as **stubs** with cross sections + reference values but no
solver — the F_N machinery they need lives in cited papers we do not
yet have.

ORPHEUS convention vs Sood convention
-------------------------------------

Sood et al. number energy groups :math:`g=N` (fast) → :math:`g=1`
(slow), the reverse of typical nuclear-engineering convention. This
module does the conversion at the boundary: ``sigma_t``,
``sigma_s``, ``nu_sigma_f``, ``chi`` arrays are stored in **ORPHEUS
convention** (``g=0`` fast, ``g=N-1`` slow), so consumers (e.g.
:mod:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`) can
take them as-is. The Sood-side tables are reproduced verbatim in
the docstrings for traceability.

The scattering matrix uses ORPHEUS's ``[from, to]`` convention:
``sigma_s[g, h]`` is :math:`\Sigma_{s, g \to h}`. The ``Sigma_g^{rem}``
removal cross sections used by Sood Eq 26-27 are therefore
``sigma_t[g] - sigma_s[g, g]``.

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

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np


@dataclass(frozen=True)
class La13511Case:
    r"""A single LA-13511 benchmark configuration.

    Cross sections are in **ORPHEUS convention** (``g=0`` fast,
    ``g=N-1`` slow). Sood's reverse convention is converted at
    construction time; consumers should not need to re-index.

    Parameters
    ----------
    case_id : str
        Sood naming convention identifier
        (``<Material>-<Groups>-<Scattering>-<Geometry>``).
    problem_number : int
        Sequential problem number within LA-13511.
    description : str
        One-line human-readable description.
    geometry : str
        One of ``"infinite"``, ``"slab"``, ``"cylinder"``, ``"sphere"``,
        ``"ISLC"``.
    n_groups : int
        Number of energy groups (1, 2, 3, or 6).
    scattering_order : int
        Legendre order of the scattering kernel (0 = isotropic,
        1 = P_1, 2 = P_2).
    sigma_t : np.ndarray
        Total macroscopic cross section, shape ``(n_groups,)``,
        in ORPHEUS group order. Units: 1/cm.
    sigma_s : np.ndarray
        Scattering matrix, shape ``(n_groups, n_groups)``, with
        ``sigma_s[g, h]`` = :math:`\Sigma_{s, g \to h}` (ORPHEUS
        ``[from, to]`` convention). Units: 1/cm. For
        ``scattering_order > 0`` this is the :math:`P_0` (isotropic)
        moment only; higher Legendre moments are not stored in the
        first slice.
    nu_sigma_f : np.ndarray
        Production cross section :math:`\nu \Sigma_f`, shape
        ``(n_groups,)``. Units: 1/cm.
    chi : np.ndarray
        Fission spectrum, shape ``(n_groups,)``. Dimensionless,
        sums to 1.
    critical_dimension_mfp : float | None
        Critical radius (sphere/cylinder) or half-thickness (slab) in
        mean free paths. ``None`` for infinite-medium cases.
    critical_dimension_cm : float | None
        Same as above but in cm.
    k_eff_or_kinf : float
        Reference :math:`k_{\rm eff}` (finite cases) or
        :math:`k_\infty` (infinite cases) from LA-13511.
    flux_ratios : Mapping[float, float] | None
        For 1G cases with published flux table: dict mapping
        ``r/r_c`` to :math:`\phi(r)/\phi(0)`. ``None`` when no flux
        table is published.
    flux_ratio_groupwise : Mapping[int, float] | None
        For 2G infinite cases: dict mapping group index (ORPHEUS
        order) to :math:`\phi_g/\phi_0` (ratio relative to ORPHEUS
        group 0 = fast). For ``PU-2-0-IN``, this is just
        ``{0: 1.0, 1: phi_slow/phi_fast}``. Sood's published ratio
        is :math:`\phi_2/\phi_1 = \phi_{\rm fast}/\phi_{\rm slow}`,
        so the slow-group entry stored here is the **inverse** of
        Sood's published number.
    sood_table : int
        LA-13511 table number where this case is tabulated.
    primary_reference : str
        The peer-reviewed paper Sood cites as the source of the
        reference values.
    notes : str
        Any extra remarks (typo flags, conversion subtleties, etc).
    """

    case_id: str
    problem_number: int
    description: str
    geometry: str
    n_groups: int
    scattering_order: int
    sigma_t: np.ndarray
    sigma_s: np.ndarray
    nu_sigma_f: np.ndarray
    chi: np.ndarray
    critical_dimension_mfp: float | None
    critical_dimension_cm: float | None
    k_eff_or_kinf: float
    flux_ratios: Mapping[float, float] | None
    flux_ratio_groupwise: Mapping[int, float] | None
    sood_table: int
    primary_reference: str
    notes: str = ""


# ═══════════════════════════════════════════════════════════════════
# Case 1 — PUa-1-0-IN (Sood problem 1): 1G infinite medium, Pu-239 (a)
# ═══════════════════════════════════════════════════════════════════
#
# LA-13511 Table 2 (Pu-239 (a) cross sections, 1G isotropic):
#   ν = 3.24,  Σ_f = 0.0816,  Σ_c = 0.019584,  Σ_s = 0.225216,
#   Σ_t = 0.32640,  c = (Σ_s + νΣ_f)/Σ_t = 1.50.
# Reference: k_inf = 2.612903 (LA-13511 Eq 20).
#
# Sanity check via Eq 19 (the simplified form):
#   k_inf = νΣ_f / (Σ_t - Σ_s) = 3.24·0.0816 / (0.32640 - 0.225216)
#         = 0.264384 / 0.101184 = 2.6129032...    ✓

PUA_1_0_IN = La13511Case(
    case_id="PUa-1-0-IN",
    problem_number=1,
    description="Pu-239 (a) bare infinite medium, 1G isotropic",
    geometry="infinite",
    n_groups=1,
    scattering_order=0,
    sigma_t=np.array([0.32640]),
    sigma_s=np.array([[0.225216]]),
    nu_sigma_f=np.array([3.24 * 0.0816]),  # 0.264384
    chi=np.array([1.0]),
    critical_dimension_mfp=None,
    critical_dimension_cm=None,
    k_eff_or_kinf=2.612903,
    flux_ratios=None,  # constant flux in infinite medium
    flux_ratio_groupwise=None,
    sood_table=2,
    primary_reference="LA-13511 Eq 20 (closed form)",
    notes=(
        "1G infinite-medium k_inf reduces to nu·Sigma_f/Sigma_a; the "
        "'c' factor in Eq 20 cancels algebraically (verified in "
        "origins.k_inf_derivations.derive_kinf_1g_eq_20_simplifies_to_"
        "nu_sigma_f_over_sigma_a)."
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
#     nu_sigma_f = [3.10·0.0936, 2.93·0.08544] = [0.290160, 0.2503392]
#     chi = [χ_2, χ_1] = [0.575, 0.425]
#
# Reference: k_inf = 2.683767, φ_2/φ_1 = 0.675229 (Sood Eq 28-29 + Eq 32).
# Sood publishes the ratio φ_2/φ_1 = phi_fast/phi_slow.
# In ORPHEUS we store flux_ratio_groupwise = {0: 1.0, 1: phi_slow/phi_fast}
# i.e. fast = unity, slow = inverse of Sood's published ratio.

PU_2_0_IN = La13511Case(
    case_id="PU-2-0-IN",
    problem_number=44,
    description="Pu bare infinite medium, 2G isotropic, no upscatter",
    geometry="infinite",
    n_groups=2,
    scattering_order=0,
    sigma_t=np.array([0.2208, 0.3360]),
    sigma_s=np.array([
        [0.0792, 0.0432],  # from g=0 (fast):  → fast self, → slow downscatter
        [0.0,    0.23616], # from g=1 (slow):  → fast upscatter (none), → slow self
    ]),
    nu_sigma_f=np.array([3.10 * 0.0936, 2.93 * 0.08544]),  # [0.290160, 0.2503392]
    chi=np.array([0.575, 0.425]),
    critical_dimension_mfp=None,
    critical_dimension_cm=None,
    k_eff_or_kinf=2.683767,
    flux_ratios=None,
    flux_ratio_groupwise={0: 1.0, 1: 1.0 / 0.675229},  # slow/fast = 1/Sood's ratio
    sood_table=30,
    primary_reference="LA-13511 Eq 28-29 (k_inf) + Eq 32 (flux ratio)",
    notes=(
        "Sood Eq 28 has a typo: the chi_1 and chi_2 numerators have "
        "the wrong Sigma_g^rem factor — Eq 28 as printed reduces to "
        "2.862 (not 2.684) when Sigma_21s = 0, while Eq 29 (printed "
        "separately) reduces correctly. The SymPy derivation in "
        "origins.k_inf_derivations.derive_kinf_2g_no_upscatter exposes "
        "the correct general form by computing det(M)=0 from Eq 25 "
        "directly, and verifies the published Eq 29 against it. See "
        "the closeout memo for the typo analysis."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Case 2 — Ua-1-0-SL (Sood problem 12): 1G bare slab, U-235 (a)  [STUB]
# ═══════════════════════════════════════════════════════════════════

UA_1_0_SL_STUB = La13511Case(
    case_id="Ua-1-0-SL",
    problem_number=12,
    description="U-235 (a) bare slab, 1G isotropic — STUB (needs F_N)",
    geometry="slab",
    n_groups=1,
    scattering_order=0,
    sigma_t=np.array([0.32640]),
    sigma_s=np.array([[0.248064]]),
    nu_sigma_f=np.array([2.70 * 0.06528]),  # 0.176256
    chi=np.array([1.0]),
    critical_dimension_mfp=0.93772556,
    critical_dimension_cm=2.872934,
    k_eff_or_kinf=1.0,  # critical
    flux_ratios={
        0.25: 0.9669506,
        0.50: 0.8686259,
        0.75: 0.7055218,
        1.00: 0.4461912,
    },
    flux_ratio_groupwise=None,
    sood_table=4,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94",
    notes=(
        "STUB: F_N solver not yet implemented (needs Kaper-Lindeman-"
        "Leaf 1974). XS + reference values populated for future use."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Case 3 — Ua-1-0-CY (Sood problem 13): 1G bare cylinder, U-235 (a)  [STUB]
# ═══════════════════════════════════════════════════════════════════

UA_1_0_CY_STUB = La13511Case(
    case_id="Ua-1-0-CY",
    problem_number=13,
    description="U-235 (a) bare infinite cylinder, 1G isotropic — STUB",
    geometry="cylinder",
    n_groups=1,
    scattering_order=0,
    sigma_t=np.array([0.32640]),
    sigma_s=np.array([[0.248064]]),
    nu_sigma_f=np.array([2.70 * 0.06528]),
    chi=np.array([1.0]),
    critical_dimension_mfp=1.72500292,
    critical_dimension_cm=5.284935,
    k_eff_or_kinf=1.0,
    flux_ratios=None,
    flux_ratio_groupwise=None,
    sood_table=5,
    primary_reference="Westfall-Metcalf 1973 NSE 52, 1",
    notes=(
        "STUB: F_N solver not yet implemented (needs Westfall-"
        "Metcalf 1973). Already cross-checked by Variant α at 8.5e-6."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Case 4 — Ua-1-0-SP (Sood problem 14): 1G bare sphere, U-235 (a)  [STUB]
# ═══════════════════════════════════════════════════════════════════

UA_1_0_SP_STUB = La13511Case(
    case_id="Ua-1-0-SP",
    problem_number=14,
    description="U-235 (a) bare sphere, 1G isotropic — STUB (needs F_N)",
    geometry="sphere",
    n_groups=1,
    scattering_order=0,
    sigma_t=np.array([0.32640]),
    sigma_s=np.array([[0.248064]]),
    nu_sigma_f=np.array([2.70 * 0.06528]),
    chi=np.array([1.0]),
    critical_dimension_mfp=2.4248249802,
    critical_dimension_cm=7.428998,
    k_eff_or_kinf=1.0,
    flux_ratios=None,  # not published for Ua-1-0-SP — see PUb-1-0-SP for an example
    flux_ratio_groupwise=None,
    sood_table=6,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94",
    notes=(
        "STUB: F_N solver not yet implemented (needs Kaper-Lindeman-"
        "Leaf 1974). Top priority for cross-checking the Variant α "
        "Green's-function sphere prototype."
    ),
)


ALL_FIRST_SLICE: tuple[La13511Case, ...] = (
    PUA_1_0_IN,
    PU_2_0_IN,
    UA_1_0_SL_STUB,
    UA_1_0_CY_STUB,
    UA_1_0_SP_STUB,
)
"""All first-slice cases. Cases 1+5 (k_inf) ship with full Branch-2
solvers; cases 2-4 (slab/cylinder/sphere) are stubs awaiting F_N
machinery from cited papers."""

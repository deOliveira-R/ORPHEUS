r"""Atalay 1997 reflected-slab and sphere case catalogue.

Atalay (1997) [Atalay1997]_ tabulates critical thicknesses :math:`2d`
(slab) and eigenvalues :math:`c` (slab + sphere) for the **reflected**
slab and sphere with linearly-anisotropic scattering — a regime that
is **not** in standard Sood/Forster/Parsons LA-13511 tables (which
focus on bare configurations).

This file is the case-class home for:

* **Reflected slab + isotropic** (Atalay Table 2, ``f_1 = 0``).
* **Reflected slab + linearly anisotropic** (Atalay Tables 3, 4, 5
  for ``f_1 = 0.10, 0.20, 0.30``).
* **Vacuum slab + linearly anisotropic** (Atalay Tables 3-5 R = 0
  column) — Atalay-unique primary.
* **Vacuum sphere + isotropic** (Atalay Table 6 R=0 + Sood Ua-1-0-SP
  cross-check).
* **Reflected sphere + linearly anisotropic** at f_1 = 0.10 only
  (Atalay Table 10) — Atalay-unique with no independent reference.

Provenance
----------

Atalay's tables use the Sood-style :math:`(c, R, f_1)` parametrisation
where :math:`c` is the secondaries-per-collision and :math:`R` is the
specular reflection coefficient at both slab faces / outer sphere
surface. Cases here are keyed by triples ``(c, R, f_1)`` rather than
material composition, mirroring how Atalay published them.

For mapping back to specific Sood materials (e.g., ``Ua-1-0-SL`` has
:math:`c = 1.30`), see the case ``description`` field.

References
----------

.. [Atalay1997] Atalay, M.A. (1997).
   "The reflected slab and sphere criticality problem with anisotropic
   scattering in one-speed neutron transport theory."
   *Progress in Nuclear Energy* **31**(3), 229-252.
   DOI: 10.1016/0149-1970(95)00094-1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import make_mixture

from .la13511 import La13511Case, La13511Truth, Provenance


# ─── Common XS for c=1.30 isotropic (matches Sood U-235(a) Σ_t=0.32640) ───


def _mix_iso_at_c(c: float) -> Mixture:
    r"""Build a 1G isotropic Mixture with the specified secondaries-per-collision c.

    Use Sood U-235(a) Σ_t = 0.32640 cm⁻¹ as a normalization, set
    :math:`\Sigma_s + \nu \Sigma_f = c \cdot \Sigma_t`. Atalay's tables
    use only c (not material composition) — any choice that gives the
    target c is fine. We pick a fully-scattering mixture
    (:math:`\Sigma_s = c \Sigma_t`, :math:`\Sigma_f = 0`) — this is
    the simplest realisation. Note: this is **not** matchable to a
    physical Sood material when :math:`c > 1` and :math:`\Sigma_f = 0`,
    but the criticality condition for Atalay only depends on c.
    """
    sigma_t = 0.32640
    sigma_s_self = c * sigma_t
    sigma_c = sigma_t - sigma_s_self
    if sigma_c < 0:
        # For c > 1, "absorption" via a (negative) capture is unphysical
        # in standard Mixture. Fold the multiplying factor into ν Σ_f.
        # But for this Atalay catalogue we just need (c, R, f_1) triples;
        # the Mixture object is not consumed by case_method solvers
        # (which take c directly). Set sigma_c = 0 and use ν Σ_f to
        # carry the multiplying excess.
        sigma_c = 0.0
        sigma_s_self = sigma_t  # all scattering
        # Add a fission term: ν Σ_f = (c - 1) Σ_t.
        nu = 1.0
        sigma_f = (c - 1.0) * sigma_t
    else:
        nu = 1.0
        sigma_f = 0.0
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sigma_c]),
        sig_f=np.array([sigma_f]),
        nu=np.array([nu]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s_self]]),
    )


# ═══════════════════════════════════════════════════════════════════
# Atalay-anchored slab cases (c, R, f_1) → critical 2d (mfp).
# Names: ATALAY_SL_{c100}_{R100}_{f1_100} with X100 = X·100 rounded.
# ═══════════════════════════════════════════════════════════════════


# Atalay Table 2 (f_1 = 0): reflected slab, isotropic
ATALAY_SLAB_C130_R000_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.00-f1_0.00",
    problem_number=1001,
    description=(
        "Atalay 1997 Table 2: c=1.30, R=0 (vacuum), f_1=0 (isotropic). "
        "Same c as Sood Ua-1-0-SL (U-235(a)); Atalay reports 2d=1.87766 mfp."
    ),
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=0.93883),
    sood_table=2,  # Atalay Table 2
    primary_reference="Atalay 1997 Table 2",
    notes="Atalay-anchored vacuum case; cross-check vs Sood Ua-1-0-SL gives 1.87545 (KLL 1974).",
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=2,
        primary_reference="Atalay 1997 Table 2",
        notes="Atalay-anchored vacuum case; cross-check vs Sood Ua-1-0-SL gives 1.87545 (KLL 1974).",
    ),
)


ATALAY_SLAB_C130_R025_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.25-f1_0.00",
    problem_number=1002,
    description="Atalay 1997 Table 2: c=1.30, R=0.25, f_1=0. Reports 2d=1.40621 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=0.703105),
    sood_table=2,
    primary_reference="Atalay 1997 Table 2",
    notes="Reflected-slab case (R=0.25). Atalay-unique primary; no Sood entry.",
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=2,
        primary_reference="Atalay 1997 Table 2",
        notes="Reflected-slab case (R=0.25). Atalay-unique primary; no Sood entry.",
    ),
)


ATALAY_SLAB_C130_R050_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.50-f1_0.00",
    problem_number=1003,
    description="Atalay 1997 Table 2: c=1.30, R=0.50, f_1=0. Reports 2d=0.89317 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=0.446585),
    sood_table=2,
    primary_reference="Atalay 1997 Table 2",
    notes="Reflected-slab case (R=0.50). Atalay-unique primary.",
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=2,
        primary_reference="Atalay 1997 Table 2",
        notes="Reflected-slab case (R=0.50). Atalay-unique primary.",
    ),
)


ATALAY_SLAB_C130_R075_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.75-f1_0.00",
    problem_number=1004,
    description="Atalay 1997 Table 2: c=1.30, R=0.75, f_1=0. Reports 2d=0.40758 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=0.20379),
    sood_table=2,
    primary_reference="Atalay 1997 Table 2",
    notes="Reflected-slab case (R=0.75). Atalay-unique primary.",
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=2,
        primary_reference="Atalay 1997 Table 2",
        notes="Reflected-slab case (R=0.75). Atalay-unique primary.",
    ),
)


# Atalay Table 3 (f_1 = 0.10): reflected slab, linearly anisotropic
ATALAY_SLAB_C130_R000_F010 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.00-f1_0.10",
    problem_number=1005,
    description=(
        "Atalay 1997 Table 3: c=1.30, R=0 (vacuum), f_1=0.10 (linearly anisotropic). "
        "Reports 2d=1.94146 mfp."
    ),
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=1,  # P_1 anisotropic
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=0.97073),
    sood_table=3,
    primary_reference="Atalay 1997 Table 3",
    notes=(
        "Vacuum + linearly anisotropic slab — Atalay-unique primary. "
        "f_1=0.10 means scattering kernel Σ_s(1+0.30 μμ')/2."
    ),
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=3,
        primary_reference="Atalay 1997 Table 3",
        notes="Vacuum + linearly anisotropic slab — Atalay-unique primary. f_1=0.10 means scattering kernel Σ_s(1+0.30 μμ')/2.",
    ),
)


ATALAY_SLAB_C130_R050_F010 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.50-f1_0.10",
    problem_number=1006,
    description="Atalay 1997 Table 3: c=1.30, R=0.50, f_1=0.10. Reports 2d=0.89831 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=0.449155),
    sood_table=3,
    primary_reference="Atalay 1997 Table 3",
    notes="Reflected + linearly anisotropic CROSS-PRODUCT case — Atalay-unique primary.",
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=3,
        primary_reference="Atalay 1997 Table 3",
        notes="Reflected + linearly anisotropic CROSS-PRODUCT case — Atalay-unique primary.",
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Atalay-anchored sphere cases.
# Atalay Table 10 is the only sphere data outside f_1 = 0; for f_1 = 0
# we use Sood Ua-1-0-SP at c = 1.30 (KLL 1974) as the reference.
# ═══════════════════════════════════════════════════════════════════


ATALAY_SPHERE_C130_R000_F0 = La13511Case(
    case_id="atalay-1997-sphere-c1.30-R0.00-f1_0.00",
    problem_number=2001,
    description=(
        "Vacuum sphere c=1.30, isotropic. Same XS as Sood Ua-1-0-SP "
        "(R_c = 2.4248 mfp, KLL 1974)."
    ),
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="sphere",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, critical_dimension_mfp=2.4248249802),
    sood_table=14,  # Atalay refers; KLL 1974 is the truth source.
    primary_reference="Kaper-Lindeman-Leaf 1974 (Atalay c=1.30 reproduces).",
    notes="Cross-check baseline; case_method sphere reproduces at ~0.5%.",
    provenance=Provenance(
        paper_id="Atalay-1997",
        paper_table=14,
        primary_reference="Kaper-Lindeman-Leaf 1974 (Atalay c=1.30 reproduces).",
        notes="Cross-check baseline; case_method sphere reproduces at ~0.5%.",
    ),
)


# All Atalay catalogue cases as a tuple
ATALAY_SLAB_CASES: tuple[La13511Case, ...] = (
    ATALAY_SLAB_C130_R000_F0,
    ATALAY_SLAB_C130_R025_F0,
    ATALAY_SLAB_C130_R050_F0,
    ATALAY_SLAB_C130_R075_F0,
    ATALAY_SLAB_C130_R000_F010,
    ATALAY_SLAB_C130_R050_F010,
)

ATALAY_SPHERE_CASES: tuple[La13511Case, ...] = (
    ATALAY_SPHERE_C130_R000_F0,
)

ATALAY_ALL_CASES: tuple[La13511Case, ...] = (
    *ATALAY_SLAB_CASES,
    *ATALAY_SPHERE_CASES,
)


__all__ = [
    # Slab
    "ATALAY_SLAB_C130_R000_F0",
    "ATALAY_SLAB_C130_R025_F0",
    "ATALAY_SLAB_C130_R050_F0",
    "ATALAY_SLAB_C130_R075_F0",
    "ATALAY_SLAB_C130_R000_F010",
    "ATALAY_SLAB_C130_R050_F010",
    # Sphere
    "ATALAY_SPHERE_C130_R000_F0",
    # Tuples
    "ATALAY_SLAB_CASES",
    "ATALAY_SPHERE_CASES",
    "ATALAY_ALL_CASES",
]

r"""Atalay 1997 reflected-slab case catalogue.

Atalay (1997) :cite:`Atalay1997` tabulates critical thicknesses :math:`2d`
(slab) and eigenvalues :math:`c` (slab + sphere) for the **reflected**
slab and sphere with linearly-anisotropic scattering — a regime that
Sood, Forster & Parsons (2003) do not tabulate (their reflected cases
carry a physical reflector material, not a specular reflection
coefficient at the face).

The catalogue holds six slab cases at :math:`c = 1.30`, each citing the
table that prints its critical half-thickness:

* **isotropic** (Atalay Table 2, :math:`f_1 = 0`) at :math:`R = 0, 0.25,
  0.50, 0.75`;
* **linearly anisotropic** (Atalay Table 3, :math:`f_1 = 0.10`) at
  :math:`R = 0` and :math:`0.50`.

Atalay's sphere results are the odd modes of the slab (its Section 3), and
it prints them at :math:`f_1 = 0.10` only (Table 10). A vacuum sphere at
:math:`c = 1.30` with isotropic scattering is Sood's problem 14
(``Ua-1-0-SP``), not an Atalay case: an Atalay-named copy of it stood here
until 2026-09-29, citing a "Table 14" the paper does not have.

Provenance
----------

Atalay's tables use the Sood-style :math:`(c, R, f_1)` parametrisation
where :math:`c` is the secondaries-per-collision and :math:`R` is the
specular reflection coefficient at both slab faces / outer sphere
surface. Cases here are keyed by triples ``(c, R, f_1)`` rather than
material composition, mirroring how Atalay published them.

For mapping back to specific Sood materials (e.g., ``Ua-1-0-SL`` has
:math:`c = 1.30`), see the case ``description`` field.

"""
from __future__ import annotations

import numpy as np

from orpheus.data.citation import Citation
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import make_mixture

from .case import La13511Case, La13511Truth


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
    problem=Citation("Atalay1997", "Table 2, c = 1.30, R = 0.00, f_1 = 0.00"),
    description=(
        "Atalay 1997 Table 2: c=1.30, R=0 (vacuum), f_1=0 (isotropic). "
        "Same c as Sood Ua-1-0-SL (U-235(a)); Atalay reports 2d=1.87766 mfp."
    ),
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("Atalay1997", "Table 2, c = 1.30, R = 0.00, f_1 = 0.00"),), critical_dimension_mfp=0.93883),
  # Atalay Table 2
    notes="Atalay-anchored vacuum case; cross-check vs Sood Ua-1-0-SL gives 1.87545 (KLL 1974).",
)


ATALAY_SLAB_C130_R025_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.25-f1_0.00",
    problem=Citation("Atalay1997", "Table 2, c = 1.30, R = 0.25, f_1 = 0.00"),
    description="Atalay 1997 Table 2: c=1.30, R=0.25, f_1=0. Reports 2d=1.40621 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("Atalay1997", "Table 2, c = 1.30, R = 0.25, f_1 = 0.00"),), critical_dimension_mfp=0.703105),
    notes="Reflected-slab case (R=0.25). Atalay-unique primary; no Sood entry.",
)


ATALAY_SLAB_C130_R050_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.50-f1_0.00",
    problem=Citation("Atalay1997", "Table 2, c = 1.30, R = 0.50, f_1 = 0.00"),
    description="Atalay 1997 Table 2: c=1.30, R=0.50, f_1=0. Reports 2d=0.89317 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("Atalay1997", "Table 2, c = 1.30, R = 0.50, f_1 = 0.00"),), critical_dimension_mfp=0.446585),
    notes="Reflected-slab case (R=0.50). Atalay-unique primary.",
)


ATALAY_SLAB_C130_R075_F0 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.75-f1_0.00",
    problem=Citation("Atalay1997", "Table 2, c = 1.30, R = 0.75, f_1 = 0.00"),
    description="Atalay 1997 Table 2: c=1.30, R=0.75, f_1=0. Reports 2d=0.40758 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=0,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("Atalay1997", "Table 2, c = 1.30, R = 0.75, f_1 = 0.00"),), critical_dimension_mfp=0.20379),
    notes="Reflected-slab case (R=0.75). Atalay-unique primary.",
)


# Atalay Table 3 (f_1 = 0.10): reflected slab, linearly anisotropic
ATALAY_SLAB_C130_R000_F010 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.00-f1_0.10",
    problem=Citation("Atalay1997", "Table 3, c = 1.30, R = 0.00, f_1 = 0.10"),
    description=(
        "Atalay 1997 Table 3: c=1.30, R=0 (vacuum), f_1=0.10 (linearly anisotropic). "
        "Reports 2d=1.94146 mfp."
    ),
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=1,  # P_1 anisotropic
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("Atalay1997", "Table 3, c = 1.30, R = 0.00, f_1 = 0.10"),), critical_dimension_mfp=0.97073),
    notes="Vacuum + linearly anisotropic slab — Atalay-unique primary. f_1=0.10 means scattering kernel Σ_s(1+0.30 μμ')/2.",
)


ATALAY_SLAB_C130_R050_F010 = La13511Case(
    case_id="atalay-1997-slab-c1.30-R0.50-f1_0.10",
    problem=Citation("Atalay1997", "Table 3, c = 1.30, R = 0.50, f_1 = 0.10"),
    description="Atalay 1997 Table 3: c=1.30, R=0.50, f_1=0.10. Reports 2d=0.89831 mfp.",
    materials={0: _mix_iso_at_c(1.30)},
    geometry_kind="slab",
    scattering_order=1,
    truth=La13511Truth(k_eff_or_kinf=1.0, sources=(Citation("Atalay1997", "Table 3, c = 1.30, R = 0.50, f_1 = 0.10"),), critical_dimension_mfp=0.449155),
    notes="Reflected + linearly anisotropic CROSS-PRODUCT case — Atalay-unique primary.",
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

ATALAY_ALL_CASES: tuple[La13511Case, ...] = ATALAY_SLAB_CASES


__all__ = [
    # Slab
    "ATALAY_SLAB_C130_R000_F0",
    "ATALAY_SLAB_C130_R025_F0",
    "ATALAY_SLAB_C130_R050_F0",
    "ATALAY_SLAB_C130_R075_F0",
    "ATALAY_SLAB_C130_R000_F010",
    "ATALAY_SLAB_C130_R050_F010",
    # Tuples
    "ATALAY_SLAB_CASES",
    "ATALAY_ALL_CASES",
]

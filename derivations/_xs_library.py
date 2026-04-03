"""Cross-section library for verification cases.

Four abstract regions (A, B, C, D) in 1, 2, and 4 energy groups.
All XS satisfy the consistency relation Σ_t = Σ_c + Σ_f + Σ_s(row sum).

Region semantics:
    A — fissile (fuel-like)
    B — non-fissile scatterer (moderator-like)
    C — non-fissile absorber (cladding-like)
    D — thin, low-opacity scatterer (gap-like)

Shared by all derivation modules to avoid duplication.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from data.macro_xs.mixture import Mixture


def make_mixture(
    sig_t: np.ndarray,
    sig_c: np.ndarray,
    sig_f: np.ndarray,
    nu: np.ndarray,
    chi: np.ndarray,
    sig_s: np.ndarray,
) -> Mixture:
    """Build a Mixture from N-group arrays."""
    ng = len(sig_t)
    eg = np.logspace(7, -3, ng + 1)
    return Mixture(
        SigC=sig_c.copy(), SigL=np.zeros(ng),
        SigF=sig_f.copy(), SigP=(nu * sig_f).copy(),
        SigT=sig_t.copy(), SigS=[csr_matrix(sig_s)],
        Sig2=csr_matrix((ng, ng)), chi=chi.copy(), eg=eg.copy(),
    )


# ═══════════════════════════════════════════════════════════════════════
# Region A — fissile
# ═══════════════════════════════════════════════════════════════════════

_A_1G = dict(
    sig_t=np.array([1.0]),
    sig_c=np.array([0.2]),
    sig_f=np.array([0.3]),
    nu=np.array([2.5]),
    chi=np.array([1.0]),
    sig_s=np.array([[0.5]]),
)

_A_2G = dict(
    sig_t=np.array([0.50, 1.00]),
    sig_c=np.array([0.01, 0.02]),
    sig_f=np.array([0.01, 0.08]),
    nu=np.array([2.50, 2.50]),
    chi=np.array([1.00, 0.00]),
    sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
)

_A_4G = dict(
    sig_c=np.array([0.01, 0.02, 0.03, 0.05]),
    sig_f=np.array([0.005, 0.01, 0.05, 0.10]),
    nu=np.array([2.80, 2.60, 2.50, 2.45]),
    chi=np.array([0.60, 0.35, 0.05, 0.00]),
    sig_s=np.array([
        [0.28, 0.08, 0.02, 0.005],
        [0.00, 0.40, 0.12, 0.06],
        [0.00, 0.00, 0.55, 0.22],
        [0.00, 0.00, 0.00, 0.90],
    ]),
)
_A_4G["sig_t"] = _A_4G["sig_c"] + _A_4G["sig_f"] + _A_4G["sig_s"].sum(axis=1)


# ═══════════════════════════════════════════════════════════════════════
# Region B — non-fissile scatterer (moderator)
# ═══════════════════════════════════════════════════════════════════════

_B_1G = dict(
    sig_t=np.array([2.0]),
    sig_c=np.array([0.1]),
    sig_f=np.array([0.0]),
    nu=np.array([0.0]),
    chi=np.array([1.0]),
    sig_s=np.array([[1.9]]),
)

_B_2G = dict(
    sig_t=np.array([0.60, 2.00]),
    sig_c=np.array([0.02, 0.05]),
    sig_f=np.array([0.00, 0.00]),
    nu=np.array([0.00, 0.00]),
    chi=np.array([1.00, 0.00]),
    sig_s=np.array([[0.40, 0.18], [0.00, 1.95]]),
)

_B_4G = dict(
    sig_c=np.array([0.02, 0.03, 0.04, 0.06]),
    sig_f=np.array([0.00, 0.00, 0.00, 0.00]),
    nu=np.array([0.00, 0.00, 0.00, 0.00]),
    chi=np.array([0.60, 0.35, 0.05, 0.00]),
    sig_s=np.array([
        [0.50, 0.15, 0.04, 0.01],
        [0.00, 0.70, 0.20, 0.08],
        [0.00, 0.00, 0.85, 0.30],
        [0.00, 0.00, 0.00, 1.80],
    ]),
)
_B_4G["sig_t"] = _B_4G["sig_c"] + _B_4G["sig_f"] + _B_4G["sig_s"].sum(axis=1)


# ═══════════════════════════════════════════════════════════════════════
# Region C — non-fissile absorber (cladding)
# ═══════════════════════════════════════════════════════════════════════

_C_1G = dict(
    sig_t=np.array([1.5]),
    sig_c=np.array([0.15]),
    sig_f=np.array([0.0]),
    nu=np.array([0.0]),
    chi=np.array([1.0]),
    sig_s=np.array([[1.35]]),
)

_C_2G = dict(
    sig_t=np.array([0.55, 1.50]),
    sig_c=np.array([0.02, 0.08]),
    sig_f=np.array([0.00, 0.00]),
    nu=np.array([0.00, 0.00]),
    chi=np.array([1.00, 0.00]),
    sig_s=np.array([[0.38, 0.15], [0.00, 1.42]]),
)

_C_4G = dict(
    sig_c=np.array([0.03, 0.04, 0.06, 0.10]),
    sig_f=np.array([0.00, 0.00, 0.00, 0.00]),
    nu=np.array([0.00, 0.00, 0.00, 0.00]),
    chi=np.array([0.60, 0.35, 0.05, 0.00]),
    sig_s=np.array([
        [0.35, 0.10, 0.03, 0.008],
        [0.00, 0.50, 0.15, 0.06],
        [0.00, 0.00, 0.60, 0.25],
        [0.00, 0.00, 0.00, 1.20],
    ]),
)
_C_4G["sig_t"] = _C_4G["sig_c"] + _C_4G["sig_f"] + _C_4G["sig_s"].sum(axis=1)


# ═══════════════════════════════════════════════════════════════════════
# Region D — thin, low-opacity scatterer (gap)
# ═══════════════════════════════════════════════════════════════════════

_D_1G = dict(
    sig_t=np.array([0.05]),
    sig_c=np.array([0.005]),
    sig_f=np.array([0.0]),
    nu=np.array([0.0]),
    chi=np.array([1.0]),
    sig_s=np.array([[0.045]]),
)

_D_2G = dict(
    sig_t=np.array([0.04, 0.06]),
    sig_c=np.array([0.003, 0.005]),
    sig_f=np.array([0.00, 0.00]),
    nu=np.array([0.00, 0.00]),
    chi=np.array([1.00, 0.00]),
    sig_s=np.array([[0.030, 0.007], [0.000, 0.055]]),
)

_D_4G = dict(
    sig_c=np.array([0.002, 0.003, 0.004, 0.006]),
    sig_f=np.array([0.00, 0.00, 0.00, 0.00]),
    nu=np.array([0.00, 0.00, 0.00, 0.00]),
    chi=np.array([0.60, 0.35, 0.05, 0.00]),
    sig_s=np.array([
        [0.030, 0.008, 0.002, 0.001],
        [0.000, 0.040, 0.010, 0.005],
        [0.000, 0.000, 0.050, 0.020],
        [0.000, 0.000, 0.000, 0.080],
    ]),
)
_D_4G["sig_t"] = _D_4G["sig_c"] + _D_4G["sig_f"] + _D_4G["sig_s"].sum(axis=1)


# ═══════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════

# Organized as XS[region][ng_key] -> dict of arrays
XS: dict[str, dict[str, dict]] = {
    "A": {"1g": _A_1G, "2g": _A_2G, "4g": _A_4G},
    "B": {"1g": _B_1G, "2g": _B_2G, "4g": _B_4G},
    "C": {"1g": _C_1G, "2g": _C_2G, "4g": _C_4G},
    "D": {"1g": _D_1G, "2g": _D_2G, "4g": _D_4G},
}

# Region layouts for each region count
LAYOUTS: dict[int, list[str]] = {
    1: ["A"],
    2: ["A", "B"],
    4: ["A", "D", "C", "B"],
}


def get_xs(region: str, ng_key: str) -> dict:
    """Get XS dict for a region and group count.

    Parameters
    ----------
    region : "A", "B", "C", or "D"
    ng_key : "1g", "2g", or "4g"
    """
    return XS[region][ng_key]


def get_mixture(region: str, ng_key: str) -> Mixture:
    """Get a Mixture for a region and group count."""
    xs = get_xs(region, ng_key)
    return make_mixture(**xs)


def get_materials(n_regions: int, ng_key: str) -> dict[int, Mixture]:
    """Get materials dict for a given number of regions and group count.

    Material IDs start from 0 (outermost) counting inward, matching the
    convention where higher IDs are more interior (fuel = highest ID).
    """
    layout = LAYOUTS[n_regions]
    # Assign IDs: A gets the highest (innermost), last in layout gets 0
    mats = {}
    for i, region in enumerate(layout):
        mat_id = len(layout) - 1 - i
        mats[mat_id] = get_mixture(region, ng_key)
    return mats


def validate_all() -> None:
    """Verify XS consistency for all regions and group counts."""
    for region, groups in XS.items():
        for ng_key, xs in groups.items():
            sig_t_check = xs["sig_c"] + xs["sig_f"] + xs["sig_s"].sum(axis=1)
            assert np.allclose(xs["sig_t"], sig_t_check), (
                f"XS inconsistency in region {region}, {ng_key}: "
                f"sig_t={xs['sig_t']} ≠ sig_c+sig_f+sig_s={sig_t_check}"
            )


# Validate on import
validate_all()

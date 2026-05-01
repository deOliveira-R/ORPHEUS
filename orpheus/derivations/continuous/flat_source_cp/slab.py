"""Semi-analytical slab collision probability eigenvalues.

Thin facade over :mod:`~orpheus.derivations.continuous.flat_source_cp.geometry` with
:data:`~.cp_geometry.SLAB` pre-selected. Derives ``k_inf`` for
{1, 2, 4} energy groups × {1, 2, 4} regions using the :math:`E_3`
exponential-integral kernel. See :doc:`/theory/peierls_unified`
§§11-17 for the three-tier integration hierarchy and the
derivation of the geometry-invariant :math:`\\Delta^{2}` operator.
"""

from __future__ import annotations

import numpy as np

from . import geometry as _cpg
from ...common.eigenvalue import kinf_from_cp
from ...common.verification_case import VerificationCase
from ...common.xs_library import LAYOUTS, get_xs, get_mixture


# ═══════════════════════════════════════════════════════════════════════
# Geometry singleton (binds the unified infrastructure)
# ═══════════════════════════════════════════════════════════════════════

GEOMETRY = _cpg.SLAB


# ═══════════════════════════════════════════════════════════════════════
# Backward-compatible slab CP matrix
# ═══════════════════════════════════════════════════════════════════════

def _slab_cp_matrix(
    sig_t_all: np.ndarray,
    t_arr: np.ndarray,
) -> np.ndarray:
    """Compute the infinite-lattice CP matrix for a slab.

    Delegates to :func:`cp_geometry.build_cp_matrix` with the
    pre-bound :data:`~.cp_geometry.SLAB` geometry. The slab's
    "volumes" are the thicknesses (1-D, V == t), and ``R_cell``
    is the total slab thickness for the white-BC surface-area
    normalisation (which returns 1 for the slab).

    Parameters
    ----------
    sig_t_all : (N_reg, ng) — total XS per region and group
    t_arr : (N_reg,) — region thicknesses

    Returns
    -------
    P_inf : (N_reg, N_reg, ng) — collision probability matrix
    """
    return _cpg.build_cp_matrix(
        GEOMETRY,
        sig_t_all=sig_t_all,
        radii_or_thicknesses=np.asarray(t_arr, dtype=float),
        volumes=np.asarray(t_arr, dtype=float),
        R_cell=float(np.sum(t_arr)),
    )


# ═══════════════════════════════════════════════════════════════════════
# Slab geometry parameters
# ═══════════════════════════════════════════════════════════════════════

# Region thicknesses (innermost to outermost)
_THICKNESSES = {
    1: [0.5],                       # A only
    2: [0.5, 0.5],                  # A + B
    4: [0.4, 0.05, 0.1, 0.45],     # A + D + C + B
}

# Material IDs per region count, matching solver convention
# (innermost = highest ID)
_MAT_IDS = {
    1: [2],          # A → fuel(2)
    2: [2, 0],       # A → fuel(2), B → cool(0)
    4: [2, 3, 1, 0], # A → fuel(2), D → gap(3), C → clad(1), B → cool(0)
}


# ═══════════════════════════════════════════════════════════════════════
# Case generation
# ═══════════════════════════════════════════════════════════════════════

def _build_case(ng_key: str, n_regions: int) -> VerificationCase:
    """Build a slab CP verification case for given groups and regions."""
    layout = LAYOUTS[n_regions]
    ng = int(ng_key[0])
    t_arr = np.array(_THICKNESSES[n_regions])

    # Collect XS per region (innermost first)
    xs_list = [get_xs(region, ng_key) for region in layout]
    sig_t_all = np.vstack([xs["sig_t"] for xs in xs_list])

    P_inf_g = _slab_cp_matrix(sig_t_all, t_arr)

    k_inf = kinf_from_cp(
        P_inf_g=P_inf_g,
        sig_t_all=sig_t_all,
        V_arr=t_arr,
        sig_s_mats=[xs["sig_s"] for xs in xs_list],
        nu_sig_f_mats=[xs["nu"] * xs["sig_f"] for xs in xs_list],
        chi_mats=[xs["chi"] for xs in xs_list],
    )

    # Build materials dict with mat_ids matching the solver convention
    mat_ids = _MAT_IDS[n_regions]
    materials = {}
    for i, region in enumerate(layout):
        materials[mat_ids[i]] = get_mixture(region, ng_key)

    # Geometry params for building SlabGeometry in solver tests
    geom_params_out = dict(
        thicknesses=t_arr.tolist(),
        mat_ids=mat_ids,
    )

    name = f"cp_slab_{ng}eg_{n_regions}rg"
    dim = n_regions * ng

    latex = (
        rf"Slab CP eigenvalue with {ng} groups, {n_regions} regions, "
        r"white boundary condition. "
        rf"The E₃-based CP matrix yields a {dim}×{dim} eigenvalue problem."
        "\n\n"
        r".. math::" "\n"
        rf"   k_\infty = {k_inf:.10f}"
    )

    labels: list[str] = ["collision-rate", "e3-def", "self-slab", "p-inf"]
    if n_regions > 1:
        labels += ["dd-slab", "dc-slab", "second-diff-general"]
    if ng == 1 and n_regions == 1:
        labels.append("one-group-kinf")
    if ng > 1:
        labels += ["matrix-eigenvalue", "mg-balance"]

    return VerificationCase(
        name=name,
        k_inf=k_inf,
        method="cp",
        geometry="slab",
        n_groups=ng,
        n_regions=n_regions,
        materials=materials,
        geom_params=geom_params_out,
        latex=latex,
        description=f"{ng}G {n_regions}-region slab CP (E₃ kernel, white BC)",
        tolerance="< 1e-6",
        vv_level="L1",
        equation_labels=tuple(labels),
    )


def all_cases() -> list[VerificationCase]:
    """Return all slab CP verification cases: {1,2,4}eg × {1,2,4}rg."""
    cases = []
    for ng_key in ["1g", "2g", "4g"]:
        for n_regions in [1, 2, 4]:
            cases.append(_build_case(ng_key, n_regions))
    return cases

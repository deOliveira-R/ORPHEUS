r"""The Case two-region slab of ``sn_slab_1eg_2rg_S8``, posed once for the gates that solve it.

One group, a fuel and a reflector region, reflective on both faces, S8
Gauss-Legendre: the problem whose singular-eigenfunction (Case) eigenvalue is
the continuous reference ``sn_slab_1eg_2rg_S8``. The L1 standoff rows
(``test_l1_standoff_slab_cylinder.py``) compare SN's k against it, and the
sweep-preconditioned Krylov gates (``tests/gates/sn/solve/
test_krylov_sweep_preconditioner.py``, issue #200) solve the same slab for
their fixed-point and rate rows.
"""
from __future__ import annotations

from orpheus.derivations.reference_values import continuous_get
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher


def case_slab_mesh(n_per: int) -> tuple[Mesh1D, dict, int]:
    r"""Build the ``sn_slab_1eg_2rg_S8`` Case singular-eigenfunction mesh.

    Returns ``(mesh, materials, N_ord)`` for direct use by ``solve_sn``.
    """
    ref = continuous_get("sn_slab_1eg_2rg_S8")
    geom = ref.problem.geometry_params
    materials = ref.problem.materials
    H_A = float(geom["fuel_height"])
    H_B = float(geom["refl_height"])
    N_ord = int(geom["n_ordinates"])
    slab = StructuredGeometry.slab(
        (0.0, H_A, H_A + H_B), (0, 1), left=BC.reflective, right=BC.reflective,
    )
    mesh = Mesher(slab).partition(CellsByCount.uniform_width(n_per)).mesh
    return mesh, materials, N_ord


def case_slab_k_ref() -> float:
    k_eff = continuous_get("sn_slab_1eg_2rg_S8").k_eff
    if k_eff is None:
        raise ValueError("sn_slab_1eg_2rg_S8 carries no eigenvalue: the slab reference answers the k question")
    return float(k_eff)

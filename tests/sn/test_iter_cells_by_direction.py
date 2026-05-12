r"""Foundation tests for ``SNMesh.iter_cells_by_direction`` + ``EquationMap.unknowns_at_cell_for_mask``.

Issue #168 Phase C — API helpers underpinning the sweep-frame apply
matvec rewrite. The matvec operates on whole ordinate subsets
simultaneously using ``outgoing_mask`` / ``incoming_mask`` and consumes
the direction-only cell ordering exposed by
:meth:`SNMesh.iter_cells_by_direction`.

These foundation tests pin the equivalence to existing iteration
patterns at bit precision so the sweep-frame matvec has a stable
substrate to build on.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    EquationMap,
    build_equation_map_cylindrical,
    build_equation_map_spherical,
)
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LevelSymmetricSN,
    ProductQuadrature,
)


# ═══════════════════════════════════════════════════════════════════════
# SNMesh.iter_cells_by_direction — equivalence to iter_cell_visits
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_iter_cells_by_direction_spherical_outward_matches_iter_cell_visits():
    """Spherical +1 direction yields same cell sequence as any μ ≥ 0 ordinate."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = GaussLegendre1D.create(8)
    sn_mesh = SNMesh(mesh, quad)

    seq_by_dir = [
        v.cell_idx for v in sn_mesh.iter_cells_by_direction(+1)
    ]
    # Every μ ≥ 0 ordinate must produce the same cell sequence.
    for n in range(quad.N):
        if quad.mu_x[n] >= 0:
            seq_by_ord = [
                v.cell_idx for v in sn_mesh.iter_cell_visits(ordinate_idx=n)
            ]
            assert np.array_equal(seq_by_dir, seq_by_ord), (
                f"Spherical +1: direction-keyed sequence != "
                f"ordinate {n} (μ={quad.mu_x[n]:.4f}) sequence"
            )


@pytest.mark.foundation
def test_iter_cells_by_direction_spherical_inward_matches_iter_cell_visits():
    """Spherical -1 direction yields same cell sequence as any μ < 0 ordinate."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = GaussLegendre1D.create(8)
    sn_mesh = SNMesh(mesh, quad)

    seq_by_dir = [
        v.cell_idx for v in sn_mesh.iter_cells_by_direction(-1)
    ]
    for n in range(quad.N):
        if quad.mu_x[n] < 0:
            seq_by_ord = [
                v.cell_idx for v in sn_mesh.iter_cell_visits(ordinate_idx=n)
            ]
            assert np.array_equal(seq_by_dir, seq_by_ord), (
                f"Spherical -1: direction-keyed sequence != "
                f"ordinate {n} (μ={quad.mu_x[n]:.4f}) sequence"
            )


@pytest.mark.foundation
def test_iter_cells_by_direction_slab_matches_iter_cell_visits():
    """Slab (1-D Cartesian) sweep direction yields the same cell sequence."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 13),
        mat_ids=np.zeros(12, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = GaussLegendre1D.create(6)
    sn_mesh = SNMesh(mesh, quad)

    for sign in (+1, -1):
        seq_by_dir = [
            v.cell_idx for v in sn_mesh.iter_cells_by_direction(sign)
        ]
        for n in range(quad.N):
            if (sign == +1 and quad.mu_x[n] >= 0) or (
                sign == -1 and quad.mu_x[n] < 0
            ):
                seq_by_ord = [
                    v.cell_idx
                    for v in sn_mesh.iter_cell_visits(ordinate_idx=n)
                ]
                assert np.array_equal(seq_by_dir, seq_by_ord)


@pytest.mark.foundation
def test_iter_cells_by_direction_cylindrical_per_level_matches():
    """Cylindrical per-level direction yields same cell sequence as level ordinates."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.5, 9),
        mat_ids=np.zeros(8, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = ProductQuadrature.create(n_mu=2, n_phi=4)
    sn_mesh = SNMesh(mesh, quad)

    level_indices = quad.level_indices
    for level_p in range(len(level_indices)):
        level_ords = np.asarray(level_indices[level_p])
        eta_at_level = quad.mu_x[level_ords]
        for sign in (+1, -1):
            # Exclude pure-azimuthal degenerate ordinates: their
            # iter_cell_visits path always iterates forward (no
            # signed direction), which would NOT match the
            # direction-keyed sweep. The degenerate branch is the
            # special case |η| < 1e-15.
            if sign == +1:
                matches = np.where(eta_at_level > 1e-15)[0]
            else:
                matches = np.where(eta_at_level < -1e-15)[0]
            if matches.size == 0:
                continue
            seq_by_dir = [
                v.cell_idx
                for v in sn_mesh.iter_cells_by_direction(sign, level_p)
            ]
            for within in matches:
                seq_by_ord = [
                    v.cell_idx
                    for v in sn_mesh.iter_cell_visits(
                        ordinate_idx=int(within),
                        mu_level_idx=level_p,
                    )
                ]
                assert np.array_equal(seq_by_dir, seq_by_ord), (
                    f"Cyl level {level_p} dir {sign} != "
                    f"within-level idx {within}"
                )


@pytest.mark.foundation
def test_iter_cells_by_direction_invalid_sign_raises():
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = GaussLegendre1D.create(4)
    sn_mesh = SNMesh(mesh, quad)
    with pytest.raises(ValueError, match="direction_sign"):
        list(sn_mesh.iter_cells_by_direction(0))
    with pytest.raises(ValueError, match="direction_sign"):
        list(sn_mesh.iter_cells_by_direction(2))


@pytest.mark.foundation
def test_iter_cells_by_direction_cylindrical_requires_level():
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = ProductQuadrature.create(n_mu=2, n_phi=4)
    sn_mesh = SNMesh(mesh, quad)
    with pytest.raises(ValueError, match="mu_level_idx"):
        list(sn_mesh.iter_cells_by_direction(+1))


# ═══════════════════════════════════════════════════════════════════════
# EquationMap.unknowns_at_cell_for_mask
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_unknowns_at_cell_for_mask_spherical_matches_brute_force():
    """For each cell + mask combo, the helper returns the same k set as a linear scan."""
    quad = GaussLegendre1D.create(8)
    nx = 10
    ng = 2
    eq_map = build_equation_map_spherical(nx, quad, ng)

    outgoing = quad.mu_x > 1e-15
    incoming = quad.mu_x < -1e-15
    for i in range(nx):
        for mask in (outgoing, incoming):
            # Brute-force reference: scan eq_map.ix/ordinate.
            ref_ks = [
                k
                for k in range(eq_map.n_eq)
                if eq_map.ix[k] == i and mask[eq_map.ordinate[k]]
            ]
            ref_ks = np.asarray(ref_ks, dtype=np.intp)
            helper_ks = eq_map.unknowns_at_cell_for_mask(i, mask)
            # The helper orders by ascending ordinate index; the
            # eq_map was built ordinate-inner, so ref_ks naturally
            # matches that order too.
            assert np.array_equal(np.sort(helper_ks), np.sort(ref_ks)), (
                f"Cell {i}: helper {helper_ks} != ref {ref_ks}"
            )


@pytest.mark.foundation
def test_unknowns_at_cell_for_mask_cylindrical_matches_brute_force():
    quad = LevelSymmetricSN.create(4)
    nx = 6
    ng = 1
    eq_map = build_equation_map_cylindrical(nx, quad, ng)

    outgoing = quad.mu_x > 1e-15
    for i in range(nx):
        ref_ks = sorted(
            k
            for k in range(eq_map.n_eq)
            if eq_map.ix[k] == i and outgoing[eq_map.ordinate[k]]
        )
        helper_ks = sorted(
            eq_map.unknowns_at_cell_for_mask(i, outgoing).tolist()
        )
        assert helper_ks == ref_ks


@pytest.mark.foundation
def test_unknowns_at_cell_for_mask_empty_mask_yields_empty():
    quad = GaussLegendre1D.create(4)
    eq_map = build_equation_map_spherical(5, quad, 1)
    empty_mask = np.zeros(quad.N, dtype=bool)
    for i in range(5):
        result = eq_map.unknowns_at_cell_for_mask(i, empty_mask)
        assert result.size == 0


@pytest.mark.foundation
def test_unknowns_at_cell_for_mask_lazy_table_caches():
    quad = GaussLegendre1D.create(4)
    eq_map = build_equation_map_spherical(5, quad, 1)
    assert eq_map._cell_ord_to_k is None  # type: ignore[attr-defined]
    mask = quad.mu_x > 0
    eq_map.unknowns_at_cell_for_mask(0, mask)
    assert eq_map._cell_ord_to_k is not None  # type: ignore[attr-defined]
    table_id = id(eq_map._cell_ord_to_k)  # type: ignore[attr-defined]
    eq_map.unknowns_at_cell_for_mask(1, mask)
    # Cached; second call reuses the table object.
    assert id(eq_map._cell_ord_to_k) == table_id  # type: ignore[attr-defined]


@pytest.mark.foundation
def test_unknowns_at_cell_for_mask_outer_boundary_inward_absent():
    r"""At i=nx-1, inward μ<0 ordinates are NOT unknowns (BC determines them).

    The eq_map for spherical skips inward-at-outer-boundary slots
    (line 399 in operator.py); the helper must return an empty array
    when queried with an incoming mask at the outer boundary.
    """
    quad = GaussLegendre1D.create(6)
    nx = 8
    eq_map = build_equation_map_spherical(nx, quad, 1)
    incoming = quad.mu_x < -1e-15
    result = eq_map.unknowns_at_cell_for_mask(nx - 1, incoming)
    assert result.size == 0

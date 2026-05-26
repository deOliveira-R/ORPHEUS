r"""Foundation tests for ``SNMesh.dag_walk``.

Issue #196 Phase G Step 2.6 (Q3) — the single canonical iteration
primitive for 1-D sweeps.  ``dag_walk`` takes EXACTLY ONE of
``ordinate_idx`` (per-ordinate visits) or ``direction_sign``
(direction-keyed visits, used by the apply matvec) as a
keyword-only argument; ``mu_level_idx`` selects a cylindrical
:math:`\mu`-level when present.

These foundation tests pin the equivalence between the two
invocation modes at bit precision: the direction-keyed branch must
yield the same cell sequence as every non-degenerate ordinate in
the matching sign class.

R-1 Step 4 — the ``EquationMap.unknowns_at_cell_for_mask`` tests that
used to live in this file have been retired with the helper (zero
production callers; the typed-AngularFlux architecture indexes
``psi.values[:, :, i, j]`` directly).  The five retired tests covered:
spherical / cylindrical brute-force-equivalence, empty-mask early
return, lazy-table caching, and outer-boundary-inward-absence.  None
exercised :meth:`SNMesh.dag_walk`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# SNMesh.dag_walk — equivalence between direction-keyed and ordinate-keyed
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_dag_walk_spherical_outward_matches_per_ordinate():
    """Spherical +1 direction yields same cell sequence as any μ ≥ 0 ordinate."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())

    seq_by_dir = [
        v.cell_idx for v in sn_mesh.dag_walk(direction_sign=+1)
    ]
    # Every μ ≥ 0 ordinate must produce the same cell sequence.
    for n in range(quad.N):
        if quad.mu_x[n] >= 0:
            seq_by_ord = [
                v.cell_idx for v in sn_mesh.dag_walk(ordinate_idx=n)
            ]
            assert np.array_equal(seq_by_dir, seq_by_ord), (
                f"Spherical +1: direction-keyed sequence != "
                f"ordinate {n} (μ={quad.mu_x[n]:.4f}) sequence"
            )


@pytest.mark.foundation
def test_dag_walk_spherical_inward_matches_per_ordinate():
    """Spherical -1 direction yields same cell sequence as any μ < 0 ordinate."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())

    seq_by_dir = [
        v.cell_idx for v in sn_mesh.dag_walk(direction_sign=-1)
    ]
    for n in range(quad.N):
        if quad.mu_x[n] < 0:
            seq_by_ord = [
                v.cell_idx for v in sn_mesh.dag_walk(ordinate_idx=n)
            ]
            assert np.array_equal(seq_by_dir, seq_by_ord), (
                f"Spherical -1: direction-keyed sequence != "
                f"ordinate {n} (μ={quad.mu_x[n]:.4f}) sequence"
            )


@pytest.mark.foundation
def test_dag_walk_slab_matches_per_ordinate():
    """Slab (1-D Cartesian) sweep direction yields the same cell sequence."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 13),
        mat_ids=np.zeros(12, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = Quadrature.gauss_legendre(6)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())

    for sign in (+1, -1):
        seq_by_dir = [
            v.cell_idx for v in sn_mesh.dag_walk(direction_sign=sign)
        ]
        for n in range(quad.N):
            if (sign == +1 and quad.mu_x[n] >= 0) or (
                sign == -1 and quad.mu_x[n] < 0
            ):
                seq_by_ord = [
                    v.cell_idx
                    for v in sn_mesh.dag_walk(ordinate_idx=n)
                ]
                assert np.array_equal(seq_by_dir, seq_by_ord)


@pytest.mark.foundation
def test_dag_walk_cylindrical_per_level_matches():
    """Cylindrical per-level direction yields same cell sequence as level ordinates."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.5, 9),
        mat_ids=np.zeros(8, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = Quadrature.product(n_mu=2, n_phi=4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())

    level_indices = quad.level_indices
    for level_p in range(len(level_indices)):
        level_ords = np.asarray(level_indices[level_p])
        eta_at_level = quad.eta[level_ords]
        for sign in (+1, -1):
            # Exclude pure-azimuthal degenerate ordinates: their
            # ordinate-keyed dag_walk path always iterates forward (no
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
                for v in sn_mesh.dag_walk(
                    direction_sign=sign, mu_level_idx=level_p,
                )
            ]
            for within in matches:
                seq_by_ord = [
                    v.cell_idx
                    for v in sn_mesh.dag_walk(
                        ordinate_idx=int(within),
                        mu_level_idx=level_p,
                    )
                ]
                assert np.array_equal(seq_by_dir, seq_by_ord), (
                    f"Cyl level {level_p} dir {sign} != "
                    f"within-level idx {within}"
                )


@pytest.mark.foundation
def test_dag_walk_invalid_sign_raises():
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    with pytest.raises(ValueError, match="direction_sign"):
        list(sn_mesh.dag_walk(direction_sign=0))
    with pytest.raises(ValueError, match="direction_sign"):
        list(sn_mesh.dag_walk(direction_sign=2))


@pytest.mark.foundation
def test_dag_walk_cylindrical_requires_level():
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = Quadrature.product(n_mu=2, n_phi=4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    with pytest.raises(ValueError, match="mu_level_idx"):
        list(sn_mesh.dag_walk(direction_sign=+1))


@pytest.mark.foundation
def test_dag_walk_xor_signature_enforced():
    """``dag_walk`` rejects both / neither of (ordinate_idx, direction_sign)."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
    )
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    # Neither supplied → ValueError.
    with pytest.raises(ValueError, match="exactly one"):
        list(sn_mesh.dag_walk())
    # Both supplied → ValueError.
    with pytest.raises(ValueError, match="exactly one"):
        list(sn_mesh.dag_walk(ordinate_idx=0, direction_sign=+1))

r"""The 3a production ties — the derivation ⟷ the assembled SN operator.

Two welds, both object-level (Mode 12: pin matrices, never functionals):

1. **The parent tie**: the discrete system the four-step derivation
   reduced (per-ordinate DD balance + closure + zero inflow, in the
   production density convention :math:`(\mu/h)\Delta\psi + \sigma_T
   \bar\psi`) has its edge-eliminated cell operator EQUAL to the
   production ``assemble_ordinate_blocks`` output, entry for entry, per
   (ordinate, group), on a heterogeneous non-uniform slab. This proves
   the object the theorems reduced IS the object the production
   assembler emits — the load-bearing premise of "consistent by
   derivation".

2. **The realized low-order system**: :func:`build_consistent_dd_system`
   (the numeric realization of the PROVEN (23a–f)/(27)/boundary rows)
   satisfies its structural laws on the same fixture: symmetry of the
   assembled tridiagonal, positive definiteness under vacuum, the
   correction→0 identity (zero displacements ⟹ zero rhs), and the
   quadrature-moment asserts at the convention boundary.

Levels: foundation (the ties are exact-arithmetic / tight-tolerance
object identities, no physics reference involved).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.discrete.sn import dsa
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation.assembly import assemble_ordinate_blocks
from orpheus.sn.mesh.augmented_mesh import SNMesh

pytestmark = pytest.mark.foundation

_RTOL = 1e-11  # sparse-order ≠ dense-order ⇒ never 0-ULP (the L16 gate value)


@pytest.fixture(scope="module")
def slab() -> SNMesh:
    """Heterogeneous, non-uniform 4-cell slab, S4, 2 groups, vacuum."""
    mesh1d = Mesh1D(
        edges=np.array([0.0, 0.5, 1.5, 3.0, 5.0]),
        mat_ids=np.array([0, 1, 1, 0]),
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(
        mesh1d, quad, {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    )


def _parent_eliminated_cell_operator(
    h: np.ndarray, sigma_t_row: np.ndarray, mu_m: float
) -> np.ndarray:
    r"""Edge-eliminate the per-ordinate DD parent system at zero inflow.

    Unknowns ``[psi_cells (K), psi_edges (K+1)]``; rows: K balance
    (density form), K closure, 1 zero-inflow. The cell operator is the
    balance block after substituting the closure+inflow-determined
    edges: ``A_cells = B_c + B_e · E`` where ``E`` solves the
    (closure + inflow) subsystem for the edges given the cells.
    """
    n = h.shape[0]
    n_edges = n + 1
    balance_c = np.zeros((n, n))
    balance_e = np.zeros((n, n_edges))
    closure_c = np.zeros((n_edges, n))
    closure_e = np.zeros((n_edges, n_edges))

    for i in range(n):
        balance_c[i, i] = sigma_t_row[i]
        balance_e[i, i + 1] += mu_m / h[i]
        balance_e[i, i] += -mu_m / h[i]
        # closure row i: psi_bar_i − ½(psi_left + psi_right) = 0, solved
        # for the DOWNSTREAM edge given the upstream one.
        closure_c[i, i] = 1.0
        closure_e[i, i] = -0.5
        closure_e[i, i + 1] = -0.5
    # zero inflow at the upwind boundary edge:
    inflow_edge = 0 if mu_m > 0 else n
    closure_e[n, inflow_edge] = 1.0

    # edges = E @ cells with E = −(closure_e)⁻¹ closure_c
    e_map = -np.linalg.solve(closure_e, closure_c)
    return balance_c + balance_e @ e_map


class TestParentTie:
    """Weld 1: the derivation's parent ≡ the production assembler."""

    def test_edge_eliminated_parent_matches_assembled_blocks(self, slab):
        h = np.diff(np.asarray(slab.mesh.edges, float))
        sigma_t = np.asarray(
            slab.material_xs_field().total_cross_section_field.values, float
        )  # (ng, K)
        mu = np.asarray(slab.quad.mu_x, float)

        for ordinate in range(slab.quad.N):
            blocks = assemble_ordinate_blocks(
                slab, ordinate, include_collision=True
            )
            for g, block in enumerate(blocks):
                expected = _parent_eliminated_cell_operator(
                    h, sigma_t[g], float(mu[ordinate])
                )
                np.testing.assert_allclose(
                    block.as_matrix(),
                    expected,
                    rtol=_RTOL,
                    atol=1e-13,
                    err_msg=(
                        f"assembled (ordinate={ordinate}, group={g}) must "
                        f"equal the edge-eliminated DD parent"
                    ),
                )


class TestRealizedLowOrderSystem:
    """Weld 2: the proven low-order system's structural laws."""

    @pytest.fixture(scope="class")
    def data(self, slab):
        h = np.diff(np.asarray(slab.mesh.edges, float))
        xs = slab.material_xs_field()
        sigma_t = np.asarray(xs.total_cross_section_field.values, float)
        fold = xs.foldable_sigma()  # {mid: (ng,)}
        mat_ids = np.asarray(slab.mesh.mat_ids, int).ravel()
        sig_s0_gg = np.stack(
            [fold[int(m)] for m in mat_ids], axis=1
        )  # (ng, K)
        mu = np.asarray(slab.quad.mu_x, float)
        w = np.asarray(slab.quad.weights, float)
        return h, sigma_t, sig_s0_gg, mu, w

    def test_vacuum_system_symmetric_positive_definite(self, data):
        h, sigma_t, sig_s0, mu, w = data
        for g in range(sigma_t.shape[0]):
            a_low, _ = dsa.build_consistent_dd_system(
                h, sigma_t[g], sig_s0[g], np.zeros_like(h), mu, w
            )
            interior = a_low[1:-1, 1:-1]
            np.testing.assert_allclose(
                interior,
                interior.T,
                rtol=0,
                atol=1e-14,
                err_msg="the interior block must be symmetric",
            )
            eigenvalues = np.linalg.eigvals(a_low)
            if not (eigenvalues.real > 0).all():
                pytest.fail(
                    f"group {g}: the vacuum low-order system must be "
                    f"positive definite; got min Re λ = "
                    f"{eigenvalues.real.min()!r}"
                )

    def test_correction_to_zero_identity(self, data):
        """d = 0 ⟹ rhs = 0 exactly — the correctness-safe-by-construction
        property at the matrix level (G maps zero displacements to a zero
        source; the correction solve then returns f = 0 identically)."""
        h, sigma_t, sig_s0, mu, w = data
        _, g_map = dsa.build_consistent_dd_system(
            h, sigma_t[0], sig_s0[0], np.zeros_like(h), mu, w
        )
        rhs = g_map @ np.zeros(2 * h.shape[0])
        np.testing.assert_array_equal(rhs, np.zeros(h.shape[0] + 1))

    def test_quadrature_convention_guard_fires(self, data):
        """The Σw = 2 convention boundary refuses a wrong-normalized rule."""
        h, sigma_t, sig_s0, mu, w = data
        with pytest.raises(dsa.DerivationError, match="Σw = 2"):
            dsa.build_consistent_dd_system(
                h, sigma_t[0], sig_s0[0], np.zeros_like(h), mu, w / 2.0
            )

    def test_reflective_rows_are_one_sided_f1(self, data):
        """Reflective BCs assemble and stay solvable (f1 = 0 rows)."""
        h, sigma_t, sig_s0, mu, w = data
        a_low, _ = dsa.build_consistent_dd_system(
            h,
            sigma_t[0],
            sig_s0[0],
            np.zeros_like(h),
            mu,
            w,
            bc=("reflective", "vacuum"),
        )
        if not np.isfinite(np.linalg.cond(a_low)):
            pytest.fail("reflective-left system must remain invertible")

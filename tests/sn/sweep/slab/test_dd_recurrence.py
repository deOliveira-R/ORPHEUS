"""L0 term-verification of the slab diamond-difference cumprod recurrence.

Split from the legacy ``tests/sn/test_cartesian.py`` (SN taxonomy reorg):
this is the per-cell DD recurrence gate (T2a sweep / slab closure), the
most direct possible check that ``DiamondDifference.update`` matches the
symbolic derivation in ``sn.balance.derive_cumprod_recurrence``. The
eigenvalue / convergence claims that shared the legacy file moved to
``tests/sn/eigenvalue/test_keff_slab.py``.
"""

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature

# The full equation-coverage list from the legacy test_cartesian module
# is preserved verbatim on every split so no verifies(...) edge is lost.
pytestmark = pytest.mark.verifies(
    "transport-cartesian",
    "dd-cartesian-1d",
    "dd-solve",
    "dd-recurrence",
    "multigroup",
    "reflective-bc",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
)


# ─── L0 term verification of the DD cumprod recurrence (ERR-025) ────

@pytest.mark.sentinel
@pytest.mark.l0
@pytest.mark.catches("ERR-025")
def test_dd_per_cell_recurrence_matches_symbolic_derivation():
    """Term-level verification that ``DiamondDifference.update``'s
    per-cell recurrence matches the symbolic derivation in
    :func:`orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence`.

    Issue #196 Step 2.5: ``_sweep_1d_cumprod`` was retired in favour
    of a unified fold over DAG visits that delegates per-cell
    algebra to :class:`DiamondDifference`.  This test directly
    invokes ``DiamondDifference.update`` on a synthetic per-cell
    input and checks the ``ψ_cell = ½(ψ_in + a·ψ_in + b·Q/W)``
    identity that gated the legacy ERR-025 coefficient drift.  The
    DD strategy's per-cell algebra IS the recurrence, so this is the
    most direct gate possible — no sweep / BC / boundary scaffolding
    in the way.
    """
    import sympy as sp
    from orpheus.derivations.discrete.sn.balance import derive_cumprod_recurrence
    from orpheus.geometry import CoordSystem, Mesh1D
    from orpheus.geometry.reduced_operator import slab_streaming
    from orpheus.sn.spatial.cell_update import CellVisit, UpstreamState
    from orpheus.sn.spatial.diamond import DiamondDifference

    # Symbolic coefficients, captured silently.
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        a_sym, b_sym = derive_cumprod_recurrence()

    mu_sym, dx_sym, Sig_t_sym, S_sym = sp.symbols(
        "mu dx Sigma_t S", positive=True
    )

    ng = 1
    sig_t_val = 1.5
    dx_val = 0.7
    Q_val = 3.0

    quad = Quadrature.gauss_legendre(4)
    edges = np.array([0.0, dx_val])
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(1, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    op = slab_streaming(mesh, quad)
    W = quad.weights.sum()
    n_half = quad.N // 2
    mu_pos = np.abs(quad.mu_x[n_half:])

    # Synthetic ψ_in per positive ordinate.
    psi_in_per_ordinate = [0.4, 0.9]
    strat = DiamondDifference()
    total_xs = np.array([sig_t_val])

    for n in range(n_half):
        mu_val = float(mu_pos[n])
        direction_idx = n_half + n
        st = op.streaming_terms(cell_idx=0, direction_idx=direction_idx)

        # The contract source is Q · V · weight_norm = Q · dx / W.
        source = np.array([Q_val * dx_val / W])
        psi_in = np.array([psi_in_per_ordinate[n]])
        upstream = UpstreamState(
            spatial_upstream=psi_in, angular_upstream=None,
        )
        visit = CellVisit(
            cell_idx=0, streaming_terms=st, face_area_downstream=1.0,
        )
        result = strat.update(visit, total_xs, source, upstream)
        cell_avg_code = float(result.cell_average_flux[0])

        # Symbolic-derived reference.
        a_num = float(a_sym.subs({mu_sym: mu_val, dx_sym: dx_val,
                                  Sig_t_sym: sig_t_val}))
        b_num = float(b_sym.subs({mu_sym: mu_val, dx_sym: dx_val,
                                  Sig_t_sym: sig_t_val,
                                  S_sym: Q_val / W}))
        psi_out_expected = a_num * psi_in[0] + b_num
        cell_avg_expected = 0.5 * (psi_in[0] + psi_out_expected)

        # Issue #196 Step 2.5: DD's unified body computes
        # ``(source + numer_upstream)/denom`` (algebraically
        # identical to the cumprod's ``a·ψ_in + 2q/denom; ½(ψ_in +
        # ψ_out)`` but ULP-different).  Re-baseline rtol=1e-13.
        assert abs(cell_avg_code - cell_avg_expected) < 1e-12, (
            f"Ordinate n={n} (μ={mu_val:.6f}): "
            f"DD gave {cell_avg_code:.10e}, "
            f"derivation gives {cell_avg_expected:.10e}, "
            f"Δ={cell_avg_code - cell_avg_expected:+.2e}. "
            "DiamondDifference does not match "
            "sn_balance.derive_cumprod_recurrence."
        )

        # Spatial WDD closure: ψ_out = 2·ψ_avg − ψ_in MUST equal the
        # symbolic ψ_out_expected.  Asserting this (NOT just cell_avg)
        # gives the sentinel sensitivity to the ``psi_spat_out`` sign /
        # factor in DiamondDifference.update — a sentinel-mutation gap
        # found in Phase S2 (the cell_avg-only assertion left the
        # spatial-closure mutants alive).  np.testing.assert_allclose
        # is a function call, so it fires even under ``-O`` (vv Mode 8).
        np.testing.assert_allclose(
            float(result.outgoing_spatial_flux[0]), psi_out_expected,
            rtol=1e-12,
            err_msg=(
                f"Ordinate n={n} (μ={mu_val:.6f}): WDD spatial closure "
                f"ψ_out=2ψ_avg−ψ_in disagrees with symbolic ψ_out."
            ),
        )

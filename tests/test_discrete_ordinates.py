"""Verify the 2D discrete ordinates solver (mesh convergence)."""

import numpy as np
import pytest

from derivations import get
from derivations.cp_slab import _XS_A, _XS_B, _make_mixture
from discrete_ordinates import DOParams, PinCellGeometry, solve_discrete_ordinates


def _slab_geom_for_do(n_fuel, n_mod, delta, ny=2):
    """Build an SN PinCellGeometry for a 2-region slab."""
    nx = n_fuel + n_mod
    vol = np.full((nx, ny), delta**2)
    vol[0, :] /= 2
    vol[-1, :] /= 2
    vol[:, 0] /= 2
    vol[:, -1] /= 2

    mat = np.zeros((nx, ny), dtype=int)
    mat[:n_fuel, :] = 2
    mat[n_fuel:, :] = 0

    return PinCellGeometry(nx=nx, ny=ny, delta=delta, mat_map=mat, volume=vol)


@pytest.mark.slow
@pytest.mark.parametrize("suffix,label", [("1g", "1G"), ("2g", "2G")])
def test_do_mesh_convergence(suffix, label):
    """2D SN solver must converge with mesh refinement."""
    fuel = _make_mixture(
        _XS_A[f"sig_t_{suffix}"], _XS_A[f"sig_c_{suffix}"],
        _XS_A[f"sig_f_{suffix}"], _XS_A[f"nu_{suffix}"],
        _XS_A[f"chi_{suffix}"], _XS_A[f"sig_s_{suffix}"],
    )
    mod = _make_mixture(
        _XS_B[f"sig_t_{suffix}"], _XS_B[f"sig_c_{suffix}"],
        _XS_B[f"sig_f_{suffix}"], _XS_B[f"nu_{suffix}"],
        _XS_B[f"chi_{suffix}"], _XS_B[f"sig_s_{suffix}"],
    )
    materials = {2: fuel, 0: mod}

    keffs = []
    deltas = [0.1, 0.05, 0.02]
    for delta in deltas:
        n_fuel = max(2, round(0.5 / delta))
        n_mod = max(2, round(0.5 / delta))
        geom = _slab_geom_for_do(n_fuel, n_mod, delta)
        result = solve_discrete_ordinates(
            materials, geom,
            params=DOParams(max_outer=300, bicgstab_tol=1e-6),
        )
        keffs.append(result.keff)

    # The sequence should converge (decreasing error between successive values)
    diffs = [abs(keffs[i] - keffs[i + 1]) for i in range(len(keffs) - 1)]
    assert diffs[-1] < diffs[0], (
        f"SN not converging: diffs={diffs}"
    )

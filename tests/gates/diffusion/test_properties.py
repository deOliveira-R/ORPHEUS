"""Structural property gates for the modern diffusion solver (L0).

Three defining properties of the converged fundamental mode:

- **Marshak vacuum semantics** — ``BC("vacuum")`` MEANS zero incoming
  partial current (:math:`J^- = 0`, albedo :math:`\\mathcal A = 0`;
  #290 ruling 3). Sharpened from the retired island test's weak
  "boundary flux < 0.3·peak" proxy to the law itself, read off the
  typed :math:`(J^+, J^-)` trace at LU exactness.
- **Positivity** — the loss operator is an M-matrix, so the
  fundamental mode is strictly positive (Perron–Frobenius).
- **Symmetry** — a uniform slab under symmetric BCs has a mirror-
  symmetric flux and equal outward net currents on both faces.
"""

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import mixture_from_diffusion_tables
from orpheus.diffusion import solve_diffusion_1d
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D

pytestmark = pytest.mark.l0  # Diffusion property checks (BC law, positivity, symmetry)


# CORE1D.m fuel column (the ruling-4 legacy encoding — see
# ``mixture_from_diffusion_tables``).
_CORE1D_FUEL = dict(
    transport=np.array([0.2181, 0.7850]),
    absorption=np.array([0.0096, 0.0959]),
    fission=np.array([0.0024, 0.0489]),
    production=np.array([0.0061, 0.1211]),
    chi=np.array([1.0, 0.0]),
    scattering=np.array([0.0160, 0.0]),
)


def _bare_slab(bc: BC, n_cells: int = 20, height: float = 50.0):
    """Uniform bare fuel slab under the same BC on both faces."""
    edges = np.linspace(0.0, height, n_cells + 1)
    mesh = Mesh1D(
        edges=edges, mat_ids=np.zeros(n_cells, dtype=int),
        bc_left=bc, bc_right=bc,
    )
    materials = {0: mixture_from_diffusion_tables(_CORE1D_FUEL)}
    return solve_diffusion_1d(
        materials, mesh, keff_tol=1e-12, flux_tol=1e-10, max_outer=3000,
    )


def test_vacuum_means_zero_incoming_current():
    """``BC("vacuum")`` is the Marshak law: J⁻ = 0 at the solution.

    The inflow rows of the composite system are constraint rows with
    zero source, so the converged trace satisfies the law to LU
    roundoff — asserted at 1e-12 of the trace scale, not a loose
    heuristic. The flux at the boundary CELLS stays strictly positive:
    vacuum is *not* the zero-flux Dirichlet idealization.
    """
    result = _bare_slab(BC("vacuum"))
    trace = result.flux.boundary
    scale = float(np.abs(trace.values).max())
    for face in ("xmin", "xmax"):
        np.testing.assert_allclose(
            trace.inflow_view(face), 0.0, atol=1e-12 * scale,
        )
    assert np.all(result.flux.interior.values[:, [0, -1]] > 0.0), (
        "Vacuum (Marshak) must NOT pin the boundary-cell flux to zero"
    )


def test_flux_positivity():
    """All flux values must be positive in the fundamental mode."""
    result = _bare_slab(BC("zero_flux"))
    phi = result.flux.interior.values
    assert np.all(phi > 0), f"Non-positive flux: min={phi.min():.6e}"


def test_flux_symmetry():
    """For a symmetric bare slab, flux must be mirror-symmetric and the
    outward net currents on the two faces must agree."""
    result = _bare_slab(BC("zero_flux"))
    phi = result.flux.interior.values  # (2, n_cells)
    for g in range(2):
        np.testing.assert_allclose(
            phi[g, :], phi[g, ::-1], rtol=1e-10,
            err_msg=f"Group {g} flux is not symmetric",
        )
    trace = result.flux.boundary
    np.testing.assert_allclose(
        trace.net_current("xmin"), trace.net_current("xmax"), rtol=1e-10,
        err_msg="Outward net currents differ across the mirror plane",
    )

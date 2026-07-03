r"""MMS gate for the diffusion loss operator — heterogeneous D, coupled groups (#93).

Method-of-Manufactured-Solutions verification of the FIXED-SOURCE
diffusion problem :math:`A\psi = q` with :math:`A = L + C - S - B`:
a smooth manufactured :math:`\phi_g(x)` is substituted into the
CONTINUOUS operator (SymPy differentiation — structurally independent
of the finite-difference assembly under test), the resulting forcing
is fed through the solver's exact resolvent, and the discrete solution
must converge to the manufactured field at :math:`\mathcal O(h^2)`.

**Ansatz design (vv Mode-7 declaration).** Issue #93's original
proposal — a single-group sine with constant :math:`D` — NULLS the two
hardest terms by construction: with :math:`D' \equiv 0` the
face-interpolated conductance never differs from the cell value, and
with one group there is no scatter coupling. This gate therefore uses:

- **heterogeneous per-group** :math:`D_g(x)` (smooth, :math:`D' \neq 0`
  across the whole domain, distinct profiles per group) — activates the
  :math:`(D\phi')'` product-rule content and the per-face conductance
  interpolation at EVERY interior face (every cell is its own material);
- **multigroup-coupled distinct shapes** (:math:`\phi_2 \not\propto
  \phi_1`) — activates the down-scatter in-scatter row with an
  :math:`\mathcal O(1)` mismatch against any transpose/sign confusion.

Terms NULLED by this ansatz, each with its covering gate elsewhere:
fission (absent — covered by the L1 eigenvalue anchors and the L2
infinite-medium gate), the non-zero-flux boundary laws (covered by the
per-law trace-semantics gates in ``test_solver.py``), curvilinear
area/volume factors (slab — pinned by the P4 hand-posed stencil gate).

**Constrained, not merely activated (vv Mode-10).** The companion
``test_mms_gate_sees_the_hard_terms`` COMMITS the evidence that the
error functional responds :math:`\mathcal O(1)` — not sub-floor — to a
corruption of each declared-activated term (heterogeneous-D flattened;
in-scatter forcing sign flipped), so a production sign/interpolation
bug in those terms cannot hide under the convergence floor.

See :ref:`the theory-page MMS section <diffusion-mms-section>` and
equation :math:numref:`diffusion-mms`.
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.diffusion import DiffusionSolver
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.material_mesh import MaterialMesh

pytestmark = [pytest.mark.l1, pytest.mark.verifies("diffusion-mms")]


# ── The manufactured problem (module constants = the SSOT the theory
#    page's diffusion-mms equation transcribes) ─────────────────────────

_L = 10.0                          # slab length (cm)
_SIG_A = np.array([0.010, 0.080])  # absorption per group
_SIG_12 = 0.015                    # down-scatter 1 → 2


def _symbolic_problem():
    r"""Manufactured :math:`(D_g, \phi_g, q_g)` with SymPy-exact forcing.

    .. math::

        q_1 = -(D_1\phi_1')' + (\Sigma_{a,1} + \Sigma_{1\to2})\,\phi_1,
        \qquad
        q_2 = -(D_2\phi_2')' + \Sigma_{a,2}\,\phi_2
              - \Sigma_{1\to2}\,\phi_1

    (the in-group scatter cancels against the collision term by the
    P4 column-sum theorem, leaving removal = absorption + out-scatter).
    """
    x = sp.symbols("x", real=True)
    L = sp.Float(_L)
    D1 = 1.2 * (1 + 0.35 * sp.sin(sp.pi * x / L))
    D2 = 0.45 * (1 - 0.25 * sp.cos(2 * sp.pi * x / L))
    phi1 = sp.sin(sp.pi * x / L)
    phi2 = sp.sin(2 * sp.pi * x / L) + 0.6 * sp.sin(sp.pi * x / L)
    q1 = -sp.diff(D1 * sp.diff(phi1, x), x) + (_SIG_A[0] + _SIG_12) * phi1
    q2 = -sp.diff(D2 * sp.diff(phi2, x), x) + _SIG_A[1] * phi2 - _SIG_12 * phi1
    return tuple(
        sp.lambdify(x, expr, "numpy") for expr in (D1, D2, phi1, phi2, q1, q2)
    )


_D1, _D2, _PHI1, _PHI2, _Q1, _Q2 = _symbolic_problem()


def _materials_from_d(d1: np.ndarray, d2: np.ndarray) -> dict:
    r"""Per-cell Mixtures with :math:`\sigma_t = 1/(3 D)` (the P1 seam
    inverted), fixed absorption, chained down-scatter, in-group
    backfill — every cell is its own material, so EVERY interior face
    is a material interface for the conductance interpolation."""
    sig_t = np.vstack([1.0 / (3.0 * d1), 1.0 / (3.0 * d2)])
    materials = {}
    for i in range(sig_t.shape[1]):
        st = sig_t[:, i]
        sig_s = np.array([
            [st[0] - _SIG_A[0] - _SIG_12, _SIG_12],
            [0.0, st[1] - _SIG_A[1]],
        ])
        materials[i] = make_mixture(
            sig_t=st.copy(), sig_c=_SIG_A.copy(), sig_f=np.zeros(2),
            nu=np.zeros(2), chi=np.zeros(2), sig_s=sig_s,
        )
    return materials


def _solve_mms(
    n_cells: int,
    *,
    flatten_d: bool = False,
    flip_inscatter_forcing: bool = False,
) -> float:
    r"""Solve :math:`A\psi = q_{\rm mms}` on ``n_cells`` and return the
    cell-width-weighted :math:`\ell^2` error against the manufactured
    flux.

    The two keyword mutations build the Mode-10 controls: ``flatten_d``
    solves with the WRONG operator (constant midslab :math:`D` while the
    forcing keeps the heterogeneous one); ``flip_inscatter_forcing``
    corrupts the coupling term's sign in the forcing.
    """
    edges = np.linspace(0.0, _L, n_cells + 1)
    z = 0.5 * (edges[:-1] + edges[1:])

    if flatten_d:
        d1 = np.full_like(z, float(_D1(_L / 2)))
        d2 = np.full_like(z, float(_D2(_L / 2)))
    else:
        d1, d2 = np.asarray(_D1(z), float), np.asarray(_D2(z), float)

    mesh = Mesh1D(
        edges=edges, mat_ids=np.arange(n_cells),
        bc_left=BC("zero_flux"), bc_right=BC("zero_flux"),
    )
    solver = DiffusionSolver(MaterialMesh(mesh, _materials_from_d(d1, d2)))

    sign = +1.0 if flip_inscatter_forcing else -1.0
    q_bulk = np.vstack([
        np.asarray(_Q1(z), float),
        # q2 = -(D2 φ2')' + Σa2 φ2 + sign·Σ12 φ1  (sign = -1 is the truth)
        np.asarray(_Q2(z), float) + (sign + 1.0) * _SIG_12 * np.asarray(_PHI1(z), float),
    ])

    n_flat = int(np.asarray(solver.template.to_flat()).size)
    q_flat = np.concatenate([q_bulk.ravel(), np.zeros(n_flat - q_bulk.size)])
    # Layout self-check: the typed composite must read our bulk block
    # back verbatim and a zero trace block (guards the bulk-first flat
    # convention this construction rides on).
    q_full = FullField.from_flat(q_flat, solver.template)
    np.testing.assert_array_equal(q_full.bulk.values, q_bulk)
    np.testing.assert_array_equal(
        q_full.boundary.values, np.zeros_like(q_full.boundary.values),
    )

    psi = solver.unflatten(solver.solve_fixed_source(q_flat, q_flat))
    phi_exact = np.vstack([
        np.asarray(_PHI1(z), float), np.asarray(_PHI2(z), float),
    ])
    dz = _L / n_cells
    return float(np.sqrt(dz * np.sum((psi.bulk.values - phi_exact) ** 2)))


def test_mms_flux_converges_second_order():
    r"""The discrete solution must converge to the manufactured flux at
    the design :math:`\mathcal O(h^2)` of the cell-centered FD scheme,
    with the D-gradient and group-coupling terms ACTIVE (Mode-7
    declaration in the module docstring)."""
    ns = [20, 40, 80, 160]
    errs = np.asarray([_solve_mms(n) for n in ns])
    orders = np.log2(errs[:-1] / errs[1:])

    assert np.all(orders > 1.8), (
        f"MMS convergence below O(h²): errors={errs}, orders={orders}"
    )
    assert errs[-1] < 1e-3, (
        f"Finest-mesh MMS error {errs[-1]:.2e} above the expected floor"
    )


def test_mms_gate_sees_the_hard_terms():
    r"""Mode-10 committed controls: the error functional must respond
    :math:`\mathcal O(1)` to a corruption of each declared-activated
    term — the terms are CONSTRAINED, not merely exercised.

    A clean n=40 solve sits on the :math:`\mathcal O(h^2)` floor; the
    flattened-D operator and the sign-flipped in-scatter forcing must
    each blow the error far above it (×30 is orders of magnitude below
    the observed separation — the assertion is a tripwire, not a fit).
    """
    err_clean = _solve_mms(40)
    err_flat_d = _solve_mms(40, flatten_d=True)
    err_flip = _solve_mms(40, flip_inscatter_forcing=True)

    assert err_flat_d > 30 * err_clean, (
        f"Flattening D(x) barely moved the MMS error "
        f"({err_flat_d:.3e} vs clean {err_clean:.3e}) — the D-gradient "
        f"term is not constrained by this gate"
    )
    assert err_flip > 30 * err_clean, (
        f"Flipping the in-scatter forcing sign barely moved the MMS "
        f"error ({err_flip:.3e} vs clean {err_clean:.3e}) — the group "
        f"coupling is not constrained by this gate"
    )

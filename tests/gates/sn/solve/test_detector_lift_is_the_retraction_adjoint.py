r"""The SN adjoint entry lifts a detector by the retraction's adjoint (#405 P1 step 6, S6.6).

A detector response :math:`\Sigma_d(\vec r, g)` is the functional
:math:`\psi \mapsto \langle \Sigma_d, R\psi\rangle`, with :math:`R` the angular
retraction (fibre integration). Its Riesz representative in phase space is
:math:`R^\dagger \Sigma_d`, the pullback
:class:`~orpheus.numerics.operator.AxisPullbackOperator`. A source rate is
lifted by the section :math:`E` instead; the two differ by the mass of the
angular measure, which is never typed (the user's ruling of 2026-10-02).

Until step 6 ``solve_sn_adjoint_fixed_source`` wrote the lift by hand,
``np.broadcast_to(sigma_d[None], …)``: R†'s content, but not R†. This file
gates the ROUTE, which the value gates cannot see (the hand broadcast and the
pullback are bit-identical, ``[M]`` 2026-10-02: φ* and ψ* ``array_equal``
before and after the re-spelling on this fixture and on a random detector):

- a counting spy on the pullback reads exactly 1 call per solve with an
  ndarray detector, and 0 with a composite (``FullField``) detector, which
  arrives already in phase space;
- a decoy pullback scaled by 2 moves φ* to exactly 2·φ* (the solve is linear
  in its right-hand side), so the entry's answer depends on the pullback.

First red (2026-10-02, the test-architect's ``route_probe.out`` T0): on the
hand broadcast the spy read 0 calls and the decoy left φ* unmoved.

The value of the lift is pinned independently by
``test_sn_adjoint_entries.py::TestSolveSnAdjointFixedSource::test_duality_cross_group_source_detector``
(the hand volume sum Σ V·Σ_d·φ and the hand-built pairing), which reddens
when the lift is the section E or the plain transpose Rᵀ.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.operator import AxisPullbackOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import _as_problem, solve_sn_adjoint_fixed_source
from orpheus.transport.full_field import FullField
from orpheus.transport.source_sinks import AngularBoundarySourceSink, AngularSourceSink

pytestmark = pytest.mark.l1


def _fixture():
    """The P1.2 duality fixture: a 2-material vacuum slab, 2 groups, GL8, and
    a thermal detector in the right region."""
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesher(StructuredGeometry.slab(
        (0.0, 1.0, 3.0, 4.0), (0, 1, 0), left=BC.vacuum, right=BC.vacuum,
    )).partition((
        CellsByCount.uniform_width(2), CellsByCount.uniform_width(4), CellsByCount.uniform_width(2),
    )).mesh
    sigma_d = np.zeros((2, 8))
    sigma_d[1, 6:8] = 0.7
    return mats, mesh, Quadrature.gauss_legendre(n_ordinates=8), sigma_d


def _phi_star(mats, mesh, quad, detector) -> np.ndarray:
    adj = solve_sn_adjoint_fixed_source(mats, mesh, quad, detector, inner_tol=1e-12, max_inner=2000)
    return np.asarray(adj.scalar_flux.values)


def test_s6_6_a_the_ndarray_detector_routes_through_the_pullback(monkeypatch) -> None:
    mats, mesh, quad, sigma_d = _fixture()
    calls = []
    honest = AxisPullbackOperator.apply

    def counting(self, x):
        calls.append(1)
        return honest(self, x)

    monkeypatch.setattr(AxisPullbackOperator, "apply", counting)
    _phi_star(mats, mesh, quad, sigma_d)
    if len(calls) != 1:
        pytest.fail(f"an ndarray detector made {len(calls)} pullback calls, not 1")

    calls.clear()
    problem = _as_problem(mesh, quad, mats)
    composite = FullField(
        interior=AngularSourceSink(
            values=np.ascontiguousarray(np.broadcast_to(sigma_d[None], (quad.N, *sigma_d.shape))),
            space=problem.angular_trial_space,
        ),
        boundary=AngularBoundarySourceSink.zeros(problem.angular_trace),
    )
    _phi_star(mats, mesh, quad, composite)
    if calls:
        pytest.fail(f"a composite detector made {len(calls)} pullback calls; it is already in phase space")


def test_s6_6_b_the_answer_depends_on_the_pullback(monkeypatch) -> None:
    mats, mesh, quad, sigma_d = _fixture()
    honest_phi = _phi_star(mats, mesh, quad, sigma_d)
    honest = AxisPullbackOperator.apply
    monkeypatch.setattr(AxisPullbackOperator, "apply", lambda self, x: 2.0 * honest(self, x))
    npt.assert_array_equal(_phi_star(mats, mesh, quad, sigma_d), 2.0 * honest_phi)

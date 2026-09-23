r"""The within-group loss representations on a 2-D Cartesian mesh: cost and agreement.

Three strategies realise the same operator :math:`(L+C)` on a multi-D
Cartesian mesh (``orpheus.sn.loss_representation``): ``ScanMarch`` (the
row-march default since the Fork-B2 ruling of #222), ``MovingFrontierWindow``
(the anti-diagonal wavefront, a selectable peer) and ``FullFieldWavefront``
(the full-field spine). This suite times the sweep :math:`(L+C)^{-1}q` and the
matvec :math:`(L+C)\psi` of each over a size, angular-order and group grid,
records the sweep's peak memory, and records beside the cost the sweep's
disagreement with the ``FullFieldWavefront`` oracle, so that a speed-up is
never read without its accuracy. The window's disagreement reads exactly 0
because the window and the oracle share their base class
(``_DAGWavefront``) and its arithmetic: that zero is one implementation
agreeing with itself, not independent evidence; ``ScanMarch``'s 1e-15 is
the comparison between two different walks.

It rests on the verification of the capability it times: the capability tier
``sweep_cartesian_2d`` (the ``cap`` marker) and, for agreement between the
strategies, ``tests/gates/sn/solve/test_scan_march_end_to_end.py``. A case is
a measurement, never a gate: `docs/theory/verification/principles.rst`,
"Where a case lives".

Provenance: the measurement behind the representation table in
``docs/theory/methods/sn/loss_representation.rst`` (#222), first a diagnostic
script (``git show f36572c8^:derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py``).
"""
from __future__ import annotations

from contextlib import contextmanager
from unittest.mock import patch

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import loss_representation as lr
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.timed_full_field import TimedFullField

CAPABILITY = "sweep_cartesian_2d"

REPRESENTATIONS = {
    "scanmarch": lr.ScanMarch,
    "window": lr.MovingFrontierWindow,
    "fullfield": lr.FullFieldWavefront,
}

#: (nx, ny, level-symmetric order, groups): mesh size, angular order and groups
#: varied one at a time about the 48 x 48, LS-4, 2-group centre.
GRID = ["24x24-LS4-2g", "48x48-LS4-2g", "96x96-LS4-2g", "48x48-LS8-2g", "48x48-LS16-2g", "48x48-LS8-4g"]


def _parse(config: str) -> tuple[int, int, int, int]:
    size, order, groups = config.split("-")
    nx, ny = (int(n) for n in size.split("x"))
    return nx, ny, int(order.removeprefix("LS")), int(groups.removesuffix("g"))


def _problem(config: str) -> tuple[SNProblem, np.ndarray, np.ndarray]:
    nx, ny, order, ng = _parse(config)
    reflective = BC("reflective")
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1), edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int), coord=CoordSystem.CARTESIAN,
        bc_xmin=reflective, bc_xmax=reflective, bc_ymin=reflective, bc_ymax=reflective,
    )
    problem = SNProblem(mesh, Quadrature.level_symmetric(order), {0: get_mixture("A", f"{ng}g")})
    rng = np.random.default_rng(7)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx, ny))
    source = rng.uniform(0.0, 2.0, size=(problem.quad.N, ng, nx, ny))
    return problem, sig_t, source


def _representation(name: str, problem: SNProblem) -> lr.LossRepresentation:
    return REPRESENTATIONS[name](problem, problem.scheme, problem.angular_closure)


class LossRepresentationKernels:
    """Sweep and matvec cost of each representation, with the sweep's agreement."""

    params = (GRID, list(REPRESENTATIONS))
    param_names = ("config", "representation")
    timeout = 300

    def setup(self, config: str, representation: str) -> None:
        self.problem, self.sig_t, self.source = _problem(config)
        self.rep = _representation(representation, self.problem)
        self.stratum = self.rep.bind_sigma(self.sig_t)
        self.psi = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, space=self.problem.full_field_space,
        )
        self.psi.interior.values[...] = np.random.default_rng(3).standard_normal(
            self.psi.interior.values.shape,
        )

    def _sweep(self) -> np.ndarray:
        angular, _ = self.rep.sweep(
            self.source, self.stratum, AngularBoundaryFlux.zeros(self.problem.angular_trace),
        )
        return angular

    def time_sweep(self, config: str, representation: str) -> None:
        self._sweep()

    def time_loss_action(self, config: str, representation: str) -> None:
        self.rep.loss_action(self.sig_t, self.psi)

    def peakmem_sweep(self, config: str, representation: str) -> None:
        self._sweep()

    def track_sweep_disagreement_with_fullfield(self, config: str, representation: str) -> float:
        """max |psi - psi_oracle| / max |psi_oracle| against ``FullFieldWavefront``."""
        oracle_rep = _representation("fullfield", self.problem)
        oracle, _ = oracle_rep.sweep(
            self.source, oracle_rep.bind_sigma(self.sig_t),
            AngularBoundaryFlux.zeros(self.problem.angular_trace),
        )
        return float(np.max(np.abs(self._sweep() - oracle)) / np.max(np.abs(oracle)))

    setattr(track_sweep_disagreement_with_fullfield, "unit", "relative")  # asv reads it


@contextmanager
def _window_forced():
    """The window wherever the selector would pick ``ScanMarch`` on a multi-D mesh.

    Refuses to report a timing for a forcing that never took effect: a forced
    leg that fell back to the default would time the default twice.
    """
    real = lr.default_for
    forced_count = 0

    def forced(problem, spatial_closure, angular_closure):
        nonlocal forced_count
        rep = real(problem, spatial_closure, angular_closure)
        if isinstance(rep, lr.ScanMarch) and not problem.is_1d:
            forced_count += 1
            return lr.MovingFrontierWindow(problem, spatial_closure, angular_closure)
        return rep

    with patch.object(lr, "default_for", forced):
        yield
    if forced_count == 0:
        raise RuntimeError("the window was never selected: the forced leg timed the default")


class FixedSourceSolve:
    """One heterogeneous fixed-source solve end to end, default against the forced window."""

    params = (["scanmarch-default", "window-forced"],)
    param_names = ("representation",)
    timeout = 600

    def setup(self, representation: str) -> None:
        nx = ny = 48
        mat = np.zeros((nx, ny), dtype=int)
        mat[: nx // 2, :] = 2
        self.mesh = Mesh2D(
            edges_x=np.linspace(0.0, 12.0, nx + 1), edges_y=np.linspace(0.0, 12.0, ny + 1),
            mat_map=mat,
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
        )
        self.quad = Quadrature.level_symmetric(8)
        self.materials = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
        self.source = np.full((self.quad.N, 2, nx, ny), 1.0 / float(self.quad.weights.sum()))

    def _solve(self, representation: str):
        def run():
            return solve_sn_fixed_source(
                self.materials, self.mesh, self.quad, self.source, scattering_order=1,
                inner_solver="source_iteration", max_inner=2000, inner_tol=1e-10,
            )
        if representation == "window-forced":
            with _window_forced():
                return run()
        return run()

    def time_solve(self, representation: str) -> None:
        self._solve(representation)

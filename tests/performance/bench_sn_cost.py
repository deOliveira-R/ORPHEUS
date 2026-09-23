r"""SN cost measurements that used to be timing gates in ``tests/gates/``.

The user's ruling (2026-09-23, plan ``.claude/plans/vv_suite_layout.md``, R6):
a wall-clock time is never a gate on a machine whose load is not
reproducible, and a slower correct answer never blocks publishing. So every
timing assertion left the gates and is measured here, where a regression is
visible in the history without failing anything. The exactly reproducible cost
checks (call counts, allocation bytes) stay gates in
``tests/gates/sn/architecture/test_composition_cost.py``.

Each case is self-contained on the public ``orpheus`` API: asv installs only
the package when it builds a commit, so a case cannot import the gates'
builders.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.loss_representation import CumprodScan, FullFieldWavefront
from orpheus.sn.problem import SNProblem
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.source_sinks import AngularSourceSink


def _graded_edges(length: float, n_cells: int, *, ratio: float) -> np.ndarray:
    widths = ratio ** np.arange(n_cells, dtype=float)
    return np.concatenate(([0.0], np.cumsum(widths))) * (length / widths.sum())


class CompositionOverhead:
    """One posed ``A.apply(x)`` against a dense contraction of the same nominal FLOP count.

    Formerly the P-2 leg of the composition-cost gate. The ratio divides out
    the host's raw numpy throughput; the call-count leg (P-1) and the
    allocation leg (P-3) remain gates. Rests on the capability
    ``sweep_cartesian_2d``. The fixture follows the gate's 32 x 40 graded,
    two-region, mixed-boundary mesh, with the shipped 2-group mixtures A and
    B (P0) in place of the gate's hand-built P1 pair.
    """

    #: Nominal FLOPs per degree of freedom for a DD cell balance: it sets the
    #: calibration size, so changing it changes what the ratio means.
    flops_per_dof = 16
    timeout = 300

    def setup(self) -> None:
        nx, ny = 32, 40
        mat_map = np.zeros((nx, ny), dtype=int)
        mat_map[nx // 2:, :] = 1
        mesh = Mesh2D(
            edges_x=_graded_edges(2.0, nx, ratio=1.02), edges_y=_graded_edges(3.0, ny, ratio=1.03),
            mat_map=mat_map, coord=CoordSystem.CARTESIAN,
            bc_xmin=BC("reflective"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("reflective"), bc_ymax=BC("vacuum"),
        )
        problem = SNProblem(
            mesh, Quadrature.level_symmetric(sn_order=4),
            {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")},
        )
        record = build_within_group_system(problem, problem.mat_xs)
        self.A = record.loss
        template = record.space.zeros()
        rng = np.random.default_rng([20260728, 7])
        self.x = type(template).from_flat(rng.standard_normal(template.to_flat().size), template)
        side = int(round(np.sqrt(self.flops_per_dof * self.x.to_flat().size)))
        self.dense = rng.standard_normal((side, side))
        self.vector = rng.standard_normal(side)
        self.A.apply(self.x)

    def time_apply(self) -> None:
        self.A.apply(self.x)

    def time_calibration_contraction(self) -> None:
        np.einsum("ij,j->i", self.dense, self.vector)

    def track_apply_over_calibration(self) -> float:
        """min over repeats of the apply time over min of the contraction time."""
        import time

        def best(fn, n):
            fn()
            out = float("inf")
            for _ in range(n):
                t0 = time.perf_counter()
                fn()
                out = min(out, time.perf_counter() - t0)
            return out

        return best(lambda: self.A.apply(self.x), 7) / best(
            lambda: np.einsum("ij,j->i", self.dense, self.vector), 25,
        )

    setattr(track_apply_over_calibration, "unit", "ratio")  # asv reads it


class SlabSweep:
    """One slab sweep, nx = 160, Gauss-Legendre S16, 4 groups (formerly a 2 ms gate)."""

    def setup(self) -> None:
        nx, ng = 160, 4
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
            bc_left=BC("vacuum"), bc_right=BC("vacuum"),
        )
        self.problem = SNProblem(mesh, Quadrature.gauss_legendre(16), {0: get_mixture("A", f"{ng}g")})
        self.rep = CumprodScan.pose(self.problem)
        self.stratum = self.rep.bind_sigma(np.ones((ng, nx)))
        self.source = AngularSourceSink.from_isotropic(np.ones((ng, nx)), self.problem).values
        self.rep.sweep(self.source, self.stratum, AngularBoundaryFlux.zeros(self.problem.angular_trace))

    def time_sweep(self) -> None:
        self.rep.sweep(self.source, self.stratum, AngularBoundaryFlux.zeros(self.problem.angular_trace))


class CumprodAgainstSpine:
    """The d = 1 cumprod scan against the full-field spine on a long chain (nx = 4096).

    Formerly a speed-up tripwire: the reason cumprod is the d = 1 default.
    """

    params = (["cumprod", "spine"],)
    param_names = ("representation",)

    def setup(self, representation: str) -> None:
        nx = 4096
        mesh = Mesh1D(
            edges=np.linspace(0.0, 10.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
            coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
        )
        self.problem = SNProblem(mesh, Quadrature.gauss_legendre(n_ordinates=8), {0: get_mixture("A", "2g")})
        rng = np.random.default_rng(99)
        sig_t = rng.uniform(0.3, 3.0, size=(self.problem.ng, nx))
        iso = rng.uniform(0.2, 1.5, size=(self.problem.ng, nx))
        self.source = AngularSourceSink.from_isotropic(iso, self.problem).values
        cls = CumprodScan if representation == "cumprod" else FullFieldWavefront
        self.rep = cls.pose(self.problem)
        self.stratum = self.rep.bind_sigma(sig_t)
        self.rep.sweep(self.source, self.stratum, AngularBoundaryFlux.zeros(self.problem.angular_trace))

    def time_sweep(self, representation: str) -> None:
        self.rep.sweep(self.source, self.stratum, AngularBoundaryFlux.zeros(self.problem.angular_trace))


class SolverComponents:
    """The solver's per-component costs on the 6 x 4 fuel|moderator 2-group Lebedev-17 case.

    Formerly ``TestPerformanceBaseline::test_profile_components``, a gate that
    printed timings and asserted nothing.
    """

    def setup(self) -> None:
        nx, ny, delta = 6, 4, 0.2
        mat = np.zeros((nx, ny), dtype=int)
        mat[:3, :] = 2
        mesh = Mesh2D(
            edges_x=np.linspace(0, nx * delta, nx + 1), edges_y=np.linspace(0, ny * delta, ny + 1),
            mat_map=mat,
        )
        problem = SNProblem(mesh, Quadrature.lebedev(order=17), {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")})
        self.solver = SNSolver(problem)
        rng = np.random.default_rng(42)
        self.phi = rng.random((self.solver.ng, *problem.spatial_shape)) + 0.1
        self.Q = rng.random((self.solver.ng, *problem.spatial_shape))
        self.fission_source = self.solver.compute_fission_source(self.phi, 1.0)

    def time_p0_scattering_source(self) -> None:
        self.solver.problem.system.factors.scattering.transfer.add_p0_source(self.Q.copy(), self.phi)

    def time_compute_keff(self) -> None:
        self.solver.compute_keff(self.phi)

    def time_solve_fixed_source_one_outer(self) -> None:
        self.solver.solve_fixed_source(self.fission_source, self.phi)

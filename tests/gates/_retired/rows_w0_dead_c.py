"""The rows of unit ``dead`` sub-unit dead-C (the sn package; specification ``ontology_w0_spec.md`` Appendix D-C.7).

22 rows of Appendix D-C.1 as 1 ``RetiredModule`` (the whole-module delete of ``orpheus/sn/loss_representation/assembly.py``, rows 85, 86,
141, 142, 143), 16 ``Moved`` rows (to the test side), 6 ``Retired`` rows (row 92 is a delegate of row 37, so ``Retired`` with a rewire,
decision D1 of the sixth pass); ``KEPT`` holds rows 83 and 100 and the deferred row 99 (decision D2, private state). Row 126 is dead-D's.

First red at the base tip: every ``Retired`` and ``Moved`` row is defined today and the ``RetiredModule`` file exists. The ``Kept`` entries
are green today and must stay green.

Exclusion of the reference leg (done-when 4 resolves a reference by import): it cannot see an INSTANCE-ATTRIBUTE read of a retired
member (``m.loss_kernel_gauge``, ``m.boundary_face_layout``, ``m.n_unknowns_flat``, ``s.compute_group_absorption_rate``,
``term.is_gain``). The only catcher of such a read is the ``AttributeError`` raised when the reading test runs (``[M]``: planting
``layout = m.boundary_face_layout`` in a copy of ``test_boundary_face_layout.py`` turned 0 rows red).
"""

from __future__ import annotations

from tests.gates._retired import Kept, Moved, Retired, RetiredModule, Row
ROWS: tuple[Row, ...] = (
    RetiredModule('orpheus.sn.loss_representation.assembly'),  # the whole-module delete of rows 85, 86, 141, 142, 143
    Moved('orpheus/mesh/axis.py', 'n_unknowns_flat', 'tests/gates/sn/primitives/test_axis_primitive.py', 'n_unknowns_flat'),  # row 37; free function n_unknowns_flat(axes, quad, ng)
    Retired('orpheus/sn/angular/closure.py', '_MMHalfGrid.downstream_per_ordinate'),  # row 80; free function downstream_per_ordinate(half_grid)  # was a property; moved to the test side 2026-11-20, retired 2026-11-21 (no consumer)
    Retired('orpheus/sn/angular/closure.py', '_MMHalfGrid.upstream'),  # row 81; free function upstream(half_grid, m); moved to the test side 2026-11-20, retired 2026-11-21 (no consumer)
    Moved('orpheus/sn/angular/closure.py', 'compute_psi_half_per_level', 'tests/gates/sn/sweep/curvilinear/_psi_half.py', 'compute_psi_half_per_level'),  # row 82; free function compute_psi_half_per_level(psi_level, tau_level, *, psi_half_seed=None)
    Retired('orpheus/sn/coupled_system.py', 'build_coupled_system'),
    Moved('orpheus/sn/loss_representation/assembly.py', 'ordinate_walk_order', 'tests/gates/sn/sweep/_ordinate_assembly.py', 'ordinate_walk_order'),  # row 85; free function ordinate_walk_order(problem, ordinate)
    Moved('orpheus/sn/loss_representation/assembly.py', 'assemble_ordinate_blocks', 'tests/gates/sn/sweep/_ordinate_assembly.py', 'assemble_ordinate_blocks'),  # row 86; free function assemble_ordinate_blocks(problem, ordinate, *, include_collision=True)
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_in_x', 'tests/gates/sn/sweep/_sweep_graph_faces.py', 'face_in_x'),  # row 87; free function face_in_x(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_out_x', 'tests/gates/sn/sweep/_sweep_graph_faces.py', 'face_out_x'),  # row 88; free function face_out_x(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_in_y', 'tests/gates/sn/sweep/_sweep_graph_faces.py', 'face_in_y'),  # row 89; free function face_in_y(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/loss_representation/sweep_graph.py', 'SweepDependencyGraph.face_out_y', 'tests/gates/sn/sweep/_sweep_graph_faces.py', 'face_out_y'),  # row 90; free function face_out_y(sweep_dependency_graph)  # was a property
    Moved('orpheus/sn/operators/loss_kernel_gauge.py', 'predicted_kernel_dimension', 'tests/gates/sn/operators/_kernel_dimension.py', 'predicted_kernel_dimension'),  # row 91; free function predicted_kernel_dimension(problem)
    Retired('orpheus/sn/problem.py', 'SNProblem.n_unknowns_flat'),  # row 92; a DELEGATE of row 37 (`SNProblem.n_unknowns_flat` is `return _axis_n_unknowns_flat(self.axes, self.quad, self.ng)`): nothing moves, because no second definition exists; the row is `Retired` with a rewire (decision D1 of the sixth pass): every `p.n_unknowns_flat` is rewritten `n_unknowns_flat(p.axes, p.quad, p.ng)`, a call of row 37's helper (D6)
    Retired('orpheus/sn/problem.py', 'SNProblem.loss_kernel_gauge'),
    Retired('orpheus/sn/problem.py', 'SNProblem.boundary_face_layout'),
    Moved('orpheus/sn/solver.py', 'boundary_vs_interior_split', 'tests/gates/sn/operators/test_typed_residual_evaluation.py', 'boundary_vs_interior_split'),  # row 95; free function boundary_vs_interior_split(residual)
    Retired('orpheus/sn/solver.py', 'SNSolver.compute_group_absorption_rate'),
    Retired('orpheus/sn/splitting.py', 'LossTerm.is_gain'),
    Moved('orpheus/sn/splitting.py', 'Splitting.law_residual', 'tests/gates/sn/architecture/_splitting_oracles.py', 'law_residual'),  # row 98; free function law_residual(splitting, state)
    Moved('orpheus/sn/loss_representation/assembly.py', '_ordinate_label', 'tests/gates/sn/sweep/_ordinate_assembly.py', '_ordinate_label'),  # row 141; free function _ordinate_label(problem, ordinate)
    Moved('orpheus/sn/loss_representation/assembly.py', '_iter_cells_in_walk_order', 'tests/gates/sn/sweep/_ordinate_assembly.py', '_iter_cells_in_walk_order'),  # row 142; free function _iter_cells_in_walk_order(graph)
    Moved('orpheus/sn/loss_representation/assembly.py', '_probe_coefficient_blocks', 'tests/gates/sn/sweep/_ordinate_assembly.py', '_probe_coefficient_blocks'),  # row 143; free function _probe_coefficient_blocks(scheme, s_axes_k, reaction_flat, *, cm, fm, n_cells)
)
KEPT: tuple[Kept, ...] = (  # Kept(file, qualname, keeper): the keeper is the wave or the consumer that keeps the row; asserted still defined
    Kept('orpheus/sn/angular/closure.py', 'MorelMontryAngularSweep.angular_endpoint_defect_per_level', "W7 (#653's radial cases un-defer part of the curvilinear build) / curvilinear SN cluster"),  # row 83
    Kept('orpheus/sn/sweep/pairing.py', 'pair_diffusion_limit_consistent', "W4, concept rows 'Spatial scheme' and 'Angular closure'"),  # row 100
    Kept('orpheus/sn/sweep/cache.py', 'CollisionCache.reset_build_count', 'deferred, private state (decision D2; section 8.3)'),  # row 99
)

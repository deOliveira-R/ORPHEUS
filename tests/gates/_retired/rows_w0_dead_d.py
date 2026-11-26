"""The rows of unit `dead`, sub-unit dead-D (transport): 19 `Retired` and 8 `Moved` rows (Appendix D-D.1 rows 79, 101-103, 105, 107-115, 117-119, 121-126, 129, 138-140) and 8 `Kept` entries
(rows 78, 104, 120, 127, 128, the deferred rows 106 and 116, and the homonym `MaterialXSField.foldable_sigma`).
Each `Moved` home is a test-side file that does not exist until the build; the rows are red at the base by design.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Moved, Retired, Row

ROWS: tuple[Row, ...] = (
    Retired('orpheus/numerics/spaces/scalar_trace_space.py', 'ScalarTraceSpace.slot_shape'),
    Retired('orpheus/transport/boundary_conditions/admission.py', 'Approximated.issues'),
    Retired('orpheus/transport/boundary_conditions/admission.py', 'FaceSite.hollow_inner'),
    Retired('orpheus/transport/boundary_conditions/binding.py', 'Binding.faces_of'),
    Moved('orpheus/transport/fields/_bases.py', 'FaceField.face_views', 'tests/gates/transport/fields/test_angular_boundary_flux.py', 'face_views'),  # row 105; free function face_views(face_field)  # was a property
    Retired('orpheus/transport/fields/angular_source_sink.py', 'AngularSourceSink.at_ordinate'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.l_block'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.isotropic_part'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.anisotropic_part'),
    Retired('orpheus/transport/fields/harmonic_moment_flux.py', 'HarmonicMomentFlux.truncate'),
    Retired('orpheus/transport/fields/harmonic_moment_source_sink.py', 'HarmonicMomentSourceSink._CELL_GROUP_NAME'),
    Moved('orpheus/transport/fields/scalar_boundary_flux.py', 'ScalarBoundaryFlux.p1_boundary_scalar_flux', 'tests/gates/transport/fields/_scalar_boundary_oracles.py', 'p1_boundary_scalar_flux'),  # row 113; free function p1_boundary_scalar_flux(scalar_boundary_flux, face)
    Retired('orpheus/transport/fields/scalar_flux.py', 'ScalarFlux.at_group'),
    Retired('orpheus/transport/kernels.py', 'FissionKernel.dyad'),
    Moved('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.fission_production_per_material', 'tests/gates/transport/test_kernels.py', 'fission_production_per_material'),  # row 117; free function fission_production_per_material(mat_xs, material_id)
    Moved('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.chi_per_material', 'tests/gates/transport/test_kernels.py', 'chi_per_material'),  # row 118; free function chi_per_material(mat_xs, material_id)
    Moved('orpheus/transport/operators/isotropic_transfer.py', 'IsotropicTransfer.dense_per_material', 'tests/gates/transport/_transfer_oracles.py', 'dense_per_material'),  # row 119; free function dense_per_material(isotropic_transfer)
    Moved('orpheus/transport/operators/transfer.py', 'LegendreMomentTransfer.on_basis', 'tests/gates/transport/_transfer_oracles.py', 'on_basis'),  # row 121; free function on_basis(transfer, basis, *, skip_l0=True), its body names `LegendreMomentTransfer`  # was a classmethod: the class parameter is dropped (decision D2)
    Moved('orpheus/transport/operators/transfer.py', 'LegendreMomentTransfer.on_frame', 'tests/gates/transport/_transfer_oracles.py', 'on_frame'),  # row 122; free function on_frame(transfer, frame, *, skip_l0=True), its body names `LegendreMomentTransfer`  # was a classmethod: the class parameter is dropped (decision D2)
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.foldable_part'),
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.residual_part'),
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.is_foldable_into_sigma_r'),
    Retired('orpheus/transport/residuals/angular_residual.py', 'AngularResidual.balance_map'),
    Moved('orpheus/transport/spatial/_ubld.py', 'D1ClosedForm.kernel_rhs', 'tests/gates/transport/spatial/test_ld_ubld_primitive.py', 'kernel_rhs'),  # row 129; free function kernel_rhs(d1_closed_form, Q_cells, psi_in)
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator.foldable_sigma'),
    Retired('orpheus/transport/operators/transfer.py', 'TransferOperator._sibling'),
    Retired('orpheus/numerics/spaces/moment_head.py', 'truncated_head'),
)

KEPT: tuple[Kept, ...] = (  # Kept(file, qualname, keeper): the keeper is the wave or the consumer that keeps the row; asserted still defined
    Kept('orpheus/transport/boundary_conditions/realization.py', 'RealizerTable.catalogue', "door: the method tables' admission matrix, asserted through it by 13 test files (5 methods); nothing in W0 opens the realizer-table cycle (3.4)"),  # row 104
    Kept('orpheus/transport/operators/multiplication_operator.py', 'MultiplicationOperator.from_mesh', 'a documented survivor pinned by a design gate (transport/fields/test_meshless_construction.py:245)'),  # row 120
    Kept('orpheus/transport/residuals/angular_residual.py', 'AngularResidual.relative_to', "W5, concept row 'Iteration record' (one estimator, one norm)"),  # row 127
    Kept('orpheus/transport/residuals/scalar_residual.py', 'ScalarResidual', "W3, concept row 'Flux and source members' (the role forest)"),  # row 128
    Kept('orpheus/transport/fields/_bases.py', 'BoundaryField.from_face_arrays', 'deferred, private state (decision D2; section 8.3)'),  # row 106
    Kept('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.n2n_matrix', 'deferred, private state (decision D2; section 8.3)'),  # row 116
    Kept('orpheus/transport/mesh/material_xs_field.py', 'MaterialXSField.foldable_sigma', 'homonym of row 138, read by sn/acceleration/dsa.py:258 and fenced by test_dsa_rate.py (ERR-070 sentinel)'),  # not a row
)

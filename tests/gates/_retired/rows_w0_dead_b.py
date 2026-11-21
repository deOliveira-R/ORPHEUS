"""The rows of unit ``dead``, sub-unit dead-B (numerics rest and eigenvalue), of the ontology W0 specification (Appendix D-B).

``Retired``: 12 rows (the definition is gone and nothing reads it). ``Moved``: 7 rows (gone from production, defined at the
test-side ``home`` under ``home_name``). ``KEPT``: 8 entries, the 7 ``keep`` rows (38, 47, 48, 53, 62, 71, 75) and the deferred
row 52 (``Field.where_largest``, private state, decision D2), each asserted still defined.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Moved, Retired, Row

ROWS: tuple[Row, ...] = (
    Retired('orpheus/data/emission_spectrum.py', 'EmissionSpectrum.is_emitting'),
    Retired('orpheus/data/macro_xs/mixture.py', 'Mixture.out_scattering_xs'),
    Retired('orpheus/diffusion/method_space.py', 'DiffusionMethodSpace.minimal'),
    Moved('orpheus/numerics/basis/overlap_basis.py', 'OverlapBasis.dominant_column', 'tests/gates/data/test_energy_grid.py', 'dominant_column'),  # row 42; free function dominant_column(overlap_basis)  # was a property
    Moved('orpheus/numerics/basis/overlap_basis.py', 'OverlapBasis.fractional_columns', 'tests/gates/numerics/_overlap_basis_oracles.py', 'fractional_columns'),  # row 43; free function fractional_columns(overlap_basis)  # was a property
    Retired('orpheus/numerics/basis/spherical_harmonic_basis.py', 'SphericalHarmonicBasis.evaluate_from_components'),
    Retired('orpheus/numerics/basis/spherical_harmonic_basis.py', 'MirrorEvenSphericalHarmonicBasis.evaluate_from_components'),
    Moved('orpheus/numerics/coefficient_layout.py', 'TensorLayout.multi_index', 'tests/gates/numerics/test_interval_legendre_basis.py', 'multi_index'),  # row 46; free function multi_index(tensor_layout, flat)
    Moved('orpheus/numerics/eigenvalue.py', 'direct_eigenvalue', 'tests/gates/numerics/_dense_eigen_oracle.py', 'direct_eigenvalue'),  # row 49; free function direct_eigenvalue(A, F, *, imag_tol=1e-9)
    Retired('orpheus/numerics/eigenvalue.py', 'rayleigh_quotient_iteration'),
    Retired('orpheus/numerics/operator.py', 'InverseMetricOperator'),
    Retired('orpheus/numerics/pencil.py', 'OperatorPencil.is_regular'),
    Retired('orpheus/numerics/pencil.py', 'OperatorPencil.rhs_rank'),
    Retired('orpheus/numerics/pencil.py', 'OperatorPencil.apply_at'),
    Moved('orpheus/numerics/space.py', 'FunctionSpace.project_live', 'tests/gates/numerics/_space_oracles.py', 'project_live'),  # row 74; free function project_live(function_space, x)
    Retired('orpheus/numerics/space.py', 'angular_flux_space'),
    Retired('orpheus/numerics/space.py', 'scalar_flux_space'),
    Moved('orpheus/numerics/eigenvalue.py', 'dominant_eigenpair', 'tests/gates/numerics/_dense_eigen_oracle.py', 'dominant_eigenpair'),  # row 130; free function dominant_eigenpair(M, *, imag_tol=1e-9)
    Moved('orpheus/numerics/eigenvalue.py', '_sign_normalised', 'tests/gates/numerics/_dense_eigen_oracle.py', '_sign_normalised'),  # row 131; free function _sign_normalised(v)
)

KEPT: tuple[Kept, ...] = (  # Kept(file, qualname, keeper): the keeper is the wave or the consumer that keeps the row; asserted still defined
    Kept('orpheus/numerics/basis/base.py', 'Basis.synthesize', "W8a, concept row 'Basis' (one table contraction)"),  # row 38
    Kept('orpheus/numerics/coupled_system.py', 'CoupledSpace.system_slices', "W1, concept row 'Direct sum of spaces' (4 offset loops become one)"),  # row 47
    Kept('orpheus/numerics/coupled_system.py', 'CoupledSpace.kernel_basis', "W1/W2: the metric-kernel family (concept 'Quotient by a metric kernel', A3)"),  # row 48
    Kept('orpheus/numerics/field.py', 'Field.cone_violations', "concept row 'Positive cone' (one, landed); W7 'Wall law on its cone'"),  # row 53
    Kept('orpheus/numerics/metric.py', 'HilbertMetric.kernel_basis', 'as row 48'),  # row 62
    Kept('orpheus/numerics/posing.py', 'ALPHA_MAP', "W5, concept row 'Question bound to a system' (spectral map derived from the parameter)"),  # row 71
    Kept('orpheus/numerics/space.py', 'FunctionSpace.kernel_basis', 'as row 48'),  # row 75
    Kept('orpheus/numerics/field.py', 'Field.where_largest', 'deferred, private state (decision D2; section 8.3)'),  # row 52
)

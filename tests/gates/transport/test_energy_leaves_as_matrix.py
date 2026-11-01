"""MOVED from tests/gates/numerics/test_matrix_inverse_operator.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
from orpheus.transport.fields.cross_section_field import CrossSectionField
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from orpheus.numerics.assembled_operator import dense
from dataclasses import replace


pytestmark = pytest.mark.foundation




_EPS = np.finfo(float).eps
















def _asymmetric_2g_carrier():
    """2G mixture with ASYMMETRIC ``SigS`` AND non-zero asymmetric ``Sig2``
    on a unit one-cell carrier (Mode-6: a symmetric matrix is
    blind to the transpose-assembly mutation; zero ``Sig2`` would vacuum
    the ``IsotropicN2N`` oracle row)."""
    from orpheus.data.macro_xs.synthetic import make_mixture
    from tests.gates.transport._carrier_helpers import unit_cell_carrier

    sig_s0 = np.array([[0.30, 0.08], [0.01, 0.45]])  # [g_from, g_to], asymmetric
    sig_2 = np.array([[0.020, 0.006], [0.001, 0.030]])  # asymmetric, non-zero
    m = make_mixture(
        sig_t=np.array([1.0, 1.5]), sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=sig_s0,
    )
    m = replace(m, SigS=[csr_matrix(sig_s0)], Sig2=[csr_matrix(sig_2)])
    return unit_cell_carrier({0: m}), sig_s0, sig_2, np.array([1.0, 1.5])






def test_as_matrix_energy_leaves_vs_storage_oracle():
    """§28.2 ASM-ORACLE — structural independence, spelled out.

    ``as_matrix(basis_shape=(ng,1))`` drives ``apply`` → the
    ``add_p0_source`` einsum kernel of either channel (per-column matvec
    accumulation over g'); ``dense_per_material()[mid]`` reads
    ``sig_s_legendre(mid)[0].T`` — a direct storage TRANSPOSE-copy. Neither
    computes the other: they agree ONLY if the apply is faithful to the
    stored cross sections — exactly the L0 claim. The mixture is ≥2G with
    ASYMMETRIC SigS/Sig2 (Mode-6: symmetric would null the transpose
    mutation).
    """
    from orpheus.transport.operators.isotropic_transfer import (
        IsotropicN2N,
        IsotropicScattering,
    )
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    hub, sig_s0, sig_2, sig_t = _asymmetric_2g_carrier()
    mat_xs = hub.mat_xs
    ng = 2
    for op in (
        IsotropicScattering.from_material_xs(
            mat_xs, space=mat_xs.space,
        ),
        IsotropicN2N.from_material_xs(
            mat_xs, space=mat_xs.space,
        ),
    ):
        got = op.as_matrix(basis_shape=(ng, 1))
        ref = op.dense_per_material()[0]  # the single material
        np.testing.assert_allclose(
            got, ref, rtol=0, atol=4 * _EPS * np.abs(ref).max(),
            err_msg=f"{type(op).__name__}: apply-to-basis ≠ storage transpose",
        )
    # The OperatorSum path — C − K_iso materialises as the fused storage oracle.
    loss = MultiplicationOperator.from_mesh(
        CrossSectionField.total_of(mat_xs), hub,
    ) - (
        IsotropicScattering.from_material_xs(
            mat_xs, space=mat_xs.space,
        )
        + IsotropicN2N.from_material_xs(
            mat_xs, space=mat_xs.space,
        )
    )
    ref_sum = np.diag(sig_t) - (sig_s0 + 2.0 * sig_2).T
    np.testing.assert_allclose(
        loss.as_matrix(basis_shape=(ng, 1)), ref_sum, atol=1e-12,
        err_msg="OperatorSum as_matrix ≠ diag(Σt) − (Σs0 + 2Σ2)ᵀ",
    )


def test_as_matrix_equals_retired_as_dense_loop():
    """HOMOG-EQUIV (#226 step 5, spec §32.2): the promoted ``as_matrix``
    reproduces the retired ``_as_dense`` apply-to-basis loop BYTE-for-byte
    on the homogeneous solver's own loss composition — the retirement is a
    pure relocation, not a re-derivation (the fuller-view-oracle
    discipline: the retired loop lives on here as the test-local
    reference, guarding a later ``as_matrix`` refactor from silently
    moving the homogeneous basis columns). Keep through the merge cycle,
    then reassess.

    DEMOTED at numerics_untangling S3 (ruling 60), stated: ``as_matrix``
    IS the apply-to-basis loop again, so its leg compares one algorithm
    with itself and pins only the column-order convention. The ``dense``
    leg (the structural assembly the solver factors) is the second
    algorithm that keeps the row a comparison (`[M]` bit-equal at S3).

    HOME delta vs spec §35.1: hosted HERE (foundation-marked) rather than
    in ``test_homogeneous.py``, whose file-level ``l1 + verifies(...)``
    marks would have written FALSE equation-TESTS edges for this pure
    software invariant (the registry warned on the marker conflict).
    """
    from orpheus.derivations import get
    from tests.gates.transport._carrier_helpers import unit_cell_carrier
    from orpheus.transport.operators.isotropic_transfer import (
        IsotropicN2N,
        IsotropicScattering,
    )
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    case = get("homo_2eg_n2n")  # asymmetric SigS + non-zero Sig2 (Mode-6)
    mix = next(iter(case.materials.values()))
    hub = unit_cell_carrier({0: mix})
    mat_xs = hub.mat_xs
    ng = mix.ng
    loss = MultiplicationOperator.from_mesh(
        CrossSectionField.total_of(mat_xs), hub,
    ) - (
        IsotropicScattering.from_material_xs(
            mat_xs, space=mat_xs.space,
        )
        + IsotropicN2N.from_material_xs(
            mat_xs, space=mat_xs.space,
        )
    )
    # local oracle = the retired _as_dense loop, verbatim
    cols = []
    for i in range(ng):
        e_i = np.zeros((ng, 1))
        e_i[i, 0] = 1.0
        cols.append(np.asarray(loss.apply(e_i)).ravel())
    oracle = np.column_stack(cols)
    np.testing.assert_array_equal(
        loss.as_matrix(basis_shape=(ng, 1)), oracle,
        err_msg="as_matrix diverged from the retired _as_dense loop",
    )
    np.testing.assert_array_equal(
        dense(loss, basis_shape=(ng, 1)), oracle,
        err_msg="dense (the assembly) diverged from the apply-to-basis loop",
    )



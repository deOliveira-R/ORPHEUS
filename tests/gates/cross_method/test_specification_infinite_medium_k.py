"""The infinite-medium k of a specification is the dense pencil's k, through the 0-D and the SN lift (#405 P1 step 8, S8.4 (ii)).

Specified by the test-architect (2026-10-02, ``.claude/plans/reference_p1_spec.md`` §1.8, S8.4) as one of the two rows of the
specification's assembly rung, an integration rung green on arrival (no first red of its own): it shows that a specification
carries what a method needs, through the test's own lift, since posing a method from a specification is P4's. Its teeth are the
reference side's activation legs and the battery arms named in the spec. Moved from
``tests/gates/specification/test_specification_assembly.py`` by the layer-DAG re-filing (Phase 1b, 2026-10-26), body unchanged:
one claim (``InfiniteMediumSpecification`` asks k-infinity of its one material) read through two methods, the homogeneous
baseline and an SN reflective slab, each against one reference independent of every operator, so the row lives with the
cross-method comparisons (``tests/gates/cross_method/``, which may import every method). The other row of the rung is
``tests/gates/transport/test_specification_material_mesh.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.geometry import MirrorDeck, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.question import Eigen
from orpheus.specification import InfiniteMediumSpecification
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import fuel

pytestmark = pytest.mark.l1

F = Channel.FISSION_EMISSION



def _k() -> CellCoefficient:
    return CellCoefficient.every(F)



def _dense_k_inf(mixture) -> float:
    """k-infinity from the raw ``Mixture`` arrays, independent of every operator:
    A = diag(Sigma_t) - Sigma_s0^T - 2 Sigma_2^T, F = chi (x) nu Sigma_f, k = rho(A^-1 F)."""
    a = np.diag(mixture.SigT) - mixture.SigS[0].toarray().T - 2.0 * mixture.Sig2[0].toarray().T
    f = np.outer(mixture.chi, mixture.SigP)
    return float(max(abs(np.linalg.eigvals(np.linalg.solve(a, f)))))



@pytest.mark.rests_on(
    "tests/gates/sn/solve/test_d3_admission.py::test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective",
    "tests/gates/derivations/test_adjoint_spectrum_reference.py::TestAdjointSpectrumReference",
)
@pytest.mark.parametrize("lift", ["homogeneous", "sn-reflective-slab"])
def test_s8_4_the_infinite_medium_k_is_the_dense_pencil_k(lift: str) -> None:
    """``InfiniteMediumSpecification(0, fuel, Eigen(every(F)))`` asks k-infinity of its one material.

    The fixture activates what a 2-group infinite medium can: upscatter and a
    non-zero (n,2n) emission (asserted, with the reference's sensitivity to the
    multiplicity 2). Tolerances: the 0-D solve is direct (``[M]`` 1 ULP from the
    dense value; 50 ULP admitted, never a LAPACK bit pin, vv anti-#38); the SN
    solve is iterative, 10 x its ``keff_tol`` (``[M]`` 1.5e-13 at 1e-10)."""
    from orpheus.homogeneous.solver import solve_homogeneous_infinite
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    spec = InfiniteMediumSpecification(0, fuel(), Eigen(_k()))
    require(isinstance(spec.question, Eigen) and spec.question.parameter == CellCoefficient({(0, F)}), "the rung does not ask k")
    mixture = spec.mixture
    require(mixture.Sig2[0].count_nonzero() > 0 and mixture.SigS[0].toarray()[1, 0] > 0, "activation: (n,2n), upscatter")
    k_ref = _dense_k_inf(mixture)
    a1 = np.diag(mixture.SigT) - mixture.SigS[0].toarray().T - mixture.Sig2[0].toarray().T
    k_mult1 = float(max(abs(np.linalg.eigvals(np.linalg.solve(a1, np.outer(mixture.chi, mixture.SigP))))))
    require(abs(k_ref - k_mult1) > 1e-2, "activation: the reference does not read the (n,2n) multiplicity")
    if lift == "homogeneous":
        k = solve_homogeneous_infinite(mixture).k_inf
        tol = 50 * np.spacing(k_ref)
    else:
        keff_tol = 1e-10
        mesh = Mesher(StructuredGeometry.from_homogeneous(0.5, MirrorDeck())).partition(CellsByCount.uniform_width(2)).mesh
        k = solve_sn(dict(spec.materials.items()), mesh, Quadrature.gauss_legendre(n_ordinates=8), keff_tol=keff_tol, inner_tol=1e-11).outcome.keff
        tol = 10 * keff_tol
    require(abs(k - k_ref) <= tol, f"{lift}: k {k!r} vs the dense pencil {k_ref!r} (|diff| {abs(k - k_ref):.3e} > {tol:.1e})")

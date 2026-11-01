"""MOVED from tests/gates/numerics/test_assembly_functor.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
from orpheus.transport.fields.cross_section_field import CrossSectionField
import pathlib
import numpy as np
import pytest
import orpheus.numerics.assembled_operator as asm
import orpheus.numerics.matrix_inverse_operator as mio
from orpheus.numerics.assembled_operator import assemble, dense
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import LinearOperator


pytestmark = [pytest.mark.foundation]


ORPHEUS = pathlib.Path(asm.__file__).resolve().parents[1]


SUM_LAW = "[Σ cᵢwᵢ] = Σ cᵢ[wᵢ]"


WORD_LAW = "[A₁…Aₙ] = [A₁]…[Aₙ]"


















def _lu_inputs(monkeypatch) -> list[np.ndarray]:
    seen: list[np.ndarray] = []
    original = mio.lu_factor

    def spy(a, *args, **kwargs):
        seen.append(np.array(a, copy=True))
        return original(a, *args, **kwargs)

    monkeypatch.setattr(mio, "lu_factor", spy)
    return seen


def _functor_calls(monkeypatch) -> list[LinearOperator]:
    seen: list[LinearOperator] = []
    original = asm.assemble

    def spy(op):
        seen.append(op)
        return original(op)

    monkeypatch.setattr(asm, "assemble", spy)
    return seen






class TestMatrixInverseFactorsDense:
    def test_g8a_the_factored_matrix_is_dense_on_a_route_discriminating_loss(
        self, monkeypatch,
    ) -> None:
        """The diffusion loss with its collision from Σ_a: a legitimate
        operator whose probe and assembly differ in the last bit (`[M]` S3
        spec finding 3: entries (2, 2), (6, 6)), so this VALUE row sees
        which route ``MatrixInverseOperator`` took."""
        from tests.gates.diffusion.test_operators import _config_setup
        from orpheus.diffusion.operators import DiffusionBoundaryOperator, LeakageOperator
        from orpheus.transport.operators.isotropic_transfer import IsotropicScattering
        from orpheus.transport.operators.lift import BulkLift
        from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

        mesh, mat_xs, template, _ = _config_setup("zero_flux/zero_flux")
        ffs = mesh.full_field_space
        loss = (
            LeakageOperator(mesh)
            + MultiplicationOperator(CrossSectionField.absorption_of(mat_xs), domain=ffs, codomain=ffs)
            - BulkLift(IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space), domain=ffs, codomain=ffs)
            - DiffusionBoundaryOperator(mesh)
        )
        flat = FlattenedOperator(loss, template)
        assembled = dense(flat)
        assert not np.array_equal(flat.as_matrix(), assembled), (
            "premise lost: probe and assembly agree bitwise, the row cannot see the route"
        )
        seen = _lu_inputs(monkeypatch)
        MatrixInverseOperator(flat)
        assert len(seen) == 1
        np.testing.assert_array_equal(seen[0], assembled)

    def test_g8b_the_diffusion_resolvent_assembles(self, monkeypatch) -> None:
        from tests.gates.diffusion.test_operators import _config_setup, _loss

        mesh, mat_xs, template, _ = _config_setup("reflective/marshak_vacuum")
        flat = FlattenedOperator(_loss(mesh, mat_xs), template)
        calls = _functor_calls(monkeypatch)
        MatrixInverseOperator(flat)
        assert [op for op in calls if op is flat], calls



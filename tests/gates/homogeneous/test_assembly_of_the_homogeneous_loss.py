"""MOVED from tests/gates/numerics/test_assembly_functor.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import pathlib
import numpy as np
import pytest
import orpheus.numerics.assembled_operator as asm
import orpheus.numerics.matrix_inverse_operator as mio
from orpheus.numerics.assembled_operator import assemble, dense
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






class TestTwoAlgorithms:
    def test_g7b_probe_and_dense_agree_on_every_homogeneous_loss(self) -> None:
        """`[M]` 8 of 8 mixtures at S3: the k∞ bound reads ``dense`` (the
        matrix the solver factors); this row keeps the probe and the
        assembly one matrix on the same population."""
        from tests.gates.homogeneous._homogeneous_population import mixture_cases
        from orpheus.homogeneous.solver import HomogeneousProblem

        cases = mixture_cases()
        assert len(cases) == 8, sorted(cases)
        for name, mix in cases.items():
            loss = HomogeneousProblem(mix).loss
            assert loss.is_assemblable, name
            np.testing.assert_array_equal(loss.as_matrix(), dense(loss), err_msg=name)


class TestMatrixInverseFactorsDense:
    def test_g8b_the_homogeneous_solver_assembles(self, monkeypatch) -> None:
        from orpheus.derivations import get
        from orpheus.homogeneous.solver import solve_homogeneous_infinite

        calls = _functor_calls(monkeypatch)
        seen = _lu_inputs(monkeypatch)
        solve_homogeneous_infinite(next(iter(get("homo_2eg_n2n").materials.values())))
        assert calls, "the homogeneous solver never reached the functor"
        assert len(seen) >= 1
        reached = list(calls)
        assert any(np.array_equal(s, dense(op)) for s in seen for op in reached), (
            "no factored matrix is the assembly of an operator the solver assembled"
        )



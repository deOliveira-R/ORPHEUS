r"""No S\ :sub:`N` entry takes a singular boundary condition, and the default filler is gone (S3.10, P1 step 3c).

The user's ruling of 2026-09-29: ``None`` retires as a boundary declaration everywhere; the declared conditions travel
as the keyword ``boundary_conditions``, keyed by tag. These two rows are the retirement's witness on the S\ :sub:`N`
entries (``retirement-audit`` D.18): a ``boundary_condition`` parameter re-added to any entry, or a ``**kwargs`` that
would swallow it in silence, reds the signature row; ``_apply_default_bcs`` re-added reds the name row.

They lived in ``tests/gates/mesh/test_mesh1d.py::TestTheDefaultsAreGone`` until the layer-DAG re-filing
(``.claude/plans/layer_dag.md``, 2026-10-26, candidate 4b of the Phase 3 sub-cone report): their claim is about
``orpheus.sn.solver``, which a mesh gate may not import. The class's third row (``realize_boundary`` constructs no
declaration by default) is a transport claim and stays there. Every test is ``foundation``; claim kind THEOREM.
"""
from __future__ import annotations

import pytest

pytestmark = pytest.mark.foundation

_SN_ENTRIES = (
    "solve_sn_fixed_source",
    "solve_sn_adjoint_fixed_source",
    "solve_sn_multiplying_source",
    "solve_sn",
    "solve_sn_adjoint",
    "_as_problem",
)


class TestTheDefaultsAreGone:
    r"""S3.10, step 3c, re-posed at wave 2 S2b: no S\ :sub:`N` entry takes a singular ``boundary_condition``, and
    ``_apply_default_bcs`` is gone."""

    @pytest.mark.parametrize("name", _SN_ENTRIES)
    def test_no_entry_takes_a_boundary_condition(self, name):
        import inspect

        import orpheus.sn.solver as solver

        parameters = inspect.signature(getattr(solver, name)).parameters
        assert "boundary_condition" not in parameters, name
        # a ``**kwargs`` would swallow the retired keyword in silence
        kinds = {p.kind for p in parameters.values()}
        assert inspect.Parameter.VAR_KEYWORD not in kinds, name

    def test_the_default_filler_is_gone(self):
        import orpheus.sn.solver as solver

        assert not hasattr(solver, "_apply_default_bcs")

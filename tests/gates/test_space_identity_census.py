"""MOVED from tests/gates/numerics/test_axis_identity.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import pytest
from orpheus.numerics import frame_faces as ff
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import DualSpace, FunctionSpace
from orpheus.numerics.spaces.full_field_space import FullFieldSpace


pytestmark = pytest.mark.foundation












def _population():
    from tests.gates.sn.regression import _generate_snapshots as G
    from orpheus.sn.solver import _as_problem
    out = []
    for c in (G._slab_3region("2g", 6), G._cylinder_3region("2g", 4, "folded_4x6")):
        p = _as_problem(c["mesh"], c["quadrature"], c["materials"], scattering_order=1, boundary_conditions=G.case_conditions(c))
        ffs = p.full_field_space
        out += [ffs, ffs.dual(), ffs.interior_space, ffs.interior_space.dual(), getattr(ffs, "trace_space", None), p.system.space,
                p.moment_space(1), p.angular_bulk_space]
        rc = getattr(p, "radial_characteristic_field_space", None)
        if rc is not None:
            out += [rc, rc.interior_space]
    q = Quadrature.folded_product(n_mu=2, n_phi=4)
    out += [ff.basis_space(ff.angular_frame(q, 2)), FunctionSpace("named", (3,))]
    return [s for s in out if isinstance(s, FunctionSpace)]


def test_K2_eq_and_hash_derive_from_the_identity_key():
    pop = _population()
    kinds = {type(s).__name__ for s in pop}
    assert len(pop) >= 12 and {"FullFieldSpace", "CoupledSpace", "TensorProductSpace", "DualSpace"} <= kinds, kinds
    for s in pop:
        assert hash(s) == hash(s.identity_key()), f"{type(s).__name__}: hash is not hash(identity_key)"
    for a in pop:
        for b in pop:
            assert (a == b) == (a.identity_key() == b.identity_key()), f"{type(a).__name__} vs {type(b).__name__}"



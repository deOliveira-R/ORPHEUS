r"""S6a on the S\ :sub:`N` carriers — the composite duals and the mint count.

The production carriers are a ``CoupledSpace`` over ``FullFieldSpace`` members
(slab, sphere, cylinder, and a cylinder whose trace metric is exactly zero on
its tangential ordinates). Their dual is the coupled space of the member duals
(biproduct duality, rulings 26 and 75); the dagger–dual square commutes bit for
bit on the Strategy's implicit operator; and the DualSpace constructions of an
adjoint solve do not grow with its iteration count. Row ids are those of the
verification spec (``scratch/numerics_untangling/packages/S6a_gate_spec.md``).
The numerics rows are ``tests/gates/numerics/test_dual_metric.py``.

``-O``-safe: every assertion is ``_require`` (a raise) or ``np.testing``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import DualSpace

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_dual_metric.py::test_S1_the_dual_space_reads_the_inverse_metric_bitwise",
        "tests/gates/numerics/test_dual_metric.py::test_K1_the_dual_is_memoised_per_space",
    ),
]

_EPS = float(np.finfo(float).eps)


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


# ═══════════════════════════════════════════════════════════════════════
# C — the composite duals (SN fixtures)
# ═══════════════════════════════════════════════════════════════════════


def _coupled_domain(M):
    from orpheus.numerics.coupled_system import CoupledSpace

    V = M.domain
    if not isinstance(V, CoupledSpace):
        pytest.fail(f"fixture: the implicit operator's domain is a {type(V).__name__}, not a CoupledSpace")
    return V


def _sn_problem(name):
    from orpheus.sn.solver import _adjoint_posing_parts, _as_problem
    from tests.gates.sn.regression import _generate_snapshots as G

    c = {
        "slab3": lambda: G._slab_3region("2g", 40),
        "sphere3": lambda: G._sphere_3region("2g", 40),
        "cyl3": lambda: G._cylinder_3region("2g", 20, "folded_4x6"),
        "cyl_degenerate": lambda: dict(G._cylinder_3region("2g", 20, "folded_4x6"),
                                       quadrature=Quadrature.folded_product(n_mu=2, n_phi=6)),
    }[name]()
    problem = _as_problem(c["mesh"], c["quadrature"], c["materials"], scattering_order=c["scattering_order"])
    M, N, F, template, _ = _adjoint_posing_parts(problem)
    return problem, M, template


def _coupled_sample(template, seed):
    from orpheus.numerics.coupled_system import CoupledField

    n = np.asarray(template.to_flat()).size
    return CoupledField.from_flat(np.random.default_rng(seed).standard_normal(n), template)


_SN = ["slab3", "sphere3", "cyl3", "cyl_degenerate"]


@pytest.mark.parametrize("name", _SN)
def test_C1_a_composite_dual_reads_the_member_inverse_metrics(name):
    """On the production carriers (``CoupledSpace`` over ``FullFieldSpace``
    members; the cylinder's trace degenerate): ``V*.apply_metric ≡
    V.apply_inverse_metric`` and the reverse, bitwise, and the dual metric MOVES
    the sample. First red TODAY (census §3(a)): ``V*`` is a ``DualSpace`` with no
    members, Euclidean (``[M]`` probe p2: the dual metric moves y by 0.0 while the
    primal's moves it by 3.3)."""
    _, M, template = _sn_problem(name)
    V = _coupled_domain(M)
    y = _coupled_sample(template, 31)
    Vd = V.dual()
    np.testing.assert_array_equal(np.asarray(Vd.apply_metric(y).to_flat()), np.asarray(V.apply_inverse_metric(y).to_flat()))
    np.testing.assert_array_equal(np.asarray(Vd.apply_inverse_metric(y).to_flat()), np.asarray(V.apply_metric(y).to_flat()))
    _require(not np.array_equal(np.asarray(Vd.apply_metric(y).to_flat()), np.asarray(y.to_flat())),
             f"[{name}] the composite dual metric is the identity")
    f, g = _coupled_sample(template, 32), _coupled_sample(template, 33)
    lhs = Vd.inner_product(f, g)
    rhs = float(np.asarray(V.apply_inverse_metric(f).to_flat()) @ np.asarray(g.to_flat()))
    _require(abs(lhs - rhs) <= 1e3 * _EPS * max(1.0, abs(rhs)), f"[{name}] ⟨f,g⟩_V* {lhs!r} vs {rhs!r}")


@pytest.mark.parametrize("name", _SN)
def test_C2_the_composite_dual_is_the_biproduct_of_member_duals(name):
    """``CoupledSpace.dual()`` is a ``CoupledSpace`` whose systems are the member
    duals (each a ``FullFieldSpace`` over the interior and trace duals); it is not
    a ``DualSpace``; ``V*.dual() is V`` through the memo; the zero factories swap
    (the cotangent zeros are V*'s zeros; accepted 2026-10-21); ``V* ≠ V``."""
    from orpheus.numerics.coupled_system import CoupledSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    _, M, _ = _sn_problem(name)
    V = _coupled_domain(M)
    Vd = V.dual()
    if not isinstance(Vd, CoupledSpace) or type(Vd) is not CoupledSpace or isinstance(Vd, DualSpace):
        pytest.fail(f"V* is a {type(Vd).__name__}, not the CoupledSpace of member duals")
    for s, sd in zip(V.systems, Vd.systems, strict=True):
        _require(sd == s.dual(), f"member dual mismatch on {s.name}")
        if isinstance(s, FullFieldSpace):
            si, st = s.interior_space, s.trace_space
            if not isinstance(sd, FullFieldSpace):
                pytest.fail(f"the full-field dual is a {type(sd).__name__}")
            _require(sd.interior_space == si.dual() and sd.trace_space == st.dual(),
                     "the full-field dual is not interior* ⊕ trace*")
    _require(Vd.dual() is V, "V*.dual() is not V")
    _require(Vd.zeros_factory is V.dual_zeros_factory and Vd.dual_zeros_factory is V.zeros_factory,
             "the zero factories did not swap")
    _require(Vd != V and V != Vd, "V* == V on a composite")


@pytest.mark.parametrize("name", _SN)
def test_C3_the_dagger_dual_square_commutes_on_the_sn_carrier(name):
    """``M.dual().H ≡ M.H.dual()`` bitwise on the Strategy's implicit operator
    (the baseline's ``ops_*`` rows). After the carve the left side runs the BASE
    sandwich over the composite dual metrics, so this is the end-to-end of C1.
    `[M]` arm ``euclid_dual`` (base sandwich, composite duals Euclidean): 0.43
    (slab), 0.94 (sphere), 0.76 (cyl), 0.74 (degenerate cyl) relative."""
    _, M, template = _sn_problem(name)
    y = _coupled_sample(template, 41)
    np.testing.assert_array_equal(np.asarray(M.dual().H.apply(y).to_flat()), np.asarray(M.H.dual().apply(y).to_flat()))


# ═══════════════════════════════════════════════════════════════════════
# K2 — the mint count of a solve does not grow with its iterations
# ═══════════════════════════════════════════════════════════════════════


def test_K2_dual_space_mints_do_not_grow_with_the_iteration_count(monkeypatch):
    """``solve_sn_adjoint`` on the 3-region slab at two tolerance settings (3 vs 5
    outers): the DualSpace constructions are EQUAL and bounded by the number of
    distinct primals asked for. `[M]` TODAY 5645 vs 10965, every one from
    ``_DualOperator.domain``/``.codomain`` (probe p1); a per-arrow cache leaves both
    unchanged (the _DualOperators are rebuilt per apply, 3383 of 3394 at
    ``OperatorProduct.dual``; probe p3). Activation: the tight solve runs more
    outers (asserted)."""
    from orpheus.sn.solver import solve_sn_adjoint
    from tests.gates.sn.regression import _generate_snapshots as G

    real_init = DualSpace.__init__
    mints, primals = [], []

    def counting(self, *a, **k):
        real_init(self, *a, **k)
        mints.append(1)
        primals.append(self.primal)

    monkeypatch.setattr(DualSpace, "__init__", counting)
    counts, outers = [], []
    for tols in ((1e-7, 1e-6, 1e-8), (1e-10, 1e-9, 1e-11)):
        c = G._slab_3region("2g", 40)
        mints.clear(); primals.clear()
        sol = solve_sn_adjoint(c["materials"], c["mesh"], c["quadrature"],
                               keff_tol=tols[0], flux_tol=tols[1], inner_tol=tols[2])
        counts.append(len(mints))
        outers.append(len(sol.outcome.trajectory))
        distinct = len({id(p) for p in primals})
        _require(len(mints) <= max(distinct, 1) * 1, f"{len(mints)} DualSpace mints for {distinct} distinct primals")
    _require(outers[1] > outers[0], f"activation: the tight solve ran {outers} outers")
    _require(counts[0] == counts[1], f"DualSpace mints grow with iterations: {counts} for outers {outers}")




def test_R2b_dual_space_of_a_composite_is_its_composite_dual():
    """``DualSpace.of(V)`` on the coupled carrier and on its full-field member
    is the composite dual (qa F2), never a ``DualSpace`` over the composite.
    First red: ``of`` wrapping the composite (the pre-S6a member-less dual)."""
    from orpheus.numerics.coupled_system import CoupledSpace
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace

    _, M, _ = _sn_problem("slab3")
    V = _coupled_domain(M)
    _require(DualSpace.of(V) is V.dual() and type(V.dual()) is CoupledSpace, "of(coupled) is not its dual")
    member = V.systems[0]
    _require(isinstance(member, FullFieldSpace) and DualSpace.of(member) is member.dual()
             and type(member.dual()) is FullFieldSpace, "of(full field) is not its dual")


@pytest.mark.parametrize("name", ["slab3", "cyl_degenerate"])
def test_R5_the_dual_of_a_sum_has_its_adjoint_by_reciprocity(name):
    """qa's row: ``S = M + N`` (the Strategy's implicit part plus its gain),
    ``⟨S.dual() f, g⟩_{V*} = ⟨f, S.dual().H g⟩_{W*}`` in the dual end metrics,
    which are seminorms on the degenerate cylinder (the pairing never reads a
    kernel component, so the comparison is in the G-seminorm, never
    componentwise), 1e-12 relative; and the square ``S.dual().H`` vs
    ``S.H.dual()`` on the same pairing. First red: a composite dual read as
    Euclidean (the census §3(a) state), arm ``composite_euclid``."""
    from orpheus.sn.solver import _adjoint_posing_parts

    problem, _, _ = _sn_problem(name)
    M, N, _F, template, _ = _adjoint_posing_parts(problem)
    S = M + N
    D = S.dual()  # W* → V*
    # Carrier classes, [M] 2026-10-21: D maps a flux-classed composite to a
    # source-classed one, and D.H (like S.H.dual()) also takes a flux-classed
    # one, so h is a flux sample. The pairings below read values only.
    f = _coupled_sample(template, 81)
    h = _coupled_sample(template, 82)
    lhs = D.codomain.inner_product(D.apply(f), h)
    a = D.H.apply(h)
    rhs = D.domain.inner_product(f, a)
    _require(abs(lhs - rhs) <= 1e-12 * max(1.0, abs(lhs), abs(rhs)),
             f"[{name}] ⟨Dᵀf, h⟩ = {lhs!r} vs ⟨f, Dᴴh⟩ = {rhs!r}")
    b = S.H.dual().apply(h)
    gap = abs(D.domain.inner_product(f, a) - D.domain.inner_product(f, b))
    _require(gap <= 1e-12 * max(1.0, abs(rhs)), f"[{name}] (S.dual()).H and (S.H).dual() differ by {gap:.2e} in the seminorm")

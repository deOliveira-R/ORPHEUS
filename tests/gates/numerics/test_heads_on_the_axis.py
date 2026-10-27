r"""S7b gates: the moment heads live on the head axis, and every basis spans the one generic coefficient axis.

Spec: ``scratch/numerics_untangling/packages/S7_gate_spec.md`` §3 rows H*, written by the test-architect
against the names of ``S7_design_S7b.md``; rulings 72, 81, 83, 84. Names:
:mod:`orpheus.numerics.coefficient_layout` (``HarmonicLayout``, ``LegendreLayout``, ``CellLayout``),
``Basis.layout``, ``Axis.layout``, ``FromBasis(basis)``, ``Axis.of_basis`` / ``FunctionSpace.of_basis``,
``moment_head_axis(space).layout``. H2b's pair (equal widths, shifted edges) read EQUAL raw and dressed on
the pre-carve and S7a trees (``S7_gates/probes/p4_h2b_shifted_edges.log``); the main agent ruled that the
cell layout carries the edges.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from orpheus.numerics import frame_faces as ff
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation


def _sh(L):
    from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
    return SphericalHarmonicBasis(L=L)


def _mirror(L):
    from orpheus.numerics.basis.spherical_harmonic_basis import MirrorEvenSphericalHarmonicBasis
    return MirrorEvenSphericalHarmonicBasis(L=L)


def _legendre(L, axis):
    from orpheus.numerics.basis.legendre_basis import LegendreBasis
    return LegendreBasis(L=L, axis=axis)


def _indicator(edges):
    from orpheus.numerics.basis.indicator_basis import IndicatorBasis
    from orpheus.geometry.manifold import RealSpace
    return IndicatorBasis((np.asarray(edges, float),), RealSpace(1))


def _cell_measure(edges, per_cell=2):
    """``per_cell`` equal-weight midpoint atoms per cell: the measure an indicator frame pairs with."""
    from orpheus.geometry.manifold import RealSpace
    from orpheus.numerics.measure import DiscreteMeasure
    e = np.asarray(edges, float)
    nodes, weights = [], []
    for lo, hi in zip(e[:-1], e[1:]):
        h = (hi - lo) / per_cell
        nodes += [lo + (k + 0.5) * h for k in range(per_cell)]
        weights += [h] * per_cell
    return DiscreteMeasure(nodes=np.array(nodes), weights=np.array(weights), support=RealSpace(1))


def _separated(a, b, what):
    assert not (a == b) and a != b and len({a, b}) == 2, f"{what}: compare EQUAL"
    assert hash(a) != hash(b), f"{what}: equal hashes"


def _one(a, b, what):
    assert a == b and hash(a) == hash(b) and len({a, b}) == 1, f"{what}: two mints of one head differ"


# ── the layout value (tier 0) ─────────────────────────────────────────────

@pytest.mark.parametrize("L", [0, 1, 2, 3])
def test_H0_the_layout_values_match_the_coefficient_tables(L):
    """Hand reference (the table's definition, not the code): rectangular SH (L+1, 2L+1), slot (0,0),
    block (l, 0:2l+1); flat Legendre (L+1,), slot (0,), block (l,); unit e0 at the slot; a block above L
    refuses; a cell layout over three cells has shape (3,); a cell layout refuses non-increasing edges."""
    from orpheus.numerics.coefficient_layout import CellLayout, HarmonicLayout, LegendreLayout
    h, g = HarmonicLayout(L), LegendreLayout(L)
    assert h.L == L and h.shape == (L + 1, 2 * L + 1) and h.isotropic_slot == (0, 0)
    assert g.L == L and g.shape == (L + 1,) and g.isotropic_slot == (0,)
    for l in range(L + 1):
        assert h.degree_block(l) == (l, slice(0, 2 * l + 1)) and g.degree_block(l) == (l,)
    for lay in (h, g):
        e0 = np.zeros(lay.shape)
        e0[lay.isotropic_slot] = 1.0
        np.testing.assert_array_equal(lay.unit, e0)
        with pytest.raises(ValueError):
            lay.degree_block(L + 1)
    assert HarmonicLayout(L) == HarmonicLayout(L) and HarmonicLayout(1) != LegendreLayout(1)
    c = CellLayout((np.array([0.0, 1.0, 2.0, 3.0]),))
    assert c.shape == (3,)
    for bad in ([2.0, 1.0, 1.0, 0.0], [0.0, np.nan, 1.0]):
        with pytest.raises(ValueError, match="strictly increasing"):
            CellLayout((np.array(bad),))


def test_H0b_each_basis_carries_its_layout():
    from orpheus.numerics.coefficient_layout import CellLayout, HarmonicLayout, LegendreLayout
    assert _sh(2).layout == HarmonicLayout(2) and _mirror(2).layout == HarmonicLayout(2)
    assert _legendre(2, "x").layout == LegendreLayout(2)
    assert _indicator([0.0, 1.0, 2.0]).layout == CellLayout((np.array([0.0, 1.0, 2.0]),))


def test_H0c_the_origin_kind_is_one_match_on_the_layout():
    """FromBasis: CellLayout → NODAL; a degree layout → Modal(layout.unit). Axis.layout is read from the origin."""
    from orpheus.numerics.axis import Axis, BasisKind, Modal
    sh = Axis.of_basis(_sh(2))
    assert isinstance(sh.kind, Modal)
    np.testing.assert_array_equal(sh.kind.unit, _sh(2).layout.unit)
    assert sh.layout == _sh(2).layout
    cell = Axis.of_basis(_indicator([0.0, 1.0, 2.0]))
    assert cell.kind is BasisKind.NODAL and cell.layout == _indicator([0.0, 1.0, 2.0]).layout
    q = Quadrature.gauss_legendre(4)
    assert Axis.of_quadrature(q).layout is None, "a sampled axis has no coefficient layout"


# ── the #631 pairs on the basis side ─────────────────────────────────────

@pytest.mark.parametrize("L", [1, 2, 3])
def test_H1_pair3_the_basis_domain_is_part_of_the_head(L):
    """#631 pair 3. SH(L) on S² vs mirror-even SH(L) on S²/σ_y: equal layout, different support."""
    _separated(FunctionSpace.of_basis(_sh(L)), FunctionSpace.of_basis(_mirror(L)), f"SH vs mirror-even, L={L}")
    _one(FunctionSpace.of_basis(_sh(L)), FunctionSpace.of_basis(_sh(L)), "control: two SH(L) mints")


def test_H2_pair4_indicator_edges_are_content():
    """#631 pair 4. Edges [0,1,2] vs [0,0.5,2] (they differ in edges AND masses), raw and dressed."""
    from orpheus.numerics.frame import GalerkinFrame
    a, b = _indicator([0.0, 1.0, 2.0]), _indicator([0.0, 0.5, 2.0])
    _separated(FunctionSpace.of_basis(a), FunctionSpace.of_basis(b), "indicator heads")
    m = _cell_measure([0.0, 1.0, 2.0])
    da, db = ff.basis_space(GalerkinFrame(a, m)), ff.basis_space(GalerkinFrame(b, m))
    _separated(da, db, "dressed indicator heads")
    (ax,) = FunctionSpace.of_basis(a).axes
    np.testing.assert_array_equal(ax.unit, np.ones(2))       # NODAL, unit ones
    (dax,) = da.axes
    np.testing.assert_array_equal(dax.unit, np.array([1.0, 1.0]))  # ruling 83: Modal(masses), masses = cell widths
    assert da.has_coordinate_cone is False


def test_H2c_the_dressed_coarse_head_unit_is_the_cell_masses():
    """Ruling 83's witness on NON-unit cells (qa F3: H2's fixture has unit widths, so the masses equal the raw
    unit of ones and ``Modal(ones)`` passed it). Edges [0, 1, 3]: masses [1, 2], read off the frame's measure."""
    from orpheus.numerics.axis import Modal
    from orpheus.numerics.frame import GalerkinFrame
    edges = [0.0, 1.0, 3.0]
    (dax,) = ff.basis_space(GalerkinFrame(_indicator(edges), _cell_measure(edges))).axes
    assert isinstance(dax.kind, Modal)
    np.testing.assert_allclose(dax.unit, np.array([1.0, 2.0]), rtol=0, atol=1e-15)


def test_H2d_a_raw_cell_head_is_measure_free_and_refuses_to_integrate():
    """The user's ruling (2026-10-21): the raw indicator head has no fine measure until a frame binds one
    (a spherical shell is 4πr²Δr, a coarse group counts its fine groups), so its axis marginal refuses,
    naming the frame; the frame's dressed head integrates (1 over [0, 3] is 3)."""
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.axis_marginal import AxisRetractionOperator
    from orpheus.numerics.frame import GalerkinFrame
    edges = [0.0, 1.0, 3.0]
    other = Axis("energy", (2,), origin=Declared(BasisKind.NODAL))
    raw = FunctionSpace.of_axes(Axis.of_basis(_indicator(edges)), other)
    with pytest.raises(TypeError, match="MEASURE-FREE"):
        AxisRetractionOperator.of(raw, "cells")
    dressed = ff.basis_space(GalerkinFrame(_indicator(edges), _cell_measure(edges)))
    R = AxisRetractionOperator.of(FunctionSpace.of_axes(*dressed.axes, other), "cells")
    # the constant function's coordinates on the dressed head are its moments, the unit [1, 2]
    one = np.broadcast_to(dressed.axes[0].unit[:, None], (2, 2))
    np.testing.assert_allclose(R.apply(one), np.array([3.0, 3.0]), rtol=0, atol=1e-14)


def test_H2b_equal_widths_shifted_edges_are_two_heads():
    """[0,1,2] vs [1,2,3]: equal widths, so equal masses, equal cell counts, equal support RealSpace(1); only
    the edges differ. A head identity that held the edges only through the masses COLLAPSES this pair. `[M]`
    pre-carve and S7a: equal raw, equal dressed (``probes/p4_h2b_shifted_edges.log``); the S7a measure spaces
    of the two frames already separate (atoms). Green because the cell layout carries the edges (the main
    agent's ruling, 2026-10-21: the cells are the partition's atoms, ruling 71)."""
    from orpheus.numerics.frame import GalerkinFrame
    a, b = _indicator([0.0, 1.0, 2.0]), _indicator([1.0, 2.0, 3.0])
    ra, rb = FunctionSpace.of_basis(a), FunctionSpace.of_basis(b)
    np.testing.assert_array_equal(ra.apply_metric(np.ones(2)), rb.apply_metric(np.ones(2)))  # PREMISE: one metric
    _separated(ra, rb, "shifted-edge indicator heads")
    ma, mb = _cell_measure([0.0, 1.0, 2.0]), _cell_measure([1.0, 2.0, 3.0])
    _separated(ff.basis_space(GalerkinFrame(a, ma)), ff.basis_space(GalerkinFrame(b, mb)),
               "shifted-edge dressed heads")
    _one(ra, FunctionSpace.of_basis(_indicator([0.0, 1.0, 2.0])), "control: two [0,1,2] mints")


def test_H3_the_retired_classes_are_gone():
    """Two retired modules are absent (``find_spec`` is ``None``); three retired classes are absent from the modules
    that survive. The positive control: ``find_spec`` finds a module that exists."""
    from orpheus.numerics import axis
    from orpheus.numerics.spaces import moment_head

    for module_name in ("orpheus.numerics.spaces.spherical_harmonic_space", "orpheus.numerics.spaces.legendre_space"):
        assert importlib.util.find_spec(module_name) is None, f"{module_name} survives S7b (rulings 72, 81)"
    for module, name in ((moment_head, "MomentHead"), (axis, "HarmonicAxis"), (axis, "LegendreAxis")):
        assert not hasattr(module, name), f"{module.__name__}.{name} survives S7b (rulings 72, 81)"
    from orpheus.numerics.basis.base import Basis
    assert "space" not in vars(Basis), "Basis.space retires (FunctionSpace.of_basis)"
    assert importlib.util.find_spec("orpheus.numerics.coefficient_layout") is not None  # positive control


@pytest.mark.parametrize("L", [0, 1, 2, 3])
def test_H4_the_door_reads_the_layout_off_the_head_axis(L):
    """moment_head_axis(space).layout is the basis's layout, on the continuum and on the dressed head."""
    from orpheus.numerics.coefficient_layout import HarmonicLayout, LegendreLayout
    from orpheus.numerics.spaces.moment_head import moment_head_axis
    assert moment_head_axis(FunctionSpace.of_basis(_sh(L))).layout == HarmonicLayout(L)
    assert moment_head_axis(FunctionSpace.of_basis(_legendre(L, "x"))).layout == LegendreLayout(L)
    head = ff.basis_space(ff.angular_frame(Quadrature.product(4, 8), L))
    assert moment_head_axis(head).layout == HarmonicLayout(L)


def test_H5_the_door_refuses_a_non_head_with_the_shared_message():
    """A sampled angular space and a cell (indicator) head both refuse: neither carries a degree layout."""
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.spaces.moment_head import moment_head_axis
    with pytest.raises(TypeError, match="is not a moment head"):
        moment_head_axis(FunctionSpace.of_axes(Axis.of_quadrature(Quadrature.gauss_legendre(4))))
    with pytest.raises(TypeError, match="is not a moment head"):
        moment_head_axis(FunctionSpace.of_basis(_indicator([0.0, 1.0, 2.0])))


def test_H6_legendre_about_two_axes_differ_by_domain_not_by_a_spent_axis_field():
    a, b = FunctionSpace.of_basis(_legendre(1, "x")), FunctionSpace.of_basis(_legendre(1, "z"))
    _separated(a, b, "Legendre about x vs about z")
    assert a.axes is not None and not hasattr(a.axes[0], "spent_axis")
    assert a.axes[0].layout == b.axes[0].layout, "PREMISE: one layout, so the support separates them"


@pytest.mark.parametrize("rule,L", [("gl8", 2), ("product48", 2), ("folded24", 2)])
def test_H7_truncated_head_is_the_frame_head_at_the_lower_order(rule, L):
    """Green on both commits (one body); S7b re-routes it through the layout."""
    from orpheus.numerics.spaces.moment_head import truncated_head
    q = {"gl8": lambda: Quadrature.gauss_legendre(8), "product48": lambda: Quadrature.product(4, 8),
         "folded24": lambda: Quadrature.folded_product(n_mu=2, n_phi=4)}[rule]()
    head = ff.basis_space(ff.angular_frame(q, L))
    for L_new in range(L + 1):
        assert truncated_head(head, L_new) == ff.basis_space(ff.angular_frame(q, L_new))


def test_H8_the_pending_concept_edges_are_closed():
    """G0: basis/ and measure.py import no space module; coefficient_layout is tier 0."""
    from tests.gates.numerics import test_concept_order as G0
    pending = G0.PENDING
    s7 = {k: v for k, v in pending.items() if "S7" in str(v)}
    assert not s7, f"S7 PENDING edges remain: {s7}"

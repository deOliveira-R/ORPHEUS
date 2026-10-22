r"""S7b gates (DRAFT, to the level the design supports): the moment heads live on the head axis.

Spec: ``S7_gate_spec.md`` §3 rows H*. Design §4, rulings 72, 81, 83. NAMES the design fixes are used as
written (``FunctionSpace.of_basis``, ``moment_head_axis``, the head axis's ``L`` / ``isotropic_slot`` /
``degree_block``); the layout VALUE's name (ruling 81) is not fixed by the design and is not read here.
Every row is RED on today's tree (the API does not exist) and its first red is its named reason only
after S7b lands; the S7b commit re-runs them and records each first red by mutation (§5 of the spec).
"""
from __future__ import annotations

import importlib

import numpy as np
import pytest

from orpheus.numerics import frame_faces as ff
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation

# S7a lands first (ruling 84); S7b, in the same merge, makes these rows green and drops the mark.
# H7 is green on both commits (one body today, DEMOTED until S7b re-routes it), so it is unmarked.
_S7B = pytest.mark.xfail(strict=True, reason="S7b: the heads on the head axis")


def _sh(L):
    from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
    return SphericalHarmonicBasis(L=L)


def _mirror(L):
    from orpheus.numerics.basis.spherical_harmonic_basis import MirrorEvenSphericalHarmonicBasis
    return MirrorEvenSphericalHarmonicBasis(L=L)


def _legendre(L, axis):
    from orpheus.numerics.basis.legendre_basis import LegendreBasis
    return LegendreBasis(L=L, axis=axis)


def _separated(a, b, what):
    assert not (a == b) and a != b and len({a, b}) == 2, f"{what}: compare EQUAL"


@_S7B
@pytest.mark.parametrize("L", [1, 2, 3])
def test_H1_pair3_the_basis_domain_is_part_of_the_head(L):
    """#631 pair 3. SH(L) on S² and mirror-even SH(L) on S²/σ_y: TODAY ``.space`` equal at L = 1, 2, 3."""
    a, b = FunctionSpace.of_basis(_sh(L)), FunctionSpace.of_basis(_mirror(L))
    _separated(a, b, f"SH vs mirror-even SH heads, L={L}")
    c = FunctionSpace.of_basis(_sh(L))
    assert a == c and hash(a) == hash(c), "control: two mints of SH(L)"


@_S7B
def test_H2_pair4_indicator_edges_are_content():
    """#631 pair 4. Edges [0,1,2] vs [0,0.5,2]: TODAY equal (axes-less nominal), also once dressed."""
    from orpheus.numerics.basis.indicator_basis import IndicatorBasis
    from orpheus.numerics.frame import GalerkinFrame
    from orpheus.numerics.manifold import RealSpace
    from orpheus.numerics.measure import DiscreteMeasure
    a = IndicatorBasis((np.array([0.0, 1.0, 2.0]),), RealSpace(1))
    b = IndicatorBasis((np.array([0.0, 0.5, 2.0]),), RealSpace(1))
    _separated(FunctionSpace.of_basis(a), FunctionSpace.of_basis(b), "indicator heads")
    m = DiscreteMeasure(nodes=np.array([0.25, 0.75, 1.25, 1.75]), weights=np.full(4, 0.5), support=RealSpace(1))
    da, db = ff.basis_space(GalerkinFrame(a, m)), ff.basis_space(GalerkinFrame(b, m))
    _separated(da, db, "dressed indicator heads")
    (ax,) = FunctionSpace.of_basis(a).axes
    assert np.array_equal(ax.unit, np.ones(2)), "the indicator axis is NODAL, unit ones (design §4)"
    (dax,) = da.axes
    # ruling 83: the dressed coarse head is Modal(masses); masses = cell measures [1, 1] for edges [0,1,2]
    np.testing.assert_array_equal(dax.unit, np.array([1.0, 1.0]))
    assert ff.basis_space(GalerkinFrame(a, m)).has_coordinate_cone is False


@_S7B
def test_H3_the_retired_classes_are_gone():
    for mod, name in (("orpheus.numerics.spaces.spherical_harmonic_space", "SphericalHarmonicSpace"),
                      ("orpheus.numerics.spaces.legendre_space", "LegendreSpace"),
                      ("orpheus.numerics.spaces.moment_head", "MomentHead"),
                      ("orpheus.numerics.axis", "HarmonicAxis"), ("orpheus.numerics.axis", "LegendreAxis")):
        try:
            m = importlib.import_module(mod)
        except ImportError:
            continue
        assert not hasattr(m, name), f"{mod}.{name} survives S7b (rulings 72, 81)"
    # positive control: the census reaches a module that exists
    assert importlib.import_module("orpheus.numerics.axis")


@_S7B
@pytest.mark.parametrize("L", [0, 1, 2, 3])
def test_H4_the_layout_surface_on_the_head_axis_matches_the_closed_form(L):
    """Hand reference (the coefficient table's definition, not the code): rectangular SH head
    (L+1, 2L+1), isotropic slot (0, 0), degree-l block (l, 0:2l+1); flat Legendre head (L+1,),
    isotropic slot (0,), block (l,)."""
    from orpheus.numerics.spaces.moment_head import moment_head_axis
    sh = moment_head_axis(FunctionSpace.of_basis(_sh(L)))
    assert sh.shape == (L + 1, 2 * L + 1) and sh.L == L and sh.isotropic_slot == (0, 0)
    for l in range(L + 1):
        assert sh.degree_block(l) == (l, slice(0, 2 * l + 1))
    lg = moment_head_axis(FunctionSpace.of_basis(_legendre(L, "x")))
    assert lg.shape == (L + 1,) and lg.L == L and lg.isotropic_slot == (0,)
    assert [lg.degree_block(l) for l in range(L + 1)] == [(l,) for l in range(L + 1)]
    with pytest.raises(ValueError):
        sh.degree_block(L + 1)


@_S7B
def test_H5_moment_head_axis_refuses_a_non_head_with_the_shared_message():
    from orpheus.numerics.spaces.moment_head import moment_head_axis
    q = Quadrature.gauss_legendre(4)
    from orpheus.numerics.axis import Axis
    with pytest.raises(TypeError, match="is not a moment head"):
        moment_head_axis(FunctionSpace.of_axes(Axis.of_quadrature(q)))


@_S7B
def test_H6_legendre_about_two_axes_differ_by_domain_not_by_a_spent_axis_field():
    """Ruling 81: spent_axis dissolves into the domain S²/O2_a, which ruling 71 makes content."""
    a, b = FunctionSpace.of_basis(_legendre(1, "x")), FunctionSpace.of_basis(_legendre(1, "z"))
    _separated(a, b, "Legendre about x vs about z")
    assert not hasattr(a.axes[0], "spent_axis")


@pytest.mark.parametrize("rule,L", [("gl8", 2), ("product48", 2), ("folded24", 2)])
def test_H7_truncated_head_is_the_frame_head_at_the_lower_order(rule, L):
    """truncated_head(space, L_new) equals the dressed head the frame mints at L_new (one body, design §4)."""
    from orpheus.numerics.spaces.moment_head import truncated_head
    q = {"gl8": lambda: Quadrature.gauss_legendre(8), "product48": lambda: Quadrature.product(4, 8),
         "folded24": lambda: Quadrature.folded_product(n_mu=2, n_phi=4)}[rule]()
    head = ff.basis_space(ff.angular_frame(q, L))
    for L_new in range(L + 1):
        assert truncated_head(head, L_new) == ff.basis_space(ff.angular_frame(q, L_new))


@_S7B
def test_H8_the_three_pending_concept_edges_are_closed():
    """G0: basis/ and measure.py import no space module (the three S7 PENDING entries)."""
    from tests.gates.numerics import test_concept_order as G0
    pending = getattr(G0, "PENDING")
    s7 = {k: v for k, v in pending.items() if "S7" in str(v)} if isinstance(pending, dict) else pending
    assert not s7, f"S7 PENDING edges remain: {s7}"

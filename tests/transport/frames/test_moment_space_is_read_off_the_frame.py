r"""#429 tracker 2.5 — the angular moment space is READ off the frame, never minted from ``L``.

The defect this gates (`[M]` 2026-09-02, the fused step's opener): the
angular moment space had EIGHT homes. The frame carried it
(``frame.basis.space``) and seven production sites re-minted it from the
integer ``L`` as ``SphericalHarmonicSpace.from_L(L)`` — the scattering
operator's :math:`\Lambda` ends (three spellings), the fission and (n,2n)
:math:`\ell = 0` ends, the moment-flux field's head, and ``truncate`` —
while two ``isinstance`` doors on :class:`HarmonicFrame` admitted one basis
class only. Every copy silently chose the full-sphere family; the day a
1-D rule binds its Legendre basis, every copy mismatches the frame at the
``(name, shape)`` composability guards, and the FIRST mismatch is at
:math:`L = 0` on every solve (fission and (n,2n) mint there).

Three gates, each with the input it rejects (``plan-authoring`` §6c):

* **the ROUTE gate** — bind a FOREIGN truncated basis (not a spherical
  harmonic subclass at all, with a RENAMED coefficient space) into the
  quadrature's frame and require every operator end and every moment
  field space to MOVE with it. A reverted producer (one that still mints
  from ``L``) fails the composability guard loudly, which is the red.
  The mutant is unconstructible before the door widened — so the door
  and the producers are ONE step (``plan-authoring`` §6b);
* **the METRIC-IDENTITY gate** — Landing A binds the basis's own
  CONTINUUM space (``basis.space``), not the frame's Parseval-dressed
  ``basis_space``: `[M]` the two are ``(name, shape)``-equal and
  metric-DIFFERENT on 33 of 33 shipped (rule, L) rows (the per-ℓ ratio is
  exactly ``[(2ℓ+1)/4π]²``), and Λ's Hilbert adjoint under the continuum
  end is its transpose exactly while the dressed end would move it on 10
  of 33 rows. A ``(name, shape)`` equality cannot see that fork, so this
  gate asserts the metric ARRAY, with the dressed space as its negative
  control;
* **the DOOR gate** — the frame asks for the
  :class:`~orpheus.numerics.basis.base.TruncatedBasis` SURFACE (a
  truncation order), typed, at the door: an indicator trial is refused
  there, and a foreign truncated basis is admitted where the old door
  refused it.

Foundation mark: software invariants (route, metric identity, typing),
mutation-proven; no physics claim rides here.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.basis.base import Basis, GramStructure, TruncatedBasis
from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.frame import GalerkinFrame
from orpheus.numerics.manifold import SPHERE, Manifold, RealSpace
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import IncompatibleOperatorComposition, OperatorProduct
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace, TensorProductSpace
from orpheus.numerics.spaces.spherical_harmonic_space import SphericalHarmonicSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.frames import HarmonicFrame
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.n2n import N2NOperator
from orpheus.transport.operators.scattering import (
    LegendreMomentScattering,
    N2NMomentOperator,
    ScatteringOperator,
)
from tests.sn._test_helpers import material_xs_from_raw, placeholder_materials

pytestmark = pytest.mark.foundation

_MUTANT_NAME = "mutated_moment_space"
_SIGS0 = np.array([[0.20, 0.00], [0.05, 0.18]])
_SIGS1 = np.array([[0.02, 0.00], [0.01, 0.015]])
_SIG2 = np.array([[0.00, 0.03], [0.01, 0.00]])


# ── the FOREIGN truncated basis: carries L and a renamed space, is NOT an SH subclass ──


@dataclass(frozen=True)
class _ForeignTruncatedBasis(Basis):
    """The same functions as the real harmonics, delegated — with a RENAMED
    coefficient space, and NOT a ``SphericalHarmonicBasis`` subclass. It
    carries a truncation order, which is the ONE thing the door asks for;
    its renamed space is ``(name, shape)``-unequal to ``from_L(L)``, which
    is what the route gate reads."""

    L: int

    def _parent(self) -> SphericalHarmonicBasis:
        return SphericalHarmonicBasis(L=self.L)

    def evaluate(self, points, /):
        return self._parent().evaluate(points)

    def synthesize(self, coefficients, table, /):
        return self._parent().synthesize(coefficients, table)

    def analyze(self, values, table, weights, /):
        return self._parent().analyze(values, table, weights)

    def analyze_transpose(self, coefficients, table, weights, /):
        return self._parent().analyze_transpose(coefficients, table, weights)

    def reconstruct(self, coefficients, table, /):
        return self._parent().reconstruct(coefficients, table)

    def reconstruct_transpose(self, values, table, /):
        return self._parent().reconstruct_transpose(values, table)

    def mass_matrix(self, measure, /):
        return self._parent().mass_matrix(measure)

    @property
    def gram_structure(self) -> GramStructure:
        return GramStructure.DIAGONAL

    @property
    def domain(self) -> Manifold:
        return SPHERE

    @property
    def space(self) -> FunctionSpace:
        return replace(SphericalHarmonicSpace.from_L(self.L), name=_MUTANT_NAME)


def _slab(nx: int = 4, n_ord: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(n_ordinates=n_ord), placeholder_materials(ng=ng))


def _mat_xs(nx: int = 4):
    cells = {0: (np.arange(nx), np.zeros(nx, dtype=int))}
    return material_xs_from_raw(
        sig_s={0: [_SIGS0, _SIGS1]}, sig2={0: _SIG2},
        cells_by_mat=cells, ng=2, nx=nx, ny=1,
    )


def _inner_factor_domain_name(kernel: OperatorProduct) -> str:
    """``frame.conjugate(X)`` is ``R ∘ (X ∘ M)``; the domain of ``X`` is the end under test."""
    inner = kernel.b
    assert isinstance(inner, OperatorProduct)
    domain = inner.a.domain
    assert domain is not None
    return domain.name


def _bind_foreign(sn: SNMesh, L: int) -> HarmonicFrame:
    """Install the foreign basis as the quadrature's frame at ``L`` through
    the production chain's own cache — so ``HarmonicFrame.for_space`` (the
    ONE spelling every consumer uses) hands the foreign frame out."""
    # the foreign frame rides the SAME measure the production frame carries
    # (today the forged 3-D padding of a 1-D rule; after #429's fix the rule's
    # own) — the test is about the SPACE, not the nodes
    honest = sn.quad.angular_frame(L)
    cache = sn.quad._angular_frames or {}
    cache[L] = GalerkinFrame(_ForeignTruncatedBasis(L=L), honest.measure)
    object.__setattr__(sn.quad, "_angular_frames", cache)
    frame = HarmonicFrame.for_space(sn.angular_bulk_space, L)
    assert frame.basis.space.name == _MUTANT_NAME, "the route through for_space did not reach the foreign frame"
    return frame


# ═══════════════════════════════════════════════════════════════════════
# A1 — the ROUTE gate
# ═══════════════════════════════════════════════════════════════════════


def test_swapping_the_frames_basis_moves_every_operator_end_and_field_space() -> None:
    """Every producer that used to mint ``SphericalHarmonicSpace.from_L(L)``
    now reads the frame's basis: with a FOREIGN basis bound, the
    scattering, fission and (n,2n) kernels COMPOSE (a reverted producer
    would fail the ``(name, shape)`` guard), the moment field's head is the
    foreign space, and the tier-2 mints take the basis, not the integer."""
    sn = _slab()
    L = 1
    _bind_foreign(sn, L)
    _bind_foreign(sn, 0)          # fission and (n,2n) mint at ℓ = 0 on every solve
    mat = _mat_xs()
    composite = sn.full_field_space

    S = ScatteringOperator.from_solver_data(mat_xs=mat, space=composite, scattering_order=L)
    lam = S._moment_scattering(skip_l0=False)
    assert lam.domain.name == _MUTANT_NAME and lam.codomain.name == _MUTANT_NAME
    assert S.full_scatter_kernel is not None      # R∘Λ∘M composes on the foreign ends
    assert S.kernel is not None

    F = FissionOperator.from_solver_data(mat_xs=mat, space=composite)
    assert _inner_factor_domain_name(F.full_fission_kernel) == _MUTANT_NAME   # R∘(F₀∘M)

    N = N2NOperator.from_solver_data(mat_xs=mat, space=composite)
    assert _inner_factor_domain_name(N.full_n2n_kernel) == _MUTANT_NAME

    field = HarmonicMomentFlux.zeros_for_mesh_and_L(sn, L)
    assert isinstance(field.space, TensorProductSpace)
    assert field.space.factors[0].name == _MUTANT_NAME
    # truncation stays in the head's OWN family, one order down — the head
    # keeps its IDENTITY (name) and only its order moves; a truncate that
    # re-minted the family from an integer would hand back the default name
    truncated = field.truncate(0).space
    assert isinstance(truncated, TensorProductSpace)
    assert truncated.factors[0].shape == (1, 1)
    assert truncated.factors[0].name == _MUTANT_NAME

    foreign = _ForeignTruncatedBasis(L=L)
    assert LegendreMomentScattering.from_material_xs(mat, foreign).domain.name == _MUTANT_NAME
    assert N2NMomentOperator.from_material_xs(mat, foreign).codomain.name == _MUTANT_NAME

    # the negative leg: an end minted from L alone does NOT compose with the
    # foreign frame — this is the red every reverted producer produces
    minted = SphericalHarmonicSpace.from_L(L)
    stale = LegendreMomentScattering(
        S.scattering, skip_l0=False, domain=minted, codomain=minted,
    )
    with pytest.raises(IncompatibleOperatorComposition, match=r"A\.domain == B\.codomain"):
        S.flux_analysis.frame.conjugate(stale)


# ═══════════════════════════════════════════════════════════════════════
# A2 — the METRIC-IDENTITY gate (and the fork control)
# ═══════════════════════════════════════════════════════════════════════


_RULES = {
    "gauss_legendre(2)": lambda: Quadrature.gauss_legendre(2),
    "gauss_legendre(8)": lambda: Quadrature.gauss_legendre(8),
    "gauss_legendre(16)": lambda: Quadrature.gauss_legendre(16),
    "level_symmetric(4)": lambda: Quadrature.level_symmetric(4),
    "level_symmetric(8)": lambda: Quadrature.level_symmetric(8),
    "lebedev(11)": lambda: Quadrature.lebedev(11),
    "lebedev(17)": lambda: Quadrature.lebedev(17),
    "product(4,6)": lambda: Quadrature.product(4, 6),
    "product(8,8)": lambda: Quadrature.product(8, 8),
    "folded_product(2,4)": lambda: Quadrature.folded_product(2, 4),
    "folded_product(4,8)": lambda: Quadrature.folded_product(4, 8),
}


@pytest.mark.parametrize("label", sorted(_RULES))
@pytest.mark.parametrize("L", [0, 1, 2])
def test_landing_a_binds_the_continuum_space_and_is_metric_identical(label: str, L: int) -> None:
    """The space the operator ends and field heads now read is the basis's
    OWN (continuum-metric) space — ``(name, shape)``-equal AND
    ``array_equal`` on the metric to the ``from_L(L)`` mint it replaces —
    and NOT the frame's Parseval-dressed ``basis_space``, which is
    ``(name, shape)``-equal and metric-different (the fork a ``==`` gate
    cannot see)."""
    q = _RULES[label]()              # built in the body, never in the parametrize list
    frame = HarmonicFrame.from_galerkin(q.angular_frame(L))
    ends = frame.basis.space
    minted = SphericalHarmonicSpace.from_L(L)
    assert ends == minted
    assert ends.inner_product_weights is not None and minted.inner_product_weights is not None
    np.testing.assert_array_equal(ends.inner_product_weights, minted.inner_product_weights)
    dressed = frame.basis_space
    assert dressed == ends                       # equality is metric-blind ...
    if dressed.metric is not None:               # ... the DENSE arm installs a matrix metric
        assert dressed.inner_product_weights is None
    else:                                        # ... the diagonal arm installs the Parseval inverse
        assert dressed.inner_product_weights is not None
        assert not np.array_equal(dressed.inner_product_weights, ends.inner_product_weights)


def test_the_operator_ends_carry_the_continuum_metric_not_the_dressed_one() -> None:
    """On a real posed composite the scattering factor's ends carry the
    basis's continuum Gram bit-for-bit (`[M]` 2026-09-02: the dressed
    metric would move ``apply_metric`` by 96–161 %)."""
    sn = _slab()
    L = 1                      # the synthetic data carries P0 and P1
    S = ScatteringOperator.from_solver_data(mat_xs=_mat_xs(), space=sn.full_field_space, scattering_order=L)
    ends = S._moment_scattering(skip_l0=False).domain
    minted = SphericalHarmonicSpace.from_L(L)
    assert ends.inner_product_weights is not None and minted.inner_product_weights is not None
    np.testing.assert_array_equal(ends.inner_product_weights, minted.inner_product_weights)
    assert ends.metric is None


# ═══════════════════════════════════════════════════════════════════════
# A3 — the DOOR gate
# ═══════════════════════════════════════════════════════════════════════


def test_the_door_asks_for_a_truncation_order_not_for_one_class() -> None:
    """An indicator trial (no truncation order) is refused at BOTH doors
    with a message naming the surface; a foreign truncated basis — which
    the old ``isinstance(basis, SphericalHarmonicBasis)`` door refused —
    is admitted, and the frame reads its order."""
    indicator = IndicatorBasis((np.array([0.0, 1.0, 2.0]),), RealSpace(1))
    measure = DiscreteMeasure(
        nodes=np.array([0.5, 1.5]), weights=np.ones(2), support=RealSpace(1),
    )
    assert not isinstance(indicator, TruncatedBasis)
    with pytest.raises(TypeError, match="truncation order"):
        HarmonicFrame(indicator, measure)
    with pytest.raises(TypeError, match="truncation order"):
        HarmonicFrame.from_galerkin(GalerkinFrame(indicator, measure))

    foreign = _ForeignTruncatedBasis(L=1)
    assert isinstance(foreign, TruncatedBasis)
    assert not isinstance(foreign, SphericalHarmonicBasis)
    frame = HarmonicFrame(foreign, Quadrature.gauss_legendre(4).measure)
    assert frame.truncation_order == 1
    upgraded = HarmonicFrame.from_galerkin(GalerkinFrame(foreign, Quadrature.gauss_legendre(4).measure))
    assert upgraded.truncation_order == 1

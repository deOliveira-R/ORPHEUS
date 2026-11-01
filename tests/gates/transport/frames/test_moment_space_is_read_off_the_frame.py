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
* **the METRIC gate** — ⛔ RE-POSED at CS4c step 6 item 6.2c-ii (ruling
  R-6.2c-1, 2026-09-08). Landing A (2026-09-02) bound the basis's own
  CONTINUUM space (``basis.space``) and recorded as its reason that *"Λ's
  Hilbert adjoint under the continuum end is its transpose exactly while
  the dressed end would move it on 10 of 33 rows (the dense-Gram rows)"*
  and that the dressed metric *"would move apply_metric by 96–161 %"*.
  `[M]` re-measured (``scratch/_step6_2c/p3_scan_161.py``,
  ``p4_lambda_adjoint.py``, ``p9_on_range.py``, 33 shipped (rule, L)
  rows): the adjoint moves on **5** rows for an ARBITRARY head draw,
  **3 of them DIAGONAL-Gram** (the mechanism is the Parseval metric's
  Moore–Penrose projection of the slots a folded rule cannot see, not
  Gram density), and on **0 of 33** for a physical moment ``φ = Mψ``;
  and NO statistic reproduces 161 % (the draw-free per-element movement
  spans 0.5 %…100 % over the 33 rows). What the continuum binding COST
  was measured too: Parseval ``‖Mψ‖² = ‖ψ‖²_W`` fails on 33 of 33 rows
  under it (ratio 3.41…157.91) and holds on 33 of 33 under the dressed
  metric. The tree therefore binds the FRAME's Parseval-dressed
  ``basis_space`` everywhere — the operator ends, the moment fields, the
  carrier's cached moment space — and the two heads are structurally
  UNEQUAL (the head is axis-built; its measure enters the identity), so
  this gate asserts the ends ARE the frame's space and are NOT the
  continuum mint, with the continuum head as its negative control;
* **the DOOR gate** — the frame asks for the
  :class:`~orpheus.numerics.basis.base.TruncatedBasis` SURFACE (a
  truncation order), typed, at the door: an indicator trial is refused
  there, and a foreign truncated basis is admitted where the old door
  refused it.

Foundation mark: software invariants (route, metric identity, typing),
mutation-proven; no physics claim rides here.

Layer DAG Phase 4 (2026-11-01): the two rows that read the SN problem's moment space moved to
tests/gates/sn/mesh/test_moment_space_is_read_off_the_frame_on_the_problem.py.
"""
from __future__ import annotations

from orpheus.numerics.basis.legendre_basis import LegendreBasis

from dataclasses import dataclass, replace

import numpy as np
import pytest

from orpheus.numerics.basis.base import Basis, GramStructure, TruncatedBasis
from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.metric import DenseMetric, FactoredMetric
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.numerics.frame import GalerkinFrame
from orpheus.geometry.manifold import Manifold, RealSpace
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.frames import HarmonicFrame
from orpheus.numerics import frame_faces

pytestmark = pytest.mark.foundation

#: The foreign twin's one content change: its analysis (and so the dressed
#: head's unit, the moments of 1) and its continuum Gram are scaled by this.
_MUTANT_SCALE = 2.0


# ── the FOREIGN truncated basis: carries L and a moved content, is NOT an SH subclass ──


@dataclass(frozen=True)
class _ForeignTruncatedBasis(Basis):
    """The same functions as the rule's OWN basis, delegated — with its
    analysis and its continuum Gram scaled by ``_MUTANT_SCALE``, and NOT a
    ``SphericalHarmonicBasis`` subclass. It carries a truncation order, which
    is the ONE thing the door asks for. Its coefficient spaces, the continuum
    one (:meth:`FunctionSpace.of_basis`) and the frame-dressed one, are
    CONTENT-unequal to the honest ones (the weights, and the dressed head's
    unit, the moments of 1), which is what the route gate reads. Until S7 the
    twin renamed its space instead; since S7 a space is its content and a
    basis names no space.

    ⛔ **The parent is a FIELD since 2026-09-02 (#429).** It used to hard-code
    ``SphericalHarmonicBasis(L)``, hence ``domain = SPHERE`` — and the frame's
    G0 arrow now refuses a basis on :math:`S^2` bound to a rule whose measure
    lives on :math:`S^2/O(2)_x`, which is the whole of ERR-080. Delegating
    the domain keeps the mutant admissible **wherever the honest basis is**,
    so the ONE thing that differs between the two runs stays the coefficient
    SPACE — which is what this module's route gate is about (``vv-principles``
    #18: a mutation that also breaks a structural law over-states its own
    coverage).
    """

    L: int
    parent: Basis

    @classmethod
    def like(cls, basis: Basis) -> "_ForeignTruncatedBasis":
        """A foreign twin of ``basis``: same functions, same domain, renamed space."""
        order = getattr(basis, "L", None)
        assert isinstance(order, int)
        return cls(L=order, parent=basis)

    def _parent(self) -> Basis:
        # the parent FOLLOWS this twin's order, so a re-mint at another order
        # (``at_order``, item 6.2c-ii — the head truncates through it) still
        # delegates to the honest family at THAT order
        return self.parent.at_order(self.L) if getattr(self.parent, "L", None) != self.L else self.parent  # type: ignore[attr-defined]

    def at_order(self, L_new: int, /) -> "_ForeignTruncatedBasis":
        return replace(self, L=L_new)

    def evaluate(self, points, /):
        return self._parent().evaluate(points)

    def _synthesize(self, coefficients, table, /):
        return self._parent().synthesize(coefficients, table)

    def _analyze(self, values, table, weights, /):
        return _MUTANT_SCALE * self._parent().analyze(values, table, weights)

    def _analyze_transpose(self, coefficients, table, weights, /):
        return self._parent().analyze_transpose(coefficients, table, weights)

    def _reconstruct(self, coefficients, table, /):
        return self._parent().reconstruct(coefficients, table)

    def _reconstruct_transpose(self, values, table, /):
        return self._parent().reconstruct_transpose(values, table)

    @property
    def gram_structure(self) -> GramStructure:
        return self._parent().gram_structure

    @property
    def domain(self) -> Manifold:
        return self._parent().domain

    @property
    def layout(self):
        return self._parent().layout

    @property
    def coefficient_weights(self):
        return _MUTANT_SCALE * self._parent().coefficient_weights


# ═══════════════════════════════════════════════════════════════════════
# A1 — the ROUTE gate
# ═══════════════════════════════════════════════════════════════════════


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


def _head_measure(space: FunctionSpace) -> np.ndarray | None:
    """The head's diagonal metric where the axis doctrine put it — the single head axis's measure."""
    assert space.axes is not None and len(space.axes) == 1, "a moment head is a single-axis space"
    return space.axes[0].weights


@pytest.mark.parametrize("label", sorted(_RULES))
@pytest.mark.parametrize("L", [0, 1, 2])
def test_the_tree_binds_the_frames_parseval_space_and_not_the_continuum_head(label: str, L: int) -> None:
    """RE-POSED at CS4c step 6 item 6.2c-ii (ruling R-6.2c-1): the space the
    operator ends and field heads read is the FRAME's Parseval-dressed
    ``basis_space`` — axis-built, the frame its head axis's generator — and
    it is structurally UNEQUAL to the basis's own continuum head (the
    ``from_L(L)`` mint, ``basis.space``), which is this gate's negative
    control: same family, same order, same name, different measure ⟹ a
    different space (the metric-blind ``(name, shape)`` seam Landing A
    leaned on is gone).

    The metric is asserted on the ARRAY per row: on a DIAGONAL-Gram row the
    head axis's measure is the Moore–Penrose reciprocal of the discrete Gram's
    diagonal (zero on dead slots); on a DENSE-Gram row (`[M]` 2 of these 33:
    ``gauss_legendre(2)`` and ``folded_product(2,4)`` at L = 2) the axis
    carries no measure and the matrix pseudo-inverse is POSITIONED on the
    space's derived metric object (item 6.2c-i)."""
    q = _RULES[label]()              # built in the body, never in the parametrize list
    frame = HarmonicFrame.from_galerkin(frame_faces.angular_frame(q, L))
    dressed = frame_faces.basis_space(frame)
    # ⛔ RE-KEYED 2026-09-02 (#429): the family to compare against is the one
    # the rule's ORBIT SPACE carries — a 1-D rule's frame binds the flat
    # Legendre head, a sphere rule's the rectangular harmonics. ⛔ RE-KEYED at
    # S7 (ruling 71, #631 pair 3): the folded rule's head is the MIRROR-EVEN
    # family on S²/σ_y, which the head axis's domain now separates from the
    # plain harmonics on S² (until S7 the two compared equal).
    from orpheus.numerics.basis.spherical_harmonic_basis import MirrorEvenSphericalHarmonicBasis

    if label.startswith("gauss_legendre"):
        continuum = FunctionSpace.of_basis(LegendreBasis(L=L, axis="x"))
    elif label.startswith("folded_product"):
        continuum = FunctionSpace.of_basis(MirrorEvenSphericalHarmonicBasis(L=L, mirror_axis=1))
    else:
        continuum = FunctionSpace.of_basis(SphericalHarmonicBasis(L=L))
    assert FunctionSpace.of_basis(frame.basis) == continuum, "the basis's own space IS the continuum mint"
    assert dressed.name == continuum.name and dressed.shape == continuum.shape
    assert dressed != continuum and continuum != dressed, (
        "the dressed head and the continuum head must be two spaces — the metric is the identity"
    )
    assert dressed.axes is not None and len(dressed.axes) == 1
    assert dressed.axes[0].origin.frame is frame, "the frame dresses the head and becomes its origin"
    assert dressed.inner_product_weights is None
    diag = np.diagonal(frame.discrete_gram).reshape(dressed.shape)
    live = diag > 0.0
    if frame.discrete_gram_structure is GramStructure.DENSE:
        assert (label, L) in {("gauss_legendre(2)", 2), ("folded_product(2,4)", 2)}
        assert _head_measure(dressed) is None, "a dense Gram leaves the axis measure-less"
        assert isinstance(dressed.metric, FactoredMetric) and len(dressed.metric.entries) == 1
        assert isinstance(dressed.metric.entries[0][1], DenseMetric)
    else:
        assert dressed.metric is None
        measure = _head_measure(dressed)
        assert measure is not None
        np.testing.assert_allclose(measure[live], 1.0 / diag[live], rtol=1e-15)
        np.testing.assert_array_equal(measure[~live], 0.0)
        # the negative control: the continuum head's measure is the OTHER array
        cont = _head_measure(continuum)
        assert cont is not None and not np.array_equal(cont, measure)


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

    # ⚠ The foreign basis is built as a twin of the RULE'S OWN basis
    # (2026-09-02, #429): the frame's G0 arrow refuses a basis on
    # :math:`S^2` bound to a 1-D rule's :math:`S^2/O(2)_x`, so a mutant
    # hard-coding the harmonics would red for the WRONG reason and this gate
    # would stop measuring the door.
    rule = Quadrature.gauss_legendre(4)
    foreign = _ForeignTruncatedBasis.like(frame_faces.angular_frame(rule, 1).basis)
    assert isinstance(foreign, TruncatedBasis)
    assert not isinstance(foreign, SphericalHarmonicBasis)
    frame = HarmonicFrame(foreign, rule.measure)
    assert frame.truncation_order == 1
    upgraded = HarmonicFrame.from_galerkin(GalerkinFrame(foreign, rule.measure))
    assert upgraded.truncation_order == 1

    # …and the same twin on a SPHERE rule, so the door gate is not keyed to
    # one family either (the foreign basis there wraps the harmonics).
    sphere = Quadrature.level_symmetric(4)
    foreign_sphere = _ForeignTruncatedBasis.like(frame_faces.angular_frame(sphere, 1).basis)
    assert HarmonicFrame(foreign_sphere, sphere.measure).truncation_order == 1

r"""Foundation tests for :attr:`DiscreteMeasure.phase` — the phase-space factor.

A measure discretises exactly one factor of transport phase space
(position × direction × energy). :attr:`phase` is a **derived** category,
read from the measure's structure (symmetry group for the angular factor,
support tag for spatial), NOT a hand-set literal. These tests pin that
derivation and the documented per-category seam (energy / untagged → raises
until its typed machinery is minted).

Tagged ``foundation`` — software invariants on the type surface.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasure, gauss_chebyshev, gauss_legendre
from orpheus.numerics.manifold import SPHERE, Quotient, RealSpace
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3


# ── angular: derived from the O(3) invariance group (the Erlangen route) ──

@pytest.mark.foundation
@pytest.mark.parametrize(
    "make",
    [
        lambda: Quadrature.lebedev(order=17).measure,        # O_h sphere cubature
        lambda: Quadrature.level_symmetric(sn_order=4).measure,  # O_h triangular
        lambda: Quadrature.gauss_legendre(8).measure,        # SO(2) slab μ-rule
        lambda: Quadrature.product(n_mu=4, n_phi=6).measure,  # SO(2) product
    ],
)
def test_angular_quadratures_have_angular_phase(make):
    """Every angular quadrature carries an ``O(3)``-subgroup invariance, from
    which ``phase == "angular"`` is derived — the symmetry IS the signature."""
    m = make()
    assert m.invariance_group is not None
    assert isinstance(m.invariance_group, SubgroupOfO3)
    assert m.phase == "angular"


@pytest.mark.foundation
def test_angular_phase_survives_stripping_the_symmetry_group():
    r"""⛔ **This gate's claim INVERTED at tracker 2.0c, and the inversion is
    the point** — do not read the old name as still describing it.

    It used to assert that ``phase`` reads the *group*, not the support: strip
    ``invariance_group`` from a Lebedev rule and ``phase`` RAISED, because the
    only angular route was ``invariance_group is not None``.  ``"S^2"`` was an
    inert string that no branch consulted.

    Typing ``support`` makes :class:`~orpheus.numerics.manifold.Sphere` a
    *reason*: the direction variable lives on the sphere, so a measure on one
    is angular whether or not the rule happens to be closed under a subgroup.
    That is what lets the shipped fold answer at all —
    ``folded_product(4, 8)`` has ``invariance_group=None`` **legitimately**
    (the mirror was quotiented away) and used to raise.

    ⟹ the honest claim is now the opposite of the old one: the manifold is the
    primary route and the group is the tie-breaker for the one genuinely
    ambiguous point set (an ``Interval``), which
    :func:`test_untagged_generic_rule_phase_raises` still pins.
    """
    base = Quadrature.lebedev(order=11).measure
    assert base.phase == "angular"
    from dataclasses import replace
    stripped = replace(base, invariance_group=None)
    assert stripped.support == SPHERE
    assert stripped.phase == "angular"


# ── spatial: derived from the support tag (a different category) ──────────

@pytest.mark.foundation
@pytest.mark.parametrize("d", [1, 2, 3])
def test_a_measure_on_R_d_has_spatial_phase(d):
    r"""Spatial measures carry NO :math:`O(3)` symmetry (a mesh is not
    rotation-invariant) — a different phase-space factor, and now recognised
    by the TYPE of the point set rather than by a tag's prefix.

    ⭐ Two things changed at tracker 2.0c, both visible in this test's own
    parameters. (1) The list was ``["spatial_R1", "spatial_R2", "cells"]``;
    ``"cells"`` was a support with **zero production producers** whose only
    producer was *this parametrize list* — a synthetic witness for an arm that
    existed to catch it (``[M]`` §V.5e(e)). It retires with the arm, per
    ``coding-standards`` ("retirement means the gate row goes too"). (2) The
    string version tested ``support.startswith("spatial")``, so a third
    dimension needed a new tag; ``RealSpace(3)`` needs nothing, which is why
    ``d = 3`` is here and passes without a line of production change.
    """
    m = DiscreteMeasure(
        nodes=np.zeros((4, d)) if d > 1 else np.linspace(0.0, 1.0, 4),
        weights=np.full(4, 0.25),
        support=RealSpace(d),
    )
    assert m.invariance_group is None
    assert m.phase == "spatial"


# ── the per-category seam: undetermined factors raise ─────────────────────

@pytest.mark.foundation
@pytest.mark.parametrize("make", [lambda: gauss_legendre(5), lambda: gauss_chebyshev(4)])
def test_untagged_generic_rule_phase_raises(make):
    """A bare generic 1-D rule (no symmetry group, no spatial tag) has no
    determined phase — a slab μ-interval is geometrically indistinguishable
    from a spatial interval, so the physical identity must be supplied. The
    energy factor lands here too, until it earns its typed support."""
    with pytest.raises(NotImplementedError, match="phase is undetermined"):
        _ = make().phase


@pytest.mark.foundation
def test_phase_is_a_closed_category_consistent_with_angular_frame():
    """``phase`` agrees with the ``angular_frame`` axis name — the angular
    quadrature's frame and its measure both say "angular".

    ⭐ The last assertion is the one this docstring always claimed and could
    not make. Until 2026-09-01 ``angular_frame`` REBUILT its measure from bare
    nodes + weights + a literal support, dropping the ``invariance_group``
    that ``phase`` derives from — so ``frame.measure.phase`` raised
    ``NotImplementedError`` on all 12 shipped rules, including this one, whose
    own measure answers ``"angular"`` correctly one line above. The assertable
    substitute was the support TAG, which is what the body settled for.

    Phase 0.1a of ``.claude/plans/angular_spaces_derived_from_symmetry.md``
    routes the rule's own measure into the frame, so the claim is now true and
    checkable. Both assertions are kept: the tag one is the weaker statement
    and remains a genuine pin on the routing.
    """
    q = Quadrature.lebedev(order=17)
    assert q.measure.phase == "angular"
    # the frame built on it is the angular frame on the same (angular) measure
    assert q.angular_frame(2).measure.support == SPHERE
    assert q.angular_frame(2).measure.phase == "angular"


# ── the fold: what the measure HAS vs what it has SPENT ───────────────────

@pytest.mark.foundation
def test_the_folded_rule_reports_its_phase_and_its_spent_group() -> None:
    r"""⛔ **A LIVE DEFECT's witness** — this raised until 2026-09-01.

    ``folded_product(4, 8)`` lives on :math:`S^2/\sigma_y` with
    ``invariance_group=None``, and the ``None`` is **correct**: the mirror was
    quotiented AWAY, so the folded nodes are genuinely not closed under it.
    The retired string dispatch had exactly three arms — ``invariance_group is
    not None``, a ``"spatial…"`` prefix, an ``"energy"`` prefix — and a
    quotient of the sphere matches none of them, so asking the shipped
    production fold which factor of phase space it discretises **raised
    ``NotImplementedError``** (``[M]`` §V.5e(f)).

    It stayed latent only because ``.phase`` had zero production consumers at
    the time — which is not a defence (``process-discipline``: "no current
    consumer" ≠ speculative), and the same property raised once before this
    campaign, on all 12 rules, from a different cause (0.1a).

    Dispatching on the manifold's TYPE fixes it structurally rather than by
    adding a fourth string prefix: a ``Quotient`` of a ``Sphere`` is the
    direction variable's home whatever symmetry the rule spent.
    """
    folded = Quadrature.folded_product(n_mu=4, n_phi=8).measure

    assert folded.phase == "angular"
    assert isinstance(folded.support, Quotient)
    assert folded.support.base == SPHERE

    # ⭐ HAS vs SPENT — the two are complementary here, and conflating them is
    # what made the string dispatch look total.
    assert folded.invariance_group is None          # spent: no longer invariant
    assert folded.quotient_group == SubgroupOfO3.Mirror("y")   # what it spent

    # `quotient_group` is DERIVED, and this is a second, independent record of
    # the same fact: `Quadrature.folded_by` is stored on the rule, not read off
    # the measure's support.
    assert folded.quotient_group == Quadrature.folded_product(
        n_mu=4, n_phi=8).folded_by

    # Negative leg: an unfolded rule has spent nothing (vv-principles #11).
    unfolded = Quadrature.product(n_mu=4, n_phi=8).measure
    assert unfolded.quotient_group is None
    assert unfolded.invariance_group is not None

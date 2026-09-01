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

from orpheus.numerics.measure import (
    DiscreteMeasure,
    SPACE_SPHERE,
    gauss_chebyshev,
    gauss_legendre,
)
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
def test_angular_phase_is_symmetry_derived_not_support_tagged():
    """The derivation reads the *group*, not the support string: a measure on
    ``S^2`` with the group set is angular, and removing the group flips it
    (no support-string fallback masks the symmetry route)."""
    base = Quadrature.lebedev(order=11).measure
    assert base.phase == "angular"
    # Strip the group: an S^2-tagged measure with no symmetry is no longer
    # auto-angular — the group, not the "S^2" string, is the discriminator.
    from dataclasses import replace
    stripped = replace(base, invariance_group=None)
    with pytest.raises(NotImplementedError):
        _ = stripped.phase


# ── spatial: derived from the support tag (a different category) ──────────

@pytest.mark.foundation
@pytest.mark.parametrize("support", ["spatial_R1", "spatial_R2", "cells"])
def test_spatial_supports_have_spatial_phase(support):
    """Spatial measures carry NO O(3) symmetry (a mesh is not rotation-
    invariant) — a different phase-space factor, recognised by its support
    tag until its typed machinery is minted."""
    m = DiscreteMeasure(
        nodes=np.linspace(0.0, 1.0, 4),
        weights=np.full(4, 0.25),
        support=support,
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
    assert q.angular_frame(2).measure.support == SPACE_SPHERE
    assert q.angular_frame(2).measure.phase == "angular"

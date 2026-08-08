r"""Cylindrical quadrature ADMISSION — ``SNMesh(CYLINDRICAL)`` constructs iff
every μ-level is CARRYING (Q5.6 step 6.3, the flip).

The admission tier sits BETWEEN two older guards, and the three fire in a
fixed order on the constructor path:

1. **trace-rank** (``numerics/spaces/angular_trace_space.build_omega_dot_n``)
   — refuses a rule with no genuine μ_x cosine (the all-tangential
   ``folded_product(n, 2)`` class; gate_design §8.1);
2. **structure-less** (``geometry/reduced_operator.cylindrical_streaming``,
   fragment ``level structure``) — refuses slab/sphere cubatures that carry
   no ``LevelStructure`` at all (gated by
   ``tests/sn/sweep/core/test_sweep_regression.py`` and
   ``tests/sn/sweep/curvilinear/test_cyl_sweep_regression.py``);
3. **carrying** (``sn/sweep/pole_angular_closure.assert_carrying_quadrature``,
   fragment ``non-carrying``) — THIS module's subject: the rule HAS levels,
   but a level's march start is an ordinate (``on_edge_node``) or its seed
   thread weight vanishes (``degenerate``), so route (a) has no honest ψ½
   state to march.

Admission is decided by STRUCTURE, not provenance: the facts are
:func:`~orpheus.sn.sweep.pole_angular_closure.march_start_structure_per_level`'s
— the same producer ``SNMesh.radial_characteristic_levels`` reads.  The
upstream quadrature-level classification battery is
``tests/sn/sweep/test_march_start_structure.py``; this module gates the
MESH's consumption of it.

Per the neighbouring ``tests/sn/mesh/`` convention the whole module is
``foundation`` and NOTHING carries a ``verifies(...)`` edge (the R12a
equation node is claimed by the classification battery, not by admission).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_mu
from orpheus.numerics.quadrature.rules_circle import (
    NODE_ALIGNED,
    STAGGERED,
    periodic_trapezoid,
)
from orpheus.numerics.quadrature.rules_product import spherical_product
from orpheus.numerics.quadrature.rules_sphere import LevelStructure
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.sweep.pole_angular_closure import (
    MarchStart,
    assert_carrying_quadrature,
    non_carrying_levels,
)
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

# The four message fragments of the admission contract (gate_design §2.2).
# Short by design — tests pin substrings, so the grep-sensitive fragment
# must be the SHORTEST distinctive one (coding-standards).
_F_NEW = "non-carrying"
_F_EDGE = "starts on an ordinate"
_F_TIE = "shared eta-minimum"
_F_FIX = "folded_product"
_F_STRUCTURELESS = "level structure"


def _cyl_mesh(nx: int = 4) -> Mesh1D:
    return Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_right=BC("reflective"),
    )


def _build(quad: Quadrature) -> SNMesh:
    return SNMesh(_cyl_mesh(), quad, placeholder_materials())


def _full_product(n_mu: int, n_phi: int, shift) -> Quadrature:
    """A FULL-circle product rule at a chosen azimuthal shift, unquotiented."""
    measure, structure = spherical_product(
        gauss_legendre_on_mu(n_mu),
        periodic_trapezoid(n_phi, shift=shift),
    )
    return Quadrature(measure=measure, level_structure=structure)


# ── positive rows (vv anti-#11: admission asserted, not only refusal) ──


@pytest.mark.parametrize("n_mu,n_phi", [(4, 8), (2, 4), (6, 12)])
def test_the_folded_product_family_constructs(n_mu, n_phi):
    """A-P1 — the carrying folded family is ADMITTED at three independent
    (n_mu, n_phi) points; the mesh derives the full multi-level carrier."""
    sn = _build(Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi))
    if sn.radial_characteristic_levels != tuple(range(n_mu)):
        pytest.fail(
            f"folded({n_mu},{n_phi}) admitted but carrier "
            f"{sn.radial_characteristic_levels} != every level position"
        )


def test_a_hand_built_sigma_y_quotient_constructs_without_the_factory_tag():
    """A-P2 — admission by STRUCTURE, not provenance: a rule bit-equal to
    ``folded_product(4,8)`` but hand-assembled with ``folded_by is None``
    is ADMITTED.  A gate keyed on the factory tag (``folded_by``) instead
    of the march-start facts reds here (mutation MA7)."""
    ref = Quadrature.folded_product(n_mu=4, n_phi=8)
    hand = Quadrature(
        measure=DiscreteMeasure(
            nodes=np.array(ref.measure.nodes, copy=True),
            weights=np.array(ref.measure.weights, copy=True),
            support="S^2",
        ),
        level_structure=LevelStructure(
            n_levels=ref.level_structure.n_levels,
            level_indices=[
                np.array(li, copy=True)
                for li in ref.level_structure.level_indices
            ],
            level_mu=np.array(ref.level_structure.level_mu, copy=True),
            polar_invariant=ref.level_structure.polar_invariant,
            azimuth=np.array(ref.level_structure.azimuth, copy=True),
            hemisphere=np.array(ref.level_structure.hemisphere, copy=True),
        ),
    )
    # The control's precondition, asserted explicitly (lessons L2): the
    # hand rule must genuinely lack the tag, or this row silently
    # degenerates into a duplicate of A-P1.
    if hand.folded_by is not None:
        pytest.fail("the hand-built rule carries a folded_by tag — the "
                    "structure-vs-provenance control is invalid.")
    np.testing.assert_array_equal(hand.measure.nodes, ref.measure.nodes)
    np.testing.assert_array_equal(hand.measure.weights, ref.measure.weights)
    sn = _build(hand)
    if sn.radial_characteristic_levels != (0, 1, 2, 3):
        pytest.fail(f"hand-built quotient admitted but carrier "
                    f"{sn.radial_characteristic_levels} != (0, 1, 2, 3)")


def test_the_spherical_arm_is_untouched():
    """A-P3 — the guard is CYLINDRICAL-scoped: the open-GL sphere (carrying
    by the same predicate) constructs exactly as before the flip."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials())
    if sn.radial_characteristic_levels != (0,):
        pytest.fail("sphere carrier moved — the admission guard leaked "
                    "outside the CYLINDRICAL arm.")


# ── negative rows (each asserts the fired fact AND the absent one) ─────


@pytest.mark.parametrize("n_mu,n_phi", [(4, 8), (2, 4)])
def test_a_full_node_aligned_product_refuses_on_an_edge_node_start(n_mu, n_phi):
    """A-N1 — the full NODE_ALIGNED product: a node ON Σ per level
    (`[M]` on_edge_node=True, degenerate=False — edge-node ONLY), so the
    refusal quotes the edge reason and NOT the tie reason."""
    with pytest.raises(ValueError) as exc:
        _build(Quadrature.product(n_mu=n_mu, n_phi=n_phi))
    msg = str(exc.value)
    assert _F_NEW in msg and _F_EDGE in msg and _F_FIX in msg, msg
    assert _F_TIE not in msg, (
        f"wrong-reason refusal: the NODE_ALIGNED product's message quotes "
        f"the tie reason — {msg}"
    )


def test_a_full_staggered_product_refuses_on_the_mirror_tie():
    """A-N2 — the UNFOLDED staggered product: mirror partners tie the
    η-minimum bit-exactly (`[M]` degenerate=True, on_edge_node=False), so
    the refusal quotes the tie reason and NOT the edge reason."""
    with pytest.raises(ValueError) as exc:
        _build(_full_product(4, 8, STAGGERED))
    msg = str(exc.value)
    assert _F_NEW in msg and _F_TIE in msg, msg
    assert _F_EDGE not in msg, (
        f"wrong-reason refusal: the staggered product's message quotes "
        f"the edge reason — {msg}"
    )


@pytest.mark.parametrize("order", [4, 18])
def test_level_symmetric_refuses_at_every_order_on_the_hemisphere_tie(order):
    """A-N3 — level-symmetric rules: hemisphere sign replication ties the
    η-minimum on every level at every order (`[M]` all four probed orders
    produce the identical ``(False, True)`` facts, so extra rows would
    exercise ONE mechanism — S2/S8 are deliberate ride-alongs, not
    overlooked).  S4 is the suite's common fixture; S18 is the family's
    #337-pinned TOP of range — if the defined range moves again, this row
    keeps the refusal battery aligned with the family's actual extent."""
    with pytest.raises(ValueError) as exc:
        _build(Quadrature.level_symmetric(order))
    msg = str(exc.value)
    assert _F_NEW in msg and _F_TIE in msg, msg
    assert _F_EDGE not in msg, msg


def test_a_folded_node_aligned_quotient_refuses_despite_carrying_the_fold_tag():
    """A-N4 — the provenance pincer's other jaw: a σ_y QUOTIENT of the
    NODE_ALIGNED product carries the factory tag yet its arcs start ON an
    ordinate — REFUSED.  With A-P2 this pins admission-by-structure from
    both sides (mutation MA7 reds exactly this pair)."""
    quad = _full_product(4, 8, NODE_ALIGNED).quotient(SubgroupOfO3.Mirror("y"))
    if quad.folded_by is None:
        pytest.fail("the quotient lost its folded_by tag — the pincer's "
                    "tagged-but-non-carrying jaw is invalid.")
    with pytest.raises(ValueError) as exc:
        _build(quad)
    msg = str(exc.value)
    assert _F_NEW in msg and _F_EDGE in msg, msg
    assert _F_TIE not in msg, msg


# ── discrimination rows — which guard owns which input ─────────────────


def test_the_structureless_guard_still_owns_the_slab_and_sphere_cubatures():
    """A-D1 — the wiring-order pin: GL / Lebedev refuse with the OLD
    ``level structure`` message and never see the new guard (the carrying
    helper sits AFTER ``cylindrical_streaming``; `[M]` the classifier
    itself does not raise on GL, so a helper wired one line earlier would
    silently take ownership and change the message).  Also asserts the
    fragment DISJOINTNESS both ways, against A-N1's captured message.

    The two existing gates this row does NOT twin (they assert "GL raises
    with `level structure`"; this row asserts "…and nothing else claims
    it"): ``tests/sn/sweep/core/test_sweep_regression.py`` and
    ``tests/sn/sweep/curvilinear/test_cyl_sweep_regression.py``.
    """
    for quad in (Quadrature.gauss_legendre(4), Quadrature.lebedev(17)):
        with pytest.raises(ValueError) as exc:
            _build(quad)
        msg = str(exc.value)
        assert _F_STRUCTURELESS in msg, msg
        assert _F_NEW not in msg, (
            f"the NEW refusal fired for a structure-less cubature — the "
            f"carrying guard was wired before cylindrical_streaming: {msg}"
        )
    # Fragment disjointness, asserted once in each direction.
    with pytest.raises(ValueError) as exc_new:
        _build(Quadrature.product(n_mu=4, n_phi=8))
    assert _F_STRUCTURELESS not in str(exc_new.value), (
        "the new admission message contains 'level structure' — the two "
        "existing structure-less gates would start matching it."
    )


def test_the_trace_rank_guard_owns_the_all_tangential_folded_rule():
    """A-D2 (gate_design §8.1) — ``folded_product(n, 2)`` folds to ONE
    pure-azimuthal angle per level (every μ_x == 0.0 bit-exact) and is
    refused by the TRACE tier's ``build_omega_dot_n`` — an EARLIER guard
    than the carrying predicate (which `[M]` reports carrying=True for
    this rule).  Admission is the tier guards' CONJUNCTION in firing
    order; this row pins the order for the all-tangential class."""
    with pytest.raises(ValueError) as exc:
        _build(Quadrature.folded_product(n_mu=2, n_phi=2))
    msg = str(exc.value)
    assert "mu_x" in msg, msg
    assert _F_NEW not in msg, (
        f"the carrying refusal fired for the all-tangential rule — the "
        f"trace-rank guard no longer owns it: {msg}"
    )


# ── invariant rows on the pure predicate (the ∀-quantifier's teeth) ────


def _synthetic_mixed() -> tuple[MarchStart, ...]:
    """(carrying, edge-node, carrying, degenerate) — `[M]` unreachable
    through every shipped factory (a mixed rule cannot be built), which
    is exactly why the pure function must be gated directly."""
    return (
        MarchStart(on_edge_node=False, degenerate=False),
        MarchStart(on_edge_node=True, degenerate=False),
        MarchStart(on_edge_node=False, degenerate=False),
        MarchStart(on_edge_node=False, degenerate=True),
    )


def test_a_mixed_rule_reports_only_the_offending_levels(monkeypatch):
    """A-Q1 — the ∀ quantifier: a mixed tuple yields exactly the offender
    positions (1, 3); and the raiser's message names the FIRST offender
    (level 1) with the fact that fired there (edge) and NOT the fact that
    fired elsewhere (the tie at position 3)."""
    assert non_carrying_levels(_synthetic_mixed()) == (1, 3)

    from orpheus.sn.sweep import pole_angular_closure as pac

    monkeypatch.setattr(
        pac, "march_start_structure_per_level",
        lambda quad, coord: _synthetic_mixed(),
    )
    with pytest.raises(ValueError) as exc:
        pac.assert_carrying_quadrature(object(), CoordSystem.CYLINDRICAL)
    msg = str(exc.value)
    assert "level 1" in msg, msg
    assert _F_EDGE in msg, msg
    assert _F_TIE not in msg, (
        f"the message quotes the tie reason, which fired at position 3, "
        f"not at the reported level 1: {msg}"
    )


def test_an_all_carrying_tuple_has_no_offenders():
    """A-Q2 — the predicate's positive fixed point: all-carrying → ()."""
    starts = tuple(
        MarchStart(on_edge_node=False, degenerate=False) for _ in range(5)
    )
    assert non_carrying_levels(starts) == ()


def test_the_offender_set_is_positions_not_a_bool():
    """A-Q3 (vv anti-#14) — the function returns the STRUCTURE the
    docstring names: a tuple of int positions, never a bool about them."""
    out = non_carrying_levels(_synthetic_mixed())
    assert isinstance(out, tuple)
    assert all(isinstance(p, int) for p in out)
    assert not isinstance(out, bool)

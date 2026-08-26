r"""The R12a march-start facts — the integer re-pose of the seed predicate (Q5.4/T26).

Two structural facts per μ-level (:class:`MarchStart`), each a bit-exact
identity on the level's realization, replace the raw first-ordinate M-M
float that used to conflate them:

* ``on_edge_node`` — the march-start edge IS an ordinate (η-minimum node
  on Σ; the ω = π most-inward direction).
* ``degenerate``  — the η-minimum is shared (a double-cover tie: mirror
  partners on a full circle, hemisphere partners on level-symmetric
  rules), killing the recurrence's (1−τ₀) thread weight.

Carrying ⟺ neither. The former predicate — ``τ_raw,0 ∈ (0,1)``
exclusive — is demoted to a THEOREM about the closure's edge arithmetic
and gated here bit-exactly per family: ``on_edge_node ⟹ τ_raw,0 == 0.0``,
``degenerate ⟹ τ_raw,0 == 1.0``, neither ⟹ strict interior, with NO
epsilon anywhere.

T26's measured hazard — odd ``n_φ`` flipping the classification on a
5.6e-16 trig-round-off gap — was an artifact of the pre-E3
``linspace``+cos azimuths. `[M]` post-E3 (roots of unity) the trichotomy
is bit-exact at every shipped ``n_φ`` including 5 and 7; the theorem
gate pins that exactness, so a regression to sloppy azimuths reds here.

The σ_y-FOLDED product rule is the forward-looking row: its levels are
arcs (T22b), the start is genuinely off-node, and the classification is
CARRYING — pinned now, consumed when the fold reaches the solve (Q5.6).
"""

from __future__ import annotations

import pytest

from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import NODE_ALIGNED, STAGGERED, Quadrature
from orpheus.sn.angular.closure import (
    march_start_structure_per_level,
    morel_montry_tau_per_level,
)
from tests.sn._test_helpers import seam_quad


# (label, quad builder, coord, expected (on_edge_node, degenerate))
_FAMILIES = [
    (
        "product_node_aligned_even_4",
        lambda: Quadrature.product(4, 4),
        CoordSystem.CYLINDRICAL,
        (True, False),
    ),
    (
        "product_node_aligned_even_8",
        lambda: Quadrature.product(4, 8),
        CoordSystem.CYLINDRICAL,
        (True, False),
    ),
    # T26's flip case: odd n_phi. Pre-E3 the tie was 5.6e-16 off exact
    # and the float predicate flipped to CARRYING; post-E3 the mirror
    # eta-tie is bit-exact and the fact is DEGENERATE, stably.
    (
        "product_node_aligned_odd_5",
        lambda: Quadrature.product(4, 5),
        CoordSystem.CYLINDRICAL,
        (False, True),
    ),
    (
        "product_node_aligned_odd_7",
        lambda: Quadrature.product(4, 7),
        CoordSystem.CYLINDRICAL,
        (False, True),
    ),
    # A FULL staggered rule is a double cover marched as one sequence:
    # the two most-negative-eta ordinates are mirror partners.
    (
        "product_staggered_full",
        lambda: seam_quad(4, 8, STAGGERED, folded=False),
        CoordSystem.CYLINDRICAL,
        (False, True),
    ),
    (
        "level_symmetric_4",
        lambda: Quadrature.level_symmetric(4),
        CoordSystem.CYLINDRICAL,
        (False, True),
    ),
    (
        "level_symmetric_8",
        lambda: Quadrature.level_symmetric(8),
        CoordSystem.CYLINDRICAL,
        (False, True),
    ),
    (
        "sphere_gl_8",
        lambda: Quadrature.gauss_legendre(8),
        CoordSystem.SPHERICAL,
        (False, False),
    ),
    # The fold's future: a folded staggered level is an ARC (T22b) with
    # a genuine off-node start — the cylinder starts carrying ψ½ state.
    (
        "folded_staggered",
        lambda: seam_quad(4, 8, STAGGERED, folded=True),
        CoordSystem.CYLINDRICAL,
        (False, False),
    ),
    # A folded NODE_ALIGNED rule keeps its Σ endpoints: start ON a node.
    (
        "folded_node_aligned",
        lambda: seam_quad(4, 8, NODE_ALIGNED, folded=True),
        CoordSystem.CYLINDRICAL,
        (True, False),
    ),
]


@pytest.mark.foundation
@pytest.mark.verifies("sn-direct-seed-r12a-predicate")
@pytest.mark.parametrize(
    ("label", "build", "coord", "expected"),
    _FAMILIES,
    ids=[f[0] for f in _FAMILIES],
)
def test_the_two_facts_classify_every_family(
    label: str, build, coord, expected: tuple[bool, bool]
) -> None:
    """Every level of every family reads the same (on_edge_node, degenerate)."""
    starts = march_start_structure_per_level(build(), coord)
    for p, start in enumerate(starts):
        assert (start.on_edge_node, start.degenerate) == expected, (
            f"{label} level {p}: got "
            f"({start.on_edge_node}, {start.degenerate}), expected {expected}"
        )
        carrying = expected == (False, False)
        assert start.consumes_independent_seed == carrying


@pytest.mark.foundation
@pytest.mark.parametrize(
    ("label", "build", "coord", "expected"),
    _FAMILIES,
    ids=[f[0] for f in _FAMILIES],
)
def test_the_tau_trichotomy_is_a_theorem_about_the_facts(
    label: str, build, coord, expected: tuple[bool, bool]
) -> None:
    r"""The facts' bit-exact τ consequence, where a τ exists at all.

    ``on_edge_node ⟹ τ_0 == 0.0`` (the start edge and the first node
    coincide bit-for-bit, so the numerator is exactly zero);
    ``degenerate ⟹ τ_0 == 1.0``; neither ⟹ strictly interior.  This
    gate is ALSO the post-E3 exactness pin: under the pre-E3
    linspace+cos azimuths the odd-n_φ rows read 0.9999999999999994
    (T26) and this comparison fails bit-exactly.

    ⛔ **NARROWED at Q5.6.4 (2026-08-11), and the narrowing is the
    point.** The trichotomy above was a theorem about the retired
    **η-midpoint** cell partition — in particular ``degenerate ⟹
    τ_0 == 1`` held *because* ``eta_edge[m+1] = (η_m + η_{m+1})/2``
    collapses to ``η_0`` when the two coincide.  The cylinder partition
    is now the midpoint in **ω**
    (:func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`),
    and a level that is not a monotone ω-arc has **no angular cells at
    all** — so for those families there is no τ to make a claim about,
    and the producer REFUSES.  That is asserted here instead of the
    trichotomy, which keeps the family coverage rather than deleting
    the rows.

    `[M]` Which arm each family lands in, on this tree:

    * **sphere GL** — cumulative-weight edges, ``μ_{1/2} = −1``, so
      ``on_edge_node`` (``μ_0 == −1``) still gives ``τ_0 == 0``
      bit-exactly.  The trichotomy is LIVE here.
    * **folded arcs** — every level is carrying (``Σ = ∅``), so only the
      strict-interior arm is reachable.  ⚠ And it is reachable
      *necessarily*, not contingently: on a monotone arc the ω-midpoint
      edges bracket their own node, so ``τ ∈ (0,1)`` is forced (see
      ``test_tau_arc_wellposedness.py::test_the_cylinder_cannot_violate_P3_once_its_arc_is_MONOTONE``).
    * **full-circle products / level_symmetric** — REFUSED at the
      partition.  These are also inadmissible for a cylindrical
      ``SNMesh`` since the 6.3 flip, so no production path loses a
      claim.
    """
    quadlike = build()
    starts = march_start_structure_per_level(quadlike, coord)
    try:
        tau_levels = morel_montry_tau_per_level(quadlike, coord)
    except ValueError as exc:
        # No angular cells ⟹ no τ ⟹ no trichotomy. Pin the refusal and
        # its REASON, so this row still fails if the producer silently
        # starts accepting a non-arc level again.
        assert "not a monotone arc in omega" in str(exc), (
            f"{label}: the producer refused, but not for the reason this "
            f"row pins (expected a non-monotone-arc refusal): {exc}"
        )
        assert coord is CoordSystem.CYLINDRICAL, (
            f"{label}: only a CYLINDRICAL level can fail to be an arc; a "
            f"spherical rule has one level and no ω. Got {coord!r}: {exc}"
        )
        return
    for p, (start, tau_level) in enumerate(
        zip(starts, tau_levels, strict=True)
    ):
        t0 = float(tau_level[0])
        if start.on_edge_node:
            assert t0 == 0.0, f"{label} level {p}: expected exact 0, got {t0!r}"
        elif start.degenerate:
            assert t0 == 1.0, f"{label} level {p}: expected exact 1, got {t0!r}"
        else:
            assert 0.0 < t0 < 1.0, (
                f"{label} level {p}: expected strict interior, got {t0!r}"
            )


@pytest.mark.foundation
def test_cartesian_is_refused() -> None:
    """Cartesian has no angular march — asking for its start is an error."""
    with pytest.raises(ValueError, match="no angular march"):
        march_start_structure_per_level(
            Quadrature.gauss_legendre(4), CoordSystem.CARTESIAN
        )

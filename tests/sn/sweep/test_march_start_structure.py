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
from orpheus.sn.sweep.pole_angular_closure import (
    march_start_structure_per_level,
    morel_montry_tau_raw_per_level,
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
    """The old encoding is the facts' bit-exact consequence — no epsilon.

    ``on_edge_node ⟹ τ_raw,0 == 0.0`` (the start edge and the first
    node coincide bit-for-bit, so the numerator is exactly zero);
    ``degenerate ⟹ τ_raw,0 == 1.0`` (η₀ == η₁ makes the midpoint edge
    exactly η₀, so numerator == denominator); neither ⟹ strictly
    interior. This gate is ALSO the post-E3 exactness pin: under the
    pre-E3 linspace+cos azimuths the odd-n_φ rows read
    0.9999999999999994 (T26) and this comparison fails bit-exactly.
    """
    quadlike = build()
    starts = march_start_structure_per_level(quadlike, coord)
    raw = morel_montry_tau_raw_per_level(quadlike, coord)
    for p, (start, tau_level) in enumerate(zip(starts, raw, strict=True)):
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

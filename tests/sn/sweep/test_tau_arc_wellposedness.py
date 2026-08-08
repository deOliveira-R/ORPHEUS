r"""τ_raw's [0,1] membership RAISES; the fold bounds τ away from the singularity (Q5.5/T27).

T27 adjudicated the cylinder's fused ``max(0.5, min(1.0, τ_raw))`` as TWO
objects welded together:

* the ``[0, 1]`` MEMBERSHIP — promoted (Q5.5) to a guard that RAISES at
  the producer: ``τ_raw ∉ [0, 1]`` means an ordinate sits outside its own
  angular cell, impossible on a monotone arc, so membership certifies the
  level's march order.  T22's ω-ordered mis-ordering produced
  ``τ_raw = 1.079`` and the absorber silently laundered it into a finite
  wrong answer 400 lines downstream; the guard stops it at source.
  ⚠ It does NOT catch the double cover — the full-circle
  ``[0, 1, 0, 1, …]`` fingerprint is entirely INSIDE ``[0, 1]``; that
  detector is the singular set Σ / the fold criterion (T24).
* the ``[½, 1]`` ABSORPTION — stays until the fold wiring (Q5.6),
  because production NODE_ALIGNED cylinders still start on an edge node
  (``τ_raw,0 = 0`` bit-exact; the recurrence divides by zero unclamped).
  On the folded arc its reason is structurally GONE, which is what the
  mechanism gates below pin.

THE MECHANISM GATE (T27's spec): at the most-activating configuration —
folded + staggered, swept to the largest supported ``n_φ`` (min τ_raw
falls monotonically toward 1/5, so large ``n_φ`` is the worst case;
τ_raw is ``n_μ``-independent, the sinθ cancels in the ratio) — assert
the MECHANISM (``Σ = ∅``, computed via ``singular_set``, never declared)
and the CONSEQUENCE (``τ_raw ⊂ [1/5, 4/5]``, strictly away from
``{0, 1}``, plus the reversal identity ``τ_m + τ_{M−1−m} = 1``).
Reddenable on ONE mutation: revert Q5.2's staggered offset to δ = 0
(patch ``_FOLD_SHIFT`` to ``NODE_ALIGNED``) and BOTH legs fail
together — Σ gains its fixed nodes AND τ_raw hits ``{0, 1}``.  That is
what attributes the pass to the mechanism rather than to luck.

`[M]` 2026-08-07 (the Q5.5 probe, n_μ = 4): folded staggered τ_raw spans
``[0.2195, 0.7805]`` at n_φ = 8 falling to ``[0.200289, 0.799711]`` at
n_φ = 64, ``|Σ| = 0`` throughout, reversal residual **0.0 bit-exact** at
every n_φ — the selection-descent fold (Q5.3: charts bit-copied, order
by restriction) upgraded the identity from the ~1e-15 of T27's
hand-built pre-Q5.3 probe to exact, so the gate asserts it with NO
epsilon.  Under δ = 0: ``|Σ| = 8``, τ_raw = ``[0.0, 1.0]``.

Every row is a claim about quadrature/closure algebra, not a transport
equation — the module is ``foundation``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import (
    NODE_ALIGNED,
    STAGGERED,
    Quadrature,
    gauss_legendre_on_mu,
    periodic_trapezoid,
    spherical_product,
)
from orpheus.numerics.symmetry import SubgroupOfO3, singular_set
from orpheus.sn.sweep.pole_angular_closure import morel_montry_tau_raw_per_level
from tests.sn._test_helpers import seam_quad

pytestmark = pytest.mark.foundation

# The T27 mechanism gates build their parent through this shift — the
# fold's design offset (Q5.2).  The teeth mutation ("revert Q5.2's
# offset to δ = 0") patches ``orpheus.numerics.quadrature.STAGGERED``
# to NODE_ALIGNED before collection, so this binding inherits it; the
# mutation must red BOTH mechanism legs together.
_FOLD_SHIFT = STAGGERED

_MIRROR_Y = SubgroupOfO3.Mirror("y")

# n_φ sweep for the mechanism gates: min τ_raw falls monotonically
# toward 1/5, so the largest supported size is the most-activating row.
_N_PHI_SWEEP = (8, 16, 32, 64)


# ── the [0, 1] membership guard (vv #11: positive + negative) ───────────

def test_shipped_families_pass_the_guard_including_the_closed_endpoints():
    """Positive: production and folded rules stay inside [0, 1] — no raise.

    The NODE_ALIGNED full circle is the sharp row: its double-cover
    ``[0, 1, 0, 1, …]`` τ pattern ATTAINS both closed endpoints
    (``0`` = edge-node start, ``1`` = η-degenerate tie — legal march
    starts per Q5.4), so this row proves the guard admits the boundary
    rather than merely the interior.
    """
    rows = [
        ("node_aligned_full_circle", lambda: Quadrature.product(4, 8),
         CoordSystem.CYLINDRICAL),
        ("folded_staggered", lambda: seam_quad(4, 8, STAGGERED, folded=True),
         CoordSystem.CYLINDRICAL),
        ("sphere_gl", lambda: Quadrature.gauss_legendre(8),
         CoordSystem.SPHERICAL),
    ]
    for label, build, coord in rows:
        raw = morel_montry_tau_raw_per_level(build(), coord)  # must not raise
        values = np.concatenate(raw)
        assert np.all((values >= 0.0) & (values <= 1.0)), (
            f"{label}: the guard passed a value outside [0, 1]: {values}"
        )
        if label == "node_aligned_full_circle":
            assert 0.0 in values and 1.0 in values, (
                f"{label}: expected BOTH closed endpoints attained "
                f"(the double-cover fingerprint), got {values}"
            )


def test_a_mis_ordered_level_is_refused_as_outside_its_cell():
    """Negative: ω-ordered members on the full circle → τ_raw = 1.079 → raise.

    The ω-chart orders the FULL circle's level non-monotonically in η
    (η = sinθ·cosω wraps), so one ordinate lands outside its own angular
    cell — T22's exact mis-ordering, which the retired absorber used to
    launder into a finite wrong answer.
    """
    measure, structure = spherical_product(
        gauss_legendre_on_mu(4), periodic_trapezoid(8, shift=NODE_ALIGNED)
    )
    mis_levels = []
    for members in structure.level_indices:
        members = np.asarray(members)
        mis_levels.append(members[np.argsort(structure.azimuth[members])])
    quad = SimpleNamespace(
        mu_x=measure.nodes[:, 0],
        mu_y=measure.nodes[:, 1],
        mu_z=measure.nodes[:, 2],
        weights=measure.weights,
        level_indices=tuple(mis_levels),
    )
    with pytest.raises(ValueError, match="outside its own angular cell"):
        morel_montry_tau_raw_per_level(quad, CoordSystem.CYLINDRICAL)


# ── the T27 mechanism gates (two legs, one mutation reds both) ──────────

@pytest.mark.parametrize("n_phi", _N_PHI_SWEEP)
def test_the_fold_mechanism_is_an_empty_singular_set(n_phi: int):
    """MECHANISM leg: the staggered parent has Σ = ∅ — computed, not declared.

    Σ = ∅ means σ_y fixes NO node, every orbit has size exactly 2, and
    the fold is a clean 2-to-1 selection whose arc has no endpoint on
    the mirror — which is WHY τ_raw stays away from the singularity
    (the consequence leg).  Under δ = 0 the parent regains its fixed
    nodes at φ ∈ {0, π} and this leg reds with |Σ| = 2·n_levels.
    """
    parent, _ = spherical_product(
        gauss_legendre_on_mu(4), periodic_trapezoid(n_phi, shift=_FOLD_SHIFT)
    )
    sigma = singular_set(parent, _MIRROR_Y)
    assert sigma.size == 0, (
        f"n_phi={n_phi}: the fold's mechanism is broken — the parent has "
        f"{sigma.size} singular node(s) under {_MIRROR_Y.name}; the arc "
        f"would start ON the mirror and τ_raw would reach the {{0, 1}} "
        f"singularity"
    )


@pytest.mark.verifies("morel-montry-folded-arc")
@pytest.mark.parametrize("n_phi", _N_PHI_SWEEP)
def test_the_folded_tau_is_bounded_with_the_reversal_identity(n_phi: int):
    """CONSEQUENCE leg: τ_raw ⊂ [1/5, 4/5] and τ_m + τ_{M−1−m} = 1 bit-exactly.

    With midpoint nodes on the arc the smallest-η node sits at
    ω = π − ε, ε = π/2n, giving τ_0 = (ε²/2)/(5ε²/2) = 1/5 in the limit,
    approached monotonically FROM INSIDE (`[M]` 0.2929 at n_φ = 4 →
    0.200289 at n_φ = 64) — so the closed [1/5, 4/5] box holds at every
    finite n_φ and the {0, 1} singularity is structurally unreachable.
    The reversal identity is bit-exact because the staggered
    roots-of-unity η-grid is bit-antisymmetric under ω → π − ω and the
    fold descends by bit-copied selection (Q5.3); the retired absorber
    would DESTROY it (clamping τ_0 to ½ breaks τ_0 + τ_{M−1} = 1).
    """
    quad = seam_quad(4, n_phi, _FOLD_SHIFT, folded=True)
    raw = morel_montry_tau_raw_per_level(quad, CoordSystem.CYLINDRICAL)
    for p, tau in enumerate(raw):
        assert np.all(tau >= 0.2) and np.all(tau <= 0.8), (
            f"n_phi={n_phi} level {p}: τ_raw left [1/5, 4/5]: {tau} — the "
            f"fold no longer bounds the closure away from the singularity"
        )
        np.testing.assert_array_equal(
            tau + tau[::-1],
            np.ones_like(tau),
            err_msg=(
                f"n_phi={n_phi} level {p}: the reversal identity "
                f"τ_m + τ_{{M−1−m}} = 1 is no longer bit-exact — the "
                f"arc lost its reversal symmetry (or the fold stopped "
                f"descending by bit-copied selection)"
            ),
        )

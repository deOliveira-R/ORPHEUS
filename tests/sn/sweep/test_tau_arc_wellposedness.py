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
from orpheus.sn.sweep.pole_angular_closure import morel_montry_tau_per_level
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

def _arc_quad(omega_per_level, sin_theta=1.0):
    """A hand-built CYLINDRICAL level from bare ω values (one level).

    Lets a row choose the arc's node placement directly, which is the
    only way to reach P3's closed endpoints now that no shipped rule
    attains them (see the two rows below).
    """
    om = np.asarray(omega_per_level, dtype=float)
    mu_z = float(np.sqrt(max(0.0, 1.0 - sin_theta ** 2)))
    return SimpleNamespace(
        mu_x=sin_theta * np.cos(om),
        mu_y=sin_theta * np.sin(om),
        mu_z=np.full(om.size, mu_z),
        weights=np.full(om.size, np.pi / om.size),
        # ⚠ an ARRAY, never a tuple: ``mu_x[(0, 1, 2, 3)]`` is numpy
        # MULTI-DIMENSIONAL indexing, not fancy indexing, and raises
        # "too many indices". Real level_indices are arrays.
        level_indices=(np.arange(om.size),),
    )


def test_shipped_families_pass_the_guard_including_the_closed_endpoints():
    r"""Positive: shipped rules stay inside [0, 1] — no raise.

    ⛔ **Re-posed at Q5.6.4 (2026-08-11).** The sharp row used to be
    ``node_aligned_full_circle`` (``Quadrature.product(4, 8)``), whose
    double-cover ``[0, 1, 0, 1, …]`` τ pattern attained BOTH closed
    endpoints.  That rule is now **refused one frame earlier**, by the
    arc guard in
    :func:`~orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level`
    — a full-circle level has ω of both signs and therefore no angular
    cells at all.  Its refusal is gated in
    ``test_mms_ordering_blindness.py::test_the_full_circle_double_cover_is_REFUSED_by_the_cell_partition``.

    So the closed-endpoint claim moves to a hand-built **node-aligned
    ARC**: put a node exactly at :math:`\omega = \pi` (i.e. ON
    :math:`\Sigma`) and :math:`\tau_0 = 0` bit-exactly, because the
    partition's lower endpoint IS :math:`\omega = \pi`.  ``0`` therefore
    remains a legal, ATTAINABLE value that the guard must admit — the
    ``on_edge_node`` march start of Q5.4.
    """
    rows = [
        ("folded_staggered", lambda: seam_quad(4, 8, STAGGERED, folded=True),
         CoordSystem.CYLINDRICAL),
        ("sphere_gl", lambda: Quadrature.gauss_legendre(8),
         CoordSystem.SPHERICAL),
        ("node_aligned_ARC (a node ON Sigma)",
         lambda: _arc_quad(np.pi - np.arange(4) * (np.pi / 4)),
         CoordSystem.CYLINDRICAL),
    ]
    attained = {}
    for label, build, coord in rows:
        tau = morel_montry_tau_per_level(build(), coord)  # must not raise
        values = np.concatenate(tau)
        attained[label] = values
        assert np.all((values >= 0.0) & (values <= 1.0)), (
            f"{label}: the guard passed a value outside [0, 1]: {values}"
        )
    endpoint_row = attained["node_aligned_ARC (a node ON Sigma)"]
    assert endpoint_row[0] == 0.0, (
        f"the closed endpoint 0 is no longer attainable: a node placed "
        f"exactly on Sigma (omega = pi) must give tau_0 = 0 bit-exactly, "
        f"got {endpoint_row}"
    )


def test_the_cylinder_cannot_violate_P3_once_its_arc_is_MONOTONE():
    r"""⭐ P3 is a THEOREM on the cylinder, and the sphere is where it still bites.

    ⛔ **Re-posed at Q5.6.4 (2026-08-11).** This row used to be
    ``test_a_mis_ordered_level_is_refused_as_outside_its_cell``: it fed
    ω-ordered members of a FULL circle to the producer and asserted the
    P3 message ``"outside its own angular cell"`` (T22's mis-ordering,
    :math:`\tau = 1.079`).  The refusal still happens — one frame
    earlier and better targeted, from the arc guard — so the old regex
    no longer matches.

    ⚠ **And the reason is structural, not cosmetic.** With the cell
    partition taken as the **ω-midpoint**, a strictly ω-monotone level
    has ``edge[m] > ω_m > edge[m+1]`` by construction, and ``cos`` is
    monotone on :math:`(0, \pi)`, so
    ``η_edge[m] < η_m < η_edge[m+1]`` — i.e. :math:`\tau \in (0,1)` is
    **forced**.  `[M]` verified over 4000 random monotone arcs:
    ``min(τ) = 4.739e-07``, ``min(1−τ) = 7.599e-10``, never outside.
    ⟹ **once the arc guard passes, cylinder-P3 CANNOT fire**, and its
    only equality case is a node ON an arc endpoint — i.e. a node on
    :math:`\Sigma`.  **Cylinder-P3 reduces to the fold criterion
    Σ = ∅.**

    So P3 keeps its teeth on the **SPHERE**, where cumulative-weight
    edges genuinely do NOT bracket their nodes: a 3-D ``level_symmetric``
    rule on the 1-D spherical arm measured
    :math:`\tau \in [-20.3, 1.13]`, 23 of 24 outside (#336).  That is
    the row below.  Stated plainly so no audit reads cylinder-P3 as
    live coverage it is not (`vv-principles` #22).
    """
    # (a) cylinder — the mis-ordering is refused, by the ARC guard.
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
    with pytest.raises(ValueError, match="monotone arc in omega"):
        morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)

    # (b) sphere — P3 itself, still falsifiable, still the only catcher.
    bad_sphere = SimpleNamespace(
        mu_x=np.array([-0.9, -0.1, 0.1, 0.9]),
        weights=np.array([0.1, 0.1, 0.1, 0.1]),   # sum 0.4 != 2 -> cells tiny
        N=4,
    )
    with pytest.raises(ValueError, match="outside its own angular cell"):
        morel_montry_tau_per_level(bad_sphere, CoordSystem.SPHERICAL)


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
    r"""CONSEQUENCE leg: τ ⊂ [1/4, 3/4], and τ_m + τ_{M−1−m} = 1 to ULP.

    ⛔ **Both numbers were RE-POSED at Q5.6.4 (2026-08-11)**, and the
    originals are kept because they are still the correct claims about
    the partition they were measured on. They read: *"τ_raw ⊂ [1/5, 4/5]
    … τ_0 = (ε²/2)/(5ε²/2) = 1/5 in the limit, `[M]` 0.2929 at n_φ = 4 →
    0.200289 at n_φ = 64"*, with the reversal identity **bit-exact**.

    Those were properties of the retired **η-midpoint (chord)**
    partition. On the ω-midpoint partition the closed form is
    :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\tan(\Delta\omega/4)`,
    so the most-inward node (:math:`\omega_0 = \pi - \Delta\omega/2`)
    gives

    > :math:`\tau_0 = \tfrac12 - \tfrac12\cot(\Delta\omega/2)
    > \tan(\Delta\omega/4) \;\longrightarrow\; \tfrac14`

    since :math:`\cot(\Delta\omega/2)\tan(\Delta\omega/4) \to \tfrac12`.
    ⟹ the box is :math:`[\tfrac14, \tfrac34]`, again approached
    monotonically FROM INSIDE, so :math:`\{0, 1\}` stays structurally
    unreachable. `[M]` measured:

    ======= ======================== ==================
    n_φ     τ range                  reversal residual
    ======= ======================== ==================
    4       [0.292893, 0.707107]     0.5 ULP
    8       [0.259892, 0.740108]     0.5 ULP
    16      [0.252425, 0.747575]     2.0 ULP
    32      [0.250603, 0.749397]     7.0 ULP
    64      [0.250151, 0.749849]    12.0 ULP
    ======= ======================== ==================

    ⚠ **The reversal identity is no longer BIT-exact, and that is a
    trade, not a regression.** The chord partition's reversal symmetry
    was exact *because* both end cells were stretched symmetrically —
    the 17.5 %-ω-width defect cancelled itself under
    :math:`\omega \to \pi - \omega`. The ω partition has the correct
    cells and pays 0.5–12 ULP from :math:`\arctan2`/:math:`\cos`
    round-off. Asserted at 64 ULP: the residual grows like the arc
    refinement, so a bit-exact assertion here would be a latent false
    red (`vv-principles` #16 — never assert tighter than the producer
    achieves).
    """
    quad = seam_quad(4, n_phi, _FOLD_SHIFT, folded=True)
    tau_levels = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)
    for p, tau in enumerate(tau_levels):
        assert np.all(tau >= 0.25) and np.all(tau <= 0.75), (
            f"n_phi={n_phi} level {p}: τ left [1/4, 3/4]: {tau} — the fold "
            f"no longer bounds the closure away from the {{0, 1}} "
            f"singularity (the limit tau_0 -> 1/4 is approached from "
            f"INSIDE, so any value below 1/4 means the arc changed)"
        )
        np.testing.assert_array_almost_equal_nulp(
            tau + tau[::-1], np.ones_like(tau), nulp=64
        )

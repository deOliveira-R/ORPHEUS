r"""The ANGULAR CELL PARTITION producer's value gate (:eq:`angular-cell-partition`).

:func:`~orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level`
is the **single producer** of a :math:`\mu`-level's angular cell
partition — the object every coefficient that references "the boundary
between cell :math:`m` and :math:`m+1`" must read.  Until this module
landed (2026-08-11) **nothing in the tree asserted what it returns**, on
either arm: the Q5.6.4 carve minted it as the one-source partition and
the only gates downstream read :math:`\tau`, which is P2 *applied* to the
partition.  A wrong partition that happens to keep :math:`\tau` in
:math:`[0,1]` was therefore invisible to every committed gate except the
two cylinder :math:`\tau` rows in
``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py``.

**References, and both are authored HERE** (`vv-principles` L11 —
structural independence; the production body must share no code path
with them):

* **SPHERE** — the cumulative-weight ladder :math:`\mu_{1/2} = -1`,
  :math:`\mu_{m+1/2} = \mu_{m-1/2} + w_m`
  (:cite:`BaileyMorelChang2010` Eq. 12; Lathrop 2000 p. 249's
  :math:`\sum \Delta\mu_m = 2` is the closing identity).  Written with
  ``np.cumsum``, whose pairwise reduction is a *different* association
  from the producer's sequential loop — which is what makes it a second
  computation rather than a transcription.
* **CYLINDER** — the analytic equispaced-arc closed form
  :math:`e_k = \sin\theta\,\cos(\pi - k\,\Delta\omega)`,
  :math:`\Delta\omega = \pi/M`, derived from the arc geometry with no
  reference to the production expression (which takes midpoints of the
  *stored* :math:`\omega` values, correct for any monotone arc).

**Negative controls** (`vv-principles` #19 — a positive reading cannot
show a gate is loaded on the structure it is credited with):

* sphere → the naive **uniform**-:math:`\Delta\mu` partition, which the
  cumulative-weight rows must DIFFER from (`[M]` by 0.15–0.21);
* cylinder → the **retired chord** (:math:`\eta`-midpoint) partition,
  which the ω-midpoint rows must DIFFER from at :math:`M \ge 3`.

⚠ **The M = 2 trap, gated explicitly below.**  ``folded_product(·, 4)``
gives :math:`M = 2` per level, where the chord and ω-midpoint partitions
coincide to 0 ULP in edge space — **no** :math:`n_\varphi = 4` fixture
can see a partition change.  Every cylinder row here therefore runs at
:math:`n_\varphi \ge 6`, and both parities of :math:`M` are covered
(:math:`n_\varphi = 6/10/18 \Rightarrow M = 3/5/9` odd;
:math:`8/16/32/64 \Rightarrow M = 4/8/16/32` even) because the arc's
middle edge behaves differently at odd :math:`M` — at odd :math:`M` an
ordinate sits exactly at :math:`\omega = \pi/2`, i.e. :math:`\eta = 0`,
which is the equality case of the orientation law below.

Every row is a claim about quadrature/closure algebra evaluated
solve-free — the module is ``foundation``.  Mode-8 discipline: every
value assertion is ``np.testing`` / ``pytest.fail`` so it fires under the
canonical ``python -O`` invocation.

Mutation-proof — measured 2026-08-11
====================================

13 in-process mutations (source-transformed from ``inspect.getsource``,
rebound across every ``sys.modules`` holder) run over a 298-row scope:
this module, ``test_tau_arc_wellposedness``, ``test_march_start_structure``,
``tests/sn/sweep/curvilinear/``, ``test_mms_ordering_blindness`` and
``test_cylinder{ical}_quadrature_admission``.  Un-mutated control:
298 passed / 0 failed.  Rows reddened per gate (``.`` = green):

======================================================  == == == == === == == == == == ==
gate                                                    MC M1 M2 M3 M3b M4 M5 M6 M7 M8 M9
======================================================  == == == == === == == == == == ==
``sphere…cumulative_weight_ladder``                      7  8  8  7   .  .  .  .  .  .  .
``sphere…spans_exactly_minus_one_to_plus_one``           .  8  8  .   .  .  .  .  .  .  .
``sphere…is_NOT_the_uniform_partition``                  4  .  .  4   .  .  .  .  .  .  .
``cylinder…equispaced_arc_closed_form``                  7  .  .  .   .  7  7  7  7  .  .
``cylinder…spans_the_level_arc_exactly``                 .  .  .  .   .  .  7  7  .  .  .
``cylinder…is_NOT_the_retired_chord_convention``         .  .  .  .   .  7  .  .  .  .  .
``the_M2_fold_is_BLIND_to_the_partition_choice``         .  .  .  .   .  .  1  1  1  .  .
``the_march_orientation_sign_is_pinned``                 6  6  6  3   .  .  6  6  6 12 12
``a_bare_level_namespace…``                              1  .  .  .   .  1  .  1  1  .  .
======================================================  == == == == === == == == == == ==

* **MC** POSITIVE CONTROL — the producer returns a uniform partition on
  both arms: 75 rows red across 9 files (46 here + in the two sibling
  modules, 29 pre-existing), so an all-blind verdict would have been a
  broken instrument (`vv-principles` #17).
* **M1** sphere seed ``-1 → 0``; **M2** ``+ w[n] → + w[n]/2``;
  **M3** accumulation replaced by a uniform step; **M4** the retired
  chord partition; **M5** ``sinθ`` dropped; **M6** the arc endpoints
  swapped (:math:`\omega_{1/2} = 0`); **M7** interior edge at
  :math:`\tfrac34\omega_m + \tfrac14\omega_{m+1}` instead of the
  midpoint; **M8** :math:`\tau \to 1-\tau`; **M9** :math:`\tau \equiv
  \tfrac12`.
* ⛔ **M3b is a REFUTED mutation candidate**: reversing the sphere's
  weight accumulation order (``w[n] → w[N-1-n]``) is **bit-identically
  invisible** — Gauss–Legendre weights are palindromic, so it reddened
  **0 of 298** rows.  No gate in this tree can catch it and none should
  claim to; it is a Mode-12 stabiliser of the rule's own symmetry.
* ⛔ **What the committed suite could NOT see before this module.** For
  the two mutations the Q5.6.4 partition carve is exposed to, the
  pre-existing catchers were **6 rows each** (of 298): M4 and M7 are
  caught only by ``test_tau_arc_wellposedness``'s four
  :math:`[\tfrac14,\tfrac34]` box rows plus the two committed cylinder
  :math:`\tau` rows.  For the ORIENTATION flip M8 the box rows are
  **provably blind** (see the orientation gate's docstring) and the
  pre-existing catchers are 6 :math:`\tau`-value rows.  Nothing asserted
  the partition itself.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.sweep.pole_angular_closure import (
    angular_cell_edges_per_level,
    morel_montry_tau_per_level,
)

pytestmark = pytest.mark.foundation

_EPS = float(np.finfo(float).eps)

# Sphere orders.  Includes ODD N (3, 5) — legal Gauss–Legendre rules that
# put a node exactly at μ = 0, the equality case of the orientation law.
_SPHERE_N = (2, 3, 4, 5, 8, 16, 32, 64)

# Cylinder parent-circle sizes.  M = n_phi/2 per level, so this covers
# M = 3, 4, 5, 8, 9, 16, 32 — BOTH parities.  n_phi = 4 (M = 2) is
# deliberately absent: it is the degenerate row, gated on its own below.
_CYL_N_PHI = (6, 8, 10, 16, 18, 32, 64)


# ── references, authored here, sharing no code path with production ────


def _sphere_reference_edges(quad) -> np.ndarray:
    r"""BMC Eq. 12 by hand: :math:`\mu_{1/2} = -1`, then cumulative weight.

    ``np.cumsum`` reduces pairwise; the producer accumulates in a
    sequential Python loop.  The two associations differ, so agreement is
    to FP-association noise rather than bit-exact — that difference IS
    the independence.
    """
    w = np.asarray(quad.weights, dtype=float)
    return np.concatenate(([-1.0], -1.0 + np.cumsum(w)))


def _sphere_edge_atol(N: int) -> float:
    r"""Pairwise-vs-sequential summation bound for :math:`N` terms on
    :math:`[-1, 1]`: :math:`\lesssim 2 N \varepsilon` per edge, doubled
    for margin.

    `[M]` 2026-08-11, measured ``max|producer − cumsum|`` (Gauss–Legendre):

    ====  ==========  ==============
    N     max gap     gap / (N·eps)
    ====  ==========  ==============
    2     ``0.0``     0.00
    8     ``2.2e-16`` 0.12
    16    ``2.2e-16`` 0.06
    32    ``1.1e-16`` 0.02
    64    ``4.4e-16`` 0.03
    128   ``8.9e-16`` 0.03
    ====  ==========  ==============

    so the bound carries ≥ 30× margin at every measured order while still
    widening correctly with ``N`` (a fixed ``atol`` would become a false
    red eventually).
    """
    return 4.0 * N * _EPS


def _arc_reference_edges(sin_theta: float, M: int) -> np.ndarray:
    r"""The analytic equispaced-arc partition, hand-derived.

    The level's arc runs :math:`\omega = \pi \to 0` in :math:`M` equal
    steps :math:`\Delta\omega = \pi/M`, so the cell edges sit at
    :math:`\omega_k = \pi - k\,\Delta\omega` and their radial cosine is
    :math:`\sin\theta\cos\omega_k`.  No reference to the producer's
    stored-:math:`\omega` midpoint expression.
    """
    return sin_theta * np.cos(np.pi - np.arange(M + 1) * (np.pi / M))


def _chord_edges(eta: np.ndarray, sin_theta: float) -> np.ndarray:
    r"""The RETIRED :math:`\eta`-midpoint (chord) partition — the
    NEGATIVE control.

    Interior edges are the midpoint of consecutive *radial cosines*
    rather than of the march variable :math:`\omega`; the endpoints stay
    pinned at :math:`\mp\sin\theta`, so the outermost cells stretch.
    This is the convention Q5.6.4 retracted, and it is
    **bit-indistinguishable from the survivor at** :math:`M = 2`.
    """
    M = eta.size
    edges = np.empty(M + 1)
    edges[0], edges[M] = -sin_theta, sin_theta
    edges[1:M] = 0.5 * (eta[:-1] + eta[1:])
    return edges


def _level_geometry(quad, level_idx) -> tuple[np.ndarray, float, int]:
    """``(eta, sin_theta, M)`` for one cylinder level."""
    eta = np.asarray(quad.mu_x)[level_idx]
    sin_theta = float(np.sqrt(1.0 - np.asarray(quad.mu_z)[level_idx[0]] ** 2))
    return eta, sin_theta, len(level_idx)


# ═══════════════════════════════════════════════════════════════════════
# SPHERE — the cumulative-weight ladder
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("angular-cell-partition")
@pytest.mark.parametrize("N", _SPHERE_N)
def test_sphere_partition_is_the_cumulative_weight_ladder(N: int):
    r"""Sphere edges == a hand-written cumulative-weight reference.

    **Pins**: the sphere branch of :eq:`angular-cell-partition` —
    :math:`\mu_{1/2} = -1` and :math:`\mu_{m+1/2} = \mu_{m-1/2} + w_m`,
    i.e. *a cell's width IS its quadrature weight*.  A seed drift, a
    missing/duplicated weight, a halved weight, a reversed accumulation,
    or a switch to any other partition (uniform, node-midpoint) moves an
    edge by O(0.1–1) and reds here.

    **Cannot catch**: (a) anything about the *choice* of the
    cumulative-weight convention — that is BMC Eq. 12 verbatim and is
    argued, not tested, in
    :func:`~orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level`;
    the row below (`…_is_NOT_the_uniform_partition`) is the only
    convention-loading evidence here.  (b) A mutation that permutes the
    *weights* palindromically — Gauss–Legendre weights are symmetric, so
    ``w[n] → w[N-1-n]`` is bit-identically invisible to this and to every
    other gate in the tree (a Mode-12 stabiliser of the rule's own
    symmetry, not a hole in this row).  (c) Anything on the cylinder arm.
    """
    quad = Quadrature.gauss_legendre(N)
    (edges,) = angular_cell_edges_per_level(quad, CoordSystem.SPHERICAL)
    reference = _sphere_reference_edges(quad)

    np.testing.assert_allclose(
        edges, reference, rtol=0.0, atol=_sphere_edge_atol(N),
        err_msg=(
            f"sphere N={N}: the partition producer departs from the "
            f"hand-written BMC Eq. 12 cumulative-weight ladder\n"
            f"producer ={edges}\nreference={reference}"
        ),
    )


@pytest.mark.verifies("angular-cell-partition")
@pytest.mark.parametrize("N", _SPHERE_N)
def test_sphere_partition_spans_exactly_minus_one_to_plus_one(N: int):
    r"""Sphere: :math:`M+1` edges, strictly increasing, from :math:`-1`
    bit-exactly, closing at :math:`+1`.

    **Pins** the structural contract the return type promises AND the
    CLOSING identity :math:`\mu_{M+1/2} = +1` — equivalently Lathrop
    2000 p. 249's :math:`\sum_m \Delta\mu_m = 2`, which is the statement
    that the cells *tile* :math:`[-1, 1]` with no gap and no overlap.
    A weight-scaling error (the classic missing/spurious factor 2, or a
    half-range rule fed to the full-range arm) leaves the ladder
    monotone and reds only HERE.

    **Cannot catch**: an error that redistributes width between cells
    while preserving the total — the closing identity is a telescoping
    sum (`vv-principles` #8), so per-cell errors that cancel are
    invisible to this row.  The value row above is what sees those.

    ⚠ The closing edge is NOT bit-exactly ``1.0``: `[M]` the residual is
    ``0.0`` at N ∈ {2,4,8,32,64} but ``+2.2e-16`` at N ∈ {12,16} and
    ``+4.4e-16`` at N = 48 (the GL weights themselves sum to
    ``2.0000000000000004`` at N ∈ {10, 64}).  Asserted at ``1e-14``;
    ``==`` would be a latent false red (`vv-principles` #16).
    """
    quad = Quadrature.gauss_legendre(N)
    (edges,) = angular_cell_edges_per_level(quad, CoordSystem.SPHERICAL)

    if edges.size != N + 1:
        pytest.fail(f"sphere N={N}: expected {N + 1} edges, got {edges.size}")
    if edges[0] != -1.0:
        pytest.fail(
            f"sphere N={N}: the march START edge must be mu = -1 "
            f"bit-exactly (the Carlson starting direction); got {edges[0]!r}"
        )
    gaps = np.diff(edges)
    if not np.all(gaps > 0.0):
        pytest.fail(
            f"sphere N={N}: the partition is not strictly increasing — "
            f"cell widths {gaps} (a cell's width IS its weight, so a "
            f"non-positive width means a non-positive quadrature weight "
            f"or a mis-ordered rule)"
        )
    np.testing.assert_allclose(
        edges[-1], 1.0, rtol=0.0, atol=1e-14,
        err_msg=(
            f"sphere N={N}: the partition does not close at mu = +1 "
            f"(sum of cell widths = {edges[-1] + 1.0!r}, must be 2 — "
            f"Lathrop 2000 p. 249); the cells do not tile [-1, 1]"
        ),
    )


@pytest.mark.parametrize("N", (4, 8, 16, 32))
def test_sphere_partition_is_NOT_the_uniform_partition(N: int):
    r"""NEGATIVE control: the sphere partition is weight-cumulative, NOT
    equal-width.

    `vv-principles` #19 — without this leg the value row above would pass
    just as happily against ``np.linspace(-1, 1, N+1)``, i.e. it would
    certify the *arithmetic* while carrying no information about the
    *convention*.  The uniform partition is the plausible-wrong
    alternative a fresh implementer reaches for ("M cells over
    :math:`[-1,1]`"), and it satisfies every structural property the row
    above checks: right length, strictly increasing, endpoints
    :math:`\mp 1`, closing identity exact.

    `[M]` 2026-08-11 max separation: 0.1521 (N=4) / 0.1764 (8) / 0.1908
    (16) / 0.2020 (32) — asserted at 0.1, ≥ 1.5× margin.

    **Cannot catch** any error smaller than the convention gap; this is a
    discriminator, not a value gate.
    """
    quad = Quadrature.gauss_legendre(N)
    (edges,) = angular_cell_edges_per_level(quad, CoordSystem.SPHERICAL)
    uniform = np.linspace(-1.0, 1.0, N + 1)
    gap = float(np.max(np.abs(edges - uniform)))
    if gap <= 0.1:
        pytest.fail(
            f"sphere N={N}: the partition is indistinguishable from the "
            f"UNIFORM equal-width partition (max separation {gap:.3e}); "
            f"this row and the cumulative-weight value row would then be "
            f"certifying nothing about the convention"
        )


# ═══════════════════════════════════════════════════════════════════════
# CYLINDER — the ω-midpoint arc partition
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies("angular-cell-partition")
@pytest.mark.parametrize("n_phi", _CYL_N_PHI)
def test_cylinder_partition_is_the_equispaced_arc_closed_form(n_phi: int):
    r"""Cylinder edges == :math:`\sin\theta\cos(\pi - k\,\Delta\omega)`.

    **Pins** the cylinder branch of :eq:`angular-cell-partition`: the
    boundary is the midpoint **in** :math:`\omega`, the variable the
    azimuthal march marches in (T22b), and its radial cosine is
    :math:`\sin\theta\cos\omega_{m\pm1/2}`.  Reds on: the retired chord
    partition, a dropped :math:`\sin\theta`, a wrong arc endpoint, a
    left/right-node edge instead of the midpoint, an :math:`\eta`-space
    cumulative-weight transplant of the sphere convention.

    **Cannot catch**: (a) a producer that is wrong only on a
    NON-equispaced arc — the reference is the equispaced closed form, and
    every shipped rule is equispaced in :math:`\varphi`, so the
    "midpoint of the *stored* :math:`\omega`" generalisation is exercised
    but not *constrained* here (`vv-principles` Mode 10); (b) the
    convention choice itself — the chord control below is that evidence;
    (c) anything at :math:`M = 2`, which is excluded from the parameter
    set (see the M = 2 control).

    ⚠ `[M]` 2026-08-11 the agreement is **M-independent** at ≤ 1.5 ULP —
    ``0.0`` / ``2.2e-16`` / ``0.0`` / ``0.0`` / ``1.1e-16`` / ``1.9e-16``
    / ``2.2e-16`` / ``3.3e-16`` / ``3.3e-16`` / ``3.3e-16`` at
    :math:`n_\varphi = 4/6/8/10/16/18/32/64/128/256` — because both sides
    are one ``cos`` of an O(1) angle.  (Contrast :math:`\tau`, which
    divides by a cell width and inherits :math:`\cot`'s condition number,
    degrading like :math:`M\varepsilon`.)  Asserted at ``8·eps``.
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    edges = angular_cell_edges_per_level(quad, CoordSystem.CYLINDRICAL)

    for p, level_idx in enumerate(quad.level_indices):
        _eta, sin_theta, M = _level_geometry(quad, level_idx)
        reference = _arc_reference_edges(sin_theta, M)
        np.testing.assert_allclose(
            edges[p], reference, rtol=0.0, atol=8.0 * _EPS,
            err_msg=(
                f"cylinder n_phi={n_phi} level {p} (M={M}, "
                f"sin(theta)={sin_theta!r}): the partition producer "
                f"departs from the analytic equispaced-arc closed form "
                f"e_k = sin(theta)*cos(pi - k*dw)\n"
                f"producer ={edges[p]}\nreference={reference}"
            ),
        )


@pytest.mark.verifies("angular-cell-partition")
@pytest.mark.parametrize("n_phi", _CYL_N_PHI)
def test_cylinder_partition_spans_the_level_arc_exactly(n_phi: int):
    r"""Cylinder: :math:`M+1` edges, strictly increasing, endpoints
    :math:`\mp\sin\theta` **bit-exactly**.

    **Pins** the closing identity of the arc: the level's cells tile the
    whole half-circle :math:`\omega \in [0, \pi]`, whose radial cosine
    runs :math:`-\sin\theta \to +\sin\theta`.  The endpoints are asserted
    with ``==``: `[M]` 2026-08-11 the residual is **exactly 0.0** at
    every :math:`n_\varphi \in \{4,\dots,256\}` because the producer
    writes :math:`\omega = \pi` and :math:`\omega = 0` as literals and
    ``cos`` is exact on both.  That is a bit-exactness the producer
    genuinely achieves, so asserting it is not over-tight
    (`vv-principles` #16) — and it is what makes a dropped
    :math:`\sin\theta`, or a first/last cell shifted by half a step,
    red HERE rather than only in the value row.

    **Cannot catch**: interior edge placement (telescoping, `vv` #8) —
    the value row above owns that.
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    edges = angular_cell_edges_per_level(quad, CoordSystem.CYLINDRICAL)

    if len(edges) != len(quad.level_indices):
        pytest.fail(
            f"n_phi={n_phi}: expected one edge array per mu-level "
            f"({len(quad.level_indices)}), got {len(edges)}"
        )
    for p, level_idx in enumerate(quad.level_indices):
        _eta, sin_theta, M = _level_geometry(quad, level_idx)
        if edges[p].size != M + 1:
            pytest.fail(
                f"n_phi={n_phi} level {p}: expected {M + 1} edges for "
                f"{M} ordinates, got {edges[p].size}"
            )
        if edges[p][0] != -sin_theta or edges[p][-1] != sin_theta:
            pytest.fail(
                f"n_phi={n_phi} level {p}: the arc endpoints must be "
                f"-/+ sin(theta) = -/+{sin_theta!r} bit-exactly (omega = "
                f"pi and 0); got {edges[p][0]!r} and {edges[p][-1]!r}"
            )
        gaps = np.diff(edges[p])
        if not np.all(gaps > 0.0):
            pytest.fail(
                f"n_phi={n_phi} level {p}: the partition is not strictly "
                f"increasing in the radial cosine — widths {gaps}; the "
                f"march is not a monotone arc"
            )


@pytest.mark.parametrize("n_phi", _CYL_N_PHI)
def test_cylinder_partition_is_NOT_the_retired_chord_convention(n_phi: int):
    r"""NEGATIVE control: the edges are ω-midpoints, NOT η-midpoints.

    `vv-principles` #19.  The chord (:math:`\eta`-midpoint) partition is
    the convention this producer used until 2026-08-11, and it satisfies
    every structural property the rows above check — right length,
    strictly increasing, endpoints :math:`\mp\sin\theta` exactly.  Its
    interior edges are exactly :math:`\cos(\Delta\omega/2)\times` the
    survivor's, so the separation is

    .. math::  \text{gap} \;\approx\;
       \sin\theta\,\bigl(1 - \cos(\Delta\omega/2)\bigr),

    which is what this row asserts (at 0.4× for margin) rather than a
    magic constant — the gap SHRINKS like :math:`M^{-2}`, so any fixed
    floor becomes a false red at high :math:`n_\varphi`.

    `[M]` 2026-08-11, measured gap ÷ predicted floor, level 0:
    0.500 (M=3) / 0.707 (4) / 0.809 (5) / 0.924 (8) / 0.940 (9) /
    0.981 (16) / 0.995 (32) / 0.999 (64) — approaching 1 from below, so
    0.4 carries ≥ 1.25× margin at the tightest row and ≥ 2× everywhere
    else.

    **Cannot catch** anything except this one convention; it is a
    discriminator.  In particular it is **structurally blind at**
    :math:`M = 2` — see the control row below, which is why
    :math:`n_\varphi = 4` is not in the parameter set.
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    edges = angular_cell_edges_per_level(quad, CoordSystem.CYLINDRICAL)

    for p, level_idx in enumerate(quad.level_indices):
        eta, sin_theta, M = _level_geometry(quad, level_idx)
        chord = _chord_edges(eta, sin_theta)
        gap = float(np.max(np.abs(edges[p] - chord)))
        floor = 0.4 * sin_theta * (1.0 - np.cos(0.5 * np.pi / M))
        if gap <= floor:
            pytest.fail(
                f"cylinder n_phi={n_phi} level {p} (M={M}): the shipped "
                f"partition is indistinguishable from the RETIRED chord "
                f"(eta-midpoint) convention — max separation {gap:.4e} "
                f"<= the predicted floor {floor:.4e}. Every value row in "
                f"this module would then pass equally against the "
                f"partition the Q5.6.4 carve replaced."
            )


def test_the_M2_fold_is_BLIND_to_the_partition_choice():
    r"""CONTROL: at :math:`M = 2` the chord and ω-midpoint partitions
    coincide to round-off — no ``folded_product(·, 4)`` fixture can see
    this carve.

    ``folded_product(n_mu, 4)`` gives :math:`M = 2` ordinates per level,
    whose single interior edge is the level's midpoint under *both*
    conventions (:math:`\omega_{1/2} = \pi/2 \Rightarrow \eta = 0`, and
    the chord midpoint of :math:`\pm\eta_0` is also 0).

    ⚠ `[M]` 2026-08-11 the agreement is **not bit-exact but is 1 ULP**:
    the ω-midpoint edge is ``sin(theta)*cos(pi/2)`` =
    ``3.11e-17``/``5.76e-17`` (``np.cos(np.pi/2)`` is ``6.12e-17``, not
    zero) while the chord's is exactly ``0.0``; in :math:`\tau` space the
    two agree to ``1.11e-16``.  Asserted at ``1e-15`` — the point is that
    the separation is **15 orders below** the ``1e-1``-scale gap the
    partition choice produces at :math:`M \ge 3`, not that the bits
    match.

    This row exists so that no future audit reads an :math:`n_\varphi=4`
    green as coverage of the partition choice (`vv-principles` #20 — a
    row that structurally cannot see the varied thing is not coverage),
    and so that the blindness is a MEASURED, committed fact rather than
    a comment.  It reds if the degeneracy ever stops holding — which
    would mean one of the two conventions moved.

    ⚠ ``folded_product(·, 4)`` is the fixture ~15 cylinder gates in this
    tree use.  Read this row before crediting any of them with partition
    coverage.
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=4)
    edges = angular_cell_edges_per_level(quad, CoordSystem.CYLINDRICAL)

    for p, level_idx in enumerate(quad.level_indices):
        eta, sin_theta, M = _level_geometry(quad, level_idx)
        if M != 2:
            pytest.fail(
                f"folded_product(4, 4) level {p} no longer has M = 2 "
                f"(got {M}); this control's whole subject is the M = 2 "
                f"degeneracy — re-derive it"
            )
        np.testing.assert_allclose(
            edges[p], _chord_edges(eta, sin_theta), rtol=0.0, atol=1e-15,
            err_msg=(
                f"level {p}: the M = 2 degeneracy no longer holds — the "
                f"omega-midpoint and chord partitions now DIFFER by more "
                f"than round-off. One of the two conventions moved; the "
                f"'n_phi = 4 is blind to the partition' statement carried "
                f"by this module and by test_tau_producer_equivalence.py "
                f"must be re-measured."
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# BOTH ARMS — the march-orientation sign
# ═══════════════════════════════════════════════════════════════════════


# ⚠ Parametrized by LABEL only, never by a pre-computed ``(cosines, tau)``
# tuple.  A ``parametrize`` argument list is evaluated at module IMPORT,
# i.e. at COLLECTION — so calling the producer there turns any production
# defect that RAISES (the P3 guard is one line downstream of the
# partition) into ``Interrupted: 1 error during collection``, which
# reports **0 failures** for the whole run.  Measured 2026-08-11 while
# mutation-testing this very module: six of thirteen mutations came back
# as a flattering "0 caught" until the build moved into the body
# (`vv-principles` #17 — the harness lies before the code does; the
# sibling shape is test-architect lesson L41e).
_ORIENTATION_LABELS = (
    "sphere_GL3", "sphere_GL4", "sphere_GL5",
    "sphere_GL8", "sphere_GL16", "sphere_GL32",
    "cyl_folded_4x6", "cyl_folded_4x8", "cyl_folded_4x10",
    "cyl_folded_4x16", "cyl_folded_4x18", "cyl_folded_4x32",
)


def _orientation_case(label: str):
    """``(cosines per level, tau per level)`` — built INSIDE the test."""
    if label.startswith("sphere_GL"):
        quad = Quadrature.gauss_legendre(int(label.removeprefix("sphere_GL")))
        return (
            (np.asarray(quad.mu_x),),
            morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL),
        )
    quad = Quadrature.folded_product(
        n_mu=4, n_phi=int(label.rsplit("x", 1)[1])
    )
    return (
        tuple(np.asarray(quad.mu_x)[i] for i in quad.level_indices),
        morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL),
    )


@pytest.mark.parametrize("label", _ORIENTATION_LABELS)
def test_the_march_orientation_sign_is_pinned(label):
    r"""⭐ The ordinate sits DOWNSTREAM of its cell centre iff its radial
    cosine is positive: :math:`(\tau_m - \tfrac12)\,\mu_m \ge 0`,
    bit-exactly, on both arms.

    **This row exists because the march ORIENTATION is otherwise
    ungated.**  The partition is traversed from the most-inward direction
    (:math:`\mu = -1` sphere; :math:`\omega = \pi` cylinder) outward, and
    :math:`\tau_m = (\mu_m - \mu_{m-1/2})/(\mu_{m+1/2} - \mu_{m-1/2})`
    measures the ordinate's barycentric position **from the upstream
    edge**.  Reverse the march, or measure from the downstream edge, and
    every :math:`\tau` becomes :math:`1 - \tau` — a change that is
    **invisible** to the three properties the committed suite asserts,
    because all three are symmetric under :math:`\tau \to 1-\tau`:

    * the P3 membership guard :math:`\tau \in [0,1]` — symmetric;
    * the fold's box :math:`\tau \in [\tfrac14, \tfrac34]`
      (``test_tau_arc_wellposedness.py``) — symmetric;
    * the reversal identity :math:`\tau_m + \tau_{M-1-m} = 1` (same
      file) — invariant, since :math:`(1-\tau_m) + (1-\tau_{M-1-m}) = 1`.

    So this is a `vv-principles` Mode-12 hole in the committed set, and
    the law here is its catcher.

    On the cylinder the closed form is
    :math:`\tau_m - \tfrac12 = \tfrac12\cot\omega_m\tan(\Delta\omega/4)`
    and :math:`\cot\omega_m` has the sign of :math:`\eta_m`; on the
    sphere the same ordering holds for Gauss–Legendre because its nodes
    and weights are symmetric about :math:`\mu = 0`.

    **The equality case is exact and is the point of covering both
    parities**: at odd :math:`M` (cylinder :math:`n_\varphi = 6/10/18`)
    and odd :math:`N` (sphere GL3/GL5) an ordinate sits exactly at
    :math:`\mu = 0`, where `[M]` :math:`|\tau - \tfrac12| \le 1.11e-16`
    (0.5 ULP).  Both legs are asserted: the product is
    :math:`\ge 0` **with no tolerance** (``-0.0 >= 0.0`` is True in
    IEEE-754, so the zero cosine passes cleanly), and at a zero cosine
    :math:`\tau` must equal :math:`\tfrac12` to 4 ULP.

    **Cannot catch**: (a) a *magnitude* error in :math:`\tau` that keeps
    the sign — the value gates in
    ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` own
    that; (b) :math:`\tau \equiv \tfrac12` (the plain angular diamond),
    which makes every product exactly zero — caught here only by the
    ACTIVATION leg below, which is why that leg is committed rather than
    left as a comment (`vv-principles` Mode 8, tautological-guard class).
    """
    cosines, taus = _orientation_case(label)
    activated = False
    for p, (mu, tau) in enumerate(zip(cosines, taus, strict=True)):
        offset = tau - 0.5
        product = offset * mu
        if not np.all(product >= 0.0):
            bad = int(np.argmin(product))
            pytest.fail(
                f"{label} level {p}: the march orientation is inverted at "
                f"ordinate {bad} — (tau - 1/2)*mu = {product[bad]!r} < 0 "
                f"(tau={tau[bad]!r}, mu={mu[bad]!r}).\n"
                f"tau measures the ordinate's position from the UPSTREAM "
                f"cell edge along an inward-to-outward march, so tau > 1/2 "
                f"exactly where the radial cosine is positive. A value "
                f"below 1/2 there means tau was measured from the "
                f"downstream edge (tau -> 1 - tau) or the level's march "
                f"order was reversed — neither of which the [0,1] guard, "
                f"the [1/4,3/4] box, nor the reversal identity can see.\n"
                f"tau    ={tau}\nmu     ={mu}\nproduct={product}"
            )
        on_axis = mu == 0.0
        if np.any(on_axis):
            np.testing.assert_allclose(
                tau[on_axis], 0.5, rtol=0.0, atol=4.0 * _EPS,
                err_msg=(
                    f"{label} level {p}: an ordinate lying exactly on the "
                    f"radial-cosine zero (omega = pi/2) must be the "
                    f"barycentre of its own cell, tau = 1/2; got "
                    f"{tau[on_axis]}"
                ),
            )
        activated = activated or bool(np.any(product > 0.0))

    if not activated:
        pytest.fail(
            f"{label}: EVERY (tau - 1/2)*mu product is zero, so the "
            f"inequality above is satisfied vacuously. Either tau has "
            f"collapsed to the constant 1/2 (the plain angular diamond — "
            f"NOT the Morel-Montry closure this module gates) or the "
            f"radial cosines are all zero."
        )


# ═══════════════════════════════════════════════════════════════════════
# CARTESIAN — no angular march, hence no partition
# ═══════════════════════════════════════════════════════════════════════


def test_cartesian_has_no_angular_cell_partition():
    """Cartesian is refused by the producer, naming the reason.

    The slab has no angular redistribution, so there is no march and no
    cell partition; the closure supplies the neutral :math:`\\tau = 1`
    by its TYPE (``IdentityAngularClosure``), never by a coord branch
    inside this producer.

    **Cannot catch** a wrong *value* on either curvilinear arm; it pins
    only that the third coordinate system is refused rather than
    silently given the sphere's or the cylinder's convention.
    """
    quad = Quadrature.gauss_legendre(8)
    with pytest.raises(ValueError, match="no angular cell partition"):
        angular_cell_edges_per_level(quad, CoordSystem.CARTESIAN)


def test_a_bare_level_namespace_is_enough_to_reach_the_producer():
    r"""Foundation: the producer reads only ``mu_x``/``mu_y``/``mu_z``/
    ``level_indices`` — no ``Quadrature`` type dependency.

    Guards the hand-built-arc idiom the neighbouring modules rely on
    (``test_tau_arc_wellposedness.py::_arc_quad``), and pins the
    ⚠ ``level_indices`` trap in one place: the entries must be numpy
    ARRAYS.  A tuple of ints triggers numpy MULTI-dimensional indexing
    (``mu_x[(0,1,2)]`` is ``mu_x[0,1,2]``) and raises ``IndexError:
    too many indices`` several frames from the cause.

    **Cannot catch** any numerical claim — it is a structural/API row.
    """
    M = 5
    omega = np.pi - (np.arange(M) + 0.5) * (np.pi / M)
    bare = SimpleNamespace(
        mu_x=np.cos(omega),
        mu_y=np.sin(omega),
        mu_z=np.zeros(M),
        weights=np.full(M, np.pi / M),
        level_indices=(np.arange(M),),   # an ARRAY, never a tuple of ints
    )
    (edges,) = angular_cell_edges_per_level(bare, CoordSystem.CYLINDRICAL)
    np.testing.assert_allclose(
        edges, _arc_reference_edges(1.0, M), rtol=0.0, atol=8.0 * _EPS,
        err_msg=(
            "a hand-built equispaced arc at sin(theta) = 1 must reproduce "
            "the analytic partition; if this reds, the hand-built-arc "
            "idiom used by the sibling gates no longer reaches the same "
            "code path production does"
        ),
    )

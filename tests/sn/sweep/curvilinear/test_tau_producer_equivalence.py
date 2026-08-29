r"""Issue #236 Phase 2 — τ producer-equivalence gate (Leg 1).

The :class:`~orpheus.sn.angular.closure.MorelMontryAngularSweep`
closure now PRODUCES the Morel–Montry angular weight τ from the quadrature
``(μ, w, levels)`` it already binds (an angular-scheme property), instead of
reading it back from the streaming-GEOMETRY factory (``reduced_operator.py``).

Step C (the geometry-τ retirement) deleted the geometry-side τ producer, so
the former ``*_equals_geometry_factory_0ulp`` legs (closure τ vs
``spherical_streaming(...).tau_mm`` / ``cylindrical_streaming(...).tau_mm_per_level``)
are now vacuous and have been removed.  The independence floor is carried
entirely by the surviving legs below, and **the two arms now use two
DIFFERENT references** — a distinction that must not be flattened again:

⛔ ``angular_differencing.morel_montry_weights`` **DELEGATED to the
production producer** (Q5.6.4) — deliberately, so a "reference" can never
drift into a second definition of the angular cell, which is exactly how
its cylinder arm had gone wrong.  It was therefore **never again an
independent reference for τ**, and it was RETIRED outright on
2026-08-12 (its body was an ``orpheus.sn`` import, an edge
``tests/test_layer_imports.py`` forbids from ``derivations/``).  Both
arms below use
hand-authored references instead:

* **SPHERE** → a hand-written cumulative-weight expression, inline.  BMC
  2010 Eq. 12/42: ``mu_edge[n+1] = mu_edge[n] + w[n]`` from ``-1``, then
  P2.  Authored here, independent of the producer.
* **CYLINDER** → the **analytic closed form**
  :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\tan(\Delta\omega/4)`,
  hand-derived from the arc geometry.  STRUCTURAL independence — it shares
  no code path with the producer, a strictly stronger footing than the
  procedural twin it replaces.

This gate proves the closure-produced τ:

* equals a structurally- (cylinder) or procedurally- (sphere) INDEPENDENT
  reference;
* on the cylinder, DIFFERS from the retired chord convention — a negative
  control (`vv-principles` #19), without which the row would pass equally
  against the partition the Q5.6.4 carve replaced and so could not be
  evidence about the partition choice;
* on Cartesian, is the neutral τ = 1.0 produced WITHOUT a geometry branch
  (the closure TYPE is the dispatch — Pattern 4).

Mode-8 discipline: every numeric assertion is ``np.testing`` /
``np.array_equal`` / ``pytest.fail`` so it FIRES under the canonical ``-O``
invocation.  The two guard-the-guard predicates use bare ``assert`` (pytest
rewrites those in test modules so they fire regardless of ``-O``); the VALUE
checks never rely on a bare assert.

⭐ **Widened and re-derived 2026-08-11.**  Both ladders grew (sphere
N ∈ {4,8,16,32,64}, cylinder n_φ ∈ {4,6,8,10,16,18,32,64} — **both
parities of M**, where the committed pair {8,16} was a single congruence
class), both flat tolerances became DERIVED functions of the mechanism
(:func:`_sphere_tau_atol` / :func:`_cyl_tau_atol` — the old
``atol=1e-13`` was a latent FALSE RED at N ≥ 64, `[M]` 2.247e-13), and
the sphere row gained the negative control it lacked.

Mutation-proof — measured 2026-08-11 (13 in-process mutations, 298-row
scope; definitions in
``tests/sn/sweep/test_angular_cell_partition.py``):

====================================  == == == == === == == == == == ==
row                                   MC M1 M2 M3 M3b M4 M5 M6 M7 M8 M9
====================================  == == == == === == == == == == ==
``sphere_tau_matches…reference``       5  5  5  5   .  .  .  .  .  5  5
``cyl_tau_equals_the_ANALYTIC…``       7  .  .  .   .  7  8  8  8  8  8
====================================  == == == == === == == == == == ==

⛔ The **cylinder** row is one of only **6 pre-existing catchers (of 298
rows)** for each of the chord-partition revert (M4), the shifted-edge
partition (M7) and the orientation flip (M8) — and before the widening
it contributed **2** of those.  That is the measurement the new
``test_angular_cell_partition.py`` module exists to change.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.angular.closure import (
    IdentityAngularClosure,
    morel_montry_tau_per_level,
)
from tests.sn._test_helpers import placeholder_materials

_EPS = float(np.finfo(float).eps)


def _sphere_tau_atol(weights: np.ndarray) -> float:
    r"""Derived tolerance for the sphere τ row: :math:`16\,\varepsilon /
    w_{\min}`.

    τ divides an edge difference by the cell width, so an edge
    discrepancy of :math:`O(\varepsilon)` (the pairwise-vs-sequential
    summation noise that IS this reference's independence) is amplified
    by :math:`1/w_m`, and the narrowest Gauss–Legendre cell sits at the
    poles.  `[M]` 2026-08-11 the measured ratio
    ``max|ref − prod| ÷ (eps/w_min)``:

    ====  ============  ==================
    N     max gap       gap ÷ (eps/w_min)
    ====  ============  ==================
    4     ``0.000e+00``  0.00
    8     ``1.776e-15``  0.81
    16    ``6.550e-15``  0.80
    32    ``6.106e-15``  0.19
    64    ``2.247e-13``  1.81
    128   ``1.880e-12``  3.81
    ====  ============  ==================

    ⟹ the factor 16 carries ≥ 4× margin at N = 128 and ≥ 8× at every
    order this module runs.  ⛔ **This replaces a flat ``atol=1e-13``**,
    which was simultaneously ~450× too loose at N = 8 and a **false red**
    from N = 64 (measured ``2.247e-13``) — the very widening its own
    docstring warned a future row would need.
    """
    return 16.0 * _EPS / float(np.min(weights))


def _cyl_tau_atol(M: int) -> float:
    r"""Derived tolerance for the cylinder τ row: :math:`40\,M\,
    \varepsilon`.

    :math:`\tau_m - \tfrac12 = \tfrac12\cot\omega_m\tan(\Delta\omega/4)`,
    and near the arc's endpoints :math:`\cot` has condition
    :math:`\csc^2\omega \approx (2M/\pi)^2` while
    :math:`\tan(\Delta\omega/4) \approx \pi/(4M)` — so an
    :math:`O(\varepsilon)` error in :math:`\omega` (from ``arctan2``)
    reaches τ multiplied by :math:`\approx M/(2\pi)`.  The agreement
    therefore DEGRADES linearly in :math:`M`; `[M]` 2026-08-11:

    =======  ====  ============  =================
    n_φ      M     max gap       gap ÷ (M·eps)
    =======  ====  ============  =================
    4        2     ``1.110e-16``  0.25
    6        3     ``1.110e-16``  0.17
    8        4     ``2.220e-16``  0.25
    10       5     ``4.441e-16``  0.40
    16       8     ``7.772e-16``  0.44
    18       9     ``5.551e-16``  0.28
    26       13    ``3.664e-15``  1.27
    32       16    ``7.438e-15``  2.09
    64       32    ``2.270e-14``  3.20
    128      64    ``1.135e-13``  7.99
    =======  ====  ============  =================

    ⟹ the factor 40 carries ≥ 5× margin at M = 64 and ≥ 12× at every
    order this module runs, while a flat ``atol`` would be a false red
    somewhere above M = 64.
    """
    return 40.0 * M * _EPS


# ═══════════════════════════════════════════════════════════════════════
# SPHERE — closure producer == structurally-independent reference (0-ULP)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("N", [4, 8, 16, 32, 64])
def test_sphere_tau_matches_independent_reference(N):
    r"""Closure-produced τ == a HAND-WRITTEN cumulative-weight reference.

    ⛔ **Reference replaced at Q5.6.4 (2026-08-11).** This row compared
    against ``contamination.morel_montry_weights``, described as "a
    DIFFERENT code path to the SAME unclamped BMC weight".  That module's
    successor **delegates to production**, so the comparison would now be
    τ against itself through a wrapper — green forever, and unable to
    detect the drift its name advertises (`coding-standards`: a rewire can
    demote a gate's claim class without touching one line of the body).

    The reference is therefore authored HERE: BMC 2010 Eq. 12 —
    :math:`\mu_{1/2} = -1`, :math:`\mu_{m+1/2} = \mu_{m-1/2} + w_m` — then
    P2 (Eq. 42).

    ⚠ **Not bit-exact, and the reason is the independence itself.**
    ``np.cumsum`` sums pairwise while the producer accumulates
    sequentially, so the two edge ladders differ by FP association — which
    is precisely what makes this a second computation rather than a copy.
    `[M]` the gap GROWS with N, because τ divides that
    :math:`O(\varepsilon)` edge noise by the cell width and the narrowest
    Gauss–Legendre cell narrows with N:

    ===== ===================  ==========
    N     max|ref − prod|      ULP
    ===== ===================  ==========
    4     ``0.000e+00``            0
    8     ``1.776e-15``           16
    16    ``6.550e-15``           59
    32    ``6.106e-15``           55
    64    ``2.247e-13``         2024
    128   ``1.880e-12``        16930
    ===== ===================  ==========

    ⛔ **The flat ``atol=1e-13`` this row used until 2026-08-11 was a
    FALSE RED from N = 64** (measured ``2.247e-13``) — exactly the
    widening its own docstring predicted a future row would need.  The
    bound is now DERIVED, :func:`_sphere_tau_atol`
    (:math:`16\varepsilon/w_{\min}`), which tracks the mechanism instead
    of a magic floor, and the ladder runs to N = 64.

    **Leg (b) — NEGATIVE control** (`vv-principles` #19), new at
    2026-08-11.  Until now this row pinned the *implementation* of the
    cumulative-weight partition and carried **zero** information about
    the *choice*: it would have passed identically had the producer used
    equal-width cells, because the reference would have been re-derived
    from … the same wrong convention only if a human wrote it that way.
    The plausible-wrong alternative — ``linspace(-1, 1, N+1)`` — is
    computed inline and τ must DIFFER from it (`[M]` by 0.159 / 0.646 /
    1.525 / 3.206 at N = 4/8/16/32).

    **Cannot catch**: the convention *being wrong in the literature's own
    terms* — BMC 2010 Eq. 12 is transcribed here as well as in
    production, so a mis-reading of the source is shared.  The corroborating
    check is Lathrop 2000's :math:`\sum\Delta\mu_m = 2`, gated as the
    closing identity in
    ``tests/sn/sweep/test_angular_cell_partition.py``.
    """
    quad = Quadrature.gauss_legendre(N)
    # (a) hand-written reference: cumulative-weight edges from −1, then P2.
    w = np.asarray(quad.weights)
    mu_edge = np.concatenate(([-1.0], -1.0 + np.cumsum(w)))
    tau_ref = (np.asarray(quad.mu_x) - mu_edge[:-1]) / np.diff(mu_edge)
    (tau_close,) = morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)

    np.testing.assert_allclose(
        tau_close, tau_ref, rtol=0.0, atol=_sphere_tau_atol(w),
        err_msg=(
            f"sphere closure τ != the hand-written cumulative-weight "
            f"reference (N={N}); both are unclamped BMC-2010 Eq.12/42 "
            f"weights and must agree to FP-association noise"
        ),
    )

    # (b) NEGATIVE control: the naive equal-width partition.
    uniform_edge = np.linspace(-1.0, 1.0, N + 1)
    tau_uniform = (
        (np.asarray(quad.mu_x) - uniform_edge[:-1]) / np.diff(uniform_edge)
    )
    gap = float(np.max(np.abs(tau_ref - tau_uniform)))
    if gap <= 0.1:
        pytest.fail(
            f"sphere N={N}: the cumulative-weight reference is "
            f"indistinguishable from the naive EQUAL-WIDTH partition "
            f"(max gap {gap:.3e}) — leg (a) cannot then be evidence "
            f"about the partition CONVENTION, only about its arithmetic"
        )


# ═══════════════════════════════════════════════════════════════════════
# CYLINDER — closure producer == clamp(independent reference) + neg control
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [4, 6, 8, 10, 16, 18, 32, 64])
def test_cyl_tau_equals_the_ANALYTIC_closed_form_not_the_chord_convention(
    n_phi,
):
    r"""Cylinder: closure τ == the ANALYTIC P2 closed form on the arc.

    ⛔ **Re-posed at Q5.6.4 (2026-08-11).** This row was
    ``test_cyl_tau_clamp_is_the_only_difference_from_reference``, and its
    whole thesis — *"closure τ == clip(reference τ_raw, ½, 1), and must
    DIFFER from the raw reference where the clamp bites"* — became
    vacuous when the ``[½, 1]`` absorber retired: there is no clamp left
    to be the only difference.  It is **kept, not deleted**, because it
    is the cylinder's only vv-L11 producer gate (two independently
    produced sides).

    ⭐ **What changed about the INDEPENDENCE, and it is an upgrade.** The
    old reference was ``contamination.morel_montry_weights`` — a second
    *procedural* implementation of the same edge recursion, i.e. only
    procedural independence (`vv-principles` L11), and one that silently
    became WRONG when the partition moved.  The reference is now
    the **analytic closed form** obtained by putting the ω-midpoint
    partition through P2 (BMC Eq. 43):

    > :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\,
    > \tan(\Delta\omega/4)`

    derived by hand from the arc geometry, sharing no code path with the
    producer.  That is *structural* independence.

    ⚠ **NEGATIVE control (`vv-principles` #19 — the positive leg alone
    cannot show a gate is loaded on the structure it is credited with).**
    The retired **chord** (η-midpoint) convention is computed inline
    below and the closed form must DIFFER from it — otherwise this row
    would pass just as happily against the partition the carve replaced,
    and would be certifying nothing about the partition choice.

    ⚠ ``angular_differencing.morel_montry_weights`` is not used here — it
    delegated to production, so comparing against it would have been
    comparing τ with itself through a wrapper (see the module-level note).
    It was retired outright on 2026-08-12.

    ⭐ **Widened at 2026-08-11** from ``n_phi ∈ {8, 16}`` (both
    :math:`M` EVEN) to ``{4, 6, 8, 10, 16, 18, 32, 64}``, i.e.
    :math:`M \in \{2,3,4,5,8,9,16,32\}` — **both parities**.  An
    all-even ladder is a single congruence class (`vv-principles` #13's
    refinement-ladder sharpening): odd :math:`M` puts an ordinate exactly
    at :math:`\omega = \pi/2`, where :math:`\cot\omega = 6.1e{-17}` and
    :math:`\tau - \tfrac12` is :math:`0.0` or :math:`\pm 1.1e{-16}` — the
    one place the closed form is evaluated at a near-cancellation, and
    the equality case of the orientation law in
    ``tests/sn/sweep/test_angular_cell_partition.py``.  The tolerance is
    now DERIVED (:func:`_cyl_tau_atol`) rather than a flat ``1e-13``,
    because the agreement degrades like :math:`M\varepsilon`.

    ⚠ **The** :math:`M = 2` **row (n_phi = 4) is a labelled NON-catcher
    for leg (b)** and skips it: there the chord and ω-midpoint partitions
    coincide to 1 ULP, so no ``folded_product(·, 4)`` fixture can see the
    partition choice at all.  It is kept for leg (a) (the closed form
    still holds) and the degeneracy itself is gated in
    ``test_angular_cell_partition.py::test_the_M2_fold_is_BLIND_to_the_partition_choice``.

    **Cannot catch**: a producer wrong only on a NON-equispaced arc (the
    reference is the equispaced closed form, and every shipped rule is
    equispaced in :math:`\varphi`); anything on the sphere arm; and — for
    leg (a) alone — a march-orientation flip is caught, but only because
    the closed form is *signed*: the three properties the fold gates
    assert (:math:`\tau \in [0,1]`, :math:`\tau \in [\tfrac14,\tfrac34]`,
    :math:`\tau_m + \tau_{M-1-m} = 1`) are ALL invariant under
    :math:`\tau \to 1-\tau`.  That Mode-12 hole has its own dedicated
    catcher in ``test_angular_cell_partition.py``.
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    tau_close = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)

    for p, level_idx in enumerate(quad.level_indices):
        eta = quad.mu_x[level_idx]
        xi = quad.mu_y[level_idx]
        sin_theta = float(np.sqrt(1.0 - quad.mu_z[level_idx[0]] ** 2))
        M = len(level_idx)
        omega = np.arctan2(xi, eta)
        d_omega = np.pi / M

        # (a) the ANALYTIC reference — hand-derived, no shared code path.
        closed_form = 0.5 + 0.5 / np.tan(omega) * np.tan(d_omega / 4.0)
        np.testing.assert_allclose(
            tau_close[p], closed_form, rtol=0.0, atol=_cyl_tau_atol(M),
            err_msg=(
                f"cylinder level {p} (n_phi={n_phi}, M={M}): production τ "
                f"departs from the analytic P2 closed form "
                f"½ + ½·cot(ω)·tan(Δω/4)\n"
                f"production ={tau_close[p]}\nclosed form={closed_form}"
            ),
        )

        # (b) NEGATIVE control: the RETIRED chord (η-midpoint) partition.
        #     Structurally VACUOUS at M = 2 (the two coincide) — skipped
        #     there rather than asserted with a threshold it cannot meet.
        if M < 3:
            continue
        chord_edge = np.empty(M + 1)
        chord_edge[0] = -sin_theta
        chord_edge[M] = sin_theta
        chord_edge[1:M] = 0.5 * (eta[:-1] + eta[1:])
        chord_tau = (eta - chord_edge[:-1]) / np.diff(chord_edge)
        gap = float(np.max(np.abs(closed_form - chord_tau)))
        # `[M]` 2026-08-11 the τ-space gap GROWS with refinement (unlike
        # the edge-space gap, which shrinks like M⁻²): 3.166e-2 (M=3) /
        # 4.459e-2 (4) / 6.050e-2 (5) / 7.514e-2 (8) / 7.694e-2 (9) /
        # 8.137e-2 (16) / 8.285e-2 (32) / 8.321e-2 (64), identical on
        # every level (τ is sinθ-independent).  1e-2 ⇒ ≥ 3× margin at the
        # tightest row and no upper-order false red.
        if gap <= 1e-2:
            pytest.fail(
                f"cylinder level {p} (n_phi={n_phi}, M={M}): the analytic "
                f"reference is indistinguishable from the RETIRED chord "
                f"convention (max gap {gap:.3e}) — this row cannot then be "
                f"evidence about the partition choice, which is the thing "
                f"the carve changed"
            )



# ═══════════════════════════════════════════════════════════════════════
# CARTESIAN — IdentityAngularClosure supplies neutral τ = 1.0 (no branch)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_identity_closure_tau_is_neutral_one():
    """Cartesian: IdentityAngularClosure supplies τ = 1.0 (no redistribution).

    POSITIVE control: the neutral element is produced WITHOUT a geometry
    branch — the closure TYPE IS the dispatch (Pattern 4).  τ = 1.0 makes
    ``c_out = α_out/1`` and ``c_in = α_in`` with α ≡ 0 (slab), so the M-M
    contribution vanishes identically.
    """
    nx = 8
    edges = np.linspace(0.0, 5.0, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    reduced = sn_mesh.reduced
    assert reduced is not None  # 1-D mesh => minted by the ctor (narrowing)
    closure = IdentityAngularClosure(
        reduced.angular, reduced.redistribution_pairing,
    )

    (tau,) = closure._tau_per_level
    np.testing.assert_array_equal(
        tau,
        np.ones_like(tau),
        err_msg=(
            "IdentityAngularClosure τ is not the neutral 1.0 — Cartesian must "
            "carry no angular redistribution"
        ),
    )
    # No geometry dispatch: morel_montry_tau_per_level raises on Cartesian;
    # Identity supplies the neutral τ by its TYPE, not by a coord branch.
    with pytest.raises(ValueError):
        morel_montry_tau_per_level(quad, CoordSystem.CARTESIAN)

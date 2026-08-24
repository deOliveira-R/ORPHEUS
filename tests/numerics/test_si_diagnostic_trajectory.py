r"""CS3 capture gate — the SI convergence-diagnostic trajectory, frozen PRE-carve.

**What this file is for.** Campaign-1 phase **CS3** (the cone overturn,
``.claude/plans/space_and_kernel_binding_campaign.md`` §4) RELOCATES the three
iterate diagnostics off the ``Displacement`` type and onto the iteration layer:

* ``Displacement.contraction_ratio`` (retired with its type at step 3)
  → the ρ trajectory computed by :class:`~orpheus.numerics.iteration.SourceIteration`,
* ``Displacement.true_error_estimate`` (retired likewise) → ``‖Δψ‖/(1−ρ)``
  on the iteration layer,
* ``Displacement.where_largest`` (retired likewise) → the per-entry
  convergence map.

At capture time (``000cf144``, pre-carve) the SI loop minted a typed
displacement every pass (then-``iteration.py:787``:
``displacement = psi - psi_prev``), found its bulk leaf
(``_flux_displacement_leaf``, then-``:421``) and called
``disp_leaf.contraction_ratio(_prev_disp_leaf)`` — which was
``self.l2 / previous.l2`` with ``Field.l2 = space.norm(values)`` on the
**INTERIOR LEAF**. This module froze what that produced on a c→1 fixture, so
the post-carve surface can be shown to reproduce it.

**Re-pointed at CS3 step 1 (2026-08-19), numbers UNMOVED.** The surface is now
the iteration record: ``record.contraction_ratios`` (derived from
``record.increment_norms``), ``record.increment_norms[-1]``,
``record.true_error_estimate()``, and the per-entry map is
:meth:`~orpheus.numerics.field.Field.where_largest` on the interior leaf of
the final increment (replayed by :func:`_replay_final_increment` through the
same production builders). **The frozen numbers below did NOT move**; the
fixture builder did not move; the control and activation legs did not move —
the value-neutrality claim this module exists to certify.

**Claim kind = RECORD** (a frozen trajectory of what the code produced on
2026-08-19 at ``000cf144``), so on its own it says *something changed*, never
*which side is right*. Two things keep it honest:

* the REFERENCE anchor for ρ lives one directory over and is untouched by this
  file — :mod:`tests.sn.solve.test_si_convergence_diagnostics` pins
  ρ ≈ c = Σ_s/Σ_t (Adams & Larsen 2002) in bands on a homogeneous slab. That
  file must be RE-POINTED by step 1, not retired: it is this pin's anchor.
* :func:`test_the_pin_discriminates_the_composite_norm_convention` is an
  in-file CONTROL proving the pin can RED — see below.

**Why these numbers can be trusted to a tight tolerance, and what that tolerance
is chosen to catch.** Measured on this fixture at ``000cf144``:

* interior-leaf **metric** norm vs interior-leaf **flat** ``np.linalg.norm``:
  at the CS3 freeze these agreed to ≤1 ULP (``2.29e-16`` max relative)
  because the then-``angular_flux`` tag space carried
  ``inner_product_weights is None`` (Euclidean) — so battery row M2 measured
  the pin BLIND to a metric↔flat swap. ⛔ SUPERSEDED at CS4b S2 (2026-08-22):
  the space is the carrier's axis-built ``angular_bulk_space`` carrying
  ``V_cell × w_n``, the two norms now differ at O(1), and the same M2
  mutation (``Field.l2 → np.linalg.norm``) REDS this module — the pin is
  metric-LOADED, which is the CS4b evidence pair (blind before, loaded
  after) the campaign's M7 battery arm records.
* interior-leaf metric vs **whole-composite** flat norm (``_l2_norm`` of the
  raw increment — the spelling the SI loop had in hand at capture time, which
  additionally ravels the boundary trace block): max relative difference
  ``4.71e-3``. **That** is the convention error this pin exists to catch, and
  it sits 9 orders above the tolerance.

⚠ **Phase-ordering hazard, measured — RESOLVED.** As written at CS3 this block
predicted: recomputing the trajectory under a physical ``V_cell × w_n`` metric
moves ρ by up to ``1.12e-3`` relative — 9 orders above ``rtol=1e-12`` — and
assigned the re-derivation to CS2 ("✅ RULED (user, 2026-08-19) — the
diagnostic is defined on the **SPACE norm**, so CS2 owns re-deriving these
frozen numbers"). The metric arrived EARLIER than the ownership sentence
assumed: **CS4b S2** (the field-layer space-source flip) installed it, so
CS4b owned and executed the re-derivation on 2026-08-22 (Q4, verification
plan §2.1 G-M1..3) — measured move ``1.117e-3``, inside the prediction; the
regeneration recipe is the provenance comment on the frozen block below.

**Measured mutation battery** (in-process pytest plugins; ``vv`` §0.5 — a gate
is not evidence until a named mutation reddens it). At ``000cf144``, 5 tests,
0.59 s green:

===  ===========================================================  ==================================
 #   mutation                                                     result
===  ===========================================================  ==================================
 M1  ρ from the WHOLE COMPOSITE's flat norm (the relocation       **2 failed** — the ρ pin AND the
     error this module exists for)                                ``‖Δψ‖`` pin
 M2  ``Field.l2`` → ``np.linalg.norm(values)`` (metric → flat)    **5 passed** — the declared ≤1 ULP
                                                                  blindness, MEASURED not assumed
 M3  one spurious leading ratio (recording-cadence drift)         **1 failed** — the length leg
 M4  ``where_largest`` ravels before locating (layout loss)       **1 failed** — the map leg only
 M5  every norm × (1 + 1e-9) — the positive control               **1 failed** — the ``‖Δψ‖`` leg ONLY
===  ===========================================================  ==================================

⛔ **M5 is a stabiliser fact about this module, stated because a reader would
otherwise over-credit the ρ pin: the ρ trajectory is INVARIANT under any factor
applied uniformly to the norm**, because ρ is a ratio and a common factor
cancels exactly (``vv`` Mode 12). That whole error class is caught only by
:func:`test_true_error_estimate_and_norm_reproduce`, which is why ``‖Δψ‖`` is
pinned in its own right and must never be folded into the ρ test.

Marks ``foundation`` — a software-invariant / wiring gate, no physics
``:label:`` claim, so no ``verifies(...)``. Every assertion is
``np.testing.assert_*`` or an explicit ``raise`` (``vv`` Mode 8: bare ``assert``
is stripped under the canonical ``python -O``).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver, _within_group_si
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import (
    AngularBoundarySourceSink,
    AngularSourceSink,
)
from orpheus.transport.timed_full_field import TimedFullField


pytestmark = pytest.mark.foundation


# ── The fixture, stated so the numbers below are re-derivable ────────────
#
# vv-principles §"Numbers carry their CONFIGURATION". 2 GROUPS and
# HETEROGENEOUS deliberately: a 1-group homogeneous slab makes a per-group
# norm and a total norm the SAME number, so a relocation that reduced over the
# wrong axis would be invisible (Mode 12 at the fixture).
_C = 0.99                    #: scattering ratio Σ_s/Σ_t in BOTH groups of BOTH zones
_NX = 40                     #: cells (2 zones of 20)
_N_ORD = 8                   #: Gauss-Legendre ordinates
_WIDTH = 10.0                #: cm — 10 mean free paths in group 0 of zone 0
_MAX_ITER = 12               #: SI passes; the budget IS the trajectory length + 1
_TOL = 1e-14                 #: far below anything reachable in 12 passes at ρ ≈ 0.97
                             #  ⟹ the stop never fires ⟹ the length is deterministic


def _mixture(sig_t, s00, s01, s11) -> Mixture:
    """2-group mixture with ASYMMETRIC scattering (downscatter only, no upscatter).

    ``SigS[g_from, g_to]`` — rows are source groups (the ORPHEUS convention,
    ERR-002). Asymmetric on purpose: a symmetric matrix is the degenerate case
    for every transpose-convention error.
    """
    zero = np.zeros(2)
    scatter = np.array([[s00, s01], [0.0, s11]])
    sig_t = np.asarray(sig_t, dtype=float)
    return Mixture(
        SigC=sig_t - scatter.sum(axis=1),   # no fission ⟹ Σ_t = Σ_c + Σ_s
        SigL=zero.copy(), SigF=zero.copy(), SigP=zero.copy(), SigT=sig_t,
        SigS=[csr_matrix(scatter)], Sig2=csr_matrix((2, 2)), chi=zero.copy(),
    )


#: zone 0 — Σ_t = (1.0, 2.0), c = 0.99 in both groups
_FUEL = _mixture([1.0, 2.0], 0.60, 0.39, 1.98)
#: zone 1 — Σ_t = (1.5, 1.0), c = 0.99 in both groups
_MODERATOR = _mixture([1.5, 1.0], 1.00, 0.485, 0.99)


# ── RE-DERIVED at CS4b S2 (2026-08-22) under the PHYSICAL space metric ──
#
# The CS3 freeze (``000cf144``, 2026-08-19, Euclidean interior norm) recorded
# ratios starting ``0.9759292618901443, 1.0358831699038207, …`` with
# ``‖Δψ‖ = 8.848958875813311`` and estimate ``308.0488097665712``. CS4b S2
# flipped the field space source to the carrier's axis-built mints, which
# carry the ``V_cell × w_n`` Hilbert metric — so ``Field.l2`` became the
# PHYSICAL norm and these diagnostics legitimately moved: max relative ρ move
# ``1.117e-3`` ([M] 2026-08-22 — the hazard block below had predicted "up to
# ``1.12e-3``" from the pre-carve probe), ``where_largest`` unchanged (the
# per-entry map reads values, not norms). Licence: the structurally
# independent ρ ≈ c = Σ_s/Σ_t anchor
# (:mod:`tests.sn.solve.test_si_convergence_diagnostics`, Adams & Larsen
# 2002) re-run FIRST under the new metric — 4 passed — never old-vs-new.
#
# Produced by the code path described in the module docstring: a ratio of the
# INTERIOR leaf's ``space.norm`` between successive SI passes. NOTE the
# trajectory crosses 1.0 in the transient — the increment GROWS for two
# passes before contracting — which is what makes it a discriminating
# fingerprint rather than a near-constant.
_CONTRACTION_RATIOS = (
    0.9748389411907307,
    1.03521237670015,
    1.0264602980995396,
    1.008420630121672,
    0.9942082644544892,
    0.9846326637492744,
    0.9785442985358258,
    0.9748296517576658,
    0.972677197134772,
    0.9715347620676633,
    0.9710349088657109,
)
#: ‖Δψ‖ of the LAST recorded increment (the interior leaf's space norm).
_LAST_DISPLACEMENT_L2 = 2.22362669469903
#: ``‖Δψ‖/(1−ρ)`` at the last recorded ρ — the c→1 false-convergence estimate.
_TRUE_ERROR_ESTIMATE = 76.76919379917577
#: ``where_largest(3)`` index tuples into the ``(n_ordinate, group, cell)`` layout.
_WHERE_LARGEST_3 = [(3, 1, 24), (4, 1, 26), (4, 1, 25)]

#: The pin's tolerance. Justified in the module docstring: 4 orders above the
#: ≤1-ULP norm-implementation floor, 9 orders below the 4.71e-3 composite-norm
#: convention error it exists to catch.
_RTOL = 1e-12


def _build_solver() -> SNSolver:
    mat_ids = np.zeros(_NX, dtype=int)
    mat_ids[_NX // 2:] = 1
    mesh = Mesh1D(
        edges=np.linspace(0.0, _WIDTH, _NX + 1), mat_ids=mat_ids,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    sn_mesh = SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=_N_ORD),
        {0: _FUEL, 1: _MODERATOR},
    )
    return SNSolver(sn_mesh, inner_solver="source_iteration", scattering_order=0)


def _run_si():
    """Run the within-group SI through the PRODUCTION single-source builders.

    ``build_within_group_system`` + ``_within_group_si`` are the same two
    builders the production fixed-source and eigenvalue paths call
    (``solver.py:1058``) — no SI construction is duplicated here.

    Returns the solve's :class:`IterationRecord` (the post-carve diagnostic
    surface) together with the composite final increment Δψ from
    :func:`_replay_final_increment` (the record holds only FLOATS, so the
    per-entry map and the boundary-block control read the replayed field).
    """
    solver = _build_solver()
    sn_mesh = solver.sn_mesh
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    si, _base, _gains, windowed = _within_group_si(
        system, sn_mesh, inner_schedule=solver.inner_schedule,
        max_iter=_MAX_ITER, tol=_TOL,
    )
    if windowed:
        raise AssertionError(
            "the 1-D slab windowed the SI iterate — windowing is 2-D Cartesian, "
            "so the interior leaf is no longer the AngularFlux this pin froze."
        )
    q_ext = TimedFullField(
        interior=AngularSourceSink.from_isotropic(
            np.ones((sn_mesh.ng, *sn_mesh.spatial_shape)), sn_mesh,
        ),
        boundary=AngularBoundarySourceSink.zeros(sn_mesh.angular_trace),
        _history=(), history_depth=2,
    )
    guess = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space,
    )
    _, record = si.solve(q_ext, initial_guess=guess)
    return record


@pytest.fixture(scope="module")
def run():
    """One solve + one replay for the whole module (`[M]` ~1.2 s total)."""
    class _Run:
        record = _run_si()
        increment = _replay_final_increment()
    return _Run


# ═══════════════════════════════════════════════════════════════════════
# ⟸ THE FOUR READS BELOW ARE WHAT STEP 1 RE-POINTS. Nothing else moves.
# ═══════════════════════════════════════════════════════════════════════


def _diagnostics(run):
    """The CS3-relocated surface, isolated to one function.

    Post step 1: ρ and the geometric-tail estimate are DERIVED on the
    :class:`~orpheus.numerics.convergence.IterationRecord` from its
    ``increment_norms`` trajectory; the per-entry map is
    :meth:`~orpheus.numerics.field.Field.where_largest` on the interior leaf
    of the replayed final increment. The RETURNED NUMBERS MUST NOT MOVE
    relative to the pre-carve freeze below.
    """
    record = run.record
    if not record.increment_norms:
        raise AssertionError(
            "increment_norms is EMPTY — the SI loop recorded no iterate "
            "increment at all (the diagnostic went silent)."
        )
    return {
        "ratios": tuple(record.contraction_ratios),
        "l2": float(record.increment_norms[-1]),
        "true_error": float(record.true_error_estimate()),
        "where_largest": list(run.increment.interior.where_largest(3)),
    }


# ── 1. The pin ──────────────────────────────────────────────────────────


def test_contraction_ratio_trajectory_reproduces(run) -> None:
    r"""The ρ trajectory reproduces the frozen pre-carve values.

    The keystone. ``vv`` claim kind = RECORD: a red says the relocation changed
    the number, and the CONTROL below says the pin is sensitive to the one
    convention change a relocation is most likely to make.

    ⛔ **What this row structurally CANNOT see** (`[M]` battery M5): a factor
    applied uniformly to every norm. ρ is a ratio, so any common scale cancels
    exactly — a 1e-9 relative perturbation of every ``l2`` leaves this row
    GREEN. :func:`test_true_error_estimate_and_norm_reproduce` is the catcher
    for that class; do not fold the two together.
    """
    got = _diagnostics(run)["ratios"]
    if len(got) != len(_CONTRACTION_RATIOS):
        raise AssertionError(
            f"trajectory LENGTH moved: {len(got)} recorded vs "
            f"{len(_CONTRACTION_RATIOS)} frozen — the recording CADENCE changed "
            f"(which pass records, or whether the zero-norm guard still skips "
            f"the first), independently of any value."
        )
    np.testing.assert_allclose(
        np.asarray(got), np.asarray(_CONTRACTION_RATIOS), rtol=_RTOL, atol=0.0,
        err_msg=(
            "the SI contraction-ratio trajectory moved. Check WHICH norm and "
            "WHICH block the relocated diagnostic uses: the frozen values are "
            "the ratio of the INTERIOR LEAF's space norm. Ratios of the whole "
            "composite's flat norm differ by up to 4.7e-3 (see the control)."
        ),
    )


def test_true_error_estimate_and_norm_reproduce(run) -> None:
    r"""``‖Δψ‖`` and ``‖Δψ‖/(1−ρ)`` reproduce the frozen values.

    The c→1 amplification is the whole point of the diagnostic (Signature 9),
    so it is pinned in its own right rather than inferred from ρ.
    """
    got = _diagnostics(run)
    np.testing.assert_allclose(
        got["l2"], _LAST_DISPLACEMENT_L2, rtol=_RTOL, atol=0.0,
        err_msg="‖Δψ‖ of the last recorded increment moved",
    )
    np.testing.assert_allclose(
        got["true_error"], _TRUE_ERROR_ESTIMATE, rtol=_RTOL, atol=0.0,
        err_msg="the true-error estimate ‖Δψ‖/(1−ρ) moved",
    )


def test_where_largest_reproduces(run) -> None:
    r"""``where_largest(3)`` reproduces the frozen index tuples.

    The convergence MAP is the diagnostic whose relocation is most likely to
    lose information: indices are into the leaf's own ``(n_ordinate, group,
    cell)`` layout, and a relocation that ravels the composite first would
    return flat or shifted indices while every NORM above still reproduced.
    """
    got = _diagnostics(run)["where_largest"]
    if [tuple(int(i) for i in t) for t in got] != _WHERE_LARGEST_3:
        raise AssertionError(
            f"where_largest(3) moved: {got} vs frozen {_WHERE_LARGEST_3} — the "
            f"per-entry convergence map lost the leaf's index layout."
        )


# ── 2. The control — the pin CAN red (vv #17 / #19) ─────────────────────


def test_the_pin_discriminates_the_composite_norm_convention(run) -> None:
    r"""The rival convention differs from the pin by ≥ 1e-3 on THIS fixture.

    ``vv`` #19: a positive reading cannot tell a loaded gate from a blind one —
    only the deliberately-wrong structure can. The rival is the natural
    relocation spelling: ``_l2_norm`` of the WHOLE composite increment,
    which additionally ravels the boundary-trace block
    (``full_field.py``'s ``to_flat`` = ``[interior.ravel(), boundary]``).

    `[M]` at ``000cf144``: max relative difference ``4.71e-3``, i.e. 9 orders
    above ``rtol = 1e-12``. If this leg ever reds, the FIXTURE has stopped
    discriminating (e.g. a vanishing boundary block of the increment — an
    all-reflective or zero-inflow config would do it) and the pin above has
    silently become blind to the convention it was written for.
    """
    if not run.record.increment_norms:
        raise AssertionError("no increment recorded — control cannot run")
    boundary_norm = float(
        np.linalg.norm(np.asarray(run.increment.boundary.values, dtype=float))
    )
    if boundary_norm <= 0.0:
        raise AssertionError(
            "the boundary increment block is ZERO on this fixture, so the "
            "interior-only and whole-composite conventions coincide and the "
            "pin cannot see the difference. Restore a config with non-trivial "
            "boundary movement (vacuum walls + an interior source)."
        )
    interior_norm = float(run.record.increment_norms[-1])
    composite_norm = float(np.hypot(interior_norm, boundary_norm))
    separation = abs(composite_norm - interior_norm) / interior_norm
    if separation < 1e-3:
        raise AssertionError(
            f"the two norm conventions differ by only {separation:.3e} on this "
            f"fixture — below the 1e-3 the pin needs to discriminate them."
        )


def _replay_final_increment():
    """Replay the fixture's SI loop and return the COMPOSITE final increment.

    The record holds only floats (deliberately — O(1) memory per solve), so
    the per-entry map leg and the boundary-block control read the increment
    FIELD from this replay, which runs the same production builders and the
    same arithmetic as :func:`_run_si`'s solve. Deterministic: same operators,
    same source, same zero guess, and the ``tol=1e-14`` stop never fires.
    """
    solver = _build_solver()
    sn_mesh = solver.sn_mesh
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    si, _base, _gains, _windowed = _within_group_si(
        system, sn_mesh, inner_schedule=solver.inner_schedule,
        max_iter=_MAX_ITER, tol=_TOL,
    )
    q_ext = TimedFullField(
        interior=AngularSourceSink.from_isotropic(
            np.ones((sn_mesh.ng, *sn_mesh.spatial_shape)), sn_mesh,
        ),
        boundary=AngularBoundarySourceSink.zeros(sn_mesh.angular_trace),
        _history=(), history_depth=2,
    )
    psi = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space,
    )
    increment = None
    for _ in range(_MAX_ITER):
        previous = psi
        rhs = q_ext
        for gain in si.gains:
            rhs = rhs + gain.apply(psi)
        psi = si.A_inv.apply(rhs, initial_guess=previous)
        increment = psi - previous
    if increment is None:                          # pragma: no cover — _MAX_ITER ≥ 1
        raise AssertionError("the replay ran zero passes")
    return increment


# ── 3. Activation — the fixture really is in the c→1 regime ─────────────


def test_the_fixture_is_in_the_c_to_one_regime(run) -> None:
    r"""ρ is near 1 and the true-error amplification is ≫ 1.

    Without this the pin could be frozen on a fast-converging fixture where
    ``1/(1−ρ) ≈ 1`` and the true-error estimate carries no information —
    reproducing a number that means nothing. The bound ``ρ < c`` is the
    Adams & Larsen (2002) contraction bound ``ρ ≈ max Σ_s/Σ_t``; the REFERENCE
    claim itself lives in
    :mod:`tests.sn.solve.test_si_convergence_diagnostics`.
    """
    got = _diagnostics(run)
    rho = got["ratios"][-1]
    if not (0.9 < rho < _C):
        raise AssertionError(
            f"ρ_last = {rho:.6f} is outside (0.9, c={_C}) — the fixture is no "
            f"longer the c→1 regime this pin was frozen in."
        )
    amplification = 1.0 / (1.0 - rho)
    if amplification < 20.0:
        raise AssertionError(
            f"1/(1−ρ) = {amplification:.2f} < 20 — the true-error amplification "
            f"is too weak for the estimate to be a meaningful pin."
        )
    np.testing.assert_allclose(
        got["true_error"], got["l2"] / (1.0 - rho), rtol=1e-12, atol=0.0,
        err_msg="the true-error estimate is not ‖Δψ‖/(1−ρ) of the SAME iterate",
    )

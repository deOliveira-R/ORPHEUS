r"""AC-b / AC-b′ — **one splitting per posed equation**, and the law
:math:`A = M - N` that makes a splitting a splitting.

Campaign: ``.claude/plans/operator_strategy_realization_campaign.md``.
Normative gate spec: ``.claude/plans/campaign_verification_plan.md`` §1.
Phase **P0** — gates only, no production change.

Why this file exists (the motivating defect, R7)
================================================

:class:`~orpheus.sn.coupled_system.WithinGroupSystem` welds two independent
stages into one record: the **posing** (``loss`` — what the physics *is*) and
the **splitting** (``implicit_operator`` / ``explicit_gains`` — which part is
solved implicitly).  Because the two share a record, there is no named
boundary at which "strategy may enter" can be asserted — and because there is
no boundary, **a second splitting grew beside the first and the first went
silently stale**:

.. code-block:: text

   inner_schedule    driver actually runs                record advertises
   ---------------   ---------------------------------   -----------------------------
   jacobi            StreamingCollisionOperator + SNBoundaryOperator   ← the same objects
   gauss_seidel      ScheduledInvertibleOperator + SNMaskedBoundaryOperator
                                                          StreamingCollisionOperator
                                                          + SNBoundaryOperator

``_select_si_splitting`` (``sn/solver.py:706``) re-derives a splitting the
record never hears about.  Nothing is *numerically* wrong today — both
splittings are consistent — but the record's claim is false, and a consumer
that trusts ``record.implicit_operator`` (the spectral gate, an admission
check, a preconditioner) reads an operator the solver does not run.

Ordering constraint **O-1**, binding
====================================

This gate is written **RED, before anything touches**
``WithinGroupSystem`` or ``_select_si_splitting``.  Fix R7 first and the fix
is unprovable — and a future re-split has no catcher.  The red row ships as
``xfail(strict=True)``: when P3/P5 lands, it XPASSes, which is a hard failure
that forces the marker's removal.  **The strict-xfail set IS the campaign's
todo list.**

Why the ``jacobi`` row must ship too
====================================

It is the **control leg**.  Without it a change that broke *both* arms would
read as "the gate was always partly red"; with it, the asymmetry is the
measurement.  Likewise the slab row below pins a trap that has already cost
this campaign one wrong answer.

Marks
=====

``foundation`` — software/architecture invariants of the operator algebra,
no theory ``:label:``, no ``verifies(...)`` (the verifies⊥level doctrine).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.solver import _within_group_krylov, _within_group_si
from tests.sn.architecture._config import (
    cart2d_seedless,
    isotropic_slab,
    random_state,
    reconstruction_residual,
    record_for,
    scattering_gain,
    seedless_implicit,
    sigma_s0_times_identity,
    slab_seedless,
    sphere_carrying,
    split_image,
    system_a,
)

pytestmark = pytest.mark.foundation

_SEED = 20260729

_R7_XFAIL = pytest.mark.xfail(
    strict=True,
    reason=(
        "R7 — the splitting is a TWIN: _select_si_splitting re-derives "
        "(ScheduledInvertibleOperator, SNMaskedBoundaryOperator) and the "
        "record still advertises (StreamingCollisionOperator, "
        "SNBoundaryOperator). Flipped by campaign P3/P5 (the partition "
        "becomes a value and the driver stops re-splitting). WHEN THIS "
        "XPASSES: P3/P5 has landed — delete this marker."
    ),
)


# ═════════════════════════════════════════════════════════════════════════
# AC-b — the driver runs the objects the record advertises
# ═════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    ("build_mesh", "inner_schedule"),
    [
        # CONTROL LEG: same geometry, same record, the other schedule.
        pytest.param(cart2d_seedless, "jacobi", id="cart2d-jacobi"),
        # THE RED: multi-D Cartesian + G-S is the only arm that re-splits.
        pytest.param(
            cart2d_seedless, "gauss_seidel",
            id="cart2d-gauss_seidel", marks=_R7_XFAIL,
        ),
        # The carrying arm returns the record's own splitting on BOTH rows —
        # inner_schedule is structurally inert there (G-S is multi-D
        # Cartesian ⟹ seedless, so it never reaches the selector).
        pytest.param(sphere_carrying, "jacobi", id="sphere-jacobi"),
        pytest.param(sphere_carrying, "gauss_seidel", id="sphere-gauss_seidel"),
    ],
)
def test_driver_consumes_the_records_own_splitting(build_mesh, inner_schedule):
    r"""The SI driver iterates the **same objects** the record advertises.

    Object identity (``is``), not value equality: two operators that happen
    to agree numerically today are still two operators, and the second one
    is where the drift lives.  This is the campaign's acceptance leg AC-b.
    """
    sn_mesh = build_mesh()
    record = record_for(sn_mesh)
    _si, driver_implicit, driver_gains, _windowed = _within_group_si(
        record, sn_mesh, inner_schedule=inner_schedule,
        max_iter=2, tol=1e-10,
    )

    if driver_implicit is not record.implicit_operator:
        pytest.fail(
            f"[{inner_schedule}] the driver's implicit operator is NOT the "
            f"record's: driver ran {type(driver_implicit).__name__}, record "
            f"advertises {type(record.implicit_operator).__name__} — the "
            f"posed record's splitting claim is stale (R7)."
        )
    advertised = tuple(map(id, record.explicit_gains))
    if tuple(map(id, driver_gains)) != advertised:
        pytest.fail(
            f"[{inner_schedule}] the driver's gains are NOT the record's: "
            f"driver ran {[type(g).__name__ for g in driver_gains]}, record "
            f"advertises {[type(g).__name__ for g in record.explicit_gains]} "
            f"— a second splitting exists beside the advertised one (R7)."
        )


@pytest.mark.parametrize(
    "build_mesh", [cart2d_seedless, sphere_carrying],
    ids=["cart2d", "sphere"],
)
def test_krylov_driver_consumes_the_records_own_splitting(build_mesh):
    r"""The Krylov driver is **green on both arms** — the twin is SI-specific.

    This is the second control leg, and it carries a positive architectural
    claim: ``_within_group_krylov`` is handed ``(record.implicit_operator,
    *record.explicit_gains)`` and stores them verbatim, so R7 is not a
    property of "within-group solves" but specifically of the SI schedule
    path.  If this row ever REDs, the defect has spread.
    """
    sn_mesh = build_mesh()
    record = record_for(sn_mesh)
    state = random_state(record, seed=_SEED)
    krylov = _within_group_krylov(
        record.implicit_operator, *record.explicit_gains,
        n_dof=int(state.to_flat().size), max_iter=2, tol=1e-10,
    )
    if krylov.A is not record.implicit_operator:
        pytest.fail(
            "the Krylov driver's operator is not the record's implicit "
            "operator — R7 has spread to the Krylov path."
        )
    if tuple(map(id, krylov.gains)) != tuple(map(id, record.explicit_gains)):
        pytest.fail(
            "the Krylov driver's gains are not the record's — R7 has spread "
            "to the Krylov path."
        )


def test_a_slab_hides_r7_the_documented_trap():
    r"""A 1-D slab reports GREEN on ``gauss_seidel`` — and that is a TRAP.

    ``_select_si_splitting`` falls back to Jacobi unless ``is_cartesian and
    not is_1d``, so a slab exercises the *control* arm under the *red* arm's
    name.  The campaign's first probe of R7 used a slab, measured ``True`` on
    both rows, and concluded there was no twin.

    Committing the trap as a named row makes that mistake unrepeatable: this
    test asserts the slab is green **and** says why, so a future reader who
    reaches for a slab fixture meets the explanation first.  It also has real
    teeth — if the 1-D fallback is ever removed, this REDs and the reader is
    told exactly which invariant moved.
    """
    sn_mesh = slab_seedless()
    record = record_for(sn_mesh)
    _si, implicit, gains, _w = _within_group_si(
        record, sn_mesh, inner_schedule="gauss_seidel", max_iter=2, tol=1e-10,
    )
    if implicit is not record.implicit_operator or tuple(map(id, gains)) != tuple(
        map(id, record.explicit_gains),
    ):
        pytest.fail(
            "the 1-D slab no longer falls back to Jacobi under "
            "inner_schedule='gauss_seidel' — the R7 observability precondition "
            "changed; re-read _select_si_splitting's geometry gate."
        )


# ═════════════════════════════════════════════════════════════════════════
# AC-b′ — the splitting LAW, with the σ_r bug as its teeth
# ═════════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    ("build_mesh", "exact"),
    [
        # Seedless: one flat operator sum — MEASURED exactly 0.0.
        pytest.param(cart2d_seedless, True, id="cart2d-seedless"),
        # Carrying: the 2x2 block grid re-associates the sums — MEASURED
        # 4.48e-17 (~2 ULP). A uniform bit-identity contract is REFUSED here
        # (vv-principles §bit-identity); the honest contract is per-arm.
        pytest.param(sphere_carrying, False, id="sphere-carrying"),
    ],
)
def test_reconstruction_identity_A_equals_M_minus_N(build_mesh, exact):
    r"""``A x == (M − ΣNᵢ) x`` — the law that makes a splitting a splitting.

    Green today on both arms (it is a **regression floor**, not a red gate).
    Its value is entirely in its teeth: the mutations below re-introduce the
    #215 σ_r defect, which shipped **46–56 % silent flux errors** and was
    gated by nothing.

    Tolerance is per-arm and MEASURED, never assumed: the seedless arm is a
    single flat operator sum and lands at exactly 0 ULP; the carrying arm's
    block grid re-associates and lands at ≤ 8 ULP (it fails at 4).
    """
    sn_mesh = build_mesh()
    record = record_for(sn_mesh)
    state = random_state(record, seed=_SEED)
    loss_image = record.loss.apply(state).to_flat()
    split = split_image(record, state)

    if exact:
        np.testing.assert_array_equal(
            loss_image, split,
            err_msg=(
                "A != M - N on the seedless arm, which is a single flat "
                "operator sum and MUST be bit-identical (measured 0.0). A "
                "non-zero defect here is an algebraic error, not FP drift."
            ),
        )
    else:
        np.testing.assert_array_almost_equal_nulp(loss_image, split, nulp=8)


@pytest.mark.parametrize(
    "build_mesh", [cart2d_seedless, sphere_carrying],
    ids=["cart2d-seedless", "sphere-carrying"],
)
def test_mutation_dropping_the_gains_reddens_the_splitting_law(build_mesh):
    r"""**M-5** — ``N := 0`` (claiming ``A = M``) REDs the law.

    This is the historical #215 bug in its purest form: the σ_r fold asserted
    the within-group loss WAS its own implicit operator.  MEASURED defect:
    1.00e-01 seedless, 1.83e-02 carrying — 15+ orders above the 0-ULP /
    8-nulp contracts the law is gated at.
    """
    sn_mesh = build_mesh()
    record = record_for(sn_mesh)
    state = random_state(record, seed=_SEED)
    defect = reconstruction_residual(record, state, gains=())
    if defect < 1e-3:
        pytest.fail(
            f"dropping every explicit gain moved the splitting law by only "
            f"{defect:.3e} — the law has no teeth on this fixture. Either "
            f"the gains no longer carry the scattering source, or the probe "
            f"state does not excite them."
        )


@pytest.mark.parametrize(
    "build_mesh", [cart2d_seedless, sphere_carrying],
    ids=["cart2d-seedless", "sphere-carrying"],
)
def test_mutation_sign_flipped_gain_reddens_the_splitting_law(build_mesh):
    r"""**M-7** — flipping ONE gain's sign REDs the law.

    Distinct from M-5: a sign flip preserves the gain's *presence* (an arity
    or None-block check still passes) and moves only the value, so it is the
    catcher for a convention drift rather than a dropped term.  MEASURED:
    3.32e-02 seedless, 3.66e-02 carrying.
    """
    sn_mesh = build_mesh()
    record = record_for(sn_mesh)
    state = random_state(record, seed=_SEED)
    flipped = (-record.explicit_gains[0], *record.explicit_gains[1:])
    defect = reconstruction_residual(record, state, gains=flipped)
    if defect < 1e-3:
        pytest.fail(
            f"flipping the first gain's sign moved the splitting law by only "
            f"{defect:.3e} — the law cannot see a gain-sign convention drift."
        )


# ═════════════════════════════════════════════════════════════════════════
# M-6 — the exact σ_r fold, with its Mode-9 control leg
# ═════════════════════════════════════════════════════════════════════════


def test_the_sigma_r_fold_is_a_splitting_only_with_its_anisotropic_remainder():
    r"""**M-6** — the #215 defect, and the degeneracy that hid it.

    The σ_r fold takes :math:`M = (L+C) - \Sigma_{s0}\,\mathbb{1}` — a
    removal that is **diagonal in angle**.  But the operator it removes from
    :math:`A` is :math:`S \supset \Sigma_{s0}P_{\text{iso}}`, the **isotropic
    projection**.  They are not the same operator, so the honest complement
    is

    .. math::

       N \;=\; M - A \;=\; S + B - \Sigma_{s0}\,\mathbb{1}
             \;=\; -\Sigma_{s0}\bigl(\mathbb{1} - P_{\text{iso}}\bigr) + B

    — the *anisotropic remainder*.  #215 shipped ``N = 0``.

    Three legs, and the third is the point:

    1. **honest** ``N`` satisfies the law (MEASURED 2.93e-17);
    2. ``N = 0`` REDs on an **anisotropic** state (MEASURED 5.43e-03
       relative, 2.63 absolute);
    3. ``N = 0`` is **invisible** on an angularly-flat state (MEASURED
       3.57e-18 — *machine zero*).

    Leg 3 is ``vv-principles`` **Mode 9** in closed form, and it is why #215
    survived: on an isotropic flux :math:`P_{\text{iso}}\psi = \psi`, so the
    two operators coincide EXACTLY and no tolerance, refinement, or regime
    change can expose the bug through an isotropic gate.  It ships as a
    permanent control: if it ever REDs, the "isotropic" state is no longer
    angularly flat and leg 2's claim to be *what caught the bug* is void.

    The fixture is vacuum-on-both-faces precisely so ``B ψ ≡ 0`` (MEASURED
    ``|Bψ|∞ = 0.0``) — the legs then isolate the projection mechanism and
    nothing else.
    """
    sn_mesh = isotropic_slab(c=0.9)
    record = record_for(sn_mesh, scattering_order=0)
    scattering = scattering_gain(record)
    boundary = record.explicit_gains[1]
    sigma_s0 = sigma_s0_times_identity(sn_mesh, scattering)
    folded_implicit = seedless_implicit(record) - sigma_s0
    honest_gains = (scattering, boundary, -sigma_s0)

    anisotropic = random_state(record, seed=4242)
    flat = random_state(record, seed=4242, angularly_flat=True)

    honest = reconstruction_residual(
        record, anisotropic, implicit=folded_implicit, gains=honest_gains,
    )
    if honest > 1e-13:
        pytest.fail(
            f"the σ_r fold's HONEST splitting violates A = M - N by "
            f"{honest:.3e} — the anisotropic remainder "
            f"-Σ_s0(I - P_iso) is no longer the correct complement."
        )

    caught = reconstruction_residual(
        record, anisotropic, implicit=folded_implicit, gains=(),
    )
    if caught < 1e-4:
        pytest.fail(
            f"the #215 bug (N := 0 under a σ_r fold) moved the splitting law "
            f"by only {caught:.3e} on an ANISOTROPIC state — this gate no "
            f"longer retro-catches a defect that shipped 46-56% flux errors."
        )

    hidden = reconstruction_residual(
        record, flat, implicit=folded_implicit, gains=(),
    )
    if hidden > 1e-12:
        pytest.fail(
            f"MODE-9 CONTROL BROKEN: the #215 bug measured {hidden:.3e} on an "
            f"angularly-flat state, but it must be machine-zero there "
            f"(P_iso psi == psi). The 'flat' probe is no longer flat, so the "
            f"anisotropic leg's claim to be what caught the bug is void."
        )


def test_the_two_sigma_s0_operators_are_indistinguishable_on_a_flat_flux():
    r"""The Mode-9 **mechanism**, asserted directly rather than inferred.

    :math:`\Sigma_{s0}P_{\text{iso}}` (what :math:`S` realizes) and
    :math:`\Sigma_{s0}\mathbb{1}` (what the σ_r sweep inverts) are different
    operators that agree **exactly** on an angularly-flat flux — MEASURED
    ``|Sψ|∞ = |Σ_s0 ψ|∞ = 2.045`` there, and materially apart otherwise.

    Gating the mechanism, not only the symptom, is what makes the M-6 control
    leg falsifiable: it pins *why* the isotropic state is blind, so a future
    reader cannot mistake leg 3's machine-zero for "the fold is fine".
    """
    sn_mesh = isotropic_slab(c=0.9)
    record = record_for(sn_mesh, scattering_order=0)
    scattering = scattering_gain(record)
    sigma_s0 = sigma_s0_times_identity(sn_mesh, scattering)

    flat = system_a(random_state(record, seed=4242, angularly_flat=True))
    anisotropic = system_a(random_state(record, seed=4242))

    np.testing.assert_allclose(
        scattering.apply(flat).to_flat(), sigma_s0.apply(flat).to_flat(),
        rtol=1e-13, atol=1e-14,
        err_msg=(
            "Sigma_s0*P_iso and Sigma_s0*I disagree on an ANGULARLY-FLAT "
            "flux, where P_iso is the identity — either the probe is not "
            "flat or S carries a channel beyond within-group P0."
        ),
    )
    spread = float(np.max(np.abs(
        scattering.apply(anisotropic).to_flat()
        - sigma_s0.apply(anisotropic).to_flat()
    )))
    if spread < 1e-3:
        pytest.fail(
            f"Sigma_s0*P_iso and Sigma_s0*I differ by only {spread:.3e} on an "
            f"ANISOTROPIC flux — the fixture is angularly degenerate and the "
            f"whole M-6 discriminator is vacuous."
        )

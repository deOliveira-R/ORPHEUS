r"""AC-b / AC-b′ — **one splitting per posed equation**, and the law
:math:`A = M - N` that makes a splitting a splitting.

Campaign: ``.claude/plans/operator_strategy_realization_campaign.md``.
Normative gate spec: ``.claude/plans/campaign_verification_plan.md`` §1.
Phase **P0** — gates only, no production change.

Why this file exists (the motivating defect, R7 — RESOLVED 2026-09-13)
======================================================================
Until the consumers campaign's step 2,
:class:`~orpheus.sn.coupled_system.WithinGroupSystem` welded two independent
stages into one record: the **posing** (``loss`` — what the physics *is*) and
the **splitting** (``implicit_operator`` / ``explicit_gains`` — which part is
solved implicitly).  Because the two shared a record, there was no named
boundary at which "strategy may enter" could be asserted — and because there
was no boundary, **a second splitting grew beside the first and the first went
silently stale**:

.. code-block:: text

   inner_schedule    driver actually ran                 record advertised
   ---------------   ---------------------------------   -----------------------------
   jacobi            StreamingCollisionOperator + SNBoundaryOperator   ← the same objects
   gauss_seidel      ScheduledInvertibleOperator + SNMaskedBoundaryOperator
                                                          StreamingCollisionOperator
                                                          + SNBoundaryOperator

``_select_si_splitting`` re-derived a splitting the record never heard
about.  Nothing was *numerically* wrong — both splittings were consistent —
but the record's claim was false, and a consumer that trusted
``record.implicit_operator`` (the spectral gate, an admission check, a
preconditioner) read an operator the solver did not run.

**The resolution (R-cc6 (i), 2026-09-13):** the record carries the loss and
its FACTORS (:class:`~orpheus.sn.coupled_system.SNLossFactors`, by role) and
no splitting; the splitting is the Strategy VALUE
:class:`~orpheus.sn.splitting.Splitting`, minted from the factors by the ONE
labelling site (:meth:`~orpheus.sn.splitting.Splitting.from_schedule`) and
consumed by the drivers AS IT IS — so "the driver runs the objects the value
advertises" holds by construction on every arm, and the strict-xfail marker
this file carried for R7 is gone (it XPASSed).  What keeps teeth is the LAW
``A = M − N`` per VALUE (AC-b′ below; the value's own
:meth:`~orpheus.sn.splitting.Splitting.law_residual`), and the pair of
schedules that must differ on the trace while agreeing on the limit
(``tests/sn/architecture/test_step2_terminal_object_anchors.py``).


Ordering constraint **O-1**, discharged
======================================
This gate was written **RED, before anything touched**
``WithinGroupSystem`` or ``_select_si_splitting`` (P0, ``9a546640``), and its
red row shipped as ``xfail(strict=True)`` so the fix could not land
unnoticed: it XPASSed at the step-2 split and the marker was deleted in the
same commit.  The AC-b rows below are now TAUTOLOGICAL by construction (the
driver is handed the value and runs it) — kept as the wiring claim
("the driver consumes the value it is handed"), with the teeth in AC-b′.

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
    sigma_s0_times_identity,
    slab_seedless,
    sphere_carrying,
    splitting_for,
    split_image,
    system_a,
)

pytestmark = pytest.mark.foundation

_SEED = 20260729


# ═════════════════════════════════════════════════════════════════════════
# AC-b — the driver runs the objects the record advertises
# ═════════════════════════════════════════════════════════════════════════


def _is_the_records_gain_rebound_on_the_iterate(driver_gain, record_gain) -> bool:
    """On a windowed (2-D Cartesian) mesh the driver's gain is the record's
    gain RE-BOUND to consume the moment iterate (CS4c step 5:
    ``on_moment_domain()`` — the same datum, the same faces, the same
    codomain; the domain's interior is the analysis face's codomain). Not a
    second splitting: nothing but the end the operand arrives on differs."""
    from orpheus.numerics.spaces.full_field_space import FullFieldSpace
    from orpheus.transport.operators.angular_lift import AngularLift

    if not isinstance(record_gain, AngularLift) or type(driver_gain) is not type(record_gain):
        return False
    # The datum + the two faces, by NAME — a rename must red this gate, not
    # degrade it (an `if hasattr` here would let `all([])` read True).
    datum_fields = ("transfer", "flux_analysis", "source_reconstruction")
    missing = [f for f in datum_fields if not hasattr(record_gain, f)]
    if missing:
        pytest.fail(
            f"the record's gain {type(record_gain).__name__} lacks "
            f"{missing} — the re-binding predicate names fields that no "
            f"longer exist; re-key it, do not weaken it."
        )
    same_datum = all(
        getattr(driver_gain, f) is getattr(record_gain, f) for f in datum_fields
    )
    trace = record_gain.domain.trace_space if isinstance(record_gain.domain, FullFieldSpace) else None
    return (
        same_datum
        and driver_gain.codomain is record_gain.codomain
        and trace is not None
        and driver_gain.domain == FullFieldSpace.from_blocks(
            record_gain.flux_analysis.codomain, trace,
        )
    )


@pytest.mark.parametrize(
    ("build_mesh", "inner_schedule"),
    [
        # CONTROL LEG: same geometry, same record, the other schedule.
        pytest.param(cart2d_seedless, "jacobi", id="cart2d-jacobi"),
        # The ex-RED: multi-D Cartesian + G-S was the only arm that re-split
        # (R7); since step 2 the value carries the G-S labelling itself.
        pytest.param(cart2d_seedless, "gauss_seidel", id="cart2d-gauss_seidel"),
        # The carrying arm labels Jacobi on BOTH rows — inner_schedule is
        # structurally inert there (G-S is multi-D Cartesian ⟹ seedless;
        # resolve_schedule falls back to Jacobi).
        pytest.param(sphere_carrying, "jacobi", id="sphere-jacobi"),
        pytest.param(sphere_carrying, "gauss_seidel", id="sphere-gauss_seidel"),
    ],
)
def test_driver_consumes_the_records_own_splitting(build_mesh, inner_schedule):
    r"""The SI driver iterates the **same objects** the Strategy value carries.

    Object identity (``is``), not value equality: two operators that happen
    to agree numerically today are still two operators, and the second one
    is where the drift lives.  This is the campaign's acceptance leg AC-b —
    since the step-2 split TRUE BY CONSTRUCTION (the driver is handed the
    value and runs its ``implicit``/``explicit``), so this row is the wiring
    claim, not the catcher; the catcher is AC-b′ (the law per value) and
    the two-schedules-differ row in the step-2 anchors.

    ⚠ The ONE sanctioned exception (CS4c step 5): on a WINDOWED mesh the
    iterate is the moment composite, and the record's angular-bound gains
    cannot read it — each binding acts through the body its ends select
    — so the driver consumes the record's own gains RE-BOUND on that end
    (``on_moment_domain()``): same datum, same faces, same codomain,
    checked field by field. The boundary gain is unchanged (it reads the
    trace) and stays the record's object.
    """
    sn_mesh = build_mesh()
    splitting = splitting_for(sn_mesh, inner_schedule)
    _si, driver_implicit, driver_gains, windowed = _within_group_si(
        splitting, sn_mesh, max_iter=2, tol=1e-10,
    )

    if driver_implicit is not splitting.implicit:
        pytest.fail(
            f"[{inner_schedule}] the driver's implicit operator is NOT the "
            f"value's: driver ran {type(driver_implicit).__name__}, the value "
            f"carries {type(splitting.implicit).__name__} — a splitting is "
            f"being re-derived behind the Strategy value (R7 again)."
        )
    advertised = tuple(map(id, splitting.explicit))
    if windowed:
        if len(driver_gains) != len(splitting.explicit) or not all(
            d is r or _is_the_records_gain_rebound_on_the_iterate(d, r)
            for d, r in zip(driver_gains, splitting.explicit)
        ):
            pytest.fail(
                f"[{inner_schedule}] windowed: the driver's gains are neither "
                f"the value's nor the value's re-bound on the moment "
                f"iterate: driver ran {[type(g).__name__ for g in driver_gains]}, "
                f"the value carries "
                f"{[type(g).__name__ for g in splitting.explicit]} (R7)."
            )
        return
    if tuple(map(id, driver_gains)) != advertised:
        pytest.fail(
            f"[{inner_schedule}] the driver's gains are NOT the value's: "
            f"driver ran {[type(g).__name__ for g in driver_gains]}, the value "
            f"carries {[type(g).__name__ for g in splitting.explicit]} "
            f"— a second splitting exists beside the value (R7 again)."
        )


@pytest.mark.parametrize(
    "build_mesh", [cart2d_seedless, sphere_carrying],
    ids=["cart2d", "sphere"],
)
def test_krylov_driver_consumes_the_records_own_splitting(build_mesh):
    r"""The Krylov driver is **green on both arms** — the twin is SI-specific.

    This is the second control leg, and it carries a positive architectural
    claim: ``_within_group_krylov`` is handed ``(splitting.implicit,
    *splitting.explicit)`` — the Jacobi value's members — and stores them
    verbatim, so R7 was never a property of "within-group solves" but
    specifically of the SI schedule path.  If this row ever REDs, the defect
    has spread.
    """
    sn_mesh = build_mesh()
    splitting = splitting_for(sn_mesh, "jacobi")
    state = random_state(splitting.system, seed=_SEED)
    krylov = _within_group_krylov(
        splitting.implicit, *splitting.explicit,
        n_dof=int(state.to_flat().size), max_iter=2, tol=1e-10,
    )
    if krylov.A is not splitting.implicit:
        pytest.fail(
            "the Krylov driver's operator is not the value's implicit "
            "operator — R7 has spread to the Krylov path."
        )
    if tuple(map(id, krylov.gains)) != tuple(map(id, splitting.explicit)):
        pytest.fail(
            "the Krylov driver's gains are not the record's — R7 has spread "
            "to the Krylov path."
        )


def test_a_slab_hides_r7_the_documented_trap():
    r"""A 1-D slab labels JACOBI under ``gauss_seidel`` — and that is a TRAP.

    :func:`~orpheus.sn.splitting.resolve_schedule` falls back to Jacobi unless
    ``is_cartesian and not is_1d``, so a slab exercises the *control* arm
    under the *other* arm's name.  The campaign's first probe of R7 used a
    slab, measured ``True`` on both rows, and concluded there was no twin.

    Committing the trap as a named row makes that mistake unrepeatable: this
    test asserts the slab is green **and** says why, so a future reader who
    reaches for a slab fixture meets the explanation first.  It also has real
    teeth — if the 1-D fallback is ever removed, this REDs and the reader is
    told exactly which invariant moved.
    """
    sn_mesh = slab_seedless()
    splitting = splitting_for(sn_mesh, "gauss_seidel")
    if splitting.schedule.is_sequenced or (
        splitting.implicit is not splitting.system.factors.streaming_collision
    ) or splitting.explicit[-1] is not splitting.system.factors.boundary:
        pytest.fail(
            "the 1-D slab no longer falls back to Jacobi under "
            "inner_schedule='gauss_seidel' — the R7 observability precondition "
            "changed; re-read resolve_schedule's geometry gate."
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
    splitting = splitting_for(sn_mesh, "jacobi")
    state = random_state(splitting.system, seed=_SEED)
    loss_image = splitting.system.loss.apply(state).to_flat()
    split = split_image(splitting, state)

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
    splitting = splitting_for(sn_mesh, "jacobi")
    state = random_state(splitting.system, seed=_SEED)
    defect = reconstruction_residual(splitting, state, gains=())
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
    splitting = splitting_for(sn_mesh, "jacobi")
    state = random_state(splitting.system, seed=_SEED)
    flipped = (-splitting.explicit[0], *splitting.explicit[1:])
    defect = reconstruction_residual(splitting, state, gains=flipped)
    if defect < 1e-3:
        pytest.fail(
            f"flipping the first gain's sign moved the splitting law by only "
            f"{defect:.3e} — the law cannot see a gain-sign convention drift."
        )


# ═════════════════════════════════════════════════════════════════════════
# M-6 — the exact σ_r fold, with its Mode-9 control leg
# ═════════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-070")
def test_the_sigma_r_fold_is_a_splitting_only_with_its_anisotropic_remainder():
    r"""**M-6** — the #215 / ERR-070 defect, and the degeneracy that hid it.

    ``catches("ERR-070")`` is earned, not asserted: the catalog entry names
    the bug as *"treats* :math:`\Sigma_{s0}P_{\text{iso}}` *as*
    :math:`\Sigma_{s0}\mathbb{1}`\ *"* with the difference operator
    :math:`\sigma_{s0}(\mathbb{1}-P_{\text{iso}})` annihilating exactly the
    isotropic subspace — and this test constructs precisely that splitting
    and measures it RED (5.43e-03) on an anisotropic state, invisible
    (3.57e-18) on a flat one.

    It is a **third, earlier** catcher, not a duplicate of the two on record.
    The existing pair are a *value* gate (the DSA fixed point shifts 43 %,
    ``test_dsa_rate.py::TestSigmaRFoldCaught``) and a *structural fence* (an
    AST sweep of the foldable accessors' consumers,
    ``TestD10RoutingSentinel``).  This one is *algebraic*: it reds the moment
    the splitting is **constructed**, before any solve runs — the cheapest of
    the three, and the layer at which the campaign intends to make the bug
    unspellable rather than merely detectable.

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
    splitting = splitting_for(sn_mesh, "jacobi", scattering_order=0)
    factors = splitting.system.factors
    scattering, boundary = factors.scattering, factors.boundary
    sigma_s0 = sigma_s0_times_identity(sn_mesh, scattering)
    folded_implicit = factors.streaming_collision - sigma_s0
    honest_gains = (scattering, boundary, -sigma_s0)

    anisotropic = random_state(splitting.system, seed=4242)
    flat = random_state(splitting.system, seed=4242, angularly_flat=True)

    honest = reconstruction_residual(
        splitting, anisotropic, implicit=folded_implicit, gains=honest_gains,
    )
    if honest > 1e-13:
        pytest.fail(
            f"the σ_r fold's HONEST splitting violates A = M - N by "
            f"{honest:.3e} — the anisotropic remainder "
            f"-Σ_s0(I - P_iso) is no longer the correct complement."
        )

    caught = reconstruction_residual(
        splitting, anisotropic, implicit=folded_implicit, gains=(),
    )
    if caught < 1e-4:
        pytest.fail(
            f"the #215 bug (N := 0 under a σ_r fold) moved the splitting law "
            f"by only {caught:.3e} on an ANISOTROPIC state — this gate no "
            f"longer retro-catches a defect that shipped 46-56% flux errors."
        )

    hidden = reconstruction_residual(
        splitting, flat, implicit=folded_implicit, gains=(),
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

    **Deliberately carries NO** ``catches("ERR-070")``.  This test asserts a
    property that is true both before and *after* the bug — it characterises
    the degeneracy, it does not detect the defect.  A marker here would be a
    phantom coverage edge: it would inflate ERR-070's catcher count with a
    test that stays green under the exact mutation it appears to cover.
    """
    sn_mesh = isotropic_slab(c=0.9)
    record = record_for(sn_mesh, scattering_order=0)
    scattering = record.factors.scattering
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

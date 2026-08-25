r"""Boundary Gauss-Seidel is a COHERENT splitting of ``A``, on a singular box.

#344's disposition turns on a distinction that looks identical from a distance.
On the default ``k_inf`` lattice ``A = L + C - S - B`` is exactly singular
(:ref:`sn-loss-kernel-gauge`), and boundary Gauss-Seidel returns a boundary
trace several percent away from the analytic answer while Jacobi returns it
exactly. Two mechanisms produce that signature:

* **manifold selection** — both schedules are honest splittings of the SAME
  ``A``, and a singular ``A`` has a solution *manifold*, so the oblique
  projector SI freezes differs by splitting and the two land on different
  members. Benign; repaired at the exit by the gauge.
* **an INCOHERENT schedule** — ``M - N != A``, e.g. ERR-056's reflect-after-the-
  FIRST-outflowing-group. A correctness bug.

`vv` Mode 9's discriminator is that an incoherent splitting moves the **bulk**
too, and that with the kernel REMOVED the two schedules must agree on the
**trace** as well. This module gates all three checks — ``M - N == A``, the
kernel-free trace agreement over four independent kernel-removal mechanisms,
and the ERR-056 positive control that gives the second one teeth — plus the two
facts that make the disposition legible: neither convergence functional can see
a ``ker A`` shift, and the driver's manifold selection is first-order in ``h``.

⛔ **Read every test name here for its TIER.** All but one measure the SI
**driver**, which does not gauge; the exception says so in its name
(:func:`test_the_PUBLIC_ENTRY_returns_the_SAME_trace_under_BOTH_schedules`) and
is the repaired counterpart of the driver row directly above it. Two of these
gates arrived from ``derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py``
(deleted in the same change that closed #344; it lives in ``git log``)
under names that read as claims about what a USER receives — they were true when
written and Step 5 of the campaign repealed them (`plan-authoring` §3, *"a fact
can die by being FIXED"*). They are re-titled here to name the tier they
measure, and the driver/entry pair is now adjacent so neither can be read
alone.

Adjacent surfaces: :mod:`tests.sn.solve.test_every_entry_gauges_its_trace`
(that the gauge FIRES at every entry, and the physical witness it repairs) and
:mod:`tests.sn.operators.test_loss_nullspace_reflective_box` (the structural
facts about :math:`\ker A` itself).

Promoted 2026-08-15 (GitHub #344); ``pyproject.toml``'s ``testpaths = ["tests"]``
never ran the diagnostic, so this is net-new executed coverage.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.solver import _balance_projection
from tests.sn._singular_loss_box import (
    REFLECTIVE as _R,
    VACUUM as _V,
    absorber,
    assemble,
    both_drivers,
    build,
    drive,
    drive_recorded,
    gains_matvec,
    isotropic_source,
    ld_reflective_box,
    loss_matvec,
    null_basis,
    scatterer,
    select_splitting,
    uniform_source_fixture,
)

pytestmark = pytest.mark.foundation

#: Smallest UNEQUAL-cell singular box — see the operator-tier module for why the
#: shrink from ``(3,4)`` is licensed (``dim ker A = 12`` either way).
_CELLS = (2, 3)


# ─────────────────────────────────────────────────────────────────────
# 1. both schedules really are splittings of the SAME operator
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    "label,cells,bcs,mixture",
    [
        ("all reflective, c=0.9", _CELLS, [(_R, _R)] * 2, scatterer(2, c=0.9)),
        ("x vacuum, c=0.5", _CELLS, [(_V, _V), (_R, _R)], scatterer(2, c=0.5)),
    ],
)
def test_both_schedules_are_splittings_of_the_SAME_A(
        label, cells, bcs, mixture):
    r"""``M - N == A`` BIT-EXACTLY, for Jacobi and for boundary Gauss-Seidel.

    The definition of a splitting, and the premise of #344's whole determinism
    verdict: if it failed, the deviation would not be manifold selection, it
    would be a wrong operator.

    The pieces of this are pinned elsewhere and the conjunction is not.
    ``test_w2_split_exactness`` gives ``B == B_lower + B_upper`` per face,
    ``test_factory_returns_reified_pair`` gives the shape of what
    :func:`~orpheus.sn.solver._select_si_splitting` returns, and
    ``test_w2_round_trip_machine_precision`` pins ``M.apply`` against
    ``M.inverse()``. This asserts the resulting identity **at the level of the
    whole operator**, exhaustively — one unit probe per degree of freedom, so
    every entry of ``M - N - A`` is checked, not a random subspace.

    ``[M]`` bit-exact ``0.0`` on 6 configurations (absorber / c=0.5 / c=0.9 /
    x-vacuum, at ``(2,3)`` and ``(3,4)``). It is exact rather than
    round-off-close because the split writes **disjoint rows**: on any row that
    carries a boundary value, exactly one of ``B_lower``/``B_upper`` is
    non-zero, so no addition is reordered.

    ⚠ **The absorber fixture the diagnostic used for a third row is dropped
    deliberately**: ``[M]`` ``|S x| = 0.000000e+00`` on it, so ``S`` is the ZERO
    operator and the row cannot see where the scattering gain lands. Both rows
    here carry ``S != 0`` (``[M]`` ``6.741510e+00`` at ``c = 0.9``).

    ``[M]`` mutation: drop ``B_upper`` from the Gauss-Seidel arm's gains —
    ``M - N`` then differs from ``A`` by ``B_upper`` and this reds while the
    per-face partition gate over in ``test_gauss_seidel_reification.py`` stays
    green, because ``B`` itself is untouched.
    """
    sn_mesh, system, template = build(cells, bcs, mixture)
    dense_a = assemble(loss_matvec(system), template)
    scattering, boundary = system.explicit_gains

    probe = np.random.default_rng(0).standard_normal(template.to_flat().size)
    activation = float(np.linalg.norm(scattering.apply(
        type(template).from_flat(probe, template)).to_flat()))
    assert activation > 1e-3, (
        f"{label}: the scattering gain is inert on this mixture "
        f"(|S x| = {activation:.3e}) — the row cannot see an error in where S "
        f"lands in the splitting"
    )

    for schedule in ("jacobi", "gauss_seidel"):
        base, gains = select_splitting(system, sn_mesh, schedule)
        dense_m = assemble(lambda x: base.apply(x), template)
        dense_n = assemble(gains_matvec(gains), template)
        assert np.array_equal(dense_m - dense_n, dense_a), (
            f"{label}/{schedule}: ||M - N - A||/||A|| = "
            f"{np.linalg.norm(dense_m - dense_n - dense_a) / np.linalg.norm(dense_a):.6e} "
            f"— the schedule is NOT a splitting of A, so #344's DETERMINISM "
            f"verdict reverses to a correctness bug. (A non-zero that is "
            f"merely round-off would still be a finding: the identity is "
            f"exact because the row supports are disjoint.)"
        )


def test_the_gauss_seidel_inverse_is_exact_on_the_DRIVER_RHS_subspace():
    r"""``M_GS.inverse()`` is a SUBSPACE inverse — exact exactly where it is used.

    ⛔ It is **not** a full-space inverse: the scheduled walk re-derives the
    outflow-definition rows, so it realizes :math:`M^{-1}` only on
    ``{y : y.outflow-rows == 0}`` (ERR-071, pinned by
    ``test_off_domain_outflow_rhs_is_not_completed_yet``). What that gate
    measures is the size of the gap; what it *states in prose* — that every
    production right-hand side lies inside the subspace — is what this gate
    ASSERTS, on the actual driver rhs ``q + S psi + B_upper psi`` for random
    ``psi``.

    That first assertion is the load-bearing one. It is the tripwire for a
    change that gives the driver rhs outflow-trace content: the inverse would
    then be applied off its domain, silently, and #344's coherence conclusion
    would no longer follow. ``[M]`` the rhs's outflow rows are ``0.0`` exactly,
    and ``max|M y - r| / max|r|`` reads ``3.6e-16 .. 5.5e-16``.

    ⚠ HAZARD this pins: a **Krylov preconditioner** built from
    ``M_GS.inverse()`` would feed it arbitrary vectors and is outside the
    contract. ``M_J.inverse()`` is exact in full space.

    ``[M]`` mutation: seed the rhs's outflow rows with the diagnostic value
    ``1.0`` — the first assertion reds, and the round-trip defect on the
    lower-coupled inflow rows rises to O(1).
    """
    sn_mesh, system, template = build(_CELLS, [(_R, _R)] * 2, absorber(2))
    n_dof = template.to_flat().size
    n_bulk = template.interior.values.size

    omega_dot_n = np.asarray(sn_mesh.angular_trace.omega_dot_n)
    layout = sn_mesh.angular_trace.layout
    outflow = np.zeros(n_dof, dtype=bool)
    for index, face in enumerate(layout.faces):
        slot = layout.faces[face]
        per_ordinate = int(np.prod(slot.shape[1:]))
        for ordinate in range(slot.shape[0]):
            if omega_dot_n[index, ordinate] > 0:
                start = n_bulk + slot.offset + ordinate * per_ordinate
                outflow[start:start + per_ordinate] = True
    assert outflow.any(), "no outflow rows — the subspace claim is vacuous"

    base, gains = select_splitting(system, sn_mesh, "gauss_seidel")
    dense_m = assemble(lambda x: base.apply(x), template)
    dense_n = assemble(gains_matvec(gains), template)
    inverse = base.inverse()
    zero = type(template).from_flat(np.zeros(n_dof), template)
    source = isotropic_source(sn_mesh, template).to_flat()

    generator = np.random.default_rng(1)
    for draw in range(3):
        rhs = source + dense_n @ generator.standard_normal(n_dof)
        worst_outflow = float(np.max(np.abs(rhs[outflow])))
        assert worst_outflow == 0.0, (
            f"draw {draw}: the driver rhs grew outflow-trace content "
            f"(max |rhs| there = {worst_outflow:.6e}) — M_GS^-1 is only an "
            f"inverse where that is zero, so the subspace argument behind "
            f"#344's coherence verdict must be re-derived"
        )
        solved = inverse.apply(
            type(template).from_flat(rhs, template), initial_guess=zero,
        ).to_flat()
        relative = float(np.max(np.abs(dense_m @ solved - rhs))
                         / np.max(np.abs(rhs)))
        assert relative < 1e-12, (
            f"draw {draw}: M_GS^-1 is not an inverse on the driver rhs "
            f"subspace ({relative:.6e})"
        )


# ─────────────────────────────────────────────────────────────────────
# 2. ⭐ the coherence gate, and the mutation that gives it teeth
# ─────────────────────────────────────────────────────────────────────
#: Four INDEPENDENT ways to make ``dim ker A == 0``, so the coherence verdict
#: cannot be an artefact of HOW the kernel was removed (`vv` §structural
#: independence). Three break the GEOMETRY, one breaks the CLOSURE.
_KERNEL_FREE = {
    "x vacuum": ((3, 4), [(_V, _V), (_R, _R)], None),
    "xmin reflective / xmax vacuum": ((3, 4), [(_R, _V), (_R, _R)], None),
    "LD on the ALL-reflective box": (_CELLS, [(_R, _R)] * 2, "LD"),
    "d3, one reflective pair": ((2, 2, 2), [(_V, _V), (_V, _V), (_R, _R)],
                                None),
}


def _kernel_free(label):
    """``(sn_mesh, system, template, null basis)`` for one removal mechanism.

    The LD row routes through :func:`~tests.sn._singular_loss_box.ld_reflective_box`
    so its ``[M]`` 10.8 s assembly is shared with the operator-tier gate that
    owns the same measurement's other claim.
    """
    cells, bcs, tag = _KERNEL_FREE[label]
    if tag == "LD":
        sn_mesh, system, template, _dense, basis, _s = ld_reflective_box(cells)
        return sn_mesh, system, template, basis
    sn_mesh, system, template = build(cells, bcs, scatterer())
    return (sn_mesh, system, template,
            null_basis(assemble(loss_matvec(system), template))[0])


@pytest.mark.parametrize("label", list(_KERNEL_FREE))
def test_kernel_free_configs_give_the_SAME_TRACE_under_both_schedules(label):
    r"""⭐ With ``dim ker A == 0`` there is nowhere to hide — trace AND bulk.

    This is the gate that separates *manifold selection* from *incoherence*.
    Remove the kernel and the two schedules must agree on **everything**; if
    boundary Gauss-Seidel were reflecting a face whose trace is inconsistent
    with the bulk state (ERR-056's failure mode), at least one of these four
    would disagree. ``[M]`` trace ``<= 3.7e-13``, bulk ``<= 2.4e-13``.

    ⭐ The **trace** comparison is what is new. The tree's two existing
    Mode-9 ERR-056 gates — ``test_w2_fixed_point_equivalence_diagonal_cubature``
    and ``test_d3_gauss_seidel_jacobi_fixed_point_invariance`` — both compare
    ``scalar_flux``/``keff``, i.e. the BULK, which is precisely the half that
    manifold selection leaves alone. A trace-level disagreement on a
    kernel-free box is invisible to them.

    The four rows remove the kernel four different ways, and the LD row is the
    only one that does it by changing the **closure** rather than the geometry
    — the others all open a vacuum face, so without it the conclusion would be
    "vacuum faces make the schedules agree". Its precondition is asserted
    (``dim ker A == 0``) rather than assumed, because
    :func:`~orpheus.sn.operators.loss_kernel_gauge.gauge_freedom` reports
    UNDETERMINED for LD at ``ndim = 3`` and a future closure could do the same
    at ``ndim = 2``.

    Teeth: :func:`test_the_err056_first_group_reflect_mutation_reddens_this`.
    """
    sn_mesh, system, template, basis = _kernel_free(label)
    assert basis.shape[1] == 0, (
        f"{label} is NOT kernel-free (dim ker A = {basis.shape[1]}) — it "
        f"cannot serve as the coherence control, because a disagreement would "
        f"then have somewhere benign to hide"
    )

    n_bulk = template.interior.values.size
    iterates = both_drivers(sn_mesh, system, template)
    difference = iterates["gauss_seidel"] - iterates["jacobi"]
    reference = iterates["jacobi"]
    trace = float(np.max(np.abs(difference[n_bulk:]))
                  / np.max(np.abs(reference[n_bulk:])))
    bulk = float(np.max(np.abs(difference[:n_bulk]))
                 / np.max(np.abs(reference[:n_bulk])))
    assert trace < 1e-10 and bulk < 1e-10, (
        f"{label}: the two schedules disagree with NO kernel to hide in — "
        f"trace {trace:.6e}, bulk {bulk:.6e}. Boundary Gauss-Seidel is "
        f"INCOHERENT and #344 reverses from a determinism disposition to a "
        f"correctness bug"
    )


@pytest.mark.catches("ERR-056")
def test_the_err056_first_group_reflect_mutation_reddens_this(monkeypatch):
    r"""⭐ The positive control (`vv` #17) for the gate above.

    Reflect each face after the **FIRST** outflowing octant group instead of the
    LAST — exactly the failure
    :meth:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule.gauss_seidel`'s
    deferred-reflect rule exists to prevent, on a ``level_symmetric`` cubature
    where faces are shared across groups. On a kernel-free configuration this
    must move both the trace and the **bulk**; if it did not, every coherence
    green above would carry no information.

    ``[M]`` trace ``1.0000e+00`` / bulk ``7.17e-01`` against the baseline
    ``7.81e-13`` / ``1.67e-13`` — twelve orders of dynamic range.

    It also carries ``catches("ERR-056")`` on its own account, and that is
    mutation-verified rather than assumed: if production ever shipped the
    first-group rule, the baseline leg (``base_trace < 1e-10``) and the
    schedules-really-differ leg would both red here. The tree's other
    ``catches("ERR-056")`` gates compare bulk quantities only.

    ⚠ The mutation goes in through ``monkeypatch`` — belt (its teardown) and
    braces (an explicit ``undo()`` inside the body, so the revert can be
    *checked* by re-solving and requiring bit-identity with the baseline). A
    leaked mutation would poison every later test in the session, and no
    assertion in those tests would attribute the damage back here.
    """
    from orpheus.sn.loss_representation import sweep_schedule as schedules

    label = "x vacuum"
    sn_mesh, system, template, _basis = _kernel_free(label)
    n_bulk = template.interior.values.size
    baseline = both_drivers(sn_mesh, system, template)
    baseline_trace = float(
        np.max(np.abs(baseline["gauss_seidel"][n_bulk:]
                      - baseline["jacobi"][n_bulk:]))
        / np.max(np.abs(baseline["jacobi"][n_bulk:]))
    )
    assert baseline_trace < 1e-10, (
        f"the UNMUTATED baseline already disagrees ({baseline_trace:.6e}) — "
        f"production may be shipping the ERR-056 rule"
    )

    def first_group_gauss_seidel(cls, ndim, octants, reflective):
        ordered, by_label = [], {}
        for entry in octants:
            sweep = schedules._octant_sweep(entry, ndim)
            if sweep.label not in by_label:
                by_label[sweep.label] = []
                ordered.append(sweep.label)
            by_label[sweep.label].append(sweep)
        first = {}
        for index, name in enumerate(ordered):
            for face in schedules._outgoing_faces(name):
                if face in reflective and face not in first:
                    first[face] = index          # ⛔ the ERR-056 mutation
        by_group = {index: [] for index in range(len(ordered))}
        for face, index in first.items():
            by_group[index].append(face)
        return cls(
            groups=tuple(
                schedules.OctantSweepGroup(
                    sweeps=tuple(by_label[name]),
                    reflect_faces=tuple(sorted(by_group[index])),
                )
                for index, name in enumerate(ordered)
            ),
            kind="gauss_seidel",
        )

    shipped = tuple(group.reflect_faces for group
                    in schedules.SweepSchedule.gauss_seidel(
                        sn_mesh.ndim, sn_mesh.quad.octants,
                        schedules.reflective_faces(sn_mesh)).groups)
    mutated = tuple(group.reflect_faces for group in first_group_gauss_seidel(
        schedules.SweepSchedule, sn_mesh.ndim, sn_mesh.quad.octants,
        schedules.reflective_faces(sn_mesh)).groups)
    assert shipped != mutated, (
        f"the ERR-056 mutation did not change the schedule ({shipped}) — the "
        f"control is inert, so it certifies nothing about the gate above"
    )

    monkeypatch.setattr(schedules.SweepSchedule, "gauss_seidel",
                        classmethod(first_group_gauss_seidel))
    broken = both_drivers(sn_mesh, system, template)
    monkeypatch.undo()

    broken_trace = float(
        np.max(np.abs(broken["gauss_seidel"][n_bulk:]
                      - broken["jacobi"][n_bulk:]))
        / np.max(np.abs(broken["jacobi"][n_bulk:])))
    broken_bulk = float(
        np.max(np.abs(broken["gauss_seidel"][:n_bulk]
                      - broken["jacobi"][:n_bulk]))
        / np.max(np.abs(broken["jacobi"][:n_bulk])))
    assert broken_trace > 1e-2 and broken_bulk > 1e-2, (
        f"the ERR-056 mutation did NOT redden the coherence comparison "
        f"(trace {broken_trace:.6e}, bulk {broken_bulk:.6e}) — the gate above "
        f"is a dud and every coherence green in this module is uninformative"
    )

    restored = both_drivers(sn_mesh, system, template)
    assert np.array_equal(restored["gauss_seidel"],
                          baseline["gauss_seidel"]), (
        "the schedule monkeypatch did not revert cleanly — every later test "
        "in this session is running against a mutated production schedule"
    )


# ─────────────────────────────────────────────────────────────────────
# 3. neither convergence functional can see the kernel
# ─────────────────────────────────────────────────────────────────────
def _balance(flat, template, sn_mesh) -> float:
    """The production balance projection's norm on a flat state."""
    return float(np.linalg.norm(np.asarray(_balance_projection(
        type(template).from_flat(flat, template), sn_mesh=sn_mesh))))


@pytest.mark.verifies("sn-loss-kernel-gauge-projection")
def test_two_cold_starts_a_KERNEL_APART_report_the_SAME_certificate():
    r"""⭐ The reason the defect went unseen: the certificates cannot show it.

    Run the production driver twice, from cold starts differing **only inside**
    :math:`\ker A`. Every convergence quantity a caller can read agrees, and the
    returned trace does not. ``[M]`` iteration count **344 both**; final
    residual ``9.028098e-14`` vs ``9.022488e-14`` (``6.2e-04`` apart, both under
    the ``1e-13`` request); the production balance projection ``2.795085``
    **bit-identically** on each; bulk agreeing to ``8.4e-16`` — and the traces
    **11.26 %** apart.

    ⚠ **This gate was re-posed after promotion, and the original could not have
    failed.** The diagnostic asserted ``||A v|| ~ 0`` for ``v`` taken from a
    dense SVD's null space — which is what an SVD null space *is*. Both
    "blindness" legs were theorems about the factorisation, not measurements of
    the solver; only the two anchor legs and the positive control had teeth.
    Here the measurand is what
    :class:`~orpheus.numerics.iteration.SourceIteration` actually reports, so a
    stopping rule that DID see the kernel would red.

    The mechanism, and why the traces stay apart rather than reconverging:
    ``v`` in :math:`\ker A = \ker(M - N)` is a **fixed direction** of the
    iteration operator :math:`G = M^{-1}N` (:math:`Gv = v`), so SI preserves the
    kernel component of its initial guess exactly and lands on
    :math:`\psi^* + v`. ``[M]`` the difference IS ``v``, to ``2.3e-14``.

    ⭐ **The control** (`vv` #19 — a positive reading alone cannot tell
    *blind* from *loaded*) is the same experiment on a kernel-FREE box, where
    :math:`G` has no fixed direction and a perturbation of the same size is
    damped: ``[M]`` the two traces agree to ``2.6e-15``. Without it, "the
    certificates agree and the traces differ" is also what a broken comparison
    produces.

    Jacobi is used deliberately: from a zero start it lands ON the canonical
    member, so the only thing separating the two runs is the injected kernel
    vector. The splitting's own manifold selection is a different claim, gated
    by :func:`test_the_UNGAUGED_SI_driver_returns_a_splitting_dependent_trace`.
    """
    sn_mesh, system, template, n_dof, n_bulk, source, _exact = (
        uniform_source_fixture(_CELLS))
    basis, _singular = null_basis(assemble(loss_matvec(system), template))
    assert basis.shape[1] > 0, "fixture is no longer singular"

    generator = np.random.default_rng(0)
    coefficients = generator.standard_normal(basis.shape[1])
    kernel_start = basis @ (coefficients / np.linalg.norm(coefficients))
    assert np.max(np.abs(kernel_start[:n_bulk])) < 1e-14, (
        "the null basis has bulk content — see the operator-tier PURE TRACE "
        "gate; this experiment would then not be a pure trace perturbation"
    )
    # Scale to a trace perturbation of ~11 %, the magnitude #344 reports.
    reference = float(np.max(np.abs(
        drive(system, sn_mesh, template, source, "jacobi")[n_bulk:])))
    kernel_start *= 0.1126 * reference / float(
        np.max(np.abs(kernel_start[n_bulk:])))

    cold, cold_record = drive_recorded(
        system, sn_mesh, template, source, "jacobi")
    shifted, shifted_record = drive_recorded(
        system, sn_mesh, template, source, "jacobi", initial=kernel_start)

    def final_residual(record):
        criterion = record.binding_criterion
        assert criterion is not None and criterion.trajectory, (
            "the driver reported no stopping trajectory — there is then no "
            "certificate to compare and this gate has no measurand"
        )
        return float(criterion.trajectory[-1])

    assert cold_record.n_iterations == shifted_record.n_iterations, (
        f"the iteration COUNT saw the kernel shift: "
        f"{cold_record.n_iterations} vs {shifted_record.n_iterations}"
    )
    # [M] 6.2e-04 apart; a functional that COULD see an 11 % trace shift moves
    # by O(1) relative, so the 5 % bar is ~80x above the noise and ~20x below
    # any real signal.
    assert final_residual(cold_record) == pytest.approx(
        final_residual(shifted_record), rel=5e-2), (
        f"the SI residual saw the kernel shift: "
        f"{final_residual(cold_record):.6e} vs "
        f"{final_residual(shifted_record):.6e}"
    )
    assert _balance(cold, template, sn_mesh) == pytest.approx(
        _balance(shifted, template, sn_mesh), abs=1e-11), (
        f"the balance projection saw the kernel shift: "
        f"{_balance(cold, template, sn_mesh):.6e} vs "
        f"{_balance(shifted, template, sn_mesh):.6e}"
    )
    bulk = float(np.max(np.abs(cold[:n_bulk] - shifted[:n_bulk]))
                 / np.max(np.abs(cold[:n_bulk])))
    assert bulk < 1e-11, (
        f"the BULK moved between the two cold starts ({bulk:.3e}) — the kernel "
        f"is pure-trace, so this would mean the perturbation left it"
    )

    trace = float(np.max(np.abs(cold[n_bulk:] - shifted[n_bulk:]))
                  / np.max(np.abs(cold[n_bulk:])))
    assert trace > 1e-2, (
        f"the two cold starts returned the SAME trace ({trace:.3e}) — every "
        f"agreement above is then vacuous. SI preserves a ker A component "
        f"exactly (Gv = v), so this can only mean the fixture stopped being "
        f"singular or the perturbation was not in ker A"
    )
    residue = float(np.linalg.norm((shifted - cold) - kernel_start)
                    / np.linalg.norm(kernel_start))
    assert residue < 1e-9, (
        f"the two answers differ by something other than the injected kernel "
        f"vector (relative residue {residue:.3e}) — the mechanism is not "
        f"manifold selection"
    )

    # ⭐ CONTROL — the same experiment where there is no manifold to select on.
    free_mesh, free_system, free_template = build(
        (3, 4), [(_V, _V), (_R, _R)], scatterer())
    free_basis, _s = null_basis(
        assemble(loss_matvec(free_system), free_template))
    assert free_basis.shape[1] == 0, "the control box is not kernel-free"
    free_source = isotropic_source(free_mesh, free_template)
    free_bulk = free_template.interior.values.size
    free_cold = drive(free_system, free_mesh, free_template, free_source,
                      "jacobi")
    perturbation = np.zeros(free_template.to_flat().size)
    perturbation[free_bulk:] = 0.1126 * float(
        np.max(np.abs(free_cold[free_bulk:]))) * np.random.default_rng(
            1).standard_normal(perturbation.size - free_bulk)
    free_shifted = drive(free_system, free_mesh, free_template, free_source,
                         "jacobi", initial=perturbation)
    free_trace = float(np.max(np.abs(free_cold[free_bulk:]
                                     - free_shifted[free_bulk:]))
                       / np.max(np.abs(free_cold[free_bulk:])))
    assert free_trace < 1e-10, (
        f"a NON-singular box also returned a start-dependent trace "
        f"({free_trace:.3e}) — the experiment is measuring something other "
        f"than the kernel, and the readings above carry no information"
    )


# ─────────────────────────────────────────────────────────────────────
# 4. the driver's manifold selection — CHARACTERIZATION, at the driver tier
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("cells", [(3, 2), (3, 3)])
def test_the_UNGAUGED_SI_driver_returns_a_splitting_dependent_trace(cells):
    r"""A splitting cannot change the equation — it selects a manifold MEMBER.

    ⛔ **Tier.** This runs the SI **driver** (:func:`tests.sn._singular_loss_box.drive`),
    which does not gauge. Every public entry does
    (:func:`~orpheus.sn.solver._exit_gauge_trace`, five call sites), so **the
    trace a user receives no longer has this property** — that is asserted in
    :mod:`tests.sn.solve.test_every_entry_gauges_its_trace`
    (``test_every_exercisable_entry_reports_a_gauge_correction``,
    ``test_the_spurious_TANGENTIAL_current_along_a_mirror_is_gone``). The
    diagnostic this came from called it *"the returned boundary trace depends on
    the SPLITTING"*, which was true when written and repealed by Step 5 of the
    campaign.

    What it still gates is the *mechanism*, and three legs are needed to
    identify it:

    * Jacobi lands on the analytic member (``[M]`` ``1.4e-12``) and
      Gauss-Seidel does not (``[M]`` ``1.05e-01`` / ``9.80e-02``);
    * the difference is **in** :math:`\ker A` — ``[M]``
      :math:`\lVert Ad\rVert/\lVert d\rVert \le 1.2\times10^{-12}`, computed
      with the production matvec, so this leg owes nothing to the gauge. Change
      the mechanism and it reds even if the magnitudes survive;
    * ⭐ the **BULK is untouched under both schedules** (``[M]``
      ``<= 1.2e-13``). That is the `vv` Mode-9 discriminator that separates this
      from ERR-056-class incoherence, which moves the bulk by 0.39–0.80. Without
      it, "the trace differs" is compatible with a real bug.

    `vv` Mode 9 says a splitting must not move the fixed POINT. Here there is no
    fixed point — there is a manifold — so that guarantee is void, and this
    records the void rather than reading it as a failure.
    """
    sn_mesh, system, template, _n, n_bulk, source, exact = (
        uniform_source_fixture(cells))
    apply_loss = loss_matvec(system)

    readings = {}
    for schedule in ("gauss_seidel", "jacobi"):
        iterate = drive(system, sn_mesh, template, source, schedule)
        difference = iterate - exact
        readings[schedule] = (
            float(np.max(np.abs(difference[n_bulk:] / exact[n_bulk:]))),
            float(np.max(np.abs(difference[:n_bulk] / exact[:n_bulk]))),
            float(np.linalg.norm(apply_loss(
                type(template).from_flat(difference, template)).to_flat())
                / max(np.linalg.norm(difference), 1e-300)),
        )
    gs_trace, gs_bulk, gs_residual = readings["gauss_seidel"]
    jacobi_trace, jacobi_bulk, _ = readings["jacobi"]

    assert jacobi_trace < 1e-10, (
        f"jacobi no longer returns the analytic trace ({jacobi_trace:.6e}) — "
        f"the contrast this row is built on is gone"
    )
    assert gs_trace > 1e-2, (
        f"the gauss_seidel driver's trace deviation is {gs_trace:.6e}; if the "
        f"SPLITTING was changed so it now lands on the canonical member, this "
        f"characterization is obsolete — delete it rather than re-baselining"
    )
    assert gs_residual < 1e-9, (
        f"the gauss_seidel deviation is NOT a ker A component "
        f"(||A d||/||d|| = {gs_residual:.3e}) — it is then not manifold "
        f"selection, and #344's disposition does not cover it"
    )
    assert gs_bulk < 1e-10 and jacobi_bulk < 1e-10, (
        f"the BULK moved between schedules (gauss_seidel {gs_bulk:.3e}, "
        f"jacobi {jacobi_bulk:.3e}) — that is the ERR-056 incoherence "
        f"signature, not manifold selection"
    )


def test_the_PUBLIC_ENTRY_returns_the_SAME_trace_under_BOTH_schedules():
    r"""⭐ The repaired side, and the counterpart the gate above needs.

    The driver's trace is splitting-dependent; **the trace a user receives is
    not**. That is the whole point of #344's Step 5, and until this row it was
    asserted only indirectly — :mod:`tests.sn.solve.test_every_entry_gauges_its_trace`
    gates that the gauge FIRES at every entry, that the spurious tangential
    current is gone and that the bulk is untouched, but nothing compared the two
    schedules' **returned** traces to each other.

    Read as a pair with
    :func:`test_the_UNGAUGED_SI_driver_returns_a_splitting_dependent_trace`:
    same box, same source, same two schedules, one tier apart. ``[M]`` the
    driver's traces are ``8.97e-02`` apart on this fixture and the entry's are
    ``2.5e-13``.

    ⚠ The activation leg is not optional. On an EVEN first axis with this
    symmetric source the kernel is never excited, so the two entries agree
    whether or not the gauge exists and the row would be inert while green
    (`vv` #13). ``gauge_correction > 1e-3`` is what proves the fixture is in the
    regime the claim is about.
    """
    import warnings

    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.operators.loss_kernel_gauge import GaugeFreedomWarning
    from orpheus.sn.solver import solve_sn_fixed_source
    from orpheus.transport.mesh.axis import AxisMesh

    cells = (3, 4)
    quad = Quadrature.level_symmetric(sn_order=4)
    entry_axes = tuple(
        AxisMesh(edges=np.linspace(0.0, extent, n + 1), bc_low=_R, bc_high=_R)
        for extent, n in zip((1.0, 2.0), cells))
    source = np.full((quad.weights.size, 2) + cells,
                     1.0 / float(quad.weights.sum()))

    returned = {}
    for schedule in ("gauss_seidel", "jacobi"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", GaugeFreedomWarning)
            solution = solve_sn_fixed_source(
                {0: absorber()}, entry_axes, quad, external_source=source,
                boundary_condition=None, inner_solver="source_iteration",
                inner_schedule=schedule, inner_tol=1e-13, max_inner=400_000,
            )
        assert solution.history is not None
        returned[schedule] = (
            np.asarray(solution.boundary_flux.values, dtype=float),
            solution.history.gauge_correction,
        )

    gs_trace, gs_correction = returned["gauss_seidel"]
    jacobi_trace, jacobi_correction = returned["jacobi"]
    assert gs_correction is not None and jacobi_correction is not None, (
        "an entry returned without gauging — the comparison below would then "
        "be about the driver, not about what a user receives"
    )
    assert gs_correction > 1e-3, (
        f"the fixture no longer excites the kernel "
        f"({gs_correction:.3e}) — the agreement below is then free, and this "
        f"row proves nothing. Check that n_x is still ODD"
    )
    agreement = float(np.max(np.abs(gs_trace - jacobi_trace))
                      / np.max(np.abs(jacobi_trace)))
    assert agreement < 1e-10, (
        f"the two schedules return DIFFERENT traces from the public entry "
        f"({agreement:.3e}) — the exit gauge is supposed to make the returned "
        f"trace a function of the problem rather than of the splitting"
    )


def test_the_UNGAUGED_gauss_seidel_driver_trace_error_is_first_order_in_h():
    r"""``err * n`` is CONSTANT — the driver's manifold selection CONVERGES.

    ⛔ **Tier**, as above: the SI driver, un-gauged. The public entries project
    this away entirely, so what is characterised here is the mechanism, not a
    user-visible error. The diagnostic's title —
    *"the gauss_seidel trace error is first order in h"* — read as a claim about
    a shipped error, and Step 5 of the campaign repealed that reading.

    Why it is worth gating at all: it is the leg that makes #344 a *determinism*
    disposition rather than a wrong-answer bug. An O(1) or non-converging
    deviation would mean boundary Gauss-Seidel converges to something that is
    not a solution; a clean O(h) says it converges to a different member of a
    manifold that is itself shrinking.

    ``[M]`` ``err * n = 0.29387244`` at **every** ``n`` in ``3, 5, 7, 9`` — 8
    significant figures over a 3x mesh range, on extents ``(1.0, 2.0)``. The
    ladder therefore starts at 3 (the diagnostic used 5, 7, 9 and cost
    ``[M]`` 19.7 s against 9.2 s here, for a *smaller* range).

    ⚠ **Two fixture facts are load-bearing and neither is optional.**

    * The ladder must stay inside ONE PARITY CLASS. At even ``n_x`` the
      deviation is ``~1e-12`` under this symmetric source, so a ``4, 8, 16``
      ladder reports "no effect" (`vv` #13's congruence-class trap). ⛔ That is
      a fact about the SOURCE's symmetry, not about the operator: ``dim ker A``
      is 12 at every parity.
    * Both axes must refine together. ``[M]`` on the ``(n, 4)`` family — only
      ``n_x`` refined — ``err * n`` reads ``0.269 / 0.307 / 0.319`` (spread
      ``1.7e-01``): that is not an ``h``-refinement at all, and using it would
      falsely refute the order.

    The ORDER is the claim and the spread pins it; the coefficient is a property
    of this box (``[M]`` ``0.311671`` on extents ``(1.0, 1.0)``) and is pinned
    loosely, as a mechanism tripwire.
    """
    ladder = (3, 5, 7)
    products = []
    for n_cells in ladder:
        sn_mesh, system, template, _n, n_bulk, source, exact = (
            uniform_source_fixture((n_cells, n_cells)))
        iterate = drive(system, sn_mesh, template, source, "gauss_seidel")
        deviation = float(np.max(np.abs(
            (iterate[n_bulk:] - exact[n_bulk:]) / exact[n_bulk:])))
        products.append(deviation * n_cells)

    assert min(products) > 1e-2, (
        f"the driver deviation collapsed ({products}) — an even n_x would do "
        f"that, and the order below would then be measured on round-off"
    )
    spread = (max(products) - min(products)) / float(np.mean(products))
    # [M] 2.124e-11 at tol=1e-13; a genuine order change reads >= 1e-1, so the
    # bar sits ~500x above the measurement and ~7 orders below a real failure.
    assert spread < 1e-8, (
        f"err * n is not constant, so the trace error is not O(h): "
        f"{products} at n = {ladder} (relative spread {spread:.3e})"
    )
    assert abs(float(np.mean(products)) - 0.29387244) < 1e-6, (
        f"the O(h) coefficient on extents (1.0, 2.0) moved to "
        f"{float(np.mean(products)):.8f} (was 0.29387244) — the schedule's "
        f"manifold selection changed; re-derive it, do not re-baseline"
    )

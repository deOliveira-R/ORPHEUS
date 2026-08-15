r"""The gauge fires at every public entry, and here is the defect it repairs.

Two claims, and they need different instruments.

**Coverage.** :func:`~orpheus.sn.solver._exit_gauge_trace` is the sibling of
``_exit_balance_defect`` with one sharpening: that one REPORTS and this one
MUTATES, so a forgotten call site does not lose a diagnostic — it silently
returns a non-physical answer. The structural guarantee a single ``Solution``
construction site would have given is unavailable (the two fixed-source arms
deliberately bypass ``_package_solution`` to keep their DG slope structure,
``solver.py:3765-3773``), so coverage is GATED here instead. The enumeration is
**derived from the module**, not hand-listed, so a new entry that forgets to
gauge cannot pass by being unknown to this file.

**The witness** (`plan-authoring` §6c — a gate that rejects nothing ships green
and unfalsifiable). ⛔ The obvious functional does NOT work: `[M]` the gauge
moves the trace by 6.1 % while the ``|Ω·n|⁰`` moment, the partial current and
the G-weighted total over Γ are **all bit-invariant** — and they must be, by
the campaign's own blindness theorem (a summed moment is mirror-EVEN, every
kernel mode is mirror-ODD). The functional that sees it is mirror-ODD: the
current **tangential** to a reflective face.

⭐ And that is not merely an observable, it is the defect. On an all-reflective
box with a uniform isotropic source the exact flux is flat, so every current is
zero everywhere. The solver already gets the NORMAL currents right (`[M]`
``~1e-15``). Un-gauged it reports a **spurious 7 % net current flowing sideways
along a mirror surface**, which cannot exist. The normal currents are therefore
a built-in negative control: they must not move.

⚠ **PARITY IS LOAD-BEARING IN EVERY FIXTURE HERE**, under the uniform
isotropic source these gates use: `[M]` ``(3,3) (3,4) (3,2) (5,4) (5,5)`` read
``4.1e-02 .. 7.8e-02`` and ``(4,3) (4,4) (4,5) (2,3) (2,2) (6,5)`` read
``2.1e-15 .. 1.0e-14``. A gate written on an even-``n_x`` mesh **with a
symmetric source** is inert. `vv` #13's congruence-class trap: a sub-agent exit
map probed only 4×4 and 6×6 and concluded the gauge could not bite at the
eigenvalue entry at all.

⛔ **But "odd ``n_x``" is NOT the rule, and stating it as one was this
campaign's own quantifier error** (`plan-authoring` §2 — the unwritten
denominator was *"…with a uniform isotropic source"*). `[M]`
``dim ker A = 12`` on **every** mesh here, odd or even: ``(2,2) (3,4) (4,4)
(5,6) (6,8)`` alike. The operator is singular whenever
:func:`~orpheus.sn.operators.loss_kernel_gauge.gauge_freedom` says so, full
stop. What parity governs is whether the ITERATION **excites** the kernel, and
that is a joint property of the mesh AND the source's symmetry:

============================  ===================  ===================
source                        even ``n_x``         odd ``n_x``
============================  ===================  ===================
uniform isotropic             `[M]` ``~1e-14``     `[M]` ``4.1e-02 ..
                              (not excited)        7.8e-02`` (excited)
anisotropic ``(1+mu_x)/W``    `[M]` ``1.03e-02 ..  excited
                              3.60e-02`` (EXCITED)
============================  ===================  ===================

⛔ **Never infer kernel-freedom from a mesh property** — ask the predicate, or
assert ``dim ker == 0``. An even-``n_x`` fixture is unexcited only for a
symmetric source; `[M]` ``(4,4)`` reads ``6.7e-14`` under a uniform source and
``1.756363e-02`` under an anisotropic one.
"""

from __future__ import annotations

import inspect
import warnings

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.face_layout import face_normal
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solver as solver_module
from orpheus.sn.operators.loss_kernel_gauge import (
    GAUGE_ESCALATION_FLAG,
    GaugeFreedomWarning,
)
from orpheus.sn.solver import (
    _as_sn_mesh,
    solve_sn,
    solve_sn_fixed_source,
)
from orpheus.transport.mesh.axis import AxisMesh

_R = BC("reflective")
_QUAD = Quadrature.level_symmetric(sn_order=4)

#: ⚠ ODD first axis — see the module docstring. Paired with the uniform
#: isotropic source these gates use, an EVEN one leaves every gate here inert
#: and green. (The operator is singular either way; parity only decides whether
#: a SYMMETRIC source excites the kernel.)
_EXCITED_CELLS = (3, 4)
_UNEXCITED_CELLS = (4, 4)


def _absorber(ng: int = 2):
    sig_t = np.linspace(0.8, 1.6, ng)
    return make_mixture(
        sig_t=sig_t, sig_c=sig_t.copy(), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=np.zeros((ng, ng)),
    )


def _reflective_mesh2d(cells):
    return Mesh2D(
        edges_x=np.linspace(0.0, 1.0, cells[0] + 1),
        edges_y=np.linspace(0.0, 2.0, cells[1] + 1),
        mat_map=np.zeros(cells, dtype=int),
        bc_xmin=_R, bc_xmax=_R, bc_ymin=_R, bc_ymax=_R,
    )


def _reflective_axes(cells):
    return tuple(
        AxisMesh(edges=np.linspace(0.0, extent, n + 1), bc_low=_R, bc_high=_R)
        for extent, n in zip((1.0, 2.0), cells)
    )


def _uniform_source(cells, ng: int = 2):
    """`Q = 1/W` per ordinate ⟹ the exact flux is FLAT, so every current is 0."""
    return np.full((_QUAD.weights.size, ng) + tuple(cells),
                   1.0 / float(_QUAD.weights.sum()))


def _currents(solution, sn_mesh) -> dict[tuple[str, str], float]:
    r"""``{(face, "normal"|"tangential"): J}`` from the RETURNED trace.

    :math:`J_b = \sum_n w_n \mu_b \psi_n` summed over the face's cells — the
    signed cosine, so it is mirror-ODD in axis ``b`` and can see the kernel.
    """
    trace = np.asarray(solution.boundary_flux.values, dtype=float)
    weights = np.asarray(sn_mesh.quad.weights, dtype=float)
    cosines = {0: np.asarray(sn_mesh.quad.mu_x, dtype=float),
               1: np.asarray(sn_mesh.quad.mu_y, dtype=float)}
    out: dict[tuple[str, str], float] = {}
    for face, slot in sn_mesh.angular_trace.layout.faces.items():
        axis, _ = face_normal(face)
        per_ordinate = slot.slice_view(trace).reshape(slot.shape[0], -1).sum(1)
        for component in (0, 1):
            role = "normal" if component == axis else "tangential"
            out[(face, role)] = float(
                (weights * cosines[component]) @ per_ordinate)
    return out


# ─────────────────────────────────────────────────────────────────────
# 1. COVERAGE — no public entry may skip the gauge, and none may hide
# ─────────────────────────────────────────────────────────────────────
#: Every public entry, with why it is or is not exercisable on a SINGULAR
#: configuration. An exemption must state a STRUCTURAL reason, not an
#: inconvenience — see :func:`test_no_public_entry_is_unaccounted_for`.
_ENTRY_LEDGER = {
    "solve_sn": "exercised — the eigenvalue entry, all-reflective by default",
    "solve_sn_fixed_source": "exercised — under a declared reflective BC",
    "solve_sn_adjoint": (
        "wired but STRUCTURALLY INERT: the adjoint routes through (L+C)^H, "
        "whose transpose solve is 1-D-scan-only (#280 Phase 2.5b), and a 1-D "
        "problem has at most ONE reflective axis pair, so gauge_freedom is "
        "never present on any configuration this entry can run"
    ),
    "solve_sn_adjoint_fixed_source": "wired but structurally inert — as above",
}


@pytest.mark.foundation
def test_no_public_entry_is_unaccounted_for():
    """The entry list is DERIVED from the module, not hand-maintained here.

    A new ``solve_sn*`` entry that forgets to gauge would otherwise pass this
    file by being unknown to it — the failure mode a hand-written list cannot
    catch, and the one that matters because the gauge MUTATES.
    """
    discovered = {
        name for name, obj in vars(solver_module).items()
        if name.startswith("solve_sn") and inspect.isfunction(obj)
        and obj.__module__ == solver_module.__name__
    }
    assert discovered, "the discovery predicate matched nothing — it is broken"
    assert discovered == set(_ENTRY_LEDGER), (
        f"public entries not accounted for in this file: "
        f"{sorted(discovered - set(_ENTRY_LEDGER))}; "
        f"ledger rows with no entry: {sorted(set(_ENTRY_LEDGER) - discovered)}"
    )


@pytest.mark.foundation
def test_every_exercisable_entry_reports_a_gauge_correction():
    """Both forward entries populate ``gauge_correction`` on a singular config.

    ``None`` here would mean the entry never called the gauge — which is
    exactly the silent-skip this gate exists to refuse. A *measured* value is
    the proof it fired, whatever its magnitude.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GaugeFreedomWarning)
        eigen = solve_sn(
            {0: get_mixture("A", "2g")},
            _reflective_mesh2d(_EXCITED_CELLS), _QUAD,
            inner_solver="source_iteration", inner_schedule="gauss_seidel",
            inner_tol=1e-13,
        )
    assert eigen.history is not None
    assert eigen.history.gauge_correction is not None, "solve_sn did not gauge"
    assert eigen.history.gauge_correction > 1e-3

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GaugeFreedomWarning)
        fixed = solve_sn_fixed_source(
            {0: _absorber()}, _reflective_axes(_EXCITED_CELLS), _QUAD,
            external_source=_uniform_source(_EXCITED_CELLS),
            boundary_condition=None, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
    assert fixed.history is not None
    assert fixed.history.gauge_correction is not None, (
        "solve_sn_fixed_source did not gauge"
    )
    assert fixed.history.gauge_correction > 1e-3


@pytest.mark.foundation
@pytest.mark.parametrize("schedule", ["gauss_seidel", "krylov"])
def test_BOTH_fixed_source_arms_gauge_not_just_one(schedule):
    """The SI and Krylov arms both project — the Krylov comment notwithstanding.

    ``solver.py`` used to describe the Krylov arm's returned boundary as *"the
    matvec's B1'' face residual"*, which would have made a reader exempt it.
    `[M]` it is a flux trace, and the two arms return **different members of
    the solution set** when the kernel is non-trivial — SI carrying ~2.9 % of
    it and Krylov ~1e-12. Exempting Krylov would leave the arms' post-gauge
    status decided by an accident of the driver.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solution = solve_sn_fixed_source(
            {0: _absorber()}, _reflective_axes(_EXCITED_CELLS), _QUAD,
            external_source=_uniform_source(_EXCITED_CELLS),
            boundary_condition=None,
            inner_solver="krylov" if schedule == "krylov"
            else "source_iteration",
            inner_schedule="gauss_seidel",
            inner_tol=1e-12, max_inner=40_000,
        )
    assert solution.history is not None
    assert solution.history.gauge_correction is not None, (
        f"the {schedule} arm returned without gauging"
    )


# ─────────────────────────────────────────────────────────────────────
# 2. ⭐ THE WITNESS — a spurious tangential current, and its control
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
def test_the_spurious_TANGENTIAL_current_along_a_mirror_is_gone():
    r"""⭐ The acceptance gate: ~7 % → round-off, with the normal currents pinned.

    The exact answer here is flat (uniform isotropic source, all-reflective),
    so **every** current component is zero. `[M]` un-gauged, the returned trace
    carries :math:`J_{\rm tangential} = 7.381060\times10^{-2}` on both y-faces
    — a net neutron current running sideways along a mirror, which is not a
    thing. Gauged, it is round-off.

    The NORMAL currents are the negative control: the solver already gets them
    right (`[M]` ``~1e-15``) and the gauge must leave them there. Without that
    leg this gate would be satisfied by any transformation that shrinks the
    trace.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GaugeFreedomWarning)
        solution = solve_sn_fixed_source(
            {0: _absorber()}, _reflective_axes(_EXCITED_CELLS), _QUAD,
            external_source=_uniform_source(_EXCITED_CELLS),
            boundary_condition=None, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
    sn_mesh = _as_sn_mesh(_reflective_axes(_EXCITED_CELLS), _QUAD,
                          {0: _absorber()})
    currents = _currents(solution, sn_mesh)

    tangential = {k: v for k, v in currents.items() if k[1] == "tangential"}
    normal = {k: v for k, v in currents.items() if k[1] == "normal"}
    assert len(tangential) == len(normal) == 4

    worst_tangential = max(abs(v) for v in tangential.values())
    assert worst_tangential < 1e-11, (
        f"a tangential current survives on a mirror face: {tangential}"
    )
    worst_normal = max(abs(v) for v in normal.values())
    assert worst_normal < 1e-11, (
        f"the NEGATIVE CONTROL moved — the normal currents were already "
        f"correct and must stay so: {normal}"
    )
    # …and the fixture still exhibits the defect it is here to catch.
    assert solution.history is not None
    assert solution.history.gauge_correction is not None
    assert solution.history.gauge_correction > 1e-3, (
        f"fixture no longer excites the kernel "
        f"({solution.history.gauge_correction:.3e}) — check n_x is still ODD; "
        f"an even first axis makes this gate inert while leaving it green"
    )


@pytest.mark.foundation
def test_the_bulk_and_keff_are_untouched():
    """The kernel is pure-trace, so nothing a user reads from the bulk moves.

    Pinned against a mesh whose parity makes the gauge inert: same ``keff`` to
    full precision. If gauging ever reached the bulk, these two would part.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GaugeFreedomWarning)
        excited = solve_sn(
            {0: get_mixture("A", "2g")},
            _reflective_mesh2d(_EXCITED_CELLS), _QUAD,
            inner_solver="source_iteration", inner_schedule="gauss_seidel",
            inner_tol=1e-13,
        )
    unexcited = solve_sn(
        {0: get_mixture("A", "2g")}, _reflective_mesh2d(_UNEXCITED_CELLS),
        _QUAD, inner_solver="source_iteration", inner_schedule="gauss_seidel",
        inner_tol=1e-13,
    )
    assert excited.keff == pytest.approx(unexcited.keff, rel=1e-12)
    # k_inf = nu*SigF / SigA is flux-shape independent on a homogeneous
    # all-reflective box — an INDEPENDENT anchor, not a self-comparison.
    assert excited.keff == pytest.approx(1.875, rel=1e-9)


@pytest.mark.foundation
@pytest.mark.verifies("sn-kernel-mirror-blindness")
def test_every_MIRROR_EVEN_functional_is_blind_to_the_gauge():
    r"""⭐ The blindness theorem, asserted — and why the obvious gate is inert.

    Every kernel mode is mirror-ODD, so any mirror-EVEN functional annihilates
    it. A face-SUMMED moment is mirror-even. ⟹ the ``|Ω·n|⁰`` moment, the
    partial current and the G-weighted total over Γ cannot move, no matter how
    large the correction — `[M]` they read ``0.0`` while the trace itself moves
    6.1 %.

    This is the theorem that makes the gauge CANONICAL (``ψ_exact`` is
    mirror-even, hence G-orthogonal to ``ker A``, hence IS the minimum-norm
    member). It is also the reason
    :func:`test_the_spurious_TANGENTIAL_current_along_a_mirror_is_gone` had to
    reach for a mirror-ODD functional: the campaign's first witness design
    asserted the ``|Ω·n|⁰`` moment and would have shipped green and
    structurally unable to fail (`plan-authoring` §6c).
    """
    axes = _reflective_axes(_EXCITED_CELLS)
    sn_mesh = _as_sn_mesh(axes, _QUAD, {0: _absorber()})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GaugeFreedomWarning)
        solution = solve_sn_fixed_source(
            {0: _absorber()}, axes, _QUAD,
            external_source=_uniform_source(_EXCITED_CELLS),
            boundary_condition=None, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
    gauged = np.asarray(solution.boundary_flux.values, dtype=float)
    gauge = sn_mesh.loss_kernel_gauge
    weights = np.asarray(_QUAD.weights, dtype=float)
    metric = np.asarray(sn_mesh.angular_trace.inner_product_weights, dtype=float)

    # A DIFFERENT member of the same solution manifold: the returned trace plus
    # an arbitrary kernel vector. Both solve the discrete system exactly, so any
    # functional that distinguishes them is one the equation does not determine.
    kernel_member = gauge.apply(np.random.default_rng(7).standard_normal(
        metric.size))
    assert np.linalg.norm(kernel_member) > 1e-6, "degenerate: nothing in ker A"
    other_member = gauged + kernel_member

    for face, slot in sn_mesh.angular_trace.layout.faces.items():
        face_index = list(sn_mesh.angular_trace.layout.faces).index(face)
        normal = np.abs(
            np.asarray(sn_mesh.angular_trace.omega_dot_n)[face_index])
        before = slot.slice_view(other_member).reshape(slot.shape[0], -1).sum(1)
        after = slot.slice_view(gauged).reshape(slot.shape[0], -1).sum(1)
        for name, weighting in (("|Omega.n|^0 moment", weights),
                                ("partial current", weights * normal)):
            lhs, rhs = float(weighting @ before), float(weighting @ after)
            assert lhs == pytest.approx(rhs, abs=1e-11), (
                f"{face} {name} moved under a change of manifold member: "
                f"{lhs:.6e} -> {rhs:.6e} — the blindness theorem says it "
                f"cannot"
            )
    total_before = float(np.sum(metric * other_member))
    total_after = float(np.sum(metric * gauged))
    assert total_before == pytest.approx(total_after, abs=1e-11)
    # …and the two members really are different, or the above is vacuous.
    assert (np.linalg.norm(other_member - gauged)
            / np.linalg.norm(gauged)) > 1e-3


@pytest.mark.foundation
def test_an_EVEN_mesh_is_excited_too_once_the_source_stops_being_symmetric():
    r"""⛔ The gate that refutes "excited iff ``n_x`` is odd".

    That rule was this campaign's own quantifier error (`plan-authoring` §2):
    it was measured over 11 meshes **all carrying a uniform isotropic source**,
    and written up as a property of the MESH. It is not.

    `[M]` ``dim ker A = 12`` at ``(2,2) (3,4) (4,4) (5,6) (6,8)`` alike — the
    operator is singular at every parity. Parity only decides whether a
    *symmetric* source excites the kernel. Break the symmetry and an even mesh
    is excited too: `[M]` ``(4,4)`` reads ``6.735136e-14`` under a uniform
    source and **``1.756363e-02``** under ``(1 + mu_x)/W``.

    This row exists so the corrected understanding is pinned rather than
    remembered — and so a reader who meets ``_UNEXCITED_CELLS`` elsewhere in
    this file cannot conclude that an even mesh is safe.
    """
    cells = _UNEXCITED_CELLS
    assert cells[0] % 2 == 0, "this row needs an EVEN first axis to mean anything"
    axes = _reflective_axes(cells)
    sn_mesh = _as_sn_mesh(axes, _QUAD, {0: _absorber()})

    # The operator does not care about the source, and never did.
    assert sn_mesh.loss_kernel_gauge.dimension > 0, (
        "the configuration must be singular for either row below to mean "
        "anything — this is the claim the parity rule obscured"
    )

    mu_x = np.asarray(_QUAD.mu_x, dtype=float)
    total_weight = float(_QUAD.weights.sum())
    anisotropic = (
        ((1.0 + mu_x) / total_weight)[:, None, None, None]
        * np.ones((mu_x.size, 2) + tuple(cells))
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GaugeFreedomWarning)
        skewed = solve_sn_fixed_source(
            {0: _absorber()}, axes, _QUAD, external_source=anisotropic,
            boundary_condition=None, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
        symmetric = solve_sn_fixed_source(
            {0: _absorber()}, axes, _QUAD,
            external_source=_uniform_source(cells),
            boundary_condition=None, inner_solver="source_iteration",
            inner_schedule="gauss_seidel", inner_tol=1e-13, max_inner=400_000,
        )
    assert skewed.history is not None and symmetric.history is not None
    assert skewed.history.gauge_correction is not None
    assert symmetric.history.gauge_correction is not None

    assert skewed.history.gauge_correction > 1e-3, (
        f"an anisotropic source on an EVEN mesh must excite the kernel; got "
        f"{skewed.history.gauge_correction:.3e}"
    )
    assert symmetric.history.gauge_correction < 1e-10, (
        f"…and the symmetric source on the SAME mesh must not — that contrast "
        f"is the whole content of this row; got "
        f"{symmetric.history.gauge_correction:.3e}"
    )


# ─────────────────────────────────────────────────────────────────────
# 3. The warning — audible when it acted, silent when it did not
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.foundation
def test_the_repair_is_AUDIBLE_and_names_the_root_fix():
    """R2b: say the fixing was applied, and that a closure change removes it.

    ⭐ The root-fix clause is ASKED of the scheme registry, not tabulated, so
    it names the closures that actually damp the face mode in THIS build.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        solve_sn(
            {0: get_mixture("A", "2g")},
            _reflective_mesh2d(_EXCITED_CELLS), _QUAD,
            inner_solver="source_iteration", inner_schedule="gauss_seidel",
            inner_tol=1e-13,
        )
    gauge_warnings = [
        w for w in caught if issubclass(w.category, GaugeFreedomWarning)]
    assert len(gauge_warnings) == 1, (
        f"expected exactly one GaugeFreedomWarning, got {len(gauge_warnings)}")
    message = str(gauge_warnings[0].message)
    assert "GAUGE-FIXED" in message
    assert "6.08%" in message, f"the magnitude is not quoted: {message}"
    # the ROOT fix, and the escalation recipe
    assert "linear_discontinuous" in message, (
        f"the root fix must name the closures that damp the mode: {message}")
    assert GAUGE_ESCALATION_FLAG in message
    assert "error::orpheus.sn.operators.loss_kernel_gauge." in message, (
        "the escalation category must be DOTTED — `-W` resolves an undotted "
        "one against builtins, so the short spelling gates nothing"
    )


@pytest.mark.foundation
def test_it_is_SILENT_when_the_answer_was_already_canonical():
    """No warning on an unexcited mesh — and the record still carries the number.

    The warning reports an ACTION, not a configuration property. All-reflective
    is the DEFAULT for the eigenvalue entry, so warning on the configuration
    would fire on the standard ``k_inf`` lattice every time.

    ⚠ The distinction is visible only because ``gauge_correction`` is still
    *measured* here (`[M]` ``~5e-15``, not ``None``): the freedom is real, the
    solve simply landed on the canonical member.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        solution = solve_sn(
            {0: get_mixture("A", "2g")},
            _reflective_mesh2d(_UNEXCITED_CELLS), _QUAD,
            inner_solver="source_iteration", inner_schedule="gauss_seidel",
            inner_tol=1e-13,
        )
    assert not [
        w for w in caught if issubclass(w.category, GaugeFreedomWarning)]
    assert solution.history is not None
    assert solution.history.gauge_correction is not None, (
        "silence must come from the gauge having done nothing, NOT from it "
        "never having run — those are different states and only one is fine"
    )
    assert solution.history.gauge_correction < 1e-10


@pytest.mark.foundation
@pytest.mark.parametrize(
    "entry", ["solve_sn", "solve_sn_fixed_source"],
)
def test_the_warning_blames_the_CALLER_not_orpheus(entry):
    r"""``stacklevel`` attribution, asserted by exact file AND line.

    ⚠ A gate of the form *"the warning is not attributed inside orpheus/"* is
    Mode-12 blind to ``stacklevel`` being too LARGE — it would pass while
    blaming the user's caller's caller. Only an exact ``(filename, lineno)``
    pin catches both directions.

    The fixed-source row is the load-bearing one: its projection happens two
    frames down in a private arm, so the warning had to be **hoisted** to the
    public entry. Emitted from the arm, ``stacklevel=3`` lands on
    ``orpheus/sn/solver.py`` — verbatim the defect #340 N4.7 measured at 2 of
    8 sites.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if entry == "solve_sn":
            expected_line = inspect.currentframe().f_lineno + 1  # type: ignore[union-attr]
            solve_sn(
                {0: get_mixture("A", "2g")},
                _reflective_mesh2d(_EXCITED_CELLS), _QUAD,
                inner_solver="source_iteration",
                inner_schedule="gauss_seidel", inner_tol=1e-13,
            )
        else:
            expected_line = inspect.currentframe().f_lineno + 1  # type: ignore[union-attr]
            solve_sn_fixed_source(
                {0: _absorber()}, _reflective_axes(_EXCITED_CELLS), _QUAD,
                external_source=_uniform_source(_EXCITED_CELLS),
                boundary_condition=None, inner_solver="source_iteration",
                inner_schedule="gauss_seidel", inner_tol=1e-13,
                max_inner=400_000,
            )
    gauge_warnings = [
        w for w in caught if issubclass(w.category, GaugeFreedomWarning)]
    assert len(gauge_warnings) == 1
    emitted = gauge_warnings[0]
    assert emitted.filename == __file__, (
        f"{entry} blamed {emitted.filename}, not the caller — stacklevel is "
        f"wrong (too small blames orpheus, too large blames the caller's "
        f"caller)"
    )
    assert emitted.lineno == expected_line, (
        f"{entry} blamed line {emitted.lineno}, expected {expected_line}"
    )

r"""The §4.6 non-vacuum MMS, driven through a **DECLARED** boundary law.

Phase **P4** of `.claude/plans/affine_boundary_source_channel.md`; design of
record `scratch/p4_mms_design.md`. Sibling to ``test_mms_prescribed_inflow.py``,
which stays — the two together are the two-route claim:

* **the sibling** hands ``solve_sn_fixed_source`` a composite
  ``q = q_bulk ⊕ q_∂`` built by ``case.prescribed_inflow(sn)``. That is the
  **channel tier**, and it deliberately bypasses the law tier
  (``docs/theory/verification/sn.rst`` says so).
* **this module** declares ``PrescribedInflow(source=<manufactured spec>)`` on
  the GEOMETRY and passes a **bulk-only** source. ``from_mesh_laws`` assembles
  ``q_∂`` from the declaration. That is the **law tier** — every step a user
  actually takes.

⭐ **Nothing above the delivery line is new.** The ansatz, the manufactured
``Q^ext``, ``phi_exact``, the SymPy provenance and the ``O(h²)`` claim are all
inherited from :func:`build_slab_2g_nonvacuum_mms_case`. P4 changes exactly one
thing — *how ``q_∂`` arrives* — which is why it is a **re-routing**, not a new
manufactured solution.

The ansatz (inherited), per group :math:`g`:

.. math::

    \psi_{n,g}(x) = \frac{A_g(x) + \mu_n B_g(x)}{W},
    \qquad A_g = c_g(a_0 + a_1\sin kx), \quad B_g = c_g b_0 \cos kx

with :math:`a_0 > 0`, which is what makes :math:`\gamma_-\psi \neq 0` at the
faces. The :math:`\mu_n B_g` term is what makes the inflow **angle-varying** —
a constant-in-angle inflow is expressible by the shipped
``ConstantInflowSource``, so it would exercise none of the machinery this phase
is about.

⚠ **The user path stops one step short of ideal, and this module says so rather
than implying otherwise.** The declaration itself is fully public since the
declaration channel landed (`985497b5`): a law object is a legal ``bcs=`` entry,
and ``solve_sn_fixed_source`` resolves it. What is still a stopgap is
:class:`_ManufacturedFaceInflow` — every constructor argument it takes is
information the shipped ``InflowSourceSpec`` Protocol cannot supply. That is
**P6's specification**, enumerated in `scratch/p4_mms_design.md` §10.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_slab_2g_nonvacuum_mms_case,
)
from orpheus.geometry.boundary import PrescribedInflow
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.source_sinks import AngularBoundarySourceSink
from tests.sn._test_helpers import volume_weighted_l2

pytestmark = pytest.mark.l1


def _require(condition: bool, message: str) -> None:
    """``-O``-survivable assertion (a bare ``assert`` is stripped under ``-O``)."""
    if not condition:
        raise AssertionError(message)


# ═══════════════════════════════════════════════════════════════════════
# The P6 stopgap — a spec that knows what the Protocol cannot tell it
# ═══════════════════════════════════════════════════════════════════════


class _ManufacturedFaceInflow:
    r""":math:`\gamma_-\psi = (A_g(x_f) + \mu_n B_g(x_f))/W` on ONE face.

    ⭐ **The P6 payoff: this is now an ORDINARY user-written source.** It takes
    the manufactured case and the face's coordinate — both known when you
    DECLARE the boundary — and nothing else. Everything it used to smuggle
    through its constructor now arrives with the space:

    ===========================  =============================================
    was a constructor argument   now read from ``space``
    ===========================  =============================================
    ``mu_inflow``                ``space.directions`` — the per-row Ω, in this
                                 space's own row order, so a permuted trace is
                                 unspellable rather than silently wrong
    ``n_ordinates``              nothing needs it: there is ONE call at ONE
                                 shape, so there is no probe to recognise
    ===========================  =============================================

    The two that remain are genuinely the author's: :math:`A_g` / :math:`B_g`
    come from the ansatz at this face's coordinate, and the :math:`1/W` is the
    ansatz's own normalisation (``φ = Σ w ψ``), which is why the Protocol
    documents the units rather than trying to supply the factor. ⚠ It is the
    FULL weight sum — ``case.quadrature.weights.sum()`` — not the
    inflow-restricted one; `[M]` those are 2.0 and 1.0 on GL-8 and the wrong
    choice is a silent ×2 in the same direction as a double delivery.

    ⛔ **What this class no longer has to do**, kept as the record of what the
    old ``evaluate(shape)`` signature cost: recognise a rank-1 realize-time
    probe by its rank, return a deliberately non-zero value into it so the
    ERR-047 certification would not be skipped, and be rebuilt per face from a
    throwaway probe ``SNMesh`` just to learn the row order.
    """

    def __init__(self, *, case, x_face: float) -> None:
        self.case = case
        self.x_face = float(x_face)
        self.W = float(np.asarray(case.quadrature.weights).sum())
        #: Every space this spec was asked for, by name — G8 reads it.
        self.spaces_seen: list[str] = []

    def _coefficients(self) -> "tuple[np.ndarray, np.ndarray]":
        x = np.array([self.x_face])
        ng = self.case.n_groups
        A = np.array([float(np.atleast_1d(self.case.A(x, g))[0])
                      for g in range(ng)])
        B = np.array([float(np.atleast_1d(self.case.B(x, g))[0])
                      for g in range(ng)])
        return A, B

    def evaluate(self, space) -> np.ndarray:
        self.spaces_seen.append(space.name)
        # γ₋ψ is the angular FLUX on the trace, so it carries the ansatz's own
        # 1/W — a DIFFERENT reason from the bulk source's /W (design §1.3/§1.4).
        mu = np.asarray(space.directions, dtype=float)
        mu = mu if mu.ndim == 1 else mu[:, 0]
        A, B = self._coefficients()
        return (A[None, :] + mu[:, None] * B[None, :]) / self.W


# ═══════════════════════════════════════════════════════════════════════
# Fixture — declare on the GEOMETRY, solve through the PUBLIC entry
# ═══════════════════════════════════════════════════════════════════════


def _face_spec(case, face: str, x_face: float) -> _ManufacturedFaceInflow:
    """No mesh argument, and no ``face`` read — the space carries both.

    ``face`` survives only as documentation of which coordinate ``x_face`` is;
    the spec never branches on it.
    """
    return _ManufacturedFaceInflow(case=case, x_face=x_face)


def _declared_solve(case, n_cells: int, inner: str = "source_iteration"):
    r"""Build → declare → solve, entirely through public surfaces.

    ⭐ **P6 deleted a step here.** This used to build a throwaway probe
    ``SNMesh`` first, purely so the spec could be handed the inflow row order —
    "a spec that received the trace would not need the probe", as the note then
    said. It receives the trace now, so the probe is gone and the fixture is
    build → declare → solve with nothing in front of it.

    Returns ``(sn, solution, specs)``.
    """
    mesh0 = case.build_mesh(n_cells)
    specs = {
        face: _face_spec(case, face, x_face)
        for face, x_face in (("xmin", 0.0), ("xmax", case.slab_length))
    }

    # ⭐ THE DECLARATION. A law object is a legal geometry-level declaration
    # since `985497b5`; `solve_sn_fixed_source` rebuilds the SNMesh from this
    # mesh and `resolve_boundary_conditions` reads the law straight off it.
    mesh = replace(
        mesh0,
        bc_left=PrescribedInflow(source=specs["xmin"]),
        bc_right=PrescribedInflow(source=specs["xmax"]),
    )

    # A BULK-ONLY source: `_build_fixed_source_rhs`'s array arm calls
    # `from_mesh_laws`, so q_∂ comes from the DECLARATION, not from the caller.
    solution = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, case.external_source(mesh),
        inner_solver=inner, max_inner=1000, inner_tol=1e-13,
    )
    return SNMesh(mesh, case.quadrature, case.materials), solution, specs


def _gamma_minus(solution, sn, face: str) -> np.ndarray:
    rows = np.asarray(sn.angular_trace.inflow_indices_for_face(face))
    return np.asarray(solution.angular_flux.boundary.face_view(face))[rows]


def _expected_trace(case, sn, face: str, x_face: float) -> np.ndarray:
    """The oracle: the CASE's own ansatz, never the spec (design §7.2)."""
    rows = np.asarray(sn.angular_trace.inflow_indices_for_face(face))
    mu = case.quadrature.mu_x[rows]
    W = float(case.quadrature.weights.sum())
    return np.stack(
        [(case.A(x_face, g) + mu * case.B(x_face, g)) / W
         for g in range(case.n_groups)],
        axis=1,
    )


def _faces(case):
    return (("xmin", 0.0), ("xmax", case.slab_length))


# ═══════════════════════════════════════════════════════════════════════
# G1 — ⭐⭐ THE KEYSTONE
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.catches("ERR-075")
@pytest.mark.verifies("bc-single-delivery", "sn-mms-nonvacuum-psi")
@pytest.mark.parametrize("inner", ["source_iteration", "krylov"])
def test_the_declared_manufactured_trace_holds_on_the_answer(inner: str) -> None:
    r"""⭐⭐ :math:`\gamma_-\psi|_f = (A_g + \mu_n B_g)/W`, per face, group, ordinate.

    The boundary condition is a **definition**, so the converged inflow trace IS
    the declared source — no reference solver, no discretization assumption. And
    because the declared trace varies in BOTH the ordinate and the group axis,
    this row pins far more than the delivery count: the row ORDER, the group
    axis, the face assignment, and the :math:`1/W`.

    ⚠ **A ``2×`` reading is NOT self-attributing.** On Gauss–Legendre
    :math:`W = 2`, so "the source was delivered twice" and "the :math:`1/W` was
    dropped" are *the same number*. Discriminate with the ordinate structure: a
    doubled delivery scales :math:`A` and :math:`B` alike, whereas a dropped
    :math:`1/W` on one term only changes the :math:`\mu`-slope relative to the
    offset.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    sn, solution, _ = _declared_solve(case, 40, inner)

    for face, x_face in _faces(case):
        got = _gamma_minus(solution, sn, face)
        want = _expected_trace(case, sn, face, x_face)
        message = (
            f"[{inner}] γ₋ψ({face}) does not equal the DECLARED manufactured "
            f"trace. max|Δ| = {float(np.max(np.abs(got - want))):.3e}. A factor "
            f"~2 is ambiguous: on Gauss-Legendre W = 2, so a doubled delivery "
            f"and a dropped 1/W are the same number — check whether the "
            f"mu-slope moved with the offset."
        )
        if inner == "source_iteration":
            # `[M]` bit-exact: SI writes q into the inflow slot and sweeps from
            # it, so the converged trace is a COPY of the declared values.
            np.testing.assert_array_equal(got, want, err_msg=message)
        else:
            # `[M]` Krylov SOLVES the trace rows, so they carry the iteration
            # residual. A ULP budget for THIS fixture at inner_tol = 1e-13 —
            # not a floating-point law; loosening the tolerance without
            # re-deriving it is a legitimate red.
            _require(
                np.max(np.abs(got - want)) < 1e-9,
                message,
            )
            np.testing.assert_array_almost_equal_nulp(got, want, nulp=256)


def test_the_declared_trace_actually_varies_in_both_axes() -> None:
    r"""Activation guard for G1 — without it the keystone could be vacuous.

    If :math:`B_g(x_f)` happened to vanish at both faces the trace would be
    constant in angle, and G1 would pass against an implementation that
    broadcast one value across :math:`\Gamma_-`. If the two groups' amplitudes
    coincided, the group axis would be equally undiscriminating.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    sn = SNMesh(case.build_mesh(20), case.quadrature, case.materials)
    for face, x_face in _faces(case):
        want = _expected_trace(case, sn, face, x_face)
        spread = float(np.max(want[:, 0]) - np.min(want[:, 0]))
        _require(
            spread > 1e-3,
            f"the declared trace on {face} is nearly constant in angle "
            f"(spread {spread:.3e}) — G1 cannot see a row-order defect",
        )
        _require(
            float(np.max(np.abs(want[:, 0] - want[:, 1]))) > 1e-3,
            f"the two groups' traces on {face} coincide — G1 cannot see a "
            f"collapsed group axis",
        )


# ═══════════════════════════════════════════════════════════════════════
# G2 / G3 — convergence, and the value it converges TO
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.verifies(
    "transport-cartesian", "dd-slab",
    "sn-mms-nonvacuum-qext", "sn-mms-nonvacuum-qext-mg",
)
def test_the_declared_route_converges_second_order() -> None:
    r"""``O(h²)`` per group, through the DECLARED channel.

    The sibling module already pins this for the hand-supplied channel; the
    claim here is that re-routing the delivery does not disturb the order. A
    mesh-independent boundary error would show up as an order FLOOR — the rate
    collapsing toward 0 as the bulk truncation shrinks past it.

    ⚠ This gate's sensitivity floor is ``~5e-4`` relative on ``q`` (design §7.1),
    four orders coarser than G1's ``~3e-12``, and in between the two a ``q``
    perturbation partially CANCELS the ``O(h²)`` truncation. G1 is the keystone
    for exactly that reason; this row corroborates consistency, it does not
    police the value.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    n_cells = [20, 40, 80]
    errors = {g: [] for g in range(case.n_groups)}
    for nc in n_cells:
        sn, solution, _ = _declared_solve(case, nc)
        mesh = sn.mesh
        for g in range(case.n_groups):
            phi = np.asarray(solution.scalar_flux.values)[g, :]
            exact = case.phi_exact(mesh.centers, g)
            errors[g].append(volume_weighted_l2(phi, exact, mesh.volumes))

    for g in range(case.n_groups):
        e = errors[g]
        orders = [np.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)]
        _require(
            all(o > 1.9 for o in orders),
            f"group {g} lost second order through the DECLARED channel: "
            f"orders={[f'{o:.3f}' for o in orders]}, errors={e}. An order "
            f"FLOOR (rate collapsing as h shrinks) means a mesh-INDEPENDENT "
            f"boundary error — i.e. the declared q is wrong by a constant.",
        )


def test_it_converges_to_the_manufactured_value_not_merely_at_a_rate() -> None:
    r"""``O(h²)`` to the WRONG limit is still ``O(h²)`` (``vv`` anti-pattern #5).

    The failure this guards is specific and plausible: if the declaration were
    silently dropped, the solve would be a clean vacuum problem — and it would
    converge at a perfect second order, to the boundary-ZERO solution. The rate
    gate above cannot tell those apart; only comparing the converged VALUE can.

    The non-vacuum activation legs make that concrete by asserting the exact
    solution really is non-zero at the faces.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    sn, solution, _ = _declared_solve(case, 80)
    for g in range(case.n_groups):
        phi = np.asarray(solution.scalar_flux.values)[g, :]
        exact = case.phi_exact(sn.mesh.centers, g)
        np.testing.assert_allclose(
            phi, exact, rtol=5e-3,
            err_msg=f"group {g} converged to the wrong limit",
        )
        for x_face in (0.0, case.slab_length):
            amplitude = float(np.abs(case.A(np.array([x_face]), g)[0]))
            _require(
                amplitude > 0.1 * float(case.c_groups[g]),
                f"A_{g}({x_face}) = {amplitude:.3e} is nearly zero, so this "
                f"fixture is effectively VACUUM at that face and the whole "
                f"module's premise is void",
            )


# ═══════════════════════════════════════════════════════════════════════
# G4 / G8 — the two routes, and the spec contract
# ═══════════════════════════════════════════════════════════════════════


def test_the_declared_and_supplied_channels_are_one_float_program() -> None:
    r"""⭐ ``from_mesh_laws(declared)`` is BIT-identical to ``case.prescribed_inflow``.

    The law tier and the channel tier must not merely agree — they must be the
    same computation, since the declaration's only job is to reach the same
    assembler. Anything less than ``array_equal`` here means the two routes have
    diverged, and a tolerance would hide exactly that.

    No solve: this compares the assembled ``q_∂`` directly, so it is the cheapest
    row in the module and the first to look at when G1 reddens.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    mesh0 = case.build_mesh(40)
    specs = {face: _face_spec(case, face, x) for face, x in _faces(case)}
    mesh = replace(
        mesh0,
        bc_left=PrescribedInflow(source=specs["xmin"]),
        bc_right=PrescribedInflow(source=specs["xmax"]),
    )
    sn = SNMesh(mesh, case.quadrature, case.materials)

    declared = AngularBoundarySourceSink.from_mesh_laws(sn)
    supplied = case.prescribed_inflow(SNMesh(mesh0, case.quadrature, case.materials))

    for face, _ in _faces(case):
        rows = np.asarray(sn.angular_trace.inflow_indices_for_face(face))
        d = np.asarray(declared.face_view(face))
        s = np.asarray(supplied.face_view(face))
        np.testing.assert_array_equal(
            d[rows], s[rows],
            err_msg=(
                f"the DECLARED and SUPPLIED q_∂ differ on {face}'s inflow rows "
                f"— the law tier and the channel tier are not one program"
            ),
        )
        off = np.ones(d.shape[0], dtype=bool)
        off[rows] = False
        _require(bool(off.any()), f"{face} has no off-inflow row; claim vacuous")
        np.testing.assert_array_equal(d[off], 0.0)


def test_the_spec_is_asked_exactly_ONCE_and_for_gamma_minus_ITSELF() -> None:
    r"""⭐ RE-POSED at P6 — this row asserted the opposite until 2026-08-06.

    It read ``test_the_spec_receives_both_shapes_and_the_probe_is_non_zero``
    and pinned that the spec is called at TWO shapes: a rank-1 ``(N,)``
    realize-time probe and the rank-2 ``(|Γ₋|, ng)`` materialisation. It also
    pinned that this spec's probe answer is non-zero, because
    ``assert_source_lives_on_incoming_trace`` returned early on a zero probe
    and a source could therefore decline the certification it was meant to
    pass.

    ⛔ **Both of those were properties of a defect, not of a contract.** Rank
    was the only signal distinguishing the probe call from the delivery call —
    "an accident, not a contract", as the design said — and a source that
    answered the probe with zeros skipped the ERR-047 check and delivered
    anyway (`[M]` zeros at ``(N,)``, ``7.0`` at ``(|Γ₋|, ng)``, delivered
    ``7.0``). P6 removed the probe. There is now ONE call, at ONE object, and
    that object is :math:`\Gamma_-(f)` itself.

    So the row now pins the *absence* of what it used to require: exactly one
    evaluation, against the inflow space of the declaring face, with the
    delivered values matching what the spec returns for that space. A future
    edit that reintroduced a second call — a probe, a shape sniff, a
    speculative evaluation — reddens here.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    mesh0 = case.build_mesh(20)
    spec = _face_spec(case, "xmin", 0.0)

    mesh = replace(mesh0, bc_left=PrescribedInflow(source=spec))
    sn = SNMesh(mesh, case.quadrature, case.materials)
    q = AngularBoundarySourceSink.from_mesh_laws(sn)

    gamma_minus = sn.angular_trace.inflow_space("xmin")
    _require(
        spec.spaces_seen == [gamma_minus.name],
        f"the spec was asked for {spec.spaces_seen}; expected exactly one "
        f"call, for {gamma_minus.name!r}. More than one means a probe or a "
        f"shape sniff is back; a different name means it was handed the wrong "
        f"face or the wrong directional tier.",
    )
    # …and what it returned for that space is what was delivered, so the single
    # call is the REAL one rather than a discarded rehearsal.
    rows = np.asarray(sn.angular_trace.inflow_indices_for_face("xmin"))
    np.testing.assert_array_equal(
        np.asarray(q.face_view("xmin"))[rows],
        spec.evaluate(gamma_minus),
        err_msg="the delivered q_∂ is not what the spec returned for Γ₋(xmin)",
    )


# ═══════════════════════════════════════════════════════════════════════
# G5 (value half) + the fixture's own structural precondition
# ═══════════════════════════════════════════════════════════════════════


def test_the_krylov_path_reproduces_the_manufactured_solution() -> None:
    r"""⭐ "Test the matvec path as well" — and it is stronger than ``A(0) = 0``.

    Reproducing the manufactured solution through GMRES with a declared law
    installed says the whole matvec is right, not merely that it is linear at
    the origin. It is also the row that could not have existed before P3: an
    affine ``A(x) = A_lin(x) − c`` breaks the Arnoldi relation
    ``A V_k = V_{k+1} H_k``, so this path RAISED ``ConvergenceCertificateError``
    (``‖Aψ − q‖/‖q‖ = 1.718``) rather than returning a wrong answer.

    One mesh, not a ladder: the convergence ORDER is a property of the spatial
    discretization and is already pinned on SI, and SI ≡ Krylov at the fixed
    point. What is Krylov-specific is that the fixed point is reached at all.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    sn, solution, _ = _declared_solve(case, 40, "krylov")
    for g in range(case.n_groups):
        phi = np.asarray(solution.scalar_flux.values)[g, :]
        exact = case.phi_exact(sn.mesh.centers, g)
        np.testing.assert_allclose(
            phi, exact, rtol=1e-2,
            err_msg=(
                f"[krylov] group {g} did not reproduce the manufactured "
                f"solution through the declared channel"
            ),
        )


def test_B_is_the_zero_morphism_when_every_face_is_prescribed() -> None:
    r"""The MMS fixture's own structural precondition — and a returning-affine tripwire.

    Both slab faces declare :class:`PrescribedInflow`, so since P3 both realize
    to the zero morphism and the assembled ``B`` is identically zero: `[M]`
    ``|B(x)|_inf = 0.0`` for a random ``x``. Every ordinate of the answer's
    inflow trace therefore comes from ``q_∂`` and from nowhere else, which is
    what lets G1 read the declared source straight off the converged flux.

    ⚠ **Deliberately NOT called a linearity gate.** On this fixture ``B(0) = 0``
    and ``B(2x) = 2B(x)`` hold because both sides are *structurally zero* — no
    input can make them red, so asserting them here would be a tautological
    companion guard (``vv`` Mode 8). Linearity is a genuine claim only where
    ``B`` is non-trivial, which is P5's prescribed + **reflective** fixture
    (`[M]` there: ``|B(x)|_inf = 1.320``, with an activation guard).

    What this row DOES catch is the affine operator returning: a ``B`` that
    emits ``q`` makes ``|B(x)|`` non-zero on exactly this fixture.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    sn, _, _ = _declared_solve(case, 20)

    from orpheus.sn.operators.boundary import SNBoundaryOperator
    from orpheus.numerics.operator import BlockRole

    B = SNBoundaryOperator(sn)
    _require(B.block_role is BlockRole.BOUNDARY, "B lost its block role")

    for face, _ in _faces(case):
        leaf = sn.bc[face]
        _require(
            leaf.block_role is BlockRole.BOUNDARY,
            f"the realized law on {face} is unstamped — before P3 the affine "
            f"operator carried no stamp, and its absence was (wrongly) believed "
            f"to fence it out of B",
        )

    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
    from orpheus.transport.timed_full_field import TimedFullField

    zero = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
    )
    rng = np.random.default_rng(20260806)
    x = replace(
        zero,
        interior=replace(
            zero.interior,
            values=rng.uniform(0.5, 2.0, zero.interior.values.shape),
        ),
        boundary=replace(
            zero.boundary,
            values=rng.uniform(0.5, 2.0, zero.boundary.values.shape),
        ),
    )
    out = np.asarray(B.apply(x).boundary.values)
    np.testing.assert_array_equal(
        out, 0.0,
        err_msg=(
            f"|B(x)|_inf = {float(np.max(np.abs(out))):.6e}, expected 0.0. With "
            f"every face prescribed, B must be identically the zero morphism — "
            f"a non-zero reading means an operator is emitting the source "
            f"again, which is ERR-075 returning."
        ),
    )

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

    ⚠ **P6 STOPGAP WITH AN OWNER.** Every constructor argument is information
    :meth:`~orpheus.geometry.boundary._source.InflowSourceSpec.evaluate` cannot
    supply, because its whole signature is a bare ``shape``:

    * ``mu_inflow`` — the per-row direction cosines **in trace order**. The spec
      is handed ``(|Γ₋|, ng)`` and must know that row ``i`` is ordinate
      ``inflow_indices_for_face(face)[i]``. Getting that wrong is silent: a
      constant-valued spec cannot detect a permuted or reversed trace.
    * ``A_g`` / ``B_g`` — evaluated at **this face's own coordinate**. The spec
      is not told which face it is on, so the coordinate is baked in and one
      spec instance is built per face.
    * ``W`` and the ``1/W`` convention — a property of the ansatz, invisible here.
    * ``n_ordinates`` — needed only to recognise the rank-1 probe (below).

    ⟹ P6's specification is "a source that receives the trace and the face",
    not "a source that receives a shape". See `scratch/p4_mms_design.md` §10.

    **The two evaluate shapes** (`[M]`, design §3):

    ``realize_boundary_law``  → ``evaluate((N,))``       — rank 1, ALL ordinates,
                                                            no group axis
    ``from_mesh_laws``        → ``evaluate((|Γ₋|, ng))`` — the real delivery

    A spec written for only one of these dies with ``IndexError`` before the MMS
    ever runs.

    ⛔ **The rank-1 probe MUST be non-zero.**
    ``assert_source_lives_on_incoming_trace`` opens with
    ``probe = source.evaluate((N,)); if not np.any(probe): return`` — so a spec
    that returns zeros there **silently skips the ERR-047 certification it is
    supposed to pass**. Measured: a spec returning zeros at rank-1 and ``7.0`` at
    rank-2 realizes cleanly. A presence predicate a source can decline is not a
    guard; that is a P6 finding in its own right, and this class works *with* the
    guard rather than around it.
    """

    def __init__(self, *, mu_inflow, A_g, B_g, W, n_ordinates) -> None:
        self.mu = np.asarray(mu_inflow, dtype=float)    # (|Γ₋|,) IN TRACE ORDER
        self.A = np.asarray(A_g, dtype=float)           # (ng,) at THIS face
        self.B = np.asarray(B_g, dtype=float)           # (ng,) at THIS face
        self.W = float(W)
        self.N = int(n_ordinates)
        #: Every shape this spec was asked for — G8 reads it.
        self.shapes_seen: list[tuple[int, ...]] = []

    def _trace(self) -> np.ndarray:
        # γ₋ψ is the angular FLUX on the trace, so it carries the ansatz's own
        # 1/W — a DIFFERENT reason from the bulk source's /W (design §1.3/§1.4).
        return (self.A[None, :] + self.mu[:, None] * self.B[None, :]) / self.W

    def evaluate(self, shape) -> np.ndarray:
        shape = tuple(int(s) for s in shape)
        self.shapes_seen.append(shape)
        if len(shape) == 1:
            _require(
                shape[0] == self.N,
                f"unexpected rank-1 probe {shape}; expected (N={self.N},). The "
                f"two evaluate shapes are (N,) at realize and (|Γ₋|, ng) at "
                f"materialisation — see this class's docstring.",
            )
            probe = np.zeros(shape, dtype=float)
            probe[: self.mu.size] = self._trace()[:, 0]
            return probe
        return self._trace()


# ═══════════════════════════════════════════════════════════════════════
# Fixture — declare on the GEOMETRY, solve through the PUBLIC entry
# ═══════════════════════════════════════════════════════════════════════


def _face_spec(case, sn, face: str, x_face: float) -> _ManufacturedFaceInflow:
    rows = np.asarray(sn.angular_trace.inflow_indices_for_face(face))
    return _ManufacturedFaceInflow(
        mu_inflow=case.quadrature.mu_x[rows],
        A_g=[float(np.atleast_1d(case.A(np.array([x_face]), g))[0])
             for g in range(case.n_groups)],
        B_g=[float(np.atleast_1d(case.B(np.array([x_face]), g))[0])
             for g in range(case.n_groups)],
        W=float(case.quadrature.weights.sum()),
        n_ordinates=case.quadrature.N,
    )


def _declared_solve(case, n_cells: int, inner: str = "source_iteration"):
    r"""Build → declare → solve, entirely through public surfaces.

    The trace is read from a **probe** ``SNMesh`` before the declaration is
    made, because the spec needs the inflow row order and only the trace knows
    it. That two-step is itself a P6 signal: a spec that received the trace
    would not need the probe. Deriving the ordering from the quadrature instead
    would be re-implementing production logic inside the test — the one thing an
    oracle must not do.

    Returns ``(sn, solution, specs)``.
    """
    mesh0 = case.build_mesh(n_cells)
    probe = SNMesh(mesh0, case.quadrature, case.materials)
    specs = {
        face: _face_spec(case, probe, face, x_face)
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
    probe = SNMesh(mesh0, case.quadrature, case.materials)
    specs = {
        face: _face_spec(case, probe, face, x)
        for face, x in _faces(case)
    }
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


def test_the_spec_receives_both_shapes_and_the_probe_is_non_zero() -> None:
    r"""⭐ The ERR-047 certification must actually RUN, not be opted out of.

    ``assert_source_lives_on_incoming_trace`` returns early on a zero rank-1
    probe, so a source can silently decline the certification it is meant to
    pass. This row pins that (a) both evaluate shapes really are requested —
    a refactor collapsing them would break every hand-rolled spec — and (b) this
    spec's probe is non-zero, so the guard fired.
    """
    case = build_slab_2g_nonvacuum_mms_case()
    mesh0 = case.build_mesh(20)
    probe_mesh = SNMesh(mesh0, case.quadrature, case.materials)
    spec = _face_spec(case, probe_mesh, "xmin", 0.0)

    mesh = replace(mesh0, bc_left=PrescribedInflow(source=spec))
    sn = SNMesh(mesh, case.quadrature, case.materials)
    AngularBoundarySourceSink.from_mesh_laws(sn)

    ranks = {len(s) for s in spec.shapes_seen}
    _require(
        ranks == {1, 2},
        f"the spec saw shapes {spec.shapes_seen} (ranks {sorted(ranks)}); both "
        f"the rank-1 (N,) realize-time probe and the rank-2 (|Γ₋|, ng) "
        f"materialisation call must occur",
    )
    n_inflow = int(np.asarray(
        sn.angular_trace.inflow_indices_for_face("xmin")).size)
    _require(
        (case.quadrature.N,) in spec.shapes_seen,
        f"no (N={case.quadrature.N},) probe in {spec.shapes_seen}",
    )
    _require(
        any(s[0] == n_inflow and len(s) == 2 for s in spec.shapes_seen),
        f"no (|Γ₋|={n_inflow}, ng) call in {spec.shapes_seen}",
    )
    _require(
        bool(np.any(spec.evaluate((case.quadrature.N,)))),
        "the rank-1 probe is all-zero, so assert_source_lives_on_incoming_trace "
        "returned early and the ERR-047 certification was SKIPPED",
    )

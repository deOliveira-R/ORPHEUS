r""":math:`\ker(L+C-S-B)` on a closed reflective box — the facts, from a DENSE SVD.

The sibling module :mod:`tests.sn.operators.test_loss_kernel_gauge` verifies the
**construction**: that :class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge`
builds the right subspace and projects onto it correctly. This module verifies
the **operator**, and the two are not the same question — the construction and
the counting law :eq:`dd-null-counting-law` share
:func:`~orpheus.sn.operators.loss_kernel_gauge._reflection_orbits` and
:func:`~orpheus.sn.operators.loss_kernel_gauge.gauge_freedom`, so a defect in
either moves *both* sides of every gate over there together. The only oracle
outside that shared root is a dense SVD of the assembled :math:`A`, and here it
is asked the four questions the gauge suite structurally cannot ask:

===============================================  ============================================
question                                          why the gauge suite cannot answer it
===============================================  ============================================
is :math:`A` really NON-singular where the        every ``_SINGULAR`` fixture over there is
predicate says so?                                singular; the ``dim == 0`` rows assert the
                                                  PREDICATE and the CONSTRUCTION, both of
                                                  which read the same two functions
is the kernel really PURE TRACE?                  the gauge is a trace-space endomorphism by
                                                  construction, so it cannot observe a bulk
                                                  component even if one existed
does the ``T`` term of                            every fixture there is ``level_symmetric``,
:eq:`dd-null-counting-law` hold?                  where :math:`T = 0` — the ``+\#\{`` tangential
                                                  trace DOFs ``\}`` term is unexercised
is :math:`\psi_{\rm exact}` the MINIMUM-norm       the canonicality theorem is the *licence*
member?                                           for the gauge; asserting it with the gauge's
                                                  own span would be circular
===============================================  ============================================

Promoted 2026-08-15 from ``derivations/diagnostics/diag_344_reflective_box_loss_nullspace.py``
(GitHub #344), which ``pyproject.toml``'s ``testpaths = ["tests"]`` never ran.
These are **characterization** gates in the `vv` sense — they pin what the
discretization does today so a change to the boundary closure or the spatial
scheme cannot move it silently — and every one of them is mutation-verified
(see the promotion report; each docstring names the mutation that reddens it).

⚠ **Fixture parity is load-bearing in two DIFFERENT ways here, and they are not
the same rule.** (a) Whether a *symmetric source* EXCITES the kernel is decided
by the first axis's parity — that rule lives in
:mod:`tests.sn.solve.test_every_entry_gauges_its_trace`, and it does not apply
to any gate in this file except
:func:`test_psi_exact_is_the_minimum_G_norm_member_of_the_manifold`, the only
one that solves. (b) Whether a *face-summed* mirror-odd functional can see a
kernel mode is decided by the TRANSVERSE cell count of that face — see
:func:`test_a_kernel_mode_carries_metric_but_no_NET_CURRENT`. ``[M]`` on
``(2,2)`` both counts are even and the tangential-current witness reads
``2.6e-15``; on ``(2,3) (3,2) (3,3) (3,4)`` it reads ``1.2e-01 .. 3.3e-01``.
**Never infer either from a mesh property without measuring it.**
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.face_layout import face_normal
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.operators.loss_kernel_gauge import predicted_kernel_dimension
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
from tests.sn._singular_loss_box import (
    LS4,
    REFLECTIVE as _R,
    VACUUM as _V,
    absorber,
    assemble,
    build,
    diamond_scattering_box,
    diamond_singular_box,
    drive,
    ld_reflective_box,
    loss_matvec,
    null_basis,
    rank_gap,
    tangential_dof_count,
    uniform_source_fixture,
)

pytestmark = pytest.mark.foundation

#: The canonical singular box for this file. ``(2,3)`` is the smallest mesh with
#: UNEQUAL cell counts — small because every gate here assembles a dense
#: :math:`A` (``n_dof = 768``, ``[M]`` 1.2 s), unequal because an x/y confusion
#: must have somewhere to show. ``[M]`` it carries ``dim ker A = 12``, the same
#: as ``(3,4)`` — the ``d = 2`` kernel dimension is mesh-INDEPENDENT
#: (:eq:`dd-null-counting-law`), which is exactly what licenses the shrink.
_CELLS = (2, 3)


@pytest.fixture(scope="module")
def singular_box():
    """``(sn_mesh, system, template, dense A, null basis, singular values)``.

    Shared through :func:`tests.sn._singular_loss_box.diamond_singular_box`, so
    the several gates that read this dense SVD pay for it once.
    """
    return diamond_singular_box(_CELLS)


# ─────────────────────────────────────────────────────────────────────
# 1. the operator-side oracle on BOTH sides of the predicate
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.verifies("dd-null-counting-law")
@pytest.mark.parametrize(
    "label,bcs,expected",
    [
        ("all reflective", [(_R, _R), (_R, _R)], 12),
        ("x vacuum", [(_V, _V), (_R, _R)], 0),
        ("y vacuum", [(_R, _R), (_V, _V)], 0),
        ("all vacuum", [(_V, _V), (_V, _V)], 0),
        ("xmin reflective / xmax vacuum", [(_R, _V), (_R, _R)], 0),
    ],
)
def test_the_singularity_needs_TWO_closed_axis_pairs(label, bcs, expected):
    r"""``dim ker A`` from a dense SVD, on both sides of the predicate.

    ⭐ **The four zero rows are the point.** Over in the gauge suite,
    ``dim == 0`` is asserted as ``LossKernelGauge.for_mesh(mesh).dimension == 0``
    and ``predicted_kernel_dimension(mesh) == 0`` — but the construction and the
    law BOTH gate on :func:`~orpheus.sn.operators.loss_kernel_gauge.gauge_freedom`
    and BOTH walk :func:`~orpheus.sn.operators.loss_kernel_gauge._reflection_orbits`,
    so they move together and agree on any wrong answer. Only an oracle outside
    that root can say the operator is *actually* non-singular there, and that
    matters in exactly one direction: a predicate that wrongly reported ABSENT
    would make the solver skip the gauge on a configuration that needs it, and
    no gate in the tree would notice.

    The last row is the sharp one: ``xmin`` reflective with ``xmax`` vacuum
    still gives the x axis a reflective law, so a predicate that counted FACES
    rather than closed PAIRS would call it two pairs. ``[M]`` widening the
    closed-pair criterion to *any* reflective face reddens **the mixed row and
    nothing else** in either module — 1 of 25.

    ✅ ``[M]`` **and that mutation is now ONE edit** — repaired 2026-08-15, in
    the commit that promoted this module. The criterion lives once, in
    :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.reflective_axes`, and
    :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.reflective_axis_pairs` is its
    ``len``; :func:`~orpheus.sn.operators.loss_kernel_gauge._reflection_orbits`
    reads the same property. ``[M]`` widening that ONE body — ``all(faces)`` to
    ``any(faces)`` — now reddens **4 of 83** gates across this module, the gauge
    suite and ``tests/sn/mesh/test_reflective_axis_pairs.py``: this row, plus
    ``test_one_vacuum_face_breaks_its_axis``,
    ``test_a_MIXED_axis_contributes_nothing`` and
    ``test_TWO_mixed_axes_are_worth_zero_not_one``. Baseline 83 passed.

    ⛔ *Kept as the record of what the promotion measured, per*
    ``plan-authoring`` *§3.* Until that repair the criterion had **two**
    line-for-line implementations — the mesh property and a private
    ``_reflective_axes`` in the gauge module, the same
    ``len(faces) == 2 and all(faces)`` test over the same ``by_axis`` loop, one
    returning the count and the other the axes. ``[M]`` widening **either alone
    was inert** (0 of 25 red) because the survivor guarded the gate; only
    widening both reddened this row. That measurement is what identified the
    twin — a mutation arm that reddens nothing is evidence about the SOURCE, not
    only about the gate (``vv`` #17).
    """
    if label == "all reflective":       # the shared build — assembled once
        sn_mesh, _system, _template, _dense, basis, singular = (
            diamond_singular_box(_CELLS))
    else:
        _mesh, system, template = build(_CELLS, bcs, absorber(2))
        sn_mesh = _mesh
        basis, singular = null_basis(assemble(loss_matvec(system), template))
    nullity = basis.shape[1]

    assert nullity == expected, (
        f"{label}: dim ker A = {nullity}, expected {expected} "
        f"(sigma_min/sigma_max = {singular[-1] / singular[0]:.3e})"
    )
    # …and the combinatorial law agrees, on BOTH sides.
    assert predicted_kernel_dimension(sn_mesh) == expected, (
        f"{label}: the counting law says "
        f"{predicted_kernel_dimension(sn_mesh)}, the SVD says {nullity}"
    )
    if expected:
        # The rank reading is not a threshold choice — it is a cliff.
        gap = rank_gap(singular, nullity)
        assert gap > 1e8, (
            f"{label}: no clean rank gap (ratio {gap:.3e}) — the integer "
            f"above is then an artefact of RANK_RTOL, not a property of A"
        )


# ─────────────────────────────────────────────────────────────────────
# 2. the kernel is a TRACE object — which is why the gauge can be one
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.verifies("dd-null-sawtooth")
def test_the_kernel_is_PURE_TRACE_and_carries_no_tangential_mass(singular_box):
    r"""No bulk component, and none on a tangential slot.

    ⭐ This is the premise that licenses
    :class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge` being an
    endomorphism of the boundary-trace space rather than of the full field, and
    **nothing else asserts it**: the gauge cannot observe a bulk component
    because it has no bulk rows, and
    ``test_the_bulk_and_keff_are_untouched`` asserts that *gauging* leaves the
    bulk alone — which is true by construction and says nothing about where
    :math:`\ker A` lives.

    The mechanism is :eq:`dd-null-sawtooth`: the mode drives
    :math:`\psi_c \equiv 0`, so the absorption term :math:`\Sigma_t V \psi_c`
    never sees it. ``[M]`` bulk null mass ``4.193e-29`` of dimension 12.

    The second leg pins the fixture's own precondition — ``level_symmetric``
    places no vanishing cosine (``[M]`` ``min |Omega.n| = 3.5002e-01``), so
    ``T`` is empty and the whole kernel is the component ``R`` the gauge
    builds. Without it the first leg would be silently reporting on a mixture
    of two different objects.

    ``[M]`` mutation: make the sawtooth sign ``+1`` at the far face regardless
    of :math:`(-1)^{n_a}` — the modes leave :math:`\ker A`, the SVD basis
    changes, and the bulk mass rises off ``1e-29``.
    """
    sn_mesh, _system, template, _dense, basis, _singular = singular_box
    assert basis.shape[1] > 0, "fixture is no longer singular — re-derive #344"
    n_bulk = template.interior.values.size

    # Basis-INDEPENDENT: diag(V V^T) is the null projector's diagonal, so this
    # is a statement about the SUBSPACE, not about the SVD's choice of basis.
    null_mass = np.einsum("ij,ij->i", basis, basis)
    bulk_mass = float(null_mass[:n_bulk].sum())
    assert bulk_mass < 1e-20, (
        f"the kernel is NOT pure-trace: bulk mass = {bulk_mass:.6e} of "
        f"dimension {basis.shape[1]} — LossKernelGauge is a trace-space "
        f"operator and would be unable to reach the rest"
    )
    assert float(null_mass[n_bulk:].sum()) == pytest.approx(
        basis.shape[1], rel=1e-9), "the null mass does not sum to the dimension"

    omega_dot_n = np.asarray(sn_mesh.angular_trace.omega_dot_n)
    assert tangential_dof_count(sn_mesh) == 0, (
        f"fixture precondition broken: level_symmetric grew a tangential "
        f"ordinate (min |Omega.n| = {np.min(np.abs(omega_dot_n)):.6e}), so "
        f"the kernel measured here is no longer pure R"
    )


# ─────────────────────────────────────────────────────────────────────
# 3. ⭐ the T term of the counting law — unexercised anywhere else
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.verifies("dd-null-counting-law")
@pytest.mark.parametrize(
    "label,quad,tangential,remainder",
    [
        ("level_symmetric(4)", LS4, 0, 12),
        ("product(4,4)", Quadrature.product(n_mu=4, n_phi=4), 160, 0),
        ("lebedev(11)", Quadrature.lebedev(order=11), 160, 18),
    ],
)
def test_the_T_plus_R_split_is_exact_per_quadrature(
        label, quad, tangential, remainder):
    r"""⭐ ``dim ker A == #{tangential trace DOFs} + R``, exactly.

    :eq:`dd-null-counting-law` ends in ``+ #{tangential trace DOFs}``, and
    **that term has no other test**: every fixture in the gauge suite's
    ``_SINGULAR`` list is ``level_symmetric``, where it is zero.
    :func:`~orpheus.sn.operators.loss_kernel_gauge.predicted_kernel_dimension`
    deliberately counts ``R`` only (its own docstring says so), so the law as
    published and the law as implemented differ by exactly this term.

    The three rows are the three structural cases, and all three are needed:

    * ``level_symmetric(4)`` — **pure R** (:math:`T = 0`). The kernel the gauge
      builds, and nothing else.
    * ``product(4,4)`` — **pure T** (:math:`R = 0`). Every ordinate has
      :math:`\mu_x = 0` or :math:`\mu_y = 0`, so no orbit has two active axes
      and the gauge correctly returns an EMPTY basis while ``A`` is 160-fold
      singular. A construction that "helpfully" spanned the tangential rows
      would divide by :math:`\sqrt{G} = 0`.
    * ``lebedev(11)`` — **both non-zero** (``160 + 18``). This is the only row
      on which ``T + R`` is an equation rather than two separate degenerate
      readings, and it is the row that makes the other two informative.

    ⛔ ``[M]`` **this gate would have caught a live documentation defect.** The
    gauge suite's module-docstring SVD table records
    ``product(4,4) | mine 224 | law 224 | svd 224`` and
    ``lebedev(11) | mine 242 | law 242 | svd 242`` at ``(3,4)``, while its own
    prose four lines below says ``product(4,4)`` is *"pure T (224 tangential,
    R = 0)"*. Measured at ``(3,4)``: ``gauge.dimension = 0`` and ``law = 0`` for
    ``product(4,4)``, ``18`` and ``18`` for ``lebedev(11)``, against
    ``svd = 224`` and ``242``. The table's ``mine``/``law`` columns recorded
    ``T + R`` where they claim ``R``. The prose was right and the table was
    wrong, and no executed row could tell them apart.

    ⚠ Cell counts are ``(2,3)``, so the two axes contribute UNEQUALLY to ``T``
    (``[M]`` 96 from the x-faces, 64 from the y-faces = 160): on a square mesh
    the two contributions are equal and a per-axis bookkeeping error cancels.
    """
    sn_mesh, system, template = build(_CELLS, [(_R, _R)] * 2, absorber(2),
                                      quad=quad)
    measured_tangential = tangential_dof_count(sn_mesh)
    dense = assemble(loss_matvec(system), template)
    basis, singular = null_basis(dense)
    nullity = basis.shape[1]

    assert measured_tangential == tangential, (
        f"{label}: {measured_tangential} tangential trace DOFs, "
        f"expected {tangential}"
    )
    assert nullity - measured_tangential == remainder, (
        f"{label}: dim ker A = {nullity}, T = {measured_tangential}, so "
        f"R = {nullity - measured_tangential} (expected {remainder})"
    )
    # R is what the gauge scopes to, and the law must report the same.
    assert predicted_kernel_dimension(sn_mesh) == remainder, (
        f"{label}: predicted_kernel_dimension = "
        f"{predicted_kernel_dimension(sn_mesh)} but the SVD's R is {remainder} "
        f"— the law and the operator disagree on the R/T boundary"
    )
    assert sn_mesh.loss_kernel_gauge.dimension == remainder, (
        f"{label}: the gauge spans {sn_mesh.loss_kernel_gauge.dimension} "
        f"dimensions but R is {remainder}"
    )
    gap = rank_gap(singular, nullity)
    assert gap > 1e8, f"{label}: no clean rank gap (ratio {gap:.3e})"


# ─────────────────────────────────────────────────────────────────────
# 4. which functionals can see a kernel mode — the R/T discriminator
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.verifies("sn-kernel-mirror-blindness")
def test_a_kernel_mode_carries_metric_but_no_NET_CURRENT(singular_box):
    r"""The three readings that decide what a gauge can even mean.

    1. :math:`\langle v, v\rangle_G > 0`. This is what separates ``R`` from
       ``T``: a tangential slot carries :math:`G = |\Omega\cdot\hat n|\,w_n = 0`
       **exactly**, so every value there has the same (zero) :math:`G`-norm and
       there is no minimum-norm representative to choose. ``R`` sits on
       current-carrying ordinates, so the minimum-:math:`G`-norm member is
       unique and the gauge is well-posed. ``[M]`` ``0.4636`` for a unit-norm
       mode.
    2. The face-summed **mirror-EVEN** moments are exactly zero —
       :eq:`sn-kernel-mirror-blindness`, on a basis that owes nothing to the
       gauge. ``[M]`` ``|Sum w psi| <= 2.1e-15`` and
       ``|Sum w |Omega.n| psi| <= 1.0e-15``.
    3. The **net current** :math:`J = \sum_n (\Omega_n\cdot\hat n) w_n \psi_n`
       is zero too (``[M]`` ``<= 2.6e-16``) — the same character-orthogonality
       identity with :math:`F = w_n \mu_a`, which is ODD in
       :math:`\mathrm{sign}\,\mu_a` while the mode's face-:math:`a` values are
       EVEN in it. This is the *operator-level* root of the negative control in
       ``test_the_spurious_TANGENTIAL_current_along_a_mirror_is_gone``: the
       normal currents cannot move because no member of the manifold carries
       one.

    ⭐ **And the activation leg, without which 2 and 3 are satisfied by the zero
    vector**: the current TANGENTIAL to a face is mirror-EVEN in both signs and
    survives — ``[M]`` ``1.9134e-01``.

    ⚠ **That leg has its own parity precondition, and it is NOT the source
    parity rule.** The face sum runs over the transverse cells against the
    checkerboard :math:`(-1)^{\sum_{c \neq a} i_c}`, so it cancels whenever the
    transverse cell count is EVEN. ``[M]`` on ``(2,2)`` — both counts even —
    the tangential current reads ``2.5650e-15`` and this leg is **inert**;
    ``(2,3) (3,2) (3,3) (3,4)`` read ``1.9134e-01 / 3.2691e-01 / 1.1985e-01 /
    2.4252e-01``. Here the x-faces carry it (transverse ``n_y = 3``, odd).

    ``[M]`` mutation: drop the transverse checkerboard from
    :func:`~orpheus.sn.operators.loss_kernel_gauge._transverse_factors` — the
    modes leave :math:`\ker A`, the SVD basis changes, and legs 2/3 red.
    """
    sn_mesh, _system, template, _dense, basis, _singular = singular_box
    n_bulk = template.interior.values.size
    coefficients = np.random.default_rng(0).standard_normal(basis.shape[1])
    mode = basis @ (coefficients / np.linalg.norm(coefficients))
    field = type(template).from_flat(mode, template)

    metric = np.asarray(sn_mesh.angular_trace.partial_current_metric)
    g_norm = float(np.sqrt(np.sum(metric * mode[n_bulk:] ** 2)))
    assert g_norm > 1e-2, (
        f"the kernel is metric-annihilated after all: <v,v>_G^(1/2) = "
        f"{g_norm:.6e} — that would put it in the tangential class T, where "
        f"no minimum-G-norm representative exists and the gauge is undefined"
    )

    weights = np.asarray(sn_mesh.quad.weights, dtype=float)
    cosines = (np.asarray(sn_mesh.quad.mu_x, dtype=float),
               np.asarray(sn_mesh.quad.mu_y, dtype=float))
    omega_dot_n = np.asarray(sn_mesh.angular_trace.omega_dot_n)
    worst_tangential = 0.0
    for index, (name, slot) in enumerate(
            sn_mesh.angular_trace.layout.faces.items()):
        axis, _sign = face_normal(name)
        per_ordinate = slot.slice_view(
            mode[n_bulk:]).reshape(slot.shape[0], -1).sum(1)
        for moment, weighting in (
                ("|Omega.n|^0", weights),
                ("partial current", weights * np.abs(omega_dot_n[index]))):
            value = float(weighting @ per_ordinate)
            assert abs(value) < 1e-11, (
                f"{name}: the mirror-EVEN {moment} moment reads {value:.6e} on "
                f"a kernel mode — eq. sn-kernel-mirror-blindness says it is "
                f"exactly zero, and the whole no-scalar-summary-reveals-it "
                f"argument rests on that"
            )
        net = float(np.max(np.abs(field.boundary.net_current(name))))
        assert net < 1e-11, (
            f"{name}: a kernel mode carries net current max|J| = {net:.6e}; "
            f"if a manifold member can carry one, the normal-current negative "
            f"control in tests/sn/solve/test_every_entry_gauges_its_trace.py "
            f"is measuring an accident"
        )
        worst_tangential = max(worst_tangential, abs(
            float((weights * cosines[1 - axis]) @ per_ordinate)))

    assert worst_tangential > 1e-2, (
        f"no face's TANGENTIAL current sees the kernel (worst "
        f"{worst_tangential:.3e}) — every assertion above is then satisfied by "
        f"the zero vector. Check the transverse cell counts: the face sum "
        f"cancels the checkerboard whenever they are all EVEN (cells={_CELLS})"
    )


# ─────────────────────────────────────────────────────────────────────
# 5. ⭐ the canonicality theorem — why min-norm is not a convention
# ─────────────────────────────────────────────────────────────────────
@pytest.mark.verifies("sn-kernel-mirror-blindness")
def test_psi_exact_is_the_minimum_G_norm_member_of_the_manifold():
    r"""⭐ The licence for the whole disposition, asserted rather than argued.

    A singular operator leaves a solution *manifold*, so "return the
    minimum-:math:`G`-norm member" is a convention unless the physical answer
    happens to BE that member. It does: :math:`\psi_{\rm exact}` is
    :math:`G`-orthogonal to :math:`\ker A` (``[M]`` ``1.551e-15`` of its own
    :math:`G`-norm), which follows from :eq:`sn-kernel-mirror-blindness` —
    :math:`\psi_{\rm exact}` is flat here, hence mirror-EVEN, hence annihilated
    by every kernel mode. **Nothing in the tree asserted this**; the gauge suite
    shows the projection recovers the exact trace, which is the *consequence*,
    not the licence.

    The second leg then closes the loop with a projector built from the **SVD**
    null basis — so both sides of the recovery are independent of
    :class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelGauge`, unlike the
    flagship in the gauge suite, where the span under test is the gauge's own.
    ``[M]`` ``1.0485e-01 -> 1.178e-12``.

    ⚠ This is the one gate in this file that SOLVES, so it is the one the
    source-excitation parity rule applies to: ``n_x = 3`` is ODD, and with the
    uniform isotropic source an even first axis leaves the Gauss-Seidel
    deviation at round-off and the second leg vacuous. The ``before`` assertion
    is the guard.

    ⛔ **What the first leg CANNOT see, measured** (`vv` Mode 12 — compute the
    stabiliser, do not assume one). Orthogonality here is a *character* identity,
    so it survives every metric that is CONSTANT ACROSS AN ORBIT'S CELLS. On this
    fixture, ``frac``:

    ===============================================  ==============
    trace metric                                      ``frac``
    ===============================================  ==============
    shipped :math:`|\Omega\cdot\hat n|\,w`             ``1.5514e-15``
    x2 everywhere                                     ``1.5137e-15``
    x(1 + ½ sign :math:`\Omega\cdot\hat n`)            ``1.5135e-15``
    x(1 + ½ sign of the PARTNER face's)               ``1.7500e-15``
    one whole face x3                                 ``1.8743e-15``
    random per-DOF                                    **``1.9377e-02``**
    ===============================================  ==============

    So this leg is **not** a gate on the metric's values — a whole battery of
    plausible metric errors is designed-green on it, and only a perturbation
    that breaks the orbit-block structure reddens it. It is a gate on the
    KERNEL's parity: a construction whose modes were not mirror-odd would fail
    it under the shipped metric. Say that rather than claiming metric
    sensitivity the row does not have.

    ``[M]`` in the promotion battery the row reds under the damped-closure
    positive control and under every arm that breaks the Gauss-Seidel solve —
    on the ACTIVATION leg and the RECOVERY leg respectively, never on ``frac``.
    """
    cells = (3, 2)
    sn_mesh, system, template, n_dof, n_bulk, source, exact = (
        uniform_source_fixture(cells))
    dense = assemble(loss_matvec(system), template)
    basis, _singular = null_basis(dense)
    assert basis.shape[1] > 0, "fixture is no longer singular"

    # G on the full field: quadrature weight x cell volume in the bulk, the
    # partial-current metric on the trace — the trace space's own atoms.
    weights = np.asarray(sn_mesh.quad.weights, dtype=float)
    volumes = np.asarray(sn_mesh.volumes, dtype=float)
    metric = np.empty(n_dof)
    bulk = np.zeros(template.interior.values.shape)
    bulk[...] = (weights.reshape((weights.size, 1) + (1,) * (bulk.ndim - 2))
                 * volumes.reshape((1, 1) + template.interior.values.shape[2:]))
    metric[:n_bulk] = bulk.ravel()
    metric[n_bulk:] = np.asarray(sn_mesh.angular_trace.partial_current_metric)

    root = np.sqrt(metric)
    orthonormal, _r = np.linalg.qr(root[:, None] * basis)

    def project(vector):
        """The G-orthogonal projection onto ker A, via the SVD basis."""
        return (orthonormal @ (orthonormal.T @ (root * vector))) / root

    in_kernel = float(np.sqrt(np.sum(metric * project(exact) ** 2))
                      / np.sqrt(np.sum(metric * exact ** 2)))
    assert in_kernel < 1e-12, (
        f"psi_exact is NOT G-orthogonal to ker A ({in_kernel:.6e}) — the "
        f"minimum-G-norm member would then be a CONVENTION rather than the "
        f"physical answer, and gauging would be picking one"
    )

    iterate = drive(system, sn_mesh, template, source, "gauss_seidel")
    before = float(np.max(np.abs(
        (iterate[n_bulk:] - exact[n_bulk:]) / exact[n_bulk:])))
    after = float(np.max(np.abs(
        ((iterate - project(iterate))[n_bulk:] - exact[n_bulk:])
        / exact[n_bulk:])))
    assert before > 1e-2, (
        f"the driver landed on the canonical member already ({before:.3e}) — "
        f"the recovery leg proves nothing. n_x must be ODD for the uniform "
        f"isotropic source to excite the kernel (cells={cells})"
    )
    assert after < 1e-10, (
        f"projecting off the SVD null basis did NOT recover the exact trace: "
        f"{before:.6e} -> {after:.6e}"
    )


# ─────────────────────────────────────────────────────────────────────
# 6. the closure, not the geometry, is what creates the kernel
# ─────────────────────────────────────────────────────────────────────
def test_a_damping_closure_leaves_the_IDENTICAL_box_non_singular():
    r"""Linear-discontinuous on the same all-reflective box: ``dim ker A == 0``.

    The gauge suite asserts this through the predicate and the construction
    (``gauge_freedom`` says DAMPED, ``dimension == 0``); this asserts it of the
    **operator**, which is the form the theory page's remedy hierarchy actually
    claims — *"switch to a closure that damps the face mode"* is a claim about
    :math:`A`, not about what a predicate reports.

    It is also the substitution that PROVES the mechanism rather than arguing
    it: same box, same quadrature, same materials, same boundary laws — only
    :eq:`dd-null-sawtooth`'s :math:`\psi_{\rm out} = -\psi_{\rm in}` is no
    longer compatible with a zero cell moment. ``[M]`` diamond gives 12 here,
    LD gives 0.

    ⚠ The LD assembly is ``[M]`` 10.8 s (three times the DOFs), and
    :func:`tests.sn.solve.test_boundary_gs_is_a_coherent_splitting.test_kernel_free_configs_give_the_SAME_TRACE_under_both_schedules`
    needs the identical object as its precondition. They share ONE assembly
    through :func:`tests.sn._singular_loss_box.ld_reflective_box` rather than
    computing it twice — the two make different claims about one measurement.
    """
    _dd_mesh, _dd_sys, _dd_tmpl, _dd_dense, dd_basis, _dd_s = (
        diamond_scattering_box(_CELLS))
    assert dd_basis.shape[1] == 12, (
        f"the diamond control is no longer singular (dim ker A = "
        f"{dd_basis.shape[1]}) — the contrast below then proves nothing"
    )

    sn_mesh, _system, _template, _dense, basis, singular = (
        ld_reflective_box(_CELLS))
    assert isinstance(sn_mesh.scheme, LinearDiscontinuous), (
        f"the LD fixture is closed by {type(sn_mesh.scheme).__name__} — this "
        f"row is then comparing diamond with diamond"
    )
    assert sn_mesh.reflective_axis_pairs == 2, (
        "the geometry conjunct must still be satisfied, or this row is about "
        "the BOX and not about the CLOSURE"
    )
    assert basis.shape[1] == 0, (
        f"linear-discontinuous left dim ker A = {basis.shape[1]} on the box "
        f"where diamond gives 12 (sigma_min/sigma_max = "
        f"{singular[-1] / singular[0]:.3e}) — the closure is the mechanism, so "
        f"this is the claim the remedy hierarchy rests on"
    )

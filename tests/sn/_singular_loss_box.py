r"""Builders for the box on which ``A = L + C - S - B`` is exactly singular.

Two test modules characterise that configuration from opposite tiers —
:mod:`tests.sn.operators.test_loss_nullspace_reflective_box` (facts about
:math:`\ker A` itself) and
:mod:`tests.sn.solve.test_boundary_gs_is_a_coherent_splitting` (facts about the
splitting and the iteration that runs on it). They need the *same* fixtures, so
the fixtures live here once (coding-elegance Pattern 2): two copies of
``_apply_A`` would be two spellings of the operator under test, and the day one
of them drifted the two modules would be characterising different objects while
both stayed green.

⚠ **This module is NOT collected by pytest, so a bare ``assert`` here would be
stripped by the canonical ``python -O`` runner** (`vv` Mode 8 — the hazard's
real domain is exactly non-collected helpers like this one). It therefore
contains **no assertions at all**: it constructs, and every claim is made in a
collected test module.

The operator is assembled from the **production** splitting — ``A =
system.implicit_operator - sum(system.explicit_gains)`` — so the object these
tests characterise is the one the SI driver iterates, not a re-derivation of it
(`vv` §structural independence: the independent oracle is LAPACK's SVD, not a
second hand-written matvec).
"""

from __future__ import annotations

from typing import Any, Callable, Optional

import numpy as np
from numpy.typing import NDArray

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import (
    SNSolver,
    _as_sn_mesh,
    _build_fixed_source_rhs,
    _select_si_splitting,
    _unwindowed_cold_start,
)
from orpheus.transport.mesh.axis import AxisMesh

VACUUM, REFLECTIVE = BC("vacuum"), BC("reflective")

#: The default rule for every fixture here. ``level_symmetric`` places its
#: cosines on shells :math:`|\mu| \ge \mu_1 > 0`, so it has **zero tangential
#: trace DOFs** and the whole kernel is the component ``R`` the gauge builds —
#: which is what makes a dimension comparison an equality rather than an
#: inequality. The tangential-bearing rules are used deliberately, and only, by
#: the ``T + R`` split gate.
LS4 = Quadrature.level_symmetric(sn_order=4)

#: Relative singular-value cut for reading a rank off a dense SVD. Not a taste
#: parameter: ``[M]`` the gap between the last retained and the first discarded
#: singular value measures ``1.9e+14`` (``level_symmetric(4)``), ``9.8e+13``
#: (``product(4,4)``) and ``4.3e+13`` (``lebedev(11)``) on the ``(2,3)`` box, so
#: any cut in ``[1e-13, 1e-2]`` returns the same integer. The gap is *asserted*,
#: not assumed — see ``test_the_T_plus_R_split_is_exact_per_quadrature``.
RANK_RTOL = 1e-10

#: Extents are deliberately UNEQUAL so a fixture's two axes are never
#: interchangeable — an x/y confusion has somewhere to show.
EXTENTS = (1.0, 2.0, 3.0)


def axes(cells, bcs):
    """``AxisMesh`` per axis, on :data:`EXTENTS`, with the given BC pairs."""
    return tuple(
        AxisMesh(edges=np.linspace(0.0, extent, n + 1), bc_low=low, bc_high=high)
        for extent, n, (low, high) in zip(EXTENTS, cells, bcs)
    )


def absorber(ng: int = 2, sig_t=(0.8, 1.6)):
    r"""A pure absorber — ``[M]`` ``S`` is the **ZERO** operator on it.

    ⚠ That is a config-blindness fact, not an aside: ``[M]``
    ``|S x| = 0.000000e+00`` for a random ``x`` here versus ``6.741510e+00`` on
    :func:`scatterer`. Any claim about where the scattering gain lands in the
    splitting is **vacuous** on this mixture, so the ``M - N == A`` gate uses
    :func:`scatterer` instead. This mixture is for the fixtures whose exact
    answer must be analytic (:func:`uniform_source_fixture`).
    """
    total = np.asarray(sig_t[:ng], dtype=float)
    return make_mixture(
        sig_t=total, sig_c=total.copy(), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=np.zeros((ng, ng)),
    )


def scatterer(ng: int = 2, c: float = 0.5, sig_t=(0.8, 1.6)):
    """A within-group scatterer with scattering ratio ``c`` — ``S != 0``."""
    total = np.asarray(sig_t[:ng], dtype=float)
    sig_s = np.zeros((ng, ng))
    for group in range(ng):
        sig_s[group, group] = c * total[group]
    return make_mixture(
        sig_t=total, sig_c=total - sig_s.sum(axis=1), sig_f=np.zeros(ng),
        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=sig_s,
    )


def build(cells, bcs, mixture, quad=LS4, scheme=None):
    """``(sn_mesh, within_group_system, cold-start template)``."""
    sn_mesh = _as_sn_mesh(axes(cells, bcs), quad, {0: mixture}, scheme=scheme)
    solver = SNSolver(sn_mesh, inner_solver="source_iteration")
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    return sn_mesh, system, _unwindowed_cold_start(sn_mesh, history_depth=0)


def loss_matvec(system) -> Callable[[Any], Any]:
    r""":math:`A x = (\text{implicit} - \sum \text{gains})\,x` — the PRODUCTION
    composition the SI driver iterates (:math:`M - N`)."""

    implicit = system.implicit_operator
    gains = list(system.explicit_gains)

    def apply(x):
        out = implicit.apply(x)
        for gain in gains:
            out = out - gain.apply(x)
        return out

    return apply


def gains_matvec(gains) -> Callable[[Any], Any]:
    r""":math:`N x = \sum_i g_i x` for a splitting's explicit gains."""

    def apply(x):
        out = None
        for gain in gains:
            contribution = gain.apply(x)
            out = contribution if out is None else out + contribution
        return out

    return apply


def assemble(apply_fn, template) -> NDArray:
    """Dense matrix of ``apply_fn`` by unit probes — exhaustive, one column each."""
    n_dof = template.to_flat().size
    dense = np.empty((n_dof, n_dof))
    unit = np.zeros(n_dof)
    for column in range(n_dof):
        unit[column] = 1.0
        dense[:, column] = apply_fn(
            type(template).from_flat(unit, template)).to_flat()
        unit[column] = 0.0
    return dense


def null_basis(dense: NDArray) -> tuple[NDArray, NDArray]:
    """``(orthonormal null basis as columns, the singular values)``."""
    _left, singular, right = np.linalg.svd(dense)
    return right[singular < RANK_RTOL * singular[0], :].T, singular


#: Assembled operators, keyed by configuration label.
#:
#: A dense ``A`` costs ``[M]`` 1.2 s on the ``(2,3)`` diamond box and **10.8 s**
#: on the linear-discontinuous one (3x the DOFs), and two modules make different
#: claims about the SAME LD configuration — the operator-tier module that a
#: damping closure leaves the box non-singular, the solve-tier module that the
#: two schedules therefore agree on it. Assembling it twice would be a
#: Pattern-2 duplicate that happens to cost 11 s; the cache makes the two
#: consumers share one computation while each states its own claim.
#:
#: ⚠ Correct under a mutation battery because a plugin installs its mutation at
#: configure time, before any test runs — one process sees one production.
_ASSEMBLED: dict[str, tuple] = {}


def assembled_loss(label: str, factory):
    """``(sn_mesh, system, template, dense A, null basis, singular values)``.

    ``factory`` is a zero-argument callable returning ``build(...)``'s triple;
    it runs at most once per ``label`` per session. Consumers must only READ
    the returned arrays.
    """
    if label not in _ASSEMBLED:
        sn_mesh, system, template = factory()
        dense = assemble(loss_matvec(system), template)
        basis, singular = null_basis(dense)
        _ASSEMBLED[label] = (sn_mesh, system, template, dense, basis, singular)
    return _ASSEMBLED[label]


#: The two configurations both modules assemble. Labels, not builders, are the
#: shared surface — a consumer names the configuration it means and cannot
#: silently ask for a slightly different one.
def diamond_singular_box(cells=(2, 3)):
    """The canonical singular box: diamond difference, all reflective, ng=2."""
    return assembled_loss(
        f"diamond all-reflective {cells} absorber",
        lambda: build(cells, [(REFLECTIVE, REFLECTIVE)] * len(cells),
                      absorber(2)),
    )


def diamond_scattering_box(cells=(2, 3)):
    """The same box with ``c = 0.5`` — the CLOSURE contrast's control.

    Same mixture as :func:`ld_reflective_box`, so the only thing that differs
    between them is the spatial closure. (The kernel is cross-section
    independent, so this is belt-and-braces — but a contrast argued as
    *"the IDENTICAL box"* should not have a second variable in it.)
    """
    return assembled_loss(
        f"diamond all-reflective {cells} c=0.5",
        lambda: build(cells, [(REFLECTIVE, REFLECTIVE)] * len(cells),
                      scatterer()),
    )


def ld_reflective_box(cells=(2, 3)):
    """The IDENTICAL box closed by linear discontinuous — non-singular."""
    from orpheus.transport.spatial.linear_discontinuous import (
        LinearDiscontinuous,
    )

    return assembled_loss(
        f"LD all-reflective {cells} c=0.5",
        lambda: build(cells, [(REFLECTIVE, REFLECTIVE)] * len(cells),
                      scatterer(), scheme=LinearDiscontinuous()),
    )


def rank_gap(singular: NDArray, nullity: int) -> float:
    """Ratio of the first discarded to the last retained singular value.

    ``inf`` when nothing was discarded — the honest reading for a
    full-rank operator, where there is no rank decision to justify.
    """
    if nullity == 0 or nullity >= singular.size:
        return float("inf")
    return float(singular[-nullity - 1] / singular[-nullity])


def tangential_dof_count(sn_mesh) -> int:
    r"""``#{trace DOFs on ordinates with :math:`\Omega\cdot\hat n = 0`}``.

    The ``T`` term of the counting law :eq:`dd-null-counting-law`. Rows and
    columns of ``A`` are identically zero there, and the trace metric
    :math:`G = |\Omega\cdot\hat n|\,w_n` vanishes, so they are a *different*
    object from the component ``R`` the gauge builds.
    """
    omega_dot_n = sn_mesh.angular_trace.omega_dot_n
    layout = sn_mesh.angular_trace.layout
    total = 0
    for index, face in enumerate(layout.faces):
        slot = layout.faces[face]
        per_ordinate = int(np.prod(slot.shape[1:]))
        total += int(np.sum(np.abs(omega_dot_n[index]) < 1e-14)) * per_ordinate
    return total


def uniform_source_fixture(cells, ng: int = 2):
    r"""All-reflective box + uniform isotropic source, whose exact answer is flat.

    :math:`\psi = Q / (\Sigma_t \sum w)` everywhere, bulk **and** trace — an
    analytic reference that owes nothing to any operator in this campaign.

    Returns ``(sn_mesh, system, template, n_dof, n_bulk, source, exact)``.
    """
    sn_mesh, system, template = build(
        cells, [(REFLECTIVE, REFLECTIVE)] * len(cells), absorber(ng),
    )
    n_dof = template.to_flat().size
    n_bulk = template.interior.values.size
    total_weight = float(sn_mesh.quad.weights.sum())
    per_group = 1.0 / (total_weight * np.asarray(
        sn_mesh.materials[0].SigT, dtype=float))

    interior = np.zeros(template.interior.values.shape)
    interior[:] = 1.0 / total_weight
    source = _build_fixed_source_rhs(np.asarray(interior, dtype=float), sn_mesh)

    exact = np.zeros(n_dof)
    block = np.zeros(template.interior.values.shape)
    for group in range(ng):
        block[:, group, ...] = per_group[group]
    exact[:n_bulk] = block.ravel()
    for slot in sn_mesh.angular_trace.layout.faces.values():
        face = np.zeros(slot.shape)
        for group in range(ng):
            face[:, group, ...] = per_group[group]
        exact[n_bulk + slot.offset:
              n_bulk + slot.offset + slot.flat_size] = face.ravel()
    return sn_mesh, system, template, n_dof, n_bulk, source, exact


def isotropic_source(sn_mesh, template):
    """A group-graded isotropic source — ``(1 + g)/W`` per ordinate."""
    total_weight = float(sn_mesh.quad.weights.sum())
    interior = np.zeros(template.interior.values.shape)
    for group in range(template.interior.values.shape[1]):
        interior[:, group, ...] = (1.0 + group) / total_weight
    return _build_fixed_source_rhs(np.asarray(interior, dtype=float), sn_mesh)


def select_splitting(system, sn_mesh, schedule: str):
    r"""``(M, (N_i,))`` for one schedule — the production selector, narrowed.

    :func:`~orpheus.sn.solver._select_si_splitting` is typed for the SEEDLESS
    record, and ``build_within_group_system`` returns the union that also covers
    the seed-carrying (coupled) arm. Production narrows the same way, with an
    explicit guard and a loud ``TypeError``
    (``orpheus/sn/solver.py``, ``_within_group_si``) — mirrored here rather
    than suppressed, so a record that ever arrives in the other shape says so
    instead of mis-binding.
    """
    from orpheus.sn.operators.streaming import StreamingCollisionOperator
    from orpheus.transport.operators.scattering import ScatteringOperator

    implicit = system.implicit_operator
    if not isinstance(implicit, StreamingCollisionOperator):
        raise TypeError(
            f"these fixtures are all SEEDLESS, so the record's implicit "
            f"operator must be the StreamingCollisionOperator; got "
            f"{type(implicit).__name__} — a seed-carrying mesh never reaches "
            f"_select_si_splitting at all"
        )
    scattering, n2n, boundary = system.explicit_gains
    if not isinstance(scattering, ScatteringOperator):
        raise TypeError(
            f"the seedless record's first gain must be the ScatteringOperator "
            f"(the builder's (S, N2N, B_a) convention); got "
            f"{type(scattering).__name__}"
        )
    implicit, boundary_gain = _select_si_splitting(
        implicit, boundary, sn_mesh, schedule)
    # The selector decides the BOUNDARY half only; the gains are named here
    # the way the driver names them (S, N₂ₙ, boundary gain — §14.1, B LAST).
    return implicit, (scattering, n2n, boundary_gain)


def drive_recorded(system, sn_mesh, template, source, schedule: str,
                   tol: float = 1e-13, initial: Optional[NDArray] = None):
    r"""``(flat iterate, IterationRecord)`` from the SI **driver**.

    ⚠ This is the driver, **not** a public entry: nothing here gauges the
    returned trace. Every claim built on it is a claim about the un-gauged
    iterate — see each caller's docstring for why that is the tier it means.
    """
    from orpheus.numerics.iteration import SourceIteration

    n_dof = template.to_flat().size
    base, gains = select_splitting(system, sn_mesh, schedule)
    iteration = SourceIteration(base.inverse(), *gains,
                                max_iter=400_000, tol=tol)
    start = np.zeros(n_dof) if initial is None else initial
    solution, record = iteration.solve(
        source, initial_guess=type(template).from_flat(start, template),
    )
    if not record.converged:
        raise RuntimeError(
            f"the {schedule} driver did not converge on this fixture — every "
            f"claim downstream of it would be about a truncated iterate"
        )
    return solution.to_flat(), record


def drive(system, sn_mesh, template, source, schedule: str,
          tol: float = 1e-13, initial: Optional[NDArray] = None) -> NDArray:
    """The flat converged iterate — :func:`drive_recorded` without the record."""
    return drive_recorded(system, sn_mesh, template, source, schedule,
                          tol=tol, initial=initial)[0]


def both_drivers(sn_mesh, system, template, tol: float = 1e-13):
    """``{schedule: flat iterate}`` from the SAME zero cold start."""
    source = isotropic_source(sn_mesh, template)
    return {
        schedule: drive(system, sn_mesh, template, source, schedule, tol=tol)
        for schedule in ("gauss_seidel", "jacobi")
    }

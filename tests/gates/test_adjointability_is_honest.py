r"""The adjointability law (ruling 27 of ``numerics_untangling.md``): the census over every consumer.

The law and its two legs are stated once, in :mod:`tests.gates.numerics._adjointability`, which also builds the
numerics rows and the leg-1 check. This file builds the consumers' rows (fission, diffusion, SN), runs leg 1 on
them, and leg 2 and the registry's reach on the union (layer DAG Phase 2, ruling 28).

The reds at landing (2026-10-18) were the fix's done-when. Six rows were
strict ``xfail`` against #614 (the diffusion operators whose composite ends
speak flat arrays). The two leg-1 rows (``diffusion_assembled_leakage``,
``diffusion_resolvent``) passed, reciprocity included, once the composite
space's metric legs took its flat coordinates (space-membership U0, #666's flat
half): ``.H`` is ``♯ ∘ dual ∘ ♭`` on the flat chart. Four leg-2 rows remain
strict ``xfail`` (:data:`_CARRIER_BLOCKED_LEG2`): their flag reads ``False``
while the assembled-matrix witness realises the transpose, so a row turns red
the day the flag is made honest. The witness route of each leg-2 row is in its parametrize id
(``<row>--<route>``) and in :data:`_WITNESS_ROUTES`.

``-O``-safe: every assertion is a ``pytest.fail`` call.

``-O``-safe: every assertion is a ``pytest.fail`` call.
"""

from __future__ import annotations

from dataclasses import field, replace
from typing import Any, Callable

import numpy as np
import pytest
from scipy import sparse

from orpheus.geometry import CoordSystem
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.numerics.assembled_operator import (
    MissingAssembly,
    SparseAssembledOperator,
    assemblable,
    assemble,
)
from orpheus.numerics.coupled_system import (
    CoupledField,
    CoupledOperator,
    CoupledSpace,
)
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import LinearOperator, TensorProductOperator, posed
from orpheus.numerics.operator import RieszLowerOperator
from orpheus.numerics.operator import RieszRaiseOperator
from tests.gates.numerics._adjointability import NUMERICS_REGISTRY, Posed, assert_doors, flat, like, rel

pytestmark = [
    pytest.mark.foundation,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_operator.py::test_adjoint_weighted_identity",
        "tests/gates/numerics/test_riesz_legs.py::"
        "test_lowering_realizes_the_metric_pairing",
        "tests/gates/numerics/test_riesz_legs.py::"
        "test_the_adjoint_composite_routes_through_the_legs",
        "tests/gates/numerics/test_operator_capability_predicates.py::"
        "test_inverse_adjoint_contract_keystone_v2",
        "tests/gates/sn/operators/test_g_adjoint_reciprocity.py::"
        "test_g_adjoint_reciprocity_full_block",
        "tests/gates/sn/operators/test_g_adjoint_reciprocity.py::"
        "test_full_loss_g_adjoint_reciprocity",
    ),
]

#: The Green splitting's convergence tolerance, and its reciprocity band
#: (10 × the tolerance, amplified by the conditioning of the het slab).
_GREEN_TOL = 1e-13
_GREEN_RTOL = 1e-9


# ═══════════════════════════════════════════════════════════════════════
# The registry
# ═══════════════════════════════════════════════════════════════════════


def _sn(case: str):
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    return G._BUILDERS[case]()


def _sn_composite(sn, seed: int = 2026):
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    return G._random_composite(sn, np.random.default_rng(seed))


def _lc(sn, sig_t):
    from orpheus.sn.operators.streaming import StreamingCollisionOperator, StreamingOperator
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    # S2 (ruling 52): ``+`` never fuses; the sweep atom is built explicitly.
    return StreamingCollisionOperator(StreamingOperator.pose(sn), MultiplicationOperator.from_mesh(
        sig_t, sn,
    ))


def _source_like(sn, sig_t, seed: int):
    r"""A source-classed composite: the loss applied to a random flux."""
    return _lc(sn, sig_t).apply(_sn_composite(sn, seed))


# ── fission: the energy dyad and its isotropic operator ──


def _fission_dyad_tensor() -> Posed:
    r"""``|χ⟩⟨νΣ_f| ⊗ I``, the factor shape of every fission operator, as the
    isotropic fission's own kernel: ends ``None``.

    The dyad's transpose is coded, so ``.dual()`` opens; the tensor is unbound
    and not a multiplier, so ``.H`` is honestly unposed. The defect of the
    first fix read the dyad's flag ``False`` here, and every fission adjoint
    built over it raised (`[M]` 2026-10-18, qa)."""
    from tests.gates.sn.operators import test_isotropic_fission as IF

    fission = IF._op()
    return Posed(fission.kernel, np.ones(fission.domain.shape))


def _isotropic_fission() -> Posed:
    from tests.gates.sn.operators import test_isotropic_fission as IF

    op = IF._op()
    return Posed(op, np.ones(op.domain.shape))


# ── diffusion: the leaves, the loss, the flat view, the resolvent ──


def _diffusion():
    from tests.gates.diffusion import test_operators as D

    mesh = D._diffusion_mesh(MirrorDeck(), WallCondition.vacuum())
    return D, mesh


def _assembled_sandwich(op: LinearOperator, template: Any) -> Callable[[Any], Any]:
    r"""The witness ``♯_V ∘ Mᵀ ∘ ♭_W`` from ``op``'s assembled matrix ``M``.

    The transpose is the sparse matrix's own (``M.T``); the metrics are the
    operator's end spaces' (``apply_metric`` / ``apply_inverse_metric``). No
    operator-level transpose is written by hand.
    """
    if not assemblable(op):
        raise MissingAssembly(f"{type(op).__name__} does not assemble")
    M = assemble(op).matrix
    V, W = op.domain, op.codomain
    if V is None or W is None:
        raise RuntimeError(f"fixture drift: {type(op).__name__} has an unbound end")

    def witness(y: Any) -> Any:
        lowered = flat(W.apply_metric(y))
        transposed = type(template).from_flat(np.asarray(M.T @ lowered), template)
        return V.apply_inverse_metric(transposed)

    return witness


def _diffusion_leaf(kind: str) -> Callable[[], Posed]:
    def build() -> Posed:
        from orpheus.diffusion.operators import (
            DiffusionBoundaryOperator,
            LeakageOperator,
        )

        D, mesh = _diffusion()
        x = D._random_flux(mesh)
        if kind == "leakage":
            op: LinearOperator = LeakageOperator(mesh)
        elif kind == "boundary":
            op = DiffusionBoundaryOperator(mesh)
        else:
            op = D._loss(mesh, mesh.mat_xs)
        return Posed(op, x, witness=_assembled_sandwich(op, x))

    return build


def _diffusion_flat_pairing(template: Any) -> Callable[[Any, Any], float]:
    space = template.space

    def pair(a: Any, b: Any) -> float:
        return float(space.inner_product(
            type(template).from_flat(np.asarray(a), template),
            type(template).from_flat(np.asarray(b), template),
        ))

    return pair


def _diffusion_flattened() -> Posed:
    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    A = D._loss(mesh, mesh.mat_xs)
    op = FlattenedOperator(A, x)
    sandwich = _assembled_sandwich(A, x)
    pair = _diffusion_flat_pairing(x)
    return Posed(
        op, x.to_flat(), pair_v=pair, pair_w=pair,
        witness=lambda y: flat(sandwich(type(x).from_flat(np.asarray(y), x))),
    )


def _diffusion_flattened_fission() -> Posed:
    r"""The diffusion pencil's flattened fission, ``FlattenedOperator(F, template)``, the production spelling.

    Its inner transposes, so the adapter does (the flat chart is an isometry of
    coordinates). Until 2026-10-31 the adapter read the base's ``False`` while
    its ``apply_transpose`` worked; the posing's annihilation check, once it
    spelled the transpose through the ``.dual()`` door, refused it
    (``tests/gates/diffusion/test_fixed_source_verdict.py`` D6 reddened). The
    flattened LOSS row beside this one stays #614-blocked: its inner does not
    transpose.
    """
    from orpheus.diffusion.solver import DiffusionSolver

    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    op = DiffusionSolver(mesh).pencil.rhs
    pair = _diffusion_flat_pairing(x)
    return Posed(op, x.to_flat(), pair_v=pair, pair_w=pair)


def _diffusion_resolvent() -> Posed:
    r"""The production diffusion resolvent, ``diffusion/solver.py``'s spelling."""
    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    op = MatrixInverseOperator(FlattenedOperator(D._loss(mesh, mesh.mat_xs), x))
    pair = _diffusion_flat_pairing(x)
    return Posed(op, x.to_flat(), pair_v=pair, pair_w=pair)


def _diffusion_assembled() -> Posed:
    from orpheus.diffusion.operators import LeakageOperator

    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    pair = _diffusion_flat_pairing(x)
    return Posed(
        assemble(LeakageOperator(mesh)), x.to_flat(), pair_v=pair, pair_w=pair,
    )


# ── SN: the leaves and composites on a 1-D and a 2-D pose ──


def _sn_streaming(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        from orpheus.sn.operators.streaming import StreamingOperator

        sn, _ = _sn(case)
        x = _sn_composite(sn)
        return Posed(StreamingOperator.pose(sn), x, y=x)

    return build


def _sn_streaming_collision(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        sn, sig_t = _sn(case)
        x = _sn_composite(sn)
        return Posed(_lc(sn, sig_t), x, y=x)

    return build


def _sn_sweep(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        sn, sig_t = _sn(case)
        return Posed(_lc(sn, sig_t).inverse(), _source_like(sn, sig_t, 7))

    return build


def _sn_boundary(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        from orpheus.sn.operators.boundary import SNBoundaryOperator

        sn, _ = _sn(case)
        x = _sn_composite(sn)
        return Posed(SNBoundaryOperator(sn), x, y=x)

    return build


def _sn_within_group_loss(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

        sn, sig_t = _sn(case)
        x = _sn_composite(sn)
        return Posed(G._loss_operator(sn, sig_t), x, y=x)

    return build


def _sn_full_loss(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

        sn, A, _ = G._FULL_LOSS_BUILDERS[case]()
        x = _sn_composite(sn)
        return Posed(A, x, y=x)

    return build


def _sn_masked_boundary_lower() -> Posed:
    r"""``B_lower`` of the Gauss–Seidel split on the reflective 2-D box.

    Witness ``(P B)ᴴ = Bᴴ P``: ``B_lower = P B`` with ``P`` the coordinate
    projection onto the emitted inflow rows, self-adjoint under the diagonal
    trace metric ``|Ω·n| w``.
    """
    from orpheus.sn.loss_representation.sweep_schedule import (
        SweepSchedule,
        reflective_faces,
    )
    from orpheus.sn.operators.boundary import SNBoundaryOperator

    sn, _ = _sn("cart2d_2g")
    B = SNBoundaryOperator(sn)
    lower = B.split(
        SweepSchedule.gauss_seidel(sn.ndim, sn.quad.octants, reflective_faces(sn)),
    ).lower
    if not any(np.size(r) for r in lower.rows.values()):
        raise RuntimeError("fixture drift: B_lower emits no rows on cart2d_2g")

    def project(y: Any) -> Any:
        values = np.zeros_like(y.boundary.values)
        kept = replace(y.boundary, values=values)
        for face, rows in lower.rows.items():
            kept.face_view(face)[rows] = y.boundary.face_view(face)[rows]
        return replace(
            y, interior=replace(y.interior, values=np.zeros_like(y.interior.values)),
            boundary=kept,
        )

    BH = B.H
    x = _sn_composite(sn)
    return Posed(lower, x, y=x, witness=lambda y: BH.apply(project(y)))


def _sn_coupled_substitution(case: str) -> Callable[[], Posed]:
    r"""The 1×1 block grid ``[[L + C]]`` and its substitution inverse."""

    def build() -> Posed:
        sn, sig_t = _sn(case)
        space = CoupledSpace.from_systems([sn.full_field_space])
        grid = CoupledOperator([[_lc(sn, sig_t)]], domain=space, codomain=space)
        x = CoupledField(systems=(_source_like(sn, sig_t, 5),))
        return Posed(grid.inverse(), x)

    return build


def _sn_inverse_of_leaf() -> Posed:
    r"""``InverseOperator(C)`` of the pointwise collision leaf.

    Witness ``(C⁻¹)ᴴ = (Cᴴ)⁻¹ = C⁻¹``: ``C`` is a real field multiplier BOUND
    to the slab's spaces, whose metrics are diagonal, so ``Cᴴ = C`` there and
    its inverse is self-adjoint. (Since the S6a review ``C`` is not
    metric-free: a field multiplier commutes only with diagonal metrics, #625;
    the precondition is now the diagonal ends, asserted below.)
    """
    from orpheus.numerics.metric import DiagonalMetric, FactoredMetric
    from orpheus.numerics.operator import InverseOperator
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    sn, sig_t = _sn("slab_2g")
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    for end in (C.domain, C.codomain):
        m = getattr(end, "_resolved_metric", "absent")
        blocks = [f for _, f in m.entries] if isinstance(m, FactoredMetric) else [m]
        if not all(f is None or isinstance(f, DiagonalMetric) for f in blocks):
            raise RuntimeError(f"fixture drift: C's end {end!r} is not diagonal-metric")
    op = C.inverse()
    if type(op) is not InverseOperator:
        raise RuntimeError(f"fixture drift: C.inverse() is {type(op).__name__}")
    return Posed(op, _source_like(sn, sig_t, 3), witness=op.apply)


def _sn_inverse_of_product() -> Posed:
    r"""``InverseOperator((L + C) @ C)``, the #285 product inverse.

    Witness ``((L+C) C)⁻ᴴ = (L+C)⁻ᴴ C⁻ᴴ = (L+C)⁻¹.H ∘ C⁻¹`` (the product law,
    the swap law on the sweep, and the pointwise law on ``C``).
    """
    from orpheus.numerics.operator import InverseOperator
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    sn, sig_t = _sn("slab_2g")
    lc = _lc(sn, sig_t)
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    op = (lc @ C).inverse()
    if type(op) is not InverseOperator:
        raise RuntimeError(f"fixture drift: the product inverse is {type(op).__name__}")
    sweep_h = lc.inverse().H
    c_inv = C.inverse()
    return Posed(
        op, _source_like(sn, sig_t, 4),
        witness=lambda y: sweep_h.apply(c_inv.apply(y)),
    )


def _sn_green() -> Posed:
    r"""``GreenOperator(L + C, S)`` (the named split, since S2) on the het 2-G slab.

    Witness ``((L+C) − S)⁻ᴴ = ((L+C)ᴴ − Sᴴ)⁻¹``: the sum law, then the Green
    of the adjoint split ``GreenOperator((L+C)ᴴ, Sᴴ)`` (its implicit part inverts by the swap law
    onto the 1-D reverse sweep). The witness is INDEPENDENT of the row's
    ``.H``: that is the base ``♯ ∘ G.dual() ∘ ♭``, the Green of the DUAL sum
    between the Riesz legs. Until the revision of 2026-10-18 the Green
    overrode ``.H`` with this very Green of the adjoint sum, and the witness
    then agreed with it to exactly 0.0, a tautology (the elegance review).
    """
    sn, sig_t, lc, S, op = _green_case()
    adjoint = GreenOperator(lc.H, S.H, max_iter=4000, tol=_GREEN_TOL)
    # Both ends are source → flux, so the adjoint's input is source-classed.
    # No row threads a seed: the Green's cold start is its first Neumann
    # partial sum, in the iterate's own class.
    return Posed(
        op, _source_like(sn, sig_t, 6), y=_source_like(sn, sig_t, 8),
        witness=adjoint.apply, rtol=_GREEN_RTOL,
    )


def _green_case():
    r"""``GreenOperator(L + C, S)`` on the het 2-G slab, with its pose."""
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    sn, _, S = G._full_loss_case(CoordSystem.CARTESIAN, 2)
    sig_t = np.asarray(CrossSectionField.total_of(sn.mat_xs).values, dtype=float)
    lc = _lc(sn, sig_t)
    return sn, sig_t, lc, S, GreenOperator(lc, S, max_iter=4000, tol=_GREEN_TOL)


def _sn_green_composed() -> Posed:
    r"""``C @ G``: a Green inside a product (qa, 2026-10-18: ``.H.apply`` on
    ``C@G``, ``G@C``, ``2·G`` and ``G+G`` raised while the Green read
    ``True``, because the Green coded no ``apply_transpose`` and its cold start
    was a source-classed zero)."""
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    sn, sig_t, _, _, G = _green_case()
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    return Posed(C @ G, _source_like(sn, sig_t, 9), y=_source_like(sn, sig_t, 10), rtol=_GREEN_RTOL)


def _sn_bulk_analysis() -> Posed:
    r"""``BulkAnalysisOperator`` = ``lift(M) ⊕ Id_Γ`` on the 2-D full-loss box.

    Witness ``(M ⊕ I)ᴴ = Mᴴ ⊕ I``: the frame analysis face's own ``.H`` on the
    moment bulk, the trace passed through.
    """
    from orpheus.sn.operators.windowing import BulkAnalysisOperator
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    sn, _, S = G._full_loss_case_cart2d()
    face = S.flux_analysis
    op = BulkAnalysisOperator(face, sn.full_field_space)
    face_h = face.H
    x = _sn_composite(sn)

    def witness(y: Any) -> Any:
        # The face's ``.H`` reads the bare moment array (see the face's row).
        pulled = np.asarray(face_h.apply(np.asarray(y.interior.values)))
        return replace(
            x, interior=replace(x.interior, values=pulled), boundary=y.boundary,
        )

    return Posed(op, x, witness=witness)


def _sn_flux_analysis_face(case: str) -> Callable[[], Posed]:
    r"""The scattering operator's minted flux-analysis face ``M`` (the face
    :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` lifts)."""

    def build() -> Posed:
        from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

        if case == "cart2d_2g":
            sn, _, S = G._full_loss_case_cart2d()
        else:
            sn, _, S = G._FULL_LOSS_BUILDERS[case]()
        from orpheus.transport.fields.angular_flux import AngularFlux

        face = S.flux_analysis
        x = _sn_composite(sn).interior
        if not isinstance(x, AngularFlux):
            raise RuntimeError(f"fixture drift: the bulk is {type(x).__name__}")
        # ``face.H`` reads a bare moment array: it refuses the face's own typed
        # codomain element (``[M]`` 2026-10-18: TypeError on cart2d, a
        # broadcast ValueError on the slab), a carrier-typing defect of the
        # metric's Riesz leg (#615), not a flag lie. ``apply`` admits only the typed leaf.
        return Posed(face, x, y=np.asarray(face.apply(x).values))

    return build


def _sn_fission(case: str) -> Callable[[], Posed]:
    def build() -> Posed:
        from orpheus.transport.operators.fission import FissionOperator

        sn, _ = _sn(case)
        x = _sn_composite(sn)
        op = FissionOperator.from_solver_data(mat_xs=sn.mat_xs, space=sn.full_field_space)
        return Posed(op, x, y=x)

    return build


def _sn_coupled_substitution_carrying() -> Posed:
    r"""The carrying sphere's implicit grid ``M`` and its block substitution:
    the radial block's transposed solve is on the diagonal (qa, 2026-10-18:
    no row reached it, so its flag reddened nothing)."""
    from orpheus.sn.splitting import Splitting, resolve_schedule
    from tests.gates.sn.architecture._config import (
        random_state,
        record_for,
        sphere_carrying,
    )

    problem = sphere_carrying()
    record = record_for(problem)
    M = Splitting.from_schedule(record, resolve_schedule(problem, "jacobi")).implicit
    # The adjoint equation's right-hand side is source-classed, as for the
    # Green row: the transposed substitution subtracts the off-diagonal
    # seeding's transpose from it, which a flux-classed sample refuses.
    return Posed(
        M.inverse(), M.apply(random_state(record, seed=31)),
        y=M.apply(random_state(record, seed=37)),
    )


def _sn_sweep_dual() -> Posed:
    r"""``((L+C)⁻¹)ᵀ``, the dual arrow of the slab sweep: its ends are the
    duals of composite (bulk ⊕ trace) spaces, and its ``.H`` is the base door
    over them. The reciprocity reads the operator's OWN end metrics, i.e. the
    composite duals' inverse forms (S6a: a composite's dual is the composite
    of its members' duals). Until S6a this row paired through a hand-built
    ``G⁻¹`` because the composite dual was Euclidean (census: the dual metric
    moved a sample by 0.0); `[M]` 2026-10-21 the row reds when the composite
    dual is restored to a member-less ``DualSpace``."""
    sn, sig_t = _sn("slab_2g")
    sweep = _lc(sn, sig_t).inverse()
    x = _sn_composite(sn)
    return Posed(
        sweep.dual(), x, y=_source_like(sn, sig_t, 12),
    )


def _sn_green_cart2d() -> Posed:
    r"""A Green over the 2-D head ``(L+C)``: the 2-D twin of :func:`_sn_green`.

    Until #616 the 2-D head had no reverse sweep and both doors refused (qa:
    no row held a Green over a 2-D head, so dropping its head clause reddened
    nothing). Since #616 the dual split's implicit part inverts by the swap
    law onto the multi-D reverse sweep, so both doors open, with the same
    independent witness as the slab row: the Green of the adjoint split
    ``GreenOperator((L+C)ᴴ, Sᴴ)``. ``y`` is source-classed, the class the
    Green's adjoint takes (a Green maps source to flux); the refusal-era row
    passed none, and the default (the forward output's class, flux) is the
    wrong class for this operator."""
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    sn, _, S = G._full_loss_case_cart2d()
    sig_t = np.asarray(CrossSectionField.total_of(sn.mat_xs).values, dtype=float)
    lc = _lc(sn, sig_t)
    op = GreenOperator(lc, S, max_iter=4000, tol=_GREEN_TOL)
    adjoint = GreenOperator(lc.H, S.H, max_iter=4000, tol=_GREEN_TOL)
    return Posed(
        op, _source_like(sn, sig_t, 13), y=_source_like(sn, sig_t, 14),
        witness=adjoint.apply, rtol=_GREEN_RTOL,
    )


def _sn_dsa_correction() -> Posed:
    from orpheus.sn.acceleration import DSACorrection

    sn, _ = _sn("slab_2g")
    x = _sn_composite(sn)
    return Posed(DSACorrection.from_problem(sn), x, y=x)


#: The registry: row id → builder. One row per operator the gate judges.
#: The consumers' rows; the numerics rows are imported from the numerics file, so the census is their union.
_CONSUMER_REGISTRY: dict[str, Callable[[], Posed]] = {
    "fission_dyad_tensor_unbound": _fission_dyad_tensor,
    "isotropic_fission_bound": _isotropic_fission,
    # diffusion
    "diffusion_leakage": _diffusion_leaf("leakage"),
    "diffusion_boundary": _diffusion_leaf("boundary"),
    "diffusion_loss": _diffusion_leaf("loss"),
    "diffusion_flattened_loss": _diffusion_flattened,
    "diffusion_flattened_fission": _diffusion_flattened_fission,
    "diffusion_assembled_leakage": _diffusion_assembled,
    "diffusion_resolvent": _diffusion_resolvent,
    # SN, 1-D and 2-D poses
    "sn_streaming_slab": _sn_streaming("slab_2g"),
    "sn_streaming_cart2d": _sn_streaming("cart2d_2g"),
    "sn_streaming_collision_slab": _sn_streaming_collision("slab_2g"),
    "sn_streaming_collision_cart2d": _sn_streaming_collision("cart2d_2g"),
    "sn_sweep_slab": _sn_sweep("slab_2g"),
    "sn_sweep_cyl_product": _sn_sweep("cyl_product_2g"),
    "sn_sweep_ld_slab": _sn_sweep("ld_slab_2g"),
    "sn_sweep_cart2d": _sn_sweep("cart2d_2g"),
    "sn_sweep_ld_2d": _sn_sweep("ld_2d_2g"),
    "sn_boundary_slab": _sn_boundary("slab_2g"),
    "sn_boundary_cart2d": _sn_boundary("cart2d_2g"),
    "sn_within_group_loss_slab": _sn_within_group_loss("slab_2g"),
    "sn_within_group_loss_cart2d": _sn_within_group_loss("cart2d_2g"),
    "sn_full_loss_slab": _sn_full_loss("slab_2g"),
    "sn_full_loss_sphere": _sn_full_loss("sphere_2g"),
    "sn_full_loss_cart2d": _sn_full_loss("cart2d_2g"),
    "sn_masked_boundary_lower_cart2d": _sn_masked_boundary_lower,
    "sn_coupled_substitution_slab": _sn_coupled_substitution("slab_2g"),
    "sn_coupled_substitution_cart2d": _sn_coupled_substitution("cart2d_2g"),
    "sn_inverse_of_leaf": _sn_inverse_of_leaf,
    "sn_inverse_of_product": _sn_inverse_of_product,
    "sn_green": _sn_green,
    "sn_green_composed": _sn_green_composed,
    "sn_green_cart2d": _sn_green_cart2d,
    "sn_sweep_slab_dual": _sn_sweep_dual,
    "sn_fission_slab": _sn_fission("slab_2g"),
    "sn_fission_cart2d": _sn_fission("cart2d_2g"),
    "sn_coupled_substitution_carrying_sphere": _sn_coupled_substitution_carrying,
    "sn_flux_analysis_face_slab": _sn_flux_analysis_face("slab_2g"),
    "sn_flux_analysis_face_cart2d": _sn_flux_analysis_face("cart2d_2g"),
    "sn_bulk_analysis_cart2d": _sn_bulk_analysis,
    "sn_dsa_correction_slab": _sn_dsa_correction,
}
#: The census's population: the numerics rows (the helper's) and the consumers' rows.
_REGISTRY: dict[str, Callable[[], Posed]] = {**NUMERICS_REGISTRY, **_CONSUMER_REGISTRY}

#: Leg 2's rows: row id → the law its witness realises. A row absent here has
#: no coded transpose route (report-only for leg 2; see the module docstring).
_WITNESS_ROUTES: dict[str, str] = {
    "diffusion_leakage": "assembled-sandwich",
    "diffusion_boundary": "assembled-sandwich",
    "diffusion_loss": "assembled-sandwich",
    "diffusion_flattened_loss": "assembled-sandwich",
    "sn_masked_boundary_lower_cart2d": "BH-P",
    "sn_inverse_of_leaf": "pointwise-self-adjoint-inverse",
    "sn_inverse_of_product": "product-swap-sweepH-Cinv",
    "sn_green": "green-of-adjoint-sum",
    "sn_green_cart2d": "green-of-adjoint-sum",
    "sn_bulk_analysis_cart2d": "MH-direct-sum-I",
}


# ═══════════════════════════════════════════════════════════════════════
# Leg 1 — the flag predicts the behaviour
# ═══════════════════════════════════════════════════════════════════════


#: Rows blocked by the diffusion carrier defect (#614): composite-ended
#: operators whose flag reads False while the assembled-matrix witness realises
#: the transpose. Strict, so a row turns red the day the flag is made honest and
#: the mark cannot go stale. (The two leg-1 rows, the flag True while ``.H.apply``
#: raised on the flat carrier, passed in space-membership U0: the composite
#: space's metric legs take its flat coordinates.)
_CARRIER_BLOCKED_LEG2 = frozenset({
    "diffusion_leakage", "diffusion_boundary", "diffusion_loss", "diffusion_flattened_loss",
})


def _rows(names, blocked, ids=None):
    marked = pytest.mark.xfail(
        strict=True,
        reason="#614: composite-ended diffusion operators speak flat arrays, so their adjoint cannot be realised",
    )
    ids = ids or list(names)
    return [
        pytest.param(name, id=i, marks=marked) if name in blocked else pytest.param(name, id=i)
        for name, i in zip(names, ids)
    ]


@pytest.mark.parametrize("row", list(_CONSUMER_REGISTRY))
@pytest.mark.catches("ERR-117")
def test_each_transpose_door_opens_iff_predicted(row: str) -> None:
    r"""``.dual()`` works ⟺ ``is_adjointable``; ``.H`` works ⟺ ``is_adjointable ∧ posed``.

    A refused door raises ``MissingAdjoint`` at construction. Where a door
    opens, the sanity leg: Euclidean reciprocity for ``.dual()``, and
    ``⟨A x, y⟩_W = ⟨x, A.H y⟩_V`` in the operator's own end metrics for ``.H``.
    """
    assert_doors(row, _CONSUMER_REGISTRY[row]())


# ═══════════════════════════════════════════════════════════════════════
# Leg 2 — a transpose the tree can realise is not hidden
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "row",
    _rows(
        list(_WITNESS_ROUTES), _CARRIER_BLOCKED_LEG2,
        [f"{row}--{route}" for row, route in _WITNESS_ROUTES.items()],
    ),
)
@pytest.mark.catches("ERR-117")
def test_a_realisable_transpose_is_advertised(row: str) -> None:
    r"""A row whose witness route (a law of the algebra) realises the adjoint
    reads ``is_adjointable is True``.

    The witness is checked first, by reciprocity against the forward; a broken
    witness fails as a FIXTURE defect, never as the flag's.
    """
    posed = _REGISTRY[row]()
    if posed.witness is None:
        pytest.fail(f"[{row}] registry drift: listed in _WITNESS_ROUTES with no witness")
    rng = np.random.default_rng(4096)
    x = like(posed.x, rng)
    y = posed.codomain_sample(rng)
    lhs = posed.inner_w(posed.act(x), y)
    rhs = posed.inner_v(x, posed.witness(y))
    if not rel(lhs, rhs) < posed.rtol:
        pytest.fail(
            f"[{row}] FIXTURE: the witness route {_WITNESS_ROUTES[row]!r} is not "
            f"the adjoint (rel {rel(lhs, rhs):.2e}); repair the witness, not the flag"
        )
    if not posed.op.is_adjointable:
        pytest.fail(
            f"[{row}] {type(posed.op).__name__}.is_adjointable=False while the "
            f"route {_WITNESS_ROUTES[row]!r} realises its adjoint (reciprocity "
            f"rel {rel(lhs, rhs):.1e})"
        )


def test_the_registry_reaches_every_lying_class() -> None:
    r"""Input count and a positive control on the registry (X1, X2): the rows
    build, and each class the census named is present at least once."""
    names = {type(_REGISTRY[r]().op).__name__ for r in _REGISTRY}
    required = {
        "SweepOperator", "RankOneOperator", "RieszLowerOperator",
        "RieszRaiseOperator", "CoupledSubstitutionOperator", "LeakageOperator",
        "DiffusionBoundaryOperator", "FlattenedOperator", "InverseOperator",
        "GreenOperator", "DSACorrection", "SNMaskedBoundaryOperator",
        "BulkAnalysisOperator", "SparseAssembledOperator", "StreamingOperator",
        "SNBoundaryOperator", "OperatorSum", "MatrixInverseOperator",
        "TensorProductOperator", "IsotropicFission", "FissionOperator",
        "OperatorProduct",
    }
    missing = required - names
    if len(_REGISTRY) < 45 or missing:
        pytest.fail(
            f"registry of {len(_REGISTRY)} rows misses {sorted(missing)} "
            f"(built classes: {sorted(names)})"
        )


def test_the_inverse_swaps_the_forwards_transpose_verbs_on_the_sn_sweep() -> None:
    r"""``(A⁻¹)ᵀ`` solves with ``Aᵀ`` on the SN slab sweep: the inverse's transposed solve IS the forward's
    transpose, and it exists iff the forward transposes. Split from numerics'
    ``test_adjointability_is_honest.py::test_the_inverse_swaps_the_forwards_transpose_verbs[sweep_slab]`` (layer DAG
    Phase 2 refresh): the numerics rows (a diagonal, a solve-only leaf) stay there; this is the consumer's row."""
    from orpheus.numerics.operator import adjointable

    sn, sig_t = _sn("slab_2g")
    inverse: Any = _lc(sn, sig_t).inverse()   # typed as the original family row (a list of tuples) typed it
    y = _sn_composite(sn)
    expected = adjointable(inverse.inner)
    if inverse.is_transpose_invertible != expected:
        pytest.fail(f"[sweep_slab] is_transpose_invertible={inverse.is_transpose_invertible}, but the forward's adjointable is {expected}")
    if not expected:
        return
    got, ref = flat(inverse.solve_transpose(y)), flat(inverse.inner.apply_transpose(y))
    if not np.array_equal(got, ref):
        pytest.fail(f"[sweep_slab] solve_transpose ≠ the forward's apply_transpose (max |Δ| {np.max(np.abs(got - ref)):.2e})")

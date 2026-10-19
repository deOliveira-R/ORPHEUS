r"""The adjointability flag is honest (ruling 27 of ``numerics_untangling.md``).

Two questions open the two transpose doors, and this module asserts each one
in both directions, over a registry of operators the tree's own fixtures build.

* ``LinearOperator.is_adjointable`` answers *a Euclidean transpose is coded*:
  it gates :meth:`~orpheus.numerics.operator.LinearOperator.dual`, the dual
  arrow :math:`A^{\mathsf T}`.
* :func:`~orpheus.numerics.operator.posed` answers *both ends' metrics give
  Riesz legs*: an unbound end, or a ``DualSpace`` end (the G² trap of ruling
  26), leaves the Hilbert adjoint unposed, unless it is metric-free.

:attr:`~orpheus.numerics.operator.LinearOperator.H` = ``♯ ∘ dual ∘ ♭`` needs both.
One flag answering both questions was the defect of the first fix: a rank-one
dyad with unbound ends read ``False`` and broke the fission adjoint built over
it (`[M]` 2026-10-18, qa: 7 of 7 fission builders, ``solve_sn_adjoint``).

Leg 1, each prediction holds (ruling 27).
    ``.dual()`` constructs and applies on a codomain sample exactly when
    ``is_adjointable``; ``.H`` constructs and applies exactly when
    ``is_adjointable`` and ``posed``; a refused door raises ``MissingAdjoint``
    at construction, never another exception. Where a door opens, the leg also
    checks its reciprocity, Euclidean for ``.dual()``
    (:math:`\langle A x, y\rangle = \langle x, A^{\mathsf T} y\rangle` on the flat
    coordinates) and in the operator's OWN end metrics for ``.H``
    (:math:`\langle A x, y\rangle_W = \langle x, A^{\mathsf H} y\rangle_V`), as a
    sanity leg. That leg is not the adjoint's correctness evidence: the
    structurally independent reciprocity gates this module ``rests_on`` are.

Leg 2, a realisable transpose is not hidden.
    Leg 1 is green on every base-``False`` operator, because ``.H`` then raises
    ``MissingAdjoint`` consistently. A flag reading ``False`` while the tree can
    already compute the transpose is the same lie in the other direction. For
    each such row the registry names a WITNESS ROUTE built from the algebra's
    own laws or from the assembled matrix, never from a second hand-written
    transpose: :math:`(A^{-1})^{\mathsf H} = (A^{\mathsf H})^{-1}`,
    :math:`(AB)^{\mathsf H} = B^{\mathsf H}A^{\mathsf H}`,
    :math:`(PB)^{\mathsf H} = B^{\mathsf H}P` for a coordinate projection ``P``,
    :math:`(M\oplus I)^{\mathsf H} = M^{\mathsf H}\oplus I`, the pointwise law
    :math:`C^{\mathsf H} = C`, and the sandwich :math:`\sharp_V M^{\mathsf T}\flat_W`
    of an assembled matrix. The witness is first checked by reciprocity
    against the forward (a broken witness is a fixture failure, reported as
    such); then the row must read ``is_adjointable is True``. A row whose
    parts have no coded transpose either (the 2-D sweep, the DSA correction)
    carries no witness: flag ``False`` is a legal fix there.

The reds at landing (2026-10-18) were the fix's done-when. Six rows remain
strict ``xfail`` against #614 (the diffusion operators whose composite ends
speak flat arrays, :data:`_CARRIER_BLOCKED_LEG1` and
:data:`_CARRIER_BLOCKED_LEG2`): strict, so a row turns red the day its carrier
is fixed. The witness route of each leg-2 row is in its parametrize id
(``<row>--<route>``) and in :data:`_WITNESS_ROUTES`.

``-O``-safe: every assertion is a ``pytest.fail`` call.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Callable

import numpy as np
import pytest
from scipy import sparse

from orpheus.geometry import BC, CoordSystem
from orpheus.numerics import frame_faces
from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.coupled_system import (
    CoupledField,
    CoupledOperator,
    CoupledSpace,
)
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.functional import InnerProductFunctional
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import (
    DiagonalOperator,
    IdentityOperator,
    LinearOperator,
    MissingAdjoint,
    MissingAssembly,
    TensorProductOperator,
    assemblable,
    outer,
    posed,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.operator import RieszLowerOperator
from orpheus.numerics.operator import RieszRaiseOperator

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

#: Sanity-leg reciprocity tolerance for a direct (non-iterative) operator:
#: the reductions are a few inner products of O(10²) terms.
_EXACT_RTOL = 1e-12
#: The Green splitting's convergence tolerance, and its reciprocity band
#: (10 × the tolerance, amplified by the conditioning of the het slab).
_GREEN_TOL = 1e-13
_GREEN_RTOL = 1e-9


# ═══════════════════════════════════════════════════════════════════════
# Carrier helpers
# ═══════════════════════════════════════════════════════════════════════


def _like(v: Any, rng: np.random.Generator) -> Any:
    r"""A random element of the same carrier (class, space, shape) as ``v``."""
    if isinstance(v, np.ndarray):
        return rng.standard_normal(v.shape)
    if isinstance(v, CoupledField):
        return replace(v, systems=tuple(_like(m, rng) for m in v.systems))
    if hasattr(v, "interior") and hasattr(v, "boundary"):
        return replace(
            v, interior=_like(v.interior, rng), boundary=_like(v.boundary, rng),
        )
    if hasattr(v, "values"):
        return replace(v, values=rng.standard_normal(np.shape(v.values)))
    raise TypeError(f"_like: no sampler for {type(v).__name__}")


def _flat(v: Any) -> np.ndarray:
    if isinstance(v, np.ndarray):
        return v.ravel()
    if hasattr(v, "to_flat"):
        return np.asarray(v.to_flat()).ravel()
    return np.asarray(v.values).ravel()


def _pairing(space: Any, a: Any, b: Any) -> float:
    r"""``⟨a, b⟩`` in ``space``'s own metric; Euclidean on an unbound end."""
    if space is None:
        return float(_flat(a) @ _flat(b))
    return float(space.inner_product(_leaf_values(a), _leaf_values(b)))


def _leaf_values(v: Any) -> Any:
    r"""A single-space field's values; a composite or an array passes through.

    A product space's ``inner_product`` reads bare arrays (``[M]`` 2026-10-18:
    the harmonic face's moment and angular ends refuse a typed leaf there),
    while a direct-sum space reads its composite.
    """
    if hasattr(v, "values") and not hasattr(v, "interior"):
        return v.values
    return v


def _rel(lhs: float, rhs: float) -> float:
    return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)


def _err(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {str(exc)[:160]}"


# ═══════════════════════════════════════════════════════════════════════
# The registry
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class _Posed:
    r"""One registry row, built: the operator, a domain sample, and how to act.

    ``pair_v`` / ``pair_w`` default to the operator's
    own end metrics (:func:`_pairing`). ``y`` is a codomain sample of the
    class the adjoint admits; it defaults to the forward's output class (an
    S\ :sub:`N` loss maps flux to source, and its adjoint takes the dual of
    flux, a flux-classed composite, so those rows pass one). ``witness`` maps a codomain sample to
    the adjoint action through a law of the algebra (leg 2), or is ``None``.
    """

    op: LinearOperator
    x: Any
    y: Any = None
    pair_v: Callable[[Any, Any], float] | None = None
    pair_w: Callable[[Any, Any], float] | None = None
    witness: Callable[[Any], Any] | None = None
    rtol: float = _EXACT_RTOL
    extra: dict = field(default_factory=dict)

    def codomain_sample(self, rng: np.random.Generator) -> Any:
        return _like(self.act(self.x) if self.y is None else self.y, rng)

    def act(self, x: Any) -> Any:
        return self.op.apply(x)

    def inner_v(self, a: Any, b: Any) -> float:
        if self.pair_v is not None:
            return self.pair_v(a, b)
        return _pairing(self.op.domain, a, b)

    def inner_w(self, a: Any, b: Any) -> float:
        if self.pair_w is not None:
            return self.pair_w(a, b)
        return _pairing(self.op.codomain, a, b)


def _sn(case: str):
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    return G._BUILDERS[case]()


def _sn_composite(sn, seed: int = 2026):
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    return G._random_composite(sn, np.random.default_rng(seed))


def _lc(sn, sig_t):
    from orpheus.sn.operators.streaming import StreamingOperator
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    return StreamingOperator.pose(sn) + MultiplicationOperator.from_mesh(
        sig_t, sn,
    )


def _source_like(sn, sig_t, seed: int):
    r"""A source-classed composite: the loss applied to a random flux."""
    return _lc(sn, sig_t).apply(_sn_composite(sn, seed))


# ── numerics: the Riesz legs, the dyad, the assembled matrix, the frame ──

_V = FunctionSpace("V_honesty", (3,), inner_product_weights=np.array([0.5, 2.0, 4.0]))
_W = FunctionSpace("W_honesty", (4,), inner_product_weights=np.array([3.0, 0.25, 1.5, 0.7]))


def _riesz_lower_primal() -> _Posed:
    return _Posed(RieszLowerOperator(_V), np.array([1.0, -2.0, 0.5]))


def _riesz_raise_primal() -> _Posed:
    # ♯ : V* → V, so its sample is a co-vector (a bare array on V*).
    return _Posed(RieszRaiseOperator(_V), np.array([0.3, 1.0, -1.5]))


def _rank_one() -> _Posed:
    D = outer(
        np.array([1.0, 2.0]),
        InnerProductFunctional(np.array([3.0, 4.0, 5.0]), axis=0),
    )
    return _Posed(D, np.array([0.2, -0.7, 1.1]))


def _sparse_assembled_bound() -> _Posed:
    M = np.array([[1.0, 0.0, 2.0], [0.0, 3.0, -1.0], [4.0, 0.5, 0.0],
                  [0.0, -2.0, 1.0]])
    op = SparseAssembledOperator(sparse.csr_array(M), domain=_V, codomain=_W)
    return _Posed(op, np.array([0.4, -1.0, 2.0]))


def _inverse_of_diagonal() -> _Posed:
    r"""``D⁻¹`` of an unbound diagonal: metric-free, so ``.H`` is posed with no
    spaces, and its transpose is the diagonal's transposed solve (qa F4: that
    solve had no catcher)."""
    return _Posed(DiagonalOperator(np.array([2.0, -0.5, 4.0])).inverse(), np.array([1.0, 3.0, -2.0]))


def _inverse_of_scaled_product() -> _Posed:
    r"""``(D @ 2D)⁻¹``: a product whose factor is a scaled operator, which has
    no native transposed solve; the product reaches it through each factor's
    canonical inverse (the elegance review: it read not adjointable while
    ``(2D)⁻¹`` read adjointable)."""
    from orpheus.numerics.operator import ScaledOperator

    D = DiagonalOperator(np.array([2.0, -0.5, 4.0]))
    return _Posed((D @ ScaledOperator(2.0, D)).inverse(), np.array([1.0, 3.0, -2.0]))


def _fission_dyad_tensor() -> _Posed:
    r"""``|χ⟩⟨νΣ_f| ⊗ I``, the factor shape of every fission operator, as the
    isotropic fission's own kernel: ends ``None``.

    The dyad's transpose is coded, so ``.dual()`` opens; the tensor is unbound
    and not a multiplier, so ``.H`` is honestly unposed. The defect of the
    first fix read the dyad's flag ``False`` here, and every fission adjoint
    built over it raised (`[M]` 2026-10-18, qa)."""
    from tests.gates.sn.operators import test_isotropic_fission as IF

    fission = IF._op()
    return _Posed(fission.kernel, np.ones(fission.domain.shape))


def _isotropic_fission() -> _Posed:
    from tests.gates.sn.operators import test_isotropic_fission as IF

    op = IF._op()
    return _Posed(op, np.ones(op.domain.shape))


def _frame_face(rule: Callable[[], Quadrature], L: int, face: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        frame = frame_faces.angular_frame(rule(), L)
        op = {
            "analysis": frame_faces.FrameAnalysis,
            "reconstruction": frame_faces.FrameReconstruction,
        }[face].of(frame)
        rng = np.random.default_rng(31)
        return _Posed(op, rng.standard_normal(op.domain.shape))

    return build


# ── diffusion: the leaves, the loss, the flat view, the resolvent ──


def _diffusion():
    from tests.gates.diffusion import test_operators as D

    mesh = D._diffusion_mesh(BC("reflective"), BC("vacuum"))
    return D, mesh


def _assembled_sandwich(op: LinearOperator, template: Any) -> Callable[[Any], Any]:
    r"""The witness ``♯_V ∘ Mᵀ ∘ ♭_W`` from ``op``'s assembled matrix ``M``.

    The transpose is the sparse matrix's own (``M.T``); the metrics are the
    operator's end spaces' (``apply_metric`` / ``apply_inverse_metric``). No
    operator-level transpose is written by hand.
    """
    if not assemblable(op):
        raise MissingAssembly(f"{type(op).__name__} does not assemble")
    M = op.assemble().matrix
    V, W = op.domain, op.codomain
    if V is None or W is None:
        raise RuntimeError(f"fixture drift: {type(op).__name__} has an unbound end")

    def witness(y: Any) -> Any:
        lowered = _flat(W.apply_metric(y))
        transposed = type(template).from_flat(np.asarray(M.T @ lowered), template)
        return V.apply_inverse_metric(transposed)

    return witness


def _diffusion_leaf(kind: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
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
        return _Posed(op, x, witness=_assembled_sandwich(op, x))

    return build


def _diffusion_flat_pairing(template: Any) -> Callable[[Any, Any], float]:
    space = template.space

    def pair(a: Any, b: Any) -> float:
        return float(space.inner_product(
            type(template).from_flat(np.asarray(a), template),
            type(template).from_flat(np.asarray(b), template),
        ))

    return pair


def _diffusion_flattened() -> _Posed:
    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    A = D._loss(mesh, mesh.mat_xs)
    op = FlattenedOperator(A, x)
    sandwich = _assembled_sandwich(A, x)
    pair = _diffusion_flat_pairing(x)
    return _Posed(
        op, x.to_flat(), pair_v=pair, pair_w=pair,
        witness=lambda y: _flat(sandwich(type(x).from_flat(np.asarray(y), x))),
    )


def _diffusion_resolvent() -> _Posed:
    r"""The production diffusion resolvent, ``diffusion/solver.py``'s spelling."""
    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    op = MatrixInverseOperator(FlattenedOperator(D._loss(mesh, mesh.mat_xs), x))
    pair = _diffusion_flat_pairing(x)
    return _Posed(op, x.to_flat(), pair_v=pair, pair_w=pair)


def _diffusion_assembled() -> _Posed:
    from orpheus.diffusion.operators import LeakageOperator

    D, mesh = _diffusion()
    x = D._random_flux(mesh)
    pair = _diffusion_flat_pairing(x)
    return _Posed(
        LeakageOperator(mesh).assemble(), x.to_flat(), pair_v=pair, pair_w=pair,
    )


# ── SN: the leaves and composites on a 1-D and a 2-D pose ──


def _sn_streaming(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        from orpheus.sn.operators.streaming import StreamingOperator

        sn, _ = _sn(case)
        x = _sn_composite(sn)
        return _Posed(StreamingOperator.pose(sn), x, y=x)

    return build


def _sn_streaming_collision(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        sn, sig_t = _sn(case)
        x = _sn_composite(sn)
        return _Posed(_lc(sn, sig_t), x, y=x)

    return build


def _sn_sweep(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        sn, sig_t = _sn(case)
        return _Posed(_lc(sn, sig_t).inverse(), _source_like(sn, sig_t, 7))

    return build


def _sn_boundary(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        from orpheus.sn.operators.boundary import SNBoundaryOperator

        sn, _ = _sn(case)
        x = _sn_composite(sn)
        return _Posed(SNBoundaryOperator(sn), x, y=x)

    return build


def _sn_within_group_loss(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

        sn, sig_t = _sn(case)
        x = _sn_composite(sn)
        return _Posed(G._loss_operator(sn, sig_t), x, y=x)

    return build


def _sn_full_loss(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

        sn, A, _ = G._FULL_LOSS_BUILDERS[case]()
        x = _sn_composite(sn)
        return _Posed(A, x, y=x)

    return build


def _sn_masked_boundary_lower() -> _Posed:
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
    return _Posed(lower, x, y=x, witness=lambda y: BH.apply(project(y)))


def _sn_coupled_substitution(case: str) -> Callable[[], _Posed]:
    r"""The 1×1 block grid ``[[L + C]]`` and its substitution inverse."""

    def build() -> _Posed:
        sn, sig_t = _sn(case)
        space = CoupledSpace.from_systems([sn.full_field_space])
        grid = CoupledOperator([[_lc(sn, sig_t)]], domain=space, codomain=space)
        x = CoupledField(systems=(_source_like(sn, sig_t, 5),))
        return _Posed(grid.inverse(), x)

    return build


def _sn_inverse_of_leaf() -> _Posed:
    r"""``InverseOperator(C)`` of the pointwise collision leaf.

    Witness ``(C⁻¹)ᴴ = (Cᴴ)⁻¹ = C⁻¹``: ``C`` is a real multiplier, metric-free
    (``Cᴴ = C`` in every diagonal metric), so its inverse is self-adjoint.
    """
    from orpheus.numerics.operator import InverseOperator
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    sn, sig_t = _sn("slab_2g")
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    if not C.is_metric_free_adjoint:
        raise RuntimeError("fixture drift: C is no longer metric-free")
    op = C.inverse()
    if type(op) is not InverseOperator:
        raise RuntimeError(f"fixture drift: C.inverse() is {type(op).__name__}")
    return _Posed(op, _source_like(sn, sig_t, 3), witness=op.apply)


def _sn_inverse_of_product() -> _Posed:
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
    return _Posed(
        op, _source_like(sn, sig_t, 4),
        witness=lambda y: sweep_h.apply(c_inv.apply(y)),
    )


def _sn_green() -> _Posed:
    r"""``GreenOperator((L + C) − S)`` on the het 2-G slab.

    Witness ``((L+C) − S)⁻ᴴ = ((L+C)ᴴ − Sᴴ)⁻¹``: the sum law, then the Green
    splitting of the adjoint sum (its leading term inverts by the swap law
    onto the 1-D reverse sweep). The witness is INDEPENDENT of the row's
    ``.H``: that is the base ``♯ ∘ G.dual() ∘ ♭``, the Green of the DUAL sum
    between the Riesz legs. Until the revision of 2026-10-18 the Green
    overrode ``.H`` with this very Green of the adjoint sum, and the witness
    then agreed with it to exactly 0.0, a tautology (the elegance review).
    """
    sn, sig_t, lc, S, op = _green_case()
    adjoint = GreenOperator(lc.H - S.H, max_iter=4000, tol=_GREEN_TOL)
    # Both ends are source → flux, so the adjoint's input is source-classed.
    # No row threads a seed: the Green's cold start is its first Neumann
    # partial sum, in the iterate's own class.
    return _Posed(
        op, _source_like(sn, sig_t, 6), y=_source_like(sn, sig_t, 8),
        witness=adjoint.apply, rtol=_GREEN_RTOL,
    )


def _green_case():
    r"""``GreenOperator((L + C) − S)`` on the het 2-G slab, with its pose."""
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    sn, _, S = G._full_loss_case(CoordSystem.CARTESIAN, 2)
    sig_t = np.asarray(sn.mat_xs.total_cross_section_field.values, dtype=float)
    lc = _lc(sn, sig_t)
    return sn, sig_t, lc, S, GreenOperator(lc - S, max_iter=4000, tol=_GREEN_TOL)


def _sn_green_composed() -> _Posed:
    r"""``C @ G``: a Green inside a product (qa, 2026-10-18: ``.H.apply`` on
    ``C@G``, ``G@C``, ``2·G`` and ``G+G`` raised while the Green read
    ``True``, because the Green coded no ``apply_transpose`` and its cold start
    was a source-classed zero)."""
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    sn, sig_t, _, _, G = _green_case()
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    return _Posed(C @ G, _source_like(sn, sig_t, 9), y=_source_like(sn, sig_t, 10), rtol=_GREEN_RTOL)


def _sn_bulk_analysis() -> _Posed:
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

    return _Posed(op, x, witness=witness)


def _sn_flux_analysis_face(case: str) -> Callable[[], _Posed]:
    r"""The scattering operator's minted flux-analysis face ``M`` (the face
    :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` lifts)."""

    def build() -> _Posed:
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
        return _Posed(face, x, y=np.asarray(face.apply(x).values))

    return build


def _sn_fission(case: str) -> Callable[[], _Posed]:
    def build() -> _Posed:
        from orpheus.transport.operators.fission import FissionOperator

        sn, _ = _sn(case)
        x = _sn_composite(sn)
        op = FissionOperator.from_solver_data(mat_xs=sn.mat_xs, space=sn.full_field_space)
        return _Posed(op, x, y=x)

    return build


def _sn_coupled_substitution_carrying() -> _Posed:
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
    return _Posed(
        M.inverse(), M.apply(random_state(record, seed=31)),
        y=M.apply(random_state(record, seed=37)),
    )


def _dual_pairing(primal: Any) -> Callable[[Any, Any], float]:
    r"""``⟨a, b⟩`` on the dual of ``primal`` in its honest metric :math:`G^{-1}`
    (a ``DualSpace`` carries :math:`G`, #617), so a dual arrow's ``.H`` is
    checked in the metric its law ``(Aᵀ)ᴴ = (Aᴴ)ᵀ`` holds in."""
    return lambda a, b: float(_flat(a) @ _flat(primal.apply_inverse_metric(b)))


def _sn_sweep_dual() -> _Posed:
    r"""``((L+C)⁻¹)ᵀ``, the dual arrow of the slab sweep: its ends are
    ``DualSpace``\ s, yet its ``.H`` is posed by the dagger–dual commutation
    (qa, 2026-10-18: 52 of 308 census rows were dual arrows whose ``.H``
    opened while the end-reading predicate said unposed)."""
    sn, sig_t = _sn("slab_2g")
    sweep = _lc(sn, sig_t).inverse()
    x = _sn_composite(sn)
    return _Posed(
        sweep.dual(), x, y=_source_like(sn, sig_t, 12),
        pair_v=_dual_pairing(sweep.codomain), pair_w=_dual_pairing(sweep.domain),
    )


def _riesz_lower_dual() -> _Posed:
    r"""``♭ᵀ``: the dual arrow of an unposed operator stays unposed."""
    return _Posed(RieszLowerOperator(_V).dual(), np.array([1.0, -2.0, 0.5]))


def _sn_green_cart2d() -> _Posed:
    r"""A Green over the 2-D head ``(L+C)``, which has no reverse sweep (#616):
    both doors refuse (qa: no row held a Green over a 2-D head, so dropping
    its head clause reddened nothing)."""
    from tests.gates.sn.operators import test_g_adjoint_reciprocity as G

    sn, _, S = G._full_loss_case_cart2d()
    sig_t = np.asarray(sn.mat_xs.total_cross_section_field.values, dtype=float)
    op = GreenOperator(_lc(sn, sig_t) - S, max_iter=4000, tol=_GREEN_TOL)
    return _Posed(op, _source_like(sn, sig_t, 13), rtol=_GREEN_RTOL)


def _sn_dsa_correction() -> _Posed:
    from orpheus.sn.acceleration import DSACorrection

    sn, _ = _sn("slab_2g")
    x = _sn_composite(sn)
    return _Posed(DSACorrection.from_problem(sn), x, y=x)


#: The registry: row id → builder. One row per operator the gate judges.
_REGISTRY: dict[str, Callable[[], _Posed]] = {
    # numerics
    "riesz_lower_primal": _riesz_lower_primal,
    "riesz_raise_primal": _riesz_raise_primal,
    "rank_one_unbound": _rank_one,
    "sparse_assembled_bound": _sparse_assembled_bound,
    "inverse_of_diagonal_unbound": _inverse_of_diagonal,
    "inverse_of_scaled_product_unbound": _inverse_of_scaled_product,
    "fission_dyad_tensor_unbound": _fission_dyad_tensor,
    "isotropic_fission_bound": _isotropic_fission,
    "frame_analysis_gl8_L2": _frame_face(lambda: Quadrature.gauss_legendre(8), 2, "analysis"),
    "frame_reconstruction_gl8_L2": _frame_face(lambda: Quadrature.gauss_legendre(8), 2, "reconstruction"),
    "frame_analysis_product48_L4": _frame_face(lambda: Quadrature.product(4, 8), 4, "analysis"),
    "frame_reconstruction_product48_L4": _frame_face(lambda: Quadrature.product(4, 8), 4, "reconstruction"),
    # diffusion
    "diffusion_leakage": _diffusion_leaf("leakage"),
    "diffusion_boundary": _diffusion_leaf("boundary"),
    "diffusion_loss": _diffusion_leaf("loss"),
    "diffusion_flattened_loss": _diffusion_flattened,
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
    "riesz_lower_primal_dual": _riesz_lower_dual,
    "sn_fission_slab": _sn_fission("slab_2g"),
    "sn_fission_cart2d": _sn_fission("cart2d_2g"),
    "sn_coupled_substitution_carrying_sphere": _sn_coupled_substitution_carrying,
    "sn_flux_analysis_face_slab": _sn_flux_analysis_face("slab_2g"),
    "sn_flux_analysis_face_cart2d": _sn_flux_analysis_face("cart2d_2g"),
    "sn_bulk_analysis_cart2d": _sn_bulk_analysis,
    "sn_dsa_correction_slab": _sn_dsa_correction,
}

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
    "sn_bulk_analysis_cart2d": "MH-direct-sum-I",
}

#: Rows whose witness route is written but cannot run today, because a PART's
#: own ``.H`` does not apply (constraint 3 of the dispatch: no coded route, so
#: report-only for leg 2). Each moves into :data:`_WITNESS_ROUTES` when the
#: named part's leg-1 row turns green.
_PENDING_WITNESS_ROUTES: dict[str, tuple[str, str]] = {}


# ═══════════════════════════════════════════════════════════════════════
# Leg 1 — the flag predicts the behaviour
# ═══════════════════════════════════════════════════════════════════════


#: Rows blocked by the diffusion carrier defect (#614): composite-ended
#: operators that speak flat arrays, so no adjoint of them can be realised
#: until the assembly homomorphism keeps their carriers typed. Strict, so a
#: row turns red the day its carrier is fixed and the mark cannot go stale.
#: Leg 1: the flag reads True while ``.H.apply`` raises on the flat carrier.
_CARRIER_BLOCKED_LEG1 = frozenset({"diffusion_assembled_leakage", "diffusion_resolvent"})
#: Leg 2: the flag reads False while the assembled-matrix witness passes.
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


def _open_door(door: str, op: LinearOperator, y: Any) -> tuple[Any, list[str]]:
    r"""Build ``op.dual()`` or ``op.H`` and apply it to ``y``; the image, or
    what refused. A construction refusal other than ``MissingAdjoint`` is
    itself a report."""
    try:
        T = op.dual() if door == "dual" else op.H
    except MissingAdjoint as exc:
        return None, [f".{door} refused: {_err(exc)}"]
    except Exception as exc:  # noqa: BLE001 — the reading IS the exception
        return None, [f".{door} raised a non-MissingAdjoint {_err(exc)}", "foreign"]
    try:
        return T.apply(y), []
    except Exception as exc:  # noqa: BLE001
        return None, [f".{door}.apply raised {_err(exc)}"]


@pytest.mark.parametrize("row", _rows(list(_REGISTRY), _CARRIER_BLOCKED_LEG1))
@pytest.mark.catches("ERR-117")
def test_each_transpose_door_opens_iff_predicted(row: str) -> None:
    r"""``.dual()`` works ⟺ ``is_adjointable``; ``.H`` works ⟺ ``is_adjointable ∧ posed``.

    A refused door raises ``MissingAdjoint`` at construction. Where a door
    opens, the sanity leg: Euclidean reciprocity for ``.dual()``, and
    ``⟨A x, y⟩_W = ⟨x, A.H y⟩_V`` in the operator's own end metrics for ``.H``.
    """
    built = _REGISTRY[row]()
    op = built.op
    rng = np.random.default_rng(1729)
    x = _like(built.x, rng)
    y = built.codomain_sample(rng)

    flag, is_posed = bool(op.is_adjointable), posed(op)
    Ax = built.act(x)
    for door, predicted in (("dual", flag), ("H", flag and is_posed)):
        image, report = _open_door(door, op, y)
        if "foreign" in report:
            pytest.fail(f"[{row}] {type(op).__name__}: {report[0]}")
        if predicted != (not report):
            detail = "; ".join(report) if report else f".{door} and its apply succeed"
            pytest.fail(
                f"[{row}] {type(op).__name__}.is_adjointable={flag}, "
                f"posed={is_posed}, so .{door} should "
                f"{'open' if predicted else 'refuse'}, but {detail}"
            )
        if not predicted:
            continue
        if door == "dual":
            lhs, rhs = float(_flat(Ax) @ _flat(y)), float(_flat(x) @ _flat(image))
            pairing = "⟨Ax,y⟩ vs ⟨x,Aᵀy⟩ (Euclidean)"
        else:
            lhs, rhs = built.inner_w(Ax, y), built.inner_v(x, image)
            pairing = "⟨Ax,y⟩_W vs ⟨x,A.H y⟩_V"
        if not _rel(lhs, rhs) < built.rtol:
            pytest.fail(
                f"[{row}] {type(op).__name__}: reciprocity {pairing}: "
                f"{lhs:.12e} vs {rhs:.12e} (rel {_rel(lhs, rhs):.2e} ≥ {built.rtol:.0e})"
            )


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
    x = _like(posed.x, rng)
    y = posed.codomain_sample(rng)
    lhs = posed.inner_w(posed.act(x), y)
    rhs = posed.inner_v(x, posed.witness(y))
    if not _rel(lhs, rhs) < posed.rtol:
        pytest.fail(
            f"[{row}] FIXTURE: the witness route {_WITNESS_ROUTES[row]!r} is not "
            f"the adjoint (rel {_rel(lhs, rhs):.2e}); repair the witness, not the flag"
        )
    if not posed.op.is_adjointable:
        pytest.fail(
            f"[{row}] {type(posed.op).__name__}.is_adjointable=False while the "
            f"route {_WITNESS_ROUTES[row]!r} realises its adjoint (reciprocity "
            f"rel {_rel(lhs, rhs):.1e})"
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


# ═══════════════════════════════════════════════════════════════════════
# The inverse family's verbs, by value (qa, 2026-10-18: four laws no
# mutation reddened, because leg 1 reads the same predicates the doors do)
# ═══════════════════════════════════════════════════════════════════════


class _SolveOnlyLeaf(LinearOperator):
    r"""An invertible leaf with no transpose: its inverse has no transposed solve."""

    def __init__(self, c: np.ndarray) -> None:
        self._c = np.asarray(c, dtype=float)

    @property
    def domain(self):
        return None

    @property
    def codomain(self):
        return None

    def apply(self, x, /):
        return self._c * x

    def solve(self, b, /):
        return b / self._c

    @property
    def is_invertible(self) -> bool:
        return True


def _inverse_family_rows():
    from orpheus.numerics.operator import InverseOperator

    sn, sig_t = _sn("slab_2g")
    D = DiagonalOperator(np.array([2.0, -0.5, 4.0]))
    return [
        ("inverse_of_diagonal", D.inverse(), np.array([1.0, 3.0, -2.0])),
        ("sweep_slab", _lc(sn, sig_t).inverse(), _sn_composite(sn)),
        ("inverse_of_solve_only_leaf", InverseOperator(_SolveOnlyLeaf(np.array([2.0, 3.0]))), np.ones(2)),
    ]


@pytest.mark.parametrize("index", range(3), ids=["inverse_of_diagonal", "sweep_slab", "inverse_of_solve_only_leaf"])
def test_the_inverse_swaps_the_forwards_transpose_verbs(index: int) -> None:
    r"""``(A⁻¹)ᵀ`` solves with ``Aᵀ``: the inverse's transposed solve IS the
    forward's transpose, and it exists iff the forward transposes."""
    from orpheus.numerics.operator import adjointable

    name, inverse, y = _inverse_family_rows()[index]
    expected = adjointable(inverse.inner)
    if inverse.is_transpose_invertible != expected:
        pytest.fail(
            f"[{name}] is_transpose_invertible={inverse.is_transpose_invertible}, "
            f"but the forward's adjointable is {expected}"
        )
    if not expected:
        return
    got, ref = _flat(inverse.solve_transpose(y)), _flat(inverse.inner.apply_transpose(y))
    if not np.array_equal(got, ref):
        pytest.fail(f"[{name}] solve_transpose ≠ the forward's apply_transpose (max |Δ| {np.max(np.abs(got - ref)):.2e})")


def test_the_inverse_of_a_multiplier_is_a_metric_free_adjoint() -> None:
    r"""``(D⁻¹)ᴴ = D⁻¹`` on unbound ends: the inverse inherits the forward's
    metric-free adjoint, so ``.H`` is posed with no spaces and equals the
    dense ``diag(1/d)`` exactly."""
    d = np.array([2.0, -0.5, 4.0])
    y = np.array([1.0, 3.0, -2.0])
    got = DiagonalOperator(d).inverse().H.apply(y)
    if not np.array_equal(got, y / d):
        pytest.fail(f"(D⁻¹)ᴴ y = {got}, expected diag(1/d) y = {y / d}")


def test_a_converged_first_partial_sum_is_returned_unswept() -> None:
    r"""A cold Green start is the first Neumann partial sum, checked before it
    is driven: when it already meets the promise, one preconditioner solve is
    the whole apply, as it was when the driver started from zero (`[M]`
    2026-10-18: 1 solve old, 2 without the check)."""
    from orpheus.numerics.operator import PermutationOperator, ScaledOperator

    calls = [0]

    class _CountingDiagonal(DiagonalOperator):
        def solve(self, b_vec):
            calls[0] += 1
            return super().solve(b_vec)

    D = _CountingDiagonal(np.array([2.0, 3.0, 4.0, 5.0]))
    gain = ScaledOperator(3e-12, PermutationOperator(np.array([1, 2, 3, 0])))
    green = GreenOperator(D - gain, max_iter=10, tol=1e-10)
    green.apply(np.array([1.0, -2.0, 0.5, 3.0]))
    if calls[0] != 1:
        pytest.fail(f"a converged first partial sum cost {calls[0]} preconditioner solves, expected 1")

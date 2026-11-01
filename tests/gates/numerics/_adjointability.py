r"""The adjointability law (ruling 27 of ``numerics_untangling.md``), shared by its two test modules.

Two questions open the two transpose doors, and this module asserts each one
in both directions, over a registry of operators the tree's own fixtures build.

* ``LinearOperator.is_adjointable`` answers *a Euclidean transpose is coded*:
  it gates :meth:`~orpheus.numerics.operator.LinearOperator.dual`, the dual
  arrow :math:`A^{\mathsf T}`.
* :func:`~orpheus.numerics.operator.posed` answers *both ends' metrics give
  Riesz legs*: an unbound end leaves the Hilbert adjoint unposed, unless it
  is metric-free. (Until S6a a ``DualSpace`` end did too, the G² trap of
  ruling 26; a dual space now carries G⁻¹ and its legs exist.)

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
    structurally independent reciprocity gates each test module ``rests_on`` are.

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

Who runs what. ``tests/gates/numerics/test_adjointability_is_honest.py`` runs leg 1 on the numerics operators
(:data:`NUMERICS_REGISTRY`, built here) and the inverse family's numerics rows; the census
``tests/gates/test_adjointability_is_honest.py`` builds the consumers' rows (fission, diffusion, SN), runs leg 1 on
them, and leg 2 and the registry's reach on the union (layer DAG Phase 2, ruling 28: a census of one layer's
contract over every consumer lives at the ``tests/gates/`` root). Not collected (the leading underscore).

``-O``-safe: every assertion is a ``pytest.fail`` call.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Callable

import numpy as np
import pytest
from scipy import sparse

from orpheus.numerics import frame_faces
from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.functional import InnerProductFunctional
from orpheus.numerics.elementary_operators import DiagonalOperator, outer
from orpheus.numerics.operator import LinearOperator, MissingAdjoint, posed
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.operator import RieszLowerOperator
from orpheus.numerics.operator import RieszRaiseOperator


#: Sanity-leg reciprocity tolerance for a direct (non-iterative) operator:
#: the reductions are a few inner products of O(10²) terms.
_EXACT_RTOL = 1e-12


def like(v: Any, rng: np.random.Generator) -> Any:
    r"""A random element of the same carrier (class, space, shape) as ``v``."""
    if isinstance(v, np.ndarray):
        return rng.standard_normal(v.shape)
    if isinstance(v, CoupledField):
        return replace(v, systems=tuple(like(m, rng) for m in v.systems))
    if hasattr(v, "interior") and hasattr(v, "boundary"):
        return replace(
            v, interior=like(v.interior, rng), boundary=like(v.boundary, rng),
        )
    if hasattr(v, "values"):
        return replace(v, values=rng.standard_normal(np.shape(v.values)))
    raise TypeError(f"like: no sampler for {type(v).__name__}")


def flat(v: Any) -> np.ndarray:
    if isinstance(v, np.ndarray):
        return v.ravel()
    if hasattr(v, "to_flat"):
        return np.asarray(v.to_flat()).ravel()
    return np.asarray(v.values).ravel()


def _pairing(space: Any, a: Any, b: Any) -> float:
    r"""``⟨a, b⟩`` in ``space``'s own metric; Euclidean on an unbound end."""
    if space is None:
        return float(flat(a) @ flat(b))
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


def rel(lhs: float, rhs: float) -> float:
    return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)


def _err(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {str(exc)[:160]}"


@dataclass(frozen=True)
class Posed:
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
        return like(self.act(self.x) if self.y is None else self.y, rng)

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


_V = FunctionSpace("V_honesty", (3,), inner_product_weights=np.array([0.5, 2.0, 4.0]))


_W = FunctionSpace("W_honesty", (4,), inner_product_weights=np.array([3.0, 0.25, 1.5, 0.7]))


def _riesz_lower_primal() -> Posed:
    return Posed(RieszLowerOperator(_V), np.array([1.0, -2.0, 0.5]))


def _riesz_raise_primal() -> Posed:
    # ♯ : V* → V, so its sample is a co-vector (a bare array on V*).
    return Posed(RieszRaiseOperator(_V), np.array([0.3, 1.0, -1.5]))


def _rank_one() -> Posed:
    D = outer(
        np.array([1.0, 2.0]),
        InnerProductFunctional(np.array([3.0, 4.0, 5.0]), axis=0),
    )
    return Posed(D, np.array([0.2, -0.7, 1.1]))


def _sparse_assembled_bound() -> Posed:
    M = np.array([[1.0, 0.0, 2.0], [0.0, 3.0, -1.0], [4.0, 0.5, 0.0],
                  [0.0, -2.0, 1.0]])
    op = SparseAssembledOperator(sparse.csr_array(M), domain=_V, codomain=_W)
    return Posed(op, np.array([0.4, -1.0, 2.0]))


def _inverse_of_diagonal() -> Posed:
    r"""``D⁻¹`` of an unbound diagonal: its transpose ``.dual()`` is the
    diagonal's transposed solve (qa F4: that solve had no catcher), and its
    ``.H`` refuses as unbound, because a field multiplier is not metric-free
    (S6a review, #625: it commutes only with diagonal metrics). Until the
    review this row's ``.H`` opened."""
    return Posed(DiagonalOperator(np.array([2.0, -0.5, 4.0]), shape=(3,), spans=(0,)).inverse(), np.array([1.0, 3.0, -2.0]))


def _inverse_of_scaled_product() -> Posed:
    r"""``(D @ 2D)⁻¹``: a product whose factor is a scaled operator, which has
    no native transposed solve; the product reaches it through each factor's
    canonical inverse (the elegance review: it read not adjointable while
    ``(2D)⁻¹`` read adjointable). Its ``.H`` refuses as unbound since the S6a
    review (a word of field multipliers is not metric-free, #625)."""
    D = DiagonalOperator(np.array([2.0, -0.5, 4.0]), shape=(3,), spans=(0,))
    # S2 (D7): the scalar leaves the word, so ``D @ (2D)`` is the one-term
    # combination 2·(D D); its inverse is ½·(D D)⁻¹, the row's subject.
    return Posed((D @ (2.0 * D)).inverse(), np.array([1.0, 3.0, -2.0]))


def _frame_face(rule: Callable[[], Quadrature], L: int, face: str) -> Callable[[], Posed]:
    def build() -> Posed:
        frame = frame_faces.angular_frame(rule(), L)
        op = {
            "analysis": frame_faces.FrameAnalysis,
            "reconstruction": frame_faces.FrameReconstruction,
        }[face].of(frame)
        rng = np.random.default_rng(31)
        return Posed(op, rng.standard_normal(op.domain.shape))

    return build


def _riesz_lower_dual() -> Posed:
    r"""``♭ᵀ``: the dual arrow of the lowering leg; posed since S6a (its
    ``.H`` opens, test_dual_metric H1)."""
    return Posed(RieszLowerOperator(_V).dual(), np.array([1.0, -2.0, 0.5]))


#: The registry: row id → builder. One row per operator the gate judges.
#: The numerics operators' rows: the honesty law's witnesses in its own layer. The consumers' rows (fission,
#: diffusion, SN) are the root census's, tests/gates/test_adjointability_is_honest.py (layer DAG Phase 2 refresh,
#: ruling 28).
NUMERICS_REGISTRY: dict[str, Callable[[], Posed]] = {
    "riesz_lower_primal": _riesz_lower_primal,
    "riesz_raise_primal": _riesz_raise_primal,
    "rank_one_unbound": _rank_one,
    "sparse_assembled_bound": _sparse_assembled_bound,
    "inverse_of_diagonal_unbound": _inverse_of_diagonal,
    "inverse_of_scaled_product_unbound": _inverse_of_scaled_product,
    "frame_analysis_gl8_L2": _frame_face(lambda: Quadrature.gauss_legendre(8), 2, "analysis"),
    "frame_reconstruction_gl8_L2": _frame_face(lambda: Quadrature.gauss_legendre(8), 2, "reconstruction"),
    "frame_analysis_product48_L4": _frame_face(lambda: Quadrature.product(4, 8), 4, "analysis"),
    "frame_reconstruction_product48_L4": _frame_face(lambda: Quadrature.product(4, 8), 4, "reconstruction"),
    "riesz_lower_primal_dual": _riesz_lower_dual,
}


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


def assert_doors(row: str, built: "Posed") -> None:
    r"""Leg 1 on one built row (shared by this file's numerics rows and the root census's consumer rows)."""
    op = built.op
    rng = np.random.default_rng(1729)
    x = like(built.x, rng)
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
            lhs, rhs = float(flat(Ax) @ flat(y)), float(flat(x) @ flat(image))
            pairing = "⟨Ax,y⟩ vs ⟨x,Aᵀy⟩ (Euclidean)"
        else:
            lhs, rhs = built.inner_w(Ax, y), built.inner_v(x, image)
            pairing = "⟨Ax,y⟩_W vs ⟨x,A.H y⟩_V"
        if not rel(lhs, rhs) < built.rtol:
            pytest.fail(
                f"[{row}] {type(op).__name__}: reciprocity {pairing}: "
                f"{lhs:.12e} vs {rhs:.12e} (rel {rel(lhs, rhs):.2e} ≥ {built.rtol:.0e})"
            )


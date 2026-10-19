r"""The Green inverse operator — the preconditioned-splitting member of the
#226 inverse family (taxonomy §12 step 4).

A general operator SUM has no direct inverse: :math:`(A + B)^{-1}` is not a
function of :math:`A^{-1}` and :math:`B^{-1}` (Sherman–Morrison–Woodbury
applies only under low-rank structure — which the bulk blocks L, C lack
and the boundary block B HAS: rank ≤ N/2 per face, Issue #300), and —
unlike the schedule-triangular
``(L+C)`` — there is no substitution order that solves it in one pass.  What
it DOES have, when its LEADING term :math:`A` is invertible, is a convergent
SPLITTING: write the sum :math:`A - B` (gains carried with their signs), and

.. math::

    (A - B)^{-1}\,q
    \;=\; \bigl(I - A^{-1}B\bigr)^{-1} A^{-1} q
    \;=\; \sum_{k=0}^{\infty} \bigl(A^{-1}B\bigr)^{k}\,A^{-1}\,q ,

the NEUMANN SERIES of the :math:`A`-preconditioned splitting — in transport,
the multiple-scattering expansion: :math:`A^{-1}q` is the uncollided flux,
:math:`(A^{-1}S)^k A^{-1} q` the k-times-rescattered contribution (Lewis &
Miller §3.2; Adams & Larsen 2002 §II).  Its partial sums are EXACTLY the
Richardson / source-iteration iterates
:math:`x_{n+1} = A^{-1}(q + B\,x_n)` started from zero, convergent iff
:math:`\rho(A^{-1}B) < 1` — for the within-group transport loss
:math:`A_{\rm loss} = (L+C) - S` this is the scattering ratio bound
:math:`\rho \le \max_{\rm cell} c = \Sigma_s/\Sigma_t < 1`, guaranteed for
any physical (absorbing) medium.

**The name is a promise (taxonomy §13).**  :class:`GreenOperator` is the
discrete Green's function of its forward sum — ``G.apply(δ_j)`` is column j
of :math:`(A-B)^{-1}`, the flux response to a unit point source — and the
name is EARNED by the distinguishing G-Neumann invariant: ``G.apply(q)``
equals the converged Neumann partial sums, which a generic ``A⁻¹`` has no
splitting to satisfy.  The family is keyed by the mathematical OBJECT, not
the algorithm: a Richardson-realized Green and a Krylov-realized Green are
the SAME Green operator (which is why ``KrylovInverseOperator`` was rejected
as a name — "Krylov" names the orthogonal realization axis).

**Green WRAPS the driver; the drivers stay public primitives (taxonomy
§11.2, Pattern 5).**  The application engine is
:class:`~orpheus.numerics.iteration.SourceIteration` — the SAME Richardson
driver the SN solver and :class:`~orpheus.numerics.iteration.KEigenvalue`
consume directly.  This class derives the splitting from the sum's structure
(leading term → preconditioner via its own structure-keyed ``.inverse()`` —
the WDD :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` for
``(L+C)``, an :class:`~orpheus.numerics.operator.InverseOperator` for a
value leaf — remaining terms → negated gains) and hands it to the driver;
it re-implements NOTHING.  A GMRES-realized Green
(:class:`~orpheus.numerics.iteration.KrylovAcceleration` as the engine) is a
future realization STRATEGY of this same object, not a sibling type.

Note the two contracts are deliberately different: the k-eigenvalue path
keeps consuming the driver DIRECTLY (a warm-started, budget-bounded inner
RELAXATION — partial convergence is by design in nested iteration), while
``GreenOperator`` promises the CONVERGED :math:`A^{-1}q` or a loud
:class:`ConvergenceFailure`.  ``K = A_loss.inverse() @ F`` is the normal
FORM of the k-problem; production power iteration replaces the exact
:math:`G` with the inexact warm-started inner — classic inexact power
iteration, spelled at the driver layer.

**The ordering ruling (taxonomy §12 step 4, MUST-land).**  Operand spelling
selects the ALGORITHM:

* ``L + C`` never reaches this module — the fusion dispatch on
  :meth:`~orpheus.sn.operators.streaming.StreamingOperator.__add__` (#261:
  "the canonical and only ordering") returns the
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`
  specialisation, whose ``.inverse()`` override (→ direct sweep) shadows
  the generic :meth:`~orpheus.numerics.operator.OperatorSum.inverse` by
  MRO (type-as-structure, §11.1).
* ``A - S`` (invertible term FIRST) derives the physical splitting:
  preconditioner :math:`A^{-1}` (the sweep), gain :math:`S` — convergence
  certified by :math:`c < 1`.
* ``(-S) + A`` (invertible term NOT first) is REFUSED at construction
  (:class:`~orpheus.numerics.operator.NotInvertible`): the left-spine
  head is the DESIGNATED preconditioner, and the canonical ordering puts
  the invertible operator first.
* ``C + L`` (a legal spelling whose leading term happens to be invertible)
  constructs — the algebra cannot know :math:`\rho(C^{-1}L) > 1` — but the
  divergent splitting raises :class:`ConvergenceFailure` LOUDLY at apply
  time.  Same math as ``L + C``, different algorithm by spelling: the
  spelling chose collision-preconditioned Richardson instead of the sweep.
  Never a silent wrong answer (Cardinal Rule 1).

**Seed semantics (§38 application context; taxonomy §11.3).**
``apply(q, initial_guess=x₀)`` threads the seed as the splitting
iteration's START — a genuine consumer (contrast
:class:`~orpheus.numerics.operator.InverseOperator`'s accept-and-ignore).
The RESULT is :math:`A^{-1}q` to the configured tolerance regardless of the
seed (the convergence test, not the start, defines the promise), so the
kwarg is context, not identity.

**Module placement (taxonomy §17 W3).**  This is a LEAF module: it imports
both :mod:`~orpheus.numerics.operator` (the algebra) and
:mod:`~orpheus.numerics.iteration` (the drivers), and
:meth:`OperatorSum.inverse` LATE-imports it — mirroring the
``streaming.py`` → ``sweep_operator.py`` precedent — so the one-way
``iteration → operator`` import direction is preserved.

**#284 scope note.**  The Richardson driver feeds the preconditioner's
sweep only source-shaped right-hand sides (:math:`q + B\,x_n` writes
bulk/inflow rows), so the sweep-inverse's source-subspace contract holds on
this path. Feeding sign-indefinite RESIDUALS into a sweep inverse is the
GMRES-preconditioning case: since #200 (2026-10-04) the SN Krylov solve
preconditions with the sweep on the full typed state, and its gates
(``tests/gates/sn/solve/test_krylov_sweep_preconditioner.py``: linearity on
full and boundary-only residuals, and the boundary round trip
:math:`A M q_b = q_b`) witness that the sweep is a linear map on residuals,
trace rows included.
"""
from __future__ import annotations

import math
import operator
from functools import reduce
from typing import Any, cast

from orpheus.numerics.iteration import SourceIteration, _l2_norm, seeded_inverse
from orpheus.numerics.operator import (
    InverseWrapMixin,
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    OperatorSum,
    ScaledOperator,
    adjointable,
)

__all__ = [
    "ConvergenceFailure",
    "GreenOperator",
]


class ConvergenceFailure(RuntimeError):
    r"""An iterative inverse application failed to reach its tolerance.

    Raised by :meth:`GreenOperator.apply` when the splitting iteration
    exhausts its TOTAL budget with the promise unmet — the TRUE relative
    residual :math:`\lVert(A-B)\psi - q\rVert/\lVert q\rVert` still at or
    above ``tol`` (spec §18.A: the promise reads the equation residual,
    never the driver's ρ-blind increment) — or when that residual is
    non-finite (a divergent split overflows the iterate, or, sharper,
    overflows the driver's stopping DENOMINATOR to ``inf`` so it
    "converges" at ``increment = finite/inf = 0.0`` onto a garbage
    iterate).  In both shapes the iterate is NOT the promised
    :math:`A^{-1}q`, and returning it silently would be a wrong answer
    (Cardinal Rule 1).  A divergent splitting usually means the sum's
    SPELLING selected the wrong preconditioner (see the module
    docstring's ordering ruling); a slowly-converging one means the
    budget (``max_iter``) is too small for the physics (:math:`c \to 1`).

    Deliberately NOT raised by the driver itself:
    :class:`~orpheus.numerics.iteration.SourceIteration` returns its final
    iterate with the residual history because nested-iteration consumers
    (the k-eigenvalue inner) legitimately run partial, budget-bounded
    solves.  The PROMISE of an exact inverse — and therefore the raise —
    belongs to the inverse OBJECT.
    """


def _left_spine_terms(op: OperatorSum) -> list[LinearOperator]:
    r"""Flatten a plain-sum tree's LEFT SPINE into its ordered terms.

    ``(A - S) - F`` (built by left-associative ``+``/``-``) is
    ``OperatorSum(OperatorSum(A, -S), -F)``; the walk returns
    ``[A, -S, -F]`` so the splitting derivation sees the terms as
    spelled.  The walk descends EXACT ``OperatorSum`` nodes only:
    subclasses — the fused
    :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` family —
    are structural LEAVES (their sum-ness is an MRO fact, their identity
    a fused operator with its own direct inverse), and a RIGHT-nested
    sum ``A + (B + C)`` stays one composite gain term (a gain only needs
    ``apply``; flattening it would buy nothing).
    """
    terms: list[LinearOperator] = []
    node: LinearOperator = op
    while type(node) is OperatorSum:
        terms.append(node.b)
        node = node.a
    terms.append(node)
    terms.reverse()
    return terms


def _negated(term: LinearOperator) -> LinearOperator:
    r"""Return :math:`-{\rm term}`, unwrapping an existing ``(-1)·X`` scaling.

    The ``A - S`` spelling reaches the sum as ``ScaledOperator(-1, S)``
    (built by ``__sub__``); the splitting's gain is :math:`S` itself —
    unwrapping hands the driver the NAMED operator (its ``repr`` reads as
    the math, and the gain applies without a pass-through wrapper).  A
    genuinely additive term gets the explicit :math:`(-1)` scaling.
    """
    if isinstance(term, ScaledOperator) and term.scalar == -1.0:
        return term.op
    return ScaledOperator(-1.0, term)


class GreenOperator(InverseWrapMixin[OperatorSum], LinearOperator):
    r"""The inverse operator :math:`(A - B)^{-1}` of a general operator sum,
    realized as the convergent :math:`A`-preconditioned splitting.

    Built by :meth:`~orpheus.numerics.operator.OperatorSum.inverse` (the
    structure-keyed factory, taxonomy §3) or directly when the iteration
    budget needs configuring — ``.inverse()`` takes no arguments, so
    consumers needing a non-default ``max_iter``/``tol`` construct
    ``GreenOperator(sum_op, max_iter=..., tol=...)`` themselves.

    The splitting is derived from the sum's structure at construction
    (never re-derived per apply): the LEFT-SPINE HEAD is the designated
    preconditioner — refused loudly if not invertible (the canonical
    ordering spells the invertible operator first) — and every remaining
    term rides as a negated gain of the wrapped
    :class:`~orpheus.numerics.iteration.SourceIteration` driver.

    :meth:`apply` returns the CONVERGED :math:`A^{-1}q` or raises
    :class:`ConvergenceFailure` — never a silent partial iterate.  The
    wrap-delegate back-half (domain↔codomain swap /
    ``solve`` = the forward sum's matvec / object-identity involution
    ``inverse() → inner``) is inherited from
    :class:`~orpheus.numerics.operator.InverseWrapMixin` — this class is
    the THIRD sibling, the one that fired the twins' extraction trigger.

    The Euclidean transpose is the Green of the dual sum,
    :math:`((A - B)^{-1})^{\mathsf T} = (A^{\mathsf T} - B^{\mathsf T})^{-1}`
    (:meth:`dual`, which :meth:`apply_transpose` applies): the terms' dual
    arrows, head first, wrapped in a ``GreenOperator`` with the same budget.
    Its preconditioner is the dual arrow's swap law
    :math:`(A^{\mathsf T})^{-1} = (A^{-1})^{\mathsf T}`, so it exists exactly
    where the head's inverse transposes, i.e. where the head is
    transpose-invertible (the SN sweep: the 1-D scan family; the multi-D
    reverse scan is #616). :attr:`is_adjointable` reads that condition and
    that every term transposes. These override the mixin's inverse law,
    since a sum has no transposed solve of its own. The Hilbert adjoint is
    the base's ``♯ ∘ dual ∘ ♭``: the law
    :math:`((A - B)^{-1})^{*} = (A^{*} - B^{*})^{-1}` holds and the honesty
    gate checks it independently (its witness is the Green of the adjoint
    sum), but it is not a second route to ``.H``.

    Parameters
    ----------
    inner : OperatorSum
        The forward sum this is the inverse of.  Its left-spine head must
        be invertible (`is_invertible`) — it becomes the splitting's
        preconditioner through its OWN structure-keyed ``.inverse()``.
    max_iter : int, optional
        TOTAL splitting-step budget across :meth:`apply`'s refinement
        loop (default ``1000``).  Exhausting it with the promise unmet
        raises :class:`ConvergenceFailure`.
    tol : float, optional
        The promise tolerance (default ``1e-8``): :meth:`apply` returns
        only when the TRUE relative residual
        :math:`\lVert(A-B)\psi - q\rVert/\lVert q\rVert` is below it
        (spec §18.A — never the driver's ρ-blind increment).  ``tol=0.0``
        is unsatisfiable for an iterative inverse and always raises.

    Raises
    ------
    NotInvertible
        At construction, when the left-spine head is not invertible (the
        canonical-ordering refusal) — or from the head's own ``.inverse()``
        builder (e.g. a singular diagonal), where the invertibility
        obligation lives.
    """

    def __init__(
        self,
        inner: OperatorSum,
        *,
        max_iter: int = 1000,
        tol: float = 1e-8,
    ) -> None:
        terms = _left_spine_terms(inner)
        #: The signed left-spine terms, head first: the splitting's operands,
        #: and the terms whose dual arrows :meth:`dual` sums.
        self._terms = terms
        #: The Green this one is the dual of, when :meth:`dual` built it:
        #: the involution (Gᵀ)ᵀ = G as an object identity.
        self._dual_of: "GreenOperator | None" = None
        leading = terms[0]
        if not leading.is_invertible:
            raise NotInvertible(
                f"GreenOperator requires an invertible LEADING term: the "
                f"left-spine head of the sum is the splitting's "
                f"preconditioner (canonical ordering — spell the invertible "
                f"operator FIRST, as in A - S); got "
                f"{type(leading).__name__} with is_invertible=False."
            )
        super().__init__(inner)
        #: The promise tolerance: :meth:`apply` returns only when the TRUE
        #: relative residual is below it.  Doubles as the driver's
        #: per-call increment stopping heuristic.
        self.tol = float(tol)
        #: TOTAL splitting-step budget across the refinement loop.
        self.max_iter = int(max_iter)
        # The splitting, derived ONCE from structure: x_{n+1} =
        # A⁻¹(q + Σ gains·x_n) with gains = the negated non-leading terms
        # — handed to the SAME Richardson driver the solvers consume
        # (§11.2: wrap, never re-implement).
        self._preconditioner = seeded_inverse(leading)
        self._gains = tuple(_negated(t) for t in terms[1:])

    def _driver(self, budget: int) -> SourceIteration[Any]:
        r"""The Richardson driver of this splitting, allowed ``budget`` passes.

        Built per refinement call so the TOTAL budget is exact: each call
        receives what the earlier calls (and a cold start's first partial
        sum) left. Construction stores references only.
        """
        return SourceIteration(self._preconditioner, *self._gains, max_iter=budget, tol=self.tol)

    def _true_relative_residual(self, psi: Any, q: Any) -> float:
        r"""The PROMISE metric :math:`\lVert(A-B)\,\psi - q\rVert / \lVert q\rVert`.

        The driver's stopping test reads the iterate INCREMENT
        :math:`\lVert\Delta\psi\rVert/\lVert\psi\rVert`, which relates to
        this true residual by the factor :math:`\rho/(1-\rho)`
        (numerical-bug-signatures Signature 9, ρ-blind stopping;
        verification spec §18.A) — near :math:`\rho \to 1` an
        increment-converged iterate can sit orders of magnitude off the
        equation.  The inverse's promise is therefore checked (and, in
        :meth:`apply`, DRIVEN) on the equation residual itself, at one
        extra forward matvec per check.
        """
        defect = self.inner.apply(psi) - q
        return _l2_norm(defect) / max(_l2_norm(q), 1e-30)

    def apply(self, q: Any, /, *, initial_guess: Any | None = None) -> Any:
        r"""Return :math:`(A - B)^{-1}\,q` with TRUE relative residual
        below ``tol``, or raise :class:`ConvergenceFailure`.

        ``initial_guess`` seeds the splitting iteration's start (§38
        application context — the converged result is seed-independent to
        ``tol`` by the promise check itself); ``None`` starts the Neumann
        series from zero, whose partial sums the G-Neumann gate pins. The
        cold start is taken as the series' first partial sum
        :math:`A^{-1}q`, which is one splitting step from zero
        (:math:`A^{-1}(q + B\,0)`) and is counted as one: the zero itself is
        never formed, because it lives in the iterate's class (a flux) while
        the only template at hand is the right-hand side's (a source), and a
        gain refuses a source-classed zero (`[M]` 2026-10-18: ``TypeError``
        from the scattering gain on every typed S\ :sub:`N` Green applied
        without a seed, and on every composite that holds one). The defect
        is the driver's (#618), and this seed retires when the driver's own
        cold start stops forming the zero.

        The promise is verified — and driven — on the equation residual,
        not the driver's increment (spec §18.A): the driver stops on
        increment < tol, which under-delivers by :math:`\rho/(1-\rho)`,
        so each increment-converged iterate is CHECKED against the true
        residual and, if above tolerance with budget remaining, handed
        back to the driver as its own warm start (the refinement loop —
        a check-only design would falsely raise for every split with
        :math:`\rho > 1/2`, i.e. most physical scattering ratios).  The
        iteration math stays entirely inside the driver; this loop is
        pure tolerance bookkeeping.
        """
        # A cold start's first partial sum is checked before it is driven,
        # as the driver checks its own first iterate; a handed-in seed is
        # driven first.
        cold = initial_guess is None
        psi: Any = self._preconditioner.apply(q) if cold else initial_guess
        steps = int(cold)
        record = None
        while True:
            if not cold:
                # Only what the earlier passes left of the TOTAL budget: the
                # promise test below raises before the budget can be spent.
                psi, record = self._driver(self.max_iter - steps).solve(
                    q, initial_guess=psi,
                )
                # The driver's OWN pass count, not its trajectory length —
                # the two differ by one for SourceIteration (#340 F10), and
                # this loop's budget is spent in passes.  ``max(_, 1)``
                # survives as the guarantee that a driver returning nothing
                # still advances the outer counter, so the while-loop cannot
                # spin forever.
                steps += max(record.n_iterations, 1)
            cold = False
            true_res = self._true_relative_residual(psi, q)
            # NaN-safe promise test (spelled "converged", never "above
            # tol"): a hard-divergent split overflows the iterate — or,
            # sharper, overflows the driver's stopping DENOMINATOR ‖ψ‖ to
            # inf one step before the numerator, so the driver "converges"
            # at increment = finite/inf = 0.0 onto a ~1e154 garbage
            # iterate (found by the divergence tooth, 2026-07-02).  The
            # true residual of that iterate is huge/non-finite, so the
            # promise test catches both shapes.
            if true_res < self.tol:
                return psi
            if steps >= self.max_iter or not math.isfinite(true_res):
                criterion = None if record is None else record.binding_criterion
                increment = (
                    float("nan")
                    if criterion is None or criterion.last is None
                    else criterion.last
                )
                raise ConvergenceFailure(
                    f"GreenOperator: splitting iteration did not reach the "
                    f"promise — TRUE relative residual ‖(A-B)ψ − q‖/‖q‖ = "
                    f"{true_res:.3e} after {steps} step(s), tol="
                    f"{self.tol:.3e} (driver increment {increment:.3e}); "
                    f"the iterate is NOT (A-B)⁻¹q. A divergent split means "
                    f"the sum's SPELLING chose the wrong preconditioner "
                    f"(the leading term; see the canonical-ordering "
                    f"ruling) — a slow one, that max_iter is too small for "
                    f"the physics."
                )

    @property
    def is_adjointable(self) -> bool:
        # The dual sum needs every term's transpose, and its head
        # (A^T)^{-1} = (A^{-1})^T needs the dual arrow's swap law.
        return all(adjointable(t) for t in self._terms) and (
            self._terms[0].dual().is_invertible
        )

    def dual(self) -> "GreenOperator":
        r"""Return :math:`((A - B)^{-1})^{\mathsf T} = (A^{\mathsf T} - B^{\mathsf T})^{-1}`, the Green of the dual sum.

        The terms' dual arrows, head first, with the same budget; the
        preconditioner is the dual arrow's swap law
        :math:`(A^{\mathsf T})^{-1} = (A^{-1})^{\mathsf T}`. The Hilbert
        adjoint is NOT overridden: ``.H`` is the base's
        ``♯ ∘ dual ∘ ♭`` over this dual, through the one door that reads
        :func:`~orpheus.numerics.operator.posed`.
        """
        if self._dual_of is not None:
            return self._dual_of  # (Gᵀ)ᵀ = G, by object identity
        if not self.is_adjointable:
            raise MissingAdjoint(
                f"GreenOperator.dual: the Green of the dual sum needs every "
                f"term to transpose and the head's inverse to transpose "
                f"((A^T)^{{-1}} = (A^{{-1}})^T); the head is "
                f"{type(self._terms[0]).__name__}."
            )
        # A Green holds two terms or more, so the reduction is a sum.
        dual_sum = cast(OperatorSum, reduce(operator.add, (t.dual() for t in self._terms)))
        transposed = GreenOperator(dual_sum, max_iter=self.max_iter, tol=self.tol)
        transposed._dual_of = self
        return transposed

    def apply_transpose(self, x: Any, /, *, initial_guess: Any | None = None) -> Any:
        r"""Return :math:`((A - B)^{-1})^{\mathsf T} x` through :meth:`dual`, converged or loud."""
        return self.dual().apply(x, initial_guess=initial_guess)

    def __repr__(self) -> str:
        return (
            f"GreenOperator({self.inner!r}, "
            f"max_iter={self.max_iter}, tol={self.tol})"
        )

r"""The Green inverse operator — the preconditioned-splitting member of the
#226 inverse family (taxonomy §12 step 4).

A general operator SUM has no direct inverse: :math:`(A + B)^{-1}` is not a
function of :math:`A^{-1}` and :math:`B^{-1}` (Sherman–Morrison–Woodbury
applies only under low-rank structure — which the bulk blocks L, C lack
and the boundary block B HAS: rank ≤ N/2 per face, Issue #300), and —
unlike the schedule-triangular
``(L+C)`` — there is no substitution order that solves it in one pass.  What
it DOES have, once it is SPLIT into an invertible implicit part :math:`M` and
a lagged part :math:`N`, is a convergent iteration: write the operator
:math:`M - N`, and

.. math::

    (M - N)^{-1}\,q
    \;=\; \bigl(I - M^{-1}N\bigr)^{-1} M^{-1} q
    \;=\; \sum_{k=0}^{\infty} \bigl(M^{-1}N\bigr)^{k}\,M^{-1}\,q ,

the NEUMANN SERIES of the :math:`M`-preconditioned splitting — in transport,
the multiple-scattering expansion: :math:`M^{-1}q` is the uncollided flux,
:math:`(M^{-1}S)^k M^{-1} q` the k-times-rescattered contribution (Lewis &
Miller §3.2; Adams & Larsen 2002 §II).  Its partial sums are EXACTLY the
Richardson / source-iteration iterates
:math:`x_{n+1} = M^{-1}(q + N\,x_n)` started from zero, convergent iff
:math:`\rho(M^{-1}N) < 1` — for the within-group transport loss
:math:`A_{\rm loss} = (L+C) - S` this is the scattering ratio bound
:math:`\rho \le \max_{\rm cell} c = \Sigma_s/\Sigma_t < 1`, guaranteed for
any physical (absorbing) medium.

**The name is a promise (taxonomy §13).**  :class:`GreenOperator` is the
discrete Green's function of its forward sum — ``G.apply(δ_j)`` is column j
of :math:`(M-N)^{-1}`, the flux response to a unit point source — and the
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
consume directly.  This class reads the splitting off the named split
(the implicit part → preconditioner via its own structure-keyed
``.inverse()`` — the WDD :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`
for ``(L+C)``, an :class:`~orpheus.numerics.operator.InverseOperator` for a
value leaf — the lagged part's terms → gains) and hands it to the driver;
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

**The split is named, never read off a spelling (numerics_untangling.md
ruling 50).**  A sum is a formal linear combination with no order
(:class:`~orpheus.numerics.operator.OperatorSum`), so it has no "leading
term" to precondition with, and a sum of two or more terms has no inverse of
its own.  The Green is built from the split itself:
``GreenOperator(M, N)`` is :math:`(M - N)^{-1}`, iterating
:math:`x \leftarrow M^{-1}(q + N x)`.  Which part is implicit is the
Strategy's choice (``Splitting.from_schedule`` labels the loss terms), not a
property of the spelling:

* ``GreenOperator(L + C, S)`` (the sweep-invertible implicit part, the
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` atom
  its owner built) is the physical splitting: preconditioner the sweep,
  gain :math:`S`, convergence certified by :math:`c < 1`.
* ``GreenOperator(M, N)`` with a non-invertible :math:`M` is REFUSED at
  construction (:class:`~orpheus.numerics.operator.NotInvertible`).
* ``GreenOperator(C, -L)`` (a legal split whose implicit part happens to be
  invertible) constructs — the algebra cannot know
  :math:`\rho(C^{-1}L) > 1` — but the divergent splitting raises
  :class:`ConvergenceFailure` LOUDLY at apply time.  Same operator as the
  sweep split, a different algorithm by the choice of :math:`M`.  Never a
  silent wrong answer (Cardinal Rule 1).

Until ruling 50 the sum's left-spine head was the designated preconditioner
(``OperatorSum.is_invertible`` read the first-spelled term, and
``OperatorSum.inverse`` built this class), so the spelling chose the
algorithm; an unordered sum cannot carry that contract.

**Seed semantics (§38 application context; taxonomy §11.3).**
``apply(q, initial_guess=x₀)`` threads the seed as the splitting
iteration's START — a genuine consumer (contrast
:class:`~orpheus.numerics.operator.InverseOperator`'s accept-and-ignore).
The RESULT is :math:`A^{-1}q` to the configured tolerance regardless of the
seed (the convergence test, not the start, defines the promise), so the
kwarg is context, not identity.

**Module placement (taxonomy §17 W3).**  This is a LEAF module: it imports
both :mod:`~orpheus.numerics.operator` (the algebra) and
:mod:`~orpheus.numerics.iteration` (the drivers), and nothing in the algebra
imports it back: since ruling 50 no sum builds its own Green, so the late
import the sum's ``inverse`` once made is gone.

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
from typing import Any

from orpheus.numerics.iteration import SourceIteration, _l2_norm, seeded_inverse
from orpheus.numerics.operator import (
    InverseWrapMixin,
    LinearOperator,
    MissingAdjoint,
    NotInvertible,
    OperatorSum,
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
    residual :math:`\lVert(M-N)\psi - q\rVert/\lVert q\rVert` still at or
    above ``tol`` (spec §18.A: the promise reads the equation residual,
    never the driver's ρ-blind increment) — or when that residual is
    non-finite (a divergent split overflows the iterate, or, sharper,
    overflows the driver's stopping DENOMINATOR to ``inf`` so it
    "converges" at ``increment = finite/inf = 0.0`` onto a garbage
    iterate).  In both shapes the iterate is NOT the promised
    :math:`A^{-1}q`, and returning it silently would be a wrong answer
    (Cardinal Rule 1).  A divergent splitting means the named implicit
    part :math:`M` is the wrong preconditioner, :math:`\rho(M^{-1}N) \ge 1`
    (see the module docstring's split table); a slowly-converging one means
    the budget (``max_iter``) is too small for the physics
    (:math:`c \to 1`).

    Deliberately NOT raised by the driver itself:
    :class:`~orpheus.numerics.iteration.SourceIteration` returns its final
    iterate with the residual history because nested-iteration consumers
    (the k-eigenvalue inner) legitimately run partial, budget-bounded
    solves.  The PROMISE of an exact inverse — and therefore the raise —
    belongs to the inverse OBJECT.
    """


class GreenOperator(InverseWrapMixin[LinearOperator], LinearOperator):
    r"""The inverse operator :math:`(M - N)^{-1}` of a split operator, realized
    as the convergent :math:`M`-preconditioned splitting.

    Built from the split: ``GreenOperator(M, N)`` for :math:`A = M - N`, with
    :math:`M` the invertible implicit part and :math:`N` the lagged part
    (ruling 50). The forward :attr:`inner` is the combination ``M - N``, so
    :meth:`inverse` returns it (the object-identity involution of the
    family).

    The splitting is fixed at construction (never re-derived per apply): the
    preconditioner is :math:`M`'s own structure-keyed ``.inverse()`` (the
    sweep for the S\ :sub:`N` implicit part, an
    :class:`~orpheus.numerics.operator.InverseOperator` for a value leaf),
    and the driver's gains are :math:`N`'s TERMS, in :math:`N`'s insertion
    order (each ``cᵢwᵢ``), so ``GreenOperator(A, S + F)`` hands the driver
    ``(S, F)`` and associates :math:`(q + S\psi) + F\psi` exactly as a hand
    ``SourceIteration(A⁻¹, S, F)`` does.

    :meth:`apply` returns the CONVERGED :math:`(M - N)^{-1}q` or raises
    :class:`ConvergenceFailure` — never a silent partial iterate.  The
    wrap-delegate back-half (domain↔codomain swap / ``solve`` = the forward's
    matvec / object-identity involution ``inverse() → inner``) is inherited
    from :class:`~orpheus.numerics.operator.InverseWrapMixin`.

    The Euclidean transpose is the Green of the dual split,
    :math:`((M - N)^{-1})^{\mathsf T} = (M^{\mathsf T} - N^{\mathsf T})^{-1}`
    (:meth:`dual`, which :meth:`apply_transpose` applies), with the same
    budget. Its preconditioner is the dual arrow's swap law
    :math:`(M^{\mathsf T})^{-1} = (M^{-1})^{\mathsf T}`, so it exists exactly
    where :math:`M`'s inverse transposes (the SN sweep: the 1-D scan family;
    the multi-D reverse scan is #616). :attr:`is_adjointable` reads that
    condition and that both parts transpose. The Hilbert adjoint is the
    base's ``♯ ∘ dual ∘ ♭``: the law
    :math:`((M - N)^{-1})^{*} = (M^{*} - N^{*})^{-1}` holds and the honesty
    gate checks it independently, but it is not a second route to ``.H``.

    Parameters
    ----------
    implicit : LinearOperator
        :math:`M`, the implicit part; must be invertible (`is_invertible`).
        It becomes the preconditioner through its OWN ``.inverse()``.
    lagged : LinearOperator
        :math:`N`, the lagged part; its terms are the driver's gains.
    max_iter : int, optional
        TOTAL splitting-step budget across :meth:`apply`'s refinement
        loop (default ``1000``).  Exhausting it with the promise unmet
        raises :class:`ConvergenceFailure`.
    tol : float, optional
        The promise tolerance (default ``1e-8``): :meth:`apply` returns
        only when the TRUE relative residual
        :math:`\lVert(M-N)\psi - q\rVert/\lVert q\rVert` is below it
        (spec §18.A — never the driver's ρ-blind increment).  ``tol=0.0``
        is unsatisfiable for an iterative inverse and always raises.

    Raises
    ------
    NotInvertible
        At construction, when :math:`M` is not invertible — or from
        :math:`M`'s own ``.inverse()`` builder (e.g. a singular diagonal),
        where the invertibility obligation lives.
    """

    def __init__(
        self,
        implicit: LinearOperator,
        lagged: LinearOperator,
        *,
        max_iter: int = 1000,
        tol: float = 1e-8,
    ) -> None:
        if not implicit.is_invertible:
            raise NotInvertible(
                f"GreenOperator requires an invertible implicit part M of the "
                f"split A = M − N; got {type(implicit).__name__} with "
                f"is_invertible=False."
            )
        super().__init__(implicit - lagged)
        #: The split :math:`A = M - N`: the implicit part and the lagged part.
        self.implicit = implicit
        self.lagged = lagged
        #: The Green this one is the dual of, when :meth:`dual` built it:
        #: the involution (Gᵀ)ᵀ = G as an object identity.
        self._dual_of: "GreenOperator | None" = None
        #: The promise tolerance: :meth:`apply` returns only when the TRUE
        #: relative residual is below it.  Doubles as the driver's
        #: per-call increment stopping heuristic.
        self.tol = float(tol)
        #: TOTAL splitting-step budget across the refinement loop.
        self.max_iter = int(max_iter)
        # The splitting, fixed ONCE: x_{n+1} = M⁻¹(q + Σ gains·x_n), the gains
        # N's terms — handed to the SAME Richardson driver the solvers consume
        # (§11.2: wrap, never re-implement).
        self._preconditioner = seeded_inverse(implicit)
        self._gains = (
            tuple(c * w if c != 1.0 else w for c, w in lagged.terms)
            if isinstance(lagged, OperatorSum)
            else (lagged,)
        )

    def _driver(self, budget: int) -> SourceIteration[Any]:
        r"""The Richardson driver of this splitting, allowed ``budget`` passes.

        Built per refinement call so the TOTAL budget is exact: each call
        receives what the earlier calls (and a cold start's first partial
        sum) left. Construction stores references only.
        """
        return SourceIteration(self._preconditioner, *self._gains, max_iter=budget, tol=self.tol)

    def _true_relative_residual(self, psi: Any, q: Any) -> float:
        r"""The PROMISE metric :math:`\lVert(M-N)\,\psi - q\rVert / \lVert q\rVert`.

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
        r"""Return :math:`(M - N)^{-1}\,q` with TRUE relative residual
        below ``tol``, or raise :class:`ConvergenceFailure`.

        ``initial_guess`` seeds the splitting iteration's start (§38
        application context — the converged result is seed-independent to
        ``tol`` by the promise check itself); ``None`` starts the Neumann
        series from zero, whose partial sums the G-Neumann gate pins. The
        cold start is taken as the series' first partial sum
        :math:`M^{-1}q`, which is one splitting step from zero
        (:math:`M^{-1}(q + N\,0)`) and is counted as one: the zero itself is
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
                    f"promise — TRUE relative residual ‖(M-N)ψ − q‖/‖q‖ = "
                    f"{true_res:.3e} after {steps} step(s), tol="
                    f"{self.tol:.3e} (driver increment {increment:.3e}); "
                    f"the iterate is NOT (M-N)⁻¹q. A divergent split means "
                    f"the implicit part M is the wrong preconditioner "
                    f"(ρ(M⁻¹N) ≥ 1) — a slow one, that max_iter is too "
                    f"small for the physics."
                )

    @property
    def is_adjointable(self) -> bool:
        # The dual split needs both parts' transposes, and its implicit part
        # (M^T)^{-1} = (M^{-1})^T needs the dual arrow's swap law.
        return (
            adjointable(self.implicit)
            and adjointable(self.lagged)
            and self.implicit.dual().is_invertible
        )

    def dual(self) -> "GreenOperator":
        r"""Return :math:`((M - N)^{-1})^{\mathsf T} = (M^{\mathsf T} - N^{\mathsf T})^{-1}`, the Green of the dual split.

        Both parts' dual arrows, with the same budget; the preconditioner is
        the dual arrow's swap law :math:`(M^{\mathsf T})^{-1} =
        (M^{-1})^{\mathsf T}`. The Hilbert adjoint is NOT overridden:
        ``.H`` is the base's ``♯ ∘ dual ∘ ♭`` over this dual, through the
        one door that reads :func:`~orpheus.numerics.operator.posed`.
        """
        if self._dual_of is not None:
            return self._dual_of  # (Gᵀ)ᵀ = G, by object identity
        if not self.is_adjointable:
            raise MissingAdjoint(
                f"GreenOperator.dual: the Green of the dual split needs both "
                f"parts to transpose and the implicit part's inverse to "
                f"transpose ((M^T)^{{-1}} = (M^{{-1}})^T); the implicit part is "
                f"{type(self.implicit).__name__}."
            )
        transposed = GreenOperator(
            self.implicit.dual(), self.lagged.dual(), max_iter=self.max_iter, tol=self.tol,
        )
        transposed._dual_of = self
        return transposed

    def apply_transpose(self, x: Any, /, *, initial_guess: Any | None = None) -> Any:
        r"""Return :math:`((M - N)^{-1})^{\mathsf T} x` through :meth:`dual`, converged or loud."""
        return self.dual().apply(x, initial_guess=initial_guess)

    def __repr__(self) -> str:
        return (
            f"GreenOperator({self.implicit!r}, {self.lagged!r}, "
            f"max_iter={self.max_iter}, tol={self.tol})"
        )

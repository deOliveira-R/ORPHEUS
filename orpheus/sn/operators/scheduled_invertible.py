r"""The reified splitting matrix :math:`M = (L+C-B_{\rm lower})` — the
schedule-folded sweep-invertible composite (#226 taxonomy step 2, §17 W2).

The boundary Gauss-Seidel is a **matrix splitting** of the
within-group loss:

.. math::

    (L+C-B) \;=\; \underbrace{(L+C-B_{\rm lower})}_{M}
              \;-\; \underbrace{B_{\rm upper}}_{N},
    \qquad
    \psi_{k+1} \;=\; M^{-1}\bigl(q + S\,\psi_k + B_{\rm upper}\,\psi_k\bigr).

⛔ It is a *splitting*, **not** a "regular splitting" — this docstring said
*regular matrix splitting* until 2026-08-09 (#341).  Varga's *regular*
requires ``M⁻¹ ≥ 0`` AND ``N ≥ 0``; the multi-D DD face-to-face
transmission has ``d−1`` eigenvalues of exactly ``−1``, so it is not a
non-negative operator, the comparison theorem is void, and boundary G-S is
NOT guaranteed to be at least as fast as Jacobi (it is measured slower on
some configurations).  Full derivation + measurements:
:ref:`sn-boundary-gs-not-regular`.

``B_lower`` is the whole-trace boundary law masked to the couplings the
octant-group schedule realizes IN-sweep (rows whose octant group is swept
strictly after their face's reflect —
:meth:`~orpheus.sn.operators.boundary.SNBoundaryOperator.split`); the
specular map has no octant-diagonal, so the split is exact and :math:`M`
stays schedule-triangular.  The scheduled walk (the SAME uniform
sweep-and-reflect loop as Jacobi, differing ONLY in the
:class:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule`) is
EXACT forward substitution for :math:`M` — so ``M.solve`` round-trips
``M.apply`` at machine precision, the W2-round-trip keystone.

**What this dissolves.**  The former ``_GaussSeidelResolvent`` paired
``apply = (L+C)ψ`` with ``solve = (L+C−B_lower)⁻¹`` — inverses of DIFFERENT
operators (the §17 falsifier-3 round-trip defect was O(1): 2.667).  Here the
forward and the inverse are the two faces of the ONE reified :math:`M`:
``apply`` is ``(L+C).apply − B_lower.apply``, the sweep atom's own action
minus the lower half, and ``inverse()`` returns a genuine
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` on ``self``.  The
lagged half rides the DRIVER as the ordinary external gain
``(S, B_upper)`` — structurally congruent with the Jacobi path's ``(S, B)``,
so :class:`~orpheus.numerics.iteration.SourceIteration` needs no case split.
The Jacobi schedule is the degenerate member (``B_lower = 0``,
``M = L+C``) — production never constructs it (the factory returns the plain
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` there).

**Construction** — explicit, by the owner that needs the scheduled sweep
(the :class:`~orpheus.sn.splitting.Splitting`, when the schedule labels
``B_lower`` implicit):

.. code-block:: python

    lower, upper = B.split(SweepSchedule.gauss_seidel(
        problem.ndim, problem.quad.octants, reflective_faces(problem),
    ))
    M = ScheduledInvertibleOperator(sco, lower)   # its expansion: L + C − B_lower
    gains = (S, upper)                            # the lagged complement, driver-side

``sco - lower`` is the plain two-term combination: ``-`` never fuses (ruling
52), so the realisation is chosen where ruling 20 puts it, at the Strategy.
The mask carries the schedule it was derived from, so a foreign schedule
cannot be paired with a mismatched mask.

**Substrate, not a twin.**  ``solve`` delegates to the ONE
``StreamingCollisionOperator._solve_timed_full_field`` body with the schedule and the
whole-trace reflect threaded through the representation's own ``sweep`` door
(S6.5 — the scheduled walk runs the operator's ONE
:class:`~orpheus.sn.loss_representation.LossRepresentation` instance and its
own interior kernel).  The former solver-side ``_solve_scheduled`` body — a
plumbing twin of the operator's solve body — retired into this delegation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Final

from orpheus.numerics.operator import FusedAtom, LinearOperator

from .boundary import SNMaskedBoundaryOperator
from .streaming import StreamingCollisionOperator
from .sweep_operator import ScheduleTriangularForward, SweepOperator

if TYPE_CHECKING:
    from orpheus.numerics.frame import FrameBase
    from orpheus.sn.loss_representation import LossRepresentation
    from orpheus.sn.loss_representation.sweep_schedule import (
        SweepSchedule,
        reflective_faces,
    )
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField
    from ..problem import SNProblem
    from .boundary import SNBoundaryOperator


__all__ = ["ScheduledInvertibleOperator"]


class ScheduledInvertibleOperator(FusedAtom["FullField"], ScheduleTriangularForward):
    r"""The sweep-invertible schedule-folded atom :math:`M = (L+C-B_{\rm lower})`.

    An atom with two named legs, :attr:`invertible` (the
    :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`
    :math:`L+C`) and :attr:`lower` (:math:`B_{\rm lower}`), whose algebraic
    content is :meth:`expansion`, the plain combination
    :math:`L + C - B_{\rm lower}` flattened over the original atoms (ruling
    52). It carries the schedule-triangular solve identity a combination
    cannot: the octant-group forward substitution IS the inverse algorithm
    for this specific operator, exactly the way the WDD sweep is for
    :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` (whose
    shape this deliberately mirrors). :meth:`apply` is
    ``invertible.apply − lower.apply``, the sweep atom's own action, so the
    forward and the scheduled sweep invert each other at machine precision.
    """

    def __init__(
        self,
        invertible: "StreamingCollisionOperator",
        lower: "SNMaskedBoundaryOperator",
    ) -> None:
        if not isinstance(invertible, StreamingCollisionOperator):
            raise TypeError(
                f"ScheduledInvertibleOperator: 'invertible' must be an "
                f"StreamingCollisionOperator; got {type(invertible).__name__}."
            )
        if not isinstance(lower, SNMaskedBoundaryOperator):
            raise TypeError(
                f"ScheduledInvertibleOperator: 'lower' must be an "
                f"SNMaskedBoundaryOperator (from SNBoundaryOperator.split); "
                f"got {type(lower).__name__}."
            )
        # Mesh-identity invariant (as for StreamingCollisionOperator): the scheduled
        # walk pairs the mask's row split with the streaming geometry.
        if invertible.problem is not lower.problem:
            raise ValueError(
                "ScheduledInvertibleOperator: the invertible composite and "
                "the boundary mask must act on the same mesh instance "
                "(mesh-identity invariant)."
            )
        #: :math:`L + C`, the sweep atom whose walk the schedule re-orders.
        self.invertible: Final = invertible
        #: :math:`B_{\rm lower}`, the boundary couplings realized in-sweep.
        self.lower: Final = lower
        # The roles of the expansion L + C − B_lower (FULL).
        self._adopt_expansion_roles()
        # The scheduled forward substitution IS this composite's inverse —

    def expansion(self) -> "LinearOperator[FullField]":
        r"""The plain combination :math:`L + C - B_{\rm lower}` this atom realises.

        Flattened over the ORIGINAL atoms (the streaming and collision legs of
        :attr:`invertible`, and :attr:`lower`), never over the sweep atom: an
        expansion is algebraic content, and one holding a realisation would
        refer back to a realisation (the S2 spec, finding 9).
        """
        return self.invertible.expansion() - self.lower

    def apply(self, psi: "FullField") -> "FullField":
        r""":math:`M\psi = (L+C)\psi - B_{\rm lower}\psi`, the sweep atom's own action minus the lower half."""
        return self.invertible.apply(psi) - self.lower.apply(psi)

    def apply_transpose(self, phi: "FullField") -> "FullField":
        r""":math:`M^{\mathsf T}\phi = (L+C)^{\mathsf T}\phi - B_{\rm lower}^{\mathsf T}\phi`."""
        return self.invertible.apply_transpose(phi) - self.lower.apply_transpose(phi)

    @property
    def schedule(self) -> "SweepSchedule":
        """The octant-order schedule — carried by the mask it derived."""
        return self.lower.schedule

    @property
    def boundary(self) -> "SNBoundaryOperator":
        """The whole-trace boundary law ``B`` the in-sweep reflect applies."""
        return self.lower.inner

    @property
    def problem(self) -> "SNProblem":
        """The shared :class:`SNProblem` (validated mesh-identity at init)."""
        return self.invertible.problem

    @property
    def loss_representation(self) -> "LossRepresentation":
        """The operator's ONE representation — the invertible operand's (S6.5)."""
        return self.invertible.loss_representation

    # ── solve: scheduled forward substitution ─────────────────────────

    @property
    def is_invertible(self) -> bool:
        # M is schedule-triangular: the octant-group forward substitution IS
        # its inverse operator — the
        # faithfulness keystone rides with StreamingCollisionOperator's.
        return True

    def inverse(self) -> "SweepOperator":
        r"""Return the inverse OPERATOR :math:`M^{-1}` — a genuine
        :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` on ``self``.

        ``M.inverse().apply(b)`` is the scheduled forward substitution,
        bit-identical to :meth:`solve` — and it round-trips ``M.apply`` at
        machine precision (the W2-round-trip keystone), which the dissolved
        ``_GaussSeidelResolvent`` could not (its forward was ``(L+C)``).
        """
        return SweepOperator(self)

    def solve(
        self,
        rhs: "FullField",
        *,
        initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""Invert :math:`M\,\psi = \text{rhs}` by octant-group forward substitution.

        The inflow rows of ``rhs.boundary`` supply the given data for the
        ``B_upper`` (lagged) couplings — the driver folds
        ``q.boundary + B_upper·ψ_k`` there, exactly as the Jacobi path folds
        ``q.boundary + B·ψ_k``.  The ``B_lower`` rows need NO seed: each is
        read only by octant groups swept after its face's reflect, which
        writes the fresh current-iterate reflection in-sweep.  This is a
        multi-D Cartesian walk — an EXACT direct inverse with no
        previous-iterate seed. ``initial_guess`` is the uniform
        :class:`~orpheus.numerics.iteration.SupportsSeededApply` keyword
        (#285), ACCEPTED and DROPPED; the warm start lives at the iteration
        layer, never on a direct sweep (#280 2.5c).
        """
        del initial_guess  # accept-and-drop: exact direct inverse, nothing to seed (#280 2.5c)
        return self._solve_timed_full_field(rhs)

    def _solve_timed_full_field(
        self,
        rhs: "FullField",
        *,
        moment_frame: "FrameBase | None" = None,
    ) -> "TimedFullField":
        r"""The uniform private solve surface (duck-shared with
        :meth:`StreamingCollisionOperator._solve_timed_full_field`) — the windowed
        product's fused evaluation calls THIS with its ``moment_frame``.

        Delegates to the invertible operand's one body with the schedule and
        the whole-trace in-place reflect threaded; the representation's
        ``sweep`` door injects its OWN interior kernel (S6.5), so the
        scheduled walk is a re-scheduling of THE operator, not a parallel
        construction.

        **Honest domain (ERR-071).** The scheduled walk realizes
        :math:`M^{-1}` exactly on the SOURCE SUBSPACE
        ``{y : y.outflow-rows = 0}`` — which contains every production rhs
        (``q + S·ψ + B_upper·ψ`` write bulk/inflow slots only) and, to FP
        dust, every ``M.apply`` image of a trace-consistent state (the
        round-trip domain).  Off that subspace the mid-march reflect
        consumes UN-restored streamed values (the outflow-defect restore
        in the delegated body fires only after the whole march), so a
        STRUCTURALLY nonzero outflow row leaves the lower-coupled inflow
        rows off by :math:`B(\text{rhs}_{\rm out})` — the full-space
        completion needs the restore interleaved per-group with the
        schedule (subtract each group's seed rows before any reflect
        reads them), a walk-internal carve deferred until a full-space
        consumer exists (e.g. a G-S-preconditioned Krylov posture).  No
        value-guard is placed here (a threshold would be arbitrary and an
        exact-zero test rejects legitimate FP-dust round-trips); the
        honest-scope witness is the W2 off-domain characterization pin
        (``tests/gates/sn/solve/test_gauss_seidel_reification.py``), and the
        production catcher for a future off-domain consumer is the
        end-of-solve convergence-claim check that caught ERR-071 itself.  The bare
        :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` sweep
        (no schedule, no mid-march reader) IS exact on the whole space.
        """
        return self.invertible._solve_timed_full_field(
            rhs,
            moment_frame=moment_frame,
            schedule=self.schedule,
            # The row-masked ADDITIVE update ``bf[lower] += (B·bf)[lower]`` —
            # completes the inhomogeneous lower row ``z_in = y_row + (Bz)_row``
            # on top of the seed, which is what makes this an EXACT inverse
            # of ``apply`` on arbitrary SOURCE-SUBSPACE rhs (the
            # W2-round-trip keystone) and leaves the upper rows carrying
            # the splitting's honest lagged seed.  See
            # :meth:`SNMaskedBoundaryOperator.reflect_rows_inplace`.
            reflect=self.lower.reflect_rows_inplace,
        )

    def __repr__(self) -> str:
        return (
            f"ScheduledInvertibleOperator({self.invertible!r} - {self.lower!r})"
        )

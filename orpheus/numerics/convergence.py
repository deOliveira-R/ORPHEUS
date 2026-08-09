r"""The shared vocabulary for *"did this iteration reach its promise?"*.

An iterative solver can stop for two structurally different reasons, and a
result object that cannot tell them apart is a trap:

* it **converged** — the stopping criterion was met, and the answer is the
  one the caller asked for;
* its **budget ran out** — the loop hit ``max_iter`` with the criterion
  unmet, and the answer is a *best-effort iterate*, mid-descent.

Both come back through the same return type, so unless the distinction is
carried explicitly a caller cannot tell a converged answer from an arbitrary
point on the way to one.

Why this module exists — the measured failure
=============================================

`[M]` 2026-08-08.  ``test_d3_pure_absorber_per_ordinate_psi_exact`` asserted
a closed-form identity to ``rtol=1e-10`` on an all-reflective 3-D box.  That
problem needs **1631** source-iteration sweeps; ``max_inner`` defaulted to
**1000**.  The solve returned the 999th iterate — honestly flagged
``history.converged = False`` — and the gate, which never read the flag,
asserted physics against it.  The error was ``3.287e-10``, and it was
*bit-identical* at every ``inner_tol`` from ``1e-9`` to ``1e-15``, because
the running residual (``1.185e-09``) never fell below even the loosest of
them: all four runs hit the same cap and returned the same bytes.

Two lessons are worth carrying forward, both counter-intuitive:

#. **Tolerance-insensitivity is the signature of a budget truncation, not of
   a discretization floor.**  If an error does not move when you tighten the
   tolerance, the first thing to check is whether the tolerance ever *bound*
   — read the iteration count against the cap before concluding anything
   about the discretization.
#. **A value gate that does not assert convergence is asserting an arbitrary
   iterate.**  It will be green or red by luck.  That gate had ridden an
   unconverged exit since it was written; it passed for months because the
   truncated error happened to land inside the tolerance, and a *correct*
   quadrature change (#337) later moved it out.

The remedy has two halves.  The **structural** half is that a loop which
knows why it stopped must say so in its return value — see
:class:`~orpheus.numerics.eigenvalue.PowerIterationOutcome`, and the
``converged`` field every solver result now derives rather than asserts.
The **loudness** half is this module's warning.

Why a warning and not an exception
==================================

Raising would be louder, and it is the wrong default here.  The project has
already taken this decision once, for the same defect class: when
``KrylovAcceleration`` stopped discarding scipy's ``info`` flag
(``orpheus/numerics/iteration.py``, D-H.1e / ERR-053), the ruling was

    *"Both surface as warnings — raising would break long-standing callers
    that tolerate slow convergence and need the best-effort iterate."*

That reasoning holds at the public solver entries too: rate studies harvest
the residual history rather than the answer, spy tests discard the result
entirely, and a diverging configuration is sometimes exactly what is being
measured.  `[M]` an audit of the tree found **zero** production callers and
**zero** ``examples/`` callers that depend on a truncated answer, but six
legitimate test consumers that do.

The warning is nonetheless *escalatable*, which is what makes it a gate
rather than a decoration.  ``pyproject.toml`` sets no ``filterwarnings`` and
no ``-W error``, so emitting it costs nothing today, and CI (or any careful
caller) turns the whole category into a hard failure with one flag::

    python -O -m pytest -W error::ConvergenceWarning

That is the same named-warning + opt-in-escalation recipe the regression
suite already uses for its bit-identity tripwire.

Related, and deliberately NOT merged with this
==============================================

* :class:`~orpheus.sn.solver.ConvergenceCertificateError` asserts a
  *different* proposition — *"a convergence claim was made and it was
  FALSE"* (the in-M lag-death class).  It stays a hard error: a false claim
  is a bug, whereas an honest best-effort answer is a legitimate result.
* :class:`~orpheus.numerics.green_operator.ConvergenceFailure` is the
  *raising* form of this module's proposition, scoped to ``GreenOperator``'s
  own promise.  It keeps its Green-operator provenance; a caller who wants
  hard failure at a solver entry escalates the warning instead.
"""

from __future__ import annotations

__all__ = ["ConvergenceWarning"]


class ConvergenceWarning(RuntimeWarning):
    r"""An iterative solve exhausted its budget; the answer is best-effort.

    Emitted by a public solver entry when the returned result carries
    ``converged = False`` — i.e. the loop hit its iteration cap with the
    stopping criterion unmet.  The result is still returned (it is the best
    iterate available), which is why this is a warning and not an exception;
    see the module docstring for the ERR-053 precedent behind that choice.

    The message states the budget that ran out, the tolerance that was not
    reached, and the last residual, so the reader can tell *how far* from
    converged the answer is — the distance between "one more sweep" and
    "diverging" is the whole diagnostic content.

    Escalate to a hard failure with ``-W error::ConvergenceWarning``.
    """

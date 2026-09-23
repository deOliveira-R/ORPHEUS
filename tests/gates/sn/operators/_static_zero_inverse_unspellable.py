# pyright: reportUnnecessaryTypeIgnoreComment=error
r"""STATIC designed-error pin: ``ZeroOperator().inverse()`` is UNSPELLABLE.

Carve P4 (taxonomy §12 step 6, verification spec §39.3) — Design C's
STRUCTURAL arm: a type with no mathematical inverse declares NO
``inverse()`` at all, so misuse is a *static* error at the call site,
never a runtime raise.  This file pins that fact with INVERTED polarity
so the checked tree stays noise-free (#226):

* today the expression below is a genuine ``reportAttributeAccessIssue``,
  so the rule-scoped ignore on it is NECESSARY and this file checks
  CLEAN under ``npx pyright``;
* if anyone ever declares ``inverse()`` on ``ZeroOperator`` — or hosts
  it on the ``LinearOperator`` base, the REJECTED demote-to-runtime arm
  of the step-6 false dichotomy — the expression starts type-checking,
  the ignore turns UNNECESSARY, and the file-level
  ``reportUnnecessaryTypeIgnoreComment=error`` pragma (line 1) REDs
  this file.

The ignore comment is therefore the ASSERTION mechanism, not a
suppression of type debt: delete it → this file reds today; break the
unspellability → this file reds then.  (Mutation-verified both ways —
M-STATIC-UNSPELLABLE in the step-6 tooth bank.)

Leading underscore → pytest never collects this module, and the
forbidden expression lives inside a never-called function so an
accidental import cannot execute it.  Runtime twin: the keystone-v2
``STRUCTURAL_ABSENT`` leg (``tests/_harness/predicates.py``) pins
``not hasattr(ZeroOperator(), "inverse")``.
"""
from orpheus.numerics.operator import ZeroOperator


def _zero_inverse_is_unspellable() -> None:
    ZeroOperator().inverse()  # pyright: ignore[reportAttributeAccessIssue]

"""A withdrawn reference generator: the value, the opt-in, and the lock.

A *withdrawal* is a maintainer ruling that a reference generator is not
fit to be believed — not research grade — and so must neither be cited
as evidence nor run by default, until the issue that records the ruling
closes. It is one value, :class:`Withdrawal` ``(reason, issue)``, read in
two places that must agree:

* **statically**, by the test harness: a test that consumes a withdrawn
  generator carries ``@pytest.mark.withdrawn(reason, issue=N)``;
  ``tests/conftest.py`` parses the marker with :meth:`Withdrawal.from_mark`
  at collection time, skips the test with its reason and issue, and
  records the withdrawal on the test's registry entry, so the V&V matrix
  and the error catalogue count a withdrawn test as neither verifying nor
  catching;
* **at run time**, by the generator itself: each withdrawn symbol is
  decorated with :func:`withdrawn_generator`, which refuses the call with
  :class:`GeneratorWithdrawn` unless the issue is lifted. The lock is what
  makes the marker placement a checked claim: an unmarked test that
  reaches a withdrawn generator is red, a module reload cannot undo the
  lock (it is part of the function object), and a subprocess cannot
  escape it (the opt-in travels in the environment).

The opt-in is ``ORPHEUS_RUN_WITHDRAWN``: a comma-separated list of issue
numbers, or ``all``. It is read at CALL time (:func:`lifted_withdrawals`),
so a variable set after import takes effect, and a subprocess inherits
it. Naming the issue keeps a lift scoped: lifting #506 lifts nothing
filed under another issue. The variable permits a run; it never changes
a value, so it is not part of any reference's cache key.

This module imports nothing outside the standard library and nothing
from ``pytest``: :meth:`Withdrawal.from_mark` reads a marker through its
two public fields (``args``, ``kwargs``), so the production lock and the
test hook share ONE definition of the value (X4) without production
depending on the test runner.

The plan of record is ``.claude/plans/reference_cache.md`` (phase P0,
the lock; phase P4, the ``Withdrawn`` state of the reference certificate
that retires it) and its specification ``.claude/plans/reference_p0_spec.md``
§2.
"""

from __future__ import annotations

import functools
import os
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Final, ParamSpec, Protocol, TypeVar, final

__all__ = [
    "ALL_WITHDRAWALS",
    "AllWithdrawals",
    "GeneratorWithdrawn",
    "RUN_WITHDRAWN_VARIABLE",
    "Withdrawal",
    "lifted_withdrawals",
    "withdrawal_of",
    "withdrawn_generator",
]

#: The environment variable that lifts a withdrawal: comma-separated issue
#: numbers (``506`` or ``506,512``), or ``all``.
RUN_WITHDRAWN_VARIABLE: Final = "ORPHEUS_RUN_WITHDRAWN"


class _MarkLike(Protocol):
    """The two fields of a ``pytest`` ``Mark`` that :meth:`Withdrawal.from_mark` reads."""

    @property
    def args(self) -> tuple[Any, ...]: ...

    @property
    def kwargs(self) -> Mapping[str, Any]: ...


@dataclass(frozen=True)
class Withdrawal:
    """A ruling that a reference generator is withdrawn, and the issue that records it.

    Parameters
    ----------
    reason
        Why the generator is withdrawn, in one sentence. Non-empty.
    issue
        The GitHub issue that records the ruling and its return criterion.
        A positive ``int`` (a ``bool`` or a numeric string is refused).

    Raises
    ------
    ValueError
        On an empty reason or an issue that is not a positive ``int``;
        the message names the offending field.
    """

    reason: str
    issue: int

    def __post_init__(self) -> None:
        if not isinstance(self.reason, str) or not self.reason.strip():
            raise ValueError(
                f"Withdrawal.reason must be a non-empty string, got {self.reason!r}"
            )
        if type(self.issue) is not int:
            raise ValueError(
                "Withdrawal.issue must be an int (the GitHub issue number), got "
                f"{type(self.issue).__name__} {self.issue!r}"
            )
        if self.issue <= 0:
            raise ValueError(
                f"Withdrawal.issue must be a positive issue number, got {self.issue}"
            )

    @classmethod
    def from_mark(cls, mark: _MarkLike, *, where: str) -> Withdrawal:
        """Parse ``@pytest.mark.withdrawn(reason, issue=N)`` into a :class:`Withdrawal`.

        ``where`` is the node id of the test carrying the marker; every
        refusal names it, so a malformed marker is a collection error that
        says which test to fix.

        Raises
        ------
        ValueError
            When the marker does not carry exactly one positional reason and
            an ``issue=`` keyword (and nothing else), or when the values fail
            the constructor's law.
        """
        spelling = "@pytest.mark.withdrawn(reason, issue=N)"
        if len(mark.args) != 1:
            raise ValueError(
                f"{where}: a withdrawn marker takes exactly one positional "
                f"argument, the reason; got {len(mark.args)} ({spelling})"
            )
        if "issue" not in mark.kwargs:
            raise ValueError(
                f"{where}: a withdrawn marker must name its issue ({spelling})"
            )
        unexpected = sorted(set(mark.kwargs) - {"issue"})
        if unexpected:
            raise ValueError(
                f"{where}: a withdrawn marker takes only issue=; got {unexpected} "
                f"({spelling})"
            )
        try:
            return cls(reason=mark.args[0], issue=mark.kwargs["issue"])
        except ValueError as exc:
            raise ValueError(f"{where}: {exc}") from exc



@final
class AllWithdrawals:
    """The lift ``ORPHEUS_RUN_WITHDRAWN=all``: it contains every issue."""

    def __contains__(self, issue: object) -> bool:
        return True

    def __repr__(self) -> str:
        return "ALL_WITHDRAWALS"


#: The one instance of :class:`AllWithdrawals`.
ALL_WITHDRAWALS: Final = AllWithdrawals()


def lifted_withdrawals() -> frozenset[int] | AllWithdrawals:
    """The issues ``ORPHEUS_RUN_WITHDRAWN`` lifts, read from the environment NOW.

    Unset or empty lifts nothing. ``all`` lifts every withdrawal. Otherwise
    the value is a comma-separated list of positive issue numbers.

    Raises
    ------
    ValueError
        On any entry that is not a positive integer: a mistyped lift is
        refused, never read as "lift nothing" (a silent no-op would leave
        the caller believing the generator ran).
    """
    raw = os.environ.get(RUN_WITHDRAWN_VARIABLE, "").strip()
    if not raw:
        return frozenset()
    if raw.lower() == "all":
        return ALL_WITHDRAWALS
    issues: set[int] = set()
    for entry in raw.split(","):
        token = entry.strip()
        if not token.isdecimal() or int(token) <= 0:
            raise ValueError(
                f"{RUN_WITHDRAWN_VARIABLE}={raw!r}: {token!r} is not a positive "
                "issue number (expected comma-separated issue numbers, or 'all')"
            )
        issues.add(int(token))
    return frozenset(issues)


class GeneratorWithdrawn(BaseException):
    """A withdrawn reference generator was called without its issue lifted.

    Derives from :class:`BaseException`, not :class:`Exception`, on
    purpose. The refusal is a policy (a lockout), not a failure of the
    computation, and no fallback may absorb it: an ``except Exception``
    around a generator call (a root-finder loop that skips a failed
    bracket, a registry walk that skips a broken producer, a test's
    ``pytest.raises(Exception)``) would otherwise turn the lock into a
    silent fallback or a green test. That is the same reason
    ``KeyboardInterrupt`` and ``SystemExit`` are ``BaseException``. The
    cost is that an ``except Exception`` cleanup block does not run for
    it; a ``finally`` block does. ``pytest`` reports it as an ordinary
    failure (it re-raises only ``Exit`` and ``KeyboardInterrupt``).

    It pickles (``__reduce__``), so a refusal raised in a worker process
    reaches the parent as itself rather than as a pickling error.

    Attributes
    ----------
    withdrawal
        The ruling that withdrew the generator.
    generator
        The qualified name of the refused callable.
    invalid_lift
        ``None`` for an ordinary refusal; otherwise the reason the
        ``ORPHEUS_RUN_WITHDRAWN`` value could not be read (a mistyped lift
        refuses the call, naming the bad value, rather than running it).
    """

    def __init__(
        self, withdrawal: Withdrawal, generator: str, invalid_lift: str | None = None
    ) -> None:
        self.withdrawal = withdrawal
        self.generator = generator
        self.invalid_lift = invalid_lift
        if invalid_lift is not None:
            message = (
                f"{generator} is withdrawn (#{withdrawal.issue}) and the lift could "
                f"not be read: {invalid_lift}"
            )
        else:
            message = (
                f"{generator} is withdrawn (#{withdrawal.issue}): {withdrawal.reason}. "
                f"A test that reaches it carries "
                f"@pytest.mark.withdrawn({withdrawal.reason!r}, issue={withdrawal.issue}); "
                f"to run it anyway set {RUN_WITHDRAWN_VARIABLE}={withdrawal.issue}."
            )
        super().__init__(message)

    def __reduce__(
        self,
    ) -> tuple[type[GeneratorWithdrawn], tuple[Withdrawal, str, str | None]]:
        return (type(self), (self.withdrawal, self.generator, self.invalid_lift))


#: The attribute a locked callable carries its :class:`Withdrawal` on; read
#: through :func:`withdrawal_of`, never spelled elsewhere.
_WITHDRAWAL_ATTRIBUTE: Final = "__withdrawal__"


def withdrawal_of(obj: object) -> Withdrawal | None:
    """The :class:`Withdrawal` a callable is locked behind, or ``None``."""
    withdrawal = getattr(obj, _WITHDRAWAL_ATTRIBUTE, None)
    return withdrawal if isinstance(withdrawal, Withdrawal) else None


_P = ParamSpec("_P")
_R = TypeVar("_R")


def withdrawn_generator(
    withdrawal: Withdrawal,
) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]:
    """Lock a reference generator behind ``withdrawal``: the call is refused unless lifted.

    The wrapped callable keeps its name, docstring, signature and
    ``__wrapped__`` (``functools.wraps``), so ``inspect.signature``, a
    registry keyed on names and a ``monkeypatch`` by attribute all see the
    generator as before. :func:`withdrawal_of` reads the lock back, so the
    set of locked symbols is itself enumerable.

    A call reads ``ORPHEUS_RUN_WITHDRAWN`` NOW (:func:`lifted_withdrawals`);
    a value that does not parse refuses the call with a
    :class:`GeneratorWithdrawn` naming the bad value.

    **ELEGANCE-DEBT[guard] #506** — a run-time refusal stands where the
    reference machinery cannot yet say "this reference is withdrawn" as a
    value; it retires at phase P4 of ``.claude/plans/reference_cache.md``,
    when the ``ReferenceCertificate`` carries the ``Withdrawn(reason,
    issue)`` state and a withdrawn generator is never called.
    """

    def lock(generator: Callable[_P, _R]) -> Callable[_P, _R]:
        name = f"{generator.__module__}.{generator.__qualname__}"

        @functools.wraps(generator)
        def locked(*args: _P.args, **kwargs: _P.kwargs) -> _R:
            try:
                lifted = lifted_withdrawals()
            except ValueError as exc:
                raise GeneratorWithdrawn(withdrawal, name, invalid_lift=str(exc)) from exc
            if withdrawal.issue not in lifted:
                raise GeneratorWithdrawn(withdrawal, name)
            return generator(*args, **kwargs)

        setattr(locked, _WITHDRAWAL_ATTRIBUTE, withdrawal)
        return locked

    return lock

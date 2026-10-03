r"""A withdrawal: a maintainer ruling that a reference must not be believed (#405 P0; moved here at P2 step 5).

A :class:`Withdrawal` ``(reason, issue)`` records that a reference (a
derivation, or a published solution in erratum) is not fit to be cited as
evidence until the issue that records the ruling closes. It is one value
read in three places that must agree (X4):

* the test harness parses ``@pytest.mark.withdrawn(reason, issue=N)`` into
  it (:meth:`Withdrawal.from_mark`) and skips the test with its reason;
* the derivations' run-time lock (``orpheus.derivations.common.withdrawal``,
  :func:`~orpheus.derivations.common.withdrawal.withdrawn_generator`)
  refuses a withdrawn generator's call unless its issue is lifted;
* a :class:`~orpheus.reference.published.PublishedSolution` and a
  :class:`~orpheus.reference.certificate.ReferenceCertificate` carry it as
  their standing, so a withdrawn reference reads as withdrawn wherever it is
  held.

It lives in the input-tier reference package so that all three share the
one class: the derivations import this package, never the converse. It
imports nothing from ``pytest``: :meth:`Withdrawal.from_mark` reads a marker
through its two public fields.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Protocol, final

from orpheus.numerics.content import ContentIdentity, content_digest

__all__ = ["Withdrawal"]


class _MarkLike(Protocol):
    """The two fields of a ``pytest`` ``Mark`` that :meth:`Withdrawal.from_mark` reads."""

    @property
    def args(self) -> tuple[Any, ...]: ...

    @property
    def kwargs(self) -> Mapping[str, Any]: ...


@final
@dataclass(frozen=True, eq=False)
class Withdrawal(ContentIdentity):
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
        content_digest(self)

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

r"""``BoundOperator`` — every transport operator carries its TWO spaces.

The CS4c binding base (design record
``.claude/plans/cs4c_binding_design.md`` §1, ruled 2026-08-30): a bound
operator is an ARROW, and an arrow has two ends. The base carries exactly
what EVERY binding shares — the user's articulation: *"StreamingOperator is
also a BoundOperator, just bound to other data"* — so it holds the two
mandatory :class:`~orpheus.numerics.space.FunctionSpace` ends and the
per-END admission helper, and nothing else. Each channel declares its own
datum fields (C: a coefficient; S: a kernel + minted faces; F: a kernel;
L, later: the closure pair) and its own extract-and-mint classmethods
(the three-tier construction discipline, §3 of the record).

Why the ends are **kw-only mandatory fields** (ruled at the step-1+2
checkpoint): the datum stays positional-first at every call site, and the
domain/codomain SWAP — the ERR-002/ERR-076 transposition family, which
type-checks and produces a well-formed *reversed* arrow that only a
reciprocity gate would catch — is unspellable-silently, because both ends
are NAMED at every exact-ctor site. Endomorphism sugar (one ``space=``
supplying both ends) lives on the TIER-2 classmethods only, never here.

Why the ends are **write-once** (:meth:`BoundOperator.__setattr__`): the
one soundness gap of realizing an abstract read-only property as a
dataclass field is mutability — ``op.domain = other`` would silently
re-bind an arrow's end, the exact illegal state this base exists to make
unrepresentable. The guard refuses re-binding after construction;
``dataclasses.replace`` (which builds a NEW instance through
``__init__``, re-running every admission) is the sanctioned respelling.

The **energy-conformity admission** runs per END through
:meth:`BoundOperator._assert_energy_extent_both_ends`: the shared guard
BODY stays single-sourced in
:mod:`~orpheus.transport.operators._energy_conformity` and each channel
still wires its own operand (data, never metadata — CS4a-R EE-3); what
the base adds is the codomain leg, collapsed onto one check when the
binding is an endomorphism. A 47g→2g condensation binding — EE-8's
recorded limit, the reason one ``(space, ng)`` check could never be the
whole story — runs both legs against its own data extent.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Generic

from orpheus.numerics.operator import Codomain, Domain, LinearOperator
from orpheus.transport.operators._energy_conformity import (
    assert_energy_extent_conforms,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = ["BoundOperator"]


@dataclass(eq=False)
class BoundOperator(
    LinearOperator[Domain, Codomain], Generic[Domain, Codomain],
):
    r"""The two-space binding base — see the module docstring.

    The two fields realize the base's abstract ``domain`` / ``codomain``
    properties as write-once data: sound at runtime (the ``__setattr__``
    guard makes each end single-assignment), invisible to pyright's
    property-vs-attribute override rule — hence the two targeted
    suppressions, each carrying this rationale.
    """

    #: The function space this operator consumes — MANDATORY, kw-only
    #: (the swap-unspellability ruling). Non-Optional by annotation: the
    #: R1 ledger row reads this declaration off the MRO.
    domain: "FunctionSpace" = field(kw_only=True)  # pyright: ignore[reportIncompatibleMethodOverride]

    #: The function space this operator produces — MANDATORY, kw-only.
    #: Distinct from ``domain`` on a genuine non-endomorphism (a
    #: condensation binding; the 2-D S witness of the step-0 census).
    codomain: "FunctionSpace" = field(kw_only=True)  # pyright: ignore[reportIncompatibleMethodOverride]

    def _assert_energy_extent_both_ends(
        self, data_ng: int, *, operator: str,
    ) -> None:
        r"""Per-END energy-conformity admission (CS4c §1).

        Runs the single-sourced guard against the DOMAIN's energy axis
        and — when the binding is not an endomorphism — against the
        CODOMAIN's too. The operand is the CHANNEL's own data extent
        (wired per site, data-not-metadata: CS4a-R EE-3); the collapse
        for endomorphisms keeps the common case one check.
        """
        assert_energy_extent_conforms(self.domain, data_ng, operator=operator)
        if self.codomain is not self.domain and self.codomain != self.domain:
            assert_energy_extent_conforms(
                self.codomain, data_ng, operator=operator,
            )


def _write_once_end(name: str) -> property:
    r"""A write-once property realizing one arrow end.

    Installed onto :class:`BoundOperator` AFTER the ``@dataclass``
    decorator has run — the one order that satisfies all three parties
    (`[M]` 2026-08-30, prototyped both ways):

    * ``dataclass`` collects ``domain``/``codomain`` from the class-body
      field sentinels and generates the kw-only mandatory ``__init__`` —
      then DELETES the no-default sentinels and re-runs
      :func:`abc.update_abstractmethods`, which would re-expose
      :class:`~orpheus.numerics.operator.LinearOperator`'s ABSTRACT
      properties and make every leaf uninstantiable;
    * the abstract property is a data descriptor, so it would also
      shadow the instance attribute ``__init__`` writes;
    * installing THIS concrete property afterwards closes both: the
      generated ``self.domain = ...`` lands in the setter below (first
      write passes, re-binding refuses — an arrow's end is WRITE-ONCE;
      ``dataclasses.replace`` is the sanctioned respelling, re-running
      every admission), and the second
      :func:`abc.update_abstractmethods` call re-derives an EMPTY
      abstract set from it.
    """
    storage = "_" + name

    def fget(self: "BoundOperator") -> "FunctionSpace":
        return getattr(self, storage)

    def fset(self: "BoundOperator", value: "FunctionSpace") -> None:
        if hasattr(self, storage):
            raise AttributeError(
                f"{type(self).__name__}.{name} is write-once — an arrow's "
                f"end cannot be re-bound after construction; build a new "
                f"binding (dataclasses.replace re-runs the admissions)."
            )
        object.__setattr__(self, storage, value)

    fget.__doc__ = (
        f"The binding's {name} — a mandatory, write-once "
        f":class:`~orpheus.numerics.space.FunctionSpace` end "
        f"(see the module docstring)."
    )
    return property(fget, fset)


BoundOperator.domain = _write_once_end("domain")  # type: ignore[assignment]
BoundOperator.codomain = _write_once_end("codomain")  # type: ignore[assignment]
abc.update_abstractmethods(BoundOperator)

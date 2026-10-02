r"""``Materials`` — stage 1 of the posing filtration: the declaration.

The posing chain (charter: ``.claude/plans/posing_filtration_charter.md``,
ratified 2026-08-25) opens with a **declaration**: the materials of THIS
problem, ``{id → Mixture}``, each entry retaining its own raw grid
provenance (:attr:`~orpheus.data.macro_xs.mixture.Mixture.eg` — no early
collapse, because a later method head may unionize; two heads over one
declaration may legitimately realize *different* energy axes). This is
NOT the library: the library is everything you own; ``Materials`` is what
this problem declares. The distinction is class semantics, not a package
boundary — the declaration lives beside the data it declares.

What this class deliberately does NOT carry:

* **No ``ng``.** A scalar group count is only well-defined in the uniform
  multigroup regime, and resolution of concrete properties is LAZY — each
  consuming stage resolves what it consumes (the mesh derives its ``ng``
  where one uniform structure is its own contract; the fine/ACE regime
  never reads a scalar ``ng`` at all, it unionizes). The future per-kind
  consistency checks (GENDF / PENDF / collapsed) land with the data-layer
  overhaul, not here.
* **No energy-axis preview.** The coherence law lives wholly in the mint
  (:meth:`~orpheus.numerics.axis.EnergyAxis.from_materials`, called at
  the consuming stage over the *reachable* set); a preview that does not
  exist cannot disagree with the law. Admission refuses only
  assignment-independent defects (the empty declaration) — a
  cross-material refusal here would let a declared-but-unassigned
  spectator flip admissibility, violating the leak principle.

The one operation the chain adds is :meth:`restrict` — the
reachable-subset constructor. Where an assignment meets the declaration
(the mesh pullback today; the d≥2 overlay object when it lands), an
assigned-but-undeclared id REFUSES in the declaration's vocabulary, and
the leak principle's "reachable materials" becomes a first-class value
the mint consumes.
"""

from __future__ import annotations

import operator
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from orpheus.numerics.content import ContentIdentity, FrozenMapping

if TYPE_CHECKING:
    from collections.abc import ItemsView, KeysView, ValuesView

    from orpheus.data.macro_xs.mixture import Mixture

__all__ = ["Materials"]


@dataclass(frozen=True, eq=False)
class Materials(ContentIdentity):
    r"""The materials of THIS problem — ``{id → Mixture}``, a declaration.

    A frozen VALUE with content identity
    (:class:`~orpheus.numerics.content.ContentIdentity`, #405 P1 step 5,
    2026-10-02): two declarations of equal mixtures under equal ids are
    equal and hash alike, whatever objects they are and in whatever order
    the ids were written, so a reference cache can key on a declaration.
    Until then it compared by identity (``eq=False``) and could not be
    pickled (it held a read-only proxy). The ids are coerced to ``int`` at
    admission (an ``np.int64`` id is the same id; a ``bool`` is refused),
    and the mapping is stored as a
    :class:`~orpheus.numerics.content.FrozenMapping`, keeping the declared
    order (iteration order is behaviour, :meth:`restrict` documents it; it
    is not content), so the declaration cannot be mutated after the fact —
    later stages may safely hold it, and it pickles.

    Parameters
    ----------
    mixtures : Mapping[int, Mixture]
        The declaration, keyed by the integer ids the geometry
        assignment will reference. Refused when empty — a problem with
        no materials has no stage 1. Per-entry grid provenance rides
        each ``Mixture`` (``eg`` or ``None``) untouched.
    """

    mixtures: "Mapping[int, Mixture]"

    def __post_init__(self) -> None:
        if not self.mixtures:
            raise ValueError(
                "Materials requires a non-empty materials declaration; "
                "a problem with no declared mixtures has no stage 1."
            )
        declared = {}
        for mid, mix in self.mixtures.items():
            if isinstance(mid, bool):
                raise TypeError(
                    f"Materials: a material id is an int, got the bool {mid!r}"
                )
            declared[operator.index(mid)] = mix
        object.__setattr__(self, "mixtures", FrozenMapping(declared))

    @classmethod
    def of(
        cls, materials: "Materials | Mapping[int, Mixture]",
    ) -> "Materials":
        """Parse-at-boundary coercion: admit a bare mapping, once.

        The consuming boundary (``MaterialMesh``) accepts either the
        declaration object or the plain ``{id: Mixture}`` dict every
        call site has always written, and normalizes HERE — downstream
        of this parse the type is ``Materials``, so the admission runs
        exactly once (coding-elegance Pattern 4: parse, don't validate).
        """
        if isinstance(materials, Materials):
            return materials
        return cls(materials)

    # ── Mapping surface ───────────────────────────────────────────────
    # The reads every existing consumer performs ([], in, iter, len,
    # keys/values/items, get) — a declaration is read like the dict it
    # replaced, so the parse-at-boundary migration is consumer-invisible.

    def __getitem__(self, mid: int) -> "Mixture":
        return self.mixtures[mid]

    def __iter__(self) -> Iterator[int]:
        return iter(self.mixtures)

    def __len__(self) -> int:
        return len(self.mixtures)

    def __contains__(self, mid: object) -> bool:
        return mid in self.mixtures

    def keys(self) -> "KeysView[int]":
        return self.mixtures.keys()

    def values(self) -> "ValuesView[Mixture]":
        return self.mixtures.values()

    def items(self) -> "ItemsView[int, Mixture]":
        return self.mixtures.items()

    def get(self, mid: int, default=None):
        return self.mixtures.get(mid, default)

    @property
    def ids(self) -> frozenset[int]:
        """The declared ids, as a set (the assignment's admissible range)."""
        return frozenset(self.mixtures)

    def __repr__(self) -> str:
        return f"Materials(ids={sorted(self.mixtures)})"

    # ── The chain's operation ─────────────────────────────────────────

    def restrict(self, ids: "Iterable[int]") -> "Materials":
        r"""The reachable-subset constructor — the assigned-undeclared guard.

        Called where an assignment meets the declaration (today: the
        ``MaterialMesh`` pullback validates its ``mat_map`` through this;
        the d≥2 overlay object when it lands). An id absent from the
        declaration refuses in the declaration's vocabulary — the
        fail-early guard the chain places at the earliest point the
        predicate is decidable. The returned subset preserves the given
        id order (first occurrence), so a sorted reachable set feeds the
        energy-arm rule deterministically.

        A declared-but-unassigned entry is simply *not selected* —
        legal and inert by the leak principle, never warned about.
        """
        wanted = [int(i) for i in ids]
        missing = sorted({i for i in wanted if i not in self.mixtures})
        if missing:
            raise ValueError(
                f"Materials.restrict: the assignment references material "
                f"ids {missing} that are NOT in the declaration "
                f"(available ids: {sorted(self.mixtures)}; used ids: "
                f"{sorted(set(wanted))})."
            )
        return Materials({i: self.mixtures[i] for i in wanted})

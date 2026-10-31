r"""The absorbing chain of the energy groups: which classes of a closed domain can never lose a neutron.

The multigroup loss :math:`A` is the generator of a continuous-time absorbing
Markov chain on the energy groups (and space and angle): a neutron in group
:math:`g` is removed (captured, fissioned, leaks) or transferred to a group
:math:`g'` it scatters into. Group the energy states by the transfer graph
(an edge :math:`g \to g'` when some present material scatters, or
(n,2n)-transfers, :math:`g` into :math:`g'`). A FINAL communicating class
:math:`C` is one no neutron leaves by transfer. If, in addition,

1. every group of :math:`C` has zero NET REMOVAL in every present material,
   :math:`\Sigma_t - \sum_{g'}\Sigma_{s0}(g \to g') - 2\sum_{g'}\Sigma_{2}(g \to g') = 0`
   (the removal the loss operator itself reads, summed exactly), and
2. every face of the domain conserves (a deck, or a condition whose return
   fraction :math:`J^-/J^+` is exactly 1; an infinite medium has no faces),

then a neutron that reaches :math:`C` is never lost, and :math:`A` is singular:
the constant (per class, in the method's representation) is in its kernel. The
converse holds on a connected domain: a class with removal anywhere, or a face
that leaks, makes the class transient. This is the STRUCTURAL singularity of
the loss, decided from declared data alone (the mixtures and the binding),
never from a float difference of an assembled operator.

What is singular is then answered per posing by the Fredholm alternative
(:mod:`orpheus.numerics.fredholm`): each method mints the kernel's bases from
its own loss on the classes named here.

The predicate reads exact zeros of declared data. A removal that is merely
small (:math:`\Sigma_a = 10^{-14}`) is regular and ill-conditioned, solved and
its conditioning recorded; a declared balance that fails at round-off
(:math:`\Sigma_t = 0.1 + 0.2`, :math:`\Sigma_s = 0.3`) is regular on the data
as stored (#665 tracks the data defect). A group with NEGATIVE net removal
(net (n,2n) production) is not absorption-free; a class mixing positive and
negative removals can be singular by a coincidence of values, which is a value
question with no structural witness, not this predicate's.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Iterable, Optional

import numpy as np

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.transport.boundary_conditions.binding import Binding

__all__ = ["absorption_free_closed_classes", "closed_faces", "final_classes"]


def final_classes(transfers: np.ndarray) -> tuple[tuple[int, ...], ...]:
    r"""The final communicating classes of a directed graph on the groups, each sorted, in order of first member.

    ``transfers[g, h]`` is true when :math:`g` feeds :math:`h`. A class is the
    set of groups mutually reachable from its members; it is final when no
    member reaches a group outside it. The transitive closure is Warshall's,
    which is exact on a boolean matrix of the group count's size.
    """
    edges = np.asarray(transfers, dtype=bool)
    n = edges.shape[0]
    reach = edges | np.eye(n, dtype=bool)
    for k in range(n):
        reach |= reach[:, [k]] & reach[[k], :]
    classes: list[tuple[int, ...]] = []
    seen: set[int] = set()
    for g in range(n):
        if g in seen:
            continue
        members = tuple(h for h in range(n) if reach[g, h] and reach[h, g])
        seen.update(members)
        if not any(reach[m, h] for m in members for h in range(n) if h not in members):
            classes.append(members)
    return tuple(classes)


def _net_removal_is_zero(mixture: "Mixture", group: int) -> bool:
    scattering = np.asarray(mixture.SigS[0].toarray(), dtype=float)[group]
    n2n = np.asarray(mixture.Sig2[0].toarray(), dtype=float)[group]
    return math.fsum([float(mixture.SigT[group]), *(-scattering), *(-2.0 * n2n)]) == 0.0


def closed_faces(binding: "Optional[Binding]") -> bool:
    r"""Whether every face conserves: a deck, or a condition with return fraction exactly 1. ``None`` is an infinite medium."""
    if binding is None:
        return True
    return all(face in binding.decks or binding.condition_of(face).return_fraction == 1.0 for face in binding.faces)


def absorption_free_closed_classes(
    mixtures: "Iterable[Mixture]", binding: "Optional[Binding]",
) -> tuple[tuple[int, ...], ...]:
    r"""The final classes that never lose a neutron on this domain; ``()`` when the loss is structurally regular.

    ``mixtures`` are the PRESENT materials (those a region of the domain
    carries); ``binding`` the domain's faces, or ``None`` for an infinite
    medium. The classes are the witness a refusal or a record names.
    """
    present = tuple(mixtures)
    if not present or not closed_faces(binding):
        return ()
    groups = int(present[0].ng)
    transfers = np.zeros((groups, groups), dtype=bool)
    for mixture in present:
        transfers |= np.asarray(mixture.SigS[0].toarray()) != 0.0
        transfers |= np.asarray(mixture.Sig2[0].toarray()) != 0.0
    np.fill_diagonal(transfers, False)
    return tuple(
        members for members in final_classes(transfers)
        if all(_net_removal_is_zero(mixture, g) for mixture in present for g in members)
    )

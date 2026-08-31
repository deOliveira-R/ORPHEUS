r"""The producer-side isotropic per-ordinate combine — ONE home.

The per-ordinate source of an isotropic scalar emission is
:math:`q_n = q^{\rm iso}/W` (:math:`W = \int_{S^2} d\Omega` — the angular
measure's total weight): the :math:`\ell=0` reconstruction of the frame,
realized on the reaction-rate fast path (no moment tensor — a
load-bearing PERF property of the SI hot loop). Every producer that
emits a per-ordinate :class:`~orpheus.transport.source_sinks.AngularSourceSink`
from an isotropic scalar source routes HERE — the :math:`/W` convention
lives once (its normalisation chain:
``docs/theory/methods/sn/slab_multigroup.rst``):

* :class:`~orpheus.transport.operators.scattering.ScatteringOperator`'s
  P0 half (combined with its :math:`\ell\ge 1` anisotropic part);
* :class:`~orpheus.transport.operators.n2n.N2NOperator`'s whole action
  (the §14.1 extraction — ORPHEUS MODELS :math:`(n,2n)` emission as
  isotropic, a :math:`P_0` truncation of the evaluated data rather than
  a property of the reaction; ``docs/theory/methods/sn/adjoint.rst``
  §sn-n2n-p0-truncation, issue #426).

Shared as a free function (the CS4c §14.1 landing): the two consumers
are different OPERATORS of one composite algebra, and the combine is the
piece of arithmetic they must never spell twice. When a third isotropic
lifted channel appears, this is the seed of the generic lift operator
(defer-until-2 — recorded, not minted).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

from orpheus.transport.source_sinks import AngularSourceSink

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.source_sinks import ScalarSourceSink

__all__ = ["assemble_per_ordinate_isotropic"]


def assemble_per_ordinate_isotropic(
    iso: "ScalarSourceSink",
    aniso: "AngularSourceSink | None",
    angular_space: "FunctionSpace",
    total_weight: float,
) -> "AngularSourceSink":
    r"""Combine an isotropic scalar source with an (optional) per-ordinate
    anisotropic part: :math:`(\text{iso}/W) + \text{aniso}`.

    Parameters
    ----------
    iso : ScalarSourceSink
        The scalar-magnitude isotropic source (P0 in-scatter, (n,2n)
        emission, …) riding the driving flux's own space.
    aniso : AngularSourceSink or None
        Per-ordinate :math:`\ell\ge 1` contribution ALREADY in
        per-ordinate magnitude, or ``None`` for a purely isotropic
        producer.
    angular_space : FunctionSpace
        The per-ordinate target space — sizes the zero accumulator when
        ``aniso`` is ``None`` (the caller holding the pose supplies it).
    total_weight : float
        :math:`W` — the binding measure's total angular weight.
    """
    aniso_part = (
        aniso if aniso is not None else AngularSourceSink.zeros(angular_space)
    )
    # The containment dunder's cross-class arm returns the LARGER
    # (angular) class — the #288 principled LSP exception the static
    # union cannot carry.
    return cast(
        "AngularSourceSink", (iso / total_weight) + aniso_part,
    )

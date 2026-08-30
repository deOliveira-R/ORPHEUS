r"""The one energy-conformity guard the CS4a bindings share.

A bound operator's data has an energy extent (``ng``); the space it is
bound to MAY carry an :class:`~orpheus.numerics.axis.EnergyAxis` stating
the posing's group count. When both are present they must agree — a 2g
coefficient bound to a 4g posing is ill-posed at construction, not at
apply time (Pattern 4: refuse where the illegal state is FORMED).

**Measured reach (vv#28 — the fraction is the contract, 2026-08-20):**
live on **4 of 13 production bindings** = 192 of 1022 shipped
constructions (18.8 %); **inert on 830** — 578 axes-less (the SN /
diffusion composite spaces carry ``axes=None``, so there is nothing to
key on) and 252 space-less — the two ``isotropic_kernel`` constructions
(``scattering.py``), which pass no ``space`` by signature. C's
Optional retired at CS4c step 2 (2026-08-30: mandatory kw-only ends on
the BoundOperator base, per-END admission through
``_assert_energy_extent_both_ends``); C already contributed 0 of the
252 (CS4a-R CEN-3 — every shipped binding passed a space). The axis-keyed strengthening for composites arrives
with CS2, when the composite spaces gain axes. Keying on the space's
*shape* instead was measured UNRUNNABLE: ``space.shape[0]`` is not
``ng`` on any composite (flat ``(64,)``; interior ordinate-first), and a
shape-keyed guard destroyed 250 of 845 battery rows when probed.

The refusal's message fragment is ``"energy extent"`` — deliberately
disjoint from the ``OperatorSum`` guards' ``"equal domains"`` /
``"A.domain == B.codomain"`` vocabulary, so a test pinning this guard
can never be intercepted by those.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.numerics.axis import EnergyAxis

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = ["assert_energy_extent_conforms"]


def assert_energy_extent_conforms(
    space: "FunctionSpace | None", ng: int, *, operator: str
) -> None:
    """Refuse a binding whose space states a DIFFERENT group count than its data.

    Fires only when there is something to key on (module docstring): a
    ``None`` space and an axes-less space both pass — presence is the
    constructor's own signature contract, and the axes-less inertness is
    declared, not accidental.
    """
    if space is None:
        return
    axes = space.axes
    if axes is None:
        return
    for axis in axes:
        if isinstance(axis, EnergyAxis):
            (extent,) = axis.shape
            if extent != int(ng):
                raise ValueError(
                    f"{operator}: energy extent mismatch — the space's "
                    f"EnergyAxis has {extent} groups while the bound data "
                    f"has ng={int(ng)}. The posing and the cross-section "
                    f"data must state the SAME group count; condense or "
                    f"re-pose before binding."
                )
            return

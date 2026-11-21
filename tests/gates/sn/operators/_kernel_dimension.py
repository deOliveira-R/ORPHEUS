r"""The counting-law oracle ``predicted_kernel_dimension`` (moved off production).

``dim ker A`` of the SN loss operator from the combinatorial counting law of
``orpheus.sn.operators.loss_kernel_gauge`` §3 — **without building any vector**. It is the
structurally independent check on the gauge: this walks the combinatorics and returns an
integer, while the gauge walks the pair generators and takes an SVD. Agreement between them
is evidence that neither has a bookkeeping error, and it is the only available oracle at
production sizes, where the dense SVD is not viable.

No production path read it; it lived beside the gauge it checks and moved to the gate
directory of the package it left. It reads three production names, the public
``gauge_freedom`` and the module-private ``_reflection_orbits`` and ``_anova_dimension``
(the separable-equation dimension of §3, which production keeps).
"""
from __future__ import annotations

from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np

from orpheus.sn.operators.loss_kernel_gauge import (
    _anova_dimension,
    _reflection_orbits,
    gauge_freedom,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from orpheus.sn.problem import SNProblem


def predicted_kernel_dimension(problem: "SNProblem") -> int:
    r"""``dim ker A`` from the counting law, counting the **R** component only.

    Matching the gauge's scope: the tangential component ``T`` contributes one dimension per
    tangential trace DOF and is not gauged.
    """
    if not gauge_freedom(problem).present:
        return 0
    spatial = tuple(int(n) for n in problem.spatial_shape)
    total = 0
    for _orbit, active in _reflection_orbits(problem):
        if len(active) < 2:
            continue
        per_orbit = 0
        for size in range(2, len(active) + 1):
            for subset in combinations(active, size):
                spectators = int(np.prod(
                    [spatial[c] for c in range(problem.ndim)
                     if c not in subset] or [1]
                ))
                per_orbit += _anova_dimension(
                    tuple(spatial[c] for c in subset)
                ) * spectators
        total += per_orbit
    return total * int(problem.ng)

r"""The SN hub's predicted kernel, typed non-Optional (the narrowing site of the retired ``SNProblem.loss_kernel_gauge``).

``OperatorPencil.kernel`` is ``PredictedKernel | None`` because a fixture or an opaque pair carries none. An SN hub's pencil
always carries one (``SNProblem.pencil`` builds it from ``predicted_kernel``), and the retired property was the single site that
narrowed it. This is that site on the test side: it reads through the door the posing layer already uses
(``orpheus.numerics.posing._predicted``), so a consumer reads ``hub_kernel(problem).projector`` with no ``assert`` of its own.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.numerics.posing import _predicted

if TYPE_CHECKING:  # pragma: no cover - typing only
    from orpheus.numerics.fredholm import PredictedKernel
    from orpheus.sn.problem import SNProblem


def hub_kernel(problem: "SNProblem") -> "PredictedKernel":
    r"""The predicted kernel of ``problem.pencil`` (its gauge is ``.projector``); refuses when the pencil carries none."""
    return _predicted(problem.pencil, "hub_kernel")

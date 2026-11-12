r"""#707 true ends: the one adapter the Step A gates reach the transport pairing through (a root helper, not judged by
the test-layer gate).

The user's ruling (2026-11-12): the pairing is a first-class object that OWNS the adjoint.
``DualPairing(V, V*)`` is the transport pairing (the angular metric ⊗ the canonical spatial contraction);
``pairing(ψ, q)`` evaluates it; ``pairing.adjoint(A)`` returns the ``A†`` with ``⟨ψ†, Aψ⟩ = ⟨A†ψ†, ψ⟩``. ``.dual()``
stays the metric-free arrow and ``.H`` refuses on V → V*.

WHERE the pairing is obtained from (a space verb such as ``V.pairing``, or a problem attribute) is not ruled. That one
spelling lives in :func:`pairing_of` and nowhere else, so the carve chooses it by editing this function only. Until
the carve, every call raises, and every gate that reaches the pairing is red for that reason.
"""

from __future__ import annotations

from typing import Any

#: The fragment the ``.H`` refusal on a V → V* operator must carry: it names the replacement (Pattern 4 corollary 3).
H_REFUSAL_FRAGMENT = "DualPairing"


def pairing_of(state_space: Any, data_space: Any) -> Any:
    """The transport pairing between a state space V and its data side V* (the carve chooses the spelling)."""
    from orpheus.numerics.pairing import DualPairing  # type: ignore[import-not-found]  # the carve's module

    return DualPairing(state_space, data_space)

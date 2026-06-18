r"""Shared test surrogate for the Morel--Montry weighted-diamond constants.

Issue #236 Phase 2 B2 (Fix 3 — dedup) — ``DiamondDifference.residual``
reads the angular-closure constants off ``CellVisit.c_in`` / ``c_out``,
sourced in PRODUCTION from ``SNMesh.pole_angular_closure`` (via
``SNMesh._make_cell_visit`` and the closure's ``c_{in,out}_per_ordinate``
accessors).

The diamond / cell-balance round-trip tests build a ``CellVisit`` from a
bare ``ReducedStreamingOperator`` or a synthetic ``StreamingTerms`` with
NO ``SNMesh`` / closure in scope, so they must stamp the fixture visit
with the SAME constants the closure would.  This surrogate is the ONE
hand-transcribed recompute of that formula:

.. math::

   c_{\rm out} = \alpha_{\rm out}/\tau, \qquad
   c_{\rm in}  = (1-\tau)/\tau\,\alpha_{\rm out} + \alpha_{\rm in}.

It is deliberately HAND-TRANSCRIBED and does NOT import the production
closure: it is the structurally-independent reference (vv L11) that keeps
the round-trip honest.  ``update`` reads the constants via
``cell_balance_terms`` (the ``DD.update`` solve path), the rewired
``residual`` reads them off the visit, and they must agree.  Slab carries
neutral ``α = 0`` / ``τ = 1`` so both are ``0.0``.

Before this module the recompute existed as two byte-identical private
copies (``_c_from_streaming_terms`` in ``test_diamond.py`` and ``_visit_c``
in ``test_cell_balance_for_streaming.py``); the production-stamp catcher
``test_cell_visit_c_stamp.py`` is a third consumer.  Unify-after-two-
instances (coding-elegance Pattern 2): ONE shared helper, all consumers
import it.
"""
from __future__ import annotations

from orpheus.geometry.reduced_operator import StreamingTerms


def c_from_streaming_terms(st: StreamingTerms) -> tuple[float, float]:
    r"""Weighted-diamond ``(c_in, c_out)`` from a per-cell ``StreamingTerms``.

    Hand-transcribed independent reference (vv L11 — NOT the production
    closure): ``c_out = α_out/τ``, ``c_in = (1-τ)/τ·α_out + α_in``.
    """
    tau, a_out, a_in = st.tau_mm, st.alpha_out, st.alpha_in
    if tau is None or a_out is None or a_in is None:
        raise ValueError(
            "c_from_streaming_terms requires a fully-populated StreamingTerms "
            "(tau_mm / alpha_out / alpha_in are None — an unbound packet)."
        )
    c_out = a_out / tau
    c_in = (1.0 - tau) / tau * a_out + a_in
    return c_in, c_out

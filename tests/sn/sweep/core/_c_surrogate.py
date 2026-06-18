r"""Shared test surrogate for the Morel--Montry weighted-diamond constants.

Issue #236 Phase 2 B2 (Fix 3 — dedup) — ``DiamondDifference.residual``
reads the angular-closure constants off ``CellVisit.c_in`` / ``c_out``
and the angular weight off ``CellVisit.tau``, all sourced in PRODUCTION
from ``SNMesh.pole_angular_closure`` (via ``SNMesh._make_cell_visit`` and
the closure's ``c_{in,out}_per_ordinate`` / ``tau_per_ordinate``
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
the round-trip honest.

Issue #236 Phase 2 Step C (the geometry-τ retirement)
-----------------------------------------------------
Before Step C the surrogate read τ/α straight off a ``StreamingTerms``
packet (``st.tau_mm`` / ``st.alpha_in`` / ``st.alpha_out``).  Step C
RETIRES those three fields from ``StreamingTerms`` (and the geometry-side
τ producer that baked them), so the surrogate can no longer take a bare
``st``.  The independent τ now comes from the structurally-independent
``morel_montry_weights`` (``contamination.py`` — a DIFFERENT code path to
the SAME BMC-2010-Eq.43 weight; vv L11) WITH the production clamp applied
(spherical UNCLAMPED, cylindrical ``clip(τ_raw, ½, 1)`` — mirroring the
streaming factories), and α comes from the SURVIVING dome arrays on the
``ReducedStreamingOperator`` (``alpha_half`` for spheres,
``alpha_per_level`` for cylinders; α is NOT retired, only its
``StreamingTerms`` packing is).

Two entry points:

* :func:`c_from_constants` — the pure hand-transcribed formula on an
  explicit ``(tau, alpha_in, alpha_out)`` triple (used by tests that
  build synthetic packets with chosen τ/α, or that resolve the triple
  via :func:`mm_constants_for_ordinate`).
* :func:`mm_constants_for_ordinate` — resolves the ``(tau, alpha_in,
  alpha_out)`` triple for one (cell, ordinate) of a geometry operator,
  sourcing τ from the INDEPENDENT reference and α from the operator's
  surviving dome.  Slab (Cartesian) is the neutral element τ = 1, α = 0.

This is the surrogate the diamond / cell-balance fixtures stamp visits
with, and the production-stamp catcher ``test_cell_visit_c_stamp.py``
cross-checks against the live ``SNMesh._make_cell_visit`` stamp.
Unify-after-two-instances (coding-elegance Pattern 2): ONE shared helper,
all consumers import it.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.discrete.sn.contamination import morel_montry_weights
from orpheus.geometry import CoordSystem
from orpheus.geometry.reduced_operator import ReducedStreamingOperator


def c_from_constants(
    tau: float, alpha_in: float, alpha_out: float,
) -> tuple[float, float]:
    r"""Weighted-diamond ``(c_in, c_out)`` from an explicit M-M triple.

    Hand-transcribed independent reference (vv L11 — NOT the production
    closure): ``c_out = α_out/τ``, ``c_in = (1-τ)/τ·α_out + α_in``.
    """
    c_out = alpha_out / tau
    c_in = (1.0 - tau) / tau * alpha_out + alpha_in
    return c_in, c_out


def mm_constants_for_ordinate(
    op: ReducedStreamingOperator,
    cell_idx: int,
    direction_idx: int,
    mu_level_idx: int | None = None,
) -> tuple[float, float, float]:
    r"""Resolve the M-M ``(tau, alpha_in, alpha_out)`` for one ordinate.

    The independent τ comes from
    :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
    (a structurally-independent code path to the BMC-2010-Eq.43 weight;
    vv L11) WITH the production clamp — spheres UNCLAMPED, cylinders
    ``clip(τ_raw, ½, 1)``.  α comes from the operator's surviving dome
    (``alpha_half`` / ``alpha_per_level``).  Slab is the neutral element
    τ = 1, α = 0.

    Parameters mirror
    :meth:`ReducedStreamingOperator.streaming_terms`: ``direction_idx``
    is the global ordinate for slab / sphere, the within-level azimuthal
    index for cylinder (with ``mu_level_idx`` selecting the μ-level).

    ``cell_idx`` is unused (the M-M constants are cell-independent) but
    accepted so call sites pass the same ordinate triple they pass to
    ``streaming_terms``.
    """
    del cell_idx  # M-M constants are independent of the radial cell

    if op.coord is CoordSystem.CARTESIAN:
        # Slab: neutral identity closure (no half-angle dome).
        return 1.0, 0.0, 0.0

    quad = op._quadrature
    if quad is None:
        raise ValueError(
            "mm_constants_for_ordinate requires the operator to carry its "
            "quadrature (op._quadrature is None)."
        )

    if op.coord is CoordSystem.SPHERICAL:
        if op.alpha_half is None:
            raise ValueError("spherical operator missing alpha_half dome.")
        tau_ref = morel_montry_weights(quad, "spherical")  # (N,), unclamped
        tau = float(tau_ref[direction_idx])
        alpha_in = float(op.alpha_half[direction_idx])
        alpha_out = float(op.alpha_half[direction_idx + 1])
        return tau, alpha_in, alpha_out

    if op.coord is CoordSystem.CYLINDRICAL:
        if mu_level_idx is None:
            raise ValueError(
                "cylindrical mm_constants_for_ordinate requires mu_level_idx."
            )
        if op.alpha_per_level is None:
            raise ValueError("cylindrical operator missing alpha_per_level.")
        # Independent reference is UNCLAMPED; the production cylinder τ is
        # clamped to [½, 1] (structural ÷0 block at the most-inward
        # ordinate; see cylindrical_streaming / #229).
        tau_raw = morel_montry_weights(quad, "cylindrical")  # list[(M,)], raw
        tau = float(np.clip(tau_raw[mu_level_idx][direction_idx], 0.5, 1.0))
        alpha_lv = op.alpha_per_level[mu_level_idx]
        alpha_in = float(alpha_lv[direction_idx])
        alpha_out = float(alpha_lv[direction_idx + 1])
        return tau, alpha_in, alpha_out

    raise ValueError(f"Unknown coord system: {op.coord!r}")

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
``st``.  τ now comes straight from the production producer
``sn.angular.closure.morel_montry_tau_per_level`` (see the correction
below), and α comes from the SURVIVING dome arrays on the
``ReducedStreamingOperator`` (``alpha_half`` for spheres,
``alpha_per_level`` for cylinders; α is NOT retired, only its
``StreamingTerms`` packing is).

⛔ **WHAT THIS SURROGATE IS AND IS NOT INDEPENDENT OF (corrected
2026-08-11, Q5.6.4).** Two claims that used to sit here are now FALSE, and
are recorded rather than quietly dropped:

1. *"WITH the production clamp applied (… cylindrical ``clip(τ_raw, ½,
   1)``)"* — the ``[½, 1]`` absorber is **RETIRED**.  Worse, this file had
   *re-implemented* it (``np.clip(tau_raw, 0.5, 1.0)``) inside the
   supposedly independent oracle, directly under a comment claiming the
   reference was UNCLAMPED — a Pattern-2 twin of the very thing it stood
   in judgement over, which no symbol grep could find (only the CONCEPT
   grep did).  Deleted.
2. *"the structurally-independent ``morel_montry_weights`` — a DIFFERENT
   code path to the SAME weight (vv L11)"* — it no longer is.  That
   function delegated to the production producer, deliberately, so
   that a "reference" can never drift into a second definition of the
   angular cell — which is exactly how its cylinder arm went wrong.
   ⛔ **It was RETIRED 2026-08-12** and this file now calls
   ``sn.angular.closure.morel_montry_tau_per_level`` directly. The
   wrapper lived in ``derivations/`` (L0), and its body WAS an
   ``orpheus.sn`` import — an edge ``tests/test_layer_imports.py``
   forbids. Nothing is lost: the delegation it performed is exactly what
   the direct call performs, and naming production at the call site makes
   the tautology recorded below visible HERE rather than one frame down.

⟹ **the τ leg of any gate built on this surrogate is now TAUTOLOGICAL**:
production τ compared with itself through a wrapper.  Do not read it as
coverage.  What this surrogate still tests honestly:

* :func:`c_from_constants` — the hand-transcribed
  ``c_out = α_out/τ``, ``c_in = (1-τ)/τ·α_out + α_in``, which shares no
  code with the closure and IS genuinely independent;
* the (cell, ordinate) → (level, within-level) **gather/permutation**,
  which is what ``test_cell_visit_c_stamp.py`` exists to catch.

τ itself is pinned independently elsewhere —
``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` (the
analytic closed form on the cylinder, a hand-written cumulative-weight
expression on the sphere).  That is the gate to consult for τ, not this
one.  ⚠ And note the c-stamp cylinder fixture is ``folded_product(2, 4)``,
i.e. M = 2 per level — `[M]` the ω-midpoint and retired chord partitions
are **bit-identical** at M = 2, so that fixture could not see a partition
change even if this leg were independent.

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

from orpheus.geometry import CoordSystem
from orpheus.geometry.reduced_operator import ReducedStreamingOperator
from orpheus.sn.angular.closure import morel_montry_tau_per_level


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

    τ comes from the production producer
    :func:`~orpheus.sn.angular.closure.morel_montry_tau_per_level`,
    unclamped on both arms.  α comes from the operator's surviving dome
    (``alpha_half`` / ``alpha_per_level``).  Slab is the neutral element
    τ = 1, α = 0.

    ⛔ This paragraph carried THREE falsehoods until 2026-08-12, each
    already refuted by this module's own header — recorded rather than
    quietly dropped, because a docstring that contradicts its own file is
    exactly what a reader trusts:

    1. it named
       ``angular_differencing.morel_montry_weights``, RETIRED that day
       (its body was an ``orpheus.sn`` import, which L0 may not make);
    2. it called that "a structurally-independent code path … vv L11" —
       header §2 says it is not, and the τ leg here is TAUTOLOGICAL;
    3. it claimed "WITH the production clamp — cylinders
       ``clip(τ_raw, ½, 1)``" — header §1 says the ``[½, 1]`` absorber is
       retired and this file's re-implementation of it was DELETED.

    Parameters mirror
    :meth:`ReducedStreamingOperator.streaming_terms`: ``direction_idx``
    is the global ordinate for slab / sphere, the within-level azimuthal
    index for cylinder (with ``mu_level_idx`` selecting the μ-level).

    ``cell_idx`` is unused (the M-M constants are cell-independent) but
    accepted so call sites pass the same ordinate triple they pass to
    ``streaming_terms``.
    """
    del cell_idx  # M-M constants are independent of the radial cell

    if op.mesh.coord is CoordSystem.CARTESIAN:
        # Slab: neutral identity closure (no half-angle dome).
        return 1.0, 0.0, 0.0

    # The measure comes from the ANGULAR factor, which is non-optional on
    # every chart -- so the former ``if quad is None: raise`` guard retired
    # with the ``_quadrature`` twin on 2026-08-26, exactly as the
    # ``op.alpha_half is None`` guard below it did at the un-weld.
    quad = op.angular.quadrature

    if op.mesh.coord is CoordSystem.SPHERICAL:
        # Production τ, named at the call site (2026-08-12): the L0 wrapper
        # ``angular_differencing.morel_montry_weights`` was retired because
        # its body WAS an ``orpheus.sn`` import, which L0 may not make. It
        # added nothing but a reshape, and the tautology it hid is stated
        # below. A test file may import ``sn`` freely.
        tau_ref = morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)[0]
        tau = float(tau_ref[direction_idx])
        # The dome is the shared ANGULAR factor (2026-08-26 un-weld).  The
        # former ``if op.alpha_half is None: raise`` guard retired with it:
        # every chart carries an angular factor now (Cartesian the neutral
        # one), so the None-ness branch was a coincidence proxy for
        # "spherical", which ``op.mesh.coord`` states directly one line
        # above.
        dome = op.angular.alpha_per_level[0]
        alpha_in = float(dome[direction_idx])
        alpha_out = float(dome[direction_idx + 1])
        return tau, alpha_in, alpha_out

    if op.mesh.coord is CoordSystem.CYLINDRICAL:
        if mu_level_idx is None:
            raise ValueError(
                "cylindrical mm_constants_for_ordinate requires mu_level_idx."
            )

        # ⛔ Q5.6.4 (2026-08-11): this line used to be
        #     tau = float(np.clip(tau_raw[...], 0.5, 1.0))
        # under a comment asserting "Independent reference is UNCLAMPED" —
        # i.e. the supposedly INDEPENDENT oracle re-implemented the very
        # [½,1] absorber it was standing in judgement over (a Pattern-2
        # twin that no symbol grep could find; only the CONCEPT grep did).
        # The absorber is retired, so the clip is deleted rather than
        # mirrored: this reference is now genuinely unclamped.
        tau_raw = morel_montry_tau_per_level(  # tuple[(M_p,), ...]
            quad, CoordSystem.CYLINDRICAL,
        )
        tau = float(tau_raw[mu_level_idx][direction_idx])
        alpha_lv = op.angular.alpha_per_level[mu_level_idx]
        alpha_in = float(alpha_lv[direction_idx])
        alpha_out = float(alpha_lv[direction_idx + 1])
        return tau, alpha_in, alpha_out

    raise ValueError(f"Unknown coord system: {op.mesh.coord!r}")

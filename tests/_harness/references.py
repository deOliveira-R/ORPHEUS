r"""Independent geometric references for mirror-gather gates.

:func:`mirror_partner_indices` is the REFERENCE side of gates that pin a
production mirror gather — realized specular operators, narrowed
permutations, hand-written specular closures. It derives the partner map
from nothing but the quadrature's direction cosines and plain numpy: no
call into :meth:`Quadrature.ordinate_permutation` /
:meth:`RigidMotion.preserves` / ``_orbit_closure``, so a defect in the
production derivation cannot cancel against it (``vv-principles`` L11
structural independence).

Until G6.3 §7d the retired ``Quadrature.reflection_index`` table played
this role: production and reference shared the partner DATUM while the
gather arithmetic stayed independent. Step 7 moved production onto
``ordinate_permutation``, so a reference that read the same source would
have demoted every one of those gates to production-vs-itself on the
partner map (the coding-standards rewire-demotion hazard). This helper
restores the stronger form: the reference now shares NOTHING with
production above the numpy line.

**Why explicit ``raise``, not ``assert``** — the same Mode-8 note as
``predicates.py``: this is a plain helper module pytest does not
AST-rewrite, so under the canonical ``python -O`` invocation a bare
``assert`` would be stripped to a no-op.
"""

from __future__ import annotations

import numpy as np

from orpheus.numerics.invariance import _embedded_nodes

#: Axis letter → direction-cosine column. Deliberately a local literal
#: rather than an import of the production convention home
#: (``face_layout.AXIS_NAMES``): the reference must stay readable as
#: self-contained numpy, and a convention drift between the two IS a
#: defect these gates should surface, not absorb.
_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def mirror_partner_indices(quad, axis: int | str) -> np.ndarray:
    r"""Partner-ordinate indices under the mirror about ``axis``, derived
    independently of production.

    Lifts the ordinates to ``(N, 3)`` orbit barycentres — direction
    cosines for a full-sphere rule, ``(μ, 0, 0)`` for a slab rule (its
    missing components are genuinely zero), and since R4 of #434
    (2026-09-03) the barycentre ``(x, 0, z)`` of a σ_y-folded rule's
    representative, which lies in the ball — negates the named component,
    and EXACT-matches each image against the node set. The
    window (``1e-12``) only absorbs trig-roundoff of the ordinate
    constructors — a sign flip itself is IEEE-exact.

    Raises when any image misses the node set or the match is not
    injective: a nearest-neighbour ``argmin`` without those two checks
    computes a RELATION, not a permutation (``vv-principles`` #14 /
    ERR-073) — the exact defect class ERR-074 catalogued on the retired
    partner table's pre-certification body.
    """
    idx = _AXIS_INDEX[axis] if isinstance(axis, str) else int(axis)
    # The canonical R^3 embedding, from its one home. This used to be
    # ``column_stack([quad.axis_cosines(a) for a in range(3)])``, which worked
    # only because the accessor silently padded a suppressed axis; phase 0.2
    # made that a refusal and surfaced the duplicate. A mirror is LINEAR, so it
    # commutes with the orbit mean — matching embedded barycentres is therefore
    # equivalent to matching orbits, and the 1-D case stays correct.
    #
    # ⚠ Independence note (``vv-principles`` #22), PER SUPPORT: production's
    # ``ordinate_permutation`` certifies through the orbit-space closure,
    # which embeds with this SAME function. So this reference no longer
    # cross-checks the EMBEDDING convention — it cross-checks the PERMUTATION,
    # which is its subject. For the nine unfolded rules the embedding is the
    # identity or a zero-pad and the sharing ends there; for a FOLD (since R4
    # of #434, 2026-09-03) the embedding routes through the entry's chart and
    # lift, the same pair the closure reads — a wider sharing, stated so that
    # a fold row is not read as an independent check of the lift. On a fold
    # the y-mirror's answer is the identity permutation, and it agrees with
    # ``DiscreteMeasure.permutation_under`` (gated in test_manifold.py's R4
    # rows); until R4 the two disagreed — this raised where the closure
    # answered.  ``_embedded_nodes`` lives in ``orpheus.numerics.invariance``
    # since R2 of #434.
    nodes = np.asarray(_embedded_nodes(quad.measure), dtype=float)
    mirrored = nodes.copy()
    mirrored[:, idx] *= -1.0
    distances = np.linalg.norm(mirrored[:, None, :] - nodes[None, :, :], axis=-1)
    partner = np.argmin(distances, axis=1)
    residual = float(np.max(distances[np.arange(len(nodes)), partner]))
    if residual > 1e-12:
        raise AssertionError(
            f"no mirror partner about axis {axis!r}: a node's mirror image "
            f"misses the node set by {residual:.3e} — the rule is not closed "
            f"under that mirror, so no reference partner map exists."
        )
    if len(set(partner.tolist())) != len(partner):
        raise AssertionError(
            f"the mirror match about axis {axis!r} is not injective — a "
            f"relation, not a permutation (ERR-073 class)."
        )
    return partner

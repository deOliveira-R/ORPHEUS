r"""The Lathrop α-defect measures the DEFECT, not the α normalization.

Written 2026-08-27 as the witness for a repair, per ``plan-authoring``
§6c (a fix lands with the case it catches).

**The defect this gate exists to catch.**
:func:`~orpheus.derivations.discrete.sn.angular_differencing.alpha_defect_beta`
returns Lathrop 2000 Eq. 25's :math:`\beta_{m+1/2} = \alpha_{m+1/2} -
(1 - \eta_{m+1/2}^2)`.  Lathrop writes Eq. 25 in the **Hébert**
convention :math:`\alpha^{H}`; :func:`alpha_dome` returns the ORPHEUS
one, :math:`\alpha^{O} = \alpha^{H}/2`.  Until 2026-08-27 the two sides
were subtracted directly, so the function returned the *normalization
difference* rather than the defect — `[M]` it converged to **0.5**
(``0.4787 → 0.4942 → 0.4985 → 0.4996 → 0.4999`` at N = 4/8/16/32/64)
and its own docstring's "zero iff τ ≡ ½" was unsatisfiable by any input.

It shipped undetected because the function had **zero consumers and
zero gates** tree-wide — a prior QA review had recorded it as
"unguarded, UNGATED" and nothing followed.  This module is that gate.

**Instrument type: CONSTRAINT** (``lessons.md`` L51).  It answers "are
the two sides of Eq. 25 in one normalization?" — pass/fail.  It is NOT
a ranker: it cannot order two angular schemes by accuracy, and must
never be cited in an accuracy argument.

**Configuration** (``plan-authoring`` §4 — a number without its fixture
is not reusable): sphere, Gauss–Legendre with :math:`\sum w = 2`, and
the **cumulative-weight** angular cell partition (BMC Eq. 12, the
sphere's own convention).  The edge set is built here from a literal
``-1 + cumsum(w)`` rather than read from production, so the two sides
of the comparison are independently produced.

⚠ **What this does NOT cover.** A second, independent hazard is
recorded in the function's docstring and is *not* gated here because it
was never measured: ``edges`` must be a partition in the **cosine**,
and nothing validates the argument's units.  Feeding radians yields
silent garbage.  Do not read this module's green as covering that.

**Mutation-verified**: setting ``_HEBERT_PER_ORPHEUS_ALPHA`` back to
``1.0`` (the pre-repair behaviour) reddens every row below.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.discrete.sn import angular_differencing as ad
from orpheus.numerics.quadrature import Quadrature

pytestmark = pytest.mark.foundation

# Ladder deliberately breaks its own arithmetic pattern (vv-principles
# #13): two odd orders and two non-powers-of-two, not 4/8/16/32.
_ORDERS = [4, 6, 10, 16, 24, 32]


def _sphere_case(n: int):
    """Ascending-μ GL rule plus its cumulative-weight edge partition."""
    quad = Quadrature.gauss_legendre(n)
    order = np.argsort(quad.mu_x)
    w = np.asarray(quad.weights, dtype=float)[order]
    edges = np.concatenate(([-1.0], -1.0 + np.cumsum(w)))
    return quad, [edges]


class TestDefectIsTheDefect:
    r"""β must vanish under refinement, not tend to the normalization."""

    @pytest.mark.parametrize("n", _ORDERS)
    def test_defect_is_small_not_one_half(self, n: int) -> None:
        """The whole discriminator, in one assertion.

        The pre-repair expression reads ≈0.48–0.50 at EVERY order; the
        repaired one is below 0.05 by N=4 and falls from there.  A
        threshold of 0.25 sits an order of magnitude clear of the
        repaired values and a factor of two clear of the broken one, so
        it cannot be satisfied by the normalization difference at any n.
        """
        quad, edges = _sphere_case(n)
        beta = ad.alpha_defect_beta(quad, "spherical", edges=edges)[0]
        worst = float(np.abs(beta).max())
        assert worst < 0.25, (
            f"S{n}: max|beta| = {worst:.4f}. A value near 0.5 means the two "
            "sides of Lathrop Eq. 25 are in different alpha normalizations "
            "(alpha^H vs alpha^O) — see this module's docstring."
        )

    def test_defect_converges_at_second_order(self) -> None:
        r"""β → 0 as :math:`O(N^{-2})`, which is the claim Eq. 25 makes.

        Convergence — not merely smallness — is what distinguishes the
        defect from any constant offset.  The pre-repair form is FLAT
        (it rises slightly toward 0.5), so this leg reddens on shape
        alone even if a threshold were relaxed.
        """
        worst = []
        for n in (8, 16, 32, 64):
            quad, edges = _sphere_case(n)
            beta = ad.alpha_defect_beta(quad, "spherical", edges=edges)[0]
            worst.append(float(np.abs(beta).max()))

        for coarse, fine in zip(worst[:-1], worst[1:]):
            assert fine < coarse, (
                f"beta must FALL under refinement; got {worst}. A flat or "
                "rising sequence is the normalization signature."
            )
        # Doubling N must buy at least 3x (second order is 4x; 3x leaves
        # room for the edge partition's own first-order content).
        ratios = [c / f for c, f in zip(worst[:-1], worst[1:])]
        assert min(ratios) > 3.0, f"ratios {ratios} are not second-order-like"

    def test_the_repair_constant_is_load_bearing(self) -> None:
        """Positive control: reverting the lift reproduces the defect.

        vv-principles #17 — a gate whose failure mode cannot be
        demonstrated is not evidence.  This drives the mutation
        IN-PROCESS via monkeypatching the module global, never by
        editing the file (``process-discipline``: a mutation left on
        disk by a killed run is committed silently).
        """
        quad, edges = _sphere_case(16)
        healthy = float(np.abs(
            ad.alpha_defect_beta(quad, "spherical", edges=edges)[0]
        ).max())

        saved = ad._HEBERT_PER_ORPHEUS_ALPHA
        try:
            ad._HEBERT_PER_ORPHEUS_ALPHA = 1.0  # the pre-repair behaviour
            broken = float(np.abs(
                ad.alpha_defect_beta(quad, "spherical", edges=edges)[0]
            ).max())
        finally:
            ad._HEBERT_PER_ORPHEUS_ALPHA = saved

        assert healthy < 0.05, f"repaired reading unexpectedly large: {healthy}"
        assert broken > 0.4, f"mutation did not bite: {broken}"
        assert ad._HEBERT_PER_ORPHEUS_ALPHA == 2.0, "restore failed"

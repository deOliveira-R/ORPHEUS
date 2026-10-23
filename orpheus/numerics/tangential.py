r"""The tangential tolerance for the direction cosine :math:`\Omega\cdot\hat n`.

A leaf of ``numerics`` (no intra-package imports) so that both the angular
trace spaces and the quadrature certification
(:func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing`) read
it without ``quadrature`` importing ``spaces`` (the module-cycle gate counts
every import kind). Re-exported by
:mod:`orpheus.numerics.spaces.angular_trace_space` under the same name.
"""

import numpy as np

__all__ = ["TANGENTIAL_EPS"]

# Tangential tolerance for the unit-vector projection ``Ω · n``. A
# safety factor (×4) over the IEEE-754 dot-product round-off bound for a
# 3-component unit-vector inner product. Empirically bit-identical to
# the legacy ``1e-15`` (operator) and ``1e-12`` (realizer) tolerances —
# and since #325's node repoint a provably-inert defensive guard (every
# shipped rule's tangential cosines are exactly 0.0, so the band
# ``(0, eps]`` is empty); see the module docstring +
# :func:`test_eps_sits_in_the_round_off_to_genuine_gap`.
# Public: the ONE tangential threshold every trace-classification
# consumer shares (the selectors below; the reflective law's
# inflow→outflow invariant, ERR-045) — a second locally-minted epsilon
# would be a twin source of truth for "what counts as grazing".
TANGENTIAL_EPS: float = 4.0 * np.finfo(np.float64).eps

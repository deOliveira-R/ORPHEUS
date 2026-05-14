r"""First-order linear-recurrence scan primitive — Blelloch 1990 §1.5.

Issue #196 Phase G Step 2.5b.  This module exposes a single free
function, :func:`ordinate_scan`, which is the canonical closed-form
solution to the affine recurrence

.. math::

    \psi[i+1] \;=\; a[i]\,\psi[i] + b[i],\qquad \psi[0] = \psi_0.

The recurrence is the per-ordinate spatial-sweep update of the
discrete S\ :sub:`N` operator: forward substitution on the block-
triangular (streaming + collision) operator.  Each per-cell update
is one element of the **pair-monoid** :math:`(a, b)` with composition

.. math::

    (\alpha_1, \beta_1) \;\oplus\; (\alpha_2, \beta_2)
        \;=\; (\alpha_1 \cdot \alpha_2,\; \alpha_2\,\beta_1 + \beta_2),

equivalent to the 2×2 lower-triangular affine matrix
:math:`M = [[a, 0], [b, 1]]` acting on the augmented vector
:math:`(\psi, 1)^\top`.  Matrix-multiplication composes them; the
scan of nx cells is the prefix-product of these matrices.

Because the 2×2 has lower-triangular structure, the matrix scan
factors into **two scalar scans**:

.. math::

    \psi[n] \;=\; \left(\prod_{i=0}^{n-1} a[i]\right)
                  \left(\psi_0 + \sum_{i=0}^{n-1}
                        \frac{b[i]}{\prod_{j=0}^{i} a[j]}\right).

This is the Blelloch §1.5 first-order linear-recurrence form.  In
numpy: ``cumprod_a · (psi_0 + cumsum(b / cumprod_a))`` — three numpy
ops.  No Python loop over cells.

The function is a **free function**, not a method on a class or on
the :class:`CellUpdate` Protocol.  The 1-D sweep
(:func:`~orpheus.sn.sweep._sweep_1d_unified`) CONSUMES ``ordinate_scan``
on the cache's ``a_attenuation`` + the per-iteration ``b`` vector
(Step 2.5c).

Pair-monoid contract
====================

The closed-form decomposition is justified by associativity:

* :math:`(M_1 \oplus M_2) \oplus M_3 = M_1 \oplus (M_2 \oplus M_3)`
  (pair-monoid associativity — the theorem).
* :math:`(1, 0) \oplus M = M = M \oplus (1, 0)` (identity element).
* Brent's theorem: associative scans admit an
  :math:`O(N/\log N)`-work parallel decomposition; the sequential
  numpy implementation is one valid backend.

The algebraic-theorem test suite at
``tests/sn/spatial/test_ordinate_scan.py`` pins these invariants —
the theorems justify the implementation, not the other way around.

References
==========

* Blelloch, G. E. (1990).  *Prefix Sums and Their Applications*.
  CMU-CS-90-190 Technical Report.  §1.5 First-Order Linear
  Recurrences.
* Sengupta, S., Harris, M., Zhang, Y., & Owens, J. D. (2007).
  *Scan Primitives for GPU Computing*.  Graphics Hardware 2007,
  pp. 97–106.  Sec. 3 segmented scan (background; not used here).
* Brent, R. P. (1974).  *The Parallel Evaluation of General
  Arithmetic Expressions*.  J. ACM 21(2), 201–206.  Brent's
  theorem — work-efficient associative-scan reduction.
"""

from __future__ import annotations

import numpy as np


def ordinate_scan(
    a: np.ndarray,
    b: np.ndarray,
    psi_0: np.ndarray,
) -> np.ndarray:
    r"""First-order linear-recurrence scan (Blelloch 1990 §1.5).

    Returns the per-cell array :math:`\psi[1\ldots n_x]` satisfying

    .. math::

        \psi[i+1] \;=\; a[i]\,\psi[i] + b[i],\qquad \psi[0] = \psi_0,

    using the closed-form decomposition

    .. math::

        \psi[n] \;=\; \left(\prod_{i=0}^{n-1} a[i]\right)
                      \left(\psi_0 \;+\; \sum_{i=0}^{n-1}
                            \frac{b[i]}{\prod_{j=0}^{i} a[j]}\right).

    Equivalent (under associativity) to the pair-monoid
    composition :math:`(a_0, b_0) \oplus (a_1, b_1) \oplus \cdots`
    with :math:`(\alpha_1, \beta_1) \oplus (\alpha_2, \beta_2) =
    (\alpha_1 \alpha_2, \alpha_2 \beta_1 + \beta_2)`.

    Parameters
    ----------
    a : ndarray, shape ``(nx, ...)``.
        Per-cell multiplier sequence.  Leading axis is the chain
        (cell-ordering) axis.
    b : ndarray, shape ``(nx, ...)``.
        Per-cell additive sequence.  Same shape as ``a``.
    psi_0 : ndarray, shape ``(...)`` or scalar.
        Initial chain state (e.g. boundary-face inflow).  Trailing
        shape must broadcast against ``a[0]``.

    Returns
    -------
    psi : ndarray, shape ``(nx, ...)``.
        Per-cell chain output.  ``psi[i]`` is the state AFTER cell
        ``i`` (i.e. :math:`\psi[i+1]` in the recurrence notation
        above): the spatial-face flux on the downstream side of
        cell ``i``.

    Notes
    -----
    Numerical regime.  The Blelloch form requires ``cumprod_a`` to
    stay finite and bounded away from zero.  Diamond-Difference SN
    sweeps produce ``a[i] = 2|μ|·A_down/denom − 1`` with ``|a| < 1``
    in well-resolved regimes — both the cumprod and the
    ``b / cumprod_a`` quotient stay in IEEE-754's well-conditioned
    band.  For ``a → 0`` or ``a → ∞`` regimes outside DD's normal
    operating envelope, consult the test catalog at
    ``tests/sn/spatial/test_ordinate_scan.py::test_ordinate_scan_small_attenuation``
    for the documented regime limits.
    """
    cumprod_a = np.cumprod(a, axis=0)
    return cumprod_a * (psi_0 + np.cumsum(b / cumprod_a, axis=0))


__all__ = ["ordinate_scan"]

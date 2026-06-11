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
(:func:`~orpheus.sn.loss_representation._sweep_1d_unified`) CONSUMES ``ordinate_scan``
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
``tests/sn/sweep/core/test_ordinate_scan.py`` pins these invariants —
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
    Numerical regime.  Two backends compute the SAME recurrence; the
    dispatch is by **whether the Blelloch closed form is finite**:

    * **Closed form finite** (the common SN sweep path).  The Blelloch
      §1.5 form ``cumprod_a · (psi_0 + cumsum(b / cumprod_a))`` is
      returned unchanged — three numpy ops, bit-identical to every
      prior release.  Diamond-Difference produces ``a[i] = 2|μ|·A_down
      /denom − 1`` with ``|a| < 1`` in well-resolved regimes, so both
      the cumprod and the ``b / cumprod_a`` quotient stay in IEEE-754's
      well-conditioned band.

    * **Closed form non-finite** — the division-free **pair-monoid**
      backend (:func:`_pair_monoid_scan`) takes over.  This arises in
      two distinct regimes, BOTH of which the previous proxy guard
      ``cumprod_a[-1] == 0`` mishandled:

      - **Exact reset** (``a[i] = 0`` exactly — e.g. the cylindrical
        pole cell whose inner radial face area vanishes at the
        ``2|μ|·A_total = ΔA_w·c_out + Σ_t·V`` resonance, Issue #209).
        At a reset the recurrence *forgets its history*
        (``psi[i+1] = b[i]``, a chain restart); the Blelloch division
        ``b / cumprod_a`` is ``b/0 = inf`` from the reset onward.

      - **Gradual underflow into the denormal band** (Issue #222,
        ERR-057).  For a long, optically-thick chain (``|a| < 1`` to
        ``nx ≳ 300`` of near-constant attenuation) the cumulative
        product decays into the IEEE-754 denormals
        (``cumprod_a ~ 1e-312``) WITHOUT hitting an exact zero, so
        ``b / cumprod_a`` overflows to ``inf`` and ``cumprod_a · inf =
        NaN``.  The old ``cumprod_a[-1] == 0`` guard tested a PROXY for
        the reset and silently MISSED this denormal band, leaking a
        NaN — a guard-predicate-incompleteness bug.

      The pair-monoid composes the affine cells ``(a, b)`` under
      ``(α₁, β₁) ⊕ (α₂, β₂) = (α₁α₂, α₂β₁ + β₂)`` via a Hillis–Steele
      log-depth scan — exact at ``a = 0`` by construction (no division
      anywhere), finite through the denormal band, and robust to
      consecutive and chain-leading resets.

    The dispatch tests the TRUE failure condition —
    ``np.all(np.isfinite(closed_form))`` over the already-computed
    closed form — NOT a proxy for one of its causes.  The check is a
    single O(N) reduction over an array the closed form already
    materialised: one extra pass, negligible against the
    cumprod/cumsum/multiply that produced it.
    """
    cumprod_a = np.cumprod(a, axis=0)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        closed_form = cumprod_a * (psi_0 + np.cumsum(b / cumprod_a, axis=0))
    if not np.all(np.isfinite(closed_form)):
        # Dispatch on the TRUE failure condition — a non-finite closed
        # form — NOT the cumprod_a[-1] == 0 proxy, which tested only one
        # of its causes (an exact reset a[i] = 0) and silently MISSED the
        # gradual denormal underflow (cumprod_a ~ 1e-312 ⇒ b / cumprod_a
        # overflows to inf ⇒ cumprod_a · inf = NaN), leaking a NaN
        # (ERR-057 / issue #222).  The division-free pair-monoid is exact
        # and finite in every one of these regimes.
        return _pair_monoid_scan(a, b, psi_0)
    return closed_form


def _pair_monoid_scan(
    a: np.ndarray,
    b: np.ndarray,
    psi_0: np.ndarray,
) -> np.ndarray:
    r"""Division-free affine-recurrence scan (Hillis–Steele 1986).

    Computes the same :math:`\psi[1\ldots n_x]` as :func:`ordinate_scan`
    but as a direct prefix scan of the **affine pair-monoid**

    .. math::

        (\alpha_1, \beta_1) \oplus (\alpha_2, \beta_2)
            \;=\; (\alpha_1\alpha_2,\; \alpha_2\,\beta_1 + \beta_2),

    rather than the Blelloch factored form.  Because the composition
    contains no reciprocal, the scan is exact when any
    :math:`\alpha = 0` (a chain reset ``psi[i+1] = b[i]``) and handles
    consecutive resets and a reset at the chain start with no
    special-casing — the reset *is* the monoid element ``(0, b)``,
    which annihilates everything to its left under ``⊕``.

    Implementation: an inclusive Hillis–Steele scan over the chain
    (leading) axis.  After :math:`\lceil\log_2 n_x\rceil` doubling
    passes, ``(alpha[i], beta[i])`` is the composition of cells
    ``0..i``; applying it to ``psi_0`` gives
    ``psi[i] = alpha[i]·psi_0 + beta[i]``.  Each pass is a single
    vectorised numpy slice-multiply-add over all trailing lanes, so
    the work is :math:`O(n_x \log n_x)` with no Python loop over
    cells.

    See Hillis & Steele (1986), *Data Parallel Algorithms*, CACM
    29(12), 1170–1183 — the inclusive parallel-prefix scan.
    """
    alpha = np.array(a, dtype=float, copy=True)
    beta = np.array(b, dtype=float, copy=True)
    nx = alpha.shape[0]
    step = 1
    while step < nx:
        # Compose cell i with the accumulated cell (i - step) sitting
        # to its left: (α_i, β_i) ⊕-after (α_{i-step}, β_{i-step}).
        alpha_left = alpha[:-step]
        beta_left = beta[:-step]
        beta[step:] = alpha[step:] * beta_left + beta[step:]
        alpha[step:] = alpha[step:] * alpha_left
        step *= 2
    return alpha * psi_0 + beta


__all__ = ["ordinate_scan"]


# ═══════════════════════════════════════════════════════════════════════
# The scan-march row primitives (relocated from sweep.py at S6.4(f) —
# the scan family's face recurrences live with the scan substrate)
# ═══════════════════════════════════════════════════════════════════════


def _x_scan_faces(
    alpha: np.ndarray,
    beta: np.ndarray,
    psi_x_in: np.ndarray,
    x_reverse: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""One x-line of the scan-march: solve the face recurrence, recover the faces.

    The shared face-scan primitive of BOTH directions (Pattern 2 — single
    source of the row-scan):

    * the forward **sweep** (:func:`_scanmarch_row`) passes the *solve*
      coefficients ``α = 2 s_x/D − 1``, ``β = 2 (Q + s_y ψ_{y,in})/D``;
    * the **matvec** twin
      (``ScanMarch.loss_action`` in ``orpheus.sn.loss_representation``, the
      apply-direction row-march S6.3 moved off the operator)
      passes the *apply* coefficients ``α = −1``, ``β = 2 ψ̄_probe`` — since ψ̄
      is KNOWN, the closure ``out_x = 2ψ̄ − in_x`` IS a first-order recurrence in
      the faces (a pure-reflection scan, ``|α| = 1`` so no underflow).

    Solves ``out_x(i) = α(i)·in_x(i) + β(i)`` (``in_x(i) = out_x(i−1)``,
    ``in_x(0) = psi_x_in``) by :func:`~orpheus.sn.spatial.scan.ordinate_scan`
    along the x-chain.

    Parameters
    ----------
    alpha, beta : np.ndarray
        ``(N_oct, ng, nx)`` scan coefficients in **mesh** x-order.
    psi_x_in : np.ndarray
        ``(N_oct, ng)`` domain x-inflow face value (the chain seed).
    x_reverse : bool
        ``True`` for a −x octant (scan high-x → low-x): the chain is reversed
        for the forward scan, the per-cell faces reversed back to mesh order.

    Returns
    -------
    in_x, out_x : np.ndarray
        ``(N_oct, ng, nx)`` upstream / downstream x-face flux per cell, mesh
        order.
    x_outflow : np.ndarray
        ``(N_oct, ng)`` domain x-outflow face value (the last swept cell's
        downstream x-face).
    """
    if x_reverse:
        alpha = alpha[..., ::-1]
        beta = beta[..., ::-1]
    # ordinate_scan scans axis 0 (the chain): move x to the front, scan, move
    # back. out_x_sweep[i] is the downstream x-face of swept-cell i.
    out_x_sweep = np.moveaxis(
        ordinate_scan(
            np.moveaxis(alpha, -1, 0), np.moveaxis(beta, -1, 0), psi_x_in,
        ),
        0, -1,
    )                                                  # (N_oct, ng, nx) sweep order
    # in_x[i] = out_x[i−1] (upstream face), in_x[0] = psi_x_in.
    in_x_sweep = np.concatenate(
        [psi_x_in[..., None], out_x_sweep[..., :-1]], axis=-1,
    )
    x_outflow = out_x_sweep[..., -1]                   # last swept cell = domain x-out
    if x_reverse:                                      # back to mesh order
        return in_x_sweep[..., ::-1], out_x_sweep[..., ::-1], x_outflow
    return in_x_sweep, out_x_sweep, x_outflow


def _scanmarch_row(
    alpha: np.ndarray,
    beta: np.ndarray,
    psi_x_in: np.ndarray,
    psi_y_in: np.ndarray,
    x_reverse: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""One y-row of the forward scan-march: x-scan, then recover :math:`\bar\psi`,
    the downstream y-face, and the domain x-outflow.

    Solves the in-row x-face recurrence via :func:`_x_scan_faces` (the *solve*
    coefficients ``α = 2 s_x/D − 1``, ``β = 2 (Q + s_y ψ_{y,in})/D``), then
    closes the diamond-difference cell with :math:`\bar\psi = \tfrac12
    (\mathrm{in}_x + \mathrm{out}_x)` and :math:`\mathrm{out}_y = 2\bar\psi -
    \psi_{y,\mathrm{in}}`.

    Returns ``(psi_avg, out_y, x_outflow)`` — the cell average and the
    downstream y-face (the next row's ``psi_y_in``) in mesh order, and the
    domain x-outflow face value.
    """
    in_x, out_x, x_outflow = _x_scan_faces(alpha, beta, psi_x_in, x_reverse)
    psi_avg = 0.5 * (in_x + out_x)
    out_y = 2.0 * psi_avg - psi_y_in
    return psi_avg, out_y, x_outflow

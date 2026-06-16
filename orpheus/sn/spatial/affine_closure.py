r"""Generic affine-scan closure operations — the COEFFICIENT MODEL (#158 Inc B).

A *consistent* affine spatial discretization (Diamond Difference, Linear
Discontinuous, Step, ...) is fully characterized, per cell, by **three
:math:`\Sigma_t`-epoch coefficients** returned by
:meth:`~orpheus.sn.spatial.cell_update.CellUpdate.affine_scan_coefficients`:

* ``a`` — the face transmission multiplier (``ψ_out = a·ψ_in + b``),
* ``inverse_denom`` — the reciprocal of the cell-balance diagonal ``1/S``,
* ``w`` — the **cell-average blend weight** (``ψ̄ = (1−w)ψ_in + w·ψ_out``).

The DAG-free **scan SOLVE** (forward substitution) applies these coefficients
to *factored* unknowns — it has no concrete inflow until the prefix scan has
run — so it consumes the two generic operations below, parameterized by the
coefficients and carrying NO scheme-specific constant.  The discretization math
therefore lives in exactly one place (the scheme's ``affine_scan_coefficients``),
never duplicated in a sweep body (Cardinal Rule 2 / the coefficient-model
litmus: *if an explicit-matrix representation would have to re-derive a
calculation the sweep does, that calculation belongs in the scheme*).

Why these forms are universal
=============================

Exactness on a spatially-constant flux (``ψ_in = ψ_out = ψ̄`` under a matched
source) forces the cell-average to be a **convex** face blend — the two weights
sum to one — so ``ψ̄ = (1−w)ψ_in + w·ψ_out`` holds for *any* consistent affine
scheme, with ``w`` the only free per-cell number.  The source emission then
follows algebraically (SymPy-verified, all schemes):

.. math::

    b   &= QV \cdot \mathrm{inverse\_denom} / w
       \qquad(\text{DD: } w=\tfrac12 \Rightarrow b = 2\,QV\cdot\mathrm{inverse\_denom})\\
    \bar\psi &= (1-w)\,\psi_{\rm in} + w\,\psi_{\rm out}
       \qquad(\text{DD: } w=\tfrac12 \Rightarrow \bar\psi = \tfrac12(\psi_{\rm in}+\psi_{\rm out}))

These two ops are in the ×V "denom" convention (``inverse_denom = 1/S`` with
``S`` the ×V cell-balance diagonal) — the convention the scan cache
(:class:`~orpheus.sn.spatial.sweep_cache.CollisionCache`) uses.  The
group-3 ≡ group-2 gate (``test_group3_equals_group2_scan_flat``) pins the SOLVE
direction (``ψ̄``/``ψ_out``) against the trusted Increment-A kernel.

The apply direction (matvec) is NOT a generic op here: with a CONCRETE probe
``ψ̄`` it rides the scheme's group-2
:meth:`~orpheus.sn.spatial.cell_update.CellUpdate.residual_kernel_batch` (the ÷V
``g=|μ|/Δ`` form — it returns the density residual AND the outgoing face in one
call, the natural apply twin of the scan solve).  EVERY affine scheme routes its
Cartesian matvec through this one kernel UNIFORMLY (#158/#240 — no per-scheme
matvec branch, no capability flag): Diamond Difference reproduces its diamond
march ``ψ_out = 2ψ̄ − ψ_in``, Linear-Discontinuous its Schur residual (its kernel
halves the scheme-agnostic ``s = 2|μ|/Δ`` to ``g = |μ|/Δ`` internally).  A future
``ExplicitMatrix`` representation would assemble the **×V** coefficients into
matrix entries (diagonal ``S = 1/inverse_denom``, upstream coupling
``(1−w+w·a)/inverse_denom`` — which equals DD's ``2|μ|`` and LD's ``|μ|(1+k)``);
that is the ×V matrix convention, NOT the ÷V ``residual_kernel_batch`` density
form (they differ by the cell volume ``V``).

.. note::

   **The scan SOLVE ops** (:func:`source_emission` / :func:`cell_average`) are,
   for ``w=½`` (Diamond Difference), **byte-identical** to the
   pre-coefficient-model closures: ``0.5*in + 0.5*out`` equals ``0.5*(in+out)``
   and ``QV·inv/0.5`` equals ``2·QV·inv`` bit-for-bit, because multiply/divide
   by ½ is an exact power-of-2 scaling that commutes with IEEE rounding.  So DD's
   SCAN snapshots stay byte-identical (the ``tests/sn/sweep/core tests/sn/solve
   -W error::DriftWarning`` gate is green).

   **The matvec APPLY is a deliberate principled-equivalence re-baseline.**  DD's
   Cartesian matvec moved off the ×V ``cell_balance`` density path onto the ÷V
   ``residual_kernel_batch`` kernel (#240, retiring the bit-identity-only
   ``matvec_via_kernel`` special case), which re-associates ~1 ULP on
   non-power-of-2 cell widths.  This is sanctioned by ``vv-principles``
   §"Bit-identity vs principled-equivalence": bit-identity is never a design
   constraint — the architecture (one uniform matvec kernel, no scheme flag) is
   the compounding asset; a regenerated ~1-ULP golden is the negligible cost.
   The same principled-equivalence (~1 ULP) governs LD's ×V scan vs ÷V kernel
   two-paths agreement.
"""

from __future__ import annotations

import numpy as np


def source_emission(
    QV: np.ndarray, inverse_denom: np.ndarray, w: np.ndarray,
) -> np.ndarray:
    r"""Affine-scan additive source coefficient ``b = QV · inverse_denom / w``.

    The per-cell source term of the recurrence ``ψ_out = a·ψ_in + b``.  ``QV`` is
    the volumetric source on the cell (already weight-normalised per ordinate,
    times cell volume; for a curvilinear sweep it carries the Morel–Montry
    angular contribution folded in by the caller).  DD's historical
    ``2·QV·inverse_denom`` is the ``w=½`` special case.
    """
    return QV * inverse_denom / w


def cell_average(
    face_in: np.ndarray, face_out: np.ndarray, w: np.ndarray,
) -> np.ndarray:
    r"""Cell-average from the face pair: ``ψ̄ = (1−w)·ψ_in + w·ψ_out``.

    The universal convex face blend (exactness-on-constants forces the weights
    to sum to one).  DD's ``½(ψ_in + ψ_out)`` is the ``w=½`` special case.
    """
    return (1.0 - w) * face_in + w * face_out


def outgoing_face_from_average(
    psi_bar: np.ndarray, face_in: np.ndarray, w: np.ndarray,
) -> np.ndarray:
    r"""Downstream face flux from the cell average: ``ψ_out = (ψ̄ − (1−w)·ψ_in)/w``.

    The inverse of :func:`cell_average` (the convex face blend
    ``ψ̄ = (1−w)·ψ_in + w·ψ_out``).  The universal affine-scheme outflow
    reconstruction: DD's diamond mean ``2ψ̄ − ψ_in`` is the ``w=½`` case;
    LD's ``(1+k)ψ̄ − k·ψ_in`` is the ``w=1/(1+k)`` case.

    .. note::

       For ``w=½`` (Diamond Difference) this is **byte-identical** to the
       inlined ``2·ψ̄ − ψ_in``: ``÷0.5`` is the exact power-of-2 ``×2`` and
       round-to-nearest commutes with exact doubling, so
       ``fl(2·(ψ̄ − 0.5·ψ_in)) == fl(2ψ̄ − ψ_in)`` bit-for-bit.  For LD's
       ``w=1/(1+k)`` it is algebraically equal to ``ψ̄ + (g/θ)(ψ̄ − ψ_in)/D₂``
       but takes a DIFFERENT reduction tree (compute ``w`` then ``÷w`` vs the
       direct ``ψ̄ + k·(…)``), so LD reconstruction is a principled
       ~1-ULP re-baseline (``vv-principles`` §"Bit-identity vs
       principled-equivalence"), not a byte-identical move.
    """
    return (psi_bar - (1.0 - w) * face_in) / w

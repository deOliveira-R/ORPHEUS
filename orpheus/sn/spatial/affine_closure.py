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
``ψ̄`` it rides a scheme-specific density residual (the natural apply twin of the
scan solve).  A scheme whose matvec IS its group-2
:meth:`~orpheus.sn.spatial.cell_update.CellUpdate.residual_kernel_batch` (the ÷V
``g=|μ|/Δ`` form — it returns the density residual AND the outgoing face in one
call) declares
:attr:`~orpheus.sn.spatial.cell_update.CellUpdateBase.matvec_via_kernel`
(e.g. Linear-Discontinuous).  Diamond Difference keeps its specialised
``cell_balance`` density path (byte-identical to the pre-#158 operator).  A
future ``ExplicitMatrix`` representation would assemble the **×V** coefficients
into matrix entries (diagonal ``S = 1/inverse_denom``, upstream coupling
``(1−w+w·a)/inverse_denom`` — which equals DD's ``2|μ|`` and LD's ``|μ|(1+k)``);
that is the ×V matrix convention, NOT the ÷V ``residual_kernel_batch`` density
form (they differ by the cell volume ``V``).

.. note::

   For ``w=½`` (Diamond Difference) the generic ops are **byte-identical** to
   the pre-coefficient-model closures: ``0.5*in + 0.5*out`` equals
   ``0.5*(in+out)`` and ``QV·inv/0.5`` equals ``2·QV·inv`` bit-for-bit, because
   multiply/divide by ½ is an exact power-of-2 scaling that commutes with IEEE
   rounding.  So DD stays byte-identical — its regression snapshots do NOT
   re-baseline (the gate ``tests/sn/sweep/core tests/sn/solve -W
   error::DriftWarning`` is green at 505/1/4).  **Principled-equivalence** (~1
   ULP, ``vv-principles`` §"Bit-identity vs principled-equivalence") applies only
   to ``w ≠ ½`` schemes — e.g. LD, whose ×V scan vs ÷V kernel two-paths agree at
   nULP, not bit-for-bit.
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

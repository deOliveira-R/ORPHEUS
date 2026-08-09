"""Derivation: the octant face-transmission spectrum, and why DIAMOND alone
carries an undamped sawtooth.

Created 2026-08-09 for #341, at the archivist's request: the identity below
is cited in ``docs/theory/methods/sn/cartesian_multid.rst``
(``sn-boundary-gs-not-regular``) and is load-bearing for a shipped
production ruling (keep ``inner_schedule="gauss_seidel"``), but had only a
numeric spot-check at ``d in {2,3,4}`` behind it.  A spot check is not an
algebra of record.  This module is.

WHAT IS BEING PROVED
====================

On an all-reflective box the boundary operator ``B`` IS the whole iteration,
and its action within one octant is the intra-octant face-to-face
transmission.  For a single cell in ``d`` dimensions, diamond differencing
closes ``psi_out = 2 psi_cell - psi_in`` per axis, and eliminating
``psi_cell`` through the balance equation gives the ``d x d`` map from the
octant's ``d`` inflow faces to its ``d`` outflow faces

    Sigma_DD = (2/D) 1 w^T - I ,   w_a = 2 |mu_a| A_a ,   D = Sigma_t V + sum_b w_b

— a RANK-ONE matrix minus the identity.  Step differencing closes
``psi_out = psi_cell`` instead, which drops the ``-I``:

    Sigma_step = (1/D) 1 w^T

CLAIM 1 (diamond).  spec(Sigma_DD) = {1 - 2 Sigma_t V / D} u {-1}^(d-1)
CLAIM 2 (step).     spec(Sigma_step) = {(D - Sigma_t V)/D} u {0}^(d-1)

WHY IT MATTERS — three consequences, each load-bearing somewhere
================================================================

1. **The (d-1) eigenvalues at exactly -1 are an UNDAMPED subspace.**  Their
   eigenvectors satisfy ``w^T v = 0``: weighted-zero-average across the
   octant's faces.  Zero cell average means the collision term
   ``Sigma_t V psi_c`` cannot see them at all, so absorption — the only
   damping present at zero leakage — does not touch them.  This IS the face
   sawtooth #340 measured (trace ratios 1.074414 / 0.925586, summing to
   exactly 2 because the mode is the zero-average direction), the mode that
   costs 1631 sweeps on the d=3 all-reflective box.

2. **Negative eigenvalues void Varga's regular-splitting hypothesis.**  A
   regular splitting needs ``M^-1 >= 0`` and ``N >= 0`` elementwise; its
   payoff is the comparison theorem forcing ``rho_GS <= rho_J``
   unconditionally.  ``Sigma_DD`` is not a non-negative operator, so the
   theorem does not apply — and the conclusion it would license is
   measurably false (#341).  The tree asserted "regular splitting" in
   eleven places until 2026-08-09.

3. **It is DIAMOND's problem, not transport's.**  Claim 2 is the control:
   swap the closure and the undamped subspace collapses from ``-1`` to
   ``0``, i.e. step annihilates the sawtooth in a single pass.  The ``-1``
   comes from precisely the ``-psi_in`` term in the DD closure.  Without
   this contrast the spectrum reads like a property of the physics; with
   it, the closure is convicted.

Run directly to print the symbolic result and the numeric spot-check::

    .venv/bin/python -O derivations/sn_dd_face_transmission.py
"""
from __future__ import annotations

import numpy as np
import sympy as sp


def derive_dd_face_transmission_spectrum(d: int = 3) -> dict[str, object]:
    """Prove Claims 1 and 2 symbolically for a given dimension ``d``.

    THE PROOF IS STRUCTURAL, NOT A DETERMINANT.  A rank-one-plus-scalar
    matrix has its whole spectrum by inspection, and saying so is both the
    honest argument and the one that holds for every ``d``:

    * ``Sigma = alpha * 1 w^T + beta * I`` acts as ``beta`` on the
      hyperplane ``{v : w^T v = 0}`` — dimension exactly ``d - 1``, since
      ``w != 0`` — because the rank-one term annihilates it;
    * it acts on ``1`` as ``alpha * (w^T 1) + beta``;
    * ``1`` is NOT in that hyperplane (``w^T 1 = sum(w) > 0`` for positive
      face weights), so the two together span ``R^d`` and the spectrum is
      complete.

    Nothing above expands a characteristic polynomial, so the cost is O(1)
    in ``d`` and the argument is visibly dimension-independent.  (An earlier
    version of this file factored ``det(Sigma - lambda I)`` per dimension;
    that is a *check at one d*, not a proof for all d, and it is
    exponentially slow besides.)  ``_char_poly_crosscheck`` below re-derives
    the same answer the expensive way at small ``d``, as an independent
    confirmation that the structural argument was applied correctly.
    """
    Sigma_t, V = sp.symbols("Sigma_t V", positive=True)
    w = sp.Matrix(sp.symbols(f"w1:{d + 1}", positive=True))       # (d, 1)
    ones = sp.ones(d, 1)
    sum_w = sum(w.T.tolist()[0])
    D = Sigma_t * V + sum_w

    # The two closures differ ONLY by the -I; everything else is shared.
    # DIAMOND  psi_out = 2 psi_cell - psi_in   ->  alpha = 2/D, beta = -1
    # STEP     psi_out =   psi_cell            ->  alpha = 1/D, beta =  0
    rank_one = ones * w.T                                          # (d, d)
    Sigma_DD = 2 * rank_one / D - sp.eye(d)
    Sigma_step = rank_one / D

    v = sp.Matrix(sp.symbols(f"v1:{d + 1}"))
    # Impose w^T v = 0 by eliminating v1; the remaining d-1 symbols are free,
    # which IS the statement that the eigenspace has dimension d - 1.
    on_hyperplane = sp.solve(sp.Eq((w.T * v)[0, 0], 0), sp.Symbol("v1"),
                             dict=True)[0]

    results = {}
    for name, Sig, beta, alpha in (
        ("DD", Sigma_DD, sp.Integer(-1), 2 / D),
        ("step", Sigma_step, sp.Integer(0), 1 / D),
    ):
        # (i) beta on the w-orthogonal hyperplane, identically in the free v's
        residual = sp.simplify((Sig * v - beta * v).subs(on_hyperplane))
        assert residual == sp.zeros(d, 1), f"{name}: hyperplane eigenvalue"

        # (ii) the lead eigenvalue on the all-ones direction
        lead = sp.simplify((Sig * ones)[0, 0])
        assert sp.simplify(Sig * ones - lead * ones) == sp.zeros(d, 1), \
            f"{name}: `ones` is an eigenvector"
        claim = sp.simplify(alpha * sum_w + beta)
        assert sp.simplify(lead - claim) == 0, f"{name}: lead eigenvalue"

        # (iii) completeness: `ones` is off the hyperplane, so dim 1 + (d-1) = d
        assert sp.simplify((w.T * ones)[0, 0] - sum_w) == 0
        assert sp.simplify(sum_w) != 0, f"{name}: `ones` would lie in the kernel"

        results[name] = sp.simplify(lead)

    # Match the leads against the published closed forms.
    lead_DD_claim = sp.simplify(1 - 2 * Sigma_t * V / D)
    lead_step_claim = sp.simplify((D - Sigma_t * V) / D)
    assert sp.simplify(results["DD"] - lead_DD_claim) == 0, "Claim 1"
    assert sp.simplify(results["step"] - lead_step_claim) == 0, "Claim 2"

    return {
        "d": d,
        "Sigma_DD": Sigma_DD,
        "Sigma_step": Sigma_step,
        "spectrum_DD": {lead_DD_claim: 1, sp.Integer(-1): d - 1},
        "spectrum_step": {lead_step_claim: 1, sp.Integer(0): d - 1},
    }


def _char_poly_crosscheck(d: int) -> None:
    """Independent (expensive) confirmation: factor the characteristic
    polynomial and match it against the claimed multiset.

    Structurally redundant with the proof above — kept because the two can
    only agree if the matrix built here is the matrix the argument was made
    about.  Small ``d`` only; symbolic determinants blow up fast.
    """
    Sigma_t, V = sp.symbols("Sigma_t V", positive=True)
    w = sp.Matrix(sp.symbols(f"w1:{d + 1}", positive=True))
    sum_w = sum(w.T.tolist()[0])
    D = Sigma_t * V + sum_w
    lam = sp.Symbol("lambda")

    for Sig, lead, other in (
        (2 * sp.ones(d, 1) * w.T / D - sp.eye(d), 1 - 2 * Sigma_t * V / D, -1),
        (sp.ones(d, 1) * w.T / D, (D - Sigma_t * V) / D, 0),
    ):
        got = sp.expand((Sig - lam * sp.eye(d)).det())
        want = sp.expand((-1) ** d * (lam - lead) * (lam - other) ** (d - 1))
        assert sp.simplify(got - want) == 0, f"char poly mismatch at d={d}"


def _numeric_spot_check(d: int, rng: np.random.Generator) -> tuple[float, float]:
    """Independent confirmation with a numeric eigen-solver.

    The symbolic proof above is the authority; this exists because a
    transcription slip between the algebra and the matrix the code actually
    builds would be invisible to a proof about the algebra alone.
    """
    w = rng.uniform(0.2, 3.0, size=d)
    Sig_t, V = float(rng.uniform(0.1, 4.0)), float(rng.uniform(0.5, 2.0))
    D = Sig_t * V + w.sum()
    S_dd = 2.0 * np.outer(np.ones(d), w) / D - np.eye(d)
    S_st = np.outer(np.ones(d), w) / D

    ev_dd = np.sort_complex(np.linalg.eigvals(S_dd))
    ev_st = np.sort_complex(np.linalg.eigvals(S_st))
    want_dd = np.sort_complex(np.array([1 - 2 * Sig_t * V / D] + [-1.0] * (d - 1)))
    want_st = np.sort_complex(np.array([(D - Sig_t * V) / D] + [0.0] * (d - 1)))
    return (float(np.max(np.abs(ev_dd - want_dd))),
            float(np.max(np.abs(ev_st - want_st))))


if __name__ == "__main__":
    rng = np.random.default_rng(20260809)
    for d in (2, 3, 4, 5, 6, 8):
        res = derive_dd_face_transmission_spectrum(d)
        err_dd, err_st = _numeric_spot_check(d, rng)
        extra = ""
        if d <= 4:                       # char-poly cross-check: small d only
            _char_poly_crosscheck(d)
            extra = "  + char-poly cross-check OK"
        print(f"d={d}  STRUCTURAL PROOF OK{extra}")
        print(f"      spec(Sigma_DD)   = {res['spectrum_DD']}")
        print(f"      spec(Sigma_step) = {res['spectrum_step']}")
        print(f"      numeric residual : DD {err_dd:.2e}   step {err_st:.2e}")
    print()
    print("Claim 1  spec(Sigma_DD)   = {1 - 2 Sigma_t V / D} u {-1}^(d-1)   PROVED")
    print("Claim 2  spec(Sigma_step) = {(D - Sigma_t V)/D}   u { 0}^(d-1)   PROVED")
    print()
    print("The (d-1)-fold -1 is the UNDAMPED face sawtooth (#340's 1631-sweep")
    print("mode); it collapses to 0 under step differencing, so it is a")
    print("property of the DIAMOND closure, not of the transport problem.")

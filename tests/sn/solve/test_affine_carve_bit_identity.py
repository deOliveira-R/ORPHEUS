r"""Bit-identity gate for the #208 affine flux-algebra carve.

The affine typing (Piece 1: ``FluxDisplacement`` + the ``flux+flux`` gate),
the SI-loop displacement retype (Piece 2), and the typed equation-residual
evaluation (Piece 3) add **ZERO numerical change** to the converged flux:

* The SI stopping norm STAYS the flat ``_l2_norm(displacement.to_flat())``
  (the displacement composite's ``to_flat()`` reads the SAME ``.values`` as the
  pre-carve ``(psi-psi_prev).to_flat()`` — see the displacement mechanism in
  ``.claude/plans/issue_208_residual_typing_closeout.md`` §3a).
* The residual evaluation is additive / diagnostic — never in the convergence
  path.

So the converged flux MUST be **byte-identical** to the pre-carve value. This
test freezes a ``sha256`` of the converged ``angular_flux.interior.values`` and
``scalar_flux.values``, generated once at the pre-carve commit ``63719a2`` and
committed WITH Piece 1.

Why a dedicated golden, not the DD regression snapshots
=======================================================

``tests/sn/regression/test_dd_regression.py`` compares iterative results at
``SAFETY × conv_tol ≈ 1e-11`` (escalatable to strict only via
``-W error::DriftWarning``), and the ``2d_2g_p1_aniso_dd_8x4_het_si`` snapshot
ALREADY pre-drifts ~6920 ULP / ~9.8e-13 (Phase-5b/5c inheritance) — so the
snapshot suite CANNOT verify this carve's stronger zero-change claim. A
``sha256`` of the raw bytes is the sharpest (sub-ULP) bit-identity assertion.

Coverage (Cardinal-6 ≥2G + heterogeneous; vv-principles §H1/H2)
===============================================================

* ``si_2d_p1_aniso_het`` — 2-D Cartesian, 2G, fuel|moderator, P1-anisotropic,
  source-iteration: the **windowed SI** path whose iterate ``bulk`` is a
  :class:`HarmonicMomentFlux` (so ``psi-psi_prev`` becomes a
  ``MomentDisplacement``). The mirror of the frozen snapshot config.
* ``krylov_2d_p1_aniso_het`` — same config, full-angular Krylov (scipy's flat
  ``b-Ax`` — exercises the typing without the SI displacement path).
* ``si_slab_2g_het`` — 1-D slab, 2G, fuel|moderator, P1, source-iteration: the
  ``AngularFlux``-bulk SI path (1-D never windows → ``AngularDisplacement``).

Marks: ``foundation`` — this is a software-invariant gate (the converged value
is inherited from the pre-carve reference), not a physics ``:label:`` claim.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D
from orpheus.geometry.mesh import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source


pytestmark = pytest.mark.foundation


# Generated at pre-carve HEAD 63719a2 (commit with Piece 1). sha256 of the
# C-contiguous float64 bytes of the converged fields.
#
# REGENERATION HISTORY:
# * 2026-06-11 (S6.9 Fork-B2, #222) — the FOUR 2-D hashes regenerated in the
#   default-flip commit: the multi-D Cartesian production default changed
#   MovingFrontierWindow → ScanMarch, a SCHEDULE change (row-march vs
#   anti-diagonal) that shifts the converged bytes at FP-association level
#   (principled-equivalent, NOT a numerics change — vv §bit-identity-vs-
#   principled).  Output-identity evidence for the flip: the G4.a/G4.b
#   Mode-9 FP-invariance gates (test_scan_march_end_to_end.py, scanmarch
#   default ≡ forced window at solver tol) + the ScanMarch G2.c nulp oracle.
#   The 1-D slab hashes are UNTOUCHED (CumprodScan stays the 1-D default —
#   the flip's blast radius pin).
# * 2026-06-15 (#240 Phase 2 Step B) — the TWO ``krylov_2d_p1_aniso_het`` hashes
#   regenerated. The Step-B carve made ``InvertibleOperator.apply`` OWN its
#   matvec (``loss_action(self.sigma)`` direct) instead of the inherited leaf
#   sum ``(loss_action(σ_t) − σ_t·ψ) + σ_t·ψ`` — the override DROPS that
#   ``−σ_t·ψ + σ_t·ψ`` round-trip, so the 2-D Cartesian matvec re-associates
#   ≤5 ULP per application; accumulated through GMRES (inner_tol=1e-12) the
#   converged Krylov bytes shift at the ~1e-12 GMRES-tolerance level. Verified
#   structurally-independent: the Krylov converged φ agrees with the SI 2-D φ
#   (which does NOT use the apply override — bit-identical to ITS golden) to
#   3.9e-12 rel (vv-principles 3-criteria: named ``loss_action`` / SI cross-
#   check / drift = GMRES-tol × ULP). The SI 2-D + slab hashes are UNTOUCHED
#   (SI rides ``solve``, not ``apply`` — the carve's apply-only blast-radius
#   pin).
# * 2026-06-16 (#240 Phase 2 Step D5a) — ALL FOUR 2-D hashes regenerated (BOTH
#   ``si_2d_p1_aniso_het`` AND ``krylov_2d_p1_aniso_het``). D5a folded the 2-D
#   row-march ``ScanMarch._sweep_interior`` (SOLVE) + ``_loss_action_interior``
#   (matvec) off inline DiamondDifference onto the scheme's coefficient model
#   (#239): the SOLVE now consumes ``scheme.cartesian_scan_coefficients`` →
#   ``(a, inverse_denom, w)`` + ``source_emission``/``cell_average`` (the
#   ``×inverse_denom`` reciprocal form, replacing the ONE remaining inline ``÷D``
#   DIVISION — so this is the same byte re-association the 1-D CumprodScan
#   already rides, NOT a numerics change); the matvec rides
#   ``scheme.residual_kernel_batch`` + ``reflect_scan_coefficients``. Both folds
#   re-associate the 2-D cell-balance reduction ~1 ULP. The SI shift accumulates
#   through source iteration (iter_count × ULP); the Krylov shift through GMRES
#   (GMRES-tol × ULP). Verified structurally-independent: the post-fold SI 2-D φ
#   (``.solve`` fold) and Krylov 2-D φ (``.apply`` fold) — DIFFERENT code paths —
#   agree to 4.2e-12 rel, AND the ScanMarch≡FullFieldWavefront oracle
#   (test_scan_march_equivalence.py) pins the value to the analytical
#   ``k_inf``/``φ=Q/Σ_t`` grounds (vv-principles §"Bit-identity vs principled-
#   equivalence" 3-criteria: named coefficient-model ops / two-path + oracle
#   cross-check / drift = iter-or-GMRES-tol × ULP). The 1-D slab hashes are
#   UNTOUCHED — D5a's blast radius is the 2-D row-march ONLY (the 1-D CumprodScan
#   scan path is byte-identical), which is exactly the D5a.3 negative-control
#   proof that the fold did not leak into the 1-D solve.
GOLDEN = {
    "si_2d_p1_aniso_het_psi_sha":
        "847f302c28cc3e096252412d7ce266f982c112b416178ddbd40be7a62d97459c",
    "si_2d_p1_aniso_het_phi_sha":
        "688e65b34c3563a54ae39bedf22ff533cd2518bdd6462752efe191912dc45713",
    "krylov_2d_p1_aniso_het_psi_sha":
        "ae9b35efb3b4f3b09411db1dcd9d1201a60985bb7cbd7087e6dddbcba0a1a196",
    "krylov_2d_p1_aniso_het_phi_sha":
        "0c8a6ce13c9fcdd831fe421e59a5d24cc6571bb5b323edf2063571bda6cca749",
    "si_slab_2g_het_psi_sha":
        "353d7db054781af44dc4682ca3330c0c7490d54185bf5a3a8a83b83b85b4b1f3",
    "si_slab_2g_het_phi_sha":
        "136a3e1779e2f280fb6dd8a99eae825430855e1fde0ff6c09ea9902e2151b028",
}


def _sha(arr: np.ndarray) -> str:
    a = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
    return hashlib.sha256(a.tobytes()).hexdigest()


def _build_2d():
    """2-D Cartesian 2G het-aniso (mirror of ``2d_2g_p1_aniso_dd_8x4_het_si``)."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    sum_w = float(quad.weights.sum())
    q_ext = np.full((quad.N, 2, nx, ny), 1.0 / sum_w)
    return {2: fuel, 0: mod}, mesh, quad, q_ext


def _build_slab():
    """1-D slab 2G het (fuel|moderator split, vacuum BC)."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    nx = 16
    mat_ids = np.zeros(nx, dtype=int)
    mat_ids[: nx // 2] = 2
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, nx + 1),
        mat_ids=mat_ids,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sum_w = float(quad.weights.sum())
    q_ext = np.full((quad.N, 2, nx), 1.0 / sum_w)
    return {2: fuel, 0: mod}, mesh, quad, q_ext


def _solve_case(case: str):
    if case == "si_2d_p1_aniso_het":
        mats, mesh, quad, q = _build_2d()
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="source_iteration", max_inner=3000, inner_tol=1e-12,
        )
    if case == "krylov_2d_p1_aniso_het":
        mats, mesh, quad, q = _build_2d()
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="krylov", max_inner=3000, inner_tol=1e-12,
        )
    if case == "si_slab_2g_het":
        mats, mesh, quad, q = _build_slab()
        return solve_sn_fixed_source(
            mats, mesh, quad, q, scattering_order=1,
            inner_solver="source_iteration", max_inner=3000, inner_tol=1e-12,
        )
    raise ValueError(f"unknown case {case!r}")


@pytest.mark.parametrize(
    "case",
    ["si_2d_p1_aniso_het", "krylov_2d_p1_aniso_het", "si_slab_2g_het"],
)
def test_converged_flux_bit_identical_after_affine_carve(case: str) -> None:
    r"""The converged flux is BIT-IDENTICAL to the pre-carve (63719a2) golden.

    The affine typing + the SI displacement retype + the residual evaluation
    add ZERO numerical change; the converged ``psi`` / ``phi`` bytes hash to the
    frozen pre-carve ``sha256``. A drift here means the carve changed the
    numerics (the DD snapshot suite, gated at ``≈1e-11`` and pre-drifting, would
    NOT catch a sub-1e-11 regression — this dedicated gate does, sub-ULP).
    """
    sol = _solve_case(case)
    if not sol.history.converged:
        raise AssertionError(f"{case}: solve did not converge")
    psi_sha = _sha(sol.angular_flux.interior.values)
    phi_sha = _sha(sol.scalar_flux.values)
    if psi_sha != GOLDEN[f"{case}_psi_sha"]:
        raise AssertionError(
            f"{case}: converged psi NOT bit-identical to pre-carve golden "
            f"(sha {psi_sha} != {GOLDEN[f'{case}_psi_sha']}) — the carve "
            f"changed the numerics."
        )
    if phi_sha != GOLDEN[f"{case}_phi_sha"]:
        raise AssertionError(
            f"{case}: converged phi NOT bit-identical to pre-carve golden "
            f"(sha {phi_sha} != {GOLDEN[f'{case}_phi_sha']})."
        )

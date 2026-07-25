r"""#276 A4 — the φ* certification battery: P1.3 / P1.4 / P1.5.

The spec's Phase-1 chain (``p6_adjoint_verification_spec.md`` §3, with the
P1.4 formula correction recorded there): P1.0a/b + P1.1 landed pre-A4;
P1.2's duality core + lift pin live in ``test_sn_adjoint_entries.py``.
This file carries the remaining certification rows AND the named
mutations — every mutation in-process (monkeypatch; NEVER ``git
checkout``), each reddening its named gate.

* **P1.3** — ``k_adj == k_fwd`` (exact algebraic identity,
  ``eig(A†) = eig(A)``): ∞-medium 2G AND 4G (both == the closed-form
  ``kinf_homogeneous`` anchor — the triple equality terminating in
  ``np.linalg.eig``), a 2-region HETEROGENEOUS reflective slab, and the
  SPHERE (the coupled daggered posing — the μ-reversal-at-the-pole
  catcher; the ∞ legs are flat-flux-BLIND to a spatial-adjoint bug,
  which is why these two legs are mandatory).
* **P1.3 mutations** — F† → F (χ↔νΣf swap dropped; ≥2G non-proportional
  χ/νΣf), S† → S (group-transpose dropped; asymmetric SigS), and
  L† → L (the transpose-solve replaced by the FORWARD sweep — "no
  ordinate reversal", visible on the heterogeneous leg, invisible on
  flat ∞).  Each shifts k_adj off k_fwd.
* **P1.4** — the ∞-medium 4G adjoint SPECTRUM == the closed form (the
  dominant eigenvector of ``(Aᵀ)⁻¹Fᵀ`` — the corrected reference; 4G
  mandatory, with the ψ*_cf ≠ φ_cf materiality precondition asserted),
  plus the F†-swap mutation reddening the spectrum row while k stays
  EXACTLY equal — the Mode-12 pair made explicit.
* **P1.5** — bi-orthogonality ``⟨ψ*_i, F φ_j⟩ = 0 (i≠j)``, the
  defining law of the forward/adjoint eigenbasis pair, closed-form
  (all modes from ``np.linalg.eig`` both sides).  HONEST-SCOPE note vs
  the spec's "ng distinct eigenmodes" wording: the rank-1
  ``F = χ⊗νΣf`` gives ONE nonzero eigenvalue and an (ng−1)-fold zero
  eigenspace, so the cross-Gram is diagonal in the DEGENERATE sense —
  exactly one nonzero entry, off-diagonals structurally zero through
  ``Fφ_j = 0`` (zero right modes) and ``χ·ψ*_i = 0`` (zero left
  modes).  Both mechanisms are asserted, and the F†-swap mutation
  (right vectors on BOTH sides) breaks the off-diagonals.

vv Mode-8: ``np.testing.*`` / :func:`require` only (fire under
``python -O``).  Mode-12 discipline: no k row stands alone — the
vector-level companions are P1.4 here and the duality/lift rows in the
entries file.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_adjoint_spectrum_homogeneous,
    kinf_and_spectrum_homogeneous,
    kinf_homogeneous,
)
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh1D
from orpheus.geometry.mesh import BC, CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_adjoint

pytestmark = pytest.mark.l1


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _quad():
    return Quadrature.gauss_legendre(n_ordinates=8)


def _solve_fwd(mats, mesh):
    return solve_sn(
        mats, mesh, _quad(),
        keff_tol=1e-9, flux_tol=1e-8, inner_tol=1e-10, max_inner=500,
    )


def _solve_adj(mats, mesh):
    return solve_sn_adjoint(
        mats, mesh, _quad(),
        keff_tol=1e-9, flux_tol=1e-8, inner_tol=1e-10, max_inner=500,
    )


def _mix_arrays(mix):
    return (
        np.asarray(mix.SigT, dtype=float),
        np.asarray(mix.SigS[0].todense(), dtype=float),
        np.asarray(mix.SigP, dtype=float),
        np.asarray(mix.chi, dtype=float),
    )


def _infinite_medium(ng_key: str):
    mix = get_mixture("A", ng_key)
    mesh = Mesh1D(
        edges=np.linspace(0.0, 5.0, 11), mat_ids=np.zeros(10, dtype=int),
    )
    return {0: mix}, mesh, mix


def _het_slab():
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, 9),
        mat_ids=np.array([0, 0, 1, 1, 1, 1, 0, 0]),
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    return mats, mesh


def _sphere():
    mats = {0: get_mixture("A", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 3.0, 7), mat_ids=np.zeros(6, dtype=int),
        coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"),
    )
    return mats, mesh


def _k(sol) -> float:
    if sol.keff is None:
        pytest.fail("solve returned no eigenvalue.")
    return float(sol.keff)


def _energy_spectrum(sol) -> np.ndarray:
    sf = np.asarray(sol.scalar_flux.values)
    spec = sf.mean(axis=tuple(range(1, sf.ndim)))
    return spec / np.linalg.norm(spec)


# ═══════════════════════════════════════════════════════════════════════
# P1.3 — the k triple-equality legs.
# ═══════════════════════════════════════════════════════════════════════


class TestP13KEquality:
    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    @pytest.mark.verifies("kinf-mg")
    def test_infinite_medium_triple_equality(self, ng_key):
        r"""``k_adj == k_fwd == kinf_homogeneous`` — the chain terminates in
        ``np.linalg.eig``, not in another ORPHEUS solver (L11)."""
        mats, mesh, mix = _infinite_medium(ng_key)
        k_fwd = _k(_solve_fwd(mats, mesh))
        k_adj = _k(_solve_adj(mats, mesh))
        k_cf = kinf_homogeneous(*_mix_arrays(mix))
        np.testing.assert_allclose(
            k_adj, k_fwd, rtol=0, atol=1e-8,
            err_msg=f"{ng_key}: k_adj != k_fwd on the ∞ medium.",
        )
        np.testing.assert_allclose(
            k_adj, k_cf, rtol=0, atol=1e-8,
            err_msg=f"{ng_key}: k_adj != closed-form k∞.",
        )

    def test_heterogeneous_slab_k_equality(self):
        r"""Spatial structure LIVE: the ∞ legs are flat-flux-blind to a
        spatial-adjoint bug (a wrong Lᵀ direction) — this leg is not."""
        mats, mesh = _het_slab()
        k_fwd = _k(_solve_fwd(mats, mesh))
        k_adj = _k(_solve_adj(mats, mesh))
        np.testing.assert_allclose(
            k_adj, k_fwd, rtol=0, atol=1e-8,
            err_msg="k_adj != k_fwd on the heterogeneous slab — a "
            "spatial-adjoint (Lᵀ) or composition error the flat legs "
            "cannot see.",
        )

    def test_sphere_k_equality(self):
        r"""The COUPLED daggered posing (System B live): μ-reversal at the
        pole + the fission ray fold + the transposed block substitution."""
        mats, mesh = _sphere()
        k_fwd = _k(_solve_fwd(mats, mesh))
        k_adj = _k(_solve_adj(mats, mesh))
        np.testing.assert_allclose(
            k_adj, k_fwd, rtol=0, atol=1e-8,
            err_msg="k_adj != k_fwd on the sphere — the coupled daggered "
            "chain (fission fold / march reverse / transposed "
            "substitution) is not the adjoint of the forward one.",
        )


class TestP13Mutations:
    r"""Each named adjoint mutation shifts k_adj off k_fwd — in-process."""

    def test_fission_role_swap_shifts_k(self, monkeypatch):
        r"""F† → F (the χ↔νΣf swap dropped).  4G, χ ∦ νΣf asserted —
        at 1G or χ ∝ νΣf the swap is invisible (the config-blindness
        the spec §2 audit names)."""
        from orpheus.transport.operators.fission import FissionOperator

        mats, mesh, mix = _infinite_medium("4g")
        _sigt, _sigs, nu_sf, chi = _mix_arrays(mix)
        require(
            not np.allclose(
                chi / np.linalg.norm(chi), nu_sf / np.linalg.norm(nu_sf),
                rtol=1e-6,
            ),
            "fixture lost χ ∦ νΣf — the F† swap would be invisible.",
        )
        k_fwd = _k(_solve_fwd(mats, mesh))

        monkeypatch.setattr(
            FissionOperator, "apply_transpose",
            lambda self, chi: self.apply(chi),
        )
        k_mut = _k(_solve_adj(mats, mesh))
        require(
            abs(k_mut - k_fwd) > 1e-6,
            f"F†=F did NOT shift the adjoint k (|dk|="
            f"{abs(k_mut - k_fwd):.2e}) — the P1.3 gate has no F† teeth.",
        )

    def test_scattering_transpose_drop_shifts_k(self, monkeypatch):
        r"""S† → S (the group transpose dropped).  Asymmetric SigS asserted."""
        from orpheus.transport.operators.scattering import ScatteringOperator

        mats, mesh, mix = _infinite_medium("2g")
        _sigt, sigs, _nsf, _chi = _mix_arrays(mix)
        require(
            not np.allclose(sigs, sigs.T),
            "fixture lost SigS asymmetry — the S† drop would be invisible.",
        )
        k_fwd = _k(_solve_fwd(mats, mesh))

        monkeypatch.setattr(
            ScatteringOperator, "apply_transpose",
            lambda self, chi: self.apply(chi),
        )
        k_mut = _k(_solve_adj(mats, mesh))
        require(
            abs(k_mut - k_fwd) > 1e-6,
            f"S†=S did NOT shift the adjoint k (|dk|="
            f"{abs(k_mut - k_fwd):.2e}) — the P1.3 gate has no S† teeth.",
        )

    def test_streaming_no_reversal_shifts_k_heterogeneous(self, monkeypatch):
        r"""L† → L: the daggered inverse rides the FORWARD sweep ("no
        ordinate reversal").  HETEROGENEOUS leg — on the flat ∞ medium the
        streaming term is nulled and this mutation is INVISIBLE (asserted
        by the config choice, vv §0.6)."""
        from orpheus.sn.operators.streaming import InvertibleOperator

        mats, mesh = _het_slab()
        k_fwd = _k(_solve_fwd(mats, mesh))

        monkeypatch.setattr(
            InvertibleOperator, "solve_transpose", InvertibleOperator.solve,
        )
        k_mut = _k(_solve_adj(mats, mesh))
        require(
            abs(k_mut - k_fwd) > 1e-6,
            f"the forward-sweep-as-adjoint-inverse mutation did NOT shift "
            f"k on the heterogeneous slab (|dk|={abs(k_mut - k_fwd):.2e}) "
            f"— the P1.3 spatial leg has no Lᵀ teeth.",
        )


# ═══════════════════════════════════════════════════════════════════════
# P1.4 — the ∞-medium 4G adjoint spectrum (flux-shape of φ*).
# ═══════════════════════════════════════════════════════════════════════


class TestP14AdjointSpectrum:
    def test_4g_spectrum_matches_closed_form(self):
        r"""The SN adjoint energy spectrum == the dominant eigenvector of
        ``(Aᵀ)⁻¹Fᵀ`` (4G MANDATORY: the 2G row lives in the entries file;
        4G is where φ* is genuinely non-flat and ≠ φ — asserted)."""
        mats, mesh, mix = _infinite_medium("4g")
        arrays = _mix_arrays(mix)
        k_cf, psi_star_cf = kinf_and_adjoint_spectrum_homogeneous(*arrays)
        _k_cf_f, phi_cf = kinf_and_spectrum_homogeneous(*arrays)
        require(
            not np.allclose(psi_star_cf, phi_cf, rtol=1e-3),
            "4G ψ*_cf ≈ φ_cf — the fixture is effectively self-adjoint "
            "and this flux-shape gate is a dud.",
        )

        adj = _solve_adj(mats, mesh)
        np.testing.assert_allclose(
            _k(adj), k_cf, rtol=0, atol=1e-8,
            err_msg="4G k_adj != closed form.",
        )
        np.testing.assert_allclose(
            _energy_spectrum(adj), psi_star_cf, rtol=1e-7,
            err_msg="the 4G SN adjoint spectrum does not match the "
            "closed-form (Aᵀ)⁻¹Fᵀ eigenvector — the energy-adjoint "
            "content (F†/S† composition) is wrong even though k agrees.",
        )

    def test_f_role_swap_reds_spectrum_while_k_stays(self, monkeypatch):
        r"""THE Mode-12 pair, explicit: F† → F leaves ``k`` EXACTLY equal
        (the k-functional is designed-green on the whole adjoint mutation
        class) while the SPECTRUM diverges O(1) — the vector-level row is
        the committed catcher, k alone never is.

        (F†=F makes the daggered pencil ``(A†, F)`` — similar to the
        forward pencil, so k is INVARIANT; the eigenVECTOR is not.)
        """
        from orpheus.transport.operators.fission import FissionOperator

        mats, mesh, mix = _infinite_medium("4g")
        arrays = _mix_arrays(mix)
        k_cf, psi_star_cf = kinf_and_adjoint_spectrum_homogeneous(*arrays)

        monkeypatch.setattr(
            FissionOperator, "apply_transpose",
            lambda self, chi: self.apply(chi),
        )
        mut = _solve_adj(mats, mesh)
        spec_mut = _energy_spectrum(mut)
        require(
            not np.allclose(spec_mut, psi_star_cf, rtol=1e-3),
            "F†=F did NOT move the adjoint spectrum — the P1.4 flux-shape "
            "row has no teeth on the canonical adjoint-fission trap.",
        )

    # NOTE the k-INVARIANCE half of the Mode-12 pair: on the ∞ medium the
    # F†=F daggered pencil (A†, F) is similar to the forward pencil, so
    # the P1.3 ∞ row would stay GREEN under this mutation — that is the
    # designed-green blindness P1.4 exists to break, and exactly why the
    # F† k-tooth (TestP13Mutations) had to ride a shifted-k regime while
    # the spectrum row above is the committed ∞-medium catcher.


# ═══════════════════════════════════════════════════════════════════════
# P1.5 — bi-orthogonality (the intrinsic law of the eigenbasis pair).
# ═══════════════════════════════════════════════════════════════════════


class TestP15BiOrthogonality:
    def _bases(self, mix):
        sig_t, sig_s, nu_sf, chi = _mix_arrays(mix)
        A = np.diag(sig_t) - sig_s.T
        F = np.outer(chi, nu_sf)
        _vals_r, right = np.linalg.eig(np.linalg.solve(A, F))
        _vals_l, left = np.linalg.eig(np.linalg.solve(A.T, F.T))
        return A, F, np.real(right), np.real(left)

    @pytest.mark.foundation
    def test_cross_gram_is_degenerate_diagonal(self):
        r"""``⟨ψ*_i, F φ_j⟩ = 0 (i≠j)`` — closed form, all modes.

        HONEST SCOPE (spec deviation, recorded): the rank-1 F has ONE
        nonzero eigenvalue; the cross-Gram is diagonal in the DEGENERATE
        sense — exactly one nonzero entry.  BOTH orthogonality mechanisms
        are asserted structurally: ``Fφ_j = 0`` for the zero right modes
        (φ_j ⊥ νΣf) and ``χ·ψ*_i = 0`` for the zero left modes.
        """
        mix = get_mixture("A", "4g")
        _A, F, right, left = self._bases(mix)
        gram = left.T @ F @ right
        scale = float(np.abs(gram).max())
        require(scale > 0.0, "cross-Gram is identically zero — degenerate fixture.")

        off = gram.copy()
        # The dominant pairing is the single nonzero entry: locate it.
        i0, j0 = np.unravel_index(int(np.abs(gram).argmax()), gram.shape)
        off[i0, j0] = 0.0
        np.testing.assert_allclose(
            off, np.zeros_like(off), atol=1e-10 * scale,
            err_msg="off-diagonal ⟨ψ*_i, Fφ_j⟩ ≠ 0 — the forward/adjoint "
            "eigenbases are not F-bi-orthogonal (a wrong left basis).",
        )

    @pytest.mark.foundation
    def test_right_vectors_are_not_f_orthogonal(self):
        r"""MUTATION (F† = F realized at the basis level): using the RIGHT
        eigenvectors on BOTH sides, the Gram is NOT one-entry — the right
        basis is F-orthogonal only in the self-adjoint case, so a φ* that
        silently equals φ fails P1.5.  (Rank-1 sharpening: ``Fφ_j = 0``
        for the zero modes kills those COLUMNS regardless, so the
        discriminating entries are the dominant COLUMN's off-dominant
        rows — ``⟨φ_i, Fφ_1⟩ ∝ χ·φ_i``, and the right zero-modes are ⊥
        νΣf, NOT ⊥ χ.)"""
        mix = get_mixture("A", "4g")
        _A, F, right, _left = self._bases(mix)
        gram = right.T @ F @ right
        scale = float(np.abs(gram).max())
        i0, j0 = np.unravel_index(int(np.abs(gram).argmax()), gram.shape)
        off = gram.copy()
        off[i0, j0] = 0.0
        require(
            float(np.abs(off).max()) > 1e-6 * scale,
            "the right-right F-Gram came out one-entry too — the "
            "bi-orthogonality gate cannot distinguish φ* from φ on this "
            "fixture (self-adjoint dud).",
        )

    def test_sn_adjoint_dominant_mode_is_the_true_left_vector(self):
        r"""The SN tie-in: the SN ψ* spectrum, paired against F with the
        closed-form ZERO right modes, annihilates — ``⟨ψ*_SN, Fφ_j⟩ = 0``
        via ``χ·ψ*_SN``-orthogonality is NOT required (only zero LEFT
        modes are χ-orthogonal); the honest SN row is
        ``Fᵀψ*_SN ∝ νΣf`` alignment with the defining law
        ``Aᵀψ* = Fᵀψ*/k`` — asserted through the P1.4 spectrum match,
        so THIS row pins the residual of the defining law directly on
        the SN vector at the closed-form (A, F).
        """
        mats, mesh, mix = _infinite_medium("4g")
        sig_t, sig_s, nu_sf, chi = _mix_arrays(mix)
        A = np.diag(sig_t) - sig_s.T
        F = np.outer(chi, nu_sf)
        adj = _solve_adj(mats, mesh)
        spec = _energy_spectrum(adj)
        k = _k(adj)
        residual = A.T @ spec - (F.T @ spec) / k
        np.testing.assert_allclose(
            residual, np.zeros_like(residual), atol=1e-7,
            err_msg="the SN adjoint spectrum does not solve the daggered "
            "eigenproblem Aᵀψ* = Fᵀψ*/k at the closed-form (A, F) — the "
            "defining-law residual is nonzero.",
        )

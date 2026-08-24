r"""#276 A4 — the φ* certification battery: P1.3 / P1.4 / P1.5.

The spec's Phase-1 chain (``p6_adjoint_verification_spec.md`` §3, with the
P1.4 formula correction recorded there): P1.0a/b + P1.1 landed pre-A4;
P1.2's duality core + lift pin live in ``test_sn_adjoint_entries.py``.
This file carries the remaining certification rows AND the named
mutations — every mutation in-process (monkeypatch; NEVER ``git
checkout``), each reddening its named gate.

* **P1.3** — ``k_adj == k_fwd`` (exact algebraic identity —
  ``eig(M†) = eig(M)`` on the iteration operator ``M = A_loss⁻¹F``,
  every factor daggered): ∞-medium 2G AND 4G (both == the closed-form
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
  plus the F†-swap mutation reddening the spectrum row — the vector
  gate's own teeth, independent of the k shift the P1.3 tooth pins on
  the same fixture.  (The Mode-12 designed-green class for k is the
  factor-ORDER/similarity family — ``eig(Mᵀ) = eig(M)``, the trap the
  P1.4 reference itself fell into — NOT leaf-transpose drops, which do
  exit the stabiliser and shift k.)
* **The A5 sphere φ*-shape row** — the coupled daggered VECTOR
  (System-A bulk + trace + System-B ray) solves the defining law
  ``A_lossᵀ(G·ψ*) = Fᵀ(G·ψ*)/k`` against a dense FORWARD-probed
  ``(A_loss, F)`` and the raw-data coupled G, on a HETEROGENEOUS
  2-region sphere (anti-vacuity: the dense pair reproduces the forward
  k first).  Teeth: the ``G_sd = V_cell → 1`` metric drop — k EXACTLY
  blind (metric-similar pencil), residual O(1) red (the ERR-067 family
  catcher) — plus the k-visible F†→F corroboration.
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

import functools
from typing import Any

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
from orpheus.sn.solver import (
    _adjoint_posing_parts,
    _as_sn_mesh,
    solve_sn,
    solve_sn_adjoint,
)
from tests.sn._test_helpers import energy_spectrum, g_coupled_diagonal



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


# ═══════════════════════════════════════════════════════════════════════
# P1.3 — the k triple-equality legs.
# ═══════════════════════════════════════════════════════════════════════


class TestP13KEquality:
    pytestmark = pytest.mark.l1

    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    @pytest.mark.verifies("kinf-mg")
    @pytest.mark.verifies("sn-adjoint-eigenproblem")
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

    @pytest.mark.verifies("sn-adjoint-eigenproblem")
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

    @pytest.mark.verifies("sn-adjoint-eigenproblem")
    def test_sphere_k_equality(self):
        r"""The COUPLED daggered posing (System B live): μ-reversal at the
        pole + the fission ray fold + the transposed block substitution.

        HONEST SCOPE: k-equality only on THIS row — the coupled daggered
        flux SHAPE is certified by :class:`TestP14SphereAdjointVector`
        (the A5 carve-in, landed): the defining-law residual against a
        dense forward-probed reference on the heterogeneous sphere.
        """
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

    pytestmark = pytest.mark.l1

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
        from orpheus.sn.operators.streaming import StreamingCollisionOperator

        mats, mesh = _het_slab()
        k_fwd = _k(_solve_fwd(mats, mesh))

        monkeypatch.setattr(
            StreamingCollisionOperator, "solve_transpose", StreamingCollisionOperator.solve,
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
    pytestmark = pytest.mark.l1

    @pytest.mark.verifies("sn-adjoint-eigenproblem")
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
            energy_spectrum(adj), psi_star_cf, rtol=1e-7,
            err_msg="the 4G SN adjoint spectrum does not match the "
            "closed-form (Aᵀ)⁻¹Fᵀ eigenvector — the energy-adjoint "
            "content (F†/S† composition) is wrong even though k agrees.",
        )

    def test_f_role_swap_reds_spectrum(self, monkeypatch):
        r"""The VECTOR row's own F† teeth: F† → F moves the adjoint
        spectrum O(1) off the closed form — evidence INDEPENDENT of the
        k shift (``test_fission_role_swap_shifts_k`` pins that half on
        the same fixture: a one-factor transpose drop is NOT a pencil
        similarity, so k moves too — 1.488 → 0.171 here, the SN-solve measurement; the 0-D closed-form proxy gives 0.153).

        The DESIGNED-GREEN k-blindness (Mode 12) lives elsewhere: on
        the factor-ORDER family — every product similar to
        ``(Aᵀ)⁻¹Fᵀ`` (``FᵀA⁻ᵀ = Mᵀ``, …) shares the k spectrum EXACTLY
        while its dominant eigenvector degenerates (rank-1 F: to ν̂Σf,
        zero A-physics).  No k row can ever catch THAT class; the
        spectrum match (``test_4g_spectrum_matches_closed_form``) + the
        reference trap-catcher are the committed catchers.
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
        spec_mut = energy_spectrum(mut)
        require(
            not np.allclose(spec_mut, psi_star_cf, rtol=1e-3),
            "F†=F did NOT move the adjoint spectrum — the P1.4 flux-shape "
            "row has no teeth on the canonical adjoint-fission trap.",
        )

    # Mode-12 accounting for the k rows (measured, 2026-07-25 sweep —
    # correcting the earlier "pencil similarity" gloss): leaf-transpose
    # DROPS (F†→F, S†→S, L†→L) are NOT k-invariant — transposing one
    # factor is not a similarity of the pencil (F†=F moves k
    # 1.488 → 0.171 on this very fixture (SN-solve measured; the 0-D
    # closed-form proxy of the same mutation gives 0.153 — cite the
    # solve, not the proxy); TestP13Mutations asserts the
    # shifts, with the regime preconditions that keep them visible).
    # The class k IS exactly blind to is the factor-ORDER/similarity
    # family (eig(Mᵀ) = eig(M)) — the trap the P1.4 reference itself
    # fell into — plus everything about the VECTOR: those are covered
    # only by the spectrum rows here and the defining-law residual
    # (TestP15).  k rows are teeth for drops, never for order/shape.


# ═══════════════════════════════════════════════════════════════════════
# The A5 sphere φ*-shape row — the coupled daggered VECTOR's independent
# reference (the carve-in the sphere k row's honest-scope note names).
# ═══════════════════════════════════════════════════════════════════════


def _het_sphere():
    r"""2-region HETEROGENEOUS sphere (A core + B shell, 2G asymmetric SigS).

    The homogeneous ``_sphere()`` weakens a vector row: forward and
    adjoint SPATIAL shapes coincide on a bare homogeneous sphere, so the
    daggered-vector materiality collapses toward angular-only.  Two
    regions make φ* differ from φ spatially AND angularly.
    """
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 3.0, 7),
        mat_ids=np.array([0, 0, 0, 1, 1, 1]),
        coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"),
    )
    return mats, mesh


def _sphere_daggered_run(sn_mesh):
    r"""Mirror :func:`solve_sn_adjoint`'s daggered chain, RAW.

    Returns ``(k_adj, psi_star)`` with ``psi_star`` the converged coupled
    ITERATE (the packaged entry cell-averages and w-reduces — this row
    needs the full System-A ⊕ System-B flat vector).  Tight knobs per the
    A5 memo §2.5 so the defining-law residual floor sits ~1e-9, two
    orders below the 1e-7 gate.
    """
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.numerics.iteration import KEigenvalue

    parts: tuple[Any, Any, Any, Any] = _adjoint_posing_parts(sn_mesh, 0)
    implicit_operator, gain, F_posed, template = parts
    if not isinstance(template, CoupledField):
        pytest.fail(
            "carrying sphere expected — this row's domain is the "
            "coupled System-A ⊕ System-B space."
        )
    ke = KEigenvalue(
        implicit_operator.H, gain.H, F_posed.H,
        max_outer=800, keff_tol=1e-10, flux_tol=1e-9,
        max_inner=800, inner_tol=1e-11,
    )
    guess = CoupledField.from_flat(
        np.ones(template.to_flat().size), template,
    )
    _out = ke.solve(initial_guess=guess)
    k_adj, psi_star = _out.keff, _out.flux_distribution
    return float(k_adj), psi_star


@functools.lru_cache(maxsize=1)
def _sphere_dense_reference():
    r"""Dense FORWARD-probed ``(A_loss, F)`` + the raw-data coupled G.

    Three structurally-independent grounds, none touching the
    ``.H``/reverse-scan machinery under test (L11 / anti-R1):
    unit-vector columns through the FORWARD ``apply`` of the posing
    parts (``A_loss = implicit_operator − gain``), numpy ``.T`` for the
    transpose, and :func:`g_coupled_diagonal`'s raw-data metric.
    Cached once — the mutation teeth judge their MUTATED iterates
    against this same unmutated reference.
    """
    from orpheus.numerics.coupled_system import CoupledField

    mats, mesh = _het_sphere()
    sn = _as_sn_mesh(mesh, _quad(), mats)
    parts: tuple[Any, Any, Any, Any] = _adjoint_posing_parts(sn, 0)
    implicit_operator, gain, F_posed, template = parts
    if not isinstance(template, CoupledField):
        pytest.fail(
            "carrying sphere expected — the dense probe spans the "
            "coupled flat space."
        )
    n = int(template.to_flat().size)
    R = np.empty((n, n))
    N_gain = np.empty((n, n))
    F_dense = np.empty((n, n))
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1.0
        x = CoupledField.from_flat(e, template)
        R[:, i] = np.asarray(implicit_operator.apply(x).to_flat(), dtype=float)
        N_gain[:, i] = np.asarray(gain.apply(x).to_flat(), dtype=float)
        F_dense[:, i] = np.asarray(F_posed.apply(x).to_flat(), dtype=float)
    return R - N_gain, F_dense, g_coupled_diagonal(sn)


def _defining_law_relative_residual(psi_star, k_adj, A_dense, F_dense, g_diag):
    r"""``‖A_lossᵀ(G·ψ*) − Fᵀ(G·ψ*)/k‖ / ‖Fᵀ(G·ψ*)/k‖``.

    ``.H = G⁻¹AᵀG`` (operator.py), so ``v := G·ψ*`` solves the
    forward-Euclidean-transpose adjoint eigenproblem
    ``A_lossᵀ v = Fᵀ v / k`` — the residual is normalization- and
    sign-free and never inverts the (grazing-singular) trace metric.
    """
    v = g_diag * np.asarray(psi_star.to_flat(), dtype=float)
    rhs = (F_dense.T @ v) / k_adj
    return float(np.linalg.norm(A_dense.T @ v - rhs) / np.linalg.norm(rhs))


class TestP14SphereAdjointVector:
    r"""The coupled sphere daggered VECTOR against an independent dense
    reference — closing the sphere k row's honest-scope gap (A5 carve-in).

    Mode-12: ``A.H = G⁻¹AᵀG`` is metric-similar to ``Aᵀ`` for ANY
    invertible G, so k is EXACTLY blind to a metric error — the G-metric
    is a free parameter of the daggered posing no eigenvalue gate can
    ever see.  THIS row is the sole catcher for a coupled metric bug;
    its primary tooth (the ``G_sd`` drop) demonstrates the blindness
    and the catch in one test.
    """

    pytestmark = pytest.mark.l1

    @pytest.mark.verifies("sn-adjoint-eigenproblem")
    def test_coupled_defining_law_residual(self):
        r"""``A_lossᵀ(G·ψ*) = Fᵀ(G·ψ*)/k`` on the heterogeneous sphere —
        the FULL coupled vector (System-A bulk + trace + System-B ray),
        with the dense probe anti-vacuity-pinned to the forward physics
        and the φ* ≠ φ materiality asserted."""
        A_dense, F_dense, g_diag = _sphere_dense_reference()
        mats, mesh = _het_sphere()

        # Anti-vacuity: the dense forward pair reproduces the forward SN
        # eigenvalue BEFORE its transpose is trusted as the reference.
        fwd = _solve_fwd(mats, mesh)
        k_dense = float(
            np.max(np.abs(np.linalg.eigvals(np.linalg.solve(A_dense, F_dense))))
        )
        np.testing.assert_allclose(
            k_dense, _k(fwd), rtol=0, atol=1e-7,
            err_msg="the dense forward-probed (A_loss, F) does not "
            "reproduce the forward SN k — the probe/chart is broken, so "
            "its transpose cannot referee the adjoint.",
        )

        sn = _as_sn_mesh(mesh, _quad(), mats)
        k_adj, psi_star = _sphere_daggered_run(sn)

        # Materiality: φ* ≢ φ (normalized) — heterogeneous + vacuum makes
        # the daggered vector genuinely non-trivial (else this row is a
        # self-adjoint dud, the P1.4 precondition discipline).
        w_n = np.asarray(sn.quad.weights, dtype=float)
        psi_a = np.asarray(psi_star.systems[0].interior.values, dtype=float)
        phi_star = np.einsum("n...,n->...", psi_a, w_n)
        phi = np.asarray(fwd.scalar_flux.values, dtype=float)
        require(
            not np.allclose(
                phi / np.linalg.norm(phi),
                phi_star / np.linalg.norm(phi_star),
                atol=1e-3,
            ),
            "φ* ≈ φ on the heterogeneous vacuum sphere — the fixture "
            "lost its adjoint materiality and this row is a dud.",
        )

        rel = _defining_law_relative_residual(
            psi_star, k_adj, A_dense, F_dense, g_diag,
        )
        np.testing.assert_allclose(
            rel, 0.0, atol=1e-7,
            err_msg="the coupled sphere daggered vector does not solve "
            "the true forward-transpose adjoint eigenproblem "
            "A_lossᵀ(Gψ*) = Fᵀ(Gψ*)/k under the raw-data G — the "
            "defining-law residual is above the iteration floor.",
        )

    def test_gsd_metric_drop_is_k_blind_but_vector_red(self):
        r"""THE Mode-12 tooth: ``G_sd = V_cell → 1`` leaves k EXACTLY
        equal (metric-similar pencil — every k row is designed-green on
        a metric error) while the residual reds O(1).  ERR-067 family:
        the ray state-metric is sphere-specific System-B content only
        this vector row can see."""
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )

        A_dense, F_dense, g_diag = _sphere_dense_reference()
        mats, mesh = _het_sphere()
        k_fwd = _k(_solve_fwd(mats, mesh))

        sn_mut = _as_sn_mesh(mesh, _quad(), mats)
        ii = sn_mut.radial_characteristic_interior_space
        bb = sn_mut.radial_characteristic_boundary_space
        cspace = sn_mut.radial_characteristic_field_space
        if ii is None or bb is None or cspace is None:
            pytest.fail("the carrying sphere must expose the ray spaces.")
        iw0 = ii.inner_product_weights
        bw0 = bb.inner_product_weights
        if iw0 is None or bw0 is None:
            pytest.fail(
                "the carrying sphere must expose explicit ray metric "
                "weights (None selects Euclidean — the ghost family)."
            )
        # Propagation probe: the mutation must reach the metric surface
        # the daggered .H reads at solve time (frozen dataclass →
        # object.__setattr__; a construction-time snapshot would leave
        # `after == before` and fail here rather than silently pass).
        tmpl = RadialCharacteristicField.flux_zeros(sn_mut.radial_characteristic_field_space)
        ones_ray = RadialCharacteristicField.from_flat(
            np.ones(int(cspace.shape[0])), tmpl,
        )
        before = float(cspace.inner_product(ones_ray, ones_ray))
        object.__setattr__(
            ii, "inner_product_weights",
            np.ones_like(np.asarray(iw0, dtype=float)),
        )
        object.__setattr__(
            bb, "inner_product_weights",
            np.ones_like(np.asarray(bw0, dtype=float)),
        )
        after = float(cspace.inner_product(ones_ray, ones_ray))
        require(
            abs(after - before) > 1e-12 * max(abs(before), 1.0),
            "the G_sd mutation did not reach the ray metric surface — "
            "the tooth would be vacuous (construction-time snapshot?).",
        )

        k_mut, psi_mut = _sphere_daggered_run(sn_mut)
        np.testing.assert_allclose(
            k_mut, k_fwd, rtol=0, atol=1e-8,
            err_msg="k moved under an INVERTIBLE metric mutation — "
            "G'⁻¹AᵀG' is similar to Aᵀ for any invertible G', so a k "
            "shift here means the mutation broke something beyond the "
            "metric (the k-blindness demonstration is invalid).",
        )
        rel_mut = _defining_law_relative_residual(
            psi_mut, k_mut, A_dense, F_dense, g_diag,
        )
        require(
            rel_mut > 1e-3,
            f"the G_sd-dropped daggered vector still satisfies the true "
            f"defining law (rel={rel_mut:.2e}) — the φ*-shape row has no "
            f"metric teeth and cannot catch the ERR-067 family.",
        )

    def test_fission_transpose_drop_reds_the_residual(self, monkeypatch):
        r"""Corroborating Mode-10 tooth (k-VISIBLE, unlike the G_sd drop):
        F† → F corrupts the daggered vector and the residual reds — the
        row catches leaf-transpose drops too, independent of the k shift
        the P1.3 teeth pin on other fixtures."""
        from orpheus.transport.operators.fission import FissionOperator

        A_dense, F_dense, g_diag = _sphere_dense_reference()
        mats, mesh = _het_sphere()

        monkeypatch.setattr(
            FissionOperator, "apply_transpose",
            lambda self, chi: self.apply(chi),
        )
        sn_mut = _as_sn_mesh(mesh, _quad(), mats)
        k_mut, psi_mut = _sphere_daggered_run(sn_mut)
        rel_mut = _defining_law_relative_residual(
            psi_mut, k_mut, A_dense, F_dense, g_diag,
        )
        require(
            rel_mut > 1e-3,
            f"F†=F left the coupled defining-law residual at the floor "
            f"(rel={rel_mut:.2e}) — the sphere vector row has no "
            f"fission-transpose teeth.",
        )


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

    @pytest.mark.l1
    @pytest.mark.verifies("sn-adjoint-eigenproblem")
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
        spec = energy_spectrum(adj)
        k = _k(adj)
        residual = A.T @ spec - (F.T @ spec) / k
        np.testing.assert_allclose(
            residual, np.zeros_like(residual), atol=1e-7,
            err_msg="the SN adjoint spectrum does not solve the daggered "
            "eigenproblem Aᵀψ* = Fᵀψ*/k at the closed-form (A, F) — the "
            "defining-law residual is nonzero.",
        )

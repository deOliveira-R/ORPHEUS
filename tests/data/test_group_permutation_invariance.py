r"""L2 — Energy-group permutation invariance of the multigroup transport solve.

What this gate proves (and why it exists)
==========================================

ORPHEUS is flipping its energy-group-index convention so that **group 0 =
FASTEST** (highest energy), ``eg`` DESCENDING, and downscatter
``SigS[g_from, g_to]`` lives in the UPPER triangle. The production GENDF
cross-section DATA today is the *opposite* (group 0 = thermal, ``eg``
ascending, downscatter LOWER triangle).

A 3-explorer audit established that **no ORPHEUS solver runs an
energy-group Gauss–Seidel sweep**: every multigroup solve builds its
in-scatter source as the FULL contraction
:math:`\sum_{g'} \Sigma_s[g', g]\,\phi_{g'}` (``SigS.T @ phi``) and its
removal as ``diag(SigT) - SigS.T - 2 Sig2.T``. Both are
**permutation-invariant** in the group index — relabelling the groups
(reverse the order) and relabelling them back leaves the converged answer
unchanged up to that relabelling. The data flip is therefore
**physics-identical**: k-eff is invariant and the flux is merely reversed
along the group axis. It is a *re-baseline*, not a correctness fix.

This file is the **empirical proof of that invariance on the CURRENT
data convention**, written BEFORE the data flip lands. If any of these
tests FAIL on current data, the "order-agnostic" premise is FALSE — a
hidden energy-order dependency (a latent bug) exists, and the flip is
NOT safe. A green run here is the load-bearing evidence that the flip is
a pure relabelling.

The group-axis crosswalk (the exact reversal)
=============================================

A :class:`~orpheus.data.macro_xs.mixture.Mixture` group-reversal flips,
along the group axis:

* **Vectors** ``(ng,)`` — reverse the single axis::

      SigT, SigC, SigF, SigP, SigL (= SigL here), nubar(folded into SigP), chi

* **Matrices** ``(ng, ng) = [g_from, g_to]`` — reverse **BOTH** axes::

      SigS[ell] (every Legendre order), Sig2

* ``eg`` ``(ng+1,)`` — reverse (energy boundaries; here synthetic XS carry
  ``eg=None`` so the field is left ``None``).

* NOT reversed: ``sig0`` (background XS, not an energy-indexed quantity).
  The synthetic :class:`Mixture` carries no ``sig0`` field (that lives on
  :class:`~orpheus.data.micro_xs.isotope.Isotope`), so there is nothing to
  guard here — confirmed by reading the dataclass.

The :func:`reverse_mixture` helper below encodes exactly this crosswalk.
It is **TEST-LOCAL** — production code is NOT modified by this gate.

Claim layer / pillar (vv-principles §1.5)
=========================================

The claim is a **symmetry / permutation-invariance** claim, NOT an
eigenvalue-*correctness* claim against an external truth value. The
"reference" for each row is the SAME solver applied to the SAME problem,
transformed by the known, exact group-reversal relabelling. This is
structurally sound because the transformation is a pure index
permutation — it is NOT another ORPHEUS solver and NOT a
procedurally-different derivation of the same identity. The symmetry is
exact in real arithmetic; the only drift is IEEE-754 non-associativity in
the reduction (well below any solver tolerance).

To *additionally* terminate the homogeneous leg in a
structurally-independent ground, the ``TestHomogeneousPermutation`` rows
cross-check the production power-iteration solver
(:func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`) against
the closed-form ``k = lambda_max(A^{-1} F)`` reference
(:func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`)
under the SAME reversal — two independent solvers, one symmetry.

Config-blindness discipline (vv-principles §0.6, L1)
====================================================

* **≥2 groups everywhere** — 1G is degenerate (``k = nuSigma_f/Sigma_a`` is
  flux-shape-independent) AND the reversal is trivial (one group does not
  reorder). Both 2G and 4G are exercised so a 2G-specific coincidence
  cannot pass silently.
* **Heterogeneous** for the SN rows — fuel (region A, fissile) + moderator
  (region B, non-fissile) so the spatial redistribution term is live and
  the reversal threads through a per-cell, per-material loop crossing a
  material boundary.
* **Asymmetric scattering** — region A 2G has ``SigS = [[0.38, 0.10],
  [0.00, 0.90]]`` (off-diagonal 0.10 vs 0.00), so the reversal is a
  *non-trivial* permutation (a symmetric matrix would make the reversal a
  no-op and the test vacuous). :func:`test_reversal_is_nontrivial` pins
  that the twin XS actually differ from the original.

The synthetic ``make_mixture`` / Sood XS are ALREADY group-0-fast and
carry no real energy grid; this gate proves the relabelling symmetry the
real-data flip relies on, using those synthetic mixtures as the vehicle.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.eigenvalue import kinf_and_spectrum_homogeneous
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.homogeneous.solver import solve_homogeneous_infinite
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn, solve_sn_fixed_source

pytestmark = [pytest.mark.l2]


# ─────────────────────────────────────────────────────────────────────
# The crosswalk helper — TEST-LOCAL group reversal (production untouched)
# ─────────────────────────────────────────────────────────────────────


def _rev_vec(v: np.ndarray) -> np.ndarray:
    """Reverse a ``(ng,)`` group-indexed vector along its single axis."""
    return np.asarray(v)[::-1].copy()


def _rev_mat(mat: csr_matrix) -> csr_matrix:
    r"""Reverse a ``(ng, ng) = [g_from, g_to]`` matrix along BOTH axes.

    Reversing both axes maps ``M[i, j] -> M[ng-1-i, ng-1-j]``, i.e. the
    scattering ``g_from -> g_to`` becomes ``(ng-1-g_from) -> (ng-1-g_to)``
    under the group relabelling ``g |-> ng-1-g``. This is the matrix half
    of the crosswalk; a single-axis reversal here would be the canonical
    convention-drift bug (Mode 6).
    """
    dense = np.asarray(mat.todense())
    return csr_matrix(dense[::-1, ::-1])


def reverse_mixture(mix: Mixture) -> Mixture:
    """Return the group-reversed twin of ``mix`` per the crosswalk.

    Reverses every group-indexed array along the group axis: vectors on
    their one axis, ``SigS[ell]`` and ``Sig2`` on BOTH axes, ``eg`` on its
    one axis (here ``None`` for synthetic XS). Leaves ``sig0`` untouched
    (not a field on :class:`Mixture`). TEST-LOCAL — does not mutate
    production code or data.
    """
    return replace(
        mix,
        SigC=_rev_vec(mix.SigC),
        SigL=_rev_vec(mix.SigL),
        SigF=_rev_vec(mix.SigF),
        SigP=_rev_vec(mix.SigP),
        SigT=_rev_vec(mix.SigT),
        SigS=[_rev_mat(s) for s in mix.SigS],
        Sig2=_rev_mat(mix.Sig2),
        chi=_rev_vec(np.asarray(mix.chi)),
        eg=(_rev_vec(mix.eg) if mix.eg is not None else None),
    )


# ─────────────────────────────────────────────────────────────────────
# Mesh / problem builders
# ─────────────────────────────────────────────────────────────────────


def _heterogeneous_slab(n_cells_per_zone: int = 4, half_width: float = 2.0) -> Mesh1D:
    """Two-zone Cartesian slab: fuel [0, w], moderator [w, 2w].

    The material boundary at ``x = w`` is what makes the reversal thread
    through a per-cell, per-material in-scatter loop (heterogeneous).
    """
    n = n_cells_per_zone
    edges = np.linspace(0.0, 2.0 * half_width, 2 * n + 1)
    mat_ids = np.array([0] * n + [1] * n, dtype=int)
    return Mesh1D(edges=edges, mat_ids=mat_ids, coord=CoordSystem.CARTESIAN)


# ─────────────────────────────────────────────────────────────────────
# Solver 1 — SN eigenvalue (heterogeneous, the spatial-redistribution path)
# ─────────────────────────────────────────────────────────────────────


class TestSNEigenvaluePermutation:
    """k-eff invariant + per-cell flux reversed under group reversal (SN)."""

    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    def test_sn_heterogeneous_eigenvalue_group_reversal(self, ng_key: str) -> None:
        r"""SN k-eff is invariant and ``phi`` reverses along the group axis.

        Heterogeneous fuel+moderator slab. Solves the eigenvalue problem
        on the original materials and on the group-reversed twin, then
        asserts:

        * ``keff_twin == keff`` to the eigenvalue's own tolerance, AND
        * ``phi_twin[g, cell] == phi[NG-1-g, cell]`` per group, per cell
          (the group axis reverses; the spatial axis is fixed).

        Catches any latent energy-ORDER dependence in the SN in-scatter
        contraction, removal-operator assembly, or fission-source build.
        A green result is the load-bearing evidence that the data flip
        leaves SN physics-identical (k invariant, flux reversed).
        """
        mat_a = get_mixture("A", ng_key)  # fissile fuel
        mat_b = get_mixture("B", ng_key)  # non-fissile moderator
        materials = {0: mat_a, 1: mat_b}
        materials_rev = {0: reverse_mixture(mat_a), 1: reverse_mixture(mat_b)}

        mesh = _heterogeneous_slab()
        quad = Quadrature.gauss_legendre(n_ordinates=8)

        sol = solve_sn(materials, mesh, quad, scattering_order=0)
        sol_rev = solve_sn(materials_rev, mesh, quad, scattering_order=0)

        assert sol.keff is not None and sol_rev.keff is not None

        # k-eff invariance: tolerance = the solver's own eigenvalue
        # convergence (keff_tol=1e-7 default) with safety headroom. The
        # observed drift is pure IEEE-754 non-associativity (~1e-15).
        np.testing.assert_allclose(
            sol_rev.keff, sol.keff, rtol=1e-6, atol=0.0,
            err_msg=(
                f"[{ng_key}] SN keff NOT group-permutation-invariant: "
                f"keff={sol.keff:.12f} vs reversed-twin={sol_rev.keff:.12f} "
                f"(diff={abs(sol.keff - sol_rev.keff):.3e}). A non-zero "
                "diff means a hidden energy-ORDER dependency — the data "
                "flip would NOT be physics-identical."
            ),
        )

        phi = sol.scalar_flux.values        # (ng, n_cells)
        phi_rev = sol_rev.scalar_flux.values
        assert phi.shape == phi_rev.shape

        # Eigenvectors carry an arbitrary scale; normalise both before the
        # shape comparison so we test the SHAPE reversal, not the scale.
        phi_n = phi / np.linalg.norm(phi)
        phi_rev_n = phi_rev / np.linalg.norm(phi_rev)
        np.testing.assert_allclose(
            phi_rev_n, phi_n[::-1, :], rtol=0.0, atol=1e-10,
            err_msg=(
                f"[{ng_key}] SN flux NOT reversed under group reversal: "
                "phi_twin[g, cell] should equal phi[NG-1-g, cell]. Max "
                f"abs diff (normalised) = "
                f"{np.max(np.abs(phi_rev_n - phi_n[::-1, :])):.3e}."
            ),
        )


# ─────────────────────────────────────────────────────────────────────
# Solver 2 — SN fixed-source (the flux-DISTRIBUTION reversal, no eigenvalue)
# ─────────────────────────────────────────────────────────────────────


class TestSNFixedSourcePermutation:
    """Fixed-source flux distribution reverses under group reversal (SN).

    Fixed-source isolates the flux SHAPE (no eigenvalue iteration). A
    uniform isotropic source is itself group-permutation-symmetric, so the
    converged flux must reverse exactly when the materials reverse. This is
    the fixed-source companion to the eigenvalue row — it pins the
    distribution, not just the dominant ratio.
    """

    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    def test_sn_fixed_source_flux_reversal(self, ng_key: str) -> None:
        r"""``phi`` from a uniform source reverses along the group axis.

        Homogeneous moderator slab (non-fissile — fixed-source is
        well-posed without an eigenvalue; an eigen-solve on a non-fissile
        mixture would be ``k = 0/abs -> nan``, a dead test — see L7). A
        uniform per-group isotropic source ``Q`` is the SAME for original
        and reversed groups (flat in g), so the converged flux must
        satisfy ``phi_twin[g] == phi[NG-1-g]`` per cell.
        """
        mat_b = get_mixture("B", ng_key)  # non-fissile scatterer
        materials = {0: mat_b}
        materials_rev = {0: reverse_mixture(mat_b)}

        ng = mat_b.ng
        # Homogeneous slab, vacuum BC (default). Uniform isotropic source,
        # group-flat so it is invariant under the group reversal.
        mesh = Mesh1D(
            edges=np.linspace(0.0, 4.0, 9),
            mat_ids=np.zeros(8, dtype=int),
            coord=CoordSystem.CARTESIAN,
        )
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        n_cells = 8
        # external_source shape for a 1-D mesh: (N, ng, *spatial) =
        # (N, ng, nx). Build a group-flat, ordinate-flat, cell-flat unit
        # source so the ONLY thing that differs between the two solves is
        # the material group order.
        n_ord = quad.n_ordinates
        q = np.ones((n_ord, ng, n_cells))

        sol = solve_sn_fixed_source(materials, mesh, quad, external_source=q)
        sol_rev = solve_sn_fixed_source(
            materials_rev, mesh, quad, external_source=q
        )

        phi = sol.scalar_flux.values        # (ng, n_cells)
        phi_rev = sol_rev.scalar_flux.values

        np.testing.assert_allclose(
            phi_rev, phi[::-1, :], rtol=1e-9, atol=1e-12,
            err_msg=(
                f"[{ng_key}] SN fixed-source flux NOT reversed under group "
                "reversal: phi_twin[g, cell] should equal phi[NG-1-g, "
                f"cell]. Max abs diff = {np.max(np.abs(phi_rev - phi[::-1, :])):.3e}."
            ),
        )


# ─────────────────────────────────────────────────────────────────────
# Solver 3 — homogeneous infinite medium (+ closed-form independent anchor)
# ─────────────────────────────────────────────────────────────────────


class TestHomogeneousPermutation:
    """k_inf invariant + spectrum reversed; cross-checked vs closed-form."""

    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    def test_homogeneous_kinf_and_spectrum_group_reversal(self, ng_key: str) -> None:
        r"""``k_inf`` invariant and the flux spectrum reverses under reversal.

        Uses the production power-iteration solver
        :func:`solve_homogeneous_infinite` on the fissile fuel mixture and
        its group-reversed twin, asserting ``k_inf`` invariance and
        ``phi_twin == phi[::-1]``.
        """
        mat_a = get_mixture("A", ng_key)
        res = solve_homogeneous_infinite(mat_a)
        res_rev = solve_homogeneous_infinite(reverse_mixture(mat_a))

        np.testing.assert_allclose(
            res_rev.k_inf, res.k_inf, rtol=1e-9, atol=0.0,
            err_msg=(
                f"[{ng_key}] homogeneous k_inf NOT group-permutation-"
                f"invariant: {res.k_inf:.12f} vs {res_rev.k_inf:.12f}."
            ),
        )

        phi = np.asarray(res.flux)
        phi_rev = np.asarray(res_rev.flux)
        phi_n = phi / np.linalg.norm(phi)
        phi_rev_n = phi_rev / np.linalg.norm(phi_rev)
        np.testing.assert_allclose(
            phi_rev_n, phi_n[::-1], rtol=0.0, atol=1e-10,
            err_msg=(
                f"[{ng_key}] homogeneous spectrum NOT reversed: "
                f"max abs diff (normalised) = "
                f"{np.max(np.abs(phi_rev_n - phi_n[::-1])):.3e}."
            ),
        )

    @pytest.mark.parametrize("ng_key", ["2g", "4g"])
    def test_kinf_matches_closed_form_under_reversal(self, ng_key: str) -> None:
        r"""Closed-form ``k = lambda_max(A^{-1}F)`` confirms the homogeneous leg.

        Structural-independence anchor (vv §1.5): the closed-form solver
        builds ``A = diag(SigT) - (SigS + 2 Sig2)^T`` and ``F = chi (x)
        nuSigma_f`` directly, independent of the power-iteration. It is run
        on BOTH the original and reversed arrays and must agree with the
        production homogeneous solver under the same reversal — two
        independent solvers, one symmetry. This terminates the homogeneous
        leg in a structurally-independent ground rather than only in the
        production solver's self-symmetry.
        """
        mix = get_mixture("A", ng_key)
        sigS0 = np.asarray(mix.SigS[0].todense())
        sig2 = np.asarray(mix.Sig2.todense())
        chi = np.asarray(mix.chi)

        k_cf, phi_cf = kinf_and_spectrum_homogeneous(
            mix.SigT, sigS0, mix.SigP, chi, sig2
        )
        # Closed-form on the reversed arrays (vectors reversed on axis 0,
        # matrices on both axes — the crosswalk applied to raw arrays).
        k_cf_rev, phi_cf_rev = kinf_and_spectrum_homogeneous(
            _rev_vec(mix.SigT),
            sigS0[::-1, ::-1],
            _rev_vec(mix.SigP),
            _rev_vec(chi),
            sig2[::-1, ::-1],
        )

        # Closed-form is itself permutation-invariant.
        np.testing.assert_allclose(
            k_cf_rev, k_cf, rtol=1e-12, atol=0.0,
            err_msg=f"[{ng_key}] closed-form k NOT reversal-invariant.",
        )
        np.testing.assert_allclose(
            phi_cf_rev, phi_cf[::-1], rtol=0.0, atol=1e-10,
            err_msg=f"[{ng_key}] closed-form spectrum NOT reversed.",
        )

        # And the production homogeneous solver matches the closed-form
        # (independent solver agreement, both under the reversal).
        res_rev = solve_homogeneous_infinite(reverse_mixture(mix))
        np.testing.assert_allclose(
            res_rev.k_inf, k_cf_rev, rtol=1e-7, atol=0.0,
            err_msg=(
                f"[{ng_key}] production homogeneous k_inf disagrees with "
                f"closed-form under reversal: {res_rev.k_inf:.12f} vs "
                f"{k_cf_rev:.12f}."
            ),
        )


# ─────────────────────────────────────────────────────────────────────
# Control — the reversal is a NON-TRIVIAL permutation (not a no-op)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("ng_key", ["2g", "4g"])
def test_reversal_is_nontrivial(ng_key: str) -> None:
    r"""The group-reversed twin XS actually DIFFER from the original.

    Guards against a vacuous gate: if the scattering matrix were
    symmetric (or every vector palindromic), the reversal would be a no-op
    and every invariance test above would pass trivially. Region A has
    asymmetric scattering (``SigS[0, 1] != SigS[1, 0]``) and a
    non-palindromic ``chi`` / ``SigP``, so the reversal genuinely permutes
    the physics. This control pins that the twin is a different mixture —
    so the invariance results are meaningful, not trivial.
    """
    mix = get_mixture("A", ng_key)
    twin = reverse_mixture(mix)

    # At least one group-indexed array must differ after reversal.
    sigt_differs = not np.array_equal(mix.SigT, twin.SigT)
    sigs_differs = not np.array_equal(
        np.asarray(mix.SigS[0].todense()), np.asarray(twin.SigS[0].todense())
    )
    chi_differs = not np.array_equal(np.asarray(mix.chi), np.asarray(twin.chi))

    # ``pytest.fail`` (a function call) rather than a bare ``assert`` so the
    # control fires under ``python -O`` regardless of invocation method —
    # this is an always-on anti-vacuity sentinel (vv Mode 8 / L4).
    if not (sigt_differs or sigs_differs or chi_differs):
        pytest.fail(
            f"[{ng_key}] group reversal is a NO-OP — the invariance gate "
            "would be vacuous. Expected asymmetric scattering / "
            "non-palindromic vectors in region A."
        )
    # Specifically: the off-diagonal scattering asymmetry must survive the
    # reversal as a genuine reordering (it is the term a convention-drift
    # bug would corrupt).
    if not sigs_differs:
        pytest.fail(
            f"[{ng_key}] SigS reversal produced an identical matrix — the "
            "scattering is unexpectedly symmetric; pick an asymmetric region."
        )

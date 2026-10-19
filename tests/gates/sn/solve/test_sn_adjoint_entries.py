r"""#276 A4 — the adjoint solver entries: ``solve_sn_adjoint`` +
``solve_sn_adjoint_fixed_source``.

The entry-level gates for the daggered posing (the posing-layer siblings
live in ``tests/gates/numerics/test_iteration.py``):

* **Eigenvalue entry** — the k triple-equality (entry == forward
  ``solve_sn`` == closed-form :math:`k_\infty`), the adjoint SPECTRUM vs
  the corrected closed form (the dominant eigenvector of
  :math:`(\mathbf{A}^T)^{-1}\mathbf{F}^T` — Mode 12: k is EXACTLY blind
  to the factor-ORDER/similarity family, ``eig(Mᵀ) = eig(M)``, and
  carries no vector information; leaf-transpose DROPS do shift k — the
  certification battery's P1.3 teeth — so the spectrum row is the
  committed catcher for order/shape, k rows for drops), the
  ∞-medium ISOTROPY of ψ*
  (angle-flatness — a cheap G-consistency signal: a metric mismatch
  between the daggered leaves distorts the angular shape first), and the
  unified-``Solution`` packaging contract.
* **Fixed-source entry** — the DUALITY identity
  :math:`\langle\psi^*, q\rangle_G = \langle q^*, \psi\rangle_G` between
  an INDEPENDENT forward ``solve_sn_fixed_source`` run and the adjoint
  entry, on the P1.2 killer config: 2-material heterogeneous VACUUM slab,
  asymmetric SigS, source and detector in DIFFERENT groups AND regions
  (fast source left, thermal detector right — the config that forces the
  downscatter chain through :math:`S^T` and makes a kernel-transpose
  error un-hideable).  The detector-response side is additionally
  cross-checked against the HAND volume sum
  :math:`\sum V\,\Sigma_d\,\varphi` — pinning the entry's angle-flat
  dual lift (no weights, no :math:`1/W`) as exactly the adjoint of the
  scalar-flux extraction.
* **Refusals** — the carrying-mesh fixed-source refusal (typed, loud —
  the unexercised coupled arm ships as a refusal, not silently) and the
  detector shape validation.

The full P1.2/P1.3/P1.4/P1.5 batteries (4G, heterogeneous + sphere k
legs, F†/S†/μ mutation teeth) are the A4-3 stage; these rows pin the
ENTRIES' construction, packaging, and first-contact physics.

vv Mode-8: ``np.testing.*`` / :func:`require` only (fire under
``python -O``).
"""
from __future__ import annotations

import numpy as np

from orpheus.numerics.outcome import SourceOutcome
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_adjoint_spectrum_homogeneous,
)
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solution import AdjointSolution, Solution
from orpheus.sn.solver import (
    solve_sn,
    solve_sn_adjoint,
    solve_sn_adjoint_fixed_source,
    solve_sn_fixed_source,
)

from tests.gates.sn._test_helpers import energy_spectrum, g_inner
from orpheus.numerics import frame_faces

pytestmark = pytest.mark.l1


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _quad():
    return Quadrature.gauss_legendre(n_ordinates=8)


def _homogeneous_2g(law: BC = BC.reflective):
    """The homogeneous 2G slab [0, 5], 10 cells, ``law`` on both faces
    (reflective: the infinite medium)."""
    mix = get_mixture("A", "2g")
    mesh = Mesher(StructuredGeometry.from_homogeneous(5.0, law)).partition(
        CellsByCount.uniform_width(10),
    ).mesh
    return {0: mix}, mesh, mix


# ═══════════════════════════════════════════════════════════════════════
# The eigenvalue entry.
# ═══════════════════════════════════════════════════════════════════════


class TestSolveSnAdjoint:
    @pytest.mark.verifies("sn-adjoint-eigenproblem")
    def test_infinite_medium_k_and_spectrum(self):
        r"""k triple-equality + the adjoint spectrum + ∞-medium isotropy.

        Mode-12 discipline: ``eig(M†) = eig(M)`` (M = A_loss⁻¹F) holds for the FULL
        dagger, so the k rows gate the posing identity but say nothing
        about the vector — and they are EXACTLY blind to the
        factor-ORDER/similarity family (the eig(Mᵀ) reference trap).
        Leaf-transpose drops DO shift k (certification P1.3 teeth); the
        spectrum row is the committed vector-level catcher, and the
        isotropy row pins the angular shape the 0-D reference cannot
        see.
        """
        materials, mesh, mix = _homogeneous_2g()
        fwd = solve_sn(
            materials, mesh, _quad(),
            keff_tol=1e-9, flux_tol=1e-8, inner_tol=1e-10, max_inner=500,
        )
        adj = solve_sn_adjoint(
            materials, mesh, _quad(),
            keff_tol=1e-9, flux_tol=1e-8, inner_tol=1e-10, max_inner=500,
        )
        k_adj, k_fwd = float(adj.outcome.keff), float(fwd.outcome.keff)

        k_cf, phi_star_cf = kinf_and_adjoint_spectrum_homogeneous(
            np.asarray(mix.SigT),
            np.asarray(mix.SigS[0].todense()),
            np.asarray(mix.SigP),
            np.asarray(mix.chi),
        )
        np.testing.assert_allclose(
            k_adj, k_fwd, rtol=0, atol=1e-8,
            err_msg="k_adjoint != k_forward — the exact algebraic identity "
            "eig(M†)=eig(M) is violated by the entry's posing.",
        )
        np.testing.assert_allclose(
            k_adj, k_cf, rtol=0, atol=1e-8,
            err_msg="k_adjoint does not match the closed-form k∞ anchor.",
        )

        # The energy spectrum of ψ* vs the corrected closed form
        # (the ONE shared reduction — tests/gates/sn/_test_helpers.py).
        spec = energy_spectrum(adj)
        np.testing.assert_allclose(
            spec, phi_star_cf, rtol=1e-8,
            err_msg="the entry's adjoint energy spectrum does not match the "
            "closed-form dominant eigenvector of (Aᵀ)⁻¹Fᵀ.",
        )
        nsf_hat = np.asarray(mix.SigP) / np.linalg.norm(np.asarray(mix.SigP))
        require(
            not np.allclose(spec, nsf_hat, rtol=1e-6),
            "the entry's spectrum equals ν̂Σf — the eig(Mᵀ) factor-order "
            "degenerate regressed somewhere in the chain.",
        )

        # ∞-medium isotropy: the true flat-problem adjoint is angle-flat.
        bulk = np.asarray(adj.angular_flux.interior.values)  # (N, ng, nx)
        require(
            bool(
                np.allclose(
                    bulk, np.broadcast_to(bulk[:1], bulk.shape),
                    rtol=1e-7, atol=1e-12 * float(np.abs(bulk).max()),
                )
            ),
            "∞-medium ψ* is NOT angle-flat — an angular-metric "
            "inconsistency between the daggered leaves (the G factors "
            "V·w_n must conjugate identically through A†, S†, F†).",
        )

    def test_solution_packaging_contract(self):
        r"""The role-typed return: AdjointSolution carrier, typed fields,
        mesh identity, history population — the A5 ruling landed (the
        role is the TYPE, no longer only the entry name)."""
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.timed_full_field import TimedFullField

        materials, mesh, _ = _homogeneous_2g()
        adj = solve_sn_adjoint(materials, mesh, _quad())
        require(
            type(adj) is AdjointSolution,
            "solve_sn_adjoint must return the role-typed AdjointSolution "
            "leaf EXACTLY (the A5 carrier ruling — not the forward "
            "Solution, not the bare base).",
        )
        require(
            adj.importance is adj.scalar_flux,
            "AdjointSolution.importance must alias scalar_flux — one "
            "storage, two vocabularies (φ* IS the importance map).",
        )
        require(
            isinstance(adj.angular_flux, TimedFullField)
            and isinstance(adj.angular_flux.interior, AngularFlux),
            "adjoint Solution.angular_flux must be the TimedFullField "
            "composite with an AngularFlux interior (duality typing: the "
            "adjoint flux is a first-class importance field).",
        )
        require(
            isinstance(adj.scalar_flux, ScalarFlux),
            "adjoint Solution.scalar_flux must be a ScalarFlux (φ*, the "
            "importance map).",
        )
        require(
            adj.angular_flux.interior.space
            == adj.angular_flux.interior.space_on(adj.problem)
            and adj.scalar_flux.space == adj.scalar_flux.space_on(adj.problem),
            "Solution space-content contract broken on the adjoint entry.",
        )
        require(
            adj.record.n_iterations >= 3 and adj.record.converged,
            "the adjoint Solution's record must carry the outer iteration.",
        )
        require(
            adj.radial_characteristic is None,
            "a Cartesian slab has no System B — radial_characteristic "
            "must be None.",
        )


# ═══════════════════════════════════════════════════════════════════════
# The fixed-source entry — the duality identity on the killer config.
# ═══════════════════════════════════════════════════════════════════════


def _het_vacuum_slab():
    r"""2-material heterogeneous VACUUM slab, asymmetric SigS (A/B 2G)."""
    mats = {0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}
    mesh = Mesher(StructuredGeometry.slab(
        (0.0, 1.0, 3.0, 4.0), (0, 1, 0), left=BC.vacuum, right=BC.vacuum,
    )).partition((
        CellsByCount.uniform_width(2),
        CellsByCount.uniform_width(4),
        CellsByCount.uniform_width(2),
    )).mesh
    return mats, mesh


class TestSolveSnAdjointFixedSource:
    def test_entry_family_role_types(self):
        r"""The A5 ruling across the ENTRY FAMILY: forward entries return
        exactly ``Solution``, adjoint entries exactly ``AdjointSolution``
        — the role is the TYPE, stamped by the entry, never a field the
        caller inspects.  (The eigenvalue-adjoint leaf is pinned in
        :meth:`TestSolveSnAdjoint.test_solution_packaging_contract`;
        this row covers the remaining three on the cheap fixture: the
        fixed-source entries on the vacuum slab, the eigenvalue entry on the
        reflective one.)"""
        materials, mesh, _ = _homogeneous_2g(BC.vacuum)
        _, reflective_mesh, _ = _homogeneous_2g(BC.reflective)
        quad = _quad()
        N, ng, nx = quad.N, 2, 10

        fwd_fs = solve_sn_fixed_source(
            materials, mesh, quad, np.ones((N, ng, nx)),
        )
        require(
            type(fwd_fs) is Solution,
            "solve_sn_fixed_source must return exactly the forward "
            "Solution leaf (A5 role axis).",
        )
        adj_fs = solve_sn_adjoint_fixed_source(
            materials, mesh, quad, np.ones((ng, nx)),
        )
        require(
            type(adj_fs) is AdjointSolution,
            "solve_sn_adjoint_fixed_source must return exactly the "
            "role-typed AdjointSolution leaf (A5 carrier ruling).",
        )
        require(
            isinstance(adj_fs.outcome, SourceOutcome) and not hasattr(adj_fs, "keff"),
            "the adjoint fixed-source kind is the outcome's TYPE (step 3: a "
            "SourceOutcome, no keff attribute at all) — the two "
            "discrimination axes are independent.",
        )
        fwd_k = solve_sn(materials, reflective_mesh, quad)
        require(
            type(fwd_k) is Solution,
            "solve_sn must return exactly the forward Solution leaf "
            "(A5 role axis).",
        )

    @pytest.mark.verifies("sn-adjoint-duality")
    def test_duality_cross_group_source_detector(self):
        r"""``⟨ψ*, q⟩_G == ⟨q*, ψ⟩_G`` — the discrete duality identity,
        cross-group AND cross-region (the P1.2 killer config).

        Forward: a FAST-group source in the LEFT region.  Adjoint: a
        THERMAL-group detector in the RIGHT region.  The response must
        flow through the downscatter chain, so an S† transpose error (the
        wrong transfer direction) breaks the identity O(1) — with a
        symmetric SigS or same-group source/detector it could pass
        falsely (the config-blindness the spec §2 audit names).

        Both pairings are evaluated with the INDEPENDENT ``g_inner``
        (the reciprocity file's hand-built G — anti-R1).  The detector
        side is ALSO cross-checked against the hand volume sum
        ``Σ V·Σ_d·φ`` — pinning the entry's detector lift, the angular
        retraction's adjoint R† (the pullback, a plain broadcast).
        """
        from orpheus.sn.problem import SNProblem
        from orpheus.transport.source_sinks import (
            AngularBoundarySourceSink,
            AngularSourceSink,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        mats, mesh = _het_vacuum_slab()
        quad = _quad()
        ng, nx = 2, 8

        # Forward: fast (g=0) iso source in the LEFT region, vacuum BCs.
        q_iso = np.zeros((ng, nx))
        q_iso[0, 0:2] = 1.3
        W = float(quad.weights.sum())
        q_per_ord = np.broadcast_to(
            (q_iso / W)[None], (quad.N, ng, nx),
        ).copy()
        fwd = solve_sn_fixed_source(
            mats, mesh, quad, q_per_ord,
            inner_tol=1e-12, max_inner=2000,
        )

        # Adjoint: thermal (g=1) detector in the RIGHT region.
        sigma_d = np.zeros((ng, nx))
        sigma_d[1, 6:8] = 0.7
        adj = solve_sn_adjoint_fixed_source(
            mats, mesh, quad, sigma_d,
            inner_tol=1e-12, max_inner=2000,
        )

        # The pairings, on the SAME independent G.
        sn = adj.problem
        require(isinstance(sn, SNProblem), "entry must return its SNProblem.")
        q_composite = TimedFullField(
            interior=AngularSourceSink(values=q_per_ord, space=sn.angular_bulk_space),
            boundary=AngularBoundarySourceSink.zeros(sn.angular_trace),
            _history=(), history_depth=2,
        )
        qstar_composite = TimedFullField(
            interior=AngularSourceSink(values=np.broadcast_to(sigma_d[None], (quad.N, ng, nx)).copy(), space=sn.angular_bulk_space),
            boundary=AngularBoundarySourceSink.zeros(sn.angular_trace),
            _history=(), history_depth=2,
        )
        # ψ composites: the Solutions' angular members carry bulk + trace.
        # NOTE the forward ran on ITS OWN SNProblem instance — rebuild the
        # pairing on the adjoint's mesh via raw values (the meshes are
        # declaration-identical; g_inner reads values only).
        lhs = g_inner(adj.angular_flux, q_composite, sn)      # ⟨ψ*, q⟩_G
        rhs = g_inner(qstar_composite, fwd.angular_flux, sn)  # ⟨q*, ψ⟩_G
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-7,
            err_msg=f"duality ⟨ψ*,q⟩_G={lhs:.10e} != ⟨q*,ψ⟩_G={rhs:.10e} — "
            "the daggered fixed-source system is not the adjoint of the "
            "forward one (S†/L†/lift error).",
        )

        # The detector-response spelling: ⟨q*, ψ⟩_G == Σ V·Σd·φ (hand).
        V = np.diff(np.asarray(mesh.edges))
        phi_fwd = np.asarray(fwd.scalar_flux.values)  # (ng, nx)
        hand_response = float(np.sum(V[None, :] * sigma_d * phi_fwd))
        np.testing.assert_allclose(
            rhs, hand_response, rtol=1e-10,
            err_msg="⟨q*, ψ⟩_G != Σ V·Σd·φ — the detector lift is NOT "
            "R†, the adjoint of the scalar-flux extraction (a stray w_n "
            "or 1/W in the detector lift).",
        )

        # Anti-vacuous: the response is genuinely nonzero (the downscatter
        # chain is live) and the adjoint flux genuinely reaches the source
        # region.
        require(
            abs(hand_response) > 1e-8,
            "the cross-group response is ~0 — the config does not "
            "exercise the downscatter chain; the duality row is vacuous.",
        )

    def test_carrying_mesh_refusal_is_typed_and_loud(self):
        r"""The daggered coupled fixed-source arm ships as a REFUSAL, not
        silently unexercised (#276 A4 scope note)."""
        mats = {0: get_mixture("A", "2g")}
        sphere = Mesher(StructuredGeometry.sphere(
            (0.0, 3.0), (0,), outer=BC.vacuum,
        )).partition(CellsByCount.uniform_width(6)).mesh
        with pytest.raises(NotImplementedError, match="(?i)coupled"):
            solve_sn_adjoint_fixed_source(
                mats, sphere, Quadrature.gauss_legendre(n_ordinates=8),
                np.ones((2, 6)),
            )

    def test_detector_shape_validation(self):
        mats, mesh = _het_vacuum_slab()
        with pytest.raises(ValueError, match="detector_response shape"):
            solve_sn_adjoint_fixed_source(
                mats, mesh, _quad(), np.ones((3, 8)),
            )


# ═══════════════════════════════════════════════════════════════════════
# The adjoint at an UNDER-RESOLVED scattering order — the rung above L ≤ 1.
# ═══════════════════════════════════════════════════════════════════════
#
# Every other reciprocity or adjoint gate under ``tests/gates/sn`` runs at
# scattering order L = 0 or 1, where every shipped rule integrates the
# products of two live harmonics exactly, so the frame square M* = R/W
# closes. Production's scattering gain is standard S_N,
# S = (1/W)·R Λ M, with S.H = G⁻¹SᵀG, G = diag(V·w_n): it is self-adjoint in
# the quadrature weights on every rule and does NOT rest on the square. A
# spelling that daggers the two frame faces separately and assumes the
# square (ERR-039's class) agrees with production wherever the square
# closes, so the L ≤ 1 gates cannot see it. These rows pin the adjoint where
# the square FAILS: level_symmetric(4) at L = 3 (the Gram goes dense) and at
# L = 5 (66 harmonics on 24 ordinates: the reconstruction face has a
# kernel), with L = 2 on the same rule as the resolved control.
# `[M]` 2026-10-18, ``scratch/adjoint_underresolved/`` (diagnosis) and
# ``scratch/adjoint_rung/`` (these rows' probes and battery).

_RUNG_LMAX = 5
_RUNG_INNER_TOL = 1e-15
#: Reciprocity bound: about four decades above the measured gap at this
#: inner tolerance, so it is not the iteration exit that the row reads.
_RUNG_RECIPROCITY_BOUND = 1e-12

#: ASYMMETRIC P0 blocks (strong 0→1 downscatter, weak upscatter) and a
#: per-element anisotropy ratio, so Σ_ℓ = Σ_0 ∘ ratio^ℓ is group-asymmetric at
#: EVERY order: a dropped group transpose is visible at each ℓ, and every
#: |Σ_ℓ|/Σ_0 = |ratio|^ℓ ≤ 1 (a physical Legendre stack).
_RUNG_P0 = {
    0: np.array([[0.50, 0.30], [0.02, 0.90]]),
    1: np.array([[0.30, 0.10], [0.01, 0.40]]),
}
_RUNG_RATIO = {
    0: np.array([[0.8, 0.6], [0.5, 0.7]]),
    1: np.array([[0.6, 0.4], [0.3, 0.5]]),
}
_RUNG_SIG_C = {0: np.array([0.05, 0.10]), 1: np.array([0.20, 0.50])}

_RECIPROCITY_ID = (
    "tests/gates/sn/solve/test_sn_adjoint_entries.py::"
    "TestSolveSnAdjointFixedSource::test_duality_cross_group_source_detector"
)
_SCATTERING_ADJOINT = "tests/gates/sn/operators/test_scattering_adjoint.py::"
_FRAME_EDGE_ID = (
    "tests/gates/numerics/test_frame.py::"
    "test_the_square_fails_where_the_order_outruns_the_rule"
)
_RUNG = (
    "tests/gates/sn/solve/test_sn_adjoint_entries.py::"
    "TestAdjointAtAnUnderResolvedOrder::"
)


def _rung_mixture(material: int, *, transposed: bool = False):
    r"""The P0–P5 anisotropic 2G mixture; ``transposed`` flips every Σ_ℓ.

    The transposed stack is the STANDARD discrete-ordinates adjoint's
    data (Lewis & Miller: the same quadrature reconstruction, the transfer
    transposed in the groups), so a forward ``S`` built on it is identity
    (b)'s right-hand side, assembled with no adjoint machinery at all.
    """
    from scipy.sparse import csr_matrix

    from orpheus.data.macro_xs.mixture import Mixture

    p0, ratio, sig_c = _RUNG_P0[material], _RUNG_RATIO[material], _RUNG_SIG_C[material]
    stack = [p0 * ratio**ell for ell in range(_RUNG_LMAX + 1)]
    if transposed:
        stack = [block.T for block in stack]
    zero = np.zeros(2)
    return Mixture(
        SigC=sig_c.copy(), SigL=zero.copy(), SigF=zero.copy(), SigP=zero.copy(),
        SigT=sig_c + stack[0].sum(axis=1),
        SigS=[csr_matrix(block) for block in stack],
        Sig2=[csr_matrix((2, 2))], chi=zero.copy(),
    )


def _rung_slab(*, transposed: bool = False):
    r"""The heterogeneous vacuum slab of :func:`_het_vacuum_slab` (A | B | A,
    8 cells) carrying the anisotropic mixtures."""
    _, mesh = _het_vacuum_slab()
    mats = {m: _rung_mixture(m, transposed=transposed) for m in (0, 1)}
    return mats, mesh


def _ls4():
    return Quadrature.level_symmetric(4)


#: (order, square closes?) on level_symmetric(4); the ids say which side
#: of the frame-square edge each row sits on.
_RUNG_ROWS = [
    pytest.param(2, True, id="LS4-L2-resolved"),
    pytest.param(3, False, id="LS4-L3-under-resolved"),
    pytest.param(5, False, id="LS4-L5-non-injective"),
]


def _assert_rung_premise(quad, order: int, square_closes: bool) -> None:
    r"""The row's premise, asserted so a change of quadrature cannot silently
    move it across the edge: the square's residual
    ``max|M* − R/W| / max|R/W|`` (draw-free, both faces as matrices) is ≤
    1e-12 on the resolved row and > 1e-3 on the others; the L = 5 row
    additionally has more harmonics than ordinates. `[M]` 2026-10-18:
    residuals 0, 0.1261, 0.6180 at L = 2, 3, 5."""
    from tests.gates.numerics.test_frame import _collapse_residual

    frame = frame_faces.angular_frame(quad, order)
    residual = _collapse_residual(frame)
    if square_closes:
        require(
            residual < 1e-12,
            f"premise lost: the frame square of level_symmetric(4) at L={order} "
            f"no longer closes (residual {residual:.3e}); this row is the "
            f"RESOLVED control — re-point it to a row where the square closes.",
        )
        return
    require(
        residual > 1e-3,
        f"premise lost: the frame square of level_symmetric(4) at L={order} "
        f"now closes (residual {residual:.3e}), so this row no longer tests "
        f"the adjoint at an under-resolved order — re-point it to an order "
        f"that outruns the rule.",
    )
    if order == _RUNG_LMAX:
        n_modes = int(np.prod(frame_faces.basis_space(frame).shape))
        require(
            n_modes > quad.N,
            f"premise lost: L={order} carries {n_modes} harmonics on {quad.N} "
            f"ordinates — the reconstruction face is no longer non-injective.",
        )


def _response_pair(quad, order: int):
    r"""``(⟨q*, φ⟩, ⟨q, φ*⟩)`` from an INDEPENDENT forward and adjoint solve:
    a fast isotropic source in the left region, a thermal detector in the
    right one (the duality row's killer configuration). Both pairings are
    the hand volume sums ``Σ V·Σ_d·φ`` and ``Σ V·q·φ*/W`` (the per-ordinate
    source is ``q/W``), so no production inner product is read."""
    mats, mesh = _rung_slab()
    ng, nx = 2, 8
    W = float(quad.weights.sum())
    q_iso = np.zeros((ng, nx))
    q_iso[0, 0:2] = 1.3
    sigma_d = np.zeros((ng, nx))
    sigma_d[1, 6:8] = 0.7
    fwd = solve_sn_fixed_source(
        mats, mesh, quad,
        np.broadcast_to((q_iso / W)[None], (quad.N, ng, nx)).copy(),
        scattering_order=order, inner_tol=_RUNG_INNER_TOL, max_inner=5000,
    )
    adj = solve_sn_adjoint_fixed_source(
        mats, mesh, quad, sigma_d,
        scattering_order=order, inner_tol=_RUNG_INNER_TOL, max_inner=5000,
    )
    V = np.diff(np.asarray(mesh.edges))
    forward = float(np.sum(V * sigma_d * np.asarray(fwd.scalar_flux.values)))
    adjoint = float(np.sum(V * q_iso * np.asarray(adj.scalar_flux.values))) / W
    return forward, adjoint


def _dense_scattering(quad, order: int, *, transposed: bool = False):
    r"""``(S, S.H, G)`` as dense matrices on the rung slab — columns swept on
    unit vectors through the production bound operator, ``G = w_n ⊗ 1 ⊗ V``
    flattened."""
    from orpheus.sn.problem import SNProblem
    from orpheus.transport.fields.angular_flux import AngularFlux
    from orpheus.transport.source_sinks import AngularSourceSink
    from tests.gates.sn.operators._composite_operand import bulk_apply

    mats, mesh = _rung_slab(transposed=transposed)
    problem = SNProblem(mesh, quad, mats, scattering_order=order)
    S = problem.system.factors.scattering
    shape = (quad.N, 2, mesh.N)
    n = int(np.prod(shape))
    forward, adjoint = [], []
    for j in range(n):
        unit = np.zeros(n)
        unit[j] = 1.0
        unit = unit.reshape(shape)
        forward.append(np.asarray(bulk_apply(
            S, AngularFlux(values=unit, space=problem.angular_bulk_space),
        ).values).ravel())
        source = AngularSourceSink(values=unit, space=problem.angular_bulk_space)
        adjoint.append(np.asarray(bulk_apply(S.H, source).values).ravel())
    V = np.diff(np.asarray(mesh.edges))
    G = np.einsum("n,g,c->ngc", np.asarray(quad.weights), np.ones(2), V).ravel()
    return np.array(forward).T, np.array(adjoint).T, G


def _relative(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _install_mixed_convention(monkeypatch, quad, mesh) -> None:
    r"""ERR-039's class, in process: S.H := M* Λᵀ M — the analysis face
    daggered on its own and the reconstruction's adjoint taken as W·M (the
    square assumed). Installed as the Euclidean transpose G·(M* Λᵀ M)·G⁻¹,
    so production's ``.H = G⁻¹SᵀG`` returns exactly the mixed spelling. It
    equals the correct (1/W)·R Λᵀ M wherever the square closes."""
    from orpheus.transport.operators.transfer import TransferOperator

    V = np.diff(np.asarray(mesh.edges))
    G = np.asarray(quad.weights)[:, None, None] * V[None, None, :]

    def mixed(self, bulk):
        analysis = frame_faces.FrameAnalysis.of(self.frame)
        transfer = self._moment_transfer(skip_l0=False)
        x = np.asarray(bulk.values) / G
        y = np.asarray(analysis.H.apply(transfer.apply_transpose(np.asarray(analysis.apply(x)))))
        return self._end.cotangent(self, G * y)

    # On the transfer binding only (scattering and (n,2n)); fission is not a
    # frame-conjugated transfer and keeps its honest transpose.
    monkeypatch.setattr(TransferOperator, "_interior_transpose", mixed)


class TestAdjointAtAnUnderResolvedOrder:
    r"""The SN adjoint is exact at a scattering order the angular rule does
    not integrate — level_symmetric(4) at L = 3 and L = 5, with L = 2 on the
    same rule as the resolved control.

    Claim layer: the operator (L0, identities (a) and (b)) and the
    fixed-source flux functional (L1, reciprocity). Activated: the ℓ ≥ 2
    moments in two groups with group-asymmetric Σ_ℓ at every order, two
    materials, vacuum leakage, cross-group and cross-region response.
    Nulled: curvature and reflection (slab, vacuum faces) — the
    curvilinear rows are the cylinder's business.

    First reds (X1), `[M]` 2026-10-18, ``scratch/adjoint_rung/``:

    * the Λᵀ → Λ drop (``TransferMaterialField.moment_source_transpose``
      replaced by ``moment_source``) reddens every row, the resolved one
      included (reciprocity rel 0.92);
    * the mixed-convention spelling (:func:`_install_mixed_convention`)
      reddens the two under-resolved rows only (reciprocity rel 4.7e-4 at
      L = 3, 3.2e-3 at L = 5; the resolved row reads 7e-16) — the defect
      every L ≤ 1 gate is blind to, and the reason this rung exists.
      :meth:`test_the_rung_sees_the_mixed_convention` pins that split.
    """

    @pytest.mark.rests_on(
        _SCATTERING_ADJOINT + "TestLambdaTranspose::test_transpose_matches_dense_per_material",
        _SCATTERING_ADJOINT + "TestFullScatterKernel::test_S_euclidean_reciprocity",
        _FRAME_EDGE_ID,
    )
    @pytest.mark.parametrize("order, square_closes", _RUNG_ROWS)
    def test_the_scattering_adjoint_is_the_metric_adjoint_and_the_transposed_transfer(
        self, order, square_closes,
    ):
        r"""(a) ``S.H = G⁻¹SᵀG`` and (b) ``S.H = S[Σ_ℓ → Σ_ℓᵀ]``, as dense
        matrices.

        (a) pins the ``.H`` route (the Euclidean transpose and the Riesz
        legs) against the FORWARD apply transposed by hand. (b) pins it
        against an independent assembly: a forward operator posed on the
        transposed data, which is the standard discrete-ordinates adjoint
        and involves no adjoint machinery. Tolerance: a direct product
        whose deepest reduction sums over the ordinates and the harmonics,
        ``8·(N + n_modes)·ε`` relative (`[M]` ≤ 3e-16 on all three rows;
        the bound is 1.7e-13 at L = 5).
        """
        quad = _ls4()
        _assert_rung_premise(quad, order, square_closes)
        S, SH, G = _dense_scattering(quad, order)
        S_transposed_data, _, _ = _dense_scattering(quad, order, transposed=True)
        n_modes = int(np.prod(frame_faces.basis_space(frame_faces.angular_frame(quad, order)).shape))
        tol = 8.0 * (quad.N + n_modes) * np.finfo(float).eps
        require(
            _relative(S, S_transposed_data) > 0.1,
            "the transposed data reproduces S — the group asymmetry the "
            "identities discriminate on is gone; the row is vacuous.",
        )
        a = _relative(SH, (S.T * G[None, :]) / G[:, None])
        require(
            a <= tol,
            f"(a) S.H != G^-1 S^T G at level_symmetric(4) L={order}: rel {a:.3e} "
            f"> {tol:.1e} — the .H route is not the metric adjoint of the "
            f"forward scattering gain.",
        )
        b = _relative(SH, S_transposed_data)
        require(
            b <= tol,
            f"(b) S.H != S[Sigma_l^T] at level_symmetric(4) L={order}: rel "
            f"{b:.3e} > {tol:.1e} — the adjoint is not the standard S_N "
            f"adjoint (same reconstruction, transfer transposed in the groups).",
        )

    @pytest.mark.verifies("sn-adjoint-duality")
    @pytest.mark.rests_on(
        _RECIPROCITY_ID,
        _RUNG + "test_the_scattering_adjoint_is_the_metric_adjoint_and_the_transposed_transfer",
    )
    @pytest.mark.parametrize("order, square_closes", _RUNG_ROWS)
    def test_reciprocity_through_the_entries(self, order, square_closes):
        r"""``⟨q*, φ⟩ = ⟨q, φ*⟩`` between an independent forward
        ``solve_sn_fixed_source`` and ``solve_sn_adjoint_fixed_source``.

        Tolerance: the solves run at inner tolerance 1e-15 and the bound is
        1e-12 relative, about four decades above the gap. `[M]` 2026-10-18,
        ``PYTHONPATH=. .venv/bin/python -O scratch/adjoint_rung/q3_tol.py``:
        3.7e-16, 1.9e-16, 5.6e-16 at L = 2, 3, 5. At inner tolerance 1e-13
        the gap was 4.6e-13: the iteration exit, not the operator
        (``scratch/adjoint_underresolved/p1t_tol.py``). Activation leg: the order is live — the
        forward response at L differs from the one at L − 1 by more than
        1e-4 relative, so the top moment carries the answer.
        """
        quad = _ls4()
        _assert_rung_premise(quad, order, square_closes)
        forward, adjoint = _response_pair(quad, order)
        lower, _ = _response_pair(quad, order - 1)
        require(
            abs(forward - lower) > 1e-4 * abs(forward),
            f"the ell={order} moments do not move the response "
            f"({lower:.10e} at L={order - 1}, {forward:.10e} at L={order}); "
            f"the row does not exercise the order it names.",
        )
        rel = abs(forward - adjoint) / abs(forward)
        require(
            rel <= _RUNG_RECIPROCITY_BOUND,
            f"reciprocity broken at level_symmetric(4) L={order}: <q*,phi>="
            f"{forward:.15e}, <q,phi*>={adjoint:.15e}, rel {rel:.3e} > "
            f"{_RUNG_RECIPROCITY_BOUND:.0e} — the adjoint entry is not the "
            f"adjoint of the forward one at this scattering order.",
        )

    @pytest.mark.rests_on(_RUNG + "test_reciprocity_through_the_entries")
    def test_the_rung_sees_the_mixed_convention(self, monkeypatch):
        r"""The loaded reading (``vv-principles`` #19): under the mixed
        spelling of :func:`_install_mixed_convention`, reciprocity BREAKS at
        L = 3 and HOLDS at L = 2 on the same rule.

        The L = 2 half is the finding the rung rests on: a resolved order is
        blind to this defect, so the under-resolved rows are its only
        catcher. `[M]` 2026-10-18: rel 4.7e-4 at L = 3, 7.4e-16 at L = 2."""
        quad = _ls4()
        _, mesh = _rung_slab()
        _install_mixed_convention(monkeypatch, quad, mesh)
        f3, a3 = _response_pair(quad, 3)
        f2, a2 = _response_pair(quad, 2)
        require(
            abs(f3 - a3) / abs(f3) > 1e-5,
            f"the mixed-convention adjoint passes reciprocity at L=3 (rel "
            f"{abs(f3 - a3) / abs(f3):.3e}) — the under-resolved row is blind "
            f"to ERR-039's class.",
        )
        require(
            abs(f2 - a2) / abs(f2) <= _RUNG_RECIPROCITY_BOUND,
            f"the mixed-convention adjoint breaks reciprocity at the RESOLVED "
            f"L=2 (rel {abs(f2 - a2) / abs(f2):.3e}): the mutant is not the "
            f"square-assuming spelling, and this row's split means nothing.",
        )

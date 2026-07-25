r"""#276 A4 — the adjoint solver entries: ``solve_sn_adjoint`` +
``solve_sn_adjoint_fixed_source``.

The entry-level gates for the daggered posing (the posing-layer siblings
live in ``tests/numerics/test_iteration.py``):

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
import pytest

from orpheus.derivations.common.eigenvalue import (
    kinf_and_adjoint_spectrum_homogeneous,
)
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solution import AdjointSolution, Solution
from orpheus.sn.solver import (
    solve_sn,
    solve_sn_adjoint,
    solve_sn_adjoint_fixed_source,
    solve_sn_fixed_source,
)

from tests.sn._test_helpers import energy_spectrum, g_inner

pytestmark = pytest.mark.l1


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _quad():
    return Quadrature.gauss_legendre(n_ordinates=8)


def _homogeneous_2g():
    mix = get_mixture("A", "2g")
    mesh = Mesh1D(
        edges=np.linspace(0.0, 5.0, 11), mat_ids=np.zeros(10, dtype=int),
    )
    return {0: mix}, mesh, mix


# ═══════════════════════════════════════════════════════════════════════
# The eigenvalue entry.
# ═══════════════════════════════════════════════════════════════════════


class TestSolveSnAdjoint:
    def test_infinite_medium_k_and_spectrum(self):
        r"""k triple-equality + the adjoint spectrum + ∞-medium isotropy.

        Mode-12 discipline: ``eig(A†) = eig(A)`` holds for the FULL
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
        if adj.keff is None or fwd.keff is None:
            pytest.fail("an entry returned no eigenvalue — the reference "
                        "or adjoint leg is broken.")
        k_adj, k_fwd = float(adj.keff), float(fwd.keff)

        k_cf, phi_star_cf = kinf_and_adjoint_spectrum_homogeneous(
            np.asarray(mix.SigT),
            np.asarray(mix.SigS[0].todense()),
            np.asarray(mix.SigP),
            np.asarray(mix.chi),
        )
        np.testing.assert_allclose(
            k_adj, k_fwd, rtol=0, atol=1e-8,
            err_msg="k_adjoint != k_forward — the exact algebraic identity "
            "eig(A†)=eig(A) is violated by the entry's posing.",
        )
        np.testing.assert_allclose(
            k_adj, k_cf, rtol=0, atol=1e-8,
            err_msg="k_adjoint does not match the closed-form k∞ anchor.",
        )

        # The energy spectrum of ψ* vs the corrected closed form
        # (the ONE shared reduction — tests/sn/_test_helpers.py).
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
            adj.angular_flux.interior.mesh is adj.mesh
            and adj.scalar_flux.mesh is adj.mesh,
            "Solution mesh-identity contract broken on the adjoint entry.",
        )
        require(
            adj.history is not None
            and adj.history.n_outer is not None
            and adj.history.n_outer >= 3
            and adj.history.converged,
            "adjoint Solution.history must carry the outer trajectory.",
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
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, 9),
        mat_ids=np.array([0, 0, 1, 1, 1, 1, 0, 0]),
        bc_left=None, bc_right=None,  # vacuum via the entry default
    )
    return mats, mesh


class TestSolveSnAdjointFixedSource:
    def test_entry_family_role_types(self):
        r"""The A5 ruling across the ENTRY FAMILY: forward entries return
        exactly ``Solution``, adjoint entries exactly ``AdjointSolution``
        — the role is the TYPE, stamped by the entry, never a field the
        caller inspects.  (The eigenvalue-adjoint leaf is pinned in
        :meth:`TestSolveSnAdjoint.test_solution_packaging_contract`;
        this row covers the remaining three on the cheap fixture.)"""
        materials, mesh, _ = _homogeneous_2g()
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
            adj_fs.keff is None and adj_fs.is_fixed_source(),
            "the adjoint fixed-source kind rides the SAME problem-kind "
            "property as the forward (keff None) — the two "
            "discrimination axes are independent.",
        )
        fwd_k = solve_sn(materials, mesh, quad)
        require(
            type(fwd_k) is Solution,
            "solve_sn must return exactly the forward Solution leaf "
            "(A5 role axis).",
        )

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
        ``Σ V·Σ_d·φ`` — pinning the entry's angle-flat dual lift.
        """
        from orpheus.sn.mesh.augmented_mesh import SNMesh
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
            boundary_condition="vacuum", inner_tol=1e-12, max_inner=2000,
        )

        # Adjoint: thermal (g=1) detector in the RIGHT region.
        sigma_d = np.zeros((ng, nx))
        sigma_d[1, 6:8] = 0.7
        adj = solve_sn_adjoint_fixed_source(
            mats, mesh, quad, sigma_d,
            boundary_condition="vacuum", inner_tol=1e-12, max_inner=2000,
        )

        # The pairings, on the SAME independent G.
        sn = adj.mesh
        require(isinstance(sn, SNMesh), "entry must return its SNMesh.")
        q_composite = TimedFullField(
            interior=AngularSourceSink.from_mesh(q_per_ord, sn),
            boundary=AngularBoundarySourceSink.zeros_on(sn),
            _history=(), history_depth=2,
        )
        qstar_composite = TimedFullField(
            interior=AngularSourceSink.from_mesh(
                np.broadcast_to(sigma_d[None], (quad.N, ng, nx)).copy(), sn,
            ),
            boundary=AngularBoundarySourceSink.zeros_on(sn),
            _history=(), history_depth=2,
        )
        # ψ composites: the Solutions' angular members carry bulk + trace.
        # NOTE the forward ran on ITS OWN SNMesh instance — rebuild the
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
            err_msg="⟨q*, ψ⟩_G != Σ V·Σd·φ — the angle-flat dual lift is "
            "NOT the adjoint of the scalar-flux extraction (a stray w_n "
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
        from orpheus.geometry.mesh import CoordSystem

        mats = {0: get_mixture("A", "2g")}
        sphere = Mesh1D(
            edges=np.linspace(0.0, 3.0, 7),
            mat_ids=np.zeros(6, dtype=int),
            coord=CoordSystem.SPHERICAL,
        )
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

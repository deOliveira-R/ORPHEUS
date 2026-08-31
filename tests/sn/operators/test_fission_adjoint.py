r"""P1.0a — adjoint fission :math:`F^\dagger` correctness (campaign #276, phase A1).

The transpose of the rank-1 fission dyad :math:`F = |\chi\rangle\langle\nu\Sigma_f|`
is the **dual dyad** :math:`F^\dagger = |\nu\Sigma_f\rangle\langle\chi|` — the
**χ↔νΣf role swap**:

.. math::

    (F^\dagger\psi^*)_g(\vec r) = \nu\Sigma_{f,g}(\vec r)\,
        \sum_{g'} \chi_{g'}(\vec r)\,\psi^*_{g'}(\vec r).

Gates (vv ``foundation`` / L0):

* **CORRECTNESS** — ``F.apply_transpose`` matches a STRUCTURALLY-INDEPENDENT
  explicit Python double-loop. The reference is the SAME hand loop the forward
  gate uses (``hand_derived_fission_emission``), with its ``chi`` and
  ``nu_sigma_f`` arguments role-swapped — because F† IS the dyad swap. It shares
  NO numpy reduction with the production ``RankOneOperator.apply_transpose`` /
  ``InnerProductFunctional`` path. 2G **and** 4G (vv L2: ≥2 groups; 1G nulls the
  group transpose).
* **RECIPROCITY** — ``⟨Fφ, ψ*⟩ == ⟨φ, F†ψ*⟩`` (full Euclidean), the
  transpose-DEFINING identity. A wrong F† (wrong axis / not actually the
  transpose) violates it.
* **ROLE-SWAP DISCRIMINATOR** (the canonical adjoint-fission trap, vv Mode-2) —
  with χ ≁ νΣf per group, swapping χ↔νΣf in F† changes the emission. This is the
  mutation precondition: if the fixture lost its asymmetry the correctness gate
  would be blind to the swap.
* **CAPABILITY + ROUTING** (vv Mode-11) — ``F`` advertises ``apply_transpose``,
  and ``F.apply_transpose`` ROUTES THROUGH the kernel's tensor-product transpose
  (the rank-1 dual-dyad primitive), not an inline reduction.

Fixture: an asymmetric ≥2-group fissile cell (mixture ``A``: χ per group ≠ νΣf
per group) in a fuel region, with a non-fissile moderator (mixture ``B``) second
region → heterogeneous νΣf field (vv: a flat field nulls redistribution-style
discrimination).

Since CS4c step 4 the fission channel is TWO bindings of one datum: the
scalar rows below gate the ENERGY binding (``IsotropicFission`` — the
solver-held ``fission_op`` the k-outer feeds), and the composite rows
gate the ANGULAR binding (``FissionOperator``, the frame's ℓ=0
conjugation, minted here exactly as the eigen-M posing mints it — its
transpose is the reversed ``full_fission_kernel`` product, so these rows
now pin factor REVERSAL rather than an inline w-spelling; the
independent references are unchanged and the divide-order difference is
ULP-tier, inside every tolerance).

vv Mode-8: every gate uses ``np.testing.*`` / :func:`require` (function calls,
fire under ``python -O``) — NEVER a bare ``assert``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import TensorProductOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField

from tests.transport._integral_kernel_helpers import (
    hand_derived_fission_emission,
    require,
)

pytestmark = pytest.mark.foundation


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


def _solver(groups):
    """Asymmetric heterogeneous fissile fixture (fuel region + moderator region).

    Mixture ``A`` (``groups``) carries an asymmetric emission spectrum χ and
    production cross section νΣf, per-group non-proportional (χ ≁ νΣf), so the
    χ↔νΣf role swap is detectable. Mixture ``B`` (non-fissile) is the second
    region → heterogeneous νΣf across the mesh.
    """
    fuel = get_mixture("A", groups)
    mod = get_mixture("B", groups)
    nx, ny = 6, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2  # fuel
    mat[3:, :] = 0  # moderator (no fission → asymmetric νΣf field)
    mesh = _uniform_2d(nx, ny, 0.2, mat)
    quad = Quadrature.lebedev(order=17)
    return SNSolver(SNMesh(mesh, quad, {2: fuel, 0: mod}))


def _angular_F(solver):
    """The ANGULAR composite binding, minted as the eigen-M posing mints it."""
    from orpheus.transport.operators.fission import FissionOperator

    return FissionOperator.from_solver_data(
        mat_xs=solver.mat_xs, space=solver.sn_mesh.full_field_space,
    )


def _shape(solver):
    return solver.ng, *solver.sn_mesh.spatial_shape


def _asymmetric_field(ng, nx, ny, seed):
    """A field distinct per group AND per cell (no symmetry to hide a swap)."""
    rng = np.random.default_rng(seed)
    return rng.uniform(0.05, 1.0, size=(ng, nx, ny))


# ═══════════════════════════════════════════════════════════════════════
# CORRECTNESS — F† matches a structurally-independent hand loop (role-swapped).
# ═══════════════════════════════════════════════════════════════════════


class TestAdjointFissionCorrectness:
    @pytest.mark.parametrize("groups,seed", [("2g", 20260628), ("4g", 20260629)])
    def test_apply_transpose_matches_hand_derived_adjoint(self, groups, seed):
        r"""``F†ψ* = νΣf·(χ·ψ*)`` equals the EXPLICIT-loop adjoint emission (0 ULP-ish).

        The reference is ``hand_derived_fission_emission(νΣf, χ, ψ*)`` — the SAME
        structurally-independent double-loop the forward gate uses, args
        role-swapped (the loop computes ``arg1_g · Σ_g' arg2_g'·field_g'``, which
        with ``arg1=νΣf, arg2=χ`` IS ``(F†ψ*)_g``). A wrong axis, a dropped
        broadcast, or a χ/νΣf role error disagrees with it.
        """
        solver = _solver(groups)
        op = solver.fission_op  # the ENERGY binding (CS4c step 4)
        ng, nx, ny = _shape(solver)
        psi_star = _asymmetric_field(ng, nx, ny, seed)

        out = op.apply_transpose(psi_star)  # bare-ndarray arm → (ng, nx, ny)
        expected = hand_derived_fission_emission(
            solver.mat_xs.fission_production,  # νΣf as the reconstruction column
            solver.mat_xs.emission_spectrum,  # χ as the contracted row
            psi_star,
        )
        np.testing.assert_allclose(
            out, expected, rtol=1e-13, atol=0.0,
            err_msg="F†.apply_transpose disagrees with the hand-derived "
            "νΣf·(χ·ψ*) adjoint emission — a sign / axis / χ↔νΣf-role bug in F†.",
        )


# ═══════════════════════════════════════════════════════════════════════
# RECIPROCITY — ⟨Fφ, ψ*⟩ == ⟨φ, F†ψ*⟩ (the transpose-defining identity).
# ═══════════════════════════════════════════════════════════════════════


class TestForwardAdjointReciprocity:
    @pytest.mark.parametrize("groups", ["2g", "4g"])
    def test_euclidean_reciprocity(self, groups):
        r"""``⟨Fφ, ψ*⟩ == ⟨φ, F†ψ*⟩`` (full Euclidean sum).

        The defining property of the transpose. ``F`` and ``F†`` both route
        through the bare-ndarray arms here. A F† that is not the genuine
        transpose of F breaks this identity.
        """
        solver = _solver(groups)
        op = solver.fission_op  # the ENERGY binding (CS4c step 4)
        ng, nx, ny = _shape(solver)
        phi = _asymmetric_field(ng, nx, ny, 11)
        psi_star = _asymmetric_field(ng, nx, ny, 12)

        lhs = float((op.apply(phi) * psi_star).sum())            # ⟨Fφ, ψ*⟩
        rhs = float((phi * op.apply_transpose(psi_star)).sum())  # ⟨φ, F†ψ*⟩
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-12,
            err_msg="forward-adjoint reciprocity ⟨Fφ,ψ*⟩ = ⟨φ,F†ψ*⟩ violated — "
            "F† is not the transpose of F.",
        )


# ═══════════════════════════════════════════════════════════════════════
# ROLE-SWAP DISCRIMINATOR — the fixture genuinely constrains χ vs νΣf (Mode-2).
# ═══════════════════════════════════════════════════════════════════════


class TestRoleSwapDiscriminator:
    @pytest.mark.parametrize("groups", ["2g", "4g"])
    def test_role_swap_changes_adjoint_emission(self, groups):
        r"""With χ ≁ νΣf, the χ↔νΣf role swap in F† changes the emission.

        The genuine adjoint-fission trap is the ROLE swap: the CORRECT F† is
        ``νΣf·(χ·ψ*)``; the WRONG one (= the forward F) is ``χ·(νΣf·ψ*)``. With
        asymmetric per-group χ AND νΣf these differ — so the correctness gate
        genuinely constrains the χ-vs-νΣf assignment. If this row shows
        equality, the fixture lost its asymmetry and the gate is blind to the
        swap.
        """
        solver = _solver(groups)
        op = solver.fission_op  # the ENERGY binding (CS4c step 4)
        ng, nx, ny = _shape(solver)
        psi_star = _asymmetric_field(ng, nx, ny, 13)
        chi = solver.mat_xs.emission_spectrum
        nu_sf = solver.mat_xs.fission_production

        correct = hand_derived_fission_emission(nu_sf, chi, psi_star)   # F† = |νΣf⟩⟨χ|
        role_swapped = hand_derived_fission_emission(chi, nu_sf, psi_star)  # the forward F
        require(
            not np.allclose(correct, role_swapped),
            "The fission fixture must be asymmetric enough that a χ↔νΣf ROLE "
            "swap changes the adjoint emission (Mode-2 discriminability). They "
            "agreed — the fixture lost its asymmetry; the F† correctness gate "
            "would be blind to the role swap.",
        )


# ═══════════════════════════════════════════════════════════════════════
# PREDICATE + ROUTING — F is adjointable; routes through the kernel.
# ═══════════════════════════════════════════════════════════════════════


class TestAdjointFissionCapabilityAndRouting:
    def test_fission_advertises_apply_transpose(self):
        op = _solver("4g").fission_op
        require(
            op.is_adjointable,
            "the fission energy binding must advertise the adjoint axis "
            "(F† via the rank-1 dyad swap).",
        )

    def test_apply_transpose_routes_through_kernel_transpose(self, monkeypatch):
        r"""Mode-11: F†.apply_transpose ENTERS the kernel's TP transpose.

        F† routes through :attr:`FissionOperator.kernel`'s
        :meth:`TensorProductOperator.apply_transpose` (the rank-1 dual dyad +
        identity), NOT an inline reduction. Wrap the TP transpose in-process and
        assert the counter advances — a routed-around inline transpose leaves it
        at 0 and reddens this gate (the strictly-stronger Mode-11 proof).
        """
        calls = {"n": 0}
        real = TensorProductOperator.apply_transpose

        def counting(self, x):
            calls["n"] += 1
            return real(self, x)

        monkeypatch.setattr(TensorProductOperator, "apply_transpose", counting)

        solver = _solver("4g")
        op = solver.fission_op
        ng, nx, ny = _shape(solver)
        op.apply_transpose(_asymmetric_field(ng, nx, ny, 14))

        require(
            calls["n"] > 0,
            "F†.apply_transpose did NOT route through the kernel's "
            "TensorProductOperator.apply_transpose — it computed an inline "
            "transpose instead of reusing the rank-1 dual-dyad primitive "
            "(Mode-11: a path that bypasses the rewired reader).",
        )


# ═══════════════════════════════════════════════════════════════════════
# COMPOSITE ARM (#276 A4) — F† on the FullField carrier: the forward's
# ``(1/W)·broadcast ∘ K ∘ (w-Σ)`` pulls back to ``(w·) ∘ Kᵀ ∘ (Σ/W)``.
# ═══════════════════════════════════════════════════════════════════════


def _composite(solver, seed):
    """Random angular FullField (bulk AND trace random, ANGLE-VARYING bulk)."""
    sn = solver.sn_mesh
    rng = np.random.default_rng(seed)
    bulk = AngularFlux(values=rng.uniform(0.05, 1.0, size=(sn.quad.N, sn.ng, *sn.spatial_shape)), space=sn.angular_bulk_space)
    trace = AngularBoundaryFlux(
        values=rng.uniform(0.05, 1.0, size=int(sn.angular_trace.layout.total_size)),
        space=sn.angular_trace,
    )
    return FullField(interior=bulk, boundary=trace)


class TestCompositeTransposeArm:
    @pytest.mark.parametrize("groups", ["2g", "4g"])
    def test_composite_pairing_identity(self, groups):
        r"""``⟨Fψ, χ⟩ == ⟨ψ, F†χ⟩`` on FullField composites (flat coordinates).

        The transpose-defining identity for the COMPOSITE arm — F's forward
        composite has the angular reduce/broadcast structure the bare arm
        lacks, so this row genuinely gates the weight-role swap
        (:meth:`test_weight_swap_discriminator`), which the bare-arm
        reciprocity above cannot see.  Precondition: the cotangent bulk must
        VARY in angle — under an angle-flat χ the weighted and unweighted
        angular reduces collapse to multiples of the SAME iso vector
        (differing only by the global ``N/W`` scale — they are NOT equal;
        for GL-8 the ratio is 4), so the swap loses its angular-SHAPE
        signature and the row's discrimination degenerates to a scale it
        cannot attribute.  The require conservatively keeps the shape
        discriminant live (the Mode-7-style blindness this pins).
        """
        solver = _solver(groups)
        op = _angular_F(solver)
        psi = _composite(solver, 31)
        chi = _composite(solver, 32)
        bulk = np.asarray(chi.interior.values)
        require(
            not np.allclose(bulk, np.broadcast_to(bulk[:1], bulk.shape)),
            "composite cotangent must vary in angle or the weight-swap "
            "error class is invisible to this pairing row.",
        )

        lhs = float(op.apply(psi).to_flat() @ chi.to_flat())
        rhs = float(psi.to_flat() @ op.apply_transpose(chi).to_flat())
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-12,
            err_msg="composite F pairing ⟨Fψ,χ⟩ = ⟨ψ,F†χ⟩ violated — the "
            "composite transpose arm is not the transpose of the composite "
            "forward arm.",
        )

    @pytest.mark.parametrize("groups", ["2g", "4g"])
    def test_composite_matches_independent_spelling(self, groups):
        r"""Composite ``F†χ`` bulk == ``w ⊗ hand-loop(νΣf, χ_spec, Σ_n χ_n/W)``;
        trace exactly zero.

        The reference is STRUCTURALLY INDEPENDENT of the production arm: the
        ordinate sum, the ``/W``, the hand-derived dual-dyad loop, and the
        explicit ``np.multiply.outer`` weight broadcast — no kernel call, no
        factory.  Pins the whole ``(w·) ∘ Kᵀ ∘ (Σ/W)`` composition including
        WHICH side carries the quadrature weights.  The zero trace is the
        pure-bulk pullback (the transpose of emitting nothing into the trace).
        """
        solver = _solver(groups)
        op = _angular_F(solver)
        chi = _composite(solver, 33)
        w = np.asarray(solver.sn_mesh.quad.weights, dtype=float)

        out = op.apply_transpose(chi)
        iso_star = np.asarray(chi.interior.values).sum(axis=0) / float(w.sum())
        expected_bulk = np.multiply.outer(
            w,
            hand_derived_fission_emission(
                solver.mat_xs.fission_production,
                solver.mat_xs.emission_spectrum,
                iso_star,
            ),
        )
        np.testing.assert_allclose(
            np.asarray(out.interior.values), expected_bulk, rtol=1e-12, atol=0.0,
            err_msg="composite F† bulk disagrees with the independent "
            "w ⊗ dual-dyad(Σχ/W) spelling — a reduce/broadcast weight-role "
            "or χ↔νΣf error in the composite arm.",
        )
        np.testing.assert_array_equal(
            np.asarray(out.boundary.values),
            np.zeros_like(np.asarray(out.boundary.values)),
            err_msg="composite F† emitted a non-zero trace — fission is pure "
            "bulk; its transpose must be too.",
        )

    def test_weight_swap_discriminator(self):
        r"""The WEIGHT-role swap (forward's weights kept on the forward sides:
        ``(1/W)·broadcast ∘ Kᵀ ∘ (w-Σ)``) differs O(1) from the true
        transpose AND breaks the pairing — the mutation this gate class
        exists to catch, verified as an explicit wrong-spelling discriminator
        (the arm is inline, so the wrong spelling is computed in-test rather
        than monkeypatched).
        """
        solver = _solver("4g")
        op = _angular_F(solver)
        psi = _composite(solver, 34)
        chi = _composite(solver, 35)
        w = np.asarray(solver.sn_mesh.quad.weights, dtype=float)
        W = float(w.sum())

        true_bulk = np.asarray(op.apply_transpose(chi).interior.values)
        # The WRONG transpose: weights NOT swapped across Kᵀ.
        wrong_iso = np.einsum(
            "n,n...->...", w, np.asarray(chi.interior.values),
        ) / W
        wrong_bulk = (
            np.asarray(op.kernel.apply_transpose(wrong_iso))[None] / W
        ) * np.ones((w.size,) + wrong_iso.shape)
        require(
            not np.allclose(true_bulk, wrong_bulk, rtol=1e-6),
            "weight-role-swapped F† coincided with the true transpose on an "
            "angle-varying cotangent — the discriminator lost its teeth.",
        )

        lhs = float(op.apply(psi).to_flat() @ chi.to_flat())
        wrong_rhs = float(
            (np.asarray(psi.interior.values) * wrong_bulk).sum()
        )
        require(
            abs(lhs - wrong_rhs) > 1e-8 * abs(lhs),
            "the pairing identity FAILED to distinguish the weight-role-"
            "swapped F† — the composite pairing row has no teeth on this "
            "mutation class.",
        )

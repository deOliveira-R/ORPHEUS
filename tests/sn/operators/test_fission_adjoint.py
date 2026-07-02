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
        op = _solver(groups).fission_op
        ng = op.chi.shape[0]
        nx, ny = op.chi.shape[1:]
        psi_star = _asymmetric_field(ng, nx, ny, seed)

        out = op.apply_transpose(psi_star)  # bare-ndarray arm → (ng, nx, ny)
        expected = hand_derived_fission_emission(
            op.mat_xs.fission_production,  # νΣf as the reconstruction column
            op.mat_xs.emission_spectrum,  # χ as the contracted row
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
        op = _solver(groups).fission_op
        ng = op.chi.shape[0]
        nx, ny = op.chi.shape[1:]
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
        op = _solver(groups).fission_op
        ng = op.chi.shape[0]
        nx, ny = op.chi.shape[1:]
        psi_star = _asymmetric_field(ng, nx, ny, 13)
        chi = op.mat_xs.emission_spectrum
        nu_sf = op.mat_xs.fission_production

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
            "FissionOperator must advertise the adjoint axis (F† via the "
            "rank-1 dyad swap).",
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

        op = _solver("4g").fission_op
        ng = op.chi.shape[0]
        nx, ny = op.chi.shape[1:]
        op.apply_transpose(_asymmetric_field(ng, nx, ny, 14))

        require(
            calls["n"] > 0,
            "F†.apply_transpose did NOT route through the kernel's "
            "TensorProductOperator.apply_transpose — it computed an inline "
            "transpose instead of reusing the rank-1 dual-dyad primitive "
            "(Mode-11: a path that bypasses the rewired reader).",
        )

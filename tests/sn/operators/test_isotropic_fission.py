r"""Intrinsic gates for :class:`IsotropicFission` — the fission ENERGY binding (CS4c step 4a).

The rank-1 dyad :math:`F = |\chi\rangle\langle\nu\Sigma_f|` on the scalar
flux — what the SN k-outer, the homogeneous :math:`K = A^{-1}F`, and the
diffusion scalar composite consume (the ANGULAR binding is
:class:`~orpheus.transport.operators.fission.FissionOperator`, gated in
``test_fission_operator.py``). Pinned in isolation on a synthetic
two-fissile-material carrier:

* **apply ≡ hand-rolled per-cell dyad** — plain ``np.dot`` loop, no shared
  code (`vv` L11); asymmetric χ AND νΣf per material so a factor swap or a
  material mixup is detectable.
* **transpose** — hand-rolled DUAL dyad + Euclidean reciprocity
  ``⟨Fφ,ψ*⟩=⟨φ,Fᵀψ*⟩`` + a deliberately-swapped flip CONTROL that must NOT
  match (the catcher's catcher — `vv` #11's negative leg for the
  factor-swap theorem).
* **admission** — per-END ng conformity; the angular-leading-bulk
  construction refusal (that consumer wants the angular binding); the
  composite-transpose refusal.
* **structure** — ``kernel`` is a cached §5.6 ``TensorProductOperator``
  (the satellite ruling: minted once per binding); ``production_rate`` is
  the typed :class:`ReactionRateFunctional`; adjointable, structurally
  non-invertible.

The scalar-COMPOSITE arm rides the production consumer's gates
(``tests/diffusion/test_operators.py`` — the arm's output blocks ride the
operand's spaces, so the real diffusion composite is the honest fixture).
vv Mode-8: ``np.testing.*`` / :func:`require` only (fire under ``-O``).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import Mesh2D
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.mesh.material_xs_field import MaterialXSField
from orpheus.transport.operators.isotropic_transfer import IsotropicFission

pytestmark = pytest.mark.foundation


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


_NG, _NX = 2, 6
# Asymmetric per-material factor pairs (χ₀ ≠ χ₁, νΣf₀ ≠ νΣf₁, and χ ∦ νΣf
# in either material — so swapping the dyad's factors moves every output).
_CHI = {0: np.array([0.9, 0.1]), 1: np.array([0.7, 0.3])}
_NU_SIG_F = {0: np.array([0.013, 0.26]), 1: np.array([0.002, 0.11])}


def _mat_xs():
    z = np.zeros(_NG)
    materials = {
        mid: Mixture(
            SigC=z.copy(), SigL=z.copy(),
            SigF=_NU_SIG_F[mid] / 2.4, SigP=_NU_SIG_F[mid].copy(),
            SigT=np.ones(_NG),
            SigS=[csr_matrix(np.zeros((_NG, _NG)))],
            Sig2=[csr_matrix(np.zeros((_NG, _NG)))],
            chi=_CHI[mid].copy(),
        )
        for mid in _CHI
    }
    mat_map = np.zeros((_NX, 1), dtype=int)
    mat_map[_NX // 2:, :] = 1
    mesh = Mesh2D(
        edges_x=np.arange(_NX + 1, dtype=float),
        edges_y=np.arange(2, dtype=float),
        mat_map=mat_map,
    )
    return MaterialXSField.from_mesh(MaterialMesh(mesh, materials))


def _op(mat_xs=None):
    mat_xs = mat_xs if mat_xs is not None else _mat_xs()
    return IsotropicFission.from_material_xs(
        mat_xs, space=mat_xs.mesh.bulk_space,
    )


def _mid_of(ix):
    return 0 if ix < _NX // 2 else 1


class TestForward:
    def test_apply_matches_hand_rolled_dyad(self):
        op = _op()
        phi = np.random.default_rng(7).uniform(0.1, 1.0, size=(_NG, _NX, 1))
        out = np.asarray(op.apply(phi))
        ref = np.empty_like(phi)
        for ix in range(_NX):
            mid = _mid_of(ix)
            rate = float(np.dot(_NU_SIG_F[mid], phi[:, ix, 0]))
            ref[:, ix, 0] = _CHI[mid] * rate
        np.testing.assert_allclose(out, ref, rtol=1e-15, atol=0)

    def test_carrier_values_are_unwrapped(self):
        op = _op()
        phi = np.random.default_rng(3).uniform(0.1, 1.0, size=(_NG, _NX, 1))

        class _Carrier:
            values = phi

        np.testing.assert_array_equal(op.apply(_Carrier()), op.apply(phi))


class TestTranspose:
    def test_transpose_matches_hand_rolled_dual_dyad(self):
        op = _op()
        psi_star = np.random.default_rng(11).uniform(
            -1.0, 1.0, size=(_NG, _NX, 1),
        )
        out = op.apply_transpose(psi_star)
        ref = np.empty_like(psi_star)
        for ix in range(_NX):
            mid = _mid_of(ix)
            importance = float(np.dot(_CHI[mid], psi_star[:, ix, 0]))
            ref[:, ix, 0] = _NU_SIG_F[mid] * importance
        np.testing.assert_allclose(out, ref, rtol=1e-15, atol=0)

    def test_flip_control_swapped_factors_do_not_match(self):
        """The catcher's catcher: a hand reference with χ↔νΣf UN-swapped
        (i.e. the forward dyad used as the transpose) must disagree —
        proves the transpose row above can red on a factor-order bug."""
        op = _op()
        psi_star = np.random.default_rng(11).uniform(
            -1.0, 1.0, size=(_NG, _NX, 1),
        )
        out = op.apply_transpose(psi_star)
        wrong = np.empty_like(psi_star)
        for ix in range(_NX):
            mid = _mid_of(ix)
            rate = float(np.dot(_NU_SIG_F[mid], psi_star[:, ix, 0]))
            wrong[:, ix, 0] = _CHI[mid] * rate
        require(
            bool(np.max(np.abs(out - wrong)) > 1e-3),
            "the flip control matched — the transpose gate cannot red",
        )

    def test_euclidean_reciprocity(self):
        op = _op()
        rng = np.random.default_rng(5)
        phi = rng.uniform(0.1, 1.0, size=(_NG, _NX, 1))
        psi_star = rng.uniform(-1.0, 1.0, size=(_NG, _NX, 1))
        lhs = float(np.sum(np.asarray(op.apply(phi)) * psi_star))
        rhs = float(np.sum(phi * op.apply_transpose(psi_star)))
        np.testing.assert_allclose(lhs, rhs, rtol=1e-13)

    def test_composite_transpose_refuses(self):
        from orpheus.transport.full_field import FullField

        op = _op()
        with pytest.raises(TypeError, match="angular binding"):
            op.apply_transpose(
                FullField.__new__(FullField)  # type: ignore[call-arg]
            )


class TestAdmission:
    def test_wrong_ng_end_refuses(self):
        from orpheus.numerics.axis import Axis, BasisKind
        from orpheus.numerics.space import FunctionSpace

        mat_xs = _mat_xs()
        wrong = FunctionSpace.of_axes(
            Axis("energy", (_NG + 1,), kind=BasisKind.NODAL),
            Axis("x", (_NX,), kind=BasisKind.NODAL),
            Axis("y", (1,), kind=BasisKind.NODAL),
        )
        with pytest.raises((TypeError, ValueError)):
            IsotropicFission.from_material_xs(mat_xs, space=wrong)

    def test_angular_leading_bulk_refuses_to_construct(self):
        from orpheus.numerics.axis import Axis, BasisKind
        from orpheus.numerics.space import FunctionSpace

        mat_xs = _mat_xs()
        angularish = FunctionSpace.of_axes(
            Axis("angle", (8,), kind=BasisKind.NODAL),
            Axis("energy", (_NG,), kind=BasisKind.NODAL),
            Axis("x", (_NX,), kind=BasisKind.NODAL),
            Axis("y", (1,), kind=BasisKind.NODAL),
        )
        with pytest.raises(TypeError, match="ANGULAR binding"):
            IsotropicFission.from_material_xs(mat_xs, space=angularish)

    def test_classmethod_equals_ctor(self):
        """G-C1's row for this class: tier-2 ≡ hand-extracted ctor on the
        same inputs (`vv` #28 — the fixture ctor provably represents the
        production factory)."""
        from orpheus.transport.material_field import FissionMaterialField

        mat_xs = _mat_xs()
        via_classmethod = IsotropicFission.from_material_xs(
            mat_xs, space=mat_xs.mesh.bulk_space,
        )
        via_ctor = IsotropicFission(
            FissionMaterialField.from_material_xs(mat_xs),
            domain=mat_xs.mesh.bulk_space,
            codomain=mat_xs.mesh.bulk_space,
        )
        phi = np.random.default_rng(2).uniform(
            0.1, 1.0, size=(_NG, _NX, 1),
        )
        np.testing.assert_array_equal(
            via_classmethod.apply(phi), via_ctor.apply(phi),
        )


class TestStructure:
    def test_kernel_is_a_cached_tensor_product(self):
        from orpheus.numerics.operator import TensorProductOperator

        op = _op()
        require(
            isinstance(op.kernel, TensorProductOperator),
            "kernel must expose the §5.6 TensorProductOperator surface",
        )
        require(op.kernel is op.kernel, "kernel must be minted ONCE")

    def test_production_rate_is_the_typed_functional(self):
        from orpheus.transport.reaction_rate_functional import (
            ReactionRateFunctional,
        )

        op = _op()
        require(
            isinstance(op.production_rate, ReactionRateFunctional),
            "production_rate must be the typed diagnostic functional",
        )
        require(
            op.production_rate is op.production_rate,
            "production_rate must be minted ONCE",
        )

    def test_predicates(self):
        op = _op()
        require(op.is_adjointable, "the dual dyad is spelled")
        require(
            not op.is_invertible,
            "a rank-1 production operator is singular",
        )

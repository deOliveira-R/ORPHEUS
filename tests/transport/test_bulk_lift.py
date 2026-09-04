r"""The lift's own laws — how a BULK action enters the composite (CS4c step 5).

:mod:`orpheus.transport.operators.lift` is the ONE home of *extension by
zero on the trace* (`[M]` 2026-09-04: nine hand spellings before it) and
of the family's two carrier parses. These rows pin the verbs and
:class:`~orpheus.transport.operators.lift.BulkLift` on their own terms
(the intrinsic-properties standard), on the diffusion scalar composite:

* **definition** — ``lift(inner).apply(ψ).interior.values`` IS
  ``inner.apply(ψ.interior.values)`` (``array_equal``: no arithmetic in
  the lift), the trace is the zero SOURCE/SINK of the operand's boundary
  class on the operand's trace space; the transpose likewise;
* **linearity** with an ACTIVATION leg (`lessons L40c` — a zero morphism
  satisfies every linearity row);
* **transpose reciprocity** ``⟨lift ψ, χ⟩ = ⟨ψ, liftᵀ χ⟩`` on raw flats;
* **capability axes** delegate to the inner; ``block_role`` is BULK;
* **admission** — the inner bound off the composite's interiors, a plain
  end, a bare array, a typed field, a composite on ANOTHER interior: each
  a typed refusal naming the operator (the ends select the carrier);
* **the assembly embed** is index-identity on the bulk block and zero on
  the trace rows/columns, and ``lift.assemble()`` acts as ``lift.apply``
  on the flat.

The inner used for the value rows is the diffusion energy binding
``IsotropicScattering + IsotropicN2N`` on the mesh's plain ``bulk_space``
— the binding the diffusion solver lifts after step 5 (R-4).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest
from scipy import sparse

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.diffusion import DiffusionMesh
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.operator import BlockRole, MissingAssembly
from orpheus.transport.fields._bases import FieldRole
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.bound_operator import BoundOperator
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.lift import (
    BulkLift,
    admit_array,
    admit_composite,
    embed_bulk_assembly,
    lift_bulk_action,
)
from orpheus.transport.source_sinks import (
    ScalarBoundarySourceSink,
    ScalarSourceSink,
)

pytestmark = pytest.mark.foundation

# An asymmetric 2-group scatter matrix (anti-#3 / anti-#4: the transpose
# convention is load-bearing on every row below).
_SIG_S = np.array([[0.38, 0.10], [0.05, 0.90]])
_SIG_T = np.array([0.60, 1.20])


def _mesh(n_cells: int = 4) -> DiffusionMesh:
    mix = make_mixture(
        sig_t=_SIG_T, sig_c=_SIG_T - _SIG_S.sum(axis=1),
        sig_f=np.zeros(2), nu=np.zeros(2), chi=np.zeros(2), sig_s=_SIG_S,
    )
    mesh1d = Mesh1D(
        edges=np.linspace(0.0, 2.0, n_cells + 1), mat_ids=np.zeros(n_cells, dtype=int),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return DiffusionMesh(mesh1d, {0: mix})


@pytest.fixture
def mesh() -> DiffusionMesh:
    return _mesh()


@pytest.fixture
def inner(mesh):
    mat_xs = mesh.material_xs_field()
    return (
        IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space)
        + IsotropicN2N.from_material_xs(mat_xs, space=mesh.bulk_space)
    )


@pytest.fixture
def lift(mesh, inner) -> BulkLift:
    space = mesh.full_field_space
    return BulkLift(inner, domain=space, codomain=space)


def _flux(mesh: DiffusionMesh, seed: int) -> FullField:
    rng = np.random.default_rng(seed)
    return FullField(
        interior=ScalarFlux(values=rng.random(mesh.bulk_space.shape) + 0.5, space=mesh.bulk_space),
        boundary=ScalarBoundaryFlux(values=rng.random(mesh.scalar_trace.shape) + 0.1, space=mesh.scalar_trace),
    )


class TestDefinition:
    def test_interior_is_the_inner_action_verbatim(self, mesh, inner, lift):
        psi = _flux(mesh, 1)
        out = lift.apply(psi)
        assert type(out.interior) is ScalarSourceSink
        np.testing.assert_array_equal(out.interior.values, inner.apply(psi.interior.values))
        assert out.interior.space is mesh.bulk_space

    def test_trace_is_the_zero_source_sink_on_the_operand_trace_space(self, mesh, lift):
        psi = _flux(mesh, 1)
        out = lift.apply(psi)
        assert type(out.boundary) is ScalarBoundarySourceSink
        assert out.boundary.space is psi.boundary.space
        assert not out.boundary.values.any()

    def test_transpose_is_the_inner_transpose_verbatim(self, mesh, inner, lift):
        chi = _flux(mesh, 2)
        out = lift.apply_transpose(chi)
        np.testing.assert_array_equal(
            out.interior.values, inner.apply_transpose(chi.interior.values),
        )
        assert type(out.boundary) is ScalarBoundarySourceSink
        assert not out.boundary.values.any()


class TestLaws:
    def test_activation(self, mesh, lift):
        assert np.abs(lift.apply(_flux(mesh, 3)).interior.values).max() > 0.0

    def test_linearity(self, mesh, lift):
        a, b = _flux(mesh, 4), _flux(mesh, 5)
        lhs = lift.apply(2.0 * a - 3.0 * b)
        rhs = 2.0 * lift.apply(a) - 3.0 * lift.apply(b)
        np.testing.assert_allclose(lhs.to_flat(), rhs.to_flat(), rtol=1e-13, atol=0.0)

    def test_transpose_reciprocity_on_raw_flats(self, mesh, lift):
        psi, chi = _flux(mesh, 6), _flux(mesh, 7)
        lhs = float(lift.apply(psi).to_flat() @ chi.to_flat())
        rhs = float(psi.to_flat() @ lift.apply_transpose(chi).to_flat())
        assert lhs != 0.0
        assert abs(lhs - rhs) <= 1e-13 * abs(lhs)

    def test_capability_axes_and_role(self, inner, lift):
        assert lift.block_role is BlockRole.BULK
        assert lift.is_adjointable is inner.is_adjointable is True
        assert lift.is_assemblable is inner.is_assemblable
        assert lift.is_invertible is False


class TestAdmission:
    def test_inner_bound_off_the_interiors_is_refused(self, mesh):
        other = _mesh(5)
        wrong = IsotropicScattering.from_material_xs(
            other.material_xs_field(), space=other.bulk_space,
        )
        space = mesh.full_field_space
        with pytest.raises(TypeError, match="bind the inner operator on the composite's interior"):
            BulkLift(wrong, domain=space, codomain=space)

    def test_a_plain_end_is_refused(self, mesh, inner):
        with pytest.raises(TypeError, match="composite FullFieldSpace"):
            BulkLift(inner, domain=mesh.bulk_space, codomain=mesh.bulk_space)

    def test_a_bare_array_is_refused_naming_the_operator(self, mesh, lift):
        with pytest.raises(TypeError, match="BulkLift.*FullField"):
            lift.apply(np.zeros(mesh.bulk_space.shape))

    def test_a_composite_on_another_interior_is_refused(self, mesh, lift):
        other = _flux(_mesh(5), 8)
        with pytest.raises(TypeError, match="body its ends select"):
            lift.apply(other)
        with pytest.raises(TypeError, match="body its ends select"):
            lift.apply_transpose(other)

    def test_admit_array_refusals_name_the_lift_and_the_values(self, mesh, inner):
        plain = IsotropicScattering.from_material_xs(
            mesh.material_xs_field(), space=mesh.bulk_space,
        )
        psi = _flux(mesh, 9)
        with pytest.raises(TypeError, match="BulkLift"):
            admit_array(plain, psi)
        with pytest.raises(TypeError, match=r"\.values"):
            admit_array(plain, psi.interior)
        with pytest.raises(TypeError, match="shape"):
            admit_array(plain, np.zeros((2, 5)))
        with pytest.raises(TypeError, match="IsotropicScattering"):
            admit_array(plain, object())
        assert admit_array(plain, psi.interior.values) is psi.interior.values

    def test_admit_composite_reads_the_named_end(self, mesh, lift):
        psi = _flux(mesh, 10)
        assert admit_composite(lift, psi, end="domain") is psi
        assert admit_composite(lift, psi, end="codomain") is psi
        with pytest.raises(TypeError, match="composite FullFieldSpace"):
            admit_composite(lift.inner, psi)  # a plain-bound operator has no interior


class TestLiftBulkAction:
    @pytest.mark.parametrize(
        "role, trace_cls", [
            (FieldRole.SOURCE_SINK, ScalarBoundarySourceSink),
            (FieldRole.FLUX, ScalarBoundaryFlux),
        ], ids=["source_sink", "flux"],
    )
    def test_the_zero_trace_takes_the_requested_role(self, mesh, role, trace_cls):
        psi = _flux(mesh, 11)
        out = lift_bulk_action(
            psi, lambda bulk: bulk.into_role(role, 2.0 * bulk.values), trace_role=role,
        )
        assert type(out.boundary) is trace_cls
        assert out.boundary.space is psi.boundary.space
        assert not out.boundary.values.any()
        np.testing.assert_array_equal(out.interior.values, 2.0 * psi.interior.values)
        assert type(out.interior) is psi.interior.role_partner(role)


@dataclass(eq=False)
class _ShapedDiagonal(BoundOperator):
    r"""A shape-preserving, assemblable bulk operator — the contract a
    lifted inner honours (bare array of the interior shape in and out)."""

    diagonal: np.ndarray

    def apply(self, x, /):
        return self.diagonal * x

    def apply_transpose(self, x, /):
        return self.diagonal * x

    @property
    def is_adjointable(self) -> bool:
        return True

    @property
    def is_assemblable(self) -> bool:
        return True

    def assemble(self) -> SparseAssembledOperator:
        return SparseAssembledOperator(
            sparse.diags_array(self.diagonal.ravel()), domain=self.domain, codomain=self.codomain,
        )


class TestAssemblyEmbed:
    def test_embed_is_index_identity_on_the_bulk_block(self, mesh):
        rng = np.random.default_rng(12)
        n = int(np.prod(mesh.bulk_space.shape))
        block = sparse.random_array((n, n), density=0.5, rng=rng)
        bulk = SparseAssembledOperator(block, domain=mesh.bulk_space, codomain=mesh.bulk_space)
        space = mesh.full_field_space
        embedded = embed_bulk_assembly(bulk, domain=space, codomain=space).as_matrix()
        n_total = int(space.shape[0])
        assert embedded.shape == (n_total, n_total) and n_total > n
        np.testing.assert_array_equal(embedded[:n, :n], bulk.as_matrix())
        assert not embedded[n:, :].any() and not embedded[:, n:].any()

    def test_a_bulk_of_the_wrong_size_is_refused(self, mesh):
        space = mesh.full_field_space
        bulk = SparseAssembledOperator(sparse.eye_array(3), domain=mesh.bulk_space, codomain=mesh.bulk_space)
        with pytest.raises(ValueError, match="not bound on the composite's interiors"):
            embed_bulk_assembly(bulk, domain=space, codomain=space)

    def test_lift_assemble_acts_as_lift_apply_on_the_flat(self, mesh):
        rng = np.random.default_rng(13)
        inner = _ShapedDiagonal(
            rng.random(mesh.bulk_space.shape) + 0.5,
            domain=mesh.bulk_space, codomain=mesh.bulk_space,
        )
        space = mesh.full_field_space
        lift = BulkLift(inner, domain=space, codomain=space)
        psi = _flux(mesh, 14)
        assert lift.is_assemblable
        np.testing.assert_array_equal(
            lift.assemble().apply(psi.to_flat()), lift.apply(psi).to_flat(),
        )

    def test_a_non_assemblable_inner_refuses_naming_the_lift(self, lift):
        assert not lift.is_assemblable  # the plain energy binding has no assembly yet (step 5 lands it)
        with pytest.raises(MissingAssembly, match="BulkLift.assemble"):
            lift.assemble()

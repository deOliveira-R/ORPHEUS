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
  the trace rows/columns, and ``assemble(lift)`` acts as ``lift.apply``
  on the flat.

The inner used for the value rows is the diffusion energy binding
``IsotropicScattering + IsotropicN2N`` on the mesh's plain ``bulk_space``
— the binding the diffusion solver lifts after step 5 (R-4).

Layer DAG Phase 4 (2026-11-01): the file builds material meshes and spells each space it reads from the raw
constructors of orpheus.transport.phase_spaces (the SN problem's and the diffusion mesh's cached mints are those
constructors on the problem's own ingredients, gated in
tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py and
tests/gates/diffusion/test_scalar_trace_is_the_raw_constructor.py).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest
from scipy import sparse

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import scalar_trace_space
from tests.gates.transport._phase_space_cases import scalar_full_field_space
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.assembled_operator import MissingAssembly, assemble
from orpheus.numerics.operator import BlockRole
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
from orpheus.transport.fields import (
    ScalarBoundarySourceSink,
    ScalarSourceSink,
)
import orpheus.numerics.assembled_operator as asm

pytestmark = pytest.mark.foundation

# An asymmetric 2-group scatter matrix (anti-#3 / anti-#4: the transpose
# convention is load-bearing on every row below).
_SIG_S = np.array([[0.38, 0.10], [0.05, 0.90]])
_SIG_T = np.array([0.60, 1.20])


def _mesh(n_cells: int = 4) -> MaterialMesh:
    mix = make_mixture(
        sig_t=_SIG_T, sig_c=_SIG_T - _SIG_S.sum(axis=1),
        sig_f=np.zeros(2), nu=np.zeros(2), chi=np.zeros(2), sig_s=_SIG_S,
    )
    geometry = StructuredGeometry.slab(
        (0.0, 2.0), (0,), left=MirrorDeck(), right=0,
    )
    mesh1d = Mesher(geometry).partition(CellsByCount.uniform_width(n_cells)).mesh
    return MaterialMesh(mesh1d, {0: mix})


@pytest.fixture
def mesh() -> MaterialMesh:
    return _mesh()


@pytest.fixture
def inner(mesh):
    mat_xs = mesh.mat_xs
    return (
        IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space)
        + IsotropicN2N.from_material_xs(mat_xs, space=mesh.bulk_space)
    )


@pytest.fixture
def lift(mesh, inner) -> BulkLift:
    space = scalar_full_field_space(mesh)
    return BulkLift(inner, domain=space, codomain=space)


def _flux(mesh: MaterialMesh, seed: int) -> FullField:
    rng = np.random.default_rng(seed)
    return FullField(
        interior=ScalarFlux(values=rng.random(mesh.bulk_space.shape) + 0.5, space=mesh.bulk_space),
        boundary=ScalarBoundaryFlux(values=rng.random(scalar_trace_space(mesh).shape) + 0.1, space=scalar_trace_space(mesh)),
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
            other.mat_xs, space=other.bulk_space,
        )
        space = scalar_full_field_space(mesh)
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

    def test_admit_array_refuses_the_composite_and_parses_the_rest(self, mesh, inner):
        """``admit_array`` keeps two refusals of its own (the composite names
        :class:`BulkLift`; a non-array, non-member by type); everything else is
        the end's one parse (D1, 2026-10-30): a typed member of the end is
        admitted with its values, a wrong-shape array refused at the exact
        shape (``ValueError``)."""
        plain = IsotropicScattering.from_material_xs(
            mesh.mat_xs, space=mesh.bulk_space,
        )
        psi = _flux(mesh, 9)
        with pytest.raises(TypeError, match="BulkLift"):
            admit_array(plain, psi)
        assert admit_array(plain, psi.interior) is psi.interior.values
        with pytest.raises(ValueError, match="exact shape"):
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

    def assembled_matrix(self, emit, /):
        # A leaf: its raw emission; the functor wraps it (ruling 62).
        return sparse.diags_array(self.diagonal.ravel())


@dataclass(eq=False)
class _NoAssembly(BoundOperator):
    """A bulk operator with no assembly mode (the identity) — the lift's
    capability delegation has something to refuse."""

    def apply(self, x, /):
        return x

    def apply_transpose(self, x, /):
        return x

    @property
    def is_adjointable(self) -> bool:
        return True


class TestAssemblyEmbed:
    def test_embed_is_index_identity_on_the_bulk_block(self, mesh):
        rng = np.random.default_rng(12)
        n = int(np.prod(mesh.bulk_space.shape))
        block = sparse.random_array((n, n), density=0.5, rng=rng)
        space = scalar_full_field_space(mesh)
        embedded = embed_bulk_assembly(block, domain=space, codomain=space).toarray()
        n_total = int(space.shape[0])
        assert embedded.shape == (n_total, n_total) and n_total > n
        np.testing.assert_array_equal(embedded[:n, :n], block.toarray())
        assert not embedded[n:, :].any() and not embedded[:, n:].any()

    def test_a_bulk_of_the_wrong_size_is_refused(self, mesh):
        space = scalar_full_field_space(mesh)
        bulk = sparse.coo_array(sparse.eye_array(3))
        # The refusal is the one matrix admission (#683), against the interiors.
        with pytest.raises(ValueError, match=r"embed_bulk_assembly: a 3×3 matrix does not map .* which needs 8×8"):
            embed_bulk_assembly(bulk, domain=space, codomain=space)

    def test_lift_assemble_acts_as_lift_apply_on_the_flat(self, mesh):
        rng = np.random.default_rng(13)
        inner = _ShapedDiagonal(
            rng.random(mesh.bulk_space.shape) + 0.5,
            domain=mesh.bulk_space, codomain=mesh.bulk_space,
        )
        space = scalar_full_field_space(mesh)
        lift = BulkLift(inner, domain=space, codomain=space)
        psi = _flux(mesh, 14)
        assert lift.is_assemblable
        np.testing.assert_array_equal(
            assemble(lift).apply(psi.to_flat()), lift.apply(psi).to_flat(),
        )

    def test_the_diffusion_binding_assembles_through_the_lift(self, mesh, inner, lift):
        """The plain energy pair assembles on its own ends (CS4c step 5) and
        the lift embeds it: the composite emission's bulk block IS the
        plain emission (``array_equal``), the trace rows/columns are zero."""
        assert inner.is_assemblable and lift.is_assemblable
        n = int(np.prod(mesh.bulk_space.shape))
        embedded = assemble(lift).as_matrix()
        np.testing.assert_array_equal(embedded[:n, :n], assemble(inner).as_matrix())
        assert not embedded[n:, :].any() and not embedded[:, n:].any()

    def test_a_non_assemblable_inner_refuses_naming_the_lift(self, mesh):
        inner = _NoAssembly(domain=mesh.bulk_space, codomain=mesh.bulk_space)
        space = scalar_full_field_space(mesh)
        lift = BulkLift(inner, domain=space, codomain=space)
        assert not lift.is_assemblable
        # The functor runs the lift's own recursion, so the refusal names
        # the inner part under the lift's law (S3 review, qa F1).
        with pytest.raises(MissingAssembly, match="BulkLift's assembly: _NoAssembly emits no assembly"):
            assemble(lift)


# ── moved here by the layer DAG Phase 2 refresh (specs_2r §5) ──


def _recording_emit():
    calls: list[object] = []

    def emit(part):
        calls.append(part)
        return asm._emit(part)
    return emit, calls


@pytest.mark.rests_on("tests/gates/numerics/test_assembly_functor.py::TestDomainCompositesRecurseThroughTheArgument::test_g2_coupled_block_law")
def test_g2_bulk_lift_embedding(mesh) -> None:
    """The lift's embedding recurses through the ``emit`` it is handed: one call, on the inner (numerics_untangling
    S3, gate G2). Moved from numerics' ``test_assembly_functor.py`` (BulkLift is transport's); the fixture is this
    file's mesh instead of the diffusion test module's config."""
    inner = _ShapedDiagonal(
        np.random.default_rng(13).random(mesh.bulk_space.shape) + 0.5,
        domain=mesh.bulk_space, codomain=mesh.bulk_space,
    )
    space = scalar_full_field_space(mesh)
    emit, calls = _recording_emit()
    BulkLift(inner, domain=space, codomain=space).assembled_matrix(emit)
    assert len(calls) == 1 and calls[0] is inner

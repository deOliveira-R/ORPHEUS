r"""The diffusion operator family — object-level gates (#290 P4).

The loss :math:`A = L + C - S - B` on the scalar composite, gated at the
OBJECT level (vv Mode 12: the k-references downstream are
invariant-functional — a resolvent factor swap or a full transpose is
spectrally invisible — so the committed catcher is the matrix itself):

* **the stencil gate** — ``[A]`` via
  :class:`~orpheus.numerics.flat_operator.FlattenedOperator` ≡ an
  independently hand-posed FD matrix (plain loops, the crosswalk
  conventions hand-transcribed — vv L11: production derives through
  operators, the test hand-lists), on a HETEROGENEOUS 2-material,
  non-uniform 4-cell, 2-group slab with ASYMMETRIC scatter (anti-#3 /
  anti-#4), across three BC configs;
* **mutation teeth** — the three named error classes RED the stencil
  gate: D-face pairing swap, :math:`\Sigma_t`-vs-:math:`\Sigma_a`
  confusion, scatter transpose (+ the boundary-closure sign flip);
* **the family laws** — column-sum conservation
  :math:`\mathbf{1}^T(C - S) = \Sigma_a`, the M-matrix sign pattern,
  per-group SPD of the bulk Schur complement, the reflective
  constant-mode annihilation :math:`(L - B_{\rm refl})\,[\phi_0,
  \phi_0/4] = 0` EXACTLY, and the dense resolvent
  (``MatrixInverseOperator`` over the flattened loss);
* **the scalar-composite substrate + shared-operator arms** minted /
  widened in P4 (``ScalarBoundarySourceSink``, the composite carrier
  ``DiffusionMesh.full_field_space``, the per-cell ``D`` gather, and
  the C/S/F scalar arms — each pinned against its own bare-kernel arm,
  the single-source cross-arm discipline).

Boundary laws enter through the MESH (#290 P7a): every config's laws
are declared as ``BC`` tags on the ``Mesh1D`` and realized at
``DiffusionMesh`` construction (``mesh.bc``); the boundary operator
``B`` reads them off the mesh. Mesh-construction laws themselves are
gated in ``test_augmented_mesh.py``.

Conventions: ``.claude/plans/diffusion_crosswalk.md``. Fixture flat
layout (ng=2, nx=4): bulk ``g*nx + i`` → indices 0–7; trace xmin slot
``(2, 2)`` C-raveled at 8 (J⁺g0, J⁺g1, J⁻g0, J⁻g1), xmax slot at 12.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.diffusion import (
    DiffusionBoundaryOperator,
    DiffusionMesh,
    LeakageOperator,
)
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh1D
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import (
    BlockRole,
    BoundaryOperator,
    FullOperator,
    NotInvertible,
)
from orpheus.transport.fields.scalar_boundary_flux import ScalarBoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.mesh.material_xs_field import MaterialXSField
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicFission,
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.lift import BulkLift
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.source_sinks import (
    ScalarBoundarySourceSink,
    ScalarSourceSink,
)

# Per-class V&V marks (the stencil gate is l0 term-verification; the
# rest foundation) — a file-level mark would conflict with the class
# mark and warn on every stencil test.


# ═══════════════════════════════════════════════════════════════════════
# Fixture: heterogeneous 2-material, non-uniform 4-cell, 2G slab with
# asymmetric scatter — every anti-degeneracy lever pulled (anti-#3 2G;
# anti-#4 heterogeneous; asymmetric Σ_s so a transpose is observable;
# non-uniform h so a face-pairing swap is observable).
# ═══════════════════════════════════════════════════════════════════════

_SIG_T_A = np.array([0.2181, 0.7850])
_SIG_T_B = np.array([0.3416, 0.9431])
_SIG_S_A = np.array([[0.1900, 0.0160], [0.0, 0.4200]])   # [g_from, g_to]
_SIG_S_B = np.array([[0.1000, 0.0020], [0.0, 0.0500]])
_EDGES = np.array([0.0, 0.5, 1.5, 3.0, 5.0])             # non-uniform
_MAT_IDS = np.array([0, 1, 1, 0])


_SIG_F_A = np.array([0.0024, 0.0489])


def _fixture_materials() -> dict[int, object]:
    # σ_c derived so the balance σ_t = σ_c + σ_f + Σ_g' σ_s(g→g') holds
    # EXACTLY — the column-sum conservation gate reads σ_a off the
    # balanced mixture.
    mix_a = make_mixture(
        sig_t=_SIG_T_A,
        sig_c=_SIG_T_A - _SIG_F_A - _SIG_S_A.sum(axis=1),
        sig_f=_SIG_F_A, nu=np.array([2.54, 2.47]),
        chi=np.array([1.0, 0.0]), sig_s=_SIG_S_A,
    )
    mix_b = make_mixture(
        sig_t=_SIG_T_B,
        sig_c=_SIG_T_B - _SIG_S_B.sum(axis=1),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=_SIG_S_B,
    )
    return {0: mix_a, 1: mix_b}


def _diffusion_mesh(
    bc_left: BC = BC("reflective"), bc_right: BC = BC("vacuum"),
) -> DiffusionMesh:
    """The fixture phase space; boundary laws declared as mesh BC tags
    and realized at construction (#290 P7a)."""
    mesh1d = Mesh1D(
        edges=_EDGES, mat_ids=_MAT_IDS,
        bc_left=bc_left, bc_right=bc_right,
    )
    return DiffusionMesh(mesh1d, _fixture_materials())


@pytest.fixture
def mesh() -> DiffusionMesh:
    return _diffusion_mesh()


@pytest.fixture
def mat_xs(mesh) -> MaterialXSField:
    return mesh.material_xs_field()


@pytest.fixture
def template(mesh) -> FullField:
    return FullField.zeros(
        interior=ScalarFlux, boundary=ScalarBoundaryFlux, space=mesh.full_field_space,
    )


@pytest.fixture
def flux(mesh) -> FullField:
    rng = np.random.default_rng(42)
    return FullField(
        interior=ScalarFlux(values=rng.random((2, 4)) + 0.5, space=mesh.bulk_space),
        boundary=ScalarBoundaryFlux(values=rng.random(mesh.scalar_trace.shape[0]) + 0.1, space=mesh.scalar_trace),
    )


def _random_flux(mesh: DiffusionMesh) -> FullField:
    """The ``flux`` fixture's composite on an arbitrary mesh instance
    (fields and operators must share the ONE mesh — identity guard)."""
    rng = np.random.default_rng(42)
    return FullField(
        interior=ScalarFlux(values=rng.random((2, 4)) + 0.5, space=mesh.bulk_space),
        boundary=ScalarBoundaryFlux(values=rng.random(mesh.scalar_trace.shape[0]) + 0.1, space=mesh.scalar_trace),
    )


def _loss(mesh, mat_xs):
    """Assemble A = L + C − S − B; B reads the mesh's own realized laws."""
    ffs = mesh.full_field_space
    L = LeakageOperator(mesh)
    C = MultiplicationOperator(mat_xs.total_cross_section_field, domain=ffs, codomain=ffs)
    # The energy binding is PLAIN-bound on the scalar bulk and LIFTED onto
    # the composite (CS4c step 5, R-4) — the diffusion solver's own spelling.
    S = BulkLift(
        IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space),
        domain=ffs, codomain=ffs,
    )
    B = DiffusionBoundaryOperator(mesh)
    return L + C - S - B


# ═══════════════════════════════════════════════════════════════════════
# The independent hand-posed FD matrix (vv L11: plain loops, crosswalk
# conventions hand-transcribed; reads ONLY the hand-held XS tables +
# the mesh's geometric bookkeeping)
# ═══════════════════════════════════════════════════════════════════════

_NG, _NX = 2, 4
_N_BULK = _NG * _NX
# Trace flat layout, hand-transcribed from the crosswalk storage table
# ((2, ng) slots C-raveled; faces in axis-ascending min/max order):
_JP = {"xmin": [8, 9], "xmax": [12, 13]}     # J⁺ columns per group
_JM = {"xmin": [10, 11], "xmax": [14, 15]}   # J⁻ columns per group
_EDGE = {"xmin": 0, "xmax": _NX - 1}


def _hand_posed_loss(albedo_by_face: dict[str, float]) -> np.ndarray:
    """[A] = [L + C − S − B] by plain loops from the fixture tables."""
    h = np.diff(_EDGES)
    volumes = h.copy()                        # slab: unit area
    sig_t = np.stack([_SIG_T_A, _SIG_T_B])    # [mat, g]
    sig_s = np.stack([_SIG_S_A, _SIG_S_B])    # [mat, g_from, g_to]
    # D = 1/(3Σ_tr); the fixture mixtures carry no P1 moment, so
    # Σ_tr = Σ_t EXACTLY (the P1 seam's isotropic limit).
    D = 1.0 / (3.0 * sig_t)                   # [mat, g]

    def b(g: int, i: int) -> int:
        return g * _NX + i

    A = np.zeros((16, 16))
    for g in range(_NG):
        # Interior-face conductances (series half-cell resistances).
        for f in range(1, _NX):               # faces between cells f-1, f
            m_l, m_r = _MAT_IDS[f - 1], _MAT_IDS[f]
            r = h[f - 1] / (2 * D[m_l, g]) + h[f] / (2 * D[m_r, g])
            gf = 1.0 / r
            for i, j in ((f - 1, f), (f, f - 1)):
                A[b(g, i), b(g, j)] -= gf / volumes[i]
                A[b(g, i), b(g, i)] += gf / volumes[i]
        # Bulk reaction terms: +σ_t (C) − σ_s0ᵀ (S).
        for i in range(_NX):
            m = _MAT_IDS[i]
            A[b(g, i), b(g, i)] += sig_t[m, g]
            for g_from in range(_NG):
                A[b(g, i), b(g_from, i)] -= sig_s[m, g_from, g]
        # Boundary faces: edge trace coupling + the two trace rows.
        for face, alb in albedo_by_face.items():
            e = _EDGE[face]
            m = _MAT_IDS[e]
            # Edge bulk row reads the net OUTWARD current (J⁺ − J⁻)/V.
            A[b(g, e), _JP[face][g]] += 1.0 / volumes[e]
            A[b(g, e), _JM[face][g]] -= 1.0 / volumes[e]
            # Outflow-definition defect row: J⁺ − c_φ φ_e − c_J J⁻.
            rho = h[e] / (2 * D[m, g])
            c_phi = 1.0 / (rho + 2.0)
            c_j = (rho - 2.0) / (rho + 2.0)
            A[_JP[face][g], _JP[face][g]] = 1.0
            A[_JP[face][g], b(g, e)] = -c_phi
            A[_JP[face][g], _JM[face][g]] = -c_j
            # Inflow row: J⁻ − 𝒜 J⁺ (L's identity − B's albedo).
            A[_JM[face][g], _JM[face][g]] = 1.0
            A[_JM[face][g], _JP[face][g]] = -alb
    return A


# Config → (BC tags declared on the mesh, expected 𝒜 per face in the
# hand-posed matrix). The laws travel through DiffusionMesh
# construction — the production path, not a parallel realization.
_CONFIGS = {
    "zero_flux/zero_flux": (
        (BC("zero_flux"), BC("zero_flux")),
        {"xmin": -1.0, "xmax": -1.0},
    ),
    "reflective/marshak_vacuum": (
        (BC("reflective"), BC("vacuum")),
        {"xmin": 1.0, "xmax": 0.0},
    ),
    "albedo(0.3)/reflective": (
        (BC("albedo", {"albedo": 0.3}), BC("reflective")),
        {"xmin": 0.3, "xmax": 1.0},
    ),
}


def _config_setup(config: str):
    """(mesh, mat_xs, template, albedos) for one BC config — the
    composite template must live on the config's own mesh instance."""
    (bc_left, bc_right), albedos = _CONFIGS[config]
    mesh = _diffusion_mesh(bc_left, bc_right)
    template = FullField.zeros(
        interior=ScalarFlux, boundary=ScalarBoundaryFlux, space=mesh.full_field_space,
    )
    return mesh, mesh.material_xs_field(), template, albedos


# ═══════════════════════════════════════════════════════════════════════
# The stencil gate (vv Mode 12 — pin the OBJECT, not its spectrum)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
class TestStencilGate:
    @pytest.mark.parametrize("config", list(_CONFIGS), ids=list(_CONFIGS))
    def test_assembled_loss_matches_hand_posed_stencil(self, config):
        mesh, mat_xs, template, albedos = _config_setup(config)
        A = _loss(mesh, mat_xs)
        produced = FlattenedOperator(A, template).as_matrix()
        expected = _hand_posed_loss(albedos)
        np.testing.assert_allclose(produced, expected, rtol=1e-13, atol=1e-16)

    def test_matrix_action_equals_typed_action(
        self, mesh, mat_xs, template, flux,
    ):
        """The functor honesty check: [A] @ x.to_flat() == A.apply(x).to_flat().
        (The fixture mesh IS the reflective/marshak_vacuum config.)"""
        A = _loss(mesh, mat_xs)
        M = FlattenedOperator(A, template).as_matrix()
        np.testing.assert_allclose(
            M @ flux.to_flat(), A.apply(flux).to_flat(),
            rtol=1e-12, atol=1e-14,
        )

    # ── Mutation teeth (the plan's three named error classes + the
    #    closure sign): each mutation must RED the stencil gate. The
    #    monkeypatch is in-process; mesh AND operators are built AFTER
    #    patching (process-discipline: never git-checkout to revert). ──

    def _stencil_delta(self) -> float:
        mesh, mat_xs, template, albedos = _config_setup("zero_flux/zero_flux")
        A = _loss(mesh, mat_xs)
        produced = FlattenedOperator(A, template).as_matrix()
        return float(np.abs(produced - _hand_posed_loss(albedos)).max())

    def test_mutation_d_face_pairing_swap_reds(self, monkeypatch):
        """Crossing the half-cell pairing (h_L with D_R) moves the
        interior conductances on the heterogeneous fixture — O(1) red."""
        import orpheus.diffusion.operators as ops

        def swapped(D, h):
            return 1.0 / (h[:-1] / (2.0 * D[:, 1:]) + h[1:] / (2.0 * D[:, :-1]))

        monkeypatch.setattr(ops, "_interior_conductance", swapped)
        assert self._stencil_delta() > 1e-3

    def test_mutation_sigma_a_for_sigma_t_reds(self):
        """Building C from Σ_a instead of Σ_t (the removal-vs-total
        confusion) breaks the in-group cancellation theorem — red."""
        mesh, mat_xs, template, albedos = _config_setup("zero_flux/zero_flux")
        ffs = mesh.full_field_space
        wrong_C = MultiplicationOperator(
            mat_xs.absorption_cross_section_field, domain=ffs, codomain=ffs,
        )
        A = (
            LeakageOperator(mesh) + wrong_C
            - BulkLift(
                IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space),
                domain=ffs, codomain=ffs,
            )
            - DiffusionBoundaryOperator(mesh)
        )
        produced = FlattenedOperator(A, template).as_matrix()
        delta = np.abs(produced - _hand_posed_loss(albedos)).max()
        assert delta > 1e-3

    def test_mutation_scatter_transpose_reds(self, monkeypatch):
        """Swapping the P0 in-scatter kernel for its transpose is
        observable on the ASYMMETRIC fixture Σ_s — red.

        CS4c step-3 re-point: the verb moved from the MaterialXSField
        facade to the kernel field (this sentinel caught the re-route —
        a mutation of the retired arm reddened nothing).
        """
        from orpheus.transport.material_field import TransferMaterialField

        monkeypatch.setattr(
            TransferMaterialField, "add_p0_source",
            TransferMaterialField.add_p0_source_transpose,
        )
        assert self._stencil_delta() > 1e-3

    def test_mutation_boundary_closure_sign_reds(self, monkeypatch):
        """Flipping the c_J sign in the P1 face closure — red."""
        import orpheus.diffusion.operators as ops

        def flipped(D_edge, h_edge):
            rho = h_edge / (2.0 * D_edge)
            return 1.0 / (rho + 2.0), -(rho - 2.0) / (rho + 2.0)

        monkeypatch.setattr(ops, "_boundary_closure", flipped)
        assert self._stencil_delta() > 1e-3


# ═══════════════════════════════════════════════════════════════════════
# Family laws (the plan's foundation gates)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestFamilyLaws:
    def test_column_sum_conservation(self, mesh, mat_xs, template):
        """1ᵀ(C − S) = Σ_a per cell per group — the in-group
        cancellation theorem: removal is DERIVED, never an input."""
        ffs = mesh.full_field_space
        CS = (
            MultiplicationOperator(mat_xs.total_cross_section_field, domain=ffs, codomain=ffs)
            - BulkLift(
                IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space),
                domain=ffs, codomain=ffs,
            )
        )
        M = FlattenedOperator(CS, template).as_matrix()
        bulk = M[:_N_BULK, :_N_BULK]
        col_sums = bulk.sum(axis=0).reshape(_NG, _NX)
        np.testing.assert_allclose(
            col_sums, mat_xs.absorption_cross_section, rtol=1e-13,
        )

    def test_m_matrix_sign_pattern(self):
        """The bulk block of the loss: positive diagonal, non-positive
        off-diagonal (spatial coupling AND group transfer) — the
        M-matrix structure behind flux positivity."""
        mesh, mat_xs, template, _ = _config_setup("zero_flux/zero_flux")
        A = _loss(mesh, mat_xs)
        bulk = FlattenedOperator(A, template).as_matrix()[:_N_BULK, :_N_BULK]
        assert np.all(np.diag(bulk) > 0.0)
        off = bulk - np.diag(np.diag(bulk))
        assert np.all(off <= 1e-15)

    @pytest.mark.parametrize(
        "config", ["zero_flux/zero_flux", "reflective/marshak_vacuum"],
    )
    def test_per_group_spd_of_bulk_schur_complement(self, config):
        """Per-group −∇·D_g∇ + Σ_r,g (the bulk Schur complement's
        group-diagonal block) is symmetric positive definite **under the
        volume metric** — the conservative divergence divides each row
        by ITS cell volume, so the raw matrix is V-similar to the
        symmetric bilinear form on a non-uniform mesh: ``diag(V) @
        block`` is the honest SPD object (exactly the composite bulk
        space's ``inner_product_weights``). The FULL multigroup loss is
        NOT symmetric (down-scatter is lower-triangular) — the SPD
        claim is per-group only."""
        mesh, mat_xs, template, _ = _config_setup(config)
        A = _loss(mesh, mat_xs)
        M = FlattenedOperator(A, template).as_matrix()
        A_bb = M[:_N_BULK, :_N_BULK]
        A_bt = M[:_N_BULK, _N_BULK:]
        A_tb = M[_N_BULK:, :_N_BULK]
        A_tt = M[_N_BULK:, _N_BULK:]
        schur = A_bb - A_bt @ np.linalg.solve(A_tt, A_tb)
        V = np.diag(np.diff(_EDGES))
        for g in range(_NG):
            block = V @ schur[g * _NX:(g + 1) * _NX, g * _NX:(g + 1) * _NX]
            np.testing.assert_allclose(block, block.T, rtol=1e-12, atol=1e-15)
            assert np.all(np.linalg.eigvalsh(0.5 * (block + block.T)) > 0.0), config

    def test_reflective_annihilates_constant_mode_exactly(self):
        """(L − B_reflective) on [φ = const, J± = φ/4] is EXACTLY zero:
        interior currents vanish, the net trace current vanishes, the
        outflow defect closes (c_φ = (1 − c_J)/4 algebraically), and
        the inflow row reads J⁻ − J⁺ = 0."""
        mesh = _diffusion_mesh(BC("reflective"), BC("reflective"))
        L = LeakageOperator(mesh)
        B = DiffusionBoundaryOperator(mesh)
        const = FullField(
            interior=ScalarFlux(values=np.full((_NG, _NX), 3.7), space=mesh.bulk_space),
            boundary=ScalarBoundaryFlux(values=np.full(mesh.scalar_trace.shape[0], 3.7 / 4.0), space=mesh.scalar_trace),
        )
        out = (L - B).apply(const)
        # Bulk is BIT-exact zero (all currents are differences of equal
        # floats); the outflow-defect row closes through the algebraic
        # identity c_φ = (1 − c_J)/4, which is exact in reals and
        # ULP-level in floats.
        np.testing.assert_array_equal(out.interior.values, 0.0)
        np.testing.assert_allclose(out.boundary.values, 0.0, atol=1e-15)

    def test_dense_resolvent_over_flattened_loss(self):
        """The #290 resolvent spelling: MatrixInverseOperator over the
        flattened loss — M-materialise both ways at machine·cond grain,
        and seed-independent apply."""
        mesh, mat_xs, template, _ = _config_setup("zero_flux/zero_flux")
        A = _loss(mesh, mat_xs)
        A_flat = FlattenedOperator(A, template)
        A_inv = MatrixInverseOperator(A_flat)
        M = A_flat.as_matrix()
        Minv = A_inv.as_matrix()
        np.testing.assert_allclose(Minv @ M, np.eye(16), atol=1e-11)
        np.testing.assert_allclose(M @ Minv, np.eye(16), atol=1e-11)
        rng = np.random.default_rng(7)
        q = rng.random(16)
        x0 = A_inv.apply(q)
        x1 = A_inv.apply(q, initial_guess=rng.random(16))
        np.testing.assert_array_equal(np.asarray(x0), np.asarray(x1))

    def test_loss_block_role_joins_to_full(self, mesh, mat_xs):
        A = _loss(mesh, mat_xs)
        assert A.block_role is BlockRole.FULL
        assert isinstance(A, FullOperator)

    def test_composition_guard_validates_shared_space(self, mesh, mat_xs):
        """Every member advertises the SAME scalar composite space, so
        the OperatorSum guard actually validates the build."""
        A = _loss(mesh, mat_xs)
        assert A.domain == mesh.full_field_space
        assert A.codomain == mesh.full_field_space


# ═══════════════════════════════════════════════════════════════════════
# The two new leaves
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestLeakageOperator:
    def test_block_role_full(self, mesh):
        assert LeakageOperator(mesh).block_role is BlockRole.FULL

    # (The multi-D refusal moved to the PHASE SPACE at #290 P7a — a
    # multi-D DiffusionMesh is unrepresentable; gated in
    # test_augmented_mesh.py.)

    def test_role_parses_and_mesh_identity(self, mesh, flux):
        L = LeakageOperator(mesh)
        # Wrong bulk role/family:
        bad_bulk = FullField(
            interior=ScalarSourceSink(values=np.ones((_NG, _NX)), space=mesh.bulk_space),
            boundary=flux.boundary,
        )
        with pytest.raises(TypeError, match="ScalarFlux"):
            L.apply(bad_bulk)
        # CS4b S3 (F2): a twin carrier's composite is content-equal — legal.
        twin = _diffusion_mesh()
        L.apply(FullField(
            interior=ScalarFlux(values=np.ones((_NG, _NX)), space=twin.bulk_space),
            boundary=ScalarBoundaryFlux.zeros(twin.scalar_trace),
        ))
        # A carrier whose volumes differ refuses (space-content invariant).
        other = DiffusionMesh(
            Mesh1D(edges=2.0 * _EDGES, mat_ids=_MAT_IDS),
            _fixture_materials(),
        )
        foreign = FullField(
            interior=ScalarFlux(values=np.ones((_NG, _NX)), space=other.bulk_space),
            boundary=ScalarBoundaryFlux.zeros(other.scalar_trace),
        )
        with pytest.raises(ValueError, match="space-content"):
            L.apply(foreign)

    def test_interior_stencil_on_uniform_slab(self):
        """Uniform-h single-material sanity: (Lφ)_i = −D(φ_{i+1} − 2φ_i
        + φ_{i−1})/h² in the interior — the classic 3-point stencil."""
        mats = _fixture_materials()
        mesh1d = Mesh1D(
            edges=np.linspace(0.0, 4.0, 5), mat_ids=np.zeros(4, dtype=int),
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
        mm = DiffusionMesh(mesh1d, {0: mats[0]})
        L = LeakageOperator(mm)
        rng = np.random.default_rng(3)
        phi = rng.random((_NG, 4))
        psi = FullField(
            interior=ScalarFlux(values=phi, space=mm.bulk_space),
            boundary=ScalarBoundaryFlux.zeros(mm.scalar_trace),
        )
        out = L.apply(psi).interior.values
        D = 1.0 / (3.0 * _SIG_T_A)
        h = 1.0
        for g in range(_NG):
            for i in (1, 2):
                expected = -D[g] * (phi[g, i + 1] - 2 * phi[g, i] + phi[g, i - 1]) / h**2
                np.testing.assert_allclose(out[g, i], expected, rtol=1e-13)

    def test_trace_rows_and_edge_coupling(self, mesh, flux):
        """The trace block hand-checked on one face: outflow defect
        J⁺ − c_φ φ_e − c_J J⁻ and the inflow identity J⁻; the edge bulk
        row reads +net/V."""
        L = LeakageOperator(mesh)
        out = L.apply(flux)
        trace_in = flux.boundary
        phi = flux.interior.values
        h = np.diff(_EDGES)
        D = 1.0 / (3.0 * np.stack([_SIG_T_A, _SIG_T_B]))  # [mat, g]
        for face, e in _EDGE.items():
            m = _MAT_IDS[e]
            rho = h[e] / (2 * D[m])
            c_phi = 1.0 / (rho + 2.0)
            c_j = (rho - 2.0) / (rho + 2.0)
            jp = trace_in.outflow_view(face)
            jm = trace_in.inflow_view(face)
            slot = out.boundary.face_view(face)
            np.testing.assert_allclose(
                slot[0], jp - c_phi * phi[:, e] - c_j * jm, rtol=1e-13,
            )
            np.testing.assert_array_equal(slot[1], jm)


@pytest.mark.foundation
class TestDiffusionBoundaryOperator:
    def test_inflow_row_carries_albedo_action_zero_elsewhere(self):
        mesh = _diffusion_mesh(
            BC("albedo", {"albedo": 0.3}), BC("zero_flux"),
        )
        flux = _random_flux(mesh)
        B = DiffusionBoundaryOperator(mesh)
        out = B.apply(flux)
        np.testing.assert_array_equal(out.interior.values, 0.0)
        for face, alb in (("xmin", 0.3), ("xmax", -1.0)):
            slot = out.boundary.face_view(face)
            np.testing.assert_array_equal(slot[0], 0.0)      # no outflow emission
            np.testing.assert_allclose(
                slot[1], alb * flux.boundary.outflow_view(face), rtol=1e-15,
            )
        assert isinstance(out.boundary, ScalarBoundarySourceSink)

    # (Face-law coverage is STRUCTURAL since #290 P7a — mesh.bc and the
    # trace share the one face_labels inventory; the positive invariant
    # is gated in test_augmented_mesh.py. The pre-P7a mismatched-dict
    # refusals guarded a state that is no longer representable.)

    def test_block_role_boundary(self, mesh):
        B = DiffusionBoundaryOperator(mesh)
        assert B.block_role is BlockRole.BOUNDARY
        assert isinstance(B, BoundaryOperator)


# ═══════════════════════════════════════════════════════════════════════
# The scalar-composite substrate + the shared-operator scalar arms
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestScalarCompositeSubstrate:
    def test_diffusion_coefficient_per_cell_gather(self, mat_xs):
        """D = 1/(3Σ_tr) gathered through mat_map — hand arithmetic
        (the fixture has no P1 moment, so Σ_tr = Σ_t exactly)."""
        D = mat_xs.diffusion_coefficient
        expected = np.empty((_NG, _NX))
        for i, m in enumerate(_MAT_IDS):
            expected[:, i] = 1.0 / (3.0 * (_SIG_T_A if m == 0 else _SIG_T_B))
        np.testing.assert_array_equal(D, expected)

    def test_full_field_space_blocks_and_metric(self, mesh):
        ffs = mesh.full_field_space
        assert ffs.shape == (_NG * _NX + 2 * 2 * _NG,)
        # CS4b S4: the interior IS the carrier's cached axis-built mint;
        # the bulk metric ACTION is cell volumes broadcast over groups
        # (asserted on the action — axis-built spaces store per-axis
        # weights, not one dense array).
        assert ffs.interior_space is mesh.bulk_space
        assert ffs.trace_space is mesh.scalar_trace
        x = np.arange(1.0, 1.0 + float(_NG * _NX)).reshape(_NG, _NX)
        np.testing.assert_array_equal(
            ffs.interior_space.apply_metric(x),
            x * np.diff(_EDGES)[None, :],
        )

    def test_scalar_boundary_source_sink_class_gate(self, mesh):
        """Same trace space, different ROLE: source ± flux is rejected
        by class identity; same-class arithmetic is closed."""
        src = ScalarBoundarySourceSink.zeros(mesh.scalar_trace)
        flx = ScalarBoundaryFlux.zeros(mesh.scalar_trace)
        with pytest.raises(TypeError):
            _ = src + flx  # type: ignore[operator]
        total = src + ScalarBoundarySourceSink.zeros(mesh.scalar_trace)
        assert isinstance(total, ScalarBoundarySourceSink)


@pytest.mark.foundation
class TestSharedOperatorScalarArms:
    def test_multiplication_scalar_apply_and_solve(self, mesh, mat_xs, flux):
        space = mesh.full_field_space
        C = MultiplicationOperator(
            mat_xs.total_cross_section_field, domain=space, codomain=space,
        )
        out = C.apply(flux)
        assert isinstance(out.interior, ScalarSourceSink)
        assert isinstance(out.boundary, ScalarBoundarySourceSink)
        np.testing.assert_array_equal(
            out.interior.values, mat_xs.total_cross_section * flux.interior.values,
        )
        np.testing.assert_array_equal(out.boundary.values, 0.0)
        # solve is the typed division back to a flux composite.
        back = C.solve(out)
        assert isinstance(back.interior, ScalarFlux)
        np.testing.assert_allclose(
            back.interior.values, flux.interior.values, rtol=1e-15,
        )

    def test_multiplication_solve_gates_spectrum_on_scalar_arm(self, mesh, flux):
        zero_coeff = MultiplicationOperator.from_mesh(
            np.zeros((_NG, _NX)), mesh,
        )
        with pytest.raises(NotInvertible):
            zero_coeff.solve(flux)

    def test_fission_scalar_composite_arm_matches_scalar_flux_arm(
        self, mesh, mat_xs, flux,
    ):
        # CS4c step 4: diffusion consumes the fission ENERGY binding
        # (IsotropicFission); since step 5 (R-4) that binding is PLAIN-bound
        # on the scalar bulk and its composite action is the LIFT's — the
        # bare arm is the same dyad, the lift performs no arithmetic.
        F = IsotropicFission.from_material_xs(mat_xs, space=mesh.bulk_space)
        composite = BulkLift(
            F, domain=mesh.full_field_space, codomain=mesh.full_field_space,
        ).apply(flux)
        direct = np.asarray(F.apply(flux.interior.values))
        assert isinstance(composite.interior, ScalarSourceSink)
        np.testing.assert_array_equal(composite.interior.values, direct)
        np.testing.assert_array_equal(composite.boundary.values, 0.0)
        assert isinstance(composite.boundary, ScalarBoundarySourceSink)

    def test_k_iso_composite_arm_matches_bare_kernel(self, mesh, mat_xs, flux):
        ffs = mesh.full_field_space
        for op in (
            IsotropicScattering.from_material_xs(
                mat_xs, space=mat_xs.mesh.bulk_space,
            ),
            IsotropicN2N.from_material_xs(
                mat_xs, space=mat_xs.mesh.bulk_space,
            ),
        ):
            composite = BulkLift(op, domain=ffs, codomain=ffs).apply(flux)
            bare = op.apply(flux.interior.values)
            assert isinstance(composite.interior, ScalarSourceSink)
            np.testing.assert_array_equal(composite.interior.values, bare)
            np.testing.assert_array_equal(composite.boundary.values, 0.0)

    def test_plain_binding_refuses_the_composite_and_the_lift_transposes_it(
        self, mesh, mat_xs, flux,
    ):
        """INVERTED at CS4c step 5. Until then this row pinned the iso
        binding's own composite-transpose REFUSAL (the `#281` note: "lands
        with the adjoint diffusion consumer"). R-4 makes the plain binding
        refuse EVERY composite — apply and transpose alike, naming the lift
        — and the lift's transpose is the plain transpose extended by zero
        on the trace, so the composite adjoint the note deferred now exists
        with no fission-, scatter- or diffusion-specific arithmetic."""
        ffs = mesh.full_field_space
        S = IsotropicScattering.from_material_xs(
            mat_xs, space=mat_xs.mesh.bulk_space,
        )
        with pytest.raises(TypeError, match="BulkLift"):
            S.apply_transpose(flux)
        with pytest.raises(TypeError, match="BulkLift"):
            S.apply(flux)
        lifted = BulkLift(S, domain=ffs, codomain=ffs).apply_transpose(flux)
        assert isinstance(lifted.interior, ScalarSourceSink)
        np.testing.assert_array_equal(
            lifted.interior.values, S.apply_transpose(flux.interior.values),
        )
        np.testing.assert_array_equal(lifted.boundary.values, 0.0)


# ═══════════════════════════════════════════════════════════════════════
# The assembly mode (stencil-assembly 2b — the diffusion family is the
# FIRST emitter consumer). Gates per the L16 spec: G1 assembled@x ≡
# apply(x) (never a scalar functional — Mode 12); G3 = THE diffusion
# family's ONE probed≡assembled pin (the fuller-view-oracle keep
# decision — forcing the RETAINED _as_matrix_by_probing pathway, since
# as_matrix itself now DELEGATES to assembly and would compare assembly
# with itself); a Mode-11 sentinel proving the resolvent path actually
# executes the delegation; and the one-source teeth (a sign flip in the
# SHARED coefficient source must red BOTH absolute gates while the
# cross-mode equivalence PERSISTS — divergence under a shared-source
# mutation is precisely the twin-path signature).
#
# The B / C / S emitters are one-source BY CONSTRUCTION (they extract
# coefficients THROUGH their own production apply/kernels — law probing,
# the coefficient array, group-impulse probing), so the mutation teeth
# target the one DIRECT emitter: LeakageOperator's conductance +
# boundary-closure sources.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l0
class TestAssemblyMode:
    _RTOL = 1e-11   # L16: CSR summation order ≠ apply order — never 0-ULP

    @pytest.mark.parametrize("config", list(_CONFIGS), ids=list(_CONFIGS))
    def test_g1_assembled_matvec_equals_typed_apply(self, config):
        """G1 on the full loss, per BC config: the assembled sparse
        matvec reproduces the typed composite action on a non-flat
        seeded random composite."""
        mesh, mat_xs, template, _ = _config_setup(config)
        A = _loss(mesh, mat_xs)
        x = _random_flux(mesh)
        np.testing.assert_allclose(
            A.assemble().apply(x.to_flat()),
            np.asarray(A.apply(x).to_flat()),
            rtol=self._RTOL, atol=1e-14,
        )

    def test_g1_per_leaf(self, mesh, mat_xs, flux):
        """G1 leaf-localized: each emitter's sparse matvec reproduces
        its own typed action (failure here names the broken emitter
        before the composite gate smears it)."""
        ffs = mesh.full_field_space
        leaves = {
            "L": LeakageOperator(mesh),
            "B": DiffusionBoundaryOperator(mesh),
            "C": MultiplicationOperator(
                mat_xs.total_cross_section_field, domain=ffs, codomain=ffs,
            ),
            # the energy bindings assemble on their plain ends; the lift
            # embeds the emission in the composite layout (G5.6 — owed for
            # exactly these two leaves: IsotropicFission has no assembly)
            "S": BulkLift(
                IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space),
                domain=ffs, codomain=ffs,
            ),
            "N2N": BulkLift(
                IsotropicN2N.from_material_xs(mat_xs, space=mesh.bulk_space),
                domain=ffs, codomain=ffs,
            ),
        }
        for name, leaf in leaves.items():
            np.testing.assert_allclose(
                leaf.assemble().apply(flux.to_flat()),
                np.asarray(leaf.apply(flux).to_flat()),
                rtol=self._RTOL, atol=1e-14,
                err_msg=f"leaf {name}: assembled matvec != typed apply",
            )

    def test_g3_probed_equals_assembled_pin(self, mesh, mat_xs, template):
        """G3 — the diffusion family's ONE permanent probed≡assembled
        oracle: the RETAINED apply-to-basis loop (forced explicitly —
        as_matrix now delegates) against the densified emission."""
        A = _loss(mesh, mat_xs)
        flat = FlattenedOperator(A, template)
        probed = flat._as_matrix_by_probing((flat.n_flat,))
        np.testing.assert_allclose(
            A.assemble().as_matrix(), probed,
            rtol=self._RTOL, atol=1e-15,
        )

    def test_resolvent_materializes_through_assembly(
        self, mesh, mat_xs, template, monkeypatch,
    ):
        """Mode-11 sentinel: the production resolvent construction
        (MatrixInverseOperator over the flattened loss) must EXECUTE the
        assembly delegation — no probing at the COMPOSITE flat dimension
        may fire. (Probes at the tiny law dimension ng ARE expected: the
        B emitter extracts each realized law's (ng × ng) matrix THROUGH
        the law's own apply — that in-emitter probing is the one-source
        discipline, not a delegation escape.) A green equivalence gate
        proves values; only this sentinel proves the consumer is
        actually on the new path."""
        from orpheus.numerics.operator import LinearOperator

        probe_dims: list[tuple[str, int]] = []
        original = LinearOperator._as_matrix_by_probing

        def counting(self, shape):
            probe_dims.append((type(self).__name__, int(np.prod(shape))))
            return original(self, shape)

        monkeypatch.setattr(
            LinearOperator, "_as_matrix_by_probing", counting,
        )
        A = _loss(mesh, mat_xs)
        n_flat = int(np.asarray(template.to_flat()).size)
        resolvent = MatrixInverseOperator(FlattenedOperator(A, template))
        composite_probes = [d for d in probe_dims if d[1] == n_flat]
        assert composite_probes == [], (
            f"the resolvent probed the composite instead of assembling: "
            f"{composite_probes}"
        )
        assert all(dim <= _NG for _, dim in probe_dims), (
            f"unexpected large probe (neither law-sized nor delegated): "
            f"{probe_dims}"
        )
        # And the inverse is real: A⁻¹(A x) round-trips.
        x = _random_flux(mesh).to_flat()
        forward = A.assemble().apply(x)
        np.testing.assert_allclose(
            np.asarray(resolvent.apply(forward)).ravel(), x,
            rtol=1e-10, atol=1e-12,
        )

    # ── One-source teeth (L16): mutate the SHARED source — both
    #    absolute gates red, the cross-mode equivalence persists. ──────

    def _mutated_deltas_and_equivalence(self) -> tuple[float, float]:
        """(assembled-vs-hand-posed delta, apply-vs-hand-posed delta)
        under the CURRENT (possibly monkeypatched) coefficient sources;
        asserts the cross-mode equivalence G1 inside."""
        mesh, mat_xs, template, albedos = _config_setup("zero_flux/zero_flux")
        A = _loss(mesh, mat_xs)
        hand = _hand_posed_loss(albedos)
        x = _random_flux(mesh)
        x_flat = x.to_flat()
        assembled = A.assemble()
        asm_delta = float(np.abs(assembled.as_matrix() - hand).max())
        apply_delta = float(
            np.abs(np.asarray(A.apply(x).to_flat()) - hand @ x_flat).max()
        )
        # The one-source signature: both modes moved TOGETHER.
        np.testing.assert_allclose(
            assembled.apply(x_flat), np.asarray(A.apply(x).to_flat()),
            rtol=self._RTOL, atol=1e-14,
        )
        return asm_delta, apply_delta

    def test_teeth_conductance_sign_flip_reds_both_modes(self, monkeypatch):
        """−g_f in the SHARED interior-conductance source: the assembled
        matrix AND the typed apply BOTH leave the hand-posed reference
        O(1), while assembled@x ≡ apply(x) persists — one coefficient
        source, no twin stencil."""
        import orpheus.diffusion.operators as ops

        original = ops._interior_conductance
        monkeypatch.setattr(
            ops, "_interior_conductance",
            lambda D, h: -original(D, h),
        )
        asm_delta, apply_delta = self._mutated_deltas_and_equivalence()
        assert asm_delta > 1e-2, "assembly blind to the shared-source flip"
        assert apply_delta > 1e-2, "apply blind to the shared-source flip"

    def test_teeth_boundary_closure_sign_flip_reds_both_modes(
        self, monkeypatch,
    ):
        """−c_J in the SHARED P1 face-closure source: same one-source
        signature on the trace-row emission family."""
        import orpheus.diffusion.operators as ops

        def flipped(D_edge, h_edge):
            rho = h_edge / (2.0 * D_edge)
            return 1.0 / (rho + 2.0), -(rho - 2.0) / (rho + 2.0)

        monkeypatch.setattr(ops, "_boundary_closure", flipped)
        asm_delta, apply_delta = self._mutated_deltas_and_equivalence()
        assert asm_delta > 1e-2, "assembly blind to the closure flip"
        assert apply_delta > 1e-2, "apply blind to the closure flip"

    def test_plain_bulk_leaves_assemble_on_their_own_layout(self, mat_xs):
        """INVERTED at CS4c step 5. Until then this row pinned the plain
        bulk bindings' assemble REFUSAL ("no composite flat layout to emit
        into"). The plain binding now emits on its OWN layout — the
        domain's C-ravel, ``n_bulk × n_bulk`` — and the composite
        embedding is the lift's job. Surviving claim: each plain leaf is
        assemblable and its emission's matvec reproduces its ``apply`` on
        a bare array of the bound shape (0 ULP for the diagonal, the L16
        CSR band for the block-diagonal energy pair)."""
        coef = mat_xs.total_cross_section_field
        x = np.random.default_rng(3).random(coef.values.shape) + 0.5
        n = x.size
        for op, rtol in (
            (MultiplicationOperator(coef, domain=coef.space, codomain=coef.space), 0.0),
            (IsotropicScattering.from_material_xs(mat_xs, space=mat_xs.mesh.bulk_space), 1e-11),
            (IsotropicN2N.from_material_xs(mat_xs, space=mat_xs.mesh.bulk_space), 1e-11),
        ):
            assert op.is_assemblable
            emitted = op.assemble()
            assert emitted.shape == (n, n)
            np.testing.assert_allclose(
                emitted.apply(x.ravel()), np.asarray(op.apply(x)).ravel(),
                rtol=rtol, atol=0.0,
            )

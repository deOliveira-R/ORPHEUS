"""Foundation tests for :class:`orpheus.transport.operators.scattering.ScatteringOperator`.

Round 1.2 of Wave D of the SN reshape campaign (Issue #162). The
operator carries the same math :class:`SNSolver` used to expose under
``_add_scattering_source``, ``_build_aniso_scattering``, and
``_add_n2n_source`` (now thin delegators); these tests pin the lifted
math at the operator level so any drift would be observable here as
well as via the underscore-prefixed delegators.

The load-bearing test is the **bit-identical extraction** suite: a
synthetic ``(psi, phi, Q)`` triple is fed through both the new
:meth:`ScatteringOperator.apply` (and the in-place helpers
:meth:`add_iso_source` / :meth:`add_n2n_source` / :meth:`build_aniso_source`)
and the explicit per-cell reference implementations from
``test_solver_components.py``. The two paths must agree to round-off,
because the operator is a structural extraction, not a re-derivation.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture, make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.operator import LinearOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import material_xs_from_raw
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation  # software-invariant tier


class _StubQuad:
    """Minimal :class:`~orpheus.numerics.quadrature.Quadrature`
    stand-in for synthetic tests.

    Issue #197 PR-TYPED-1 — :class:`ScatteringOperator` now consumes
    a quadrature handle (instead of separate ``N`` / ``weights`` / ``Y``
    constructor args).  Tests that exercise the per-material dispatch
    in isolation use this stub to avoid building a full quadrature.
    """

    def __init__(self, *, N: int, weights: np.ndarray) -> None:
        self.N = N
        self.weights = np.asarray(weights)

    def spherical_harmonics(self, L: int) -> np.ndarray:
        # Only called when scattering_order > 0; synthetic tests
        # below pass scattering_order=0 or =1 with anisotropy unused.
        return np.zeros((self.N, L + 1, 2 * L + 1))


def _uniform_2d(nx, ny, delta, mat_map):
    """Helper: build a uniform Mesh2D."""
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def solver_2g_p0():
    """SNSolver fixture, 2-group, P0 only (no anisotropic data)."""
    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    materials = {2: fuel, 0: mod}

    nx, ny = 6, 4
    delta = 0.2
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0

    mesh = _uniform_2d(nx, ny, delta, mat)
    quad = Quadrature.lebedev(order=17)
    sn_mesh = SNMesh(mesh, quad, materials)
    solver = SNSolver(sn_mesh)
    return solver


@pytest.fixture
def solver_2g_p0_n2n():
    """SNSolver fixture, 2-group, P0 only, with a NON-ZERO (n,2n) matrix.

    #269 (Mode-10 cure): the (n,2n)-term tests below must NOT ride a
    library mixture, whose ``Sig2`` is all-zero — that nulls the exact
    term they verify (``2·Σ_2n``), so a sign/factor mutation in the
    (n,2n) path cannot move the measured value (vacuous green).  This
    fixture builds the SAME 2-group, P0-only, fuel/moderator geometry as
    :func:`solver_2g_p0` but injects an asymmetric ``Sig2`` into the fuel
    material so the term is genuinely activated AND constrained.  ``Sig2``
    is cross-group only (the asymmetry catches a ``g_from↔g_to`` swap).
    """
    from scipy.sparse import csr_matrix

    fuel = get_mixture("A", "2g")
    mod = get_mixture("B", "2g")
    # Inject a non-zero, asymmetric (n,2n) matrix into the fuel mixture.
    # SigT must absorb the (n,2n) reaction XS so the mixture stays balanced.
    sig2 = np.array([[0.0, 0.03], [0.01, 0.0]])
    fuel.Sig2 = csr_matrix(sig2)
    fuel.SigT = np.asarray(fuel.SigT) + sig2.sum(axis=1)
    materials = {2: fuel, 0: mod}

    nx, ny = 6, 4
    delta = 0.2
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2
    mat[3:, :] = 0

    mesh = _uniform_2d(nx, ny, delta, mat)
    quad = Quadrature.lebedev(order=17)
    sn_mesh = SNMesh(mesh, quad, materials)
    return SNSolver(sn_mesh)


# ── Reference implementations (per-cell loops, known correct) ─────────


def _ref_iso_scatter_inplace(solver, Q, phi):
    """Reference per-cell P0 in-scatter (bit-identical to the legacy code).

    Issue #196 PR-INDEX-4: ``Q`` / ``phi`` are principled
    ``(ng, nx, ny)``.  Per-cell update reads ``(ng,)`` slices over the
    spatial pair.
    """
    out = Q.copy()
    nx, ny = solver.sn_mesh.spatial_shape
    for ix in range(nx):
        for iy in range(ny):
            mid = int(solver.sn_mesh.mat_map[ix, iy])
            out[:, ix, iy] += {mid: solver.mat_xs.sig_s_legendre(mid)[0] for mid in solver.mat_xs.materials}[mid].T @ phi[:, ix, iy]
    return out


def _ref_n2n_inplace(solver, Q, phi):
    """Reference per-cell (n,2n) source (bit-identical to the legacy code).

    Issue #196 PR-INDEX-4 — principled ``(ng, nx, ny)`` (see
    :func:`_ref_iso_scatter_inplace`).
    """
    out = Q.copy()
    nx, ny = solver.sn_mesh.spatial_shape
    for ix in range(nx):
        for iy in range(ny):
            mid = int(solver.sn_mesh.mat_map[ix, iy])
            out[:, ix, iy] += 2.0 * ({mid: solver.mat_xs.n2n_matrix(mid) for mid in solver.mat_xs.materials}[mid].T @ phi[:, ix, iy])
    return out


# ──────────────────────────────────────────────────────────────────────
# Protocol contract
# ──────────────────────────────────────────────────────────────────────


class TestProtocolCompliance:
    """ScatteringOperator must satisfy the LinearOperator Protocol."""

    def test_implements_linear_operator(self, solver_2g_p0):
        """isinstance LinearOperator (runtime-checkable Protocol)."""
        assert isinstance(solver_2g_p0.scattering_op, LinearOperator)

    def test_predicates_adjointable_not_invertible(self, solver_2g_p0):
        """``is_adjointable`` True, ``is_invertible`` False — the adjoint S† is free
        via full_scatter_kernel (#276 A2b / #118); still no useful inverse."""
        op = solver_2g_p0.scattering_op
        assert op.is_adjointable and not op.is_invertible

    def test_apply_accepts_psi_shape(self, solver_2g_p0):
        """apply(psi) must accept typed AngularFlux ``(N, ng, nx, ny)`` (D-I.2)."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi_values = np.ones((N, op.ng, *op.spatial_shape))
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p0.sn_mesh)
        out = op.apply(psi)
        assert out.values.shape == psi.values.shape


# ──────────────────────────────────────────────────────────────────────
# Bit-identical extraction (load-bearing)
# ──────────────────────────────────────────────────────────────────────


class TestBitIdenticalExtractionP0:
    """The lifted math must match the legacy reference per-cell code."""

    @pytest.mark.sentinel
    def test_add_iso_source_matches_reference(self, solver_2g_p0):
        """ScatteringOperator.add_iso_source = the per-cell reference.

        Issue #196 PR-INDEX-4: principled ``(ng, nx, ny)`` end-to-end.
        """
        np.random.seed(42)
        (nx, ny), ng = solver_2g_p0.sn_mesh.spatial_shape, solver_2g_p0.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q = np.random.rand(ng, nx, ny)

        expected = _ref_iso_scatter_inplace(solver_2g_p0, Q, phi)

        Q_actual = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q_actual, phi)

        np.testing.assert_allclose(Q_actual, expected, rtol=1e-13)

    def test_add_n2n_source_matches_reference(self, solver_2g_p0_n2n):
        """ScatteringOperator.add_n2n_source = the per-cell reference.

        #269 (Mode-10 cure): rides ``solver_2g_p0_n2n`` (NON-zero,
        asymmetric ``Sig2`` in the fuel) rather than the library
        ``solver_2g_p0`` (``Sig2 = 0``), so the ``2·Σ_2n`` term is
        genuinely constrained — a sign/factor mutation in
        ``add_n2n_source`` reddens this test (see the in-process
        monkeypatch proof in #269 closeout).  The reference
        :func:`_ref_n2n_inplace` is a structurally-independent per-cell
        loop (explicit ``2·Σ_2nᵀ@φ``), not the SUT's reduction.
        """
        solver = solver_2g_p0_n2n
        np.random.seed(123)
        (nx, ny), ng = solver.sn_mesh.spatial_shape, solver.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q = np.random.rand(ng, nx, ny)

        expected = _ref_n2n_inplace(solver, Q, phi)

        Q_actual = Q.copy()
        solver.scattering_op.add_n2n_source(Q_actual, phi)

        np.testing.assert_allclose(Q_actual, expected, rtol=1e-13)

    @pytest.mark.sentinel
    @pytest.mark.parametrize("trailing", [(), (4,)], ids=["scalar", "LD-2^d=4"])
    def test_isotropic_kernel_bit_identical_to_legacy_verbs(
        self, solver_2g_p0_n2n, trailing,
    ):
        r"""#276 P2 — the production ``isotropic_kernel`` (:math:`\Sigma_{s,0} +
        2\Sigma_{2n}`, the :class:`~orpheus.numerics.operator.OperatorSum` that
        :meth:`~orpheus.transport.operators.scattering.ScatteringOperator._assemble_per_ordinate_source`
        now routes the SN forward isotropic source through) is **bit-identical**
        (0-ULP) to the legacy ``add_iso_source`` then ``add_n2n_source`` in-place
        accumulation.

        Both reach the SAME per-material ``mat_xs`` verbs, so the K_iso routing
        introduces zero numerical change — this is the 0-ULP-inheritance anchor
        for the #276 P2 forward re-expression (the structural CORRECTNESS of the
        verbs themselves is pinned independently by
        :meth:`test_add_iso_source_matches_reference` /
        :meth:`test_add_n2n_source_matches_reference` against the per-cell
        reference loops).  Non-zero (n,2n) fixture (#269) + scalar AND LD
        (trailing :math:`2^d` spectator axis, #240 D5b-S3).  ``-O``-safe
        (``np.testing``).
        """
        op = solver_2g_p0_n2n.scattering_op
        rng = np.random.default_rng(0)
        (nx, ny), ng = solver_2g_p0_n2n.sn_mesh.spatial_shape, op.ng
        phi = rng.uniform(0.1, 1.0, size=(ng, nx, ny, *trailing))

        # Production path: the cached OperatorSum the forward now uses.
        got = op.isotropic_kernel.apply(phi)

        # Legacy path: zeros → add_iso_source → add_n2n_source (raw-in, in place).
        ref = np.zeros_like(phi)
        op.add_iso_source(ref, phi)
        op.add_n2n_source(ref, phi)

        np.testing.assert_array_equal(
            got, ref,
            err_msg="isotropic_kernel.apply (#276 P2 production iso path) must "
            "equal the legacy add_iso_source + add_n2n_source accumulation (0-ULP).",
        )

    def test_zero_flux_zero_addition(self, solver_2g_p0):
        """φ = 0 => ScatteringOperator adds zero (linearity guard)."""
        (nx, ny), ng = solver_2g_p0.sn_mesh.spatial_shape, solver_2g_p0.ng
        Q = np.ones((ng, nx, ny))
        phi = np.zeros_like(Q)
        Q_before = Q.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q, phi)
        np.testing.assert_array_equal(Q, Q_before)

    def test_delegators_match_operator_directly(self, solver_2g_p0):
        """SNSolver._add_scattering_source delegates to op.add_iso_source bit-identically.

        Issue #196 PR-INDEX-5: delegator's PUBLIC contract is now
        principled ``(ng, nx, ny)``.  No bridge.
        """
        np.random.seed(7)
        (nx, ny), ng = solver_2g_p0.sn_mesh.spatial_shape, solver_2g_p0.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q_init = np.random.rand(ng, nx, ny)

        Q_via_delegator = Q_init.copy()
        solver_2g_p0._add_scattering_source(Q_via_delegator, phi)

        Q_via_operator = Q_init.copy()
        solver_2g_p0.scattering_op.add_iso_source(Q_via_operator, phi)

        np.testing.assert_array_equal(Q_via_delegator, Q_via_operator)

    def test_delegator_n2n_matches_operator_directly(self, solver_2g_p0_n2n):
        """SNSolver._add_n2n_source delegates to op.add_n2n_source bit-identically.

        Issue #196 PR-INDEX-5: principled ``(ng, nx, ny)`` end-to-end.
        #269 (Mode-10 cure): rides the NON-zero ``Sig2`` fixture so the
        delegated quantity is a genuinely non-trivial (n,2n) source — a
        delegator that dropped/altered the (n,2n) term would now diverge
        from the operator's output (library ``Sig2=0`` made both sides
        the un-mutated input, hiding any term drift).
        """
        solver = solver_2g_p0_n2n
        np.random.seed(11)
        (nx, ny), ng = solver.sn_mesh.spatial_shape, solver.ng
        phi = np.random.rand(ng, nx, ny) + 0.1
        Q_init = np.random.rand(ng, nx, ny)

        Q_via_delegator = Q_init.copy()
        solver._add_n2n_source(Q_via_delegator, phi)

        Q_via_operator = Q_init.copy()
        solver.scattering_op.add_n2n_source(Q_via_operator, phi)

        np.testing.assert_array_equal(Q_via_delegator, Q_via_operator)


# ──────────────────────────────────────────────────────────────────────
# Pℓ Galerkin reconstruction
# ──────────────────────────────────────────────────────────────────────


class TestAnisotropicScatteringExtraction:
    """build_aniso_source must reproduce the legacy _build_aniso_scattering."""

    @pytest.fixture
    def solver_2g_p1(self):
        """4-group with P1 anisotropic scattering data."""
        # Use 421-group library which carries P1; if not available, skip.
        try:
            mix = get_mixture("A", "4g")
        except Exception:
            pytest.skip("4g library unavailable")
        if len(mix.SigS) < 2:
            pytest.skip("No P1 data in 4g library")

        mesh = _uniform_2d(2, 2, 0.5, np.zeros((2, 2), dtype=int))
        quad = Quadrature.lebedev(order=17)
        return SNSolver(SNMesh(mesh, quad, {0: mix}), scattering_order=1)

    def test_returns_none_for_p0(self, solver_2g_p0):
        """L=0 => Pℓ contribution is None (signal: no aniso source needed)."""
        N = solver_2g_p0.quad.N
        (nx, ny), ng = solver_2g_p0.sn_mesh.spatial_shape, solver_2g_p0.ng
        # D-I.2: typed AngularFlux carrier.
        psi_values = np.ones((N, ng, nx, ny))
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p0.sn_mesh)
        out = solver_2g_p0.scattering_op.build_aniso_source(psi)
        assert out is None

    def test_returns_none_for_no_angular_flux(self, solver_2g_p0):
        """First-iteration sentinel: psi=None => return None."""
        out = solver_2g_p0.scattering_op.build_aniso_source(None)
        assert out is None

    def test_isotropic_flux_zero_aniso_source(self, solver_2g_p1):
        """Isotropic ψ_n = const for every ordinate => P1+ Galerkin moments = 0."""
        op = solver_2g_p1.scattering_op
        N = op.n_ordinates
        psi_iso_values = np.ones((N, op.ng, *op.spatial_shape))
        psi_iso = AngularFlux.from_mesh(psi_iso_values, solver_2g_p1.sn_mesh)
        Q_aniso = op.build_aniso_source(psi_iso)
        assert Q_aniso is not None
        np.testing.assert_allclose(Q_aniso.values, 0, atol=1e-12)

    def test_delegator_matches_operator(self, solver_2g_p1):
        """SNSolver._build_aniso_scattering delegates bit-identically.

        D-I.2 — the delegator preserves a bare-ndarray external
        contract for backward compat (solver.py:1203 caller); the
        operator now consumes typed :class:`AngularFlux` and returns
        :class:`AngularSourceSink`.  Compare ``out_via_delegator``
        (bare ndarray) to ``out_via_operator.values`` (typed
        :class:`AngularSourceSink` unwrapped).
        """
        op = solver_2g_p1.scattering_op
        N = op.n_ordinates
        np.random.seed(42)
        psi_values = np.random.rand(N, op.ng, *op.spatial_shape) + 0.1
        psi_typed = AngularFlux.from_mesh(psi_values, solver_2g_p1.sn_mesh)

        out_via_delegator = solver_2g_p1._build_aniso_scattering(psi_values)
        out_via_operator = op.build_aniso_source(psi_typed)
        np.testing.assert_array_equal(out_via_delegator, out_via_operator.values)


# ──────────────────────────────────────────────────────────────────────
# apply() — the LinearOperator surface
# ──────────────────────────────────────────────────────────────────────


class TestApplySemantics:
    """apply(psi) returns the per-ordinate scattering source.

    Combines P0 in-scatter + (n,2n) (broadcast across N) + Pℓ (genuine
    per-ordinate) into a single ``(N, ng, nx, ny)`` array (principled
    storage; see :ref:`theory-sn-index-convention`).
    """

    def test_apply_isotropic_flux_p0_only(self, solver_2g_p0):
        """For P0-only solver, apply(ψ) = (P0 in-scatter + (n,2n))(φ) / W broadcast.

        R-1 Step 4 A1 — ``ScatteringOperator.apply`` returns per-ordinate
        density at the producer boundary (the ``/sum_w`` projection
        lives at the producer per Pattern 7).  Pre-A1 this returned
        iso magnitude broadcast across N; post-A1 it returns
        ``Q_iso / sum_w`` broadcast — the per-ordinate value each
        ordinate sees in the per-ordinate transport equation
        ``(Ω·∇ + Σ_t) ψ_n = Q_iso/W + …``.

        D-I.2: typed AngularFlux carrier → AngularSourceSink output.
        """
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        (nx, ny), ng = op.spatial_shape, op.ng

        np.random.seed(5)
        psi_values = np.random.rand(N, ng, nx, ny) + 0.1
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p0.sn_mesh)

        # Compute scalar flux the same way apply() does internally.
        phi = np.einsum('n,ngxy->gxy', op.weights, psi_values)

        # Reference: compute Q_iso explicitly, then project to
        # per-ordinate via /sum_w (R-1 Step 4 A1).
        Q_iso = np.zeros((ng, nx, ny))
        op.add_iso_source(Q_iso, phi)
        op.add_n2n_source(Q_iso, phi)
        sum_w = float(op.weights.sum())
        expected = np.broadcast_to(
            (Q_iso / sum_w)[None, :, :, :], psi_values.shape,
        )

        actual = op.apply(psi)
        np.testing.assert_allclose(actual.values, expected, rtol=1e-13)

    def test_apply_zero_psi_returns_zero(self, solver_2g_p0):
        """ψ = 0 => S·ψ = 0 (linearity guard)."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi_values = np.zeros((N, op.ng, *op.spatial_shape))
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p0.sn_mesh)
        out = op.apply(psi)
        np.testing.assert_array_equal(out.values, np.zeros_like(psi_values))

    def test_apply_linearity(self, solver_2g_p0):
        """S is linear — homogeneity and DIRECT additivity.

        Since campaign 1 CS3 (flux lives in V) the textbook laws are
        directly spellable: ``S(c·ψ) = c·S(ψ)`` and
        ``S(ψ₁+ψ₂) = S(ψ₁)+S(ψ₂)``. The additivity row alone reds an
        affine S (the pre-CS3 blend spelling
        ``S(ψ₁+λ(ψ₂⊖ψ₁)) = (1−λ)S(ψ₁)+λS(ψ₂)`` could not — affine maps
        preserve affine combinations; the sharpness argument is in
        ``test_declared_law_is_linear.py``). ``op.apply`` stays on flux
        states (its domain — S guards it, rejecting a non-flux input)."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        (nx, ny), ng = op.spatial_shape, op.ng
        m = solver_2g_p0.sn_mesh

        np.random.seed(13)
        psi1 = AngularFlux.from_mesh(np.random.rand(N, ng, nx, ny) + 0.1, m)
        psi2 = AngularFlux.from_mesh(np.random.rand(N, ng, nx, ny) + 0.1, m)
        c = 2.5

        np.testing.assert_allclose(
            op.apply(c * psi1).values, (c * op.apply(psi1)).values,
            rtol=1e-12, atol=1e-13,
        )
        lhs = op.apply(psi1 + psi2)
        rhs = op.apply(psi1) + op.apply(psi2)
        np.testing.assert_allclose(lhs.values, rhs.values, rtol=1e-12, atol=1e-13)


# ──────────────────────────────────────────────────────────────────────
# Producer-side normalisation invariant (R-1 Step 4 A1 — Pattern 7)
# ──────────────────────────────────────────────────────────────────────


class TestProducerSideNormalisation:
    """Producer-side ``/sum_w`` invariant on the typed ``apply`` boundary.

    R-1 Step 4 A1 lifted the per-ordinate projection from the sweep
    interior to the producer boundary (Pattern 7 per
    ``coding-elegance`` SKILL.md).  This test class pins the algebraic
    identity that makes the projection ``sum_w``-independent for uniform
    input: under a uniform AngularFlux ``ψ_n = c`` (so the iso magnitude
    is ``φ = c · sum_w``), the producer's apply returns per-ordinate
    density ``Q_n = c · Σ_{g'} (Σ_{s,0}[g'→g] + 2·Σ_{2n}[g'→g])`` —
    explicitly *without* the ``sum_w`` factor.  If any future refactor
    re-introduces a sweep-internal ``/W`` (or drops the producer's
    ``/sum_w``), this test fails.
    """

    @pytest.mark.l0
    @pytest.mark.verifies("matrix-eigenvalue")
    def test_typed_apply_returns_per_ordinate_already_normalised(
        self, solver_2g_p0,
    ):
        r"""Uniform :math:`\psi_n = c` ⇒ producer-side per-ord output is
        :math:`Q_n = c \sum_{g'}(\Sigma_{s,0}[g'\to g] + 2\Sigma_{2n}[g'\to g])`.

        Algebra:

        * :math:`\phi_g = \sum_n w_n\,c = c \cdot \mathrm{sum\_w}` (iso scalar).
        * Iso magnitude :math:`Q_{\rm iso}[g] = \sum_{g'} \big(\Sigma_{s,0}[g'\to g] + 2\Sigma_{2n}[g'\to g]\big) \cdot c \cdot \mathrm{sum\_w}`.
        * Producer-side :math:`1/W`: :math:`Q_n[g] = Q_{\rm iso}[g]/\mathrm{sum\_w} = c \sum_{g'}(\Sigma_{s,0}[g'\to g] + 2\Sigma_{2n}[g'\to g])`.

        The ``sum_w`` factor cancels by construction.  This is the
        load-bearing producer-side identity that Pattern 7 introduces.
        """
        solver = solver_2g_p0
        op = solver.scattering_op
        N = op.n_ordinates
        (nx, ny), ng = op.spatial_shape, op.ng

        c = 0.37
        psi_values = np.full((N, ng, nx, ny), c)
        psi = AngularFlux.from_mesh(psi_values, solver.sn_mesh)

        # Reference: compute Σ_{g'}(Σ_{s,0}[g'→g] + 2·Σ_{2n}[g'→g]) at each
        # cell from the cell's material.  This is the per-ord magnitude
        # the producer must emit (sum_w-independent).
        expected = np.zeros((N, ng, nx, ny))
        for ix in range(nx):
            for iy in range(ny):
                mid = int(solver.sn_mesh.mat_map[ix, iy])
                sig_s0 = solver.mat_xs.sig_s_legendre(mid)[0]   # (ng, ng) — [g'→g]
                sig_2n = solver.mat_xs.n2n_matrix(mid)          # (ng, ng) — [g'→g]
                # Σ_{g'} (Σ_s0[g'→g] + 2·Σ_2n[g'→g]) per target group.
                # sig_s0.T @ ones gives column sums (over g') indexed by g.
                col_sum = (sig_s0.T + 2.0 * sig_2n.T) @ np.ones(ng)
                expected[:, :, ix, iy] = c * col_sum[None, :]

        actual = op.apply(psi)
        np.testing.assert_allclose(actual.values, expected, rtol=1e-13, atol=1e-13)


# ──────────────────────────────────────────────────────────────────────
# D-H.2-C1 — Composite TimedFullField invariants (volumetric scattering).
#
# Per Option β3 (Issue #208), scattering is volumetric — the output
# boundary is the implicit-zero :class:`AngularBoundaryFlux`.  Bulk follows
# the full :math:`P_\ell` Galerkin path identical to the legacy
# :class:`AngularFlux` branch.  The parity tests vs. legacy retired
# with D-H.2-C1 (the legacy class itself retires in C5; both branches
# share the same Galerkin kernel — composite-branch exercise alone is
# sufficient).
# ──────────────────────────────────────────────────────────────────────


class TestCompositeInvariants:
    """Composite :class:`TimedFullField` variant: bulk-only scattering."""

    def test_returns_timeless_full_field(self, solver_2g_p0):
        """Composite input → TIMELESS composite output (#257 S8a base arrow)."""
        from dataclasses import replace

        from orpheus.transport.fields.angular_flux import (
            AngularFlux,
        )
        from orpheus.transport.full_field import FullField
        from orpheus.transport.source_sinks import AngularSourceSink
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = solver_2g_p0.sn_mesh
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        np.random.seed(41)
        bulk_values = np.random.rand(*state.interior.values.shape) + 0.1
        state = replace(state, interior=replace(state.interior, values=bulk_values))

        out = solver_2g_p0.scattering_op.apply(state)

        # #257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``,
        # so the output is the TIMELESS FullField (history-free).
        assert isinstance(out, FullField)
        assert not isinstance(out, TimedFullField)
        assert isinstance(out.interior, AngularSourceSink)
        assert out.interior.mesh is sn_mesh

    def test_implicit_zero_boundary(self, solver_2g_p0):
        """Scattering is volumetric — boundary member is all zeros."""
        from dataclasses import replace

        sn_mesh = solver_2g_p0.sn_mesh
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh)
        np.random.seed(42)
        bulk_values = np.random.rand(*state.interior.values.shape) + 0.1
        state = replace(state, interior=replace(state.interior, values=bulk_values))

        out = solver_2g_p0.scattering_op.apply(state)

        # Implicit-zero boundary (Option β3 / Wave O #208).
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    def test_zero_bulk_zero_output(self, solver_2g_p0):
        """ψ = 0 ⇒ S·ψ = 0 (linearity guard at composite layer)."""
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=solver_2g_p0.sn_mesh)
        out = solver_2g_p0.scattering_op.apply(state)
        np.testing.assert_array_equal(out.interior.values, 0.0)
        np.testing.assert_array_equal(out.boundary.values, 0.0)

    def test_output_is_timeless_full_field(self, solver_2g_p0):
        """#257 S8a — the matvec leaf is a base arrow ``FullField -> FullField``.

        The output is the TIMELESS FullField (history-free) regardless of the
        input iterate's ``history_depth`` (was: the old convention stamped
        ``history_depth`` onto the output — re-pointed).
        """
        from orpheus.transport.full_field import FullField
        from orpheus.transport.timed_full_field import TimedFullField

        sn_mesh = solver_2g_p0.sn_mesh
        for depth in (0, 1, 2, 4):
            state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh, history_depth=depth)
            out = solver_2g_p0.scattering_op.apply(state)
            assert isinstance(out, FullField)
            assert not isinstance(out, TimedFullField)


# ──────────────────────────────────────────────────────────────────────
# Algebraic identities (P0 + (n,2n))
# ──────────────────────────────────────────────────────────────────────


class TestP0AlgebraicIdentities:
    """Hand-checkable cases for the P0 + (n,2n) algebra."""

    def test_p0_uniform_flux_homogeneous(self):
        """Homogeneous medium, uniform φ_g = 1: Q_iso[g] = Σ_g' Σ_s0[g'->g].

        Issue #196 PR-INDEX-4: principled ``(ng, nx, ny)`` ψ / Q.
        """
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]),
            nu=np.array([2.5, 2.5]),
            chi=np.array([1.0, 0.0]),
            sig_s=np.array([[0.38, 0.10], [0.00, 0.90]]),
        )
        nx, ny = 2, 2
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}))
        op = solver.scattering_op

        phi = np.ones((op.ng, nx, ny))
        Q = np.zeros_like(phi)
        op.add_iso_source(Q, phi)

        # Hand-computed: Q[g] = Σ_g' σ_s0[g'->g] · φ[g'] = column-sum · 1.
        sig_s0_dense = np.array(mix.SigS[0].todense())
        # Convention: ORPHEUS ``SigS[l]`` matrix entry ``[g_from, g_to]``.
        # phi @ sig_s0 sums over g_from for each g_to.
        expected_per_cell = np.ones(op.ng) @ sig_s0_dense
        for ix in range(nx):
            for iy in range(ny):
                np.testing.assert_allclose(Q[:, ix, iy], expected_per_cell, rtol=1e-14)

    def test_n2n_doubling_factor(self):
        """For a pure-(n,2n) mixture (Σ_s0 = 0), Q = 2·φ·Σ_2n."""
        # Build a synthetic mixture with zero P0 scatter and known sig2.
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.0, 0.0]),  # no fission
            nu=np.array([0.0, 0.0]),
            chi=np.zeros(2),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
            sig_s=np.array([[0.0, 0.0], [0.0, 0.0]]),  # zero P0
        )
        # Inject a non-zero (n,2n) matrix manually after construction.
        from scipy.sparse import csr_matrix
        sig2_test = np.array([[0.0, 0.05], [0.0, 0.0]])
        mix.Sig2 = csr_matrix(sig2_test)

        nx, ny = 2, 2
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}))
        op = solver.scattering_op

        np.random.seed(31)
        # Issue #196 PR-INDEX-4: principled (ng, nx, ny).
        phi = np.random.rand(op.ng, nx, ny) + 0.1
        Q = np.zeros_like(phi)
        op.add_iso_source(Q, phi)
        # P0 contribution should be zero
        np.testing.assert_allclose(Q, 0, atol=1e-15)

        # (n,2n) contribution
        op.add_n2n_source(Q, phi)
        # Hand-computed: Q[g, ix, iy] = 2 · sum_g' phi[g', ix, iy] · sig2[g'->g]
        for ix in range(nx):
            for iy in range(ny):
                expected = 2.0 * phi[:, ix, iy] @ sig2_test
                np.testing.assert_allclose(Q[:, ix, iy], expected, rtol=1e-14)


# ──────────────────────────────────────────────────────────────────────
# Foldable / residual split — Phase G Step 3+4.a (Issue #196)
# ──────────────────────────────────────────────────────────────────────


@pytest.fixture
def solver_2g_p1_n2n():
    """2G solver with non-trivial cross-group P0 AND a Pℓ≥1 channel AND
    a non-zero (n,2n) matrix. Stresses every channel of the residual."""
    # Asymmetric P0 matrix: non-trivial diagonal AND off-diagonal entries.
    p0 = np.array([[0.38, 0.10], [0.05, 0.90]])
    # Non-trivial P1 block — Pℓ≥1 is unconditionally residual.
    p1 = np.array([[0.02, 0.01], [0.00, 0.04]])
    from scipy.sparse import csr_matrix
    mix = make_mixture(
        sig_t=np.array([0.5, 1.0]),
        sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.01, 0.08]),
        nu=np.array([2.5, 2.5]),
        chi=np.array([1.0, 0.0]),
        sig_s=p0,
    )
    # Append a P1 block manually; ORPHEUS's SigS is a list[csr_matrix]
    # indexed by Legendre order.
    mix.SigS = [csr_matrix(p0), csr_matrix(p1)]
    # Inject (n,2n) — non-zero on a cross-group entry only (the brief
    # explicitly notes diagonal sig2 entries are rare but legal).
    mix.Sig2 = csr_matrix(np.array([[0.0, 0.03], [0.01, 0.0]]))

    nx, ny = 3, 2
    mesh = _uniform_2d(nx, ny, 0.4, np.zeros((nx, ny), dtype=int))
    quad = Quadrature.lebedev(order=17)
    return SNSolver(SNMesh(mesh, quad, {0: mix}), scattering_order=1)


class TestFoldablePart:
    """``foldable_part()`` carries ONLY the P0 within-group diagonal."""

    def test_returns_scattering_operator_instance(self, solver_2g_p0):
        """Mechanism criterion 1 — sibling class, not a new class."""
        S = solver_2g_p0.scattering_op
        assert isinstance(S.foldable_part(), ScatteringOperator)

    def test_scattering_order_is_zero(self, solver_2g_p0):
        """Mechanism criterion 2 — no Pℓ structure in foldable."""
        S = solver_2g_p0.scattering_op
        assert S.foldable_part().scattering_order == 0

    def test_Y_is_None(self, solver_2g_p0):
        """Mechanism criterion 3 — no spherical harmonics for ℓ=0."""
        S = solver_2g_p0.scattering_op
        assert S.foldable_part().Y is None

    def test_Y_is_None_even_for_p1_source(self, solver_2g_p1_n2n):
        """Even when S carries P1+ data, the foldable sibling has no Y."""
        S = solver_2g_p1_n2n.scattering_op
        assert S.foldable_part().Y is None

    def test_sig_s_is_diagonal_only(self, solver_2g_p1_n2n):
        """Mechanism criterion 4a — sig_s[mid][0] is diagonal-only."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig_s:
            mat = S_fold.sig_s[mid][0]
            expected = np.diag(np.diag(S.sig_s[mid][0]))
            np.testing.assert_array_equal(mat, expected)
            # Off-diagonal is literally zero, not just small.
            off_diag = mat - np.diag(np.diag(mat))
            assert np.all(off_diag == 0.0)

    def test_sig_s0_matches_sig_s_l0(self, solver_2g_p1_n2n):
        """Mechanism criterion 4b — sig_s0 == sig_s[mid][0]."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig_s:
            np.testing.assert_array_equal(
                S_fold.sig_s0[mid], S_fold.sig_s[mid][0]
            )

    def test_sig_s_has_length_one(self, solver_2g_p1_n2n):
        """Mechanism criterion 4c — no Pℓ≥1 entries in foldable."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig_s:
            assert len(S_fold.sig_s[mid]) == 1

    def test_sig2_is_zero_matrix(self, solver_2g_p1_n2n):
        """Mechanism criterion 4d — (n,2n) belongs to residual unconditionally."""
        S = solver_2g_p1_n2n.scattering_op
        S_fold = S.foldable_part()
        for mid in S.sig2:
            assert S_fold.sig2[mid].shape == S.sig2[mid].shape
            assert S_fold.sig2[mid].dtype == S.sig2[mid].dtype
            np.testing.assert_array_equal(
                S_fold.sig2[mid], np.zeros_like(S.sig2[mid])
            )

    def test_does_not_mutate_parent_sig_s(self, solver_2g_p1_n2n):
        """Anti-rec 4 — split returns new arrays; parent unchanged."""
        S = solver_2g_p1_n2n.scattering_op
        # Snapshot every parent array.
        before = {mid: [m.copy() for m in S.sig_s[mid]] for mid in S.sig_s}
        before_sig2 = {mid: S.sig2[mid].copy() for mid in S.sig2}
        _ = S.foldable_part()
        # Parent unchanged.
        for mid in S.sig_s:
            for l, m in enumerate(S.sig_s[mid]):
                np.testing.assert_array_equal(m, before[mid][l])
        for mid in S.sig2:
            np.testing.assert_array_equal(S.sig2[mid], before_sig2[mid])


class TestResidualPart:
    """``residual_part()`` carries everything but P0 within-group diagonal."""

    def test_returns_scattering_operator_instance(self, solver_2g_p0):
        """Mechanism criterion 1 — sibling class."""
        S = solver_2g_p0.scattering_op
        assert isinstance(S.residual_part(), ScatteringOperator)

    def test_sig_s_l0_diagonal_zeroed(self, solver_2g_p1_n2n):
        """Mechanism criterion 5a — cross-group only on P0."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        for mid in S.sig_s:
            expected = S.sig_s[mid][0] - np.diag(np.diag(S.sig_s[mid][0]))
            np.testing.assert_array_equal(S_res.sig_s[mid][0], expected)
            # The diagonal IS zero, not just close.
            diag = np.diag(S_res.sig_s[mid][0])
            assert np.all(diag == 0.0)

    def test_sig_s0_matches_diagonal_zeroed(self, solver_2g_p1_n2n):
        """Mechanism criterion 5b — sig_s0 alias of sig_s[mid][0]."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        for mid in S.sig_s:
            np.testing.assert_array_equal(
                S_res.sig_s0[mid], S_res.sig_s[mid][0]
            )

    def test_pl_ge_1_carried_verbatim(self, solver_2g_p1_n2n):
        """Mechanism criterion 5c — Pℓ≥1 blocks unchanged."""
        S = solver_2g_p1_n2n.scattering_op
        assert S.scattering_order >= 1, "fixture must carry P1+ data"
        S_res = S.residual_part()
        for mid in S.sig_s:
            for l in range(1, S.scattering_order + 1):
                np.testing.assert_array_equal(
                    S_res.sig_s[mid][l], S.sig_s[mid][l]
                )

    def test_sig2_carried_verbatim(self, solver_2g_p1_n2n):
        """Mechanism criterion 5d — (n,2n) unconditionally residual."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        for mid in S.sig2:
            np.testing.assert_array_equal(S_res.sig2[mid], S.sig2[mid])

    def test_scattering_order_preserved(self, solver_2g_p1_n2n):
        """Mechanism criterion 5e — Pℓ structure preserved."""
        S = solver_2g_p1_n2n.scattering_op
        assert S.residual_part().scattering_order == S.scattering_order

    def test_Y_is_self_Y(self, solver_2g_p1_n2n):
        """Mechanism criterion 5f — precomputed harmonics reusable."""
        S = solver_2g_p1_n2n.scattering_op
        S_res = S.residual_part()
        # Either same object (preferred) or array-equal.
        assert S_res.Y is S.Y or np.array_equal(S_res.Y, S.Y)

    def test_Y_None_for_p0_solver(self, solver_2g_p0):
        """If S has no harmonics (L=0), residual has none either."""
        S = solver_2g_p0.scattering_op
        assert S.residual_part().Y is None

    def test_does_not_mutate_parent_sig_s(self, solver_2g_p1_n2n):
        """Anti-rec 4 — split returns new arrays; parent unchanged."""
        S = solver_2g_p1_n2n.scattering_op
        before = {mid: [m.copy() for m in S.sig_s[mid]] for mid in S.sig_s}
        _ = S.residual_part()
        for mid in S.sig_s:
            for l, m in enumerate(S.sig_s[mid]):
                np.testing.assert_array_equal(m, before[mid][l])


class TestFoldableSigma:
    """``foldable_sigma()`` returns the per-material (ng,) σ_{s,0}^{g→g}."""

    def test_returns_dict_of_ndarrays(self, solver_2g_p1_n2n):
        """Mechanism criterion 6a — dict[int, ndarray]."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        assert isinstance(result, dict)
        for mid, arr in result.items():
            assert isinstance(mid, int)
            assert isinstance(arr, np.ndarray)

    def test_shape_is_ng(self, solver_2g_p1_n2n):
        """Mechanism criterion 6b — each value is (ng,)."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        for arr in result.values():
            assert arr.shape == (S.ng,)

    def test_values_are_diagonal_of_sig_s0(self, solver_2g_p1_n2n):
        """Mechanism criterion 6c — equals np.diag(sig_s[mid][0])."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        for mid, arr in result.items():
            np.testing.assert_array_equal(arr, np.diag(S.sig_s[mid][0]))

    def test_returned_arrays_are_copies(self, solver_2g_p1_n2n):
        """Mutating the returned dict's values must not affect ``self``."""
        S = solver_2g_p1_n2n.scattering_op
        result = S.foldable_sigma()
        # Snapshot parent diagonal.
        before = {mid: np.diag(S.sig_s[mid][0]).copy() for mid in S.sig_s}
        # Mutate the returned arrays.
        for arr in result.values():
            arr[:] = -999.0
        # Parent unchanged.
        for mid in S.sig_s:
            np.testing.assert_array_equal(
                np.diag(S.sig_s[mid][0]), before[mid]
            )


class TestAlgebraicIdentity:
    """The load-bearing contract:
    ``S.apply(ψ) ≈ S.foldable_part().apply(ψ) + S.residual_part().apply(ψ)``
    at ``rtol=1e-14`` (FP-non-associativity precision).

    Covers P0-only, Pℓ≥1, non-zero (n,2n), and cross-group + diagonal
    coupling — the four cases enumerated in the brief's criterion 7.
    """

    def _check_identity(self, op, psi):
        """``psi`` is a typed :class:`AngularFlux`; ``apply`` returns
        :class:`AngularSourceSink`.  Compare via ``.values``."""
        full = op.apply(psi)
        split_sum = op.foldable_part().apply(psi) + op.residual_part().apply(psi)
        np.testing.assert_allclose(full.values, split_sum.values, rtol=1e-14, atol=1e-15)

    def test_identity_p0_only_random_psi(self, solver_2g_p0):
        """Case 1 — scattering_order == 0 only (no Pℓ).

        D-I.2: typed AngularFlux carrier.
        """
        op = solver_2g_p0.scattering_op
        assert op.scattering_order == 0
        N = op.n_ordinates
        np.random.seed(42)
        psi_values = np.random.rand(N, op.ng, *op.spatial_shape) + 0.1
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p0.sn_mesh)
        self._check_identity(op, psi)

    def test_identity_p0_only_uniform_psi(self, solver_2g_p0):
        """Case 1b — uniform ψ probes the diagonal isolation path."""
        op = solver_2g_p0.scattering_op
        N = op.n_ordinates
        psi_values = np.ones((N, op.ng, *op.spatial_shape))
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p0.sn_mesh)
        self._check_identity(op, psi)

    def test_identity_with_pl_ge_1(self, solver_2g_p1_n2n):
        """Case 2 — scattering_order >= 1 (with non-zero P1 block)."""
        op = solver_2g_p1_n2n.scattering_op
        assert op.scattering_order >= 1
        N = op.n_ordinates
        np.random.seed(101)
        psi_values = np.random.rand(N, op.ng, *op.spatial_shape) + 0.1
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p1_n2n.sn_mesh)
        self._check_identity(op, psi)

    def test_identity_with_nonzero_n2n(self, solver_2g_p1_n2n):
        """Case 3 — non-zero (n,2n) coupling."""
        op = solver_2g_p1_n2n.scattering_op
        # Fixture explicitly sets (n,2n) cross-group entries.
        any_nonzero_n2n = any(
            np.any(op.sig2[mid] != 0.0) for mid in op.sig2
        )
        assert any_nonzero_n2n, "fixture must carry non-zero (n,2n)"
        N = op.n_ordinates
        np.random.seed(202)
        psi_values = np.random.rand(N, op.ng, *op.spatial_shape) + 0.1
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p1_n2n.sn_mesh)
        self._check_identity(op, psi)

    def test_identity_multigroup_cross_group_plus_diagonal(self, solver_2g_p1_n2n):
        """Case 4 — non-trivial cross-group P0 + diagonal coupling."""
        op = solver_2g_p1_n2n.scattering_op
        # Fixture's P0 matrix has both diagonal AND off-diagonal entries.
        for mid in op.sig_s:
            p0 = op.sig_s[mid][0]
            diag = np.diag(p0)
            off = p0 - np.diag(diag)
            assert np.any(diag != 0.0)
            assert np.any(off != 0.0)
        N = op.n_ordinates
        np.random.seed(303)
        psi_values = np.random.rand(N, op.ng, *op.spatial_shape) + 0.1
        psi = AngularFlux.from_mesh(psi_values, solver_2g_p1_n2n.sn_mesh)
        self._check_identity(op, psi)

    def test_residual_zero_when_p0_diagonal_only_no_n2n(self):
        """Pure-diagonal P0 with no (n,2n) and no Pℓ≥1 ⇒ residual.apply(ψ)=0
        and full == foldable.apply(ψ) by construction."""
        from scipy.sparse import csr_matrix
        # Strictly diagonal P0, zero (n,2n).
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.0, 0.0]),
            nu=np.array([0.0, 0.0]),
            chi=np.zeros(2),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
            sig_s=np.diag([0.3, 0.8]),
        )
        mix.Sig2 = csr_matrix(np.zeros((2, 2)))
        nx, ny = 2, 2
        mesh = _uniform_2d(nx, ny, 0.5, np.zeros((nx, ny), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver = SNSolver(SNMesh(mesh, quad, {0: mix}))
        op = solver.scattering_op

        N = op.n_ordinates
        np.random.seed(404)
        # D-I.2: typed AngularFlux carrier.
        psi_values = np.random.rand(N, op.ng, nx, ny) + 0.1
        psi = AngularFlux.from_mesh(psi_values, solver.sn_mesh)
        full = op.apply(psi)
        residual_part = op.residual_part().apply(psi)
        np.testing.assert_allclose(residual_part.values, 0.0, atol=1e-15)
        # And full ≡ foldable up to FP-non-associativity.
        foldable_part = op.foldable_part().apply(psi)
        np.testing.assert_allclose(full.values, foldable_part.values, rtol=1e-14, atol=1e-15)


class TestPurity:
    """``foldable_part()`` / ``residual_part()`` are pure functions —
    calling twice returns instances with equal per-material arrays
    (mechanism criterion 8)."""

    def test_foldable_part_pure(self, solver_2g_p1_n2n):
        S = solver_2g_p1_n2n.scattering_op
        a, b = S.foldable_part(), S.foldable_part()
        assert a.scattering_order == b.scattering_order
        assert a.Y is None and b.Y is None
        for mid in S.sig_s:
            np.testing.assert_array_equal(a.sig_s[mid][0], b.sig_s[mid][0])
            np.testing.assert_array_equal(a.sig_s0[mid], b.sig_s0[mid])
            np.testing.assert_array_equal(a.sig2[mid], b.sig2[mid])

    def test_residual_part_pure(self, solver_2g_p1_n2n):
        S = solver_2g_p1_n2n.scattering_op
        a, b = S.residual_part(), S.residual_part()
        assert a.scattering_order == b.scattering_order
        for mid in S.sig_s:
            for l in range(S.scattering_order + 1):
                np.testing.assert_array_equal(
                    a.sig_s[mid][l], b.sig_s[mid][l]
                )
            np.testing.assert_array_equal(a.sig2[mid], b.sig2[mid])

    def test_foldable_sigma_pure(self, solver_2g_p1_n2n):
        S = solver_2g_p1_n2n.scattering_op
        a, b = S.foldable_sigma(), S.foldable_sigma()
        assert set(a.keys()) == set(b.keys())
        for mid in a:
            np.testing.assert_array_equal(a[mid], b[mid])


# ──────────────────────────────────────────────────────────────────────
# is_foldable_into_sigma_r — Phase G Step 3+4.b.i (Issue #196)
# ──────────────────────────────────────────────────────────────────────


class TestIsFoldableIntoSigmaR:
    """``S.is_foldable_into_sigma_r()`` returns True iff S carries only
    diagonal P0 + zero sig2.

    Consumed by substep 3+4.b.ii's ``OperatorSum.solve`` fusion hook
    to detect "this S is the foldable_part — fuse into σ_r and route
    to the within-group sweep". A STRUCTURAL predicate on the
    operator's data, not an identity claim about its action.
    """

    def test_full_scattering_returns_false(self, solver_2g_p1_n2n):
        """Full S with non-zero off-diagonal P0 + non-zero P1 + non-zero
        sig2 → NOT foldable."""
        S = solver_2g_p1_n2n.scattering_op
        # Sanity: the fixture's S has all three non-foldable channels.
        assert S.scattering_order >= 1
        assert S.is_foldable_into_sigma_r() is False

    def test_foldable_part_roundtrip_is_true(self, solver_2g_p1_n2n):
        """``S.foldable_part().is_foldable_into_sigma_r() == True``.

        The load-bearing round-trip: the operator constructed by
        ``foldable_part()`` IS, by definition, the foldable part of
        itself.
        """
        S = solver_2g_p1_n2n.scattering_op
        foldable = S.foldable_part()
        assert foldable.is_foldable_into_sigma_r() is True

    def test_residual_part_returns_false(self, solver_2g_p1_n2n):
        """``S.residual_part().is_foldable_into_sigma_r() == False``.

        The residual carries the cross-group off-diagonal P0
        unconditionally (every multi-group system has at least one
        cross-group entry) — so the diagonal-only check fails.
        """
        S = solver_2g_p1_n2n.scattering_op
        residual = S.residual_part()
        assert residual.is_foldable_into_sigma_r() is False

    def test_p0_only_diagonal_returns_true(self):
        """P0-only ScatteringOperator with diagonal sig_s + zero sig2 →
        True (positive control).

        Build a synthetic ScatteringOperator directly (bypassing
        SNSolver) to isolate the predicate from any fixture setup.
        """
        ng = 2
        p0_diag = np.diag([0.38, 0.90])
        mat_xs = material_xs_from_raw(
            sig_s={0: [p0_diag]},
            sig2={0: np.zeros((ng, ng))},
            cells_by_mat={0: (
                np.array([0, 0, 1, 1]),
                np.array([0, 1, 0, 1]),
            )},
            ng=ng, nx=2, ny=2,
        )
        S = ScatteringOperator(
            mat_xs=mat_xs,
            quadrature=_StubQuad(N=12, weights=np.ones(12) / 12.0),
            scattering_order=0,
        )
        assert S.is_foldable_into_sigma_r() is True

    def test_p0_with_off_diagonal_returns_false(self):
        """scattering_order=0 with non-diagonal P0 → False.

        Off-diagonal P0 is cross-group scattering — couples distinct
        energy groups and cannot collapse into a per-cell scalar.
        """
        ng = 2
        # Non-diagonal P0 — non-zero off-diagonal entry.
        p0 = np.array([[0.38, 0.10], [0.00, 0.90]])
        mat_xs = material_xs_from_raw(
            sig_s={0: [p0]},
            sig2={0: np.zeros((ng, ng))},
            cells_by_mat={0: (
                np.array([0, 0, 1, 1]),
                np.array([0, 1, 0, 1]),
            )},
            ng=ng, nx=2, ny=2,
        )
        S = ScatteringOperator(
            mat_xs=mat_xs,
            quadrature=_StubQuad(N=12, weights=np.ones(12) / 12.0),
            scattering_order=0,
        )
        assert S.is_foldable_into_sigma_r() is False

    def test_p0_diagonal_with_nonzero_sig2_returns_false(self):
        """scattering_order=0 with diagonal P0 BUT non-zero sig2 → False.

        (n,2n) doubling is unconditionally residual: folding into a
        "removal" cross-section is conceptually wrong because (n,2n)
        emits two neutrons per absorption.
        """
        ng = 2
        p0_diag = np.diag([0.38, 0.90])
        sig2 = np.array([[0.0, 0.05], [0.0, 0.0]])
        mat_xs = material_xs_from_raw(
            sig_s={0: [p0_diag]},
            sig2={0: sig2},
            cells_by_mat={0: (
                np.array([0, 0, 1, 1]),
                np.array([0, 1, 0, 1]),
            )},
            ng=ng, nx=2, ny=2,
        )
        S = ScatteringOperator(
            mat_xs=mat_xs,
            quadrature=_StubQuad(N=12, weights=np.ones(12) / 12.0),
            scattering_order=0,
        )
        assert S.is_foldable_into_sigma_r() is False

    def test_scattering_order_ge_1_returns_false_even_with_diagonal_p0(
        self,
    ):
        """scattering_order >= 1 → False even if P0 is diagonal.

        Pℓ ≥ 1 is direction-dependent (Y_ℓ^m(Ω_n)) — unconditionally
        residual. The presence of ANY Pℓ ≥ 1 channel disqualifies the
        operator from foldability.
        """
        ng = 2
        p0_diag = np.diag([0.38, 0.90])
        p1 = np.array([[0.02, 0.00], [0.00, 0.04]])
        mat_xs = material_xs_from_raw(
            sig_s={0: [p0_diag, p1]},
            sig2={0: np.zeros((ng, ng))},
            cells_by_mat={0: (
                np.array([0, 0, 1, 1]),
                np.array([0, 1, 0, 1]),
            )},
            ng=ng, nx=2, ny=2,
        )
        S = ScatteringOperator(
            mat_xs=mat_xs,
            quadrature=_StubQuad(N=12, weights=np.ones(12) / 12.0),
            scattering_order=1,
        )
        assert S.is_foldable_into_sigma_r() is False


# ──────────────────────────────────────────────────────────────────────
# Wave T step T.3 — per-ℓ kernel structure tests (substep T.3b).
#
# T.3 design context (test-architect spec §6 Q6, user resolution):
# the §15.2 form `Σ_ℓ (Σ_{s,ℓ} ⊗ A_ℓ ⊗ G_ℓ)` does NOT satisfy the
# disjoint-axes TensorProductOperator contract because the per-material
# per-ℓ einsum couples group + spatial axes via cells_by_material
# indexing.  Math-honest fallback: kernel is an OperatorSum of
# per-ℓ summands (each a custom LinearOperator, NOT a TP).  T.3 is
# therefore NOT the first SOTP production consumer — T.4 (streaming)
# inherits that role.
# ──────────────────────────────────────────────────────────────────────


class TestAnisoMomentSourcePath:
    """Anisotropic in-scatter via the moment→source map ``R·Λ_{ℓ≥1}``.

    Phase 5a (angular-windowing) retired the per-ℓ
    ``_PerLegendreOrderScattering`` kernel — which recomputed ``M`` for
    every Legendre order — in favour of the shared ``R·Λ`` reconstruction.
    The two aniso paths are:

    * :meth:`ScatteringOperator.build_aniso_source` — the full-angular path,
      ``(1/W)·kernel`` where ``kernel = frame.conjugate(Λ) = R∘Λ∘M``;
    * the windowed moment-iterate ``apply`` arm — whose iterate bulk IS
      ``φ`` (the 2-D Cartesian angular-windowing SI iterate), so ``M`` is
      already done; the P4 carve made it the EXPLICIT typed grid path
      ``Λ : HarmonicMomentFlux → HarmonicMomentSourceSink`` then
      ``frame.reconstruct : HarmonicMomentSourceSink → AngularSourceSink``
      (the role-changing edge materialised as a typed carrier).

    The ndarray ``R∘Λ`` reference
    :meth:`ScatteringOperator._aniso_source_from_moment_values`
    (= ``frame.reconstruct_after(Λ)``) remains the 0-ULP scattering-kernel
    crosscheck's oracle — same ``Λ`` kernel + ``R`` face as both production
    paths. ``build_aniso_source``'s numerical correctness is pinned by the
    pre-T.3 bit-identical snapshot below (a structurally-independent
    reference captured BEFORE the kernel ever existed).  This class adds
    the load-bearing Phase-5a guard: the moment ``apply`` arm reproduces
    the full-angular arm bit-for-bit, plus a sentinel that the windowed arm
    actually executes the typed ``frame.reconstruct``.
    """

    @pytest.fixture
    def op_p1(self, solver_2g_p1_n2n):
        """ScatteringOperator with P1 aniso (asymmetric SigS + n2n)."""
        return solver_2g_p1_n2n.scattering_op

    @pytest.fixture
    def op_p0(self, solver_2g_p0):
        """ScatteringOperator with scattering_order=0 (P0 only)."""
        return solver_2g_p0.scattering_op

    def test_moment_apply_arm_bit_identical_to_angular_arm(
        self, op_p1, solver_2g_p1_n2n,
    ):
        """Phase 5a load-bearing guard: ``S.apply(HarmonicMomentFlux)``
        reproduces ``S.apply(AngularFlux)`` BIT-FOR-BIT when the moments
        are ``φ = Mψ``.

        This is the angular-windowing carve's correctness core —
        windowing the within-group SI iterate from a full per-ordinate
        :class:`AngularFlux` down to :class:`HarmonicMomentFlux` moments
        loses NO scattering-source information (``S`` is a pure function
        of the moments).  The P1 asymmetric-``SigS`` + (n,2n) fixture
        activates the ℓ=0 (iso + cross-group + n2n) AND ℓ≥1 (aniso)
        paths, so a windowing that dropped a moment, swapped an index, or
        drifted a convention would fail here (vv-principles: ≥2 groups,
        anisotropic, asymmetric scatter — Modes 2/6 + the dropped-moment
        trap).
        """
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )

        op = op_p1
        psi = self._reproduce_psi(solver_2g_p1_n2n, seed=7)

        # The full-angular arm S consumes today.
        src_angular = op.apply(psi)

        # The windowed path: project φ = Mψ via the frame analysis face
        # (the SAME frame the operator uses), then the moment arm.
        moments = op.frame.analysis.apply(psi.values)
        phi_field = HarmonicMomentFlux.from_mesh_and_L(
            moments, psi.mesh, op.scattering_order,
        )
        src_moments = op.apply(phi_field)

        # Bit-for-bit (de-risk proven: same M, Y_0^0 = 1 ⇒ ℓ=0 IS the
        # scalar flux, shared R·Λ map for ℓ≥1).
        np.testing.assert_array_equal(src_moments.values, src_angular.values)
        # Non-degeneracy: ℓ≥1 genuinely carries signal (else the moment
        # arm collapses to the P0-only arm and the guard is vacuous).
        assert op.scattering_order >= 1
        assert np.any(moments[1:] != 0.0)

    def test_windowed_arm_executes_typed_role_changing_edge(
        self, op_p1, solver_2g_p1_n2n, monkeypatch,
    ):
        """Mode-11 sentinel (vv-principles): the windowed
        ``apply(HarmonicMomentFlux)`` arm actually executes the EXPLICIT typed
        grid path — ``Λ`` constructs a
        :class:`HarmonicMomentSourceSink` (the role-changing edge) which the
        frame's ``R`` then reconstructs. The Phase-5a value guard above proves
        the NUMBERS are right but cannot tell whether the rewired typed line
        ran or a bypass produced the same value; this counter-spy proves it.

        Spy point: ``HarmonicMomentSourceSink.from_mesh_and_L`` is constructed
        ONLY by ``Λ``'s typed (flux → source) arm. If the windowed apply
        bypassed the typed grid (e.g. fell back to the ndarray
        ``reconstruct_after`` reference) this counter would stay at zero.
        """
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )
        from orpheus.transport.source_sinks import HarmonicMomentSourceSink

        calls = {"n": 0}
        original = HarmonicMomentSourceSink.from_mesh_and_L.__func__

        def counting(cls, *args, **kwargs):
            calls["n"] += 1
            return original(cls, *args, **kwargs)

        monkeypatch.setattr(
            HarmonicMomentSourceSink, "from_mesh_and_L", classmethod(counting),
        )

        op = op_p1
        psi = self._reproduce_psi(solver_2g_p1_n2n, seed=7)
        moments = op.frame.analysis.apply(psi.values)
        phi_field = HarmonicMomentFlux.from_mesh_and_L(
            moments, psi.mesh, op.scattering_order,
        )
        _ = op.apply(phi_field)

        assert op.scattering_order >= 1
        assert calls["n"] >= 1, (
            "windowed apply(HarmonicMomentFlux) did not construct a "
            "HarmonicMomentSourceSink — the explicit typed role-changing "
            "edge (Λ: flux → source) was bypassed."
        )

    def _load_snapshot(self):
        from tests.sn._test_helpers import SN_TESTS_ROOT

        snapshot_path = (
            SN_TESTS_ROOT / "_fixtures" / "wave_t_t3" / "pre_t3_snapshots.npz"
        )
        return np.load(snapshot_path)

    def _reproduce_psi(self, solver, seed: int):
        """Mirror `_capture_pre_t3_snapshots.py::_make_psi`."""
        from orpheus.transport.fields.angular_flux import AngularFlux

        N = solver.quad.N
        ng = solver.ng
        nx, ny = solver.sn_mesh.spatial_shape
        rng = np.random.default_rng(seed)
        psi_values = rng.uniform(0.05, 1.0, size=(N, ng, nx, ny))
        return AngularFlux.from_mesh(psi_values, solver.sn_mesh)

    def _reproduce_phi(self, solver, seed: int):
        """Mirror `_capture_pre_t3_snapshots.py::_make_phi`."""
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        ng = solver.ng
        nx, ny = solver.sn_mesh.spatial_shape
        rng = np.random.default_rng(seed)
        phi_values = rng.uniform(0.05, 1.0, size=(ng, nx, ny))
        return ScalarFlux.from_mesh(phi_values, solver.sn_mesh)

    def test_apply_angular_flux_bit_identical_to_pre_t3_snapshot(
        self, op_p1, solver_2g_p1_n2n,
    ):
        """L1-1 per spec §3 — `apply(AngularFlux)` per-ordinate output
        matches the pre-T.3 captured snapshot within
        `nulp ≤ 4·scattering_order`.

        Post-T.3c the AngularFlux arm inherits the kernel-routed
        numerics via its call to `build_aniso_source`.  P0 + (n,2n)
        contribution is bit-identical (unchanged code path).  The
        per-ℓ aniso contribution may drift by `(L+1) × ULP` per the
        principled-equivalence three-criteria gate; the
        `(iso/sum_w + aniso)` combination at the apply boundary
        preserves the drift bound.
        """
        psi = self._reproduce_psi(solver_2g_p1_n2n, seed=20260530)
        out_post_t3 = op_p1.apply(psi).values
        expected = self._load_snapshot()["p1_apply_angular_flux"]
        nulp_bound = max(4, 4 * op_p1.scattering_order)
        np.testing.assert_array_almost_equal_nulp(
            out_post_t3, expected, nulp=nulp_bound,
        )

    def test_apply_scalar_flux_bit_identical_to_pre_t3_snapshot(
        self, op_p1, solver_2g_p1_n2n,
    ):
        """L1-2 per spec §3 — `apply(ScalarFlux)` iso scalar output
        matches the pre-T.3 captured snapshot **bit-identically**.

        The ScalarFlux arm is P0 + (n,2n) only — it does NOT call
        `build_aniso_source` and therefore NOT the kernel.  T.3
        leaves this path untouched; `np.array_equal` is the
        appropriate gate (no FP reduction reorder).
        """
        phi = self._reproduce_phi(solver_2g_p1_n2n, seed=20260530 + 1)
        out_post_t3 = op_p1.apply(phi).values
        expected = self._load_snapshot()["p1_apply_scalar_flux"]
        np.testing.assert_array_equal(out_post_t3, expected)

    def test_apply_timed_full_field_bit_identical_to_pre_t3_snapshot(
        self, op_p1, solver_2g_p1_n2n,
    ):
        """L1-3 per spec §3 — `apply(TimedFullField)` bulk + boundary
        output matches the pre-T.3 captured snapshot.

        Bulk: nulp ≤ 4·order (inherits from AngularFlux-style
        kernel-routed numerics).  Boundary: bit-identical (the
        implicit-zero AngularBoundaryFlux from Option β3 is unchanged
        across T.3).
        """
        from dataclasses import replace

        psi = self._reproduce_psi(solver_2g_p1_n2n, seed=20260530)
        state = TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=solver_2g_p1_n2n.sn_mesh)
        state = replace(state, interior=replace(state.interior, values=psi.values))

        out_post_t3 = op_p1.apply(state)
        snapshots = self._load_snapshot()

        # Bulk: principled-equivalence relaxation.
        nulp_bound = max(4, 4 * op_p1.scattering_order)
        np.testing.assert_array_almost_equal_nulp(
            out_post_t3.interior.values,
            snapshots["p1_apply_timed_full_field_bulk"],
            nulp=nulp_bound,
        )
        # Boundary: bit-identical (implicit zero, untouched by T.3).
        np.testing.assert_array_equal(
            out_post_t3.boundary.values,
            snapshots["p1_apply_timed_full_field_boundary"],
        )

    def test_build_aniso_source_bit_identical_to_pre_t3_snapshot(
        self, op_p1, solver_2g_p1_n2n,
    ):
        """L1-4 per spec §3 — `build_aniso_source` output matches the
        pre-T.3 captured snapshot within `nulp ≤ 4·scattering_order`.

        Pre-T.3 the body inlined `R(Λ(M(psi))) / sum_w`.  Post-Phase-5a
        the body projects `φ = M(psi)` once, then applies the shared
        moment→source map `_aniso_source_from_moment_values` (= `R·Λ`)
        and the `/ sum_w` boundary normalisation — the SAME composition,
        with the per-ℓ `_PerLegendreOrderScattering` kernel (which rebuilt
        `M`/`R` per ℓ) retired.  The reduction tree may differ from the
        capture at the `Σ_ℓ` outer sum; drift bounded by `(L+1) × ULP`
        per `vv-principles` §"Bit-identity vs principled-equivalence"
        three-criteria gate.

        Snapshot lives at
        `tests/sn/_fixtures/wave_t_t3/pre_t3_snapshots.npz`, captured
        in T.3a (commit `ed05ea3`) on the same fixture
        (`solver_2g_p1_n2n`) with the same seed (20260530).
        """
        from orpheus.transport.fields.angular_flux import AngularFlux

        # Reproduce the fixture's seed + psi construction from
        # `_capture_pre_t3_snapshots.py::_make_psi(solver, seed=20260530)`.
        op = op_p1
        sn_mesh = solver_2g_p1_n2n.sn_mesh
        rng = np.random.default_rng(20260530)
        N = op.n_ordinates
        psi_values = rng.uniform(
            0.05, 1.0, size=(N, op.ng, *op.spatial_shape),
        )
        psi = AngularFlux.from_mesh(psi_values, sn_mesh)

        # Post-T.3c output via the kernel-routed `build_aniso_source`.
        out_post_t3 = op.build_aniso_source(psi).values

        # Pre-T.3 snapshot.
        from tests.sn._test_helpers import SN_TESTS_ROOT

        snapshot_path = (
            SN_TESTS_ROOT / "_fixtures" / "wave_t_t3" / "pre_t3_snapshots.npz"
        )
        snapshots = np.load(snapshot_path)
        expected = snapshots["p1_build_aniso_source"]

        nulp_bound = max(4, 4 * op.scattering_order)
        np.testing.assert_array_almost_equal_nulp(
            out_post_t3, expected, nulp=nulp_bound,
        )

    def test_per_material_einsum_invariance_p1(self, op_p1, solver_2g_p1_n2n):
        """L6-1 per spec §3 — `MaterialXSField.apply_legendre_scattering_moments`
        output is bit-identical at P=1 to the pre-T.3 snapshot.

        T.3 does NOT touch `material_xs_field.py:515-572` — this test
        defends against an unintentional modernisation while in the
        file.  The per-material per-ℓ einsum is the leaf primitive;
        no FP reduction reorder; `np.array_equal` is the appropriate
        gate.
        """
        from orpheus.transport.fields.angular_flux import AngularFlux

        # Reproduce the snapshot script's psi (seed=20260530).
        psi_p1 = self._reproduce_psi(solver_2g_p1_n2n, seed=20260530)
        L = 1
        moments_values = op_p1.frame.analysis.apply(psi_p1.values)

        # apply_legendre_scattering_moments inline (mirror snapshot
        # capture script).  skip_l0=False → full block coverage.
        out = op_p1.mat_xs.apply_legendre_scattering_moments(
            moments_values, L=L, skip_l0=False,
        )
        expected = self._load_snapshot()["p1_apply_legendre_scattering_moments"]
        np.testing.assert_array_equal(out, expected)

    def test_per_material_einsum_invariance_p3(self):
        """L6-2 per spec §3 — same as L6-1 but at P=3, exercising the
        higher-order ℓ loop body.

        Builds an independent P3 solver (mirroring
        `_capture_pre_t3_snapshots.py::build_p3_solver`) to reach the
        captured snapshot.
        """
        from orpheus.transport.fields.angular_flux import AngularFlux
        from scipy.sparse import csr_matrix

        p0 = np.array([[0.38, 0.10], [0.05, 0.90]])
        p1 = np.array([[0.02, 0.01], [0.00, 0.04]])
        p2 = np.array([[0.005, 0.002], [0.000, 0.010]])
        p3 = np.array([[0.001, 0.0005], [0.000, 0.002]])
        mix = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.01, 0.08]),
            nu=np.array([2.5, 2.5]),
            chi=np.array([1.0, 0.0]),
            sig_s=p0,
        )
        mix.SigS = [csr_matrix(p0), csr_matrix(p1), csr_matrix(p2), csr_matrix(p3)]
        mix.Sig2 = csr_matrix(np.array([[0.0, 0.03], [0.01, 0.0]]))

        nx, ny = 3, 2
        mesh = _uniform_2d(nx, ny, 0.4, np.zeros((nx, ny), dtype=int))
        quad = Quadrature.lebedev(order=17)
        solver_p3 = SNSolver(SNMesh(mesh, quad, {0: mix}), scattering_order=3)
        op_p3 = solver_p3.scattering_op

        rng = np.random.default_rng(20260530 + 2)
        psi_p3 = AngularFlux.from_mesh(
            rng.uniform(0.05, 1.0, size=(quad.N, 2, nx, ny)), solver_p3.sn_mesh,
        )
        L = 3
        moments_values = op_p3.frame.analysis.apply(psi_p3.values)
        out = op_p3.mat_xs.apply_legendre_scattering_moments(
            moments_values, L=L, skip_l0=False,
        )
        expected = self._load_snapshot()["p3_apply_legendre_scattering_moments"]
        np.testing.assert_array_equal(out, expected)

"""Foundation tests for the AngularClosureBase ABC + concrete strategies.

These tests pin the **software contract** of the angular-closure
strategy family.  They are foundation-tagged because the claims are
about the strategy API (ABC conformance, registry self-registration,
hand-calc closure algebra, α-recursion identities) rather than
transport-equation identities — those are verified transitively via the
operator-level tests at
:file:`tests/sn/test_streaming_operator.py` /
:file:`tests/sn/test_streaming_operator_decomposition.py` and the
curvilinear MMS suite.

PR-TYPED-6c Step 7 (2026-05-18) retired
``LegacyTauSymmetricInterpolation`` (pre-Phase-B inlined form) and
``BaileyFlatFluxRedist`` (Phase B ablation strategy for back-bisection
against the legacy form).  After PR-TYPED-6.5's default switch to
``MorelMontryAngularSweep`` for curvilinear + ``IdentityAngularClosure``
for Cartesian, neither legacy strategy had any production consumer;
the unified matvec body's typed accessor (``closure.cell_contribution``)
makes them unimplementable as full Protocol members without compatibility
shims.  This file's coverage now focuses on the two surviving strategies.

Issue #248 (2026-06-18) retired the orphaned ``@runtime_checkable
PoleAngularClosure`` **Protocol** (Issue #236 Phase 2 B2 retyped every
production consumer onto the ``AngularClosureBase`` **ABC** and made
the strategy methods abstract on it — the Protocol was left orphaned and
divergent).  Conformance is now ABC inheritance, and the hand-calc
redistribution algebra (formerly pinned through the dead legacy
``__call__`` bundle interface) is re-pinned through the LIVE production
surface: the half-angle ψ-thread via the module-level
:func:`~orpheus.sn.angular.closure.compute_psi_half_per_level`
(the same recurrence kernel ``_psi_half_grid_single_level`` the matvec's
:meth:`~MorelMontryAngularSweep.precompute_psi_state` consumes) composed
with the geometry redistribution fold ``(ΔA/w)/V·(α_{m+1/2}ψ_{m+1/2} −
α_{m-1/2}ψ_{m-1/2})`` the test owns explicitly.

C5 (2026-07-03) retired the unbound ``MorelMontryAngularSweep()`` legacy
mode — construction tests bind to the tiny-sphere SNMesh helper (the
family's ``cls(sn_mesh)`` contract), and the hand-calc algebra tests
call the pure module-level surface (no instance at all).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.angular.closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    AngularClosureBase,
)
from tests.sn._test_helpers import (
    make_tiny_spherical_sn_mesh,
    redistribution_via_live_path,
)

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Protocol conformance
# ═══════════════════════════════════════════════════════════════════════


def _mm_from(sn_mesh) -> MorelMontryAngularSweep:
    """Build the closure from a mesh's two TENSOR FACTORS.

    The family's contract is ``cls(angular, pairing)`` (the un-weld arc's
    Phase B) — the mesh is no longer an operand.  These protocol/repr rows
    only need *a* constructed closure, so they spell the two-factor call
    once here rather than at each site.
    """
    return MorelMontryAngularSweep(
        sn_mesh.reduced.angular, sn_mesh.reduced.redistribution_pairing,
    )

class TestProtocolConformance:
    """The :class:`AngularClosureBase` ABC contract is honoured by the
    surviving strategies.  Issue #248 retired the orphaned
    ``@runtime_checkable PoleAngularClosure`` Protocol; conformance is now
    ABC inheritance (the production matvec/sweep type against the ABC, and
    the strategy methods ``precompute_psi_state`` / ``cell_contribution`` /
    ``angular_adjoint`` are abstract on it — Issue #236 Phase 2 B2)."""

    def test_morel_montry_satisfies_protocol(self) -> None:
        """``MorelMontryAngularSweep`` is a ``AngularClosureBase`` subtype."""
        closure = _mm_from(make_tiny_spherical_sn_mesh())
        assert isinstance(closure, AngularClosureBase)

    def test_is_linear_class_attr_advertised(self) -> None:
        """``MorelMontryAngularSweep`` advertises ``is_linear = True``."""
        assert MorelMontryAngularSweep.is_linear is True

    def test_morel_montry_repr(self) -> None:
        # The "()" repr is contractual — the bound mesh is an
        # implementation detail not surfaced in the repr.
        closure = _mm_from(make_tiny_spherical_sn_mesh())
        assert repr(closure) == "MorelMontryAngularSweep()"


# ═══════════════════════════════════════════════════════════════════════
# Registry / self-registration
# ═══════════════════════════════════════════════════════════════════════


class TestRegistry:
    """RegistryMixin self-registration for the AngularClosure family.

    Mirrors the BoundaryFaceFlux / DiscretizationScheme registry contracts —
    each strategy self-registers at import time and is name-keyed
    addressable via ``AngularClosureBase.create(...)``.
    """

    def test_morel_montry_registered(self) -> None:
        assert "morel_montry_angular_sweep" in AngularClosureBase.registry
        assert (
            AngularClosureBase.registry["morel_montry_angular_sweep"]
            is MorelMontryAngularSweep
        )

    def test_create_morel_montry_returns_instance(self) -> None:
        # ``create`` forwards kwargs to the concrete ``__init__`` — the
        # family's ``cls(angular, pairing)`` two-tensor-factor contract rides
        # through the registry (C5: no unbound construction; the un-weld
        # arc's Phase B replaced the mesh operand with the two factors).
        sn_mesh = make_tiny_spherical_sn_mesh()
        instance = AngularClosureBase.create(
            "morel_montry_angular_sweep",
            angular=sn_mesh.reduced.angular,
            pairing=sn_mesh.reduced.redistribution_pairing,
        )
        assert isinstance(instance, MorelMontryAngularSweep)

    def test_create_unknown_key_raises(self) -> None:
        with pytest.raises(
            KeyError, match="unknown AngularClosureBase key",
        ):
            AngularClosureBase.create("unknown_strategy")


# ═══════════════════════════════════════════════════════════════════════
# α-recursion identity (Hébert Eq. 3.423-3.424 in ORPHEUS normalization)
# ═══════════════════════════════════════════════════════════════════════


class TestAlphaRecursionIdentities:
    """Verify the α-dome recursion identities Phase B's pole closure
    consumes match Hébert §3.9.4 (in the ORPHEUS factor-of-2-absorbed
    normalization).

    ORPHEUS recurrence: ``α[n+1] = α[n] - w[n] * μ[n]`` with α[0] = 0.
    Hébert recurrence: ``α^H_{n+1/2} = α^H_{n-1/2} - 2·𝒲_n·μ_n`` with
    α^H_{1/2} = 0.  ORPHEUS α corresponds to Hébert α/2 (with the
    factor of 2 absorbed into the redistribution divisor); the
    recurrences are equivalent under that re-scaling.
    """

    def test_alpha_endpoints_zero_gauss_legendre(self) -> None:
        """α_{1/2} = α_{N+1/2} = 0 by Hébert Eq. 3.423 + GL antisymmetry."""
        from orpheus.numerics.quadrature import Quadrature
        for N in (4, 8, 16):
            quad = Quadrature.gauss_legendre(n_ordinates=N)
            mu = quad.mu_x
            w = quad.weights
            alpha = np.zeros(N + 1)
            for n in range(N):
                alpha[n + 1] = alpha[n] - w[n] * mu[n]
            assert alpha[0] == 0.0, (
                f"α_{{1/2}} must be zero (Hébert Eq. 3.423); got {alpha[0]}"
            )
            assert abs(alpha[N]) < 1e-13, (
                f"α_{{N+1/2}} must be zero by GL antisymmetry; "
                f"got {alpha[N]:.3e} at N={N}"
            )

    def test_alpha_recurrence_signed_step(self) -> None:
        """``α[n+1] - α[n] = -w[n]·μ[n]`` (the recurrence rule)."""
        from orpheus.numerics.quadrature import Quadrature
        N = 8
        quad = Quadrature.gauss_legendre(n_ordinates=N)
        mu = quad.mu_x
        w = quad.weights
        alpha = np.zeros(N + 1)
        for n in range(N):
            alpha[n + 1] = alpha[n] - w[n] * mu[n]
        for n in range(N):
            np.testing.assert_allclose(
                alpha[n + 1] - alpha[n], -w[n] * mu[n], rtol=1e-13,
                err_msg=f"α-recurrence step n={n} violated",
            )

    def test_alpha_dome_non_negative(self) -> None:
        """α is non-negative across the GL quadrature (the 'dome' property)."""
        from orpheus.numerics.quadrature import Quadrature
        for N in (4, 8, 16):
            quad = Quadrature.gauss_legendre(n_ordinates=N)
            mu = quad.mu_x
            w = quad.weights
            alpha = np.zeros(N + 1)
            for n in range(N):
                alpha[n + 1] = alpha[n] - w[n] * mu[n]
            assert np.all(alpha >= -1e-13), (
                f"α-dome non-negativity violated at N={N}; got {alpha}"
            )


# ═══════════════════════════════════════════════════════════════════════
# Hand-calc fixtures: 2-ordinate sphere
# ═══════════════════════════════════════════════════════════════════════
#
# Smallest non-trivial sphere case for verbatim algebra tests.  N=2 GL
# gives μ = ±1/√3, w = 1; α[0] = 0, α[1] = -w·μ_0 = 1/√3, α[2] = α[1] -
# w·μ_1 = 0 (closes by GL antisymmetry).  τ_mm clamped to 1/2 (the M-M
# clamp's lower bound).


@pytest.fixture
def two_ordinate_sphere_fixture():
    """Tiny 2-ordinate sphere fixture.

    Returns ``(psi_cells, alpha, dAw, tau, V)``:

    * ``psi_cells`` shape ``(1, 2, 1)`` — single cell, 2 ordinates, 1 group.
    * ``alpha`` shape ``(3,)`` — Hébert/GL recurrence values.
    * ``dAw`` shape ``(1, 2)`` — geometry factor.
    * ``tau`` shape ``(2,)`` — M-M clamp.
    * ``V`` shape ``(1,)`` — volume.
    """
    # GL N=2: μ = ±1/√3, w = 1.0, α[0]=0, α[1]=1/√3, α[2]=0.
    inv_sqrt3 = 1.0 / np.sqrt(3.0)
    alpha = np.array([0.0, inv_sqrt3, 0.0])
    # Choose ψ values that vary in angle: ψ_0 = 1, ψ_1 = 3.
    psi_cells = np.array([[[1.0], [3.0]]])  # (ng=1, N=2, nx=1)
    dAw = np.array([[1.0, 1.0]])             # (nx=1, N=2)
    tau = np.array([0.5, 0.5])               # (N=2,) clamped at lower bound
    V = np.array([1.0])                      # (nx=1,)
    return psi_cells, alpha, dAw, tau, V


class TestMorelMontryHandCalc:
    """Verbatim Hébert §3.9.4 algebra on a 2-ordinate sphere fixture.

    Issue #248 — re-pinned onto the LIVE production path (the half-angle
    ψ-thread via :func:`~orpheus.sn.angular.closure.compute_psi_half_per_level`
    plus the explicit geometry redistribution fold) after the dead legacy
    ``__call__`` bundle interface was retired.  The expected redistribution
    values (sphere ``R_0 = 2/√3``, ``R_1 = -2/√3``) are unchanged; only the
    surface that produces them moved from the orphaned bundle path to the
    matvec's own recurrence kernel.
    """

    def test_morel_montry_2ord_recurrence(
        self, two_ordinate_sphere_fixture,
    ) -> None:
        r"""Per-cell M-M weighted DD recurrence on N=2 sphere.

        Hand-calc with τ = 1/2 (pure DD limit):

        * Seed: ψ_face_{1/2} = 0.
        * Ordinate 0 (μ_0 = -1/√3, ψ_0 = 1):
            ψ_face_{3/2} = (1 - 0.5·0)/0.5 = 2.
            R_0 = (ΔA/w)·(α_{3/2}·2 - α_{1/2}·0) / V
                = 1·(1/√3·2 - 0·0) / 1 = 2/√3.
        * Ordinate 1 (μ_1 = +1/√3, ψ_1 = 3):
            ψ_face_{5/2} = (3 - 0.5·2)/0.5 = 4.
            R_1 = (ΔA/w)·(α_{5/2}·4 - α_{3/2}·2) / V
                = 1·(0·4 - 1/√3·2) / 1 = -2/√3.

        The ψ-thread (``ψ_face_{3/2} = 2``, ``ψ_face_{5/2} = 4``) is pinned
        by the live
        :func:`~orpheus.sn.angular.closure.compute_psi_half_per_level`
        recurrence; the α-weighted geometry fold gives ``R_0``/``R_1``.
        """
        psi, alpha, dAw, tau, V = two_ordinate_sphere_fixture
        result = redistribution_via_live_path(psi, alpha, dAw, tau, V)
        inv_sqrt3 = 1.0 / np.sqrt(3.0)
        expected = np.array([[[2.0 * inv_sqrt3], [-2.0 * inv_sqrt3]]])
        np.testing.assert_allclose(result, expected, rtol=1e-14)

    def test_morel_montry_recurrence_seed_alpha_half_zero_unused(
        self, two_ordinate_sphere_fixture,
    ) -> None:
        """At ordinate 0 the term ``α_{1/2}·ψ_{1/2,i}`` should vanish
        (α_{1/2} = 0), so the seed value is irrelevant for the n=0
        contribution.  Verify algebraically: α_{1/2} = alpha[0] = 0 ⇒ the
        term is zero regardless of seed."""
        psi, alpha, dAw, tau, V = two_ordinate_sphere_fixture
        # The closure's recurrence seeds at ψ_{1/2}=0 (Carlson convention);
        # the seed enters R_0 only via α_{1/2}·ψ_{1/2}, and α_{1/2}=0.
        assert alpha[0] == 0.0
        # Running the live path gives the expected n=0 redistribution
        # (which only depends on α_{3/2}·ψ_face_{3/2}).
        result = redistribution_via_live_path(psi, alpha, dAw, tau, V)
        inv_sqrt3 = 1.0 / np.sqrt(3.0)
        # Expected n=0: α[1]·ψ_face_{3/2} = (1/√3)·2 = 2/√3.
        np.testing.assert_allclose(
            result[0, 0, 0], 2.0 * inv_sqrt3, rtol=1e-14,
        )


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical multi-level sanity check
# ═══════════════════════════════════════════════════════════════════════


class TestCylindricalLevelDispatch:
    """The cylindrical multi-level redistribution: each μ-level runs the
    per-level recurrence independently (its own α-dome / ΔA/w / τ) and
    scatters back to global ordinate slots.

    Issue #248 — the dead legacy ``__call__(... level_indices=...)`` bundle
    is gone; the per-level dispatch is reconstructed here through the LIVE
    :func:`~orpheus.sn.angular.closure.compute_psi_half_per_level` recurrence
    (one call per level on the level's ordinate slice) plus the explicit
    geometry fold.  This mirrors how production drives cylindrical: the
    matvec's :meth:`~MorelMontryAngularSweep.precompute_psi_state` loops
    :attr:`level_indices` calling :meth:`_psi_half_grid_for_level` — the
    same single-level recurrence — per μ-level.  The hand-calc
    redistribution values (``R = ±1.0`` on level 0, ``R = ±6.3`` on
    level 1) are unchanged.
    """

    def test_cylindrical_dispatch_two_levels(self) -> None:
        """Two-level cylindrical fixture; the per-level recurrence runs
        independently and scatters back to global ordinate slots."""
        # 4 global ordinates split into 2 levels of 2 each.
        ng, nx = 1, 1
        N_global = 4
        psi_global = np.array([[[1.0], [2.0], [3.0], [5.0]]])  # (1, 4, 1)
        # Level 0 carries ordinates 0, 1; Level 1 carries 2, 3.
        level_indices = [np.array([0, 1]), np.array([2, 3])]
        alpha_per_level = [
            np.array([0.0, 0.5, 0.0]),  # M=2 dome on level 0
            np.array([0.0, 0.7, 0.0]),  # M=2 dome on level 1
        ]
        dAw_per_level = [
            np.array([[1.0, 1.0]]),
            np.array([[1.5, 1.5]]),
        ]
        tau_per_level = [
            np.array([0.5, 0.5]),
            np.array([0.5, 0.5]),
        ]
        V = np.array([1.0])

        result = np.zeros((ng, N_global, nx))
        for p, level_idx in enumerate(level_indices):
            psi_level = psi_global[:, level_idx, :]  # (ng, M_p, nx)
            result[:, level_idx, :] = redistribution_via_live_path(
                psi_level,
                alpha_per_level[p], dAw_per_level[p],
                tau_per_level[p], V,
            )
        assert result.shape == (ng, N_global, nx)

        # Hand-calc level 0 (ψ_0=1, ψ_1=2, α_dome=[0, 0.5, 0], τ=0.5):
        # ψ_face_{1/2} = 0
        # ψ_face_{3/2} = (1 - 0)/0.5 = 2
        # R_0 = 1·(0.5·2 - 0·0)/1 = 1.0
        # ψ_face_{5/2} = (2 - 0.5·2)/0.5 = 2
        # R_1 = 1·(0·2 - 0.5·2)/1 = -1.0
        np.testing.assert_allclose(result[0, 0, 0], 1.0, rtol=1e-14)
        np.testing.assert_allclose(result[0, 1, 0], -1.0, rtol=1e-14)

        # Hand-calc level 1 (ψ_2=3, ψ_3=5, α_dome=[0, 0.7, 0], dAw=1.5,
        # τ=0.5):
        # ψ_face_{1/2} = 0
        # ψ_face_{3/2} = (3 - 0)/0.5 = 6
        # R_2 = 1.5·(0.7·6 - 0·0)/1 = 6.3
        # ψ_face_{5/2} = (5 - 0.5·6)/0.5 = 4
        # R_3 = 1.5·(0·4 - 0.7·6)/1 = -6.3
        np.testing.assert_allclose(result[0, 2, 0], 6.3, rtol=1e-14)
        np.testing.assert_allclose(result[0, 3, 0], -6.3, rtol=1e-14)


# ═══════════════════════════════════════════════════════════════════════
# Linearity (a precondition for use as a linear-operator component)
# ═══════════════════════════════════════════════════════════════════════


def test_strategy_is_linear_in_psi():
    """The M-M redistribution is linear in ``psi_cells`` — the algebra
    behind the advertised ``MorelMontryAngularSweep.is_linear = True``.
    Verify ``f(αψ + βφ) = α·f(ψ) + β·f(φ)``.

    Issue #248 — driven through the LIVE redistribution path (the
    module-level
    :func:`~orpheus.sn.angular.closure.compute_psi_half_per_level`
    recurrence composed with the geometry fold) instead of the retired
    ``__call__`` bundle.  Linearity is the load-bearing property for the
    matvec's consumption: both the recurrence and the α-weighted fold are
    linear in ``psi``, so the redistribution is too."""
    rng = np.random.default_rng(seed=42)
    ng, N, nx = 2, 4, 5

    psi1 = rng.standard_normal((ng, N, nx))
    psi2 = rng.standard_normal((ng, N, nx))
    # Random but valid α-dome (non-negative, endpoints zero).
    alpha = np.zeros(N + 1)
    for n in range(N):
        alpha[n + 1] = alpha[n] + abs(rng.standard_normal())
    alpha -= alpha[N] * np.linspace(0, 1, N + 1)  # rebalance to close at 0
    # Make sure endpoints are zero (the dome property).
    alpha[0] = 0.0
    alpha[N] = 0.0
    # Ensure non-negative dome.
    alpha = np.abs(alpha)
    dAw = np.abs(rng.standard_normal((nx, N))) + 0.1
    tau = np.full(N, 0.5)
    V = np.abs(rng.standard_normal(nx)) + 0.1

    a, b = 1.7, -2.3
    lhs = redistribution_via_live_path(a * psi1 + b * psi2, alpha, dAw, tau, V)
    rhs = (
        a * redistribution_via_live_path(psi1, alpha, dAw, tau, V)
        + b * redistribution_via_live_path(psi2, alpha, dAw, tau, V)
    )
    np.testing.assert_allclose(
        lhs, rhs, rtol=1e-13, atol=1e-14,
        err_msg="M-M redistribution non-linearity detected.",
    )


class TestPairingContract:
    r"""The closure takes the two TENSOR FACTORS, and refuses a pairing it
    cannot yet honour.

    The redistribution operator factors as
    :math:`R_{\rm spatial} \otimes A_{\rm angular}`, so a member of this
    family is constructed from exactly those two objects
    (``cls(angular, pairing)``, the un-weld arc's Phase B).  The pairing's axes
    ``(nx, n_mom, n_thread)`` are real — ``n_mom`` is the spatial scheme's
    moment count, ``n_thread`` is what the angular device propagates — but
    the CELL SOLVE for ``n_mom > 1`` does not exist yet (Issue #158's
    curvilinear linear-discontinuous arm), so the constructor must refuse
    rather than silently divide by a scalar.

    ⭐ **These rows are why the guard is not unfalsifiable.** No shipped
    scheme produces a multi-moment pairing, but the pairing is a plain array, so
    a hand-built one is a genuine witness — not a mutation of the code
    under test (``vv-principles`` #17's granularity trap: each refusal arm
    is exercised separately, and each reddens on its own).

    Both refused shapes are the ones the literature actually proposes:
    ``(nx, 2, 2)`` is Adams--Martin 1992 App. A (the angular closure
    applied per spatial moment), ``(nx, 2, 1)`` is ONETRAN (Hill 1975
    Eq. 32, the angular index closed on the cell AVERAGE only).
    """

    @staticmethod
    def _angular():
        sn = make_tiny_spherical_sn_mesh()
        return sn.reduced.angular, sn.reduced.redistribution_pairing

    @pytest.mark.foundation
    def test_the_shipped_gram_is_single_moment_and_admits(self) -> None:
        """The positive leg — without it the refusals below prove nothing."""
        angular, pairing = self._angular()
        if pairing.shape[1:] != (1, 1):
            pytest.fail(
                f"the shipped pairing should be single-moment; got {pairing.shape}"
            )
        closure = MorelMontryAngularSweep(angular, pairing)
        if closure is None:  # pragma: no cover — construction must succeed
            pytest.fail("the single-moment pairing must be admitted")

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mom,n_thread,family", [
        (2, 2, "Adams-Martin: closed per spatial moment"),
        (2, 1, "ONETRAN: closed on the cell average only"),
    ])
    def test_multi_moment_gram_refuses_naming_the_missing_solve(
        self, n_mom: int, n_thread: int, family: str,
    ) -> None:
        angular, pairing = self._angular()
        nx = pairing.shape[0]
        with pytest.raises(NotImplementedError, match="158"):
            MorelMontryAngularSweep(angular, np.zeros((nx, n_mom, n_thread)))

    @pytest.mark.foundation
    def test_a_gram_without_the_moment_axes_refuses(self) -> None:
        """The pre-Phase-B spelling (a bare ``(nx,)`` ΔA) is not a pairing."""
        angular, pairing = self._angular()
        with pytest.raises(ValueError, match="n_mom"):
            MorelMontryAngularSweep(angular, np.zeros(pairing.shape[0]))

    @pytest.mark.foundation
    def test_the_identity_member_takes_the_same_two_operands(self) -> None:
        """The dispatch site constructs whichever class was selected without
        knowing which, so the contract must be uniform — and Cartesian's
        factors are both the NEUTRAL element."""
        angular, pairing = self._angular()
        closure = IdentityAngularClosure(angular, pairing)
        if repr(closure) != "IdentityAngularClosure()":
            pytest.fail(
                "the identity member's repr must not carry a mesh — it no "
                f"longer holds one; got {repr(closure)}"
            )


# ═══════════════════════════════════════════════════════════════════════
# P4.9a — the minted scan constants + the two-representation weld
# ═══════════════════════════════════════════════════════════════════════


class TestMintedScanConstants:
    """The closure MINTS its scan-normal representation; the two forms weld.

    P4.9a (plan phase; ``scratch/p4_9a_verification_plan.md`` §5): the march
    has ONE owner spelling (:func:`march_psi_half_step`) and one minted
    scan-normal form (``tau_inv`` / ``march_a_in_coeff``).  They are
    algebraically identical and NOT bitwise equal — [M] 2026-08-28, real
    fp(4, 6) τ: the STABLE discriminator is absolute, max |Δ| = 1.776e-15
    exact across seeds; bit-equal-fraction and max-ULP are draw-dependent
    (46–51 %, 1e2–1e5 ULP over 200 seeds) — so this class WELDS them
    within an ABSOLUTE band instead of pretending they are one (a nulp
    band on this pair certifies a seed, not the weld: near-zero outputs
    blow up relative ULP without bound).

    In-family value comparisons (``advance_psi_half`` vs the batch kernel)
    are deliberately ABSENT: after the single-sourcing both route through
    one body, so such a gate would be ``f(x) == f(x)``.  The external
    anchors are the CYL_DEG snapshots and the M1 owner-mutation arm.
    """

    @staticmethod
    def _mm_cylinder():
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from tests.sn._test_helpers import placeholder_materials

        mesh = Mesh1D(
            edges=np.linspace(0.01, 2.0, 9),
            mat_ids=np.zeros(8, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.folded_product(n_mu=4, n_phi=6)
        sn = SNMesh(mesh, quad, placeholder_materials(ng=2))
        closure = sn.angular_closure
        assert isinstance(closure, MorelMontryAngularSweep)
        return closure

    def test_minted_constants_pin_their_spelling(self) -> None:
        """[foundation] ``1.0/τ`` and ``(1−τ)/τ`` — bit-exact, spelled THUS.

        Positive leg: ``array_equal`` against the defining expressions on
        the closure's own τ (pins the SPELLING, not an independent value —
        the M7 mutation arm is this gate's teeth).  Discrimination leg:
        the algebraically-equal respelling ``tau_inv − 1.0`` must differ
        ([M] 1–2 ULP on all four measured configs) — proving the pin can
        move under the realistic defect.
        """
        closure = self._mm_cylinder()
        tau = closure.tau_per_ordinate
        assert np.array_equal(closure.tau_inv_per_ordinate, 1.0 / tau)
        assert np.array_equal(
            closure.march_a_in_coeff_per_ordinate, (1.0 - tau) / tau,
        )
        respelled = closure.tau_inv_per_ordinate - 1.0
        assert not np.array_equal(
            closure.march_a_in_coeff_per_ordinate, respelled,
        ), "the discrimination leg lost its subject: (1-τ)/τ == τ⁻¹−1 bitwise"

    def test_march_step_and_scan_normal_form_weld_within_band(self) -> None:
        """[foundation] The two representations agree within an ABSOLUTE
        band, and are NOT bit-equal — the honest weld.

        [M] the band: max |Δ| = 1.776e-15 on real fp(4, 6) τ, EXACT
        across 200 seeds (the archivist's reproduction, 2026-08-28);
        4e-15 gives ~2.2× headroom.  Deliberately NOT a nulp band —
        max-ULP on this pair ranges 1e2–1e5 with the draw (near-zero
        outputs), so a nulp assertion would pin a seed, not the weld.
        """
        closure = self._mm_cylinder()
        tau = closure.tau_per_ordinate
        tau_inv = closure.tau_inv_per_ordinate
        a_in = closure.march_a_in_coeff_per_ordinate
        rng = np.random.default_rng(20260828)
        psi_avg = rng.uniform(0.1, 5.0, size=(64, tau.size))
        psi_in = rng.uniform(0.1, 5.0, size=(64, tau.size))
        form_a = np.empty_like(psi_avg)
        for n in range(tau.size):
            form_a[:, n] = closure.advance_psi_half(
                psi_avg[:, n], psi_in[:, n], ordinate=n,
            )
        form_b = tau_inv * psi_avg - a_in * psi_in
        max_abs_delta = float(np.max(np.abs(form_a - form_b)))
        assert max_abs_delta <= 4.0e-15, (
            f"the two M-M representations drifted past the weld band: "
            f"max |Δ| = {max_abs_delta:.3e} > 4.0e-15"
        )
        assert not np.array_equal(form_a, form_b), (
            "the weld's premise failed: the two forms read bit-equal on "
            "this τ — the band assertion is then vacuous, re-derive it"
        )

    def test_identity_closure_constants_are_neutral(self) -> None:
        """[foundation] Identity: τ⁻¹ ≡ 1, coeff ≡ 0, advance ≡ ψ̄ exactly."""
        sn = make_tiny_spherical_sn_mesh()
        closure = IdentityAngularClosure(
            sn.reduced.angular, sn.reduced.redistribution_pairing,
        )
        n = closure.tau_per_ordinate.size
        assert np.array_equal(closure.tau_inv_per_ordinate, np.ones(n))
        assert np.array_equal(
            closure.march_a_in_coeff_per_ordinate, np.zeros(n),
        )
        psi_avg = np.array([1.3, -0.7, 2.9])
        psi_in = np.array([9.9, 9.9, 9.9])
        out = closure.advance_psi_half(psi_avg, psi_in, ordinate=n - 1)
        assert np.array_equal(out, psi_avg)

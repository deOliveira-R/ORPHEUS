"""Foundation tests for the PoleAngularClosureBase ABC + concrete strategies.

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
production consumer onto the ``PoleAngularClosureBase`` **ABC** and made
the strategy methods abstract on it — the Protocol was left orphaned and
divergent).  Conformance is now ABC inheritance, and the hand-calc
redistribution algebra (formerly pinned through the dead legacy
``__call__`` bundle interface) is re-pinned through the LIVE production
surface: the half-angle ψ-thread via the module-level
:func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
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

from orpheus.sn.sweep.pole_angular_closure import (
    MorelMontryAngularSweep,
    PoleAngularClosureBase,
)
from tests.sn._test_helpers import (
    make_tiny_spherical_sn_mesh,
    redistribution_via_live_path,
)

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# Protocol conformance
# ═══════════════════════════════════════════════════════════════════════


class TestProtocolConformance:
    """The :class:`PoleAngularClosureBase` ABC contract is honoured by the
    surviving strategies.  Issue #248 retired the orphaned
    ``@runtime_checkable PoleAngularClosure`` Protocol; conformance is now
    ABC inheritance (the production matvec/sweep type against the ABC, and
    the strategy methods ``precompute_psi_state`` / ``cell_contribution`` /
    ``angular_adjoint`` are abstract on it — Issue #236 Phase 2 B2)."""

    def test_morel_montry_satisfies_protocol(self) -> None:
        """``MorelMontryAngularSweep`` is a ``PoleAngularClosureBase`` subtype."""
        closure = MorelMontryAngularSweep(make_tiny_spherical_sn_mesh())
        assert isinstance(closure, PoleAngularClosureBase)

    def test_is_linear_class_attr_advertised(self) -> None:
        """``MorelMontryAngularSweep`` advertises ``is_linear = True``."""
        assert MorelMontryAngularSweep.is_linear is True

    def test_morel_montry_repr(self) -> None:
        # The "()" repr is contractual — the bound mesh is an
        # implementation detail not surfaced in the repr.
        closure = MorelMontryAngularSweep(make_tiny_spherical_sn_mesh())
        assert repr(closure) == "MorelMontryAngularSweep()"


# ═══════════════════════════════════════════════════════════════════════
# Registry / self-registration
# ═══════════════════════════════════════════════════════════════════════


class TestRegistry:
    """RegistryMixin self-registration for the PoleAngularClosure family.

    Mirrors the BoundaryFaceFlux / DiscretizationScheme registry contracts —
    each strategy self-registers at import time and is name-keyed
    addressable via ``PoleAngularClosureBase.create(...)``.
    """

    def test_morel_montry_registered(self) -> None:
        assert "morel_montry_angular_sweep" in PoleAngularClosureBase.registry
        assert (
            PoleAngularClosureBase.registry["morel_montry_angular_sweep"]
            is MorelMontryAngularSweep
        )

    def test_create_morel_montry_returns_instance(self) -> None:
        # ``create`` forwards kwargs to the concrete ``__init__`` — the
        # family's ``cls(sn_mesh)`` construction contract rides through
        # the registry (C5: no unbound construction).
        instance = PoleAngularClosureBase.create(
            "morel_montry_angular_sweep",
            sn_mesh=make_tiny_spherical_sn_mesh(),
        )
        assert isinstance(instance, MorelMontryAngularSweep)

    def test_create_unknown_key_raises(self) -> None:
        with pytest.raises(
            KeyError, match="unknown PoleAngularClosureBase key",
        ):
            PoleAngularClosureBase.create("unknown_strategy")


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
    ψ-thread via :meth:`MorelMontryAngularSweep.compute_psi_half_per_level`
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
        :func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
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
    :meth:`MorelMontryAngularSweep.compute_psi_half_per_level` recurrence
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
    :func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
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

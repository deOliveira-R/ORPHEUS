r"""#257 S6 — Spec B: the fission ``production_rate`` reproduces ``F.apply`` bit-identically.

S6 adds a ``production_rate`` property to :class:`FissionOperator` →
the :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
over ``νΣf`` (the S2 accessor ``mat_xs.fission_production_field``; the
original S5 ``ProductionRateFunctional`` was retired — the reaction-rate
functional generalises it). With ``kernel`` (already present since Wave T)
the operator becomes an ``IntegralKernelOperator``. The matvec arms are
UNCHANGED.

This is Frame 3's "first test" — the *semantic composition* equals the
*fused realization*. The Frame-3 decomposition is
``F = M_χ ∘ ProductionRate ∘ M_νΣf``: the production-rate functional
contracts the group axis to the per-cell density ``p(r)``, and χ
re-broadcasts it across groups. The fused realization is
:meth:`FissionOperator.apply`, whose kernel reduces to
``χ · (νΣf · φ).sum(axis=0, keepdims=True)`` (the ``RankOneOperator.apply``
``inner`` line, ``operator.py:1776``). So:

    χ · production_rate.evaluate(φ)  ==  F.apply(φ)   (0 ULP)

because ``production_rate.evaluate(φ)`` IS the ``inner`` line byte-for-byte
(same numpy primitive, same axis, same ``keepdims=True``), and the χ
broadcast reproduces ``RankOneOperator``'s ``left * inner``.

vv claim layer (1.5 gate): this file carries TWO distinct claim types,
explicitly demarcated.

* B.1 CORRECTNESS — ``F.apply(φ)`` matches a STRUCTURALLY-INDEPENDENT
  hand-derived Python double-loop (``hand_derived_fission_emission``,
  shares NO numpy reduction with the production path). This pins that the
  fused realization computes the right physics. (foundation / value claim.)
* B.2 EQUIVALENCE / de-risk — the NEW ``production_rate`` property,
  re-broadcast by χ, reproduces ``F.apply(φ)`` at 0 ULP. Clearly
  demarcated as equivalence (the new property reproduces the existing
  realization), NOT a second correctness claim — L11: the structurally
  independent reference for the physics is B.1's hand loop, not this
  same-array comparison.

vv Mode-2 (νΣf↔φ / wrong axis): the inputs are a ≥3-group fissile cell,
heterogeneous, with asymmetric νΣf AND χ per group, so a swap or a
wrong-axis contraction is detectable (a νΣf↔φ swap produces a different
contraction; a wrong χ broadcast produces a different group profile).

vv Mode-11: B.2 reads ``op.production_rate`` OFF the live operator and
compares against ``F.apply`` — so mutating the production_rate property
(or the χ binding) reddens the gate. It does NOT route around the new
property via a sibling path.

vv Mode-8: every gate uses ``np.testing.*`` / ``require`` (function
calls, fire under ``python -O``) — NEVER a bare ``assert``.

``foundation`` — software invariants on the operator's type surface plus
an L0 value-correctness check against a hand-derived reference. No
theory-page ``:label:`` is claimed.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.scalar_flux import ScalarFlux

from tests.transport._integral_kernel_helpers import (
    hand_derived_fission_emission,
    require,
    require_production_rate_property,
)

pytestmark = pytest.mark.foundation


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def solver_4g():
    """≥3-group fissile, heterogeneous fixture.

    Mixture ``A`` 4g carries an asymmetric production cross section
    ``νΣf = [0.014, 0.026, 0.125, 0.245]`` and an asymmetric emission
    spectrum ``χ = [0.6, 0.35, 0.05, 0]`` — per-group asymmetric, so a
    νΣf↔φ swap or a wrong χ broadcast (Mode 2) is detectable. Mixture
    ``B`` 4g (non-fissile moderator) is the second region → heterogeneous
    νΣf field (H2: a flat field would null redistribution-style
    discrimination).
    """
    fuel = get_mixture("A", "4g")
    mod = get_mixture("B", "4g")
    nx, ny = 6, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:3, :] = 2  # fuel
    mat[3:, :] = 0  # moderator (no fission → asymmetric across the mesh)
    mesh = _uniform_2d(nx, ny, 0.2, mat)
    quad = Quadrature.lebedev(order=17)
    return SNSolver(SNMesh(mesh, quad, {2: fuel, 0: mod}))


def _asymmetric_phi(ng, nx, ny, seed=20260620):
    """A flux distinct per group AND per cell (no symmetry to hide a swap)."""
    rng = np.random.default_rng(seed)
    return rng.uniform(0.05, 1.0, size=(ng, nx, ny))


# ═══════════════════════════════════════════════════════════════════════
# B.1 — CORRECTNESS: F.apply matches a structurally-independent hand loop.
# ═══════════════════════════════════════════════════════════════════════


class TestFissionApplyCorrectness:
    def test_apply_matches_hand_derived_emission(self, solver_4g):
        """``F.apply(φ)`` equals the EXPLICIT-loop fission emission (0 ULP).

        The reference shares NO numpy reduction primitive with the
        production ``RankOneOperator.apply`` path — it is a per-cell,
        per-group Python double-loop. A νΣf↔φ swap, a wrong contraction
        axis, or a dropped χ broadcast disagrees with it.
        """
        op = solver_4g.fission_op
        ng = solver_4g.ng
        nx, ny = solver_4g.sn_mesh.spatial_shape
        phi = _asymmetric_phi(ng, nx, ny)

        out = op.apply(phi)  # bare-ndarray arm → (ng, nx, ny)
        expected = hand_derived_fission_emission(
            op.mat_xs.emission_spectrum, op.mat_xs.fission_production, phi,
        )
        np.testing.assert_allclose(
            out, expected, rtol=1e-13, atol=0.0,
            err_msg="FissionOperator.apply disagrees with the hand-derived "
            "χ·Σνσf·φ emission — a sign/axis/broadcast bug in the fused "
            "fission realization.",
        )

    def test_hand_loop_is_role_swap_sensitive(self, solver_4g):
        """Discriminator: the asymmetric inputs make a χ↔νΣf ROLE swap detectable.

        The production contraction ``Σ_g' νΣf_g'·φ_g'`` is symmetric in the
        PRODUCT, so swapping the two *contracted* arrays (νΣf↔φ) is a
        no-op — that is NOT the fission failure mode. The genuine Mode-2
        trap is the ROLE swap: using νΣf as the broadcast spectrum and χ
        as a contracted array (``χ·φ`` contracted, ``νΣf`` broadcast). With
        asymmetric per-group χ AND νΣf this produces a DIFFERENT emission,
        so the fixture genuinely constrains the χ-vs-νΣf assignment. If
        this row shows equality, the fixture lost its asymmetry and B.1 is
        blind to the role swap.
        """
        op = solver_4g.fission_op
        ng = solver_4g.ng
        nx, ny = solver_4g.sn_mesh.spatial_shape
        phi = _asymmetric_phi(ng, nx, ny)
        chi = op.mat_xs.emission_spectrum
        nu_sf = op.mat_xs.fission_production

        straight = hand_derived_fission_emission(chi, nu_sf, phi)
        # ROLE swap: broadcast by νΣf, contract χ·φ. χ ≠ νΣf per group →
        # a different per-group emission profile.
        role_swapped = hand_derived_fission_emission(nu_sf, chi, phi)
        require(
            not np.allclose(straight, role_swapped),
            "The fission fixture must be asymmetric enough that a χ↔νΣf "
            "ROLE swap changes the emission (Mode 2 discriminability). They "
            "agreed — the fixture lost its asymmetry; B.1 would be blind to "
            "the role swap.",
        )


# ═══════════════════════════════════════════════════════════════════════
# B.2 — EQUIVALENCE / de-risk (NOT correctness): the NEW production_rate
# property, re-broadcast by χ, reproduces F.apply BIT-IDENTICALLY.
# ═══════════════════════════════════════════════════════════════════════


class TestProductionRateReproducesApply:
    def test_chi_times_production_rate_is_fission_apply(self, solver_4g):
        """``χ · production_rate.evaluate(φ) == F.apply(φ)`` at 0 ULP.

        EQUIVALENCE / de-risk (NOT a second correctness claim — the
        structurally independent physics reference is B.1's hand loop).
        ``production_rate.evaluate(φ)`` IS the ``inner`` reduction line of
        ``RankOneOperator.apply`` (same numpy primitive, same axis,
        ``keepdims=True``), and ``χ · inner`` is the ``left * inner``
        broadcast — so the S6 composition is bit-identical to the legacy
        rank-1 fission by construction.

        Mode-11: reads ``op.production_rate`` OFF the live operator (the
        NEW property) and compares against the live ``F.apply`` — mutating
        the property reddens this gate.
        """
        op = solver_4g.fission_op
        ng = solver_4g.ng
        nx, ny = solver_4g.sn_mesh.spatial_shape
        phi = _asymmetric_phi(ng, nx, ny)

        pr = require_production_rate_property(op)  # NEW S6 member; skip if PRE-IMPL
        density = np.asarray(pr.evaluate(phi))  # (1, nx, ny) keepdims
        # χ broadcast reproduces RankOneOperator's `left * inner`.
        chi = op.mat_xs.emission_spectrum  # (ng, nx, ny)
        composed = chi * density  # (ng, nx, ny)

        fused = op.apply(phi)  # the unchanged matvec arm → (ng, nx, ny)
        np.testing.assert_array_equal(
            composed, fused,
            err_msg="χ · production_rate.evaluate(φ) must reproduce "
            "FissionOperator.apply(φ) BIT-IDENTICALLY (0 ULP) — the S6 "
            "composition F = M_χ ∘ ProductionRate ∘ M_νΣf shares the rank-1 "
            "`inner` reduction with the fused realization. A mismatch means "
            "the new production_rate property does not reproduce the arm.",
        )

    def test_production_rate_carries_nu_sigma_f(self, solver_4g):
        """The ``production_rate`` functional contracts ``νΣf`` (not Σt / Σa).

        Reads the functional's value at a known flux and confirms the
        density equals ``Σ_g νΣf_g · φ_g`` (the production cross section),
        NOT some other cross section — pins that S6 wired the
        ``mat_xs.fission_production_field`` accessor into the property.
        """
        op = solver_4g.fission_op
        ng = solver_4g.ng
        nx, ny = solver_4g.sn_mesh.spatial_shape
        phi = _asymmetric_phi(ng, nx, ny)

        pr = require_production_rate_property(op)
        density = np.asarray(pr.evaluate(phi)).reshape(nx, ny)
        expected = (op.mat_xs.fission_production * phi).sum(axis=0)
        np.testing.assert_array_equal(
            density, expected,
            err_msg="production_rate.evaluate must contract νΣf against φ "
            "over groups (the fission_production_field accessor). It did "
            "not match Σ_g νΣf_g·φ_g — wrong cross section wired in.",
        )


# ═══════════════════════════════════════════════════════════════════════
# B.4 — Mode-11 sentinel: F.apply ROUTES THROUGH the functional's evaluate.
# ═══════════════════════════════════════════════════════════════════════


class TestFissionApplyRoutesThroughFunctional:
    """The matvec ENTERS ``ReactionRateFunctional.evaluate`` — not a green twin.

    B.2's value-equivalence (``χ·pr.evaluate == F.apply``) holds even if the
    matvec computed its OWN inline reduction and never touched the functional —
    both sides equal ``χ·(νΣf·φ).sum``. This sentinel closes that Mode-11 gap:
    it wraps ``ReactionRateFunctional.evaluate`` in-process and asserts the
    counter advanced after ``F.apply``, proving the production-rate co-vector
    IS the contraction the kernel performs (the dyad's row-factor), not a
    parallel description. An inline-reduction regression leaves the counter at
    0 and reddens this gate (the strictly-stronger Mode-11 proof over B.2).
    """

    def test_apply_enters_reaction_rate_functional_evaluate(self, solver_4g, monkeypatch):
        from orpheus.transport.reaction_rate_functional import ReactionRateFunctional

        calls = {"n": 0}
        real_evaluate = ReactionRateFunctional.evaluate

        def counting_evaluate(self, x, /):
            calls["n"] += 1
            return real_evaluate(self, x)

        monkeypatch.setattr(ReactionRateFunctional, "evaluate", counting_evaluate)

        op = solver_4g.fission_op
        ng = solver_4g.ng
        nx, ny = solver_4g.sn_mesh.spatial_shape
        op.apply(_asymmetric_phi(ng, nx, ny))

        require(
            calls["n"] > 0,
            "F.apply did NOT route through ReactionRateFunctional.evaluate — "
            "the fission matvec computed an inline reduction instead of calling "
            "the dyad's row-factor (Mode-11: a green twin routing around the "
            "rewired reader).",
        )


# ═══════════════════════════════════════════════════════════════════════
# B.3 — the existing kernel gate MUST stay green unchanged. (Stated here as
# a pointer; the gate itself lives in test_fission_operator.py.)
# ═══════════════════════════════════════════════════════════════════════
#
# ``tests/sn/operators/test_fission_operator.py::TestRankOneTensorProductKernel``
# pins ``kernel.ops[0].reconstruction is emission_spectrum`` (reference
# identity) and ``kernel.apply ≡ apply`` (single source of truth). This carve
# must NOT perturb it — the kernel is the OTHER half of the
# IntegralKernelOperator surface and is already green.


# ═══════════════════════════════════════════════════════════════════════
# W-F live-arm sentinel — Mode-11: the K-eigenvalue loop EXECUTES the
# fission bare-``np.ndarray`` dispatch arm. Pins the load-bearing arm so
# the W-F dead-arm retirement cannot delete it.
# ═══════════════════════════════════════════════════════════════════════


class TestFissionNdarrayArmIsKEigenvalueLive:
    r"""The bare-``np.ndarray`` fission ``apply`` arm is on the K-loop call graph.

    W-E resolution (``d30d4a6``): the K-eigenvalue outer loop
    :func:`orpheus.numerics.eigenvalue.power_iteration` feeds a bare
    :class:`numpy.ndarray` flux to
    :meth:`~orpheus.sn.solver.SNSolver.compute_fission_source`, which calls
    ``self.fission_op.apply(flux_distribution) / keff`` (``sn/solver.py``).
    So the bare-ndarray dispatch arm (``fission.py`` —
    ``@_apply_impl.register def _(self, phi_arr: np.ndarray)``) is the LIVE
    arm at the outer-iteration boundary, NOT dead weight.

    This sentinel is the W-F safety net: it proves — by a real keff solve
    with an in-process counter on the *registered ndarray implementation
    function itself* — that the arm is genuinely executed, so the W-F
    dead-arm retirement keeps it (and a future regression that routes the
    K-loop around it reddens here).

    vv Mode-11 (the strictly-stronger proof): a green keff value alone does
    NOT prove the ndarray arm ran — a refactor could route fission through a
    typed carrier and still converge. The sentinel WRAPS the production
    reader (the ndarray arm) in-process and asserts the counter advanced —
    the routed-around path cannot fake the wrap (``vv-principles`` Mode-11
    "pytest-plugin sentinel that WRAPS the internal call").

    vv Mode-8: the gate uses ``require`` (a ``pytest.fail`` call) — fires
    under ``python -O``; NEVER a bare ``assert``.

    Why wrap the *registry* function and not ``F.apply``: wrapping the outer
    :func:`functools.singledispatchmethod` callable DEFEATS type-based
    dispatch (the wrapper's ``__class__`` is seen, the input falls to the
    base ``TypeError`` arm) — the exact Mode-11 hazard in miniature. The
    counter must wrap the ``np.ndarray``-registered leaf, reached via the
    descriptor's ``dispatcher.registry``.
    """

    @pytest.mark.sentinel
    def test_keff_solve_executes_fission_ndarray_arm(self, solver_4g):
        from orpheus.numerics.eigenvalue import power_iteration
        from orpheus.transport.operators.fission import FissionOperator

        # Reach the np.ndarray-registered implementation leaf (NOT the
        # dispatcher) and wrap it with a counter. The descriptor lives on the
        # class __dict__; `.dispatcher` is the underlying functools
        # singledispatch carrying the read-only `.registry` mappingproxy.
        descriptor = FissionOperator.__dict__["_apply_impl"]
        registry = descriptor.dispatcher.registry
        require(
            np.ndarray in registry,
            "FissionOperator._apply_impl has NO np.ndarray-registered arm — "
            "either the arm was already retired (the K-loop will TypeError) "
            "or the dispatch handle moved; W-F must keep this arm.",
        )
        ndarray_impl = registry[np.ndarray]

        calls = {"n": 0}

        def counting_ndarray_arm(self, phi_arr):
            calls["n"] += 1
            return ndarray_impl(self, phi_arr)

        # The registry is a read-only mappingproxy — the ONLY supported
        # mutation is `singledispatchmethod.register`. Re-register the
        # counting wrapper for np.ndarray, run the solve, then restore the
        # original leaf in `finally` (manual revert because the registry is
        # class-global and there is no monkeypatch hook for it; NEVER leave
        # the live dispatch table mutated — Mode-11 probe hygiene).
        descriptor.register(np.ndarray, counting_ndarray_arm)
        try:
            keff, history, _flux = power_iteration(solver_4g, max_iter=60)
        finally:
            descriptor.register(np.ndarray, ndarray_impl)

        require(
            calls["n"] > 0,
            "The K-eigenvalue loop did NOT execute the fission bare-ndarray "
            "apply arm — power_iteration converged WITHOUT routing fission "
            "through `F.apply(np.ndarray)`. Either the live arm was deleted "
            "(W-F over-retirement) or the K-loop now feeds a typed carrier. "
            "The ndarray arm is the load-bearing outer-iteration boundary "
            "(W-E `d30d4a6`); W-F must not retire it (Mode-11).",
        )
        # Corroborating sanity: the solve actually ran (so the counter result
        # is meaningful, not a zero-iteration vacuum).
        require(
            len(history) > 0 and np.isfinite(keff),
            f"power_iteration returned no usable trajectory "
            f"(keff={keff}, n_outer={len(history)}) — the sentinel's "
            f"counter reading is not anchored to a real solve.",
        )

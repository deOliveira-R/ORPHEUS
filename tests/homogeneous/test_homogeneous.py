"""Verify the infinite-medium eigenvalue solver against SymPy analytical solutions."""

import numpy as np
import pytest

import orpheus.numerics.eigenvalue as _eig
from orpheus.derivations import get
from orpheus.homogeneous.solver import solve_homogeneous_infinite

# File-level verifies marker: every test in this file exercises the
# homogeneous eigenvalue chain end-to-end by asserting k_inf matches
# the SymPy-derived analytical reference to 1e-12. That tolerance is
# tight enough to pin every step of the derivation — if any of the
# labelled equations below were implemented incorrectly the k mismatch
# would be far larger than 1e-12, so a passing test is equation-level
# (L1) verification for every link in the chain.
#
# Declared explicitly here (rather than inherited from
# VerificationCase) so the Nexus AST pass picks it up via decorator
# parsing and writes TESTS edges.
#
# The 2G labels (two-group-*) and the power-iteration step labels
# (fission-source, fixed-source-solve, keff-update) are all exercised
# by the homo_2eg / homo_4eg cases — the analytical k_inf is derived
# symbolically via exactly those equations, so a solver k that matches
# to 1e-12 implies every link in the chain is correct. The absorption-
# xs label is the derived property used inside keff-update.
pytestmark = [pytest.mark.l1, pytest.mark.verifies(
    "one-group-kinf",
    "inf-hom-balance",
    "matrix-eigenvalue",
    "removal-matrix",
    "fission-matrix",
    "mg-balance",
    # B.1 additions (issue #87): the full 2G analytical chain and the
    # power-iteration step labels, all verified end-to-end by the
    # homo_2eg and homo_4eg parametric cases.
    "two-group-A",
    "two-group-F",
    "two-group-Ainv",
    "two-group-M",
    "two-group-charpoly",
    "two-group-roots",
    "fission-source",
    "fixed-source-solve",
    "keff-update",
    "absorption-xs",
)]


@pytest.mark.parametrize("case_name", [
    "homo_1eg",
    "homo_2eg",
    "homo_4eg",
    # Non-trivial, asymmetric (n,2n): de-vacuums the n2n-in-A convention.
    # Every other case has Sig2=0, so the 2·Σ₂ᵀ loss term is never
    # exercised; here it moves k_inf ~0.6 (a drop/double-count reds).
    "homo_2eg_n2n",
])
def test_kinf_exact(case_name):
    """Eigenvalue must match analytical solution to machine precision.

    The refounded solver assembles A = C − K_iso from the transport
    operators (a different FP reduction tree than the oracle's fused
    ``(Σ_s + 2Σ_2)ᵀ``), so the tolerance is now principled-equivalence
    (FP-non-associativity, ~1 ULP), not bit-identity — still ≪ 1e-12.
    """
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    result = solve_homogeneous_infinite(mix)
    assert abs(result.k_inf - case.k_inf) < 1e-12, (
        f"k_inf mismatch: solver={result.k_inf:.10f} "
        f"analytical={case.k_inf:.10f}"
    )


@pytest.mark.verifies("normalisation")
def test_post_solve_production_rate_is_100():
    """L1: post-convergence flux is normalised to 100 n/cm^3/s production.

    After :func:`solve_homogeneous_infinite` solves, the flux is rescaled
    so the **fission** production rate

    .. math::

       \\nu\\Sigma_f \\cdot \\boldsymbol{\\phi} = 100

    (see Eq. ``normalisation`` in docs/theory/homogeneous.rst).
    Production is :math:`\\nu\\Sigma_f` only — the (n,2n) reaction is a
    loss-side transfer folded into the loss matrix as
    :math:`2\\Sigma_2^T`, NOT a production channel (it does NOT enter
    the numerator).  The ``homo_2eg_n2n`` case (non-zero, asymmetric
    :math:`\\Sigma_2`) makes this non-vacuous: under the retired
    ``(\\Sigma_p + 2\\cdot\\text{colsum}(\\Sigma_2))`` formula the
    production would not equal 100.  The 1G case is a degenerate
    one-scalar normalisation a bug could accidentally satisfy, so it is
    excluded.
    """
    for case_name in ("homo_2eg", "homo_4eg", "homo_2eg_n2n"):
        case = get(case_name)
        mix = next(iter(case.materials.values()))
        result = solve_homogeneous_infinite(mix)

        # Production = νΣ_f @ φ  (fission only; n2n lives in A, not F).
        production = mix.SigP @ result.flux

        assert abs(production - 100.0) < 1e-9, (
            f"{case_name}: production rate = {production:.6e}, "
            f"expected 100.0 (normalisation constraint)"
        )


# ── Operator-algebra assembly: A-level oracle + Mode-11 liveness (#276) ──


@pytest.mark.verifies("removal-matrix")
def test_assemble_loss_matrix_matches_fused_oracle():
    """The operator-composed A = C − K_iso (apply-to-basis) matches the fused
    ``diag(Σ_t) − (Σ_s0 + 2Σ_2)ᵀ`` on the non-trivial-(n,2n) case.

    A SHARP procedural pin at the A level — it localises a sign/term/omission
    bug in the operator assembly faster than the end-to-end eig. It shares
    ``mat_xs`` data with the fused form (so it is NOT structurally
    independent), and therefore PAIRS with the SymPy ``case.k_inf`` anchor
    in :func:`test_kinf_exact` rather than replacing it.
    """
    from orpheus.homogeneous.solver import _assemble_loss_matrix
    from orpheus.transport.mesh.material_mesh import MaterialMesh

    case = get("homo_2eg_n2n")
    mix = next(iter(case.materials.values()))
    mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()

    A = _assemble_loss_matrix(mat_xs)
    sig_t = mat_xs.total_cross_section[:, 0]
    sig_s0 = mat_xs.sig_s_legendre(0)[0]  # (ng, ng), [g_from, g_to]
    sig_2 = mat_xs.n2n_matrix(0)
    A_fused = np.diag(sig_t) - (sig_s0 + 2.0 * sig_2).T
    np.testing.assert_allclose(A, A_fused, atol=1e-12, rtol=0)


def test_kinf_gate_executes_the_bare_multiplication_arm(monkeypatch):
    """Mode-11: the homogeneous k∞ gate actually EXECUTES the new
    ``MultiplicationOperator`` bare-ndarray arm.

    The apply-to-basis A-assembly routes the collision diagonal C = M[Σ_t]
    through the bare arm. Perturbing ONLY that arm (×1.5 on ndarray input)
    moves k_inf O(1) — proving the arm is on the gate's call graph and
    load-bearing, not a vacuous green. (``-O``-safe: the monkeypatch is an
    in-process attribute swap, reverted by the fixture; never a
    ``git checkout``.)
    """
    from orpheus.transport.operators.multiplication_operator import (
        MultiplicationOperator,
    )

    raw = MultiplicationOperator.__dict__["apply"]  # the singledispatchmethod

    def perturbed(self, x):
        out = raw.__get__(self, type(self))(x)
        if isinstance(x, np.ndarray):  # corrupt ONLY the meshless collision arm
            return out * 1.5
        return out

    monkeypatch.setattr(MultiplicationOperator, "apply", perturbed)

    case = get("homo_2eg_n2n")
    mix = next(iter(case.materials.values()))
    result = solve_homogeneous_infinite(mix)
    assert abs(result.k_inf - case.k_inf) > 1e-3, (
        f"perturbing the bare M[Σ_t] arm left k_inf at {result.k_inf:.6f} "
        f"(oracle {case.k_inf:.6f}) — the homogeneous gate does NOT execute "
        f"the bare arm (Mode-11 vacuous-green)"
    )


# ── #276 P4-D: solve_homogeneous_infinite rewired to DirectEigenvalue ──
#
# These gates pin the rewire of the inline eig (solver.py ~182-187) onto the
# new DirectEigenvalue primitive and the rerouting of the reaction rates
# through IntegratedReactionRate. Expected equivalence classes:
#
#   * k_inf  — BIT-IDENTICAL. The eig computation is structurally unchanged
#     (same ``eig(solve(A, F))``, argmax, real(), sign flip), so routing it
#     through the verified primitive must not move a single bit. The
#     structurally-INDEPENDENT correctness anchor stays ``test_kinf_exact``
#     (SymPy ``case.k_inf``, 1e-12) — DirectEigenvalue must NOT be wired into
#     that oracle path.
#   * rates / flux via IntegratedReactionRate — BIT-IDENTICAL on the shipped
#     cases (ng ≤ 4): V_cell = 1 and the short Σ_g νΣf_g φ_g reduction matches
#     ``νΣf @ φ`` exactly. (A larger group structure would relax this to ≤ few
#     ULP — a reduction-tree change per vv-principles bit-identity criterion 3.)


def _require(condition: bool, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare assert) for the liveness sentinel."""
    if not condition:
        pytest.fail(message)


def test_direct_eigenvalue_is_on_the_homogeneous_call_path(monkeypatch):
    """Mode-11: after the rewire, ``solve_homogeneous_infinite`` actually CALLS
    the direct-eigenvalue primitive — not a still-live inline ``eig``.

    In-process WRAP sentinel (the gold-standard Mode-11 proof): spies the
    direct-eigenvalue symbol in the SOLVER's own namespace (so a
    ``from ... import DirectEigenvalue`` local binding is the thing wrapped) and
    asserts the counter fires during the solve. A routed-around inline ``eig``
    cannot bump it. SKIPS pre-impl (symbol not yet imported into the solver).
    """
    import orpheus.homogeneous.solver as hsolver

    # The rewire SHOULD import the primitive at module level (top of solver.py)
    # so the solver's own binding is wrappable here; ``hasattr(hsolver, ...)`` is
    # the rewired-yet? proxy. We also wrap the definition site ``_eig`` to catch
    # an in-function ``from ... import`` (belt-and-suspenders). SKIP only when
    # the solver's namespace carries no such symbol (not yet rewired).
    target = next(
        (n for n in ("DirectEigenvalue", "direct_eigenvalue") if hasattr(hsolver, n)),
        None,
    )
    if target is None:
        pytest.skip(
            "#276 P4-D PRE-IMPL: solve_homogeneous_infinite not yet rewired to "
            "DirectEigenvalue (symbol absent from its namespace)."
        )

    calls: list[int] = []

    def _wrap(ns: object, name: str) -> None:
        original = getattr(ns, name)

        def spy(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        monkeypatch.setattr(ns, name, spy)

    _wrap(hsolver, target)
    if hasattr(_eig, target):  # also wrap the definition site (in-function import)
        _wrap(_eig, target)

    case = get("homo_2eg_n2n")
    mix = next(iter(case.materials.values()))
    result = solve_homogeneous_infinite(mix)
    _require(
        len(calls) >= 1,
        "solve_homogeneous_infinite did NOT call the direct-eigenvalue primitive "
        "— the rewire is claimed but the inline eig is still live (Mode-11 "
        "vacuous green).",
    )
    np.testing.assert_allclose(result.k_inf, case.k_inf, atol=1e-12, rtol=0)


def test_kinf_is_the_direct_eigenvalue_of_the_assembled_pair():
    """Bit-identity: solver ``k_inf`` == the primitive's ``k`` on the SAME (A, F).

    Builds A, F the way the solver does and asserts ``solve_homogeneous_infinite``'s
    ``k_inf`` is BYTE-equal to ``DirectEigenvalue(A, F).solve()[0]`` — proving
    (a) the rewired solver's k_inf IS the primitive's output (it routes THROUGH
    it, no drift) and (b) the primitive reproduces the inline eig bit-for-bit.
    PAIRS with ``test_kinf_exact`` (the structurally-independent SymPy anchor);
    this pin localises a rewire regression to the primitive boundary. SKIPS
    pre-impl.
    """
    DE = getattr(_eig, "DirectEigenvalue", None)
    fn = getattr(_eig, "direct_eigenvalue", None)
    if DE is None and fn is None:
        pytest.skip(
            "#276 P4-D PRE-IMPL: DirectEigenvalue / direct_eigenvalue not yet on "
            "orpheus.numerics.eigenvalue."
        )

    from orpheus.homogeneous.solver import _as_dense, _assemble_loss_matrix
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.operators.fission import FissionOperator

    case = get("homo_2eg_n2n")
    mix = next(iter(case.materials.values()))
    mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()
    A = _assemble_loss_matrix(mat_xs)
    F = _as_dense(FissionOperator.from_solver_data(mat_xs=mat_xs), mix.ng)
    if DE is not None:
        k_prim = DE(A, F).solve()[0]
    else:
        assert fn is not None  # guaranteed non-None by the pre-impl skip above
        k_prim = fn(A, F)[0]

    result = solve_homogeneous_infinite(mix)
    _require(
        result.k_inf == k_prim,
        f"solver k_inf={result.k_inf!r} != primitive k={k_prim!r} — the rewire "
        f"is not bit-identical (or the solver routes around the primitive).",
    )


def test_rates_via_integrated_reaction_rate_are_bit_identical():
    r"""Bit-identity of the rate rerouting: ``IntegratedReactionRate(νΣf).evaluate(φ)``
    == ``νΣf @ φ`` on the meshless unit-volume cell, for every shipped case.

    ``V_cell = 1`` and the short ``Σ_g νΣf_g φ_g`` reduction matches the dot
    bit-for-bit for ng ∈ {1, 2, 4}. Pins the production-rate rerouting
    independently of the eig; runs green TODAY (both paths exist) and guards
    the rewire. (No PRE-IMPL skip — it protects the rate-side claim regardless
    of the eigenvalue-side rewire state.)
    """
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.transport.reaction_rate_functional import IntegratedReactionRate

    for case_name in ("homo_1eg", "homo_2eg", "homo_4eg", "homo_2eg_n2n"):
        case = get(case_name)
        mix = next(iter(case.materials.values()))
        mat_xs = MaterialMesh.from_materials({0: mix}).material_xs_field()
        ng = mix.ng
        phi = solve_homogeneous_infinite(mix).flux

        legacy_prod = float(mat_xs.fission_production[:, 0] @ phi)
        irr_prod = float(
            IntegratedReactionRate(mat_xs.fission_production_field).evaluate(
                phi.reshape(ng, 1)
            )
        )
        _require(
            legacy_prod == irr_prod,
            f"{case_name}: IntegratedReactionRate production {irr_prod!r} != "
            f"νΣf@φ {legacy_prod!r} — the rate rerouting is not bit-identical.",
        )

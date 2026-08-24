r"""GreenOperator L1 SN gates — ``A_loss = (L+C) − S`` on a heterogeneous 2G
VACUUM slab (taxonomy §12 step 4; verification spec PART III §21.2/§22/§24.4).

**What is pinned.**  The generic sum inverse on the REAL transport operators:
``((L+C) − S).inverse()`` is a
:class:`~orpheus.numerics.green_operator.GreenOperator` whose derived
splitting is EXACTLY the canonical within-group source iteration —
preconditioner = the WDD sweep ``(L+C)⁻¹``, gain = ``S`` — proven by (a) the
G-Neumann multiple-scattering expansion :math:`\sum_k ((L+C)^{-1}S)^k
(L+C)^{-1} q` converging to ``green.apply(q)`` with the geometric tail ratio
of the physical scattering ratio (§17 falsifier-5, run as a permanent gate),
(b) an operator-level MANUFACTURED anchor (exact, structurally independent
of the iteration), and (c) BIT-IDENTITY of the auto-derived driver against
the hand-assembled ``SourceIteration(sweep, S)``.

**Config (Mode-9 discipline, spec §21.2):** heterogeneous (two materials,
distinct σ_t and ASYMMETRIC 2G scattering blocks) VACUUM slab — a reflective
isotropic box would null the streaming↔scattering redistribution the Neumann
series expands, and 1G is blind to a scattering-matrix transpose (Mode 6).

**The manufactured anchor + the #284 source subspace.**  Sweep-realized
inverses are exact on the SOURCE subspace (rhs out-rows zero); a dense LU of
the materialized sum would disagree on arbitrary trace rows.  The anchor
therefore manufactures at the operator level: ``x_tc = (L+C)⁻¹(random)`` is
trace-CONSISTENT, so ``q = A_loss.apply(x_tc)`` lies in the source subspace
and the exact solution IS ``x_tc`` — no dense oracle needed, no trace
caveat.  ``q`` is honestly SOURCE-typed (``AngularSourceSink`` bulk), the
production rhs convention.

**The typed cold start (shared driver wart, documented).**  On role-typed
SN carriers ``initial_guess=None`` is unusable — the driver's
``_zeros_like(q)`` inherits ``q``'s SOURCE role and the gain ``S`` rejects a
source-typed iterate (cross-class arithmetic is forbidden).  Production
always warm-starts with a flux-typed state; these gates do the same
(``TimedFullField.zeros(interior=AngularFlux, …)``).  This is the DRIVER's
pre-existing typing boundary, inherited by Green — not a step-4 novelty.

**Ordering-ruling edges (§22.1):** ``(L+C)`` → ``SweepOperator`` (the MRO
shadow — the fused type owns the DIRECT inverse); ``A_loss`` → Green,
converges; ``C + L`` (legal spelling, invertible leading collision term) →
Green constructs, its collision-preconditioned Richardson DIVERGES →
:class:`ConvergenceFailure`, never a silent wrong answer; ``(−S) + (L+C)``
→ refused at construction naming the canonical ordering.

**Teeth (spec §25, mutation-verified 2026-07-02):** M-GRN-FLATTEN (flatten
through the ``StreamingCollisionOperator`` subclass → the leading term becomes the
bare non-invertible ``L`` → the ``A_loss`` gates RED at construction);
M-GRN-SIGN/SWAP RED the anchor + Neumann; M-GRN-TOL REDs the ``C+L`` raise.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.green_operator import ConvergenceFailure, GreenOperator
from orpheus.numerics.iteration import SourceIteration, seeded_inverse
from orpheus.numerics.operator import (
    NotInvertible,
    OperatorSum,
    ScaledOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingCollisionOperator, StreamingOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn.operators.test_removal_form_matvec_sweep import _random_state

pytestmark = pytest.mark.foundation

#: §20.0 driver-tol rule (iterative inverse — SAFETY × tol, not nulp).
_SAFETY = 10.0
_TOL = 1e-8


def _mix(sig_t, p0):
    """2G mixture with an ASYMMETRIC P0 block (Mode-6: transpose-detectable)."""
    m = make_mixture(
        sig_t=np.array(sig_t), sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=np.array(p0),
    )
    m.SigS = [csr_matrix(np.array(p0))]
    m.Sig2 = csr_matrix(np.zeros((2, 2)))
    return m


def _het_scattering_slab() -> SNMesh:
    """Two-material 2G VACUUM slab, GL-4 — het σ_t AND het asymmetric SigS."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=0, outer_thickness_cm=1.0),
            Region(mat_id=1, outer_thickness_cm=2.0),
        ),
        bcs=(BC("vacuum"), BC("vacuum")),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=3), RegionMesh(n_cells=3)),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=4),
        {0: _mix([1.0, 1.5], [[0.38, 0.10], [0.05, 0.60]]),
         1: _mix([1.2, 1.8], [[0.55, 0.03], [0.12, 0.40]])},
    )


def _operators():
    """``(sn, L+C fused, S, A_loss = (L+C) − S)`` on the het slab."""
    sn = _het_scattering_slab()
    mat_xs = sn.material_xs_field()
    lc = StreamingOperator(sn) + MultiplicationOperator.from_mesh(
        mat_xs.total_cross_section, sn,
    )
    S = ScatteringOperator(
        mat_xs=mat_xs, quadrature=sn.quad, scattering_order=0,
    )
    return sn, lc, S, lc - S


def _flux_zeros(sn: SNMesh) -> TimedFullField:
    """The flux-typed cold start (the production warm-start convention)."""
    return TimedFullField.zeros(interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn.full_field_space)


def _manufactured(sn, lc, a_loss):
    """Trace-consistent manufactured pair ``(x_tc, q = A_loss x_tc)``.

    ``x_tc`` is a sweep output (trace-consistent by the shed), so ``q``
    lies in the sweep's source subspace (#284) and the EXACT solution of
    ``A_loss ψ = q`` is ``x_tc`` — the structurally-independent anchor.
    """
    x_tc = lc.inverse().apply(_random_state(sn, seed=11))
    return x_tc, a_loss.apply(x_tc)


def test_g_neumann_l1_expansion_and_manufactured_anchor():
    """G-Neumann on the REAL operators (§21.2): the multiple-scattering
    partial sums ``Σ_k ((L+C)⁻¹S)^k (L+C)⁻¹ q`` converge to
    ``green.apply(q)`` AND both to the exact manufactured ``x_tc``; the
    tail decay is geometric with ratio < 1 (the het scattering ratio —
    0.4015 on this fixture), STABILIZED (a non-geometric tail would mean
    the splitting is not the physical one)."""
    sn, lc, S, a_loss = _operators()
    x_tc, q = _manufactured(sn, lc, a_loss)
    green = a_loss.inverse()
    assert isinstance(green, GreenOperator)
    psi = green.apply(q, initial_guess=_flux_zeros(sn))

    # rtol + atol: the promise is a NORM statement (true residual), so
    # small-magnitude elements carry an absolute band of the same size —
    # element-wise rtol alone over-demands on near-zero flux entries.
    np.testing.assert_allclose(
        psi.interior.values, x_tc.interior.values,
        rtol=_SAFETY * _TOL, atol=_SAFETY * _TOL,
        err_msg="Green ≠ the exact manufactured solution (anchor)",
    )

    lc_inv = lc.inverse()
    term = lc_inv.apply(q)
    acc = term.interior.values.copy()
    norms = [float(np.linalg.norm(term.interior.values))]
    for _k in range(1, 40):
        term = lc_inv.apply(S.apply(term))     # ((L+C)⁻¹S)^k (L+C)⁻¹ q
        acc = acc + term.interior.values
        norms.append(float(np.linalg.norm(term.interior.values)))

    np.testing.assert_allclose(
        acc, psi.interior.values, rtol=_SAFETY * _TOL, atol=_SAFETY * _TOL,
        err_msg="Neumann partial sums ≠ Green application (G-Neumann L1)",
    )
    tail = [norms[i + 1] / norms[i] for i in (30, 34, 38)]
    if not all(r < 0.6 for r in tail):
        pytest.fail(f"Neumann tail not contracting below the physical c: {tail}")
    if not max(tail) - min(tail) < 1e-3:
        pytest.fail(f"Neumann tail not geometric (splitting suspect): {tail}")


def test_green_driver_bit_identical_to_hand_source_iteration():
    """The flatten derivation reproduces the CANONICAL within-group SI —
    ``green._driver`` ≡ hand-built ``SourceIteration(sweep, S)`` BIT-
    identically.  Pinned at the DRIVER grain: ``green.apply`` adds the
    true-residual refinement layer on top (spec §18.A), so the promise
    surface may legitimately run extra steps; the derivation — which
    preconditioner, which gains, in what form — is what this gate owns
    (its anchor for "was right" is the manufactured gate above)."""
    sn, lc, S, a_loss = _operators()
    _, q = _manufactured(sn, lc, a_loss)
    green = GreenOperator(a_loss, max_iter=1000, tol=_TOL)
    hand = SourceIteration(
        seeded_inverse(lc), S, max_iter=1000, tol=_TOL,
    )
    got, _ = green._driver.solve(q, initial_guess=_flux_zeros(sn))
    ref, _ = hand.solve(q, initial_guess=_flux_zeros(sn))
    np.testing.assert_array_equal(
        got.interior.values, ref.interior.values,
        err_msg="auto-derived splitting ≠ the canonical hand-built SI",
    )


def test_ordering_dispatch_sweep_vs_green():
    """§22.1/§24.4 — the operand spelling selects the algorithm, keyed by
    TYPE: the fused ``(L+C)`` owns the DIRECT sweep (its ``.inverse()``
    override shadows the generic sum's by MRO); the plain sums route to
    the iterative Green."""
    _, lc, _s_unused, a_loss = _operators()
    del _s_unused
    assert isinstance(lc, StreamingCollisionOperator)
    assert type(lc.inverse()) is SweepOperator      # the MRO shadow holds
    assert type(a_loss) is OperatorSum              # −S broke the fusion
    assert type(a_loss.inverse()) is GreenOperator  # generic arm → Green
    assert a_loss.inverse().inverse() is a_loss     # involution, object id


def test_minus_s_plus_lc_refuses_with_canonical_ordering_message():
    """§22.3 — ``(−S) + (L+C)`` spells the non-invertible term first; the
    factory refuses at construction naming the canonical ordering (the
    §18.B contract pin: ``is_invertible`` reads the LEADING term)."""
    _, lc, S, _ = _operators()
    backwards = ScaledOperator(-1.0, S) + lc
    assert backwards.is_invertible is False
    with pytest.raises(NotInvertible, match="canonical ordering"):
        backwards.inverse()


def test_c_plus_l_trap_constructs_then_fails_loudly():
    """§22.2, the order-dependent strategy trap on the REAL operators:
    ``C + L`` is the same math as ``L + C`` but the spelling misses the
    fusion dispatch (it lives on ``StreamingOperator.__add__``, #261) —
    the plain sum's leading term ``C`` IS invertible, so a Green
    CONSTRUCTS (the algebra cannot know ρ(C⁻¹L) > 1)… and its
    collision-preconditioned Richardson fails LOUDLY at apply.  Never a
    silent wrong answer — the ruling's whole point."""
    sn, lc, S, a_loss = _operators()
    mat_xs = sn.material_xs_field()
    C = MultiplicationOperator.from_mesh(mat_xs.total_cross_section, sn)
    L = StreamingOperator(sn)
    cl = C + L
    assert type(cl) is OperatorSum          # no fusion on this spelling
    assert cl.is_invertible is True         # leading C invertible — constructs
    assert isinstance(cl.inverse(), GreenOperator)
    _, q = _manufactured(sn, lc, a_loss)
    with pytest.raises(ConvergenceFailure):
        GreenOperator(cl, max_iter=8, tol=_TOL).apply(
            q, initial_guess=_flux_zeros(sn),
        )

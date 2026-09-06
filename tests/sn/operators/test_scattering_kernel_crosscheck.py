r"""#257 S6 / CS4c step 5 — the FUSED angular route ≡ the TYPED moment route.

:attr:`ScatteringOperator.kernel` is the :math:`\ell\ge1` redistribution as a
typed ``OperatorProduct``. On the ANGULAR end it is the frame conjugation
:math:`R\circ\Lambda\circ M` — one fused ndarray chain; on the MOMENT end
(``S.on_moment_domain()``, the windowed SI driver's binding) the production
body is the explicit TYPED grid path: :math:`\Lambda` maps flux moments to
source moments (the role-changing edge, in the signature), then the minted
source-reconstruction FACE synthesises the per-ordinate source. Two different
spellings of one operator; this file cross-checks them.

⛔ **What this file used to compare, and why it changed.** Until CS4c step 5
the second side was the private helper ``_aniso_source_from_moment_values``,
reached by handing a moment iterate to the ANGULAR-bound operator, which
dispatched on the carrier's class per call (`[M]` 143 such feeds per windowed
solve — the shipped non-endomorphism the step-0 census measured). That helper
is RETIRED: *each binding acts through the body its ends select*, so the moment
operand now rides an operator bound on the moment end, and the comparison moved
UP a tier — from a private helper against a private chain, to the two BOUND
operators' own public actions. The claim did not weaken; the second side is now
the production route rather than a fragment of it.

    S.apply(FullField(ψ, 0))          ==  S_w.apply(FullField(M·ψ, 0))    (0 ULP)
    S.kernel.apply(ψ.values)          ==  S_w.kernel.apply((M·ψ).values)  (0 ULP)

`[M]` 2026-09-04, this fixture, 200 seeds (``default_rng(s)``, standard-normal
ψ): **200/200 ``array_equal``, max |Δ| = 0.0** over the whole sweep — so
``array_equal`` is a property of the FIXTURE, not of one draw (`vv` anti-#31),
and the two routes genuinely share :math:`\Lambda` and the frame's :math:`R`.

CRITICAL DEMARCATION (L11): this is the EQUIVALENCE / de-risk leg — the two
spellings compose the same numpy chain. It is NOT a correctness claim about
the scattering *physics*. The structurally-independent reference for that is
the aniso MMS gate
(``tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py``, both
l0 ``verifies("pn-scatter", "flux-moments")`` cases) — that gate stays green
and is the L1 backing.

The config is ANISOTROPIC (≥P1) + heterogeneous + multi-group so :math:`\Lambda`'s
ℓ≥1 blocks are genuinely exercised (a P0-only config would null the aniso path
and make the cross-check vacuous — guarded by the non-degeneracy row, which
also asserts the two ends' ℓ≥1 halves are non-zero: a zero morphism satisfies
every equality row).

vv Mode-11: both sides are read OFF live operators (``S`` and the sibling
``S.on_moment_domain()`` returns) — mutating either body reddens the gate; no
row routes around a production property.

vv Mode-8: every gate uses ``np.testing.*`` / ``require`` (function calls, fire
under ``python -O``) — NEVER a bare ``assert``.

``foundation`` — software invariants on the operator's type surface plus a
bit-identity equivalence check. No theory-page ``:label:`` is claimed (the
physics ``:label:`` s ``pn-scatter`` / ``flux-moments`` are pinned by the MMS
gate, NOT here).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_flux import AngularFlux

from tests.sn.operators._composite_operand import bulk_apply, zero_trace_composite
from tests.transport._integral_kernel_helpers import (
    require,
    require_scattering_kernel_property,
)

pytestmark = pytest.mark.foundation


def _uniform_2d(nx, ny, delta, mat_map):
    return Mesh2D(
        edges_x=np.linspace(0, nx * delta, nx + 1),
        edges_y=np.linspace(0, ny * delta, ny + 1),
        mat_map=np.asarray(mat_map, dtype=int),
    )


@pytest.fixture
def solver_p1_het():
    """ANISOTROPIC (P1) + heterogeneous + 2G scattering fixture.

    Two materials with DIFFERENT asymmetric P0 + P1 blocks so the moment
    chain ``R∘Λ∘M`` exercises ℓ≥1 with a genuinely heterogeneous,
    asymmetric ``SigS`` (H2 + Mode 6). Returns the SOLVER (the operator is
    ``solver.scattering_op``) so the typed :class:`AngularFlux` is built on
    the solver's own :class:`SNMesh`, exactly like the existing
    ``test_scattering_operator.py`` fixtures.
    """
    p0_a = np.array([[0.38, 0.10], [0.05, 0.90]])
    p1_a = np.array([[0.02, 0.01], [0.00, 0.04]])
    p0_b = np.array([[0.55, 0.03], [0.12, 0.40]])
    p1_b = np.array([[0.06, 0.02], [0.01, 0.03]])

    def _mix(p0, p1):
        m = make_mixture(
            sig_t=np.array([0.5, 1.0]),
            sig_c=np.array([0.01, 0.02]),
            sig_f=np.array([0.0, 0.0]),
            nu=np.array([0.0, 0.0]),
            chi=np.zeros(2),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
            sig_s=p0,
        )
        m.SigS = [csr_matrix(p0), csr_matrix(p1)]
        m.Sig2 = [csr_matrix(np.array([[0.0, 0.03], [0.01, 0.0]]))]
        return m

    nx, ny = 4, 3
    mat = np.zeros((nx, ny), dtype=int)
    mat[:2, :] = 0
    mat[2:, :] = 1
    mesh = _uniform_2d(nx, ny, 0.4, mat)
    quad = Quadrature.lebedev(order=17)
    sn_mesh = SNMesh(mesh, quad, {0: _mix(p0_a, p1_a), 1: _mix(p0_b, p1_b)})
    return SNSolver(sn_mesh, scattering_order=1)


def _aniso_psi(solver, seed=20260620):
    """A non-isotropic angular flux (so ℓ≥1 moments are non-zero)."""
    N = solver.quad.N
    ng = solver.ng
    nx, ny = solver.sn_mesh.spatial_shape
    rng = np.random.default_rng(seed)
    psi_values = rng.uniform(0.05, 1.0, size=(N, ng, nx, ny))
    return AngularFlux(values=psi_values, space=solver.sn_mesh.angular_bulk_space)


# ═══════════════════════════════════════════════════════════════════════
# C — EQUIVALENCE / de-risk: the NEW kernel reproduces the EXISTING
# R∘Λ∘M aniso realization byte-for-byte.
# ═══════════════════════════════════════════════════════════════════════


class TestFusedAngularRouteEqualsTypedMomentRoute:
    r"""The two ends of the retained analysis face compute ONE operator.

    ``S`` is bound on the angular composite (its domain interior is the
    analysis face's DOMAIN, so the operand is :math:`\psi`);
    ``S.on_moment_domain()`` is the SAME datum and the SAME faces bound to
    consume :math:`M\psi` (the analysis face's CODOMAIN). Their bodies are
    selected at construction and are different code — the fused
    ``frame.conjugate`` product vs the typed
    :math:`\Lambda`-then-reconstruction-face grid path — so agreement is a
    claim, not a tautology.
    """

    def test_full_action_bit_identical_across_the_two_ends(self, solver_p1_het):
        r"""``S.apply(ψ ⊕ 0) == S_w.apply(Mψ ⊕ 0)`` on the interior, 0 ULP.

        THE crosscheck, at the operator tier: the angular end's fused
        :math:`R\Lambda M\psi/W + E\phi/W` against the moment end's typed
        :math:`R(\Lambda\, M\psi)/W + E\phi_{\ell=0}/W`. Two selected bodies,
        one datum, one frame. `[M]` 200/200 seeds bit-equal (module
        docstring), so ``array_equal`` is the honest contract.

        Replaces ``test_kernel_apply_equals_existing_R_Lambda_M_chain``,
        whose second side (``_aniso_source_from_moment_values``) is retired.
        """
        op = solver_p1_het.scattering_op
        sn_mesh = solver_p1_het.sn_mesh
        psi = _aniso_psi(solver_p1_het)
        moments = op.flux_analysis.apply(psi)          # M·ψ, TYPED

        via_angular = bulk_apply(op, psi).values
        via_moment = op.on_moment_domain().apply(
            zero_trace_composite(moments, op.domain.trace_space),
        ).interior.values

        np.testing.assert_array_equal(
            np.asarray(via_moment), np.asarray(via_angular),
            err_msg="the moment-domain sibling's typed route must reproduce the "
            "angular binding's fused R∘Λ∘M chain BIT-IDENTICALLY (0 ULP) — the "
            "two bodies share Λ and the frame's R, and the ℓ=0 halves are the "
            "same scalar flux (pinned separately below). A mismatch means the "
            "moment end's typed grid path drifted from the validated aniso "
            "realization.",
        )

    def test_kernel_property_agrees_at_the_operator_ends(self, solver_p1_het):
        r"""``S.kernel.apply(ψ) == S_w.kernel.apply(M·ψ)`` — the ℓ≥1 halves alone.

        The row above compares the WHOLE emission; this one isolates the
        :math:`\ell\ge1` redistribution, which is where the two spellings
        actually differ (``frame.conjugate(Λ)`` = :math:`R\circ(\Lambda\circ M)`
        on the angular end, ``frame.reconstruct_after(Λ)`` = :math:`R\circ\Lambda`
        on the moment end — :math:`M` is skipped there because the operand is
        already projected; re-projecting would double-project). Fed the
        corresponding operands the two must coincide.
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)
        kernel = require_scattering_kernel_property(op)
        moments = op.flux_analysis.apply(psi)

        np.testing.assert_array_equal(
            np.asarray(op.on_moment_domain().kernel.apply(moments.values)),
            np.asarray(kernel.apply(psi.values)),
            err_msg="the moment end's R∘Λ fed M·ψ must equal the angular end's "
            "R∘Λ∘M fed ψ (0 ULP) — same Λ, same R, and M is exactly the "
            "projection the moment operand already carries.",
        )

    def test_the_two_ends_read_the_same_scalar_flux(self, solver_p1_het):
        r"""The DECOMPOSITION leg: the ℓ=0 halves are the same :math:`\phi`.

        Without this row the full-action equality above is a claim about a
        SUM and could hold by two errors cancelling. The angular end reduces
        :math:`\phi = \int\psi\,d\Omega` (``AngularFlux.integrate_angular``);
        the moment end reads the :math:`\ell = 0` slot of :math:`M\psi`
        (:math:`Y_0^0 = 1`). They are different reductions of the same data,
        and `[M]` they agree bit-for-bit — which is what makes the row above
        an :math:`\ell\ge1` statement.
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)
        moments = op.flux_analysis.apply(psi)
        scalar_space = op.on_moment_domain().isotropic_energy.domain
        np.testing.assert_array_equal(
            np.asarray(moments.scalar_flux(space=scalar_space).values),
            np.asarray(psi.integrate_angular().values),
            err_msg="the moment end's ℓ=0 slot must BE the angular end's "
            "∫ψ dΩ — if these differ the full-action equality is a "
            "coincidence of two compensating errors.",
        )

    def test_kernel_apply_shape_matches_per_ordinate_source(self, solver_p1_het):
        """The kernel output is per-ordinate ``(N, ng, nx, ny)``, pre-``1/W``.

        Confirms the ``kernel`` reproduces the ``R∘Λ∘M`` chain's shape (a
        per-ordinate angular field), NOT a moment tensor — the kernel is the
        moment→source composition, and the ``/W`` producer-side normalisation
        lives OUTSIDE it (lesson L18). A shape mismatch signals the kernel
        composed the wrong factor order.
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)
        kernel = require_scattering_kernel_property(op)
        out = np.asarray(kernel.apply(psi.values))
        N = solver_p1_het.quad.N
        ng = solver_p1_het.ng
        nx, ny = solver_p1_het.sn_mesh.spatial_shape
        require(
            out.shape == (N, ng, nx, ny),
            f"S.kernel.apply output must be per-ordinate (N, ng, nx, ny) = "
            f"{(N, ng, nx, ny)} (the R∘Λ∘M moment→source field, pre-1/W); "
            f"got {out.shape}.",
        )

    def test_aniso_path_is_non_degenerate(self, solver_p1_het):
        """Non-degeneracy guard: ℓ≥1 moments genuinely carry signal.

        If ``scattering_order < 1`` or the ℓ≥1 moments were all zero, the
        cross-check would collapse to comparing two zero arrays (vacuous).
        This row asserts the config exercises the aniso path so the
        equivalence gates have teeth — and, since step 5, that the ℓ≥1
        EMISSION itself is non-zero on both ends (a zero morphism satisfies
        every equality row above with both sides structurally zero).
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)
        moments = op.flux_analysis.apply(psi)
        require(
            op.legendre_order >= 1,
            "Cross-check config must be ANISOTROPIC (scattering_order ≥ 1) "
            "so Λ's ℓ≥1 blocks are exercised — else the kernel cross-check "
            "is vacuous.",
        )
        require(
            bool(np.any(np.asarray(moments.values)[1:] != 0.0)),
            "ℓ≥1 moments must be non-zero (the flux must be non-isotropic) "
            "so the R∘Λ∘M chain genuinely runs the aniso reconstruction. "
            "They were all zero — the cross-check would compare two zeros.",
        )
        for name, emitted in (
            ("angular", np.asarray(op.kernel.apply(psi.values))),
            (
                "moment",
                np.asarray(op.on_moment_domain().kernel.apply(moments.values)),
            ),
        ):
            require(
                float(np.max(np.abs(emitted))) > 1e-6,
                f"the {name} end's ℓ≥1 emission is ≈ 0 — the equality rows "
                f"above would hold with both sides structurally zero.",
            )


# ═══════════════════════════════════════════════════════════════════════
# C′ — the kernel is the ℓ≥1 aniso redistribution ONLY, NOT the full S.
#      Pins the documented partial-kernel semantics (qa-flagged blind spot).
# ═══════════════════════════════════════════════════════════════════════


class TestScatteringKernelIsAnisoSubcomponent:
    """``S.kernel`` is the nonlocal ℓ≥1 redistribution — NOT the full S.

    The §5.6 ``kernel`` exposes the genuinely-nonlocal-in-angle ``R∘Λ∘M``
    redistribution (ℓ≥1, ``skip_l0=True``, pre-``1/W``). The full
    :meth:`ScatteringOperator.apply` ALSO adds the isotropic ℓ=0 P0
    in-scatter, the (n,2n) up-scatter, and the producer-side ``1/W``
    normalisation. A future consumer who reads ``kernel`` as the COMPLETE
    scattering operator would silently lose those (~95% of the source on
    a typical config). This gate pins the documented partial-kernel
    semantics so that misreading reddens here — and complements the
    fission asymmetry (``FissionOperator.kernel`` IS the full F, but
    ``ScatteringOperator.kernel`` is a strict sub-component).
    """

    def test_kernel_is_strict_subcomponent_of_full_apply(self, solver_p1_het):
        """``S.kernel.apply(ψ) != S.apply(ψ)`` — the kernel omits P0/n2n/(1/W).

        Not an FP-level difference: the kernel (ℓ≥1 aniso, pre-``1/W``)
        and the full per-ordinate source differ substantially. ``require``
        (Mode-8 ``-O``-safe) asserts they do NOT coincide.
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)
        kernel = require_scattering_kernel_property(op)

        aniso_pre_w = np.asarray(kernel.apply(psi.values))   # ℓ≥1, pre-1/W
        # CS4c step 5: the gain is composite-bound; the bulk emission rides a
        # zero-trace composite.
        full = np.asarray(                                   # iso+aniso, /W
            bulk_apply(op, psi).values,
        )

        require(
            aniso_pre_w.shape == full.shape,
            f"shapes must match for the sub-component comparison; got "
            f"{aniso_pre_w.shape} vs {full.shape}.",
        )
        require(
            not np.allclose(aniso_pre_w, full),
            "S.kernel.apply (the ℓ≥1 aniso redistribution, pre-1/W) must NOT "
            "equal the full S.apply (which adds the ℓ=0 P0 iso in-scatter, the "
            "(n,2n) up-scatter, and the 1/W normalisation). They coincided — "
            "a consumer reading `kernel` as the COMPLETE scattering operator "
            "would silently drop those contributions. The kernel is the §5.6 "
            "nonlocal-in-angle sub-component, NOT the full operator.",
        )

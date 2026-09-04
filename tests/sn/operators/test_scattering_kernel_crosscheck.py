r"""#257 S6 — Spec C: the scattering ``kernel`` reproduces the aniso ``R∘Λ∘M`` path.

S6 adds a ``kernel`` property to :class:`ScatteringOperator` → a typed
``OperatorProduct(R, OperatorProduct(Λ, M))`` reproducing the EXISTING
``_aniso_source_from_moment_values`` chain ``R(Λ(M·ψ))`` (M = the
:attr:`frame`'s analysis face, R = its reconstruction face, Λ =
``LegendreMomentScattering``; ``OperatorProduct.apply(x) = a.apply(b.apply(x))``,
``operator.py:826``). With ``kernel`` the operator becomes an
``IntegralKernelOperator``. The 5 ``@singledispatchmethod`` arms are
UNCHANGED.

The new ``kernel`` reproduces the EXISTING realization byte-for-byte —
same einsums, same order (``M`` then ``Λ`` then ``R``, the exact chain
``build_aniso_source`` runs after its own ``M`` projection). So:

    S.kernel.apply(ψ.values)  ==  _aniso_source_from_moment_values(M.apply(ψ.values))   (0 ULP)

where ``M = S.frame.analysis`` (the SAME analysis face the kernel composes)
and ``_aniso_source_from_moment_values`` is the SHARED ``R∘Λ_{ℓ≥1}`` map
(``skip_l0=True``).

CRITICAL DEMARCATION (L11): this is the EQUIVALENCE / de-risk leg — the
new ``kernel`` reproduces the EXISTING aniso realization. It is NOT a
correctness claim about the scattering *physics*. The
structurally-independent reference for the scattering physics is the
EXISTING aniso MMS gate
(``tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py``,
both l0 ``verifies("pn-scatter", "flux-moments")`` cases) — that gate
stays green and is the L1 backing; this file only proves the new typed
``kernel`` composes the same numpy chain the validated path already runs.

The config is ANISOTROPIC (≥P1) + heterogeneous + multi-group so Λ's
ℓ≥1 blocks are genuinely exercised (a P0-only config would null the
aniso path and make the cross-check vacuous — guarded by the
non-degeneracy assert).

vv Mode-11: the cross-check reads ``S.kernel`` OFF the live operator (the
NEW property) and compares against the live ``_aniso_source_from_moment_values``
+ ``S.frame.analysis`` chain — mutating the kernel property reddens the
gate. It does NOT route around the new property.

vv Mode-8: every gate uses ``np.testing.*`` / ``require`` (function
calls, fire under ``python -O``) — NEVER a bare ``assert``.

``foundation`` — software invariants on the operator's type surface plus
a bit-identity equivalence check. No theory-page ``:label:`` is claimed
(the physics ``:label:`` s ``pn-scatter`` / ``flux-moments`` are pinned by
the MMS gate, NOT here).
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


class TestScatteringKernelReproducesAnisoPath:
    def test_kernel_apply_equals_existing_R_Lambda_M_chain(self, solver_p1_het):
        """``S.kernel.apply(ψ.values) == _aniso_source_from_moment_values(M·ψ)`` (0 ULP).

        EQUIVALENCE / de-risk (NOT correctness — the physics reference is
        the aniso MMS gate, see module docstring). The new typed
        ``kernel`` (``OperatorProduct(R, OperatorProduct(Λ, M))``)
        composes the SAME einsums in the SAME order the existing
        ``build_aniso_source`` runs after its own ``M`` projection.
        Mode-11: reads ``S.kernel`` OFF the live operator.
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)

        kernel = require_scattering_kernel_property(op)  # NEW S6 member
        via_kernel = kernel.apply(psi.values)

        # The EXISTING chain: M projection (the frame's analysis face — the
        # SAME object the kernel composes) then the shared R∘Λ_{ℓ≥1} map.
        moments = op.frame.analysis.apply(psi.values)
        via_existing = op._aniso_source_from_moment_values(moments)

        np.testing.assert_array_equal(
            via_kernel, via_existing,
            err_msg="S.kernel.apply(ψ.values) must reproduce the EXISTING "
            "_aniso_source_from_moment_values(M·ψ) chain BIT-IDENTICALLY "
            "(0 ULP) — the typed kernel composes the same R∘Λ∘M einsums in "
            "the same order. A mismatch means the new kernel does not "
            "reproduce the validated aniso realization.",
        )

    def test_kernel_apply_shape_matches_per_ordinate_source(self, solver_p1_het):
        """The kernel output is per-ordinate ``(N, ng, nx, ny)``, pre-``1/W``.

        Confirms the new ``kernel`` reproduces the ``R∘Λ∘M`` chain's shape
        (a per-ordinate angular field), NOT a moment tensor — the kernel
        is the moment→source composition, and the ``/W`` producer-side
        normalisation lives OUTSIDE it (lesson L18). A shape mismatch
        signals the kernel composed the wrong factor order.
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
        equivalence gate has teeth.
        """
        op = solver_p1_het.scattering_op
        psi = _aniso_psi(solver_p1_het)
        moments = op.frame.analysis.apply(psi.values)
        require(
            op.scattering_order >= 1,
            "Cross-check config must be ANISOTROPIC (scattering_order ≥ 1) "
            "so Λ's ℓ≥1 blocks are exercised — else the kernel cross-check "
            "is vacuous.",
        )
        require(
            bool(np.any(moments[1:] != 0.0)),
            "ℓ≥1 moments must be non-zero (the flux must be non-isotropic) "
            "so the R∘Λ∘M chain genuinely runs the aniso reconstruction. "
            "They were all zero — the cross-check would compare two zeros.",
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
        full = np.asarray(op.apply(psi).values)              # iso+aniso+n2n, /W

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

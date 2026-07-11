r"""The ψ½ Ray-Characteristic coupled system — regression floor + per-step invariants.

**Campaign home.** This module is the verification substrate for the
*coupled block operator* campaign (`.claude/plans/coupled_block_operator_campaign.md`):
the within-group augmented SN system posed as a 2×2 coupled block operator

.. code-block:: text

    [ A_AA   A_AB ] [ transport ]      A_AA = L + C − S − B     (System A)
    [ A_BA   A_BB ] [ ray       ]      A_BB = RayOp             (System B: the ψ½ ray)

over two systems — System A (the transport bulk⊕trace `FullField`) and System B
(the ψ½ **Radial-Characteristic** starting-direction flux, a codim-1 `FaceField`
carrying two boundary conditions: r=R Dirichlet + r=0 pole-reflection). The
seed IS System B's boundary condition; `A_BA` (bulk→ray) is the Schur-fold
coupling currently welded un-named into the S/F scattering/fission arms — a
coupling *block* never posed in the algebra. The campaign names all four blocks
as operators; "where the block goes" (folded into the resolvent, an explicit
block operator, or a DSA preconditioner) then becomes a *composition choice* the
machinery supports — that flexibility IS operator algebra realized.

**:class:`TestRegressionFloor` — landed FIRST, before any production change**, so
every subsequent campaign step (System B posing, the `A_BA`/`A_AB` un-welds, the
`CoupledOperator` assembly, the coupled solve) diffs against a *pinned* baseline
of the measured block structure. Promoted verbatim from the numerics-investigator
design diagnostics ``derivations/diagnostics/diag_coupled_0{1,2}_*.py``
(2026-07-07). The six pins and their measured values on the seed-carrying
vacuum sphere (GL S4, nx=5, 2G, c=0.4):

======================================================  ==================================
pin                                                      measured (the baseline this floor holds)
======================================================  ==================================
loss ``(L+C)`` is block-TRIANGULAR in the ray            ``A_sb = A_st = 0`` exact; ``A_bs = 7.505``, ``A_ss = 5.000``
bulk→ray ``A_BA`` lives in the LAGGED scattering gain     ``S_sb = 0.183``; ``S_bs = 0`` exact (ψ½ zero moment weight)
outer-SI splitting rate is bounded by ``c``               ``ρ(M⁻¹N) = 0.371 < c = 0.4``
the folded ray seed is a DIRECT (nilpotent) solve         ``ρ(lag) = 0`` — no bulk→seed back-edge
the welded sweep is the EXACT inverse of ``(L+C)``        ``‖solve(apply(ψ)) − ψ‖ = 3.5e-16``
extraction is PRINCIPLED-equivalent, not bit-identical    ``‖welded − dense_LU‖ = 5.5e-16`` (distinct reduction trees)
======================================================  ==================================

The last row is the oracle every EXTRACT step of the campaign pins against
(``coding-elegance``/``vv-principles`` §"Bit-identity vs principled-equivalence":
principled-equivalence + the invariant test is the bar, never byte-identity).

**Runtime discipline (vv Mode 8).** Every gate raises via :func:`pytest.fail`
(a function call — fires under the canonical ``python -O``), never a bare
``assert`` (which ``-O`` strips outside pytest-rewritten test bodies). The
``-s`` print lines echo the measured values for the design record.

References: GH #284 (the triangular sweep = forward substitution), #282
(route-a direct ψ½ seed), #280 (the walk unification).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from numpy.typing import NDArray

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import (
    RadialCharacteristicBoundaryOperator,
    SNBoundaryOperator,
    _RayBoundaryFullFieldGain,
)
from orpheus.sn.operators.streaming import StreamingOperator
import orpheus.sn.loss_representation as _lr_mod
import orpheus.sn.operators.radial_characteristic as _rc_mod
from orpheus.sn.operators.radial_characteristic import (
    RadialCharacteristicEmission,
    RadialCharacteristicOperator,
    RadialCharacteristicReconstruction,
    RadialCharacteristicSeeding,
    _RayEmissionFullFieldGain,
)
import orpheus.numerics.spaces.radial_characteristic_space as _rcs_mod
# Campaign step 4c (THE LIFT): the Fold ``RadialCharacteristicReconstruction``
# migrated transport → sn (it is a factor of the sn coupling operator
# ``RadialCharacteristicEmission``); ``_rcr_mod`` is the reconstruction/fold home
# (now the same module as ``_rc_mod``) — kept as a distinct alias for the spies
# that patch the fold as the reconstruction sees it.
import orpheus.sn.operators.radial_characteristic as _rcr_mod
import orpheus.sn.solver as _solver_mod
from orpheus.sn.solver import SNSolver, _lagged_gains, _within_group_triple
from orpheus.sn.spatial.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
    carlson_inward_sweep_transpose,
    radial_characteristic_forward_residual,
)
from orpheus.numerics.operator import SystemRole, _join_system_roles
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.full_field import FullField
from orpheus.transport.radial_characteristic_composite import (
    RadialCharacteristicComposite,
)
from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)
from orpheus.transport.source_sinks.angular_source_sink import AngularSourceSink
from orpheus.transport.source_sinks.radial_characteristic_source_sink import (
    RadialCharacteristicSourceSink,
)
from orpheus.derivations.common.xs_library import make_mixture

pytestmark = pytest.mark.foundation


# ── compact dense-probe helpers (mirror test_radial_characteristic_metric) ──


def _mixture(sig_t: float, sig_s: float, ng: int):
    """A group-graded diagonal-scatter mixture (asymmetric in g — vv L2/anti-#2)."""
    st = np.array([sig_t * (1.0 + 0.4 * g) for g in range(ng)])
    ss = np.diag([sig_s * (1.0 + 0.4 * g) for g in range(ng)])
    return make_mixture(sig_t=st, sig_c=st - ss.sum(axis=0), sig_f=np.zeros(ng),
                        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=ss)


def _sphere(nx: int = 5, ng: int = 2, sigma: float = 1.0, c: float = 0.4,
            bc: str = "vacuum"):
    """Seed-carrying sphere (GL S4) with scattering ratio ``c`` and outer law ``bc``.

    Default vacuum (the regression floor). ``bc="reflective"`` gives a
    seed-carrying sphere whose outer face drives ``B_b``'s non-trivial specular
    corner swap — the vacuum floor never exercises that arm (``_reflect_corner``
    returns zeros for vacuum), so the ``B_b`` gates below need it.
    """
    mesh = Mesh1D(edges=np.linspace(0.0, 4.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
                  coord=CoordSystem.SPHERICAL, bc_right=BC(bc))
    return SNMesh(mesh, Quadrature.gauss_legendre(4), {0: _mixture(sigma, c * sigma, ng)})


def _loss(sn, slope: float = 0.4):
    r"""The within-group loss operator ``L + C`` with per-group ``σ_t = 1 + slope·g``.

    ``slope`` parametrizes the two diagnostics' constructions into one helper:
    the block-structure pins (triangularity, nilpotent lag) used ``slope = 0.4``
    (matching the sphere mixture's own ``σ_t``); the round-trip pins (welded =
    exact inverse, extract = principled-equiv) used ``slope = 0.3`` (a distinct
    invertible ``σ_t`` — the round-trip holds for *any* invertible ``(L+C)``).
    """
    sig_t = np.stack(
        [np.full(sn.spatial_shape, 1.0 + slope * g) for g in range(sn.ng)], 0)
    return StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)


def _template(sn):
    """A zero :class:`FullField` (bulk ⊕ trace ⊕ ψ½ seed) to seed dense probes."""
    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    n_tr = int(sn.angular_trace.layout.total_size)
    n_sd = sn.radial_characteristic_space.shape[0]
    return FullField(
        interior=AngularFlux.from_mesh(np.zeros((N, ng, nx)), sn),
        boundary=AngularBoundaryFlux(values=np.zeros(n_tr), space=sn.angular_trace, mesh=sn),
        radial_characteristic=RadialCharacteristicFlux(
            values=np.zeros(n_sd), space=sn.radial_characteristic_space, mesh=sn))


def _dense(fn, tpl):
    """Materialise ``fn`` as a dense matrix by probing the ``FullField`` basis."""
    n = tpl.to_flat().size
    M = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1.0
        M[:, j] = fn(FullField.from_flat(e, tpl)).to_flat()
    return M


def _blocks(sn):
    """The (bulk, trace, seed) row/col slices of the flat ``FullField`` layout."""
    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    nb = N * ng * nx
    nt = int(sn.angular_trace.layout.total_size)
    ns = sn.radial_characteristic_space.shape[0]
    return slice(0, nb), slice(nb, nb + nt), slice(nb + nt, nb + nt + ns), (nb, nt, ns)


def _bn(M, r, c):
    """The max-abs of block ``M[r, c]`` (a block-norm for triangularity probes)."""
    return float(np.max(np.abs(M[r, c])))


class TestRegressionFloor:
    r"""The pinned baseline of the ψ½ coupled-block structure — landed BEFORE any
    production change so every campaign step diffs against it.

    Each test carries a mutation tooth in its docstring: the specific structural
    corruption (a bulk→seed back-edge in the sweep, a lost scattering source, a
    lagged direct seed) that moves the measured block-norm O(1) and reddens the
    gate. Promoted from ``diag_coupled_01_psi_half_block_structure`` (pins 1–4)
    and ``diag_coupled_02_wrap_vs_extract`` (pins 5–6).
    """

    def test_loss_operator_is_block_triangular_in_the_ray(self):
        r"""Within ``(L+C)`` the seed rows are self-contained (``A_sb = A_st = 0``)
        while the seed feeds the bulk (``A_bs ≠ 0``) — the #284 direct-solve
        certificate. Mutation tooth: a bulk→seed coupling leaking into ``(L+C)``
        (e.g. the ray reading a bulk moment in-sweep) makes ``A_sb`` jump O(1)."""
        sn = _sphere()
        LC = _loss(sn, slope=0.4)
        b, t, s, _ = _blocks(sn)
        A = _dense(LC.apply, _template(sn))
        a_sb, a_st = _bn(A, s, b), _bn(A, s, t)
        a_bs, a_ss = _bn(A, b, s), _bn(A, s, s)
        print(f"  (L+C): A_sb={a_sb:.2e} A_st={a_st:.2e} | A_bs={a_bs:.3f} A_ss={a_ss:.3f}")
        a_ba = max(a_sb, a_st)
        if not (a_ba < 1e-12):
            pytest.fail(f"A_BA(within L+C)={a_ba:.3e} ≠ 0 — the ray is NOT block-triangular; "
                        f"the sweep is no longer exact forward substitution (#284 broken).")
        if not (a_bs > 1.0 and a_ss > 1.0):
            pytest.fail(f"seed→bulk feed A_bs={a_bs:.3e} or self-block A_ss={a_ss:.3e} vanished "
                        f"— the ψ½ coupling is degenerate (route-a seed not wired).")

    def test_bulk_to_ray_coupling_lives_in_the_lagged_A_BA_gain(self):
        r"""The full ``A_BA`` (bulk→ray) is carried by the LIFTED
        :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicEmission`
        gain (``A_BA_sb ≠ 0``, the ray scattering source) — NOT by the
        model-generic ``S``, which is now PURE BULK (``S_sb = 0``, campaign step
        4c THE LIFT re-pointed this pin from S to A_BA). The ray still cannot feed
        the scalar flux (``S_bs = 0``, ψ½ has zero moment weight). So the coupled
        system's off-diagonal ``A_BA`` is an OUTER-iterated block riding its OWN
        lagged gain, never a within-sweep one.

        Mutation tooth: an ``A_BA`` that emitted into the bulk (a double-count with
        S's bulk) would make ``A_BA_bb ≠ 0``; an S that re-grew the ray arm would
        make ``S_sb ≠ 0`` (the pre-lift regression this flip retired).

        B.2b re-point: the FullField-densifiable face is the TRANSIENT gain
        adapter (the block itself is FullField → RadialCharacteristicComposite —
        un-densifiable on one template); on the BLOCK the ``A_BA_bb = 0``
        disjointness is STRUCTURAL (no bulk slot in the codomain), so the
        measured ``aba_bb`` row now pins the ADAPTER's embed."""
        sn = _sphere()
        _, S, _ = _within_group_triple(SNSolver(sn))
        A_BA = _RayEmissionFullFieldGain(
            RadialCharacteristicEmission(sn, S.isotropic_kernel))
        b, _, s, _ = _blocks(sn)
        tpl = _template(sn)
        Sd = _dense(S.apply, tpl)
        Ad = _dense(A_BA.apply, tpl)
        s_sb, s_bs = _bn(Sd, s, b), _bn(Sd, b, s)
        aba_sb, aba_bb = _bn(Ad, s, b), _bn(Ad, b, b)
        print(f"  S(pure bulk): S_sb={s_sb:.2e} S_bs={s_bs:.2e} | "
              f"A_BA: A_BA_sb(bulk→ray)={aba_sb:.3f} A_BA_bb(no bulk-out)={aba_bb:.2e}")
        # The coupling lives in A_BA, not S.
        if not (aba_sb > 1e-6):
            pytest.fail(f"A_BA_sb={aba_sb:.3e} ≈ 0 — the lifted A_BA does NOT source the ray; "
                        f"the bulk→ray coupling is missing (the ψ½ gain is unwired).")
        # S is now pure bulk (the LIFT dropped the ray side-channel).
        if not (s_sb < 1e-12):
            pytest.fail(f"S_sb={s_sb:.3e} ≠ 0 — the model-generic scattering gain still sources "
                        f"the ray; the LIFT did not make S pure bulk (a re-grown ray arm).")
        # A_BA is disjoint from the bulk (no double-count with S's bulk).
        if not (aba_bb < 1e-12):
            pytest.fail(f"A_BA_bb={aba_bb:.3e} ≠ 0 — A_BA emits into the bulk (a double-count "
                        f"with S's bulk); it must write ONLY the ray (disjoint direct sum).")
        # The ray still cannot feed the scalar flux (ψ½ zero moment weight).
        if not (s_bs < 1e-12):
            pytest.fail(f"S_bs={s_bs:.3e} ≠ 0 — the ray leaks into the scalar flux; ψ½ must have "
                        f"ZERO moment weight (the ghost-metric physics is violated).")

    def test_outer_si_splitting_rate_is_bounded_by_scattering_ratio(self):
        r"""``A = (L+C) − S − B = M − N`` with ``M = (L+C)``, ``N = S + B``: the
        outer source iteration's spectral radius ``ρ(M⁻¹N) ≤ c`` (Adams–Larsen;
        strictly below c for vacuum leakage). This is the convergence rate of the
        ONE genuinely-iterated coupling in the ψ½ system."""
        c = 0.4
        sn = _sphere(c=c)
        LC, S, B = _within_group_triple(SNSolver(sn))
        tpl = _template(sn)
        M = _dense(LC.apply, tpl)
        N = _dense(S.apply, tpl) + _dense(B.apply, tpl)
        rho = float(np.max(np.abs(np.linalg.eigvals(np.linalg.solve(M, N)))))
        print(f"  ρ(M⁻¹N) = {rho:.4f}   (c={c}; below c for vacuum leakage)")
        if not (0.0 < rho < c + 1e-6):
            pytest.fail(f"ρ(M⁻¹N)={rho:.4f} not in (0, c={c}] — the within-group SI splitting "
                        f"rate is not bounded by the scattering ratio (splitting mis-posed).")

    def test_folded_ray_seed_is_a_direct_solve_zero_spectral_radius(self):
        r"""FOLDED (route a): the ray→bulk seed is solved in-sweep (forward
        substitution). Lagging it instead (moving ``A_bs`` to the iteration's N
        side) yields a NILPOTENT iteration matrix (``ρ = 0``, converges in 2 steps)
        — because ``A_sb = 0`` there is no back-edge. So the fold buys a direct
        ρ=0 solve; lagging the *direct* seed is merely wasteful (2 passes for the
        same answer), NOT unstable. (The historical #282 EDGE-EXTRAPOLATION seed —
        ψ½ = E·ψ_bulk, a genuine cycle — diverges ρ≈70; that is why route-a was
        needed. Documented, not reconstructed here.)"""
        sn = _sphere()
        LC = _loss(sn, slope=0.4)
        _, _, _, (nb, nt, ns) = _blocks(sn)
        M = _dense(LC.apply, _template(sn))
        ab = np.r_[np.arange(nb), np.arange(nb, nb + nt)]      # bulk+trace = System A
        sd = np.arange(nb + nt, nb + nt + ns)                  # seed = System B
        perm = np.r_[ab, sd]
        Mp = M[np.ix_(perm, perm)]
        nA = len(ab)
        M_lag = Mp.copy(); M_lag[:nA, nA:] = 0.0               # lag the seed→bulk feed
        N_lag = np.zeros_like(Mp); N_lag[:nA, nA:] = Mp[:nA, nA:]
        rho_lag = float(np.max(np.abs(np.linalg.eigvals(np.linalg.solve(M_lag, N_lag)))))
        print(f"  ρ(lag route-a's direct seed) = {rho_lag:.3e}  (nilpotent → 2 steps)")
        if not (rho_lag < 1e-10):
            pytest.fail(f"ρ_lag={rho_lag:.3e} ≠ 0 — a bulk→seed back-edge appeared (A_sb≠0), so "
                        f"lagging the seed is no longer nilpotent (the triangular structure broke).")

    def test_welded_sweep_is_exact_direct_inverse(self):
        r"""``(L+C).solve((L+C).apply(ψ)) ≈ ψ`` at machine precision — the welded
        sphere sweep IS the exact direct inverse (route-a). A WRAP inherits this
        bit-for-bit. Mutation tooth: a re-introduced ψ½ seed lag makes the sphere
        round-trip blow up (pre-route-a residual was O(1e5))."""
        sn = _sphere()
        LC = _loss(sn, slope=0.3)
        tpl = _template(sn)
        nb = sn.quad.N * sn.ng * sn.nx
        rng = np.random.default_rng(7)
        psi0 = np.zeros(tpl.to_flat().size)
        psi0[:nb] = rng.standard_normal(nb)                 # random physical bulk
        psi0_ff = FullField.from_flat(psi0, tpl)
        back = LC.solve(LC.apply(psi0_ff)).to_flat()
        rel = np.max(np.abs(back[:nb] - psi0[:nb])) / (np.max(np.abs(psi0[:nb])) + 1e-300)
        print(f"  ||solve(apply(ψ)) − ψ||_bulk_rel = {rel:.3e}")
        if not (rel < 1e-12):
            pytest.fail(f"round-trip {rel:.3e} — the welded sweep is NOT the exact inverse of "
                        f"(L+C); the direct ψ½ solve regressed (a WRAP would inherit the error).")

    def test_extract_to_dense_is_principled_equivalent_not_bit_identical(self):
        r"""The welded WDD sweep and a LAPACK LU of the SAME assembled ``(L+C)``
        agree to a few ULP (~1e-15) — principled-equivalent, different reduction
        trees. This is the numerical cost of EXTRACTION: the answer is preserved to
        machine precision, but bit-identity is lost. WRAP (same code) keeps
        bit-identity; EXTRACT trades it for ~1e-15 drift. This row is the oracle
        every EXTRACT step of the campaign pins against."""
        sn = _sphere()
        LC = _loss(sn, slope=0.3)
        tpl = _template(sn)
        nb = sn.quad.N * sn.ng * sn.nx
        M = _dense(LC.apply, tpl)                            # the EXTRACTED explicit matrix
        rng = np.random.default_rng(11)
        psi0 = np.zeros(tpl.to_flat().size)
        psi0[:nb] = rng.standard_normal(nb)
        q = LC.apply(FullField.from_flat(psi0, tpl))
        psi_weld = LC.solve(q).to_flat()                    # welded sweep (forward substitution)
        psi_dense = np.linalg.solve(M, q.to_flat())         # extracted dense LU
        diff = np.max(np.abs(psi_weld[:nb] - psi_dense[:nb])) / (np.max(np.abs(psi_dense[:nb])) + 1e-300)
        print(f"  ||welded − dense_LU||_bulk_rel = {diff:.3e}  (principled ~1e-15, not 0)")
        if not (diff < 1e-11):
            pytest.fail(f"welded vs dense LU differ by {diff:.3e} — an EXTRACTED block solve does "
                        f"not even reach principled-equivalence; the extraction dropped the row "
                        f"contract of the sweep (naive dense M⁻¹ ignores inflow/seed rows).")


# ── Step-1b helpers: the boundary un-weld B = B_a + B_b ───────────────────


def _seed_composite(sn, seed_values: NDArray) -> FullField:
    """A composite with zero bulk, the given trace, and the given ψ½ seed."""
    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    n_tr = int(sn.angular_trace.layout.total_size)
    return FullField(
        interior=AngularFlux.from_mesh(np.zeros((N, ng, nx)), sn),
        boundary=AngularBoundaryFlux(values=np.zeros(n_tr), space=sn.angular_trace, mesh=sn),
        radial_characteristic=RadialCharacteristicFlux(
            values=seed_values, space=sn.radial_characteristic_space, mesh=sn),
    )


def _random_composite(sn, rng) -> FullField:
    """A composite with random bulk, trace, and ψ½ seed (all blocks non-zero)."""
    N, nx, ng = sn.quad.N, sn.nx, sn.ng
    n_tr = int(sn.angular_trace.layout.total_size)
    ns = sn.radial_characteristic_space.shape[0]
    return FullField(
        interior=AngularFlux.from_mesh(rng.standard_normal((N, ng, nx)), sn),
        boundary=AngularBoundaryFlux(
            values=rng.standard_normal(n_tr), space=sn.angular_trace, mesh=sn),
        radial_characteristic=RadialCharacteristicFlux(
            values=rng.standard_normal(ns), space=sn.radial_characteristic_space, mesh=sn),
    )


def _ray_composite(sn, seed_values: NDArray) -> RadialCharacteristicComposite:
    """System B's own carrier with the unified seed layout split onto its members
    (the B.2b block probes — the bridge round-trip is the exact re-labeling)."""
    return RadialCharacteristicComposite.from_unified(
        RadialCharacteristicFlux(
            values=np.asarray(seed_values, dtype=float),
            space=sn.radial_characteristic_space, mesh=sn))


def _dense_ray(fn, sn) -> NDArray:
    """Densify a System-B block (composite → composite) by probing the composite
    basis in the UNIFIED layout — bit-comparable with the pre-B.2b ``_dense_seed``
    matrices (the role-preserving bridge is a pure gather/scatter)."""
    ns = sn.radial_characteristic_space.shape[0]
    M = np.zeros((ns, ns))
    for j in range(ns):
        e = np.zeros(ns)
        e[j] = 1.0
        M[:, j] = fn(_ray_composite(sn, e)).to_unified().values
    return M


def _v_cell_seed(sn) -> NDArray:
    """The ``G_sd = V_cell`` seed metric (production ``inner_product_weights``)."""
    return np.asarray(
        sn.radial_characteristic_space.inner_product_weights, dtype=float)


def _g_recip(fwd: NDArray, T: NDArray, g: NDArray, rng) -> float:
    """The metric-reciprocity defect ``|⟨fwd x, y⟩_g − ⟨x, T y⟩_g| / norm``."""
    x = rng.standard_normal(g.size)
    y = rng.standard_normal(g.size)
    lhs = float((fwd @ x) @ (g * y))
    rhs = float(x @ (g * (T @ y)))
    return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-30)


class TestBoundaryUnweld:
    r"""``B = B_a + B_b`` — the boundary un-weld is a DISJOINT direct sum.

    ``B_a`` (:class:`SNBoundaryOperator`, System A's trace boundary) touches only
    the trace and emits a **present-zero** ray block; ``B_b``
    (:class:`RadialCharacteristicBoundaryOperator`, System B's ray corner) touches
    only the ray and emits a present-zero trace. Their sum reconstructs the whole
    augmented boundary bit-identically (RULING P1). The reflective sphere is used
    so BOTH arms are non-trivial (vacuum would leave B_b zero).
    """

    def test_b_a_touches_only_the_trace_present_zero_ray(self):
        r"""``B_a`` reflects the trace but emits a PRESENT-ZERO ray (not ``None``
        — the mixed-presence law would raise under ``B_a + B_b``). Mutation tooth:
        a ``B_a`` that re-grew the ray arm would emit a non-zero ray block."""
        sn = _sphere(bc="reflective")
        out = SNBoundaryOperator(sn).apply(_random_composite(sn, np.random.default_rng(1)))
        if out.radial_characteristic is None:
            pytest.fail("B_a emitted a None ray on a seed-carrying composite — the "
                        "mixed-presence law will raise under B_a + B_b (must be present-zero).")
        np.testing.assert_array_equal(
            out.radial_characteristic.values, 0.0,
            err_msg="B_a emitted a NON-ZERO ray block — it re-grew the ray arm "
                    "that belongs to B_b (the un-weld leaked).")
        if not np.max(np.abs(out.boundary.values)) > 0.0:
            pytest.fail("B_a emitted a zero trace on a reflective sphere — it is not "
                        "reflecting the trace (System A boundary is dead).")

    def test_b_b_touches_only_the_ray_present_zero_trace(self):
        r"""Through the TRANSIENT gain adapter, ``B_b`` emits a PRESENT-ZERO
        trace and zero bulk (the pre-B.2b byte shape the ``B_a + B_b`` sum
        needs). On the BLOCK itself "B_b touches the trace" is UNSPELLABLE
        since the B.2b re-type (its codomain has no trace slot — Pattern 4);
        the container-type rows pin that. Mutation tooth: an adapter leaking a
        trace action emits a non-zero trace."""
        sn = _sphere(bc="reflective")
        psi = _random_composite(sn, np.random.default_rng(2))
        block = RadialCharacteristicBoundaryOperator(sn)
        out = _RayBoundaryFullFieldGain(block).apply(psi)
        np.testing.assert_array_equal(
            out.boundary.values, 0.0,
            err_msg="B_b's adapter emitted a NON-ZERO trace block — it leaked a "
                    "trace action that belongs to B_a.")
        if not np.max(np.abs(out.radial_characteristic.values)) > 0.0:
            pytest.fail("B_b emitted a zero ray corner on a reflective sphere — the "
                        "System B boundary arm is dead.")
        # The BLOCK's structural half: composite in/out, SOURCE members out.
        block_out = block.apply(
            RadialCharacteristicComposite.from_unified(psi.radial_characteristic))
        if type(block_out) is not RadialCharacteristicComposite:
            pytest.fail(f"B_b block returned {type(block_out).__name__}")
        if type(block_out.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"B_b interior member is {type(block_out.interior).__name__} "
                        f"— the block must emit the SOURCE pair.")
        if type(block_out.boundary) is not RadialCharacteristicBoundarySourceSink:
            pytest.fail(f"B_b boundary member is {type(block_out.boundary).__name__} "
                        f"— the block must emit the SOURCE pair.")
        # Adapter ≡ block, byte-for-byte on the ray (the embed is a re-label).
        np.testing.assert_array_equal(
            out.radial_characteristic.values, block_out.to_unified().values,
            err_msg="adapter ray ≠ block ray — the transient embed moved values.")

    def test_sum_reconstructs_both_blocks_disjointly(self):
        r"""``(B_a + B_b).apply`` = ``B_a``'s trace ⊕ ``B_b``'s ray, byte-for-byte
        per block — the disjoint direct sum reconstructs the whole boundary."""
        sn = _sphere(bc="reflective")
        psi = _random_composite(sn, np.random.default_rng(3))
        B_a = SNBoundaryOperator(sn)
        # B.2b: the production sum rides the transient adapter (the raw block
        # declares System B's composite spaces, which the OperatorSum guard
        # correctly refuses to sum with B_a's FullField spaces).
        B_b = _RayBoundaryFullFieldGain(RadialCharacteristicBoundaryOperator(sn))
        out = (B_a + B_b).apply(psi)
        out_a, out_b = B_a.apply(psi), B_b.apply(psi)
        # The composite trace is exactly B_a's (B_b contributes present-zero).
        np.testing.assert_array_equal(out.boundary.values, out_a.boundary.values)
        # The composite ray is exactly B_b's (B_a contributes present-zero).
        np.testing.assert_array_equal(
            out.radial_characteristic.values, out_b.radial_characteristic.values)

    def test_seedless_mesh_has_no_b_b_and_b_is_b_a_alone(self):
        r"""B.2b re-point of the old None-pass-through row: a seedless mesh has
        no System B, so ``B_b`` is UNCONSTRUCTABLE there (Pattern 4 — the old
        "B_b passes None through" behavior is now unspellable), and the
        production boundary is ``B_a`` ALONE (``_within_group_triple``'s
        presence branch). The ray of ``B_a``'s output stays ``None``."""
        slab = SNMesh(
            Mesh1D(edges=np.linspace(0.0, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CARTESIAN, bc_right=BC("reflective"),
                   bc_left=BC("reflective")),
            Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, 2)})
        with pytest.raises(ValueError, match="carries no ψ½ ray"):
            RadialCharacteristicBoundaryOperator(slab)
        # The production boundary on a seedless mesh is B_a alone.
        _, _, B = _within_group_triple(SNSolver(slab))
        if not isinstance(B, SNBoundaryOperator):
            pytest.fail(f"seedless _within_group_triple boundary is "
                        f"{type(B).__name__} — must be B_a alone (no B_b arm).")
        N, nx, ng = slab.quad.N, slab.nx, slab.ng
        n_tr = int(slab.angular_trace.layout.total_size)
        psi = FullField(
            interior=AngularFlux.from_mesh(np.random.default_rng(4).standard_normal((N, ng, nx)), slab),
            boundary=AngularBoundaryFlux(
                values=np.random.default_rng(5).standard_normal(n_tr),
                space=slab.angular_trace, mesh=slab),
            radial_characteristic=None,
        )
        out = B.apply(psi)
        if out.radial_characteristic is not None:
            pytest.fail("seedless B_a emitted a non-None ray on a seedless composite.")


class TestB_b_RayBoundary:
    r"""``B_b`` — the ψ½ ray-corner boundary law (RULINGS P1 + P2).

    The specular corner swap, its Euclidean transpose (the mirror), and — the
    load-bearing gate — the ``G_sd = V_cell`` reciprocity that keeps Mode-12
    closed (P2: the corner gauge is symmetric, so Euclidean = Hilbert). Reflective
    sphere for the non-trivial arm; vacuum for the null control.
    """

    def test_reflective_corner_swap_forward(self):
        r"""Forward: ``out.corner(level, −1) = seed.corner(level, +1)`` per level;
        the cells and the +1 corner stay zero (B_b touches only the inflow row).
        B.2b: probed on System B's OWN carrier; the unified view of the output
        (the role-preserving bridge) keeps the pre-B.2b assertions bitwise."""
        sn = _sphere(bc="reflective")
        space = sn.radial_characteristic_space
        seed_vals = np.random.default_rng(6).standard_normal(space.shape[0])
        out = RadialCharacteristicBoundaryOperator(sn).apply(_ray_composite(sn, seed_vals))
        ov = out.to_unified().values
        for level in space.levels:
            np.testing.assert_array_equal(
                space.corner_view(ov, level, -1), space.corner_view(seed_vals, level, +1),
                err_msg=f"level {level}: corner(−1) ≠ seed.corner(+1) (specular swap wrong).")
            np.testing.assert_array_equal(
                space.corner_view(ov, level, +1), 0.0,
                err_msg=f"level {level}: the +1 corner is non-zero (B_b touched a non-inflow row).")
            np.testing.assert_array_equal(
                space.cells_view(ov, level, -1), 0.0,
                err_msg=f"level {level}: the cells leg is non-zero (B_b is corner-only).")

    def test_transpose_is_exact_euclidean_mirror(self):
        r"""``dense(B_b.apply_transpose) ≡ dense(B_b.apply).T`` (0 ULP). A
        same-direction transpose would equal ``dense(apply)`` (not its transpose)
        and — since the swap matrix is non-symmetric — red this gate."""
        sn = _sphere(bc="reflective")
        B_b = RadialCharacteristicBoundaryOperator(sn)
        fwd = _dense_ray(B_b.apply, sn)
        T = _dense_ray(B_b.apply_transpose, sn)
        np.testing.assert_array_equal(
            T, fwd.T, err_msg="B_bᵀ ≠ (B_b)ᵀ — the transpose is not the Euclidean mirror.")

    def test_euclidean_transpose_is_the_vcell_hilbert_adjoint(self):
        r"""Mode-12 CLOSURE (RULING P2): ``⟨B_b x, y⟩_{G_sd} = ⟨x, B_bᵀ y⟩_{G_sd}``
        under ``G_sd = V_cell``. Euclidean IS the Hilbert adjoint because the
        corner gauge is symmetric (``g₊ = g₋ = V(R)``). CONTROL = 0 + two teeth
        (a wrong-direction transpose, an asymmetric gauge) prove it is not vacuous
        — a future asymmetric gauge that reopened Mode-12 would red the gate."""
        sn = _sphere(bc="reflective")
        B_b = RadialCharacteristicBoundaryOperator(sn)
        fwd = _dense_ray(B_b.apply, sn)
        T = _dense_ray(B_b.apply_transpose, sn)
        # The unified-layout G_sd carries the SAME numbers as System B's
        # composite member space (interior V_cell ⊕ boundary V(R)) — the b2
        # member-wise ≡ direct gates pin that equivalence (G-b3.4).
        g = _v_cell_seed(sn)
        rng = np.random.default_rng(8)
        ctrl = _g_recip(fwd, T, g, rng)
        print(f"  B_b G_sd-reciprocity: control={ctrl:.2e}")
        if not (ctrl < 1e-12):
            pytest.fail(f"CONTROL defect {ctrl:.3e} ≠ 0 — the Euclidean transpose is NOT the "
                        f"V_cell Hilbert adjoint; a Euclidean block adjoint on System B has "
                        f"reopened Mode-12 (the corner gauge is not symmetric).")
        # Tooth a: a same-direction (wrong) transpose breaks reciprocity.
        tooth_a = _g_recip(fwd, fwd, g, rng)
        # Tooth b: an asymmetric corner gauge breaks reciprocity even with the correct T.
        g_bad = g.copy()
        for level in sn.radial_characteristic_space.levels:
            sn.radial_characteristic_space.corner_view(g_bad, level, +1)[:] *= 2.0
        tooth_b = _g_recip(fwd, T, g_bad, rng)
        print(f"    teeth: wrong-transpose={tooth_a:.2f}  gauge-asymmetry={tooth_b:.2f}")
        if not (tooth_a > 1e-3 and tooth_b > 1e-3):
            pytest.fail(f"reciprocity gate is VACUOUS: wrong-transpose tooth {tooth_a:.3e} or "
                        f"gauge-asymmetry tooth {tooth_b:.3e} did not red (Mode-12 gate toothless).")

    def test_vacuum_outer_emits_zero_corner(self):
        r"""``kind == "vacuum"`` → ``B_b`` emits an all-zero ray block (no
        re-emission at the outer ray). Positive law (anti-#11)."""
        sn = _sphere(bc="vacuum")
        seed_vals = np.random.default_rng(10).standard_normal(
            sn.radial_characteristic_space.shape[0])
        out = RadialCharacteristicBoundaryOperator(sn).apply(_ray_composite(sn, seed_vals))
        np.testing.assert_array_equal(
            out.to_unified().values, 0.0,
            err_msg="vacuum B_b emitted a non-zero corner (it did the reflective swap).")

    def test_unruled_outer_law_is_loud_deferred(self, monkeypatch):
        r"""``kind ∈ {white, albedo, periodic}`` → ``NotImplementedError`` with the
        specific message (NEGATIVE law, anti-#11: a bare ``raises`` false-greens on
        a downstream crash). Monkeypatch the xmax law kind (no white-sphere mesh
        needed) — auto-reverts, never a git checkout."""
        sn = _sphere(bc="reflective")
        B_b = RadialCharacteristicBoundaryOperator(sn)  # construct BEFORE the patch
        monkeypatch.setattr(sn.bc["xmax"], "kind", "white")
        seed_vals = np.random.default_rng(12).standard_normal(
            sn.radial_characteristic_space.shape[0])
        with pytest.raises(NotImplementedError, match="no ruled corner action yet"):
            B_b.apply(_ray_composite(sn, seed_vals))

    def test_is_adjointable_is_per_leaf(self):
        r"""``B_b.is_adjointable`` is the OUTER ray-face law's, not the whole-trace
        intersection: reflective + vacuum → True; the loud-deferred set → False."""
        if not RadialCharacteristicBoundaryOperator(_sphere(bc="reflective")).is_adjointable:
            pytest.fail("reflective B_b is not adjointable (the involution should be).")
        if not RadialCharacteristicBoundaryOperator(_sphere(bc="vacuum")).is_adjointable:
            pytest.fail("vacuum B_b is not adjointable (the zero map should be).")


class TestSplitInteraction:
    r"""The schedule ``split()`` lives on ``B_a`` alone; ``B_b`` is schedule-atomic.

    RULING P1 corollary: a grading is a refinement of ONE system's boundary block,
    never the composite. Because ``B_a`` sheds the ray arm, its masked halves emit
    ZERO ray — so ``B_lower + B_upper + B_b`` carries the ray corner exactly ONCE
    (the latent double-count the un-weld closes; it would go live on curvilinear
    multi-D, #22). Verified on the seed-carrying sphere's degenerate split (the
    masked ``_apply_faces`` path is what matters — it emits present-zero ray)."""

    def test_split_masked_halves_emit_zero_ray(self):
        from orpheus.sn.loss_representation.sweep_schedule import SweepSchedule

        sn = _sphere(bc="reflective")
        B_a = SNBoundaryOperator(sn)
        parts = B_a.split(SweepSchedule.gauss_seidel(sn))
        psi = _random_composite(sn, np.random.default_rng(13))
        for name, half in (("lower", parts.lower), ("upper", parts.upper)):
            out = half.apply(psi)
            if out.radial_characteristic is None:
                pytest.fail(f"B_{name} emitted a None ray (mixed-presence raises under the sum).")
            np.testing.assert_array_equal(
                out.radial_characteristic.values, 0.0,
                err_msg=f"masked B_{name} emitted a NON-ZERO ray corner — the split doubles "
                        f"the ray (the latent bug; the grading must live on B_a's TRACE only).")


# ── Step-1c helpers: A_BB = RadialCharacteristicOperator (the ψ½ radial BVP) ──


def _graded_sphere(nx: int, ng: int = 2, p: float = 1.5, R: float = 4.0,
                   bc: str = "vacuum"):
    r"""A seed-carrying sphere on a power-graded radial mesh ``r_j = R·(j/nx)^p``.

    The genuinely non-uniform ``dr`` the reversal / index-drift gates need: on a
    uniform mesh ``dr[::-1] == dr`` and ``dr[k−1] == dr[k]``, so those gates are
    vv Mode-5 vacuous — the grading breaks the blind spot (§0.6)."""
    edges = R * (np.arange(nx + 1, dtype=float) / nx) ** p
    mesh = Mesh1D(edges=edges, mat_ids=np.zeros(nx, dtype=int),
                  coord=CoordSystem.SPHERICAL, bc_right=BC(bc))
    return SNMesh(mesh, Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, ng)})


def _ray_sigma(sn, slope: float = 0.3) -> CrossSectionField:
    r"""Heterogeneous per-group per-cell ``σ_t`` as a typed, mesh-bound
    ``CrossSectionField`` on ``sn`` (``.values`` are ``(ng, nx)``, varying in
    BOTH group AND cell so an index / group-axis bug in the march is not nulled —
    anti-#2 asymmetry). The operator's ``C_ray`` collision coefficient; the typed
    field carries ``.mesh`` for the operator's mesh-identity guard."""
    nx, ng = sn.nx, sn.ng
    raw = np.stack([1.0 + slope * g + 0.15 * np.arange(nx) for g in range(ng)], 0)
    return CrossSectionField.from_mesh(raw, sn)


def _ray_source(sn, rng) -> RadialCharacteristicSourceSink:
    """A random q½ source block on ``sn``'s ray carrier (all slots non-zero)."""
    space = sn.radial_characteristic_space
    return RadialCharacteristicSourceSink(
        values=rng.standard_normal(space.shape[0]), space=space, mesh=sn)


def _ray_cotangent(sn, rng) -> RadialCharacteristicFlux:
    """A random flux-space cotangent (the solve's codomain) on ``sn``'s carrier."""
    space = sn.radial_characteristic_space
    return RadialCharacteristicFlux(
        values=rng.standard_normal(space.shape[0]), space=space, mesh=sn)


def _member(unified) -> RadialCharacteristicComposite:
    """Bridge a unified ψ½ leaf into System B's member composite — the B.2c
    re-typed block I/O (role-preserving + bitwise, the B.1d licence). The
    oracles stay UNIFIED (engine layout); only the operator calls bridge."""
    return RadialCharacteristicComposite.from_unified(unified)


def _two_leg_reference(op, source) -> RadialCharacteristicFlux:
    r"""Replicate ``A_BB.solve``'s two-leg march with the REAL engine — the WRAP
    oracle for the bit-identity gate. Calls ``carlson_inward_sweep_from_source``
    (the test-module imported name, UNPATCHED when the operator's module attr is
    spied) so ``.solve`` (patched with a delegating spy) and this reference
    compute with the SAME engine → any divergence is a WRAP bug, not FP."""
    space = source.space
    sigma = op.total_cross_section.values
    dr = np.asarray(op.sn_mesh.axis_widths[0])
    out = RadialCharacteristicFlux.zeros_on(op.sn_mesh)
    buf, sv = out.values, source.values
    for lv in space.levels:
        q_minus = space.cells_view(sv, lv, -1)
        q_plus = space.cells_view(sv, lv, +1)
        corner_in = space.corner_view(sv, lv, -1)
        cells_minus, pole_face = carlson_inward_sweep_from_source(
            q_minus, sigma, dr, corner_in)
        cells_plus_rev, corner_out = carlson_inward_sweep_from_source(
            q_plus[:, ::-1], sigma[:, ::-1], dr[::-1], pole_face)
        space.cells_view(buf, lv, -1)[...] = cells_minus
        space.corner_view(buf, lv, -1)[...] = corner_in
        space.cells_view(buf, lv, +1)[...] = cells_plus_rev[:, ::-1]
        space.corner_view(buf, lv, +1)[...] = corner_out
    return out


def _install_engine_spy(monkeypatch) -> list[dict]:
    r"""Mode-11 sentinel: wrap ``carlson_inward_sweep_from_source`` in the
    OPERATOR's module namespace, recording ``(args, result)`` per call and
    delegating to the real engine. Returns the calls list (2 per level: inward
    then outward). Proves ``.solve`` EXECUTES the production engine (a divergent
    inlined copy would leave the list empty)."""
    calls: list[dict] = []
    real = carlson_inward_sweep_from_source  # the unpatched module-top import

    def spy(Q_bar, sigma_t, dr, bc_outer_value):
        result = real(Q_bar, sigma_t, dr, bc_outer_value)
        calls.append({
            "Q": np.asarray(Q_bar).copy(),
            "sigma": np.asarray(sigma_t).copy(),
            "dr": np.asarray(dr).copy(),
            "bc": np.asarray(bc_outer_value).copy(),
            "cells": result[0].copy(),
            "exit_face": result[1].copy(),
        })
        return result

    monkeypatch.setattr(_rc_mod, "carlson_inward_sweep_from_source", spy)
    return calls


def _euclid_adjoint_defect(op, u, v) -> float:
    r"""The relative Euclidean reciprocity defect
    ``|⟨solve(u), v⟩ − ⟨u, solve_transpose(v)⟩| / (|·| + |·|)``.

    Plain dot products (NOT the ``G_sd`` metric): ``solve_transpose`` is the
    ISOLATED EUCLIDEAN adjoint of the resolvent (the pure ray-block transpose —
    operator docstring), so its consistency partner is the Euclidean inner
    product, not the ``V_cell`` Hilbert adjoint (which is realized once at the
    composite, L19). ``u``/``v`` are unified; the calls bridge (B.2c I/O)."""
    lhs = float(op.solve(_member(u)).to_unified().values @ v.values)
    rhs = float(u.values @ op.solve_transpose(_member(v)).to_unified().values)
    return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)


class TestA_BB_RadialBVP:
    r"""``A_BB`` = :class:`RadialCharacteristicOperator` — the ψ½ radial two-point
    BVP resolvent (campaign step 1c).

    Foundation invariants of the full operator — the forward ``apply`` /
    ``apply_transpose``, the resolvent ``solve`` / ``solve_transpose``, and the
    operator-returning ``inverse()`` (step 1c posed the resolvent; step 4b
    completes the forward via the shared kernel).
    The convergence-ORDER claim lives in the sibling L1 module
    ``test_ray_operator.py`` (don't conflate foundation + verifies, L9). Every
    value row is ≥2G; the sphere-GL S4 carrier is the ONLY seed-carrying member
    (cylinder/slab are the non-carrying CONTROL — the constructor rejects them).

    Runtime: gates raise via :func:`pytest.fail` / ``np.testing.assert_*`` (fire
    under ``python -O``), never a bare ``assert`` (vv Mode 8).
    """

    # ── solve / solve_transpose adjoint consistency (Euclidean) ────────────

    def test_adjoint_consistency_euclidean(self):
        r"""``⟨solve(u), v⟩ = ⟨u, solve_transpose(v)⟩`` (Euclidean) to < 1e-11 on
        heterogeneous σ, ≥2 random draws — the PRIMARY solve/solve_transpose
        consistency gate (the resolvent adjoint is the reverse-mode transpose of
        the two-leg march). The source ``μ = +1`` corner cotangent is EXACTLY
        zero (the q½ fold writes only cells + the ``μ = -1`` corner, so the
        outflow-corner source slot is unused — R13)."""
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        space = sn.radial_characteristic_space
        for seed in (1, 2, 3):                       # ≥2 draws
            rng = np.random.default_rng(seed)
            u, v = _ray_source(sn, rng), _ray_cotangent(sn, rng)
            defect = _euclid_adjoint_defect(op, u, v)
            if not defect < 1e-11:
                pytest.fail(
                    f"seed {seed}: Euclidean adjoint defect {defect:.3e} ≥ 1e-11 — "
                    f"solve_transpose is NOT the transpose of solve (the reverse "
                    f"march or its leg chaining is mis-wired).")
            src_bar = op.solve_transpose(_member(v)).to_unified()
            for lv in space.levels:
                np.testing.assert_array_equal(
                    space.corner_view(src_bar.values, lv, +1), 0.0,
                    err_msg=f"seed {seed} level {lv}: the source μ=+1 corner "
                            f"cotangent is non-zero (it must stay 0 — the q½ fold "
                            f"never writes that slot).")

    def test_adjoint_sign_flip_tooth(self, monkeypatch):
        r"""TOOTH for the adjoint gate: a sign flip in the reverse-mode
        recurrence (``carlson_inward_sweep_transpose``'s incoming face-cotangent
        ``-f_bar → +f_bar``) breaks reciprocity — the defect jumps to O(1).
        Proves the < 1e-11 gate above is not vacuously green."""

        def transpose_sign_flip(cells_bar, final_face_bar, sigma_t, dr):
            ng, nx = cells_bar.shape
            Q_bar = np.zeros((ng, nx), dtype=cells_bar.dtype)
            f_bar = final_face_bar.copy()
            for k in range(nx):
                denom = dr[k] * sigma_t[:, k] + 2.0
                c_bar = cells_bar[:, k] + 2.0 * f_bar
                f_in_bar = +f_bar                    # SIGN FLIP: production is -f_bar
                Q_bar[:, k] = (dr[k] / denom) * c_bar
                f_in_bar = f_in_bar + (2.0 / denom) * c_bar
                f_bar = f_in_bar
            return Q_bar, f_bar

        monkeypatch.setattr(
            _rc_mod, "carlson_inward_sweep_transpose", transpose_sign_flip)
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        rng = np.random.default_rng(1)
        defect = _euclid_adjoint_defect(op, _ray_source(sn, rng),
                                        _ray_cotangent(sn, rng))
        if not defect > 1e-3:
            pytest.fail(
                f"the transpose sign-flip left the adjoint defect at {defect:.3e} "
                f"— the adjoint-consistency gate has no teeth.")

    # ── WRAP bit-identity via a Mode-11 call-counter sentinel ──────────────

    def test_wrap_executes_engine_bit_identical(self, monkeypatch):
        r"""``solve`` WRAPS the production engine: the Mode-11 sentinel counts
        exactly ``2·n_levels`` calls to ``carlson_inward_sweep_from_source`` (2
        legs/level) AND the result is bit-identical (``array_equal``) to an
        independent two-leg reference on the SAME engine — a divergent inlined
        copy would leave the counter at 0 (Cardinal Rule 2 single source)."""
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        source = _ray_source(sn, np.random.default_rng(4))
        reference = _two_leg_reference(op, source)   # real engine, before the spy
        calls = _install_engine_spy(monkeypatch)
        flux = op.solve(_member(source)).to_unified()
        n_levels = len(sn.radial_characteristic_space.levels)
        if len(calls) != 2 * n_levels:
            pytest.fail(
                f"solve called the engine {len(calls)}× (expected 2·n_levels = "
                f"{2 * n_levels}) — it is NOT the two-leg WRAP (a divergent copy?).")
        np.testing.assert_array_equal(
            flux.values, reference.values,
            err_msg="solve is not bit-identical to the two-leg engine reference — "
                    "the WRAP diverged from the production march.")

    def test_pole_continuation_threads_exit_to_entry(self, monkeypatch):
        r"""Pole continuation ``ψ½⁺(0) = ψ½⁻(0)``: per level the OUTWARD leg's
        entry face (call #2's ``bc_outer_value``) EQUALS the INWARD leg's exit
        face (call #1's ``phi_face_final``) — the inward exit IS the outward
        entry (internal to the march, R13). The exit face is asserted non-trivial
        so the gate is not vacuously satisfied by zeros."""
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        calls = _install_engine_spy(monkeypatch)
        op.solve(_member(_ray_source(sn, np.random.default_rng(5))))
        for i in range(0, len(calls), 2):
            inward, outward = calls[i], calls[i + 1]
            np.testing.assert_array_equal(
                outward["bc"], inward["exit_face"],
                err_msg="the outward leg's entry face ≠ the inward leg's exit face "
                        "— pole continuation ψ½⁺(0)=ψ½⁻(0) is broken.")
            if not np.max(np.abs(inward["exit_face"])) > 0.0:
                pytest.fail("the inward exit (pole) face is identically zero — the "
                            "pole-continuation gate is vacuous for this source.")

    def test_outward_leg_marches_reversed_data(self, monkeypatch):
        r"""The 2.5a discipline — orientation is carried by the DATA, never a
        flag: the OUTWARD (+1) leg rides the same engine on the ``[:, ::-1]`` /
        ``[::-1]`` reversed level data, the INWARD (-1) leg on forward data. Run
        on a GRADED mesh so ``dr[::-1] ≠ dr`` — the reversal is a genuine
        constraint (on a uniform mesh it is vv Mode-5 vacuous; the non-vacuity
        check enforces that). If the operator dropped a reversal, call #2 would
        carry forward data and these equalities RED."""
        sn = _graded_sphere(nx=8)
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        space = sn.radial_characteristic_space
        dr = np.asarray(sn.axis_widths[0])
        sigma = op.total_cross_section.values
        source = _ray_source(sn, np.random.default_rng(6))
        sv = source.values
        # Non-vacuity (Mode 5): on this graded mesh reversed ≠ forward.
        if np.array_equal(dr[::-1], dr):
            pytest.fail("dr is uniform — the reversal gate is Mode-5 vacuous; the "
                        "graded mesh must give dr[::-1] ≠ dr.")
        calls = _install_engine_spy(monkeypatch)
        op.solve(_member(source))
        for idx, lv in enumerate(space.levels):
            inward, outward = calls[2 * idx], calls[2 * idx + 1]
            np.testing.assert_array_equal(
                inward["dr"], dr,
                err_msg=f"level {lv}: the inward leg did not march FORWARD widths.")
            np.testing.assert_array_equal(
                inward["Q"], space.cells_view(sv, lv, -1),
                err_msg=f"level {lv}: the inward leg read the wrong source cells.")
            np.testing.assert_array_equal(
                outward["dr"], dr[::-1],
                err_msg=f"level {lv}: the outward leg did not march REVERSED widths "
                        f"(the 2.5a data-carried orientation broke).")
            np.testing.assert_array_equal(
                outward["Q"], space.cells_view(sv, lv, +1)[:, ::-1],
                err_msg=f"level {lv}: the outward leg did not read REVERSED source cells.")
            np.testing.assert_array_equal(
                outward["sigma"], sigma[:, ::-1],
                err_msg=f"level {lv}: the outward leg did not read REVERSED σ_t.")

    # ── r = R Dirichlet propagation ────────────────────────────────────────

    def test_r_R_dirichlet_propagates_into_interior(self):
        r"""A nonzero ``r = R`` inflow corner (μ=−1) vs zero — same cells source —
        changes the INTERIOR cells (by the ``e^{−σ(R−r)}`` envelope), not merely
        the boundary. Two solves differing ONLY in ``corner_in`` must differ in
        the interior; equal interiors would mean the Dirichlet datum is ignored."""
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        space = sn.radial_characteristic_space
        s0 = RadialCharacteristicSourceSink(
            values=np.zeros(space.shape[0]), space=space, mesh=sn)
        s1 = RadialCharacteristicSourceSink(
            values=np.zeros(space.shape[0]), space=space, mesh=sn)
        for s in (s0, s1):
            for lv in space.levels:
                space.cells_view(s.values, lv, -1)[...] = 0.5     # identical cells
        for lv in space.levels:
            space.corner_view(s1.values, lv, -1)[...] = 3.0       # nonzero inflow only
        a = op.solve(_member(s0)).interior.cells(0, -1)
        b = op.solve(_member(s1)).interior.cells(0, -1)
        interior_diff = float(np.max(np.abs(a[:, :-1] - b[:, :-1])))  # exclude outer cell
        if not interior_diff > 1e-6:
            pytest.fail(
                f"the interior cells are unchanged ({interior_diff:.3e}) when the "
                f"r=R inflow corner goes 0 → 3 — the Dirichlet datum does NOT "
                f"propagate inward (the corner is being ignored).")

    def test_dirichlet_bc_ignore_tooth(self, monkeypatch):
        r"""TOOTH for the Dirichlet-propagation gate: an engine that ignores its
        ``bc_outer_value`` (always enters at 0) makes the two solves' interiors
        IDENTICAL — the interior difference collapses to 0. Proves the gate
        above catches a dropped inflow datum."""

        def ignore_bc(Q_bar, sigma_t, dr, bc_outer_value):
            return carlson_inward_sweep_from_source(
                Q_bar, sigma_t, dr, np.zeros_like(bc_outer_value))

        monkeypatch.setattr(
            _rc_mod, "carlson_inward_sweep_from_source", ignore_bc)
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        space = sn.radial_characteristic_space
        s0 = RadialCharacteristicSourceSink(
            values=np.zeros(space.shape[0]), space=space, mesh=sn)
        s1 = RadialCharacteristicSourceSink(
            values=np.zeros(space.shape[0]), space=space, mesh=sn)
        for s in (s0, s1):
            space.cells_view(s.values, 0, -1)[...] = 0.5
        space.corner_view(s1.values, 0, -1)[...] = 3.0
        interior_diff = float(np.max(np.abs(
            op.solve(_member(s0)).interior.cells(0, -1)[:, :-1]
            - op.solve(_member(s1)).interior.cells(0, -1)[:, :-1])))
        if not interior_diff < 1e-14:
            pytest.fail(
                f"the bc-ignoring engine still produced an interior difference "
                f"({interior_diff:.3e}) — the Dirichlet-propagation gate's tooth "
                f"does not bite.")

    # ── Fixed-source Q/Σ equilibrium (conservation + spatial distribution) ─

    def test_fixed_source_equilibrium_Q_over_sigma(self):
        r"""The single most powerful curvilinear diagnostic: uniform source at
        equilibrium ``q̄ = σ·C`` with the consistent inflow ``φ_R = C`` → every
        cell of BOTH legs sits at ``C = q̄/σ`` (the flat identity
        ``(Δr·σ·C + 2C)/(Δr·σ + 2) = C``, self-similar through the pole). ≥2G with
        DISTINCT per-group ``C`` and heterogeneous per-cell σ, so a missing ``Δr``
        / factor / group-axis bug would break the equilibrium."""
        sn = _sphere(nx=6)
        sigma = _ray_sigma(sn)                       # het in g AND cell
        sig = sigma.values                           # (ng, nx) for the σ·C source
        op = RadialCharacteristicOperator(sn, sigma)
        space = sn.radial_characteristic_space
        C = np.array([0.5, 1.3])                     # distinct per-group equilibrium
        src = RadialCharacteristicSourceSink(
            values=np.zeros(space.shape[0]), space=space, mesh=sn)
        for g in range(sn.ng):
            for sign in (-1, +1):                    # BOTH legs' cells source = σ·C
                space.cells_view(src.values, 0, sign)[g, :] = sig[g] * C[g]
            space.corner_view(src.values, 0, -1)[g] = C[g]   # consistent inflow
        flux = op.solve(_member(src))
        expected = np.broadcast_to(C[:, None], (sn.ng, sn.nx))
        for sign in (-1, +1):
            np.testing.assert_allclose(
                flux.interior.cells(0, sign), expected, atol=1e-13,
                err_msg=f"leg {sign}: the equilibrium flux ≠ q̄/σ = C — the "
                        f"fixed-source Q/Σ balance (conservation + spatial "
                        f"distribution) is broken.")

    # ── Constructor negative gates (the non-carrying CONTROL) ──────────────

    def test_constructor_rejects_non_carrying_foreign_mesh_and_nonpositive(self):
        r"""The constructor's guards, NET-NEW teeth (L4). Three illegal states,
        each with ``match=`` the SPECIFIC message (a downstream crash would
        false-green a bare ``raises``):

        * **non-carrying CONTROL** — a cylinder (level-symmetric) and a slab
          (Cartesian) have ``radial_characteristic_space is None`` → the seedless
          guard fires (before σ_t is even read);
        * **foreign-mesh σ_t** — a ``CrossSectionField`` on a DIFFERENT sphere
          (same ``(ng, nx)``, different Δr) is refused by the mesh-identity
          invariant. THIS is the Pattern-4 illegal state the typed, mesh-bound
          coefficient closes — a bare ``(ng, nx)`` ndarray could not catch it (it
          carries no mesh), so the operator would silently march this mesh's Δr
          against a foreign σ_t;
        * **σ_t ≤ 0** — the DD-denominator ``Δr·σ + 2`` guard.

        Positive control: a σ_t on THIS mesh constructs cleanly."""
        # Non-carrying CONTROL — a cylinder (needs a level-structured quadrature).
        cyl = SNMesh(
            Mesh1D(edges=np.linspace(0.05, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CYLINDRICAL, bc_right=BC("vacuum")),
            Quadrature.level_symmetric(4), {0: _mixture(1.0, 0.4, 2)})
        slab = SNMesh(
            Mesh1D(edges=np.linspace(0.0, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CARTESIAN, bc_right=BC("reflective"),
                   bc_left=BC("reflective")),
            Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, 2)})
        if cyl.radial_characteristic_space is not None:
            pytest.fail("the cylinder carries a ray space — CONTROL invalid.")
        if slab.radial_characteristic_space is not None:
            pytest.fail("the slab carries a ray space — CONTROL invalid.")
        # The seedless guard fires before σ_t is read (a valid field on the mesh).
        with pytest.raises(ValueError, match="carries no starting-direction ray"):
            RadialCharacteristicOperator(
                cyl, CrossSectionField.from_mesh(np.ones((2, 5)), cyl))
        with pytest.raises(ValueError, match="carries no starting-direction ray"):
            RadialCharacteristicOperator(
                slab, CrossSectionField.from_mesh(np.ones((2, 5)), slab))
        # Positive control — a σ_t on THIS mesh constructs cleanly.
        sn = _sphere(nx=6)
        RadialCharacteristicOperator(
            sn, CrossSectionField.from_mesh(np.ones((sn.ng, sn.nx)), sn))
        # THE Pattern-4 closure: a σ_t bound to a DIFFERENT sphere (graded — a
        # genuinely different Δr) is refused. The typed coefficient makes the
        # foreign-mesh march unconstructable; a bare ndarray could not.
        foreign_mesh = _graded_sphere(nx=6)
        foreign_sigma = CrossSectionField.from_mesh(
            np.ones((sn.ng, sn.nx)), foreign_mesh)
        with pytest.raises(ValueError, match="mesh-identity invariant"):
            RadialCharacteristicOperator(sn, foreign_sigma)
        # σ_t ≤ 0 → the DD-denominator guard.
        bad = np.ones((sn.ng, sn.nx))
        bad[1, 2] = 0.0
        with pytest.raises(ValueError, match="strictly positive"):
            RadialCharacteristicOperator(sn, CrossSectionField.from_mesh(bad, sn))


def _install_forward_spy(monkeypatch) -> list[int]:
    r"""Mode-11 anti-twin sentinel: wrap the SHARED forward kernel
    ``radial_characteristic_forward_residual`` in BOTH namespaces that import
    it — the operator (``_rc_mod``) and the fused ``(L+C)`` walk (``_lr_mod``).
    Every entry appends to the shared list, so a caller that inlined a
    divergent copy would leave its delta at 0 (L16: only-new-reds ⟹ twin)."""
    calls: list[int] = []
    real = radial_characteristic_forward_residual

    def spy(values, space, sigma, dr):
        calls.append(1)
        return real(values, space, sigma, dr)

    monkeypatch.setattr(_rc_mod, "radial_characteristic_forward_residual", spy)
    monkeypatch.setattr(_lr_mod, "radial_characteristic_forward_residual", spy)
    return calls


class TestA_BB_Forward:
    r"""``A_BB`` step-4b — the forward ``apply`` / ``apply_transpose`` /
    ``inverse()`` that complete the operator, single-sourced with the ``(L+C)``
    walk (the user's "extract the shared kernel now" ruling — no forward twin).

    Sphere-GL S4 carrier, ≥2G, graded σ. ``solve∘apply`` is principled-equiv at
    ~FP ULP for the cells (the forward's ``2/Δr`` and the march's ``Δr·σ+2``
    reassociate, L7) and BIT-EXACT ``0.0`` on the μ=+1 outflow corner.
    """

    def test_apply_is_the_exact_march_inverse(self):
        # solve∘apply=id on the CONSISTENT subspace ψ0 = solve(q0) (the +1
        # outflow corner is a FREE datum apply MEASURES / solve OVERWRITES — an
        # arbitrary ψ falsely reds, refutation R2). Cells at rtol
        # (principled-equiv, R1); the apply∘solve +1 corner closes bit-exact 0.
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        space = sn.radial_characteristic_space
        for seed in range(2):
            rng = np.random.default_rng(seed)
            q0 = _member(_ray_source(sn, rng))
            psi0 = op.solve(q0)
            psi1 = op.solve(op.apply(psi0))
            np.testing.assert_allclose(
                psi1.to_unified().values, psi0.to_unified().values,
                rtol=1e-11, atol=1e-13)
            qr = op.apply(op.solve(q0))
            for p in space.levels:
                corner = qr.boundary.corner(p, +1)
                np.testing.assert_array_equal(corner, np.zeros_like(corner))

    def test_apply_routes_through_the_shared_forward_kernel(self, monkeypatch):
        # Mode-11 anti-twin (operator side): A_BB.apply MUST enter the shared
        # radial_characteristic_forward_residual, not an inlined copy.
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        calls = _install_forward_spy(monkeypatch)
        op.apply(_member(_ray_cotangent(sn, np.random.default_rng(1))))
        if len(calls) == 0:
            pytest.fail(
                "A_BB.apply did NOT route through the shared forward kernel "
                "(a divergent inlined copy — twin).")

    def test_walk_forward_routes_through_the_shared_kernel(self, monkeypatch):
        # Mode-11 anti-twin (walk side): (L+C).apply's ψ½ rows MUST enter the
        # SAME shared kernel — the single-source proof (both callers, one body).
        sn = _sphere()
        LC = _loss(sn, slope=0.4)
        calls = _install_forward_spy(monkeypatch)
        LC.apply(_random_composite(sn, np.random.default_rng(2)))
        if len(calls) == 0:
            pytest.fail(
                "(L+C).apply did NOT route through the shared forward kernel "
                "(a leftover inline copy survives — twin).")

    def test_apply_transpose_is_the_euclidean_adjoint(self):
        # ⟨apply(u), v⟩ = ⟨u, apply_transpose(v)⟩ — plain flat dot (the metric
        # Hilbert adjoint .H is realized at the composite, L19/R4). Graded het σ.
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        for seed in range(2):
            rng = np.random.default_rng(seed + 10)
            u = _ray_cotangent(sn, rng)
            v = _ray_cotangent(sn, rng)
            lhs = float(op.apply(_member(u)).to_unified().values @ v.values)
            rhs = float(
                u.values @ op.apply_transpose(_member(v)).to_unified().values)
            defect = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
            if defect > 1e-11:
                pytest.fail(
                    f"forward Euclidean adjoint defect {defect:.2e} > 1e-11")

    def test_inverse_is_the_march_involution(self):
        # inverse() delegates: apply → inner.solve, solve → inner.apply;
        # inverse().inverse() is self (mixin identity); predicates report True.
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        rng = np.random.default_rng(3)
        q0 = _member(_ray_source(sn, rng))
        psi0 = op.solve(q0)
        inv = op.inverse()
        np.testing.assert_array_equal(
            inv.apply(q0).to_unified().values, op.solve(q0).to_unified().values)
        np.testing.assert_array_equal(
            inv.solve(psi0).to_unified().values,
            op.apply(psi0).to_unified().values)
        if inv.inverse() is not op:
            pytest.fail("inverse().inverse() is not the original operator.")
        if not (op.is_invertible and op.is_adjointable):
            pytest.fail(
                "A_BB must report is_invertible and is_adjointable True at 4b.")

    def test_b2c_member_composite_block_boundary(self):
        r"""B.2c re-type (G-c1.1): the four action surfaces speak System B's
        member composite. Containers (source composite out of apply /
        apply_transpose / solve_transpose; flux composite out of solve);
        declared domain/codomain asserted by object IDENTITY, not ``==`` — the
        unified space collides with the composite space on ``(name, shape)``
        (memo F2), so ``==`` is Mode-12-blind to a block left typed on the
        unified carrier; and the block-boundary refusals (foreign carrier /
        foreign mesh / the solve source-role parse)."""
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        rng = np.random.default_rng(7)
        src = _member(_ray_source(sn, rng))
        cot = _member(_ray_cotangent(sn, rng))
        if op.domain is not sn.radial_characteristic_composite_space:
            pytest.fail("A_BB.domain is not THE composite member-space object "
                        "(F2: == cannot see a unified-typed block).")
        if op.codomain is not sn.radial_characteristic_composite_space:
            pytest.fail("A_BB.codomain is not THE composite member-space object.")
        for name, out in (
            ("apply", op.apply(cot)),
            ("apply_transpose", op.apply_transpose(cot)),
            ("solve_transpose", op.solve_transpose(cot)),
        ):
            if type(out) is not RadialCharacteristicComposite:
                pytest.fail(f"{name} did not emit the member composite; got "
                            f"{type(out).__name__}.")
            if type(out.interior) is not RadialCharacteristicInteriorSourceSink:
                pytest.fail(f"{name} did not emit SOURCE members; got "
                            f"{type(out.interior).__name__}.")
        out_solve = op.solve(src)
        if type(out_solve) is not RadialCharacteristicComposite:
            pytest.fail("solve did not emit the member composite.")
        if type(out_solve.interior) is not RadialCharacteristicInteriorFlux:
            pytest.fail(f"solve did not emit FLUX members; got "
                        f"{type(out_solve.interior).__name__}.")
        # Block-boundary refusals (parse-don't-validate, match= the parse's
        # own message — a downstream crash must not false-green these).
        with pytest.raises(TypeError, match="System B's member carrier"):
            op.apply(_ray_cotangent(sn, rng))          # a unified leaf
        with pytest.raises(TypeError, match="SOURCE members"):
            op.solve(cot)                              # flux into the resolvent
        with pytest.raises(ValueError, match="mesh-identity invariant"):
            op.apply(_member(_ray_cotangent(_sphere(), rng)))

    def test_b2c_member_composite_bridge_teeth(self, monkeypatch):
        r"""TEETH for the G-c1.1 rows: (a) an apply that SKIPS the
        ``from_unified`` re-split (returns the unified engine leaf) reds the
        container assertion; (b) a value-corrupting bridge (the re-split
        doubling the interior member) moves the apply value — the value rows
        are not blind to a broken bridge. In-process monkeypatch, both legs."""
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        cot = _member(_ray_cotangent(sn, np.random.default_rng(8)))
        reference = op.apply(cot).to_unified().values      # unpatched control
        # (a) bridge-drop → the container gate reds.
        with monkeypatch.context() as m:
            m.setattr(RadialCharacteristicComposite, "from_unified",
                      classmethod(lambda cls, unified: unified))
            out = op.apply(cot)
            if isinstance(out, RadialCharacteristicComposite):
                pytest.fail("the bridge-drop tooth is mis-wired — apply still "
                            "emitted a composite under a dropped re-split.")
        # (b) value corruption → the value rows red.
        real_from_unified = RadialCharacteristicComposite.from_unified.__func__

        def corrupt(cls, unified):
            out = real_from_unified(cls, unified)
            return replace(out, interior=out.interior * 2.0)

        with monkeypatch.context() as m:
            m.setattr(RadialCharacteristicComposite, "from_unified",
                      classmethod(corrupt))
            corrupted = op.apply(cot).to_unified().values
        if not np.max(np.abs(corrupted - reference)) > 0.0:
            pytest.fail("the value-corruption tooth left apply unchanged — the "
                        "value rows have no teeth against a broken bridge.")


# ── Step-2 helpers: A_BA the ψ½ Schur fold (bulk → ray q½ source) ──────────
#
# A_BA folds a bulk isotropic cell-emission q₀ (the ℓ=0 K_iso·φ₀ for S, the
# χ·νΣf·φ for F) onto the ray q½ source at the closed rays μ=±1, per carrying
# radial level. It lives ENTIRELY in the lagged S + F/k gain, OUTSIDE the
# resolvent (it cannot touch the #284 forward-substitution certificate). The
# Step-2 un-weld routes the three hand-rolled S/F fold arms (scattering.py S-fwd
# + S-adj, fission.py F-fwd) through ONE single source — the operator
# RadialCharacteristicReconstruction, which wraps the fold helper
# fold_moments_to_radial_characteristic; these gates pin that. (from_angular_
# source keeps its own per-level analysis — a distinct operation, not a twin.)
#
# The fold math (single source ``fold_moments_to_radial_characteristic``):
#   Q̄(μ=±1) = Σ_ℓ (2ℓ+1)/2 · Q_ℓ · P_ℓ(±1) = Σ_ℓ (2ℓ+1)/2 · Q_ℓ · (±1)^ℓ.
# At ℓ=0 it collapses to ½·Q₀ (both signs — P₀≡1); the PRODUCTION S/F arms feed
# ℓ=0 ONLY, so an S/F-only gate is BLIND to P_ℓ(±1) for ℓ≥1 (refutation #3) —
# the contract gate manufactures ℓ≥1 to activate the ``sign^ℓ`` line.


def _fissile_mixture(sig_t: float, sig_s: float, ng: int):
    """A group-graded FISSILE mixture (asymmetric in g).

    F's seed emission is ``χ·νΣf·φ`` — identically zero on a non-fissile
    mixture (``_mixture`` sets ``sig_f = nu = chi = 0``), which makes the F
    routing / bit-identity rows VACUOUS (0 == 0). vv refutation #4: split the
    fixed-source ``Q/Σ`` scattering config (non-fissile ``_mixture``) from the
    fissile config (here) so the F emission is genuinely nonzero. χ births all
    fission neutrons in group 0 (``chi=[1,0,…]``) — an asymmetric emission.
    """
    st = np.array([sig_t * (1.0 + 0.4 * g) for g in range(ng)])
    ss = np.diag([sig_s * (1.0 + 0.4 * g) for g in range(ng)])
    sf = np.array([0.05 + 0.1 * g for g in range(ng)])
    nu = np.full(ng, 2.4)
    chi = np.zeros(ng)
    chi[0] = 1.0
    return make_mixture(sig_t=st, sig_c=st - ss.sum(axis=0) - sf, sig_f=sf,
                        nu=nu, chi=chi, sig_s=ss)


def _fissile_sphere(nx: int = 5, ng: int = 2, sigma: float = 1.0, c: float = 0.4):
    """A seed-carrying FISSILE sphere (GL S4) — the F-arm carrier (non-vacuous emission)."""
    mesh = Mesh1D(edges=np.linspace(0.0, 4.0, nx + 1), mat_ids=np.zeros(nx, dtype=int),
                  coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"))
    return SNMesh(mesh, Quadrature.gauss_legendre(4),
                  {0: _fissile_mixture(sigma, c * sigma, ng)})


def _s_emission(S, psi: FullField) -> NDArray:
    """The ℓ=0 iso cell-emission ``q₀ = K_iso·φ₀`` that the S seed arm folds — the
    A_BA *input* (``(ng, nx)``). Computed via the operator's own isotropic kernel
    (the emission is the bulk-scattering job, verified elsewhere; A_BA's job is the
    FOLD of this emission, so it is the correct oracle input)."""
    phi0 = psi.interior.integrate_angular().values          # (ng, nx)
    return np.asarray(S.isotropic_kernel.apply(phi0))


def _f_emission(F, psi: FullField) -> NDArray:
    """The ℓ=0 iso fission cell-emission ``χ·νΣf·φ`` that the F seed folds — the
    A_BA_fission *input* (``(ng, nx)``). The eigenvalue outer loop computes exactly
    this (as ``fission_source`` from the SCALAR flux, ``F.kernel ∘ integrate``), so
    the migrated F seed :func:`~orpheus.sn.solver._radial_characteristic_fission_seed`
    only needs the FOLD of this emission (``A_BA_fission = Fold ∘ F.kernel``,
    factored)."""
    return np.asarray(F.apply(psi.interior.integrate_angular()).values)


def _ba_oldloop_reference(emission: NDArray, sn) -> NDArray:
    r"""The EXACT pre-un-weld hand-rolled fold loop (the S/F seed arms before
    Step 2 routed them through RadialCharacteristicReconstruction): the
    bit-identity ORACLE. For each carrying (level, sign)
    the cells are ``fold(emission[None], sign)`` and the corners stay zero. The
    Step-2 un-weld MUST reproduce this byte-for-byte (``np.array_equal``),
    inheriting verification from the ℓ-fold contract gate (vv §Bit-identity: free
    verification by inheritance from a verified reference)."""
    space = sn.radial_characteristic_space
    vals = np.zeros(space.shape[0])
    for lv in space.levels:
        for sign in (-1, +1):
            space.cells_view(vals, lv, sign)[:] = (
                _rcs_mod.fold_moments_to_radial_characteristic(emission[None], sign))
    return vals


def _fold_transpose_reference(y: NDArray, sign: int, n_moments: int,
                              *, coeff0: float | None = None) -> NDArray:
    r"""The Euclidean transpose of the ℓ-fold: ``moments_bar[ℓ] = coeff[ℓ]·y``
    with ``coeff[ℓ] = ((2ℓ+1)/2)·sign^ℓ`` — the contract the production
    fold-transpose (the S-adjoint's single-sourced ``0.5``) must satisfy. The
    ``coeff0`` override is the tooth: a ``0.6`` at ℓ=0 (≠ the fold's ``0.5``)
    breaks the adjoint identity."""
    ell = np.arange(n_moments)
    coeff = ((2.0 * ell + 1.0) / 2.0) * np.float64(sign) ** ell
    if coeff0 is not None:
        coeff = coeff.copy()
        coeff[0] = coeff0
    return coeff[:, None, None] * y[None]


def _install_fold_spy(monkeypatch) -> dict:
    r"""Mode-11 sentinel: WRAP the shared Legendre fold as the RECONSTRUCTION
    operator sees it. Post-un-weld the S/F seed arms route through
    ``RadialCharacteristicReconstruction.apply``, which module-level-imports
    ``fold_moments_to_radial_characteristic`` — so the object it fetches at call
    time is the reconstruction module's global (``_rcr_mod``), NOT the numerics
    source module. Wrapping it here counts every fold the reconstruction runs; an
    S/F arm that bypassed the reconstruction (re-inlined the numerics fold in its
    own namespace) would leave this counter at 0 (Cardinal Rule 2 single source;
    Mode 11)."""
    calls: dict = {"n": 0, "signs": []}
    real = _rcr_mod.fold_moments_to_radial_characteristic

    def spy(moments, sign):
        calls["n"] += 1
        calls["signs"].append(sign)
        return real(moments, sign)

    monkeypatch.setattr(_rcr_mod, "fold_moments_to_radial_characteristic", spy)
    return calls


def _apply_A_BA(emission: NDArray, sn) -> NDArray:
    r"""BIND POINT — how the extracted single-source A_BA fold is invoked.

    Input: the ℓ=0 iso cell-emission ``(ng, nx)``. Output: the folded ray-source
    ``RadialCharacteristicSourceSink.values`` (cells at μ=±1 per carrying level,
    corners zero). The Step-2 un-weld gives this ONE surface (the loop sites
    2/3/4 inline today). Its production shape is being decided in parallel — the
    main agent flips the ``# BIND:`` line to the chosen surface:
    """
    # Bound (Step 2, operator shape): the extracted single-source A_BA fold is
    # RadialCharacteristicReconstruction.apply — the emission is the ℓ=0 moment
    # (a unit ℓ axis, n_moments=1).
    return _rcr_mod.RadialCharacteristicReconstruction(sn).apply(
        emission[None]).values


def _apply_A_BA_transpose(seed_cotangent: NDArray, sn) -> NDArray:
    r"""BIND POINT — the A_BA Euclidean transpose (ray-cotangent → emission-cotangent),
    IF the chosen A_BA shape exposes a transpose surface. If the user picks a
    factory-only shape (no transpose surface), this gate is subsumed by the
    fold-helper-transpose contract (``test_fold_transpose_euclidean_contract``)
    and the operator-level consistency gate — flag the binding, do not force it.
    """
    # Bound (Step 2, operator shape): wrap the flat ray cotangent as a
    # RadialCharacteristicFlux and pull it back through the reconstruction's
    # Euclidean transpose → the (n_moments=1, ng, nx) bulk-moment cotangent.
    field = RadialCharacteristicFlux(
        values=np.asarray(seed_cotangent, dtype=float),
        space=sn.radial_characteristic_space, mesh=sn)
    return _rcr_mod.RadialCharacteristicReconstruction(sn).apply_transpose(field)


def _wrap_extracted_A_BA(monkeypatch) -> dict:
    r"""Mode-11 sentinel on the EXTRACTED single-source A_BA surface (post-un-weld).
    The wrap MUST sit on the SAME object the production S/F seed arms construct —
    the factory or the operator the un-weld routes them through. xfail until the
    bind is chosen; flip the ``# BIND:`` monkeypatch target:
    """
    calls: dict = {"n": 0}
    # Bound (Step 2, operator shape): wrap the operator method on the CLASS, so
    # every S/F seed arm that constructs a RadialCharacteristicReconstruction and
    # calls ``.apply`` is counted — immune to import binding (it patches the
    # class object all instances share).
    real = _rcr_mod.RadialCharacteristicReconstruction.apply

    def spy(self, moments, /):
        calls["n"] += 1
        return real(self, moments)

    monkeypatch.setattr(
        _rcr_mod.RadialCharacteristicReconstruction, "apply", spy)
    return calls


class TestA_BA_SchurFold:
    r"""``A_BA`` — the ψ½ Schur fold FACTOR (the ``RadialCharacteristicReconstruction``
    single source the ℓ-moment emission folds through at μ = ±1).

    Post-LIFT (campaign step 4c, commit 1) this class pins the FOLD FACTOR + the
    non-carrying control: the fold contract (``P_ℓ(±1) = (±1)^ℓ`` on a manufactured
    anisotropic input), the fold-transpose Euclidean contract, the fold-factor
    surface value + its Euclidean transpose, and the cyl/slab non-carrying control.
    The FULL coupling operator ``RadialCharacteristicEmission`` (``Fold ∘ K_iso ∘
    integrate``), the S/F→pure-bulk lift, and the driver routing are pinned in
    :class:`TestCoupledLift` (the step-2 "S/F EMIT the ray" gates are retired
    there — S/F are now pure bulk).

    Carrying member = **sphere-GL S4 ONLY** (the only geometry that carries a ψ½
    level, R12a; 1 level → 2 fold calls/arm). cylinder/slab are the non-carrying
    CONTROL. Every value row is ≥2G (1G is degenerate, vv anti-#3). Runtime: gates
    raise via ``pytest.fail`` / ``np.testing.assert_*`` (fire under ``python -O``),
    never a bare ``assert`` (vv Mode 8).
    """

    # ── Gate 1: the fold contract on a MANUFACTURED anisotropic input ──────

    def test_fold_contract_anisotropic_activates_p_ell(self):
        r"""Load-bearing refutation #3. Manufacture ``moments`` with ℓ=0 AND ℓ=1
        (≥2G, distinct per group) and assert the closed form
        ``Q̄(+1) = ½Q₀ + (3/2)Q₁``, ``Q̄(−1) = ½Q₀ − (3/2)Q₁`` — the ``sign^ℓ``
        line's P₁(±1) = ±1 asymmetry.

        Tooth (in-process, local): a fold that DROPS ``sign^ℓ`` (``coeff`` = the
        same ``(2ℓ+1)/2`` for both signs) reds the anisotropic assertion by
        ``3·|Q₁|`` (measured ≈ 2.7). NECESSITY: the SAME mutated fold on an
        ℓ=0-only input stays green (P₀≡1, ``sign^0=1`` always) — so the production
        S/F arms, which feed ℓ=0 ONLY, are STRUCTURALLY blind to this bug; the
        anisotropic input is what earns the coverage (§0.6 iso-snapshot blindness).
        """
        Q0 = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])       # (ng=2, nx=3)
        Q1 = np.array([[0.7, -0.3, 0.2], [0.1, 0.9, -0.5]])     # distinct per group
        mom = np.stack([Q0, Q1], axis=0)                        # (L=2, ng, nx)
        fold = _rcs_mod.fold_moments_to_radial_characteristic
        np.testing.assert_allclose(
            fold(mom, +1), 0.5 * Q0 + 1.5 * Q1, rtol=0.0, atol=1e-13,
            err_msg="Q̄(+1) ≠ ½Q₀ + (3/2)Q₁ — the (2ℓ+1)/2·(+1)^ℓ fold is wrong.")
        np.testing.assert_allclose(
            fold(mom, -1), 0.5 * Q0 - 1.5 * Q1, rtol=0.0, atol=1e-13,
            err_msg="Q̄(−1) ≠ ½Q₀ − (3/2)Q₁ — the P₁(−1) = −1 sign is dropped/flipped.")

        # Tooth + necessity: a fold that drops sign^ℓ.
        def _fold_drop_sign(moments, sign):
            moments = np.asarray(moments)
            ell = np.arange(moments.shape[0])
            coeff = (2.0 * ell + 1.0) / 2.0                     # sign^ℓ DROPPED
            return np.tensordot(coeff, moments, axes=(0, 0))

        red_aniso = float(np.max(np.abs(_fold_drop_sign(mom, -1) - fold(mom, -1))))
        red_iso = float(np.max(np.abs(
            _fold_drop_sign(Q0[None], -1) - fold(Q0[None], -1))))
        print(f"  [G1] drop-sign^ℓ: anisotropic RED={red_aniso:.3f}  iso-only={red_iso:.1e}")
        if not red_aniso > 1e-3:
            pytest.fail(f"the drop-sign^ℓ mutation did not red the anisotropic contract "
                        f"({red_aniso:.3e}) — the ℓ≥1 gate is toothless.")
        if not red_iso < 1e-15:
            pytest.fail(f"the iso-only input MOVED under the mutation ({red_iso:.3e}) — "
                        f"the necessity claim is false (iso should be blind to sign^ℓ).")

    # NOTE (campaign step 4c, THE LIFT): the step-2 gates that asserted S/F EMIT
    # the ray via the fold (``test_scattering_and_fission_both_route_through_the_
    # shared_fold``, ``test_scattering_seed_is_the_half_emission_fold``,
    # ``test_fission_seed_is_the_half_emission_fold_fissile``) are RETIRED — S/F
    # are now PURE BULK and the emission moved to the ``RadialCharacteristicEmission``
    # gain (``A_BA``). Their successor coverage lives in :class:`TestCoupledLift`
    # (L1-FWD / L2 / L3), which pins that A_BA emits the fold and S/F carry no ray.

    # ── Gate 5a (LIVE): the fold-helper Euclidean transpose contract ───────

    def test_fold_transpose_euclidean_contract(self):
        r"""The transpose the S-adjoint arm's hard-coded ``0.5`` must be
        single-sourced through. The Euclidean adjoint identity
        ``⟨fold(m, sign), y⟩ = ⟨m, fold_transpose(y, sign)⟩`` with
        ``fold_transpose(y, sign)[ℓ] = ((2ℓ+1)/2)·sign^ℓ · y``, on a MANUFACTURED
        anisotropic ``m`` (ℓ=0 AND ℓ=1) and random ``y``.

        Teeth: (a) a ``0.6`` ℓ=0 coefficient in the reference transpose (≠ the
        fold's ``0.5`` — the scattering.py:1846 hard-code) breaks the identity
        (measured ≈ 0.02–0.04); (b) dropping ``sign^ℓ`` in the transpose breaks
        the sign=−1 leg (the P₁(−1) transpose consistency)."""
        rng = np.random.default_rng(0)
        m = rng.standard_normal((2, 2, 3))                      # (L=2, ng, nx) anisotropic
        y = rng.standard_normal((2, 3))
        fold = _rcs_mod.fold_moments_to_radial_characteristic
        for sign in (-1, +1):
            folded = fold(m, sign)
            lhs = float(np.sum(folded * y))
            rhs = float(np.sum(m * _fold_transpose_reference(y, sign, 2)))
            np.testing.assert_allclose(
                lhs, rhs, rtol=1e-12, atol=1e-12,
                err_msg=f"sign {sign}: ⟨fold(m),y⟩ ≠ ⟨m, fold_transpose(y)⟩ — the "
                        f"fold-transpose is not the Euclidean adjoint of the fold.")
            # Tooth (a): 0.6 ≠ 0.5 at ℓ=0.
            rhs_06 = float(np.sum(m * _fold_transpose_reference(y, sign, 2, coeff0=0.6)))
            d06 = abs(lhs - rhs_06) / (abs(lhs) + abs(rhs_06) + 1e-300)
            if not d06 > 1e-3:
                pytest.fail(f"sign {sign}: the 0.6 ℓ=0-coefficient tooth did not red "
                            f"({d06:.3e}) — the transpose contract is toothless to the "
                            f"scattering.py:1846 hard-coded 0.5.")

        # Tooth (b): the sign^ℓ transpose consistency (P₁(−1) = −1). A transpose
        # that drops sign^ℓ agrees at sign=+1 but breaks at sign=−1.
        ell = np.arange(2)
        for sign in (-1, +1):
            folded = fold(m, sign)
            lhs = float(np.sum(folded * y))
            no_sign = ((2.0 * ell + 1.0) / 2.0)[:, None, None] * y[None]   # sign^ℓ dropped
            d = abs(lhs - float(np.sum(m * no_sign))) / (abs(lhs) + 1e-300)
            print(f"  [G5a] sign={sign}: drop-sign^ℓ transpose defect = {d:.3f}")
            if sign == -1 and not d > 1e-3:
                pytest.fail(f"the drop-sign^ℓ transpose stayed green at sign=−1 "
                            f"({d:.3e}) — the P₁(−1) transpose consistency is unpinned.")

    # NOTE (campaign step 4c, THE LIFT): the step-2 gate
    # ``test_scattering_seed_arm_euclidean_transpose_consistency`` (S's hand-rolled
    # adjoint 0.5 == the forward fold's transpose) is RETIRED — the S-adjoint no
    # longer carries the seed pullback (S is pure bulk). Its stronger successor is
    # :meth:`TestCoupledLift.test_A_BA_scatter_carries_the_seed_pullback_S_carries_none`
    # (L1-ADJ), which pins the pullback ``w·K_isoᵀ(Reconstructionᵀ χ_seed)`` in
    # ``A_BA.apply_transpose`` and its Euclidean fwd↔adj consistency.

    # ── Gate 6 (LIVE): the non-carrying CONTROL (no ray, no fold) ──────────

    def test_non_carrying_control_no_ray_no_fold(self, monkeypatch):
        r"""cylinder + slab are the non-carrying CONTROL (``radial_characteristic_
        space is None`` — refutation #6, NOT "other geometries"): the S/F seed arm
        emits a ``None`` ray and NO fold is invoked. Feed a bulk-only composite
        (the arm is gated on ``psi.radial_characteristic is not None``), assert the
        output ray is ``None`` and the fold spy counter stays 0."""
        cyl = SNMesh(
            Mesh1D(edges=np.linspace(0.05, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CYLINDRICAL, bc_right=BC("vacuum")),
            Quadrature.level_symmetric(4), {0: _mixture(1.0, 0.4, 2)})
        slab = SNMesh(
            Mesh1D(edges=np.linspace(0.0, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CARTESIAN, bc_right=BC("reflective"),
                   bc_left=BC("reflective")),
            Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, 2)})
        for tag, sn in (("cylinder", cyl), ("slab", slab)):
            if sn.radial_characteristic_space is not None:
                pytest.fail(f"{tag} carries a ray space — the non-carrying CONTROL is invalid.")
            N, nx, ng = sn.quad.N, sn.nx, sn.ng
            n_tr = int(sn.angular_trace.layout.total_size)
            bulk_only = FullField(
                interior=AngularFlux.from_mesh(
                    np.random.default_rng(50).standard_normal((N, ng, nx)), sn),
                boundary=AngularBoundaryFlux(values=np.zeros(n_tr), space=sn.angular_trace, mesh=sn),
                radial_characteristic=None)
            calls = _install_fold_spy(monkeypatch)
            out = SNSolver(sn).scattering_op.apply(bulk_only)
            if out.radial_characteristic is not None:
                pytest.fail(f"{tag}: S emitted a non-None ray on a non-carrying mesh "
                            f"(the seed arm must be skipped when the space is None).")
            if calls["n"] != 0:
                pytest.fail(f"{tag}: the fold was invoked {calls['n']}× on a non-carrying "
                            f"mesh — A_BA must not fire without a ray carrier.")

    # ── the fold FACTOR surface (a factor of A_BA, post-lift LIVE) ─────────
    #
    # ``_apply_A_BA`` / ``_apply_A_BA_transpose`` exercise the fold FACTOR
    # (``RadialCharacteristicReconstruction`` — emission → ray) that the coupling
    # operator ``RadialCharacteristicEmission`` composes with ``K_iso ∘ integrate``.
    # The FULL production A_BA (integrate → K_iso → fold) and the driver routing
    # are pinned in :class:`TestCoupledLift` (L1-FWD / L2 / L4-S); these two rows
    # pin just the fold factor's value + Euclidean transpose.

    def test_A_BA_fold_factor_folds_half_emission(self):
        r"""The A_BA fold FACTOR surface
        (``RadialCharacteristicReconstruction.apply``) folds an emission to the
        ray q½ source, matching the closed-form ½·emission loop
        (``_ba_oldloop_reference``)."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        psi = _random_composite(sn, np.random.default_rng(60))
        emission = _s_emission(S, psi)
        got = _apply_A_BA(emission, sn)
        np.testing.assert_array_equal(
            got, _ba_oldloop_reference(emission, sn),
            err_msg="the A_BA fold factor surface ≠ the documented fold loop.")

    def test_A_BA_fold_factor_transpose_euclidean_contract(self):
        r"""The fold FACTOR exposes ``.apply_transpose``, which satisfies the
        Euclidean adjoint contract against the forward fold
        (``⟨fold·emission, y⟩ = ⟨emission, foldᵀ·y⟩``)."""
        sn = _sphere()
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(80)
        emission = rng.standard_normal((sn.ng, sn.nx))
        y = rng.standard_normal(space.shape[0])
        fwd = _apply_A_BA(emission, sn)
        bwd = _apply_A_BA_transpose(y, sn)
        lhs = float(fwd @ y)
        rhs = float(emission.ravel() @ np.asarray(bwd).ravel())
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-11, atol=1e-12,
            err_msg="A_BA.apply_transpose is not the Euclidean adjoint of A_BA.apply.")


# ── Step-4c helpers: THE LIFT — S/F → pure bulk, A_BA born driver-side ──────
#
# The scatter LIFT (campaign step 4c, commit 1): the model-generic S / F gains
# drop their hand-rolled ψ½ seed arm and become PURE BULK; the bulk→ray emission
# is posed as the first-class coupling operator A_BA =
# RadialCharacteristicEmission (= Fold ∘ K_iso ∘ integrate), which the SI/Krylov
# driver lags as its OWN gain (the Wave-O #208 pattern that separated B from S).
# S_bulk ⊕ A_BA is bit-identical to the old monolithic S.apply.


def _a_ba_scatter(sn) -> RadialCharacteristicEmission:
    r"""The production scattering bulk→ray coupling A_BA — the SI driver's OWN
    lagged gain (``RadialCharacteristicEmission`` over S's isotropic kernel, the
    SAME shared object the bulk scatter gain uses — single-sourced emission)."""
    return RadialCharacteristicEmission(sn, SNSolver(sn).scattering_op.isotropic_kernel)


def _pullback_reconstruction(sn, S, chi_seed_values: NDArray) -> NDArray:
    r"""The seed pullback ``w·K_isoᵀ(Reconstructionᵀ χ_seed)`` rebuilt from its
    NAMED factors — the structural decomposition
    ``A_BAᵀ = (∫dμ)ᵀ ∘ K_isoᵀ ∘ Foldᵀ`` the S-adjoint carried inline before the
    LIFT. Shape ``(N, ng, nx)``. (This shares A_BA's fold + kernel objects, so the
    ``== A_BA.apply_transpose.interior`` check is INHERITANCE — necessary-not-
    sufficient; the load-bearing structural cross-check is the fwd↔adj Euclidean
    reciprocity, leg (c) of L1-ADJ.)"""
    space = sn.radial_characteristic_space
    fold = RadialCharacteristicReconstruction(sn, n_moments=1)
    chi_ray = RadialCharacteristicFlux(
        values=np.asarray(chi_seed_values, dtype=float), space=space, mesh=sn)
    m_bar = fold.apply_transpose(chi_ray)                       # Foldᵀ: (1, ng, nx)
    phi0_bar = np.asarray(S.isotropic_kernel.apply_transpose(m_bar[0]))  # K_isoᵀ: (ng, nx)
    w = np.asarray(sn.quad.weights, dtype=float)               # (∫dμ)ᵀ = ×w_n
    return w.reshape((w.size, 1, 1)) * phi0_bar[None]          # (N, ng, nx)


def _fold_half_to(coeff0: float):
    r"""A fold that overrides the ℓ=0 coefficient (½ → ``coeff0``) — the L2 value
    tooth (the un-weld's single-sourcing eliminated the hand-rolled ``0.5``)."""

    def _fold(moments, sign):
        moments = np.asarray(moments)
        ell = np.arange(moments.shape[0])
        coeff = ((2.0 * ell + 1.0) / 2.0) * np.float64(sign) ** ell
        coeff = coeff.copy()
        coeff[0] = coeff0
        return np.tensordot(coeff, moments, axes=(0, 0))

    return _fold


class TestCoupledLift:
    r"""Step 4c, commit 1 — THE SCATTER LIFT (S/F → pure bulk, A_BA born
    driver-side, BIT-IDENTICAL). The successor of the retired step-2 "S/F emit the
    ray" gates.

    Carrying member = **sphere-GL S4 ONLY** (R12a); ≥2G every value row; nonzero
    seed + bulk where the term needs activating (§0.6). Runtime: gates raise via
    :func:`pytest.fail` / ``np.testing.assert_*`` (fire under ``python -O``), never
    a bare ``assert`` (vv Mode 8). Every gate carries a mutation-verified tooth
    (a companion ``*_has_teeth`` / ``*_pins_the_object`` row).

    Commit 1 = the SCATTER lift (S/F ``.apply`` pure bulk, A_BA_scatter its own
    driver gain). **Commit 2 = the FISSION outer-seam migration**: the eigenvalue F
    ray seed moved from the per-ordinate ``from_isotropic → from_angular_source``
    round-trip to the DIRECT moments-fold
    :func:`~orpheus.sn.solver._radial_characteristic_fission_seed` (=
    ``A_BA_fission = Fold ∘ F.kernel``, factored — the outer loop already applied
    ``F.kernel ∘ integrate`` to build ``fission_source``, so only the Fold remains).
    The **L1-F** value gate + the **L4-F** Mode-11 sentinel pin it.
    (``from_angular_source`` STAYS for its other consumers — the final total-source
    reconstruction and the fixed-source external source.)
    """

    # ── L1-FWD: S / F pure bulk, A_BA emits (LIFT deliverable 1, forward) ───

    def test_L1_fwd_S_and_F_apply_are_pure_bulk_A_BA_scatter_emits(self):
        r"""S / F ``.apply`` are PURE BULK — ray present-zero, bulk UNCHANGED
        (``S.apply(psi).interior == S.apply(psi.interior)``, the FullField arm's bulk IS
        the model-generic scatter — the LIFT touched no bulk). The lifted A_BA
        gain carries the emission (ray nonzero, bulk present-zero — the disjoint
        direct sum ``S_bulk ⊕ A_BA``). Positive (A_BA emits) + pure-bulk
        (S/F ray present-zero); ≥2G, nonzero seed + bulk."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        psi = _random_composite(sn, np.random.default_rng(100))
        s_out = S.apply(psi)
        if s_out.radial_characteristic is None:
            pytest.fail("S emitted a None ray on a carrying composite — the composite "
                        "presence law raises under S ⊕ A_BA (must be present-zero).")
        np.testing.assert_array_equal(
            s_out.radial_characteristic.values, 0.0,
            err_msg="S.apply ray ≠ present-zero — the model-generic scatter still "
                    "emits the ψ½ ray (the LIFT did not make S pure bulk).")
        np.testing.assert_array_equal(
            s_out.interior.values, S.apply(psi.interior).values,
            err_msg="S.apply(FullField).interior ≠ S.apply(bulk) — the LIFT altered the "
                    "model-generic scatter bulk (it must touch ONLY the ray).")
        # F pure bulk (fissile sphere; the F-fwd ray arm is dead — the fission
        # emission rides from_angular_source / commit 2, not F.apply).
        snf = _fissile_sphere()
        f_out = SNSolver(snf).fission_op.apply(
            _random_composite(snf, np.random.default_rng(101)))
        if f_out.radial_characteristic is not None:
            np.testing.assert_array_equal(
                f_out.radial_characteristic.values, 0.0,
                err_msg="F.apply ray ≠ present-zero — the F-fwd ray arm is not dropped.")
        # A_BA carries the emission on System B's OWN carrier (B.2b re-type):
        # the "bulk present-zero" disjointness is now STRUCTURAL — the codomain
        # has no bulk slot (Pattern 4) — so the pins are the container types
        # (G-b3.1 (ii)) + the nonzero folded cells; the boundary member is a
        # REAL zero (the fold writes cells only).
        a_out = _a_ba_scatter(sn).apply(psi)
        if type(a_out) is not RadialCharacteristicComposite:
            pytest.fail(f"A_BA.apply returned {type(a_out).__name__} — the block "
                        f"must emit System B's RadialCharacteristicComposite.")
        if type(a_out.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"A_BA interior member is {type(a_out.interior).__name__} "
                        f"— the emission must carry the SOURCE pair.")
        if type(a_out.boundary) is not RadialCharacteristicBoundarySourceSink:
            pytest.fail(f"A_BA boundary member is {type(a_out.boundary).__name__} "
                        f"— the emission must carry the SOURCE pair.")
        if not np.max(np.abs(a_out.to_unified().values)) > 1e-6:
            pytest.fail("A_BA.apply ray ≈ 0 — the lifted gain does not carry the "
                        "bulk→ray emission (the coupling is unwired).")
        np.testing.assert_array_equal(
            a_out.boundary.values, 0.0,
            err_msg="A_BA.apply boundary member ≠ 0 — the fold writes CELLS only; "
                    "the corner datum is the boundary's job.")

    # ── L1-ADJ: the DECISIVE lost-pullback catcher (LIFT deliv 1, adjoint) ──

    def test_L1_adj_A_BA_scatter_carries_the_seed_pullback_S_carries_none(self):
        r"""THE DECISIVE pullback catcher (refutation R2). The S-adjoint's
        ``w·K_isoᵀ(Reconstructionᵀ χ_seed)`` bulk pullback moved to
        ``A_BA.apply_transpose``; ``S.apply_transpose`` is now pure bulk. On a
        NONZERO seed cotangent χ (the previously-nulled input — a present-zero χ
        gives ``Reconstructionᵀ(0) = 0`` so every ``.H`` reciprocity gate is BLIND
        to a lost pullback; ONLY a nonzero χ catches it):

        (a) ``A_BA.apply_transpose(χ).interior`` == the pullback (lives in A_BA), and
        is NONZERO (non-vacuity);
        (b) ``S.apply_transpose(χ_seed).interior`` == 0 (S dropped it — pure bulk);
        (c) structurally-independent fwd↔adj Euclidean reciprocity
        ``⟨A_BA·ψ, χ_seed⟩ = ⟨ψ, A_BAᵀ·χ⟩`` (a corrupted pullback lands on ONE
        side — the load-bearing cross-check; (a) is only an INHERITANCE decomposition).

        Tooth: :meth:`test_L1_adj_pullback_catcher_has_teeth`."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        A_BA = _a_ba_scatter(sn)
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(110)
        chi_seed = rng.standard_normal(space.shape[0])
        # B.2b: the block's cotangent is System B's OWN carrier (a FullField
        # seed-only cotangent is now unspellable at the block boundary).
        chi_cot = _ray_composite(sn, chi_seed)
        # (a) the pullback lives in A_BA (== the named-factor reconstruction),
        # and the output is the forward-looking 2-block System-A shape
        # (radial_characteristic is None — G-b3.2's shape pin).
        adj_out = A_BA.apply_transpose(chi_cot)
        if adj_out.radial_characteristic is not None:
            pytest.fail("A_BA.apply_transpose returned a 3-slot System-A cotangent "
                        "— the block speaks the 2-block shape (ray=None); the "
                        "present-zero re-pad is the ADAPTER's job.")
        adj_bulk = adj_out.interior.values
        np.testing.assert_array_equal(
            adj_bulk, _pullback_reconstruction(sn, S, chi_seed),
            err_msg="A_BA.apply_transpose.interior ≠ w·K_isoᵀ(Reconstructionᵀ χ_seed) — "
                    "the lifted seed pullback is wrong.")
        if not np.max(np.abs(adj_bulk)) > 1e-6:
            pytest.fail("A_BA.apply_transpose bulk ≈ 0 on a nonzero seed cotangent — "
                        "the pullback catcher is vacuous.")
        # (b) S dropped the pullback (pure-bulk transpose: seed-only χ → zero
        # bulk). S is FullField-typed — feed it the seed-only FullField.
        s_adj = S.apply_transpose(_seed_composite(sn, chi_seed))
        np.testing.assert_array_equal(
            s_adj.interior.values, 0.0,
            err_msg="S.apply_transpose(seed-only χ).interior ≠ 0 — S still carries the "
                    "seed pullback (the LIFT did not move it to A_BA).")
        # (c) structurally-independent Euclidean fwd↔adj reciprocity, NONZERO seed.
        psi = _random_composite(sn, rng)
        lhs = float(A_BA.apply(psi).to_unified().values @ chi_seed)
        rhs = float(psi.interior.values.ravel() @ adj_bulk.ravel())
        defect = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
        if not defect < 1e-11:
            pytest.fail(f"A_BA Euclidean fwd↔adj reciprocity defect {defect:.3e} ≥ "
                        f"1e-11 — apply_transpose is not the transpose of apply.")

    def test_L1_adj_pullback_catcher_has_teeth(self, monkeypatch):
        r"""TOOTH for L1-ADJ: a corrupted ``A_BA.apply_transpose`` (drop the bulk
        pullback — return present-zero bulk) breaks the fwd↔adj Euclidean
        reciprocity O(1). Proves L1-ADJ's leg (c) catches a lost pullback (and, via
        the composite ``.H``, the reciprocity-leaf sphere row too — see
        ``test_g_adjoint_reciprocity``). Monkeypatch-revert; ``-O``-safe."""
        sn = _sphere()
        A_BA = _a_ba_scatter(sn)
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(111)
        chi_seed = rng.standard_normal(space.shape[0])
        psi = _random_composite(sn, rng)

        from orpheus.transport.full_field import FullField as _FF
        from orpheus.transport.source_sinks import (
            AngularBoundarySourceSink, AngularSourceSink,
        )

        def _drop_pullback(self, cotangent, /):
            # A zero bulk (the pullback dropped) in the b3 block shape
            # (2-block System-A cotangent, ray=None) — the symmetric-drop bug.
            return _FF(
                interior=AngularSourceSink.zeros_on(self.sn_mesh),
                boundary=AngularBoundarySourceSink.zeros_on(self.sn_mesh),
                radial_characteristic=None)

        monkeypatch.setattr(RadialCharacteristicEmission, "apply_transpose", _drop_pullback)
        adj_bulk = A_BA.apply_transpose(_ray_composite(sn, chi_seed)).interior.values
        lhs = float(A_BA.apply(psi).to_unified().values @ chi_seed)
        rhs = float(psi.interior.values.ravel() @ adj_bulk.ravel())
        defect = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
        print(f"  [L1-ADJ tooth] dropped-pullback reciprocity defect = {defect:.3f}")
        if not defect > 1e-3:
            pytest.fail(f"the dropped-pullback mutation left the reciprocity defect at "
                        f"{defect:.3e} — L1-ADJ has no teeth to a lost pullback.")

    # ── L2: A_BA = Fold ∘ K_iso ∘ integrate, single-source fold (deliv 2) ───

    def test_L2_A_BA_scatter_routes_through_the_shared_fold(self, monkeypatch):
        r"""A_BA folds through the SINGLE-SOURCE reconstruction: a Mode-11 wrap on
        the shared ``fold_moments_to_radial_characteristic`` fires EXACTLY
        ``2·n_levels`` (2 signs/level; sphere-GL S4 → 1 level → 2), and the ray
        value is the documented fold loop (``_ba_oldloop_reference``,
        ``array_equal``, inheriting Gate 1's value). The lifted S no longer folds
        (0× — the emission left it).

        Tooth (½ → 0.6 fold coefficient): :meth:`test_L2_fold_value_has_teeth`."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        A_BA = _a_ba_scatter(sn)
        psi = _random_composite(sn, np.random.default_rng(120))
        emission = _s_emission(S, psi)
        n_levels = len(sn.radial_characteristic_space.levels)   # == 1 (sphere-GL S4)
        fold_calls = _install_fold_spy(monkeypatch)
        recon_calls = _wrap_extracted_A_BA(monkeypatch)
        ray = A_BA.apply(psi).to_unified().values
        # (i) the extracted single-source RECONSTRUCTION operator is on the call
        # path (once per A_BA.apply — a re-inlined fold copy leaves this at 0).
        if recon_calls["n"] != 1:
            pytest.fail(f"A_BA routed through the reconstruction operator {recon_calls['n']}× "
                        f"(expected 1) — it does not use the extracted single-source surface.")
        # (ii) the shared inner fold fires 2·n_levels (2 signs/level).
        if fold_calls["n"] != 2 * n_levels:
            pytest.fail(f"A_BA folded {fold_calls['n']}× (expected 2·n_levels = "
                        f"{2 * n_levels}) — it does not route through the single-source fold.")
        np.testing.assert_array_equal(
            ray, _ba_oldloop_reference(emission, sn),
            err_msg="A_BA.apply ray ≠ the documented fold loop (the fold value diverged).")
        # S no longer folds (pure bulk).
        fold_calls["n"] = 0
        S.apply(psi)
        if fold_calls["n"] != 0:
            pytest.fail(f"S folded {fold_calls['n']}× — the lifted S must be pure bulk "
                        f"(the emission moved to A_BA).")

    def test_L2_fold_value_has_teeth(self, monkeypatch):
        r"""TOOTH for L2: a ½ → 0.6 fold coefficient (as the RECONSTRUCTION sees it)
        moves ``A_BA.apply`` off the documented loop (which uses the numerics fold,
        unpatched) — the ``array_equal`` reds. Proves the fold VALUE is pinned."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        A_BA = _a_ba_scatter(sn)
        psi = _random_composite(sn, np.random.default_rng(121))
        emission = _s_emission(S, psi)
        monkeypatch.setattr(_rcr_mod, "fold_moments_to_radial_characteristic", _fold_half_to(0.6))
        ray = A_BA.apply(psi).to_unified().values
        oracle = _ba_oldloop_reference(emission, sn)   # numerics fold, unpatched (0.5)
        red = float(np.max(np.abs(ray - oracle)))
        print(f"  [L2 tooth] ½→0.6 fold: |A_BA.ray − oracle| = {red:.4f}")
        if not red > 1e-3:
            pytest.fail(f"the ½→0.6 fold mutation left A_BA.ray on the oracle "
                        f"({red:.3e}) — the L2 fold-value gate is toothless.")

    # ── L3: S_bulk ⊕ A_BA_ray reconstructs the monolith (deliv 3, Mode-12) ──

    def test_L3_S_bulk_plus_A_BA_ray_reconstructs_the_monolith(self):
        r"""The disjoint direct sum ``S_bulk ⊕ A_BA_ray`` reconstructs the OLD
        monolithic ``S.apply`` — pin the OBJECT (the exact ray PLACEMENT, never a
        keff/sum proxy; Mode-12 R4). S contributes ONLY the bulk (ray present-zero),
        A_BA ONLY the ray (bulk present-zero); the ray placement is byte-for-byte
        the documented monolith fold (``_ba_oldloop_reference``, ``array_equal``).

        Tooth (a ray permutation preserving the level sum):
        :meth:`test_L3_ray_placement_pins_the_object`."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        A_BA = _a_ba_scatter(sn)
        psi = _random_composite(sn, np.random.default_rng(130))
        emission = _s_emission(S, psi)
        s_out, a_out = S.apply(psi), A_BA.apply(psi)
        # Disjoint direct sum: S no ray; A_BA's "no bulk" is STRUCTURAL since
        # B.2b (the block codomain has no bulk slot — Pattern 4).
        np.testing.assert_array_equal(
            s_out.radial_characteristic.values, 0.0,
            err_msg="S contributes a nonzero ray — the direct sum is not disjoint.")
        # The reconstructed monolith ray == the documented fold loop (exact placement).
        monolith_ray = _ba_oldloop_reference(emission, sn)
        if not np.max(np.abs(monolith_ray)) > 1e-6:
            pytest.fail("the reconstructed monolith ray ≈ 0 — L3 is vacuous.")
        np.testing.assert_array_equal(
            a_out.to_unified().values, monolith_ray,
            err_msg="S_bulk ⊕ A_BA_ray ≠ the monolithic S.apply ray (the ray "
                    "placement drifted — a permutation the OBJECT pin catches).")

    def test_L3_ray_placement_pins_the_object(self, monkeypatch):
        r"""TOOTH for L3 (Mode-12): a fold that PERMUTES the cells within a level
        (a radial roll — preserving the per-level SUM) reds the ``array_equal`` ray
        placement while a sum proxy would stay green. Proves L3 pins the OBJECT."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        A_BA = _a_ba_scatter(sn)
        psi = _random_composite(sn, np.random.default_rng(131))
        emission = _s_emission(S, psi)
        real_fold = _rcs_mod.fold_moments_to_radial_characteristic

        def _fold_rolled(moments, sign):
            return np.roll(real_fold(moments, sign), 1, axis=-1)  # permute cells; sum preserved

        monkeypatch.setattr(_rcr_mod, "fold_moments_to_radial_characteristic", _fold_rolled)
        ray = A_BA.apply(psi).to_unified().values
        monolith_ray = _ba_oldloop_reference(emission, sn)   # numerics fold, unpatched
        # The OBJECT (placement) differs …
        placement_red = float(np.max(np.abs(ray - monolith_ray)))
        # … while the per-(level,sign) SUM proxy is BLIND to the permutation.
        space = sn.radial_characteristic_space
        sum_ray = sum(float(space.cells_view(ray, lv, s).sum())
                      for lv in space.levels for s in (-1, +1))
        sum_ref = sum(float(space.cells_view(monolith_ray, lv, s).sum())
                      for lv in space.levels for s in (-1, +1))
        sum_gap = abs(sum_ray - sum_ref) / (abs(sum_ref) + 1e-300)
        print(f"  [L3 tooth] permutation: placement RED={placement_red:.4f}  sum-proxy gap={sum_gap:.1e}")
        if not placement_red > 1e-3:
            pytest.fail(f"the ray permutation did not red the OBJECT pin "
                        f"({placement_red:.3e}) — L3 is a sum proxy, not an OBJECT pin.")
        if not sum_gap < 1e-9:
            pytest.fail(f"the permutation did NOT preserve the sum ({sum_gap:.3e}) — "
                        f"the Mode-12 necessity (a sum proxy would be blind) is not shown.")

    # ── L3.5 (G-b3.3): the transient adapter — byte-identity + DELEGATION ───

    def test_L35_adapter_is_byte_identical_and_delegates(self, monkeypatch):
        r"""G-b3.3 — the two facts that keep production honest through B.2b:

        (i) BYTE-IDENTITY: the FullField-gain adapter reproduces the pre-B.2b
        embedded output exactly — ray == the documented fold loop
        (``_ba_oldloop_reference``), bulk/boundary present-zero (whole-flat
        ``array_equal``: value AND placement).

        (ii) DELEGATION (THE SENTINEL-TEETH CHECK): a counting spy on
        ``RadialCharacteristicEmission.apply`` (the CLASS method) fires
        EXACTLY once per ``adapter.apply`` — proving the adapter DELEGATES to
        the wrapped block rather than inlining a fold copy. This is the
        load-bearing assumption that keeps the L4-S Mode-11 sentinel
        non-vacuous through the adapter (an inlining adapter would leave the
        L4-S spy at 0 while production stayed green)."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        adapter = _RayEmissionFullFieldGain(_a_ba_scatter(sn))
        psi = _random_composite(sn, np.random.default_rng(140))
        emission = _s_emission(S, psi)
        out = adapter.apply(psi)
        np.testing.assert_array_equal(
            out.radial_characteristic.values, _ba_oldloop_reference(emission, sn),
            err_msg="adapter ray ≠ the documented fold loop — the transient embed "
                    "is not byte-identical to the pre-B.2b output.")
        np.testing.assert_array_equal(
            out.interior.values, 0.0,
            err_msg="adapter bulk ≠ present-zero — the embed leaked into the bulk.")
        np.testing.assert_array_equal(
            out.boundary.values, 0.0,
            err_msg="adapter trace ≠ present-zero — the embed leaked into the trace.")
        # (ii) the delegation proof — one isolated call, counter == 1.
        counter = {"n": 0}
        real = RadialCharacteristicEmission.apply

        def spy(self, p, /):
            counter["n"] += 1
            return real(self, p)

        monkeypatch.setattr(RadialCharacteristicEmission, "apply", spy)
        adapter.apply(psi)
        if counter["n"] != 1:
            pytest.fail(f"adapter.apply fired the block's class method {counter['n']}× "
                        f"(expected exactly 1) — the adapter does not DELEGATE, so "
                        f"the L4-S sentinel is Mode-11-vacuous through it.")

    # ── L4-S: the driver routes A_BA as its OWN lagged gain (deliv 3/d) ─────

    def test_L4S_driver_routes_A_BA_scatter_as_its_own_gain(self, monkeypatch):
        r"""THE DECISIVE own-slot driver-routing sentinel (refutation R5: a green
        bulk gate measures the UNCHANGED sibling; only wrapping the NEW A_BA.apply
        proves the driver rewired). Wrap ``RadialCharacteristicEmission.apply`` and
        run a REAL within-group sphere ``solve_fixed_source``: the counter fires
        (> 0). Structural pair: ``_lagged_gains`` carries an A_BA on the sphere and
        NOT on a seedless slab.

        Tooth (drop A_BA from ``_lagged_gains``): :meth:`test_L4S_sentinel_has_teeth`."""
        sn = _sphere()
        _, S, B = _within_group_triple(SNSolver(sn))
        # B.2b: the gain slot carries the ADAPTER wrapping the block — the
        # wraps-predicate reads the wrapped ``_emission``.
        if not any(isinstance(getattr(g, "_emission", None), RadialCharacteristicEmission)
                   for g in _lagged_gains(S, B, sn)):
            pytest.fail("_lagged_gains(sphere) carries no adapter-wrapped "
                        "RadialCharacteristicEmission — A_BA is not wired as its "
                        "own lagged gain.")
        slab = SNMesh(
            Mesh1D(edges=np.linspace(0.0, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CARTESIAN, bc_right=BC("reflective"),
                   bc_left=BC("reflective")),
            Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, 2)})
        _, Ss, Bs = _within_group_triple(SNSolver(slab))
        if any(isinstance(getattr(g, "_emission", None), RadialCharacteristicEmission)
               for g in _lagged_gains(Ss, Bs, slab)):
            pytest.fail("_lagged_gains(slab) carries an A_BA — a seedless mesh has no "
                        "bulk→ray coupling.")
        # Mode-11 sentinel: a REAL within-group sphere solve APPLIES A_BA.
        counter = {"n": 0}
        real = RadialCharacteristicEmission.apply

        def spy(self, psi, /):
            counter["n"] += 1
            return real(self, psi)

        monkeypatch.setattr(RadialCharacteristicEmission, "apply", spy)
        SNSolver(sn).solve_fixed_source(
            np.ones((sn.ng, sn.nx)), np.ones((sn.ng, sn.nx)))
        if counter["n"] <= 0:
            pytest.fail("a real within-group sphere solve did NOT apply A_BA — the SI "
                        "driver does not route the lifted gain (the rewire is missing).")

    def test_L4S_sentinel_has_teeth(self, monkeypatch):
        r"""TOOTH for L4-S: a driver that forgot to widen its gains (``_lagged_gains``
        returns just ``(S, B)`` — the pre-lift shape) leaves the A_BA.apply counter
        at 0. Proves the L4-S sentinel catches an unwired driver (Mode 11)."""
        sn = _sphere()

        def _no_a_ba(S, B, sn_mesh):
            return (S, B)                                    # the un-widened gain tuple

        monkeypatch.setattr(_solver_mod, "_lagged_gains", _no_a_ba)
        counter = {"n": 0}
        real = RadialCharacteristicEmission.apply

        def spy(self, psi, /):
            counter["n"] += 1
            return real(self, psi)

        monkeypatch.setattr(RadialCharacteristicEmission, "apply", spy)
        SNSolver(sn).solve_fixed_source(
            np.ones((sn.ng, sn.nx)), np.ones((sn.ng, sn.nx)))
        print(f"  [L4-S tooth] un-widened _lagged_gains: A_BA.apply calls = {counter['n']}")
        if counter["n"] != 0:
            pytest.fail(f"the un-widened driver STILL applied A_BA {counter['n']}× — "
                        f"the L4-S sentinel would not catch a missing rewire.")

    # ── L1-F: the eigenvalue F seed IS the direct moments-fold (commit 2) ───

    def test_L1F_fission_seed_is_the_moments_fold(self):
        r"""[L1-F value gate, commit 2] The eigenvalue fission ray seed
        (:func:`~orpheus.sn.solver._radial_characteristic_fission_seed`) IS the
        direct ℓ=0 moments-fold of the fission emission ``χ·νΣf·φ`` — pinning the
        migration ``A_BA_fission = Fold ∘ F.kernel`` (factored: the outer loop
        already applied ``F.kernel ∘ integrate`` to build ``fission_source``, so
        only the Fold remains).

        (a) BIT-IDENTICAL to the documented fold loop (``_ba_oldloop_reference`` —
        both fold ``emission[None]`` through ``RadialCharacteristicReconstruction``);
        (b) PRINCIPLED-EQUIV (~ULP) to the RETIRED per-ordinate route
        ``from_isotropic → from_angular_source`` — documenting the migration is a
        genuine ~ULP RE-baseline (the removed round-trip's per-ordinate ``·w``
        reassociates; vv §Bit-identity, criterion-3 FP-non-associativity), NOT
        bit-identical. FISSILE sphere, ≥2G (refutation #4: a non-fissile mixture
        makes the emission — hence the seed — identically zero, a VACUOUS gate).

        Tooth (½ → 0.6 fold): :meth:`test_L1F_fission_seed_value_has_teeth`."""
        from orpheus.sn.solver import (
            _radial_characteristic_fission_seed,
            _radial_characteristic_source_from_per_ordinate,
        )
        snf = _fissile_sphere()
        F = SNSolver(snf).fission_op
        psi = _random_composite(snf, np.random.default_rng(140))
        emission = _f_emission(F, psi)
        if not np.max(np.abs(emission)) > 1e-6:
            pytest.fail("the fission emission ≈ 0 — the L1-F value gate is VACUOUS "
                        "(the mixture is not actually fissile / νΣf = 0; vv #4).")
        got = _radial_characteristic_fission_seed(emission, snf).values
        # (a) bit-identical to the direct moments-fold loop.
        np.testing.assert_array_equal(
            got, _ba_oldloop_reference(emission, snf),
            err_msg="_radial_characteristic_fission_seed ≠ the direct moments-fold "
                    "(_ba_oldloop_reference) — the F seed is not the ℓ=0 fold.")
        # (b) principled-equiv (~ULP) to the RETIRED per-ordinate round-trip.
        old = _radial_characteristic_source_from_per_ordinate(
            AngularSourceSink.from_isotropic(emission, snf).values, snf).values
        # Budget: measured ~9-16 ULP on this config (maxabs 5.6e-16, maxrel 1.9e-15);
        # nulp=32 gives headroom for the removed per-ordinate ·w reassociation (the
        # migration is principled-equiv — vv §Bit-identity criterion 3, NOT byte-id).
        np.testing.assert_array_almost_equal_nulp(got, old, nulp=32)

    def test_L1F_fission_seed_value_has_teeth(self, monkeypatch):
        r"""TOOTH for L1-F (a): a ½ → 0.6 fold coefficient (as the RECONSTRUCTION
        sees it) moves ``_radial_characteristic_fission_seed`` off the documented
        loop (``_ba_oldloop_reference`` uses the numerics fold, unpatched) — the
        ``array_equal`` reds. Proves the F seed's fold VALUE is pinned (a
        dropped/wrong ½ coefficient in the migrated fold is caught)."""
        from orpheus.sn.solver import _radial_characteristic_fission_seed
        snf = _fissile_sphere()
        F = SNSolver(snf).fission_op
        psi = _random_composite(snf, np.random.default_rng(141))
        emission = _f_emission(F, psi)
        monkeypatch.setattr(
            _rcr_mod, "fold_moments_to_radial_characteristic", _fold_half_to(0.6))
        got = _radial_characteristic_fission_seed(emission, snf).values
        oracle = _ba_oldloop_reference(emission, snf)   # numerics fold, unpatched (0.5)
        red = float(np.max(np.abs(got - oracle)))
        print(f"  [L1-F tooth] ½→0.6 fold: |seed − oracle| = {red:.4f}")
        if not red > 1e-3:
            pytest.fail(f"the ½→0.6 fold mutation left the F seed on the oracle "
                        f"({red:.3e}) — the L1-F value gate is toothless.")

    # ── L4-F: the OUTER eigenvalue loop routes F through the fold (deliv/d) ──

    def test_L4F_outer_fission_seed_routes_through_the_moments_fold(self, monkeypatch):
        r"""L4-F — THE DECISIVE Mode-11 sentinel (refutation R5): a GREEN eigenvalue
        solve is BLIND to a leftover ``from_angular_source`` route on the F seed;
        only instrumenting the fold ON the fission path proves the OUTER
        fission-source loop routes through the migrated moments-fold.

        ⚠ REFINEMENT of the brief's literal design — a GLOBAL counter on
        ``RadialCharacteristicReconstruction.apply`` is Mode-11-CONTAMINATED: the
        SCATTER ``A_BA`` (:class:`RadialCharacteristicEmission`, a within-group
        lagged gain) folds through the SAME reconstruction EVERY SI iteration
        (measured: 322 scatter folds even with the F seed reverted), so a global
        ``n > 0`` fires from scatter alone and reverting the F seed would NOT red it.
        The fission-SPECIFIC catcher wraps the seam
        ``_radial_characteristic_fission_seed`` and measures the fold-count DELTA
        *attributable to it* — isolating the fission fold from the scatter fold.

        Structural pair (HAZARD 5): ``_lagged_gains`` carries the SCATTER ``A_BA``
        (over ``S.isotropic_kernel``), NOT a fission fold — F is the OUTER
        ``q_ext``, never a within-group gain.

        Tooth (revert the migration): :meth:`test_L4F_sentinel_has_teeth`."""
        snf = _fissile_sphere()
        # Global fold counter (scatter + fission both fold through Reconstruction).
        fold = {"n": 0}
        real_fold = RadialCharacteristicReconstruction.apply

        def fold_spy(self, moments, /):
            fold["n"] += 1
            return real_fold(self, moments)

        monkeypatch.setattr(RadialCharacteristicReconstruction, "apply", fold_spy)
        # Fission-SPECIFIC: the fold delta DURING _radial_characteristic_fission_seed.
        seam = {"n": 0, "fold_delta": 0}
        real_seed = _solver_mod._radial_characteristic_fission_seed

        def seam_spy(fission_source, sn_mesh):
            before = fold["n"]
            out = real_seed(fission_source, sn_mesh)
            seam["n"] += 1
            seam["fold_delta"] += fold["n"] - before
            return out

        monkeypatch.setattr(
            _solver_mod, "_radial_characteristic_fission_seed", seam_spy)
        # A REAL fissile-sphere eigenvalue solve (default source_iteration → the
        # eig-SI fission-seed site, solver.py:1453).
        _solver_mod.solve_sn(snf.materials, snf.mesh, snf.quad)
        # (i) the OUTER eigenvalue loop calls the fission seam.
        if seam["n"] <= 0:
            pytest.fail("the eigenvalue outer loop never called "
                        "_radial_characteristic_fission_seed — the F seed is not on "
                        "the eigenvalue path (Mode-11: green keff, uncaught).")
        # (ii) each fission-seam call routes through the moments-fold (delta > 0).
        if seam["fold_delta"] <= 0:
            pytest.fail(f"the fission seed fired {seam['n']}× but folded "
                        f"{seam['fold_delta']}× through RadialCharacteristic"
                        f"Reconstruction — a leftover from_angular_source route "
                        f"survives on the F seed (Mode-11).")
        print(f"  [L4-F] eigenvalue solve: seam_n={seam['n']} "
              f"seam_fold_delta={seam['fold_delta']} global_fold={fold['n']}")
        # Structural pair (HAZARD 5): the within-group gains carry the SCATTER A_BA
        # (over S.isotropic_kernel), NOT a fission fold — F is the outer q_ext.
        _, S, B = _within_group_triple(SNSolver(snf))
        gains = _lagged_gains(S, B, snf)
        emissions = [
            g._emission for g in gains
            if isinstance(getattr(g, "_emission", None), RadialCharacteristicEmission)
        ]
        if len(emissions) != 1 or emissions[0].emission_kernel is not S.isotropic_kernel:
            pytest.fail("_lagged_gains does not carry EXACTLY the scatter A_BA (over "
                        "S.isotropic_kernel) — the F fold must be the OUTER q_ext "
                        "seam, never a within-group gain (HAZARD 5).")

    def test_L4F_sentinel_has_teeth(self, monkeypatch):
        r"""TOOTH for L4-F: reverting the migration — pointing
        ``_radial_characteristic_fission_seed`` back to the RETIRED per-ordinate
        BYPASS (``from_isotropic → from_angular_source``, which never touches
        ``RadialCharacteristicReconstruction``) — leaves the fission-path
        ``fold_delta`` at 0, EVEN THOUGH the scatter A_BA keeps the GLOBAL fold
        counter > 0 (measured 322). Proves BOTH (a) the L4-F ``fold_delta`` catcher
        reds a reverted migration, and (b) a GLOBAL fold counter would NOT (it stays
        > 0 from scatter) — the fission-specific delta is WHY the sentinel has teeth
        (Mode 11). This tooth bakes the refutation of the brief's literal design in."""
        from orpheus.sn.solver import _radial_characteristic_source_from_per_ordinate
        snf = _fissile_sphere()
        fold = {"n": 0}
        real_fold = RadialCharacteristicReconstruction.apply

        def fold_spy(self, moments, /):
            fold["n"] += 1
            return real_fold(self, moments)

        monkeypatch.setattr(RadialCharacteristicReconstruction, "apply", fold_spy)
        # The reverted seam: compute the seed via the RETIRED per-ordinate route
        # (no Reconstruction) but record the fold delta the sentinel measures.
        seam = {"n": 0, "fold_delta": 0}

        def bypass_seam(fission_source, sn_mesh):
            before = fold["n"]
            out = _radial_characteristic_source_from_per_ordinate(
                AngularSourceSink.from_isotropic(fission_source, sn_mesh).values,
                sn_mesh)
            seam["n"] += 1
            seam["fold_delta"] += fold["n"] - before
            return out

        monkeypatch.setattr(
            _solver_mod, "_radial_characteristic_fission_seed", bypass_seam)
        _solver_mod.solve_sn(snf.materials, snf.mesh, snf.quad)
        print(f"  [L4-F tooth] reverted seed: global_fold={fold['n']} (scatter) "
              f"seam_n={seam['n']} seam_fold_delta={seam['fold_delta']}")
        # The seam still fires (the bypass IS called) …
        if seam["n"] <= 0:
            pytest.fail("the reverted seam never fired — the tooth is mis-wired.")
        # … but the fission path no longer folds — the fold_delta catcher reds.
        if seam["fold_delta"] != 0:
            pytest.fail(f"the reverted (bypass) fission seed STILL folded "
                        f"{seam['fold_delta']}× — the L4-F fold_delta catcher would "
                        f"not red a reverted migration.")
        # The Mode-11 lesson: a GLOBAL counter would MISS this (scatter keeps it > 0).
        if not fold["n"] > 0:
            pytest.fail(f"the global fold counter is {fold['n']} — the scatter A_BA "
                        f"should keep it > 0 even with F reverted (the Mode-11 point "
                        f"a global counter is blind to; the fission-specific delta "
                        f"is why the sentinel has teeth).")


# ── Step-3 helpers: A_AB the ψ½ seed injection (ray → bulk) ────────────────
#
# A_AB = RadialCharacteristicSeeding: the ray ψ½ seed injected into the bulk
# Morel–Montry angular recurrence. It is CELL-LOCAL ANGULAR (the seed at cell i
# feeds cell i's ordinate recurrence; NO spatial coupling — the radial march is
# A_BB's job), so — unlike A_BB's spatially-woven forward matvec — BOTH
# directions realize HERE as thin WRAPs of the single-sourced closure methods
# (precompute_psi_state / cell_contribution / angular_adjoint). σ-INDEPENDENT:
# with the bulk zeroed the collision/streaming terms drop out, so A_AB needs no
# σ_t (the constructor takes sn_mesh only). The forward .apply's contribution to
# (L+C).apply is isolated by LINEARITY (interior=0, boundary=0 → only the seed's
# angular numerator survives); the transpose is the seed_cells_bar term the
# in-sweep reverse adds on cells(p,-1). A_sb=0 (block-triangular) and A_bs≈7.5
# (this coupling's magnitude) are already pinned by TestRegressionFloor — not
# re-tested here. The sphere carries ONE level (R12a), so the per-level loop is
# length 1: a multi-carrying-level indexing bug is UNTESTABLE with current
# geometry (cylinder is non-carrying) — an inherited blind spot, noted not faked.


def _bulk_composite(sn, bulk_values: NDArray) -> FullField:
    """A composite with the given bulk, zero trace, and a zero ψ½ ray — so
    ``(L+C).apply_transpose`` isolates A_AB's ``seed_cells_bar``: the ray=0 nulls
    the A_BB self-block (``_seed_rows_transpose`` on a zero χ_seed → 0), leaving
    only the M-M thread cotangent on ``cells(p,-1)``."""
    n_tr = int(sn.angular_trace.layout.total_size)
    ns = sn.radial_characteristic_space.shape[0]
    return FullField(
        interior=AngularFlux.from_mesh(bulk_values, sn),
        boundary=AngularBoundaryFlux(
            values=np.zeros(n_tr), space=sn.angular_trace, mesh=sn),
        radial_characteristic=RadialCharacteristicFlux(
            values=np.zeros(ns), space=sn.radial_characteristic_space, mesh=sn),
    )


def _seed_flux(sn, rng) -> RadialCharacteristicFlux:
    """A random ψ½ ray seed — the ``A_AB.apply`` input."""
    ns = sn.radial_characteristic_space.shape[0]
    return RadialCharacteristicFlux(
        values=rng.standard_normal(ns), space=sn.radial_characteristic_space, mesh=sn)


def _bulk_cotangent(sn, rng) -> AngularSourceSink:
    """A random bulk-residual cotangent — the ``A_AB.apply_transpose`` input."""
    return AngularSourceSink.from_mesh(
        rng.standard_normal((sn.quad.N, sn.ng, sn.nx)), sn)


def _install_closure_spy(monkeypatch, sn, method_name: str) -> list[dict]:
    r"""Mode-11 sentinel: WRAP a ``MorelMontryAngularSweep`` method (the shared
    M-M closure kernel A_AB routes through) on the closure CLASS, recording each
    call's ``(args, kwargs)`` and delegating to the real method. Class-level so
    it is robust to the closure's storage layout (an instance-attr patch could
    trip ``__slots__``); the test uses ONE closure, so no cross-instance leak.
    Proves ``apply`` / ``apply_transpose`` EXECUTE the production kernel — a
    divergent inlined copy would leave the list empty (Cardinal Rule 2)."""
    cls = type(sn.pole_angular_closure)
    real = getattr(cls, method_name)

    def spy(self, *args, **kwargs):
        calls.append({"args": args, "kwargs": kwargs})
        return real(self, *args, **kwargs)

    calls: list[dict] = []
    monkeypatch.setattr(cls, method_name, spy)
    return calls


class TestA_AB_SeedInjection:
    r"""``A_AB`` = :class:`RadialCharacteristicSeeding` — the ray→bulk ψ½ seed
    injection (campaign step 3).

    The off-diagonal ``(transport, ray)`` coupling: the ψ½ ray seeds the bulk
    Morel–Montry angular recurrence. CELL-LOCAL ANGULAR (no spatial coupling),
    so both ``apply`` (ray → bulk residual) and ``apply_transpose`` (bulk
    cotangent → ray seed cotangent) are realized as WRAPs of the single-sourced
    closure methods and ``is_adjointable = True`` — the same both-directions
    completeness ``A_BB`` reached at step 4b (its forward is the radial march,
    single-sourced with the walk via ``radial_characteristic_forward_residual``).

    Sphere-GL S4 is the ONLY carrying member (cylinder/slab are the non-carrying
    CONTROL — the constructor rejects them). Every value row is ≥2G. Gates raise
    via :func:`pytest.fail` / ``np.testing.assert_*`` (fire under ``python -O``),
    never a bare ``assert`` (vv Mode 8).

    L13/Mode-11 caveat: the bit-identity gates (:meth:`test_apply_matches_the_in_sweep_seed_injection`,
    :meth:`test_apply_transpose_is_the_in_sweep_seed_cells_bar`) route both
    ``A_AB`` and the ``(L+C)`` reference through the SAME closure methods, so
    they INHERIT bit-identity and are blind to a bug inside a shared method. The
    correctness cross-check is :meth:`test_euclidean_adjoint_consistency`
    (forward ↔ transpose — a shared-method sign bug lands on one side and breaks
    reciprocity)."""

    # ── Forward — the seed injection ≡ the in-sweep contribution ───────────

    def test_apply_matches_the_in_sweep_seed_injection(self, monkeypatch):
        r"""``A_AB.apply(s)`` ≡ the seed's contribution to ``(L+C).apply``
        (interior=0, boundary=0 isolate ``A_bs`` by linearity), BIT-IDENTICALLY
        (``array_equal`` — with ψ_cell=0 the σ-diagonal cancels, so this is
        0-ULP, not principled-equiv). A Mode-11 sentinel proves ``apply``
        EXECUTES the shared closure kernel with the bulk ZEROED and the seed
        passed as ``radial_characteristic``. σ-independence is asserted
        positively (the reference is identical at two ``σ_t`` slopes). ≥2G."""
        sn = _sphere()
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(30)
        sv = rng.standard_normal(space.shape[0])
        seed = RadialCharacteristicFlux(values=sv, space=space, mesh=sn)
        # Reference (real methods, computed BEFORE the spy): the seed's
        # contribution to (L+C).apply — interior=0/boundary=0 ⇒ A_AA·0 = 0, so the
        # bulk output is exactly A_AB·s.
        reference = _loss(sn).apply(_seed_composite(sn, sv)).interior
        ref_other_sigma = _loss(sn, slope=0.9).apply(_seed_composite(sn, sv)).interior
        np.testing.assert_array_equal(
            reference.values, ref_other_sigma.values,
            err_msg="the seed→bulk contribution changed with σ_t — A_AB is not "
                    "σ-independent (the isolation-by-linearity premise is wrong).")
        # Mode-11 sentinel on the shared closure kernel.
        pre = _install_closure_spy(monkeypatch, sn, "precompute_psi_state")
        cc = _install_closure_spy(monkeypatch, sn, "cell_contribution")
        out = RadialCharacteristicSeeding(sn).apply(_member(seed))
        if len(pre) != 1:
            pytest.fail(
                f"precompute_psi_state called {len(pre)}× (expected 1) — A_AB."
                f"apply is not the single-precompute WRAP.")
        psi_view_arg = np.asarray(pre[0]["args"][0])
        if np.max(np.abs(psi_view_arg)) != 0.0:
            pytest.fail(
                "A_AB.apply did NOT zero the bulk psi_view — A_AA's angular "
                "redistribution would leak into the isolated coupling.")
        # B.2c: the engine receives the BRIDGED unified seed (a fresh
        # to_unified object, so the pre-re-type object-identity pin becomes a
        # role + value-fidelity pin — the bridge is bitwise).
        bridged = pre[0]["kwargs"].get("radial_characteristic")
        if type(bridged) is not RadialCharacteristicFlux:
            pytest.fail(
                "A_AB.apply did not pass a role-preserved ψ½ FLUX as "
                f"radial_characteristic; got {type(bridged).__name__}.")
        np.testing.assert_array_equal(
            bridged.values, sv,
            err_msg="the bridged seed is not value-faithful to the input "
                    "composite (the to_unified bridge must be bitwise).")
        if len(cc) < sn.nx:
            pytest.fail(
                f"cell_contribution called {len(cc)}× (< nx = {sn.nx}) — the "
                f"cell-local angular injection did not visit every cell.")
        np.testing.assert_array_equal(
            out.interior.values, reference.values,
            err_msg="A_AB.apply is not bit-identical to the in-sweep injection.")

    # ── Transpose — the seed_cells_bar term ≡ the in-sweep reverse ─────────

    def test_apply_transpose_is_the_in_sweep_seed_cells_bar(self, monkeypatch):
        r"""``A_AB.apply_transpose(v).cells(p,-1)`` ≡ the ``seed_cells_bar`` term
        the in-sweep reverse adds — ``(L+C).apply_transpose(interior=v, ray=0)`` on
        the ray block (ray=0 nulls the A_BB self-block, isolating the M-M thread
        cotangent), BIT-IDENTICALLY. The ``+1`` leg and both corners stay EXACTLY
        0 (the forward writes only the inward leg). A Mode-11 sentinel proves
        ``angular_adjoint`` runs exactly once. ≥2G."""
        sn = _sphere()
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(31)
        vv = rng.standard_normal((sn.quad.N, sn.ng, sn.nx))
        reference = _loss(sn).apply_transpose(
            _bulk_composite(sn, vv)).radial_characteristic
        aa = _install_closure_spy(monkeypatch, sn, "angular_adjoint")
        # B.2c: the cotangent arrives as System A's FullField (the codomain);
        # only its interior member is read (trace/ray structurally discarded).
        out = RadialCharacteristicSeeding(sn).apply_transpose(
            _bulk_composite(sn, vv))
        if len(aa) != 1:
            pytest.fail(
                f"angular_adjoint called {len(aa)}× (expected 1) — A_AB."
                f"apply_transpose is not the single-adjoint WRAP.")
        for p in space.levels:
            np.testing.assert_array_equal(
                out.interior.cells(p, -1), reference.cells(p, -1),
                err_msg=f"level {p}: apply_transpose cells(-1) ≠ the in-sweep "
                        f"seed_cells_bar.")
            np.testing.assert_array_equal(
                out.interior.cells(p, +1), 0.0,
                err_msg=f"level {p}: apply_transpose wrote the +1 leg (must be 0).")
            np.testing.assert_array_equal(
                out.boundary.corner(p, -1), 0.0,
                err_msg=f"level {p}: apply_transpose wrote the -1 corner (be 0).")
            np.testing.assert_array_equal(
                out.boundary.corner(p, +1), 0.0,
                err_msg=f"level {p}: apply_transpose wrote the +1 corner (be 0).")

    # ── Euclidean adjoint consistency — THE correctness cross-check ────────

    def test_euclidean_adjoint_consistency(self):
        r"""``⟨A_AB·u, v⟩ = ⟨u, A_ABᵀ·v⟩`` (Euclidean, plain dot — NOT the
        ``V_cell`` metric, which is the COMPOSITE Hilbert adjoint realized once
        at the coupled operator, L19) to < 1e-11, ≥3 draws. THE load-bearing
        correctness gate: it compares ``apply`` (precompute + cell_contribution)
        to ``apply_transpose`` (angular_adjoint) — two separately-implemented
        duals — so a sign/wiring bug in EITHER shared method lands on ONE side
        and breaks reciprocity (unlike the bit-identity gates, which route both
        sides through the same method and are blind to a shared-method bug).
        ≥2G."""
        sn = _sphere()
        op = RadialCharacteristicSeeding(sn)
        for seed in (1, 2, 3):
            rng = np.random.default_rng(seed)
            u, v = _seed_flux(sn, rng), _bulk_cotangent(sn, rng)
            # B.2c I/O: apply's trace/ray output slots are present-zero, so
            # pairing the interior member against the bulk cotangent IS the
            # full Euclidean dot; the transpose reads a bulk-only FullField.
            lhs = float(
                op.apply(_member(u)).interior.values.ravel() @ v.values.ravel())
            rhs = float(
                u.values.ravel()
                @ op.apply_transpose(
                    _bulk_composite(sn, v.values)).to_unified().values.ravel())
            defect = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
            if not defect < 1e-11:
                pytest.fail(
                    f"seed {seed}: Euclidean adjoint defect {defect:.3e} ≥ 1e-11 "
                    f"— apply_transpose is not the transpose of apply.")

    def test_euclidean_adjoint_consistency_tooth(self, monkeypatch):
        r"""TOOTH for the adjoint gate: flipping the sign of ``cell_contribution``'s
        angular numerator flips ``apply`` (the forward side) but NOT
        ``apply_transpose`` (which routes through ``angular_adjoint``), so
        reciprocity breaks and the defect jumps to O(1). Proves the < 1e-11 gate
        has teeth AND that a shared-method bug DOES surface on the adjoint
        cross-check (the L13 escape hatch the bit-identity gates lack)."""
        sn = _sphere()
        cls = type(sn.pole_angular_closure)
        real = cls.cell_contribution

        def flip(self, *args, **kwargs):
            denom, upstream = real(self, *args, **kwargs)
            return denom, -upstream

        monkeypatch.setattr(cls, "cell_contribution", flip)
        op = RadialCharacteristicSeeding(sn)
        rng = np.random.default_rng(1)
        u, v = _seed_flux(sn, rng), _bulk_cotangent(sn, rng)
        lhs = float(
            op.apply(_member(u)).interior.values.ravel() @ v.values.ravel())
        rhs = float(
            u.values.ravel()
            @ op.apply_transpose(
                _bulk_composite(sn, v.values)).to_unified().values.ravel())
        defect = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
        if not defect > 1e-3:
            pytest.fail(
                f"the cell_contribution sign flip left the adjoint defect at "
                f"{defect:.3e} — the consistency gate has no teeth.")

    # ── Seed-consumed asymmetry — reads ONLY the inward leg ────────────────

    def test_apply_reads_only_the_inward_leg(self):
        r"""``A_AB.apply`` reads ONLY the inward ``cells(p,-1)`` leg of the seed
        (the recurrence seed): two seeds sharing that leg but differing in the
        ``+1`` leg and both corners give IDENTICAL bulk output. Non-vacuity
        (Mode-5): the output is asserted non-trivial (a zero output would satisfy
        the identity vacuously)."""
        sn = _sphere()
        op = RadialCharacteristicSeeding(sn)
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(32)
        full = rng.standard_normal(space.shape[0])
        only_minus = np.zeros_like(full)
        for p in space.levels:
            space.cells_view(only_minus, p, -1)[...] = space.cells_view(full, p, -1)
        if not np.max(np.abs(only_minus)) > 0.0:
            pytest.fail("the inward -1 leg is zero — the asymmetry gate is vacuous.")
        out_full = op.apply(_member(
            RadialCharacteristicFlux(values=full, space=space, mesh=sn)))
        out_minus = op.apply(_member(
            RadialCharacteristicFlux(values=only_minus, space=space, mesh=sn)))
        if not np.max(np.abs(out_full.interior.values)) > 1e-6:
            pytest.fail("A_AB.apply output is ~0 — the asymmetry gate is vacuous.")
        np.testing.assert_array_equal(
            out_full.interior.values, out_minus.interior.values,
            err_msg="A_AB.apply changed when only the +1 leg / corners changed — "
                    "it reads more than the inward starting-direction leg.")

    # ── Non-carrying CONTROL + mesh-identity ───────────────────────────────

    def test_constructor_and_mesh_identity_reject_non_carrying(self):
        r"""The guards, NET-NEW teeth (L4). Non-carrying CONTROL — a cylinder
        (level-symmetric) and a slab (Cartesian) have
        ``radial_characteristic_space is None`` → the seedless guard fires with
        ``match=`` the specific message. Positive control — the sphere
        constructs. Mesh-identity (Pattern 4) — ``apply`` / ``apply_transpose``
        refuse a field on a DIFFERENT sphere."""
        cyl = SNMesh(
            Mesh1D(edges=np.linspace(0.05, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CYLINDRICAL, bc_right=BC("vacuum")),
            Quadrature.level_symmetric(4), {0: _mixture(1.0, 0.4, 2)})
        slab = SNMesh(
            Mesh1D(edges=np.linspace(0.0, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CARTESIAN, bc_right=BC("reflective"),
                   bc_left=BC("reflective")),
            Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, 2)})
        if cyl.radial_characteristic_space is not None:
            pytest.fail("the cylinder carries a ray space — CONTROL invalid.")
        if slab.radial_characteristic_space is not None:
            pytest.fail("the slab carries a ray space — CONTROL invalid.")
        with pytest.raises(ValueError, match="carries no starting-direction ray"):
            RadialCharacteristicSeeding(cyl)
        with pytest.raises(ValueError, match="carries no starting-direction ray"):
            RadialCharacteristicSeeding(slab)
        # Positive control + the mesh-identity guard (a field on ANOTHER sphere).
        sn = _sphere()
        op = RadialCharacteristicSeeding(sn)
        other = _sphere()
        with pytest.raises(ValueError, match="mesh-identity invariant"):
            op.apply(_member(_seed_flux(other, np.random.default_rng(9))))
        with pytest.raises(ValueError, match="mesh-identity invariant"):
            op.apply_transpose(_bulk_composite(
                other,
                np.random.default_rng(9).standard_normal(
                    (other.quad.N, other.ng, other.nx))))
        # B.2c block-boundary refusals: a unified leaf / a bulk field are not
        # the typed carriers (parse-don't-validate at the block boundary).
        with pytest.raises(TypeError, match="System B's member carrier"):
            op.apply(_seed_flux(sn, np.random.default_rng(9)))
        with pytest.raises(TypeError, match="expected a FullField"):
            op.apply_transpose(_bulk_cotangent(sn, np.random.default_rng(9)))

    def test_b2c_grid_entry_containers(self):
        r"""B.2c re-type (G-c1.3): A_AB's grid-entry carriers. ``apply`` emits
        System A's FullField — the interior member carries the bulk term, the
        trace and ψ½ slots are PRESENT-ZERO (the transitional 3-block presence
        convention; ray → ``None`` at the B.2d eviction); ``apply_transpose``
        emits System B's composite (source members). Declared domain/codomain
        pinned by IDENTITY (memo F2 — the composite space object, and
        System A's cached full_field_space)."""
        sn = _sphere()
        op = RadialCharacteristicSeeding(sn)
        rng = np.random.default_rng(11)
        if op.domain is not sn.radial_characteristic_composite_space:
            pytest.fail("A_AB.domain is not THE composite member-space object.")
        if op.codomain is not sn.full_field_space:
            pytest.fail("A_AB.codomain is not THE full_field_space object.")
        out = op.apply(_member(_seed_flux(sn, rng)))
        if type(out) is not FullField:
            pytest.fail(f"apply did not emit a FullField; got {type(out).__name__}.")
        if type(out.interior) is not AngularSourceSink:
            pytest.fail(f"apply's interior is not an AngularSourceSink; got "
                        f"{type(out.interior).__name__}.")
        np.testing.assert_array_equal(
            out.boundary.values, 0.0,
            err_msg="apply's trace slot is not present-zero (A_AB writes only "
                    "the interior).")
        if out.radial_characteristic is None:
            pytest.fail("apply's ψ½ slot is None — the transitional 3-block "
                        "presence convention requires present-zero until B.2d.")
        np.testing.assert_array_equal(
            out.radial_characteristic.values, 0.0,
            err_msg="apply's ψ½ slot is not present-ZERO.")
        out_t = op.apply_transpose(_bulk_composite(
            sn, rng.standard_normal((sn.quad.N, sn.ng, sn.nx))))
        if type(out_t) is not RadialCharacteristicComposite:
            pytest.fail(f"apply_transpose did not emit the member composite; "
                        f"got {type(out_t).__name__}.")
        if type(out_t.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"apply_transpose did not emit SOURCE members; got "
                        f"{type(out_t.interior).__name__}.")


class TestSystemRoleLattice:
    r"""4a — the ``SystemRole {A, B, COUPLED}`` two-system role lattice.

    The COARSE two-system partition that makes System B first-class in the
    operator algebra (orthogonal to :class:`~orpheus.numerics.operator.BlockRole`,
    the within-System-A bulk↔boundary refinement): the self-block ``A_BB`` and
    the ray boundary ``B_b`` are System B; the off-diagonal couplings ``A_AB`` /
    ``A_BA`` span both systems (COUPLED); every model-generic System-A leaf stays
    unclassified (``None``). The join is the two-system analogue of the
    block-role union — ``A ⊔ B = COUPLED``. Foundation: a software-invariant
    gate (no ``verifies`` — it pins no equation).
    """

    def test_join_is_the_two_system_union(self):
        # The defining law A ⊔ B = COUPLED, its symmetry, idempotence, COUPLED
        # absorption, and the conservative None propagation (an operator outside
        # the two-system decomposition stays outside under a sum).
        A, B, C = SystemRole.A, SystemRole.B, SystemRole.COUPLED
        assert _join_system_roles(A, A) is A
        assert _join_system_roles(B, B) is B
        assert _join_system_roles(C, C) is C
        assert _join_system_roles(A, B) is C
        assert _join_system_roles(B, A) is C          # symmetric
        assert _join_system_roles(A, C) is C
        assert _join_system_roles(C, B) is C          # COUPLED absorbs
        assert _join_system_roles(None, A) is None
        assert _join_system_roles(A, None) is None
        assert _join_system_roles(None, None) is None

    def test_psi_half_blocks_carry_their_system_role(self):
        # The four ψ½ blocks are stamped at the class level (the classification
        # is a class attribute — readable without instantiation).
        assert RadialCharacteristicOperator.system_role is SystemRole.B          # A_BB
        assert RadialCharacteristicSeeding.system_role is SystemRole.COUPLED     # A_AB
        assert RadialCharacteristicBoundaryOperator.system_role is SystemRole.B  # B_b
        assert (
            _rcr_mod.RadialCharacteristicReconstruction.system_role
            is SystemRole.COUPLED  # A_BA
        )

    def test_model_generic_operators_stay_unclassified(self):
        # The CONTROL: System-A model-generic leaves carry NO intrinsic
        # two-system membership — an SN context composes them into System A, but
        # they belong to no system by construction (the honest None default).
        assert SNBoundaryOperator.system_role is None      # B_a (System A trace boundary)
        assert MultiplicationOperator.system_role is None  # C (collision)
        assert ScatteringOperator.system_role is None      # S
        assert FissionOperator.system_role is None         # F

    def test_role_propagates_through_the_composers(self):
        # The derivation fires through the composers exactly as block_role does:
        # OperatorSum joins its summands, the G-adjoint preserves. A_AB is
        # σ-independent, so a sphere instance is cheap; it carries COUPLED, so
        # both a sum with itself and its adjoint stay COUPLED.
        a_ab = RadialCharacteristicSeeding(_sphere())
        assert (a_ab + a_ab).system_role is SystemRole.COUPLED   # OperatorSum join
        assert (2.0 * a_ab).system_role is SystemRole.COUPLED    # ScaledOperator passthrough
        assert a_ab.H.system_role is SystemRole.COUPLED          # _AdjointOperator passthrough

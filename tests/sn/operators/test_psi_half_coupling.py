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

import numpy as np
import pytest
from numpy.typing import NDArray

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import (
    RadialCharacteristicBoundaryOperator,
    SNBoundaryOperator,
)
from orpheus.sn.operators.streaming import StreamingOperator
import orpheus.sn.loss_representation as _lr_mod
import orpheus.sn.operators.radial_characteristic as _rc_mod
from orpheus.sn.operators.radial_characteristic import (
    RadialCharacteristicOperator,
    RadialCharacteristicSeeding,
)
import orpheus.numerics.spaces.radial_characteristic_space as _rcs_mod
import orpheus.transport.operators.radial_characteristic_reconstruction as _rcr_mod
from orpheus.sn.solver import SNSolver, _within_group_triple
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
from orpheus.transport.full_field import FullField
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
        bulk=AngularFlux.from_mesh(np.zeros((N, ng, nx)), sn),
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

    def test_bulk_to_ray_coupling_lives_in_the_lagged_scattering_gain(self):
        r"""The full ``A_BA`` (bulk→ray) is carried ONLY by the LAGGED scattering
        gain ``S`` (``S_sb ≠ 0``, the ray scattering source), and the ray cannot
        feed the scalar flux (``S_bs = 0``, ψ½ has zero moment weight). So the
        coupled system's off-diagonal ``A_BA`` is an OUTER-iterated block, never a
        within-sweep one."""
        sn = _sphere()
        _, S, _ = _within_group_triple(SNSolver(sn))
        b, _, s, _ = _blocks(sn)
        Sd = _dense(S.apply, _template(sn))
        s_sb, s_bs, s_bb = _bn(Sd, s, b), _bn(Sd, b, s), _bn(Sd, b, b)
        print(f"  S: S_sb(bulk→ray src)={s_sb:.3e}  S_bs(ray→scalar)={s_bs:.2e}  S_bb={s_bb:.3e}")
        if not (s_sb > 1e-6):
            pytest.fail(f"S_sb={s_sb:.3e} ≈ 0 — scattering does NOT source the ray; the bulk→ray "
                        f"coupling A_BA is missing (the ψ½ within-group source is unwired).")
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
        bulk=AngularFlux.from_mesh(np.zeros((N, ng, nx)), sn),
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
        bulk=AngularFlux.from_mesh(rng.standard_normal((N, ng, nx)), sn),
        boundary=AngularBoundaryFlux(
            values=rng.standard_normal(n_tr), space=sn.angular_trace, mesh=sn),
        radial_characteristic=RadialCharacteristicFlux(
            values=rng.standard_normal(ns), space=sn.radial_characteristic_space, mesh=sn),
    )


def _dense_seed(fn, sn) -> NDArray:
    """Densify a seed-block operator (``FullField -> FullField``) by probing the
    ψ½ basis — the ``(n_sd, n_sd)`` matrix of its radial_characteristic block."""
    ns = sn.radial_characteristic_space.shape[0]
    M = np.zeros((ns, ns))
    for j in range(ns):
        e = np.zeros(ns)
        e[j] = 1.0
        out = fn(_seed_composite(sn, e))
        M[:, j] = out.radial_characteristic.values
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
        r"""``B_b`` reflects the ray corner but emits a PRESENT-ZERO trace.
        Mutation tooth: a ``B_b`` leaking a trace action emits a non-zero trace."""
        sn = _sphere(bc="reflective")
        out = RadialCharacteristicBoundaryOperator(sn).apply(
            _random_composite(sn, np.random.default_rng(2)))
        np.testing.assert_array_equal(
            out.boundary.values, 0.0,
            err_msg="B_b emitted a NON-ZERO trace block — it leaked a trace action "
                    "that belongs to B_a.")
        if not np.max(np.abs(out.radial_characteristic.values)) > 0.0:
            pytest.fail("B_b emitted a zero ray corner on a reflective sphere — the "
                        "System B boundary arm is dead.")

    def test_sum_reconstructs_both_blocks_disjointly(self):
        r"""``(B_a + B_b).apply`` = ``B_a``'s trace ⊕ ``B_b``'s ray, byte-for-byte
        per block — the disjoint direct sum reconstructs the whole boundary."""
        sn = _sphere(bc="reflective")
        psi = _random_composite(sn, np.random.default_rng(3))
        B_a = SNBoundaryOperator(sn)
        B_b = RadialCharacteristicBoundaryOperator(sn)
        out = (B_a + B_b).apply(psi)
        out_a, out_b = B_a.apply(psi), B_b.apply(psi)
        # The composite trace is exactly B_a's (B_b contributes present-zero).
        np.testing.assert_array_equal(out.boundary.values, out_a.boundary.values)
        # The composite ray is exactly B_b's (B_a contributes present-zero).
        np.testing.assert_array_equal(
            out.radial_characteristic.values, out_b.radial_characteristic.values)

    def test_seedless_composite_sum_is_still_b_a(self):
        r"""On a seedless composite ``B_b`` is a no-op (ray ``None``); ``B_a + B_b``
        equals ``B_a`` on the trace, ray stays ``None`` (no mixed-presence raise).
        A slab (Cartesian) is genuinely seedless — a sphere carries the pole ray."""
        slab = SNMesh(
            Mesh1D(edges=np.linspace(0.0, 4.0, 6), mat_ids=np.zeros(5, dtype=int),
                   coord=CoordSystem.CARTESIAN, bc_right=BC("reflective"),
                   bc_left=BC("reflective")),
            Quadrature.gauss_legendre(4), {0: _mixture(1.0, 0.4, 2)})
        N, nx, ng = slab.quad.N, slab.nx, slab.ng
        n_tr = int(slab.angular_trace.layout.total_size)
        psi = FullField(
            bulk=AngularFlux.from_mesh(np.random.default_rng(4).standard_normal((N, ng, nx)), slab),
            boundary=AngularBoundaryFlux(
                values=np.random.default_rng(5).standard_normal(n_tr),
                space=slab.angular_trace, mesh=slab),
            radial_characteristic=None,
        )
        B_a = SNBoundaryOperator(slab)
        B_b = RadialCharacteristicBoundaryOperator(slab)
        out = (B_a + B_b).apply(psi)
        if out.radial_characteristic is not None:
            pytest.fail("seedless B_a + B_b emitted a non-None ray — B_b must pass "
                        "None through on a seedless composite.")
        np.testing.assert_array_equal(out.boundary.values, B_a.apply(psi).boundary.values)


class TestB_b_RayBoundary:
    r"""``B_b`` — the ψ½ ray-corner boundary law (RULINGS P1 + P2).

    The specular corner swap, its Euclidean transpose (the mirror), and — the
    load-bearing gate — the ``G_sd = V_cell`` reciprocity that keeps Mode-12
    closed (P2: the corner gauge is symmetric, so Euclidean = Hilbert). Reflective
    sphere for the non-trivial arm; vacuum for the null control.
    """

    def test_reflective_corner_swap_forward(self):
        r"""Forward: ``out.corner(level, −1) = seed.corner(level, +1)`` per level;
        the cells and the +1 corner stay zero (B_b touches only the inflow row)."""
        sn = _sphere(bc="reflective")
        space = sn.radial_characteristic_space
        seed_vals = np.random.default_rng(6).standard_normal(space.shape[0])
        out = RadialCharacteristicBoundaryOperator(sn).apply(_seed_composite(sn, seed_vals))
        ov = out.radial_characteristic.values
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
        fwd = _dense_seed(B_b.apply, sn)
        T = _dense_seed(B_b.apply_transpose, sn)
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
        fwd = _dense_seed(B_b.apply, sn)
        T = _dense_seed(B_b.apply_transpose, sn)
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
        out = RadialCharacteristicBoundaryOperator(sn).apply(_seed_composite(sn, seed_vals))
        np.testing.assert_array_equal(
            out.radial_characteristic.values, 0.0,
            err_msg="vacuum B_b emitted a non-zero corner (it did the reflective swap).")

    def test_unruled_outer_law_is_loud_deferred(self, monkeypatch):
        r"""``kind ∈ {white, albedo, periodic}`` → ``NotImplementedError`` with the
        specific message (NEGATIVE law, anti-#11: a bare ``raises`` false-greens on
        a downstream crash). Monkeypatch the xmax law kind (no white-sphere mesh
        needed) — auto-reverts, never a git checkout."""
        sn = _sphere(bc="reflective")
        monkeypatch.setattr(sn.bc["xmax"], "kind", "white")
        seed_vals = np.random.default_rng(12).standard_normal(
            sn.radial_characteristic_space.shape[0])
        with pytest.raises(NotImplementedError, match="no ruled corner action yet"):
            RadialCharacteristicBoundaryOperator(sn).apply(_seed_composite(sn, seed_vals))

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
    composite, L19)."""
    lhs = float(op.solve(u).values @ v.values)
    rhs = float(u.values @ op.solve_transpose(v).values)
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
            src_bar = op.solve_transpose(v)
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
        flux = op.solve(source)
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
        op.solve(_ray_source(sn, np.random.default_rng(5)))
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
        op.solve(source)
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
        a = op.solve(s0).cells(0, -1)
        b = op.solve(s1).cells(0, -1)
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
            op.solve(s0).cells(0, -1)[:, :-1] - op.solve(s1).cells(0, -1)[:, :-1])))
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
        flux = op.solve(src)
        expected = np.broadcast_to(C[:, None], (sn.ng, sn.nx))
        for sign in (-1, +1):
            np.testing.assert_allclose(
                flux.cells(0, sign), expected, atol=1e-13,
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
            q0 = _ray_source(sn, rng)
            psi0 = op.solve(q0)
            psi1 = op.solve(op.apply(psi0))
            np.testing.assert_allclose(
                psi1.values, psi0.values, rtol=1e-11, atol=1e-13)
            qr = op.apply(op.solve(q0))
            for p in space.levels:
                corner = space.corner_view(qr.values, p, +1)
                np.testing.assert_array_equal(corner, np.zeros_like(corner))

    def test_apply_routes_through_the_shared_forward_kernel(self, monkeypatch):
        # Mode-11 anti-twin (operator side): A_BB.apply MUST enter the shared
        # radial_characteristic_forward_residual, not an inlined copy.
        sn = _sphere()
        op = RadialCharacteristicOperator(sn, _ray_sigma(sn))
        calls = _install_forward_spy(monkeypatch)
        op.apply(_ray_cotangent(sn, np.random.default_rng(1)))
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
            lhs = float(op.apply(u).values @ v.values)
            rhs = float(u.values @ op.apply_transpose(v).values)
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
        q0 = _ray_source(sn, rng)
        psi0 = op.solve(q0)
        inv = op.inverse()
        np.testing.assert_array_equal(inv.apply(q0).values, op.solve(q0).values)
        np.testing.assert_array_equal(inv.solve(psi0).values, op.apply(psi0).values)
        if inv.inverse() is not op:
            pytest.fail("inverse().inverse() is not the original operator.")
        if not (op.is_invertible and op.is_adjointable):
            pytest.fail(
                "A_BB must report is_invertible and is_adjointable True at 4b.")


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
    phi0 = psi.bulk.integrate_angular().values          # (ng, nx)
    return np.asarray(S.isotropic_kernel.apply(phi0))


def _f_emission(F, psi: FullField) -> NDArray:
    """The ℓ=0 iso fission cell-emission ``χ·νΣf·φ`` that the F seed arm folds."""
    return np.asarray(F.apply(psi.bulk.integrate_angular()).values)


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
    r"""``A_BA`` — the ψ½ Schur fold (bulk → ray q½ source), the coupling Step 2
    un-welds from the five hand-rolled fold sites into ONE single source.

    Carrying member = **sphere-GL S4 ONLY** (the only geometry that carries a ψ½
    level, R12a; 1 level → 2 fold calls/arm). cylinder/slab are the non-carrying
    CONTROL. Every value row is ≥2G (1G is degenerate, vv anti-#3). Runtime: gates
    raise via ``pytest.fail`` / ``np.testing.assert_*`` (fire under ``python -O``),
    never a bare ``assert`` (vv Mode 8).

    Live TODAY (green, teeth mutation-verified in-process): the fold contract, the
    fold-transpose contract, the operator seed-arm transpose consistency, the S/F
    closed-form bit-identity, the shared-fold Mode-11 routing, the non-carrying
    control. xfail-skeleton (flip the ``# BIND`` in ``_apply_A_BA`` /
    ``_wrap_extracted_A_BA`` once the un-weld lands): the direct-A_BA-surface rows.
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

    # ── Gate 2a (LIVE): S and F both route through the ONE shared fold ──────

    def test_scattering_and_fission_both_route_through_the_shared_fold(self, monkeypatch):
        r"""Mode-11 CENTERPIECE (single-source routing). A wrap-COUNTER on the
        shared ``fold_moments_to_radial_characteristic`` (the SAME object the S/F
        seed arms local-import) must be entered by BOTH the S forward seed arm AND
        the F forward seed arm — a green-but-unrouted arm (a divergent inlined copy
        of the P_ℓ math) leaves the counter at 0 (Mode 11). Each arm folds
        ``2·n_levels`` times (2 signs/level; sphere-GL S4 → 1 level → 2 calls).

        POST-UN-WELD: re-point the wrap to the EXTRACTED single-source surface
        (``_wrap_extracted_A_BA`` — the ``from_moments`` factory / ``A_BA.apply``);
        the count then proves S and F route through the ONE extracted loop, not
        merely the shared inner math. (That sharper gate is the xfail row below.)
        """
        n_levels = 1  # sphere-GL S4 carries exactly one ψ½ level (levels == (0,))

        sn_s = _sphere()
        S = SNSolver(sn_s).scattering_op
        calls = _install_fold_spy(monkeypatch)
        S.apply(_random_composite(sn_s, np.random.default_rng(30)))
        s_calls = calls["n"]
        if s_calls != 2 * n_levels:
            pytest.fail(f"S folded {s_calls}× (expected 2·n_levels = {2 * n_levels}) — "
                        f"the S seed arm does not route through the shared fold.")

        sn_f = _fissile_sphere()
        F = SNSolver(sn_f).fission_op
        calls["n"] = 0
        F.apply(_random_composite(sn_f, np.random.default_rng(31)))
        f_calls = calls["n"]
        print(f"  [G2a] S fold-calls={s_calls}  F fold-calls={f_calls}  (both = 2·n_levels)")
        if f_calls != 2 * n_levels:
            pytest.fail(f"F folded {f_calls}× (expected 2·n_levels = {2 * n_levels}) — "
                        f"the F seed arm does not route through the shared fold "
                        f"(or the fissile emission is zero — check the mixture).")

    # ── Gate 3 (LIVE): the S / F seed output IS the ½·emission fold ─────────

    def test_scattering_seed_is_the_half_emission_fold(self):
        r"""The S forward seed output equals the ℓ=0 fold of its iso emission:
        byte-for-byte the documented old loop (``np.array_equal`` — the un-weld's
        bit-identity contract, inheriting from Gate 1) AND, structurally
        INDEPENDENT of the fold, the closed form ``½·q₀`` per (level, sign) with
        zero corners (P₀≡1; scattering is volumetric). Non-fissile ``Q/Σ``
        scattering config (refutation #4 — no ``nan`` k here)."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        space = sn.radial_characteristic_space
        psi = _random_composite(sn, np.random.default_rng(21))
        emission = _s_emission(S, psi)
        if not np.max(np.abs(emission)) > 1e-6:
            pytest.fail("scattering iso emission ≈ 0 — the S seed gate is vacuous.")
        seed = S.apply(psi).radial_characteristic
        if seed is None:
            pytest.fail("S emitted a None ray on a seed-carrying sphere composite.")
        # (a) bit-identity vs the documented hand-rolled loop (survives the un-weld).
        np.testing.assert_array_equal(
            seed.values, _ba_oldloop_reference(emission, sn),
            err_msg="S seed ≠ the pre-un-weld hand-rolled fold loop (the un-weld "
                    "diverged from the documented site-3 loop).")
        # (b) structurally-independent closed form: ½·emission cells, zero corners.
        for lv in space.levels:
            for sign in (-1, +1):
                np.testing.assert_allclose(
                    space.cells_view(seed.values, lv, sign), 0.5 * emission,
                    rtol=1e-13, atol=1e-14,
                    err_msg=f"level {lv} sign {sign}: S seed cells ≠ ½·q₀ (the ℓ=0 fold).")
                np.testing.assert_array_equal(
                    space.corner_view(seed.values, lv, sign), 0.0,
                    err_msg=f"level {lv} sign {sign}: S seed corner ≠ 0 (fold writes cells only).")

    def test_fission_seed_is_the_half_emission_fold_fissile(self):
        r"""The F forward seed output equals ``½·(χ·νΣf·φ)`` per (level, sign),
        zero corners — same bit-identity + closed-form pair as S. FISSILE mixture
        (refutation #4: a non-fissile mixture makes the emission — hence the seed —
        identically zero, a VACUOUS gate; the non-vacuity guard asserts a genuine
        nonzero emission)."""
        sn = _fissile_sphere()
        F = SNSolver(sn).fission_op
        space = sn.radial_characteristic_space
        psi = _random_composite(sn, np.random.default_rng(22))
        emission = _f_emission(F, psi)
        if not np.max(np.abs(emission)) > 1e-6:
            pytest.fail("fission iso emission ≈ 0 — the F seed gate is VACUOUS "
                        "(mixture is not actually fissile / νΣf = 0).")
        seed = F.apply(psi).radial_characteristic
        if seed is None:
            pytest.fail("F emitted a None ray on a seed-carrying fissile sphere composite.")
        np.testing.assert_array_equal(
            seed.values, _ba_oldloop_reference(emission, sn),
            err_msg="F seed ≠ the pre-un-weld hand-rolled fold loop (the F seed "
                    "arm before the Step-2 un-weld).")
        for lv in space.levels:
            for sign in (-1, +1):
                np.testing.assert_allclose(
                    space.cells_view(seed.values, lv, sign), 0.5 * emission,
                    rtol=1e-13, atol=1e-14,
                    err_msg=f"level {lv} sign {sign}: F seed cells ≠ ½·fission_iso.")
                np.testing.assert_array_equal(
                    space.corner_view(seed.values, lv, sign), 0.0,
                    err_msg=f"level {lv} sign {sign}: F seed corner ≠ 0.")

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

    # ── Gate 5b (LIVE): the S seed-arm forward/adjoint transpose consistency ─

    def test_scattering_seed_arm_euclidean_transpose_consistency(self, monkeypatch):
        r"""The S adjoint seed arm (site 5, hard-coded ``0.5``) is the EXACT
        Euclidean transpose of the S forward seed arm:
        ``⟨A_BA·φ, χ̄⟩ = ⟨φ, A_BAᵀ·χ̄⟩`` (< 1e-11), with the adjoint's OUTPUT ray
        block present-zero (``∂S/∂ψ½ = 0``). This is the gate that pins the
        hand-rolled ``0.5`` is CONSISTENT with the forward fold — the risk the
        un-weld's single-sourcing eliminates by construction.

        Tooth: monkeypatch the FORWARD fold to ``0.6`` at ℓ=0 (the forward
        local-imports it; the adjoint hand-rolls ``0.5``, unaffected) → the
        forward/adjoint coefficients DISAGREE → the identity reds (measured ≈ 0.09).
        A shared-coefficient value is invisible to this consistency gate (both
        legs scale together) — the fold's VALUE is pinned by Gate 1/3, this gate
        pins the fwd↔adj CONSISTENCY. The tooth survives the un-weld (the fold and
        the fold-transpose stay distinct entry points)."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        space = sn.radial_characteristic_space
        N, nx, ng = sn.quad.N, sn.nx, sn.ng
        n_tr = int(sn.angular_trace.layout.total_size)
        n_sd = space.shape[0]
        rng = np.random.default_rng(40)

        def euclid_defect() -> float:
            phi_ff = FullField(
                bulk=AngularFlux.from_mesh(rng.standard_normal((N, ng, nx)), sn),
                boundary=AngularBoundaryFlux(values=np.zeros(n_tr), space=sn.angular_trace, mesh=sn),
                radial_characteristic=RadialCharacteristicFlux(
                    values=np.zeros(n_sd), space=space, mesh=sn))
            chi = rng.standard_normal(n_sd)
            chi_ff = FullField(
                bulk=AngularFlux.from_mesh(np.zeros((N, ng, nx)), sn),
                boundary=AngularBoundaryFlux(values=np.zeros(n_tr), space=sn.angular_trace, mesh=sn),
                radial_characteristic=RadialCharacteristicFlux(values=chi, space=space, mesh=sn))
            seed_out = S.apply(phi_ff).radial_characteristic.values         # A_BA·φ
            adj = S.apply_transpose(chi_ff)
            # ∂S/∂ψ½ = 0: the adjoint's ray block is present-zero.
            if adj.radial_characteristic is not None:
                np.testing.assert_array_equal(
                    adj.radial_characteristic.values, 0.0,
                    err_msg="S adjoint ray block ≠ present-zero (∂S/∂ψ½ must be 0).")
            bulk_out = adj.bulk.values                                       # A_BAᵀ·χ̄
            lhs = float(seed_out @ chi)
            rhs = float(phi_ff.bulk.values.ravel() @ bulk_out.ravel())
            return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)

        control = euclid_defect()
        if not control < 1e-11:
            pytest.fail(f"seed-arm Euclidean transpose defect {control:.3e} ≥ 1e-11 — "
                        f"the S adjoint seed arm's 0.5 is NOT the transpose of the "
                        f"forward fold (forward/adjoint coefficients disagree).")

        # Tooth: force a forward/adjoint coefficient disagreement (forward → 0.6).
        def _fold_06(moments, sign):
            moments = np.asarray(moments)
            ell = np.arange(moments.shape[0])
            coeff = ((2.0 * ell + 1.0) / 2.0) * np.float64(sign) ** ell
            coeff = coeff.copy()
            coeff[0] = 0.6
            return np.tensordot(coeff, moments, axes=(0, 0))

        # Patch the fold as the RECONSTRUCTION sees it (forward → 0.6); the
        # adjoint uses the SEPARATE fold-transpose (unpatched, 0.5) → the fwd/adj
        # coefficients disagree → the identity reds.
        monkeypatch.setattr(_rcr_mod, "fold_moments_to_radial_characteristic", _fold_06)
        tooth = euclid_defect()
        print(f"  [G5b] control={control:.2e}  fwd-fold=0.6 tooth={tooth:.3f}")
        if not tooth > 1e-3:
            pytest.fail(f"the forward-fold=0.6 mismatch left the transpose defect at "
                        f"{tooth:.3e} — the fwd↔adj consistency gate has no teeth.")

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
                bulk=AngularFlux.from_mesh(
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

    # ── xfail-skeleton: the direct extracted-A_BA surface (BIND undecided) ──

    def test_A_BA_direct_surface_folds_half_emission(self):
        r"""The extracted single-source A_BA surface
        (``RadialCharacteristicReconstruction.apply``) folds an emission to the
        ray q½ source, matching the closed-form ½·emission loop
        (``_ba_oldloop_reference``). LIVE post-un-weld (Step 2)."""
        sn = _sphere()
        S = SNSolver(sn).scattering_op
        psi = _random_composite(sn, np.random.default_rng(60))
        emission = _s_emission(S, psi)
        got = _apply_A_BA(emission, sn)
        np.testing.assert_array_equal(
            got, _ba_oldloop_reference(emission, sn),
            err_msg="the extracted A_BA surface ≠ the documented fold loop.")

    def test_scattering_and_fission_route_through_extracted_A_BA(self, monkeypatch):
        r"""Mode-11, SHARPENED: S and F both route through the EXTRACTED
        single-source surface (``RadialCharacteristicReconstruction.apply``), not
        merely the shared inner fold. A green-but-unrouted arm (re-inlining the
        loop) leaves the counter at 0 and reds. LIVE post-un-weld (Step 2)."""
        sn_s = _sphere()
        S = SNSolver(sn_s).scattering_op
        calls = _wrap_extracted_A_BA(monkeypatch)
        S.apply(_random_composite(sn_s, np.random.default_rng(70)))
        if calls["n"] <= 0:
            pytest.fail("S did not route through the extracted A_BA surface (Mode 11).")
        sn_f = _fissile_sphere()
        F = SNSolver(sn_f).fission_op
        calls["n"] = 0
        F.apply(_random_composite(sn_f, np.random.default_rng(71)))
        if calls["n"] <= 0:
            pytest.fail("F did not route through the extracted A_BA surface (Mode 11).")

    def test_A_BA_transpose_surface_euclidean_contract(self):
        r"""The operator-shape A_BA exposes ``.apply_transpose``, which satisfies
        the Euclidean adjoint contract against the forward fold
        (``⟨A_BA·emission, y⟩ = ⟨emission, A_BAᵀ·y⟩``). LIVE post-un-weld (Step 2 —
        the operator shape carries a real transpose surface)."""
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
# (L+C).apply is isolated by LINEARITY (bulk=0, boundary=0 → only the seed's
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
        bulk=AngularFlux.from_mesh(bulk_values, sn),
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
        (bulk=0, boundary=0 isolate ``A_bs`` by linearity), BIT-IDENTICALLY
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
        # contribution to (L+C).apply — bulk=0/boundary=0 ⇒ A_AA·0 = 0, so the
        # bulk output is exactly A_AB·s.
        reference = _loss(sn).apply(_seed_composite(sn, sv)).bulk
        ref_other_sigma = _loss(sn, slope=0.9).apply(_seed_composite(sn, sv)).bulk
        np.testing.assert_array_equal(
            reference.values, ref_other_sigma.values,
            err_msg="the seed→bulk contribution changed with σ_t — A_AB is not "
                    "σ-independent (the isolation-by-linearity premise is wrong).")
        # Mode-11 sentinel on the shared closure kernel.
        pre = _install_closure_spy(monkeypatch, sn, "precompute_psi_state")
        cc = _install_closure_spy(monkeypatch, sn, "cell_contribution")
        out = RadialCharacteristicSeeding(sn).apply(seed)
        if len(pre) != 1:
            pytest.fail(
                f"precompute_psi_state called {len(pre)}× (expected 1) — A_AB."
                f"apply is not the single-precompute WRAP.")
        psi_view_arg = np.asarray(pre[0]["args"][0])
        if np.max(np.abs(psi_view_arg)) != 0.0:
            pytest.fail(
                "A_AB.apply did NOT zero the bulk psi_view — A_AA's angular "
                "redistribution would leak into the isolated coupling.")
        if pre[0]["kwargs"].get("radial_characteristic") is not seed:
            pytest.fail(
                "A_AB.apply did not pass the seed as radial_characteristic.")
        if len(cc) < sn.nx:
            pytest.fail(
                f"cell_contribution called {len(cc)}× (< nx = {sn.nx}) — the "
                f"cell-local angular injection did not visit every cell.")
        np.testing.assert_array_equal(
            out.values, reference.values,
            err_msg="A_AB.apply is not bit-identical to the in-sweep injection.")

    # ── Transpose — the seed_cells_bar term ≡ the in-sweep reverse ─────────

    def test_apply_transpose_is_the_in_sweep_seed_cells_bar(self, monkeypatch):
        r"""``A_AB.apply_transpose(v).cells(p,-1)`` ≡ the ``seed_cells_bar`` term
        the in-sweep reverse adds — ``(L+C).apply_transpose(bulk=v, ray=0)`` on
        the ray block (ray=0 nulls the A_BB self-block, isolating the M-M thread
        cotangent), BIT-IDENTICALLY. The ``+1`` leg and both corners stay EXACTLY
        0 (the forward writes only the inward leg). A Mode-11 sentinel proves
        ``angular_adjoint`` runs exactly once. ≥2G."""
        sn = _sphere()
        space = sn.radial_characteristic_space
        rng = np.random.default_rng(31)
        vv = rng.standard_normal((sn.quad.N, sn.ng, sn.nx))
        v = AngularSourceSink.from_mesh(vv, sn)
        reference = _loss(sn).apply_transpose(
            _bulk_composite(sn, vv)).radial_characteristic
        aa = _install_closure_spy(monkeypatch, sn, "angular_adjoint")
        out = RadialCharacteristicSeeding(sn).apply_transpose(v)
        if len(aa) != 1:
            pytest.fail(
                f"angular_adjoint called {len(aa)}× (expected 1) — A_AB."
                f"apply_transpose is not the single-adjoint WRAP.")
        for p in space.levels:
            np.testing.assert_array_equal(
                out.cells(p, -1), reference.cells(p, -1),
                err_msg=f"level {p}: apply_transpose cells(-1) ≠ the in-sweep "
                        f"seed_cells_bar.")
            np.testing.assert_array_equal(
                space.cells_view(out.values, p, +1), 0.0,
                err_msg=f"level {p}: apply_transpose wrote the +1 leg (must be 0).")
            np.testing.assert_array_equal(
                space.corner_view(out.values, p, -1), 0.0,
                err_msg=f"level {p}: apply_transpose wrote the -1 corner (be 0).")
            np.testing.assert_array_equal(
                space.corner_view(out.values, p, +1), 0.0,
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
            lhs = float(op.apply(u).values.ravel() @ v.values.ravel())
            rhs = float(u.values.ravel() @ op.apply_transpose(v).values.ravel())
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
        lhs = float(op.apply(u).values.ravel() @ v.values.ravel())
        rhs = float(u.values.ravel() @ op.apply_transpose(v).values.ravel())
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
        out_full = op.apply(
            RadialCharacteristicFlux(values=full, space=space, mesh=sn))
        out_minus = op.apply(
            RadialCharacteristicFlux(values=only_minus, space=space, mesh=sn))
        if not np.max(np.abs(out_full.values)) > 1e-6:
            pytest.fail("A_AB.apply output is ~0 — the asymmetry gate is vacuous.")
        np.testing.assert_array_equal(
            out_full.values, out_minus.values,
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
            op.apply(_seed_flux(other, np.random.default_rng(9)))
        with pytest.raises(ValueError, match="mesh-identity invariant"):
            op.apply_transpose(_bulk_cotangent(other, np.random.default_rng(9)))


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

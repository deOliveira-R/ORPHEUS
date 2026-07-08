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
from orpheus.sn.solver import SNSolver, _within_group_triple
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.full_field import FullField
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

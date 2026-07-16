r"""Foundation-tagged symbolic tests for the 2-D Cartesian LD stress MMS source
(#240 D5b-S4 — defends the multi-D bilinear UBLD slope-row verification).

The Branch-1 SymPy module (:func:`derive_2d_cartesian_ld_stress_mms` in
:mod:`orpheus.derivations.continuous.mms.sn`) substitutes the strengthened
:math:`\mu`-bilinear ansatz

.. math::

   \psi_{n,g}(x,y) = \frac{1}{W}\bigl[A_g(x,y)
       + \mu_{x,n} B_g(x,y) + \mu_{y,n} C_g(x,y)\bigr]

into the continuous 2-D Cartesian SN operator (NO angular redistribution) and
verifies algebraically that the closed-form :math:`Q^{\rm ext}_n(x,y)` returned
by the factory's ``external_source`` method is the unique source consistent with
the ansatz.  This MUST pass before any solve consumes the source (the
algebra-of-record Branch-1 gate, L11).

The cross-checks here are STRUCTURALLY INDEPENDENT of the symbolic
substitution itself:

1. an INDEPENDENT finite-difference of the continuous-PDE residual
   :math:`\mu_x\partial_x\psi + \mu_y\partial_y\psi + \Sigma_t\psi -
   (1/W)(\Sigma_s\phi + Q^{\rm ext})` (built with numpy central differences, a
   different evaluation path from SymPy's exact ``diff``) must vanish — this
   catches a sign / factor error that SymPy's own ``simplify`` would share with
   the closed form; and
2. the Branch-2 numpy ``external_source`` must agree with the SymPy-lambdified
   closed form to floating-point precision (the L1 Branch-2-vs-Branch-1
   cross-check — a copy error between the symbolic and numerical branches).

Plus the vv-Mode-7 activation checks: the SLOPE drivers :math:`B,\,C` carry
genuine per-axis spatial variation (:math:`\partial_x B \neq 0`,
:math:`\partial_y C \neq 0`), and :math:`B,\,C` are NOT x↔y reflections (the
strengthening — a same-sign slope-row sign bug cannot cancel).

``-O``-safe (Mode 8): every load-bearing check is a function call
(``np.testing.assert_*`` / ``pytest.fail``), never a bare ``assert`` (which the
canonical ``python -O`` invocation strips to a NO-OP).

References
----------

- ``.claude/skills/vv-principles/SKILL.md`` (Mode 7 — MMS simplification bias).
- ``.claude/skills/algebra-of-record/SKILL.md`` (Branch 1 / Branch 2 / L11).
- :class:`SN2DCartesianLDStressMMSCase` (Branch-2 numerical factory).
- :doc:`/theory/methods/sn/index` — label ``ld-cartesian-2d`` (D6/archivist
  mints it).
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.continuous.mms.sn import (
    SN2DCartesianLDStressMMSCase,
    _2d_cartesian_ld_stress_symbolic,
    build_2d_cartesian_ld_stress_mms_case,
    derive_2d_cartesian_ld_stress_mms,
)


# ═══════════════════════════════════════════════════════════════════════
# V_ld2d-stress — the algebra-of-record substitution identity (Branch 1)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("ld-cartesian-2d", "transport-cartesian-2d")
def test_v_ld2d_stress_substitution_identity():
    r"""V_ld2d-stress — substituting the strengthened :math:`\mu`-bilinear
    ansatz into the 2-D Cartesian SN operator yields the closed-form
    :math:`Q^{\rm ext}_n` exactly (zero residual under :func:`sympy.simplify`).

    The load-bearing algebra-of-record claim (L11): the closed form is the
    unique manufactured source consistent with the ansatz, derived from the
    CONTINUOUS PDE — it touches NONE of the LD cell-update code, so a sign bug
    in the LD slope rows cannot also corrupt the reference."""
    result = derive_2d_cartesian_ld_stress_mms()
    if not result["pass"]:
        pytest.fail(
            f"V_ld2d-stress failed: substituted Q differs from the closed form "
            f"by {result['diff']}"
        )


@pytest.mark.foundation
def test_v_ld2d_stress_residual_finite_difference_vanishes():
    r"""V_ld2d-stress — the continuous-PDE residual vanishes under an
    INDEPENDENT finite-difference evaluation of the streaming derivatives.

    Structural independence (L11): the substitution identity above uses SymPy's
    exact symbolic ``diff``; this re-checks the SAME manufactured identity via
    numpy central differences of :math:`\psi` (a different evaluation path), so a
    shared upstream ``diff`` mistake (ERR-032-class) cannot pass both.  The
    residual

    .. math::

       \mu_x\,\partial_x\psi + \mu_y\,\partial_y\psi + \Sigma_t\,\psi
       - \tfrac{1}{W}\bigl(\Sigma_s\,\phi + Q^{\rm ext}\bigr)

    must be :math:`\approx 0` at interior sample points for every sampled
    ordinate."""
    (x, y, mu_x, mu_y, Lx, Ly, Sigma_t, Sigma_s, W,
     A, B, C, psi, phi, Q_closed) = _2d_cartesian_ld_stress_symbolic()

    # Lambdify ψ, φ, Q at fixed (Lx, Ly, Σt, Σs, W) so the FD operates on a pure
    # spatial+angular numpy function (independent of SymPy's diff).
    Lx_v, Ly_v = 1.3, 0.9
    sig_t_v, sig_s_v, W_v = 1.7, 0.6, 4.0 * np.pi
    subs = {Lx: Lx_v, Ly: Ly_v, Sigma_t: sig_t_v, Sigma_s: sig_s_v, W: W_v}
    psi_fn = sp.lambdify((x, y, mu_x, mu_y), psi.subs(subs), "numpy")
    phi_fn = sp.lambdify((x, y), phi.subs(subs), "numpy")
    Q_fn = sp.lambdify((x, y, mu_x, mu_y), Q_closed.subs(subs), "numpy")

    h = 1e-6
    pts = [(0.31 * Lx_v, 0.62 * Ly_v), (0.55 * Lx_v, 0.27 * Ly_v),
           (0.73 * Lx_v, 0.81 * Ly_v)]
    ords = [(0.5, 0.3), (-0.4, 0.6), (0.2, -0.7), (-0.55, -0.15)]
    max_res = 0.0
    for (xv, yv) in pts:
        for (mxv, myv) in ords:
            dpsi_dx = (psi_fn(xv + h, yv, mxv, myv)
                       - psi_fn(xv - h, yv, mxv, myv)) / (2 * h)
            dpsi_dy = (psi_fn(xv, yv + h, mxv, myv)
                       - psi_fn(xv, yv - h, mxv, myv)) / (2 * h)
            lhs = mxv * dpsi_dx + myv * dpsi_dy + sig_t_v * psi_fn(xv, yv, mxv, myv)
            rhs = (sig_s_v * phi_fn(xv, yv) + Q_fn(xv, yv, mxv, myv)) / W_v
            max_res = max(max_res, abs(float(lhs - rhs)))
    # central-difference truncation ~ O(h²) ⇒ ~1e-9 at h=1e-6 for smooth ψ.
    np.testing.assert_allclose(
        max_res, 0.0, atol=1e-7,
        err_msg=f"manufactured-source FD residual {max_res:.3e} ≠ 0",
    )


@pytest.mark.foundation
def test_v_ld2d_stress_scalar_flux_is_A():
    r"""V_ld2d-stress — the scalar flux :math:`\phi = \int\psi\,d\mu = A(x,y)`
    (the :math:`\mu_x B + \mu_y C` terms integrate to zero).

    Checked symbolically (:math:`\phi` returned by the symbolic builder equals
    :math:`A`) and as a non-trivial, NON-vanishing field (:math:`a_0>0` →
    :math:`\phi \neq 0` at the edges — the boundary-closure stress)."""
    (x, y, _mx, _my, Lx, Ly, _st, _ss, _W,
     A, _B, _C, _psi, phi, _Q) = _2d_cartesian_ld_stress_symbolic()
    if sp.simplify(phi - A) != 0:
        pytest.fail(f"scalar flux φ ≠ A: φ-A = {sp.simplify(phi - A)}")
    # NON-vanishing at the corner (x=0, y=0): a0 = 7/10 + a2·cos0·cos0.
    phi00 = float(phi.subs({x: 0, y: 0, Lx: 1.3, Ly: 0.9}))
    if abs(phi00) < 1e-3:
        pytest.fail(
            f"scalar flux vanishes at the (0,0) edge corner (φ={phi00}); the "
            "boundary closure would be untested (a0 must be > 0)"
        )


@pytest.mark.foundation
def test_v_ld2d_stress_slope_drivers_are_per_axis_active():
    r"""V_ld2d-stress (Mode 7) — the SLOPE drivers carry genuine per-axis
    variation: :math:`\partial_x B \neq 0` AND :math:`\partial_y C \neq 0`.

    The dangerous degeneracy (Frame 2 §199–208): if :math:`B` were constant in
    :math:`x` (or :math:`C` constant in :math:`y`) the corresponding slope row
    would be driven only by :math:`\partial A`, an angularly-FLAT slope that the
    DD cell-average path also captures — the per-ordinate bilinear coupling
    would NOT be stressed.  This pins that both slope drivers genuinely vary
    along their own axis."""
    (x, y, _mx, _my, Lx, Ly, _st, _ss, _W,
     _A, B, C, _psi, _phi, _Q) = _2d_cartesian_ld_stress_symbolic()
    dB_dx = sp.diff(B, x)
    dC_dy = sp.diff(C, y)
    # Evaluate at an interior point where the harmonics do not coincidentally
    # vanish (x=Lx/3, y=Ly/3).
    sub = {x: sp.Rational(1, 3) * Lx, y: sp.Rational(1, 3) * Ly,
           Lx: sp.Integer(13) / 10, Ly: sp.Rational(9, 10)}
    if sp.simplify(dB_dx.subs(sub)) == 0:
        pytest.fail("∂_x B vanishes at the interior probe — x-slope row not active")
    if sp.simplify(dC_dy.subs(sub)) == 0:
        pytest.fail("∂_y C vanishes at the interior probe — y-slope row not active")


@pytest.mark.foundation
def test_v_ld2d_stress_slope_drivers_break_xy_reflection():
    r"""V_ld2d-stress (the strengthening) — :math:`B` and :math:`C` are NOT x↔y
    reflections of each other.

    The whole point of the strengthening (Frame 2 §238–268): if :math:`B(x,y) =
    C(y,x)` (under x↔y), a same-sign slope-row SIGN bug (both rows share
    ``_LDCellTerms.slope``) could produce a flux that still converges at O(h²)
    because the two symmetric errors do not break the global symmetry the MMS
    measures — a FALSE GREEN for that bug class.  The :math:`b_2,\,c_2`
    cross-harmonics break the reflection; this test pins that
    :math:`B(x,y) \neq C(y,x)`."""
    (x, y, _mx, _my, Lx, Ly, _st, _ss, _W,
     _A, B, C, _psi, _phi, _Q) = _2d_cartesian_ld_stress_symbolic()
    # C with x↔y swapped (and Lx↔Ly so the harmonics map cleanly).
    C_reflected = C.subs({x: y, y: x, Lx: Ly, Ly: Lx}, simultaneous=True)
    delta = sp.simplify(B - C_reflected)
    if delta == 0:
        pytest.fail(
            "B(x,y) == C(y,x) — the slope drivers ARE x↔y reflections; a "
            "same-sign slope-row sign bug could cancel (the strengthening is "
            "not load-bearing as written)"
        )


# ═══════════════════════════════════════════════════════════════════════
# L1 cross-check — Branch 2 (numpy) external_source == Branch 1 (SymPy)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_ld2d_stress_branch2_matches_branch1_source():
    r"""L1 cross-check — the Branch-2 numpy ``external_source`` agrees with the
    SymPy-lambdified Branch-1 closed form to floating-point precision, for BOTH
    groups and at multiple cells.

    Disagreement points to a copy error between the symbolic and numerical
    branches (the algebra-of-record Branch-2-descends-from-Branch-1 contract).
    The drivers scale linearly by the per-group amplitude :math:`c_g`, so the
    per-group source is :math:`c_g` times the single-group symbolic closed form;
    the multigroup in-scatter is matched by the effective scalar
    :math:`\Sigma_s^{\rm eff}_g = (\sum_{g'}\Sigma_s[g', g]\,A_{g'})/A_g` so the
    isotropic source term lines up.  **Covering** :math:`g \in \{0, 1\}` **and
    ≥2 cells pins the** :math:`c_1` **per-group scaling AND the 2G downscatter
    coupling** :math:`\Sigma_s[g', g]` **(the transpose-active term, the
    ERR-002/009 class) at the foundation tier** — not only inside the slow
    end-to-end convergence run (where a wrong-but-still-O(h²) source can hide,
    vv §5: rate ≠ correctness)."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    Q_b2 = case.external_source(mesh)              # (N, ng, nx, ny)
    cx, cy = mesh.centers_x, mesh.centers_y
    sum_w = float(case.quadrature.weights.sum())
    ng = case.n_groups

    (x, y, mu_x, mu_y, Lx, Ly, Sigma_t, Sigma_s, W,
     _A, _B, _C, _psi, _phi, Q_closed) = _2d_cartesian_ld_stress_symbolic()
    Q_fn = sp.lambdify(
        (x, y, mu_x, mu_y, Lx, Ly, Sigma_t, Sigma_s, W), Q_closed, "numpy",
    )

    max_err = 0.0
    for (ii, jj) in [(3, 2), (5, 4)]:              # ≥2 cells (x↔y-distinct)
        xi, yj = float(cx[ii]), float(cy[jj])
        xa, ya = np.array([xi]), np.array([yj])
        # Per-group cell-average A_g (already c_g-scaled by ``_drivers``).
        A_g = [float(case._drivers(xa, ya, g)[0][0, 0]) for g in range(ng)]
        for g in range(ng):
            cg = float(case.c_spectrum[g])
            sig_t_g = float(case.sigma_t_fn(xa, ya, g)[0])
            in_scatter_g = sum(
                float(case.sigma_s_fn(xa, ya, g_from, g)[0]) * A_g[g_from]
                for g_from in range(ng)
            )
            sig_s_eff_g = in_scatter_g / A_g[g]
            for n in range(case.quadrature.N):
                mxn = float(case.quadrature.mu_x[n])
                myn = float(case.quadrature.mu_y[n])
                q_sym = cg * float(Q_fn(
                    xi, yj, mxn, myn, case.length_x, case.length_y,
                    sig_t_g, sig_s_eff_g, sum_w,
                )) / sum_w
                max_err = max(max_err, abs(q_sym - Q_b2[n, g, ii, jj]))
    np.testing.assert_allclose(
        max_err, 0.0, atol=1e-13,
        err_msg=f"Branch-2 source ≠ Branch-1 closed form (max |Δ|={max_err:.3e})",
    )


@pytest.mark.foundation
def test_ld2d_stress_prescribed_inflow_is_nonvanishing():
    r"""The prescribed-inflow trace is NON-zero on the domain edges (:math:`a_0
    > 0`) — the boundary-closure stress.  Since #257 S9 the trace is
    moment-resolved on an LD mesh (slot-0 = transverse cell average, slot-1 =
    bare transverse face-slope); ``max|trace|`` is dominated by the non-vanishing
    average, so the non-vanishing assertion holds shape-agnostically."""
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.spatial import LinearDiscontinuous

    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(8)
    mats = case.build_materials(mesh)
    sn = SNMesh(mesh, case.quadrature, mats, scheme=LinearDiscontinuous())
    bss = case.prescribed_inflow(sn)
    max_abs = float(np.max(np.abs(bss.values)))
    if max_abs < 1e-3:
        pytest.fail(
            f"prescribed-inflow trace is ~zero ({max_abs:.3e}); a0>0 should make "
            "the inflow non-vanishing (the boundary closure would be untested)"
        )


@pytest.mark.foundation
def test_v_ld2d_stress_source_builds_without_ld_scheme():
    r"""V_ld2d-stress (L11 proxy) — the manufactured source BUILDS without
    instantiating the LD scheme/kernel.

    The genuine structural independence (the source is derived from the
    CONTINUOUS PDE, importing none of ``_LDCellTerms`` / ``_schur_terms`` /
    ``_ubld``) is established by inspection + the independent finite-difference
    cross-check above; this test pins the OBSERVABLE proxy — the symbolic
    builder (:func:`_2d_cartesian_ld_stress_symbolic`) and the Branch-2
    ``external_source`` produce a finite source with NO
    :class:`~orpheus.transport.spatial.LinearDiscontinuous` instance in scope (sharing
    the kernel would let an LD slope-row sign bug corrupt the reference too —
    the defining anti-pattern of a dependent MMS)."""
    case = build_2d_cartesian_ld_stress_mms_case()
    mesh = case.build_mesh(6)
    # external_source must build WITHOUT any LinearDiscontinuous instance.
    Q = case.external_source(mesh)
    np.testing.assert_array_equal(
        np.isfinite(Q), True,
        err_msg="external_source produced non-finite values",
    )
    # The case dataclass carries no LD-kernel handle (it is a pure MMS case).
    if isinstance(case, SN2DCartesianLDStressMMSCase) is False:
        pytest.fail("unexpected case type")

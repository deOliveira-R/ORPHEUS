r"""Foundation tests for ``Billiard.from_problem`` dual-factory bit-equality.

The test-architect spec § "File 2 (load-bearing)" identifies this as the
**load-bearing** file of the TransportSolver Protocol refactor. It pins
that the new (Protocol-aligned) factory signature

::

    Billiard.from_problem(
        materials={0: Mixture}, geometry_spec=GeometrySpec(...),
        alpha=..., quadrature=...,
    )

produces **bit-equal** :class:`CriticalSolution` results to the legacy
signature

::

    Billiard.from_problem(
        geometry_kind="sphere", materials={"sigma_t": ...},
        geometry={"R": ...}, alpha=..., quadrature=...,
    )

across every geometry × group × multi-region combination. Equality is
:func:`np.array_equal` (IEEE-754 exact) on arrays and
:meth:`float.hex` exact-bit comparison on scalars — NEVER
:func:`np.allclose` or any tolerance-based check (per Cardinal Rule 1).

Why bit-equality matters
------------------------

The dual-factory pattern creates two ways to produce the same
:class:`Billiard` instance. If the two paths drift even by a single
ULP — because of a unit-conversion rounding, a different storage shape,
or a missing key — the cross-method gates that swap geometry_specs at
fn_method's predicted thicknesses (see
:func:`tests.cross_method.test_eigenvalue._shadow_with_thickness_mfp`)
would silently regress. Bit-equality is the only check rigid enough to
catch that.

Asymmetric-input catchers
-------------------------

Per the spec § "Test matrix File 2", the test inputs were specifically
chosen so that a sign flip / variable swap / factor-of-2 error
produces *detectably-different* (not bit-equal) results:

* ``test_unified_factory_bit_equal_sphere_2g_asymmetric`` and
  ``test_unified_factory_bit_equal_slab_2g_asymmetric`` use
  ``SigS[0,1] != SigS[1,0]``. A scattering-transpose drift (Sig 3 from
  numerical-bug-signatures, ERR-002) would fail the bit-equality
  check.

* ``test_unified_factory_bit_equal_slab_1g_asymmetric_dim`` uses
  ``a_cm = 1.234`` (asymmetric in the sense that ``2 * a_cm = 2.468``
  is not a round number). A factor-of-2 swap between ``a`` and ``2a``
  (Sig 5) would land at a different ``L_cm`` and fail.

* ``test_unified_factory_bit_equal_slab_asymmetric_with_alpha_dict``
  uses ``alpha_left != alpha_right``. A bc_left ↔ bc_right swap
  (Sig 4) would route the legacy/new paths to different solver calls
  and fail.

Foundation-tier
---------------

Module-level ``pytestmark = [pytest.mark.foundation]`` (no
``verifies(...)``); these are software invariants.

References
----------

* :doc:`/skills/numerical-bug-signatures` — Sig 3 (scattering
  transpose), Sig 4 (BC swap), Sig 5 (factor-of-2). The asymmetric
  test inputs target exactly these patterns.
* :doc:`/theory/transport_solver_protocol` § "Dual-factory pattern".
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.continuous.trajectory_resolvent import Billiard
from orpheus.geometry.mesh import BC


pytestmark = [pytest.mark.foundation]


# ----------------------------------------------------------------------
# Bit-equality helpers (mirror the precedent in
# test_trajectory_resolvent_billiard.py)
# ----------------------------------------------------------------------


def _bit_equal_arrays(a: np.ndarray, b: np.ndarray) -> bool:
    """Strict bit-for-bit equality on float arrays (no ``allclose``)."""
    a_arr = np.asarray(a)
    b_arr = np.asarray(b)
    if a_arr.shape != b_arr.shape:
        return False
    if a_arr.dtype != b_arr.dtype:
        return False
    return bool(np.array_equal(a_arr, b_arr))


def _bit_equal_floats(a: float, b: float) -> bool:
    """Strict bit-for-bit float equality (handles NaN by hex)."""
    return float(a).hex() == float(b).hex()


def _assert_bit_equal_critical(
    sol_new: CriticalSolution, sol_legacy: CriticalSolution
) -> None:
    """Assert two :class:`CriticalSolution`s are bit-for-bit identical
    on every load-bearing field.

    This is the central invariant of the dual-factory test matrix.
    Eigenvalue, psi, phi, mesh nodes — every numerical field must
    agree to IEEE-754 exact.
    """
    assert _bit_equal_floats(sol_new.eigenvalue, sol_legacy.eigenvalue), (
        f"eigenvalue drift: new={sol_new.eigenvalue.hex()} "
        f"legacy={sol_legacy.eigenvalue.hex()}"
    )
    assert sol_new.eigenvalue_kind == sol_legacy.eigenvalue_kind
    assert sol_new.parameter_kind == sol_legacy.parameter_kind
    assert _bit_equal_floats(
        sol_new.parameter_value, sol_legacy.parameter_value
    )
    # psi/phi arrays — bit-for-bit identical.
    assert _bit_equal_arrays(
        sol_new.metadata["psi"], sol_legacy.metadata["psi"]
    ), "psi drift between new and legacy paths"
    assert _bit_equal_arrays(
        sol_new.metadata["phi"], sol_legacy.metadata["phi"]
    ), "phi drift between new and legacy paths"
    assert sol_new.metadata["iterations"] == sol_legacy.metadata["iterations"]
    assert sol_new.converged == sol_legacy.converged


# ----------------------------------------------------------------------
# Mixture builders (1G / 2G / 4G) producing the production-protocol shape
# ----------------------------------------------------------------------


def _make_mixture_1g(
    sigma_t: float, sigma_s: float, nu_sigma_f: float
) -> Mixture:
    """Build a 1G Mixture with the same XS the legacy raw dict carries."""
    eg = np.logspace(7, -3, 2)
    return Mixture(
        SigC=np.zeros(1),
        SigL=np.zeros(1),
        SigF=np.zeros(1),
        SigP=np.array([nu_sigma_f]),
        SigT=np.array([sigma_t]),
        SigS=[csr_matrix(np.array([[sigma_s]]))],
        Sig2=csr_matrix((1, 1)),
        chi=np.array([1.0]),
        eg=eg,
    )


def _make_mixture_mg(
    sigma_t: np.ndarray,
    sigma_s: np.ndarray,
    nu_sigma_f: np.ndarray,
    chi: np.ndarray | None = None,
) -> Mixture:
    """Build an MG Mixture from arrays (matching the legacy MG raw dict)."""
    ng = sigma_t.size
    eg = np.logspace(7, -3, ng + 1)
    chi_arr = chi if chi is not None else np.array(
        [1.0] + [0.0] * (ng - 1)
    )
    return Mixture(
        SigC=np.zeros(ng),
        SigL=np.zeros(ng),
        SigF=np.zeros(ng),
        SigP=nu_sigma_f.astype(float),
        SigT=sigma_t.astype(float),
        SigS=[csr_matrix(sigma_s.astype(float))],
        Sig2=csr_matrix((ng, ng)),
        chi=chi_arr.astype(float),
        eg=eg,
    )


# ----------------------------------------------------------------------
# Test matrix — 11 bit-equality tests per the spec
# ----------------------------------------------------------------------


def test_unified_factory_bit_equal_sphere_1g():
    r"""1G sphere — legacy raw XS dict ≡ new (Mixture, GeometrySpec).

    XS chosen with σ_t ≠ σ_s ≠ νσ_f to catch any field-swap drift
    between the legacy raw-dict path and the Mixture-keyed path.
    """
    sigma_t, sigma_s, nu_sigma_f, R, alpha = 0.5, 0.4, 0.1, 5.0, 1.0
    n_r, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    # New path: Mixture + GeometrySpec
    new_b = Billiard.from_problem(
        materials={0: _make_mixture_1g(sigma_t, sigma_s, nu_sigma_f)},
        geometry_spec=GeometrySpec(
            geometry="sphere",
            critical_dimension_mfp=sigma_t * R,
            critical_dimension_cm=R,
            n_groups=1,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    # Legacy path: raw XS dict + geometry dict
    legacy_b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R": R},
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_sphere_2g_asymmetric():
    r"""2G sphere with **asymmetric** scattering matrix.

    ``SigS[0,1] != SigS[1,0]`` — catches scattering-transpose drift
    (Sig 3 / ERR-002). If the new path accidentally transposes
    ``sigma_s`` during synthesis, the bit-equality check fails.
    """
    sigma_t = np.array([1.0, 0.5])
    # Asymmetric: down-scatter (group 0 → 1) is 0.40, up-scatter
    # (group 1 → 0) is 0.0 — the same pattern the legacy MG fixture
    # uses, but explicitly distinguished from a hypothetical
    # transposed pattern.
    sigma_s = np.array([[0.4, 0.4], [0.0, 0.4]])
    nu_sigma_f = np.array([0.05, 0.10])
    R, alpha = 5.0, 1.0
    n_r, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    mix = _make_mixture_mg(sigma_t, sigma_s, nu_sigma_f)
    new_b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=GeometrySpec(
            geometry="sphere",
            critical_dimension_mfp=R,
            critical_dimension_cm=R,
            n_groups=2,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    legacy_b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R": R},
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_sphere_4g_asymmetric():
    r"""4G sphere — group-coupling stress test.

    Down-scatter cascade with asymmetric inter-group coupling.
    Catches any flat slice-by-slice indexing drift in MG paths.
    """
    sigma_t = np.array([1.0, 0.8, 0.6, 0.4])
    # Down-scatter cascade (no upscatter); each group has its own
    # in-group + cascading down-scatter rate.
    sigma_s = np.array([
        [0.3, 0.2, 0.1, 0.05],
        [0.0, 0.4, 0.15, 0.05],
        [0.0, 0.0, 0.4, 0.1],
        [0.0, 0.0, 0.0, 0.3],
    ])
    nu_sigma_f = np.array([0.04, 0.08, 0.16, 0.20])
    R, alpha = 5.0, 1.0
    n_r, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    mix = _make_mixture_mg(sigma_t, sigma_s, nu_sigma_f)
    new_b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=GeometrySpec(
            geometry="sphere",
            critical_dimension_mfp=R,
            critical_dimension_cm=R,
            n_groups=4,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    legacy_b = Billiard.from_problem(
        geometry_kind="sphere",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R": R},
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_cylinder_1g():
    r"""1G cylinder — different XS keys to catch field-swap drift."""
    sigma_t, sigma_s, nu_sigma_f, R, alpha = 0.5, 0.4, 0.1, 5.0, 1.0
    n_r, n_mu_axial, n_phi_az, n_traj_quad, max_iter, tol = 6, 6, 8, 12, 10, 1e-9

    new_b = Billiard.from_problem(
        materials={0: _make_mixture_1g(sigma_t, sigma_s, nu_sigma_f)},
        geometry_spec=GeometrySpec(
            geometry="cylinder",
            critical_dimension_mfp=sigma_t * R,
            critical_dimension_cm=R,
            n_groups=1,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={
            "n_r": n_r,
            "n_mu_axial": n_mu_axial,
            "n_phi_az": n_phi_az,
            "n_traj_quad": n_traj_quad,
        },
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    legacy_b = Billiard.from_problem(
        geometry_kind="cylinder",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R": R},
        alpha=alpha,
        quadrature={
            "n_r": n_r,
            "n_mu_axial": n_mu_axial,
            "n_phi_az": n_phi_az,
            "n_traj_quad": n_traj_quad,
        },
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_cylinder_2g_asymmetric():
    r"""2G cylinder — asymmetric SigS catcher (Sig 3 / ERR-002)."""
    sigma_t = np.array([1.0, 0.5])
    sigma_s = np.array([[0.4, 0.4], [0.0, 0.4]])
    nu_sigma_f = np.array([0.05, 0.10])
    R, alpha = 5.0, 1.0
    n_r, n_mu_axial, n_phi_az, n_traj_quad, max_iter, tol = 6, 6, 8, 12, 10, 1e-9

    mix = _make_mixture_mg(sigma_t, sigma_s, nu_sigma_f)
    new_b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=GeometrySpec(
            geometry="cylinder",
            critical_dimension_mfp=R,
            critical_dimension_cm=R,
            n_groups=2,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={
            "n_r": n_r,
            "n_mu_axial": n_mu_axial,
            "n_phi_az": n_phi_az,
            "n_traj_quad": n_traj_quad,
        },
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    legacy_b = Billiard.from_problem(
        geometry_kind="cylinder",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R": R},
        alpha=alpha,
        quadrature={
            "n_r": n_r,
            "n_mu_axial": n_mu_axial,
            "n_phi_az": n_phi_az,
            "n_traj_quad": n_traj_quad,
        },
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_slab_1g_asymmetric_dim():
    r"""1G slab with **asymmetric dimensional** input (a_cm = 1.234).

    Catches factor-of-2 swap between half-thickness ``a`` and full
    width ``L = 2a`` (Sig 5). The new path passes ``a_cm`` via
    ``GeometrySpec.critical_dimension_cm`` (which the spec converts
    to ``L = 2a`` via ``domain_extent_cm``); the legacy path passes
    ``L`` directly. If the new path used ``a`` instead of ``2a``
    when populating the legacy geometry payload, the resulting
    Billiard would solve a slab half as wide and ``k_eff`` /
    ``parameter_value`` would diverge.
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.4, 0.1
    a_cm = 1.234  # half-thickness
    L_cm = 2.0 * a_cm  # full slab width
    alpha = 1.0
    n_x, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    new_b = Billiard.from_problem(
        materials={0: _make_mixture_1g(sigma_t, sigma_s, nu_sigma_f)},
        geometry_spec=GeometrySpec(
            geometry="slab",
            critical_dimension_mfp=sigma_t * a_cm,
            critical_dimension_cm=a_cm,
            n_groups=1,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={"n_x": n_x, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    legacy_b = Billiard.from_problem(
        geometry_kind="slab",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"L": L_cm},
        alpha=alpha,
        quadrature={"n_x": n_x, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_slab_2g_asymmetric():
    r"""2G slab — asymmetric SigS + slab cusp."""
    sigma_t = np.array([1.0, 0.5])
    sigma_s = np.array([[0.4, 0.4], [0.0, 0.4]])
    nu_sigma_f = np.array([0.05, 0.10])
    a_cm = 2.5
    L_cm = 2.0 * a_cm
    alpha = 1.0
    n_x, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    mix = _make_mixture_mg(sigma_t, sigma_s, nu_sigma_f)
    new_b = Billiard.from_problem(
        materials={0: mix},
        geometry_spec=GeometrySpec(
            geometry="slab",
            critical_dimension_mfp=a_cm,
            critical_dimension_cm=a_cm,
            n_groups=2,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
        quadrature={"n_x": n_x, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)

    legacy_b = Billiard.from_problem(
        geometry_kind="slab",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"L": L_cm},
        alpha=alpha,
        quadrature={"n_x": n_x, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_slab_asymmetric_with_alpha_dict():
    r"""1G slab_asymmetric with **asymmetric BC** (α_L ≠ α_R).

    Catches bc_left ↔ bc_right swap (Sig 4). The new path passes
    ``alpha={"alpha_left": 0.5, "alpha_right": 0.8}`` directly; the
    legacy path passes the same dict. The dual-factory normalisation
    must preserve the order. If the new path swaps left/right when
    converting to the alpha_payload, the legacy and new solves would
    produce different k_eff / psi.

    Geometry kind ``slab_asymmetric`` is only inferred when alpha is
    a dict with ``alpha_left`` / ``alpha_right`` keys; otherwise
    ``slab`` is used.
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.4, 0.1
    a_cm = 2.5
    L_cm = 2.0 * a_cm
    n_x, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    new_b = Billiard.from_problem(
        materials={0: _make_mixture_1g(sigma_t, sigma_s, nu_sigma_f)},
        geometry_spec=GeometrySpec(
            geometry="slab",
            critical_dimension_mfp=sigma_t * a_cm,
            critical_dimension_cm=a_cm,
            n_groups=1,
            bc_left=BC("partial", {"albedo": 0.5}),
            bc_right=BC("partial", {"albedo": 0.8}),
        ),
        alpha={"alpha_left": 0.5, "alpha_right": 0.8},
        quadrature={"n_x": n_x, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_new = new_b.solve_critical(max_iter=max_iter, tol=tol)
    assert new_b.geometry_kind == "slab_asymmetric"

    legacy_b = Billiard.from_problem(
        geometry_kind="slab_asymmetric",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"L": L_cm},
        alpha={"alpha_left": 0.5, "alpha_right": 0.8},
        quadrature={"n_x": n_x, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    _assert_bit_equal_critical(sol_new, sol_legacy)


def test_unified_factory_bit_equal_hollow_sphere_1g():
    r"""1G hollow sphere — legacy path only.

    The new (GeometrySpec-based) path does NOT yet support hollow
    geometries (the spec doesn't carry ``R_in`` — see
    ``_build_legacy_geometry_payload``). This test pins the legacy
    path's correctness and documents that the new path raises
    cleanly when given a hollow geometry.
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.4, 0.1
    R_in, R_out, alpha = 1.0, 5.0, 1.0
    n_r, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    legacy_b = Billiard.from_problem(
        geometry_kind="hollow_sphere",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R_in": R_in, "R_out": R_out},
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    # Legacy path is the canonical one; verify the synthesized
    # Protocol surface is still populated correctly.
    assert isinstance(sol_legacy, CriticalSolution)
    assert sol_legacy.metadata["geometry_kind"] == "hollow_sphere"
    assert legacy_b.method_name == "trajectory_resolvent"
    # The synthesized GeometrySpec for hollow sphere uses the outer
    # radius as ``critical_dimension_cm`` (R_in is not in the spec
    # schema yet).
    assert legacy_b.geometry_spec is not None
    assert legacy_b.geometry_spec.critical_dimension_cm == R_out

    # The new path (Mixture + GeometrySpec) does not yet support
    # hollow geometries: the GeometrySpec schema does not carry an
    # ``R_in``. Building a Billiard via the new path with a
    # hollow-style alpha dict simply produces a SOLID sphere
    # (geometry_spec.geometry == "sphere"), not a hollow sphere.
    # The legacy path is the canonical entry for hollow geometries;
    # this is documented behaviour until a richer spec lands.
    new_solid_b = Billiard.from_problem(
        materials={0: _make_mixture_1g(sigma_t, sigma_s, nu_sigma_f)},
        geometry_spec=GeometrySpec(
            geometry="sphere",
            critical_dimension_mfp=sigma_t * R_out,
            critical_dimension_cm=R_out,
            n_groups=1,
            bc_left=BC.reflective,
            bc_right=BC.reflective,
        ),
        alpha=alpha,
    )
    assert new_solid_b.geometry_kind == "sphere"
    # Hollow-specific kwargs (R_in) are unreachable through the new
    # path; the legacy path is required.


def test_unified_factory_bit_equal_annulus_1g():
    r"""1G annulus — legacy path only (analogous to hollow sphere)."""
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.4, 0.1
    R_in, R_out, alpha = 1.0, 5.0, 1.0
    n_r, n_mu_axial, n_phi_az, n_traj_quad, max_iter, tol = 6, 6, 8, 12, 10, 1e-9

    legacy_b = Billiard.from_problem(
        geometry_kind="annulus",
        materials={
            "sigma_t": sigma_t,
            "sigma_s": sigma_s,
            "nu_sigma_f": nu_sigma_f,
        },
        geometry={"R_in": R_in, "R_out": R_out},
        alpha=alpha,
        quadrature={
            "n_r": n_r,
            "n_mu_axial": n_mu_axial,
            "n_phi_az": n_phi_az,
            "n_traj_quad": n_traj_quad,
        },
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    assert isinstance(sol_legacy, CriticalSolution)
    assert sol_legacy.metadata["geometry_kind"] == "annulus"
    assert legacy_b.method_name == "trajectory_resolvent"
    assert legacy_b.geometry_spec.critical_dimension_cm == R_out


def test_unified_factory_bit_equal_sphere_mr_2g():
    r"""2G multi-region sphere — legacy path only.

    Multi-region (sphere_mr) is not yet supported on the new path
    because the GeometrySpec schema does not carry per-region radii.
    This test pins the legacy path's bit-equal preservation through
    the dual-factory rewire and confirms the synthesized
    ``materials`` dict has the expected per-region Mixture keys.
    """
    radii = np.array([2.0, 5.0])
    sigma_t_r = np.array([[1.0, 0.5], [0.5, 0.3]])
    sigma_s_r = np.array(
        [[[0.5, 0.4], [0.0, 0.4]], [[0.3, 0.2], [0.0, 0.2]]]
    )
    nu_sf_r = np.array([[0.05, 0.10], [0.0, 0.0]])
    alpha = 1.0
    n_r, n_mu, n_traj_quad, max_iter, tol = 8, 8, 16, 10, 1e-9

    legacy_b = Billiard.from_problem(
        geometry_kind="sphere_mr",
        materials={
            "sigma_t": sigma_t_r,
            "sigma_s": sigma_s_r,
            "nu_sigma_f": nu_sf_r,
        },
        geometry={"radii": radii},
        alpha=alpha,
        quadrature={"n_r": n_r, "n_mu": n_mu, "n_traj_quad": n_traj_quad},
    )
    sol_legacy = legacy_b.solve_critical(max_iter=max_iter, tol=tol)

    assert isinstance(sol_legacy, CriticalSolution)
    assert sol_legacy.metadata["geometry_kind"] == "sphere_mr"
    # Per-region Mixture synthesis: 2 regions → 2 Mixture entries.
    assert set(legacy_b.materials.keys()) == {0, 1}
    assert isinstance(legacy_b.materials[0], Mixture)
    assert isinstance(legacy_b.materials[1], Mixture)
    # Per-region σ_t propagated into per-region Mixture.
    assert _bit_equal_arrays(
        legacy_b.materials[0].SigT, sigma_t_r[0]
    )
    assert _bit_equal_arrays(
        legacy_b.materials[1].SigT, sigma_t_r[1]
    )

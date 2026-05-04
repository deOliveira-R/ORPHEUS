r"""Foundation tests for ``MomentSpace.from_problem`` Protocol conformance.

The test-architect spec § "File 3" identified that
:class:`~orpheus.derivations.continuous.fn_method.moment_space.MomentSpace`
already accepts the production-protocol input shape — a
``dict[int, Mixture]`` + :class:`GeometrySpec` pair — but the
foundation gate referenced in MomentSpace's docstring
(``tests.derivations.test_moment_space``) was never written. This
file fills that gap.

What this file pins
-------------------

1. **Bit-equality with the function-level API** — direct calls to
   :func:`solve_fn_slab_bare_critical` /
   :func:`solve_fn_sphere_bare_critical` /
   :func:`compute_kinf_*` produce the same float results as
   ``MomentSpace.from_problem(...).solve_critical()``. Verified via
   :meth:`float.hex` exact-bit comparison.

2. **Keyword + positional argument routes match** — the factory
   accepts ``materials`` and ``geometry`` either positionally or by
   keyword; both routes produce equivalent results.

3. **Protocol surface (``materials`` / ``geometry_spec`` /
   ``method_name``)** — round-trip preservation of the
   production-protocol input.

4. **Asymmetric chi MG case** (Sig 3 catcher) — a 2G mixture with
   chi-asymmetric fission spectrum produces the same ``k_inf`` as
   :func:`compute_kinf_2g_general` direct evaluation.

Foundation-tier
---------------

Module-level ``pytestmark = [pytest.mark.foundation]`` (no
``verifies(...)``); these are software invariants.

References
----------

* :class:`MomentSpace` docstring § "Notes on bit-equality" — names
  this file as the foundation gate that pins the bit-equality
  invariant.
* :doc:`/skills/algebra-of-record` § "Branch 1 / Branch 2 separation"
  — bit-equality preservation is what guarantees the F_N moment
  space's Branch-1 ⊥ Branch-2 cross-checks remain valid.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.geometry_spec import GeometrySpec
from orpheus.derivations.common.solution_types import CriticalSolution
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.fn_method.multi_group.k_inf import (
    compute_kinf_1g,
    compute_kinf_2g_general,
)
from orpheus.derivations.continuous.fn_method.slab import (
    solve_fn_slab_bare_critical,
)
from orpheus.derivations.continuous.fn_method.sphere import (
    solve_fn_sphere_bare_critical,
)
from orpheus.geometry.mesh import BC


pytestmark = [pytest.mark.foundation]


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _bit_equal_floats(a: float, b: float) -> bool:
    """Strict bit-for-bit float equality."""
    return float(a).hex() == float(b).hex()


def _make_1g_mixture(c: float, sigma_t: float = 1.0) -> Mixture:
    r"""Build a 1G Mixture with mean-secondaries c.

    Splits ``c·sigma_t`` into ``sigma_s + nu_sigma_f`` such that
    ``(sigma_s + nu_sigma_f) / sigma_t == c`` to **machine precision**
    (NOT just ``c·sigma_t·0.7 + c·sigma_t·0.3``, which drifts by 1 ULP
    from c due to FP addition rounding). To preserve bit-equality with
    direct function-level F_N calls (which take ``c`` as a scalar
    input), the mixture's reconstructed c MUST be IEEE-754 identical
    to the input.
    """
    nu_sigma_f = 0.3 * c * sigma_t
    sigma_s = c * sigma_t - nu_sigma_f  # exact: x - y for FP-rounded y
    sig_c = sigma_t - sigma_s - nu_sigma_f  # absorption
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sig_c]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
    )


# ----------------------------------------------------------------------
# 1. Construction routes — keyword vs positional
# ----------------------------------------------------------------------


def test_moment_space_from_problem_keyword_only_compat():
    r"""Keyword and positional argument routes produce equivalent objects.

    The factory's signature is
    ``from_problem(materials, geometry, *, fn_order=..., flux_reconstruction=...)``.
    Calling with ``materials=`` / ``geometry=`` keywords vs. positional
    args must yield equivalent MomentSpace instances (same field
    contents, same Protocol surface).
    """
    mix = _make_1g_mixture(c=1.30)
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=0.93772556,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )

    ms_kw = MomentSpace.from_problem(
        materials={0: mix}, geometry=spec, fn_order=8,
    )
    ms_pos = MomentSpace.from_problem(
        {0: mix}, spec, fn_order=8,
    )

    assert ms_kw.fn_order == ms_pos.fn_order
    assert ms_kw.geometry is ms_pos.geometry  # same instance
    assert ms_kw.materials[0] is ms_pos.materials[0]
    assert ms_kw.method_name == ms_pos.method_name == "fn_method"


# ----------------------------------------------------------------------
# 2. Bit-equality with function-level API (slab + sphere)
# ----------------------------------------------------------------------


def test_moment_space_solve_critical_slab_bit_equal_function_api():
    r"""MomentSpace slab ≡ ``solve_fn_slab_bare_critical``.

    Pins the class facade's bit-equal preservation against direct
    function-level calls. The class is documented as a thin facade;
    this test enforces that.
    """
    c = 1.30
    fn_order = 10

    mix = _make_1g_mixture(c=c)
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.94,
        critical_dimension_cm=0.94,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=spec, fn_order=fn_order,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol_class = ms.solve_critical()
        res_func = solve_fn_slab_bare_critical(c=c, n_modes=fn_order)

    assert isinstance(sol_class, CriticalSolution)
    assert sol_class.eigenvalue_kind == "k_eff"
    assert sol_class.parameter_kind == "half_thickness_mfp"
    assert _bit_equal_floats(
        sol_class.parameter_value, res_func.a_critical_mfp
    ), (
        f"slab bit-equality drift: class={sol_class.parameter_value!r} "
        f"function={res_func.a_critical_mfp!r}"
    )


def test_moment_space_solve_critical_sphere_bit_equal_function_api():
    r"""MomentSpace sphere ≡ ``solve_fn_sphere_bare_critical``."""
    c = 1.30
    fn_order = 10

    mix = _make_1g_mixture(c=c)
    spec = GeometrySpec(
        geometry="sphere",
        critical_dimension_mfp=2.4,
        critical_dimension_cm=2.4,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    )
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=spec, fn_order=fn_order,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol_class = ms.solve_critical()
        res_func = solve_fn_sphere_bare_critical(c=c, n_modes=fn_order)

    assert sol_class.parameter_kind == "radius_mfp"
    assert _bit_equal_floats(
        sol_class.parameter_value, res_func.R_critical_mfp
    ), (
        f"sphere bit-equality drift: class="
        f"{sol_class.parameter_value!r} function="
        f"{res_func.R_critical_mfp!r}"
    )


# ----------------------------------------------------------------------
# 3. Bit-equality on infinite-medium k_inf (1G + 2G asymmetric chi)
# ----------------------------------------------------------------------


def test_moment_space_solve_critical_infinite_1g():
    r"""MomentSpace infinite 1G ≡ :func:`compute_kinf_1g`."""
    sigma_t, sigma_s, nu_sigma_f = 1.0, 0.7, 0.4
    mix = make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sigma_t - sigma_s - nu_sigma_f]),
        sig_f=np.array([nu_sigma_f / 2.0]),
        nu=np.array([2.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s]]),
    )
    spec = GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    )
    ms = MomentSpace.from_problem(materials={0: mix}, geometry=spec)
    sol = ms.solve_critical()

    expected = compute_kinf_1g(sigma_t, sigma_s, nu_sigma_f)
    assert sol.eigenvalue_kind == "k_inf"
    assert _bit_equal_floats(sol.eigenvalue, expected), (
        f"1G k_inf drift: ms={sol.eigenvalue!r} expected={expected!r}"
    )


def test_moment_space_solve_critical_infinite_2g_asymmetric_chi():
    r"""2G k_inf with **chi-asymmetric** fission spectrum.

    Stress test for the 2G branch (Sig 3 catch — group-coupling
    drift). The chi vector is non-trivial (``chi[0] = 0.6, chi[1] =
    0.4``) so a swap of group indices in ``compute_kinf_2g_general``
    would land at a different k_inf.
    """
    sigma_t = np.array([1.0, 0.5])
    sigma_s = np.array([[0.4, 0.4], [0.0, 0.4]])
    nu_sigma_f = np.array([0.05, 0.10])
    chi = np.array([0.6, 0.4])

    mix = Mixture(
        SigC=np.zeros(2),
        SigL=np.zeros(2),
        SigF=np.zeros(2),
        SigP=nu_sigma_f.copy(),
        SigT=sigma_t.copy(),
        SigS=[csr_matrix(sigma_s)],
        Sig2=csr_matrix((2, 2)),
        chi=chi.copy(),
        eg=np.logspace(7, -3, 3),
    )
    spec = GeometrySpec(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    )
    ms = MomentSpace.from_problem(materials={0: mix}, geometry=spec)
    sol = ms.solve_critical()

    sigma_s_arr = np.asarray(sigma_s)
    expected = compute_kinf_2g_general(
        sigma_t.reshape(2),
        sigma_s_arr.reshape(2, 2),
        nu_sigma_f.reshape(2),
        chi.reshape(2),
    )
    assert sol.eigenvalue_kind == "k_inf"
    assert _bit_equal_floats(sol.eigenvalue, expected), (
        f"2G chi-asymmetric k_inf drift: ms={sol.eigenvalue!r} "
        f"expected={expected!r}"
    )


# ----------------------------------------------------------------------
# 4. Protocol surface round-trip
# ----------------------------------------------------------------------


def test_moment_space_method_name_and_properties():
    r"""``method_name`` + Protocol-surface round-trip.

    The :class:`TransportSolver` Protocol requires:

    * ``materials`` returns the production-protocol ``dict[int, Mixture]``.
    * ``geometry_spec`` returns the :class:`GeometrySpec` instance.
    * ``method_name`` returns ``"fn_method"``.

    This test pins all three after construction.
    """
    mix = _make_1g_mixture(c=1.30)
    spec = GeometrySpec(
        geometry="slab",
        critical_dimension_mfp=0.94,
        critical_dimension_cm=0.94,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    )
    ms = MomentSpace.from_problem(
        materials={0: mix}, geometry=spec, fn_order=8,
    )

    # method_name
    assert ms.method_name == "fn_method"

    # materials — production-protocol shape preserved.
    assert isinstance(ms.materials, dict)
    assert set(ms.materials.keys()) == {0}
    assert ms.materials[0] is mix

    # geometry_spec — Protocol alias for the .geometry field.
    assert ms.geometry_spec is spec
    assert ms.geometry_spec is ms.geometry  # alias is identity-stable
    assert isinstance(ms.geometry_spec, GeometrySpec)
    assert ms.geometry_spec.geometry == "slab"

r"""V&V test suite for the spherical-harmonic space / basis / projection algebra.

The five identities from §P1.5 of the moment-space + layering plan
become executable tests here. ALL tests in this file carry
``@pytest.mark.catches("ERR-039")`` — they pin the post-fix contract
that the (R, Π^T, Π^*) operators are separately typed with
mathematically distinct semantics.

The ``@pytest.mark.verifies(<label>)`` markers reference equation
labels that ship in the Sphinx ``docs/theory/spherical_harmonics.rst``
page under P1.6. The labels follow the ``sh-`` prefix established by
the test-architect's verification plan so that Phase 2's ``dual-`` /
``tensor-`` / ``sum-`` labels sit alongside in the same namespace.

Pillars (per ``vv-principles`` skill):

* ``test_basis_mass_matrix_against_lebedev`` — semi-analytical (the
  Lebedev rule integrates spherical harmonics exactly up to its
  declared degree; the resulting Gram matrix is the discrete
  semi-analytical reference for the continuous identity
  :math:`\langle Y_\ell^m, Y_{\ell'}^{m'}\rangle_{L^2(S^2)} =
  (4\pi/(2\ell+1)) \delta`).
* ``test_R_equals_2l_plus_1_times_S0`` / ``test_pi_R_is_4pi_identity``
  / ``test_H_equals_g_C_times_S0`` / ``test_T_carries_w_n`` —
  closed-form (the identities are algebraic, not numerical).
* ``test_*_codomain_*`` / ``test_*_roundtrip`` / ``test_*_equality_*``
  — type-system / construction-API checks (L0 software invariants).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.basis import SphericalHarmonicBasis
from orpheus.numerics.projection import (
    HarmonicMomentReconstruction,
    MomentProjection,
)
from orpheus.numerics.quadrature import (
    lebedev_sphere,
    product_mu_phi,
)
from orpheus.numerics.spaces import SphericalHarmonicSpace


# ─────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────


@pytest.fixture
def lebedev_L_pair():
    """Pair a Lebedev rule with an L that the rule integrates exactly.

    The Lebedev rule of order ``order`` integrates
    :math:`Y_\\ell^m Y_{\\ell'}^{m'}` exactly for
    :math:`\\ell + \\ell' \\le \\mathrm{order}`. For L=3, need
    order >= 6 — pick 13 for headroom.
    """
    measure = lebedev_sphere(13)
    L = 3
    return measure, L


def _mask_non_existent_m(c: np.ndarray, L: int) -> np.ndarray:
    """Zero out the |m| > l padding entries of a (L+1, 2L+1, ...) array.

    The :math:`(L+1, 2L+1)` storage shape leaves padded slots that the
    addition-theorem identity assumes to be zero; tests use this helper
    to construct "band-limited" inputs explicitly.
    """
    out = c.copy()
    for l_idx in range(L + 1):
        out[l_idx, 2 * l_idx + 1 :] = 0.0
    return out


# ─────────────────────────────────────────────────────────────────────
# B.1 — the five P1.5 identities (one test each)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("sh-space-metric")
def test_space_inner_product_weights_equal_4pi_over_2l_plus_1():
    r"""``SphericalHarmonicSpace.from_L(L).inner_product_weights[\ell] == 4\pi/(2\ell+1)``.

    The Gram-matrix diagonal :math:`g_C` lives in exactly one place —
    on the space — and the padded ``(L+1, 2L+1)`` layout matches the
    :class:`~orpheus.sn.harmonic_moment_field.HarmonicMomentField`
    storage convention (row :math:`\ell` carries
    :math:`4\pi/(2\ell+1)` in the :math:`2\ell+1` valid slots, zero in
    the :math:`|m|>\ell` padding).
    """
    L = 4
    space = SphericalHarmonicSpace.from_L(L)
    expected_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    weights = space.inner_product_weights
    assert weights.shape == (L + 1, 2 * L + 1)
    for ell in range(L + 1):
        np.testing.assert_allclose(
            weights[ell, : 2 * ell + 1],
            expected_per_ell[ell],
            rtol=1e-15,
        )
        np.testing.assert_array_equal(weights[ell, 2 * ell + 1 :], 0.0)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("real-sh-discrete-orthogonality")
def test_basis_mass_matrix_against_lebedev(lebedev_L_pair):
    r"""``SphericalHarmonicBasis.mass_matrix(lebedev_measure)`` ≈ the theoretical metric.

    Pins the SH convention against a structurally-independent
    discrete-orthogonality computation: the Lebedev rule of degree
    :math:`\ge 2L` integrates :math:`Y_\ell^m Y_{\ell'}^{m'}` exactly,
    giving the Gram matrix :math:`(4\pi/(2\ell+1))\delta_{\ell\ell'}
    \delta_{m m'}` on the diagonal of the ``(L+1, 2L+1, L+1, 2L+1)``
    4-tensor and zero (to FP roundoff) on the off-diagonal.
    """
    measure, L = lebedev_L_pair
    basis = SphericalHarmonicBasis(L=L)
    G = basis.mass_matrix(measure)  # (L+1, 2L+1, L+1, 2L+1)
    assert G.shape == (L + 1, 2 * L + 1, L + 1, 2 * L + 1)

    expected_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    # Diagonal: G[l, l+m, l, l+m] == 4π/(2l+1) for |m| <= l.
    for ell in range(L + 1):
        for m_off in range(2 * ell + 1):
            actual = G[ell, m_off, ell, m_off]
            np.testing.assert_allclose(
                actual, expected_per_ell[ell],
                rtol=1e-12,
                err_msg=f"ell={ell}, m_off={m_off}",
            )

    # Off-diagonal (ell != ell' or m != m'): ≈ 0 to quadrature
    # precision.
    for ell in range(L + 1):
        for ell_p in range(L + 1):
            for m_off in range(2 * ell + 1):
                for m_off_p in range(2 * ell_p + 1):
                    if (ell, m_off) == (ell_p, m_off_p):
                        continue
                    np.testing.assert_allclose(
                        G[ell, m_off, ell_p, m_off_p],
                        0.0,
                        atol=1e-12,
                    )


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("sh-addition-theorem-reconstruction")
def test_R_equals_2l_plus_1_times_S0(lebedev_L_pair):
    r"""``R.apply(c) == (2\ell+1) \cdot S_0(c)`` for random band-limited c.

    The addition-theorem reconstruction :math:`R` differs from the
    naked synthesis :math:`S_0` by the per-:math:`\ell` factor
    :math:`(2\ell+1)`. This pins the ERR-039 distinction at the
    operator-construction level: ``R`` is built via the canonical
    :meth:`~HarmonicMomentReconstruction.from_spherical_harmonic_space`
    constructor that sources :math:`(2\ell+1)` from
    :attr:`SphericalHarmonicSpace.addition_theorem_factor`.
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    R = HarmonicMomentReconstruction.from_spherical_harmonic_space(space, Y)

    rng = np.random.default_rng(seed=2026)
    c = _mask_non_existent_m(rng.standard_normal((L + 1, 2 * L + 1)), L)

    actual = R.apply(c)
    # (2ℓ+1) · S_0(c) per the addition-theorem formula.
    expected = (
        basis.synthesize(c * basis.addition_theorem_factor[:, None], measure.nodes)
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-14)

    # Structurally-independent cross-check on unit vectors: for
    # c = e_{ℓ₀, m_off₀} (single nonzero entry at one (ℓ, m) slot),
    # the einsum collapses to a single multiplication per ordinate
    # (no accumulation), so the (2ℓ+1) literal carried by R can be
    # read off bit-identically as the per-ordinate scalar factor on
    # the column ``Y[:, ℓ, m_off]``.  No FP-non-associativity room
    # — this is the structural independence per `lessons-L11`.
    for ell_0 in range(L + 1):
        for m_off_0 in range(2 * ell_0 + 1):
            e = np.zeros((L + 1, 2 * L + 1))
            e[ell_0, m_off_0] = 1.0
            R_e = R.apply(e)
            expected_e = (2.0 * ell_0 + 1.0) * Y[:, ell_0, m_off_0]
            np.testing.assert_array_equal(
                R_e, expected_e,
                err_msg=f"R(e_{{ell={ell_0}, m_off={m_off_0}}}) "
                        f"should equal (2ℓ+1) · Y[:, ℓ, m_off] bit-identically",
            )


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("pi-r-equals-4pi-i")
def test_pi_R_is_4pi_identity_on_band_limited(lebedev_L_pair):
    r""":math:`\Pi \cdot R = 4\pi \cdot I` on the band-limited coefficient space.

    Sibling of
    ``tests/numerics/test_projection_operators.py::TestGalerkinIdempotencyOnLebedev::test_pi_R_is_identity_on_band_limited``;
    this version constructs the operators through the new
    :class:`SphericalHarmonicSpace` API and pins the genuine
    :math:`\Pi R = 4\pi I` identity (NOT the broken ``Π R == I`` that
    the retired :meth:`assert_galerkin_idempotency` was checking — see
    P1.6).
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )
    R = HarmonicMomentReconstruction.from_spherical_harmonic_space(space, Y)

    rng = np.random.default_rng(seed=42)
    c = _mask_non_existent_m(rng.standard_normal((L + 1, 2 * L + 1)), L)

    out = M.apply(R.apply(c))
    np.testing.assert_allclose(out, 4.0 * np.pi * c, rtol=1e-10, atol=1e-12)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("hilbert-adjoint-equals-metric-times-S0")
def test_H_equals_g_C_times_S0(lebedev_L_pair):
    r"""``M.H`` computed generically equals :math:`g_C \cdot S_0(c)`.

    The :attr:`MomentProjection.H` property routes through the
    generic :class:`~orpheus.numerics.operator._AdjointOperator`
    wrapper, which composes :math:`(1/w_V) \cdot \Pi^\top(w_W \cdot c)`
    using :attr:`MomentProjection.domain`'s quadrature weights as
    :math:`w_V` and :attr:`MomentProjection.codomain`'s SH-Gram
    diagonal as :math:`w_W`. The result is the metric-weighted naked
    synthesis :math:`g_C \cdot S_0(c)` — the W-weighted Hilbert
    adjoint :math:`\Pi^*`.

    This is the ERR-039 endpoint: the metric, the transpose, and the
    Hilbert adjoint are SEPARATELY TYPED and their composition falls
    out of the generic machinery without prose warnings.
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )

    rng = np.random.default_rng(seed=99)
    c = _mask_non_existent_m(rng.standard_normal((L + 1, 2 * L + 1)), L)

    actual = M.H.apply(c)
    # Expected: g_C · S_0(c) with g_C^ℓ = 4π/(2ℓ+1).
    g_C_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    c_scaled = c * g_C_per_ell[:, None]
    expected = basis.synthesize(c_scaled, measure.nodes)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# B.2 — constructor / API surface
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("real-sh-discrete-orthogonality")
@pytest.mark.parametrize(
    "quadrature_factory",
    [
        pytest.param(lambda: lebedev_sphere(7), id="lebedev_7"),
        pytest.param(lambda: lebedev_sphere(13), id="lebedev_13"),
        pytest.param(
            lambda: product_mu_phi(6, 8)[0],
            id="product_mu_phi_6x8",
        ),
        pytest.param(
            lambda: product_mu_phi(8, 16)[0],
            id="product_mu_phi_8x16",
        ),
    ],
)
def test_mass_matrix_under_multiple_quadratures(quadrature_factory):
    r""":meth:`SphericalHarmonicBasis.mass_matrix` is exact across SH-degree-sufficient quadratures.

    For :math:`L=2` and a quadrature that integrates the
    :math:`Y_\ell^m Y_{\ell'}^{m'}` products up to total degree
    :math:`\ge 2L = 4`, the discrete Gram diagonal equals the
    theoretical :math:`4\pi/(2\ell+1)` to FP precision. Lebedev rules
    at order :math:`\ge 2L` and product (Gauss-Legendre × equispaced
    Chebyshev-equivalent) rules with sufficient
    :math:`(n_\mu, n_\phi)` both satisfy this.

    Level-symmetric :math:`S_N` rules are deliberately omitted: they
    are optimised for moment integration in transport (the integrals
    that arise in the SN transport equation when reconstructing the
    angular-flux moments), NOT for arbitrary
    :math:`Y_\ell^m Y_{\ell'}^{m'}` products beyond degree 1. At
    :math:`L=2`, LS_8 has a 24% diagonal error and no LS order makes
    it exact — the rule's design objective is structurally different.
    """
    L = 2
    basis = SphericalHarmonicBasis(L=L)
    measure = quadrature_factory()
    G = basis.mass_matrix(measure)
    expected_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    for ell in range(L + 1):
        for m_off in range(2 * ell + 1):
            actual = G[ell, m_off, ell, m_off]
            np.testing.assert_allclose(
                actual, expected_per_ell[ell],
                rtol=1e-12,
                err_msg=f"ell={ell}, m_off={m_off}",
            )


@pytest.mark.l0
def test_moment_projection_codomain_is_spherical_harmonic_space():
    r""":attr:`MomentProjection.codomain` returns a typed :class:`SphericalHarmonicSpace`.

    Type-level guarantee that ``M.H`` composition via the generic
    adjoint machinery finds the SH metric correctly. The equality
    convention is ``(name, shape)``: two SphericalHarmonicSpaces of
    matching :math:`L` compare equal.
    """
    measure = lebedev_sphere(7)
    L = 2
    basis = SphericalHarmonicBasis(L=L)
    space = SphericalHarmonicSpace.from_L(L)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=basis.evaluate(measure.nodes),
    )
    cod = M.codomain
    assert isinstance(cod, SphericalHarmonicSpace)
    assert cod.L == L
    assert cod.shape == (L + 1, 2 * L + 1)
    assert cod == SphericalHarmonicSpace.from_L(L)

    # ``range`` aliases ``codomain`` (pre-P3.5 framework name read by
    # the generic _AdjointOperator).
    assert M.range == M.codomain


@pytest.mark.l0
def test_from_spherical_harmonic_space_roundtrip():
    r"""Construction via the new classmethod is bit-identical to the legacy path.

    Both :meth:`MomentProjection.from_spherical_harmonic_space` and
    :meth:`HarmonicMomentReconstruction.from_spherical_harmonic_space`
    produce operators whose :meth:`apply` outputs are bit-identical to
    the legacy ``(weights, Y, L)`` constructor and ``from_Y(Y)``
    respectively. The migration is API-only; no numerical change.
    """
    measure = lebedev_sphere(13)
    L = 3
    basis = SphericalHarmonicBasis(L=L)
    space = SphericalHarmonicSpace.from_L(L)
    Y = basis.evaluate(measure.nodes)

    # New API
    M_new = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )
    R_new = HarmonicMomentReconstruction.from_spherical_harmonic_space(space, Y)

    # Legacy API (preserved as shim per plan §P1.3)
    M_legacy = MomentProjection(weights=measure.weights, Y=Y, L=L)
    R_legacy = HarmonicMomentReconstruction.from_Y(Y)

    rng = np.random.default_rng(seed=12345)
    psi = rng.standard_normal(measure.n_points)
    c = _mask_non_existent_m(rng.standard_normal((L + 1, 2 * L + 1)), L)

    np.testing.assert_array_equal(M_new.apply(psi), M_legacy.apply(psi))
    np.testing.assert_array_equal(R_new.apply(c), R_legacy.apply(c))


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
@pytest.mark.verifies("moment-projection-transpose-T", "hilbert-adjoint-equals-metric-times-S0")
def test_T_carries_w_n_and_H_carries_g_C(lebedev_L_pair):
    r"""``M.apply_transpose`` carries :math:`w_n`; ``M.H`` carries :math:`g_C`.

    Direct pin of the post-P1.4 contract: the two operators differ by
    the per-:math:`\ell` factor :math:`g_C / w_n` (which is not a
    proper scalar — it lives in different axes), and ERR-039's
    original confusion no longer arises because the two are typed
    distinctly:

    .. code-block:: python

        M.apply_transpose(c)  → w_n · S_0(c)      # representation transpose
        M.H.apply(c)          → g_C · S_0(c)      # Hilbert adjoint
    """
    measure, L = lebedev_L_pair
    space = SphericalHarmonicSpace.from_L(L)
    basis = SphericalHarmonicBasis(L=L)
    Y = basis.evaluate(measure.nodes)
    M = MomentProjection.from_spherical_harmonic_space(
        space, weights=measure.weights, Y=Y,
    )

    rng = np.random.default_rng(seed=7)
    psi = rng.standard_normal(measure.n_points)
    c = _mask_non_existent_m(rng.standard_normal((L + 1, 2 * L + 1)), L)

    # ⟨Π ψ, c⟩  — Euclidean on coefficient space
    Mpsi = M.apply(psi)
    lhs_euclidean = float(np.sum(Mpsi * c))

    # ⟨ψ, Π^T c⟩_V_Euclidean — Π^T already carries w_n
    rhs_T = float(np.sum(psi * M.apply_transpose(c)))
    np.testing.assert_allclose(lhs_euclidean, rhs_T, rtol=1e-12, atol=1e-14)

    # ⟨Π ψ, c⟩_C  — coefficient inner product with g_C metric
    g_C_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
    c_in_C = c * g_C_per_ell[:, None]
    lhs_metric = float(np.sum(Mpsi * c_in_C))

    # ⟨ψ, Π^* c⟩_V_W  — angular inner product with W metric
    H_c = M.H.apply(c)
    rhs_H = float(np.sum(measure.weights * psi * H_c))
    np.testing.assert_allclose(lhs_metric, rhs_H, rtol=1e-12, atol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# B.3 — forward-compatibility test that paves the way for Phase 2
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_spherical_harmonic_space_equality_by_name_shape():
    r"""SphericalHarmonicSpace equality follows :class:`FunctionSpace`'s ``(name, shape)`` convention.

    Two ``SphericalHarmonicSpace.from_L(L)`` instances with the same
    :math:`L` produce equal :class:`SphericalHarmonicSpace` objects
    even when their ``inner_product_weights`` arrays are distinct
    ``ndarray`` allocations. This convention paves the way for Phase 2's
    :class:`DualSpace` /
    :class:`TensorProductSpace` /
    :class:`DirectSumSpace` constructors to compose
    :class:`SphericalHarmonicSpace` instances without spurious
    inequalities from FP-level metric differences.
    """
    a = SphericalHarmonicSpace.from_L(3)
    b = SphericalHarmonicSpace.from_L(3)
    c = SphericalHarmonicSpace.from_L(4)

    assert a == b
    assert a != c
    assert hash(a) == hash(b)
    assert hash(a) != hash(c)

    # Cross-class equality with a bare FunctionSpace carrying the same
    # (name, shape) — supports the "equal-shape spaces compare equal"
    # invariant under DualSpace / TensorProductSpace introduction in
    # Phase 2.
    from orpheus.numerics.space import FunctionSpace
    bare = FunctionSpace(name="spherical_harmonic_space", shape=(4, 7))
    assert a == bare

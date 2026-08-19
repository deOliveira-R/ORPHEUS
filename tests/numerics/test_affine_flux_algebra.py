r"""Foundation tests for the #208/#201 affine flux algebra.

Flux states form an **affine space** :math:`A` over a distinct **difference
vector space** :math:`V` (the displacements). This module pins the torsor
algebra across all four flux-state leaves (``AngularFlux``, ``ScalarFlux``,
``HarmonicMomentFlux``, ``AngularBoundaryFlux``):

* :math:`\psi_2 \ominus \psi_1 \to \Delta\psi \in V` (the mint),
* :math:`\psi \oplus \Delta\psi \to \psi' \in A` (the torsor action),
* :math:`\psi_1 \oplus \psi_2 \to \bot` (no origin — the #201 gate),
* :math:`\sum_i\lambda_i\psi_i,\ \sum\lambda_i=1 \to \psi` (affine combination),

(The displacement diagnostics that lived here migrated at campaign 1 CS3,
2026-08-19: ρ / true-error to :mod:`tests.numerics.test_iteration_record`
[the record derives them from ``increment_norms``], the per-entry map to
:mod:`tests.numerics.test_field` [``Field.where_largest``].)

These are ``foundation`` tests — software invariants, not physics ``:label:``
claims. The structural ground is numpy ``+``/``-`` on the same buffer (no
solver). Per ``-O`` discipline every assertion is ``np.testing.assert_*`` or an
explicit ``raise`` (bare ``assert`` is stripped under ``-O``).

See ``.claude/agent-memory/test-architect/issue_208_affine_algebra_verification_spec.md``
and ``.claude/plans/issue_208_residual_typing_closeout.md`` §3a.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.displacements import Displacement
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from tests.sn._test_helpers import placeholder_materials


pytestmark = pytest.mark.foundation

_MOMENT_L = 1
_LEAVES = ["angular", "scalar", "moment", "boundary"]


def _mesh() -> SNMesh:
    m = Mesh1D(
        edges=np.linspace(0.0, 2.0, 5), mat_ids=np.zeros(4, dtype=int),
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    return SNMesh(m, Quadrature.gauss_legendre(n_ordinates=4), placeholder_materials())


@pytest.fixture
def mesh() -> SNMesh:
    return _mesh()


@pytest.fixture
def mesh2() -> SNMesh:
    return _mesh()  # distinct instance, structurally identical


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(208)


def _make_flux(leaf: str, m: SNMesh, rng: np.random.Generator):
    """A flux leaf with DISTINCT random values (NOT flat — Mode-9)."""
    if leaf == "angular":
        return AngularFlux.from_mesh(
            rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape)), m,
        )
    if leaf == "scalar":
        return ScalarFlux.from_mesh(rng.standard_normal((m.ng, *m.spatial_shape)), m)
    if leaf == "moment":
        # HarmonicMomentFlux is rank-d like the other leaves: its space is
        # (L+1, 2L+1, ng, *spatial) — rank-1 (nx,) on a 1-D mesh, no phantom ny.
        shape = (_MOMENT_L + 1, 2 * _MOMENT_L + 1, m.ng, *m.spatial_shape)
        return HarmonicMomentFlux.from_mesh_and_L(
            rng.standard_normal(shape), m, _MOMENT_L,
        )
    if leaf == "boundary":
        zero = AngularBoundaryFlux.zeros_on(m)
        return AngularBoundaryFlux.from_mesh(rng.standard_normal(zero.values.shape), m)
    raise ValueError(leaf)


def _displacement_type(flux):
    """Rename-safe indirection: the class ``flux ⊖ flux`` mints."""
    return type(flux - flux)


# ── 1. Torsor round-trip + vector-space axioms (the keystone) ────────────


@pytest.mark.parametrize("leaf", _LEAVES)
def test_mint_stores_raw_difference_and_torsor_round_trips(leaf, mesh, rng):
    r"""ψ₂ ⊖ ψ₁ stores the RAW difference (bit-exact); ψ₁ ⊕ (ψ₂ ⊖ ψ₁) recovers
    ψ₂ to FP rounding.

    Two claims. (1) the mint is the literal single subtraction
    ``d.values == ψ₂.values − ψ₁.values`` (bit-exact) and the torsor is the
    literal single addition ``(ψ₁ ⊕ d).values == ψ₁.values + d.values``
    (bit-exact) — these are the affine algebra's defining operations. (2) the
    composed round-trip ``ψ₁ ⊕ (ψ₂ ⊖ ψ₁) ≈ ψ₂`` is exact ONLY up to IEEE-754
    rounding (``a + (b − a) ≠ b`` bit-for-bit for arbitrary a, b — the
    subtraction rounds), so it is asserted to a few ULP, NOT bit-exact. The
    result type MUST be the displacement (mint) then the flux (action)."""
    psi1, psi2 = _make_flux(leaf, mesh, rng), _make_flux(leaf, mesh, rng)
    d = psi2 - psi1
    if type(d) is type(psi1):
        raise AssertionError(
            f"{leaf}: ψ₂⊖ψ₁ returned the FLUX type {type(d).__name__}, not a "
            f"displacement — the affine-axiom violation #208 fixes."
        )
    if not isinstance(d, Displacement):
        raise AssertionError(f"{leaf}: ψ₂⊖ψ₁ is not a Displacement")
    # (1) mint stores the raw difference; torsor adds the raw buffer (bit-exact)
    np.testing.assert_array_equal(d.values, psi2.values - psi1.values)
    recovered = psi1 + d
    if type(recovered) is not type(psi1):
        raise AssertionError(
            f"{leaf}: ψ₁⊕d returned {type(recovered).__name__}, not the flux "
            f"type {type(psi1).__name__} — torsor action broken."
        )
    np.testing.assert_array_equal(recovered.values, psi1.values + d.values)
    # (2) the COMPOSED round-trip is exact up to FP rounding (NOT bit-exact)
    np.testing.assert_array_almost_equal_nulp(recovered.values, psi2.values, nulp=8)


@pytest.mark.parametrize("leaf", _LEAVES)
def test_displacement_telescoping(leaf, mesh, rng):
    r"""(ψ₂ ⊖ ψ₁) + (ψ₃ ⊖ ψ₂) == (ψ₃ ⊖ ψ₁) — path-independence of V.

    Bit-exact: both sides reduce to ``ψ₃.values − ψ₁.values`` once the
    intermediate ``ψ₂`` cancels (the left side is
    ``(ψ₂−ψ₁) + (ψ₃−ψ₂)``; FP-associativity makes this NOT generally
    bit-equal to ``ψ₃−ψ₁``, so assert to a few ULP)."""
    psi1, psi2, psi3 = (_make_flux(leaf, mesh, rng) for _ in range(3))
    lhs = (psi2 - psi1) + (psi3 - psi2)
    rhs = psi3 - psi1
    if type(lhs) is not _displacement_type(psi1):
        raise AssertionError(f"{leaf}: disp + disp did not return a displacement")
    np.testing.assert_array_almost_equal_nulp(lhs.values, rhs.values, nulp=8)


@pytest.mark.parametrize("leaf", _LEAVES)
def test_displacement_scalar_action_and_zero(leaf, mesh, rng):
    r"""2·d == d + d (scalar action, bit-exact) AND ψ ⊕ 0_disp == ψ."""
    psi1, psi2 = _make_flux(leaf, mesh, rng), _make_flux(leaf, mesh, rng)
    d = psi2 - psi1
    np.testing.assert_array_equal((2.0 * d).values, (d + d).values)
    zero = psi1 - psi1  # the zero displacement
    recovered = psi1 + zero
    if type(recovered) is not type(psi1):
        raise AssertionError(f"{leaf}: ψ ⊕ 0_disp did not return a flux")
    np.testing.assert_array_equal(recovered.values, psi1.values)


# ── 2. The gate — flux+flux forbidden, legal ops NOT forbidden (L11) ─────


@pytest.mark.parametrize("leaf", _LEAVES)
def test_flux_add_flux_forbidden_but_torsor_allowed(leaf, mesh, rng):
    r"""L11 paired: ψ ⊕ ψ → TypeError (no origin); ψ ⊕ displacement → flux."""
    psi1, psi2 = _make_flux(leaf, mesh, rng), _make_flux(leaf, mesh, rng)
    with pytest.raises(TypeError) as ei:
        _ = psi1 + psi2
    if "affine_combination" not in str(ei.value):
        raise AssertionError(
            f"{leaf}: flux+flux TypeError must name the legal alternative "
            f"(affine_combination); got: {str(ei.value)!r}"
        )
    out = psi1 + (psi2 - psi1)  # POSITIVE — torsor MUST NOT raise
    if type(out) is not type(psi1):
        raise AssertionError(f"{leaf}: torsor action did not return a flux")


@pytest.mark.parametrize("leaf", _LEAVES)
def test_displacement_add_displacement_returns_displacement(leaf, mesh, rng):
    r"""disp ⊕ disp → disp (V is closed under +); no raise."""
    psi1, psi2, psi3 = (_make_flux(leaf, mesh, rng) for _ in range(3))
    s = (psi2 - psi1) + (psi3 - psi2)
    if type(s) is not _displacement_type(psi1):
        raise AssertionError(f"{leaf}: disp+disp did not return a displacement")


@pytest.mark.parametrize("leaf", _LEAVES)
def test_affine_combination_partition_of_unity(leaf, mesh, rng):
    r"""L11 paired: Σλ=1 mints a flux equal to the weighted average AND the
    torsor relaxation form; Σλ≠1 RAISES."""
    psi1, psi2 = _make_flux(leaf, mesh, rng), _make_flux(leaf, mesh, rng)
    omega = 0.7  # 0.7 + 0.3 == 1.0 EXACTLY in float64
    combo = type(psi1).affine_combination([(omega, psi2), (1.0 - omega, psi1)])
    if type(combo) is not type(psi1):
        raise AssertionError(f"{leaf}: affine_combination(Σλ=1) did not return a flux")
    expected = omega * psi2.values + (1.0 - omega) * psi1.values
    np.testing.assert_array_equal(combo.values, expected)
    # The torsor relaxation form is the SAME blend in a different evaluation
    # order (ψ₁ + ω·Δψ vs (1−ω)ψ₁ + ω·ψ₂) — equal up to FP reorder. Relative
    # tolerance (NOT element-wise ULP, which inflates on the moment field's
    # near-zero ℓ≥1 / zero-padded entries — the Phase-5c small-magnitude
    # artifact). vv-principles: drift is (reduction depth)×ULP, FP-explainable.
    relax = psi1 + omega * (psi2 - psi1)
    np.testing.assert_allclose(relax.values, expected, rtol=1e-12, atol=1e-13)
    with pytest.raises(ValueError) as ei:  # Σλ = 1.4 ≠ 1
        _ = type(psi1).affine_combination([(0.7, psi2), (0.7, psi1)])
    if "1" not in str(ei.value):
        raise AssertionError(
            f"{leaf}: Σλ≠1 message should name the Σλ=1 constraint; "
            f"got {str(ei.value)!r}"
        )


def test_subtraction_typed_and_guards_intact(mesh, mesh2, rng):
    r"""__sub__ returns the typed displacement; cross-class and cross-mesh
    guards STILL raise (the carve must not loosen them)."""
    psi1 = _make_flux("angular", mesh, rng)
    psi2 = _make_flux("angular", mesh, rng)
    d = psi1 - psi2  # POSITIVE: works
    if type(d) is type(psi1):
        raise AssertionError("__sub__ returned a flux, not a displacement")
    from orpheus.transport.source_sinks.angular_source_sink import AngularSourceSink
    src = AngularSourceSink.zeros_on(mesh)
    with pytest.raises(TypeError):  # NEGATIVE cross-class
        _ = psi1 - src  # type: ignore[operator]
    other = AngularFlux.zeros_on(mesh2)
    with pytest.raises(ValueError):  # NEGATIVE cross-mesh
        _ = psi1 - other

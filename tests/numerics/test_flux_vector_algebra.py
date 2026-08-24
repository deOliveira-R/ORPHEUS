r"""Foundation tests for the flux VECTOR algebra — flux lives in V (CS3).

**The ruling this module pins** (campaign 1 CS3, 2026-08-19,
``.claude/plans/space_and_kernel_binding_campaign.md`` §0/§4): flux lives in
the positive cone :math:`K \subset V` of an ordered vector space. The additive
algebra of every flux leaf is therefore the plain vector-space algebra of
:class:`~orpheus.numerics.field.Field`:

* :math:`\psi_1 + \psi_2 \to \psi \in V` — LEGAL, the headline change,
* :math:`\psi_2 - \psi_1 \to` the SAME leaf type, SIGNED (differences leave K),
* the zero flux is the additive identity (V has an origin),
* the fiber discipline SURVIVES: cross-class and cross-mesh arithmetic still
  refuse, via the retained ``_check_partner`` chain (class identity + space
  equality + MESH identity, ``transport/fields/_bases.py``) — "fluxes of
  different problems don't mix" was always about the fiber, never about
  affine structure.

**History** (past tense, kept per ``coding-standards``): until 2026-08-19 this
module was ``test_affine_flux_algebra.py`` and pinned the #208/#201 affine
TORSOR algebra — ``flux ⊖ flux`` minted a ``Displacement``, ``flux + flux``
raised ``TypeError``, and ``affine_combination`` (Σλ=1) was the only multi-flux
blend. That ontology was OVERTURNED by user ruling ("we should not be bound by
past mistakes"); the relaxation blend it canonicalised is carried forward below
as plain algebra (`test_relaxation_blend_is_plain_algebra` — Q4 ruling: the
named-operation ceremony dissolves, the CONTENT stays pinned).

These are ``foundation`` tests — software invariants, not physics ``:label:``
claims. The structural ground is numpy ``+``/``-`` on the same buffers (no
solver). Per ``-O`` discipline every assertion is ``np.testing.assert_*`` or an
explicit ``raise`` (bare ``assert`` is stripped under ``-O``).

Leaf coverage note: this battery parameterizes the four historical leaves
(angular / scalar / moment / angular-boundary). The other three flux leaves
(``ScalarBoundaryFlux``, the two RadialCharacteristic leaves) need heavier
mesh fixtures and have the flip pinned in their own modules
(``tests/transport/fields/test_scalar_boundary_flux.py``,
``tests/transport/test_radial_characteristic_field.py``,
``tests/sn/mesh/test_radial_characteristic_split_leaves.py``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
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


def _stretched_mesh() -> SNMesh:
    """Doubled width, same shape — the cell VOLUMES differ, so the carrier
    mints an UNEQUAL space (the F2 content discriminator)."""
    m = Mesh1D(
        edges=np.linspace(0.0, 4.0, 5), mat_ids=np.zeros(4, dtype=int),
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


def _make_flux(leaf: str, m: SNMesh, rng: np.random.Generator) -> Field:
    """A flux leaf with DISTINCT random values (NOT flat — Mode-9)."""
    if leaf == "angular":
        return AngularFlux(values=rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape)), space=m.angular_bulk_space)
    if leaf == "scalar":
        return ScalarFlux(values=rng.standard_normal((m.ng, *m.spatial_shape)), space=m.bulk_space)
    if leaf == "moment":
        # HarmonicMomentFlux is rank-d like the other leaves: its space is
        # (L+1, 2L+1, ng, *spatial) — rank-1 (nx,) on a 1-D mesh, no phantom ny.
        shape = (_MOMENT_L + 1, 2 * _MOMENT_L + 1, m.ng, *m.spatial_shape)
        return HarmonicMomentFlux.from_mesh_and_L(
            rng.standard_normal(shape), m, _MOMENT_L,
        )
    if leaf == "boundary":
        zero = AngularBoundaryFlux.zeros(m.angular_trace)
        return AngularBoundaryFlux(values=rng.standard_normal(zero.values.shape), space=m.angular_trace)
    raise ValueError(leaf)


def _zeros_like_flux(leaf: str, m: SNMesh) -> Field:
    """The zero flux of the leaf — a legal, freely-constructed origin."""
    if leaf == "angular":
        return AngularFlux.zeros(m.angular_bulk_space)
    if leaf == "scalar":
        return ScalarFlux(values=np.zeros((m.ng, *m.spatial_shape)), space=m.bulk_space)
    if leaf == "moment":
        shape = (_MOMENT_L + 1, 2 * _MOMENT_L + 1, m.ng, *m.spatial_shape)
        return HarmonicMomentFlux.from_mesh_and_L(np.zeros(shape), m, _MOMENT_L)
    if leaf == "boundary":
        return AngularBoundaryFlux.zeros(m.angular_trace)
    raise ValueError(leaf)


# ── 1. flux + flux is LEGAL and arithmetically correct (the headline) ────


@pytest.mark.parametrize("leaf", _LEAVES)
def test_flux_add_flux_is_legal_and_exact(leaf, mesh, rng):
    r"""ψ₁ + ψ₂ returns the SAME leaf type carrying the bit-exact numpy sum.

    The whole ontology change in one row: until CS3 this raised ``TypeError``
    ("an affine space has no origin"). Commutativity is bit-exact (numpy `+`
    commutes elementwise); associativity only to FP rounding.

    Mutation (verified at the carve): ``Field.__add__`` body → subtraction
    reds the value legs here while the TYPE leg stays green — which is why
    the type leg alone is not the gate.
    """
    psi1, psi2, psi3 = (_make_flux(leaf, mesh, rng) for _ in range(3))
    s = psi1 + psi2
    if type(s) is not type(psi1):
        raise AssertionError(
            f"{leaf}: ψ₁+ψ₂ returned {type(s).__name__}, not the leaf type "
            f"{type(psi1).__name__} — the sum must stay in V's fiber."
        )
    np.testing.assert_array_equal(s.values, psi1.values + psi2.values)
    np.testing.assert_array_equal((psi2 + psi1).values, s.values)  # commutes
    # associativity to FP reorder — relative tolerance, not element-wise
    # ULP: near-cancelling sum entries inflate ULP distance (the same
    # small-magnitude artifact the relaxation row documents).
    np.testing.assert_allclose(
        ((psi1 + psi2) + psi3).values, (psi1 + (psi2 + psi3)).values,
        rtol=1e-12, atol=1e-13,
    )


@pytest.mark.parametrize("leaf", _LEAVES)
def test_zero_flux_is_the_additive_identity(leaf, mesh, rng):
    r"""ψ + 0 == ψ where 0 is a FLUX — the origin exists.

    The exact statement the retired "no natural zero" doctrine denied. The
    result must be a NEW object with equal values (``replace``-constructed),
    not ``ψ`` itself — a ``return self`` short-circuit would break the copy
    contract while passing a value-only assertion.
    """
    psi = _make_flux(leaf, mesh, rng)
    zero = _zeros_like_flux(leaf, mesh)
    out = psi + zero
    if out is psi:
        raise AssertionError(f"{leaf}: ψ + 0 returned ψ itself, not a copy")
    if type(out) is not type(psi):
        raise AssertionError(f"{leaf}: ψ + 0 changed the leaf type")
    np.testing.assert_array_equal(out.values, psi.values)


@pytest.mark.parametrize("leaf", _LEAVES)
def test_relaxation_blend_is_plain_algebra(leaf, mesh, rng):
    r"""0.7·ψ₂ + 0.3·ψ₁ — the relaxation blend, as ordinary V arithmetic.

    Carries forward the only real content of the retired
    ``affine_combination`` (Q4 ruling: the Σλ=1 ceremony dissolves; the blend
    is spellable everywhere). The torsor spelling ψ₁ + ω·(ψ₂ − ψ₁) is the
    SAME blend in a different evaluation order — equal to FP reorder
    (relative tolerance, NOT element-wise ULP, which inflates on the moment
    leaf's near-zero ℓ≥1 entries).
    """
    psi1, psi2 = _make_flux(leaf, mesh, rng), _make_flux(leaf, mesh, rng)
    omega = 0.7  # 0.7 + 0.3 == 1.0 exactly in float64
    blend = omega * psi2 + (1.0 - omega) * psi1
    if type(blend) is not type(psi1):
        raise AssertionError(f"{leaf}: the blend left the leaf type")
    expected = omega * psi2.values + (1.0 - omega) * psi1.values
    np.testing.assert_array_equal(blend.values, expected)
    relax = psi1 + omega * (psi2 - psi1)
    np.testing.assert_allclose(relax.values, expected, rtol=1e-12, atol=1e-13)


# ── 2. flux − flux is SAME-TYPED and SIGNED ──────────────────────────────


@pytest.mark.parametrize("leaf", _LEAVES)
def test_flux_difference_same_typed_and_signed(leaf, mesh, rng):
    r"""ψ₂ − ψ₁ returns the SAME leaf type, stores the raw signed difference.

    The inverted assertion — the retired battery RAISED if the difference
    kept the flux type. The signed leg is the cone's own boundary: a
    difference of cone members leaves K (V is signed, K ⊂ V strictly), which
    is exactly why cone membership must be an element PREDICATE and never a
    constructor invariant (§0 ruling 1).
    """
    psi1, psi2 = _make_flux(leaf, mesh, rng), _make_flux(leaf, mesh, rng)
    d = psi2 - psi1
    if type(d) is not type(psi1):
        raise AssertionError(
            f"{leaf}: ψ₂−ψ₁ returned {type(d).__name__}, not the leaf type — "
            f"the difference is a signed element of the SAME V."
        )
    np.testing.assert_array_equal(d.values, psi2.values - psi1.values)
    if not np.any(d.values < 0.0):
        raise AssertionError(
            f"{leaf}: the random difference has no negative entry — the "
            f"SIGNED claim was not exercised (degenerate fixture)."
        )


@pytest.mark.parametrize("leaf", _LEAVES)
def test_difference_round_trip_and_telescoping(leaf, mesh, rng):
    r"""ψ₁ + (ψ₂ − ψ₁) ≈ ψ₂ and (ψ₂−ψ₁) + (ψ₃−ψ₂) ≈ ψ₃−ψ₁ (to a few ULP).

    The update-step round-trip and path-independence, now statements about
    ONE type. ``a + (b − a) ≠ b`` bit-for-bit under IEEE-754 (the
    subtraction rounds), so both are asserted to ``nulp=8``, exactly as the
    retired torsor spelling was.
    """
    psi1, psi2, psi3 = (_make_flux(leaf, mesh, rng) for _ in range(3))
    np.testing.assert_array_almost_equal_nulp(
        (psi1 + (psi2 - psi1)).values, psi2.values, nulp=8,
    )
    np.testing.assert_array_almost_equal_nulp(
        ((psi2 - psi1) + (psi3 - psi2)).values, (psi3 - psi1).values, nulp=8,
    )


# ── 3. The fiber discipline SURVIVES (the charter's owed negative) ───────


def test_fiber_guard_cross_mesh_refuses(mesh, mesh2, rng):
    r"""The fiber discipline is space CONTENT (CS4b S3 re-derivation).

    This row was the CS3 negative control for the retained MESH arm of
    ``BulkField._check_partner``; its own docstring instructed re-derivation
    when identity semantics moved. CS4b S3 retired the provenance arm (the
    F2 ruling): the fiber IS the space now, so twin carriers — distinct
    instances, EQUAL content — legitimately mix (the correctly-blind leg),
    and a carrier whose cell volumes differ refuses on the base gate's
    space-content arm (the discriminator leg). Deleting the SPACE arm of
    ``Field._check_partner`` reds the refusal leg while the positive legs
    stay green (battery M6).
    """
    psi_a = _make_flux("angular", mesh, rng)
    psi_a2 = _make_flux("angular", mesh, rng)
    psi_b = _make_flux("angular", mesh2, rng)
    if mesh is mesh2:
        raise AssertionError("fixture defect: the two meshes are one object")
    if not (psi_a.space == psi_b.space):
        raise AssertionError(
            "twin carriers minted UNEQUAL spaces — the content identity has "
            "changed; re-derive this row against the new semantics."
        )
    out_twin = psi_a + psi_b  # twin content — legal since the F2 re-key
    if type(out_twin) is not type(psi_a):
        raise AssertionError("twin-carrier add did not return the leaf type")
    psi_c = _make_flux("angular", _stretched_mesh(), rng)
    with pytest.raises(ValueError, match="equal space"):  # NEGATIVE leg
        _ = psi_a + psi_c
    with pytest.raises(ValueError, match="equal space"):
        _ = psi_a - psi_c
    out = psi_a + psi_a2  # POSITIVE — same fiber adds (vv #11 pairing)
    if type(out) is not type(psi_a):
        raise AssertionError("same-mesh add did not return the leaf type")


def test_cross_class_still_refuses(mesh, rng):
    r"""AngularFlux ± AngularSourceSink → TypeError (Layer 1 survives).

    Same units never grant meaning: the state/rate-density distinction is a
    CLASS distinction and outlives the affine gate (it is what the
    ``test_typed_source_sinks`` contrast row re-poses onto).
    """
    psi = _make_flux("angular", mesh, rng)
    from orpheus.transport.source_sinks.angular_source_sink import AngularSourceSink

    src = AngularSourceSink.zeros(mesh.angular_bulk_space)
    with pytest.raises(TypeError):
        _ = psi + src  # type: ignore[operator]
    with pytest.raises(TypeError):
        _ = psi - src  # type: ignore[operator]


# ── 4. Scalar algebra unchanged (the regression floor) ───────────────────


@pytest.mark.parametrize("leaf", _LEAVES)
def test_scalar_algebra_unchanged(leaf, mesh, rng):
    r"""c·ψ, ψ/k, −ψ return the leaf type with exact values — the operations
    that were always legal (incl. the eigenvalue normalisation ψ/k that
    ``sn/solver.py``'s adjoint comment wrongly claimed the old doctrine
    forbade — #353's neighbourhood)."""
    psi = _make_flux(leaf, mesh, rng)
    for out, expected in [
        (2.5 * psi, 2.5 * psi.values),
        (psi / 1.25, psi.values / 1.25),
        (-psi, -psi.values),
    ]:
        if type(out) is not type(psi):
            raise AssertionError(f"{leaf}: scalar algebra changed the type")
        np.testing.assert_array_equal(out.values, expected)

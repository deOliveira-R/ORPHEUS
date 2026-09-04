r"""Gates for :mod:`orpheus.transport.material_field` — the per-material
kernel fields (CS4c step 3a, the O-6 landing).

Three tiers, deliberately separated:

* **Admission** — the pairing's construction contract, positive AND
  negative legs (`vv` #11).
* **Independent-reference gates** (permanent) — every verb against a
  HAND-ROLLED per-cell Python loop (plain ``@`` matmuls, no einsum, no
  shared dispatch code): the arms' arithmetic pinned by a structurally
  independent realization (`vv` L11).  The (n,2n) multiplicity is spelled
  as a hand-authored literal ``2.0`` here — the external written-down pin
  that keeps ``N2NKernel.multiplicity`` anchored after XD-2 removes every
  production literal (`coding-standards`, single-sourcing clause).
* The transitional facade battery (each verb bit-identical to the
  ``MaterialXSField.apply_*`` arm it replaced) RETIRED WITH THE ARMS at
  step 3c, as scheduled at its birth; the permanent coverage is the
  independent-reference tier.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from tests.sn._test_helpers import material_xs_from_raw
from orpheus.transport.kernels import N2NKernel, ScatteringKernel
from orpheus.transport.material_field import (
    MaterialField,
    N2NMaterialField,
    ScatteringMaterialField,
)

pytestmark = pytest.mark.foundation


# ── the asymmetric two-material fixture (L27: group-flip detectable) ────

_SIGS_A = [
    np.array([[0.38, 0.10], [0.05, 0.90]]),   # ℓ=0
    np.array([[0.20, -0.04], [0.02, 0.30]]),  # ℓ=1
    np.array([[0.06, 0.01], [-0.01, 0.09]]),  # ℓ=2
]
_SIGS_B = [
    np.array([[0.55, 0.03], [0.12, 0.40]]),
    np.array([[0.15, 0.05], [-0.03, 0.22]]),
    np.array([[0.04, -0.02], [0.02, 0.05]]),
]
_SIG2_A = np.array([[0.00, 0.03], [0.01, 0.00]])
_SIG2_B = np.array([[0.02, 0.00], [0.00, 0.04]])
_NX, _NG, _L = 6, 2, 2


def _mat_xs(nx=_NX):
    half = nx // 2
    ix, iy = np.arange(nx), np.zeros(nx, dtype=int)
    return material_xs_from_raw(
        sig_s={0: _SIGS_A, 1: _SIGS_B},
        sig2={0: _SIG2_A, 1: _SIG2_B},
        cells_by_mat={0: (ix[:half], iy[:half]), 1: (ix[half:], iy[half:])},
        ng=_NG, nx=nx, ny=1,
    )


def _fields(mat_xs=None):
    mat_xs = mat_xs if mat_xs is not None else _mat_xs()
    return (
        ScatteringMaterialField.from_material_xs(mat_xs),
        N2NMaterialField.from_material_xs(mat_xs),
        mat_xs,
    )


def _phi(sm=0, seed=1, nx=_NX):
    rng = np.random.default_rng(seed)
    shape = (_NG, nx, 1) if sm == 0 else (_NG, nx, 1, sm)
    return rng.uniform(0.1, 1.0, size=shape)


def _moments(sm=0, seed=3, nx=_NX):
    rng = np.random.default_rng(seed)
    shape = (
        (_L + 1, 2 * _L + 1, _NG, nx, 1)
        if sm == 0
        else (_L + 1, 2 * _L + 1, _NG, nx, 1, sm)
    )
    return rng.uniform(-1.0, 1.0, size=shape)


def _head(rank: int = 2):
    r"""The angular HEAD to hand a moment verb (#429 tracker 2.5).

    A moments tensor's leading axes ARE the head's, and which head a
    consumer holds is read off its frame — the rectangular ``(L+1, 2L+1)``
    real harmonics on a sphere rule, the FLAT ``(L+1,)`` Legendre family on a
    1-D rule's :math:`S^2/O(2)_x`. The verbs take it as a parameter because
    the per-degree contraction spec DEPENDS on it: `[M]` on a flat head the
    old einsum ``mfc...,fg->mgc...`` would have contracted the GROUP axis as
    if it were :math:`m` — silently, with no shape error.

    Both ranks are exercised in this module (``rank=2`` on the default
    fixtures, ``rank=1`` in :class:`TestFlatAngularHead`), so neither
    layout is certified by the other.
    """
    from orpheus.numerics.spaces.legendre_space import LegendreSpace
    from orpheus.numerics.spaces.spherical_harmonic_space import (
        SphericalHarmonicSpace,
    )

    return (
        SphericalHarmonicSpace.from_L(_L) if rank == 2
        else LegendreSpace.from_L(_L, "x")
    )


def _flat_moments(sm=0, seed=3, nx=_NX):
    """Moments in the FLAT (Legendre) head's layout — one coefficient per degree."""
    rng = np.random.default_rng(seed)
    shape = (_L + 1, _NG, nx, 1) if sm == 0 else (_L + 1, _NG, nx, 1, sm)
    return rng.uniform(-1.0, 1.0, size=shape)


def _sig_of(mid):
    return {0: _SIGS_A, 1: _SIGS_B}[mid], {0: _SIG2_A, 1: _SIG2_B}[mid]


def _mid_of_cell(ix):
    return 0 if ix < _NX // 2 else 1


# ═══════════════════════════════════════════════════════════════════════
# Admission — positive + negative legs (vv #11)
# ═══════════════════════════════════════════════════════════════════════


class TestAdmission:
    def test_a_correct_field_constructs(self):
        sf, nf, _ = _fields()
        if sf.ng != _NG or sf.order != _L or nf.ng != _NG:
            pytest.fail("correct field must construct with derived ng/order")

    def test_empty_map_refuses(self):
        with pytest.raises(ValueError, match="at least one material"):
            MaterialField(per_material={}, cells_by_material={})

    def test_layout_naming_kernel_less_material_refuses(self):
        k = ScatteringKernel(moments=(np.eye(_NG),))
        with pytest.raises(ValueError, match="carry no kernel"):
            MaterialField(
                per_material={0: k},
                cells_by_material={
                    0: (np.array([0]),), 7: (np.array([1]),),
                },
            )

    def test_mixed_ng_refuses(self):
        with pytest.raises(ValueError, match="uniform group structure"):
            MaterialField(
                per_material={
                    0: ScatteringKernel(moments=(np.eye(2),)),
                    1: ScatteringKernel(moments=(np.eye(3),)),
                },
                cells_by_material={0: (np.array([0]),)},
            )

    def test_mixed_order_refuses_on_order_read(self):
        f = MaterialField(  # base admits it — order is scattering-tier
            per_material={
                0: ScatteringKernel(moments=(np.eye(2),)),
                1: ScatteringKernel(moments=(np.eye(2), np.zeros((2, 2)))),
            },
            cells_by_material={0: (np.array([0]),)},
        )
        sf = ScatteringMaterialField(
            per_material=dict(f.per_material),
            cells_by_material=dict(f.cells_by_material),
        )
        with pytest.raises(ValueError, match="mixed truncation orders"):
            _ = sf.order

    def test_mappings_are_read_only(self):
        sf, _, _ = _fields()
        with pytest.raises(TypeError):
            sf.per_material[0] = None  # type: ignore[index]
        with pytest.raises(TypeError):
            sf.cells_by_material[0] = ()  # type: ignore[index]

    def test_moment_verb_refuses_order_mismatch(self):
        sf, _, _ = _fields()
        bad = np.zeros((_L + 2, 2 * (_L + 1) + 1, _NG, _NX, 1))
        with pytest.raises(ValueError, match="truncate the field"):
            sf.moment_source(bad, skip_l0=True, head=_head())


# ═══════════════════════════════════════════════════════════════════════
# Truncation + layout sharing
# ═══════════════════════════════════════════════════════════════════════


class TestTruncationAndSharing:
    def test_truncated_is_the_sub_stack_and_shares_the_layout(self):
        sf, _, _ = _fields()
        p0 = sf.truncated(0)
        if p0.order != 0:
            pytest.fail("truncated(0) must carry order 0")
        for mid in sf.per_material:
            np.testing.assert_array_equal(
                p0.per_material[mid].p0, sf.per_material[mid].p0,
            )
            if sf.cells_by_material[mid] is not p0.cells_by_material[mid]:
                pytest.fail(
                    "truncation must SHARE the layout index arrays, not copy"
                )

    def test_truncated_above_stored_order_refuses(self):
        sf, _, _ = _fields()
        with pytest.raises(ValueError, match="not invented"):
            sf.truncated(_L + 1)

    def test_two_fields_share_the_mesh_layout_object(self):
        mat_xs = _mat_xs()
        sf = ScatteringMaterialField.from_material_xs(mat_xs)
        nf = N2NMaterialField.from_material_xs(mat_xs)
        for mid in sf.cells_by_material:
            if sf.cells_by_material[mid] is not nf.cells_by_material[mid]:
                pytest.fail(
                    "both channel fields must share the mesh's ONE cached "
                    "where-partition (no per-field recomputation)"
                )

    def test_mesh_layout_is_cached(self):
        mat_xs = _mat_xs()
        if mat_xs.mesh.cells_by_material is not mat_xs.mesh.cells_by_material:
            pytest.fail("MaterialMesh.cells_by_material must be cached")


# ═══════════════════════════════════════════════════════════════════════
# Independent-reference gates (permanent; hand-rolled per-cell loops)
# ═══════════════════════════════════════════════════════════════════════


class TestIndependentReference:
    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    def test_p0_source(self, sm):
        sf, _, _ = _fields()
        phi = _phi(sm)
        Q = np.zeros_like(phi)
        sf.add_p0_source(Q, phi)
        ref = np.zeros_like(phi)
        for ix in range(_NX):
            s0 = _sig_of(_mid_of_cell(ix))[0][0]
            ref[:, ix, 0] += np.tensordot(s0.T, phi[:, ix, 0], axes=1)
        np.testing.assert_allclose(Q, ref, rtol=0, atol=1e-15)

    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    def test_p0_source_transpose(self, sm):
        sf, _, _ = _fields()
        chi = _phi(sm, seed=2)
        Q = np.zeros_like(chi)
        sf.add_p0_source_transpose(Q, chi)
        ref = np.zeros_like(chi)
        for ix in range(_NX):
            s0 = _sig_of(_mid_of_cell(ix))[0][0]
            ref[:, ix, 0] += np.tensordot(s0, chi[:, ix, 0], axes=1)
        np.testing.assert_allclose(Q, ref, rtol=0, atol=1e-15)

    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    @pytest.mark.parametrize("skip_l0", [True, False], ids=["l>=1", "l>=0"])
    def test_moment_source(self, sm, skip_l0):
        sf, _, _ = _fields()
        mom = _moments(sm)
        out = sf.moment_source(mom, skip_l0=skip_l0, head=_head())
        ref = np.zeros_like(mom)
        for ix in range(_NX):
            sig_s = _sig_of(_mid_of_cell(ix))[0]
            for l in range(1 if skip_l0 else 0, _L + 1):
                for m in range(2 * l + 1):
                    ref[l, m, :, ix, 0] = np.tensordot(
                        sig_s[l].T, mom[l, m, :, ix, 0], axes=1,
                    )
        np.testing.assert_allclose(out, ref, rtol=0, atol=1e-15)

    @pytest.mark.parametrize("skip_l0", [True, False], ids=["l>=1", "l>=0"])
    def test_moment_source_transpose(self, skip_l0):
        sf, _, _ = _fields()
        mom = _moments(seed=5)
        out = sf.moment_source_transpose(mom, skip_l0=skip_l0, head=_head())
        ref = np.zeros_like(mom)
        for ix in range(_NX):
            sig_s = _sig_of(_mid_of_cell(ix))[0]
            for l in range(1 if skip_l0 else 0, _L + 1):
                for m in range(2 * l + 1):
                    ref[l, m, :, ix, 0] = np.tensordot(
                        sig_s[l], mom[l, m, :, ix, 0], axes=1,
                    )
        np.testing.assert_allclose(out, ref, rtol=0, atol=1e-15)

    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    def test_n2n_emission(self, sm):
        _, nf, _ = _fields()
        phi = _phi(sm)
        Q = np.zeros_like(phi)
        nf.add_emission(Q, phi)
        ref = np.zeros_like(phi)
        for ix in range(_NX):
            sig2 = _sig_of(_mid_of_cell(ix))[1]
            # 2.0 = the hand-authored multiplicity pin (see module docstring)
            ref[:, ix, 0] += 2.0 * np.tensordot(sig2.T, phi[:, ix, 0], axes=1)
        np.testing.assert_allclose(Q, ref, rtol=0, atol=1e-15)

    def test_n2n_emission_transpose(self):
        _, nf, _ = _fields()
        chi = _phi(seed=7)
        Q = np.zeros_like(chi)
        nf.add_emission_transpose(Q, chi)
        ref = np.zeros_like(chi)
        for ix in range(_NX):
            sig2 = _sig_of(_mid_of_cell(ix))[1]
            ref[:, ix, 0] += 2.0 * np.tensordot(sig2, chi[:, ix, 0], axes=1)
        np.testing.assert_allclose(Q, ref, rtol=0, atol=1e-15)

    def test_n2n_moment_emission_touches_only_l0(self):
        _, nf, _ = _fields()
        mom = _moments(seed=9)
        out = nf.moment_emission(mom, head=_head())
        if not np.array_equal(out[1:], np.zeros_like(out[1:])):
            pytest.fail(
                "the (n,2n) field carries ONE matrix (ORPHEUS's P0 model of the "
                "channel — NOT a property of the reaction, see #426), so every "
                "ℓ≥1 block must be zero"
            )
        if not np.array_equal(out[0, 1:], np.zeros_like(out[0, 1:])):
            pytest.fail("(n,2n) writes only the (ℓ=0, m=0) slot")
        ref = np.zeros_like(mom[0, 0])
        for ix in range(_NX):
            sig2 = _sig_of(_mid_of_cell(ix))[1]
            ref[:, ix, 0] = 2.0 * np.tensordot(sig2.T, mom[0, 0, :, ix, 0], axes=1)
        np.testing.assert_allclose(out[0, 0], ref, rtol=0, atol=1e-15)

    def test_n2n_moment_emission_transpose(self):
        _, nf, _ = _fields()
        mom = _moments(seed=11)
        out = nf.moment_emission_transpose(mom, head=_head())
        ref = np.zeros_like(mom[0, 0])
        for ix in range(_NX):
            sig2 = _sig_of(_mid_of_cell(ix))[1]
            ref[:, ix, 0] = 2.0 * np.tensordot(sig2, mom[0, 0, :, ix, 0], axes=1)
        np.testing.assert_allclose(out[0, 0], ref, rtol=0, atol=1e-15)

    def test_group_rate(self):
        _, nf, _ = _fields()
        phi = _phi()
        vol = np.linspace(0.5, 2.0, _NX).reshape(_NX, 1)
        rate = np.zeros(_NG)
        nf.add_to_group_rate(rate, phi, vol)
        ref = np.zeros(_NG)
        for ix in range(_NX):
            sig2 = _sig_of(_mid_of_cell(ix))[1]
            for g in range(_NG):
                for gp in range(_NG):
                    ref[g] += vol[ix, 0] * 2.0 * phi[gp, ix, 0] * sig2[gp, g]
        np.testing.assert_allclose(rate, ref, rtol=1e-13)

    def test_reciprocity_forward_vs_transpose(self):
        r"""``⟨Λφ, χ⟩ = ⟨φ, Λᵀχ⟩`` over the full moment contraction — the
        forward/transpose pair cross-checks itself with no reference."""
        sf, _, _ = _fields()
        a, b = _moments(seed=13), _moments(seed=14)
        lhs = float(np.sum(sf.moment_source(a, skip_l0=False, head=_head()) * b))
        rhs = float(
            np.sum(a * sf.moment_source_transpose(b, skip_l0=False, head=_head()))
        )
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12)

# (The transitional facade battery — every verb bit-identical to the
# MaterialXSField arm it replaced — RETIRED WITH THE ARMS at step 3c,
# as its own docstring scheduled. Permanent coverage: the hand-rolled
# independent references above.)


# ═══════════════════════════════════════════════════════════════════════
# FissionMaterialField (CS4c step 4a) — validated pairs + gather verbs
# ═══════════════════════════════════════════════════════════════════════

# Two DISTINCT producing materials (asymmetric χ AND νΣf so a factor swap
# or a material-id mixup is detectable), one non-producing (the null-χ
# branch of the simplex law rides through the field path).
_CHI = {0: np.array([0.9, 0.1]), 1: np.array([0.7, 0.3]),
        2: np.array([0.0, 0.0])}
_NU_SIG_F = {0: np.array([0.013, 0.26]), 1: np.array([0.002, 0.11]),
             2: np.array([0.0, 0.0])}


def _fission_field(nx=_NX):
    from orpheus.transport.kernels import FissionKernel
    from orpheus.transport.material_field import FissionMaterialField

    third = nx // 3
    ix, iy = np.arange(nx), np.zeros(nx, dtype=int)
    return FissionMaterialField(
        per_material={
            mid: FissionKernel(chi=_CHI[mid], nu_sig_f=_NU_SIG_F[mid])
            for mid in _CHI
        },
        cells_by_material={
            0: (ix[:third], iy[:third]),
            1: (ix[third:2 * third], iy[third:2 * third]),
            2: (ix[2 * third:], iy[2 * third:]),
        },
    )


class TestFissionAdmission:
    def test_a_correct_field_constructs(self):
        ff = _fission_field()
        if ff.ng != _NG:
            pytest.fail("correct fission field must construct with ng=2")

    def test_simplex_law_runs_per_material_through_the_field_path(self):
        """A producing material with a non-simplex χ cannot enter the
        field — the law fires in :class:`FissionKernel`'s ctor, which is
        the ONLY door (Pattern 4: the field holds validated kernels, so
        an invalid spectrum is not a value a field can carry)."""
        from orpheus.transport.kernels import FissionKernel
        from orpheus.transport.material_field import FissionMaterialField

        with pytest.raises(ValueError):
            FissionMaterialField(
                per_material={
                    0: FissionKernel(
                        chi=np.array([0.5, 0.1]),         # not a simplex
                        nu_sig_f=np.array([0.1, 0.2]),    # producing
                    ),
                },
                cells_by_material={0: (np.arange(_NX), np.zeros(_NX, int))},
            )

    def test_gather_arity_mismatch_refuses(self):
        with pytest.raises(ValueError, match="spatial axes"):
            _fission_field().gather_chi((_NX,))

    def test_gather_out_of_bounds_fails_loudly(self):
        """A shape smaller than the layout's reach raises (numpy's
        IndexError) — a wrong-space gather can never silently truncate."""
        with pytest.raises(IndexError):
            _fission_field().gather_chi((2, 1))

    def test_gathered_factors_are_write_protected(self):
        chi = _fission_field().gather_chi((_NX, 1))
        with pytest.raises(ValueError):
            chi[0, 0, 0] = 5.0


class TestFissionIndependentReference:
    """Gathers against a hand-rolled per-cell loop (`vv` L11: plain
    indexing, no shared dispatch code)."""

    def test_gather_chi_and_nu_sig_f(self):
        ff = _fission_field()
        chi = ff.gather_chi((_NX, 1))
        nu = ff.gather_nu_sig_f((_NX, 1))
        third = _NX // 3
        for ix in range(_NX):
            mid = 0 if ix < third else (1 if ix < 2 * third else 2)
            np.testing.assert_array_equal(chi[:, ix, 0], _CHI[mid])
            np.testing.assert_array_equal(nu[:, ix, 0], _NU_SIG_F[mid])

    def test_gather_is_bit_identical_to_the_facade_views(self):
        """The gather ≡ the facade's dense per-cell views (both are pure
        index gathers of the same per-material vectors — this row pins
        that no arithmetic crept into either route). Lives until F-1
        retires the facade views; then the hand-rolled row above is the
        surviving reference."""
        from scipy.sparse import csr_matrix

        from orpheus.data.macro_xs.mixture import Mixture
        from orpheus.geometry import Mesh2D
        from orpheus.transport.material_field import FissionMaterialField
        from orpheus.transport.mesh.material_mesh import MaterialMesh
        from orpheus.transport.mesh.material_xs_field import MaterialXSField

        z = np.zeros(_NG)
        materials = {
            mid: Mixture(
                SigC=z.copy(), SigL=z.copy(),
                SigF=_NU_SIG_F[mid] / 2.4, SigP=_NU_SIG_F[mid].copy(),
                SigT=np.ones(_NG),
                SigS=[csr_matrix(np.zeros((_NG, _NG)))],
                Sig2=[csr_matrix(np.zeros((_NG, _NG)))],
                chi=_CHI[mid].copy(),
            )
            for mid in _CHI
        }
        third = _NX // 3
        mat_map = np.zeros((_NX, 1), dtype=int)
        mat_map[third:2 * third, :] = 1
        mat_map[2 * third:, :] = 2
        mesh = Mesh2D(
            edges_x=np.arange(_NX + 1, dtype=float),
            edges_y=np.arange(2, dtype=float),
            mat_map=mat_map,
        )
        mat_xs = MaterialXSField.from_mesh(MaterialMesh(mesh, materials))
        ff = FissionMaterialField.from_material_xs(mat_xs)
        np.testing.assert_array_equal(
            ff.gather_chi((_NX, 1)), mat_xs.emission_spectrum,
        )
        np.testing.assert_array_equal(
            ff.gather_nu_sig_f((_NX, 1)), mat_xs.fission_production,
        )


# ═══════════════════════════════════════════════════════════════════════
# The FLAT angular head (#429 tracker 2.5) — the verbs' new dispatch axis
# ═══════════════════════════════════════════════════════════════════════


class TestFlatAngularHead:
    r"""A rank-1 head is DISPATCHED, and the failure it replaces was SILENT.

    ⛔ **§6c.** ``head=`` is a new parameter whose whole purpose is to select a
    per-degree contraction. Every row above hands the RECTANGULAR head, so
    without this class the rank-1 arm would land with no input that reaches
    it — a gate structurally unable to fail. The rows here construct that
    input and pin the arm against a hand-rolled reference.

    ⭐ **Why the old spelling was silent rather than loud.** The contraction
    was ``mfc...,fg->mgc...``: on a flat ``(L+1, ng, nx, …)`` tensor the
    einsum's ``m`` binds the DEGREE axis and ``f`` binds the GROUP axis, so it
    contracts the group index against the degree index and returns a
    well-shaped array of the wrong quantity — the same failure family as
    ``values[0, 0]`` on a flat moment field.
    """

    @pytest.mark.parametrize("sm", [0, 4], ids=["scalar", "LD-2^d=4"])
    @pytest.mark.parametrize("skip_l0", [True, False], ids=["l>=1", "l>=0"])
    def test_moment_source_on_a_flat_head(self, sm, skip_l0):
        """The independent reference, one degree per coefficient."""
        sf, _, _ = _fields()
        mom = _flat_moments(sm)
        out = sf.moment_source(mom, skip_l0=skip_l0, head=_head(rank=1))
        ref = np.zeros_like(mom)
        for ix in range(_NX):
            sig_s = _sig_of(_mid_of_cell(ix))[0]
            for l in range(1 if skip_l0 else 0, _L + 1):
                ref[l, :, ix, 0] = np.tensordot(
                    sig_s[l].T, mom[l, :, ix, 0], axes=1,
                )
        np.testing.assert_allclose(out, ref, rtol=0, atol=1e-15)

    def test_moment_source_transpose_on_a_flat_head(self):
        sf, _, _ = _fields()
        mom = _flat_moments(seed=5)
        out = sf.moment_source_transpose(mom, skip_l0=False, head=_head(rank=1))
        ref = np.zeros_like(mom)
        for ix in range(_NX):
            sig_s = _sig_of(_mid_of_cell(ix))[0]
            for l in range(_L + 1):
                ref[l, :, ix, 0] = np.tensordot(
                    sig_s[l], mom[l, :, ix, 0], axes=1,
                )
        np.testing.assert_allclose(out, ref, rtol=0, atol=1e-15)

    def test_the_rectangular_spec_would_be_SILENT_on_a_flat_tensor(self):
        r"""The demonstration: the old per-degree loop RETURNS on a flat tensor.

        Not a claim about production — production now dispatches on the head —
        but the record of why no test could have caught this by running. The
        retired loop sliced ``moments[l, :2l+1]`` and contracted
        ``"mfc...,fg->mgc..."``. On a flat ``(L+1, ng, nx, 1)`` tensor axis 1
        is the GROUP axis, so ``[:2l+1]`` selects groups and the einsum's
        ``m`` binds what is left of them: a well-shaped array of the wrong
        quantity, no exception anywhere.
        """
        sf, _, _ = _fields()
        mom = _flat_moments(seed=11)
        honest = sf.moment_source(mom, skip_l0=False, head=_head(rank=1))

        # the RETIRED loop, transcribed and applied to the flat tensor
        wrong = np.zeros_like(mom)
        for kernel_idx, ix in enumerate(range(_NX)):
            sig_s = _sig_of(_mid_of_cell(ix))[0]
            for l in range(_L + 1):
                view = mom[l, : 2 * l + 1][..., ix, 0]      # selects GROUPS
                if view.shape[0] == 0:
                    continue
                contracted = np.einsum(
                    "mf,fg->mg", view.reshape(view.shape[0], -1)[:, :_NG],
                    sig_s[l],
                )
                wrong[l, : contracted.shape[0], ix, 0] = contracted[:, 0]

        assert wrong.shape == honest.shape, (
            "the retired spelling is SHAPE-compatible on a flat tensor — that "
            "is the finding, and it is why the defect was silent"
        )
        assert not np.allclose(wrong, honest), (
            "the retired spelling must disagree with the head-dispatched one, "
            "else this arm carries no claim"
        )

    def test_a_head_of_an_unshipped_rank_is_refused_naming_the_shipped_ones(self):
        """The dispatch names its own domain rather than falling through."""

        @dataclass(frozen=True)
        class _RankThreeHead:
            """A minimal ``MomentHead`` of an unshipped rank — a stub, not a space."""

            L: int
            shape: tuple[int, ...]
            name: str = "rank3_head"
            isotropic_slot: tuple[int, ...] = (0, 0, 0)

            def degree_block(self, l, /):
                return (l, slice(None), slice(None))

            def truncated(self, L_new, /):
                raise NotImplementedError

        sf, _, _ = _fields()
        head = _RankThreeHead(L=_L, shape=(_L + 1, 2 * _L + 1, 1))
        mom = np.zeros((_L + 1, 2 * _L + 1, 1, _NG, _NX, 1))
        with pytest.raises(NotImplementedError) as excinfo:
            sf.moment_source(mom, skip_l0=False, head=head)
        message = str(excinfo.value)
        assert "rank 3" in message
        assert "rank 2" in message and "rank 1" in message, (
            f"the refusal must name the shipped ranks: {message}"
        )

    def test_a_moments_tensor_whose_leading_axes_are_not_the_head_is_refused(self):
        r"""The head is a CONTRACT on the operand, not a hint.

        ⚠ **The guard is one-sided, and the limitation is recorded rather
        than asserted away.** It compares ``moments.shape[:rank]`` against the
        head's shape, so it catches a FLAT tensor offered to a RECTANGULAR
        head (``(3, 2)`` vs ``(3, 5)`` — the direction below) and does NOT
        catch the reverse when the leading extent happens to agree: a
        rectangular ``(3, 5, ng, …)`` tensor with a flat ``(3,)`` head passes
        it and dies further in with ``IndexError: index 2 is out of bounds for
        axis 1``. Unreachable from production (every call site hands the
        operator's OWN head, minted from the same frame as the moments), so it
        is a message-quality gap, not a correctness one — see the return.
        """
        sf, _, _ = _fields()
        with pytest.raises(ValueError, match="lead with the angular head"):
            sf.moment_source(_flat_moments(), skip_l0=False, head=_head())

        # the reverse direction — a RECTANGULAR tensor offered with a FLAT
        # head — is refused at the SAME guard since the group axis is checked
        # too (until 2026-09-02 it passed the guard and died in the loop with
        # an IndexError, reading its m-axis as the groups)
        with pytest.raises(ValueError, match="lead with the angular head"):
            sf.moment_source(_moments(), skip_l0=False, head=_head(rank=1))

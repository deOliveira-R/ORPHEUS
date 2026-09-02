r"""G-C1 rows for the S-family tier-2 classmethods (CS4c step 3; the
three-tier discipline §3 — EVERY extract-and-mint classmethod owes
``classmethod-built ≡ ctor-built`` on the same inputs, closing vv#28's
simple-ctor-vs-composite-factory blindness for this family).

Rows: ``ScatteringOperator.from_solver_data``,
``LegendreMomentScattering.from_material_xs``,
``N2NMomentOperator.from_material_xs``, the iso pair's
``from_material_xs``, ``N2NOperator.from_solver_data``, and one
parametrized row over the three kernel ``from_mixture``s (the datum
tier — proportionate to its one-``todense()`` exposure; the existing
``test_kernels.py`` value pins are cross-source gates, not this shape).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import SphericalHarmonicSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.frames.harmonic_frame import HarmonicFrame
from orpheus.transport.kernels import FissionKernel, N2NKernel, ScatteringKernel
from orpheus.transport.material_field import (
    N2NMaterialField,
    ScatteringMaterialField,
)
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.n2n import N2NOperator
from orpheus.transport.operators.scattering import (
    LegendreMomentScattering,
    N2NMomentOperator,
    ScatteringOperator,
)

from tests.sn._test_helpers import material_xs_from_raw
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis

pytestmark = pytest.mark.foundation

_SIGS = [
    np.array([[0.38, 0.10], [0.05, 0.90]]),
    np.array([[0.20, -0.04], [0.02, 0.30]]),
]
_SIG2 = np.array([[0.00, 0.03], [0.01, 0.00]])
_NX = 4


def _mat_xs():
    return material_xs_from_raw(
        sig_s={0: _SIGS}, sig2={0: _SIG2},
        cells_by_mat={0: (np.arange(_NX), np.zeros(_NX, dtype=int))},
        ng=2, nx=_NX, ny=1,
    )


def _sn(mat_xs):
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(n_ordinates=4), mat_xs.materials)


def _fields_equal(a, b):
    if set(a.per_material) != set(b.per_material):
        return False
    for mid in a.per_material:
        ka, kb = a.per_material[mid], b.per_material[mid]
        def _arrs(k):
            # per-channel datum arrays: Legendre stack (S), reaction
            # matrix (N2N), or the factor pair (Fission).
            if hasattr(k, "moments"):
                return tuple(k.moments)
            if hasattr(k, "matrix"):
                return (k.matrix,)
            return (k.chi, k.nu_sig_f)

        arrs_a = _arrs(ka)
        arrs_b = _arrs(kb)
        if len(arrs_a) != len(arrs_b):
            return False
        for x, y in zip(arrs_a, arrs_b):
            if not np.array_equal(x, y):
                return False
    return True


class TestSFamilyTierTwoEquivalence:
    def test_scattering_from_solver_data_equals_the_exact_ctor(self):
        mat_xs = _mat_xs()
        sn = _sn(mat_xs)
        space = sn.full_field_space
        rich = ScatteringOperator.from_solver_data(
            mat_xs=mat_xs, scattering_order=1, space=space,
        )
        interior = space.interior_space
        frame = HarmonicFrame.for_space(interior, 1)
        exact = ScatteringOperator(
            ScatteringMaterialField.from_material_xs(mat_xs).truncated(1),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space,
            codomain=space,
        )
        if not (rich.domain is exact.domain and rich.codomain is exact.codomain):
            pytest.fail("ends must be the SAME space object")
        if rich.flux_analysis.frame is not exact.flux_analysis.frame:
            pytest.fail("both routes must land on the ONE interned frame")
        if not _fields_equal(rich.scattering, exact.scattering):
            pytest.fail("the extracted kernel field drifted from the ctor's")

    def test_moment_pair_from_material_xs_equals_the_exact_ctor(self):
        mat_xs = _mat_xs()
        sh = SphericalHarmonicSpace.from_L(1)
        rich = LegendreMomentScattering.from_material_xs(mat_xs, SphericalHarmonicBasis(L=1), skip_l0=False)
        exact = LegendreMomentScattering(
            ScatteringMaterialField.from_material_xs(mat_xs).truncated(1),
            skip_l0=False, domain=sh, codomain=sh,
        )
        if not (rich.domain == exact.domain and rich.codomain == exact.codomain):
            pytest.fail("Λ ends drifted")
        if rich.skip_l0 != exact.skip_l0 or not _fields_equal(
            rich.scattering, exact.scattering,
        ):
            pytest.fail("Λ datum drifted")

        rich_n = N2NMomentOperator.from_material_xs(mat_xs, SphericalHarmonicBasis(L=1))
        exact_n = N2NMomentOperator(
            N2NMaterialField.from_material_xs(mat_xs), domain=sh, codomain=sh,
        )
        if not (
            rich_n.domain == exact_n.domain
            and _fields_equal(rich_n.n2n, exact_n.n2n)
        ):
            pytest.fail("N2N moment datum/ends drifted")

    def test_iso_pair_from_material_xs_equals_the_exact_ctor(self):
        mat_xs = _mat_xs()
        space = mat_xs.mesh.bulk_space
        rich = IsotropicScattering.from_material_xs(mat_xs, space=space)
        exact = IsotropicScattering(
            ScatteringMaterialField.from_material_xs(mat_xs).truncated(0),
            domain=space, codomain=space,
        )
        if not (
            rich.domain is exact.domain
            and _fields_equal(rich.scattering, exact.scattering)
        ):
            pytest.fail("IsoS classmethod drifted from the ctor")
        rich_n = IsotropicN2N.from_material_xs(mat_xs, space=space)
        exact_n = IsotropicN2N(
            N2NMaterialField.from_material_xs(mat_xs),
            domain=space, codomain=space,
        )
        if not (
            rich_n.domain is exact_n.domain
            and _fields_equal(rich_n.n2n, exact_n.n2n)
        ):
            pytest.fail("IsoN2N classmethod drifted from the ctor")

    def test_n2n_operator_from_solver_data_equals_the_exact_ctor(self):
        from orpheus.numerics.space import FunctionSpace

        mat_xs = _mat_xs()
        sn = _sn(mat_xs)
        space = sn.full_field_space
        rich = N2NOperator.from_solver_data(mat_xs=mat_xs, space=space)
        interior = space.interior_space
        assert interior is not None and interior.axes is not None
        scalar = FunctionSpace.of_axes(*interior.axes[1:])
        # CS4c step-4 harmonization: the retained state is the L=0 frame
        # (hub-interned — the exact ctor reaches the SAME object, so the
        # equivalence is an identity, stronger than the old array_equal
        # on a weights copy).
        from orpheus.transport.frames.harmonic_frame import HarmonicFrame

        exact = N2NOperator(
            IsotropicN2N(
                N2NMaterialField.from_material_xs(mat_xs),
                domain=scalar, codomain=scalar,
            ),
            frame=HarmonicFrame.for_space(interior, 0),
            domain=space, codomain=space,
        )
        if not (rich.domain is exact.domain and rich.codomain is exact.codomain):
            pytest.fail("N2N ends drifted")
        if rich.frame is not exact.frame:
            pytest.fail("N2N frame drifted from the interned hub route")
        np.testing.assert_array_equal(
            np.asarray(rich.frame.measure.weights), np.asarray(sn.quad.weights),
        )
        if not (
            rich.energy.domain == exact.energy.domain
            and _fields_equal(rich.energy.n2n, exact.energy.n2n)
        ):
            pytest.fail("N2N energy binding drifted from the ctor route")

    @pytest.mark.parametrize("cls,fields", [
        (ScatteringKernel, lambda mix: {"moments": tuple(
            np.asarray(s.todense()) for s in mix.SigS)}),
        (N2NKernel, lambda mix: {"matrix": np.asarray(mix.Sig2.todense())}),
        (FissionKernel, lambda mix: {"chi": mix.chi, "nu_sig_f": mix.SigP}),
    ], ids=["scattering", "n2n", "fission"])
    def test_kernel_from_mixture_equals_the_exact_ctor(self, cls, fields):
        """The datum tier's single parametrized row (plan §4: the
        extraction is one ``todense()``; the existing value pins are
        cross-source gates, not this shape)."""
        mix = _mat_xs().materials[0]
        rich = cls.from_mixture(mix)
        exact = cls(**fields(mix))
        for name in fields(mix):
            a, b = getattr(rich, name), getattr(exact, name)
            if isinstance(a, tuple):
                for x, y in zip(a, b):
                    np.testing.assert_array_equal(x, y)
                np.testing.assert_equal(len(a), len(b))
            else:
                np.testing.assert_array_equal(a, b)


def test_fission_operator_from_solver_data_equals_the_exact_ctor():
    """G-C1 row for the ANGULAR fission binding (CS4c step 4): tier-2 ≡
    the exact ctor on the same inputs — energy binding, hub-interned L=0
    frame (an IDENTITY, the interning theorem), and both ends."""
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.frames.harmonic_frame import HarmonicFrame
    from orpheus.transport.material_field import FissionMaterialField
    from orpheus.transport.operators.fission import FissionOperator
    from orpheus.transport.operators.isotropic_scattering import (
        IsotropicFission,
    )

    mat_xs = _mat_xs()
    sn = _sn(mat_xs)
    space = sn.full_field_space
    rich = FissionOperator.from_solver_data(mat_xs=mat_xs, space=space)
    interior = space.interior_space
    assert interior is not None and interior.axes is not None
    scalar = FunctionSpace.of_axes(*interior.axes[1:])
    exact = FissionOperator(
        IsotropicFission(
            FissionMaterialField.from_material_xs(mat_xs),
            domain=scalar, codomain=scalar,
        ),
        frame=HarmonicFrame.for_space(interior, 0),
        domain=space, codomain=space,
    )
    if not (rich.domain is exact.domain and rich.codomain is exact.codomain):
        pytest.fail("F ends drifted")
    if rich.frame is not exact.frame:
        pytest.fail("F frame drifted from the interned hub route")
    if not _fields_equal(rich.energy.fission, exact.energy.fission):
        pytest.fail("F energy binding drifted from the ctor route")

r"""G-C1 rows for the S-family tier-2 classmethods (CS4c step 3; the
three-tier discipline §3 — EVERY extract-and-mint classmethod owes
``classmethod-built ≡ ctor-built`` on the same inputs, closing vv#28's
simple-ctor-vs-composite-factory blindness for this family).

Rows: ``ScatteringOperator.from_solver_data``,
``LegendreMomentTransfer.from_field`` (both channels), the iso pair's
``from_material_xs``, ``N2NOperator.from_solver_data`` (since #426 step 2
at the solve's order — the SAME interned frame as S), and one parametrized
row over the three kernel channel constructors (the datum tier —
proportionate to its one-``todense()`` exposure; the existing
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
from orpheus.transport.kernels import N2N_MULTIPLICITY, FissionKernel, TransferKernel
from orpheus.transport.material_field import TransferMaterialField
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.n2n import N2NOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.operators.transfer import LegendreMomentTransfer

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
            # per-channel datum arrays: the Legendre stack + yield (a
            # transfer channel) or the factor pair (fission).
            if hasattr(k, "moments"):
                return (*k.moments, np.asarray(k.multiplicity))
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
        assert interior is not None
        frame = HarmonicFrame.for_space(interior, 1)
        exact = ScatteringOperator(
            TransferMaterialField.scattering(mat_xs).at_order(1),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space,
            codomain=space,
        )
        if not (rich.domain is exact.domain and rich.codomain is exact.codomain):
            pytest.fail("ends must be the SAME space object")
        if rich.flux_analysis.frame is not exact.flux_analysis.frame:
            pytest.fail("both routes must land on the ONE interned frame")
        if not _fields_equal(rich.transfer, exact.transfer):
            pytest.fail("the extracted kernel field drifted from the ctor's")

    def test_moment_pair_from_material_xs_equals_the_exact_ctor(self):
        mat_xs = _mat_xs()
        sh = SphericalHarmonicSpace.from_L(1)
        rich = LegendreMomentTransfer.on_basis(
            TransferMaterialField.scattering(mat_xs), SphericalHarmonicBasis(L=1), skip_l0=False,
        )
        exact = LegendreMomentTransfer(
            TransferMaterialField.scattering(mat_xs).at_order(1),
            skip_l0=False, domain=sh, codomain=sh,
        )
        if not (rich.domain == exact.domain and rich.codomain == exact.codomain):
            pytest.fail("Λ ends drifted")
        if rich.skip_l0 != exact.skip_l0 or not _fields_equal(
            rich.transfer, exact.transfer,
        ):
            pytest.fail("Λ datum drifted")

        rich_n = LegendreMomentTransfer.on_basis(
            TransferMaterialField.n2n(mat_xs), SphericalHarmonicBasis(L=1), skip_l0=False,
        )
        exact_n = LegendreMomentTransfer(
            TransferMaterialField.n2n(mat_xs).at_order(1),
            skip_l0=False, domain=sh, codomain=sh,
        )
        if not (
            rich_n.domain == exact_n.domain
            and rich_n.skip_l0 == exact_n.skip_l0
            and _fields_equal(rich_n.transfer, exact_n.transfer)
        ):
            pytest.fail("N2N moment datum/ends drifted")

    def test_iso_pair_from_material_xs_equals_the_exact_ctor(self):
        mat_xs = _mat_xs()
        space = mat_xs.mesh.bulk_space
        rich = IsotropicScattering.from_material_xs(mat_xs, space=space)
        exact = IsotropicScattering(
            TransferMaterialField.scattering(mat_xs).at_order(0),
            domain=space, codomain=space,
        )
        if not (
            rich.domain is exact.domain
            and _fields_equal(rich.transfer, exact.transfer)
        ):
            pytest.fail("IsoS classmethod drifted from the ctor")
        rich_n = IsotropicN2N.from_material_xs(mat_xs, space=space)
        exact_n = IsotropicN2N(
            TransferMaterialField.n2n(mat_xs).at_order(0),
            domain=space, codomain=space,
        )
        if not (
            rich_n.domain is exact_n.domain
            and _fields_equal(rich_n.transfer, exact_n.transfer)
        ):
            pytest.fail("IsoN2N classmethod drifted from the ctor")

    def test_n2n_operator_from_solver_data_equals_the_exact_ctor(self):
        r"""#426 step 2: the (n,2n) binding is minted at the SOLVE's order on
        the SAME interned frame S is — the exact ctor is the core's (a field
        at that order + the two faces), and the P0 energy binding is DERIVED
        (``isotropic_energy``), of the role's own class."""
        mat_xs = _mat_xs()
        sn = _sn(mat_xs)
        space = sn.full_field_space
        rich = N2NOperator.from_solver_data(
            mat_xs=mat_xs, scattering_order=1, space=space,
        )
        interior = space.interior_space
        assert interior is not None
        frame = HarmonicFrame.for_space(interior, 1)
        exact = N2NOperator(
            TransferMaterialField.n2n(mat_xs).at_order(1),
            flux_analysis=frame.flux_analysis_on(interior),
            source_reconstruction=frame.source_reconstruction_on(interior),
            domain=space,
            codomain=space,
        )
        if not (rich.domain is exact.domain and rich.codomain is exact.codomain):
            pytest.fail("N2N ends drifted")
        if rich.frame is not exact.frame:
            pytest.fail("N2N frame drifted from the interned hub route")
        S = ScatteringOperator.from_solver_data(
            mat_xs=mat_xs, scattering_order=1, space=space,
        )
        if rich.frame is not S.frame:
            pytest.fail("the two transfer terms must share ONE interned frame")
        np.testing.assert_array_equal(
            np.asarray(rich.frame.measure.weights), np.asarray(sn.quad.weights),
        )
        if not _fields_equal(rich.transfer, exact.transfer):
            pytest.fail("the extracted (n,2n) kernel field drifted from the ctor's")
        if not (
            isinstance(rich.isotropic_energy, IsotropicN2N)
            and rich.isotropic_energy.transfer.order == 0
            and _fields_equal(
                rich.isotropic_energy.transfer, exact.transfer.at_order(0),
            )
        ):
            pytest.fail("the (n,2n) P0 energy binding drifted from the role's")

    @pytest.mark.parametrize("cls,ctor,fields", [
        (TransferKernel, TransferKernel.scattering, lambda mix: {
            "moments": tuple(np.asarray(s.todense()) for s in mix.SigS),
            "multiplicity": 1,
        }),
        (TransferKernel, TransferKernel.n2n, lambda mix: {
            "moments": tuple(np.asarray(s.todense()) for s in mix.Sig2),
            "multiplicity": N2N_MULTIPLICITY,
        }),
        (FissionKernel, FissionKernel.from_mixture,
         lambda mix: {"chi": mix.chi, "nu_sig_f": mix.SigP}),
    ], ids=["scattering", "n2n", "fission"])
    def test_kernel_from_mixture_equals_the_exact_ctor(self, cls, ctor, fields):
        """The datum tier's single parametrized row (plan §4: the
        extraction is one ``todense()`` per order plus the channel's yield;
        the existing value pins are cross-source gates, not this shape)."""
        mix = _mat_xs().materials[0]
        rich = ctor(mix)
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
    from orpheus.transport.operators.isotropic_transfer import (
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

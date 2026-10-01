r"""S1 gates for the carrier's angular-bulk space mint (campaign 1 CS4b).

G1.1–G1.5 of the CS4b verification plan
(``scratch/cs4b_verification_plan.md`` §11, step S1) plus the scheme-side
``moment_axis`` admission pair. The step is provably behaviour-neutral —
nothing consumes :attr:`SNProblem.angular_bulk_space` yet — so every gate here
is either a RECORD of the mint's content, a LAW comparing it against the
SHIPPED dense composite interior (the §6c witness that exists today), or an
ADMISSION with both legs (vv #11).

Conventions gated (CS4b crosswalk B1/B5, ``.claude/plans/cs4b_crosswalk.md``):

* the axis order is ``(angular, energy, spatial)``, matching the bulk tensor
  ``(N, ng, *spatial)``;
* the energy and spatial arms are ``bulk_space``'s axes REUSED VERBATIM
  (object identity — the energy-arm rule is spelled once);
* the Gram of the axis product equals the hand-built dense ``V·w_n``
  oracle on both the DD and the LD arm (LD composes the scheme-owned MODAL
  ``moment_axis`` carrying ``moment_mass_diagonal``), and since S2b the
  composite's interior IS the cached mint (identity on DD, ``==`` on LD);
* the derived space NAME is never pinned (R4: CS2's typed axis subclasses
  change the digest — every assertion here is per-axis content or relative
  ``is``/``==``);
* cone predicates: nodal bulk families answer ``True``, the harmonic-moment
  and trace families ``None``, the LD moment-tailed product ``False`` (the
  Q6-ratified routing base for ``cone_violations``; the exact vertex test
  is #400).

Fixture: the verification plan's §2 configuration — a NON-uniform 5-cell
slab (uniform volumes would collapse ``V`` to a scalar and blind half the
metric claims), ``gauss_legendre(4)``, ``ng = 2``, vacuum/vacuum.
"""

from __future__ import annotations

from fractions import Fraction

from tests._harness.float_bounds import gamma

import numpy as np
import pytest

from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellEdges, Mesher
from orpheus.mesh import AxisCoord, AxisMesh, RadialAxisMesh
from orpheus.numerics.axis import BasisKind, EnergyAxis
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.spatial import LinearDiscontinuous
from orpheus.transport.spatial.diamond import DiamondDifference
from tests.gates.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

#: NON-uniform edges — ``V = [0.2, 0.3, 0.4, 0.7, 1.4]``, a genuine vector.
_EDGES = np.array([0.0, 0.2, 0.5, 0.9, 1.6, 3.0])
_NG = 2


def _slab(*, scheme=None, ng: int = _NG) -> SNProblem:
    geometry = StructuredGeometry.slab(
        (_EDGES[0], _EDGES[-1]), (0,), left=BC("vacuum"), right=BC("vacuum"),
    )
    mesh = Mesher(geometry).partition(CellEdges(_EDGES)).mesh
    kwargs = {} if scheme is None else {"scheme": scheme}
    return SNProblem(
        mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=ng), **kwargs
    )


class TestG11AxisTuple:
    """G1.1 — the axis tuple IS (angular w_n, energy, spatial V). RECORD."""

    def test_axes_are_angular_energy_spatial_with_the_carrier_measures(self):
        sn = _slab()
        space = sn.angular_bulk_space
        assert space.axes is not None and len(space.axes) == 3
        angular, energy, spatial = space.axes

        assert angular.label == "angular"
        assert angular.shape == (sn.quad.N,)
        assert angular.kind is BasisKind.NODAL
        assert angular.weights is not None
        assert np.array_equal(angular.weights, sn.quad.weights)

        assert isinstance(energy, EnergyAxis)
        assert energy.shape == (sn.ng,)

        assert spatial.label == "spatial"
        assert spatial.shape == sn.spatial_shape
        assert spatial.kind is BasisKind.NODAL
        assert spatial.weights is not None
        assert np.array_equal(spatial.weights, sn.volumes)

        assert space.shape == (sn.quad.N, sn.ng, *sn.spatial_shape)

    def test_energy_and_spatial_arms_are_bulk_space_axes_verbatim(self):
        """The scalar arms are REUSED objects, not respelled twins — the
        energy-arm rule (``EnergyAxis.from_materials``) is spelled exactly
        once, in ``bulk_space`` (Pattern 2)."""
        sn = _slab()
        scalar_axes = sn.bulk_space.axes
        assert scalar_axes is not None
        assert sn.angular_bulk_space.axes is not None
        assert sn.angular_bulk_space.axes[1] is scalar_axes[0]
        assert sn.angular_bulk_space.axes[2] is scalar_axes[1]


class TestG12Cache:
    """G1.2 — the mint is CACHED. LAW (the is/== asymmetry IS the gate)."""

    def test_same_carrier_reads_the_same_instance(self):
        sn = _slab()
        assert sn.angular_bulk_space is sn.angular_bulk_space

    def test_twin_carriers_mint_equal_but_distinct_spaces(self):
        a, b = _slab(), _slab()
        assert a.angular_bulk_space == b.angular_bulk_space
        assert a.angular_bulk_space is not b.angular_bulk_space


class TestG13GramEquivalenceDD:
    """G1.3 (re-scoped at S2b) — the axis product's Gram equals the
    HAND-BUILT dense ``G_bulk = V·w_n``, and the composite interior IS the
    cached mint.

    Until the Q2 re-point this compared against
    ``full_field_space.interior_space``'s own dense array — the shipped
    §6c witness. S2b re-pointed that interior AT ``angular_bulk_space``,
    which made the comparison tautological (single-sourcing demotes every
    gate that compared the copies), so the dense side moved IN-TEST: a
    fuller-view oracle built from raw mesh data, independent of both
    production spellings. The identity row is the unification claim
    itself."""

    def test_composite_interior_is_the_cached_mint(self):
        sn = _slab()
        assert sn.full_field_space.interior_space is sn.angular_bulk_space

    @staticmethod
    def _exact(a) -> np.ndarray:
        return np.vectorize(lambda v: Fraction(float(v)), otypes=[object])(np.asarray(a, dtype=float))

    def test_all_three_metric_faces_agree_with_the_dense_oracle(self):
        r"""G1.3 — every metric face of the axis-built bulk space is :math:`G = V\cdot w_n`, within a DERIVED bound.

        The claim is the identity of the Gram, :math:`G_{n,g,i} = w_n V_i`; the
        two spellings (the axis product and the hand-densified oracle) are
        BOTH asserted against the EXACT bilinear form of the float inputs,
        computed in rational arithmetic — so neither spelling is the reference
        and the row is independent of how either associates.

        **Bounds, derived** (Higham, *Accuracy and Stability*, Lemma 3.1 and
        eq. (4.4); the standard model per operation, any association):

        * scalar :math:`\langle x, y\rangle_G = \sum_k w_n V_i x_k y_k` — ``n``
          terms, each a product of 4 factors (3 roundings in any order), summed
          in any order (at most ``n − 1`` additions per term):
          :math:`|\widehat{\langle x,y\rangle} - \langle x,y\rangle| \le
          \gamma_{n+2}\sum_k |w_n V_i x_k y_k|`;
        * norm :math:`\sqrt{\langle x, x\rangle_G}` — all terms positive, one
          more rounding for the square root, and
          :math:`|\sqrt{1+\theta}-1| \le |\theta|`: relative error
          :math:`\le \gamma_{n+3}`, asserted squared in exact arithmetic;
        * vector faces, per element: :math:`G x` is a 3-factor product, and
          :math:`G^{-1}x` is :math:`x/(wV)` or :math:`(x/w)/V` — two roundings
          either way: relative :math:`\gamma_2`.

        Until 2026-10 the scalar faces asserted ``==`` between the two
        spellings — two reduction orders of one inner product, equal only by the
        draw (red in 3 of 10 ±3-ULP perturbations of the GL-4 rule). The
        non-vacuity leg below shows the bound still separates the Gram from the
        two named wrong ones (the angular weight or the cell volume dropped) by
        orders of magnitude.
        """
        sn = _slab()
        axis_built = sn.angular_bulk_space
        w = np.asarray(sn.quad.weights, dtype=float)
        V = np.asarray(sn.volumes, dtype=float)
        dense = FunctionSpace(
            name="dense_oracle",
            shape=axis_built.shape,
            inner_product_weights=w.reshape(-1, 1, 1) * V.reshape(1, 1, -1),
        )

        rng = np.random.default_rng(0)
        x = rng.standard_normal(dense.shape)
        y = rng.standard_normal(dense.shape)

        W = self._exact(w).reshape(-1, 1, 1)
        VV = self._exact(V).reshape(1, 1, -1)
        X, Y = self._exact(x), self._exact(y)
        terms = W * VV * X * Y
        exact_ip = terms.sum()
        abs_sum = np.abs(terms).sum()
        n = terms.size
        ip_bound = gamma(n + 2) * abs_sum
        sq = (W * VV * X * X).sum()
        norm_rel = gamma(n + 3)

        for name, space in (("axis-built", axis_built), ("dense oracle", dense)):
            ip = Fraction(float(space.inner_product(x, y)))
            if abs(ip - exact_ip) > ip_bound:
                pytest.fail(
                    f"{name}: <x,y>_G = {float(ip)!r} is {float(abs(ip - exact_ip)):.3e} from "
                    f"the exact {float(exact_ip)!r}, outside gamma_(n+2)·Σ|terms| = {float(ip_bound):.3e}"
                )
            nrm = Fraction(float(space.norm(x)))
            if abs(nrm * nrm - sq) > sq * ((1 + norm_rel) ** 2 - 1):
                pytest.fail(f"{name}: ||x||_G = {float(nrm)!r} is not sqrt of the exact form within gamma_(n+3)")
            for face, got, want in (
                ("apply_metric", space.apply_metric(x), W * VV * X),
                ("apply_inverse_metric", space.apply_inverse_metric(x), X / (W * VV)),
            ):
                got_x = self._exact(np.broadcast_to(np.asarray(got, dtype=float), x.shape))
                err = np.abs(got_x - want)
                if np.any(err > gamma(2) * np.abs(want)):
                    k = int(np.argmax(err > gamma(2) * np.abs(want)))
                    pytest.fail(f"{name}.{face}: element {k} outside gamma_2 of the exact G-action")

        # Non-vacuity: the bound separates G = w·V from the two named wrong Grams.
        for wrong_name, wrong in (
            ("angular weight dropped", (VV * X * Y).sum()),
            ("cell volume dropped", (W * X * Y).sum()),
        ):
            if abs(wrong - exact_ip) <= 1000 * ip_bound:
                pytest.fail(f"the scalar bound cannot tell G from the Gram with the {wrong_name}")


class TestG14GramEquivalenceLD:
    """G1.4 (re-scoped at S2b) — the LD arm: the Gram carries the scheme's
    moment mass on the trailing 2^d axis; the axis form reproduces the
    hand-built dense oracle, and the composite interior is the widened
    product. LAW.

    [M] R9 measured its draw's inner product bit-identical; that was the
    draw's luck, not a law — the two spellings associate the weight
    products differently, and on THIS fixture's ``rng(0)`` draw the
    near-cancelling bilinear form lands 6 ULP apart (measured 2026-08-22).
    A hand-picked ``nulp ≤ 64`` stood here until 2026-10-01; it was itself a
    rounding accident (red under ±3-ULP jitter of the rule). Both spellings
    are now checked against the exact rational form within derived bounds
    (γ_{n+3} for the scalar face, γ_{n+4} for the norm, γ_3 for the vector
    faces), and a non-vacuity leg shows the bound still tells the Gram from
    versions with the moment mass, the angular weight or the cell volume
    dropped."""

    def test_ld_composite_interior_is_the_widened_product(self):
        """The unification claim on the LD arm: the composite's interior
        IS the cached trial mint (CS4b S5 upgraded this from ``==`` to
        ``is`` — the widening composition moved from an inline
        ``of_axes`` here into :attr:`SNProblem.angular_trial_space`, so the
        composite, the trial property, and every LD allocation share ONE
        instance), and that mint equals the moment-widened product of
        the cached base."""
        sn = _slab(scheme=LinearDiscontinuous())
        assert sn.angular_bulk_space.axes is not None
        widened = FunctionSpace.of_axes(
            *sn.angular_bulk_space.axes, sn.scheme.moment_axis(sn.axes)
        )
        assert sn.full_field_space.interior_space == widened
        assert sn.full_field_space.interior_space is sn.angular_trial_space

    def test_all_three_metric_faces_agree_on_the_ld_interior(self):
        r"""G1.4 — every metric face of the LD interior is :math:`G = w_n V_i m_j`, within a DERIVED bound.

        The claim is the identity of the widened Gram,
        :math:`G_{n,g,i,j} = w_n V_i m_j` with :math:`m` the scheme's moment
        mass. The two spellings (the axis product widened by the scheme's
        modal ``moment_axis``, and the hand-densified oracle built from raw
        mesh and scheme data) are BOTH asserted against the EXACT bilinear
        form of the float inputs, computed in rational arithmetic, so neither
        spelling is the reference and the row is independent of how either
        associates. The method is G1.3's; the only change is the fifth factor.

        **Bounds, derived** (Higham, *Accuracy and Stability*, Lemma 3.1 and
        eq. (4.4); the standard model per operation, any association):

        * scalar :math:`\langle x, y\rangle_G = \sum_k w_n V_i m_j x_k y_k` —
          ``n`` terms, each a product of 5 factors (4 roundings in any order),
          summed in any order: :math:`\le \gamma_{n+3}\sum_k |\text{term}_k|`;
        * norm — positive terms, one more rounding for the square root:
          relative :math:`\le \gamma_{n+4}`, asserted squared in exact
          arithmetic;
        * vector faces, per element: :math:`G x` is a 4-factor product and
          :math:`G^{-1}x` is :math:`x` divided by three factors in some
          bracketing — three roundings either way: relative :math:`\gamma_3`.

        Until 2026-10 the scalar face asserted ``nulp <= 64`` between the two
        spellings and the other faces ``nulp <= 4``: hand-picked counts on a
        cancellation-conditioned sum, with no derivation behind either number.
        The non-vacuity leg shows the bound still separates the Gram from the
        three named wrong ones (the moment mass, the angular weight or the cell
        volume dropped) by more than three orders of magnitude.
        """
        sn = _slab(scheme=LinearDiscontinuous())
        base = sn.angular_bulk_space
        assert base.axes is not None
        widened = FunctionSpace.of_axes(
            *base.axes, sn.scheme.moment_axis(sn.axes)
        )
        # The oracle: G_bulk = V·w_n ⊗ moment_mass, densified BY HAND from
        # raw mesh + scheme data (the retired production spelling, now the
        # test-side fuller-view reference).
        w = np.asarray(sn.quad.weights, dtype=float)
        V = np.asarray(sn.volumes, dtype=float)
        mass = np.asarray(sn.scheme.moment_mass_diagonal(sn.axes), dtype=float)
        g = (w.reshape(-1, 1, 1) * V.reshape(1, 1, -1))[..., None] * mass
        dense = FunctionSpace(
            name="dense_oracle_ld",
            shape=widened.shape,
            inner_product_weights=g,
        )
        if widened.shape != dense.shape:
            pytest.fail(f"widened shape {widened.shape} != oracle shape {dense.shape}")

        rng = np.random.default_rng(0)
        x = rng.standard_normal(dense.shape)
        y = rng.standard_normal(dense.shape)

        exact = TestG13GramEquivalenceDD._exact
        W = exact(w).reshape(-1, 1, 1, 1)
        VV = exact(V).reshape(1, 1, -1, 1)
        M = exact(mass).reshape(1, 1, 1, -1)
        X, Y = exact(x), exact(y)
        G = W * VV * M
        terms = G * X * Y
        exact_ip = terms.sum()
        n = terms.size
        ip_bound = gamma(n + 3) * np.abs(terms).sum()
        sq = (G * X * X).sum()
        norm_rel = gamma(n + 4)

        for name, space in (("axis-built", widened), ("dense oracle", dense)):
            ip = Fraction(float(space.inner_product(x, y)))
            if abs(ip - exact_ip) > ip_bound:
                pytest.fail(
                    f"{name}: <x,y>_G = {float(ip)!r} is {float(abs(ip - exact_ip)):.3e} from "
                    f"the exact {float(exact_ip)!r}, outside gamma_(n+3)·Σ|terms| = {float(ip_bound):.3e}"
                )
            nrm = Fraction(float(space.norm(x)))
            if abs(nrm * nrm - sq) > sq * ((1 + norm_rel) ** 2 - 1):
                pytest.fail(f"{name}: ||x||_G = {float(nrm)!r} is not sqrt of the exact form within gamma_(n+4)")
            for face, got, want in (
                ("apply_metric", space.apply_metric(x), G * X),
                ("apply_inverse_metric", space.apply_inverse_metric(x), X / G),
            ):
                got_x = exact(np.broadcast_to(np.asarray(got, dtype=float), x.shape))
                outside = np.abs(got_x - want) > gamma(3) * np.abs(want)
                if np.any(outside):
                    k = int(np.argmax(outside))
                    pytest.fail(f"{name}.{face}: element {k} outside gamma_3 of the exact G-action")

        # Non-vacuity: the bound separates G = w·V·m from the three named wrong Grams.
        for wrong_name, wrong in (
            ("moment mass dropped", (W * VV * X * Y).sum()),
            ("angular weight dropped", (VV * M * X * Y).sum()),
            ("cell volume dropped", (W * M * X * Y).sum()),
        ):
            if abs(wrong - exact_ip) <= 1000 * ip_bound:
                pytest.fail(f"the scalar bound cannot tell G from the Gram with the {wrong_name}")

def _cart_axes():
    # A minimal 1-D Cartesian axis tuple (P4.6: the family consumes axes).
    return (AxisMesh(edges=np.array([0.0, 1.0]),
                     bc_low=BC.reflective, bc_high=BC.reflective),)


def _radial_axes(kind: AxisCoord):
    # A minimal 1-D radial axis tuple of the given kind (P4.6).
    return (RadialAxisMesh(edges=np.array([0.0, 1.0]), coord=kind,
                           bc_outer=BC.reflective),)


class TestMomentAxisAdmission:
    """The scheme-side mint's ADMISSION pair (vv #11: both legs)."""

    def test_ld_mints_the_modal_mass_axis(self):
        scheme = LinearDiscontinuous()
        axis = scheme.moment_axis(_cart_axes())
        assert axis.label == "spatial_moment"
        assert axis.shape == (2,)
        assert axis.kind is BasisKind.MODAL
        assert axis.weights is not None
        assert np.array_equal(
            axis.weights,
            scheme.moment_mass_diagonal(_cart_axes()),
        )

    def test_slopeless_closure_refuses(self):
        with pytest.raises(NotImplementedError, match="no moment axis"):
            DiamondDifference().moment_axis(_cart_axes())

    # ── The CHART admission (2026-08-26).  Third arm of the same pair:
    # a multi-moment mass is defined on a Cartesian chart and is NOT
    # expressible on a curved one, so the producer must refuse rather
    # than hand back the slab's diagonal.
    #
    # ⭐ §6c — THE WITNESS IS CONSTRUCTIBLE, and that is the point of this
    # class of gate.  Before the guard, `SNProblem(Mesh1D(coord=SPHERICAL),
    # gauss_legendre(4), ..., scheme=LinearDiscontinuous())` BUILT, and its
    # moment weights measured [1., 0.33333333] -- bit-identical to a slab's,
    # on both the sphere AND the cylinder.  The wrong value was being
    # installed on two shipped charts, silently.  A gate that only mutated
    # the SUT would have proved nothing about that.

    @pytest.mark.parametrize(
        "kind", [AxisCoord.RADIAL_SPHERICAL, AxisCoord.RADIAL_CYLINDRICAL]
    )
    def test_curvilinear_multi_moment_mass_is_refused_not_slab_defaulted(
        self, kind: AxisCoord,
    ) -> None:
        """LD on a curved chart REFUSES; before the guard it returned the slab's.

        The true ``M/V`` there is cell-dependent AND non-diagonal (a
        spherical pole cell wants ``[[1, 0.5], [0.5, 0.4]]``), which a
        per-axis ``Axis`` weight vector cannot express -- so no honest
        value exists to return.  The MACHINERY half of the old two-blocker
        wording (#409, the non-Hadamard metric) was discharged by P7 (the
        dense-metric family); the refusal stands on the VALUE alone --
        #158's cell solve is what gives a chosen ``G`` a consumer.
        """
        with pytest.raises(NotImplementedError, match="no moment mass"):
            LinearDiscontinuous().moment_mass_diagonal(_radial_axes(kind))
        with pytest.raises(NotImplementedError, match="no moment mass"):
            LinearDiscontinuous().moment_axis(_radial_axes(kind))

    @pytest.mark.parametrize(
        "kind", [AxisCoord.RADIAL_SPHERICAL, AxisCoord.RADIAL_CYLINDRICAL]
    )
    def test_the_moment_mass_refusal_names_only_the_value_blocker(
        self, kind: AxisCoord,
    ) -> None:
        """E1 (P7 S4): the refusal stands on ONE blocker after the
        dense-metric family landed.

        The message still refuses (the pinned ``no moment mass`` fragment
        survives), still names #158 (the value's missing consumer), and
        no longer names #409 — the machinery half was discharged by P7.
        The absence assert is the half a ``match=`` cannot pin: it stops
        the two-blocker wording drifting back while reading as a mere
        rewording.
        """
        with pytest.raises(NotImplementedError, match="no moment mass") as exc:
            LinearDiscontinuous().moment_mass_diagonal(_radial_axes(kind))
        message = str(exc.value)
        assert "158" in message, "the VALUE blocker (#158) must stay named"
        assert "409" not in message, (
            "the discharged MACHINERY blocker (#409) must not be re-cited"
        )

    def test_mixed_axes_refuse_and_name_only_the_curved_kind(self):
        """P4.6's granularity witness: a mixed (z, r)-style tuple refuses,
        NAMING only the radial axis kind — the per-axis question the
        whole-mesh enum structurally could not pose (its projection
        refuses mixed multi-axis tuples outright, ``coord_system`` at
        ``transport/mesh/axis.py``).  No mesh ctor builds this today;
        the bare-axes spelling is the constructible witness (§6c).
        """
        axes = (
            _cart_axes()[0],
            _radial_axes(AxisCoord.RADIAL_CYLINDRICAL)[0],
        )
        with pytest.raises(
            NotImplementedError, match="mass on a radial_cylindrical axis",
        ):
            LinearDiscontinuous().moment_mass_diagonal(axes)

    @pytest.mark.parametrize(
        "kind", [AxisCoord.RADIAL_SPHERICAL, AxisCoord.RADIAL_CYLINDRICAL]
    )
    def test_slopeless_mass_is_admitted_on_a_curved_chart(
        self, kind: AxisCoord,
    ) -> None:
        """The width control: the guard must not be too WIDE.

        A single-moment scheme's cell-average mass is :math:`V/V = 1`
        whatever the measure, so DD is unaffected by the chart.  Without
        this leg the refusal above is compatible with a guard that simply
        rejects every curvilinear chart.
        """
        assert np.array_equal(
            DiamondDifference().moment_mass_diagonal(_radial_axes(kind)),
            np.ones(1),
        )


class TestAngularTrialSpace:
    """``SNProblem.angular_trial_space`` — the ONE widening mint (CS4b S5).

    The property replaces the retired ``spatial_moments=`` factory int:
    a call site widens by SELECTING this mint instead of threading the
    scheme's basis size through a factory parameter. Claims: the
    slopeless identity (LAW — the two mints collapse to one instance),
    the LD structure (RECORD — base axes + the scheme's moment axis),
    and the single-source identity with the composite interior (LAW).
    """

    def test_slopeless_trial_space_IS_the_bulk_space(self):
        """DD (width 1): not merely ``==`` — the SAME cached instance,
        so slopeless consumers pay nothing and the mints cannot drift."""
        sn = _slab(scheme=DiamondDifference())
        assert sn.angular_trial_space is sn.angular_bulk_space

    def test_ld_trial_space_appends_the_scheme_moment_axis(self):
        sn = _slab(scheme=LinearDiscontinuous())
        base = sn.angular_bulk_space
        trial = sn.angular_trial_space
        assert base.axes is not None and trial.axes is not None
        # The base's axes verbatim, then the scheme's own moment axis.
        assert trial.axes[: len(base.axes)] == base.axes
        (tail,) = trial.axes[len(base.axes) :]
        assert tail == sn.scheme.moment_axis(sn.axes)
        assert trial.shape == (*base.shape, 2)

    def test_ld_trial_space_is_cached_and_single_sourced(self):
        sn = _slab(scheme=LinearDiscontinuous())
        assert sn.angular_trial_space is sn.angular_trial_space
        assert sn.full_field_space.interior_space is sn.angular_trial_space

    def test_field_allocation_rides_the_trial_mint(self):
        """The S5 end-state, at both widths: a field allocated on the
        trial mint IS an element of it — the DD leg doubles as the
        collapse witness (its trial mint IS the bulk instance). Until
        the sugar retired (S5.4) this row was the BRIDGE gate: ``[M]``
        2026-08-24, ``zeros_on(mesh, spatial_moments=…)``'s derived
        space was ``is``-identical to this mint at DD and ``==`` at LD,
        proving the ~700-site migration a pure re-spelling."""
        from orpheus.transport.fields.angular_flux import AngularFlux

        dd = _slab(scheme=DiamondDifference())
        assert AngularFlux.zeros(dd.angular_trial_space).space is dd.angular_bulk_space
        ld = _slab(scheme=LinearDiscontinuous())
        psi = AngularFlux.zeros(ld.angular_trial_space)
        assert psi.space is ld.angular_trial_space
        assert psi.values.shape == (*ld.angular_bulk_space.shape, 2)


class TestG15ConePredicates:
    """G1.5 — the cone predicates, stated. RECORD ([M] verification plan
    Finding 8 + R9)."""

    def test_nodal_bulk_families_answer_true(self):
        sn = _slab()
        assert sn.angular_bulk_space.has_coordinate_cone is True
        assert sn.bulk_space.has_coordinate_cone is True

    def test_the_harmonic_moment_family_answers_false_and_the_trace_family_none(self):
        """The moment head is a MODAL axis since CS4c step 6 item 6.2c-ii
        (a spectral coefficient may be negative for a positive function), so
        the moment product's cone answer is a definite ``False`` — the typed
        refusal :meth:`Field.cone_violations` turns it into — where the
        axes-less head answered ``None`` (unanswerable). The trace family is
        still name-built and still answers ``None``."""
        sn = _slab()
        moment_space = HarmonicMomentFlux.zeros_for_problem_and_L(sn, 1).space
        assert moment_space.has_coordinate_cone is False
        assert sn.angular_trace.has_coordinate_cone is None

    def test_the_ld_moment_tail_is_modal_so_the_cone_reads_false(self):
        """The Q6-ratified routing base: a moment-tailed LD bulk space
        answers ``False`` (signed slope coefficients are legal on a
        positive function), which ``Field.cone_violations`` turns into the
        typed refusal. The exact modal test (the vertex theorem) is #400."""
        sn = _slab(scheme=LinearDiscontinuous())
        assert sn.angular_bulk_space.axes is not None
        widened = FunctionSpace.of_axes(
            *sn.angular_bulk_space.axes, sn.scheme.moment_axis(sn.axes)
        )
        assert widened.has_coordinate_cone is False

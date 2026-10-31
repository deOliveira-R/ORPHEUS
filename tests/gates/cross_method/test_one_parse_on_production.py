r"""One parse at the domain, on production objects: the 20 rows that were silent under ``python -O`` refused on
every run, and every admission spelling reaching the ONE primitive (U1b, plan ``space_membership_covector_gauge.md``).

It lives in ``cross_method`` because its rows span the homogeneous, SN and diffusion objects; the numerics-level
rows it rests on are ``tests/gates/numerics/test_coordinate_space.py``.

Every refusal is a :class:`ValueError` (``NotInCoordinateSpace`` subclasses it alone), matched on the primitive's
fragment: ``exact shape`` for a scalar verb or an engine minted with its domain shape, ``leading shape`` for an
action ``A ⊗ I_P``.

First reds: every ``TestTheSilentRowsAreRefused`` row RETURNS under ``python -O`` at ``6aa26d00``, the last commit
before U1b (`[M]` 2026-10-30, this file run against a ``git archive 6aa26d00`` copy outside the repository: 20 of 20
``DID NOT RAISE``, the 2 honest rows green there too). The routing gate's first red is that a site keeping its own
comparison returns under the patched primitive (``DID NOT RAISE _Routed``).
"""
from __future__ import annotations

import sys

import numpy as np
import pytest

pytestmark = [
    pytest.mark.l0,
    pytest.mark.rests_on(
        "tests/gates/numerics/test_coordinate_space.py::TestThePrimitive",
        "tests/gates/numerics/test_coordinate_space.py::TestTheSpaceVerbsParseOnEveryRun",
        "tests/gates/numerics/test_coordinate_space.py::TestAnEngineIsMintedWithItsDomainShape",
    ),
]

EXACT, LEADING = "exact shape", "leading shape"


# ── the fixtures (the explorer's probe builders, kept identical so the first reds are the probe's) ────────

@pytest.fixture(scope="module")
def homog():
    from orpheus.data.macro_xs.synthetic import get_mixture
    from orpheus.homogeneous.solver import HomogeneousProblem
    return HomogeneousProblem(get_mixture("A", "2g"))


@pytest.fixture(scope="module")
def sn():
    from orpheus.sn import solver as sns
    from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
    from tests.gates.sn.architecture import _sn_consumers_helpers as S
    from tests.gates.transport import _carry_cases as carry
    mesh, quad, conds = S.build("slab", ("vacuum",))
    prob = sns._as_problem(mesh, quad, carry.materials(), boundary_conditions=conds)
    solver = sns.SNSolver(prob, inner_solver="source_iteration")
    phi = np.asarray(solver.initial_flux_distribution(), dtype=float)
    phi = phi * (1 + np.random.default_rng(0).random(phi.shape))
    collision = MultiplicationOperator.from_mesh(np.asarray(prob.mat_xs.total_cross_section), prob)
    return dict(prob=prob, solver=solver, phi=phi, bulk=prob.bulk_space, engine=collision.engine,
                n_ord=len(np.asarray(quad.weights)), mesh=mesh, quad=quad, conds=conds, carry=carry)


@pytest.fixture(scope="module")
def sn_solution(sn):
    from orpheus.sn.solver import solve_sn
    return solve_sn(sn["carry"].materials(), sn["mesh"], sn["quad"], boundary_conditions=sn["conds"],
                    max_outer=30, keff_tol=1e-5, flux_tol=1e-4)


@pytest.fixture(scope="module")
def diffusion_bulk():
    from orpheus.diffusion.solver import DiffusionSolver
    from tests.gates.diffusion import _wave2_helpers as H
    ds = DiffusionSolver(H.diffusion_mesh(*H.build("slab", ("vacuum",))))
    psi = ds.unflatten(ds.initial_flux_distribution())
    return psi.interior.space, np.asarray(psi.interior.values)


# ── the 20 rows silent under -O at 6aa26d00 ─────────────────────────────────────────────────────────────────

class TestTheSilentRowsAreRefused:
    # (i) outer products: 5 rows
    @pytest.mark.catches("ERR-122")
    def test_the_reaction_rate_refuses_a_flat_carrier(self, homog):
        """ERR-122's own input, refused under ``-O`` (it was the debug-only row of
        ``test_reaction_rate_on_the_pose.py``, #669). `[M]` at 6aa26d00 it read 0.225 (0.095)."""
        with pytest.raises(ValueError, match=EXACT):
            homog.production_rate.evaluate(np.array([0.6, 0.4]))

    def test_the_covector_refuses_a_flat_carrier(self, homog):          # [M] 0.225 (0.095)
        with pytest.raises(ValueError, match=EXACT):
            homog.production_rate.covector.evaluate(np.array([0.6, 0.4]))

    def test_the_dyad_functional_refuses_a_flat_carrier(self, homog):
        """The rank-one fission dyad's ``InnerProductFunctional`` (the IPF raw residue). `[M]` (1, 2), sum 0.225."""
        ipf = homog.production.production_rate      # a ReactionRateFunctional, an InnerProductFunctional
        with pytest.raises(ValueError, match=EXACT):
            ipf.evaluate(np.array([0.6, 0.4]))

    def test_the_homogeneous_inner_product_refuses_a_flat_operand(self, homog):     # [M] 1.0 (0.52)
        col = np.array([[0.6], [0.4]])
        with pytest.raises(ValueError, match=EXACT):
            homog.space.inner_product(col, col.ravel())

    def test_the_homogeneous_norm_refuses_a_square_operand(self, homog):            # [M] 1.0198 (0.7211)
        with pytest.raises(ValueError, match=EXACT):
            homog.space.norm(np.ones((2, 2)))

    # (ii) implicit pullbacks x ⊗ 1: 11 rows
    @pytest.mark.parametrize("verb", ["inner_product", "integral", "flat"])
    def test_the_sn_bulk_scalar_verbs_refuse_a_column(self, sn, verb):
        """`[M]` inner_product 15.58 (13.91); integral 10.48 (8.93); flat = flat(φ⊗1)."""
        B, phi = sn["bulk"], sn["phi"]
        col = phi[:, :1]
        call = {"inner_product": lambda: B.inner_product(phi, col), "integral": lambda: B.integral.evaluate(col),
                "flat": lambda: B.flat(col)}[verb]
        with pytest.raises(ValueError, match=EXACT):
            call()

    @pytest.mark.parametrize("verb", ["apply_metric", "apply_inverse_metric"])
    def test_the_sn_bulk_metric_legs_refuse_a_column(self, sn, verb):     # [M] G(φ⊗1), G⁺(φ⊗1)
        with pytest.raises(ValueError, match=LEADING):
            getattr(sn["bulk"], verb)(sn["phi"][:, :1])

    def test_the_diffusion_bulk_inner_product_refuses_a_column(self, diffusion_bulk):   # [M] 13.41 (15.19)
        space, x = diffusion_bulk
        with pytest.raises(ValueError, match=EXACT):
            space.inner_product(x, x[:, :1])

    @pytest.mark.parametrize("verb", ["apply", "solve"])
    def test_the_sn_collision_engine_refuses_a_spatial_column(self, sn, verb):
        """Silent with asserts on as well (`[M]` probe_debug.log): the rank check admits (8, 2, 1) for (8, 2, 12)."""
        x = np.ones((sn["n_ord"], 2, 1))
        with pytest.raises(ValueError, match=EXACT):
            getattr(sn["engine"], verb)(x)

    def test_the_sn_absorption_rate_refuses_a_column(self, sn):          # [M] Σ 0.516 (0.424), asserts on too
        with pytest.raises(ValueError, match=EXACT):
            sn["solver"].compute_group_absorption_rate(sn["phi"][:, :1])

    def test_the_reaction_rate_density_refuses_a_column_cross_section(self, sn, sn_solution):   # [M] Σ 4.10 (3.35)
        xs = np.asarray(sn["prob"].mat_xs.absorption_cross_section)
        with pytest.raises(ValueError, match=EXACT):
            sn_solution.reaction_rate_density(xs[:, :1])

    def test_the_frame_analysis_face_refuses_one_ordinate_for_eight(self, sn):
        """`[M]` Σ 26.42 (24.72): one ordinate taken for all 8, silent with asserts on."""
        from orpheus.numerics.frame_faces import FrameAnalysis
        M = FrameAnalysis.of(sn["prob"].system.factors.fission.source_reconstruction.frame)
        with pytest.raises(ValueError, match=LEADING):
            M.apply(np.ones((1, 2, 12)))

    # (iii) non-members out: 4 rows
    def test_the_sn_metric_leg_refuses_a_one_group_row(self, sn):          # [M] out (1, 12)
        with pytest.raises(ValueError, match=LEADING):
            sn["bulk"].apply_metric(sn["phi"][:1])

    def test_the_sn_collision_engine_refuses_one_ordinate(self, sn):      # [M] out (1, 2, 12), asserts on too
        with pytest.raises(ValueError, match=EXACT):
            sn["engine"].apply(np.ones((1, 2, 12)))

    def test_the_homogeneous_engine_refuses_a_payload(self, homog):        # [M] out (2, 3), asserts on too
        with pytest.raises(ValueError, match=EXACT):
            homog.collision.engine.apply(np.ones((2, 3)))

    def test_the_volume_measure_refuses_a_row(self, sn):                    # [M] out (12,), asserts on too
        vm = sn["mesh"].volume_measure
        with pytest.raises(ValueError, match=LEADING):
            vm.integrate(np.ones((1, vm.n_points)))


class TestTheRouteOnlySitesRefuse:
    """Two sites the route gate alone guarded (qa, U1b review: their mutation arms reddened only the route row). The
    production rate integrates through the volume measure, which refused before U1b; ``SourcePosing`` compared spaces
    only for a typed source. A refusal row each, so the site's own parse has a witness that is not its owner tag."""

    def test_the_sn_production_rate_refuses_a_column(self, sn):
        with pytest.raises(ValueError, match=EXACT):
            sn["solver"].compute_group_production_rate(sn["phi"][:, :1])

    def test_a_bare_source_is_parsed_at_the_codomain(self, homog):
        from orpheus.numerics.pencil import OperatorPencil
        from orpheus.numerics.posing import SourcePosing
        prompt = OperatorPencil(homog.loss, homog.production).at(1.0)
        with pytest.raises(ValueError, match=EXACT):
            SourcePosing(prompt, np.ones(2))


class TestTheHonestCallsAreUnchanged:
    """The positive leg of the rows above: the correct operand reads the value it read at ``6aa26d00`` bit for bit
    (the U1b head capture; the parse returns its operand itself, so no product or reduction tree moves)."""

    def test_the_reaction_rate(self, homog):
        assert homog.production_rate.evaluate(np.array([[0.6], [0.4]])) == 0.09500000000000001   # [M] head capture

    def test_the_engine_keeps_its_payload_free_product(self, sn):
        x = np.random.default_rng(1).random((sn["n_ord"], 2, 12))
        np.testing.assert_array_equal(sn["engine"].apply(x), np.expand_dims(sn["engine"].coefficient, (0,)) * x)

    def test_a_site_holding_its_space_admits_a_typed_member(self, sn):
        """D1 at a bound site: the axis marginal asks its own space, so a member of that space is read by its
        coordinates, equal to the bare call (`[M]` 2026-10-30, at 4882293a: refused, ``values.shape ()``, because
        the shape-only parse saw the member as a 0-d object array)."""
        from orpheus.numerics.axis_marginal import AxisRetractionOperator
        from orpheus.numerics.space import Member
        retraction = AxisRetractionOperator.of(sn["prob"].angular_bulk_space, "angular")
        x = np.random.default_rng(2).random(retraction.domain.shape)
        np.testing.assert_array_equal(retraction.apply(Member(values=x, space=retraction.domain)), retraction.apply(x))


# ── every admission spelling reaches the ONE primitive (the mutation of the one verb reddens them all) ────

class _Routed(Exception):
    """Raised by the patched primitive: a site that propagates it routed its shape leg through the primitive."""


def _patch_primitive(monkeypatch):
    from orpheus.numerics import coordinate_space as cs
    originals = {cs.coordinates, cs.leading_coordinates}

    def routed(*_a, owner="", **_k):
        raise _Routed(owner)

    n = 0
    for mod in list(sys.modules.values()):
        for name in ("coordinates", "leading_coordinates"):
            if mod is not None and getattr(mod, name, None) in originals:
                monkeypatch.setattr(mod, name, routed)
                n += 1
    assert n >= 2, f"the primitive was rebound in {n} bindings: the arm did not install"
    return n


#: The correct-operand call of each admission spelling and each parse site, resolved by NAME at run time (a
#: rebinding arm cannot reach a row parametrized over objects built at collection, lessons L100). 25 independent
#: sites: the space's two verbs, ``Member``, ``admit_array``, the axis marginal, the objects below the space, the
#: transport and SN operators, and the solver and solution observables. ``harmonic_frame._admit`` is not a row: it
#: is typed-only, with no shape leg of its own, so it could see the primitive only through the carrier's
#: ``Member.__post_init__``, a reuse of that row (U1b spec §3.4).
_SITES = [
    "FunctionSpace.coordinates_of", "FunctionSpace.leading_coordinates_of", "Member.__post_init__",
    "admit_array", "AxisRetractionOperator.apply", "FissionMomentOperator._admit",
    "DiagonalOperator.apply", "DiagonalMetric.apply_block", "DenseMetric.apply_block", "FactoredMetric.apply",
    "DiscreteMeasure.integrate", "InnerProductFunctional.evaluate", "Basis.reconstruct", "Basis.analyze",
    "FrameAnalysis.apply", "FrameReconstruction.apply", "CrossGramInverse.apply",
    "SNSolver.compute_group_absorption_rate", "SNSolver.compute_group_production_rate",
    "Solution.reaction_rate_density", "PartialCurrentOperator.apply", "IsotropicEmissionOperator.apply_transpose",
    "MaterialXSField.project_through_bilinear", "SourcePosing.__post_init__",
    "RadialCharacteristicReconstruction.apply",
]

#: The two composite parses compare block spaces by content and carry no shape leg: a composite is not a flat
#: member, and routing its identity through one parse is #666, not U1b. Declared, measured exceptions.
_COMPOSITE_SITES = ["FullField.require_member", "RadialCharacteristicField.require_member"]


def _route(call, monkeypatch):
    call()  # honest once, unpatched: every lazy member (a cached flux, a cached field) is built before the patch
    _patch_primitive(monkeypatch)
    with pytest.raises(_Routed) as routed:
        call()
    return routed.value.args[0]


@pytest.mark.parametrize("site", _SITES)
def test_every_admission_reaches_the_one_primitive(site, monkeypatch, homog, sn, sn_solution, sphere):
    """The mutation of the one verb reddens every spelling: if a site keeps its own shape comparison, the patched
    primitive is never called and the site returns (red here). Correct operands only, so no site can raise for its
    own reason first; every object is built, and the call run once honestly, BEFORE the patch, so a constructor's or
    a lazy member's parse cannot be credited to the verb under test. The FIRST parse reached must be the site's own
    (the ``owner`` it passes): without that, a site keeping its own comparison is still credited by any parse
    downstream of it (`[M]` 2026-10-30, against the tree before U1b step 7: five rows of this step's sites green that
    way, the group rates through the measure, the moment operator through the dyad). A site that asks a space's verb
    names ``repr(space)`` as its owner, so for those rows the owner check proves "the first parse is at this space",
    not "at this site" (qa, U1b review): the site itself is witnessed by its refusal row above (the group rates, the
    reaction-rate density, ``SourcePosing``, the axis retraction's typed-member row)."""
    call, expected = _site_calls(site, homog, sn, sn_solution, sphere)
    owner = _route(call, monkeypatch)
    assert owner == expected, (
        f"{site}: the first parse reached is {owner!r}, not the site's own {expected!r} (a downstream parse "
        f"credited to this site)"
    )


@pytest.mark.xfail(strict=True, raises=pytest.fail.Exception,
                   reason="#666: a composite's identity is not yet one parse; its require_member compares block "
                          "spaces and has no shape leg (an XPASS means #666 landed: move the row into _SITES)")
@pytest.mark.parametrize("site", _COMPOSITE_SITES)
def test_the_composite_parses_reach_the_one_primitive(site, monkeypatch, sphere):
    """`[M]` 2026-10-30: both return a correct composite under the patched primitive (``DID NOT RAISE _Routed``)."""
    _route(_composite_calls(site, sphere), monkeypatch)


def test_the_declared_exceptions_are_strict():
    """The composite rows' xfail must be strict, or an XPASS after #666 would pass in silence."""
    marks = [m for m in test_the_composite_parses_reach_the_one_primitive.pytestmark if m.name == "xfail"]
    assert len(marks) == 1 and marks[0].kwargs.get("strict") is True


@pytest.fixture(scope="module")
def sphere():
    """A seed-carrying sphere (GL S4, 5 cells, 2 groups) and its ψ½ fold, from the ψ½ coupling gates."""
    from tests.gates.sn.operators import test_psi_half_coupling as H
    problem = H._sphere()
    return dict(problem=problem, fold=H._rc_fold(problem), H=H)


def _composite_calls(site, sphere):
    from orpheus.transport.full_field import FullField
    from orpheus.transport.radial_characteristic_field import RadialCharacteristicField
    problem, H = sphere["problem"], sphere["H"]
    space = problem.radial_characteristic_field_space
    if site == "FullField.require_member":
        x = H._zero_composite(problem)
        return lambda: FullField.require_member(x, mesh=problem, context="routing gate")
    x = RadialCharacteristicField.flux_zeros(space)
    return lambda: RadialCharacteristicField.require_member(x, space=space, context="routing gate")


def _site_calls(site, homog, sn, sn_solution, sphere):
    """The correct-operand call for one site and the ``owner`` its own parse names, every object built here (before
    the patch). Builders: the explorer's
    probe; the parse domain of each site is the U1b spec's per-object table."""
    from orpheus.numerics.axis_marginal import AxisRetractionOperator
    from orpheus.numerics.basis.legendre_basis import LegendreBasis
    from orpheus.numerics.frame import GalerkinFrame
    from orpheus.numerics.frame_faces import CrossGramInverse, FrameAnalysis, FrameReconstruction
    from orpheus.numerics.metric import DenseMetric, DiagonalMetric, FactoredMetric
    from orpheus.numerics.pencil import OperatorPencil
    from orpheus.numerics.posing import SourcePosing
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import Member
    from orpheus.sn.boundary.angular import IsotropicEmissionOperator, PartialCurrentOperator
    from orpheus.transport.operators.fission import FissionMomentOperator
    prob, B, phi = sn["prob"], sn["bulk"], sn["phi"]
    q4 = Quadrature.gauss_legendre(n_ordinates=4)
    lb = LegendreBasis(2)
    gf = GalerkinFrame(lb, q4.measure)
    fr = prob.system.factors.fission.source_reconstruction.frame
    diagonal, dense = DiagonalMetric(np.ones(12)), DenseMetric(np.eye(2))
    factored = FactoredMetric((((2,), None), ((12,), DiagonalMetric(np.ones(12)))))
    analysis, reconstruction, cross_gram = FrameAnalysis.of(fr), FrameReconstruction.of(fr), CrossGramInverse.of(gf)
    retraction = AxisRetractionOperator.of(prob.angular_bulk_space, "angular")
    moment_op = _find(prob.system.factors.fission.full_fission_kernel, FissionMomentOperator)
    assert moment_op is not None, "the fission kernel holds no FissionMomentOperator: the builder lost its subject"
    moments = np.ones(tuple(moment_op.domain.shape) + (2, 12))
    current = PartialCurrentOperator(np.array([0.5, 0.2, 0.3]))
    emission = IsotropicEmissionOperator(1.0, 3)
    fine = prob.volume_measure
    trial, fine_fields = sn["mesh"].indicator_basis(), np.ones((fine.n_points, 2))
    prompt = OperatorPencil(homog.loss, homog.production).at(1.0)
    fold = sphere["fold"]
    ball = sphere["problem"]
    table = {
        "FunctionSpace.coordinates_of": (lambda: B.coordinates_of(phi), repr(B)),
        "FunctionSpace.leading_coordinates_of": (lambda: B.leading_coordinates_of(phi), repr(B)),
        "Member.__post_init__": (lambda: Member(values=phi, space=B), "Member"),
        "admit_array": (lambda: homog.production.apply(np.array([[0.6], [0.4]])), repr(homog.production.domain)),
        "AxisRetractionOperator.apply": (lambda: retraction.apply(np.ones(retraction.domain.shape)), repr(retraction.domain)),
        "FissionMomentOperator._admit": (lambda: moment_op.apply(moments), "FissionMomentOperator"),
        "DiagonalOperator.apply": (lambda: sn["engine"].apply(np.ones((sn["n_ord"], 2, 12))), "DiagonalOperator.apply"),
        "InnerProductFunctional.evaluate": (lambda: homog.production.production_rate.evaluate(np.array([[0.6], [0.4]])), "InnerProductFunctional.evaluate"),
        "DiscreteMeasure.integrate": (lambda: fine.integrate(np.ones(fine.n_points)), "DiscreteMeasure.integrate"),
        "SNSolver.compute_group_absorption_rate": (lambda: sn["solver"].compute_group_absorption_rate(phi), repr(prob.bulk_space)),
        "SNSolver.compute_group_production_rate": (lambda: sn["solver"].compute_group_production_rate(phi), repr(prob.bulk_space)),
        "Solution.reaction_rate_density": (lambda: sn_solution.reaction_rate_density(
            np.asarray(prob.mat_xs.absorption_cross_section)), repr(sn_solution.scalar_flux.space)),
        "DiagonalMetric.apply_block": (lambda: diagonal.apply_block(
            np.ones((2, 12)), start=1, block_shape=(12,), inverse=False), "DiagonalMetric.apply_block"),
        "DenseMetric.apply_block": (lambda: dense.apply_block(np.ones(2), start=0, block_shape=(2,), inverse=False), "DenseMetric.apply_block"),
        "FactoredMetric.apply": (lambda: factored.apply(phi), "FactoredMetric.apply"),
        "Basis.reconstruct": (lambda: lb.reconstruct(np.ones(3), gf.table), "LegendreBasis.reconstruct"),
        "Basis.analyze": (lambda: lb.analyze(np.ones(4), gf.table, q4.measure.weights), "LegendreBasis.analyze"),
        "FrameAnalysis.apply": (lambda: analysis.apply(np.ones((sn["n_ord"], 2, 12))), repr(analysis.domain)),
        "FrameReconstruction.apply": (lambda: reconstruction.apply(np.ones((1, 2, 12))), repr(reconstruction.domain)),
        "CrossGramInverse.apply": (lambda: cross_gram.apply(np.ones(3)), repr(cross_gram.domain)),
        "PartialCurrentOperator.apply": (lambda: current.apply(np.ones((3, 2))), "PartialCurrentOperator.apply (|Γ₊|)"),
        "IsotropicEmissionOperator.apply_transpose": (lambda: emission.apply_transpose(np.ones((3, 2))), "IsotropicEmissionOperator.apply_transpose (|Γ₋|)"),
        "MaterialXSField.project_through_bilinear": (lambda: prob.mat_xs.project_through_bilinear(
            trial, fine, phi=fine_fields, phi_star=fine_fields, rho=fine_fields), "MaterialXSField.project_through_bilinear (phi)"),
        "SourcePosing.__post_init__": (lambda: SourcePosing(prompt, np.ones((2, 1))), repr(prompt.codomain)),
        "RadialCharacteristicReconstruction.apply": (lambda: fold.apply(np.ones((1, ball.ng, ball.nx))), "RadialCharacteristicReconstruction.apply (n_moments, ng, nx)"),
    }
    return table[site]


def _find(obj, T, depth=0, seen=None):
    seen = seen if seen is not None else set()
    if depth > 7 or id(obj) in seen:
        return None
    seen.add(id(obj))
    if isinstance(obj, T):
        return obj
    for o in list(getattr(obj, "__dict__", {}).values()) + (list(obj) if isinstance(obj, (tuple, list)) else []):
        r = _find(o, T, depth + 1, seen)
        if r is not None:
            return r
    return None

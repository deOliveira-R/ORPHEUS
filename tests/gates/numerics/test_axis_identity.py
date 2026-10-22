r"""S7a gates: an axis is identified by what fixes its space, and a space by one key.

Spec: ``scratch/numerics_untangling/packages/S7_gate_spec.md`` (rows I*, O*, Q*, K*, R58, LB*),
written by the test-architect before the carve; each row was measured red for its named reason
on the pre-carve tree (``e5e47220``) and its mutation arm is named in the spec (§2). Issue #631;
rulings 58, 63, 64, 71, 73, 82, 84 of ``.claude/plans/numerics_untangling.md``.
"""
from __future__ import annotations

import dataclasses
import importlib
import inspect
import pkgutil

import numpy as np
import pytest

import orpheus
from orpheus.numerics import axis as axis_mod
from orpheus.numerics import frame_faces as ff
from orpheus.numerics.axis import Declared, Axis, BasisKind, Modal
from orpheus.values.content import content_digest
from orpheus.numerics.frame import GalerkinFrame, PetrovGalerkinFrame
from orpheus.geometry.manifold import Quotient
from orpheus.numerics.metric import DenseMetric
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.quadrature.rules_1d import gauss_legendre_on_polar_orbit
from orpheus.numerics.space import DualSpace, FunctionSpace
from orpheus.numerics.spaces.full_field_space import FullFieldSpace

pytestmark = pytest.mark.foundation



def _separated(a, b, what):
    """The separation legs (lessons M-identity (b)): ==, !=, both orders, and the container."""
    assert not (a == b) and not (b == a), f"{what}: the two compare EQUAL"
    assert a != b and b != a, f"{what}: != disagrees with =="
    assert len({a, b}) == 2, f"{what}: a set holds ONE of them (hash and == collide)"
    assert {a: 1}.get(b) is None, f"{what}: a dict keyed by one answers for the other"
    # Content-digest hashes: a collision on these specific pairs would itself be a defect.
    assert hash(a) != hash(b), f"{what}: equal hashes"


def _one(a, b, what):
    assert a == b and b == a, f"{what}: two spellings of one space compare unequal"
    assert hash(a) == hash(b), f"{what}: equal spaces with unequal hashes"
    assert len({a, b}) == 1, f"{what}: a set holds both"


# ═════════════════════════════════════════════════════════════════════════
# I — the #631 pairs, each beside the control pair that must stay equal
# ═════════════════════════════════════════════════════════════════════════

def test_I1_measure_spaces_with_different_weights_are_different_spaces():
    """Pair 1 (S7a). LS(4) and product(3, 8): 24 nodes each, different weights. TODAY the inline
    ``measure.space`` is nominal ``(L2[S^2], (24,))`` and they compare equal (census: 21 of 105 pairs)."""
    a, b = Quadrature.level_symmetric(4).measure, Quadrature.product(3, 8).measure
    assert a.n_points == b.n_points == 24
    assert not np.array_equal(np.sort(a.weights), np.sort(b.weights)), "PREMISE: the weights differ"
    _separated(FunctionSpace.of_measure(a, "angular"), FunctionSpace.of_measure(b, "angular"), "LS(4) vs product(3,8) measure spaces")


def test_I1b_equal_atoms_and_support_with_different_weights_are_different_spaces():
    """The weights-only discriminator (I1's pair also differs in its nodes, so a weights-blind identity
    keeps I1 green): one node set on S², the weights doubled."""
    from orpheus.numerics.measure import DiscreteMeasure
    m = Quadrature.level_symmetric(4).measure
    m2 = DiscreteMeasure(nodes=m.nodes, weights=2.0 * m.weights, support=m.support)
    _separated(FunctionSpace.of_measure(m, "angular"), FunctionSpace.of_measure(m2, "angular"), "one node set, weights x2")
    _separated(Axis.of_measure(m, "angular"), Axis.of_measure(m2, "angular"), "one node set, weights x2 (axes)")


def test_I1c_two_mints_of_one_measure_are_one_space():
    """Control for I1: a fresh mint of the same rule is the same space (content, not instance)."""
    _one(FunctionSpace.of_measure(Quadrature.level_symmetric(4).measure, "angular"),
         FunctionSpace.of_measure(Quadrature.level_symmetric(4).measure, "angular"), "two LS(4) mints")


def test_I1d_the_quadrature_axis_and_the_measure_axis_of_one_rule_are_one_axis():
    """Control: ``FromQuadrature(q)`` and ``FromMeasure(q.measure)`` under one label are equal (the arm is
    provenance; only support and atoms are content). Green today, must stay green."""
    q = Quadrature.product(4, 8)
    _one(Axis.of_quadrature(q), Axis.of_measure(q.measure, "angular"), "of_quadrature vs of_measure")
    _one(FunctionSpace.of_axes(Axis.of_quadrature(q)), FunctionSpace.of_measure(q.measure, "angular"),
         "the SN angular space vs the measure space of the same rule (design §6)")


def test_I2_angular_axes_with_different_directions_are_different_axes():
    """Pair 2 (S7a). LS(4) vs product(2, 12): equal weights (all 4π/24), different directions.
    TODAY ``Axis.of_quadrature`` of the two is equal (the production angular axis)."""
    ls, pr = Quadrature.level_symmetric(4), Quadrature.product(2, 12)
    assert np.array_equal(ls.measure.weights, pr.measure.weights), "PREMISE: equal weights"
    assert not np.array_equal(np.sort(ls.measure.nodes, axis=0), np.sort(pr.measure.nodes, axis=0)), \
        "PREMISE: different direction sets"
    _separated(Axis.of_quadrature(ls), Axis.of_quadrature(pr), "LS(4) vs product(2,12) axes")
    _separated(FunctionSpace.of_axes(Axis.of_quadrature(ls)), FunctionSpace.of_axes(Axis.of_quadrature(pr)),
               "LS(4) vs product(2,12) spaces")


def test_I2b_the_same_mu_rule_about_two_poles_gives_two_axes():
    """Ruling 71's refutation of "nodes and label suffice": GL(8) on S²/O2_x vs S²/O2_z share nodes and
    weights; only the support separates them. TODAY ``Axis.of_measure`` of the two is equal."""
    gx, gz = gauss_legendre_on_polar_orbit(8, "x"), gauss_legendre_on_polar_orbit(8, "z")
    assert np.array_equal(gx.nodes, gz.nodes) and np.array_equal(gx.weights, gz.weights), "PREMISE"
    _separated(Axis.of_measure(gx, "angular"), Axis.of_measure(gz, "angular"), "GL8 about x vs about z")
    _one(Axis.of_measure(gx, "angular"), Axis.of_measure(gauss_legendre_on_polar_orbit(8, "x"), "angular"),
         "control: two GL8-about-x mints")


def _dense_pair():
    unit = Modal(np.array([1.0, 0.0, 0.0]))
    ax = Axis("head", (3,), origin=Declared(unit))
    G = np.array([[2.0, 0.3, 0.0], [0.3, 1.0, 0.2], [0.0, 0.2, 3.0]])
    base = FunctionSpace.of_axes(ax)
    return base, G


def test_I5_a_dense_overlay_is_part_of_the_identity():
    """Pair 5 (S7a, design §1.5). Equal axes, different positioned dense forms. TODAY equal (the overlay
    lives outside the identity); the metric demonstrably differs (activation)."""
    base, G = _dense_pair()
    s1, s2 = base._with_overlay((DenseMetric(G),)), base._with_overlay((DenseMetric(2.0 * G),))
    x = np.ones(3)
    assert not np.array_equal(s1.apply_metric(x), s2.apply_metric(x)), "ACTIVATION: the metrics differ"
    _separated(s1, s2, "two dense overlays on equal axes")
    _one(s1, base._with_overlay((DenseMetric(G.copy()),)), "control: one dense form, two mints")
    _separated(s1, base, "a dense overlay vs the bare head")


def test_I5c_a_dense_form_is_its_matrix_not_its_inverse_face():
    """The inverse face is derived from the matrix, so it is not content: the production dense head
    (whose frame supplies the exact inverse face) and the same axes re-dressed with ``DenseMetric(G)``
    (which pseudo-inverts) are ONE space. [M] before the fix (the S7 archivist, 2026-10-22): unequal,
    with array_equal matrices, the digests differing in the inverse face's round-off."""
    f = ff.angular_frame(Quadrature.folded_product(n_mu=2, n_phi=4), 2)
    bs = ff.basis_space(f)
    form = bs._overlay_forms()[0]
    assert isinstance(form, DenseMetric)
    G = np.array(form.matrix)
    rebuilt = FunctionSpace.of_axes(*bs.axes)._with_overlay((DenseMetric(G.copy()),))
    assert not np.array_equal(np.asarray(form.inverse_matrix), np.asarray(rebuilt._overlay_forms()[0].inverse_matrix)), \
        "PREMISE: the two inverse faces differ in round-off"
    _one(bs, rebuilt, "one dense form, its inverse supplied vs derived")


def test_I5b_the_production_dense_head_carries_its_overlay_in_its_identity():
    """Pair 5 on a production head: ``basis_space`` of a DENSE-Gram frame (folded_product(2, 4), L = 2),
    re-dressed with twice its form, is a different space; re-dressed with its own form, the same."""
    f = ff.angular_frame(Quadrature.folded_product(n_mu=2, n_phi=4), 2)
    bs = ff.basis_space(f)
    forms = bs._overlay_forms()
    dense = [k for k, form in enumerate(forms) if isinstance(form, DenseMetric)]
    assert dense == [0], f"PREMISE: the head carries one dense form on axis 0, got {forms}"
    form = forms[0]
    assert isinstance(form, DenseMetric)
    G = form.matrix
    assert bs.axes is not None
    twice = FunctionSpace.of_axes(*bs.axes)._with_overlay((DenseMetric(2.0 * np.asarray(G)),))
    _separated(bs, twice, "the dressed head vs the head with twice its Gram form")


def test_I6_F7_a_composite_dual_is_not_a_composite_of_named_spaces():
    """Pair 6 (S7a, F7 of the S6a review). TODAY ``from_blocks(a, b).dual() == from_blocks(a*, b*)`` with
    equal hashes: the composite folds member NAMES, which carry no covariance."""
    a, b = FunctionSpace("a", (2,)), FunctionSpace("b", (3,))
    lhs = FullFieldSpace.from_blocks(a, b).dual()
    rhs = FullFieldSpace.from_blocks(FunctionSpace("a*", (2,)), FunctionSpace("b*", (3,)))
    _separated(lhs, rhs, "F7: composite dual vs composite of starred names")
    _one(lhs, FullFieldSpace.from_blocks(FunctionSpace("a", (2,)), FunctionSpace("b", (3,))).dual(),
         "control: two composite duals of equal primals")
    _one(DualSpace.of(a), DualSpace.of(FunctionSpace("a", (2,))), "control: DualSpace of equal primals")


# ═════════════════════════════════════════════════════════════════════════
# LB — the measure space's label (orchestrator's ruling inside 63)
# ═════════════════════════════════════════════════════════════════════════

def test_LB_ii_the_frame_measure_space_is_the_sn_angular_space():
    q = Quadrature.level_symmetric(6)
    _one(ff.measure_space(ff.angular_frame(q, 2)), FunctionSpace.of_axes(Axis.of_quadrature(q)),
         "measure_space(angular_frame(q, L)) vs of_axes(Axis.of_quadrature(q))")


def test_LB_iii_a_spatial_frame_measure_axis_is_labelled_spatial():
    from orpheus.numerics.basis.indicator_basis import IndicatorBasis
    from orpheus.geometry.manifold import RealSpace
    from orpheus.numerics.measure import DiscreteMeasure
    basis = IndicatorBasis((np.array([0.0, 1.0, 2.0]),), RealSpace(1))
    m = DiscreteMeasure(nodes=np.array([0.25, 0.75, 1.25, 1.75]), weights=np.full(4, 0.5), support=RealSpace(1))
    vs = ff.measure_space(GalerkinFrame(basis, m))
    assert vs.axes is not None, "the measure space is axis-built (S7a, ruling 63)"
    assert [ax.label for ax in vs.axes] == ["spatial"]


def _gauge_production():
    """The production gauge on the reflective 3x4 box (test_loss_kernel_gauge's fixture)."""
    from orpheus.sn.operators.loss_kernel_gauge import LossKernelGauge
    from tests.gates.sn.operators.test_loss_kernel_gauge import _R, _mesh
    return LossKernelGauge.for_mesh(_mesh((3, 4), [(_R, _R)] * 2))


def test_LB_i_iv_the_gauge_gather_codomain_is_the_projector_domain():
    """(i) one space by construction; (iv) labelled with the index set's own name. TODAY the codomain
    is the axes-less inline space, so (iv) reds on ``axes is None``."""
    blocks = _gauge_production()._blocks if hasattr(_gauge_production(), "_blocks") else _gauge_production().blocks
    assert blocks, "PREMISE: the gauge has blocks"
    for blk in blocks:
        assert blk.gather.codomain is not None and blk.gather.codomain == blk.projector.domain, "(i) gather codomain != projector domain"
        codomain = blk.gather.codomain
        assert codomain is not None
        axes = codomain.axes
        assert axes is not None and len(axes) == 1, "(iv) the index-set measure space is axis-built"
        assert axes[0].label == blk.basis.domain.name, f"(iv) label {axes[0].label!r} != {blk.basis.domain.name!r}"


def test_LB_v_axis_label_refuses_a_support_with_no_role():
    """(v) ``DiscreteMeasure.axis_label`` refuses a support that is none of sphere/orbit space,
    RealSpace, EnergyGroups, IndexSet, naming the support. A bare Interval with no O(3) tag is the
    witness (a rule on ``Interval(-1, 1)`` spelled as a spatial interval)."""
    from orpheus.geometry.manifold import Interval
    from orpheus.numerics.measure import DiscreteMeasure
    m = DiscreteMeasure(nodes=np.array([-0.5, 0.5]), weights=np.array([1.0, 1.0]), support=Interval(-1.0, 1.0))
    if m.invariance_group is not None:
        pytest.fail("PREMISE: the witness carries no invariance tag (else phase answers angular)")
    with pytest.raises(NotImplementedError, match=r"Interval|\[-1"):
        _ = m.axis_label


# ═════════════════════════════════════════════════════════════════════════
# O — the Origin closed sum (ruling 64 with Declared, ruling 82)
# ═════════════════════════════════════════════════════════════════════════

def test_O1_an_axis_takes_no_kind_beside_its_origin():
    """``Axis(kind=NODAL, generator=basis)`` is unspellable: the signature has ``origin`` and no
    ``kind`` / ``generator``, and passing ``kind=`` is a TypeError. TODAY ``kind`` is a parameter."""
    params = inspect.signature(Axis).parameters
    assert "origin" in params and "kind" not in params and "generator" not in params, list(params)
    with pytest.raises(TypeError):
        Axis("x", (3,), kind=BasisKind.NODAL)  # type: ignore[call-arg]
    assert not hasattr(Axis, "generator_as"), "generator_as retires with the union (design §2)"


def test_O2_each_arm_implies_its_kind_and_exposes_its_content():
    D, FM, FQ, FF = (getattr(axis_mod, n) for n in ("Declared", "FromMeasure", "FromQuadrature", "FromFrame"))
    q = Quadrature.product(4, 8)
    assert FQ(q).kind is BasisKind.NODAL and FM(q.measure).kind is BasisKind.NODAL
    assert FQ(q).support == q.measure.support and np.array_equal(FQ(q).nodes, q.measure.nodes)
    assert D(BasisKind.NODAL).support is None and D(BasisKind.NODAL).nodes is None
    f = ff.angular_frame(q, 2)
    # FromFrame's unit is Phi^T w, read independently off the frame's trial table and the measure.
    unit = np.einsum("n...,n->...", np.asarray(f.table), np.asarray(q.measure.weights))  # Phi^T w
    k = FF(f).kind
    assert isinstance(k, Modal)
    # 8 eps x the reduction depth (32 nodes): [M] today's dressed unit sits 4.4e-16 from this reference
    np.testing.assert_allclose(k.unit.ravel(), unit.ravel(), rtol=0, atol=1e-13)


def test_O3_the_direction_reads_match_FromQuadrature_and_nothing_else():
    """AST: the set of ``case FromQuadrature(...)`` patterns is the declared consumer set (re-poses
    test_axis_generator G5_6b). Positive control: the set is non-empty and holds for_space."""
    import ast
    import pathlib
    root = pathlib.Path(orpheus.__file__).parent
    found = set()
    for f in root.rglob("*.py"):
        tree = ast.parse(f.read_text())
        for fn in ast.walk(tree):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for n in ast.walk(fn):
                if isinstance(n, ast.MatchClass) and getattr(n.cls, "id", getattr(n.cls, "attr", "")) == "FromQuadrature":
                    found.add((str(f.relative_to(root.parent)), fn.name))
    declared = {
        ("orpheus/sn/mesh/reduced_operator.py", "streaming_terms"),
        ("orpheus/sn/angular/closure.py", "__init__"),
        ("orpheus/transport/frames/harmonic_frame.py", "for_space"),
        ("orpheus/transport/frames/harmonic_frame.py", "moment_space_on"),
    }
    assert found, "CENSUS BROKEN: no `case FromQuadrature(` anywhere"
    assert found == declared, f"extra {sorted(found - declared)}; missing {sorted(declared - found)}"


@pytest.mark.parametrize("arm", ["declared", "from_measure"])
def test_O4_for_space_refuses_a_non_quadrature_axis_naming_itself(arm):
    """The per-consumer refusal, widened to every non-FromQuadrature arm (improves
    test_frame_mints_refuse_a_moment_space.py, which builds the generator-less arm only)."""
    from orpheus.transport.frames import HarmonicFrame
    q = Quadrature.gauss_legendre(4)
    ax = (Axis("angular", (q.n_ordinates,), np.asarray(q.weights, float), origin=Declared(BasisKind.NODAL))
          if arm == "declared" else Axis.of_measure(q.measure, "angular"))
    with pytest.raises(ValueError, match="HarmonicFrame.for_space needs the generating Quadrature"):
        HarmonicFrame.for_space(FunctionSpace.of_axes(ax), 1)


def test_O5_the_energy_axis_is_Declared_nodal_and_takes_no_origin():
    """Ruling 82 + Q2: EnergyAxis fixes Declared(NODAL) and refuses another origin."""
    from orpheus.numerics.axis import EnergyAxis
    ax = EnergyAxis.synthetic(2)
    assert isinstance(ax.origin, axis_mod.Declared) and ax.origin.kind is BasisKind.NODAL
    with pytest.raises(TypeError):  # Q2 ruling: EnergyAxis takes no origin= argument at all
        EnergyAxis("energy", (2,), origin=axis_mod.Declared(BasisKind.NODAL))  # type: ignore[call-arg]


# ═════════════════════════════════════════════════════════════════════════
# Q — the support is content-encodable (Quotient)
# ═════════════════════════════════════════════════════════════════════════

def test_Q1_a_quotient_has_a_content_digest_and_its_sympy_fields_are_not_content():
    sx, sz = gauss_legendre_on_polar_orbit(8, "x").support, gauss_legendre_on_polar_orbit(8, "z").support
    assert isinstance(sx, Quotient)
    assert content_digest(sx) != content_digest(sz), "S²/O2_x and S²/O2_z must differ as content"
    assert content_digest(sx) == content_digest(gauss_legendre_on_polar_orbit(8, "x").support)
    compare = {f.name: f.compare for f in dataclasses.fields(Quotient)}
    for derived in ("generators", "syzygy", "gram", "det_gram"):
        assert compare[derived] is False, f"Quotient.{derived} is derived from (base, by) and must not be compared"
    assert compare["derived_by"] is True and compare["base"] is True and compare["by"] is True


# ═════════════════════════════════════════════════════════════════════════
# K — one identity key (design §5)
# ═════════════════════════════════════════════════════════════════════════

def _all_space_classes():
    for m in pkgutil.walk_packages(orpheus.__path__, "orpheus."):
        if ".derivations" in m.name or "_traced_memo_boot" in m.name:
            continue
        importlib.import_module(m.name)

    def subs(c):
        for s in c.__subclasses__():
            yield s
            yield from subs(s)
    return sorted(set(subs(FunctionSpace)), key=lambda c: c.__qualname__)


def test_K1_equality_and_hash_are_defined_once_on_the_base():
    """TODAY 12 of 12 subclasses define their own __eq__/__hash__ (8 delegate, 4 dataclass-generated:
    AngularTraceSpace, AngularFaceTraceSpace, the two radial-characteristic spaces)."""
    # the production family only: a test module may define a toy subclass
    classes = [c for c in _all_space_classes() if c.__module__.startswith("orpheus.")]
    names = {c.__qualname__ for c in classes}
    assert len(classes) >= 10 and {"AngularTraceSpace", "DualSpace", "FullFieldSpace"} <= names, \
        f"CENSUS BROKEN: {sorted(names)}"
    own = sorted(c.__qualname__ for c in classes if "__eq__" in vars(c) or "__hash__" in vars(c))
    assert own == [], f"{len(own)} of {len(classes)} space classes define their own __eq__/__hash__: {own}"
    assert "identity_key" in vars(FunctionSpace)


def _population():
    from tests.gates.sn.regression import _generate_snapshots as G
    from orpheus.sn.solver import _as_problem
    out = []
    for c in (G._slab_3region("2g", 6), G._cylinder_3region("2g", 4, "folded_4x6")):
        p = _as_problem(c["mesh"], c["quadrature"], c["materials"], scattering_order=1)
        ffs = p.full_field_space
        out += [ffs, ffs.dual(), ffs.interior_space, ffs.interior_space.dual(), getattr(ffs, "trace_space", None), p.system.space,
                p.moment_space(1), p.angular_bulk_space]
        rc = getattr(p, "radial_characteristic_field_space", None)
        if rc is not None:
            out += [rc, rc.interior_space]
    q = Quadrature.folded_product(n_mu=2, n_phi=4)
    out += [ff.basis_space(ff.angular_frame(q, 2)), FunctionSpace("named", (3,))]
    return [s for s in out if isinstance(s, FunctionSpace)]


def test_K2_eq_and_hash_derive_from_the_identity_key():
    pop = _population()
    kinds = {type(s).__name__ for s in pop}
    assert len(pop) >= 12 and {"FullFieldSpace", "CoupledSpace", "TensorProductSpace", "DualSpace"} <= kinds, kinds
    for s in pop:
        assert hash(s) == hash(s.identity_key()), f"{type(s).__name__}: hash is not hash(identity_key)"
    for a in pop:
        for b in pop:
            assert (a == b) == (a.identity_key() == b.identity_key()), f"{type(a).__name__} vs {type(b).__name__}"


# ═════════════════════════════════════════════════════════════════════════
# R58 — the Galerkin test-space predicate is `frame.test is frame.basis`
# ═════════════════════════════════════════════════════════════════════════

def test_R58_a_petrov_galerkin_frame_whose_test_is_its_basis_gets_the_dressed_test_space():
    from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
    b, m = SphericalHarmonicBasis(L=2), Quadrature.product(4, 8).measure
    pg = PetrovGalerkinFrame(b, m, b)
    _one(ff.test_space(pg), ff.basis_space(pg), "test_space(PG(b, m, b)) vs basis_space")
    assert not (ff.test_space(pg) == FunctionSpace.of_basis(b)), "the pre-S7 undressed answer"
    g = GalerkinFrame(b, m)
    _one(ff.test_space(g), ff.basis_space(g), "control: the Galerkin frame")

"""Foundation tests for ``orpheus.numerics.quadrature.registry``.

Verifies the five-stage selection chain (domain → symmetry → V →
structural → minimum points) end-to-end on the four canonical geometry
tags shipped today (``slab``, ``sphere``, ``cylinder``,
``cartesian2d``), plus the explainability log and the negative paths
(no-rule-fits, unknown flags, unknown geometry).

Stages 0 and 1 are the two halves of :class:`AngularSymmetry`: the
continuous isotropy the dimensional reduction SPENDS (which fixes the
angular domain) and the discrete residual still OWED (which the rule
must realize as an ordinate permutation). Before 2026-08-02 the table
recorded only the spent half and compared it against a rule's declared
group — a gate no discrete azimuthal rule can pass, which is why it
only ever passed on a false declaration.

The selection chain is a *software* invariant — it does not verify a
physics equation in a Sphinx theory page — so the tests are tagged
``foundation`` rather than carrying ``verifies(...)`` markers. The
underlying primitive contracts (Gauss-Legendre exactness, Lebedev
:math:`O_h`-invariance, level-symmetric :math:`S_N` polar-level
structure) carry their own L1 tests in
``tests/numerics/test_rules_*.py``.
"""

from __future__ import annotations

import pytest

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import (
    GEOMETRY_ANGULAR_SYMMETRY,
    AngularSymmetry,
    QuadratureSelectionError,
    QuadratureSpec,
    SelectionLog,
    quadrature_registry,
    select_quadrature,
)
from orpheus.numerics.quadrature.registry import (
    _gl1d_invert,
    _lebedev_invert,
    _ls_sn_invert,
    _product_invert,
)
from orpheus.numerics.symmetry import SubgroupOfO3


# ---------------------------------------------------------------------------
# Registry sanity
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_registry_population() -> None:
    """The registry must contain the four shipped rules."""
    names = {spec.name for spec in quadrature_registry}
    assert names == {
        "GaussLegendre1D",
        "LebedevSphere",
        "LevelSymmetricSN",
        "ProductQuadrature",
    }


@pytest.mark.foundation
def test_registry_specs_are_well_formed() -> None:
    """Every spec must have callable inversion / node-count.

    A spec deliberately carries NO invariance group: a rule's symmetry
    is parameter-dependent, so a parameter-free field cannot state it
    truthfully. The selector computes it from the built measure.
    """
    for spec in quadrature_registry:
        assert isinstance(spec, QuadratureSpec)
        assert callable(spec.factory)
        assert callable(spec.degree_of_exactness_for)
        assert callable(spec.expected_node_count)
        assert not hasattr(spec, "invariance_group")
        assert isinstance(spec.parameters, dict)
        assert all(isinstance(t, type) for t in spec.parameters.values())
        # Structural-flags accessor must surface the four flags.
        flags = spec.structural_flags()
        assert set(flags) == {
            "positive_weights",
            "axis_aligned",
            "level_structured",
            "half_range_clean",
        }


@pytest.mark.foundation
def test_geometry_angular_symmetry_table() -> None:
    """Each geometry declares both halves of its angular symmetry.

    The SPENT continuous half fixes the angular domain; the OWED
    discrete half is what a quadrature must realize as an ordinate
    permutation. Recording only the spent half (as this table did
    until 2026-08-02) makes the selection gate unsatisfiable — no
    finite point set on S^2 is SO(2)-closed.
    """
    assert GEOMETRY_ANGULAR_SYMMETRY == {
        "slab": AngularSymmetry(
            continuous_isotropy=SubgroupOfO3.SO2,
            discrete_residual=SubgroupOfO3.Mirror("x"),
        ),
        "sphere": AngularSymmetry(
            continuous_isotropy=SubgroupOfO3.SO2,
            discrete_residual=SubgroupOfO3.Mirror("x"),
        ),
        "cylinder": AngularSymmetry(
            continuous_isotropy=SubgroupOfO3.Trivial,
            discrete_residual=SubgroupOfO3.Dnh(2),
        ),
        "cartesian2d": AngularSymmetry(
            continuous_isotropy=SubgroupOfO3.Trivial,
            discrete_residual=SubgroupOfO3.Dnh(2),
        ),
    }


# ---------------------------------------------------------------------------
# Inversion helpers
# ---------------------------------------------------------------------------
#
# The four ``_*_invert`` helpers must return the smallest parameter set
# that satisfies the target-degree constraint. Verify each independently
# so the upstream selector tests can rely on the contract.


@pytest.mark.foundation
@pytest.mark.parametrize(
    "target,expected_n",
    [
        (0, 2),    # ceil(1/2)=1, rounded up to even -> 2 (deg 3)
        (1, 2),
        (3, 2),
        (4, 4),    # ceil(5/2)=3 -> 4 (deg 7)
        (15, 8),   # ceil(16/2)=8 (even already) -> deg 15
        (16, 10),  # ceil(17/2)=9 -> 10 (even) -> deg 19
    ],
)
def test_gl1d_invert(target: int, expected_n: int) -> None:
    params = _gl1d_invert(target)
    assert params == {"n": expected_n}


@pytest.mark.foundation
def test_gl1d_invert_negative_returns_none() -> None:
    assert _gl1d_invert(-1) is None


@pytest.mark.foundation
@pytest.mark.parametrize(
    "target,expected_order",
    [
        (3, 3),
        (4, 5),
        (5, 5),
        (15, 15),
        (32, 35),  # 31 < 32 < 35; pick 35
        (47, 47),
    ],
)
def test_lebedev_invert(target: int, expected_order: int) -> None:
    params = _lebedev_invert(target)
    assert params == {"order": expected_order}


@pytest.mark.foundation
def test_lebedev_invert_too_high_returns_none() -> None:
    """target_degree above the highest tabulated Lebedev order has no
    solution — the inverter MUST return ``None`` (not raise)."""
    assert _lebedev_invert(48) is None
    assert _lebedev_invert(100) is None


@pytest.mark.foundation
@pytest.mark.parametrize(
    "target,expected_N",
    [
        (1, 2),    # min even N=2 (deg 1)
        (3, 4),    # N=4 (deg 3)
        (4, 6),    # need N>=5; round up to even -> 6 (deg 5)
        (15, 16),  # N=16 (deg 15)
    ],
)
def test_ls_sn_invert(target: int, expected_N: int) -> None:
    params = _ls_sn_invert(target)
    assert params == {"sn_order": expected_N}


@pytest.mark.foundation
def test_product_invert() -> None:
    params = _product_invert(4)
    # n_mu = ceil(5/2) = 3; n_phi = 5
    assert params == {"n_mu": 3, "n_phi": 5}

    params = _product_invert(15)
    # n_mu = ceil(16/2) = 8; n_phi = 16
    assert params == {"n_mu": 8, "n_phi": 16}


# ---------------------------------------------------------------------------
# Selector — happy paths
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_select_slab_returns_gauss_legendre() -> None:
    """A slab discretises S^2/SO(2) = [-1,1], so GL1D is the only rule
    on the right domain — the three S^2 rules fail stage 0."""
    measure, log = select_quadrature("slab", target_degree=15)

    assert isinstance(measure, DiscreteMeasure)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"
    assert log.chosen_parameters == {"n": 8}
    assert measure.n_points == 8
    assert measure.degree_of_exactness == 15
    assert measure.invariance_group == SubgroupOfO3.SO2
    # Every S^2 rule is rejected at the DOMAIN stage.
    rejected_names = {name for name, _ in log.rejected}
    assert rejected_names == {
        "LebedevSphere",
        "LevelSymmetricSN",
        "ProductQuadrature",
    }


@pytest.mark.foundation
def test_select_sphere_returns_gauss_legendre() -> None:
    """1-D radial spherical SN reduces to GL on :math:`\\mu_r` — the
    same spent/owed split as the slab, so the same rule wins.

    The spent half is not free here: its fiber action reappears in the
    curvilinear sweep as the angular-redistribution alpha term."""
    measure, log = select_quadrature("sphere", target_degree=15)

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"
    assert log.chosen_parameters == {"n": 8}
    assert measure.degree_of_exactness == 15


@pytest.mark.foundation
def test_select_cylinder_with_level_structured_returns_product() -> None:
    """Cylindrical SN sweep needs polar-level structure, and the
    product rule supplies it — *when its azimuthal count is even*.

    ``target_degree=5`` inverts to ``n_phi = 6``. D_2h needs mirror
    planes at 0 and pi/2, and D_{n h} carries them at k*pi/n, so
    ``n_phi`` even is exactly the admissibility condition.
    """
    measure, log = select_quadrature(
        "cylinder", target_degree=5, level_structured=True
    )

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "ProductQuadrature"
    assert log.chosen_parameters == {"n_mu": 3, "n_phi": 6}
    # GL1D is rejected at the DOMAIN stage now, not the structural one:
    # a mu-marginal cannot carry a cylinder's two angular DOF at all,
    # which is a stronger and earlier objection than "no levels".
    gl_reason = next(
        reason for name, reason in log.rejected
        if name == "GaussLegendre1D"
    )
    assert "domain mismatch" in gl_reason


@pytest.mark.foundation
def test_select_cylinder_rejects_odd_azimuthal_product_rule() -> None:
    """ERR-042, structurally: an odd-``n_phi`` product rule is NOT
    invariant under the cylinder's owed D_2h, and the selector must
    refuse it rather than hand back an asymmetric rule.

    ``target_degree=4`` inverts to ``n_phi = 5``. The sigma_x mirror
    sends phi -> 180 - phi, mapping the 0-degree node to 180 degrees,
    which is not a node of a 5-fold grid. Before 2026-08-02 this rule
    was ADMITTED here, on the strength of a false ``SO(2)`` tag.

    The selector degrades gracefully: it falls back to a rule that
    genuinely carries the symmetry, rather than failing.
    """
    _, log = select_quadrature(
        "cylinder", target_degree=4, level_structured=True
    )

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "LevelSymmetricSN"

    product_reason = next(
        reason for name, reason in log.rejected
        if name == "ProductQuadrature"
    )
    assert "symmetry mismatch" in product_reason
    assert "D_2h" in product_reason

    # And the discriminator is parity, not degree: one step up in
    # target degree makes n_phi even and the same rule admissible.
    _, log_even = select_quadrature(
        "cylinder", target_degree=5, level_structured=True
    )
    assert log_even.chosen_spec is not None
    assert log_even.chosen_spec.name == "ProductQuadrature"


@pytest.mark.foundation
def test_select_cartesian2d_prefers_lebedev_over_ls_sn() -> None:
    """2-D Cartesian owes :math:`D_{2h}`. Lebedev, LS_N and the
    even-``n_phi`` product rule all satisfy it; Lebedev wins on
    minimum-points (14 nodes for order 5 vs 48 for LS_6, 18 for the
    product rule).

    The :class:`SelectionLog` must record this preference: the losers
    on cost are NOT in the rejected list.

    Only GL1D is rejected, and on DOMAIN — a mu-marginal cannot serve
    a geometry that keeps both angular degrees of freedom. Before
    2026-08-02 the product rule was rejected here too, because the
    table over-claimed :math:`O_h`: that demands the x<->z exchange,
    never a symmetry of a z-uniform problem.
    """
    measure, log = select_quadrature("cartesian2d", target_degree=5)

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "LebedevSphere"
    assert log.chosen_parameters == {"order": 5}
    assert measure.n_points == 14

    rejected_names = {name for name, _ in log.rejected}
    assert rejected_names == {"GaussLegendre1D"}
    gl_reason = next(
        reason for name, reason in log.rejected
        if name == "GaussLegendre1D"
    )
    assert "domain mismatch" in gl_reason


# ---------------------------------------------------------------------------
# Selector — explainability log
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_log_records_geometry_metadata() -> None:
    """The :class:`SelectionLog` must carry the geometry, group, target
    degree, and requested flags from the call so consumers can render
    the decision."""
    measure, log = select_quadrature(
        "cylinder", target_degree=4, level_structured=True
    )
    assert log.geometry == "cylinder"
    assert log.angular_symmetry == AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.Trivial,
        discrete_residual=SubgroupOfO3.Dnh(2),
    )
    # Both halves are reconstructible from the log, so a reader can
    # replay stages 0 and 1 without re-deriving the decomposition.
    assert log.angular_symmetry.support == "S^2"
    assert log.target_degree == 4
    assert log.requested_flags == {"level_structured": True}


@pytest.mark.foundation
def test_log_summary_string_is_informative() -> None:
    """``log.summary()`` is the one-line debug string used by SN
    diagnostic output — must include geometry, target degree, chosen
    rule, and parameters."""
    _, log = select_quadrature("slab", target_degree=15)
    summary = log.summary()
    assert "slab" in summary
    assert "15" in summary
    assert "GaussLegendre1D" in summary
    assert "n" in summary  # parameter key


@pytest.mark.foundation
def test_log_rejected_list_carries_reasons() -> None:
    """Each rejection must come with a reason string identifying the
    failing stage (domain / symmetry / V / structural)."""
    _, log = select_quadrature("slab", target_degree=15)
    assert log.rejected, "expected the S^2 rules to be rejected for a slab"
    for name, reason in log.rejected:
        assert isinstance(name, str)
        assert isinstance(reason, str)
        # Every slab rejection is at the DOMAIN stage: Lebedev, LS_N
        # and the product rule all live on S^2, while a slab
        # discretises the quotient S^2/SO(2) = [-1,1].
        assert "domain mismatch" in reason
        assert "[-1,1]" in reason


@pytest.mark.foundation
def test_dont_care_flags_dont_filter() -> None:
    """Passing a flag with value ``False`` must NOT filter rules — only
    ``True`` constrains."""
    # No flag request: any structural set is accepted.
    measure_a, _ = select_quadrature("slab", target_degree=15)
    # axis_aligned=False (don't care): same answer.
    measure_b, _ = select_quadrature(
        "slab", target_degree=15, axis_aligned=False
    )
    assert measure_a.n_points == measure_b.n_points


# ---------------------------------------------------------------------------
# Selector — negative paths
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_unknown_geometry_raises() -> None:
    """Unknown geometry tags raise :class:`KeyError` immediately —
    this is a programming error, not a selection failure."""
    with pytest.raises(KeyError, match="hexagonal"):
        select_quadrature("hexagonal", target_degree=5)


@pytest.mark.foundation
def test_unknown_structural_flag_raises() -> None:
    """A typo in a structural flag must fail loudly, not be silently
    ignored."""
    with pytest.raises(KeyError, match="level_strucutred"):
        select_quadrature(
            "cylinder", target_degree=4, level_strucutred=True
        )


@pytest.mark.foundation
def test_no_rule_fits_raises_with_log() -> None:
    """Asking for a target degree no rule can supply must raise
    :class:`QuadratureSelectionError` carrying a populated
    :class:`SelectionLog`.

    Forcing this state with the full registry is awkward because LS_N
    inversion accepts arbitrarily large :math:`N`. Instead we override
    the registry to contain only Lebedev (which tops out at order 47)
    and ask for ``target_degree=100``, which forces a V-stage rejection.
    """
    only_lebedev = [
        s for s in quadrature_registry if s.name == "LebedevSphere"
    ]
    with pytest.raises(QuadratureSelectionError) as excinfo:
        select_quadrature(
            "cartesian2d", target_degree=100, registry=only_lebedev
        )

    assert excinfo.value.log is not None
    assert excinfo.value.log.chosen_spec is None
    rejected_names = {name for name, _ in excinfo.value.log.rejected}
    assert "LebedevSphere" in rejected_names
    leb_reason = next(
        reason for name, reason in excinfo.value.log.rejected
        if name == "LebedevSphere"
    )
    assert "V mismatch" in leb_reason or "target_degree" in leb_reason


@pytest.mark.foundation
def test_slab_cannot_supply_polar_levels_and_says_so() -> None:
    """Asking a slab for ``level_structured=True`` is a category error,
    and the log must be able to say why.

    A slab discretises the quotient S^2/SO(2) = [-1,1]. Its angular
    variable IS the polar axis, so there are no per-mu polar *levels*
    to expose — a level structure is a fibration of S^2 over mu, and
    the slab has already taken that quotient.

    Every rule is therefore rejected, and the two rejection reasons
    are the two different objections: the S^2 rules fail on DOMAIN,
    and the one domain-admissible rule (GL1D) fails on the structural
    flag. Before 2026-08-02 this returned ProductQuadrature — an S^2
    rule handed to a 1-D solver, admitted only because the geometry
    table compared spent-group to spent-group.
    """
    with pytest.raises(QuadratureSelectionError) as excinfo:
        select_quadrature(
            "slab",
            target_degree=4,
            level_structured=True,
            half_range_clean=True,
        )

    log = excinfo.value.log
    assert log.chosen_spec is None
    reasons = dict(log.rejected)
    assert set(reasons) == {
        "GaussLegendre1D",
        "LebedevSphere",
        "LevelSymmetricSN",
        "ProductQuadrature",
    }
    assert "structural mismatch" in reasons["GaussLegendre1D"]
    assert "level_structured" in reasons["GaussLegendre1D"]
    for name in ("LebedevSphere", "LevelSymmetricSN", "ProductQuadrature"):
        assert "domain mismatch" in reasons[name]

    # GL1D rejected at structural (no level structure).
    gl_reason = next(
        reason for name, reason in log.rejected
        if name == "GaussLegendre1D"
    )
    assert "structural" in gl_reason.lower()


@pytest.mark.foundation
def test_truly_incompatible_flags_raises() -> None:
    """Cylinder + ``axis_aligned=True`` + ``level_structured=True``:
    no rule on S^2 carries both flags, so the request has no candidate
    and must raise with a populated log."""
    with pytest.raises(QuadratureSelectionError) as excinfo:
        select_quadrature(
            "cylinder",
            target_degree=4,
            level_structured=True,
            axis_aligned=True,
        )
    log = excinfo.value.log
    assert log.chosen_spec is None
    # Every spec must appear in the rejected list.
    rejected_names = {name for name, _ in log.rejected}
    assert rejected_names == {
        "GaussLegendre1D",
        "LebedevSphere",
        "LevelSymmetricSN",
        "ProductQuadrature",
    }


# ---------------------------------------------------------------------------
# Selector — registry override
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_registry_override_failure_path() -> None:
    """Registry override that excludes the only matching rule must
    raise :class:`QuadratureSelectionError`."""
    only_gl = [s for s in quadrature_registry if s.name == "GaussLegendre1D"]
    with pytest.raises(QuadratureSelectionError):
        select_quadrature("cartesian2d", target_degree=5, registry=only_gl)


@pytest.mark.foundation
def test_registry_override_success_path() -> None:
    """Registry override containing only one matching rule must succeed
    and return that rule."""
    only_lebedev = [
        s for s in quadrature_registry if s.name == "LebedevSphere"
    ]
    measure, log = select_quadrature(
        "cartesian2d", target_degree=5, registry=only_lebedev
    )
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "LebedevSphere"
    assert measure.invariance_group == SubgroupOfO3.OctahedralOh


# ---------------------------------------------------------------------------
# AngularSymmetry — the type's own laws
# ---------------------------------------------------------------------------
#
# The two halves are not interchangeable, and the failure they prevent is
# specific: recording the SPENT half where the OWED one belongs produces a
# gate that is unsatisfiable by construction. These gates pin both the
# derivation (support from the spent group) and the discrimination (which
# rules the owed group actually selects).


@pytest.mark.foundation
def test_support_is_derived_from_the_spent_group_not_declared() -> None:
    """The angular domain IS the quotient S^2/G^0.

    Deriving it rather than storing a second column is what keeps the
    domain and the spent group from drifting apart — they are one fact,
    so there is no state in which they can disagree.
    """
    assert AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.SO2,
        discrete_residual=SubgroupOfO3.Mirror("z"),
    ).support == "[-1,1]"

    assert AngularSymmetry(
        continuous_isotropy=SubgroupOfO3.Trivial,
        discrete_residual=SubgroupOfO3.Dnh(2),
    ).support == "S^2"

    # An unmapped quotient must refuse loudly rather than guess: a wrong
    # domain silently admits a rule of the wrong dimensionality.
    with pytest.raises(NotImplementedError, match="S\\^2/"):
        _ = AngularSymmetry(
            continuous_isotropy=SubgroupOfO3.OctahedralOh,
            discrete_residual=SubgroupOfO3.Mirror("z"),
        ).support


@pytest.mark.foundation
def test_the_two_stages_are_independent_and_both_load_bearing() -> None:
    """Neither conjunct subsumes the other — drop either and the gate
    admits something wrong.

    Measured, not assumed (the T18 lesson: before folding predicates
    onto one primitive, measure what each actually selects):

    * symmetry alone admits a Lebedev rule for a slab, because an
      O_h-invariant rule certainly satisfies Z_2;
    * domain alone admits an odd-``n_phi`` product rule for a cylinder,
      because it does live on S^2.
    """
    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    cylinder = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]

    lebedev = next(
        s for s in quadrature_registry if s.name == "LebedevSphere"
    ).build({"order": 5})
    odd_product = next(
        s for s in quadrature_registry if s.name == "ProductQuadrature"
    ).build({"n_mu": 3, "n_phi": 5})

    # Symmetry alone would let an S^2 rule serve a 1-D slab.
    assert slab.admits_symmetry(lebedev)
    assert not slab.admits_domain(lebedev)

    # Domain alone would let an asymmetric rule serve a cylinder.
    assert cylinder.admits_domain(odd_product)
    assert not cylinder.admits_symmetry(odd_product)


@pytest.mark.foundation
def test_owed_symmetry_selects_by_azimuthal_parity() -> None:
    """D_2h needs mirror planes at 0 AND pi/2; D_{n h} carries planes at
    k*pi/n, so pi/2 is present iff n is even.

    This is the whole content of the cylinder's admissibility rule, and
    it is ERR-042 stated structurally instead of as a hand-written
    guard.
    """
    cylinder = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
    product = next(
        s for s in quadrature_registry if s.name == "ProductQuadrature"
    )
    for n_phi in range(2, 10):
        measure = product.build({"n_mu": 3, "n_phi": n_phi})
        assert cylinder.admits_symmetry(measure) == (n_phi % 2 == 0), (
            f"n_phi={n_phi}: parity is the discriminator"
        )


@pytest.mark.foundation
def test_selector_asks_the_nodes_not_the_declared_tag() -> None:
    """Stage 1 must not route through ``measure.invariance_group``.

    A declaration may be true-but-not-maximal, and one shipped rule is:
    ``gauss_legendre_on_mu`` declares SO(2), the group its domain was
    quotiented BY. Z_2 is a reflection, hence not a subgroup of the
    rotation group SO(2) — so a lattice query against the declared tag
    would reject Gauss-Legendre for a slab, which is the one rule a
    slab must accept. Asking the nodes directly cannot fail that way.
    """
    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    gl = next(
        s for s in quadrature_registry if s.name == "GaussLegendre1D"
    ).build({"n": 8})

    # The lattice route is the trap...
    assert not slab.discrete_residual.is_subgroup_of(gl.invariance_group)
    # ...and the nodes give the right answer.
    assert slab.admits_symmetry(gl)

    # End-to-end: the slab does get Gauss-Legendre.
    _, log = select_quadrature("slab", target_degree=15)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"

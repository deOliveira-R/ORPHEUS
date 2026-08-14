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
    """Every spec must have a callable factory and inversion.

    A spec deliberately carries NO invariance group: a rule's symmetry
    is parameter-dependent, so a parameter-free field cannot state it
    truthfully. The selector computes it from the built measure.

    ⭐ It carries no ``expected_node_count`` either, for the same reason
    one level down. The cost is a property of the rule the selector has
    *already built* (stages 0 and 1 force the instantiation), so a
    formula-shaped second source could only ever drift from it. `[M]`
    2026-08-14 the two agreed on all 25 shipped configurations — which is
    what makes a twin dangerous rather than safe — and the first family
    that would have broken it already exists: ``folded_product``
    quotients by a mirror, so ``n_mu * n_phi`` over-counts it by 2x.

    Both absence assertions below are the load-bearing half of this
    gate: they are what stops either field being re-added by a
    contributor who reads the selector and reaches for a cheap tag.
    """
    for spec in quadrature_registry:
        assert isinstance(spec, QuadratureSpec)
        assert callable(spec.factory)
        assert callable(spec.degree_of_exactness_for)
        assert not hasattr(spec, "invariance_group")
        assert not hasattr(spec, "expected_node_count")
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
        (1, 2),     # min even N=2
        (3, 4),     # N=4
        (4, 6),     # need N>=5; round up to even -> 6
        (11, 12),   # N=12 — see the docstring: NOT the family's last order
    ],
)
def test_ls_sn_invert(target: int, expected_N: int) -> None:
    r"""The ``N - 1`` inversion, pinned at four targets.

    ⛔ The last row's comment read "N=12 — the LAST order the family can
    serve" until 2026-08-14. `[M]` the family serves through :math:`S_{18}`
    (:math:`S_{20}` is the smallest even order whose per-orbit weight solve
    goes negative), and the very next test in this file
    (``test_ls_sn_invert_serves_the_moment_matched_frontier``) asserts
    :math:`S_{14}` and :math:`S_{16}` — so the refutation was already sitting
    twenty lines below the claim. ``12`` was the pre-#337 convention seed's
    frontier; #337 moved it and this comment did not follow.

    ⚠ These rows pin the inversion, **not** the cheapest order that meets the
    target. `[M]` the two differ at 6 of the 9 buildable orders because the
    realized degree is build-measured and exceeds :math:`N - 1`: target 4 is
    met by :math:`S_4` (realized degree 5, 24 nodes) and this returns
    :math:`S_6` (48 nodes). See :func:`_ls_sn_invert`'s docstring — the
    over-shoot is a known cost defect, and these rows will need re-posing when
    it is fixed.
    """
    params = _ls_sn_invert(target)
    assert params == {"sn_order": expected_N}


@pytest.mark.foundation
@pytest.mark.parametrize(
    "target,expected", [(13, {"sn_order": 14}), (15, {"sn_order": 16})]
)
def test_ls_sn_invert_serves_the_moment_matched_frontier(
    target: int, expected: "dict[str, int]",
) -> None:
    r"""⭐ RE-POSED at #337 — exactly as the #327 row's own docstring
    predicted: "if a future node choice pushes the frontier up, the
    inverter follows it and these targets start being served — at which
    point this row reddens and says so, which is the correct outcome."

    The moment-matched seed (#337) pushed the positivity frontier from
    :math:`S_{12}` to :math:`S_{18}`, so targets 13 and 15 are served
    again (:math:`S_{14}` achieves degree 15 and :math:`S_{16}` degree
    15 — the ``n_min = target + 1`` inverse stays SAFE, if non-tight,
    because ``deg(N) >= N - 1`` at every buildable order). The inverter
    discovered this by ATTEMPTING the construction — no literal moved.
    """
    assert _ls_sn_invert(target) == expected


@pytest.mark.foundation
def test_ls_sn_invert_REFUSES_above_the_positivity_frontier() -> None:
    r"""Above :math:`S_{18}` the solve has no POSITIVE solution
    (`[M]` #337: −2.191e-4 at :math:`S_{20}`, 50-digit-confirmed), so
    the inverter returns ``None`` — the selector's documented "cannot
    reach target_degree with any supported parameters" channel. The
    inverter discovers this by attempting the construction, never by
    comparing against a literal frontier.
    """
    assert _ls_sn_invert(21) is None


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
    assert measure.invariance_group == SubgroupOfO3.Mirror("x")
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
@pytest.mark.parametrize(
    "geometry", ["slab", "sphere", "cylinder", "cartesian2d"],
)
def test_the_winner_really_is_the_argmin_over_the_survivors(
    geometry: str,
) -> None:
    r"""Stage 4 states :math:`\arg\min\{n(Q)\}` — pin the ``argmin``.

    The four ``test_select_*`` gates pin stage 4 by EXAMPLE (Lebedev's 14
    nodes beat 48 and 18 on ``cartesian2d``). This pins the *criterion*:
    across every geometry, no rule that survived stages 0-3 has fewer
    nodes than the one that was returned.

    The survivor set is read from the log's own rejection list rather than
    re-derived here, so this does not re-implement the filter — it checks
    only the minimisation, which is the one part of the criterion no
    per-example gate states in general.

    Independent of the cost SOURCE by construction: it would redden for a
    reversed sort, an off-by-one in the tie-break, or a candidate tuple
    whose cost slot drifted away from the measure it accompanies.

    ⚠ **2 of these 4 rows cannot catch an argmin bug, and say so here
    rather than being counted as coverage** (``vv-principles`` #20). `[M]`
    2026-08-14 at ``target_degree=5``: ``slab`` and ``sphere`` leave a
    **singleton** survivor set (only ``GaussLegendre1D``, 4 nodes — the
    three :math:`S^2` rules are rejected on domain), so no ordering over it
    can be wrong and the loop's assertion is trivially true. Only
    ``cylinder`` and ``cartesian2d`` have a real choice
    (``LebedevSphere``\ =14, ``ProductQuadrature``\ =18,
    ``LevelSymmetricSN``\ =48; winner 14).
    ⟹ mutation-verified at **2 rows**, not 4: `[M]` replacing the stage-4
    sort with ``reverse=True`` reddens 3 gates (both real rows here, plus
    ``test_select_cartesian2d_prefers_lebedev_over_ls_sn``) and leaves the
    two singleton rows green. The singleton rows are kept anyway — they
    still assert that every geometry admits *something* and that each
    survivor's inversion is non-``None`` — but they are not evidence about
    the minimisation.
    """
    measure, log = select_quadrature(geometry, target_degree=5)
    rejected_names = {name for name, _ in log.rejected}
    survivors = [s for s in quadrature_registry if s.name not in rejected_names]
    assert survivors, f"{geometry}: nothing survived, so there is no argmin"

    for spec in survivors:
        params = spec.degree_of_exactness_for(5)
        assert params is not None, (
            f"{spec.name} survived stage 2 but its inversion returns None"
        )
        assert measure.n_points <= spec.build(params).n_points, (
            f"{geometry}: {log.chosen_spec.name if log.chosen_spec else '?'} "
            f"was returned with {measure.n_points} nodes, but {spec.name} "
            f"survived every stage with fewer"
        )


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
    stages = set()
    for name, reason in log.rejected:
        assert isinstance(name, str)
        assert isinstance(reason, str)
        # The CLAIM: every rejection names the stage that rejected it. That is
        # what this row is for, and it holds for all of them.
        assert any(
            marker in reason
            for marker in ("domain mismatch", "V mismatch", "symmetry",
                           "structural")
        ), f"{name} was rejected without naming a stage: {reason!r}"
        stages.add("V" if "V mismatch" in reason else "domain")

    # ⭐ RE-POSED at #327. This used to assert ``"domain mismatch" in reason``
    # for EVERY rejection — true only incidentally, because every rule could
    # reach every target degree, so the V stage never rejected anything and
    # domain was always the first refusal. The level-symmetric family is now
    # BOUNDED (no positive solution above S12), so at ``target_degree=15`` it
    # is rejected at V — which runs FIRST by design, since V fixes the
    # parameters the domain check needs nodes for.
    #
    # The domain claim is kept where it still applies: any rule that CAN reach
    # the target still has to be rejected for living on S^2.
    for name, reason in log.rejected:
        if "V mismatch" not in reason:
            assert "domain mismatch" in reason and "[-1,1]" in reason, (
                f"{name} reached the domain stage but its reason does not name "
                f"the S^2 vs [-1,1] mismatch: {reason!r}"
            )
    assert "domain" in stages, (
        "every rule was rejected at V, so this row no longer exercises the "
        "domain stage at all — pick a target_degree the S^2 rules can reach"
    )


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
    r"""Stage 1 must not route through ``measure.invariance_group``.

    **Re-posed 2026-08-03, because its original witness was DISSOLVED by a
    fix — and that is the more interesting story.** It used to read: a
    declaration may be true-but-not-maximal, and ``gauss_legendre_on_mu``
    declares :math:`SO(2)` (the group its domain was quotiented BY), so a
    lattice query against the declared tag would reject Gauss-Legendre for
    a slab — the one rule a slab must accept.

    That trap existed *because the declaration was wrong*. The tag has
    since been corrected to :math:`\sigma_x` — the symmetry the measure
    actually has, rather than the one its reduction spent — so the lattice
    route and the node route now agree, and the original assertion
    (``not Mirror('x').is_subgroup_of(Mirror('x'))``) is simply false.

    A gate whose witness can be removed by fixing an unrelated bug is
    pinning a COINCIDENCE, not a mechanism. So this now injects a
    deliberately weakened-but-TRUE declaration and shows stage 1 is
    unmoved. That cannot be dissolved: it is the "compute the answer,
    never read a declaration — even a true one" rule, and the injected tag
    IS true, merely not maximal.
    """
    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    gl = next(
        s for s in quadrature_registry if s.name == "GaussLegendre1D"
    ).build({"n": 8})

    # The honest declaration and the nodes now agree — no trap left here.
    declared = gl.invariance_group
    assert declared is not None
    assert slab.discrete_residual.is_subgroup_of(declared)
    assert slab.admits_symmetry(gl)

    # Weaken the DECLARATION to something true-but-not-maximal. Every
    # measure is Trivial-invariant, so this is not a lie — it is exactly
    # the "true but under-claiming" case a declaration is allowed to be.
    understated = gl.with_metadata(invariance_group=SubgroupOfO3.Trivial)
    # A lattice query against the tag would now REJECT the one rule a slab
    # must accept...
    weakened = understated.invariance_group
    assert weakened is not None
    assert not slab.discrete_residual.is_subgroup_of(weakened)
    # ...while stage 1, which asks the NODES, is unmoved.
    assert slab.admits_symmetry(understated)

    # End-to-end: the slab does get Gauss-Legendre.
    _, log = select_quadrature("slab", target_degree=15)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"

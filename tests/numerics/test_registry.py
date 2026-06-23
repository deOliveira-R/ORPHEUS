"""Foundation tests for ``orpheus.numerics.quadrature.registry``.

Verifies the four-stage selection chain (G → V → structural → minimum
points) end-to-end on the four canonical geometry tags shipped today
(``slab``, ``sphere``, ``cylinder``, ``cartesian2d``), plus the
explainability log and the negative paths (no-rule-fits, unknown flags,
unknown geometry).

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
    GEOMETRY_GROUPS,
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
    """Every spec must have callable inversion / node-count and a
    valid :class:`SubgroupOfO3` invariance group."""
    for spec in quadrature_registry:
        assert isinstance(spec, QuadratureSpec)
        assert callable(spec.factory)
        assert callable(spec.degree_of_exactness_for)
        assert callable(spec.expected_node_count)
        assert isinstance(spec.invariance_group, SubgroupOfO3)
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
def test_geometry_groups_table() -> None:
    """The four geometries the selector handles must each carry a
    well-defined :class:`SubgroupOfO3`."""
    assert GEOMETRY_GROUPS == {
        "slab": SubgroupOfO3.SO2,
        "sphere": SubgroupOfO3.SO2,
        "cylinder": SubgroupOfO3.SO2,
        "cartesian2d": SubgroupOfO3.OctahedralOh,
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
    """slab geometry is :math:`SO(2)`-invariant; only :math:`SO(2)`-tagged
    rules pass the G stage. GL1D wins on minimum-points (8 << 16
    GL × n_phi for the product rule)."""
    measure, log = select_quadrature("slab", target_degree=15)

    assert isinstance(measure, DiscreteMeasure)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"
    assert log.chosen_parameters == {"n": 8}
    assert measure.n_points == 8
    assert measure.degree_of_exactness == 15
    assert measure.invariance_group == SubgroupOfO3.SO2
    # Lebedev and LS_N must have been rejected at the G stage.
    rejected_names = {name for name, _ in log.rejected}
    assert "LebedevSphere" in rejected_names
    assert "LevelSymmetricSN" in rejected_names


@pytest.mark.foundation
def test_select_sphere_returns_gauss_legendre() -> None:
    """1-D radial spherical SN reduces to GL on :math:`\\mu_r` — the
    geometry's group is :math:`SO(2)`, matching the slab path."""
    measure, log = select_quadrature("sphere", target_degree=15)

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"
    assert log.chosen_parameters == {"n": 8}
    assert measure.degree_of_exactness == 15


@pytest.mark.foundation
def test_select_cylinder_with_level_structured_returns_product() -> None:
    """Cylindrical SN sweep needs polar-level structure (Bailey 2009
    Eq. 50). Among :math:`SO(2)`-tagged rules, only ProductQuadrature
    has ``level_structured=True``."""
    measure, log = select_quadrature(
        "cylinder", target_degree=4, level_structured=True
    )

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "ProductQuadrature"
    assert log.chosen_parameters == {"n_mu": 3, "n_phi": 5}
    # GL1D should be rejected at the structural stage (no level
    # structure), and Oh-invariant rules at the G stage.
    rejected_names = {name for name, _ in log.rejected}
    assert "GaussLegendre1D" in rejected_names
    # Verify the rejection reason for GL1D mentions "structural".
    gl_reason = next(
        reason for name, reason in log.rejected
        if name == "GaussLegendre1D"
    )
    assert "structural" in gl_reason.lower()


@pytest.mark.foundation
def test_select_cartesian2d_prefers_lebedev_over_ls_sn() -> None:
    """2-D Cartesian carries :math:`O_h`. Both Lebedev and LS_N pass G
    and V; Lebedev wins on minimum-points (14 nodes for order 5 vs.
    48 nodes for LS_6).

    The :class:`SelectionLog` must record this preference: LS_N is NOT
    in the rejected list (it is a valid candidate, just more
    expensive).
    """
    measure, log = select_quadrature("cartesian2d", target_degree=5)

    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "LebedevSphere"
    assert log.chosen_parameters == {"order": 5}
    assert measure.n_points == 14

    # The :math:`SO(2)`-tagged rules (GL1D, Product) must be rejected
    # at G. LS_N must NOT be rejected — it was a valid candidate that
    # lost on cost.
    rejected_names = {name for name, _ in log.rejected}
    assert "GaussLegendre1D" in rejected_names
    assert "ProductQuadrature" in rejected_names
    assert "LevelSymmetricSN" not in rejected_names


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
    assert log.geometry_group == SubgroupOfO3.SO2
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
    failing stage (G / V / structural)."""
    _, log = select_quadrature("slab", target_degree=15)
    for name, reason in log.rejected:
        assert isinstance(name, str)
        assert isinstance(reason, str)
        # All rejections for slab are at G (Lebedev, LS_N).
        assert "G mismatch" in reason


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
def test_incompatible_structural_flags_explainable() -> None:
    """Asking for ``half_range_clean=True`` together with
    ``level_structured=True`` on a slab geometry: the slab's
    :math:`SO(2)` group rejects both Oh-tagged rules at G; the
    remaining :math:`SO(2)`-tagged rules are GL1D (no level structure)
    and ProductQuadrature (does have it). The selector picks Product;
    log shows GL1D rejected at structural, Lebedev/LS_N rejected at G.
    """
    measure, log = select_quadrature(
        "slab",
        target_degree=4,
        level_structured=True,
        half_range_clean=True,
    )
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "ProductQuadrature"

    # GL1D rejected at structural (no level structure).
    gl_reason = next(
        reason for name, reason in log.rejected
        if name == "GaussLegendre1D"
    )
    assert "structural" in gl_reason.lower()


@pytest.mark.foundation
def test_truly_incompatible_flags_raises() -> None:
    """Cylinder + ``axis_aligned=True``: among :math:`SO(2)`-tagged
    rules, GL1D has axis_aligned=True but lacks level structure (which
    the cylinder usually wants); Product has level_structured=True
    but axis_aligned=False. Asking for BOTH at once leaves no
    candidate and must raise."""
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

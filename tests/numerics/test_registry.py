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

from dataclasses import replace

import numpy as np
import pytest

from orpheus.numerics.exactness import UNIFORM_ON_SPHERE
from orpheus.numerics.generating_measure import CHEBYSHEV_T, LEGENDRE
from orpheus.numerics.manifold import COSINE_INTERVAL, SPHERE, Quotient
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
from orpheus.numerics.quadrature.rules_sphere import level_symmetric_sn
from orpheus.numerics.quadrature.registry import (
    _gl1d_invert,
    _lebedev_invert,
    _lebedev_orders,
    _ls_sn_invert,
    _product_invert,
    _registry,
)
from orpheus.numerics.symmetry import SubgroupOfO3

# 2.2b gates (test-architect, 2026-09-02) — the draft's imports, verbatim
import numpy as np
import pytest

from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.manifold import (
    COSINE_INTERVAL,
    SPHERE,
    Ball,
    Quotient,
    barycentre,
    quotient_onto,
)
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.registry import GEOMETRY_ANGULAR_SYMMETRY
from orpheus.numerics.quadrature.rules_1d import (
    gauss_legendre_on_mu,
    gauss_legendre_on_polar_orbit,
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

    The SPENT stabiliser fixes the angular domain; the UNSPENT finite
    symmetry — what the solution still has in the local frame — is the fold
    licence (R3 of #434, 2026-09-03); the OWED closure is what a quadrature
    must realize as an ordinate permutation for a reflecting face. Recording
    only the spent half (as this table did until 2026-08-02) makes the
    selection gate unsatisfiable — no finite point set on S^2 is
    SO(2)-closed; recording owed without unspent (until R3) let the face
    closure be read as a fold licence.
    """
    assert GEOMETRY_ANGULAR_SYMMETRY == {
        # The spent group names the SAME axis as the owed mirror's normal:
        # the polar marginal embeds along x, and a slab spends O(2)_x.
        "slab": AngularSymmetry(
            spent=SubgroupOfO3.O2("x"),
            unspent=SubgroupOfO3.Trivial,
            owed=SubgroupOfO3.Mirror("x"),
        ),
        "sphere": AngularSymmetry(
            spent=SubgroupOfO3.O2("x"),
            unspent=SubgroupOfO3.Trivial,
            owed=SubgroupOfO3.Mirror("x"),
        ),
        # z-uniform (even in the axial cosine mu: sigma_z) AND azimuthally
        # symmetric (even in the azimuthal cosine xi: sigma_y) — D_1h.
        "cylinder": AngularSymmetry(
            spent=SubgroupOfO3.Trivial,
            unspent=SubgroupOfO3.Dnh(1),
            owed=SubgroupOfO3.Dnh(2),
        ),
        # z-uniform only: even in mu_z, and in nothing else (D1 of #434).
        "cartesian2d": AngularSymmetry(
            spent=SubgroupOfO3.Trivial,
            unspent=SubgroupOfO3.Mirror("z"),
            owed=SubgroupOfO3.Dnh(2),
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
    r"""A target above every order SciPy serves has no solution — the
    inverter MUST return ``None`` (not raise).

    The threshold is **derived** from the discovered set, so this proves a
    refusal whatever SciPy's table is, instead of pinning a number that
    expires the next time SciPy grows one.

    ⛔ This asserted ``_lebedev_invert(48) is None`` and
    ``_lebedev_invert(100) is None`` until 2026-08-14, against a frozen
    18-tuple topping at order 47. `[M]` SciPy 1.17.1 serves 32 orders to
    **131**, so both were pinning the defect: at HEAD 48 resolves to order
    53 and 100 to order 101, both of which build. The gate was *correct
    about the code and wrong about the world* — exactly what a stale
    tabulation buys, and the reason the set is now asked rather than
    stored.
    """
    top = max(_lebedev_orders())
    assert _lebedev_invert(top) == {"order": top}, (
        "the largest discovered order must itself be reachable"
    )
    assert _lebedev_invert(top + 1) is None
    assert _lebedev_invert(top + 1000) is None


#: The slab's angular domain — what a rule must DECLARE to clear stage 0
#: since tracker 2.4 (2026-09-01). A chart-level rule on ``[-1,1]`` is
#: refused there, so every test-local 1-D spec below reads its rule on the
#: orbit space exactly as the shipped ``GaussLegendre1D`` entry does.
_SLAB_ORBIT_SPACE = SPHERE.quotient(SubgroupOfO3.O2("x"))


def _on_slab(measure: DiscreteMeasure) -> DiscreteMeasure:
    """A 1-D rule on the chart, read on the slab's orbit space (with the
    residual :math:`\\sigma_x` the shipped adopter also re-tags)."""
    return measure.on_orbit_space(_SLAB_ORBIT_SPACE).with_metadata(
        invariance_group=SubgroupOfO3.Mirror("x"),
    )


def _gauss_chebyshev_spec() -> QuadratureSpec:
    """A Gauss-Chebyshev entry, shaped exactly like the shipped GL one.

    Deliberately NOT registered in production: it could never be selected
    for any transport geometry (they all integrate unweighted), so a
    permanent entry would be a rule that builds on every call and can never
    win. It lives here, where the selector's own ``registry=`` override
    puts it in front of the stage that must refuse it.
    """
    return QuadratureSpec(
        name="GaussChebyshev1D",
        factory=lambda n: _on_slab(CHEBYSHEV_T.gauss(n)),
        parameters={"n": int},
        degree_of_exactness_for=_gl1d_invert,   # same deg = 2n - 1
        positive_weights=True,
        axis_aligned=True,
        level_structured=False,
        half_range_clean=True,
    )


@pytest.mark.foundation
def test_gauss_chebyshev_clears_every_stage_except_the_reference() -> None:
    r"""The reference conjunct is LOAD-BEARING: nothing else catches this.

    Following ``test_the_two_stages_are_independent_and_both_load_bearing``:
    a conjunct earns its place by exhibiting an input the *other* conjuncts
    admit. `[M]` 2026-08-14, ``CHEBYSHEV_T.gauss(4)`` against ``slab``:

    ==================  =======================================  ========
    stage               reading                                  verdict
    ==================  =======================================  ========
    0 domain            ``support == S^2/O2_x`` (declared;      **admits**
                        the chart-level rule is refused here
                        since tracker 2.4)
    1 symmetry          nodes invariant under the owed           **admits**
                        :math:`\sigma_x` (computed, not
                        declared)
    2 degree            ``2n - 1 = 7 >= 5``                      **admits**
    3 structural        positive weights, axis-aligned,          **admits**
                        half-range clean
    2' **reference**    exact against ``chebyshev_t``, query     **REFUSES**
                        wants ``legendre``
    ==================  =======================================  ========

    So without this conjunct the rule is admissible, and list position
    alone decides whether it wins.

    What it would cost: `[M]` at :math:`n = 4` both rules advertise
    algebraic degree 7, and on :math:`q(x) = x^6` — comfortably inside that
    degree — Legendre returns :math:`0.285714` (exact :math:`2/7`) while
    Chebyshev returns :math:`0.981748`. An error of **0.696**, roughly
    :math:`3.4\times` the true value, from a rule that is not broken and is
    delivering its full advertised accuracy against the measure it was
    built for.
    """
    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    gc = _gauss_chebyshev_spec().build({"n": 4})

    assert slab.admits_domain(gc), "premise: stage 0 must admit it"
    assert slab.admits_symmetry(gc), "premise: stage 1 must admit it"
    assert gc.exactness is not None
    assert gc.exactness.degree >= 5, "premise: stage 2 must admit it"
    assert bool((gc.weights > 0).all()), "premise: stage 3 must admit it"

    # ...and the reference is the one thing that differs.
    assert gc.exactness.reference != slab.reference
    assert gc.exactness.reference.name == "chebyshev_t"
    assert slab.reference.name == "legendre"

    # The quantified consequence, recomputed here rather than quoted.
    q = LEGENDRE.gauss(4)
    exact_unweighted = 2.0 / 7.0
    legendre_value = float(np.sum(q.weights * q.nodes**6))
    chebyshev_value = float(np.sum(gc.weights * gc.nodes**6))
    assert legendre_value == pytest.approx(exact_unweighted, abs=1e-14)
    assert abs(chebyshev_value - exact_unweighted) == pytest.approx(
        0.696, abs=1e-3
    )


@pytest.mark.foundation
def test_a_rule_exact_against_the_wrong_measure_is_refused() -> None:
    """The selector must reject Gauss-Chebyshev for an unweighted query.

    Two registries, because the two failure modes are different: alone it
    must be REFUSED (not merely out-competed), and beside Gauss-Legendre it
    must still be refused rather than winning on a tie-break.
    """
    gc_only = [_gauss_chebyshev_spec()]
    with pytest.raises(QuadratureSelectionError) as excinfo:
        select_quadrature("slab", target_degree=5, registry=gc_only)

    reason = next(
        r for name, r in excinfo.value.log.rejected
        if name == "GaussChebyshev1D"
    )
    assert "V mismatch" in reason
    assert "chebyshev_t" in reason and "legendre" in reason

    # Chebyshev FIRST, so registry order would hand it the win if the
    # conjunct were absent (ties break on order, and both cost 4 nodes).
    both = [_gauss_chebyshev_spec()] + [
        s for s in quadrature_registry if s.name == "GaussLegendre1D"
    ]
    measure, log = select_quadrature("slab", target_degree=5, registry=both)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"
    assert "GaussChebyshev1D" in {name for name, _ in log.rejected}
    assert measure.exactness is not None
    assert measure.exactness.reference == GEOMETRY_ANGULAR_SYMMETRY[
        "slab"
    ].reference


@pytest.mark.foundation
def test_a_rule_with_no_exactness_claim_at_all_is_refused() -> None:
    r"""No claim is not a weak claim — it is no certification whatsoever.

    A measure may carry ``exactness=None``. Nothing then establishes that
    it integrates *anything* exactly, so it cannot satisfy a V constraint
    at any target degree, including 0.

    ⚠ **This arm is reachable, but not by any rule the registry ships**, so
    it needs a constructed witness or it would be an unfalsifiable guard —
    `[M]` 2026-08-14, disabling it in-process left all 54 gates green
    before this test existed, which is precisely the "cannot fail, wearing
    an authoritative name" shape.

    The real future occupant is ``Quadrature.folded_product``: `[M]` it is
    the one shipped angular family whose measure carries
    ``exactness=None`` (and ``invariance_group=None``, on support
    ``'S^2/sigma_y'``). It cannot be registered today — its quotient
    support matches no geometry, so stage 0 would refuse it first — but
    when that is resolved this arm is what stops it being selected on a
    degree it never claimed.
    """
    claimless = QuadratureSpec(
        name="ClaimlessGL",
        factory=lambda n: replace(_on_slab(LEGENDRE.gauss(n)), exactness=None),
        parameters={"n": int},
        degree_of_exactness_for=_gl1d_invert,
        positive_weights=True,
        axis_aligned=True,
        level_structured=False,
        half_range_clean=True,
    )
    # It clears every other stage: same nodes, same support, same symmetry.
    probe = claimless.build({"n": 4})
    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    assert slab.admits_domain(probe) and slab.admits_symmetry(probe)
    assert probe.exactness is None

    for target in (0, 5):
        with pytest.raises(QuadratureSelectionError) as excinfo:
            select_quadrature(
                "slab", target_degree=target, registry=[claimless]
            )
        reason = next(
            r for name, r in excinfo.value.log.rejected
            if name == "ClaimlessGL"
        )
        assert "V mismatch" in reason
        assert "no exactness claim" in reason


@pytest.mark.foundation
def test_an_inversion_that_over_promises_is_caught_not_trusted() -> None:
    """A spec cannot talk its way past the V stage.

    ``degree_of_exactness_for`` is a per-spec callable the selector used to
    take on faith: it returned parameters, and the rule built from them was
    assumed to meet the target. It is a *search hint*, and the authority is
    the claim the built rule carries.

    Here a spec's inversion always answers ``n = 1`` — Gauss-Legendre at one
    node is exact to degree 1 — while claiming to serve any target. The
    selector must refuse it rather than return a rule three degrees short.

    ⭐ This is the arm that would have caught ``_ls_sn_invert`` had its
    ``N - 1`` inversion erred the other way. `[M]` it over-shoots (returns
    an order 2 higher than needed at 6 of 9 buildable orders), which is safe
    and merely expensive; the same class of staleness pointing the other
    direction returns a rule that silently misses the requested accuracy,
    and nothing downstream re-checks.
    """
    lying = QuadratureSpec(
        name="OverPromisingGL",
        factory=lambda n: _on_slab(LEGENDRE.gauss(n)),
        parameters={"n": int},
        degree_of_exactness_for=lambda _target: {"n": 1},
        positive_weights=True,
        axis_aligned=True,
        level_structured=False,
        half_range_clean=True,
    )
    with pytest.raises(QuadratureSelectionError) as excinfo:
        select_quadrature("slab", target_degree=5, registry=[lying])

    reason = next(
        r for name, r in excinfo.value.log.rejected
        if name == "OverPromisingGL"
    )
    assert "V mismatch" in reason
    assert "degree 1" in reason and "target_degree=5" in reason

    # And the honest twin passes, so the gate is not just refusing everything.
    honest = replace(lying, name="HonestGL", degree_of_exactness_for=_gl1d_invert)
    measure, log = select_quadrature("slab", target_degree=5, registry=[honest])
    assert log.chosen_spec is not None and log.chosen_spec.name == "HonestGL"
    assert measure.exactness is not None
    assert measure.exactness.degree >= 5


@pytest.mark.foundation
def test_every_registered_rule_speaks_one_of_the_two_reference_measures(
) -> None:
    r"""Close the reference vocabulary, and gate the closure.

    The conjunct compares references with ``==``, and that is sound here
    only because the registry's references are drawn from the two
    canonical constants a geometry can ask for. It is **not** sound in
    general: `[M]` 2026-08-14 the tree spells "Lebesgue on :math:`[-1,1]`"
    three mutually-unequal ways — ``LEGENDRE`` (``'legendre'``),
    ``LEGENDRE.on(-1, 1)`` (``'legendre_on[-1.0,1.0]'``, because the
    identity case does not canonicalise) and the anonymous
    ``UniformMeasure`` that ``equispaced`` builds
    (``'uniform([-1.0,1.0])'``). All three denote the same measure; none
    compares equal to another, and cross-class ``__eq__`` returns ``False``
    before it looks at a single field.

    ⟹ this gate makes the assumption explicit and enforced instead of
    accidental. The day a rule arrives spelling its reference the third
    way, it reddens **here, at registration**, rather than being silently
    mis-selected somewhere downstream. Repairing ``ReferenceMeasure``
    equality across realizations is the real fix and is filed separately;
    until then this is the fence around the shortcut.
    """
    admissible = {LEGENDRE, UNIFORM_ON_SPHERE}
    for spec in quadrature_registry:
        params = spec.degree_of_exactness_for(5)
        assert params is not None
        claim = spec.build(params).exactness
        assert claim is not None, (
            f"{spec.name} carries no exactness claim, so no query can ever "
            f"establish that its degree is against the right measure"
        )
        assert claim.reference in admissible, (
            f"{spec.name} is exact against {claim.reference.name!r}, which is "
            f"outside the closed vocabulary {{legendre, uniform(S^2)}} the "
            f"reference conjunct compares by equality. Either it is a rule "
            f"no geometry can ask for, or it spells a canonical measure a "
            f"second way -- and the second case is a silent mis-selection"
        )


@pytest.mark.foundation
def test_the_registry_is_frozen_and_cannot_carry_a_duplicate_name() -> None:
    r"""Two rules cannot share a name, and the tuple cannot be appended to.

    ``name`` is not merely a label: :class:`SelectionLog` reports every
    rejection under it, and the theory page teaches
    ``dict(log.rejected)["ProductQuadrature"]`` as the way to read one.
    ``dict`` keeps the LAST value for a repeated key, so two rules sharing
    a name make one rejection **silently vanish** from that view while the
    log still lists both — a disappearing diagnostic, with nothing warning.

    `[M]` 2026-08-14 ``quadrature_registry`` was also the only module-level
    *list*-shaped registry in ``orpheus/``; every sibling
    (``LOSS_REPRESENTATIONS`` and friends) is a ``tuple``.

    The refusal is at **import**, not here — this test only proves it
    fires. A construction-time invariant that lives in a test can be
    skipped; one that raises on import cannot.
    """
    assert isinstance(quadrature_registry, tuple)
    with pytest.raises(AttributeError):
        quadrature_registry.append(_gauss_chebyshev_spec())  # type: ignore[attr-defined]

    names = [s.name for s in quadrature_registry]
    assert len(names) == len(set(names))

    twin = replace(quadrature_registry[0])
    with pytest.raises(ValueError, match="duplicate spec name"):
        _registry(*quadrature_registry, twin)

    # ...and the honest case is accepted, so the guard is not refusing all.
    assert _registry(*quadrature_registry) == quadrature_registry


@pytest.mark.foundation
def test_the_discovered_lebedev_set_matches_what_scipy_advertises() -> None:
    r"""Two independent channels of SciPy's own truth must agree.

    :func:`_lebedev_orders` discovers the admissible set by **constructing**
    — it asks for each odd order and keeps the ones that succeed. SciPy
    also **advertises** the set, in the ``NotImplementedError`` it raises
    for an unavailable order. Those are different code paths reading
    different data, so agreement is evidence and disagreement localises the
    fault:

    * discovery is a strict subset ⟹ the search window is truncating;
    * discovery is a strict superset ⟹ an order builds that SciPy does not
      list, i.e. our wrapper is reaching something undocumented.

    ⭐ This is the gate that makes the discovery trustworthy, and it is
    why the refusal gates above are allowed to derive their thresholds from
    :func:`_lebedev_orders` without becoming self-referential: those check
    the *inverter* against the set, this checks the *set* against SciPy.
    Without it, a discovery that silently returned a truncated tuple would
    keep every other Lebedev gate green — which is precisely the failure
    the frozen 18-tuple shipped from 2026-05-06 to 2026-08-14.

    ⭐ **This gate and the ceiling guard inside** :func:`_lebedev_orders`
    **catch disjoint halves of the same defect, and neither covers both.**
    `[M]` 2026-08-14, shrinking the search window in-process:

    ==============  ==================  ==================================
    ceiling         guard               this gate
    ==============  ==================  ==================================
    47 (a served    **FIRES** — the     never reached
    order)          window is provably
                    binding
    49, 51 (in the  **silent** — the    **BITES** — discovers max 47
    47→53 gap)      last order found    against an advertised 131
                    is 47, not the
                    ceiling
    ==============  ==================  ==================================

    So the guard is sound but incomplete: it can only notice a window that
    ends exactly on a served order. A ceiling that happens to land in one
    of SciPy's gaps truncates silently, and only comparing against the
    advertised list catches it. Do not retire either as redundant.
    """
    import re

    from scipy.integrate import lebedev_rule

    with pytest.raises(NotImplementedError) as excinfo:
        lebedev_rule(4)  # 4 is even, so never available in any SciPy
    advertised = {int(m) for m in re.findall(r"\d+", str(excinfo.value))}
    advertised.discard(4)  # the order we asked for, echoed in the message

    discovered = set(_lebedev_orders())
    assert discovered == advertised, (
        f"discovery and SciPy's advertised list disagree.\n"
        f"  only discovered: {sorted(discovered - advertised)}\n"
        f"  only advertised: {sorted(advertised - discovered)}"
    )
    assert discovered, "SciPy advertised no Lebedev orders at all"


@pytest.mark.foundation
@pytest.mark.parametrize(
    "target,expected_N",
    [
        (1, 2),     # deg(S_2) = 3
        (3, 2),     # deg(S_2) = 3 — exactly meets it
        (4, 4),     # deg(S_4) = 5
        (11, 10),   # deg(S_10) = 11 — S_12 also gives 11, and costs more
    ],
)
def test_ls_sn_invert(target: int, expected_N: int) -> None:
    r"""The inversion returns the CHEAPEST order that meets the target.

    ⛔ These rows read ``(3, 4)``, ``(4, 6)``, ``(11, 12)`` until
    2026-08-14, pinning an inversion that literally inverted
    :math:`\deg = N - 1`. That formula is a true **lower** bound, so the
    answers were safe — and `[M]` too expensive at 6 of the 9 buildable
    orders, because the realized degree is build-measured and exceeds
    :math:`N - 1`. Target 3 is met by :math:`S_2` (8 nodes) and the old
    inversion returned :math:`S_4` (24) — **3×**; target 11 by
    :math:`S_{10}` (120) against :math:`S_{12}` (168).

    Stage 4 ranks by node count, so the over-shoot was not cosmetic: it
    priced the family above its true cost and could lose it a
    minimum-points tie-break it should win.

    ⛔ The last row's comment also read "N=12 — the LAST order the family
    can serve". `[M]` the family serves through :math:`S_{18}`
    (:math:`S_{20}` is the smallest even order whose per-orbit weight solve
    goes negative), and the very next test in this file asserted
    :math:`S_{14}` — the refutation was sitting twenty lines below the
    claim. ``12`` was the pre-#337 convention seed's frontier.
    """
    params = _ls_sn_invert(target)
    assert params == {"sn_order": expected_N}


@pytest.mark.foundation
@pytest.mark.parametrize("target", list(range(0, 18)))
def test_ls_sn_invert_is_minimal_and_never_short(target: int) -> None:
    r"""The two halves of "cheapest that meets it", checked directly.

    Rather than pinning a table (which is what went stale), assert the
    *property*: whatever order comes back must (a) actually deliver the
    target degree, and (b) have no cheaper buildable order that also does.
    Node count :math:`N(N+2)` is increasing in :math:`N`, so "cheaper"
    is "smaller :math:`N`".

    ⭐ This is the gate that would have caught the over-shoot, and it is
    written against the *contract* rather than against measured values so
    that a future seed change moves it automatically — the same reason
    :func:`_ls_sn_invert` asks the family for the frontier instead of
    tabulating it.
    """
    params = _ls_sn_invert(target)
    if params is None:
        pytest.skip(f"family cannot serve degree {target}")
    chosen = params["sn_order"]

    measure, _ = level_symmetric_sn(chosen)
    assert measure.degree_of_exactness is not None
    assert measure.degree_of_exactness >= target, "returned order is SHORT"

    for smaller in range(2, chosen, 2):
        cheaper, _ = level_symmetric_sn(smaller)
        assert (
            cheaper.degree_of_exactness is None
            or cheaper.degree_of_exactness < target
        ), (
            f"S_{smaller} ({smaller * (smaller + 2)} nodes) also reaches "
            f"degree {target}, so S_{chosen} "
            f"({chosen * (chosen + 2)} nodes) is not minimal"
        )


@pytest.mark.foundation
@pytest.mark.parametrize(
    "target,expected", [(13, {"sn_order": 14}), (15, {"sn_order": 14})]
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
    again. The inverter discovered this by ATTEMPTING the construction —
    no literal moved.

    ⭐ **Both targets are now served by the same order**, and that is the
    point of the second re-pose. :math:`S_{14}` achieves degree **15**, so
    it meets 13 *and* 15; the old ``n_min = target + 1`` inversion sent
    target 15 to :math:`S_{16}` (288 nodes) purely because it inverted a
    lower bound instead of asking what the rule delivers. `[M]`
    :math:`S_{16}` also achieves only degree 15 — so the old answer paid
    **64 extra nodes for nothing at all**, not even a degree it could use.

    ⛔ The expected value for target 15 read ``{"sn_order": 16}`` until
    2026-08-14, and the docstring justified it as "SAFE, if non-tight".
    Safe it was; the non-tightness was a live cost defect, because stage 4
    ranks by node count.
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
        spent=SubgroupOfO3.Trivial,
        unspent=SubgroupOfO3.Dnh(1),
        owed=SubgroupOfO3.Dnh(2),
    )
    # All three entries are reconstructible from the log, so a reader can
    # replay stages 0 and 1 without re-deriving the ledger.
    assert log.angular_symmetry.support == SPHERE
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
            assert "domain mismatch" in reason and "S^2/O2_x" in reason, (
                f"{name} reached the domain stage but its reason does not name "
                f"the S^2 vs S^2/O2_x mismatch: {reason!r}"
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
    the registry to contain only Lebedev and ask for a degree beyond
    everything the installed SciPy serves, which forces a V-stage
    rejection.

    ⛔ This asked for ``target_degree=100`` against a Lebedev "which tops
    out at order 47" until 2026-08-14. `[M]` SciPy 1.17.1 serves order
    **101**, so the selection succeeded and this gate stopped proving a
    refusal — it would have gone red on the fix, for the right reason.
    The target is now **derived** from the discovered set, so the gate
    keeps proving a refusal when SciPy's table grows instead of pinning
    a number that expires.
    """
    beyond_scipy = max(_lebedev_orders()) + 2
    only_lebedev = [
        s for s in quadrature_registry if s.name == "LebedevSphere"
    ]
    with pytest.raises(QuadratureSelectionError) as excinfo:
        select_quadrature(
            "cartesian2d", target_degree=beyond_scipy, registry=only_lebedev
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
    """The angular domain IS the orbit space S^2/K of the SPENT stabiliser.

    Deriving it rather than storing a second column is what keeps the
    domain and the spent group from drifting apart — they are one fact,
    so there is no state in which they can disagree.
    """
    # The orbit space itself, not its chart: until tracker 2.4 this row read
    # ``== COSINE_INTERVAL``, which is the chart of S^2/SO(2) and was also
    # exactly what a SPATIAL rule on [-1,1] declared.
    for axis in ("x", "y", "z"):
        assert AngularSymmetry(
            spent=SubgroupOfO3.O2(axis),
            unspent=SubgroupOfO3.Trivial,
            owed=SubgroupOfO3.Mirror(axis),
        ).support == SPHERE.quotient(SubgroupOfO3.O2(axis))

    assert AngularSymmetry(
        spent=SubgroupOfO3.Trivial,
        unspent=SubgroupOfO3.Mirror("z"),
        owed=SubgroupOfO3.Dnh(2),
    ).support == SPHERE

    # An unmapped quotient must refuse loudly rather than guess: a wrong
    # domain silently admits a rule of the wrong dimensionality.
    with pytest.raises(NotImplementedError, match="S\\^2/"):
        _ = AngularSymmetry(
            spent=SubgroupOfO3.OctahedralOh,
            unspent=SubgroupOfO3.Trivial,
            owed=SubgroupOfO3.Mirror("z"),
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
    (``not Mirror('x').contains(Mirror('x'))``) is simply false.

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
    assert declared.contains(slab.owed)
    assert slab.admits_symmetry(gl)

    # Weaken the DECLARATION to something true-but-not-maximal. Every
    # measure is Trivial-invariant, so this is not a lie — it is exactly
    # the "true but under-claiming" case a declaration is allowed to be.
    understated = gl.with_metadata(invariance_group=SubgroupOfO3.Trivial)
    # A lattice query against the tag would now REJECT the one rule a slab
    # must accept...
    weakened = understated.invariance_group
    assert weakened is not None
    assert not weakened.contains(slab.owed)
    # ...while stage 1, which asks the NODES, is unmoved.
    assert slab.admits_symmetry(understated)

    # End-to-end: the slab does get Gauss-Legendre.
    _, log = select_quadrature("slab", target_degree=15)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"


@pytest.mark.foundation
def test_the_reference_is_READ_off_the_catalogue_entry() -> None:
    r"""Pattern 2 (#429 tracker 3.1, 2026-09-02): the hat-box answer has ONE
    home, and it is the orbit space's own derivation.

    ``AngularSymmetry.reference`` used to TABULATE ``LEGENDRE`` for any
    rotation axis; it now reads :attr:`Quotient.reference` off the entry
    :attr:`support` derives — the same collapse ``support`` underwent at
    tracker 2.4. Three legs: the axial geometry reads the entry's field by
    IDENTITY; the geometry that spends nothing has the bare sphere as its
    domain, so its reference is the sphere's own Lebesgue measure (a
    property of the base, not a catalogue read — user-ruled); and the
    arm for an entry whose pushforward no shipped realization spells
    RAISES naming the missing work (``plan-authoring`` §6c: the gate lands
    with its witness — a constructed geometry spending a mirror).
    """
    from orpheus.numerics.manifold import Quotient

    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    assert isinstance(slab.support, Quotient)
    assert slab.reference is slab.support.reference
    assert slab.reference is LEGENDRE

    cylinder = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
    assert cylinder.support is SPHERE
    assert cylinder.reference is UNIFORM_ON_SPHERE

    spends_a_mirror = AngularSymmetry(
        spent=SubgroupOfO3.Mirror("y"),
        unspent=SubgroupOfO3.Trivial,
        owed=SubgroupOfO3.Trivial,
    )
    mirror_entry = SPHERE.quotient(SubgroupOfO3.Mirror("y"))
    assert spends_a_mirror.support == mirror_entry
    assert mirror_entry.reference is None
    with pytest.raises(NotImplementedError, match="no shipped ReferenceMeasure"):
        _ = spends_a_mirror.reference


# ---------------------------------------------------------------------------
# The 1-D geometries spend the full STABILISER — #432, 2026-09-02
# ---------------------------------------------------------------------------


@pytest.mark.foundation
def test_the_one_dimensional_geometries_spend_the_STABILISER_not_its_rotation_half() -> None:
    r"""``slab`` and ``sphere`` spend :math:`O(2)_x`, not :math:`SO(2)_x`.

    A slab is invariant under every rotation about its normal **and** under
    :math:`y \to -y` — the whole stabiliser of :math:`\hat e_x` — which is
    why its flux depends on :math:`\Omega` only through :math:`\mu`. Spending
    the rotation half alone would understate :math:`G^0`, and the recorded
    residual :math:`\sigma_x` would then not be :math:`G/G^0`.

    ⭐ The row that gives this content rather than restating the table: the
    residual :math:`\sigma_x` must lie **outside** the spent group, or
    "residual" is the wrong word. `[M]` 2026-09-02 :math:`\sigma_x \notin
    O(2)_x` (it flips the axis) while :math:`\sigma_y, \sigma_z \in O(2)_x`
    (they contain it) — so the recorded pair really is a complement, and
    exactly one of the three coordinate mirrors could have played that role.

    The 2-D/3-D geometries are the CONTROL: they spend ``Trivial`` and
    discretise the bare sphere, so a change that gave every geometry an axial
    group would redden here.
    """
    for geometry in ("slab", "sphere"):
        symmetry = GEOMETRY_ANGULAR_SYMMETRY[geometry]
        assert symmetry.spent == SubgroupOfO3.O2("x"), geometry
        assert symmetry.spent != SubgroupOfO3.SO2("x")
        assert symmetry.owed == SubgroupOfO3.Mirror("x")
        # nothing finite is left unspent by a 1-D reduction
        assert symmetry.unspent == SubgroupOfO3.Trivial
        # the owed closure is genuinely OUTSIDE the spent group
        assert not symmetry.spent.contains(symmetry.owed)
        # ... and the other two mirrors are INSIDE it, so "outside" is a
        # discriminating statement about sigma_x and not about mirrors.
        for inside in ("y", "z"):
            assert symmetry.spent.contains(SubgroupOfO3.Mirror(inside))

    for geometry in ("cylinder", "cartesian2d"):
        symmetry = GEOMETRY_ANGULAR_SYMMETRY[geometry]
        assert symmetry.spent == SubgroupOfO3.Trivial, geometry
        assert symmetry.support == SPHERE


@pytest.mark.foundation
@pytest.mark.parametrize("geometry", ["slab", "sphere"])
def test_the_support_is_DERIVED_as_the_stabiliser_orbit_space(geometry: str) -> None:
    r"""``AngularSymmetry.support`` is ``SPHERE.quotient(G^0)`` — the
    catalogue entry, derived, never a hand-written chart.

    Since #432 that entry is ``S^2/O2_x``, and the row asserts the whole
    tuple the derivation produces: the entry itself, its ``by``, its name,
    its realization and its pushforward reference. The reference matters
    independently — a degree of exactness on this orbit space is against
    Lebesgue on :math:`[-1,1]` (Archimedes' hat box), which is what makes
    ``LEGENDRE`` the right generating measure for the selector's stage-2
    comparison.
    """
    symmetry = GEOMETRY_ANGULAR_SYMMETRY[geometry]
    entry = SPHERE.quotient(SubgroupOfO3.O2("x"))
    support = symmetry.support
    # the narrowing IS part of the claim: ``support`` is typed ``Manifold``,
    # and the row asserts it is specifically a catalogued QUOTIENT.
    assert isinstance(support, Quotient)
    assert support == entry
    assert support.by == SubgroupOfO3.O2("x")
    assert support.name == "S^2/O2_x"
    assert support.realization == COSINE_INTERVAL
    assert symmetry.reference is LEGENDRE
    # the chart is NOT the orbit space (tracker 2.4's ruling, still binding)
    assert symmetry.support != COSINE_INTERVAL


@pytest.mark.foundation
@pytest.mark.parametrize("geometry", ["slab", "sphere"])
def test_selection_still_lands_on_GaussLegendre1D_over_the_stabiliser_orbit_space(
    geometry: str,
) -> None:
    r"""End-to-end: the renamed orbit space must not move the selector's
    answer, and the three :math:`S^2` rules must still be rejected on DOMAIN
    with a message naming the new entry.

    `[M]` 2026-09-02 at ``target_degree=5``: ``GaussLegendre1D`` with
    ``n = 4``, measure on ``S^2/O2_x`` carrying ``invariance_group ==
    sigma_x``, and all three sphere rules rejected as *"geometry '…'
    discretises S^2/O2_x, but the rule's nodes live on S^2"*.

    ⚠ The rejection REASON is asserted because it is the one place the
    entry's NAME reaches a user-facing string; a stage-0 comparison that had
    kept comparing charts would still choose GL1D and would print the old
    name, so the choice alone is not a discriminating observation.
    """
    measure, log = select_quadrature(geometry, target_degree=5)
    assert log.chosen_spec is not None
    assert log.chosen_spec.name == "GaussLegendre1D"
    assert log.chosen_parameters == {"n": 4}
    assert measure.support == SPHERE.quotient(SubgroupOfO3.O2("x"))
    assert measure.invariance_group == SubgroupOfO3.Mirror("x")

    rejected = dict(log.rejected)
    assert set(rejected) == {
        "LebedevSphere", "LevelSymmetricSN", "ProductQuadrature",
    }
    for name, reason in rejected.items():
        assert "S^2/O2_x" in reason, f"{name}: {reason}"
        assert "domain mismatch" in reason, f"{name}: {reason}"


# =============================================================================
# 2.2b — the Γ-slot: gates drafted by the test-architect (2026-09-02)
# =============================================================================

_g22b_AXES = ("x", "y", "z")


def _g22b_trivially_quotiented_product() -> DiscreteMeasure:
    from dataclasses import replace

    return replace(
        Quadrature.product(4, 8).measure,
        support=SPHERE.quotient(SubgroupOfO3.Trivial),
    )




def _g22b_rot(axis: str, theta: float) -> RigidMotion:
    return RigidMotion.rotation_about_axis(axis=np.eye(3)["xyz".index(axis)], angle=theta)


def _g22b_mirror(axis: str) -> RigidMotion:
    return RigidMotion.reflection(normal=np.eye(3)["xyz".index(axis)])


#: Pairwise-incommensurate angles.  A right-angle sample of a continuous
#: family generates C_4 and certifies what it never tested (ERR-072,
#: ``vv-principles`` #13); the CONTROL below is allowed to sample because it
#: can only REFUTE, and it is compared against a criterion that never does.
_g22b_INCOMMENSURATE = (1.0, float(np.sqrt(2.0)), float(np.e), 2.5, float(np.sqrt(7.0)))


# =============================================================================


# G12 — tests/numerics/test_registry.py
# =============================================================================
@pytest.mark.catches("ERR-081")
class TestStageZeroIsTheDescentArrowPlusTheUnspentSymmetry:
    r"""(vii) Stage 0 through the arrow, and the coverage test.

    `[M]` 2026-09-02, over 4 geometries × 7 rules = 28 rows, the arrow carve
    moved **4**: cylinder and cartesian2d each admitted the σ_y fold AND a σ_z
    fold, where equality refused both.  R3 of #434 (2026-09-03) then moved
    exactly **one** of those back: the σ_y fold on ``cartesian2d`` is REFUSED,
    because the fold spent a symmetry a z-uniform Cartesian solution does
    not have (its only unspent symmetry is σ_z) — the owed closure had been
    read as a fold licence (D1).

    Two legs, and BOTH have a live witness at landing (``plan-authoring``
    §6c):

    * the **descent arrow** ``S²/K → X`` must exist — remove it and the slab
      admits a σ_y fold it has no arrow onto;
    * the **coverage leg** ``H ⊆ Γ·K`` (:meth:`SubgroupOfO3.is_subset_of_product`)
      — remove it and the arrow ``S² → S²/O2_x`` admits a 1-D POLAR rule for
      a 2-D geometry, and cartesian2d admits the σ_y fold again.  It is
      TOTAL: the slab's own rule passes it as ``O2_x ⊆ {e}·O2_x``, so there
      is no equality short circuit to keep in front of it (until R3 there
      was one, because the predecessor ``owed.contains(H)`` refused the
      slab's own rule — ``σ_x ⊉ O(2)_x`` — and needed the identity case
      short-circuited).

    ⚠ **And the honest scope of the coverage leg at the selector tier.**
    `[M]` at ``select_quadrature`` it changes the REASON and never the
    selection: ``GaussLegendre1D`` IS registered and IS a fold (``S^2/O2_x``,
    ``H = O(2)_x``), and the shipped log for ``select_quadrature("cylinder",
    5)`` refuses it at stage 0 by the coverage clause (*"a fold by O2_x …
    unspent D_1h, spent Trivial"*); neuter the leg in-process and the refusal
    moves to stage 2's V conjunct while the choice, ``LebedevSphere(5)``,
    stays — 0 of 48 (geometry × degree) choices move.  So the leg's gate
    belongs HERE, at the ``admits_domain`` tier, and no end-to-end selector
    row may be credited to it.  (Until 2026-09-03 this paragraph said
    "nothing registered is a fold" and named stage 2 as the refuser — both
    false, the archivist measured.)
    """

    @pytest.mark.foundation
    def test_the_cylinder_admits_the_shipped_fold_and_the_slab_does_not(self):
        cyl = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
        slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        fold = Quadrature.folded_product(4, 8).measure
        assert cyl.admits_domain(fold)
        assert quotient_onto(cyl.support, fold.support) is not None
        assert fold.support.by.is_subset_of_product(gamma=cyl.unspent, kappa=cyl.spent)
        # the slab has no arrow at all: S^2/O2_x -> S^2/sigma_y does not exist
        assert quotient_onto(slab.support, fold.support) is None
        assert not slab.admits_domain(fold)

    @pytest.mark.foundation
    def test_a_z_fold_is_admitted_iff_its_mirror_lies_in_the_unspent_symmetry(self):
        """`[M]` ``σ_z ⊆ D_1h·{e}`` on the cylinder and ``σ_z ⊆ σ_z·{e}`` on
        cartesian2d, so a z-fold is admitted for both; the leg is what makes
        that a LATTICE verdict rather than a spelling coincidence — and the
        σ_y fold, ``σ_y ⊄ σ_z·{e}``, is refused for cartesian2d by the same
        leg (the D1 witness)."""
        zfold = Quadrature.product(4, 8).measure.quotient(SubgroupOfO3.Mirror("z"))
        assert zfold.support == SPHERE.quotient(SubgroupOfO3.Mirror("z"))
        yfold = Quadrature.folded_product(4, 8).measure
        for geometry in ("cylinder", "cartesian2d"):
            sym = GEOMETRY_ANGULAR_SYMMETRY[geometry]
            assert SubgroupOfO3.Mirror("z").is_subset_of_product(gamma=sym.unspent, kappa=sym.spent)
            assert sym.admits_domain(zfold), geometry
        plane = GEOMETRY_ANGULAR_SYMMETRY["cartesian2d"]
        assert not SubgroupOfO3.Mirror("y").is_subset_of_product(gamma=plane.unspent, kappa=plane.spent)
        assert not plane.admits_domain(yfold)                      # D1 closed
        assert GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain(yfold)
        # and the reason names the fold licence, not the arrow (the arrow
        # exists: S^2 -> S^2/sigma_y is the entry's own quotient map)
        reason = plane.domain_refusal(yfold)
        assert reason is not None and "unspent sigma_z" in reason and "fold by sigma_y" in reason
        assert "descent arrow" not in reason
        assert plane.domain_refusal(zfold) is None

    @pytest.mark.foundation
    def test_a_ONE_DIMENSIONAL_rule_is_refused_by_the_unspent_symmetry_alone(self):
        r"""The arrow ``S² → S²/O2_x`` EXISTS (the entry's own quotient map),
        so the arrow is not what refuses a polar rule for a 2-D geometry —
        ``O(2)_x ⊄ Γ·{e}`` is: its identity component ``SO(2)_x`` is not
        inside ``{e}``, the first conjunct of the coverage theorem."""
        for geometry in ("cylinder", "cartesian2d"):
            sym = GEOMETRY_ANGULAR_SYMMETRY[geometry]
            slab_rule = Quadrature.gauss_legendre(8).measure
            assert quotient_onto(sym.support, slab_rule.support) is not None
            assert not slab_rule.support.by.is_subset_of_product(gamma=sym.unspent, kappa=sym.spent)
            assert not sym.admits_domain(slab_rule)

    @pytest.mark.foundation
    def test_the_slab_still_admits_its_own_rule_through_the_identity_arrow(self):
        r"""⛔ The row a naive spelling breaks.  The slab OWES ``σ_x`` and
        ``σ_x ⊉ O(2)_x`` (`[M]` ``False``), so reading the owed closure as the
        fold licence refuses the slab's own Gauss-Legendre rule — which is
        why, until R3 of #434, the identity case had to short-circuit in
        front of that leg.  The coverage test needs no short circuit: the
        slab's own rule is ``O(2)_x ⊆ {e}·O(2)_x``, total and true."""
        slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        rule = Quadrature.gauss_legendre(8).measure
        assert rule.support == slab.support
        assert not slab.owed.contains(rule.support.by)               # the trap
        assert rule.support.by.is_subset_of_product(gamma=slab.unspent, kappa=slab.spent)
        assert slab.admits_domain(rule)
        for geometry in ("cylinder", "cartesian2d"):
            sym = GEOMETRY_ANGULAR_SYMMETRY[geometry]
            for factory in (
                lambda: Quadrature.product(4, 8).measure,
                lambda: Quadrature.level_symmetric(4).measure,
                lambda: Quadrature.lebedev(5).measure,
            ):
                assert sym.admits_domain(factory())

    @pytest.mark.foundation
    def test_the_rejection_message_names_the_one_failing_clause(self):
        """The 1-D rule on the cylinder HAS a descent arrow (``S² → S²/O2_x`` is
        the entry's own quotient map); what refuses it is the fold licence —
        ``O(2)_x ⊄ D_1h·{e}`` — and the reason says so, and ONLY so.  Until R3
        of #434 the message was a disjunction naming both causes on every
        refusal (`[M]` 14 of 17 shipped stage-0 refusals fail the arrow only,
        3 the licence only, 0 both) and this row was named for the arrow.
        """
        from orpheus.numerics.quadrature.registry import select_quadrature

        _, log = select_quadrature("cylinder", 5)
        reasons = [r for name, r in log.rejected if name == "GaussLegendre1D"]
        assert reasons, "GaussLegendre1D was not rejected for the cylinder"
        assert "domain mismatch" in reasons[0]
        assert "S^2/O2_x" in reasons[0]
        assert "fold by O2_x" in reasons[0] and "unspent D_1h" in reasons[0]
        assert "arrow" not in reasons[0]
        # ... and a rule with NO arrow is told exactly that, and nothing else
        _, log = select_quadrature("slab", 5)
        reasons = [r for name, r in log.rejected if name == "LebedevSphere"]
        assert reasons and "no descent arrow" in reasons[0]
        assert "unspent" not in reasons[0]


# =============================================================================


# G13 — tests/numerics/test_registry.py
# =============================================================================
class TestStageOneOnAFoldAsksTheOrbitSpaceNotTheRepresentatives:
    r"""(viii) Stage 1 on the fold — the positive leg and a real negative one.

    `[M]` before the carve ``cylinder.admits_symmetry(folded_product(4, 8))``
    was ``False``: Γ = D_2h was asked of the AMBIENT representatives and σ_y
    maps a y ≥ 0 node to its absent mate.  After, it is ``True`` — and the
    negative leg below proves the row is not simply "everything passes now".
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 8), (6, 8)])
    def test_the_cylinder_owed_symmetry_is_realized_on_the_fold(self, n_mu, n_phi):
        cyl = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
        assert cyl.admits_symmetry(Quadrature.folded_product(n_mu, n_phi).measure)

    @pytest.mark.foundation
    def test_a_fold_whose_representatives_are_not_closed_is_refused(self):
        cyl = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
        m = Quadrature.folded_product(4, 8).measure
        nodes = np.asarray(m.nodes, dtype=float)
        weights = np.asarray(m.weights, dtype=float)
        bumped = weights.copy()
        bumped[0] *= 1.5
        assert not cyl.admits_symmetry(
            DiscreteMeasure(nodes=nodes, weights=bumped, support=m.support)
        )
        assert not cyl.admits_symmetry(
            DiscreteMeasure(nodes=nodes[1:], weights=weights[1:], support=m.support)
        )

    @pytest.mark.foundation
    def test_both_stages_are_needed_for_the_fold(self):
        """Stage 0 and stage 1 are independent on the fold too: the slab
        admits neither, and a perturbed fold clears stage 0 and fails
        stage 1."""
        cyl = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
        slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        m = Quadrature.folded_product(4, 8).measure
        assert cyl.admits_domain(m) and cyl.admits_symmetry(m)
        assert not slab.admits_domain(m)
        bumped = np.asarray(m.weights, dtype=float).copy()
        bumped[0] *= 1.5
        broken = DiscreteMeasure(
            nodes=np.asarray(m.nodes, dtype=float), weights=bumped, support=m.support
        )
        assert cyl.admits_domain(broken) and not cyl.admits_symmetry(broken)


# =============================================================================


# G14 — tests/numerics/test_registry.py
# =============================================================================
class TestTheSelectorIsUNMOVEDByTheCarve:
    r"""(ix) ``select_quadrature`` end-to-end.

    ``folded_product`` is NOT registered (the design records the structural
    reasons), so **no end-to-end selection exercises the fold** and the carve
    cannot be credited with an end-to-end row.  What the selector owes
    instead is a REGRESSION statement: `[M]` 4 geometries × 4 degrees = 16
    rows, the chosen spec, its parameters and its point count are identical
    before and after the carve.
    """

    _EXPECTED = {
        ("slab", 1): ("GaussLegendre1D", {"n": 2}, 2),
        ("slab", 3): ("GaussLegendre1D", {"n": 2}, 2),
        ("slab", 5): ("GaussLegendre1D", {"n": 4}, 4),
        ("slab", 7): ("GaussLegendre1D", {"n": 4}, 4),
        ("sphere", 1): ("GaussLegendre1D", {"n": 2}, 2),
        ("sphere", 3): ("GaussLegendre1D", {"n": 2}, 2),
        ("sphere", 5): ("GaussLegendre1D", {"n": 4}, 4),
        ("sphere", 7): ("GaussLegendre1D", {"n": 4}, 4),
        ("cylinder", 1): ("ProductQuadrature", {"n_mu": 1, "n_phi": 2}, 2),
        ("cylinder", 3): ("LebedevSphere", {"order": 3}, 6),
        ("cylinder", 5): ("LebedevSphere", {"order": 5}, 14),
        ("cylinder", 7): ("LebedevSphere", {"order": 7}, 26),
        ("cartesian2d", 1): ("ProductQuadrature", {"n_mu": 1, "n_phi": 2}, 2),
        ("cartesian2d", 3): ("LebedevSphere", {"order": 3}, 6),
        ("cartesian2d", 5): ("LebedevSphere", {"order": 5}, 14),
        ("cartesian2d", 7): ("LebedevSphere", {"order": 7}, 26),
    }

    @pytest.mark.foundation
    @pytest.mark.parametrize("key", sorted(_EXPECTED))
    def test_the_chosen_rule_is_unmoved(self, key):
        from orpheus.numerics.quadrature.registry import select_quadrature

        geometry, degree = key
        name, params, n_points = self._EXPECTED[key]
        measure, log = select_quadrature(geometry, degree)
        assert log.chosen_spec.name == name
        assert log.chosen_parameters == params
        assert measure.n_points == n_points

    @pytest.mark.foundation
    def test_the_fold_is_not_registered_so_no_selection_reaches_it(self):
        from orpheus.numerics.quadrature.registry import quadrature_registry

        assert "folded" not in {s.name.lower() for s in quadrature_registry}
        supports = set()
        for spec in quadrature_registry:
            params = spec.degree_of_exactness_for(5)
            assert params is not None
            supports.add(spec.build(params).support.name)
        assert "S^2/sigma_y" not in supports

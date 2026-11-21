"""Foundation tests for ``orpheus.numerics.quadrature.registry``.

Verifies the catalogue of rules (the inversions of the degree formulas, the
frozen registry) and the geometry-to-symmetry table on the four canonical
geometry tags shipped today (``slab``, ``sphere``, ``cylinder``,
``cartesian2d``). The selector that once walked the catalogue retired on
2026-11-20 (wave W0, unit ``dead``); its law, the admission of a rule by a
geometry's symmetry, lives in the legs below.

The domain and symmetry legs are the two halves of :class:`AngularSymmetry`: the
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
``tests/gates/numerics/test_rules_*.py``.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.numerics.exactness import UNIFORM_ON_SPHERE
from orpheus.numerics.generating_measure import CHEBYSHEV_T, LEGENDRE, pushforward_reference
from orpheus.geometry.manifold import COSINE_INTERVAL, SPHERE, Quotient
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import (
    Quadrature,
    GEOMETRY_ANGULAR_SYMMETRY,
    AngularSymmetry,
    QuadratureSpec,
    quadrature_registry,
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
from orpheus.geometry.symmetry import SubgroupOfO3

# 2.2b gates (test-architect, 2026-09-02) — the draft's imports, verbatim
import numpy as np
import pytest

from orpheus.geometry.transformation import RigidMotion
from orpheus.geometry.manifold import (
    COSINE_INTERVAL,
    SPHERE,
    Ball,
    Quotient,
    quotient_onto,
)
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.registry import GEOMETRY_ANGULAR_SYMMETRY
from orpheus.numerics.quadrature.rules_1d import (
    gauss_legendre_on_mu,
    gauss_legendre_on_polar_orbit,
)
from orpheus.geometry.symmetry import SubgroupOfO3


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
    truthfully. The admission legs compute it from the built measure.

    ⭐ It carries no ``expected_node_count`` either, for the same reason
    one level down. The cost is a property of the rule the selector (retired on
    2026-11-20) had *already built* (its stages 0 and 1 forced the instantiation), so a
    formula-shaped second source could only ever drift from it. `[M]`
    2026-08-14 the two agreed on all 25 shipped configurations — which is
    what makes a twin dangerous rather than safe — and the first family
    that would have broken it already exists: ``folded_product``
    quotients by a mirror, so ``n_mu * n_phi`` over-counts it by 2x.

    Both absence assertions below are the load-bearing half of this
    gate: they are what stops either field being re-added by a
    contributor who reaches for a cheap tag.
    """
    for spec in quadrature_registry:
        assert isinstance(spec, QuadratureSpec)
        assert callable(spec.factory)
        assert callable(spec.degree_of_exactness_for)
        assert not hasattr(spec, "invariance_group")
        assert not hasattr(spec, "expected_node_count")
        assert not hasattr(spec, "parameters")  # the classmethod's signature is the one list


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
    assert params == {"n_ordinates": expected_n}


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
        factory=lambda n_ordinates: Quadrature(measure=_on_slab(CHEBYSHEV_T.gauss(n_ordinates)), level_structure=None),
        degree_of_exactness_for=_gl1d_invert,   # same deg = 2n - 1
        axis_aligned=True,
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
    gc = _gauss_chebyshev_spec().build({"n_ordinates": 4}).measure

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
        claim = spec.build(params).measure.exactness
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

    ``name`` is not merely a label: it is the key a spec is addressed by.
    ``dict`` keeps the LAST value for a repeated key, so two rules sharing
    a name make one of them **silently vanish** from a view keyed on it,
    with nothing warning.

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
    ).build({"order": 5}).measure
    odd_product = next(
        s for s in quadrature_registry if s.name == "ProductQuadrature"
    ).build({"n_mu": 3, "n_phi": 5}).measure

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
        measure = product.build({"n_mu": 3, "n_phi": n_phi}).measure
        assert cylinder.admits_symmetry(measure) == (n_phi % 2 == 0), (
            f"n_phi={n_phi}: parity is the discriminator"
        )


@pytest.mark.foundation
def test_the_symmetry_leg_asks_the_nodes_not_the_declared_tag() -> None:
    r"""The symmetry leg (``admits_symmetry``) must not route through ``measure.invariance_group``.

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
    ).build({"n_ordinates": 8}).measure

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



@pytest.mark.foundation
def test_the_reference_is_READ_off_the_catalogue_entry() -> None:
    r"""Pattern 2 (#429 tracker 3.1, 2026-09-02): the hat-box answer has ONE
    home, and it is the orbit space's own derivation.

    ``AngularSymmetry.reference`` used to TABULATE ``LEGENDRE`` for any
    rotation axis; it now reads
    :func:`~orpheus.numerics.generating_measure.pushforward_reference` of the
    orbit space :attr:`support` derives — the same collapse ``support``
    underwent at tracker 2.4. Three legs: the axial geometry reads that answer by
    IDENTITY; the geometry that spends nothing has the bare sphere as its
    domain, so its reference is the sphere's own Lebesgue measure (a
    property of the base, not a catalogue read — user-ruled); and the
    arm for an entry whose pushforward no shipped realization spells
    RAISES naming the missing work (``plan-authoring`` §6c: the gate lands
    with its witness — a constructed geometry spending a mirror).
    """
    from orpheus.geometry.manifold import Quotient

    slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
    assert isinstance(slab.support, Quotient)
    assert slab.reference is pushforward_reference(slab.support)
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
    assert pushforward_reference(mirror_entry) is None
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
def test_the_stabiliser_orbit_space_admits_GaussLegendre1D_and_refuses_the_three_S2_rules(
    geometry: str,
) -> None:
    r"""The domain leg: ``Quadrature.gauss_legendre(4)`` (the rule the retired
    selector chose at ``target_degree=5``) lives on ``S^2/O2_x`` and is
    admitted; the three :math:`S^2` rules are refused on DOMAIN with a
    message naming the entry.

    `[M]` 2026-09-02 at ``target_degree=5``: ``GaussLegendre1D`` with
    ``n = 4``, measure on ``S^2/O2_x`` carrying ``invariance_group ==
    sigma_x``, and all three sphere rules refused as *"geometry '…'
    discretises S^2/O2_x, but the rule's nodes live on S^2"*.

    ⚠ The refusal REASON is asserted because it is the one place the
    entry's NAME reaches a user-facing string; a domain comparison that had
    kept comparing charts would still admit GL1D and would print the old
    name, so the admission alone is not a discriminating observation.
    """
    symmetry = GEOMETRY_ANGULAR_SYMMETRY[geometry]
    measure = Quadrature.gauss_legendre(4).measure
    assert symmetry.admits_domain(measure)
    assert measure.support == SPHERE.quotient(SubgroupOfO3.O2("x"))
    assert measure.invariance_group == SubgroupOfO3.Mirror("x")

    for name, rule in (
        ("LebedevSphere", Quadrature.lebedev(5)),
        ("LevelSymmetricSN", Quadrature.level_symmetric(4)),
        ("ProductQuadrature", Quadrature.product(n_mu=3, n_phi=5)),
    ):
        reason = symmetry.domain_refusal(rule.measure)
        assert reason is not None, name
        assert "S^2/O2_x" in reason, f"{name}: {reason}"


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


# G12 — tests/gates/numerics/test_registry.py
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

    ⚠ **And the honest scope of the coverage leg at the selector tier (the
    selector retired on 2026-11-20).** `[M]` at the selector it changed the
    REASON and never the selection: ``GaussLegendre1D`` IS registered and IS
    a fold (``S^2/O2_x``, ``H = O(2)_x``), and the shipped log for the
    cylinder at degree 5 refused it at stage 0 by the coverage clause (*"a
    fold by O2_x … unspent D_1h, spent Trivial"*); neutering the leg
    in-process moved the refusal to stage 2's V conjunct while the choice,
    ``LebedevSphere(5)``, stayed — 0 of 48 (geometry × degree) choices moved.
    So the leg's gate belongs HERE, at the ``admits_domain`` tier, and no
    end-to-end row may be credited to it.  (Until 2026-09-03 this paragraph said
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
        reason = GEOMETRY_ANGULAR_SYMMETRY["cylinder"].domain_refusal(
            Quadrature.gauss_legendre(8).measure
        )
        assert reason is not None, "GaussLegendre1D was not refused for the cylinder"
        assert "S^2/O2_x" in reason
        assert "fold by O2_x" in reason and "unspent D_1h" in reason
        assert "arrow" not in reason
        # ... and a rule with NO arrow is told exactly that, and nothing else
        reason = GEOMETRY_ANGULAR_SYMMETRY["slab"].domain_refusal(Quadrature.lebedev(5).measure)
        assert reason is not None and "no descent arrow" in reason
        assert "unspent" not in reason


# =============================================================================


# G13 — tests/gates/numerics/test_registry.py
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


# G14 — tests/gates/numerics/test_registry.py
# =============================================================================
class TestTheFoldIsNotRegistered:
    r"""(ix) ``folded_product`` is NOT registered (the design records the
    structural reasons), so **no shipped catalogue entry reaches the fold**
    and the carve cannot be credited with an end-to-end row.  (The selector
    that once chose among the entries, and whose 16-row regression table
    pinned that choice, retired on 2026-11-20.)
    """

    @pytest.mark.foundation
    def test_the_fold_is_not_registered_so_no_entry_reaches_it(self):
        from orpheus.numerics.quadrature.registry import quadrature_registry

        assert "folded" not in {s.name.lower() for s in quadrature_registry}
        supports = set()
        for spec in quadrature_registry:
            params = spec.degree_of_exactness_for(5)
            assert params is not None
            supports.add(spec.build(params).measure.support.name)
        assert "S^2/sigma_y" not in supports


@pytest.mark.foundation
def test_lebedev_order_13_is_not_positive_weighted() -> None:
    """The fact ERR-118 recorded, read off the built rule: Lebedev order 13 carries a non-positive weight (3 of 32 orders)."""
    assert not np.all(Quadrature.lebedev(order=13).weights > 0.0)


# =============================================================================
# Re-filed from tests/gates/geometry/test_manifold.py (layer DAG ruling 23, 2026-10-25)
# =============================================================================


@pytest.mark.foundation
class TestQuotient:
    """The numerics half of geometry's ``TestQuotient``.

    Re-filed 2026-10-25 from ``tests/gates/geometry/test_manifold.py`` (layer-DAG
    ruling 23): the row below asserts a numerics object's verb, so it lives here; the class keeps its geometry name so the node id
    changes only in its path.  The manifold rows of the class stay in geometry.
    """

    def test_the_derived_orbit_space_agrees_with_the_hand_written_table(self):
        r"""⭐ Two independent producers of one orbit space must agree.

        ``AngularSymmetry.support`` (``quadrature/registry.py``) was a
        **hand-written table** — ``SO2 -> COSINE_INTERVAL``,
        ``Trivial -> SPHERE``, raise otherwise.  :meth:`GSpace.quotient`
        **derives** the same orbit space from invariant theory (the
        Procesi--Schwarz matrix of the invariant ring).  Neither read the
        other, so this was a genuine two-implementation agreement and the
        only check that either was right about :math:`S^2/G^0`.

        ⛔ **RE-SCOPED 2026-09-01 by tracker 2.4, which single-sourced the
        SO(2) row.** The registry now DERIVES its domain through
        ``SPHERE.quotient(spent)`` — so for the axial row ``theirs`` is the
        very object ``mine`` was built from, and asking whether they agree
        is asking one call whether it equals itself (``coding-standards``:
        single-sourcing a duplicate demotes every gate that compared its
        copies, and the fix stays; the gate's DESCRIPTION is what moves).
        What survives, and is asserted below: (i) the hand-written
        ``expected`` column is still an INDEPENDENT pin on the derivation's
        realization — the one input that could make it fail is a wrong
        Procesi–Schwarz derivation; (ii) the ``Trivial`` row is still two
        producers, because the registry answers ``SPHERE`` explicitly while
        the catalogue re-derives ``M/{e}``; (iii) the axial row's registry
        answer must be the ORBIT SPACE, carrying the spent group, not its
        chart — the reading a spatial rule could share.

        ⭐ **Strengthened 2026-09-01 by the very retype it was written to
        pin** (tracker 2.0c).  It used to compare ``mine.name == theirs``
        because ``theirs`` was a *string* — the strongest claim the string
        vocabulary admitted.  Both sides are now ``Manifold``, so the gate
        asserts **object equality**: not merely that two producers spell the
        orbit space the same way, but that they produce the same point set.
        A name gate is satisfied by any self-consistent lie; this one is not.
        (``coding-standards``' mirror clause: a retirement can silently
        PROMOTE a gate's claim class, and the docstring must move with it or
        it goes on advertising the weaker claim.)

        ⚠ If this ever reddens, do NOT re-baseline: one of the two is wrong
        about a quotient, and which one is a mathematical question, not a
        test-maintenance one.
        """
        from orpheus.numerics.quadrature.registry import AngularSymmetry

        rows = {
            SubgroupOfO3.O2("x"): COSINE_INTERVAL,
            SubgroupOfO3.O2("z"): COSINE_INTERVAL,
            SubgroupOfO3.Trivial: SPHERE,
        }
        for group, expected in rows.items():
            derived = SPHERE.quotient(group)
            mine = derived.realization
            theirs = AngularSymmetry(
                spent=group, unspent=SubgroupOfO3.Trivial, owed=SubgroupOfO3.Trivial,
            ).support
            # (i) the independent pin on the derivation
            assert mine == expected, f"{group.name}: mine gave {mine.name}"
            if group == SubgroupOfO3.Trivial:
                # (ii) still two producers: an explicit SPHERE vs M/{e}
                assert theirs == SPHERE and mine == theirs
            else:
                # (iii) the registry hands out the orbit space itself
                assert isinstance(theirs, Quotient) and theirs == derived, (
                    f"AngularSymmetry.support for S^2/{group.name} is "
                    f"{theirs!r}, not the derived orbit space"
                )
                assert theirs.by == group and theirs.realization == expected
